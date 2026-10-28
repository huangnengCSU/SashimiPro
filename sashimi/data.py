"""Reading alignments: per-base coverage and splice junctions for a region.

Coordinates exposed by this module are 1-based and inclusive, like the
region strings users type ("chr1:1,001-9,000").  A junction (intron) is
stored as (first intron base, last intron base), the STAR SJ.out convention.
"""

import os
import re
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import pysam

LIBTYPES = ["NONE", "SENSE", "ANTISENSE", "MATE1_SENSE", "MATE2_SENSE"]

# CIGAR operations (pysam codes)
_M, _I, _D, _N, _S, _H, _P, _EQ, _X = range(9)
_REF_CONSUMING = (_M, _D, _N, _EQ, _X)


def parse_region(text):
    """'chr1:1,001-9,000' -> ('chr1', 1001, 9000).  Raises ValueError."""
    m = re.fullmatch(r"\s*([^:\s]+):([\d,]+)\s*-\s*([\d,]+)\s*", text or "")
    if not m:
        raise ValueError(f"Cannot parse region '{text}' (expected chr:start-end)")
    chrom, start, end = m.group(1), int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))
    if start < 1 or end < start:
        raise ValueError(f"Invalid region '{text}'")
    return chrom, start, end


def format_region(chrom, start, end):
    return f"{chrom}:{start:,}-{end:,}"


# --------------------------------------------------------------------------
# Sample sheet
# --------------------------------------------------------------------------

def _is_alignment_path(s):
    return s.lower().endswith((".bam", ".cram", ".sam"))


def read_samples(paths):
    """Build the sample table from a TSV sample sheet or a list of BAM files.

    The TSV follows ggsashimi: col1 = id, col2 = bam path (relative paths are
    resolved against the TSV's directory), further columns are free metadata.
    A header line is optional; it is detected when column 2 is not a
    .bam/.cram path.  Without a header, metadata columns are named col3, col4...
    """
    if isinstance(paths, str):
        paths = [paths]
    if len(paths) == 1 and not _is_alignment_path(paths[0]):
        sheet = paths[0]
        with open(sheet, encoding="utf-8") as fh:
            lines = [l.rstrip("\n\r") for l in fh if l.strip() and not l.startswith("#")]
        if not lines:
            raise ValueError(f"Sample sheet {sheet} is empty")
        # Tab-separated is canonical; fall back to any whitespace when no tabs.
        if any("\t" in l for l in lines):
            rows = [l.split("\t") for l in lines]
        else:
            rows = [l.split() for l in lines]
        if all(len(r) < 2 for r in rows):
            raise ValueError(f"Sample sheet {sheet} needs at least two columns (id, bam path)")
        if len(rows[0]) > 1 and not _is_alignment_path(rows[0][1]):
            header, rows = rows[0], rows[1:]
        else:
            header = ["id", "bam"] + [f"col{i + 1}" for i in range(2, max(map(len, rows)))]
        header = [h.strip() or f"col{i + 1}" for i, h in enumerate(header)]
        header[0], header[1] = "id", "bam"
        width = len(header)
        rows = [(r + [""] * width)[:width] for r in rows]
        df = pd.DataFrame(rows, columns=header)
        missing = df.loc[df["bam"].str.strip() == "", "id"].tolist()
        if missing:
            raise ValueError(f"Sample sheet {sheet}: no bam path for {missing}")
        base = os.path.dirname(os.path.abspath(sheet))
        df["bam"] = [p if os.path.isabs(p) else os.path.join(base, p) for p in df["bam"]]
    else:
        ids = [re.sub(r"\.(bam|cram|sam)$", "", os.path.basename(p), flags=re.I) for p in paths]
        df = pd.DataFrame({"id": ids, "bam": [os.path.abspath(p) for p in paths]})
    if df["id"].duplicated().any():
        dup = df.loc[df["id"].duplicated(), "id"].tolist()
        raise ValueError(f"Duplicated sample ids: {dup}")
    return df.reset_index(drop=True)


# --------------------------------------------------------------------------
# Alignment counting
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class ReadFilter:
    min_mapq: int = 0
    skip_duplicates: bool = False
    skip_secondary: bool = True
    skip_supplementary: bool = False
    skip_qcfail: bool = True


@dataclass
class SampleSignal:
    """Coverage and junctions of one sample, per strand ('+', '-')."""
    coverage: dict                      # strand -> np.ndarray[int32], index 0 == region start
    junctions: dict                     # strand -> {(intron_start, intron_end): count}
    total_mapped: int                   # mapped reads in the whole BAM (from the index)
    n_reads: int = 0                    # reads counted in the region


def _read_strand(read, libtype):
    if libtype == "NONE":
        return "+"
    flip = libtype == "ANTISENSE" \
        or (libtype == "MATE1_SENSE" and read.is_read2) \
        or (libtype == "MATE2_SENSE" and read.is_read1)
    return "-" if read.is_reverse != flip else "+"


def count_region(bam_path, chrom, start, end, libtype="NONE", read_filter=ReadFilter(), reference=None):
    """Per-base coverage and junctions of the reads overlapping chrom:start-end."""
    strands = ["+"] if libtype == "NONE" else ["+", "-"]
    length = end - start + 1
    origin = start - 1                  # 0-based coordinate of array index 0
    block_starts = {s: [] for s in strands}
    block_ends = {s: [] for s in strands}
    junctions = {s: {} for s in strands}
    rf = read_filter
    n = 0

    with pysam.AlignmentFile(bam_path, reference_filename=reference) as bam:
        if chrom not in bam.references:
            alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
            if alt not in bam.references:
                raise ValueError(f"{os.path.basename(bam_path)}: contig '{chrom}' not in BAM header")
            chrom = alt
        try:
            total = bam.mapped
        except ValueError:
            total = 0
        for read in bam.fetch(chrom, start - 1, end):
            if read.is_unmapped or read.mapping_quality < rf.min_mapq:
                continue
            if (rf.skip_secondary and read.is_secondary) or (rf.skip_supplementary and read.is_supplementary) \
                    or (rf.skip_duplicates and read.is_duplicate) or (rf.skip_qcfail and read.is_qcfail):
                continue
            n += 1
            s = _read_strand(read, libtype)
            # Aligned blocks (M/=/X) are what we count as coverage; D and N are gaps.
            bs, be = block_starts[s], block_ends[s]
            for b0, b1 in read.get_blocks():
                bs.append(b0)
                be.append(b1)
            pos = read.reference_start
            jd = junctions[s]
            for op, ln in read.cigartuples:
                if op == _N:
                    key = (pos + 1, pos + ln)       # 1-based first/last intron base
                    if key[0] > start and key[1] < end:
                        jd[key] = jd.get(key, 0) + 1
                if op in _REF_CONSUMING:
                    pos += ln

    coverage = {}
    for s in strands:
        diff = np.zeros(length + 1, dtype=np.int64)
        if block_starts[s]:
            b0 = np.clip(np.asarray(block_starts[s]) - origin, 0, length)
            b1 = np.clip(np.asarray(block_ends[s]) - origin, 0, length)
            np.add.at(diff, b0, 1)
            np.add.at(diff, b1, -1)
        coverage[s] = np.cumsum(diff[:-1]).astype(np.int32)
        junctions[s] = OrderedDict(sorted(junctions[s].items()))
    return SampleSignal(coverage, junctions, total, n)


@dataclass
class RegionData:
    chrom: str
    start: int
    end: int
    libtype: str
    samples: pd.DataFrame               # sample sheet rows that were loaded
    signals: dict                       # sample id -> SampleSignal
    warnings: list = field(default_factory=list)

    @property
    def strands(self):
        return ["+"] if self.libtype == "NONE" else ["+", "-"]


_CACHE = OrderedDict()
_CACHE_SIZE = 256


def _cached_count(bam, chrom, start, end, libtype, rf, reference):
    try:
        mtime = os.path.getmtime(bam)
    except OSError:
        mtime = None
    key = (bam, mtime, chrom, start, end, libtype, rf, reference)
    if key in _CACHE:
        _CACHE.move_to_end(key)
        return _CACHE[key]
    sig = count_region(bam, chrom, start, end, libtype, rf, reference)
    _CACHE[key] = sig
    while len(_CACHE) > _CACHE_SIZE:
        _CACHE.popitem(last=False)
    return sig


def load_region(samples, chrom, start, end, libtype="NONE", read_filter=ReadFilter(), reference=None, threads=4):
    """Count every sample of the sheet over the region (in parallel, cached)."""
    if libtype not in LIBTYPES:
        raise ValueError(f"Unknown library type {libtype}; choose from {LIBTYPES}")
    warnings, keep = [], []
    for i, row in samples.iterrows():
        if os.path.isfile(row["bam"]):
            keep.append(i)
        else:
            warnings.append(f"{row['id']}: file not found ({row['bam']})")
    samples = samples.loc[keep].reset_index(drop=True)
    if samples.empty:
        raise ValueError("No readable BAM files. " + "; ".join(warnings))

    def work(bam):
        return _cached_count(bam, chrom, start, end, libtype, read_filter, reference)

    with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
        results = list(ex.map(work, samples["bam"]))
    signals = OrderedDict()
    for sid, sig in zip(samples["id"], results):
        if sig.n_reads == 0:
            warnings.append(f"{sid}: no reads in region")
        signals[sid] = sig
    return RegionData(chrom, start, end, libtype, samples, signals, warnings)


def junction_table(region, min_count=1):
    """All junctions of all samples as a tidy DataFrame (for BED/TSV export)."""
    rows = []
    for sid, sig in region.signals.items():
        for strand, jd in sig.junctions.items():
            for (a, b), c in jd.items():
                if c >= min_count:
                    rows.append((region.chrom, a, b, sid, c, strand if region.libtype != "NONE" else "."))
    return pd.DataFrame(rows, columns=["chrom", "intron_start", "intron_end", "sample", "count", "strand"])
