"""GTF annotation: loaded once, then queried per region and by gene name."""

import gzip
import hashlib
import os
import re

import numpy as np
import pandas as pd

_ATTRS = ["gene_id", "gene_name", "transcript_id", "transcript_name", "transcript_type", "transcript_biotype"]
_CACHE_DIR = os.path.join(os.path.expanduser("~"), ".cache", "sashimi_pro")


def _cache_path(path):
    st = os.stat(path)
    key = f"{os.path.abspath(path)}|{st.st_size}|{st.st_mtime}"
    return os.path.join(_CACHE_DIR, hashlib.md5(key.encode()).hexdigest() + ".pkl")


def _parse_gtf(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 6, 8],
                         names=["chrom", "feature", "start", "end", "strand", "attr"],
                         dtype={0: str, 2: str, 3: np.int64, 4: np.int64, 6: str, 8: str})
    df = df[df["feature"].isin(["exon", "CDS"])].reset_index(drop=True)
    for a in _ATTRS:
        df[a] = df["attr"].str.extract(rf'{a} "([^"]*)"', expand=False)
    if df["transcript_id"].isna().all():
        raise ValueError("GTF has no 'transcript_id' attributes on exon lines")
    df["gene_name"] = df["gene_name"].fillna(df["gene_id"])
    df["transcript_name"] = df["transcript_name"].fillna(df["transcript_id"])
    df["biotype"] = df["transcript_type"].combine_first(df["transcript_biotype"]).fillna("").astype(str)
    df = df.drop(columns=["attr", "transcript_type", "transcript_biotype"])
    df["chrom"] = df["chrom"].astype("category")
    df["feature"] = df["feature"].astype("category")
    return df


class Annotation:
    """In-memory exon/CDS table.  Large GTFs are parsed once and pickled to ~/.cache."""

    def __init__(self, path):
        self.path = path
        cache = _cache_path(path)
        if os.path.exists(cache):
            self.df = pd.read_pickle(cache)
        else:
            self.df = _parse_gtf(path)
            try:
                os.makedirs(_CACHE_DIR, exist_ok=True)
                self.df.to_pickle(cache)
            except OSError:
                pass
        ex = self.df[self.df["feature"] == "exon"]
        self.genes = (ex.groupby(["gene_name", "chrom"], observed=True)
                        .agg(start=("start", "min"), end=("end", "max"), strand=("strand", "first"))
                        .reset_index())
        self._upper = self.genes["gene_name"].str.upper()
        self.chroms = set(self.df["chrom"].cat.categories)

    def _chrom(self, chrom):
        if chrom in self.chroms:
            return chrom
        alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        return alt if alt in self.chroms else chrom

    def find_gene(self, name):
        """Gene name / id -> (chrom, start, end) or None."""
        hits = self.genes[self._upper == name.strip().upper()]
        if hits.empty:
            gid = self.df.loc[self.df["gene_id"].str.split(".").str[0] == name.strip().split(".")[0], "gene_name"]
            if gid.empty:
                return None
            hits = self.genes[self.genes["gene_name"] == gid.iloc[0]]
        h = hits.iloc[0]
        return str(h["chrom"]), int(h["start"]), int(h["end"])

    def region(self, chrom, start, end):
        """Transcripts overlapping the region.

        Returns a list of dicts (sorted by gene, then start) with keys
        id, name, gene, strand, start, end, biotype, exons [(s, e)], cds [(s, e)].
        """
        chrom = self._chrom(chrom)
        df = self.df[(self.df["chrom"] == chrom)]
        tx_ids = df.loc[(df["end"] >= start) & (df["start"] <= end), "transcript_id"].unique()
        sub = df[df["transcript_id"].isin(tx_ids)]
        out = []
        for tid, g in sub.groupby("transcript_id", sort=False):
            exons = sorted(zip(g.loc[g["feature"] == "exon", "start"], g.loc[g["feature"] == "exon", "end"]))
            if not exons:
                continue
            cds = sorted(zip(g.loc[g["feature"] == "CDS", "start"], g.loc[g["feature"] == "CDS", "end"]))
            r = g.iloc[0]
            out.append(dict(id=tid, name=r["transcript_name"], gene=r["gene_name"], strand=r["strand"],
                            start=exons[0][0], end=exons[-1][1], biotype=r["biotype"],
                            exons=[(int(a), int(b)) for a, b in exons], cds=[(int(a), int(b)) for a, b in cds]))
        out.sort(key=lambda t: (t["gene"], t["start"], -(t["end"] - t["start"])))
        return out


def collapse_genes(transcripts):
    """Merge the transcripts of each gene into one 'meta-transcript'."""
    by_gene = {}
    for t in transcripts:
        by_gene.setdefault(t["gene"], []).append(t)
    out = []
    for gene, ts in by_gene.items():
        merged = {}
        for key in ("exons", "cds"):
            iv = sorted(e for t in ts for e in t[key])
            m = []
            for a, b in iv:
                if m and a <= m[-1][1] + 1:
                    m[-1] = (m[-1][0], max(m[-1][1], b))
                else:
                    m.append((a, b))
            merged[key] = m
        out.append(dict(id=gene, name=gene, gene=gene, strand=ts[0]["strand"],
                        start=min(t["start"] for t in ts), end=max(t["end"] for t in ts),
                        biotype="", **merged))
    return out


def annotated_introns(transcripts):
    """Set of (intron_start, intron_end) implied by the transcripts' exons."""
    s = set()
    for t in transcripts:
        ex = t["exons"]
        for (a0, b0), (a1, b1) in zip(ex, ex[1:]):
            s.add((b0 + 1, a1 - 1))
    return s
