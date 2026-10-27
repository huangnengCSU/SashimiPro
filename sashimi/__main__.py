"""Command line: `python -m sashimi app ...` (interactive) or `python -m sashimi plot ...` (static file)."""

import argparse
import json
import os
import sys

from . import __version__
from .data import LIBTYPES, ReadFilter, format_region, junction_table, load_region, parse_region, read_samples
from .plot import PALETTES


def _common(sp):
    sp.add_argument("-b", "--bam", nargs="+", help="Sample sheet TSV (id, bam path, extra columns) or BAM file(s)")
    sp.add_argument("-c", "--coordinates", "--region", dest="region", help="chr:start-end, or a gene name with -g")
    sp.add_argument("-g", "--gtf", help="GTF annotation (.gtf or .gtf.gz)")
    sp.add_argument("-s", "--strand", dest="libtype", default="NONE", choices=LIBTYPES,
                    help="Library strandedness [%(default)s]")
    sp.add_argument("--min-mapq", type=int, default=0, help="Minimum mapping quality [%(default)s]")
    sp.add_argument("--skip-duplicates", action="store_true", help="Ignore reads flagged as duplicates")
    sp.add_argument("--params", help="Settings JSON saved from the app (plot parameters, region, strandedness)")


def _column(sheet_cols, value):
    """Column given by name or by 1-based index (ggsashimi style)."""
    if value is None:
        return None
    if value.isdigit():
        i = int(value) - 1
        if not 0 <= i < len(sheet_cols):
            sys.exit(f"ERROR: column index {value} out of range (sheet has {len(sheet_cols)} columns)")
        return sheet_cols[i]
    if value not in sheet_cols:
        sys.exit(f"ERROR: column '{value}' not in sample sheet {list(sheet_cols)}")
    return value


def main(argv=None):
    ap = argparse.ArgumentParser(prog="sashimi", description="Interactive sashimi plots (pure Python).")
    ap.add_argument("--version", action="version", version=f"SashimiPro {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("app", help="Start the interactive browser app")
    _common(a)
    a.add_argument("--host", default="127.0.0.1", help="Use 0.0.0.0 to serve on the network [%(default)s]")
    a.add_argument("--port", type=int, default=8050)
    a.add_argument("--debug", action="store_true")

    p = sub.add_parser("plot", help="Render a static figure (pdf/svg/png), like ggsashimi")
    _common(p)
    p.add_argument("-o", "--out", default="sashimi.pdf", help="Output file; format from extension [%(default)s]")
    p.add_argument("-S", "--out-strand", default="both", choices=["both", "plus", "minus"])
    p.add_argument("-M", "--min-coverage", type=float, help="Minimum reads for a junction to be drawn")
    p.add_argument("-O", "--overlay", help="Column (name or 1-based index) to group/overlay samples by")
    p.add_argument("-A", "--aggr", choices=["mean", "median", "sum", "mean_j", "median_j"],
                   help="Aggregate group coverage + junctions; *_j keeps overlaid coverage, aggregates junctions only")
    p.add_argument("-C", "--color-factor", help="Column (name or 1-based index) with color levels")
    p.add_argument("-L", "--labels", help="Column (name or 1-based index) with track labels")
    p.add_argument("-P", "--palette", help=f"Palette name {list(PALETTES)} or a file with one color per line")
    p.add_argument("--alpha", type=float)
    p.add_argument("--shrink", action="store_true", default=None)
    p.add_argument("--fix-y-scale", action="store_true", default=None)
    p.add_argument("--cpm", action="store_true", default=None, help="Normalise to reads per million mapped")
    p.add_argument("--height", type=int, help="Track height in px")
    p.add_argument("--width", type=int, help="Figure width in px")
    p.add_argument("--base-size", type=float, help="Font size")
    p.add_argument("--highlight", help="Introns to highlight, e.g. 'chr7:128,940,001-128,940,500; 128941001-128941200'")
    p.add_argument("--scale", type=float, default=2, help="Raster scale factor for png [%(default)s]")
    p.add_argument("-j", "--junctions-tsv", help="Also write all junction counts to this TSV")

    args = ap.parse_args(argv)
    settings = {}
    if args.params:
        with open(args.params) as fh:
            settings = json.load(fh)
    region = args.region or settings.get("region")
    libtype = args.libtype if args.libtype != "NONE" or not settings.get("libtype") else settings["libtype"]

    if args.cmd == "app":
        from .app import run
        init = dict(bam=" ".join(args.bam or []), gtf=args.gtf or "", region=region or "", libtype=libtype,
                    mapq=args.min_mapq, params=settings.get("params"))
        print(f"SashimiPro {__version__}: open http://{args.host if args.host != '0.0.0.0' else 'localhost'}:{args.port}")
        run(init, host=args.host, port=args.port, debug=args.debug)
        return

    # ---- static plot ----
    from .annotation import Annotation
    from .plot import build_figure

    if not args.bam or not region:
        ap.error("plot needs -b/--bam and -c/--coordinates (or --params with a region)")
    samples = read_samples(args.bam)
    ann = Annotation(args.gtf) if args.gtf else None
    try:
        chrom, start, end = parse_region(region)
    except ValueError:
        hit = ann.find_gene(region) if ann else None
        if not hit:
            sys.exit(f"ERROR: cannot parse region / find gene '{region}'")
        chrom, start, end = hit
        flank = max(50, int((end - start) * 0.1))
        start, end = max(1, start - flank), end + flank

    rf = ReadFilter(min_mapq=args.min_mapq, skip_duplicates=args.skip_duplicates)
    data = load_region(samples, chrom, start, end, libtype, rf)
    for w in data.warnings:
        print("WARN:", w, file=sys.stderr)
    transcripts = ann.region(chrom, start, end) if ann else None

    params = dict(settings.get("params") or {})
    cols = list(samples.columns)
    if args.overlay:
        params["group_by"] = _column(cols, args.overlay)
    if args.color_factor:
        params["color_by"] = _column(cols, args.color_factor)
    if args.labels:
        params["label_by"] = _column(cols, args.labels)
    if args.aggr:
        base = args.aggr.replace("_j", "")
        params["group_mode"] = "overlay" if args.aggr.endswith("_j") else base
        params["junction_agg"] = base
    elif args.overlay and "junction_agg" not in params:
        params["junction_agg"] = "each"
    if args.palette:
        if os.path.isfile(args.palette):
            with open(args.palette) as fh:
                params["custom_colors"] = ",".join(l.split("\t")[0].strip() for l in fh if l.strip())
        else:
            params["palette"] = args.palette
    for k, v in (("min_junction", args.min_coverage), ("alpha", args.alpha), ("shrink", args.shrink),
                 ("fix_y", args.fix_y_scale), ("track_height", args.height), ("width", args.width),
                 ("font_size", args.base_size)):
        if v is not None:
            params[k] = v
    if args.cpm:
        params["normalize"] = "cpm"
    if args.highlight:
        params["highlight"] = args.highlight
    params.setdefault("width", 1000)

    root, ext = os.path.splitext(args.out)
    fmt = ext.lstrip(".").lower() or "pdf"
    fmt = {"jpg": "jpeg", "tif": "tiff"}.get(fmt, fmt)
    single = f"{root}.{fmt}"
    if libtype == "NONE":
        outs = [("combined", single)]
    elif args.out_strand != "both":
        outs = [({"plus": "+", "minus": "-"}[args.out_strand], single)]
    elif "strand" in params:            # strand chosen in the app's saved settings
        outs = [(params["strand"], single)]
    else:                               # stranded library: one file per strand, like ggsashimi
        outs = [("+", f"{root}_plus.{fmt}"), ("-", f"{root}_minus.{fmt}")]
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    for strand, path in outs:
        fig = build_figure(data, dict(params, strand=strand), transcripts, interactive=False)
        fig.write_image(path, scale=args.scale if fmt in ("png", "jpeg") else 1)
        print(f"Wrote {path}  ({format_region(chrom, start, end)})")
    if args.junctions_tsv:
        junction_table(data).to_csv(args.junctions_tsv, sep="\t", index=False)
        print(f"Wrote {args.junctions_tsv}")


if __name__ == "__main__":
    main()
