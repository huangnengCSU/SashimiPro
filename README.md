# SashimiPro

Interactive sashimi plots for RNA-seq splicing, in pure Python (pysam + Plotly/Dash) with no R dependency.
Load your BAMs once, then adjust grouping, colors, junction filters, intron shrinking and more in the browser; the plot updates instantly.

📖 User manual (every parameter explained, plus recipes): [English](docs/manual_en.md) · [中文](docs/manual_zh.md)

## Features

- **Interactive browser app**: every display setting redraws the plot instantly and keeps your zoom; BAMs are only re-read when the region or read filters change.
- **Navigate by gene name** (with a GTF), pan left/right, zoom out, or reload exactly the range you zoomed into.
- **Group samples** by any column of the sample sheet: overlay replicates, or show the mean / median / sum per group.
- **Highlight introns of interest**: click an arc (or type coordinates) to emphasize it across all tracks.
- **Annotated vs novel junctions** from the GTF (dotted lines or filtering), CDS/UTR-aware transcript models, collapsed gene view.
- **CPM normalization**, MAPQ / duplicate / secondary-alignment filters, strand-specific libraries.
- **Hover** for coverage values and junction coordinates, length and read counts.
- **Export** PDF / SVG / PNG figures, per-sample junction counts (TSV), and all plot settings (JSON) for reproducible batch plotting.
- Coverage is identical, base by base, to `samtools depth`.

## Installation

Requires Python ≥ 3.9. Installing provides the `sashimi` command (`sashimipro` is an equivalent alias):

```bash
git clone https://github.com/huangnengCSU/SashimiPro.git
cd SashimiPro
pip install -e .
sashimi --version
```

Dependencies (pysam, numpy, pandas, plotly, dash, kaleido) are installed automatically. To keep them out of your main environment, install into a virtual environment first:

```bash
python -m venv .venv          # add --system-site-packages to reuse an existing pysam (e.g. from conda)
.venv/bin/pip install -e .
```

PDF/SVG/PNG export uses [kaleido](https://github.com/plotly/Kaleido), which needs Google Chrome or Chromium installed.

## Quick start: interactive app

```bash
sashimi app -b samples.tsv -g genes.gtf -c chr7:128,938,246-128,942,070 -s MATE2_SENSE
# then open http://127.0.0.1:8050
```

- Press **Load** (or Enter in the region box) to read the BAMs; all other settings apply immediately.
- The region box also accepts a gene name or Ensembl gene id (requires a GTF); ±10% flanks are added.
- ◀ ▶ move half a window; **zoom out** widens the region 3×; drag a box on the plot and press **zoom in** to reload just that range.
- Click a junction arc to highlight that intron in every track; click again to remove it.

On a remote server or HPC node, start with `--host 0.0.0.0 --port 8050`, then forward the port locally with `ssh -L 8050:localhost:8050 user@server` and open the same address.

## Sample sheet

A tab-separated file in the same format as ggsashimi: column 1 = sample id, column 2 = BAM path (relative paths are resolved against the sheet's folder), then any metadata columns. The header line is optional (without one, metadata columns are named `col3`, `col4`, …).

```
id	bam	condition	replicate
ctrl_1	bams/ctrl_1.bam	ctrl	1
ctrl_2	bams/ctrl_2.bam	ctrl	2
kd_1	bams/kd_1.bam	kd	1
kd_2	bams/kd_2.bam	kd	2
```

In the app, choose any column to group, color or label tracks. BAM files can also be given directly: `-b a.bam b.bam c.bam`. BAMs must be coordinate-sorted and indexed.

## Static figures and batch mode

`sashimi plot` renders a figure without the browser and accepts ggsashimi's main options:

```bash
sashimi plot -b samples.tsv -g genes.gtf -c chr7:128938246-128942070 \
    -O condition -C condition -A mean_j --shrink -M 5 -o fig.pdf

# Reuse the settings saved from the app, so batch figures match what you saw:
sashimi plot -b samples.tsv -g genes.gtf -c IRF5 --params sashimi.settings.json -o IRF5.pdf

# Highlight introns of interest:
sashimi plot -b samples.tsv -c chr7:128938246-128942070 --highlight "chr7:128,940,876-128,941,278" -o fig.svg
```

With many samples, add `--tracks-per-file 25` to split the figure into several readable files (`fig_part01.png`, …) that share the same scales.

`-O/-C/-L` take a column name or a ggsashimi-style 1-based column index. For stranded libraries (`-s` other than `NONE`) without `-S`, separate `_plus` and `_minus` files are written, as in ggsashimi. Run `sashimi plot -h` for all options.

## Library strandedness

| `-s` | Use for |
|---|---|
| `NONE` | Unstranded libraries (default) |
| `SENSE` / `ANTISENSE` | Single-end stranded libraries (forward / reverse, e.g. dUTP = `ANTISENSE`) |
| `MATE1_SENSE` | Paired-end, read 1 on the transcript strand (Salmon `ISF`, featureCounts `-s 1`) |
| `MATE2_SENSE` | Paired-end, read 2 on the transcript strand: Illumina TruSeq Stranded / dUTP (Salmon `ISR`, featureCounts `-s 2`) |

## Compared with ggsashimi

| | ggsashimi | SashimiPro |
|---|---|---|
| Dependencies | Python + R (ggplot2, data.table, gtable, …) | Python only |
| Reading BAMs | Per-base Python loop | NumPy difference arrays |
| Changing a parameter | Re-run the whole script | Instant redraw in the browser |
| Additional features | — | Gene-name lookup, CPM, read filters, intron highlighting, annotated vs novel junctions, CDS/UTR display, collapsed genes, adjustable shrink strength, hover details |

Notes: as in ggsashimi, only junctions with both ends inside the region are drawn. Arc widths scale with the square root of the read count, relative to the largest junction in the figure.

## Acknowledgements

SashimiPro is an independent reimplementation inspired by [ggsashimi](https://github.com/guigolab/ggsashimi) (Garrido-Martín et al., *PLoS Comput Biol* 2018), whose sample-sheet format and command-line options it follows.

## License

[MIT](LICENSE)
