# SashimiPro User Manual

[中文版](manual_zh.md)

## Contents

1. [Concepts](#1-concepts)
2. [Starting the app](#2-starting-the-app)
3. [Layout](#3-layout)
4. [Sidebar parameters](#4-sidebar-parameters)
   - [Data](#data)
   - [Samples & colors](#samples--colors)
   - [Junctions](#junctions)
   - [Highlight introns](#highlight-introns)
   - [Coverage & axes](#coverage--axes)
   - [Size & text](#size--text)
   - [Annotation](#annotation)
   - [Export](#export)
5. [Interacting with the plot](#5-interacting-with-the-plot)
6. [Recipes](#6-recipes)
7. [Troubleshooting](#7-troubleshooting)

---

## 1. Concepts

A sashimi plot draws one **track** per sample (or per group of samples). Each track has:

- **Coverage**: the filled area, i.e. how many reads cover each base.
- **Junction arcs**: arcs joining the two ends of an intron, meaning reads were spliced from one exon directly to the other. The number on an arc is the number of supporting reads; thicker arcs have more reads.
- At the bottom, the **annotation track**: transcript models from the GTF (thick boxes = CDS, thin boxes = UTR, line = intron, arrowheads = direction of transcription).

**Two kinds of parameters:**

- **Ones that re-read the BAMs**: sample files, GTF, region, strandedness and read filters. After changing them press **Load** (or Enter in the Region box).
- **Everything else** only changes how the plot is drawn. The plot updates instantly, without reloading, and keeps your current zoom.

## 2. Starting the app

```bash
sashimi app -b samples.tsv -g genes.gtf -c chr7:128,938,246-128,942,070 -s MATE2_SENSE
```

Open <http://127.0.0.1:8050> in a browser. Command-line options only pre-fill the sidebar; everything can be changed in the app.

On a server: `sashimi app ... --host 0.0.0.0 --port 8050`, then run `ssh -L 8050:localhost:8050 user@server` locally and open the same address.

## 3. Layout

- **Left sidebar**: parameters, in collapsible sections (click a section title to open/close it).
- **Right**: the plot. While data is loading, a blue progress bar and "Loading BAMs…" appear at the top and the buttons are disabled.

## 4. Sidebar parameters

### Data

| Parameter | Meaning |
|---|---|
| **Samples** | Path to a sample sheet (TSV), or one or more BAM paths separated by spaces. Sheet format below. Paths must not contain spaces. |
| **Annotation GTF** | Gene annotation (`.gtf` or `.gtf.gz`), optional. Needed for the annotation track, gene-name lookup and annotated/novel junction calls. A large GTF is slow the first time (up to about a minute), then cached in `~/.cache/sashimi_pro/`. |
| **Region or gene name** | Region such as `chr7:128,938,246-128,942,070` (commas optional, 1-based, same as IGV). You can also type a gene name or gene id (e.g. `IRF5`, `ENSG00000128604`): the whole gene plus 10% flanks on each side is loaded. Maximum 2 Mb. |
| **Library strandedness** | How reads relate to the transcript strand (table below). Decides whether each read counts towards the plus or minus strand. A wrong choice mixes up the strands, but total coverage in "both" mode is unaffected. |
| **min MAPQ** | Reads with mapping quality below this are ignored. In STAR, unique = 255 and multi-mapped = 0/1/3; use 255 (STAR) or ~10 (HISAT2 etc.) for unique reads only. Default 0 = no filter. |
| **skip duplicates** | Ignore reads flagged as PCR duplicates (the BAM must have been marked with Picard / samtools markdup). Off by default. |
| **skip secondary** | Ignore secondary alignments (non-primary hits of multi-mapped reads) so a read is not counted twice. On by default. |

**Choosing the library strandedness:**

| Option | Meaning | Typical case |
|---|---|---|
| `NONE` | Unstranded | When unsure |
| `SENSE` | Every read is on the transcript's strand | Single-end, forward-stranded |
| `ANTISENSE` | Every read is on the opposite strand | Single-end dUTP |
| `MATE1_SENSE` | read1 sense, read2 antisense | Paired-end; Salmon `ISF`, featureCounts `-s 1` |
| `MATE2_SENSE` | read2 sense, read1 antisense | **Paired-end Illumina TruSeq Stranded / dUTP (most common)**; Salmon `ISR`, featureCounts `-s 2` |

**Sample sheet (TSV)**: column 1 = sample id, column 2 = BAM path (relative paths are resolved against the TSV's folder), then any number of metadata columns (condition, batch, replicate…) that can be used for grouping, colouring and labels. The first line may be a header:

```
id	bam	condition	replicate
ctrl_1	bams/ctrl_1.bam	ctrl	1
ctrl_2	bams/ctrl_2.bam	ctrl	2
kd_1	bams/kd_1.bam	kd	1
```

Without a header, metadata columns are named `col3`, `col4`, …

**Buttons:**

| Button | Action |
|---|---|
| **Load** | Read the BAMs with the current settings. |
| **◀ / ▶** | Move half a window left/right and reload. |
| **zoom out** | Make the region 3× wider (one window width added on each side) and reload. |
| **zoom in** | First drag a box on the plot, then click: reloads only the selected range. |

After loading, the text under the buttons shows the number of samples, the region and the number of transcripts. Problems (missing file, a sample with no reads in the region) show as ⚠ warnings; errors show in red.

### Samples & colors

| Parameter | Meaning |
|---|---|
| **Samples shown** | Which samples to show. Click × to remove one; pick from the dropdown to add it back. Track order follows the sample sheet. |
| **Group (overlay) by column** | Merge samples sharing a value in this column into one track; e.g. with `condition`, all ctrl samples form one track and all kd samples another. Empty = one track per sample. |
| **Within a group** | How samples in a group are combined: `overlay` = each sample drawn semi-transparent with a thin outline, on top of each other, so you can see variation within the group (when a group shares one colour, samples get lighter/darker shades of it; near-identical samples will still look like a single layer); `mean`; `median`; `sum`. |
| **Color by column** | Column that decides colours. Each distinct value gets a colour, in order of first appearance in the sheet. Empty = everything grey-blue, and the Palette has no effect. |
| **Label column (ungrouped)** | Column shown as the track name (top-left) when not grouping; default `id`. When grouping, the label is the group name and size, e.g. `ctrl (n=2)`. |
| **Palette** | The list of colours handed out to the values of the Color column (only used when Color by column is set). `Okabe-Ito` (default, colour-blind safe, good for papers), `Set1` (bright), `Dark2` (muted), `D3`, `Tableau`, `ggsashimi` (same as the old tool: red, green, blue, black). |
| **Custom colors** | Your own colours, comma separated, e.g. `#d62728, steelblue, grey`; hex codes or colour names. Overrides the Palette, in the same value order. |

**Tip**: setting Group and Color to the same column (e.g. both `condition`) is the usual choice. In `overlay` mode, a Color column different from the Group column (e.g. Group = `condition`, Color = `replicate`) gives each sample within a group its own colour.

### Junctions

| Parameter | Meaning |
|---|---|
| **Minimum reads to draw** | Junctions with fewer reads are not drawn; use it to hide low-level noise. It is compared with the **displayed value**: after group aggregation, and after CPM scaling if enabled. |
| **Aggregate counts in a group** | How junction counts are combined within a group: `mean`, `median`, `sum`, or `list each` (show every sample's count, e.g. `355,346`; width and filtering use the mean). No effect without grouping. |
| **Arc placement** | `alternate` = arcs alternate above/below in genomic order (default, reduces overlap); `above` = all above; `below` = all below. Arcs above start at the coverage height of the flanking exons; arcs below start at 0. |
| **Max arc width** | Width in pixels of the thickest arc. The junction with most reads in the whole figure gets this width; others scale with the square root of their read count. |
| **show counts** | Print the read count on each arc. |
| **only annotated junctions** | Draw only junctions present in the GTF (needs a GTF). |
| **dotted line for novel junctions** | Draw junctions absent from the GTF with a dotted line, to spot new splicing events (needs a GTF). |

Only junctions with both ends inside the loaded region are drawn.

### Highlight introns

Emphasise the arcs of one or more introns of interest in every track.

| Parameter | Meaning |
|---|---|
| **Introns to highlight** | Intron coordinates such as `chr7:128,940,876-128,941,278`; separate several with `;`. Coordinates are the first and last intron base (1-based, as in STAR `SJ.out.tab` and the hover text); the chromosome can be omitted. **Easier: click an arc on the plot** to add its intron to the box; click again to remove it. |
| **color** | Highlight colour, default red `#d62728`; hex or colour name. |
| **match ± bp** | Matching tolerance. Tools sometimes differ by 1 in coordinate conventions (e.g. BED is 0-based); 1 covers that. Default 0 = exact. |
| **fade the other arcs** | Make the other arcs pale so the highlighted ones stand out. |
| **shade the intron across tracks** | Draw a light background band over the intron in all tracks. |
| **Clear** | Remove all highlights. |

Highlighted arcs are thicker, in the highlight colour, with bold counts. **Highlighted introns ignore the Minimum reads and only-annotated filters**, so they are drawn even with few reads; a track has no arc only if that sample has no such junction at all.

Command line: `sashimi plot ... --highlight "chr7:128,940,876-128,941,278"`.

### Coverage & axes

| Parameter | Meaning |
|---|---|
| **Strand** | Which strand's signal to show: `both` = combined, `plus`, `minus`. Only meaningful for stranded libraries (strandedness ≠ NONE). For a gene, choosing its own strand removes antisense transcription. |
| **Scale** | `raw reads`, or `per million mapped` = CPM, divided by the number of mapped reads in that BAM (in millions). Use CPM to compare samples with different sequencing depth, and with `mean` grouping. |
| **same y scale for all tracks** | Use one y-axis range for every track so heights are directly comparable. Otherwise each track is scaled to its own maximum. |
| **shrink introns** | Draw introns shorter so exons are easier to see (introns are usually much longer). Axis labels stay in real coordinates, but ticks are no longer evenly spaced. |
| **Shrink strength** | An intron of length L is drawn with length L^x. Smaller x = stronger shrinking; 1 = none; default 0.7. |
| **Fill opacity** | Opacity of the coverage fill (0.1 faint – 1 solid). Lower it in `overlay` mode to see overlapping samples. |

### Size & text

| Parameter | Meaning |
|---|---|
| **Track height** | Height of each sample track, in pixels. |
| **Width** | Figure width in pixels. 0 = fit the browser window; exports then use 1000 px. Set a fixed width for figures so every export is identical. |
| **Font size** | Base font size; labels, axes and junction counts scale with it. |
| **show title** | Show the title (region, strand, CPM, shrinking). |

### Annotation

Requires a GTF.

| Parameter | Meaning |
|---|---|
| **Show** | `transcripts` = one row per transcript; `collapsed genes` = one row per gene (union of all its exons), handy for genes with many isoforms; `none` = hide the annotation. |
| **Max transcripts** | In `transcripts` mode, show at most this many (sorted by gene and start position, first N kept). |
| **Pick transcripts** | Choose exactly which transcripts to show (type to search). Overrides Max transcripts. Clear it to go back to the first N. |

Hovering an exon shows the transcript name, gene, transcript id, biotype and strand.

### Export

| Parameter / button | Action |
|---|---|
| **File name** | Name of the downloaded file, without extension. |
| **Format dropdown** | `pdf` / `svg` (vector, editable in Illustrator/Inkscape; recommended for papers) / `png`. |
| **PNG resolution ×** | PNG pixel multiplier. The figure is Width × total height in pixels; × 2 doubles both (e.g. 1000 px wide becomes 2000 px), so it is sharper and larger but the layout is unchanged. PNG only; PDF/SVG are vector and ignore it. Default 2, range 1–6. PNGs are limited to about 60,000 px per side; beyond that the multiplier is lowered automatically with a notice (otherwise the file would be corrupt). |
| **Tracks per file** | Maximum tracks per file, default 25. With more tracks, Figure downloads a ZIP of `name_part01`, `_part02`, … (e.g. 205 samples → 9 images). All files share the same x range, arc-width scale and (with same y scale) y axis, so they can be compared directly. 0 = everything in one file. |
| **Figure** | Export the whole figure with the current parameters (the full loaded region, not the current zoom). |
| **Junctions** | Export raw read counts of every junction in every sample (TSV: chromosome, intron start/end, sample, count, strand; plus an `annotated` column with a GTF). Ignores display filters such as Minimum reads, so it is suitable for statistics. |
| **Settings** | Save all plot parameters to a JSON file. |
| **Drop / click to load a settings JSON** | Load a saved settings file and restore all plot parameters (region and data are not changed). |

A saved settings JSON also works with the command line, so batch figures match what you saw in the app:

```bash
sashimi plot -b samples.tsv -g genes.gtf --params sashimi.settings.json -o fig.pdf
```

## 5. Interacting with the plot

- **Hover**: on coverage, position and value; on arcs, intron coordinates, length, read count and whether the junction is annotated.
- **Drag a box** to zoom (display only, no reloading); **double-click** to reset.
- After zooming, click **zoom in** in the sidebar to reload just that range (finer detail, junctions re-filtered).
- Changing parameters keeps the zoom; changing the region or toggling shrink resets it.
- Camera icon in the plot's top-right toolbar: download the **current view** as SVG.

## 6. Recipes

**Compare splicing between two conditions (e.g. exon skipping)**
1. Group by = `condition`, Color by = `condition`.
2. Within a group = `mean`, Scale = `per million mapped`.
3. Aggregate counts = `mean` (or `list each` to see each replicate).
4. Tick same y scale and shrink introns.

**Check replicate consistency**: Group by = `condition`, Within a group = `overlay`, Fill opacity 0.3–0.4, Aggregate counts = `list each`.

**Follow one intron across groups**: click its arc (or type it under Highlight introns); other arcs fade. Combine with Group by and same y scale to compare read counts between groups.

**Find novel splicing**: provide a GTF, tick dotted line for novel junctions, hover the dotted arcs for coordinates.

**Too many arcs**: raise Minimum reads, set Arc placement to `above`, or untick show counts.

**Publication figure**: fixed Width (e.g. 1000), Font size 12–14, export SVG/PDF. Also click Settings to save the parameters and reuse them for other genes.

## 7. Troubleshooting

**The plot is empty / only shows a message.** Make sure Samples and Region are filled in, press Load, and read the status text under the buttons.

**A sample has no signal.** The status shows a "no reads in region" warning. Check the chromosome name (`chr7` and `7` are tried automatically) and the genome build (hg19 vs hg38).

**Plus and minus strands look swapped.** The strandedness is wrong; paired-end dUTP libraries need `MATE2_SENSE`. Press Load after changing it.

**Changing the Palette does nothing.** Set Color by column first.

**Gene name not found.** A GTF is required; names are case-insensitive, and Ensembl gene ids work too.

**With hundreds of samples, the exported PNG is unreadable or won't open.** Hundreds of tracks make a very tall image (205 tracks ≈ 24,000 px), which image viewers shrink to fit the window, so it looks blurry; a large multiplier can also exceed the browser's pixel limit and corrupt the file. Use **Tracks per file** to split the export (25 tracks per image by default), or group / select fewer samples. On the command line: `--tracks-per-file 25`.
