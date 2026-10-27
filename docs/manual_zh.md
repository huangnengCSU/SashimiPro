# SashimiPro 用户手册

[English version](manual_en.md)

## 目录

1. [基本概念](#1-基本概念)
2. [启动](#2-启动)
3. [界面总览](#3-界面总览)
4. [左侧参数逐项说明](#4-左侧参数逐项说明)
   - [Data（数据）](#data数据)
   - [Samples & colors（样本与颜色）](#samples--colors样本与颜色)
   - [Junctions（剪接位点）](#junctions剪接位点)
   - [Highlight introns（高亮内含子）](#highlight-introns高亮内含子)
   - [Coverage & axes（覆盖度与坐标轴）](#coverage--axes覆盖度与坐标轴)
   - [Size & text（尺寸与文字）](#size--text尺寸与文字)
   - [Annotation（注释）](#annotation注释)
   - [Export（导出）](#export导出)
5. [图上的交互操作](#5-图上的交互操作)
6. [常见任务](#6-常见任务)
7. [常见问题](#7-常见问题)

---

## 1. 基本概念

sashimi plot 每个样本（或每组样本）画一个**轨道（track）**，每个轨道包括：

- **覆盖度（coverage）**：填充的面积图，表示每个碱基上有多少 reads 覆盖。
- **剪接弧线（junction arc）**：连接内含子两端的弧线，表示有 reads 从一个外显子直接跳到另一个外显子。弧线上的数字是支持这个剪接的 reads 数，线越粗 reads 越多。
- 最下方是 **注释轨道**：GTF 里的转录本结构（粗块 = CDS，细块 = UTR，线 = 内含子，箭头 = 转录方向）。

**两类参数：**

- **需要重新读 BAM 的**：样本文件、GTF、区域、链特异性、read 过滤。改完要按 **Load**（或在 Region 框按回车）。
- **其他所有参数**：只改变画法，修改后图会立刻更新，不用重新加载，而且会保持当前的缩放。

## 2. 启动

```bash
sashimi app -b samples.tsv -g genes.gtf -c chr7:128,938,246-128,942,070 -s MATE2_SENSE
```

浏览器打开 <http://127.0.0.1:8050>。命令行参数只是预先填好左侧的输入框，启动后都可以在界面里修改。

在服务器上运行：`sashimi app ... --host 0.0.0.0 --port 8050`，本地执行 `ssh -L 8050:localhost:8050 用户名@服务器` 后打开同样的地址。

## 3. 界面总览

- **左侧栏**：参数，分成几个可以折叠的分组（点分组标题展开/收起）。
- **右侧**：图。加载数据时顶部会出现蓝色进度条和 “Loading BAMs…” 提示，此时按钮暂时不可用。

## 4. 左侧参数逐项说明

### Data（数据）

| 参数 | 含义 |
|---|---|
| **Samples** | 样本表（TSV）路径，或者一个或多个 BAM 文件路径（空格分隔）。样本表格式见下文。路径里不能有空格。 |
| **Annotation GTF** | 基因注释文件（`.gtf` 或 `.gtf.gz`），可不填。填了才有注释轨道、才能用基因名定位、才能区分已注释/新的 junction。大文件第一次读取较慢（可能一分钟左右），之后有缓存（`~/.cache/sashimi_pro/`）。 |
| **Region or gene name** | 要看的区域，格式 `chr7:128,938,246-128,942,070`（逗号可有可无，1-based 坐标，与 IGV 一致）。也可以直接写基因名或 gene id（如 `IRF5`、`ENSG00000128604`），会自动取整个基因并两边各加 10% 的侧翼。最大 2 Mb。 |
| **Library strandedness** | 文库的链特异性，见下表。决定每条 read 被算到正链还是负链。选错了链分不开，但 “both” 模式下的总覆盖度不受影响。 |
| **min MAPQ** | 最低比对质量。比对质量低于此值的 read 不计入。STAR 中唯一比对 = 255，多重比对 = 0/1/3；想只看唯一比对可设为 255（STAR）或 10 左右（HISAT2 等）。默认 0 = 不过滤。 |
| **skip duplicates** | 忽略被标记为 PCR 重复的 reads（需要 BAM 已经用 Picard/samtools markdup 标记过）。默认关闭。 |
| **skip secondary** | 忽略 secondary alignment（多重比对中的非主比对），避免同一条 read 被重复计数。默认开启。 |

**Library strandedness 怎么选：**

| 选项 | 含义 | 常见情况 |
|---|---|---|
| `NONE` | 非链特异性，不区分正负链 | 不确定时选这个 |
| `SENSE` | 所有 reads 与转录本同向 | 单端、正向链特异性文库 |
| `ANTISENSE` | 所有 reads 与转录本反向 | 单端 dUTP 文库 |
| `MATE1_SENSE` | read1 与转录本同向、read2 反向 | 双端，Salmon `ISF`，featureCounts `-s 1` |
| `MATE2_SENSE` | read2 与转录本同向、read1 反向 | **双端 Illumina TruSeq Stranded / dUTP（最常见）**，Salmon `ISR`，featureCounts `-s 2` |

**样本表（TSV）格式**：第 1 列样本 id，第 2 列 BAM 路径（相对路径以 TSV 所在目录为准），之后可以有任意多列元数据（条件、批次、重复号……），用来分组、上色、做标签。第一行可以是表头：

```
id	bam	condition	replicate
ctrl_1	bams/ctrl_1.bam	ctrl	1
ctrl_2	bams/ctrl_2.bam	ctrl	2
kd_1	bams/kd_1.bam	kd	1
```

没有表头时，元数据列自动命名为 `col3`、`col4`…

**按钮：**

| 按钮 | 作用 |
|---|---|
| **Load** | 按当前设置读取 BAM。 |
| **◀ / ▶** | 向左/右移动半个窗口并重新加载。 |
| **zoom out** | 把范围扩大到原来的 3 倍（两边各扩一个窗口宽度）并重新加载。 |
| **zoom in** | 先在图上拖动鼠标框选一段，再点它：只重新加载框选的区间。 |

加载完成后，按钮下面会显示样本数、区域、转录本数；有问题（文件找不到、某样本在该区域没有 reads）会显示 ⚠ 警告，出错显示红字。

### Samples & colors（样本与颜色）

| 参数 | 含义 |
|---|---|
| **Samples shown** | 显示哪些样本。点 × 去掉某个样本，点下拉框加回来。轨道顺序 = 样本表里的顺序。 |
| **Group (overlay) by column** | 按样本表的某一列分组，同一组的样本合并到一个轨道。例如选 `condition`，所有 ctrl 合成一个轨道、所有 kd 合成一个轨道。不选 = 每个样本一个轨道。 |
| **Within a group** | 一组内多个样本怎么合并：`overlay` = 每个样本半透明叠在一起并带细轮廓线，能看到组内差异（同组颜色相同时自动用该颜色的深浅不同色调区分样本；样本数据很接近时重叠部分看起来几乎是一层，这是正常的）；`mean` = 平均值；`median` = 中位数；`sum` = 相加。 |
| **Color by column** | 按哪一列上色。这一列的每个不同取值一种颜色（按在样本表里出现的顺序分配）。不选 = 全部灰蓝色，下面的 Palette 不起作用。 |
| **Label column (ungrouped)** | 不分组时，轨道左上角显示哪一列作为名字（默认 `id`）。分组时显示的是组名和样本数，如 `ctrl (n=2)`。 |
| **Palette** | 调色板，即上色时依次使用的一组颜色（只在设置了 Color by column 时生效）。`Okabe-Ito`（默认，色盲友好，适合论文）、`Set1`（鲜艳）、`Dark2`（偏暗沉稳）、`D3`、`Tableau`、`ggsashimi`（与旧版一致：红绿蓝黑）。 |
| **Custom colors** | 自定义颜色，逗号分隔，如 `#d62728, steelblue, grey`，十六进制或颜色名都可以。填了以后覆盖 Palette，按取值出现顺序使用。 |

**提示**：Group 和 Color 选同一列（例如都选 `condition`）是最常见的用法。在 `overlay` 模式下，如果 Color by 选了和 Group 不同的列（比如 Group = `condition`，Color = `replicate`），同一组内的不同样本会用不同颜色叠加。

### Junctions（剪接位点）

| 参数 | 含义 |
|---|---|
| **Minimum reads to draw** | reads 数低于此值的 junction 不画。用来去掉低丰度噪音。注意比较的是**最终显示的数值**：分组后是聚合后的值，开了 CPM 是标准化后的值。 |
| **Aggregate counts in a group** | 分组时 junction 的 reads 数怎么合并：`mean`、`median`、`sum`，或 `list each`（把每个样本的数都列出来，如 `355,346`；线宽和过滤按平均值）。不分组时不起作用。 |
| **Arc placement** | 弧线画在哪里：`alternate` = 按位置顺序上下交替（默认，避免重叠），`above` = 全在上方，`below` = 全在下方。上方弧线从两端外显子的覆盖度高度起画；下方弧线从 0 起画。 |
| **Max arc width** | 最粗弧线的线宽（像素）。全图 reads 最多的 junction 用这个线宽，其他按 reads 数的平方根成比例变细。 |
| **show counts** | 在弧线上显示 reads 数。 |
| **only annotated junctions** | 只画 GTF 里已有的剪接（需要 GTF）。 |
| **dotted line for novel junctions** | GTF 里没有的“新”剪接用点线画，方便发现新的剪接事件（需要 GTF）。 |

只画两端都在当前区域内的 junction；跨出区域边界的不画。

### Highlight introns（高亮内含子）

对某个 intron 特别感兴趣时，可以把它的弧线在所有轨道里突出显示。

| 参数 | 含义 |
|---|---|
| **Introns to highlight** | 要高亮的内含子坐标，如 `chr7:128,940,876-128,941,278`；多个用 `;` 分隔。坐标是内含子的第一个和最后一个碱基（1-based，与 STAR `SJ.out.tab` 和悬停信息中的一致），染色体名可省略。**更简单的方法：直接在图上点击弧线**，就会把该内含子加入输入框；再点一次取消。 |
| **color** | 高亮颜色，默认红色 `#d62728`，可填十六进制或颜色名。 |
| **match ± bp** | 匹配容差。不同工具的坐标约定可能差 1（例如 BED 是 0-based），设为 1 可以兼容。默认 0 = 精确匹配。 |
| **fade the other arcs** | 把其他弧线变淡，让高亮的更突出。 |
| **shade the intron across tracks** | 在所有轨道中给该内含子区域加一条浅色背景带。 |
| **Clear** | 清空所有高亮。 |

高亮的弧线更粗、数字加粗并用高亮颜色显示。**高亮的 intron 不受 Minimum reads 和 only annotated 的过滤**，即使 reads 很少也会画出来；如果某个样本里根本没有这个剪接，该轨道上就没有弧线。

命令行出图：`sashimi plot ... --highlight "chr7:128,940,876-128,941,278"`。

### Coverage & axes（覆盖度与坐标轴）

| 参数 | 含义 |
|---|---|
| **Strand** | 显示哪条链的信号：`both` = 两条链合在一起，`plus` = 正链，`minus` = 负链。只对链特异性文库有意义（Library strandedness ≠ NONE）；看基因时选与基因同向的链可以去掉反义链转录的干扰。 |
| **Scale** | `raw reads` = 原始 read 数；`per million mapped` = CPM，除以该 BAM 全部比对上的 reads 数（百万）。测序深度不同的样本相互比较、或做 `mean` 分组时建议用 CPM。 |
| **same y scale for all tracks** | 所有轨道使用相同的 y 轴范围，方便直接比较高度。不勾选时每个轨道按自己的最大值缩放。 |
| **shrink introns** | 压缩内含子显示长度，让外显子更显眼（内含子通常比外显子长得多）。坐标轴标签仍是真实坐标，但刻度不再等距。 |
| **Shrink strength** | 压缩强度：长度为 L 的内含子显示为 L^x。x 越小压缩越狠；1 = 不压缩；默认 0.7。 |
| **Fill opacity** | 覆盖度填充的不透明度（0.1 很淡 – 1 不透明）。`overlay` 模式下调低可以看清重叠部分。 |

### Size & text（尺寸与文字）

| 参数 | 含义 |
|---|---|
| **Track height** | 每个样本轨道的高度（像素）。 |
| **Width** | 图的宽度（像素）。0 = 自动适应浏览器窗口；导出时如果是 0 则按 1000 像素。投稿时可以固定一个宽度，保证每次导出一致。 |
| **Font size** | 基础字号，标签、坐标轴、junction 数字都按它缩放。 |
| **show title** | 是否显示顶部标题（区域、链、是否 CPM、是否压缩）。 |

### Annotation（注释）

需要填写 GTF。

| 参数 | 含义 |
|---|---|
| **Show** | `transcripts` = 每个转录本一行；`collapsed genes` = 每个基因合并成一行（所有转录本的外显子取并集），适合转录本很多的基因；`none` = 不显示注释。 |
| **Max transcripts** | `transcripts` 模式下最多显示几个转录本（按基因、起始位置排序后取前 N 个）。 |
| **Pick transcripts** | 手动挑选要显示的转录本，可以输入名字搜索。选了以后忽略 Max transcripts。清空 = 回到前 N 个。 |

鼠标悬停在外显子上会显示转录本名、基因名、转录本 id、biotype 和链。

### Export（导出）

| 参数 / 按钮 | 作用 |
|---|---|
| **File name** | 下载文件的名字（不含扩展名）。 |
| **格式下拉框** | `pdf` / `svg`（矢量图，可在 Illustrator/Inkscape 里编辑，投稿推荐）/ `png`。 |
| **PNG resolution ×** | PNG 的分辨率倍数。图的大小按 Width × 总高度（像素）算，× 2 就是长宽各乘 2（如 1000 像素宽变成 2000 像素），图更清晰、文件更大；排版尺寸不变。只对 PNG 有效，PDF/SVG 是矢量图，忽略此项。默认 2，范围 1–6。PNG 每边最多约 60,000 像素，超过时会自动降低倍数并提示（否则文件会损坏）。 |
| **Tracks per file** | 每个文件最多放几个轨道，默认 25。轨道数超过它时，Figure 会下载一个 ZIP，里面是 `名字_part01`、`_part02`… 多个文件（例如 205 个样本 → 9 张图）。各文件的 x 范围、弧线粗细标准和（勾选 same y scale 时的）y 轴完全一致，可以直接对比。0 = 全部放在一个文件里。 |
| **Figure** | 按当前参数导出整张图（整个加载的区域，不是当前缩放的视野）。 |
| **Junctions** | 导出当前区域所有样本、所有 junction 的原始 reads 数（TSV：染色体、内含子起止、样本、reads 数、链；有 GTF 时还有 annotated 列）。不受 Minimum reads 等显示参数影响，适合做后续统计。 |
| **Settings** | 把当前所有画图参数保存为 JSON 文件。 |
| **Drop / click to load a settings JSON** | 读入之前保存的参数文件，恢复所有画图参数（不会改变区域和数据）。 |

保存的 settings JSON 也可以给命令行批量出图使用，保证和界面中看到的完全一致：

```bash
sashimi plot -b samples.tsv -g genes.gtf --params sashimi.settings.json -o fig.pdf
```

## 5. 图上的交互操作

- **悬停**：在覆盖度上显示位置和数值；在弧线上显示内含子坐标、长度、reads 数，以及是否为注释中的剪接。
- **拖动框选**：放大（只放大显示，不重新读数据）；**双击**：恢复。
- 放大后点左侧 **zoom in**：按放大的范围重新加载（更精细，junction 也重新筛选）。
- 改变参数时会保持当前缩放；换区域或开关 shrink 时会重置。
- 图右上角工具栏的相机图标：把**当前视野**下载为 SVG。

## 6. 常见任务

**比较两组条件的可变剪接（如 exon skipping）**
1. Group by = `condition`，Color by = `condition`。
2. Within a group = `mean`，Scale = `per million mapped`。
3. Aggregate counts = `mean`（或 `list each` 看每个重复）。
4. 勾选 same y scale、shrink introns。

**看每个重复的一致性**：Group by = `condition`，Within a group = `overlay`，Fill opacity 调到 0.3–0.4，Aggregate counts = `list each`。

**盯住某个 intron 比较各组**：在图上点击它的弧线（或在 Highlight introns 里输入坐标），其他弧线会变淡；配合 Group by 和 same y scale 比较各组的 reads 数。

**找新的剪接事件**：填 GTF，勾选 dotted line for novel junctions，悬停在点线上看坐标。

**弧线太多太乱**：调高 Minimum reads；或 Arc placement 换成 `above`；或取消 show counts。

**出投稿图**：Width 设为固定值（如 1000），Font size 12–14，导出 SVG/PDF；同时点 Settings 保存参数，方便日后用同样的参数出其他基因的图。

## 7. 常见问题

**图是空的 / 只有一行提示？** 确认 Samples 和 Region 都已填写，然后点 Load。看按钮下方的状态文字有没有报错。

**某个样本没有信号？** 状态中会有 “no reads in region” 警告。检查染色体名（`chr7` 与 `7` 会自动互相尝试）和坐标版本（hg19/hg38）。

**plus / minus 链信号看起来反了？** Library strandedness 选错了，双端 dUTP 文库应选 `MATE2_SENSE`；改完要按 Load。

**修改 Palette 没反应？** 需要先设置 Color by column。

**基因名找不到？** 需要填写 GTF；名字不区分大小写，也可以用 Ensembl gene id。

**样本很多（上百个）时导出的 PNG 看不清 / 打不开？** 一张图里放几百个轨道会非常高（205 个轨道约 24,000 像素），看图软件会把整张图缩到窗口大小，所以看起来模糊；倍数太大还会超过浏览器的像素上限导致文件损坏。用 **Tracks per file** 分成多个文件（默认 25 个轨道一张），或者先分组 / 只选部分样本。命令行用 `--tracks-per-file 25`。
