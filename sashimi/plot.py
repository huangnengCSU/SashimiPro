"""Build a sashimi plot (Plotly figure) from loaded RegionData + display params.

Everything here is cheap compared to reading the BAMs, so the app calls
build_figure() on every parameter change.
"""

import math
import re

import numpy as np
import plotly.graph_objects as go
from plotly.colors import qualitative
from plotly.subplots import make_subplots

from .annotation import annotated_introns, collapse_genes

PALETTES = {
    "Okabe-Ito": ["#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7", "#F0E442", "#000000"],
    "Set1": qualitative.Set1,
    "Dark2": qualitative.Dark2,
    "D3": qualitative.D3,
    "Tableau": qualitative.T10,
    "ggsashimi": ["#ff0000", "#00ff00", "#0000ff", "#000000"],
}
SINGLE_COLOR = "#5b6b7f"
ANN_COLOR = "#2f3b4c"

DEFAULTS = dict(
    # which signal
    strand="combined",          # '+', '-', 'combined'
    samples=None,               # list of sample ids to show (None = all)
    normalize="raw",            # 'raw' or 'cpm' (per million mapped reads)
    # grouping / colors / labels (column names of the sample sheet)
    group_by=None,
    group_mode="overlay",       # 'overlay', 'mean', 'median', 'sum'
    color_by=None,
    label_by="id",
    palette="Okabe-Ito",
    custom_colors="",           # comma separated; overrides palette
    # junctions
    min_junction=1.0,
    junction_agg="mean",        # 'mean', 'median', 'sum', 'each'
    annotated_only=False,
    novel_dashed=False,
    show_labels=True,
    arc_layout="alternate",     # 'alternate', 'above', 'below'
    arc_width=4.0,
    # highlighted introns
    highlight="",               # "chr7:128,940,001-128,940,500; 128941001-128941200"
    highlight_color="#d62728",
    highlight_tol=0,            # match introns within +-N bp (coordinate conventions differ by 1 between tools)
    dim_others=True,
    highlight_band=True,
    # axes
    fix_y=False,
    shrink=False,
    shrink_power=0.7,
    # look
    alpha=0.6,
    track_height=110,
    width=1000,
    font_size=12,
    show_title=True,
    # annotation track
    ann_mode="transcripts",     # 'transcripts', 'genes', 'none'
    max_transcripts=15,
    transcripts=None,           # explicit list of transcript ids to show
)

MAX_POINTS = 4000               # per coverage curve


def merged_params(params):
    p = dict(DEFAULTS)
    p.update({k: v for k, v in (params or {}).items() if k in DEFAULTS})
    return p


# --------------------------------------------------------------------------
# Coordinate shrinking
# --------------------------------------------------------------------------

def intersect_introns(introns):
    """Collapse overlapping introns to their common (intersected) part."""
    out = []
    for a, b in sorted(introns):
        if out and a <= out[-1][1]:
            out[-1] = (max(out[-1][0], a), min(out[-1][1], b))
        else:
            out.append((a, b))
    return out


class CoordMap:
    """Piecewise-linear map real -> display coordinates; each shrunk interval of
    length L is drawn with length L**power.  Identity when no intervals."""

    def __init__(self, start, end, intervals=(), power=0.7):
        real, disp = [start - 0.5], [start - 0.5]
        for a, b in intervals:
            a, b = a - 0.5, b + 0.5
            if a <= real[-1] or b >= end + 0.5:
                continue
            disp.append(disp[-1] + (a - real[-1]))
            real.append(a)
            disp.append(disp[-1] + (b - a) ** power)
            real.append(b)
        disp.append(disp[-1] + (end + 0.5 - real[-1]))
        real.append(end + 0.5)
        self.real, self.disp = np.array(real), np.array(disp)
        self.identity = len(real) == 2

    def fwd(self, x):
        return x if self.identity else np.interp(x, self.real, self.disp)

    def inv(self, x):
        return x if self.identity else np.interp(x, self.disp, self.real)


def parse_introns(text):
    """'chr7:1,001-2,000; 3001-4000' -> [(1001, 2000), (3001, 4000)].  Chromosome names are ignored."""
    out = []
    for m in re.finditer(r"(\d[\d,]*)\s*-\s*(\d[\d,]*)", text or ""):
        a, b = (int(x.strip(",").replace(",", "")) for x in m.groups())
        out.append((min(a, b), max(a, b)))
    return out


def format_introns(chrom, introns):
    prefix = f"{chrom}:" if chrom else ""
    return "; ".join(f"{prefix}{a:,}-{b:,}" for a, b in introns)


def nice_step(span, n):
    raw = span / max(n, 1)
    mag = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1
    for m in (1, 2, 2.5, 5, 10):
        if raw <= m * mag:
            return m * mag
    return 10 * mag


def nice_ticks(vmax, n=3):
    if vmax <= 0:
        return [0]
    step = nice_step(vmax, n)
    return list(np.arange(0, vmax + step * 0.01, step))


def fmt_count(v, normalized):
    if not normalized or v >= 100:
        return f"{v:,.0f}"
    return f"{v:.3g}"


# --------------------------------------------------------------------------
# Tracks
# --------------------------------------------------------------------------

def _palette(p):
    custom = [c.strip() for c in (p["custom_colors"] or "").replace(";", ",").split(",") if c.strip()]
    return custom or list(PALETTES.get(p["palette"], PALETTES["Okabe-Ito"]))


def _to_rgb(color):
    """'#rrggbb' / 'rgb(r,g,b)' / colour name -> (r, g, b) in 0-255, or None if unknown."""
    c = color.strip()
    if c.startswith("#") and len(c) in (4, 7):
        h = c[1:] if len(c) == 7 else "".join(ch * 2 for ch in c[1:])
        return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))
    if c.startswith("rgb"):
        return tuple(int(float(v)) for v in c[c.index("(") + 1:c.index(")")].split(",")[:3])
    try:
        from matplotlib.colors import to_rgb
        return tuple(int(round(v * 255)) for v in to_rgb(c))
    except Exception:
        return None


def _shades(color, n):
    """n distinguishable shades of one colour (darker -> lighter), for overlaid replicates."""
    rgb = _to_rgb(color)
    if n <= 1 or rgb is None:
        return [color] * n
    out = []
    for f in np.linspace(-0.35, 0.45, n):     # <0 mixes towards black, >0 towards white
        target = 0 if f < 0 else 255
        mixed = [round(v + (target - v) * abs(f)) for v in rgb]
        out.append("#%02x%02x%02x" % tuple(mixed))
    return out


def _sample_signal(sig, strand):
    # Unstranded data only has '+', so any strand request falls back to all reads.
    if strand == "combined" or strand not in sig.coverage:
        cov =sum(v.astype(np.float64) for v in sig.coverage.values())
        jd = {}
        for d in sig.junctions.values():
            for k, c in d.items():
                jd[k] = jd.get(k, 0) + c
        return cov, jd
    return sig.coverage[strand].astype(np.float64), dict(sig.junctions[strand])


def build_tracks(region, p, known_introns=None):
    """Turn per-sample signal into the list of drawn tracks."""
    sheet = region.samples.set_index("id", drop=False)
    ids = [s for s in region.signals if p["samples"] is None or s in p["samples"]]
    normalized = p["normalize"] == "cpm"
    palette = _palette(p)
    color_levels = {}

    def color_of(level):
        if level not in color_levels:
            color_levels[level] = palette[len(color_levels) % len(palette)]
        return color_levels[level]

    def col(sid, name):
        if not name or name not in sheet.columns:
            return None
        return str(sheet.at[sid, name])

    wanted = parse_introns(p["highlight"])
    tol = int(p["highlight_tol"] or 0)

    def highlighted(k):
        return any(abs(k[0] - a) <= tol and abs(k[1] - b) <= tol for a, b in wanted)

    per_sample = {}
    for sid in ids:
        sig = region.signals[sid]
        cov, jd = _sample_signal(sig, p["strand"])
        f = 1e6 / sig.total_mapped if normalized and sig.total_mapped else 1.0
        per_sample[sid] = (cov * f, {k: v * f for k, v in jd.items()})

    groups = {}
    for sid in ids:
        g = col(sid, p["group_by"]) if p["group_by"] else sid
        groups.setdefault(g, []).append(sid)

    tracks = []
    for g, members in groups.items():
        levels = {col(s, p["color_by"]) for s in members} if p["color_by"] else {None}
        if not p["color_by"]:
            tcolor = SINGLE_COLOR
        else:
            tcolor = color_of(levels.pop() if len(levels) == 1 else g)
        if p["group_by"]:
            label = f"{g} (n={len(members)})" if len(members) > 1 else g
        else:
            label = col(g, p["label_by"]) or g

        covs = np.vstack([per_sample[s][0] for s in members])
        mode = p["group_mode"] if len(members) > 1 else "mean"
        if mode == "overlay":
            colors = [color_of(col(s, p["color_by"])) if p["color_by"] else tcolor for s in members]
            if len(set(colors)) == 1:       # same colour for the whole group: use shades to tell samples apart
                colors = _shades(colors[0], len(members))
            curves = [(s, per_sample[s][0], c) for s, c in zip(members, colors)]
            disp = covs.max(axis=0)
        else:
            agg = {"mean": np.mean, "median": np.median, "sum": np.sum}[mode]
            disp = agg(covs, axis=0)
            curves = [(g, disp, tcolor)]

        keys = sorted({k for s in members for k in per_sample[s][1]})
        junctions = {}
        for k in keys:
            vals = [per_sample[s][1].get(k, 0.0) for s in members]
            if p["junction_agg"] == "each" and len(vals) > 1:
                v, text = float(np.mean(vals)), ",".join(fmt_count(x, normalized) for x in vals)
            else:
                fn = {"mean": np.mean, "median": np.median, "sum": np.sum, "each": np.mean}[p["junction_agg"]]
                v = float(fn(vals))
                text = fmt_count(v, normalized)
            hl = highlighted(k)             # highlighted introns bypass the display filters
            if v < p["min_junction"] and not hl:
                continue
            annotated = known_introns is not None and k in known_introns
            if p["annotated_only"] and known_introns is not None and not annotated and not hl:
                continue
            junctions[k] = dict(value=v, text=text, annotated=annotated, highlight=hl)
        tracks.append(dict(name=g, label=label, color=tcolor, curves=curves, disp=disp,
                           junctions=junctions, members=members))
    return tracks


# --------------------------------------------------------------------------
# Drawing helpers
# --------------------------------------------------------------------------

def _curve_xy(y, start, cmap):
    """Step-curve points of a per-base signal, run-length or max-binned."""
    n = len(y)
    edges = cmap.fwd(np.arange(start, start + n + 1) - 0.5)
    change = np.flatnonzero(np.diff(y)) + 1
    if len(change) < MAX_POINTS:
        idx = np.concatenate([[0], change])
        xs, ys = edges[idx], y[idx]
    else:
        centers = (edges[:-1] + edges[1:]) / 2
        bins = np.linspace(edges[0], edges[-1], MAX_POINTS + 1)[:-1]
        idx = np.unique(np.clip(np.searchsorted(centers, bins), 0, n - 1))
        xs, ys = edges[idx], np.maximum.reduceat(y, idx)
    xs = np.append(xs, edges[-1])
    ys = np.append(ys, ys[-1])
    real = np.round(cmap.inv(xs + 0.5)).astype(int)
    return xs, ys, real


def _bezier(x0, y0, x1, y1, c, n=32):
    t = np.linspace(0, 1, n)
    u = 1 - t
    xs = u**3 * x0 + 3 * u**2 * t * x0 + 3 * u * t**2 * x1 + t**3 * x1
    ys = u**3 * y0 + 3 * u**2 * t * c + 3 * u * t**2 * c + t**3 * y1
    return xs, ys


def _arcs(track, region, cmap, p, vmax):
    """Geometry of every junction arc of a track."""
    start = region.start
    disp = track["disp"]
    ymax = float(disp.max()) if len(disp) else 0.0
    ymax = ymax or 1.0
    items = sorted(track["junctions"].items())
    sides = []
    for i, _ in enumerate(items):
        if p["arc_layout"] == "above":
            sides.append(1)
        elif p["arc_layout"] == "below":
            sides.append(-1)
        else:
            sides.append(1 if i % 2 == 0 else -1)
    spans = {s: sorted(b - a for (a, b), _ in [it for it, sd in zip(items, sides) if sd == s]) for s in (1, -1)}

    arcs = []
    for ((a, b), j), side in zip(items, sides):
        x0, x1 = float(cmap.fwd(a - 0.5)), float(cmap.fwd(b + 0.5))
        sp = spans[side]
        frac = sp.index(b - a) / (len(sp) - 1) if len(sp) > 1 else 0.5
        if side > 0:
            y0 = float(disp[a - 1 - start]) if a - 1 >= start else 0.0
            y1 = float(disp[b + 1 - start]) if b + 1 - start < len(disp) else 0.0
            peak = max(y0, y1) + ymax * (0.12 + 0.28 * frac)
        else:
            y0 = y1 = 0.0
            peak = -ymax * (0.15 + 0.25 * frac)
        c = (peak - (y0 + y1) / 8) / 0.75
        xs, ys = _bezier(x0, y0, x1, y1, c)
        w = 0.6 + (p["arc_width"] - 0.6) * math.sqrt(j["value"] / vmax) if vmax > 0 else 1
        arcs.append(dict(key=(a, b), xs=xs, ys=ys, peak=peak, xmid=(x0 + x1) / 2,
                         width=round(w * 2) / 2, j=j))
    return arcs


def _annotation_rows(transcripts, p):
    if p["ann_mode"] == "genes":
        rows = collapse_genes(transcripts)
    else:
        rows = transcripts
        if p["transcripts"]:
            chosen = set(p["transcripts"])
            rows = [t for t in rows if t["id"] in chosen]
        else:
            rows = rows[: int(p["max_transcripts"])]
    return rows


def _draw_annotation(fig, row, rows, region, cmap, p):
    start, end = region.start, region.end
    lines_x, lines_y, arr = [], [], {"+": ([], []), "-": ([], [])}
    view = cmap.fwd(end + 0.5) - cmap.fwd(start - 0.5)
    spacing = view / 60
    for i, t in enumerate(rows):
        y = -i
        s, e = max(t["start"], start), min(t["end"], end)
        if s > e:
            continue
        xs, xe = float(cmap.fwd(s - 0.5)), float(cmap.fwd(e + 0.5))
        lines_x += [xs, xe, None]
        lines_y += [y, y, None]
        if t["strand"] in arr:
            pts = np.arange(xs + spacing / 2, xe, spacing)
            arr[t["strand"]][0].extend(pts)
            arr[t["strand"]][1].extend([y] * len(pts))

    fig.add_trace(go.Scatter(x=lines_x, y=lines_y, mode="lines", line=dict(color=ANN_COLOR, width=1),
                             hoverinfo="skip", showlegend=False), row=row, col=1)
    for strand, sym in (("+", "triangle-right"), ("-", "triangle-left")):
        if arr[strand][0]:
            fig.add_trace(go.Scatter(x=arr[strand][0], y=arr[strand][1], mode="markers",
                                     marker=dict(symbol=sym, size=5, color=ANN_COLOR),
                                     hoverinfo="skip", showlegend=False), row=row, col=1)

    for i, t in enumerate(rows):
        y = -i
        px, py = [], []
        thin = 0.12 if t["cds"] else 0.22
        for parts, h in ((t["exons"], thin), (t["cds"], 0.22)):
            for a, b in parts:
                a, b = max(a, start), min(b, end)
                if a > b:
                    continue
                x0, x1 = float(cmap.fwd(a - 0.5)), float(cmap.fwd(b + 0.5))
                px += [x0, x1, x1, x0, x0, None]
                py += [y - h, y - h, y + h, y + h, y - h, None]
        if not px:
            continue
        info = f"{t['name']}" + (f" | {t['gene']}" if t["gene"] != t["name"] else "") \
            + (f" | {t['id']}" if t["id"] not in (t["name"], t["gene"]) else "") \
            + (f" | {t['biotype']}" if t["biotype"] else "") + f" ({t['strand']})"
        fig.add_trace(go.Scatter(x=px, y=py, mode="lines", fill="toself", fillcolor=ANN_COLOR,
                                 line=dict(color=ANN_COLOR, width=0.5), hoveron="fills",
                                 hoverinfo="text", text=info, name=t["name"], showlegend=False),
                      row=row, col=1)
        lx = float(cmap.fwd(max(t["start"], start) - 0.5))
        fig.add_annotation(x=lx, y=y + 0.24, text=t["name"], showarrow=False, xanchor="left", yanchor="bottom",
                           font=dict(size=p["font_size"] * 0.8, color=ANN_COLOR), row=row, col=1)
    fig.update_yaxes(range=[-len(rows) + 0.6, 0.85], visible=False, row=row, col=1)


# --------------------------------------------------------------------------
# Main entry
# --------------------------------------------------------------------------

def build_figure(region, params=None, transcripts=None, interactive=True, return_cmap=False):
    """region: data.RegionData; transcripts: Annotation.region(...) output or None.

    With return_cmap, also returns the CoordMap (display <-> real coordinates)."""
    p = merged_params(params)
    known = annotated_introns(transcripts) if transcripts else None
    tracks = build_tracks(region, p, known)
    if not tracks:
        fig = go.Figure(layout=dict(template="simple_white", height=250, xaxis_visible=False, yaxis_visible=False))
        fig.add_annotation(text="No samples selected", showarrow=False, x=0.5, y=0.5, xref="paper", yref="paper",
                           font=dict(size=16, color="#889"))
        return (fig, CoordMap(region.start, region.end)) if return_cmap else fig
    normalized = p["normalize"] == "cpm"
    start, end = region.start, region.end

    intervals = ()
    if p["shrink"]:
        intervals = intersect_introns({k for t in tracks for k in t["junctions"]})
    cmap = CoordMap(start, end, intervals, float(p["shrink_power"]))

    ann_rows = _annotation_rows(transcripts, p) if transcripts and p["ann_mode"] != "none" else []
    n = len(tracks)
    ann_h = (len(ann_rows) * (p["font_size"] * 2 + 10) + 16) if ann_rows else 0
    heights = [p["track_height"]] * n + ([ann_h] if ann_rows else [])
    gap = 6
    plot_h = sum(heights) + gap * (len(heights) - 1)
    margin = dict(l=70, r=20, t=50 if p["show_title"] else 15, b=45)
    fig = make_subplots(rows=len(heights), cols=1, shared_xaxes=True, row_heights=heights,
                        vertical_spacing=gap / plot_h if len(heights) > 1 else 0)

    vmax = max((j["value"] for t in tracks for j in t["junctions"].values()), default=1.0)
    any_hl = any(j["highlight"] for t in tracks for j in t["junctions"].values())
    hcolor = p["highlight_color"] or DEFAULTS["highlight_color"]
    dim = 0.3 if any_hl and p["dim_others"] else 1.0
    halo = "1px 1px 1px #fff, -1px -1px 1px #fff, 1px -1px 1px #fff, -1px 1px 1px #fff"
    font = dict(size=p["font_size"])
    track_ranges = []
    for r, t in enumerate(tracks, start=1):
        overlaid = len(t["curves"]) > 1
        for name, y, color in t["curves"]:
            xs, ys, real = _curve_xy(y, start, cmap)
            fig.add_trace(go.Scatter(
                x=xs, y=ys, mode="lines", line=dict(width=0, shape="hv", color=color),
                fill="tozeroy", fillcolor=color, opacity=p["alpha"], name=name, customdata=real,
                hovertemplate=f"<b>{name}</b><br>{region.chrom}:%{{customdata:,}}<br>%{{y:.4~g}}<extra></extra>",
                showlegend=False), row=r, col=1)
            if overlaid:                    # opaque outline so each overlaid sample stays visible
                fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(width=1, shape="hv", color=color),
                                         hoverinfo="skip", showlegend=False), row=r, col=1)

        arcs = _arcs(t, region, cmap, p, vmax)
        buckets = {}
        for a in arcs:
            dash = "dot" if p["novel_dashed"] and known is not None and not a["j"]["annotated"] else "solid"
            buckets.setdefault((a["j"]["highlight"], a["width"], dash), []).append(a)
        for (hl, w, dash), items in sorted(buckets.items()):    # highlighted arcs drawn last (on top)
            xs, ys, cd = [], [], []
            for a in items:
                (ja, jb), j = a["key"], a["j"]
                tag = "" if known is None else ("annotated" if j["annotated"] else "<b>novel</b>")
                info = (f"<b>{t['label']}</b><br>intron {region.chrom}:{ja:,}-{jb:,} ({jb - ja + 1:,} bp)"
                        f"<br>reads: {j['text']}" + (f"<br>{tag}" if tag else "")
                        + ("<br>★ highlighted" if hl else "") + "<br><i>click to toggle highlight</i>")
                xs += list(a["xs"]) + [None]
                ys += list(a["ys"]) + [None]
                cd += [[info, f"{ja}-{jb}"]] * len(a["xs"]) + [None]
            fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", customdata=cd,
                                     line=dict(color=hcolor if hl else t["color"], width=w + 2 if hl else w, dash=dash),
                                     opacity=1.0 if hl else dim, hovertemplate="%{customdata[0]}<extra></extra>",
                                     showlegend=False), row=r, col=1)
        if p["show_labels"] and arcs:
            for hl in (False, True):
                sel = [a for a in arcs if a["j"]["highlight"] == hl]
                if not sel:
                    continue
                fig.add_trace(go.Scatter(
                    x=[a["xmid"] for a in sel], y=[a["peak"] for a in sel], mode="text",
                    text=[f"<b>{a['j']['text']}</b>" if hl else a["j"]["text"] for a in sel],
                    textfont=dict(size=p["font_size"] * (0.9 if hl else 0.75), color=hcolor if hl else "#222",
                                  shadow=halo),
                    opacity=1.0 if hl else dim, hoverinfo="skip", showlegend=False), row=r, col=1)

        ymax = float(t["disp"].max()) if len(t["disp"]) else 0.0
        top = max([ymax] + [a["peak"] for a in arcs]) or 1.0
        bottom = min([0.0] + [a["peak"] for a in arcs])
        track_ranges.append((ymax, top, bottom))
        fig.add_annotation(text=f"<b>{t['label']}</b>", x=0.005, y=0.98, xref="x domain", yref="y domain",
                           xanchor="left", yanchor="top", showarrow=False,
                           font=dict(size=p["font_size"], color=t["color"] if p["color_by"] else "#222"),
                           row=r, col=1)

    for r, (ymax, top, bottom) in enumerate(track_ranges, start=1):
        if p["fix_y"]:
            ymax = max(x[0] for x in track_ranges)
            top = max(x[1] for x in track_ranges)
            bottom = min(x[2] for x in track_ranges)
        pad = (top - bottom) * 0.08
        ticks = nice_ticks(ymax)
        fig.update_yaxes(range=[bottom - (pad if bottom < 0 else 0), top + pad], tickvals=ticks,
                         ticktext=[fmt_count(v, normalized) for v in ticks],
                         zeroline=True, zerolinecolor="#999", zerolinewidth=1, showgrid=False,
                         ticks="outside", tickfont=font, fixedrange=False, row=r, col=1)

    if ann_rows:
        _draw_annotation(fig, n + 1, ann_rows, region, cmap, p)

    if p["highlight_band"]:
        for a, b in parse_introns(p["highlight"]):
            a, b = max(a, start), min(b, end)
            if a <= b:
                fig.add_vrect(x0=float(cmap.fwd(a - 0.5)), x1=float(cmap.fwd(b + 0.5)), row="all", col=1,
                              fillcolor=hcolor, opacity=0.08, layer="below", line_width=0)

    x0, x1 = float(cmap.fwd(start - 0.5)), float(cmap.fwd(end + 0.5))
    xaxis = dict(range=[x0, x1], showgrid=False, ticks="outside", tickfont=font, tickformat=",d")
    if not cmap.identity:
        k = 7
        pos = np.linspace(x0, x1, k + 2)[1:-1]
        real = cmap.inv(pos)
        step = nice_step(end - start, 40)
        real = np.unique(np.round(real / step) * step)
        xaxis.update(tickvals=list(cmap.fwd(real)), ticktext=[f"{int(v):,}" for v in real])
    fig.update_xaxes(**xaxis)
    fig.update_xaxes(showticklabels=True, row=len(heights), col=1)

    strand_txt = {"+": "plus strand", "-": "minus strand", "combined": ""}.get(p["strand"], "")
    if region.libtype == "NONE":
        strand_txt = ""
    title = f"{region.chrom}:{start:,}-{end:,}" + (f" · {strand_txt}" if strand_txt else "") \
        + (" · CPM" if normalized else "") + (" · introns shrunk" if not cmap.identity else "")
    fig.update_layout(
        height=plot_h + margin["t"] + margin["b"],
        width=None if interactive and not p["width"] else p["width"],
        margin=margin, template="simple_white", font=dict(size=p["font_size"], family="Arial, Helvetica, sans-serif"),
        title=dict(text=title, x=0.01, font=dict(size=p["font_size"] * 1.15)) if p["show_title"] else None,
        hovermode="closest", dragmode="zoom", showlegend=False,
        plot_bgcolor="white", paper_bgcolor="white",
    )
    return (fig, cmap) if return_cmap else fig
