"""Interactive sashimi plot app (Dash).

BAMs are read only when the region / library settings change ("Load");
every other control re-renders the figure from the in-memory data.
"""

import base64
import json
import os
import re
import traceback

from dash import Dash, Input, Output, State, ctx, dcc, html, no_update

from . import __version__
from .annotation import Annotation, annotated_introns
from .data import LIBTYPES, ReadFilter, format_region, junction_table, load_region, parse_region, read_samples
from .plot import DEFAULTS, PALETTES, build_figure, format_introns, merged_params, parse_introns

# Server-side state (the app is meant to be run locally by one user).
STATE = {"region": None, "transcripts": None, "annotation": None, "ann_path": None, "cmap": None}


# --------------------------------------------------------------------------
# Small widget helpers
# --------------------------------------------------------------------------

def row(label, *children, hint=None):
    kids = [html.Label(label)] if label else []
    kids += list(children)
    if hint:
        kids.append(html.Div(hint, className="hint"))
    return html.Div(kids, className="row")


def switch(id_, label, on):
    return dcc.Checklist(id=id_, options=[{"label": label, "value": "on"}], value=["on"] if on else [],
                         className="check", inputStyle={"marginRight": "4px"})


def radio(id_, options, value):
    opts = [{"label": l, "value": v} for v, l in options]
    return dcc.RadioItems(id=id_, options=opts, value=value, inline=True, className="radio")


def slider(id_, lo, hi, step, value, marks=None):
    return dcc.Slider(id=id_, min=lo, max=hi, step=step, value=value, marks=marks,
                      tooltip={"placement": "bottom", "always_visible": False}, updatemode="mouseup")


# Plot parameters <-> controls.  kind: 'value' (as is), 'bool' (Checklist), 'num'
PARAM_CONTROLS = [
    ("samples", "p-samples", "value"),
    ("group_by", "p-group_by", "value"),
    ("group_mode", "p-group_mode", "value"),
    ("color_by", "p-color_by", "value"),
    ("label_by", "p-label_by", "value"),
    ("palette", "p-palette", "value"),
    ("custom_colors", "p-custom_colors", "value"),
    ("strand", "p-strand", "value"),
    ("normalize", "p-normalize", "value"),
    ("min_junction", "p-min_junction", "num"),
    ("junction_agg", "p-junction_agg", "value"),
    ("show_labels", "p-show_labels", "bool"),
    ("annotated_only", "p-annotated_only", "bool"),
    ("novel_dashed", "p-novel_dashed", "bool"),
    ("arc_layout", "p-arc_layout", "value"),
    ("arc_width", "p-arc_width", "num"),
    ("highlight", "p-highlight", "value"),
    ("highlight_color", "p-highlight_color", "value"),
    ("highlight_tol", "p-highlight_tol", "num"),
    ("dim_others", "p-dim_others", "bool"),
    ("highlight_band", "p-highlight_band", "bool"),
    ("fix_y", "p-fix_y", "bool"),
    ("shrink", "p-shrink", "bool"),
    ("shrink_power", "p-shrink_power", "num"),
    ("alpha", "p-alpha", "num"),
    ("track_height", "p-track_height", "num"),
    ("width", "p-width", "num"),
    ("font_size", "p-font_size", "num"),
    ("show_title", "p-show_title", "bool"),
    ("ann_mode", "p-ann_mode", "value"),
    ("max_transcripts", "p-max_transcripts", "num"),
    ("transcripts", "p-transcripts", "value"),
]


def collect_params(values):
    p = {}
    for (name, _, kind), v in zip(PARAM_CONTROLS, values):
        if kind == "bool":
            v = bool(v)
        elif kind == "num":
            v = DEFAULTS[name] if v in (None, "") else float(v)
        p[name] = v
    p["highlight"] = p["highlight"] or ""
    for k in ("track_height", "width", "max_transcripts", "highlight_tol"):
        p[k] = int(p[k] or 0)
    if not p["samples"]:
        p["samples"] = []
    if not p["transcripts"]:
        p["transcripts"] = None
    return p


def params_to_controls(p):
    out = []
    for name, _, kind in PARAM_CONTROLS:
        v = p.get(name, DEFAULTS[name])
        if name == "samples" and not v:
            out.append(no_update)       # settings without a sample list keep the current selection
            continue
        out.append((["on"] if v else []) if kind == "bool" else v)
    return out


def placeholder_figure(text):
    return {"layout": {"xaxis": {"visible": False}, "yaxis": {"visible": False}, "height": 300,
                       "plot_bgcolor": "white", "paper_bgcolor": "white",
                       "annotations": [{"text": text, "showarrow": False, "xref": "paper", "yref": "paper",
                                        "x": 0.5, "y": 0.5, "font": {"size": 16, "color": "#889"}}]}}


# --------------------------------------------------------------------------
# Layout
# --------------------------------------------------------------------------

def _sheet_columns(bam):
    try:
        cols = [c for c in read_samples([os.path.expanduser(x) for x in bam.split()]).columns if c != "bam"]
    except Exception:
        cols = ["id"]
    return cols, [c for c in cols if c != "id"]


def make_layout(init):
    d = merged_params(init.get("params"))
    cols, meta = _sheet_columns(init.get("bam") or "")
    sidebar = html.Div(className="sidebar", children=[
        html.Div(["SashimiPro", html.Span(f" v{__version__}")], className="brand"),

        html.Details(open=True, children=[
            html.Summary("Data"),
            row("Samples (TSV sheet or BAM paths, space separated)",
                dcc.Input(id="in-bam", type="text", value=init.get("bam", ""), debounce=True)),
            row("Annotation GTF (optional)", dcc.Input(id="in-gtf", type="text", value=init.get("gtf", ""), debounce=True)),
            row("Region or gene name", dcc.Input(id="in-region", type="text", value=init.get("region", ""),
                                                 debounce=True, placeholder="chr1:1,000-9,000 or GAPDH"),
                hint="Press Enter or Load. Genes get ±10% flanks."),
            row("Library strandedness", dcc.Dropdown(id="in-libtype", options=LIBTYPES, value=init.get("libtype", "NONE"),
                                                     clearable=False)),
            row("Read filters", html.Div(className="inline", children=[
                html.Div([html.Span("min MAPQ ", className="lbl"),
                          dcc.Input(id="in-mapq", type="number", value=init.get("mapq", 0), min=0)]),
                html.Div([switch("in-dup", "skip duplicates", False), switch("in-sec", "skip secondary", True)]),
            ])),
            html.Div(className="inline", children=[
                html.Button("Load", id="btn-load", className="primary"),
                html.Button("◀", id="btn-left", title="Move left by half a window", className="nav"),
                html.Button("▶", id="btn-right", title="Move right by half a window", className="nav"),
                html.Button("zoom out", id="btn-out", title="Zoom out 3×"),
                html.Button("zoom in", id="btn-view", title="Reload just the range you zoomed into on the plot"),
            ]),
            html.Div(id="status", className="status"),
        ]),

        html.Details(open=True, children=[
            html.Summary("Samples & colors"),
            row("Samples shown", dcc.Dropdown(id="p-samples", multi=True, options=[], value=[])),
            row("Group (overlay) by column", dcc.Dropdown(id="p-group_by", options=meta, value=d["group_by"])),
            row("Within a group", radio("p-group_mode", [("overlay", "overlay"), ("mean", "mean"), ("median", "median"),
                                                         ("sum", "sum")], d["group_mode"])),
            row("Color by column", dcc.Dropdown(id="p-color_by", options=meta, value=d["color_by"])),
            row("Label column (ungrouped)", dcc.Dropdown(id="p-label_by", options=cols, value=d["label_by"], clearable=False)),
            row("Palette", dcc.Dropdown(id="p-palette", options=list(PALETTES), value=d["palette"], clearable=False)),
            row("Custom colors", dcc.Input(id="p-custom_colors", type="text", value=d["custom_colors"], debounce=True,
                                           placeholder="#d62728, steelblue, ..."), hint="Overrides the palette, in level order."),
        ]),

        html.Details(open=True, children=[
            html.Summary("Junctions"),
            row("Minimum reads to draw (after aggregation / CPM)",
                dcc.Input(id="p-min_junction", type="number", value=d["min_junction"], min=0, step=1, debounce=True)),
            row("Aggregate counts in a group", radio("p-junction_agg", [("mean", "mean"), ("median", "median"),
                                                                       ("sum", "sum"), ("each", "list each")],
                                                     d["junction_agg"])),
            row("Arc placement", radio("p-arc_layout", [("alternate", "alternate"), ("above", "above"), ("below", "below")],
                                       d["arc_layout"])),
            row("Max arc width", slider("p-arc_width", 1, 12, 0.5, d["arc_width"], {1: "1", 6: "6", 12: "12"})),
            switch("p-show_labels", "show counts", d["show_labels"]),
            switch("p-annotated_only", "only annotated junctions (needs GTF)", d["annotated_only"]),
            switch("p-novel_dashed", "dotted line for novel junctions", d["novel_dashed"]),
        ]),

        html.Details(open=True, children=[
            html.Summary("Highlight introns"),
            row("Introns to highlight", dcc.Input(id="p-highlight", type="text", value=d["highlight"], debounce=True,
                                                  placeholder="chr7:128,940,001-128,940,500; ..."),
                hint="Separate several with ';'. Or click an arc on the plot to add / remove it."),
            html.Div(className="inline", children=[
                html.Div([html.Span("color", className="lbl"),
                          dcc.Input(id="p-highlight_color", type="text", value=d["highlight_color"], debounce=True)]),
                html.Div([html.Span("match ± bp", className="lbl"),
                          dcc.Input(id="p-highlight_tol", type="number", value=d["highlight_tol"], min=0, step=1,
                                    debounce=True)]),
            ]),
            switch("p-dim_others", "fade the other arcs", d["dim_others"]),
            switch("p-highlight_band", "shade the intron across tracks", d["highlight_band"]),
            html.Button("Clear", id="btn-hl-clear", style={"marginTop": "4px"}),
        ]),

        html.Details(open=True, children=[
            html.Summary("Coverage & axes"),
            row("Strand", radio("p-strand", [("combined", "both"), ("+", "plus"), ("-", "minus")], d["strand"])),
            row("Scale", radio("p-normalize", [("raw", "raw reads"), ("cpm", "per million mapped")], d["normalize"])),
            switch("p-fix_y", "same y scale for all tracks", d["fix_y"]),
            switch("p-shrink", "shrink introns", d["shrink"]),
            row("Shrink strength (intron length ^ x)", slider("p-shrink_power", 0.3, 1.0, 0.05, d["shrink_power"],
                                                                {0.3: "0.3", 0.7: "0.7", 1: "1"})),
            row("Fill opacity", slider("p-alpha", 0.1, 1.0, 0.05, d["alpha"], {0.1: "0.1", 0.5: "0.5", 1: "1"})),
        ]),

        html.Details(open=False, children=[
            html.Summary("Size & text"),
            row("Track height (px)", slider("p-track_height", 50, 300, 10, d["track_height"], {50: "50", 150: "150", 300: "300"})),
            row("Width (px, 0 = fit window)", dcc.Input(id="p-width", type="number", value=init.get("width", 0), min=0,
                                                        step=50, debounce=True)),
            row("Font size", slider("p-font_size", 8, 22, 1, d["font_size"], {8: "8", 14: "14", 22: "22"})),
            switch("p-show_title", "show title", d["show_title"]),
        ]),

        html.Details(open=False, children=[
            html.Summary("Annotation"),
            row("Show", radio("p-ann_mode", [("transcripts", "transcripts"), ("genes", "collapsed genes"),
                                             ("none", "none")], d["ann_mode"])),
            row("Max transcripts", dcc.Input(id="p-max_transcripts", type="number", value=d["max_transcripts"], min=1,
                                             debounce=True)),
            row("Pick transcripts", dcc.Dropdown(id="p-transcripts", multi=True, options=[], value=d["transcripts"],
                                                 placeholder="(first N)")),
        ]),

        html.Details(open=True, children=[
            html.Summary("Export"),
            row("File name", dcc.Input(id="ex-name", type="text", value="sashimi", debounce=True)),
            html.Div(className="inline", children=[
                html.Div([html.Span("Format", className="lbl"),
                          dcc.Dropdown(id="ex-format", options=["pdf", "svg", "png"], value="pdf", clearable=False)]),
                html.Div([html.Span("PNG resolution ×", className="lbl"),
                          dcc.Input(id="ex-scale", type="number", value=2, min=1, max=6, step=1)]),
            ]),
            html.Div("× only affects PNG: 2 = twice the pixels in each direction (sharper).", className="hint"),
            html.Div(className="inline", style={"marginTop": "6px"}, children=[
                html.Button("Figure", id="btn-export"),
                html.Button("Junctions", id="btn-junc", title="All junction counts per sample (TSV)"),
                html.Button("Settings", id="btn-settings", title="Save all plot parameters as JSON"),
            ]),
            dcc.Upload(id="up-settings", children=html.Div("Drop / click to load a settings JSON", className="hint"),
                       style={"border": "1px dashed #c8cdd3", "borderRadius": "4px", "padding": "6px",
                              "textAlign": "center", "marginTop": "8px", "cursor": "pointer"}),
            html.Div(id="export-status", className="status"),
            dcc.Download(id="dl"),
        ]),
    ])
    main = html.Div(className="main", children=[
        html.Div(id="load-banner", className="load-banner", style={"display": "none"}, children=[
            html.Div(className="load-bar"),
            html.Div("Loading BAMs… (the first load of a large GTF can take a minute)", className="load-text"),
        ]),
        dcc.Loading(dcc.Graph(id="graph", figure=placeholder_figure("Loading…" if init.get("region") else
                                                                    "Enter samples and a region, then press Load"),
                              config={
            "displaylogo": False, "scrollZoom": False,
            "toImageButtonOptions": {"format": "svg", "filename": "sashimi"},
            "modeBarButtonsToRemove": ["lasso2d", "select2d"],
        }), type="dot", delay_show=400),
    ])
    return html.Div([dcc.Store(id="loaded"), html.Div([sidebar, main], className="layout")])


# --------------------------------------------------------------------------
# App
# --------------------------------------------------------------------------

def _resolve_region(text):
    try:
        return parse_region(text)
    except ValueError:
        ann = STATE["annotation"]
        hit = ann.find_gene(text) if ann and text and text.strip() else None
        if not hit:
            raise ValueError(f"'{text}' is neither a region (chr:start-end) nor a gene in the GTF")
        chrom, s, e = hit
        flank = max(50, int((e - s) * 0.1))
        return chrom, max(1, s - flank), e + flank


def _view_range(relayout):
    """Real coordinates of the current zoom (from Plotly relayoutData)."""
    if not relayout or STATE["cmap"] is None:
        return None
    lo = hi = None
    for k, v in relayout.items():
        if re.fullmatch(r"xaxis\d*\.range\[0\]", k):
            lo = v
        elif re.fullmatch(r"xaxis\d*\.range\[1\]", k):
            hi = v
        elif re.fullmatch(r"xaxis\d*\.range", k):
            lo, hi = v
    if lo is None or hi is None:
        return None
    cm = STATE["cmap"]
    return int(round(float(cm.inv(lo)) + 0.5)), int(round(float(cm.inv(hi)) - 0.5))


def create_app(init=None):
    init = init or {}
    app = Dash(__name__, title="SashimiPro", suppress_callback_exceptions=True)
    app.layout = make_layout(init)

    @app.callback(
        Output("loaded", "data"), Output("status", "children"), Output("status", "className"),
        Output("in-region", "value"),
        Output("p-samples", "options"), Output("p-samples", "value"),
        Output("p-group_by", "options"), Output("p-color_by", "options"), Output("p-label_by", "options"),
        Output("p-label_by", "value"), Output("p-transcripts", "options"),
        Input("btn-load", "n_clicks"), Input("in-region", "value"), Input("btn-left", "n_clicks"),
        Input("btn-right", "n_clicks"), Input("btn-out", "n_clicks"), Input("btn-view", "n_clicks"),
        State("in-bam", "value"), State("in-gtf", "value"), State("in-libtype", "value"), State("in-mapq", "value"),
        State("in-dup", "value"), State("in-sec", "value"), State("graph", "relayoutData"),
        State("p-samples", "value"), State("p-label_by", "value"), State("loaded", "data"),
        running=[
            (Output("load-banner", "style"), {"display": "block"}, {"display": "none"}),
            *[(Output(b, "disabled"), True, False) for b in ("btn-load", "btn-left", "btn-right", "btn-out", "btn-view")],
        ],
    )
    def load(_n, region_text, _l, _r, _o, _v, bam, gtf, libtype, mapq, dup, sec, relayout, chosen, label_by, loaded):
        trig = ctx.triggered_id
        try:
            if not bam or not bam.strip():
                return no_update, "Enter a sample sheet or BAM file(s), a region, then Load.", "status", *[no_update] * 8
            gtf = (gtf or "").strip()
            if gtf != STATE["ann_path"]:
                STATE["annotation"] = Annotation(os.path.expanduser(gtf)) if gtf else None
                STATE["ann_path"] = gtf
            if not region_text:
                return no_update, "Enter a region or gene name.", "status", *[no_update] * 8
            chrom, start, end = _resolve_region(region_text)
            span = end - start + 1
            if trig == "btn-left":
                start, end = max(1, start - span // 2), max(1, start - span // 2) + span - 1
            elif trig == "btn-right":
                start, end = start + span // 2, end + span // 2
            elif trig == "btn-out":
                start, end = max(1, start - span), end + span
            elif trig == "btn-view":
                vr = _view_range(relayout)
                if vr:
                    start, end = max(start, vr[0]), min(end, vr[1])
            if end - start > 2_000_000:
                raise ValueError(f"Region is {end - start + 1:,} bp; please stay below 2 Mb")

            samples = read_samples([os.path.expanduser(x) for x in bam.split()])
            rf = ReadFilter(min_mapq=int(mapq or 0), skip_duplicates=bool(dup), skip_secondary=bool(sec))
            region = load_region(samples, chrom, start, end, libtype, rf)
            ann = STATE["annotation"]
            transcripts = ann.region(chrom, start, end) if ann else None
            STATE["region"], STATE["transcripts"] = region, transcripts

            ids = list(region.signals)
            # keep the user's sample selection when only the region changed
            same = loaded and loaded.get("ids") == ids
            value = [s for s in (chosen or []) if s in ids] if same else ids
            value = value or ids
            cols = [c for c in region.samples.columns if c not in ("bam",)]
            meta = [c for c in cols if c != "id"]
            tx_opts = [{"label": f"{t['name']} ({t['gene']})", "value": t["id"]} for t in (transcripts or [])]
            msg = (f"Loaded {len(ids)} sample(s), {format_region(chrom, start, end)} "
                   f"({end - start + 1:,} bp)" + (f", {len(transcripts)} transcripts" if transcripts is not None else ""))
            if region.warnings:
                msg += "\n⚠ " + "\n⚠ ".join(region.warnings)
            data = {"ids": ids, "region": format_region(chrom, start, end), "libtype": libtype,
                    "token": (loaded or {}).get("token", 0) + 1}
            label_by = label_by if label_by in cols else "id"
            return (data, msg, "status", format_region(chrom, start, end), ids, value, meta, meta, cols, label_by,
                    tx_opts)
        except Exception as e:  # show any loading error in the sidebar
            traceback.print_exc()
            return no_update, f"Error: {e}", "status err", *[no_update] * 8

    @app.callback(
        Output("graph", "figure"),
        Input("loaded", "data"), *[Input(cid, "value") for _, cid, _ in PARAM_CONTROLS],
    )
    def render(loaded, *values):
        region = STATE["region"]
        if not loaded or region is None:
            return placeholder_figure("Enter samples and a region, then press Load")
        p = collect_params(values)
        fig, STATE["cmap"] = build_figure(region, p, STATE["transcripts"], return_cmap=True)
        # keep the zoom while tweaking parameters; reset when region or x scaling changes
        fig.update_layout(uirevision=f"{loaded['region']}|{p['shrink']}|{p['shrink_power']}")
        return fig

    @app.callback(
        Output("dl", "data"), Output("export-status", "children"),
        Input("btn-export", "n_clicks"), Input("btn-junc", "n_clicks"), Input("btn-settings", "n_clicks"),
        State("ex-name", "value"), State("ex-format", "value"), State("ex-scale", "value"),
        *[State(cid, "value") for _, cid, _ in PARAM_CONTROLS],
        prevent_initial_call=True,
    )
    def export(_f, _j, _s, name, fmt, scale, *values):
        name = (name or "sashimi").strip()
        p = collect_params(values)
        region = STATE["region"]
        try:
            if ctx.triggered_id == "btn-settings":
                region_txt = format_region(region.chrom, region.start, region.end) if region else None
                payload = {"region": region_txt, "libtype": region.libtype if region else None, "params": p}
                return dcc.send_string(json.dumps(payload, indent=2), f"{name}.settings.json"), "Settings saved."
            if region is None:
                return no_update, "Nothing loaded yet."
            if ctx.triggered_id == "btn-junc":
                df = junction_table(region)
                if STATE["transcripts"] is not None:
                    known = annotated_introns(STATE["transcripts"])
                    df["annotated"] = [(a, b) in known for a, b in zip(df["intron_start"], df["intron_end"])]
                return dcc.send_data_frame(df.to_csv, f"{name}.junctions.tsv", sep="\t", index=False), \
                    f"{len(df)} junction rows exported."
            if not p["width"]:
                p["width"] = 1000
            fig = build_figure(region, p, STATE["transcripts"], interactive=False)
            img = fig.to_image(format=fmt, scale=float(scale or 1) if fmt == "png" else 1)
            return dcc.send_bytes(img, f"{name}.{fmt}"), f"Saved {name}.{fmt}"
        except Exception as e:
            traceback.print_exc()
            return no_update, f"Export failed: {e}"

    @app.callback(
        Output("p-highlight", "value", allow_duplicate=True),
        Input("graph", "clickData"), Input("btn-hl-clear", "n_clicks"),
        State("p-highlight", "value"),
        prevent_initial_call=True,
    )
    def toggle_highlight(click, _clear, current):
        if ctx.triggered_id == "btn-hl-clear":
            return ""
        try:
            key = click["points"][0]["customdata"][1]       # junction arcs carry [hover text, "start-end"]
        except (TypeError, KeyError, IndexError):
            return no_update                                # clicked coverage / annotation: ignore
        a, b = (int(x) for x in key.split("-"))
        introns = parse_introns(current)
        introns = [iv for iv in introns if iv != (a, b)] if (a, b) in introns else introns + [(a, b)]
        region = STATE["region"]
        return format_introns(region.chrom if region else "", introns)

    @app.callback(
        *[Output(cid, "value", allow_duplicate=True) for _, cid, _ in PARAM_CONTROLS],
        Output("export-status", "children", allow_duplicate=True),
        Input("up-settings", "contents"),
        prevent_initial_call=True,
    )
    def restore(contents):
        try:
            payload = json.loads(base64.b64decode(contents.split(",", 1)[1]))
            p = merged_params(payload.get("params", payload))
            return *params_to_controls(p), "Settings applied."
        except Exception as e:
            return *[no_update] * len(PARAM_CONTROLS), f"Could not read settings: {e}"

    return app


def run(init, host="127.0.0.1", port=8050, debug=False):
    app = create_app(init)
    app.run(host=host, port=port, debug=debug)
