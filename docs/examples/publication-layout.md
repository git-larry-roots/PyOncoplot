# Compact publication layout (synthetic cohort)

Run from the repository root with the existing environment:

~~~bash
PYTHONPATH=src .venv/bin/python docs/examples/publication_layout.py
~~~

The [complete executable example](publication_layout.py) writes
.venv/issue-2-publication/publication-layout.{png,pdf,svg}. No real patient
information is used. All styling is passed to oncoplot(..., backend="matplotlib",
options={...}); no returned axes are edited.

The example combines aligned top event counts, right gene-frequency bars,
compact bottom categorical and numeric tracks, sample IDs under the last track,
complete categorical legends and a horizontal numeric colorbar. The matrix takes
most of the width. EMPTY remains visible; hatched unresolved/not-assessed cells
are not negative calls or mutation events.

## Public controls

| Option | Meaning |
| --- | --- |
| `compact_layout=True` | Measure legend/text extents and reserve a bottom band; pack full categorical legends and horizontal numeric colorbars together when they fit. Requires bottom legend positions (or hidden legends). |
| `sample_id_padding=6` | Tick-label separation in points. Existing show_sample_ids, sample_id_position="bottom", sample_id_angle, and font_size_samples apply; bottom IDs automatically move below bottom metadata, including in legacy layout. |
| `tmb_bar_color="#476B8A"` | Override marginal bar color, including stacked segments. Omit to retain mutation/custom tmb_palette colors. |
| `tmb_bar_width=.72` | Width in sample units; empty spacing is 1 - width. Samples and centers never move. |
| `tmb_tick_interval=1` | Positive integer major-tick interval in displayed units. Use log10_transform_tmb=False for raw integer event counts. |
| `tmb_y_label="Recorded\nevents"` | Custom multiline label; enable existing show_tmb_y_label. |
| `gene_bar_label_format="count_percent"` | count, percent (existing default), or count_percent, e.g. 12 (24%); enable show_gene_bar_labels. Existing gene_bar_label_padding is a fraction of the bar scale and gene_bar_label_nudge is in data units. Compact layout fits full text inside the marginal panel. |
| `metadata_track_heights={"arm": .24, "purity": .32}` | Individual positive relative track heights, keyed by original column name. Unspecified tracks retain height 1. Heights also determine metadata panel weight relative to the matrix. |
| `metadata_track_spacing=.08` | Gap between tracks, in the same relative units. Existing buffer_metadata controls the panel-to-matrix gap. |
| `metadata_numeric_scales={"purity": {"min": 0, "max": 1, "ticks": [0, .5, 1], "missing_color": "#CCCCCC"}}` | Explicit display scale, colorbar ticks and missing fill. Each key is optional. Omitted bounds use observed range (0–1 for all missing). Out-of-range data saturate at endpoints; original values are untouched. Numeric bar tracks also clamp display height to the range. Effective scales/ticks are validated. |
| `tile_edgecolor="white", tile_linewidth=.3` | Cell border color and thickness, independent from row separators. |
| `row_separator_color="#666666", row_separator_linewidth=0` | Row separator color/width; zero removes separators. Pathway outlines remain independently styled by existing pathway options. |
| `title_gap_points=24` | Gap between overall title lower edge and plotting area upper edge; title can be cropped during assembly. |
| `legend_label_overrides={"Missense_Mutation": "Missense"}` | Exact original legend value to display string, before prettification/truncation; source data/order are unchanged. State keys may use the complete label, e.g. Cell state: unresolved. |
| `legend_title_overrides={"cohort_arm": "Arm", "Mutation Type": ""}` | Exact original column/heading to display string. Also applies to metadata track labels. Empty strings suppress selected headings; existing show_legend_titles=False suppresses generic and metadata headings globally. |

Existing metadata_legend_position="bottom" and
metadata_legend_orientation_heatmap="horizontal" already select horizontal
colorbars; no duplicate orientation API was added. Existing
metadata_legend_nrow/metadata_legend_ncol, legend font sizes, key sizes and
manual legend_offsets remain available in legacy layout. Compact mode uses
measured automatic column packing rather than manual row/column counts, and
rejects manual offsets and manual metadata legend row/column counts. Increase width/height or reduce fonts if the requested
layout cannot fit. Compact mode does not yet support continuous main-grid
variant-value colorbars; that combination explicitly raises NotImplementedError,
while legacy layout continues to render those grids.

## Scientific semantics and backends

Issue #1 existing preserve_empty_genes, cell_states, cell_state_palette,
gene_denominators, and row_summaries are reused, not duplicated. Percentage
labels require explicit denominators when cell states are supplied. Count-only
labels do not invent denominators. Layout changes never alter sample ordering,
event membership, mutation counts, or declared percentage denominators. External
summary provenance remains available in prepared_data.row_summary_external even when display headings are customized.

Matplotlib supports all publication controls above. Plotly supports display
label/title overrides, as well as its existing legend orientation and empty/state
controls. New grid color controls are Matplotlib-only because Plotly heatmaps
and scatter tiles do not provide equivalent borders. Nondefault
values for the other newly added publication controls raise a descriptive
NotImplementedError on Plotly rather than silently doing nothing. Plotly static
exports continue to depend on optional Kaleido; PNG/PDF/SVG acceptance for
this example is verified with Matplotlib, not a claimed Plotly export guarantee.

## Publication refinements (issue #5)

Run the [complete external-summary example](publication_refinements.py):

~~~bash
PYTHONPATH=src .venv/bin/python docs/examples/publication_refinements.py
~~~

It writes `.venv/issue-5-publication/publication-refinements.{png,pdf,svg}`.
It demonstrates supplied aggregate counts/denominators, independently colored
summary bars, gray hatching, neutral missing fill without a marker, multiline
metadata labels, shorter panels, and a shared bottom legend/colorbar shelf.
All data are synthetic; it never modifies returned axes.

| Public control | Meaning |
| --- | --- |
| `gene_bar_color="#476B8A"` | One shared color for computed **and** supplied row-summary bars (all stacked segments); does not change tile or mutation legend colors. `None` retains palette defaults. |
| `gene_bar_subplot_title="", gene_bar_x_label=""` | Suppress heading and axis label independently. Any string replaces display text verbatim; `None` retains automatic defaults, including External summary for supplied aggregates. Provenance and counts remain unchanged. |
| `cell_state_palette={"unresolved": {"color": "#F5F5F5", "hatch": "///", "hatch_color": "#888888", "edgecolor": "#BBBBBB", "linewidth": .3}}` | Extend the existing state-to-color mapping with optional style dictionaries. `color` is background fill; hatch ink and border are independent. String colors retain legacy styling. Missing style fields inherit background/tile border options and the default `//` hatch. Empty hatch disables hatching, so choose a distinct fill/border if doing so. |
| `metadata_numeric_scales={"purity": {"show_missing_marker": False, "missing_color": "#EEEEEE", "colorbar": {"length": 95, "thickness": 7, "label_pad": 2}}}` | Extend the existing column-keyed numeric mapping. Marker visibility defaults to True for both numeric bars and heatmaps; original missing values remain missing. Colorbar sizes and label spacing are in typographic points, positive sizes/nonnegative spacing. Length/thickness apply to horizontal and vertical metadata colorbars. Omitted sizing uses historical defaults. |
| `legend_title_overrides={"cohort_arm": "Study\narm", "purity": "Tumor\npurity"}` | Existing source-column-keyed mapping supports multiline track labels and legend titles; no renamed input columns are necessary. |
| `tmb_height_ratio=.07, metadata_height_ratio=.14` | Existing ratios now scale Matplotlib historical panel weights relative to the main matrix, anchored at defaults `.14` and `.18`. Halving either halves the requested weight relative to the main matrix; compact mode may enlarge panels to fit measured labels, ticks and gaps. Noncompact defaults retain historical geometry. Track heights/spacing still determine internal proportions. |

Compact layout measures complete categorical legends and numeric colorbars,
including ticks and multiline labels, and shares a row when space permits.
It wraps to additional shelves when needed and rejects impossible canvas sizes
rather than clipping. Explicit colorbar sizes do not shrink automatically.
Plotly explicitly rejects state style dictionaries and nondefault gene-bar color,
gene-bar axis label, and numeric scale extensions; ordinary string state colors,
summary heading overrides, metadata display labels, and existing ratios remain
supported. No new dependencies are required.

## Matrix-aligned compact layout (issue 7)

The complete standalone [compact example](compact_alignment.py) adds counted
legend display strings, a numeric colorbar title independent of categorical
headings, explicit multiline marginal labels and a 4.33-inch figure. See
[compact options](../options-reference.md#compact-alignment-and-short-figures-matplotlib)
for automatic measured spacing, sizing limits and backend scope. No returned
axes need editing; source categories and external-cohort provenance stay intact.
