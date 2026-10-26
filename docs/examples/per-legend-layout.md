# Per-legend columns and short panels

Run the standalone [synthetic public-API example](per_legend_layout.py):

    PYTHONPATH=src .venv/bin/python docs/examples/per_legend_layout.py

It writes PNG, PDF and SVG to .venv/issue-13-per-legend/ without modifying returned
axes. It retains sample order, categorical order and independent row-summary
counts/denominators. The source columns deliberately differ from display titles.

The example uses metadata_legend_ncol={"z_response": 1, "cohort_arm": 2} and
metadata_legend_nrow={"a_site": 2}. The two-entry response legend is one column,
the five-entry arm legend uses two columns, and the three-entry site legend uses
two columns to limit rows to two. Entries fill down each column, then right.
An absent mapping entry stays automatic; integers still apply globally.

The histogram uses tmb_height_ratio=.015 and tmb_y_label_angle=90 for a multiline
label. Metadata uses metadata_height_ratio=.025 with thin weighted tracks.
Measured text minima keep labels readable without reducing fonts; requested
ratios are relative weights, not hard pixel heights. Rotation can impose a real
minimum height (a long vertical label cannot fit into an arbitrarily short bar).
For unrotated labels, adjacent tick spacing is measured rather than summing all
tick heights, avoiding unnecessary histogram inflation.

Bottom legends share measured shelves, wrap in caller order when needed, and
center-align within shelves; numeric bars remain right-aligned. The automatic
fitter reserves real text clearance and raises clear errors for impossible sizes.
Use larger width/height if the full requested columns cannot fit.
