# Compact legend columns and histogram label rotation

Run the standalone [synthetic public-API example](legend_columns_rotation.py):

```bash
PYTHONPATH=src .venv/bin/python docs/examples/legend_columns_rotation.py
```

It writes PNG, PDF and SVG for both 0° and 90° labels to
`.venv/issue-11-columns-rotation/`, on the original canvas (no tight cropping).
Three categorical metadata groups use `metadata_legend_ncol=2`, including five
study arms split into columns of three and two entries. The numeric purity
colorbar is measured and right-aligned automatically. Categories read down each
column, then right, matching standard Matplotlib legends.

Use `metadata_legend_nrow` instead to cap rows. With both options supplied the
effective column count is `max(ncol, ceil(entries / nrow))`. Overwide explicit
columns fail descriptively; increase width or request fewer columns. Omitting
both retains automatic layout.

`tmb_y_label_angle=90` gives the standard rotated histogram y-label; `0` makes
it horizontal. Explicit `\n` line breaks and `font_size_tmb_axis` are preserved.
The default `None` retains existing automatic behavior: horizontal multiline
labels in compact mode, otherwise standard 90° labels. Angles are finite degrees
using Matplotlib rotation conventions, like `sample_id_angle`. A rotated label
may need a taller histogram; compact layout reserves measured text space rather
than shrinking fonts. Plotly rejects explicit angles with `NotImplementedError`.

The example includes explicitly ordered samples, empty genes, unresolved and
not-assessed cell states, and external counts with independent denominators.
These styling controls do not alter assignments, counts, or denominators.
