# Final compact layout

Run the standalone [synthetic example](final_compact_layout.py):

```sh
PYTHONPATH=src .venv/bin/python docs/examples/final_compact_layout.py
```


It writes PNG, PDF and SVG to .venv/issue-9-compact/, using only public API
options (no returned-axes edits). All events and external summaries are synthetic.

Metadata tracks retain the exact metadata_cols order, including numeric tracks.
Compact categorical legends retain that relative order, read left-to-right then
top-to-bottom when wrapping. Numeric legends form a separate right-aligned group
in their relative metadata_cols order; their titles, ticks and scales remain tied
to the source column. A numeric legend is therefore not inserted between two
categorical groups even when its track is between them. Measured full legend
extents are vertically centered on a common shelf when they fit; otherwise shelves
wrap without truncation or font shrinking.

Reuse tmb_height_ratio, metadata_height_ratio, metadata_track_heights and
metadata_track_spacing to request short panels. Compact multiline histogram labels
are horizontal beside the axis, avoiding an unnecessarily tall rotated label.
The renderer reserves the physical clearance needed by ticks and multiline track
labels; impossible sizes still fail explicitly rather than clip. Track weights are
relative, not promises of absolute pixel heights.

The only new option is gene_bar_label_position="column" (Matplotlib only).
Its default "bar_end" preserves legacy placement. The shared column starts beyond
the longest rendered bar, with existing padding and nudge applied uniformly.
Count/percent/count_percent formatting is unchanged: computed labels use observed
events and the prepared sample denominator; external labels use the supplied counts
and explicit per-gene denominators. Placement never recomputes events, normalizes
external summaries, or changes provenance. Existing percent bar semantics are
unchanged (computed mutation composition vs external count/denominator).

This example deliberately requests non-alphabetic metadata order, several
categorical tracks, a numeric track between them, separate categorical-title
visibility and numeric show_title=True, gray hatch cell states, thin tracks,
multiline labels and an independent cohort summary. At width=2400 it fits one
bottom shelf; width=1600 demonstrates measured wrapping at the same font sizes.
