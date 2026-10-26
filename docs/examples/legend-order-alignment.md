# Independent legend order and alignment

Run the [complete public API synthetic example](legend_order_alignment.py):


```sh
PYTHONPATH=src .venv/bin/python docs/examples/legend_order_alignment.py
```

It exports PNG, PDF and SVG files for horizontal and rotated multiline histogram
labels beneath .venv/issue-15-layout/. Metadata tracks stay response, measurement,
treatment, while bottom legends display treatment, measurement, response.


```python
options = {
    "compact_layout": True,
    "metadata_legend_position": "bottom",
    "metadata_legend_order": ["treatment", "measurement", "response"],
    "metadata_legend_alignment": {
        "treatment": "top", "measurement": "top", "response": "top",
    },
    "metadata_legend_ncol": {"treatment": 1, "response": 2},
    "legend_gap_points": 1,
    "tmb_height_ratio": .025,
    "tmb_y_label": "Recorded\nevents",
    "tmb_y_label_angle": 90,
}
```

All three full legend bounds are top-aligned despite their different heights.
Options affect layout only, not sample order, assignments, counts, denominators
or provenance. Unspecified alignment entries remain centered. Automatic defaults
retain historical placement. These controls require compact Matplotlib layout;
unsupported renderers/layouts report an explicit error rather than ignoring them.
