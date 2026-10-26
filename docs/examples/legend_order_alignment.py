"""Public API example: independent legend order, mixed heights and tight gaps.

Run with PYTHONPATH=src .venv/bin/python docs/examples/legend_order_alignment.py.
All data are synthetic; exports are written beneath .venv/issue-15-layout/.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import pandas as pd
from pyoncoplot import oncoplot


def publication_example(*, overrides=None):
    samples = [f"S{i}" for i in range(12)]
    events = pd.DataFrame({"gene": ["G1", "G1", "G2", "G3"],
                           "sample": samples[:4], "event": ["Missense"] * 4})
    metadata = pd.DataFrame({"sample": samples,
                             "treatment": list("ABCABCABCABC"),
                             "response": ["Stable", "Response"] * 6,
                             "measurement": [i / 11 for i in range(12)]})
    options = dict(compact_layout=True, width=1800, height=850, font_family="DejaVu Sans",
                   metadata_legend_position="bottom", show_sample_ids=True,
                   metadata_legend_order=["treatment", "measurement", "response"],
                   metadata_legend_alignment={key: "top" for key in ("treatment", "response", "measurement")},
                   metadata_legend_ncol={"treatment": 1, "response": 2},
                   legend_title_overrides={"treatment": "Study\ntreatment", "measurement": "Observed\nmeasurement"},
                   legend_gap_points=1, show_tmb_y_label=True, tmb_y_label="Recorded\nevents",
                   tmb_y_label_angle=90, tmb_height_ratio=.025, tmb_tick_interval=1,
                   log10_transform_tmb=False, metadata_height_ratio=.09,
                   title_text="Synthetic cohort — independent legend layout")
    options.update(overrides or {})
    return oncoplot(events, gene_col="gene", sample_col="sample", mutation_type_col="event",
                    metadata=metadata, metadata_sample_col="sample",
                    metadata_cols=["response", "measurement", "treatment"],
                    metadata_require_mutations=False, sample_order=samples, show_all_samples=True,
                    backend="matplotlib", draw_tmb_bar=True, options=options)


def main():
    output = Path(".venv/issue-15-layout")
    output.mkdir(parents=True, exist_ok=True)
    for angle in (0, 90):
        result = publication_example(overrides={"tmb_y_label_angle": angle})
        for suffix in ("png", "pdf", "svg"):
            path = output / f"legend-layout-{angle}.{suffix}"
            result.save(path, dpi=150, bbox_inches=None)
            print(path)


if __name__ == "__main__":
    main()
