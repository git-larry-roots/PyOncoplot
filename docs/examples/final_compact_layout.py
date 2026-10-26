"""Synthetic publication layout: run PYTHONPATH=src .venv/bin/python docs/examples/final_compact_layout.py.

Artifacts go to .venv/issue-9-compact/ and contain no patient data.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import pandas as pd
from pyoncoplot import oncoplot


def publication_example(*, width=2400, overrides=None):
    samples = [f"S{i:02}" for i in range(1, 13)]
    events = pd.DataFrame({
        "gene": ["TP53", "TP53", "KRAS", "KRAS", "PIK3CA", "PTEN", "RB1"],
        "sample": [samples[i] for i in [0, 3, 1, 4, 2, 7, 9]],
        "event": ["Missense_Mutation"] * 5 + ["Frame_Shift_Del", "Nonsense_Mutation"],
    })
    metadata = pd.DataFrame({"sample": samples,
                             "cohort_arm": ["A", "B", "C", "D"] * 3,
                             "z_response": ["Stable", "Responding"] * 6,
                             "a_site": ["Lung", "Liver", "Other"] * 4,
                             "purity": [.2, .3, .5, None, .7, .6, .8, .4, .9, .5, .6, .7]})
    genes = ["TP53", "KRAS", "PIK3CA", "PTEN", "RB1", "EMPTY"]
    states = pd.DataFrame({"gene": ["TP53", "EMPTY"], "sample": ["S02", "S12"],
                           "state": ["unresolved", "not_assessed"]})
    options = {
        "width": width, "gene_bar_color": "#476B8A", "gene_bar_subplot_title": "",
        "gene_bar_x_label": "Independent cohort\nAffected (%)", "gene_bar_mode": "percent",
        "tmb_height_ratio": .055, "metadata_height_ratio": .10,
        "legend_title_overrides": {"z_response": "Treatment\nresponse", "a_site": "Sample\nsite", "cohort_arm": "Study\narm", "purity": "Tumor\npurity"},
        "legend_label_overrides": {"Missense_Mutation": "Missense (n=5)", "Frame_Shift_Del": "Frameshift (n=1)",
                                   "Nonsense_Mutation": "Nonsense (n=1)", "Cell state: unresolved": "Unassigned",
                                   "Cell state: not_assessed": "Not assessed",
                                   **{arm: f"{arm} (n=3)" for arm in "ABCD"}},
        "metadata_numeric_scales": {"purity": {"min": 0, "max": 1, "ticks": [0, .5, 1],
            "missing_color": "#EEEEEE", "show_missing_marker": False, "show_title": True,
            "colorbar": {"length": 95, "thickness": 7, "label_pad": 2}}},
    }
    options.update(overrides or {})
    summaries = pd.DataFrame({"gene": genes, "count": [9, 7, 5, 4, 3, 1], "denominator": [20] * 6})
    result = oncoplot(
        events, gene_col="gene", sample_col="sample", mutation_type_col="event",
        include_genes=genes, preserve_empty_genes=True, cell_states=states,
        row_summaries=summaries,
        cell_state_palette={state: {"color": "#F5F5F5", "hatch": "///", "hatch_color": "#888888",
                                     "edgecolor": "#BBBBBB", "linewidth": .3}
                            for state in ("unresolved", "not_assessed")}, metadata=metadata, metadata_sample_col="sample",
        metadata_cols=["z_response", "purity", "cohort_arm", "a_site"], metadata_require_mutations=False,
        sample_order=samples, show_all_samples=True,
        backend="matplotlib", draw_gene_bar=True, draw_tmb_bar=True,
        options={
            "compact_layout": True, "height": 600, "font_family": "DejaVu Sans",
            "metadata_position": "bottom", "metadata_legend_position": "bottom",
            "metadata_legend_orientation_heatmap": "horizontal", "show_sample_ids": True,
            "sample_id_padding": 1, "legend_gap_points": 3, "sample_id_angle": 90,
            "log10_transform_tmb": False, "tmb_bar_color": "#476B8A", "tmb_bar_width": .72,
            "tmb_tick_interval": 1, "show_tmb_y_label": True, "tmb_y_label": "Recorded\nevents",
            "show_gene_bar_labels": True, "gene_bar_label_position": "column", "gene_bar_label_format": "count_percent",
            "gene_bar_label_padding": .08, "gene_bar_width_ratio": .17,
            "metadata_track_heights": {"cohort_arm": .24, "purity": .24, "z_response": .24, "a_site": .24}, "metadata_track_spacing": .08,
            "tile_linewidth": .3, "tile_edgecolor": "white", "row_separator_linewidth": 0,
            "title_text": "Synthetic cohort — compact alignment", "title_gap_points": 10,
            "show_legend_titles": False,
            **options,
        },
    )
    return result


def main():
    output = Path(".venv/issue-9-compact")
    output.mkdir(parents=True, exist_ok=True)
    result = publication_example()
    for suffix in ("png", "pdf", "svg"):
        path = output / f"final-compact-layout.{suffix}"
        result.save(path, dpi=150)
        print(path)


if __name__ == "__main__":
    main()
