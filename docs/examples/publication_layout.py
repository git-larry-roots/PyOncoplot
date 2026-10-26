"""Synthetic publication layout: run PYTHONPATH=src .venv/bin/python docs/examples/publication_layout.py.

Artifacts go to .venv/issue-2-publication/ and contain no patient data.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import pandas as pd
from pyoncoplot import oncoplot


def publication_example():
    samples = [f"S{i:02}" for i in range(1, 13)]
    events = pd.DataFrame({
        "gene": ["TP53", "TP53", "KRAS", "KRAS", "PIK3CA", "PTEN", "RB1"],
        "sample": [samples[i] for i in [0, 3, 1, 4, 2, 7, 9]],
        "event": ["Missense_Mutation"] * 5 + ["Frame_Shift_Del", "Nonsense_Mutation"],
    })
    metadata = pd.DataFrame({"sample": samples,
                             "cohort_arm": ["A", "B", "C", "D"] * 3,
                             "purity": [.2, .3, .5, None, .7, .6, .8, .4, .9, .5, .6, .7]})
    genes = ["TP53", "KRAS", "PIK3CA", "PTEN", "RB1", "EMPTY"]
    states = pd.DataFrame({"gene": ["TP53", "EMPTY"], "sample": ["S02", "S12"],
                           "state": ["unresolved", "not_assessed"]})
    result = oncoplot(
        events, gene_col="gene", sample_col="sample", mutation_type_col="event",
        include_genes=genes, preserve_empty_genes=True, cell_states=states,
        gene_denominators={g: 12 for g in genes}, metadata=metadata, metadata_sample_col="sample",
        metadata_cols=["cohort_arm", "purity"], metadata_require_mutations=False,
        sample_order=samples, show_all_samples=True,
        backend="matplotlib", draw_gene_bar=True, draw_tmb_bar=True,
        options={
            "compact_layout": True, "width": 1200, "height": 800, "font_family": "DejaVu Sans",
            "metadata_position": "bottom", "metadata_legend_position": "bottom",
            "metadata_legend_orientation_heatmap": "horizontal", "show_sample_ids": True,
            "sample_id_padding": 6, "sample_id_angle": 90,
            "log10_transform_tmb": False, "tmb_bar_color": "#476B8A", "tmb_bar_width": .72,
            "tmb_tick_interval": 1, "show_tmb_y_label": True, "tmb_y_label": "Recorded\nevents",
            "show_gene_bar_labels": True, "gene_bar_label_format": "count_percent",
            "gene_bar_label_padding": .08, "gene_bar_width_ratio": .17,
            "metadata_track_heights": {"cohort_arm": .24, "purity": .32}, "metadata_track_spacing": .08,
            "metadata_numeric_scales": {"purity": {"min": 0, "max": 1, "ticks": [0, .5, 1], "missing_color": "#CCCCCC"}},
            "tile_linewidth": .3, "tile_edgecolor": "white", "row_separator_linewidth": 0,
            "title_text": "Synthetic cohort — publication layout", "title_gap_points": 24,
            "legend_label_overrides": {"Missense_Mutation": "Missense", "Frame_Shift_Del": "Frameshift", "Nonsense_Mutation": "Nonsense"},
            "legend_title_overrides": {"cohort_arm": "Arm", "purity": "Purity"}, "show_legend_titles": False,
        },
    )
    return result


def main():
    output = Path(".venv/issue-2-publication")
    output.mkdir(parents=True, exist_ok=True)
    result = publication_example()
    for suffix in ("png", "pdf", "svg"):
        path = output / f"publication-layout.{suffix}"
        result.save(path, dpi=150)
        print(path)


if __name__ == "__main__":
    main()
