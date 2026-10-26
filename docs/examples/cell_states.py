"""Run from the repository root: PYTHONPATH=src python docs/examples/cell_states.py.

Writes a static PNG and interactive HTML under .venv/ (review artifacts only).
"""
from pathlib import Path
import pandas as pd
from pyoncoplot import oncoplot


def main():
    events = pd.DataFrame({
        "gene": ["TP53", "KRAS"], "sample": ["S1", "S2"],
        "mutation_type": ["Missense_Mutation", "Missense_Mutation"],
    })
    states = pd.DataFrame({
        "gene": ["TP53", "EMPTY"], "sample": ["S2", "S3"],
        "state": ["unresolved", "not_assessed"],
    })
    # Deliberately out of display order; this aggregate cohort is independent.
    summaries = pd.DataFrame({
        "gene": ["KRAS", "EMPTY", "TP53"], "count": [3, 2, 7],
        "denominator": [10, 8, 20],
    })
    output = Path(".venv/issue-1-smoke")
    output.mkdir(parents=True, exist_ok=True)
    for backend, suffix in [("matplotlib", "png"), ("plotly", "html")]:
        result = oncoplot(
            events, gene_col="gene", sample_col="sample", mutation_type_col="mutation_type",
            include_genes=["EMPTY", "TP53", "KRAS"], preserve_empty_genes=True,
            cell_states=states, row_summaries=summaries, show_all_samples=True,
            backend=backend, draw_gene_bar=True, draw_tmb_bar=True,
            options={"gene_bar_mode": "percent", "show_gene_bar_labels": True,
                     "title_text": "Recorded events, non-event states, independent aggregates",
                     "log10_transform_tmb": False, "width": 1200, "height": 650},
        )
        result.save(output / f"cell-states.{suffix}")
        print(output / f"cell-states.{suffix}")


if __name__ == "__main__":
    main()
