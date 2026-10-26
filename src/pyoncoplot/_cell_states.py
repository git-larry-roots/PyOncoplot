"""Opt-in non-event observations and independent row summaries.

These tables never enter event collapse, ranking, TMB, or mutation-type counts.
"""
from collections.abc import Mapping

import numpy as np
import pandas as pd


def validate_cell_states(states, palette, data, gene_col, sample_col):
    from ._data import _validate_id_column

    colors = {"unresolved": "#E69F00", "not_assessed": "#999999"}
    if palette is not None:
        if not isinstance(palette, Mapping):
            raise TypeError("cell_state_palette must be a mapping of state to color.")
        colors.update(palette)
    if states is None:
        return pd.DataFrame(columns=["Gene", "Sample", "State"]), colors
    if not isinstance(states, pd.DataFrame):
        raise TypeError("cell_states must be a pandas DataFrame.")
    for col in (gene_col, sample_col, "state"):
        _validate_id_column(states, col, "Cell state")
    out = states[[gene_col, sample_col, "state"]].copy().astype(str)
    out.columns = ["Gene", "Sample", "State"]
    if out.duplicated(["Gene", "Sample"]).any():
        raise ValueError("cell_states must have one state per gene/sample pair.")
    if not set(out.State).issubset(colors):
        raise ValueError("cell_state_palette must define every supplied state.")
    events = data[[gene_col, sample_col]].astype(str).drop_duplicates()
    events.columns = ["Gene", "Sample"]
    if not out.merge(events, on=["Gene", "Sample"]).empty:
        raise ValueError("A cell cannot contain both a mutation and a non-mutation state.")
    return out, colors


def prepare_row_summary(external, denominators, gene_col, genes, tiles, states_supplied):
    from ._data import _validate_id_column, check_valid_dataframe_columns

    if external is None and not states_supplied and denominators is None:
        return None  # Preserve legacy percentage/composition semantics exactly.
    if external is not None:
        if denominators is not None:
            raise ValueError("Use row_summaries denominator column, not gene_denominators, for external summaries.")
        if not isinstance(external, pd.DataFrame):
            raise TypeError("row_summaries must be a pandas DataFrame.")
        _validate_id_column(external, gene_col, "Row summary")
        check_valid_dataframe_columns(external, ["count"])
        out = external.rename(columns={gene_col: "Gene", "count": "Count", "denominator": "Denominator"}).copy()
        out["Gene"] = out["Gene"].astype(str)
        if out.Gene.duplicated().any():
            raise ValueError("row_summaries identifiers must be unique.")
        if not set(genes).issubset(set(out.Gene)):
            raise ValueError("row_summaries must cover every displayed gene; missing summaries are not zero.")
        if "Denominator" not in out:
            out["Denominator"] = np.nan
        out = out.set_index("Gene").reindex(genes)[["Count", "Denominator"]]
    else:
        out = tiles.groupby("Gene", observed=False).size().reindex(genes, fill_value=0).to_frame("Count")
        if denominators is not None and not isinstance(denominators, Mapping):
            raise TypeError("gene_denominators must be a mapping by gene identifier.")
        out["Denominator"] = pd.Series({str(k): v for k, v in (denominators or {}).items()}, dtype=object).reindex(genes)
    for column in ("Count", "Denominator"):
        original = out[column]
        numeric = pd.to_numeric(original, errors="coerce")
        bad = original.notna() & (~np.isfinite(numeric) | numeric.isna())
        if bad.any() or (numeric.dropna() < (1 if column == "Denominator" else 0)).any():
            raise ValueError("Summary counts must be finite/nonnegative and denominators finite/positive.")
        if column == "Count" and (numeric.isna().any() or (numeric % 1 != 0).any()):
            raise ValueError("Summary counts must be nonnegative integers.")
        out[column] = numeric.astype(float)
    if (out.Count > out.Denominator).any():
        raise ValueError("Summary count cannot exceed its declared denominator.")
    out.index.name = "Gene"
    return out


def state_grid(prepared):
    """Expand state observations across tracks without fabricating mutation tiles."""
    states = prepared.cell_states
    if states is None or states.empty:
        return pd.DataFrame()
    return states.merge(prepared.main_grid_rows[["Gene", "RowIndex"]], on="Gene")


def summary_groups(prepared, mutation_levels):
    """Return actual counts, per-row denominators, and stacked bar counts."""
    tiles = prepared.tiles
    totals = tiles.groupby("Gene", observed=False).size().reindex(prepared.genes, fill_value=0)
    groups = [(level, tiles[tiles.MutationType.astype(str) == level]) for level in mutation_levels]
    unspecified = tiles[tiles.MutationType.isna()]
    if not unspecified.empty:
        groups.append((np.nan, unspecified))
    groups = [(level, group.groupby("Gene", observed=False).size().reindex(prepared.genes, fill_value=0))
              for level, group in groups]
    denominators = pd.Series(max(prepared.total_samples, 1), index=prepared.genes)
    if prepared.row_summary is not None:
        denominators = prepared.row_summary.Denominator
        if prepared.row_summary_external:
            totals = prepared.row_summary.Count
            groups = [("External summary", totals)]
    return totals, denominators, groups


def configure_summary_rendering(prepared, options, draw_gene_bar):
    """Apply the same denominator guard and source label to both entrypoints."""
    if not draw_gene_bar or prepared.row_summary is None:
        return options
    if options.show_gene_bar_labels or options.gene_bar_mode == "percent":
        if prepared.row_summary.Denominator.isna().any():
            raise ValueError("Percentage summaries require explicit denominators for every displayed gene.")
    if prepared.row_summary_external:
        from dataclasses import replace
        options = replace(options)
        title = options.gene_bar_subplot_title
        options.gene_bar_subplot_title = "External summary" + (f": {title}" if title else "")
    return options
