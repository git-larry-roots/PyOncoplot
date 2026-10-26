"""Synthetic acceptance, invalid-input, and default regressions for issue #1."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal, assert_series_equal

from pyoncoplot import oncoplot, prepare_oncoplot_data


def events():
    return pd.DataFrame({"gene": ["A", "A", "B"], "sample": ["S1", "S1", "S2"],
                         "type": ["Missense_Mutation"] * 3, "vaf": [.2, .4, .3]})


def states():
    return pd.DataFrame({"gene": ["A", "EMPTY"], "sample": ["S2", "S3"],
                         "state": ["unresolved", "not_assessed"]})


def params():
    return dict(gene_col="gene", sample_col="sample", mutation_type_col="type",
                include_genes=["EMPTY", "A", "B"], preserve_empty_genes=True,
                show_all_samples=True)


def test_states_do_not_enter_event_summaries():
    baseline = prepare_oncoplot_data(events(), **params(), prepare_tmb=True)
    marked = prepare_oncoplot_data(events(), **params(), cell_states=states(), prepare_tmb=True)
    assert marked.genes == ["EMPTY", "A", "B"]
    assert marked.samples == ["S1", "S2", "S3"]
    assert_frame_equal(baseline.tiles, marked.tiles, check_categorical=False)
    assert_frame_equal(baseline.mutation_counts, marked.mutation_counts)
    assert_series_equal(baseline.tmb_totals, marked.tmb_totals.reindex(baseline.samples))
    assert marked.tmb_totals["S3"] == 0
    assert "unresolved" not in marked.mutation_type_levels
    assert marked.row_summary.Denominator.isna().all()
    assert marked.tiles.MutationCount.sum() == 3


@pytest.mark.parametrize("backend", ["plotly", "matplotlib"])
@pytest.mark.parametrize("expanded", [False, True])
def test_render_acceptance(backend, expanded):
    external = pd.DataFrame({"gene": ["B", "A", "EMPTY"], "count": [3, 7, 2], "denominator": [10, 20, 8]})
    extra = {"variant_value_cols": ["vaf"]} if expanded else {}
    result = oncoplot(events(), **params(), cell_states=states(), row_summaries=external,
                      backend=backend, draw_gene_bar=True, draw_tmb_bar=True,
                      options={"gene_bar_mode": "percent", "show_gene_bar_labels": True}, **extra)
    prepared = result.prepared_data
    assert prepared.row_summary.Count.tolist() == [2, 7, 3]
    assert prepared.row_summary.Denominator.tolist() == [8, 20, 10]
    assert len(prepared.tiles) == 2
    if backend == "plotly":
        bars = [t for t in result.figure.data if t.type == "bar" and t.orientation == "h"]
        assert len(bars) == 1
        assert list(bars[0].x) == [25, 35, 30]
        assert any(t.name == "Cell state: unresolved" and t.showlegend for t in result.figure.data)
        labels = [text for t in result.figure.data if t.type == "scatter" and t.mode == "text" for text in t.text]
        assert set(["25%", "35%", "30%"]).issubset(labels)
        assert "External summary" in result.figure.to_json()
    else:
        bar_axis = next(ax for ax in result.figure.axes if "External summary" in ax.get_xlabel())
        assert [p.get_width() for p in bar_axis.patches] == [25, 35, 30]
        assert [t.get_text() for t in bar_axis.texts] == ["25%", "35%", "30%"]
        assert any(p.get_hatch() == "//" for ax in result.figure.axes for p in ax.patches)
        assert any("Cell state: unresolved" == t.get_text() for legend in result.figure.legends for t in legend.get_texts())
        plt.close(result.figure)


@pytest.mark.parametrize("backend", ["plotly", "matplotlib"])
def test_states_require_explicit_denominators_only_for_percentages(backend):
    for options in ({"show_gene_bar_labels": True}, {"gene_bar_mode": "percent"}):
        with pytest.raises(ValueError, match="explicit denominators"):
            oncoplot(events(), **params(), cell_states=states(), draw_gene_bar=True, backend=backend, options=options)
    result = oncoplot(events(), **params(), cell_states=states(), draw_gene_bar=True, backend=backend)
    if backend == "matplotlib":
        plt.close(result.figure)
    result = oncoplot(events(), **params(), cell_states=states(), gene_denominators={"A": 1, "B": 2, "EMPTY": 1},
                      draw_gene_bar=True, backend=backend, options={"gene_bar_mode": "percent"})
    assert result.prepared_data.row_summary.Count.tolist() == [0, 1, 1]
    if backend == "matplotlib":
        plt.close(result.figure)


@pytest.mark.parametrize("bad", [
    pd.DataFrame({"gene": ["A"], "sample": ["S1"], "state": ["unresolved"]}),
    pd.concat([states(), states()]),
    states().assign(state="unknown"),
    states().assign(sample=None),
])
def test_invalid_states(bad):
    with pytest.raises(ValueError):
        prepare_oncoplot_data(events(), **params(), cell_states=bad)


@pytest.mark.parametrize("counts,denominators", [([-1, 0, 0], [1, 1, 1]), ([.5, 0, 0], [1, 1, 1]),
    ([np.inf, 0, 0], [1, 1, 1]), ([2, 0, 0], [1, 1, 1]), ([0, 0, 0], [0, 1, 1]),
    ([0, 0, 0], [np.inf, 1, 1]), ([None, 0, 0], [1, 1, 1])])
def test_invalid_external_values(counts, denominators):
    summary = pd.DataFrame({"gene": ["A", "B", "EMPTY"], "count": counts, "denominator": denominators})
    with pytest.raises(ValueError):
        prepare_oncoplot_data(events(), **params(), row_summaries=summary)


def test_external_missing_duplicate_and_count_only():
    for genes in (["A", "B"], ["A", "B", "EMPTY", "A"]):
        summary = pd.DataFrame({"gene": genes, "count": 0})
        with pytest.raises(ValueError):
            prepare_oncoplot_data(events(), **params(), row_summaries=summary)
    summary = pd.DataFrame({"gene": ["A", "B", "EMPTY"], "count": [7, 2, 0]})
    result = oncoplot(events(), **params(), row_summaries=summary, draw_gene_bar=True)
    assert result.prepared_data.row_summary.Count.tolist() == [0, 7, 2]
    with pytest.raises(ValueError, match="explicit denominators"):
        oncoplot(events(), **params(), row_summaries=summary, draw_gene_bar=True, options={"show_gene_bar_labels": True})


def test_default_regression_and_preserved_order_over_pathways():
    legacy = prepare_oncoplot_data(events(), gene_col="gene", sample_col="sample", include_genes=["EMPTY", "B", "A"])
    assert legacy.genes == ["B", "A"]
    assert legacy.row_summary is None
    pathway = pd.DataFrame({"gene": ["A", "B", "EMPTY"], "pathway": ["X", "Y", "X"]})
    preserved = prepare_oncoplot_data(events(), **params(), pathway=pathway)
    assert preserved.genes == ["EMPTY", "A", "B"]


def test_yaml_table_loading(tmp_path):
    import yaml
    events().to_csv(tmp_path / "events.tsv", sep="\t", index=False)
    states().to_csv(tmp_path / "states.csv", index=False)
    pd.DataFrame({"gene": ["A", "B", "EMPTY"], "count": [1, 2, 0], "denominator": [5, 5, 5]}).to_csv(tmp_path / "summary.csv", index=False)
    config = {**params(), "data": "events.tsv", "cell_states": {"path": "states.csv"}, "row_summaries": "summary.csv", "draw_gene_bar": True}
    (tmp_path / "params.yaml").write_text(yaml.safe_dump(config))
    result = oncoplot(params=tmp_path / "params.yaml")
    assert result.prepared_data.row_summary.Count.tolist() == [0, 1, 2]


@pytest.mark.parametrize("backend", ["plotly", "matplotlib"])
def test_all_empty_rows_with_state_only_cohort(backend):
    empty = events().iloc[:0]
    result = oncoplot(empty, **params(), cell_states=states(), backend=backend,
                      draw_gene_bar=True, draw_tmb_bar=True,
                      options={"log10_transform_tmb": False})
    assert result.prepared_data.genes == ["EMPTY", "A", "B"]
    assert result.prepared_data.tiles.empty
    assert result.prepared_data.tmb_totals.sum() == 0
    if backend == "matplotlib":
        plt.close(result.figure)


def test_counts_hover_never_infers_unresolved_denominator():
    result = oncoplot(events(), **params(), cell_states=states(), draw_gene_bar=True)
    bars = [t for t in result.figure.data if t.type == "bar" and t.orientation == "h"]
    assert bars
    for bar in bars:
        assert all("percentage unavailable" in hover and "%" not in hover for hover in bar.hovertext)


def test_filters_custom_states_and_input_immutability():
    original_events, original_states = events(), states()
    e, st = original_events.copy(deep=True), original_states.copy(deep=True)
    st["state"] = "pending"
    prepared = prepare_oncoplot_data(e, **params(), cell_states=st,
        cell_state_palette={"pending": "purple"}, filter_samples_by_isin_lists={"sample": ["S1", "S2"]})
    assert set(prepared.cell_states.Sample) == {"S2"}
    assert_frame_equal(e, original_events)
    assert_frame_equal(original_states, states())
    assert_frame_equal(st, original_states.assign(state="pending"))


@pytest.mark.parametrize("included", [None, [], "A", ["A", "A"], [None], [""]])
def test_preserve_requires_explicit_unique_rows(included):
    with pytest.raises(ValueError):
        prepare_oncoplot_data(events(), gene_col="gene", sample_col="sample", preserve_empty_genes=True, include_genes=included)


@pytest.mark.parametrize("backend", ["plotly", "matplotlib"])
def test_direct_renderer_denominators_and_external_label(backend):
    from pyoncoplot._plotly import render_plotly_oncoplot
    from pyoncoplot._matplotlib import render_matplotlib_oncoplot
    renderer = render_plotly_oncoplot if backend == "plotly" else render_matplotlib_oncoplot
    prepared = prepare_oncoplot_data(events(), **params(), cell_states=states())
    with pytest.raises(ValueError, match="explicit denominators"):
        renderer(prepared, palette={}, draw_gene_bar=True, options={"gene_bar_mode": "percent"})
    summary = pd.DataFrame({"gene": ["A", "B", "EMPTY"], "count": [5, 2, 1]})
    prepared = prepare_oncoplot_data(events(), **params(), row_summaries=summary)
    fig = renderer(prepared, palette={}, draw_gene_bar=True, options={"show_gene_bar_axis": False})
    if backend == "matplotlib":
        assert any(ax.get_title() == "External summary" for ax in fig.axes)
        plt.close(fig)
    else:
        assert any(a.text == "External summary" for a in fig.layout.annotations)


@pytest.mark.parametrize("backend", ["plotly", "matplotlib"])
def test_preserved_rows_without_events_or_states(backend):
    from pyoncoplot import OncoplotOptions
    options = OncoplotOptions(show_gene_bar_axis=False)
    meta = pd.DataFrame({"sample": ["N1", "N2"]})
    summary = pd.DataFrame({"gene": ["X", "Y"], "count": [0, 4]})
    for _ in range(2):
        result = oncoplot(events().iloc[:0], gene_col="gene", sample_col="sample",
            include_genes=["Y", "X"], preserve_empty_genes=True,
            metadata=meta, metadata_require_mutations=False, show_all_samples=True,
            row_summaries=summary, options=options, draw_gene_bar=True, backend=backend)
        assert result.prepared_data.genes == ["Y", "X"]
        assert result.prepared_data.samples == ["N1", "N2"]
        assert result.prepared_data.tiles.empty
        assert options.gene_bar_subplot_title is None
        if backend == "matplotlib":
            plt.close(result.figure)


def test_unselected_state_rows_do_not_add_subjects():
    prepared = prepare_oncoplot_data(events(), gene_col="gene", sample_col="sample",
        include_genes=["A"], cell_states=states().iloc[[1]], show_all_samples=True)
    assert "S3" not in prepared.samples
    assert prepared.cell_states.empty
