"""Publication layout tests use only synthetic, explicitly named samples."""
import importlib.util
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba, to_hex
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
import pytest
from pyoncoplot import OncoplotOptions, oncoplot

spec = importlib.util.spec_from_file_location("publication_layout", Path(__file__).parents[1] / "docs/examples/publication_layout.py")
example = importlib.util.module_from_spec(spec)
spec.loader.exec_module(example)


def axes_by_name(figure):
    return {ax.get_label(): ax for ax in figure.axes if ax.get_label()}


def test_publication_geometry_and_exports(tmp_path):
    result = example.publication_example()
    fig = result.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    axes = axes_by_name(fig)
    main, tmb, metadata, bars = (axes[k] for k in ("main", "tmb", "metadata", "gene_bar"))
    # All sample centers and edges must be identical in display coordinates.
    for axis in (tmb, metadata):
        np.testing.assert_allclose(axis.transData.transform([(0, 0), (12, 0)])[:, 0],
                                   main.transData.transform([(0, 0), (12, 0)])[:, 0], atol=1e-8)
    np.testing.assert_allclose([main.get_position().y0, main.get_position().y1],
                               [bars.get_position().y0, bars.get_position().y1])
    assert main.get_position().width > .60
    assert main.get_position().height > metadata.get_position().height * 4
    assert not main.get_xticklabels()
    assert [t.get_text() for t in metadata.get_xticklabels()] == result.prepared_data.samples
    labels = [t.get_window_extent(renderer) for t in metadata.get_xticklabels()]
    assert max(b.y1 for b in labels) < metadata.bbox.y0
    assert all(not a.overlaps(b) for a, b in zip(labels, labels[1:]))
    legend_boxes = [l.get_window_extent(renderer) for l in fig.legends]
    colorbars = [ax for ax in fig.axes if ax not in axes.values()]
    assert len(colorbars) == 1
    colorbar = colorbars[0]
    assert colorbar.get_position().width > colorbar.get_position().height * 3
    np.testing.assert_allclose(colorbar.get_xticks(), [0, .5, 1])
    decoration_boxes = legend_boxes + [colorbar.get_tightbbox(renderer)]
    assert all(not a.overlaps(b) for i, a in enumerate(decoration_boxes) for b in decoration_boxes[i+1:])
    assert max(b.y1 for b in decoration_boxes) < min(b.y0 for b in labels)
    # No clipping in the original figure canvas, even without bbox_inches='tight'.
    for box in decoration_boxes + labels + [ax.get_tightbbox(renderer) for ax in axes.values()]:
        assert box.x0 >= -1 and box.y0 >= -1
        assert box.x1 <= fig.bbox.x1 + 1 and box.y1 <= fig.bbox.y1 + 1
    title_box = fig._suptitle.get_window_extent(renderer)
    assert title_box.y0 - tmb.bbox.y1 == pytest.approx(24 / 72 * fig.dpi)
    assert tmb.get_ylabel() == "Recorded\nevents"
    assert all(p.get_width() == pytest.approx(.72) for p in tmb.patches)
    assert all(to_rgba(p.get_facecolor()) == to_rgba("#476B8A") for p in tmb.patches)
    assert all(float(t).is_integer() for t in tmb.get_yticks())
    assert [t.get_text() for t in bars.texts][:2] == ["2 (17%)", "2 (17%)"]
    assert all(t.get_window_extent(renderer).x1 <= bars.bbox.x1 for t in bars.texts)
    assert any(p.get_hatch() == "//" for p in main.patches)
    assert any(t.get_text() == "D" for l in fig.legends for t in l.get_texts())  # no truncation to three levels
    assert all(l.get_title().get_text() == "" for l in fig.legends)
    missing = [p for p in metadata.patches if p.get_x() == 3 and p.get_y() > .3][0]
    assert missing.get_facecolor() == to_rgba("#CCCCCC")
    assert metadata.patches[0].get_height() == pytest.approx(.24)
    assert missing.get_height() == pytest.approx(.32)
    assert missing.get_y() == pytest.approx(.32)
    for suffix in ("png", "pdf", "svg"):
        path = tmp_path / f"publication.{suffix}"
        result.save(path)
        assert path.stat().st_size > 1000
    plt.close(fig)


def small(**kwargs):
    data = pd.DataFrame({"gene": ["A", "A", "B"], "sample": ["S2", "S2", "S1"], "type": ["x"] * 3})
    return oncoplot(data, gene_col="gene", sample_col="sample", mutation_type_col="type",
                    include_genes=["EMPTY", "B", "A"], preserve_empty_genes=True,
                    sample_order=["S2", "S1"], draw_gene_bar=True, draw_tmb_bar=True,
                    metadata=pd.DataFrame({"sample": ["S2", "S1"], "value": [3., 7.]}),
                    metadata_sample_col="sample", metadata_cols=["value"], **kwargs)


def test_scientific_invariance():
    plain = small(backend="matplotlib")
    styled = small(backend="matplotlib", options={"compact_layout": True, "metadata_legend_position": "bottom",
                   "metadata_numeric_scales": {"value": {"min": 0, "max": 10}}, "show_sample_ids": True,
                   "gene_bar_label_format": "count_percent", "show_gene_bar_labels": True,
                   "metadata_track_heights": {"value": .25}, "row_separator_linewidth": 0})
    a, b = plain.prepared_data, styled.prepared_data
    assert a.samples == b.samples == ["S2", "S1"]
    assert a.genes == b.genes == ["EMPTY", "B", "A"]
    for name in ("tiles", "mutation_counts", "tmb", "metadata", "main_grid_rows"):
        assert_frame_equal(getattr(a, name), getattr(b, name))
    for name in ("main", "tmb", "gene_bar"):
        p = axes_by_name(plain.figure)[name].patches
        q = axes_by_name(styled.figure)[name].patches
        assert [(x.get_x(), x.get_y(), x.get_width(), x.get_height()) for x in p] == [
            (x.get_x(), x.get_y(), x.get_width(), x.get_height()) for x in q]
    plt.close(plain.figure)
    plt.close(styled.figure)


@pytest.mark.parametrize("options", [
    {"tmb_bar_width": 0}, {"tmb_bar_width": 1.1}, {"tmb_bar_width": np.nan},
    {"tmb_tick_interval": .5}, {"tmb_tick_interval": 0}, {"sample_id_padding": -1},
    {"title_gap_points": np.inf}, {"metadata_track_heights": {"value": 0}},
    {"metadata_track_spacing": np.nan}, {"tile_linewidth": -1}, {"row_separator_linewidth": -1},
    {"tile_edgecolor": "not-a-color"}, {"gene_bar_label_format": "both"},
    {"metadata_numeric_scales": {"value": {"min": 1, "max": 0}}},
    {"metadata_numeric_scales": {"value": {"ticks": [0, np.nan]}}},
    {"metadata_numeric_scales": {"value": {"ticks": [1, 0]}}},
    {"metadata_numeric_scales": {"value": {"missing_color": "bad"}}},
    {"metadata_numeric_scales": {"value": {"typo": 1}}},
    {"legend_label_overrides": {"x": 4}},
])
def test_invalid_options(options):
    with pytest.raises(ValueError):
        OncoplotOptions(**options)


@pytest.mark.parametrize("options", [
    {"metadata_track_heights": {"unknown": 1}},
    {"metadata_numeric_scales": {"value": {"min": 10}}},
    {"metadata_numeric_scales": {"value": {"min": 0, "max": 10, "ticks": [0, 11]}}},
    {"compact_layout": True, "legend_offsets": {"mutation": {"x": .1}}},
    {"compact_layout": True},
    {"compact_layout": True, "metadata_legend_position": "bottom", "metadata_legend_ncol": 2},
    {"compact_layout": True, "metadata_legend_position": "bottom", "font_size_legend_text": 200},
])
def test_invalid_effective_layout(options):
    with pytest.raises(ValueError):
        small(backend="matplotlib", options=options)
    plt.close("all")


def test_plotly_parity_and_explicit_limitations():
    result = small(backend="plotly", options={"legend_label_overrides": {"x": "Display X"},
                   "legend_title_overrides": {"value": "Display value"}})
    assert "Display X" in result.figure.to_json()
    assert "Display value" in result.figure.to_json()
    for options in ({"compact_layout": True}, {"tmb_bar_width": .5}, {"tile_edgecolor": "red"}, {"row_separator_color": "blue"},
                    {"metadata_numeric_scales": {"value": {"min": 0}}},
                    {"gene_bar_label_format": "count_percent"}):
        with pytest.raises(NotImplementedError, match="matplotlib"):
            small(backend="plotly", options=options)


def test_count_labels_without_denominators_and_independent_grid():
    data = pd.DataFrame({"gene": ["A"], "sample": ["S1"], "type": ["x"]})
    states = pd.DataFrame({"gene": ["A"], "sample": ["S2"], "state": ["unresolved"]})
    result = oncoplot(data, gene_col="gene", sample_col="sample", mutation_type_col="type",
                      cell_states=states, draw_gene_bar=True, backend="matplotlib",
                      options={"show_gene_bar_labels": True, "gene_bar_label_format": "count",
                               "tile_edgecolor": "red", "tile_linewidth": .2,
                               "row_separator_color": "blue", "row_separator_linewidth": 0})
    axes = axes_by_name(result.figure)
    assert [t.get_text() for t in axes["gene_bar"].texts] == ["1"]
    assert all(p.get_edgecolor() == to_rgba("red") for p in axes["main"].patches)
    assert all(p.get_linewidth() == .2 for p in axes["main"].patches)
    assert all(all(w == 0 for w in c.get_linewidths()) for c in axes["main"].collections)
    assert result.prepared_data.row_summary.Denominator.isna().all()
    plt.close(result.figure)


@pytest.mark.parametrize("plot_type", ["bar", "heatmap"])
def test_numeric_scale_endpoints_ticks_and_missing(plot_type):
    data = pd.DataFrame({"gene": ["A"], "sample": ["S1"], "type": ["x"]})
    metadata = pd.DataFrame({"sample": ["S1", "S2", "S3", "S4"], "value": [-1., 2., np.nan, .5]})
    result = oncoplot(data, gene_col="gene", sample_col="sample", mutation_type_col="type",
                      metadata=metadata, metadata_sample_col="sample", metadata_cols=["value"],
                      metadata_require_mutations=False, show_all_samples=True, backend="matplotlib",
                      metadata_palette={"value": "viridis"},
                      options={"metadata_numeric_plot_type": plot_type, "metadata_track_heights": {"value": .3},
                               "metadata_numeric_scales": {"value": {"min": 0, "max": 1, "ticks": [0, .25, 1], "missing_color": "pink"}}})
    fig = result.figure
    metadata_axis = axes_by_name(fig)["metadata"]
    colorbar = next(ax for ax in fig.axes if ax.get_label() not in {"metadata", "main"})
    np.testing.assert_allclose(colorbar.get_yticks(), [0, .25, 1])
    patches = metadata_axis.patches[1:] if plot_type == "bar" else metadata_axis.patches
    from matplotlib import colormaps
    assert patches[0].get_facecolor() == to_rgba(to_hex(colormaps["viridis"](0.)))
    assert patches[1].get_facecolor() == to_rgba(to_hex(colormaps["viridis"](1.)))
    assert patches[2].get_facecolor() == to_rgba("pink")
    assert all(0 <= p.get_height() <= .3 for p in patches)
    assert result.prepared_data.metadata.value.dropna().tolist() == [-1., 2., .5]
    plt.close(fig)
