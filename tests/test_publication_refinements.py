"""Issue 5: public controls separate scientific inputs from presentation."""
import importlib.util
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import numpy as np
import pytest
from pandas.testing import assert_frame_equal
from pyoncoplot import OncoplotOptions, oncoplot
from test_publication_layout import axes_by_name, small

spec = importlib.util.spec_from_file_location("refinements", Path(__file__).parents[1] / "docs/examples/publication_refinements.py")
example = importlib.util.module_from_spec(spec)
spec.loader.exec_module(example)


@pytest.mark.parametrize("width,shared", [(1800, True), (1200, False)])
def test_shared_band_and_wrap_geometry(width, shared, tmp_path):
    result = example.publication_example(width=width)
    fig = result.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    axes = axes_by_name(fig)
    cb = next(ax for ax in fig.axes if ax not in axes.values())
    boxes = [leg.get_window_extent(renderer) for leg in fig.legends]
    colorbox = cb.get_tightbbox(renderer)
    assert any(box.y0 < colorbox.y1 and colorbox.y0 < box.y1 for box in boxes) == shared
    boxes.append(colorbox)
    assert all(not a.overlaps(b) for i, a in enumerate(boxes) for b in boxes[i+1:])
    assert all(b.x0 >= 0 and b.y0 >= 0 and b.x1 <= fig.bbox.x1 and b.y1 <= fig.bbox.y1 for b in boxes)
    labels = [t.get_window_extent(renderer) for t in axes["metadata"].get_xticklabels()]
    assert max(b.y1 for b in boxes) < min(b.y0 for b in labels)
    for name in ("tmb", "metadata"):
        np.testing.assert_allclose(axes[name].transData.transform([(0, 0), (12, 0)])[:, 0],
                                   axes["main"].transData.transform([(0, 0), (12, 0)])[:, 0])
    assert axes["gene_bar"].get_title() == ""
    assert axes["gene_bar"].get_xlabel() == "Affected (%)"
    assert all(p.get_facecolor() == to_rgba("#476B8A") for p in axes["gene_bar"].patches)
    assert [t.get_text() for t in axes["metadata"].get_yticklabels()] == ["Study\narm", "Tumor\npurity"]
    assert not any(t.get_text() == "!" for t in axes["metadata"].texts)
    hatch = next(p for p in axes["main"].patches if p.get_hatch())
    assert hatch.get_hatch() == "///"
    assert hatch.get_edgecolor() == to_rgba("#888888")
    assert hatch.get_facecolor() == to_rgba("#F5F5F5")
    assert cb.get_position().width * fig.get_figwidth() * 72 == pytest.approx(95)
    assert cb.get_position().height * fig.get_figheight() * 72 == pytest.approx(7)
    assert cb.xaxis.labelpad == 2
    for ext in ("png", "pdf", "svg"):
        result.save(tmp_path / f"refined-{width}.{ext}")
    plt.close(fig)


def test_invariance_and_ratio_defaults():
    plain = example.publication_example(overrides={"tmb_height_ratio": .14, "metadata_height_ratio": .18})
    styled = example.publication_example()
    a, b = plain.prepared_data, styled.prepared_data
    assert a.samples == b.samples and a.genes == b.genes
    assert a.row_summary_external and b.row_summary_external
    for key in ("row_summary", "tiles", "mutation_counts", "tmb", "metadata", "cell_states", "main_grid_rows"):
        assert_frame_equal(getattr(a, key), getattr(b, key))
    p, q = axes_by_name(plain.figure), axes_by_name(styled.figure)
    assert q["main"].get_position().height > p["main"].get_position().height
    # Histogram weights also respect the measured multiline-label minimum.
    assert q["tmb"].get_position().height/q["main"].get_position().height + 1e-12 >= (
        .5*p["tmb"].get_position().height/p["main"].get_position().height)
    # Metadata is now subject to measured multiline-label minimums. Assert the
    # requested weight remains a lower bound and the complete labels do not overlap.
    assert q["metadata"].get_position().height/q["main"].get_position().height >= (
        .14/.18*p["metadata"].get_position().height/p["main"].get_position().height)
    styled.figure.canvas.draw()
    boxes = [t.get_window_extent(styled.figure.canvas.get_renderer()) for t in q["metadata"].get_yticklabels()]
    assert all(not a.overlaps(b) for i,a in enumerate(boxes) for b in boxes[i+1:])
    default = small(backend="matplotlib")
    axes = axes_by_name(default.figure)
    assert axes["tmb"].get_position().height/axes["main"].get_position().height == pytest.approx(1.2/3)
    assert axes["metadata"].get_position().height/axes["main"].get_position().height == pytest.approx(1/3)
    colored = small(backend="matplotlib", options={"gene_bar_color": "pink"})
    for pbar, qbar in zip(axes["gene_bar"].patches, axes_by_name(colored.figure)["gene_bar"].patches):
        assert pbar.get_width() == qbar.get_width()
        assert qbar.get_facecolor() == to_rgba("pink")
    plt.close("all")


@pytest.mark.parametrize("plot_type", ["bar", "heatmap"])
@pytest.mark.parametrize("show", [True, False])
def test_missing_marker_toggle(plot_type, show):
    result = example.publication_example(overrides={"metadata_numeric_plot_type": plot_type,
        "metadata_numeric_scales": {"purity": {"show_missing_marker": show}}})
    assert any(t.get_text() == "!" for t in axes_by_name(result.figure)["metadata"].texts) == show
    plt.close(result.figure)


@pytest.mark.parametrize("override", [{"gene_bar_subplot_title": "Custom", "gene_bar_x_label": "Custom axis"},
                                       {"gene_bar_subplot_title": "", "gene_bar_x_label": ""}])
def test_summary_display_override(override):
    result = example.publication_example(overrides=override)
    axes = axes_by_name(result.figure)
    assert axes["gene_bar"].get_title() == override["gene_bar_subplot_title"]
    assert axes["gene_bar"].get_xlabel() == override["gene_bar_x_label"]
    assert result.prepared_data.row_summary_external
    assert result.prepared_data.row_summary.Count.tolist() == [9, 7, 5, 4, 3, 1]
    plt.close(result.figure)


@pytest.mark.parametrize("scale", [{"show_missing_marker": "no"}, {"colorbar": {"length": 0}},
    {"colorbar": {"thickness": np.nan}}, {"colorbar": {"label_pad": -1}}, {"colorbar": {"typo": 1}}])
def test_invalid_numeric_controls(scale):
    with pytest.raises(ValueError):
        OncoplotOptions(metadata_numeric_scales={"value": scale})


@pytest.mark.parametrize("style", [{"hatch": "?"}, {"color": "bad"}, {"linewidth": -1}, {"typo": 1}])
def test_invalid_state_controls(style):
    from pyoncoplot import prepare_oncoplot_data
    import pandas as pd
    with pytest.raises(ValueError):
        prepare_oncoplot_data(pd.DataFrame({"gene": ["A"], "sample": ["S1"]}),
                              gene_col="gene", sample_col="sample", cell_state_palette={"unresolved": style})


@pytest.mark.parametrize("options", [{"gene_bar_color": "red"}, {"gene_bar_x_label": ""},
    {"metadata_numeric_scales": {"value": {"show_missing_marker": False}}}])
def test_plotly_explicit_limitations(options):
    with pytest.raises(NotImplementedError, match="matplotlib"):
        small(backend="plotly", options=options)


def test_long_multiline_titles_and_ticks():
    result = example.publication_example(width=1200, overrides={
        "show_legend_titles": True,
        "legend_title_overrides": {"purity": "A longer numeric\nmetadata heading", "cohort_arm": "Study\narm"},
        "metadata_numeric_scales": {"purity": {"min": -10000, "max": 10000,
            "ticks": [-10000, 0, 10000], "colorbar": {"length": 130, "label_pad": 8}}}})
    fig = result.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    axes = axes_by_name(fig)
    boxes = [l.get_window_extent(renderer) for l in fig.legends]
    boxes += [ax.get_tightbbox(renderer) for ax in fig.axes if ax not in axes.values()]
    assert all(not a.overlaps(b) for i, a in enumerate(boxes) for b in boxes[i+1:])
    plt.close(fig)


def test_plotly_rejects_state_style_mapping():
    from pyoncoplot import prepare_oncoplot_data
    from pyoncoplot._plotly import render_plotly_oncoplot
    import pandas as pd
    prepared = prepare_oncoplot_data(pd.DataFrame({"gene": ["A"], "sample": ["S1"]}),
        gene_col="gene", sample_col="sample", cell_state_palette={"unresolved": {"color": "gray"}})
    with pytest.raises(NotImplementedError, match="matplotlib"):
        render_plotly_oncoplot(prepared, palette={}, options={})


def test_vertical_colorbar_sizes_and_impossible_compact_length():
    result = small(backend="matplotlib", options={"show_legend_titles": True,
        "metadata_numeric_scales": {"value": {"colorbar": {"length": 75, "thickness": 5, "label_pad": 3}}}})
    axes = axes_by_name(result.figure)
    cb = next(ax for ax in result.figure.axes if ax not in axes.values())
    assert cb.get_position().height * result.figure.get_figheight() * 72 == pytest.approx(75)
    assert cb.get_position().width * result.figure.get_figwidth() * 72 == pytest.approx(5)
    plt.close(result.figure)
    with pytest.raises(ValueError, match="Colorbar labels exceed figure width"):
        example.publication_example(overrides={"metadata_numeric_scales": {"purity": {"colorbar": {"length": 5000}}}})
    plt.close("all")


def test_multiple_numeric_colorbars_pack_without_overlap():
    import pandas as pd
    result = oncoplot(pd.DataFrame({"gene": ["A"], "sample": ["S1"]}),
        gene_col="gene", sample_col="sample", backend="matplotlib",
        metadata=pd.DataFrame({"sample": ["S1", "S2"], "x": [0., 1.], "y": [3., 6.]}),
        metadata_sample_col="sample", metadata_cols=["x", "y"], metadata_require_mutations=False,
        show_all_samples=True, options={"compact_layout": True, "metadata_legend_position": "bottom",
            "legend_title_overrides": {"x": "First\nmeasurement", "y": "Second\nmeasurement"}})
    fig = result.figure
    fig.canvas.draw()
    axes = axes_by_name(fig)
    cb = [ax for ax in fig.axes if ax not in axes.values()]
    assert len(cb) == 2
    boxes = [ax.get_tightbbox(fig.canvas.get_renderer()) for ax in cb]
    assert not boxes[0].overlaps(boxes[1])
    plt.close(fig)


def test_new_options_preserve_existing_positional_field_order():
    from dataclasses import fields
    names = [field.name for field in fields(OncoplotOptions)]
    assert names[102:104] == ["gene_bar_color", "gene_bar_x_label"]
    assert names[104] == "legend_gap_points"
    assert names[names.index("sample_id_padding") + 1] == "tmb_bar_color"
