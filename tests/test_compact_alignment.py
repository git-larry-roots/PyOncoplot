"""Issue 7: measured matrix-aligned shelves and independent label styling."""
from dataclasses import fields
import importlib.util
from pathlib import Path
import subprocess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest
from pandas.testing import assert_frame_equal
from pyoncoplot import OncoplotOptions
from test_publication_layout import axes_by_name, small

spec = importlib.util.spec_from_file_location("compact", Path(__file__).parents[1] / "docs/examples/compact_alignment.py")
example = importlib.util.module_from_spec(spec)
spec.loader.exec_module(example)


@pytest.mark.parametrize("width,shared", [(2200, True), (1600, False)])
def test_matrix_band_and_exports(width, shared, tmp_path):
    result = example.publication_example(width=width, overrides={"height": 560})
    fig = result.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    axes = axes_by_name(fig)
    main = axes["main"].bbox
    cb = next(a for a in fig.axes if a not in axes.values())
    legends = [l.get_window_extent(renderer) for l in fig.legends]
    colorbox = cb.get_tightbbox(renderer)
    assert any(b.overlaps(colorbox) for b in legends) is False
    assert any(b.y0 < colorbox.y1 and colorbox.y0 < b.y1 for b in legends) == shared
    boxes = legends + [colorbox]
    assert min(b.x0 for b in boxes) == pytest.approx(main.x0)
    assert all(main.x0 - 1e-6 <= b.x0 and b.x1 <= main.x1 + 1e-6 for b in boxes)
    assert cb.get_xlabel() == "Tumor\npurity"
    assert all(l.get_title().get_text() == "" for l in fig.legends)
    texts = [t.get_text() for l in fig.legends for t in l.get_texts()]
    assert all(f"{arm} (n=3)" in texts for arm in "ABCD")
    assert "Missense (n=5)" in texts
    assert axes["tmb"].get_ylabel() == "Recorded\nevents"
    assert axes["gene_bar"].get_xlabel() == "Independent cohort\nAffected (%)"
    for ax in axes.values():
        assert not ax.xaxis.label.get_wrap() and not ax.yaxis.label.get_wrap() and not ax.title.get_wrap()
        boxes.append(ax.get_tightbbox(renderer))
    assert all(b.x0 >= 0 and b.y0 >= 0 and b.x1 <= fig.bbox.x1 and b.y1 <= fig.bbox.y1 for b in boxes)
    labels = [t.get_window_extent(renderer) for t in axes["metadata"].get_yticklabels()]
    assert all(not a.overlaps(b) for i,a in enumerate(labels) for b in labels[i+1:])
    samples = [t.get_window_extent(renderer) for t in axes["metadata"].get_xticklabels()]
    assert max(b.y1 for b in legends + [colorbox]) < min(b.y0 for b in samples)
    for name in ("tmb", "metadata"):
        np.testing.assert_allclose(axes[name].bbox.intervalx, main.intervalx)
    np.testing.assert_allclose(axes["gene_bar"].bbox.intervaly, main.intervaly)
    for ext, signature in (("png", b"\x89PNG"), ("pdf", b"%PDF"), ("svg", b"<?xml")):
        path = tmp_path / f"compact.{ext}"
        result.save(path, bbox_inches=None)
        assert path.read_bytes().startswith(signature)
        assert path.stat().st_size > 1000
    plt.close(fig)


def test_short_figure_and_scientific_invariance():
    plain = example.publication_example(overrides={"height": 800, "legend_label_overrides": {}, "legend_gap_points": 6})
    styled = example.publication_example()
    assert styled.figure.get_figheight() == pytest.approx(520/120)
    assert axes_by_name(styled.figure)["main"].bbox.height < axes_by_name(plain.figure)["main"].bbox.height
    a,b = plain.prepared_data, styled.prepared_data
    assert a.samples == b.samples and a.genes == b.genes
    assert a.row_summary_external and b.row_summary_external
    for key in ("row_summary", "tiles", "mutation_counts", "tmb", "metadata", "cell_states", "main_grid_rows"):
        assert_frame_equal(getattr(a,key), getattr(b,key))
    assert b.row_summary.Count.tolist() == [9,7,5,4,3,1]
    assert b.row_summary.Denominator.tolist() == [20]*6
    assert b.metadata.cohort_arm.tolist() == list("ABCD")*3
    plt.close("all")


@pytest.mark.parametrize("show,global_show", [(True,False),(False,True)])
def test_numeric_title_override_both_orientations(show, global_show):
    result = small(backend="matplotlib", options={"show_legend_titles": global_show,
        "metadata_numeric_scales": {"value": {"show_title": show}}})
    cb = next(a for a in result.figure.axes if a not in axes_by_name(result.figure).values())
    assert bool(cb.get_title()) == show
    plt.close(result.figure)
    result = example.publication_example(overrides={"height": 700, "show_legend_titles": global_show,
        "metadata_numeric_scales": {"purity": {"show_title": show}}})
    cb = next(a for a in result.figure.axes if a not in axes_by_name(result.figure).values())
    assert bool(cb.get_xlabel()) == show
    plt.close(result.figure)


@pytest.mark.parametrize("options", [{"legend_gap_points": -1}, {"legend_gap_points": np.nan}, {"legend_gap_points": True},
    {"metadata_numeric_scales": {"value": {"show_title": "yes"}}}])
def test_invalid_controls(options):
    with pytest.raises(ValueError):
        OncoplotOptions(**options)


@pytest.mark.parametrize("options", [{"legend_gap_points": 2}, {"metadata_numeric_scales": {"value": {"show_title": True}}}])
def test_plotly_explicit_errors(options):
    with pytest.raises(NotImplementedError, match="matplotlib"):
        small(backend="plotly", options=options)


def test_complete_positional_prefix():
    source = subprocess.check_output(["git", "show", "f3542c9:src/pyoncoplot/_options.py"], text=True)
    import ast
    baseline = [n.target.id for n in ast.parse(source).body[-1].body if isinstance(n, ast.AnnAssign)]
    assert [f.name for f in fields(OncoplotOptions)][:len(baseline)] == baseline


def test_oversized_labels_fail_instead_of_clipping():
    with pytest.raises(ValueError, match="exceed figure width"):
        example.publication_example(overrides={"legend_label_overrides": {"A": "Long "*300}})
    plt.close("all")


def test_long_newlines_and_every_visible_text():
    result = example.publication_example(overrides={
        "height": 700,
        "title_text": "Explicit plot title\nkeeps exactly two lines",
        "metadata_numeric_scales": {"purity": {"min": -10000, "max": 10000,
            "ticks": [-10000, 0, 10000], "show_title": True,
            "colorbar": {"length": 130, "label_pad": 8}}},
        "legend_title_overrides": {"cohort_arm": "Study\narm", "purity": "Long numeric\nheading preserved"},
    })
    fig = result.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes = []
    for ax in fig.axes:
        xticks = [t for v,t in zip(ax.get_xticks(), ax.get_xticklabels()) if min(ax.get_xlim()) <= v <= max(ax.get_xlim())]
        yticks = [t for v,t in zip(ax.get_yticks(), ax.get_yticklabels()) if min(ax.get_ylim()) <= v <= max(ax.get_ylim())]
        for t in [ax.xaxis.label, ax.yaxis.label, ax.title] + xticks + yticks + list(ax.texts):
            if t.get_visible() and t.get_text():
                boxes.append(t.get_window_extent(renderer))
    assert all(not a.overlaps(b) for i,a in enumerate(boxes) for b in boxes[i+1:])
    assert all(b.x0 >= 0 and b.y0 >= 0 and b.x1 <= fig.bbox.x1 and b.y1 <= fig.bbox.y1 for b in boxes)
    assert fig._suptitle.get_text() == "Explicit plot title\nkeeps exactly two lines"
    assert not fig._suptitle.get_wrap()
    plt.close(fig)


def test_direct_renderers_and_legacy_size_defaults():
    from pyoncoplot._matplotlib import render_matplotlib_oncoplot
    from pyoncoplot._plotly import render_plotly_oncoplot
    baseline = small(backend="matplotlib", options={"height": 480})
    assert baseline.figure.get_figheight() == 5
    prepared = baseline.prepared_data
    with pytest.raises(NotImplementedError, match="matplotlib"):
        render_plotly_oncoplot(prepared, palette={}, options={"legend_gap_points": 1})
    result = render_matplotlib_oncoplot(prepared, palette={}, options={
        "compact_layout": True, "metadata_legend_position": "bottom", "height": 480,
        "metadata_numeric_scales": {"value": {"show_title": True}}, "show_legend_titles": False})
    assert result.get_figheight() == 4
    cb = next(ax for ax in result.axes if ax.get_label() not in {"main", "metadata"})
    assert cb.get_xlabel()
    plt.close("all")


def test_too_short_matrix_rejected():
    with pytest.raises(ValueError, match="increase width/height"):
        example.publication_example(overrides={"height": 300})
    plt.close("all")
