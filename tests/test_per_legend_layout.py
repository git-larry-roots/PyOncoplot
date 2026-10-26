"""Per-source constraints: geometry and scientific invariance on synthetic inputs."""
from collections import ChainMap
from dataclasses import fields
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from pandas.testing import assert_frame_equal
import pytest
from pyoncoplot import OncoplotOptions, load_oncoplot_params
from pyoncoplot._matplotlib import render_matplotlib_oncoplot
from pyoncoplot._plotly import render_plotly_oncoplot
from test_publication_layout import axes_by_name

spec = importlib.util.spec_from_file_location("per_legend", Path(__file__).parents[1] / "docs/examples/per_legend_layout.py")
example = importlib.util.module_from_spec(spec)
spec.loader.exec_module(example)


def rendered_columns(legend, renderer):
    boxes = [t.get_window_extent(renderer) for t in legend.get_texts()]
    assert sorted(range(len(boxes)), key=lambda i: (round(boxes[i].x0, 5), -boxes[i].y0)) == list(range(len(boxes)))
    return len({round(b.x0, 5) for b in boxes})


@pytest.mark.parametrize("constraints,expected", [
    ({}, [1, 2, 2]),
    ({"metadata_legend_ncol": {}, "metadata_legend_nrow": {}}, [2, 5, 3]),
    ({"metadata_legend_ncol": {"cohort_arm": 2}, "metadata_legend_nrow": {}}, [2, 2, 3]),
    ({"metadata_legend_ncol": 2, "metadata_legend_nrow": {"cohort_arm": 2}}, [2, 3, 2]),
    ({"metadata_legend_ncol": {"cohort_arm": 1}, "metadata_legend_nrow": 2}, [1, 3, 2]),
    ({"metadata_legend_ncol": ChainMap({"z_response": 1, "cohort_arm": 2}), "metadata_legend_nrow": {"a_site": 2}}, [1, 2, 2]),
])
def test_independent_constraints(constraints, expected):
    r = example.publication_example(overrides=constraints)
    assert [rendered_columns(l, r.figure.canvas.get_renderer()) for l in r.figure.legends[1:]] == expected
    plt.close("all")


@pytest.mark.parametrize("width,angle", [(2400, 0), (2400, 90), (1600, 90)])
def test_short_thin_geometry_exports_invariance(width, angle, tmp_path):
    r = example.publication_example(width=width, angle=angle)
    base = example.publication_example(overrides={"metadata_legend_ncol": None, "metadata_legend_nrow": None,
                                               "tmb_height_ratio": .14, "metadata_height_ratio": .35})
    fig = r.figure
    renderer = fig.canvas.get_renderer()
    axes = axes_by_name(fig)
    tmb, meta, main = (axes[k] for k in ("tmb", "metadata", "main"))
    assert tmb.bbox.height < axes_by_name(base.figure)["tmb"].bbox.height
    assert meta.bbox.height < axes_by_name(base.figure)["metadata"].bbox.height
    assert tmb.yaxis.label.get_fontsize() == 10 and meta.get_yticklabels()[0].get_fontsize() == 10
    assert tmb.yaxis.label.get_window_extent(renderer).height <= tmb.bbox.height
    # The compact fitter stays close to physical text minima, not default panel sizes.
    assert tmb.bbox.height < tmb.yaxis.label.get_window_extent(renderer).height + 8
    label_boxes = [t.get_window_extent(renderer) for t in meta.get_yticklabels()]
    assert all(2 <= a.y0 - b.y1 < 4 for a, b in zip(label_boxes, label_boxes[1:]))
    for ax in (tmb, meta):
        np.testing.assert_allclose(ax.transData.transform([(0, 0), (12, 0)])[:, 0],
                                   main.transData.transform([(0, 0), (12, 0)])[:, 0])
    cb = next(ax for ax in fig.axes if ax not in axes.values())
    assert cb.get_tightbbox(renderer).x1 == pytest.approx(main.bbox.x1)
    cats = fig.legends[1:]
    assert [[t.get_text() for t in l.get_texts()] for l in cats] == [
        ["Stable", "Responding"], [f"Arm {v}" for v in "ABCDE"], ["Lung", "Liver", "Other"]]
    boxes = [l.get_window_extent(renderer) for l in fig.legends] + [cb.get_tightbbox(renderer)]
    samples = [t.get_window_extent(renderer) for t in meta.get_xticklabels()]
    for i, b in enumerate(boxes):
        assert all(not b.overlaps(other) for other in boxes[i+1:] + samples)
    for b in boxes + [ax.get_tightbbox(renderer) for ax in axes.values()]:
        assert 0 <= b.x0 < b.x1 <= fig.bbox.x1 and 0 <= b.y0 < b.y1 <= fig.bbox.y1
    if width == 2400:
        centers = [(l.get_window_extent(renderer).y0+l.get_window_extent(renderer).y1)/2 for l in cats]
        assert centers == pytest.approx([centers[0]] * 3)
    a, b = r.prepared_data, base.prepared_data
    assert (a.samples, a.genes, a.total_samples, a.row_summary_external) == (b.samples, b.genes, b.total_samples, b.row_summary_external)
    baseline_metadata = axes_by_name(base.figure)["metadata"]
    np.testing.assert_allclose([p.get_facecolor() for p in meta.patches],
                               [p.get_facecolor() for p in baseline_metadata.patches])
    for name in ("tiles", "metadata", "mutation_counts", "tmb", "row_summary", "cell_states", "main_grid_rows"):
        assert_frame_equal(getattr(a, name), getattr(b, name))
    for suffix, signature in (("png", b"\x89PNG"), ("pdf", b"%PDF"), ("svg", b"<?xml")):
        path = tmp_path / f"per-legend-{width}-{angle}.{suffix}"
        r.save(path, bbox_inches=None)
        assert path.read_bytes().startswith(signature)
    plt.close("all")


@pytest.mark.parametrize("name", ["metadata_legend_ncol", "metadata_legend_nrow"])
@pytest.mark.parametrize("value", [0, -1, True, False, 2., 1.5, "2", None, float("nan"), float("inf")])
def test_invalid_mapping_counts(name, value):
    with pytest.raises(ValueError, match=name):
        OncoplotOptions(**{name: {"cohort_arm": value}})


@pytest.mark.parametrize("name", ["metadata_legend_ncol", "metadata_legend_nrow"])
@pytest.mark.parametrize("key", [1, True, None])
def test_non_string_keys(name, key):
    with pytest.raises(ValueError, match="keys.*strings"):
        OncoplotOptions(**{name: {key: 2}})


@pytest.mark.parametrize("renderer", [render_matplotlib_oncoplot, render_plotly_oncoplot])
@pytest.mark.parametrize("name", ["metadata_legend_ncol", "metadata_legend_nrow"])
@pytest.mark.parametrize("key", ["Study\narm", "metadata:cohort_arm", "missing"])
def test_unknown_source_names_even_hidden(renderer, name, key):
    prepared = example.publication_example().prepared_data
    with pytest.raises(ValueError, match=name + ".*unknown metadata"):
        renderer(prepared, palette={}, options={name: {key: 2}, "show_metadata_legends": False})
    plt.close("all")


@pytest.mark.parametrize("renderer", [render_matplotlib_oncoplot, render_plotly_oncoplot])
@pytest.mark.parametrize("value", [{}, {"cohort_arm": 2}])
def test_unsupported_backend_layout_is_explicit(renderer, value):
    prepared = example.publication_example().prepared_data
    with pytest.raises(NotImplementedError, match="metadata_legend_ncol/nrow mappings require"):
        renderer(prepared, palette={}, options={"metadata_legend_ncol": value})
    plt.close("all")


def test_yaml_source_constraints_and_numeric_policy(tmp_path):
    path = tmp_path / "layout.yaml"
    path.write_text("""options:
  metadata_legend_ncol: {z_response: 1, cohort_arm: 2}
  metadata_legend_nrow: {a_site: 2}
""")
    options = load_oncoplot_params(path)["options"]
    result = example.publication_example(overrides=options)
    assert [rendered_columns(l, result.figure.canvas.get_renderer()) for l in result.figure.legends[1:]] == [1, 2, 2]
    # Numeric source names are valid but do not turn colorbars into categorical legends.
    numeric = example.publication_example(overrides={"metadata_legend_ncol": {"purity": 2}, "metadata_legend_nrow": {}})
    assert len(numeric.figure.legends) == 4
    assert len(numeric.figure.axes) == 5
    plt.close("all")


def test_positional_prefix_is_fixed_not_git_dependent():
    names = tuple(f.name for f in fields(OncoplotOptions))
    assert names[:8] == ("width", "height", "selection_type", "tmb_height_ratio", "gene_bar_width_ratio",
                         "metadata_height_ratio", "metadata_position", "buffer_metadata")
    assert names[25:27] == ("metadata_legend_nrow", "metadata_legend_ncol")
    assert names[105:107] == ("gene_bar_label_position", "tmb_y_label_angle")


@pytest.mark.parametrize("payload,match", [
    ('{cohort_arm: true}', 'positive integer'),
    ('{cohort_arm: .inf}', 'positive integer'),
    ('{Study: 2}', 'unknown metadata'),
    ('{1: 2}', 'keys.*strings'),
])
def test_invalid_yaml_options(tmp_path, payload, match):
    path = tmp_path / "invalid.yaml"
    path.write_text("options:" + chr(10) + "  metadata_legend_ncol: " + payload + chr(10))
    with pytest.raises(ValueError, match=match):
        options = load_oncoplot_params(path)["options"]
        example.publication_example(overrides=options)
    plt.close("all")


def test_independent_display_layout_preserves_external_summaries_and_states():
    base = example.publication_example()
    changed = example.publication_example(overrides={
        "metadata_legend_order": ["cohort_arm", "purity", "z_response", "a_site"],
        "metadata_legend_alignment": {"cohort_arm": "top", "purity": "top", "z_response": "top", "a_site": "top"},
        "legend_gap_points": 1})
    a, b = changed.prepared_data, base.prepared_data
    assert (a.samples, a.genes, a.total_samples, a.row_summary_external, a.metadata_cols) == (b.samples, b.genes, b.total_samples, b.row_summary_external, b.metadata_cols)
    assert a.row_summary_external
    for name in ("tiles", "metadata", "mutation_counts", "tmb", "row_summary", "cell_states", "main_grid_rows"):
        assert_frame_equal(getattr(a, name), getattr(b, name))
    for name in ("main", "metadata", "gene_bar", "tmb"):
        old = axes_by_name(base.figure)[name]
        new = axes_by_name(changed.figure)[name]
        np.testing.assert_allclose([p.get_facecolor() for p in old.patches], [p.get_facecolor() for p in new.patches])
        np.testing.assert_allclose([p.get_width() for p in old.patches], [p.get_width() for p in new.patches])
        np.testing.assert_allclose([p.get_height() for p in old.patches], [p.get_height() for p in new.patches])
    plt.close("all")
