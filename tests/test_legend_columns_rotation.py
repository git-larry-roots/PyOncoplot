"""Issue 11: measured compact metadata columns and histogram label rotation."""
from dataclasses import fields
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
from pandas.testing import assert_frame_equal
import pytest
from pyoncoplot import OncoplotOptions
from test_publication_layout import axes_by_name, small

spec = importlib.util.spec_from_file_location("columns_rotation", Path(__file__).parents[1] / "docs/examples/legend_columns_rotation.py")
example = importlib.util.module_from_spec(spec)
spec.loader.exec_module(example)


@pytest.mark.parametrize("width,angle", [(2400, 0), (2400, 90), (1600, 0), (1600, 90)])
def test_public_geometry_exports_and_invariance(width, angle, tmp_path):
    result = example.publication_example(width=width, angle=angle)
    base = example.publication_example(overrides={"metadata_legend_ncol": None, "tmb_y_label_angle": None})
    fig = result.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    axes = axes_by_name(fig)
    cats = fig.legends[1:]
    assert [[t.get_text() for t in l.get_texts()] for l in cats] == [
        ["Stable", "Responding"], [f"Arm {v}" for v in "ABCDE"], ["Lung", "Liver", "Other"]]
    # Matplotlib's standard column-major order: down each column, then right.
    arm = cats[1]
    positions = [t.get_window_extent(renderer) for t in arm.get_texts()]
    assert sorted(range(5), key=lambda i: (round(positions[i].x0, 5), -positions[i].y0)) == list(range(5))
    assert len(set(round(b.x0, 5) for b in positions)) == 2
    assert positions[0].y0 > positions[1].y0 > positions[2].y0
    assert positions[3].y0 > positions[4].y0
    boxes = [l.get_window_extent(renderer) for l in cats]
    assert sorted(range(3), key=lambda i: (-round((boxes[i].y0+boxes[i].y1)/2, 4), boxes[i].x0)) == [0, 1, 2]
    cb = next(a for a in fig.axes if a not in axes.values())
    assert cb.get_xticks().tolist() == [0, .5, 1]
    assert cb.collections[-1].norm.vmin == 0 and cb.collections[-1].norm.vmax == 1
    assert cb.get_tightbbox(renderer).x1 == pytest.approx(axes["main"].bbox.x1)
    if width == 2400:
        assert all((b.y0+b.y1)/2 == pytest.approx((boxes[0].y0+boxes[0].y1)/2) for b in boxes)
    else:
        assert boxes[0].y0 > cb.get_tightbbox(renderer).y1
    decorations = [l.get_window_extent(renderer) for l in fig.legends] + [cb.get_tightbbox(renderer)]
    sample_boxes = [t.get_window_extent(renderer) for t in axes["metadata"].get_xticklabels()]
    for i, box in enumerate(decorations):
        assert all(not box.overlaps(b) for b in decorations[i+1:] + sample_boxes)
    for box in decorations + [ax.get_tightbbox(renderer) for ax in axes.values()]:
        assert box.x0 >= 0 and box.y0 >= 0
        assert box.x1 <= fig.bbox.x1 and box.y1 <= fig.bbox.y1
    label = axes["tmb"].yaxis.label
    assert label.get_text() == "Recorded\nevents" and label.get_rotation() == angle
    assert label.get_fontsize() == 10
    assert label.get_window_extent(renderer).height <= axes["tmb"].bbox.height
    assert [p.get_width() for p in axes["gene_bar"].patches] == pytest.approx([45,35,25,20,15,5])
    a, b = result.prepared_data, base.prepared_data
    assert a.samples == b.samples and a.genes == b.genes and a.total_samples == b.total_samples
    assert a.row_summary_external and b.row_summary_external
    for name in ("tiles", "metadata", "mutation_counts", "tmb", "row_summary", "cell_states", "main_grid_rows"):
        assert_frame_equal(getattr(a, name), getattr(b, name))
    for suffix, signature in (("png", b"\x89PNG"), ("pdf", b"%PDF"), ("svg", b"<?xml")):
        path = tmp_path / f"columns-{angle}.{suffix}"
        result.save(path, bbox_inches=None)
        assert path.read_bytes().startswith(signature) and path.stat().st_size > 1000
    plt.close("all")


@pytest.mark.parametrize("constraints,columns", [({"metadata_legend_ncol": 1},1),
    ({"metadata_legend_ncol": 2},2), ({"metadata_legend_nrow": 2},3),
    ({"metadata_legend_nrow": 3},2), ({"metadata_legend_ncol": 2, "metadata_legend_nrow": 2},3)])
def test_row_column_constraints(constraints, columns):
    result = example.publication_example(overrides={"metadata_legend_ncol": None, **constraints})
    legend = result.figure.legends[2]
    renderer = result.figure.canvas.get_renderer()
    boxes = [t.get_window_extent(renderer) for t in legend.get_texts()]
    assert len(set(round(b.x0, 5) for b in boxes)) == columns
    # The mutation legend remains automatic, unaffected by metadata controls.
    assert len(set(round(t.get_window_extent(renderer).y0, 5) for t in result.figure.legends[0].get_texts())) == 1
    plt.close("all")


@pytest.mark.parametrize("name", ["metadata_legend_ncol", "metadata_legend_nrow"])
@pytest.mark.parametrize("value", [0,-1,True,False,1.5,2.,"2",float("nan"),float("inf")])
def test_invalid_constraints(name, value):
    with pytest.raises(ValueError, match=name):
        OncoplotOptions(**{name:value})


@pytest.mark.parametrize("value", [True,False,"90",float("nan"),float("inf"),-float("inf")])
def test_invalid_angle(value):
    with pytest.raises(ValueError, match="tmb_y_label_angle"):
        OncoplotOptions(tmb_y_label_angle=value)


def test_explicit_overwide_columns_fail_not_repacked():
    options = {"legend_label_overrides": {v: "Category " + v * 25 for v in "ABCDE"}, "metadata_legend_ncol": 5}
    with pytest.raises(ValueError, match="exceed figure width"):
        example.publication_example(width=1600, overrides=options)
    # Exactly the same labels fit when the caller asks for fewer columns.
    example.publication_example(width=1600, overrides={**options, "metadata_legend_ncol": 2, "height": 1200})
    plt.close("all")


@pytest.mark.parametrize("compact", [True,False])
@pytest.mark.parametrize("angle", [0,90])
def test_direct_renderers_and_defaults(compact, angle):
    from pyoncoplot._matplotlib import render_matplotlib_oncoplot
    from pyoncoplot._plotly import render_plotly_oncoplot
    prepared = small(backend="matplotlib").prepared_data
    fig = render_matplotlib_oncoplot(prepared, palette={}, draw_tmb_bar=True,
        options={"show_tmb_y_label": True, "tmb_y_label": "First\nsecond", "tmb_y_label_angle": angle,
                 "compact_layout":compact, "metadata_legend_position":"bottom"})
    label = axes_by_name(fig)["tmb"].yaxis.label
    assert label.get_rotation() == angle and label.get_text() == "First\nsecond"
    with pytest.raises(NotImplementedError, match="tmb_y_label_angle"):
        render_plotly_oncoplot(prepared, palette={}, options={"tmb_y_label_angle":angle})
    assert OncoplotOptions().tmb_y_label_angle is None
    names = tuple(f.name for f in fields(OncoplotOptions))
    assert names[105:] == ("gene_bar_label_position", "tmb_y_label_angle")
    plt.close("all")
