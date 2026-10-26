"""Issue 9: source-ordered shelves, compact proportions and shared labels."""
from dataclasses import fields
import importlib.util
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pytest
from pandas.testing import assert_frame_equal
from pyoncoplot import OncoplotOptions
from test_publication_layout import axes_by_name, small

spec = importlib.util.spec_from_file_location("final_compact", Path(__file__).parents[1] / "docs/examples/final_compact_layout.py")
example = importlib.util.module_from_spec(spec)
spec.loader.exec_module(example)


@pytest.mark.parametrize("width,shared", [(2400, True), (1600, False)])
def test_source_order_geometry_exports(width, shared, tmp_path):
    result = example.publication_example(width=width)
    fig = result.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    axes = axes_by_name(fig)
    cats = fig.legends[1:]
    expected = [["Stable", "Responding"], [f"{v} (n=3)" for v in "ABCD"], ["Lung", "Liver", "Other"]]
    assert [[t.get_text() for t in l.get_texts()] for l in cats] == expected
    # Check *rendered* reading order, including when the last group wraps.
    boxes = [l.get_window_extent(renderer) for l in cats]
    order = sorted(range(len(boxes)), key=lambda i: (-round((boxes[i].y0+boxes[i].y1)/2, 4), boxes[i].x0))
    assert order == [0, 1, 2]
    cb = next(a for a in fig.axes if a not in axes.values())
    cbbox = cb.get_tightbbox(renderer)
    assert cb.get_xlabel() == "Tumor\npurity"
    assert cb.get_xticks().tolist() == [0, .5, 1]
    assert cbbox.x1 == pytest.approx(axes["main"].bbox.x1)
    if shared:
        for box in boxes:
            assert (box.y0+box.y1)/2 == pytest.approx((cbbox.y0+cbbox.y1)/2)
        assert boxes[-1].x1 < cbbox.x0
    else:
        assert any(box.y0 > cbbox.y1 for box in boxes)
    assert [t.get_text() for t in axes["metadata"].get_yticklabels()] == ["Treatment\nresponse", "Tumor\npurity", "Study\narm", "Sample\nsite"]
    assert result.prepared_data.metadata_cols == ["z_response", "purity", "cohort_arm", "a_site"]
    assert axes["tmb"].bbox.height < axes["main"].bbox.height * .25
    if shared:
        assert axes["metadata"].bbox.height / 4 < axes["main"].bbox.height / 6
    # Narrow wrapping consumes height: retain measured multiline clearance.
    assert axes["metadata"].bbox.height < axes["main"].bbox.height
    assert axes["tmb"].yaxis.label.get_rotation() == 0
    assert axes["tmb"].yaxis.label.get_fontsize() == 10
    for name in ("metadata", "tmb"):
        np.testing.assert_allclose(axes[name].bbox.intervalx, axes["main"].bbox.intervalx)
    np.testing.assert_allclose(axes["gene_bar"].bbox.intervaly, axes["main"].bbox.intervaly)
    for suffix, signature in (("png", b"\x89PNG"), ("pdf", b"%PDF"), ("svg", b"<?xml")):
        path = tmp_path / f"final.{suffix}"
        result.save(path, bbox_inches=None)
        assert path.read_bytes().startswith(signature) and path.stat().st_size > 1000
    plt.close(fig)


@pytest.mark.parametrize("mode", ["counts", "percent"])
def test_shared_labels_external_semantics_and_invariance(mode):
    base = example.publication_example(overrides={"gene_bar_mode": mode, "gene_bar_label_position": "bar_end"})
    column = example.publication_example(overrides={"gene_bar_mode": mode})
    a, b = base.prepared_data, column.prepared_data
    assert a.samples == b.samples and a.genes == b.genes
    assert a.row_summary_external and b.row_summary_external
    for name in ("tiles", "metadata", "mutation_counts", "tmb", "row_summary", "cell_states", "main_grid_rows"):
        assert_frame_equal(getattr(a, name), getattr(b, name))
    bars = axes_by_name(column.figure)["gene_bar"]
    xs = [t.get_window_extent(column.figure.canvas.get_renderer()).x0 for t in bars.texts]
    np.testing.assert_allclose(xs, xs[0])
    assert [t.get_text() for t in bars.texts] == ["9 (45%)", "7 (35%)", "5 (25%)", "4 (20%)", "3 (15%)", "1 (5%)"]
    assert len(set(t.get_position()[0] for t in axes_by_name(base.figure)["gene_bar"].texts)) > 1
    assert [p.get_width() for p in bars.patches] == pytest.approx([45,35,25,20,15,5] if mode == "percent" else [9,7,5,4,3,1])
    plt.close("all")


def test_computed_column_and_direct_backend_scope():
    from pyoncoplot._matplotlib import render_matplotlib_oncoplot
    from pyoncoplot._plotly import render_plotly_oncoplot
    base = small(backend="matplotlib")
    prepared = base.prepared_data
    fig = render_matplotlib_oncoplot(prepared, palette={}, draw_gene_bar=True,
        options={"show_gene_bar_labels": True, "gene_bar_label_position": "column", "gene_bar_label_format": "count_percent"})
    bars = axes_by_name(fig)["gene_bar"]
    assert len(bars.texts) == len(prepared.genes) > 0
    totals = prepared.tiles.groupby("Gene", observed=False).size().reindex(prepared.genes, fill_value=0)
    assert [t.get_text() for t in bars.texts] == [f"{n} ({round(n/prepared.total_samples*100):g}%)" for n in totals]
    assert len(set(t.get_position()[0] for t in bars.texts)) == 1
    with pytest.raises(NotImplementedError, match="matplotlib"):
        render_plotly_oncoplot(prepared, palette={}, options={"gene_bar_label_position": "column"})
    assert OncoplotOptions().gene_bar_label_position == "bar_end"
    plt.close("all")


@pytest.mark.parametrize("value", ["fixed", "", None, 1])
def test_invalid_label_position(value):
    with pytest.raises(ValueError, match="gene_bar_label_position"):
        OncoplotOptions(gene_bar_label_position=value)


def test_complete_issue7_positional_prefix():
    # Reuse the *literal*, fully enumerated previous snapshot, then include the
    # field appended in issue 7. Locate by name, never by AST statement index.
    import ast
    source = ast.parse((Path(__file__).parent / "test_compact_alignment.py").read_text())
    function = next(n for n in source.body if isinstance(n, ast.FunctionDef) and n.name == "test_complete_positional_prefix")
    assignment = next(n for n in function.body if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "baseline" for t in n.targets))
    baseline = ast.literal_eval(assignment.value) + ("legend_gap_points",)
    assert len(baseline) == 105
    names = tuple(f.name for f in fields(OncoplotOptions))
    assert names[:105] == baseline
    assert names[105:] == ("gene_bar_label_position",)


def test_multiple_numeric_source_order_and_titles():
    import pandas as pd
    from pyoncoplot import oncoplot
    events = pd.DataFrame({"gene": ["G", "H"], "sample": ["S1", "S2"]})
    metadata = pd.DataFrame({"sample": ["S1", "S2"], "z_cat": ["Z1", "Z2"],
        "z_number": [10.,20.], "a_cat": ["A1", "A2"], "a_number": [100.,200.]})
    result = oncoplot(events, gene_col="gene", sample_col="sample", metadata=metadata,
        metadata_sample_col="sample", metadata_cols=["z_cat", "z_number", "a_cat", "a_number"],
        backend="matplotlib", options={"compact_layout": True, "width": 2400,
        "metadata_legend_position": "bottom", "font_family": "DejaVu Sans",
        "show_legend_titles": True, "legend_title_overrides": {"z_number": "First scale", "a_number": "Second scale"},
        "metadata_numeric_scales": {"z_number": {"ticks": [10,20]}, "a_number": {"ticks": [100,200]}}})
    fig = result.figure
    axes = axes_by_name(fig)
    colorbars = [a for a in fig.axes if a not in axes.values()]
    assert [a.get_xlabel() for a in colorbars] == ["First scale", "Second scale"]
    assert [a.get_xticks().tolist() for a in colorbars] == [[10,20], [100,200]]
    assert colorbars[0].bbox.x1 < colorbars[1].bbox.x0
    assert [l.get_title().get_text() for l in fig.legends][-2:] == ["Z Cat", "A Cat"]
    assert [t.get_text() for t in fig.legends[-2].get_texts()] == ["Z1", "Z2"]
    plt.close(fig)
