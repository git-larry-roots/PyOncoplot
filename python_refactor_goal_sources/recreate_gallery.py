"""Recreate the direct public ``oncoplot()`` gallery."""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path
from typing import Optional, Sequence

from PIL import Image, ImageDraw
from pyoncoplot import oncoplot

try:
    from python_refactor_goal_sources.make_oncoplots import (
        CLEAN_OUT,
        COMPARISON_OUT,
        CONFIG_PATH,
        GALLERY_CONFIG,
        GENERATED_ROOT,
        GOAL_PLOTS,
        INPUTS,
    )
except ModuleNotFoundError:  # pragma: no cover - supports ``python path/to/script.py``.
    from make_oncoplots import (  # type: ignore
        CLEAN_OUT,
        COMPARISON_OUT,
        CONFIG_PATH,
        GALLERY_CONFIG,
        GENERATED_ROOT,
        GOAL_PLOTS,
        INPUTS,
    )


def _render_oncoplot_preset(name: str, out_dir: Path) -> Path:
    run_config = GALLERY_CONFIG["plot_runs"].get(name)
    if run_config is None:
        raise ValueError(f"Unknown gallery preset {name!r}.")
    if run_config.get("style") != "clean" or run_config["renderer"] != "oncoplot":
        raise ValueError(f"Preset {name!r} is not a direct oncoplot gallery preset.")

    out_dir.mkdir(parents=True, exist_ok=True)
    output_path = out_dir / run_config["output_name"]
    save = dict(run_config["params"]["save"])
    save["path"] = output_path
    result = oncoplot(
        params=CONFIG_PATH,
        params_key=f"gallery_params.plot_runs.{name}.params.oncoplot",
        save=save,
    )
    if result.backend == "matplotlib":
        import matplotlib.pyplot as plt

        plt.close(result.figure)
    return output_path


def _render_comparison_sheet(name: str, out_dir: Path) -> Path:
    comparison_config = GALLERY_CONFIG.get("comparison_runs", {})
    if name not in comparison_config:
        raise ValueError("Comparison sheets are only available for configured BRCA gallery presets.")

    run_config = comparison_config[name]
    clean_name = run_config["clean_preset"]
    clean_run = GALLERY_CONFIG["plot_runs"][clean_name]
    out_dir.mkdir(parents=True, exist_ok=True)
    output_path = out_dir / run_config["output_name"]

    with tempfile.TemporaryDirectory(prefix="pyoncoplot_compare_") as temporary_dir:
        work_dir = Path(temporary_dir)
        clean_path = _render_oncoplot_preset(clean_name, work_dir / "clean")
        original_path = GOAL_PLOTS / clean_run["goal_plot"]
        labels = ["Original reference", "Clean generated"]
        paths = [original_path, clean_path]
        thumb_width = 620
        thumb_height = 360
        label_height = 38
        canvas = Image.new("RGB", (thumb_width * len(paths), thumb_height + label_height), "white")
        draw = ImageDraw.Draw(canvas)
        for index, (label, path) in enumerate(zip(labels, paths)):
            if path.exists():
                with Image.open(path) as source_image:
                    image = source_image.convert("RGB")
                image.thumbnail((thumb_width - 20, thumb_height - 20), Image.Resampling.LANCZOS)
                x = index * thumb_width + (thumb_width - image.width) // 2
                y = label_height + (thumb_height - image.height) // 2
                canvas.paste(image, (x, y))
            draw.text((index * thumb_width + 12, 10), label, fill="black")
        canvas.save(output_path)
    return output_path


def render_preset(name: str, out_dir: Optional[Path] = None, style: str = "clean") -> Path:
    if style == "comparison":
        return _render_comparison_sheet(name, out_dir or COMPARISON_OUT)
    if style != "clean":
        raise ValueError("style must be one of: clean, comparison.")
    return _render_oncoplot_preset(name, out_dir or CLEAN_OUT)


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--style", choices=["clean", "comparison", "both"], default="clean")
    parser.add_argument("--preset", action="append", help="Preset name to render. May be supplied more than once.")
    parser.add_argument("--out-dir", type=Path, help="Override the output directory for the selected style.")
    args = parser.parse_args(argv)

    if args.style in {"clean", "both"}:
        clean_out = args.out_dir if args.style == "clean" and args.out_dir else CLEAN_OUT
        runs = GALLERY_CONFIG["plot_runs"]
        selected_names = args.preset or [
            name for name, run in runs.items() if run.get("style") == "clean" and run.get("run", True)
        ]
        for name in selected_names:
            render_preset(name, clean_out, style="clean")
        print(f"Wrote clean gallery images to {clean_out}")

    if args.style in {"comparison", "both"}:
        comparison_out = args.out_dir if args.style == "comparison" and args.out_dir else COMPARISON_OUT
        comparison_names = args.preset or ["brca_large"]
        for name in comparison_names:
            render_preset(name, comparison_out, style="comparison")
        print(f"Wrote comparison sheets to {comparison_out}")


if __name__ == "__main__":
    main()
