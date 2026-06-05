# Python Refactor Goal Sources

This directory is the canonical source area for reproducing the retained
training and gallery plots used by `pyoncoplot`.

## Folder Layout

```text
python_refactor_goal_sources/
  config.yaml                      config-driven named gallery runs
  make_oncoplots.py                direct public oncoplot() gallery renderer
  recreate_gallery.py              gallery orchestration CLI
  generate_synthetic_inputs.py     deterministic synthetic input generator for non-fuc fixtures
  fuc_sources/                     copied fuc example scripts, rebuild helper, and manifest
  source_info_for_training.md      this source map
  goal_plots/                      numbered source/reference plots
  syntheitic_goal_data/            deterministic TSV/JSON inputs
  generated_plots/                 tracked generated outputs for gallery review
  ggoncoplot/                      cloned R reference implementation
```

The folder names `python_refactor_goal_sources` and `syntheitic_goal_data` keep
their current spelling intentionally.

## Source Rules

- Active goal plots are direct public `oncoplot()` targets.
- `goal_plots/` stores numbered source/reference material.
- Generated images must never be written to `goal_plots/`.
- `generated_plots/` stores reviewed generated outputs and is tracked when plots
  are intentionally updated.
- All active source reference plots are PNG files named `goal_plot_NN.png`.
- Clean generated counterparts are configured in `config.yaml` as
  `gen.goal_plot_01.png` through `gen.goal_plot_13.png`.

## Source Sections

Goal plots are compactly numbered by retained source family:

| Range | Section | Description |
| --- | --- | --- |
| `goal_plot_01.png` - `goal_plot_06.png` | ggoncoplot/R | ggoncoplot README and paper oncoplot examples |
| `goal_plot_07.png` - `goal_plot_09.png` | other R-based paper | Published oncoplot-style figures from non-ggoncoplot papers |
| `goal_plot_10.png` - `goal_plot_13.png` | Python/fuc AML | fuc tutorial AML oncoplot examples |

## Generated Counterparts

| Goal plot | Clean generated plot | Source section | Renderer intent |
| --- | --- | --- | --- |
| `goal_plot_01.png` | `gen.goal_plot_01.png` | ggoncoplot/R | large BRCA oncoplot |
| `goal_plot_02.png` | `gen.goal_plot_02.png` | ggoncoplot/R | small README oncoplot |
| `goal_plot_03.png` | `gen.goal_plot_03.png` | ggoncoplot/R | README basic oncoplot |
| `goal_plot_04.png` | `gen.goal_plot_04.png` | ggoncoplot/R | README oncoplot with marginal bars |
| `goal_plot_05.png` | `gen.goal_plot_05.png` | ggoncoplot/R | README oncoplot with clinical metadata |
| `goal_plot_06.png` | `gen.goal_plot_06.png` | ggoncoplot/R | compact paper-style GBM oncoplot |
| `goal_plot_07.png` | `gen.goal_plot_07.png` | other R-based paper | compact BRCA oncoplot |
| `goal_plot_08.png` | `gen.goal_plot_08.png` | other R-based paper | compact CSSC alteration matrix |
| `goal_plot_09.png` | `gen.goal_plot_09.png` | other R-based paper | GBM clinical/molecular panel |
| `goal_plot_10.png` | `gen.goal_plot_10.png` | Python/fuc AML | basic AML oncoplot |
| `goal_plot_11.png` | `gen.goal_plot_11.png` | Python/fuc AML | AML metadata oncoplot, unsorted |
| `goal_plot_12.png` | `gen.goal_plot_12.png` | Python/fuc AML | AML metadata oncoplot, sorted |
| `goal_plot_13.png` | `gen.goal_plot_13.png` | Python/fuc AML | AML metadata oncoplot, survival-filtered |

## Goal Plot Source Table

| Goal plot | Source status | Original local source path | Original filename | Source project/page | Source URL | License/attribution | Dimensions | SHA256 | Import action | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `goal_plot_01.png` | resolved | `/ggoncoplot/paper/oncoplot.png` | `73b90cf5-8b0b-4787-b38e-f0b1a2d6.png`; matches `oncoplot.png` | ggoncoplot paper Figure 1 | https://github.com/selkamand/ggoncoplot/tree/main/paper | ggoncoplot; MIT + file LICENSE; Sam El-Kamand, Julian M. W. Quinn, Mark J. Cowley | 3600 x 1800 | `a062b5f1868a18eb5754b149c91f0609b1c127d4207de4f0e7dfab48dadd85c1` | preserved existing PNG | Exact SHA match to nested `paper/oncoplot.png`; duplicate was not re-imported. |
| `goal_plot_02.png` | resolved | `/ggoncoplot/man/figures/README-example-1.png` | `README-example-1.png` | ggoncoplot GitHub man/figures | https://github.com/selkamand/ggoncoplot/tree/main/man/figures | ggoncoplot; MIT + file LICENSE; Sam El-Kamand | 672 x 480 | `4f7bd78b8c0b7b7c92f1fd6397e3d0bb2b92098a0314721bcd9ebd6495c33451` | copied PNG | README example figure. |
| `goal_plot_03.png` | resolved | `/ggoncoplot/man/figures/README-unnamed-chunk-2-1.png` | `README-unnamed-chunk-2-1.png` | ggoncoplot README/pkgdown basic example | https://selkamand.github.io/ggoncoplot/ | ggoncoplot; MIT + file LICENSE; Sam El-Kamand | 7200 x 3000 | `630e24bdd8e29193901dea8db2a4691fe0281043f11f2cc88264d9246c494a87` | copied PNG | Basic static example from package README. |
| `goal_plot_04.png` | resolved | `/ggoncoplot/man/figures/README-unnamed-chunk-3-1.png` | `README-unnamed-chunk-3-1.png` | ggoncoplot README/pkgdown marginal plots example | https://selkamand.github.io/ggoncoplot/ | ggoncoplot; MIT + file LICENSE; Sam El-Kamand | 7200 x 3600 | `40aebfb68a06b726dd0618c097bad0951a04c881d50f95b4821c9bdb093dd51c` | copied PNG | README marginal bar plots example. |
| `goal_plot_05.png` | resolved | `/ggoncoplot/man/figures/README-unnamed-chunk-4-1.png` | `README-unnamed-chunk-4-1.png` | ggoncoplot README/pkgdown clinical metadata example | https://selkamand.github.io/ggoncoplot/ | ggoncoplot; MIT + file LICENSE; Sam El-Kamand | 7200 x 3600 | `aac710c8a364b6b27fecedf483c9122fc293238bb127030fa1d8fed6deae5b2e` | copied PNG | README clinical metadata example. |
| `goal_plot_06.png` | resolved | `/ggoncoplot/paper/oncoplot_gbm.png` | `oncoplot_gbm.png` | ggoncoplot paper assets | https://github.com/selkamand/ggoncoplot/tree/main/paper | ggoncoplot; MIT + file LICENSE; Sam El-Kamand, Julian M. W. Quinn, Mark J. Cowley | 864 x 432 | `569397c0aab009dd92dabe077cac8de87d38e1ef3ee3c1807f3a1bd27726ac11` | copied PNG | GBM oncoplot paper asset. |
| `goal_plot_07.png` | resolved | `previous example_oncoplots/9c4c2e0f-74b4-4992-82be-57048160.png` | `9c4c2e0f-74b4-4992-82be-57048160.png` | Zhang et al. Annals of Translational Medicine 2019, Figure 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC6526269/ | Figure 1 from "Characterization of frequently mutated cancer genes in Chinese breast tumors: a comparison of Chinese and TCGA cohorts"; DOI 10.21037/atm.2019.04.23 | 850 x 683 | `be239584d040d28222dc7606c337ffdc118b90cdb26365ea57a77966ad6f18a8` | preserved existing PNG | Figure caption: "The mutational landscape of 33 cancer genes in Chinese breast tumors (n=305)." |
| `goal_plot_08.png` | resolved | `previous example_oncoplots/7a1ffb71-a4d9-4f42-8cf8-d992b037.png` | `7a1ffb71-a4d9-4f42-8cf8-d992b037.png` | Thind et al. Frontiers in Oncology 2022, Figure 2G | https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2022.919118/full | Figure panel from "Whole genome analysis reveals the genomic complexity in metastatic cutaneous squamous cell carcinoma"; Frontiers open-access CC BY; DOI 10.3389/fonc.2022.919118 | 1400 x 700 | `3c02828a7f23f6364828aebc42f25b9ce68825cbcb1ce88f05373abcc54573ad` | preserved existing PNG | Figure 2G is described as detailed sample-level SNV and variant-type information for top altered genes in metastatic CSCC. |
| `goal_plot_09.png` | resolved | `previous example_oncoplots/how-to-go-about-making-a-plot-like-this-v0-7ixn5xszvh3d1.webp` | `how-to-go-about-making-a-plot-like-this-v0-7ixn5xszvh3d1.webp` | Hoogstrate et al. Cancer Cell 2023 | https://www.sciencedirect.com/science/article/pii/S1535610823000478 | Oncoprint/clinical molecular panel from "Transcriptome analysis reveals tumor microenvironment changes in glioblastoma"; Cancer Cell / Elsevier open-access CC BY; DOI 10.1016/j.ccell.2023.02.019 | 1080 x 436 | `f4563ea506d1c5dd27de72b77c6ed7b20a6faa5d306d14a2173ea906b000f99f` | preserved existing converted PNG | User-identified source; article metadata also indexed by PubMed PMID 36898379 and Cardiff ORCA. |
| `goal_plot_10.png` | resolved | `previous example_oncoplots/oncoplot.png` | `oncoplot.png` | fuc tutorials, Create customized oncoplots | https://sbslee-fuc.readthedocs.io/en/latest/tutorials.html | fuc documentation by Seung-been Steven Lee; original license not bundled here | 1080 x 720 | `97e920006300cc26d5d301018863b4bd482f008061164fd832f4336972381710` | preserved existing PNG | Standard oncoplot tutorial image. |
| `goal_plot_11.png` | resolved | `previous example_oncoplots/customized_oncoplot_1.png` | `customized_oncoplot_1.png` | fuc tutorials, Create customized oncoplots | https://sbslee-fuc.readthedocs.io/en/latest/tutorials.html | fuc documentation by Seung-been Steven Lee; original license not bundled here | 1080 x 720 | `1589f9f9cc3899fb88738ea401f560c69b8f6ea648c44551f07bf7cd31885359` | preserved existing PNG | Customized oncoplot with annotations. |
| `goal_plot_12.png` | resolved | `previous example_oncoplots/customized_oncoplot_2.png` | `customized_oncoplot_2.png` | fuc tutorials, Create customized oncoplots | https://sbslee-fuc.readthedocs.io/en/latest/tutorials.html | fuc documentation by Seung-been Steven Lee; original license not bundled here | 1080 x 720 | `c838ebc206baa7ee3548472884fa1ce7c7baecbac1115803d757713cb8f33660` | preserved existing PNG | Customized oncoplot with sorted annotations. |
| `goal_plot_13.png` | resolved | `python_refactor_goal_sources/fuc_sources/customized_oncoplot_3.py` | `customized_oncoplot_3.png` | fuc tutorials, Create customized oncoplots | https://sbslee-fuc.readthedocs.io/en/latest/tutorials.html | fuc documentation by Seung-been Steven Lee; original license not bundled here | 1080 x 720 | `6e6f19675709e97c08f5c9e1326f96046775f73ab9b170c5dee7f5d2a9d97f85` | renamed user-provided PNG to numbered source plot | Customized oncoplot filtered to `Overall_Survival_Status == 1`. |

## Gallery Inputs

| Family | Files |
| --- | --- |
| AML/fuc | `aml_mutations.tsv`, `aml_metadata.tsv`, `aml_tmb.tsv`, `aml_palette.json`, `aml_gallery_params.json` |
| BRCA | `brca_mutations.tsv`, `brca_metadata.tsv`, `brca_tmb.tsv`, `brca_palette.json` |
| CSSC | `cssc_mutations.tsv`, `cssc_tmb.tsv`, `cssc_palette.json` |
| GBM | `gbm_clinical_tracks.tsv`, `gbm_events.tsv`, `gbm_palette.json` |
| ggoncoplot README | `ggoncoplot_readme_mutations.tsv`, `ggoncoplot_readme_metadata.tsv`, `ggoncoplot_readme_tmb.tsv`, `ggoncoplot_readme_palette.json` |

## Regeneration Commands

Regenerate deterministic synthetic input data for the non-fuc examples:

```bash
python3 python_refactor_goal_sources/generate_synthetic_inputs.py
```

Regenerate compact fuc-backed AML fixtures after downloading the upstream
`fuc-data` files listed in `fuc_sources/manifest.json`:

```bash
python3 python_refactor_goal_sources/fuc_sources/rebuild_fuc_fixtures.py --source-dir /path/to/fuc-data
```

Render clean generated plots:

```bash
python3 python_refactor_goal_sources/recreate_gallery.py --style clean
```

Render a BRCA comparison sheet:

```bash
python3 python_refactor_goal_sources/recreate_gallery.py --style comparison --preset brca_large
```
