# Replication Package — IOAE performance assesment

This replication package accompanies:

**Evaluating the predictive power of Mexico’s Timely Economic Activity Indicator: Real and pseudo real time performance**  
Francisco Corona, Jesús López-Pérez, Edgar René Benavidez-Maruri, and Francisco J. López-Ballesteros  
Manuscript under review

## Overview

The goal of this project is to reproduce the analyses and outputs reported in the manuscript using R. Dependencies are installed automatically by the replication pipeline to facilitate execution across different machines.

The replication pipeline produces:
- Intermediate objects and model outputs stored as `.rds` files under `outputs/rds/`
- Final figures saved as `.png` files under `outputs/figures/`
- An optional HTML viewer (`outputs/outputs.html`) rendered from `analysis/04_outputs.Rmd` for convenient inspection of final outputs

## Quick start

```r
source("analysis/run_all.R")
````

## Requirements

* **R 4.5.1**
* An internet connection is required the first time you run the pipeline (dependencies are installed automatically).
* The pipeline installs required R packages from **CRAN** and installs `nowcasting` from **GitHub** (`nmecsys/nowcasting`).
* **RMarkdown** is optional and only required to render the outputs viewer document (`analysis/04_outputs.Rmd`).
* **System dependency note (seasonal adjustment):** this project uses `seasonal::seas()` for seasonal adjustment. Depending on the operating system and local configuration, `seasonal` may require the X-13ARIMA-SEATS binaries.

## Distribution

This replication package is distributed as a compressed archive (e.g., `.zip`/`.tar.gz`). After extraction, run the workflow from the project root directory following the instructions below.

## Archive contents

The archive is expected to contain at least: `requirements.R`, `ioae_assesmente.Rproj`, `analysis/`, `R/`, and `data/` (if data are included), as well as this `README.md`.

## Setup

1. Download the replication package as a compressed archive (e.g., `.zip`/`.tar.gz`) and extract it to a local folder.

2. Open an R session (e.g., RStudio or the R console).

   **Recommended (RStudio):** open the project by double-clicking `ioae_assesment.Rproj` (located in the extracted project folder). This helps ensure paths are set correctly.

3. Set the **R working directory** to the **project root** (the folder that contains `requirements.R`, `analysis/`, and `data/`).

   You can verify you are in the correct folder by running:

```r
file.exists("requirements.R") && dir.exists("analysis")
```

4. Run the full replication pipeline:

```r
source("analysis/run_all.R")
```

## Optional: outputs viewer (`analysis/04_outputs.Rmd`)

The file `analysis/04_outputs.Rmd` is provided only as a convenient way to view the final tables and figures. It is **not required** to reproduce the underlying results.

### Prerequisites

Rendering this `.Rmd` requires **Pandoc >= 1.12.3**. The simplest way to ensure this requirement is met is to use a **recent version of RStudio** (which typically bundles a compatible Pandoc).

You can check Pandoc availability from R:

```r
rmarkdown::pandoc_available("1.12.3")
rmarkdown::pandoc_version()
```

### Recommended workflow (RStudio)

1. Open the project by double-clicking `ioae_assesmente.Rproj` (located in the extracted project folder).
2. Run the replication pipeline first (this creates `outputs/` and its contents):

```r
source("analysis/run_all.R")
```

3. Confirm that outputs exist, then render the viewer:

```r
stopifnot(dir.exists(here::here("outputs")))

rmarkdown::render(
  input = here::here("analysis", "04_outputs.Rmd"),
  output_file = "outputs.html",
  output_dir = here::here("outputs"),
  quiet = TRUE
)
```

The rendered HTML file will be saved as `outputs/outputs.html`.

## Project structure

The replication package is organized as follows:

* `analysis/`
  Analysis scripts and the main entry point:

  * `run_all.R` (master script; runs the full pipeline and installs dependencies)
  * `01_preparing.R` (data loading + seasonal adjustment + transformations)
  * `02_models.R` (nowcasting scenarios + model evaluation + RDS outputs)
  * `03_figures_tables.R` (generates final figures and tables)
  * `04_outputs.Rmd` (optional viewer document; renders `outputs/outputs.html`)

* `data/`
  Raw input time-series data (CSV) used by the pipeline.

* `R/`
  Project functions sourced by the pipeline.

* `outputs/`
  All generated artifacts (figures, stored objects, and the HTML viewer).

## Replication workflow

Running `analysis/run_all.R` executes the following steps:

1. **Prepare inputs (`analysis/01_preparing.R`)**

   * Loads raw time-series inputs from `data/`
   * Applies seasonal adjustment where required (`seasonal::seas()`)
   * Applies transformations and constructs prepared inputs used downstream

2. **Modeling and evaluation (`analysis/02_models.R`)**

   * Executes the scenario loop (e.g., varying model/hyperparameter settings and windowing strategy)
   * Produces nowcasts and computes evaluation metrics
   * Writes intermediate and final results to disk as `.rds` objects

3. **Figures and tables (`analysis/03_figures_tables.R`)**

   * Generates final figures and tables based on stored `.rds` results
   * Writes final figure files (`.png`) to `outputs/figures/`
   * Writes table objects (`.rds`) to `outputs/rds/prepared_tables/`

4. **Outputs viewer (optional) (`analysis/04_outputs.Rmd`)**

   * Renders `outputs/outputs.html`, a convenient way to view final tables and figures

## Expected outputs

After a successful run, the following artifacts will be created under `outputs/`:

### Figures (`.png`)

Saved under `outputs/figures/`, including (non-exhaustive):

* Main figures: `figure_2.png`, `figure_3.png`, `figure_5.png`, `figure_6.png`, `figure_7.png`
* Black-and-white versions: `figure_2_bw.png`, `figure_3_bw.png`, `figure_5_bw.png`, `figure_6_bw.png`, `figure_7_bw.png`
* Appendix figures: `figure_app_1.png`–`figure_app_4.png` and corresponding `_bw` versions

### Saved objects (`.rds`)

Saved under `outputs/rds/` and organized by purpose:

* `outputs/rds/prepared_inputs/`

  * `prepared_inputs.rds` (consolidated prepared inputs used downstream)

* `outputs/rds/models_objects/`

  * Model and evaluation objects (e.g., factor counts, error objects, MCS objects, p-values, test outputs)

* `outputs/rds/figures_inputs/`

  * Inputs required to build final figures (e.g., `f5.rds`, `f6.rds`, `f7.rds`, `TS.rds`)
  * Performance matrices (e.g., `mat_mae_*`, `mat_rmse_*`)

* `outputs/rds/prepared_tables/`

  * Tables stored as R objects (e.g., `t1.rds`–`t4.rds`, `t_dm_rmse.rds`, `t_mcs_rmse.rds`, `t_n_factors.rds`, `t_spa_test.rds`)

### HTML viewer (optional)

* `outputs/outputs.html` (only if `analysis/04_outputs.Rmd` is rendered)

## Reproducibility notes

* Dependencies are installed automatically by the pipeline (via `requirements.R`). For a clean run, start from a fresh R session.
* Some steps involve time-series transformations and seasonal adjustment. Results should be deterministic given a fixed codebase and the seed set by the pipeline.
* The replication workflow writes intermediate objects to `outputs/rds/` and final figures to `outputs/figures/`.

## Troubleshooting

### Dependency installation (CRAN / GitHub)

The pipeline installs dependencies automatically when you run `analysis/run_all.R`. If package installation fails:

1. Make sure you have an active internet connection.
2. Re-run the pipeline in a fresh R session.
3. If you are behind a corporate proxy/firewall, ensure that CRAN and GitHub are accessible from your network.
4. If the error message indicates missing system build tools (compilation), install the required tools for your operating system (e.g., Rtools on Windows; Xcode Command Line Tools on macOS), then re-run.

#### GitHub package `nowcasting`

This project installs `nowcasting` from GitHub (`nmecsys/nowcasting`) because it is not available for some R versions on CRAN. If GitHub installation fails, verify that access to GitHub is allowed on your network and try again.

### `seasonal::seas()` / X-13ARIMA-SEATS issues

This project uses `seasonal::seas()` for seasonal adjustment. On some machines, `seas()` may fail if the required X-13ARIMA-SEATS binaries are not available or not properly configured.

If you encounter an error at the seasonal-adjustment step:

1. Confirm that the `seasonal` package installed successfully (the pipeline installs it automatically).
2. Re-run the pipeline to verify the error is reproducible.
3. If the error persists, install/configure X-13ARIMA-SEATS for your operating system (per the `seasonal` package guidance) and re-run the pipeline.

If you are running on a restricted/corporate environment where installing external binaries is not possible, please contact the maintainer (see Contact below) for alternatives.

### Rendering the outputs document (`analysis/04_outputs.Rmd`)

Rendering is **optional** and may fail on machines with an outdated Pandoc installation. The core replication outputs (figures and `.rds` objects under `outputs/`) are still produced by the pipeline even if the HTML viewer cannot be rendered.

If you want to render `analysis/04_outputs.Rmd`, make sure you are using **Pandoc >= 1.12.3** (a recent RStudio version typically satisfies this), then run:

```r
rmarkdown::render(
  input = here::here("analysis", "04_outputs.Rmd"),
  output_file = "outputs.html",
  output_dir = here::here("outputs")
)
```

### Working directory / project root

If the pipeline fails because it cannot find files (e.g., under `data/` or `analysis/`), ensure your **R working directory** is set to the project root (the folder that contains `requirements.R`). You can check:

```r
file.exists("requirements.R") && dir.exists("analysis") && dir.exists("data")
```

## Citation

If you use this replication package, please cite the associated manuscript:

Corona, F., López-Pérez, J., Benavidez-Maruri, E. R., & López-Ballesteros, F. J. (Manuscript under review). *Evaluating the predictive power of Mexico’s Timely Economic Activity Indicator: Real and pseudo real time performance*.

## License

All rights reserved.

## Contact

For questions, issues, or replication support, please contact:

[franciscoj.corona@inegi.org.mx](mailto:franciscoj.corona@inegi.org.mx)

