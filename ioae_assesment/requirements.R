# requirements.R
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Helper: install CRAN pkg if missing
ensure_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Helper: install GitHub pkg if missing
ensure_github <- function(pkg, repo_spec) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    ensure_cran("pak")
    pak::pkg_install(repo_spec)
  }
}

# --- CRAN dependencies (unique list) ---
cran_pkgs <- c(
  "remotes",
  "here",
  "seasonal",
  "imputeTS",
  "glmnet",
  "forecast",
  "caret",
  "nnfor",
  "vars",
  "tidyverse",
  "murphydiagram",
  "MCS",
  "gt",
  "zoo",
  "glue",
  "rmarkdown"
)

for (p in cran_pkgs) ensure_cran(p)

# --- GitHub dependency: nowcasting (pinned) ---
# Repo: https://github.com/nmecsys/nowcasting
# Pinning avoids "moving target" installs
remotes::install_github("nmecsys/nowcasting")

# Sanity checks
stopifnot(requireNamespace("nowcasting", quietly = TRUE))
