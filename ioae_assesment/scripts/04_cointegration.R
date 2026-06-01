# =============================================================================
# 04_cointegration.R
# Cointegration Analysis: log(IGAE) — log(IMSS), Mexico 2004M1 onward
# =============================================================================
# Contents
#   0.  Libraries & data loading
#   1.  Integration order table (ADF — formal report version)
#   2.  Block-exogeneity Granger causality  (Toda-Yamamoto, VAR(p+1))
#   3.  Engle-Granger residual-based cointegration test
#   4.  Johansen FIML cointegration test    (trace + max-lambda)
#   5.  VECM: cointegrating vector beta + loading matrix alpha
#   6.  Long-run restrictions on beta       (unit-elasticity LR test)
#   7.  Weak exogeneity of log(IMSS)        (alpha restriction)
#   8.  Diagnostic plots                    (series + EG residuals)
# =============================================================================
#
# Data sources (read-only — no modifications to 01/02/03 scripts):
#   ../outputs/lmn/lmn_uncertainty_indices.csv   MU, FU (from 01_jln)
#   ../data/DFM/Datos.csv                        raw panel (IGAE, IMSS, etc.)
#   ../data/imef_recession_dates.xlsx            IMEF recession calendar
# =============================================================================

rm(list = ls())

##### 0. LIBRARIES & DATA LOADING #############################################

suppressMessages({
  library(urca)      # ur.df, ca.jo, cajorls, blrtest, alrtest
  library(vars)      # VAR, causality
  library(readr)     # read_csv
  library(readxl)    # read_xlsx
  library(dplyr)     # data manipulation
  library(tidyr)     # pivot_longer
  library(zoo)       # as.yearmon
  library(ggplot2)   # plots
  library(xtable)    # LaTeX table export
})

# ── Uncertainty indices ───────────────────────────────────────────────────────
ui <- read_csv("../outputs/lmn/lmn_uncertainty_indices.csv",
               show_col_types = FALSE) 
ui <- ui |> dplyr::select(date, MU, FU)

# ── Raw panel: IGAE, IMSS, PEDIDOS_MANU, IAI (same source as 01_jln) ─────────
raw <- read.csv("../data/DFM/Datos.csv", row.names = 1)
raw <- raw[rownames(raw) >= "2004/01",
           c("IGAE", "IMSS", "PEDIDOS_MANU", "IAI")]

dates_raw <- as.Date(as.yearmon(rownames(raw), "%Y/%m"))

df_raw <- as.data.frame(raw) |>
  mutate(
    date       = dates_raw,
    log_igae   = log(IGAE),
    log_imss   = log(IMSS),
    log_iai    = log(IAI),
    manuf_orders   = as.numeric(scale(PEDIDOS_MANU))
  )

# ── Merge on date → common sample ─────────────────────────────────────────────
# Drop trailing rows where IGAE is not yet published (NA in raw levels).
# IMSS and other series may be released earlier than IGAE.
df <- inner_join(df_raw, ui, by = "date") |>
  arrange(date) |>
  filter(!is.na(log_igae), !is.na(log_imss))

cat(sprintf("Common sample: %s to %s  (T = %d obs)\n",
            min(df$date), max(df$date), nrow(df)))

# ── Structural break dummies (impulse — same dates as in 01_jln) ──────────────
D_2009 <- as.integer(df$date == as.Date("2009-05-01"))
D_2020 <- as.integer(df$date == as.Date("2020-04-01"))

# ── Objects used throughout ───────────────────────────────────────────────────
# Unit root diagnostics (ADF + KPSS) classify all four series as I(1):
#   log(IGAE): ADF over-rejects due to COVID break; confirmed I(1) via ZA
#   log(IMSS): ADF(trend) does not reject; clearly I(1)
#   MU, FU   : KPSS(Level) and KPSS(Trend) both reject stationarity;
#              ADF on first differences strongly rejects → I(1) empirically
# → Expand the cointegrating system to all four variables (Option B).
# → Only structural break dummies remain as exogenous I(0) controls.

Y      <- cbind(log_igae = df$log_igae,
                log_imss = df$log_imss,
                MU       = df$MU,
                FU       = df$FU)                # T × 4  I(1) system

X_exog <- cbind(D_2009 = D_2009,
                D_2020 = D_2020)                 # T × 2  I(0) exogenous only

p_base <- 4    # VAR lag order (from main SVAR analysis)

# ── IMEF recession dates (for shading in plots) ───────────────────────────────
imef <- read_xlsx("../data/imef_recession_dates.xlsx") |>
  filter(phase == "recession", in_sample %in% c("Yes", "Partial")) |>
  mutate(
    start = as.Date(start),
    end   = if_else(is.na(end) | end == "", Sys.Date(), as.Date(end))
  )

rec_df <- data.frame(xmin = imef$start, xmax = imef$end,
                     ymin = -Inf,       ymax = Inf)

# ── Output directory ──────────────────────────────────────────────────────────
dir.create("../outputs/coint", recursive = TRUE, showWarnings = FALSE)

cat("Data loaded successfully.\n\n")


##### 1. INTEGRATION ORDER TABLE (ADF — FORMAL) ###############################
# Critical values (MacKinnon 1996, T ≈ 250):
#   with trend  (tau3): 1% = -4.04 | 5% = -3.45 | 10% = -3.15
#   with drift  (tau2): 1% = -3.46 | 5% = -2.88 | 10% = -2.57
#   with none   (tau1): 1% = -2.60 | 5% = -1.95 | 10% = -1.61
#
# NOTE on integration orders:
#
# log(IGAE): standard ADF over-rejects (tau3 = -4.79) due to a large level
#   shift. ZA (model="both") endogenously finds a break at 2019M2 (pre-COVID
#   Mexican recession) with tau_ZA = -5.46, clearing the 5% CV (-5.08) but
#   NOT the 1% CV (-5.57). Treated as I(1): borderline ZA, economic theory,
#   and first-difference ADF (tau2 = -12.6) all support this.
#
# log(IMSS): ZA tau = -3.18 (break ≈ 2018M11); cannot reject unit root at
#   any level. Clearly I(1).
#
# MU:  ZA tau = -4.50 (break ≈ 2019M1);  > 10% CV (-4.82) — I(1) confirmed.
# FU:  ZA tau = -4.33 (break ≈ 2010M2);  > 10% CV (-4.82) — I(1) confirmed.
#   (Break at 2010M2 captures post-GFC financial recovery.)
#   KPSS(Level) and KPSS(Trend) both reject stationarity; ADF on differences
#   strongly rejects. High persistence typical of GARCH volatility processes.

adf_row <- function(x, name, type_lev = "trend", override_order = NULL) {
  lev <- ur.df(na.omit(x),       type = type_lev, selectlags = "AIC")
  dif <- ur.df(diff(na.omit(x)), type = "drift",  selectlags = "AIC")
  cv5_lev <- if (type_lev == "trend") -3.45 else if (type_lev == "drift") -2.88 else -1.95
  auto_order <- if (lev@teststat[1] > cv5_lev && dif@teststat[1] < -2.88)
    "I(1)" else "I(0)"
  order <- if (!is.null(override_order)) override_order else auto_order
  data.frame(
    Series  = name,
    k_lev   = lev@lags,
    tau_lev = round(lev@teststat[1], 3),
    k_dif   = dif@lags,
    tau_dif = round(dif@teststat[1], 3),
    Order   = order,
    stringsAsFactors = FALSE
  )
}

int_table <- do.call(rbind, list(
  adf_row(df$log_igae, "log(IGAE)", type_lev = "trend", override_order = "I(1)$^\\dag$"),
  adf_row(df$log_imss, "log(IMSS)", type_lev = "trend"),
  adf_row(df$MU,       "MU",        type_lev = "drift", override_order = "I(1)$^\\ddag$"),
  adf_row(df$FU,       "FU",        type_lev = "drift", override_order = "I(1)$^\\ddag$")
))

cat("======= 1. INTEGRATION ORDER TABLE =======\n")
print(int_table, row.names = FALSE)

# Export to LaTeX
int_xtab <- xtable(
  int_table,
  caption = paste0(
    "ADF unit-root tests. $k$ = lag length (AIC). ",
    "log(IGAE) and log(IMSS) tested with trend ($\\tau_3$); ",
    "MU and FU tested with drift ($\\tau_2$). ",
    "Critical values (5\\%): $\\tau_3=-3.45$, $\\tau_2=-2.88$ (MacKinnon 1996). ",
    "$^\\dag$~Standard ADF over-rejects for log(IGAE) ($\\tau_3=-4.79$) due to a ",
    "level shift at the 2019M2 pre-COVID Mexican recession. ",
    "Zivot-Andrews (1992, model = both): $\\tau_{ZA}=-5.46$ at 2019M2, ",
    "clearing the 5\\% CV ($-5.08$) but not the 1\\% CV ($-5.57$); ",
    "treated as I(1) on economic grounds and by preponderance of evidence. ",
    "log(IMSS): $\\tau_{ZA}=-3.18$ (break 2018M11) confirms I(1). ",
    "$^\\ddag$~KPSS(Level) and KPSS(Trend) both reject stationarity ($p<0.01$). ",
    "ZA: $\\tau_{ZA}=-4.50$ for MU (break 2019M1) and $\\tau_{ZA}=-4.33$ for FU ",
    "(break 2010M2), both above the 10\\% CV ($-4.82$), confirming I(1). ",
    "High persistence of GARCH-based volatility is discussed in ",
    "Bollerslev et al.\\ (1994)."
  ),
  label  = "tab:adf_coint",
  digits = c(0, 0, 0, 3, 0, 3, 0)
)
colnames(int_xtab) <- c("Series", "$k$", "$\\tau$ (levels)",
                          "$k$", "$\\tau$ ($\\Delta$)", "Order")
print(int_xtab,
      file                   = "../outputs/coint/tab_adf_coint.tex",
      sanitize.text.function = identity,
      include.rownames       = FALSE,
      booktabs               = TRUE,
      caption.placement      = "top")
cat("Saved: outputs/coint/tab_adf_coint.tex\n\n")

library(tseries)
kpss.test(df$MU, null = "Level")   # H0: stationary (level)
kpss.test(df$FU, null = "Level")   # H0: stationary (level)


##### 2. BLOCK-EXOGENEITY GRANGER CAUSALITY (TODA-YAMAMOTO) ###################
# Toda & Yamamoto (1995): fit VAR(p + d_max) in levels; test only first p lags.
# p_base = 4, d_max = 1 → VAR(5). Wald statistic ~ chi^2 asymptotically.

d_max   <- 1
var_ty  <- VAR(
  cbind(log_igae = df$log_igae,
        log_imss = df$log_imss,
        MU       = df$MU,
        FU       = df$FU),
  p      = p_base + d_max,     # VAR(5) for T-Y
  type   = "const",
  exogen = cbind(D_2009 = D_2009,
                 D_2020 = D_2020)
)

cat("======= 2. BLOCK-EXOGENEITY GRANGER CAUSALITY (TODA-YAMAMOTO, VAR(5)) =======\n")

# (a) {MU, FU} → {log_igae, log_imss}
cat("\n[a] H0: {MU, FU} do NOT Granger-cause {log_igae, log_imss}\n")
gc_mu_fu  <- causality(var_ty, cause = c("MU", "FU"))
print(gc_mu_fu$Granger)
cat("Instantaneous causality:\n"); print(gc_mu_fu$Instant)

# (b) log_igae → log_imss (within I(1) block)
cat("\n[b] H0: log_igae does NOT Granger-cause log_imss / MU / FU\n")
gc_igae   <- causality(var_ty, cause = "log_igae")
print(gc_igae$Granger)

# (c) log_imss → log_igae (within I(1) block)
cat("\n[c] H0: log_imss does NOT Granger-cause log_igae / MU / FU\n")
gc_imss   <- causality(var_ty, cause = "log_imss")
print(gc_imss$Granger)

cat("\n")

# ── Export Granger causality table ────────────────────────────────────────────
fmt_p <- function(p) {
  stars <- ifelse(p < 0.01, "$^{***}$", ifelse(p < 0.05, "$^{**}$",
           ifelse(p < 0.10, "$^{*}$",   "")))
  paste0(formatC(p, format = "f", digits = 4), stars)
}

granger_tab <- data.frame(
  Hypothesis = c(
    "$\\{MU,FU\\} \\not\\rightarrow \\{\\log\\text{IGAE},\\log\\text{IMSS}\\}$",
    "$\\log\\text{IGAE} \\not\\rightarrow \\{\\log\\text{IMSS}, MU, FU\\}$",
    "$\\log\\text{IMSS} \\not\\rightarrow \\{\\log\\text{IGAE}, MU, FU\\}$"
  ),
  Test      = c("$F$", "$F$", "$F$"),
  Statistic = c(
    formatC(gc_mu_fu$Granger$statistic, format = "f", digits = 4),
    formatC(gc_igae$Granger$statistic,  format = "f", digits = 4),
    formatC(gc_imss$Granger$statistic,  format = "f", digits = 4)
  ),
  df1 = c(
    gc_mu_fu$Granger$parameter["df1"],
    gc_igae$Granger$parameter["df1"],
    gc_imss$Granger$parameter["df1"]
  ),
  df2 = c(
    gc_mu_fu$Granger$parameter["df2"],
    gc_igae$Granger$parameter["df2"],
    gc_imss$Granger$parameter["df2"]
  ),
  p_value = c(
    fmt_p(gc_mu_fu$Granger$p.value),
    fmt_p(gc_igae$Granger$p.value),
    fmt_p(gc_imss$Granger$p.value)
  ),
  stringsAsFactors = FALSE
)

# Instantaneous causality row (chi-squared test)
instant_row <- data.frame(
  Hypothesis = "$\\{MU,FU\\}$ no instantaneous causality with $\\{\\log\\text{IGAE},\\log\\text{IMSS}\\}$",
  Test       = "$\\chi^2$",
  Statistic  = formatC(gc_mu_fu$Instant$statistic, format = "f", digits = 4),
  df1        = gc_mu_fu$Instant$parameter,
  df2        = NA_real_,
  p_value    = fmt_p(gc_mu_fu$Instant$p.value),
  stringsAsFactors = FALSE
)

granger_full <- rbind(granger_tab)

dir.create("../outputs/coint", recursive = TRUE, showWarnings = FALSE)

gc_xtab <- xtable(
  granger_full,
  caption = paste0(
    "Block-exogeneity Granger causality tests. ",
    "Toda-Yamamoto (1995) procedure: VAR(", p_base + d_max, ") in levels; ",
    "exogenous: $D_{2009}$, $D_{2020}$. ",
    "Granger tests use an $F$-statistic; instantaneous causality uses $\\chi^2$. ",
    "Significance: $^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."
  ),
  label  = "tab:granger"
)
colnames(gc_xtab) <- c("$H_0$", "Test", "Statistic", "$df_1$", "$df_2$", "$p$-value")

print(gc_xtab,
      file                   = "../outputs/coint/tab_granger.tex",
      sanitize.text.function = identity,
      include.rownames       = FALSE,
      booktabs               = TRUE,
      caption.placement      = "top",
      NA.string              = "---")
cat("Saved: outputs/coint/tab_granger.tex\n\n")


##### 3. ENGLE-GRANGER RESIDUAL-BASED COINTEGRATION TEST ######################
# Full 4-variable system: normalise on log(IGAE).
# OLS regression of log(IGAE) on log(IMSS), MU, FU → test residuals for I(0).
# MacKinnon (1991) response-surface CVs for K = 4 variables (3 regressors):
#   1% = -4.64   5% = -4.10   10% = -3.81
# Note: EG tests for at least one CI vector but cannot determine rank r.
#       Johansen (§4) provides the full rank determination.

eg_ols <- lm(log_igae ~ log_imss + MU + FU, data = df)
eg_res <- residuals(eg_ols)

eg_adf <- ur.df(eg_res, type = "none", selectlags = "AIC")

cat("======= 3. ENGLE-GRANGER RESIDUAL-BASED TEST (K=4) =======\n")
cat(sprintf("Cointegrating regression:\n"))
cat(sprintf("  log(IGAE) = %.4f + %.4f*log(IMSS) + %.4f*MU + %.4f*FU\n",
            coef(eg_ols)[1], coef(eg_ols)[2],
            coef(eg_ols)[3], coef(eg_ols)[4]))
cat(sprintf("ADF on residuals: tau = %.4f  (lags = %d)\n",
            eg_adf@teststat[1], eg_adf@lags))
cat("MacKinnon (1991) CVs [K=4]:  1% = -4.64  |  5% = -4.10  |  10% = -3.81\n")
cat(sprintf("Decision (5%%): %s\n\n",
            ifelse(eg_adf@teststat[1] < -4.10,
                   "REJECT H0 — at least one CI vector found",
                   "FAIL TO REJECT H0 — no cointegration")))


##### 4. JOHANSEN FIML COINTEGRATION TEST #####################################
# 4-variable system: Y = (log_igae, log_imss, MU, FU)'.
# ecdet = "const": constant restricted to CI space.
# K = p_base = 4 (levels VAR lag order).
# dumvar: break dummies only (D_2009, D_2020).
# Tests r = 0, 1, 2, 3 sequentially; stop at first non-rejection.

jo_trace <- ca.jo(Y, type = "trace", ecdet = "const", K = p_base,
                  dumvar = X_exog)
jo_eigen <- ca.jo(Y, type = "eigen", ecdet = "const", K = p_base,
                  dumvar = X_exog)

cat("======= 4. JOHANSEN FIML (K=4 variables) =======\n")
cat("\n--- Trace test ---\n");        print(summary(jo_trace))
cat("\n--- Max-eigenvalue test ---\n"); print(summary(jo_eigen))

# Compact table: 4 hypotheses, both tests side by side
# jo_trace@teststat is ordered r=K-1,...,0 (bottom to top); reverse for display
n_hyp   <- length(jo_trace@teststat)
h0_labs <- paste0("$r \\leq ", seq(0, n_hyp - 1), "$")
h0_labs[1] <- "$r = 0$"

jo_tab <- data.frame(
  H0      = h0_labs,
  tr_stat = round(rev(jo_trace@teststat), 3),
  tr_cv5  = rev(jo_trace@cval[, "5pct"]),
  me_stat = round(rev(jo_eigen@teststat), 3),
  me_cv5  = rev(jo_eigen@cval[, "5pct"]),
  stringsAsFactors = FALSE
)

cat("\n--- Compact table ---\n")
print(jo_tab, row.names = FALSE)

# Export to LaTeX
jo_xtab <- xtable(
  jo_tab,
  caption = paste0("Johansen cointegration tests. VAR lag $K=4$, ",
                   "constant restricted to CI space (\\texttt{ecdet=''const''}). ",
                   "Exogenous: MU, FU, PEDIDOS\\_MANU, $D_{2009}$, $D_{2020}$. ",
                   "Critical values from Osterwald-Lenum (1992)."),
  label   = "tab:johansen",
  digits  = c(0, 0, 3, 2, 3, 2)
)
colnames(jo_xtab) <- c("$H_0$",
                        "$\\lambda_{\\text{trace}}$", "CV 5\\%",
                        "$\\lambda_{\\max}$",          "CV 5\\%")
print(jo_xtab,
      file                   = "../outputs/coint/tab_johansen.tex",
      sanitize.text.function = identity,
      include.rownames       = FALSE,
      booktabs               = TRUE,
      caption.placement      = "top")
cat("Saved: outputs/coint/tab_johansen.tex\n\n")


##### 5. VECM: COINTEGRATING VECTOR + LOADING MATRIX ##########################
# Set r from Johansen results above. Update after inspecting §4 output.
# jo_trace@V is (K+1) × K with K=4; rows = [log_igae, log_imss, MU, FU, const]
# Columns = eigenvectors (CI vectors), ordered by eigenvalue size.

r    <- 1    # <-- UPDATE after running §4 if Johansen suggests r ≠ 1
vecm <- cajorls(jo_trace, r = r)

# --- Cointegrating vector(s) beta ---
# Normalise each CI vector on log_igae (first row)
beta_raw <- jo_trace@V[, 1:r, drop = FALSE]
beta_hat <- sweep(beta_raw, 2, beta_raw[1, ], FUN = "/")
rownames(beta_hat) <- c("log_igae", "log_imss", "MU", "FU", "const")

cat("======= 5. VECM (r =", r, ") =======\n")
cat("\n--- Cointegrating vector(s) beta (normalised on log_igae) ---\n")
print(round(beta_hat, 4))
cat(sprintf("\n  Okun elasticity (beta_IMSS): %.4f\n",  beta_hat["log_imss", 1]))
cat(sprintf("  MU long-run coefficient:      %.4f\n",  beta_hat["MU", 1]))
cat(sprintf("  FU long-run coefficient:      %.4f\n",  beta_hat["FU", 1]))

# --- Loading matrix alpha (4 × r) ---
cat("\n--- Loading matrix alpha (speed-of-adjustment) ---\n")
ect_rows  <- paste0("ect", seq_len(r))
alpha_hat <- vecm$rlm$coefficients[ect_rows, , drop = FALSE]
print(round(alpha_hat, 4))
cat("  Negative alpha_igae => IGAE error-corrects toward long-run equilibrium.\n")
cat("  Insignificant alpha_X => variable X is weakly exogenous.\n")

# --- Full VECM summary ---
cat("\n--- VECM full equation summaries ---\n")
vecm_summ <- summary(vecm$rlm)
print(vecm_summ)

# Helper: significance stars
stars <- function(p) ifelse(p < 0.01, "$^{***}$",
                     ifelse(p < 0.05, "$^{**}$",
                     ifelse(p < 0.10, "$^{*}$", "")))

# --- Save full output to text ---
sink("../outputs/coint/vecm_summary.txt")
cat("VECM — cajorls(jo_trace, r =", r, ")\n")
cat("System: log_igae, log_imss, MU, FU\n\n")
cat("Cointegrating vector(s) beta (normalised on log_igae):\n")
print(round(beta_hat, 4))
cat("\nLoading matrix alpha:\n")
print(round(alpha_hat, 4))
cat("\nFull coefficient table:\n")
print(vecm_summ)
sink()
cat("Saved: outputs/coint/vecm_summary.txt\n\n")


##### 6. LONG-RUN RESTRICTIONS ON BETA ########################################
# K+1 = 5 rows in V (4 variables + restricted constant).
# Two economically motivated restrictions on the first CI vector:
#
# Restriction A — Exclusion of MU from the long-run relation
#   H0: beta_MU = 0  (uncertainty does not enter the CI vector;
#       it only shifts short-run dynamics)
#   H = 5×4 matrix: all columns of the 5-dim identity EXCEPT row 3 (MU)
#
# Restriction B — Exclusion of FU from the long-run relation
#   H0: beta_FU = 0  (symmetric argument for financial uncertainty)
#
# Restriction C — Unit Okun elasticity: beta_IMSS = -1

cat("======= 6. LONG-RUN RESTRICTIONS =======\n")

# ── Restriction A: beta_MU = 0 ────────────────────────────────────────────────
# beta = H*phi;  H is 5×4 (remove row 3 = MU from free parameters)
# Free: log_igae, log_imss, FU, const  → s=4, df=(5-4)*r=1
H_noMU <- cbind(c(1,0,0,0,0),   # log_igae free
                c(0,1,0,0,0),   # log_imss free
                c(0,0,0,1,0),   # FU free
                c(0,0,0,0,1))   # const free
lr_noMU <- blrtest(jo_trace, H = H_noMU, r = r)
cat("\n[A] H0: beta_MU = 0 (MU excluded from CI vector)\n")
print(summary(lr_noMU))

# ── Restriction B: beta_FU = 0 ────────────────────────────────────────────────
# Free: log_igae, log_imss, MU, const → s=4, df=1
H_noFU <- cbind(c(1,0,0,0,0),
                c(0,1,0,0,0),
                c(0,0,1,0,0),
                c(0,0,0,0,1))
lr_noFU <- blrtest(jo_trace, H = H_noFU, r = r)
cat("\n[B] H0: beta_FU = 0 (FU excluded from CI vector)\n")
print(summary(lr_noFU))

# ── Restriction C: unit Okun elasticity beta_IMSS = -1 ───────────────────────
# Free: log_igae, FU, MU, const; beta_imss fixed at -1
# H = 5×4: col1=(1,-1,0,0,0)', col2=(0,0,1,0,0)', col3=(0,0,0,1,0)', col4=(0,0,0,0,1)'
H_unit <- cbind(c(1,-1,0,0,0),
                c(0, 0,1,0,0),
                c(0, 0,0,1,0),
                c(0, 0,0,0,1))
lr_unit <- blrtest(jo_trace, H = H_unit, r = r)
cat("\n[C] H0: beta_IMSS = -1 (unit Okun elasticity)\n")
cat(sprintf("    Unrestricted beta_IMSS = %.4f\n", beta_hat["log_imss", 1]))
print(summary(lr_unit))


##### 7. WEAK EXOGENEITY TESTS #################################################
# alrtest: H0: alpha ∈ col(A), i.e. some rows of alpha are zero.
# K=4, so A is 4×s; restricting one equation sets s=3, df=(4-3)*r=r.
#
# Test each variable for weak exogeneity in turn.
# A variable is weakly exogenous if its alpha row = 0:
#   it does not respond to deviations from the long-run equilibrium.

cat("======= 7. WEAK EXOGENEITY TESTS =======\n")

vars_K  <- c("log_igae", "log_imss", "MU", "FU")
for (i in seq_along(vars_K)) {
  # A = identity(4) with row i removed → s = 3 free rows
  A_i <- diag(4)[, -i, drop = FALSE]
  we  <- alrtest(jo_trace, A = A_i, r = r)
  cat(sprintf("\n  H0: alpha_%s = 0 (weakly exogenous)\n", vars_K[i]))
  cat(sprintf("  LR = %.4f,  df = %d,  p = %.4f\n",
              we@teststat, we@df, pchisq(we@teststat, we@df, lower.tail = FALSE)))
}


##### 8. DIAGNOSTIC PLOTS #####################################################

# ── Fig C1: all four I(1) series standardised — two panels ───────────────────
# Panel A: real activity (log_igae, log_imss)
# Panel B: uncertainty (MU, FU)

df_real <- df |>
  mutate(across(c(log_igae, log_imss), ~ (. - mean(.)) / sd(.))) |>
  select(date, log_igae, log_imss) |>
  pivot_longer(-date, names_to = "series", values_to = "value") |>
  mutate(series = recode(series, log_igae = "log(IGAE)", log_imss = "log(IMSS)"))

df_unc <- df |>
  mutate(across(c(MU, FU), ~ (. - mean(.)) / sd(.))) |>
  select(date, MU, FU) |>
  pivot_longer(-date, names_to = "series", values_to = "value")

theme_coint <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_real <- ggplot(df_real, aes(x = date, y = value, colour = series)) +
  geom_rect(data = rec_df,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            inherit.aes = FALSE, fill = "grey70", alpha = 0.3) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c("log(IGAE)" = "steelblue",
                                  "log(IMSS)" = "firebrick")) +
  labs(title = "Real activity", x = NULL,
       y = "Standardised level", colour = NULL) +
  theme_coint

p_unc <- ggplot(df_unc, aes(x = date, y = value, colour = series)) +
  geom_rect(data = rec_df,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            inherit.aes = FALSE, fill = "grey70", alpha = 0.3) +
  geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c("MU" = "darkorange", "FU" = "purple4")) +
  labs(title = "Uncertainty indices", x = NULL,
       y = "Standardised level", colour = NULL) +
  theme_coint

library(patchwork)
p_series <- (p_real / p_unc) +
  plot_annotation(
    title    = "I(1) system: log(IGAE), log(IMSS), MU, FU — Mexico",
    subtitle = "Standardised; shaded = IMEF recession periods",
    theme    = theme(plot.title = element_text(face = "bold"))
  )

ggsave("../outputs/coint/fig_coint_series.png",
       p_series, width = 10, height = 6, dpi = 150)
cat("Saved: outputs/coint/fig_coint_series.png\n")

# ── Fig C2: Engle-Granger residuals + zero line + recession shading ───────────
eg_res_df <- data.frame(date = df$date, residual = eg_res)

p_eg <- ggplot(eg_res_df, aes(x = date, y = residual)) +
  geom_rect(data = rec_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey70", alpha = 0.3) +
  geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed", linewidth = 0.6) +
  geom_line(colour = "steelblue", linewidth = 0.8) +
  labs(
    title    = "Engle-Granger cointegrating residuals",
    subtitle = expression(hat(z)[t] == log(IGAE)[t] - hat(alpha) - hat(beta) * log(IMSS)[t]),
    x        = NULL,
    y        = expression(hat(z)[t])
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("../outputs/coint/fig_coint_residuals.png",
       p_eg, width = 10, height = 3.5, dpi = 150)
cat("Saved: outputs/coint/fig_coint_residuals.png\n")

# ── Fig C3: ECT (VECM equilibrium error) over time ───────────────────────────
# ECT = beta' * [Y; 1] — 4 variables + restricted constant
# beta_hat has rows: log_igae, log_imss, MU, FU, const
ect_vals <- Y %*% beta_hat[1:4, 1] + beta_hat["const", 1]
ect_df   <- data.frame(date = df$date, ect = as.numeric(ect_vals))

p_ect <- ggplot(ect_df, aes(x = date, y = ect)) +
  geom_rect(data = rec_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey70", alpha = 0.3) +
  geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed", linewidth = 0.6) +
  geom_line(colour = "darkorange", linewidth = 0.8) +
  labs(
    title    = "VECM error-correction term (ECT)",
    subtitle = expression(ECT[t] == hat(beta)^"'" * Y[t] + hat(mu)),
    x        = NULL,
    y        = expression(ECT[t])
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("../outputs/coint/fig_coint_ect.png",
       p_ect, width = 10, height = 3.5, dpi = 150)
cat("Saved: outputs/coint/fig_coint_ect.png\n")

cat("\n===== 04_cointegration.R complete =====\n")
cat("Outputs written to: outputs/coint/\n")
cat("  tab_adf_coint.tex       — ADF integration order table (4 variables)\n")
cat("  tab_granger.tex         — Toda-Yamamoto block-exogeneity table\n")
cat("  tab_johansen.tex        — Johansen trace + max-lambda table\n")
cat("  vecm_summary.txt        — Full VECM coefficient output\n")
cat("  fig_coint_series.png    — Real activity + uncertainty time series\n")
cat("  fig_coint_residuals.png — EG residuals (4-variable)\n")
cat("  fig_coint_ect.png       — VECM error-correction term\n")
