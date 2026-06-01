rm(list= ls())
# =============================================================================
# LMN Uncertainty Indices for the Mexican Economy
# =============================================================================
# Built on top of Corona et al. (IOAE) pipeline:
# cloned from https://github.com/FranciscoL-B/ioae_assesment/tree/main
#   - Stationarity transformations: their best_trans() / cpm() logic
#   - Seasonal adjustment:          seasonal::seas()
#   - DFM factor extraction:        nowcasting::nowcast(), method = "2s"
#   - Factor count:                 nowcasting::bai_ng()
#
# Added for JLN / LMN:
#   - GARCH(1,1) on idiosyncratic residuals per variable (rugarch)
#   - Aggregation into MU and FU (Option B: all vars in factor step,
#     only macro + financial in uncertainty aggregation)
#   - Proxy SVAR with Cholesky identification (vars)
#
# Variable grouping (Option B):
#   Factor extraction  : ALL variables (macro, financial, external,
#                        non-traditional)
#   MU aggregation     : macro variables only
#   FU aggregation     : financial variables only
#   Non-traditional    : enter factors only; excluded from MU / FU
#   External           : enter factors + assigned to MU
# =============================================================================

### first time install nowcasting and dependencies
#remotes::install_github("cran/matlab")
#remotes::install_github("cran/nowcasting")


##### 0. LIBRARIES #####
suppressMessages({
  library(here)        # portable paths
  library(seasonal)    # seas() – same as Corona et al.
  library(imputeTS)    # na_locf()
  library(nowcasting)  # nowcast(), bai_ng()
  library(vars)        # VAR(), irf(), fevd()
  library(rugarch)     # ugarchspec() / ugarchfit() for GARCH
  library(tidyverse)   # data wrangling + ggplot2
  library(zoo)         # na.locf(), as.yearmon()
  library(tseries)     # adf.test()
  library(patchwork)   # plot composition
  library(lubridate)
  library(urca)
  library(readxl)
  library(mFilter)
})

set.seed(12345)   # mirror Corona et al.


##### 1. USER-BUILT HELPER FUNCTIONS  #####
# These replicate / extend the user-built functions referenced in the
# Corona et al. scripts.
# load functions
source(paste0("../R","/functions.R"))


imef <- read_xlsx("../data/imef_recession_dates.xlsx") |>
  filter(phase == "recession", in_sample %in% c("Yes", "Partial")) |>
  mutate(
    start = as.Date(start),
    end   = if_else(is.na(end) | end == "", Sys.Date(), as.Date(end))
  )

##### 2. PARAMETERS  #####
start          <- "2004/01"          # 
start_lmn      <- "2004/01"          # first date kept for uncertainty output
year_s         <- as.numeric(substr(start, 1, 4))

# Variable groups for Option B 
# Adjust these vectors if the catalog uses different short-names.
MACRO_VARS <- c(
  # Mexican real activity
  "IGAE", "T_EMP_MAN", "ISBSVM", "IAI", "X", "M", "TDU", "REMESAS",
  "ANTAD", "GASOLINAS", "PROD_VEH", "INPC", "CONF_CONS", "CONF_MAN",
  "PEDIDOS_MANU", "CONF_COM", "CONF_SERV", "IMSS", "IMEF",
  "OCUP_HOT",
  # External / US macro (assigned to MU )
  "U_US", "TOTAL_NONFARM", "TOTAL_CONSTR", "TOTAL_MANUF",
  "TOTAL_SERV", "MANUF_USA", "IPI_EUA"
)

FINANCIAL_VARS <- c(
  "IPC", "SP_500", "TC", "TIIE_28", "SPEI", "TARJETAS", "PRECIO", "M4"
)

# 37 vars = 27 + 8 + 2

# Non-traditional → enter factor extraction only, excluded from MU/FU
NONTRADITIONAL_VARS <- c("MOVILIDAD", "SATELITES")


##### 3. Univariate analysis #####
##### 3.0 LOAD RAW DATA AND CATALOG 
cat_dfm     <- read.csv(paste0("../data", "/DFM", "/Catalogo.csv"))
db_dfm_raw  <- read.csv(paste0("../data", "/DFM", "/Datos.csv"),  row.names = 1)

dates_dfm   <- rownames(db_dfm_raw)

# Active variables for IGAE (mirrors Corona et al. line 30-31)
vari        <- "IGAE"

# Target: IGAE in monthly growth-rate form
igae_raw    <- db_dfm_raw[, "IGAE"]
names(igae_raw) <- row.names(db_dfm_raw)
igae_mv     <- cpm(igae_raw, 1)
names(igae_mv) <- dates_dfm[-length(dates_dfm)]
idx_start   <- which(dates_dfm == start)
igae_mv <- igae_mv[idx_start:length(igae_mv)]

igae_yoy <- cpm(igae_raw,12)
names(igae_yoy) <- dates_dfm[-(1:12)]
igae_yoy <- igae_yoy[idx_start:length(igae_yoy)]


# recessions
rect_left  <- as.numeric(format(imef$start, "%Y")) +
  (as.numeric(format(imef$start, "%m")) - 1) / 12
rect_right <- as.numeric(format(imef$end,   "%Y")) +
  (as.numeric(format(imef$end,   "%m")) - 1) / 12


par(mfrow = c(1, 1), mar = c(1, 1, 1, 2))
png(paste('../outputs/eda/','igae_lev.png', sep = ""),
    width = 2145, height = 1155, res = 300)

# 1 blank canvas
igae_ts <- ts(igae_raw[idx_start:length(igae_raw)],
              start = c(2004, 1), frequency = 12)
plot(igae_ts, type = "n",
     xlab = "", ylab = "Index")

# 2. Recession shading — behind everything
rect(xleft   = rect_left,
     xright  = rect_right,
     ybottom = par("usr")[3],
     ytop    = par("usr")[4],
     col     = adjustcolor("grey60", alpha.f = 0.30),
     border  = NA)

# 3. Data on top
lines(igae_ts, col = "blue")
abline(h = 0, lty = 2, col = "grey40")
dev.off()

png(paste('../outputs/eda/','igae_mv.png', sep = ""),
    width = 2145, height = 1155, res = 300)

igae_yoy_ts <- ts(igae_yoy, start = c(2004,1), frequency = 12)
plot(igae_yoy_ts,
     xlab = "", ylab = "% yearly variation")
rect(xleft   = rect_left,
     xright  = rect_right,
     ybottom = par("usr")[3],
     ytop    = par("usr")[4],
     col     = adjustcolor("grey60", alpha.f = 0.30),
     border  = NA)

lines(igae_yoy_ts, col = 'blue')
abline(h = 0, lty = 2)
dev.off()


png(paste('../outputs/eda/','igae_yoy.png', sep = ""),
    width = 2145, height = 1155, res = 300)

igae_mv_ts <- ts(igae_mv, start=c(2004,01), frequency = 12)
plot(igae_mv_ts,
     xlab = 'Time', ylab = '% monthly variation')

rect(xleft   = rect_left,
     xright  = rect_right,
     ybottom = par("usr")[3],
     ytop    = par("usr")[4],
     col     = adjustcolor("grey60", alpha.f = 0.30),
     border  = NA)

lines(igae_mv_ts, col = 'blue')
abline(h = 0, lty =2)

dev.off()

par(mfrow = c(1,1))



##### 3.1 ADF tests for univariate data
# ── 1. SERIES ─────────────────────────────────────────────────────────────────
igae_lev   <- ts(na.omit(igae_raw[idx_start:length(igae_raw)]),
                 start = c(year_s, 1), frequency = 12)
igae_mv_ts <- ts(na.omit(igae_mv),
                 start = c(year_s, 1), frequency = 12)
igae_yoy_ts <- ts(na.omit(igae_yoy_ts),
                  start = c(year_s,1), frequency = 12)
# ── 2. ADF  (H0: unit root) ───────────────────────────────────────────────────
# levels:    include trend + drift  (series has upward trend)
# variation: include drift only     (mean-zero after differencing)
adf_lev <- ur.df(igae_lev,   type = "trend", selectlags = "AIC")
adf_mv  <- ur.df(igae_mv_ts, type = "drift", selectlags = "AIC")
adf_yoy <- ur.df(igae_yoy_ts, type = 'drift', selectlags = 'AIC')
cat("── ADF: IGAE levels ──────────────────────────────\n"); summary(adf_lev)
cat("── ADF: IGAE monthly variation ───────────────────\n"); summary(adf_mv)
cat("── ADF: IGAE YoY variation ───────────────────\n"); summary(adf_yoy)


# ── 3. PP  (H0: unit root, non-parametric correction) ────────────────────────
pp_lev <- ur.pp(igae_lev,   type = "Z-tau", model = "trend")
pp_mv  <- ur.pp(igae_mv_ts, type = "Z-tau", model = "constant")
pp_yoy <- ur.pp(igae_yoy_ts, type = "Z-tau", model = 'constant')
cat("── PP: IGAE levels ───────────────────────────────\n"); summary(pp_lev)
cat("── PP: IGAE monthly variation ────────────────────\n"); summary(pp_mv)
cat("── PP: IGAE YoY variation ────────────────────\n"); summary(pp_yoy)

# ── 4. KPSS  (H0: stationary — reverse null) ─────────────────────────────────
kpss_lev <- ur.kpss(igae_lev,   type = "tau")   # trend-stationary H0
kpss_mv  <- ur.kpss(igae_mv_ts, type = "mu")    # level-stationary H0
kpss_yoy  <- ur.kpss(igae_yoy_ts, type = "mu")    # level-stationary H0

cat("── KPSS: IGAE levels ─────────────────────────────\n"); summary(kpss_lev)
cat("── KPSS: IGAE monthly variation ──────────────────\n"); summary(kpss_mv)
cat("── KPSS: IGAE YoY variation ──────────────────\n"); summary(kpss_yoy)

# ── 5. SUMMARY TABLE ──────────────────────────────────────────────────────────
extract_adf  <- function(u) c(stat = u@teststat[1], cv1 = u@cval[1,1],
                              cv5 = u@cval[1,2], cv10 = u@cval[1,3])
extract_pp   <- function(u) c(stat = u@teststat,    cv1 = u@cval[1],
                              cv5 = u@cval[2],   cv10 = u@cval[3])
extract_kpss <- function(u) c(stat = u@teststat,    cv1 = u@cval[4],
                              cv5 = u@cval[3],   cv10 = u@cval[2])

results <- rbind(
  data.frame(series = "Levels",           test = "ADF",  t(extract_adf(adf_lev))),
  data.frame(series = "Monthly variation",test = "ADF",  t(extract_adf(adf_mv))),
  data.frame(series = "Yearly variation",test = "ADF",  t(extract_adf(adf_yoy))),
  data.frame(series = "Levels",           test = "PP",   t(extract_pp(pp_lev))),
  data.frame(series = "Monthly variation",test = "PP",   t(extract_pp(pp_mv))),
  data.frame(series = "Yearly variation",test = "PP",   t(extract_pp(pp_yoy))),
  data.frame(series = "Levels",           test = "KPSS", t(extract_kpss(kpss_lev))),
  data.frame(series = "Monthly variation",test = "KPSS", t(extract_kpss(kpss_mv))),
  data.frame(series = "Yearly variation",test = "KPSS", t(extract_kpss(kpss_yoy)))
)

cat("\n── Unit root tests — IGAE ──────────────────────────────────\n")
cat("   ADF/PP  H0: unit root    → reject if stat < cv (left tail)\n")
cat("   KPSS    H0: stationary   → reject if stat > cv (right tail)\n\n")
results_print <- results
results_print[, sapply(results_print, is.numeric)] <- round(
  results_print[, sapply(results_print, is.numeric)], 3
)
print(results_print, row.names = FALSE)

##### 3.2 Structural breaks 

library(strucchange)   # Chow, CUSUM, Bai-Perron
# urca already loaded  # ur.za (Zivot-Andrews)

# ── helper: date → observation index ─────────────────────────────────────────
obs_idx <- function(ts_obj, yr, mo)
  which(round(time(ts_obj), 6) == round(yr + (mo - 1) / 12, 6))

# ── 1. CHOW TEST at known Mexican recession dates ─────────────────────────────
# H0: no structural break at the specified point
chow_mv_2009  <- sctest(igae_mv_ts ~ 1, type = "Chow",
                        point = obs_idx(igae_mv_ts, 2009, 5))
chow_mv_2020  <- sctest(igae_mv_ts ~ 1, type = "Chow",
                        point = obs_idx(igae_mv_ts, 2020, 5))
chow_lev_2009 <- sctest(igae_lev   ~ 1, type = "Chow",
                        point = obs_idx(igae_lev, 2009, 5))
chow_lev_2020 <- sctest(igae_lev   ~ 1, type = "Chow",
                        point = obs_idx(igae_lev, 2020, 5))

chow_res <- data.frame(
  Series = c("Levels", "Levels", "Monthly var.", "Monthly var."),
  Break  = c("2009/05", "2020/05", "2009/05", "2020/05"),
  F_stat = round(c(chow_lev_2009$statistic,  chow_lev_2020$statistic,
                   chow_mv_2009$statistic,   chow_mv_2020$statistic), 3),
  p_value = round(c(chow_lev_2009$p.value,   chow_lev_2020$p.value,
                    chow_mv_2009$p.value,    chow_mv_2020$p.value), 4)
)
cat("\n── Chow tests ──\n"); print(chow_res, row.names = FALSE)


# # ── 2. CUSUM TEST (OLS-CUSUM — global parameter stability) ───────────────────
# cusum_lev <- efp(igae_lev   ~ 1, type = "OLS-CUSUM")
# cusum_mv  <- efp(igae_mv_ts ~ 1, type = "OLS-CUSUM")
# 
# cat("\n── CUSUM: levels ──\n");       print(sctest(cusum_lev))
# cat("── CUSUM: monthly variation ──\n"); print(sctest(cusum_mv))
# 
# par(mfrow = c(1, 2))
# plot(cusum_lev, main = "OLS-CUSUM: IGAE levels",
#      ylab = "Empirical fluctuation process")
# abline(v = c(2009 + 4/12, 2020 + 4/12), col = "red", lty = 2)
# 
# plot(cusum_mv, main = "OLS-CUSUM: IGAE monthly variation",
#      ylab = "")
# abline(v = c(2009 + 4/12, 2020 + 4/12), col = "red", lty = 2)
# par(mfrow = c(1, 1))

##### 3.3 BAI-PERRON: data-driven break detection ────────────────────────────────
# h = 0.10 → minimum segment length = 10% of sample (~24 months)
bp_lev <- breakpoints(igae_lev   ~ 1, h = 0.10)
bp_mv  <- breakpoints(igae_mv_ts ~ 1, h = 0.10)

cat("\n── Bai-Perron: levels ──\n");           print(summary(bp_lev))
cat("── Bai-Perron: monthly variation ──\n"); print(summary(bp_mv))

# Detected break dates
cat("\nDetected break dates — levels:\n")
print(time(igae_lev)[bp_lev$breakpoints])
cat("Detected break dates — monthly variation:\n")
print(time(igae_mv_ts)[bp_mv$breakpoints])


##### 3.4 ZIVOT-ANDREWS: unit root allowing for one endogenous break ─────────────
# Complements ADF — if ZA rejects, series is I(0) around a broken trend/mean
za_lev <- ur.za(igae_lev,   model = "both")        # break in trend + intercept
za_mv  <- ur.za(igae_mv_ts, model = "intercept")   # break in mean only

cat("\n── Zivot-Andrews: levels ──\n");           summary(za_lev)
cat("── Zivot-Andrews: monthly variation ──\n"); summary(za_mv)


### 3.5 Permanent + Transitory decomposition


# ── 3.3 Permanent/Transitory Decomposition — Beveridge-Nelson and HP ──────────


# Log-level IGAE (reuse igae_lev from section 3.1)
igae_log_ts <- log(igae_lev)

# ── Step 1: ARIMA model selection ─────────────────────────────────────────────
fit210 <- arima(igae_log_ts, order = c(2,1,0))
fit110 <- arima(igae_log_ts, order = c(1,1,0))
fit212 <- arima(igae_log_ts, order = c(2,1,2))

ic_table <- data.frame(
  Model = c("ARIMA(2,1,0)", "ARIMA(1,1,0)", "ARIMA(2,1,2)"),
  N = c(fit210$nobs,fit110$nobs,fit212$nobs),
  Loglik = c(fit210$loglik,fit110$loglik,fit212$loglik),
  AIC   = c(AIC(fit210), AIC(fit110), AIC(fit212)),
  BIC   = c(BIC(fit210), BIC(fit110), BIC(fit212))
  )
cat("\n── ARIMA model selection ──────────────────────────────────────\n")
print(ic_table, row.names = FALSE)

# ── Step 2: Extract ARIMA(2,1,0) parameters ───────────────────────────────────
coefs  <- coef(fit210)
phi1   <- coefs["ar1"]
phi2   <- coefs["ar2"]
theta1 <- 0
theta2 <- 0
mu     <- ifelse("intercept" %in% names(coefs), coefs["intercept"], 0)
denom  <- 1 - phi1 - phi2

cat(sprintf("\nARIMA(2,1,0): phi1=%.4f  phi2=%.4f  theta1=%.4f  theta2=%.4f  mu=%.6f\n",
            phi1, phi2, theta1, theta2, mu))

# Ljung-Box white-noise test on residuals (equivalent to Stata wntestq)
lb <- Box.test(residuals(fit212), lag = 12, type = "Ljung-Box")
cat(sprintf("Ljung-Box Q(12): stat=%.3f  p=%.4f\n", lb$statistic, lb$p.value))

# fail to reject white noise by a wide margin

# ── Step 3: BN transitory component ───────────────────────────────────────────
# BN cycle for ARIMA(2,1,0):
#   C_t = [(φ1+φ2)/D]*(Δy_t−μ) + [φ2/D]*(Δy_{t-1}−μ)
#        + [(θ1+θ2)/D]*ε_t      + [θ2/D]*ε_{t-1}
# where D = 1 − φ1 − φ2

n      <- length(igae_log_ts)
resids <- as.numeric(residuals(fit210))

# Pad diff to length n so index t maps to Δy_t = y_t − y_{t-1}
dlgy   <- c(NA, as.numeric(diff(igae_log_ts)))

C_bn   <- numeric(n)   # initialised to zero

for (t in 3:n) {
  if (anyNA(c(resids[t], resids[t-1], dlgy[t], dlgy[t-1]))) next
  C_bn[t] <- ((phi1 + phi2) / denom) * (dlgy[t]   - mu) +
    (phi2            / denom) * (dlgy[t-1] - mu) +
    ((theta1 + theta2) / denom) * resids[t]      +
    (theta2           / denom) * resids[t-1]
}

trend_bn_ts <- ts(as.numeric(igae_log_ts) - C_bn,
                  start = start(igae_log_ts), frequency = 12)
C_bn_ts     <- ts(C_bn, start = start(igae_log_ts), frequency = 12)

# ── Step 4: HP filter (λ = 14400 for monthly) ─────────────
hp_res   <- hpfilter(igae_log_ts, freq = 14400)
trend_hp <- hp_res$trend
C_hp     <- hp_res$cycle

# ── Step 5: Plot — cyclical components ────────────────────────────────────────
png("../outputs/eda/igae_cycle.png", width = 1950, height = 1050, res = 300)
par(mar = c(4, 4, 3, 1))
ylim_c <- range(c(C_bn, C_hp), na.rm = TRUE)
plot(C_hp, type = "n", ylim = ylim_c,
     xlab = "Time", ylab = "Log deviation from trend",
     main = "Cyclical Components of IGAE (log)")
rect(xleft = rect_left, xright = rect_right,
     ybottom = par("usr")[3], ytop = par("usr")[4],
     col = adjustcolor("grey60", alpha.f = 0.25), border = NA)
abline(h = 0, col = "grey40", lwd = 0.8)
lines(C_hp,     col = "steelblue", lwd = 1.5)
lines(C_bn_ts,  col = "firebrick", lwd = 1.5, lty = 2)
legend("bottomleft", legend = c("HP cycle", "BN cycle"),
       col = c("steelblue", "firebrick"), lty = c(1,2), lwd = 1.5, bty = "n")
dev.off()

# ── Step 6: Plot — permanent components ───────────────────────────────────────
png("../outputs/eda/igae_trend.png", width = 1950, height = 1050, res = 300)
par(mar = c(4, 4, 3, 1))
plot(igae_log_ts, col = "grey40", lwd = 1,
     xlab = "Time", ylab = "Log IGAE",
     main = "Permanent Components of IGAE (log)")
rect(xleft = rect_left, xright = rect_right,
     ybottom = par("usr")[3], ytop = par("usr")[4],
     col = adjustcolor("grey60", alpha.f = 0.25), border = NA)
lines(trend_hp,    col = "steelblue", lwd = 1.8)
lines(trend_bn_ts, col = "firebrick", lwd = 1.8, lty = 2)
legend("topleft",
       legend = c("Actual ln(IGAE)", "HP trend", "BN trend"),
       col    = c("grey40", "steelblue", "firebrick"),
       lty    = c(1, 1, 2), lwd = c(1, 1.8, 1.8), bty = "n")
dev.off()



##### 4. Multivariate analysis #####

##### continue with analysis
act         <- as.character(cat_dfm[cat_dfm[, vari] == 1, "Short"])
db_dfm_raw  <- db_dfm_raw[, act, drop = FALSE]

# Keep from 'start' onward
idx_start   <- which(dates_dfm == start)
db_dfm_raw  <- db_dfm_raw[idx_start:nrow(db_dfm_raw), , drop = FALSE]
dates_dfm   <- rownames(db_dfm_raw)

# Forward-fill NAs (same as na_locf in Corona et al.)
db_dfm_raw  <- na_locf(db_dfm_raw)
# Backward-fill remaining leading NAs (variables that start after 'start',
# e.g. GASOLINAS ~2014, MOVILIDAD/SATELITES ~2020).
db_dfm_raw  <- as.data.frame(
  zoo::na.locf(as.matrix(db_dfm_raw), fromLast = TRUE, na.rm = FALSE)
)

cat(sprintf("Data loaded: %d observations × %d variables (%s – %s)\n",
            nrow(db_dfm_raw), ncol(db_dfm_raw),
            head(dates_dfm, 1), tail(dates_dfm, 1)))


##### 4.1 SEASONAL ADJUSTMENT  ───────────────────────────────────────────────────
# Mirrors 01_preparing.R lines 88-105 exactly.

sa_vars <- as.character(cat_dfm[cat_dfm[, "SA"] == 1, "Short"])
sa_vars <- intersect(sa_vars, colnames(db_dfm_raw))

for (v in sa_vars) {
  date_v  <- rownames(db_dfm_raw)[!is.na(db_dfm_raw[, v])]
  ts_v    <- ts(db_dfm_raw[, v],
                start     = c(as.integer(substr(date_v[1], 1, 4)),
                              as.integer(substr(date_v[1], 6, 7))),
                frequency = 12)
  sa_v    <- tryCatch(seas(ts_v)$series$s11, error = function(e) NULL)
  if (!is.null(sa_v)) {
    db_dfm_raw[!is.na(ts_v), v] <- as.numeric(sa_v)
    db_dfm_raw[db_dfm_raw[, v] < 0 & !is.na(db_dfm_raw[, v]), v] <- 0.01
  }
}

cat(sprintf("Seasonal adjustment applied to %d variables.\n", length(sa_vars)))


##### 4.5. OPTIMAL STATIONARITY TRANSFORMATIONS  ──────────────────────────────────
# Mirrors 01_preparing.R lines 107-154.
# For each variable, best_trans() chooses the transformation that maximises
# its correlation with monthly IGAE growth (d = 1).

db_dfm_trans  <- matrix(NA_real_, nrow(db_dfm_raw), ncol(db_dfm_raw),
                        dimnames = dimnames(db_dfm_raw))

type_trans <- data.frame(
  Transformation = character(ncol(db_dfm_raw)),
  Correlation    = numeric(ncol(db_dfm_raw)),
  row.names      = colnames(db_dfm_raw),
  stringsAsFactors = FALSE
)

for (i in seq_len(ncol(db_dfm_raw))) {
  var_i <- colnames(db_dfm_raw)[i]
  
  rel_i <- as.character(
    cat_dfm[cat_dfm[, "Short"] == var_i, "Rel", drop = TRUE]
  )
  if (length(rel_i) != 1 || is.na(rel_i) || !rel_i %in% c("I", "P", "N"))
    rel_i <- "I"
  
  td <- best_trans(
    y       = igae_mv,
    x       = db_dfm_raw[, var_i],
    names_y = names(igae_mv),        # dates_dfm[-1]  — now has names
    names_x = dates_dfm,             # full date vector
    rel     = rel_i
  )
  
  db_dfm_trans[, var_i]               <- td$x_trans
  type_trans[var_i, "Transformation"] <- td$trans
  type_trans[var_i, "Correlation"]    <- td$rho[td$trans]
}

# Drop the first 12 rows (NAs from annual transformation)
db_dfm_trans <- db_dfm_trans[13:nrow(db_dfm_trans), , drop = FALSE]
dates_trans  <- rownames(db_dfm_trans)

cat("\nTransformation summary:\n")
print(table(type_trans$Transformation))


# plot igae_mv against each of the variables 
dates_trans  <- rownames(db_dfm_trans)

# start date of the analysis window (after 12-row annual-transformation drop)
yr_s <- as.integer(substr(dates_trans[1], 1, 4))   # 2005
mo_s <- as.integer(substr(dates_trans[1], 6, 7))   # 1

igae_aligned <- igae_raw[names(igae_raw) %in% dates_trans]
igae_ts <- ts(scale(igae_aligned), start = c(yr_s, mo_s), frequency = 12)

# ── helper: plot one group ────────────────────────────────────────────────────
igae_ts <- ts(scale(igae_raw), start = c(2004, 1), frequency = 12)

var_labels <- setNames(cat_dfm$Description, cat_dfm$Short)

plot_group <- function(group_vars, group_label) {
  
  idx <- which(colnames(db_dfm_trans) %in% group_vars)
  n   <- length(idx)
  nc  <- ceiling(sqrt(n))
  nr  <- ceiling(n / nc)
  
  par(mfrow = c(nr, nc), mar = c(2, 2, 1.5, 0.5), oma = c(0, 0, 2, 0))
  
  for (i in idx) {
    yr_s <- as.integer(substr(dates_trans[1], 1, 4))
    mo_s <- as.integer(substr(dates_trans[1], 6, 7))
    var_ts <- ts(db_dfm_trans[, i], start = c(yr_s, mo_s), frequency = 12)
    trans  <- type_trans$Transformation[i]          # positional — same order as db_dfm_trans
    
    short_name <- colnames(db_dfm_trans)[i]
    eng_label  <- var_labels[short_name] 
    
    if (is.na(eng_label)) eng_label <- short_name
    
    plot(var_ts,
         xlab = "", ylab = "",
         main = paste(eng_label, "-", trans),
         col  = "blue", cex.main = 0.8)
    lines(igae_ts * sd(db_dfm_trans[, i], na.rm = TRUE),
          col = "red", lty = 2, lwd = 0.8)
    abline(h = 0, col = "grey50")
  }
  
}

macro_vars_present <- intersect(MACRO_VARS, colnames(db_dfm_trans))
png(paste('../outputs/eda/','macro_dfm_vars.png', sep = ""),
    width = 2145, height = 2145, res = 300)
plot_group(macro_vars_present, "Macro variables (blue) vs. IGAE growth (red dashed)")
dev.off()

fin_vars_present <- intersect(FINANCIAL_VARS, colnames(db_dfm_trans))
png(paste('../outputs/eda/','fin_dfm_vars.png', sep = ""),
    width = 2145, height = 1155, res = 300)
plot_group(fin_vars_present, "Financial variables (blue) vs. IGAE growth (red dashed)")
dev.off()

##### 5. FACTOR EXTRACTION #####
##### 5.1 FACTOR COUNT SELECTION  ────────────────────────────────────────────────
# Use bai_ng() from the nowcasting package — same criterion as Corona et al.
# We also add Onatski as a robustness check.

# Scale the panel before factor extraction (same as scale(Y) in Corona et al.)
X_scaled    <- scale(db_dfm_trans)

# Bai-Ng (2002) IC criterion
bn_result   <- bai_ng(X_scaled, demean = 2)
rhat        <- bn_result$rhat[1]   # IC_p1 — conservative, use IC_p2 for robustness
cat(sprintf("\nBai-Ng factor selection: r = %d (IC_p1)\n", rhat))
cat(sprintf("                         r = %d (IC_p2)\n", bn_result$rhat[2]))

# Optional Onatski check (uncomment if nowcasting version supports it)
# onat <- tryCatch(onatski2010(X_scaled), error = function(e) NULL)
# if (!is.null(onat)) cat(sprintf("Onatski (2010) ED:        r = %d\n", onat["ed"]))


##### 5.2. DFM EXTRACTION  ────────────────────────────────────────────────────────
# Mirrors Corona et al. 02_models.R lines 83-119.
# We use nowcast() with method = "2s" on the FULL panel (all variable groups
# per Option B), then extract:
#   fhat_    = common factors  (T × r)
#   Phat_    = loadings        (N × r)
#   e_idio   = idiosyncratic component per variable (T × N)


# ── Build db_now: YFCST aligned to dates_trans with one trailing NA ──────────
# igae_mv[13:] skips the 12 periods consumed by annual transforms, matching
# X_scaled row 1. Appending NA gives the single trailing NA nowcast() needs.
yfcst_in <- c(igae_mv[13:length(igae_mv)], NA)

db_now   <- ts(
  cbind(yfcst_in, X_scaled),
  start     = c(year_s + 1, 1),
  frequency = 12
)
colnames(db_now) <- c("YFCST", colnames(X_scaled))

freq_now <- rep(12, ncol(db_now))

# ── Verify before calling nowcast ────────────────────────────────────
cat("db_now dim:      ", dim(db_now), "\n")
cat("YFCST NAs:       ", sum(is.na(db_now[, "YFCST"])), "\n")   # should be 1
cat("freq_now length: ", length(freq_now), "\n")                 # must = ncol(db_now)
cat("rhat:            ", rhat, "\n")

# ── Call nowcast ──────────────────────────────────────────────────────
nowcast_obj <- nowcast(
  YFCST ~.,
  data      = db_now,
  r         = rhat,
  p         = 1,
  q         = rhat,
  method    = "2s",
  frequency = freq_now
)


# Loadings (N × r)
Phat_ <- nowcast_obj$factors$eigen$vectors[, 1:rhat, drop = FALSE]
rownames(Phat_) <- colnames(db_now)[-1]

# Common factors (T × r)
fhat_ <- as.matrix(nowcast_obj$factors$dynamic_factors)

# Sign normalisation — align factor with IGAE growth (mirrors lines 113-116)
igae_ts <- db_now[, "YFCST"]
if (cor(igae_ts, fhat_[, 1], use = "pairwise.complete.obs") < 0) {
  Phat_ <- -Phat_
  fhat_ <- -fhat_
}

# Common component and idiosyncratic residuals (mirrors lines 118-119)
T_obs    <- nrow(X_scaled)
X_trim   <- X_scaled[1:nrow(fhat_), , drop = FALSE]
fhat_tr  <- fhat_[1:nrow(X_trim), , drop = FALSE]
FP       <- fhat_tr %*% t(Phat_)          # T × N  common component
e_idio   <- X_trim - FP                   # T × N  idiosyncratic component

dates_fhat <- dates_trans[1:nrow(fhat_tr)]

cat(sprintf("DFM extracted: %d factors, %d observations, %d variables\n",
            rhat, nrow(fhat_tr), ncol(e_idio)))

##### 5.3 DFM FACTOR DIAGNOSTICS ────────────────────────────────────────────────

### percentage of variation attributed to idiosyncratic error
# ── Variance decomposition: common vs idiosyncratic ───────────────────────────

r2_by_var <- 1 - apply(e_idio,  2, var, na.rm = TRUE) /
  apply(X_trim,  2, var, na.rm = TRUE)

idio_share <- data.frame(
  variable      = colnames(e_idio),
  r2_common     = round(r2_by_var, 3),          # share explained by factors
  idio_share    = round(1 - r2_by_var, 3)       # share in idiosyncratic
) |> arrange(desc(idio_share))

cat("\n── Variance decomposition (common factors vs idiosyncratic) ──\n")
print(idio_share, row.names = FALSE)
cat(sprintf("\nMedian idiosyncratic share across variables: %.1f%%\n",
            median(idio_share$idio_share) * 100))

# Three figures:
#   fig7_factors_ts.png   — factor time series (small multiples)
#   fig8_loadings_heatmap.png — variable-factor loadings heat map
#   fig9_variance_contrib.png — cumulative variance explained by factor

# Build tidy factor data frame
fhat_df <- as.data.frame(fhat_tr)
colnames(fhat_df) <- paste0("F", seq_len(ncol(fhat_tr)))
fhat_df$date <- as.Date(paste0(dates_fhat, "/01"), format = "%Y/%m/%d")

# ── Fig 7: Factor time series ─────────────────────────────────────────────────
fhat_long <- fhat_df |>
  pivot_longer(-date, names_to = "factor", values_to = "value") |>
  mutate(factor = factor(factor, levels = paste0("F", seq_len(rhat))))

p_factors <- ggplot(fhat_long, aes(x = date, y = value)) +
  geom_line(colour = "#1f4e79", linewidth = 0.55) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  facet_wrap(~ factor, ncol = 4, scales = "free_y") +
  labs(
    title    = "DFM Common Factors — Time Series",
    subtitle = sprintf("%d factors", rhat),
    x        = NULL,
    y        = "Factor value (standardised units)"
  ) +
  theme_bw(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(size = 7),   # panel labels
    axis.text        = element_text(size = 6)    # tick labels
        )

ggsave(paste0("../outputs", "/lmn", "/fig7_factors_ts.png"),
       p_factors, width = 9.5, height = 5, dpi = 200)
cat("\nFig 7 saved: fig7_factors_ts.png\n")

##### 5.4 DFM FACTOR INTERPRETATION 
##### Loadings heatmap ───────────────────────────────────────────────────
loadings_df <- as.data.frame(Phat_)
colnames(loadings_df) <- paste0("F", seq_len(rhat))
loadings_df$variable <- rownames(Phat_)

# Order variables by group then idiosyncratic share (for readability)
var_order <- idio_share$variable[idio_share$variable %in% rownames(Phat_)]
if (length(var_order) == 0) var_order <- rownames(Phat_)

loadings_long <- loadings_df |>
  pivot_longer(-variable, names_to = "factor", values_to = "loading") |>
  mutate(
    # map short name → English description; fall back to short name if missing
    variable = var_labels[variable] |> (\(x) ifelse(is.na(x), variable, x))(),
    variable = factor(variable, levels = rev(var_labels[var_order] |>
                                               (\(x) ifelse(is.na(x), var_order, x))())),
    factor   = factor(factor, levels = paste0("F", seq_len(rhat)))
  )

p_heat <- ggplot(loadings_long, aes(x = factor, y = variable, fill = loading)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = round(loading, 2)), size = 2.2, colour = "black") +
  scale_fill_gradient2(low = "#c55a11", mid = "white", high = "#1f4e79",
                       midpoint = 0, name = "Loading") +
  labs(
    title    = "DFM Factor Loadings",
    subtitle = "Variables ordered by idiosyncratic share (high → low)",
    x        = "Factor",
    y        = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_text(size = 7),
        panel.grid  = element_blank())

ggsave(paste0("../outputs", "/lmn", "/fig8_loadings_heatmap.png"),
       p_heat, width = 5.9, height = 5, dpi = 200)
cat("Fig 8 saved: fig8_loadings_heatmap.png\n")

##### 5.5 Variance contribution per factor ───────────────────────────────────
# Share of total panel variance explained by each factor individually
# Var explained by F_k = var(fhat_tr[,k] * Phat_[,k]') summed over variables
total_var <- sum(apply(X_trim, 2, var, na.rm = TRUE))
var_by_factor <- sapply(seq_len(rhat), function(k) {
  comp_k <- outer(fhat_tr[, k], Phat_[, k])   # T × N rank-1 component
  sum(apply(comp_k, 2, var, na.rm = TRUE))
})
pct_by_factor <- var_by_factor / total_var * 100

scree_df <- data.frame(
  factor    = paste0("F", seq_len(rhat)),
  pct       = pct_by_factor,
  cumulative = cumsum(pct_by_factor)
) |> mutate(factor = factor(factor, levels = paste0("F", seq_len(rhat))))

p_scree <- ggplot(scree_df, aes(x = factor)) +
  geom_col(aes(y = pct), fill = "#1f4e79", alpha = 0.75, width = 0.6) +
  geom_line(aes(y = cumulative, group = 1), colour = "#c55a11",
            linewidth = 0.8) +
  geom_point(aes(y = cumulative), colour = "#c55a11", size = 2.5) +
  scale_y_continuous(
    name     = "Individual contribution (%)",
    sec.axis = sec_axis(~ ., name = "Cumulative (%)")
  ) +
  labs(
    title    = "DFM: Variance Contribution by Factor",
    subtitle = "Bars = individual; Orange line = cumulative",
    x        = "Factor"
  ) +
  theme_bw(base_size = 10)

ggsave(paste0("../outputs", "/lmn", "/fig9_variance_contrib.png"),
       p_scree, width = 7, height = 4.5, dpi = 150)
cat("Fig 9 saved: fig9_variance_contrib.png\n")

##### 5.5 Top-5 loadings table per factor ──────────────────────────────────────────
cat("\n── Top-5 loadings per factor (absolute value) ──\n")
for (k in seq_len(rhat)) {
  top5 <- sort(abs(Phat_[, k]), decreasing = TRUE)[1:5]
  top5_names <- names(top5)
  top5_vals  <- round(Phat_[top5_names, k], 3)
  cat(sprintf("F%d: %s\n", k,
              paste(sprintf("%s(%.3f)", top5_names, top5_vals),
                    collapse = ", ")))
}

# Print scree table for notes
cat("\n── Variance contribution by factor ──\n")
print(scree_df |> mutate(across(where(is.numeric), \(x) round(x, 2))),
      row.names = FALSE)


##### 6. UNCERTAINTY COMPUTATION — JLN STEP #####
# For each variable i:
#   1. Idiosyncratic residual e_it = x_it - F_t' * lambda_i  (already in e_idio)
#   2. GARCH(1,1) on e_it  →  conditional std dev  sigma_it
# This is the JLN "genuine unpredictability" measure.

cat("\nFitting GARCH(1,1) on idiosyncratic components...\n")

sigma_mat <- matrix(NA_real_, nrow(e_idio), ncol(e_idio),
                    dimnames = dimnames(e_idio))

for (j in seq_len(ncol(e_idio))) {
  if (j %% 10 == 0)
    cat(sprintf("  Variable %d / %d\r", j, ncol(e_idio)))
  sigma_mat[, j] <- fit_garch11(e_idio[, j])
}
cat(sprintf("  Variable %d / %d — done.\n", ncol(e_idio), ncol(e_idio)))

sigma_df <- as.data.frame(sigma_mat)
rownames(sigma_df) <- dates_fhat


##### 6.1 AGGREGATE MU AND FU  ───────────────────────────────────────────────────
# Option B: average conditional vols only over their respective groups.
# Non-traditional variables have sigma estimated but are NOT averaged in.

macro_in_data  <- intersect(MACRO_VARS,        colnames(sigma_df))
fin_in_data    <- intersect(FINANCIAL_VARS,    colnames(sigma_df))
nontr_in_data  <- intersect(NONTRADITIONAL_VARS, colnames(sigma_df))

MU_raw <- rowMeans(sigma_df[, macro_in_data, drop = FALSE], na.rm = TRUE)
FU_raw <- rowMeans(sigma_df[, fin_in_data,   drop = FALSE], na.rm = TRUE)

# Standardise for comparability (zero mean, unit variance)
MU <- (MU_raw - mean(MU_raw, na.rm = TRUE)) / sd(MU_raw, na.rm = TRUE)
FU <- (FU_raw - mean(FU_raw, na.rm = TRUE)) / sd(FU_raw, na.rm = TRUE)

uncertainty_indices <- tibble(
  date = as.Date(as.yearmon(dates_fhat, "%Y/%m")),
  MU   = MU,
  FU   = FU,
  MU_raw = MU_raw,
  FU_raw = FU_raw
)

cat(sprintf("\nMU: mean = %.3f | sd = %.3f | n_series = %d\n",
            mean(MU, na.rm = TRUE), sd(MU, na.rm = TRUE), length(macro_in_data)))
cat(sprintf("FU: mean = %.3f | sd = %.3f | n_series = %d\n",
            mean(FU, na.rm = TRUE), sd(FU, na.rm = TRUE), length(fin_in_data)))
cat(sprintf("MU–FU correlation: %.3f\n", cor(MU, FU, use = "complete.obs")))


# 6.4 Plot MU/FU with IMEF business cycle dating 


plot_df <- tibble(
  date = as.Date(as.yearmon(dates_fhat, "%Y/%m")),   # adjust to your date vector
  MU   = MU,
  FU   = FU,
  
) |> pivot_longer(c(MU, FU), names_to = "index", values_to = "value")

mu_fu_imef_plot <- ggplot(plot_df, aes(x = date, y = value, colour = index)) +
  geom_rect(
    data = imef,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "grey80", alpha = 0.4
  ) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  scale_colour_manual(values = c(MU = "#1f4e79", FU = "#c55a11")) +
  labs(
    title    = "Macro and Financial Uncertainty Indices — Mexico",
    subtitle = "Shaded: IMEF recession dates",
    x = NULL, y = "Standardised index", colour = NULL
  ) +
  theme_bw(base_size = 10)


##### 7. PROXY SVAR — LMN STEP  #####
# Identification: Cholesky ordering
#   Ordering A (LMN baseline): MU first  → macro uncertainty is the
#                               primitive shock; FU responds to MU
#   Ordering B (alternative):  FU first  → financial uncertainty is
#                               the primitive shock
#
# The divergence in the IGAE impulse responses across orderings is the
# core LMN diagnostic.

# 7.0 Align IGAE growth with uncertainty dates
igae_growth_svar <- igae_mv[dates_fhat]      # names already correct from line 99

# 7.1 VAR system: [MU_t, FU_t, igae_growth_t]
svar_df <- tibble(
  date        = uncertainty_indices$date,
  MU          = MU,
  FU          = FU,
  igae_growth = as.numeric(igae_growth_svar)
) |> drop_na()

# 7.2 time series format
var_data_ts <- ts(svar_df[, c("MU", "igae_growth", "FU")],
                  start = c(year(min(svar_df$date)),
                            month(min(svar_df$date))),
                  frequency = 12)

# 7.3 structural break 
ts_start <- start(var_data_ts)          # e.g. c(2005, 1)
n        <- nrow(as.data.frame(var_data_ts))

# Row index for any (year, month) within the ts
row_of <- function(yr, mo) (yr - ts_start[1]) * 12 + (mo - ts_start[2]) + 1

D_2009 <- integer(n);  D_2009[row_of(2009, 5)] <- 1
D_2020 <- integer(n);  D_2020[row_of(2020, 4)] <- 1

exog_mat <- cbind(D_2009, D_2020)

# # VAR with break dummies
# var_fit_breaks <- VAR(
#   var_data_ts[, c("MU", "igae_growth", "FU")],
#   p      = p_opt,
#   type   = "const",
#   exogen = exog_mat
# )

# 7.4 lag selection with dummies return 
lag_sel  <- VARselect(var_data_ts, lag.max = 6, type = "const",
                      exogen = exog_mat)
p_aic    <- lag_sel$selection["AIC(n)"]
p_sc     <- lag_sel$selection["SC(n)"]
p_opt    <- max(1, min(p_sc, 4))    # SC, capped at 4

cat(sprintf("\nVAR lag order: AIC = %d, SC = %d → using p = %d\n",
            p_aic, p_sc, p_opt))


# 7.5 LIKELIHOOD RATIO TEST for lag order ───────────────────────────────────────
# H0: p = p_r  (restricted)  vs  H1: p = p_u  (unrestricted)
# LR = T * (log|Σ_r| - log|Σ_u|) ~ χ²(K² × (p_u - p_r))
# Both models must share the same estimation sample → trim to p_max lags

K     <- ncol(var_data_ts)   # 3 variables
p_max <- 4                   # ceiling for comparison

# Re-estimate all models on the same sample (rows p_max+1 : T)
common_data <- window(var_data_ts,
                      start = time(var_data_ts)[p_max + 1])

exog_common <- tail(exog_mat, nrow(as.data.frame(common_data)))

var_p1 <- VAR(common_data[, c("MU", "igae_growth", "FU")], p = 1, 
              type = "const", exogen = exog_common)
var_p2 <- VAR(common_data[, c("MU", "igae_growth", "FU")], p = 2, 
              type = "const", exogen = exog_common)
var_p3 <- VAR(common_data[, c("MU", "igae_growth", "FU")], p = 3, 
              type = "const", exogen = exog_common)
var_p4 <- VAR(common_data[, c("MU", "igae_growth", "FU")], p = 4, 
              type = "const", exogen = exog_common)
var_p5 <- VAR(common_data[, c("MU", "igae_growth", "FU")], p = 5, 
              type = "const", exogen = exog_common)


lr_test <- function(var_r, var_u) {
  e_r <- residuals(var_r)
  e_u <- residuals(var_u)
  n   <- nrow(e_u)                           # common sample size
  e_r <- tail(e_r, n)                        # align: drop extra rows from restricted
  Sigma_r <- crossprod(e_r) / n
  Sigma_u <- crossprod(e_u) / n
  lr  <- n * (log(det(Sigma_r)) - log(det(Sigma_u)))
  df  <- K^2 * (var_u$p - var_r$p)
  data.frame(
    H0       = paste0("p = ", var_r$p),
    H1       = paste0("p = ", var_u$p),
    LR_stat  = round(lr, 3),
    df       = df,
    p_value  = round(1 - pchisq(lr, df), 4)
  )
}

lr_results <- rbind(
  lr_test(var_p1, var_p2),   # sequential: is p=2 needed?
  lr_test(var_p2, var_p3),   # sequential: is p=3 needed?
  lr_test(var_p3, var_p4),   # sequential: is p=4 needed?
  lr_test(var_p3, var_p5),   # sequential: is p=4 needed?
  lr_test(var_p1, var_p4)    # joint:      is p=4 better than p=1?
)


cat("\n── LR tests for VAR lag order ──────────────────────────────\n")
cat("   χ²(df), same estimation sample for all models\n\n")
print(lr_results, row.names = FALSE)

p_opt <- 4 # based on LR tests

# 7.6 var with breaks, first ordering
# X_t = (MU, igae_growth, FU): baseline
# macro uncertainty is contemporaneously exogenous;
# financial uncertainty responds within-period to MU and output.
var_A_fit_breaks <- VAR(
  var_data_ts[,c("MU", "igae_growth", "FU")],
  p      = p_opt,
  type   = "const",
  exogen = exog_mat
)

irf_A  <- irf(var_A_fit_breaks, impulse = c("MU", "FU"), response = "igae_growth",
              n.ahead = 24, ortho = TRUE, boot = TRUE,
              ci = 0.90, runs = 500)

fevd_A <- fevd(var_A_fit_breaks, n.ahead = 24)


# 7.7 Ordering B: FU → MU → IGAE  ──────────────────────────────────────────────
var_B_fit_breaks <- VAR(var_data_ts[, c("FU", "MU", "igae_growth")],
             p = p_opt, type = "const", exogen = exog_mat)

irf_B  <- irf(var_B_fit_breaks, impulse = c("FU", "MU"), response = "igae_growth",
              n.ahead = 24, ortho = TRUE, boot = TRUE,
              ci = 0.90, runs = 500)

fevd_B <- fevd(var_B_fit_breaks, n.ahead = 24)


# 7.8. VARIANCE DECOMPOSITION SUMMARY  ───────────────────────────────────────
h12 <- 12   # report FEVD at 12-month horizon

fevd_A_igae <- as.data.frame(fevd_A$igae_growth) |>
  mutate(horizon = row_number()) |>
  filter(horizon == h12) |>
  dplyr::select(-horizon)

fevd_B_igae <- as.data.frame(fevd_B$igae_growth) |>
  mutate(horizon = row_number()) |>
  filter(horizon == h12) |>
  dplyr::select(-horizon)

cat(sprintf("\n── FEVD of IGAE growth at %d-month horizon ──\n", h12))
cat("Ordering A (MU first):\n")
print(round(fevd_A_igae, 3))
cat("Ordering B (FU first):\n")
print(round(fevd_B_igae, 3))

# 7.9 Augmented VAR with IPI_EUA ─────────────────────────────────────────────────
# Pull IPI_EUA from the transformed panel, aligned to uncertainty dates
ipi_dates  <- match(dates_fhat, rownames(db_dfm_trans))
ipi_eua_v  <- as.numeric(scale(db_dfm_trans[, "IPI_EUA"]))[ipi_dates]

svar_ext_df <- svar_df |>
  mutate(ipi_eua = ipi_eua_v[match(format(date, "%Y/%m"), dates_fhat)]) |>
  drop_na()

var_data_ext <- ts(
  svar_ext_df[, c("ipi_eua", "MU", "FU", "igae_growth")],
  start     = c(year(min(svar_ext_df$date)),
                month(min(svar_ext_df$date))),
  frequency = 12
)

var_ext_breaks  <- VAR(var_data_ext, p = p_opt, type = "const", exogen = exog_mat)
fevd_ext <- fevd(var_ext_breaks, n.ahead = 24)

as.data.frame(fevd_ext$igae_growth) |>
  mutate(horizon = row_number()) |>
  filter(horizon %in% c(1, 6, 12, 24)) |>
  dplyr::select(horizon, ipi_eua, MU, FU, igae_growth)


###### 8 Save all figures #####

# Colour palette consistent with LMN literature
col_mu <- "#1f4e79"
col_fu <- "#c55a11"

# ── Figure 1: MU and FU indices + IGAE growth ────────────────────────────────

# Join IGAE growth onto the uncertainty_indices date frame
p1_df <- uncertainty_indices |>
  left_join(svar_df |> dplyr::select(date, igae_growth), by = "date")

# Scale factor: map IGAE growth range onto MU/FU range
scale_factor_p1 <- max(abs(p1_df$MU), abs(p1_df$FU), na.rm = TRUE) /
  max(abs(p1_df$igae_growth),          na.rm = TRUE)

p1 <- p1_df |>
  ggplot(aes(x = date)) +
  # IGAE growth bars (background, secondary axis)
  geom_rect(
    data = imef |> filter(in_sample %in% c("Yes", "Partial")),
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "grey70", alpha = 0.25
  ) +
  geom_col(aes(y = igae_growth * scale_factor_p1),
           fill = "#70ad47", alpha = 0.35, width = 25) +
  # MU and FU lines (primary axis)
  geom_line(aes(y = MU, colour = "MU", linetype = "MU"), linewidth = 0.8) +
  geom_line(aes(y = FU, colour = "FU", linetype = "FU"), linewidth = 0.8) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
  scale_colour_manual(values = c(MU = col_mu, FU = col_fu)) +
  scale_linetype_manual(values = c(MU = "solid", FU = "dashed")) +
  scale_y_continuous(
    name     = "Standard deviations (MU / FU)",
    sec.axis = sec_axis(~ . / scale_factor_p1,
                        name = "IGAE monthly growth (%)")
  ) +
  labs(title    = "Macro and Financial Uncertainty Indices",
       subtitle = "Bars: IGAE monthly growth  |  Lines: MU (solid) and FU (dashed)",
       x = NULL, colour = NULL, linetype = NULL) +
  theme_bw(base_size = 11) +
  theme(
    legend.position        = "top",
    axis.title.y.right     = element_text(colour = "#70ad47")
  )

# ── Figure 2: MU – FU spread ─────────────────────────────────────────────────
p2 <- uncertainty_indices |>
  mutate(spread = MU - FU) |>
  ggplot(aes(x = date, y = spread)) +
  geom_ribbon(aes(ymin = pmin(spread, 0), ymax = 0),
              fill = col_fu, alpha = 0.45) +
  geom_ribbon(aes(ymin = 0, ymax = pmax(spread, 0)),
              fill = col_mu, alpha = 0.45) +
  geom_line(linewidth = 0.5, colour = "grey30") +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
  labs(title = "MU – FU Spread",
       subtitle  = paste0("Blue = macro-driven uncertainty dominates  |  ",
                          "Orange = finance-driven uncertainty dominates"),
       x = NULL, y = "MU − FU") +
  theme_bw(base_size = 11)

# ── Figure 3: IGAE growth and MU ─────────────────────────────────────────────
scale_factor <- max(abs(svar_df$igae_growth), na.rm = TRUE) /
  max(abs(svar_df$MU),          na.rm = TRUE)

p3 <- svar_df |>
  ggplot(aes(x = date)) +
  geom_col(aes(y = igae_growth), fill = "#70ad47", alpha = 0.6, width = 25) +
  geom_line(aes(y = MU * scale_factor), colour = col_mu, linewidth = 0.9) +
  scale_y_continuous(
    name = "IGAE monthly growth (pp)",
    sec.axis = sec_axis(~ . / scale_factor, name = "MU (standardised)")
  ) +
  labs(title = "IGAE Growth vs. Macro Uncertainty",
       x = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.title.y.right = element_text(colour = col_mu))

# ── Figure 4: IRFs under both orderings ──────────────────────────────────────
irf_to_tibble <- function(irf_obj, impulse_var, ordering_label) {
  resp  <- irf_obj$irf[[impulse_var]]
  lower <- irf_obj$Lower[[impulse_var]]
  upper <- irf_obj$Upper[[impulse_var]]
  tibble(
    horizon  = 0:(nrow(resp) - 1),
    estimate = resp[, "igae_growth"],
    lower    = lower[, "igae_growth"],
    upper    = upper[, "igae_growth"],
    impulse  = impulse_var,
    ordering = ordering_label
  )
}

irf_df <- bind_rows(
  irf_to_tibble(irf_A, "MU", "A: MU first (LMN baseline)"),
  irf_to_tibble(irf_A, "FU", "A: MU first (LMN baseline)"),
  irf_to_tibble(irf_B, "FU", "B: FU first (alternative)"),
  irf_to_tibble(irf_B, "MU", "B: FU first (alternative)")
)

p4 <- irf_df |>
  mutate(impulse = factor(impulse,
                          levels = c("MU", "FU"),
                          labels = c("MU shock → IGAE",
                                     "FU shock → IGAE"))) |>
  ggplot(aes(x = horizon, y = estimate, colour = impulse, fill = impulse)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15,
              colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  facet_wrap(~ ordering, ncol = 2) +
  scale_colour_manual(values = c("MU shock → IGAE" = col_mu,
                                 "FU shock → IGAE" = col_fu)) +
  scale_fill_manual(values   = c("MU shock → IGAE" = col_mu,
                                 "FU shock → IGAE" = col_fu)) +
  labs(title    = "IRF: Response of IGAE Growth to Uncertainty Shocks",
       subtitle = "90% bootstrap confidence bands | 500 replications",
       x = "Months after shock", y = "IGAE response",
       colour = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top",
        strip.background = element_rect(fill = "grey92"))

# ── Figure 5: FEVD at all horizons ───────────────────────────────────────────
fevd_long <- bind_rows(
  as.data.frame(fevd_A$igae_growth) |>
    mutate(horizon = row_number(), ordering = "A: MU first"),
  as.data.frame(fevd_B$igae_growth) |>
    mutate(horizon = row_number(), ordering = "B: FU first")
) |>
  pivot_longer(-c(horizon, ordering), names_to = "shock",
               values_to = "share")

p5 <- fevd_long |>
  filter(shock %in% c("MU", "FU")) |>
  ggplot(aes(x = horizon, y = share, colour = shock, linetype = ordering)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c(MU = col_mu, FU = col_fu)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title    = "FEVD: Share of IGAE Variance Explained by MU and FU",
       x = "Horizon (months)", y = "Fraction of forecast error variance",
       colour = NULL, linetype = "Ordering") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

EPU <- read_xlsx("../data/Mexico_Policy_Uncertainty_Data.xlsx")
EPU <- EPU[,c('Date', 'EPU_MX')]


# ── EPU comparison plot ────────────────────────────────────────────────────────
# Standardise EPU_MX to zero mean / unit variance so it is on the same scale
# as MU and FU, then merge on month-year key

EPU <- EPU |>
  mutate(
    date   = as.Date(as.yearmon(Date, "%Y/%m")),  # adjust format if needed
    EPU_MX = as.numeric(EPU_MX),
    EPU_MX = (EPU_MX - mean(EPU_MX, na.rm = TRUE)) / sd(EPU_MX, na.rm = TRUE)
  ) |>
  dplyr::select(date, EPU_MX)

# Build comparison data frame (inner join on date)
compare_df <- tibble(
  date = as.Date(as.yearmon(dates_fhat, "%Y/%m")),
  MU   = MU,
  FU   = FU
) |>
  inner_join(EPU, by = "date") |>
  pivot_longer(c(MU, FU, EPU_MX), names_to = "index", values_to = "value")

# Plot
epu_compare_plot <- ggplot(compare_df, aes(x = date, y = value, colour = index)) +
  geom_rect(
    data = imef,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "grey80", alpha = 0.4
  ) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  scale_colour_manual(
    values = c(MU = "#1f4e79", FU = "#c55a11", EPU_MX = "#538135"),
    labels = c(MU = "MU (this paper)", FU = "FU (this paper)", EPU_MX = "EPU (Baker et al.)")
  ) +
  labs(
    title    = "Uncertainty Indices — Mexico: LMN vs EPU",
    subtitle = "All series standardised (z-score) | Shaded: IMEF recessions",
    x = NULL, y = "Standard deviations", colour = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

png("../outputs/lmn/fig_epu_compare.png", width = 1950, height = 1050, res = 300)
print(epu_compare_plot)
dev.off()




dir.create(paste0("../outputs", "/lmn"), recursive = TRUE, showWarnings = FALSE)

ggsave(paste0("../outputs", "/lmn", "/fig1_mu_fu_indices.png"),
       p1, width = 10, height = 4, dpi = 150)

ggsave(paste0("../outputs", "/lmn", "/fig2_mu_fu_spread.png"),
       p2, width = 10, height = 4, dpi = 150)

ggsave(paste0("../outputs", "/lmn", "/fig3_igae_vs_mu.png"),
       p3, width = 10, height = 4, dpi = 150)

ggsave(paste0("../outputs", "/lmn", "/fig4_irfs.png"),
       p4, width = 12, height = 5, dpi = 150)

ggsave(paste0("../outputs", "/lmn", "/fig5_fevd.png"),
       p5, width = 9, height = 4, dpi = 150)



cat("\nFigures saved to outputs/lmn/\n")


###### 9. EXPORT RESULTS  #####
dir.create(paste0("../outputs", "/rds", "/lmn"), recursive = TRUE,
           showWarnings = FALSE)

lmn_results <- list(
  uncertainty_indices = uncertainty_indices,
  p1_df               = p1_df,
  igae_yoy            = igae_yoy,
  sigma_mat           = sigma_df,
  type_trans          = type_trans,
  fhat                = fhat_tr,
  Phat                = Phat_,
  e_idio              = e_idio,
  var_fit_A           = var_A_fit_breaks,
  var_fit_B           = var_B_fit_breaks,
  irf_A               = irf_A,
  irf_B               = irf_B,
  fevd_A              = fevd_A,
  fevd_B              = fevd_B,
  group_assignments   = list(
    macro         = macro_in_data,
    financial     = fin_in_data,
    nontraditional = nontr_in_data
  )
)

saveRDS(lmn_results,
        paste0("../outputs", "/rds", "/lmn", "/lmn_results.rds"))

write_csv(uncertainty_indices,
          paste0("../outputs", "/lmn", "/lmn_uncertainty_indices.csv"))

cat("\nResults saved:\n")
cat("  outputs/rds/lmn/lmn_results.rds\n")
cat("  outputs/lmn/lmn_uncertainty_indices.csv\n")
cat("\n=== Done ===\n")


# ── Unit root tests: MU and FU ──────────────────────────────────────────────
adf_mu  <- ur.df(ts(MU,  frequency = 12), type = "drift", selectlags = "AIC")
adf_fu  <- ur.df(ts(FU,  frequency = 12), type = "drift", selectlags = "AIC")
kpss_mu <- ur.kpss(ts(MU, frequency = 12), type = "mu")
kpss_fu <- ur.kpss(ts(FU, frequency = 12), type = "mu")

cat("ADF MU: ");  cat(round(adf_mu@teststat[1], 3),  "\n")
cat("ADF FU: ");  cat(round(adf_fu@teststat[1], 3),  "\n")
cat("KPSS MU: "); cat(round(kpss_mu@teststat, 3),    "\n")
cat("KPSS FU: "); cat(round(kpss_fu@teststat, 3),    "\n")


# Check GARCH persistence (alpha + beta) for all variables
# persistence <- sapply(sigma_df, function(fit)
#   sum(coef(fit)[c("alpha1", "beta1")])
# )
# cat("Max GARCH persistence:", round(max(persistence), 4), "\n")
# # Should be < 1 for all; any value ≥ 0.99 warrants a note
