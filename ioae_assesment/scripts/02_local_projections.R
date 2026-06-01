rm(list = ls())
# =============================================================================
# Local Projections (Jordà 2005) — Uncertainty Shocks → IGAE Growth
# =============================================================================
# Identification: Cholesky, LMN ordering  X_t = (MU_t, IGAE_t, FU_t)'
#
#   MU  first  → macro uncertainty is contemporaneously exogenous;
#                a MU shock can shift output and financial uncertainty
#                within the same period.
#   FU  last   → financial uncertainty is contemporaneously endogenous;
#                it absorbs MU and output shocks on impact.
#
# NOTE: the main script (01_jln_lmn_methodology.R) uses (MU, FU, IGAE)
#       in Ordering A, which is NOT the LMN baseline. Use the VAR estimated
#       here for any inference about causal ordering.
#
# LP advantages over VAR-based IRFs:
#   - No mis-specification accumulation across horizons
#   - Confidence bands via Newey-West HAC (bandwidth = horizon)
#   - Robust to lag-length misspecification
# =============================================================================

# ── 0. LIBRARIES ──────────────────────────────────────────────────────────────
suppressMessages({
  library(here)
  library(tidyverse)
  library(vars)
  library(sandwich)
  library(lubridate)
  library(patchwork)
})

# ── 1. LOAD RESULTS FROM MAIN SCRIPT ──────────────────────────────────────────
res   <- readRDS(paste0("../outputs", "/rds", "/lmn", "/lmn_results.rds"))
p_opt <- res$var_fit_A$p                       # lag order from main script

# Recover original VAR data (T × 3: MU, FU, igae_growth)
var_df <- as.data.frame(res$var_fit_A$y)

cat(sprintf("VAR data: %d obs, lag order p = %d\n", nrow(var_df), p_opt))


# ── 2. RE-ESTIMATE VAR WITH CORRECT LMN ORDERING ──────────────────────────────
# X_t = (MU, igae_growth, FU)
var_lmn <- VAR(
  var_df[, c("MU", "igae_growth", "FU")],
  p    = p_opt,
  type = "const"
)

cat("VAR (LMN ordering) estimated.\n")


# ── 3. CHOLESKY STRUCTURAL SHOCKS ─────────────────────────────────────────────
eta   <- residuals(var_lmn)                    # T_eff × 3 reduced-form residuals u_t = A(L)X_t
Sigma <- crossprod(eta) / nrow(eta)            # Σ_u = E[u_t u_t'] — residual covariance
P     <- t(chol(Sigma))                        # Σ_u = PP', P lower-triangular
eps   <- t(solve(P) %*% t(eta))                # ε_t = P⁻¹ u_t — orthogonal structural 
colnames(eps) <- colnames(eta)

# Sanity check — off-diagonals should be near zero
cat("\nStructural shock covariance (should be ~identity):\n")
print(round(cov(eps), 3))


# ── 4. LP HELPER FUNCTION ─────────────────────────────────────────────────────
# For each horizon h, estimates:
#   y_{t+h} = alpha_h + beta_h * shock_t + gamma_h * controls_t + u_{t+h}
#
# Standard errors: Newey-West HAC with bandwidth = max(1, h) to account for
# the MA(h-1) serial correlation induced by the overlapping projection window.

run_lp <- function(y, shock, controls, H = 24, level = 0.90) {
  T_lp <- length(y)
  z    <- qnorm(1 - (1 - level) / 2)

  map_dfr(0:H, function(h) {
    idx_t  <- seq_len(T_lp - h)
    idx_th <- idx_t + h

    df_h <- data.frame(
      y_h   = y[idx_th],
      shock = shock[idx_t],
      controls[idx_t, , drop = FALSE]
    )

    # IGAE_{t+h} = α_h + β_h · ε_t + γ_h · controls_t + error_{t+h}
    fit  <- lm(y_h ~ ., data = df_h)
    V_nw <- NeweyWest(fit, lag = max(1, h), prewhite = FALSE, adjust = TRUE)
    b    <- coef(fit)["shock"]
    se   <- sqrt(V_nw["shock", "shock"])

    tibble(
      horizon  = h,
      estimate = b,
      lower    = b - z * se,
      upper    = b + z * se
    )
  })
}


# ── 5. BUILD LP INPUTS ────────────────────────────────────────────────────────
T_eff <- nrow(eta)                             # effective sample after lag trimming
n_full <- nrow(var_df)

# Response: igae_growth aligned to structural shock sample
y_lp <- tail(var_df$igae_growth, T_eff)

# Controls: p lags of all three variables
lag_mat <- do.call(cbind, lapply(seq_len(p_opt), function(j) {
  blk <- var_df[(p_opt - j + 1):(n_full - j), , drop = FALSE]
  setNames(blk, paste0(colnames(var_df), "_lag", j))
}))


# ── 6. RUN LOCAL PROJECTIONS ──────────────────────────────────────────────────
H <- 24

cat("\nRunning LP for MU shock...\n")
lp_mu <- run_lp(y_lp, eps[, "MU"], lag_mat, H) |>
  mutate(impulse = "MU shock → IGAE", method = "LP")

cat("Running LP for FU shock...\n")
lp_fu <- run_lp(y_lp, eps[, "FU"], lag_mat, H) |>
  mutate(impulse = "FU shock → IGAE", method = "LP")


# ── 7. VAR IRFs FOR COMPARISON ────────────────────────────────────────────────
cat("Bootstrapping VAR IRFs (500 runs)...\n")
irf_lmn <- irf(
  var_lmn,
  impulse  = c("MU", "FU"),
  response = "igae_growth",
  n.ahead  = H,
  ortho    = TRUE,
  boot     = TRUE,
  ci       = 0.90,
  runs     = 500
)

irf_to_tbl <- function(irf_obj, impulse_var) {
  tibble(
    horizon  = 0:H,
    estimate = as.numeric(irf_obj$irf[[impulse_var]][,   "igae_growth"]),
    lower    = as.numeric(irf_obj$Lower[[impulse_var]][, "igae_growth"]),
    upper    = as.numeric(irf_obj$Upper[[impulse_var]][, "igae_growth"]),
    impulse  = paste0(impulse_var, " shock → IGAE"),
    method   = "VAR (Cholesky)"
  )
}

var_mu <- irf_to_tbl(irf_lmn, "MU")
var_fu <- irf_to_tbl(irf_lmn, "FU")


# ── 8. PRINT SUMMARIES ────────────────────────────────────────────────────────
fmt <- function(df) {
  df |>
    filter(horizon %in% c(0, 3, 6, 12, 18, 24)) |>
    dplyr::select(horizon, estimate, lower, upper) |>
    mutate(across(where(is.numeric), \(x) round(x, 4)))
}

cat("\n── LP IRF: MU shock → IGAE (90% NW CI) ──\n");  print(fmt(lp_mu))
cat("\n── LP IRF: FU shock → IGAE (90% NW CI) ──\n");  print(fmt(lp_fu))

cat("\n── VAR FEVD of IGAE (LMN ordering, h = 12) ──\n")
fevd_lmn <- fevd(var_lmn, n.ahead = 24)
print(round(as.data.frame(fevd_lmn$igae_growth)[12, ], 4))


# ── 9. PLOT ───────────────────────────────────────────────────────────────────
col_mu <- "#1f4e79"
col_fu <- "#c55a11"

plot_df <- bind_rows(lp_mu, lp_fu, var_mu, var_fu) |>
  mutate(
    method  = factor(method,  levels = c("LP", "VAR (Cholesky)")),
    impulse = factor(impulse, levels = c("MU shock → IGAE",
                                         "FU shock → IGAE"))
  )

p_lp <- ggplot(plot_df,
               aes(x = horizon, y = estimate,
                   colour = impulse, fill = impulse,
                   linetype = method)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  facet_wrap(~ impulse, ncol = 2) +
  scale_colour_manual(values = c("MU shock → IGAE" = col_mu,
                                 "FU shock → IGAE" = col_fu)) +
  scale_fill_manual(values   = c("MU shock → IGAE" = col_mu,
                                 "FU shock → IGAE" = col_fu)) +
  scale_linetype_manual(values = c("LP"             = "solid",
                                   "VAR (Cholesky)" = "dashed")) +
  labs(
    title    = "IRF: IGAE Response to Uncertainty Shocks",
    subtitle = paste0("Solid: Local Projections (90% Newey-West CI)  |  ",
                      "Dashed: VAR Cholesky (90% bootstrap)\n",
                      "Ordering: X_t = (MU, IGAE, FU)'  — LMN baseline"),
    x        = "Months after shock",
    y        = "IGAE growth response",
    colour   = NULL, fill = NULL, linetype = "Method"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "top",
    strip.background = element_rect(fill = "grey92")
  )

dir.create(paste0("../outputs", "/lmn"), recursive = TRUE, showWarnings = FALSE)
ggsave(paste0("../outputs", "/lmn", "/fig6_lp_irf.png"),
       p_lp, width = 12, height = 5, dpi = 150)

cat("\nFigure saved: outputs/lmn/fig6_lp_irf.png\n")
cat("=== Done ===\n")
