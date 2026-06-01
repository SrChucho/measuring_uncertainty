rm(list = ls())
# =============================================================================
# GaR (Adrian et al. 2019 - ABG) — Uncertainty Indicators → IGAE Growth by quantile
# =============================================================================


# 0. LIBRARIES 
suppressMessages({
  library(here)
  library(tidyverse)
  library(vars)
  library(sandwich)
  library(lubridate)
  library(patchwork)
  library(quantreg)
  library(readxl)
})

# 1. LOAD RESULTS FROM MAIN SCRIPT 
res   <- readRDS(paste0("../outputs", "/rds", "/lmn", "/lmn_results.rds"))

# # simple quantile regression
# models <- list()
# qq <- seq(.05,.95, by = 0.05) 
# QR <- qq |> map(~ rq(igae_growth ~ MU + FU, data = res$p1_df, tau = .x))
# names(QR) <- qq

# dimensions
h <- c(1,3,6,12) # 4 horizons
qq <- seq(.05,.95, by = 0.05) # 19 quantiles
n_row <- nrow(res$p1_df)


# df_h <- data.frame(
#   Yh = res$p1_df$igae_growth[((h+1):n_row)],
#   res$p1_df[1:(n_row-h) ,c("MU","FU")]
# )
# 
# GaR <- qq |> map(~ rq(Yh ~ MU + FU, data = df_h, tau = .x))


igae_yoy  <- res$igae_yoy
n_row     <- length(igae_yoy)       


horizons <- c(1,3,6,12)
# this part follows QRboot function in AAB Replication package
df_list <- horizons |> 
  set_names(paste0("h_", horizons)) |> 
  map(function(h) {
    usable_rows <- n_row - h
    
    data.frame(
      Yh = igae_yoy[(h + 1):n_row],
      ylag = igae_yoy[1:usable_rows],
      res$p1_df[1:usable_rows, c("MU", "FU")]
    ) |> 
      dplyr::filter(!is.na(Yh))
  })

# Fit quantile regression for each forecasting horizon/quantile pair
GaR_models <- df_list |>
  map(function(df_h) {
    qq |>
      set_names(paste0("tau_", qq)) |>
     map(~ rq(Yh ~ ylag + MU + FU, data = df_h, tau = .x))  
  })
  


## coef analysis (sidequest)
h1coef <- sapply(GaR_models[[1]], coef)
h3coef <- sapply(GaR_models[[2]], coef)
h6coef <- sapply(GaR_models[[3]], coef)
h12coef <- sapply(GaR_models[[4]], coef)

coef_list <- list(h1coef, h3coef, h6coef, h12coef)

par(mfrow = c(2, 2), mar = c(4.5, 4, 3, 1), oma = c(2, 0, 3, 0))
for(i in 1:4){
  hh <- horizons[i]
  mat <- coef_list[[i]]

  y_range <- range(c(mat["MU", ], mat["FU", ]), na.rm = TRUE)
  
  plot(qq, mat["MU",], type = "l", col = 'red', pch = 1,
       lty = 2, lwd = 2, ylim= y_range,  ylab = "Coef", xlab = "Quantile",
       main = paste("Horizon h = ",hh))
  lines(qq, mat["FU",], col = 'blue', lty = 2, lwd = 2)
  
  grid(nx = NULL, ny = NULL, col = "gray90", lty = "solid")

}
mtext("Quantile Regression Coefficients across Forecasting Horizons", outer = TRUE, side = 3, line = 0.5, font = 2, cex = 1.2)
legend("bottom", legend = c("MU", "FU"), 
       col = c("red", "blue"), lty = 2, lwd = 2, 
       horiz = TRUE, bty = "n", cex = 1.1, xpd = TRUE, inset = c(0, -4.5))


# fitter  quantiles
fore_quantile1 <- sapply(GaR_models[[1]], predict)
fore_quantile3 <- sapply(GaR_models[[2]], predict)
fore_quantile6 <- sapply(GaR_models[[3]], predict)
fore_quantile12 <- sapply(GaR_models[[4]], predict)


fore_quantile_list <- list(
  '1' = fore_quantile1, 
  '3' = fore_quantile3, 
  '6' = fore_quantile6, 
  '12' = fore_quantile12
)


# Pad and align each horizon's fitted data frame
# each horizon-h matrix has (nrow-h) rows
# need to pad (add NA's) rows at the begining to align them all

aligned_horizons_list <- fore_quantile_list |> 
  imap(function(df_fitted, horizon_name) {
    
    num_padding <- n_row - nrow(df_fitted)
    

    padding_df <- as.data.frame(matrix(NA, nrow = num_padding, 
                                       ncol = ncol(df_fitted)))
    colnames(padding_df) <- colnames(df_fitted)
    
    aligned_df <- rbind(padding_df, df_fitted)
    
    # Rename columns to include the horizon prefix (e.g., h_1_tau_0.1)
    colnames(aligned_df) <- paste0("h",horizon_name, "_", colnames(aligned_df))
    
    return(aligned_df)
  })

# Combine all horizons side-by-side and bind them to your original dates
final_plot_df <- bind_cols(aligned_horizons_list) |> 
  mutate(
    Date = res$p1_df$date,                
    Actual_Growth = igae_yoy 
  ) |> 
  dplyr::select(Date, Actual_Growth, everything()) # Bring identifiers to the front



imef <- read_xlsx("../data/imef_recession_dates.xlsx") |>
  filter(phase == "recession", in_sample %in% c("Yes", "Partial")) |>
  mutate(
    start = as.Date(start),
    end   = if_else(is.na(end) | end == "", Sys.Date(), as.Date(end))
  )

plot_gar <- function(h_label, tau_lo, tau_med, tau_hi, title,
                     recessions = NULL,
                     ylim_range = NULL) {
  
  par(mfrow = c(1, 1), mar = c(6, 4, 4, 2))
  
  # ── 1. Auto ylim from quantile range (not raw IGAE which can be extreme) ────
  if (is.null(ylim_range)) {
    ylim_range <- range(
      c(final_plot_df[[tau_lo]], final_plot_df[[tau_hi]],
        final_plot_df$Actual_Growth),
      na.rm = TRUE
    )
    ylim_range <- ylim_range * 1.05   # 5% padding
  }
  
  # ── 2. Set up axes ───────────────────────────────────────────────────────────
  plot(final_plot_df$Date, final_plot_df$Actual_Growth,
       type = "n", ylim = ylim_range,
       xlab = "Time", ylab = "YoY IGAE growth (%)", main = title)
  
  # ── 3. Recession shading ─────────────────────────────────────────────────────
  if (!is.null(recessions) && nrow(recessions) > 0) {
    rect(xleft   = recessions$start,
         xright  = recessions$end,
         ybottom = par("usr")[3],
         ytop    = par("usr")[4],
         col     = adjustcolor("grey60", alpha.f = 0.30),
         border  = NA)
  }
  
  # ── 4. Grid ──────────────────────────────────────────────────────────────────
  grid(nx = NULL, ny = NULL, col = "gray93", lty = "solid")
  
  # ── 5. Zero reference line ───────────────────────────────────────────────────
  abline(h = 0, col = "grey30", lwd = 0.9, lty = 1)
  
  # ── 6. Fan ribbon: shaded area between 10th and 90th percentile ─────────────
  valid <- !is.na(final_plot_df[[tau_lo]]) & !is.na(final_plot_df[[tau_hi]])
  d_v   <- final_plot_df$Date[valid]
  lo_v  <- final_plot_df[[tau_lo]][valid]
  hi_v  <- final_plot_df[[tau_hi]][valid]
  
  polygon(c(d_v, rev(d_v)),
          c(hi_v, rev(lo_v)),
          col    = adjustcolor("steelblue", alpha.f = 0.13),
          border = NA)
  
  # ── 7. Data lines ────────────────────────────────────────────────────────────
  lines(final_plot_df$Date, final_plot_df$Actual_Growth,
        col = "grey40", lwd = 1.5, lty = 1)
  lines(final_plot_df$Date, final_plot_df[[tau_hi]],
        col = "steelblue", lwd = 1.2, lty = 2)
  lines(final_plot_df$Date, final_plot_df[[tau_med]],
        col = "black",     lwd = 1.2, lty = 2)
  lines(final_plot_df$Date, final_plot_df[[tau_lo]],
        col = "firebrick", lwd = 2.0, lty = 1)
  
  # ── 8. Legend ────────────────────────────────────────────────────────────────
  legend("bottom",
         legend = c("IGAE yoy", "5th-p (GaR)", "50th-p", "95th-p"),
         lty    = c(1, 1, 2, 2),
         lwd    = c(1.5, 2.0, 1.2, 1.2),
         col    = c("grey40", "firebrick", "black", "steelblue"),
         horiz  = TRUE, bty = "n",
         inset  = c(0, -0.22), xpd = TRUE)
}


# h = 1
png(paste('../outputs/abg/','fan_plot_h1.png', sep = ""),
    width = 1080, height = 780, res = 150)
plot_gar("1",  "h1_tau_0.05",  "h1_tau_0.5",  "h1_tau_0.95",
         "Growth-at-Risk (h = 1 month)",  recessions = imef)
dev.off()

# h = 3
png(paste('../outputs/abg/','fan_plot_h3.png', sep = ""),
    width = 1080, height = 780, res = 150)
plot_gar("3", "h3_tau_0.05", "h3_tau_0.5", "h3_tau_0.95",
         "Growth-at-Risk (h = 3 months)", recessions = imef)
dev.off()

# h = 6
png(paste('../outputs/abg/','fan_plot_h6.png', sep = ""),
    width = 1080, height = 780, res = 150)
plot_gar("6", "h6_tau_0.05", "h6_tau_0.5", "h6_tau_0.95",
         "Growth-at-Risk (h = 6 months)", recessions = imef)
dev.off()

# h = 12
png(paste('../outputs/abg/','fan_plot_h12.png', sep = ""),
    width = 1080, height = 780, res = 150)
plot_gar("12", "h12_tau_0.05", "h12_tau_0.5", "h12_tau_0.95",
         "Growth-at-Risk (h = 12 months)",recessions = imef)
dev.off()


