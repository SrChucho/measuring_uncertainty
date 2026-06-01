library(glmnet) |> suppressMessages()
library(forecast) |> suppressMessages()
library(nowcasting) |> suppressMessages()
library(caret) |> suppressMessages()
library(nnfor) |> suppressMessages()
library(vars) |> suppressMessages()
library(tidyverse) |> suppressMessages()
library(murphydiagram) |> suppressMessages()
library(MCS) |> suppressMessages()

# -------------------------------
# preparing output structure
prepared_inputs <- readRDS(here::here("outputs", "rds", "prepared_inputs", "prepared_inputs.rds"))
list2env(prepared_inputs, envir = parent.frame())

dir.create(here::here("outputs", "rds", "figures_inputs"),
           recursive = TRUE, showWarnings = FALSE)

dir.create(here::here("outputs", "rds", "models_objects"),
           recursive = TRUE, showWarnings = FALSE)
# -------------------------------


# -------------------------------
# parameters
pibo <- TRUE
machine_learning <- TRUE
lasso_reg <- TRUE
if(rolling) w <- "roll" else w <- "expa"
# -------------------------------

  ########################### current nowcasting ####################################
  
  # nowcasting horizon
  H <- sum(is.na(yfcst))
  
  # using PIBO
  yfcst_o <- gdp_o(yfcst_levels, yfcst, pib, d, vari = vari_o)
  
  
  if(!is.null(yfcst_o)){
    Y <- cbind(db_dfm_p, yfcst_o)
    colnames(Y) <- c(colnames(db_dfm_p), vari)
  }else{
    Y <- db_dfm_p
  }
  
  db_now <- ts(cbind(yfcst, scale(Y)), start = year_s, frequency = 12)
  colnames(db_now) <- c("YFCST", colnames(Y))
  
  ############ number of factors #####
  n_fts_Bai_Ng <- t(sapply(1:Ht, function(j) bai_ng(Y[1:(Ty-Ht+j),], demean = 2)$rhat))
  colnames(n_fts_Bai_Ng) <- paste0("Bai_Ng_", colnames(n_fts_Bai_Ng))
  n_fts_Bai_Ng <- as_tibble(n_fts_Bai_Ng)
  
  n_fts_Onatski <- sapply(1:Ht, function(j) onatski2010(Y[1:(Ty-Ht+j),])["ed"])
  n_fts_Onatski <- tibble(Onatski_ed = n_fts_Onatski)
  
  n_fts_Ahn_Hor <- t(sapply(1:Ht, function(j) ratio.test(Y[1:(Ty-Ht+j),])[c("ker", "kgr")]))
  colnames(n_fts_Ahn_Hor) <- paste0("Ahn_Horenstein_", colnames(n_fts_Ahn_Hor))
  n_fts_Ahn_Hor <- as_tibble(n_fts_Ahn_Hor)
  
  df_n_factors <- tibble(
    Date = zoo::as.Date.yearmon((time(db_dfm_p)))[(Ty-Ht+1):(Ty)]
  ) |>
    bind_cols(
      n_fts_Bai_Ng,
      n_fts_Onatski,
      n_fts_Ahn_Hor
    )
  
  
  if(rhat == 1 & w == "expa") {
    cat("\nNumber of factors determination under different criteria, November 2018–November 2024\n")
    
    df_n_factors|> print(n = 200)
  
    df_n_factors|>
      saveRDS(here::here("outputs/rds/models_objects/df_n_factors.rds"))
  
  }
  
  freq_now <- rep(12, ncol(db_now))
  
  nowcast_current <- nowcast(YFCST ~., data = db_now, r = rhat, p = 1, 
                             q = rhat, method = "2s", frequency =  freq_now)
  
  # objects and nowcastings
  fhat <- as.matrix(nowcast_current$factors$dynamic_factors)
  model <- Arima(yfcst[1:Ty], c(0,0,0), xreg = fhat[1:Ty,])
  nowcast <- forecast(model, xreg = fhat[-(1:Ty),])
  
  confidence <- confidence_dfm(nowcast_current, rhat = rhat)
  
  Phat <- confidence$Phat[[1]]
  Names <- cat_dfm[cat_dfm[,vari] == 1, "Description"]
  Names[Names == "Trend U Manuf"] <- "Trend L Manuf"
  
  ################## pooled and ADF tests #######
  db_now_ <- db_now[1:(nrow(db_now)-H),]
  
  # regressors
  X_ <- db_now_[,-1]
  
  # loadings
  Phat_ <- nowcast_current$factors$eigen$vectors[,1:rhat,drop=FALSE]
  rownames(Phat_) <- colnames(db_now_)[-1]
  
  fhat_ <- as.matrix(nowcast_current$factors$dynamic_factors) 
  fhat_ <- as.matrix(fhat_[1:nrow(X_),])
  
  # factor rotation
  if(cor(db_now_[,1], fhat_[,1], use ="pairwise.complete.obs") < 0){
    Phat_ <- -Phat_
    fhat_ <- -fhat_
  }
  
  FP <- fhat_ %*% t(Phat_)
  e_idio <- X_ - FP
  
  pt_comm <- pooled.test(FP) |> suppressWarnings()
  pt_idio <- pooled.test(e_idio) |> suppressWarnings()
  
  
  
  if(rhat == 1 & w == "expa") {
    
    ##### factor ADF test #####
    lags <- sapply(0:7, function(x) adf(fhat, "const", k = x)$bic)
    
    f_adf <- adf(fhat, "const", k = which(lags == min(lags))-1)
    ###########################
    
    df_fct_adf_test_pv <- tibble(
      DF_statistic = f_adf$statistic |> round(3),
      p_value = f_adf$p.value |> round(3)
    )
    
    df_fct_adf_test_pv |>
      saveRDS("outputs/rds/models_objects/df_fct_adf_test_pv.rds")
    
    cat("\nP-value of common factor ADF test\n")
    
    df_fct_adf_test_pv |> print()
    
    df_pooled_tests_pv <- tibble(
      P_statistic = pt_idio$stats[1] |> round(3),
      p_value = pt_idio$stats[2] |> round(3)
    )
  
    df_pooled_tests_pv |> 
      saveRDS("outputs/rds/models_objects/df_pooled_tests_pv.rds")
  
    cat("\nP-value of idyocincratic component pooled test\n")
    
    df_pooled_tests_pv |> print()
    
  }
  
  
  #################
  
  if(pibo)
    Names <- c(Names, "PIBO")
  
  rownames(Phat) <- Names
  
  ######## Variables and correlations ####
  TS <- data.frame(Names)
  saveRDS(TS, here::here("outputs/rds/figures_inputs/TS.rds"))
  
  rownames(type_trans) <- Names[-length(Names)]
  
  df_type_trans <- type_trans |>
    as_tibble(
      rownames = "Variable"
    ) |> 
    mutate(
      Transformation = case_when(
        Transformation == "monthly" ~ "MV",
        Transformation == "annual" ~ "AV",
        TRUE ~ "Level"
      )
    )
  
  if(rhat == 1 & w == "expa") {
    cat("\nTable 1: IOAE macroeconomic and financial time series\n")
    
    df_type_trans|> print(n = 200)
    cat("\n")
  
  df_type_trans|>
    saveRDS(here::here("outputs/rds/models_objects/Trans.rds"))
  
  }
  
  ######################### comparative performance assessment ##############################
  set.seed(12345)
  ####
  avai_t2 <- read.csv(here::here("Data/DFM/Catalogo.csv"))
  v_na_t2 <- avai_t2$Short[avai_t2[,"Avai_T2"] != 1]
  ####
  
  # NA indicator
  na_ind <- is.na(db_dfm_raw[(nrow(db_dfm_raw)+1-H):
                               (nrow(db_dfm_raw)),,drop=FALSE])
    
  # T+H sequence
  out_sample <- matrix(0, Ht, H)
  for(h in (H-1) : 0) 
    out_sample[,H-h] <- (Ty-Ht-h+1):(Ty-h)
    
  ####
  yfcst_dates <- names(db_yfcst[!is.na(db_yfcst)])
  dates_matrix <- matrix(0, nrow = nrow(out_sample), ncol = 2)
  for(i in 1:nrow(out_sample)) {
    dates_matrix[i, 1] <- yfcst_dates[out_sample[i, 1]]
    dates_matrix[i, 2] <- yfcst_dates[out_sample[i, 2]]
  }
  ####
    
  # models_test
  if(rhat == 1 & w == "expa") {
    dfm_mae <- arima_mae <- arimax_mae <- mlp_mae <- lasso_mae <- favar_mae <- 
      matrix(NA, Ht, H)  
  } else {
    dfm_mae <- matrix(NA, Ht, H)
  }

  fhat_list <- empirical_mae <- loadings_dfm <- nowcasts_models <- 
    coeff_lasso_sig <- list()
    
  for(h in 1 : Ht){ 
    
    cat(glue::glue("\nNowcast in progress: {h}/{Ht}\r"))
    flush.console()

    fhat_models <- matrix(NA, Ty, rhat)
      
    # historical data set
    into <- 1:(out_sample[h,H])
      
    if(rolling)
      into <- h:(out_sample[h,H])
      
    yfcst_h <- yfcst[into]
    yfcst_h_levels <- yfcst_levels[into]
    
    if(rolling){
      yfcst_h[out_sample[1,]] <- NA
      yfcst_h_levels[out_sample[1,]] <- NA
    }else{
      yfcst_h[out_sample[h,]] <- NA
      yfcst_h_levels[out_sample[h,]] <- NA
    }
      
    Y_h <- db_dfm_p[into,]
      
    if(rolling){
      google_db <- google[h:nrow(google),]
    }else{
      google_db <- google
    }
      
    # current NAs
    if(rolling){
      for(i in 1 : H) 
        Y_h[out_sample[1,][i],][na_ind[i,]] <- NA
    }else{
      for(i in 1 : H) 
        Y_h[out_sample[h,][i],][na_ind[i,]] <- NA
    }
      
    Y_h <- na.locf(Y_h)
      
    if(pibo){ 
      # PIBO
      yfcst_o_h <- gdp_o(yfcst_h_levels, yfcst_h, pib, d, vari = vari_o)
        
      if(rolling){
        if(!is.null(yfcst_o_h)){
          Y_h <- cbind(Y_h, yfcst_o_h[h:length(yfcst_o_h)])
            colnames(Y_h) <- c(colnames(db_dfm_p), vari)
        }    
      }else{
         if(!is.null(yfcst_o_h)){
          Y_h <- cbind(Y_h, yfcst_o_h)
          colnames(Y_h) <- c(colnames(db_dfm_p), vari)
        }
      }
    }
      
    Ty_h <- sum(!is.na(yfcst_h))
      
    # generating the cross-section time series validation
    time_slices_h <- createTimeSlices(1:Ty_h, initialWindow = Ty_h-Ht, 
                                      horizon = H, fixedWindow = FALSE)
      
    if(rolling){
      time_slices_h <- createTimeSlices(1:Ty_h, initialWindow = Ty_h-Ht, 
                                        horizon = H, fixedWindow = TRUE)
    }
      
    train_control_h <- trainControl(method = "cv", index = time_slices_h$train, 
                                    indexOut = time_slices_h$test)
      
    lasso_google_h <- caret::train(x = as.data.frame(google[1:Ty_h, ]), 
                                   y=yfcst_h[1:Ty_h],  
                                   method = "glmnet", 
                                   trControl = train_control_h, tuneLength = 10, 
                                   tuneGrid = expand.grid(alpha = 1, 
                                                          lambda = seq(0.0001, 0.1, length.out = 100)))
      
    # best lambda
    best_lambda_google_h <- lasso_google_h$bestTune["lambda"]
      
    # best regression
    best_lasso_google_h <- glmnet(x = as.data.frame(google_db[1:Ty_h, ]), 
                                  y=yfcst_h[1:Ty_h], alpha = 1, 
                                  lambda = as.numeric(best_lambda_google_h))
      
    coef_google <- rownames(coef(best_lasso_google_h))[
      which(abs(coef(best_lasso_google_h)) > 0)][-1]
      
    if(length(coef_google) > 0)
      Y_h <- cbind(Y_h, google_db[1:(Ty_h+H),coef_google,drop=FALSE])
      
    # nowcasting model
    db_now_h <- ts(cbind(yfcst_h, scale(Y_h)), start = year_s, frequency = 12)
    colnames(db_now_h) <- c("YFCST", colnames(Y_h))
    
    freq_now <- rep(12, ncol(db_now_h))
      
    nowcast_object <- nowcast(YFCST ~., data = db_now_h, r = rhat, p = 1, 
                                q = rhat, method = "2s", frequency =  freq_now)
    
    # loadings
    Phat <- nowcast_object$factors$eigen$vectors[,1:rhat,drop=FALSE]
    rownames(Phat) <- colnames(Y_h)
    #rownames(Phat) <- colnames(Y_h)[!(colnames(Y_h) %in% v_na_t2)]
      
    fhat <- as.matrix(nowcast_object$factors$dynamic_factors)
      
    # factor rotation
    if(cor(yfcst_h, fhat[,1], use ="pairwise.complete.obs") < 0){
      Phat <- -Phat
      fhat <- -fhat
    }
      
    # factors and loadings
    loadings_dfm[[h]] <- Phat
      
    fhat_models[1:nrow(fhat), ] <- fhat
    fhat_list[[h]] <- fhat_models
      
    # nowcasting estimation
    fhat_h <- as.matrix(nowcast_object$factors$dynamic_factors)
    model_h <- Arima(yfcst_h[1:Ty_h], c(0,0,0), xreg = fhat_h[1:Ty_h,])
    nowcast_h <- forecast(model_h, xreg = fhat_h[-(1:Ty_h),])
      
    # nowcasting object
    nowcast_dfm <- matrix(cbind(nowcast_h$mean, nowcast_h$lower[, "95%"],
                                nowcast_h$upper[, "95%"]), nrow = H, ncol = 3)
      
    colnames(nowcast_dfm) <- c("mean", "lower", "upper")
    rownames(nowcast_dfm) <- names(yfcst_h_levels)[is.na(yfcst_h_levels)]
      
    nowcasts_dfm <- lin_function(yfcst_h_levels, nowcast_dfm, d = d)
      
    nowcasts_dfm$levels <- cbind(nowcasts_dfm$levels, 
                                 yfcst_levels[names(yfcst_h_levels)
                                                [is.na(yfcst_h_levels)]])
      
    nowcasts_dfm$monthly <- cbind(nowcasts_dfm$monthly, 
                                  yfcst[names(yfcst_h_levels)
                                          [is.na(yfcst_h_levels)]])
      
    nowcasts_dfm$annual <- cbind(nowcasts_dfm$annual, 
                                 yfcst_annual[names(yfcst_h_levels)
                                                [is.na(yfcst_h_levels)]])
      
    colnames(nowcasts_dfm$levels) <- colnames(nowcasts_dfm$monthly) <-
      colnames(nowcasts_dfm$annual) <- c("mean", "lower", "upper", "observed")
      
    nowcasts_models[[h]] <- nowcasts_dfm
      
    # dfm error
    dfm_mae[h, ] <- nowcast_h$mean - yfcst[out_sample[h,]]
  
    
    if(rhat == 1 & w == "expa") {    
    
      # arima
      arima_mae[h,] <- forecast(auto.arima(yfcst_h[1:Ty_h]), h = H)$mean -
        yfcst[out_sample[h,]]      
      
      # arimax
      arimax_mae[h,] <- forecast(auto.arima(yfcst_h[1:Ty_h], 
                                            xreg = fhat_h[1:Ty_h,]), xreg = fhat_h[-(1:Ty_h),])$mean -
        yfcst[out_sample[h,]]     
      
      # neural network
      if(machine_learning){
        y_into <- ts(yfcst_h[1:Ty_h], start = year_s, frequency = 12)
        fhat_into <- as.matrix(fhat_h[1:Ty_h,])
        
        mlp_m <- mlp(y_into, xreg = fhat_into)
        
        fhat_out <- as.matrix(fhat_h)
        
        fore_mlp <- forecast(mlp_m, h = H, xreg = fhat_out)$mean
        
        mlp_mae[h,] <- fore_mlp - yfcst[out_sample[h,]]
      }
      
      # empirical model
      var_h <- colnames(Y_h)[which(apply(Y_h, 2, sd) != 0) & 
                               cor(yfcst_h, Y_h, use = "pairwise.complete.obs") != 1]
      
      Y_h_o <- Y_h[,var_h]
      
      nowcasts_var_h <- matrix(0, H, length(var_h))
      colnames(nowcasts_var_h) <- var_h
      
      for(i in 1 : length(var_h)){
        model_emp_h <-  try(Arima(yfcst_h[1:Ty_h], c(0,0,0), 
                                  xreg = Y_h_o[1:Ty_h, i]), silent = TRUE)
        
        
        nowcasts_var_h[,i] <- forecast(model_emp_h, xreg = Y_h_o[-(1:Ty_h),i])$mean 
      }
      
      empirical_mae[[h]] <- nowcasts_var_h - yfcst[out_sample[h,]]
      
      # model
      if(lasso_reg){
        # LASSO regression
        
        # training the model
        lasso_model_h <- caret::train(x = as.data.frame(db_now_h[1:Ty_h,-1]), 
                                      y=db_now_h[1:Ty_h,1],  
                                      method = "glmnet", 
                                      trControl = train_control_h, tuneLength = 10, 
                                      tuneGrid = expand.grid(alpha = 1, 
                                                             lambda = seq(0.0001, 0.1, length.out = 100)))
        
        # best lambda
        best_lambda_h <- lasso_model_h$bestTune["lambda"]
        
        # best regression
        lasso_best_h <- glmnet(y = c(db_now_h[1:Ty_h,1]), 
                               x = db_now_h[1:Ty_h,-1], alpha = 1, 
                               lambda = as.numeric(best_lambda_h))
        
        
        coef_lasso <- rownames(coef(lasso_best_h))[which(abs(coef(lasso_best_h)) 
                                                         > 0)][-1]
        
        
        coeff_lasso_sig[[h]] <- coef_lasso
        
        # prediction
        pred <- predict(lasso_best_h, 
                        as.matrix(db_now_h[(Ty_h+1):(Ty_h+H),-1]))
        
        lasso_mae[h,] <- pred - yfcst[out_sample[h,]]
      }
      
      # FAVAR [Corona et al. (2017) with contemporaneous effects]
      dfm <- cbind(yfcst_h[1:Ty_h], fhat_h[1:Ty_h,])
      colnames(dfm) <- c("YFCST", paste("fhat", 1: rhat, sep = ""))
      
      exogen <- as.matrix(fhat_h)
      colnames(exogen) <- paste("exo",1:rhat, sep = "")
      
      p <- VARselect(dfm, exogen = exogen[1:Ty_h,,drop=FALSE])$selection["SC(n)"]
      var_1 <- VAR(dfm, p = p, exogen =  exogen[1:Ty_h,,drop=FALSE]) 
      
      restricted <- matrix(0, rhat + 1, 1 + (rhat+1)*p + rhat)
      colnames(restricted) <- rownames(coef(var_1)$YFCST)
      
      for(j in 1 : rhat){
        restricted[j+1, paste("fhat",j,".l",1:p, sep = "")] <- 1 
      }
      restricted[1,] <- 1
      restricted[,"const"] <- 1
      
      var_1r <- restrict(var_1, method = "manual", resmat = restricted)
      
      favar_predict <- predict(var_1r, n.ahead = H, 
                               dumvar = exogen[-(1:Ty_h),,drop=FALSE])[[1]]$YFCST[,"fcst"]
      
      favar_mae[h,] <- favar_predict - yfcst[out_sample[h,]]
      
    }
  }  

  if(rhat == 1 & w == "expa") {
      
    colnames(Y_h)[!is.element(colnames(Y_h), coef_lasso)]
      
    # selecting models
    K_names <- sapply(1:ncol(db_dfm_p), 
                        function(x) colnames(empirical_mae[[x]])[1:ncol(db_dfm_p)])
    K <- apply(K_names, 2, length)
    K_w <- which(K == min(K))[1] 
    k <- K[K_w]
    k_names <- K_names[,K_w]
    
    ####
    v_avai_t2 <- !(k_names %in% v_na_t2)
    k_names <- k_names[v_avai_t2]
    ####
    
    
    ######## Getting errors ########
    errors_list <- list()
    for(i in 1 : length(k_names)){
      empirical_errors <- matrix(0, Ht, H)
      for(h in 1 : Ht){
        empirical_errors[h,] <- empirical_mae[[h]][,k_names[i]]
      }
      errors_list[[i]] <- empirical_errors
    }
    names(errors_list) <- k_names
    
    # generates dfm errors data frame
    f_err_df <-function(v_err, metric = c("MAE", "RMSE"), rhat, w) {
      tibble(
        metric = rep(metric,2),
        h = c(1:2),
        r = rep(rhat, 2),
        window = rep(w, 2),
        value = v_err
      )
    }
    
    # MAE
    mae_empirical <- sapply(1:length(k_names), 
                            function(x) colMeans(abs(errors_list[[x]])))
    colnames(mae_empirical) <- k_names
  
  }
  
  # errors
  colMeans(abs(dfm_mae))
  colMeans(abs(mlp_mae))
  colMeans(abs(arima_mae))
  colMeans(abs(arimax_mae))
  colMeans(abs(lasso_mae))
  colMeans(abs(favar_mae))
  mae_empirical
  
  v_mae_dfm <- colMeans(abs(dfm_mae))
  df_mae_dfm <- v_mae_dfm |>
    f_err_df("MAE", rhat, w)
  
  # RMSE
  rmse_empirical <- sapply(1:length(k_names), 
                           function(x) sqrt(colMeans(errors_list[[x]]^2)))
  colnames(rmse_empirical) <- k_names
  
  #errors
  sqrt(colMeans(dfm_mae^2))
  sqrt(colMeans(mlp_mae^2))
  sqrt(colMeans(arima_mae^2))
  sqrt(colMeans(arimax_mae^2))
  sqrt(colMeans(lasso_mae^2))
  sqrt(colMeans(favar_mae^2))
  rmse_empirical
  
  v_rmse_dfm <- sqrt(colMeans(dfm_mae^2))
  df_rmse_dfm <- v_rmse_dfm |>
    f_err_df("RMSE", rhat, w)
  
  ## for selected dates 
  months <- paste("0", 3:6, sep = "")
  dates_sel <- paste("2020/", months, sep = "")
  
  # indicator of months
  ind_dates <- sapply(1:H, function(h) !is.element(dates_matrix[,h], 
                                                   dates_sel))
  
  # MAE for selected dates
  mae_empirical_sel <- sapply(1:length(k_names), function(i) 
    sapply(1:H, function(h) 
      mean(abs(errors_list[[i]][,h][ind_dates[,h]]))))
  
  # errors
  sapply(1:H, function(h) mean(abs(dfm_mae[,h][ind_dates[,h]])))
  sapply(1:H, function(h) mean(abs(mlp_mae[,h][ind_dates[,h]])))
  sapply(1:H, function(h) mean(abs(arima_mae[,h][ind_dates[,h]])))
  sapply(1:H, function(h) mean(abs(arimax_mae[,h][ind_dates[,h]])))
  sapply(1:H, function(h) mean(abs(lasso_mae[,h][ind_dates[,h]])))
  sapply(1:H, function(h) mean(abs(favar_mae[,h][ind_dates[,h]])))
  mae_empirical_sel
  
  v_mae_dfm_sel <- sapply(1:H, function(h) mean(abs(dfm_mae[,h][ind_dates[,h]])))
  df_mae_dfm_sel <- v_mae_dfm_sel |>
    f_err_df("MAE", rhat, w)
  
  # RMSE for selected dates
  rmse_empirical_sel <- sapply(1:length(k_names), function(i) 
    sapply(1:H, function(h) 
      sqrt(mean((errors_list[[i]][,h][ind_dates[,h]])^2))))
  colnames(rmse_empirical_sel) <- colnames(mae_empirical_sel) <- k_names
  
  # errors
  sapply(1:H, function(h) sqrt(mean(dfm_mae[,h][ind_dates[,h]]^2)))
  sapply(1:H, function(h) sqrt(mean(mlp_mae[,h][ind_dates[,h]]^2)))
  sapply(1:H, function(h) sqrt(mean(arima_mae[,h][ind_dates[,h]]^2)))
  sapply(1:H, function(h) sqrt(mean(arimax_mae[,h][ind_dates[,h]]^2)))
  sapply(1:H, function(h) sqrt(mean(lasso_mae[,h][ind_dates[,h]]^2)))
  sapply(1:H, function(h) sqrt(mean(favar_mae[,h][ind_dates[,h]]^2)))
  rmse_empirical_sel
  
  v_rmse_dfm_sel <- sapply(1:H, function(h) sqrt(mean(dfm_mae[,h][ind_dates[,h]]^2)))
  df_rmse_dfm_sel <- v_rmse_dfm_sel |>
    f_err_df("RMSE", rhat, w)
  
  l_df_dfm_err <- list(df_mae_dfm, df_rmse_dfm, df_mae_dfm_sel, df_rmse_dfm_sel)
  
  tryCatch(
    {
      df_dfm_err_ <- readRDS(here::here("outputs/rds/models_objects/df_dfm_err.rds"))
      df_dfm_err <- map_df(l_df_dfm_err, ~ bind_rows(.x))
      df_dfm_err <- bind_rows(df_dfm_err_, df_dfm_err)
      if(nrow(df_dfm_err) <= (H * length(l_df_dfm_err) * length(rhats) * length(rolling_cases)))
      df_dfm_err |> 
        saveRDS(here::here("outputs/rds/models_objects/df_dfm_err.rds"))
    },
    error = function(e) {
      df_dfm_err <- map_df(l_df_dfm_err, ~ bind_rows(.x))
      df_dfm_err |>
        saveRDS(here::here("outputs/rds/models_objects/df_dfm_err.rds"))
    }
  )
  
  df_dfm_err_p <- l_df_dfm_err |>
    map_df(~ bind_rows(.x)) |>
    mutate(
      sample = c(rep("full", 4), rep("covid-19 excluded", 4))
    )
  
  cat("\n\nDFM erros under current scenario\n")
  
  print(df_dfm_err_p)
  
  cat("\n")

  
  ######################################## DM test & MCS procedure ##############################
  ###### Evaluating errors ######
  if(rhat == 1 & w == "expa") {
    
    set.seed(12345)
    #### DM test for MAE ####
    dm_empirical_mae <- dm_lasso_mae <- dm_favar_mae <- dm_mlp_mae <- 
      dm_arima_mae <- list()
    
    for(h in 1 : H){
      
      dm_reg_mae <- list()
      for(i in 1 : length(k_names)){
        dm_reg_mae[[i]] <- dm_test_hac(e1 = dfm_mae[,h],loss= "ae", 
                                       e2 = errors_list[[i]][,h], 
                                       h = h, alternative = "less", small.sample = TRUE)
      }
      names(dm_reg_mae) <- k_names
      
      dm_empirical_mae[[h]] <- dm_reg_mae
      
      dm_lasso_mae[[h]] <- dm_test_hac(e1 = dfm_mae[,h],loss= "ae", e2 = lasso_mae[,h], 
                                       h = h, alternative = "less", small.sample = TRUE)  
      
      dm_favar_mae[[h]] <- dm_test_hac(e1 = dfm_mae[,h],loss= "ae", e2 = favar_mae[,h],
                                       h = h, alternative = "less", small.sample = TRUE)

      dm_mlp_mae[[h]] <- dm_test_hac(e1 = dfm_mae[,h],loss= "ae", e2 = mlp_mae[,h], 
                                     h = h, alternative = "less", small.sample = TRUE)  
      
      dm_arima_mae[[h]] <- dm_test_hac(e1 = dfm_mae[,h],loss= "ae", e2 = arima_mae[,h], 
                                       h = h, alternative = "less", small.sample = TRUE)  
    }
    
    dm_mae <- matrix(0, length(k_names) + 4, 4)
    rownames(dm_mae) <- c(k_names, "LASSO", "FAVAR", "MLP", "ARIMA")
    
    colnames(dm_mae) <- c(paste("DM h = ", 1:H, sep = ""), 
                          paste( "p-value h = ", 1:H, sep = ""))
    
    
    dm_mae <- dm_mae[,c(1,3,2,4)]
    
    for(h in 1 : H){
      for(i in 1 : length(k_names)){
        dm_mae[k_names[i], paste("DM h = ", h, sep = "")] <- 
          dm_empirical_mae[[h]][[i]]$dm_hac
        
        dm_mae[k_names[i], paste("p-value h = ", h, sep = "")] <- 
          dm_empirical_mae[[h]][[i]]$p_value_asymptotic
        
      }
      
      dm_mae["LASSO", paste("DM h = ", h, sep = "")] <- dm_lasso_mae[[h]]$dm_hac
      dm_mae["FAVAR", paste("DM h = ", h, sep = "")] <- dm_favar_mae[[h]]$dm_hac
      dm_mae["MLP", paste("DM h = ", h, sep = "")] <- dm_mlp_mae[[h]]$dm_hac
      dm_mae["ARIMA", paste("DM h = ", h, sep = "")] <- dm_arima_mae[[h]]$dm_hac
      
      dm_mae["LASSO", paste("p-value h = ", h, sep = "")] <- dm_lasso_mae[[h]]$p_value_asymptotic
      dm_mae["FAVAR", paste("p-value h = ", h, sep = "")] <- dm_favar_mae[[h]]$p_value_asymptotic
      dm_mae["MLP", paste("p-value h = ", h, sep = "")] <- dm_mlp_mae[[h]]$p_value_asymptotic
      dm_mae["ARIMA", paste("p-value h = ", h, sep = "")] <- dm_arima_mae[[h]]$p_value_asymptotic
      
    }  
    
    ## DM test outputs
    df_dm_mae <- round(dm_mae, 4) |>
      as_tibble(rownames = "Short") |>
      left_join(
        {
          cat_dfm[,c("Short", "Description")] |>
            as_tibble()
        },
        by = join_by(Short)
      ) |> 
      mutate(
        Description = if_else(is.na(Description), Short, Description)
      ) |>
      select(-Short) |>
      relocate(Description) |>
      select(`Variable / Model` = Description, everything())
    
    df_dm_mae <- df_dm_mae |>
      slice(1:(nrow(df_dm_mae) - 4)) |>
      arrange(`DM h = 1`) |>
      bind_rows(
        df_dm_mae |>
          slice((nrow(df_dm_mae) - 3):nrow(df_dm_mae))
      ) 
    
    df_dm_mae|> 
      saveRDS(here::here("outputs/rds/models_objects/dm_mae.rds"))
    
    cat("\nTable 2: Diebold–Mariano tests with HAC correction: IOAE versus alternative predictors and models\n")
    
    df_dm_mae |> print(n = 200)
    
    cat("\n")
    
    #### SPA tests ####
    hz <- c(1, 2)
    
    errors_alt <- matrix(0, nrow(dfm_mae), nrow(dm_mae))
    colnames(errors_alt) <- rownames(dm_mae)
    
    # The function creates a data frame with the statistics of the tests 
    f_spa <- function(h) {
      for(i in 1 : length(errors_list))
        errors_alt[,i] <- errors_list[[colnames(errors_alt)[i]]][,h]
      
      errors_alt[,c("LASSO", "FAVAR", "MLP", "ARIMA")] <- cbind(lasso_mae[,h], 
                                                                favar_mae[,h], 
                                                                mlp_mae[,h],
                                                                arima_mae[,h])
      
      spa_mae <- spa_test(abs(dfm_mae[,h]), abs(errors_alt), verbose = FALSE)
      spa_rmse <- spa_test(dfm_mae[,h]^2, errors_alt^2, verbose = FALSE)
      
      tibble(
        Metric = c("MAE", "RMSE"),
        h = c(h, h),
        Statistic = c(spa_mae$statistic, spa_rmse$statistic),
        p_value = c(spa_mae$p_value, spa_rmse$p_value)
      )
    }
    
    df_spa_test <- map_df(hz, ~ f_spa(.x)) |>
      arrange(Metric) 
    
    df_spa_test|>
      saveRDS(here::here("outputs/rds/models_objects/spa_test.rds"))
    
    cat("\nSPA test using the IOAE as the benchmark\n")
    
    print(df_spa_test)
    
    cat("\n")
    
    #### MCS for MAE ####
    MCS_mae <- list()
    
    for(h in 1 : H){ 
      
      cat(str_glue("\nIn progress: {h} of {H} MCS procedures of MAE\r"))
      flush.console()
      
      mae_h <- abs(cbind(dfm_mae[,h], sapply(1:length(errors_list), 
                                             function(x) errors_list[[x]][,h]), 
                         lasso_mae[,h], favar_mae[,h], mlp_mae[,h], arima_mae[,h]))
      
      colnames(mae_h) <- c("DFM", k_names,
                           "LASSO", "FAVAR", "MLP", "ARIMA")
      
      MCS_mae[[h]] <- MCSprocedure(Loss=mae_h, alpha=0.1, B=5000, statistic='TR', 
                                   cl=NULL, verbose = FALSE)
    }
    
    mcs_mae <- matrix(NA, nrow(dm_mae) + 1, 4)
    rownames(mcs_mae) <- c("DFM", rownames(dm_mae))
    
    
    colnames(mcs_mae) <- c(paste("Rank h = ", 1:H, sep = ""), 
                           paste( "TR h = ", 1:H, sep = ""))
    
    mcs_mae <- mcs_mae[,c(1,3,2,4)]
    
    for(h in 1 : H){
      mcs_mae[rownames( MCS_mae[[h]]@show), c(paste("Rank h = ", h, sep = ""), 
                                              paste( "TR h = ", h, sep = ""))] <- MCS_mae[[h]]@show[,c("Rank_R", "v_R")]
      
    }
    
    # Eliminated variables by procedure
    # Generic variables mapping
    var_map <- as_tibble(cat_dfm[, c("Short", "Description")]) |>
      bind_rows(
        tibble(
          Short = c("DFM", "LASSO", "FAVAR", "MLP", "ARIMA"),
          Description = c("DFM", "LASSO", "FAVAR", "MLP", "ARIMA")
        )
      )
    
    # function for mapping variables
    f_v_Description <- function(Short) {
      var_map <- cat_dfm[, c("Short", "Description")]
      rbind(
        data.frame(
            Short = c("DFM", "LASSO", "FAVAR", "MLP", "ARIMA"),
            Description = c("DFM", "LASSO", "FAVAR", "MLP", "ARIMA")
        )
      )
      var_map[,"Description"][which(var_map[,"Short"] %in% Short)]
    }
    
    # List of eliminated variables by procedure
    l_eliminated <- map(hz, ~ {
      if(length(rownames(mcs_mae)) == length(MCS_mae[[.x]]@Info$model.names)) {
        c("None")
      } else{
        v_shorts <- rownames(mcs_mae)[!(rownames(mcs_mae) %in% MCS_mae[[.x]]@Info$model.names)]
        f_v_Description(v_shorts)  
      }
    })
    
    l_eliminated |> saveRDS("outputs/rds/models_objects/l_elim_mae.rds")
    
    ## Procedure outputs
    df_mcs_mae <- round(mcs_mae, 4) |>
      as_tibble(
        rownames = "Short" 
      ) |>
      left_join(
        var_map,
        by = join_by(Short)
      ) |>
      mutate(
        Description = if_else(is.na(Description), Short, Description)
      ) |>
      select(`Variable / Model` = Description, everything(), -Short) 
    
    df_mcs_mae|>
      saveRDS(here::here("outputs/rds/models_objects/mcs_mae.rds"))
    
    cat("\n\nTable 3: MCS results at the 10% significance level based on the range statistic TR\n")
    
    df_mcs_mae|> print(n = 200)
    
    df_mcs_mae_pv <- tibble(
      h = 1:2,
      p_value = sapply(1:H, function(h) MCS_mae[[h]]@Info$mcs_pvalue) 
    )
    
    df_mcs_mae_pv |>
      saveRDS(here::here("outputs/rds/models_objects/mcs_mae_pvalue.rds"))
    
    cat("\nMCS p-value\n")
    
    print(df_mcs_mae_pv)
    
    #### DM test for RMSE ####
    dm_empirical_mse <- dm_lasso_mse <- dm_favar_mse <- dm_mlp_mse <- 
      dm_arima_mse <- list()
    
    for(h in 1 : H){
      
      dm_reg_mse <- list()
      for(i in 1 : length(k_names)){
        dm_reg_mse[[i]] <- dm_test_hac(e1 = dfm_mae[,h],loss= "se", 
                                       e2 = errors_list[[i]][,h], 
                                       h = h, alternative = "less", small.sample = TRUE)
      }
      names(dm_reg_mse) <- k_names
      
      dm_empirical_mse[[h]] <- dm_reg_mse
      
      dm_lasso_mse[[h]] <- dm_test_hac(e1 = dfm_mae[,h],loss= "se", e2 = lasso_mae[,h], 
                                       h = h, alternative = "less", small.sample = TRUE)  
      
      dm_favar_mse[[h]] <- dm_test_hac(e1 = dfm_mae[,h],loss= "se", e2 = favar_mae[,h],
                                       h = h, alternative = "less", small.sample = TRUE)
      # dm_favar_mse[[h]] <- dm_test_hac(e1 = dfm_mae[,h],loss= "se", e2 = favar_mae_fc[,h], 
      #                                  h = h, alternative = "less", small.sample = TRUE)
      
      dm_mlp_mse[[h]] <- dm_test_hac(e1 = dfm_mae[,h],loss= "se", e2 = mlp_mae[,h], 
                                     h = h, alternative = "less", small.sample = TRUE)  
      
      dm_arima_mse[[h]] <- dm_test_hac(e1 = dfm_mae[,h],loss= "se", e2 = arima_mae[,h], 
                                       h = h, alternative = "less", small.sample = TRUE)  
    }
    
    dm_mse <- matrix(0, length(k_names) + 4, 4)
    rownames(dm_mse) <- c(k_names, "LASSO", "FAVAR", "MLP", "ARIMA")
    
    colnames(dm_mse) <- colnames(dm_mae)
    
    
    for(h in 1 : H){
      for(i in 1 : length(k_names)){
        dm_mse[k_names[i], paste("DM h = ", h, sep = "")] <- 
          dm_empirical_mse[[h]][[i]]$dm_hac
        
        dm_mse[k_names[i], paste("p-value h = ", h, sep = "")] <- 
          dm_empirical_mse[[h]][[i]]$p_value_asymptotic
        
      }
      
      dm_mse["LASSO", paste("DM h = ", h, sep = "")] <- dm_lasso_mse[[h]]$dm_hac
      dm_mse["FAVAR", paste("DM h = ", h, sep = "")] <- dm_favar_mse[[h]]$dm_hac
      dm_mse["MLP", paste("DM h = ", h, sep = "")] <- dm_mlp_mse[[h]]$dm_hac
      dm_mse["ARIMA", paste("DM h = ", h, sep = "")] <- dm_arima_mse[[h]]$dm_hac
      
      dm_mse["LASSO", paste("p-value h = ", h, sep = "")] <- dm_lasso_mse[[h]]$p_value_asymptotic
      dm_mse["FAVAR", paste("p-value h = ", h, sep = "")] <- dm_favar_mse[[h]]$p_value_asymptotic
      dm_mse["MLP", paste("p-value h = ", h, sep = "")] <- dm_mlp_mse[[h]]$p_value_asymptotic
      dm_mse["ARIMA", paste("p-value h = ", h, sep = "")] <- dm_arima_mse[[h]]$p_value_asymptotic
      
    }  
    
    df_dm_rmse <- round(dm_mse, 4) |>
      as_tibble(rownames = "Short") |>
      left_join(
        {
          cat_dfm[,c("Short", "Description")] |>
            as_tibble()
        },
        by = join_by(Short)
      ) |> 
      mutate(
        Description = if_else(is.na(Description), Short, Description)
      ) |>
      select(-Short) |>
      relocate(Description) |>
      select(`Variable / Model` = Description, everything())
    
    df_dm_rmse <- df_dm_rmse |>
      slice(1:(nrow(df_dm_rmse) - 4)) |>
      arrange(`DM h = 1`) |>
      bind_rows(
        df_dm_rmse |>
          slice((nrow(df_dm_rmse) - 3):nrow(df_dm_rmse))
      )
    
    df_dm_rmse|> 
      saveRDS(here::here("outputs/rds/models_objects/dm_rmse.rds"))
    
    cat("\nDM-HAC test, RMSE as the loss function\n")
    
    df_dm_rmse|> print(n = 200)
    
    cat("\n")
    
    #### MCS for RMSE ####
    MCS_mse <- list()
    for(h in 1 : H){ 
      
      cat(str_glue("\nExecuting {h} of {H} MCS procedures of RMSE\r"))
      flush.console()

      mse_h <- cbind(dfm_mae[,h], sapply(1:length(errors_list), 
                                         function(x) errors_list[[x]][,h]), lasso_mae[,h], 
                     favar_mae[,h], mlp_mae[,h], arima_mae[,h])^2
      
      colnames(mse_h) <-colnames(mae_h)
      
      MCS_mse[[h]] <- MCSprocedure(Loss=mse_h, alpha=0.1, B=5000, statistic='TR', 
                                   cl=NULL, verbose = FALSE)
    }
    
    mcs_mse <- matrix(NA, nrow(dm_mse) + 1, 4)
    rownames(mcs_mse) <- rownames(mcs_mae)
    
    colnames(mcs_mse) <- c(paste("Rank h = ", 1:H, sep = ""), 
                           paste( "TR h = ", 1:H, sep = ""))
    
    mcs_mse <- mcs_mse[,c(1,3,2,4)]
    
    for(h in 1 : H){
      mcs_mse[rownames(MCS_mse[[h]]@show), c(paste("Rank h = ", h, sep = ""), 
                                             paste( "TR h = ", h, sep = ""))] <- MCS_mse[[h]]@show[,c("Rank_R", "v_R")]
      
    }
    
    # List of eliminated variables by procedure
    l_eliminated <- map(hz, ~ {
      if(length(rownames(mcs_mse)) == length(MCS_mse[[.x]]@Info$model.names)) {
        c("None")
      } else{
        v_shorts <- rownames(mcs_mse)[!(rownames(mcs_mse) %in% MCS_mse[[.x]]@Info$model.names)]
        f_v_Description(v_shorts)  
      }
    })
    
    l_eliminated |> saveRDS("outputs/rds/models_objects/l_elim_rmse.rds")
    
    ## Procedure outputs
    df_mcs_rmse <- round(mcs_mse, 4) |>
      as_tibble(
        rownames = "Short" 
      ) |>
      left_join(
        var_map,
        by = join_by(Short)
      ) |>
      mutate(
        Description = if_else(is.na(Description), Short, Description)
      ) |>
      select(`Variable / Model` = Description, everything(), -Short)
    
    df_mcs_rmse |>
      saveRDS(here::here("outputs/rds/models_objects/mcs_rmse.rds"))
    
    cat("\n\nMCS results, RMSE as the loss function\n")
    
    df_mcs_rmse |> print(n = 200)
    
    df_mcs_rmse_pv <- tibble(
      h = 1:2,
      p_value = sapply(1:H, function(h) MCS_mse[[h]]@Info$mcs_pvalue) 
    ) 
    
    df_mcs_rmse_pv |>
      saveRDS(here::here("outputs/rds/models_objects/mcs_rmse_pvalue.rds"))
    
    cat("\nMCS p-value\n")
    
    print(df_mcs_rmse_pv)
    
    cat("\n")
  }
  
  
  
  ######################################## getting figures inputs ##############################
  
  if(rhat == 1 & w == "expa") {
    
    #### loading factor
    j <- 1
    Phat_t <- t(sapply(1:Ht, function(x) loadings_dfm[[x]][1:ncol(db_dfm_p), j]))
    rownames(Phat_t) <- names(yfcst_levels)[out_sample[,1]-1]
    
    phat_names <- rep(NA, Ht)
    phat_names[seq(1, Ht, 3)] <- rownames(Phat_t)[seq(1, Ht, 3)]
    
    variables <- cat_dfm[,"Description"]
    colnames(Phat_t) <- variables
    
    ord_phat <- order(abs(Phat_t[Ht,]), decreasing = TRUE)
    Phat_t <- Phat_t[,ord_phat]
    
    
    #### input figure 7
    types <- c("annual", "monthly", "levels")
    
    for(type in types) {
      for(x in 1:Ht) {
        if(x == 1 & type == "annual"){ df_mae <- as_tibble(nowcasts_models[[x]][type][[1]]) |>
          mutate(horizon = 1:2,
                 type = type,
                 period = row.names(nowcasts_models[[x]][type][[1]]))
        } else {
          df_mae <- df_mae |>
            bind_rows(
              as_tibble(nowcasts_models[[x]][type][[1]]) |>
                mutate(horizon = 1:2,
                       type = type,
                       period = row.names(nowcasts_models[[x]][type][[1]]))
            )
        }
      }
    }
    
    df_mae |> saveRDS(here::here("outputs/rds/figures_inputs/f7.rds"))
    
    
    #### input figure 6
    fhat <- -nowcast_object$factors$dynamic_factors

    #tibble(fhat_1) |>
    tibble(fhat) |>
      mutate(period = time(fhat),
             period = as.Date(period)) |>
      bind_cols(
        tibble(nowcast_object$yfcst[,"y"])
      ) |>
      #rename(igae = 3, fhat = fhat_1) |>
      rename(igae = 3, fhat = fhat) |>
      pivot_longer(cols = c(fhat, igae)) |> 
      saveRDS(here::here("outputs/rds/figures_inputs/f6.rds"))
  
  
    #### input figure 5
    confidence$Phat[[1]] |>
      cbind(TS |> dplyr::select(variable = Names)) |>
      as_tibble() |> 
      saveRDS(here::here("outputs/rds/figures_inputs/f5.rds"))
  
    #### input appendix figures
    # MAE cum
    for(h in 1 : H){
      mat_lasso <- ts(apply(abs(cbind(dfm_mae[,h], lasso_mae[,h])), 2, cumsum)/(1:Ht),
                      start = c(2018, 9+h), frequency = 12)
      
      saveRDS(mat_lasso, here::here(paste0("outputs/rds/figures_inputs/mat_mae_lasso_1",h,".rds")))
      
    }  
    
    for(h in 1 : H){
      mat_favar <- ts(apply(abs(cbind(dfm_mae[,h], favar_mae[,h])), 2, cumsum)/(1:Ht),
                      start = c(2018, 9+h), frequency = 12)
      
      saveRDS(mat_favar, here::here(paste0("outputs/rds/figures_inputs/mat_mae_favar_1",h,".rds")))
      
    }
    
    for(h in 1 : H){
      mat_mlp <- ts(apply(abs(cbind(dfm_mae[,h], mlp_mae[,h])), 2, cumsum)/(1:Ht),
                    start = c(2018, 9+h), frequency = 12)
      
      saveRDS(mat_mlp, here::here(paste0("outputs/rds/figures_inputs/mat_mae_mlp_1",h,".rds")))
      
    }
    
    for(h in 1 : H){
      mat_arima <- ts(apply(abs(cbind(dfm_mae[,h], arima_mae[,h])), 2, cumsum)/(1:Ht),
                      start = c(2018, 9+h), frequency = 12)
      
      saveRDS(mat_arima, here::here(paste0("outputs/rds/figures_inputs/mat_mae_arima_1",h,".rds")))
      
    }
    
    for(h in 1 : H){ 
      
      mat_lasso <- matrix(NA, Ht, 2)
      mat_lasso[ind_dates[,h], 1] <- cumsum(abs(dfm_mae[ind_dates[,h],h]))/
        (1:(Ht-sum(!ind_dates[,h])))
      mat_lasso[ind_dates[,h], 2] <- cumsum(abs(lasso_mae[ind_dates[,h],h]))/
        (1:(Ht-sum(!ind_dates[,h])))
      mat_lasso <- ts(mat_lasso, start = c(2018, 9+h), frequency = 12)
      
      saveRDS(mat_lasso, here::here(paste0("outputs/rds/figures_inputs/mat_mae_lasso_2",h,".rds")))
    }  
    
    for(h in 1 : H){
      mat_favar <- matrix(NA, Ht, 2)
      mat_favar[ind_dates[,h], 1] <- cumsum(abs(dfm_mae[ind_dates[,h],h]))/
        (1:(Ht-sum(!ind_dates[,h])))
      mat_favar[ind_dates[,h], 2] <- cumsum(abs(favar_mae[ind_dates[,h],h]))/
        (1:(Ht-sum(!ind_dates[,h])))
      mat_favar <- ts(mat_favar, start = c(2018, 9+h), frequency = 12)
      
      saveRDS(mat_favar, here::here(paste0("outputs/rds/figures_inputs/mat_mae_favar_2",h,".rds")))
    }
    
    for(h in 1 : H){
      mat_mlp <- matrix(NA, Ht, 2)
      mat_mlp[ind_dates[,h], 1] <- cumsum(abs(dfm_mae[ind_dates[,h],h]))/
        (1:(Ht-sum(!ind_dates[,h])))
      mat_mlp[ind_dates[,h], 2] <- cumsum(abs(mlp_mae[ind_dates[,h],h]))/
        (1:(Ht-sum(!ind_dates[,h])))
      mat_mlp <- ts(mat_mlp, start = c(2018, 9+h), frequency = 12)
      
      saveRDS(mat_mlp, here::here(paste0("outputs/rds/figures_inputs/mat_mae_mlp_2",h,".rds")))
    }
    
    for(h in 1 : H){
      mat_arima <- matrix(NA, Ht, 2)
      mat_arima[ind_dates[,h], 1] <- cumsum(abs(dfm_mae[ind_dates[,h],h]))/
        (1:(Ht-sum(!ind_dates[,h])))
      mat_arima[ind_dates[,h], 2] <- cumsum(abs(arima_mae[ind_dates[,h],h]))/
        (1:(Ht-sum(!ind_dates[,h])))
      mat_arima <- ts(mat_arima, start = c(2018, 9+h), frequency = 12)
      
      saveRDS(mat_arima, here::here(paste0("outputs/rds/figures_inputs/mat_mae_arima_2",h,".rds")))
    }
    
    
    # RMSE cum
    for(h in 1 : H){
      mat_lasso <- sqrt(ts(apply(cbind(dfm_mae[,h], lasso_mae[,h])^2, 2, cumsum)/(1:Ht),
                           start = c(2018, 9+h), frequency = 12))
      
      saveRDS(mat_lasso, here::here(paste0("outputs/rds/figures_inputs/mat_rmse_lasso_1",h,".rds")))
    }  
    
    for(h in 1 : H){
      mat_favar <- sqrt(ts(apply(cbind(dfm_mae[,h], favar_mae[,h])^2, 2, cumsum)/(1:Ht),
                           start = c(2018, 9+h), frequency = 12))
      
      saveRDS(mat_favar, here::here(paste0("outputs/rds/figures_inputs/mat_rmse_favar_1",h,".rds")))
    }
    
    for(h in 1 : H){
      mat_mlp <- sqrt(ts(apply(cbind(dfm_mae[,h], mlp_mae[,h])^2, 2, cumsum)/(1:Ht),
                         start = c(2018, 9+h), frequency = 12))
      
      saveRDS(mat_mlp, here::here(paste0("outputs/rds/figures_inputs/mat_rmse_mlp_1",h,".rds")))
    }
    for(h in 1 : H){
      mat_arima <- sqrt(ts(apply(cbind(dfm_mae[,h], arima_mae[,h])^2, 2, cumsum)/(1:Ht),
                           start = c(2018, 9+h), frequency = 12))
      
      saveRDS(mat_arima, here::here(paste0("outputs/rds/figures_inputs/mat_rmse_arima_1",h,".rds")))
    }
    
    for(h in 1 : H){ 
      mat_lasso <- matrix(NA, Ht, 2)
      mat_lasso[ind_dates[,h], 1] <- sqrt(cumsum(dfm_mae[ind_dates[,h],h]^2)/
                                            (1:(Ht-sum(!ind_dates[,h]))))
      mat_lasso[ind_dates[,h], 2] <- sqrt(cumsum(lasso_mae[ind_dates[,h],h]^2)/
                                            (1:(Ht-sum(!ind_dates[,h]))))
      mat_lasso <- ts(mat_lasso, start = c(2018, 9+h), frequency = 12)
      
      saveRDS(mat_lasso, here::here(paste0("outputs/rds/figures_inputs/mat_rmse_lasso_2",h,".rds")))
    }  
    
    for(h in 1 : H){
      mat_favar <- matrix(NA, Ht, 2)
      mat_favar[ind_dates[,h], 1] <- sqrt(cumsum(dfm_mae[ind_dates[,h],h]^2)/
                                            (1:(Ht-sum(!ind_dates[,h]))))
      mat_favar[ind_dates[,h], 2] <- sqrt(cumsum(favar_mae[ind_dates[,h],h]^2)/
                                            (1:(Ht-sum(!ind_dates[,h]))))
      mat_favar <- ts(mat_favar, start = c(2018, 9+h), frequency = 12)
      
      saveRDS(mat_favar, here::here(paste0("outputs/rds/figures_inputs/mat_rmse_favar_2",h,".rds")))
    }
    
    for(h in 1 : H){
      mat_mlp <- matrix(NA, Ht, 2)
      mat_mlp[ind_dates[,h], 1] <- sqrt(cumsum(dfm_mae[ind_dates[,h],h]^2)/
                                          (1:(Ht-sum(!ind_dates[,h]))))
      mat_mlp[ind_dates[,h], 2] <- sqrt(cumsum(mlp_mae[ind_dates[,h],h]^2)/
                                          (1:(Ht-sum(!ind_dates[,h]))))
      mat_mlp <- ts(mat_mlp, start = c(2018, 9+h), frequency = 12)
      
      saveRDS(mat_mlp, here::here(paste0("outputs/rds/figures_inputs/mat_rmse_mlp_2",h,".rds")))
    }
    for(h in 1 : H){
      mat_arima <- matrix(NA, Ht, 2)
      mat_arima[ind_dates[,h], 1] <- sqrt(cumsum(dfm_mae[ind_dates[,h],h]^2)/
                                            (1:(Ht-sum(!ind_dates[,h]))))
      mat_arima[ind_dates[,h], 2] <- sqrt(cumsum(arima_mae[ind_dates[,h],h]^2)/
                                            (1:(Ht-sum(!ind_dates[,h]))))
      mat_arima <- ts(mat_arima, start = c(2018, 9+h), frequency = 12)
      
      saveRDS(mat_arima, here::here(paste0("outputs/rds/figures_inputs/mat_rmse_arima_2",h,".rds")))
    }
    
  }
