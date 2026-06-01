# -------------------------------
# adf test
# -------------------------------
adf <- function(x, type = c("none", "const", "trend"),
                alternative = c("stationary", "explosive"),
                k = trunc((length(x) - 1)^(1/3))){
  
  if (NCOL(x) > 1)
    stop("x is not a vector or univariate time series")
  if (any(is.na(x)))
    stop("NAs in x")
  if (k < 0)
    stop("k negative")
  
  type <- match.arg(type)
  alternative <- match.arg(alternative)
  DNAME <- deparse(substitute(x))
  
  k <- k + 1
  x <- as.vector(x, mode = "double")
  y <- diff(x)
  n <- length(y)
  z <- embed(y, k)
  yt <- z[, 1]
  xt1 <- x[k:n]
  tt <- k:n
  
  if(type == "none"){
    if(k > 1){
      yt1 <- z[, 2:k]
      res <- lm(yt ~ xt1 + yt1 - 1)
    }else{
      res <- lm(yt ~ xt1 - 1)
    }
    res.sum <- summary(res)
    STAT <- res.sum$coefficients[1, 1]/res.sum$coefficients[1,
                                                            2]
    Table <- adfTable(trend = "nc", statistic = "t")
  }
  
  if(type == "const"){
    if(k > 1){
      yt1 <- z[, 2:k]
      res <- lm(yt ~ xt1 + 1 + yt1)
    }else{
      res <- lm(yt ~ xt1 + 1)
    }
    res.sum <- summary(res)
    STAT <- res.sum$coefficients[2, 1]/res.sum$coefficients[2,
                                                            2]
    Table <- adfTable(trend = "c", statistic = "t")
  }
  
  if(type == "trend"){
    if (k > 1) {
      yt1 <- z[, 2:k]
      res <- lm(yt ~ xt1 + 1 + tt + yt1)
    }else{
      res <- lm(yt ~ xt1 + 1 + tt)
    }
    res.sum <- summary(res)
    STAT <- res.sum$coefficients[2, 1]/res.sum$coefficients[2,
                                                            2]
    Table <- adfTable(trend = "ct", statistic = "t")
  }
  bic <- BIC(res)
  table <- Table$z
  tablen <- dim(table)[2]
  tableT <- Table$x
  tablep <- Table$y
  tableipl <- numeric(tablen)
  for (i in (1:tablen)) tableipl[i] <- approx(tableT, table[,
                                                            i], n, rule = 2)$y
  interpol <- approx(tableipl, tablep, STAT, rule = 2)$y
  if (is.na(approx(tableipl, tablep, STAT, rule = 1)$y))
    if (interpol == min(tablep))
      warning("p-value smaller than printed p-value")
  else warning("p-value greater than printed p-value")
  if (alternative == "stationary")
    PVAL <- interpol
  else if (alternative == "explosive")
    PVAL <- 1 - interpol
  else stop("irregular alternative")
  PARAMETER <- k - 1
  METHOD <- "Augmented Dickey-Fuller Test"
  names(STAT) <- "Dickey-Fuller"
  names(PARAMETER) <- "Lag order"
  structure(list(bic = bic, statistic = STAT, parameter = PARAMETER,
                 alternative = alternative, p.value = PVAL, method = METHOD,
                 data.name = DNAME), class = "htest")
}


# -------------------------------
# adf tables
# -------------------------------
adfTable <- function(trend = c("nc", "c", "ct"), statistic = c("t", "n"),
                     includeInf = TRUE)
{
  # A function implemented by Diethelm Wuertz
  
  # Description:
  #   Tables critical values for augmented Dickey-Fuller test.
  
  # Note:
  #   x=-3:0; y=0:3; z=outer(x,y,"*"); rownames(z)=x; colnames(z)=y; z
  
  # Examples:
  #   adfTable()
  
  # FUNCTION:
  
  # Match Arguments:
  type = trend = match.arg(trend)
  statistic = match.arg(statistic)
  
  # Tables:
  if (statistic == "t") {
    # Hamilton Table B.6 - OLS t-Statistic
    if (type == "nc") {
      table = cbind(
        c(-2.66, -2.26, -1.95, -1.60, +0.92, +1.33, +1.70, +2.16),
        c(-2.62, -2.25, -1.95, -1.61, +0.91, +1.31, +1.66, +2.08),
        c(-2.60, -2.24, -1.95, -1.61, +0.90, +1.29, +1.64, +2.03),
        c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.29, +1.63, +2.01),
        c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.28, +1.62, +2.00),
        c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.28, +1.62, +2.00))
    } else if (type == "c") {
      table = cbind(
        c(-3.75, -3.33, -3.00, -2.63, -0.37, +0.00, +0.34, +0.72),
        c(-3.58, -3.22, -2.93, -2.60, -0.40, -0.03, +0.29, +0.66),
        c(-3.51, -3.17, -2.89, -2.58, -0.42, -0.05, +0.26, +0.63),
        c(-3.46, -3.14, -2.88, -2.57, -0.42, -0.06, +0.24, +0.62),
        c(-3.44, -3.13, -2.87, -2.57, -0.43, -0.07, +0.24, +0.61),
        c(-3.43, -3.12, -2.86, -2.57, -0.44, -0.07, +0.23, +0.60))
    } else if (type == "ct") {
      table = cbind(
        c(-4.38, -3.95, -3.60, -3.24, -1.14, -0.80, -0.50, -0.15),
        c(-4.15, -3.80, -3.50, -3.18, -1.19, -0.87, -0.58, -0.24),
        c(-4.04, -3.73, -3.45, -3.15, -1.22, -0.90, -0.62, -0.28),
        c(-3.99, -3.69, -3.43, -3.13, -1.23, -0.92, -0.64, -0.31),
        c(-3.98, -3.68, -3.42, -3.13, -1.24, -0.93, -0.65, -0.32),
        c(-3.96, -3.66, -3.41, -3.12, -1.25, -0.94, -0.66, -0.33))
    } else {
      stop("Invalid type specified")
    }
  } else if (statistic == "z" || statistic == "n") {
    # Hamilton Table B.5 - Based on OLS Autoregressive Coefficient
    if (type == "nc") {
      table = cbind(
        c(-11.9,  -9.3,  -7.3,  -5.3, +1.01, +1.40, +1.79, +2.28),
        c(-12.9,  -9.9,  -7.7,  -5.5, +0.97, +1.35, +1.70, +2.16),
        c(-13.3, -10.2,  -7.9,  -5.6, +0.95, +1.31, +1.65, +2.09),
        c(-13.6, -10.3,  -8.0,  -5.7, +0.93, +1.28, +1.62, +2.04),
        c(-13.7, -10.4,  -8.0,  -5.7, +0.93, +1.28, +1.61, +2.04),
        c(-13.8, -10.5,  -8.1,  -5.7, +0.93, +1.28, +1.60, +2.03))
    } else if (type == "c") {
      table = cbind(
        c(-17.2, -14.6, -12.5, -10.2, -0.76, +0.01, +0.65, +1.40),
        c(-18.9, -15.7, -13.3, -10.7, -0.81, -0.07, +0.53, +1.22),
        c(-19.8, -16.3, -13.7, -11.0, -0.83, -0.10, +0.47, +1.14),
        c(-20.3, -16.6, -14.0, -11.2, -0.84, -0.12, +0.43, +1.09),
        c(-20.5, -16.8, -14.0, -11.2, -0.84, -0.13, +0.42, +1.06),
        c(-20.7, -16.9, -14.1, -11.3, -0.85, -0.13, +0.41, +1.04))
    } else if (type == "ct") {
      table = cbind(
        c(-22.5, -19.9, -17.9, -15.6, -3.66, -2.51, -1.53, -0.43),
        c(-25.7, -22.4, -19.8, -16.8, -3.71, -2.60, -1.66, -0.65),
        c(-27.4, -23.6, -20.7, -17.5, -3.74, -2.62, -1.73, -0.75),
        c(-28.4, -24.4, -21.3, -18.0, -3.75, -2.64, -1.78, -0.82),
        c(-28.9, -24.8, -21.5, -18.1, -3.76, -2.65, -1.78, -0.84),
        c(-29.5, -25.1, -21.8, -18.3, -3.77, -2.66, -1.79, -0.87))
    } else {
      stop("Invalid type specified")
    }
  } else {
    stop("Invalid statistic specified")
  }
  
  # Transpose:
  Table = t(table)
  colnames(Table) = c("0.010", "0.025", "0.050", "0.100", "0.900",
                      "0.950", "0.975", "0.990")
  rownames(Table) = c(" 25", " 50", "100", "250", "500", "Inf")
  ans = list(
    x = as.numeric(rownames(Table)),
    y = as.numeric(colnames(Table)),
    z = Table)
  class(ans) = "gridData"
  
  # Exclude Inf:
  if (!includeInf) {
    nX = length(ans$x)
    ans$x = ans$x[-nX]
    ans$z = ans$z[-nX, ]
  }
  
  # Add Control:
  attr(ans, "control") <-
    c(table = "adf", trend = trend, statistic = statistic)
  
  # Return Value:
  ans
}


# -------------------------------
# @ pooled test
# e is the idiosyncratic errors
# type is the specification in the unit root test
# lag_max is the maximum number of lags to testing
# -------------------------------
pooled.test <- function(e, type = "const", lag_max = 7){
  N <- ncol(e)
  bics <- matrix(0, lag_max + 1, N)
  for(i in 1 : N){
    bics[,i] <- sapply(0:lag_max, function(x) adf(e[,i], type = type, 
                                                  k = x)$bic)
  }
  lags <- apply(bics, 2, function(x) which(x == min(x)) - 1)
  adf_ehat <- sapply(1:N, function(x) adf(e[,x], type = type,
                                          k = lags[x])$p.value)
  P <- (-2*sum(log(adf_ehat)) - 2*N)/sqrt(4*N)
  p.value <- dnorm(P)
  #p.value <- 1 - pnorm(P)
  stats <- round(c(P, p.value), 4)
  names(stats) <- c("P", "p.value")
  result <- list(adf_ehat, stats)
  names(result) <- c("ADF", "stats")
  return(result)
}

  
# -------------------------------
# cambios porcentuales
# -------------------------------
cpm <- function(x, q = 1){

  index <- matrix(0, length(x) - q, 2)
  index[,1] <- 1:nrow(index)
  index[,2] <- (q+1):(nrow(index)+q)

  xm <- rep(0, nrow(index))
  for(i in 1 : length(xm))
    xm[i] <- x[index[i,2]]/x[index[i,1]]*100-100

  return(xm)
}


# -------------------------------
# transformation
# -------------------------------
best_trans <- function(y, x, names_y, names_x, rel){
  
  x_trans <- rep(NA, length(names_x))
  names(x_trans) <- names_x
  
  l <- x
  names(l) <- names_x
  
  mv <- cpm(x)
  names(mv) <- names_x[-1]
    
  av <- cpm(x, 12)
  names(av) <- names_x[-(1:12)]
  
  rho <- c(cor(y, l[names_y], use = "pairwise.complete.obs"),
          cor(y, mv[names_y], use = "pairwise.complete.obs"),
          cor(y, av[names_y], use = "pairwise.complete.obs"))
  names(rho) <- c("linear", "monthly", "annual")
  
  if(rel == "I")
    max_rho <- names(rho)[which(abs(rho) == max(abs(rho)))]
  
   if(rel == "P")
    max_rho <- names(rho)[which(rho == max(rho))]
  
  if(rel == "N")
    max_rho <- names(rho)[which(rho == min(rho))]
 
  
  if(max_rho == "linear")
    x_trans[names(l)] <- l
  
  if(max_rho == "monthly")
    x_trans[names(mv)] <- mv  
  
  if(max_rho == "annual")
    x_trans[names(av)] <- av
  
  result <- list(rho, max_rho, x_trans)
  names(result) <- c("rho", "trans", "x_trans")
  
  return(result)
}


# -------------------------------
# funcion para obtener nowcasts alternativos
# -------------------------------
lin_function <- function(y, x, d = 12){
  
  # y <- db_yfcst_raw[,vari] 
  # x <-  now_test$forecasts_dfm
  # d <- 1
  
  y <- y[!is.na(y)]
  n <- length(y)
  
  H <- nrow(x)
  
  now_lev <- now_mv <- now_av <- matrix(0, H, 3)
  colnames(now_lev) <- colnames(now_mv) <- colnames(now_av) <- colnames(x)
  
  # niveles
  if(is.null(d)){
    now_lev <- x
    y <- c(y, now_lev[,"mean"])
    now_mv[,"mean"] <- cpm(c(y, now_lev[,"mean"]), 1)[(n):(n+H-1)]
    
    #for(i in 1 : H){
      init <- y[n]
      
      for(i in 1 : H){
        now_mv[i, "lower"] <- cpm(c(init, now_lev[i,"lower"]), 1)
        now_mv[i, "upper"] <- cpm(c(init, now_lev[i,"upper"]), 1)
        
        init <- now_lev[i,"mean"]
      }
    #}
    
    for(i in 1 : H){
      init <- y[n-12+i]
  
      now_av[i, "mean"] <- cpm(c(init, now_lev[i,"mean"]), 1)
      now_av[i, "lower"] <- cpm(c(init, now_lev[i,"lower"]), 1)
      now_av[i, "upper"] <- cpm(c(init, now_lev[i,"upper"]), 1)
    }
  }else{
    
    # variaciones mensuales
    if(d == 1){
      now_mv <- x
      
      #for(i in 1 : H){
        init <- y[n]
        
        for(i in 1 : H){ 
          now_lev[i, "mean"] <- init*(1 + now_mv[i, "mean"]/100)
          now_lev[i, "lower"] <- init*(1 + now_mv[i, "lower"]/100)
          now_lev[i, "upper"] <- init*(1 + now_mv[i, "upper"]/100)
          
          init <- now_lev[i,"mean"]
          
          y[n+i] <- init
        }
      #}
      
      for(i in 1 : H){
        init <- y[n-12+i]
    
        now_av[i, "mean"] <- cpm(c(init, now_lev[i,"mean"]), 1)
        now_av[i, "lower"] <- cpm(c(init, now_lev[i,"lower"]), 1)
        now_av[i, "upper"] <- cpm(c(init, now_lev[i,"upper"]), 1)
        
      }
      
    }
    
    # variaciones anuales
    if(d == 12){
      now_av <- x
      
      for(i in 1 : H){ 
        init <- y[n-12+i]
  
        now_lev[i, "mean"] <- init*(1 + now_av[i, "mean"]/100)
        now_lev[i, "lower"] <- init*(1 + now_av[i, "lower"]/100)
        now_lev[i, "upper"] <- init*(1 + now_av[i, "upper"]/100)
        
        y[n+i] <- now_lev[i, "mean"]
      }
      
      now_mv[,"mean"] <- cpm(c(y, now_lev[,"mean"]), 1)[(n):(n+H-1)]
      #for(i in 1 : H){
        init <- y[n]
        
        for(i in 1 : H){
          now_mv[i, "lower"] <- cpm(c(init, now_lev[i,"lower"]), 1)
          now_mv[i, "upper"] <- cpm(c(init, now_lev[i,"upper"]), 1)
          
          init <- now_lev[i,"mean"]
        }
      #}
    }
  }
  
  rownames(now_lev) <- rownames(now_mv) <- rownames(now_av) <- rownames(x)
      
  result <- list(now_lev, now_mv, now_av)
  names(result) <- c("levels", "monthly", "annual")
  
  return(result)
}


# -------------------------------
# intervalos de confianza en el dfm
# -------------------------------
confidence_dfm <- function(dfm, alpha = 0.05, rhat = 1, qs = 2){
  
  # matrices
  y <- dfm$yfcst[,"y"]
  Ty <- min(which(is.na(y))) - 1
  X_s <- scale(dfm$xfcst[1:Ty,])
  
  y <- dfm$yfcst[,"y"]
  
  # number of variables
  N <- ncol(X_s)
  
  # @ PC factor extraction
  
  # eigen decomposition
  ed <- eigen(t(X_s)%*%X_s)
  
  # eigenvectors
  L <- ed$vectors[,1:rhat,drop=FALSE]
  
  # factors and loadings
  Phat <- sqrt(N)*L 
  fhat <- X_s%*%Phat/N
  
  rownames(Phat) <- colnames(X_s)
  colnames(Phat) <- colnames(fhat) <- paste("f", 1:rhat, sep ="")
  
  # idiosyncratic noises
  chi <- fhat%*%t(Phat)
  ehat <- X_s - chi
  
  # rotating the factors and loadings
  rho <- cor(fhat[,1], y[1:Ty], use = "complete.obs")
  
  if(rho < 0){
    fhat <- -fhat
    Phat <- -Phat
    rho <- -rho
  }
  
  # @ static factor confidence interval 
  cil <- ciu <- St <- matrix(0, Ty, rhat)
  
  # estimates
  z <- qnorm(1-alpha/2, 0, 1)
  er <- eigen(t(X_s)%*%X_s/(N*Ty))
  V <- diag(er$values[1:rhat],rhat,rhat)
  Vinv <- solve(V)
  
  # asymptotic variance for factors
  for (t in 1:Ty) {
    Gamma <- matrix(0, rhat, rhat)
    
    for (i in 1:N) {
      Gamma.i <- matrix(Phat[i,])%*%t(matrix(Phat[i,]))*(ehat[t,i]^2)
      Gamma <- Gamma + Gamma.i
    }
    Gamma <- (1/N)*(Gamma)
    St[t,] <- sqrt(diag(Vinv%*%Gamma%*%Vinv))
    
    # confidence Intervals
    cil[t,] <- fhat[t,] - z*N^(-1/2)*St[t,]
    ciu[t,] <- fhat[t,] + z*N^(-1/2)*St[t,]
  }
  
  # asymptotic variance for loadings
  Phi_list <- list()
  for(i in 1 : N){
    Phi <- matrix(0, rhat, rhat)
    for (t in 1:Ty) {
      Phi.i <- matrix(fhat[t,])%*%t(matrix(fhat[t,]))*(ehat[t,i]^2)
      Phi <- Phi + Phi.i
    }
    Phi <- Phi/Ty
    Phi_list[[i]] <- Phi
  }
  
  Divq <- list()
  for(v in 1 : qs){
    Div_list <- list()
    for(i in 1 : N){
      Div <- matrix(0, rhat, rhat)
      for(t in (v+1):Ty){
        Div.i <- fhat[t,]*fhat[t-v,]*ehat[t,i]*ehat[t-v,i]
        Div <- Div + Div.i
      }
      Div <- (1-v/(qs+1))*(Div + t(Div))
      Div_list[[i]] <- Div/Ty
    }
    Divq[[v]] <- Div_list
  }
  
  D_list <- list()
  for(i in 1 : N)
    D_list[[i]] <- Divq[[1]][[i]] + Divq[[2]][[i]]
  
  se <- matrix(0, N, rhat)
  for(i in 1 : N){
    se[i,] <- sqrt(diag((Phi_list[[i]] + D_list[[i]])/Ty))
  }
  
  # intervals
  Fhat_list <- Phat_list <- list()
  
  for(i in 1 : rhat){
    f_mat <- cbind(cil[,i], fhat[,i], ciu[,i])
    colnames(f_mat) <- c("lower", "fhat", "upper")
    
    Fhat_list[[i]] <- f_mat
    
    li <- Phat[,i] - z*se[,i]
    ls <- Phat[,i] + z*se[,i]
    
    P_mat <- cbind(li, Phat[,i], ls)
    colnames(P_mat) <- c("lower", "Phat", "upper")
    
    Phat_list[[i]] <- P_mat
  }
  
  result <- list(Fhat_list, Phat_list)
  names(result) <- c("Fhat", "Phat")
  
  return(result)
}


# -------------------------------
# Diebold-Mariano con HAC, bootstrap y small-sample correction (HLN, 1997)
# -------------------------------
dm_test_hac <- function(e1, e2, h = 1,
                        loss = c("se","ae","custom"),
                        d_custom = NULL,
                        alternative = c("two.sided","less","greater"),
                        lag = NULL,              # si NULL usa max(0, h-1)
                        prewhite = FALSE,        # pre-blanqueo en NeweyWest
                        bootstrap = FALSE,
                        B = 2000,
                        block_length = NULL,     # si NULL usa ceiling(sqrt(n))
                        seed = NULL,
                        small.sample = FALSE) {  # NUEVO argumento
  # Dependencias
  if (!requireNamespace("sandwich", quietly = TRUE)) {
    stop("Necesitas el paquete 'sandwich'. install.packages('sandwich')")
  }
  if (!requireNamespace("lmtest", quietly = TRUE)) {
    stop("Necesitas el paquete 'lmtest'. install.packages('lmtest')")
  }
  
  loss <- match.arg(loss)
  alternative <- match.arg(alternative)
  
  # Alinear y limpiar NAs
  e1 <- as.numeric(e1); e2 <- as.numeric(e2)
  if (length(e1) != length(e2)) stop("e1 y e2 deben tener la misma longitud.")
  idx <- is.finite(e1) & is.finite(e2)
  e1 <- e1[idx]; e2 <- e2[idx]
  n <- length(e1)
  if (n < 10) stop("Muy pocas observaciones tras remover NAs.")
  
  # Serie de pérdidas diferenciales d_t = L1 - L2
  d <- switch(loss,
              "se"     = e1^2 - e2^2,
              "ae"     = abs(e1) - abs(e2),
              "custom" = {
                if (is.null(d_custom) || length(d_custom) != n)
                  stop("Para loss='custom', provee d_custom de longitud n.")
                as.numeric(d_custom)
              })
  
  dbar <- mean(d)
  
  # Regresión contra constante y varianza HAC para la media
  fit <- stats::lm(d ~ 1)  # intercept-only
  if (is.null(lag)) lag <- max(0, h - 1)  # recomendación clásica para h-step-ahead
  V <- sandwich::NeweyWest(fit, lag = lag, prewhite = prewhite, adjust = TRUE)
  # Estadístico t de la constante (media de d)
  ct <- lmtest::coeftest(fit, vcov = V)
  t_stat <- as.numeric(ct[1, "t value"])
  se_hat <- as.numeric(ct[1, "Std. Error"])
  
  # --- Small-sample correction HLN (Harvey, Leybourne & Newbold 1997)
  if (small.sample) {
    # factor de ajuste
    hln_factor <- sqrt((n + 1 - 2*h + n/h) / n)
    t_stat <- t_stat * hln_factor
  }
  
  # p-valor asintótico (normal) o HLN (t-Student)
  if (small.sample) {
    p_asym <- switch(alternative,
                     "two.sided" = 2 * stats::pt(-abs(t_stat), df = n-1),
                     "greater"   = 1 - stats::pt(t_stat, df = n-1),
                     "less"      = stats::pt(t_stat, df = n-1))
  } else {
    p_asym <- switch(alternative,
                     "two.sided" = 2 * stats::pnorm(-abs(t_stat)),
                     "greater"   = 1 - stats::pnorm(t_stat),
                     "less"      = stats::pnorm(t_stat))
  }
  
  # Bootstrap de bloques (opcional) para p-valor
  p_boot <- NA_real_
  if (bootstrap) {
    if (!is.null(seed)) set.seed(seed)
    if (is.null(block_length)) block_length <- ceiling(sqrt(n))
    
    # Función para generar una serie bootstrap por bloques móviles
    moving_block_resample <- function(x, L) {
      n <- length(x)
      nb <- ceiling(n / L)           # número de bloques
      starts <- sample.int(n - L + 1, nb, replace = TRUE)
      x_star <- unlist(lapply(starts, function(s) x[s:(s + L - 1)]))
      x_star[1:n]
    }
    
    # Estadístico centrado bajo H0: E[d]=0
    t_boot <- numeric(B)
    for (b in seq_len(B)) {
      d_star <- moving_block_resample(d - dbar, block_length)  # centrado
      dbar_star <- mean(d_star)
      t_boot[b] <- dbar_star / se_hat
    }
    # p-valor bootstrap según alternativa
    p_boot <- switch(alternative,
                     "two.sided" = mean(abs(t_boot) >= abs(t_stat)),
                     "greater"   = mean(t_boot >= t_stat),
                     "less"      = mean(t_boot <= t_stat))
  }
  
  out <- list(
    method = paste0("Diebold–Mariano con HAC (Newey–West)", 
                    if (small.sample) " + small-sample HLN" else ""),
    alternative = alternative,
    h = h,
    lag = lag,
    prewhite = prewhite,
    loss = loss,
    n = n,
    estimate = c(mean_diff_loss = dbar),
    statistic = c(t = t_stat),
    p_value_asymptotic = p_asym,
    p_value_bootstrap = p_boot,
    se = se_hat,
    kernel = "Bartlett (Newey–West)",
    bootstrap = bootstrap,
    dm_hac = c(mean_diff_loss = dbar) / se_hat,
    B = if (bootstrap) B else NA_integer_,
    block_length = if (bootstrap) block_length else NA_integer_
  )
  class(out) <- "htest"
  return(out)
}


# -------------------------------
# Bai and Ng (2002) information criteria
# Y is the time series
# demean is the type of transformation
# k_max is the maximum number of factors to be tested
# -------------------------------
bai_ng <- function(Y, demean = 0, k_max = 8){
  # N and T
  n <- ncol(Y)
  tt <- nrow(Y)
  
  # minimum between N and T
  CNT2 <- min(c(n,tt))
  
  # information criteria
  IC <- matrix(0, k_max + 1, 3)
  colnames(IC) <- c("IC1", "IC2", "IC3")
  
  # estimating
  for(k in 0 : k_max){
    if(k > 0){
      pc_est <- pc(Y, k, demean)
      Phat <- pc_est$Phat
      Fhat <- pc_est$Fhat
      
      ehat <- Y - Fhat%*%t(Phat)
    }else{
      ehat <- Y
    }
    
    # variance of the errors
    Vfk <- mean(colMeans(ehat^2))
    
    # criteria
    IC[k+1, "IC1"] <- log(Vfk) + (k*((n+tt)/(n*tt)))*log((n*tt)/(n+tt))
    IC[k+1, "IC2"] <- log(Vfk) + (k*((n+tt)/(n*tt)))*log(CNT2)
    IC[k+1, "IC3"] <- log(Vfk) + (k*((log(CNT2)/CNT2)))
  }
  
  # r hat
  rhat <- apply(IC, 2, function(x) which(x == min(x)) - 1)
  
  # result
  result <- list(IC, rhat)
  names(result) <- c("IC", "rhat")
  
  return(result)
}


# -------------------------------
# number of factors procedure of Onatski (2010)
# X is the database
# demean is the type of transformation
# -------------------------------
onatski2010 <- function(X, demean = 2){  
  
  # raw data
  if(demean == 0)
    X <- X
  
  # demean data
  if(demean == 1)
    X <- scale(X, scale = FALSE)
  
  # standardized data
  if(demean == 2)
    X <- scale(X, scale = TRUE)
  
  # transpose
  X <- t(X)
  
  # T and N
  T <- ncol(X)
  N <- nrow(X)
  
  # r max (Onatski WP)
  rmax <- round(1.55*min(T^(2/5), N^(2/5)))
  
  # j initial
  j <- rmax + 1
  
  # eigen value decomposition
  XXt <- X%*%t(X)
  SVD <- eigen(XXt/T)
  
  # lambda
  lambda <- SVD$values
  
  # inital values
  rhat <- Inf
  conver <- Inf
  i <- 1
  Conver <- Inf
  
  while(length(rhat) < 4 | conver > 0 | conver < 0){
    # iterative
    
    # betahat
    betahat <- 2*abs(coef(lm(lambda[j:(j+4)] ~ c(((j - 1):(j + 3))^(2/3))))[2])
    
    # criterion (page 1008)
    cond <- which(lambda[-length(lambda)] - lambda[-1] >= betahat)
    
    if(length(cond) != 0){
      rHat <- max(cond)
      rhat[i] <- rHat  
      
      # repeat 2 and 3
      j <- rhat[i] + 1
      # convergence                               
      if(length(rhat) > 4)
        conver <- sum(diff(rhat[(length(rhat)-4):length(rhat)]))
    }else{
      rhat <- rep(0, i)
      conver <- 0
    }    
    
    i <- i + 1
  }    
  
  # rhat
  r.hat <- c(rhat[length(rhat)], betahat)
  names(r.hat) <- c("ed", "betahat")
  
  # return
  return(r.hat)
}


# -------------------------------
# eigenvalue ratio test
# X is the database
# kmax is the maximum number to test
# demean is the type of transformation
# -------------------------------
ratio.test <- function(X, kmax = 8, demean = 2){
  
  # raw data
  if(demean == 0)
    X <- X
  
  # demean data
  if(demean == 1)
    X <- scale(X, scale = FALSE)
  
  # standardized data
  if(demean == 2)
    X <- scale(X, scale = TRUE)
  
  # T and N
  T <- nrow(X)
  N <- ncol(X)                      
  
  # extracting the eigenvalues
  if(T > N){
    mu <- eigen(t(X)%*%X/(N*T))$values
  }else{
    mu <- eigen(X%*%t(X)/(N*T))$values
  }
  
  # Vfk
  Vfk <- c()
  for(i in 1 : c(kmax + 1))        
    Vfk[i] <- sum(mu[i:min(N,T)])
  
  # mock eigenvalue
  mu0 <- sum(mu[1:(min(N,T))])/(log(min(N,T))*N)
  mu <- c(mu0, mu)
  
  # test ER(k)
  ER <- mu[1:c(kmax+2)][-c(kmax+2)]/mu[1:c(kmax+2)][-1]
  if(is.null(Vfk)){
    # test GR(k)
    GR <- NULL
  }else{
    # test GR(k)
    GR <- log(1 + mu[1:c(kmax+2)][-c(kmax+2)]/c(Vfk))/
      log(1 + mu[1:c(kmax+2)][-1]/c(Vfk))
  }
  # estimators
  ker <- which(ER == max(ER)) - 1
  kgr <- which(GR == max(GR)) - 1
  
  # result
  result <- c(ker, kgr, mu)
  names(result) <- c("ker", "kgr", paste("mu_",0:min(N,T), sep = ""))
  # return
  return(result)
}


# -------------------------------
# principal components factor extraction
# Y is the matrix of time series
# r is the number of factors
# demean is the type of transformation
# -------------------------------
pc <- function(Y, r = 1, demean = 2){
  # number of time series and observations
  N <- ncol(Y)
  Ty <- nrow(Y)

  # transformation
  if(demean == 0)
    Y <- Y

  if(demean == 1)
    Y <- scale(Y, center = TRUE)

  if(demean == 2)
    Y <- scale(Y)

  # eigenvalue decomposition
  ed <- eigen(t(Y)%*%Y)

  # loading matrix
  Phat <- sqrt(N)*ed$vectors[,1:r,drop=FALSE]
  rownames(Phat) <- colnames(Y)

  # factor estimates
  Fhat <- Y%*%Phat/N

  # idiosyncratic errors
  ehat <- Y - Fhat%*%t(Phat)

  colnames(Phat) <- colnames(Fhat) <- paste("Factor",1:r, sep = "")

  # result
  result <- list(ed, Y, Phat, Fhat, ehat)
  names(result) <- c("ed", "Y", "Phat", "Fhat", "ehat")

  return(result)
}

# -------------------------------
# SPA test (Hansen, 2005) - simple implementation with block bootstrap
# Benchmark vs many alternatives
# -------------------------------

spa_test <- function(loss_bench, loss_alt,
                     B = 2000,
                     block_length = NULL,
                     studentize = TRUE,
                     seed = 123,
                     verbose = TRUE) {
  # loss_bench: numeric vector length T (benchmark losses)
  # loss_alt:   matrix/data.frame T x K of losses for alternative models
  # Returns: list with observed statistic, p-value, and per-model means
  
  set.seed(seed)
  
  loss_bench <- as.numeric(loss_bench)
  loss_alt <- as.matrix(loss_alt)
  
  if (nrow(loss_alt) != length(loss_bench)) {
    stop("loss_alt must have same number of rows as length(loss_bench).")
  }
  
  Tn <- length(loss_bench)
  K  <- ncol(loss_alt)
  
  # Differential losses: positive => alternative beats benchmark
  d <- matrix(rep(loss_bench, K), ncol = K) - loss_alt  # T x K
  
  # Drop rows with any NA
  ok <- complete.cases(d)
  d <- d[ok, , drop = FALSE]
  Tn <- nrow(d)
  
  if (Tn < 30) stop("Too few complete observations after NA removal.")
  
  # Default block length (rule-of-thumb)
  if (is.null(block_length)) {
    block_length <- max(5, floor(Tn^(1/3)))
  }
  block_length <- as.integer(block_length)
  
  # Helper: circular block bootstrap indices
  boot_indices <- function(Tn, L) {
    # number of blocks
    nb <- ceiling(Tn / L)
    starts <- sample.int(Tn, nb, replace = TRUE)
    idx <- unlist(lapply(starts, function(s) {
      ((s - 1 + 0:(L - 1)) %% Tn) + 1
    }))
    idx[1:Tn]
  }
  
  # Observed means
  dbar <- colMeans(d)
  
  # Observed test statistic: max over k of t-stat (only positive evidence matters)
  if (studentize) {
    # simple sd using sample sd; bootstrap handles dependence
    sd_k <- apply(d, 2, sd)
    sd_k[sd_k == 0] <- NA_real_
    t_obs_k <- sqrt(Tn) * dbar / sd_k
    t_obs <- max(t_obs_k, na.rm = TRUE)
  } else {
    t_obs_k <- sqrt(Tn) * dbar
    t_obs <- max(t_obs_k)
  }
  t_obs <- max(0, t_obs)  # one-sided
  
  # Centering choice (SPA-style): ensure null is benchmark not inferior
  # We use "positive-part" centering: dtilde = d - pmax(dbar, 0)
  # This is a common practical variant to avoid over-rejection.
  d_center <- sweep(d, 2, pmax(dbar, 0), FUN = "-")
  
  # Bootstrap distribution of max statistic under H0
  t_boot <- numeric(B)
  
  if (verbose) {
    message(sprintf("SPA bootstrap: T=%d, K=%d, B=%d, block_length=%d, studentize=%s",
                    Tn, K, B, block_length, studentize))
  }
  
  for (b in seq_len(B)) {
    idx <- boot_indices(Tn, block_length)
    db <- d_center[idx, , drop = FALSE]
    db_bar <- colMeans(db)
    
    if (studentize) {
      sd_b <- apply(db, 2, sd)
      sd_b[sd_b == 0] <- NA_real_
      t_b_k <- sqrt(Tn) * db_bar / sd_b
      t_boot[b] <- max(t_b_k, na.rm = TRUE)
    } else {
      t_b_k <- sqrt(Tn) * db_bar
      t_boot[b] <- max(t_b_k)
    }
    t_boot[b] <- max(0, t_boot[b])
  }
  
  # One-sided p-value
  pval <- mean(t_boot >= t_obs)
  
  out <- list(
    statistic = t_obs,
    p_value = pval,
    dbar = dbar,
    per_model_stat = t_obs_k,
    block_length = block_length,
    B = B,
    studentize = studentize
  )
  
  class(out) <- "spa_test"
  out
}

print.spa_test <- function(x, ...) {
  cat("SPA test (benchmark vs alternatives)\n")
  cat(sprintf("Statistic: %.4f\n", x$statistic))
  cat(sprintf("p-value:   %.4f\n", x$p_value))
  cat(sprintf("B=%d, block_length=%d, studentize=%s\n",
              x$B, x$block_length, x$studentize))
  invisible(x)
}


# -------------------------------
# function that incorporates the PIBO data to use it as a covariate in the DFM
# -------------------------------
gdp_o <- function(yfcst_levels, yfcst, pib, d, vari){
  
  dates_yfcst <- names(yfcst_levels)
  
  date_yfcst <- names(which(is.na(yfcst_levels))[1])
  year_yfcst <- substring(date_yfcst, 1, 4) 
  month_yfcst <- substring(date_yfcst, 6, 7) 
  
  if(month_yfcst == "03")
    quarter <- paste(year_yfcst, "/01", sep = "")
  
  if(month_yfcst == "06")
    quarter <- paste(year_yfcst, "/02", sep = "")
  
  if(month_yfcst == "09")
    quarter <- paste(year_yfcst, "/03", sep = "")
  
  if(month_yfcst == "12")
    quarter <- paste(year_yfcst, "/04", sep = "")
  
  if(month_yfcst == "03" | month_yfcst == "06" | month_yfcst == "09"| 
     month_yfcst == "12"){
    yfcst_o_levels <- yfcst_levels
    yfcst_o <- yfcst
    
    current <- yfcst_levels[(which(is.na(yfcst_levels))[1]-2):
                              (which(is.na(yfcst_levels))[1]-1)]
    
    past <- yfcst_levels[which(names(current)[1] == dates_yfcst):
                           (which(names(current)[2] == dates_yfcst)+1)-12]
    
    yfcst_o_levels[names(which(is.na(yfcst_levels)))] <- 
      (((pib[quarter, vari] + 100)/100)*mean(past) - sum(current)/3)*3
    
    timely <- cpm(yfcst_o_levels, d)
    names(timely) <- dates_yfcst[-(1:d)]
    
    yfcst_o[which(is.na(yfcst))] <- timely[names(which(is.na(yfcst_levels)))] 
  }else{
    yfcst_o <- NULL
  }
  
  return(yfcst_o)
} 

#' fit_garch11 — Fit GARCH(1,1) to a residual series; fall back to
#'               rolling std if it fails.
#' @param e  numeric vector of residuals
#' @return   numeric vector of conditional standard deviations
fit_garch11 <- function(e) {
  e <- as.numeric(e)
  e[is.nan(e) | is.infinite(e)] <- NA
  if (sum(!is.na(e)) < 20) return(rep(NA_real_, length(e)))
  
  spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = "norm"
  )
  fit <- tryCatch(
    ugarchfit(spec, data = na.omit(e) * 100, solver = "hybrid",
              solver.control = list(trace = 0)),
    error = function(err) NULL
  )
  if (is.null(fit)) {
    # fallback: rolling 6-month std
    cond_vol <- zoo::rollapply(e, width = 6, FUN = sd, na.rm = TRUE,
                               fill = NA, align = "right")
    return(as.numeric(cond_vol))
  }
  # map back to original scale
  cond_vol_full <- rep(NA_real_, length(e))
  not_na        <- which(!is.na(e))
  cond_vol_full[not_na] <- as.numeric(sigma(fit)) / 100
  cond_vol_full
}

