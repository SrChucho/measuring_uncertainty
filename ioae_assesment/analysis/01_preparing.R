library(seasonal) |> suppressMessages()
library(imputeTS) |> suppressMessages()

# -------------------------------
# seed
set.seed(12345)
# -------------------------------

# -------------------------------
# parameters
vari <- "IGAE"
start <- "2004/01"
d <- 1
start_nowcasting <- "2019/01"
# -------------------------------

# data: catalog
cat_dfm <- read.csv(here::here("Data/DFM/Catalogo.csv"))
cat_yfcst <- read.csv(here::here("Data/YFCST/Catalogo.csv"))
cat_gt <- read.csv(here::here("Data/Google Trends/Catalogo.csv"))

# data: read time series
db_dfm_raw <- read.csv(here::here("Data/DFM/Datos.csv"), row.names = 1)
db_yfcst_raw <- read.csv(here::here("Data/YFCST/Datos.csv"), row.names = 1)
db_yfcst_raw$IGAE <- trunc(db_yfcst_raw$IGAE * 1e12) / 1e12
pib <- read.csv(here::here("Data/DFM/PIB.csv"), row.names = 1)
google <- read.csv(here::here("Data/Google Trends/Google_Trends.csv"), row.names = 1)

# active time series
act <- as.character(cat_dfm[which(cat_dfm[,vari] == 1), "Short"])
db_dfm_raw <- db_dfm_raw[,act]

Tf <- nrow(db_dfm_raw)
#db_dfm_raw[(Tf-1):Tf,]

# dates
dates_dfm <- rownames(db_dfm_raw)
dates_yfcst <- rownames(db_yfcst_raw)

Ht <- length(which(dates_yfcst == start_nowcasting):length(dates_yfcst))

# raw data set beginnig with the parameter start
db_dfm_raw_s_fit <- db_dfm_raw[which(start == dates_dfm):length(dates_dfm),]
db_dfm_raw_s_fit <- na_locf(db_dfm_raw_s_fit)

# T raw
Tr <- sum(!is.na(db_yfcst_raw[,vari]))

# models specification
if(is.null(d)){
  d_alt_1 <- 1
  d_alt_2 <- 12
  
  alt_1_real <- cpm(db_yfcst_raw[1:Tr,vari], d_alt_1)
  alt_2_real <- cpm(db_yfcst_raw[1:Tr,vari], d_alt_2)
  
  name <- "LEV"
}else{
  
  if(d == 1){
    d_alt_1 <- 0
    d_alt_2 <- 12
    
    alt_1_real <- db_yfcst_raw[1:Tr,vari]
    alt_2_real <- cpm(db_yfcst_raw[1:Tr,vari], d_alt_2)
    
    name <- "MV"
  }
  
  if(d == 12){
    d_alt_1 <- 0
    d_alt_2 <- 1
    
    alt_1_real <- db_yfcst_raw[1:Tr,vari]
    alt_2_real <- cpm(db_yfcst_raw[1:Tr,vari], d_alt_2)
    
    name <- "AV"
  }
}

# start year
year_s <- as.numeric(substring(start, 1, 4)[1])

# catalog with active variables
cat_dfm <- cat_dfm[cat_dfm[,vari] == 1,]

# seasonal adjustment of series that require it
sa_dfm <- as.character(cat_dfm[which(cat_dfm[,"SA"] == 1), "Short"])

for(i in 1 : length(sa_dfm)){ 
  
  # date
  date_i <- rownames(db_dfm_raw)[!is.na(db_dfm_raw[,sa_dfm[i]])]
  
  # time series
  ts_dfm_i <- ts(db_dfm_raw[,sa_dfm[i]],  
                 start = c(as.numeric(substring(date_i, 1, 4))[1], as.numeric(substring(date_i, 6, 7)))[1], 
                 frequency = 12)
  
  # seasonal adjustment
  db_dfm_raw[!is.na(ts_dfm_i),sa_dfm[i]] <- seas(ts_dfm_i)$series$s11
  
  # restriction of negatives
  db_dfm_raw[which(db_dfm_raw[,sa_dfm[i]] < 0), sa_dfm[i]] <- 0.01 
}

# transformation according to the data
type_trans <- as.data.frame(matrix(0, ncol(db_dfm_raw), 2))
colnames(type_trans) <- c("Transformation", "Correlation")
rownames(type_trans) <- colnames(db_dfm_raw)

if(!is.null(d)){
  if(d != 0){
    db_yfcst_trans <- apply(db_yfcst_raw, 2, function(x) cpm(x, d))
    rownames(db_yfcst_trans) <- dates_yfcst[-(1:d)]
  }
  if(d == 0){
    db_yfcst_trans <- db_yfcst_raw
    rownames(db_yfcst_trans) <- dates_yfcst
  }
  db_dfm_trans <- matrix(NA, nrow(db_dfm_raw), ncol(db_dfm_raw))
  rownames(db_dfm_trans) <- rownames(db_dfm_raw)
  colnames(db_dfm_trans) <- colnames(db_dfm_raw)
  
  for(i in 1 : ncol(db_dfm_trans)){ 
    transform_data <- best_trans(db_yfcst_trans[,vari], db_dfm_raw[,i],
                                 rownames(db_yfcst_trans), rownames(db_dfm_raw), 
                                 as.character(cat_yfcst[which(cat_yfcst[,"Short"] 
                                                              == vari), colnames(db_dfm_raw)[i]]))
    
    
    db_dfm_trans[,i] <-  transform_data$x_trans
    type_trans[i, "Transformation"] <- transform_data$trans
    type_trans[i, "Correlation"] <- transform_data$rho[transform_data$trans]
  }
}else{
  db_yfcst_trans <- db_yfcst_raw
  
  db_dfm_trans <- matrix(NA, nrow(db_dfm_raw), ncol(db_dfm_raw))
  rownames(db_dfm_trans) <- rownames(db_dfm_raw)
  colnames(db_dfm_trans) <- colnames(db_dfm_raw)
  
  for(i in 1 : ncol(db_dfm_trans)){ 
    transform_data <- best_trans(db_yfcst_trans[,vari], db_dfm_raw[,i],
                                 rownames(db_yfcst_trans), rownames(db_dfm_raw), 
                                 as.character(cat_yfcst[which(cat_yfcst[,"Short"] 
                                                              == vari), colnames(db_dfm_raw)[i]]))
    
    
    db_dfm_trans[,i] <-  transform_data$x_trans
    type_trans[i, "Transformation"] <- transform_data$trans
    type_trans[i, "Correlation"] <- transform_data$rho[transform_data$trans]
  }
}

# data initial date
start_ycst <- which(rownames(db_yfcst_trans) == start)
start_dfm <- which(rownames(db_dfm_trans) == start)

db_yfcst <- ts(db_yfcst_trans[start_ycst:nrow(db_yfcst_trans),], start = year_s,
               frequency = 12)
db_dfm <- ts(db_dfm_trans[start_dfm:nrow(db_dfm_trans),], start = year_s,
             frequency = 12)

# model parameters
Ty <- length(db_yfcst[!is.na(db_yfcst)])

# data set
db_dfm_p <- db_dfm
db_dfm_p[1:Ty, ] <- na_locf(db_dfm_p[1:Ty, ]) 

# incorporation of the variable to be predicted 
yfcst <- db_yfcst

yfcst_levels <- db_yfcst_raw[,vari]
names(yfcst_levels) <- dates_yfcst

yfcst_levels <- yfcst_levels[which(dates_yfcst == start):length(yfcst_levels)]
names(yfcst) <- names(yfcst_levels)

yfcst_annual <- cpm(yfcst_levels, 12)
names(yfcst_annual) <- names(yfcst)[-(1:12)]

vari_o <- vari

## writing prepared inputs
dir.create(here::here("outputs", "rds", "prepared_inputs"),
           recursive = TRUE, showWarnings = FALSE)

prepared_inputs <- list(
  # Scalars / params
  vari = vari, vari_o = vari_o, d = d, year_s = year_s, Ht = Ht, Ty = Ty,
  
  # Catalog
  cat_dfm = cat_dfm,
  
  # Core data
  db_dfm_raw = db_dfm_raw,
  db_dfm_p   = db_dfm_p,
  db_yfcst   = db_yfcst,
  yfcst      = yfcst,
  yfcst_levels = yfcst_levels,
  yfcst_annual = yfcst_annual,
  pib = pib,
  google = google,
  
  # Aux
  type_trans = type_trans
)

saveRDS(prepared_inputs, here::here("outputs", "rds", "prepared_inputs", "prepared_inputs.rds"))
