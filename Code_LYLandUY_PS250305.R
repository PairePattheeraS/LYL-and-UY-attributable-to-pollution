# This script calculates Life Years Lost (LYL) and Unhealthy Years (UY) attributable to pollution
# using mortality data, morbidity data, life table and population size from IHME.
# ------------------------------------------------------
# Data source: IHME (https://vizhub.healthdata.org/gbd-results/)
# Author: Pattheera (Paire) Somboonsin
# Last Updated: 05 March 2025
# The calculation details follow Somboonsin et al. (2024). DOI: 10.1016/j.atmosenv.2024.120763.
# ------------------------------------------------------

# Set working directory (Change this path based on your system)
setwd("")

# Load required libraries
library(tidyverse)
library(dplyr)

# ==============================
# FUNCTION: Compute Life Table
# ==============================
# This function estimates life expectancy and related indicators based on age- and sex-specific mortality rates (qx) provided by IHME.
LifeTableQx<-function(data,Sex){
  qx <-data
  N<-length(qx)
  AgeI<-c(1, 4, rep(5, 19))
  
  ax <- AgeI / 2
  
  if (Sex == "Male") {
    ax[1] <- ifelse(qx[1] < 0.107, 0.045 + 2.684 * qx[1], 0.330)
    ax[2] <- ifelse(qx[2] < 0.107, 1.651 - 2.816 * qx[2], 1.352)
  } else if (Sex == "Female") {
    ax[1] <- ifelse(qx[1] < 0.107, 0.053 + 2.800 * qx[1], 0.350)
    ax[2] <- ifelse(qx[2] < 0.107, 1.522 - 1.518 * qx[2], 1.361)
  } 
  
  mx<-qx/(AgeI-qx*(AgeI-ax))
  px<-1-qx
  
  lx<-100000
  dx<-lx[1]-lx[2]
  
  for (i in 1:(length(px)-1)){
    lx<-c(lx,lx[i]*px[i])
    dx<-c(dx,lx[i]-lx[i+1])}
  
  dx<-c(dx[2:21],lx[21])
  
  Lx<-(AgeI[1]*lx[2])+(ax[1]*dx[1])
  for (i in 1:(length(px)-1)){
    Lx<-c(Lx,(AgeI[i]*lx[i+1])+(ax[i]*dx[i]))}
  Lx<-Lx[2:21]
  Lx[21]<-lx[21]/mx[21]
  Tx<-rev(cumsum(rev(Lx[1:21])))
  
  ex<-Tx/lx
  
  Age<-c(0,1,seq(from=5, to=95, by=5))
  ALL<-data.frame(Age,qx,mx,lx,dx,Lx,Tx,ex,ax)
  return(ALL)
}

Calculate_a0 <- function(mx,Sex) {
  a0 <- numeric(length(mx))
  if (Sex == "Male") {
    a0[1] <- ifelse(mx[1] < 0.107, 0.045 + 2.684 * mx[1], 0.330)
    a0[2] <- ifelse(mx[2] < 0.107, 1.651 - 2.816 * mx[2], 1.352)
  } else if (Sex == "Female") {
    a0[1] <- ifelse(mx[1] < 0.107, 0.053 + 2.800 * mx[1], 0.350)
    a0[2] <- ifelse(mx[2] < 0.107, 1.522 - 1.518 * mx[2], 1.361)
  }
  
  return(a0) }

lifetable.mx<-function(mx,Sex){
  AGEL<-95
  AGEF<-0
  N<-length(mx)
  AgeI<-c(1,4,rep(5,19))
  a0<-Calculate_a0(mx,Sex)[1:2]
  ax<-c(a0,rep(2.5,(N-2)))
  if(mx[N]>0){ax[N]<-1/mx[N]}
  qx<-(AgeI*mx)/(1+(AgeI-ax)*mx)
  qx[N]<-1             
  
  px<-1-qx
  
  lx<-100000
  
  for(y in 1:(N-1)){          
    lx[y+1]<-lx[y]*px[y]
  }
  
  dx<-lx*qx
  dx[N]<-lx[N]
  
  Lx<-AgeI*lx+(ax-AgeI)*dx
  Lx[N]<-lx[N]*ax[N]                 
  
  
  Tx<-c()
  for(y in 1:N){
    Tx[y]<-sum(Lx[y:N])
  }
  
  ex<-Tx/lx
  Age<-c(0,1,seq(from=5, to=95, by=5))
  ALL<-cbind(Age,AgeI,ax,mx,qx,lx,dx,Lx,Tx,ex,px)
  return(ALL)
}

# ========================
# FUNCTION: Calculate LYL
# ========================

RxiMatrix<-function(PP,Cum){
  NumC<-dim(PP)[2]
  NumR<-dim(PP)[1]
  
  G<-PP
  FD<-colSums(PP)
  FD<-as.numeric(FD)
  
  FRx3<-t(matrix(rep(FD,(NumR)),NumC))
  FRx<-G/FRx3
  return(FRx)
}

LostYears2<-function(FLT,B){
  # Now we calculate the total number of lost years up to certain age
  AGEL<-95
  AGEF<-0
  N<-dim(FLT)[1]
  AgeI<-c(1,4,rep(5,19))
  lx<-as.numeric(FLT$lx)/100000
  dx<-c(lx[-N]-lx[-1],lx[N])
  ax<-as.numeric(FLT$ax)
  Lx<-AgeI[-N]*lx[-1]+ax[-N]*dx[-N]
  qx <-as.numeric(FLT$qx)
  if(N==21){ Lx[N]<-lx[N]*ax[N]}
  mx<-dx/Lx
  Rxi<-RxiMatrix(B,0)
  fxi<-RxiMatrix(B,1)
  
  Nrow<-dim(Rxi)[1]
  Ncol<-dim(Rxi)[2]
  ## person years and person lost
  ex<-sum(Lx[-(N)])/lx[1]
  LYL<-AgeI[-(N)]-Lx[-(N)]
  #ex+sum(LYL)
  #sum(LYL)
  ## we use the life table functions to separate the person years and 
  ## person lost
  #LYL<-AgeI-Lx
  
  LYLi<-(AgeI[1]-ax[1])*dx[1]*Rxi[,1]
  LYLi<-cbind(LYLi,AgeI[2]*(1-lx[2])*LYLi/sum(LYLi)+(AgeI[2]-ax[2])*dx[2]*Rxi[,2])
  
  for (y in 3:(length(LYL)-1)){
    LYLi<-cbind(LYLi,AgeI[y]*(1-lx[y])*rowSums(LYLi[,1:(y-1)])/(sum(rowSums(LYLi[,1:(y-1)])))+(AgeI[y]-ax[y])*dx[y]*Rxi[,y])
  }
  
  # age-specific LYLs
  Ncol<-Ncol-1
  LYLi2<-c()
  LYLi3<-matrix(0,Nrow,Ncol)
  for (y in 1:(Ncol-1)){
    L<-matrix(0,Nrow,Ncol)
    L[,y]<-(AgeI[y]-ax[y])*dx[y]*Rxi[,y]
    L[,(y+1):Ncol]<-matrix(rep(AgeI[(y+1):Ncol],each=Nrow),Nrow)*matrix(rep(dx[y]*Rxi[,y],length((y+1):Ncol)),Nrow)
    LYLi2<- cbind(LYLi2,rowSums(L))
    LYLi3<- (LYLi3+L)
  }
  y<-Ncol
  LYLi2<- cbind(LYLi2,(AgeI[y]-ax[y])*dx[y]*Rxi[,y])
  LYLi3[,y]<-LYLi3[,y]+LYLi2[,y]
  
  return(LYLi3)
  
}

# =================================================
# FUNCTION: Calculate Confidence Intervals for LYL
# =================================================

LTci <- function(dx,mx,Sex) {
  Age <-c(0,1,seq(5,95,by=5))
  AGEL<-95
  AGEF<-0
  Aage<-Age
  
  NN<-which(Aage == AGEL)
  nn<-which(Aage == AGEF)
  
  m <- length(mx)
  Ntil <- round(dx/mx)
  Y <- suppressWarnings(matrix(rbinom(m * 1000,
                                      Ntil,
                                      mx),
                               m, 1000))
  MX <- Y/Ntil
  fun.ex <- function(mx) {
    return((lifetable.mx(mx,Sex)[1,9] - lifetable.mx(mx,Sex)[NN,9])/100000)
  }
  fun.lyl <- function(mx) {
    return(AGEL-AGEF-(lifetable.mx(mx,Sex)[1,9] - lifetable.mx(mx,Sex)[NN,9])/100000)
  }
  exsim.ex <- apply(MX, 2, fun.ex)
  exsim.lyl <- apply(MX, 2, fun.lyl)
  
  ## confidence interval
  CI.ex <- quantile(exsim.ex,
                    probs = c((1-0.95)/2,0.5,
                              1 - (1-0.95)/2))
  CI.lyl <- quantile(exsim.lyl,
                     probs = c((1-0.95)/2,0.5,
                               1 - (1-0.95)/2))
  ## output
  out <- data.frame(
    CIex=CI.ex,
    CIlyl=CI.lyl
  )
  return(out)
}

LTcii2 <- function(dx,mx,Sex,B) {
  Age <-c(0,1,seq(5,95,by=5))
  AGEL<-95
  AGEF<-0
  Aage<-Age
  N<-length(mx)
  NN<-which(Aage == AGEL)
  
  mxi <- RxiMatrix(B,0)
  mxi[is.na(mxi)]<-0
  Nmxi<-dim(mxi)[1]  
  Ncol<-dim(mxi)[2]
  
  m <- N
  Ntil <- round(dx/mx)
  Y <- suppressWarnings(matrix(rbinom(m * 1000,
                                      Ntil,
                                      mx),
                               m, 1000))
  MX <- Y/Ntil
  MX[N,]<-1
  
  BM<-list() 
  BM2<-BM
  n=1000
  
  for (t in 1:Ncol){
    xx<-rmultinom(n, size =1000, prob = mxi[,t])/1000
    BM[[t]]<-lapply(seq_len(ncol(xx)), function(i) xx[,i])
  }
  BM2<- Map(cbind,BM[[1]])
  for (t in 2:Ncol){
    BM2<- Map(cbind,BM2,BM[[t]])
  }
  
  fun.lyli3 <- function(mx, BM3) {
    LT_sim <- lifetable.mx(mx, Sex)
    LT_sim <-as.data.frame(LT_sim)
    lost_years <- LostYears2(LT_sim, B = BM3)
    if (is.null(dim(lost_years))) {
      stop("LostYears2 returned a result with unexpected dimensions.")
    }
    return(rowSums(lost_years[, c(1:(NN - 1))]))
  }
  
  MX2<-lapply(seq_len(ncol(MX)), function(x) MX[,x])
  
  exsim.lyli3 <-mapply(fun.lyli3, mx=MX2, BM3=BM2)
  
  Mat <- matrix(unlist(exsim.lyli3), Nmxi)
  
  ## confidence interval
  CIlyl1<-c()
  for (i in 1:3){
    CI.lyli3 <- quantile(Mat[i,],
                         probs = c((1-0.95)/2,.5,
                                   1 - (1-0.95)/2),na.rm = TRUE)
    CIlyl1<-rbind(CIlyl1,CI.lyli3)}
  
  ## output
  out <- data.frame(CIlyl1)
  
  return(out)
}

# ================================================================================
# ANALYSE: Life Years Lost (LYL) for Different Countries, Sex and Pollution types
# ================================================================================
# Define study parameters
sex_labels<- c("Female","Male")
Sex_code <-c("F","M")
sex_code <-c("f","m")

Country_code <-c("CHN","SLB")
Country <-c("China","Solomon Islands")

Year <-2019

# Initialise empty dataframes
DTf_merge <-data.frame()
DT_LT_merge <-data.frame()

# Loop over countries and sexes
for (ct in 1:length(Country)){
  for (se in 1:length(sex_labels)){
    for (ye in 1:length(Year)){
      LT <- read.table(paste0("Data for LYL/LT", Sex_code[se], "_", Country_code[ct], Year, ".csv"), header = TRUE, fill = TRUE, sep = ",")
      LT1 <- LifeTableQx(LT$qx, sex_labels[se])  # Compute life table
      # Load cause-specific death data
      PP1 <- read.table(paste0("Data for LYL/",Country_code[ct], "_", tolower(Sex_code[se]), Year, ".csv"), header = TRUE, fill = TRUE, sep = ",")
      PP<-PP1[c(3:5),-(1:6)]
      
      Age <-c(0,1,seq(5,95,by=5))
      FLT <-LT1
      B<-PP
      
      dx<-LT1[,5]
      mx<-LT1[,3]
      ax<-LT1[,9]
      qx<-LT1[,2]
      
      print(paste("Processing year", Year[ye], "for", Country[ct], sex_labels[se]))
      print(LTci(dx,mx,sex_labels[se]))
      print(LTcii2(dx,mx,sex_labels[se],B))
      
      Result_CI <-rbind(as.matrix(t(LTci(dx, mx, sex_labels[se]))), as.matrix(LTcii2(dx, mx, sex_labels[se], B)))
      rownames(Result_CI) <- NULL
      
      DTf <-cbind(Country[ct],Country_code[ct],sex_labels[se],Year[ye],Cate=c("Life Expectancy","LYL_total","LYL_APM","LYLHAP","LYL_Others"),Result_CI)
      colnames(DTf) <-c("Country","Country_code","Sex","Year","Cate","LYL25","LYL50","LYL97")
      DTf <- as.data.frame(DTf, stringsAsFactors = FALSE)
      DTf$LYL_CI <- paste0(
        format(round(as.numeric(DTf$LYL50), 2), nsmall = 2), " (", 
        format(round(as.numeric(DTf$LYL25), 2), nsmall = 2), ",", 
        format(round(as.numeric(DTf$LYL97), 2), nsmall = 2), ")"
      )
      
      DT_LT_merge <-rbind(DT_LT_merge,FLT)
      DTf_merge <-rbind(DTf_merge,DTf)
    }
  }
}

# =======================
# FUNCTION: Calculate UY
# =======================
lifetable_qx_morbid <- function(x, qx, PM, APM, HAP, sex, ax = NULL,last_ax=2.5) {

  m <- length(x)
  n <- c(diff(x), 5)

  # Estimate ax values
  ax <- numeric(m)
  if (sex == "M") {
    ax[1] <- ifelse(qx[1] < 0.107, 0.045 + 2.684 * qx[1], 0.330)
    ax[2] <- ifelse(qx[2] < 0.107, 1.651 - 2.816 * qx[2], 1.352)
  } else if (sex == "F") {
    ax[1] <- ifelse(qx[1] < 0.107, 0.053 + 2.800 * qx[1], 0.350)
    ax[2] <- ifelse(qx[2] < 0.107, 1.522 - 1.518 * qx[2], 1.361)
  }
  
  ax[3:m] <- n[3:m] / 2  # Assign midpoints for other age groups
  
  px  <- 1-qx
  lx  <- cumprod(c(1,px))*100000
  dx  <- -diff(lx)
  # Lx  <- n*lx[-1] + ax*dx
  Lx  <- n*lx[-1] + ax*dx
  lx <- lx[-(m+1)]
  
  # Ensure no NA or Inf values
  Lx <- replace(Lx, is.na(Lx) | is.infinite(Lx), 0)
  
  # Adjust for pollution factors
  Lx_PM <- (1 - PM) * Lx
  Lx_APM <- (1 - APM) * Lx
  Lx_HAP <- (1 - HAP) * Lx
  
  # Handle missing values in pollution-adjusted Lx
  Lx_PM <- replace(Lx_PM, is.na(Lx_PM) | is.infinite(Lx_PM), 0)
  Lx_APM <- replace(Lx_APM, is.na(Lx_APM) | is.infinite(Lx_APM), 0)
  Lx_HAP <- replace(Lx_HAP, is.na(Lx_HAP) | is.infinite(Lx_HAP), 0)
  
  # Compute Tx (sum of future years lived)
  Tx <- rev(cumsum(rev(Lx)))
  Tx_PM <- rev(cumsum(rev(Lx_PM)))
  Tx_APM <- rev(cumsum(rev(Lx_APM)))
  Tx_HAP <- rev(cumsum(rev(Lx_HAP)))
  
  # Ensure Tx and lx are of equal length
  if (length(Tx) != length(lx)) {
    stop("Error: Length mismatch between Tx and lx.")
  }
  
  ex <- Tx / lx
  ex_PM <- Tx_PM / lx
  ex_APM <- Tx_APM / lx
  ex_HAP <- Tx_HAP / lx
  
  # Compute unhealthy years (UY)
  UY_PM <- ex - ex_PM
  UY_APM <- ex - ex_APM
  UY_HAP <- ex - ex_HAP
  HLE_All <- ex - (UY_PM + UY_APM + UY_HAP)
  
  # Ensure all vectors have the same length before returning the data frame
  final_length <- min(length(x), length(n), length(ax), length(qx), length(px), 
                      length(lx), length(dx), length(Lx), length(Tx), 
                      length(ex), length(ex_PM), length(ex_APM), length(ex_HAP), 
                      length(HLE_All), length(UY_PM), length(UY_APM), length(UY_HAP))
  
  return.df <- data.frame(x, n, ax, qx, px, lx, dx, Lx, Tx, ex, ex_PM, ex_APM, ex_HAP, HLE_All, UY_PM,UY_APM, UY_HAP)
  return(return.df)
}

# =================================================
# FUNCTION: Calculate Confidence Intervals for UY
# =================================================
CIex_morbid <- function(x, qx, Dx, PM, APM, HAP, sex, ax = NULL,
                        which_x=0, ns=1000, level=0.95) {
  
  m <- length(x)
  Ntil <- round(Dx / qx)
  
  # Generate simulated probabilities
  Y <- suppressWarnings(matrix(rbinom(m*ns,
                                      Ntil,
                                      qx),
                               m,ns))
  ## simulated probabilities
  QX <- Y/Ntil
  
  
  YPM <- suppressWarnings(matrix(rbinom(m*ns,
                                        Ntil,
                                        PM),
                                 m,ns))
  
  YAPM <- suppressWarnings(matrix(rbinom(m*ns,
                                         Ntil,
                                         APM),
                                  m,ns))
  
  YHAP <- suppressWarnings(matrix(rbinom(m*ns,
                                         Ntil,
                                         HAP),
                                  m,ns))
  
  QXPM <- YPM/Ntil
  QXAPM <- YAPM/Ntil
  QXHAP <- YHAP/Ntil
  
  wh <- which(x==which_x)
  
  exsimA <- rep(0,ns)
  for(s in 1:ns){
    exsimA[s] <-lifetable_qx_morbid(x, qx=QX[,s], PM=QXPM[,s], APM=QXAPM[,s], HAP=QXHAP[,s], sex= sex,
                                    last_ax=last_ax)$HLE_All[wh]
  }
  exsimB <- rep(0,ns)
  for(s in 1:ns){
    exsimB[s] <-lifetable_qx_morbid(x, qx=QX[,s], PM=QXPM[,s], APM=QXAPM[,s], HAP=QXHAP[,s], sex= sex,
                                    last_ax=last_ax)$UY_PM[wh]
  }
  exsimC <- rep(0,ns)
  for(s in 1:ns){
    exsimC[s] <-lifetable_qx_morbid(x, qx=QX[,s], PM=QXPM[,s], APM=QXAPM[,s], HAP=QXHAP[,s], sex= sex,
                                    last_ax=last_ax)$UY_APM[wh]
  }
  exsimD <- rep(0,ns)
  for(s in 1:ns){
    exsimD[s] <-lifetable_qx_morbid(x, qx=QX[,s], PM=QXPM[,s], APM=QXAPM[,s], HAP=QXHAP[,s], sex= sex,
                                    last_ax=last_ax)$UY_HAP[wh]
  }
  
  exsimE <- rep(0,ns)
  for(s in 1:ns){
    exsimE[s] <-lifetable_qx_morbid(x, qx=QX[,s], PM=QXPM[,s], APM=QXAPM[,s], HAP=QXHAP[,s], sex= sex,
                                    last_ax=last_ax)$ex[wh]
  }
  
  # Compute confidence intervals
  compute_CI <- function(sim) {
    quantile(sim, probs = c((1 - level) / 2, 1 - (1 - level) / 2))
  }
  
  CIHealthyAll <- compute_CI(exsimA)
  CIUY_PM <- compute_CI(exsimB)
  CIUY_APM <- compute_CI(exsimC)
  CIUY_HAP <- compute_CI(exsimD)
  CIex <- compute_CI(exsimE)
  
  # Output results as a data frame
  output_UY <-data.frame(
    CIHealthyAll = c(CIHealthyAll[1], CIHealthyAll[2]),
    CIUY_PM = c(CIUY_PM[1], CIUY_PM[2]),
    CIUY_APM = c(CIUY_APM[1], CIUY_APM[2]),
    CIUY_HAP = c(CIUY_HAP[1], CIUY_HAP[2]),
    CIex = c(CIex[1], CIex[2]))
  
  return(output_UY)
  
}

# ===============================================================================
# ANALYSE: Unhealthy Years (UY) for Different Countries, Sex and Pollution types
# ===============================================================================
# Define study parameters
country_codes <- c("SLB", "CHN")
countries <- c("Solomon Islands", "China")
sex_labels <- c("Female", "Male")
sex_codes <- c("F", "M")
sex_codes_lower <- c("f", "m")
pollutants <- c("PM", "aPM", "HAP")
years <- c(2019)

# Initialise an empty data frame to store results
UY_merge <- data.frame()

# Loop through country, sex, and year
for (country_idx in seq_along(country_codes)) {
  for (sex_idx in seq_along(sex_codes)) {
    for (year in years) {
      
      # Read input data
      pop_data <- read.csv(paste0("Data for UY/",country_codes[country_idx], "_POP", sex_codes[sex_idx], year, ".csv"), header = TRUE, fill = TRUE)
      yld_all <- read.csv(paste0("Data for UY/",country_codes[country_idx], "_YLD", sex_codes_lower[sex_idx], year, ".csv"), header = TRUE, fill = TRUE)
      LT <- read.csv(paste0("Data for UY/","LT",sex_codes[sex_idx],"_",country_codes[country_idx],year,".csv",sep=""),header=TRUE,fill=TRUE,sep=",")
      
      yld_pollute <- list()
      for (p in seq_along(pollutants)) {
        yld_pollute[[pollutants[p]]] <- read.csv(paste0("Data for UY/",country_codes[country_idx], "_YLD", sex_codes_lower[sex_idx], year, "_", pollutants[p], ".csv"), header = TRUE, fill = TRUE)
      }
      
      death_data <- read.csv(paste0("Data for UY/",country_codes[country_idx], "_DX", sex_codes_lower[sex_idx], year, ".csv"), header = TRUE, fill = TRUE)
      
      # Create morbidity dataset
      morbidity_data <- data.frame(
        Pop = pop_data$val,
        YLD_All = yld_all$val,
        PM = yld_pollute[["PM"]]$val,
        aPM = yld_pollute[["aPM"]]$val,
        HAP = yld_pollute[["HAP"]]$val
      )
      
      # Compute proportion of healthy years
      health_ratios <- morbidity_data %>%
        mutate(
          All = 1 - (YLD_All / Pop),
          PM = 1 - (PM / Pop),
          aPM = 1 - (aPM / Pop),
          HAP = 1 - (HAP / Pop))
      
      # Compute life tables
      Lh <- LT$Lx * health_ratios[, c("All", "PM", "aPM", "HAP")]
      Tx <- as.data.frame(lapply(Lh, function(x) rev(cumsum(rev(x)))))
      eh <- Tx / LT$lx
      
      # Compute unhealthy years (UY)
      UY <- LT$ex - eh
      colnames(UY) <- c("UY_All", "UY_PM", "UY_APM", "UY_HAP")
      
      # Compute disability burden
      MorbData <- morbidity_data %>%
        mutate(
          All = YLD_All / Pop,
          PM = PM / Pop,
          APM = aPM / Pop,
          HAP = HAP / Pop
        ) %>%
        select(All, PM, APM, HAP)
      
      # Bind life table data
      Databind <- cbind(LT, MorbData, Dx = death_data[10])
      Databind$Country_code <- country_codes[country_idx]
      Databind$Sex <- sex_codes[sex_idx]
      
      # Print progress
      print(paste0("Processing year ", year, " for ", countries[country_idx], " (", sex_labels[sex_idx], ")"))
      
      # Compute life table with morbidity
      morbidDT <- lifetable_qx_morbid(
        x = Databind$Age, 
        qx = Databind$qx, 
        PM = Databind$PM, 
        APM = Databind$APM, 
        HAP = Databind$HAP, 
        sex = Databind$Sex[1]
      )
      
      #print(morbidDT)
      
      # Compute confidence intervals
      CIdata <- CIex_morbid(
        x = Databind$Age, qx = Databind$qx, Dx = Databind$val, 
        PM = Databind$PM, APM = Databind$APM, HAP = Databind$HAP, 
        sex = Databind$Sex[1], which_x = 0
      )
      
      CIdata[3, ] <- c(morbidDT$HLE_All[1], morbidDT$UY_PM[1], morbidDT$UY_APM[1], morbidDT$UY_HAP[1], morbidDT$ex[1])
      
      # Reshape data for final output
      UY_data <- as.data.frame(t(CIdata))
      colnames(UY_data) <- c("UY_low", "UY_high", "UY_mid")
      UY_data$Cate <- c("Healthy Life Expectancy", "UY_PM", "UY_APM", "UY_HAP", "Life Expectancy")
      rownames(UY_data) <-NULL
      
      # Format CI column
      UY_data$UY_CI <- paste0(
        format(round(UY_data$UY_mid, 2), nsmall = 2), " (",
        format(round(UY_data$UY_low, 2), nsmall = 2), ", ",
        format(round(UY_data$UY_high, 2), nsmall = 2), ")"
      )
      
      # Add metadata columns
      UY_data <- UY_data %>%
        mutate(
          Country_code = country_codes[country_idx],
          Country = countries[country_idx],
          Year = year,
          Sex = sex_labels[sex_idx]
        ) %>%
        select(Country, Country_code, Sex, Year, Cate, UY_low, UY_mid, UY_high, UY_CI)
      
      # Print and merge results
      print(UY_data)
      UY_merge <- bind_rows(UY_merge, UY_data)
    }
  }
}

