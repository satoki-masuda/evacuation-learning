### Evacuation simulation
### @author:Satoki Masuda
### encoding: UTF-8
### 2022/1/10

rm(list=ls(all=TRUE)) #delete all variables
library(tictoc)
library(Matrix)
library(dplyr)
library(compiler)
library(parallel)
library(ggplot2)
library(ggrepel)

createEmptyDf = function(nrow, ncol, colnames = c()){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}

### Settings ------------------
set.seed(123)
#モンテカルロシミュレーションの回数
dim<-1000

# Parameter
par <- c(-2.227076383, 0.712642215, -1.292325166, 6.33181656, -0.556804972, -0.51096569, -0.169770418, 1.285324338, 2.496594752, 0.917186604, -1.219138708, -1.543252011, 1.0)

project <- "YOUR FOLDER PATH"
#data for individual LOS
inputTSNWpath <- paste(project, "/input/TSNW/", sep = "")
personal <- read.csv(paste(project, "/input/PersonalInfo.csv", sep=""), fileEncoding = "Shift-jis")

estimates1 <- read.csv(paste(project, "/result/1_safe_evac.csv", sep=""))$estimates
estimates2 <- read.csv(paste(project, "/result/2_safe_nevac.csv", sep=""))$estimates
estimates3 <- read.csv(paste(project, "/result/3_danger_evac.csv", sep=""))$estimates
estimates4 <- read.csv(paste(project, "/result/4_danger_nevac.csv", sep=""))$estimates

learning_hessian1 <- read.csv(paste(project, "/result/hessian1_safe_evac.csv", sep=""))
learning_hessian2 <- read.csv(paste(project, "/result/hessian2_safe_nevac.csv", sep=""))
learning_hessian3 <- read.csv(paste(project, "/result/hessian3_danger_evac.csv", sep=""))
learning_hessian4 <- read.csv(paste(project, "/result/hessian4_danger_nevac.csv", sep=""))

setwd(project)

theta <- par[length(par)] #discount factor
li_theta <- c(1, theta^4, theta^2, theta, 1)


# Probability func --------------
prob <- function(par, learning_par, tsnw, g1_drill, g1_cong, g2_cong, g3_drill, g4_map){
  
  #time-space NW node data
  node_space <- sort(unique(tsnw$a_space))
  N <- length(node_space)
  node_TSN <- c(sort(unique(c(tsnw$k, tsnw$a))), 99999) #absorbing state
  Nt <- length(node_TSN)
  
  #link to absorbing state
  absLink <- createEmptyDf(nrow = N, ncol = ncol(tsnw), colnames = colnames(tsnw)) # absorbing link
  absLink[,] <- 0
  absLink$k <- node_space + 40000
  absLink[,c("a", "a_space")] <- 99999
  absLink$k_space <- node_space
  tsnw <- rbind(tsnw, absLink)
  tsnw$linkID <- 1:nrow(tsnw)
  L <- nrow(tsnw)
  
  #adjacent matrix I of time-structured NW
  I <- array(0, dim = c(Nt, Nt))
  for(i in 1:nrow(tsnw)){
    I[(1:Nt)[node_TSN == tsnw$k[i]],(1:Nt)[node_TSN == tsnw$a[i]]] <- 1
  }
  

  #------------------- Variables -------------------

  #calculate utility with exponential form (note!)
  tsnw[tsnw$distance>0, "distance"] <- log(10*tsnw[tsnw$distance>0, "distance"])/10 #距離は0.1-900kmくらいの間でばらついてる
  tsnw[, c("depth", "home_depth")] <- tsnw[, c("depth", "home_depth")]/10
  tsnw[, "subj_prob"] <- tsnw[, "subj_prob"]/100
  
  #for value function of the last node
  #depth * non-home dummy
  dst_depth <- array(exp(0.0), dim = c(Nt,Nt))
  #school * non-home dummy
  dst_school <- dst_depth
  #living in inundation area dummy * under 2nd floor dummy * home dummy
  dst_home_low <- dst_depth
  #ASC
  dst_home <- dst_depth
  #interest
  dst_interest <- dst_depth
  #evacuation order * non-home dummy
  order <- dst_depth
  #evacuation order * low * non-home dummy
  order_low <- dst_depth
  subj_prob1 <- dst_depth
  
  #distance * non-home dummy
  distance_1 <- dst_depth
  distance_2 <- dst_depth
  distance_3 <- dst_depth
  distance_4 <- dst_depth
  
  #後で使う
  congestion <- dst_depth
  flood_time <- dst_depth
  drill <- dst_depth
  
  for (i in 1:L){
    k_i <- (1:Nt)[node_TSN == tsnw$k[i]]
    a_i <- (1:Nt)[node_TSN == tsnw$a[i]]
    
    if (tsnw$a[i]==99999){
      if (tsnw[tsnw$a==tsnw$k[i], "home"][1]==0){
        dst_depth[k_i,a_i] <- exp(tsnw[tsnw$a==tsnw$k[i], "depth"][1])
        dst_school[k_i,a_i] <- exp((tsnw[tsnw$a==tsnw$k[i], "a_type"][1]=="primary_school") | (tsnw[tsnw$a==tsnw$k[i], "a_type"][1]=="school"))
        if (tsnw$distance[i]<1){
          congestion[k_i, a_i] <- exp(1)
        }
      }
      else{
        dst_home_low[k_i, a_i] <- exp((tsnw[1, "home_depth"]>0) * (tsnw[1, "height"]<=2))
        dst_home[k_i, a_i] <- exp(1)
        dst_interest[k_i, a_i] <- exp(6 - tsnw[1, "bousai_interest"])
        drill[k_i, a_i] <- exp(1)
      }
      flood_time[k_i,a_i] <- exp(tsnw[(tsnw$a==tsnw$k[i]), "depth"][1]>0)
    }
    
    order[k_i,a_i] <- exp((tsnw[tsnw$a_space==tsnw$k_space[i], "home"][1]) * (tsnw$home[i]==0) *
                            ((tsnw$koiki_order[i]) + (tsnw$level4[i])))
    order_low[k_i,a_i] <- exp((tsnw[tsnw$a_space==tsnw$k_space[i], "home"][1]) * (tsnw$home[i]==0) *
                                (tsnw$home_depth[i]>0) * (tsnw$height[i]<=2) *
                                ((tsnw$koiki_order[i]) + (tsnw$level4[i])))
    
    distance_1[k_i,a_i] <- exp(tsnw$distance[i] * (tsnw$home[i]==0) * (tsnw$a[i]<2000))
    distance_2[k_i,a_i] <- exp(tsnw$distance[i] * (tsnw$home[i]==0) * ((tsnw$a[i]>2000)&(tsnw$a[i]<3000)))
    distance_3[k_i,a_i] <- exp(tsnw$distance[i] * (tsnw$home[i]==0) * ((tsnw$a[i]>3000)&(tsnw$a[i]<4000)))
    distance_4[k_i,a_i] <- exp(tsnw$distance[i] * (tsnw$home[i]==0) * ((tsnw$a[i]>4000)&(tsnw$a[i]<5000)))
    subj_prob1[k_i,a_i] <- exp((tsnw$subj_prob[i]) * (tsnw$home[i]==1))
  }
  

  # --------- log-likelihood function ---------------------
  # utility function matrix M
  M <- array(0, dim = c(Nt, Nt))
  M <- I*
    (dst_depth^par[1])*(dst_school^par[2])*
    (dst_home_low^par[3])*(dst_home^par[4])*(dst_interest^par[5])*(subj_prob1^par[6])*
    (order^par[7])*(order_low^par[8])*
    (distance_1^par[9])*(distance_2^par[10])*(distance_3^par[11])*(distance_4^par[12])
  
  
  # value function V=log(z)
  # calculate z(t)=Mz(t+1)+b from t=T until t=0 backward
  
  z <- array(1, dim = c(1, Nt)) #exp(Vd). value function of the absorbing step
  
  N_ts <- c(sum(node_TSN<10000), #num of nodes at time step 0 (home)
            sum((node_TSN>10000)&(node_TSN<20000)), #num of nodes at time step 1
            sum((node_TSN>20000)&(node_TSN<30000)), #num of nodes at time step 2
            sum((node_TSN>30000)&(node_TSN<40000)), #num of nodes at time step 3
            sum((node_TSN>40000)&(node_TSN<50000)), #num of nodes at time step 4
            sum(node_TSN>90000)) #num of absorbing nodes
  
  idx_1 <- Nt - N_ts[6] + 1
  idx_2 <- Nt

  for (i in 5:1){ #num of step (4) + time step 0 (1)
    z_i <- matrix(M[(idx_1 - N_ts[i]):(idx_2 - N_ts[i+1]), (idx_1):(idx_2)], c(N_ts[i], N_ts[i+1])) %*% (z[(idx_1):(idx_2)]^(li_theta[i]))
    idx_1 <- idx_1 - N_ts[i]
    idx_2 <- idx_2 - N_ts[i+1]
    z[idx_1:idx_2] <- (z_i == 0)*1 + (z_i != 0)*z_i
  } 
  
  V <- z#log(z) #delete exp()
  
  #calculate route choice probability
  p <- 0.0*dst_depth
  for(i in 1:nrow(tsnw)){
    k_i <- (1:Nt)[node_TSN == tsnw$k[i]]
    a_i <- (1:Nt)[node_TSN == tsnw$a[i]]
    theta_i <- (tsnw$a[i]%/%10000==1)*li_theta[1]+(tsnw$a[i]%/%10000==2)*li_theta[2]+(tsnw$a[i]%/%10000==3)*li_theta[3]+(tsnw$a[i]%/%10000==4)*li_theta[4]+(tsnw$a[i]==99999)*li_theta[5]
    p[k_i,a_i] <- (M[k_i,a_i] * V[a_i]^theta_i) / V[k_i]
  }
  
  #to absorbing node
  tsnw[(tsnw$a==99999)&(tsnw$k == tsnw[(tsnw$a %/% 10000==4)&(tsnw$prev_choice==1),"a"]),"prev_choice"] <- 1
  
  # calculate the probability of non-evacuation
  p_nev <- numeric(4)
  link_home <- tsnw[tsnw$home==1,c("k","a")]
  for (i in 1:nrow(link_home)){
    if (i==1){
      p_nev[i] <- p[(1:Nt)[node_TSN == link_home$k[i]],(1:Nt)[node_TSN == link_home$a[i]]]
    }else{
      p_nev[i] <- p_nev[i-1] * p[(1:Nt)[node_TSN == link_home$k[i]],(1:Nt)[node_TSN == link_home$a[i]]]
    }
  }
  
  
  #------------------ 情報を出した場合の効用関数の更新---------------------------------
  tsnw$prev_choice_info <- 0
  prev <- array(exp(0.0), dim = c(Nt,Nt))
  
  df_prev <- tsnw[tsnw$prev_choice==1,]
  for (i in 1:nrow(df_prev)){
    k_i <- (1:Nt)[node_TSN == df_prev$k[i]]
    a_i <- (1:Nt)[node_TSN == df_prev$a[i]]
    prev[k_i,a_i] <- exp(1)
  }
  
  if ((g1_drill==TRUE) & (tsnw$home_depth[1]==0) & (p_nev[4]<0.5)){#安全避難
    M <- M * (drill ^ learning_par[2])*(prev^learning_par[7])#避難訓練
  }else if ((g1_cong==TRUE) & (tsnw$home_depth[1]==0) & (p_nev[4]<0.5)){#安全避難
    M <- M * (congestion ^ learning_par[3])*(prev^learning_par[7])#混雑
  }else if ((g2_cong==TRUE) & (tsnw$home_depth[1]==0) & (p_nev[4]>0.5)){#安全非避難
    M <- M * (congestion ^ learning_par[3])*(prev^learning_par[7]) #混雑
  }else if ((g3_drill==TRUE) & (tsnw$home_depth[1]>0) & (p_nev[4]<0.5)){#危険避難
    M <- M * (drill ^ learning_par[2])*(prev^learning_par[7]) #避難訓練
  }else if ((g4_map==TRUE) & (tsnw$home_depth[1]>0) & (p_nev[4]>0.5)){#危険非避難
    M <- M * (flood_time ^ learning_par[1])*(prev^learning_par[7])#ハザードマップ
  }else{
    M <- M #情報を出しても変わらない
  }
  
  z <- array(1, dim = c(1, Nt)) #exp(Vd). value function of the absorbing step
  
  idx_1 <- Nt - N_ts[6] + 1
  idx_2 <- Nt
  
  for (i in 5:1){ #num of step (4) + time step 0 (1)
    z_i <- matrix(M[(idx_1 - N_ts[i]):(idx_2 - N_ts[i+1]), (idx_1):(idx_2)], c(N_ts[i], N_ts[i+1])) %*% (z[(idx_1):(idx_2)] ^(li_theta[i]))
    idx_1 <- idx_1 - N_ts[i]
    idx_2 <- idx_2 - N_ts[i+1]
    z[idx_1:idx_2] <- (z_i == 0)*1 + (z_i != 0)*z_i
  }
  
  V <- z#log(z) #delete exp()
  
  #calculate route choice probability
  p <- 0.0*dst_depth
  for(i in 1:nrow(tsnw)){
    k_i <- (1:Nt)[node_TSN == tsnw$k[i]]
    a_i <- (1:Nt)[node_TSN == tsnw$a[i]]
    theta_i <- (tsnw$a[i]%/%10000==1)*li_theta[1]+(tsnw$a[i]%/%10000==2)*li_theta[2]+(tsnw$a[i]%/%10000==3)*li_theta[3]+(tsnw$a[i]%/%10000==4)*li_theta[4]+(tsnw$a[i]==99999)*li_theta[5]
    p[k_i,a_i] <- (M[k_i,a_i] * V[a_i]^theta_i) / V[k_i]
  }
  
  # calculate the probability of non-evacuation
  p_nev <- numeric(4)
  link_home <- tsnw[tsnw$home==1,c("k","a")]
  for (i in 1:nrow(link_home)){
    if (i==1){
      p_nev[i] <- p[(1:Nt)[node_TSN == link_home$k[i]],(1:Nt)[node_TSN == link_home$a[i]]]
    }else{
      p_nev[i] <- p_nev[i-1] * p[(1:Nt)[node_TSN == link_home$k[i]],(1:Nt)[node_TSN == link_home$a[i]]]
    }
  }
  
  
  result <- list(p_nev)
  return(result)
}


# Calculation of evacuation rate -------
df1 <- data.frame(timing=c(-48, -24, -12, -6), rate = rep(0,4), low = rep(0,4), high = rep(0,4), group = as.factor("g1_without"))
TSNW_file <- list.files(inputTSNWpath, pattern="wave1")[which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1))]
n_data <- length(TSNW_file)
for (i_data in 1:n_data){
  TSNW <- read.csv(paste(inputTSNWpath, TSNW_file[i_data], sep = ""))
  result <- prob(par=par, learning_par = rep(0,7), tsnw=TSNW, g1_drill=F, g1_cong=F, g2_cong=F, g3_drill=F, g4_map=F)
  df1[,"rate"] <- df1[,"rate"] + c(1-result[[1]])
}
df1[,"rate"] <- df1[,"rate"]/n_data

df3 <- data.frame(timing=c(-48, -24, -12, -6), low = rep(0,4), high = rep(0,4), rate = rep(0,4), group = as.factor("g3_without"))
TSNW_file <- list.files(inputTSNWpath, pattern="wave1")[which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1))]
n_data <- length(TSNW_file)
for (i_data in 1:n_data){
  TSNW <- read.csv(paste(inputTSNWpath, TSNW_file[i_data], sep = ""))
  result <- prob(par=par, learning_par = rep(0,7), tsnw=TSNW, g1_drill=F, g1_cong=F, g2_cong=F, g3_drill=F, g4_map=F)
  df3[,"rate"] <- df3[,"rate"] + c(1-result[[1]])
}
df3[,"rate"] <- df3[,"rate"]/n_data

df4 <- data.frame(timing=c(-48, -24, -12, -6), low = rep(0,4), high = rep(0,4), rate = rep(0,4), group = as.factor("g4_without"))
TSNW_file <- list.files(inputTSNWpath, pattern="wave1")[which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0))]
n_data <- length(TSNW_file)
for (i_data in 1:n_data){
  TSNW <- read.csv(paste(inputTSNWpath, TSNW_file[i_data], sep = ""))
  result <- prob(par=par, learning_par = rep(0,7), tsnw=TSNW, g1_drill=F, g1_cong=F, g2_cong=F, g3_drill=F, g4_map=F)
  df4[,"rate"] <- df4[,"rate"] + c(1-result[[1]])
}
df4[,"rate"] <- df4[,"rate"]/n_data



df5 <- data.frame(timing=c(-48, -24, -12, -6), low = rep(0,4), high = rep(0,4), rate = rep(0,4), group = as.factor("g1_drill"))
TSNW_file <- list.files(inputTSNWpath, pattern="wave1")[which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1))]
n_data <- length(TSNW_file)

hes1.col <- t(chol(learning_hessian1))
#パラメータ推定値が多変量正規分布に従うと仮定して、その分布からパラメータをサンプリング
par_combination <- matrix(rnorm(dim*length(estimates1), 0, 1), nrow = dim, ncol=length(estimates1)) %*% hes1.col + matrix(rep(estimates1, dim), nrow=dim, ncol=length(estimates1), byrow=T)
#学習直後の効果を見るので，記憶率は一定
par_combination[,4] <- estimates1[4]
par_combination[,5] <- estimates1[5]
par_combination[,6] <- estimates1[6]
rate_list <- matrix(rep(0, dim*4), nrow=dim, ncol=4)

for (par_i in 1:dim){
  if (par_i %% 100 == 0){
    print(par_i)
  }
  learning_par <- par_combination[par_i,]
  
  for (i_data in 1:n_data){
    TSNW <- read.csv(paste(inputTSNWpath, TSNW_file[i_data], sep = ""))
    result <- prob(par=par, learning_par = learning_par, tsnw=TSNW, g1_drill=T, g1_cong=F, g2_cong=F, g3_drill=F, g4_map=F)
    rate_list[par_i,] <- rate_list[par_i,] + c(1-result[[1]])
  }
  rate_list[par_i,] <- rate_list[par_i,]/n_data
}
df5[,"rate"] <- apply(rate_list, 2, mean)#列方向に平均をとる
df5[,"low"] <- apply(rate_list, 2, function(x) quantile(x, 0.025))#列方向に2.5%点をとる
df5[,"high"] <- apply(rate_list, 2, function(x) quantile(x, 0.975))#列方向に97.5%点をとる



df6 <- data.frame(timing=c(-48, -24, -12, -6), low = rep(0,4), high = rep(0,4), rate = rep(0,4), group = as.factor("g1_cong"))
TSNW_file <- list.files(inputTSNWpath, pattern="wave1")[which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1))]
n_data <- length(TSNW_file)

hes1.col <- t(chol(learning_hessian1))
#パラメータ推定値が多変量正規分布に従うと仮定して、その分布からパラメータをサンプリング
par_combination <- matrix(rnorm(dim*length(estimates1), 0, 1), nrow = dim, ncol=length(estimates1)) %*% hes1.col + matrix(rep(estimates1, dim), nrow=dim, ncol=length(estimates1), byrow=T)
#学習直後の効果を見るので，記憶率は一定
par_combination[,4] <- estimates1[4]
par_combination[,5] <- estimates1[5]
par_combination[,6] <- estimates1[6]
rate_list <- matrix(rep(0, dim*4), nrow=dim, ncol=4)

for (par_i in 1:dim){
  if (par_i %% 100 == 0){
    print(par_i)
  }
  learning_par <- par_combination[par_i,]
  
  for (i_data in 1:n_data){
    TSNW <- read.csv(paste(inputTSNWpath, TSNW_file[i_data], sep = ""))
    result <- prob(par=par, learning_par = learning_par, tsnw=TSNW, g1_drill=F, g1_cong=T, g2_cong=F, g3_drill=F, g4_map=F)
    rate_list[par_i,] <- rate_list[par_i,] + c(1-result[[1]])
  }
  rate_list[par_i,] <- rate_list[par_i,]/n_data
}
df6[,"rate"] <- apply(rate_list, 2, mean)#列方向に平均をとる
df6[,"low"] <- apply(rate_list, 2, function(x) quantile(x, 0.025))#列方向に2.5%点をとる
df6[,"high"] <- apply(rate_list, 2, function(x) quantile(x, 0.975))#列方向に97.5%点をとる

df8 <- data.frame(timing=c(-48, -24, -12, -6), low = rep(0,4), high = rep(0,4), rate = rep(0,4), group = as.factor("g3_drill"))
TSNW_file <- list.files(inputTSNWpath, pattern="wave1")[which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1))]
n_data <- length(TSNW_file)
hes3.col <- t(chol(learning_hessian3))
#パラメータ推定値が多変量正規分布に従うと仮定して、その分布からパラメータをサンプリング
par_combination <- matrix(rnorm(dim*length(estimates3), 0, 1), nrow = dim, ncol=length(estimates3)) %*% hes3.col + matrix(rep(estimates3, dim), nrow=dim, ncol=length(estimates3), byrow=T)
#学習直後の効果を見るので，記憶率は一定
par_combination[,4] <- estimates3[4]
par_combination[,5] <- estimates3[5]
par_combination[,6] <- estimates3[6]
rate_list <- matrix(rep(0, dim*4), nrow=dim, ncol=4)

for (par_i in 1:dim){
  if (par_i %% 100 == 0){
    print(par_i)
  }
  learning_par <- par_combination[par_i,]
  
  for (i_data in 1:n_data){
    TSNW <- read.csv(paste(inputTSNWpath, TSNW_file[i_data], sep = ""))
    result <- prob(par=par, learning_par = learning_par, tsnw=TSNW, g1_drill=F, g1_cong=F, g2_cong=F, g3_drill=T, g4_map=F)
    rate_list[par_i,] <- rate_list[par_i,] + c(1-result[[1]])
  }
  rate_list[par_i,] <- rate_list[par_i,]/n_data
}
df8[,"rate"] <- apply(rate_list, 2, mean)#列方向に平均をとる
df8[,"low"] <- apply(rate_list, 2, function(x) quantile(x, 0.025))#列方向に2.5%点をとる
df8[,"high"] <- apply(rate_list, 2, function(x) quantile(x, 0.975))#列方向に97.5%点をとる

df9 <- data.frame(timing=c(-48, -24, -12, -6), low = rep(0,4), high = rep(0,4), rate = rep(0,4), group = as.factor("g4_map"))
TSNW_file <- list.files(inputTSNWpath, pattern="wave1")[which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0))]
n_data <- length(TSNW_file)
hes4.col <- t(chol(learning_hessian4))
#パラメータ推定値が多変量正規分布に従うと仮定して、その分布からパラメータをサンプリング
par_combination <- matrix(rnorm(dim*length(estimates4), 0, 1), nrow = dim, ncol=length(estimates4)) %*% hes4.col + matrix(rep(estimates4, dim), nrow=dim, ncol=length(estimates4), byrow=T)
#学習直後の効果を見るので，記憶率は一定
par_combination[,4] <- estimates4[4]
par_combination[,5] <- estimates4[5]
par_combination[,6] <- estimates4[6]
rate_list <- matrix(rep(0, dim*4), nrow=dim, ncol=4)

for (par_i in 1:dim){
  if (par_i %% 100 == 0){
    print(par_i)
  }
  learning_par <- par_combination[par_i,]
  
  for (i_data in 1:n_data){
    TSNW <- read.csv(paste(inputTSNWpath, TSNW_file[i_data], sep = ""))
    result <- prob(par=par, learning_par = learning_par, tsnw=TSNW, g1_drill=F, g1_cong=F, g2_cong=F, g3_drill=F, g4_map=T)
    rate_list[par_i,] <- rate_list[par_i,] + c(1-result[[1]])
  }
  rate_list[par_i,] <- rate_list[par_i,]/n_data
}
df9[,"rate"] <- apply(rate_list, 2, mean)#列方向に平均をとる
df9[,"low"] <- apply(rate_list, 2, function(x) quantile(x, 0.025))#列方向に2.5%点をとる
df9[,"high"] <- apply(rate_list, 2, function(x) quantile(x, 0.975))#列方向に97.5%点をとる




# plot evacuation rate -----------------
df.merged <- rbind(df5, df1)
rate_change <- ggplot(df.merged, aes(x = timing, y = rate, 
                                     color = group, linetype = group, shape = group)) +
  geom_errorbar(aes(ymin = low, ymax = high),  width=1.5, position = position_dodge(1)) +
  geom_line(linewidth=0.7, position = position_dodge(1)) +
  geom_point(size=3, position = position_dodge(1)) +
  geom_text_repel(data = subset(df.merged, group == "g1_without"), aes(label = round(rate,3)), 
                  nudge_x = c(rep(3,4)), nudge_y = c(rep(-0.02,3),-0.03),  size = 6, show.legend = FALSE,  segment.size = 0, segment.color = "white") +
  geom_text_repel(data = subset(df.merged, group == "g1_drill"), aes(label = round(rate,3)), 
                  nudge_x = c(0,rep(-3.5,3)), nudge_y = c(0.05, rep(0.03,3)),  size = 6, show.legend = FALSE,  segment.size = 0, segment.color = "white") +
  geom_hline(yintercept = 0, color = "black", linewidth=1) +  # y軸の原点に線を描く
  geom_vline(xintercept = -50, color = "black",linewidth=1) +  # x軸の原点に線を描く
  scale_x_continuous(limits = c(-50, -1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.75), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-50,-1), ylim = c(0, 0.75)) +
  scale_color_manual(values = c("g1_drill" = "blue", "g1_without" = "black"),
                     labels = c("Evacuation drill", "Without")) +
  scale_linetype_manual(values = c("g1_drill" = "solid", "g1_without" = "dashed"),
                        labels = c("Evacuation drill", "Without")) +
  scale_shape_manual(values = c("g1_drill" = 15, "g1_without" = 16),
                     labels = c("Evacuation drill", "Without")) +
  theme(plot.title = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = c(0.1, 0.9),
        legend.justification = c(0.1, 0.9),
        legend.key.size = unit(1.7, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black",hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 12, color = "black",hjust = 1, vjust = 1),
        axis.line = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Time before typhoon landfall [hour]", y = "Cumulative evacuation probability", title = "Group 1")
rate_change
ggsave(file=paste(project, "/result/picture/rate_change_group1_drill.pdf", sep = ""), plot = rate_change)#, width = 12, height = 12, units = "cm")


df.merged <- rbind(df6, df1)
rate_change <- ggplot(df.merged, aes(x = timing, y = rate, 
                                     color = group, linetype = group, shape = group)) +
  geom_errorbar(aes(ymin = low, ymax = high),  width=1.5, position = position_dodge(1)) +
  geom_line(linewidth=0.7, position = position_dodge(1)) +
  geom_point(size=3, position = position_dodge(1)) +
  geom_text_repel(data = subset(df.merged, group == "g1_without"), aes(label = round(rate,3)), 
                  nudge_x = c(rep(3,4)), nudge_y = c(rep(-0.02,3),-0.03),  size = 6, show.legend = FALSE,  segment.size = 0, segment.color = "white") +
  geom_text_repel(data = subset(df.merged, group == "g1_cong"), aes(label = round(rate,3)), 
                  nudge_x = c(0,rep(-3.5,3)), nudge_y = c(0.05, rep(0.03,3)),  size = 6, show.legend = FALSE,  segment.size = 0, segment.color = "white") +
  geom_hline(yintercept = 0, color = "black", linewidth=1) +  # y軸の原点に線を描く
  geom_vline(xintercept = -50, color = "black",linewidth=1) +  # x軸の原点に線を描く
  scale_x_continuous(limits = c(-50, -1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.75), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-50,-1), ylim = c(0, 0.75)) +
  scale_color_manual(values = c("g1_cong" = "orange", "g1_without" = "black"),
                     labels = c("Info. on others", "Without")) +
  scale_linetype_manual(values = c("g1_cong" = "solid", "g1_without" = "dashed"),
                        labels = c("Info. on others", "Without")) +
  scale_shape_manual(values = c("g1_cong" = 17, "g1_without" = 16),
                     labels = c("Info. on others", "Without")) +
  theme(plot.title = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = c(0.1, 0.9),
        legend.justification = c(0.1, 0.9),
        legend.key.size = unit(1.7, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black",hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 12, color = "black",hjust = 1, vjust = 1),
        axis.line = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Time before typhoon landfall [hour]", y = "Cumulative evacuation probability", title = "Group 1")
rate_change
ggsave(file=paste(project, "/result/picture/rate_change_group1_cong.pdf", sep = ""), plot = rate_change)#, width = 12, height = 12, units = "cm")


df.merged <- rbind(df8, df3)
rate_change <- ggplot(df.merged, aes(x = timing, y = rate, 
                                     color = group, linetype = group, shape = group)) +
  geom_errorbar(aes(ymin = low, ymax = high),  width=1.5, position = position_dodge(1)) +
  geom_line(linewidth=0.7, position = position_dodge(1)) +
  geom_point(size=3, position = position_dodge(1)) +
  geom_text_repel(data = subset(df.merged, group == "g3_without"), aes(label = round(rate,3)), 
                  nudge_x = c(rep(3,4)), nudge_y = c(rep(-0.02,3),-0.03),  size = 6, show.legend = FALSE,  segment.size = 0, segment.color = "white") +
  geom_text_repel(data = subset(df.merged, group == "g3_drill"), aes(label = round(rate,3)), 
                  nudge_x = c(0,rep(-3.5,3)), nudge_y = c(0.05, rep(0.03,3)),  size = 6, show.legend = FALSE,  segment.size = 0, segment.color = "white") +
  geom_hline(yintercept = 0, color = "black", linewidth=1) +  # y軸の原点に線を描く
  geom_vline(xintercept = -50, color = "black",linewidth=1) +  # x軸の原点に線を描く
  scale_x_continuous(limits = c(-50, -1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.75), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-50,-1), ylim = c(0, 0.75)) +
  scale_color_manual(values = c("g3_drill" = "blue", "g3_without" = "black"),
                     labels = c("Evacuation drill", "Without")) +
  scale_linetype_manual(values = c("g3_drill" = "solid", "g3_without" = "dashed"),
                        labels = c("Evacuation drill", "Without")) +
  scale_shape_manual(values = c("g3_drill" = 15, "g3_without" = 16),
                     labels = c("Evacuation drill", "Without")) +
  theme(plot.title = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = c(0.1, 0.9),
        legend.justification = c(0.1, 0.9),
        legend.key.size = unit(1.7, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black",hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 12, color = "black",hjust = 1, vjust = 1),
        axis.line = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Time before typhoon landfall [hour]", y = "Cumulative evacuation probability", title = "Group 3")
rate_change
ggsave(file=paste(project, "/result/picture/rate_change_group3.pdf", sep = ""), plot = rate_change)#, width = 12, height = 12, units = "cm")


df.merged <- rbind(df9, df4)
rate_change <- ggplot(df.merged, aes(x = timing, y = rate, 
                                     color = group, linetype = group, shape = group)) +
  geom_errorbar(aes(ymin = low, ymax = high),  width=1.5, position = position_dodge(1)) +
  geom_line(linewidth=0.7, position = position_dodge(1)) +
  geom_point(size=3, position = position_dodge(1)) +
  geom_text_repel(data = subset(df.merged, group == "g4_without"), aes(label = round(rate,3)), 
                  nudge_x = c(rep(3,4)), nudge_y = c(rep(-0.02,3),-0.03),  size = 6, show.legend = FALSE,  segment.size = 0, segment.color = "white") +
  geom_text_repel(data = subset(df.merged, group == "g4_map"), aes(label = round(rate,3)), 
                  nudge_x = c(0,rep(-3.5,3)), nudge_y = c(0.05, rep(0.03,3)),  size = 6, show.legend = FALSE,  segment.size = 0, segment.color = "white") +
  geom_hline(yintercept = 0, color = "black", linewidth=1) +  # y軸の原点に線を描く
  geom_vline(xintercept = -50, color = "black",linewidth=1) +  # x軸の原点に線を描く
  scale_x_continuous(limits = c(-50, -1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.75), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-50,-1), ylim = c(0, 0.75)) +
  scale_color_manual(values = c("g4_map" = "red", "g4_without" = "black"),
                     labels = c("Hazard map", "Without")) +
  scale_linetype_manual(values = c("g4_map" = "solid", "g4_without" = "dashed"),
                        labels = c("Hazard map", "Without")) +
  scale_shape_manual(values = c("g4_map" = 18, "g4_without" = 16),
                     labels = c("Hazard map", "Without")) +
  theme(plot.title = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = c(0.1, 0.9),
        legend.justification = c(0.1, 0.9),
        legend.key.size = unit(1.7, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black",hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 12, color = "black",hjust = 1, vjust = 1),
        axis.line = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Time before typhoon landfall [hour]", y = "Cumulative evacuation probability", title = "Group 4")
rate_change
ggsave(file=paste(project, "/result/picture/rate_change_group4.pdf", sep = ""), plot = rate_change)#, width = 12, height = 12, units = "cm")


df.merged <- rbind(df1, df3, df4, df5, df6, df8, df9)
write.csv(df.merged, paste(project, "/result/simulation_confidence_interval.csv", sep = ""), row.names = FALSE)
