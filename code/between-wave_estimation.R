### Information learning and forgetting model
### @author:Satoki Masuda
### encoding: UTF-8
### 2022/11/15

rm(list=ls(all=TRUE)) #delte all variables
library(tictoc)
library(Matrix)

createEmptyDf = function(nrow, ncol, colnames = c()){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}

################
### Settings ###
################

#Estimated intial (wave1) parameter
x <- c(-2.227076383, 0.712642215, -1.292325166, 6.33181656, -0.556804972, -0.51096569, -0.169770418, 1.285324338, 2.496594752, 0.917186604, -1.219138708, -1.543252011, 1.0)

project <- "YOUR FOLDER PATH"
name <- "1_safe_evac"
inputpath <- paste(project, "/input", sep = "")

inputlinkpath <- paste(inputpath, "/TSNW/", sep="")
inputroutepath <- paste(inputpath, "/route/user/", sep="")

personal <- read.csv(paste(inputpath, "/PersonalInfo.csv", sep=""), fileEncoding = "Shift-jis")

setwd(project)

####################
### data reading ###
####################
#route data
route_file2 <- list.files(inputroutepath, pattern="wave2")[which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1))]
route_file3 <- list.files(inputroutepath, pattern="wave3")[which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1))]
route_file4 <- list.files(inputroutepath, pattern="wave4")[which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1))]
#link data
link_file2 <- list.files(inputlinkpath, pattern="wave2")[which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1))]
link_file3 <- list.files(inputlinkpath, pattern="wave3")[which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1))]
link_file4 <- list.files(inputlinkpath, pattern="wave4")[which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1))]

n_data2 <- length(route_file2)
n_data3 <- length(route_file3)
n_data4 <- length(route_file4)


#######################
### estimation part ###
#######################

#wave2-wave4の尤度関数．wave1の推定パラメータに変化分を足し算．追加された報酬関数のパラメータを求める．
fr_LL <- function(r2){
  LL <- 0
  print(r2)
  theta <- x[length(x)] #discount factor
  
  for (i_data in 1:n_data2){
    if (is.na(LL)){
      print(route_file2[i_data - 1])
      break
    }
    
    link <- read.csv(paste(inputlinkpath, link_file2[i_data], sep = ""))
    route <- read.csv(paste(inputroutepath, route_file2[i_data], sep = ""))
    
    #node data
    node_space <- sort(unique(link$a_space))
    n <- length(node_space)
    node_TSN <- c(sort(unique(c(link$k, link$a))), 99999) #absorbing state
    N <- length(node_TSN)
    
    #link data
    absLink <- createEmptyDf(nrow = n, ncol = ncol(link), colnames = colnames(link)) # absorbing link
    absLink[,] <- 0
    absLink$linkID <- 1:n
    absLink$k <- node_space + 40000
    absLink$a <- 99999
    absLink$k_space <- node_space
    absLink$a_space <- 99999
    link <- rbind(link, absLink)
    link$linkID <- 1:nrow(link)
    L <- nrow(link)
    
    #adjacent matrix I of time-structured NW
    I <- array(0, dim = c(N, N))
    for(i in 1:nrow(link)){
      k_i <- (1:N)[node_TSN == link$k[i]]
      a_i <- (1:N)[node_TSN == link$a[i]]
      I[k_i,a_i] <- 1
    }
    
    
    #################
    ### Variables ###
    #################
    #calculate utility with exponential form (note!)
    
    link[link$distance>0, "distance"] <- log(10*link[link$distance>0, "distance"])/10 #距離は0.1-900kmくらいの間でばらついてる
    
    #for value function of the last node
    #depth * non-home dummy
    dst_depth <- array(exp(0.0), dim = c(N,N))
    #depth * hazardmap * non-home dummy
    dst_depth_map <- dst_depth
    #school * non-home dummy
    dst_school <- dst_depth
    
    #initial subj. depth of home * home dummy
    dst_home_sbj <- dst_depth
    #living in inundation area dummy * under 2nd floor dummy * home dummy
    dst_home_low <- dst_depth
    #care dummy * home dummy
    dst_home_care <- dst_depth
    #ASC
    dst_home <- dst_depth
    #interest
    dst_interest <- dst_depth
    
    link[, "depth"] <- link[, "depth"]/10
    link[, "home_depth"] <- link[, "home_depth"]/10
    link[, "init_subj_depth"] <- link[, "init_subj_depth"]/10
    link[, "subj_prob1"] <- link[, "subj_prob1"]/100
    
    #evacuation order * home dummy
    order <- dst_depth
    
    #evacuation order * home dummy
    order_low <- dst_depth
    
    #distance * non-home dummy
    distance_1 <- dst_depth
    distance_2 <- dst_depth
    distance_3 <- dst_depth
    distance_4 <- dst_depth
    
    #reward
    congestion2 <- dst_depth
    dst_congestion2 <- dst_depth
    drill2 <- dst_depth
    #toyosu2 <- dst_depth
    subj_prob1 <- dst_depth
    prev <- dst_depth
    
    for (i in 1:L){
      k_i <- (1:N)[node_TSN == link$k[i]]
      a_i <- (1:N)[node_TSN == link$a[i]]
      
      if (link$a[i]==99999){
        if (link[link$a==link$k[i], "home"][1]==0){
          dst_depth[k_i,a_i] <- exp(link[link$a==link$k[i], "depth"][1])
          dst_school[k_i,a_i] <- exp((link[link$a==link$k[i], "a_type"][1]=="primary_school") | (link[link$a==link$k[i], "a_type"][1]=="school"))
          dst_congestion2[k_i, a_i] <- exp(sum(link[link$a_space==link$k_space[i], "congestion2"])>0)
        }
        else{
          dst_home_low[k_i, a_i] <- exp((link[1, "home_depth"]>0) * (link[1, "height"]<=2))
          dst_home[k_i, a_i] <- exp(1)
          dst_interest[k_i, a_i] <- exp(6 - link[1, "bousai_interest"])
          drill2[k_i,a_i] <- exp(link[1, "drill2"])
          }
      }
      
      order[k_i,a_i] <- exp((link[link$a_space==link$k_space[i], "home"][1]) * (link$home[i]==0) *
                              ((link$koiki_order[i]) + (link$level4[i])))
      order_low[k_i,a_i] <- exp((link[link$a_space==link$k_space[i], "home"][1]) * (link$home[i]==0) *
                                  (link$home_depth[i]>0) * (link$height[i]<=2) *
                                  ((link$koiki_order[i]) + (link$level4[i])))
      
      distance_1[k_i,a_i] <- exp(link$distance[i] * (link$home[i]==0) * (link$a[i]<20000))
      distance_2[k_i,a_i] <- exp(link$distance[i] * (link$home[i]==0) * ((link$a[i]>20000)&(link$a[i]<30000)))
      distance_3[k_i,a_i] <- exp(link$distance[i] * (link$home[i]==0) * ((link$a[i]>30000)&(link$a[i]<40000)))
      distance_4[k_i,a_i] <- exp(link$distance[i] * (link$home[i]==0) * ((link$a[i]>40000)&(link$a[i]<50000)))
      
      subj_prob1[k_i,a_i] <- exp((link$subj_prob1[i]) * (link$home[i]==1))
      #主観確率はwave1のものを使う。wave2以降に学習により主観確率が変化して行動が変容するのは学習効果に含むこととする。
      prev[k_i,a_i] <- exp(link$prev[i])
    }
    
    
    ###############################
    ### log-likelihood function ###
    ###############################
    
    # utility function matrix M for each class
    M_1 <- array(0, dim = c(N, N))
    M_1 <- I*
      (dst_depth^x[1])*(dst_school^x[2])*
      (dst_home_low^x[3])*(dst_home^x[4])*(dst_interest^x[5])*(subj_prob1^x[6])*
      (order^x[7])*(order_low^x[8])*
      (distance_1^x[9])*(distance_2^x[10])*(distance_3^x[11])*(distance_4^x[12])*
      ((drill2^r2[2])^(r2[5])^1)*
      ((dst_congestion2^r2[3])^(r2[6]^1))*
      (prev^r2[7])
    
    # value function V=log(z)
    # calculate z(t)=Mz(t+1)+b from t=T until t=0 backward
    
    #exp(Vd). value function of the absorbing step for each class
    z_1 <- array(1, dim = c(1, N))
    
    N_ts <- c(sum(node_TSN<10000), #num of nodes at time step 0 (home)
              sum((node_TSN>10000)&(node_TSN<20000)), #num of nodes at time step 1
              sum((node_TSN>20000)&(node_TSN<30000)), #num of nodes at time step 2
              sum((node_TSN>30000)&(node_TSN<40000)), #num of nodes at time step 3
              sum((node_TSN>40000)&(node_TSN<50000)), #num of nodes at time step 4
              sum(node_TSN>90000)) #num of absorbing nodes
    
    idx_1 <- N - N_ts[6] + 1
    idx_2 <- N
    
    li_theta <- c(1, theta^4, theta^2, theta, 1)
    
    for (i in 5:1){ #num of step (4) + time step 0 (1)
      z_i1 <- matrix(M_1[(idx_1 - N_ts[i]):(idx_2 - N_ts[i+1]), (idx_1):(idx_2)], c(N_ts[i], N_ts[i+1])) %*% (z_1[(idx_1):(idx_2)] ^(li_theta[i]))
      idx_1 <- idx_1 - N_ts[i]
      idx_2 <- idx_2 - N_ts[i+1]
      z_1[idx_1:idx_2] <- (z_i1 == 0)*1 + (z_i1 != 0)*z_i1
    }
    
    V_1 <- z_1 #exp(V^d) for calculation of denominator of choice probability
    
    #calculate log-likelihood
    k_n <- route[,"k"]
    a_n <- route[,"a"]
    hhd <- length(k_n) #num of links in the route
    
    for (h in 1:hhd){
      k <- (1:N)[node_TSN == k_n[h]]
      a <- (1:N)[node_TSN == a_n[h]]
      LL <- LL + log((M_1[k,a] * V_1[a]^(li_theta[h])) / V_1[k])
    }
  }
  
  
  for (i_data in 1:n_data3){
    if (is.na(LL)){
      print(route_file3[i_data - 1])
      break
    }
    
    link <- read.csv(paste(inputlinkpath, link_file3[i_data], sep = ""))
    route <- read.csv(paste(inputroutepath, route_file3[i_data], sep = ""))
    
    #node data
    node_space <- sort(unique(link$a_space))
    n <- length(node_space)
    node_TSN <- c(sort(unique(c(link$k, link$a))), 99999) #absorbing state
    N <- length(node_TSN)
    
    #link data
    absLink <- createEmptyDf(nrow = n, ncol = ncol(link), colnames = colnames(link)) # absorbing link
    absLink[,] <- 0
    absLink$linkID <- 1:n
    absLink$k <- node_space + 40000
    absLink$a <- 99999
    absLink$k_space <- node_space
    absLink$a_space <- 99999
    link <- rbind(link, absLink)
    link$linkID <- 1:nrow(link)
    L <- nrow(link)
    
    #adjacent matrix I of time-structured NW
    I <- array(0, dim = c(N, N))
    for(i in 1:nrow(link)){
      k_i <- (1:N)[node_TSN == link$k[i]]
      a_i <- (1:N)[node_TSN == link$a[i]]
      I[k_i,a_i] <- 1
    }
    
    
    #################
    ### Variables ###
    #################
    #calculate utility with exponential form (note!)
    
    link[link$distance>0, "distance"] <- log(10*link[link$distance>0, "distance"])/10 #距離は0.1-900kmくらいの間でばらついてる
    
    #for value function of the last node
    #depth * non-home dummy
    dst_depth <- array(exp(0.0), dim = c(N,N))
    #depth * hazardmap * non-home dummy
    dst_depth_map <- dst_depth
    #school * non-home dummy
    dst_school <- dst_depth
    
    #initial subj. depth of home * home dummy
    dst_home_sbj <- dst_depth
    #living in inundation area dummy * under 2nd floor dummy * home dummy
    dst_home_low <- dst_depth
    #care dummy * home dummy
    dst_home_care <- dst_depth
    #ASC
    dst_home <- dst_depth
    #interest
    dst_interest <- dst_depth
    
    link[, "depth"] <- link[, "depth"]/10
    link[, "home_depth"] <- link[, "home_depth"]/10
    link[, "init_subj_depth"] <- link[, "init_subj_depth"]/10
    link[, "subj_prob1"] <- link[, "subj_prob1"]/100
    
    #evacuation order * home dummy
    order <- dst_depth
    #evacuation order * home dummy
    order_low <- dst_depth
    
    #distance * non-home dummy
    distance_1 <- dst_depth
    distance_2 <- dst_depth
    distance_3 <- dst_depth
    distance_4 <- dst_depth
    
    #reward
    congestion2 <- dst_depth
    dst_congestion2 <- dst_depth
    drill2 <- dst_depth
    #toyosu2 <- dst_depth
    congestion3 <- dst_depth
    dst_congestion3 <- dst_depth
    flood_time3 <- dst_depth
    drill3 <- dst_depth
    #toyosu3 <- dst_depth
    subj_prob1 <- dst_depth
    prev <- dst_depth
    
    for (i in 1:L){
      k_i <- (1:N)[node_TSN == link$k[i]]
      a_i <- (1:N)[node_TSN == link$a[i]]
      
      if (link$a[i]==99999){
        if (link[link$a==link$k[i], "home"][1]==0){
          dst_depth[k_i,a_i] <- exp(link[link$a==link$k[i], "depth"][1])
          dst_school[k_i,a_i] <- exp((link[link$a==link$k[i], "a_type"][1]=="primary_school") | (link[link$a==link$k[i], "a_type"][1]=="school"))
          dst_congestion2[k_i, a_i] <- exp(sum(link[link$a_space==link$k_space[i], "congestion2"])>0)
          dst_congestion3[k_i, a_i] <- exp(sum(link[link$a_space==link$k_space[i], "congestion3"])>0)
        }
        else{
          dst_home_low[k_i, a_i] <- exp((link[1, "home_depth"]>0) * (link[1, "height"]<=2))
          dst_home[k_i, a_i] <- exp(1)
          dst_interest[k_i, a_i] <- exp(6 - link[1, "bousai_interest"])
          drill2[k_i,a_i] <- exp(link[1, "drill2"])
          drill3[k_i,a_i] <- exp(link[1, "drill3"])
        }
        #if flood info is given, flood_time=1
        flood_time3[k_i,a_i] <- exp(link[1, "flood_info3"]  * (link[(link$a==link$k[i]), "duration"][1]>0))
      }
      
      order[k_i,a_i] <- exp((link[link$a_space==link$k_space[i], "home"][1]) * (link$home[i]==0) *
                              ((link$koiki_order[i]) + (link$level4[i])))
      order_low[k_i,a_i] <- exp((link[link$a_space==link$k_space[i], "home"][1]) * (link$home[i]==0) *
                                  (link$home_depth[i]>0) * (link$height[i]<=2) *
                                  ((link$koiki_order[i]) + (link$level4[i])))
      
      distance_1[k_i,a_i] <- exp(link$distance[i] * (link$home[i]==0) * (link$a[i]<20000))
      distance_2[k_i,a_i] <- exp(link$distance[i] * (link$home[i]==0) * ((link$a[i]>20000)&(link$a[i]<30000)))
      distance_3[k_i,a_i] <- exp(link$distance[i] * (link$home[i]==0) * ((link$a[i]>30000)&(link$a[i]<40000)))
      distance_4[k_i,a_i] <- exp(link$distance[i] * (link$home[i]==0) * ((link$a[i]>40000)&(link$a[i]<50000)))
      
      subj_prob1[k_i,a_i] <- exp((link$subj_prob1[i]) * (link$home[i]==1))
      prev[k_i,a_i] <- exp(link$prev[i])
    }
    
    
    ###############################
    ### log-likelihood function ###
    ###############################
    
    # utility function matrix M for each class
    M_1 <- array(0, dim = c(N, N))
    M_1 <- I*
      (dst_depth^x[1])*(dst_school^x[2])*
      (dst_home_low^x[3])*(dst_home^x[4])*(dst_interest^x[5])*(subj_prob1^x[6])*
      (order^x[7])*(order_low^x[8])*
      (distance_1^x[9])*(distance_2^x[10])*(distance_3^x[11])*(distance_4^x[12])*
      (prev^r2[7])*
      ((flood_time3^r2[1])^(r2[4]^1))*
      ((drill3^r2[2])^(r2[5]^1))*
      ((dst_congestion3^r2[3])^(r2[6]^1))*
         ((drill2^r2[2])^(r2[5]^4))*
         ((dst_congestion2^r2[3])^(r2[6]^4))
    
    # value function V=log(z)
    # calculate z(t)=Mz(t+1)+b from t=T until t=0 backward
    
    #exp(Vd). value function of the absorbing step for each class
    z_1 <- array(1, dim = c(1, N))
    
    N_ts <- c(sum(node_TSN<10000), #num of nodes at time step 0 (home)
              sum((node_TSN>10000)&(node_TSN<20000)), #num of nodes at time step 1
              sum((node_TSN>20000)&(node_TSN<30000)), #num of nodes at time step 2
              sum((node_TSN>30000)&(node_TSN<40000)), #num of nodes at time step 3
              sum((node_TSN>40000)&(node_TSN<50000)), #num of nodes at time step 4
              sum(node_TSN>90000)) #num of absorbing nodes
    
    idx_1 <- N - N_ts[6] + 1
    idx_2 <- N
    
    li_theta <- c(1, theta^4, theta^2, theta, 1)
    
    for (i in 5:1){ #num of step (4) + time step 0 (1)
      z_i1 <- matrix(M_1[(idx_1 - N_ts[i]):(idx_2 - N_ts[i+1]), (idx_1):(idx_2)], c(N_ts[i], N_ts[i+1])) %*% (z_1[(idx_1):(idx_2)] ^(li_theta[i]))
      idx_1 <- idx_1 - N_ts[i]
      idx_2 <- idx_2 - N_ts[i+1]
      z_1[idx_1:idx_2] <- (z_i1 == 0)*1 + (z_i1 != 0)*z_i1
    }
    
    V_1 <- z_1 #exp(V^d) for calculation of denominator of choice probability
    
    #calculate log-likelihood
    k_n <- route[,"k"]
    a_n <- route[,"a"]
    hhd <- length(k_n) #num of links in the route
    
    for (h in 1:hhd){
      k <- (1:N)[node_TSN == k_n[h]]
      a <- (1:N)[node_TSN == a_n[h]]
      LL <- LL + log((M_1[k,a] * V_1[a]^(li_theta[h])) / V_1[k])
    }
  }
  
  
  for (i_data in 1:n_data4){
    if (is.na(LL)){
      print(route_file4[i_data - 1])
      break
    }
    
    link <- read.csv(paste(inputlinkpath, link_file4[i_data], sep = ""))
    route <- read.csv(paste(inputroutepath, route_file4[i_data], sep = ""))
    
    #node data
    node_space <- sort(unique(link$a_space))
    n <- length(node_space)
    node_TSN <- c(sort(unique(c(link$k, link$a))), 99999) #absorbing state
    N <- length(node_TSN)
    
    #link data
    absLink <- createEmptyDf(nrow = n, ncol = ncol(link), colnames = colnames(link)) # absorbing link
    absLink[,] <- 0
    absLink$linkID <- 1:n
    absLink$k <- node_space + 40000
    absLink$a <- 99999
    absLink$k_space <- node_space
    absLink$a_space <- 99999
    link <- rbind(link, absLink)
    link$linkID <- 1:nrow(link)
    L <- nrow(link)
    
    #adjacent matrix I of time-structured NW
    I <- array(0, dim = c(N, N))
    for(i in 1:nrow(link)){
      k_i <- (1:N)[node_TSN == link$k[i]]
      a_i <- (1:N)[node_TSN == link$a[i]]
      I[k_i,a_i] <- 1
    }
    
    
    #################
    ### Variables ###
    #################
    #calculate utility with exponential form (note!)
    
    link[link$distance>0, "distance"] <- log(10*link[link$distance>0, "distance"])/10 #距離は0.1-900kmくらいの間でばらついてる
    
    #for value function of the last node
    #depth * non-home dummy
    dst_depth <- array(exp(0.0), dim = c(N,N))
    #depth * hazardmap * non-home dummy
    dst_depth_map <- dst_depth
    #school * non-home dummy
    dst_school <- dst_depth
    
    #initial subj. depth of home * home dummy
    dst_home_sbj <- dst_depth
    #living in inundation area dummy * under 2nd floor dummy * home dummy
    dst_home_low <- dst_depth
    #care dummy * home dummy
    dst_home_care <- dst_depth
    #ASC
    dst_home <- dst_depth
    #interest
    dst_interest <- dst_depth
    
    link[, "depth"] <- link[, "depth"]/10
    link[, "home_depth"] <- link[, "home_depth"]/10
    link[, "init_subj_depth"] <- link[, "init_subj_depth"]/10
    link[, "subj_prob1"] <- link[, "subj_prob1"]/100
    
    #evacuation order * home dummy
    order <- dst_depth
    #evacuation order * home dummy
    order_low <- dst_depth
    
    #distance * non-home dummy
    distance_1 <- dst_depth
    distance_2 <- dst_depth
    distance_3 <- dst_depth
    distance_4 <- dst_depth
    
    #reward
    congestion2 <- dst_depth
    dst_congestion2 <- dst_depth
    drill2 <- dst_depth
    #toyosu2 <- dst_depth
    congestion3 <- dst_depth
    dst_congestion3 <- dst_depth
    flood_time3 <- dst_depth
    drill3 <- dst_depth
    #toyosu3 <- dst_depth
    congestion4 <- dst_depth
    dst_congestion4 <- dst_depth
    flood_time4 <- dst_depth
    drill4 <- dst_depth
    #toyosu4 <- dst_depth
    subj_prob1 <- dst_depth
    prev <- dst_depth
    
    for (i in 1:L){
      k_i <- (1:N)[node_TSN == link$k[i]]
      a_i <- (1:N)[node_TSN == link$a[i]]
      
      if (link$a[i]==99999){
        if (link[link$a==link$k[i], "home"][1]==0){
          dst_depth[k_i,a_i] <- exp(link[link$a==link$k[i], "depth"][1])
          dst_school[k_i,a_i] <- exp((link[link$a==link$k[i], "a_type"][1]=="primary_school") | (link[link$a==link$k[i], "a_type"][1]=="school"))
          dst_congestion2[k_i, a_i] <- exp(sum(link[link$a_space==link$k_space[i], "congestion2"])>0)
          dst_congestion3[k_i, a_i] <- exp(sum(link[link$a_space==link$k_space[i], "congestion3"])>0)
          dst_congestion4[k_i, a_i] <- exp(sum(link[link$a_space==link$k_space[i], "congestion4"])>0)
        }
        else{
          dst_home_low[k_i, a_i] <- exp((link[1, "home_depth"]>0) * (link[1, "height"]<=2))
          dst_home[k_i, a_i] <- exp(1)
          dst_interest[k_i, a_i] <- exp(6 - link[1, "bousai_interest"])
          drill2[k_i,a_i] <- exp(link[1, "drill2"])
          drill3[k_i,a_i] <- exp(link[1, "drill3"])
          drill4[k_i,a_i] <- exp(link[1, "drill4"])
        }
        #if flood info is given, flood_time=1
        flood_time3[k_i,a_i] <- exp(link[1, "flood_info3"]  * (link[(link$a==link$k[i]), "duration"][1]>0))
        flood_time4[k_i,a_i] <- exp(link[1, "flood_info4"]  * (link[(link$a==link$k[i]), "duration"][1]>0))
      }
      
      order[k_i,a_i] <- exp((link[link$a_space==link$k_space[i], "home"][1]) * (link$home[i]==0) *
                              ((link$koiki_order[i]) + (link$level4[i])))
      order_low[k_i,a_i] <- exp((link[link$a_space==link$k_space[i], "home"][1]) * (link$home[i]==0) *
                                  (link$home_depth[i]>0) * (link$height[i]<=2) *
                                  ((link$koiki_order[i]) + (link$level4[i])))
      
      distance_1[k_i,a_i] <- exp(link$distance[i] * (link$home[i]==0) * (link$a[i]<20000))
      distance_2[k_i,a_i] <- exp(link$distance[i] * (link$home[i]==0) * ((link$a[i]>20000)&(link$a[i]<30000)))
      distance_3[k_i,a_i] <- exp(link$distance[i] * (link$home[i]==0) * ((link$a[i]>30000)&(link$a[i]<40000)))
      distance_4[k_i,a_i] <- exp(link$distance[i] * (link$home[i]==0) * ((link$a[i]>40000)&(link$a[i]<50000)))
      subj_prob1[k_i,a_i] <- exp((link$subj_prob1[i]) * (link$home[i]==1))
      prev[k_i,a_i] <- exp(link$prev[i])
    }
    
    
    ###############################
    ### log-likelihood function ###
    ###############################
    
    # utility function matrix M for each class
    M_1 <- array(0, dim = c(N, N))
    M_1 <- I*
      (dst_depth^x[1])*(dst_school^x[2])*
      (dst_home_low^x[3])*(dst_home^x[4])*(dst_interest^x[5])*(subj_prob1^x[6])*
      (order^x[7])*(order_low^x[8])*
      (distance_1^x[9])*(distance_2^x[10])*(distance_3^x[11])*(distance_4^x[12])*
      (prev^r2[7])*
      ((flood_time4^r2[1])^(r2[4]^2))*
      ((drill4^r2[2])^(r2[5]^2))*
      ((dst_congestion4^r2[3])^(r2[6]^2))*
         ((flood_time3^r2[1])^(r2[4]^5))*
         ((drill3^r2[2])^(r2[5]^5))*
         ((dst_congestion3^r2[3])^(r2[6]^5))*
            ((drill2^r2[2])^(r2[5]^8))*
            ((dst_congestion2^r2[3])^(r2[6]^8))
    
    # value function V=log(z)
    # calculate z(t)=Mz(t+1)+b from t=T until t=0 backward
    
    #exp(Vd). value function of the absorbing step for each class
    z_1 <- array(1, dim = c(1, N))
    
    N_ts <- c(sum(node_TSN<10000), #num of nodes at time step 0 (home)
              sum((node_TSN>10000)&(node_TSN<20000)), #num of nodes at time step 1
              sum((node_TSN>20000)&(node_TSN<30000)), #num of nodes at time step 2
              sum((node_TSN>30000)&(node_TSN<40000)), #num of nodes at time step 3
              sum((node_TSN>40000)&(node_TSN<50000)), #num of nodes at time step 4
              sum(node_TSN>90000)) #num of absorbing nodes
    
    idx_1 <- N - N_ts[6] + 1
    idx_2 <- N
    
    li_theta <- c(1, theta^4, theta^2, theta, 1)
    
    for (i in 5:1){ #num of step (4) + time step 0 (1)
      z_i1 <- matrix(M_1[(idx_1 - N_ts[i]):(idx_2 - N_ts[i+1]), (idx_1):(idx_2)], c(N_ts[i], N_ts[i+1])) %*% (z_1[(idx_1):(idx_2)] ^(li_theta[i]))
      idx_1 <- idx_1 - N_ts[i]
      idx_2 <- idx_2 - N_ts[i+1]
      z_1[idx_1:idx_2] <- (z_i1 == 0)*1 + (z_i1 != 0)*z_i1
    }
    
    V_1 <- z_1 #exp(V^d) for calculation of denominator of choice probability
    
    #calculate log-likelihood
    k_n <- route[,"k"]
    a_n <- route[,"a"]
    hhd <- length(k_n) #num of links in the route
    
    for (h in 1:hhd){
      k <- (1:N)[node_TSN == k_n[h]]
      a <- (1:N)[node_TSN == a_n[h]]
      LL <- LL + log((M_1[k,a] * V_1[a]^(li_theta[h])) / V_1[k])
    }
  }
  
  
  print(LL)
  return(LL)
}




# Estimation parameter r2
binit <- c(numeric(3),1,1,1,0)

# -------------- estimation start! -------------- #
tic()

lower <-  c(rep(-Inf, 3),rep(0,3), -Inf)
upper <-  c(rep(Inf, 3),rep(1,3), Inf)
res <- optim(binit, fr_LL, method = "L-BFGS-B", lower = lower, upper = upper, hessian = TRUE, control=list(fnscale=-1))
b <- res$par
hhh <- res$hessian
tval <- b/sqrt(abs(diag(solve(hhh)))) #-(solve(hhh))が分散共分散行列
LL <- res$value
Lc <- fr_LL(binit)

# 結果の表示
print(res)
print(Lc) 
print(LL)
print((Lc-LL)/Lc) 
print((Lc-(LL-length(b)))/Lc) 
print(b)
print(tval)

# -------------- estimation finished! -------------- #
toc()


params <- data.frame(estimates=b,t.value=tval,initial_likelihood=c(Lc,rep("",length(b)-1)), final_likelihood=c(LL,rep("",length(b)-1)), likelihood_ratio=c((Lc-LL)/Lc,rep("",length(b)-1)), fixed_likelihood_ratio=c((Lc-(LL-length(b)))/Lc,rep("",length(b)-1)))
write.csv(params, paste(project, "/result/",name,".csv", sep = ""))
#分散共分散行列行列
write.csv(-(solve(hhh)), paste(project, "/result/hessian",name,".csv", sep = ""), row.names=FALSE)
