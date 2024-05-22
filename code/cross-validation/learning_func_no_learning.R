### Information learning model
### @author:Satoki Masuda
### encoding: UTF-8
### 2024/5/2

library(tictoc)
library(Matrix)

createEmptyDf = function(nrow, ncol, colnames = c()){
  data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
}

learning_model_estimation <- function(route_file, link_file){
  
  #######################
  ### estimation part ###
  #######################
  n_data <- length(route_file)
  
  #wave2-wave4の尤度関数．wave1の推定パラメータに変化分を足し算．追加された報酬関数のパラメータを求める．
  fr <- function(x){
    LL <- 0
    #print(x)
    theta <- 1.0
    #theta <- x[length(x)] #discount factor
    
    for (i_data in 1:n_data){
      if (is.na(LL)){
        print(route_file[i_data - 1])
        break
      }
      
      link <- read.csv(paste(inputlinkpath, link_file[i_data], sep = ""))
      route <- read.csv(paste(inputroutepath, route_file[i_data], sep = ""))
      
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
      #subjective probability
      subj_prob <- dst_depth
      
      link[, "depth"] <- link[, "depth"]/10
      link[, "home_depth"] <- link[, "home_depth"]/10
      link[, "init_subj_depth"] <- link[, "init_subj_depth"]/10
      link[, "subj_prob"] <- link[, "subj_prob"]/100
      
      #evacuation order * non-home dummy
      order <- dst_depth
      #evacuation order * low * non-home dummy
      order_low <- dst_depth
      
      distance_1 <- dst_depth
      distance_2 <- dst_depth
      distance_3 <- dst_depth
      distance_4 <- dst_depth
      
      for (i in 1:L){
        k_i <- (1:N)[node_TSN == link$k[i]]
        a_i <- (1:N)[node_TSN == link$a[i]]
        
        if (link$a[i]==99999){
          if (link[link$a==link$k[i], "home"][1]==0){
            dst_depth[k_i,a_i] <- exp(link[link$a==link$k[i], "depth"][1])
            dst_school[k_i,a_i] <- exp((link[link$a==link$k[i], "a_type"][1]=="primary_school") | (link[link$a==link$k[i], "a_type"][1]=="school"))
          }
          else{
            dst_home_low[k_i, a_i] <- exp((link[1, "home_depth"]>0) * (link[1, "height"]<=2))
            dst_home[k_i, a_i] <- exp(1)
            dst_interest[k_i, a_i] <- exp(6 - link[1, "bousai_interest"])
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
        subj_prob[k_i,a_i] <- exp((link$subj_prob[i]) * (link$home[i]==1))
      }
      
      
      ###############################
      ### log-likelihood function ###
      ###############################
      # utility function matrix M
      M <- array(0, dim = c(N, N))
      #for all data
      M <- I*
        (dst_depth^x[1])*(dst_school^x[2])*
        (dst_home_low^x[3])*(dst_home^x[4])*(dst_interest^x[5]) *(subj_prob^x[6])*
        (order^x[7])*(order_low^x[8])*
        (distance_1^x[9])*(distance_2^x[10])*(distance_3^x[11])*(distance_4^x[12])
      
      # value function V=log(z)
      # calculate z(t)=Mz(t+1)+b from t=T until t=0 backward
      
      z <- array(1, dim = c(1, N)) #exp(Vd). value function of the absorbing step
      
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
        z_i <- matrix(M[(idx_1 - N_ts[i]):(idx_2 - N_ts[i+1]), (idx_1):(idx_2)], c(N_ts[i], N_ts[i+1])) %*% (z[(idx_1):(idx_2)] ^(li_theta[i]))
        idx_1 <- idx_1 - N_ts[i]
        idx_2 <- idx_2 - N_ts[i+1]
        z[idx_1:idx_2] <- (z_i == 0)*1 + (z_i != 0)*z_i
      } 
      
      V <- log(z) #delete exp()
      
      #calculate log-likelihood
      k_n <- route[,"k"]
      a_n <- route[,"a"]
      hhd <- length(k_n) #num of links in the route
      
      for (h in 1:hhd){
        k <- (1:N)[node_TSN == k_n[h]]
        a <- (1:N)[node_TSN == a_n[h]]
        LL <- LL + (log(M[k,a]) + (li_theta[h])*V[a] - V[k])
      }
    }
    #print(LL)
    return(LL)
  }
  
  
  # -------------- estimation start! -------------- #
  
  binit <- numeric(12)
  lower <-  c(rep(-Inf, length(binit)))
  upper <-  c(rep(Inf, length(binit)))
  # lower <-  c(rep(-Inf, length(binit)-1),0)
  # upper <-  c(rep(Inf, length(binit)-1),1)
  res <- optim(binit, fr, method = "L-BFGS-B", lower = lower, upper = upper, hessian = TRUE, control=list(fnscale=-1))
  #res <- optim(binit, fr, method = "L-BFGS-B", lower = c(-10, -10, 0), upper = c(10, 10, 1), hessian = TRUE, control=list(fnscale=-1))
  b <- res$par
  # hhh <- res$hessian
  # tval <- b/sqrt(abs(diag(solve(hhh))))
  # LL <- res$value
  # Lc <- fr_LL(binit)
  
  return(b)
}





learning_model_prediction <- function(route_file, link_file, param){
  n_data <- length(route_file)
  n_sample <- n_data/3 #1人につき2-4の3wave分あるので
  
  #wave2-wave4の尤度関数．wave1の推定パラメータに変化分を足し算．追加された報酬関数のパラメータを求める．
  fr <- function(x){
    LL <- 0
    #print(x)
    theta <- 1.0#discount factor
    
    for (i_data in 1:n_data){
      if (is.na(LL)){
        print(route_file[i_data - 1])
        break
      }
      
      link <- read.csv(paste(inputlinkpath, link_file[i_data], sep = ""))
      route <- read.csv(paste(inputroutepath, route_file[i_data], sep = ""))
      
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
      #subjective probability
      subj_prob <- dst_depth
      
      link[, "depth"] <- link[, "depth"]/10
      link[, "home_depth"] <- link[, "home_depth"]/10
      link[, "init_subj_depth"] <- link[, "init_subj_depth"]/10
      link[, "subj_prob"] <- link[, "subj_prob"]/100
      
      #evacuation order * non-home dummy
      order <- dst_depth
      #evacuation order * low * non-home dummy
      order_low <- dst_depth
      
      distance_1 <- dst_depth
      distance_2 <- dst_depth
      distance_3 <- dst_depth
      distance_4 <- dst_depth
      
      for (i in 1:L){
        k_i <- (1:N)[node_TSN == link$k[i]]
        a_i <- (1:N)[node_TSN == link$a[i]]
        
        if (link$a[i]==99999){
          if (link[link$a==link$k[i], "home"][1]==0){
            dst_depth[k_i,a_i] <- exp(link[link$a==link$k[i], "depth"][1])
            dst_school[k_i,a_i] <- exp((link[link$a==link$k[i], "a_type"][1]=="primary_school") | (link[link$a==link$k[i], "a_type"][1]=="school"))
          }
          else{
            dst_home_low[k_i, a_i] <- exp((link[1, "home_depth"]>0) * (link[1, "height"]<=2))
            dst_home[k_i, a_i] <- exp(1)
            dst_interest[k_i, a_i] <- exp(6 - link[1, "bousai_interest"])
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
        subj_prob[k_i,a_i] <- exp((link$subj_prob[i]) * (link$home[i]==1))
      }
      
      
      ###############################
      ### log-likelihood function ###
      ###############################
      # utility function matrix M
      M <- array(0, dim = c(N, N))
      #for all data
      M <- I*
        (dst_depth^x[1])*(dst_school^x[2])*
        (dst_home_low^x[3])*(dst_home^x[4])*(dst_interest^x[5]) *(subj_prob^x[6])*
        (order^x[7])*(order_low^x[8])*
        (distance_1^x[9])*(distance_2^x[10])*(distance_3^x[11])*(distance_4^x[12])
      
      # value function V=log(z)
      # calculate z(t)=Mz(t+1)+b from t=T until t=0 backward
      
      z <- array(1, dim = c(1, N)) #exp(Vd). value function of the absorbing step
      
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
        z_i <- matrix(M[(idx_1 - N_ts[i]):(idx_2 - N_ts[i+1]), (idx_1):(idx_2)], c(N_ts[i], N_ts[i+1])) %*% (z[(idx_1):(idx_2)] ^(li_theta[i]))
        idx_1 <- idx_1 - N_ts[i]
        idx_2 <- idx_2 - N_ts[i+1]
        z[idx_1:idx_2] <- (z_i == 0)*1 + (z_i != 0)*z_i
      } 
      
      V <- log(z) #delete exp()
      
      #calculate log-likelihood
      k_n <- route[,"k"]
      a_n <- route[,"a"]
      hhd <- length(k_n) #num of links in the route
      
      for (h in 1:hhd){
        k <- (1:N)[node_TSN == k_n[h]]
        a <- (1:N)[node_TSN == a_n[h]]
        LL <- LL + (log(M[k,a]) + (li_theta[h])*V[a] - V[k])
      }
    }
    #print(LL)
    return(LL)
  }
  
  
  return(fr(param)/n_sample)
}