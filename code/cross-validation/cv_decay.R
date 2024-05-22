rm(list=ls(all=TRUE)) #delte all variables
library(parallel)

source("learning_func_decay.R")

seed_value <- 100
n_core <- 10
set.seed(seed_value)  # 結果の再現性のため
kfold <- 10

################
### Settings ###
################
#Estimated intial (wave1) parameter
x <- c(-2.227076383, 0.712642215, -1.292325166, 6.33181656, -0.556804972, -0.51096569, -0.169770418, 1.285324338, 2.496594752, 0.917186604, -1.219138708, -1.543252011, 1.0)

project <- "YOUR FOLDER PATH"
inputpath <- paste(project, "/input", sep = "")

inputlinkpath <- paste(inputpath, "/TSNW/", sep="")
inputroutepath <- paste(inputpath, "/route/user/", sep="")

personal <- read.csv(paste(inputpath, "/PersonalInfo.csv", sep=""), fileEncoding = "Shift-jis")



len_data <- nrow(personal)
colnames <- c("LL_g1", "LL_g2", "LL_g3", "LL_g4", "LL")
results <- data.frame(matrix(NA, nrow = kfold, ncol = 5, dimnames = list( c(), colnames )))

cross_validation <- function(ifold) {
  # データをランダムにシャッフル
  shuffled_indices <- sample(len_data)
  # トレーニングセットとテストセットに分割
  split_index <- ceiling(0.8 * len_data)
  train_indices <- shuffled_indices[1:split_index]
  test_indices <- shuffled_indices[(split_index + 1):len_data]
  
  ### training data reading 1 ###
  #route data
  route_file2 <- sort(list.files(inputroutepath, pattern="wave2"))[intersect(train_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  route_file3 <- sort(list.files(inputroutepath, pattern="wave3"))[intersect(train_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  route_file4 <- sort(list.files(inputroutepath, pattern="wave4"))[intersect(train_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  #link data
  link_file2 <- sort(list.files(inputlinkpath, pattern="wave2"))[intersect(train_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  link_file3 <- sort(list.files(inputlinkpath, pattern="wave3"))[intersect(train_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  link_file4 <- sort(list.files(inputlinkpath, pattern="wave4"))[intersect(train_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  list_route <- list(route_file2, route_file3, route_file4)
  list_link <- list(link_file2, link_file3, link_file4)
  learning_param1 <- learning_model_estimation(list_route, list_link, x)
  
  ### training data reading 2 ###
  #route data
  route_file2 <- sort(list.files(inputroutepath, pattern="wave2"))[intersect(train_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  route_file3 <- sort(list.files(inputroutepath, pattern="wave3"))[intersect(train_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  route_file4 <- sort(list.files(inputroutepath, pattern="wave4"))[intersect(train_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  #link data
  link_file2 <- sort(list.files(inputlinkpath, pattern="wave2"))[intersect(train_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  link_file3 <- sort(list.files(inputlinkpath, pattern="wave3"))[intersect(train_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  link_file4 <- sort(list.files(inputlinkpath, pattern="wave4"))[intersect(train_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  list_route <- list(route_file2, route_file3, route_file4)
  list_link <- list(link_file2, link_file3, link_file4)
  learning_param2 <- learning_model_estimation(list_route, list_link, x)
  
  ### training data reading 3 ###
  #route data
  route_file2 <- sort(list.files(inputroutepath, pattern="wave2"))[intersect(train_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  route_file3 <- sort(list.files(inputroutepath, pattern="wave3"))[intersect(train_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  route_file4 <- sort(list.files(inputroutepath, pattern="wave4"))[intersect(train_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  #link data
  link_file2 <- sort(list.files(inputlinkpath, pattern="wave2"))[intersect(train_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  link_file3 <- sort(list.files(inputlinkpath, pattern="wave3"))[intersect(train_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  link_file4 <- sort(list.files(inputlinkpath, pattern="wave4"))[intersect(train_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  list_route <- list(route_file2, route_file3, route_file4)
  list_link <- list(link_file2, link_file3, link_file4)
  learning_param3 <- learning_model_estimation(list_route, list_link, x)
  
  ### training data reading 4 ###
  #route data
  route_file2 <- sort(list.files(inputroutepath, pattern="wave2"))[intersect(train_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  route_file3 <- sort(list.files(inputroutepath, pattern="wave3"))[intersect(train_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  route_file4 <- sort(list.files(inputroutepath, pattern="wave4"))[intersect(train_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  #link data
  link_file2 <- sort(list.files(inputlinkpath, pattern="wave2"))[intersect(train_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  link_file3 <- sort(list.files(inputlinkpath, pattern="wave3"))[intersect(train_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  link_file4 <- sort(list.files(inputlinkpath, pattern="wave4"))[intersect(train_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  list_route <- list(route_file2, route_file3, route_file4)
  list_link <- list(link_file2, link_file3, link_file4)
  learning_param4 <- learning_model_estimation(list_route, list_link, x)
  
  
  
  
  # Validation. Log-likelihoodの計算
  ### test data reading 1 ###
  #route data
  route_file2 <- sort(list.files(inputroutepath, pattern="wave2"))[intersect(test_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  route_file3 <- sort(list.files(inputroutepath, pattern="wave3"))[intersect(test_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  route_file4 <- sort(list.files(inputroutepath, pattern="wave4"))[intersect(test_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  #link data
  link_file2 <- sort(list.files(inputlinkpath, pattern="wave2"))[intersect(test_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  link_file3 <- sort(list.files(inputlinkpath, pattern="wave3"))[intersect(test_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  link_file4 <- sort(list.files(inputlinkpath, pattern="wave4"))[intersect(test_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  list_route <- list(route_file2, route_file3, route_file4)
  list_link <- list(link_file2, link_file3, link_file4)
  n_g1 <- length(list_route[[1]])
  
  predictions1 <- learning_model_prediction(list_route, list_link, x, learning_param1)
  #predictions1 <- 1
  
  ### test data reading 2 ###
  #route data
  route_file2 <- sort(list.files(inputroutepath, pattern="wave2"))[intersect(test_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  route_file3 <- sort(list.files(inputroutepath, pattern="wave3"))[intersect(test_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  route_file4 <- sort(list.files(inputroutepath, pattern="wave4"))[intersect(test_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  #link data
  link_file2 <- sort(list.files(inputlinkpath, pattern="wave2"))[intersect(test_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  link_file3 <- sort(list.files(inputlinkpath, pattern="wave3"))[intersect(test_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  link_file4 <- sort(list.files(inputlinkpath, pattern="wave4"))[intersect(test_indices, which(is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  list_route <- list(route_file2, route_file3, route_file4)
  list_link <- list(link_file2, link_file3, link_file4)
  n_g2 <- length(list_route[[1]])
  
  predictions2 <- learning_model_prediction(list_route, list_link, x, learning_param2)
  #predictions2 <- 2
  
  ### test data reading 3 ###
  #route data
  route_file2 <- sort(list.files(inputroutepath, pattern="wave2"))[intersect(test_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  route_file3 <- sort(list.files(inputroutepath, pattern="wave3"))[intersect(test_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  route_file4 <- sort(list.files(inputroutepath, pattern="wave4"))[intersect(test_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  #link data
  link_file2 <- sort(list.files(inputlinkpath, pattern="wave2"))[intersect(test_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  link_file3 <- sort(list.files(inputlinkpath, pattern="wave3"))[intersect(test_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  link_file4 <- sort(list.files(inputlinkpath, pattern="wave4"))[intersect(test_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==1)))]
  list_route <- list(route_file2, route_file3, route_file4)
  list_link <- list(link_file2, link_file3, link_file4)
  n_g3 <- length(list_route[[1]])
  
  predictions3 <- learning_model_prediction(list_route, list_link, x, learning_param3)
  #predictions3 <- 3
  
  ### test data reading 4 ###
  #route data
  route_file2 <- sort(list.files(inputroutepath, pattern="wave2"))[intersect(test_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  route_file3 <- sort(list.files(inputroutepath, pattern="wave3"))[intersect(test_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  route_file4 <- sort(list.files(inputroutepath, pattern="wave4"))[intersect(test_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  #link data
  link_file2 <- sort(list.files(inputlinkpath, pattern="wave2"))[intersect(test_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  link_file3 <- sort(list.files(inputlinkpath, pattern="wave3"))[intersect(test_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  link_file4 <- sort(list.files(inputlinkpath, pattern="wave4"))[intersect(test_indices, which(!is.na(personal$Depth_nodeファイルの値使う)&(personal$sp1避難有無==0)))]
  list_route <- list(route_file2, route_file3, route_file4)
  list_link <- list(link_file2, link_file3, link_file4)
  n_g4 <- length(list_route[[1]])
  
  predictions4 <- learning_model_prediction(list_route, list_link, x, learning_param4)
  #predictions4 <- 4
  
  # 結果の評価
  print(paste("ifold:", ifold))
  print("learning_param1:")
  print(learning_param1)
  print("learning_param2:")
  print(learning_param2)
  print("learning_param3:")
  print(learning_param3)
  print("learning_param4:")
  print(learning_param4)
  print("predictions1:")
  print(predictions1)
  print("predictions2:")
  print(predictions2)
  print("predictions3:")
  print(predictions3)
  print("predictions4:")
  print(predictions4)
  return(list(predictions1, 
              predictions2, 
              predictions3, 
              predictions4, 
              (n_g1*predictions1 + n_g2*predictions2 + n_g3*predictions3 + n_g4*predictions4)/(n_g1+n_g2+n_g3+n_g4) #全体の平均
              )
         )
}


parallels <- mclapply(1:kfold, function(ifold) { #並列化処理
  cross_validation(ifold)
}, mc.cores = n_core)

for (ifold in 1:kfold) {
  results[ifold, "LL_g1"] <- parallels[[ifold]][1]
  results[ifold, "LL_g2"] <- parallels[[ifold]][2]
  results[ifold, "LL_g3"] <- parallels[[ifold]][3]
  results[ifold, "LL_g4"] <- parallels[[ifold]][4]
  results[ifold, "LL"] <- parallels[[ifold]][5]
}

#resultをcsvに出力
write.csv(results, paste(project, "/result/","learning_cv",seed_value,".csv", sep = ""))
mean_LL <- mean(results$LL)
print(paste("Average LL:", mean_LL))
