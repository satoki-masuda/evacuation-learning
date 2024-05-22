rm(list=ls(all=TRUE)) #delte all variables
library(parallel)

source("learning_func_no_learning.R")

seed_value <- 100
n_core <- 10
set.seed(seed_value)  # 結果の再現性のため
kfold <- 10

################
### Settings ###
################
project <- "YOUR FOLDER PATH"
inputpath <- paste(project, "/input", sep = "")

inputlinkpath <- paste(inputpath, "/TSNW/", sep="")
inputroutepath <- paste(inputpath, "/route/user/", sep="")

personal <- read.csv(paste(inputpath, "/PersonalInfo.csv", sep=""), fileEncoding = "Shift-jis")

len_data <- nrow(personal)
colnames <- c("LL")
results <- data.frame(matrix(NA, nrow = kfold, ncol = 1, dimnames = list( c(), colnames )))

cross_validation <- function(ifold) {
  # データをランダムにシャッフル
  shuffled_indices <- sample(len_data)
  # トレーニングセットとテストセットに分割
  split_index <- ceiling(0.8 * len_data)
  train_indices <- shuffled_indices[1:split_index]
  test_indices <- shuffled_indices[(split_index + 1):len_data]
  
  #route data
  route_file1 <- sort(list.files(inputroutepath, pattern="wave1"))[train_indices]
  route_file2 <- sort(list.files(inputroutepath, pattern="wave2"))[train_indices]
  route_file3 <- sort(list.files(inputroutepath, pattern="wave3"))[train_indices]
  route_file4 <- sort(list.files(inputroutepath, pattern="wave4"))[train_indices]
  route_file <- c(route_file1, route_file2, route_file3, route_file4)
  #link data
  link_file1 <- sort(list.files(inputlinkpath, pattern="wave1"))[train_indices]
  link_file2 <- sort(list.files(inputlinkpath, pattern="wave2"))[train_indices]
  link_file3 <- sort(list.files(inputlinkpath, pattern="wave3"))[train_indices]
  link_file4 <- sort(list.files(inputlinkpath, pattern="wave4"))[train_indices]
  link_file <- c(link_file1, link_file2, link_file3, link_file4)
  
  learning_param <- learning_model_estimation(route_file, link_file)
  
  
  # Validation. Log-likelihoodの計算
  #route data
  #route_file1 <- sort(list.files(inputroutepath, pattern="wave1"))[test_indices]
  route_file2 <- sort(list.files(inputroutepath, pattern="wave2"))[test_indices]
  route_file3 <- sort(list.files(inputroutepath, pattern="wave3"))[test_indices]
  route_file4 <- sort(list.files(inputroutepath, pattern="wave4"))[test_indices]
  route_file <- c(route_file2, route_file3, route_file4)
  #link data
  #link_file1 <- sort(list.files(inputlinkpath, pattern="wave1"))[test_indices]
  link_file2 <- sort(list.files(inputlinkpath, pattern="wave2"))[test_indices]
  link_file3 <- sort(list.files(inputlinkpath, pattern="wave3"))[test_indices]
  link_file4 <- sort(list.files(inputlinkpath, pattern="wave4"))[test_indices]
  link_file <- c(link_file2, link_file3, link_file4)
  
  predictions <- learning_model_prediction(route_file, link_file, learning_param)
  
  # 結果の評価
  print(paste("ifold:", ifold))
  print("learning_param:")
  print(learning_param)
  print("predictions:")
  print(predictions)
  return(predictions)
}


parallels <- mclapply(1:kfold, function(ifold) { #並列化処理
  cross_validation(ifold)
}, mc.cores = n_core)

for (ifold in 1:kfold) {
  results[ifold, "LL"] <- parallels[[ifold]][1]
}

#resultをcsvに出力
write.csv(results, paste(project, "/result/","learning_cv_no_learning",seed_value,".csv", sep = ""))
mean_LL <- mean(results$LL)
print(paste("Average LL:", mean_LL))
