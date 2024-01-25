library(doParallel)

registerDoParallel(cores = 40)

results_list <- foreach(i = 1:3000) %dopar% {
  
  mylm <- lm(response_df[[i]]$V1 ~ predictor1 + predictor2)
  
  result <- broom::glance(mylm)
  
  c(result)
}


outsideLoopHolder <- c()

for(i in 1:3000){
  
  mylm <- lm(response_df[[i]]$V1 ~ predictor1 + predictor2)
  
  result <- broom::glance(mylm)
  
  outsideLoopHolder <- rbind(outsideLoopHolder, result)
}


