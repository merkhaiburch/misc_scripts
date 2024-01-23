# dopar (doparallel) on a linear model in R
# https://stackoverflow.com/questions/55506788/how-to-transform-a-for-loop-in-a-foreach-loop-in-r
# coefficients: https://stackoverflow.com/questions/36820274/applying-foreach-library-to-a-multi-regression-analysis-in-r
# https://stackoverflow.com/questions/19791609/saving-multiple-outputs-of-foreach-dopar-loop



## Example data 
outcome_list <- list(as.data.frame(cbind(rnorm(32), dataframe_id = c(1))),
                     as.data.frame(cbind(rnorm(32), dataframe_id = c(2))),
                     as.data.frame(cbind(rnorm(32), dataframe_id = c(3))))

## Parallel code
library(doParallel)
registerDoParallel(cl <- makeCluster(3))
results_list <- foreach(i = 1:3) %dopar% {
  
  mylm <- lm(outcome_list[[i]]$V1 ~ mtcars$mpg + mtcars$cyl)
  gof <- broom::glance(mylm)
  betas <- broom::tidy(mylm)
  
  c(outcome_list[[i]]$V2[1], 
    betas$estimate[1],
    betas$estimate[2], betas$p.value[2], 
    betas$estimate[3], betas$p.value[3],
    gof$p.value, gof$r.squared, gof$AIC,
    gof$BIC)
}
stopCluster(cl)

results_df <- setNames(as.data.frame(do.call("rbind", results_list)),
                       c("dataframe_id", "intercept", "b_mpg", "p_mpg", 
                         "b_disp", "p_disp", "p.model", "AIC", "BIC"))


# Old version ----------------------------------------------------------
mtcars <- mtcars #I will use the explanatory variables from here
x <- list()
results_df <- as.data.frame(cbind(dataframe_id = c(0), intercept = c(0),
                                  b_mpg = c(0), p_mpg = c(0),
                                  b_cyl = c(0), p_cyl = c(0),
                                  p.model = c(0), AIC = c(0),
                                  BIC = c(0)))
i <- 1
for(i in 1:3){
  x[[i]] <- lm(outcome_list[[i]]$V1 ~ mtcars$mpg + mtcars$cyl)
  gof <- broom::glance(x[[i]])
  betas <- broom::tidy(x[[i]])
  results_df <- rbind(results_df, c(outcome_list[[i]]$V2[1], 
                                    betas$estimate[1],
                                    betas$estimate[2], betas$p.value[2], 
                                    betas$estimate[3], betas$p.value[3],
                                    gof$p.value, gof$r.squared, gof$AIC,
                                    gof$BIC))
  
  if(i %% i == 0){
    message(paste(i, "of 3")) # To know if my machine has not crashed
    x <- list() # To keep RAM clean of useless data
  }
  gc()
}

results_df <- results_df[-1, ]

