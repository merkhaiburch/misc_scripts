# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-08-11
#
# Description 
#   - Linear Regression for Prediction
# ---------------------------------------------------------------

# Create a data set using the following code:
set.seed(1, sample.kind="Rounding") 
n <- 100
Sigma <- 9*matrix(c(1.0, 0.5, 0.5, 1.0), 2, 2)
dat <- MASS::mvrnorm(n = 100, c(69, 69), Sigma) %>%
  data.frame() %>% setNames(c("x", "y"))

# Build 100 linear models using the data from above and calculate
# - mean & standard deviation.
library(caret)
set.seed(1, sample.kind="Rounding")
run_model <- function(dat) {
  # Separate into test and training
  test_index <- caret::createDataPartition(dat$y, times = 1, p = 0.5, list = FALSE)
  test_set <- dat[test_index, ]
  train_set <- dat[-test_index, ]
  
  # Run model
  model <- lm(y~x, data = train_set)
  model_pred <- predict(model, test_set)
  caret::RMSE(model_pred, test_set$y)
}

temp <- replicate(100, run_model(dat))

# Calculate mean and sd
mean(temp)
sd(temp)


# Now we will repeat the exercise above but using larger datasets. 
# Write a function that takes a size n, then 
# (1) builds a dataset using the code provided at the top of Q1 but 
#     with n observations instead of 100 and without the set.seed(1), 
# (2) runs the replicate() loop that you wrote to answer Q1, which builds 
#   100 linear models and returns a vector of RMSEs, and 
# (3) calculates the mean and standard deviation of the 100 RMSEs.
q2 <- function(n_size) {
  # Make dummy dataset
  n <- n_size
  Sigma <- 9*matrix(c(1.0, 0.5, 0.5, 1.0), 2, 2)
  dat_large <- MASS::mvrnorm(n = n_size, c(69, 69), Sigma) %>%
    data.frame() %>% setNames(c("x", "y"))
  
  # Run replicate loop
  temp <- replicate(100, run_model(dat_large))
  
  return(c(nrow(dat_large), mean(temp), sd(temp)))
}
 # Set seed
set.seed(1, sample.kind="Rounding")

# Use sapply to apply function to get results
n <- c(100, 500, 1000, 5000, 10000)
sapply(n, q2)


# Q4
# Now repeat the exercise from Q1, this time making the correlation 
# between x and y larger, as in the following code:
set.seed(1, sample.kind="Rounding")
n <- 100
Sigma <- 9*matrix(c(1.0, 0.95, 0.95, 1.0), 2, 2)
dat_2 <- MASS::mvrnorm(n = 100, c(69, 69), Sigma) %>%
  data.frame() %>% setNames(c("x", "y"))
set.seed(1, sample.kind="Rounding")
temp <- replicate(100, run_model(dat_2))

# Calculate mean and sd
mean(temp)
sd(temp)

# Q5
set.seed(1, sample.kind="Rounding")
Sigma <- matrix(c(1.0, 0.75, 0.75, 0.75, 1.0, 0.25, 0.75, 0.25, 1.0), 3, 3)
dat_3 <- MASS::mvrnorm(n = 100, c(0, 0, 0), Sigma) %>%
  data.frame() %>% setNames(c("y", "x_1", "x_2"))
cor(dat_3)

set.seed(1, sample.kind="Rounding")

test_index <- caret::createDataPartition(dat_3$y, times = 1, p = 0.5, list = FALSE)
test_set <- dat_3[test_index, ]
train_set <- dat_3[-test_index, ]

model <- lm(y~x_1, data = train_set)
model_pred <- predict(model, test_set)
caret::RMSE(model_pred, test_set$y)

model <- lm(y~x_2, data = train_set)
model_pred <- predict(model, test_set)
caret::RMSE(model_pred, test_set$y)

model <- lm(y~x_1+x_2, data = train_set)
model_pred <- predict(model, test_set)
caret::RMSE(model_pred, test_set$y)

# Logistic regression
set.seed(2, sample.kind="Rounding") #if you are using R 3.6 or later
make_data <- function(n = 1000, p = 0.5, 
                      mu_0 = 0, mu_1 = 2, 
                      sigma_0 = 1,  sigma_1 = 1){
  
  y <- rbinom(n, 1, p)
  f_0 <- rnorm(n, mu_0, sigma_0)
  f_1 <- rnorm(n, mu_1, sigma_1)
  x <- ifelse(y == 1, f_1, f_0)
  
  test_index <- createDataPartition(y, times = 1, p = 0.5, list = FALSE)
  
  list(train = data.frame(x = x, y = as.factor(y)) %>% slice(-test_index),
       test = data.frame(x = x, y = as.factor(y)) %>% slice(test_index))
}
dat <- make_data()
dat$train %>% ggplot(aes(x, color = y)) + geom_density()

# Set seed to 1
set.seed(1, sample.kind="Rounding")
# Use make_data to generate 25 different datasets
mu_1_25 <- seq(0, 3, len=25)
q1 <- function(fds) {
  temp <- make_data(mu_1 = fds)
  glm_fit <- temp$train %>% 
    glm(y ~ x, data=., family = "binomial")
  p_hat_logit <- predict(glm_fit, newdata = temp$test, type = "response")
  y_hat_logit <- ifelse(p_hat_logit > 0.5, "1", "0") %>% factor
  confusionMatrix(y_hat_logit, temp$test$y)$overall[["Accuracy"]]
}
lala <- sapply(mu_1_25, q1)
plot(lala~mu_1_25)

