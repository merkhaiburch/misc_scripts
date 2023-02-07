# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-08-19 
# Updated... 2022-08-19

# Description 
# Test out ridge regression on the mtcars dataset
# ---------------------------------------------------------------

# Load in packages
library(glmnet)

# Turn mtcars into df with row names as a column3
mtcars2 <- mtcars
mtcars2$car <- rownames(mtcars2) 
rownames(mtcars2) <- NULL

# Split data into training and testing setes
train = caret::predict(mtcars2, train[,cols])
test = predict(pre_proc_val, test[,cols])

# Subset mtcars
head(mtcars)

# Trun data into matrix
x = as.matrix(mtcars)
y_train = train$mpg

x_test = as.matrix(test_dummies)
y_test = test$unemploy

lambdas <- 10^seq(2, -3, by = -.1)
ridge_reg = glmnet(x, y_train, nlambda = 25, alpha = 0, family = 'gaussian', lambda = lambdas)

summary(ridge_reg)