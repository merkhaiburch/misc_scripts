# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-03-01 
# Updated... 2022-03-01

# Description 
# Testing case specific random models on test data
# ---------------------------------------------------------------

# Helpful literature
# http://rstudio-pubs-static.s3.amazonaws.com/269829_8285925c922e445097f47925b112841f.html
# https://www.tandfonline.com/doi/full/10.1080/10618600.2014.983641
# https://uc-r.github.io/random_forests

# Load in helpful packages
library(dplyr)
library(ranger)
library(caret)
library(ggplot2)


# ------------------------------------------
# Split into training and testing sets
# ------------------------------------------

index <- caret::createDataPartition(iris$Species, p=0.80, list=FALSE)
testset <- iris[-index,] # testing 20%
trainset <- iris[index,] #training 80%


# ------------------------------------------
# Use Iris dataset for regular random forest
# ------------------------------------------

# Regular random forest
iris_rf <- ranger::ranger(Sepal.Length ~ ., data = trainset, importance = 'impurity')

# Extract relative importance of variables 
importance_iris_rf <- data.frame(iris_rf$variable.importance/max(iris_rf$variable.importance))
importance_iris_rf$features <- rownames(importance_iris_rf)
colnames(importance_iris_rf)[1] <- "relative_importance"

# Plot features in order of importance
importance_iris_rf %>% 
  dplyr::arrange(desc(relative_importance)) %>%
  ggplot(aes(x= reorder(features, relative_importance), y = relative_importance)) +
  geom_col() +
  coord_flip()

# Predict values
iris_rf_values <- predict(iris_rf, testset)$predictions


# ------------------------------------------------
# Use Iris dataset for case-specifc random forest
# ------------------------------------------------

# Case specific random forest
iris_case <- ranger::csrf(Sepal.Length ~ ., 
                          training_data = trainset,
                          test_data = testset,
                          params1 = list(num.trees = 50, mtry = 4), 
                          params2 = list(num.trees = 5))


# ------------------------------------------------
# Compare the distributions of predicted values
# ------------------------------------------------

# Compare the distributions of regular random forest and case specific
model_comparisons <- cbind(iris[,5], data.frame(iris_rf_values), data.frame(iris_case))
colnames(model_comparisons) <- c("species", "rf", "csrf")

# Correlate values
cor(model_comparisons$rf, model_comparisons$csrf)

# Plot predictions against each other
ggplot(aes(x = rf, y = csrf), data = model_comparisons) +
  geom_point() +
  geom_smooth()


