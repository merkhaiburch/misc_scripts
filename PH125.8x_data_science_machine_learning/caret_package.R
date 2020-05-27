# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-04-28
#
# Description 
#   - Caret package, training and test sets, and overall accuracy
# ---------------------------------------------------------------

# --------------------------------------------------------------
# Caret package, training and test sets, and overall accuracy
# --------------------------------------------------------------

# Load packages
library(tidyverse)
library(caret)
library(dslabs)
data(heights)

# define the outcome and predictors
y <- heights$sex
x <- heights$height

# generate training and test sets
set.seed(2007)
test_index <- createDataPartition(y, times = 1, p = 0.5, list = FALSE)
test_set <- heights[test_index, ]
train_set <- heights[-test_index, ]

# guess the outcome
y_hat <- sample(c("Male", "Female"), length(test_index), replace = TRUE)
y_hat <- sample(c("Male", "Female"), length(test_index), replace = TRUE) %>% 
  factor(levels = levels(test_set$sex))

# compute accuracy
mean(y_hat == test_set$sex)
heights %>% group_by(sex) %>% summarize(mean(height), sd(height))
y_hat <- ifelse(x > 62, "Male", "Female") %>% factor(levels = levels(test_set$sex))
mean(y == y_hat)

# examine the accuracy of 10 cutoffs
cutoff <- seq(61, 70)
accuracy <- map_dbl(cutoff, function(x){
  y_hat <- ifelse(train_set$height > x, "Male", "Female") %>% 
    factor(levels = levels(test_set$sex))
  mean(y_hat == train_set$sex)
})
max(accuracy)
plot(cutoff, accuracy)
best_cutoff <- cutoff[which.max(accuracy)]
best_cutoff
y_hat <- ifelse(test_set$height > best_cutoff, "Male", "Female") %>% 
  factor(levels = levels(test_set$sex))
y_hat <- factor(y_hat)
mean(y_hat == test_set$sex)


# --------------------------------------------------------------
# Basics of Evaluating Machine Learning Algorithms
# --------------------------------------------------------------

# tabulate each combination of prediction and actual value
table(predicted = y_hat, actual = test_set$sex)
test_set %>% 
  mutate(y_hat = y_hat) %>%
  group_by(sex) %>% 
  summarize(accuracy = mean(y_hat == sex))
prev <- mean(y == "Male")

confusionMatrix(data = y_hat, reference = test_set$sex)


# --------------------------------------------------------------
# Balanced Accuracy and F1 Score
# --------------------------------------------------------------

# maximize F-score
cutoff <- seq(61, 70)
F_1 <- map_dbl(cutoff, function(x){
  y_hat <- ifelse(train_set$height > x, "Male", "Female") %>% 
    factor(levels = levels(test_set$sex))
  F_meas(data = y_hat, reference = factor(train_set$sex))
})
max(F_1)

best_cutoff <- cutoff[which.max(F_1)]
y_hat <- ifelse(test_set$height > best_cutoff, "Male", "Female") %>% 
  factor(levels = levels(test_set$sex))
sensitivity(data = y_hat, reference = test_set$sex)
specificity(data = y_hat, reference = test_set$sex)


# --------------------------------------------------------------
# Roc and precision recall curves
# --------------------------------------------------------------

p <- 0.9
n <- length(test_index)
y_hat <- sample(c("Male", "Female"), n, replace = TRUE, prob=c(p, 1-p)) %>% 
  factor(levels = levels(test_set$sex))
mean(y_hat == test_set$sex)

# ROC curve
probs <- seq(0, 1, length.out = 10)
guessing <- map_df(probs, function(p){
  y_hat <- 
    sample(c("Male", "Female"), n, replace = TRUE, prob=c(p, 1-p)) %>% 
    factor(levels = c("Female", "Male"))
  list(method = "Guessing",
       FPR = 1 - specificity(y_hat, test_set$sex),
       TPR = sensitivity(y_hat, test_set$sex))
})
guessing %>% qplot(FPR, TPR, data =., xlab = "1 - Specificity", ylab = "Sensitivity")

cutoffs <- c(50, seq(60, 75), 80)
height_cutoff <- map_df(cutoffs, function(x){
  y_hat <- ifelse(test_set$height > x, "Male", "Female") %>% 
    factor(levels = c("Female", "Male"))
  list(method = "Height cutoff",
       FPR = 1-specificity(y_hat, test_set$sex),
       TPR = sensitivity(y_hat, test_set$sex))
})

# plot both curves together
bind_rows(guessing, height_cutoff) %>%
  ggplot(aes(FPR, TPR, color = method)) +
  geom_line() +
  geom_point() +
  xlab("1 - Specificity") +
  ylab("Sensitivity")

library(ggrepel)
map_df(cutoffs, function(x){
  y_hat <- ifelse(test_set$height > x, "Male", "Female") %>% 
    factor(levels = c("Female", "Male"))
  list(method = "Height cutoff",
       cutoff = x, 
       FPR = 1-specificity(y_hat, test_set$sex),
       TPR = sensitivity(y_hat, test_set$sex))
}) %>%
  ggplot(aes(FPR, TPR, label = cutoff)) +
  geom_line() +
  geom_point() +
  geom_text_repel(nudge_x = 0.01, nudge_y = -0.01)

# plot precision against recall
guessing <- map_df(probs, function(p){
  y_hat <- sample(c("Male", "Female"), length(test_index), 
                  replace = TRUE, prob=c(p, 1-p)) %>% 
    factor(levels = c("Female", "Male"))
  list(method = "Guess",
       recall = sensitivity(y_hat, test_set$sex),
       precision = precision(y_hat, test_set$sex))
})

height_cutoff <- map_df(cutoffs, function(x){
  y_hat <- ifelse(test_set$height > x, "Male", "Female") %>% 
    factor(levels = c("Female", "Male"))
  list(method = "Height cutoff",
       recall = sensitivity(y_hat, test_set$sex),
       precision = precision(y_hat, test_set$sex))
})

bind_rows(guessing, height_cutoff) %>%
  ggplot(aes(recall, precision, color = method)) +
  geom_line() +
  geom_point()
guessing <- map_df(probs, function(p){
  y_hat <- sample(c("Male", "Female"), length(test_index), replace = TRUE, 
                  prob=c(p, 1-p)) %>% 
    factor(levels = c("Male", "Female"))
  list(method = "Guess",
       recall = sensitivity(y_hat, relevel(test_set$sex, "Male", "Female")),
       precision = precision(y_hat, relevel(test_set$sex, "Male", "Female")))
})

height_cutoff <- map_df(cutoffs, function(x){
  y_hat <- ifelse(test_set$height > x, "Male", "Female") %>% 
    factor(levels = c("Male", "Female"))
  list(method = "Height cutoff",
       recall = sensitivity(y_hat, relevel(test_set$sex, "Male", "Female")),
       precision = precision(y_hat, relevel(test_set$sex, "Male", "Female")))
})
bind_rows(guessing, height_cutoff) %>%
  ggplot(aes(recall, precision, color = method)) +
  geom_line() +
  geom_point()


# Comprehension Check: Practice with Machine Learning, Part 1
library(dslabs)
library(dplyr)
library(lubridate)
data(reported_heights)

dat <- mutate(reported_heights, date_time = ymd_hms(time_stamp)) %>%
  filter(date_time >= make_date(2016, 01, 25) & date_time < make_date(2016, 02, 1)) %>%
  mutate(type = ifelse(day(date_time) == 25 & hour(date_time) == 8 & between(minute(date_time), 15, 30), "inclass","online")) %>%
  select(sex, type)

dat$sex <- factor(dat$sex, c("Female", "Male"))
x <- dat$type

table(dat)

# Use type (online/in person) to predict sex
y_hat <- ifelse(dat$type == "online", "Male", "Female") %>% factor(levels = levels(as.factor(dat$sex)))
mean(y == y_hat)

# Write a line of code using the table() function to show the confusion matrix between y_hat and y
# table(y_hat, y)
table(predicted = y_hat, actual = dat$sex)

# Test the sensitivity of the prediction
caret::sensitivity(data = y_hat, reference = dat$sex)

# Test the specificity of tge prection
caret::specificity(data = y_hat, reference = dat$sex)

# What is the prevalance of females in the dataset (in the entire confusion matrix)
confusionMatrix(data = y_hat, reference = dat$sex)


# ------------------------------
# Comprehension Check: Practice with Machine Learning, Part 2
# ------------------------------

# Use the iris dataset, remove the setosa species and foucs on
# versicolor and virginica iris species

library(caret)
data(iris)
iris <- iris[-which(iris$Species=='setosa'),]
y <- iris$Species

# create an even split of the data into train and test partitions using 
#    createDataPartition() from the caret package.
set.seed(2, sample.kind="Rounding")
test_index <- createDataPartition(y, times = 1, p = 0.5, list = FALSE)
test <- iris[test_index,]
train <- iris[-test_index,]

# Next we will figure out the singular feature in the dataset that yields the 
#   greatest overall accuracy when predicting species. You can use the code 
#   from the introduction and from Q7 to start your analysis.
# Using only the train iris dataset, for each feature, perform a simple search to 
#   find the cutoff that produces the highest accuracy, predicting virginica if greater 
#   than the cutoff and versicolor otherwise. Use the seq function over the range of each 
#   feature by intervals of 0.1 for this search.

# # guess the outcome
# y_hat <- sample(c("versicolor", "virginica"), length(test), replace = TRUE)
# y_hat <- sample(c("versicolor","virginica"), length(test), replace = TRUE) %>% 
#   factor(levels = levels(test$Species))
# 
# # compute accuracy
# mean(y_hat == test$Species)

# Set up variables
x_dat <- train$Petal.Length

# examine the accuracy of multiple cutoffs (ONLY PRODUCES THE RIGHT ANSWER WHEN SPECIES LEVELS ARE IN THIS ORDER)
cutoff <- seq(min(x_dat), max(x_dat), by = 0.1)
accuracy <- sapply(cutoff, function(x){
  y_hat <- ifelse(train$Petal.Length > x, "virginica", "versicolor") 
  mean(y_hat == train$Species)
})
max(accuracy)
plot(cutoff, accuracy)
best_cutoff <- cutoff[which.max(accuracy)]
best_cutoff

# From the answers
foo <- function(x){
  rangedValues <- seq(range(x)[1],range(x)[2],by=0.1)
  print(rangedValues)
  sapply(rangedValues,function(i){
    y_hat <- ifelse(x>i,'virginica','versicolor')
    mean(y_hat==train$Species)
  })
}
predictions <- apply(test[,-5],2,foo)
sapply(predictions,max)	

# Use the best cutoff value from the training data to calculate accuracy in the test data
y_hat <- ifelse(test$Petal.Length > 4, "versicolor","virginica")
mean(y_hat == train$Species)


# CORRECT ANSWER
var1 <- 3
predictions <- foo(train[,var1])
rangedValues <- seq(range(train[,var1])[1],range(train[,var1])[2],by=0.1)
cutoffs <-rangedValues[which(predictions==max(predictions))]

y_hat <- ifelse(test[,var1]>cutoffs[1],'virginica','versicolor')
mean(y_hat==test$Species)

# Use the test data as the training data 
# From the answers
foo <- function(x){
  rangedValues <- seq(range(x)[1],range(x)[2],by=0.1)
  print(rangedValues)
  sapply(rangedValues,function(i){
    y_hat <- ifelse(x>i,'virginica','versicolor')
    mean(y_hat==test$Species)
  })
}
predictions <- apply(test[,-5],2,foo)
sapply(predictions,max)	

# Use the best cutoff value from the training data to calculate accuracy in the test data
y_hat <- ifelse(test$Petal.Length > 4, "versicolor","virginica")
mean(y_hat == test$Species)


# Plot the data
plot(iris,pch=21,bg=iris$Species)

# OPtimize test data to be the new training data
foo <- function(x){
  rangedValues <- seq(range(x)[1],range(x)[2],by=0.1)
  sapply(rangedValues,function(i){
    y_hat <- ifelse(x>i,'virginica','versicolor')
    mean(y_hat==test$Species)
  })
}
predictions <- apply(test[,-5],2,foo)
sapply(predictions,max)	

var1 <- 4
predictions <- foo(test[,var1])
rangedValues <- seq(range(test[,var1])[1],range(test[,var1])[2],by=0.1)
cutoffs <-rangedValues[which(predictions==max(predictions))]

y_hat <- ifelse(train[,var1]>cutoffs[1],'virginica','versicolor')
mean(y_hat==train$Species)