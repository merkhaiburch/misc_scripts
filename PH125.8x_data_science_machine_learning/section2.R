# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-08-11
#
# Description 
#   - Conditional expectations and loss function
# ---------------------------------------------------------------

# We are now going to write code to compute conditional probabilities for being male 
# in the heights dataset. Round the heights to the closest inch. Plot the estimated 
# conditional probability  𝑃(𝑥)=Pr(Male|height=𝑥)  for each  𝑥 .
library(dslabs)
data("heights")
heights %>% 
  mutate(height = round(height)) %>%
  group_by(height) %>%
  summarize(p = mean(sex == "Male")) %>%
qplot(height, p, data =.)

#Q7
# Cut height into quantiles
ps <- seq(0, 1, 0.1)
heights %>% 
  mutate(g = cut(male, quantile(height, ps), include.lowest = TRUE)) %>%
  group_by(g) %>%
  summarize(p = mean(sex == "Male"), height = mean(height)) %>%
  qplot(height, p, data =.)

# Q8
# You can generate data from a bivariate normal distrubution using the 
# MASS package using the following code:
Sigma <- 9*matrix(c(1,0.5,0.5,1), 2, 2)
dat <- MASS::mvrnorm(n = 10000, c(69, 69), Sigma) %>%
  data.frame() %>% setNames(c("x", "y"))
plot(dat)

# Using an approach similar to that used in the previous exercise, 
# let's estimate the conditional expectations and make a plot. Part 
# of the code has again been provided for you:
ps <- seq(0, 1, 0.1)
dat %>% 
  mutate(g = cut(x, quantile(x, ps), include.lowest = TRUE)) %>%
  group_by(g) %>%
  summarize(y = mean(y), x = mean(x)) %>%
	qplot(x, y, data =.)

