# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-12-19
# Updated... 2023-12-19
#
# Description:
# Testing exponential decay models in R with the helper package nlstools
# ------------------------------------------------------------------------------

library(nlstools)

# Generate example data
set.seed(123)
distance <- seq(0, 100, by = 10)
expression <- 100 * exp(-0.05 * distance) + rnorm(length(distance), sd = 10)
data <- data.frame(Distance = distance, Expression = expression)

# Plot the data
plot(data$Distance, data$Expression, main = "Non-linear Effect Model",
     xlab = "Distance from Gene", ylab = "Gene Expression")

#  fitting nonlinear regression models requires the provision of starting values for model parameters. 
# A poor choice of starting values may cause non-convergence or convergence to an unwanted local 
# (rather than global) minimum when trying to minimize the least-squares criterion.
#
# nlstools provides the graphical function preview(), which can be used to assess the suitability 
# of the chosen starting values, prior to fitting the model. 
#
# 1) First specify model function
# 2) use this formula as first argument of the function preview(), 
# 3) supply the name of your dataset as second argument
# 4) finally provide a list of values for the model parameters as third argument. 
# An additional argument variable can be used to specify which independent variable is plotted against 
# the dependent variable (column index of the original dataset; default is 1) when more than one independent variable is modeled.
# Expression = a * exp(b * Distance) where a = initial amount and b = growth rate
# For an expoential growth model remove the negative in front of b
formulaExp <- as.formula(Expression ~ a * exp(-b * Distance))
preview(formulaExp, data = data, start = list(a = max(data$Expression), b = 0.05)) 

# Fit a non-linear model given paramets visualized above
model <- nls(Expression ~ a * exp(-b * Distance), data = data, start = list(a = 100, b = 0.05))

# Add the fitted curve to the plot
curve(predict(model, newdata = data.frame(Distance = x)), add = TRUE, col = "red")

# Display the model summary
summary(model)

# Look at overview of model
nlstools::overview(model)

# Look at plot
nlstools::plotfit(model, smooth = F)
nlstools::plotfit(model, smooth = T)

# Assess the goodness of fit through residuals
# Do things fall along the diagnoal of the QQ plot?
model_res <- nlsResiduals(model)
plot(model_res)

# Test residuals
# the null hypothesis of normal distribution could not be rejected (Shapiro-Wilk test: p = if non-significant) 
# and there was also no indication of autocorrelation (runs test: p = if non-significant).
test.nlsResiduals(model_res)


# Test non-linear decay with multiple distances --------------------------------

# Generate example data for multiple distance types
set.seed(123)
distances_A <- seq(0, 100, by = 10)
distances_B <- seq(0, 150, by = 15)
distances_C <- seq(0, 80, by = 8)
expression_data <- matrix(100 * exp(-0.05 * distances) + matrix(rnorm(length(distances) * num_samples, sd = 10), ncol = num_samples), ncol = num_samples)

expression_data <- data.frame(
  Distance_A = rep(distances_A, each = num_samples),
  Distance_B = rep(distances_B, each = num_samples),
  Distance_C = rep(distances_C, each = num_samples),
  Expression = c(expression_data)
)

# Plot the data
plot(expression_data$Distance_A, expression_data$Expression, 
     main = "Exponential Decay Model for Multiple Distances",
     xlab = "Distance A", ylab = "Gene Expression", pch = 16, col = "blue")

# Fit non-linear decay model for multiple distances
multi_distance_model <- nls(Expression ~ a * exp(-b * Distance_A) + c * exp(-d * Distance_B) + e * exp(-f * Distance_C),
                            data = expression_data, start = list(a = 100, b = 0.05, c = 100, d = 0.05, e = 100, f = 0.05))

# Add fitted curve to the plot
curve(predict(multi_distance_model, 
              newdata = data.frame(Distance_A = x, Distance_B = mean(distances_B), Distance_C = mean(distances_C))),
      add = TRUE, col = "red")

# Display the model summary
summary(multi_distance_model)


