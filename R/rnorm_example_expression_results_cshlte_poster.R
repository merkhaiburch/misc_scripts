# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-10-04
# Updated... 2022-10-04
#
# Description:
# Simulate data and plot for test poster 
# ------------------------------------------------------------------------------

# Load packages
library(dplyr)
library(ggplot2)

# Simulate data
flag <- rnorm(400, mean = 0.2, sd = 0.1) %>% data.frame()
tip <- rnorm(400, mean = 0.4, sd = 0.1) %>% data.frame()
base <- rnorm(400, mean = 0.6, sd = 0.1) %>% data.frame()
grow <- rnorm(400, mean = 0.8, sd = 0.1) %>% data.frame()
expression <- rbind(flag, tip, base, grow)

# Add identifier
expression$Tissue <- c(rep("Flag Leaf", 400),
                       rep("Leaf Tip", 400),
                       rep("Leaf Base", 400),
                       rep("Growing Point", 400))
colnames(expression)[1] <- "Value"

# Plot results 
crPal <- croix::croix_palette(name = "mov_edward_scissorhands", n = 5)
axis_text_size <- 22
plot1 <- ggplot(expression, aes(x = Value, fill = Tissue)) +
  geom_histogram(position="identity", alpha=0.5, bins = 50) +
  scale_colour_manual(values=crPal) +
  scale_fill_manual(values=crPal) +
  labs(x = "Mixed Model P-value",
       y = "Count") +
  theme_bw() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size))  +
  theme(legend.position="bottom",
        legend.key.size = unit(1.5, 'cm'),
        legend.text = element_text(size=axis_text_size),
        legend.title = element_text(size=axis_text_size))
plot1

# save
ggsave("~/Box Sync/Cornell_PhD/presentations/2022_CSHL_TE/example_expression.png", 
       plot1, height = 10, width = 11.5, units = "in")
