# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-08-20 
# Updated... 2022-08-20

# Description 
# The data taken from Supp. tables D1, D9, and D10
# Science, 2021
# ---------------------------------------------------------------

# Load packages
library(dplyr)
library(ggplot2)

# Color blind pallete
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Load in their yield data
sdata <- read.csv("/Users/merri/Downloads/yield_numbers.csv", header = TRUE)

# Calculate median by line and year
sdata |> group_by(Line, year) |> summarise(median = median(their_t_ha))

# Change the order of the WT line
level_order <- factor(sdata$Line, level = c("WT","YZ-19-8A","ND-18-56",
                                            "ND-18-34A", "ND-17-20",
                                            "YZ-26-1C", "YZ-19-21",
                                            "ND-18-44","ND-18-12B"))

# Plot yield boxplots both years
yield_plot <- ggplot(sdata, aes(x = level_order, y = their_t_ha, fill = as.factor(year))) +
  geom_boxplot() +
  theme_bw() +
  geom_hline(yintercept = 4.70, color = "black", size = .6) +
  xlab("Line") +
  ylab("Yield (t/ha)") +
  scale_fill_manual(values=c("#56B4E9","#CC79A7")) +
  ggtitle("Science soy yield data from Supp. D1, D9, and D10")

# Save to file
ggsave(yield_plot, "/Users/merri/Downloads/science_soy_yield_2020_2021.png",
       width = 5, height = 6, units = "in")


# Run lm -----------------------------------------------------------------------

# Across year model for lines represented in both years
yield_complete <- sdata |> 
  filter(Line %in% c("WT", "ND-18-34A", "YZ-19-21","YZ-19-21"))
summary(lm(their_t_ha ~ Line + year, data = yield_complete))


# 2020 model
yield_2020 <- sdata |> filter(year == 2020)
summary(lm(their_t_ha~Line, data = yield_2020))

# 2021 model
yield_2021 <- sdata |> filter(year == 2021)
summary(lm(their_t_ha~Line, data = yield_2021))







