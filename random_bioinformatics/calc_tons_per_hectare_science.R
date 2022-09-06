library(dplyr)

# Calculate tons per hectare from science paper

# scale up plots
# 20 rows*10 plants per row*.233 interspace rows *.167 interspace plants=7.78 meters^2/plot
(plot <- 20*10*.233*.167)

# square meters in a hectare
# 1 hectare=(100 meters)^2=10000 meters^2
hectare <- 10000

# number of their plots that could fit in a hectare
# of their plots that could fit in a ha: (10,000m^2)/(7.78 m^2/plot)=1284 plots/ha
(their_plot_hectares <- hectare/plot)

# number of their hectares times their yield per hectare
# 1284 plots/ha*1084g=1.39 T/ha
yield <- 1084 # in grams
(their_plot_hectares*yield)/1000000# 1,392,922



# Confirm yield numbers from second science paper
# 9 rows
# 0.75 m between rows
# 17 plants per row
# 3.5 cm spacing between rows
(plot <- 9*17*0.035*.75)

# square meters in a hectare
# 1 hectare=(100 meters)^2=10000 meters^2
hectare <- 10000

# number of their plots that could fit into a hectare
(their_plot_hectares <- hectare/plot)

# Load in their yield data, look at mean of line by year
sdata <- read.csv("~/Desktop/yield_numbers.csv", header = TRUE)
sdata %>% 
  group_by(Line, year) %>% 
  summarise(mean = mean(num_seeds_div_plant_g))

# Calculate our version of t/ha
sdata$our_t_ha <- (their_plot_hectares * sdata$num_seeds_div_plant_g)/1000000

# Plot their yield vs our yield
library(ggplot2)
a <- ggplot(sdata, aes(x = their_t_ha, y = our_t_ha)) +
  geom_point()+
  geom_smooth(method=lm) +
  facet_grid(~year, scales = "free")
ggsave("~/Desktop/yield_across_lines.png", a)

ggplot(sdata, aes(x = their_t_ha, y = our_t_ha, color = Line)) +
  geom_point()+
  geom_smooth(method=lm) +
  facet_grid(~year, scales = "free")


temp <- sdata %>% filter(year == 2021)
ggplot(temp, aes(x = their_t_ha, y = our_t_ha, color = Line)) +
  geom_point()+
  geom_smooth(method=lm) +
  facet_grid(~year, scales = "free")


