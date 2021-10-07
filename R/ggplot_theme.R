# Select color palette
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Set plot options
axis_text_size <- 8

# To use for fills, add
scale_fill_manual(values=cbPalette)
head(iris)
# ggplot theme
ggplot(iris, aes(x = Sepal.Length, y = Petal.Length, fill = Species)) +
  geom_point() +
  geom_smooth(method = lm) + 
  labs(x = "", 
       y = "", 
       title = "") +
  annotate(geom="text", label = "", x = 5, y = 15) +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size), 
        legend.position="bottom", 
        legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.3, 'cm')) +
  ylim(0, 10) +
  scale_fill_manual(values=cbbPalette)

