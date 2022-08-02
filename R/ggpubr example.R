library(ggpubr)
library(ggplot2)
iris <- iris


a <- ggplot(data=iris, aes(x = Sepal.Length, y = Sepal.Width)) + 
  geom_point(aes(color=Species, shape=Species)) +
  xlab("Sepal Length") +  ylab("Sepal Width") +
  ggtitle("Sepal Length-Width")

b <-  ggplot(data=iris, aes(x=Species, y=Sepal.Length)) + 
  geom_boxplot(aes(fill=Species)) + 
  ylab("Sepal Length") + ggtitle("Iris Boxplot") +
  stat_summary(fun=mean, geom="point", shape=5, size=4) 

c <- ggplot(data=iris, aes(x=Sepal.Length, y=Sepal.Width, color=Species)) + 
  geom_point(aes(shape=Species), size=1.5) + xlab("Sepal Length") + 
  ylab("Sepal Width") + 
  ggtitle("Scatterplot with smoothers") + 
  geom_smooth(method="lm")

d <- ggplot(data=iris, aes(Sepal.Length, y=Sepal.Width, color=Species)) + 
  geom_point(aes(shape=Species), size=1.5) + 
  geom_smooth(method="lm") +
  xlab("Sepal Length") + 
  ylab("Sepal Width") +
  facet_grid(. ~ Species)

e <- ggplot(data=iris, aes(x = Sepal.Length))+ 
  stat_density(aes(ymax = ..density..,  ymin = -..density.., 
                       fill = Species, color = Species), 
                   geom = "ribbon", position = "identity") +
  facet_grid(. ~ Species) + coord_flip() + xlab("Sepal Length") 


# Look at image together---------------------------------------------
ggpubr::ggarrange(a,b,c,d,e,
                  nrow = 3, ncol = 2,
                  common.legend = TRUE,
                  legend = "bottom",
                  align = "v")

# Save to file-------------------------------------------------------
ggsave(filename = "~/git_projects/pleiotropy/images/proportion_hits_original_vs_mean_permuted_log.png",
       plot = ggpubr::ggarrange(pn, pg, m, e,
                                nrow = 2, ncol = 2, 
                                common.legend = TRUE, 
                                legend = "bottom", 
                                align = "v"),
       width = 7.5,
       height = 6,
       units = "in",
       dpi = "retina"
)

