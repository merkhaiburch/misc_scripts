# Garbage individual TE moodel plots. Thrown out because plotting at the 
# superfamily and class level only really reveals their annotation method
# and not necessarily biology


# Plot all families with no filtering ------------------------------------------
noFilterPlot <- ggplot(allModelEffects, aes(x = Term, y = estimate, fill = Model)) +
  geom_boxplot() +
  xlab("Distance from Gene") +
  ylab("Effect") +
  ggtitle(paste0("Unfiltered TE set with: ", length(unique(allModelEffects$family)), " families")) +
  theme_bw() +
  facet_wrap(vars(Classification), scales = "free_y") +
  theme(axis.text=element_text(size=axis_text_size),
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size),
        legend.position="bottom",
        legend.text = element_text(size = axis_text_size),
        legend.title=element_text(size = legend_text_size),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_colour_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette)
# ggsave("/home/mbb262/git_projects/te_ase_nam/figs/noFilter_jointSingleModelEffects.png",
#        plot = noFilterPlot,
#        width = 10,
#        height = 7.5,
#        units = "in",
#        dpi = "retina")



# Filter out problem family and re-plot ----------------------------------------
# The menace family: FAM14838
temp <- sub_allModelEffects %>% filter(family != "FAM14838")
menace <- ggplot(temp, aes(x = Term, y = estimate, fill = Model)) +
  geom_boxplot() +
  xlab("Distance from Gene") +
  ylab("Effect") +
  ggtitle(paste0("Filtered TE set with: ", length(unique(temp$family)), " families - FAM14838 removed")) +
  theme_bw() +
  facet_wrap(vars(Classification), scales = "free_y") +
  theme(axis.text=element_text(size=axis_text_size),
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size),
        legend.position="bottom",
        legend.text = element_text(size = axis_text_size),
        legend.title=element_text(size = legend_text_size),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_colour_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette)
# ggsave("/home/mbb262/git_projects/te_ase_nam/figs/filter2_jointSingleModelEffects.png",
#        plot = menace,
#        width = 10,
#        height = 7.5,
#        units = "in",
#        dpi = "retina")



# Plot class specific effects --------------------------------------------------
# Not recommended by Michelle as this mostly reflects annotation styles and
# not actual TE biology

temp <- sub_allModelEffects %>% filter(family != "FAM14838")
classy <- temp %>%
  left_join(temp %>% group_by(Classification) %>% summarise(N=n())) %>%
  mutate(Label=paste0(Classification,' (n = ',N,')')) %>%
  ggplot(., mapping = aes(x = Term, y = estimate, fill = Model)) +
  geom_boxplot() +
  xlab("Distance from Gene") +
  ylab("Effect") +
  ggtitle(paste0("Filtered TE set with: ", length(unique(temp$family)), " families - FAM14838 removed")) +
  theme_bw() +
  facet_wrap(vars(Label), scales = "free_y") +
  theme(axis.text=element_text(size=axis_text_size),
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size),
        legend.position="bottom",
        legend.text = element_text(size = axis_text_size),
        legend.title=element_text(size = legend_text_size),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_colour_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette)
# ggsave("/home/mbb262/git_projects/te_ase_nam/figs/filter2_classification_jointSingleModelEffects.png",
#        plot = classy,
#        width = 17,
#        height = 10,
#        units = "in",
#        dpi = "retina")

# Plot R2 against TE family copy number and split by class ---------------------

# Make a threshold line
temp <- quantile(subjoint$adj_r2, .99)
temp

# Facet results out by classification
facetsuperfam <- ggplot(subjoint, aes(x = bp, y = adj_r2)) +
  geom_point() +
  xlab("Total TE Family Base Pairs") +
  ylab("Adjusted R2") +
  facet_wrap(vars(Classification), scales = "free_x") +
  geom_hline(yintercept = temp) +
  theme_bw() +
  theme(axis.text=element_text(size=axis_text_size),
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size),
        legend.position="bottom",
        legend.text = element_text(size = axis_text_size),
        legend.title=element_text(size = legend_text_size),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  scale_colour_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette)
facetsuperfam

ggsave("/home/mbb262/git_projects/te_ase_nam/figs/facet_joint_modelR2_filteredByTEPrevalance_totalTEBP.png",
       plot = facetsuperfam,
       width = 8.5,
       height = 7.5, 
       units = "in",
       dpi = "retina")


# Hybrid ------------------------------------------------------------------------------

# Plot class specific effects --------------------------------------------------
# Not recommended by Michelle as this mostly reflects annotation styles and
# not actual TE biology

temp <- sub_allModelEffects %>% filter(family != "FAM14838")
classy <- temp %>% 
  left_join(temp %>% group_by(Classification) %>% summarise(N=n())) %>%
  mutate(Label=paste0(Classification,' (n = ',N,')')) %>%
  ggplot(., mapping = aes(x = Term, y = estimate, fill = Model)) +
  geom_boxplot() +
  xlab("Distance from Gene") +
  ylab("Effect") +
  ggtitle(paste0("Filtered TE set with: ", length(unique(temp$family)), " families - FAM14838 removed")) +
  theme_bw() +
  facet_wrap(vars(Label), scales = "free_y") +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = axis_text_size),
        legend.title=element_text(size = legend_text_size),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_colour_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette)
# ggsave("/home/mbb262/git_projects/te_ase_nam/figs/hybrid_filter2_classification_jointSingleModelEffects.png",
#        plot = classy,
#        width = 17,
#        height = 10, 
#        units = "in",
#        dpi = "retina")
