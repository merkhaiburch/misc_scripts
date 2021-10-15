# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-08-26
# Updated... 2021-10-15
# Description 
#   - calculates Adjusted Entropy and fits it to gerp counts
# ---------------------------------------------------------------

# Followed this tutorial to pass variables into functions:
# https://aosmith.rbind.io/2019/06/24/function-for-model-fitting/

# -----------------------------------------------
#     Plot all gerp levels with points and a line
# -----------------------------------------------


# ----------------------------
# nam physiological - filtered
all_gerp_levels1 <-ggplot(npf_plot, aes(x = gerp_count, y = fitted_entropy, color = Percentile)) +
  geom_point(alpha = .3) +
  geom_smooth(method = lm) + 
  labs(x = "Log Count of GERP Sites", 
       y = "Adjusted Entropy", 
       title = "NAM Physiological Filtered") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size), 
        legend.position="none") +
  ylim(-6, 7.5)

# ----------------------------
# nam physiological - random
all_gerp_levels12 <-ggplot(npr_plot, aes(x = gerp_count, y = fitted_entropy, color = Percentile)) +
  geom_point(alpha = .3) +
  geom_smooth(method = lm) + 
  labs(x = "Log Count of GERP Sites", 
       y = "Adjusted Entropy", 
       title = "NAM Physiological random") +
  scale_color_manual(values=cbbbPalette) +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size), 
        legend.position="none") +
  ylim(-6, 7.5)


# ----------------------------------
# goodman physiological - filtered
all_gerp_levels2 <- ggplot(gpf_plot, aes(x = gerp_count, y = fitted_entropy, color = Percentile)) +
  geom_point(alpha = .3) +
  geom_point(alpha = .5) +
  geom_smooth(method = lm) + 
  labs(x = "Log Count of GERP Sites", 
       y = "Adjusted Entropy", 
       title = "Goodman Physiological Filtered") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size), 
        legend.position="none") +
  ylim(-6, 7.5) +
  scale_color_manual(values=cbbbPalette)

# ----------------------------------
# goodman physiological - random
all_gerp_levels22 <- ggplot(gpr_plot, aes(x = gerp_count, y = fitted_entropy, color = Percentile)) +
  geom_point(alpha = .3) +
  geom_point(alpha = .5) +
  geom_smooth(method = lm) + 
  labs(x = "Log Count of GERP Sites", 
       y = "Adjusted Entropy", 
       title = "Goodman Physiological Random") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size), 
        legend.position="none") +
  ylim(-6, 7.5) +
  scale_color_manual(values=cbbbPalette)

# ------------------------------
# goodman metabolite - filtered
all_gerp_levels3 <- ggplot(gmf_plot, aes(x = gerp_count, y = fitted_entropy, color = Percentile)) +
  geom_point(alpha = .3) +
  geom_point(alpha = .5) +
  geom_smooth(method = lm) + 
  labs(x = "Log Count of GERP Sites", 
       y = "Adjusted Entropy", 
       title = "Goodman Metabolite Filtered") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size), 
        legend.position="none") +
  ylim(-6, 7.5) +
  scale_color_manual(values=cbbbPalette)

# ------------------------------
# goodman metabolite - random
all_gerp_levels32 <- ggplot(gmr_plot, aes(x = gerp_count, y = fitted_entropy, color = Percentile)) +
  geom_point(alpha = .3) +
  geom_point(alpha = .5) +
  geom_smooth(method = lm) + 
  labs(x = "Log Count of GERP Sites", 
       y = "Adjusted Entropy", 
       title = "Goodman Metabolite Random") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size), 
        legend.position="none") +
  ylim(-6, 7.5) +
  scale_color_manual(values=cbbbPalette)


# ---------------------------------
# goodman expression - filtered
all_gerp_levels4 <- ggplot(gef_plot, aes(x = gerp_count, y = fitted_entropy, color = Percentile)) +
  geom_point(alpha = .3) +
  geom_point(alpha = .5) +
  geom_smooth(method = lm) + 
  labs(x = "Log Count of GERP Sites", 
       y = "Adjusted Entropy", 
       title = "Goodman Expression Filtered") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size), 
        legend.position="bottom", 
        legend.title = element_text(size = axis_text_size_legend), 
        legend.text = element_text(size = axis_text_size_legend)) +
  scale_color_manual(values=cbbbPalette)+
  ylim(-6, 7.5) 

# ---------------------------------
# goodman expression - random
all_gerp_levels42 <- ggplot(ger_plot, aes(x = gerp_count, y = fitted_entropy, color = Percentile)) +
  geom_point(alpha = .3) +
  geom_point(alpha = .5) +
  geom_smooth(method = lm) + 
  labs(x = "Log Count of GERP Sites", 
       y = "Adjusted Entropy", 
       title = "Goodman Expression Random") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size), 
        legend.position="bottom", 
        legend.title = element_text(size = axis_text_size_legend), 
        legend.text = element_text(size = axis_text_size_legend)) +
  scale_color_manual(values=cbbbPalette) +
  ylim(-6, 7.5)

ggsave("/workdir/mbb262/interval_data/adjusted_entropy_all_gerp_levels.png",
       plot = ((all_gerp_levels1 + all_gerp_levels12)/(all_gerp_levels2 + all_gerp_levels22)/(all_gerp_levels3 + all_gerp_levels32)/(all_gerp_levels4 + all_gerp_levels42)),
       width = 9,
       height = 16,
       units = "in")



# -----------------------------------------------
#  Plot all percentiles without total gerp sites 
# points and line
# -----------------------------------------------



# ----------------------------
# nam physiological - filtered
gerp_levels_minusTotalSites1 <-ggplot(npf_plot_noTotalSites, aes(x = gerp_count, y = fitted_entropy, color = Percentile)) +
  geom_point(alpha = .3) +
  geom_smooth(method = lm) + 
  labs(x = "Log Count of GERP Sites", 
       y = "Adjusted Entropy", 
       title = "NAM Physiological Filtered") +
  scale_color_manual(values=cbbbPalette) +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size), 
        legend.position="none") +
  ylim(-6, 7.5)

# ----------------------------
# nam physiological - random
gerp_levels_minusTotalSites12 <-ggplot(npr_plot_noTotalSites, aes(x = gerp_count, y = fitted_entropy, color = Percentile)) +
  geom_point(alpha = .3) +
  geom_smooth(method = lm) + 
  labs(x = "Log Count of GERP Sites", 
       y = "Adjusted Entropy", 
       title = "NAM Physiological Random") +
  scale_color_manual(values=cbbbPalette) +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size), 
        legend.position="none") +
  ylim(-6, 7.5)


# ----------------------------------
# goodman physiological - filtered
gerp_levels_minusTotalSites2 <- ggplot(gpf_plot_noTotalSites, aes(x = gerp_count, y = fitted_entropy, color = Percentile)) +
  geom_point(alpha = .3) +
  geom_point(alpha = .5) +
  geom_smooth(method = lm) + 
  labs(x = "Log Count of GERP Sites", 
       y = "Adjusted Entropy", 
       title = "Goodman Physiological Filtered") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size), 
        legend.position="none") +
  ylim(-6, 7.5) +
  scale_color_manual(values=cbbbPalette)

# ----------------------------------
# goodman physiological - random
gerp_levels_minusTotalSites22 <- ggplot(gpr_plot_noTotalSites, aes(x = gerp_count, y = fitted_entropy, color = Percentile)) +
  geom_point(alpha = .3) +
  geom_point(alpha = .5) +
  geom_smooth(method = lm) + 
  labs(x = "Log Count of GERP Sites", 
       y = "Adjusted Entropy", 
       title = "Goodman Physiological Random") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size), 
        legend.position="none") +
  ylim(-6, 7.5) +
  scale_color_manual(values=cbbbPalette)

# ------------------------------
# goodman metabolite - filtered
gerp_levels_minusTotalSites3 <- ggplot(gmf_plot_noTotalSites, aes(x = gerp_count, y = fitted_entropy, color = Percentile)) +
  geom_point(alpha = .3) +
  geom_point(alpha = .5) +
  geom_smooth(method = lm) + 
  labs(x = "Log Count of GERP Sites", 
       y = "Adjusted Entropy", 
       title = "Goodman Metabolite Filtered") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size), 
        legend.position="none") +
  ylim(-6, 7.5) +
  scale_color_manual(values=cbbbPalette)

# ------------------------------
# goodman metabolite - random
gerp_levels_minusTotalSites32 <- ggplot(gmr_plot_noTotalSites, aes(x = gerp_count, y = fitted_entropy, color = Percentile)) +
  geom_point(alpha = .3) +
  geom_point(alpha = .5) +
  geom_smooth(method = lm) + 
  labs(x = "Log Count of GERP Sites", 
       y = "Adjusted Entropy", 
       title = "Goodman Metabolite Random") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size), 
        legend.position="none") +
  ylim(-6, 7.5) +
  scale_color_manual(values=cbbbPalette)


# ---------------------------------
# goodman expression - filtered
gerp_levels_minusTotalSites4 <- ggplot(gef_plot_noTotalSites, aes(x = gerp_count, y = fitted_entropy, color = Percentile)) +
  geom_point(alpha = .3) +
  geom_point(alpha = .5) +
  geom_smooth(method = lm) + 
  labs(x = "Log Count of GERP Sites", 
       y = "Adjusted Entropy", 
       title = "Goodman Expression Filtered") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size), 
        legend.position="bottom", 
        legend.title = element_text(size = axis_text_size_legend), 
        legend.text = element_text(size = axis_text_size_legend)) +
  scale_color_manual(values=cbbbPalette) +
  ylim(-6, 7.5)

# ---------------------------------
# goodman expression - filtered
gerp_levels_minusTotalSites42 <- ggplot(ger_plot_noTotalSites, aes(x = gerp_count, y = fitted_entropy, color = Percentile)) +
  geom_point(alpha = .3) +
  geom_point(alpha = .5) +
  geom_smooth(method = lm) + 
  labs(x = "Log Count of GERP Sites", 
       y = "Adjusted Entropy", 
       title = "Goodman Expression Random") +
  theme_classic() +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size), 
        legend.position="bottom", 
        legend.title = element_text(size = axis_text_size_legend), 
        legend.text = element_text(size = axis_text_size_legend)) +
  scale_color_manual(values=cbbbPalette) +
  ylim(-6, 7.5)


ggsave("/workdir/mbb262/interval_data/adjusted_entropy_gerp_levels_minusTotalSites.png",
       plot = ((gerp_levels_minusTotalSites1 + gerp_levels_minusTotalSites12)/(gerp_levels_minusTotalSites2 + gerp_levels_minusTotalSites22)/(gerp_levels_minusTotalSites3 + gerp_levels_minusTotalSites32)/(gerp_levels_minusTotalSites4 + gerp_levels_minusTotalSites42)),
       width = 9,
       height = 16, 
       units = "in")
