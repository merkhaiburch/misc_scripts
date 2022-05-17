# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-05-14
# Updated... 2022-05-14

# Description 
# Create a helper function to plot GWAS results  in a manhattan
# plot and qq plot
# Manhattan plot:
# https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
#
# QQ Plot
# https://www.broadinstitute.org/files/shared/diabetes/scandinavs/qqplot.R
# ---------------------------------------------------------------

plot_gwas_qq <- function(gwas_df, plot_title, plot_id) {
  
  # Format the manhattan plot data
  don <- gwas_df %>% 
    group_by(Chr) %>% 
    summarise(chr_len=max(Pos)) %>% 
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    left_join(gwas_df, ., by=c("Chr"="Chr")) %>%
    arrange(Chr, Pos) %>%
    mutate(BPcum=Pos+tot) %>%
    dplyr::filter(-log10(p)>1.0)
  
  # Create different colors for chromosomes
  don <- don %>% mutate(chr_color = ifelse(Chr%%2, "a" ,"b"))
  
  # Set centers of chroms
  axis_set <- don %>% 
    group_by(Chr) %>% 
    summarize(center = mean(BPcum))
  
  # Make the Manhattan plot on the gwasResults dataset
  a <- ggplot(don, aes(x = BPcum, y = -log10(p))) +
    geom_point(aes(colour = chr_color), size = 1) +
    scale_color_manual(values=c("grey28","grey61")) +
    geom_hline(yintercept = -log10(0.00001), color = "black", linetype = "solid") +
    labs(x = NULL, y = "-log10(p)", title = "y ~ SNP + management") + 
    ggtitle(plot_title) +
    theme(legend.position="none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    scale_x_continuous(label = axis_set$Chr, breaks = axis_set$center)
  
  # # Make QQ plot
  #https://www.broadinstitute.org/files/shared/diabetes/scandinavs/qqplot.R
  observed <- sort(gwas_df$p)
  lobs <- -(log10(observed))
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  qq_data <- data.frame(expected, lexp, lobs)
  b <- ggplot(qq_data, aes(x = lexp, y = lobs)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0) +
    ylim(0,15) +
    ggtitle("QQ Plot") +
    labs(x = "-log10(Expected)", y = "-log10(Observed)") 
  
  # export the plot
  ggsave(paste0("~/Box Sync/Cornell_PhD/labProjects/debugging/", plot_id, ".png"),
         plot = ggpubr::ggarrange(a, b, nrow = 1, ncol = 2),
         width = 10,
         height = 5.5, 
         units = "in",
         dpi = "retina")
}
