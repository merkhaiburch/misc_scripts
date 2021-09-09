# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-01-08 
#
# Description 
#   - Standard function to plot a ggplot Manhattan plot
# https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
# ---------------------------------------------------------------

# Load Packages
library(dplyr)

# Make a function
manhattan_plot <- function(gwas_result, sig_threshold_1, sig_threshold_2, ylim, title){
  # Sort dataframe by chromosome
  gwas_result <- gwas_result[order(gwas_result$Chr),]
  
  # Make x axis values
  nCHR <- length(unique(gwas_result$Chr))
  gwas_result$BPcum <- NA
  s <- 0
  nbp <- c()
  for (i in unique(gwas_result$Chr)){
    nbp[i] <- max(gwas_result[gwas_result$Chr == i,]$Pos)
    gwas_result[gwas_result$Chr == i,"BPcum"] <- gwas_result[gwas_result$Chr == i,"Pos"] + s
    s <- s + nbp[i]
  }
  
  # Format x axis to plot by chromosome
  axis.set <- gwas_result %>% 
    group_by(Chr) %>% 
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  
  # Plot everything out
  plot <-  ggplot(gwas_result, aes(x = BPcum, y = -log10(p), 
                                   color = as.factor(Chr), size = -log10(p))) +
    geom_point(size = 1) +
    geom_hline(yintercept = sig_threshold_1, color = "grey40", linetype = "dashed") +
    geom_hline(yintercept = sig_threshold_2, color = "black", linetype = "dashed") +
    scale_x_continuous(label = axis.set$Chr, breaks = axis.set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
    # scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, y = "-log10(p)", title = title) + 
    # ggtitle(title) +
    theme_minimal() +
    theme( 
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size = 9, vjust = 0.5),
      axis.text.y = element_text(size = 9, vjust = 0.5))
  return(plot)
}



