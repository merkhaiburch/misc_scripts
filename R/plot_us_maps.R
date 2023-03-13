# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-03-09
# Updated... 2023-03-09
#
# Description:
# Make a US map plot that colors where I have lived
# https://cran.r-project.org/web/packages/usmap/vignettes/mapping.html
# ------------------------------------------------------------------------------



# Load package
library(usmap)
library(ggplot2)

# Plot the presence absence variation of me
mer_states <- data.frame(state = c("HI", "AZ", "SD", "CA", "ID", "NV", "OR", 
                                   "WA", "MT", "WY", "CO", "NM", "UT", "ND", 
                                   "NE", "KA", "OK", "TX", "KS"),
                         mer_presabs = c(1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))

# Select certain states
a <- plot_usmap(include = c("CA", "ID", "NV", "OR", "WA", "HI", 
                       "MT", "WY", "CO", "NM", "UT", "AZ", "ND",
                       "SD", "NE", "KA", "OK", "TX", "KS"),
           data = mer_states,
           values = "mer_presabs",
           labels = TRUE) +
  scale_fill_continuous(low = "white", high = "#009E73", label = scales::comma) + 
  theme(legend.position = "null")

ggsave(filename = "~/Downloads/states.png", plot = a)
