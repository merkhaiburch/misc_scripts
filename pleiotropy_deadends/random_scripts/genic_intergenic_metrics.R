# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-09-27 
#
# Description 
#   - Get metics for genic and intergenic ranges
#   - Number of ranges, size, the split between genic & intergenic
# ---------------------------------------------------------------

# Helpful packages
library(magrittr)
library(dplyr)

# load in ranges
ranges <- read.csv("git_projects/pleiotropy/data/genic_intergenic_intervals_b73_v4.49.csv")

# Number of ranges --> 
dim(ranges)

# ranges length
ranges$length <- ranges$end - ranges$start
summary(ranges$length)

# count of type of range
ranges$range_type <- gsub("_.*", "", ranges$rr_id)
ranges %>% 
  group_by(range_type) %>% 
  summarize(n())
  
ranges %>% 
  group_by(range_type) %>% 
  summarize(mean_size = mean(length, na.rm = TRUE))
