# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-04-18
# Updated... 2022-04-18
#
# Description 
# Calculating the number of intervals 'under the curve'
# Code modified from: 24_plot_pleiotropy.R
# ---------------------------------------------------------------

# Load in packages
library(magrittr)
library(dplyr)
library(ggplot2)


# -------------------------------
# Data gathering & formatting
# -------------------------------

# read in data aggregatied by aggregate_interval.data.R
data_dir <- "/Volumes/merrittData1/pleiotropy/interval_data"
data_dir <- "~/Downloads/"

# Load in melted and unmelted data
all_interval_data_melted <- data.table::fread(paste0(data_dir, "/all_interval_data_melted.csv"))

# subset to just columns and rows needed for plotting
all_interval_data_melted <- all_interval_data_melted %>% 
  filter(is_adjusted == "Unadjusted")


# ---------------------------------------------------------------------------
# Calculate the difference in proportion of trait counts between the 
# permuted and non-permuted data. Aggregating across genic/intergenic intervals
# ---------------------------------------------------------------------------

# Create a function that counts the number of intervals with 0, 1, >2, and >0 difference trait counts
curve_counter <- function(df, name_of_pop) {
  # Count the number of intervals with X number of traits mapping to them
  interval_proportion <- df %>% group_by(pleiotropy_score, perm_type) %>% summarise(n = n()) 
  
  # Count the total number of intervals
  lala <- df %>% group_by(perm_type) %>% summarise(n = n())
  
  # Combine files
  prop_intervals <- merge(interval_proportion, lala, by = c("perm_type"))
  
  # Calculate the proportion of X intervals in a count / total number of intervals
  prop_intervals$proportion_pleio <- prop_intervals$n.x/prop_intervals$n.y

  # Stupid idea to make data for lab meeting, split  oof and perm data, merge together
  og <- prop_intervals %>% filter(perm_type == "Original Data")
  perm <- prop_intervals %>% filter(perm_type == "Permuted")
  diff_intervals <- merge(og, perm, by = "pleiotropy_score", all = TRUE)

  # Do some formatting
  diff_intervals <- diff_intervals %>% mutate_at(c(5,9), ~replace(., is.na(.), 0)) # replace NA with zero (in instances where permuted data has no counts at value of #unique)
  
  # Calculate the difference in the proportion of intervals
  diff_intervals$difference <- diff_intervals$proportion_pleio.x-diff_intervals$proportion_pleio.y

  # Check sum of og data proportions, should equal one
  print(sum(diff_intervals$proportion_pleio.x))

  # number of intervals at certain counts
  a <- prop_intervals %>% 
    filter(pleiotropy_score == 0 & perm_type == "Original Data") %>% 
    select(proportion_pleio, n.x) # 0 traits
  
  # intervals with 2-25 traits
  b <- diff_intervals %>% 
    select(difference, proportion_pleio.x, pleiotropy_score, n.x.x) %>%
    filter(difference > 0 & pleiotropy_score >= 2 & pleiotropy_score <= 25) %>% 
    summarise(sum = sum(proportion_pleio.x), n = sum(n.x.x))
  
  # Intervals with 26-50 traits
  c <- diff_intervals %>% 
    select(difference, proportion_pleio.x, pleiotropy_score, n.x.x) %>%
    filter(difference > 0 & pleiotropy_score >= 26 & pleiotropy_score <= 50) %>% 
    summarise(sum = sum(proportion_pleio.x), n = sum(n.x.x))
  
  # intervals with 51-100
  d <- diff_intervals %>% 
    select(difference, proportion_pleio.x, pleiotropy_score, n.x.x) %>%
    filter(difference > 0 & pleiotropy_score >= 51 & pleiotropy_score <= 100) %>% 
    summarise(sum = sum(proportion_pleio.x), n = sum(n.x.x))
  
  # all intervals with positive difference
  e <- diff_intervals %>% 
    select(difference, proportion_pleio.x, pleiotropy_score, n.x.x) %>%
    filter(difference > 0) %>% 
    summarise(sum = sum(proportion_pleio.x), n = sum(n.x.x)) 

  # aggregate metrics
  props <- cbind(name_of_pop, a, b, c, d, e) 
  colnames(props) <- c("population", "zero", "zero_n",
                       "two_twentyfive", "two_twentyfive_n",
                       "twentysix_fifty", "twentysix_fifty_n", 
                       "fiftyone_onehundred", "fiftyone_onehundred_n",
                       "all_above_0_difference", "all_above_0_difference_n")
  
  return(props)
}

# Use function ---------------------------------------------

# Field NAM
nf <- all_interval_data_melted %>% filter(Population == "NAM")
nf_res <- curve_counter(nf, "NAM Field")

# Field Goodman
gf <- all_interval_data_melted %>% filter(Population == "Goodman" & trait_type == "Physiological")
gf_res <- curve_counter(gf, "Goodman Field")

# Metabolite
gm <- all_interval_data_melted %>% filter(Population == "Goodman" & trait_type == "Metabolite")
gm_res <- curve_counter(gm, "Goodman Mass Features")

# Expression
ge <- all_interval_data_melted %>% filter(Population == "Goodman" & trait_type == "Expression")
ge_res <- curve_counter(ge, "Goodman Expression")

# Combine all
tmp <- rbind(nf_res, gf_res, gm_res, ge_res)

# Write to file
write.csv(tmp, "~/git_projects/pleiotropy/data/interval_count.csv", quote = FALSE, row.names = FALSE)

