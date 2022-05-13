# ---------------------------------------------------------------------------
# 1) Plot un-adjusted, splitting by genic or intergenic intervals
# log(Observed data/mean permuted data count)
# ---------------------------------------------------------------------------
all_interval_data <- data.table::fread(paste0(data_dir, "/all_interval_data_unique_filter.csv"))


# NOT LOG TRANSFORMED ------------------------------------------------------------------------------

# Physiological NAM-------------------------------------------------
temp <- all_interval_data %>% select("rr_id", "range_type", "start", "seqid", "nam_filtered_all_uniqueCount_2", 
                                     "npp_unique_perm_1","npp_unique_perm_2","npp_unique_perm_3",
                                     "npp_unique_perm_4","npp_unique_perm_5","npp_unique_perm_6",
                                     "npp_unique_perm_7","npp_unique_perm_8","npp_unique_perm_9",
                                     "npp_unique_perm_10")
dim(temp)

# Take mean over all permutations
temp$perm_mean <- rowMeans(temp[,6:15])

# Calculate the log(observed/permuted mean)
temp$scaled_pleio_mean <- (temp[,5] + 1)/(temp$perm_mean + 1)

temp %>% group_by(range_type) %>% summarise(m = mean(scaled_pleio_mean), med = median(scaled_pleio_mean), n = n())

# plot
pn_nolog <- ggplot(data=temp, aes(x = scaled_pleio_mean,  fill = range_type)) +
  geom_histogram(size=0.5, bins = 100, alpha = 0.4) +
  xlab("(Observed + 1)/(Mean Permuted Trait Count + 1)") +
  ylab("Count of Intervals")+
  ggtitle("NAM Field") +
  theme_classic()+
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  scale_colour_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette) +
  labs(tag = "A") +
  geom_vline(xintercept = 1)

pn_man

ggplot(temp, aes(x = start, y = scaled_pleio_mean, color = as.factor(seqid)))+
  geom_point() +
  geom_hline(yintercept =1) +
  theme_classic() +
  facet_grid(.~seqid, scales = "free") +
  xlab("Position") +
  ylab("Number of Traits Mapped to Interval")+
  ggtitle("NAM Field: Observed Data") +
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size),
        axis.text.x = element_text(angle = 45, hjust=1),
        plot.tag=element_text(size=tag_size),
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm')) +
  labs(tag = "A") 

# Physiological Goodman---------------------------------------------
temp <- all_interval_data %>% select("rr_id", "range_type", "goodman_filtered_physiological_uniqueCount_2", 
                                     "gpp_unique_perm_1","gpp_unique_perm_2","gpp_unique_perm_3",
                                     "gpp_unique_perm_4","gpp_unique_perm_5","gpp_unique_perm_6",
                                     "gpp_unique_perm_7","gpp_unique_perm_8","gpp_unique_perm_9",
                                     "gpp_unique_perm_10")
dim(temp)

# Take mean over all permutations
temp$perm_mean <- rowMeans(temp[,6:15])

# Calculate the log(observed/permuted mean)
temp$scaled_pleio_mean <- (temp[,5] + 1)/(temp$perm_mean + 1)

# plot
pg_nolog <- ggplot(data=temp, aes(x = scaled_pleio_mean, fill = range_type)) +
  geom_histogram(size=0.5, bins = 100, alpha = 0.4) +
  xlab("(Observed + 1)/(Mean Permuted Trait Count + 1)") +
  ylab("Count of Intervals")+
  ggtitle("Goodman Field") +
  theme_classic()+
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  scale_colour_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette) +
  labs(tag = "B") +
  geom_vline(xintercept = 1)

# Metabolite-------------------------------------------------------
temp <- all_interval_data %>% select("rr_id", "range_type", "goodman_filtered_metabolite_uniqueCount", 
                                     "gpm_perm_1","gpm_perm_2","gpm_perm_3",
                                     "gpm_perm_4","gpm_perm_5","gpm_perm_6",
                                     "gpm_perm_7","gpm_perm_8","gpm_perm_9",
                                     "gpm_perm_10")
dim(temp)

# Take mean over all permutations
temp$perm_mean <- rowMeans(temp[,6:15])

# Calculate the log(observed/permuted mean)
temp$scaled_pleio_mean <- (temp[,5] + 1)/(temp$perm_mean + 1)

# plot
m_nolog <- ggplot(data=temp, aes(x = scaled_pleio_mean, fill = range_type)) +
  geom_histogram(size=0.5, bins = 100, alpha = 0.4) +
  xlab("(Observed + 1)/(Mean Permuted Trait Count + 1)") +
  ylab("Count of Intervals")+
  ggtitle("Goodman Mass Features") +
  theme_classic()+
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  scale_colour_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette) +
  labs(tag = "C") +
  geom_vline(xintercept = 1)

# Expression--------------------------------------------------------
temp <- all_interval_data %>% select("rr_id", "range_type", "goodman_filtered_expression_uniqueCount", 
                                     "gpe_perm_1","gpe_perm_2","gpe_perm_3",
                                     "gpe_perm_4","gpe_perm_5")
dim(temp)

# Take mean over all permutations
temp$perm_mean <- rowMeans(temp[,6:15])

# Calculate the log(observed/permuted mean)
temp$scaled_pleio_mean <- (temp[,5] + 1)/(temp$perm_mean + 1)

# plot
e_nolog <- ggplot(data=temp, aes(x = scaled_pleio_mean, fill = range_type)) +
  geom_histogram(size=0.5, bins = 100, alpha = 0.4) +
  xlab("(Observed + 1)/(Mean Permuted Trait Count + 1)") +
  ylab("Count of Intervals")+
  ggtitle("Goodman Expression") +
  theme_classic()+
  theme(axis.text=element_text(size=axis_text_size), 
        axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=title_text_size), 
        legend.position="bottom", 
        legend.text = element_text(size = legend_text_size),
        legend.title=element_blank(),
        legend.key.size = unit(legend_shape_size, 'cm'),
        plot.tag=element_text(size=tag_size)) +
  scale_colour_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette) +
  labs(tag = "D")  +
  geom_vline(xintercept = 1)

# Look at image together---------------------------------------------
ggpubr::ggarrange(pn_nolog, pg_nolog, m_nolog, e_nolog,
                  nrow = 2, ncol = 2,
                  common.legend = TRUE,
                  legend = "bottom",
                  align = "v")

# Save to file-------------------------------------------------------
ggsave(filename = "~/git_projects/pleiotropy/images/.png",
       plot = ggpubr::ggarrange(pn_nolog, pg_nolog, m_nolog, e_nolog,
                                nrow = 2, ncol = 2,
                                common.legend = TRUE,
                                legend = "bottom",
                                align = "v"),
       width = 7.5,
       height = 6,
       units = "in",
       dpi = "retina"
)


