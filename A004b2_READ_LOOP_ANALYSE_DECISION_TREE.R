rm(list=ls())
library(dplyr)
library(tidyr)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)
library(zoo)
library(slider)
library(purrr)
library(tibble)
library(RcppRoll)


args = commandArgs(trailingOnly = FALSE)
script_path = sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0) {
  # Running via Rscript
  script_dir = dirname(normalizePath(script_path))
} else {
  # Running interactively (e.g. in RStudio), fallback to current working dir or a known path
  script_dir = dirname(rstudioapi::getActiveDocumentContext()$path)  # for RStudio
  # or script_dir = getwd()
}
setwd(script_dir)

source('~/Dropbox/Treg_problem_local/A004b1_READ_LOOP_MERGE_CORES.R')

all_data   = readRDS('all_data_merged.rds')
all_params = readRDS('all_params_merged.rds')

all_data_epithelial = all_data[c("t","realization_ind","sterile","allow_tregs_to_do_their_job",
                                 "epithelial_healthy","epithelial_inj_1","epithelial_inj_2",
                                 "epithelial_inj_3","epithelial_inj_4","epithelial_inj_5","pathogen")]



all_data_epithelial_ss = all_data_epithelial %>% filter(t>2000)
# Compute rolling standard deviation (window = windowsize time steps)
windowsize = 100
all_data_stationarity = all_data_epithelial_ss %>%
  group_by(realization_ind, sterile, allow_tregs_to_do_their_job) %>%
  arrange(t) %>%
  mutate(
    rolling_sd_healthy = roll_sd(epithelial_healthy, n = windowsize, fill = NA, align = "right"),
    rolling_sd_inj_1   = roll_sd(epithelial_inj_1,   n = windowsize, fill = NA, align = "right"),
    rolling_sd_inj_2   = roll_sd(epithelial_inj_2,   n = windowsize, fill = NA, align = "right"),
    rolling_sd_inj_3   = roll_sd(epithelial_inj_3,   n = windowsize, fill = NA, align = "right"),
    rolling_sd_inj_4   = roll_sd(epithelial_inj_4,   n = windowsize, fill = NA, align = "right"),
    rolling_sd_inj_5   = roll_sd(epithelial_inj_5,   n = windowsize, fill = NA, align = "right"),
    rolling_sd_pat     = roll_sd(pathogen, n = windowsize, fill = NA, align = "right")
  )

# For each realization, get the first time after which all 4 scenarios have stable epithelial_healthy (e.g., rolling SD < 1)
threshold_sd = 1

# Flag time points where all epithelial types are stable
all_data_stationarity = all_data_stationarity %>%
  mutate(stable_all = rolling_sd_healthy < threshold_sd &
           rolling_sd_inj_1 < threshold_sd &
           rolling_sd_inj_2 < threshold_sd &
           rolling_sd_inj_3 < threshold_sd &
           rolling_sd_inj_4 < threshold_sd &
           rolling_sd_inj_5 < threshold_sd &
           rolling_sd_pat < threshold_sd)

# Get first stable time per group
first_stable_time = all_data_stationarity %>%
  filter(stable_all) %>%
  group_by(realization_ind, sterile, allow_tregs_to_do_their_job) %>%
  summarise(first_stable_t = min(t), .groups = "drop")

# Now aggregate to the latest across 4 scenarios per realization
stationary_start = first_stable_time %>%
  group_by(realization_ind) %>%
  summarise(stationary_start = max(first_stable_t), .groups = "drop")
# Join to get stationary start per realization
all_data_poststationary = all_data %>%
  inner_join(stationary_start, by = "realization_ind") %>%
  filter(t >= stationary_start)

# saveRDS(all_data_stationarity,'all_data_stationarity.rds')
# saveRDS(all_data_poststationary,'all_data_poststationary.rds')
# 
# all_data_stationarity   = readRDS('all_data_stationarity.rds')
# all_data_poststationary = readRDS('all_data_poststationary.rds')

# keep_objects = c("all_data_poststationary") # otherwise memory issues 
# rm(list = setdiff(ls(), keep_objects))

df_sterile_1 = all_data_poststationary %>% filter(sterile == 1)
df_sterile_0 = all_data_poststationary %>% filter(sterile == 0)

# Per realization_ind: run Wilcoxon test (TRUE vs FALSE)
wilcox_results_sterile_1 = df_sterile_1 %>%
  dplyr::group_by(realization_ind) %>%
  group_split() %>%
  map_df(function(group_df) {
    df_true  = group_df %>% filter(allow_tregs_to_do_their_job == TRUE) %>% pull(epithelial_healthy)
    df_false = group_df %>% filter(allow_tregs_to_do_their_job == FALSE) %>% pull(epithelial_healthy)
    
    # Match by timepoints if you want true paired test
    common_t = intersect(group_df$t[group_df$allow_tregs_to_do_their_job == TRUE],
                          group_df$t[group_df$allow_tregs_to_do_their_job == FALSE])
    
    df_true  = group_df %>% filter(allow_tregs_to_do_their_job == TRUE, t %in% common_t) %>% arrange(t) %>% pull(epithelial_healthy)
    df_false = group_df %>% filter(allow_tregs_to_do_their_job == FALSE, t %in% common_t) %>% arrange(t) %>% pull(epithelial_healthy)
    
    tibble(
      realization_ind = unique(group_df$realization_ind),
      p_value = if(length(df_true) > 0 && length(df_false) > 0) {
        wilcox.test(df_true, df_false, paired = TRUE)$p.value
      } else { NA_real_ },
      median_diff = median(df_true) - median(df_false)
    )
  })

# Per realization_ind: run Wilcoxon test (TRUE vs FALSE)
wilcox_results_sterile_0 = df_sterile_0 %>%
  dplyr::group_by(realization_ind) %>%
  group_split() %>%
  map_df(function(group_df) {
    df_true  = group_df %>% filter(allow_tregs_to_do_their_job == TRUE) %>% pull(epithelial_healthy)
    df_false = group_df %>% filter(allow_tregs_to_do_their_job == FALSE) %>% pull(epithelial_healthy)
    
    # Match by timepoints if you want true paired test
    common_t = intersect(group_df$t[group_df$allow_tregs_to_do_their_job == TRUE],
                          group_df$t[group_df$allow_tregs_to_do_their_job == FALSE])
    
    df_true  = group_df %>% filter(allow_tregs_to_do_their_job == TRUE, t %in% common_t) %>% arrange(t) %>% pull(epithelial_healthy)
    df_false = group_df %>% filter(allow_tregs_to_do_their_job == FALSE, t %in% common_t) %>% arrange(t) %>% pull(epithelial_healthy)
    
    tibble(
      realization_ind = unique(group_df$realization_ind),
      p_value = if(length(df_true) > 0 && length(df_false) > 0) {
        wilcox.test(df_true, df_false, paired = TRUE)$p.value
      } else { NA_real_ },
      median_diff = median(df_true) - median(df_false)
    )
  })

colnames(wilcox_results_sterile_0)[2:3] = c('p_value_sterile_0','median_diff_sterile_0')
colnames(wilcox_results_sterile_1)[2:3] = c('p_value_sterile_1','median_diff_sterile_1')
wilcox_results_merged                   = merge(wilcox_results_sterile_0,wilcox_results_sterile_1,by='realization_ind')

prop_healthy_better_all_sterile_0 = all_data_poststationary %>%
  filter(sterile==0 & allow_tregs_to_do_their_job==T) %>%
  mutate(healthy_greater_than_all_injured =
           epithelial_healthy > epithelial_inj_1 &
           epithelial_healthy > epithelial_inj_2 &
           epithelial_healthy > epithelial_inj_3 &
           epithelial_healthy > epithelial_inj_4 &
           epithelial_healthy > epithelial_inj_5) %>%
  group_by(realization_ind, sterile, allow_tregs_to_do_their_job) %>%
  summarise(prop_healthy_better_all = mean(healthy_greater_than_all_injured),
            .groups = "drop")
prop_healthy_better_all_sterile_0 = prop_healthy_better_all_sterile_0[c('realization_ind','prop_healthy_better_all')]
colnames(prop_healthy_better_all_sterile_0)[2] = 'prop_healthy_better_sterile_0'

prop_healthy_better_all_sterile_1 = all_data_poststationary %>%
  filter(sterile==1 & allow_tregs_to_do_their_job==T) %>%
  mutate(healthy_greater_than_all_injured =
           epithelial_healthy > epithelial_inj_1 &
           epithelial_healthy > epithelial_inj_2 &
           epithelial_healthy > epithelial_inj_3 &
           epithelial_healthy > epithelial_inj_4 &
           epithelial_healthy > epithelial_inj_5) %>%
  group_by(realization_ind, sterile, allow_tregs_to_do_their_job) %>%
  summarise(prop_healthy_better_all = mean(healthy_greater_than_all_injured),
            .groups = "drop")
prop_healthy_better_all_sterile_1 = prop_healthy_better_all_sterile_1[c('realization_ind','prop_healthy_better_all')]
colnames(prop_healthy_better_all_sterile_1)[2] = 'prop_healthy_better_sterile_1'

prop_healthy_better   = merge(prop_healthy_better_all_sterile_1, prop_healthy_better_all_sterile_0, by='realization_ind')
wilcox_results_merged = merge(wilcox_results_merged, prop_healthy_better, by='realization_ind')

grid_size   = 25
effect_size = grid_size*0.125
prop_cutoff = 0.75

wilcox_results_merged_tregs_useful_for_sterile_0    = wilcox_results_merged %>% filter(median_diff_sterile_0>effect_size & p_value_sterile_0<0.01 & prop_healthy_better_sterile_0>prop_cutoff)
wilcox_results_merged_tregs_useful_for_sterile_1    = wilcox_results_merged %>% filter(median_diff_sterile_1>effect_size & p_value_sterile_1<0.01 & prop_healthy_better_sterile_1>prop_cutoff)
wilcox_results_merged_tregs_useful_for_sterile_both = wilcox_results_merged %>% filter((median_diff_sterile_1>effect_size & p_value_sterile_1<0.01 & prop_healthy_better_sterile_1>prop_cutoff) 
                                                                                       & (median_diff_sterile_0>effect_size & p_value_sterile_0<0.01 & prop_healthy_better_sterile_0>prop_cutoff))


100*dim(wilcox_results_merged_tregs_useful_for_sterile_0)[1]/dim(wilcox_results_merged)[1]
100*dim(wilcox_results_merged_tregs_useful_for_sterile_1)[1]/dim(wilcox_results_merged)[1]
100*dim(wilcox_results_merged_tregs_useful_for_sterile_both)[1]/dim(wilcox_results_merged)[1]

# any case that NOT having tregs is better? No, yay :)
wilcox_results_merged_tregs_bad_for_sterile_0    = wilcox_results_merged %>% filter(median_diff_sterile_0<(-1*effect_size) & p_value_sterile_0<0.01 & prop_healthy_better_sterile_0<prop_cutoff)
wilcox_results_merged_tregs_bad_for_sterile_1    = wilcox_results_merged %>% filter(median_diff_sterile_1<(-1*effect_size) & p_value_sterile_1<0.01 & prop_healthy_better_sterile_1<prop_cutoff)
wilcox_results_merged_tregs_bad_for_sterile_both = wilcox_results_merged %>% filter((median_diff_sterile_1<(-1*effect_size) & p_value_sterile_1<0.01 & prop_healthy_better_sterile_1<prop_cutoff) 
                                                                                       & (median_diff_sterile_0<(-1*effect_size) & p_value_sterile_0<0.01 & prop_healthy_better_sterile_0<prop_cutoff))


100*dim(wilcox_results_merged_tregs_bad_for_sterile_0)[1]/dim(wilcox_results_merged)[1]
100*dim(wilcox_results_merged_tregs_bad_for_sterile_1)[1]/dim(wilcox_results_merged)[1]
100*dim(wilcox_results_merged_tregs_bad_for_sterile_both)[1]/dim(wilcox_results_merged)[1]


### also filter according to pathogens, only for sterile=0 of course
pathogen_threshold = 50 # can't have a peak more than this

wilcox_results_sterile_0_pat = df_sterile_0 %>%
  group_by(realization_ind) %>%
  group_split() %>%
  map_df(function(group_df) {
    
    df_true  = group_df %>%
      filter(allow_tregs_to_do_their_job == TRUE)
    
    df_false = group_df %>%
      filter(allow_tregs_to_do_their_job == FALSE)
    
    # Match common timepoints
    common_t = intersect(df_true$t, df_false$t)
    
    df_true_matched = df_true %>%
      filter(t %in% common_t) %>%
      arrange(t)
    
    df_false_matched = df_false %>%
      filter(t %in% common_t) %>%
      arrange(t)
    
    # If either group is empty or mismatched, return NA
    if (nrow(df_true_matched) == 0 || nrow(df_false_matched) == 0 ||
        nrow(df_true_matched) != nrow(df_false_matched)) {
      return(tibble(
        realization_ind = unique(group_df$realization_ind),
        p_value = NA_real_,
        median_diff = NA_real_,
        all_pathogen_below_threshold = NA
      ))
    }
    
    # Extract pathogen vectors
    pat_true = df_true_matched$pathogen
    pat_false = df_false_matched$pathogen
    
    tibble(
      realization_ind = unique(group_df$realization_ind),
      p_value = wilcox.test(pat_true, pat_false, paired = TRUE)$p.value,
      median_diff = median(pat_true) - median(pat_false),
      all_pathogen_below_threshold = all(pat_true < pathogen_threshold)
    )
  })


# filter also according to pathogen
wilcox_results_merged_tregs_useful_for_sterile_0_pat = wilcox_results_sterile_0_pat %>% filter(median_diff<0 & p_value<0.01 & all_pathogen_below_threshold==T)

# check for all indexes?
usefull_for_everything = intersect(unique(wilcox_results_merged_tregs_useful_for_sterile_both$realization_ind),unique(wilcox_results_merged_tregs_useful_for_sterile_0_pat$realization_ind))

source("./PLOT_FUNCTIONS.R")
dir_name      = './random_simulations'
subdir_name   = '/figures'
dir.create(paste0(dir_name,subdir_name))


for (ind in usefull_for_everything){

  print(ind)
  realization_data = all_data %>% filter(realization_ind==ind)
  
  # first_stable_time_use = first_stable_time %>% filter(realization_ind == ind)
  # realization_data      = merge(realization_data, first_stable_time_use, by=c('realization_ind','sterile','allow_tregs_to_do_their_job'))
  
  p_epithelium = plot_faceted(realization_data, 
                              c("epithelial_healthy", paste0("epithelial_inj_", 1:5)),
                              "Epithelial Cell Dynamics")
  # p_resting    = plot_faceted(realization_data,
  #                             c("phagocyte_M0", "treg_resting"),
  #                             "M0 Phagocytes & Resting Tregs")
  # p_M1         = plot_faceted(realization_data,
  #                             paste0("phagocyte_M1_L_", 0:5),
  #                             "M1 Phagocyte Levels")
  # p_reg_active = plot_faceted(realization_data,
  #                             c(paste0("phagocyte_M2_L_", 0:5), "treg_active"),
                              # "M2 Phagocytes & Active Tregs")
  p_microbes   = plot_faceted(realization_data, c("commensal", "pathogen"),
                              "Microbe Dynamics")
  # p_cum_com    = plot_faceted(realization_data,
  #                             c("C_ROS", "C_M0", "C_M1", "C_M2"),
  #                             "Cumulative Commensal Deaths")
  # p_cum_path    = plot_faceted(realization_data,
  #                              c("P_ROS", "P_M0", "P_M1", "P_M2"),
  #                              "Cumulative Pathogen Deaths")
  
  # Save all 7 plots with realization_ind in filename
  ggsave(paste0(dir_name,subdir_name,"/epithelium_", ind,"_useful_both.png"), p_epithelium, dpi = 300, width = 16, height = 10, bg = 'white')
  # ggsave(paste0(dir_name,subdir_name,"/resting_tregs_", ind, "_useful_both.png"), p_resting, dpi = 300, width = 8, height = 5, bg = 'white')
  # ggsave(paste0(dir_name,subdir_name,"/phagocyte_M1_", ind, "_useful_both.png"), p_M1, dpi = 300, width = 8, height = 5, bg = 'white')
  # ggsave(paste0(dir_name,subdir_name,"/phagocyte_M2_activeTreg_", ind, "_useful_both.png"), p_reg_active, dpi = 300, width = 8, height = 5, bg = 'white')
  ggsave(paste0(dir_name,subdir_name,"/microbes_", ind, "_useful_both.png"), p_microbes, dpi = 300, width = 8, height = 5, bg = 'white')
  # ggsave(paste0(dir_name,subdir_name,"/cumulative_commensal_deaths_", ind, "_useful_both.png"), p_cum_com, dpi = 300, width = 8, height = 5, bg = 'white')
  # ggsave(paste0(dir_name,subdir_name,"/cumulative_pathogen_deaths_", ind, "_useful_both.png"), p_cum_path, dpi = 300, width = 8, height = 5, bg = 'white')
}


for (ind in unique(wilcox_results_merged_tregs_bad_for_sterile_0$realization_ind)){
  
  print(ind)
  realization_data = all_data %>% filter(realization_ind==ind)

  p_epithelium = plot_faceted(realization_data, 
                              c("epithelial_healthy", paste0("epithelial_inj_", 1:5)),
                              "Epithelial Cell Dynamics")
  p_microbes   = plot_faceted(realization_data, c("commensal", "pathogen"),
                              "Microbe Dynamics")
  # Save all 7 plots with realization_ind in filename
  ggsave(paste0(dir_name,subdir_name,"/epithelium_", ind,"_bad_for_sterile_0.png"), p_epithelium, dpi = 300, width = 16, height = 10, bg = 'white')
  ggsave(paste0(dir_name,subdir_name,"/microbes_", ind, "_bad_for_sterile_0.png"), p_microbes, dpi = 300, width = 8, height = 5, bg = 'white')
}

for (ind in unique(wilcox_results_merged_tregs_bad_for_sterile_1$realization_ind)){
  
  print(ind)
  realization_data = all_data %>% filter(realization_ind==ind)
  
  p_epithelium = plot_faceted(realization_data, 
                              c("epithelial_healthy", paste0("epithelial_inj_", 1:5)),
                              "Epithelial Cell Dynamics")
  p_microbes   = plot_faceted(realization_data, c("commensal", "pathogen"),
                              "Microbe Dynamics")
  # Save all 7 plots with realization_ind in filename
  ggsave(paste0(dir_name,subdir_name,"/epithelium_", ind,"_bad_for_sterile_1.png"), p_epithelium, dpi = 300, width = 16, height = 10, bg = 'white')
  ggsave(paste0(dir_name,subdir_name,"/microbes_", ind, "_bad_for_sterile_1.png"), p_microbes, dpi = 300, width = 8, height = 5, bg = 'white')
}

### What are the parameters?

all_params = readRDS('all_params_merged.rds')

all_params_useful = all_params %>% filter(realization_ind %in% usefull_for_everything)

##### Check PCA analysis?

library(rpart)
library(rpart.plot)

all_params$useful = ifelse(all_params$realization_ind %in% usefull_for_everything, "Useful", "Not Useful")
all_params$useful = factor(all_params$useful, levels = c("Useful", "Not Useful"))

# #Train tree
# #Config below
# #Treat Useful/Not Useful fairly (equal weight)
# #Grow a detailed tree (low pruning threshold)
tree_model = rpart(
  useful ~ ., 
  data = all_params %>% dplyr::select(-realization_ind),
  method = "class",
  parms = list(prior = c(Useful = 0.5, `Not Useful` = 0.5)),
  control = rpart.control(minsplit = 2, cp = 0.0001, maxdepth = 10)
)

# Plot tree, Save as PNG
png(paste0(dir_name,subdir_name,"/decision_tree.png"), width = 1200, height = 600)
# rpart.plot(tree_model)
# rpart.plot(tree_model, extra = 104, type = 3, under = TRUE, tweak = 1.2)
rpart.plot(tree_model,
           extra = 104,
           type = 3,
           under = TRUE,
           tweak = 1.2,
           box.palette = c("Blues", "Greens"),  # left for Not Useful, right for Useful
           fallen.leaves = TRUE)
dev.off()

