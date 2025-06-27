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
library(tseries)

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
source("~/Dropbox/Treg_problem_local/MISC_FUNCTIONS.R")

args_trailing = commandArgs(trailingOnly = TRUE)
if (length(args_trailing) >= 1) {
  core_index = as.integer(args_trailing[1])
} else {
  stop("Please provide a core_index as the first argument when running the script.")
}
#################### READ ####################################################################################################
all_data_merge   = c()
all_params_merge = c()
dir_name         = './random_simulations'

parameters_df = readRDS(paste0(dir_name,'/parameters_df_',core_index,'.rds'))
#List all files
all_files     = list.files(dir_name, pattern = "\\.rds$", full.names = FALSE)
# Extract the index and case
# Extract index and case from filenames
file_info = strcapture(
  pattern = paste0("dataframe_([0-9]+)_([0-9]+)\\","_",core_index,".rds"),
  x = all_files,
  proto = data.frame(index = integer(), case = integer())
)
# Add filenames back to the data frame
file_info$filename = all_files
# Filter complete sets
complete_sets = file_info %>%
  dplyr::group_by(index) %>%
  dplyr::filter(n() == 4 & all(1:4 %in% case)) %>%
  dplyr::ungroup()

# Read only those files that are part of complete sets
file_paths = file.path(dir_name, complete_sets$filename)
all_data   = do.call(rbind, lapply(file_paths, readRDS))

parameters_df_in_data = parameters_df[sort(unique(all_data$realization_ind)), ]
parameters_df_in_data = parameters_df_in_data %>% mutate(realization_ind = row_number())

all_data$realization_ind              = all_data$realization_ind+1000*core_index # core index to realization index, realizations do not go above 999 per core
parameters_df_in_data$realization_ind = parameters_df_in_data$realization_ind+1000*core_index  # core index to realization index

colnames_insert = c('epithelial_healthy','epithelial_inj_1','epithelial_inj_2','epithelial_inj_3','epithelial_inj_4','epithelial_inj_5',
                    'phagocyte_M0','phagocyte_M1_L_0','phagocyte_M1_L_1','phagocyte_M1_L_2','phagocyte_M1_L_3','phagocyte_M1_L_4','phagocyte_M1_L_5',
                    'phagocyte_M2_L_0','phagocyte_M2_L_1','phagocyte_M2_L_2','phagocyte_M2_L_3','phagocyte_M2_L_4','phagocyte_M2_L_5',
                    'commensal','pathogen','treg_resting','treg_active','C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')

colnames(all_data)[7:37] = colnames_insert

all_data_merge   = rbind(all_data_merge, all_data)
all_params_merge = rbind(all_params_merge, parameters_df_in_data)

saveRDS(all_data_merge,paste0('all_data_merged_',core_index,'.rds'))
saveRDS(all_params_merge,paste0('all_params_merged',core_index,'.rds'))
#################### READ ####################################################################################################
all_data_merge    = readRDS(paste0('all_data_merged_',core_index,'.rds'))
all_params_merged = readRDS(paste0('all_params_merged',core_index,'.rds'))

all_data_epithelial = all_data_merge[c("t","realization_ind","sterile","allow_tregs_to_do_their_job",
                                       "epithelial_healthy","epithelial_inj_1","epithelial_inj_2",
                                       "epithelial_inj_3","epithelial_inj_4","epithelial_inj_5","pathogen")]

tregs_effect_epithelium_sterile_0 = matrix(0, nrow = length(unique(all_data_epithelial$realization_ind)), ncol=4) # realization_ind, better, worse, no diff
tregs_effect_epithelium_sterile_1 = matrix(0, nrow = length(unique(all_data_epithelial$realization_ind)), ncol=4) # realization_ind, better, worse, no diff
tregs_effect_pathogens_sterile_0  = matrix(0, nrow = length(unique(all_data_epithelial$realization_ind)), ncol=4) # realization_ind, better, worse, no diff

row_idx            = 1
grid_size          = 25 
max_level_injury   = 5
mean_level_injury  = mean(c(1,2,3,4,5))
pathogen_threshold = 25
max_t              = max(all_data_merge$t)
check_last         = 1000

for (ind in unique(all_data_epithelial$realization_ind)){
  print(c(core_index,row_idx,length(unique(all_data_epithelial$realization_ind))))
  all_data_epithelial_use    = all_data_epithelial %>% filter(realization_ind==ind)
  all_params_merged_use      = all_params_merged %>% filter(realization_ind==ind)
  
  all_data_epithelial_use = all_data_epithelial_use %>% dplyr::rowwise() %>% dplyr::mutate(injury_score = 5*epithelial_inj_5+
                                                                                             4*epithelial_inj_4+
                                                                                             3*epithelial_inj_3+
                                                                                             2*epithelial_inj_2+
                                                                                             1*epithelial_inj_1)
  
  data_sterile_0_tregs_true  = all_data_epithelial_use %>% filter(sterile==0 & allow_tregs_to_do_their_job==T)
  data_sterile_0_tregs_false = all_data_epithelial_use %>% filter(sterile==0 & allow_tregs_to_do_their_job==F)
  data_sterile_1_tregs_true  = all_data_epithelial_use %>% filter(sterile==1 & allow_tregs_to_do_their_job==T)
  data_sterile_1_tregs_false = all_data_epithelial_use %>% filter(sterile==1 & allow_tregs_to_do_their_job==F)
  
  e_injured_sterile_0_tregs_true  = data_sterile_0_tregs_true$injury_score
  e_injured_sterile_0_tregs_false = data_sterile_0_tregs_false$injury_score
  
  pathogens_sterile_0_tregs_true  = data_sterile_0_tregs_true$pathogen
  pathogens_sterile_0_tregs_false = data_sterile_0_tregs_false$pathogen

  e_injured_sterile_1_tregs_true  = data_sterile_1_tregs_true$injury_score
  e_injured_sterile_1_tregs_false = data_sterile_1_tregs_false$injury_score
  
  pathogens_sterile_1_tregs_true  = data_sterile_1_tregs_true$pathogen
  pathogens_sterile_1_tregs_false = data_sterile_1_tregs_false$pathogen
  
  ######### STERILE 0 EPITHELIUM ##################  
  # determine steady time through epithelial cells
  # t_sterile_0_tregs_true  = get_steady_time(e_injured_sterile_0_tregs_true, 500, 100, 0.5)
  # t_sterile_0_tregs_false = get_steady_time(e_injured_sterile_0_tregs_false, 500, 100, 0.5)
  
  t_sterile_0_tregs_true  = max_t-check_last+1
  t_sterile_0_tregs_false = max_t-check_last+1
  
  steady_e_injured_sterile_0_tregs_true  = e_injured_sterile_0_tregs_true[t_sterile_0_tregs_true:max_t]
  steady_e_injured_sterile_0_tregs_false = e_injured_sterile_0_tregs_false[t_sterile_0_tregs_false:max_t]
  
  # # steady_e_injured_sterile_0_tregs_true = zoo::rollmean(steady_e_injured_sterile_0_tregs_true, k = 100, fill = NA, align = "center")
  # # steady_e_injured_sterile_0_tregs_false = zoo::rollmean(steady_e_injured_sterile_0_tregs_false, k = 100, fill = NA, align = "center")
  
  # df_t = data.frame(
  #   value = c(steady_e_injured_sterile_0_tregs_true, steady_e_injured_sterile_0_tregs_false),
  #   group = c(
  #     rep("tregs_true", length(steady_e_injured_sterile_0_tregs_true)),
  #     rep("tregs_false", length(steady_e_injured_sterile_0_tregs_false))
  #   )
  # )
  # 
  # ggplot(df_t, aes(x = value, fill = group)) +
  #   geom_density(alpha = 0.5) +
  #   labs(title = "Density Plot", x = "Value") +
  #   theme_minimal()
   
  # Compute Cohen's d
  m1 = mean(steady_e_injured_sterile_0_tregs_true)
  m2 = mean(steady_e_injured_sterile_0_tregs_false)
  s1 = sd(steady_e_injured_sterile_0_tregs_true)
  s2 = sd(steady_e_injured_sterile_0_tregs_false)
  n1 = length(steady_e_injured_sterile_0_tregs_true)
  n2 = length(steady_e_injured_sterile_0_tregs_false)
  s_pooled = sqrt(((n1 - 1)*s1^2 + (n2 - 1)*s2^2) / (n1 + n2 - 2))
  d = (m1 - m2) / s_pooled
  
  # Interpret effect size
  if (!is.na(d) & d < -0.5 & abs((m1 - m2))>mean_level_injury*grid_size*0.25) {
    tregs_effect_epithelium_sterile_0[row_idx,] = c(ind,1,0,0)
  } else if (!is.na(d) & d > +0.5 & abs((m1 - m2))>mean_level_injury*grid_size*0.25) {
    tregs_effect_epithelium_sterile_0[row_idx,] = c(ind,0,1,0)
  } else {
    tregs_effect_epithelium_sterile_0[row_idx,] = c(ind,0,0,1)
  }
  ######### STERILE 0 EPITHELIUM ##################  
  
  ######### STERILE 1 EPITHELIUM ##################  
  # determine steady time through epithelial cells
  # t_sterile_1_tregs_true  = get_steady_time(e_injured_sterile_1_tregs_true, 500, 100, 0.5)
  # t_sterile_1_tregs_false = get_steady_time(e_injured_sterile_1_tregs_false, 500, 100, 0.5)
  
  t_sterile_1_tregs_true  = max_t-check_last+1
  t_sterile_1_tregs_false = max_t-check_last+1
  
  steady_e_injured_sterile_1_tregs_true  = e_injured_sterile_1_tregs_true[t_sterile_1_tregs_true:max_t]
  steady_e_injured_sterile_1_tregs_false = e_injured_sterile_1_tregs_false[t_sterile_1_tregs_false:max_t]
  
  m1 = mean(steady_e_injured_sterile_1_tregs_true)
  m2 = mean(steady_e_injured_sterile_1_tregs_false)
  s1 = sd(steady_e_injured_sterile_1_tregs_true)
  s2 = sd(steady_e_injured_sterile_1_tregs_false)
  n1 = length(steady_e_injured_sterile_1_tregs_true)
  n2 = length(steady_e_injured_sterile_1_tregs_false)
  s_pooled = sqrt(((n1 - 1)*s1^2 + (n2 - 1)*s2^2) / (n1 + n2 - 2))
  d = (m1 - m2) / s_pooled
  
  if (!is.na(d) & d < -0.5 & abs((m1 - m2))>mean_level_injury*grid_size*0.25) {
    tregs_effect_epithelium_sterile_1[row_idx,] = c(ind,1,0,0)
  } else if (!is.na(d) & d > +0.5 & abs((m1 - m2))>mean_level_injury*grid_size*0.25) {
    tregs_effect_epithelium_sterile_1[row_idx,] = c(ind,0,1,0)
  } else {
    tregs_effect_epithelium_sterile_1[row_idx,] = c(ind,0,0,1)
  }
  ######### STERILE 1 EPITHELIUM ##################  
  
  ######### STERILE 0 PATHOGEN ##################  
  # use steady time determined through epithelial cells
  steady_pathogens_sterile_0_tregs_true  = pathogens_sterile_0_tregs_true[t_sterile_0_tregs_true:max_t]
  steady_pathogens_sterile_0_tregs_false = pathogens_sterile_0_tregs_false[t_sterile_0_tregs_false:max_t]
  
  # Compute Cohen's d (negative d means tregs helped reduce pathogens)
  m1 = mean(steady_pathogens_sterile_0_tregs_true)
  m2 = mean(steady_pathogens_sterile_0_tregs_false)
  s1 = sd(steady_pathogens_sterile_0_tregs_true)
  s2 = sd(steady_pathogens_sterile_0_tregs_false)
  n1 = length(steady_pathogens_sterile_0_tregs_true)
  n2 = length(steady_pathogens_sterile_0_tregs_false)
  s_pooled = sqrt(((n1 - 1)*s1^2 + (n2 - 1)*s2^2) / (n1 + n2 - 2))
  d = (m1 - m2) / s_pooled
  
  
  df_t = data.frame(
    value = c(steady_pathogens_sterile_0_tregs_true, steady_pathogens_sterile_0_tregs_false),
    group = c(
      rep("tregs_true", length(steady_pathogens_sterile_0_tregs_true)),
      rep("tregs_false", length(steady_pathogens_sterile_0_tregs_false))
    )
  )

  ggplot(df_t, aes(x = value, fill = group)) +
    geom_density(alpha = 0.5) +
    labs(title = "Density Plot", x = "Value") +
    theme_minimal()
  
  
  # Interpret direction: LOWER pathogen = BETTER = d < -0.5
  if (!is.na(d) & d < -0.5 & m1<m2) {
    tregs_effect_pathogens_sterile_0[row_idx,] = c(ind,1,0,0)  # tregs helped reduce pathogens
  } else if (!is.na(d) & d > 0.5 & m1>m2) {
    tregs_effect_pathogens_sterile_0[row_idx,] = c(ind,0,1,0)  # tregs worsened pathogen burden
  } else {
    tregs_effect_pathogens_sterile_0[row_idx,] = c(ind,0,0,1)  # no strong effect
  }
  ######### STERILE 0 PATHOGEN ##################  
  
  ######### INCREASE ROW INDEX ##################  
  row_idx = row_idx + 1
  ######### INCREASE ROW INDEX ##################  
  
}

############# STERILE 0
####### TREGS *USEFUL* FOR EPITHELIUM HEALTH *AND CAN* ELIMINATE THE PATHOGEN, STERILE 0
vec_1 = tregs_effect_epithelium_sterile_0[which(tregs_effect_epithelium_sterile_0[,2]==1),1]
vec_2 = tregs_effect_pathogens_sterile_0[which(tregs_effect_pathogens_sterile_0[,2]==1),1]
tregs_useful_e_yes_p_yes_sterile_0 = sort(intersect(vec_1, vec_2))

####### TREGS *USEFUL* FOR EPITHELIUM HEALTH *BUT CAN'T* ELIMINATE THE PATHOGEN, STERILE 0
vec_1 = tregs_effect_epithelium_sterile_0[which(tregs_effect_epithelium_sterile_0[,2]==1),1]
vec_2 = tregs_effect_pathogens_sterile_0[which(tregs_effect_pathogens_sterile_0[,3]==1),1]
tregs_useful_e_yes_p_no_sterile_0 = sort(intersect(vec_1, vec_2))

####### TREGS *HARMFUL* FOR EPITHELIUM HEALTH *BUT CAN* ELIMINATE THE PATHOGEN, STERILE 0
vec_1 = tregs_effect_epithelium_sterile_0[which(tregs_effect_epithelium_sterile_0[,3]==1),1]
vec_2 = tregs_effect_pathogens_sterile_0[which(tregs_effect_pathogens_sterile_0[,2]==1),1]
tregs_useful_e_no_p_yes_sterile_0 = sort(intersect(vec_1, vec_2))

####### TREGS *HARMFUL* FOR EPITHELIUM HEALTH *AND CAN'T* ELIMINATE THE PATHOGEN, STERILE 0
vec_1 = tregs_effect_epithelium_sterile_0[which(tregs_effect_epithelium_sterile_0[,3]==1),1]
vec_2 = tregs_effect_pathogens_sterile_0[which(tregs_effect_pathogens_sterile_0[,3]==1),1]
tregs_useful_e_no_p_no_sterile_0 = sort(intersect(vec_1, vec_2))

####### TREGS *DON'T MATTER* FOR EPITHELIUM HEALTH *AND CAN'T* ELIMINATE THE PATHOGEN, STERILE 0
vec_1 = tregs_effect_epithelium_sterile_0[which(tregs_effect_epithelium_sterile_0[,4]==1),1]
vec_2 = tregs_effect_pathogens_sterile_0[which(tregs_effect_pathogens_sterile_0[,3]==1),1]
tregs_useful_e_dm_p_no_sterile_0 = sort(intersect(vec_1, vec_2))

####### TREGS *DON'T MATTER* FOR EPITHELIUM HEALTH *BUT CAN* ELIMINATE THE PATHOGEN, STERILE 0
vec_1 = tregs_effect_epithelium_sterile_0[which(tregs_effect_epithelium_sterile_0[,4]==1),1]
vec_2 = tregs_effect_pathogens_sterile_0[which(tregs_effect_pathogens_sterile_0[,2]==1),1]
tregs_useful_e_dm_p_yes_sterile_0 = sort(intersect(vec_1, vec_2))

####### TREGS *USEFUL* FOR EPITHELIUM HEALTH *BUT DON'T MATTER* ELIMINATE THE PATHOGEN, STERILE 0
vec_1 = tregs_effect_epithelium_sterile_0[which(tregs_effect_epithelium_sterile_0[,2]==1),1]
vec_2 = tregs_effect_pathogens_sterile_0[which(tregs_effect_pathogens_sterile_0[,4]==1),1]
tregs_useful_e_yes_p_dm_sterile_0 = sort(intersect(vec_1, vec_2))

####### TREGS *HARMFUL* FOR EPITHELIUM HEALTH *BUT DON'T MATTER* ELIMINATE THE PATHOGEN, STERILE 0
vec_1 = tregs_effect_epithelium_sterile_0[which(tregs_effect_epithelium_sterile_0[,3]==1),1]
vec_2 = tregs_effect_pathogens_sterile_0[which(tregs_effect_pathogens_sterile_0[,4]==1),1]
tregs_useful_e_no_p_dm_sterile_0 = sort(intersect(vec_1, vec_2))

####### TREGS *DON'T MATTER* FOR EPITHELIUM HEALTH *AND DON'T MATTER* ELIMINATE THE PATHOGEN, STERILE 0
vec_1 = tregs_effect_epithelium_sterile_0[which(tregs_effect_epithelium_sterile_0[,4]==1),1]
vec_2 = tregs_effect_pathogens_sterile_0[which(tregs_effect_pathogens_sterile_0[,4]==1),1]
tregs_useful_e_dm_p_dm_sterile_0 = sort(intersect(vec_1, vec_2))

tregs_useful_e_yes_p_yes_sterile_0
saveRDS(tregs_useful_e_yes_p_yes_sterile_0, paste0('tregs_useful_e_yes_p_yes_sterile_0_',core_index,'.rds'))

tregs_useful_e_yes_p_no_sterile_0
saveRDS(tregs_useful_e_yes_p_no_sterile_0, paste0('tregs_useful_e_yes_p_no_sterile_0_',core_index,'.rds'))

tregs_useful_e_yes_p_dm_sterile_0
saveRDS(tregs_useful_e_yes_p_dm_sterile_0, paste0('tregs_useful_e_yes_p_dm_sterile_0_',core_index,'.rds'))

tregs_useful_e_no_p_yes_sterile_0
saveRDS(tregs_useful_e_no_p_yes_sterile_0, paste0('tregs_useful_e_no_p_yes_sterile_0_',core_index,'.rds'))

tregs_useful_e_no_p_no_sterile_0
saveRDS(tregs_useful_e_no_p_no_sterile_0, paste0('tregs_useful_e_no_p_no_sterile_0_',core_index,'.rds'))

tregs_useful_e_no_p_dm_sterile_0
saveRDS(tregs_useful_e_no_p_dm_sterile_0, paste0('tregs_useful_e_no_p_dm_sterile_0_',core_index,'.rds'))

tregs_useful_e_dm_p_yes_sterile_0
saveRDS(tregs_useful_e_dm_p_yes_sterile_0, paste0('tregs_useful_e_dm_p_yes_sterile_0_',core_index,'.rds'))

tregs_useful_e_dm_p_no_sterile_0
saveRDS(tregs_useful_e_dm_p_no_sterile_0, paste0('tregs_useful_e_dm_p_no_sterile_0_',core_index,'.rds'))

tregs_useful_e_dm_p_dm_sterile_0
saveRDS(tregs_useful_e_dm_p_dm_sterile_0, paste0('tregs_useful_e_dm_p_dm_sterile_0_',core_index,'.rds'))


############# STERILE 1
####### TREGS *USEFUL* FOR EPITHELIUM HEALTH, STERILE 1
vec_1 = tregs_effect_epithelium_sterile_1[which(tregs_effect_epithelium_sterile_1[,2]==1),1]
tregs_useful_e_yes_p_dm_sterile_1 = sort(vec_1)

####### TREGS *HARMFUL* FOR EPITHELIUM HEALTH *BUT DON'T MATTER* ELIMINATE THE PATHOGEN, STERILE 1
vec_1 = tregs_effect_epithelium_sterile_1[which(tregs_effect_epithelium_sterile_1[,3]==1),1]
tregs_useful_e_no_p_dm_sterile_1 = sort(vec_1)

####### TREGS *DON'T MATTER* FOR EPITHELIUM HEALTH *AND DON'T MATTER* ELIMINATE THE PATHOGEN, STERILE 1
vec_1 = tregs_effect_epithelium_sterile_1[which(tregs_effect_epithelium_sterile_1[,4]==1),1]
tregs_useful_e_dm_p_dm_sterile_1 = sort(vec_1)

tregs_useful_e_yes_p_dm_sterile_1
saveRDS(tregs_useful_e_yes_p_dm_sterile_1, paste0('tregs_useful_e_yes_p_dm_sterile_1_',core_index,'.rds'))

tregs_useful_e_no_p_dm_sterile_1
saveRDS(tregs_useful_e_no_p_dm_sterile_1, paste0('tregs_useful_e_no_p_dm_sterile_1_',core_index,'.rds'))

tregs_useful_e_dm_p_dm_sterile_1
saveRDS(tregs_useful_e_dm_p_dm_sterile_1, paste0('tregs_useful_e_dm_p_dm_sterile_1_',core_index,'.rds'))






