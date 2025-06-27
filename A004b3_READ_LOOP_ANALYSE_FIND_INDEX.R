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

source('~/Dropbox/Treg_problem_local/A004b1_READ_LOOP_MERGE_CORES.R')

all_data   = readRDS('all_data_merged.rds')
all_params = readRDS('all_params_merged.rds')

all_data_epithelial = all_data[c("t","realization_ind","sterile","allow_tregs_to_do_their_job",
                                 "epithelial_healthy","epithelial_inj_1","epithelial_inj_2",
                                 "epithelial_inj_3","epithelial_inj_4","epithelial_inj_5","pathogen")]

# #### SUBSAMPLE FOR POC
# inds_pick = c(1301, 1504, 2037, 2286, 2418, 3054, 3076, 3077, 3329, 3671, 4118, 4175, sample(4000: max(all_data_epithelial$realization_ind),8))
# all_data_epithelial = all_data_epithelial %>% filter(realization_ind %in% inds_pick)

max_t = 5000
tregs_effect_epithelium_sterile_0 = matrix(0, nrow = length(unique(all_data_epithelial$realization_ind)), ncol=4) # realization_ind, better, worse, no diff
tregs_effect_epithelium_sterile_1 = matrix(0, nrow = length(unique(all_data_epithelial$realization_ind)), ncol=4) # realization_ind, better, worse, no diff
tregs_effect_pathogens_sterile_0  = matrix(0, nrow = length(unique(all_data_epithelial$realization_ind)), ncol=4) # realization_ind, better, worse, no diff

row_idx = 1
for (ind in unique(all_data_epithelial$realization_ind)){
  print(c(row_idx,length(unique(all_data_epithelial$realization_ind))))
  all_data_epithelial_use    = all_data_epithelial %>% filter(realization_ind==ind)
  
  data_sterile_0_tregs_true  = all_data_epithelial_use %>% filter(sterile==0 & allow_tregs_to_do_their_job==T)
  data_sterile_0_tregs_false = all_data_epithelial_use %>% filter(sterile==0 & allow_tregs_to_do_their_job==F)
  data_sterile_1_tregs_true  = all_data_epithelial_use %>% filter(sterile==1 & allow_tregs_to_do_their_job==T)
  data_sterile_1_tregs_false = all_data_epithelial_use %>% filter(sterile==1 & allow_tregs_to_do_their_job==F)

  e_healthy_sterile_0_tregs_true  = data_sterile_0_tregs_true$epithelial_healthy
  e_healthy_sterile_0_tregs_false = data_sterile_0_tregs_false$epithelial_healthy
  pathogens_sterile_0_tregs_true  = data_sterile_0_tregs_true$pathogen
  pathogens_sterile_0_tregs_false = data_sterile_0_tregs_false$pathogen
  
  e_healthy_sterile_1_tregs_true  = data_sterile_1_tregs_true$epithelial_healthy
  e_healthy_sterile_1_tregs_false = data_sterile_1_tregs_false$epithelial_healthy
  pathogens_sterile_1_tregs_true  = data_sterile_1_tregs_true$pathogen
  pathogens_sterile_1_tregs_false = data_sterile_1_tregs_false$pathogen
  
  ######### STERILE 0 EPITHELIUM ##################  
  # determine steady time through epithelial cells
  t_sterile_0_tregs_true  = get_steady_time(e_healthy_sterile_0_tregs_true, 250, 100, 0.5)
  t_sterile_0_tregs_false = get_steady_time(e_healthy_sterile_0_tregs_false, 250, 100, 0.5)
  
  steady_e_healthy_sterile_0_tregs_true  = e_healthy_sterile_0_tregs_true[t_sterile_0_tregs_true:max_t]
  steady_e_healthy_sterile_0_tregs_false = e_healthy_sterile_0_tregs_false[t_sterile_0_tregs_false:max_t]
  
  p_better = t.test(steady_e_healthy_sterile_0_tregs_true, steady_e_healthy_sterile_0_tregs_false, alternative = "greater")$p.value
  p_worse  = t.test(steady_e_healthy_sterile_0_tregs_true, steady_e_healthy_sterile_0_tregs_false, alternative = "less")$p.value
  
  # Interpret result
  if (!is.na(p_better) & p_better < 0.05) {
    tregs_effect_epithelium_sterile_0[row_idx,] = c(ind,1,0,0)
  } else if (!is.na(p_worse) & p_worse < 0.05) {
    tregs_effect_epithelium_sterile_0[row_idx,] = c(ind,0,1,0)
  } else {
    tregs_effect_epithelium_sterile_0[row_idx,] = c(ind,0,0,1)
  }
  ######### STERILE 0 EPITHELIUM ##################  
  
  ######### STERILE 1 EPITHELIUM ##################  
  # determine steady time through epithelial cells
  t_sterile_1_tregs_true  = get_steady_time(e_healthy_sterile_1_tregs_true, 250, 100, 0.5)
  t_sterile_1_tregs_false = get_steady_time(e_healthy_sterile_1_tregs_false, 250, 100, 0.5)
  
  steady_e_healthy_sterile_1_tregs_true  = e_healthy_sterile_1_tregs_true[t_sterile_1_tregs_true:max_t]
  steady_e_healthy_sterile_1_tregs_false = e_healthy_sterile_1_tregs_false[t_sterile_1_tregs_false:max_t]
  
  p_better = t.test(steady_e_healthy_sterile_1_tregs_true, steady_e_healthy_sterile_1_tregs_false, alternative = "greater")$p.value
  p_worse  = t.test(steady_e_healthy_sterile_1_tregs_true, steady_e_healthy_sterile_1_tregs_false, alternative = "less")$p.value
  
  # Interpret result
  if (!is.na(p_better) & p_better < 0.05) {
    tregs_effect_epithelium_sterile_1[row_idx,] = c(ind,1,0,0)
  } else if (!is.na(p_worse) & p_worse < 0.05) {
    tregs_effect_epithelium_sterile_1[row_idx,] = c(ind,0,1,0)
  } else {
    tregs_effect_epithelium_sterile_1[row_idx,] = c(ind,0,0,1)
  }
  ######### STERILE 1 EPITHELIUM ##################  

  ######### STERILE 0 PATHOGEN ##################  
  # use steady time determined through epithelial cells
  steady_pathogens_sterile_0_tregs_true  = pathogens_sterile_0_tregs_true[t_sterile_0_tregs_true:max_t]
  steady_pathogens_sterile_0_tregs_false = pathogens_sterile_0_tregs_false[t_sterile_0_tregs_false:max_t]
  
  p_better = t.test(steady_pathogens_sterile_0_tregs_true, steady_pathogens_sterile_0_tregs_false, alternative = "less")$p.value
  p_worse  = t.test(steady_pathogens_sterile_0_tregs_true, steady_pathogens_sterile_0_tregs_false, alternative = "greater")$p.value
  
  # Interpret result
  if (!is.na(p_better) & p_better < 0.05) {
    tregs_effect_pathogens_sterile_0[row_idx,] = c(ind,1,0,0)
  } else if (!is.na(p_worse) & p_worse < 0.05) {
    tregs_effect_pathogens_sterile_0[row_idx,] = c(ind,0,1,0)
  } else {
    tregs_effect_pathogens_sterile_0[row_idx,] = c(ind,0,0,1)
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
saveRDS(tregs_useful_e_yes_p_yes_sterile_0, 'tregs_useful_e_yes_p_yes_sterile_0.rds')

tregs_useful_e_yes_p_no_sterile_0
saveRDS(tregs_useful_e_yes_p_no_sterile_0, 'tregs_useful_e_yes_p_no_sterile_0.rds')

tregs_useful_e_yes_p_dm_sterile_0
saveRDS(tregs_useful_e_yes_p_dm_sterile_0, 'tregs_useful_e_yes_p_dm_sterile_0.rds')

tregs_useful_e_no_p_yes_sterile_0
saveRDS(tregs_useful_e_no_p_yes_sterile_0, 'tregs_useful_e_no_p_yes_sterile_0.rds')

tregs_useful_e_no_p_no_sterile_0
saveRDS(tregs_useful_e_no_p_no_sterile_0, 'tregs_useful_e_no_p_no_sterile_0.rds')

tregs_useful_e_no_p_dm_sterile_0
saveRDS(tregs_useful_e_no_p_dm_sterile_0, 'tregs_useful_e_no_p_dm_sterile_0.rds')

tregs_useful_e_dm_p_yes_sterile_0
saveRDS(tregs_useful_e_dm_p_yes_sterile_0, 'tregs_useful_e_dm_p_yes_sterile_0.rds')

tregs_useful_e_dm_p_no_sterile_0
saveRDS(tregs_useful_e_dm_p_no_sterile_0, 'tregs_useful_e_dm_p_no_sterile_0.rds')

tregs_useful_e_dm_p_dm_sterile_0
saveRDS(tregs_useful_e_dm_p_dm_sterile_0, 'tregs_useful_e_dm_p_dm_sterile_0.rds')


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
saveRDS(tregs_useful_e_yes_p_dm_sterile_1, 'tregs_useful_e_yes_p_dm_sterile_1.rds')

tregs_useful_e_no_p_dm_sterile_1
saveRDS(tregs_useful_e_no_p_dm_sterile_1, 'tregs_useful_e_no_p_dm_sterile_1.rds')

tregs_useful_e_dm_p_dm_sterile_1
saveRDS(tregs_useful_e_dm_p_dm_sterile_1, 'tregs_useful_e_dm_p_dm_sterile_1.rds')






