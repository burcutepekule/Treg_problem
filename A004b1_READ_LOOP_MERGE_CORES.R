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

all_data_merge   = c()
all_params_merge = c()
for (prefix  in seq(1,4)){
  print(prefix)
  keep_objects = c("all_data_merge","all_params_merge","prefix")
  rm(list = setdiff(ls(), keep_objects))
  
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
  
  # source("./PLOT_FUNCTIONS.R")
  dir_name      = './random_simulations'
  # subdir_name   = '/figures'
  # dir.create(paste0(dir_name,subdir_name))
  
  parameters_df = readRDS(paste0(dir_name,'/parameters_df_',prefix,'.rds'))
  #List all files
  all_files     = list.files(dir_name, pattern = "\\.rds$", full.names = FALSE)
  # Extract the index and case
  # Extract index and case from filenames
  file_info = strcapture(
    pattern = paste0("dataframe_([0-9]+)_([0-9]+)\\","_",prefix,".rds"),
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
  
  all_data$realization_ind              = all_data$realization_ind+1000*prefix # core index to realization index, realizations do not go above 999 per core
  parameters_df_in_data$realization_ind = parameters_df_in_data$realization_ind+1000*prefix  # core index to realization index
  

  colnames_insert = c('epithelial_healthy','epithelial_inj_1','epithelial_inj_2','epithelial_inj_3','epithelial_inj_4','epithelial_inj_5',
                      'phagocyte_M0','phagocyte_M1_L_0','phagocyte_M1_L_1','phagocyte_M1_L_2','phagocyte_M1_L_3','phagocyte_M1_L_4','phagocyte_M1_L_5',
                      'phagocyte_M2_L_0','phagocyte_M2_L_1','phagocyte_M2_L_2','phagocyte_M2_L_3','phagocyte_M2_L_4','phagocyte_M2_L_5',
                      'commensal','pathogen','treg_resting','treg_active','C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')
  
  colnames(all_data)[7:37] = colnames_insert

  all_data_merge   = rbind(all_data_merge, all_data)
  all_params_merge = rbind(all_params_merge, parameters_df_in_data)
}

dim(all_data_merge)
dim(all_params_merge)

saveRDS(all_data_merge,'all_data_merged.rds')
saveRDS(all_params_merge,'all_params_merged.rds')
