rm(list=ls())
library(dplyr)
library(tidyr)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)

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

source("./BASE_FUNCTIONS.R")
source("./PLOT_FUNCTIONS.R")
dir_name      = './random_simulations'
subdir_name   = '/figures'
dir.create(paste0(dir_name,subdir_name))

parameters_df = readRDS(paste0(dir_name,'/parameters_df.rds'))
#List all files
all_files     = list.files(dir_name, pattern = "\\.rds$", full.names = FALSE)
# Extract the index and case
# Extract index and case from filenames
file_info = strcapture(
  pattern = "dataframe_([0-9]+)_([0-9]+)\\.rds",
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
colnames(all_data)

parameters_df_in_data = parameters_df[sort(unique(all_data$realization_ind)), ]
parameters_df_in_data = parameters_df_in_data %>% mutate(realization_ind = row_number())
rm(parameters_df) # remove, too big?

colnames_insert = c('epithelial_healthy','epithelial_inj_1','epithelial_inj_2','epithelial_inj_3','epithelial_inj_4','epithelial_inj_5',
                    'phagocyte_M0','phagocyte_M1_L_0','phagocyte_M1_L_1','phagocyte_M1_L_2','phagocyte_M1_L_3','phagocyte_M1_L_4','phagocyte_M1_L_5',
                    'phagocyte_M2_L_0','phagocyte_M2_L_1','phagocyte_M2_L_2','phagocyte_M2_L_3','phagocyte_M2_L_4','phagocyte_M2_L_5',
                    'commensal','pathogen','treg_resting','treg_active','C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')

colnames(all_data)[6:36] = colnames_insert

# all_data_complete = merge(all_data, parameters_df_in_data, by='realization_ind')

all_figure_files = list.files(paste0(dir_name,subdir_name), pattern = "\\.png$", full.names = FALSE)
file_figure_info = strcapture(
  pattern = ".*_([0-9]+)\\.png",
  x = all_figure_files,
  proto = data.frame(index = integer())
)

max_fig_index = max(file_figure_info)
if(is.infinite(max_fig_index)){
  max_fig_index = 0
}

for (ind in (1+max_fig_index):max(sort(unique(all_data$realization_ind)))){
  print(ind)
  realization_data = all_data %>% filter(realization_ind==ind)
  
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
  #                             "M2 Phagocytes & Active Tregs")
  # p_microbes   = plot_faceted(realization_data, c("commensal", "pathogen"),
  #                             "Microbe Dynamics")
  # p_cum_com    = plot_faceted(realization_data, 
  #                             c("C_ROS", "C_M0", "C_M1", "C_M2"),
  #                             "Cumulative Commensal Deaths")
  # p_cum_path    = plot_faceted(realization_data, 
  #                              c("P_ROS", "P_M0", "P_M1", "P_M2"),
  #                              "Cumulative Pathogen Deaths")
  
  # Save all 7 plots with realization_ind in filename
  ggsave(paste0(dir_name,subdir_name,"/epithelium_", ind, ".png"), p_epithelium, dpi = 300, width = 8, height = 5)
  # ggsave(paste0(dir_name,subdir_name,"/resting_tregs_", ind, ".png"), p_resting, dpi = 300, width = 8, height = 5)
  # ggsave(paste0(dir_name,subdir_name,"/phagocyte_M1_", ind, ".png"), p_M1, dpi = 300, width = 8, height = 5)
  # ggsave(paste0(dir_name,subdir_name,"/phagocyte_M2_activeTreg_", ind, ".png"), p_reg_active, dpi = 300, width = 8, height = 5)
  # ggsave(paste0(dir_name,subdir_name,"/microbes_", ind, ".png"), p_microbes, dpi = 300, width = 8, height = 5)
  # ggsave(paste0(dir_name,subdir_name,"/cumulative_commensal_deaths_", ind, ".png"), p_cum_com, dpi = 300, width = 8, height = 5)
  # ggsave(paste0(dir_name,subdir_name,"/cumulative_pathogen_deaths_", ind, ".png"), p_cum_path, dpi = 300, width = 8, height = 5)
}
