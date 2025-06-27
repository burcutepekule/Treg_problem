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
library(MASS)


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

# Initialize merged storage
all_data_merged_keep   = c()
all_params_merged_keep = c()

# Initialize Tregs index collections
tregs_useful_e_yes_p_yes_sterile_0_keep = c()
tregs_useful_e_yes_p_no_sterile_0_keep  = c()
tregs_useful_e_yes_p_dm_sterile_0_keep = c()
tregs_useful_e_no_p_yes_sterile_0_keep = c()
tregs_useful_e_no_p_no_sterile_0_keep  = c()
tregs_useful_e_no_p_dm_sterile_0_keep  = c()
tregs_useful_e_dm_p_yes_sterile_0_keep = c()
tregs_useful_e_dm_p_no_sterile_0_keep  = c()
tregs_useful_e_dm_p_dm_sterile_0_keep  = c()

tregs_useful_e_yes_p_dm_sterile_1_keep = c()
tregs_useful_e_no_p_dm_sterile_1_keep  = c()
tregs_useful_e_dm_p_dm_sterile_1_keep  = c()

# Loop to read and collect
for (core_index in seq(1, 9)) {
  print(core_index)
  # Read merged simulation data and parameters
  all_data_merged   = readRDS(paste0('all_data_merged_', core_index, '.rds'))
  all_params_merged = readRDS(paste0('all_params_merged', core_index, '.rds'))
  
  # Read Tregs classifications - sterile 0
  tregs_useful_e_yes_p_yes_sterile_0 = readRDS(paste0('tregs_useful_e_yes_p_yes_sterile_0_', core_index, '.rds'))
  tregs_useful_e_yes_p_no_sterile_0  = readRDS(paste0('tregs_useful_e_yes_p_no_sterile_0_', core_index, '.rds'))
  tregs_useful_e_yes_p_dm_sterile_0  = readRDS(paste0('tregs_useful_e_yes_p_dm_sterile_0_', core_index, '.rds'))
  tregs_useful_e_no_p_yes_sterile_0  = readRDS(paste0('tregs_useful_e_no_p_yes_sterile_0_', core_index, '.rds'))
  tregs_useful_e_no_p_no_sterile_0   = readRDS(paste0('tregs_useful_e_no_p_no_sterile_0_', core_index, '.rds'))
  tregs_useful_e_no_p_dm_sterile_0   = readRDS(paste0('tregs_useful_e_no_p_dm_sterile_0_', core_index, '.rds'))
  tregs_useful_e_dm_p_yes_sterile_0  = readRDS(paste0('tregs_useful_e_dm_p_yes_sterile_0_', core_index, '.rds'))
  tregs_useful_e_dm_p_no_sterile_0   = readRDS(paste0('tregs_useful_e_dm_p_no_sterile_0_', core_index, '.rds'))
  tregs_useful_e_dm_p_dm_sterile_0   = readRDS(paste0('tregs_useful_e_dm_p_dm_sterile_0_', core_index, '.rds'))
  
  # Read Tregs classifications - sterile 1
  tregs_useful_e_yes_p_dm_sterile_1 = readRDS(paste0('tregs_useful_e_yes_p_dm_sterile_1_', core_index, '.rds'))
  tregs_useful_e_no_p_dm_sterile_1  = readRDS(paste0('tregs_useful_e_no_p_dm_sterile_1_', core_index, '.rds'))
  tregs_useful_e_dm_p_dm_sterile_1  = readRDS(paste0('tregs_useful_e_dm_p_dm_sterile_1_', core_index, '.rds'))
  
  # Merge main data
  all_data_merged_keep   = rbind(all_data_merged_keep, all_data_merged)
  all_params_merged_keep = rbind(all_params_merged_keep, all_params_merged)
  
  # Merge Tregs indices
  tregs_useful_e_yes_p_yes_sterile_0_keep = c(tregs_useful_e_yes_p_yes_sterile_0_keep, tregs_useful_e_yes_p_yes_sterile_0)
  tregs_useful_e_yes_p_no_sterile_0_keep  = c(tregs_useful_e_yes_p_no_sterile_0_keep,  tregs_useful_e_yes_p_no_sterile_0)
  tregs_useful_e_yes_p_dm_sterile_0_keep  = c(tregs_useful_e_yes_p_dm_sterile_0_keep,  tregs_useful_e_yes_p_dm_sterile_0)
  tregs_useful_e_no_p_yes_sterile_0_keep  = c(tregs_useful_e_no_p_yes_sterile_0_keep,  tregs_useful_e_no_p_yes_sterile_0)
  tregs_useful_e_no_p_no_sterile_0_keep   = c(tregs_useful_e_no_p_no_sterile_0_keep,   tregs_useful_e_no_p_no_sterile_0)
  tregs_useful_e_no_p_dm_sterile_0_keep   = c(tregs_useful_e_no_p_dm_sterile_0_keep,   tregs_useful_e_no_p_dm_sterile_0)
  tregs_useful_e_dm_p_yes_sterile_0_keep  = c(tregs_useful_e_dm_p_yes_sterile_0_keep,  tregs_useful_e_dm_p_yes_sterile_0)
  tregs_useful_e_dm_p_no_sterile_0_keep   = c(tregs_useful_e_dm_p_no_sterile_0_keep,   tregs_useful_e_dm_p_no_sterile_0)
  tregs_useful_e_dm_p_dm_sterile_0_keep   = c(tregs_useful_e_dm_p_dm_sterile_0_keep,   tregs_useful_e_dm_p_dm_sterile_0)
  
  tregs_useful_e_yes_p_dm_sterile_1_keep = c(tregs_useful_e_yes_p_dm_sterile_1_keep, tregs_useful_e_yes_p_dm_sterile_1)
  tregs_useful_e_no_p_dm_sterile_1_keep  = c(tregs_useful_e_no_p_dm_sterile_1_keep,  tregs_useful_e_no_p_dm_sterile_1)
  tregs_useful_e_dm_p_dm_sterile_1_keep  = c(tregs_useful_e_dm_p_dm_sterile_1_keep,  tregs_useful_e_dm_p_dm_sterile_1)
}

total_realizations = length(unique(all_data_merged_keep$realization_ind))

# e_yes, sterile 0
p_e_yes_p_yes_sterile_0 = 100 * length(tregs_useful_e_yes_p_yes_sterile_0_keep) / total_realizations  # 2.32
p_e_yes_p_no_sterile_0  = 100 * length(tregs_useful_e_yes_p_no_sterile_0_keep)  / total_realizations  # 0
p_e_yes_p_dm_sterile_0  = 100 * length(tregs_useful_e_yes_p_dm_sterile_0_keep)  / total_realizations  # 9.30
print(c(p_e_yes_p_yes_sterile_0, p_e_yes_p_no_sterile_0, p_e_yes_p_dm_sterile_0, p_e_yes_p_yes_sterile_0 + p_e_yes_p_no_sterile_0 + p_e_yes_p_dm_sterile_0))

# e_no, sterile 0
p_e_no_p_yes_sterile_0 = 100 * length(tregs_useful_e_no_p_yes_sterile_0_keep) / total_realizations  # 0
p_e_no_p_no_sterile_0  = 100 * length(tregs_useful_e_no_p_no_sterile_0_keep)  / total_realizations  # 0
p_e_no_p_dm_sterile_0  = 100 * length(tregs_useful_e_no_p_dm_sterile_0_keep)  / total_realizations  # 0
print(c(p_e_no_p_yes_sterile_0, p_e_no_p_no_sterile_0, p_e_no_p_dm_sterile_0, p_e_no_p_yes_sterile_0 + p_e_no_p_no_sterile_0 + p_e_no_p_dm_sterile_0))

# e_dm, sterile 0
p_e_dm_p_yes_sterile_0 = 100 * length(tregs_useful_e_dm_p_yes_sterile_0_keep) / total_realizations  # 2.32
p_e_dm_p_no_sterile_0  = 100 * length(tregs_useful_e_dm_p_no_sterile_0_keep)  / total_realizations  # 2.32
p_e_dm_p_dm_sterile_0  = 100 * length(tregs_useful_e_dm_p_dm_sterile_0_keep)  / total_realizations  # 83.72
print(c(p_e_dm_p_yes_sterile_0, p_e_dm_p_no_sterile_0, p_e_dm_p_dm_sterile_0, p_e_dm_p_yes_sterile_0 + p_e_dm_p_no_sterile_0 + p_e_dm_p_dm_sterile_0))

# pathogenic injury (sterile == 1)
p_e_yes_p_dm_sterile_1 = 100 * length(tregs_useful_e_yes_p_dm_sterile_1_keep) / total_realizations  # 2.32
p_e_no_p_dm_sterile_1  = 100 * length(tregs_useful_e_no_p_dm_sterile_1_keep)  / total_realizations  # 2.32
p_e_dm_p_dm_sterile_1  = 100 * length(tregs_useful_e_dm_p_dm_sterile_1_keep)  / total_realizations  # 95.3
print(c(p_e_yes_p_dm_sterile_1, p_e_no_p_dm_sterile_1, p_e_dm_p_dm_sterile_1, p_e_yes_p_dm_sterile_1 + p_e_no_p_dm_sterile_1 + p_e_dm_p_dm_sterile_1))


### FOR STERILE CASES, 55% of cases are useful for epithelium, 0.83% harmful for epithelium, 44.17% don't matter. 


source("./PLOT_FUNCTIONS.R")
dir_name      = './random_simulations'
subdir_name   = '/figures'
dir.create(paste0(dir_name,subdir_name))

plot_ts = 0

if(plot_ts==1){
  subsubdir_name = '/tregs_useful_e_yes_p_dm_sterile_1_keep'
  dir.create(paste0(dir_name,subdir_name,subsubdir_name))
  for (ind in tregs_useful_e_yes_p_dm_sterile_1_keep){
    print(ind)
    realization_data = all_data_merged_keep %>% filter(realization_ind==ind)
    p_epithelium = plot_faceted(realization_data, c("epithelial_healthy", paste0("epithelial_inj_", 1:5)),"Epithelial Cell Dynamics")
    p_microbes   = plot_faceted(realization_data, c("commensal", "pathogen"),"Microbe Dynamics")
    ggsave(paste0(dir_name,subdir_name,subsubdir_name,"/epithelium_", ind,"_.png"), p_epithelium, dpi = 300, width = 16, height = 10, bg = 'white')
    ggsave(paste0(dir_name,subdir_name,subsubdir_name,"/microbes_", ind, "_.png"), p_microbes, dpi = 300, width = 8, height = 5, bg = 'white')
  }
  
  subsubdir_name = '/tregs_useful_e_no_p_dm_sterile_1_keep'
  dir.create(paste0(dir_name,subdir_name,subsubdir_name))
  for (ind in tregs_useful_e_no_p_dm_sterile_1_keep){
    print(ind)
    realization_data = all_data_merged_keep %>% filter(realization_ind==ind)
    p_epithelium = plot_faceted(realization_data, c("epithelial_healthy", paste0("epithelial_inj_", 1:5)),"Epithelial Cell Dynamics")
    p_microbes   = plot_faceted(realization_data, c("commensal", "pathogen"),"Microbe Dynamics")
    ggsave(paste0(dir_name,subdir_name,subsubdir_name,"/epithelium_", ind,"_.png"), p_epithelium, dpi = 300, width = 16, height = 10, bg = 'white')
    ggsave(paste0(dir_name,subdir_name,subsubdir_name,"/microbes_", ind, "_.png"), p_microbes, dpi = 300, width = 8, height = 5, bg = 'white')
  }
  
  subsubdir_name = '/tregs_useful_e_dm_p_dm_sterile_1_keep'
  dir.create(paste0(dir_name,subdir_name,subsubdir_name))
  for (ind in tregs_useful_e_dm_p_dm_sterile_1_keep){
    print(ind)
    realization_data = all_data_merged_keep %>% filter(realization_ind==ind)
    p_epithelium = plot_faceted(realization_data, c("epithelial_healthy", paste0("epithelial_inj_", 1:5)),"Epithelial Cell Dynamics")
    p_microbes   = plot_faceted(realization_data, c("commensal", "pathogen"),"Microbe Dynamics")
    ggsave(paste0(dir_name,subdir_name,subsubdir_name,"/epithelium_", ind,"_.png"), p_epithelium, dpi = 300, width = 16, height = 10, bg = 'white')
    ggsave(paste0(dir_name,subdir_name,subsubdir_name,"/microbes_", ind, "_.png"), p_microbes, dpi = 300, width = 8, height = 5, bg = 'white')
  }
}

all_params_merged_keep_0 = all_params_merged_keep %>% filter(realization_ind %in% tregs_useful_e_yes_p_yes_sterile_0_keep)
if(dim(all_params_merged_keep_0)[1]>0){
  all_params_merged_keep_0$class_1  = 1
  all_params_merged_keep_0$class_2  = 1
  all_params_merged_keep_0$class_3  = 1 # one is useful, other is either useful or dm
  all_params_merged_keep_0$class    = 11
}

all_params_merged_keep_1 = all_params_merged_keep %>% filter(realization_ind %in% tregs_useful_e_yes_p_no_sterile_0_keep)
if(dim(all_params_merged_keep_1)[1]>0){
  all_params_merged_keep_1$class_1  = 1
  all_params_merged_keep_1$class_2  = 0
  all_params_merged_keep_1$class_3  = 0 # at least one is harmful
  all_params_merged_keep_1$class    = 10
}

all_params_merged_keep_2 = all_params_merged_keep %>% filter(realization_ind %in% tregs_useful_e_yes_p_dm_sterile_0_keep)
if(dim(all_params_merged_keep_2)[1]>0){
  all_params_merged_keep_2$class_1  = 1
  all_params_merged_keep_2$class_2  = 2
  all_params_merged_keep_2$class_3  = 1 # one is useful, other is either useful or dm
  all_params_merged_keep_2$class    = 12
}

all_params_merged_keep_3 = all_params_merged_keep %>% filter(realization_ind %in% tregs_useful_e_no_p_yes_sterile_0_keep)
if(dim(all_params_merged_keep_3)[1]>0){
  all_params_merged_keep_3$class_1  = 0
  all_params_merged_keep_3$class_2  = 1
  all_params_merged_keep_3$class_3  = 0 # at least one is harmful
  all_params_merged_keep_3$class    = 1
}

all_params_merged_keep_4 = all_params_merged_keep %>% filter(realization_ind %in% tregs_useful_e_no_p_no_sterile_0_keep)
if(dim(all_params_merged_keep_4)[1]>0){
  all_params_merged_keep_4$class_1  = 0
  all_params_merged_keep_4$class_2  = 0
  all_params_merged_keep_4$class_3  = 0 # at least one is harmful
  all_params_merged_keep_4$class    = 0
}

all_params_merged_keep_5 = all_params_merged_keep %>% filter(realization_ind %in% tregs_useful_e_no_p_dm_sterile_0_keep)
if(dim(all_params_merged_keep_5)[1]>0){
  all_params_merged_keep_5$class_1  = 0
  all_params_merged_keep_5$class_2  = 2
  all_params_merged_keep_5$class_3  = 0 # at least one is harmful
  all_params_merged_keep_5$class    = 2
}

all_params_merged_keep_6 = all_params_merged_keep %>% filter(realization_ind %in% tregs_useful_e_dm_p_yes_sterile_0_keep)
if(dim(all_params_merged_keep_6)[1]>0){
  all_params_merged_keep_6$class_1  = 2
  all_params_merged_keep_6$class_2  = 1
  all_params_merged_keep_6$class_3  = 1 # one is useful, other is either useful or dm
  all_params_merged_keep_6$class    = 21
}

all_params_merged_keep_7 = all_params_merged_keep %>% filter(realization_ind %in% tregs_useful_e_dm_p_no_sterile_0_keep)
if(dim(all_params_merged_keep_7)[1]>0){
  all_params_merged_keep_7$class_1  = 2
  all_params_merged_keep_7$class_2  = 0
  all_params_merged_keep_7$class_3  = 0 # at least one is harmful
  all_params_merged_keep_7$class    = 20
}

all_params_merged_keep_8 = all_params_merged_keep %>% filter(realization_ind %in% tregs_useful_e_dm_p_dm_sterile_0_keep)
if(dim(all_params_merged_keep_8)[1]>0){
  all_params_merged_keep_8$class_1  = 2
  all_params_merged_keep_8$class_2  = 2
  all_params_merged_keep_8$class_3  = 2 # it doesn't matter at all 
  all_params_merged_keep_8$class    = 22
}

sterile_0_df = rbind(all_params_merged_keep_0, all_params_merged_keep_1, all_params_merged_keep_2, all_params_merged_keep_3, all_params_merged_keep_4, 
                     all_params_merged_keep_5, all_params_merged_keep_6, all_params_merged_keep_7, all_params_merged_keep_8)

sterile_0_df$class   = as.factor(sterile_0_df$class)
sterile_0_df$class_1 = as.factor(sterile_0_df$class_1)
sterile_0_df$class_2 = as.factor(sterile_0_df$class_2)
sterile_0_df$class_3 = as.factor(sterile_0_df$class_3)

features = sterile_0_df %>% dplyr::select(-realization_ind, -class, -class_1, -class_2, -class_3, -epith_recovery_chance)

# lda_model = lda(class ~ ., data = cbind(class = sterile_0_df$class, features))
lda_model = lda(class ~ ., data = cbind(class = sterile_0_df$class_3, features))

print(lda_model)

lda_pred = predict(lda_model)
sterile_0_df$LD1 = lda_pred$x[,1]
sterile_0_df$LD2 = lda_pred$x[,2]  # if 3+ classes

lda_scores = as.data.frame(lda_pred$x)
# lda_scores$class = sterile_0_df$class
lda_scores$class = sterile_0_df$class_3

loadings = as.data.frame(lda_model$scaling)
loadings$feature = rownames(loadings)
loadings$LD1 = loadings$LD1 / max(abs(loadings$LD1)) * max(lda_scores$LD1)
loadings$LD2 = loadings$LD2 / max(abs(loadings$LD2)) * max(lda_scores$LD2)

# Compute vector length (Euclidean norm) of each feature's contribution
loadings$importance = sqrt(loadings$LD1^2 + loadings$LD2^2)

# Keep top N features (e.g., top 8)
top_features = loadings %>% filter(importance > 0.25*max(loadings$importance))

lda_scores$class_label = factor(lda_scores$class,
                                levels = c(0, 1, 2),
                                labels = c("At least one is HARMFUL", "At least one is HELPFUL, other DM", "Nothing matters"))

custom_colors <- c(
  "At least one is HARMFUL" = "red",   # purple
  "At least one is HELPFUL, other DM" = "blue",
  "Nothing matters" = "#AFAAB9"   # purple
)

p_lda_sterile_0 = ggplot(lda_scores, aes(x = LD1, y = LD2, color = class_label)) +
  geom_point(alpha = 0.4, size = 2) +
  stat_ellipse(type = "norm", level = 0.95, size = 1, linetype = "dashed") +  # 95% confidence ellipses
  geom_segment(data = top_features,
               aes(x = 0, y = 0, xend = LD1, yend = LD2),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black",
               linewidth = 1.1)+  # thicker arrows
  geom_text(data = top_features,
            aes(x = LD1, y = LD2, label = feature),
            color = "black", size = 4.5, fontface = "bold", hjust = 1.1, vjust = 1.1) +
  scale_color_manual(values = custom_colors) +
  labs(title = "LDA projection of top features, pathogenic injury",
       x = "LD1", y = "LD2", color = "Treg effect class") +
  xlim(-6,6)+
  theme_minimal(base_size = 14)
p_lda_sterile_0

filname = paste0(dir_name,subdir_name,"/LDA_STERILE_0.png")
ggsave(filname, p_lda_sterile_0, width = 16, height = 12, dpi = 300, bg = 'white')


####### STERILE
all_params_merged_keep_9  = all_params_merged_keep %>% filter(realization_ind %in% tregs_useful_e_yes_p_dm_sterile_1_keep)
all_params_merged_keep_9$class  = 1

all_params_merged_keep_10 = all_params_merged_keep %>% filter(realization_ind %in% tregs_useful_e_no_p_dm_sterile_1_keep)
all_params_merged_keep_10$class = 0

all_params_merged_keep_11 = all_params_merged_keep %>% filter(realization_ind %in% tregs_useful_e_dm_p_dm_sterile_1_keep)
all_params_merged_keep_11$class = 2


sterile_1_df = rbind(all_params_merged_keep_9,all_params_merged_keep_10,all_params_merged_keep_11)
sterile_1_df$class = as.factor(sterile_1_df$class)

features = sterile_1_df %>% dplyr::select(-realization_ind, -class, -epith_recovery_chance)
lda_model = lda(class ~ ., data = cbind(class = sterile_1_df$class, features))
print(lda_model)

lda_pred = predict(lda_model)
sterile_1_df$LD1 = lda_pred$x[,1]
sterile_1_df$LD2 = lda_pred$x[,2]  # if 3+ classes

lda_scores = as.data.frame(lda_pred$x)
lda_scores$class = sterile_1_df$class

loadings = as.data.frame(lda_model$scaling)
loadings$feature = rownames(loadings)
loadings$LD1 = loadings$LD1 / max(abs(loadings$LD1)) * max(lda_scores$LD1)
loadings$LD2 = loadings$LD2 / max(abs(loadings$LD2)) * max(lda_scores$LD2)

# Compute vector length (Euclidean norm) of each feature's contribution
loadings$importance = sqrt(loadings$LD1^2 + loadings$LD2^2)

# Keep top N features (e.g., top 8)
top_features = loadings %>% filter(importance > 0.25*max(loadings$importance))

lda_scores$class_label = factor(lda_scores$class,
                                levels = c(0, 1, 2),
                                labels = c("Tregs Harmful", "Tregs Helpful", "Tregs Irrelevant"))

custom_colors <- c(
  "Tregs Helpful"    = "blue",  # green
  "Tregs Harmful"    = "red",  # orange
  "Tregs Irrelevant" = "#AFAAB9"   # purple
)

p_lda_sterile_1 = ggplot(lda_scores, aes(x = LD1, y = LD2, color = class_label)) +
  geom_point(alpha = 0.4, size = 2) +
  stat_ellipse(type = "norm", level = 0.75, size = 1, linetype = "dashed") +  # 95% confidence ellipses
  geom_segment(data = top_features,
               aes(x = 0, y = 0, xend = LD1, yend = LD2),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black",
               linewidth = 1.1)+  # thicker arrows
  geom_text(data = top_features,
            aes(x = LD1, y = LD2, label = feature),
            color = "black", size = 4.5, fontface = "bold", hjust = 1.1, vjust = 1.1) +
  scale_color_manual(values = custom_colors) +
  labs(title = "LDA projection of top features, sterile injury",
       x = "LD1", y = "LD2", color = "Treg effect class") +
  # xlim(-6,6)+
  theme_minimal(base_size = 14)
p_lda_sterile_1

filname = paste0(dir_name,subdir_name,"/LDA_STERILE_1.png")
ggsave(filname, p_lda_sterile_1, width = 16, height = 12, dpi = 300, bg = 'white')
