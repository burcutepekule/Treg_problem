# Clear workspace at the very beginning
rm(list=ls())

# Load required libraries
library(parallel)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)
library(Rcpp)
library(av)

# Set working directory and source functions
args = commandArgs(trailingOnly = FALSE)
script_path = sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0) {
  script_dir = dirname(normalizePath(script_path))
} else {
  script_dir = dirname(rstudioapi::getActiveDocumentContext()$path)
}
setwd(script_dir)

source("./FAST_FUNCTIONS.R")
source("./PLOT_FUNCTIONS.R")

# Create output directory
dir_name = './frames'
dir.create(dir_name)

# Set seed and parameters
seed_in = 3
set.seed(seed_in)
placeholder = nullGrob()

# Define all simulation parameters
t_max           = 500
plot_on         = 1
plot_every      = 1

rate_leak_commensal_injury     = 0.50
rate_leak_commensal_baseline   = 0.05
rat_com_pat_threshold          = 0.95
th_ROS_microbe                 = 0.10
th_ROS_epith_recover           = th_ROS_microbe+0.20
active_age_limit               = 5
epith_recovery_chance          = 0.25
treg_vicinity_effect           = 1

ros_decay   = 0.50
DAMPs_decay = 0.25
SAMPs_decay = 0.25

grid_size            = 25
n_phagocytes         = 200
n_tregs              = 200
cc_phagocyte         = 5
rounds_active        = 5
digestion_time       = 3

injury_percentage = 60
injury_site       = get_middle_percent(seq(1,grid_size), injury_percentage)

diffusion_speed_DAMPs = 0.1
diffusion_speed_SAMPs = 0.1
diffusion_speed_ROS   = 0.1

add_ROS   = .5
add_DAMPs = .5
add_SAMPs = .5

max_cell_value_ROS   = 1
max_cell_value_DAMPs = 1
max_cell_value_SAMPs = 1

lim_ROS    = max_cell_value_ROS
lim_DAMP   = max_cell_value_DAMPs
lim_SAMP   = max_cell_value_SAMPs

act_radius_ROS             = 1
act_radius_treg            = 1
act_radius_DAMPs           = 1
act_radius_SAMPs           = 1

activation_threshold_DAMPs = 0.05
activation_threshold_SAMPs = 0.05

activity_engulf_M0_baseline    = 0.25 
activity_engulf_M1_baseline    = 0.45
activity_engulf_M2_baseline    = 0.45 
activity_engulf_max            = 0.99

activity_ROS_M0_baseline  = 0.00
activity_ROS_M1_baseline  = 0.15
activity_ROS_M2_baseline  = 0.00
activity_ROS_max          = 0.99

activity_engulf_M1_step   = (activity_engulf_max-activity_engulf_M1_baseline)/cc_phagocyte
activity_engulf_M2_step   = (activity_engulf_max-activity_engulf_M2_baseline)/cc_phagocyte
activity_ROS_M1_step      = (activity_ROS_max-activity_ROS_M1_baseline)/cc_phagocyte

# Create scenarios
scenarios_df = expand.grid(
  sterile = c(0, 1),
  allow_tregs_to_do_their_job = c(FALSE, TRUE),
  allow_tregs_to_suppress_cognate = c(FALSE),
  randomize_tregs = c(0,1)
)

# Function to run a single scenario
run_scenario <- function(scenario_ind) {
  
  get_8n_avg_signal_fast <- function(x, y, act_radius_signal, signal_matrix) {
    loc  = c(x, y)
    x_coordinates = (loc[1]-act_radius_signal):(loc[1]+act_radius_signal)
    x_coordinates = x_coordinates[x_coordinates>0 & x_coordinates<=grid_size]
    y_coordinates = (loc[2]-act_radius_signal):(loc[2]+act_radius_signal)
    y_coordinates = y_coordinates[y_coordinates>0 & y_coordinates<=grid_size]
    dval = signal_matrix[y_coordinates, x_coordinates]
    return(mean(dval))
  }
  
  # Fast shift_insert for matrix operations
  shift_insert_fast <- function(vec, insert_vals) {
    n_insert <- length(insert_vals)
    n_vec <- length(vec)
    
    if(n_insert >= n_vec) {
      return(insert_vals[1:n_vec])
    } else {
      return(c(insert_vals, vec[1:(n_vec - n_insert)]))
    }
  }
  
  # Optimized version of iszero_coordinates
  iszero_coordinates <- function(x) {
    # Initialize with default sampling: -1, 0, 1
    y <- sample(c(-1, 0, 1), length(x), replace = TRUE)
    # Replace where x == 0 with sample from -1 or 1
    zero_idx <- which(x == 0)
    y[zero_idx] <- sample(c(-1, 1), length(zero_idx), replace = TRUE)
    
    return(y)
  }
  
  diffuse_matrix <- function(mat, D, max_cell_value) {
    nr <- nrow(mat)
    nc <- ncol(mat)
    
    # Pad the matrix with zeros around the edges
    padded <- matrix(0, nrow = nr + 2, ncol = nc + 2)
    padded[2:(nr + 1), 2:(nc + 1)] <- mat
    
    # # Compute 4-neighbor diffusion
    # laplacian <- padded[1:nr,   2:(nc+1)] +  # up
    #   padded[3:(nr+2), 2:(nc+1)] +  # down
    #   padded[2:(nr+1), 1:nc] +      # left
    #   padded[2:(nr+1), 3:(nc+2)] -  # right
    #   4 * mat                      # center
    
    # Compute 8-neighbor Laplacian (Moore neighborhood)
    laplacian <- (
      padded[1:nr,     1:nc    ] +  # top-left
        padded[1:nr,     2:(nc+1)] +  # top
        padded[1:nr,     3:(nc+2)] +  # top-right
        padded[2:(nr+1), 1:nc    ] +  # left
        padded[2:(nr+1), 3:(nc+2)] +  # right
        padded[3:(nr+2), 1:nc    ] +  # bottom-left
        padded[3:(nr+2), 2:(nc+1)] +  # bottom
        padded[3:(nr+2), 3:(nc+2)]    # bottom-right
      - 8 * mat                   # center subtraction
    )
    
    mat_new <- mat + D * laplacian
    mat_new <- matrix(pmin(max_cell_value, mat_new), nrow = nrow(mat), ncol = ncol(mat))
    
    return(mat_new)
  }
  shift_insert <- function(current_registry, new_elements_vector) {
    combined_registry <- c(new_elements_vector, current_registry)
    target_length <- length(current_registry)
    result_registry <- combined_registry[1:target_length]
    
    return(result_registry)
  }
  
  get_middle_percent <- function(seq_vector, percent) {
    n_total <- length(seq_vector)
    n_select <- ceiling(n_total * percent / 100)
    
    # Calculate start and end index for middle values
    mid <- floor(n_total / 2)
    half_window <- floor(n_select / 2)
    
    start_idx <- max(1, mid - half_window + 1)
    end_idx <- min(n_total, start_idx + n_select - 1)
    
    return(seq_vector[start_idx:end_idx])
  }
  
  
  
  vectors_to_treg_df <- function() {
    data.frame(
      x = treg_x,
      y = treg_y,
      active_age = treg_active_age,
      phenotype = treg_phenotype,
      activity_SAMPs_binary = treg_activity_SAMPs_binary,
      id = treg_id
    )
  }
  
  matrix_to_pathogen_df <- function() {
    if(nrow(pathogen_coords) == 0) {
      return(data.frame(x = numeric(0), y = numeric(0), id = numeric(0)))
    }
    data.frame(
      x = pathogen_coords[, "x"],
      y = pathogen_coords[, "y"], 
      id = pathogen_coords[, "id"]
    )
  }
  
  matrix_to_commensal_df <- function() {
    if(nrow(commensal_coords) == 0) {
      return(data.frame(x = numeric(0), y = numeric(0), id = numeric(0)))
    }
    data.frame(
      x = commensal_coords[, "x"],
      y = commensal_coords[, "y"],
      id = commensal_coords[, "id"]
    )
  }
  
  agent_colors = c(
    epithelial_healthy = "#B0E2FF",
    epithelial_inj_1   = "#8CB4E5",
    epithelial_inj_2   = "#6987CC",
    epithelial_inj_3   = "#465AB2",
    epithelial_inj_4   = "#232D99",
    epithelial_inj_5   = "#000080",
    phagocyte_M0       = "grey70",
    phagocyte_M1_L_0   = "#F8C8E8",
    phagocyte_M1_L_1   = "#F397D6",
    phagocyte_M1_L_2   = "#E754C4",
    phagocyte_M1_L_3   = "#D12CA0",
    phagocyte_M1_L_4   = "#A5177A",
    phagocyte_M1_L_5   = "#6B0C4F",
    phagocyte_M2_L_0   = "#CDEFE3", 
    phagocyte_M2_L_1   = "#97D6BC",  
    phagocyte_M2_L_2   = "#61BD96",  
    phagocyte_M2_L_3   = "#3BA578",  
    phagocyte_M2_L_4   = "#2E8B57",  
    phagocyte_M2_L_5   = "#1F5C3B",  
    phagocyte_M0   = "grey70",
    phagocyte_M1   = "#E754C4",
    phagocyte_M2   = "#3BA578",
    treg_resting = "#D8BFD8",
    treg_active  = "#967BB6",
    commensal    = "turquoise2",
    pathogen     = "firebrick1",
    C_ROS = "black",
    C_M0  = "grey70",
    C_M1  = "#D12CA0",
    C_M2  = "#3BA578",
    P_ROS = "black",
    P_M0  = "grey70",
    P_M1  = "#D12CA0",
    P_M2  = "#3BA578"
  )
  
  
  plot_faceted = function(data, variables, title) {
    data_long = data %>%
      dplyr::select(t, sterile, allow_tregs_to_do_their_job, all_of(variables)) %>%
      pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
    
    p=ggplot(data_long, aes(x = t, y = value, color = variable)) +
      geom_line(alpha = 1, linewidth = 1) +
      facet_grid(sterile ~ allow_tregs_to_do_their_job, labeller = label_both) +
      scale_color_manual(values = agent_colors) +
      theme_minimal() +
      labs(title = title, x = "Time", y = "Count", color = "Agent")
    
    return(p)
  }
  
  plot_faceted_8 = function(data, variables, title) {
    data_long = data %>%
      dplyr::select(t, sterile, allow_tregs_to_do_their_job, randomize_tregs, all_of(variables)) %>%
      pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
    
    p = ggplot(data_long, aes(x = t, y = value, color = variable)) +
      geom_line(alpha = 1, linewidth = 1) +
      facet_grid(sterile ~ allow_tregs_to_do_their_job + randomize_tregs, labeller = label_both) +
      scale_color_manual(values = agent_colors) +
      theme_minimal() +
      labs(title = title, x = "Time", y = "Count", color = "Agent")
    
    return(p)
  }
  
  plot_grid_DAMPs <- function() {
    
    # Create full grid background (invisible or white)
    full_grid <- expand.grid(x = 1:grid_size, y = 1:grid_size)
    
    # --- Convert DAMP matrix to long format ---
    damps_df <- as.data.frame(DAMPs)
    colnames(damps_df) <- paste0("x", 1:ncol(damps_df))
    damps_df <- damps_df %>%
      mutate(y = nrow(damps_df):1) %>%
      pivot_longer(cols = starts_with("x"), names_to = "x", values_to = "value") %>%
      mutate(x = as.integer(gsub("x", "", x)))
    
    # --- Epithelial layer: blue/orange by health ---
    epithelial_df <- epithelium %>%
      mutate(
        type = ifelse(level_injury == 0, "epithelial_healthy",
                      paste0("epithelial_inj_", level_injury)),
        fill = agent_colors[type],
        y = grid_size+ 1 - y
      ) %>%
      select(x, y, fill)
    
    # --- Plot ---
    p_damps = ggplot() +
      geom_tile(data = full_grid, aes(x = x, y = y), fill = "white", color = NA, width = 1, height = 1) +
      # DAMP heatmap with colorbar
      geom_tile(data = damps_df, aes(x = x, y = y, fill = value)) +
      scale_fill_gradient(low = "white", high = "black", name = "DAMPs Density", limits = c(0, lim_DAMP)) +
      
      # Epithelium on top with inline fill colors (no legend)
      geom_tile(data = epithelial_df, aes(x = x, y = y), fill = epithelial_df$fill, width = 1, height = 1) +
      
      coord_fixed() +
      # labs(title = "DAMP Heatmap with Epithelium (No Legend for Epithelium)") +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank()
      )
    return(p_damps)
    
  }
  
  plot_grid_SAMPs <- function() {
    
    # Create full grid background (invisible or white)
    full_grid <- expand.grid(x = 1:grid_size, y = 1:grid_size)
    
    # --- Convert DAMP matrix to long format ---
    samps_df <- as.data.frame(SAMPs)
    colnames(samps_df) <- paste0("x", 1:ncol(samps_df))
    samps_df <- samps_df %>%
      mutate(y = nrow(samps_df):1) %>%
      pivot_longer(cols = starts_with("x"), names_to = "x", values_to = "value") %>%
      mutate(x = as.integer(gsub("x", "", x)))
    
    # --- Epithelial layer: blue/orange by health ---
    epithelial_df <- epithelium %>%
      mutate(
        type = ifelse(level_injury == 0, "epithelial_healthy",
                      paste0("epithelial_inj_", level_injury)),
        fill = agent_colors[type],
        y = grid_size+ 1 - y
      ) %>%
      select(x, y, fill)
    
    # --- Plot ---
    p_samps = ggplot() +
      geom_tile(data = full_grid, aes(x = x, y = y), fill = "white", color = NA, width = 1, height = 1) +
      # DAMP heatmap with colorbar
      geom_tile(data = samps_df, aes(x = x, y = y, fill = value)) +
      scale_fill_gradient(low = "white", high = "black", name = "SAMPs Density", limits = c(0, lim_DAMP)) +
      
      # Epithelium on top with inline fill colors (no legend)
      geom_tile(data = epithelial_df, aes(x = x, y = y), fill = epithelial_df$fill, width = 1, height = 1) +
      
      coord_fixed() +
      # labs(title = "DAMP Heatmap with Epithelium (No Legend for Epithelium)") +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank()
      )
    return(p_samps)
    
  }
  
  plot_grid_ROS <- function() {
    
    # Create full grid background (invisible or white)
    full_grid <- expand.grid(x = 1:grid_size, y = 1:grid_size)
    
    # --- Convert DAMP matrix to long format ---
    ros_df <- as.data.frame(ROS)
    colnames(ros_df) <- paste0("x", 1:ncol(ros_df))
    ros_df <- ros_df %>%
      mutate(y = nrow(ros_df):1) %>%
      pivot_longer(cols = starts_with("x"), names_to = "x", values_to = "value") %>%
      mutate(x = as.integer(gsub("x", "", x)))
    
    # --- Epithelial layer: blue/orange by health ---
    epithelial_df <- epithelium %>%
      mutate(
        type = ifelse(level_injury == 0, "epithelial_healthy",
                      paste0("epithelial_inj_", level_injury)),
        fill = agent_colors[type],
        y = grid_size+ 1 - y
      ) %>%
      select(x, y, fill)
    
    # --- Plot ---
    p_ros = ggplot() +
      geom_tile(data = full_grid, aes(x = x, y = y), fill = "white", color = NA, width = 1, height = 1) +
      # DAMP heatmap with colorbar
      geom_tile(data = ros_df, aes(x = x, y = y, fill = value)) +
      scale_fill_gradient(low = "white", high = "black", name = "ROS Density", limits = c(0, lim_ROS)) +
      
      # Epithelium on top with inline fill colors (no legend)
      geom_tile(data = epithelial_df, aes(x = x, y = y), fill = epithelial_df$fill, width = 1, height = 1) +
      
      coord_fixed() +
      # labs(title = "ROS Heatmap with Epithelium (No Legend for Epithelium)") +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank()
      )
    return(p_ros)
  }
  
  plot_grid_antiinf <- function() {
    
    # Create a base dataframe for the grid
    grid <- expand.grid(x = 1:grid_size, y = 0:grid_size)
    
    # 1. Epithelial cell state (y = 0 only)
    # 1. Epithelial cell state (y = 0 only)
    epithelial_layer <- epithelium %>%
      mutate(type = ifelse(level_injury == 0, "epithelial_healthy",
                           paste0("epithelial_inj_", level_injury)))
    
    # 1. Phagocyte state 
    phagocytes_plot <- phagocytes %>%
      mutate(type = ifelse(phenotype==0, "phagocyte_M0", 
                           ifelse(phenotype==1, "phagocyte_M1","phagocyte_M2")))
    
    phagocytes_plot = phagocytes_plot %>% filter(type %in% c("phagocyte_M2"))
    
    tregs_plot <- tregs %>%
      mutate(type = ifelse(phenotype==0, "treg_resting","treg_active"))
    
    tregs_plot = tregs_plot %>% filter(type=="treg_active")
    
    # Create full grid background (invisible or white)
    full_grid <- expand.grid(x = 1:grid_size, y = 1:grid_size)
    
    # 2. Combine all agents into one dataframe with their type
    all_types <- c(
      "epithelial_healthy", "epithelial_inj_1","epithelial_inj_2","epithelial_inj_3","epithelial_inj_4","epithelial_inj_5",
      "phagocyte_M2","treg_active")
    
    agent_plot_df <- bind_rows(
      epithelial_layer %>% select(x, y, type),
      phagocytes_plot %>% select(x, y, type),
      tregs_plot  %>% select(x, y, type)
    ) %>%
      mutate(type = factor(type, levels = all_types))
    
    p_lym <- ggplot() +
      geom_tile(data = full_grid, aes(x = x, y = y), fill = "white", color = NA, width = 1, height = 1) +
      geom_tile(data = agent_plot_df, aes(x = x, y = y, fill = type), width = 1, height = 1) +
      scale_fill_manual(
        values = agent_colors,
        name = "Cell Type",
        drop = FALSE  # <-- ensures unused levels are shown
      ) +
      coord_fixed(ratio = 1) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, grid_size+1)) +
      scale_y_reverse(expand = c(0, 0)) +  # Flip Y so y=0 is on top
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "right"
      )
    
    return(p_lym)
  }
  
  plot_grid_phagocyte_M1 <- function() {
    
    # Create a base dataframe for the grid
    grid <- expand.grid(x = 1:grid_size, y = 0:grid_size)
    
    # 1. Epithelial cell state (y = 0 only)
    epithelial_layer <- epithelium %>%
      mutate(type = ifelse(level_injury == 0, "epithelial_healthy",
                           paste0("epithelial_inj_", level_injury)))
    # 1. Phagocyte state
    phagocytes_plot <- phagocytes %>%
      mutate(type = ifelse(phenotype==0, "phagocyte_M0", 
                           ifelse(phenotype==1, "phagocyte_M1","phagocyte_M2")))
    
    phagocytes_plot = phagocytes_plot %>% filter(type=="phagocyte_M1")
    # Create full grid background (invisible or white)
    full_grid <- expand.grid(x = 1:grid_size, y = 1:grid_size)
    
    phagocytes_plot$level <- sapply(phagocytes_plot$bacteria_registry, sum)
    
    if(dim(phagocytes_plot)[1]>0){
      phagocytes_plot$type <- paste0("phagocyte_M1_L_", phagocytes_plot$level)
    }
    
    # Combine all agents into one dataframe with their type
    agent_plot_df <- bind_rows(
      epithelial_layer %>% select(x, y, type),
      phagocytes_plot %>% select(x, y, type)
    )
    
    # 2. Combine all agents into one dataframe with their type
    # all_types <- c(
    #   "epithelial_healthy", "epithelial_unhealthy","phagocyte_M1")
    
    all_types <- c(
      "epithelial_healthy", "epithelial_inj_1","epithelial_inj_2","epithelial_inj_3","epithelial_inj_4","epithelial_inj_5",
      "phagocyte_M1_L_0","phagocyte_M1_L_1","phagocyte_M1_L_2","phagocyte_M1_L_3","phagocyte_M1_L_4","phagocyte_M1_L_5")
    
    agent_plot_df <- bind_rows(
      epithelial_layer %>% select(x, y, type),
      phagocytes_plot %>% select(x, y, type)
    ) %>%
      mutate(type = factor(type, levels = all_types))
    
    p_lym <- ggplot() +
      geom_tile(data = full_grid, aes(x = x, y = y), fill = "white", color = NA, width = 1, height = 1) +
      geom_tile(data = agent_plot_df, aes(x = x, y = y, fill = type), width = 1, height = 1) +
      scale_fill_manual(
        values = agent_colors,
        name = "Cell Type",
        drop = FALSE  # <-- ensures unused levels are shown
      ) +
      coord_fixed(ratio = 1) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, grid_size+1)) +
      scale_y_reverse(expand = c(0, 0)) +  # Flip Y so y=0 is on top
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "right"
      )
    
    return(p_lym)
  }
  
  plot_grid_mat = function(mat, grid_size, lim_mat){
    # Create full grid background (invisible or white)
    full_grid <- expand.grid(x = 1:grid_size, y = 1:grid_size)
    
    # --- Convert DAMP matrix to long forDAMPs ---
    mat_df <- as.data.frame(mat)
    colnames(mat_df) <- paste0("x", 1:ncol(mat_df))
    mat_df <- mat_df %>%
      mutate(y = nrow(mat_df):1) %>%
      pivot_longer(cols = starts_with("x"), names_to = "x", values_to = "value") %>%
      mutate(x = as.integer(gsub("x", "", x)))
    
    # --- Plot ---
    p = ggplot() +
      geom_tile(data = full_grid, aes(x = x, y = y), fill = "white", color = NA, width = 1, height = 1) +
      geom_tile(data = mat_df, aes(x = x, y = y, fill = value)) +
      scale_fill_gradient(low = "white", high = "black", name = "DAMP Density", limits = c(0, lim_mat)) +
      coord_fixed() +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank()
      )
    return(p)
  }
  
  plot_microbes_nomem = function() {
    # Count microbes
    num_pat = nrow(pathogens_lp)
    num_com = nrow(commensals_lp)
    
    # Create a new row of data
    new_data = data.frame(
      time = t,
      count = c(num_pat, num_com),
      type = c("pathogen", "commensal")
    )
    
    microbe_colors <- agent_colors[c("pathogen", "commensal")]
    
    
    # Horizontal bar plot
    p = ggplot(new_data, aes(x = type, y = count, fill = type)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = microbe_colors, name = "Microbe Type") +
      coord_flip() +  # Make bars horizontal
      ylim(0, max(200, num_pat, num_com)) +
      # labs(x = "", y = "Microbe Count", title = paste("Time =", t)) +
      theme_minimal() +
      theme(legend.position = "none")
    
    
    return(p)
  }
  
  plot_phagocyte_nomem = function() {
    # Count lymph
    # 1. Phagocyte state 
    phagocytes_plot <- phagocytes %>%
      mutate(type = ifelse(phenotype==0, "phagocyte_M0", 
                           ifelse(phenotype==1, "phagocyte_M1","phagocyte_M2")))
    
    num_M0 = dim(phagocytes_plot %>% filter(type=="phagocyte_M0"))[1]
    num_M1 = dim(phagocytes_plot %>% filter(type=="phagocyte_M1"))[1]
    num_M2 = dim(phagocytes_plot %>% filter(type=="phagocyte_M2"))[1]
    
    # Create a new row of data
    new_data = data.frame(
      time = t,
      count = c(num_M0, num_M1, num_M2),
      type = c("phagocyte_M0", "phagocyte_M1", "phagocyte_M2")
    )
    
    cell_colors <- agent_colors[c("phagocyte_M0", "phagocyte_M1", "phagocyte_M2")]
    
    
    # Horizontal bar plot
    p = ggplot(new_data, aes(x = type, y = count, fill = type)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = cell_colors, name = "Cell Type") +
      coord_flip() +  # Make bars horizontal
      ylim(0, max(30, num_M0, num_M1, num_M2)) +
      theme_minimal() +
      theme(legend.position = "none")
    return(p)
  }
  
  plot_treg_nomem = function() {
    # Count lymph
    # Treg state 
    tregs_plot <- tregs %>%
      mutate(type = ifelse(phenotype==0, "treg_resting","treg_active"))
    
    num_treg_0 = dim(tregs_plot %>% filter(type=="treg_resting"))[1]
    num_treg_1 = dim(tregs_plot %>% filter(type=="treg_active"))[1]
    
    # Create a new row of data
    new_data = data.frame(
      time = t,
      count = c(num_treg_0, num_treg_1),
      type = c("treg_resting", "treg_active")
    )
    
    cell_colors <- agent_colors[c("treg_resting", "treg_active")]
    
    
    # Horizontal bar plot
    p = ggplot(new_data, aes(x = type, y = count, fill = type)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = cell_colors, name = "Cell Type") +
      coord_flip() +  # Make bars horizontal
      ylim(0, max(30, num_treg_0, num_treg_1)) +
      theme_minimal() +
      theme(legend.position = "none")
    return(p)
  }
  
  plot_cumdeath_nomem = function() {
    # Count microbes
    
    # Create a new row of data
    new_data = data.frame(
      time = t,
      count = c(pathogens_killed_by_ROS, pathogens_killed_by_Mac,commensals_killed_by_ROS, commensals_killed_by_Mac),
      type = c("P_ROS","P_M0","P_M1","P_M2","C_ROS","C_M0","C_M1","C_M2")
    )
    
    microbe_colors <- agent_colors[c("P_ROS","P_M0","P_M1","P_M2","C_ROS","C_M0","C_M1","C_M2")]
    
    
    # Horizontal bar plot
    p = ggplot(new_data, aes(x = type, y = count, fill = type)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = microbe_colors, name = "Microbe Type") +
      coord_flip() +  # Make bars horizontal
      ylim(0, max(200, pathogens_killed_by_ROS, pathogens_killed_by_Mac,commensals_killed_by_ROS, commensals_killed_by_Mac)) +
      labs(x = "", y = "Cum. death") +
      theme_minimal() +
      theme(legend.position = "none")
    
    
    return(p)
  }
  
  plot_grid_pathogens <- function() {
    
    # Create a base dataframe for the grid
    grid <- expand.grid(x = 1:grid_size, y = 0:grid_size)
    
    # 1. Epithelial cell state (y = 0 only)
    # 1. Epithelial cell state (y = 0 only)
    epithelial_layer <- epithelium %>%
      mutate(type = ifelse(level_injury == 0, "epithelial_healthy",
                           paste0("epithelial_inj_", level_injury)))
    
    # Create full grid background (invisible or white)
    full_grid <- expand.grid(x = 1:grid_size, y = 1:grid_size)
    
    # 2. Combine all agents into one dataframe with their type
    all_types <- c(
      "epithelial_healthy", "epithelial_inj_1","epithelial_inj_2","epithelial_inj_3","epithelial_inj_4","epithelial_inj_5",
      "pathogen"
    )
    
    pathogen_density <- pathogens_lp %>%
      dplyr::group_by(x, y) %>%
      dplyr::summarise(count = n(), .groups = "drop") %>%
      mutate(
        type = "pathogen"
      )
    
    # Check if there is variation in counts
    if (nrow(pathogen_density) > 0) {
      if (length(unique(pathogen_density$count)) == 1) {
        # All counts are the same: use fixed alpha
        pathogen_density$alpha_val <- 0.8
      } else {
        # Vary alpha by density
        pathogen_density$alpha_val <- pmin(1, pathogen_density$count / max(pathogen_density$count))
      }
    }
    
    epithelial_layer_plot <- epithelial_layer %>%
      select(x, y, type) %>%
      mutate(alpha_val = 1)
    
    agent_plot_df <- bind_rows(
      epithelial_layer_plot,
      pathogen_density
    ) %>%
      mutate(type = factor(type, levels = all_types))
    
    
    p_mic <- ggplot() +
      geom_tile(data = full_grid, aes(x = x, y = y), fill = "white", color = NA, width = 1, height = 1) +
      geom_tile(data = agent_plot_df,
                aes(x = x, y = y, fill = type, alpha = alpha_val),
                width = 1, height = 1) +
      scale_fill_manual(
        values = agent_colors,
        name = "Cell Type",
        drop = FALSE
      ) +
      scale_alpha_continuous(range = c(0.5, 1), guide = FALSE) +
      coord_fixed(ratio = 1) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, grid_size + 1)) +
      scale_y_reverse(expand = c(0, 0)) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "right"
      )
    
    return(p_mic)
  }
  
  plot_grid_commensals <- function() {
    
    # Create a base dataframe for the grid
    grid <- expand.grid(x = 1:grid_size, y = 0:grid_size)
    
    # 1. Epithelial cell state (y = 0 only)
    epithelial_layer <- epithelium %>%
      mutate(type = ifelse(level_injury == 0, "epithelial_healthy",
                           paste0("epithelial_inj_", level_injury)))
    
    # Create full grid background (invisible or white)
    full_grid <- expand.grid(x = 1:grid_size, y = 1:grid_size)
    
    # 2. Combine all agents into one dataframe with their type
    all_types <- c(
      "epithelial_healthy", "epithelial_inj_1","epithelial_inj_2","epithelial_inj_3","epithelial_inj_4","epithelial_inj_5",
      "commensal"
    )
    
    commensal_density <- commensals_lp %>%
      dplyr::group_by(x, y) %>%
      dplyr::summarise(count = n(), .groups = "drop") %>%
      mutate(
        type = "commensal"
      )
    
    # Check if there is variation in counts
    if (nrow(commensal_density) > 0) {
      if (length(unique(commensal_density$count)) == 1) {
        # All counts are the same: use fixed alpha
        commensal_density$alpha_val <- 0.8
      } else {
        # Vary alpha by density
        commensal_density$alpha_val <- pmin(1, commensal_density$count / max(commensal_density$count))
      }
    }
    
    epithelial_layer_plot <- epithelial_layer %>%
      select(x, y, type) %>%
      mutate(alpha_val = 1)
    
    agent_plot_df <- bind_rows(
      epithelial_layer_plot,
      commensal_density
    ) %>%
      mutate(type = factor(type, levels = all_types))
    
    
    p_mic <- ggplot() +
      geom_tile(data = full_grid, aes(x = x, y = y), fill = "white", color = NA, width = 1, height = 1) +
      geom_tile(data = agent_plot_df,
                aes(x = x, y = y, fill = type, alpha = alpha_val),
                width = 1, height = 1) +
      scale_fill_manual(
        values = agent_colors,
        name = "Cell Type",
        drop = FALSE
      ) +
      scale_alpha_continuous(range = c(0.5, 1), guide = FALSE) +
      coord_fixed(ratio = 1) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, grid_size + 1)) +
      scale_y_reverse(expand = c(0, 0)) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "right"
      )
    
    return(p_mic)
  }
  
  plot_grid_resting <- function() {
    
    # Create a base dataframe for the grid
    grid <- expand.grid(x = 1:grid_size, y = 0:grid_size)
    
    # 1. Epithelial cell state (y = 0 only)
    # 1. Epithelial cell state (y = 0 only)
    epithelial_layer <- epithelium %>%
      mutate(type = ifelse(level_injury == 0, "epithelial_healthy",
                           paste0("epithelial_inj_", level_injury)))
    # 1. Phagocyte state 
    phagocytes_plot <- phagocytes %>%
      mutate(type = ifelse(phenotype==0, "phagocyte_M0", 
                           ifelse(phenotype==1, "phagocyte_M1","phagocyte_M2")))
    
    phagocytes_plot = phagocytes_plot %>% filter(type=="phagocyte_M0")
    
    tregs_plot <- tregs %>%
      mutate(type = ifelse(phenotype==0, "treg_resting","treg_active"))
    tregs_plot = tregs_plot %>% filter(type=="treg_resting")
    
    # Create full grid background (invisible or white)
    full_grid <- expand.grid(x = 1:grid_size, y = 1:grid_size)
    
    # 2. Combine all agents into one dataframe with their type
    all_types <- c(
      "epithelial_healthy", "epithelial_inj_1","epithelial_inj_2","epithelial_inj_3","epithelial_inj_4","epithelial_inj_5",
      "phagocyte_M0","treg_resting"
    )
    
    agent_plot_df <- bind_rows(
      epithelial_layer %>% select(x, y, type),
      phagocytes_plot %>% select(x, y, type),
      tregs_plot %>% select(x, y, type)
    ) %>%
      mutate(type = factor(type, levels = all_types))
    
    
    
    p_lym <- ggplot() +
      geom_tile(data = full_grid, aes(x = x, y = y), fill = "white", color = NA, width = 1, height = 1) +
      geom_tile(data = agent_plot_df, aes(x = x, y = y, fill = type), width = 1, height = 1) +
      scale_fill_manual(
        values = agent_colors,
        name = "Cell Type",
        drop = FALSE  # <-- ensures unused levels are shown
      ) +
      coord_fixed(ratio = 1) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, grid_size+1)) +
      scale_y_reverse(expand = c(0, 0)) +  # Flip Y so y=0 is on top
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "right"
      )
    
    return(p_lym)
  }
  
  plot_simtime_simple = function(){
    
    p_DAMPs  = plot_grid_DAMPs()
    p_SAMPs  = plot_grid_SAMPs()
    p_ROS    = plot_grid_ROS()
    p_d      = plot_grid_antiinf()
    p_a      = plot_grid_phagocyte_M1()
    p_lymp   = plot_grid_resting()
    p_com    = plot_grid_commensals()
    p_pat    = plot_grid_pathogens()
    
    p_mic_bar  = plot_microbes_nomem()
    p_lym_bar  = plot_phagocyte_nomem()
    p_treg_bar = plot_treg_nomem()
    p_cumdeath = plot_cumdeath_nomem()
    
    if(sterile==0){
      injury_type= 'pathogenic'
    }else{
      injury_type= 'sterile'
    }
    
    tit_add = paste0('Sim. time : ',t,', Injury: ',injury_type,', Tregs allowed : ',allow_tregs_to_do_their_job)
    row_0 = plot_grid(p_mic_bar,p_lym_bar,p_treg_bar, p_cumdeath, ncol = 4, rel_widths =  c(1,1,1,1))
    row_1 = plot_grid(p_com,p_pat,p_a,p_d,ncol = 4, rel_widths =  c(1,1,1.1,1.05))
    row_2 = plot_grid(p_DAMPs,p_SAMPs,p_ROS,p_lymp,ncol = 4, rel_widths =  c(1,1,0.98,1.1))
    
    combined_plot = plot_grid(row_0,row_1,row_2,align='v',nrow = 3, rel_heights =  c(0.35,1,1))
    # combined_plot = plot_grid(row_1,row_2,align='v',nrow = 2, rel_heights =  c(1,1))
    
    # Add a title
    combined_plot = plot_grid(
      ggdraw() + draw_label(tit_add, size = 14, hjust = 0.5),
      combined_plot,
      ncol = 1,
      rel_heights = c(0.05, 1)  # Adjust title vs. plot height ratio
    )
    
    # print(combined_plot)
    return(combined_plot)
  }
  ###############
  # Get scenario parameters
  sterile                         = scenarios_df[scenario_ind,]$sterile
  allow_tregs_to_do_their_job     = scenarios_df[scenario_ind,]$allow_tregs_to_do_their_job
  allow_tregs_to_suppress_cognate = scenarios_df[scenario_ind,]$allow_tregs_to_suppress_cognate
  randomize_tregs                 = scenarios_df[scenario_ind,]$randomize_tregs
  
  if(sterile==0){
    rate_leak_pathogen_injury      = 0.75
  }else{
    rate_leak_pathogen_injury      = 0.00
  }
  n_commensals_lp = 20
  n_pathogens_lp  = round(rate_leak_pathogen_injury*length(injury_site)) 
  
  # Initialize fields
  DAMPs  = matrix(0, grid_size, grid_size)
  SAMPs  = matrix(0, grid_size, grid_size)
  ROS    = matrix(0, grid_size, grid_size)
  
  epithelium  = data.frame(x = seq(1,grid_size,1),
                           y = rep(0, grid_size),
                           level_injury = 0, id = seq(1,grid_size))
  
  epithelium[injury_site,]$level_injury = 1
  max_level_injury = 5
  
  # Initialize phagocytes
  phagocyte_x = sample(1:grid_size, n_phagocytes, TRUE)
  phagocyte_y = sample(2:grid_size, n_phagocytes, TRUE) 
  phagocyte_pathogens_engulfed  = rep(0, n_phagocytes)
  phagocyte_commensals_engulfed = rep(0, n_phagocytes)
  phagocyte_num_times_activated = rep(0, n_phagocytes)
  phagocyte_phenotype       = rep(0, n_phagocytes)
  phagocyte_activity_ROS    = rep(activity_ROS_M0_baseline, n_phagocytes)
  phagocyte_activity_engulf = rep(activity_engulf_M0_baseline, n_phagocytes)
  phagocyte_active_age      = rep(0, n_phagocytes)
  
  phagocyte_bacteria_registry = matrix(0, nrow = n_phagocytes, ncol = cc_phagocyte)
  
  # Initialize tregs
  treg_x          = sample(1:grid_size, n_tregs, TRUE)
  treg_y          = sample(2:grid_size, n_tregs, TRUE)
  treg_active_age = rep(0, n_tregs)
  treg_phenotype  = rep(0, n_tregs)
  treg_activity_SAMPs_binary = rep(0, n_tregs)
  
  # Initialize pathogens
  if(n_pathogens_lp == 0) {
    pathogen_coords <- matrix(numeric(0), ncol = 3)
    colnames(pathogen_coords) <- c("x", "y", "id")
  } else {
    pathogen_coords <- matrix(c(
      sample(injury_site, n_pathogens_lp, TRUE),
      rep(1, n_pathogens_lp),
      seq(1, n_pathogens_lp)
    ), ncol = 3)
    colnames(pathogen_coords) <- c("x", "y", "id")
  }
  
  # Initialize commensals
  commensal_coords <- matrix(c(
    sample(1:grid_size, n_commensals_lp, TRUE),
    sample(1:grid_size, n_commensals_lp, TRUE),
    seq(1, n_commensals_lp)
  ), ncol = 3)
  colnames(commensal_coords) <- c("x", "y", "id")
  
  last_id_pathogen  = n_pathogens_lp
  last_id_commensal = n_commensals_lp
  
  pathogens_killed_by_ROS = 0
  pathogens_killed_by_Mac = rep(0,3)
  
  commensals_killed_by_ROS = 0
  commensals_killed_by_Mac = rep(0,3)
  
  # Initialize tracking matrices
  epithelium_longitudinal  = matrix(0,nrow=t_max, ncol=(max_level_injury+1))
  macrophages_longitudinal = matrix(0,nrow=t_max, ncol=1+2*(cc_phagocyte+1))
  microbes_longitudinal    = matrix(0,nrow=t_max, ncol=2)
  tregs_longitudinal       = matrix(0,nrow=t_max, ncol=2)
  microbes_cumdeath_longitudinal = matrix(0,nrow=t_max, ncol=2*4)
  
  # Main simulation loop
  for (t in 1:t_max) {
    print(c(scenario_ind,t))
    
    # Update injury site
    injury_site_updated = which(epithelium$level_injury>0)
    
    # Update DAMPs
    DAMPs[1,] = DAMPs[1,] + epithelium$level_injury*add_DAMPs
    
    # Update SAMPs based on activated tregs
    active_tregs = which(treg_phenotype == 1)
    if(length(active_tregs) > 0) {
      for(i in active_tregs) {
        SAMPs[treg_y[i], treg_x[i]] = SAMPs[treg_y[i], treg_x[i]] + treg_activity_SAMPs_binary[i] * add_SAMPs
      }
    }
    
    # Update ROS based on M1 phagocytes  
    M1_phagocytes = which(phagocyte_phenotype == 1)
    if(length(M1_phagocytes) > 0) {
      for(i in M1_phagocytes) {
        ROS[phagocyte_y[i], phagocyte_x[i]] = ROS[phagocyte_y[i], phagocyte_x[i]] + phagocyte_activity_ROS[i] * add_ROS
      }
    }
    
    # Diffuse & decay DAMPs, SAMPs, and ROS
    DAMPs   = diffuse_matrix(DAMPs, diffusion_speed_DAMPs, max_cell_value_DAMPs)
    SAMPs   = diffuse_matrix(SAMPs, diffusion_speed_SAMPs, max_cell_value_SAMPs)
    ROS     = diffuse_matrix(ROS, diffusion_speed_ROS, max_cell_value_ROS)
    DAMPs   = DAMPs - DAMPs_decay*DAMPs
    SAMPs   = SAMPs - SAMPs_decay*SAMPs
    ROS     = ROS - ros_decay*ROS
    
    # Move pathogens and commensals randomly
    if(nrow(pathogen_coords) > 0) {
      dy = ifelse(pathogen_coords[, "y"] == 1,
                  sample(c(1), size = nrow(pathogen_coords), replace = TRUE),
                  sample(c(-1, 0, 1), size = nrow(pathogen_coords), replace = TRUE))
      
      dx = iszero_coordinates(dy)
      
      pathogen_coords[, "x"] = pmin(pmax(pathogen_coords[, "x"] + dx, 1), grid_size)
      pathogen_coords[, "y"] = pmin(pmax(pathogen_coords[, "y"] + dy, 1), grid_size)
    }
    if(nrow(commensal_coords) > 0) {
      dy = ifelse(commensal_coords[, "y"] == 1,
                  sample(c(1), size = nrow(commensal_coords), replace = TRUE),
                  sample(c(-1, 0, 1), size = nrow(commensal_coords), replace = TRUE))
      
      dx = iszero_coordinates(dy)
      
      commensal_coords[, "x"] = pmin(pmax(commensal_coords[, "x"] + dx, 1), grid_size)
      commensal_coords[, "y"] = pmin(pmax(commensal_coords[, "y"] + dy, 1), grid_size)   
    }
    
    # Move phagocytes and tregs based on DAMPs gradient
    density_matrix_tregs      = DAMPs
    if(randomize_tregs==1){
      density_matrix_tregs = matrix(0,grid_size,grid_size)
    }
    density_matrix_phagocytes = DAMPs
    
    all_equal_treg            = all(density_matrix_tregs == density_matrix_tregs[1, 1])
    all_equal_phagocytes      = all(density_matrix_phagocytes == density_matrix_phagocytes[1, 1])
    
    if(!all_equal_treg){ 
      for(i in 1:length(treg_x)) {
        x = treg_x[i]
        y = treg_y[i]
        
        x_range = max(1, x - 1):min(grid_size, x + 1)
        y_range = max(1, y - 1):min(grid_size, y + 1)
        
        neighbors_x = rep(x_range, each = length(y_range))
        neighbors_y = rep(y_range, times = length(x_range))
        
        neighbor_densities = density_matrix_tregs[cbind(neighbors_y, neighbors_x)]
        
        total = sum(neighbor_densities)
        if(total > 0) {
          probs = neighbor_densities / total
        } else {
          probs = rep(1 / length(neighbor_densities), length(neighbor_densities))
        }
        
        chosen_idx = sample(1:length(neighbors_x), 1, prob = probs)
        
        treg_x[i] = neighbors_x[chosen_idx]
        treg_y[i] = neighbors_y[chosen_idx]
      }
    } else {
      dy_treg = ifelse(treg_y == 1,
                       sample(c(1), size = length(treg_y), replace = TRUE),
                       sample(c(-1, 0, 1), size = length(treg_y), replace = TRUE))
      dx_treg = iszero_coordinates(dy_treg)
      treg_x = pmin(pmax(treg_x + dx_treg, 1), grid_size)
      treg_y = pmin(pmax(treg_y + dy_treg, 1), grid_size)
    }
    
    if(!all_equal_phagocytes){ 
      for(i in 1:length(phagocyte_x)){
        x = phagocyte_x[i]
        y = phagocyte_y[i]
        
        x_range = max(1, x - 1):min(grid_size, x + 1)
        y_range = max(1, y - 1):min(grid_size, y + 1)
        
        neighbors_x = rep(x_range, each = length(y_range))
        neighbors_y = rep(y_range, times = length(x_range))
        
        neighbor_densities = density_matrix_phagocytes[cbind(neighbors_y, neighbors_x)]
        
        total = sum(neighbor_densities)
        if(total > 0) {
          probs = neighbor_densities / total
        } else {
          probs = rep(1 / length(neighbor_densities), length(neighbor_densities))
        }
        
        chosen_idx = sample(1:length(neighbors_x), 1, prob = probs)
        
        phagocyte_x[i] = neighbors_x[chosen_idx]
        phagocyte_y[i] = neighbors_y[chosen_idx]
      }
    } else {
      dy_phagocyte = ifelse(phagocyte_y == 1,
                            sample(c(1), size = length(phagocyte_y), replace = TRUE),
                            sample(c(-1, 0, 1), size = length(phagocyte_y), replace = TRUE))
      dx_phagocyte = iszero_coordinates(dy_phagocyte)
      phagocyte_x = pmin(pmax(phagocyte_x + dx_phagocyte, 1), grid_size)
      phagocyte_y = pmin(pmax(phagocyte_y + dy_phagocyte, 1), grid_size)
    }
    
    # Add new microbes based on the injured epithelium
    n_pathogens_lp_new = round(mean(epithelium$level_injury) * rate_leak_pathogen_injury * length(injury_site_updated))
    if(n_pathogens_lp_new > 0) {
      new_pathogen_coords = matrix(c(
        sample(1:grid_size, n_pathogens_lp_new, replace = TRUE, prob = epithelium$level_injury),
        rep(1, n_pathogens_lp_new),
        last_id_pathogen + seq(1, n_pathogens_lp_new)
      ), ncol = 3)
      colnames(new_pathogen_coords) = c("x", "y", "id")
      
      pathogen_coords = rbind(pathogen_coords, new_pathogen_coords)
      last_id_pathogen = last_id_pathogen + n_pathogens_lp_new
    }
    
    n_commensals_lp_new_injury   = round(mean(epithelium$level_injury) * rate_leak_commensal_injury * length(injury_site_updated))
    n_commensals_lp_new_baseline = round(rate_leak_commensal_baseline * grid_size)
    
    total_new_commensals = n_commensals_lp_new_baseline + n_commensals_lp_new_injury
    if(total_new_commensals > 0) {
      baseline_x = sample(1:grid_size, n_commensals_lp_new_baseline, TRUE)
      injury_x = if(n_commensals_lp_new_injury > 0) {
        sample(1:grid_size, n_commensals_lp_new_injury, TRUE, prob = epithelium$level_injury)
      } else { numeric(0) }
      
      new_commensal_coords = matrix(c(
        c(baseline_x, injury_x),
        rep(1, total_new_commensals),
        last_id_commensal + seq(1, total_new_commensals)
      ), ncol = 3)
      colnames(new_commensal_coords) = c("x", "y", "id")
      
      commensal_coords = rbind(commensal_coords, new_commensal_coords)
      last_id_commensal = last_id_commensal + total_new_commensals
    }
    
    # Plotting
    if(plot_on==1 & (t%%plot_every == 0 | t==1)){
      
      
      # Convert back to dataframes when needed for analysis/plotting
      bacteria_registry_list <- vector("list", length(phagocyte_x))
      for(i in 1:length(phagocyte_x)) {
        bacteria_registry_list[[i]] <- phagocyte_bacteria_registry[i, ]
      }
      
      # Create the complete dataframe with bacteria registry
      phagocytes <- data.frame(
        x = phagocyte_x,
        y = phagocyte_y,
        pathogens_engulfed = phagocyte_pathogens_engulfed,
        commensals_engulfed = phagocyte_commensals_engulfed,
        num_times_activated = phagocyte_num_times_activated,
        phenotype = phagocyte_phenotype,
        activity_ROS = phagocyte_activity_ROS,
        activity_engulf = phagocyte_activity_engulf,
        active_age = phagocyte_active_age
      )
      
      # Add the bacteria registry as a list column
      phagocytes$bacteria_registry <- bacteria_registry_list
      
      tregs=data.frame(
        x = treg_x,
        y = treg_y,
        active_age = treg_active_age,
        phenotype = treg_phenotype,
        activity_SAMPs_binary = treg_activity_SAMPs_binary
      )
      
      
      if(nrow(pathogen_coords) == 0) {
        pathogens_lp = data.frame(x = numeric(0), y = numeric(0), id = numeric(0))
      }else{
        pathogens_lp =   data.frame(
          x = pathogen_coords[, "x"],
          y = pathogen_coords[, "y"]
        )
      }
      
      if(nrow(commensal_coords) == 0) {
        commensals_lp = data.frame(x = numeric(0), y = numeric(0), id = numeric(0))
      }else{
        commensals_lp =   data.frame(
          x = commensal_coords[, "x"],
          y = commensal_coords[, "y"]
        )
      }

      p = plot_simtime_simple()
      ggsave(
        paste0(dir_name,"/frame_seed_",seed_in,"_STERILE_",sterile,"_TREGS_",allow_tregs_to_do_their_job,"_trnd_",randomize_tregs,"_",t,".png"),
        plot = p,
        width = 12,
        height = 10,
        dpi = 600,
        bg = "white"
      )
    }
    
    # Update phagocyte phenotypes
    M0_indices = which(phagocyte_phenotype == 0)
    M1_indices = which(phagocyte_phenotype == 1)
    M2_indices = which(phagocyte_phenotype == 2)
    
    if(t %% digestion_time == 0) {
      phagocyte_bacteria_registry = cbind(
        matrix(0, nrow = nrow(phagocyte_bacteria_registry), ncol = 1),
        phagocyte_bacteria_registry[, -ncol(phagocyte_bacteria_registry)]
      )
    }
    
    if(length(M0_indices) > 0) {
      for(i in M0_indices) {
        avg_DAMPs = get_8n_avg_signal_fast(phagocyte_x[i], phagocyte_y[i], act_radius_DAMPs, DAMPs)
        avg_SAMPs = get_8n_avg_signal_fast(phagocyte_x[i], phagocyte_y[i], act_radius_SAMPs, SAMPs)
        
        if(avg_DAMPs >= activation_threshold_DAMPs && avg_DAMPs > avg_SAMPs) {
          phagocyte_phenotype[i] = 1
          phagocyte_active_age[i] = 1
          bacteria_count = sum(phagocyte_bacteria_registry[i, ])
          phagocyte_activity_ROS[i] = activity_ROS_M1_baseline + activity_ROS_M1_step * bacteria_count
          phagocyte_activity_engulf[i] = activity_engulf_M1_baseline + activity_engulf_M1_step * bacteria_count
        } else if(avg_SAMPs >= activation_threshold_SAMPs && avg_SAMPs > avg_DAMPs) {
          phagocyte_phenotype[i] = 2
          phagocyte_active_age[i] = 1
          bacteria_count = sum(phagocyte_bacteria_registry[i, ])
          phagocyte_activity_ROS[i] = activity_ROS_M2_baseline
          phagocyte_activity_engulf[i] = activity_engulf_M2_baseline + activity_engulf_M2_step * bacteria_count
        }
      }
    }
    
    active_indices = c(M1_indices, M2_indices)
    if(length(active_indices) > 0) {
      phagocyte_active_age[active_indices] = phagocyte_active_age[active_indices] + 1
      
      old_enough = phagocyte_active_age[active_indices] >= active_age_limit
      candidates = active_indices[old_enough]
      
      for(i in candidates) {
        avg_DAMPs = get_8n_avg_signal_fast(phagocyte_x[i], phagocyte_y[i], act_radius_DAMPs, DAMPs)
        avg_SAMPs = get_8n_avg_signal_fast(phagocyte_x[i], phagocyte_y[i], act_radius_SAMPs, SAMPs)
        
        bacteria_count = sum(phagocyte_bacteria_registry[i, ])
        
        if(avg_DAMPs >= activation_threshold_DAMPs && avg_DAMPs > avg_SAMPs) {
          phagocyte_phenotype[i] = 1
          phagocyte_active_age[i] = 1
          phagocyte_activity_ROS[i] = activity_ROS_M1_baseline + activity_ROS_M1_step * bacteria_count
          phagocyte_activity_engulf[i] = activity_engulf_M1_baseline + activity_engulf_M1_step * bacteria_count
        } else if(avg_SAMPs >= activation_threshold_SAMPs && avg_SAMPs > avg_DAMPs) {
          phagocyte_phenotype[i] = 2
          phagocyte_active_age[i] = 1
          phagocyte_activity_ROS[i] = activity_ROS_M2_baseline
          phagocyte_activity_engulf[i] = activity_engulf_M2_baseline + activity_engulf_M2_step * bacteria_count
        } else if(avg_SAMPs < activation_threshold_SAMPs && avg_DAMPs < activation_threshold_DAMPs) {
          phagocyte_phenotype[i] = 0
          phagocyte_active_age[i] = 0
          phagocyte_activity_ROS[i] = activity_ROS_M0_baseline
          phagocyte_activity_engulf[i] = activity_engulf_M0_baseline
        }
      }
    }
    
    # Treg active age updated
    active_treg_indices = which(treg_phenotype == 1)
    if(length(active_treg_indices) > 0) {
      old_tregs = active_treg_indices[treg_active_age[active_treg_indices] >= active_age_limit]
      young_tregs = active_treg_indices[treg_active_age[active_treg_indices] < active_age_limit]
      
      if(length(young_tregs) > 0) {
        treg_active_age[young_tregs] = treg_active_age[young_tregs] + 1
      }
      
      if(length(old_tregs) > 0) {
        treg_phenotype[old_tregs] = 0
        treg_active_age[old_tregs] = 0
        treg_activity_SAMPs_binary[old_tregs] = 0
      }
    }
    
    # Engulfment process
    phagocyte_positions = paste(phagocyte_x, phagocyte_y, sep = "_")
    
    for(i in 1:length(phagocyte_x)) {
      px = phagocyte_x[i]
      py = phagocyte_y[i]
      
      if(nrow(pathogen_coords) > 0) {
        pathogen_overlap = (pathogen_coords[, "x"] == px) & (pathogen_coords[, "y"] == py)
        pathogen_indices = which(pathogen_overlap)
        
        if(length(pathogen_indices) > 0) {
          engulf_success = runif(length(pathogen_indices)) < phagocyte_activity_engulf[i]
          indices_to_engulf = pathogen_indices[engulf_success]
          
          if(length(indices_to_engulf) > 0) {
            phagocyte_pathogens_engulfed[i] = phagocyte_pathogens_engulfed[i] + length(indices_to_engulf)
            
            pathogen_coords = pathogen_coords[-indices_to_engulf, , drop = FALSE]
            
            phagocyte_bacteria_registry[i, ] = shift_insert_fast(phagocyte_bacteria_registry[i, ], 
                                                                 rep(1, length(indices_to_engulf)))
            
            phagocyte_phenotype_index = phagocyte_phenotype[i] + 1
            pathogens_killed_by_Mac[phagocyte_phenotype_index] = pathogens_killed_by_Mac[phagocyte_phenotype_index] + length(indices_to_engulf)
          }
        }
      }
      
      if(nrow(commensal_coords) > 0) {
        commensal_overlap = (commensal_coords[, "x"] == px) & (commensal_coords[, "y"] == py)
        commensal_indices = which(commensal_overlap)
        
        if(length(commensal_indices) > 0) {
          engulf_success = runif(length(commensal_indices)) < phagocyte_activity_engulf[i]
          indices_to_engulf = commensal_indices[engulf_success]
          
          if(length(indices_to_engulf) > 0) {
            phagocyte_commensals_engulfed[i] = phagocyte_commensals_engulfed[i] + length(indices_to_engulf)
            commensal_coords = commensal_coords[-indices_to_engulf, , drop = FALSE]
            
            phagocyte_bacteria_registry[i, ] = shift_insert_fast(phagocyte_bacteria_registry[i, ], 
                                                                 rep(1, length(indices_to_engulf)))
            
            phagocyte_phenotype_index = phagocyte_phenotype[i] + 1
            commensals_killed_by_Mac[phagocyte_phenotype_index] = commensals_killed_by_Mac[phagocyte_phenotype_index] + length(indices_to_engulf)
          }
        }
      }
    }
    
    # Treg activation & effector actions
    if(allow_tregs_to_do_their_job) {
      M1_phagocyte_indices = which(phagocyte_phenotype == 1)
      
      if(length(M1_phagocyte_indices) > 0) {
        for(i in M1_phagocyte_indices) {
          px = phagocyte_x[i]
          py = phagocyte_y[i]
          
          treg_distances_x = abs(treg_x - px)
          treg_distances_y = abs(treg_y - py)
          nearby_treg_indices = which(treg_distances_x <= treg_vicinity_effect & 
                                        treg_distances_y <= treg_vicinity_effect)
          
          if(length(nearby_treg_indices) > 0) {
            num_pat_antigens = phagocyte_pathogens_engulfed[i]
            num_com_antigens = phagocyte_commensals_engulfed[i]
            
            if((num_pat_antigens + num_com_antigens) > 0) {
              rat_com_pat = num_com_antigens / (num_pat_antigens + num_com_antigens)
              
              if(rat_com_pat > rat_com_pat_threshold) {
                treg_phenotype[nearby_treg_indices] = 1
                treg_activity_SAMPs_binary[nearby_treg_indices] = 1
                treg_active_age[nearby_treg_indices] = 1
                
                if(allow_tregs_to_suppress_cognate) {
                  phagocyte_phenotype[i] = 2
                  phagocyte_active_age[i] = 1
                  bacteria_count = sum(phagocyte_bacteria_registry[i, ])
                  phagocyte_activity_ROS[i] = activity_ROS_M2_baseline
                  phagocyte_activity_engulf[i] = activity_engulf_M2_baseline + 
                    activity_engulf_M2_step * bacteria_count
                }
              }
            }
          }
        }
      }
    }
    
    # Kill pathogens with ROS
    if(nrow(pathogen_coords) > 0) {
      pathogen_avg_ROS = numeric(nrow(pathogen_coords))
      
      for(i in 1:nrow(pathogen_coords)) {
        pathogen_avg_ROS[i] = get_8n_avg_signal_fast(pathogen_coords[i, "x"], 
                                                     pathogen_coords[i, "y"], 
                                                     act_radius_ROS, ROS)
      }
      
      pathogens_to_kill = which(pathogen_avg_ROS > th_ROS_microbe)
      
      if(length(pathogens_to_kill) > 0) {
        pathogen_coords = pathogen_coords[-pathogens_to_kill, , drop = FALSE]
        pathogens_killed_by_ROS = pathogens_killed_by_ROS + length(pathogens_to_kill)
      }
    }
    
    # Kill commensals with ROS
    if(nrow(commensal_coords) > 0) {
      commensal_avg_ROS = numeric(nrow(commensal_coords))
      
      for(i in 1:nrow(commensal_coords)) {
        commensal_avg_ROS[i] = get_8n_avg_signal_fast(commensal_coords[i, "x"], 
                                                      commensal_coords[i, "y"], 
                                                      act_radius_ROS, ROS)
      }
      
      commensals_to_kill = which(commensal_avg_ROS > th_ROS_microbe)
      
      if(length(commensals_to_kill) > 0) {
        commensal_coords = commensal_coords[-commensals_to_kill, , drop = FALSE]
        commensals_killed_by_ROS = commensals_killed_by_ROS + length(commensals_to_kill)
      }
    }
    
    # Calculate pathogen counts touching epithelium
    pathogen_epithelium_counts = rep(0, grid_size)
    if(nrow(pathogen_coords) > 0) {
      epithelium_pathogens = pathogen_coords[pathogen_coords[, "y"] == 1, , drop = FALSE]
      if(nrow(epithelium_pathogens) > 0) {
        pathogen_epithelium_counts = tabulate(epithelium_pathogens[, "x"], nbins = grid_size)
      }
    }
    
    # Calculate commensal counts touching epithelium
    commensal_epithelium_counts = rep(0, grid_size)
    if(nrow(commensal_coords) > 0) {
      epithelium_commensals = commensal_coords[commensal_coords[, "y"] == 1, , drop = FALSE]
      if(nrow(epithelium_commensals) > 0) {
        commensal_epithelium_counts = tabulate(epithelium_commensals[, "x"], nbins = grid_size)
      }
    }
    
    # Update epithelium
    for(i in 1:nrow(epithelium)) {
      px = epithelium$x[i]
      
      x_coordinates = pmax(1, pmin(grid_size, (px - act_radius_ROS):(px + act_radius_ROS)))
      ros_values = ROS[1, x_coordinates]
      mean_ros = mean(ros_values)
      
      count_pathogens = pathogen_epithelium_counts[px]
      epithelium$level_injury[i] = epithelium$level_injury[i] + round(log(count_pathogens + 1))
      
      count_commensals = commensal_epithelium_counts[px]
      epithelium$level_injury[i] = epithelium$level_injury[i] + round(log(count_commensals + 1))
      
      if(mean_ros > th_ROS_epith_recover) {
        epithelium$level_injury[i] = epithelium$level_injury[i] + 1
      }
      
      epithelium$level_injury[i] = min(epithelium$level_injury[i], max_level_injury)
      
      if(epithelium$level_injury[i] > 0 && runif(1) < epith_recovery_chance) {
        epithelium$level_injury[i] = max(0, epithelium$level_injury[i] - 1)
      }
    }
    
    # Save longitudinal data
    epithelium_longitudinal[t, ] = as.numeric(table(factor(epithelium$level_injury, levels = 0:5)))
    
    phagocyte_counts = c(
      sum(phagocyte_phenotype == 0),
      tabulate(phagocyte_active_age[phagocyte_phenotype == 1] + 1, 6),
      tabulate(phagocyte_active_age[phagocyte_phenotype == 2] + 1, 6)
    )
    macrophages_longitudinal[t, ] = phagocyte_counts
    
    microbes_longitudinal[t, ] = c(nrow(commensal_coords), nrow(pathogen_coords))
    tregs_longitudinal[t, ] = c(sum(treg_phenotype == 0), sum(treg_phenotype == 1))
    microbes_cumdeath_longitudinal[t, ] = c(commensals_killed_by_ROS, commensals_killed_by_Mac, 
                                            pathogens_killed_by_ROS, pathogens_killed_by_Mac)
  }
  
  # Prepare results
  longitudinal_df = data.frame(
    epithelium_longitudinal,
    macrophages_longitudinal,
    microbes_longitudinal,
    tregs_longitudinal,
    microbes_cumdeath_longitudinal
  )
  
  longitudinal_df$t = 1:t_max
  longitudinal_df$sterile = sterile
  longitudinal_df$allow_tregs_to_do_their_job = allow_tregs_to_do_their_job
  longitudinal_df$allow_tregs_to_suppress_cognate = allow_tregs_to_suppress_cognate
  longitudinal_df$randomize_tregs = randomize_tregs
  
  # Create videos
  pattern = paste0("^frame_seed_\\d+_STERILE_", sterile, "_TREGS_", allow_tregs_to_do_their_job, "_trnd_", randomize_tregs, "_\\d+\\.png$")
  png_files = list.files(dir_name, full.names = TRUE, pattern = pattern)
  
  if(length(png_files) > 0) {
    png_files = png_files[order(as.numeric(gsub(".*_(\\d+)\\.png$", "\\1", png_files)))]
    
    video_out = paste0(dir_name,"/simulation_sterile_", sterile, "_tregs_", allow_tregs_to_do_their_job, "_trnd_", randomize_tregs, ".mp4")

    av_encode_video(
      input = png_files,
      output = video_out,
      framerate = 5,
      vfilter = "scale=1000:-2",
      codec = "libx264"
    )
    
    gif_out = paste0(dir_name, "/simulation_sterile_", sterile, "_tregs_", allow_tregs_to_do_their_job, "_trnd_", randomize_tregs, ".gif")
    
    av_encode_video(
      input = png_files,
      output = gif_out,
      framerate = 5,
      vfilter = "scale=1000:-2,fps=2"
    )

  }
  
  return(longitudinal_df)
}

# Set up parallel processing
num_cores = detectCores() - 4
print(paste("Using", num_cores, "cores"))

cl = makeCluster(num_cores)

# Load libraries on workers
clusterEvalQ(cl, {
  library(dplyr)
  library(tidyr)
  library(cowplot)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(Rcpp)
  library(av)
})

# Set working directory on workers
clusterEvalQ(cl, {
  setwd(getwd())
})

# Export ALL variables to workers
clusterExport(cl, ls(all.names = TRUE))

# Set random seeds
clusterSetRNGStream(cl, seed_in)

# Run scenarios in parallel
start_time = Sys.time()
print("Starting parallel execution...")

results_list = parLapply(cl, 1:nrow(scenarios_df), run_scenario)

end_time = Sys.time()
print(paste("Parallel execution time:", end_time - start_time))

# Clean up cluster
stopCluster(cl)

# Combine results
longitudinal_df_keep = do.call(rbind, results_list)

print("All scenarios completed!")
print(paste("Final dataset dimensions:", paste(dim(longitudinal_df_keep), collapse = " x ")))

# Save results
write.csv(longitudinal_df_keep, "simulation_results.csv", row.names = FALSE)
print("Results saved to simulation_results.csv")