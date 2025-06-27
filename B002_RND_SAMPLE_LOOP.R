rm(list=ls())
library(dplyr)
library(tidyr)
library(Rcpp)

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

args_trailing = commandArgs(trailingOnly = TRUE)
if (length(args_trailing) >= 1) {
  core_index = as.integer(args_trailing[1])
} else {
  stop("Please provide a core_index as the first argument when running the script.")
}
set.seed(core_index)

source("./FAST_FUNCTIONS.R")

# start_time = proc.time()

dir_name = './random_simulations'
dir.create(dir_name)

##### NEED TO IMPLEMENT AGE FOR ALL CELLS! INCLUDING EPITHELIUM!
#### HOW TO ATTRACT T CELLS TO DCS?

scenarios_df = expand.grid(
  sterile = c(0, 1),
  allow_tregs_to_do_their_job = c(FALSE, TRUE),
  allow_tregs_to_suppress_cognate = c(FALSE),
  randomize_tregs = 0
)

t_max            = 2000
num_realizations = 1000 

th_ROS_microbe_vec        = runif(num_realizations, min=0, max=1)
th_ROS_epith_recover_vec  = runif(num_realizations, min=0, max=1)
epith_recovery_chance_vec = runif(num_realizations, min=0, max=1)
rat_com_pat_threshold_vec = runif(num_realizations, min=0.5, max=1) # at least above half?

#### THESE BOTH NEED TO BE BELOW 0.125! (BECAUSE 8 NEIGHBORS AND DIFFUSES **ALL** TO 8 NEIGHBORS NOTHING LEFT IN THE CENTER SINCE 0.125*8 =1)
diffusion_speed_DAMPs_vec = runif(num_realizations, min=0, max=0.12)
diffusion_speed_SAMPs_vec = runif(num_realizations, min=0, max=0.12)
diffusion_speed_ROS_vec   = runif(num_realizations, min=0, max=0.12)

add_ROS_vec   = runif(num_realizations, min=0, max=1)
add_DAMPs_vec = runif(num_realizations, min=0, max=1)
add_SAMPs_vec = runif(num_realizations, min=0, max=1)

ros_decay_vec   = runif(num_realizations, min=0, max=1)
DAMPs_decay_vec = runif(num_realizations, min=0, max=1)
SAMPs_decay_vec = runif(num_realizations, min=0, max=1)

activation_threshold_DAMPs_vec = runif(num_realizations, min=0, max=1)
activation_threshold_SAMPs_vec = runif(num_realizations, min=0, max=1)

# engulfing and ROS rates - baseline is max 50% of the capacity
activity_engulf_M0_baseline_vec    = runif(num_realizations, min=0, max=0.5)
activity_engulf_M1_baseline_vec    = runif(num_realizations, min=0, max=0.5)
activity_engulf_M2_baseline_vec    = runif(num_realizations, min=0, max=0.5)
activity_ROS_M1_baseline_vec       = runif(num_realizations, min=0, max=0.5)

rate_leak_commensal_injury_vec     = runif(num_realizations, min=0.5, max=1)
rate_leak_pathogen_injury_vec      = runif(num_realizations, min=0.5, max=1)
rate_leak_commensal_baseline_vec   = runif(num_realizations, min=0, max=0.25)

active_age_limit_vec  = sample(seq(3, 10, 1), num_realizations, replace=T)
# digestion_time_vec    = sample(seq(3, 10, 1), num_realizations, replace=T)

# Combine into dataframe
parameters_df = data.frame(
  th_ROS_microbe               = th_ROS_microbe_vec,
  th_ROS_epith_recover         = th_ROS_epith_recover_vec,
  epith_recovery_chance        = epith_recovery_chance_vec,
  rat_com_pat_threshold        = rat_com_pat_threshold_vec,
  diffusion_speed_DAMPs        = diffusion_speed_DAMPs_vec,
  diffusion_speed_SAMPs        = diffusion_speed_SAMPs_vec,
  diffusion_speed_ROS          = diffusion_speed_ROS_vec,
  add_ROS                      = add_ROS_vec,
  add_DAMPs                    = add_DAMPs_vec,
  add_SAMPs                    = add_SAMPs_vec,
  ros_decay                    = ros_decay_vec,
  DAMPs_decay                  = DAMPs_decay_vec,
  SAMPs_decay                  = SAMPs_decay_vec,
  activation_threshold_DAMPs   = activation_threshold_DAMPs_vec,
  activation_threshold_SAMPs   = activation_threshold_SAMPs_vec,
  activity_engulf_M0_baseline  = activity_engulf_M0_baseline_vec,
  activity_engulf_M1_baseline  = activity_engulf_M1_baseline_vec,
  activity_engulf_M2_baseline  = activity_engulf_M2_baseline_vec,
  activity_ROS_M1_baseline     = activity_ROS_M1_baseline_vec,
  rate_leak_commensal_injury   = rate_leak_commensal_injury_vec,
  rate_leak_pathogen_injury    = rate_leak_pathogen_injury_vec,
  active_age_limit             = active_age_limit_vec,
  # digestion_time               = digestion_time_vec,
  rate_leak_commensal_baseline = rate_leak_commensal_baseline_vec
)
#########################################
treg_vicinity_effect = 1 # if 0, that means has to be at the very same pixel
grid_size            = 25
n_phagocytes         = 200
n_tregs              = 200
injury_percentage    = 60
injury_site = get_middle_percent(seq(1,grid_size), injury_percentage)

max_cell_value_ROS   = 1
max_cell_value_DAMPs = 1
max_cell_value_SAMPs = 1

# for plotting purposes
lim_ROS    = max_cell_value_ROS
lim_DAMP   = max_cell_value_DAMPs
lim_SAMP   = max_cell_value_SAMPs

act_radius_ROS             = 1 # ROS RADIUS, SMALLER THE LOCAL
act_radius_treg            = 1 # treg radius
act_radius_DAMPs           = 1 # DAMP RADIUS, SMALLER THE LOCAL
act_radius_SAMPs           = 1 # SAMP RADIUS, SMALLER THE LOCAL

activity_engulf_max       = 0.99
activity_ROS_M0_baseline  = 0.00 
activity_ROS_M2_baseline  = 0.00 
activity_ROS_max          = 0.99

max_level_injury          = 5
cc_phagocyte              = 5
digestion_time            = 3

## logistic function
k_in  = 0.044
x0_in = 50

saveRDS(parameters_df,paste0(dir_name,'/parameters_df_',core_index,'.rds'))

for (realization_ind in 1:num_realizations){
  params                        = parameters_df[realization_ind, ]
  th_ROS_microbe                = params$th_ROS_microbe  # above this, will kill the microbe
  th_ROS_epith_recover          = params$th_ROS_epith_recover  # above this, will kill the epithelial cell
  epith_recovery_chance         = params$epith_recovery_chance 
  rat_com_pat_threshold         = params$rat_com_pat_threshold
  diffusion_speed_DAMPs         = params$diffusion_speed_DAMPs
  diffusion_speed_SAMPs         = params$diffusion_speed_SAMPs
  diffusion_speed_ROS           = params$diffusion_speed_ROS
  add_ROS                       = params$add_ROS
  add_DAMPs                     = params$add_DAMPs
  add_SAMPs                     = params$add_SAMPs
  ros_decay                     = params$ros_decay
  DAMPs_decay                   = params$DAMPs_decay
  SAMPs_decay                   = params$SAMPs_decay
  activation_threshold_DAMPs    = params$activation_threshold_DAMPs
  activation_threshold_SAMPs    = params$activation_threshold_SAMPs
  activity_engulf_M0_baseline   = params$activity_engulf_M0_baseline
  activity_engulf_M1_baseline   = params$activity_engulf_M1_baseline
  activity_engulf_M2_baseline   = params$activity_engulf_M2_baseline
  activity_ROS_M1_baseline      = params$activity_ROS_M1_baseline
  rate_leak_commensal_injury    = params$rate_leak_commensal_injury
  # digestion_time                = params$digestion_time
  active_age_limit              = params$active_age_limit
  rate_leak_commensal_baseline  = params$rate_leak_commensal_baseline
  
  for (scenario_ind in 1:nrow(scenarios_df)){
    sterile                         = scenarios_df[scenario_ind,]$sterile
    allow_tregs_to_do_their_job     = scenarios_df[scenario_ind,]$allow_tregs_to_do_their_job
    allow_tregs_to_suppress_cognate = scenarios_df[scenario_ind,]$allow_tregs_to_suppress_cognate
    randomize_tregs                 = scenarios_df[scenario_ind,]$randomize_tregs
    
    if(sterile==0){
      rate_leak_pathogen_injury     = params$rate_leak_pathogen_injury
    }else{
      rate_leak_pathogen_injury      = 0.00
    }
    
    n_commensals_lp = 20
    n_pathogens_lp  = round(rate_leak_pathogen_injury*length(injury_site)) 
    
    # Initialize fields
    DAMPs  = matrix(0, grid_size, grid_size)
    SAMPs  = matrix(0, grid_size, grid_size)
    ROS    = matrix(0, grid_size, grid_size)
    
    activity_engulf_M1_step   = (activity_engulf_max-activity_engulf_M1_baseline)/cc_phagocyte
    activity_engulf_M2_step   = (activity_engulf_max-activity_engulf_M2_baseline)/cc_phagocyte
    activity_ROS_M1_step      = (activity_ROS_max-activity_ROS_M1_baseline)/cc_phagocyte
    
    epithelium  = data.frame(x = seq(1,grid_size,1),
                             y = rep(0, grid_size),
                             level_injury = 0, id = seq(1,grid_size))
    
    epithelium[injury_site,]$level_injury = 1 # start with 1
    
    #phagocytes
    phagocyte_x = sample(1:grid_size, n_phagocytes, TRUE)
    phagocyte_y = sample(2:grid_size, n_phagocytes, TRUE) 
    phagocyte_pathogens_engulfed  = rep(0, n_phagocytes)
    phagocyte_commensals_engulfed = rep(0, n_phagocytes)
    phagocyte_num_times_activated = rep(0, n_phagocytes)
    phagocyte_phenotype       = rep(0, n_phagocytes)  # 0=M0, 1=M1, 2=M2
    phagocyte_activity_ROS    = rep(activity_ROS_M0_baseline, n_phagocytes)
    phagocyte_activity_engulf = rep(activity_engulf_M0_baseline, n_phagocytes)
    phagocyte_active_age      = rep(0, n_phagocytes)
    
    # For bacteria_registry, use a matrix instead of list of vectors
    phagocyte_bacteria_registry = matrix(0, nrow = n_phagocytes, ncol = cc_phagocyte)
    
    # tregs
    treg_x          = sample(1:grid_size, n_tregs, TRUE)
    treg_y          = sample(2:grid_size, n_tregs, TRUE)
    treg_active_age = rep(0, n_tregs)
    treg_phenotype  = rep(0, n_tregs)  # 0=resting, 1=activated
    treg_activity_SAMPs_binary = rep(0, n_tregs)
    
    if(n_pathogens_lp == 0) {
      pathogen_coords = matrix(numeric(0), ncol = 3)
      colnames(pathogen_coords) = c("x", "y", "id")
    } else {
      pathogen_coords = matrix(c(
        sample(injury_site, n_pathogens_lp, TRUE),
        rep(1, n_pathogens_lp),
        seq(1, n_pathogens_lp)
      ), ncol = 3)
      colnames(pathogen_coords) = c("x", "y", "id")
    }
    
    # Commensals - use matrix
    commensal_coords = matrix(c(
      sample(1:grid_size, n_commensals_lp, TRUE),
      sample(1:grid_size, n_commensals_lp, TRUE),
      seq(1, n_commensals_lp)
    ), ncol = 3)
    colnames(commensal_coords) = c("x", "y", "id")
    
    last_id_pathogen  = n_pathogens_lp
    last_id_commensal = n_commensals_lp
    
    p_prev_mic = NULL
    p_prev_lym = NULL
    
    pathogens_killed_by_ROS = 0
    pathogens_killed_by_Mac = rep(0,3)
    
    commensals_killed_by_ROS = 0
    commensals_killed_by_Mac = rep(0,3)
    
    ### KEEP THE NUMBERS
    epithelium_longitudinal  = matrix(0,nrow=t_max, ncol=(max_level_injury+1))
    macrophages_longitudinal = matrix(0,nrow=t_max, ncol=1+2*(cc_phagocyte+1)) # M0, M1, M2 levels
    microbes_longitudinal    = matrix(0,nrow=t_max, ncol=2) # first column commensals, second column pathogens
    tregs_longitudinal       = matrix(0,nrow=t_max, ncol=2) # two phenotypes, active and resting
    microbes_cumdeath_longitudinal = matrix(0,nrow=t_max, ncol=2*4) # first 4 columns commensals, second 4 columns pathogens - death by ROS, M0, M1, M2

    for (t in 1:t_max) {
      
      # Update injury site
      injury_site_updated = which(epithelium$level_injury>0)
      
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
      
      # Move pathogens and commensals randomly (optimized for matrices)
      if(nrow(pathogen_coords) > 0) {
        # Initialize dy vector
        dy = ifelse(pathogen_coords[, "y"] == 1,
                    sample(c(1), size = nrow(pathogen_coords), replace = TRUE),
                    sample(c(-1, 0, 1), size = nrow(pathogen_coords), replace = TRUE))
        
        # dx sampled conditionally to avoid (0,0)
        dx = iszero_coordinates(dy)
        
        # Update positions
        pathogen_coords[, "x"] = pmin(pmax(pathogen_coords[, "x"] + dx, 1), grid_size)
        pathogen_coords[, "y"] = pmin(pmax(pathogen_coords[, "y"] + dy, 1), grid_size)
      }
      if(nrow(commensal_coords) > 0) {
        # Initialize dy vector
        dy = ifelse(commensal_coords[, "y"] == 1,
                    sample(c(1), size = nrow(commensal_coords), replace = TRUE),
                    sample(c(-1, 0, 1), size = nrow(commensal_coords), replace = TRUE))
        
        # dx sampled conditionally to avoid (0,0)
        dx = iszero_coordinates(dy)
        
        # Update positions
        commensal_coords[, "x"] = pmin(pmax(commensal_coords[, "x"] + dx, 1), grid_size)
        commensal_coords[, "y"] = pmin(pmax(commensal_coords[, "y"] + dy, 1), grid_size)   
      }
      
      # Pre-calculate pathogen counts touching epithelium (y=1) for each x position
      pathogen_epithelium_counts = rep(0, grid_size)
      if(nrow(pathogen_coords) > 0) {
        epithelium_pathogens = pathogen_coords[pathogen_coords[, "y"] == 1, , drop = FALSE]
        if(nrow(epithelium_pathogens) > 0) {
          pathogen_epithelium_counts = tabulate(epithelium_pathogens[, "x"], nbins = grid_size)
        }
      }
      
      # Pre-calculate commensal counts touching epithelium (y=1) for each x position
      commensal_epithelium_counts = rep(0, grid_size)
      if(nrow(commensal_coords) > 0) {
        epithelium_commensals = commensal_coords[commensal_coords[, "y"] == 1, , drop = FALSE]
        if(nrow(epithelium_commensals) > 0) {
          commensal_epithelium_counts = tabulate(epithelium_commensals[, "x"], nbins = grid_size)
        }
      }
      
      # Update DAMPs
      DAMPs[1,] = DAMPs[1,] + epithelium$level_injury*add_DAMPs
      DAMPs[1,] = DAMPs[1,] + 1*logistic_scaled_0_to_5_quantized(pathogen_epithelium_counts+commensal_epithelium_counts)*add_DAMPs
      
      # Diffuse & decay DAMPs, SAMPs, and ROS.
      DAMPs   = diffuse_matrix(DAMPs, diffusion_speed_DAMPs, max_cell_value_DAMPs)
      SAMPs   = diffuse_matrix(SAMPs, diffusion_speed_SAMPs, max_cell_value_SAMPs)
      ROS     = diffuse_matrix(ROS, diffusion_speed_ROS, max_cell_value_ROS)
      DAMPs   = DAMPs - DAMPs_decay*DAMPs
      SAMPs   = SAMPs - SAMPs_decay*SAMPs
      ROS     = ROS - ros_decay*ROS
      
      # Move phagocytes and tregs based on DAMPs gradient
      density_matrix = DAMPs # here you can choose another density matrix if you like
      all_equal      = all(density_matrix == density_matrix[1, 1])
      
      if(!all_equal){ 
        # Move tregs
        ##### move_agents_vectors_gradient  ######################################## 
        for(i in 1:length(treg_x)) {
          x = treg_x[i]
          y = treg_y[i]
          
          # Get 3x3 neighborhood (with boundary check) - same as original
          x_range = max(1, x - 1):min(grid_size, x + 1)
          y_range = max(1, y - 1):min(grid_size, y + 1)
          
          # Create all neighbor combinations
          neighbors_x = rep(x_range, each = length(y_range))
          neighbors_y = rep(y_range, times = length(x_range))
          
          # Get density values for those cells (vectorized)
          neighbor_densities = density_matrix[cbind(neighbors_y, neighbors_x)]
          
          # Normalize to get probabilities
          total = sum(neighbor_densities)
          if(total > 0) {
            probs = neighbor_densities / total
          } else {
            probs = rep(1 / length(neighbor_densities), length(neighbor_densities))
          }
          
          # Sample a move weighted by local density
          chosen_idx = sample(1:length(neighbors_x), 1, prob = probs)
          
          # Update position (modifies in place using)
          treg_x[i] = neighbors_x[chosen_idx]
          treg_y[i] = neighbors_y[chosen_idx]
        }
        ###############################################################################
        
        # Move all phagocytes
        for(i in 1:length(phagocyte_x)){
          x = phagocyte_x[i]
          y = phagocyte_y[i]
          
          # Get 3x3 neighborhood (with boundary check) - same as original
          x_range = max(1, x - 1):min(grid_size, x + 1)
          y_range = max(1, y - 1):min(grid_size, y + 1)
          
          # Create all neighbor combinations
          neighbors_x = rep(x_range, each = length(y_range))
          neighbors_y = rep(y_range, times = length(x_range))
          
          # Get density values for those cells (vectorized)
          neighbor_densities = density_matrix[cbind(neighbors_y, neighbors_x)]
          
          # Normalize to get probabilities
          total = sum(neighbor_densities)
          if(total > 0) {
            probs = neighbor_densities / total
          } else {
            probs = rep(1 / length(neighbor_densities), length(neighbor_densities))
          }
          
          # Sample a move weighted by local density
          chosen_idx = sample(1:length(neighbors_x), 1, prob = probs)
          
          # Update position
          phagocyte_x[i] = neighbors_x[chosen_idx]
          phagocyte_y[i] = neighbors_y[chosen_idx]
        }
        
      } else {
        # Random movement when no gradient 
        # tregs
        dy_treg = ifelse(treg_y == 1,
                         sample(c(1), size = length(treg_y), replace = TRUE),
                         sample(c(-1, 0, 1), size = length(treg_y), replace = TRUE))
        dx_treg = iszero_coordinates(dy_treg)# dx sampled conditionally to avoid (0,0)
        # Update positions with boundary constraints 
        treg_x = pmin(pmax(treg_x + dx_treg, 1), grid_size)
        treg_y = pmin(pmax(treg_y + dy_treg, 1), grid_size)
        
        # All phagocytes
        dy_phagocyte = ifelse(phagocyte_y == 1,
                              sample(c(1), size = length(phagocyte_y), replace = TRUE),
                              sample(c(-1, 0, 1), size = length(phagocyte_y), replace = TRUE))
        dx_phagocyte = iszero_coordinates(dy_phagocyte)# dx sampled conditionally to avoid (0,0)
        # Update positions with boundary constraints 
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
        # Baseline commensals
        baseline_x = sample(1:grid_size, n_commensals_lp_new_baseline, TRUE)
        # Injury-site commensals  
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

      # Update the phagocyte phenotypes
      # Vectorized operations where possible
      M0_indices = which(phagocyte_phenotype == 0)
      M1_indices = which(phagocyte_phenotype == 1)
      M2_indices = which(phagocyte_phenotype == 2)
      
      # Registry shifting every digestion_time steps (vectorized)
      if(t %% digestion_time == 0) {
        # Shift all registries at once using matrix operations
        phagocyte_bacteria_registry = cbind(
          matrix(0, nrow = nrow(phagocyte_bacteria_registry), ncol = 1),
          phagocyte_bacteria_registry[, -ncol(phagocyte_bacteria_registry)]
        )
      }
      
      # Process M0 phagocytes (candidates for activation)
      if(length(M0_indices) > 0) {
        for(i in M0_indices) {
          # Get signals (this function call is unavoidable but now faster)
          avg_DAMPs = get_8n_avg_signal_fast(phagocyte_x[i], phagocyte_y[i], act_radius_DAMPs, DAMPs)
          avg_SAMPs = get_8n_avg_signal_fast(phagocyte_x[i], phagocyte_y[i], act_radius_SAMPs, SAMPs)
          
          # Activation logic
          if(avg_DAMPs >= activation_threshold_DAMPs && avg_DAMPs > avg_SAMPs) {
            phagocyte_phenotype[i] = 1  # M1
            phagocyte_active_age[i] = 1
            bacteria_count = sum(phagocyte_bacteria_registry[i, ])
            phagocyte_activity_ROS[i] = activity_ROS_M1_baseline + activity_ROS_M1_step * bacteria_count
            phagocyte_activity_engulf[i] = activity_engulf_M1_baseline + activity_engulf_M1_step * bacteria_count
          } else if(avg_SAMPs >= activation_threshold_SAMPs && avg_SAMPs > avg_DAMPs) {
            phagocyte_phenotype[i] = 2  # M2
            phagocyte_active_age[i] = 1
            bacteria_count = sum(phagocyte_bacteria_registry[i, ])
            phagocyte_activity_ROS[i] = activity_ROS_M2_baseline
            phagocyte_activity_engulf[i] = activity_engulf_M2_baseline + activity_engulf_M2_step * bacteria_count
          }
        }
      }
      
      # Process M1/M2 phagocytes
      active_indices = c(M1_indices, M2_indices)
      if(length(active_indices) > 0) {
        # Vectorized age increment
        phagocyte_active_age[active_indices] = phagocyte_active_age[active_indices] + 1
        
        # Check which ones are old enough to potentially change
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
            phagocyte_phenotype[i] = 0  # Back to M0
            phagocyte_active_age[i] = 0
            phagocyte_activity_ROS[i] = activity_ROS_M0_baseline
            phagocyte_activity_engulf[i] = activity_engulf_M0_baseline
          }
        }
      }
      
      # Treg active age updated
      active_treg_indices = which(treg_phenotype == 1)
      if(length(active_treg_indices) > 0) {
        # Vectorized age check
        old_tregs = active_treg_indices[treg_active_age[active_treg_indices] >= active_age_limit]
        young_tregs = active_treg_indices[treg_active_age[active_treg_indices] < active_age_limit]
        
        # Age the young ones
        if(length(young_tregs) > 0) {
          treg_active_age[young_tregs] = treg_active_age[young_tregs] + 1
        }
        
        # Shut down the old ones
        if(length(old_tregs) > 0) {
          treg_phenotype[old_tregs] = 0
          treg_active_age[old_tregs] = 0
          treg_activity_SAMPs_binary[old_tregs] = 0
        }
      }
      
      # Engulfment process
      # Pre-calculate all phagocyte positions for faster lookup
      phagocyte_positions = paste(phagocyte_x, phagocyte_y, sep = "_")
      
      for(i in 1:length(phagocyte_x)) {
        px = phagocyte_x[i]
        py = phagocyte_y[i]
        
        # Fast pathogen overlap check using vectorized operations
        if(nrow(pathogen_coords) > 0) {
          pathogen_overlap = (pathogen_coords[, "x"] == px) & (pathogen_coords[, "y"] == py)
          pathogen_indices = which(pathogen_overlap)
          
          if(length(pathogen_indices) > 0) {
            # Stochastic engulfment
            engulf_success = runif(length(pathogen_indices)) < phagocyte_activity_engulf[i]
            indices_to_engulf = pathogen_indices[engulf_success]
            
            if(length(indices_to_engulf) > 0) {
              phagocyte_pathogens_engulfed[i] = phagocyte_pathogens_engulfed[i] + length(indices_to_engulf)
              
              # Remove engulfed pathogens (create new matrix without these rows)
              pathogen_coords = pathogen_coords[-indices_to_engulf, , drop = FALSE]
              
              # Update bacteria registry efficiently
              phagocyte_bacteria_registry[i, ] = shift_insert_fast(phagocyte_bacteria_registry[i, ], 
                                                                   rep(1, length(indices_to_engulf)))
              
              # Update kill count
              phagocyte_phenotype_index = phagocyte_phenotype[i] + 1
              pathogens_killed_by_Mac[phagocyte_phenotype_index] = pathogens_killed_by_Mac[phagocyte_phenotype_index] + length(indices_to_engulf)
            }
          }
        }
        
        # Similar logic for commensals
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
        # Get all M1 phagocytes for treg interactions
        M1_phagocyte_indices = which(phagocyte_phenotype == 1)
        
        if(length(M1_phagocyte_indices) > 0) {
          for(i in M1_phagocyte_indices) {
            px = phagocyte_x[i]
            py = phagocyte_y[i]
            
            # Check if a Treg is in the vicinity (vectorized)
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
                  # TREG ACTIVATION (vectorized for all nearby tregs)
                  treg_phenotype[nearby_treg_indices] = 1
                  treg_activity_SAMPs_binary[nearby_treg_indices] = 1
                  treg_active_age[nearby_treg_indices] = 1
                  
                  if(allow_tregs_to_suppress_cognate) {
                    # Suppress this M1 phagocyte
                    phagocyte_phenotype[i] = 2  # Convert to M2
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
        # Vectorized ROS exposure calculation
        pathogen_avg_ROS = numeric(nrow(pathogen_coords))
        
        for(i in 1:nrow(pathogen_coords)) {
          pathogen_avg_ROS[i] = get_8n_avg_signal_fast(pathogen_coords[i, "x"], 
                                                       pathogen_coords[i, "y"], 
                                                       act_radius_ROS, ROS)
        }
        
        # Vectorized death determination
        pathogens_to_kill = which(pathogen_avg_ROS > th_ROS_microbe)
        
        if(length(pathogens_to_kill) > 0) {
          pathogen_coords = pathogen_coords[-pathogens_to_kill, , drop = FALSE]
          pathogens_killed_by_ROS = pathogens_killed_by_ROS + length(pathogens_to_kill)
        }
      }
      
      # Kill commensals with ROS
      if(nrow(commensal_coords) > 0) {
        # Vectorized ROS exposure calculation
        commensal_avg_ROS = numeric(nrow(commensal_coords))
        
        for(i in 1:nrow(commensal_coords)) {
          commensal_avg_ROS[i] = get_8n_avg_signal_fast(commensal_coords[i, "x"], 
                                                        commensal_coords[i, "y"], 
                                                        act_radius_ROS, ROS)
        }
        
        # Vectorized death determination
        commensals_to_kill = which(commensal_avg_ROS > th_ROS_microbe)
        
        if(length(commensals_to_kill) > 0) {
          commensal_coords = commensal_coords[-commensals_to_kill, , drop = FALSE]
          commensals_killed_by_ROS = commensals_killed_by_ROS + length(commensals_to_kill)
        }
      }
      
      # Injur epithelium with pathogens / commensals
      for(i in 1:nrow(epithelium)) {
        px = epithelium$x[i]
        
        # Get ROS values in vicinity
        x_coordinates = pmax(1, pmin(grid_size, (px - act_radius_ROS):(px + act_radius_ROS)))
        ros_values = ROS[1, x_coordinates]  # Row right below epithelium
        mean_ros = mean(ros_values)
        
        # Increase level_injury based on pathogen count
        count_pathogens = pathogen_epithelium_counts[px]
        # epithelium$level_injury[i] = epithelium$level_injury[i] + round(log(count_pathogens + 1))
        epithelium$level_injury[i] = epithelium$level_injury[i] + logistic_scaled_0_to_5_quantized(count_pathogens)
        
        # RULE 2: Increase level_injury based on ROS
        if(mean_ros > th_ROS_epith_recover) {
          epithelium$level_injury[i] = epithelium$level_injury[i] + 1
        }
        
        # Apply maximum injury constraint
        epithelium$level_injury[i] = min(epithelium$level_injury[i], max_level_injury)
        
        # RECOVERY: Stochastic recovery when injured
        if(epithelium$level_injury[i] > 0 && runif(1) < epith_recovery_chance) {
          epithelium$level_injury[i] = max(0, epithelium$level_injury[i] - 1)
        }
      }
      
      # Save abundances
      epithelium_longitudinal[t, ] = as.numeric(table(factor(epithelium$level_injury, levels = 0:5)))
      
      # Phenotype counting
      phagocyte_counts = c(
        sum(phagocyte_phenotype == 0),  # M0
        tabulate(phagocyte_active_age[phagocyte_phenotype == 1] + 1, 6),  # M1 by level
        tabulate(phagocyte_active_age[phagocyte_phenotype == 2] + 1, 6)   # M2 by level  
      )
      macrophages_longitudinal[t, ] = phagocyte_counts
      
      microbes_longitudinal[t, ] = c(nrow(commensal_coords), nrow(pathogen_coords))
      tregs_longitudinal[t, ] = c(sum(treg_phenotype == 0), sum(treg_phenotype == 1))
      microbes_cumdeath_longitudinal[t, ] = c(commensals_killed_by_ROS, commensals_killed_by_Mac, 
                                              pathogens_killed_by_ROS, pathogens_killed_by_Mac)
    }
    
    # Combine all matrices into one
    longitudinal_df = data.frame(
      epithelium_longitudinal,
      macrophages_longitudinal,
      microbes_longitudinal,
      tregs_longitudinal,
      microbes_cumdeath_longitudinal
    )
    
    # Add time column
    longitudinal_df$t = 1:t_max
    
    # Add metadata columns (same value repeated for each row)
    longitudinal_df$realization_ind               = realization_ind
    longitudinal_df$sterile                       = sterile
    longitudinal_df$allow_tregs_to_do_their_job   = allow_tregs_to_do_their_job
    longitudinal_df$allow_tregs_to_suppress_cognate = allow_tregs_to_suppress_cognate
    longitudinal_df$randomize_tregs                 = randomize_tregs
    
    # Optional: Reorder columns to have metadata first
    longitudinal_df = longitudinal_df %>%
      select(t, realization_ind, sterile, allow_tregs_to_do_their_job, allow_tregs_to_suppress_cognate, randomize_tregs, everything())
    
    saveRDS(longitudinal_df, paste0(dir_name,'/dataframe_',realization_ind,'_',scenario_ind,'_',core_index,'.rds'))
    print(paste0(dir_name,'/dataframe_',realization_ind,'_',scenario_ind,'_',core_index,'.rds SAVED'))
  }
}

