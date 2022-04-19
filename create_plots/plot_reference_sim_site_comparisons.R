# plot_reference_sim_site_comparisons.R
# April 2022

# coordinate plotting of all site-specific simulation results against associated reference dataset, across all included validation relationships



library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)
library(RColorBrewer)


####################################################################################################
##################### setup: filepaths, info on simulations run, source scripts ####################
####################################################################################################
simulation_coordinator_path = "/Users/moniqueam/Documents/malaria-model_validation/simulation_inputs/simulation_coordinator.csv"
base_script_plot_filepath = "/Users/moniqueam/Documents/malaria-model_validation/create_plots"
base_reference_filepath = "/Users/moniqueam/Documents/malaria-model_validation/reference_datasets"
simulation_output_filepath = "/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/simulation_output/s220418"
plot_output_filepath = "/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/simulation_output/_plots"


source(file.path(base_script_plot_filepath, 'helper_functions_par_dens.R'))
source(file.path(base_script_plot_filepath, 'helper_functions_age_inc.R'))


####################################################################################################
############################ read in data and create plots #########################################
####################################################################################################

coord_csv = read.csv(simulation_coordinator_path)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#                      age - parasite density                     #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
par_dens_sites = coord_csv$site[intersect(which(!is.na(coord_csv$site)), which(coord_csv$age_parasite_density==1))]

# determine which of the parasite-density sites have the relevant simulation output
available_sites = c()
for (ii in 1:length(par_dens_sites)){
  if (file.exists(paste0(simulation_output_filepath, '/', par_dens_sites[ii], '/parasite_densities_by_age_month.csv'))){
    available_sites = c(available_sites, par_dens_sites[ii])
  }
}

for (ss in 1:length(available_sites)){
  cur_site = available_sites[ss]
  sim_df = read.csv(paste0(simulation_output_filepath, '/', cur_site, '/parasite_densities_by_age_month.csv'))
  age_agg_sim_df = get_age_bin_averages(sim_df)
  
  filepath_ref = paste0(base_reference_filepath, '/', coord_csv$age_parasite_density_ref[which(coord_csv$site == cur_site)])
  ref_df = read.csv(filepath_ref)
  ref_df = ref_df[tolower(ref_df$Site) == tolower(cur_site),]
  
  gg_plots = plot_par_dens_ref_sim_comparison(age_agg_sim_df, ref_df)
  ggsave(filename=paste0(plot_output_filepath, '/site_compare_par_dens_age_', cur_site, '.png'), plot=gg_plots[[2]])
}




# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#                         age - incidence                         #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
age_inc_sites = coord_csv$site[intersect(which(!is.na(coord_csv$site)), which(coord_csv$age_incidence==1))]

# determine which of the age-incidence sites have the relevant simulation output
available_sites = c()
for (ii in 1:length(age_inc_sites)){
  if (file.exists(paste0(simulation_output_filepath, '/', age_inc_sites[ii], '/inc_prev_data_final.csv'))){
    available_sites = c(available_sites, age_inc_sites[ii])
  }
}

# aggregate all age-incidence simulation data into one dataframe and reference data into a second dataframe (for the relevant sites)
for (ss in 1:length(available_sites)){
  cur_site = available_sites[ss]
  sim_df_cur = read.csv(paste0(simulation_output_filepath, '/', cur_site, '/inc_prev_data_final.csv'))
  upper_ages = sort(unique(sim_df_cur$Age))
  sim_df_cur$mean_age = sapply(sim_df_cur$Age, get_mean_from_upper_age, upper_ages = upper_ages)
  sim_df_cur$p_detect_case = coord_csv$p_detect_case[which(coord_csv$site == cur_site)]
  
  filepath_ref = paste0(base_reference_filepath, '/', coord_csv$age_incidence_ref[which(coord_csv$site == cur_site)])
  ref_df_cur = read.csv(filepath_ref)
  ref_df_cur = ref_df_cur[which(tolower(ref_df_cur$Site) == tolower(cur_site)),]
  
  if(ss == 1){
    sim_df = sim_df_cur
    ref_df = ref_df_cur
  } else{
    sim_df = rbind(sim_df, sim_df_cur)
    ref_df = rbind(ref_df, ref_df_cur)
  }
}
gg_plot = plot_inc_ref_sim_comparison(sim_df, ref_df)
ggsave(filename=paste0(plot_output_filepath, '/site_compare_incidence_age.pdf'), plot=gg_plot)



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#                         age - prevalence                        #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

# TODO: finish age-prevalence section - NOT DONE!!

age_prev_sites = coord_csv$site[intersect(which(!is.na(coord_csv$site)), which(coord_csv$age_prevalence==1))]

# determine which of the age-incidence sites have the relevant simulation output
available_sites = c()
for (ii in 1:length(age_prev_sites)){
  if (file.exists(paste0(simulation_output_filepath, '/', age_prev_sites[ii], '/inc_prev_data_final.csv'))){
    available_sites = c(available_sites, age_prev_sites[ii])
  }
}

# aggregate all age-incidence simulation data into one dataframe and reference data into a second dataframe (for the relevant sites)
for (ss in 1:length(available_sites)){
  cur_site = available_sites[ss]
  sim_df_cur = read.csv(paste0(simulation_output_filepath, '/', cur_site, '/inc_prev_data_final.csv'))
  upper_ages = sort(unique(sim_df_cur$Age))
  sim_df_cur$mean_age = sapply(sim_df_cur$Age, get_mean_from_upper_age, upper_ages = upper_ages)
  sim_df_cur$p_detect_case = coord_csv$p_detect_case[which(coord_csv$site == cur_site)]
  
  filepath_ref = paste0(base_reference_filepath, '/', coord_csv$age_prevalence_ref[which(coord_csv$site == cur_site)])
  ref_df_cur = read.csv(filepath_ref)
  ref_df_cur = ref_df_cur[which(tolower(ref_df_cur$Site) == tolower(cur_site)),]
  if('agebin' %in% colnames(ref_df_cur)){
    upper_ages = sort(unique(ref_df_cur$agebin))
    ref_df_cur$mean_age = sapply(ref_df_cur$agebin, get_mean_from_upper_age, upper_ages = upper_ages)
  }
  ref_df_cur[,c("Site", "mean_age", 'month', 'START_YEAR', 'PR')]
  
  if(ss == 1){
    sim_df = sim_df_cur
    ref_df = ref_df_cur
  } else{
    sim_df = rbind(sim_df, sim_df_cur)
    ref_df = rbind(ref_df, ref_df_cur)
  }
}
gg_plot = plot_inc_ref_sim_comparison(sim_df, ref_df)
ggsave(filename=paste0(plot_output_filepath, '/site_compare_incidence_age.pdf'), plot=gg_plot)

