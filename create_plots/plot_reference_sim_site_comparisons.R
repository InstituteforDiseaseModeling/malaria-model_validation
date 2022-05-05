# plot_reference_sim_site_comparisons.R
# April 2022

# coordinate plotting of all site-specific simulation results against associated reference dataset, across all included validation relationships



library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)


####################################################################################################
##################### setup: filepaths, info on simulations run, source scripts ####################
####################################################################################################
simulation_coordinator_path = "/Users/moniqueam/Documents/malaria-model_validation/simulation_inputs/simulation_coordinator.csv"
base_script_plot_filepath = "/Users/moniqueam/Documents/malaria-model_validation/create_plots"
base_reference_filepath = "/Users/moniqueam/Documents/malaria-model_validation/reference_datasets"
simulation_output_filepath = "/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/simulation_output/s220425"
plot_output_filepath = paste0(simulation_output_filepath, "/_plots")


source(file.path(base_script_plot_filepath, 'helper_functions_par_dens.R'))
source(file.path(base_script_plot_filepath, 'helper_functions_age_inc_prev.R'))
source(file.path(base_script_plot_filepath, 'helper_functions_infection_duration.R'))
source(file.path(base_script_plot_filepath, 'helper_functions_infectiousness.R'))
source(file.path(base_script_plot_filepath, 'calculate_comparison_metrics.R'))


####################################################################################################
############################ read in data and create plots #########################################
####################################################################################################
coord_csv = read.csv(simulation_coordinator_path)
if(!dir.exists(plot_output_filepath)) dir.create(plot_output_filepath)



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

generate_parasite_density_outputs(available_sites, coord_csv, simulation_output_filepath, base_reference_filepath, plot_output_filepath)

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
ggsave(filename=paste0(plot_output_filepath, '/site_compare_incidence_age.png'), plot=gg_plot)

# additional quantitative comparisons and metrics
combined_df = prepare_inc_df(sim_df, ref_df)

correlation_output = corr_ref_sim_points(combined_df)
correlation_df = correlation_output[[2]]
slope_correlation_output = corr_ref_deriv_sim_points(combined_df)
slope_correlation_df = slope_correlation_output[[2]]
ggsave(filename=paste0(plot_output_filepath, '/scatter_regression_incidence_age.png'), plot=correlation_output[[1]])
ggsave(filename=paste0(plot_output_filepath, '/scatter_regression_slope_incidence_age.png'), plot=slope_correlation_output[[1]])

mean_diff_df = calc_mean_rel_diff(combined_df)
correlation_df$corr_slope = correlation_df$slope
correlation_df$corr_r_squared = correlation_df$r.squared
slope_correlation_df$derivative_corr_slope = slope_correlation_df$slope
slope_correlation_df$derivative_corr_r_squared = slope_correlation_df$r.squared
quantitative_comparison_df = merge(mean_diff_df, correlation_df[,c('Site', 'corr_slope','corr_r_squared')], all=TRUE)
quantitative_comparison_df = merge(mean_diff_df, slope_correlation_df[,c('Site', 'derivative_corr_slope','derivative_corr_r_squared')], all=TRUE)

write.csv(quantitative_comparison_df, paste0(plot_output_filepath, '/comparison_metric_table_incidence_age.csv'))



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#                         age - prevalence                        #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

age_prev_sites = coord_csv$site[intersect(which(!is.na(coord_csv$site)), which(coord_csv$age_prevalence==1))]

# determine which of the age-incidence sites have the relevant simulation output
available_sites = c()
for (ii in 1:length(age_prev_sites)){
  if (file.exists(paste0(simulation_output_filepath, '/', age_prev_sites[ii], '/prev_inc_by_age_month.csv'))){
    available_sites = c(available_sites, age_prev_sites[ii])
  }
}

# aggregate all age-incidence simulation data into one dataframe and reference data into a second dataframe (for the relevant sites)
for (ss in 1:length(available_sites)){
  cur_site = available_sites[ss]
  
  # read in and format reference data for this site
  filepath_ref = paste0(base_reference_filepath, '/', coord_csv$age_prevalence_ref[which(coord_csv$site == cur_site)])
  ref_df_cur = read.csv(filepath_ref)
  ref_df_cur = ref_df_cur[which(tolower(ref_df_cur$Site) == tolower(cur_site)),]
  if('agebin' %in% colnames(ref_df_cur)){
    upper_ages = sort(unique(ref_df_cur$agebin))
    ref_df_cur$mean_age = sapply(ref_df_cur$agebin, get_mean_from_upper_age, upper_ages = upper_ages)
  } else if (('PR_LAR' %in% colnames(ref_df_cur)) & ('PR_UAR' %in% colnames(ref_df_cur))){
    ref_df_cur$mean_age = (ref_df_cur$PR_LAR + ref_df_cur$PR_UAR)/2
  }
  colnames(ref_df_cur)[colnames(ref_df_cur) == 'PR_MONTH'] = 'month'
  colnames(ref_df_cur)[colnames(ref_df_cur) == 'PR'] = 'prevalence'
  colnames(ref_df_cur)[colnames(ref_df_cur) == 'N'] = 'total_sampled'
  colnames(ref_df_cur)[colnames(ref_df_cur) == 'N_POS'] = 'num_pos'
  colnames(ref_df_cur)[colnames(ref_df_cur) == 'PR_YEAR'] = 'year'
  ref_df_cur = ref_df_cur[,c("Site", "mean_age", 'month', 'total_sampled', 'num_pos', 'prevalence', 'year')]
  # remove reference rows without prevalence values
  ref_df_cur = ref_df_cur[!is.na(ref_df_cur$prevalence),]
  
  # read in and format data from simulations to match reference dataset
  sim_df_cur = read.csv(paste0(simulation_output_filepath, '/', cur_site, '/prev_inc_by_age_month.csv'))
  colnames(sim_df_cur)[colnames(sim_df_cur) == 'PfPR'] = 'prevalence'
  upper_ages = sort(unique(sim_df_cur$agebin))
  sim_df_cur$mean_age = sapply(sim_df_cur$agebin, get_mean_from_upper_age, upper_ages = upper_ages)
  # determine whether reference is for a single month or averaged over multiple months - subset from simulations to match
  if((length(unique(ref_df_cur$month))==1) & (is.character(ref_df_cur$month[1]))){  # check whether multiple months are listed in a character string
    included_months = as.numeric(unlist(strsplit(ref_df_cur$month[1],",")))
    sim_df_cur = sim_df_cur[sim_df_cur$month %in% included_months,]
    if(length(included_months)>1){
      sim_df_cur$month = 'multiple'
      ref_df_cur$month = 'multiple'
    }
  }else if (all(is.numeric(ref_df_cur$month))){
    included_months = unique(ref_df_cur$month)
    sim_df_cur = sim_df_cur[sim_df_cur$month %in% included_months,]
  } else{
    warning(paste0('The month format in the ', cur_site, ' reference dataset was not recognized.'))
  }
  sim_df_cur = sim_df_cur[,c("Site", "mean_age", 'month', 'prevalence', 'year', 'Run_Number')]
  
  # merge into dataset with all sites
  if(ss == 1){
    sim_df = sim_df_cur
    ref_df = ref_df_cur
  } else{
    sim_df = rbind(sim_df, sim_df_cur)
    ref_df = rbind(ref_df, ref_df_cur)
  }
}
gg_plot = plot_prev_ref_sim_comparison(sim_df, ref_df)
ggsave(filename=paste0(plot_output_filepath, '/site_compare_prevalence_age.png'), plot=gg_plot)






# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#                    age - infection duration                     #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# set positive threshold density for sampled parasites in simulation output (to match PCR threshold in reference)
pos_thresh_dens = 0.5  # (39=min(ref_df$DENSITY[ref_df$DENSITY>0], na.rm=TRUE) - 1)
# specify binning for duration of infection
duration_bins=c(seq(0,350,50), 500)


age_duration_sites = coord_csv$site[intersect(which(!is.na(coord_csv$site)), which(coord_csv$age_parasite_density==1))]

# determine which of the parasite-density sites have the relevant simulation output
available_sites = c()
for (ii in 1:length(age_duration_sites)){
  if (file.exists(paste0(simulation_output_filepath, '/', age_duration_sites[ii], '/patient_reports.csv'))){
    available_sites = c(available_sites, age_duration_sites[ii])
  }
}

for (ss in 1:length(available_sites)){
  cur_site = available_sites[ss]
  sim_dir = paste0(simulation_output_filepath, '/', cur_site)

  filepath_ref = paste0(base_reference_filepath, '/', coord_csv$infection_duration_ref[which(coord_csv$site == cur_site)])
  ref_df = read.csv(filepath_ref)
  ref_df = ref_df[tolower(ref_df$site) == tolower(cur_site),]
  
  gg_plots = plot_duration_ref_sim_comparison(sim_dir, ref_df)
  
  ggsave(filename=paste0(plot_output_filepath, '/site_compare_infect_duration_', cur_site, '.png'), plot=gg_plots[[1]])
  ggsave(filename=paste0(plot_output_filepath, '/site_compare_age_infect_duration_', cur_site, '.png'), plot=gg_plots[[2]])
  ggsave(filename=paste0(plot_output_filepath, '/site_compare_duration_measures_', cur_site, '.png'), plot=gg_plots[[3]])
}









# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#                   infectiouness to vectors                        #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

infectiousness_sites = coord_csv$site[intersect(which(!is.na(coord_csv$site)), which(coord_csv$infectiousness_to_mosquitos==1))]

# determine which of the parasite-density sites have the relevant simulation output
available_sites = c()
for (ii in 1:length(infectiousness_sites)){
  if (file.exists(paste0(simulation_output_filepath, '/', infectiousness_sites[ii], '/infectiousness_by_age_density_month.csv'))){
    available_sites = c(available_sites, infectiousness_sites[ii])
  }
}

for (ss in 1:length(available_sites)){
  cur_site = available_sites[ss]
  sim_df = read.csv(paste0(simulation_output_filepath, '/', cur_site, '/infectiousness_by_age_density_month.csv'))

  filepath_ref = paste0(base_reference_filepath, '/', coord_csv$infectiousness_to_mosquitos_ref[which(coord_csv$site == cur_site)])
  ref_df = read.csv(filepath_ref)
  ref_df = ref_df[tolower(ref_df$site) == tolower(cur_site),]
  
  gg_plot = plot_infectiousness_ref_sim_comparison(sim_df, ref_df)
  ggsave(filename=paste0(plot_output_filepath, '/site_compare_infectiousness_', cur_site, '.png'), plot=gg_plot)
}


