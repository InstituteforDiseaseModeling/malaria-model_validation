# helpers_coordinate_each_relationship.R

# For each validation relationship, this script contains a function that coordinates:
#     1) the reformatting and alignment of simulation and reference datasets into data frames that are used downstream, 
#     2) the quantitative comparisons between reference and new simulation results and between benchmark and new simulation results, and 
#     3) the generation of plotting outputs. 
# The main function for each relationship saves all plots and csvs to the output directory for the report-generating script to use. 


library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)
library(ggpubr)


################################################################################# 
########################### Incidence by age  ##################################
################################################################################# 

generate_age_incidence_outputs = function(coord_csv, simulation_output_filepath, base_reference_filepath, plot_output_filepath, benchmark_simulation_filepath=NA){

  # determine which of the age-incidence sites have the relevant simulation output
  age_inc_sites = coord_csv$site[intersect(which(!is.na(coord_csv$site)), which(coord_csv$age_incidence==1))]
  available_sites = c()
  for (ii in 1:length(age_inc_sites)){
    if (file.exists(paste0(simulation_output_filepath, '/', age_inc_sites[ii], '/inc_prev_data_final.csv'))){
      available_sites = c(available_sites, age_inc_sites[ii])
    }
  }
  
  # aggregate all age-incidence simulation data into one dataframe and reference data into a second dataframe (for the relevant sites)
  #   if a benchmark simulation directory was specified and the site exists there, create a third dataframe
  sim_df = data.frame()
  bench_df = data.frame()
  ref_df = data.frame()
  for (ss in 1:length(available_sites)){
    # simulations currently being evaluated
    cur_site = available_sites[ss]
    sim_df_cur = read.csv(paste0(simulation_output_filepath, '/', cur_site, '/inc_prev_data_final.csv'))
    upper_ages = sort(unique(sim_df_cur$Age))
    sim_df_cur$mean_age = sapply(sim_df_cur$Age, get_mean_from_upper_age, upper_ages = upper_ages)
    sim_df_cur$p_detect_case = coord_csv$p_detect_case[which(coord_csv$site == cur_site)]
    
    # simulations used as benchmark
    if(!is.na(benchmark_simulation_filepath) & file.exists(paste0(benchmark_simulation_filepath, '/', cur_site, '/inc_prev_data_final.csv'))){
      bench_df_cur = read.csv(paste0(benchmark_simulation_filepath, '/', cur_site, '/inc_prev_data_final.csv'))
      upper_ages = sort(unique(bench_df_cur$Age))
      bench_df_cur$mean_age = sapply(bench_df_cur$Age, get_mean_from_upper_age, upper_ages = upper_ages)
      bench_df_cur$p_detect_case = coord_csv$p_detect_case[which(coord_csv$site == cur_site)]
    } else{
      bench_df_cur = data.frame()
    }
    
    # reference data
    filepath_ref = paste0(base_reference_filepath, '/', coord_csv$age_incidence_ref[which(coord_csv$site == cur_site)])
    ref_df_cur = read.csv(filepath_ref)
    ref_df_cur = ref_df_cur[which(tolower(ref_df_cur$Site) == tolower(cur_site)),]
    
    # add site into larger dataframe
    sim_df = rbind(sim_df, sim_df_cur)
    ref_df = rbind(ref_df, ref_df_cur)
    bench_df = rbind(bench_df, bench_df_cur)
  }
  gg_plot = plot_inc_ref_sim_comparison(sim_df, ref_df, bench_df)
  ggsave(filename=paste0(plot_output_filepath, '/site_compare_incidence_age.png'), plot=gg_plot, height=2*ceiling(length(available_sites)/4), width=7.5, units='in')
  
  
  # additional quantitative comparisons and metrics between simulation and reference data
  combined_df = prepare_inc_df(sim_df, ref_df, bench_df)
  
  # correlations between new simulation and reference dataset values
  correlation_output = corr_ref_sim_points(combined_df)
  correlation_df = correlation_output[[2]]
  slope_correlation_output = corr_ref_deriv_sim_points(combined_df)
  slope_correlation_df = slope_correlation_output[[2]]
  combined_df_with_slopes = slope_correlation_output[[3]]
  correlation_plots = ggarrange(correlation_output[[1]], slope_correlation_output[[1]], nrow=1, ncol=2, 
                                common.legend=TRUE)  #, legend.grob=get_legend(correlation_output[[1]], position = 'bottom'))
  ggsave(filename=paste0(plot_output_filepath, '/scatter_regression_incidence_age.png'), plot=correlation_plots, height=4.5, width=8, units='in')
  
  # metrics comparing simulation to reference VALUE
  mean_diff_df = calc_mean_rel_diff(combined_df)
  correlation_df$corr_slope = correlation_df$slope
  correlation_df$corr_r_squared = correlation_df$r.squared
  # metrics comparing simulation to reference SLOPE
  mean_slope_diff_df = calc_mean_rel_slope_diff(combined_df=combined_df_with_slopes)
  slope_correlation_df$derivative_corr_slope = slope_correlation_df$slope
  slope_correlation_df$derivative_corr_r_squared = slope_correlation_df$r.squared
  # combine metrics and save csv
  quantitative_comparison_df = merge(mean_diff_df, correlation_df[,c('Site', 'corr_slope','corr_r_squared')], all=TRUE)
  # quantitative_comparison_slope_df = merge(mean_slope_diff_df, slope_correlation_df[,c('Site', 'derivative_corr_r_squared')], all=TRUE)
  # quantitative_comparison_df = merge(quantitative_comparison_df, quantitative_comparison_slope_df, all=TRUE)
  quantitative_comparison_df = merge(quantitative_comparison_df, mean_slope_diff_df, all=TRUE)
  write.csv(quantitative_comparison_df, paste0(plot_output_filepath, '/comparison_metric_table_incidence_age.csv'), row.names=FALSE)
  
  
  # compare simulation and benchmark simulation results
  # TODO: add likelihood component
  if('benchmark_value' %in% colnames(combined_df)){
    compare_benchmarks_output = compare_benchmark(combined_df)
    ggsave(filename=paste0(plot_output_filepath, '/scatter_benchmark_incidence_age.png'), plot=compare_benchmarks_output, height=4.5, width=8, units='in')
  }
}



################################################################################# 
########################### Prevalence by age  ##################################
################################################################################# 

generate_age_prevalence_outputs = function(coord_csv, simulation_output_filepath, base_reference_filepath, plot_output_filepath, benchmark_simulation_filepath=NA){
  
  age_prev_sites = coord_csv$site[intersect(which(!is.na(coord_csv$site)), which(coord_csv$age_prevalence==1))]
  
  # determine which of the age-incidence sites have the relevant simulation output
  available_sites = c()
  for (ii in 1:length(age_prev_sites)){
    if (file.exists(paste0(simulation_output_filepath, '/', age_prev_sites[ii], '/prev_inc_by_age_month.csv'))){
      available_sites = c(available_sites, age_prev_sites[ii])
    }
  }
  
  # aggregate all age-prevalence simulation data into one dataframe and reference data into a second dataframe (for the relevant sites)
  #   if a benchmark simulation directory was specified and the site exists there, create a third dataframe
  sim_df = data.frame()
  bench_df = data.frame()
  ref_df = data.frame()
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
    # remove rows with population = 0 (this is relevant for cohort simulations where there is only one age group in each year and all other age groups are zero)
    sim_df_cur = sim_df_cur[sim_df_cur$Pop >0,]
    
    if(!is.na(benchmark_simulation_filepath)){
      # read in and format data from benchmark simulations to match reference dataset
      bench_df_cur = read.csv(paste0(benchmark_simulation_filepath, '/', cur_site, '/prev_inc_by_age_month.csv'))
      colnames(bench_df_cur)[colnames(bench_df_cur) == 'PfPR'] = 'prevalence'
      upper_ages = sort(unique(bench_df_cur$agebin))
      bench_df_cur$mean_age = sapply(bench_df_cur$agebin, get_mean_from_upper_age, upper_ages = upper_ages)
      # remove rows with population = 0 (this is relevant for cohort simulations where there is only one age group in each year and all other age groups are zero)
      bench_df_cur = bench_df_cur[bench_df_cur$Pop >0,]
    } else{
      bench_df_cur = data.frame()
    }
    
    # determine whether reference is for a single month or averaged over multiple months - subset from simulations to match
    if((length(unique(ref_df_cur$month))==1) & (is.character(ref_df_cur$month[1]))){  # check whether multiple months are listed in a character string
      included_months = as.numeric(unlist(strsplit(ref_df_cur$month[1],",")))
      sim_df_cur = sim_df_cur[sim_df_cur$month %in% included_months,]
      bench_df_cur = bench_df_cur[bench_df_cur$month %in% included_months,]
      if(length(included_months)>1){
        sim_df_cur$month = 'multiple'
        ref_df_cur$month = 'multiple'
        if(nrow(bench_df_cur)>0) bench_df_cur$month = 'multiple'
      }
    }else if (all(is.numeric(ref_df_cur$month))){
      included_months = unique(ref_df_cur$month)
      sim_df_cur = sim_df_cur[sim_df_cur$month %in% included_months,]
      bench_df_cur = bench_df_cur[bench_df_cur$month %in% included_months,]
    } else{
      warning(paste0('The month format in the ', cur_site, ' reference dataset was not recognized.'))
    }
    sim_df_cur = sim_df_cur[,c("Site", "mean_age", 'month', 'prevalence', 'year', 'Run_Number')]
    if(nrow(bench_df_cur)>0) bench_df_cur = bench_df_cur[,c("Site", "mean_age", 'month', 'prevalence', 'year', 'Run_Number')]
    
    # add site into larger dataframe
    sim_df = rbind(sim_df, sim_df_cur)
    ref_df = rbind(ref_df, ref_df_cur)
    bench_df = rbind(bench_df, bench_df_cur)
  }
  gg_plot = plot_prev_ref_sim_comparison(sim_df, ref_df, bench_df)
  ggsave(filename=paste0(plot_output_filepath, '/site_compare_prevalence_age.png'), plot=gg_plot, height=9, width=8, units='in')
  
  
  # additional quantitative comparisons and metrics
  combined_df = prepare_prev_df(sim_df, ref_df, bench_df)
  
  # correlations between new simulation and reference dataset values
  correlation_output = corr_ref_sim_points(combined_df)
  correlation_df = correlation_output[[2]]
  slope_correlation_output = corr_ref_deriv_sim_points(combined_df)
  slope_correlation_df = slope_correlation_output[[2]]
  combined_df_with_slopes = slope_correlation_output[[3]]
  correlation_plots = ggarrange(correlation_output[[1]], slope_correlation_output[[1]], nrow=1, ncol=2, 
                                common.legend=TRUE)  #, legend.grob=get_legend(correlation_output[[1]], position = 'bottom'))
  ggsave(filename=paste0(plot_output_filepath, '/scatter_regression_prevalence_age.png'), plot=correlation_plots, height=4.5, width=8, units='in')
  
  # metrics comparing simulation to reference VALUE
  mean_diff_df = calc_mean_rel_diff(combined_df)
  correlation_df$corr_slope = correlation_df$slope
  correlation_df$corr_r_squared = correlation_df$r.squared
  # metrics comparing simulation to reference SLOPE
  mean_slope_diff_df = calc_mean_rel_slope_diff(combined_df=combined_df_with_slopes)
  slope_correlation_df$derivative_corr_slope = slope_correlation_df$slope
  slope_correlation_df$derivative_corr_r_squared = slope_correlation_df$r.squared
  # combine metrics and save csv
  quantitative_comparison_df = merge(mean_diff_df, correlation_df[,c('Site', 'corr_slope','corr_r_squared')], all=TRUE)
  # quantitative_comparison_slope_df = merge(mean_slope_diff_df, slope_correlation_df[,c('Site', 'derivative_corr_r_squared')], all=TRUE)
  # quantitative_comparison_df = merge(quantitative_comparison_df, quantitative_comparison_slope_df, all=TRUE)
  quantitative_comparison_df = merge(quantitative_comparison_df, mean_slope_diff_df, all=TRUE)
  write.csv(quantitative_comparison_df, paste0(plot_output_filepath, '/comparison_metric_table_prevalence_age.csv'), row.names=FALSE)
  
  
  # compare simulation and benchmark simulation results
  if('benchmark_value' %in% colnames(combined_df)){
    compare_benchmarks_output = compare_benchmark(combined_df)
    ggsave(filename=paste0(plot_output_filepath, '/scatter_benchmark_prevalence_age.png'), plot=compare_benchmarks_output, height=4.5, width=8, units='in')
    new_sim_loglik = get_prev_likelihood(combined_df, sim_column='simulation_value')
    bench_sim_loglik = get_prev_likelihood(combined_df, sim_column='benchmark_value')
    colnames(new_sim_loglik)[which(colnames(new_sim_loglik)=='loglikelihood')] = 'loglikelihood_new_sim'
    colnames(bench_sim_loglik)[which(colnames(bench_sim_loglik)=='loglikelihood')] = 'loglikelihood_benchmark_sim'
    loglikelihood_comparison = merge(new_sim_loglik, bench_sim_loglik)
    write.csv(loglikelihood_comparison, paste0(plot_output_filepath, '/loglikelihood_benchmark_comparison_prevalence_age.csv'), row.names=FALSE)
  }
}





################################################################################# 
########################### Parasite density by age  ##################################
################################################################################# 


generate_parasite_density_outputs = function(coord_csv, simulation_output_filepath, base_reference_filepath, plot_output_filepath, benchmark_simulation_filepath=NA){
  
  par_dens_sites = coord_csv$site[intersect(which(!is.na(coord_csv$site)), which(coord_csv$age_parasite_density==1))]
  # determine which of the parasite-density sites have the relevant simulation output
  available_sites = c()
  for (ii in 1:length(par_dens_sites)){
    if (file.exists(paste0(simulation_output_filepath, '/', par_dens_sites[ii], '/parasite_densities_by_age_month.csv'))){
      available_sites = c(available_sites, par_dens_sites[ii])
    }
  }
  
  # iterate through sites, grabbing relevant reference and simulation data to plot; also combine data into a dataframe containing all sites
  all_sim_sites = data.frame()
  all_ref_sites = data.frame()
  all_bench_sites = data.frame()
  for (ss in 1:length(available_sites)){
    cur_site = available_sites[ss]
    sim_df = read.csv(paste0(simulation_output_filepath, '/', cur_site, '/parasite_densities_by_age_month.csv'))
    age_agg_sim_df = get_age_bin_averages(sim_df)
    
    filepath_ref = paste0(base_reference_filepath, '/', coord_csv$age_parasite_density_ref[which(coord_csv$site == cur_site)])
    ref_df = read.csv(filepath_ref)
    ref_df = ref_df[tolower(ref_df$Site) == tolower(cur_site),]
    ref_df$Site = tolower(ref_df$Site)
    
    if(!is.na(benchmark_simulation_filepath)){
      bench_df = read.csv(paste0(benchmark_simulation_filepath, '/', cur_site, '/parasite_densities_by_age_month.csv'))
      age_agg_bench_df = get_age_bin_averages(sim_df=bench_df)
    } else age_agg_bench_df = data.frame()
    
    gg_plots = plot_par_dens_ref_sim_comparison(age_agg_sim_df, ref_df, age_agg_bench_df)
    gg_plots[[2]] = gg_plots[[2]] + ggtitle(available_sites[ss])
    gg_plots[[3]] = gg_plots[[3]] + ggtitle(available_sites[ss])
    ggsave(filename=paste0(plot_output_filepath, '/site_compare_par_dens_age_', cur_site, '.png'), plot=gg_plots[[2]], width=8, height=6, units='in')
    ggsave(filename=paste0(plot_output_filepath, '/site_compare_gamet_dens_age_', cur_site, '.png'), plot=gg_plots[[3]], width=8, height=6, units='in')
    
    all_sim_sites = merge(all_sim_sites, age_agg_sim_df, all=TRUE)
    all_ref_sites = merge(all_ref_sites, ref_df, all=TRUE)
    all_bench_sites = merge(all_bench_sites, age_agg_bench_df, all=TRUE)
  }
  
  loglik_df = get_dens_likelihood(sim_df=all_sim_sites, ref_df=all_ref_sites)
  write.csv(paste0(plot_output_filepath, '/loglikelihoods_par_dens.csv'))
}




################################################################################# 
###########################  Infectiouss to vectors  ##################################
################################################################################# 

generate_infectiousness_outputs = function(coord_csv, simulation_output_filepath, base_reference_filepath, plot_output_filepath, benchmark_simulation_filepath=NA){
  
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
    gg_plot = gg_plot + ggtitle(available_sites[ss])
    ggsave(filename=paste0(plot_output_filepath, '/site_compare_infectiousness_', cur_site, '.png'), plot=gg_plot, width=7.5, height=6, units='in')
  }
  
}




################################################################################# 
###########################  Duration of infection  ##################################
################################################################################# 

generate_age_infection_duration_outputs = function(coord_csv, simulation_output_filepath, base_reference_filepath, plot_output_filepath, pos_thresh_dens=0.5, duration_bins=c(seq(0,350,50), 500), benchmark_simulation_filepath=NA){
  
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
    
    ggsave(filename=paste0(plot_output_filepath, '/site_compare_infect_duration_', cur_site, '.png'), plot=gg_plots[[1]], height=4, width=8, units='in')
    ggsave(filename=paste0(plot_output_filepath, '/site_compare_age_infect_duration_', cur_site, '.png'), plot=gg_plots[[2]], height=5, width=8, units='in')
    ggsave(filename=paste0(plot_output_filepath, '/site_compare_duration_measures_', cur_site, '.png'), plot=gg_plots[[3]], height=4, width=8, units='in')
  }
}





################################################################################# 
###########################    ##################################
################################################################################# 





