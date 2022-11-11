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
  #' From simulation output and matched reference data, create plots and quantitative comparisons for all sites associated with the incidence-by-age validation relationship. 
  #' 
  #' @param coord_csv A dataframe detailing the sites simulated for each validation relationship and the corresponding reference dataset
  #' @param simulation_output_filepath The filepath where simulation output is located
  #' @param base_reference_filepath The filepath where reference datasets are located
  #' @param plot_output_filepath The filepath to the directory where plots should be created
  #' @param benchmark_simulation_filepath The filepath where benchmark simulation output is located. If NA, no comparisons are made against benchmark simulations

  # get formatted dataframe with reference and simulation incidence data from all relevant sites
  combined_df = prepare_inc_df(coord_csv, simulation_output_filepath, base_reference_filepath, benchmark_simulation_filepath)

  # create plots comparing reference and simulation outputs
  gg_plot = plot_inc_ref_sim_comparison(combined_df)
  ggsave(filename=paste0(plot_output_filepath, '/site_compare_incidence_age.png'), plot=gg_plot, height=2*ceiling(length(unique(combined_df$Site))/4), width=7.5, units='in')
  
  
  # additional quantitative comparisons and metrics between simulation and reference data
  
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
  if('benchmark' %in% colnames(combined_df)){
    compare_benchmarks_output = compare_benchmark(combined_df)
    ggsave(filename=paste0(plot_output_filepath, '/scatter_benchmark_incidence_age.png'), plot=compare_benchmarks_output, height=4.5, width=8, units='in')
    mean_diff_df_bench = calc_mean_rel_diff(combined_df, sim_colname='benchmark')
    
    add_to_summary_table(combined_df=combined_df, plot_output_filepath=plot_output_filepath, validation_relationship_name='age_incidence')
  }
}



################################################################################# 
########################### Prevalence by age  ##################################
################################################################################# 

generate_age_prevalence_outputs = function(coord_csv, simulation_output_filepath, base_reference_filepath, plot_output_filepath, benchmark_simulation_filepath=NA){
  #' From simulation output and matched reference data, create plots and quantitative comparisons for all sites associated with the prevalence-by-age validation relationship. 
  #' 
  #' @param coord_csv A dataframe detailing the sites simulated for each validation relationship and the corresponding reference dataset
  #' @param simulation_output_filepath The filepath where simulation output is located
  #' @param base_reference_filepath The filepath where reference datasets are located
  #' @param plot_output_filepath The filepath to the directory where plots should be created
  #' @param benchmark_simulation_filepath The filepath where benchmark simulation output is located. If NA, no comparisons are made against benchmark simulations

  # get formatted dataframe with reference and simulation prevalence data from all relevant sites
  combined_df = prepare_prev_df(coord_csv, simulation_output_filepath, base_reference_filepath, benchmark_simulation_filepath)

  # create plots comparing reference and simulation outputs
  gg_plot = plot_prev_ref_sim_comparison(combined_df)
  ggsave(filename=paste0(plot_output_filepath, '/site_compare_prevalence_age.png'), plot=gg_plot, height=9, width=8, units='in')
  
  
  # additional quantitative comparisons and metrics
  
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
  if('benchmark' %in% colnames(combined_df)){
    compare_benchmarks_output = compare_benchmark(combined_df)
    ggsave(filename=paste0(plot_output_filepath, '/scatter_benchmark_prevalence_age.png'), plot=compare_benchmarks_output, height=4.5, width=8, units='in')
    new_sim_loglik = get_prev_loglikelihood(combined_df, sim_column='simulation')
    bench_sim_loglik = get_prev_loglikelihood(combined_df, sim_column='benchmark')
    colnames(new_sim_loglik)[which(colnames(new_sim_loglik)=='loglikelihood')] = 'loglikelihood_new_sim'
    colnames(bench_sim_loglik)[which(colnames(bench_sim_loglik)=='loglikelihood')] = 'loglikelihood_benchmark_sim'
    loglikelihood_comparison = merge(new_sim_loglik, bench_sim_loglik)
    write.csv(loglikelihood_comparison, paste0(plot_output_filepath, '/loglikelihood_prevalence_age.csv'), row.names=FALSE)
    
    add_to_summary_table(combined_df=combined_df, plot_output_filepath=plot_output_filepath, validation_relationship_name='age_prevalence')
  }
}





################################################################################# 
########################### Parasite density by age  ##################################
################################################################################# 


generate_parasite_density_outputs = function(coord_csv, simulation_output_filepath, base_reference_filepath, plot_output_filepath, benchmark_simulation_filepath=NA){
  #' From simulation output and matched reference data, create plots and quantitative comparisons for all sites associated with the parasite density-by-age validation relationship. 
  #' 
  #' @param coord_csv A dataframe detailing the sites simulated for each validation relationship and the corresponding reference dataset
  #' @param simulation_output_filepath The filepath where simulation output is located
  #' @param base_reference_filepath The filepath where reference datasets are located
  #' @param plot_output_filepath The filepath to the directory where plots should be created
  #' @param benchmark_simulation_filepath The filepath where benchmark simulation output is located. If NA, no comparisons are made against benchmark simulations

  # get formatted dataframe with reference and simulation prevalence data from all relevant sites
  combined_dfs = prepare_dens_df(coord_csv, simulation_output_filepath, base_reference_filepath, benchmark_simulation_filepath)
  combined_df_asex = combined_dfs[[1]]
  combined_df_gamet = combined_dfs[[2]]

  # asexual parasite density
  plot_output = plot_par_dens_ref_sim_comparison(combined_df = combined_df_asex)
  gg_barplot = plot_output[[1]]
  line_plot_list = plot_output[[2]]
  all_sites = plot_output[[3]]
  ggsave(filename=paste0(plot_output_filepath, '/site_compare_barplot_asex_dens_age_.png'), plot=gg_barplot, width=10, height=20, units='in')
  for(ss in 1:length(all_sites)){
    ggsave(filename=paste0(plot_output_filepath, '/site_compare_asex_dens_age_', all_sites[ss], '.png'), plot=line_plot_list[[ss]], width=8, height=6, units='in')
  }  
  
  # gametocyte density
  plot_output = plot_par_dens_ref_sim_comparison(combined_df = combined_df_gamet)
  gg_barplot = plot_output[[1]]
  line_plot_list = plot_output[[2]]
  all_sites = plot_output[[3]]
  ggsave(filename=paste0(plot_output_filepath, '/site_compare_barplot_gamet_dens_age_.png'), plot=gg_barplot, width=5, height=15, units='in')
  for(ss in 1:length(all_sites)){
    ggsave(filename=paste0(plot_output_filepath, '/site_compare_gamet_dens_age_', all_sites[ss], '.png'), plot=line_plot_list[[ss]], width=8, height=6, units='in')
  }  
  
  
  
  
  if('benchmark' %in% colnames(combined_df_asex)){
    compare_benchmarks_output = compare_benchmark(combined_df_asex)
    ggsave(filename=paste0(plot_output_filepath, '/scatter_benchmark_asex_dens.png'), plot=compare_benchmarks_output, height=4.5, width=8, units='in')
    compare_benchmarks_output = compare_benchmark(combined_df_gamet)
    ggsave(filename=paste0(plot_output_filepath, '/scatter_benchmark_gamet_dens.png'), plot=compare_benchmarks_output, height=4.5, width=8, units='in')
    
    loglik_df_asex = get_dens_loglikelihood(combined_df=combined_df_asex, sim_column='simulation')
    colnames(loglik_df_asex)[colnames(loglik_df_asex)=='loglikelihood'] = 'loglike_asex'
    loglik_df_asex_bench = get_dens_loglikelihood(combined_df=combined_df_asex, sim_column='benchmark')
    colnames(loglik_df_asex_bench)[colnames(loglik_df_asex_bench)=='loglikelihood'] = 'benchmark_loglike_asex'
    loglik_df_asex = merge(loglik_df_asex, loglik_df_asex_bench, all=TRUE)
    
    loglik_df_gamet = get_dens_loglikelihood(combined_df=combined_df_gamet, sim_column='simulation')
    colnames(loglik_df_gamet)[colnames(loglik_df_gamet)=='loglikelihood'] = 'loglike_gamet'
    loglik_df_gamet_bench = get_dens_loglikelihood(combined_df=combined_df_gamet, sim_column='benchmark')
    colnames(loglik_df_gamet_bench)[colnames(loglik_df_gamet_bench)=='loglikelihood'] = 'benchmark_loglike_gamet'
    loglik_df_gamet = merge(loglik_df_gamet, loglik_df_gamet_bench, all=TRUE)

    loglik_df = merge(loglik_df_asex, loglik_df_gamet, all=TRUE)
    write.csv(loglik_df, paste0(plot_output_filepath, '/loglikelihoods_par_dens.csv'))
    
    add_to_summary_table(combined_df=combined_df_asex, plot_output_filepath=plot_output_filepath, validation_relationship_name='asexual_par_dens')
    add_to_summary_table(combined_df=combined_df_gamet, plot_output_filepath=plot_output_filepath, validation_relationship_name='gamet_par_dens')
    
  }
}




################################################################################# 
###########################  Infectiousness to vectors  ##################################
################################################################################# 

generate_infectiousness_outputs = function(coord_csv, simulation_output_filepath, base_reference_filepath, plot_output_filepath, benchmark_simulation_filepath=NA){
  #' From simulation output and matched reference data, create plots and quantitative comparisons for all sites associated with the infectiousness-to-vectors validation relationship. 
  #' 
  #' @param coord_csv A dataframe detailing the sites simulated for each validation relationship and the corresponding reference dataset
  #' @param simulation_output_filepath The filepath where simulation output is located
  #' @param base_reference_filepath The filepath where reference datasets are located
  #' @param plot_output_filepath The filepath to the directory where plots should be created
  #' @param benchmark_simulation_filepath The filepath where benchmark simulation output is located. If NA, no comparisons are made against benchmark simulations

  combined_df = prepare_infect_df(coord_csv, simulation_output_filepath, base_reference_filepath, benchmark_simulation_filepath)

  plot_output = plot_infectiousness_ref_sim_comparison(combined_df)
  plot_list = plot_output[[1]]
  all_sites = plot_output[[2]]
  for(ss in 1:length(all_sites)){
    ggsave(filename=paste0(plot_output_filepath, '/site_compare_infectiousness_', all_sites[ss], '.png'), plot=plot_list[[ss]], width=7.5, height=6, units='in')
  } 
  
  if('benchmark' %in% colnames(combined_df)){
    compare_benchmarks_output = compare_benchmark(combined_df)
    ggsave(filename=paste0(plot_output_filepath, '/scatter_benchmark_infectiousness.png'), plot=compare_benchmarks_output, height=4.5, width=8, units='in')
    
    # TODO: add likelihood and other quantitative comparisons
    add_to_summary_table(combined_df=combined_df, plot_output_filepath=plot_output_filepath, validation_relationship_name='infectiousness')
    
  }
}




################################################################################# 
###########################  Duration of infection  ##################################
################################################################################# 

generate_age_infection_duration_outputs = function(coord_csv, simulation_output_filepath, base_reference_filepath, plot_output_filepath, pos_thresh_dens=0.5, duration_bins=c(seq(0,350,50), 500), benchmark_simulation_filepath=NA){
  #' From simulation output and matched reference data, create plots and quantitative comparisons for all sites associated with the duration-of-infection validation relationship. 
  #' 
  #' @param coord_csv A dataframe detailing the sites simulated for each validation relationship and the corresponding reference dataset
  #' @param simulation_output_filepath The filepath where simulation output is located
  #' @param base_reference_filepath The filepath where reference datasets are located
  #' @param plot_output_filepath The filepath to the directory where plots should be created
  #' @param pos_thresh_dens A number giving the minimum true asexual parasite density a simulated individual must have to be considered positive
  #' @param duration_bins A monotonically-increasing vector of numbers giving the plotted bin breaks for the duration (in days) individuals remain infected
  #' @param benchmark_simulation_filepath The filepath where benchmark simulation output is located. If NA, no comparisons are made against benchmark simulations
  
  
  # TODO: add benchmark simulation support, add quantitative comparisons
  
  
  # determine which of the infectiousness sites have the relevant simulation output
  available_sites = get_available_sites_for_relationship(coord_csv, simulation_output_filepath, relationship_name='infection_duration', relationship_sim_filename='patient_reports.csv')

  for (ss in 1:length(available_sites)){
    cur_site = available_sites[ss]

    filepath_ref = paste0(base_reference_filepath, '/', coord_csv$infection_duration_ref[which(coord_csv$site == cur_site)])
    ref_df = read.csv(filepath_ref)
    ref_df = ref_df[tolower(ref_df$site) == tolower(cur_site),]
    ref_df$date = as.Date(ref_df$date)
    
    sim_dir = paste0(simulation_output_filepath, '/', cur_site)
    sim_data = get_sim_survey(sim_dir=sim_dir, ref_df=ref_df)
    
    # create and save comparison plots
    gg1 = plot_infection_duration_dist(ref_df=ref_df, sim_data=sim_data, pos_thresh_dens=pos_thresh_dens, duration_bins=duration_bins)
    gg2 = plot_infection_duration_dist_by_age(ref_df=ref_df, sim_data=sim_data, pos_thresh_dens=pos_thresh_dens, duration_bins=duration_bins)
    gg3 = create_barplot_frac_comparison(ref_df=ref_df, sim_data=sim_data, pos_thresh_dens=pos_thresh_dens)
    
    ggsave(filename=paste0(plot_output_filepath, '/site_compare_infect_duration_', cur_site, '.png'), plot=gg1, height=4, width=8, units='in')
    ggsave(filename=paste0(plot_output_filepath, '/site_compare_infect_duration_age_', cur_site, '.png'), plot=gg2, height=5, width=8, units='in')
    ggsave(filename=paste0(plot_output_filepath, '/site_compare_infect_duration_measures_', cur_site, '.png'), plot=gg3, height=4, width=8, units='in')
  }
}







