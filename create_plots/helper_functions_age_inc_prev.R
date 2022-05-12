# helper_functions_age_inc.R

# plot the match between the reference datasets and matched simulations

library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)
library(ggpubr)



############################ helper functions #########################################
get_substr = function(site_name_str, index){
  strsplit(site_name_str, "_")[[1]][index]
}
get_mean_from_upper_age = function(cur_age, upper_ages){
  mean_ages = (c(0, upper_ages[1:(length(upper_ages)-1)]) + upper_ages) / 2
  return(mean_ages[which(upper_ages == cur_age)])
}



########################## plot age-incidence comparisons with reference ####################
plot_inc_ref_sim_comparison = function(sim_df, ref_df){
  # scale down incidence in simulation according to probability of detecting a case
  sim_df$Incidence = sim_df$Incidence * sim_df$p_detect_case

  # set up reference and simulation dataset columns
  ref_df$Incidence = ref_df$INC / 1000
  ref_df$mean_age = (ref_df$INC_LAR + ref_df$INC_UAR)/2
  ref_df = data.frame('Incidence'=ref_df$Incidence, 'mean_age'=ref_df$mean_age, 'Site'=ref_df$Site, 'Pop_size'=ref_df$POP, 'year'=ref_df$START_YEAR)
  sim_df = data.frame('Incidence'=sim_df$Incidence, 'mean_age'=sim_df$mean_age, 'Site'=sim_df$Site, 'Pop_size'=NA, 'year'=NA)
  ref_df$source = 'reference'
  sim_df$source = 'simulation'

  df_combined = rbind(sim_df, ref_df)
  df_combined$source = factor(df_combined$source, levels=c('reference', 'simulation'))
  
  gg = ggplot(df_combined, aes(x=mean_age, y=Incidence, color=source, group=year)) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = c("reference" = rgb(169/255,23/255,23/255, alpha=0.8),
                                  "simulation"=rgb(0/255,124/255,180/255, alpha=0.8))) +
    xlab('age (midpoint of age bin)') +
    ylab('incidence') +
    facet_wrap('Site', ncol=4) +
    theme_bw()  
  
  return(gg)
}






########################## plot age-prevalence comparisons with reference without sweep background ####################
plot_prev_ref_sim_comparison = function(sim_df, ref_df){
  # get simulation average across seeds
  sim_df = sim_df %>% group_by(Site, mean_age, month) %>%
    summarise(prev_sd = sd(prevalence),
              prevalence = mean(prevalence))
  
  ref_df$Site = tolower(ref_df$Site)
  sim_df$Site = tolower(sim_df$Site)
  ref_df$source = 'reference'
  sim_df$source = 'simulation'
  df_combined = merge(sim_df, ref_df, all=TRUE)
  df_combined$source = factor(df_combined$source, levels=c('reference', 'simulation'))
  df_combined$site_month = paste0(df_combined$Site, '_month', df_combined$month)
  
  
  gg = ggplot(df_combined, aes(x=mean_age, y=prevalence, color=source, group=year)) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = c("reference" = rgb(169/255,23/255,23/255, alpha=0.8),
                                  "simulation"=rgb(0/255,124/255,180/255, alpha=0.8))) +
    xlab('age (midpoint of age bin)') +
    ylab('prevalence') +
    facet_wrap('site_month', ncol=4) +
    theme_bw()  
  
  return(gg)
}



########################### main coordinator function for incidence  ##################################


generate_age_incidence_outputs = function(coord_csv, simulation_output_filepath, base_reference_filepath, plot_output_filepath){
  
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
  ggsave(filename=paste0(plot_output_filepath, '/site_compare_incidence_age.png'), plot=gg_plot, height=2*ceiling(length(available_sites)/4), width=7.5, units='in')
  
  # additional quantitative comparisons and metrics
  combined_df = prepare_inc_df(sim_df, ref_df)
  
  correlation_output = corr_ref_sim_points(combined_df)
  correlation_df = correlation_output[[2]]
  slope_correlation_output = corr_ref_deriv_sim_points(combined_df)
  slope_correlation_df = slope_correlation_output[[2]]
  combined_df_with_slopes = slope_correlation_output[[3]]
  correlation_plots = ggarrange(correlation_output[[1]], slope_correlation_output[[1]], nrow=1, ncol=2, 
                                common.legend=TRUE)  #, legend.grob=get_legend(correlation_output[[1]], position = 'bottom'))
  ggsave(filename=paste0(plot_output_filepath, '/scatter_regression_incidence_age.png'), plot=correlation_plots, height=4.5, width=8, units='in')
  
  # ggsave(filename=paste0(plot_output_filepath, '/scatter_regression_incidence_age.png'), plot=correlation_output[[1]])
  # ggsave(filename=paste0(plot_output_filepath, '/scatter_regression_slope_incidence_age.png'), plot=slope_correlation_output[[1]])
  
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
}



########################### main coordinator function for prevalence  ##################################

generate_age_prevalence_outputs = function(coord_csv, simulation_output_filepath, base_reference_filepath, plot_output_filepath){
  
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
    # remove rows with population = 0 (this is relevant for cohort simulations where there is only one age group in each year and all other age groups are zero)
    sim_df_cur = sim_df_cur[sim_df_cur$Pop >0,]
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
  ggsave(filename=paste0(plot_output_filepath, '/site_compare_prevalence_age.png'), plot=gg_plot, height=9, width=8, units='in')
  
  
  # additional quantitative comparisons and metrics
  combined_df = prepare_prev_df(sim_df, ref_df)
  
  correlation_output = corr_ref_sim_points(combined_df)
  correlation_df = correlation_output[[2]]
  slope_correlation_output = corr_ref_deriv_sim_points(combined_df)
  slope_correlation_df = slope_correlation_output[[2]]
  combined_df_with_slopes = slope_correlation_output[[3]]
  correlation_plots = ggarrange(correlation_output[[1]], slope_correlation_output[[1]], nrow=1, ncol=2, 
                                common.legend=TRUE)  #, legend.grob=get_legend(correlation_output[[1]], position = 'bottom'))
  ggsave(filename=paste0(plot_output_filepath, '/scatter_regression_prevalence_age.png'), plot=correlation_plots, height=4.5, width=8, units='in')
  
  # ggsave(filename=paste0(plot_output_filepath, '/scatter_regression_incidence_age.png'), plot=correlation_output[[1]])
  # ggsave(filename=paste0(plot_output_filepath, '/scatter_regression_slope_incidence_age.png'), plot=slope_correlation_output[[1]])
  
  
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
  
}







# # set up simulation output columns
# sim_site_df = data.frame('REF_GROUP'=NA, 'Incidence'=NA, 'mean_age'=NA, 'ACD_LOCATION'=NA)
# for (ss in sites){
#   filepath_sim = file.path(paste0(base_filepath_sim_output, '/site_', ss, '/summary_data_final.csv'))
#   df = read.csv(filepath_sim)
#   upper_ages = sort(unique(df$Age))
#   df$mean_age = sapply(df$Age, get_mean_from_upper_age, upper_ages = upper_ages)
#   
#   df$ACD_LOCATION = ref_df$ACD_LOCATION[ref_df$REF_GROUP == as.numeric(ss)][1]
#   
# }
# sim_site_df = sim_site_df[-1,]
# df_combined$site_name = paste0(df_combined$ACD_LOCATION, ' (id #', df_combined$REF_GROUP, ')')
# df_combined_clean = df_combined[!(df_combined$site_name %in% paste0('Dielmo (id #', duplicate_sites, ')')),]
# 
# # format for comparison with Cameron
# gg = ggplot(df_combined, aes(x=mean_age, y=Incidence, color=source)) +
#   geom_line() +
#   geom_point() +
#   scale_color_manual(values = c("reference" = rgb(0.3, 0.3, 0.3),
#                                 "simulation"=rgb(0.7, 0.9, 0.1))) +
#   scale_x_continuous(breaks=c(0, 5, 15, 45, 90), trans='sqrt') +
#   coord_cartesian(xlim=c(-0.1, 90)) +  # 50
#   facet_wrap('site_name', ncol=3, scales='free') +
#   theme_bw()+
#   theme(panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       strip.background = element_blank(),
#       panel.border = element_rect(colour = "black", fill = NA))
# ggsave(file.path(base_filepath_sim_output, '_plots', paste0('compareSimRef_age_inc_CamFormat.png')), gg, width=6.5, height=4.75, units='in')
# 
# 
# # zoom into ages <20
# gg = ggplot(df_combined_clean, aes(x=mean_age, y=Incidence, color=source)) +
#   geom_line() +
#   geom_point() +
#   scale_color_manual(values = c("reference" = "red",
#                                 "simulation"="blue")) +
#   coord_cartesian(xlim=c(0, 20)) +
#   xlab('age (midpoint of age bin)') +
#   ylab('incidence') +
#   facet_wrap('site_name', ncol=4) +
#   theme_bw()
# ggsave(file.path(base_filepath_sim_output, '_plots', paste0('compareSimRef_age_inc_zoomed.png')), gg, width=9, height=5, units='in')



