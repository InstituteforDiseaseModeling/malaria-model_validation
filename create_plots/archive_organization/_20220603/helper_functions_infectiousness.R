# helper_functions_infectiousness.R
library(ggplot2)




########################### create infectiousness plots  ##################################

plot_infectiousness_ref_sim_comparison = function(sim_df, ref_df){
  months_of_year = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
  
  # remove simulation rows with zero pop
  sim_df = sim_df[sim_df$Pop>0,]
  
  # subset simulation to months in reference df
  ref_months = unique(ref_df$month)
  sim_df = sim_df[sim_df$month %in% ref_months,]
  
  
  # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
  # translate the frequencies from simulation output into fraction of individuals in each 
  #     {age bin, month, densitybin, run number} group are in each infectiousness bin
  sim_df$infectiousness_bin_count = sim_df$infectiousness_bin_freq * sim_df$Pop  # note, since this is an average over the reporting period, these may not be whole numbers
  
  # check that people are never infectious while they have no parasites - send a warning if not
  if(nrow(sim_df[sim_df$infectiousness_bin == 0 & sim_df$infectiousness_bin > 0 & sim_df$infectiousness_bin_count > 0,])>0){
    warning('Some individuals without parasites are reporting that they are infectious to mosquitoes... this suggests a bug.')
  }
  
  # aggregate individuals within an age bin (across years)
  sim_df_agg1 = sim_df %>% group_by(Site, month, Run_Number, agebin, densitybin, infectiousness_bin) %>%
    summarise(infectiousness_bin_count = sum(infectiousness_bin_count), 
              Pop = sum(Pop))
  # get total within each group (for denominator)
  sim_df_group_total = sim_df_agg1 %>% group_by(Site, month, Run_Number, agebin, densitybin) %>%
    summarise(infectiousness_group_sum = sum(infectiousness_bin_count))
  # calculate proportion of all individuals in a group fell in each infectiousness bin
  sim_df_agg1 = merge(sim_df_agg1, sim_df_group_total, all=TRUE)
  sim_df_agg1$group_infectiousness_freq = sim_df_agg1$infectiousness_bin_count / sim_df_agg1$infectiousness_group_sum
  # # check that sum within a group is 1
    # sim_subset = sim_df_agg1[sim_df_agg1$month==1 & sim_df_agg1$agebin==5 & sim_df_agg1$Run_Number == 0 & sim_df_agg1$densitybin==500,]
    # sum(sim_subset$group_infectiousness_freq)
  
  # get simulation average across seeds and within age groups (assuming equal population sizes for all ages in bin)
  sim_df_agg2 = sim_df_agg1 %>% group_by(Site, month, agebin, densitybin, infectiousness_bin) %>%
    summarise(infect_sd = sd(group_infectiousness_freq),
              infectiousness_bin_freq = mean(group_infectiousness_freq))
  # check that the sum within all groups is 1
  sim_df_check = sim_df_agg2 %>% group_by(Site, month, agebin, densitybin) %>%
    summarise(infectiousness_dens_bin_sum = sum(infectiousness_bin_freq))
  if(!all(sim_df_check$infectiousness_dens_bin_sum-0.001 < 1) | !all(sim_df_check$infectiousness_dens_bin_sum+0.001 > 1)) warning("The sum of infectiousness bin frequencies is not 1 for at least one group. Recommend checking for bugs.")
  
  
  # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
  # standardize column names and merge simulation and reference data frames
  colnames(sim_df_agg2)[colnames(sim_df_agg2) == 'Site'] = 'site'
  colnames(sim_df_agg2)[colnames(sim_df_agg2) == 'infectiousness_bin'] = 'fraction_infected_bin'
  colnames(ref_df)[colnames(ref_df) == 'freq_frac_infect'] = 'infectiousness_bin_freq'

  sim_df_agg2$source='simulation'
  ref_df$source='reference'  
  
  combined_df = merge(sim_df_agg2, ref_df, all=TRUE)
  combined_df$source = factor(combined_df$source, levels=c('reference', 'simulation'))
  
  # change facet values to intuitive labels
  combined_df$month = months_of_year[combined_df$month]
  combined_df$month = factor(combined_df$month, levels=months_of_year)
  all_age_bins = sort(unique(combined_df$agebin))
  age_bin_labels = paste0('<=', all_age_bins[1], ' years')
  for(aa in 1:(length(all_age_bins)-1)){
    age_bin_labels = c(age_bin_labels, paste0(all_age_bins[aa], '-', all_age_bins[aa+1], ' years'))
  }
  combined_df$agebin_index = match(combined_df$agebin, all_age_bins)
  combined_df$agebin = age_bin_labels[combined_df$agebin_index]
  combined_df$agebin = factor(combined_df$agebin, levels = age_bin_labels)
  
  # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
  # create plot
  gg = ggplot(data=combined_df, aes(y=fraction_infected_bin, x=densitybin, size=infectiousness_bin_freq, color=source, fill=source))+
    geom_point(alpha=0.5) +
    scale_x_log10() +
    ylab('percent of mosquitoes infected upon feeding') + #('percent of individuals (in each age-density group) who fall in each infectiousness bin') +
    xlab('gametocyte density') +
    labs(size = 'fraction of individuals') +
    scale_color_manual(values = c("reference" = rgb(169/255,23/255,23/255, alpha=0.8),
                                  "simulation"=rgb(0/255,124/255,180/255, alpha=0.8))) +
    scale_fill_manual(values = c("reference" = rgb(169/255,23/255,23/255, alpha=0.8),
                                 "simulation"=rgb(0/255,124/255,180/255, alpha=0.8))) +
    facet_grid(agebin~month)
  return(gg)
}




########################### main coordinator function  ##################################

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
