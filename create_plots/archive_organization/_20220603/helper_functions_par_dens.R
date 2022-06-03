# plot_par_dens_comparison.R


library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)
library(RColorBrewer)


############################ helper function #########################################

# get average parasite densities in each age bin, weighting all ages in bin equally (e.g., not weighted by population size)
get_age_bin_averages = function(sim_df){
  age_bins = unique(sim_df$agebin)
  # remove rows where there are zero people of the measured age bin in the simulation
  sim_df = sim_df[sim_df$Pop > 0,]
  # get the simulation mean across runs
  if (all(sim_df$Pop ==sim_df$Pop[1])){
    # get average across all years in age bins and across simulation run seeds
    age_agg_sim_df = sim_df %>% group_by(month, agebin, densitybin, Site) %>%
      summarise(asexual_par_dens_freq = mean(asexual_par_dens_freq), 
                gametocyte_dens_freq = mean(gametocyte_dens_freq), 
                Pop = mean(Pop))
  } else{
    warning("Different population sizes found across years within an age group... need to set up weighted averaging of parasite densities.") 
    # If this warning is triggered, use the population size to calculate the number of individuals in each density bin, then aggregate the sum, then divide by the total aggregated population
  }
  return(age_agg_sim_df)
}



########################## plot parasite density comparisons with reference ####################

# stacked barplots of parasite density bins by age
plot_par_dens_ref_sim_comparison = function(age_agg_sim_df, ref_df, age_agg_bench_df){
  months_of_year = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
  
  # subset simulation output to months in reference dataset
  months = sort(unique(ref_df$month))
  cur_df = age_agg_sim_df[age_agg_sim_df$month %in% months,]
  
  
  # if the maximum reference density bin is < (maximum simulation density bin / 1000), aggregate all simulation densities >= max ref bin into the max ref bin
  #   the final bin will be all densities equal to or above that value
  max_ref_dens = max(ref_df$densitybin, na.rm=TRUE)
  if(max_ref_dens < (max(cur_df$densitybin, na.rm=TRUE)/1000)){
    # get sum of frequencies within higher bins
    all_higher_dens = cur_df[cur_df$densitybin >= max_ref_dens,]
    sim_agg_higher_dens = all_higher_dens %>% group_by(month, agebin, Site) %>%
      summarise(densitybin = min(densitybin),
                asexual_par_dens_freq = sum(asexual_par_dens_freq),
                gametocyte_dens_freq = sum(gametocyte_dens_freq),
                Pop = mean(Pop))
    # remove higher density bins from df
    cur_df_lower = cur_df[cur_df$densitybin < max_ref_dens,]
    # add back in the aggregated frequencies across higher density bins
    cur_df = merge(cur_df_lower, sim_agg_higher_dens, all=TRUE)
  }
  
  # add zeros for unobserved reference densities up to max_ref_dens
  all_zeros_df = cur_df[,c('month', 'agebin', 'densitybin', 'Site')]
  ref_df = merge(ref_df, all_zeros_df, all=TRUE)
  ref_df[is.na(ref_df)] = 0

  
  # combine reference and simulation dataframes
  cur_df$source = 'simulation'
  ref_df$source = 'reference'
  combined_df0 = full_join(cur_df, ref_df)
  

  # = = = = = = = = = #
  # stacked barplots
  # = = = = = = = = = #
  # change type to factors for barplot groupings
  combined_df = combined_df0
  combined_df$densitybin = factor(combined_df$densitybin, levels=sort(unique(combined_df$densitybin)))
  combined_df$agebin = factor(combined_df$agebin, levels=sort(unique(combined_df$agebin)))
  
  # colors
  num_colors = ifelse(length(unique(combined_df$densitybin)) %% 2 ==0, length(unique(combined_df$densitybin))+1, length(unique(combined_df$densitybin)))
  colors = brewer.pal(n=num_colors, name='BrBG')
  names(colors) = sort(unique(combined_df$densitybin))
  # plot
  gg1=ggplot(combined_df, aes(fill=densitybin, y=asexual_par_dens_freq, x=agebin)) + 
    geom_bar(position="stack", stat="identity") + 
    # scale_fill_brewer(palette = "BrBG") +
    scale_fill_manual(values=colors, limits=names(colors)) +
    facet_grid(month~source)
  
  
  
  # = = = = = = = = = #
  # grid of line plots
  # = = = = = = = = = #
  
  # calculate reference error bounds using Jerrerys interval
  ci_width = 0.95
  alpha = 1-ci_width
  combined_df0$min_asex = NA
  combined_df0$max_asex = NA
  combined_df0$min_gamet = NA
  combined_df0$max_gamet = NA
  for(rr in 1:nrow(combined_df0)){
    if(combined_df0$source[rr] == 'reference'){
      if((combined_df0$count_asex[rr]>0) & (combined_df0$count_asex[rr]<combined_df0$bin_total_asex[rr])){
        combined_df0$min_asex[rr] = qbeta(p=(alpha/2), shape1=(combined_df0$count_asex[rr]+0.5), shape2=(combined_df0$bin_total_asex[rr] - combined_df0$count_asex[rr] + 0.5))
        combined_df0$max_asex[rr] = qbeta(p=(1-alpha/2), shape1=(combined_df0$count_asex[rr]+0.5), shape2=(combined_df0$bin_total_asex[rr] - combined_df0$count_asex[rr] + 0.5))
      }
      if((combined_df0$count_gamet[rr]>0) & (combined_df0$count_gamet[rr]<combined_df0$bin_total_gamet[rr])){
        combined_df0$min_gamet[rr] = qbeta(p=(alpha/2), shape1=(combined_df0$count_gamet[rr]+0.5), shape2=(combined_df0$bin_total_gamet[rr] - combined_df0$count_gamet[rr] + 0.5))
        combined_df0$max_gamet[rr] = qbeta(p=(1-alpha/2), shape1=(combined_df0$count_gamet[rr]+0.5), shape2=(combined_df0$bin_total_gamet[rr] - combined_df0$count_gamet[rr] + 0.5))
      }
    }
  }
  
  # change facet values to intuitive labels
  combined_df0$month = months_of_year[combined_df0$month]
  combined_df0$month = factor(combined_df0$month, levels=months_of_year)
  all_age_bins = sort(unique(combined_df0$agebin))
  age_bin_labels = paste0('<=', all_age_bins[1], ' years')
  for(aa in 1:(length(all_age_bins)-1)){
    age_bin_labels = c(age_bin_labels, paste0(all_age_bins[aa], '-', all_age_bins[aa+1], ' years'))
  }
  combined_df0$agebin_index = match(combined_df0$agebin, all_age_bins)
  combined_df0$agebin = age_bin_labels[combined_df0$agebin_index]
  combined_df0$agebin = factor(combined_df0$agebin, levels = age_bin_labels)

  # plot asexual densities
  gg2=ggplot(combined_df0, aes(x=densitybin, y=asexual_par_dens_freq, color=source), alpha=0.8) + 
    geom_line(size=2) + 
    geom_point() +
    scale_x_continuous(trans='log10') +
    geom_errorbar(aes(ymin=min_asex, ymax=max_asex), width=0.2) +
    theme_bw() +
    ylab('fraction of population') +
    xlab('asexual parasite density bin') +
    scale_color_manual(values = c("reference" = rgb(169/255,23/255,23/255, alpha=0.8),
                                  "simulation"=rgb(0/255,124/255,180/255, alpha=0.8))) +
    # scale_fill_brewer(palette = "BrBG") +
    # scale_fill_manual(values=colors, limits=names(colors)) +
    facet_grid(agebin~month)
  
  
  # plot gametocyte densities
  gg3=ggplot(combined_df0, aes(x=densitybin, y=gametocyte_dens_freq, color=source)) + 
    geom_line(size=2) + 
    geom_point() +
    scale_x_continuous(trans='log10') +
    geom_errorbar(aes(ymin=min_gamet, ymax=max_gamet), width=0.2) +
    theme_bw() +
    ylab('fraction of population') +
    xlab('gametocyte density bin') +
    scale_color_manual(values = c("reference" = rgb(169/255,23/255,23/255, alpha=0.8),
                                  "simulation"=rgb(0/255,124/255,180/255, alpha=0.8))) +
    facet_grid(agebin~month)

  return(list(gg1,gg2, gg3))
}




########################### main coordinator function  ##################################

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


