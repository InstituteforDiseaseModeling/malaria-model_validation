# plot_par_dens_comparison.R


library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)
library(RColorBrewer)


############################ setup: filepaths and info on simulations run #################################
filepath_ref = "/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/reference_data/par_dens_extract_Laye_2007.csv"
base_filepath_sim_output = '/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/simulation_output/par_dens_age'


# # sweep of simulations across seasonalities, EIRs, and CMs
# filepath_sim = file.path(base_filepath_sim_output, 'sweepSeasonEIRCM', 'summary_data_final.csv')

# coordinator csv - get site names for parasite density simulations
simulation_coordinator_path = "/Users/moniqueam/Documents/malaria-model_validation/simulation_inputs/simulation_coordinator.csv"



############################ helper functions #########################################

# get average parasite densities in each age bin, weighting all ages in bin equally (e.g., not weighted by population size)
get_age_bin_averages = function(sim_df){
  age_bins = unique(sim_df$agebin)
  # remove rows where there are zero people of the measured age bin in the simulation
  sim_df = sim_df[sim_df$Pop > 0,]
  # get average across all years in age bins and across simulation run seeds
  age_agg_sim_df = sim_df %>% group_by(month, agebin, densitybin, Site) %>%
    summarise(asexual_par_dens = mean(asexual_par_dens), 
              gametocyte_dens = mean(gametocyte_dens), 
              Pop = mean(Pop))
  return(age_agg_sim_df)
}



# stacked barplots of parasite density bins by age
plot_ref_sim_comparison = function(age_agg_sim_df, ref_df){
  
  # subset simulation output to months in reference dataset
  months = sort(unique(ref_df$month))
  cur_df = age_agg_sim_df[age_agg_sim_df$month %in% months,]

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
  ggplot(combined_df, aes(fill=densitybin, y=asexual_par_dens, x=agebin)) + 
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

  # plot
  ggplot(combined_df0, aes(x=densitybin, y=asexual_par_dens, color=source)) + 
    geom_line(size=2) + 
    geom_point() +
    scale_x_continuous(trans='log10') +
    geom_errorbar(aes(ymin=min_asex, ymax=max_asex), width=0.2) +
    # scale_fill_brewer(palette = "BrBG") +
    # scale_fill_manual(values=colors, limits=names(colors)) +
    facet_grid(agebin~month)
  

  
  # = = = = = = = = = #
  # plot cumulative densities
  # = = = = = = = = = #
  
  
}



############################ read in data and create plots #########################################

coord_csv = read.csv(simulation_coordinator_path)
par_dens_sites = coord_csv$site[intersect(which(!is.na(coord_csv$site)), which(coord_csv$age_parasite_density==1))]

# for (ss in 1:length(par_dens_sites)){}
cur_site = "Laye_2007"
sim_df = read.csv(paste0(base_filepath_sim_output, '/', cur_site, '/parasite_densities_by_age_month.csv'))
age_agg_sim_df = get_age_bin_averages(sim_df)
ref_df = read.csv(filepath_ref)
