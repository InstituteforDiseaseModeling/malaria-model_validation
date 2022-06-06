# helpers_reformat_sim_ref_dfs.R


# This script contains functions used to process, align, and format reference and simulation data. 
#    There is one main function per validation relationship, as well as a number of shared functions 
#    that are called for common data manipulations. The main function for a validation relationship 
#    takes the reference data, new simulation outputs, and benchmark simulation outputs (optional) 
#    as inputs and returns data frames that are used in downstream plotting and comparisons. 



library(stringr)
library(dplyr)
library(data.table)
library(lubridate)
library(broom)
library(tidyverse)


############################ helper functions #########################################
get_substr = function(site_name_str, index){
  strsplit(site_name_str, "_")[[1]][index]
}
get_mean_from_upper_age = function(cur_age, upper_ages){
  mean_ages = (c(0, upper_ages[1:(length(upper_ages)-1)]) + upper_ages) / 2
  return(mean_ages[which(upper_ages == cur_age)])
}



# get average parasite densities in each age bin, weighting all ages in bin equally (e.g., not weighted by population size)
get_age_bin_averages = function(sim_df){
  # remove rows where there are zero people of the measured age bin in the simulation
  sim_df = sim_df[sim_df$Pop > 0,]
  # get the simulation mean across runs
  if (all(sim_df$Pop ==sim_df$Pop[1])){
    # get average across all years in age bins and across simulation run seeds
    age_agg_sim_df = sim_df %>% group_by(month, mean_age, agebin, densitybin, Site) %>%
      summarise(asexual_par_dens_freq = mean(asexual_par_dens_freq), 
                gametocyte_dens_freq = mean(gametocyte_dens_freq), 
                Pop = mean(Pop))
  } else{
    warning("Different population sizes found across years within an age group... need to set up weighted averaging of parasite densities.") 
    # If this warning is triggered, use the population size to calculate the number of individuals in each density bin, then aggregate the sum, then divide by the total aggregated population
  }
  return(age_agg_sim_df)
}




match_sim_ref_ages = function(ref_df, sim_df, bench_df=data.frame()){
  
  sites = intersect(unique(sim_df$Site), unique(ref_df$Site))
  # check that ages match between reference and simulation. if there is a small difference (<1 year, update simulation)
  for(ss in sites){
    ages_ref = sort(unique(ref_df$mean_age[ref_df$Site == ss]))
    ages_sim = sort(unique(sim_df$mean_age[sim_df$Site == ss]))
    missing_ages_ref = ages_ref[which(!(ages_ref %in% ages_sim))]
    missing_ages_sim = ages_sim[which(!(ages_sim %in% ages_ref))]
    
    if(nrow(bench_df)>0 && !all(ages_sim == sort(unique(bench_df$mean_age[bench_df$Site==ss])))){
      warning(paste0("The age bins used in the benchmarking simulation are different from those used in the new simulation for site: ", ss, ". This benchmark simulation will be excluded."))
      bench_df = bench_df[bench_df$Site != ss,]
    }
    
    if(!all(ages_ref <= ages_sim+0.1) | !all(ages_ref >= ages_sim-0.1)){
      print(paste0('Imperfect age match between reference and simulations for site: ', ss))
      print('...  Mismatched reference / simulation ages are:')
      print(paste0('     ', missing_ages_ref, ' / ', missing_ages_sim, ','))
      print('... For age thresholds that differ by less than a year, replacing simulation age with reference age.')
      
      # check whether the missing ages are simply off by <1 year. If so, replace simulation age with nearby reference age
      for(mm in missing_ages_ref){
        sim_replace_age = missing_ages_sim[which(abs(missing_ages_sim - mm)<1)]
        if(length(sim_replace_age)==1){
          sim_df$mean_age[sim_df$Site == ss & sim_df$mean_age == sim_replace_age] = mm
          if(nrow(bench_df>0) && any(bench_df$Site ==ss)){
            bench_df$mean_age[bench_df$Site == ss & bench_df$mean_age == sim_replace_age] = mm
          }
        }
      }
      
      # update sim ages
      ages_sim = sort(unique(sim_df$mean_age[sim_df$Site == ss]))
      
      if(!all(ages_ref <= ages_sim+0.1) | !all(ages_ref >= ages_sim-0.1)){
        print('...After adjustment, there remains an imperfect match between reference and simulation age bins.')
        print(paste0('      Reference has ', length(ages_ref), ' age groups and simulation has ', length(ages_sim), ' age groups.'))
      } else{
        print('... All age bins now match.')
      }
    }
  }
  return(list(sim_df, bench_df))
}


get_available_sites_for_relationship = function(coord_csv, simulation_output_filepath, relationship_name, relationship_sim_filename){
  # determine which of the simulation sites both are indicated by the coordinator csv to be included in this validation relationship and have the relevant simulation output
  
  coord_sites = coord_csv$site[intersect(which(!is.na(coord_csv$site)), which(coord_csv[[relationship_name]]==1))]
  available_sites = c()
  for (ii in 1:length(coord_sites)){
    if (file.exists(paste0(simulation_output_filepath, '/', coord_sites[ii], '/', relationship_sim_filename))){
      available_sites = c(available_sites, coord_sites[ii])
    }
  }
  return(available_sites)
}



# if the maximum reference density bin is < (maximum simulation density bin / max_magnitude_difference), 
#    aggregate all simulation densities >= max ref bin into the max ref bin. The new final density bin will be all densities equal to or above that value
combine_higher_dens_freqs = function(sim_df_cur, max_ref_dens, max_magnitude_difference=100){
  if(max_ref_dens < (max(sim_df_cur$densitybin, na.rm=TRUE)/max_magnitude_difference)){
    # get sum of frequencies within higher bins
    all_higher_dens = sim_df_cur[sim_df_cur$densitybin >= max_ref_dens,]
    sim_agg_higher_dens = all_higher_dens %>% group_by(month, agebin, Site) %>%
      summarise(densitybin = min(densitybin),
                asexual_par_dens_freq = sum(asexual_par_dens_freq),
                gametocyte_dens_freq = sum(gametocyte_dens_freq),
                Pop = mean(Pop))
    # remove higher density bins from df
    cur_df_lower = sim_df_cur[sim_df_cur$densitybin < max_ref_dens,]
    # add back in the aggregated frequencies across higher density bins
    sim_df_cur = merge(cur_df_lower, sim_agg_higher_dens, all=TRUE)
  }
  return(sim_df_cur)
}





get_fraction_in_infectious_bin = function(sim_df){
  # translate the frequencies from simulation output into fraction of individuals in each 
  #     {age bin, month, densitybin, run number} group are in each infectiousness bin
  # return mean values across simulation seeds and within age groups (Assuming equal population sizes and weighting for all ages in an age bin)
  
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
  
  return(sim_df_agg2)
}









prepare_inc_df = function(coord_csv, simulation_output_filepath, base_reference_filepath, benchmark_simulation_filepath=NA){
  # read in, format, and align simulation and matched reference data into a single dataframe for incidence results across all sites
  
  
  # determine which of the age-incidence sites have the relevant simulation output
  available_sites = get_available_sites_for_relationship(coord_csv, simulation_output_filepath, relationship_name='age_incidence', relationship_sim_filename='inc_prev_data_final.csv')

  # iterate through sites, aggregating all age-incidence simulation data into one dataframe and reference data into a second dataframe (for the relevant sites)
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
  
  
  # scale down incidence in simulation according to probability of detecting a case in the reference setting
  sim_df$Incidence = sim_df$Incidence * sim_df$p_detect_case
  
  # set up reference and simulation dataset columns
  ref_df$Incidence = ref_df$INC / 1000
  ref_df$mean_age = (ref_df$INC_LAR + ref_df$INC_UAR)/2
  ref_df = data.frame('reference'=ref_df$Incidence, 'mean_age'=ref_df$mean_age, 'Site'=ref_df$Site, 'ref_pop_size'=ref_df$POP, 'ref_year'=ref_df$START_YEAR)
  sim_df = data.frame('simulation'=sim_df$Incidence, 'mean_age'=sim_df$mean_age, 'Site'=sim_df$Site)
  # format benchmark simulations
  if(nrow(bench_df)>0){
    # scale down incidence in simulation according to probability of detecting a case
    bench_df$Incidence = bench_df$Incidence * bench_df$p_detect_case
    bench_df = data.frame('benchmark'=bench_df$Incidence, 'mean_age'=bench_df$mean_age, 'Site'=bench_df$Site)
  }
  
  # check that ages match between reference and simulation. if there is a small difference (<1 year), update simulation to match reference age value after giving warning
  age_matched_dfs = match_sim_ref_ages(ref_df, sim_df, bench_df)
  sim_df = age_matched_dfs[[1]]
  bench_df = age_matched_dfs[[2]]
  
  combined_df = merge(ref_df, sim_df, all=TRUE)
  combined_df = merge(combined_df, bench_df, all=TRUE)
  combined_df$metric = 'incidence'

  return(combined_df)
}



# prepare dataframe with simulation and reference data formatted together
prepare_prev_df = function(coord_csv, simulation_output_filepath, base_reference_filepath, benchmark_simulation_filepath=NA){
  
  
  # determine which of the age-prevalence sites have the relevant simulation output
  available_sites = get_available_sites_for_relationship(coord_csv, simulation_output_filepath, relationship_name='age_prevalence', relationship_sim_filename='prev_inc_by_age_month.csv')

  
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
  
  
  # format reference data
  ref_df$Site = tolower(ref_df$Site)
  ref_df = data.frame('reference'=ref_df$prevalence, 'mean_age'=ref_df$mean_age, 'Site'=ref_df$Site, 'month'=ref_df$month,
                      'site_month'=paste0(ref_df$Site, '_month', ref_df$month), 
                      'total_sampled'=ref_df$total_sampled, 'num_pos'=ref_df$num_pos, 'ref_year' = ref_df$year)
  
  # format new simulation output
  # get simulation average across seeds
  sim_df = sim_df %>% group_by(Site, mean_age, month) %>%
    summarise(prevalence = mean(prevalence))
  sim_df$Site = tolower(sim_df$Site)
  sim_df = data.frame('simulation'=sim_df$prevalence, 'mean_age'=sim_df$mean_age, 'Site'=sim_df$Site, 'month'=sim_df$month,
                      'site_month'=paste0(sim_df$Site, '_month', sim_df$month))
    
  # format benchmark simulations
  if(nrow(bench_df)>0){
    # get simulation average across seeds
    bench_df = bench_df %>% group_by(Site, mean_age, month) %>%
      summarise(prevalence = mean(prevalence))
    bench_df$Site = tolower(bench_df$Site)
    bench_df = data.frame('benchmark'=bench_df$prevalence, 'mean_age'=bench_df$mean_age, 'Site'=bench_df$Site, 'month'=bench_df$month,
                          'site_month'=paste0(bench_df$Site, '_month', bench_df$month))
  }
  
  # check that ages match between reference and simulation. if there is a small difference (<1 year, update simulation)
  age_matched_dfs = match_sim_ref_ages(ref_df, sim_df, bench_df)
  sim_df = age_matched_dfs[[1]]
  bench_df = age_matched_dfs[[2]]
  
  combined_df = merge(ref_df, sim_df, all=TRUE)
  combined_df = merge(combined_df, bench_df, all=TRUE)
  combined_df$metric = 'prevalence'
  return(combined_df)

}






prepare_dens_df = function(coord_csv, simulation_output_filepath, base_reference_filepath, benchmark_simulation_filepath=NA){
  # read in, format, and align simulation and matched reference data into a single dataframe for parasite density results across all sites

  # determine which of the parasite density sites have the relevant simulation output
  available_sites = get_available_sites_for_relationship(coord_csv, simulation_output_filepath, relationship_name='age_parasite_density', relationship_sim_filename='parasite_densities_by_age_month.csv')
  
  # iterate through sites, grabbing relevant reference and simulation data to plot; combine data into a dataframe containing all sites
  sim_df = data.frame()
  bench_df = data.frame()
  ref_df = data.frame()
  for (ss in 1:length(available_sites)){
    cur_site = available_sites[ss]
    sim_df_cur = read.csv(paste0(simulation_output_filepath, '/', cur_site, '/parasite_densities_by_age_month.csv'))
    upper_ages = sort(unique(sim_df_cur$agebin))
    sim_df_cur$mean_age = sapply(sim_df_cur$agebin, get_mean_from_upper_age, upper_ages = upper_ages)
    age_agg_sim_df = get_age_bin_averages(sim_df_cur)
    
    filepath_ref = paste0(base_reference_filepath, '/', coord_csv$age_parasite_density_ref[which(coord_csv$site == cur_site)])
    ref_df_cur = read.csv(filepath_ref)
    ref_df_cur = ref_df_cur[tolower(ref_df_cur$Site) == tolower(cur_site),]
    ref_df_cur$Site = tolower(ref_df_cur$Site)
    upper_ages = sort(unique(ref_df_cur$agebin))
    ref_df_cur$mean_age = sapply(ref_df_cur$agebin, get_mean_from_upper_age, upper_ages = upper_ages)
    
    if(!is.na(benchmark_simulation_filepath)){
      bench_df_cur = read.csv(paste0(benchmark_simulation_filepath, '/', cur_site, '/parasite_densities_by_age_month.csv'))
      upper_ages = sort(unique(bench_df_cur$agebin))
      bench_df_cur$mean_age = sapply(bench_df_cur$agebin, get_mean_from_upper_age, upper_ages = upper_ages)
      age_agg_bench_df = get_age_bin_averages(sim_df=bench_df_cur)
      
      if(!all(sort(unique(age_agg_sim_df$mean_age)) == sort(unique(age_agg_bench_df$mean_age)))) {
        warning(paste0('New and benchmark simulation age bins are not the same for site:', cur_site, '. Removing benchmark sims.'))
        age_agg_bench_df = data.frame()
      }
      if(!all(sort(unique(age_agg_sim_df$densitybin)) == sort(unique(age_agg_bench_df$densitybin)))) {
        warning(paste0('New and benchmark simulation parasite density bins are not the same for site:', cur_site, '. Removing benchmark sims.'))
        age_agg_bench_df = data.frame()
      }
    } else age_agg_bench_df = data.frame()
    
    
    # subset simulation output to months in reference dataset
    months = sort(unique(ref_df_cur$month))
    sim_df_cur = age_agg_sim_df[age_agg_sim_df$month %in% months,]
    if(nrow(age_agg_bench_df)>0) {
      bench_df_cur = age_agg_bench_df[age_agg_bench_df$month %in% months,]
    } else bench_df_cur = data.frame()
    
    # if the maximum reference density bin is < (maximum simulation density bin / max_magnitude_difference), aggregate all simulation densities >= max ref bin into the max ref bin
    #   the final bin will be all densities equal to or above that value
    max_ref_dens = max(ref_df_cur$densitybin, na.rm=TRUE)
    sim_df_cur = combine_higher_dens_freqs(sim_df_cur, max_ref_dens, max_magnitude_difference=100)
    if(nrow(bench_df_cur)>0) bench_df_cur = combine_higher_dens_freqs(bench_df_cur, max_ref_dens, max_magnitude_difference=100)

    
    # add zeros for unobserved reference densities up to max_ref_dens
    all_zeros_df = sim_df_cur[,c('month', 'mean_age', 'agebin', 'densitybin', 'Site')]
    ref_df_cur = merge(ref_df_cur, all_zeros_df, all=TRUE)
    ref_df_cur[is.na(ref_df_cur)] = 0

    # add site into larger dataframe
    if(nrow(ref_df)>0){
      ref_df = merge(ref_df, ref_df_cur, all=TRUE)
    } else ref_df = ref_df_cur
    sim_df = rbind(sim_df, sim_df_cur)
    bench_df = rbind(bench_df, bench_df_cur)
  }

  # check that ages match between reference and simulation. if there is a small difference (<1 year, update simulation)
  age_matched_dfs = match_sim_ref_ages(ref_df, sim_df, bench_df)
  sim_df = age_matched_dfs[[1]]
  bench_df = age_matched_dfs[[2]]  
  
  # format reference data
  ref_df_asex = data.frame('reference'=ref_df$asexual_par_dens_freq, 'mean_age'=ref_df$mean_age, 'agebin'=ref_df$agebin, 'densitybin'=ref_df$densitybin, 'Site'=ref_df$Site, 'month'=ref_df$month,
                      'site_month'=paste0(ref_df$Site, '_month', ref_df$month), 
                      'ref_total'=ref_df$bin_total_asex, 'ref_bin_count'=ref_df$count_asex)
  ref_df_gamet = data.frame('reference'=ref_df$gametocyte_dens_freq, 'mean_age'=ref_df$mean_age, 'agebin'=ref_df$agebin, 'densitybin'=ref_df$densitybin, 'Site'=ref_df$Site, 'month'=ref_df$month,
                           'site_month'=paste0(ref_df$Site, '_month', ref_df$month), 
                           'ref_total'=ref_df$bin_total_gamet, 'ref_bin_count'=ref_df$count_gamet)
  
  # format new simulation output
  sim_df_asex = data.frame('simulation'=sim_df$asexual_par_dens_freq, 'mean_age'=sim_df$mean_age, 'agebin'=sim_df$agebin, 'densitybin'=sim_df$densitybin, 'Site'=sim_df$Site, 'month'=sim_df$month,
                      'site_month'=paste0(sim_df$Site, '_month', sim_df$month))
  sim_df_gamet = data.frame('simulation'=sim_df$gametocyte_dens_freq, 'mean_age'=sim_df$mean_age, 'agebin'=sim_df$agebin, 'densitybin'=sim_df$densitybin, 'Site'=sim_df$Site, 'month'=sim_df$month,
                           'site_month'=paste0(sim_df$Site, '_month', sim_df$month))
  
  # format benchmark simulations
  if(nrow(bench_df)>0){
    bench_df_asex = data.frame('benchmark'=bench_df$asexual_par_dens_freq, 'mean_age'=bench_df$mean_age, 'agebin'=bench_df$agebin, 'densitybin'=bench_df$densitybin, 'Site'=bench_df$Site, 'month'=bench_df$month,
                          'site_month'=paste0(bench_df$Site, '_month', bench_df$month))
    bench_df_gamet = data.frame('benchmark'=bench_df$gametocyte_dens_freq, 'mean_age'=bench_df$mean_age, 'agebin'=bench_df$agebin, 'densitybin'=bench_df$densitybin, 'Site'=bench_df$Site, 'month'=bench_df$month,
                               'site_month'=paste0(bench_df$Site, '_month', bench_df$month))
  }
  

  # combine reference and simulation dataframes
  combined_df_asex = merge(ref_df_asex, sim_df_asex, all=TRUE)
  combined_df_asex = merge(combined_df_asex, bench_df_asex, all=TRUE)
  combined_df_asex$metric = 'asexual_density'
  
  combined_df_gamet = merge(ref_df_gamet, sim_df_gamet, all=TRUE)
  combined_df_gamet = merge(combined_df_gamet, bench_df_gamet, all=TRUE)
  combined_df_gamet$metric = 'gametocyte_density'
  
  return(list(combined_df_asex, combined_df_gamet))
  
}






prepare_infect_df = function(coord_csv, simulation_output_filepath, base_reference_filepath, benchmark_simulation_filepath=NA){
  
  # determine which of the infectiousness sites have the relevant simulation output
  available_sites = get_available_sites_for_relationship(coord_csv, simulation_output_filepath, relationship_name='infectiousness_to_mosquitos', relationship_sim_filename='infectiousness_by_age_density_month.csv')
  
  # iterate through sites, grabbing relevant reference and simulation data to plot; combine data into a dataframe containing all sites
  sim_df = data.frame()
  bench_df = data.frame()
  ref_df = data.frame()
  for (ss in 1:length(available_sites)){
    cur_site = available_sites[ss]
    filepath_ref = paste0(base_reference_filepath, '/', coord_csv$infectiousness_to_mosquitos_ref[which(coord_csv$site == cur_site)])
    ref_df_cur = read.csv(filepath_ref)
    ref_df_cur = ref_df_cur[tolower(ref_df_cur$site) == tolower(cur_site),]
    ref_months = unique(ref_df_cur$month)
    
    sim_df_cur = read.csv(paste0(simulation_output_filepath, '/', cur_site, '/infectiousness_by_age_density_month.csv'))
    # remove simulation rows with zero pop
    sim_df_cur = sim_df_cur[sim_df_cur$Pop>0,]    
    # subset simulation to months in reference df
    sim_df_cur = sim_df_cur[sim_df_cur$month %in% ref_months,]
    # get mean (across simulation seeds and ages within a bin) fraction of individuals in each  {age bin, month, densitybin, run number} group that fall in each infectiousness bin
    sim_df_agg2 = get_fraction_in_infectious_bin(sim_df_cur)
    
    if(!is.na(benchmark_simulation_filepath)){
      bench_df_cur = read.csv(paste0(benchmark_simulation_filepath, '/', cur_site, '/infectiousness_by_age_density_month.csv'))
      # remove simulation rows with zero pop
      bench_df_cur = bench_df_cur[bench_df_cur$Pop>0,]    
      # subset simulation to months in reference df
      bench_df_cur = bench_df_cur[bench_df_cur$month %in% ref_months,]
      # get mean (across simulation seeds and ages within a bin) fraction of individuals in each  {age bin, month, densitybin, run number} group that fall in each infectiousness bin
      bench_df_agg2 = get_fraction_in_infectious_bin(bench_df_cur)
      
      if(!all(sort(unique(sim_df_agg2$agebin)) == sort(unique(bench_df_agg2$agebin)))) {
        warning(paste0('New and benchmark simulation age bins are not the same for site:', cur_site, '. Removing benchmark sims.'))
        bench_df_agg2 = data.frame()
      }
      if(!all(sort(unique(sim_df_agg2$densitybin)) == sort(unique(bench_df_agg2$densitybin)))) {
        warning(paste0('New and benchmark simulation parasite density bins are not the same for site:', cur_site, '. Removing benchmark sims.'))
        bench_df_agg2 = data.frame()
      }
    } else bench_df_agg2 = data.frame()
    
    
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
    # standardize column names and merge simulation and reference data frames
    colnames(ref_df_cur)[colnames(ref_df_cur) == 'site'] = 'Site'
    ref_df_cur = data.frame('reference'=ref_df_cur$freq_frac_infect, 'agebin'=ref_df_cur$agebin, 'densitybin'=ref_df_cur$densitybin, 'fraction_infected_bin'=ref_df_cur$fraction_infected_bin,
                        'Site'=ref_df_cur$Site, 'month'=ref_df_cur$month, 'site_month'=paste0(ref_df_cur$Site, '_month', ref_df_cur$month), 
                        'ref_total'=ref_df_cur$num_in_group, 'ref_bin_count'=ref_df_cur$count)
    sim_df_cur = data.frame('simulation'=sim_df_agg2$infectiousness_bin_freq, 'agebin'=sim_df_agg2$agebin, 'densitybin'=sim_df_agg2$densitybin, 'fraction_infected_bin'=sim_df_agg2$infectiousness_bin,
                        'Site'=sim_df_agg2$Site, 'month'=sim_df_agg2$month,
                        'site_month'=paste0(sim_df_agg2$Site, '_month', sim_df_agg2$month))

    if(nrow(bench_df_agg2)>0){
      bench_df_cur = data.frame('benchmark'=bench_df_agg2$infectiousness_bin_freq, 'agebin'=bench_df_agg2$agebin, 'densitybin'=bench_df_agg2$densitybin, 'fraction_infected_bin'=bench_df_agg2$infectiousness_bin,
                          'Site'=bench_df_agg2$Site, 'month'=bench_df_agg2$month,
                          'site_month'=paste0(bench_df_agg2$Site, '_month', bench_df_agg2$month))
    }
    # add site into larger dataframe
    sim_df = rbind(sim_df, sim_df_cur)
    ref_df = rbind(ref_df, ref_df_cur)
    bench_df = rbind(bench_df, bench_df_cur)
  }
  combined_df = merge(sim_df, ref_df, all=TRUE)
  if(nrow(bench_df)>0) combined_df = merge(combined_df, bench_df, all=TRUE)
  combined_df$metric = 'infectiousness'
  return(combined_df)
}













# infection duration functions

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# get subset of simulation dataset to match reference dataset for dates and ages of individuals of sampled individuals
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
get_sim_survey = function(sim_dir, ref_df, seeds=NA){
  sampled_sim_filename = 'sim_duration_survey_sampling.csv'
  if(file.exists(paste0(sim_dir, '/', sampled_sim_filename))){
    sim_subset_full = read.csv(paste0(sim_dir, '/', sampled_sim_filename))
  } else{
    # get first year of sampling in reference dataset. the simulation will be referenced from the first day of that year
    first_ref_date = as.Date(paste0(year(min(ref_df$date, na.rm=TRUE)), '-01-01'))
    indIDs = unique(ref_df$SID)
    
    sim_full = fread(paste0(sim_dir, '/patient_reports.csv'))
    sim_full$date = first_ref_date + sim_full$simday
    sim_full$age = sim_full$age/365
    
    if(is.na(seeds)){
      seeds = unique(sim_full$Run_Number)
    }
    for(seed in sort(seeds)){
      print(paste0('Currently on seed ', seed))
      sim = sim_full[which(sim_full$Run_Number == seed),]  # subset to desired run
      included_ids = c() # track which individuals have already been included from the simulation (to avoid double-sampling simulation individuals)
      sim_subset = data.table()
      for(ii in 1:length(indIDs)){
        if (ii %% 50 == 0) print(paste0('Currently on individual ', ii, ' out of ', length(indIDs)))
        ref_df_cur = ref_df[ref_df$SID == indIDs[ii],]
        ref_df_cur = ref_df_cur[order(ref_df_cur$date),]
        # find a matching individual
        age_cur = ref_df_cur$age[1]
        day_cur = ref_df_cur$date[1]
        
        # use age-specific matches
        id_candidates = sim$id[intersect(which(sim$date == day_cur), which(round(sim$age) ==round(age_cur)))]
        id_candidates = id_candidates[!(id_candidates %in% included_ids)]
        # if no perfect age-match remain, expand year-range until finding a match
        if(length(id_candidates) == 0){
          year_range = 0
          while(length(id_candidates) == 0 & year_range < 100){
            year_range = year_range + 5
            id_candidates = sim$id[intersect(which(sim$date == day_cur), which(round(sim$age) %in% seq((round(age_cur)-year_range), (round(age_cur)+year_range))))]
            id_candidates = id_candidates[!(id_candidates %in% included_ids)]
          }
          if(length(id_candidates) == 0){
            print(paste0('Problem: no age-matched simulation individual found for reference id: ', indIDs[ii]))
          } else{
            print(paste0('No exact age match remaining for reference id: ', indIDs[ii], '. Used simulation individual within ', year_range, ' years.'))
          }
        }
        
        id_sim_cur = sample(id_candidates, size=1)
        included_ids = c(included_ids, id_sim_cur)
        
        # keep the same simulation dates as the reference samples for this individual
        sim_subset_cur = sim[intersect(which(sim$id == id_sim_cur), which(sim$date %in% ref_df_cur$date)),]
        sim_subset = rbind(sim_subset, sim_subset_cur)
      }
      sim_subset$seed = seed
      if (seed == sort(seeds)[1]){
        sim_subset_full = sim_subset
      } else{
        sim_subset_full = rbind(sim_subset_full, sim_subset)
      }
    }
    
    # rename simulation columns to match reference data
    colnames(sim_subset_full)[colnames(sim_subset_full) == 'id'] = 'SID'
    colnames(sim_subset_full)[colnames(sim_subset_full) == 'true_asexual_parasites'] = 'DENSITY'
    
    write.csv(sim_subset_full, paste0(sim_dir, '/', sampled_sim_filename), row.names=FALSE)
  }
  return(sim_subset_full)
}







