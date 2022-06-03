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


# infection duration helper function

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












library(broom)
library(tidyverse)


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


# prepare dataframe with simulation and reference data formatted together
prepare_inc_df = function(sim_df, ref_df, bench_df = data.frame()){
  # scale down incidence in simulation according to probability of detecting a case
  sim_df$Incidence = sim_df$Incidence * sim_df$p_detect_case
  
  # set up reference and simulation dataset columns
  ref_df$Incidence = ref_df$INC / 1000
  ref_df$mean_age = (ref_df$INC_LAR + ref_df$INC_UAR)/2
  ref_df = data.frame('reference_value'=ref_df$Incidence, 'mean_age'=ref_df$mean_age, 'Site'=ref_df$Site) #, 'ref_pop_size'=ref_df$POP, 'ref_year'=ref_df$START_YEAR)
  sim_df = data.frame('simulation_value'=sim_df$Incidence, 'mean_age'=sim_df$mean_age, 'Site'=sim_df$Site)
  
  # format benchmark simulations
  if(nrow(bench_df)>0){
    # scale down incidence in simulation according to probability of detecting a case
    bench_df$Incidence = bench_df$Incidence * bench_df$p_detect_case
    bench_df = data.frame('benchmark_value'=bench_df$Incidence, 'mean_age'=bench_df$mean_age, 'Site'=bench_df$Site)
  }
  # check that ages match between reference and simulation. if there is a small difference (<1 year, update simulation)
  age_matched_dfs = match_sim_ref_ages(ref_df, sim_df, bench_df)
  sim_df = age_matched_dfs[[1]]
  bench_df = age_matched_dfs[[2]]
  
  combined_df = merge(ref_df, sim_df, all=TRUE)
  combined_df = merge(combined_df, bench_df, all=TRUE)
  combined_df$metric = 'incidence'
  return(combined_df)
}



# prepare dataframe with simulation and reference data formatted together
prepare_prev_df = function(sim_df, ref_df, bench_df=data.frame()){
  # get simulation average across seeds
  sim_df = sim_df %>% group_by(Site, mean_age, month) %>%
    summarise(prevalence = mean(prevalence))
  
  ref_df$Site = tolower(ref_df$Site)
  sim_df$Site = tolower(sim_df$Site)
  
  ref_df = data.frame('reference_value'=ref_df$prevalence, 'mean_age'=ref_df$mean_age, 'Site'=ref_df$Site, 'month'=ref_df$month,
                      'site_month'=paste0(ref_df$Site, '_month', ref_df$month), 'total_sampled'=ref_df$total_sampled, 'num_pos'=ref_df$num_pos)
  sim_df = data.frame('simulation_value'=sim_df$prevalence, 'mean_age'=sim_df$mean_age, 'Site'=sim_df$Site, 'month'=sim_df$month,
                      'site_month'=paste0(sim_df$Site, '_month', sim_df$month))
  
  # format benchmark simulations
  if(nrow(bench_df)>0){
    # get simulation average across seeds
    bench_df = bench_df %>% group_by(Site, mean_age, month) %>%
      summarise(prevalence = mean(prevalence))
    bench_df$Site = tolower(bench_df$Site)
    bench_df = data.frame('benchmark_value'=bench_df$prevalence, 'mean_age'=bench_df$mean_age, 'Site'=bench_df$Site, 'month'=bench_df$month,
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






