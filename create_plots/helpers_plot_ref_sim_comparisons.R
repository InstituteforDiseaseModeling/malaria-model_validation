# helpers_plot_ref_sim_comparisons.R

# These functions generate plots comparing between reference data and simulated outputs (and also between new and benchmark simulations). 



library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)
library(ggpubr)
library(stringr)
library(data.table)
library(RcmdrMisc)
library(lubridate)





########################## plot age-incidence comparisons with reference ####################
plot_inc_ref_sim_comparison = function(sim_df, ref_df, bench_df=data.frame()){
  # scale down incidence in simulation according to probability of detecting a case
  sim_df$Incidence = sim_df$Incidence * sim_df$p_detect_case
  
  # set up reference and simulation dataset columns
  ref_df$Incidence = ref_df$INC / 1000
  ref_df$mean_age = (ref_df$INC_LAR + ref_df$INC_UAR)/2
  ref_df = data.frame('Incidence'=ref_df$Incidence, 'mean_age'=ref_df$mean_age, 'Site'=ref_df$Site, 'Pop_size'=ref_df$POP, 'year'=ref_df$START_YEAR)
  sim_df = data.frame('Incidence'=sim_df$Incidence, 'mean_age'=sim_df$mean_age, 'Site'=sim_df$Site, 'Pop_size'=NA, 'year'=NA)
  ref_df$source = 'reference'
  sim_df$source = 'simulation'
  
  if(nrow(bench_df)>0){
    # scale down incidence in simulation according to probability of detecting a case
    bench_df$Incidence = bench_df$Incidence * bench_df$p_detect_case
    bench_df = data.frame('Incidence'=bench_df$Incidence, 'mean_age'=bench_df$mean_age, 'Site'=bench_df$Site, 'Pop_size'=NA, 'year'=NA)
    bench_df$source = 'benchmark'
  }
  
  df_combined = rbind(sim_df, ref_df, bench_df)
  df_combined$source = factor(df_combined$source, levels=c('reference', 'simulation', 'benchmark'))
  
  gg = ggplot(df_combined, aes(x=mean_age, y=Incidence, color=source, shape=source, group=year)) +
    geom_line() +
    geom_point(aes(size=source)) +
    scale_color_manual(values = c("reference" = rgb(169/255,23/255,23/255, alpha=0.8),
                                  "simulation"=rgb(0/255,124/255,180/255, alpha=0.8),
                                  "benchmark"=rgb(0,0,0,0.8))) +
    scale_shape_manual(values=c("reference" = 16,
                                "simulation"=16,
                                "benchmark"=1)) +
    scale_size_manual(values=c("reference" = 1.5,
                               "simulation"=1.5,
                               "benchmark"=2)) +
    xlab('age (midpoint of age bin)') +
    ylab('incidence') +
    facet_wrap('Site', ncol=4) +
    theme_bw()  
  
  return(gg)
}






########################## plot age-prevalence comparisons with reference without sweep background ####################
plot_prev_ref_sim_comparison = function(sim_df, ref_df, bench_df=data.frame()){
  # get simulation average across seeds
  sim_df = sim_df %>% group_by(Site, mean_age, month) %>%
    summarise(prev_sd = sd(prevalence),
              prevalence = mean(prevalence))
  
  ref_df$Site = tolower(ref_df$Site)
  sim_df$Site = tolower(sim_df$Site)
  ref_df$source = 'reference'
  sim_df$source = 'simulation'
  df_combined = merge(sim_df, ref_df, all=TRUE)
  
  if(nrow(bench_df)>0){
    # get simulation average across seeds
    bench_df = bench_df %>% group_by(Site, mean_age, month) %>%
      summarise(prev_sd = sd(prevalence),
                prevalence = mean(prevalence))
    bench_df$Site = tolower(bench_df$Site)
    bench_df$source = 'benchmark'
    df_combined = merge(df_combined, bench_df, all=TRUE)
  }
  df_combined$source = factor(df_combined$source, levels=c('reference', 'simulation', 'benchmark'))
  df_combined$site_month = paste0(df_combined$Site, '_month', df_combined$month)
  
  gg = ggplot(df_combined, aes(x=mean_age, y=prevalence, color=source, shape=source, group=year)) +
    geom_line() +
    geom_point(aes(size=source)) +
    scale_color_manual(values = c("reference" = rgb(169/255,23/255,23/255, alpha=0.8),
                                  "simulation"=rgb(0/255,124/255,180/255, alpha=0.8),
                                  "benchmark"=rgb(0,0,0,0.8))) +
    scale_shape_manual(values=c("reference" = 16,
                                "simulation"=16,
                                "benchmark"=1)) +
    scale_size_manual(values=c("reference" = 1.5,
                               "simulation"=1.5,
                               "benchmark"=2)) +
    xlab('age (midpoint of age bin)') +
    ylab('prevalence') +
    facet_wrap('site_month', ncol=4) +
    theme_bw()  
  
  return(gg)
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












# create scatter plots with new versus benchmark simulation output and get comparison of likelihoods for each simulation set and each site.
compare_benchmark = function(combined_df){
  combined_df = combined_df[!is.na(combined_df$benchmark_value),]
  metric = combined_df$metric[1]
  if('site_month' %in% colnames(combined_df)) combined_df$Site = combined_df$site_month
  min_value = min(c(combined_df$benchmark_value, combined_df$simulation_value), na.rm=TRUE)
  max_value = max(c(combined_df$benchmark_value, combined_df$simulation_value), na.rm=TRUE)
  gg = ggplot(combined_df, aes(x=benchmark_value, y=simulation_value, color=Site, fill=Site))+
    geom_point(size=2) + 
    ylab(paste0('new simulation ', metric)) + 
    xlab(paste0('benchmark sim ', metric)) + 
    ggtitle('Benchmark versus new sim values') +
    coord_fixed(ratio=1, xlim=c(min_value, max_value), ylim=c(min_value, max_value))+
    geom_abline(slope=1,intercept=0, color='grey', alpha=0.5) + 
    theme_classic() +
    theme(plot.title = element_text(size=12))
  return(gg)
}












############################# infection duration functions ################################33



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# calculate probability of going from negative to positive or from positive to negative between sample dates
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
get_frac_state_swaps = function(data){
  
  # brute force approach iterating through people and dates
  indIDs = unique(data$SID)
  sum_denom_pos = 0
  sum_turn_neg = 0
  sum_denom_neg = 0
  sum_turn_pos = 0
  
  for(ii in 1:length(indIDs)){
    data_cur = data[data$SID == indIDs[ii],]
    data_cur = data_cur[order(data_cur$date),]
    # get indices of positive tests and of negative tests
    ind_pos = which(data_cur$DENSITY > pos_thresh_dens)
    ind_neg = which(data_cur$DENSITY <= pos_thresh_dens)
    
    # denominators for each (number of each type that had an observation after then (i.e., they could have been observed to change))
    last_obs_pos = (data_cur$DENSITY[nrow(data_cur)] > pos_thresh_dens)
    sum_denom_pos = sum_denom_pos + ifelse(last_obs_pos, length(ind_pos)-1, length(ind_pos))
    sum_denom_neg = sum_denom_neg + ifelse(last_obs_pos, length(ind_neg), length(ind_neg)-1)
    
    # find how many tests change from neg to pos or pos to neg across timesteps
    sum_turn_neg = sum_turn_neg + sum((ind_pos + 1) %in% ind_neg)
    sum_turn_pos = sum_turn_pos + sum((ind_neg + 1) %in% ind_pos)
  }
  
  frac_pos_turn_neg_next_time = sum_turn_neg / sum_denom_pos
  frac_neg_turn_pos_next_time = sum_turn_pos / sum_denom_neg
  
  return(c(frac_pos_turn_neg_next_time, frac_neg_turn_pos_next_time))
}


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# calculate infection duration (time between consequtive positive samples) and whether or not the infection observation was censored
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
get_time_pos = function(data){
  
  # brute force approach iterating through people and dates
  indIDs = unique(data$SID)
  days_positive = c()  # number of days between the first and last sample date in a run of positive samples
  sample_censored = c()  # Boolean indicating whether the positivity duration might have been longer unobserved
  sample_ages = c()
  sample_seeds = c()
  cur_index = 1
  
  if(!('seed' %in% colnames(data))) data$seed = 0
  for (ss in unique(data$seed)){
    data_ss = data[data$seed == ss,]
    for(ii in 1:length(indIDs)){
      data_cur = data_ss[data_ss$SID == indIDs[ii],]
      data_cur = data_cur[order(data_cur$date),]
      cur_age = mean(data_cur$age)
      # which samples had positive tests
      pos_bool = (data_cur$DENSITY > pos_thresh_dens)
      # get the length of spans of positive and of negative tests in a row
      num_in_a_row = rle(pos_bool)
      start_index_num_in_a_row = c(1,1+cumsum(num_in_a_row$lengths))[-(1+length(num_in_a_row$lengths))]
      
      # for each span of positives, determine the length before turning negative as well as whether it was censored or not
      # if it includes the first or last sample, it's censored (it's at least that many days but may have been longer)
      num_in_a_row$lengths[num_in_a_row$values]
      
      # iterate through periods of positivity
      for(tt in which(num_in_a_row$values)){
        # indicate whether this period of positivity was censored
        sample_censored[cur_index] = (tt==1 | tt==length(num_in_a_row$lengths))
        sample_ages[cur_index] = cur_age
        sample_seeds[cur_index] = ss
        
        # get the first and last positive sample day observed for this period of positivity
        first_index = start_index_num_in_a_row[tt]
        last_index = first_index + num_in_a_row$lengths[tt] - 1
        first_positive_sample_day = data_cur$date[first_index]
        last_positive_sample_day = data_cur$date[last_index]
        
        days_positive[cur_index] = as.numeric(difftime(last_positive_sample_day, first_positive_sample_day, units = c("days")))
        
        cur_index = cur_index + 1
      }
    }
  }
  
  return(data.frame('days_positive' = days_positive, 'censored' = sample_censored, 'age' = sample_ages, 'seed' = sample_seeds))
}





# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# functions to get binned durations of infections, with quantile ranges across sim seeds, optionally faceted by censorship and age
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
bin_durations_all_seeds = function(days_positive, duration_bins){
  # create data frame of binned days positive for each seed in days_positive
  bin_w_seeds = data.frame(bin_mid=(duration_bins[-length(duration_bins)]+duration_bins[-1])/2,
                           bin_min=duration_bins[-length(duration_bins)],
                           bin_max=duration_bins[-1])
  for(seed in unique(days_positive$seed)){
    days_positive_cur = days_positive[days_positive$seed == seed,]
    binned_counts = hist(x=days_positive_cur$days_positive, breaks=duration_bins, plot=FALSE)
    binned_dens = binned_counts$counts / sum(binned_counts$counts)
    bin_w_seeds[[paste0('seed_', seed)]] = as.numeric(binned_dens)
  }
  
  if(length(which(grepl('seed_', colnames(bin_w_seeds))))>1){
    bin_w_seeds$density = apply(bin_w_seeds[,which(grepl('seed_', colnames(bin_w_seeds)))], 1, mean)
    bin_w_seeds$quant_low = apply(bin_w_seeds[,which(grepl('seed_', colnames(bin_w_seeds)))], 1, quantile, probs=0.0)  # 0.1
    bin_w_seeds$quant_high = apply(bin_w_seeds[,which(grepl('seed_', colnames(bin_w_seeds)))], 1, quantile, probs=1)  # 0.9
  } else if(length(which(grepl('seed_', colnames(bin_w_seeds)))) == 1){
    bin_w_seeds$density = bin_w_seeds[,which(grepl('seed_', colnames(bin_w_seeds)))]
    bin_w_seeds$quant_low = NA
    bin_w_seeds$quant_high = NA
  } else{ # no data present for this binning
    bin_w_seeds$density = NA
    bin_w_seeds$quant_low = NA
    bin_w_seeds$quant_high = NA
  }
  
  bin_df = bin_w_seeds[,c('bin_mid', 'density', 'quant_low', 'quant_high', 'bin_min', 'bin_max')]
  return(bin_df)
}


get_age_group = function(cur_age, age_bin_lower, age_bin_labels){
  return(age_bin_labels[max(which(cur_age>age_bin_lower))])
}


get_duration_bins = function(ind_data, duration_bins, facet_censored=TRUE, facet_age=TRUE, age_bin_lower=c(0, 5, 10, 20, 100)){
  if(facet_age) age_bin_labels = c(paste0(age_bin_lower[1:(length(age_bin_lower)-1)], '-', age_bin_lower[2:length(age_bin_lower)]), paste0('>', age_bin_lower[length(age_bin_lower)]))
  
  # get data frame of the days each observed infection lasted, also recording age, whether infection was censored, and seed
  days_positive = get_time_pos(data=ind_data)
  if(facet_age){
    days_positive$age_group = sapply(days_positive$age, get_age_group, age_bin_lower=age_bin_lower, age_bin_labels=age_bin_labels)
  }
  
  # determine how dataset should be divided up among factors/facets
  bin_df = data.frame()
  if(facet_censored & facet_age){
    for(v1 in unique(days_positive$censored)){
      for(v2 in unique(days_positive$age_group)){
        bins_cur = bin_durations_all_seeds(days_positive[intersect(which(days_positive$censored == v1), which(days_positive$age_group == v2)),], duration_bins=duration_bins)
        bins_cur$censored = v1
        bins_cur$age_group = v2
        bin_df = rbind(bin_df, bins_cur)
      }
    }
  } else if(facet_censored){
    for(v1 in unique(days_positive$censored)){
      bins_cur = bin_durations_all_seeds(days_positive[which(days_positive$censored == v1),], duration_bins=duration_bins)
      bins_cur$censored = v1
      bins_cur$age_group = 'combined'
      bin_df = rbind(bin_df, bins_cur)
    }
  } else if(facet_age){
    for(v2 in unique(days_positive$age_group)){
      bins_cur = bin_durations_all_seeds(days_positive[which(days_positive$age_group == v2),], duration_bins=duration_bins)
      bins_cur$censored = 'combined'
      bins_cur$age_group = v2
      bin_df = rbind(bin_df, bins_cur)
    }
  } else{
    bins_cur = bin_durations_all_seeds(days_positive, duration_bins=duration_bins)
    bins_cur$censored = 'combined'
    bins_cur$age_group = 'combined'
    bin_df = bins_cur
  }
  
  if(facet_age){
    bin_df$age_group = factor(bin_df$age_group, levels=age_bin_labels)
  }
  
  return(bin_df)
}



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# create barplot comparing sim and ref for fraction of samples positive and fractions of samples switching from neg--> pos or pos--> neg
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
create_barplot_frac_comparison = function(ref_df, sim_data, pos_thresh_dens, sim_dir){
  frac_swap_ref = get_frac_state_swaps(data=ref_df)
  frac_swap_sim = get_frac_state_swaps(data=sim_data)
  
  ref_df$pos = ref_df$DENSITY > pos_thresh_dens
  frac_samples_pos_ref = sum(ref_df$pos, na.rm=TRUE) / sum((ref_df$DENSITY > -1), na.rm=TRUE)
  
  sim_data$pos = sim_data$DENSITY > pos_thresh_dens
  frac_samples_pos_sim = sum(sim_data$pos, na.rm=TRUE) / sum((sim_data$DENSITY > -1), na.rm=TRUE)
  
  df = data.frame(source=c(rep('reference',3), rep('simulation', 3)), 
                  measure=rep(c('negative to positive next time', 'positive to negative next time', 'total samples positive'),2), 
                  value=c(frac_swap_ref[2], frac_swap_ref[1], frac_samples_pos_ref, frac_swap_sim[2], frac_swap_sim[1], frac_samples_pos_sim)) 
  
  gg = ggplot(df, aes(x=measure, y=value, fill=source, color=source)) + 
    geom_bar(stat="identity", position = "dodge", alpha=.3) +
    scale_color_manual(values = c("reference" = rgb(169/255,23/255,23/255, alpha=0.8),
                                  "simulation"=rgb(0/255,124/255,180/255, alpha=0.8))) +
    scale_fill_manual(values = c("reference" = rgb(169/255,23/255,23/255, alpha=0.8),
                                 "simulation"=rgb(0/255,124/255,180/255, alpha=0.8))) +
    labs(title='general dataset properties', y='fraction of samples', x='')
  
  return(gg)
}

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# create plots comparing the distributions of infection lengths
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

plot_infection_duration_dist = function(ref_df, sim_data, pos_thresh_dens, sim_dir, duration_bins=c(seq(0,350,50), 500)){
  
  # get densities for each duration bin
  ref_bin_df = get_duration_bins(ind_data=ref_df, duration_bins=duration_bins, facet_censored=TRUE, facet_age=FALSE)
  ref_bin_df$dataset = 'reference'
  ref_bin_df$censor_type = 'start & finish observed'
  ref_bin_df$censor_type[ref_bin_df$censored] = 'censored'
  
  # combine multiple simulation seeds to show mean, 10%, 90% quantile values across seed distributions
  sim_bin_df = get_duration_bins(ind_data=sim_data, duration_bins=duration_bins, facet_censored=TRUE, facet_age=FALSE)
  sim_bin_df$dataset = 'simulation'
  sim_bin_df$censor_type = 'start & finish observed'
  sim_bin_df$censor_type[sim_bin_df$censored] = 'censored'
  
  # combine reference and simulation datasets
  bin_df = rbind(ref_bin_df, sim_bin_df)
  
  # max_y=0.8
  gg = ggplot(bin_df, aes(x=bin_mid, y=density, color=dataset, fill=dataset)) +
    geom_bar(stat="identity",position = "identity", alpha=.3)+
    geom_errorbar(aes(ymin=quant_low, ymax=quant_high), width=.2,
                  position=position_dodge(.9)) +
    scale_color_manual(values = c("reference" = rgb(169/255,23/255,23/255, alpha=0.8),
                                  "simulation"=rgb(0/255,124/255,180/255, alpha=0.8))) +
    scale_fill_manual(values = c("reference" = rgb(169/255,23/255,23/255, alpha=0.8),
                                 "simulation"=rgb(0/255,124/255,180/255, alpha=0.8))) +
    # ylim(NA,max_y) +
    labs(title=paste0('infection duration'), x='infection duration (days)') + 
    facet_wrap(scales='fixed', facets='censor_type')
  return(gg)
  # ggsave(filename=paste0(sim_dir, '/duration_comparison_xEIR', xEIR_val, '_densThresh', pos_thresh_dens, '.png'), plot=gg, width=6, height=3, units='in')
}



plot_infection_duration_dist_ref_only = function(ref_df, sim_dir, duration_bins=c(seq(0,350,50), 500)){
  
  # get densities for each duration bin
  ref_bin_df = get_duration_bins(ind_data=ref_df, duration_bins=duration_bins, facet_censored=TRUE, facet_age=FALSE)
  ref_bin_df$dataset = 'reference'
  ref_bin_df$censor_type = 'start & finish observed'
  ref_bin_df$censor_type[ref_bin_df$censored] = 'censored'
  bin_df = ref_bin_df
  
  # max_y=0.8
  gg = ggplot(bin_df, aes(x=bin_mid, y=density, color=dataset, fill=dataset)) +
    geom_bar(stat="identity",position = "identity", alpha=.3)+
    scale_color_manual(values = c("reference" = rgb(169/255,23/255,23/255, alpha=0.8),
                                  "simulation"=rgb(0/255,124/255,180/255, alpha=0.8))) +
    scale_fill_manual(values = c("reference" = rgb(169/255,23/255,23/255, alpha=0.8),
                                 "simulation"=rgb(0/255,124/255,180/255, alpha=0.8))) +
    # ylim(NA,max_y) +
    labs(title=paste0('infection duration - reference'), x='infection duration (days)') + 
    facet_wrap(scales='fixed', facets='censor_type')
  
  ggsave(filename=paste0(sim_dir, '/duration_allAges_reference.png'), plot=gg, width=6, height=3, units='in')
}



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# create plots comparing the distributions of infection lengths, faceted by age group
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

plot_infection_duration_dist_by_age = function(ref_df, sim_data, pos_thresh_dens, sim_dir, age_bin_lower=c(0, 5, 10, 20, 100), duration_bins=c(seq(0,350,50), 500)){
  
  # get densities for each duration bin
  ref_bin_df = get_duration_bins(ind_data=ref_df, duration_bins=duration_bins, facet_censored=TRUE, facet_age=TRUE)
  ref_bin_df$dataset = 'reference'
  ref_bin_df$censor_type = 'start & finish observed'
  ref_bin_df$censor_type[ref_bin_df$censored] = 'censored'
  
  # combine multiple simulation seeds to show mean, 10%, 90% quantile values across seed distributions
  sim_bin_df = get_duration_bins(ind_data=sim_data, duration_bins=duration_bins, facet_censored=TRUE, facet_age=TRUE)
  sim_bin_df$dataset = 'simulation'
  sim_bin_df$censor_type = 'start & finish observed'
  sim_bin_df$censor_type[sim_bin_df$censored] = 'censored'
  
  # combine reference and simulation datasets
  bin_df = rbind(ref_bin_df, sim_bin_df)
  
  # max_y=0.9
  gg = ggplot(bin_df, aes(x=bin_mid, y=density, color=dataset, fill=dataset)) +
    geom_bar(stat="identity",position = "identity", alpha=.3)+
    geom_errorbar(aes(ymin=quant_low, ymax=quant_high), width=.2,
                  position=position_dodge(.9)) +
    scale_color_manual(values = c("reference" = rgb(169/255,23/255,23/255, alpha=0.8),
                                  "simulation"=rgb(0/255,124/255,180/255, alpha=0.8))) +
    scale_fill_manual(values = c("reference" = rgb(169/255,23/255,23/255, alpha=0.8),
                                 "simulation"=rgb(0/255,124/255,180/255, alpha=0.8))) +
    # ylim(NA,max_y) +
    labs(title=paste0('infection duration'), x='infection duration (days)') + 
    facet_grid(scales='fixed', facets=c('censor_type', 'age_group'))
  return(gg)
  # ggsave(filename=paste0(sim_dir, '/duration_comparison_by_agebin_xEIR', xEIR_val, '_densThresh', pos_thresh_dens, '.png'), plot=gg, width=8, height=4.5, units='in')
}


plot_infection_duration_dist_by_age_ref_only = function(ref_df, sim_dir, age_bin_lower=c(0, 5, 10, 20, 100), duration_bins=c(seq(0,350,50), 500)){
  
  # get densities for each duration bin
  ref_bin_df = get_duration_bins(ind_data=ref_df, duration_bins=duration_bins, facet_censored=TRUE, facet_age=TRUE)
  ref_bin_df$dataset = 'reference'
  ref_bin_df$censor_type = 'start & finish observed'
  ref_bin_df$censor_type[ref_bin_df$censored] = 'censored'
  bin_df = ref_bin_df
  
  # max_y=0.9
  gg = ggplot(bin_df, aes(x=bin_mid, y=density, color=dataset, fill=dataset)) +
    geom_bar(stat="identity", position = "identity", alpha=.3)+
    scale_color_manual(values = c("reference" = rgb(169/255,23/255,23/255, alpha=0.8),
                                  "simulation"=rgb(0/255,124/255,180/255, alpha=0.8))) +
    scale_fill_manual(values = c("reference" = rgb(169/255,23/255,23/255, alpha=0.8),
                                 "simulation"=rgb(0/255,124/255,180/255, alpha=0.8))) +
    # ylim(NA,max_y) +
    labs(title=paste0('infection duration - reference'), x='infection duration (days)') + 
    facet_grid(scales='fixed', facets=c('censor_type', 'age_group'))
  
  ggsave(filename=paste0(sim_dir, '/duration_by_agebin_reference.png'), plot=gg, width=8, height=4.5, units='in')
}







#########################################################################
# Create all reference-simulation comparison plots: infection duration
#########################################################################

plot_duration_ref_sim_comparison = function(sim_dir, ref_df){
  ref_df$date = as.Date(ref_df$date)
  sim_data = get_sim_survey(sim_dir=sim_dir, ref_df=ref_df)
  
  # create comparison plots
  gg1 = plot_infection_duration_dist(ref_df=ref_df, sim_data=sim_data, pos_thresh_dens=pos_thresh_dens, sim_dir=sim_dir, duration_bins=duration_bins)
  gg2 = plot_infection_duration_dist_by_age(ref_df=ref_df, sim_data=sim_data, pos_thresh_dens=pos_thresh_dens, sim_dir=sim_dir, duration_bins=duration_bins)
  gg3 = create_barplot_frac_comparison(ref_df=ref_df, sim_data=sim_data, pos_thresh_dens=pos_thresh_dens, sim_dir=sim_dir)
  
  return(list(gg1, gg2, gg3))
}










