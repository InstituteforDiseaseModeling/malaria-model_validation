#####################################################################################
# compare_ref_sim_infection_duration_Navrongo
# contact: mambrose
# Jan 2022
#
# Goal: simulate transmission in EMOD and sample to match reference dataset
#   Compare:
#   - probability a person goes from negative to positive between sample dates
#   - probability a person goes from positive to negative between sample dates
#   - fraction of samples positive
#   - among all individuals starting positive, how long until negative?
#   - among all individuals turning positive during study, how long until negative?
#####################################################################################


library(stringr)
library(ggplot2)
library(dplyr)
library(data.table)
library(RcmdrMisc)

#########################################################################
# setup
#########################################################################

dropbox_filepath = "C:/Users/moniqueam/Dropbox (IDM)/Malaria Team Folder"
sim_dir = 'C:/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/simulation_output/duration_of_infection/model_validation_duration_Navrongo_30y_xEIR_v2'
xEIR_vals = c(0.1, 0.5, 0.8, 1, 1.2)
seeds = 0:9  # c(0)
# set positive threshold density for sampled parasites in simulation output (to match PCR threshold in reference)
pos_thresh_dens = 0.5  # (39=min(ref_data$DENSITY[ref_data$DENSITY>0], na.rm=TRUE) - 1)
# specify binning for duration of infection
duration_bins=c(seq(0,350,50), 500)


########################################################################
# functions
########################################################################

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# get subset of simulation dataset to match reference dataset for dates and ages of individuals of sampled individuals
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
get_sim_survey = function(sim_dir, xEIR_val, ref_data, seeds=NA){
  sampled_sim_filename = paste0('sampled_sim_xEIR', xEIR_val, '.csv')
  if(file.exists(paste0(sim_dir, '/', sampled_sim_filename))){
    sim_subset_full = read.csv(paste0(sim_dir, '/', sampled_sim_filename))
  } else{
    sim_full = fread(paste0(sim_dir, '/patient_reports.csv'))
    sim_full$Site[sim_full$Site=='Navrongo_x01'] = 'Navrongo_x010'  # TODO: remove this for new experiments - currently included to fix mis-labeled version
    sim_full$x_EIR = as.numeric(sub(".*_x", "", sim_full$Site))/10
    sim_full = sim_full[which(sim_full$x_EIR == xEIR_val),]  # subset to desired xEIR
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
        ref_data_cur = ref_data[ref_data$SID == indIDs[ii],]
        ref_data_cur = ref_data_cur[order(ref_data_cur$date),]
        # find a matching individual
        age_cur = ref_data_cur$age[1]
        day_cur = ref_data_cur$date[1]
        
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
        sim_subset_cur = sim[intersect(which(sim$id == id_sim_cur), which(sim$date %in% ref_data_cur$date)),]
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
  } else{
    bin_w_seeds$density = bin_w_seeds[,which(grepl('seed_', colnames(bin_w_seeds)))]
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
create_barplot_frac_comparison = function(ref_data, sim_subset, xEIR_val, seed, pos_thresh_dens, sim_dir){
  frac_swap_ref = get_frac_state_swaps(data=ref_data)
  frac_swap_sim = get_frac_state_swaps(data=sim_subset)
  
  ref_data$pos = ref_data$DENSITY > pos_thresh_dens
  frac_samples_pos_ref = sum(ref_data$pos, na.rm=TRUE) / sum((ref_data$DENSITY > -1), na.rm=TRUE)
  
  sim_subset$pos = sim_subset$DENSITY > pos_thresh_dens
  frac_samples_pos_sim = sum(sim_subset$pos, na.rm=TRUE) / sum((sim_subset$DENSITY > -1), na.rm=TRUE)
  
  png(filename = paste0(sim_dir, '/fraction_comparison_xEIR', xEIR_val, '_seed', seed, '_densThresh', pos_thresh_dens, '.png'), width=7, height=5, units='in', res=900)
  x=barplot(matrix(c( frac_swap_ref[2], frac_swap_ref[1], frac_samples_pos_ref, frac_swap_sim[2], frac_swap_sim[1], frac_samples_pos_sim), nrow=2, byrow=TRUE), 
            # beside=T, col=c(rgb(0.2,0.2,0.2), rgb(0.9,0.9,0.9)), ylim=c(0,0.6))
            beside=T, col=c(rgb(0.9,0.4,0.45,0.5), rgb(0.3,0.9,0.9,0.5)), border=c(rgb(0.9,0.4,0.45), rgb(0.3,0.9,0.9)), ylim=c(0,0.6))
  text(cex=1, x=x[1,]+0.5, y=-0.1, xpd=TRUE, adj=0.5, srt=0, c('fraction of negatives\n turn positive next time', 'fraction of positives\n turn negative next time', 'fraction of samples\n positive'))
  legend('topleft', c('reference', 'simulation'), fill=c(rgb(0.9,0.4,0.45,0.5), rgb(0.3,0.9,0.9,0.5)), bty='n')
  dev.off()
}

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# create plots comparing the distributions of infection lengths
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# plot_infection_duration_dist = function(ref_data, sim_subset, xEIR_val, seed, pos_thresh_dens, sim_dir){
#   days_positive_ref = get_time_pos(data=ref_data)
#   days_positive_ref$dataset = 'reference'
#   days_positive_sim = get_time_pos(data=sim_subset[sim_subset$seed==seed,])
#   days_positive_sim$dataset = 'simulation'
#   
#   # combine multiple simulation seeds to show mean, min, and max values
#   
#   # combine reference and simulation datasets
#   days_positive = rbind(days_positive_ref, days_positive_sim)
#   days_positive$censor_type = 'start & finish observed'
#   days_positive$censor_type[days_positive$censored] = 'censored'
#   
#   max_x = 450
#   max_y=0.016
#   gg = ggplot(days_positive, aes(x=days_positive, color=dataset, fill=dataset)) +
#     xlim(NA, max_x) +
#     ylim(NA,max_y) +
#     geom_histogram(aes(y=..density..), alpha=0.25, position="identity", binwidth=50)+
#     labs(title=paste0('infection duration - sim xEIR=', xEIR_val), x='infection duration (days)') + 
#     facet_wrap(scales='fixed', facets='censor_type')
#   
#   ggsave(filename=paste0(sim_dir, '/duration_comparison_xEIR', xEIR_val, '_densThresh', pos_thresh_dens, '.png'), plot=gg, width=6, height=3, units='in')
# }
# 
# 
# plot_infection_duration_dist_ref_only = function(ref_data, sim_dir){
#   days_positive_ref = get_time_pos(data=ref_data)
#   days_positive_ref$dataset = 'reference'
#   days_positive = days_positive_ref
#   days_positive$censor_type = 'start & finish observed'
#   days_positive$censor_type[days_positive$censored] = 'censored'
#   
#   max_x = 450
#   max_y=0.016
#   gg = ggplot(days_positive, aes(x=days_positive, color=dataset, fill=dataset)) +
#     xlim(NA, max_x) +
#     ylim(NA,max_y) +
#     geom_histogram(aes(y=..density..), alpha=0.25, position="identity", binwidth=50)+
#     labs(title=paste0('infection duration - reference dataset'), x='infection duration (days)') + 
#     facet_wrap(scales='fixed', facets='censor_type')
#   
#   ggsave(filename=paste0(sim_dir, '/duration_allAges_reference.png'), plot=gg, width=6, height=3, units='in')
# }




plot_infection_duration_dist = function(ref_data, sim_data, xEIR_val, pos_thresh_dens, sim_dir, duration_bins=c(seq(0,350,50), 500)){

  # get densities for each duration bin
  ref_bin_df = get_duration_bins(ind_data=ref_data, duration_bins=duration_bins, facet_censored=TRUE, facet_age=FALSE)
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
    # ylim(NA,max_y) +
    labs(title=paste0('infection duration - sim xEIR=', xEIR_val), x='infection duration (days)') + 
    facet_wrap(scales='fixed', facets='censor_type')
  ggsave(filename=paste0(sim_dir, '/duration_comparison_xEIR', xEIR_val, '_densThresh', pos_thresh_dens, '.png'), plot=gg, width=6, height=3, units='in')
}



plot_infection_duration_dist_ref_only = function(ref_data, sim_dir, duration_bins=c(seq(0,350,50), 500)){
  
  # get densities for each duration bin
  ref_bin_df = get_duration_bins(ind_data=ref_data, duration_bins=duration_bins, facet_censored=TRUE, facet_age=FALSE)
  ref_bin_df$dataset = 'reference'
  ref_bin_df$censor_type = 'start & finish observed'
  ref_bin_df$censor_type[ref_bin_df$censored] = 'censored'
  bin_df = ref_bin_df
  
  # max_y=0.8
  gg = ggplot(bin_df, aes(x=bin_mid, y=density, color=dataset, fill=dataset)) +
    geom_bar(stat="identity",position = "identity", alpha=.3)+
    # ylim(NA,max_y) +
    labs(title=paste0('infection duration - reference'), x='infection duration (days)') + 
    facet_wrap(scales='fixed', facets='censor_type')
  
  ggsave(filename=paste0(sim_dir, '/duration_allAges_reference.png'), plot=gg, width=6, height=3, units='in')
}



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# create plots comparing the distributions of infection lengths, faceted by age group
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# plot_infection_duration_dist_by_age = function(ref_data, sim_subset, xEIR_val, seed, pos_thresh_dens, sim_dir, age_bin_lower=c(0, 5, 10, 20, 100)){
#   days_positive_ref = get_time_pos(data=ref_data)
#   days_positive_ref$dataset = 'reference'
#   days_positive_sim = get_time_pos(data=sim_subset)
#   days_positive_sim$dataset = 'simulation'
#   days_positive = rbind(days_positive_ref, days_positive_sim)
#   days_positive$censor_type = 'start & finish observed'
#   days_positive$censor_type[days_positive$censored] = 'censored'
#   
#   max_x = 450
#   max_y=0.02
#   
#   # plot subset by age
#   age_bin_labels = c(paste0(age_bin_lower[1:(length(age_bin_lower)-1)], '-', age_bin_lower[2:length(age_bin_lower)]), paste0('>', age_bin_lower[length(age_bin_lower)]))
#   get_age_group = function(cur_age, age_bin_lower, age_bin_labels){
#     return(age_bin_labels[max(which(cur_age>age_bin_lower))])
#   }
#   days_positive$age_group = sapply(days_positive$age, get_age_group, age_bin_lower=age_bin_lower, age_bin_labels=age_bin_labels)
#   days_positive$age_group = factor(days_positive$age_group, levels=age_bin_labels)
#   
#   gg=ggplot(days_positive, aes(x=days_positive, color=dataset, fill=dataset)) +
#     xlim(NA, max_x) +
#     ylim(NA,max_y) +
#     geom_histogram(aes(y=..density..), alpha=0.25, position="identity", binwidth=50)+
#     labs(title=paste0('infection duration - sim xEIR=', xEIR_val), x='infection duration (days)') + 
#     facet_grid(scales='fixed', facets=c('censor_type', 'age_group'))
#   
#   ggsave(filename=paste0(sim_dir, '/duration_comparison_by_agebin_xEIR', xEIR_val, '_densThresh', pos_thresh_dens, '.png'), plot=gg, width=8, height=4.5, units='in')
# }
# 
# 
# plot_infection_duration_dist_by_age_ref_only = function(ref_data, sim_dir, age_bin_lower=c(0, 5, 10, 20, 100)){
#   days_positive_ref = get_time_pos(data=ref_data)
#   days_positive_ref$dataset = 'reference'
#   days_positive = days_positive_ref
#   days_positive$censor_type = 'start & finish observed'
#   days_positive$censor_type[days_positive$censored] = 'censored'
#   
#   max_x = 450
#   max_y=0.02
#   
#   # plot subset by age
#   age_bin_labels = c(paste0(age_bin_lower[1:(length(age_bin_lower)-1)], '-', age_bin_lower[2:length(age_bin_lower)]), paste0('>', age_bin_lower[length(age_bin_lower)]))
#   get_age_group = function(cur_age, age_bin_lower, age_bin_labels){
#     return(age_bin_labels[max(which(cur_age>age_bin_lower))])
#   }
#   days_positive$age_group = sapply(days_positive$age, get_age_group, age_bin_lower=age_bin_lower, age_bin_labels=age_bin_labels)
#   days_positive$age_group = factor(days_positive$age_group, levels=age_bin_labels)
#   
#   gg=ggplot(days_positive, aes(x=days_positive, color=dataset, fill=dataset)) +
#     xlim(NA, max_x) +
#     ylim(NA,max_y) +
#     geom_histogram(aes(y=..density..), alpha=0.25, position="identity", binwidth=50)+
#     labs(title=paste0('infection duration - reference dataset'), x='infection duration (days)') + 
#     facet_grid(scales='fixed', facets=c('censor_type', 'age_group'))
#   ggsave(filename=paste0(sim_dir, '/duration_by_agebin_reference.png'), plot=gg, width=8, height=4.5, units='in')
# }



plot_infection_duration_dist_by_age = function(ref_data, sim_data, xEIR_val, pos_thresh_dens, sim_dir, age_bin_lower=c(0, 5, 10, 20, 100), duration_bins=c(seq(0,350,50), 500)){
  
  # get densities for each duration bin
  ref_bin_df = get_duration_bins(ind_data=ref_data, duration_bins=duration_bins, facet_censored=TRUE, facet_age=TRUE)
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
    # ylim(NA,max_y) +
    labs(title=paste0('infection duration - sim xEIR=', xEIR_val), x='infection duration (days)') + 
    facet_grid(scales='fixed', facets=c('censor_type', 'age_group'))
  ggsave(filename=paste0(sim_dir, '/duration_comparison_by_agebin_xEIR', xEIR_val, '_densThresh', pos_thresh_dens, '.png'), plot=gg, width=8, height=4.5, units='in')
}


plot_infection_duration_dist_by_age_ref_only = function(ref_data, sim_dir, age_bin_lower=c(0, 5, 10, 20, 100), duration_bins=c(seq(0,350,50), 500)){
  
  # get densities for each duration bin
  ref_bin_df = get_duration_bins(ind_data=ref_data, duration_bins=duration_bins, facet_censored=TRUE, facet_age=TRUE)
  ref_bin_df$dataset = 'reference'
  ref_bin_df$censor_type = 'start & finish observed'
  ref_bin_df$censor_type[ref_bin_df$censored] = 'censored'
  bin_df = ref_bin_df
  
  # max_y=0.9
  gg = ggplot(bin_df, aes(x=bin_mid, y=density, color=dataset, fill=dataset)) +
    geom_bar(stat="identity", position = "identity", alpha=.3)+
    # ylim(NA,max_y) +
    labs(title=paste0('infection duration - reference'), x='infection duration (days)') + 
    facet_grid(scales='fixed', facets=c('censor_type', 'age_group'))
  
  ggsave(filename=paste0(sim_dir, '/duration_by_agebin_reference.png'), plot=gg, width=8, height=4.5, units='in')
}





#########################################################################
# Load and format reference dataset: Navrongo, Ghana
#########################################################################
ref_data = read.csv(paste0(dropbox_filepath, "/data/Ghana/Navrongo/slide.csv"))
ref_data$date = as.Date(ref_data$DATE)
# remove 1990s date
ref_data = ref_data[ref_data$date > as.Date('1999-01-01'),]

# get birthdays from other dataset and use to calculate age at time of sample
ref_data_birthdays = read.csv(paste0(dropbox_filepath, "/data/Ghana/Navrongo/person.csv"))
ref_data_birthdays = ref_data_birthdays[,c('DOB', 'SID')]
ref_data = merge(ref_data, ref_data_birthdays, by='SID')
ref_data$DOB = as.Date(ref_data$DOB)
ref_data$age = as.numeric((ref_data$date - ref_data$DOB)/365)
indIDs = unique(ref_data$SID)

# get first year of sampling in reference dataset. the simulation will be referenced from the first day of that year
first_ref_date = as.Date(paste0(year(min(ref_data$date)), '-01-01'))

plot_infection_duration_dist_ref_only(ref_data=ref_data, sim_dir=sim_dir, duration_bins=duration_bins)
plot_infection_duration_dist_by_age_ref_only(ref_data=ref_data, sim_dir=sim_dir, duration_bins=duration_bins)
# # plot longitudinal parasite densities
# ggplot(ref_data, aes(x=date, y=DENSITY, color=SID))+
#   geom_line()+
#   scale_y_continuous(trans='log10') +
#   theme(legend.position="none")



#########################################################################
# Simulation output: iterate through simulation scenarios, sample to match reference dataset, create comparison plots
#########################################################################
xEIR_val = 0.1
for(xEIR_val in xEIR_vals){
  sim_data = get_sim_survey(sim_dir=sim_dir, xEIR_val=xEIR_val, ref_data=ref_data)
  # create comparison plots
  plot_infection_duration_dist(ref_data=ref_data, sim_data=sim_data, xEIR_val=xEIR_val, pos_thresh_dens=pos_thresh_dens, sim_dir=sim_dir, duration_bins=duration_bins)
  plot_infection_duration_dist_by_age(ref_data=ref_data, sim_data=sim_data, xEIR_val=xEIR_val, pos_thresh_dens=pos_thresh_dens, sim_dir=sim_dir, duration_bins=duration_bins)
  create_barplot_frac_comparison(ref_data=ref_data, sim_subset=sim_subset, xEIR_val=xEIR_val, seed=seed, pos_thresh_dens=pos_thresh_dens, sim_dir=sim_dir)
      
}





