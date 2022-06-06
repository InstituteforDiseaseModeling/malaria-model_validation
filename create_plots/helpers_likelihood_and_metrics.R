# helpers_loglikelihood_and_metrics.R

#  These functions calculate likelihoods or other comparison metrics evaluating how well simulation and reference data agree.
#    The main functions return a data frame where each row corresponds to a site (or site-month) and columns contain the quantitative measure. 
#    Some functions also return a ggplot corresponding to the quantitative comparison.

# The loglikelihood evaluations are approximate and generally assume that the mean simulated value is the 'true' population value and ask how 
#    likely it was to observe the reference dataset (given the study sample size).

library(broom)
library(tidyverse)



####################################################################################################### 
# loglikelihood functions for each validation relationship
####################################################################################################### 

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# prevalence
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
get_prev_likelihood = function(combined_df, sim_column='simulation'){
  combined_df$prob_pos_sim = combined_df[,sim_column]
  # only include reference sites where the sample sizes were reported
  combined_df = combined_df[!is.na(combined_df$total_sampled),]
  sites = unique(combined_df$site_month)
  
  # likelihood approximated for each age group as probability of obtaining observed num_pos from total_sampled if simulation is true prevalence
  loglik_by_site = rep(NA, length(sites))
  for(ss in 1:length(sites)){
    cur_df = combined_df[combined_df$site_month==sites[ss],]
    # get product of probabilities for each age group. If any age groups don't match, entire site is NA.
    if(!any(is.na(cur_df$prob_pos_sim))){
      loglik_total=0
      for(rr in 1:nrow(cur_df)){  # iterate through age groups (each row corresponds to a different age group)
        loglik_total = loglik_total + log(dbinom(x=cur_df$num_pos[rr], size=cur_df$total_sampled[rr], prob=cur_df$prob_pos_sim[rr]))
      }
      loglik_by_site[ss] = loglik_total
    }
  }
  loglik_df = data.frame('site_month' = sites, 'loglikelihood'=loglik_by_site)
  
  return(loglik_df)
}





# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# parasite density
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
get_dens_likelihood = function(combined_df, sim_column='simulation'){

  # remove rows where there is no reference data (sometimes there are different density bins in different sites, so rows with NA are created for the 'missing' bins - but check that there aren't values in the simulation either)
  ref_rows_na = which(is.na(combined_df$ref_bin_count))
  if(length(ref_rows_na)>0){
    sim_rows_0 = which(combined_df[[sim_column]] < 0.0001)
    if(all(ref_rows_na %in% sim_rows_0)){
      combined_df = combined_df[-ref_rows_na,]
    } else{
      warning('There may be a mismatch in the age bins from the reference data and simulation data for at least one site. No rows were removed.')
    }
  }
  
  # remove rows without simulation output
  combined_df = combined_df[!is.na(combined_df[[sim_column]]),]
  
  # multinomial draw: likelihood of getting the observed distribution of parasite densities assuming the simulations show the true population-level frequencies
  # iterate through groups of month-age-site
  loglik_df = data.frame('site_month'= c(), 'loglikelihood' = c())
  site_months = unique(combined_df$site_month)
  for (ss in site_months){
    loglikelihood = 0
    cur_agebins = unique(combined_df$agebin[combined_df$site_month == ss])
    for (aa in cur_agebins){
      cur_df = combined_df[combined_df$site_month == ss & combined_df$agebin == aa,]
      # check that the sum of counts matches the sum column
      if(sum(cur_df$ref_bin_count) == cur_df$ref_total[1] & length(unique(cur_df$ref_total))==1){
        loglikelihood = loglikelihood + log(dmultinom(x=cur_df$ref_bin_count, prob=cur_df[[sim_column]]))
      } else{
        warning(paste0('Either the sum of individuals across bins in the reference dataset does not match the reported total number of individuals included or different density bins are used in the reference and simulation. This site-age is being skipped: ', ss, ' - ', aa))
      }
    }
    loglik_df = rbind(loglik_df, data.frame('site_month'= ss, 'loglikelihood' = loglikelihood))
  }
  return(loglik_df)
}







####################################################################################################### 
# other quantitative comparison metrics
####################################################################################################### 

# mean relative difference between reference and matched simulation values across ages
#   for a given age, calculate (reference - simulation)/reference. Report the mean across all ages
calc_mean_rel_diff = function(combined_df){
  if('site_month' %in% colnames(combined_df)) combined_df$Site = combined_df$site_month
  combined_df$rel_diff = abs((combined_df$reference - combined_df$simulation) / combined_df$reference)
  combined_df$abs_diff = abs(combined_df$reference - combined_df$simulation)
  mean_diff_df = combined_df %>% group_by(Site) %>%
    summarise(mean_rel_diff = mean(rel_diff),
              mean_abs_diff = mean(abs_diff))
  return(mean_diff_df)
}


calc_mean_rel_slope_diff = function(combined_df){
  if('site_month' %in% colnames(combined_df)) combined_df$Site = combined_df$site_month
  combined_df$rel_diff = abs((combined_df$ref_slope_to_next - combined_df$sim_slope_to_next) / combined_df$ref_slope_to_next)
  combined_df$abs_diff = abs(combined_df$ref_slope_to_next - combined_df$sim_slope_to_next)
  mean_slope_diff_df = combined_df %>% group_by(Site) %>%
    summarise(mean_rel_slope_diff = mean(rel_diff, na.rm=TRUE),
              mean_abs_slpe_diff = mean(abs_diff, na.rm=TRUE))
  return(mean_slope_diff_df)
}


# correlation between reference and matched simulation data points
corr_ref_sim_points = function(combined_df){
  metric = combined_df$metric[1]
  if('site_month' %in% colnames(combined_df)) combined_df$Site = combined_df$site_month
  min_value = min(c(combined_df$reference, combined_df$simulation), na.rm=TRUE)
  max_value = max(c(combined_df$reference, combined_df$simulation), na.rm=TRUE)
  gg = ggplot(combined_df, aes(x=reference, y=simulation, color=Site, fill=Site))+
    geom_point(size=2) + 
    ylab(paste0('simulation ', metric)) + 
    xlab(paste0('reference ', metric)) + 
    ggtitle('Ref versus sim values') +
    coord_fixed(ratio=1, xlim=c(min_value, max_value), ylim=c(min_value, max_value))+
    geom_abline(slope=1,intercept=0, color='grey', alpha=0.5) + 
    geom_smooth(method = "lm", fill = NA, se=FALSE, alpha=0.5, size=0.5) +
    theme_classic() +
    theme(plot.title = element_text(size=12))
  
  
  lm_fit = combined_df %>% nest(data = -Site) %>% mutate(model = map(data, ~lm(simulation ~ reference, data = .)), tidied = map(model, tidy)) %>% unnest(tidied)
  lm_info = combined_df %>% nest(data = -Site) %>% mutate(model = map(data, ~lm(simulation ~ reference, data = .)), tidied = map(model, glance)) %>% unnest(tidied)
  
  lm_summary = merge(lm_fit[c('Site', 'term', 'estimate')], lm_info[c('Site', 'r.squared', 'p.value', 'nobs')], by='Site', all=TRUE)
  lm_summary = lm_summary[lm_summary$term!='(Intercept)',]
  colnames(lm_summary)[colnames(lm_summary)=='estimate'] = 'slope'
  return(list(gg, lm_summary))
}


# correlation between derivatives moving between ages
corr_ref_deriv_sim_points = function(combined_df){
  metric = combined_df$metric[1]
  # calculate the slope when moving between age groups
  combined_df$sim_slope_to_next = NA
  combined_df$ref_slope_to_next = NA
  sites = unique(combined_df$Site)
  for(ss in sites){
    cur_df = combined_df[combined_df$Site == ss,]
    cur_ages = sort(unique(cur_df$mean_age))
    for(aa in 1:(length(cur_ages)-1)){
      combined_df_row = intersect(which(combined_df$Site == ss), which(combined_df$mean_age == cur_ages[aa]))
      # get the slope between the value for this age and the next-largest age group
      sim_val_cur = cur_df$simulation[cur_df$mean_age == cur_ages[aa]]
      sim_val_next = cur_df$simulation[cur_df$mean_age == cur_ages[aa+1]]
      sim_slope = (sim_val_next - sim_val_cur) / (cur_ages[aa+1] - cur_ages[aa])
      combined_df$sim_slope_to_next[combined_df_row] = sim_slope
      
      ref_val_cur = cur_df$reference[cur_df$mean_age == cur_ages[aa]]
      ref_val_next = cur_df$reference[cur_df$mean_age == cur_ages[aa+1]]
      ref_slope = (ref_val_next - ref_val_cur) / (cur_ages[aa+1] - cur_ages[aa])
      combined_df$ref_slope_to_next[combined_df_row] = ref_slope
    }
  }
  min_value = min(c(combined_df$ref_slope_to_next, combined_df$sim_slope_to_next), na.rm=TRUE)
  max_value = max(c(combined_df$ref_slope_to_next, combined_df$sim_slope_to_next), na.rm=TRUE)
  gg = ggplot(combined_df, aes(x=ref_slope_to_next, y=sim_slope_to_next, color=Site, fill=Site))+
    geom_point(size=2) + 
    ylab(paste0('simulation slopes for age-', metric)) + 
    xlab(paste0('reference slopes for age-', metric)) + 
    ggtitle('Ref versus sim slopes') +
    coord_fixed(ratio=1, xlim=c(min_value, max_value), ylim=c(min_value, max_value))+
    geom_abline(slope=1,intercept=0, color='grey', alpha=0.5) + 
    # geom_smooth(method = "lm", fill = NA, se=FALSE, alpha=0.5, size=0.5) +
    theme_classic() +
    theme(plot.title = element_text(size=12))
  
  
  lm_fit = combined_df %>% nest(data = -Site) %>% mutate(model = map(data, ~lm(sim_slope_to_next ~ ref_slope_to_next, data = .)), tidied = map(model, tidy)) %>% unnest(tidied)
  lm_info = combined_df %>% nest(data = -Site) %>% mutate(model = map(data, ~lm(sim_slope_to_next ~ ref_slope_to_next, data = .)), tidied = map(model, glance)) %>% unnest(tidied)
  
  lm_summary = merge(lm_fit[c('Site', 'term', 'estimate')], lm_info[c('Site', 'r.squared', 'p.value', 'nobs')], by='Site', all=TRUE)
  lm_summary = lm_summary[lm_summary$term!='(Intercept)',]
  colnames(lm_summary)[colnames(lm_summary)=='estimate'] = 'slope'
  return(list(gg, lm_summary, combined_df))
}





