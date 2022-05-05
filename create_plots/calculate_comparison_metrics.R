# calculate_comparison_metrics.R

# given simulation and reference datasets, calculate some metrics that quantify 
#   whether they are showing similar/different patterns


library(broom)
library(tidyverse)


# prepare dataframe with simulation and reference data formatted together
prepare_inc_df = function(sim_df, ref_df){
  # scale down incidence in simulation according to probability of detecting a case
  sim_df$Incidence = sim_df$Incidence * sim_df$p_detect_case
  
  # set up reference and simulation dataset columns
  ref_df$Incidence = ref_df$INC / 1000
  ref_df$mean_age = (ref_df$INC_LAR + ref_df$INC_UAR)/2
  ref_df = data.frame('reference_value'=ref_df$Incidence, 'mean_age'=ref_df$mean_age, 'Site'=ref_df$Site, 'ref_pop_size'=ref_df$POP, 'ref_year'=ref_df$START_YEAR)
  sim_df = data.frame('simulation_value'=sim_df$Incidence, 'mean_age'=sim_df$mean_age, 'Site'=sim_df$Site)

  sites = intersect(unique(sim_df$Site), unique(ref_df$Site))
  # check that ages match between reference and simulation. if there is a small difference (<1 year, update simulation)
  for(ss in sites){
    ages_ref = sort(unique(ref_df$mean_age[ref_df$Site == ss]))
    ages_sim = sort(unique(sim_df$mean_age[sim_df$Site == ss]))
    missing_ages_ref = ages_ref[which(!(ages_ref %in% ages_sim))]
    missing_ages_sim = ages_sim[which(!(ages_sim %in% ages_ref))]
    
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
  
  combined_df = merge(ref_df, sim_df, all=TRUE)
  combined_df$metric = 'incidence'
  return(combined_df)
}



# mean relative difference between reference and matched simulation values across ages
#   for a given age, calculate (reference - simulation)/reference. Report the mean across all ages
calc_mean_rel_diff = function(combined_df){
  combined_df$rel_diff = (combined_df$reference_value - combined_df$simulation_value) / combined_df$reference_value
  combined_df$abs_diff = abs(combined_df$reference_value - combined_df$simulation_value)
  mean_diff_df = combined_df %>% group_by(Site) %>%
    summarise(mean_rel_diff = mean(rel_diff),
              mean_abs_diff = mean(abs_diff))
  return(mean_diff_df)
}


# correlation between reference and matched simulation data points
corr_ref_sim_points = function(combined_df){
  metric = combined_df$metric[1]
  gg = ggplot(combined_df, aes(x=reference_value, y=simulation_value, color=Site, fill=Site))+
    geom_point(size=2) + 
    ylab(paste0('simulation ', metric)) + 
    xlab(paste0('reference ', metric)) + 
    geom_abline(slope=1,intercept=0, color='grey', alpha=0.5) + 
    geom_smooth(method = "lm", fill = NA, se=FALSE, alpha=0.5, size=0.5) +
    theme_classic()
  

  lm_fit = combined_df %>% nest(data = -Site) %>% mutate(model = map(data, ~lm(simulation_value ~ reference_value, data = .)), tidied = map(model, tidy)) %>% unnest(tidied)
  lm_info = combined_df %>% nest(data = -Site) %>% mutate(model = map(data, ~lm(simulation_value ~ reference_value, data = .)), tidied = map(model, glance)) %>% unnest(tidied)
  
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
      sim_val_cur = cur_df$simulation_value[cur_df$mean_age == cur_ages[aa]]
      sim_val_next = cur_df$simulation_value[cur_df$mean_age == cur_ages[aa+1]]
      sim_slope = (sim_val_next - sim_val_cur) / (cur_ages[aa+1] - cur_ages[aa])
      combined_df$sim_slope_to_next[combined_df_row] = sim_slope
      
      ref_val_cur = cur_df$reference_value[cur_df$mean_age == cur_ages[aa]]
      ref_val_next = cur_df$reference_value[cur_df$mean_age == cur_ages[aa+1]]
      ref_slope = (ref_val_next - ref_val_cur) / (cur_ages[aa+1] - cur_ages[aa])
      combined_df$ref_slope_to_next[combined_df_row] = ref_slope
    }
  }
  
  gg = ggplot(combined_df, aes(x=ref_slope_to_next, y=sim_slope_to_next, color=Site, fill=Site))+
    geom_point(size=2) + 
    ylab(paste0('simulation slopes for age-', metric)) + 
    xlab(paste0('reference slopes for age-', metric)) + 
    geom_abline(slope=1,intercept=0, color='grey', alpha=0.5) + 
    geom_smooth(method = "lm", fill = NA, se=FALSE, alpha=0.5, size=0.5) +
    theme_classic()
  

  lm_fit = combined_df %>% nest(data = -Site) %>% mutate(model = map(data, ~lm(sim_slope_to_next ~ ref_slope_to_next, data = .)), tidied = map(model, tidy)) %>% unnest(tidied)
  lm_info = combined_df %>% nest(data = -Site) %>% mutate(model = map(data, ~lm(sim_slope_to_next ~ ref_slope_to_next, data = .)), tidied = map(model, glance)) %>% unnest(tidied)
  
  lm_summary = merge(lm_fit[c('Site', 'term', 'estimate')], lm_info[c('Site', 'r.squared', 'p.value', 'nobs')], by='Site', all=TRUE)
  lm_summary = lm_summary[lm_summary$term!='(Intercept)',]
  colnames(lm_summary)[colnames(lm_summary)=='estimate'] = 'slope'
  return(list(gg, lm_summary))
}



