# helper_functions_age_inc.R

# plot the match between the reference datasets and matched simulations

library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)




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
    scale_color_manual(values = c("reference" = "red",
                                  "simulation"="blue")) +
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
    scale_color_manual(values = c("reference" = "red",
                                  "simulation"="blue")) +
    xlab('age (midpoint of age bin)') +
    ylab('prevalence') +
    facet_wrap('site_month', ncol=4) +
    theme_bw()  
  
  return(gg)
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



