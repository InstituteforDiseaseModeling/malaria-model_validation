# plot_validation_age_prev_sim_sweep_ref.R

# Note: see related script plot_validation_age_inc_sim_sweep_ref.R for additional plot options

library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)


############################ setup: filepaths and info on simulations run #################################
filepath_ref = "/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/reference_data/Cameron_age_incidence_eir.csv"
filepath_ref2 = "/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/reference_data/age_prev_additional_data.csv"
base_filepath_sim_output = '/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/simulation_output/ewan_age_incidence/updated_eirs'


# sweep of simulations across seasonalities, EIRs, and CMs
filepath_sim = file.path(base_filepath_sim_output, 'sweepSeasonEIRCM', 'inc_prev_data_final.csv')


############################ helper functions #########################################
get_substr = function(site_name_str, index){
  strsplit(site_name_str, "_")[[1]][index]
}
get_mean_from_upper_age = function(cur_age, upper_ages){
  mean_ages = (c(0, upper_ages[1:(length(upper_ages)-1)]) + upper_ages) / 2
  return(mean_ages[which(upper_ages == cur_age)])
}



# get the sweep simulation seasonality and cm that best matches the site simulation for this site
get_matching_sweep_info = function(df, df_sweep, sites_cm_seek, study_site_closest_seasonality){
  df$CM = NA
  df$season = NA
  sweep_cms = unique(df_sweep$CM)
  # iterate through study sites in df
  study_sites = unique(df$Site)
  for (ss in study_sites){
    cur_cm = sum(sites_cm_seek[[paste0('X', ss)]])
    closest_cm = sweep_cms[which.min(abs(sweep_cms-cur_cm))]
    cur_season = study_site_closest_seasonality[[ss]]
    df$CM[df$Site==ss] = closest_cm
    df$season[df$Site==ss] = cur_season
  }
  df$season = factor(df$season, levels = levels(df_sweep$season))
  return(df)
}


format_ref_data_cameron = function(ref_filename){
  # Cameron et al. dataset
  ref_df = read.csv(ref_filename)
  ref_df = ref_df[!is.na(ref_df$PR_LAR),]
  ref_df = ref_df[!is.na(ref_df$REF_GROUP),]
  ref_df$EIR_FOR_SIM = ref_df$EIR_FOR_SIM_V3
  ref_df$EIR_FOR_SIM = as.numeric(ref_df$EIR_FOR_SIM)
  ref_df$SEASON_FOR_SIM = factor(ref_df$SEASON_FOR_SIM, levels=c('flatSeason', 'midSeason', 'highSeason'))
  ref_df$mean_age = (ref_df$PR_LAR + ref_df$PR_UAR) / 2
  
  return(ref_df[,c("REF_GROUP", "ACD_LOCATION", "COUNTRY", 'mean_age', 'SEASON_FOR_SIM', 'START_YEAR', 'PR')])
}

format_ref_data_other = function(ref_filename){
  # non-Cameron et al. dataset
  ref_df = read.csv(ref_filename)
  ref_df = ref_df[!is.na(ref_df$mean_age),]
  ref_df = ref_df[!is.na(ref_df$REF_GROUP),]
  ref_df$SEASON_FOR_SIM = factor(ref_df$SEASON_FOR_SIM, levels=c('flatSeason', 'midSeason', 'highSeason'))

  return(ref_df[,c("REF_GROUP", "ACD_LOCATION", "COUNTRY", 'mean_age', 'SEASON_FOR_SIM', 'START_YEAR', 'PR')])
}

load_all_ref_dfs = function(ref_filename_cameron, ref_filename_other){
  df1 = format_ref_data_cameron(ref_filename_cameron)
  df2 = format_ref_data_other(ref_filename_other)
  ref_df = rbind(df1, df2)
  return(ref_df)
}

###################### plot simulation sweeps ##############################

df_sweep = read.csv(filepath_sim)
df_sweep$season = sapply(df_sweep$Site, get_substr, index=1)
df_sweep$season = factor(df_sweep$season, levels=c('flatSeason', 'midSeason', 'highSeason'))
df_sweep$EIR = as.numeric(sapply(df_sweep$Site, get_substr, index=2))
df_sweep$CM = as.numeric(sapply(df_sweep$Site, get_substr, index=3))
upper_sweep_ages = unique(df_sweep$Age)
df_sweep$mean_age = sapply(df_sweep$Age, get_mean_from_upper_age, upper_ages = upper_sweep_ages)

gg= ggplot(df_sweep, aes(x=mean_age, y=Prevalence)) + 
  geom_line(aes(color=EIR, group=interaction(CM, EIR, season))) +
  scale_colour_gradient(
    low = rgb(0.2,0.2,0),
    high = rgb(0,0.8,0.9),
    trans='log',
    label = function(x) sprintf("%.1f", x)
  ) +
  coord_cartesian(xlim=c(0, 55)) +
  xlab('age') +
  ylab('prevalence') + 
  facet_grid(CM ~ season) +
  theme_bw()
ggsave(file.path(base_filepath_sim_output, '_plots_prevalence', paste0('sweepSeasonEIRCM.png')), gg, width=7, height=4, units='in')



########################## plot reference data against simulation sweep ####################
cm_match=35
match_seasonality = FALSE
matchedString = ''
if (match_seasonality) matchedString='_matchSeason'
ref_df = load_all_ref_dfs(filepath_ref, filepath_ref2)
ss_plot_list = unique(ref_df$REF_GROUP)
gg_list = list()
ii=1
for (ss in ss_plot_list){
  if (ss %in% ref_df$REF_GROUP){
    site_name = gsub('\\?', '_', ref_df$ACD_LOCATION[ref_df$REF_GROUP == as.numeric(ss)][1])
    ref_cur = data.frame('Prevalence'=ref_df$PR[ref_df$REF_GROUP == as.numeric(ss)], 'mean_age'=ref_df$mean_age[which(ref_df$REF_GROUP == as.numeric(ss))])
    ref_cur = distinct(ref_cur)
    if(match_seasonality){
      cur_season = ref_df$SEASON_FOR_SIM[ref_df$REF_GROUP == as.numeric(ss)][1]
    } else cur_season='midSeason'
    
    # check that there are at least three data points
    if(nrow(ref_cur)>2){
      gg= ggplot(df_sweep[intersect(which(df_sweep$CM==cm_match), which(df_sweep$season==cur_season)),], aes(x=mean_age, y=Prevalence)) + # , linetype="season"
        geom_line(aes(color=EIR, group=interaction(EIR))) +
        scale_colour_gradient(
          low = rgb(0.2,0.2,0),
          high = rgb(0,0.8,0.9),
          trans='log',
          label = function(x) sprintf("%.1f", x)
        ) +
        coord_cartesian(xlim=c(0, 55)) +
        xlab('age') +
        ylab('prevalence') + 
        theme_bw()+
        theme(legend.position='none')
      
      gg2 = gg + 
        geom_line(data=ref_cur, aes(x=mean_age, y=Prevalence), color='red', size=1.1, alpha=0.25) + 
        geom_point(data=ref_cur, aes(x=mean_age, y=Prevalence), color='red', size=2, alpha=1) +
        annotate("text",  x=Inf, y = Inf, label = site_name, vjust=1.2, hjust=1.1, color=rgb(0.5,0.5,0.5))
      gg_list[[ii]] = gg2
      ii=ii+1
    }
  }
 }

do.call("grid.arrange", c(gg_list, ncol=ceiling(sqrt(length(gg_list)))))
dev.copy(png, file.path(base_filepath_sim_output, '_plots_prevalence',paste0("sweepSeasonEIRCM_allRefSites", matchedString, ".png")), width=16, height=9, units='in', res=900)
dev.off()




########################## plot all reference data ####################
ref_df = load_all_ref_dfs(filepath_ref, filepath_ref2)
ss_plot_list = unique(ref_df$REF_GROUP)  # 10:max(ref_df$REF_GROUP, na.rm=TRUE)
# get mean of age bin and assign peak prevalence group

# group by peak incidence among any sampled age group
ref_df$prev_group = NA
prev_groups = c(0, .25, .5, 0.75)
prev_group_levels = paste0('max prev: ', c('0-0.25', '0.25-0.5', '0.5-0.75', '>0.75'))
for (ss in ss_plot_list){
  if (ss %in% ref_df$REF_GROUP){
    peak_prev = max(ref_df$PR[which(ref_df$REF_GROUP == as.numeric(ss))], na.rm=TRUE)
    cur_prev_group = prev_group_levels[max(which(prev_groups<peak_prev))]
    ref_df$prev_group[which(ref_df$REF_GROUP == as.numeric(ss))] = cur_prev_group
  }
}
ref_df$prev_group = factor(ref_df$prev_group, levels=prev_group_levels)

gg=ggplot(ref_df, aes(x=mean_age, y=PR, color=as.factor(COUNTRY), group=interaction(REF_GROUP, SEASON_FOR_SIM)), alpha=0.5) +
  geom_line() +
  geom_point(size=1.5, alpha=0.3) +
  xlab('age (midpoint of age bin)') +
  ylab('prevalence (microscopy)') + 
  theme_bw() +
  facet_grid(vars(prev_group), vars(SEASON_FOR_SIM))
ggsave(file.path(base_filepath_sim_output, '_plots_prevalence', paste0('reference_age_prev_allSites.png')), gg, width=9, height=5.5, units='in')














################################################################################################################
# Garki dataset - plot reference from https://github.com/pselvaraj87/Malaria-GarkiDB
################################################################################################################

# for months, Assign each month to a season based on the Garki project timeline. Link: http://garkiproject.nd.edu/garki-timeline.html
    # seasonType = {'1': 'DC', '2': 'DC', '3': 'DC', '4': 'DC', '5': 'DH', '6': 'DH', '7': 'DH', '8': 'W', '9': 'W',
    #                '10': 'W', '11': 'W', '0': 'W'}


library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)
library(data.table)

############################ setup: filepaths and info on simulations run #################################
filepath_ref = "/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/reference_data/Garki_df.csv"
included_villages = c('Rafin Marke', 'Sugungum', 'Matsari')
included_villages = c('Kwaru', 'Baribari', 'Tafin Sale')


ref_df = fread(filepath_ref)
ref_df$par_pos = ref_df$Parasitemia>0
ref_df = ref_df[ref_df$Village %in% included_villages,]

# order factors for seasons
season_levels = paste0(rep(c('DC', 'DH', 'W'),5), rep(2:7, each=3))
ref_df$Seasons = factor(ref_df$Seasons, levels=season_levels)


ref_sum = ref_df %>% group_by(`Age Bins`, Village, Seasons) %>%
  summarise(total_sampled = n(),
            num_pos = sum(par_pos))
ref_sum$prevalence = ref_sum$num_pos / ref_sum$total_sampled

ggplot(ref_sum, aes(x=`Age Bins`, y=prevalence)) +
  geom_point() + 
  geom_line() + 
  facet_grid(Village~Seasons)




age_bin_max = rep(NA, length(unique(ref_df$`Age Bins`)))
for (ii in 1:length(age_bin_max)){
  age_bin_max[ii] = max(ref_df$Age[ref_df$`Age Bins`==ii])
}
age_bin_max = round(age_bin_max, 0)