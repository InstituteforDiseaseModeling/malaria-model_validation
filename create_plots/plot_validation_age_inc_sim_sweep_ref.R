# plot_validation_age_inc_sim_sweep_ref.R

# first, plot the simulation relationships between age and incidence for a sweep of EIRs, seasonality patterns, and CM rates
# second, plot the match between the reference datasets, matched simulations, and sweep of EIRs

library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)


############################ setup: filepaths and info on simulations run #################################
filepath_ref = "/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/reference_data/Cameron_age_incidence_eir.csv"
base_filepath_sim_output = '/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/simulation_output/ewan_age_incidence/updated_eirs'
# base_filepath_sim_output = '/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/simulation_output/ewan_age_incidence/original_eirs'
duplicate_sites = c(6, 7)  # remove these sites to avoid duplicating reference data

# sweep of simulations across seasonalities, EIRs, and CMs
filepath_sim = file.path(base_filepath_sim_output, 'sweepSeasonEIRCM', 'summary_data_final.csv')
eirs = c(1, 3, 5, 10, 20, 30, 50, 100, 200, 400)


# simulation on sites matched with reference
sites = c('1', '2', '3', '4', '5', '6', '7', '8', '9')

# monthly EIRs used for each site in simulations
study_site_monthly_EIRs = data.frame(
  'month' = 1:12,
  '1'= 55*c(0.07, 0.08, 0.04, 0.02, 0.03, 0.04, 0.22, 0.13, 0.18, 0.09, 0.07, 0.05),   # Asembo Bay - EIR =  56
  '2'= c(3.6, 1.4, 0.0, 2.8, 0.4, 10.3, 10.6, 3.1, 0.4, 0.0, 0.7, 0.8),               # Chonyi - EIR =  34
  '3'= 10/5.6*c(0.6, 0.25, 0.0, 0.45, 0.05, 1.7, 1.75, 0.5, 0.05, 0.0, 0.1, 0.15),     # Ngerenya - EIR = 10
  '4'= 200/159.8*c(10.4, 13, 6, 2.6, 4, 6, 35, 21, 28, 15, 10.4, 8.4),                 # Dielmo - EIR = 200
  '5'= 155/159.8*c(10.4, 13, 6, 2.6, 4, 6, 35, 21, 28, 15, 10.4, 8.4),                 # Dielmo - EIR = 155
  '6'= 155/159.8*c(10.4, 13, 6, 2.6, 4, 6, 35, 21, 28, 15, 10.4, 8.4),                 # Dielmo - EIR = 155
  '7'= 101.2/159.8*c(10.4, 13, 6, 2.6, 4, 6, 35, 21, 28, 15, 10.4, 8.4),               # Dielmo - EIR = 101.2
  '8'= c(0.39, 0.19, 0.77, 0, 0, 0, 6.4, 2.2, 4.7, 3.9, 0.87, 0.58),                  # Ndiop - EIR = 20
  '9'= c(0.39, 0.19, 0.77, 0, 0, 0, 6.4, 2.2, 4.7, 3.9, 0.87, 0.58)                   # Ndiop - EIR = 20
)
sim_site_annual_eir = colSums(study_site_monthly_EIRs)[-1]
sites_cm_seek = data.frame(
  '1'= 0.35,
  '2'= 0.35,
  '3'= 0.35,
  '4'= 0.5,
  '5'= 0.5,
  '6'= 0.5,
  '7'= 0.5,
  '8'= 0.5,
  '9'= 0.5
)
# guess at which seasonality from sweep best matches each site
study_site_closest_seasonality = data.frame(
  '1'='midSeason',
  '2'='highSeason',
  '3'='highSeason',
  '4'='midSeason',
  '5'='midSeason',
  '6'='midSeason',
  '7'='midSeason',
  '8'='highSeason',
  '9'='highSeason'
)




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
    cur_cm = 100*(sites_cm_seek[[paste0('X', ss)]])
    closest_cm = sweep_cms[which.min(abs(sweep_cms-cur_cm))]
    cur_season = study_site_closest_seasonality[[ss]]
    df$CM[df$Site==ss] = closest_cm
    df$season[df$Site==ss] = cur_season
  }
  df$season = factor(df$season, levels = levels(df_sweep$season))
  return(df)
}


###################### plot simulation sweeps ##############################

df_sweep = read.csv(filepath_sim)
df_sweep$season = sapply(df_sweep$Site, get_substr, index=1)
df_sweep$season = factor(df_sweep$season, levels=c('flatSeason', 'midSeason', 'highSeason'))
df_sweep$EIR = as.numeric(sapply(df_sweep$Site, get_substr, index=2))
df_sweep$CM = as.numeric(sapply(df_sweep$Site, get_substr, index=3))
upper_sweep_ages = unique(df_sweep$Age)
df_sweep$mean_age = sapply(df_sweep$Age, get_mean_from_upper_age, upper_ages = upper_sweep_ages)

gg = ggplot(df_sweep, aes(x=mean_age, y=Incidence)) + 
  geom_line(aes(color=EIR, group=interaction(CM, EIR, season))) +
  scale_colour_gradient(
    low = rgb(0.2,0.2,0),
    high = rgb(0,0.8,0.9),
    trans='log',
    label = function(x) sprintf("%.1f", x)
  ) +
  facet_grid(CM ~ season) +
  coord_cartesian(xlim=c(0, 55)) +
  ylab('incidence') +
  xlab('age')+
  theme_bw()
ggsave(file.path(base_filepath_sim_output, '_plots', paste0('sweepSeasonEIRCM.png')), gg, width=7, height=4, units='in')



########################## plot comparisons with reference on background of sweep ####################
ref_df = read.csv(filepath_ref)
for (ss in sites){
  filepath_sim = file.path(paste0(base_filepath_sim_output, '/site_', ss, '/summary_data_final.csv'))
  df = read.csv(filepath_sim)
  upper_ages = sort(unique(df$Age))
  df$mean_age = sapply(df$Age, get_mean_from_upper_age, upper_ages = upper_ages)
  # get the sweep simulation that best matches the site simulation for this site
  df = get_matching_sweep_info(df, df_sweep, sites_cm_seek, study_site_closest_seasonality)
  
  reference = ref_df$INC[which(ref_df$REF_GROUP == as.numeric(ss))]
  reference = reference / 1000
  site_name = ref_df$ACD_LOCATION[ref_df$REF_GROUP == as.numeric(ss)][1]


  ref_cur = data.frame('Incidence'=reference, 'Age'=upper_ages)
  ref_cur$mean_age = sapply(ref_cur$Age, get_mean_from_upper_age, upper_ages = upper_ages)
  # repeat reference data for all CMs and seasonalities
  # ref_cur = merge(ref_cur, distinct(df_sweep[c('season', 'CM')]))

  
  gg2 = gg + 
    geom_line(data=df, aes(x=mean_age, y=Incidence), color='blue', size=1.1, alpha=0.5) + 
    geom_point(data=df, aes(x=mean_age, y=Incidence), color='blue', size=1.5, alpha=1) +
    geom_line(data=ref_cur, aes(x=mean_age, y=Incidence), color='red', size=1.1, alpha=0.75) + 
    geom_point(data=ref_cur, aes(x=mean_age, y=Incidence), color='red', size=1.5, alpha=1) 
      
  ggsave(file.path(base_filepath_sim_output, '_plots', paste0('sweepSeasonEIRCM_site', ss, '.png')), gg2, width=7, height=4, units='in')

}


########################## plot comparisons with reference without sweep background ####################
# set up reference dataset columns
ref_df = read.csv(filepath_ref)
ref_df = ref_df[which(ref_df$REF_GROUP %in% sites),]
ref_df$Incidence = ref_df$INC / 1000
ref_df$mean_age = NA
for (ss in sites){
  if (ss %in% ref_df$REF_GROUP){
    ref_df$mean_age[which(ref_df$REF_GROUP == as.numeric(ss))] = (ref_df$INC_LAR[which(ref_df$REF_GROUP == as.numeric(ss))] + ref_df$INC_UAR[which(ref_df$REF_GROUP == as.numeric(ss))]) / 2
  }
}
ref_df = data.frame('REF_GROUP'=ref_df$REF_GROUP, 'Incidence'=ref_df$Incidence, 'mean_age'=ref_df$mean_age, 'ACD_LOCATION'=ref_df$ACD_LOCATION)
ref_df$source = 'reference'

# set up simulation output columns
sim_site_df = data.frame('REF_GROUP'=NA, 'Incidence'=NA, 'mean_age'=NA, 'ACD_LOCATION'=NA)
for (ss in sites){
  filepath_sim = file.path(paste0(base_filepath_sim_output, '/site_', ss, '/summary_data_final.csv'))
  df = read.csv(filepath_sim)
  upper_ages = sort(unique(df$Age))
  df$mean_age = sapply(df$Age, get_mean_from_upper_age, upper_ages = upper_ages)
  
  df$ACD_LOCATION = ref_df$ACD_LOCATION[ref_df$REF_GROUP == as.numeric(ss)][1]
  
  sim_site_df = rbind(sim_site_df, data.frame('REF_GROUP'=df$Site, 'Incidence'=df$Incidence, 'mean_age'=df$mean_age, 'ACD_LOCATION'=df$ACD_LOCATION))
}
sim_site_df = sim_site_df[-1,]
sim_site_df$source = 'simulation'

df_combined = rbind(sim_site_df, ref_df)
df_combined$source = factor(df_combined$source, levels=c('reference', 'simulation'))
df_combined$site_name = paste0(df_combined$ACD_LOCATION, ' (id #', df_combined$REF_GROUP, ')')
df_combined_clean = df_combined[!(df_combined$site_name %in% paste0('Dielmo (id #', duplicate_sites, ')')),]

gg = ggplot(df_combined_clean, aes(x=mean_age, y=Incidence, color=source)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("reference" = "red",
                                "simulation"="blue")) +
  xlab('age (midpoint of age bin)') +
  ylab('incidence') +
  facet_wrap('site_name', ncol=4) +
  theme_bw()  
ggsave(file.path(base_filepath_sim_output, '_plots', paste0('compareSimRef_age_inc.png')), gg, width=9, height=5, units='in')

# format for comparison with Cameron
gg = ggplot(df_combined, aes(x=mean_age, y=Incidence, color=source)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("reference" = rgb(0.3, 0.3, 0.3),
                                "simulation"=rgb(0.7, 0.9, 0.1))) +
  scale_x_continuous(breaks=c(0, 5, 15, 45, 90), trans='sqrt') +
  coord_cartesian(xlim=c(-0.1, 90)) +  # 50
  facet_wrap('site_name', ncol=3, scales='free') +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA))
ggsave(file.path(base_filepath_sim_output, '_plots', paste0('compareSimRef_age_inc_CamFormat.png')), gg, width=6.5, height=4.75, units='in')


# zoom into ages <20
gg = ggplot(df_combined_clean, aes(x=mean_age, y=Incidence, color=source)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("reference" = "red",
                                "simulation"="blue")) +
  coord_cartesian(xlim=c(0, 20)) +
  xlab('age (midpoint of age bin)') +
  ylab('incidence') +
  facet_wrap('site_name', ncol=4) +
  theme_bw()
ggsave(file.path(base_filepath_sim_output, '_plots', paste0('compareSimRef_age_inc_zoomed.png')), gg, width=9, height=5, units='in')





########################## plot additional reference data not matched with simulation ####################
cm_match = 50
match_seasonality = FALSE
matchedString = ''
if (match_seasonality) matchedString='_matchSeason'
ref_df = read.csv(filepath_ref)
ref_df = ref_df[!is.na(ref_df$REF_GROUP),]
gg_list = list()
ii=1
ss_plot_list = unique(ref_df$REF_GROUP)  # 10:max(ref_df$REF_GROUP, na.rm=TRUE)
# remove the two Dielmo sites whose data is already included
ss_plot_list = ss_plot_list[!(ss_plot_list %in% duplicate_sites)]
# ss_plot_list = c(14, 15, 16, 17, 18, 21, 22, 24, 25, 27, 28, 29, 30, 31, 33, 34, 36, 37, 38)
for (ss in ss_plot_list){
  if (ss %in% ref_df$REF_GROUP){
    reference = ref_df$INC[which(ref_df$REF_GROUP == as.numeric(ss))]
    reference = reference / 1000
    site_name = gsub('\\?', '_', ref_df$ACD_LOCATION[ref_df$REF_GROUP == as.numeric(ss)][1])
    if(match_seasonality){
      cur_season = ref_df$SEASON_FOR_SIM[ref_df$REF_GROUP == as.numeric(ss)][1]
    } else cur_season='midSeason'
    
    mean_ages = (ref_df$INC_LAR[which(ref_df$REF_GROUP == as.numeric(ss))] + ref_df$INC_UAR[which(ref_df$REF_GROUP == as.numeric(ss))]) / 2
    ref_cur = data.frame('Incidence'=reference, 'mean_age'=mean_ages)
    ref_cur = distinct(ref_cur)
    
    # check that there are at least three data points that span at least 5 years
    if((nrow(ref_cur)>2) & ((max(ref_cur$mean_age)-min(ref_cur$mean_age))>5) & (max(ref_cur$Incidence)>0.5)){
      gg= ggplot(df_sweep[intersect(which(df_sweep$CM==cm_match), which(df_sweep$season==cur_season)),], aes(x=mean_age, y=Incidence)) + # , linetype="season"
        geom_line(aes(color=EIR, group=interaction(EIR))) +
        scale_colour_gradient(
          low = rgb(0.2,0.2,0),
          high = rgb(0,0.8,0.9),
          trans='log',
          label = function(x) sprintf("%.1f", x)
        ) +
        coord_cartesian(xlim=c(0, 55)) +
        xlab('age') +
        ylab('incidence') + 
        theme_bw()+
        theme(legend.position='none')
  
      gg2 = gg + 
        geom_line(data=ref_cur, aes(x=mean_age, y=Incidence), color='red', size=1.1, alpha=0.75) + 
        geom_point(data=ref_cur, aes(x=mean_age, y=Incidence), color='red', size=2, alpha=1)  +
        annotate("text",  x=Inf, y = Inf, label = site_name, vjust=1.2, hjust=1.1, color=rgb(0.5,0.5,0.5))
      gg_list[[ii]] = gg2
      ii=ii+1
    }
   }
 }

do.call("grid.arrange", c(gg_list, ncol=ceiling(sqrt(length(gg_list)))))
dev.copy(png, file.path(base_filepath_sim_output, '_plots',paste0("sweepSeasonEIRCM_otherSites", matchedString, ".png")), width=16, height=9, units='in', res=900)
# dev.copy(png, file.path(base_filepath_sim_output, '_plots',"sweepSeasonEIRCM_allRefSites", matchedString, ".png"), width=15, height=10, units='in', res=900)
dev.off()




########################## plot all reference data ####################
ref_df = read.csv(filepath_ref)
ref_df$EIR_FOR_SIM = ref_df$EIR_FOR_SIM_V3
ref_df = ref_df[!is.na(ref_df$REF_GROUP),]
ref_df$Incidence = ref_df$INC / 1000
ref_df$EIR_FOR_SIM = as.numeric(ref_df$EIR_FOR_SIM)
ref_df$SEASON_FOR_SIM = factor(ref_df$SEASON_FOR_SIM, levels=c('flatSeason', 'midSeason', 'highSeason'))
ref_df$mean_age = (ref_df$INC_LAR + ref_df$INC_UAR) / 2
ref_df = ref_df[!(ref_df$REF_GROUP %in% duplicate_sites), ]
ss_plot_list = unique(ref_df$REF_GROUP)  # 10:max(ref_df$REF_GROUP, na.rm=TRUE)

# gg=ggplot(ref_df, aes(x=mean_age, y=Incidence, color=EIR_FOR_SIM, group=interaction(REF_GROUP, EIR_FOR_SIM, SEASON_FOR_SIM)), alpha=0.5) +
#   geom_line() +
#   geom_point(size=1.5, alpha=0.3) +
#   scale_colour_gradient(
#     low = rgb(0.2,0.2,0),
#     high = rgb(0,0.8,0.9),
#     trans='log',
#     label = function(x) sprintf("%.1f", x)
#   ) +
#   xlab('age') +
#   ylab('incidence') + 
#   theme_bw() +
#   facet_grid('SEASON_FOR_SIM')#+
# ggsave(file.path(base_filepath_sim_output, '_plots', paste0('reference_age_inc_allSites.png')), gg, width=5, height=6, units='in')


# group by mean incidence among those under 10
ref_df$inc_group = NA
inc_groups = c(0, 1, 2, 3, 4)
inc_group_levels = paste0('mean U10 inc: ', c('0-1', '1-2', '2-3', '3-4', '>4'))
for (ss in ss_plot_list){
  if (ss %in% ref_df$REF_GROUP){
    mean_u10_inc = mean(ref_df$Incidence[intersect(which(ref_df$REF_GROUP == as.numeric(ss)), which(ref_df$INC_UAR<=10))], na.rm=TRUE)
    cur_inc_group = inc_group_levels[max(which(inc_groups<mean_u10_inc))]
    ref_df$inc_group[which(ref_df$REF_GROUP == as.numeric(ss))] = cur_inc_group
  }
}
ref_df$inc_group = factor(ref_df$inc_group, levels=inc_group_levels)
gg=ggplot(ref_df, aes(x=mean_age, y=Incidence, color=EIR_FOR_SIM, group=interaction(REF_GROUP, EIR_FOR_SIM, SEASON_FOR_SIM)), alpha=0.5) +
  geom_line() +
  geom_point(size=1.5, alpha=0.3) +
  scale_colour_gradient(
    low = rgb(0.2,0.2,0),
    high = rgb(0,0.8,0.9),
    trans='log',
    label = function(x) sprintf("%.1f", x)
  ) +
  xlab('age') +
  ylab('incidence') + 
  theme_bw() +
  facet_grid(vars(inc_group), vars(SEASON_FOR_SIM))
  # facet_grid(vars(SEASON_FOR_SIM), vars(inc_group))
ggsave(file.path(base_filepath_sim_output, '_plots', paste0('reference_age_inc_allSites_facetMeanU10Inc.png')), gg, width=8, height=5, units='in')


# group by peak incidence among any sampled age group
ref_df$inc_group = NA
inc_groups = c(0, 1.5, 3, 5)
inc_group_levels = paste0('max inc: ', c('<1.5', '1.5-3', '3-5', '>5'))
for (ss in ss_plot_list){
  if (ss %in% ref_df$REF_GROUP){
    peak_inc = max(ref_df$Incidence[which(ref_df$REF_GROUP == as.numeric(ss))], na.rm=TRUE)
    cur_inc_group = inc_group_levels[max(which(inc_groups<peak_inc))]
    ref_df$inc_group[which(ref_df$REF_GROUP == as.numeric(ss))] = cur_inc_group
  }
}
ref_df$inc_group = factor(ref_df$inc_group, levels=inc_group_levels)
# gg=ggplot(ref_df, aes(x=mean_age, y=Incidence, color=EIR_FOR_SIM, group=interaction(REF_GROUP, EIR_FOR_SIM, SEASON_FOR_SIM)), alpha=0.5) +
gg=ggplot(ref_df, aes(x=mean_age, y=Incidence, color=as.factor(COUNTRY), group=interaction(REF_GROUP, EIR_FOR_SIM, SEASON_FOR_SIM)), alpha=0.5) +
  geom_line() +
  geom_point(size=1.5, alpha=0.3) +
  # scale_colour_gradient(
  #   low = rgb(0.2,0.2,0),
  #   high = rgb(0,0.8,0.9),
  #   trans='log',
  #   label = function(x) sprintf("%.1f", x)
  # ) +
  xlab('age (midpoint of age bin)') +
  ylab('incidence') + 
  theme_bw() +
  # theme(legend.position='none') + 
  facet_grid(vars(inc_group), vars(SEASON_FOR_SIM))
# facet_grid(vars(SEASON_FOR_SIM), vars(inc_group))
ggsave(file.path(base_filepath_sim_output, '_plots', paste0('reference_age_inc_allSites_facetMaxInc.png')), gg, width=9, height=5.5, units='in')






########### plot monthly EIRs for sites and sweeps  #################

# monthly EIRs used in simulation sweeps
sweep_monthly_EIRs = data.frame('month' = 1:12)
# explore results of different EIRs and seasonalities on age-incidence relationship
highSeason = c(1, 1, 1, 1, 1, 3, 27, 70, 130, 57, 29, 1) / 322  # Dapelogo
midSeason = c(1, 1, 1, 1, 1, 7, 8, 9, 5, 4, 6, 1) / 45  # Laye
flatSeason = rep(1.0 / 12, 12)
seasonalities = list(highSeason, midSeason, flatSeason)
seasonality_names = c('highSeason', 'midSeason', 'flatSeason')
for (ss in 1:length(seasonality_names)){
  for (ee in eirs){
    sweep_monthly_EIRs[[paste0(seasonality_names[ss],'_', ee)]] = ee * seasonalities[[ss]]
  }
}




# comparison of simulation sweeps and site-specific EIRs
plot_sweep_value = TRUE
par(mfrow=c(3,3))
for(ss in 1:9){
  cur_eir = sum(study_site_monthly_EIRs[[paste0('X', ss)]])
  closest_eir = eirs[which.min(abs(eirs-cur_eir))]
  plot(study_site_monthly_EIRs[[paste0('X', ss)]], type='l', lwd=2, main=paste0(ss, ' (EIR=',cur_eir,')'), xlab='month', ylab='EIR', bty='L')
  if (plot_sweep_value){
    lines(sweep_monthly_EIRs[[paste0(seasonality_names[1],'_', closest_eir)]], col='blue')
    lines(sweep_monthly_EIRs[[paste0(seasonality_names[2],'_', closest_eir)]], col='green')
    lines(sweep_monthly_EIRs[[paste0(seasonality_names[3],'_', closest_eir)]], col='red')
    if(ss==1) legend('topleft', c('high', 'mid', 'flat'), col=c('blue', 'green', 'red'), lwd=1, bty='n')
  }
}
par(mfrow=c(1,1))



# plot comparisons of all sweep and site monthly EIRs
study_site_monthly_EIRs_wide = melt(study_site_monthly_EIRs, id.vars=c('month'))
study_site_monthly_EIRs_wide$group = 'ref_site'
sweep_monthly_EIRs_wide = melt(sweep_monthly_EIRs, id.vars=c('month'))
sweep_monthly_EIRs_wide$group = 'sweep'
sweep_monthly_EIRs_wide$season = sapply(as.character(sweep_monthly_EIRs_wide$variable), get_substr, index=1)
sweep_monthly_EIRs_wide$season = factor(sweep_monthly_EIRs_wide$season, levels=c('flatSeason', 'midSeason', 'highSeason'))
sweep_monthly_EIRs_wide$EIR = as.numeric(sapply(as.character(sweep_monthly_EIRs_wide$variable), get_substr, index=2))



df_EIRs = rbind(study_site_monthly_EIRs_wide, sweep_monthly_EIRs_wide)

ggplot(df_EIRs, aes(x=month, y=value, color=variable, linetype=group)) + 
  geom_line()

ggplot(study_site_monthly_EIRs_wide, aes(x=month, y=value, color=variable)) + 
  geom_line()

gg = ggplot(sweep_monthly_EIRs_wide, aes(x=month, y=value, color=EIR, group=interaction(EIR, season))) + 
  geom_line()+
  scale_colour_gradient(
    low = rgb(0.2,0.2,0),
    high = rgb(0,0.8,0.9),
    trans='log',
    label = function(x) sprintf("%.1f", x)
  ) +
  ylab('monthly EIR') +
  facet_grid(~season) +
  theme_bw()
ggsave(file.path(base_filepath_sim_output, '_plots', paste0('season_EIR_sweep_monthlyEIRs.png')), gg, width=6.5, height=2.5, units='in')
