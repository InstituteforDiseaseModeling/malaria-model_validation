# plot_validation_age_inc_sim_sweep_ref.R

# first, plot the simulation relationships between age and incidence for a sweep of EIRs, seasonality patterns, and CM rates
# second, plot the match between the reference datasets, matched simulations, and sweep of EIRs

library(ggplot2)
library(reshape2)



num_seeds = 2
base_filepath_sim_output = 'c:/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/simulation_output/ewan_age_incidence'
filepath_sim = file.path(base_filepath_sim_output, 'sweepSeasonEIRCM', 'summary_data_final.csv')

get_substr = function(site_name_str, index){
  strsplit(site_name_str, "_")[[1]][index]
}
get_mean_from_upper_age = function(cur_age, upper_ages){
  mean_ages = (c(0, upper_ages[1:(length(upper_ages)-1)]) + upper_ages) / 2
  return(mean_ages[which(upper_ages == cur_age)])
}

df_sweep = read.csv(filepath_sim)
df_sweep$season = sapply(df_sweep$Site, get_substr, index=1)
df_sweep$season = factor(df_sweep$season, levels=c('flatSeason', 'midSeason', 'highSeason'))
df_sweep$EIR = as.numeric(sapply(df_sweep$Site, get_substr, index=2))
df_sweep$CM = as.numeric(sapply(df_sweep$Site, get_substr, index=3))
upper_sweep_ages = unique(df_sweep$Age)
df_sweep$mean_age = sapply(df_sweep$Age, get_mean_from_upper_age, upper_ages = upper_sweep_ages)

gg = ggplot(df_sweep, aes(x=mean_age, y=Incidence, color=EIR, group=interaction(CM, EIR, season))) + # , linetype="season"
  geom_line() +
  scale_colour_gradient(
    low = rgb(0.2,0.2,0),
    high = rgb(0,0.8,0.9),
    trans='log'
  ) +
  facet_grid(CM ~ season) +
  theme_bw()
ggsave(file.path(base_filepath_sim_output, '_plots', paste0('sweepSeasonEIRCM_v2.png')), gg)


sites = c('1', '2', '3', '4', '5', '7', '8', '9')   # '6',
EIR_sites = c()
seasonality_sites = c()
num_seeds = 2
base_filepath_sim_output = '/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/simulation_output/ewan_age_incidence'
filepath_ref = "/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/reference_data/Cameron_age_incidence_eir.csv"
ref_df = read.csv(filepath_ref)

for (ss in sites){
  filepath_sim = file.path(paste0(base_filepath_sim_output, '/site_', ss, '/summary_data_final.csv'))
  df = read.csv(filepath_sim)
  upper_ages = sort(unique(df$Age))
  df$mean_age = sapply(df$Age, get_mean_from_upper_age, upper_ages = upper_ages)
  
  reference = ref_df$INC[which(ref_df$REF_GROUP == as.numeric(ss))]
  reference = reference / 1000
  site_name = ref_df$ACD_LOCATION[ref_df$REF_GROUP == as.numeric(ss)][1]
  # reference = [3.2, 5, 6.1, 4.75, 3.1, 2.75, 2.7, 1.9, 0.12, 0.8, 0.5, 0.25, 0.1, 0.2, 0.4,
  #              0.3, 0.2, 0.2, 0.2, 0.15, 0.15, 0.15]
  # sample_size = [55, 60, 55, 50, 50, 38, 38, 38, 38, 38, 26, 26, 26, 26, 26, 110, 75, 75, 150, 70, 70, 90]
  # error = np.array(reference)/np.sqrt(np.array(sample_size))

  ref_cur = data.frame('Incidence'=reference, 'Age'=upper_ages)
  ref_cur$mean_age = sapply(ref_cur$Age, get_mean_from_upper_age, upper_ages = upper_ages)
  # ggplot() + 
  #   geom_line(data=df, aes(x=Age, y=Incidence), color='blue', size=3) + 
  #   geom_line(data=ref_cur, aes(x=Age, y=Incidence), color='red', size=3) + 
  #   geom_point(data=ref_cur, aes(x=Age, y=Incidence), color='red', size=3) 
    
    
    
  
    
gg=  ggplot(df_sweep, aes(x=mean_age, y=Incidence)) + # , linetype="season"
    geom_line(aes(color=EIR, group=interaction(CM, EIR, season))) +
    scale_colour_gradient(
      low = rgb(0.2,0.2,0),
      high = rgb(0,0.8,0.9),
      trans='log'
    ) +
    facet_grid(CM ~ season) +
    theme_bw()

gg2 = gg + 
  geom_line(data=df, aes(x=mean_age, y=Incidence), color='blue', size=3) + 
  geom_line(data=ref_cur, aes(x=mean_age, y=Incidence), color='red', size=3) + 
  geom_point(data=ref_cur, aes(x=mean_age, y=Incidence), color='red', size=3) 
    
ggsave(file.path(base_filepath_sim_output, '_plots', paste0('sweepSeasonEIRCM_site', ss, '.png')), gg2)

}









# plot monthly EIRs for sites and seasonalities

# monthly EIRs used in simulation sweeps
sweep_monthly_EIRs = data.frame('month' = 1:12)
# explore results of different EIRs and seasonalities on age-incidence relationship
highSeason = c(1, 1, 1, 1, 1, 3, 27, 70, 130, 57, 29, 1) / 322  # Dapelogo
midSeason = c(1, 1, 1, 1, 1, 7, 8, 9, 5, 4, 6, 1) / 45  # Laye
flatSeason = rep(1.0 / 12, 12)
seasonalities = list(highSeason, midSeason, flatSeason)
seasonality_names = c('highSeason', 'midSeason', 'flatSeason')
eirs = c(10, 20, 30, 50, 100, 200)
for (ss in 1:length(seasonality_names)){
  for (ee in eirs){
    sweep_monthly_EIRs[[paste0(seasonality_names[ss],'_', ee)]] = ee * seasonalities[[ss]]
  }
}


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

# comparison of simulation sweeps and site-specific EIRs
par(mfrow=c(3,3))
for(ss in 1:9){
  cur_eir = sum(study_site_monthly_EIRs[[paste0('X', ss)]])
  closest_eir = eirs[which.min(abs(eirs-cur_eir))]
  plot(study_site_monthly_EIRs[[paste0('X', ss)]], type='l', lwd=2, main=ss, xlab='month', ylab='EIR', bty='L')
  lines(sweep_monthly_EIRs[[paste0(seasonality_names[1],'_', closest_eir)]], col='blue')
  lines(sweep_monthly_EIRs[[paste0(seasonality_names[2],'_', closest_eir)]], col='green')
  lines(sweep_monthly_EIRs[[paste0(seasonality_names[3],'_', closest_eir)]], col='red')
  if(ss==1) legend('topleft', c('high', 'mid', 'flat'), col=c('blue', 'green', 'red'), lwd=1, bty='n')
}
par(mfrow=c(1,1))

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


# plot comparisons of all sweep and site monthly EIRs
study_site_monthly_EIRs_wide = melt(study_site_monthly_EIRs, id.vars=c('month'))
study_site_monthly_EIRs_wide$group = 'ref_site'
sweep_monthly_EIRs_wide = melt(sweep_monthly_EIRs, id.vars=c('month'))
sweep_monthly_EIRs_wide$group = 'sweep'

df_EIRs = rbind(study_site_monthly_EIRs_wide, sweep_monthly_EIRs_wide)

ggplot(df_EIRs, aes(x=month, y=value, color=variable, linetype=group)) + 
  geom_line()

ggplot(study_site_monthly_EIRs_wide, aes(x=month, y=value, color=variable)) + 
  geom_line()

ggplot(sweep_monthly_EIRs_wide, aes(x=month, y=value, color=variable)) + 
  geom_line()
