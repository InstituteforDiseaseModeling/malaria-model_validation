
library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)
library(data.table)

############################ setup: filepaths and info on simulations run #################################
filepath_ref = "/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/reference_data/Garki_df.csv"
included_villages = c('Rafin Marke', 'Sugungum', 'Matsari')  # set as NA to keep all villages in dataset
output_dir = '/Users/moniqueam/Documents/malaria-model_validation/reference_datasets'
parasite_density_bins = c(0, 50, 500, 5000, 5000000)

# convert the numbers in the Garki reference to parasite density estimates
num_fields = 200
total_blood_volume = 0.5
positive_field_fraction_bins = (1-exp(-total_blood_volume/num_fields*parasite_density_bins))


ref_df = fread(filepath_ref)
ref_df$par_pos = ref_df$Parasitemia>0
if (!is.na(included_villages[1])){
  ref_df = ref_df[ref_df$Village %in% included_villages,]
}

age_bin_max = rep(NA, length(unique(ref_df$`Age Bins`)))
for (ii in 1:length(age_bin_max)){
  age_bin_max[ii] = max(ref_df$Age[ref_df$`Age Bins`==ii])
}
age_bin_max = ceiling(age_bin_max)
ref_df$age_bin_year = age_bin_max[ref_df$`Age Bins`]

# order factors for seasons and assign months
#  month-season relationships are: {'1': 'DC', '2': 'DC', '3': 'DC', '4': 'DC', '5': 'DH', '6': 'DH', '7': 'DH', '8': 'W', '9': 'W', '10': 'W', '11': 'W', '12': 'W'}
season_types = c('DC', 'DH', 'W')  
season_corr_months = c(2, 7, 10)
include_years = 2
all_years = 2:7
season_levels = paste0(rep(season_types, length(all_years)), rep(all_years, each=3))
ref_df$Seasons = factor(ref_df$Seasons, levels=season_levels)



# ============================================================ #
#                get prevalence by age bin
# ============================================================ #

ref_sum = ref_df %>% group_by(age_bin_year, Village, Seasons) %>%
  summarise(total_sampled = n(),
            num_pos = sum(par_pos))
ref_sum$prevalence = ref_sum$num_pos / ref_sum$total_sampled

ggplot(ref_sum, aes(x=age_bin_year, y=prevalence)) +
  geom_point() + 
  geom_line() + 
  facet_grid(Village~Seasons)

# add site name used for simulations
ref_sum$Site = paste0(ref_sum$Village, '_1970')
ref_sum$Site = gsub(' ', '_', ref_sum$Site)
colnames(ref_sum)[which(colnames(ref_sum) =='age_bin_year')] = 'agebin'

# subset to included year and switch to using month to match simulations
ref_sum$year = as.numeric(gsub("\\D", "", ref_sum$Seasons))
ref_sum$month0 = gsub("([0-9]+).*", "", ref_sum$Seasons)
ref_sum$month = season_corr_months[match(ref_sum$month0, season_types)]
ref_sum_subset = ref_sum[ref_sum$year %in% include_years,]

write.csv(ref_sum_subset, paste0(output_dir,'/garki_prev_by_age_bin.csv'), row.names = FALSE)




# ============================================================ #
#            get parasite density bins by age bin
# ============================================================ #

# add columns for parasite density bins
ref_df$asexual_bin_num = findInterval(ref_df$Parasitemia, c(0, 10^-6 + positive_field_fraction_bins))
# ref_df$asexual_bin_num = findInterval(ref_df$Parasitemia, c(-10^-6, positive_field_fraction_bins), left.open=TRUE)
ref_df$asexual_field_bin = positive_field_fraction_bins[ref_df$asexual_bin_num]
ref_df$asexual_dens_bin = parasite_density_bins[ref_df$asexual_bin_num]
ref_df$gamet_bin_num = findInterval(ref_df$Gametocytemia, c(0, 10^-6 + positive_field_fraction_bins))
ref_df$gamet_dens_bin = parasite_density_bins[ref_df$gamet_bin_num]

ref_asex_sum = ref_df %>% group_by(age_bin_year, asexual_dens_bin, Village, Seasons) %>%
  summarise(count_asex = n())
ref_gamet_sum = ref_df %>% group_by(age_bin_year, gamet_dens_bin, Village, Seasons) %>%
  summarise(count_gamet = n())
ref_age_totals = ref_df %>% group_by(age_bin_year, Village, Seasons) %>%
  summarise(num_in_age = n())
ref_asex_sum = merge(ref_asex_sum, ref_age_totals)
ref_asex_sum$asexual_par_dens_freq = ref_asex_sum$count_asex / ref_asex_sum$num_in_age
colnames(ref_asex_sum)[which(colnames(ref_asex_sum) =='asexual_dens_bin')] = 'densitybin'
colnames(ref_asex_sum)[which(colnames(ref_asex_sum) =='num_in_age')] = 'bin_total_asex'

ref_gamet_sum = merge(ref_gamet_sum, ref_age_totals)
ref_gamet_sum$gametocyte_dens_freq = ref_gamet_sum$count_gamet / ref_gamet_sum$num_in_age
colnames(ref_gamet_sum)[which(colnames(ref_gamet_sum) =='gamet_dens_bin')] = 'densitybin'
colnames(ref_gamet_sum)[which(colnames(ref_gamet_sum) =='num_in_age')] = 'bin_total_gamet'

ref_dens_agg = merge(ref_asex_sum, ref_gamet_sum)


# add site name used for simulations
ref_dens_agg$Site = paste0(ref_dens_agg$Village, '_1970')
ref_dens_agg$Site = gsub(' ', '_', ref_dens_agg$Site)
colnames(ref_dens_agg)[which(colnames(ref_dens_agg) =='age_bin_year')] = 'agebin'

# subset to included year and switch to using month to match simulations
ref_dens_agg$year = as.numeric(gsub("\\D", "", ref_dens_agg$Seasons))
ref_dens_agg$month0 = gsub("([0-9]+).*", "", ref_dens_agg$Seasons)
ref_dens_agg$month = season_corr_months[match(ref_dens_agg$month0, season_types)]
ref_dens_agg_subset = ref_dens_agg[ref_dens_agg$year %in% include_years,]

write.csv(ref_dens_agg_subset, paste0(output_dir,'/garki_par_dens_by_age_bin.csv'), row.names = FALSE)



