# extract infectiousness data from python hard-coded file


###############################################################################
# function to convert nested lists into dataframe
#   function expects as input:
#      - raw_ref_data: nested list containing the number of individuals in each {month, age-bin, parasitemia-bin} to infect a certain fraction of mosquitos
#          In the list, the outermost list level is month, then age, then parasitemia bin, then innermost is fraction_infected_bin, with the value being the number of individuals in that bin
#      - months, age_bins, parasitemia_bins, fraction_infected_bins: vectors containing the bin values for each factor (dimensions must correspond to list)
create_df_from_nested_lists = function(raw_ref_data, months, age_bins, parasitemia_bins, fraction_infected_bins){
  df = data.frame()
  for(mm in 1:length(months)){
    month = months[mm]
    for(aa in 1:length(age_bins)){
      age = age_bins[aa]
      for(pp in 1:length(parasitemia_bins)){
        par_dens = parasitemia_bins[pp]
        
        # get a total of the number of people in each month-age-par_dens group (summed across fraction_infected_bins)
        num_in_group = 0
        for(ff in 1:length(fraction_infected_bins)){
          frac_infected_bin = fraction_infected_bins[ff]
          count =  raw_ref_data[[mm]][[aa]][[pp]][[ff]]
          num_in_group = num_in_group + count
          df = rbind(df, data.frame(month=month, agebin=age, densitybin=par_dens, fraction_infected_bin=frac_infected_bin, count=count, num_in_group=NA))
        }
        df$num_in_group[(df$month==month) & (df$agebin==age) & (df$densitybin==par_dens)] = num_in_group
      }
    }
  }
  # calculate the fraction of individuals in each month-age-par_dens group that fell in each fraction-infected bin
  df$freq_frac_infect = df$count / df$num_in_group
  return(df)
}



###############################################################################
output_dir = '/Users/moniqueam/Documents/malaria-model_validation/reference_datasets'
# Dapelogo and Laye - smeared infectiousness by smeared gametocytemia and age bin
months = c(7, 9, 1)
age_bins = c(5, 15, 80)  # (, 5] (5, 15] (15, ]
parasitemia_bins = c(0, 0.5, 5, 50, 500, 5000, 50000, 500000)  # (, 0] (0, 50] ... (50000, ]
fraction_infected_bins = c(0, 5, 20, 50, 80, 100)

###############################################################################
# Dapelogo data were grabbed from https://github.com/InstituteforDiseaseModeling/idmtools-calibra/blob/0d8aa8f6ef867fe9d2823dc9536eceec3759bdd7/malaria/study_sites/dapelogo_infectiousness_calib_site.py
#      replace innermost "[" with "c(", all other "[" with "list(" and all "]" with ")"
# outermost list level is month, then age, then parasitemia bin, then innermost is fraction_infected_bin
dapelogo_raw_ref_data = list(
  # month = 7
  list(list(c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0),
       c(1, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0),
       c(1, 1, 0, 2, 0, 0), c(0, 1, 0, 1, 0, 0),
       c(0, 0, 0, 1, 1, 0),
       c(0, 0, 0, 0, 0, 0)),
      list(c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0),
       c(0, 0, 0, 0, 0, 0), c(0, 0, 1, 0, 0, 0),
       c(0, 2, 1, 0, 1, 0), c(1, 1, 3, 1, 0, 0),
       c(0, 0, 0, 0, 0, 0),
       c(0, 0, 0, 0, 0, 0)),
      list(c(7, 5, 0, 0, 0, 0), c(3, 0, 0, 0, 0, 0),
       c(2, 0, 0, 0, 0, 0), c(2, 1, 0, 0, 0, 0),
       c(2, 1, 0, 0, 0, 0), c(2, 2, 0, 0, 0, 0),
       c(0, 0, 0, 0, 0, 0),
       c(0, 0, 0, 0, 0, 0))),
 # month = 9
 list(list(c(0, 0, 1, 0, 0, 0), c(1, 0, 0, 0, 0, 0),
       c(0, 0, 0, 0, 0, 0), c(1, 1, 0, 0, 0, 0),
       c(1, 0, 1, 0, 0, 0), c(0, 0, 0, 2, 0, 0),
       c(0, 0, 0, 0, 0, 0),
       c(0, 0, 0, 0, 0, 0)),
      list(c(4, 0, 1, 0, 0, 0), c(0, 0, 0, 0, 0, 0),
       c(2, 0, 0, 0, 0, 0), c(0, 0, 1, 0, 2, 1),
       c(0, 0, 0, 0, 0, 0), c(1, 0, 0, 0, 1, 1),
       c(0, 0, 0, 0, 0, 0),
       c(0, 0, 0, 0, 0, 0)),
      list(c(13, 1, 0, 0, 0, 0),
       c(4, 0, 1, 0, 0, 0), c(1, 0, 0, 0, 0, 0),
       c(0, 0, 0, 0, 1, 0), c(1, 0, 1, 0, 0, 0),
       c(1, 0, 0, 1, 0, 0), c(0, 0, 0, 0, 0, 0),
       c(0, 0, 0, 0, 0, 0))),
 # month = 1
  list(list(c(2, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0),
       c(1, 0, 0, 0, 0, 0), c(1, 0, 0, 1, 0, 0),
       c(1, 0, 0, 1, 0, 0), c(0, 0, 0, 1, 0, 1),
       c(0, 0, 1, 0, 0, 0),
       c(0, 0, 0, 0, 0, 0)),
      list(c(2, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0),
       c(1, 0, 0, 1, 0, 0), c(3, 0, 0, 0, 0, 0),
       c(2, 0, 2, 0, 0, 0), c(1, 1, 0, 0, 0, 0),
       c(1, 0, 0, 0, 0, 0),
       c(0, 0, 0, 0, 0, 0)),
      list(c(13, 1, 0, 0, 0, 0),
       c(2, 0, 0, 0, 0, 0), c(1, 1, 0, 0, 0, 0),
       c(3, 0, 0, 0, 0, 0), c(4, 0, 0, 0, 0, 0),
       c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0),
       c(0, 0, 0, 0, 0, 0))))

df_dapelogo = create_df_from_nested_lists(dapelogo_raw_ref_data, months, age_bins, parasitemia_bins, fraction_infected_bins)
df_dapelogo$site = 'dapelogo_2007'
write.csv(df_dapelogo, paste0(output_dir, '/dapelogo_infectiousness.csv'))



###############################################################################
# Laye data were grabbed from https://github.com/InstituteforDiseaseModeling/idmtools-calibra/blob/0d8aa8f6ef867fe9d2823dc9536eceec3759bdd7/malaria/study_sites/laye_infectiousness_calib_site.py
#      replace innermost "[" with "c(", all other "[" with "list(" and all "]" with ")"
# outermost list level is month, then age, then parasitemia bin, then innermost is fraction_infected_bin
laye_raw_ref_data = 
  list(
    # month = 7
    list(
      list(c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0),
           c(1, 0, 2, 0, 1, 1), c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0)),
      list(c(1, 2, 0, 0, 0, 0), c(1, 0, 1, 1, 0, 0), c(1, 2, 0, 0, 0, 0), c(1, 1, 0, 0, 1, 0), c(3, 0, 2, 3, 0, 0),
           c(0, 0, 1, 0, 0, 0), c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0)),
      list(c(15, 1, 0, 0, 0, 0), c(1, 1, 0, 0, 0, 0), c(1, 0, 0, 0, 0, 0), c(2, 0, 0, 0, 0, 0), c(3, 0, 0, 1, 0, 1),
           c(0, 0, 1, 0, 0, 0), c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0))),
    # month = 9
    list(
      list(c(1, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0), c(1, 0, 0, 0, 0, 0),
           c(0, 0, 1, 0, 0, 0), c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0)),
      list(c(2, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0), c(4, 0, 0, 0, 0, 0), c(4, 0, 1, 2, 1, 0),
           c(1, 1, 0, 1, 1, 0), c(1, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0)),
      list(c(7, 0, 0, 0, 0, 0), c(3, 0, 0, 0, 0, 0), c(1, 0, 1, 0, 0, 0), c(4, 0, 1, 0, 0, 0), c(4, 1, 0, 0, 0, 0),
           c(1, 1, 0, 0, 1, 0), c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0))),
    # month = 1
    list(
      list(c(1, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0),
           c(0, 0, 0, 1, 0, 0), c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0)),
      list(c(7, 0, 0, 0, 0, 0), c(1, 0, 0, 0, 0, 0), c(3, 1, 1, 0, 0, 0), c(3, 0, 0, 0, 0, 0), c(2, 1, 0, 0, 0, 0),
           c(1, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0)),
      list(c(13, 0, 0, 0, 0, 0), c(3, 0, 0, 0, 0, 0), c(2, 0, 0, 0, 0, 0), c(6, 0, 0, 0, 0, 0), c(2, 1, 0, 0, 0, 0),
           c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0)))
  )

df_laye = create_df_from_nested_lists(laye_raw_ref_data, months, age_bins, parasitemia_bins, fraction_infected_bins)
df_laye$site = 'laye_2007'
write.csv(df_laye, paste0(output_dir, '/laye_infectiousness.csv'))

