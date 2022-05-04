# reformat reference duration-of-infection dataset from Navrongo, Ghana 


dropbox_filepath = "C:/Users/moniqueam/Dropbox (IDM)/Malaria Team Folder"
output_dir = '/Users/moniqueam/Documents/malaria-model_validation/reference_datasets'

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
ref_data$age = round(ref_data$age,1)

ref_data = ref_data[,c('SID', 'date', 'DENSITY', 'age')]
ref_data$site = 'navrongo_2000'

write.csv(ref_data, paste0(output_dir, '/Navrongo_infection_duration.csv'), row.names=FALSE)


