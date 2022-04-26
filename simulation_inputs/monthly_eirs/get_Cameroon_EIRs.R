# get monthly eir from sampled days using spline versus using average of sample days in each month

library(stats)
library(lubridate)

############# set filepaths and annual EIRs #################
data_base_filepath = "/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/references/EIR_data_extract/"
# filename = "Meunier_data_extract_Ebolakounou.csv"
filename = "Meunier_data_extract_Koundou.csv"
annual_eir_koundou = 176.1
annual_eir_ebolakounou = 17.7

############# approach 1 for getting monthly EIR (linear interpolation, followed by averaging within month) #################
day_eir = read.csv(paste0(data_base_filepath, filename))
day_eir$day = as.Date(day_eir$day, format='%m/%d/%Y')
day_eir$nightly_EIR = sapply(day_eir$nightly_EIR, max, 0.1)
# tack on an extra couple year repeats to get some form of continuity
day_eir2 = day_eir
day_eir2$day = day_eir2$day + 365
day_eir = rbind(day_eir, day_eir2)
day_eir2$day = day_eir2$day + 365
day_eir = rbind(day_eir, day_eir2)
day_eir$doy = as.numeric(day_eir$day - min(day_eir$day))

start_doy = as.numeric(as.Date('1998-01-01') - min(day_eir$day))
end_doy = start_doy + 365
approx_month_days = round(seq(start_doy, end_doy, 30.4))
approx_eir = rep(NA, 12)

ap = approxfun(day_eir$doy, day_eir$nightly_EIR, method='linear')
for(ii in 1:12){
  approx_eir[ii] = mean(ap(approx_month_days[ii]:(approx_month_days[ii+1]-1)))
}


############# approach 2 for getting monthly EIR (average two samples within month) #################
day_eir = read.csv(paste0(data_base_filepath, filename))
day_eir$day = as.Date(day_eir$day, format='%m/%d/%Y')
day_eir$nightly_EIR = sapply(day_eir$nightly_EIR, max, 0.1)
approx_eir_v2 = rep(NA, 12)
for(mm in 1:12){
  approx_eir_v2[mm] = mean(day_eir$nightly_EIR[month(day_eir$day) == mm])
}
plot(day_eir$day, day_eir$nightly_EIR, type='p', pch=20, cex=2)
for(mm in 1:12){
  lines(x=c(as.Date(paste0('1997-', mm, '-01')), as.Date(paste0('1997-', mm, '-28'))), y=rep(approx_eir_v2[mm],2), col='purple')
  lines(x=c(as.Date(paste0('1998-', mm, '-01')), as.Date(paste0('1998-', mm, '-28'))), y=rep(approx_eir_v2[mm],2), col='purple')
  lines(x=c(as.Date(paste0('1997-', mm, '-01')), as.Date(paste0('1997-', mm, '-28'))), y=rep(approx_eir[mm],2), col='green')
  lines(x=c(as.Date(paste0('1998-', mm, '-01')), as.Date(paste0('1998-', mm, '-28'))), y=rep(approx_eir[mm],2), col='green')
}
legend('topleft', c('mean of two sample dates in month', 'linear interpolation, then mean of all days in month'), lwd=1, col=c('purple', 'green'), bty='n')


############# rescale monthly EIRs to get appropriate annual EIR #################
rescaled_koundou = approx_eir_v2 / sum(approx_eir_v2) * annual_eir_koundou
print(data.frame(month=1:12, eir=rescaled_koundou))
plot(rescaled_koundou, type='b')

rescaled_ebolakounou = approx_eir_v2 / sum(approx_eir_v2) * annual_eir_ebolakounou
print(data.frame(month=1:12, eir=rescaled_ebolakounou))

