import pandas as pd
import numpy as np
import os
import simulations.manifest as manifest


# commission simulations with a sweep across different EIRs, CM rates, and seasonalities
highSeason = [xx/322 for xx in [1, 1, 1, 1, 1, 3, 27, 70, 130, 57, 29, 1]]  # Dapelogo
midSeason = [xx/45 for xx in [1, 1, 1, 1, 1, 7, 8, 9, 5, 4, 6, 1]]  # Laye
flatSeason = [1.0 / 12] * 12
seasonalities = [highSeason, midSeason, flatSeason]
seasonality_names = ['highSeason', 'midSeason', 'flatSeason']
eirs = [1, 3, 5, 10, 20, 30, 50, 100, 200, 400]
cms = [0, 0.35, 0.5]
simulation_duration = 365*65

all_monthly_EIRs = pd.DataFrame({'month': np.arange(1, 13, 1)})
cm_filepath = [None]*(len(eirs)*len(seasonalities)*len(cms))
site = [None]*(len(eirs)*len(seasonalities)*len(cms))
index = 0
for ss in range(len(seasonality_names)):
    for ee in eirs:
        for cm in cms:
            # site name
            site[index] ='%s_%i_%i' % (seasonality_names[ss], ee, cm*100)
            # CM filename
            if cm > 0:
                cm_filepath[index] = 'case_management/constant_%i.csv' % round(cm*100)
            else:
                cm_filepath[index] = None
            # create EIR input file with EIRs for all sites and seasonalities
            all_monthly_EIRs[site[index]] = [ee * xx for xx in seasonalities[ss]]
            index = index+1
sweep_coord_df = pd.DataFrame({
    'site': site,
    'simulation_duration': [simulation_duration] * len(site),
    'enable_vital_dynamics': [0] * len(site),
    'CM_filepath': cm_filepath,
    'include_AnnualMalariaSummaryReport': ['summary_report_age_bins/age_bins_set4.csv'] * len(site),
    'include_MonthlyMalariaSummaryReport': [False] * len(site),
    'monthly_summary_report_age_bins_filepath': [None] * len(site),
    'include_MalariaPatientReport': [False] * len(site),
    'NMF_filepath': ['nonmalarial_fevers/nmf_rates_generic.csv'] * len(site),
    'EIR_filepath': ['monthly_eirs/eir_by_sweep_site.csv'] * len(site)
})
sweep_coord_df.to_csv(os.path.join(manifest.input_files_path, 'sweep_sim_coordinator.csv'))
all_monthly_EIRs.to_csv(os.path.join(manifest.input_files_path, 'monthly_eirs/eir_by_sweep_site.csv'))


# create CM input files - assume U5 and adult coverages are the same and that severe coverage is 2x U5 coverage with max of 0.8
for cm in cms:
    if cm > 0:
        filename = 'case_management/constant_%i.csv' % round(cm*100)
        df = pd.DataFrame(data={'simday': [0], 'duration': [-1], 'U5_coverage': [cm], 'adult_coverage': [cm], 'severe_coverage': [min((cm*2), 0.8)]})
        df.to_csv(os.path.join(manifest.input_files_path, filename))



