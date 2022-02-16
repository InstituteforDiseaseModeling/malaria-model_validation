import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

font = {'family': 'normal', 'size': 16}
matplotlib.rc('font', **font)

sites = ['1', '2', '3', '4', '5', '7', '8', '9']   # '6',
num_seeds = 2
base_filepath_sim_output = '/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/simulation_output/ewan_age_incidence'
filepath_ref = "/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/reference_data/Cameron_age_incidence_eir.csv"
ref_df = pd.read_csv(filepath_ref)

for ss in sites:
    filepath_sim = os.path.join(base_filepath_sim_output, 'site_%s' % ss, 'summary_data_final.csv')
    df = pd.read_csv(filepath_sim)

    reference = ref_df.INC[ref_df.REF_GROUP == int(ss)].tolist()
    reference = [xx / 1000 for xx in reference]
    site_name = ref_df.ACD_LOCATION[ref_df.REF_GROUP == int(ss)].tolist()[0]
    # reference = [3.2, 5, 6.1, 4.75, 3.1, 2.75, 2.7, 1.9, 0.12, 0.8, 0.5, 0.25, 0.1, 0.2, 0.4,
    #              0.3, 0.2, 0.2, 0.2, 0.15, 0.15, 0.15]
    # sample_size = [55, 60, 55, 50, 50, 38, 38, 38, 38, 38, 26, 26, 26, 26, 26, 110, 75, 75, 150, 70, 70, 90]
    # error = np.array(reference)/np.sqrt(np.array(sample_size))
    ages = np.array(df['Age'].unique())

    for channel in ['Incidence']:  #'Prevalence',
        plt.figure()
        plt.plot(df['Age'], df[channel], color='xkcd:sky blue', linewidth=2, label='Simulation')
        plt.fill_between(df['Age'], df[channel] - df[channel+'_std']/np.sqrt(num_seeds),
                         df[channel] + df[channel+'_std']/np.sqrt(num_seeds),
                         color='xkcd:sky blue', alpha=0.5)

        # plot reference
        plt.plot(ages, reference, marker='o', markersize=18, color='r', label='Reference')
        # plt.errorbar(ages, reference, yerr=error, marker='o', markersize=8, color='r', label='Reference')
        plt.legend()
        plt.title('%s (#%s)' % (site_name, ss))
        plt.ylabel('Annual  Incidence')
        plt.xlabel('Age')
        plt.savefig(os.path.join(base_filepath_sim_output, '_plots', 'site_%s' % ss))
        plt.show()
