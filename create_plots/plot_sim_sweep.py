import os
import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
import matplotlib
from plotnine import ggplot, aes, geom_line, facets

font = {'family': 'normal', 'size': 16}
matplotlib.rc('font', **font)

num_seeds = 2
base_filepath_sim_output = '/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/simulation_output/ewan_age_incidence'
filepath_sim = os.path.join(base_filepath_sim_output, 'sweepSeasonEIRCM', 'summary_data_final.csv')
df = pd.read_csv(filepath_sim)
df["season"] = [x.split("_")[0] for x in df.Site]
df["EIR"] = [(x.split("_")[1]) for x in df.Site]
df["CM"] = [(x.split("_")[2]) for x in df.Site]


gg = (
        ggplot(df) +
        aes(x='Age', y='Incidence', color='EIR') +  # , linetype="season"
        geom_line() +
        facets.facet_grid("season ~ CM")
        # + facets.facet_wrap(facets="season", nrow=1)
)

gg.save(os.path.join(base_filepath_sim_output, '_plots', 'sweepSeasonEIRCM.png'))
# for ss in sites:
#
#     reference = ref_df.INC[ref_df.REF_GROUP == int(ss)].tolist()
#     reference = [xx / 1000 for xx in reference]
#     site_name = ref_df.ACD_LOCATION[ref_df.REF_GROUP == int(ss)].tolist()[0]
#     # reference = [3.2, 5, 6.1, 4.75, 3.1, 2.75, 2.7, 1.9, 0.12, 0.8, 0.5, 0.25, 0.1, 0.2, 0.4,
#     #              0.3, 0.2, 0.2, 0.2, 0.15, 0.15, 0.15]
#     # sample_size = [55, 60, 55, 50, 50, 38, 38, 38, 38, 38, 26, 26, 26, 26, 26, 110, 75, 75, 150, 70, 70, 90]
#     # error = np.array(reference)/np.sqrt(np.array(sample_size))
#     ages = np.array(df['Age'].unique())
#
#     for channel in ['Incidence']:  #'Prevalence',
#         plt.figure()
#         plt.plot(df['Age'], df[channel], color='xkcd:sky blue', linewidth=2, label='Simulation')
#         plt.fill_between(df['Age'], df[channel] - df[channel+'_std']/np.sqrt(num_seeds),
#                          df[channel] + df[channel+'_std']/np.sqrt(num_seeds),
#                          color='xkcd:sky blue', alpha=0.5)
#
#         # plot reference
#         plt.plot(ages, reference, marker='o', markersize=18, color='r', label='Reference')
#         # plt.errorbar(ages, reference, yerr=error, marker='o', markersize=8, color='r', label='Reference')
#         plt.legend()
#         plt.title('%s (#%s)' % (site_name, ss))
#         plt.ylabel('Annual  Incidence')
#         plt.xlabel('Age')
#         plt.savefig(os.path.join(base_filepath_sim_output, '_plots', 'site_%s' % ss))
#         plt.show()
