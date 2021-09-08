import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

font = {'family': 'normal', 'size': 16}

matplotlib.rc('font', **font)

filepath = '/Users/pselvaraj/Github/malaria-model_validation/ewan_sim_match/sim_data/summary_data_final.csv'
df = pd.read_csv(filepath)
num_seeds = 50

reference = [3.2, 5, 6.1, 4.75, 3.1, 2.75, 2.7, 1.9, 0.12, 0.8, 0.5, 0.25, 0.1, 0.2, 0.4,
             0.3, 0.2, 0.2, 0.2, 0.15, 0.15, 0.15]
sample_size = [55, 60, 55, 50, 50, 38, 38, 38, 38, 38, 26, 26, 26, 26, 26, 110, 75, 75, 150, 70, 70, 90]
error = np.array(reference)/np.sqrt(np.array(sample_size))
ages = np.array(df['Age'].unique())

for channel in ['Prevalence', 'Incidence']:
    plt.figure()
    plt.plot(df['Age'], df[channel], color='xkcd:sky blue', linewidth=2, label='Simulation')
    plt.fill_between(df['Age'], df[channel] - df[channel+'_std']/np.sqrt(num_seeds),
                     df[channel] + df[channel+'_std']/np.sqrt(num_seeds),
                     color='xkcd:sky blue', alpha=0.5)

    # plot reference
    # plt.plot(ages, reference, marker='o', markersize=18, color='r')
    plt.errorbar(ages, reference, yerr=error, marker='o', markersize=8, color='r', label='Reference')
    plt.legend()
    plt.title('Dielmo')
    plt.ylabel('Annual  Incidence')
    plt.xlabel('Age')
    plt.savefig('/Users/pselvaraj/Github/malaria-model_validation/ewan_sim_match/figures/Dielmo.png')
    plt.show()