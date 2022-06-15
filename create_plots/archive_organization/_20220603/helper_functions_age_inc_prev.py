from plotnine import *
import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype


def get_substr(site_name_str, index):
    return site_name_str.split("_")[1][index]


def get_mean_from_upper_age(cur_age, upper_ages):
    if not upper_ages or cur_age not in upper_ages:
        return None
    else:
        mean_ages = [upper_ages[0] / 2]
        for i in range(len(upper_ages) - 1):
            mean_ages.append((upper_ages[i] + upper_ages[i + 1]) / 2)
        return mean_ages[upper_ages.index(cur_age)]


def plot_inc_ref_sim_comparison(sim_df, ref_df):
    """
    Plot age-incidence comparisons with reference
    Args:
        sim_df ():
        ref_df ():

    Returns:

    """
    sim_df['Incidence'] = sim_df['Incidence'] * sim_df['p_detect_case']

    # set up reference and simulation dataset columns
    ref_df['Incidence'] = ref_df['INC'] / 1000
    ref_df['mean_age'] = (ref_df['INC_LAR'] + ref_df['INC_UAR']) / 2
    ref_df = pd.DataFrame(data={'Incidence': ref_df['Incidence'],
                                'mean_age': ref_df['mean_age'],
                                'Site': ref_df['Site'],
                                'Pop_size': ref_df['POP'],
                                'year': ref_df['START_YEAR']})
    sim_df = pd.DataFrame(data={'Incidence': sim_df['Incidence'],
                                'mean_age': sim_df['mean_age'],
                                'Site': sim_df['Site'],
                                'Pop_size': np.nan,
                                'year': np.nan})

    ref_df['source'] = 'reference'
    sim_df['source'] = 'simulation'

    df_combined = pd.concat(sim_df, ref_df)

    source_cat = CategoricalDtype(categories=['reference', 'simulation'], ordered=True)
    df_combined['source'] = df_combined['source'].astype(source_cat)

    gg = (ggplot(df_combined, aes(x='mean_age', y='Incidence', color='source', group='year'))
          + geom_line()
          + geom_point()
          + scale_color_manual(values={"reference": 'red',
                                       "simulation": 'blue'})
          + xlab('age (midpoint of age bin)')
          + ylab('incidence')
          + facet_wrap('Site', ncol=4)
          + theme_bw())
    return gg
