from plotnine import *
import numpy as np
import pandas as pd
from scipy.stats import beta
from pandas.api.types import CategoricalDtype


def get_age_bin_averages(sim_df):
    """
    Get average parasite densities in each age bin, weighting all ages in bin equally (e.g., not weighted by 
    population size)
    Args:
        sim_df (): 

    Returns: age_agg_sim_df

    """
    # age_bins = sim_df['agebin'].unique()
    # remove rows where there are zero people of the measured age bin in the simulation
    sim_df = sim_df[sim_df['Pop'] > 0]
    # get average across all years in age bins and across simulation run seeds
    age_agg_sim_df = sim_df.group_by(['month', 'agebin', 'densitybin', 'Site']).agg(
        asexual_par_dens_freq=('asexual_par_dens_freq', np.mean),
        gametocyte_dens_freq=('gametocyte_dens_freq', np.mean),
        Pop=('Pop', np.mean)
    )
    return age_agg_sim_df


def plot_par_dens_ref_sim_comparison(age_agg_sim_df, ref_df):
    """
    Plot parasite density comparisons with reference
    Stacked barplots of parasite density bins by age
    Args:
        age_agg_sim_df ():
        ref_df ():

    Returns:

    """
    months_of_year = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    # subset simulation output to months in reference dataset
    months = sorted(ref_df['month'].unique())
    cur_df = age_agg_sim_df[age_agg_sim_df['month'].isin(months)]

    # if the maximum reference density bin is < (maximum simulation density bin / 1000),
    # aggregate all simulation densities >= max ref bin into the max ref bin
    #   the final bin will be all densities equal to or above that value
    max_ref_dens = ref_df['densitybin'].dropna().max()
    max_cur_dens = cur_df['densitybin'].dropna().max()
    if max_ref_dens < (max_cur_dens / 1000):
        # get sum of frequencies within higher bins
        all_higher_dens = cur_df[cur_df['densitybin'] >= max_ref_dens]
        sim_agg_higher_dens = all_higher_dens.group_by(['month', 'agebin', 'Site']).agg(
            densitybin=('densitybin', np.min),
            asexual_par_dens_freq=('asexual_par_dens_freq', np.sum),
            gametocyte_dens_freq=('gametocyte_dens_freq', np.sum),
            Pop=('Pop', np.mean))
        # remove higher density bins from df
        cur_df_lower = cur_df[cur_df['densitybin'] < max_ref_dens]
        # add back in the aggregated frequencies across higher density bins
        cur_df = pd.merge(cur_df_lower, sim_agg_higher_dens, how="outer")

    # add zeros for unobserved reference densities up to max_ref_dens
    all_zeros_df = cur_df[['month', 'agebin', 'densitybin', 'Site']]
    ref_df = pd.merge(ref_df, all_zeros_df, how="outer")
    ref_df.fillna(0, inplace=True)

    # combine reference and simulation dataframes
    cur_df['source'] = 'simulation'
    ref_df['source'] = 'reference'
    combined_df0 = pd.concat([cur_df, ref_df], join='outer')

    # = = = = = = = = = #
    # stacked barplots
    # = = = = = = = = = #
    # change type to factors for barplot groupings
    combined_df = combined_df0
    convert_dict = {'densitybin': 'category',
                    'agebin': 'category'}
    combined_df = combined_df.astype(convert_dict)

    # colors
    len_density_bin = len(combined_df['densitybin'].unique())
    num_colors = len_density_bin + 1 if len_density_bin % 2 == 0 else len_density_bin
    #colors = brewer.pal(n=num_colors, name='BrBG')
    #names(colors) = sorted(combined_df['densitybin'].unique())
    # plot
    gg1 = ( ggplot(combined_df, aes(x='agebin', y='asexual_par_dens_freq', fill='densitybin'))
        + geom_bar(position="stack", stat="identity")
        + scale_fill_brewer(palette="BrBG")
        # scale_fill_manual(values=colors, limits=names(colors)) +
        + facet_grid('month~source')
      )

    # = = = = = = = = = #
    # grid of line plots
    # = = = = = = = = = #

    # calculate reference error bounds using Jerrerys interval
    ci_width = 0.95
    alpha = 1 - ci_width
    combined_df0['min_asex'] = np.nan
    combined_df0['max_asex'] = np.nan
    combined_df0['min_gamet'] = np.nan
    combined_df0['max_gamet'] = np.nan
    for rr in range(len(combined_df0.index)):

        if combined_df0['source'].iloc[rr] == 'reference':
            if ((combined_df0['count_asex'].iloc[rr] > 0) &
                    (combined_df0['count_asex'].iloc[rr] < combined_df0['bin_total_asex'].iloc[rr])):
                combined_df0['min_asex'].ilo[rr] = beta.ppf(
                    p=alpha / 2,
                    a=combined_df0['count_asex'].iloc[rr]+0.5,
                    b=combined_df0['bin_total_asex'].iloc[rr] - combined_df0['count_asex'][rr] + 0.5)

                combined_df0['max_asex'].iloc[rr] = beta.ppf(
                    p=1-alpha / 2,
                    a=combined_df0['count_asex'].iloc[rr]+0.5,
                    b=combined_df0['bin_total_asex'].iloc[rr] - combined_df0['count_asex'].iloc[rr] + 0.5)

            if ((combined_df0['count_gamet'].iloc[rr] > 0) &
                (combined_df0['count_gamet'].iloc[rr] < combined_df0['bin_total_gamet'].iloc[rr])):
                combined_df0['min_gamet'].iloc[rr] = beta.ppf(
                    p=alpha / 2,
                    a=combined_df0['count_gamet'].iloc[rr]+0.5,
                    b=combined_df0['bin_total_gamet'].iloc[rr] - combined_df0['count_gamet'].iloc[rr] + 0.5)

                combined_df0['max_gamet'].iloc[rr] = beta.ppf(
                    p=1-alpha / 2,
                    a=combined_df0['count_gamet'].iloc[rr]+0.5,
                    b=combined_df0['bin_total_gamet'].iloc[rr] - combined_df0['count_gamet'].iloc[rr] + 0.5)

    # change facet values to intuitive labels
    combined_df0['month'] = months_of_year[combined_df0['month']]
    month_cat = CategoricalDtype(categories=months_of_year, ordered=True)
    combined_df0['month'] = combined_df0['month'].astype(month_cat)
    all_age_bins = sorted(combined_df0['agebin'].unique())
    age_bin_labels = ['<=' + all_age_bins[1] + " years"]
    for aa in range(len(all_age_bins) - 1):
        age_bin_labels.append(all_age_bins[aa] + '-' + all_age_bins[aa+1] + ' years')

    combined_df0['agebin_index'] = combined_df0['agebin'].isin(all_age_bins)
    combined_df0['agebin'] = age_bin_labels[combined_df0['agebin'].isin(all_age_bins)]
    age_bin_labels_cat = CategoricalDtype(categories=age_bin_labels, ordered=True)
    combined_df0['agebin'] = combined_df0['agebin'].astype(age_bin_labels_cat)

    # plot asexual densities
    gg2 = ( ggplot(combined_df0, aes(x="densitybin", y='asexual_par_dens_freq', color='source'), alpha=0.8)
    + geom_line(size=2)
    + geom_point()
    + scale_x_continuous(trans='log10')
    + geom_errorbar(aes(ymin='min_asex', ymax='max_asex'), width=0.2)
    + theme_bw()
    + ylab('fraction of population')
    + xlab('asexual parasite density bin')
    + scale_color_manual(values={"reference": 'red',
                                 "simulation": 'blue'})
    + facet_grid('agebin~month')
    # scale_fill_brewer(palette = "BrBG") +
    # scale_fill_manual(values=colors, limits=names(colors)) +
    )

    # plot gametocyte densities
    gg3 = ( ggplot(combined_df0, aes(x='densitybin', y='gametocyte_dens_freq', color='source'))
            + geom_line(size=2)
            + geom_point()
            + scale_x_continuous(trans='log10')
            + geom_errorbar(aes(ymin='min_gamet', ymax='max_gamet'), width=0.2)
            + theme_bw()
            + ylab('fraction of population')
            + xlab('gametocyte density bin')
            + scale_color_manual(values={"reference": 'red',
                                 "simulation": 'blue'})
            + facet_grid('agebin~month')
    )

    return list(gg1, gg2, gg3)
