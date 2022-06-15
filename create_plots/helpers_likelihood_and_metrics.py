# helpers_loglikelihood_and_metrics.py
#
#  These functions calculate likelihoods or other comparison metrics evaluating how well simulation and reference data
#       agree.
#  The main functions return a data frame where each row corresponds to a site (or site-month) and columns contain
#       the quantitative measure.
#  Some functions also return a ggplot corresponding to the quantitative comparison.
#
#  The loglikelihood evaluations are approximate and generally assume that the mean simulated value is the 'true'
#       population value and ask how likely it was to observe the reference dataset (given the study sample size).

import pandas as pd
from scipy import stats
import warnings
import numpy as np


# region: loglikelihood functions for each validation relationship
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# prevalence
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
def get_prev_loglikelihood(combined_df, sim_column='simulation'):
    """
    Calculate an approximate likelihood for the simulation parameters for each site. This is estimated as the product,
    across age groups, of the probability of observing the reference values if the simulation means represented the
    true population mean
    Args:
        combined_df (): A dataframe containing both the reference and matched simulation output
        sim_column (): The name of the column of combined_df to use as the simulation output
    Returns: A dataframe of loglikelihoods where each row corresponds to a site-month

    """

    combined_df['prob_pos_sim'] = combined_df[sim_column]
    # only include reference sites where the sample sizes were reported
    combined_df = combined_df.dropna(subset=['total_sampled'])
    sites = combined_df['site_month'].unique()

    # likelihood approximated for each age group as probability of obtaining observed num_pos from total_sampled if simulation is true prevalence
    loglik_by_site = [None] * len(sites)
    for ss in range(len(sites)):
        cur_df = combined_df[combined_df['site_month'] == sites[ss]]
        # get product of probabilities for each age group. If any age groups don't match, entire site is Nan.
        # todo: need code review for this if condition
        # if (!any( is.na(cur_df$prob_pos_sim))){
        if not cur_df['prob_pos_sim'].isna().values.any():
            loglik_total=0
            for rr in range(len(cur_df)): # iterate through age groups (each row corresponds to a different age group)
                # todo: need code review for this math
                # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.binom.html
                # loglik_total = loglik_total + math.log(dbinom(x=cur_df$num_pos[rr], size = cur_df$total_sampled[rr],
                #                                           prob = cur_df$prob_pos_sim[rr]))
                loglik_total = loglik_total + stats.binom.logpmf(k=cur_df['num_pos'].iloc[rr],
                                                                 n=cur_df['total_sampled'].iloc[rr],
                                                                 p=cur_df['prob_pos_sim'].iloc[rr])

            loglik_by_site[ss] = loglik_total

    loglik_df = pd.DataFrame({'site_month': sites, 'loglikelihood': loglik_by_site})

    return loglik_df


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# parasite density
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
def get_dens_loglikelihood(combined_df, sim_column='simulation'):
    """
    Calculate an approximate likelihood for the simulation parameters for each site. This is estimated as the product,
    across age groups, of the probability of observing the reference values if the simulation means represented the
    true population mean
    Args:
        combined_df (): A dataframe containing both the reference and matched simulation output
        sim_column (): The name of the column of combined_df to use as the simulation output

    Returns: A dataframe of loglikelihoods where each row corresponds to a site-month

    """

    # remove rows where there is no reference data (sometimes there are different density bins in different sites, so
    # rows with NA are created for the 'missing' bins - but check that there aren't values in the simulation either)
    # todo: need code review
    ref_rows_na = combined_df['ref_bin_count'].isna()
    if len(ref_rows_na):
        sim_rows_0 = combined_df[sim_column] < 0.0001
        if ref_rows_na.equals(sim_rows_0):
            combined_df = combined_df[~ref_rows_na]
        else:
            warnings.warn("There may be a mismatch in the age bins from the reference data and simulation data for at "
                          "least one site. No rows were removed.")

    # remove rows without simulation output
    combined_df = combined_df[~combined_df[sim_column].isna()]

    # assuming reference data comes from a multinomial draw: likelihood of getting the observed distribution of
    # parasite densities assuming the simulations show the true population-level frequencies
    # iterate through groups of month-age-site
    loglik_df = pd.DataFrame(columns =['site_month', 'loglikelihood'])
    site_months = combined_df['site_month'].unique()
    for ss in site_months:
        loglikelihood = 0
        cur_agebins = combined_df[combined_df['site_month'] == ss]['agebin'].unique()
        for aa in cur_agebins:
            cur_df = combined_df[(combined_df['site_month'] == ss) & (combined_df['agebin'] == aa)]
            # check that the sum of counts matches the sum column
            if sum(cur_df['ref_bin_count']) == cur_df['ref_total'].iloc[0] and len(cur_df['ref_total'].unique()) == 1:
                # todo: need code review for this math
                # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.multinomial.html
                # in Python:
                # stats.multinomial.logpmf([3, 4], n=7, p=[.3, .7])
                # Out[29]: array(-1.48327013)
                # in R:
                # log(dmultinom(c(3, 4), prob=c(0.3, 0.7)))
                # [1] - 1.48327
                # loglikelihood = loglikelihood + log(dmultinom(x=cur_df$ref_bin_count, prob=cur_df[[sim_column]]))
                loglikelihood = loglikelihood + stats.multinomial.logpmf(x=cur_df['ref_bin_count'],
                                                                         n = sum(cur_df['ref_bin_count']),
                                                                         p=cur_df[sim_column])
            else:
                warnings.warn(f'Either the sum of individuals across bins in the reference dataset does not match the '
                              f'reported total number of individuals included or different density bins are used in '
                              f'the reference and simulation. This site-age is being skipped: {ss} - {aa}')

        loglik_df = pd.concat([loglik_df, pd.DataFrame({'site_month': ss, 'loglikelihood': loglikelihood})])

    return loglik_df
# endregion


# region other quantitative comparison metrics
def calc_mean_rel_diff(combined_df):
    """
    Calculate the mean relative difference between reference and matched simulation values across ages.
    For a given age, calculate (reference - simulation)/reference. Report the mean across all ages.
    Args:
        combined_df (): A dataframe containing both the reference and matched simulation output

    Returns: A dataframe with the mean absolute difference and the mean magnitude of the relative difference between
             all pairs of reference and simulation values

    """
    if 'site_month' in combined_df.columns:
        combined_df['Site'] = combined_df['site_month']
    combined_df['rel_diff'] = abs((combined_df['reference'] - combined_df['simulation']) / combined_df['reference'])
    combined_df['abs_diff'] = abs(combined_df['reference'] - combined_df['simulation'])
    # todo: need code review
    # mean_diff_df = combined_df % > % group_by(Site) % > %
    # summarise(mean_rel_diff=mean(rel_diff),
    #           mean_abs_diff=mean(abs_diff))
    mean_diff_df = combined_df.group_by(['Site']).agg(
        mean_rel_diff=('rel_diff', np.mean),
        mean_abs_diff=('abs_diff', np.mean)
    )

    return mean_diff_df


def calc_mean_rel_slope_diff(combined_df):
    """
     Calculate the mean relative difference between the slopes in the reference and matched simulation dataset, moving
     from the youngest to the oldest age bin.
     For a given age, calculate (reference slope - simulation slope)/reference slope. Report the mean across all ages.
    Args:
        combined_df (): A dataframe containing both the reference and matched simulation output

    Returns: A dataframe with the mean absolute difference and the mean magnitude of the relative difference between
             all pairs of reference and simulation values

    """
    if 'site_month' in combined_df.columns:
        combined_df['Site'] = combined_df['site_month']
    combined_df['rel_diff'] = abs((combined_df['ref_slope_to_next'] - combined_df['sim_slope_to_next'])
                                  / combined_df['ref_slope_to_next'])
    combined_df['abs_diff'] = abs(combined_df['ref_slope_to_next'] - combined_df['sim_slope_to_next'])
    # todo: should we replace all numpy.mean with numpy.nanmean() to ignore the nan values
    mean_slope_diff_df = combined_df.group_by(['Site']).agg(
        mean_rel_slope_diff=('rel_diff', np.nanmean),
        mean_abs_slpe_diff=('abs_diff', np.nanmean)
    )

    return mean_slope_diff_df
# endregion
