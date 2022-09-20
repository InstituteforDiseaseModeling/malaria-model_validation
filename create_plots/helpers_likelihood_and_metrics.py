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

import os
import pandas as pd
from scipy import stats
import warnings
import numpy as np
from plotnine import ggplot, aes, geom_point, xlab, ylab, coord_fixed, geom_abline, theme_classic, themes, \
    ggtitle, geom_smooth
import collections


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
    if any(ref_rows_na):
        sim_rows_0 = sum(combined_df[sim_column] < 0.0001)
        if ref_rows_na.equals(sum(sim_rows_0)):
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
                import numpy as np
                from scipy.special import xlogy, gammaln
                import math

                x = np.array(list(cur_df['ref_bin_count']))
                n = sum(cur_df['ref_bin_count'])
                p = np.array(list(cur_df[sim_column]))

                # stats.multinomial.logpmf(x=x, n=n, p =p)
                # state.multinomial sometimes return nan value, replace it with the following line
                try:
                    result = math.log(np.exp(gammaln(n + 1) + np.sum(xlogy(x, p) - gammaln(x + 1), axis=-1)))
                    loglikelihood += result
                except Exception as ex:
                    print(x, n, p)
                    print(ex)
                    loglikelihood = -np.inf
            else:
                warnings.warn(f'Either the sum of individuals across bins in the reference dataset does not match the '
                              f'reported total number of individuals included or different density bins are used in '
                              f'the reference and simulation. This site-age is being skipped: {ss} - {aa}')

        loglik_df = pd.concat([loglik_df, pd.DataFrame({'site_month': [ss], 'loglikelihood': [loglikelihood]})])

    return loglik_df
# endregion


# region other quantitative comparison metrics
def calc_mean_rel_diff(combined_df, sim_colname='simulation'):
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
    combined_df['rel_diff'] = abs((combined_df['reference'] - combined_df[sim_colname]) / combined_df['reference'])
    combined_df['abs_diff'] = abs(combined_df['reference'] - combined_df[sim_colname])
    # todo: need code review
    # mean_diff_df = combined_df % > % group_by(Site) % > %
    # summarise(mean_rel_diff=mean(rel_diff),
    #           mean_abs_diff=mean(abs_diff))
    mean_diff_df = combined_df.groupby(['Site']).agg(
        mean_rel_diff=('rel_diff', np.nanmean),
        mean_abs_diff=('abs_diff', np.nanmean)
    ).reset_index()

    mean_all_sites = pd.DataFrame({'Site': ['all_sites'],
                                   'mean_rel_diff': [np.nanmean(combined_df['rel_diff'])],
                                   'mean_abs_diff': [np.nanmean(combined_df['abs_diff'])]})
    mean_diff_df = pd.concat([mean_diff_df, mean_all_sites])

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
    mean_slope_diff_df = combined_df.groupby(['Site']).agg(
        mean_rel_slope_diff=('rel_diff', np.nanmean),
        mean_abs_slpe_diff=('abs_diff', np.nanmean)
    ).reset_index()

    return mean_slope_diff_df


def corr_ref_sim_points(combined_df):
    """
    Calculate the correlation between reference and matched simulation data points.
    Args:
        combined_df (): A dataframe containing both the reference and matched simulation output

    Returns: A list with two elements:
            1) A gg scatterplot showing the reference values against the simulation values
            2) A dataframe summarizing the linear regression results

    """

    metric = combined_df['metric'].iloc[0]
    if 'site_month' in combined_df.columns:
        combined_df['Site'] = combined_df['site_month']
    min_value = min(combined_df['reference'].min(skipna=True), combined_df['simulation'].min(skipna=True))
    max_value = max(combined_df['reference'].max(skipna=True), combined_df['simulation'].max(skipna=True))
    gg = (ggplot(combined_df, aes(x='reference', y='simulation', color='Site', fill='Site'))
          + geom_point(size=2)
          + ylab('simulation ' + str(metric))
          + xlab('reference ' + str(metric))
          + ggtitle('Ref versus sim values')
          + coord_fixed(ratio=1, xlim=(min_value, max_value), ylim=(min_value, max_value))
          + geom_abline(slope=1, intercept=0, color='grey', alpha=0.5)
          + geom_smooth(method="lm", fill=None, se=False, alpha=0.5, size=0.5)
          + theme_classic()
          + themes.theme(plot_title=themes.element_text(size=12)))

    # create data frame with information about linear regression and correlation for each Site
    lm_summary = collections.defaultdict(list)
    groups = combined_df.groupby("Site")
    for name, group in groups:
        # remove rows without both simulation and reference values
        group = group.dropna(axis=0, how='any', subset=['simulation', 'reference'])
        slope, intercept, r_value, p_value, std_err = stats.linregress(group['reference'], group['simulation'])
        lm_summary['Site'].append(name)
        lm_summary['slope'].append(slope)
        lm_summary['r.squared'].append(r_value**2)
        lm_summary['p.value'].append(p_value)
        lm_summary['nobs'].append(len(group))
    lm_summary = pd.DataFrame(lm_summary)
    return gg, lm_summary


def corr_ref_deriv_sim_points(combined_df):
    """
    Calculate the correlation between reference and matched simulation slopes (derivatives) when moving from the
    youngest to the oldest age group.
    Args:
        combined_df (): A dataframe containing both the reference and matched simulation output

    Returns: A list with two elements:
            1) A gg scatterplot showing the reference slopes against the simulation slopes
            2) A dataframe summarizing the linear regression results
            3) The combined_df dataframe, with columns added giving the simulation and reference slopes

    """

    metric = combined_df['metric'].iloc[0]
    # calculate the slope when moving between age groups
    combined_df['sim_slope_to_next'] = np.nan
    combined_df['ref_slope_to_next'] = np.nan
    sites = combined_df['Site'].unique()
    for ss in sites:
        cur_df = combined_df[combined_df['Site'] == ss]
        cur_ages = sorted(cur_df['mean_age'].unique())
        for aa in range(len(cur_ages) - 1):
            # todo: need code review
            combined_df_row = (combined_df['Site'] == ss) & (combined_df['mean_age'] == cur_ages[aa])
            # get the slope between the value for this age and the next-largest age group
            sim_val_cur = cur_df[cur_df['mean_age'] == cur_ages[aa]]['simulation'].iloc[0]
            sim_val_next = cur_df[cur_df['mean_age'] == cur_ages[aa + 1]]['simulation'].iloc[0]
            sim_slope = (sim_val_next - sim_val_cur) / (cur_ages[aa + 1] - cur_ages[aa])
            combined_df['sim_slope_to_next'][combined_df_row] = sim_slope

            ref_val_cur = cur_df[cur_df['mean_age'] == cur_ages[aa]]['reference'].iloc[0]
            ref_val_next = cur_df[cur_df['mean_age'] == cur_ages[aa + 1]]['reference'].iloc[0]
            ref_slope = (ref_val_next - ref_val_cur) / (cur_ages[aa + 1] - cur_ages[aa])
            combined_df['ref_slope_to_next'][combined_df_row] = ref_slope

    min_value = min(combined_df['ref_slope_to_next'].min(skipna=True), combined_df['sim_slope_to_next'].min(skipna=True))
    max_value = max(combined_df['ref_slope_to_next'].max(skipna=True), combined_df['sim_slope_to_next'].max(skipna=True))

    gg = (ggplot(combined_df, aes(x='ref_slope_to_next', y='sim_slope_to_next', color='Site', fill='Site'))
          + geom_point(size=2)
          + ylab('simulation slopes for age-' + str(metric))
          + xlab('reference slopes for age-' + str(metric))
          + ggtitle('Ref versus sim slopes')
          + coord_fixed(ratio=1, xlim=(min_value, max_value), ylim=(min_value, max_value))
          + geom_abline(slope=1, intercept=0, color='grey', alpha=0.5)
          # + geom_smooth(method="lm", fill=None, se=False, alpha=0.5, size=0.5)
          + theme_classic()
          + themes.theme(plot_title=themes.element_text(size=12)))

    # create data frame with information about linear regression and correlation for each Site
    lm_summary = collections.defaultdict(list)
    groups = combined_df.groupby("Site")
    for name, group in groups:
        # remove rows without both simulation and reference values
        group = group.dropna(axis=0, how='any', subset=['sim_slope_to_next', 'ref_slope_to_next'])
        slope, intercept, r_value, p_value, std_err = stats.linregress(group['ref_slope_to_next'], group['sim_slope_to_next'])
        lm_summary['Site'].append(name)
        lm_summary['slope'].append(slope)
        lm_summary['r.squared'].append(r_value ** 2)
        lm_summary['p.value'].append(p_value)
        lm_summary['nobs'].append(len(group))
    lm_summary = pd.DataFrame(lm_summary)
    # lm_summary = lm_summary[lm_summary$term != '(Intercept)',]
    # colnames(lm_summary)[colnames(lm_summary) == 'estimate'] = 'slope'
    return gg, lm_summary, combined_df
# endregion


# region summarize new simulation success against benchmark simulation
def add_to_summary_table(combined_df, plot_output_filepath, validation_relationship_name,
                         rel_change_threshold=0.1):

    mean_diff_df_new = calc_mean_rel_diff(combined_df, sim_colname='simulation')
    mean_diff_df_bench = calc_mean_rel_diff(combined_df, sim_colname='benchmark')
    mean_diff_df_bench.rename(columns={'mean_rel_diff': 'mean_rel_diff_bench',
                                       'mean_abs_diff': 'mean_abs_diff_bench'},
                              inplace=True)
    mean_diff_df = pd.merge(mean_diff_df_new, mean_diff_df_bench, on='Site',how='outer')

    # determine which sites improved, got worse, or stayed close to the same for absolute difference
    mean_diff_df['change_abs_diff'] = mean_diff_df['mean_abs_diff_bench'] - mean_diff_df['mean_abs_diff']
    mean_diff_df = mean_diff_df[mean_diff_df['change_abs_diff'].notnull()]
    #todo: need code review in the following line:
    # R code: mean_diff_df$abs_diff_changed = abs(mean_diff_df$change_abs_diff)/mean_diff_df$mean_abs_diff_bench > rel_change_threshold
    mean_diff_df['abs_diff_changed'] = abs(mean_diff_df['change_abs_diff']) / mean_diff_df['mean_abs_diff_bench'] > rel_change_threshold
    mean_diff_df['change_type'] = 'better'
    mean_diff_df.loc[mean_diff_df['change_abs_diff'] < 0, 'change_type'] = 'worse'
    mean_diff_df.loc[~mean_diff_df['abs_diff_changed'], 'change_type'] = 'similar'

    # save results as row in dataframe
    summary_df = pd.DataFrame({'validation_relationship': validation_relationship_name,
                               'abs_diff_new': mean_diff_df[mean_diff_df['Site'] == 'all_sites']['mean_abs_diff'],
                               'abs_diff_bench': mean_diff_df[mean_diff_df['Site'] == 'all_sites']['mean_abs_diff_bench'],
                               # 'ave_rel_diff_new_sim'=mean_diff_df$mean_rel_diff[mean_diff_df$Site == 'all_sites'],
                               # 'ave_rel_diff_bench_sim'=mean_diff_df$mean_rel_diff_bench[mean_diff_df$Site == 'all_sites'],
                               'num_sites_better': len(mean_diff_df[(mean_diff_df['Site'] != 'all_sites') & (mean_diff_df['change_type'] == 'better')]),
                               'num_sites_similar': len(mean_diff_df[(mean_diff_df['Site'] != 'all_sites') & (mean_diff_df['change_type'] == 'similar')]),
                               'num_sites_worse': len(mean_diff_df[(mean_diff_df['Site'] != 'all_sites') & (mean_diff_df['change_type'] == 'worse')])})

    # write to csv, adding as new row if csv already exists but relationship doesn't
    summary_filepath = os.path.join(plot_output_filepath, 'summary_table_sim_benchmark.csv')
    if os.path.exists(summary_filepath):
        # add rows to existing csv
        existing_summaries = pd.read_csv(summary_filepath)
        if validation_relationship_name in existing_summaries['validation_relationship']:
            existing_summaries[existing_summaries['validation_relationship'] == validation_relationship_name] = \
                summary_df.iloc[0]
            summary_df = existing_summaries
        else:
            summary_df = pd.concat([existing_summaries, summary_df])

    summary_df.to_csv(summary_filepath, header=True, index=False)
# endregion
