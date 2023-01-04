# helpers_reformat_sim_ref_dfs.R
#
# This script contains functions used to process, align, and format reference and simulation data.
#    There is one main function per validation relationship, as well as a number of shared functions
#    that are called for common data manipulations. The main function for a validation relationship
#    takes the reference data, new simulation outputs, and benchmark simulation outputs (optional)
#    as inputs and returns data frames that are used in downstream plotting and comparisons.

import numpy as np
import warnings
import pandas as pd
import math
import os
import datetime
import random


# region: helper functions
def get_mean_from_upper_age(cur_age, upper_ages):
    """
    Using the upper bounds of a set of age bins, return the mean age of an individuals in a particular bin
    Args:
        cur_age (): The upper age bound of the current age bin
        upper_ages (): The upper age bounds of all age bins

    Returns: The mean age of someone in the current age bin

    """
    if not upper_ages or cur_age not in upper_ages:
        return None
    else:
        mean_ages = [upper_ages[0] / 2]
        for i in range(len(upper_ages) - 1):
            mean_ages.append((upper_ages[i] + upper_ages[i + 1]) / 2)
        return mean_ages[upper_ages.index(cur_age)]


def get_age_bin_averages(sim_df):
    """
    get average fraction of individuals in each age bin that fall into each parasite density bin, weighting all ages in
     bin equally (e.g., not weighted by population size)
    Args:
        sim_df (): A dataframe with simulation output where each row corresponds to an unique combination of {age bin,
        parasite density bin, month, year, site, and run seed}

    Returns: A dataframe with the average gametocyte and asexual parasite density within all groupings

    """

    # remove rows where there are zero people of the measured age bin in the simulation
    sim_df = sim_df[sim_df['Pop'] > 0]
    # get the simulation mean across runs
    if (sim_df['Pop'] == sim_df['Pop'].iloc[1]).all():
        # get average across all years in age bins and across simulation run seeds
        age_agg_sim_df = sim_df.groupby(['month', 'mean_age', 'agebin', 'densitybin', 'Site']).agg(
            asexual_par_dens_freq=('asexual_par_dens_freq', np.nanmean),
            gametocyte_dens_freq=('gametocyte_dens_freq', np.nanmean),
            Pop=('Pop', np.nanmean)).reset_index()
        return age_agg_sim_df
    else:
        warnings.warn("Different population sizes found across years within an age group... "
                      "need to set up weighted averaging of parasite densities.")
        # If this warning is triggered, use the population size to calculate the number of individuals in each density
        # bin, then aggregate the sum, then divide by the total aggregated population
        # todo: NYI
        return None


def match_sim_ref_ages(ref_df, sim_df, bench_df=pd.DataFrame()):
    """
    Check that ages match between reference and simulation. if there is a small difference (<1 year), update simulation to use same ages as reference.
    Args:
        ref_df (): A dataframe with the reference values. Has column mean_age.
        sim_df (): A dataframe with the simulation values. Has column mean_age.
        bench_df (): A dataframe with the benchmark simulation values. If included, must have column mean_age.

    Returns: A list where the first element is the updated (main) simulation dataframe and the second element is the benchmark simulation dataframe

    """

    def intersection(lst1, lst2):
        return list(set(lst1) & set(lst2))
    sites = intersection(sim_df['Site'].unique(), ref_df['Site'].unique())
    for ss in sites:
        ages_ref = sorted(ref_df[ref_df['Site'] == ss]['mean_age'].unique())
        ages_sim = sorted(sim_df[sim_df['Site'] == ss]['mean_age'].unique())

        missing_ages_ref = [x for x in ages_ref if x not in ages_sim]
        missing_ages_sim = [x for x in ages_sim if x not in ages_ref]

        if (len(bench_df) > 0) and ages_sim != sorted(bench_df[bench_df['Site'] == ss]['mean_age'].unique()):
            warnings.warn(f"The age bins used in the benchmarking simulation are different from those used in the new "
                          f"simulation for site: {ss}. This benchmark simulation will be excluded.")
            bench_df = bench_df[bench_df['Site'] != ss]

        if (not all(age_ref <= age_sim + 0.1 for age_ref, age_sim in zip(ages_ref, ages_sim))
                or (not all(age_ref >= age_sim - 0.1 for age_ref, age_sim in zip(ages_ref, ages_sim)))):
            print(f'Imperfect age match between reference and simulations for site: {ss}')
            print('...  Mismatched reference / simulation ages are:')
            for missing_age_ref, missing_age_sim in zip(missing_ages_ref, missing_ages_sim):
                print(f'     {missing_age_ref} / {missing_age_sim},')
            print('... For age thresholds that differ by less than a year, replacing simulation age with reference age.')

            # check whether the missing ages are simply off by <1 year. If so, replace simulation age with nearby reference age
            for mm in missing_ages_ref:
                sim_replace_age = [x for x in missing_ages_sim if math.fabs(x - mm) < 1]
                if len(sim_replace_age) == 1:
                    sim_df.loc[(sim_df['Site'] == ss) & (sim_df['mean_age'] == sim_replace_age[0]), 'mean_age'] = mm
                    if len(bench_df) > 0 and any(bench_df['Site'] == ss):
                        bench_df.loc[(bench_df['Site'] == ss) & (bench_df['mean_age'] == sim_replace_age[0]), 'mean_age'] = mm

            # update sim ages
            ages_sim = sorted(sim_df[sim_df['Site'] == ss]['mean_age'].unique())

            if (not all(age_ref <= age_sim + 0.1 for age_ref, age_sim in zip(ages_ref, ages_sim))
                    or (not all(age_ref >= age_sim - 0.1 for age_ref, age_sim in zip(ages_ref, ages_sim)))):
                print('...After adjustment, there remains an imperfect match between reference and simulation age bins.')
                print(f'      Reference has {len(ages_ref)} age groups and simulation has {len(ages_sim)} age groups.')
            else:
                print('... All age bins now match.')

    return sim_df, bench_df


def get_available_sites_for_relationship(coord_csv, simulation_output_filepath, relationship_name,
                                         relationship_sim_filename):
    """
    Determine which of the simulation sites both are indicated by the coordinator csv to be included in this validation
     relationship and have the relevant simulation output
    Args:
        coord_csv (): A dataframe detailing the sites simulated for each validation relationship and the corresponding
                      reference dataset
        simulation_output_filepath (): The filepath where simulation output is located
        relationship_name (): The column name from coord_csv corresponding to the current validation relationship
                              (values in this column indicate whether or not a site is used for that relationship)
        relationship_sim_filename (): The name of the simulation output file used for this validation relationship

    Returns: A vector of the simulation site names that should be included for this validation relationship

    """
    coord_sites = coord_csv[(~coord_csv['site'].isna()) & (coord_csv[relationship_name] == 1)]['site']
    available_sites = list()
    for ii in range(len(coord_sites)):
        if os.path.isfile(os.path.join(simulation_output_filepath,  coord_sites.iloc[ii], relationship_sim_filename)):
            available_sites.append(coord_sites.iloc[ii])

    return available_sites


def combine_higher_dens_freqs(sim_df_cur, max_ref_dens, max_magnitude_difference=100):
    """
    Aggregate simulation parasite density bins that are substantially above the maximum reference bin.
    Specifically, if the maximum reference density bin is < (maximum simulation density bin / max_magnitude_difference),
    aggregate all simulation densities >= max ref bin into the max ref bin. The new final density bin will be all
    densities equal to or above that value
    Args:
        sim_df_cur (): A dataframe with simulation results in the original density bins
        max_ref_dens (): A number indicating the maximum parasite density bin included in the reference dataset
        max_magnitude_difference (): The maximum allowable difference between the maximum reference and simulation
                                     density bins before the largest simulation bins are aggregated

    Returns: The updated simulation dataframe with either the original or aggregated density bins

    """
    if max_ref_dens < sim_df_cur['densitybin'].max(skipna=True) / max_magnitude_difference:
        # get sum of frequencies within higher bins
        all_higher_dens = sim_df_cur[sim_df_cur['densitybin'] >= max_ref_dens]
        sim_agg_higher_dens = all_higher_dens.group_by(['month', 'agebin', 'Site']).agg(
            densitybin=('densitybin', np.nanmin),
            asexual_par_dens_freq=('asexual_par_dens_freq', np.nansum),
            gametocyte_dens_freq=('gametocyte_dens_freq', np.nansum),
            Pop=('Pop', np.nanmean)
        )
        # remove higher density bins from df
        cur_df_lower = sim_df_cur[sim_df_cur['densitybin'] < max_ref_dens]
        # add back in the aggregated frequencies across higher density bins
        sim_df_cur = pd.merge(cur_df_lower, sim_agg_higher_dens, how='outer')

    return sim_df_cur


def get_fraction_in_infectious_bin(sim_df):
    """
    Translate the infectiousness values from simulation output into the fraction of individuals from each {age bin,
    month, densitybin, run number} group that are in each infectiousness bin.
    Calculate the average value across seeds and years.
    Args:
        sim_df (): A dataframe containing the (unaggregated) simulation output

    Returns: A dataframe including the mean bin frequency values across simulation seeds and within age groups
    (assuming equal population sizes and weighting for all ages in an age bin)

    """
    # note, since this is an average over the reporting period, these may not be whole numbers
    sim_df['infectiousness_bin_count'] = sim_df['infectiousness_bin_freq'] * sim_df['Pop']

    # check that people are never infectious while they have no parasites - send a warning if not
    if len(sim_df[(sim_df['infectiousness_bin'] == 0)
                  & (sim_df['infectiousness_bin'] > 0)
                  & (sim_df['infectiousness_bin_count'] > 0)]) > 0:
        warnings.warn('Some individuals without parasites are reporting that they are infectious to mosquitoes... '
                      'this suggests a bug.')

    # aggregate individuals within an age bin (across years)
    sim_df_agg1 = sim_df.groupby(['Site', 'month', 'Run_Number', 'agebin', 'densitybin', 'infectiousness_bin']).agg(
        infectiousness_bin_count=('infectiousness_bin_count', np.nansum),
        Pop=('Pop', np.nansum)
    ).reset_index()

    # get total within each group (for denominator)
    sim_df_group_total = sim_df_agg1.groupby(['Site', 'month', 'Run_Number', 'agebin', 'densitybin']).agg(
        infectiousness_group_sum=('infectiousness_bin_count', np.nansum),
        Pop=('Pop', np.nansum)
    ).reset_index()
    # calculate proportion of all individuals in a group fell in each infectiousness bin
    sim_df_agg1 = pd.merge(sim_df_agg1, sim_df_group_total, on=['Site', 'month', 'Run_Number', 'agebin', 'densitybin'],
                    how='outer')
    sim_df_agg1['group_infectiousness_freq'] = sim_df_agg1['infectiousness_bin_count'] / \
                                               sim_df_agg1['infectiousness_group_sum']
    # # check that sum within a group is 1
    # sim_subset = sim_df_agg1[sim_df_agg1$month==1 & sim_df_agg1$agebin==5 & sim_df_agg1$Run_Number == 0 & sim_df_agg1$densitybin==500,]
    # sum(sim_subset$group_infectiousness_freq)

    # get simulation average across seeds and within age groups (assuming equal population sizes for all ages in bin)
    sim_df_agg2 = sim_df_agg1.groupby(['Site', 'month', 'agebin', 'densitybin', 'infectiousness_bin']).agg(
        infect_sd=('group_infectiousness_freq', np.nanstd),
        infectiousness_bin_freq=('group_infectiousness_freq', np.nanmean)
    ).reset_index()

    # check that the sum within all groups is 1
    # todo: this variable is defined but not used
    sim_df_check = sim_df_agg2.groupby(['Site', 'month', 'agebin', 'densitybin']).agg(
        infectiousness_dens_bin_sum=('infectiousness_bin_freq', np.nansum)
    ).reset_index()

    if not all(sim_df_check['infectiousness_dens_bin_sum'] - 0.001 < 1) or \
            not all(sim_df_check['infectiousness_dens_bin_sum'] + 0.001 > 1):
        warnings.warn(
            "The sum of infectiousness bin frequencies is not 1 for at least one group. Recommend checking for bugs.")

    return sim_df_agg2

# endregion

# region: main reformatting functions
def prepare_inc_df(coord_csv, simulation_output_filepath, base_reference_filepath, benchmark_simulation_filepath=None):
    """
    Read in, align, and combine reference and simulation data for all sites associated with the incidence-by-age
    validation relationship.
    Args:
        coord_csv (): A dataframe detailing the sites simulated for each validation relationship and the corresponding
                      reference dataset
        simulation_output_filepath (): The filepath where simulation output is located
        base_reference_filepath (): The filepath where reference datasets are located
        benchmark_simulation_filepath (): The filepath where benchmark simulation output is located. If NA,
                                          no comparisons are made against benchmark simulations

    Returns: A dataframe containing the combined reference and simulation data for this validation relationship

    """
    # determine which of the age-incidence sites have the relevant simulation output
    available_sites = get_available_sites_for_relationship(coord_csv, simulation_output_filepath,
                                                           relationship_name='age_incidence',
                                                           relationship_sim_filename='inc_prev_data_final.csv')

    # iterate through sites, aggregating all age-incidence simulation data into one dataframe and reference data into a second dataframe (for the relevant sites)
    #   if a benchmark simulation directory was specified and the site exists there, create a third dataframe
    sim_df = pd.DataFrame()
    bench_df = pd.DataFrame()
    ref_df = pd.DataFrame()
    for ss in range(len(available_sites)):
        # simulations currently being evaluated
        cur_site = available_sites[ss]
        # todo: duplicate code that can be moved to a common function
        sim_df_cur = pd.read_csv(os.path.join(simulation_output_filepath, cur_site, 'inc_prev_data_final.csv'))
        upper_ages = sorted(sim_df_cur['Age'].unique())
        sim_df_cur['mean_age'] = sim_df_cur['Age'].apply(get_mean_from_upper_age, upper_ages=upper_ages)
        sim_df_cur['p_detect_case'] = coord_csv[coord_csv['site'] == cur_site]['p_detect_case'].iloc[0]

        # simulations used as benchmark
        if (benchmark_simulation_filepath is not None) and \
                (os.path.isfile(os.path.join(benchmark_simulation_filepath, cur_site, 'inc_prev_data_final.csv'))):
            bench_df_cur = pd.read_csv(os.path.join(benchmark_simulation_filepath, cur_site, 'inc_prev_data_final.csv'))
            upper_ages = sorted(bench_df_cur['Age'].unique())
            bench_df_cur['mean_age'] = bench_df_cur['Age'].apply(get_mean_from_upper_age, upper_ages=upper_ages)
            bench_df_cur['p_detect_case'] = coord_csv[coord_csv['site'] == cur_site]['p_detect_case'].iloc[0]
        else:
            bench_df_cur = pd.DataFrame()

        # reference data
        filepath_ref = os.path.join(base_reference_filepath,
                                    coord_csv[coord_csv['site'] == cur_site]['age_incidence_ref'].iloc[0])
        ref_df_cur = pd.read_csv(filepath_ref)
        ref_df_cur = ref_df_cur[ref_df_cur['Site'].str.lower() == cur_site.lower()]

        # add site into larger dataframe
        sim_df = pd.concat([sim_df, sim_df_cur])
        ref_df = pd.concat([ref_df, ref_df_cur])
        bench_df = pd.concat([bench_df, bench_df_cur])

    # scale down incidence in simulation according to probability of detecting a case in the reference setting
    sim_df['Incidence'] = sim_df['Incidence'] * sim_df['p_detect_case']

    # set up reference and simulation dataset columns
    ref_df['Incidence'] = ref_df['INC'] / 1000
    ref_df['mean_age'] = (ref_df['INC_LAR'] + ref_df['INC_UAR']) / 2
    ref_df = pd.DataFrame({'reference': ref_df['Incidence'],
                           'mean_age': ref_df['mean_age'],
                           'Site': ref_df['Site'],
                           'ref_pop_size': ref_df['POP'],
                           'ref_year': ref_df['START_YEAR']})
    sim_df = pd.DataFrame({'simulation': sim_df['Incidence'],
                           'mean_age': sim_df['mean_age'],
                           'Site': sim_df['Site']})
    # format benchmark simulations
    if len(bench_df) > 0:
        # scale down incidence in simulation according to probability of detecting a case
        bench_df['Incidence'] = bench_df['Incidence'] * bench_df['p_detect_case']
        bench_df = pd.DataFrame({'benchmark': bench_df['Incidence'],
                                 'mean_age': bench_df['mean_age'],
                                 'Site': bench_df['Site']})

    # check that ages match between reference and simulation. if there is a small difference (<1 year), update simulation to match reference age value after giving warning
    sim_df, bench_df = match_sim_ref_ages(ref_df, sim_df, bench_df)

    combined_df = pd.merge(ref_df, sim_df, how='outer')
    combined_df = pd.merge(combined_df, bench_df, how='outer')
    combined_df['metric'] = 'incidence'

    return combined_df


# prepare dataframe with simulation and reference data formatted together
def prepare_prev_df(coord_csv, simulation_output_filepath, base_reference_filepath, benchmark_simulation_filepath=None):
    """
    Read in, align, and combine reference and simulation data for all sites associated with the prevalence-by-age
    validation relationship.
    Args:
        coord_csv (): A dataframe detailing the sites simulated for each validation relationship and the corresponding
                     reference dataset
        simulation_output_filepath (): The filepath where simulation output is located
        base_reference_filepath (): The filepath where reference datasets are located
        benchmark_simulation_filepath (): The filepath where benchmark simulation output is located. If NA, no
                                          comparisons are made against benchmark simulations

    Returns: A dataframe containing the combined reference and simulation data for this validation relationship

    """

    # determine which of the age-prevalence sites have the relevant simulation output
    available_sites = get_available_sites_for_relationship(coord_csv, simulation_output_filepath,
                                                           relationship_name='age_prevalence',
                                                           relationship_sim_filename='prev_inc_by_age_month.csv')

    # aggregate all age-prevalence simulation data into one dataframe and reference data into a second dataframe (for the relevant sites)
    #   if a benchmark simulation directory was specified and the site exists there, create a third dataframe
    sim_df = pd.DataFrame()
    bench_df = pd.DataFrame()
    ref_df = pd.DataFrame()
    for ss in range(len(available_sites)):
        # simulations currently being evaluated
        cur_site = available_sites[ss]

        # read in and format reference data for this site
        filepath_ref = os.path.join(base_reference_filepath,
                                    coord_csv[coord_csv['site'] == cur_site]['age_prevalence_ref'].iloc[0])
        ref_df_cur = pd.read_csv(filepath_ref)
        ref_df_cur = ref_df_cur[ref_df_cur['Site'].str.lower() == cur_site.lower()]

        if 'agebin' in ref_df_cur.columns:
            upper_ages = sorted(ref_df_cur['agebin'].unique())
            ref_df_cur['mean_age'] = ref_df_cur['agebin'].apply(get_mean_from_upper_age, upper_ages=upper_ages)
        elif ('PR_LAR' in ref_df_cur.columns) and ('PR_UAR' in ref_df_cur.columns):
            ref_df_cur['mean_age'] = (ref_df_cur['PR_LAR'] + ref_df_cur['PR_UAR']) / 2
        name_dict = {'PR_MONTH': 'month',
                     'PR': 'prevalence',
                     'N': 'total_sampled',
                     'N_POS': 'num_pos',
                     'PR_YEAR': 'year'}

        ref_df_cur.rename(columns=name_dict, inplace=True)

        ref_df_cur = ref_df_cur[["Site", "mean_age", 'month', 'total_sampled', 'num_pos', 'prevalence', 'year']]
        # remove reference rows without prevalence values
        ref_df_cur = ref_df_cur[~ref_df_cur['prevalence'].isna()]

        # read in and format data from simulations to match reference dataset
        sim_df_cur = pd.read_csv(os.path.join(simulation_output_filepath, cur_site, 'prev_inc_by_age_month.csv'))
        sim_df_cur.rename(columns={'PfPR': 'prevalence'}, inplace=True)

        upper_ages = sorted(sim_df_cur['agebin'].unique())
        sim_df_cur['mean_age'] = sim_df_cur['agebin'].apply(get_mean_from_upper_age, upper_ages=upper_ages)
        # remove rows with population = 0 (this is relevant for cohort simulations where there is only one age group
        # in each year and all other age groups are zero)
        sim_df_cur = sim_df_cur[sim_df_cur['Pop'] > 0]
        sim_df_cur['month'] = sim_df_cur['month'].astype(str)

        if benchmark_simulation_filepath is not None:
            # read in and format data from benchmark simulations to match reference dataset
            bench_df_cur = pd.read_csv(os.path.join(benchmark_simulation_filepath, cur_site, 'prev_inc_by_age_month.csv'))
            bench_df_cur.rename(columns={'PfPR': 'prevalence'}, inplace=True)

            upper_ages = sorted(bench_df_cur['agebin'].unique())
            bench_df_cur['mean_age'] = bench_df_cur['agebin'].apply(get_mean_from_upper_age, upper_ages=upper_ages)
            # remove rows with population = 0 (this is relevant for cohort simulations where there is only one age
            # group in each year and all other age groups are zero)
            bench_df_cur = bench_df_cur[bench_df_cur['Pop'] > 0]
            bench_df_cur['month'] = bench_df_cur['month'].astype(str)

        else:
            bench_df_cur = pd.DataFrame()

        # determine whether reference is for a single month or averaged over multiple months - subset from simulations to match
        if (len(ref_df_cur['month'].unique()) == 1) and \
                (not str(ref_df_cur['month'].iloc[0]).isdigit()):  # check whether multiple months are listed in a character string
            # todo: need code review
            # included_months = as.numeric(unlist(strsplit(ref_df_cur$month[1], ",")))
            included_months = ref_df_cur['month'].iloc[0].split(",")
            included_months = [str(int(x)) for x in included_months]

            sim_df_cur = sim_df_cur[sim_df_cur['month'].astype(str).isin(included_months)]
            bench_df_cur = bench_df_cur[bench_df_cur['month'].astype(str).isin(included_months)]
            if len(included_months) > 1:
                sim_df_cur['month'] = 'multiple'
                ref_df_cur['month'] = 'multiple'
            if len(bench_df_cur) > 0:
                bench_df_cur['month'] = 'multiple'

        elif all(ref_df_cur['month'].astype(str).str.isnumeric()):
            included_months = ref_df_cur['month'].unique()
            included_months = [str(int(month)) for month in included_months]
            sim_df_cur = sim_df_cur[sim_df_cur['month'].astype(str).isin(included_months)]
            bench_df_cur = bench_df_cur[bench_df_cur['month'].astype(str).isin(included_months)]
        else:
            warnings.warn(f'The month format in the {cur_site} reference dataset was not recognized.')

        sim_df_cur = sim_df_cur[["Site", "mean_age", 'month', 'prevalence', 'year', 'Run_Number']]
        if len(bench_df_cur) > 0:
            bench_df_cur = bench_df_cur[["Site", "mean_age", 'month', 'prevalence', 'year', 'Run_Number']]

        # add site into larger dataframe
        sim_df = pd.concat([sim_df, sim_df_cur])
        ref_df = pd.concat([ref_df, ref_df_cur])
        bench_df = pd.concat([bench_df, bench_df_cur])

    # format reference data
    ref_df['Site'] = ref_df['Site'].str.lower()
    ref_df = pd.DataFrame({'reference': ref_df['prevalence'],
                           'mean_age': ref_df['mean_age'],
                           'Site': ref_df['Site'],
                           'month': ref_df['month'],
                           'site_month': ref_df['Site'] + '_month' + ref_df['month'].astype('str'),
                           'total_sampled': ref_df['total_sampled'],
                           'num_pos': ref_df['num_pos'],
                           'ref_year': ref_df['year']})

    # format new simulation output
    # get simulation average across seeds
    # todo: add a common method to create sim_df and bench_df
    sim_df = sim_df.groupby(['Site', 'mean_age', 'month']).agg(
        prevalence=('prevalence', np.nanmean)
    ).reset_index()
    sim_df['Site'] = sim_df['Site'].str.lower()
    sim_df = pd.DataFrame({'simulation':  sim_df['prevalence'],
                           'mean_age': sim_df['mean_age'],
                           'Site': sim_df['Site'],
                           'month': sim_df['month'],
                           'site_month': sim_df['Site'] + '_month' + sim_df['month'].astype('str')})

    # format benchmark simulations
    if len(bench_df) > 0:
        # get simulation average across seeds
        bench_df = bench_df.groupby(['Site', 'mean_age', 'month']).agg(
            prevalence=('prevalence', np.nanmean)
        ).reset_index()
        bench_df['Site'] = bench_df['Site'].str.lower()
        bench_df = pd.DataFrame({'benchmark': bench_df['prevalence'],
                                 'mean_age': bench_df['mean_age'],
                                 'Site': bench_df['Site'],
                                 'month': bench_df['month'],
                                 'site_month': bench_df['Site'] + '_month' + bench_df['month'].astype('str')})

    # check that ages match between reference and simulation. if there is a small difference (<1 year, update simulation)
    sim_df, bench_df = match_sim_ref_ages(ref_df, sim_df, bench_df)

    combined_df = pd.merge(ref_df, sim_df, on=['mean_age', 'site_month'], how='outer')
    combined_df = pd.merge(bench_df, combined_df, on=['mean_age', 'site_month'], how='outer')
    combined_df['metric'] = 'prevalence'
    return combined_df


def prepare_dens_df(coord_csv, simulation_output_filepath, base_reference_filepath, benchmark_simulation_filepath=None):
    """
    Read in, align, and combine reference and simulation data for all sites associated with the parasite density-by-age
    validation relationship.
    Args:
        coord_csv (): A dataframe detailing the sites simulated for each validation relationship and the corresponding
                      reference dataset
        simulation_output_filepath (): The filepath where simulation output is located
        base_reference_filepath (): The filepath where reference datasets are located
        benchmark_simulation_filepath (): The filepath where benchmark simulation output is located. If NA, no
                                          comparisons are made against benchmark simulations

    Returns: A dataframe containing the combined reference and simulation data for this validation relationship.

    """

    # determine which of the parasite density sites have the relevant simulation output
    available_sites = get_available_sites_for_relationship(coord_csv, simulation_output_filepath,
                                                           relationship_name='age_parasite_density',
                                                           relationship_sim_filename='parasite_densities_by_age_month.csv')

    if not available_sites:
        return pd.DataFrame(), pd.DataFrame()
    # iterate through sites, grabbing relevant reference and simulation data to plot; combine data into a dataframe containing all sites
    sim_df = pd.DataFrame()
    bench_df = pd.DataFrame()
    ref_df = pd.DataFrame()
    for ss in range(len(available_sites)):
        # simulations currently being evaluated
        cur_site = available_sites[ss]

        # read in and format reference data for this site
        # todo: write a common method to generate age_agg_df for sim, ref and benchmark data
        filepath_sim = os.path.join(simulation_output_filepath, cur_site, 'parasite_densities_by_age_month.csv')
        sim_df_cur = pd.read_csv(filepath_sim)
        upper_ages = sorted(sim_df_cur['agebin'].unique())
        sim_df_cur['mean_age'] =sim_df_cur['agebin'].apply(get_mean_from_upper_age, upper_ages=upper_ages)
        age_agg_sim_df = get_age_bin_averages(sim_df_cur)

        filepath_ref = os.path.join(base_reference_filepath,
                                    coord_csv[coord_csv['site'] == cur_site]['age_parasite_density_ref'].iloc[0])
        ref_df_cur = pd.read_csv(filepath_ref)
        ref_df_cur = ref_df_cur[ref_df_cur['Site'].str.lower() == cur_site.lower()]
        ref_df_cur['Site'] =ref_df_cur['Site'].str.lower()
        upper_ages = sorted(ref_df_cur['agebin'].unique())
        ref_df_cur['mean_age'] =ref_df_cur['agebin'].apply(get_mean_from_upper_age, upper_ages=upper_ages)

        if benchmark_simulation_filepath is not None:
            bench_df_cur = pd.read_csv(os.path.join(benchmark_simulation_filepath, cur_site,
                                                    'parasite_densities_by_age_month.csv'))
            upper_ages = sorted(bench_df_cur['agebin'].unique())
            bench_df_cur['mean_age'] = bench_df_cur['agebin'].apply(get_mean_from_upper_age, upper_ages=upper_ages)
            age_agg_bench_df = get_age_bin_averages(bench_df_cur)

            if sorted(age_agg_sim_df['mean_age'].unique()) != sorted(age_agg_bench_df['mean_age'].unique()):
                warnings.warn(f'New and benchmark simulation age bins are not the same for site: {cur_site}. '
                              f'Removing benchmark sims.')
                age_agg_bench_df = pd.DataFrame()
            if sorted(age_agg_sim_df['densitybin'].unique()) != sorted(age_agg_bench_df['densitybin'].unique()):
                warnings.warn(f'New and benchmark simulation parasite density bins are not the same for site: {cur_site}. '
                              f'Removing benchmark sims.')
                age_agg_bench_df = pd.DataFrame()
        else:
            age_agg_bench_df = pd.DataFrame()


        # subset simulation output to months in reference dataset
        months = sorted(ref_df_cur['month'].unique())
        sim_df_cur = age_agg_sim_df[age_agg_sim_df['month'].isin(months)]
        if len(age_agg_bench_df) > 0:
            bench_df_cur = age_agg_bench_df[age_agg_bench_df['month'].isin(months)]
        else:
            bench_df_cur = pd.DataFrame()

        # if the maximum reference density bin is < (maximum simulation density bin / max_magnitude_difference), aggregate all simulation densities >= max ref bin into the max ref bin
        #   the final bin will be all densities equal to or above that value
        max_ref_dens = ref_df_cur['densitybin'].max(skipna=True)
        sim_df_cur = combine_higher_dens_freqs(sim_df_cur, max_ref_dens, max_magnitude_difference=100)
        if len(bench_df_cur) > 0:
            bench_df_cur = combine_higher_dens_freqs(bench_df_cur, max_ref_dens, max_magnitude_difference=100)


        # add zeros for unobserved reference densities up to max_ref_dens
        all_zeros_df = sim_df_cur[['month', 'mean_age', 'agebin', 'densitybin', 'Site']]
        ref_df_cur = pd.merge(ref_df_cur, all_zeros_df, how='outer')
        ref_df_cur.fillna(0, inplace=True)

        # add site into larger dataframe
        if len(ref_df) > 0:
            ref_df = pd.merge(ref_df, ref_df_cur, how='outer')
        else:
            ref_df = ref_df_cur
        sim_df = pd.concat([sim_df, sim_df_cur])
        bench_df = pd.concat([bench_df, bench_df_cur])

    # check that ages match between reference and simulation. if there is a small difference (<1 year, update simulation)
    sim_df, bench_df = match_sim_ref_ages(ref_df, sim_df, bench_df)

    # format reference data
    ref_df_asex = pd.DataFrame({'reference': ref_df['asexual_par_dens_freq'],
                                'mean_age': ref_df['mean_age'],
                                'agebin': ref_df['agebin'],
                                'densitybin': ref_df['densitybin'],
                                'Site': ref_df['Site'],
                                'month': ref_df['month'],
                                'site_month': ref_df['Site'] + '_month' + ref_df['month'].astype('str'),
                                'ref_total': ref_df['bin_total_asex'],
                                'ref_bin_count': ref_df['count_asex']})
    ref_df_gamet = pd.DataFrame({'reference': ref_df['gametocyte_dens_freq'],
                                'mean_age': ref_df['mean_age'],
                                'agebin': ref_df['agebin'],
                                'densitybin': ref_df['densitybin'],
                                'Site': ref_df['Site'],
                                'month': ref_df['month'],
                                'site_month': ref_df['Site'] + '_month' + ref_df['month'].astype('str'),
                                'ref_total': ref_df['bin_total_gamet'],
                                'ref_bin_count': ref_df['count_gamet']})

    # format new simulation output
    sim_df_asex = pd.DataFrame({'simulation': sim_df['asexual_par_dens_freq'],
                                'mean_age': sim_df['mean_age'],
                                'agebin': sim_df['agebin'],
                                'densitybin': sim_df['densitybin'],
                                'Site': sim_df['Site'],
                                'month': sim_df['month'],
                                'site_month': sim_df['Site'] + '_month' + sim_df['month'].astype('str')})
    sim_df_gamet = pd.DataFrame({'simulation': sim_df['gametocyte_dens_freq'],
                                'mean_age': sim_df['mean_age'],
                                'agebin': sim_df['agebin'],
                                'densitybin': sim_df['densitybin'],
                                'Site': sim_df['Site'],
                                'month': sim_df['month'],
                                'site_month': sim_df['Site'] + '_month' + sim_df['month'].astype('str')})

    # format benchmark simulations
    if len(bench_df) > 0:
        bench_df_asex = pd.DataFrame({'benchmark': bench_df['asexual_par_dens_freq'],
                                    'mean_age': bench_df['mean_age'],
                                    'agebin': bench_df['agebin'],
                                    'densitybin': bench_df['densitybin'],
                                    'Site': bench_df['Site'],
                                    'month': bench_df['month'],
                                    'site_month': bench_df['Site'] + '_month' + bench_df['month'].astype('str')})
        bench_df_gamet = pd.DataFrame({'benchmark': bench_df['gametocyte_dens_freq'],
                                      'mean_age': bench_df['mean_age'],
                                      'agebin': bench_df['agebin'],
                                      'densitybin': bench_df['densitybin'],
                                      'Site': bench_df['Site'],
                                      'month': bench_df['month'],
                                      'site_month': bench_df['Site'] + '_month' + bench_df['month'].astype('str')})

    # combine reference and simulation dataframes
    combined_df_asex = pd.merge(ref_df_asex, sim_df_asex, how='outer')
    if len(bench_df) > 0:
        combined_df_asex = pd.merge(combined_df_asex, bench_df_asex, how='outer')
    combined_df_asex['metric'] = 'asexual_density'

    combined_df_gamet = pd.merge(ref_df_gamet, sim_df_gamet, how='outer')
    if len(bench_df) > 0:
        combined_df_gamet = pd.merge(combined_df_gamet, bench_df_gamet, how='outer')
    combined_df_gamet['metric'] = 'gametocyte_density'

    return combined_df_asex, combined_df_gamet


def prepare_infect_df(coord_csv, simulation_output_filepath, base_reference_filepath, benchmark_simulation_filepath=None):
    """
    Read in, align, and combine reference and simulation data for all sites associated with the
    infectiousness-to-mosquitos validation relationship.
    Args:
        coord_csv ():
        simulation_output_filepath ():
        base_reference_filepath ():
        benchmark_simulation_filepath ():

    Returns: A dataframe containing the combined reference and simulation data for this validation relationship

    """

    # determine which of the infectiousness sites have the relevant simulation output
    available_sites = get_available_sites_for_relationship(coord_csv, simulation_output_filepath,
                                                           relationship_name='infectiousness_to_mosquitos',
                                                           relationship_sim_filename='infectiousness_by_age_density_month.csv')

    if not available_sites:
        return pd.DataFrame()
    # iterate through sites, grabbing relevant reference and simulation data to plot; combine data into a dataframe containing all sites
    sim_df = pd.DataFrame()
    bench_df = pd.DataFrame()
    ref_df = pd.DataFrame()
    for ss in range(len(available_sites)):
        # simulations currently being evaluated
        cur_site = available_sites[ss]

        filepath_ref = os.path.join(base_reference_filepath,
                                    coord_csv[coord_csv['site'] == cur_site]['infectiousness_to_mosquitos_ref'].iloc[0])
        ref_df_cur = pd.read_csv(filepath_ref)
        ref_df_cur = ref_df_cur[ref_df_cur['site'].str.lower() == str(cur_site).lower()]
        ref_months = ref_df_cur['month'].unique()

        filepath_sim = os.path.join(simulation_output_filepath, cur_site, 'infectiousness_by_age_density_month.csv')
        sim_df_cur = pd.read_csv(filepath_sim)
        # remove simulation rows with zero pop
        sim_df_cur = sim_df_cur[sim_df_cur['Pop'] > 0]
        # subset simulation to months in reference df
        sim_df_cur = sim_df_cur[sim_df_cur['month'].isin(ref_months)]
        # get mean (across simulation seeds and ages within a bin) fraction of individuals in each  {age bin, month,
        # densitybin, run number} group that fall in each infectiousness bin
        sim_df_agg2 = get_fraction_in_infectious_bin(sim_df_cur)

        if benchmark_simulation_filepath is not None:
            bench_df_cur = pd.read_csv(os.path.join(benchmark_simulation_filepath, cur_site,
                                                    'infectiousness_by_age_density_month.csv'))
            # remove simulation rows with zero pop
            bench_df_cur = bench_df_cur[bench_df_cur['Pop'] > 0]
            # subset simulation to months in reference df
            bench_df_cur = bench_df_cur[bench_df_cur['month'].isin(ref_months)]
            # get mean (across simulation seeds and ages within a bin) fraction of individuals in each  {age bin, month, densitybin, run number} group that fall in each infectiousness bin
            bench_df_agg2 = get_fraction_in_infectious_bin(bench_df_cur)

            if sorted(sim_df_agg2['agebin'].unique()) != sorted(bench_df_agg2['agebin'].unique()):
                warnings.warn(f'New and benchmark simulation age bins are not the same for site: {cur_site}. '
                              f'Removing benchmark sims.')
                bench_df_agg2 = pd.DataFrame()
            if sorted(sim_df_agg2['densitybin'].unique()) != sorted(bench_df_agg2['densitybin'].unique()):
                warnings.warn(f'New and benchmark simulation parasite density bins are not the same for site: '
                              f'{cur_site}. Removing benchmark sims.')
                bench_df_agg2 = pd.DataFrame()
        else:
            bench_df_agg2 = pd.DataFrame()

        # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
        # standardize column names and merge simulation and reference data frames
        ref_df_cur.rename(columns={'site': 'Site'}, inplace=True)
        ref_df_cur = pd.DataFrame({'reference': ref_df_cur['freq_frac_infect'],
                                   'agebin': ref_df_cur['agebin'],
                                   'densitybin': ref_df_cur['densitybin'],
                                   'fraction_infected_bin': ref_df_cur['fraction_infected_bin'],
                                   'Site': ref_df_cur['Site'],
                                   'month': ref_df_cur['month'],
                                   'site_month': ref_df_cur['Site'] + '_month' + ref_df_cur['month'].astype('str'),
                                   'ref_total': ref_df_cur['num_in_group'],
                                   'ref_bin_count': ref_df_cur['count']})
        sim_df_cur = pd.DataFrame({'simulation': sim_df_agg2['infectiousness_bin_freq'],
                                   'agebin': sim_df_agg2['agebin'],
                                   'densitybin': sim_df_agg2['densitybin'],
                                   'fraction_infected_bin': sim_df_agg2['infectiousness_bin'],
                                   'Site': sim_df_agg2['Site'],
                                   'month': sim_df_agg2['month'],
                                   'site_month': sim_df_agg2['Site'] + '_month' + sim_df_agg2['month'].astype('str')})

        if len(bench_df_agg2) > 0:
            bench_df_cur = pd.DataFrame({'benchmark': bench_df_agg2['infectiousness_bin_freq'],
                                         'agebin': bench_df_agg2['agebin'],
                                         'densitybin': bench_df_agg2['densitybin'],
                                         'fraction_infected_bin': bench_df_agg2['infectiousness_bin'],
                                         'Site': bench_df_agg2['Site'],
                                         'month': bench_df_agg2['month'],
                                         'site_month': bench_df_agg2['Site'] + '_month' + bench_df_agg2['month'].astype('str')})

        # add site into larger dataframe
        sim_df = pd.concat([sim_df, sim_df_cur])
        ref_df = pd.concat([ref_df, ref_df_cur])
        # todo: I added this if condition since I think there is case that the bench_df_cur is not defined
        if len(bench_df_cur) > 0:
            bench_df = pd.concat([bench_df, bench_df_cur])

    combined_df = pd.merge(sim_df, ref_df, how='outer')
    if len(bench_df) > 0:
        combined_df = pd.merge(combined_df, bench_df, how='outer')
    combined_df['metric'] = 'infectiousness'
    return combined_df
# endregion


# region: infection duration sampling function
# Note: this function does not follow the same pattern as the functions for the other validation relationships
def get_sim_survey(sim_dir, ref_df, seeds=None):
    """
    Subsample from the simulation output to match the survey that generated the reference dataset (i.e., match the
    dates and ages of sampled individuals)
    If they have not already been generated, the subsampled results (one subsampling per simulation seed) are saved to
     a csv
    Args:
        sim_dir (): The filepath to the directory where the simulation results are saved
        ref_df (): A dataframe containing the reference data, which is used to determine which individuals from the
            simulations are kept in the subsampled results
        seeds (): The subset of simulaton run seeds to include. If NA, all seeds are used.

    Returns: A dataframe containing simulation survey results matching the reference dataset (includes a set of
            reference-matched rows for each of the simulation seeds)

    """

    sampled_sim_filename = 'sim_duration_survey_sampling.csv'
    file_path = os.path.join(sim_dir, sampled_sim_filename)
    if os.path.isfile(file_path):
        sim_subset_full = pd.read_csv(file_path)
    else:
        # get first year of sampling in reference dataset. the simulation will be referenced from the first day of that year
        first_ref_date = datetime.date(datetime.datetime.strptime(str(ref_df['date'].dropna().min()), "%Y-%m-%d %H:%M:%S").year, 1, 1)
        indIDs = ref_df['SID'].unique()
        ref_df['date'] = ref_df['date'].apply(lambda x: x.date())

        patient_report_path = os.path.join(sim_dir, 'patient_reports.csv')
        sim_full = pd.read_csv(patient_report_path)

        if seeds is None:
            seeds = sim_full['Run_Number'].unique()

        sim_subset_full = pd.DataFrame()
        for seed in sorted(seeds):
            print('Currently on seed ' + str(seed))
            sim = sim_full[sim_full['Run_Number'] == seed]  # subset to desired run
            sim['date'] = [first_ref_date + datetime.timedelta(days=int(simday)) for simday in sim['simday']]
            sim['age'] = sim['age'] / 365
            # track which individuals have already been included from the simulation (
            # to avoid double-sampling simulation individuals)
            included_ids = set()
            sim_subset = pd.DataFrame()
            for ii in range(len(indIDs)):
                if ii % 50 == 0:
                    print('Currently on individual ' + str(ii) + ' out of ', len(indIDs))
                ref_df_cur = ref_df[ref_df['SID'] == indIDs[ii]]
                ref_df_cur = ref_df_cur.sort_values(by='date')
                # find a matching individual
                age_cur = ref_df_cur['age'].iloc[0]
                day_cur = ref_df_cur['date'].iloc[0]

                # use age-specific matches
                id_candidates = sim[(sim['date'] == day_cur) & (round(sim['age']) == round(age_cur))]['id'].tolist()
                id_candidates = [idx for idx in id_candidates if idx not in included_ids]
                # if no perfect age-match remain, expand year-range until finding a match
                if len(id_candidates) == 0:
                    year_range = 0
                    while len(id_candidates) == 0 & year_range < 100:
                        year_range = year_range + 5
                        # todo: need some code review on the math
                        # r code
                        # id_candidates = sim$id[intersect(which(sim$date == day_cur), which(round(sim$age) % in % seq((round(age_cur)-year_range), (round(age_cur)+year_range))))]
                        id_candidates = sim[(sim['date'] == day_cur)
                                            & (sim['age'].round().isin(range((round(age_cur)-year_range),
                                                                              round(age_cur)+year_range)))]['id'].tolist()
                        id_candidates = [idx for idx in id_candidates if idx not in included_ids]

                    if len(id_candidates) == 0:
                        print('Problem: no age-matched simulation individual found for reference id: ' + indIDs[ii])
                    else:
                        print('No exact age match remaining for reference id: ' + indIDs[ii]
                              + '. Used simulation individual within ', year_range, ' years.')

                id_sim_cur = random.sample(id_candidates, 1)[0]  # todo: should we remove this id after drawing?
                included_ids.add(id_sim_cur)

                # keep the same simulation dates as the reference samples for this individual
                sim_subset_cur = sim[(sim['id'] == id_sim_cur) & (sim['date'].isin(ref_df_cur['date']))]
                sim_subset = pd.concat([sim_subset, sim_subset_cur])

            sim_subset['seed'] = seed
            if sim_subset_full.empty:
                sim_subset_full = sim_subset
            else:
                sim_subset_full = pd.concat([sim_subset_full, sim_subset])

        # rename simulation columns to match reference data
        sim_subset_full.rename(columns={'id': 'SID', 'true_asexual_parasites': 'DENSITY'}, inplace=True)
        sim_subset_full.to_csv(file_path, index=False)

    return sim_subset_full
