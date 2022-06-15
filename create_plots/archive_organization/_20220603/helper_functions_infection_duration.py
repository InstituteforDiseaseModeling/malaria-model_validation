#####################################################################################
# compare_ref_sim_infection_duration_Navrongo
# contact: mambrose
# Jan 2022
#
# Goal: simulate transmission in EMOD and sample to match reference dataset
#   Compare:
#   - probability a person goes from negative to positive between sample dates
#   - probability a person goes from positive to negative between sample dates
#   - fraction of samples positive
#   - among all individuals starting positive, how long until negative?
#   - among all individuals turning positive during study, how long until negative?
#####################################################################################

import numpy as np
import pandas as pd
import os.path as path
import datetime
import math
import random


def get_sim_survey(sim_dir, ref_df, seeds=np.nan):
    """
    Get subset of simulation dataset to match reference dataset for dates and ages of individuals of sampled individuals
    Args:
        sim_dir ():
        ref_df ():
        seeds ():

    Returns:

    """
    sampled_sim_filename = 'sim_duration_survey_sampling.csv'
    file_path = path.join(sim_dir, sampled_sim_filename)
    if path.isfile(file_path):
        sim_subset_full = pd.read_csv(file_path)
    else:
        # Get first year of sampling in reference dataset. the simulation will be referenced from the first day of
        # that year
        first_ref_date = datetime.date(datetime.datetime.strptime(ref_df['date'].dropna().min(), "%Y-%m-%d").year, 1, 1)
        indIDs = ref_df['SID'].unique()

        patient_report_path = path.join(sim_dir, 'patient_reports.csv')

        sim_full = pd.read_csv(patient_report_path)
        sim_full['date'] = first_ref_date + sim_full['simday']
        sim_full['age'] = sim_full['age'] / 365

        if math.isnan(seeds):
            seeds = sim_full['Run_Number'].unique()

        for seed in sorted(seeds):
            print('Currently on seed ' + seed)
            sim = sim_full[sim_full['Run_Number'] == seed]  # subset to desired run
            # track which individuals have already been included from the simulation (
            # to avoid double-sampling simulation individuals)
            included_ids = list()
            sim_subset = pd.DataFrame()
            for ii in range(len(indIDs)):
                if ii % 50 == 0:
                    print('Currently on individual ' + ii + ' out of ', len(indIDs))
                ref_df_cur = ref_df[ref_df['SID'] == indIDs[ii]]
                ref_df_cur = ref_df_cur.sort_values(by='date')
                # find a matching individual
                age_cur = ref_df_cur['age'].iloc[0]
                day_cur = ref_df_cur['date'].iloc[0]

                # use age-specific matches
                id_candidates = sim[(sim['date'] == day_cur) & (round(sim['age'] == round(age_cur)))]['id']
                id_candidates = [idx not in included_ids for idx in id_candidates]
                # if no perfect age-match remain, expand year-range until finding a match
                if len(id_candidates) == 0:
                    year_range = 0
                    while len(id_candidates) == 0 & year_range < 100:
                        year_range = year_range + 5
                        id_candidates = sim[(sim['date'] == day_cur)
                                            & (sim['age'].round().isin(range((round(age_cur)-year_range),
                                                                             round(age_cur)+year_range)))]['id']
                        id_candidates = [idx not in included_ids for idx in id_candidates]

                    if len(id_candidates) == 0:
                        print('Problem: no age-matched simulation individual found for reference id: ' + indIDs[ii])
                    else:
                        print('No exact age match remaining for reference id: ' + indIDs[ii]
                              + '. Used simulation individual within ', year_range, ' years.')

                id_sim_cur = random.sample(id_candidates, 1) # todo: should we remove this id after drawing?
                included_ids.extend(id_sim_cur)

                # keep the same simulation dates as the reference samples for this individual
                sim_subset_cur = sim[(sim['id'] == id_sim_cur) & (sim['date'].isin(ref_df_cur['date']))]
                sim_subset = pd.concat(sim_subset, sim_subset_cur)

            sim_subset['seed'] = seed
            if seed == sorted(seeds)[1]:
                sim_subset_full = sim_subset
            else:
                sim_subset_full = pd.concat(sim_subset_full, sim_subset)

        # rename simulation columns to match reference data
        sim_subset_full.rename(columns={'id': 'SID', 'true_asexual_parasites': 'DENSITY'}, inplace=True)
        sim_subset_full.to_csv(file_path, index=False)

    return sim_subset_full


# def get_frac_state_swaps(data):
#     """
#     Calculate probability of going from negative to positive or from positive to negative between sample dates
#     Args:
#         data ():
#
#     Returns:
#
#     """
#     # brute force approach iterating through people and dates
#     indIDs = data['SID'].unique()
#     sum_denom_pos = 0
#     sum_turn_neg = 0
#     sum_denom_neg = 0
#     sum_turn_pos = 0
#
#     for ii in range(len(indIDs)):
#         data_cur = data[data['SID'] == indIDs[ii]]
#         data_cur = data_cur.sort_values(by=['date'])
#         # get indices of positive tests and of negative tests
#         ind_pos = which(data_cur$DENSITY > pos_thresh_dens)
#         ind_neg = which(data_cur$DENSITY <= pos_thresh_dens)
#
#         # denominators for each (number of each type that had an observation after then (i.e., they could have been observed to change))
#         last_obs_pos = (data_cur$DENSITY[nrow(data_cur)] > pos_thresh_dens)
#         sum_denom_pos = sum_denom_pos + ifelse(last_obs_pos, length(ind_pos)-1, length(ind_pos))
#         sum_denom_neg = sum_denom_neg + ifelse(last_obs_pos, length(ind_neg), length(ind_neg)-1)
#
#         # find how many tests change from neg to pos or pos to neg across timesteps
#         sum_turn_neg = sum_turn_neg + sum((ind_pos + 1) % in % ind_neg)
#         sum_turn_pos = sum_turn_pos + sum((ind_neg + 1) % in % ind_pos)
#
#
#     frac_pos_turn_neg_next_time = sum_turn_neg / sum_denom_pos
#     frac_neg_turn_pos_next_time = sum_turn_pos / sum_denom_neg
#
#     return (c(frac_pos_turn_neg_next_time, frac_neg_turn_pos_next_time))
