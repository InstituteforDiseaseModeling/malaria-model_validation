import os
import warnings
import pandas as pd
from functools import partial
import emod_api.demographics.Demographics as Demographics
import params as parameters
from typing import Dict, Any
from idmtools.entities.simulation import Simulation
from emodpy_malaria.interventions.diag_survey import add_diagnostic_survey

from emodpy_malaria.reporters.builtin import add_malaria_summary_report, MalariaPatientJSONReport
from emodpy_malaria import malaria_config as malconf
from emodpy_malaria.interventions.drug_campaign import add_drug_campaign
from emodpy_malaria.interventions.treatment_seeking import add_treatment_seeking
from emodpy_malaria.interventions.inputeir import InputEIR
from emod_api.interventions.common import BroadcastEvent
# import emod_api.campaign as camp
import simulations.manifest as manifest


def update_sim_random_seed(simulation, value):
    simulation.task.config.parameters.Run_Number = value
    return {"Run_Number": value}


def mAb_vs_EIR(EIR):
    # Rough cut at function from eyeballing a few BinnedReport outputs parsed into antibody fractions
    mAb = 0.9 * (1e-4*EIR*EIR + 0.7*EIR) / ( 0.7*EIR + 2 )
    return min(mAb, 1.0)


def update_mab(simulation, value):
    simulation.task.config.parameters.Maternal_Antibody_Protection *= value
    return None


def set_param_fn(config):
    """
    This function is a callback that is passed to emod-api.config to set parameters The Right Way.
    """
    config = malconf.set_team_defaults(config, manifest)
    # config = set_config.set_config(config)

    config.parameters.Base_Rainfall = 150
    config.parameters.Climate_Model = "CLIMATE_CONSTANT"
    config.parameters.Enable_Disease_Mortality = 0
    config.parameters.Enable_Vector_Species_Report = 0
    config.parameters.Enable_Vital_Dynamics = 0
    config.parameters.Demographics_Filenames = ['demographics.json']
    config.parameters.Simulation_Duration = 70*365
    config.parameters.Enable_Initial_Prevalence = 1
    config.parameters.Age_Initialization_Distribution_Type = 'DISTRIBUTION_SIMPLE'
    config.parameters.Vector_Species_Params = []
    config.parameters.Start_Time = 0
    # config.parameters.pop("Serialized_Population_Filenames")

    return config


# def update_camp_type(simulation, site):
#     # simulation.task.config.parameters.Run_Number = value
#     build_camp_partial = partial(build_camp, site=site)
#     simulation.task.create_campaign_from_callback(build_camp_partial)
#
#     update_mab(simulation, mAb_vs_EIR(sum(study_site_monthly_EIRs[site])))
#
#     return {"Site": site}


def set_simulation_scenario(simulation, site):
    # get information on this simulation setup from coordinator csv
    coord_df = pd.read_csv(manifest.simulation_coordinator_path)
    coord_df = coord_df.set_index('site')

    # === set up config === #
    # simulation duration
    simulation_duration = coord_df.at[site, 'simulation_duration'].tolist()
    simulation.task.config.parameters.Simulation_Duration = simulation_duration
    # set whether there are births and deaths
    simulation.task.config.parameters.Enable_Vital_Dynamics = coord_df.at[site, 'enable_vital_dynamics'].tolist()
    # maternal antibodies - use first 12 months of data frame to get annual EIR from monthly eir
    monthly_eirs = pd.read_csv(os.path.join(manifest.input_files_path, coord_df.at[site, 'EIR_filepath']))
    update_mab(simulation, mAb_vs_EIR(sum(monthly_eirs.loc[monthly_eirs.index[0:12], site])))

    # === set up campaigns === #
    build_camp_partial = partial(build_camp, site=site, coord_df=coord_df)
    simulation.task.create_campaign_from_callback(build_camp_partial)

    # === set up reporters === #
    if coord_df.at[site, 'include_AnnualMalariaSummaryReport']:
        if (not pd.isna(coord_df.at[site, 'annual_summary_report_age_bins_filepath'])) and (not (coord_df.at[site, 'annual_summary_report_age_bins_filepath'] == '')):
            summary_report_age_bins_df = pd.read_csv(os.path.join(manifest.input_files_path, coord_df.at[site, 'annual_summary_report_age_bins_filepath']))
            summary_report_age_bins = summary_report_age_bins_df['age'].tolist()
        else:
            summary_report_age_bins = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 19, 39, 59, 85]

        add_malaria_summary_report(simulation.task, manifest=manifest, start_day=0, duration_days=1000000,
                                   reporting_interval=365, age_bins=summary_report_age_bins,
                                   infectiousness_bins=[0, 100], max_number_reports=2000,
                                   parasitemia_bins=[0, 50, 500, 5000, 5000000], report_description='Annual_Report')

    if coord_df.at[site, 'include_MonthlyMalariaSummaryReport']:
        if (not pd.isna(coord_df.at[site, 'monthly_summary_report_age_bins_filepath'])) and (not (coord_df.at[site, 'monthly_summary_report_age_bins_filepath'] == '')):
            summary_report_age_bins_df = pd.read_csv(os.path.join(manifest.input_files_path, coord_df.at[site, 'monthly_summary_report_age_bins_filepath']))
            summary_report_age_bins = summary_report_age_bins_df['age'].tolist()
        else:
            summary_report_age_bins = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 19, 39, 59, 85]

        for yy in range(round(simulation_duration/365)):
            add_malaria_summary_report(simulation.task, manifest=manifest, start_day=(yy * 365), duration_days=365,
                                       reporting_interval=30, age_bins=summary_report_age_bins,
                                       infectiousness_bins=[0, 100], max_number_reports=1000,
                                       parasitemia_bins=[0, 50, 500, 5000, 50000, 5000000], report_description='Monthly_Report_%i' % yy)

    if coord_df.at[site, 'include_MalariaPatientReport']:
        patient_report = MalariaPatientJSONReport()  # Create the reporter
        patient_report.config(ptr_config_builder, manifest)  # Config the reporter
        simulation.task.reporters.add_reporter(patient_report)  # Add the reporter

    if coord_df.at[site, 'include_parDensSurveys']:  # surveys added as a campaign in build_camp()
        simulation.task.config.parameters.Report_Event_Recorder = 1
        simulation.task.config.parameters.Report_Event_Recorder_Events = ['parasites_on_survey_day']
        simulation.task.config.parameters.Custom_Individual_Events = ['parasites_on_survey_day']

    return {"Site": site}


def build_standard_campaign_object(manifest):
    import emod_api.campaign as campaign
    campaign.set_schema(manifest.schema_file)
    return campaign


# def build_camp(site, cross_sectional_surveys=False, survey_days=None):
def build_camp(site, coord_df):
    """
    Build a campaign input file for the DTK using emod_api.
    Right now this function creates the file and returns the filename. If calling code just needs an asset that's fine.
    """
    # create campaign object
    camp = build_standard_campaign_object(manifest)

    # === EIR === #

    # set monthly eir for site - TODO - change to daily EIR
    monthly_eirs = pd.read_csv(os.path.join(manifest.input_files_path, coord_df.at[site, 'EIR_filepath']))
    # TODO - currently recycles first 12 values; should update to use multiple years if provided
    camp.add(InputEIR(camp, monthly_eir=monthly_eirs.loc[monthly_eirs.index[0:12], site].tolist(),
             start_day=0, age_dependence="SURFACE_AREA_DEPENDENT"))


    # === INTERVENTIONS === #

    # health-seeking
    if (not pd.isna(coord_df.at[site, 'CM_filepath'])) and (not (coord_df.at[site, 'CM_filepath'] == '')):
        hs_df = pd.read_csv(os.path.join(manifest.input_files_path, coord_df.at[site, 'CM_filepath']))
    else:
        hs_df = pd.DataFrame()
    # NMFs
    if (not pd.isna(coord_df.at[site, 'NMF_filepath'])) and (not (coord_df.at[site, 'NMF_filepath'] == '')):
        nmf_df = pd.read_csv(os.path.join(manifest.input_files_path, coord_df.at[site, 'NMF_filepath']))
    else:
        nmf_df = pd.DataFrame()

    if not hs_df.empty:
        # case management for malaria
        add_hfca_hs(camp, hs_df)
        # case management for NMFs
        add_nmf_hs(camp, hs_df, nmf_df)


    # === SURVEYS === #

    # add parasite density surveys among individuals with parasitemia
    if coord_df.at[site, 'include_parDensSurveys']:
        # adding schema file, so it can be looked up when creating the campaigns
        camp.schema_path = manifest.schema_file
        survey_days = pd.read_csv(os.path.join(manifest.input_files_path, coord_df.at[site, 'survey_days_filepath'])).loc['days']
        add_broadcasting_survey(camp, survey_days=survey_days)

    return camp


def add_hfca_hs(camp, hs_df):
    for r, row in hs_df.iterrows() :
        add_hs_from_file(camp, row)


def add_hs_from_file(camp, row):
    hs_child = row['U5_coverage']
    hs_adult = row['adult_coverage']
    severe_cases = row['severe_coverage']
    start_day = row['simday']
    duration = row['duration']

    add_treatment_seeking(camp, start_day=start_day,
                          targets=[{'trigger': 'NewClinicalCase', 'coverage': hs_child, 'agemin': 0, 'agemax': 5,
                                   'seek': 1, 'rate': 0.3},
                                   {'trigger': 'NewClinicalCase', 'coverage': hs_adult, 'agemin': 5, 'agemax': 100,
                                    'seek': 1, 'rate': 0.3},],
                          drug=['Artemether', 'Lumefantrine'], duration=duration)
    add_treatment_seeking(camp, start_day=start_day,
                          targets=[{'trigger': 'NewSevereCase', 'coverage': severe_cases, 'seek': 1, 'rate': 0.5}], #change by adding column and reviewing literature
                          drug=['Artemether', 'Lumefantrine'], duration=duration)  # , broadcast_event_name='Received_Severe_Treatment')


def add_nmf_hs(camp, hs_df, nmf_df):
    # if no NMF rate is specified, assume all age groups have 0.0038 probability each day
    if nmf_df.empty:
        nmf_df = pd.DataFrame({'U5_nmf': [0.0038], 'adult_nmf': [0.0038]})
    elif nmf_df.shape[0] != 1:
        warnings.warn('The NMF dataframe has more than one row. Only values in the first row will be used.')
    nmf_row = nmf_df.iloc[0]

    # apply the health-seeking rate for clinical malaria to NMFs
    for r, row in hs_df.iterrows():
        add_nmf_hs_from_file(camp, row, nmf_row)


def add_nmf_hs_from_file(camp, row, nmf_row):
    hs_child = row['U5_coverage']
    hs_adult = row['adult_coverage']
    start_day = row['simday']
    duration = row['duration']
    if start_day == 0:  # due to dtk diagnosis/treatment configuration, a start day of 0 is not supported
        start_day = 1  # start looking for NMFs on day 1 (not day 0) of simulation
        if duration > 1:
            duration = duration - 1
    nmf_child = nmf_row['U5_nmf']
    nmf_adult = nmf_row['adult_nmf']

    add_drug_campaign(camp, 'MSAT', 'AL', start_days=[start_day],
                      target_group={'agemin': 0, 'agemax': 5},
                      coverage=nmf_child * hs_child,
                      repetitions=duration, tsteps_btwn_repetitions=1,
                      diagnostic_type='PF_HRP2', diagnostic_threshold=5,
                      receiving_drugs_event_name='Received_NMF_Treatment')
    add_drug_campaign(camp, 'MSAT', 'AL', start_days=[start_day],
                      target_group={'agemin': 5, 'agemax': 120},
                      coverage=nmf_adult * hs_adult,
                      repetitions=duration, tsteps_btwn_repetitions=1,
                      diagnostic_type='PF_HRP2', diagnostic_threshold=5,
                      receiving_drugs_event_name='Received_NMF_Treatment')


def ptr_config_builder(params):
    return params


def add_broadcasting_survey(camp, survey_days, include_neg_broadcast=False):
    pos_diag_cfg = BroadcastEvent(camp=camp, Event_Trigger='parasites_on_survey_day')
    neg_diag_cfg = BroadcastEvent(camp=camp, Event_Trigger='negative_test_on_survey_day')
    for survey_day in survey_days:
        if include_neg_broadcast:
            add_diagnostic_survey(camp=camp, start_day=survey_day,
                                  diagnostic_type='TRUE_PARASITE_DENSITY', diagnostic_threshold=0,
                                  positive_diagnosis_configs=[pos_diag_cfg], negative_diagnosis_configs=[neg_diag_cfg])
        else:
            add_diagnostic_survey(camp=camp, start_day=survey_day,
                                  diagnostic_type='TRUE_PARASITE_DENSITY', diagnostic_threshold=0,
                                  positive_diagnosis_configs=[pos_diag_cfg])


def build_demog():
    """
    Build a demographics input file for the DTK using emod_api.
    Right now this function creates the file and returns the filename. If calling code just needs an asset that's fine.
    Also right now this function takes care of the config updates that are required as a result of specific demog settings. We do NOT want the emodpy-disease developers to have to know that. It needs to be done automatically in emod-api as much as possible.
    TBD: Pass the config (or a 'pointer' thereto) to the demog functions or to the demog class/module.

    """
    demog = Demographics.from_file(manifest)

    return demog


