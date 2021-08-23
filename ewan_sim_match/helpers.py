import ewan_sim_match.manifest as manifest

import emod_api.demographics.Demographics as Demographics

from emodpy_malaria import config as malconf
from emodpy_malaria.interventions.treatment_seeking import add
from emodpy_malaria.interventions.inputeir import InputEIR

import emod_api.campaign as camp


def update_sim_random_seed(simulation, value):
    simulation.task.config.parameters.Run_Number = value
    return {"Run_Number": value}


def mAb_vs_EIR(EIR):
    # Rough cut at function from eyeballing a few BinnedReport outputs parsed into antibody fractions
    mAb = 0.9 * (1e-4*EIR*EIR + 0.7*EIR) / ( 0.7*EIR + 2 )
    return min(mAb, 1.0)


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
    config.parameters.Simulation_Duration = 2*365
    config.parameters.Enable_Initial_Prevalence = 1
    # config.parameters.pop("Serialized_Population_Filenames")

    return config


def build_camp():
    """
    Build a campaign input file for the DTK using emod_api.
    Right now this function creates the file and returns the filename. If calling code just needs an asset that's fine.
    """

    # This isn't desirable. Need to think about right way to provide schema (once)
    camp.schema_path = manifest.schema_file

    # print( f"Telling emod-api to use {manifest.schema_file} as schema." )
    add(camp, targets=[{"trigger": "NewClinicalCase", "coverage": 1, "seek": 0.5, "rate": 0.3}],
        drug=['Artemether', 'Lumefantrine'], start_day=0, broadcast_event_name='Received_Treatment')

    camp.add(InputEIR(camp, monthly_eir=[10.4, 13, 6, 2.6, 4, 6, 35, 21, 28, 15, 10.4, 8.4],
             start_day=0, age_dependence="SURFACE_AREA_DEPENDENT"))

    return camp


def build_demog():
    """
    Build a demographics input file for the DTK using emod_api.
    Right now this function creates the file and returns the filename. If calling code just needs an asset that's fine.
    Also right now this function takes care of the config updates that are required as a result of specific demog settings. We do NOT want the emodpy-disease developers to have to know that. It needs to be done automatically in emod-api as much as possible.
    TBD: Pass the config (or a 'pointer' thereto) to the demog functions or to the demog class/module.

    """
    demog = Demographics.from_file("/ewan_sim_match/input_files/demographics.json")

    return demog


def msr_config_builder(params):
    params.Age_Bins = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 20, 25, 30, 40, 50, 60, 100]
    params.Duration_Days = 1000000
    params.Infectiousness_Bins = [0, 100]
    params.Max_Number_Reports = 2000
    params.Nodeset_Config = {"class": "NodeSetAll"}
    params.Parasitemia_Bins = [0, 50, 500, 5000, 5000000]
    params.Report_Description = 'Annual_Report'
    params.Reporting_Interval = 365
    params.Start_Day = 0

    return params

