from functools import partial
from typing import Dict, Any
from idmtools.entities.simulation import Simulation
from emodpy_malaria.interventions.diag_survey import add_diagnostic_survey

from emodpy_malaria import malaria_config as malconf
from emodpy_malaria.interventions.treatment_seeking import add_treatment_seeking
from emodpy_malaria.interventions.inputeir import InputEIR
from emod_api.interventions.common import BroadcastEvent
# import emod_api.campaign as camp
import general.manifest as manifest
from general.site_eir_values import study_site_monthly_EIRs, mAb_vs_EIR, sites_with_interventions, sites_cm_seek


def update_sim_random_seed(simulation, value):
    simulation.task.config.parameters.Run_Number = value
    return {"Run_Number": value}


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



# def update_camp_type(simulation, site, cross_sectional_surveys=False, survey_days=None):
#     # simulation.task.config.parameters.Run_Number = value
#     build_camp_partial = partial(build_camp, site=site, cross_sectional_surveys=cross_sectional_surveys, survey_days=survey_days)
#     simulation.task.create_campaign_from_callback(build_camp_partial)
#
#     update_mab(simulation, mAb_vs_EIR(sum(study_site_monthly_EIRs[site])))
#
#     return {"Site": site}


def update_camp_type(simulation, site):
    # simulation.task.config.parameters.Run_Number = value
    build_camp_partial = partial(build_camp, site=site)
    simulation.task.create_campaign_from_callback(build_camp_partial)

    update_mab(simulation, mAb_vs_EIR(sum(study_site_monthly_EIRs[site])))

    return {"Site": site}




def build_standard_campaign_object(manifest):
    import emod_api.campaign as campaign
    campaign.set_schema(manifest.schema_file)
    return campaign

# def build_camp(site, cross_sectional_surveys=False, survey_days=None):
def build_camp(site):
    """
    Build a campaign input file for the DTK using emod_api.
    Right now this function creates the file and returns the filename. If calling code just needs an asset that's fine.
    """

    # This isn't desirable. Need to think about right way to provide schema (once)
    # camp.schema_path = manifest.schema_file
    camp = build_standard_campaign_object(manifest)

    if site not in study_site_monthly_EIRs.keys():
        raise Exception("Don't know how to configure site: %s " % site)

    # print( f"Telling emod-api to use {manifest.schema_file} as schema." )
    if site in sites_with_interventions:
        add_treatment_seeking(camp, targets=[{"trigger": "NewClinicalCase", "coverage": 1, "seek": sites_cm_seek[site], "rate": 0.3}],
            drug=['Artemether', 'Lumefantrine'], start_day=0, broadcast_event_name='Received_Treatment')

    camp.add(InputEIR(camp, monthly_eir=study_site_monthly_EIRs[site],
             start_day=0, age_dependence="SURFACE_AREA_DEPENDENT"))

    # if cross_sectional_surveys:
    #     # adding schema file, so it can be looked up when creating the campaigns
    #     camp.schema_path = manifest.schema_file
    #     add_broadcasting_survey(camp, survey_days=survey_days)

    return camp


def ptr_config_builder(params):
    return params

def set_param(simulation: Simulation, param: str, value: Any) -> Dict[str, Any]:
    """
    Set specific parameter value
    Args:
        simulation: idmtools Simulation
        param: parameter
        value: new value

    Returns: dict

    """
    return simulation.task.set_parameter(param, value)


def add_broadcasting_survey(camp, survey_days):
    pos_diag_cfg = BroadcastEvent(campaign=camp, Event_Trigger='parasites_on_survey_day')
    # neg_diag_cfg = BroadcastEvent(campaign, Event_Trigger='negative_test_on_survey_day')
    for survey_day in survey_days:
        add_diagnostic_survey(campaign=camp, start_day=survey_day,
                              diagnostic_type='TRUE_PARASITE_DENSITY', diagnostic_threshold=0,
                              positive_diagnosis_configs=[pos_diag_cfg])  #, negative_diagnosis_configs=[neg_diag_cfg])

    return

