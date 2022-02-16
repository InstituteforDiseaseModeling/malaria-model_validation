import emod_api.demographics.Demographics as Demographics
import manifest_asset as manifest_asset
import params as parameters

def build_demog():
    """
    Build a demographics input file for the DTK using emod_api.
    Right now this function creates the file and returns the filename. If calling code just needs an asset that's fine.
    Also right now this function takes care of the config updates that are required as a result of specific demog settings. We do NOT want the emodpy-disease developers to have to know that. It needs to be done automatically in emod-api as much as possible.
    TBD: Pass the config (or a 'pointer' thereto) to the demog functions or to the demog class/module.

    """
    demog = Demographics.from_file(manifest_asset)

    return demog



def msr_config_builder(params):
    params.Age_Bins = parameters.summary_report_age_bins
    params.Duration_Days = 1000000
    params.Infectiousness_Bins = [0, 100]
    params.Max_Number_Reports = 2000
    params.Nodeset_Config = {"class": "NodeSetAll"}
    params.Parasitemia_Bins = [0, 50, 500, 5000, 5000000]
    params.Report_Description = 'Annual_Report'
    params.Reporting_Interval = 365
    params.Start_Day = 0

    return params

