# use the coordinator csv to determine which analyzers should be run for each site
import pandas as pd
import simulations.manifest as manifest

from idmtools.core.platform_factory import Platform
from idmtools.analysis.platform_anaylsis import PlatformAnalysis
from simulations.analyzers.AnnualSummaryReportAnalyzer import AnnualSummaryReportAnalyzer
from simulations.analyzers.ParDensAgeAnalyzer import ParDensAgeAnalyzer
from simulations.analyzers.InfectiousnessByParDensAgeAnalyzer import InfectiousnessByParDensAgeAnalyzer
from simulations.analyzers.PatientReportAnalyzer import PatientAnalyzer
from simulations.analyzers.MonthlySummaryReportAnalyzer import MonthlySummaryReportAnalyzer


# TODO: replace manual specification of exp ids once we are using a single thread to track experiments
exp_name_id = {
    # 'validation_navrongo_2000': '05741b4b-85c2-ec11-a9f6-9440c9be2c51',
    # 'validation_sotuba_1999': 'ab993a49-84c2-ec11-a9f6-9440c9be2c51',
    # 'validation_dongubougou_1999': 'fb3f7607-84c2-ec11-a9f6-9440c9be2c51',
    'validation_dapelogo_2007': '592cb8ef-30c8-ec11-a9f8-b88303911bc1',
    # 'validation_sugungum_1970': '37cf618d-78c2-ec11-a9f6-9440c9be2c51',
    # 'validation_ebolakounou_1997': 'de036687-78c2-ec11-a9f6-9440c9be2c51',
    # 'validation_dielmo_1990': '7eeed80f-f7c7-ec11-a9f8-b88303911bc1',
    # 'validation_chonyi_1999': '19324580-78c2-ec11-a9f6-9440c9be2c51',
    # 'validation_ndiop_1993': 'a35829dc-e3c7-ec11-a9f8-b88303911bc1',
    # 'validation_ngerenya_1999': '810f4f7c-0ec0-ec11-a9f6-9440c9be2c51',
    # 'validation_matsari_1970': '26ebf0ce-70bf-ec11-a9f6-9440c9be2c51',
    # 'validation_rafin_marke_1970': '83c662c8-70bf-ec11-a9f6-9440c9be2c51',
    # 'validation_laye_2007': 'd75c7309-f7c7-ec11-a9f8-b88303911bc1',
}

if __name__ == "__main__":
    platform = Platform('CALCULON')
    coord_df = pd.read_csv(manifest.simulation_coordinator_path)
    coord_df = coord_df.set_index('site')

    for expt_name, id in exp_name_id.items():
        site = expt_name.replace('validation_', '')
        report_start_day = int(coord_df.at[site, 'report_start_day'])
        # determine the analyzers to run for each site
        analyzers = []
        analyzer_args = []
        if coord_df.at[site, 'include_MonthlyMalariaSummaryReport']:
            if coord_df.at[site, 'age_parasite_density']:
                analyzers.append(ParDensAgeAnalyzer)
                analyzer_args.append({'expt_name': site,
                                      'sweep_variables': ['Run_Number', 'Site'],
                                      'start_year': int(report_start_day/365),
                                      'end_year': int(coord_df.at[site, 'simulation_duration']/365)})
            if coord_df.at[site, 'infectiousness_to_mosquitos']:
                analyzers.append(InfectiousnessByParDensAgeAnalyzer)
                analyzer_args.append({'expt_name': site,
                                      'sweep_variables': ['Run_Number', 'Site'],
                                      'start_year': int(report_start_day/365),
                                      'end_year': int(coord_df.at[site, 'simulation_duration']/365)})
            if coord_df.at[site, 'age_prevalence']:
                analyzers.append(MonthlySummaryReportAnalyzer)
                analyzer_args.append({'expt_name': site,
                                      'sweep_variables': ['Run_Number', 'Site'],
                                      'start_year': int(report_start_day / 365),
                                      'end_year': int(coord_df.at[site, 'simulation_duration']/365)})
        if coord_df.at[site, 'include_AnnualMalariaSummaryReport']:
            analyzers.append(AnnualSummaryReportAnalyzer)
            analyzer_args.append({'expt_name': site,
                                  'sweep_variables': ['Run_Number', 'Site']})
        if coord_df.at[site, 'include_MalariaPatientReport']:  # infection duration
            analyzers.append(PatientAnalyzer)
            analyzer_args.append({'expt_name': site,
                                  'start_report_day': report_start_day})

        analysis = PlatformAnalysis(platform=platform, experiment_ids=[id],
                                    analyzers=analyzers,
                                    analyzers_args=analyzer_args,
                                    analysis_name=site)

        analysis.analyze(check_status=True)
        wi = analysis.get_work_item()
        print(wi)
