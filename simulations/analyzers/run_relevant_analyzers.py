# use the coordinator csv to determine which analyzers should be run for each site
import pandas as pd
import simulations.manifest as manifest

from idmtools.core.platform_factory import Platform
from idmtools.analysis.platform_anaylsis import PlatformAnalysis
from simulations.analyzers.SummaryReportAnalyzer import SummaryReportAnalyzer
from simulations.analyzers.ParDensAgeAnalyzer import ParDensAgeAnalyzer
from simulations.analyzers.PatientReportAnalyzer import PatientAnalyzer

exp_name_id = {
    'validation_chonyi_1999': '251af5cc-5bbf-ec11-a9f6-9440c9be2c51',
    'validation_ngerenya_1999': 'cd3ba6e2-5bbf-ec11-a9f6-9440c9be2c51',
    'validation_dielmo_1990': 'd843d0e8-5bbf-ec11-a9f6-9440c9be2c51',
    'validation_ndiop_1993': 'b248b3f0-5bbf-ec11-a9f6-9440c9be2c51',
    'validation_matsari_1970': '26ebf0ce-70bf-ec11-a9f6-9440c9be2c51',
    'validation_rafin_marke_1970': '83c662c8-70bf-ec11-a9f6-9440c9be2c51',
    'validation_sugungum_1970': 'f20c83bb-70bf-ec11-a9f6-9440c9be2c51',
    'validation_navrongo_2000_dur': 'ce0827ae-70bf-ec11-a9f6-9440c9be2c51',
    'validation_laye_2007': 'cce16da7-70bf-ec11-a9f6-9440c9be2c51',
    'validation_dapelogo_2007': '464baba0-70bf-ec11-a9f6-9440c9be2c51',
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
            analyzers.append(ParDensAgeAnalyzer)
            analyzer_args.append({'expt_name': site,
                                  'sweep_variables': ['Run_Number', 'Site'],
                                  'end_year': int(coord_df.at[site, 'simulation_duration']/365)})
        if coord_df.at[site, 'include_AnnualMalariaSummaryReport']:
            analyzers.append(SummaryReportAnalyzer)
            analyzer_args.append({'expt_name': site,
                                  'sweep_variables': ['Run_Number', 'Site']})
        if coord_df.at[site, 'include_MalariaPatientReport']:
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
