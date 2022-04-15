# use the coordinator csv to determine which analyzers should be run for each site
import pandas as pd
import simulations.manifest as manifest

from idmtools.core.platform_factory import Platform
from idmtools.analysis.platform_anaylsis import PlatformAnalysis
from simulations.analyzers.SummaryReportAnalyzer import SummaryReportAnalyzer
from simulations.analyzers.ParDensAgeAnalyzer import ParDensAgeAnalyzer
from simulations.analyzers.PatientReportAnalyzer import PatientAnalyzer

exp_name_id = {
    'sugungum_1970': 'b7126585-30b6-ec11-a9f6-9440c9be2c51'
}
start_report_day = 365 * 25

if __name__ == "__main__":
    platform = Platform('CALCULON')
    coord_df = pd.read_csv(manifest.simulation_coordinator_path)
    coord_df = coord_df.set_index('site')

    for expt_name, id in exp_name_id.items():
        site = expt_name
        # determine the analyzers to run for each site
        analyzers = []
        analyzer_args = []
        if coord_df.at[site, 'include_MonthlyMalariaSummaryReport']:
            analyzers.append(ParDensAgeAnalyzer)
            analyzer_args.append({'expt_name': expt_name,
                                  'sweep_variables': ['Run_Number', 'Site'],
                                  'end_year': int(coord_df.at[site, 'simulation_duration']/365)})
        if coord_df.at[site, 'include_AnnualMalariaSummaryReport']:
            analyzers.append(SummaryReportAnalyzer)
            analyzer_args.append({'expt_name': expt_name,
                                  'sweep_variables': ['Run_Number', 'Site']})
        if coord_df.at[site, 'include_MalariaPatientReport']:
            analyzers.append(PatientAnalyzer)
            analyzer_args.append({'working_dir': 'outputs',
                                  'dir_name': expt_name,
                                  'start_report_day': start_report_day})

        analysis = PlatformAnalysis(platform=platform, experiment_ids=[id],
                                    analyzers=analyzers,
                                    analyzers_args=analyzer_args,
                                    analysis_name=expt_name)

        analysis.analyze(check_status=True)
        wi = analysis.get_work_item()
        print(wi)
