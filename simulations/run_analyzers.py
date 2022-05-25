import argparse
import simulations.params as params
import simulations.manifest as manifest
from simulations.helpers import get_comps_id_filename, load_coordinator_df, get_suite_id

from idmtools.core.platform_factory import Platform
from idmtools.analysis.platform_anaylsis import PlatformAnalysis
from simulations.analyzers.AnnualSummaryReportAnalyzer import AnnualSummaryReportAnalyzer
from simulations.analyzers.ParDensAgeAnalyzer import ParDensAgeAnalyzer
from simulations.analyzers.InfectiousnessByParDensAgeAnalyzer import InfectiousnessByParDensAgeAnalyzer
from simulations.analyzers.PatientReportAnalyzer import PatientAnalyzer
from simulations.analyzers.MonthlySummaryReportAnalyzer import MonthlySummaryReportAnalyzer
from simulations.wait_for_experiment import check_experiment


def run_analyzers(site, characteristic=False):
    platform = Platform(manifest.platform_name)
    comps_id_file = get_comps_id_filename(site=site)
    with open(comps_id_file, 'r') as id_file:
        exp_id = id_file.readline()
    # Wait for experiment to be done
    if check_experiment(site, platform):
        coord_df = load_coordinator_df(characteristic=characteristic, set_index=True)

        # for expt_name, id in exp_name_id.items():
        # site = expt_name.replace('validation_', '')
        report_start_day = int(coord_df.at[site, 'report_start_day'])
        # determine the analyzers to run for each site
        analyzers = []
        analyzer_args = []
        if coord_df.at[site, 'include_MonthlyMalariaSummaryReport']:
            if coord_df.at[site, 'age_parasite_density']:
                analyzers.append(ParDensAgeAnalyzer)
                analyzer_args.append({'expt_name': site,
                                      'sweep_variables': ['Run_Number', 'Site'],
                                      'start_year': int(report_start_day / 365),
                                      'end_year': int(coord_df.at[site, 'simulation_duration'] / 365)})
            if coord_df.at[site, 'infectiousness_to_mosquitos']:
                analyzers.append(InfectiousnessByParDensAgeAnalyzer)
                analyzer_args.append({'expt_name': site,
                                      'sweep_variables': ['Run_Number', 'Site'],
                                      'start_year': int(report_start_day / 365),
                                      'end_year': int(coord_df.at[site, 'simulation_duration'] / 365)})
            if coord_df.at[site, 'age_prevalence']:
                analyzers.append(MonthlySummaryReportAnalyzer)
                analyzer_args.append({'expt_name': site,
                                      'sweep_variables': ['Run_Number', 'Site'],
                                      'start_year': int(report_start_day / 365),
                                      'end_year': int(coord_df.at[site, 'simulation_duration'] / 365)})
        if coord_df.at[site, 'include_AnnualMalariaSummaryReport']:
            analyzers.append(AnnualSummaryReportAnalyzer)
            analyzer_args.append({'expt_name': site,
                                  'sweep_variables': ['Run_Number', 'Site']})
        if coord_df.at[site, 'include_MalariaPatientReport']:  # infection duration
            analyzers.append(PatientAnalyzer)
            analyzer_args.append({'expt_name': site,
                                  'start_report_day': report_start_day})

        analysis = PlatformAnalysis(platform=platform, experiment_ids=[exp_id],
                                    analyzers=analyzers,
                                    analyzers_args=analyzer_args,
                                    analysis_name=site)

        suite_id = get_suite_id()
        analysis.tags = {'Suite': suite_id}
        analysis.analyze(check_status=True)

        wi = analysis.get_work_item()
        analyzers_id_file = get_comps_id_filename(site=site, level=2)

        if wi.succeeded:
            print(f"Analyzer work item {wi.uid} succeeded.\n")
            with open(analyzers_id_file, 'w') as id_file:
                id_file.write(wi.uid.hex)
        else:
            print(f"Analyzer work item {wi.uid} failed.")

        return wi.succeeded, wi.uid
    else:
        return False, exp_id


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process site name')
    parser.add_argument('--site', '-s', type=str, help='site name',
                        default=params.sites[0])  # not sure if we want to make this required argument
    args = parser.parse_args()
    run_analyzers(args.site)
