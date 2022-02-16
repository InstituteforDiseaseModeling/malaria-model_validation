import os
import sys
from idmtools.analysis.analyze_manager import AnalyzeManager
from idmtools.core import ItemType
from idmtools.core.platform_factory import Platform
from idmtools.assets import AssetCollection
from idmtools.analysis.platform_anaylsis import PlatformAnalysis
from infection_duration_navrongo.manifest_asset import analyzed_ouptut_path
from general.analyzers.PatientReportAnalyzer import PatientAnalyzer

use_ssmt = True

if __name__ == "__main__":
    # Set the platform where you want to run your analysis
    # In this case we are running in BELEGOST since the Work Item we are analyzing was run on COMPS
    with Platform('Calculon') as platform:

        # Set the experiment you want to analyze
        experiments = {
            'model_validation_duration_Navrongo_30y_xEIR_v2': '14dd55a8-1e83-ec11-a9f3-9440c9be2c51'}
        start_report_day = 365*25

        if use_ssmt:
            for dirname, exp_id in experiments.items():
                # Initialize the SSMT PlatformAnalysis class with all parameters
                # first add files needed for prevalence_analyzer.py to assetcollection in ssmt
                ac = AssetCollection()

                analysis = PlatformAnalysis(platform=platform,
                                            experiment_ids=[exp_id],
                                            analyzers=[PatientAnalyzer],
                                            analyzers_args=[{'working_dir': 'outputs',
                                                             'dir_name': dirname,
                                                             'start_report_day': start_report_day}],
                                            analysis_name=os.path.split(sys.argv[0])[1],
                                            tags={'WorkItem type': 'Docker'},
                                            asset_files=ac,
                                            wait_till_done=True)
                # Run analysis on COMPS
                analysis.analyze(check_status=True)
        else:
            for dirname, exp_id in experiments.items():
                # Initialize the analyser class with the path of the output csv file
                analyzers = [PatientAnalyzer(working_dir=analyzed_ouptut_path, dir_name=dirname, start_report_day=start_report_day)]

                # Create AnalyzerManager with required parameters
                manager = AnalyzeManager(ids=[(exp_id, ItemType.EXPERIMENT)],
                                         analyzers=analyzers)
                manager.analyze()
