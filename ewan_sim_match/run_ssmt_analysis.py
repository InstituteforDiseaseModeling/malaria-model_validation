from general.analyzers.SummaryReportAnalyzer import SummaryReportAnalyzer
from idmtools.core.platform_factory import Platform
from idmtools.analysis.platform_anaylsis import PlatformAnalysis

# site = '9'
# analysis_name = "ssmt_model_validation_ewan_%s" % site
# analysis_name = "ssmt_ewan_age_inc_EIR_season_cm_sweep"

exp_name_id = {
    'ssmt_model_validation_ewan_1': '911e2edc-2b8a-ec11-a9f3-9440c9be2c51',
    'ssmt_model_validation_ewan_2': 'bba50600-2d8a-ec11-a9f3-9440c9be2c51',
    'ssmt_model_validation_ewan_3': '6df8920d-2d8a-ec11-a9f3-9440c9be2c51',
    'ssmt_model_validation_ewan_4': '84f8920d-2d8a-ec11-a9f3-9440c9be2c51',
    'ssmt_model_validation_ewan_5': 'fb4ee016-2d8a-ec11-a9f3-9440c9be2c51',
    'ssmt_model_validation_ewan_7': '7bd60a1f-2d8a-ec11-a9f3-9440c9be2c51',
    'ssmt_model_validation_ewan_8': '7ed60a1f-2d8a-ec11-a9f3-9440c9be2c51',
    'ssmt_model_validation_ewan_9': 'b2d62b27-2d8a-ec11-a9f3-9440c9be2c51',
    'ssmt_ewan_age_inc_EIR_season_cm_sweep': '6c1728b9-988b-ec11-a9f3-9440c9be2c51',
}



if __name__ == "__main__":
    platform = Platform('CALCULON')
    for analysis_name, id in exp_name_id.items():
        analysis = PlatformAnalysis(platform=platform, experiment_ids=[id],
                                    analyzers=[SummaryReportAnalyzer],
                                    analyzers_args=[],
                                    analysis_name=analysis_name)

        analysis.analyze(check_status=True)
        wi = analysis.get_work_item()
        print(wi)