from general.analyzers.SummaryReportAnalyzer import SummaryReportAnalyzer
from idmtools.core.platform_factory import Platform
from idmtools.analysis.platform_anaylsis import PlatformAnalysis

exp_name_id = {
    # 'ssmt_model_validation_ewan_1': '5e5950c0-1790-ec11-a9f3-9440c9be2c51',
    # 'ssmt_model_validation_ewan_2': '615950c0-1790-ec11-a9f3-9440c9be2c51',
    # 'ssmt_model_validation_ewan_3': '116f7acb-1790-ec11-a9f3-9440c9be2c51',
    # 'ssmt_model_validation_ewan_4': '146f7acb-1790-ec11-a9f3-9440c9be2c51',
    # 'ssmt_model_validation_ewan_5': '4a3b76d2-1790-ec11-a9f3-9440c9be2c51',
    # 'ssmt_model_validation_ewan_6': '5f3b76d2-1790-ec11-a9f3-9440c9be2c51',
    # 'ssmt_model_validation_ewan_7': '9807f9df-1790-ec11-a9f3-9440c9be2c51',
    # 'ssmt_model_validation_ewan_8': '9b07f9df-1790-ec11-a9f3-9440c9be2c51',
    # 'ssmt_model_validation_ewan_9': 'bb656ee7-1790-ec11-a9f3-9440c9be2c51',
    # 'ssmt_ewan_age_inc_EIR_season_cm_sweep': '3a730310-1690-ec11-a9f3-9440c9be2c51',

    # 'ssmt_model_validation_ewan_eirv2_1': 'a3ac395a-4b93-ec11-a9f3-9440c9be2c51',
    # 'ssmt_model_validation_ewan_eirv2_2': 'adacb279-4b93-ec11-a9f3-9440c9be2c51',
    # 'ssmt_model_validation_ewan_eirv2_3': '30080580-4b93-ec11-a9f3-9440c9be2c51',
    # 'ssmt_model_validation_ewan_eirv2_4': '77da5b90-4b93-ec11-a9f3-9440c9be2c51',
    # 'ssmt_model_validation_ewan_eirv2_5': '7669a19f-4b93-ec11-a9f3-9440c9be2c51',
    # 'ssmt_model_validation_ewan_eirv2_6': '7969a19f-4b93-ec11-a9f3-9440c9be2c51',
    # 'ssmt_model_validation_ewan_eirv2_7': '5ae999a7-4b93-ec11-a9f3-9440c9be2c51',
    # 'ssmt_model_validation_ewan_eirv2_8': '81e999a7-4b93-ec11-a9f3-9440c9be2c51',
    # 'ssmt_model_validation_ewan_eirv2_9': '2ae22baf-4b93-ec11-a9f3-9440c9be2c51',
    # 'ssmt_ewan_age_inc_EIR_season_cm_sweep': '17f93f8d-7e96-ec11-a9f5-9440c9be2c51',

    # 'ssmt_model_validation_ewan_eirv3_2': '032e0114-1f94-ec11-a9f3-9440c9be2c51',
    # 'ssmt_model_validation_ewan_eirv3_4': '1628f919-1f94-ec11-a9f3-9440c9be2c51',
    # 'ssmt_model_validation_ewan_eirv3_5': '1928f919-1f94-ec11-a9f3-9440c9be2c51',
    # 'ssmt_model_validation_ewan_eirv3_7': '85c17021-1f94-ec11-a9f3-9440c9be2c51',
    'ssmt_model_validation_ewan_eirv3_6': '9a2ac1b1-d998-ec11-a9f5-9440c9be2c51'
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
