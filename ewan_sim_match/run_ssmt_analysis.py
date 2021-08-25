from ewan_sim_match.analyzers.SummaryReportAnalyzer import SummaryReportAnalyzer
from idmtools.core.platform_factory import Platform
from idmtools.analysis.platform_anaylsis import PlatformAnalysis

if __name__ == "__main__":
    platform = Platform('CALCULON')
    analysis = PlatformAnalysis(platform=platform, experiment_ids=["ab7487a1-6405-ec11-a9ed-b88303911bc1"],
                                analyzers=[SummaryReportAnalyzer], # VectorGeneticsAnalyzer],
                                analyzers_args=[],
                                analysis_name="Model Validation")

    analysis.analyze(check_status=True)
    wi = analysis.get_work_item()
    print(wi)