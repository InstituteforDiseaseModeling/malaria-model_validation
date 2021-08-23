from analyzers.InsetChartAnalyzer import InsetChartAnalyzer
from analyzers.VectorGeneticsAnalyzer import VectorGeneticsAnalyzer
from idmtools.core.platform_factory import Platform
from idmtools.analysis.platform_anaylsis import PlatformAnalysis

if __name__ == "__main__":
    platform = Platform('CALCULON')
    analysis = PlatformAnalysis(platform=platform, experiment_ids=["483e76e6-66eb-eb11-a9ed-b88303911bc1"],
                                analyzers=[VectorGeneticsAnalyzer, InsetChartAnalyzer], # VectorGeneticsAnalyzer],
                                analyzers_args=[],
                                analysis_name="Sporozoite Delay")

    analysis.analyze(check_status=True)
    wi = analysis.get_work_item()
    print(wi)