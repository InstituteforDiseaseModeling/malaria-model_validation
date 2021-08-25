import json
import os
import pandas as pd
import numpy as np
from typing import Dict, Any, Union
from idmtools.entities.ianalyzer import IAnalyzer as BaseAnalyzer

import matplotlib as mpl
from idmtools.entities.iworkflow_item import IWorkflowItem
from idmtools.entities.simulation import Simulation

from logging import getLogger

from idmtools.analysis.analyze_manager import AnalyzeManager
from idmtools.core import ItemType
from idmtools.core.platform_factory import Platform

mpl.use('Agg')


class SummaryReportAnalyzer(BaseAnalyzer):

    def __init__(self, title='idm', tags=['Run_Number']):
        super().__init__(filenames=["output\\MalariaSummaryReport_Annual_Report.json"])
        self.tags = tags
        print(title)

    def initialize(self):
        """
        Initialize our Analyzer. At the moment, this just creates our output folder
        Returns:
        """
        if not os.path.exists(os.path.join(self.working_dir, "output")):
            os.mkdir(os.path.join(self.working_dir, "output"))

    def map(self, data: Dict[str, Any], item: Union[IWorkflowItem, Simulation]) -> Any:
        """
        Extracts the Statistical Population, Data channel from InsetChart.
        Called for Each WorkItem/Simulation.
        Args:
            data: Data mapping str to content of file
            item: Item to Extract Data from(Usually a Simulation)
        Returns:
        """
        datatemp = data[self.filenames[0]]

        age_bins = datatemp['Metadata']['Age Bins']
        prevalence = datatemp['DataByTimeAndAgeBins']['PfPR by Age Bin']
        prevalence = np.array(np.array([i for i in prevalence]))
        prevalence[prevalence == 0] = np.nan
        prevalence = np.nanmean(prevalence, axis=1)

        incidence = datatemp['DataByTimeAndAgeBins']['Annual Clinical Incidence by Age Bin']
        incidence = np.array(np.array([i for i in incidence]))
        incidence[incidence == 0] = np.nan
        incidence = np.nanmean(incidence, axis=1)

        df = pd.DataFrame(list(zip(age_bins, prevalence, incidence)),
                          columns=['Age', 'Prevalence', 'Incidence'])

        return df

    def reduce(self, all_data: Dict[Union[IWorkflowItem, Simulation], Any]) -> Any:
        """
        Create the Final Population JSON and Plot
        Args:
            all_data: Populate data from all the Simulations
        Returns:
            None
        """
        output_dir = os.path.join(self.working_dir, "output")
        df_final = pd.DataFrame()
        for s, v in all_data.items():
            dftemp = v[0].copy()
            for t in self.tags:
                dftemp[t] = [s.tags[t]]*len(v[0])
            df_final = pd.concat([df_final, dftemp])
        df_final.to_csv(os.path.join(output_dir, "summary_data_full.csv"))
        # df_incidence.to_csv(os.path.join(output_dir, "sporozoite_reduction_incidence_full.csv"))

        groupby_tags = self.tags
        groupby_tags.remove('Run_Number')
        df_summarized = df_final[['Age', 'Prevalence', 'Incidence']].apply(
             np.mean).reset_index()
        df_summarized_std = df_final[['Prevalence', 'Incidence']].apply(
            np.std)
        for c in ['Prevalence', 'Incidence']:
            df_summarized[c + '_std'] = list(df_summarized_std[c])

        df_summarized.to_csv(os.path.join(output_dir, "summary_data_final.csv"))


if __name__ == '__main__':

    # Set the platform where you want to run your analysis
    # In this case we are running in BELEGOST since the Work Item we are analyzing was run on COMPS
    logger = getLogger()
    with Platform('CALCULON') as platform:

        # Initialize the analyser class with the path of the output csv file
        analyzers = [SummaryReportAnalyzer()]

        # Set the experiment id you want to analyze
        experiment_id = 'ab7487a1-6405-ec11-a9ed-b88303911bc1'

        # Specify the id Type, in this case an Experiment on COMPS
        manager = AnalyzeManager(partial_analyze_ok=True, ids=[(experiment_id, ItemType.EXPERIMENT)],
                                 analyzers=analyzers)
        manager.analyze()
