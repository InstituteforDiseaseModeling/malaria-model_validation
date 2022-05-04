import os
import warnings
import pandas as pd
import numpy as np
from typing import Dict, Any, Union
from logging import getLogger
from idmtools.core.platform_factory import Platform
from idmtools.analysis.analyze_manager import AnalyzeManager
from idmtools.core import ItemType
from idmtools.entities.ianalyzer import IAnalyzer as BaseAnalyzer


class InfectiousnessByParDensAgeAnalyzer(BaseAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir=".", start_year=0, end_year=65):
        super(InfectiousnessByParDensAgeAnalyzer, self).__init__(working_dir=working_dir,
                                                 filenames=["output/MalariaSummaryReport_Infectiousness_Monthly_Report_%d.json" % x
                                                            for x in range(start_year, end_year)])
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year
        self.working_dir = working_dir

    def initialize(self):
        """
        Initialize our Analyzer. At the moment, this just creates our output folder
        Returns:
        """
        if not os.path.exists(os.path.join(self.working_dir, self.expt_name)):
            os.mkdir(os.path.join(self.working_dir, self.expt_name))


    def map(self, data, simulation):

        agebins = data[self.filenames[0]]['Metadata']['Age Bins']
        gam_bins = data[self.filenames[0]]['Metadata']['Gametocytemia Bins']
        frac_infected_bins = data[self.filenames[0]]['Metadata']['Infectiousness Bins']

        adf = pd.DataFrame()
        for year, fname in zip(range(self.start_year, self.end_year), self.filenames):  # , 2020
            for age in list(range(0, len(agebins))):
                for dens in list(range(0, len(gam_bins))):
                    for infect in list(range(0, len(frac_infected_bins))):

                        d = data[fname]['DataByTimeAndInfectiousnessBinsAndPfPRBinsAndAgeBins']['Smeared Infectiousness by smeared Gametocytemia and Age Bin'][:12]
                        infect_bin_frac = [x[infect][dens][age] for x in d]
                        d = data[fname]['DataByTimeAndAgeBins']['Average Population by Age Bin'][:12]
                        pop = [x[age] for x in d]
                        simdata = pd.DataFrame({'month': range(1, 13),
                                                'infectiousness_bin_freq': infect_bin_frac,
                                                'Pop': pop,
                                                })
                        simdata['agebin'] = agebins[age]
                        simdata['densitybin'] = gam_bins[dens]
                        simdata['infectiousness_bin'] = frac_infected_bins[infect]
                        simdata['year'] = year
                        adf = pd.concat([adf, simdata])

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                adf[sweep_var] = simulation.tags[sweep_var]
        return adf


    def reduce(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv((os.path.join(self.working_dir, self.expt_name, 'infectiousness_by_age_density_month.csv')),
                   index=False)


if __name__ == '__main__':

    # Set the experiment id you want to analyze
    experiment_id = 'b7126585-30b6-ec11-a9f6-9440c9be2c51'
    end_year = 65

    # Set the platform where you want to run your analysis
    # In this case we are running in BELEGOST since the Work Item we are analyzing was run on COMPS
    logger = getLogger()
    with Platform('CALCULON') as platform:

        # Initialize the analyser class with the path of the output csv file
        analyzers = [InfectiousnessByParDensAgeAnalyzer(expt_name='Sugungum_1970',
                                        sweep_variables=['Run_Number', 'Site'],
                                        end_year=end_year)]

        # Specify the id Type, in this case an Experiment on COMPS
        manager = AnalyzeManager(partial_analyze_ok=True, ids=[(experiment_id, ItemType.EXPERIMENT)],
                                 analyzers=analyzers)
        manager.analyze()
