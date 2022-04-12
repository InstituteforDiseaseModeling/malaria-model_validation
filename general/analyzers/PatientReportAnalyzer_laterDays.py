import os
import pandas as pd
import numpy as np
from idmtools.entities import IAnalyzer
from idmtools.entities.simulation import Simulation

from logging import getLogger
logger = getLogger()


class PatientAnalyzer(IAnalyzer):

    def __init__(self, dir_name, working_dir='.', start_report_day=0):
        super(PatientAnalyzer, self).__init__(working_dir=working_dir,
                                              filenames=['output/MalariaPatientReport.json']
                                              )

        self.dir_name = dir_name
        self.sweep_variables = ['Run_Number', 'x_Temp_LH_values', 'Site']
        self.channels = ['true_gametocytes', 'true_asexual_parasites', 'temps']
        self.fields = ['id', 'initial_age'] + self.channels
        self.output_fname = os.path.join(self.working_dir, self.dir_name, "patient_reports.csv")
        self.start_report_day = start_report_day

        # make sure output folder exists
        os.makedirs(os.path.join(self.working_dir, self.dir_name), exist_ok=True)

    def map(self, data, simulation: Simulation):
        patients = data[self.filenames[0]]["patient_array"]
        ntsteps = data[self.filenames[0]]["ntsteps"]
        npatients = len(patients)
        logger.info('Parsing channels %s for %d patients and %d time steps',
             self.channels, npatients, ntsteps-self.start_report_day)

        # fill in long-format dataframe with patient information for each simulation day. columns are:
        #    0) Timestep of simulation
        #    1) Patient ID
        #    2) Patient age
        #    3-X) Channel values for that patient at that timestep
        patient_info_temp = np.empty((ntsteps-self.start_report_day, 3+len(self.channels)), dtype='float')
        patient_info_temp[:, 0] = np.arange(start=self.start_report_day, stop=ntsteps, step=1)  # day of simulation

        birthdays = [p['birthday'] for p in patients]
        initial_ages = [p['initial_age'] for p in patients]
        patient_df_list = []
        for i, p in enumerate(patients):
            patient_info = patient_info_temp.copy()
            patient_info[:, 1] = np.repeat(p['id'], ntsteps-self.start_report_day)  # patient ID
            pos_first_data = max(birthdays[i]-self.start_report_day, 0)  # get row when patient's data starts (corresponds to day of simulation when patient appears)
            # calculate age of patient on each day of simulation
            patient_info[pos_first_data:, 2] = np.arange(start=(initial_ages[i]+max(self.start_report_day - max(birthdays[i], 0), 0)),
                                                         stop=(initial_ages[i] + ntsteps - self.start_report_day + pos_first_data))
            # iterate through channels, adding to data frame
            for c, channel in enumerate(self.channels):
                patient_info[pos_first_data:, 3+c] = p[channel]
            patient_df_list.append(pd.DataFrame(data=patient_info, columns=['simday', 'id', 'age'] + self.channels))
        patient_df = pd.concat(patient_df_list)

        # remove any days before the start_report_day
        # patient_df['simday'] = patient_df['simday'] - self.start_report_day
        # patient_df = patient_df.loc[patient_df['simday'] >= 0]

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                patient_df[sweep_var] = simulation.tags[sweep_var]

        return patient_df

    def reduce(self, all_data):
        data_sets_per_experiment = {}

        for simulation, associated_data in all_data.items():
            experiment_name = simulation.experiment.name
            if experiment_name not in data_sets_per_experiment:
                data_sets_per_experiment[experiment_name] = []

            data_sets_per_experiment[experiment_name].append(associated_data)

        for experiment_name, data_sets in data_sets_per_experiment.items():
            d = pd.concat(data_sets).reset_index(drop=True)
            d.to_csv(self.output_fname, index=False)


if __name__ == '__main__':
    from idmtools.analysis.analyze_manager import AnalyzeManager
    from idmtools.core import ItemType
    from idmtools.core.platform_factory import Platform
    from general.manifest import analyzed_ouptut_path

    # Set the platform where you want to run your analysis
    # In this case we are running in BELEGOST since the Work Item we are analyzing was run on COMPS
    with Platform('Calculon') as platform:

        # Set the experiment you want to analyze
        experiments = {'test_infectDuration_Navrongo2': 'aaa3fc9f-787a-ec11-a9f3-9440c9be2c51'}

        for dirname, exp_id in experiments.items():
            # Initialize the analyser class with the path of the output csv file
            analyzers = [PatientAnalyzer(working_dir=analyzed_ouptut_path, dir_name=dirname, start_report_day=365*85)]

            # Create AnalyzerManager with required parameters
            manager = AnalyzeManager(ids=[(exp_id, ItemType.EXPERIMENT)],
                                     analyzers=analyzers)
            manager.analyze()


