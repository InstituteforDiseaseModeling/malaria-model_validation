# create input csv with simulation days repeated across years
import os
import pandas as pd
import simulations.manifest as manifest

target_filename = 'bfa_2007.csv'  # save output in a csv with this name (saved in the survey_days directory)
day_of_year_values = [21, 195, 264]  # days of the year when the survey should occur
num_years = 66  # number of years over which the surveys should repeat

all_days = [ii + (365 * xx) for xx in range(num_years) for ii in day_of_year_values]

df = pd.DataFrame(data=all_days, columns=['days'])
df.to_csv(os.path.join(manifest.input_files_path, 'survey_days', target_filename))

