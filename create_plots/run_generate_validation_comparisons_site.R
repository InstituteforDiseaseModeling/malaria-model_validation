# run_generate_validation_comparisons_site.R
# April 2022

# This is the main run script that sets parameters and calls the core function to create validation comparisons for each validation relationship. All site-specific simulation results are plotted and compared against their associated reference dataset.


library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)


####################################################################################################
##################### setup: filepaths, info on simulations run, source scripts ####################
####################################################################################################
# Note: these first three filephaths should be grabbed from the manifest once we move the scripts to python
simulation_output_filepath = "/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/simulation_output/s220506"
benchmark_simulation_filepath = "/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/EMOD_validation_recalibration/simulation_output/s220506"
validation_repo_filepath = "/Users/moniqueam/Documents/malaria-model_validation/"


simulation_coordinator_path = paste0(validation_repo_filepath, "/simulation_inputs/simulation_coordinator.csv")
base_script_plot_filepath = paste0(validation_repo_filepath, "/create_plots")
base_reference_filepath = paste0(validation_repo_filepath, "/reference_datasets")
# plot_output_filepath = paste0(simulation_output_filepath, "/_plots")
plot_output_filepath = paste0(validation_repo_filepath, "/report/_plots")


source(file.path(base_script_plot_filepath, 'helpers_coordinate_each_relationship.R'))
source(file.path(base_script_plot_filepath, 'helpers_reformat_sim_ref_dfs.R'))
source(file.path(base_script_plot_filepath, 'helpers_likelihood_and_metrics.R'))
source(file.path(base_script_plot_filepath, 'helpers_plot_ref_sim_comparisons.R'))


####################################################################################################
############################ read in data and create plots #########################################
####################################################################################################
coord_csv = read.csv(simulation_coordinator_path)
if(!dir.exists(plot_output_filepath)) dir.create(plot_output_filepath)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#                         age - incidence                         #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
generate_age_incidence_outputs(coord_csv, simulation_output_filepath, base_reference_filepath, plot_output_filepath, benchmark_simulation_filepath=benchmark_simulation_filepath)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#                         age - prevalence                        #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
generate_age_prevalence_outputs(coord_csv, simulation_output_filepath, base_reference_filepath, plot_output_filepath, benchmark_simulation_filepath=benchmark_simulation_filepath)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#                      age - parasite density                     #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
generate_parasite_density_outputs(coord_csv, simulation_output_filepath, base_reference_filepath, plot_output_filepath, benchmark_simulation_filepath=benchmark_simulation_filepath)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#                   infectiousness to vectors                        #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
generate_infectiousness_outputs(coord_csv, simulation_output_filepath, base_reference_filepath, plot_output_filepath, benchmark_simulation_filepath=benchmark_simulation_filepath)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#                    age - infection duration                     #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# set positive threshold density for sampled parasites in simulation output (to match PCR threshold in reference)
pos_thresh_dens = 0.5  # (39=min(ref_df$DENSITY[ref_df$DENSITY>0], na.rm=TRUE) - 1)
# specify binning for duration of infection
duration_bins=c(seq(0,350,50), 500)
generate_age_infection_duration_outputs(coord_csv, simulation_output_filepath, base_reference_filepath, plot_output_filepath, pos_thresh_dens, duration_bins, benchmark_simulation_filepath=benchmark_simulation_filepath)




