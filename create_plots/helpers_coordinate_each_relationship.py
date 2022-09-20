# helpers_coordinate_each_relationship.py

# For each validation relationship, this script contains a function that coordinates:
#     1) the reformatting and alignment of simulation and reference datasets into data frames that are used downstream,
#     2) the quantitative comparisons between reference and new simulation results and between benchmark and new
#        simulation results, and
#     3) the generation of plotting outputs.
# The main function for each relationship saves all plots and csvs to the output directory for the report-generating
# script to use.


import pandas as pd
import os

from create_plots.helpers_reformat_sim_ref_dfs import prepare_inc_df, prepare_prev_df, prepare_dens_df, \
    prepare_infect_df, get_available_sites_for_relationship, get_sim_survey
from create_plots.helpers_plot_ref_sim_comparisons import plot_inc_ref_sim_comparison, plot_prev_ref_sim_comparison, \
    compare_benchmark, plot_par_dens_ref_sim_comparison, plot_infectiousness_ref_sim_comparison, \
    plot_infection_duration_dist, plot_infection_duration_dist_by_age, create_barplot_frac_comparison
from create_plots.helpers_likelihood_and_metrics import calc_mean_rel_diff, calc_mean_rel_slope_diff, \
    get_prev_loglikelihood, get_dens_loglikelihood, corr_ref_sim_points, corr_ref_deriv_sim_points, add_to_summary_table


# todo: create one base generate output function for all
# Incidence by age
def generate_age_incidence_outputs(coord_csv, simulation_output_filepath, base_reference_filepath, plot_output_filepath,
                                   benchmark_simulation_filepath=None):
    """
    From simulation output and matched reference data, create plots and quantitative comparisons for all sites
    associated with the incidence-by-age validation relationship.
    Args:
        coord_csv (): A dataframe detailing the sites simulated for each validation relationship and the corresponding
                      reference dataset
        simulation_output_filepath (): The filepath where simulation output is located
        base_reference_filepath (): The filepath where reference datasets are located
        plot_output_filepath (): The filepath to the directory where plots should be created
        benchmark_simulation_filepath (): The filepath where benchmark simulation output is located. If None, no
                                          comparisons are made against benchmark simulations

    Returns:

    """
    # get formatted dataframe with reference and simulation incidence data from all relevant sites
    combined_df = prepare_inc_df(coord_csv, simulation_output_filepath, base_reference_filepath,
                                 benchmark_simulation_filepath)

    # create plots comparing reference and simulation outputs
    gg_plot = plot_inc_ref_sim_comparison(combined_df)
    gg_plot.save(filename=os.path.join(plot_output_filepath, 'site_compare_incidence_age.png'))
                 # height=2 * math.ceil(len(combined_df['Site'].unique()) / 4), width=7.5, units='in')

    # additional quantitative comparisons and metrics between simulation and reference data
    # correlations between new simulation and reference dataset values
    correlation_output = corr_ref_sim_points(combined_df)
    correlation_output[0].save(filename=os.path.join(plot_output_filepath, 'corr_ref_sim_points_incidence_age.png'), height=9, width=8,
                 units='in')
    correlation_df = correlation_output[1]
    slope_correlation_output = corr_ref_deriv_sim_points(combined_df)
    slope_correlation_df = slope_correlation_output[1]
    combined_df_with_slopes = slope_correlation_output[2]
    # todo: no equivalent for ggarrange in Python
    # correlation_plots = ggarrange(correlation_output[0], slope_correlation_output[0], nrow=1, ncol=2,
    #                               common.legend = TRUE)  # , legend.grob=get_legend(correlation_output[[1]], position = 'bottom'))
    # correlation_plots.save(filename=os.path.join(plot_output_filepath, 'scatter_regression_incidence_age.png'),
    #                        height=4.5, width=8, units='in')
    correlation_output[0].save(filename=os.path.join(plot_output_filepath, 'scatter_regression_incidence_age_correlation.png'),
                               height=4.5, width=8, units='in')
    slope_correlation_output[0].save(
        filename=os.path.join(plot_output_filepath, 'scatter_regression_incidence_age_slope_correlation.png'),
        height=4.5, width=8, units='in')

    # metrics comparing simulation to reference VALUE
    mean_diff_df = calc_mean_rel_diff(combined_df)
    correlation_df['corr_slope'] = correlation_df['slope']
    correlation_df['corr_r_squared'] = correlation_df['r.squared']
    # metrics comparing simulation to reference SLOPE
    mean_slope_diff_df = calc_mean_rel_slope_diff(combined_df=combined_df_with_slopes)
    slope_correlation_df['derivative_corr_slope'] = slope_correlation_df['slope']
    slope_correlation_df['derivative_corr_r_squared'] = slope_correlation_df['r.squared']
    # combine metrics and save csv
    quantitative_comparison_df = pd.merge(mean_diff_df, correlation_df[['Site', 'corr_slope', 'corr_r_squared']],
                                          how="outer")
    # quantitative_comparison_slope_df = merge(mean_slope_diff_df, slope_correlation_df[,c('Site', 'derivative_corr_r_squared')], all=TRUE)
    # quantitative_comparison_df = merge(quantitative_comparison_df, quantitative_comparison_slope_df, all=TRUE)
    quantitative_comparison_df = pd.merge(quantitative_comparison_df, mean_slope_diff_df, how="outer")
    quantitative_comparison_df.to_csv(os.path.join(plot_output_filepath, 'comparison_metric_table_incidence_age.csv'),
                                      index=False)

    # compare simulation and benchmark simulation results
    if 'benchmark' in combined_df.columns:
        compare_benchmarks_output = compare_benchmark(combined_df)
        compare_benchmarks_output.save(filename=os.path.join(plot_output_filepath,
                                                             'scatter_benchmark_incidence_age.png'),
                                       height=4.5, width=8, units='in')
        mean_diff_df_bench = calc_mean_rel_diff(combined_df, sim_colname='benchmark')

        add_to_summary_table(combined_df=combined_df, plot_output_filepath=plot_output_filepath,
                             validation_relationship_name='age_incidence')


# Prevalence by age
def generate_age_prevalence_outputs(coord_csv, simulation_output_filepath, base_reference_filepath,
                                    plot_output_filepath, benchmark_simulation_filepath=None):
    """
    From simulation output and matched reference data, create plots and quantitative comparisons for all sites
    associated with the prevalence-by-age validation relationship.
    Args:
        coord_csv (): A dataframe detailing the sites simulated for each validation relationship and the corresponding
                      reference dataset
        simulation_output_filepath (): The filepath where simulation output is located
        base_reference_filepath (): The filepath where reference datasets are located
        plot_output_filepath (): The filepath to the directory where plots should be created
        benchmark_simulation_filepath (): The filepath where benchmark simulation output is located. If None, no
                                          comparisons are made against benchmark simulations

    Returns:

    """

    combined_df = prepare_prev_df(coord_csv, simulation_output_filepath, base_reference_filepath,
                                  benchmark_simulation_filepath)

    # create plots comparing reference and simulation outputs
    gg_plot = plot_prev_ref_sim_comparison(combined_df)
    gg_plot.save(filename=os.path.join(plot_output_filepath, 'site_compare_prevalence_age.png'), height=9, width=8,
                 units='in')

    # additional quantitative comparisons and metrics
    # correlations between new simulation and reference dataset values
    correlation_output = corr_ref_sim_points(combined_df)
    correlation_output[0].save(filename=os.path.join(plot_output_filepath, 'corr_ref_sim_points_prevalence_age.png'), height=9, width=8,
                                units='in')
    correlation_df = correlation_output[1]
    slope_correlation_output = corr_ref_deriv_sim_points(combined_df)
    slope_correlation_df = slope_correlation_output[1]
    combined_df_with_slopes = slope_correlation_output[2]
    # todo: no equivalent for ggarrange in Python
    # correlation_plots = ggarrange(correlation_output[0], slope_correlation_output[0], nrow=1, ncol=2,
    #                               common.legend = TRUE)  # , legend.grob=get_legend(correlation_output[[1]], position = 'bottom'))
    # correlation_plots.save(filename=os.path.join(plot_output_filepath, 'scatter_regression_prevalence_age.png'),
    #                        height=4.5, width=8, units='in')

    # metrics comparing simulation to reference VALUE
    mean_diff_df = calc_mean_rel_diff(combined_df)
    correlation_df['corr_slope'] = correlation_df['slope']
    correlation_df['corr_r_squared'] = correlation_df['r.squared']
    # metrics comparing simulation to reference SLOPE
    mean_slope_diff_df = calc_mean_rel_slope_diff(combined_df=combined_df_with_slopes)
    slope_correlation_df['derivative_corr_slope'] = slope_correlation_df['slope']
    slope_correlation_df['derivative_corr_r_squared'] = slope_correlation_df['r.squared']
    # combine metrics and save csv
    quantitative_comparison_df = pd.merge(mean_diff_df, correlation_df[['Site', 'corr_slope', 'corr_r_squared']],
                                          how="outer")
    # quantitative_comparison_slope_df = merge(mean_slope_diff_df, slope_correlation_df[,c('Site', 'derivative_corr_r_squared')], all=TRUE)
    # quantitative_comparison_df = merge(quantitative_comparison_df, quantitative_comparison_slope_df, all=TRUE)
    quantitative_comparison_df = pd.merge(quantitative_comparison_df, mean_slope_diff_df, how="outer")
    quantitative_comparison_df.to_csv(os.path.join(plot_output_filepath, 'comparison_metric_table_prevalence_age.csv'),
                                      index=False)

    # compare simulation and benchmark simulation results
    if 'benchmark' in combined_df.columns:
        compare_benchmarks_output = compare_benchmark(combined_df)
        compare_benchmarks_output.save(filename=os.path.join(plot_output_filepath,
                                                             'scatter_benchmark_prevalence_age.png'),
                                       height=4.5, width=8, units='in')
        # add likelihood component
        new_sim_loglik = get_prev_loglikelihood(combined_df, sim_column='simulation')
        bench_sim_loglik = get_prev_loglikelihood(combined_df, sim_column='benchmark')
        new_sim_loglik.rename(columns={"loglikelihood": "loglikelihood_new_sim"}, inplace=True)
        bench_sim_loglik.rename(columns={"loglikelihood": "loglikelihood_benchmark_sim"}, inplace=True)
        loglikelihood_comparison = pd.merge(new_sim_loglik, bench_sim_loglik)
        loglikelihood_comparison.to_csv(os.path.join(plot_output_filepath, 'loglikelihood_prevalence_age.csv'),
                                        index=False)
        add_to_summary_table(combined_df=combined_df, plot_output_filepath=plot_output_filepath,
                             validation_relationship_name='age_prevalence')


# Parasite density by age
def generate_parasite_density_outputs(coord_csv, simulation_output_filepath, base_reference_filepath,
                                      plot_output_filepath, benchmark_simulation_filepath=None):
    """
    From simulation output and matched reference data, create plots and quantitative comparisons for all sites
    associated with the parasite density-by-age validation relationship.
    Args:
        coord_csv (): A dataframe detailing the sites simulated for each validation relationship and the corresponding
                      reference dataset
        simulation_output_filepath (): The filepath where simulation output is located
        base_reference_filepath (): The filepath where reference datasets are located
        plot_output_filepath (): The filepath to the directory where plots should be created
        benchmark_simulation_filepath (): The filepath where benchmark simulation output is located. If None, no
                                          comparisons are made against benchmark simulations
    Returns:

    """

    # get formatted dataframe with reference and simulation prevalence data from all relevant sites
    combined_dfs = prepare_dens_df(coord_csv, simulation_output_filepath, base_reference_filepath,
                                   benchmark_simulation_filepath)
    combined_df_asex = combined_dfs[0]
    combined_df_gamet = combined_dfs[1]

    # todo: combine these 2 plotting block in to one function
    # asexual parasite density
    plot_output = plot_par_dens_ref_sim_comparison(combined_df=combined_df_asex)
    gg_barplot = plot_output[0]
    line_plot_list = plot_output[1]
    all_sites = plot_output[2]
    gg_barplot.save(filename=os.path.join(plot_output_filepath, 'site_compare_barplot_asex_dens_age.png'), width=10,
                    height=20, units='in')
    for ss in range(len(all_sites)):
        line_plot_list[ss].save(filename=os.path.join(plot_output_filepath,
                                                      'site_compare_asex_dens_age_' + all_sites[ss] + '.png'),
                                width=8, height=6, units='in')

    # gametocyte density
    plot_output = plot_par_dens_ref_sim_comparison(combined_df = combined_df_gamet)
    gg_barplot = plot_output[0]
    line_plot_list = plot_output[1]
    all_sites = plot_output[2]
    gg_barplot.save(filename=os.path.join(plot_output_filepath, 'site_compare_barplot_gamet_dens_age.png'), width=5,
                    height=15, units='in')
    for ss in range(len(all_sites)):
        line_plot_list[ss].save(filename=os.path.join(plot_output_filepath,
                                                      'site_compare_gamet_dens_age_' + all_sites[ss] + '.png'),
                                width=8, height=6, units='in')

    # compare simulation and benchmark simulation results
    if 'benchmark' in combined_df_asex.columns:
        compare_benchmarks_output = compare_benchmark(combined_df_asex)
        compare_benchmarks_output.save(filename=os.path.join(plot_output_filepath,
                                                             'scatter_benchmark_asex_dens.png'),
                                       height=4.5, width=8, units='in')
        compare_benchmarks_output = compare_benchmark(combined_df_gamet)
        compare_benchmarks_output.save(filename=os.path.join(plot_output_filepath,
                                                             'scatter_benchmark_gamet_dens.png'),
                                       height=4.5, width=8, units='in')
        # add likelihood component
        loglik_df_asex = get_dens_loglikelihood(combined_df=combined_df_asex, sim_column='simulation')
        loglik_df_asex.rename(columns={"loglikelihood": "loglike_asex"}, inplace=True)
        loglik_df_asex_bench = get_dens_loglikelihood(combined_df=combined_df_asex, sim_column='benchmark')
        loglik_df_asex_bench.rename(columns={"loglikelihood": "benchmark_loglike_asex"}, inplace=True)
        # todo: need review, this dataframe is created but not used. Should we use it in line 234(loglik_df = pd.merge(...))
        loglikelihood_comparison = pd.merge(loglik_df_asex, loglik_df_asex_bench, how="outer")

        loglik_df_gamet = get_dens_loglikelihood(combined_df=combined_df_gamet, sim_column='simulation')
        loglik_df_gamet.rename(columns={"loglikelihood": "loglike_gamet"}, inplace=True)
        loglik_df_gamet_bench = get_dens_loglikelihood(combined_df=combined_df_gamet, sim_column='benchmark')
        loglik_df_gamet_bench.rename(columns={"loglikelihood": "benchmark_loglike_gamet"}, inplace=True)
        loglik_df_gamet = pd.merge(loglik_df_gamet, loglik_df_gamet_bench, how="outer")

        loglik_df = pd.merge(loglik_df_asex, loglik_df_gamet, how="outer")
        loglik_df.to_csv(os.path.join(plot_output_filepath, 'loglikelihoods_par_dens.csv'),
                                        index=False)
        add_to_summary_table(combined_df=combined_df_asex, plot_output_filepath=plot_output_filepath,
                             validation_relationship_name='asexual_par_dens')
        add_to_summary_table(combined_df=combined_df_gamet, plot_output_filepath=plot_output_filepath,
                             validation_relationship_name='gamet_par_dens')


# Infectiousness to vectors
def generate_infectiousness_outputs(coord_csv, simulation_output_filepath, base_reference_filepath,
                                      plot_output_filepath, benchmark_simulation_filepath=None):

    """
    From simulation output and matched reference data, create plots and quantitative comparisons for all sites
    associated with the infectiousness-to-vectors validation relationship.
    Args:
        coord_csv (): A dataframe detailing the sites simulated for each validation relationship and the corresponding
                      reference dataset
        simulation_output_filepath (): The filepath where simulation output is located
        base_reference_filepath (): The filepath where reference datasets are located
        plot_output_filepath (): The filepath to the directory where plots should be created
        benchmark_simulation_filepath (): The filepath where benchmark simulation output is located. If None, no
                                          comparisons are made against benchmark simulations

    Returns:

    """

    combined_df = prepare_infect_df(coord_csv, simulation_output_filepath, base_reference_filepath,
                                    benchmark_simulation_filepath)

    plot_output = plot_infectiousness_ref_sim_comparison(combined_df)
    plot_list = plot_output[0]
    all_sites = plot_output[1]
    for ss in range(len(all_sites)):
        plot_list[ss].save(filename=os.path.join(plot_output_filepath,
                                                 'site_compare_infectiousness_' + all_sites[ss] + '.png'),
                           width=7.5, height=6, units='in')

    # compare simulation and benchmark simulation results
    if 'benchmark' in combined_df.columns:
        compare_benchmarks_output = compare_benchmark(combined_df)
        compare_benchmarks_output.save(filename=os.path.join(plot_output_filepath,
                                                             'scatter_benchmark_infectiousness.png'),
                                       height=4.5, width=8, units='in')
        # todo: add likelihood and other quantitative comparisons
        add_to_summary_table(combined_df=combined_df, plot_output_filepath=plot_output_filepath,
                             validation_relationship_name='infectiousness')


# Duration of infection
def generate_age_infection_duration_outputs(coord_csv, simulation_output_filepath, base_reference_filepath,
                                            plot_output_filepath, pos_thresh_dens=0.5, duration_bins=None,
                                            benchmark_simulation_filepath=None):
    """
    From simulation output and matched reference data, create plots and quantitative comparisons for all sites
    associated with the duration-of-infection validation relationship.
    Args:
        coord_csv (): A dataframe detailing the sites simulated for each validation relationship and the corresponding
                      reference dataset
        simulation_output_filepath (): The filepath where simulation output is located
        base_reference_filepath (): The filepath where reference datasets are located
        plot_output_filepath (): The filepath to the directory where plots should be created
        pos_thresh_dens (): A number giving the minimum true asexual parasite density a simulated individual must have
                            to be considered positive
        duration_bins (): A monotonically-increasing vector of numbers giving the plotted bin breaks for the duration
                          (in days) individuals remain infected
        benchmark_simulation_filepath (): The filepath where benchmark simulation output is located. If None, no
                                          comparisons are made against benchmark simulations


    Returns:

    """
    # TODO: add benchmark simulation support, add quantitative comparisons
    if not duration_bins:
        duration_bins = list(range(0, 400, 50))
        duration_bins.append(500)

    # determine which of the infectiousness sites have the relevant simulation output
    available_sites = get_available_sites_for_relationship(coord_csv, simulation_output_filepath,
                                                           relationship_name='infection_duration',
                                                           relationship_sim_filename='patient_reports.csv')

    for ss in range(len(available_sites)):
        cur_site = available_sites[ss]

        filepath_ref = os.path.join(base_reference_filepath, coord_csv[coord_csv['site'] == cur_site][
            'infection_duration_ref'].iloc[0])
        ref_df = pd.read_csv(filepath_ref)
        ref_df = ref_df[ref_df['site'].str.lower() == cur_site.lower()]
        ref_df['date'] = pd.to_datetime(ref_df['date'])

        sim_dir = os.path.join(simulation_output_filepath, cur_site)
        sim_data = get_sim_survey(sim_dir=sim_dir, ref_df=ref_df)

        # create and save comparison plots
        gg1 = plot_infection_duration_dist(ref_df=ref_df, sim_data=sim_data, pos_thresh_dens=pos_thresh_dens,
                                           duration_bins=duration_bins)
        gg2 = plot_infection_duration_dist_by_age(ref_df=ref_df, sim_data=sim_data, pos_thresh_dens=pos_thresh_dens,
                                                  duration_bins=duration_bins)
        gg3 = create_barplot_frac_comparison(ref_df=ref_df, sim_data=sim_data, pos_thresh_dens=pos_thresh_dens)

        gg1.save(filename=os.path.join(plot_output_filepath, 'site_compare_infect_duration_' + cur_site + '.png'),
                 height=4, width=8, units='in')
        gg2.save(filename=os.path.join(plot_output_filepath, 'site_compare_infect_duration_age_' + cur_site + '.png'),
                 height=5, width=8, units='in')
        gg3.save(filename=os.path.join(plot_output_filepath, 'site_compare_infect_duration_measures_' + cur_site +'.png'),
                 height=4, width=8, units='in')
