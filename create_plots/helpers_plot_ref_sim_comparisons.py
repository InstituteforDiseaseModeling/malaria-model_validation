# helpers_plot_ref_sim_comparisons.py
# These functions generate plots comparing between reference data and simulated outputs (and also between new and
# benchmark simulations).


from plotnine import ggplot, aes, geom_bar, scale_fill_brewer, facet_grid, geom_line, geom_point, geom_errorbar, \
    theme_bw, xlab, ylab, scale_color_manual, scale_fill_manual, coord_fixed, geom_abline, theme_classic, themes, \
    facet_wrap, scale_shape_manual, scale_size_manual, scale_x_log10, ggtitle, labs, position_dodge, scale_x_continuous
import numpy as np
import pandas as pd
from scipy.stats import beta
from pandas.api.types import CategoricalDtype
from datetime import datetime


# todo: not sure how to define color with rbg numbers. using the builtin colors for now
color_manual = {"reference": 'red',  # (169/255,23/255,23/255, 0.8),
                "simulation": 'blue',  # (0/255,124/255,180/255, 0.8),
                "benchmark": 'black'}  # (0 ,0 , 0, 0.8)}
color_manual_2 = color_manual.copy()
color_manual_2.pop("benchmark")
# todo: See matplotlib.markers. for list of all possible shapes.
# https://plotnine.readthedocs.io/en/stable/_modules/plotnine/scales/scale_shape.html
# https://matplotlib.org/stable/api/markers_api.html
shape_manual = {"reference": "P",  # plus #16,
                "simulation": "D",  # diamond #16,
                "benchmark":  6}  # careup #1}
size_manual = {"reference": 1.5,
               "simulation": 1.5,
               "benchmark": 3}
fill_manual = {"reference": 'red',
               "simulation": 'blue',
               "benchmark": '#ffffff00'}  # this should set a transparent fill color for benchmark


# compare new and benchmark simulation values
def compare_benchmark(combined_df):
    """
    Create scatter plots with new versus benchmark simulation output
    Args:
        combined_df (): A dataframe with both the simulation and matched benchmark values

    Returns: A gg scatterplot

    """
    # filter out nan values in dataframe
    combined_df = combined_df[~combined_df['benchmark'].isna()]
    metric = combined_df['metric'].iloc[0]
    if 'site_month' in combined_df.columns:
        combined_df['Site'] = combined_df['site_month']
    min_value = min(np.nanmin(combined_df['benchmark']), np.nanmin(combined_df['simulation']))
    max_value = max(np.nanmax(combined_df['benchmark']), np.nanmax(combined_df['simulation']))
    gg = (ggplot(combined_df, aes(x='benchmark', y='simulation', color='Site', fill='Site'))
          + geom_point(size=2)
          + ylab(f'new simulation {metric}')
          + xlab(f'benchmark sim {metric}')
          + ggtitle('Benchmark versus new sim values')
          + coord_fixed(ratio=1, xlim=(min_value, max_value), ylim=(min_value, max_value))
          + geom_abline(slope=1, intercept=0, color='grey', alpha=0.5)
          + theme_classic()
          + themes.theme(plot_title=themes.element_text(size=12)))
    return gg


def plot_ref_sim_comparison(combined_df, data_column_name):
    """
    Create a panel of line plots (one for each site-month) showing the data-by-age relationship seen in the
    reference and simulation datasets
    Args:
        combined_df (): A dataframe with the reference and simulation values
        data_column_name (): string for column name that contains the data

    Returns: A panel of ggplots

    """
    # convert dataframe to long format
    # todo: need code review
    # R code:
    # combined_df_long = pivot_longer(data=combined_df, cols=c('reference', 'simulation', 'benchmark'), names_to='source',
    #                                 values_to='incidence')
    # todo: the dataframe combined_df only contains one column called simulation, one reference and one benchmark,
    # pd.wide_to_long will not work is the stubnames are exactly as the column names
    # ValueError: stubname can't be identical to a column name.
    # rename the columns to add data_column_name as pre-fix
    """
    combined_rename_df = combined_df.rename(
    columns=lambda column: f"{data_column_name}_{column}" if column in ['reference', 'simulation', 'benchmark'] else column)
    combined_rename_df["id"] = combined_rename_df.index
    print(combined_rename_df.index.name)
    print(combined_rename_df.columns)

    combined_df_long = pd.wide_to_long(combined_rename_df, stubnames=f'{data_column_name}_',
                                       j='source', i=['id'])
    # this returns an empty dataframe
    """
    # todo: try with melt function instead
    # it seems to work
    id_vars = [column for column in combined_df.columns if column not in ['reference', 'simulation', 'benchmark']]
    combined_df_long = pd.melt(combined_df.reset_index(),
                               id_vars=id_vars,
                               #id_vars=['mean_age', 'Site', 'ref_pop_size', 'ref_year', 'metric'],
                               value_vars=['reference', 'simulation', 'benchmark'],
                               var_name='source', value_name=data_column_name
                               )
    combined_df_long = combined_df_long[combined_df_long[data_column_name].notnull()]
    combined_df_long = combined_df_long[combined_df_long['ref_year'].notnull()]
    facet_wrap_col = 'Site' if data_column_name == 'incidence' else 'site_month'
    gg = (ggplot(combined_df_long,
                 aes(x='mean_age', y=data_column_name, color='source', shape='source', group='ref_year'))
          # todo: need to find equivalent of interaction() in Python
          # + geom_line(aes(group=interaction('source', 'ref_year')))
          + geom_line(aes(group='source'))
          + geom_point(aes(size='source'), alpha=0.5)
          + scale_color_manual(values=color_manual)
          + scale_shape_manual(values=shape_manual)
          + scale_size_manual(values=size_manual)
          + scale_fill_manual(values=fill_manual)
          + xlab('age (midpoint of age bin)')
          + ylab(data_column_name)
          + facet_wrap(facet_wrap_col, ncol=4)
          + theme_bw()
          )

    return gg


# plot age-incidence comparisons with reference
def plot_inc_ref_sim_comparison(combined_df):
    """
        Create a panel of line plots (one for each site-month) showing the incidence-by-age relationship seen in the
        reference and simulation datasets
        Args:
            combined_df (): A dataframe with the reference and simulation values

        Returns: A panel of ggplots

        """
    return plot_ref_sim_comparison(combined_df, 'incidence')


# plot age-prevalence comparisons with reference
def plot_prev_ref_sim_comparison(combined_df):
    """
    Create a panel of line plots (one for each site-month) showing the prevalence-by-age relationship seen in the
    reference and simulation datasets
    Args:
        combined_df (): A dataframe with the reference and simulation values

    Returns: A panel of ggplots

    """
    return plot_ref_sim_comparison(combined_df, 'prevalence')


# plot parasite density comparisons with reference
def plot_par_dens_ref_sim_comparison(combined_df):
    """
    Create a plots (one for each site-month) showing the parasite density-by-age relationship seen in the reference
    and simulation datasets
    Args:
        combined_df (): A dataframe with the reference and simulation values

    Returns: A list with three elements:
           - 1) A panel of gg barplots showing the parasite density frequencies across age groups for all sites.
           - 2) A list of ggplots comparing the reference and simulation density frequencies. Each plot corresponds to a site.
           - 3) A vector of site names, with the same ordering as the list of plots (element 2)
    """
    months_of_year = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    # convert dataframe to long format
    # combined_df_long = pd.wide_to_long(df=combined_df, stubnames=['reference', 'simulation', 'benchmark'],
    #                                   j='source', i='density_frequency')
    # pd.wide_to_long doesn't allow stubname this is identical to a column name, try with pd.melt

    id_vars = [column for column in combined_df.columns if column not in ['reference', 'simulation', 'benchmark']]
    combined_df_long = pd.melt(combined_df.reset_index(),
                               id_vars=id_vars,
                               value_vars=['reference', 'simulation', 'benchmark'],
                               var_name='source', value_name='density_frequency'
                               )

    # = = = = = = = = = #
    # stacked barplots
    # = = = = = = = = = #
    # change type to factors for barplot groupings
    convert_dict = {'densitybin': 'category',
                    'mean_age': 'category'}
    combined_df_long = combined_df_long.astype(convert_dict)

    # colors
    # todo: BrBG only have 11 colors, not sure if it's enough, using scale_fill_brewer(palette = "BrBG") in ggplot for now
    # from bokeh.palettes import BrBG
    # len_density_bin = len(combined_df_long['densitybin'].unique())
    # num_colors = len_density_bin + 1 if len_density_bin % 2 == 0 else len_density_bin
    # if 3<= num_colors <= 11:
    #     colors = BrBG[num_colors]
    # names(colors) = sorted(combined_df['densitybin'].unique())

    # plot
    gg1 = (ggplot(combined_df_long, aes(fill='densitybin', y='density_frequency', x='mean_age'))
           + geom_bar(position="stack", stat="identity")
           # + scale_fill_manual(values=colors, limits=names(colors))
           + scale_fill_brewer(type='div', palette="BrBG")
           + facet_grid('site_month~source'))

    # = = = = = = = = = = = = = = = = = = #
    # grid of line plots - one plot panel per site
    # = = = = = = = = = = = = = = = = = = #
    line_plot_list = list()
    all_sites = combined_df_long['Site'].unique()
    for ss in range(len(all_sites)):
        cur_site = all_sites[ss]
        combined_df = combined_df_long[combined_df_long['Site'] == cur_site]

        # calculate reference error bounds using Jerrerys interval
        ci_width = 0.95
        alpha = 1 - ci_width
        combined_df['min_ref'] = np.nan
        combined_df['max_ref'] = np.nan
        eligible_rows =((combined_df['ref_bin_count'] > 0) &
                                     (combined_df['ref_bin_count'] < combined_df['ref_total']) &
                                     (combined_df['source'] == 'reference'))
        combined_df['min_ref'][eligible_rows] = beta.ppf(q=alpha / 2,
                                                         a=combined_df['ref_bin_count'][eligible_rows] + 0.5,
                                                         b=combined_df['ref_total'][eligible_rows] - combined_df['ref_bin_count'][eligible_rows] + 0.5)
        combined_df['max_ref'][eligible_rows] = beta.ppf(q=1 - alpha / 2,
                                                         a=combined_df['ref_bin_count'][eligible_rows] + 0.5,
                                                         b=combined_df['ref_total'][eligible_rows] - combined_df['ref_bin_count'][eligible_rows] + 0.5)

        # change facet values to intuitive labels
        combined_df['month'] = combined_df['month'].apply(lambda x: months_of_year[x-1])
        month_cat = CategoricalDtype(categories=months_of_year, ordered=True)
        combined_df['month'] = combined_df['month'].astype(month_cat)
        all_age_bins = sorted(combined_df['agebin'].unique())
        age_bin_labels = ['<=' + str(all_age_bins[0]) + " years"]
        for aa in range(len(all_age_bins) - 1):
            age_bin_labels.append(str(all_age_bins[aa]) + '-' + str(all_age_bins[aa + 1]) + ' years')

        combined_df['agebin_index'] = combined_df['agebin'].map(dict(zip(all_age_bins, range(len(all_age_bins)))))
        combined_df['agebin'] = combined_df['agebin_index'].apply(lambda x: age_bin_labels[x])
        age_bin_labels_cat = CategoricalDtype(categories=age_bin_labels, ordered=True)
        combined_df['agebin'] = combined_df['agebin'].astype(age_bin_labels_cat)

        # plot lineplot of simulation and reference densities
        gg2 = (ggplot(combined_df, aes(x='densitybin', y='density_frequency', color='source'))
               + geom_line(aes(group='source'), size=1)
               + geom_point(aes(size='source', shape='source'))
               + geom_errorbar(aes(ymin='min_ref', ymax='max_ref'), width=0.2)
               # + scale_x_continuous(trans='log10')
               + theme_bw()
               + ylab('fraction of population')
               + xlab('parasite density bin')
               + themes.theme(axis_text_x=themes.element_text(angle=45))
               + ggtitle(cur_site)
               + scale_color_manual(values=color_manual)
               + scale_shape_manual(values=shape_manual)
               + scale_size_manual(values=size_manual)
               + facet_grid('month~agebin')) # 'agebin~month' will throw error when try to say this plot: IndexError: index 0 is out of bounds for axis 0 with size 0

        line_plot_list.append(gg2)

    return gg1, line_plot_list, all_sites


#create infectiousness plots
def plot_infectiousness_ref_sim_comparison(combined_df):
    """
    Create a plots (one for each site-month) showing the infectiousness-to-vectors by age and parasite density
    relationship seen in the reference and simulation datasets
    Args:
        combined_df (): A dataframe with the reference and simulation values

    Returns: A list with two elements:
            - 1) A list of ggplots comparing the reference and simulation density frequencies. Each plot corresponds to a site.
            - 2) A vector of site names, with the same ordering as the list of plots (element 2)

    """
    months_of_year = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    # convert dataframe to long format
    id_vars = [column for column in combined_df.columns if column not in ['reference', 'simulation', 'benchmark']]
    combined_df_long = pd.melt(combined_df.reset_index(),
                               id_vars=id_vars,
                               value_vars=['reference', 'simulation', 'benchmark'],
                               var_name='source', value_name='infectiousness_bin_freq'
                               )

    line_plot_list = list()
    all_sites = combined_df_long['Site'].unique()
    for ss in range(len(all_sites)):
        cur_site = all_sites[ss]
        combined_df = combined_df_long[combined_df_long['Site'] == cur_site]

        # change facet values to intuitive labels
        combined_df['month'] = combined_df['month'].apply(lambda x: months_of_year[x-1])
        month_cat = CategoricalDtype(categories=months_of_year, ordered=True)
        combined_df['month'] = combined_df['month'].astype(month_cat)
        all_age_bins = sorted(combined_df['agebin'].unique())
        age_bin_labels = ['<=' + str(all_age_bins[0]) + " years"]
        for aa in range(len(all_age_bins) - 1):
            age_bin_labels.append(str(all_age_bins[aa]) + '-' + str(all_age_bins[aa + 1]) + ' years')

        combined_df['agebin_index'] = combined_df['agebin'].map(dict(zip(all_age_bins, range(len(all_age_bins)))))
        combined_df['agebin'] = combined_df['agebin_index'].apply(lambda x: age_bin_labels[x])
        age_bin_labels_cat = CategoricalDtype(categories=age_bin_labels, ordered=True)
        combined_df['agebin'] = combined_df['agebin'].astype(age_bin_labels_cat)

        # plot lineplot of simulation and reference densities
        gg2 = (ggplot(combined_df, aes(x='densitybin', y='fraction_infected_bin', color='source',
                                       size='infectiousness_bin_freq', fill='source'))
               + geom_point(alpha=0.5) #pch=21
               + scale_x_log10()
               + ylab('percent of mosquitoes infected upon feeding')
               + xlab('gametocyte density')
                # todo: confirm what does labs(size=) do in R
               # + labs(size='fraction of individuals')
               + ggtitle(cur_site)
               + scale_color_manual(values=color_manual)
               + scale_fill_manual(values=fill_manual)
               + facet_grid('month~agebin'))

        line_plot_list.append(gg2)

    return line_plot_list, all_sites


# infection duration plotting functions
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# helper functions for pre-plotting data processing
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
def get_frac_state_swaps(data, pos_thresh_dens):
    """
    Calculate probability of going from negative to positive or from positive to negative between sample dates (among
    surveyed individuals)
    Args:
        data (): A dataframe of test results for all individuals and survey days
        pos_thresh_dens ():

    Returns: A vector where the first element is the fraction of the time a positive test was followed by a negative
    result and the second element is the fraction of time a negative test was followed by a positive test

    """
    # brute force approach iterating through people and dates
    indIDs = data['SID'].unique()
    sum_denom_pos = 0
    sum_turn_neg = 0
    sum_denom_neg = 0
    sum_turn_pos = 0

    for ii in range(len(indIDs)):
        data_cur = data[data['SID'] == indIDs[ii]]
        data_cur = data_cur.sort_values(by=['date'])
        # todo: need code review
        # get tests results as a series of bool values
        ind_pos = data_cur['DENSITY'] > pos_thresh_dens
        pos_count = sum(ind_pos)
        neg_count = len(ind_pos) - pos_count

        # denominators for each (number of each type that had an observation after then (i.e., they could have been observed to change))
        last_obs_pos = (data_cur['DENSITY'].iloc[len(data_cur) - 1] > pos_thresh_dens)
        sum_denom_pos = sum_denom_pos + pos_count - 1 if last_obs_pos else sum_denom_pos + pos_count
        sum_denom_neg = sum_denom_neg + neg_count if last_obs_pos else sum_denom_neg + neg_count -1

        # find how many tests change from neg to pos or pos to neg across timesteps
        for i in range(len(ind_pos) - 1):
            if ind_pos.iloc[i] and not ind_pos.iloc[i + 1]: # change from positive to negative
                sum_turn_neg += 1
            elif not ind_pos.iloc[i] and ind_pos.iloc[i + 1]: # change from negative to positive
                sum_turn_pos += 1

    frac_pos_turn_neg_next_time = sum_turn_neg / sum_denom_pos
    frac_neg_turn_pos_next_time = sum_turn_pos / sum_denom_neg

    return [frac_pos_turn_neg_next_time, frac_neg_turn_pos_next_time]


def get_time_pos(data, pos_thresh_dens):
    """
    Calculate infection duration (time between consecutive positive samples) and whether or not the infection
    observation was censored
    Args:
        data (): A dataframe of test results for all individuals and survey days
        pos_thresh_dens (): A number giving the minimum true asexual parasite density a simulated individual must have
                            to be considered positive

    Returns: A dataframe where each row corresponds to a stretch of time when an individual has uninterrupted positive tests

    """
    # brute force approach iterating through people and dates
    indIDs = data['SID'].unique()
    days_positive = list()  # number of days between the first and last sample date in a run of positive samples
    sample_censored = list()  # Boolean indicating whether the positivity duration might have been longer unobserved
    sample_ages = list()
    sample_seeds = list()

    if not ('seed' in data.columns):
        data['seed'] = 0
    for ss in data['seed'].unique():
        data_ss = data[data['seed'] == ss]
        for ii in range(len(indIDs)):
            data_cur = data_ss[data_ss['SID'] == indIDs[ii]]
            data_cur = data_cur.sort_values(by='date')
            cur_age = data_cur['age'].mean()
            # which samples had positive tests
            pos_bool = data_cur['DENSITY'] > pos_thresh_dens # todo, unresolved reference, maybe add it as an argument

            # todo: need code review
            def rle(series):
                """
                implement the rle(Run Length Encoding) function in Python for pandas series
                Args:
                    series (): A Pandas series

                Returns: a dictionary contains 2 keys lengths and values:
                        Value of key lengths is an integer list containing the length of each run.
                        Value of key values is a list of the same length as lengths with the corresponding values.

                """
                i = 0
                res_dict = {'lengths': [], 'values': []}
                while i <= len(series) - 1:
                    count = 1
                    value = series.iloc[i]
                    j = i
                    while j < len(series) - 1:
                        if series.iloc[j] == series.iloc[j + 1]:
                            count = count + 1
                            j = j + 1
                        else:
                            break
                    res_dict['lengths'].append(count)
                    res_dict['values'].append(value)
                    i = j + 1
                return res_dict

            # get the length of spans of positive and of negative tests in a row
            num_in_a_row = rle(pos_bool)
            start_index_num_in_a_row = [0] + list(np.cumsum(num_in_a_row['lengths']))[:-1]

            # for each span of positives, determine the length before turning negative as well as whether it was censored or not
            # if it includes the first or last sample, it's censored (it's at least that many days but may have been longer)
            # todo: figure out what the next line is for
            # num_in_a_row['lengths'][num_in_a_row['values']

            # iterate through periods of positivity
            pos_index = np.where(np.array(num_in_a_row['values']))[0]
            for tt in pos_index:
                # indicate whether this period of positivity was censored
                # todo: need code review: the 4 lists are empty lists, should use append() instead of assigning by index
                sample_censored.append(tt == 0 or tt == len(num_in_a_row['lengths']))
                sample_ages.append(cur_age)
                sample_seeds.append(ss)

                # get the first and last positive sample day observed for this period of positivity
                first_index = start_index_num_in_a_row[tt]
                last_index = first_index + num_in_a_row['lengths'][tt] - 1  # todo: why -1
                first_positive_sample_day = data_cur['date'].iloc[first_index]
                last_positive_sample_day = data_cur['date'].iloc[last_index]
                if isinstance(last_positive_sample_day, str):
                    last_positive_sample_day = datetime.strptime(last_positive_sample_day, '%Y-%m-%d')
                if isinstance(first_positive_sample_day, str):
                    first_positive_sample_day = datetime.strptime(first_positive_sample_day, '%Y-%m-%d')
                days_positive.append((last_positive_sample_day - first_positive_sample_day).days)

    return pd.DataFrame({'days_positive': days_positive,
                         'censored': sample_censored,
                         'age': sample_ages,
                         'seed': sample_seeds})


def bin_durations_all_seeds(days_positive, duration_bins):
    """
    Get binned durations of infections for a specific group of infections (e.g., age group and censored/non-censored),
    with quantile ranges across sim seeds
    Args:
        days_positive (): A data frame where each row corresponds to a stretch of time when an individual has
                          uninterrupted positive tests
        duration_bins (): A monotonically-increasing vector of numbers giving the plotted bin breaks for the duration
                          (in days) individuals remain infected

    Returns: A dataframe with the fraction of infections that fall in each of the duration bins (average value across
             seeds), along with the extreme quantiles among all seeds

    """
    # create data frame of binned days positive for each seed in days_positive
    bin_min = duration_bins[:-1]
    bin_max = duration_bins[1:]
    bin_mid = [sum(x) / 2 for x in zip(bin_min, bin_max)]
    bin_w_seeds = pd.DataFrame({'bin_mid': bin_mid,
                                'bin_min': bin_min,
                                'bin_max': bin_max})
    for seed in days_positive['seed'].unique():
        days_positive_cur = days_positive[days_positive['seed'] == seed]
        # todo: need code review np.histogram() vs hist()
        binned_counts = np.histogram(np.array(days_positive_cur['days_positive']), bins=duration_bins, density=True)
        binned_dens = list(binned_counts[0])
        bin_w_seeds[f'seed_{seed}'] = binned_dens

    # todo, need code review
    num_seed_columns = ['seed_' in i for i in bin_w_seeds.columns]
    seed_columns = [i for i in bin_w_seeds.columns if 'seed_' in i]
    if any(num_seed_columns):
        bin_w_seeds['density'] = bin_w_seeds[seed_columns].mean(axis=1)
        bin_w_seeds['quant_low'] = bin_w_seeds[seed_columns].quantile(0.0, axis=1)  # 0.1
        bin_w_seeds['quant_high'] = bin_w_seeds[seed_columns].quantile(1.0, axis=1)    # 0.9
    elif sum(num_seed_columns) == 1:
        bin_w_seeds['density'] = bin_w_seeds[seed_columns]
        bin_w_seeds['quant_low'] = np.nan
        bin_w_seeds['quant_high'] = np.nan
    else:  # no data present for this binning
        bin_w_seeds['density'] = np.nan
        bin_w_seeds['quant_low'] = np.nan
        bin_w_seeds['quant_high'] = np.nan

    bin_df = bin_w_seeds[['bin_mid', 'density', 'quant_low', 'quant_high', 'bin_min', 'bin_max']]
    return bin_df


def get_age_group(cur_age, age_bin_lower, age_bin_labels):
    """
    Get the appropriate age bin label for a particular age in years
    Args:
        cur_age (): A numeric value representing an individual's age in years
        age_bin_lower (): A vector of the lower age ranges of each age bin
        age_bin_labels (): A vector of age bin labels, with index-wise correspondance with age_bin_lower

    Returns: The age-bin label corresponding to an age in years

    """
    res_labels = [age_bin_labels[i] for i in range(len(age_bin_lower)) if cur_age > age_bin_lower[i]]
    return res_labels[-1] if len(res_labels) else age_bin_labels[0]


def get_duration_bins(ind_data, duration_bins, pos_thresh_dens, facet_censored=True, facet_age=True, age_bin_lower=None):
    """
    Get fraction of infections that fall in each duration bin for an entire dataset, with quantile ranges across sim
    seeds, optionally faceted by censorship and age
    Args:
        ind_data (): A dataframe with all individual infection duration data
        duration_bins (): A monotonically-increasing vector of numbers giving the plotted bin breaks for the duration
                         (in days) individuals remain infected
        pos_thresh_dens (): A number giving the minimum true asexual parasite density a simulated individual must have
                            to be considered positive
        facet_censored (): A boolean value indicating whether results should be partitioned by whether or not the
                           infection duration observation was censored
        facet_age (): A boolean value indicating whether results should be partitioned by age group
        age_bin_lower (): A vector of the lower age ranges of each age bin to use if faceting by age group

    Returns: A dataframe with the fraction of infections that fall in each of the duration bins (average value across
            seeds), along with the extreme quantiles among all seeds.

    """
    if age_bin_lower is None:
        age_bin_lower = [0, 5, 10, 20, 100]

    if facet_age:
        age_bin_labels = [f'{age_bin_lower[i]}-{age_bin_lower[i+1]}' for i in range(len(age_bin_lower) - 1)]
        age_bin_labels.append(f'>{age_bin_lower[-1]}')

    # get data frame of the days each observed infection lasted, also recording age, whether infection was censored, and seed
    days_positive = get_time_pos(data=ind_data, pos_thresh_dens=pos_thresh_dens)
    if facet_age:
        days_positive['age_group'] = days_positive['age'].apply(get_age_group,
                                                                age_bin_lower=age_bin_lower,
                                                                age_bin_labels=age_bin_labels)

    # determine how dataset should be divided up among factors/facets
    bin_df = pd.DataFrame()
    if facet_censored and facet_age:
        for v1 in days_positive['censored'].unique():
            for v2 in days_positive['age_group'].unique():
                bins_cur = bin_durations_all_seeds(days_positive[(days_positive['censored'] == v1) & (days_positive['age_group'] == v2)],
                                                   duration_bins=duration_bins)
                bins_cur['censored'] = v1
                bins_cur['age_group'] = v2
                bin_df = pd.concat([bin_df, bins_cur])
    elif facet_censored:
        for v1 in days_positive['censored'].unique():
            bins_cur = bin_durations_all_seeds(days_positive[days_positive['censored'] == v1], duration_bins=duration_bins)
            bins_cur['censored'] = v1
            bins_cur['age_group'] = 'combined'
            bin_df = pd.concat([bin_df, bins_cur])
    elif facet_age:
        for v2 in days_positive['age_group'].unique():
            bins_cur = bin_durations_all_seeds(days_positive[days_positive['age_group'] == v2], duration_bins=duration_bins)
            bins_cur['censored'] = 'combined'
            bins_cur['age_group'] = v2
            bin_df = pd.concat([bin_df, bins_cur])
    else:
        bins_cur = bin_durations_all_seeds(days_positive, duration_bins=duration_bins)
        bins_cur['censored'] = 'combined'
        bins_cur['age_group'] = 'combined'
        bin_df = bins_cur

    if facet_age:
        age_bin_labels_cat = CategoricalDtype(categories=age_bin_labels, ordered=True)
        bin_df['age_group'] = bin_df['age_group'].astype(age_bin_labels_cat)

    return bin_df


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# plotting functions
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
def create_barplot_frac_comparison(ref_df, sim_data, pos_thresh_dens):
    """
    Create a gg barplot comparing the reference dataset and matching subsampled simulations for fraction of samples
    positive and fractions of samples switching from neg--> pos or pos--> neg
    Args:
        ref_df (): A dataframe with the reference values.
        sim_data (): A dataframe with the subsampled simulation values (may include results from multiple seeds)
        pos_thresh_dens (): A number giving the minimum true asexual parasite density a simulated individual must have
                            to be considered positive

    Returns: A gg barplot

    """

    frac_swap_ref = get_frac_state_swaps(data=ref_df, pos_thresh_dens=pos_thresh_dens)
    frac_swap_sim = get_frac_state_swaps(data=sim_data, pos_thresh_dens=pos_thresh_dens)

    ref_df['pos'] = ref_df['DENSITY'] > pos_thresh_dens
    frac_samples_pos_ref = sum(ref_df['pos']) / sum((ref_df['DENSITY'] > -1))

    sim_data['pos'] = sim_data['DENSITY'] > pos_thresh_dens
    frac_samples_pos_sim = sum(sim_data['pos']) / sum(sim_data['DENSITY'] > -1)

    df = pd.DataFrame({'source': ['reference'] * 3 + ['simulation'] * 3,
                       'measure': ['negative to positive next time',
                                   'positive to negative next time',
                                   'total samples positive'] * 2,
                       'value': [frac_swap_ref[1],
                                 frac_swap_ref[0],
                                 frac_samples_pos_ref,
                                 frac_swap_sim[1],
                                 frac_swap_sim[0],
                                 frac_samples_pos_sim]})

    gg = (ggplot(df, aes(x='measure', y='value', fill='source', color='source'))
          + geom_bar(stat="identity", position="dodge", alpha=.3)
          + scale_color_manual(values=color_manual_2)
          + labs(title='general dataset properties', y='fraction of samples', x=''))

    return gg


def plot_infection_duration_dist(ref_df, sim_data, pos_thresh_dens, duration_bins=None):
    """
    Create a panel of gg barplots comparing the distributions of infection lengths in simulation versus reference datasets
    Args:
        ref_df (): A dataframe with the reference values.
        sim_data (): A dataframe with the subsampled simulation values (may include results from multiple seeds)
        pos_thresh_dens (): A number giving the minimum true asexual parasite density a simulated individual must have
                            to be considered positive
        duration_bins (): A monotonically-increasing vector of numbers giving the plotted bin breaks for the duration
                          (in days) individuals remain infected

    Returns: A panel of gg barplots, faceted by censorship status

    """
    if not duration_bins:
        duration_bins = list(range(0, 400, 50))
        duration_bins.append(500)

    # get densities for each duration bin
    ref_bin_df = get_duration_bins(ind_data=ref_df, duration_bins=duration_bins, facet_censored=True, facet_age=False,
                                   pos_thresh_dens=pos_thresh_dens)
    ref_bin_df['dataset'] = 'reference'
    ref_bin_df['censor_type'] = 'start & finish observed'
    # todo: need code review, not sure if the Python code is correct or not
    # ref_bin_df$censor_type[ref_bin_df$censored] = 'censored'
    ref_bin_df['censor_type'][ref_bin_df['censored']] = 'censored'

    # combine multiple simulation seeds to show mean, 10%, 90% quantile values across seed distributions
    sim_bin_df = get_duration_bins(ind_data=sim_data, duration_bins=duration_bins, facet_censored=True, facet_age=False,
                                   pos_thresh_dens=pos_thresh_dens)
    sim_bin_df['dataset'] = 'simulation'
    sim_bin_df['censor_type'] = 'start & finish observed'
    sim_bin_df['censor_type'][sim_bin_df['censored']] = 'censored'

    # combine reference and simulation datasets
    bin_df = pd.concat([ref_bin_df, sim_bin_df])

    # max_y=0.8
    gg = (ggplot(bin_df, aes(x='bin_mid', y='density', color='dataset', fill='dataset'))
          + geom_bar(stat="identity", position="identity", alpha=.3)
          + geom_errorbar(aes(ymin='quant_low', ymax='quant_high'), width=.2, position=position_dodge(.9))
          + scale_color_manual(values=color_manual_2)
          + scale_fill_manual(values=color_manual_2)
          + labs(title='infection duration', x='infection duration (days)')
            # ylim(NA,max_y) +
          + facet_wrap(scales='fixed', facets='censor_type'))
    return gg


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# create plots comparing the distributions of infection lengths, faceted by age group
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
def plot_infection_duration_dist_by_age(ref_df, sim_data, pos_thresh_dens, age_bin_lower=None,
                                        duration_bins=None):
    """
    Create a panel of gg barplots comparing the distributions of infection lengths in simulation versus reference
    datasets
    Args:
        ref_df (): A dataframe with the reference values.
        sim_data (): A dataframe with the subsampled simulation values (may include results from multiple seeds)
        pos_thresh_dens (): A number giving the minimum true asexual parasite density a simulated individual must have
                            to be considered positive
        age_bin_lower (): A vector of the lower age ranges of each age bin to use if faceting by age group
        duration_bins (): A monotonically-increasing vector of numbers giving the plotted bin breaks for the duration
                         (in days) individuals remain infected

    Returns: A panel of gg barplots, faceted by censorship status and age bin

    """
    if age_bin_lower is None:
        age_bin_lower = [0, 5, 10, 20, 100]
    if not duration_bins:
        duration_bins = list(range(0, 400, 50))
        duration_bins.append(500)

    # get densities for each duration bin
    ref_bin_df = get_duration_bins(ind_data=ref_df, duration_bins=duration_bins, facet_censored=True, facet_age=True,
                                   pos_thresh_dens=pos_thresh_dens, age_bin_lower=age_bin_lower)
    ref_bin_df['dataset'] = 'reference'
    ref_bin_df['censor_type'] = 'start & finish observed'
    ref_bin_df['censor_type'][ref_bin_df['censored']] = 'censored'

    # combine multiple simulation seeds to show mean, 10%, 90% quantile values across seed distributions
    sim_bin_df = get_duration_bins(ind_data=sim_data, duration_bins=duration_bins, facet_censored=True, facet_age=True,
                                   pos_thresh_dens=pos_thresh_dens, age_bin_lower=age_bin_lower)
    sim_bin_df['dataset'] = 'simulation'
    sim_bin_df['censor_type'] = 'start & finish observed'
    sim_bin_df['censor_type'][sim_bin_df['censored']] = 'censored'

    # combine reference and simulation datasets
    bin_df = pd.concat([ref_bin_df, sim_bin_df])

    # max_y=0.9
    gg = (ggplot(bin_df, aes(x='bin_mid', y='density', color='dataset', fill='dataset'))
          + geom_bar(stat="identity", position="identity", alpha=.3)
          + geom_errorbar(aes(ymin='quant_low', ymax='quant_high'), width=.2, position=position_dodge(.9))
          + scale_color_manual(values=color_manual_2)
          + scale_fill_manual(values=color_manual_2)
          + labs(title='infection duration', x='infection duration (days)')
          # ylim(NA,max_y) +
          + facet_grid(scales='fixed', facets=['censor_type', 'age_group']))

    return gg
