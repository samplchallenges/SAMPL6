#!/usr/bin/env python

# =============================================================================
# GLOBAL IMPORTS
# =============================================================================
import os
import glob
import io
import collections

import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import scipy.stats

# =============================================================================
# CONSTANTS
# =============================================================================

# Paths to input data.
LOGP_SUBMISSIONS_DIR_PATH = './logP_predictions'
EXPERIMENTAL_DATA_FILE_PATH = '../experimental_data/logP_experimental_values.csv'

# =============================================================================
# STATS FUNCTIONS
# =============================================================================

def r2(data):
    x, y = data.T
    slope, intercept, r_value, p_value, stderr = scipy.stats.linregress(x, y)
    return r_value**2


def slope(data):
    x, y = data.T
    slope, intercept, r_value, p_value, stderr = scipy.stats.linregress(x, y)
    return slope


def me(data):
    x, y = data.T
    error = np.array(x) - np.array(y)
    return error.mean()


def mae(data):
    x, y = data.T
    error = np.abs(np.array(x) - np.array(y))
    return error.mean()


def rmse(data):
    x, y = data.T
    error = np.array(x) - np.array(y)
    rmse = np.sqrt((error**2).mean())
    return rmse


def compute_bootstrap_statistics(samples, stats_funcs, percentile=0.95, n_bootstrap_samples=1000):
    """Compute bootstrap confidence interval for the given statistics functions."""
    # Handle case where only a single function is passed.
    try:
        len(stats_funcs)
    except TypeError:
        stats_funcs = [stats_funcs]

    # Compute mean statistics.
    statistics = [stats_func(samples) for stats_func in stats_funcs]

    # Generate bootstrap statistics.
    bootstrap_samples_statistics = np.zeros((len(statistics), n_bootstrap_samples))
    for bootstrap_sample_idx in range(n_bootstrap_samples):
        samples_indices = np.random.randint(low=0, high=len(samples), size=len(samples))
        for stats_func_idx, stats_func in enumerate(stats_funcs):
            bootstrap_samples_statistics[stats_func_idx][bootstrap_sample_idx] = stats_func(samples[samples_indices])

    # Compute confidence intervals.
    percentile_index = int(np.floor(n_bootstrap_samples * (1 - percentile) / 2)) - 1
    bootstrap_statistics = []
    for stats_func_idx, samples_statistics in enumerate(bootstrap_samples_statistics):
        samples_statistics.sort()
        stat_lower_percentile = samples_statistics[percentile_index]
        stat_higher_percentile = samples_statistics[-percentile_index+1]
        confidence_interval = (stat_lower_percentile, stat_higher_percentile)
        bootstrap_statistics.append([statistics[stats_func_idx], confidence_interval, samples_statistics])

    return bootstrap_statistics

# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

def plot_correlation(x, y, data, title=None, color=None, kind='joint', ax=None):
    # Extract only logP values.
    data = data[[x, y]]

    # Find extreme values to make axes equal.
    min_limit = np.ceil(min(data.min()) - 1)
    max_limit = np.floor(max(data.max()) + 1)
    axes_limits = np.array([min_limit, max_limit])

    if kind == 'joint':
        grid = sns.jointplot(x=x, y=y, data=data,
                             kind='reg', joint_kws={'ci': None}, stat_func=None,
                             xlim=axes_limits, ylim=axes_limits, color=color)
        ax = grid.ax_joint
        grid.fig.subplots_adjust(top=0.95)
        grid.fig.suptitle(title)
    elif kind == 'reg':
        ax = sns.regplot(x=x, y=y, data=data, color=color, ax=ax)
        ax.set_title(title)

    # Add diagonal line.
    ax.plot(axes_limits, axes_limits, ls='--', c='black', alpha=0.8, lw=0.7)

    # Add shaded area for 0.5-1 logP error.
    palette = sns.color_palette('BuGn_r')
    ax.fill_between(axes_limits, axes_limits - 0.5, axes_limits + 0.5, alpha=0.2, color=palette[2])
    ax.fill_between(axes_limits, axes_limits - 1, axes_limits + 1, alpha=0.2, color=palette[3])


def plot_correlation_with_SEM(x_lab, y_lab, x_err_lab, y_err_lab, data, title=None, color=None, ax=None):
    # Extract only logP values.
    x_error = data.loc[:, x_err_lab]
    y_error = data.loc[:, y_err_lab]
    x_values = data.loc[:, x_lab]
    y_values = data.loc[:, y_lab]
    data = data[[x_lab, y_lab]]

    # Find extreme values to make axes equal.
    min_limit = np.ceil(min(data.min()) - 1)
    max_limit = np.floor(max(data.max()) + 1)
    axes_limits = np.array([min_limit, max_limit])

    # Color
    current_palette = sns.color_palette()
    sns_blue = current_palette[0]

    # Plot
    plt.figure(figsize=(6, 6))
    grid = sns.regplot(x=x_values, y=y_values, data=data, color=color, ci=None)
    plt.errorbar(x=x_values, y=y_values, xerr=x_error, yerr=y_error, fmt="o", ecolor=sns_blue, capthick='2',
                 label='SEM', alpha=0.75)
    plt.axis("equal")

    if len(title) > 70:
        plt.title(title[:70]+"...")
    else:
        plt.title(title)

    # Add diagonal line.
    grid.plot(axes_limits, axes_limits, ls='--', c='black', alpha=0.8, lw=0.7)

    # Add shaded area for 0.5-1 logP error.
    palette = sns.color_palette('BuGn_r')
    grid.fill_between(axes_limits, axes_limits - 0.5, axes_limits + 0.5, alpha=0.2, color=palette[2])
    grid.fill_between(axes_limits, axes_limits - 1, axes_limits + 1, alpha=0.2, color=palette[3])

    plt.xlim(axes_limits)
    plt.ylim(axes_limits)


def barplot_with_CI_errorbars(df, x_label, y_label, y_lower_label, y_upper_label, figsize=False):
    """Creates bar plot of a given dataframe with asymmetric error bars for y axis.

    Args:
        df: Pandas Dataframe that should have columns with columnnames specified in other arguments.
        x_label: str, column name of x axis categories
        y_label: str, column name of y axis values
        y_lower_label: str, column name of lower error values of y axis
        y_upper_label: str, column name of upper error values of y axis
        figsize: tuple, size in inches. Default value is False.

    """
    # Column names for new columns for delta y_err which is calculated as | y_err - y |
    delta_lower_yerr_label = "$\Delta$" + y_lower_label
    delta_upper_yerr_label = "$\Delta$" + y_upper_label
    data = df  # Pandas DataFrame
    data.loc[:,delta_lower_yerr_label] = data.loc[:,y_label] - data.loc[:,y_lower_label]
    data.loc[:,delta_upper_yerr_label] = data.loc[:,y_upper_label] - data.loc[:,y_label]

    # Color
    current_palette = sns.color_palette()
    sns_color = current_palette[2]

    # Plot style
    plt.close()
    plt.style.use(["seaborn-talk", "seaborn-whitegrid"])
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 16
    #plt.tight_layout()

    # If figsize is specified
    if figsize != False:
        plt.figure(figsize=figsize)

    # Plot
    x = range(len(data[y_label]))
    y = data[y_label]
    plt.bar(x, y)
    plt.xticks(x, data[x_label], rotation=90)
    plt.errorbar(x, y, yerr=(data[delta_lower_yerr_label], data[delta_upper_yerr_label]),
                 fmt="none", ecolor=sns_color, capsize=3, capthick=True)
    plt.xlabel(x_label)
    plt.ylabel(y_label)


def barplot_with_CI_errorbars_colored_by_label(df, x_label, y_label, y_lower_label, y_upper_label, color_label, figsize=False):
    """Creates bar plot of a given dataframe with asymmetric error bars for y axis.

        Args:
            df: Pandas Dataframe that should have columns with columnnames specified in other arguments.
            x_label: str, column name of x axis categories
            y_label: str, column name of y axis values
            y_lower_label: str, column name of lower error values of y axis
            y_upper_label: str, column name of upper error values of y axis
            color_label: str, column name of label that will determine the color of bars
            figsize: tuple, size in inches. Default value is False.

        """
    # Column names for new columns for delta y_err which is calculated as | y_err - y |
    delta_lower_yerr_label = "$\Delta$" + y_lower_label
    delta_upper_yerr_label = "$\Delta$" + y_upper_label
    data = df  # Pandas DataFrame
    data.loc[:, delta_lower_yerr_label] = data.loc[:, y_label] - data.loc[:, y_lower_label]
    data.loc[:, delta_upper_yerr_label] = data.loc[:, y_upper_label] - data.loc[:, y_label]

    # Color
    current_palette = sns.color_palette()
    # Error bar color
    sns_color = current_palette[2]
    # Bar colors
    category_list = ["Physical", "Empirical", "Mixed", "Other"]
    bar_color_dict = {}
    for i, cat in enumerate(category_list):
        bar_color_dict[cat] = current_palette[i]
    print("bar_color_dict:\n", bar_color_dict)


    # Plot style
    plt.close()
    plt.style.use(["seaborn-talk", "seaborn-whitegrid"])
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 16
    # plt.tight_layout()

    # If figsize is specified
    if figsize != False:
        plt.figure(figsize=figsize)

    # Plot
    x = range(len(data[y_label]))
    y = data[y_label]
    #barlist = plt.bar(x, y)
    fig, ax = plt.subplots(figsize=figsize)
    barlist = ax.bar(x, y)

    plt.xticks(x, data[x_label], rotation=90)
    plt.errorbar(x, y, yerr=(data[delta_lower_yerr_label], data[delta_upper_yerr_label]),
                 fmt="none", ecolor='gray', capsize=3, elinewidth=2, capthick=True)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    # Reset color of bars ased on color label
    print("data.columns:\n",data.columns)
    for i, c_label in enumerate(data.loc[:, color_label]):
        barlist[i].set_color(bar_color_dict[c_label])

    # create legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color=bar_color_dict["Physical"], lw=5),
                    Line2D([0], [0], color=bar_color_dict["Empirical"], lw=5),
                    Line2D([0], [0], color=bar_color_dict["Mixed"], lw=5),
                    Line2D([0], [0], color=bar_color_dict["Other"], lw=5)]
    ax.legend(custom_lines, category_list)


def barplot(df, x_label, y_label, title):
    """Creates bar plot of a given dataframe.

    Args:
        df: Pandas Dataframe that should have columns with columnnames specified in other arguments.
        x_label: str, column name of x axis categories
        y_label: str, column name of y axis values
        title: str, the title of the plot

    """
    # Plot style
    plt.close()
    plt.style.use(["seaborn-talk", "seaborn-whitegrid"])
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 16
    #plt.tight_layout()

    # Plot
    data = df
    x = range(len(data[y_label]))
    y = data[y_label]
    plt.bar(x, y)
    plt.xticks(x, data[x_label], rotation=90)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if len(title) > 70:
        plt.title(title[:70]+"...")
    else:
        plt.title(title)
    plt.tight_layout()



# =============================================================================
# UTILITY CLASSES
# =============================================================================

class IgnoredSubmissionError(Exception):
    """Exception used to signal a submission that must be ignored."""
    pass


class BadFormatError(Exception):
    """Exception used to signal a submission with unexpected formatting."""
    pass


class SamplSubmission:
    """A generic SAMPL submission.
    Parameters
    ----------
    file_path : str
        The path to the submission file.
    Raises
    ------
    IgnoredSubmission
        If the submission ID is among the ignored submissions.
    """
    # The D3R challenge IDs that are handled by this class.
    CHALLENGE_IDS = {1559}

    # The IDs of the submissions used for testing the validation.
    TEST_SUBMISSIONS = {}

    # Section of the submission file.
    SECTIONS = {}

    # Sections in CSV format with columns names.
    CSV_SECTIONS = {}

    def __init__(self, file_path, user_map):
        file_name = os.path.splitext(os.path.basename(file_path))[0]
        print(file_name)
        file_data = file_name.split('-')

        # Check if this is a deleted submission.
        if file_data[0] == 'DELETED':
            raise IgnoredSubmissionError('This submission was deleted.')

        # Check if this is a test submission.
        self.receipt_id = file_data[0]
        if self.receipt_id in self.TEST_SUBMISSIONS:
            raise IgnoredSubmissionError('This submission has been used for tests.')

        # Check this is the correct challenge.
        self.challenge_id = int(file_data[1])
        assert self.challenge_id in self.CHALLENGE_IDS

        # Store user map information.
        user_map_record = user_map[user_map.receipt_id == self.receipt_id]
        assert len(user_map_record) == 1
        user_map_record = user_map_record.iloc[0]

        self.id = user_map_record.id
        self.participant = user_map_record.firstname + ' ' + user_map_record.lastname
        self.participant_id = user_map_record.uid
        #self.participant_email = user_map_record.email
        #assert self.challenge_id == user_map_record.component

    @classmethod
    def _read_lines(cls, file_path):
        """Generator to read the file and discard blank lines and comments."""
        with open(file_path, 'r', encoding='utf-8-sig') as f:
            for line in f:
                # Strip whitespaces.
                line = line.strip()
                # Don't return blank lines and comments.
                if line != '' and line[0] != '#':
                    yield line

    @classmethod
    def _load_sections(cls, file_path):
        """Load the data in the file and separate it by sections."""
        sections = {}
        current_section = None
        for line in cls._read_lines(file_path):
            # Check if this is a new section.
            if line[:-1] in cls.SECTIONS:
                current_section = line[:-1]
            else:
                if current_section is None:
                    import pdb
                    pdb.set_trace()
                try:
                    sections[current_section].append(line)
                except KeyError:
                    sections[current_section] = [line]

        # Check that all the sections have been loaded.
        found_sections = set(sections.keys())
        if found_sections != cls.SECTIONS:
            raise BadFormatError('Missing sections: {}.'.format(found_sections - cls.SECTIONS))

        # Create a Pandas dataframe from the CSV format.
        for section_name in cls.CSV_SECTIONS:
            csv_str = io.StringIO('\n'.join(sections[section_name]))
            columns = cls.CSV_SECTIONS[section_name]
            id_column = columns[0]
            section = pd.read_csv(csv_str, index_col=id_column, names=columns, skipinitialspace=True)
            #section = pd.read_csv(csv_str, names=columns, skipinitialspace=True)
            sections[section_name] = section
        return sections

    @classmethod
    def _create_comparison_dataframe(cls, column_name, submission_data, experimental_data):
        """Create a single dataframe with submission and experimental data."""
        # Filter only the systems IDs in this submissions.


        experimental_data = experimental_data[experimental_data.index.isin(submission_data.index)] # match by column index
        # Fix the names of the columns for labelling.
        submission_series = submission_data[column_name]
        submission_series.name += ' (calc)'
        experimental_series = experimental_data[column_name]
        experimental_series.name += ' (expt)'

        # Concatenate the two columns into a single dataframe.
        return pd.concat([experimental_series, submission_series], axis=1)

# =============================================================================
# LOGP PREDICTION CHALLENGE
# =============================================================================

class logPSubmission(SamplSubmission):
    """A submission for logP challenge.

    Parameters
    ----------
    file_path : str
        The path to the submission file

    Raises
    ------
    IgnoredSubmission
        If the submission ID is among the ignored submissions.

    """

    # The D3R challenge IDs that are handled by this class.
    CHALLANGE_IDS = {1559}

    # The IDs of the submissions that will be ignored in the analysis.
    TEST_SUBMISSIONS = {}

    # Section of the submission file.
    SECTIONS = {'Predictions', 'Name', 'Software', 'Category', 'Method'}

    # Sections in CSV format with columns names.
    CSV_SECTIONS = {'Predictions': ("Molecule ID", "logP mean", "logP SEM", "logP model uncertainty")}


    def __init__(self, file_path, user_map):
        super().__init__(file_path, user_map)

        file_name = os.path.splitext(os.path.basename(file_path))[0]
        file_data = file_name.split('-')

        # Check if this is a type III submission
        self.submission_type = file_data[2]
        assert self.submission_type in ['logP']

        self.file_name, self.index = file_data[3:]
        self.index = int(self.index)

        # Load predictions.
        sections = self._load_sections(file_path)  # From parent-class.
        self.data = sections['Predictions']  # This is a pandas DataFrame.
        self.name = sections['Name'][0]
        self.category = sections['Category'][0] # New section for logP challenge.

    def compute_logP_statistics(self, experimental_data, stats_funcs):
        data = self._create_comparison_dataframe('logP mean', self.data, experimental_data)

        # Create lists of stats functions to pass to compute_bootstrap_statistics.
        stats_funcs_names, stats_funcs = zip(*stats_funcs.items())
        bootstrap_statistics = compute_bootstrap_statistics(data.as_matrix(), stats_funcs, n_bootstrap_samples=1000)

        # Return statistics as dict preserving the order.
        return collections.OrderedDict((stats_funcs_names[i], bootstrap_statistics[i])
                                      for i in range(len(stats_funcs)))

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================


def load_submissions(directory_path, user_map):
    submissions = []
    for file_path in glob.glob(os.path.join(directory_path, '*.csv')):
        try:
            submission = logPSubmission(file_path, user_map)

        except IgnoredSubmissionError:
            continue
        submissions.append(submission)
    return submissions



class logPSubmissionCollection:
    """A collection of logP submissions."""

    LOGP_CORRELATION_PLOT_BY_METHOD_PATH_DIR = 'logPCorrelationPlots'
    LOGP_CORRELATION_PLOT_WITH_SEM_BY_METHOD_PATH_DIR = 'logPCorrelationPlotsWithSEM'
    LOGP_CORRELATION_PLOT_BY_LOGP_PATH_DIR = 'error_for_each_logP.pdf'
    ABSOLUTE_ERROR_VS_LOGP_PLOT_PATH_DIR = 'AbsoluteErrorPlots'


    def __init__(self, submissions, experimental_data, output_directory_path, logP_submission_collection_file_path):


        # Check if submission collection file already exists.
        if os.path.isfile(logP_submission_collection_file_path):
            print("Analysis will be done using the existing logP_submission_collection.csv file.")

            self.data = pd.read_csv(logP_submission_collection_file_path)
            print("\n SubmissionCollection: \n")
            print(self.data)

            # Populate submission.data dataframes parsing sections of collection file.
            for submission in submissions:
                data = []

                receipt_ID = submission.receipt_id
                df_collection_of_each_submission = self.data.loc[self.data["receipt_id"] == receipt_ID ]

                # Transform into Pandas DataFrame.
                submission.data = pd.DataFrame()
                submission.data["logP mean"] = df_collection_of_each_submission["logP (calc)"]
                submission.data["logP SEM"] = df_collection_of_each_submission["logP SEM (calc)"]
                submission.data["Molecule ID"] = df_collection_of_each_submission["Molecule ID"]

                submission.data.set_index("Molecule ID", inplace=True)

            # Transform into Pandas DataFrame.
            self.output_directory_path = output_directory_path

        else: # Build collection dataframe from the beginning.
            # Build full logP collection table.
            data = []

            # Submissions for logP.
            for submission in submissions:

                print("submission.data:\n", submission.data)

                for mol_ID, series in submission.data.iterrows():
                    #print("mol_ID:", mol_ID)
                    #print("series:\n", series)

                    #mol_ID = series[1]["Molecule ID"]

                    #pKa_mean_exp = experimental_data.loc[experimental_data["pKa ID"] == pKa_ID, 'pKa mean'].values[0]
                    logP_mean_exp = experimental_data.loc[mol_ID, 'logP mean']
                    logP_SEM_exp = experimental_data.loc[mol_ID, 'logP SEM']

                    #pKa_mean_pred = submission.data.loc[submission.data["pKa ID"] == pKa_ID, 'pKa mean'].values[0]
                    logP_mean_pred = submission.data.loc[mol_ID, "logP mean"]
                    logP_SEM_pred = submission.data.loc[mol_ID, "logP SEM"]

                    data.append({
                        'receipt_id': submission.receipt_id,
                        'participant': submission.participant,
                        'name': submission.name,
                        'category': submission.category,
                        'Molecule ID': mol_ID,
                        'logP (calc)': logP_mean_pred,
                        'logP SEM (calc)': logP_SEM_pred,
                        'logP (exp)': logP_mean_exp,
                        'logP SEM (exp)': logP_SEM_exp,
                        '$\Delta$logP error (calc - exp)': logP_mean_pred - logP_mean_exp
                    })

            # Transform into Pandas DataFrame.
            self.data = pd.DataFrame(data=data)
            self.output_directory_path = output_directory_path

            print("\n SubmissionCollection: \n")
            print(self.data)

            # Create general output directory.
            os.makedirs(self.output_directory_path, exist_ok=True)

            # Save collection.data dataframe in a CSV file.
            self.data.to_csv(logP_submission_collection_file_path)

    def generate_correlation_plots(self):
        # logP correlation plots.
        output_dir_path = os.path.join(self.output_directory_path,
                                       self.LOGP_CORRELATION_PLOT_BY_METHOD_PATH_DIR)
        os.makedirs(output_dir_path, exist_ok=True)
        for receipt_id in self.data.receipt_id.unique():
            data = self.data[self.data.receipt_id == receipt_id]
            title = '{} ({})'.format(receipt_id, data.name.unique()[0])

            plt.close('all')
            plot_correlation(x='logP (exp)', y='logP (calc)',
                             data=data, title=title, kind='joint')
            plt.tight_layout()
            # plt.show()
            output_path = os.path.join(output_dir_path, '{}.pdf'.format(receipt_id))
            plt.savefig(output_path)

    def generate_correlation_plots_with_SEM(self):
        # logP correlation plots.
        output_dir_path = os.path.join(self.output_directory_path,
                                       self.LOGP_CORRELATION_PLOT_WITH_SEM_BY_METHOD_PATH_DIR)
        os.makedirs(output_dir_path, exist_ok=True)
        for receipt_id in self.data.receipt_id.unique():
            data = self.data[self.data.receipt_id == receipt_id]
            title = '{} ({})'.format(receipt_id, data.name.unique()[0])

            plt.close('all')
            plot_correlation_with_SEM(x_lab='logP (exp)', y_lab='logP (calc)',
                                      x_err_lab='logP SEM (exp)', y_err_lab='logP SEM (calc)',
                                      data=data, title=title)
            plt.tight_layout()
            # plt.show()
            output_path = os.path.join(output_dir_path, '{}.pdf'.format(receipt_id))
            plt.savefig(output_path)

    def generate_molecules_plot(self):
        # Correlation plot by molecules.
        plt.close('all')
        data_ordered_by_mol_ID = self.data.sort_values(["Molecule ID"], ascending=["True"])
        sns.set(rc={'figure.figsize': (8.27,11.7)})
        sns.violinplot(y='Molecule ID', x='$\Delta$logP error (calc - exp)', data=data_ordered_by_mol_ID,
                           inner='point', linewidth=1, width=1.2)
        plt.tight_layout()
        # plt.show()
        plt.savefig(os.path.join(self.output_directory_path, self.LOGP_CORRELATION_PLOT_BY_LOGP_PATH_DIR))

    def generate_absolute_error_vs_molecule_ID_plot(self):
        """
        For each method a bar plot is generated so that absolute errors of each molecule can be compared.
        """
        # Setup output directory
        output_dir_path = os.path.join(self.output_directory_path,
                                       self.ABSOLUTE_ERROR_VS_LOGP_PLOT_PATH_DIR)
        os.makedirs(output_dir_path, exist_ok=True)

        # Calculate absolute errors.
        self.data["absolute error"] = np.NaN
        self.data.loc[:, "absolute error"] = np.absolute(self.data.loc[:, "$\Delta$logP error (calc - exp)"])

        # Create a separate plot for each submission.
        for receipt_id in self.data.receipt_id.unique():
            data = self.data[self.data.receipt_id == receipt_id]
            title = '{} ({})'.format(receipt_id, data.name.unique()[0])

            plt.close('all')
            barplot(df=data, x_label="Molecule ID", y_label="absolute error", title=title)
            output_path = os.path.join(output_dir_path, '{}.pdf'.format(receipt_id))
            plt.savefig(output_path)


def generate_statistics_tables(submissions, stats_funcs, directory_path, file_base_name,
                                sort_stat=None, ordering_functions=None,
                                latex_header_conversions=None):
    stats_names = list(stats_funcs.keys())
    ci_suffixes = ('', '_lower_bound', '_upper_bound')

    # Collect the records for the DataFrames.
    statistics_csv = []
    statistics_latex = []
    statistics_plot = []

    for i, submission in enumerate(submissions):
        receipt_id = submission.receipt_id
        category = submission.category

        print('\rGenerating bootstrap statistics for submission {} ({}/{})'
                  ''.format(receipt_id, i + 1, len(submissions)), end='')

        bootstrap_statistics = submission.compute_logP_statistics(experimental_data, stats_funcs)

        record_csv = {}
        record_latex = {}
        for stats_name, (stats, (lower_bound, upper_bound), bootstrap_samples) in bootstrap_statistics.items():
            # For CSV and JSON we put confidence interval in separate columns.
            for suffix, info in zip(ci_suffixes, [stats, lower_bound, upper_bound]):
                record_csv[stats_name + suffix] = info

            # For the PDF, print bootstrap CI in the same column.
            stats_name_latex = latex_header_conversions.get(stats_name, stats_name)
            record_latex[stats_name_latex] = '{:.3f} [{:.3f}, {:.3f}]'.format(stats, lower_bound, upper_bound)

            # For the violin plot, we need all the bootstrap statistics series.
            for bootstrap_sample in bootstrap_samples:
                statistics_plot.append(dict(ID=receipt_id, name=submission.name, category=category,
                                            statistics=stats_name_latex, value=bootstrap_sample))

        statistics_csv.append({'ID': receipt_id, 'name': submission.name, 'category': category, **record_csv})
        escaped_name = submission.name.replace('_', '\_')
        statistics_latex.append({'ID': receipt_id, 'name': escaped_name, 'category': category, **record_latex})
    print()
    print("statistics_csv:\n",statistics_csv)

    # Convert dictionary to Dataframe to create tables/plots easily.
    statistics_csv = pd.DataFrame(statistics_csv)
    #print("statistics_csv:\n", statistics_csv)
    #print("statistics_csv.columns:\n", statistics_csv.columns)
    statistics_csv.set_index('ID', inplace=True)
    statistics_latex = pd.DataFrame(statistics_latex)
    statistics_plot = pd.DataFrame(statistics_plot)

    # Sort by the given statistics.
    if sort_stat is not None:
        statistics_csv.sort_values(by=sort_stat, inplace=True)
        statistics_latex.sort_values(by=latex_header_conversions.get(sort_stat, sort_stat),
                                     inplace=True)

    # Reorder columns that were scrambled by going through a dictionaries.
    stats_names_csv = [name + suffix for name in stats_names for suffix in ci_suffixes]
    #print("stats_names_csv:", stats_names_csv)
    stats_names_latex = [latex_header_conversions.get(name, name) for name in stats_names]
    #print("stats_names_latex:", stats_names_latex)
    statistics_csv = statistics_csv[['name', "category"] + stats_names_csv]
    statistics_latex = statistics_latex[['ID', 'name'] + stats_names_latex]

    # Create CSV and JSON tables (correct LaTex syntax in column names).
    os.makedirs(directory_path)
    file_base_path = os.path.join(directory_path, file_base_name)
    with open(file_base_path + '.csv', 'w') as f:
        statistics_csv.to_csv(f)
    with open(file_base_path + '.json', 'w') as f:
        statistics_csv.to_json(f, orient='index')

    # Create LaTex table.
    latex_directory_path = os.path.join(directory_path, file_base_name + 'LaTex')
    os.makedirs(latex_directory_path, exist_ok=True)
    with open(os.path.join(latex_directory_path, file_base_name + '.tex'), 'w') as f:
        f.write('\\documentclass{article}\n'
                '\\usepackage[a4paper,margin=0.005in,tmargin=0.5in,landscape]{geometry}\n'
                '\\usepackage{booktabs}\n'
                '\\usepackage{longtable}\n'
                '\\pagenumbering{gobble}\n'
                '\\begin{document}\n'
                '\\begin{center}\n')
        statistics_latex.to_latex(f, column_format='|ccccccc|', escape=False, index=False, longtable=True)
        f.write('\end{center}\n' 
                '\nNotes\n\n'
                '- Mean and 95\% confidence intervals of statistic values were calculated by bootstrapping.\n\n'
                'pKas of the rest of the molecules in these submissions were blindly predicted before experimental data was released.\n\n'
                #'- Some logP predictions were submitted after the submission deadline to be used as a reference method.\n\n'
                '\end{document}\n')

    # Violin plots by statistics across submissions.
    plt.close('all')
    fig, axes = plt.subplots(ncols=len(stats_names), figsize=(12, 0.375 * len(submissions)))
    for ax, stats_name in zip(axes, stats_names):
        stats_name_latex = latex_header_conversions.get(stats_name, stats_name)
        data = statistics_plot[statistics_plot.statistics == stats_name_latex]
        # Plot ordering submission by statistics.
        ordering_function = ordering_functions.get(stats_name, lambda x: x)
        order = sorted(statistics_csv[stats_name].items(), key=lambda x: ordering_function(x[1]))
        order = [receipt_id for receipt_id, value in order]
        sns.violinplot(x='value', y='ID', data=data, ax=ax,
                        order=order, palette='PuBuGn_r', linewidth=0.5, width=1)
        ax.set_xlabel(stats_name_latex)
        ax.set_ylabel('')
        sns.set_style("whitegrid")
    plt.tight_layout()
    # plt.show()
    plt.savefig(file_base_path + '_bootstrap_distributions.pdf')

def generate_performance_comparison_plots(statistics_filename, directory_path):
        # Read statistics table
        statistics_file_path = os.path.join(directory_path, statistics_filename)
        df_statistics = pd.read_csv(statistics_file_path)
        #print("\n df_statistics \n", df_statistics)

        # RMSE comparison plot
        barplot_with_CI_errorbars(df=df_statistics, x_label="ID", y_label="RMSE", y_lower_label="RMSE_lower_bound",
                                  y_upper_label="RMSE_upper_bound", figsize=(22,10))
        plt.savefig(directory_path + "/RMSE_vs_method_plot.pdf")

        # RMSE comparison plot with each category colored separately
        barplot_with_CI_errorbars_colored_by_label(df=df_statistics, x_label="ID", y_label="RMSE",
                                  y_lower_label="RMSE_lower_bound",
                                  y_upper_label="RMSE_upper_bound", color_label = "category", figsize=(22,10))
        plt.ylim(0.0, 7.0)
        plt.savefig(directory_path + "/RMSE_vs_method_plot_colored_by_method_category.pdf")

        # MAE comparison plot
        # Reorder based on MAE value
        df_statistics_MAE = df_statistics.sort_values(by="MAE", inplace=False)

        barplot_with_CI_errorbars(df=df_statistics_MAE, x_label="ID", y_label="MAE", y_lower_label="MAE_lower_bound",
                                  y_upper_label="MAE_upper_bound", figsize=(22,10))
        plt.savefig(directory_path + "/MAE_vs_method_plot.pdf")

        # MAE comparison plot with each category colored separately
        barplot_with_CI_errorbars_colored_by_label(df=df_statistics_MAE, x_label="ID", y_label="MAE",
                                                   y_lower_label="MAE_lower_bound",
                                                   y_upper_label="MAE_upper_bound", color_label="category",
                                                   figsize=(22, 10))
        plt.ylim(0.0, 7.0)
        plt.savefig(directory_path + "/MAE_vs_method_plot_colored_by_method_category.pdf")


        # Plot RMSE and MAE comparison plots for each category separately
        category_list = ["Physical","Empirical", "Mixed", "Other"]

        for category in category_list:
            print("category: ",category)
            #print("df_statistics.columns:\n", df_statistics.columns)

            # Take subsection of dataframe for each category
            df_statistics_1category = df_statistics.loc[df_statistics['category'] == category]
            df_statistics_MAE_1category = df_statistics_MAE.loc[df_statistics_MAE['category'] == category]


            # RMSE comparison plot for each category
            barplot_with_CI_errorbars(df=df_statistics_1category, x_label="ID", y_label="RMSE", y_lower_label="RMSE_lower_bound",
                                      y_upper_label="RMSE_upper_bound", figsize=(12, 10))
            plt.title("Method category: {}".format(category), fontdict={'fontsize': 22})
            plt.ylim(0.0,7.0)
            plt.savefig(directory_path + "/RMSE_vs_method_plot_for_{}_category.pdf".format(category))

            # MAE comparison plot for each category
            barplot_with_CI_errorbars(df=df_statistics_MAE_1category, x_label="ID", y_label="MAE",
                                      y_lower_label="MAE_lower_bound",
                                      y_upper_label="MAE_upper_bound", figsize=(12, 10))
            plt.title("Method category: {}".format(category), fontdict={'fontsize': 22})
            plt.ylim(0.0, 7.0)
            plt.savefig(directory_path + "/MAE_vs_method_plot_for_{}_category.pdf".format(category))




# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':

    sns.set_style('whitegrid')
    sns.set_context('paper')

    # Read experimental data.
    with open(EXPERIMENTAL_DATA_FILE_PATH, 'r') as f:
        # experimental_data = pd.read_json(f, orient='index')
        names = ('Molecule ID', 'logP mean', 'logP SEM',
                 'Assay Type', 'Experimental ID', 'Isomeric SMILES')
        experimental_data = pd.read_csv(f, names=names, skiprows=1)

    # Convert numeric values to dtype float.
    for col in experimental_data.columns[1:7]:
        experimental_data[col] = pd.to_numeric(experimental_data[col], errors='coerce')


    experimental_data.set_index("Molecule ID", inplace=True)
    experimental_data["Molecule ID"] = experimental_data.index
    print("Experimental data: \n", experimental_data)

    # Import user map.
    with open('../predictions/SAMPL6-user-map-logP.csv', 'r') as f:
        user_map = pd.read_csv(f)

    # Configuration: statistics to compute.
    stats_funcs = collections.OrderedDict([
        ('RMSE', rmse),
        ('MAE', mae),
        ('ME', me),
        ('R2', r2),
        ('m', slope),
    ])
    ordering_functions = {
        'ME': lambda x: abs(x),
        'R2': lambda x: -x,
        'm': lambda x: abs(1 - x),
    }
    latex_header_conversions = {
        'R2': 'R$^2$',
        'RMSE': 'RMSE',
        'MAE': 'MAE',
        'ME': 'ME',
    }

    # Load submissions data.
    submissions_logP = load_submissions(LOGP_SUBMISSIONS_DIR_PATH, user_map)

    # Perform the analysis using the different algorithms for matching predictions to experiment


    output_directory_path='./analysis_outputs'
    logP_submission_collection_file_path = '{}/logP_submission_collection.csv'.format(output_directory_path)

    collection_logP= logPSubmissionCollection(submissions_logP, experimental_data,
                                             output_directory_path, logP_submission_collection_file_path)

    # Generate plots and tables.
    for collection in [collection_logP]:
        collection.generate_correlation_plots()
        collection.generate_correlation_plots_with_SEM()
        collection.generate_molecules_plot()
        collection.generate_absolute_error_vs_molecule_ID_plot()

    import shutil

    if os.path.isdir('{}/StatisticsTables'.format(output_directory_path)):
        shutil.rmtree('{}/StatisticsTables'.format(output_directory_path))


    for submissions, type in zip([submissions_logP], ['logP']):
        generate_statistics_tables(submissions, stats_funcs, directory_path=output_directory_path + '/StatisticsTables',
                                    file_base_name='statistics', sort_stat='RMSE',
                                    ordering_functions=ordering_functions,
                                    latex_header_conversions=latex_header_conversions)

    # Generate RMSE and MAE comparison plots.
    statistics_directory_path = os.path.join(output_directory_path, "StatisticsTables")
    generate_performance_comparison_plots(statistics_filename="statistics.csv", directory_path=statistics_directory_path)
