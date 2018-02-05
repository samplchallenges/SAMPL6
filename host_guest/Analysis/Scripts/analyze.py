#!/usr/bin/env python

# =============================================================================
# GLOBAL IMPORTS
# =============================================================================

import os
import io
import glob
import math
import collections

import numpy as np
import scipy.stats
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


# =============================================================================
# CONSTANTS
# =============================================================================

# Paths to input data.
HOST_GUEST_OA_SUBMISSIONS_DIR_PATH = '../Submissions/973/'
HOST_GUEST_CB_SUBMISSIONS_DIR_PATH = '../Submissions/974/'
EXPERIMENTAL_DATA_FILE_PATH = '../ExperimentalMeasurements/experimental_measurements.csv'


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


def compute_bootstrap_statistics(samples, stats_funcs, percentile=0.95, n_bootstrap_samples=10000):
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
    # Extract only free energies.
    data = data[[x, y]]

    # Find extreme values to make axes equal.
    min_limit = np.ceil(min(data.min()) - 2)
    max_limit = np.floor(max(data.max()) + 2)
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

    # Add shaded area for 1-2 kcal/mol error.
    palette = sns.color_palette('BuGn_r')
    ax.fill_between(axes_limits, axes_limits - 1, axes_limits + 1, alpha=0.2, color=palette[2])
    ax.fill_between(axes_limits, axes_limits - 2, axes_limits + 2, alpha=0.2, color=palette[3])


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
    CHALLENGE_IDS = {}

    # The IDs of the submissions used for testing the validation.
    TEST_SUBMISSIONS = {}

    # Section of the submission file.
    SECTIONS = {}

    # Sections in CSV format with columns names.
    CSV_SECTIONS = {}

    def __init__(self, file_path, user_map):
        file_name = os.path.splitext(os.path.basename(file_path))[0]
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
        self.participant_email = user_map_record.email
        assert self.challenge_id == user_map_record.component

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
            section = pd.read_csv(csv_str, index_col=id_column, names=columns)
            sections[section_name] = section
        return sections

    @classmethod
    def _create_comparison_dataframe(cls, column_name, submission_data, experimental_data):
        """Create a single dataframe with submission and experimental data."""
        # Filter only the systems IDs in this submissions.
        experimental_data = experimental_data[experimental_data.index.isin(submission_data.index)]
        # Fix the names of the columns for labelling.
        submission_series = submission_data[column_name]
        submission_series.name += ' (calc)'
        experimental_series = experimental_data[column_name]
        experimental_series.name += ' (expt)'
        # Concatenate the two columns into a single dataframe.
        return pd.concat([submission_series, experimental_series], axis=1)


# =============================================================================
# MAIN CHALLENGE HOST-GUEST SUBMISSION
# =============================================================================

class HostGuestSubmission(SamplSubmission):
    """A submission for the main host-guest challenge.

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
    CHALLENGE_IDS = {973, 974}

    # The IDs of the submissions used for testing the validation.
    TEST_SUBMISSIONS = {'egte7', '6odqw', 'q6igp', '2fpb5',
                        'dofcg', 'iw8mj', '5dbnp', 'exb60'}

    # Section of the submission file.
    SECTIONS = {'Predictions', 'Name', 'Software', 'Method'}

    # Sections in CSV format with columns names.
    CSV_SECTIONS = {'Predictions': ('System ID', '$\Delta$G', 'SEM $\Delta$G', 'd$\Delta$G',
                                    '$\Delta$H', 'SEM $\Delta$H', 'd$\Delta$H')}

    def __init__(self, file_path, user_map):
        super().__init__(file_path, user_map)

        file_name = os.path.splitext(os.path.basename(file_path))[0]
        file_data = file_name.split('-')

        # Check this is a known host, and tansform into uppercase.
        self.host_name = file_data[2].upper()
        assert self.host_name in ['OA', 'TEMOA', 'CB8']

        self.file_name, self.index = file_data[3:]
        self.index = int(self.index)

        # Load predictions.
        sections = self._load_sections(file_path)  # From parent-class.
        self.data = sections['Predictions']  # This is a pandas DataFrame.
        self.name = sections['Name'][0]

    def compute_free_energy_statistics(self, experimental_data, stats_funcs):
        data = self._create_comparison_dataframe('$\Delta$G', self.data, experimental_data)

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
    for file_path in glob.glob(os.path.join(directory_path, '*.txt')):
        try:
            submission = HostGuestSubmission(file_path, user_map)
        except IgnoredSubmissionError:
            continue
        submissions.append(submission)
    return submissions


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
        print('\rGenerating bootstrap statistics for submission {} ({}/{})'
              ''.format(receipt_id, i+1, len(submissions)), end='')

        bootstrap_statistics = submission.compute_free_energy_statistics(experimental_data, stats_funcs)

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
                statistics_plot.append(dict(ID=receipt_id, name=submission.name,
                                            statistics=stats_name_latex, value=bootstrap_sample))

        statistics_csv.append({'ID': receipt_id, 'name': submission.name, **record_csv})
        escaped_name = submission.name.replace('_', '\_')
        statistics_latex.append({'ID': receipt_id, 'name': escaped_name, **record_latex})
    print()

    # Convert dictionary to Dataframe to create tables/plots easily.
    statistics_csv = pd.DataFrame(statistics_csv)
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
    stats_names_latex = [latex_header_conversions.get(name, name) for name in stats_names]
    statistics_csv = statistics_csv[['name'] + stats_names_csv]
    statistics_latex = statistics_latex[['ID', 'name'] + stats_names_latex]

    # Create CSV and JSON tables (correct LaTex syntax in column names).
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
                '\\usepackage[a4paper,margin=0.4in,tmargin=0.5in,landscape]{geometry}\n'
                '\\usepackage{booktabs}\n'
                '\\pagenumbering{gobble}\n'
                '\\begin{document}\n'
                '\\begin{center}\n')
        statistics_latex.to_latex(f, column_format='|ccccccc|', escape=False)
        f.write('\end{center}\n'
                '\end{document}\n')

    # Violin plots by statistics across submissions.
    plt.close('all')
    fig, axes = plt.subplots(ncols=len(stats_names), figsize=(12, 0.375*len(submissions)))
    for ax, stats_name in zip(axes, stats_names):
        stats_name_latex = latex_header_conversions.get(stats_name, stats_name)
        data = statistics_plot[statistics_plot.statistics == stats_name_latex]
        # Plot ordering submission by statistics.
        ordering_function = ordering_functions.get(stats_name, lambda x: x)
        order = sorted(statistics_csv[stats_name].items(), key=lambda x: ordering_function(x[1]))
        order = [receipt_id for receipt_id, value in order]
        sns.violinplot(x='value', y='ID', data=data, ax=ax,
                       order=order, palette='PuBuGn_r')
        ax.set_xlabel(stats_name_latex)
        ax.set_ylabel('')
    plt.tight_layout()
    # plt.show()
    plt.savefig(file_base_path + '_bootstrap.pdf')


class HostGuestSubmissionCollection:
    """A collection of HostGuestSubmissions."""

    SUBMISSION_CORRELATION_PLOT_DIR = 'SubmissionsCorrelationPlots'
    MOLECULE_CORRELATION_PLOT_PATH = 'molecules_error.pdf'

    def __init__(self, submissions, experimental_data, output_directory_path):
        # Build full free energy table.
        data = []

        # Submissions free energies.
        for submission in submissions:
            for system_id, free_energy_calc in submission.data['$\Delta$G'].items():
                free_energy_expt = experimental_data.loc[system_id, '$\Delta$G']
                data.append({
                    'receipt_id': submission.receipt_id,
                    'participant': submission.participant,
                    'name': submission.name,
                    'system_id': system_id,
                    '$\Delta$G (calc)': free_energy_calc,
                    '$\Delta$G (expt)': free_energy_expt,
                    '$\Delta\Delta$G': free_energy_calc - free_energy_expt
                })

        # Transform into Pandas DataFrame.
        self.data = pd.DataFrame(data=data)
        self.output_directory_path = output_directory_path

        # Create general output directory.
        os.makedirs(self.output_directory_path, exist_ok=True)

    def generate_correlation_plots(self):
        # Correlation plot by submission.
        output_dir_path = os.path.join(self.output_directory_path,
                                       self.SUBMISSION_CORRELATION_PLOT_DIR)
        os.makedirs(output_dir_path, exist_ok=True)
        for receipt_id in self.data.receipt_id.unique():
            data = self.data[self.data.receipt_id == receipt_id]
            title = '{} ({})'.format(receipt_id, data.name.unique()[0])

            plt.close('all')
            plot_correlation(x='$\Delta$G (expt)', y='$\Delta$G (calc)',
                             data=data, title=title, kind='joint')
            plt.tight_layout()
            # plt.show()
            output_path = os.path.join(output_dir_path, '{}.pdf'.format(receipt_id))
            plt.savefig(output_path)

    def generate_molecules_plot(self):
        # Correlation plot by molecules.
        plt.close('all')
        sns.violinplot(y='system_id', x='$\Delta\Delta$G', data=self.data, inner='point')
        plt.tight_layout()
        # plt.show()
        plt.savefig(os.path.join(self.output_directory_path, self.MOLECULE_CORRELATION_PLOT_PATH))


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    # TODO:     I had to fix the index CB8-G12a -> CB8-G12 to make the analysis work
    # TODO:     ../Submissions/974/tb3ck-974-CB8-WGatMSU-1.txt: has an extra - in CB8-G6 enthalpy

    sns.set_style('whitegrid')
    sns.set_context('paper')

    # Read experimental data.
    with open(EXPERIMENTAL_DATA_FILE_PATH, 'r') as f:
        # experimental_data = pd.read_json(f, orient='index')
        names = ('System ID', 'name', 'SMILES', 'Ka', 'dKa', '$\Delta$H', 'd$\Delta$H',
                 'T$\Delta$S', 'dT$\Delta$S', 'n', '$\Delta$G', 'd$\Delta$G')
        experimental_data = pd.read_csv(f, sep=';', names=names, index_col='System ID', skiprows=1)

    # Convert numeric values to dtype float.
    # experimental_data = experimental_data.convert_objects(convert_numeric=True)
    for col in experimental_data.columns[3:]:
        experimental_data[col] = pd.to_numeric(experimental_data[col], errors='coerce')

    # Import user map.
    with open('../Submissions/SAMPL6_user_map.csv', 'r') as f:
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
        'RMSE': 'RMSE (kcal/mol)',
        'MAE': 'MAE (kcal/mol)',
        'ME': 'ME (kcal/mol)',
    }

    # Load submissions data. We do OA and TEMOA together.
    submissions_cb = load_submissions(HOST_GUEST_CB_SUBMISSIONS_DIR_PATH, user_map)
    collection_cb = HostGuestSubmissionCollection(submissions_cb, experimental_data,
                                                  output_directory_path='CB8')

    submissions_oa = load_submissions(HOST_GUEST_OA_SUBMISSIONS_DIR_PATH, user_map)
    collection_oa = HostGuestSubmissionCollection(submissions_oa, experimental_data,
                                                  output_directory_path='OA')

    for collection in [collection_cb, collection_oa]:
        collection.generate_correlation_plots()
        collection.generate_molecules_plot()

    # Generate plots and tables.
    for submissions, host in zip([submissions_cb, submissions_oa], ['CB8', 'OA']):
        generate_statistics_tables(submissions, stats_funcs, directory_path='StatisticsTables',
                                   file_base_name=host, sort_stat='RMSE', ordering_functions=ordering_functions,
                                   latex_header_conversions=latex_header_conversions)
