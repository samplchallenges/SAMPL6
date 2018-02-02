#!/usr/bin/env python

# =============================================================================
# GLOBAL IMPORTS
# =============================================================================

import os
import io
import glob
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
        bootstrap_statistics.append([statistics[stats_func_idx], confidence_interval])

    return bootstrap_statistics


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

    def __init__(self, file_path):
        file_name = os.path.splitext(os.path.basename(file_path))[0]
        file_data = file_name.split('-')

        # Check if this is a deleted submission.
        if file_data[0] == 'DELETED':
            raise IgnoredSubmissionError('This submission was deleted.')

        # Check if this is a test submission.
        self.submission_id = file_data[0]
        if self.submission_id in self.TEST_SUBMISSIONS:
            raise IgnoredSubmissionError('This submission has been used for tests.')

        # Check this is the correct challenge.
        self.challenge_id = int(file_data[1])
        assert self.challenge_id in self.CHALLENGE_IDS

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

    def __init__(self, file_path):
        super().__init__(file_path)

        file_name = os.path.splitext(os.path.basename(file_path))[0]
        file_data = file_name.split('-')

        # Check this is a known host, and tansform into uppercase.
        self.host_name = file_data[2].upper()
        assert self.host_name in ['OA', 'TEMOA', 'CB8']

        self.file_name, self.index = file_data[3:]
        self.index = int(self.index)

        # Load predictions.
        sections = self._load_sections(file_path)  # From parent-class.
        self.predictions = sections['Predictions']  # This is a pandas DataFrame.
        self.name = sections['Name'][0]

    def compute_free_energy_statistics(self, experimental_data, stats_funcs):
        data = self._create_comparison_dataframe('$\Delta$G', self.predictions, experimental_data)
        # Create lists of stats functions to pass to function.
        stats_funcs_names, stats_funcs = zip(*stats_funcs.items())
        bootstrap_statistics = compute_bootstrap_statistics(data.as_matrix(), stats_funcs)
        return {stats_funcs_names[i]: bootstrap_statistics[i] for i in range(len(stats_funcs))}

    def plot_free_energy_correlation(self, experimental_data, axes_limits='equal'):
        """Generate a correlation plot of the free energies."""
        # Create dataframe for seaborn.
        data = self._create_comparison_dataframe('$\Delta$G', self.predictions, experimental_data)

        # Take care of equal axes limits.
        if axes_limits == 'equal':
            # Find extreme limits.
            min_limit = np.ceil(min(data.min()) - 2)
            max_limit = np.floor(max(data.max()) + 2)
            axes_limits = np.array([min_limit, max_limit])

        grid = sns.jointplot(x='$\Delta$G (expt)', y='$\Delta$G (calc)', data=data,
                             kind='reg', joint_kws={'ci': None}, stat_func=None,
                             xlim=axes_limits, ylim=axes_limits)

        # Title.
        grid.fig.subplots_adjust(top=0.95)
        grid.fig.suptitle('{} ({})'.format(self.name, self.submission_id))

        # Add diagonal line.
        grid.ax_joint.plot(axes_limits, axes_limits, ls='--', c='black', alpha=0.8, lw=0.7)

        # Add shaded area for 1-2 kcal/mol error.
        palette = sns.color_palette('BuGn_r')
        grid.ax_joint.fill_between(axes_limits, axes_limits - 1, axes_limits + 1, alpha=0.2, color=palette[2])
        grid.ax_joint.fill_between(axes_limits, axes_limits - 2, axes_limits + 2, alpha=0.2, color=palette[3])


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def load_submissions(output_directory_path):
    submissions = []
    for file_path in glob.glob(os.path.join(output_directory_path, '*.txt')):
        try:
            submission = HostGuestSubmission(file_path)
        except IgnoredSubmissionError:
            continue
        submissions.append(submission)
    return submissions


def generate_correlation_plots(submissions, output_directory_path):
    os.makedirs(output_directory_path, exist_ok=True)
    for submission in submissions:
        plt.close('all')
        submission.plot_free_energy_correlation(experimental_data)
        # plt.show()
        submission_id = submission.submission_id
        plt.savefig(os.path.join(output_directory_path, '{}.pdf'.format(submission_id)))


def generate_statistics_tables(submissions, stats_funcs, directory_path, file_base_name):
    statistics = {}
    for i, submission in enumerate(submissions):
        submission_id = submission.submission_id
        print('\rGenerating bootstrap statistics for submission {} ({}/{})'
              ''.format(submission_id, i+1, len(submissions)), end='')

        statistics[submission_id] = collections.OrderedDict()
        bootstrap_statistics = submission.compute_free_energy_statistics(experimental_data, stats_funcs)

        # Separate statistics and its confidence interval into separate columns.
        for stats_name, (stats, confidence_interval) in bootstrap_statistics.items():
            statistics[submission_id][stats_name] = stats
            statistics[submission_id][stats_name + '_lower_bound'] = confidence_interval[0]
            statistics[submission_id][stats_name + '_upper_bound'] = confidence_interval[1]

        if len(statistics) == 2:
            break
    print()

    # Convert to DataFrame.
    statistics = pd.DataFrame.from_dict(statistics, orient='index')
    statistics.index.names = ['ID']

    # Create LaTex table.
    latex_directory_path = os.path.join(directory_path, 'LaTex')
    os.makedirs(latex_directory_path, exist_ok=True)
    with open(os.path.join(latex_directory_path, file_base_name + '.tex'), 'w') as f:
        f.write('\\documentclass{article}\n'
                '\\usepackage{booktabs}\n'
                '\\pagenumbering{gobble}\n'
                '\\begin{document}\n'
                '\\begin{center}\n')
        statistics.to_latex(f)
        f.write('\end{center}\n'
                '\end{document}\n')

    # Create CSV and JSON tables (correct LaTex syntax in column names).
    statistics.columns = [name.replace('$^2$', '2') for name in statistics.columns]
    file_base_path = os.path.join(directory_path, file_base_name)
    with open(file_base_path + '.csv', 'w') as f:
        statistics.to_csv(f, header=statistics.columns)
    with open(file_base_path + '.json', 'w') as f:
        statistics.to_json(f, orient='index')


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    sns.set_style('whitegrid')
    sns.set_context('talk')

    # Read experimental data.
    with open(EXPERIMENTAL_DATA_FILE_PATH, 'r') as f:
        # experimental_data = pd.read_json(f, orient='index')
        names = ('System ID', 'name', 'SMILES', 'Ka', 'dKa', '$\Delta$H', 'd$\Delta$H',
                 'T$\Delta$S', 'dT$\Delta$S', 'n', '$\Delta$G', 'd$\Delta$G')
        experimental_data = pd.read_csv(f, sep=';', names=names, index_col='System ID', skiprows=0)

    # Convert numeric values to dtype float.
    # experimental_data = experimental_data.convert_objects(convert_numeric=True)
    for col in experimental_data.columns[3:]:
        experimental_data[col] = pd.to_numeric(experimental_data[col], errors='coerce')

    # TODO:
    # TODO:     I had to fix the index CB8-G12a -> CB8-G12 to make the analysis work
    # TODO:     ../Submissions/974/tb3ck-974-CB8-WGatMSU-1.txt: has an extra - in CB8-G6 enthalpy

    # Load submissions data.
    submissions = load_submissions(HOST_GUEST_CB_SUBMISSIONS_DIR_PATH)

    # Correlation plots.
    # generate_correlation_plots(submissions, 'CorrelationPlots/CB8/')

    # Statistics tables.
    stats_funcs = collections.OrderedDict([
        ('R$^2$', r2),
        ('RMSE', rmse),
        ('m', slope),
    ])
    generate_statistics_tables(submissions, stats_funcs, file_base_path='StatisticsTables/CB8')


    # Bootstrap distributions plots.
