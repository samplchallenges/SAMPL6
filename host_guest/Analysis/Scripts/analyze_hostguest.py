#!/usr/bin/env python

# =============================================================================
# GLOBAL IMPORTS
# =============================================================================

import os
import copy
import collections

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from pkganalysis.submission import (SamplSubmission, IgnoredSubmissionError,
                                    load_submissions, plot_correlation)
from pkganalysis.stats import (compute_bootstrap_statistics, rmse, mae,
                               me, r2, slope, kendall_tau)


# =============================================================================
# CONSTANTS
# =============================================================================

# Paths to input data.
HOST_GUEST_OA_SUBMISSIONS_DIR_PATH = '../SubmissionsDoNotUpload/973/'
HOST_GUEST_CB_SUBMISSIONS_DIR_PATH = '../SubmissionsDoNotUpload/974/'
EXPERIMENTAL_DATA_FILE_PATH = '../ExperimentalMeasurements/experimental_measurements.csv'


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
                        'dofcg', 'iw8mj', '5dbnp', 'exb60',
                        'gjdze', 'hkgxh', 'd7xde'}

    # Section of the submission file.
    SECTIONS = {'Predictions', 'Name', 'Software', 'Method'}

    # Sections in CSV format with kwargs to pass to pandas.read_csv().
    CSV_SECTIONS = {
        'Predictions': {'names': ('System ID', '$\Delta$G', 'SEM $\Delta$G', 'd$\Delta$G',
                                  '$\Delta$H', 'SEM $\Delta$H', 'd$\Delta$H'),
                        'index_col': 'System ID'}
    }

    RENAME_METHODS = {
        'Alchemical Free Energy Calculations': 'SOMD/AM1BCC-GAFF-TIP3P/MBAR/A'
    }

    def __init__(self, file_path, user_map):
        super().__init__(file_path, user_map)

        file_name = os.path.splitext(os.path.basename(file_path))[0]
        file_data = file_name.split('-')

        # Check this is a known host, and tansform into uppercase.
        self.host_name = file_data[2].upper()
        assert self.host_name in ['OA', 'TEMOA', 'CB8']

        self.file_name, self.index = file_data[3:]

        # Load predictions.
        sections = self._load_sections(file_path)  # From parent-class.
        self.data = sections['Predictions']  # This is a pandas DataFrame.
        try:
            self.name = self.RENAME_METHODS[sections['Name'][0]]
        except KeyError:
            self.name = sections['Name'][0]

        # Add host name column to predictions.
        self.data['host_name'] = self.host_name

    def __add__(self, other):
        """Merge the data of the two submission."""
        merged_submission = copy.deepcopy(self)
        merged_submission.receipt_id = '{} + {}'.format(*sorted([self.receipt_id, other.receipt_id]))
        merged_submission.index = '{} + {}'.format(*sorted([self.index, other.index]))
        merged_submission.host_name = '{} + {}'.format(*sorted([self.host_name, other.host_name]))

        # Check if this is already a merged submission.
        if isinstance(merged_submission.file_name, list):
            merged_submission.file_name = sorted([*merged_submission.file_name, other.file_name])
        else:
            merged_submission.file_name = sorted([merged_submission.file_name, other.file_name])
        merged_submission.data = pd.concat([merged_submission.data, other.data])
        return merged_submission


# =============================================================================
# NULL METHODS
# =============================================================================




# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

class HostGuestSubmissionCollection:
    """A collection of HostGuestSubmissions."""

    FREE_ENERGY_CORRELATION_PLOT_DIR = 'FreeEnergyCorrelationPlots'
    ENTHALPIES_CORRELATION_PLOT_DIR = 'EnthalpiesCorrelationPlots'
    MOLECULE_CORRELATION_PLOT_PATH = 'molecules_error.pdf'

    def __init__(self, submissions, experimental_data, output_directory_path):
        # Build full free energy table.
        data = []

        def assign_method_class(name):
            # Roughly groups similar methods.
            method_classes = {
                'AMOEBA/BAR/Tinker': 'Alchemical/AMOEBA',
                'FS-DAM/GAFF2/TIP3P': 'Alchemical/Nonequilibrium'
            }
            if name in method_classes:
                return method_classes[name]
            if name.startswith('BSSE') or name.startswith('DFT') or name.startswith('SQM'):
                return 'QM'
            if 'MMPBSA' in name:
                return 'MMPBSA'
            if 'Umbrella Sampling' in name or (name.startswith('US') and name[-1] == '1'):
                return 'Umbrella Sampling'
            if name.startswith('US'):
                n_submission = int(name.split('_')[1])
                if n_submission == 1:
                    return 'Umbrella Sampling'
                elif n_submission == 2:
                    return 'Umbrella Sampling/Fitted'
                elif n_submission in {3, 5, 8, 10, 13, 15, 18, 20}:
                    return 'Movable Type'
                else:
                    return 'Movable Type/Fitted'
            if 'Force-Matching' in name:
                return 'Alchemical/Force-Matching'
            if name.startswith('SOMD'):
                if name[-1] == 'D':
                    return 'Alchemical/Fitted'
                else:
                    return 'Alchemical'
            if name.startswith('FEP'):
                return 'Relative/Alchemical'
            return name


        # Submissions free energies and enthalpies.
        for submission in submissions:
            for system_id, series in submission.data[['$\Delta$G', '$\Delta$H']].iterrows():
                free_energy_expt = experimental_data.loc[system_id, '$\Delta$G']
                enthalpy_expt = experimental_data.loc[system_id, '$\Delta$H']
                free_energy_calc = series['$\Delta$G']
                enthalpy_calc = series['$\Delta$H']
                data.append({
                    'receipt_id': submission.receipt_id,
                    'participant': submission.participant,
                    'name': submission.name,
                    'method': assign_method_class(submission.name),
                    'system_id': system_id,
                    'host_name': submission.data.loc[system_id, 'host_name'],
                    '$\Delta$G (calc) [kcal/mol]': free_energy_calc,
                    '$\Delta$G (expt) [kcal/mol]': free_energy_expt,
                    '$\Delta\Delta$G error (calc - expt)  [kcal/mol]': free_energy_calc - free_energy_expt,
                    '$\Delta$H (calc) [kcal/mol]': enthalpy_calc,
                    '$\Delta$H (expt) [kcal/mol]': enthalpy_expt,
                    '$\Delta\Delta$H error (calc - expt)  [kcal/mol]': enthalpy_calc - enthalpy_expt
                })

        # Transform into Pandas DataFrame.
        self.data = pd.DataFrame(data=data)
        self.output_directory_path = output_directory_path

        # Create general output directory.
        os.makedirs(self.output_directory_path, exist_ok=True)

    def generate_correlation_plots(self):
        # Free energy correlation plots.
        self._generate_correlation_plots(x='$\Delta$G (expt) [kcal/mol]', y='$\Delta$G (calc) [kcal/mol]',
                                         directory_path=self.FREE_ENERGY_CORRELATION_PLOT_DIR)
        # Enthalpies correlation plots.
        self._generate_correlation_plots(x='$\Delta$H (expt) [kcal/mol]', y='$\Delta$H (calc) [kcal/mol]',
                                         directory_path=self.ENTHALPIES_CORRELATION_PLOT_DIR)

    def _generate_correlation_plots(self, x, y, directory_path):
        output_dir_path = os.path.join(self.output_directory_path, directory_path)
        os.makedirs(output_dir_path, exist_ok=True)
        for receipt_id in self.data.receipt_id.unique():
            data = self.data[self.data.receipt_id == receipt_id]

            # If this is a merged submission, we need a hue.
            host_names = data.host_name.unique()
            if len(host_names) > 1:
                hue = 'host_name'
                title = '{} ({})'.format(data.name.unique()[0], receipt_id)
            else:
                hue = None
                title = '{} - {} ({})'.format(data.name.unique()[0], host_names[0], receipt_id)

            # Check if enthalpies were computed.
            if data[y].isnull().any():
                continue

            plt.close('all')
            plot_correlation(x=x, y=y, data=data, title=title, hue=hue)
            plt.tight_layout()
            # plt.show()
            output_path = os.path.join(output_dir_path, '{}.pdf'.format(receipt_id))
            plt.savefig(output_path)


    def generate_molecules_plot(self):
        # Correlation plot by molecules.
        plt.close('all')
        n_rows = len(self.data.system_id.unique())
        fig, ax = plt.subplots(figsize=(6, 0.4*n_rows))
        sns.violinplot(y='system_id', x='$\Delta\Delta$G error (calc - expt)  [kcal/mol]',
                       data=self.data, linewidth=1.0, inner='point', ax=ax)
        plt.tight_layout()
        # plt.show()
        plt.savefig(os.path.join(self.output_directory_path, self.MOLECULE_CORRELATION_PLOT_PATH))

    def generate_statistics_tables(self, stats_funcs, subdirectory_path, groupby,
                                   extra_fields=None, hue=None, sort_stat=None,
                                   ordering_functions=None, latex_header_conversions=None,
                                   stats_limits=None, caption=''):
        if stats_limits is None:
            stats_limits = {}
        if extra_fields is None:
            extra_fields = []

        def escape(s):
            return s.replace('_', '\_')

        extra_fields_latex = [escape(extra_field) for extra_field in extra_fields]

        file_base_name = 'statistics'
        directory_path = os.path.join(self.output_directory_path, subdirectory_path)

        stats_names, stats_funcs = zip(*stats_funcs.items())
        ci_suffixes = ('', '_lower_bound', '_upper_bound')

        # Collect the records for the DataFrames.
        statistics_csv = []
        statistics_latex = []
        statistics_plot = []

        groups = self.data[groupby].unique()
        for i, group in enumerate(groups):
            print('\rGenerating bootstrap statistics for {} {} ({}/{})'
                  ''.format(groupby, group, i+1, len(groups)), end='')

            # Select the group.
            data = self.data[self.data[groupby] == group]

            # Isolate the extra field.
            group_fields = {}
            latex_group_fields = {}
            for extra_field, extra_field_latex in zip(extra_fields, extra_fields_latex):
                assert len(data[extra_field].unique()) == 1
                extra_field_value = data[extra_field].values[0]
                group_fields[extra_field] = extra_field_value
                latex_group_fields[extra_field_latex] = escape(extra_field_value)

            hue_fields = {}
            if hue is not None:
                assert len(data[hue].unique()) == 1
                hue_fields = {hue: data[hue].values[0]}

            # Compute bootstrap statistics.
            data = data[['$\Delta$G (expt) [kcal/mol]', '$\Delta$G (calc) [kcal/mol]']]
            bootstrap_statistics = compute_bootstrap_statistics(data.as_matrix(), stats_funcs)
            # Convert to dict preserving the order of statistics.
            bootstrap_statistics = collections.OrderedDict((stats_names[i], bootstrap_statistics[i])
                                                           for i in range(len(stats_funcs)))

            record_csv = {}
            record_latex = {}
            for stats_name, (stats, (lower_bound, upper_bound), bootstrap_samples) in bootstrap_statistics.items():
                # For CSV and JSON we put confidence interval in separate columns.
                for suffix, info in zip(ci_suffixes, [stats, lower_bound, upper_bound]):
                    record_csv[stats_name + suffix] = info

                # For the PDF, print bootstrap CI in the same column.
                stats_name_latex = latex_header_conversions.get(stats_name, stats_name)
                record_latex[stats_name_latex] = '{:.2f} [{:.2f}, {:.2f}]'.format(stats, lower_bound, upper_bound)

                # For the violin plot, we need all the bootstrap statistics series.
                for bootstrap_sample in bootstrap_samples:
                    statistics_plot.append(dict(ID=group, statistics=stats_name_latex,
                                                value=bootstrap_sample,
                                                **group_fields, **hue_fields))

            statistics_csv.append({'ID': group, **group_fields, **record_csv})
            statistics_latex.append({'ID': escape(group), **latex_group_fields,
                                     **record_latex})
        print()

        # Convert dictionary to Dataframe to create tables/plots easily.
        statistics_csv = pd.DataFrame(statistics_csv)
        statistics_csv.set_index('ID', inplace=True)
        statistics_latex = pd.DataFrame(statistics_latex)
        statistics_plot = pd.DataFrame(statistics_plot)

        # Sort by the given statistics.
        if sort_stat is not None:
            ordering_function = ordering_functions.get(sort_stat, lambda x: x)
            order = sorted(statistics_csv[sort_stat].items(), key=lambda x: ordering_function(x[1]))
            order = [k for k, value in order]
            statistics_csv = statistics_csv.reindex(order)
            latex_order = [escape(k) for k in order]
            statistics_latex.ID = statistics_latex.ID.astype('category')
            statistics_latex.ID.cat.set_categories(latex_order, inplace=True)
            statistics_latex.sort_values(by='ID', inplace=True)

        # Reorder columns that were scrambled by going through a dictionaries.
        stats_names_csv = [name + suffix for name in stats_names for suffix in ci_suffixes]
        stats_names_latex = [latex_header_conversions.get(name, name) for name in stats_names]
        statistics_csv = statistics_csv[extra_fields + stats_names_csv]
        statistics_latex = statistics_latex[['ID'] + extra_fields_latex + stats_names_latex]

        # Create CSV and JSON tables (correct LaTex syntax in column names).
        os.makedirs(directory_path, exist_ok=True)
        file_base_path = os.path.join(directory_path, file_base_name)
        with open(file_base_path + '.csv', 'w') as f:
            statistics_csv.to_csv(f)
        with open(file_base_path + '.json', 'w') as f:
            statistics_csv.to_json(f, orient='index')

        # Create LaTex table.
        latex_directory_path = os.path.join(directory_path, file_base_name + 'LaTex')
        os.makedirs(latex_directory_path, exist_ok=True)
        with open(os.path.join(latex_directory_path, file_base_name + '.tex'), 'w') as f:
            f.write('\\documentclass[8pt]{article}\n'
                    '\\usepackage[a4paper,margin=0.2in,tmargin=0.5in,bmargin=0.5in,landscape]{geometry}\n'
                    '\\usepackage{booktabs}\n'
                    '\\usepackage{longtable}\n'
                    '\\pagenumbering{gobble}\n'
                    '\\begin{document}\n'
                    '\\begin{center}\n'
                    '\\begin{footnotesize}\n')
            statistics_latex.to_latex(f, column_format='|' + 'c'*(2 + len(stats_funcs)) + '|',
                                      escape=False, index=False, longtable=True, bold_rows=True)
            f.write('\end{footnotesize}\n'
                    '\end{center}\n')
            f.write(caption + '\n')
            f.write('\end{document}\n')

        # Violin plots by statistics across submissions.
        plt.close('all')
        n_cols, n_rows = len(stats_names), len(groups)

        fig, axes = plt.subplots(ncols=n_cols, figsize=(6*n_cols, 1.7*n_rows))
        for ax, stats_name in zip(axes, stats_names):
            stats_name_latex = latex_header_conversions.get(stats_name, stats_name)
            data = statistics_plot[statistics_plot.statistics == stats_name_latex]

            # Plot ordering submission by statistics.
            ordering_function = ordering_functions.get(stats_name, lambda x: x)
            order = sorted(statistics_csv[stats_name].items(), key=lambda x: ordering_function(x[1]))
            order = [receipt_id for receipt_id, value in order]
            sns.violinplot(x='value', y='ID', data=data, linewidth=1.0,
                           ax=ax, hue=hue, order=order, scale='width')

            # Configure axes.
            if hue is not None:
                ax.legend_.remove()
            if stats_name in stats_limits:
                ax.set_xlim(stats_limits[stats_name])
            ax.set_yticklabels(ax.get_yticklabels(), rotation=45)
            ax.set_xlabel(stats_name_latex)
            ax.set_ylabel('')

            handles, labels = ax.get_legend_handles_labels()
        if hue is not None:
            fig.legend(handles, labels, loc='lower right', ncol=n_cols)
        plt.tight_layout()
        # plt.show()
        plt.savefig(file_base_path + '_bootstrap_distributions.pdf')


# =============================================================================
# MERGE OA/TEMOA SUBMISSIONS
# =============================================================================

def merge_submissions(submissions, discard_not_matched=True):
    # Find all host names.
    host_names = set([submission.host_name for submission in submissions])

    # Find submissions that have the same name.
    submissions_by_name = {}
    for submission in submissions:
        try:
            submissions_by_name[submission.name].append(submission)
        except KeyError:
            submissions_by_name[submission.name] = [submission]

    # Merge OA/TEMOA submissions that use the same method into a single submission object.
    merged_submissions = []
    for method_name, method_submissions in submissions_by_name.items():
        # Check that the submissions come from the same participant,
        # and that it's the same method applied to different systems.
        assert len(method_submissions) <= len(host_names)
        assert len(set([submission.participant for submission in method_submissions])) == 1
        assert len(set([submission.host_name for submission in method_submissions])) == len(method_submissions)

        # Discard methods that were run on only a subset of the hosts.
        if len(method_submissions) == len(host_names) or not discard_not_matched:
            if len(method_submissions) == 1:
                merged_submissions.append(method_submissions[0])
            else:
                merged_submissions.append(sum(method_submissions[1:], method_submissions[0]))

    return merged_submissions


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    # TODO:     I had to fix the index CB8-G12a -> CB8-G12 in experimental_data to make the analysis work
    # TODO:     ../Submissions/974/tb3ck-974-CB8-WGatMSU-1.txt: has an extra - in CB8-G6 enthalpy
    # TODO:     ../Submissions/974/d7xde-974-CB8-NHLBI-2.txt was ignored as it is identical to 6jsye-974-CB8-NHLBI-2.txt (from two different people!)

    sns.set_style('whitegrid')
    sns.set_context('notebook')

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
    with open('../SubmissionsDoNotUpload/SAMPL6_user_map.csv', 'r') as f:
        user_map = pd.read_csv(f)

    # Configuration: statistics to compute.
    stats_funcs = collections.OrderedDict([
        ('RMSE', rmse),
        ('MAE', mae),
        ('ME', me),
        ('R2', r2),
        ('m', slope),
        ('kendall_tau', kendall_tau)
    ])
    ordering_functions = {
        'ME': lambda x: abs(x),
        'R2': lambda x: -x,
        'm': lambda x: abs(1 - x),
        'kendall_tau': lambda x: -x
    }
    latex_header_conversions = {
        'R2': 'R$^2$',
        'RMSE': 'RMSE (kcal/mol)',
        'MAE': 'MAE (kcal/mol)',
        'ME': 'ME (kcal/mol)',
        'kendall_tau': '$\\tau$',
    }
    stats_limits = {
        'RMSE': (0, 10),
        'MAE': (0, 10),
        'ME': (-5, 5),
        'R2': (0, 1),
        'm': (-3, 3),
        'kendall_tau': (-1, 1),
    }

    # Statistics by molecule.
    stats_funcs_molecules = collections.OrderedDict([
        ('RMSE', rmse),
        ('MAE', mae),
        ('ME', me),
    ])

    # Load submissions data. We do OA and TEMOA together.
    submissions_cb = load_submissions(HostGuestSubmission, HOST_GUEST_CB_SUBMISSIONS_DIR_PATH, user_map)
    submissions_oa_temoa = load_submissions(HostGuestSubmission, HOST_GUEST_OA_SUBMISSIONS_DIR_PATH, user_map)

    # Separate OA submissions from TEMOA submissions.
    submissions_oa = [submission for submission in submissions_oa_temoa if submission.host_name == 'OA']
    submissions_temoa = [submission for submission in submissions_oa_temoa if submission.host_name == 'TEMOA']

    # Merge methods that were run on both OA and TEMOA.
    submissions_oa_temoa = merge_submissions(submissions_oa_temoa)

    # Merge all methods that were run on all hosts.
    submissions_cb_oa_temoa = merge_submissions(submissions_oa_temoa + submissions_cb)

    # Merge all methods to obtain molecule statistics.
    submissions_all = merge_submissions(submissions_oa + submissions_temoa + submissions_cb,
                                        discard_not_matched=False)

    # Create submission collections
    collection_cb = HostGuestSubmissionCollection(submissions_cb, experimental_data,
                                                  output_directory_path='../CB8')
    collection_oa = HostGuestSubmissionCollection(submissions_oa, experimental_data,
                                                  output_directory_path='../OA')
    collection_temoa = HostGuestSubmissionCollection(submissions_temoa, experimental_data,
                                                     output_directory_path='../TEMOA')
    collection_oa_temoa = HostGuestSubmissionCollection(submissions_oa_temoa, experimental_data,
                                                        output_directory_path='../OA-TEMOA')
    collection_cb_oa_temoa = HostGuestSubmissionCollection(submissions_cb_oa_temoa, experimental_data,
                                                           output_directory_path='../CB8-OA-TEMOA')
    collection_all = HostGuestSubmissionCollection(submissions_all, experimental_data,
                                                   output_directory_path='../MoleculesStatistics')

    # Generate correlation plots and statistics.
    for collection in [collection_cb, collection_oa, collection_temoa,
                       collection_oa_temoa, collection_cb_oa_temoa]:
        sns.set_context('notebook')
        collection.generate_correlation_plots()

        sns.set_context('talk')
        caption = ''
        if any('NB001' in receipt_id for receipt_id in collection.data.receipt_id.unique()):
            caption += ('* NB001 was not submitted before the deadline because of a technical issue, '
                        'and it was received after the experimental results were published.')
        collection.generate_statistics_tables(stats_funcs, 'StatisticsTables', groupby='name',
                                              extra_fields=['receipt_id'], hue='method',
                                              sort_stat='RMSE', ordering_functions=ordering_functions,
                                              latex_header_conversions=latex_header_conversions,
                                              caption=caption)  # stats_limits=stats_limits,

    # Generate molecule statistics and plots. Remove null methods.
    collection_all.data.drop(collection_all.data['name'].str.startswith('NULL'))
    collection_all.generate_molecules_plot()
    collection_all.generate_statistics_tables(stats_funcs_molecules, 'StatisticsTables', groupby='system_id',
                                              hue='host_name', sort_stat='MAE', ordering_functions=ordering_functions,
                                              latex_header_conversions=latex_header_conversions)
