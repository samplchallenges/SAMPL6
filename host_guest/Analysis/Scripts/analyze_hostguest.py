#!/usr/bin/env python

# =============================================================================
# GLOBAL IMPORTS
# =============================================================================

import os
import copy
import collections
import pickle

import numpy as np
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

    _ROW_HEIGHT = 0.6

    def __init__(self, submissions, experimental_data, output_directory_path):
        # Build full free energy table.
        data = []

        # Submissions free energies and enthalpies.
        for submission in submissions:
            for system_id, series in submission.data[['$\Delta$G', 'SEM $\Delta$G', '$\Delta$H']].iterrows():
                free_energy_expt = experimental_data.loc[system_id, '$\Delta$G']
                enthalpy_expt = experimental_data.loc[system_id, '$\Delta$H']
                free_energy_calc = series['$\Delta$G']
                free_energy_calc_sem = series['SEM $\Delta$G']
                enthalpy_calc = series['$\Delta$H']
                data.append({
                    'receipt_id': submission.receipt_id,
                    'participant': submission.participant,
                    'name': submission.name,
                    'method': self._assign_paper_method_name(submission.name),
                    'system_id': system_id,
                    'host_name': submission.data.loc[system_id, 'host_name'],
                    '$\Delta$G (calc) [kcal/mol]': free_energy_calc,
                    'SEM $\Delta$G (calc) [kcal/mol]': free_energy_calc_sem,
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

    @staticmethod
    def _assign_method_class(name):
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
                return 'Alchemical/GAFF/Fitted'
            else:
                return 'Alchemical/GAFF'
        if 'GAFF' in name:
            return 'Alchemical/GAFF'
        if name.startswith('FEP'):
            return 'Alchemical/Relative'
        return name

    @staticmethod
    def _assign_paper_method_name(name):
        # Convert from submission method name to the name used in the paper.
        method_names = {
            'DDM/GAFF/AM1-BCC/TIP3P': 'FEP-GAFF',
            'HREM/BAR/RESP/Force-Matching/TIP3P': 'ForceMatch',
            'DDM/Force-Matching/FEP/HREM/MBAR': 'ForceMatch-QMMM',
            'BSSE-corrected RI-B3PW91 (SMD)/CBS': 'DFT(B3PW91)',
            'BSSE-corrected RI-B3PW91-D3 (SMD)/CBS': 'DFT(B3PW91)-D3',
            'DFT-opt': 'DFT(TPSS)-D3',
            'FEP-MM': 'RelFEP-GAFF2',
            'FEP-QM/MM': 'RelFEP-GAFF2-QMMM',
            'FS-DAM/GAFF2/TIP3P': 'FastSwitch',
            'EKEN-DIAZ/MD/MMPBSA': 'MMPBSA',
            'SQM-opt': 'SQM(PM6-DH+)',
            'AMOEBA/BAR/Tinker': 'Tinker-BAR',
            'Umbrella Sampling/TIP3P': 'US-CGenFF',
            'US/PMF/MT/MD_1': 'US-GAFF',
            'US/PMF/MT/MD_2': 'US-GAFF-C',
        }
        try:
            return method_names[name]
        except KeyError:
            pass

        if name.startswith('SOMD/AM1BCC-GAFF-TIP3P'):
            paper_name = 'SOMD-' + name[-1]
            if 'NOBUFFER' in name:
                paper_name += '-nobuffer'
        if name.startswith('US/PMF/MT/MD'):
            submission_number = name.rsplit('_', 1)[-1]
            paper_name = 'MovTyp-' + submission_number
        if 'NULL' in name:
            paper_name = name
        return paper_name

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
                                   extra_fields=None, sort_stat=None,
                                   ordering_functions=None, latex_header_conversions=None,
                                   caption=''):
        """Generate statistics tables in CSV, JSON, and LaTex format.

        Parameters
        ----------
        groupby : str
            The name of the data column to be used to compute the statistics.
            For example, 'name' to obtain statistics about individual methods,
            'system_id' to compute statistics by molecules.
        ordering_functions : dict
            Dictionary statistic_name -> ordering_function(stats), where
            ordering_function determines how to rank the the groups by
            statistics.
        """
        if extra_fields is None:
            extra_fields = []

        def escape(s):
            return s.replace('_', '\_')

        extra_fields_latex = [escape(extra_field) for extra_field in extra_fields]

        file_base_name = 'statistics'
        directory_path = os.path.join(self.output_directory_path, subdirectory_path)

        stats_names, stats_funcs = zip(*stats_funcs.items())
        ci_suffixes = ('', '_lower_bound', '_upper_bound')

        # Compute or read the bootstrap statistics from the cache.
        cache_file_path = os.path.join(self.output_directory_path, 'bootstrap_distributions.p')
        all_bootstrap_statistics = self._get_bootstrap_statistics(groupby, stats_names, stats_funcs,
                                                                  cache_file_path=cache_file_path)

        # Collect the records for the DataFrames.
        statistics_csv = []
        statistics_latex = []

        groups = self.data[groupby].unique()
        for i, group in enumerate(groups):
            print('\rGenerating bootstrap statistics tables for {} {} ({}/{})'
                  ''.format(groupby, group, i+1, len(groups)), end='')

            # Isolate bootstrap statistics.
            bootstrap_statistics = all_bootstrap_statistics[group]

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

            record_csv = {}
            record_latex = {}
            for stats_name, (stats, (lower_bound, upper_bound), bootstrap_samples) in bootstrap_statistics.items():
                # For CSV and JSON we put confidence interval in separate columns.
                for suffix, info in zip(ci_suffixes, [stats, lower_bound, upper_bound]):
                    record_csv[stats_name + suffix] = info

                # For the PDF, print bootstrap CI in the same column.
                stats_name_latex = latex_header_conversions.get(stats_name, stats_name)
                record_latex[stats_name_latex] = '{:.2f} [{:.2f}, {:.2f}]'.format(stats, lower_bound, upper_bound)

            statistics_csv.append({'ID': group, **group_fields, **record_csv})
            statistics_latex.append({'ID': escape(group), **latex_group_fields,
                                     **record_latex})
        print()

        # Convert dictionary to Dataframe to create tables/plots easily.
        statistics_csv = pd.DataFrame(statistics_csv)
        statistics_csv.set_index('ID', inplace=True)
        statistics_latex = pd.DataFrame(statistics_latex)

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

    def plot_bootstrap_distributions(self, stats_funcs, subdirectory_path, groupby,
                                     ordering_functions=None, latex_header_conversions=None,
                                     stats_limits=None, **violinplot_kwargs):
        """Generate a violin plot for the bootstrap distribution of the statistics.

        Parameters
        ----------
        subdirectory_path : str
            The path of the output directory relative to self.output_directory_path.
        groupby : str
            The name of the data column to be used to compute the statistics.
            For example, 'name' to obtain statistics about individual methods,
            'system_id' to compute statistics by molecules.
        ordering_functions : dict
            Dictionary statistic_name -> ordering_function(stats), where
            ordering_function determines how to rank the the groups by
            statistics.

        """
        if stats_limits is None:
            stats_limits = {}

        directory_path = os.path.join(self.output_directory_path, subdirectory_path)
        stats_names, stats_funcs = zip(*stats_funcs.items())

        # Compute or read the bootstrap statistics from the cache.
        ordering_data, statistics_plot = self._get_bootstrap_distribution_plot_data(
            groupby, stats_names, stats_funcs)

        # Create CSV and JSON tables (correct LaTex syntax in column names).
        os.makedirs(directory_path, exist_ok=True)

        # Violin plots by statistics across submissions.
        n_groups = len(ordering_data.index)
        for stats_name in stats_names:
            plt.close('all')
            fig, ax = plt.subplots(figsize=(8, self._ROW_HEIGHT*n_groups))

            stats_name_latex = latex_header_conversions.get(stats_name, stats_name)
            data = statistics_plot[statistics_plot.statistics_name == stats_name_latex]

            # Determine the order in which to display groups.
            ordering_function = ordering_functions.get(stats_name, lambda x: x)
            order = sorted(ordering_data[stats_name].items(), key=lambda x: ordering_function(x[1]))
            order = [group for group, stats in order]

            # Plot boot strap distributions.
            sns.violinplot(x='value', y='ID', data=data, linewidth=1.0,
                           order=order, scale='width', ax=ax, **violinplot_kwargs)

            self._modify_violinplot(ax, stats_name)

            # Configure axes.
            if stats_name in stats_limits:
                ax.set_xlim(stats_limits[stats_name])
            ax.set_xlabel(stats_name_latex)
            ax.set_ylabel('')

            # if hue is not None:
            #     ax.legend_.remove()
            # handles, labels = ax.get_legend_handles_labels()
            # if hue is not None:
            #     fig.legend(handles, labels, loc='lower right')
            plt.tight_layout()
            # plt.show()
            plt.savefig(os.path.join(directory_path, stats_name) + '_bootstrap_distributions.pdf')
            # plt.savefig(os.path.join(directory_path, stats_name) + '_bootstrap_distributions.png', dpi=300)

    def _get_bootstrap_distribution_plot_data(self, groupby, stats_names, stats_funcs):
        """Return the dataframes with the statistics necessary to plot.

        Returns
        -------
        ordering_data : pandas.Dataframe
            A dataframe containing all the statistics that can be used to decide
            the order of the bootstrap distributions to plot.
        statistics_plot_data : pandas.Dataframe
            A dataframe containing all the bootstrap samples.
        """
        cache_file_path = os.path.join(self.output_directory_path, 'bootstrap_distributions.p')
        all_bootstrap_statistics = self._get_bootstrap_statistics(groupby, stats_names, stats_funcs,
                                                                  cache_file_path=cache_file_path)

        # Collect the records for the DataFrames.
        ordering_data = []
        statistics_plot = []

        groups = self.data[groupby].unique()
        for i, group in enumerate(groups):
            print('\rCollecting bootstrap statistics for {} {} ({}/{})'
                  ''.format(groupby, group, i+1, len(groups)), end='')

            # Isolate bootstrap statistics.
            bootstrap_statistics = all_bootstrap_statistics[group]

            ordering_data_record = {}
            for stats_name, (stats, (lower_bound, upper_bound), bootstrap_samples) in bootstrap_statistics.items():
                stats_name_latex = latex_header_conversions.get(stats_name, stats_name)
                # Keep the data used for the order.
                # TODO check primary_hue
                ordering_data_record[stats_name] = stats
                # For the violin plot, we need all the bootstrap statistics series.
                for bootstrap_sample in bootstrap_samples:
                    statistics_plot.append(dict(ID=group, statistics_name=stats_name_latex,
                                                value=bootstrap_sample))
            ordering_data.append(dict(ID=group, **ordering_data_record))
        print()

        ordering_data = pd.DataFrame(ordering_data)
        ordering_data.set_index('ID', inplace=True)
        statistics_plot = pd.DataFrame(statistics_plot)
        return ordering_data, statistics_plot

    def _get_bootstrap_statistics(self, groupby, stats_names, stats_funcs, cache_file_path):
        """Generate the bootstrap distributions of all groups and cache them.

        If cached values are found on disk, the distributions are not recomputed.

        Returns
        -------
        all_bootstrap_statistics : collections.OrderedDict
            group -> {stats_name -> (statistics, confidence_interval, bootstrap_samples)}
            confidence_interval is a pair (lower_bound, upper_bound), and bootstrap_samples
            are the (ordered) bootstrap statistics used to compute the confidence interval.
        """
        # Identify all the groups (e.g. methods/molecules).
        groups = self.data[groupby].unique()

        # Initialize returned value. The OrderedDict maintains the order of statistics.
        all_bootstrap_statistics = collections.OrderedDict([(name, None) for name in stats_names])
        all_bootstrap_statistics = collections.OrderedDict(
            [(group, copy.deepcopy(all_bootstrap_statistics)) for group in groups]
        )

        # Load the statistics that we have already computed.
        try:
            with open(cache_file_path, 'rb') as f:
                print('Loading cached bootstrap distributions from {}'.format(cache_file_path))
                cached_bootstrap_statistics = pickle.load(f)
        except FileNotFoundError:
            cached_bootstrap_statistics = None

        cache_updated = False
        for i, (group, group_bootstrap_statistics) in enumerate(all_bootstrap_statistics.items()):
            # Check which statistics we still need to compute for this group.
            if cached_bootstrap_statistics is not None:
                group_stats_names = []
                group_stats_funcs = []
                for stats_name, stats_func in zip(stats_names, stats_funcs):
                    try:
                        all_bootstrap_statistics[group][stats_name] = cached_bootstrap_statistics[group][stats_name]
                    except KeyError:
                        try:
                            method_name = self._assign_paper_method_name(group)
                            all_bootstrap_statistics[group][stats_name] = cached_bootstrap_statistics[method_name][stats_name]
                        except KeyError:
                            group_stats_names.append(stats_name)
                            group_stats_funcs.append(stats_func)
            else:
                # Compute everything.
                group_stats_names = stats_names
                group_stats_funcs = stats_funcs

            if len(group_stats_names) == 0:
                continue
            cache_updated = True  # Update the cache on disk later.

            print('\rGenerating bootstrap statistics for {} {} ({}/{})'
                  ''.format(groupby, group, i+1, len(groups)), end='')

            # Select the group data.
            data = self.data[self.data[groupby] == group]

            # Check if SEMs for the free energies are reported.
            sems = data['SEM $\Delta$G (calc) [kcal/mol]'].values
            if np.any(np.isnan(sems)):
                sems = None
            else:  # Add a column of SEMs = 0.0 for the experimental values.
                sems = np.array([(0.0, sem) for sem in sems])

            # Compute bootstrap statistics.
            data = data[['$\Delta$G (expt) [kcal/mol]', '$\Delta$G (calc) [kcal/mol]']]
            new_bootstrap_statistics = compute_bootstrap_statistics(data.as_matrix(), group_stats_funcs, sems=sems,
                                                                    n_bootstrap_samples=100000)

            # Update the returned value with the statistics just computed.
            new_boostrap_statistics = {group_stats_names[i]: new_bootstrap_statistics[i]
                                       for i in range(len(group_stats_funcs))}
            group_bootstrap_statistics.update(new_boostrap_statistics)

        # Cache the computed statistics on disk. Create output directory if necessary.
        if cache_updated:
            os.makedirs(os.path.dirname(cache_file_path), exist_ok=True)
            with open(cache_file_path, 'wb') as f:
                pickle.dump(all_bootstrap_statistics, f)

        return all_bootstrap_statistics

    def _modify_violinplot(self, ax, stats_name):
        pass

# =============================================================================
# MERGE SUBMISSIONS AND COLLECTIONS
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


class SplitBootstrapSubmissionCollection(HostGuestSubmissionCollection):
    """Two collections merged, which allows to plot split bootstrap distributions.

    The violin plots are ordered according to collection1.

    Parameters
    ----------
    collection1 : HostGuestSubmissionCollection
        The first collection.
    collection2 : HostGuestSubmissionCollection
        The second collection.
    hue : str
        The name of the hue parameter in seaborn violinplot.
    collection1_hue : str
        The hue value for collection1.
    collection2_hue : str
        The hue value for collection2.
    output_directory_path : str
        The path of the directory where to save the results of the analysis.
    """

    def __init__(self, collection1, collection2, hue, collection1_hue, collection2_hue, output_directory_path):
        self.hue = hue
        self.collections = collections.OrderedDict([
            (collection1_hue, collection1),
            (collection2_hue, collection2)
        ])

        # Create output directory if needed.
        self.output_directory_path = output_directory_path
        os.makedirs(self.output_directory_path, exist_ok=True)

    def _get_bootstrap_distribution_plot_data(self, groupby, stats_names, stats_funcs):
        """Return the dataframes with the statistics necessary to plot.

        Returns
        -------
        ordering_data : pandas.Dataframe
            A dataframe containing all the statistics that can be used to decide
            the order of the bootstrap distributions to plot.
        statistics_plot_data : pandas.Dataframe
            A dataframe containing all the bootstrap samples.
        """
        # Collect the records for the DataFrames.
        ordering_data = []
        statistics_plot = []
        # Cache the statistics of the collections as we'll need them in _modify_violinplot().
        self._collections_statistics = []

        # Gather data for both collection distributions.
        for collection_hue, collection in self.collections.items():
            cache_file_path = os.path.join(collection.output_directory_path, 'bootstrap_distributions.p')
            all_bootstrap_statistics = collection._get_bootstrap_statistics(groupby, stats_names, stats_funcs,
                                                                            cache_file_path=cache_file_path)

            groups = collection.data[groupby].unique()

            for group_idx, group in enumerate(groups):
                # TODO REMOVE ME
                if 'SOMD' not in group:
                    continue
                print('\rCollecting bootstrap statistics for {} {} ({}/{})'
                      ''.format(groupby, group, group_idx+1, len(groups)), end='')
                # Isolate bootstrap statistics.
                bootstrap_statistics = all_bootstrap_statistics[group]

                ordering_data_record = {}
                for stats_name, (stats, (lower_bound, upper_bound), bootstrap_samples) in bootstrap_statistics.items():
                    stats_name_latex = latex_header_conversions.get(stats_name, stats_name)
                    # Keep the data used for the order.
                    ordering_data_record[stats_name] = stats
                    # For the violin plot, we need all the bootstrap statistics series.
                    for bootstrap_sample in bootstrap_samples:
                        statistics_plot.append(dict(ID=group, statistics_name=stats_name_latex,
                                                    value=bootstrap_sample, **{self.hue: collection_hue}))

                # We order using collection1 data.
                if collection_hue == list(self.collections.keys())[0]:
                    ordering_data.append(dict(ID=group, **ordering_data_record))
                self._collections_statistics.append(dict(ID=group, **ordering_data_record,
                                                         **{self.hue: collection_hue}))

            del all_bootstrap_statistics  # This can be huge when the number of bootstrap cycles is high.
            print()

        ordering_data = pd.DataFrame(ordering_data)
        ordering_data.set_index('ID', inplace=True)
        statistics_plot = pd.DataFrame(statistics_plot)
        self._collections_statistics = pd.DataFrame(self._collections_statistics)
        return ordering_data, statistics_plot

    def _modify_violinplot(self, ax, stats_name):
        """Add the real statistics to the violin plot."""
        for group_idx, group in enumerate(ax.get_yticklabels()):
            for collection_idx, (collection_hue, collection) in enumerate(self.collections.items()):
                collection_data = self._collections_statistics[(self._collections_statistics.ID == group.get_text()) &
                                                               (self._collections_statistics[self.hue] == collection_hue)]
                try:
                    collection_stats = collection_data[stats_name].values[0]
                except IndexError:
                    # The group doesn't have a statistics for the hue.
                    continue
                sign = -1 if collection_idx == 0 else 1
                y = np.linspace(group_idx, group_idx + sign * self._ROW_HEIGHT / 1.8)
                x = [collection_stats for _ in y]
                ax.plot(x, y, c='black', lw=3, alpha=0.6)


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
        # 'RMSE': (0, 10),
        # 'MAE': (0, 10),
        # 'ME': (-5, 5),
        'R2': (0, 1),
        # 'm': (-3, 3),
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

    # Create a CB8 collection excluding the bonus challenges.
    def remove_bonus(submission_collection_data):
        return submission_collection_data[(submission_collection_data.system_id != 'CB8-G11') &
                                          (submission_collection_data.system_id != 'CB8-G12') &
                                          (submission_collection_data.system_id != 'CB8-G13')]
    collection_cb_no_bonus = HostGuestSubmissionCollection(submissions_cb, experimental_data,
                                                           output_directory_path='../CB8-NOBONUS')
    collection_cb_no_bonus.data = remove_bonus(collection_cb_no_bonus.data)

    # Generate correlation plots and statistics.
    for collection in [collection_cb, collection_cb_no_bonus, collection_oa, collection_temoa,
                       collection_oa_temoa, collection_cb_oa_temoa]:
        sns.set_context('notebook')
        collection.generate_correlation_plots()

        sns.set_context('talk')
        caption = ''
        if any('NB001' in receipt_id for receipt_id in collection.data.receipt_id.unique()):
            caption += ('* NB001 was not submitted before the deadline because of a technical issue, '
                        'and it was received after the experimental results were published.')
        collection.generate_statistics_tables(stats_funcs, subdirectory_path='StatisticsTables',
                                              groupby='name', extra_fields=['receipt_id'],
                                              sort_stat='RMSE', ordering_functions=ordering_functions,
                                              latex_header_conversions=latex_header_conversions,
                                              caption=caption)
        collection.plot_bootstrap_distributions(stats_funcs, subdirectory_path='StatisticsPlots',
                                                groupby='name', ordering_functions=ordering_functions,
                                                latex_header_conversions=latex_header_conversions,
                                                stats_limits=stats_limits)

    # Generate molecule statistics and plots. Remove null methods.
    collection_all.data.drop(collection_all.data['name'].str.startswith('NULL'))
    collection_all.generate_molecules_plot()
    collection_all.generate_statistics_tables(stats_funcs_molecules, 'StatisticsTables', groupby='system_id',
                                              sort_stat='MAE', ordering_functions=ordering_functions,
                                              latex_header_conversions=latex_header_conversions)
    collection_all.plot_bootstrap_distributions(stats_funcs_molecules, subdirectory_path='StatisticsPlots',
                                                groupby='system_id', ordering_functions=ordering_functions,
                                                latex_header_conversions=latex_header_conversions,
                                                stats_limits=stats_limits)

    # # Violin plots of bootstrap distribution for two collections.
    # collection = SplitBootstrapSubmissionCollection(collection_oa_temoa, collection_cb,
    #                                                 hue='dataset', collection1_hue='OA/TEMOA', collection2_hue='CB8',
    #                                                 output_directory_path='../MergedCB8-OA')
    # collection.plot_bootstrap_distributions(stats_funcs, 'StatisticsPlots', groupby='method',
    #                                         ordering_functions=ordering_functions, stats_limits=stats_limits,
    #                                         latex_header_conversions=latex_header_conversions,
    #                                         hue='dataset', split=True)
