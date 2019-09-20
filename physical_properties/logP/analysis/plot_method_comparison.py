#!/usr/bin/env python

# =============================================================================
# GLOBAL IMPORTS
# =============================================================================
import os
import glob
import io
import collections
import pickle

import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import scipy.stats

## For QQ-Plots and Error Slope Calc
"""import scipy.stats
import scipy.integrate
import matplotlib.patches as patches
from pylab import rcParams
from operator import itemgetter, attrgetter"""

#for method comparison plot
import itertools


# =============================================================================
# CONSTANTS
# =============================================================================

# Paths to input data.
LOGP_SUBMISSIONS_DIR_PATH = './logP_predictions'
EXPERIMENTAL_DATA_FILE_PATH = '../experimental_data/logP_experimental_values.csv'


# =============================================================================
# STATS FUNCTIONS METHOD COMPARISON PLOT
# =============================================================================


def plot_correlation_with_SEM_Method_Comparison(x_lab, y_lab, x_err_lab, y_err_lab, dataX, dataY, receipt_id_1, receipt_id_2, title=None, color=None, ax=None):
    # Extract only logP values.
    x_error = dataX.loc[:, x_err_lab]
    y_error = dataY.loc[:, y_err_lab]
    x_values = dataX.loc[:, x_lab]
    y_values = dataY.loc[:, y_lab]
    #combine data for use below
    x = dataX.filter(['Molecule ID','logP (calc)'], axis=1)
    y = dataY.filter(['Molecule ID','logP (calc)'], axis=1)
    combined = pd.merge(x, y, on='Molecule ID')
    combined = combined[['logP (calc)_x', 'logP (calc)_y']]

    # Find extreme values to make axes equal.
    min_limit = np.ceil(min(combined.min()) - 1)
    max_limit = np.floor(max(combined.max()) + 1)
    axes_limits = np.array([min_limit, max_limit])

    # Color
    current_palette = sns.color_palette()
    sns_blue = current_palette[0]

    # Plot
    plt.figure(figsize=(6, 6))
    grid = sns.regplot(x=x_values, y=y_values, data=combined, color=color, ci=None)

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

    grid.set_xlabel("{} Calculated logP".format(receipt_id_1))
    grid.set_ylabel("{} Calculated logP".format(receipt_id_2))

    plt.xlim(axes_limits)
    plt.ylim(axes_limits)




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

    # The IDs of submissions used for reference calculations
    REF_SUBMISSIONS = ['REF01', 'REF02', 'REF03', 'REF04', 'REF05', 'REF06', 'REF07', 'REF08',
                       'REF09', 'REF10', 'REF11', 'REF12', 'REF13']


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

        # Check if this is a reference submission
        self.reference_submission = False
        if self.receipt_id in self.REF_SUBMISSIONS:
            self.reference_submission = True

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



# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================


def load_submissions(directory_path, user_map):
    """Load submissions from a specified directory using a specified user map.
    Optional argument:
        ref_ids: List specifying submission IDs (alphanumeric, typically) of
        reference submissions which are to be ignored/analyzed separately.
    Returns: submissions
    """
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

    LOGP_CORRELATION_PLOT_WITH_SEM_METHOD_COMPARISON_PATH_DIR = 'logPCorrelationPlotsWithSEMMethodComparison'


    def __init__(self, submissions, output_directory_path, logP_submission_collection_file_path, ignore_refcalcs = True):


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
                if submission.reference_submission and ignore_refcalcs:
                    continue
                df_collection_of_each_submission = self.data.loc[self.data["receipt_id"] == receipt_ID ]

                # Transform into Pandas DataFrame.
                submission.data = pd.DataFrame()
                submission.data["logP mean"] = df_collection_of_each_submission["logP (calc)"]
                submission.data["logP SEM"] = df_collection_of_each_submission["logP SEM (calc)"]
                submission.data["Molecule ID"] = df_collection_of_each_submission["Molecule ID"]
                submission.data["logP model uncertainty"] = df_collection_of_each_submission["logP model uncertainty"]

                submission.data.set_index("Molecule ID", inplace=True)

            # Transform into Pandas DataFrame.
            #Transform into Pandas DataFrame.
            self.data = pd.DataFrame(data=data)
            self.output_directory_path = output_directory_path
            print("\n SubmissionCollection: \n")
            print(self.data)

            # Create general output directory.
            os.makedirs(self.output_directory_path, exist_ok=True)

            # Save collection.data dataframe in a CSV file.
            self.data.to_csv(logP_submission_collection_file_path)


        '''else: # Build collection dataframe from the beginning.
            # Build full logP collection table.
            data = []

            # Submissions for logP.
            for submission in submissions:
                if submission.reference_submission and ignore_refcalcs:
                    continue
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
                    logP_model_uncertainty =  submission.data.loc[mol_ID, "logP model uncertainty"]

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
                        '$\Delta$logP error (calc - exp)': logP_mean_pred - logP_mean_exp,
                        'logP model uncertainty': logP_model_uncertainty
                    })

            # Transform into Pandas DataFrame.
            self.data = pd.DataFrame(data=data)
            self.output_directory_path = output_directory_path'''





    def generate_correlation_plots_with_SEM_method_comparison(self):
        # logP correlation plots.
        output_dir_path = os.path.join(self.output_directory_path,
                                       self.LOGP_CORRELATION_PLOT_WITH_SEM_METHOD_COMPARISON_PATH_DIR)
        os.makedirs(output_dir_path, exist_ok=True)


        for category_name, category_df in self.data.groupby('category'):
            #Make folder for each methods category
            save_path = os.path.join(output_dir_path, category_name)
            os.makedirs(save_path, exist_ok=True)
            # Paths for the different comparisons
            Reference_Reference_Comparison_PATH = os.path.join(save_path, "Reference_Reference_Comparison")
            Participant_Reference_Comparison_PATH = os.path.join(save_path, "Participant_Reference_Comparison")
            Participant_Participant_Comparison_PATH = os.path.join(save_path, "Participant_Participant_Comparison")

            all_ids = category_df.receipt_id.unique()
            #Make unique combinations of all the submissions
            unique_combinations = list(itertools.combinations(all_ids, 2))
            for receipt_id_1, receipt_id_2 in unique_combinations:
                dataX = self.data[self.data.receipt_id == receipt_id_1]
                dataY = self.data[self.data.receipt_id == receipt_id_2]
                title = '{} ({}) VS {} ({})'.format(receipt_id_1, dataX.name.unique()[0],
                                                    receipt_id_2, dataY.name.unique()[0])

                plt.close('all')
                plot_correlation_with_SEM_Method_Comparison(x_lab='logP (calc)', y_lab='logP (calc)',
                                                            x_err_lab='logP SEM (calc)', y_err_lab='logP SEM (calc)',
                                                            dataX=dataX, dataY=dataY, title=title,
                                                            receipt_id_1=receipt_id_1, receipt_id_2=receipt_id_2)
                plt.tight_layout()

                #Place method comparison plots in seperate folders
                #Not all method categories will have comparisons to reference calculations
                count_reference = sum(ID in [receipt_id_1, receipt_id_2] for ID in SamplSubmission.REF_SUBMISSIONS)
                if count_reference == 2:
                    os.makedirs(Reference_Reference_Comparison_PATH, exist_ok=True)
                    output_path = os.path.join(Reference_Reference_Comparison_PATH, '{}-{}.pdf'.format(receipt_id_1, receipt_id_2))
                    plt.savefig(output_path)
                elif count_reference == 1:
                    os.makedirs(Participant_Reference_Comparison_PATH, exist_ok=True)
                    output_path = os.path.join(Participant_Reference_Comparison_PATH, '{}-{}.pdf'.format(receipt_id_1, receipt_id_2))
                    plt.savefig(output_path)
                else:
                    os.makedirs(Participant_Participant_Comparison_PATH, exist_ok=True)
                    output_path = os.path.join(Participant_Participant_Comparison_PATH, '{}-{}.pdf'.format(receipt_id_1, receipt_id_2))
                    plt.savefig(output_path)



'''
def generate_statistics_tables(submissions, stats_funcs, directory_path, file_base_name,
                                sort_stat=None, ordering_functions=None,
                                latex_header_conversions=None, ignore_refcalcs = True):
    stats_names = list(stats_funcs.keys())
    ci_suffixes = ('', '_lower_bound', '_upper_bound')

    # Collect the records for the DataFrames.
    statistics_csv = []
    statistics_latex = []
    statistics_plot = []

    # Collect the records for QQ Plot
    # Dictionary of receipt ID: [X, Y, error_slope]
    QQplot_dict = {}

    for i, submission in enumerate(submissions):
        receipt_id = submission.receipt_id
        category = submission.category
        # Pull submission type
        type = 'Standard'
        if submission.reference_submission:
            type = 'Reference'

        # Ignore reference calculation, if applicable
        if submission.reference_submission and ignore_refcalcs:
            continue

        print('\rGenerating bootstrap statistics for submission {} ({}/{})'
                  ''.format(receipt_id, i + 1, len(submissions)), end='')

        bootstrap_statistics = submission.compute_logP_statistics(experimental_data, stats_funcs)

        # Compute error slope
        error_slope_bootstrap_statistics, QQplot_data = submission.compute_logP_model_uncertainty_statistics(experimental_data)
        #print("error_slope_bootstrap_statistics:\n")
        #print(error_slope_bootstrap_statistics)

        # Add data to to QQplot dictionary
        QQplot_dict.update({receipt_id : QQplot_data})

        # Add error slope and CI to bootstrap_statistics
        bootstrap_statistics.update({'ES' : error_slope_bootstrap_statistics })
        #print("bootstrap_statistics:\n", bootstrap_statistics)

        # Organize data to construct CSV and PDF versions of statistics tables
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
                statistics_plot.append(dict(ID=receipt_id, name=submission.name, category=category,
                                            statistics=stats_name_latex, value=bootstrap_sample))

        statistics_csv.append({'ID': receipt_id, 'name': submission.name, 'category': category, 'type': type, **record_csv})
        escaped_name = submission.name.replace('_', '\_')
        statistics_latex.append({'ID': receipt_id, 'name': escaped_name, 'category': category, 'type':type, **record_latex})
    print()
    print("statistics_csv:\n",statistics_csv)
    print()

    # Write QQplot_dict to a JSON file for plotting later
    #print("QQplot_dict:\n", QQplot_dict)
    QQplot_directory_path = os.path.join(output_directory_path, "QQPlots")
    os.makedirs(QQplot_directory_path, exist_ok=True)
    QQplot_dict_filename = os.path.join(QQplot_directory_path, 'QQplot_dict.pickle')

    with open(QQplot_dict_filename, 'wb') as outfile:
        pickle.dump(QQplot_dict, outfile)


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
    #print("stats_names_csv:", stats_names_csv)
    stats_names_latex = [latex_header_conversions.get(name, name) for name in stats_names]
    #print("stats_names_latex:", stats_names_latex)
    statistics_csv = statistics_csv[['name', "category", "type"] + stats_names_csv + ["ES", "ES_lower_bound", "ES_upper_bound"] ]
    statistics_latex = statistics_latex[['ID', 'name'] + stats_names_latex + ["ES"]] ## Add error slope(ES)

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
                '\\usepackage[a4paper,margin=0.005in,tmargin=0.5in,lmargin=0.5in,rmargin=0.5in,landscape]{geometry}\n'
                '\\usepackage{booktabs}\n'
                '\\usepackage{longtable}\n'
                '\\pagenumbering{gobble}\n'
                '\\begin{document}\n'
                '\\begin{center}\n'
                '\\scriptsize\n')
        statistics_latex.to_latex(f, column_format='|ccccccccc|', escape=False, index=False, longtable=True)
        f.write('\end{center}\n'
                '\nNotes\n\n'
                '- RMSE: Root mean square error\n\n'
                '- MAE: Mean absolute error\n\n'
                '- ME: Mean error\n\n'
                '- R2: R-squared, square of Pearson correlation coefficient\n\n'
                '- m: slope of the line fit to predicted vs experimental logP values\n\n'
                '- $\\tau$:  Kendall rank correlation coefficient\n\n'
                '- ES: error slope calculated from the QQ Plots of model uncertainty predictions\n\n'
                '- Mean and 95\% confidence intervals of RMSE, MAE, ME, R2, and m were calculated by bootstrapping with 10000 samples.\n\n'
                '- 95\% confidence intervals of ES were calculated by bootstrapping with 1000 samples.'
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
'''




# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':

    sns.set_style('whitegrid')
    sns.set_context('paper')

    '''# Read experimental data.
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
    print("Experimental data: \n", experimental_data)'''

    # Import user map.
    with open('../predictions/SAMPL6-user-map-logP.csv', 'r') as f:
        user_map = pd.read_csv(f)

    # Configuration: statistics to compute.
    '''stats_funcs = collections.OrderedDict([
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
        'RMSE': 'RMSE',
        'MAE': 'MAE',
        'ME': 'ME',
        'kendall_tau': '$\\tau$'
    }'''


    #==========================================================================================
    # Run analysis including reference calculations
    #==========================================================================================
    # Load submissions data.
    submissions_logP = load_submissions(LOGP_SUBMISSIONS_DIR_PATH, user_map)
    # Perform the analysis

    output_directory_path='./analysis_outputs_withrefs'
    logP_submission_collection_file_path = '{}/logP_submission_collection.csv'.format(output_directory_path)

    collection_logP= logPSubmissionCollection(submissions_logP, output_directory_path, logP_submission_collection_file_path, ignore_refcalcs = False)

    # Generate plots and tables.
    for collection in [collection_logP]:
        collection.generate_correlation_plots_with_SEM_method_comparison()
