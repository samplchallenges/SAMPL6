#!/usr/bin/env python

# =============================================================================
# GLOBAL IMPORTS
# =============================================================================

import os
import csv
import json
import collections

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from pkganalysis.submission import (IgnoredSubmissionError, load_submissions)
from pkganalysis.sampling import (SamplingSubmission, YankSamplingAnalysis,
                                  YANK_N_ITERATIONS, DG_KEY, export_dictionary)


# =============================================================================
# CONSTANTS
# =============================================================================

# Paths to input data.
SAMPLING_SUBMISSIONS_DIR_PATH = '../SubmissionsDoNotUpload/975/'
YANK_ANALYSIS_DIR_PATH = 'YankAnalysis/Sampling/'
SAMPLING_ANALYSIS_DIR_PATH = '../SAMPLing/'
SAMPLING_DATA_DIR_PATH = os.path.join(SAMPLING_ANALYSIS_DIR_PATH, 'Data')
SAMPLING_PLOT_DIR_PATH = os.path.join(SAMPLING_ANALYSIS_DIR_PATH, 'Plots')

# All system ids.
SYSTEM_IDS = [
    'CB8-G3-0', 'CB8-G3-1', 'CB8-G3-2', 'CB8-G3-3', 'CB8-G3-4',
    'OA-G3-0', 'OA-G3-1', 'OA-G3-2', 'OA-G3-3', 'OA-G3-4',
    'OA-G6-0', 'OA-G6-1', 'OA-G6-2', 'OA-G6-3', 'OA-G6-4'
]


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def load_yank_analysis():
    """Load the YANK analysis in a single dataframe."""
    yank_free_energies = {}
    for system_id in SYSTEM_IDS:
        file_path = os.path.join(YANK_ANALYSIS_DIR_PATH, 'yank-{}.json'.format(system_id))
        with open(file_path, 'r') as f:
            yank_free_energies[system_id] = json.load(f)
    return yank_free_energies


def export_submissions(submissions, reference_free_energies):
    """Export the submission data to CSV and JSON format."""
    for submission in submissions:
        mean_free_energies = submission.mean_free_energies()

        # Export data of mean trajectory and confidence intervals.
        exported_data = {}
        for system_name in mean_free_energies['System name'].unique():
            system_name_data = mean_free_energies[mean_free_energies['System name'] == system_name]

            # Obtain free energies and bias.
            free_energies = system_name_data[DG_KEY].values
            free_energies_ci = system_name_data['$\Delta$G CI'].values
            reference_diff = free_energies - reference_free_energies.loc[system_name, '$\Delta$G [kcal/mol]']

            exported_data[system_name + '-mean'] = collections.OrderedDict([
                ('DG', free_energies.tolist()),
                ('DG_CI', free_energies_ci.tolist()),
                ('reference_difference', reference_diff.tolist()),
                ('n_energy_evaluations', system_name_data['N energy evaluations'].values.tolist()),
            ])

        # Export.
        file_base_path = os.path.join(SAMPLING_DATA_DIR_PATH, submission.receipt_id)
        export_dictionary(exported_data, file_base_path)


def plot_mean_free_energy(mean_data, ax, x='Simulation percentage', **kwargs):
    """Plot mean trajectory with confidence intervals."""
    ci_key = '$\Delta$G CI'

    # Plot the mean free energy trajectory.
    ax.plot(mean_data[x].values, mean_data[DG_KEY].values, **kwargs)

    # Plot mean trajectory confidence intervals.
    mean_dg = mean_data[DG_KEY].values
    sem_dg = mean_data[ci_key].values
    ax.fill_between(mean_data[x].values, mean_dg + sem_dg, mean_dg - sem_dg, alpha=0.65)
    return ax


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    sns.set_style('whitegrid')
    sns.set_context('talk')

    # Read reference values.
    yank_analysis = YankSamplingAnalysis(YANK_ANALYSIS_DIR_PATH)

    # Obtain free energies and final reference values.
    reference_free_energies = yank_analysis.free_energies_from_iteration(YANK_N_ITERATIONS, mean_trajectory=True)
    reference_free_energies = reference_free_energies[reference_free_energies['Simulation percentage'] == 100]
    reference_free_energies.set_index('System name', inplace=True)

    # Import user map.
    with open('../SubmissionsDoNotUpload/SAMPL6_user_map.csv', 'r') as f:
        user_map = pd.read_csv(f)

    # Load submissions data. We do OA and TEMOA together.
    submissions = load_submissions(SamplingSubmission, SAMPLING_SUBMISSIONS_DIR_PATH, user_map)

    # Export YANK analysis and submissions to CSV/JSON tables.
    yank_analysis.export(os.path.join(SAMPLING_DATA_DIR_PATH, 'reference_free_energies'))
    export_submissions(submissions, reference_free_energies)

    # Plot submission data.
    os.makedirs(SAMPLING_PLOT_DIR_PATH, exist_ok=True)
    palette_mean = sns.color_palette('dark')
    for submission in submissions:
        mean_free_energies = submission.mean_free_energies()
        for system_name in submission.data['System name'].unique():
            # Select the data for only this host-guest system.
            data = submission.data[submission.data['System name'] == system_name]
            mean_data = mean_free_energies[mean_free_energies['System name'] == system_name]

            # Get the corresponding YANK free energies.
            # The cost for the same system is the same for all replicates.
            n_energy_evaluations = int(submission.cost.loc[system_name + '-0', 'N energy evaluations'])
            yank_mean_data = yank_analysis.free_energies_from_energy_evaluations(n_energy_evaluations,
                                                                                 system_name=system_name,
                                                                                 mean_trajectory=True)

            fig, ax = plt.subplots()

            # Plot the 5 replicates trajectories.
            sns.tsplot(data=data, time='Simulation percentage', value=DG_KEY,
                       unit='System name', condition='System ID', color='pastel', ax=ax)

            # Plot the submission mean trajetory with CI.
            plot_mean_free_energy(mean_data, ax=ax, color=palette_mean[0],
                                  label='Mean $\Delta$G')

            # Plot YANK mean trajectory with CI.
            plot_mean_free_energy(yank_mean_data, ax=ax, color=palette_mean[2],
                                  label='Ref mean $\Delta$G')

            # Plot reference value.
            ref_free_energy = reference_free_energies.loc[system_name, DG_KEY]
            ax.plot(mean_data['Simulation percentage'], [ref_free_energy for _ in range(100)],
                    color=palette_mean[1], ls='--', label='Ref final $\Delta$G')

            # Set axis limits/titles.
            ax.set_ylim((-20, 4))
            ax.set_title('{} - {} ({})'.format(system_name, submission.name, submission.receipt_id))
            ax.set_xlabel(ax.get_xlabel() + ' (N energy evaluations: {:,})'.format(n_energy_evaluations))
            ax.legend(ncol=2)
            plt.tight_layout()
            output_path_dir = os.path.join(SAMPLING_PLOT_DIR_PATH,
                                           '{}-{}.pdf'.format(submission.receipt_id, system_name))
            # plt.show()
            plt.savefig(output_path_dir)

    # Plot YANK replicates isolated.
    for system_name in reference_free_energies.index:
        # Select the data for only this host-guest system.
        yank_data = yank_analysis.system_free_energies(system_name)
        yank_mean_data = yank_analysis.system_free_energies(system_name, mean_trajectory=True)

        fig, ax = plt.subplots()

        # Plot the 5 replicates trajectories.
        sns.tsplot(data=yank_data, time='N energy evaluations', value=DG_KEY,
                   unit='System name', condition='System ID', color='pastel', ax=ax)

        # Plot the submission mean trajectory with CI.
        plot_mean_free_energy(yank_mean_data, ax=ax, color=palette_mean[0],
                              x='N energy evaluations', label='Mean $\Delta$G')

        # Set axis limits/titles.
        ax.set_title('{} - Reference'.format(system_name))
        ax.legend(ncol=2)
        plt.tight_layout()
        output_path_dir = os.path.join(SAMPLING_PLOT_DIR_PATH,
                                       'reference-{}.pdf'.format(system_name))
        # plt.show()
        plt.savefig(output_path_dir)



