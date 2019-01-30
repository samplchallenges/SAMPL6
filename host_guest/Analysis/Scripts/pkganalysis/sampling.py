#!/usr/bin/env python

"""Utility classes and functions for SAMPLing analysis."""


# =============================================================================
# GLOBAL IMPORTS
# =============================================================================

import os
import csv
import glob
import json
import operator
import collections

import numpy as np
import pandas as pd

from .submission import SamplSubmission
from .stats import mean_confidence_interval


# =============================================================================
# CONSTANTS
# =============================================================================

N_REPLICATES = 5

# All YANK calculations and all phases have run for the same number of iterations/steps.
YANK_N_ITERATIONS = 40000
N_STEPS_PER_ITERATIONS = 500

# Number of states in the alchemical protocol (complex + solvent) for each system.
YANK_N_STATES = {
    'CB8-G3': 69 + 62,
    'OA-G3': 59 + 54,
    'OA-G6': 55 + 52
}

DG_KEY = '$\Delta$G [kcal/mol]'
DDG_KEY = 'd$\Delta$G [kcal/mol]'


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def export_dictionary(data_dict, file_base_path):
    """Export the data in CSV and JSON format.

    Parameters
    ----------
    data_dict : dict
        Dictionary system ID -> {key: list of values}.
    file_base_path : str
        The extension-less path to the output files.

    """
    os.makedirs(os.path.dirname(file_base_path), exist_ok=True)

    # Export to JSON as it is.
    with open(file_base_path + '.json', 'w') as f:
        json.dump(data_dict, f, indent=4, sort_keys=True)

    # Flatten for CSV format.
    csv_data = []
    for system_id in sorted(data_dict):
        for key, values in data_dict[system_id].items():
            csv_data.append([system_id + '-' + key] + values)

    # Export CSV file.
    with open(file_base_path + '.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerows(csv_data)


def compute_system_name_mean_free_energies(data, reference_free_energies=None, extra_fields=()):
    """Compute mean free energies and CI for each system name in data."""
    output_data = []
    system_names = data['System name'].unique()
    for system_name in system_names:
        system_name_data = data[data['System name'] == system_name]

        # Obtain number of energy evaluations for data point.
        energy_evaluations = system_name_data['N energy evaluations'].unique().tolist()
        tot_energy_evaluations = max(energy_evaluations)

        # Obtain final free energy to be used for computing the RMSE of the estimator.
        if reference_free_energies is None:
            last_data_point = system_name_data[system_name_data['N energy evaluations'] == tot_energy_evaluations]
            reference_free_energy, _ = mean_confidence_interval(last_data_point[DG_KEY])
        else:
            reference_free_energy = reference_free_energies[system_name]

        # Add records.
        for idx, n_energy_evaluations in enumerate(energy_evaluations):
            time_point_data = system_name_data[system_name_data['N energy evaluations'] == n_energy_evaluations]
            free_energies = time_point_data[DG_KEY]
            f, ci = mean_confidence_interval(free_energies.values, confidence=0.95)
            std = np.std(free_energies.values)
            rmse = np.sqrt(1/(len(free_energies.values)) * np.sum((free_energies.values - reference_free_energy)**2))
            bias = f - reference_free_energy

            extra_values = {}
            for field in extra_fields:
                assert len(time_point_data[field].unique()) == 1
                extra_values[field] = time_point_data[field].values[0]

            # Keep "Simulation percentage" an integer if possible.
            if len(energy_evaluations) == 100:
                simulation_percentage = energy_evaluations.index(n_energy_evaluations) + 1
            else:
                simulation_percentage = n_energy_evaluations * 100 / tot_energy_evaluations

            output_data.append({
                'System name': system_name,
                DG_KEY: f,
                'std': std,
                'RMSE': rmse,
                'bias': bias,
                '$\Delta$G CI': ci,
                'Simulation percentage': simulation_percentage,
                'N energy evaluations': n_energy_evaluations,
                **extra_values
            })
    columns_order = ['System name',  DG_KEY, 'std', 'RMSE', 'bias', '$\Delta$G CI', 'Simulation percentage',
                     'N energy evaluations'] + list(extra_fields)
    output_data = pd.DataFrame(output_data, columns=columns_order)
    return output_data


# =============================================================================
# CONVERSION FROM ENERGY EVALUATIONS/CPU TIME TO ITERATIONS
# =============================================================================

def energy_evaluations_from_iterations(system_name, n_iterations):
    """Compute the number of energy evaluations necessary to run N iteration"""
    n_states = YANK_N_STATES[system_name]
    md_energy_evaluations = n_states * N_STEPS_PER_ITERATIONS
    mc_energy_evaluations = 2 * 2 * n_states  # Rotation and translation, initial and final energy.
    energy_matrix_evaluations = n_states * n_states  # Technically, we compute only the changed force groups.
    energy_evaluations_per_iteration = md_energy_evaluations + mc_energy_evaluations + energy_matrix_evaluations
    return energy_evaluations_per_iteration * n_iterations


def energy_evaluations_iteration_cutoffs(tot_energy_evaluations, system_name):
    """Compute the 100 YANK iterations to use for analysis from the total energy evaluations to consider."""
    # Find the number of energy evaluations per iteration.
    energy_evaluations_per_iteration = energy_evaluations_from_iterations(system_name, n_iterations=1)
    # Find total number of iterations.
    last_iteration = tot_energy_evaluations / energy_evaluations_per_iteration
    return get_iteration_cutoffs(last_iteration)


def cpu_time_iteration_cutoffs(tot_time, system_id, yank_cpu_times):
    """Compute the 100 YANK iterations to use for analysis from the total wall-clock time to consider."""
    # Find average time per iteration.
    yank_time_per_iteration = yank_cpu_times[system_id] / YANK_N_ITERATIONS

    # Find total number of iterations.
    last_iteration = tot_time / yank_time_per_iteration
    return get_iteration_cutoffs(last_iteration)


def get_iteration_cutoffs(last_iteration):
    """Return 100 equally spaced iterations.

    Parameters
    ----------
    last_iteration : float
        The approximate total number of iterations considered.
        This will be rounded to the nearest integer.

    Returns
    -------
    iteration_cutoffs : list of int
        A list of 100 equally-spaced iterations.

    """
    # Find all iterations cutoff.
    first_iteration = last_iteration / 100
    iteration_cutoffs = np.linspace(first_iteration, last_iteration, num=100, endpoint=True)

    # Convert to list of integers.
    iteration_cutoffs = np.rint(iteration_cutoffs).astype(int).tolist()

    assert iteration_cutoffs[-1] <= YANK_N_ITERATIONS
    assert iteration_cutoffs[0] != 0
    return iteration_cutoffs


# =============================================================================
# SAMPLing HOST-GUEST SUBMISSION
# =============================================================================

class SamplingSubmission(SamplSubmission):
    """A submission for the SAMPLing challenge.

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
    CHALLENGE_IDS = {975}

    # The IDs of the submissions used for testing the validation.
    TEST_SUBMISSIONS = {'3rhz6'}

    # Section of the submission file.
    SECTIONS = {'Predictions', 'Cost', 'Name', 'Software', 'TechnicalDetails', 'Method'}

    # Sections in CSV format with kwargs to pass to pandas.read_csv().
    CSV_SECTIONS = {
        # The predictions table is transposed in __init__.
        'Predictions': {'header': None, 'index_col': 0},
        'Cost': {'names': ('System ID', 'N energy evaluations', 'Wall clock time', 'CPU time'),
                 'index_col': 'System ID'}
    }

    def __init__(self, file_path, user_map):
        super().__init__(file_path, user_map)

        file_name = os.path.splitext(os.path.basename(file_path))[0]
        file_data = file_name.split('-')

        # Check this is a known host, and tansform into uppercase.
        self.type = file_data[2]
        assert self.type in ['relative', 'absolute']

        self.file_name, self.index = file_data[3:]
        self.index = int(self.index)

        # Load predictions.
        sections = self._load_sections(file_path)  # From parent-class.
        self.name = sections['Name'][0]
        self.cost = sections['Cost']  # This is a pandas DataFrame.
        self.paper_name = self._assign_paper_name()

        # Reformat predictions to create main data table.
        predictions = sections['Predictions']  # This is a pandas DataFrame.
        data = []
        for system_id, row in predictions.iterrows():
            row = row.values
            free_energies = row[list(range(0, 200, 2))].tolist()
            err_free_energies = row[list(range(1, 200, 2))].tolist()

            tot_energy_evals = self.cost.loc[system_id, 'N energy evaluations']
            energy_evaluations = np.linspace(tot_energy_evals / 100, tot_energy_evals,
                                             num=100, endpoint=True)
            energy_evaluations = np.rint(energy_evaluations).astype(int).tolist()

            if not np.isnan(self.cost.loc[system_id, 'CPU time']):
                tot_time = self.cost.loc[system_id, 'CPU time']
            else:
                tot_time = self.cost.loc[system_id, 'Wall clock time']
            cpu_times = np.linspace(tot_time / 100, tot_time, num=100, endpoint=True)

            # Data row.
            for timepoint_idx in range(100):
                data.append({
                    'System ID': system_id,
                    'System name': system_id[:-2],  # Remove replicate ID.
                    'N energy evaluations': energy_evaluations[timepoint_idx],
                    'CPU time [s]': cpu_times[timepoint_idx],
                    'CPU time [h]': cpu_times[timepoint_idx] / 3600,
                    'CPU time [d]': cpu_times[timepoint_idx] / 3600 / 24,
                    DG_KEY: free_energies[timepoint_idx],
                    DDG_KEY: err_free_energies[timepoint_idx],
                    'Simulation percentage': timepoint_idx + 1,
                })
        self.data = pd.DataFrame(data)

    def mean_free_energies(self):
        """Return a dataframe with mean free energies and 95% t-based confidence intervals."""
        return compute_system_name_mean_free_energies(self.data)

    def _assign_paper_name(self):
        name_table = {
            'APR/pAPRika': 'AMBER/APR',
            'Langevin/Virtual Bond/TI': 'AMBER/TI',
            'WExploreRateRatio': 'OpenMM/WExplore',
            'SOMD/AM1BCC-GAFF-TIP3P/MBAR/C': 'OpenMM/SOMD',
            'Expanded-ensemble/MBAR': 'GROMACS/EE'
        }
        if self.receipt_id == 'NB007':
            return 'GROMACS/CT-NS'
        elif self.receipt_id == 'NB008':
            return 'GROMACS/CT-NS-long'
        if self.receipt_id == 'NB009':
            return 'GROMACS/Jarz-F'
        elif self.receipt_id == 'NB010':
            return 'GROMACS/Jarz-R'
        if self.receipt_id == 'NB011':
            return 'GROMACS/Jarz-F-Gauss'
        elif self.receipt_id == 'NB012':
            return 'GROMACS/Jarz-R-Gauss'
        return name_table[self.name]


# =============================================================================
# YANK ANALYSIS
# =============================================================================

class YankSamplingAnalysis:
    """Utility class to easily access the results of the YANK analysis."""

    def __init__(self, directory_path):
        # Read in analysis files.
        self._yank_free_energies = {}
        pattern_file_path = os.path.join(directory_path, 'yank-*.json')
        for file_path in glob.glob(pattern_file_path):
            file_name = os.path.splitext(os.path.basename(file_path))[0]
            _, system_id = file_name.split('-', 1)

            with open(file_path, 'r') as f:
                # Cast string keys to int.
                self._yank_free_energies[system_id] = {int(k): v for k, v in json.load(f).items()}
        assert len(self._yank_free_energies) == 15

        yank_cpu_times_file_path = os.path.join(directory_path, 'yank_cpu_times.json')
        with open(yank_cpu_times_file_path, 'r') as f:
            self._yank_cpu_times = json.load(f)

        # Attributes exposed by a SamplingSubmission.
        self.name = 'OpenMM/HREX'
        self.paper_name = 'OpenMM/HREX'
        self.receipt_id = 'REF'
        self.file_name = 'YANK'

    @property
    def system_names(self):
        return {system_id[:-2] for system_id in self._yank_free_energies}

    def export(self, file_base_path):
        """Export the YANK analysis into CSV and JSON formats.

        Parameters
        ----------
        file_base_path : str
            The extension-less path where to save the file.
        """
        exported_data = collections.OrderedDict()

        # Export data of 5 replicates
        for system_id in sorted(self._yank_free_energies):
            iterations = sorted(self._yank_free_energies[system_id])
            system_id_data = self._get_free_energies_from_iterations(iterations, [system_id], mean_trajectory=False)
            exported_data[system_id] = collections.OrderedDict([
                ('DG', system_id_data[DG_KEY].values.tolist()),
                ('dDG', system_id_data[DDG_KEY].values.tolist()),
                ('hrex_iterations', system_id_data['HREX iteration'].values.tolist()),
                ('n_energy_evaluations', system_id_data['N energy evaluations'].values.tolist()),
                ('cpu_times', system_id_data['CPU time [s]'].values.tolist()),
            ])

        # Export data of mean trajectory and confidence intervals.
        for system_name in self.system_names:
            system_name_data = self.get_system_free_energies(system_name, mean_trajectory=True)
            exported_data[system_name + '-mean'] = collections.OrderedDict([
                ('DG', system_name_data[DG_KEY].values.tolist()),
                ('DG_CI', system_name_data['$\Delta$G CI'].values.tolist()),
                ('hrex_iterations', system_name_data['HREX iteration'].values.tolist()),
                ('n_energy_evaluations', system_name_data['N energy evaluations'].values.tolist()),
            ])

        # Export.
        export_dictionary(exported_data, file_base_path)

    def get_system_iterations(self, system_id):
        return sorted(self._yank_free_energies[system_id])

    def get_system_free_energies(self, system_name, mean_trajectory=False):
        """Get all the free energies from the system name as a Dataframe."""
        system_ids = sorted([k for k in self._yank_free_energies if k[:-2] == system_name])

        # Find all iterations in common among the system_ids.
        system_name_common_iterations = None
        for system_id in system_ids:
            iterations = sorted(self._yank_free_energies[system_id])
            if system_name_common_iterations is None:
                system_name_common_iterations = set(iterations)
            else:
                system_name_common_iterations.intersection_update(set(iterations))

        # Retrieve the dataframe.
        return self._get_free_energies_from_iterations(sorted(system_name_common_iterations),
                                                       system_ids, mean_trajectory)

    def get_reference_free_energies(self):
        """Get the mean free energy estimate of the last iteration."""
        # We can't use mean_trajectory=True because otherwise
        # _compute_mean_trajectory() will go under infinite recursion.
        last_iteration_free_energies = self._get_free_energies_from_iterations(
            [YANK_N_ITERATIONS], system_ids=[], mean_trajectory=False)

        reference_free_energies = {}
        for system_name in self.system_names:
            data = last_iteration_free_energies[last_iteration_free_energies['System name'] == system_name]
            reference_free_energies[system_name] = np.mean(data[DG_KEY].values)
        return reference_free_energies

    def get_free_energies_from_total_time(self, tot_time, system_id):
        """Get 100 equally-spaced free energies and errors covering tot_time as a DataFrame."""
        iterations = cpu_time_iteration_cutoffs(tot_time, system_id, self._yank_cpu_times)
        return self._get_free_energies_from_iterations(iterations, [system_id])

    def get_free_energies_from_energy_evaluations(self, n_energy_evaluations, system_id=None,
                                                  system_name=None, mean_trajectory=False):
        """Get 100 equally-spaced free energies and errors covering n_energy_evaluations as a DataFrame."""
        assert operator.xor(system_id is not None, system_name is not None)

        if system_name is None:
            system_name = system_id[:-2]
        iterations = energy_evaluations_iteration_cutoffs(n_energy_evaluations, system_name)

        if system_id is not None:
            system_ids = [system_id]
        else:
            system_ids = sorted([k for k in self._yank_free_energies if k[:-2] == system_name])
        return self._get_free_energies_from_iterations(iterations, system_ids, mean_trajectory)

    def get_free_energies_from_iteration(self, final_iteration, system_id=None,
                                         system_name=None, mean_trajectory=False):
        """Get 100 equally-spaced free energies and errors covering final_iterations as a DataFrame."""
        iterations = get_iteration_cutoffs(final_iteration)
        if system_id is None and system_name is None:
            system_ids = []
        elif system_id is not None:
            system_ids = [system_id]
        else:
            system_ids = sorted([k for k in self._yank_free_energies if k[:-2] == system_name])
        return self._get_free_energies_from_iterations(iterations, system_ids, mean_trajectory)

    def _get_free_energies_from_iterations(self, iterations, system_ids, mean_trajectory=False):
        # Handle default argument.
        if len(system_ids) == 0:
            # Pick everything.
            system_ids = sorted(self._yank_free_energies.keys())

        # Create dataframe.
        dataframe = []
        for system_id in system_ids:
            system_name = system_id[:-2]

            for iteration_idx, iteration in enumerate(iterations):
                f, df = self._yank_free_energies[system_id][iteration]
                n_energy_evaluations = energy_evaluations_from_iterations(system_name, iteration)
                cpu_time = self._yank_cpu_times[system_id] / YANK_N_ITERATIONS * iteration

                dataframe.append({
                    'System ID': system_id,
                    'System name': system_name,
                    DG_KEY: f,
                    DDG_KEY: df,
                    'Simulation percentage': iteration_idx + 1,
                    'HREX iteration': iteration,
                    'N energy evaluations': n_energy_evaluations,
                    'CPU time [s]': cpu_time,
                    'CPU time [h]': cpu_time / 3600,
                    'CPU time [d]': cpu_time / 3600 / 24,
                })

        # Convert to Pandas DataFrame.
        columns_order = ['System ID', 'System name', 'Simulation percentage',
                         DG_KEY, DDG_KEY, 'HREX iteration', 'N energy evaluations',
                         'CPU time [s]', 'CPU time [h]', 'CPU time [d]']
        dataframe = pd.DataFrame(dataframe, columns=columns_order)

        # Create average +- CI trajectory if requested.
        if mean_trajectory:
            dataframe = self._compute_mean_trajectory(dataframe)
        return dataframe

    def _compute_mean_trajectory(self, data):
        """Compute average free energy and t-based CI for each iteration."""
        reference_free_energies = self.get_reference_free_energies()
        return compute_system_name_mean_free_energies(data, reference_free_energies=reference_free_energies,
                                                      extra_fields=['HREX iteration'])
