#!/usr/bin/env python

# =============================================================================
# GLOBAL IMPORTS
# =============================================================================

import collections
import copy
import itertools
import json
import math
import os

import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sns
from matplotlib import pyplot as plt
from pkganalysis.sampling import (SamplingSubmission, YankSamplingAnalysis,
                                  YANK_N_ITERATIONS, DG_KEY, DDG_KEY, export_dictionary)
from pkganalysis.submission import (load_submissions)

# =============================================================================
# CONSTANTS
# =============================================================================

YANK_METHOD_PAPER_NAME = 'OpenMM/HREX'

# Paths to input data.
SAMPLING_SUBMISSIONS_DIR_PATH = '../SubmissionsDoNotUpload/975/'
YANK_ANALYSIS_DIR_PATH = 'YankAnalysis/Sampling/'
SAMPLING_ANALYSIS_DIR_PATH = '../SAMPLing/'
SAMPLING_DATA_DIR_PATH = os.path.join(SAMPLING_ANALYSIS_DIR_PATH, 'Data')
SAMPLING_PLOT_DIR_PATH = os.path.join(SAMPLING_ANALYSIS_DIR_PATH, 'Plots')
SAMPLING_PAPER_DIR_PATH = os.path.join(SAMPLING_ANALYSIS_DIR_PATH, 'PaperImages')

# All system ids.
SYSTEM_IDS = [
    'CB8-G3-0', 'CB8-G3-1', 'CB8-G3-2', 'CB8-G3-3', 'CB8-G3-4',
    'OA-G3-0', 'OA-G3-1', 'OA-G3-2', 'OA-G3-3', 'OA-G3-4',
    'OA-G6-0', 'OA-G6-1', 'OA-G6-2', 'OA-G6-3', 'OA-G6-4'
]

# Kelly's colors for maximum contrast.
#               "gray95",  "gray13",  "gold2",   "plum4",   "darkorange1", "lightskyblue2", "firebrick", "burlywood3", "gray51", "springgreen4", "lightpink2", "deepskyblue4", "lightsalmon2", "mediumpurple4", "orange", "maroon", "yellow3", "brown4", "yellow4", "sienna4", "chocolate", "gray19"
KELLY_COLORS = ['#F2F3F4', '#222222', '#F3C300', '#875692', '#F38400',     '#A1CAF1',       '#BE0032',  '#C2B280',   '#848482', '#008856',      '#E68FAC',    '#0067A5',     '#F99379',      '#604E97',      '#F6A600', '#B3446C', '#DCD300', '#882D17', '#8DB600', '#654522', '#E25822', '#2B3D26']
TAB10_COLORS = sns.color_palette('tab10')
# Index of Kelly's colors associated to each submission.
SUBMISSION_COLORS = {
    'AMBER/APR': KELLY_COLORS[11],
    'OpenMM/WExplore': KELLY_COLORS[7],
    'OpenMM/SOMD': KELLY_COLORS[4],
    'GROMACS/EE': KELLY_COLORS[3],
    'GROMACS/EE-fullequil': KELLY_COLORS[10],
    YANK_METHOD_PAPER_NAME: KELLY_COLORS[9],
    'GROMACS/CT-NS-long': KELLY_COLORS[6],
    'GROMACS/CT-NS': KELLY_COLORS[1],
    'GROMACS/NS-Jarz-F': TAB10_COLORS[0],
    'GROMACS/NS-Jarz-R': TAB10_COLORS[1],
    'GROMACS/NS-Gauss-F': TAB10_COLORS[2],
    'GROMACS/NS-Gauss-R': TAB10_COLORS[4],
}

N_ENERGY_EVALUATIONS_SCALE = 1e6

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def reduce_to_first_significant_digit(quantity, uncertainty):
    """Truncate a quantity to the first significant digit of its uncertainty."""
    first_significant_digit = math.floor(math.log10(abs(uncertainty)))
    quantity = round(quantity, -first_significant_digit)
    uncertainty = round(uncertainty, -first_significant_digit)
    return quantity, uncertainty


def load_yank_analysis():
    """Load the YANK analysis in a single dataframe."""
    yank_free_energies = {}
    for system_id in SYSTEM_IDS:
        file_path = os.path.join(YANK_ANALYSIS_DIR_PATH, 'yank-{}.json'.format(system_id))
        with open(file_path, 'r') as f:
            yank_free_energies[system_id] = json.load(f)
    return yank_free_energies


def fit_efficiency(mean_data, find_best_fit=True):
    """Compute the efficiency by fitting the model and using only the asymptotic data.

    We fit using the simulation percentage as the independent value
    because it is less prone to overflowing during fitting. We then
    return the efficiency in units of (kcal/mol)**2/n_energy_evaluations.
    """
    from scipy.optimize import curve_fit

    def model(x, log_efficiency):
        return np.exp(log_efficiency) / x

    vars = mean_data['std'].values**2

    cost = mean_data['Simulation percentage'].values
    # cost = mean_data['N energy evaluations'].values / 1e7
    if find_best_fit:
        # Find fit with best error up to discarding 70% of calculation.
        max_discarded = math.floor(0.5*len(cost))
    else:
        # Use all the data.
        max_discarded = 1

    # Fit.
    fits = []
    for n_discarded in range(max_discarded):
        cost_fit = cost[n_discarded:]
        vars_fit = vars[n_discarded:]
        fit = curve_fit(model, cost_fit, vars_fit, p0=[0.0])
        fits.append((np.exp(fit[0]), fit[1]))

    # Find the fit with the minimum error.
    n_discarded = fits.index(min(fits, key=lambda x: x[1]))
    # Convert efficiency / simulation_percentage to efficiency / n_energy_evaluations
    efficiency = fits[n_discarded][0][0] / 100 * mean_data['N energy evaluations'].values[-1]
    # efficiency = fits[n_discarded][0][0] * 1e7
    return efficiency, n_discarded


def export_submissions(submissions, reference_free_energies):
    """Export the submission data to CSV and JSON format."""
    for submission in submissions:
        exported_data = {}

        # Export data of the 5 independent replicates.
        for system_id in sorted(submission.data['System ID'].unique()):
            system_id_data = submission.data[submission.data['System ID'] == system_id]
            exported_data[system_id] = collections.OrderedDict([
                ('DG', system_id_data[DG_KEY].values.tolist()),
                ('dDG', system_id_data[DDG_KEY].values.tolist()),
                ('cpu_times', system_id_data['CPU time [s]'].values.tolist()),
                ('n_energy_evaluations', system_id_data['N energy evaluations'].values.tolist()),
            ])

        # Export data of mean trajectory and confidence intervals.
        mean_free_energies = submission.mean_free_energies()
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


# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

def plot_mean_free_energy(mean_data, ax, x='Simulation percentage',
                          color_mean=None, color_ci=None, zorder=None,
                          start=None, stride=1, scale_n_energy_evaluations=True,
                          plot_ci=True, **plot_kwargs):
    """Plot mean trajectory with confidence intervals."""
    ci_key = '$\Delta$G CI'

    if start is None:
        # Discard the first datapoint which are 0.0 (i.e. no estimate).
        start = np.nonzero(mean_data[DG_KEY].values)[0][0]

    if x == 'N energy evaluations' and scale_n_energy_evaluations:
        # Plot in millions of energy evaluations.
        scale = N_ENERGY_EVALUATIONS_SCALE
    else:
        scale = 1

    x = mean_data[x].values[start::stride] / scale
    mean_dg = mean_data[DG_KEY].values[start::stride]
    sem_dg = mean_data[ci_key].values[start::stride]

    # Plot mean trajectory confidence intervals.
    if plot_ci:
        ax.fill_between(x, mean_dg + sem_dg, mean_dg - sem_dg, alpha=0.15, color=color_ci, zorder=zorder)

    # Plot the mean free energy trajectory.
    if zorder is not None:
        # Push the CI shaded area in the background so that the trajectories are always visible.
        zorder += 20
    ax.plot(x, mean_dg, color=color_mean, alpha=1.0, zorder=zorder, **plot_kwargs)
    return ax


def plot_mean_data(mean_data, axes, color=None, label=None, x='N energy evaluations',
                   zorder=None, plot_std=True, plot_bias=True, plot_ci=True):
    """Plot free energy, variance and bias as a function of the cost in three different axes."""
    # Do not plot the part of data without index.
    first_nonzero_idx = np.nonzero(mean_data[DG_KEY].values)[0][0]

    # If the x-axis is the number of energy/force evaluations, plot it in units of millions.
    if x == 'N energy evaluations':
        scale = N_ENERGY_EVALUATIONS_SCALE
    else:
        scale = 1

    # Plot the submission mean trajectory with CI.
    plot_mean_free_energy(mean_data, x=x, ax=axes[0],
                          color_mean=color, color_ci=color, zorder=zorder,
                          start=first_nonzero_idx, label=label, plot_ci=plot_ci)

    # Plot standard deviation of the trajectories.
    if plot_std:
        axes[1].plot(mean_data[x].values[first_nonzero_idx:] / scale,
                     mean_data['std'].values[first_nonzero_idx:], color=color, alpha=0.8,
                     zorder=zorder, label=label)
    if plot_bias:
        axes[2].plot(mean_data[x].values[first_nonzero_idx:] / scale,
                     mean_data['bias'].values[first_nonzero_idx:], color=color, alpha=0.8,
                     zorder=zorder, label=label)


def align_yaxis(ax1, v1, ax2, v2):
    """Adjust ax2 ylimit so that v2 in in the twin ax2 is aligned to v1 in ax1.

    From https://stackoverflow.com/questions/10481990/matplotlib-axis-with-two-scales-shared-origin .

    """
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny+dy, maxy+dy)


# =============================================================================
# FIGURE 2
# =============================================================================

def plot_submissions_trajectory(submissions, yank_analysis, axes, y_limits=None,
                                plot_std=True, plot_bias=True, plot_bias_to_reference=False,
                                system_names=None):
    """Plot free energy trajectories, std, and bias of the given submissions."""
    if system_names is None:
        system_names = ['CB8-G3', 'OA-G3', 'OA-G6']
    n_systems = len(system_names)
    max_n_energy_evaluations = {system_name: 0 for system_name in system_names}
    min_n_energy_evaluations = {system_name: np.inf for system_name in system_names}

    # Handle default arguments.
    if y_limits is None:
        # 3 by 3 matrix of y limits for the plots.
        y_limits = [[None for _ in range(n_systems)] for _ in range(n_systems)]
    # We need a 2D array of axes for the code to work even if we're not plotting std or bias.
    try:
        axes_shape = len(axes.shape)
    except AttributeError:
        axes = np.array([[axes]])
    else:
        if axes_shape == 1:
            axes = np.array([axes])

    # Build a dictionary mapping submissions and system names to their mean data.
    all_mean_data = {}
    for submission in submissions:
        # We always want to print in order
        all_mean_data[submission.paper_name] = {}
        mean_free_energies = submission.mean_free_energies()

        for system_name in system_names:
            # CB8-G3 calculations for GROMACS/EE did not converge.
            if submission.name == 'Expanded-ensemble/MBAR' and system_name == 'CB8-G3':
                continue

            # Add mean free energies for this system.
            system_mean_data = mean_free_energies[mean_free_energies['System name'] == system_name]
            n_energy_evaluations = system_mean_data['N energy evaluations'].values[-1]
            # # TODO REMOVE ME
            # if submission.paper_name == 'GROMACS/CT-NS' and 'OA' in system_name:
            #     max_n_energy_evaluations[system_name] = n_energy_evaluations
            #     continue
            # elif 'APR' not in submission.paper_name:
            #     continue

            all_mean_data[submission.paper_name][system_name] = system_mean_data

            # Keep track of the maximum and minimum number of energy evaluations,
            # which will be used to determine how to truncate the plotted reference
            # data and determine the zorder of the trajectories respectively.
            max_n_energy_evaluations[system_name] = max(max_n_energy_evaluations[system_name],
                                                        n_energy_evaluations)
            min_n_energy_evaluations[system_name] = min(min_n_energy_evaluations[system_name],
                                                        n_energy_evaluations)

    # Add also reference YANK calculations if provided.
    if yank_analysis is not None:
        all_mean_data[YANK_METHOD_PAPER_NAME] = {}
        for system_name in system_names:
            system_mean_data = yank_analysis.get_free_energies_from_energy_evaluations(
                max_n_energy_evaluations[system_name], system_name=system_name, mean_trajectory=True)
            all_mean_data[YANK_METHOD_PAPER_NAME][system_name] = system_mean_data

    # Create a table mapping submissions and system name to the zorder used
    # to plot the free energy trajectory so that smaller shaded areas are on
    # top of bigger ones.
    # First find the average CI for all methods up to min_n_energy_evaluations.
    methods_cis = {name: {} for name in system_names}
    for method_name, method_mean_data in all_mean_data.items():
        for system_name, system_mean_data in method_mean_data.items():
            # Find index of all energy evaluations < min_n_energy_evaluations.
            n_energy_evaluations = system_mean_data['N energy evaluations'].values
            last_idx = np.searchsorted(n_energy_evaluations, min_n_energy_evaluations[system_name], side='right')
            cis = system_mean_data['$\Delta$G CI'].values[:last_idx]
            methods_cis[system_name][method_name] = np.mean(cis)

    # For each system, order methods from smallest CI (plot on top) to greatest CI (background).
    zorders = {name: {} for name in system_names}
    for system_name, system_cis in methods_cis.items():
        ordered_methods = sorted(system_cis.keys(), key=lambda method_name: system_cis[method_name])
        for zorder, method_name in enumerate(ordered_methods):
            zorders[system_name][method_name] = zorder

    # The columns are in order CB8-G3, OA-G3, and OA-G6.
    system_columns = {'CB8-G3': 0, 'OA-G3': 1, 'OA-G6': 2}

    # Plot submissions in alphabetical order to order he legend labels.
    for method_name in sorted(all_mean_data.keys()):
        submission_mean_data = all_mean_data[method_name]
        submission_color = SUBMISSION_COLORS[method_name]

        # Plot free energy trajectories.
        for system_name, mean_data in submission_mean_data.items():
            ax_idx = system_columns[system_name]

            # The OA prediction of the NS short protocol are the same of the long protocol submission file.
            if method_name == 'GROMACS/CT-NS-long' and system_name != 'CB8-G3':
                # Just add the label.
                axes[0][ax_idx].plot([], color=submission_color, label=method_name)
                continue

            # Update maximum number of energy evaluations.
            n_energy_evaluations = mean_data['N energy evaluations'].values[-1]
            max_n_energy_evaluations[system_name] = max(max_n_energy_evaluations[system_name],
                                                        n_energy_evaluations)

            # Determine zorder and plot.
            zorder = zorders[system_name][method_name]
            plot_mean_data(mean_data, axes[:,ax_idx], color=submission_color,
                           zorder=zorder, label=method_name,
                           plot_std=plot_std, plot_bias=plot_bias)

            # Plot adding full cost of Wang-Landau equilibration.
            if 'EE' in method_name:
                first_nonzero_idx = np.nonzero(mean_data[DG_KEY].values)[0][0]
                calibration_cost = mean_data['N energy evaluations'].values[first_nonzero_idx] * 4
                mean_data['N energy evaluations'] += calibration_cost
                label = method_name + '-fullequil'
                plot_mean_data(mean_data, axes[:,ax_idx], color=SUBMISSION_COLORS[label],
                               zorder=zorder, label=label, plot_std=plot_std, plot_bias=plot_bias)

    # Fix labels.
    axes[0][0].set_ylabel('$\Delta$G [kcal/mol]')
    if plot_std:
        axes[1][0].set_ylabel('std($\Delta$G) [kcal/mol]')
    if plot_bias:
        axes[2][0].set_ylabel('bias [kcal/mol]')
    central_column_idx = int(len(axes[0])/2)
    axes[-1][central_column_idx].set_xlabel('number of energy/force evaluations [10$^6$]')

    # Fix axes limits.
    for ax_idx, system_name in enumerate(system_names):
        for row_idx in range(len(axes)):
            ax = axes[row_idx][ax_idx]
            # Set the x-axis limits.
            ax.set_xlim((0, max_n_energy_evaluations[system_name]/N_ENERGY_EVALUATIONS_SCALE))
            # Keep the x-axis label only at the bottom row.
            if row_idx != len(axes)-1:
                ax.xaxis.set_ticklabels([])
            y_lim = y_limits[row_idx][ax_idx]
            if y_lim is not None:
                ax.set_ylim(y_lim)

        # Set the system name in the title.
        axes[0][ax_idx].set_title(system_name)

    # Create a bias axis AFTER the ylim has been set.
    if yank_analysis is not None and plot_bias_to_reference:
        for ax_idx, (system_name, ax) in enumerate(zip(system_names, axes[0])):
            yank_full_mean_data = yank_analysis.get_system_free_energies(system_name, mean_trajectory=True)
            ref_free_energy = yank_full_mean_data[DG_KEY].values[-1]
            with sns.axes_style('white'):
                ax2 = ax.twinx()
                # Plot a vertical line to fix the scale.
                vertical_line = np.linspace(*ax.get_ylim()) - ref_free_energy
                ax2.plot([50] * len(vertical_line), vertical_line, alpha=0.0001)
                ax2.grid(alpha=0.5, linestyle='dashed', zorder=0)
                # We add the bias y-label only on the rightmost Axis.
                if ax_idx == n_systems - 1:
                    ax2.set_ylabel('Bias to reference [kcal/mol]')
                # Set the 0 of the twin axis to the YANK reference free energy.
                align_yaxis(ax, ref_free_energy, ax2, 0.0)


def plot_all_entries_trajectory(submissions, yank_analysis, zoomed=False):
    """Plot free energy trajectories, std, and bias of the challenge entries."""
    # Create a figure with 3 columns (one for each system) and 2 rows.
    # The first row contains the free energy trajectory and CI, the second
    # a plot of the estimator variance, and the third the bias to the
    # asymptotic value.
    if zoomed:
        # figsize = (7.25, 6.2)  # Without WExplorer
        figsize = (7.25, 7.5)
    else:
        figsize = (7.25, 7.5)  # With WExplorer
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=figsize)

    # Optionally, remove WExplore.
    if zoomed:
        submissions = [s for s in submissions if s.name not in ['WExploreRateRatio']]

    if zoomed:
        # Y-axis limits when WExplore calculations are excluded.
        y_limits = [
            [(-15, -10), (-9, -4), (-9, -4)],
            [(0, 2), (0, 0.8), (0, 0.8)],
            [(-3, 1), (-0.6, 0.6), (-0.6, 0.6)],
        ]
    else:
        # Y-axis limits when WExplore calculations are included.
        y_limits = [
            [(-17, -9), (-13, -5), (-13, -5)],
            [(0, 2), (0, 1.75), (0, 1.75)],
            [(-4, 4), (-0.6, 0.6), (-0.6, 0.6)],
        ]

    plot_submissions_trajectory(submissions, yank_analysis, axes, y_limits=y_limits)

    # Show/save figure.
    if zoomed:
        plt.tight_layout(h_pad=0.2, rect=[0.0, 0.00, 1.0, 0.92], w_pad=0.0)  # Without WExplorer
    else:
        plt.tight_layout(h_pad=0.2, rect=[0.0, 0.00, 1.0, 0.92])  # With WExplorer

    # Plot legend.
    if zoomed:
        # bbox_to_anchor = (2.52, 1.55)  # Without WExplorer.
        bbox_to_anchor = (2.4, 1.48)
    else:
        bbox_to_anchor = (2.4, 1.48)  # With WExplorer.
    axes[0][1].legend(loc='upper right', bbox_to_anchor=bbox_to_anchor,
                      fancybox=True, ncol=4)
    plt.subplots_adjust(wspace=0.35)
    # plt.show()
    if zoomed:
        file_name = 'Figure2-free_energy_trajectories_zoomed'
    else:
        file_name = 'Figure2-free_energy_trajectories'
    figure_dir_path = os.path.join(SAMPLING_PAPER_DIR_PATH, 'Figure2-free_energy_trajectories')
    os.makedirs(figure_dir_path, exist_ok=True)
    output_base_path = os.path.join(figure_dir_path, file_name)
    plt.savefig(output_base_path + '.pdf')
    # plt.savefig(output_base_path + '.png', dpi=500)


def plot_all_nonequilibrium_switching(submissions):
    """Plot free energy trajectories, std, and bias of the nonequilibrium-switching calculations."""
    # Create a figure with 3 columns (one for each system) and 2 rows.
    # The first row contains the free energy trajectory and CI, the second
    # a plot of the estimator variance, and the third the bias to the
    # asymptotic value.
    figsize = (7.25, 3.5)  # With WExplorer
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=figsize)

    # Select nonequilibrium-switching calculations with estimators.
    submissions = [s for s in submissions if 'NS' in s.paper_name]

    # Y-axis limits.
    y_limits = [
        [(-20, 5), (-40, 0), (-40, 0)]
    ]

    plot_submissions_trajectory(submissions, yank_analysis=None, axes=axes,
                                y_limits=y_limits, plot_std=False, plot_bias=False)

    # Show/save figure.
    plt.tight_layout(pad=0.0, rect=[0.0, 0.00, 1.0, 0.85])

    # Plot legend.
    axes[0].legend(loc='upper left', bbox_to_anchor=(0.1, 1.3),
                   fancybox=True, ncol=3)
    plt.subplots_adjust(wspace=0.35)

    # plt.show()
    figure_dir_path = os.path.join(SAMPLING_PAPER_DIR_PATH, 'Figure3-nonequilibrium_comparison')
    os.makedirs(figure_dir_path, exist_ok=True)
    output_base_path = os.path.join(figure_dir_path, 'Figure3-nonequilibrium_comparison')
    plt.savefig(output_base_path + '.pdf')
    # plt.savefig(output_base_path + '.png', dpi=500)


# =============================================================================
# FIGURE 3
# =============================================================================

# Directories containing the volume information of YANK and GROMACS/EE.
BAROSTAT_DATA_DIR_PATH = os.path.join('..', 'SAMPLing', 'Data', 'BarostatData')
YANK_VOLUMES_DIR_PATH = os.path.join(BAROSTAT_DATA_DIR_PATH, 'YankVolumes')
EE_VOLUMES_DIR_PATH = os.path.join(BAROSTAT_DATA_DIR_PATH, 'EEVolumes')


def plot_volume_distributions(axes, plot_predicted=False):
    """Plot the volume distributions obtained with Monte Carlo and Berendsen barostat."""
    import scipy.stats
    import scipy.integrate
    from simtk import unit

    # Load data.
    yank_volumes = collections.OrderedDict([
        (1, np.load(os.path.join(YANK_VOLUMES_DIR_PATH, 'volumes_pressure100.npy'))),
        (100, np.load(os.path.join(YANK_VOLUMES_DIR_PATH, 'volumes_pressure10000.npy'))),
    ])

    ee_volumes = collections.OrderedDict([
        (1, np.load(os.path.join(EE_VOLUMES_DIR_PATH, '1atm_vanilla.npy'))),
        (100, np.load(os.path.join(EE_VOLUMES_DIR_PATH, '100atm_vanilla.npy'))),
    ])

    titles = ['Monte Carlo barostat', 'Berendsen barostat']
    for ax, volume_trajectories, title in zip(axes, [yank_volumes, ee_volumes], titles):
        for pressure, trajectory in volume_trajectories.items():
            label = '$\\rho$(V|{}atm)'.format(pressure)
            print('{}: mean={:.3f}nm^3, var={:.3f}'.format(label, np.mean(trajectory),
                                                           np.var(trajectory)))
            ax = sns.distplot(trajectory, label=label, hist=True, ax=ax)

        if plot_predicted:
            # Plot predicted distribution.
            beta = 1.0 / (unit.BOLTZMANN_CONSTANT_kB * 298.15*unit.kelvin)
            p1 = 1.0 * unit.atmosphere
            p2 = 100.0 * unit.atmosphere
            volumes = np.linspace(78.0, 82.0, num=200)
            fit = scipy.stats.norm

            # Fit the original distribution.
            original_pressure, new_pressure = list(volume_trajectories.keys())
            original_trajectory = list(volume_trajectories.values())[0]
            fit_parameters = fit.fit(original_trajectory)

            # Find normalizing constant predicted distribution.
            predicted_distribution = lambda v: np.exp(-beta*(p2 - p1)*v*unit.nanometer**3) * fit.pdf([v], *fit_parameters)
            normalizing_factor = scipy.integrate.quad(predicted_distribution, volumes[0], volumes[-1])[0]
            predicted = np.array([predicted_distribution(v) / normalizing_factor for v in volumes])

            # Set the scale.
            label = '$\\rho$(V|{}atm)$\cdot e^{{\\beta ({}atm - {}atm) V}}$'.format(original_pressure, new_pressure, original_pressure)
            ax.plot(volumes, predicted, label=label)
            # ax.plot(volumes, [fit.pdf([v], *fit_parameters) for v in volumes], label='original')
            ax.set_ylabel('density')

        ax.set_title(title + ' volume distribution')
        ax.legend(fontsize='xx-small')
        ax.set_xlim((78.5, 82.0))
        ax.set_xlabel('Volume [nm^3]')


# Directory with the restraint information.
RESTRAINT_DATA_DIR_PATH = os.path.join('YankAnalysis', 'RestraintAnalysis')

# The state index of the discharged state with LJ interactions intact.
DISCHARGED_STATE = {
    'CB8-G3': 25,
    'OA-G3': 32,
    'OA-G6': 29
}

# The final free energy predictions without restraint unbiasing.
BIASED_FREE_ENERGIES = {
    'CB8-G3-0': -10.643,
    'CB8-G3-1': -10.533,
    'CB8-G3-2': -10.463,
    'CB8-G3-3': None,  # TODO: Run the biased analysis
    'CB8-G3-4': -10.324,
    'OA-G3-0': -5.476,
    'OA-G3-1': -5.588,
    'OA-G3-2': -5.486,
    'OA-G3-3': -5.510,
    'OA-G3-4': -5.497,
    'OA-G6-0': -5.669,
    'OA-G6-1': -5.665,
    'OA-G6-2': -5.767,
    'OA-G6-3': -5.737,
    'OA-G6-4': -5.788,
}


def plot_restraint_distance_distribution(system_id, ax, kde=True, iteration_set=None):
    """Plot the distribution of restraint distances at bound, discharged, and decoupled states.

    Return the 99.99-percentile restraint radius that was used as a cutoff during analysis.
    """
    n_iterations = YANK_N_ITERATIONS + 1  # Count also iteration 0.
    system_name = system_id[:-2]
    discharged_state_idx = DISCHARGED_STATE[system_name]

    # Load all distances cached during the analysis.
    cache_dir_path = os.path.join('pkganalysis', 'cache', system_id.replace('-', ''))
    cached_distances_file_path = os.path.join(cache_dir_path, 'restraint_distances_cache.npz')
    distances_kn = np.load(cached_distances_file_path)['arr_0']
    # Distances are in nm but we plot in Angstrom.
    distances_kn *= 10
    n_states = int(len(distances_kn) / n_iterations)

    # Use the same colors that are used in the water analysis figures.
    color_palette = sns.color_palette('viridis', n_colors=n_states)
    color_palette = [color_palette[i] for i in (0, discharged_state_idx, -1)]

    # Isolate distances in the bound, discharged (only LJ), and decoupled state.
    distances_kn_bound = distances_kn[:n_iterations]
    distances_kn_discharged = distances_kn[(discharged_state_idx-1)*n_iterations:discharged_state_idx*n_iterations]
    distances_kn_decoupled = distances_kn[(n_states-1)*n_iterations:]

    # Filter iterations.
    if iteration_set is not None:
        distances_kn_bound = distances_kn_bound[iteration_set]
        distances_kn_discharged = distances_kn_discharged[iteration_set]
        distances_kn_decoupled = distances_kn_decoupled[iteration_set]
    assert len(distances_kn_bound) == len(distances_kn_decoupled)

    # Plot the distributions.
    # sns.distplot(distances_kn, ax=ax, kde=True, label='all states')
    sns.distplot(distances_kn_bound, ax=ax, kde=kde, label='bound', color=color_palette[0])
    sns.distplot(distances_kn_discharged, ax=ax, kde=kde, label='discharged', color=color_palette[1])
    sns.distplot(distances_kn_decoupled, ax=ax, kde=kde, label='decoupled', color=color_palette[2])

    # Plot the threshold used for analysis, computed as the
    # 99.99-percentile of all distances in the bound state.
    distance_cutoff = np.percentile(a=distances_kn_bound, q=99.99)
    limits = ax.get_ylim()
    ax.plot([distance_cutoff for _ in range(100)],
            np.linspace(limits[0], limits[1]/2, num=100), color='black')

    return distance_cutoff


def plot_restraint_profile(system_id, ax, restraint_cutoff):
    """Plot the free energy as a function of the restraint cutoff."""
    # Load the free energy profile for this system.
    restraint_profile_file_path = os.path.join(RESTRAINT_DATA_DIR_PATH,
                                               system_id.replace('-', '') + '.json')
    with open(restraint_profile_file_path, 'r') as f:
        free_energies_profile = json.load(f)

    # Reorder the free energies by increasing cutoff and convert str keys to floats.
    free_energies_profile = [(float(d), f) for d, f in free_energies_profile.items()]
    free_energies_profile = sorted(free_energies_profile, key=lambda x: x[0])
    distance_cutoffs, free_energies = list(zip(*free_energies_profile))
    f, df = list(zip(*free_energies))

    # Convert string to floats.
    distance_cutoffs = [float(c) for c in distance_cutoffs]

    # Plot profile.
    ax.errorbar(x=distance_cutoffs, y=f, yerr=df, label='after reweighting')
    # Plot biased free energy
    biased_f = BIASED_FREE_ENERGIES[system_id]
    x = np.linspace(*ax.get_xlim())
    ax.plot(x, [biased_f for _ in x], label='before reweighting')

    # Plot restraint distance cutoff.
    limits = ax.get_ylim()
    x = [restraint_cutoff for _ in range(100)]
    y = np.linspace(limits[0], limits[1], num=100)
    ax.plot(x, y, color='black')


def plot_restraint_analysis(system_id, axes):
    """Plot distribution of restraint distances and free energy profile on two axes."""
    # Histograms of restraint distances/energies.
    ax = axes[0]
    kde = True
    restraint_cutoff = plot_restraint_distance_distribution(system_id, ax, kde=kde)
    # Set restraint distance distribution lables and titles.
    ax.set_title('Harmonic restraint radius distribution')
    if kde is False:
        ax.set_ylabel('Number of samples')
    else:
        ax.set_ylabel('density')
    ax.legend(loc='upper right', fontsize='xx-small')
    ax.set_xlabel('Restraint radius [A]')

    # Free energy as a function of restraint distance.
    ax = axes[1]
    ax.set_title('$\Delta G$ as a function of restraint radius cutoff')
    plot_restraint_profile(system_id, ax, restraint_cutoff)
    # Labels and legend.
    ax.set_xlabel('Restraint radius cutoff [A]')
    ax.set_ylabel('$\Delta G$ [kcal/mol]')
    ax.legend(fontsize='xx-small')


def plot_restraint_and_barostat_analysis():
    """Plot the Figure showing info for the restraint and barostat analysis."""
    import seaborn as sns
    from matplotlib import pyplot as plt
    sns.set_style('whitegrid')
    sns.set_context('paper')

    # Create two columns, each of them share the x-axis.
    fig = plt.figure(figsize=(7.25, 5))
    # Restraint distribution axes.
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(223, sharex=ax1)
    barostat_axes = [ax1, ax2]
    # Volume distribution axes.
    ax3 = fig.add_subplot(222)
    ax4 = fig.add_subplot(224, sharex=ax3)
    restraint_axes = [ax3, ax4]

    # Plot barostat analysis.
    plot_volume_distributions(barostat_axes, plot_predicted=True)

    # Plot restraint analysis.
    system_id = 'OA-G3-0'
    plot_restraint_analysis(system_id, restraint_axes)
    # Configure axes.
    restraint_axes[0].set_xlim((0, 10.045))

    for ax in restraint_axes + barostat_axes:
        ax.tick_params(axis='x', which='major', pad=0.2)
        ax.tick_params(axis='y', which='major', pad=0.2)
    plt.tight_layout(pad=0.5)

    # plt.show()
    output_file_path = os.path.join(SAMPLING_PAPER_DIR_PATH, 'Figure3-restraint_barostat',
                                    'restraint_barostat.pdf')
    os.makedirs(os.path.dirname(output_file_path), exist_ok=True)
    plt.savefig(output_file_path)


# =============================================================================
# FIGURE 4
# =============================================================================

def plot_yank_system_bias(system_name, data_dir_paths, axes, shift_to_origin=True, plot_std=True):
    """Plot the YANK free energy trajectoies when discarding initial samples for a single system."""
    color_palette = sns.color_palette('viridis', n_colors=len(data_dir_paths)+1)

    # Plot trajectories with truncated data.
    all_iterations = set()
    for data_idx, data_dir_path in enumerate(data_dir_paths):
        yank_analysis = YankSamplingAnalysis(data_dir_path)

        # In the YankAnalysis folder, each analysis starting from
        # iteration N is in the folder "iterN/".
        last_dir_name = os.path.basename(os.path.normpath(data_dir_path))
        label = last_dir_name[4:]
        # First color is for the full data.
        color = color_palette[data_idx+1]

        # Collect all iterations that we'll plot for the full data.
        mean_data = yank_analysis.get_system_free_energies(system_name, mean_trajectory=True)
        all_iterations.update(mean_data['HREX iteration'].values.tolist())

        # Simulate plotting starting from the origin.
        if shift_to_origin:
            mean_data['HREX iteration'] -= mean_data['HREX iteration'].values[0]

        plot_mean_data(mean_data, axes, x='HREX iteration', color=color,
                       label=label, plot_std=plot_std, plot_bias=False, plot_ci=False)

    # Plot trajectory with full data.
    color = color_palette[0]

    # Plot an early iteration and all the iterations analyzed for the bias.
    yank_analysis = YankSamplingAnalysis(YANK_ANALYSIS_DIR_PATH)
    system_ids = [system_name + '-' + str(i) for i in range(5)]
    first_iteration = yank_analysis.get_system_iterations(system_ids[0])[2]
    iterations = [first_iteration] + sorted(all_iterations)
    mean_data = yank_analysis._get_free_energies_from_iterations(
        iterations, system_ids, mean_trajectory=True)

    # Simulate plotting starting from the origin.
    if shift_to_origin:
        mean_data['HREX iteration'] -= mean_data['HREX iteration'].values[0]

    # Simulate ploatting starting from the origin.
    plot_mean_data(mean_data, axes, x='HREX iteration', color=color,
                   label='0', plot_std=plot_std, plot_bias=False, plot_ci=False)
    axes[0].set_title(system_name)


def plot_yank_bias(plot_std=True, figure_dir_path=None):
    """Plot YANK free energy trajectories when discarding initial samples."""
    # In the first column, plot the "unshifted" trajectory of CB8-G3,
    # with all sub-trajectories shifted to the origin. In the second
    # and third columns, plot the trajectories of CB8-G3 and OA-G3
    # with all sub-trajectories shifted to the origin.
    what_to_plot = [
        ('CB8-G3', False),
        ('CB8-G3', True),
        ('OA-G3', True),
        # ('OA-G3', False),
        # ('OA-G6', False),
    ]

    if plot_std:
        n_rows = 2
    else:
        n_rows = 1
    n_cols = len(what_to_plot)
    fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(7.25, 4.6))

    # The loops are based on a two dimensional array of axes.
    if n_rows == 1:
        axes = np.array([axes])

    # Sort paths by how many samples they have.
    data_dir_paths = ['YankAnalysis/BiasAnalysis/iter{}/'.format(i) for i in [1000, 2000, 4000, 8000, 16000, 24000]]

    for column_idx, (system_name, shift_to_origin) in enumerate(what_to_plot):
        plot_yank_system_bias(system_name, data_dir_paths, axes[:,column_idx],
                              shift_to_origin=shift_to_origin, plot_std=plot_std)
        title = system_name + ' (shifted)' if shift_to_origin else system_name
        axes[0,column_idx].set_title(title)

    # Fix axes limits and labels.
    ylimits = {
        'CB8-G3': (-12.5, -10.5),
        'OA-G3': (-8, -6),
        'OA-G6': (-8, -6)
    }
    for column_idx, (system_name, _) in enumerate(what_to_plot):
        axes[0][column_idx].set_ylim(ylimits[system_name])
    if plot_std:
        for column_idx in range(3):
            axes[1][column_idx].set_ylim((0, 0.6))

    for row_idx, ax_idx in itertools.product(range(n_rows), range(n_cols)):
        # Control the number of ticks for the x axis.
        axes[row_idx][ax_idx].locator_params(axis='x', nbins=4)
        # Set x limits for number of iterations.
        axes[row_idx][ax_idx].set_xlim((0, YANK_N_ITERATIONS))
    # Remove ticks labels that are shared with the last row.
    for row_idx, ax_idx in itertools.product(range(n_rows-1), range(n_cols)):
        axes[row_idx][ax_idx].set_xticklabels([])

    # Set axes labels.
    axes[0][0].set_ylabel('$\Delta$G [kcal/mol]')
    if plot_std:
        axes[1][0].set_ylabel('std($\Delta$G) [kcal/mol]')
    axes[-1][1].set_xlabel('HREX iteration')

    plt.tight_layout(h_pad=0.1, rect=[0.0, 0.00, 1.0, 0.92])

    handles, labels = axes[0][0].get_legend_handles_labels()
    handles = [handles[-1]] + handles[:-1]
    labels = [labels[-1]] + labels[:-1]
    bbox_to_anchor = (-0.1, 1.45)
    axes[0][0].legend(handles, labels, loc='upper left', bbox_to_anchor=bbox_to_anchor,
                      title='n discarded initial iterations', ncol=len(data_dir_paths)+1,
                      fancybox=True)

    # plt.show()
    if figure_dir_path is None:
        figure_dir_path = os.path.join(SAMPLING_PAPER_DIR_PATH, 'Figure4-bias_hrex')
    os.makedirs(figure_dir_path, exist_ok=True)
    output_file_path = os.path.join(figure_dir_path, 'Figure4-bias_hrex')
    plt.savefig(output_file_path + '.pdf')
    # plt.savefig(output_file_path + '.png', dpi=600)


def plot_equilibration_methods():
    """Plot the different trajectories obtained through the reduced potential and instantaneous work equilibration."""
    sns.set_context('paper', font_scale=1.2)

    yank_analysis_potential = YankSamplingAnalysis(YANK_ANALYSIS_DIR_PATH)
    yank_analysis_work = YankSamplingAnalysis('YankAnalysis/Sampling_instantaneouswork/')

    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(7.25, 8))

    ylimits = {
        'CB8-G3': (-12.5, -10.5),
        'OA-G3': (-8, -6),
        'OA-G6': (-8, -6)
    }

    palette = sns.color_palette('tab10', n_colors=2)

    for col_idx, system_name in enumerate(['CB8-G3', 'OA-G3', 'OA-G6']):
        for i, (label, yank_analysis) in enumerate([
            ('potential', yank_analysis_potential),
            ('average instantaneous work', yank_analysis_work)
        ]):
            # mean_data = yank_analysis.get_system_free_energies(system_name, mean_trajectory=True)
            mean_data = yank_analysis.get_free_energies_from_iteration(
                YANK_N_ITERATIONS, system_name=system_name, mean_trajectory=True)
            plot_mean_data(mean_data, axes[:,col_idx], x='HREX iteration', color=palette[i], label=label,
                           plot_std=True, plot_bias=True, plot_ci=True)

            # Set labels and limits.
            ax = axes[0,col_idx]
            ax.set_title(system_name)
            ax.set_ylim(ylimits[system_name])

    axes[0][0].set_ylabel('$\Delta$G [kcal/mol]')
    axes[1][0].set_ylabel('std($\Delta$G) [kcal/mol]')
    axes[2][0].set_ylabel('bias [kcal/mol]')
    axes[-1,1].set_xlabel('N HREX iterations')

    axes[0,0].legend(loc='upper right', bbox_to_anchor=(2.4, 1.48),
                      fancybox=True, ncol=2)

    # plt.show()
    figure_dir_path = os.path.join(SAMPLING_PAPER_DIR_PATH, 'Figure4-bias_hrex')
    os.makedirs(figure_dir_path, exist_ok=True)
    output_file_path = os.path.join(figure_dir_path, 'Figure4-instantaneous_work_equil')
    plt.savefig(output_file_path + '.pdf')
    # plt.savefig(output_file_path + '.png', dpi=600)

# =============================================================================
# RELATIVE EFFICIENCY ANALYSIS
# =============================================================================

def get_relative_efficiency_input(submission, yank_analysis, system_name):
    """Prepare the data to compute the mean relative efficiencies for this system."""
    # For GROMACS/EE-fullquil we need to account for the extra equilibration
    # cost and shift all energy evaluation to the right.
    if submission.paper_name == 'GROMACS/EE-fullequil':
        mean_free_energies = submission.mean_free_energies()
        mean_data = mean_free_energies[mean_free_energies['System name'] == system_name]
        first_shifted = mean_data['N energy evaluations'].values[0]
        last_shifted = mean_data['N energy evaluations'].values[-1]
        calibration_cost = first_shifted*100/99 - last_shifted/99
    else:
        calibration_cost = 0

    # Isolate the data for the system.
    data_sub = submission.data[submission.data['System name'] == system_name]
    n_energy_evaluations = max(data_sub['N energy evaluations'])

    data_ref = yank_analysis.get_free_energies_from_energy_evaluations(
        n_energy_evaluations, system_name=system_name, mean_trajectory=False,
        start=calibration_cost)

    # Obtain the free energies for the submission.
    n_replicates = 5
    free_energy_sub = np.empty(shape=(n_replicates, 100))
    free_energy_ref = np.empty(shape=(n_replicates, 100))

    for data, free_energy in [
        (data_sub, free_energy_sub),
        (data_ref, free_energy_ref),
    ]:
        for i in range(n_replicates):
            system_id = system_name + '-' + str(i)
            system_id_data = data[data['System ID'] == system_id]
            free_energy[i] = system_id_data[DG_KEY].values

    # Discard the initial frames of WExplorer and GROMACS/EE that don't have predictions.
    from pkganalysis.efficiency import discard_initial_zeros
    free_energy_ref, free_energy_sub = discard_initial_zeros(free_energy_ref, free_energy_sub)

    # Determine the actual asymptotic free energy of YANK.
    asymptotic_free_energy_ref = yank_analysis.get_reference_free_energies()[system_name]
    return free_energy_ref, free_energy_sub, asymptotic_free_energy_ref

def compute_all_relative_efficiencies(
        free_energy_A, free_energy_B, ci, n_bootstrap_samples,
        asymptotic_free_energy_A=None, asymptotic_free_energy_B=None
):
    from pkganalysis.efficiency import EfficiencyAnalysis

    analysis = EfficiencyAnalysis(free_energy_A, free_energy_B,
                                  asymptotic_free_energy_A,
                                  asymptotic_free_energy_B)

    std_rel_eff = analysis.compute_std_relative_efficiency(
        confidence_interval=ci, n_bootstrap_samples=n_bootstrap_samples)
    abs_bias_rel_eff = analysis.compute_abs_bias_relative_efficiency(
        confidence_interval=ci, n_bootstrap_samples=n_bootstrap_samples)
    rmse_rel_eff = analysis.compute_rmse_relative_efficiency(
        confidence_interval=ci, n_bootstrap_samples=n_bootstrap_samples)

    if ci is None:
        rel_eff = [std_rel_eff, abs_bias_rel_eff, rmse_rel_eff]
        return rel_eff
    else:
        rel_eff = [std_rel_eff[0], abs_bias_rel_eff[0], rmse_rel_eff[0]]
        cis = [std_rel_eff[1], abs_bias_rel_eff[1], rmse_rel_eff[1]]
        return rel_eff, cis


def plot_relative_efficiencies(submissions, yank_analysis, ci=0.95, n_bootstrap_samples=1000,
                               same_plot=False, step_cumulative=2):
    sns.set_context('paper')

    statistic_names = ['std', 'absolute bias', 'RMSE']

    # Create output directory.
    figure_dir_path = os.path.join(SAMPLING_PAPER_DIR_PATH, 'SI_Figure2-efficiencies')
    os.makedirs(figure_dir_path, exist_ok=True)

    # Check if we need all the efficiencies in the same plot or not.
    if same_plot:
        fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(7.25, 8))
        # Keep track of data range by statistic.
        statistic_ranges = {name: [np.inf, 0] for name in statistic_names}
    # Keep track of n_energy_evaluations by column.
    max_n_energy_evaluations = [0 for _ in range(3)]

    for submission in submissions:
        if submission.paper_name in {'OpenMM/WExplore'}:
            continue
        # if submission.paper_name in {'AMBER/APR', 'GROMACS/CT-NS', 'GROMACS/CT-NS-long',
        #                              'GROMACS/EE', 'OpenMM/SOMD'}:
        #     continue
        # if submission.paper_name not in {'AMBER/APR'}:
        #     continue
        print(submission.paper_name)

        system_names = submission.data['System name'].unique()

        # Create figure.
        if not same_plot:
            # For GROMACS/EE, there are no submissions for CB8-G3.
            if 'GROMACS/EE' in submission.paper_name:
                 system_names = system_names[~(system_names == 'CB8-G3')]

            fig, axes = plt.subplots(nrows=3, ncols=len(system_names),
                                     figsize=(7.25, 8))
            statistic_ranges = {name: [np.inf, 0] for name in statistic_names}

        for col_idx, system_name in enumerate(system_names):
            # For GROMACS/EE, there are no submissions for CB8-G3.
            if 'GROMACS/EE' in submission.paper_name and system_name == 'CB8-G3':
                continue

            # Get input for EfficiencyAnalysis.
            free_energy_ref, free_energy_sub, asymptotic_free_energy_ref = get_relative_efficiency_input(
                submission, yank_analysis, system_name)

            # Get the relative efficiencies.
            rel_eff = compute_all_relative_efficiencies(
                free_energy_ref, free_energy_sub, ci, n_bootstrap_samples,
                asymptotic_free_energy_A=asymptotic_free_energy_ref
            )
            if ci is not None:
                rel_eff, cis = rel_eff  # Unpack confidence intervals.

            # Use the same asymptotic free energies to compute the absolute bias
            # relative efficiency as a function of the simulation length.
            asymptotic_free_energy_sub = free_energy_sub.mean(axis=0)[-1]

            # # Print relative efficiencies.
            # print(system_name, ci)
            # if ci is not None:
            #     for rel_eff, bounds in zip(rel_eff, cis):
            #         print('\t', rel_eff, bounds.tolist())
            # else:
            #     for rel_eff in rel_eff:
            #         print('\t', rel_eff)

            # Compute mean efficiencies as a function of the length of the simulation.
            n_costs = free_energy_ref.shape[1]
            n_rel_eff = int(n_costs / step_cumulative)
            relative_efficiencies = np.empty(shape=(3, n_rel_eff))
            low_bounds = np.empty(shape=(3, n_rel_eff))
            high_bounds = np.empty(shape=(3, n_rel_eff))

            for i, c in enumerate(range(step_cumulative-1, n_costs, step_cumulative)):
                c1 = c + 1

                rel_eff = compute_all_relative_efficiencies(
                    free_energy_ref[:,:c1], free_energy_sub[:,:c1],
                    ci, n_bootstrap_samples,
                    asymptotic_free_energy_A=asymptotic_free_energy_ref,
                    asymptotic_free_energy_B=asymptotic_free_energy_sub
                )
                if ci is not None:
                    rel_eff, cis = rel_eff  # Unpack confidence intervals.

                # Update CI lower and upper bound.
                relative_efficiencies[:,i] = rel_eff
                if ci is not None:
                    low_bounds[:,i]  = [x[0] for x in cis]
                    high_bounds[:,i]  = [x[1] for x in cis]

            # Plot.
            color = SUBMISSION_COLORS[submission.paper_name]

            # Get number of energy evaluations.
            mean_data = submission.mean_free_energies(system_name=system_name)
            # Check how many initial iteration have been discarded.
            discarded_iterations = 100 - n_costs
            n_energy_evaluations = mean_data['N energy evaluations'].values[
                                        discarded_iterations+1::step_cumulative] / 1e6
            for row_idx, rel_eff in enumerate(relative_efficiencies):
                ax = axes[row_idx][col_idx]
                ax.plot(n_energy_evaluations, rel_eff, color=color, label=submission.paper_name)

                # Plot back line at 0.
                ax.plot(n_energy_evaluations, [0 for _ in n_energy_evaluations], color='black', ls='--')

                # Update data range.
                statistic_range = statistic_ranges[statistic_names[row_idx]]
                # if ci is None:
                #     min_rel_eff = min(rel_eff)
                #     max_rel_eff = max(rel_eff)
                # else:
                #     min_rel_eff = min(*rel_eff, *low_bounds[row_idx])
                #     max_rel_eff = max(*rel_eff, *high_bounds[row_idx])
                statistic_range[0] = min(statistic_range[0], min(rel_eff))
                statistic_range[1] = max(statistic_range[1], max(rel_eff))

            # Update x-axis range.
            if same_plot:
                max_n_energy_evaluations[col_idx] = max(max_n_energy_evaluations[col_idx],
                                                        n_energy_evaluations[-1])
            else:
                for row_idx in range(len(statistic_names)):
                    axes[row_idx][col_idx].set_xlim((0, n_energy_evaluations[-1]))

            if ci is not None:
                # Plot confidence intervals.
                for row_idx, (low_bound_c, high_bound_c) in enumerate(zip(low_bounds, high_bounds)):
                    ax = axes[row_idx][col_idx]
                    ax.fill_between(n_energy_evaluations, low_bound_c, high_bound_c,
                                    alpha=0.35, color='gray')

            # We do this multiple times unnecessarily if same_plot is True, but the code is simpler.
            for col_idx, system_name in enumerate(system_names):
                axes[0][col_idx].set_title(system_name)
            for row_idx, statistic_name in enumerate(statistic_names):
                axes[row_idx][0].set_ylabel(statistic_name + ' rel eff')
                for col_idx in range(len(system_names)):
                    if same_plot:
                        extra_space = 0.1
                    else:
                        # Make space for confidence intervals.
                        extra_space = 1
                    ylimits = (statistic_ranges[statistic_name][0] - extra_space,
                               statistic_ranges[statistic_name][1] + extra_space)
                    axes[row_idx][col_idx].set_ylim(ylimits)
                    axes[row_idx][col_idx].tick_params(axis='y', which='major', pad=0.1)
            axes[-1][1].set_xlabel('Calculation length')

        # Set labels and axes limits.
        if not same_plot:
            fig.suptitle(submission.paper_name)

            output_file_base_name = '{}-{}'.format(submission.receipt_id, submission.file_name)
            output_file_base_path = os.path.join(figure_dir_path, output_file_base_name)
            plt.savefig(output_file_base_path + '.pdf')
            # plt.savefig(output_file_base_path + '.png', dpi=600)
            # plt.show()

    if same_plot:
        for row_idx in range(len(statistic_names)):
            for col_idx in range(len(system_names)):
                axes[row_idx][col_idx].set_xlim((0, max_n_energy_evaluations[col_idx]))
        axes[0][1].legend(loc='upper right', bbox_to_anchor=(2.0, 1.48),
                          fancybox=True, ncol=3)

        output_file_base_path = os.path.join(figure_dir_path, 'relative-efficiencies')
        plt.savefig(output_file_base_path + '.pdf')
        # plt.savefig(output_file_base_path + '.png', dpi=600)
        # plt.show()


def plot_absolute_efficiencies(submissions, yank_analysis, ci=0.95, n_bootstrap_samples=1000):
    sns.set_context('paper')

    # Keep track of data range by statistic.
    statistic_names = ['std', 'absolute bias', 'RMSE']

    # Keep track of maximum number of energy evaluations
    # to determine plotting range for YANK.
    system_names = ['CB8-G3', 'OA-G3', 'OA-G6']
    max_n_energy_eval = {name: 0 for name in system_names}

    # Create figure.
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(7.25, 8))

    for submission in submissions + [yank_analysis]:
        if 'WExplore' in submission.paper_name:
            continue
        print(submission.paper_name)

        # Obtain std, bias, and RMSE of the 5 trajectories.
        # If this is a YANK analysis, we get it later specifically for the system.
        if not isinstance(submission, YankSamplingAnalysis):
            mean_free_energies = submission.mean_free_energies()

        for col_idx, system_name in enumerate(system_names):
            # GROMACS/EE doesn't have submissions for CB8-G3.
            if 'GROMACS/EE' in submission.paper_name and system_name == 'CB8-G3':
                continue

            # Select the submission data for only this host-guest system.
            if isinstance(submission, YankSamplingAnalysis):
                line_style = '--'
                mean_data = submission.get_free_energies_from_energy_evaluations(
                    max_n_energy_eval[system_name], system_name=system_name, mean_trajectory=True)
            else:
                line_style = '-'
                mean_data = mean_free_energies[mean_free_energies['System name'] == system_name]

            # Update maximum number of energy evaluations.
            n_energy_evaluations = mean_data['N energy evaluations'].values
            max_n_energy_eval[system_name] = max(max_n_energy_eval[system_name], n_energy_evaluations[-1])

            # Discard initial computational costs for which there's no data.
            first_nonzero_idx = np.nonzero(mean_data[DG_KEY])[0][0]
            n_energy_evaluations = n_energy_evaluations[first_nonzero_idx:]

            # Compute cumulative total std, abs_bias, and RMSE.
            scale_energy_evaluations = 1e6
            norm_factor = (n_energy_evaluations - n_energy_evaluations[0])[1:] / scale_energy_evaluations
            avg_std = sp.integrate.cumtrapz(mean_data['std'].values[first_nonzero_idx:]) / norm_factor
            avg_abs_bias = sp.integrate.cumtrapz(np.abs(mean_data['bias'].values[first_nonzero_idx:])) / norm_factor
            avg_rmse = sp.integrate.cumtrapz(mean_data['RMSE'].values[first_nonzero_idx:]) / norm_factor

            # Plot total statistics as a function of the energy evaluations.
            # Discard first energy evaluation as cumtrapz doesn't return a result for it.
            for row_idx, avg_stats in enumerate([avg_std, avg_abs_bias, avg_rmse]):
                ax = axes[row_idx, col_idx]
                color = SUBMISSION_COLORS[submission.paper_name]
                ax.plot(n_energy_evaluations[1:] / scale_energy_evaluations, avg_stats,
                        color=color, label=submission.paper_name, ls=line_style)

                # Set x axis.
                ax.set_xlim((0, n_energy_evaluations[-1] / scale_energy_evaluations))

    # Set labels and axes limits.
    y_limits = {
        'std': (0, 0.4),
        'absolute bias': (0, 0.3),
        'RMSE': (0, 0.4)
    }
    for col_idx, system_name in enumerate(system_names):
        axes[0][col_idx].set_title(system_name)
        # Set y limits (shared for each row).
        for row_idx, statistic_name in enumerate(statistic_names):
            axes[row_idx][col_idx].set_ylim(y_limits[statistic_name])
            axes[row_idx][col_idx].tick_params(axis='y', which='major', pad=0.1)

    # # Remove shared ticks.
    # for row_idx in range(len(statistic_names)):
    #     for col_idx in range(len(system_names)):
    #         if col_idx > 0:
    #             axes[row_idx][col_idx].set_yticklabels([])
    #         if row_idx < len(statistic_names)-1:
    #             axes[row_idx][col_idx].set_xticklabels([])

    for row_idx, statistic_name in enumerate(statistic_names):
        axes[row_idx][0].set_ylabel('mean ' + statistic_name + ' [kcal/mol]')
    axes[-1][1].set_xlabel('N energy evaluations [M]')

    axes[0][1].legend(loc='upper right', bbox_to_anchor=(2.0, 1.48),
                      fancybox=True, ncol=3)

    figure_dir_path = os.path.join(SAMPLING_PAPER_DIR_PATH, 'SI_Figure2-efficiencies')
    os.makedirs(figure_dir_path, exist_ok=True)
    output_file_base_path = os.path.join(figure_dir_path, 'absolute-efficiencies')
    plt.savefig(output_file_base_path + '.pdf')
    # plt.savefig(output_file_base_path + '.png', dpi=600)
    # plt.show()


def print_relative_efficiency_table(
        submissions, yank_analysis, ci=0.95,
        n_bootstrap_samples=100,
        print_bias_corrected=False
):
    """Create a table with standard deviation, absolute bias, and RMSE relative efficiency."""
    methods = []

    # Initialize the table to be converted into a Pandas dataframe.
    system_names = ['CB8-G3', 'OA-G3', 'OA-G6']
    statistic_names = [r'$e_{\text{std}}$', r'$e_{|\text{bias}|}$', r'$e_{\text{RMSD}}$']
    column_names  = ['\\makecell{$\Delta$ G \\\\ $[$kcal/mol$]$}', '\\makecell{n eval \\\\ $[$M$]$}'] + statistic_names
    # Add columns.
    efficiency_table = collections.OrderedDict()
    for system_name, column_name in itertools.product(system_names, column_names):
        efficiency_table[(system_name, column_name)] = []

    for submission in submissions:
        # Collect method's names in the given order.
        methods.append(submission.paper_name)

        mean_free_energies = submission.mean_free_energies()
        for system_name in system_names:
            # CB8-G3 calculations for GROMACS/EE did not converge yet, and the
            # long protocol in CS-NS calculations have been run only on CB8-G3.
            if ((submission.name == 'Expanded-ensemble/MBAR' and system_name == 'CB8-G3') or
                    (submission.paper_name == 'GROMACS/CT-NS-long' and system_name != 'CB8-G3')):
                relative_efficiencies, relative_efficiencies_corrected = np.full((2, 3), fill_value=np.nan)
                dg = ''
                n_force_eval = ''
            else:
                # Get input for EfficiencyAnalysis.
                free_energy_ref, free_energy_sub, asymptotic_free_energy_ref = get_relative_efficiency_input(
                    submission, yank_analysis, system_name)

                # Get the relative efficiencies.
                relative_efficiencies, cis = compute_all_relative_efficiencies(
                    free_energy_ref, free_energy_sub, ci, n_bootstrap_samples,
                    asymptotic_free_energy_A=asymptotic_free_energy_ref
                )

                # Recompute relative efficiencies assuming that YANK converged.
                if print_bias_corrected:
                    relative_efficiencies_corrected, cis_corrected = compute_all_relative_efficiencies(
                        free_energy_ref, free_energy_sub, ci, n_bootstrap_samples)

                # Select the data for only this host-guest system.
                mean_data_sub = mean_free_energies[mean_free_energies['System name'] == system_name]

                # Get the final free energy and number of energy/force evaluations.
                dg = mean_data_sub[DG_KEY].values[-1]
                dg_CI = mean_data_sub['$\Delta$G CI'].values[-1]  # Confidence interval.
                dg, dg_CI = reduce_to_first_significant_digit(dg, dg_CI)
                n_force_eval = mean_data_sub['N energy evaluations'].values[-1]
                # Convert to string format.
                dg = '{} $\\pm$ {}'.format(dg, dg_CI)
                n_force_eval = str(int(round(n_force_eval / 1e6)))

            # Add free energy and cost entries.
            efficiency_table[(system_name, column_names[0])].append(dg)
            efficiency_table[(system_name, column_names[1])].append(n_force_eval)

            # Add efficiency entries for the table.
            for statistic_idx, statistic_name in enumerate(statistic_names):
                # Gather the format arguments.
                rel_effs = [relative_efficiencies[statistic_idx], cis[statistic_idx][0], cis[statistic_idx][1]]
                if print_bias_corrected:
                    rel_effs.append(relative_efficiencies_corrected[statistic_idx])
                    # Comment this if we don't want to print CIs for the corrected estimate.
                    rel_effs.extend([cis_corrected[statistic_idx][0], cis_corrected[statistic_idx][1]])

                # Print significant digits.
                efficiencies_format = []
                for e_idx in range(0, len(rel_effs), 3):
                    rel_eff, low_bound, high_bound = rel_effs[e_idx:e_idx+3]
                    if high_bound - rel_eff < 0.1 or rel_eff - low_bound < 0.1:
                        fmt = '{:2.2f}'
                    else:
                        fmt = '{:2.1f}'
                    # Print lower and higher bound as sub and superscripts of the estimate.
                    efficiencies_format.append(fmt + '$_{{\raisem{{2pt}}{{' + fmt + '}}}}^{{\mathstrut ' + fmt + '}}$')

                if np.isnan(rel_effs[0]):
                    data_entry = ''
                # Standard deviation efficiency is not affected by the bias.
                elif print_bias_corrected and ('std' not in statistic_name):
                    data_entry = efficiencies_format[0] + ' (' + efficiencies_format[1] + ')'
                    data_entry = data_entry.format(*rel_effs)
                else:
                    data_entry = efficiencies_format[0].format(*rel_effs[:3])

                # Remove the minus sign from "-0".
                data_entry = data_entry.replace('-0.0', '0.0')
                data_entry = data_entry.replace('-0.00', '0.00')
                efficiency_table[(system_name, statistic_name)].append(data_entry)

    # Add row for reference calculation.
    methods.append(YANK_METHOD_PAPER_NAME)

    # Add free energy and cost entries.
    for system_name in system_names:
        yank_mean_data = yank_analysis.get_free_energies_from_iteration(
            YANK_N_ITERATIONS, system_name=system_name, mean_trajectory=True)
        dg = yank_mean_data[DG_KEY].values[-1]
        dg_CI = yank_mean_data['$\Delta$G CI'].values[-1]  # Confidence interval.
        dg, dg_CI = reduce_to_first_significant_digit(dg, dg_CI)
        n_force_eval = yank_mean_data['N energy evaluations'].values[-1]
        n_force_eval = str(int(round(n_force_eval / 1e6)))
        efficiency_table[(system_name, column_names[0])].append('{} $\\pm$ {}'.format(dg, dg_CI))
        efficiency_table[(system_name, column_names[1])].append(n_force_eval)

    # All efficiencies are relative to YANK so they're all 1.
    for system_name, statistic_name in itertools.product(system_names, statistic_names):
        efficiency_table[(system_name, statistic_name)].append('0.0')

    # Convert to Pandas Dataframe.
    efficiency_table = pd.DataFrame(efficiency_table)
    # Set the method's names as index column.
    efficiency_table = efficiency_table.assign(Method=methods)
    efficiency_table.set_index(keys='Method', inplace=True)

    # Print table.
    column_format = 'lccccc|ccccc|ccccc'
    efficiency_table_latex = efficiency_table.to_latex(column_format=column_format, multicolumn_format='c',
                                                       escape=False)

    # Make header and reference method bold.
    textbf = lambda s: '\\textbf{' + s + '}'
    efficiency_table_latex = efficiency_table_latex.replace(YANK_METHOD_PAPER_NAME, textbf(YANK_METHOD_PAPER_NAME))
    efficiency_table_latex = efficiency_table_latex.replace('Method', textbf('Method'))
    for system_name in system_names:
        efficiency_table_latex = efficiency_table_latex.replace(system_name, textbf(system_name))
    for column_name in column_names:
        efficiency_table_latex = efficiency_table_latex.replace(column_name, textbf(column_name))
    print(efficiency_table_latex)


# =============================================================================
# SUPPORTING INFORMATION FIGURES
# =============================================================================

def plot_single_trajectories_figures(axes, system_data, system_mean_data,
                                     reference_system_mean_data=None,
                                     plot_errors=True, plot_methods_uncertainties=True):
    """Plot individual free energy trajectories and standard deviations for a single method and system."""
    system_name = system_data['System name'].unique()[0]
    palette_mean = sns.color_palette('pastel')
    submission_mean_color = 'black'
    reference_mean_color = palette_mean[9]

    # Plot the method uncertainties of the single replicate trajectories.
    # First scale the number of energy evaluations.
    system_data.loc[:,'N energy evaluations'] /= N_ENERGY_EVALUATIONS_SCALE

    # Plot the 5 replicates individual trajectories.
    # First remove the initial predictions that are 0.0 (i.e. there is no estimate).
    ax = axes[0]
    system_data = system_data[system_data[DG_KEY] != 0.0]
    sns.lineplot(data=system_data, x='N energy evaluations', y=DG_KEY,
                 hue='System ID', palette='bright', ax=ax, alpha=0.6)

    # Plot the submission mean trajectory with CI.
    plot_mean_free_energy(system_mean_data, x='N energy evaluations',  ax=ax,
                          color_mean=submission_mean_color, plot_ci=False,
                          color_ci=submission_mean_color, label='Best estimate',
                          scale_n_energy_evaluations=True)

    # Plot YANK mean trajectory with CI.
    if reference_system_mean_data is not None:
        plot_mean_free_energy(reference_system_mean_data, x='N energy evaluations', ax=ax,
                              color_mean=reference_mean_color, plot_ci=False,
                              color_ci=reference_mean_color, label='Reference estimate',
                              scale_n_energy_evaluations=True)

    ax.set_title(system_name)
    # Add the y-label only on the leftmost Axis.
    if system_name != 'CB8-G3':
        ax.set_ylabel('')
    # Remove the legend for now, which will be added at the end after tighting up the plot.
    ax.get_legend().remove()

    # Create a bias axis.
    if reference_system_mean_data is not None:
        ref_free_energy = reference_free_energies.loc[system_name, DG_KEY]
        with sns.axes_style('white'):
            ax2 = ax.twinx()
            # Plot a vertical line to make the scale.
            vertical_line = np.linspace(*ax.get_ylim()) - ref_free_energy
            ax2.plot([50] * len(vertical_line), vertical_line, alpha=0.0001)
            ax2.grid(alpha=0.5, linestyle='dashed', zorder=0)
            # We add the bias y-label only on the rightmost Axis.
            if system_name == 'OA-G6':
                ax2.set_ylabel('Bias to reference [kcal/mol]')
            # Set the 0 of the twin axis to the YANK reference free energy.
            align_yaxis(ax, ref_free_energy, ax2, 0.0)

    if plot_errors:
        # The x-axis is shared between the 2 rows so we can plot the ticks only in the bottom one.
        ax.xaxis.set_ticklabels([])
        ax.set_xlabel('')

        ax = axes[1]

        # WExplore uses the mean of the 5 replicates to estimate the
        # uncertainty so it doesn't add information.
        if plot_methods_uncertainties:
            sns.lineplot(data=system_data, x='N energy evaluations', y=DDG_KEY,
                         hue='System ID', palette='bright', ax=ax, alpha=0.6)

            # The legend is added later at the top.
            ax.get_legend().remove()

        # Plot the standard deviation of the free energy trajectories.
        # submission_std = system_mean_data['std']
        submission_std = system_mean_data['unbiased_std']
        # cost = system_mean_data['Simulation percentage'].values
        cost = system_mean_data['N energy evaluations'].values / N_ENERGY_EVALUATIONS_SCALE
        ax.plot(cost, submission_std, color=submission_mean_color)

        # Plot confidence interval around standard deviation.
        submission_std_low_ci = system_mean_data['unbiased_std_low_CI'].values
        submission_std_up_ci = system_mean_data['unbiased_std_up_CI'].values
        ax.fill_between(cost, submission_std_low_ci, submission_std_up_ci, alpha=0.35, color='gray')

        if reference_system_mean_data is not None:
            # reference_std = reference_system_mean_data['std']
            reference_std = reference_system_mean_data['unbiased_std']
            ax.plot(cost, reference_std, color=reference_mean_color)

        # Only the central plot shows the x-label.
        ax.set_xlabel('')
        # Add the y-label only on the leftmost Axis.
        if system_name != 'CB8-G3':
            ax.set_ylabel('')
        else:
            ax.set_ylabel('std($\Delta$G) [kcal/mol]')

    # Set x limits.
    for ax in axes:
        ax.set_xlim((0, max(system_data['N energy evaluations'])))


def plot_all_single_trajectories_figures(submissions, yank_analysis, plot_errors=True, output_path_dir=None):
    """Individual plots for each method with the 5 individual free energy and uncertainty trajectories."""
    sns.set_context('paper')

    if output_path_dir is None:
        output_path_dir = os.path.join(SAMPLING_PAPER_DIR_PATH, 'SI_Figure1-individual-trajectories/')
    os.makedirs(output_path_dir, exist_ok=True)

    # -------------------- #
    # Plot submission data #
    # -------------------- #

    # Remove nonequilibrium-switching calculations with single-direction estimators.
    submissions = [s for s in submissions if ('Jarz' not in s.paper_name and 'Gauss' not in s.paper_name)]

    for submission in submissions + [yank_analysis]:
        # CB8-G3 calculations for GROMACS/EE did not converge yet.
        if submission.name == 'Expanded-ensemble/MBAR':
            submission.data = submission.data[submission.data['System name'] != 'CB8-G3']
        # WExplore uses the mean of the 5 replicates to estimate the
        # uncertainty so it doesn't add information.
        if 'WExplore' in submission.paper_name:
            plot_methods_uncertainties = False
        else:
            plot_methods_uncertainties = True

        if not isinstance(submission, YankSamplingAnalysis):
            mean_free_energies = submission.mean_free_energies()
            unique_system_names = submission.data['System name'].unique()
        else:
            unique_system_names = sorted(submission.system_names)

        # Create a figure with 3 axes (one for each system).
        n_systems = len(unique_system_names)
        if plot_errors:
            # The second row will plot the errors.
            fig, axes = plt.subplots(nrows=2, ncols=n_systems, figsize=(7.25, 4.8))
            trajectory_axes = axes[0]
        else:
            fig, axes = plt.subplots(nrows=1, ncols=n_systems, figsize=(7.25, 2.4))
            trajectory_axes = axes

        # Set figure title.
        fig.suptitle(submission.paper_name)

        # Determine range of data across systems.
        min_DG = np.inf
        max_DG = -np.inf
        min_dDG = np.inf
        max_dDG = -np.inf

        # for system_name in unique_system_names:
        for ax_idx, system_name in enumerate(unique_system_names):

            if isinstance(submission, YankSamplingAnalysis):
                data = submission.get_free_energies_from_iteration(final_iteration=YANK_N_ITERATIONS,
                                                                   system_name=system_name)
                mean_data = submission.get_free_energies_from_iteration(final_iteration=YANK_N_ITERATIONS,
                                                                        system_name=system_name,
                                                                        mean_trajectory=True)
            else:
                # Select the data for only this host-guest system.
                data = submission.data[submission.data['System name'] == system_name]
                mean_data = mean_free_energies[mean_free_energies['System name'] == system_name]

            plot_single_trajectories_figures(axes[:,ax_idx], data, mean_data, plot_errors=plot_errors,
                                             reference_system_mean_data=None,
                                             plot_methods_uncertainties=plot_methods_uncertainties)

            # Collect max and min data to determine axes range.
            min_DG = min(min_DG, min(data[DG_KEY]), min(mean_data[DG_KEY]))
            max_DG = max(max_DG, max(data[DG_KEY]), max(mean_data[DG_KEY]))
            min_dDG = min(min_dDG, min(data[DDG_KEY]), min(mean_data['std']))
            max_dDG = max(max_dDG, max(data[DDG_KEY]), max(mean_data['std']))

        # Set limits.
        for i in range(len(unique_system_names)):
            axes[0][i].set_ylim((min_DG, max_DG))
            axes[1][i].set_ylim((min_dDG, max_dDG))
            # Keep ticks only in external plots.
            axes[0][i].set_xticklabels([])
        for i in range(1, len(unique_system_names)):
            axes[0][i].set_yticklabels([])
            axes[1][i].set_yticklabels([])

        # The x-label is shown only in the central plot.
        axes[-1][1].set_xlabel('N energy evaluations  [10$^6$]')

        plt.tight_layout(pad=0.2, rect=[0.0, 0.0, 1.0, 0.85])

        # Create legend.
        # The first handle/label is the legend title "System ID" so we get rid of it.
        handles, labels = trajectory_axes[0].get_legend_handles_labels()
        labels = ['replicate ' + str(i) for i in range(5)] + labels[6:]
        bbox_to_anchor = (-0.1, 1.35)
        trajectory_axes[0].legend(handles=handles[1:], labels=labels, loc='upper left',
                                  bbox_to_anchor=bbox_to_anchor, ncol=6, fancybox=True,
                                  labelspacing=0.8, handletextpad=0.5, columnspacing=1.2)
        # Save figure.
        output_file_name = '{}-{}'.format(submission.receipt_id, submission.file_name)
        plt.savefig(os.path.join(output_path_dir, output_file_name + '.pdf'))
        # plt.savefig(os.path.join(output_path_dir, output_file_name + '.png'), dpi=300)
        # plt.show()


def plot_hrex_stat_ineff_trajectories():
    """Individual plots for HREX with the 5 individual free energy and uncertainty trajectories
    as a function of the statistical inefficiency."""
    sns.set_context('paper')

    # Limits of y-axis (free energies, uncertainties) by system.
    y_limits = {
        'CB8-G3': [(-14, -10), (0, 2)],
        'OA-G3': [(-9, -5), (0, 1.5)],
        'OA-G6': [(-9, -5), (0, 1.5)],
    }

    # Create output dir.
    output_path_dir = os.path.join(SAMPLING_PAPER_DIR_PATH, 'SI_Figure-statistical-inefficiency')
    os.makedirs(output_path_dir, exist_ok=True)

    # Read the data, which is organized by statistical inefficiency.
    # We'll then plot by system.
    yank_analysis_by_statineff = collections.OrderedDict()
    for stat_ineff in ['2', '5', '10', '20', '50', '100']:
        data_dir_path = os.path.join('YankAnalysis', 'CorrelationAnalysis', 'statineff-{}'.format(stat_ineff))
        yank_analysis = YankSamplingAnalysis(data_dir_path)
        yank_analysis_by_statineff[stat_ineff] = yank_analysis

    # Plot by system.
    for system_name in ['CB8-G3', 'OA-G3', 'OA-G6']:
        fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(7.25, 9.8))

        # Set figure title.
        fig.suptitle('HREX uncertainty predictions as a function of\n'
                     'statistical inefficiency for {}'.format(system_name))

        # for system_name in unique_system_names:
        for stat_ineff_idx, stat_ineff in enumerate(yank_analysis_by_statineff):
            yank_analysis = yank_analysis_by_statineff[stat_ineff]
            data = yank_analysis.get_free_energies_from_iteration(final_iteration=YANK_N_ITERATIONS,
                                                               system_name=system_name)
            mean_data = yank_analysis.get_free_energies_from_iteration(final_iteration=YANK_N_ITERATIONS,
                                                                    system_name=system_name,
                                                                    mean_trajectory=True)

            # Plot on the correct axis.
            DG_row = 2*int(stat_ineff_idx / 3)
            col = stat_ineff_idx % 3
            stat_ineff_axes = axes[DG_row:DG_row+2, col]
            plot_single_trajectories_figures(stat_ineff_axes, data, mean_data, plot_errors=True,
                                             reference_system_mean_data=None,
                                             plot_methods_uncertainties=True)

            # Set titles and limits.
            title = 'Statistical inefficiency: {} ps'.format(stat_ineff)
            if DG_row > 0:
                title = '\n' + title
            stat_ineff_axes[0].set_title(title, fontweight='bold')
            stat_ineff_axes[0].set_ylim(y_limits[system_name][0])
            stat_ineff_axes[1].set_ylim(y_limits[system_name][1])

        # Keep ticks only in external plots.
        for row_idx in range(axes.shape[0]):
            for col_idx in range(axes.shape[1]):
                if row_idx != len(axes[0]) - 1:
                    axes[row_idx][col_idx].set_xticklabels([])
                if col_idx != 0:
                    axes[row_idx][col_idx].set_ylabel('')
                    axes[row_idx][col_idx].set_yticklabels([])
        # Set x label.
        axes[-1][1].set_xlabel('N energy evaluations  [10$^6$]')

        plt.tight_layout(pad=0.0, rect=[0.0, 0.0, 1.0, 0.88])

        # Create legend.
        # The first handle/label is the legend title "System ID" so we get rid of it.
        handles, labels = axes[0][0].get_legend_handles_labels()
        labels = ['replicate ' + str(i) for i in range(5)] + labels[6:]
        bbox_to_anchor = (-0.1, 1.35)
        axes[0][0].legend(handles=handles[1:], labels=labels, loc='upper left',
                          bbox_to_anchor=bbox_to_anchor, ncol=6, fancybox=True,
                          labelspacing=0.8, handletextpad=0.5, columnspacing=1.2)

        # Save figure.
        output_file_name = 'statineff-{}'.format(system_name)
        # plt.savefig(os.path.join(output_path_dir, output_file_name + '.pdf'))
        plt.savefig(os.path.join(output_path_dir, output_file_name + '.png'), dpi=300)
        # plt.show()


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    sns.set_style('whitegrid')
    sns.set_context('paper')

    # Read reference values.
    yank_analysis = YankSamplingAnalysis(YANK_ANALYSIS_DIR_PATH)

    # Obtain free energies and final reference values.
    mean_reference_free_energies = yank_analysis.get_free_energies_from_iteration(YANK_N_ITERATIONS, mean_trajectory=True)
    reference_free_energies = mean_reference_free_energies[mean_reference_free_energies['Simulation percentage'] == 100]
    reference_free_energies.set_index('System name', inplace=True)

    # Compute efficiency of reference.
    reference_efficiencies = {}
    for system_name in mean_reference_free_energies['System name'].unique():
        mean_data = mean_reference_free_energies[mean_reference_free_energies ['System name'] == system_name]
        reference_efficiencies[system_name], n_discarded = fit_efficiency(mean_data)

    # Import user map.
    with open('../SubmissionsDoNotUpload/SAMPL6_user_map.csv', 'r') as f:
        user_map = pd.read_csv(f)

    # Load submissions data. We do OA and TEMOA together.
    all_submissions = load_submissions(SamplingSubmission, SAMPLING_SUBMISSIONS_DIR_PATH, user_map)
    # Remove AMBER/TI.
    all_submissions = [s for s in all_submissions if s.name not in ['Langevin/Virtual Bond/TI']]

    # Create an extra submission for GROMACS/EE where the full cost of equilibration has been taken into account.
    gromacs_ee_submission = copy.deepcopy([s for s in all_submissions if s.paper_name == 'GROMACS/EE'][0])
    gromacs_ee_submission.paper_name = 'GROMACS/EE-fullequil'
    gromacs_ee_submission.file_name = 'EENVT-fullequil'
    data = gromacs_ee_submission.data  # Shortcut.
    mean_free_energies = gromacs_ee_submission.mean_free_energies()
    for system_name in ['OA-G3', 'OA-G6']:
        mean_data = mean_free_energies[mean_free_energies['System name'] == system_name]
        first_nonzero_idx = np.nonzero(mean_data[DG_KEY].values)[0][0]
        full_equilibration_cost = mean_data['N energy evaluations'].values[first_nonzero_idx] * 4
        for i in data[data['System name'] == system_name].index:
            data.at[i, 'N energy evaluations'] += full_equilibration_cost

    all_submissions.append(gromacs_ee_submission)

    # Sort the submissions to have all pot and tables in the same order.
    all_submissions = sorted(all_submissions, key=lambda s: s.paper_name)

    # Separate the main submissions from the data about nonequilibrium estimators.
    main_submissions = [s for s in all_submissions if not ('Jarz' in s.paper_name or 'Gauss' in s.paper_name)]
    noneq_submissions = [s for s in all_submissions if 'NS' in s.paper_name]

    # # Export YANK analysis and submissions to CSV/JSON tables.
    # yank_analysis.export(os.path.join(SAMPLING_DATA_DIR_PATH, 'reference_free_energies'))
    # for s in main_submissions:
    #     file_base_path = os.path.join(SAMPLING_DATA_DIR_PATH, s.receipt_id + '-reference')
    #     yank_analysis.export_by_submission(file_base_path, s)
    # export_submissions(all_submissions, reference_free_energies)

    # Create figure with free energy, standard deviation, and bias as a function of computational cost.
    # plot_all_entries_trajectory(main_submissions, yank_analysis, zoomed=False)
    # plot_all_entries_trajectory(main_submissions, yank_analysis, zoomed=True)

    # Create results and efficiency table.
    plot_relative_efficiencies(main_submissions, yank_analysis)
    # plot_relative_efficiencies(main_submissions, yank_analysis, ci=None, same_plot=True)
    # plot_absolute_efficiencies(main_submissions, yank_analysis)
    # print_relative_efficiency_table(main_submissions, yank_analysis, print_bias_corrected=False)

    # Plot nonequilibrium-switching single-direction estimator.
    # plot_all_nonequilibrium_switching(noneq_submissions)

    # Plot sensitivity analysis figure.
    # plot_restraint_and_barostat_analysis()

    # Plot figure for HREX bias analysis.
    # plot_yank_bias()
    # plot_equilibration_methods()

    # Plot individual trajectories.
    # plot_all_single_trajectories_figures(all_submissions, yank_analysis)

    # Plot statistical inefficiency analysis.
    # plot_hrex_stat_ineff_trajectories()
