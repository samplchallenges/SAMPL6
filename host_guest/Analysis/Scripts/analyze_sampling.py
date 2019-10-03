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
from pkganalysis.stats import mean_confidence_interval
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
    'AMBER/APR': 'dodgerblue',#KELLY_COLORS[11],
    'OpenMM/REVO': 'gold', #KELLY_COLORS[7],
    'OpenMM/SOMD': KELLY_COLORS[4],
    'GROMACS/EE': 'darkviolet', #KELLY_COLORS[3],
    'GROMACS/EE-fullequil': KELLY_COLORS[10],
    YANK_METHOD_PAPER_NAME: '#4ECC41', #'limegreen', #KELLY_COLORS[9],
    'GROMACS/NS-DS/SB-long': KELLY_COLORS[6],
    'GROMACS/NS-DS/SB': KELLY_COLORS[1],
    'GROMACS/NS-Jarz-F': TAB10_COLORS[0],
    'GROMACS/NS-Jarz-R': TAB10_COLORS[1],
    'GROMACS/NS-Gauss-F': TAB10_COLORS[2],
    'GROMACS/NS-Gauss-R': TAB10_COLORS[4],
    'NAMD/BAR': 'saddlebrown'
}

SUBMISSION_LINE_STYLES = {
    'AMBER/APR': '--',
    'OpenMM/REVO': '-',
    'OpenMM/SOMD': '-',
    'GROMACS/EE': '-',
    'GROMACS/EE-fullequil': '-',
    YANK_METHOD_PAPER_NAME: '-',
    'GROMACS/NS-DS/SB-long': '-',
    'GROMACS/NS-DS/SB': '-',
    'GROMACS/NS-Jarz-F': '-',
    'GROMACS/NS-Jarz-R': '-',
    'GROMACS/NS-Gauss-F': '-',
    'GROMACS/NS-Gauss-R': '-',
    'NAMD/BAR': '--',
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


def plot_mean_data(mean_data, axes, color=None, ls=None, label=None, x='N energy evaluations',
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
                          color_mean=color, color_ci=color, ls=ls, zorder=zorder,
                          start=first_nonzero_idx, label=label, plot_ci=plot_ci)

    # Plot standard deviation of the trajectories.
    if plot_std:
        axes[1].plot(mean_data[x].values[first_nonzero_idx:] / scale,
                     mean_data['std'].values[first_nonzero_idx:], color=color, alpha=0.8,
                     ls=ls, zorder=zorder, label=label)
    if plot_bias:
        axes[2].plot(mean_data[x].values[first_nonzero_idx:] / scale,
                     mean_data['bias'].values[first_nonzero_idx:], color=color, alpha=0.8,
                     ls=ls, zorder=zorder, label=label)


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
# FIGURE 1 - SAMPLING CHALLENGE OVERVIEW
# =============================================================================

def plot_example_bias_variance(yank_analysis, type='mixed', cost='generic',
                               max_n_eval_percentage=1.0,
                               mixed_proportion=0.5,
                               model_free_energy=None,
                               plot_experimental_value=False):
    """Free energy trajectories used to visualize bias and variance on the plots.

    This is used to illustrate how bias and uncertainty are intended in the paper.

    Parameters
    ----------
    type : str, optional
        Can be 'single' (plot only CB8-G3-1), 'all' (plot all system IDs of CB8-G3),
        'mean' (plot mean trajectory and uncertainties), and 'mixed (first part is
        all system IDs and second part is mean trajectory and uncertainties).
    cost : str, optional
        Can be 'generic' (no label on x-axis) or 'neval' (x-axis in number of
        energy evaluations).
    mixed_proportion : float, optional
        The proportion of all System IDs and mean trajectories in mixed-type plots.
    """
    # sns.set_context('paper', font_scale=1.6)
    sns.set_style('white')
    sns.set_context('paper', font_scale=1.0)

    # Load the data
    n_iterations = 40000
    cb8_data = yank_analysis.get_free_energies_from_iteration(n_iterations, system_name='CB8-G3', mean_trajectory=False)
    cb8_data_mean = yank_analysis.get_free_energies_from_iteration(n_iterations, system_name='CB8-G3', mean_trajectory=True)
    max_n_eval = max(cb8_data_mean['N energy evaluations'])
    max_n_eval_scaled = int(max_n_eval / N_ENERGY_EVALUATIONS_SCALE)
    max_displayed_n_eval = next(x for x in cb8_data_mean['N energy evaluations'] if x >= max_n_eval * max_n_eval_percentage)
    max_displayed_n_eval_scaled = int(max_displayed_n_eval / N_ENERGY_EVALUATIONS_SCALE)

    # Determine the asymptotic free energy if not given.
    if model_free_energy is None:
        model_free_energy = cb8_data_mean[DG_KEY].values[-1]

    # Scale the number of energy evaluations.
    cb8_data.loc[:,'N energy evaluations'] /= N_ENERGY_EVALUATIONS_SCALE

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(2.5, 1.8))

    if type == 'single':
        # Plot only CB8-G3-1.
        cb8_data_1 = cb8_data[cb8_data['System ID'] == 'CB8-G3-1']
        sns.lineplot(data=cb8_data_1, x='N energy evaluations', y=DG_KEY,
                     hue='System ID', palette='bright', ax=ax, alpha=0.6)
    elif type == 'all':
        # Plot the 5 replicates individual trajectories.
        sns.lineplot(data=cb8_data, x='N energy evaluations', y=DG_KEY,
                     hue='System ID', palette='bright', ax=ax, alpha=0.6)
    elif type == 'mean':
        # Plot the submission mean trajectory with CI.
        plot_mean_free_energy(cb8_data_mean, x='N energy evaluations',  ax=ax,
                              color_mean='black', plot_ci=True,
                              color_ci='black',
                              scale_n_energy_evaluations=True)
    elif type == 'mixed':
        # Plot all System IDs for the first half and mean/uncertainty in second half.
        half_n_eval = max_displayed_n_eval_scaled * mixed_proportion
        cb8_data_first_half = cb8_data[cb8_data['N energy evaluations'] <= half_n_eval + max_n_eval_scaled / 100]
        sns.lineplot(data=cb8_data_first_half, x='N energy evaluations', y=DG_KEY,
                     hue='System ID', palette='bright', ax=ax, alpha=0.6)

        cb8_data_second_half = cb8_data_mean[cb8_data_mean['N energy evaluations'] >= half_n_eval * N_ENERGY_EVALUATIONS_SCALE]
        plot_mean_free_energy(cb8_data_second_half, x='N energy evaluations',  ax=ax,
                              color_mean='black', plot_ci=True,
                              color_ci=(0.3, 0.3, 0.3), scale_n_energy_evaluations=True,
                              ls='--')

    try:
        ax.get_legend().remove()
    except AttributeError:
        pass

    # Set limits
    x_lim = (0, max_displayed_n_eval_scaled)
    ax.set_xlim(x_lim)
    y_lim = (-12.5, -10.5)
    ax.set_ylim(y_lim)

    # Plot model and experiment indication. Both values are not real data, just an example.
    model_free_energy = -10.75
    final_prediction = cb8_data_mean[cb8_data_mean['N energy evaluations'] == max_displayed_n_eval][DG_KEY].values[0]
    ax.plot(x_lim, [model_free_energy]*2, color='gray', ls='--')
    ax.text(x_lim[-1]+(max_n_eval_scaled*max_n_eval_percentage)/100, model_free_energy, r'$\Delta$G$_{\theta}$')
    ax.text(x_lim[-1]+(max_n_eval_scaled*max_n_eval_percentage)/100, final_prediction - 0.13, r'$\overline{\Delta G}$')

    # Plot experimental value horizontal line only for generic plot.
    if plot_experimental_value:
        experiment_dg = -11.75
        plt.plot(x_lim, [experiment_dg]*2, color='black')

    if cost == 'neval':
        ax.set_xlabel('N force/energy evaluations')
    else:
        ax.set_xlabel('Computational cost', labelpad=-5)
    ax.set_ylabel('$\Delta$G', labelpad=-5)
    ax.set_yticklabels([])
    ax.set_xticklabels([])

    plt.tight_layout(pad=0.1, rect=[0.0, 0.0, 0.90, 1.0])

    # Save file.
    figure_dir_path = os.path.join(SAMPLING_PAPER_DIR_PATH, 'Figure 1 - host-guest')
    os.makedirs(figure_dir_path, exist_ok=True)
    output_base_path = os.path.join(figure_dir_path, 'example_trajectories')
    plt.savefig(output_base_path + '.pdf')


# =============================================================================
# FIGURE 2 - MEAN ERROR AND RELATIVE EFFICIENCY CARTOON
# =============================================================================

def plot_mean_error_cartoon():
    """Plot the cartoon used to explain mean error and relative efficiency.

    This is used as an example to clarify some gotchas with the difinition
    of efficiency.
    """
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    sns.set_context('paper')
    sns.set_style('white')

    def err_decay_func_square(decay_coeff, c):
        return decay_coeff / np.sqrt(c)

    def mean_error_square(decay_coeff, c_min, c_max):
        return 2 * decay_coeff * (np.sqrt(c_max) - np.sqrt(c_min)) / (c_max - c_min)

    def err_decay_func_B(decay_coeff, c):
        return decay_coeff / c**(5/6)

    def mean_error_B(decay_coeff, c_min, c_max):
        return 6 * decay_coeff * (c_max**(1/6) - c_min**(1/6)) / (c_max - c_min)

    decay_coeffs = {
        'A': 1.0,
        'B': 2.5,
        'Z': 1.5,
    }
    c_ranges = collections.OrderedDict([
        ("A'", np.arange(1, 4.5, 0.1)),
        ("A''", np.arange(3, 6, 0.1)),
        ("B", np.arange(2, 6.5, 0.1)),
        ("Z", np.arange(1, 6.5, 0.1)),
    ])

    # Determine colors colors.
    colors = {m: 'C'+str(i) for i, m in enumerate(sorted(c_ranges))}

    # Plot the error trajectories.
    fig, ax = plt.subplots(figsize=(3.5, 2.6))

    # method_names = ["B", "Z", "A'", "A''"]
    method_names = ["Z", "A'", "A''"]
    for method_name in method_names:
        color = colors[method_name]
        c_range = c_ranges[method_name]
        decay_coeff = decay_coeffs[method_name[0]]

        if method_name == 'B':
            err_decay_func = err_decay_func_B
        else:
            err_decay_func = err_decay_func_square
        err = err_decay_func(decay_coeff, c_range)

        # Plot error area.
        ax.plot(c_range, err, color=color, label=method_name, zorder=1)
        ax.fill_between(c_range, err, 0, color=color, alpha=0.5, zorder=0)

        # Add method label.
        c_method_label_idx = int(len(c_range) / 8)
        ax.text(c_range[c_method_label_idx], err[c_method_label_idx]+0.01, method_name, fontsize=12)

        if method_name[0] == 'A':
            # Plot mean error.
            c_min, c_max = min(c_range), max(c_range)
            mean_err = mean_error_square(decay_coeff, c_min, c_max)
            # Start mean error horizontal line from the error curve.
            c_mean = (decay_coeff / mean_err)**2
            ax.plot([0, c_mean], [mean_err, mean_err], color='black', ls='--', alpha=0.8, zorder=1)

            # Add label mean error.
            # ax.text(1.05, mean_err+0.025, '$\mathbb{E}[RMSE_{' + method_name + '}]$', fontsize=9)
            ax.text(-0.3, mean_err+0.025, '$\mathbb{E}[RMSE_{' + method_name + '}]$', fontsize=9)

            # Add c_min/max labels.
            ax.text(c_min-0.4, -0.1, 'c$_{min,' + method_name + '}$', fontsize=9)
            ax.text(c_max-0.4, -0.1, 'c$_{max,' + method_name + '}$', fontsize=9)

    # Configure axes.
    ax.set_xlim(1, 6.4)
    ax.set_ylim(0, 2)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_ylabel('$RMSE(\Delta G)$')
    ax.set_xlabel('computational cost')
    # Pull axes labels closest to axes.
    ax.tick_params(axis='x', which='major', pad=2.0)
    ax.yaxis.set_label_coords(0.0, 0.65)

    # Plot the relative efficiencies in an inset plot.
    ax_ins = inset_axes(ax, width='100%', height='100%', bbox_to_anchor=[145, 115, 90, 50])

    # Compute relative efficiencies with respect to Z.

    relative_efficiencies = collections.OrderedDict()
    for method_name in [name for name in method_names if name != 'Z']:
        c_min, c_max = min(c_ranges[method_name]), max(c_ranges[method_name])
        if method_name == 'B':
            mean_error_func = mean_error_B
        else:
            mean_error_func = mean_error_square
        mean_err_method = mean_error_func(decay_coeffs[method_name[0]], c_min, c_max)
        mean_err_Z = mean_error_square(decay_coeffs['Z'], c_min, c_max)
        relative_efficiencies[method_name] = -np.log(mean_err_method/mean_err_Z)

    # Plot horizontal bar plot with all efficiencies.
    labels, rel_effs = zip(*relative_efficiencies.items())
    bar_colors = [colors[m] for m in labels]
    labels = [l + '/Z' for l in labels]
    # labels = ['$e_{err,' + str(l) + '/Z}$' for l in labels]
    ax_ins.barh(y=labels, width=rel_effs, color=bar_colors, alpha=0.85)
    ax_ins.set_title('relative efficiency', pad=2.5)

    # plt.tight_layout(rect=[0.0, 0.0, 1.0, 1.0])
    plt.tight_layout(rect=[0.1, 0.0, 1.0, 1.0])

    # Pull axes labels closest to axes.
    ax_ins.set_xticks([0.0])
    ax_ins.grid(axis='x')
    ax_ins.tick_params(axis='x', which='major', pad=0.0)
    ax_ins.tick_params(axis='y', which='major', pad=0.0)

    output_dir_path = os.path.join(SAMPLING_PAPER_DIR_PATH, 'Figure2-efficiency_cartoon')
    os.makedirs(output_dir_path, exist_ok=True)
    plt.savefig(os.path.join(output_dir_path, 'error_trajectories.pdf'))


# =============================================================================
# FIGURE 3 - FREE ENERGY TRAJECTORIES
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
        submission_ls = SUBMISSION_LINE_STYLES[method_name]

        # Plot free energy trajectories.
        for system_name, mean_data in submission_mean_data.items():
            ax_idx = system_columns[system_name]

            # The OA prediction of the NS short protocol are the same of the long protocol submission file.
            if method_name == 'GROMACS/NS-DS/SB-long' and system_name != 'CB8-G3':
                # Just add the label.
                axes[0][ax_idx].plot([], color=submission_color, ls=submission_ls, label=method_name)
                continue

            # Update maximum number of energy evaluations.
            n_energy_evaluations = mean_data['N energy evaluations'].values[-1]
            max_n_energy_evaluations[system_name] = max(max_n_energy_evaluations[system_name],
                                                        n_energy_evaluations)

            # Determine zorder and plot.
            zorder = zorders[system_name][method_name]
            plot_mean_data(mean_data, axes[:,ax_idx], color=submission_color,
                           ls=submission_ls, zorder=zorder, label=method_name,
                           plot_std=plot_std, plot_bias=plot_bias)

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
    sns.set_style('whitegrid')
    sns.set_context('paper')

    # Create a figure with 3 columns (one for each system) and 2 rows.
    # The first row contains the free energy trajectory and CI, the second
    # a plot of the estimator variance, and the third the bias to the
    # asymptotic value.
    if zoomed:
        figsize = (7.25, 7.0)  # Without REVO
    else:
        figsize = (7.25, 7.0)  # With REVO
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=figsize)

    # Optionally, remove REVO.
    if zoomed:
        submissions = [s for s in submissions if s.name not in ['WExploreRateRatio']]

    if zoomed:
        # Y-axis limits when REVO calculations are excluded.
        y_limits = [
            [(-15, -10), (-9, -4), (-9, -4)],
            [(0, 2), (0, 0.8), (0, 0.8)],
            [(-3, 1), (-0.6, 0.6), (-0.6, 0.6)],
        ]
    else:
        # Y-axis limits when REVO calculations are included.
        y_limits = [
            [(-17, -9), (-13, -5), (-13, -5)],
            [(0, 2), (0, 1.75), (0, 1.75)],
            [(-4, 4), (-0.6, 0.6), (-0.6, 0.6)],
        ]

    plot_submissions_trajectory(submissions, yank_analysis, axes, y_limits=y_limits)

    # Show/save figure.
    if zoomed:
        plt.tight_layout(h_pad=0.2, rect=[0.0, 0.00, 1.0, 0.92], w_pad=0.0)  # Without REVO
    else:
        plt.tight_layout(h_pad=0.2, rect=[0.0, 0.00, 1.0, 0.92])  # With REVO

    # Plot legend.
    if zoomed:
        # bbox_to_anchor = (2.52, 1.55)  # Without REVO.
        bbox_to_anchor = (2.4, 1.48)
    else:
        bbox_to_anchor = (2.4, 1.48)  # With REVO.
    axes[0][1].legend(loc='upper right', bbox_to_anchor=bbox_to_anchor,
                      fancybox=True, ncol=4)
    plt.subplots_adjust(wspace=0.35)
    # plt.show()
    if zoomed:
        file_name = 'Figure3-free_energy_trajectories_zoomed'
    else:
        file_name = 'Figure3-free_energy_trajectories'
    figure_dir_path = os.path.join(SAMPLING_PAPER_DIR_PATH, 'Figure3-free_energy_trajectories')
    os.makedirs(figure_dir_path, exist_ok=True)
    output_base_path = os.path.join(figure_dir_path, file_name)
    plt.savefig(output_base_path + '.pdf')
    # plt.savefig(output_base_path + '.png', dpi=500)


# =============================================================================
# FIGURE 4 - NONEQUILIBRIUM SWITCHING ESTIMATOR COMPARISON
# =============================================================================

def plot_all_nonequilibrium_switching(submissions):
    """Plot free energy trajectories, std, and bias of the nonequilibrium-switching calculations."""
    # Create a figure with 3 columns (one for each system) and 2 rows.
    # The first row contains the free energy trajectory and CI, the second
    # a plot of the estimator variance, and the third the bias to the
    # asymptotic value.
    figsize = (7.25, 3.5)
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
    legend = axes[0].legend(loc='upper left', bbox_to_anchor=(0.6, 1.3),
                            fancybox=True, ncol=3)
    # Change legend labels to refer to estimator used rather than overall method ID.
    legend_labels_map = {
        'GROMACS/NS-DS/SB-long': 'BAR-long',
        'GROMACS/NS-DS/SB': 'BAR',
        'GROMACS/NS-Jarz-F': 'Jarzynski-Forward',
        'GROMACS/NS-Jarz-R': 'Jarzynski-Reverse',
        'GROMACS/NS-Gauss-F': 'Gaussian-Forward',
        'GROMACS/NS-Gauss-R': 'Gaussian-Reverse',
    }
    for text in legend.get_texts():
        text.set_text(legend_labels_map[text.get_text()])

    plt.subplots_adjust(wspace=0.35)

    # plt.show()
    figure_dir_path = os.path.join(SAMPLING_PAPER_DIR_PATH, 'Figure4-nonequilibrium_comparison')
    os.makedirs(figure_dir_path, exist_ok=True)
    output_base_path = os.path.join(figure_dir_path, 'Figure4-nonequilibrium_comparison')
    plt.savefig(output_base_path + '.pdf')
    # plt.savefig(output_base_path + '.png', dpi=500)


# =============================================================================
# FIGURE 5 - BAROSTAT AND RESTRAINT
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
    mc_volumes = collections.OrderedDict([
        (1, np.load(os.path.join(YANK_VOLUMES_DIR_PATH, 'volumes_pressure100.npy'))),
        (100, np.load(os.path.join(YANK_VOLUMES_DIR_PATH, 'volumes_pressure10000.npy'))),
    ])
    mc_volumes_hrex = collections.OrderedDict([
        (1, np.load(os.path.join(YANK_VOLUMES_DIR_PATH, 'hrex_state_volumes_state0.npy'))),
        (58, np.load(os.path.join(YANK_VOLUMES_DIR_PATH, 'hrex_state_volumes_state58.npy'))),
    ])

    b_volumes = collections.OrderedDict([
        (1, np.load(os.path.join(EE_VOLUMES_DIR_PATH, '1atm_vanilla.npy'))),
        (100, np.load(os.path.join(EE_VOLUMES_DIR_PATH, '100atm_vanilla.npy'))),
    ])
    b_volumes_ee = collections.OrderedDict([
        (1, np.load(os.path.join(EE_VOLUMES_DIR_PATH, '1atm_expanded.npy'))),
        (100, np.load(os.path.join(EE_VOLUMES_DIR_PATH, '100atm_expanded.npy'))),
    ])

    # Print some statistics for each distribution.
    for volume_trajectories, label in [(mc_volumes, 'MC-MD  '),
                                       (mc_volumes_hrex, 'MC-HREX'),
                                       (b_volumes, 'BB-MD    '),
                                       (b_volumes_ee, 'BB-EE    ')]:
        for pressure, trajectory in volume_trajectories.items():
            n = len(trajectory)
            t_stat = 2.326  # 98% CI

            mean = np.mean(trajectory)
            sem = scipy.stats.sem(trajectory)
            mean_ci = t_stat * sem

            var = np.var(trajectory, ddof=1)
            # Standard error of variance if volume is gaussianly distributed
            sev = var * np.sqrt(2 / (n-1))
            var_ci = t_stat * sev

            skew = scipy.stats.skew(trajectory)
            # Standard error of skewness if volume is gaussianly distributed
            ses = np.sqrt( 6*n*(n-1) / ((n-2)*(n+1)*(n+3)) )
            skew_ci = t_stat * ses
            print('{}-{} (n={}): mean={:.3f} +- {:.3f}nm^3\t\tvar={:.3f} +- {:.3f}\tskew={:.3f} +- {:.3f}'.format(
                pressure, label, n, mean, mean_ci, var, var_ci, skew, skew_ci))

    # Plot the 1atm vs 100atm comparison.
    barostats = ['B', 'MC']
    for ax, volume_trajectories, barostat in zip(axes, [b_volumes, mc_volumes], barostats):
        barostat += ',MD'
        barostat = 'MD'

        for pressure, trajectory in volume_trajectories.items():
            label = '$\\rho_{{\mathrm{{{}}}}}$(V|{}atm)'.format(barostat, pressure)
            ax = sns.distplot(trajectory, label=label, hist=False, ax=ax)

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
            label = '$\\rho_{{\mathrm{{{}}}}}$(V|{}atm)$\cdot e^{{\\beta ({}atm - {}atm) V}}$'.format(barostat, original_pressure, new_pressure, original_pressure)
            ax.plot(volumes, predicted, ls='--', label=label)
            # ax.plot(volumes, [fit.pdf([v], *fit_parameters) for v in volumes], label='original')

    # Plot comparison MD vs expanded ensemble and HREX volumes.
    for ax_idx, (trajectory, label) in enumerate([
        (b_volumes_ee[1], 'B,EE'), (mc_volumes_hrex[1], 'MC,HREX')
    ]):
        label = 'E'
        ax = axes[ax_idx]
        label = '$\\rho_{{\mathrm{{{}}}}}$(V|1atm)'.format(label)
        sns.distplot(trajectory, label=label, hist=False, ax=ax)

    # Set titles and configure axes.
    axes[0].set_title('Berendsen barostat volume distribution', pad=2.0)
    axes[1].set_title('Monte Carlo barostat volume distribution', pad=2.0)
    for ax_idx in range(len(axes)):
        axes[ax_idx].set_xlim((78.8, 81.2))
        axes[ax_idx].set_ylim((0.0, 6.0))
        axes[ax_idx].set_ylabel('density')
    axes[0].set_xlabel('', labelpad=0.3)
    axes[1].set_xlabel('Volume [nm^3]', labelpad=0.3)

    # Create single legend for both MC and B barostat axes.
    bbox_to_anchor = (-0.1, -0.15)
    axes[0].legend(fontsize='xx-small', loc='upper left', bbox_to_anchor=bbox_to_anchor, ncol=4,
                   fancybox=True, labelspacing=0.7, handletextpad=0.4, columnspacing=1.1,)
    # axes[0].get_legend().remove()
    axes[1].get_legend().remove()

    plt.tight_layout(pad=0, rect=[0.0, 0.0, 1.0, 1.0])


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
    ax.set_title('Restrained ligand-receptor distance', pad=2.0)
    if kde is False:
        ax.set_ylabel('Number of samples')
    else:
        ax.set_ylabel('density')
    ax.legend(loc='upper right', fontsize='x-small')
    ax.set_xlabel('Restrained distance [$\mathrm{\AA}$]', labelpad=0.3)

    # Free energy as a function of restraint distance.
    ax = axes[1]
    ax.set_title('$\Delta G$ as a function of restraint radius cutoff', pad=2.0 )
    plot_restraint_profile(system_id, ax, restraint_cutoff)
    # Labels and legend.
    ax.set_xlabel('Restraint radius cutoff [$\mathrm{\AA}$]', labelpad=0.3)
    ax.set_ylabel('$\Delta G$ [kcal/mol]')
    ax.legend(fontsize='x-small')


def plot_restraint_and_barostat_analysis():
    """Plot the Figure showing info for the restraint and barostat analysis."""
    import seaborn as sns
    from matplotlib import pyplot as plt
    sns.set_style('whitegrid')
    sns.set_context('paper', font_scale=1.0)

    # Create two columns, each of them share the x-axis.
    fig = plt.figure(figsize=(7.25, 4))
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
    restraint_axes[1].set_ylim((-7, -3.9))

    for ax in restraint_axes + barostat_axes:
        ax.tick_params(axis='x', which='major', pad=0.1)
        ax.tick_params(axis='y', which='major', pad=0.1)
    plt.tight_layout(pad=0.3)

    # plt.show()
    output_file_path = os.path.join(SAMPLING_PAPER_DIR_PATH, 'Figure5-restraint_barostat',
                                    'restraint_barostat.pdf')
    os.makedirs(os.path.dirname(output_file_path), exist_ok=True)
    plt.savefig(output_file_path)


# =============================================================================
# FIGURE 6 - HREX INITIAL BIAS
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
        # ('CB8-G3', True),
        ('OA-G3', False),
        # ('OA-G3', False),
        ('OA-G6', False),
    ]

    if plot_std:
        n_rows = 2
    else:
        n_rows = 1
    n_cols = len(what_to_plot)
    fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(7.25, 4.0))

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

    # If there is an odd number of columns print x label only on the central one.
    if n_cols % 2 == 1:
        axes[-1][1].set_xlabel('HREX iteration')
    else:
        for ax in axes[-1]:
            ax.set_xlabel('HREX iteration')

    plt.tight_layout(h_pad=0.1, rect=[0.0, 0.00, 1.0, 0.91])

    handles, labels = axes[0][0].get_legend_handles_labels()
    handles = [handles[-1]] + handles[:-1]
    labels = [labels[-1]] + labels[:-1]
    bbox_to_anchor = (0.4, 1.53)
    axes[0][0].legend(handles, labels, loc='upper left', bbox_to_anchor=bbox_to_anchor,
                      title='number of discarded initial iterations', ncol=len(data_dir_paths)+1,
                      fancybox=True, labelspacing=0.8, handletextpad=0.5, columnspacing=1.2,
                      fontsize='small')

    # plt.show()
    if figure_dir_path is None:
        figure_dir_path = os.path.join(SAMPLING_PAPER_DIR_PATH, 'Figure6-bias_hrex')
    os.makedirs(figure_dir_path, exist_ok=True)
    output_file_path = os.path.join(figure_dir_path, 'Figure6-bias_hrex')
    plt.savefig(output_file_path + '.pdf')
    # plt.savefig(output_file_path + '.png', dpi=600)


# =============================================================================
# SUPPORTING INFORMATION - EXAMPLE OF HREX BIAS
# =============================================================================

def simulate_correlation_samples():
    """Simulation of bias from same initial configuration.

    There are 3 states as different harmonic oscillators, but all
    or almost all the samples come from the first (bound) state to
    simulate what happens when they don't decorrelate fast enough.
    The hypothesis is that most is that starting from the bound
    state causes the initial free energy to be artificially negative
    if the correlation times are long.

    The second (discharged) state is just a shifted harmonic oscillator
    (same free energy as bound state). The third (unbound) is shifted
    and has much higher entropy.

    """
    from numpy.random import normal
    from pymbar import MBAR

    def harmonic_oscillator_free_energy(sigma):
        """Analytical expression for the free energy of a harmonic oscillator."""
        #return - np.log(2 * np.pi * sigma**2) * 3.0 / 2.0  # 3D oscillator
        return - np.log(np.sqrt(2 * np.pi) * sigma)

    def harmonic_oscillator_potential(x, loc, std):
        """Compute potential of the given positions given location
        and standard deviation of the Gaussian distribution.

        Potentials are returned in units of kT.
        """
        spring_constant = 1 / std**2
        return spring_constant / 2.0 * (x - loc)**2

    def print_free_energies(Deltaf_ij, dDeltaf_ij):
        mbar_str = ', '.join(['{:.4f} +- {:.4f}'.format(f, df) for f, df in zip(Deltaf_ij[:,0], dDeltaf_ij[:,0])])
        print('MBAR      :', mbar_str)
        analytical_str = ', '.join(['{:.4f}          '.format(f) for f in analytical_Deltaf])
        print('Analytical:', analytical_str)

    def compute_mbar_free_energy(all_samples, shifts, stds, analytical_f):
        n_states = len(all_samples)

        # u_kn[k,n] is the reduced potential energy n-th sample evaluated at state k.
        u_kn = np.empty(shape=(n_states, n_states*n_samples))

        # Convert samples to potentials.
        for k in range(n_states):
            for sampled_k, samples in enumerate(all_samples):
                start = sampled_k * n_samples
                end = (sampled_k + 1) * n_samples
                u_kn[k,start:end] = harmonic_oscillator_potential(samples, loc=shifts[k], std=stds[k])

        # Compute MBAR free energy.
        N_k = np.array([n_samples] * n_states)
        mbar = MBAR(u_kn, N_k=N_k, initial_f_k=analytical_f)
        Deltaf_ij, dDeltaf_ij, _ = mbar.getFreeEnergyDifferences()
        return Deltaf_ij, dDeltaf_ij


    # Determine standard deviation and shift of the harmonic distributions.
    n_samples = 5000000
    stds = np.array([2.0, 2.0, 5.0])
    shifts = np.array([0.0, 2.0, 2.0])
    print('\nspring constants:', 1 / stds**2)

    # Compute analytical free energy.
    analytical_f = np.array([harmonic_oscillator_free_energy(s) for s in stds])
    analytical_Deltaf = np.array([analytical_f[0] - analytical_f[i] for i in range(len(stds))])

    # FIRST TEST.
    # Sample from all states and verify that MBAR free energy is correct.
    # -------------------------------------------------------------------
    all_samples = [normal(loc=l, scale=s, size=n_samples) for l, s in zip(shifts, stds)]
    Deltaf_ij, dDeltaf_ij = compute_mbar_free_energy(all_samples, shifts, stds, analytical_f)
    print()
    print_free_energies(Deltaf_ij, dDeltaf_ij)

    # SECOND TEST.
    # Check if the bias is not due to lack of overlap. If we sample only the end states the estimate should be correct.
    # -----------------------------------------------------------------------------------------------------------------
    for i in range(1, len(all_samples)):
        all_samples_bar = [all_samples[0], all_samples[i]]
        shifts_bar = [shifts[0], shifts[i]]
        stds_bar = [stds[0], stds[i]]
        analytical_f_bar = [analytical_f[0], analytical_f[i]]
        Deltaf_ij, dDeltaf_ij = compute_mbar_free_energy(all_samples_bar, shifts_bar, stds_bar, analytical_f_bar)
        print('\nBAR_{}0'.format(i))
        print_free_energies(Deltaf_ij, dDeltaf_ij)

    # THIRD TEST.
    # Now sample from only the bound state to see how the free energy changes.
    # ------------------------------------------------------------------------
    all_samples[1:] = [normal(loc=shifts[0], scale=stds[0], size=n_samples) for _ in range(len(stds)-1)]
    Deltaf_ij, dDeltaf_ij = compute_mbar_free_energy(all_samples, shifts, stds, analytical_f)
    print()
    print_free_energies(Deltaf_ij, dDeltaf_ij)

    # FOURTH TEST.
    # Now let the unbound state decorrelate fast (i.e. sample from its own distribution).
    # -----------------------------------------------------------------------------------
    all_samples[-1] = normal(loc=shifts[-1], scale=stds[-1], size=n_samples)
    Deltaf_ij, dDeltaf_ij = compute_mbar_free_energy(all_samples, shifts, stds, analytical_f)
    print()
    print_free_energies(Deltaf_ij, dDeltaf_ij)

    # RESULT: SUCCESS!!!


# =============================================================================
# SUPPORTING INFORMATION - COMPLEX/SOLVENT and ENTROPY/ENTHALPY DECOMPOSITION
# =============================================================================

def _mean_data_decomposition(data):
    # Convert into a numpy array to take the mean.
    # Convert None (not supported by numpy) into nans.
    try:
        # This may fail if we have computed different iterations for each.
        data = np.array(data, dtype=np.float)
    except ValueError:
        data_lengths = [len(x) for x in data]
        print('Warning: Truncating data of shape {}'.format(data_lengths))
        min_length = min(data_lengths)
        data = [x[:min_length] for x in data]
        data = np.array(data, dtype=np.float)
    # Compute std and mean along the trajectory ignoring NaNs.
    return np.nanmean(data, axis=0), np.nanstd(data, axis=0)


def _plot_phase_decomposition(ax, phase_free_energies):
    # Shortcuts.
    data = phase_free_energies
    label = '$\Delta$G'

    # Plot each phase data on a separate axis to make the comparison on different order of magnitudes easier.
    # Receipt with three axes: https://matplotlib.org/3.1.0/gallery/ticks_and_spines/multiple_yaxis_with_spines.html
    phase_axes = {
        'complex': ax.twinx(),
        'solvent': ax.twinx()
    }
    phase_colors = {
        'complex': 'C1',
        'solvent': 'C0',
    }
    for ax_name in sorted(phase_axes):
        phase_axes[ax_name].set_ylabel(label + ' ' + ax_name + ' [kcal/mol]',
                                       color=phase_colors[ax_name])
    phase_axes[ax_name].spines["right"].set_position(("axes", 1.2))

    # Compute total free energy summing complex and solvent for all replicates.
    total_mean = [np.array(data['solvent'][i]) + np.array(data['complex'][i]) for i in range(5)]
    total_mean, total_std = _mean_data_decomposition(total_mean)

    # Compute and plot the phase free energy.
    for phase_name in ['complex', 'solvent']:
        color = phase_colors[phase_name]

        # Convert into a numpy array to take the mean.
        # Convert None (not supported by numpy) into nans.
        data[phase_name], std = _mean_data_decomposition(data[phase_name])

        # Plot each phase data on a separate axis to make the comparison easier.
        phase_axes[phase_name].plot(data[phase_name], ls='-', color=color,
                                    label=label + ' ' + phase_name)
        # Plot uncertainties.
        phase_axes[phase_name].fill_between(x=list(range(len(std))), y1=data[phase_name]-std,
                                            y2=data[phase_name]+std, color=color, alpha=0.7)


    # Plot total free energy.
    # total = data['solvent'] + data['complex']
    # ax.plot(total, color='black', label=label+' total')
    ax.plot(total_mean, color='black', label=label+' total')
    ax.fill_between(x=list(range(len(total_std))), y1=total_mean-total_std,
                    y2=total_mean+total_std, color='black', alpha=0.7)
    ax.set_ylabel(label + ' total [kcal/mol]')
    ax.set_xlabel('simulation percentage')

    # Make the range of all y axes the same.
    ax.set_ylim((-21, -18))
    phase_axes['complex'].set_ylim((-151.0, -148.0))
    phase_axes['solvent'].set_ylim((129.0, 132.0))


def _plot_entropy_enthalpy_decomposition(ax, phase_free_energies, phase_enthalpy):
    # Analyze only the complex.
    phase_name = 'complex'

    # Plot each phase data on a separate axis to make the comparison on different order of magnitudes easier.
    # Receipt with three axes: https://matplotlib.org/3.1.0/gallery/ticks_and_spines/multiple_yaxis_with_spines.html
    axes = {
        '$\Delta$G': ax,
        '$\Delta$H': ax.twinx(),
        '-T$\Delta$S': ax.twinx(),
    }
    colors = {
        '$\Delta$G': 'black',
        '$\Delta$H': 'C1',
        '-T$\Delta$S': 'C0',
    }
    for ax_name in sorted(axes):
        axes[ax_name].set_ylabel(ax_name + ' ' + phase_name + ' [kcal/mol]', color=colors[ax_name])
    axes[ax_name].spines["right"].set_position(("axes", 1.2))

    # Variable used to propagate entropy decomposition.
    entropy_std = []

    # Plot the total average free energy and enthalpy and for each phase.
    for data, label in [(phase_free_energies, '$\Delta$G'),
                        (phase_enthalpy, '$\Delta$H')]:
        color = colors[label]

        # Convert into a numpy array to take the mean.
        # Convert None (not supported by numpy) into nans.
        data[phase_name], std = _mean_data_decomposition(data[phase_name])
        ns_replica = np.arange(0.0, 40.0, 40/len(std))

        # Plot each phase data on a separate axis to make the comparison easier.
        axes[label].plot(ns_replica, data[phase_name], ls='-', color=color, label=label+' '+phase_name)
        # Plot uncertainties.
        axes[label].fill_between(x=ns_replica, y1=data[phase_name]-std,
                                 y2=data[phase_name]+std, color=color, alpha=0.7)

        # Propagate uncertainty.
        if len(entropy_std) == 0:
            entropy_std = std**2
        else:
            entropy_std += std**2
    entropy_std = np.sqrt(entropy_std)

    # Plot also entropies.
    label = '-T$\Delta$S'
    color = colors[label]
    entropy = phase_free_energies[phase_name] - phase_enthalpy[phase_name]
    axes[label].plot(ns_replica, entropy, ls='-', color=color, label=label+' '+phase_name)
    # Plot uncertainties.
    axes[label].fill_between(x=ns_replica, y1=entropy-entropy_std,
                             y2=entropy+entropy_std, color=color, alpha=0.7)

    ax.set_xlabel('ns/replica')


def plot_decomposition(system_name, starting_iteration, type, output_file_path):
    """
    Decomposition of the free energy trajectory in complex/solvent phase or entropy/enthalpy.

    Parameters
    ----------
    type : str
        Can be 'entropy-enthalpy' or 'phase'.
    """
    data_file_pattern = 'YankAnalysis/BiasAnalysis/iter{}/fe-decomposition-{}-{{}}.json'.format(
        starting_iteration, system_name)

    n_replicates = 5
    phase_free_energies = {'complex': [[] for _ in range(n_replicates)],
                           'solvent': [[] for _ in range(n_replicates)]}
    phase_enthalpy = copy.deepcopy(phase_free_energies)

    for replicate_idx in range(n_replicates):
        # Read decomposition data.
        decomposition_data_file_path = data_file_pattern.format(replicate_idx)
        with open(decomposition_data_file_path, 'r') as f:
            decomposition_data = json.load(f)

        # Read free energy and enthalpy at each iteration.
        sorted_decomposition_data = sorted(decomposition_data, key=lambda x: int(x.split('-')[1]))
        for phase_iter in sorted_decomposition_data:
            decomposition = decomposition_data[phase_iter]
            phase_name, iteration = phase_iter.split('-')

            # Correct sign consistent with thermodynamic cycle.
            if phase_name == 'complex':
                sign = -1
            else:
                sign = 1

            corrected_free_energy = sign * (decomposition['DeltaF'] + decomposition['DeltaF_standard_state_correction'])
            phase_free_energies[phase_name][replicate_idx].append(corrected_free_energy)

            # Multiplication works only if enthalpy is not None.
            if decomposition['DeltaH'] is not None:
                decomposition['DeltaH'] *= sign
            phase_enthalpy[phase_name][replicate_idx].append(decomposition['DeltaH'])

    # Create figure.
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.25, 4.6))
    if type == 'entropy-enthalpy':
        _plot_entropy_enthalpy_decomposition(ax, phase_free_energies, phase_enthalpy)
    else:
        _plot_phase_decomposition(ax, phase_free_energies)

        # # Plot total free energy.
        # total = data['solvent'] + data['complex']
        # ax.plot(total, color=color, label=label)
        # totals.append(total)

    # Plot also entropies.
    # ax.plot(totals[0] - totals[1], color='blue', label='-T$\Delta$S')

    # ax.set_ylim((-20, -18))
    # phase_axes['complex'].set_ylim((-153, -148))
    # phase_axes['solvent'].set_ylim((128, 133))
    # ax.set_ylim((-23, -18))
    # phase_axes['complex'].set_ylim((30, 45))
    # phase_axes['solvent'].set_ylim((-55, -40))

    # ax.legend()

    plt.tight_layout()
    if output_file_path is not None:
        os.makedirs(os.path.dirname(output_file_path), exist_ok=True)
        plt.savefig(output_file_path)
    else:
        plt.show()


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

    # Discard the initial frames of REVO and GROMACS/EE that don't have predictions.
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
    sns.set_style('whitegrid')
    sns.set_context('paper')

    statistic_names = ['std', 'absolute bias', 'RMSE']

    # Create output directory.
    figure_dir_path = os.path.join(SAMPLING_PAPER_DIR_PATH, 'SI_Figure-efficiencies')
    os.makedirs(figure_dir_path, exist_ok=True)

    # Check if we need all the efficiencies in the same plot or not.
    if same_plot:
        fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(7.25, 8))
        # Keep track of data range by statistic.
        statistic_ranges = {name: [np.inf, 0] for name in statistic_names}
    # Keep track of n_energy_evaluations by column.
    max_n_energy_evaluations = [0 for _ in range(3)]

    for submission in submissions:
        if submission.paper_name in {'OpenMM/REVO'}:
            continue
        # if submission.paper_name in {'AMBER/APR', 'GROMACS/NS-DS/SB', 'GROMACS/NS-DS/SB-long',
        #                              'NAMD/BAR', 'GROMACS/EE', 'GROMACS/EE-fullequil', 'OpenMM/SOMD'}:
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
            color = SUBMISSION_COLORS[submission.paper_name]

            # For GROMACS/EE, there are no submissions for CB8-G3.
            if 'GROMACS/EE' in submission.paper_name and system_name == 'CB8-G3':
                continue
            # For GROMACS/NS-DS/SB-long there are no new submissions for OAs.
            if 'GROMACS/NS-DS/SB-long' in submission.paper_name and system_name != 'CB8-G3':
                # Just add the label.
                axes[0][col_idx].plot([], color=color, label=submission.paper_name)
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
            axes[-1][1].set_xlabel('Number of force/energy evaluations [10$^6$]')

        # Set labels and axes limits.
        if not same_plot:
            fig.suptitle(submission.paper_name)

            output_file_base_name = 'releff-{}-{}'.format(submission.file_name, submission.receipt_id)
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
    sns.set_style('whitegrid')
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
        if 'REVO' in submission.paper_name:
            continue
        print(submission.paper_name)

        # Obtain std, bias, and RMSE of the 5 trajectories.
        # If this is a YANK analysis, we get it later specifically for the system.
        if not isinstance(submission, YankSamplingAnalysis):
            mean_free_energies = submission.mean_free_energies()

        color = SUBMISSION_COLORS[submission.paper_name]

        for col_idx, system_name in enumerate(system_names):
            # GROMACS/EE doesn't have submissions for CB8-G3.
            if 'GROMACS/EE' in submission.paper_name and system_name == 'CB8-G3':
                continue
            # For GROMACS/NS-DS/SB-long there are no new submissions for OAs.
            if 'GROMACS/NS-DS/SB-long' in submission.paper_name and 'OA' in system_name:
                # Just add the label.
                axes[0][col_idx].plot([], color=color, label=submission.paper_name)
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

    figure_dir_path = os.path.join(SAMPLING_PAPER_DIR_PATH, 'SI_Figure-efficiencies')
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
    statistic_names = [r'$e_{\mathrm{std}}$', r'$e_{|\mathrm{bias}|}$', r'$e_{\mathrm{RMSD}}$']
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
                    (submission.paper_name == 'GROMACS/NS-DS/SB-long' and system_name != 'CB8-G3')):
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


def print_nonequilibrium_relative_efficiencies(nonequilibrium_submissions):
    """Print relative efficiencies w.r.t. for the nonequilibrium estimators table."""
    system_names = ['CB8-G3', 'OA-G3', 'OA-G6']

    def _get_free_energy_array(submission, system_name, step=1, max_c=100, get_asymptotic=False):
        n_replicates = 5
        system_data = submission.data[submission.data['System name'] == system_name]
        free_energy_array = np.empty(shape=(n_replicates, int(max_c/step)))
        for i in range(n_replicates):
            system_id = system_name + '-' + str(i)
            system_id_data = system_data[system_data['System ID'] == system_id]
            free_energy_array[i] = system_id_data[DG_KEY].values[:max_c:step]

        if get_asymptotic:
            mean_free_energies = submission.mean_free_energies()
            asymptotic = mean_free_energies[mean_free_energies['System name'] == system_name][DG_KEY].values[-1]
            return free_energy_array, asymptotic
        return free_energy_array

    # Use GROMACS/NS-DS/SB-long as reference method.
    reference_submission = [s for s in nonequilibrium_submissions if s.paper_name == 'GROMACS/NS-DS/SB-long'][0]
    # Also remove the other BAR submission.
    nonequilibrium_submissions = [s for s in nonequilibrium_submissions if 'GROMACS/NS-DS/SB' not in s.paper_name]

    # Get only the first 50 as the 1-directional estimators only have half the cost.
    free_energy_ref = {}
    asymptotic_ref = {}
    for system_name in system_names:
        DG, asympt = _get_free_energy_array(reference_submission, system_name, max_c=50, get_asymptotic=True)
        free_energy_ref[system_name] = DG
        asymptotic_ref[system_name] = asympt

    for submission in nonequilibrium_submissions:
        print(submission.paper_name, end='')
        for system_name in system_names:
            free_energy_sub = _get_free_energy_array(submission, system_name, step=2)
            rel_eff, cis = compute_all_relative_efficiencies(
                free_energy_ref[system_name], free_energy_sub, ci=0.95, n_bootstrap_samples=1000,
                asymptotic_free_energy_A=asymptotic_ref[system_name],
                asymptotic_free_energy_B=asymptotic_ref[system_name]
            )
            for i, stat_name in enumerate(['std', 'bias', 'RMSE']):
                print(r' & {:.1f}$_{{\raisem{{2pt}}{{{:.1f}}}}}^{{\mathstrut {:.1f}}}$'.format(rel_eff[i], cis[i][0], cis[i][1]), end='')
        print(r' \\')


def print_final_prediction_table(submissions, yank_analysis):
    """Plot the table containing the fina binding free energy predictions for all replicates."""
    for submission in submissions + [yank_analysis]:
        # GROMACS/EE-fullequil predictions are identical to GROMACS/EE
        if submission.paper_name == 'GROMACS/EE-fullequil':
            continue

        if isinstance(submission, YankSamplingAnalysis):
            submission_data = yank_analysis.get_free_energies_from_iteration(final_iteration=YANK_N_ITERATIONS)
        else:
            submission_data = submission.data
        submission_data = submission_data[submission_data['Simulation percentage'] == 100]

        row_str = submission.paper_name + '  &  '
        submission_final_DGs = []
        for system_id in submission_data['System ID'].unique():
            # GROMACS/EE doesn't have predictions for CB8-G3, and the
            # GROMACS/NS-DS/SB-long protocol was applied only to CB8-G3.
            if (('GROMACS/EE' in submission.paper_name and 'CB8-G3' in system_id) or
                        (submission.paper_name == 'GROMACS/NS-DS/SB-long' and 'OA' in system_id)):
                submission_final_DGs.append('')
                continue

            dg = submission_data.loc[submission_data['System ID'] == system_id, DG_KEY].values[0]
            ddg = submission_data.loc[submission_data['System ID'] == system_id, DDG_KEY].values[0]
            dg, ddg = reduce_to_first_significant_digit(dg, ddg)
            submission_final_DGs.append(r'{} $\pm$ {}'.format(dg, ddg))
        row_str += '  &  '.join(submission_final_DGs) + r'  \\'
        print(row_str)


# =============================================================================
# SUPPORTING INFORMATION - SINGLE TRAJECTORIES
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

        # REVO uses the mean of the 5 replicates to estimate the
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
    sns.set_style('whitegrid')
    sns.set_context('paper')

    if output_path_dir is None:
        output_path_dir = os.path.join(SAMPLING_PAPER_DIR_PATH, 'SI_Figure-individual-trajectories/')
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
        # REVO uses the mean of the 5 replicates to estimate the
        # uncertainty so it doesn't add information.
        if 'REVO' in submission.paper_name:
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
        output_file_name = 'replicates-{}-{}'.format(submission.file_name, submission.receipt_id)
        plt.savefig(os.path.join(output_path_dir, output_file_name + '.pdf'))
        # plt.savefig(os.path.join(output_path_dir, output_file_name + '.png'), dpi=300)
        # plt.show()


# =============================================================================
# SUPPORTING INFORMATION - HREX/MBAR STATISTICAL INEFFICIENCY ANALYSIS
# =============================================================================

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
    for stat_ineff in ['5', '10', '20', '50', '100', '200']:
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
            stat_ineff_axes[0].set_ylabel('$\Delta$G [kcal/mol]')
            stat_ineff_axes[1].set_ylabel('std($\Delta$G) [kcal/mol]')

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
        bbox_to_anchor = (0.05, 1.35)
        axes[0][0].legend(handles=handles[1:], labels=labels, loc='upper left',
                          bbox_to_anchor=bbox_to_anchor, ncol=6, fancybox=True,
                          labelspacing=0.8, handletextpad=0.5, columnspacing=1.2)

        # Save figure.
        output_file_name = 'statineff-{}'.format(system_name)
        plt.savefig(os.path.join(output_path_dir, output_file_name + '.pdf'))
        # plt.savefig(os.path.join(output_path_dir, output_file_name + '.png'), dpi=300)
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

    # Export YANK analysis and submissions to CSV/JSON tables.
    yank_analysis.export(os.path.join(SAMPLING_DATA_DIR_PATH, 'reference_free_energies'))
    for s in main_submissions:
        file_base_path = os.path.join(SAMPLING_DATA_DIR_PATH, s.receipt_id + '-reference')
        yank_analysis.export_by_submission(file_base_path, s)
    export_submissions(all_submissions, reference_free_energies)

    # Create example trajectory for the figure describing the challenge process.
    plot_example_bias_variance(yank_analysis, max_n_eval_percentage=0.4, mixed_proportion=0.3)

    # Cartoon explaining mean error and relative efficiency.
    plot_mean_error_cartoon()

    # Create figure with free energy, standard deviation, and bias as a function of computational cost.
    plot_all_entries_trajectory(main_submissions, yank_analysis, zoomed=False)
    plot_all_entries_trajectory(main_submissions, yank_analysis, zoomed=True)

    # Create results and efficiency table.
    print_relative_efficiency_table(main_submissions, yank_analysis, print_bias_corrected=False)

    # Plot nonequilibrium-switching single-direction estimator.
    plot_all_nonequilibrium_switching(noneq_submissions)

    # Plot sensitivity analysis figure.
    plot_restraint_and_barostat_analysis()

    # Plot figure for HREX bias analysis.
    plot_yank_bias()


    # Supporting information
    # ----------------------

    # Absolute/relative efficiency as a function of the computational cost.
    plot_relative_efficiencies(main_submissions, yank_analysis)
    plot_relative_efficiencies(main_submissions, yank_analysis, ci=None, same_plot=True)
    plot_absolute_efficiencies(main_submissions, yank_analysis)

    # Relative efficiency for uni/bi-directional estimators.
    print_nonequilibrium_relative_efficiencies(noneq_submissions)

    # Plot replicate predictions table.
    print_final_prediction_table(all_submissions, yank_analysis)

    # Plot individual trajectories.
    plot_all_single_trajectories_figures(all_submissions, yank_analysis)

    # Plot statistical inefficiency analysis.
    plot_hrex_stat_ineff_trajectories()

    # Supporting information for bias section.
    output_dir_path = os.path.join(SAMPLING_PAPER_DIR_PATH, 'SI_Figure-bias_hrex')
    plot_decomposition('CB8-G3', starting_iteration=5, type='phase',
                       output_file_path=output_dir_path + '/free-energy-phase-decomposition.pdf'))
    plot_decomposition('CB8-G3', starting_iteration=5, type='entropy-enthalpy',
                       output_file_path=output_dir_path + '/free-energy-entropy-decomposition.pdf')
