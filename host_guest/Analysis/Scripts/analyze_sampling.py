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

# All system ids.
SYSTEM_IDS = [
    'CB8-G3-0', 'CB8-G3-1', 'CB8-G3-2', 'CB8-G3-3', 'CB8-G3-4',
    'OA-G3-0', 'OA-G3-1', 'OA-G3-2', 'OA-G3-3', 'OA-G3-4',
    'OA-G6-0', 'OA-G6-1', 'OA-G6-2', 'OA-G6-3', 'OA-G6-4'
]

# Kelly's colors for maximum contrast.
#               "gray95",  "gray13",  "gold2",   "plum4",   "darkorange1", "lightskyblue2", "firebrick", "burlywood3", "gray51", "springgreen4", "lightpink2", "deepskyblue4", "lightsalmon2", "mediumpurple4", "orange", "maroon", "yellow3", "brown4", "yellow4", "sienna4", "chocolate", "gray19"
KELLY_COLORS = ['#F2F3F4', '#222222', '#F3C300', '#875692', '#F38400',     '#A1CAF1',       '#BE0032',  '#C2B280',   '#848482', '#008856',      '#E68FAC',    '#0067A5',     '#F99379',      '#604E97',      '#F6A600', '#B3446C', '#DCD300', '#882D17', '#8DB600', '#654522', '#E25822', '#2B3D26']
# Index of Kelly's colors associated to each submission.
SUBMISSION_COLORS = {
    'AMBER/APR': 11,
    'OpenMM/WExplore': 7,
    'OpenMM/SOMD': 4,
    'GROMACS/EE': 3,
    'GROMACS/EE-fullequil': 10,
    YANK_METHOD_PAPER_NAME: 9,
    'GROMACS/CT-NS-long': 6,
    'GROMACS/CT-NS': 1
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

    # Plot the mean free energy trajectory.
    ax.plot(x, mean_dg, color=color_mean, alpha=0.8, zorder=zorder, **plot_kwargs)
    # Plot mean trajectory confidence intervals.
    if plot_ci:
        ax.fill_between(x, mean_dg + sem_dg, mean_dg - sem_dg, alpha=0.3, color=color_ci, zorder=zorder)
    return ax


def plot_mean_data(mean_data, axes, color, label, x='N energy evaluations',
                   zorder=None, plot_bias=True, plot_ci=True):
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


def plot_all_entries_trajectory(submissions, yank_analysis, zoomed=False):
    """Plot the free energy trajectory of all the entries with CIs and first-hitting time."""
    system_names = ['CB8-G3', 'OA-G3', 'OA-G6']
    n_systems = len(system_names)

    # Create a figure with 3 columns (one for each system) and 2 rows.
    # The first row contains the free energy trajectory and CI, the second
    # a plot of the estimator variance, and the third the bias to the
    # asymptotic value.
    if zoomed:
        figsize = (7.25, 6.2)  # Without WExplorer
    else:
        figsize = (7.25, 8)  # With WExplorer
    fig, axes = plt.subplots(nrows=3, ncols=n_systems, figsize=figsize)

    # Fix axes limits based on maximum number of energy evaluations among the submissions.
    max_n_energy_evaluations = {system_name: 0 for system_name in system_names}

    # Optionally, remove WExplore.
    if zoomed:
        submissions = [s for s in submissions if s.name not in ['WExploreRateRatio']]

    for submission_idx, submission in enumerate(submissions):
        mean_free_energies = submission.mean_free_energies()
        unique_system_names = submission.data['System name'].unique()
        submission_color = KELLY_COLORS[SUBMISSION_COLORS[submission.paper_name]]

        # Plot free energy trajectories.
        for ax_idx, system_name in enumerate(unique_system_names):
            # CB8-G3 calculations for GROMACS/EE did not converge yet.
            if submission.name == 'Expanded-ensemble/MBAR' and system_name == 'CB8-G3':
                continue
            # The OA prediction of the NS short protocol are the same of the long protocol submission file.
            if submission.paper_name == 'GROMACS/CT-NS-long' and system_name != 'CB8-G3':
                # Just add the label.
                axes[0][ax_idx].plot([], color=submission_color, label=submission.paper_name)
                continue

            # Select the data for only this host-guest system.
            mean_data = mean_free_energies[mean_free_energies['System name'] == system_name]

            # Update maximum number of energy evaluations.
            n_energy_evaluations = mean_data['N energy evaluations'].values[-1]
            max_n_energy_evaluations[system_name] = max(max_n_energy_evaluations[system_name],
                                                        n_energy_evaluations)

            # TODO: Set zorder according to ranking of total variance so
            # TODO:     that smaller shaded area are on top of bigger shaded areas.
            # TODO:     Keep zorder=0 for the black color though as it sticks out.
            plot_mean_data(mean_data, axes[:,ax_idx], color=submission_color,
                           zorder=1, label=submission.paper_name)

            # Plot adding full cost of Wang-Landau equilibration.
            if 'Expanded' in submission.name:
                first_nonzero_idx = np.nonzero(mean_data[DG_KEY].values)[0][0]
                calibration_cost = mean_data['N energy evaluations'].values[first_nonzero_idx] * 4
                mean_data['N energy evaluations'] += calibration_cost
                label = submission.paper_name + '-fullequil'
                plot_mean_data(mean_data, axes[:,ax_idx], color=KELLY_COLORS[SUBMISSION_COLORS[label]],
                               zorder=1, label=label)

    # Plot the YANK mean free energies.
    for ax_idx, system_name in enumerate(system_names):
        label = YANK_METHOD_PAPER_NAME
        submission_color = KELLY_COLORS[SUBMISSION_COLORS[label]]
        yank_mean_data = yank_analysis.get_free_energies_from_energy_evaluations(max_n_energy_evaluations[system_name],
                                                                                 system_name=system_name,
                                                                                 mean_trajectory=True)
        # Black is much more visible than the other colors so let's push it on the back.
        plot_mean_data(yank_mean_data, axes[:,ax_idx], color=submission_color, zorder=0, label=label)

    # Fix labels.
    axes[0][0].set_ylabel('$\Delta$G [kcal/mol]')
    axes[1][0].set_ylabel('std($\Delta$G) [kcal/mol]')
    axes[2][0].set_ylabel('bias [kcal/mol]')
    axes[-1][1].set_xlabel('number of energy/force evaluations [10$^6$]')

    if zoomed:
        # Y-axis limits when WExplore calculations are excluded.
        y_limits = [
            [(-16, -9), (-8, -4), (-8, -4)],
            [(0, 2), (0, 1), (0, 1)],
            [(-3, 1), (-0.6, 0.6), (-0.6, 0.6)],
        ]
    else:
        # Y-axis limits when WExplore calculations are included.
        y_limits = [
            [(-17, -9), (-12.5, -5), (-12.5, -5)],
            [(0, 2), (0, 1.75), (0, 1.75)],
            [(-4, 4), (-0.6, 0.6), (-0.6, 0.6)],
        ]
    # Fix axes limits.
    for ax_idx, system_name in enumerate(system_names):
        for row_idx in range(len(axes)):
            ax = axes[row_idx][ax_idx]
            # Set the x-axis limits.
            ax.set_xlim((0, max_n_energy_evaluations[system_name]/N_ENERGY_EVALUATIONS_SCALE))
            # Keep the x-axis label only at the bottom row.
            if row_idx != 2:
                ax.xaxis.set_ticklabels([])
            y_lim = y_limits[row_idx][ax_idx]
            if y_lim is not None:
                ax.set_ylim(y_lim)

        # Set the system name in the title.
        axes[0][ax_idx].set_title(system_name)

    # Create a bias axis AFTER the ylim has been set.
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

    # Show/save figure.
    if zoomed:
        plt.tight_layout(h_pad=0.0, rect=[0.0, 0.00, 1.0, 0.92], w_pad=0.0)  # Without WExplorer
    else:
        plt.tight_layout(h_pad=0.0, rect=[0.0, 0.00, 1.0, 0.92])  # With WExplorer

    # Plot legend.
    if zoomed:
        bbox_to_anchor = (2.52, 1.55)  # Without WExplorer.
    else:
        bbox_to_anchor = (2.62, 1.48)  # With WExplorer.
    axes[0][1].legend(loc='upper right', bbox_to_anchor=bbox_to_anchor,
                      fancybox=True, ncol=4)#int(len(submissions)/2)+1)
    plt.subplots_adjust(wspace=0.35)
    # plt.show()
    if zoomed:
        file_name = 'Figure2-free_energy_trajectories_zoomed.pdf'
    else:
        file_name = 'Figure2-free_energy_trajectories.pdf'
    output_file_path = os.path.join('../SAMPLing/PaperImages', file_name)
    plt.savefig(output_file_path)


def plot_yank_system_bias(system_name, data_dir_paths, axes, shift_to_origin=True):
    """Plot the YANK free energy trajectoies when discarding initial samples for a single system."""
    color_palette = sns.color_palette('viridis', n_colors=len(data_dir_paths)+1)

    # Plot trajectories with truncated data.
    all_iterations = set()
    for data_idx, data_dir_path in enumerate(data_dir_paths):
        yank_analysis = YankSamplingAnalysis(data_dir_path)

        # In the YankAnalysis folder, each analysis starting from
        # iteration N is in the folder "BiasAnalysis-N/".
        label = data_dir_path.split('-')[-1][:-1]
        # First color is for the full data.
        color = color_palette[data_idx+1]

        # Collect all iterations that we'll plot for the full data.
        mean_data = yank_analysis.get_system_free_energies(system_name, mean_trajectory=True)
        all_iterations.update(mean_data['HREX iteration'].values.tolist())

        # Simulate plotting starting from the origin.
        if shift_to_origin:
            mean_data['HREX iteration'] -= mean_data['HREX iteration'].values[0]

        plot_mean_data(mean_data, axes, x='HREX iteration', color=color,
                       label=label, plot_bias=False, plot_ci=False)

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
                   label='0', plot_bias=False, plot_ci=False)
    axes[0].set_title(system_name)


def plot_yank_bias():
    """Plot YANK free energy trajectories when discarding initial samples."""
    system_names = ['CB8-G3', 'OA-G3', 'OA-G6']
    n_rows = 2
    n_cols = len(system_names)
    fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(7.25, 4.6))

    # The loops are based on a two dimensional array of axes.
    if n_rows == 1:
        axes = np.array([axes])

    # Sort paths by how many samples they have.
    data_dir_paths = ['YankAnalysis/BiasAnalysis-{}/'.format(i) for i in [1000, 2000, 4000, 8000, 16000]]

    # In the first column, plot the "unshifted" trajectory of CB8-G3,
    # with all sub-trajectories shifted to the origin.
    plot_yank_system_bias('CB8-G3', data_dir_paths, axes[:,0], shift_to_origin=False)
    # In the second and third columns, plot the trajectories of CB8-G3
    # and OA-G3 with all sub-trajectories shifted to the origin.
    plot_yank_system_bias('CB8-G3', data_dir_paths, axes[:,1], shift_to_origin=True)
    axes[0,0].set_title('CB8-G3 (shifted)')
    plot_yank_system_bias('OA-G3', data_dir_paths, axes[:,2], shift_to_origin=True)
    axes[0,0].set_title('OA-G3 (shifted)')

    # Fix axes limits and labels.
    ylimits = {
        'CB8-G3': (-12.5, -10.5),
        'OA-G3': (-8, -6.5),
        'OA-G6': (-8, -6.5)
    }
    for ax_idx, system_name in zip(range(3), ['CB8-G3', 'CB8-G3','OA-G3']):
        axes[0][ax_idx].set_ylim(ylimits[system_name])
    for ax_idx in range(3):
        axes[1][ax_idx].set_ylim((0, 0.3))

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
    axes[1][0].set_ylabel('std($\Delta$G) [kcal/mol]')
    axes[-1][1].set_xlabel('HREX iteration')

    plt.tight_layout(h_pad=0.1, rect=[0.0, 0.00, 1.0, 0.92])

    handles, labels = axes[0][0].get_legend_handles_labels()
    handles = [handles[-1]] + handles[:-1]
    labels = [labels[-1]] + labels[:-1]
    bbox_to_anchor = (0.2, 1.45)
    axes[0][0].legend(handles, labels, loc='upper left', bbox_to_anchor=bbox_to_anchor,
                      title='n discarded iterations', ncol=len(data_dir_paths)+1, fancybox=True)

    # plt.show()
    output_file_path = os.path.join('../SAMPLing/PaperImages', 'Figure3-bias_hrex.pdf')
    plt.savefig(output_file_path)


def compute_geometric_mean_relative_efficiencies(mean_data, reference_mean_data):
    """Compute the relative std, absolute bias, and RMSD efficiency for the data."""
    # Discard the initial frames of WExplorer and GROMACS/EE without a prediction.
    first_nonzero_idx = np.nonzero(mean_data[DG_KEY].values)[0][0]
    var = mean_data['std'].values[first_nonzero_idx:]**2
    reference_var = reference_mean_data['std'].values[first_nonzero_idx:]**2
    bias = mean_data['bias'].values[first_nonzero_idx:]
    reference_bias = reference_mean_data['bias'].values[first_nonzero_idx:]

    var_efficiencies = reference_var / var
    msd_efficiencies = (reference_var + reference_bias**2) / (var + bias**2)

    # For the bias, discard all zero elements that cause the geometric mean to be 0.
    # for b in [bias, reference_bias]:
    nonzero_indices = ~np.isclose(bias, 0.0)
    bias = bias[nonzero_indices]
    reference_bias = reference_bias[nonzero_indices]

    nonzero_indices = ~np.isclose(reference_bias, 0.0)
    bias = bias[nonzero_indices]
    reference_bias = reference_bias[nonzero_indices]

    abs_bias_efficiencies = np.abs(reference_bias) / np.abs(bias)

    # Scale the energy evaluations to make the numerical integration more stable.
    mean_std_efficiency = np.sqrt(sp.stats.mstats.gmean(var_efficiencies))
    mean_rmsd_efficiency = np.sqrt(sp.stats.mstats.gmean(msd_efficiencies))
    mean_abs_bias_efficiency = sp.stats.mstats.gmean(abs_bias_efficiencies)

    return mean_std_efficiency, mean_abs_bias_efficiency, mean_rmsd_efficiency


def compute_all_mean_relative_efficiencies(mean_data, reference_mean_data):
    # Compute the total std, bias and RMSD of the submission.
    relative_efficiencies = np.array(compute_geometric_mean_relative_efficiencies(mean_data, reference_mean_data))
    # Compute reference total statistics assuming that the reference calculation has converged.
    reference_mean_data = copy.deepcopy(reference_mean_data)
    reference_mean_data['bias'] -= reference_mean_data['bias'].values[-1]
    corrected_relative_efficiencies =  np.array(compute_geometric_mean_relative_efficiencies(mean_data, reference_mean_data))
    return relative_efficiencies, corrected_relative_efficiencies


def print_relative_efficiency_table(submissions, yank_analysis):
    """Create a table with total standard deviation, absolute bias, and error."""
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
                relative_efficiencies, corrected_relative_efficiencies = np.full((2, 3), fill_value=np.nan)
                dg = ''
                n_force_eval = ''
            else:
                # Select the data for only this host-guest system.
                mean_data = mean_free_energies[mean_free_energies['System name'] == system_name]
                # Select the corresponding iterations for the YANK calculation.
                n_energy_evaluations = mean_data['N energy evaluations'].values[-1]
                reference_mean_data = yank_analysis.get_free_energies_from_energy_evaluations(
                    n_energy_evaluations, system_name=system_name, mean_trajectory=True)
                # Compute relative std, bias, and RMSD efficiencies.
                relative_efficiencies, corrected_relative_efficiencies = compute_all_mean_relative_efficiencies(mean_data, reference_mean_data)

                # Get the final free energy and number of energy/force evaluations.
                dg = mean_data[DG_KEY].values[-1]
                dg_CI = mean_data['$\Delta$G CI'].values[-1]  # Confidence interval.
                dg, dg_CI = reduce_to_first_significant_digit(dg, dg_CI)
                n_force_eval = mean_data['N energy evaluations'].values[-1]
                # Convert to string format.
                dg = '{} $\\pm$ {}'.format(dg, dg_CI)
                n_force_eval = str(int(round(n_force_eval / 1e6)))

            # Add free energy and cost entries.
            efficiency_table[(system_name, column_names[0])].append(dg)
            efficiency_table[(system_name, column_names[1])].append(n_force_eval)

            # Add efficiency entries for the table.
            for statistic_idx, statistic_name in enumerate(statistic_names):
                relative_efficiency = relative_efficiencies[statistic_idx]
                corrected_relative_efficiency = corrected_relative_efficiencies[statistic_idx]
                # Print significant digits.
                efficiencies_format = []
                for e in [relative_efficiency, corrected_relative_efficiency]:
                    efficiencies_format.append('{:.2f}' if e < 0.09 else '{:.1f}')

                if np.isnan(relative_efficiency):
                    data_entry = ''
                elif 'std' not in statistic_name:
                    data_entry = efficiencies_format[0] + ' (' + efficiencies_format[1] + ')'
                    data_entry = data_entry.format(relative_efficiency, corrected_relative_efficiency)
                else:
                    # Standard deviation efficiency is not affected by the bias.
                    data_entry = efficiencies_format[0].format(relative_efficiency)
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
        efficiency_table[(system_name, statistic_name)].append('1.0')

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
    palette_mean = sns.color_palette('dark')
    submission_mean_color = 8
    reference_mean_color = 9

    # Plot the 5 replicates individual trajectories.
    # First remove the initial predictions that are 0.0 (i.e. there is no estimate).
    ax = axes[0]
    system_data = system_data[system_data[DG_KEY] != 0.0]
    sns.lineplot(data=system_data, x='N energy evaluations', y=DG_KEY,
                 hue='System ID', palette='bright', ax=ax, zorder=5, alpha=0.6)

    # Plot the submission mean trajectory with CI.
    plot_mean_free_energy(system_mean_data, x='N energy evaluations',  ax=ax,
                          color_mean=palette_mean[submission_mean_color],
                          color_ci=palette_mean[submission_mean_color], label='Mean $\Delta$G',
                          scale_n_energy_evaluations=False)

    # Plot YANK mean trajectory with CI.
    if reference_system_mean_data is not None:
        plot_mean_free_energy(reference_system_mean_data, x='N energy evaluations', ax=ax,
                              color_mean=palette_mean[reference_mean_color],
                              color_ci=palette_mean[reference_mean_color], label='Ref mean $\Delta$G',
                              scale_n_energy_evaluations=False)

    ax.set_title(system_name)
    # Add the y-label only on the leftmost Axis.
    if system_name != 'CB8-G3':
        ax.set_ylabel('')
    # Remove the legend for now, which will be added at the end after tighting up the plot.
    ax.get_legend().remove()

    # Create a bias axis.
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

        # Plot the standard deviation of the free energy trajectories.
        submission_std = system_mean_data['std']
        cost = system_mean_data['N energy evaluations'].values / N_ENERGY_EVALUATIONS_SCALE
        # cost = system_mean_data['Simulation percentage'].values
        ax.plot(cost, submission_std, color=palette_mean[submission_mean_color])
        if reference_system_mean_data is not None:
            reference_std = reference_system_mean_data['std']
            ax.plot(cost, reference_std, color=palette_mean[reference_mean_color])

        # WExplore uses the mean of the 5 replicates to estimate the
        # uncertainty so it doens't add information.
        if plot_methods_uncertainties:
            # Plot the method uncertainties of the single replicate trajectories.
            # First scale the number of energy evaluations.
            system_data.loc[:,'N energy evaluations'] /= N_ENERGY_EVALUATIONS_SCALE
            sns.lineplot(data=system_data, x='N energy evaluations', y=DDG_KEY,
                         hue='System ID', palette='bright', ax=ax, alpha=0.6)

            # The legend is added later at the top.
            ax.get_legend().remove()

        # Set y-axis limits.
        if reference_system_mean_data is not None:
            ax.set_ylim((0, max(max(submission_std), max(reference_std))))
        else:
            # Ignore the first data points.
            ax.set_ylim((0, max(submission_std)))
        # Only the central plot shows the x-label.
        ax.set_xlabel('')
        # Add the y-label only on the leftmost Axis.
        if system_name != 'CB8-G3':
            ax.set_ylabel('')
        else:
            ax.set_ylabel(DDG_KEY)

    # The x-label is shown only in the central plot.
    if system_name == 'OA-G3':
        ax.set_xlabel('N energy evaluations  [10$^6$]')


def plot_all_single_trajectories_figures(submissions, yank_analysis, plot_errors=True):
    """Individual plots for each method with the 5 individual free energy and uncertainty trajectories."""
    sns.set_context('paper')

    output_path_dir = '../SAMPLing/PaperImages/SI_Figure1-individual-trajectories/'
    os.makedirs(output_path_dir, exist_ok=True)

    # -------------------- #
    # Plot submission data #
    # -------------------- #

    for submission in submissions + [yank_analysis]:
        # CB8-G3 calculations for GROMACS/EE did not converge yet.
        if submission.name == 'Expanded-ensemble/MBAR':
            submission.data = submission.data[submission.data['System name'] != 'CB8-G3']
        # WExplore uses the mean of the 5 replicates to estimate the
        # uncertainty so it doens't add information.
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

        # for system_name in unique_system_names:
        for ax_idx, system_name in enumerate(unique_system_names):

            if isinstance(submission, YankSamplingAnalysis):
                data = submission.get_free_energies_from_iteration(final_iteration=YANK_N_ITERATIONS,
                                                                   system_name=system_name)
                mean_data = submission.get_free_energies_from_iteration(final_iteration=YANK_N_ITERATIONS,
                                                                        system_name=system_name,
                                                                        mean_trajectory=True)
                reference_mean_data = None
            else:
                # Select the data for only this host-guest system.
                data = submission.data[submission.data['System name'] == system_name]
                mean_data = mean_free_energies[mean_free_energies['System name'] == system_name]

                # Get the corresponding YANK free energies.
                # The cost for the same system is the same for all replicates.
                n_energy_evaluations = int(submission.cost.loc[system_name + '-0', 'N energy evaluations'])
                reference_mean_data = yank_analysis.get_free_energies_from_energy_evaluations(n_energy_evaluations,
                                                                                              system_name=system_name,
                                                                                              mean_trajectory=True)

            plot_single_trajectories_figures(axes[:,ax_idx], data, mean_data, plot_errors=plot_errors,
                                             reference_system_mean_data=reference_mean_data,
                                             plot_methods_uncertainties=plot_methods_uncertainties)

        plt.tight_layout(w_pad=0.7, rect=[0.0, 0.0, 1.0, 0.85])
        # Create legend.
        bbox_to_anchor = (0.35, 1.5)
        # The first handle/label is the legend title "System ID" so we get rid of it.
        handles, labels = trajectory_axes[0].get_legend_handles_labels()
        labels = ['replicate ' + str(i) for i in range(5)] + labels[6:]
        trajectory_axes[0].legend(handles=handles[1:], labels=labels, loc='upper left',
                                  bbox_to_anchor=bbox_to_anchor, ncol=4, fancybox=True)

        # Save figure.
        output_file_name = '{}-{}.pdf'.format(submission.receipt_id, submission.file_name)
        plt.savefig(os.path.join(output_path_dir, output_file_name))
        # plt.show()


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    sns.set_style('whitegrid')
    sns.set_context('paper')

    # Flag that controls whether to plot the trajectory of uncertainties and std.
    PLOT_ERRORS = True

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
    submissions = load_submissions(SamplingSubmission, SAMPLING_SUBMISSIONS_DIR_PATH, user_map)
    # Remove AMBER/TI.
    submissions = [s for s in submissions if s.name not in ['Langevin/Virtual Bond/TI']]

    # # Create an extra submission for GROMACS/EE where the full cost of equilibration has been taken into account.
    # gromacs_ee_submission = copy.deepcopy([s for s in submissions if s.paper_name == 'GROMACS/EE'][0])
    # gromacs_ee_submission.paper_name = 'GROMACS/EE-fullequil'
    # data = gromacs_ee_submission.data  # Shortcut.
    # mean_free_energies = gromacs_ee_submission.mean_free_energies()
    # for system_name in ['OA-G3', 'OA-G6']:
    #     mean_data = mean_free_energies[mean_free_energies['System name'] == system_name]
    #     first_nonzero_idx = np.nonzero(mean_data[DG_KEY].values)[0][0]
    #     full_equilibration_cost = mean_data['N energy evaluations'].values[first_nonzero_idx] * 4
    #     for i in (data['System name'] == system_name).index:
    #         data.at[i, 'N energy evaluations'] += full_equilibration_cost
    # submissions.append(gromacs_ee_submission)

    # Sort the submissions to have all pot and tables in the same order.
    submissions = sorted(submissions, key=lambda s: s.paper_name)

    # # Export YANK analysis and submissions to CSV/JSON tables.
    # yank_analysis.export(os.path.join(SAMPLING_DATA_DIR_PATH, 'reference_free_energies'))
    # export_submissions(submissions, reference_free_energies)

    # Create figure with free energy, standard deviation, and bias as a function of computational cost.
    plot_all_entries_trajectory(submissions, yank_analysis, zoomed=False)

    # Create results and efficiency table.
    # print_relative_efficiency_table(submissions, yank_analysis)

    # Plot figure for HREX bias analysis.
    # plot_yank_bias()

    # Plot individual trajectories.
    # plot_all_single_trajectories_figures(submissions, yank_analysis)
    import sys; sys.exit()

    # # TODO REMOVE ME: CODE FOR SINGLE PLOT FOR SLIDES
    # # Plot submission data.
    # os.makedirs(SAMPLING_PLOT_DIR_PATH, exist_ok=True)
    # palette_mean = sns.color_palette('dark')
    # for submission in submissions:
    #     if 'SOMD' not in submission.name:
    #         continue
    #
    #     mean_free_energies = submission.mean_free_energies()
    #     for system_name in submission.data['System name'].unique():
    #         # Select the data for only this host-guest system.
    #         data = submission.data[submission.data['System name'] == system_name]
    #         mean_data = mean_free_energies[mean_free_energies['System name'] == system_name]
    #
    #         # Get the corresponding YANK free energies.
    #         # The cost for the same system is the same for all replicates.
    #         n_energy_evaluations = int(submission.cost.loc[system_name + '-0', 'N energy evaluations'])
    #         yank_mean_data = yank_analysis.get_free_energies_from_energy_evaluations(n_energy_evaluations,
    #                                                                              system_name=system_name,
    #                                                                              mean_trajectory=True)
    #
    #         fig, ax = plt.subplots(figsize=(7.5, 6.5))
    #
    #         # Plot the 5 replicates trajectories.
    #         sns.tsplot(data=data, time='Simulation percentage', value=DG_KEY,
    #                    unit='System name', condition='System ID', color='pastel', ax=ax)
    #
    #         # Plot the submission mean trajetory with CI.
    #         plot_mean_free_energy(mean_data, ax=ax, color_mean=palette_mean[0],
    #                               label='Mean $\Delta$G')
    #
    #         # Plot YANK mean trajectory with CI.
    #         plot_mean_free_energy(yank_mean_data, ax=ax, color_mean=palette_mean[2],
    #                               label='Ref mean $\Delta$G')
    #
    #         # Create a bias axis.
    #         ref_free_energy = reference_free_energies.loc[system_name, DG_KEY]
    #         with sns.axes_style('white'):
    #             ax2 = ax.twinx()
    #             # Plot a vertical line to make the scale.
    #             vertical_line = np.linspace(*ax.get_ylim()) - ref_free_energy
    #             ax2.plot([50] * len(vertical_line), vertical_line, alpha=0.0001)
    #             ax2.grid(alpha=0.5, linestyle='dashed', zorder=0)
    #             # We add the bias y-label only on the rightmost Axis.
    #             ax2.set_ylabel('Bias to reference [kcal/mol]')
    #
    #         # Set axis limits/titles.
    #         ax.set_ylim((-20, 4))
    #         ax.set_title('{} - {} ({})'.format(system_name, submission.name, submission.receipt_id))
    #         ax.set_xlabel(ax.get_xlabel() + ' (N energy evaluations: {:,})'.format(n_energy_evaluations))
    #         ax.legend(ncol=2)
    #         plt.tight_layout()
    #         output_path_dir = os.path.join(SAMPLING_PLOT_DIR_PATH,
    #                                        '{}-{}.pdf'.format(submission.receipt_id, system_name))
    #         # plt.show()
    #         plt.savefig(output_path_dir)
    # # TODO END REMOVE ME: CODE FOR SINGLE PLOT FOR SLIDES

    # # TODO COMMENT FOR SEPARATE PLOTS
    # n_systems = 3
    # if PLOT_ERRORS:
    #     raise NotImplementedError('Plot errors of all 3 systems in single figure is not supported.')
    # else:
    #     fig, axes = plt.subplots(nrows=1, ncols=n_systems, figsize=(6*n_systems, 6))
    #
    # # Plot YANK replicates isolated.
    # # for system_name in reference_free_energies.index:
    # system_names = ['CB8-G3', 'OA-G3', 'OA-G6']
    # for ax_idx, (system_name, ax) in enumerate(zip(system_names[:n_systems], axes)):
    #     # Select the data for only this host-guest system.
    #     yank_data = yank_analysis.get_system_free_energies(system_name)
    #     yank_mean_data = yank_analysis.get_system_free_energies(system_name, mean_trajectory=True)
    #
    #     # # TODO DECOMMENT FOR SEPARATE PLOTS
    #     # if PLOT_ERRORS:
    #     #     fig, axes = plt.subplots(nrows=2, figsize=(7, 12), sharex=True)
    #     #     ax = axes[0]
    #     # else:
    #     #     fig, ax = plt.subplots(nrows=1, figsize=(7, 6))
    #
    #     # Plot the 5 replicates trajectories.
    #     sns.tsplot(data=yank_data, time='N energy evaluations', value=DG_KEY,
    #                unit='System name', condition='System ID', color='pastel', ax=ax)
    #
    #     # Plot the submission mean trajectory with CI.
    #     plot_mean_free_energy(yank_mean_data, ax=ax, color_mean=palette_mean[0],
    #                           x='N energy evaluations', label='Mean $\Delta$G')
    #
    #     # Plot EE NPT values.
    #     ref2_dg = None
    #     # if system_name == 'OA-G3':
    #     #     ref2_dg = -6.0057999999999989
    #     #     ref2_ci = 0.21634797116823681
    #     # elif system_name == 'OA-G6':
    #     #     ref2_dg = -6.8739999999999997
    #     #     ref2_ci = 0.28182900923056731
    #     # else:
    #     #     ref2_dg = None
    #     if ref2_dg is not None:
    #         n_energy_evaluations = yank_mean_data['N energy evaluations'].values
    #         ax.plot(n_energy_evaluations, [ref2_dg for _ in n_energy_evaluations],
    #                 color=palette_mean[2], label='Ref2 mean $\Delta$G')
    #         ax.fill_between(n_energy_evaluations, ref2_dg + ref2_ci, ref2_dg - ref2_ci, alpha=0.65)
    #
    #     # Set axis limits/titles.
    #     ax.set_ylim((-13, -6))
    #     # Set axis title.
    #     # TODO DECOMMENT TITLE FOR SINGLE PLOT
    #     # ax.set_title('{} - Reference'.format(system_name))
    #     ax.legend(loc='lower right')
    #
    #     # TODO REMOVE AXIS IDX CONDITION FOR SINGLE PLOT
    #     if ax_idx != 0:
    #         ax.set_ylabel('')
    #         ax.yaxis.set_ticklabels([])
    #     else:
    #         ax.set_ylabel('Reference\n' + ax.get_ylabel())
    #
    #     # Plot uncertainties and standard deviation of the 5 trajectories.
    #     if PLOT_ERRORS:
    #         # Remove x label of old axis.
    #         ax.set_xlabel('')
    #         # Plot the standard deviation and trajectory uncertainties.
    #         ax = axes[1]
    #         sns.tsplot(data=yank_data, time='N energy evaluations', value='d$\Delta$G [kcal/mol]',
    #                    unit='System name', condition='System ID', color='pastel', ax=ax)
    #         ax.plot(yank_mean_data['N energy evaluations'].values, yank_mean_data['std'].values, label='OA-G3 std')
    #         ax.legend()
    #
    #     ax.legend(ncol=2)
    #
    # # TODO INDENT AND SWITCH OUTPUT PATH FOR SEPARATE PLOTS
    # plt.tight_layout(pad=0.1)
    # # plt.show()
    # # output_path_dir = os.path.join(SAMPLING_PLOT_DIR_PATH,
    # #                                'reference-{}.pdf'.format(system_name))
    # output_path_dir = os.path.join(SAMPLING_PLOT_DIR_PATH,
    #                                'reference.pdf'.format(system_name))
    # plt.savefig(output_path_dir)


    # =============================================================================
    # FIGURES GENERATED FOR THE PAPER
    # =============================================================================

    # sns.set_context('paper')
    #
    # final_free_energies = {
    #     'YANK NPT': {
    #         'OA-G3': [-6.693039053, -6.712935302, -6.718535019, -6.67622114, -6.716658183],
    #         'OA-G6': [-7.142365148, -7.174144518, -7.182525502, -7.137949615, -7.246234684]
    #     },
    #     'YANK NVT production PME': {
    #         'OA-G3': [-6.765, -6.652, -6.695, -6.875, -6.682],
    #         'OA-G6': [-7.132, -7.284, -7.093, -7.109, -7.114]
    #     },
    #     'YANK NVT matched PME': {
    #         'OA-G3': [-6.639, -6.700, -6.556, -6.667, -6.654],
    #         'OA-G6': [-7.133, -7.072, -7.245, -7.007, -7.241]
    #     },
    #     'EE NPT': {
    #         'OA-G3': [-5.950, -5.841, -6.183, -5.856, -6.199],
    #         'OA-G6': [-6.693, -6.712, -7.254, -6.887, -6.824]
    #     },
    #     'EE NVT production PME': {
    #         'OA-G3': [-6.440, -6.767, -6.591, -6.522, -6.556],
    #         'OA-G6': [-6.962, -6.978, -6.974, -7.203, -6.796]
    #     },
    #     'EE NVT matched PME': {
    #         'OA-G3': [-6.645, -6.723, -6.484, -6.611, -6.583],
    #         'OA-G6': [-6.676, -6.909, -6.971, -6.955, -7.052]
    #     }
    # }
    #
    # # Convert to Pandas Dataframe.
    # replicate_free_energies = []
    # for calculation_name, calculation_data in final_free_energies.items():
    #     method, ensemble = calculation_name.split(' ', 1)
    #     for system_id, data in calculation_data.items():
    #         # dg, ddg = mean_confidence_interval(data)
    #         for data_point in data:
    #             replicate_free_energies.append({
    #                 'method': method,
    #                 'system': '{} - {}'.format(system_id, ensemble),
    #                 '-'+DG_KEY: - data_point,
    #                 # '$\Delta$G CI': ddg
    #             })
    # # Reorder bars.
    # replicate_free_energies = sorted(replicate_free_energies, key=lambda x: x['system'])
    # replicate_free_energies = pd.DataFrame(replicate_free_energies)

    # # Compute average free energies and CI.
    # from pkganalysis.stats import mean_confidence_interval
    # mean_free_energies = {calculation_name: {system_id: mean_confidence_interval(data)
    #                                          for system_id, data in calculation_data.items()}
    #                       for calculation_name, calculation_data in final_free_energies.items()}
    # from pprint import pprint
    # pprint(mean_free_energies)


    # BAR PLOT WITH NPT AND NVT RESULTS
    # ====================================

    # # Barplot.
    # ax = sns.barplot(data=replicate_free_energies, x='system', y='-'+DG_KEY, hue='method', ci='sd')
    # ax.legend(loc='upper left')
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=40)
    # plt.tight_layout()
    # # plt.show()
    # plt.savefig('NVT_NPT_free_energies.png', dpi=300)




    # =============================================================================
    # FIGURES GENERATED FOR THE TALK
    # =============================================================================
    #
    # for submission in submissions:
    #     mean_free_energies = submission.mean_free_energies()
    #
    #     # fig, axes = plt.subplots(ncols=3)
    #     for i, system_name in enumerate(submission.data['System name'].unique()):
    #         # Select the data for only this host-guest system.
    #         data = submission.data[submission.data['System name'] == system_name]
    #         mean_data = mean_free_energies[mean_free_energies['System name'] == system_name]
    #
    #         # Get the corresponding YANK free energies.
    #         # The cost for the same system is the same for all replicates.
    #         n_energy_evaluations = int(submission.cost.loc[system_name + '-0', 'N energy evaluations'])
    #         yank_mean_data = yank_analysis.get_free_energies_from_energy_evaluations(n_energy_evaluations,
    #                                                                              system_name=system_name,
    #                                                                              mean_trajectory=True)
    #
    #         fig, ax = plt.subplots(figsize=(7, 6))
    #         # ax = axes[i]
    #
    #         # Plot the 5 replicates trajectories.
    #         sns.tsplot(data=data, time='Simulation percentage', value=DG_KEY,
    #                    unit='System name', condition='System ID', color='pastel', ax=ax)
    #
    #         # Plot the submission mean trajetory with CI.
    #         plot_mean_free_energy(mean_data, ax=ax, color_mean=palette_mean[0],
    #                               label='Mean $\Delta$G')
    #
    #         # Plot YANK mean trajectory with CI.
    #         plot_mean_free_energy(yank_mean_data, ax=ax, color_mean=palette_mean[2],
    #                               label='Ref mean $\Delta$G')
    #
    #         # Plot reference value.
    #         ref_free_energy = reference_free_energies.loc[system_name, DG_KEY]
    #         ax.plot(mean_data['Simulation percentage'], [ref_free_energy for _ in range(100)],
    #                 color=palette_mean[1], ls='--', label='Ref final $\Delta$G')
    #
    #         # Set axis limits/titles.
    #         ax.set_ylim((-20, 4))
    #         ax.set_title(system_name)
    #         ax.set_xlabel(ax.get_xlabel() + ' (N energy evaluations: {:,})'.format(n_energy_evaluations))
    #         ax.legend(ncol=2)
    #         plt.tight_layout()
    #         output_path_dir = os.path.join(SAMPLING_PLOT_DIR_PATH,
    #                                        '{}-{}.pdf'.format(submission.receipt_id, system_name))
    #         # plt.show()
    #         plt.savefig(output_path_dir)



    # # fig, axes = plt.subplots(ncols=3, sharey=True, figsize=(15, 5))
    # # for ax, system_name in zip(axes, reference_free_energies.index):
    # for i, system_name in enumerate(reference_free_energies.index):
    #     # Select the data for only this host-guest system.
    #     yank_data = yank_analysis.get_system_free_energies(system_name)
    #     yank_mean_data = yank_analysis.get_system_free_energies(system_name, mean_trajectory=True)
    #
    #     fig, ax = plt.subplots(figsize=(6, 6))
    #
    #     # Plot the 5 replicates trajectories.
    #     sns.tsplot(data=yank_data, time='N energy evaluations', value=DG_KEY,
    #                unit='System name', condition='System ID', color='pastel', ax=ax)
    #     #
    #     # if i != 0:
    #     #     ax.set_ylabel('')
    #     #     ax.set_yticklabels([])
    #     # elif i != 1:
    #     #     ax.set_xlabel('')
    #     if i != 0:
    #         ax.set_ylim((-10, -5))
    #
    #     # Plot the submission mean trajectory with CI.
    #     plot_mean_free_energy(yank_mean_data, ax=ax, color_mean=palette_mean[0],
    #                           x='N energy evaluations', label='Mean $\Delta$G')
    #
    #     # Print second reference value.
    #     # if i != 0:
    #     #     mean_dg, sem_dg = travis_mean[system_name]
    #     #     print('{}: {} \pm {}'.format(system_name, mean_dg, sem_dg))
    #     #     n_energy_evaluations = yank_mean_data['N energy evaluations'].values
    #     #     ax.plot(n_energy_evaluations, [mean_dg for _ in n_energy_evaluations],
    #     #             color=palette_mean[2], label='Ref2 mean $\Delta$G')
    #     #     ax.fill_between(n_energy_evaluations, mean_dg + sem_dg, mean_dg - sem_dg, alpha=0.65)
    #
    #     # Set axis limits/titles.
    #     ax.set_title('{} - Reference'.format(system_name))
    #     # ax.set_ylim((-15, -6))
    #     ax.legend(ncol=2)
    #     plt.tight_layout()
    #     output_path_dir = os.path.join(SAMPLING_PLOT_DIR_PATH,
    #                                    'reference-{}.pdf'.format(system_name))
    #     # plt.show()
    #     plt.savefig(output_path_dir)


