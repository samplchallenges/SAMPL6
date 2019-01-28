import copy
import glob
import itertools
import os

import numpy as np
import pymbar
import scipy.spatial


# =============================================================================
# GLOBALS
# =============================================================================

ANALYSIS_DIR_PATH = os.path.join('YankAnalysis', 'WaterAnalysis')
PLOTS_DIR_PATH = os.path.join('..', 'SAMPLing', 'Plots', 'WaterAnalysis')
PAPER_IMAGES_DIR_PATH = os.path.join('..', 'SAMPLing', 'PaperImages')

# SYSTEM_IDS = ['CB8-G3-0', 'CB8-G3-1', 'CB8-G3-2', 'CB8-G3-3', 'CB8-G3-4',
#               'OA-G3-0', 'OA-G3-1', 'OA-G3-2', 'OA-G3-3', 'OA-G3-4',
#               'OA-G6-0', 'OA-G6-1', 'OA-G6-2', 'OA-G6-3', 'OA-G6-4']
SYSTEM_IDS = ['CB8-G3-0']


BINDING_SITE_SELECTION = {
    'CB8': 'resname HST and mass > 1.5',

    # All the heavy atoms except for the the 4 propionic acids at the bottom and the 4 flaps at the top.
    # PyMol selection syntax only basket heavy atoms: select oa, (id 1-40) or (id 57-88) or id 90 or (id 93-95) or (id 98-99) or id 101 or (id 104-105) or id 123 or id 125 or id 128
    'OA': '(index 0 to 39) or (index 56 to 87) or index 89 or (index 92 to 94) or (index 97 to 98) or index 100 or (index 103 to 104) or index 122 or index 124 or index 127'
    # PyMol selection syntax basket+flaps heavy atoms: select oa, (id 1-40) or (id 57-111) or (id 115-121) or (id 123-128)
    # 'OA': '(index 0 to 39) or (index 56 to 110) or (index 114 to 120) or (index 122 to 127)'
}

SYSTEM_N_STATES = {
    'CB8-G3': 69,
    'OA-G3': 59,
    'OA-G6': 55
}


# =============================================================================
# DEBUG FUNCTIONS
# =============================================================================

def compute_n_bound_waters(nc_file_path, binding_site_selection, state_idx=None, replica_idx=None):
    from yank.analyze import extract_trajectory
    # Extract the trajectory of the state/replica.
    trajectory = extract_trajectory(nc_file_path, state_index=state_idx, replica_index=replica_idx,
                                    keep_solvent=True, discard_equilibration=False, image_molecules=True)

    # Map atom_index: molecule_index
    atom_to_molecule = {}
    molecule_to_atoms = trajectory.topology.find_molecules()
    for molecule_idx, molecule_atoms in enumerate(molecule_to_atoms):
        for atom in molecule_atoms:
            atom_to_molecule[atom.index] = molecule_idx

    # Compute number of waters in the binding site.
    binding_site_atom_indices = trajectory.topology.select(binding_site_selection)
    water_atom_indices = trajectory.topology.select('water')

    # Create PDB of binding site only to make sure the selection is correct.
    if state_idx == 0:
        binding_site_traj = trajectory.slice(0)
        binding_site_traj = binding_site_traj.atom_slice(binding_site_atom_indices)
        file_name = os.path.dirname(nc_file_path).split('/')[-1] + '-binding-site.pdb'
        binding_site_traj.save(os.path.join(ANALYSIS_DIR_PATH, file_name))

    n_waters = np.zeros(trajectory.n_frames, dtype=np.int64)
    for frame_idx, frame_positions in enumerate(trajectory.xyz):
        binding_site_positions = frame_positions[binding_site_atom_indices]
        waters_positions = frame_positions[water_atom_indices]

        # Obtain a Delaunay triangulation describing the binding site volume.
        triangulation = scipy.spatial.Delaunay(binding_site_positions)
        # Array of boolean, one for each atom.
        are_water_atoms_bound = triangulation.find_simplex(waters_positions) > 0

        # Convert is_water_atom_bound to is_water_molecule_bound.
        bound_water_molecules = set()
        for water_atom_idx, is_water_atom_bound in zip(water_atom_indices, are_water_atoms_bound):
            if is_water_atom_bound:
                bound_water_molecules.add(atom_to_molecule[water_atom_idx])
        n_waters[frame_idx] = len(bound_water_molecules)
    return n_waters


# =============================================================================
# EXTRACT DATA
# =============================================================================

def get_complex_nc_file_path(system_id):
    nc_file_path_pattern = os.path.join('..', 'SAMPLing_NPT_waters', 'experiments', 'experiment-{}', '{}', 'complex.nc')
    host_name, guest_name, replicate_idx = system_id.split('-')
    host_guest_name = host_name + guest_name
    return nc_file_path_pattern.format(host_guest_name, host_guest_name + replicate_idx)


def extract_replica_state_indices(system_id):
    from yank.repex import Reporter
    nc_file_path = get_complex_nc_file_path(system_id)
    reporter = Reporter(nc_file_path, open_mode='r')
    try:
        replica_state_indices = reporter.read_replica_thermodynamic_states()
    finally:
        reporter.close()
    output_file_path = os.path.join(ANALYSIS_DIR_PATH, '{}-replica-state-indices.npy'.format(system_id))
    np.save(output_file_path, replica_state_indices)


def compute_all_n_bound_waters(job_id):
    # Prepare all jobs arguments.
    jobs = []
    for system_id in SYSTEM_IDS:
        host_name, guest_name, replicate_idx = system_id.split('-')
        job = {
            'system_id': system_id,
            'nc_file_path': get_complex_nc_file_path(system_id),
            'binding_site_selection': BINDING_SITE_SELECTION[host_name]
        }
        n_states = SYSTEM_N_STATES[host_name + '-' + guest_name]
        for state_idx in range(n_states):
            job_state = [copy.deepcopy(job), copy.deepcopy(job)]
            job_state[0]['state_idx'] = state_idx
            job_state[1]['replica_idx'] = state_idx
            jobs.extend(job_state)

    if job_id is None:
        print('jobid must be specified. Total number of jobs is {}'.format(len(jobs)))
        return

    # Create output directory.
    job_kwargs = jobs[job_id]
    system_id = job_kwargs.pop('system_id')
    os.makedirs(os.path.join(ANALYSIS_DIR_PATH, system_id), exist_ok=True)
    try:
        output_file_name = 'replica{}.npy'.format(job_kwargs['replica_idx'])
    except KeyError:
        output_file_name = 'state{}.npy'.format(job_kwargs['state_idx'])

    # Compute and save the number of waters.
    n_bound_waters = compute_n_bound_waters(**job_kwargs)
    np.save(os.path.join(ANALYSIS_DIR_PATH, system_id, output_file_name), n_bound_waters)


# =============================================================================
# FIGURES
# =============================================================================

def get_n_states(system_id):
    return len(list(glob.glob(os.path.join(ANALYSIS_DIR_PATH, system_id, 'state*.npy'))))


def iterate_bound_waters(system_id, state=True):
    n_states = get_n_states(system_id)

    # Obtain all the saved n_water_bound trajectories for the system.
    trajectory_type = 'state' if state else 'replica'
    file_path_pattern = os.path.join(ANALYSIS_DIR_PATH, system_id, trajectory_type + '{}.npy')
    for state_idx in range(n_states):
        n_bound_waters_file_path = file_path_pattern.format(state_idx)
        temp = np.load(n_bound_waters_file_path)
        yield temp


def plot_n_bound_waters_trajectory(system_id):
    # Plot the bound waters trajectories in two columns (first for replicas, second for states)
    n_states = get_n_states(system_id)
    fig, axes = plt.subplots(nrows=n_states+1, ncols=2, figsize=(8, 1.2*n_states))

    # Obtain all the saved n_water_bound trajectories for the system.
    for col_idx, trajectory_type in enumerate(['replica', 'state']):
        state_flag = True if trajectory_type == 'state' else False
        all_n_bound_waters = []
        for state_idx, n_bound_waters in enumerate(iterate_bound_waters(system_id, state=state_flag)):
            ax = axes[state_idx+1][col_idx]
            ax.plot(n_bound_waters)
            # Add info about statistical inefficiency.
            g_t = pymbar.timeseries.statisticalInefficiency(n_bound_waters)
            ax.set_title('{} {}\nstatistical inefficiency {:.2f} ps'.format(trajectory_type, state_idx, g_t))
            ax.set_xlabel('iteration')
            ax.set_ylabel('n bound waters')
            ax.set_ylim((0, 15))

            all_n_bound_waters.append(n_bound_waters)

        # Plot super-replica correlation function.
        ax = axes[0][col_idx]
        g_t, C_t = pymbar.timeseries.statisticalInefficiencyMultiple(all_n_bound_waters, fast=False,
                                                                     return_correlation_function=True)
        t = [x[0] for x in C_t]
        C_t = [x[1] for x in C_t]
        ax.plot(t, C_t)
        ax.set_title('statistical inefficiency {:.2f} ps'.format(g_t))
        ax.set_xlabel('t')
        ax.set_ylabel('correlation')

    fig.suptitle('System: ' + system_id, y=0.993)
    fig.tight_layout(rect=[0, 0, 1, 0.99])
    # plt.show()
    plt.savefig(os.path.join(PLOTS_DIR_PATH, system_id + '-traj.pdf'))


def plot_replicas_statistical_inefficiencies(system_ids):
    n_cols = len(system_ids)
    fig, axes = plt.subplots(nrows=2, ncols=len(system_ids), figsize=(4*n_cols, 4))

    for system_idx, system_id in enumerate(system_ids):
        statistical_inefficiencies = []
        all_n_bound_waters = []
        for n_bound_waters in iterate_bound_waters(system_id, state=False):
            g_t = pymbar.timeseries.statisticalInefficiency(n_bound_waters)
            statistical_inefficiencies.append(g_t)
            all_n_bound_waters.append(n_bound_waters)

        ax = axes[0][system_idx]
        sns.distplot(statistical_inefficiencies, kde=False, norm_hist=True, ax=ax)
        ax.set_title('System: {}'.format(system_id))
        ax.set_xlabel('statistical inefficiency [ps]')

        # Plot super-replica correlation function.
        ax = axes[1][system_idx]
        g_t, C_t = pymbar.timeseries.statisticalInefficiencyMultiple(all_n_bound_waters,
                                                                     return_correlation_function=True)
        t, C_t = list(zip(*C_t))
        ax.plot(t, C_t)
        ax.set_xlabel('t [ps]')
        ax.set_ylabel('Correlation')
        ax.set_title('Super-replica statistical inefficiency {:.2f} ps'.format(g_t))
    fig.suptitle('Water wetting/dewetting replica statistical inefficiency distribution')
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    # plt.show()
    plt.savefig(os.path.join(PLOTS_DIR_PATH, 'statistical-inefficiencies.pdf'))


def plot_distribution_n_bound_waters(system_id):
    n_states = get_n_states(system_id)
    n_cols = 4
    n_rows = int(np.ceil(n_states / n_cols))

    for trajectory_type in ['replica', 'state']:
        fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(12, 2*n_rows))

        state_flag = True if trajectory_type == 'state' else False
        n_bound_waters = [traj for traj in iterate_bound_waters(system_id, state=state_flag)]
        max_n_waters = max(max(traj) for traj in n_bound_waters)
        bins = np.linspace(-0.5, max_n_waters+0.5, max_n_waters+2)
        all_n_bound_waters = np.concatenate(n_bound_waters)

        for traj_idx, (ax, traj) in enumerate(zip(axes.flatten(), n_bound_waters)):
            sns.distplot(all_n_bound_waters, bins=bins, kde=False, norm_hist=True,
                         hist_kws={'edgecolor': 'w', 'lw': 2}, label='total distribution', ax=ax)
            sns.distplot(traj, bins=bins, kde=False, norm_hist=True,
                         hist_kws={'edgecolor': 'w', 'lw': 2}, label=trajectory_type+' distribution', ax=ax)
            g_t = pymbar.timeseries.statisticalInefficiency(traj)

            title = '{} {}'.format(trajectory_type, traj_idx)
            if not state_flag:
                title += '\nstatistical inefficiency {:.2f} ps'.format(g_t)
            ax.set_title(title)
            ax.set_xlabel('n bound waters')
            ax.legend()

        fig.suptitle('System: ' + system_id, y=0.993)
        fig.tight_layout(rect=[0, 0, 1, 0.99])
        # plt.show()
        plt.savefig(os.path.join(PLOTS_DIR_PATH, '{}-{}-dist.pdf'.format(system_id, trajectory_type)))


# =============================================================================
# PAPER FIGURES
# =============================================================================

def plot_replica_trajectories(replicas_state_indices, all_n_bound_waters, replica_idx,
                              xlabel, ax, step=1, state_index_ticks=True):
    """
    Plot the trajectory of replica state indices and the bound waters for the given
    replica on the same axis.
    """
    # Handle the step. The bound waters trajectories skip iteration 0.
    t = np.array(range(1, len(all_n_bound_waters[replica_idx]), step))  # In picoseconds.
    replica_bound_waters = all_n_bound_waters[replica_idx][t-1]
    replica_state_indices = replicas_state_indices[replica_idx][t]

    # Plot n bound water trajectories.
    ax.plot(t, replica_bound_waters, linestyle='None', marker='.', ms=0.8,
            color='black', alpha=0.7, label='n bound waters')
    ax.set_ylabel('n bound waters')
    ax.set_ylim((-0.2, np.max(all_n_bound_waters)+0.5))
    ax.set_xlim((0, len(all_n_bound_waters[replica_idx])))

    # Create state indices trajectories.
    ax2 = ax.twinx()
    scatter = ax2.scatter(t, replica_state_indices, marker='.', s=0.8,
                          vmin=0, vmax=np.max(replicas_state_indices),
                          label='state index', alpha=0.7, cmap='viridis',
                          c=replica_state_indices)
    if state_index_ticks:
        ax2.set_ylabel('state index')
    else:
        # ax.tick_params(axis='y', which='both', left=False, labelleft=False)
        ax2.set_yticks([])
    ax2.set_ylim((0, np.max(replicas_state_indices)))

    # Keep the number of bound waters trajectory in the front.
    ax.set_zorder(10)
    ax.patch.set_visible(False)

    # Add labels.
    ax.set_title('CB8-G3-0 replica ' + str(replica_idx), pad=1.8)
    if xlabel is not None:
        ax.set_xlabel(xlabel)

    return scatter


def plot_paper_trajectories(replica_indices):
    """
    Plot the part of Figure 5 containing the state trajectories and the
    bound water trajectories and correlation functions.
    """
    n_rows = len(replica_indices) + 1
    n_cols = 1
    fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(7.25/5*2.9, 3), sharex=True)

    correlation_palette = itertools.cycle(sns.color_palette("Paired"))

    # Compute statistical inefficiencies for all systems.
    for system_id in ['CB8-G3-0', 'OA-G3-0', 'OA-G6-0']:

        # Read in all state indices and water trajectories.
        file_path = os.path.join(ANALYSIS_DIR_PATH, system_id + '-replica-state-indices.npy')
        replica_state_indices = np.transpose(np.load(file_path))
        all_n_bound_waters = np.array(list(iterate_bound_waters(system_id, state=False)))

        # Plot the state indices and waters bound trajectories in the first rows.
        if system_id == 'CB8-G3-0':
            for row_idx, replica_idx in enumerate(replica_indices):
                # Collect the scatter plot to build a color bar for later.
                scatter = plot_replica_trajectories(replica_state_indices, all_n_bound_waters,
                                                    replica_idx, xlabel=None, ax=axes[row_idx],
                                                    step=1, state_index_ticks=False)
            fig.suptitle('Examples of replica trajectories', fontsize='medium')
            max_state_colorbar = np.max(replica_state_indices)

        # Plot the correlation functions.
        ax = axes[len(replica_indices)]
        for label, trajectories in [('state indices', replica_state_indices),
                                    ('n bound waters', all_n_bound_waters)]:
            # Compute the statistical inefficiency and the correlation function.
            g_t, C_t = pymbar.timeseries.statisticalInefficiencyMultiple(trajectories, fast=True,
                                                                         return_correlation_function=True)
            t = [x[0] for x in C_t]
            C_t = [x[1] for x in C_t]
            ax.plot(t, C_t, label='{} {} ({:.1f}ps)'.format(system_id, label, g_t),
                    color=next(correlation_palette), alpha=0.7)

        # Set labels.
        ax.set_ylim((0, 1))
        ax.set_title('Super replica correlation functions', pad=1.8)
        ax.set_xlabel('t [ps]')
        ax.set_ylabel('correlation')
        ax.tick_params(axis='x', which='major', pad=0.2)

    axes[-1].legend(fontsize='xx-small')
    fig.tight_layout(pad=0.1, rect=[0, 0, 1, 0.92])

    # Make space for the colorbar.
    fig.subplots_adjust(left=0.24)
    # Create new axis for the colorbar.
    colorbar_ax = fig.add_axes([0.0, 0.105, 0.05, 0.78])
    # Add the colorbar for the first scatter. Set the ticks
    # for the bound, discharged and decoupled states.
    colorbar = fig.colorbar(scatter, cax=colorbar_ax, ticks=[0, 25, max_state_colorbar])
    colorbar.ax.set_yticklabels(['bound', 'discharged', 'decoupled'], fontsize='x-small')
    colorbar.ax.tick_params(axis='y', which='major', pad=0.2)

    # plt.show()
    # plt.savefig(os.path.join(PAPER_IMAGES_DIR_PATH, 'Figure5-waters', 'trajectories.pdf'))
    plt.savefig(os.path.join(PAPER_IMAGES_DIR_PATH, 'Figure5-waters', 'trajectories.png'), dpi=600)


def plot_paper_distributions(system_id, skip):
    n_bound_waters = list(iterate_bound_waters(system_id, state=True))
    n_states = len(n_bound_waters)
    max_n_waters = max(max(traj) for traj in n_bound_waters)
    state_indices = list(range(0, len(n_bound_waters), skip))

    # Create bins sequence to feed to matplotlib.hist().
    bins = [-0.5+i for i in range(max_n_waters+2)]

    # Make sure the decoupled state is in.
    if state_indices[-1] != n_states - 1:
        state_indices.append(n_states - 1)

    # Select color.
    color_palette = itertools.cycle(sns.color_palette('viridis', n_colors=len(state_indices)))

    fig, ax = plt.subplots(figsize=(7.25/5*2, 3))
    for state_idx in state_indices:
        state_bound_waters = n_bound_waters[state_idx]
        if state_idx == 0:
            label = 'Bound'
        elif state_idx == n_states - 1:
            label = 'Decoupled'
        else:
            label = 'State {}'.format(state_idx+1)
        state_color = next(color_palette)
        sns.distplot(state_bound_waters, kde=True, hist=True,
                     kde_kws={'bw': 0.4, 'lw': 0.8, 'shade': False, 'color': state_color, 'label': label},
                     bins=bins, hist_kws={'color': state_color},
                     ax=ax)

    ax.set_xlim(-0.5, max_n_waters)
    ax.set_xticks(list(range(max_n_waters+1)))
    ax.tick_params(axis='x', which='major', pad=0.2)
    ax.tick_params(axis='y', which='major', pad=0.2)
    ax.set_title('Bound water molecules distribution\nby state in {}'.format(system_id))
    ax.set_ylabel('density')
    ax.set_xlabel('n bound waters')
    # ax.legend(fontsize='x-small', ncol=2, loc='upper right', bbox_to_anchor=(1.015, 1.015), fancybox=False)
    ax.get_legend().remove()

    fig.tight_layout(pad=0.1)
    # plt.show()
    fig.savefig(os.path.join(PAPER_IMAGES_DIR_PATH, 'Figure5-waters', '{}-distribution.pdf'.format(system_id)))
    # plt.savefig(os.path.join(PAPER_IMAGES_DIR_PATH, 'Figure5-waters', '{}-dist.png'.format(system_id), dpi=300))


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    # import argparse
    # parser = argparse.ArgumentParser(description='Compute the number of waters bound in time.')
    # parser.add_argument('-i', '--jobid', metavar='jobid', type=int, help='JOB ID (1-based).')
    # args = parser.parse_args()
    # if not args.jobid:
    #     job_id = None
    # else:
    #     job_id = args.jobid - 1
    # compute_all_n_bound_waters(job_id)

    # From here on, parallelization is not supported.
    # -----------------------------------------------

    # for system_id in ['CB8-G3-0', 'OA-G3-0', 'OA-G6-0']:
    #     extract_replica_state_indices(system_id)

    from matplotlib import pyplot as plt
    import seaborn as sns
    sns.set_context('paper', font_scale=0.7)
    sns.set_style('whitegrid')
    system_ids = ['CB8-G3-0', 'OA-G3-0', 'OA-G6-0']
    # system_ids = ['CB8-G3-0']

    # # Analysis plot.
    # # TODO: Make these available as SI figures?
    # plot_replicas_statistical_inefficiencies(system_ids)
    # # for system_id in system_ids:
    # #     plot_n_bound_waters_trajectory(system_id)
    # #     plot_distribution_n_bound_waters(system_id)

    # Paper figures.
    # plot_paper_trajectories(replica_indices=[1, 5])
    for system_id in system_ids:
        # plot_paper_trajectories_old(system_id, replica_indices=[1, 5, 10])
        plot_paper_distributions(system_id, skip=4)
