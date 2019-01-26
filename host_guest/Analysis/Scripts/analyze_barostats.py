import os

import netCDF4
import numpy as np
from simtk import unit


DATA_DIR_PATH = os.path.join('..', 'SAMPLing', 'Data', 'BarostatData')
YANK_VOLUMES_DIR_PATH = os.path.join(DATA_DIR_PATH, 'YankVolumes')
EE_VOLUMES_DIR_PATH = os.path.join(DATA_DIR_PATH, 'EEVolumes')


# =============================================================================
# EXTRACT/GENERATE THE DATA
# =============================================================================

def read_hrex_volumes(nc_file_path, replica_index=0, state_index=0):
    nc_file = netCDF4.Dataset(nc_file_path, 'r')
    try:
        # Read replica volumes.
        replica_volumes = nc_file.variables['volumes'][:, replica_index]

        # Read state volumes.
        replica_state_indices = nc_file.variables['states'][:]
        n_iterations = len(replica_state_indices)
        state_volumes = np.zeros(n_iterations)
        for iteration, iteration_replica_state_indices in enumerate(replica_state_indices):
            replica_index = np.where(iteration_replica_state_indices == state_index)[0][0]
            state_volumes[iteration] = nc_file.variables['volumes'][iteration, replica_index]
    finally:
        nc_file.close()
    return replica_volumes, state_volumes


def run_simulations(nc_file_path, pressure):
    n_iterations = 50000

    # Read the state and the MCMC move used for the simulation.
    from yank.repex import Reporter
    reporter = Reporter(nc_file_path, open_mode='r')
    try:
        print('Reading from NetCDF file...')
        thermodynamic_state = reporter.read_thermodynamic_states()[0][0]
        sampler_state = reporter.read_sampler_states(iteration=0, analysis_particles_only=False)[0]
        mcmc_move = reporter.read_mcmc_moves()[0]
    finally:
        reporter.close()

    # Run the simulations at the set pressure.
    print('Running at pressure {}...'.format(pressure))

    thermodynamic_state.pressure = pressure
    volumes = np.zeros(n_iterations)
    for iteration in range(n_iterations):
        mcmc_move.apply(thermodynamic_state, sampler_state)
        volumes[iteration] = sampler_state.volume / unit.nanometer**3

        # Feedback.
        if iteration + 1 % 1000 == 0:
            print('Completed {}/{} iterations'.format(iteration + 1, n_iterations))

    # Store volumes.
    os.makedirs(YANK_VOLUMES_DIR_PATH)
    pressure_str = '{:.2f}'.format(pressure / unit.atmosphere).replace('.', '')
    file_path = os.path.join(YANK_VOLUMES_DIR_PATH, 'volumes_pressure{}.npy'.format(pressure_str))
    np.save(file_path, volumes)


if __name__ == '__main__':
    nc_file_path = '../SAMPLing/experiments/experiment-OAG3/OAG30/complex.nc'

    import argparse
    parser = argparse.ArgumentParser(description='Run simulation at given pressure.')
    parser.add_argument('-p', '--pressure', metavar='pressure', type=float, help='Pressure in atm.')
    pressure = parser.parse_args().pressure
    # Run new simulations.
    run_simulations(nc_file_path, pressure*unit.atmosphere)

    # Read and store HREX volumes.
    hrex_replica_volumes, hrex_state_volumes = read_hrex_volumes(nc_file_path)
    np.save('hrex_replica_volumes.npy', hrex_replica_volumes)
    np.save('hrex_state_volumes.npy', hrex_state_volumes)

