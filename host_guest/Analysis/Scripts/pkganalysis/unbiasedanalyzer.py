#!/usr/local/bin/env python


# =============================================================================
# GLOBAL IMPORTS
# =============================================================================

import os
import copy
import scipy
import logging
import operator
import collections

import yaml
import numpy as np

import mdtraj
import pymbar
from simtk import openmm, unit as unit
import openmmtools as mmtools
import yank.restraints
from yank import analyze, repex, utils

logger = logging.getLogger(__name__)


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def get_restraint_force(system):
    """Extract the radially symmetric restraint Custom(Centroid)BondForce of the system."""
    restraint_force = None
    for i, force in enumerate(system.getForces()):
        if force.__class__.__name__ in ['CustomCentroidBondForce', 'CustomBondForce']:
            if force.getGlobalParameterName(0) == 'lambda_restraints':
                restraint_force = copy.deepcopy(force)
                break
    return restraint_force


def set_restrained_particles(restraint_force, particles1, particles2):
    try:
        # CustomCentroidBondForce
        restraint_force.setGroupParameters(0, list(particles1))
        restraint_force.setGroupParameters(1, list(particles2))
    except AttributeError:
        # CustomBondForce
        _, _, bond_parameters = restraint_force.getBondParameters(0)
        restraint_force.setBondParameters(0, particles1[0], particles2[0], bond_parameters)


def compute_centroid_distance(positions_group1, positions_group2, weights_group1, weights_group2):
    """Compute the distance between the centers of mass of the two groups.

    The two positions given must have the same units.

    Parameters
    ----------
    positions_group1 : numpy.array
        The positions of the particles in the first CustomCentroidBondForce group.
    positions_group2 : numpy.array
        The positions of the particles in the second CustomCentroidBondForce group.
    weights_group1 : list of float
        The mass of the particle in the first CustomCentroidBondForce group.
    weights_group2 : list of float
        The mass of the particles in the second CustomCentroidBondForce group.

    """
    assert len(positions_group1) == len(weights_group1)
    assert len(positions_group2) == len(weights_group2)
    # Compute center of mass for each group.
    com_group1 = np.average(positions_group1, axis=0, weights=weights_group1)
    com_group2 = np.average(positions_group2, axis=0, weights=weights_group2)
    # Compute distance between centers of mass.
    distance = np.linalg.norm(com_group1 - com_group2)
    return distance


# =============================================================================
# ANALYZER TO UNBIAS THE EFFECT OF THE RESTRAINT
# =============================================================================

# Dictionary cached_value: list of cache_values to invalidate.
_CACHE_VALUES_DEPENDENCIES = dict(
    n_discarded_initial_iterations=['equilibration_data'],
    fixed_statistical_inefficiency=['equilibration_data'],
    min_n_iterations=['equilibration_data'],
    max_n_iterations=['equilibration_data'],
    equilibration_data=['state_indices_kn', 'uncorrelated_u_kn', 'uncorrelated_N_k'],
    restraint_energy_cutoff=['unbiased_uncorrelated_u_kn', 'unbiased_uncorrelated_N_k'],
    restraint_distance_cutoff=['unbiased_uncorrelated_u_kn', 'unbiased_uncorrelated_N_k'],
    state_indices_kn=['unbiased_uncorrelated_u_kn', 'unbiased_uncorrelated_N_k'],
    uncorrelated_u_kn=['unbiased_uncorrelated_u_kn', 'unbiased_uncorrelated_N_k'],
    uncorrelated_N_k=['unbiased_uncorrelated_u_kn', 'unbiased_uncorrelated_N_k'],
    unbiased_uncorrelated_u_kn=['mbar'],
    unbiased_uncorrelated_N_k=['mbar'],
    mbar=['observables']
)


class UnbiasedAnalyzer(analyze.ReplicaExchangeAnalyzer):

    def __init__(self, *args, restraint_energy_cutoff=None, restraint_distance_cutoff=None,
                 min_n_iterations=0, max_n_iterations=None, n_discarded_initial_iterations=None,
                 fixed_statistical_inefficiency=None, **kwargs):
        # Cached values that are read directly from the Reporter.
        self._n_iterations = None
        self._n_replicas = None
        self._restraint_data = None
        self._kT = None
        self._energy_matrix = None
        self._unsampled_energy_matrix = None

        # This cache should be always set with _update_cache().
        self._cache = {}

        super().__init__(*args, **kwargs)

        if max_n_iterations is not None:
            self.max_n_iterations = max_n_iterations
        self.min_n_iterations = min_n_iterations
        self.n_discarded_initial_iterations = n_discarded_initial_iterations
        self.fixed_statistical_inefficiency = fixed_statistical_inefficiency
        self.restraint_energy_cutoff = restraint_energy_cutoff
        self.restraint_distance_cutoff = restraint_distance_cutoff

    @property
    def n_iterations(self):
        """int: The total number of iterations of the phase."""
        if self._n_iterations is None:
            # The + 1 accounts for iteration 0.
            self._n_iterations = self._reporter.read_last_iteration(full_iteration=False)
        return self._n_iterations

    @property
    def n_replicas(self):
        """int; Number of replicas."""
        if self._n_replicas is None:
            replica_state_indices = self._reporter.read_replica_thermodynamic_states(iteration=0)
            _, self._n_replicas = replica_state_indices.shape
        return self._n_replicas

    class _CachedProperty(object):
        """Descriptor of a cached value."""
        def __init__(self, name, validator=None, default_func=None, check_changes=False):
            assert name in _CACHE_VALUES_DEPENDENCIES
            assert name != 'observables'  # Reserved name.
            self._name = name
            self._validator = validator
            self._default_func = default_func
            self._check_changes = check_changes

        def __get__(self, instance, owner_class=None):
            try:
                return instance._cache[self._name]
            except KeyError:
                if self._default_func is not None:
                    return self._default_func(instance)
                raise

        def __set__(self, instance, new_value):
            if self._validator is not None:
                new_value = self._validator(instance, new_value)
            instance._update_cache(self._name, new_value, self._check_changes)

    @staticmethod
    def _max_n_iterations_validator(instance, new_value):
        if new_value > instance.n_iterations:
            new_value = instance.n_iterations
        return new_value

    n_discarded_initial_iterations = _CachedProperty('n_discarded_initial_iterations', check_changes=True)
    min_n_iterations = _CachedProperty('min_n_iterations', check_changes=True)
    max_n_iterations = _CachedProperty('max_n_iterations', validator=_max_n_iterations_validator.__func__,
                                       default_func=lambda instance: instance.n_iterations,
                                       check_changes=True)
    fixed_statistical_inefficiency = _CachedProperty('fixed_statistical_inefficiency', check_changes=True)
    restraint_energy_cutoff = _CachedProperty('restraint_energy_cutoff')
    restraint_distance_cutoff = _CachedProperty('restraint_distance_cutoff')

    def _automatic_equilibration_detection(self, input_data):
        """Perform automatic equilibration detection and set equilibration_data."""
        u_n = self.get_timeseries(input_data[:, :, self.min_n_iterations:self.max_n_iterations+1])

        # Discard equilibration samples.
        # TODO: if we include u_n[0] (the energy right after minimization) in the equilibration detection,
        # TODO:         then number_equilibrated is 0. Find a better way than just discarding first frame.
        if self.min_n_iterations == 0:
            u_n = u_n[1:]
        equilibration_data = list(analyze.get_equilibration_data(u_n))

        # Discard also minimization frame or minimum number of frames.
        if self.min_n_iterations == 0:
            equilibration_data[0] += 1
        else:
            equilibration_data[0] += self.min_n_iterations

        self._equilibration_data = tuple(equilibration_data)
        return self._equilibration_data

    def _get_equilibration_data_auto(self, input_data=None):
        """Modify original to discard minimization frame, log the number of discarded equilibration and keep max_n_iterations into account."""
        if input_data is None:
            input_data, _ = self.get_states_energies()

        # Check if we need to perform automatic equilibration detection or if the
        # user has already provided a number of initial iterations to discard.
        if self.n_discarded_initial_iterations is None and self.fixed_statistical_inefficiency is None:
            # Automatic equilibration detection.
            self._automatic_equilibration_detection(input_data)
        elif self.fixed_statistical_inefficiency is not None:
            assert self.min_n_iterations == 0
            statistical_inefficiency = self.fixed_statistical_inefficiency
            # Discard 2 times the statistical efficiency if possible, otherwise don't discard anything.
            n_equilibration_iterations = int(2 * statistical_inefficiency)
            if n_equilibration_iterations < self.max_n_iterations:
                logger.info('Could not discard 2*statistical_inefficiency initial iterations. '
                            'The entire data will be used to compute the free energy.')
                n_equilibration_iterations = 0
            n_effective_samples = (self.max_n_iterations - n_equilibration_iterations + 1) / statistical_inefficiency
            self._equilibration_data = (n_equilibration_iterations, statistical_inefficiency, n_effective_samples)
        else:
            u_n = self.get_timeseries(input_data[:, :, self.n_discarded_initial_iterations:self.max_n_iterations+1])
            # We set n_equilibration_iterations = n_discarded_initial_iterations
            # and we just compute the statistical inefficiency and the number of
            # effective samples starting from this iteration.
            n_equilibration_iterations = self.n_discarded_initial_iterations
            statistical_inefficiency = pymbar.timeseries.statisticalInefficiency(u_n)
            n_effective_samples = len(u_n) / statistical_inefficiency
            self._equilibration_data = (n_equilibration_iterations, statistical_inefficiency, n_effective_samples)

        logger.debug('Equilibration data: {}'.format(self._equilibration_data))
        return self._equilibration_data

    _equilibration_data = _CachedProperty('equilibration_data', check_changes=True,
                                          default_func=lambda instance: instance._get_equilibration_data_auto())

    @property
    def n_equilibration_iterations(self):
        return self._equilibration_data[0]

    @property
    def statistical_inefficiency(self):
        return self._equilibration_data[1]

    @property
    def _uncorrelated_iterations(self):
        """list of int: the indices of the uncorrelated iterations."""
        equilibrium_iterations = np.array(range(self.n_equilibration_iterations, self.n_iterations + 1))
        uncorrelated_iterations_indices = pymbar.timeseries.subsampleCorrelatedData(equilibrium_iterations,
                                                                                    self.statistical_inefficiency)
        return equilibrium_iterations[uncorrelated_iterations_indices]

    @staticmethod
    def _state_indices_kn_default_func(instance):
        uncorrelated_iterations = instance._uncorrelated_iterations  # Shortcut.
        replica_state_indices = instance._reporter.read_replica_thermodynamic_states()
        n_correlated_iterations, instance._n_replicas = replica_state_indices.shape

        # Initialize output array.
        n_frames = instance.n_replicas * len(uncorrelated_iterations)
        state_indices_kn = np.zeros(n_frames, dtype=np.int32)

        # Map kn columns to the sta
        for iteration_idx, iteration in enumerate(uncorrelated_iterations):
            for replica_idx in range(instance.n_replicas):
                # Deconvolute index.
                state_idx = replica_state_indices[iteration, replica_idx]
                frame_idx = state_idx*len(uncorrelated_iterations) + iteration_idx
                # Set output array.
                state_indices_kn[frame_idx] = state_idx
        instance._state_indices_kn = state_indices_kn
        return state_indices_kn

    _state_indices_kn = _CachedProperty('state_indices_kn', default_func=_state_indices_kn_default_func.__func__)
    _uncorrelated_u_kn = _CachedProperty('uncorrelated_u_kn', default_func=lambda x: None)
    _uncorrelated_N_k = _CachedProperty('uncorrelated_N_k', default_func=lambda x: None)
    _unbiased_uncorrelated_u_kn = _CachedProperty('unbiased_uncorrelated_u_kn', default_func=lambda x: None)
    _unbiased_uncorrelated_N_k = _CachedProperty('unbiased_uncorrelated_N_k', default_func=lambda x: None)
    _mbar = _CachedProperty('mbar', default_func=lambda x: None)

    def _update_cache(self, key, new_value, check_changes=False):
        """Update the cache and invalidate values that depend on it.

        Parameters
        ----------
        key : str
            The name of the value to update.
        value : object
            The new value of the key.
        check_changes : bool, optional
            If True and the new value is equal to the current one,
            the dependent cache values are not invalidated.

        """
        invalidate_cache = True
        try:
            old_value = self._cache[key]
        except KeyError:
            pass
        else:
            if check_changes and old_value == new_value:
                invalidate_cache = False
        # Update value and invalidate the cache.
        self._cache[key] = new_value
        if invalidate_cache:
            self._invalidate_cache_values(key)

    def _invalidate_cache_values(self, key):
        """Invalidate all the cache dependencies of key.

        Parameters
        ----------
        key : str
            The name of the cached whose dependencies must be invalidated.

        """
        for k in _CACHE_VALUES_DEPENDENCIES[key]:
            # Invalidate observables that are in a separate cache.
            if k == 'observables':
                for observable in self.observables:
                    self._computed_observables[observable] = None
            else:
                # Invalidate dependencies of k.
                self._invalidate_cache_values(k)
                # Remove k.
                self._cache.pop(k, None)

    def _read_thermodynamic_states(self):
        """Read thermodynamic states and caches useful info in the meantime."""
        thermodynamic_states, unsampled_states = self._reporter.read_thermodynamic_states()
        # TODO should we read all temperatures and let kT property depend on reference_states?
        self._kT = unsampled_states[0].kT  # Avoid reading TS again when we need kT.
        return thermodynamic_states, unsampled_states

    def get_states_energies(self):
        """Modify original method to cache energies."""
        if self._energy_matrix is None:
            self._energy_matrix, self._unsampled_energy_matrix = super().get_states_energies()
        return self._energy_matrix, self._unsampled_energy_matrix

    def _compute_free_energy(self):
        """Modify original method to cache the last free energy to speed up next MBAR."""
        super()._compute_free_energy()
        self._extra_analysis_kwargs['initial_f_k'] = self.mbar.f_k

    def _create_mbar_from_scratch(self):
        if self._uncorrelated_u_kn is None:
            # Compute everything.
            super()._create_mbar_from_scratch()
        else:
            # Check if we need to compute the unbiased energies.
            if self._unbiased_uncorrelated_u_kn is None:
                self._compute_unbiased_mbar_data()
            self._create_mbar(self._unbiased_uncorrelated_u_kn, self._unbiased_uncorrelated_N_k)

    def _prepare_mbar_input_data(self, sampled_energy_matrix, unsampled_energy_matrix):
        """Convert the sampled and unsampled energy matrices into MBAR ready data"""
        u_kln, N_k = super()._prepare_mbar_input_data(sampled_energy_matrix, unsampled_energy_matrix)
        u_kn = pymbar.utils.kln_to_kn(u_kln, N_k)
        self._uncorrelated_u_kn = u_kn
        self._uncorrelated_N_k = N_k
        self._compute_unbiased_mbar_data()
        return self._unbiased_uncorrelated_u_kn, self._unbiased_uncorrelated_N_k

    def _compute_unbiased_mbar_data(self):
        """Unbias the restraint, apply energy/distance cutoff and cut to max_n_iterations."""
        # Don't modify the cached uncorrelated energies.
        u_kn = copy.deepcopy(self._uncorrelated_u_kn)
        N_k = copy.deepcopy(self._uncorrelated_N_k)
        n_uncorrelated_iterations_kn = u_kn.shape[1]

        # Check if we need to unbias the restraint.
        try:
            restraint_data = self._get_restraint_data()
        except TypeError as e:
            logger.info(str(e) + ' The restraint will not be unbiased.')
            unbias_restraint = False
        else:
            unbias_restraint = True

        is_cutoff_distance = unbias_restraint and self.restraint_distance_cutoff is not None
        is_cutoff_energy = unbias_restraint and self.restraint_energy_cutoff is not None
        if unbias_restraint and self.restraint_distance_cutoff is None and self.restraint_energy_cutoff is None:
            logger.info('No cutoff specified. The restraint will not be unbiased.')
            unbias_restraint = False

        # If we don't need to unbias the restraint and max_n_iterations is not set, there's nothing to do.
        if not unbias_restraint and self.max_n_iterations == self.n_iterations:
            self._unbiased_uncorrelated_u_kn = u_kn
            self._unbiased_uncorrelated_N_k = N_k
            return

        # If we need to unbias the restraint, compute the restraint
        # energies/distances and fix the standard state correction.
        if unbias_restraint:
            reduced_system, restraint_force = restraint_data[:2]
            particle_indices_group1, particle_indices_group2 = restraint_data[2:4]
            weights_group1, weights_group2 = restraint_data[4:]

            # Compute restraint energies/distances.
            distances_kn, energies_kn = self._compute_restrain_energies(
                particle_indices_group1, particle_indices_group2, weights_group1, weights_group2,
                reduced_system, compute_distances=is_cutoff_distance)

            # Convert energies to kT unit for comparison to energy cutoff.
            energies_kn = energies_kn / self.kT
            logger.debug('Restraint energy mean: {} kT; std: {} kT'
                         ''.format(np.mean(energies_kn), np.std(energies_kn, ddof=1)))
            assert len(energies_kn) == n_uncorrelated_iterations_kn

        # We need to take into account the initial unsampled states to index correctly N_k.
        state_idx_shift = 0
        while N_k[state_idx_shift] == 0:
            state_idx_shift +=1

        # Determine which samples are outside the cutoffs or have to be truncated.
        uncorrelated_iterations = self._uncorrelated_iterations
        columns_to_keep = []
        for iteration_kn in range(n_uncorrelated_iterations_kn):
            # Convert kn iteration to the actual iteration number.
            uncorrelated_iteration_idx = iteration_kn % len(uncorrelated_iterations)
            iteration = uncorrelated_iterations[uncorrelated_iteration_idx]
            if (iteration > self.max_n_iterations or
                    (is_cutoff_energy and energies_kn[iteration_kn] > self.restraint_energy_cutoff) or
                    (is_cutoff_distance and distances_kn[iteration_kn] > self.restraint_distance_cutoff)):
                # Update the number of samples generated from its state.
                state_idx = self._state_indices_kn[iteration_kn]
                N_k[state_idx + state_idx_shift] -= 1
            else:
                columns_to_keep.append(iteration_kn)

        # Drop all columns that exceed the cutoff(s).
        n_discarded = n_uncorrelated_iterations_kn - len(columns_to_keep)
        logger.debug('Discarding {}/{} samples outside the cutoffs (max_n_iterations: {}, '
                     'restraint_distance_cutoff: {}, restraint_energy_cutoff: {}).'.format(
            n_discarded, n_uncorrelated_iterations_kn, self.max_n_iterations,
            self.restraint_distance_cutoff, self.restraint_energy_cutoff))
        u_kn = u_kn[:, columns_to_keep]

        # Add new end states that don't include the restraint.
        if unbias_restraint:
            energies_kn = energies_kn[columns_to_keep]
            n_states, n_iterations = u_kn.shape
            n_states_new = n_states + 2
            N_k_new = np.zeros(n_states_new, N_k.dtype)
            u_kn_new = np.zeros((n_states_new, n_iterations), u_kn.dtype)
            u_kn_new[0, :] = u_kn[0] - energies_kn
            u_kn_new[-1, :] = u_kn[-1] - energies_kn
            # Copy old values.
            N_k_new[1:-1] = N_k
            u_kn_new[1:-1, :] = u_kn
        else:
            u_kn_new = u_kn
            N_k_new = N_k

        # Cache new values.
        self._unbiased_uncorrelated_u_kn = u_kn_new
        self._unbiased_uncorrelated_N_k = N_k_new

    def _get_restraint_data(self):
        """Return the two unsampled states and a reduced version of them containing only the restraint force."""
        # Check cached value.
        if self._restraint_data is not None:
            return copy.deepcopy(self._restraint_data)

        # Isolate the end states.
        sampled_states, unsampled_states = self._read_thermodynamic_states()
        if len(unsampled_states) == 0:
            end_states = [sampled_states[0], sampled_states[-1]]
        else:
            end_states = unsampled_states

        # Isolate restraint force.
        system = end_states[0].system
        restraint_force = get_restraint_force(system)

        # Check this is a radially symmetric restraint and it was turned on at the end states.
        if restraint_force is None:
            raise TypeError('Cannot find a radially symmetric restraint.')
        if end_states[0].lambda_restraints != 1.0 or end_states[-1].lambda_restraints != 1.0:
            raise TypeError('Cannot unbias a restraint that is turned off at one of the end states.')

        # Log bond parameters.
        bond_parameters = restraint_force.getBondParameters(0)[-1]
        try:  # FlatBottom
            logger.debug('Bond parameters: K={}, r0={}'.format(*bond_parameters))
        except IndexError:  # Harmonic
            logger.debug('Bond parameters: K={}'.format(*bond_parameters))

        # Obtain restraint's particle indices to compute restraint distance.
        try:
            # CustomCentroidBondForce
            particle_indices_group1, weights_group1 = restraint_force.getGroupParameters(0)
            particle_indices_group2, weights_group2 = restraint_force.getGroupParameters(1)
            assert len(weights_group1) == 0  # Use masses to compute centroid.
            assert len(weights_group2) == 0  # Use masses to compute centroid.
        except AttributeError:
            # CustomBondForce
            particle_indices_group1, particle_indices_group2, _ = restraint_force.getBondParameters(0)
            particle_indices_group1 = [particle_indices_group1]  # Convert to list.
            particle_indices_group2 = [particle_indices_group2]  # Convert to list.

        # Convert tuples of np.integers to lists of ints.
        particle_indices_group1 = [int(i) for i in sorted(particle_indices_group1)]
        particle_indices_group2 = [int(i) for i in sorted(particle_indices_group2)]
        logger.debug('receptor restrained atoms: {}'.format(particle_indices_group1))
        logger.debug('ligand restrained atoms: {}'.format(particle_indices_group2))

        # Create new system with only solute and restraint forces.
        reduced_system = openmm.System()
        for particle_indices_group in [particle_indices_group1, particle_indices_group2]:
            for i in particle_indices_group:
                reduced_system.addParticle(system.getParticleMass(i))
        reduced_system.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())

        # Compute weights restrained particles.
        weights_group1 = [system.getParticleMass(i) for i in particle_indices_group1]
        weights_group2 = [system.getParticleMass(i) for i in particle_indices_group2]

        # Adapt the restraint force atom indices to the reduced system.
        assert max(particle_indices_group1) < min(particle_indices_group2)
        n_atoms_group1 = len(particle_indices_group1)
        tot_n_atoms = n_atoms_group1 + len(particle_indices_group2)
        restraint_force = copy.deepcopy(restraint_force)
        set_restrained_particles(restraint_force, particles1=range(n_atoms_group1),
                                 particles2=range(n_atoms_group1, tot_n_atoms))
        reduced_system.addForce(restraint_force)

        self._restraint_data = (reduced_system, restraint_force,
                                particle_indices_group1, particle_indices_group2,
                                weights_group1, weights_group2)
        return copy.deepcopy(self._restraint_data)

    def _compute_restrain_energies(self, particle_indices_group1, particle_indices_group2,
                                   weights_group1, weights_group2, reduced_system,
                                   compute_distances=False):
        """Compute the restrain distances for the given iterations.

        Parameters
        ----------
        particle_indices_group1 : list of int
            The particle indices of the first CustomCentroidBondForce group.
        particle_indices_group2 : list of int
            The particle indices of the second CustomCentroidBondForce group.
        weights_group1 : list of float
            The mass of the particle in the first CustomCentroidBondForce group.
        weights_group2 : list of float
            The mass of the particles in the second CustomCentroidBondForce group.

        Returns
        -------
        restrain_distances_kn : np.array
            The restrain distances.
        """
        ENERGY_UNIT = unit.kilojoules_per_mole
        MDTRAJ_DISTANCE_UNIT = unit.nanometers
        uncorrelated_iterations = self._uncorrelated_iterations  # Shortcut.
        uncorrelated_iterations_set = set(uncorrelated_iterations)

        energies_kn = None
        distances_kn = None

        def uncorrelate(array_kn):
            iterations_kn_to_keep = []
            for iteration_kn in range(len(array_kn)):
                # Convert kn iteration to the actual iteration number.
                iteration = iteration_kn % (self.n_iterations + 1)  # Include minimization iteration.
                # Discard equilibration and correlated iterations.
                if iteration in uncorrelated_iterations_set:
                    iterations_kn_to_keep.append(iteration_kn)
            # Return uncorrelated version of the array.
            return copy.deepcopy(array_kn[iterations_kn_to_keep])

        # Check cached values.
        # dir_path = os.path.dirname(self._reporter.filepath)
        dir_path = os.path.basename(os.path.dirname(self._reporter.filepath))  # TODO REMOVE ME
        dir_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'cache', dir_path)  # TODO REMOVE ME
        logger.debug('Searching for restraint cache in {}'.format(dir_path))
        os.makedirs(dir_path, exist_ok=True)  # TODO MOVE ME BELOW

        restraint_energies_file_path = os.path.join(dir_path, 'restraint_energies_cache.npz')
        restraint_distances_file_path = os.path.join(dir_path, 'restraint_distances_cache.npz')
        if compute_distances and os.path.exists(restraint_distances_file_path):
            distances_kn = uncorrelate(np.load(restraint_distances_file_path)['arr_0'])
            compute_distances = False
        if os.path.exists(restraint_energies_file_path):
            energies_kn = uncorrelate(np.load(restraint_energies_file_path)['arr_0'])

        if energies_kn is not None and (not compute_distances or distances_kn is not None):
            # Set MDTraj units if distances array has been loaded.
            if distances_kn is not None:
                distances_kn = distances_kn * MDTRAJ_DISTANCE_UNIT
            return distances_kn, energies_kn * ENERGY_UNIT

        # subset_particles_indices = list(self._reporter.analysis_particle_indices)
        subset_particles_indices = particle_indices_group1 + particle_indices_group2
        replica_state_indices = self._reporter.read_replica_thermodynamic_states()
        n_correlated_iterations, self._n_replicas = replica_state_indices.shape

        # Create output arrays. We unfold the replicas the same way
        # it is done during the kln_to_kn conversion.
        n_frames_kn = self.n_replicas * (self.n_iterations + 1)  # The +1 is for the minimization data.
        energies_kn = np.zeros(n_frames_kn, dtype=np.float64) * ENERGY_UNIT

        if compute_distances:
            n_atoms = len(subset_particles_indices)
            n_atoms_group1 = len(particle_indices_group1)
            traj_particle_indices_group1 = list(range(n_atoms_group1))
            traj_particle_indices_group2 = list(range(n_atoms_group1, n_atoms))

            # Create topology of the restrained atoms.
            serialized_topography = self._reporter.read_dict('metadata/topography')
            topology = mmtools.utils.deserialize(serialized_topography).topology
            topology = topology.subset(subset_particles_indices)

            distances_kn = np.zeros(n_frames_kn, dtype=np.float32)
            # Initialize trajectory object needed for imaging molecules.
            trajectory = mdtraj.Trajectory(xyz=np.zeros((n_atoms, 3)), topology=topology)

        # Create context used to compute the energies.
        integrator = openmm.VerletIntegrator(1.0*unit.femtosecond)
        context = openmm.Context(reduced_system, integrator)

        # TODO: Better handling of file cache to compute
        # TODO:     energies/distances only at uncorrelated iterations.
        # iterations_to_compute = uncorrelated_iterations
        iterations_to_compute = list(range(self.n_iterations + 1))

        # Pre-computing distances.
        logger.debug('Computing restraint energies...')
        for iteration_idx, iteration in enumerate(iterations_to_compute):
            if iteration % 5000 == 0:
                logger.debug('Computed restraint energies of {}/{}...'.format(iteration, self.n_iterations))

            # Obtain solute only sampler states.
            sampler_states = self._reporter.read_sampler_states(iteration=iteration,
                                                                analysis_particles_only=True)

            for replica_idx, sampler_state in enumerate(sampler_states):
                # Deconvolute index.
                state_idx = replica_state_indices[iteration, replica_idx]
                frame_idx = state_idx*len(iterations_to_compute) + iteration_idx

                sliced_sampler_state = sampler_state[subset_particles_indices]
                sliced_sampler_state.apply_to_context(context)
                potential_energy = context.getState(getEnergy=True).getPotentialEnergy()
                energies_kn[frame_idx] = potential_energy

                if compute_distances:
                    # Update trajectory positions/box vectors.
                    trajectory.xyz = (sampler_state.positions[subset_particles_indices] / unit.nanometers).astype(np.float32)
                    trajectory.unitcell_vectors = np.array([sampler_state.box_vectors / unit.nanometers], dtype=np.float32)
                    trajectory.image_molecules(inplace=True, make_whole=False)
                    positions_group1 = trajectory.xyz[0][traj_particle_indices_group1]
                    positions_group2 = trajectory.xyz[0][traj_particle_indices_group2]

                    # Set output arrays.
                    distances_kn[frame_idx] = compute_centroid_distance(positions_group1, positions_group2,
                                                                        weights_group1, weights_group2)

        # Save cache files and set MDTraj units for distances.
        np.savez_compressed(restraint_energies_file_path, energies_kn / ENERGY_UNIT)
        if compute_distances:
            np.savez_compressed(restraint_distances_file_path, distances_kn)
            distances_kn = distances_kn * MDTRAJ_DISTANCE_UNIT
        return uncorrelate(distances_kn), uncorrelate(energies_kn)

    def _compute_standard_state_correction(self, restraint_force, unbiased=False, max_distance=None):
        """Compute the standard state correction."""
        # TODO refactor: the redundant code with yank.restraints.RadiallySymmetricRestraint._standard_state_correction
        r_min = 0 * unit.nanometers
        if self.restraint_distance_cutoff is not None:
            r_max = self.restraint_distance_cutoff
        elif max_distance is not None:
            r_max = max_distance
        else:
            r_max = 100 * unit.nanometers

        # Create a System object containing two particles connected by the reference force
        system = openmm.System()
        system.addParticle(1.0 * unit.amu)
        system.addParticle(1.0 * unit.amu)
        force = copy.deepcopy(restraint_force)
        set_restrained_particles(force, particles1=[0], particles2=[1])
        # Disable the PBC if on for this approximation of the analytical solution
        force.setUsesPeriodicBoundaryConditions(False)
        system.addForce(force)

        # Create a Reference context to evaluate energies on the CPU.
        integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
        platform = openmm.Platform.getPlatformByName('Reference')
        context = openmm.Context(system, integrator, platform)

        # Set default positions.
        positions = unit.Quantity(np.zeros([2,3]), unit.nanometers)
        context.setPositions(positions)

        # Create a function to compute integrand as a function of interparticle separation.
        beta = 1 / self.kT

        def integrand(r):
            """
            Parameters
            ----------
            r : float
                Inter-particle separation in nanometers

            Returns
            -------
            dI : float
               Contribution to integrand (in nm^2).

            """
            positions[1, 0] = r * unit.nanometers
            context.setPositions(positions)
            state = context.getState(getEnergy=True)
            potential = beta * state.getPotentialEnergy()  # In kT.

            if (self.restraint_energy_cutoff is not None and
                        potential > self.restraint_energy_cutoff):
                return 0.0
            elif unbiased:
                potential = 0.0

            dI = 4.0 * np.pi * r**2 * np.exp(-potential)
            return dI

        # Integrate shell volume.
        shell_volume, shell_volume_error = scipy.integrate.quad(lambda r: integrand(r), r_min / unit.nanometers,
                                                                r_max / unit.nanometers) * unit.nanometers**3

        # Compute standard-state volume for a single molecule in a box of
        # size (1 L) / (avogadros number). Should also generate constant V0.
        liter = 1000.0 * unit.centimeters**3  # one liter
        standard_state_volume = liter / (unit.AVOGADRO_CONSTANT_NA*unit.mole)  # standard state volume

        # Compute standard state correction for releasing shell restraints into standard-state box (in units of kT).
        DeltaG = - np.log(standard_state_volume / shell_volume)

        # Return standard state correction (in kT).
        return DeltaG

    def get_standard_state_correction(self):
        """
        Compute the standard state correction free energy associated with the Phase.

        This usually is just a stored variable, but it may need other calculations.

        Returns
        -------
        standard_state_correction : float
            Free energy contribution from the standard_state_correction

        """
        if self._computed_observables['standard_state_correction'] is not None:
            return self._computed_observables['standard_state_correction']

        # Determine if we need to recompute the standard state correction.
        compute_ssc = True
        try:
            reduced_system, restraint_force = self._get_restraint_data()[:2]
        except TypeError:
            compute_ssc = False

        if compute_ssc:
            # TODO: This code is redundant with yank.py.
            # TODO: Compute average box volume here?
            box_vectors = reduced_system.getDefaultPeriodicBoxVectors()
            box_volume = mmtools.states._box_vectors_volume(box_vectors)
            ssc = - np.log(yank.restraints.V0 / box_volume)
            if self.restraint_distance_cutoff is not None or self.restraint_energy_cutoff is not None:
                max_dimension = np.max(unit.Quantity(box_vectors) / unit.nanometers) * unit.nanometers
                ssc_cutoff = self._compute_standard_state_correction(restraint_force, unbiased=True,
                                                                     max_distance=max_dimension)
                # The restraint volume can't be bigger than the box volume.
                if ssc_cutoff < ssc:
                    ssc = ssc_cutoff

            # Update observable.
            self._computed_observables['standard_state_correction'] = ssc
            logger.debug('Computed a new standard state correction of {} kT'.format(ssc))

        # This reads the SSC from the reporter if compute_ssc is False.
        return super().get_standard_state_correction()


# =============================================================================
# ANALYZER FOR AUTOCORRELATION ANALYSIS WITH INSTANTANEOUS WORK TIMESERIES
# =============================================================================

class InstantaneousWorkAnalyzer(UnbiasedAnalyzer):

    def _automatic_equilibration_detection(self, input_data):
        """Perform automatic equilibration detection using the instantaneous work values."""
        # Compute timeseries of average instantaneous work w.r.t. neighbor states.
        n_iterations = self.max_n_iterations + 1 - self.min_n_iterations
        work_timeseries = np.zeros([n_iterations], dtype=input_data.dtype)

        for iteration in range(self.min_n_iterations, self.max_n_iterations+1):
            iteration_energies = input_data[:, :, iteration]
            sampled_energies = np.diagonal(iteration_energies)
            # Use forward/backward difference for first and last state.
            forward_work = np.diagonal(iteration_energies, offset=1) - sampled_energies[:-1]
            backward_work = np.diagonal(iteration_energies, offset=-1) - sampled_energies[1:]
            avg_work = (forward_work[1:] - backward_work[:-1]) / 2
            work_timeseries[iteration] = forward_work[0] + np.sum(avg_work) - backward_work[-1]

        # Always throw away first frame (right after minimization).
        if self.min_n_iterations == 0:
            work_timeseries = work_timeseries[1:]
        equilibration_data = list(pymbar.timeseries.detectEquilibration(work_timeseries))
        if self.min_n_iterations == 0:
            equilibration_data[0] += 1
        else:
            equilibration_data[0] += self.min_n_iterations

        self._equilibration_data = tuple(equilibration_data)
        return self._equilibration_data


# =============================================================================
# ANALYZER TO ANALYZE WITH BAR
# =============================================================================

class BARAnalyzer(UnbiasedAnalyzer):
    """An UnbiasedAnalyzer that uses BAR instead of MBAR."""

    def _create_mbar(self, energy_matrix, samples_per_state):
        logger.debug("Initializing BAR...")
        self.mbar = MultiBAR(energy_matrix, samples_per_state)
        logger.debug('Done.')
        return self.mbar


class MultiBAR:
    """This class implement BAR analysis between neighbor states with the MBAR object API."""

    def __init__(self, u_kn, N_k, **kwargs):
        from pymbar import BAR

        # This attribute is used in the MultiStateSamplerAnalyzer.
        self.N_k = N_k

        # Determine number of states and iterations.
        # TODO: this is wrong. The only way to associate each sample
        #   to each state is to check N_k. For example, the first N_k[0]
        #   iterations are from state 0.
        self.n_states = u_kn.shape[0]
        n_iterations = u_kn.shape[1] / self.n_states
        assert float(n_iterations).is_integer(), (u_kn.shape, self.n_states, n_iterations, N_k.shape)
        self.n_states = int(self.n_states)
        n_iterations = int(n_iterations)

        # Compute the pairwise free energies with BAR.
        Deltaf_i = np.zeros(self.n_states-1)
        dDeltaf_i = np.zeros(self.n_states-1)
        for iteration in range(n_iterations):
            for state_i in range(self.n_states-1):
                # u_kn[k, s*n+i] is the reduced potential of a i-th configuration sampled from
                # state s and evaluated at state k, where n is the total number of iterations.
                state_j = state_i + 1

                # Here, by u_ij are all the energies sampled from state i and evaluated at state j.
                u_ii = np.array([u_kn[state_i, state_i*n_iterations + n] for n in range(n_iterations)])
                u_ij = np.array([u_kn[state_j, state_i*n_iterations + n] for n in range(n_iterations)])
                u_jj = np.array([u_kn[state_j, state_j*n_iterations + n] for n in range(n_iterations)])
                u_ji = np.array([u_kn[state_i, state_j*n_iterations + n] for n in range(n_iterations)])

                # Compute forward and reverse work and compute free energy BAR.
                w_F = u_ij - u_ii
                w_R = u_ji - u_jj
                Deltaf_i[state_i], dDeltaf_i[state_i] = BAR(w_F, w_R)

        # Transform Deltaf_ij in matrix form to make it consistent with MBAR.
        self.Deltaf_ij = np.zeros((self.n_states, self.n_states))
        self.dDeltaf_ij = np.zeros((self.n_states, self.n_states))
        for state_i in range(self.n_states):
            # Deltaf and uncertainty i == j is 0.
            for state_j in range(state_i+1, self.n_states):
                self.Deltaf_ij[state_i,state_j] = np.sum(Deltaf_i[state_i:state_j])
                self.dDeltaf_ij[state_i,state_j] = np.sqrt(np.sum(dDeltaf_i[state_i:state_j]**2))

        # Complete the rest of the antisymmetric matrix.
        self.Deltaf_ij -= self.Deltaf_ij.transpose()
        self.dDeltaf_ij += self.dDeltaf_ij.transpose()


    def getFreeEnergyDifferences(self):
        return self.Deltaf_ij, self.dDeltaf_ij


# =============================================================================
# SHORTCUT ANALYSIS FUNCTIONS
# =============================================================================

def analyze_phase(analyzer, return_enthalpy=False):
    data = dict()
    Deltaf_ij, dDeltaf_ij = analyzer.get_free_energy()
    data['DeltaF'] = Deltaf_ij[analyzer.reference_states[0], analyzer.reference_states[1]]
    data['dDeltaF'] = dDeltaf_ij[analyzer.reference_states[0], analyzer.reference_states[1]]
    data['DeltaF_standard_state_correction'] = analyzer.get_standard_state_correction()

    if return_enthalpy:
        try:
            DeltaH_ij, dDeltaH_ij = analyzer.get_enthalpy()
        except:
            data['DeltaH'] = None
            data['dDeltaH'] = None
        else:
            data['DeltaH'] = DeltaH_ij[analyzer.reference_states[0], analyzer.reference_states[1]]
            data['dDeltaH'] = dDeltaH_ij[analyzer.reference_states[0], analyzer.reference_states[1]]
    return data


def analyze_directory(source_directory, energy_cutoffs=None, distance_cutoffs=None,
                      iteration_cutoffs=None, n_discarded_initial_iterations=None,
                      fixed_statistical_inefficiency=None,
                      solvent_df=None, solvent_ddf=None, analyzer_class=None,
                      return_decomposition=False):
    """Run analysis on a set of cutoffs.

    Parameters
    ----------
    return_decomposition : bool, optional, default=False
        If True, not only all the free energies, but also their decomposition
        by phase and enthalpy are returned.

    Returns
    -------
    all_free_energies : List[Tuple[Quantity, Quantity]]
        The list of (free energies, uncertainty) pairs for all cutoffs.
    all_decompositions : List[Dict[str, Dict]]
        The list of free energy decompositions for all cutoffs indexed
        by a string including the phase name and the cutoff.

    See Also
    --------
    yank.analyze.analyze_phase
        For the data structure returned by the decomposition analysis.

    """
    if analyzer_class is None:
        analyzer_class = UnbiasedAnalyzer

    # Validate the input arguments. A Quantity is an iterable even
    # if it has only a single element so we need to check its value.
    is_distance_cutoff_iterable = (distance_cutoffs is not None and
                                   isinstance(distance_cutoffs._value, collections.Iterable))
    is_interval_cutoff_iterable = isinstance(iteration_cutoffs, collections.Iterable)
    if not operator.xor(
            isinstance(energy_cutoffs, collections.Iterable),
            operator.xor(is_interval_cutoff_iterable, is_distance_cutoff_iterable)
    ):
        raise ValueError('Only one between energy_cutoffs, distance_cutoffs, '
                         'and profile_convergence can be set.')

    # Check the variable that we want to profile.
    if isinstance(energy_cutoffs, collections.Iterable):
        cutoffs = energy_cutoffs
        cutoff_attribute = 'restraint_energy_cutoff'
    elif is_distance_cutoff_iterable:
        cutoffs = distance_cutoffs
        cutoff_attribute = 'restraint_distance_cutoff'
    else:
        cutoffs = iteration_cutoffs
        cutoff_attribute = 'max_n_iterations'

    analysis_script_path = os.path.join(source_directory, 'analysis.yaml')
    with open(analysis_script_path, 'r') as f:
        analysis = yaml.load(f)

    def cutoff_phase_name(phase_name, cutoff):
        if (is_interval_cutoff_iterable or
                        phase_name == 'complex' and 'restraint' in cutoff_attribute):
            return phase_name + '-' + str(cutoff)
        return phase_name

    data = collections.OrderedDict()
    calculation_type = ''
    for phase_name, sign in analysis:
        # Attempt to guess type of calculation
        if 'complex' in phase_name:
            calculation_type = ' of binding'
        elif 'solvent1' in phase_name:
            calculation_type = ' of solvation'

        # Avoid recomputing the solvent phase if known.
        if (phase_name == 'solvent' and not is_interval_cutoff_iterable and
                    solvent_df is not None and solvent_ddf is not None):
            data[phase_name] = {}
            data[phase_name]['DeltaF'] = solvent_df
            data[phase_name]['dDeltaF'] = solvent_ddf
            data[phase_name]['DeltaF_standard_state_correction'] = 0.0
            continue

        # Create kwargs to feed to unbiased analyzer.
        analyzer_kwargs = {'n_discarded_initial_iterations': n_discarded_initial_iterations,
                           'fixed_statistical_inefficiency': fixed_statistical_inefficiency}
        if phase_name == 'complex':
            if 'energy' not in cutoff_attribute:
                analyzer_kwargs['restraint_energy_cutoff'] = energy_cutoffs
            if 'distance' not in cutoff_attribute:
                analyzer_kwargs['restraint_distance_cutoff'] = distance_cutoffs

        phase_path = os.path.join(source_directory, phase_name + '.nc')
        reporter = repex.Reporter(phase_path, open_mode='r')
        phase = analyzer_class(reporter, **analyzer_kwargs)
        kT = phase.kT

        # For the restraint cutoff, we want to analyze exclusively the complex phase.
        if (is_interval_cutoff_iterable or
                (phase_name == 'complex' and 'restraint' in cutoff_attribute)):
            for cutoff_idx, cutoff in enumerate(cutoffs):
                logger.info('Analyzing cutoff {} ({}/{}) of {}'.format(cutoff, cutoff_idx+1,
                                                                       len(cutoffs), phase_name))
                setattr(phase, cutoff_attribute, cutoff)
                data[cutoff_phase_name(phase_name, cutoff)] = analyze_phase(phase, return_enthalpy=return_decomposition)
        else:
            data[phase_name] = analyze_phase(phase, return_enthalpy=return_decomposition)

    kT_to_kcalmol = kT / unit.kilocalories_per_mole

    # Compute free energy and enthalpy for all cutoffs.
    all_free_energies = []
    all_sscs = []
    for cutoff in cutoffs:
        DeltaF = 0.0
        dDeltaF = 0.0
        phase_names = []
        for phase_name, sign in analysis:
            phase_name = cutoff_phase_name(phase_name, cutoff)
            phase_names.append(phase_name)
            DeltaF -= sign * (data[phase_name]['DeltaF'] + data[phase_name]['DeltaF_standard_state_correction'])
            dDeltaF += data[phase_name]['dDeltaF']**2
        dDeltaF = np.sqrt(dDeltaF)

        all_free_energies.append((DeltaF * kT_to_kcalmol, dDeltaF * kT_to_kcalmol))
        all_sscs.append(data[cutoff_phase_name('complex', cutoff)]['DeltaF_standard_state_correction'] * kT_to_kcalmol)

        # Print energies
        logger.info('')
        logger.info('Reporting free energy for cutoff: {}'.format(cutoff))
        logger.info('-------------------------------------')
        logger.info('Free energy{:<13}: {:9.3f} +- {:.3f} kT ({:.3f} +- {:.3f} kcal/mol)'.format(
            calculation_type, DeltaF, dDeltaF, DeltaF * kT / unit.kilocalories_per_mole,
            dDeltaF * kT / unit.kilocalories_per_mole))

        for phase in phase_names:
            logger.info('DeltaG {:<17}: {:9.3f} +- {:.3f} kT'.format(phase, data[phase]['DeltaF'],
                                                                     data[phase]['dDeltaF']))
            if data[phase]['DeltaF_standard_state_correction'] != 0.0:
                logger.info('DeltaG {:<17}: {:18.3f} kT'.format('restraint',
                                                                data[phase]['DeltaF_standard_state_correction']))

    logger.info('all_free_energies={}'.format(all_free_energies))
    logger.info('sscs={}'.format(all_sscs))

    if return_decomposition:
        return all_free_energies, data
    else:
        return all_free_energies


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Analyzer to unbias radially symmetric restraints.')
    parser.add_argument('-s', '--store', metavar='store', type=str, help='Storage directory for NetCDF data files.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False)
    parser.add_argument('--energy-cutoff', metavar='energy_cutoff', type=float,
                        default=None, help='Energy cutoff in kT units.')
    parser.add_argument('--distance-cutoff', metavar='distance_cutoff', type=float,
                        default=None, help='Distance cutoff in Angstroms.')
    args = parser.parse_args()

    if args.verbose:
        utils.config_root_logger(verbose=True)
    else:
        utils.config_root_logger(verbose=False)

    if args.energy_cutoff:
        energy_cutoff = args.energy_cutoff
    else:
        energy_cutoff = None

    if args.distance_cutoff:
        distance_cutoff = args.distance_cutoff * unit.angstroms
    else:
        distance_cutoff = None

    analyze_directory(args.store, energy_cutoffs=energy_cutoff, distance_cutoffs=distance_cutoff)
