#!/usr/bin/env python

"""Functions for the efficiency statistical analysis."""


# =============================================================================
# GLOBAL IMPORTS
# =============================================================================

import copy

import numpy as np
import scipy as sp
import arch.bootstrap


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def discard_initial_zeros(free_energy_A, free_energy_B):
    """Discard the initial data points for which there's no data.

    Some submissions don't have estimates for the first part of
    the trajectory. This function removes the initial data points
    with mean free energy 0.0.

    Parameters
    ----------
    free_energy_A : numpy.ndarray
        free_energy_A[r][c] is the r-th replicate of the free energy estimate
        computed by method A at the c-th computational cost.
    free_energy_B : numpy.ndarray
        free_energy_B[r][c] is the r-th replicate of the free energy estimate
        computed by method B at the c-th computational cost.

    Returns
    -------
    free_energy_A : numpy.ndarray
    free_energy_B : numpy.ndarray
        A copy of the free energy trajectories after discarding the
        initial computational costs without data.

    """
    # Compute the mean free energies of the two methods.
    mean_c_A = free_energy_A.mean(axis=0)
    mean_c_B = free_energy_B.mean(axis=0)

    # Find the first indices for which both methods
    # have a non-zero free energy estimate.
    first_nonzero_idx_A = np.nonzero(mean_c_A)[0][0]
    first_nonzero_idx_B = np.nonzero(mean_c_B)[0][0]
    first_nonzero_idx = max(first_nonzero_idx_A, first_nonzero_idx_B)

    # Discard the initial computational costs.
    free_energy_A = copy.deepcopy(free_energy_A[:,first_nonzero_idx:])
    free_energy_B = copy.deepcopy(free_energy_B[:,first_nonzero_idx:])

    return free_energy_A, free_energy_B


# =============================================================================
# MAIN CLASS
# =============================================================================

class EfficiencyAnalysis:
    """Utility class for the calculation of relative efficiency statistics.

    Parameters
    ----------
    free_energy_A : numpy.ndarray
        free_energy_A[r][c] is the r-th replicate of the free energy estimate
        computed by method A at the c-th computational cost.
    free_energy_B : numpy.ndarray
        free_energy_B[r][c] is the r-th replicate of the free energy estimate
        computed by method B at the c-th computational cost.
    asymptotic_free_energy_A : float, optional
        If given, this will be used as the asymptotic free energy of
        method A to compute the bias. Otherwise, this is estimated
        as free_energy_A.mean(0)[-1].
    asymptotic_free_energy_B : float, optional
        If given, this will be used as the asymptotic free energy of
        method B to compute the bias. Otherwise, this is estimated
        as free_energy_B.mean(0)[-1].
    model : 'normal' or None
        If 'normal', the distribution of free_energy_X[:,c] at a fixed
        computational cost will be assumed to be normal, and parametric
        bootstrapping will be used. This has the disadvantage of ignoring
        eventual correlations between consecutive data points of the
        free energy trajectories, but it may model the tails of the
        distribution better when the number of independent trajectories
        is low.

    """

    def __init__(
            self, free_energy_A, free_energy_B,
            asymptotic_free_energy_A=None,
            asymptotic_free_energy_B=None,
            model=None,
    ):
        # Check that all means and stds have the same number of data points.
        shapes = {x.shape for x in [free_energy_A, free_energy_B]}
        if len(shapes) != 1:
            raise ValueError('free_energy_A and free_energy_B must have the same shape')

        if model not in {None, self._NORMAL_MODEL}:
            raise ValueError('model must be None or {}'.format(self._NORMAL_MODEL))

        self._free_energy_A = free_energy_A
        self._free_energy_B = free_energy_B
        self._asymptotic_free_energy_A = asymptotic_free_energy_A
        self._asymptotic_free_energy_B = asymptotic_free_energy_B
        self._model = model

    @property
    def n_replicates_A(self):
        """The number of replicate free energy trajectories for method A."""
        return self._free_energy_A.shape[0]

    @property
    def n_replicates_B(self):
        """The number of replicate free energy trajectories for method B."""
        return self._free_energy_B.shape[0]

    @property
    def asymptotic_free_energy_A(self):
        """The asymptotic free energy for method A."""
        return _estimate_asymptotic_free_energy(
            self._free_energy_A, self._asymptotic_free_energy_A)

    @property
    def asymptotic_free_energy_B(self):
        """The asymptotic free energy for method B."""
        return _estimate_asymptotic_free_energy(
            self._free_energy_B, self._asymptotic_free_energy_B)

    @property
    def mean_c_A(self):
        """mean_c_A[c] is the mean of the free energy estimates computed
        by method A at the c-th computational cost."""
        return self._free_energy_A.mean(axis=0)

    @property
    def mean_c_B(self):
        """mean_c_B[c] is the mean of the free energy estimates computed
        by method B at the c-th computational cost."""
        return self._free_energy_B.mean(axis=0)

    @property
    def std_c_A(self):
        """std_c_A[c] is the standard deviation of the free energy estimate
        computed by method A at the c-th computational cost."""
        if self._model == self._NORMAL_MODEL:
            return _normal_unbiased_std(self._free_energy_A)
        return self._free_energy_A.std(axis=0, ddof=1)

    @property
    def std_c_B(self):
        """std_c_B[c] is the standard deviation of the free energy estimate
        computed by method B at the c-th computational cost."""
        if self._model == self._NORMAL_MODEL:
            return _normal_unbiased_std(self._free_energy_B)
        return self._free_energy_B.std(axis=0, ddof=1)

    @property
    def var_c_A(self):
        """var_c_A[c] is the variance of the free energy estimate
        computed by method A at the c-th computational cost."""
        return self._free_energy_A.var(axis=0, ddof=1)

    @property
    def var_c_B(self):
        """var_c_B[c] is the variance of the free energy estimate
        computed by method B at the c-th computational cost."""
        return self._free_energy_B.var(axis=0, ddof=1)

    @property
    def bias_c_A(self):
        """bias_c_A[c] is the bias of the free energy estimates computed
        by method A at the c-th computational cost."""
        return _bias(self.mean_c_A, self.asymptotic_free_energy_A)

    @property
    def bias_c_B(self):
        """bias_c_B[c] is the bias of the free energy estimates computed
        by method B at the c-th computational cost."""
        return _bias(self.mean_c_B, self.asymptotic_free_energy_B)

    def compute_std_relative_efficiency(
            self,
            confidence_interval=None,
            n_bootstrap_samples=10000,
    ):
        std_relative_efficiency = _std_relative_efficiency(
            self._free_energy_A, self._free_energy_B)

        if self._model == self._NORMAL_MODEL:
            bootstrap_func = _generate_normal_std_rel_eff_sample_arch
            sampling = 'parametric'
        else:
            bootstrap_func = _std_relative_efficiency
            sampling = 'nonparametric'

        if confidence_interval is not None:
            ci = self._compute_rel_eff_ci(
                bootstrap_func, sampling, confidence_interval,
                n_bootstrap_samples, include_asymptotic=False
            )
            return std_relative_efficiency, ci
        return std_relative_efficiency

    def compute_abs_bias_relative_efficiency(
            self,
            confidence_interval=None,
            n_bootstrap_samples=10000,
    ):
        abs_bias_relative_efficiency = _abs_bias_relative_efficiency(
            self._free_energy_A, self._free_energy_B,
            self.asymptotic_free_energy_A,
            self.asymptotic_free_energy_B
        )

        if self._model == self._NORMAL_MODEL:
            bootstrap_func = _generate_normal_abs_bias_rel_eff_sample_arch
            sampling = 'parametric'
        else:
            bootstrap_func = _abs_bias_relative_efficiency
            sampling = 'nonparametric'

        if confidence_interval is not None:
            ci = self._compute_rel_eff_ci(
                bootstrap_func, sampling, confidence_interval,
                n_bootstrap_samples, include_asymptotic=True
            )
            return abs_bias_relative_efficiency, ci
        return abs_bias_relative_efficiency

    def compute_rmse_relative_efficiency(
            self,
            confidence_interval=None,
            n_bootstrap_samples=10000,
    ):
        rmse_bias_relative_efficiency = _rmse_relative_efficiency(
            self._free_energy_A, self._free_energy_B,
            self.asymptotic_free_energy_A,
            self.asymptotic_free_energy_B
        )

        if self._model == self._NORMAL_MODEL:
            bootstrap_func = _generate_normal_rmse_rel_eff_sample_arch
            sampling = 'parametric'
        else:
            bootstrap_func = _rmse_relative_efficiency
            sampling = 'nonparametric'

        if confidence_interval is not None:
            ci = self._compute_rel_eff_ci(
                bootstrap_func, sampling, confidence_interval,
                n_bootstrap_samples, include_asymptotic=True
            )
            return rmse_bias_relative_efficiency, ci
        return rmse_bias_relative_efficiency

    def _compute_rel_eff_ci(
            self, arch_func, sampling, confidence_interval,
            n_bootstrap_samples, include_asymptotic
    ):
        """Shortcut to compute a CI with arch.bootstrap."""
        bs = _IIDBootstrapNotEqual(self._free_energy_A, self._free_energy_B)

        # Configure extra keyword arguments for arch_func.
        kwargs = {}
        if self._model == self._NORMAL_MODEL:
            kwargs['params_cache'] = {}
        if include_asymptotic:
            kwargs['asymptotic_free_energy_A'] = self.asymptotic_free_energy_A
            kwargs['asymptotic_free_energy_B'] = self.asymptotic_free_energy_B

        if len(kwargs) == 0:
            kwargs = None

        ci = bs.conf_int(
            arch_func, reps=n_bootstrap_samples, method='bca',
            size=confidence_interval, sampling=sampling,
            extra_kwargs=kwargs
        )
        # Convert shape from (2,1) to (2,)
        assert ci.shape == (2, 1)
        return np.array([ci[0][0], ci[1][0]])

    # Class constants
    _NORMAL_MODEL = 'normal'


# =============================================================================
# Basic statistics utilities
# =============================================================================

def _normal_unbiased_std(data):
    """Return the unbiased estimate of the standard deviation assuming normal distribution.

    The sqrt of the Bessel-corrected variance is a biased estimate of
    the std. If the underlying data has a normal distribution, we have
    an analytical expression for the unbiased estimate.

    Parameters
    ----------
    data : np.ndarray
        data[i][j] is the i-th replicate of the j-th measurement (e.g.,
        the i-th replicate of the free energy trajectory at the j-th
        computational cost).

    Returns
    -------
    unbiased_std : np.ndarray
        unbiased_std[j] is the unbiased estimate of the standard
        deviation for the j-th measurement.

    See Also
    --------
    http://web.eecs.umich.edu/~fessler/papers/files/tr/stderr.pdf.

    """
    n_replicates = data.shape[0]

    bessel_std = np.std(data, ddof=1, axis=0)

    # Compute the factor used to unbias the Bessel-corrected sample std
    # Use the log gamma function for numerical stability.
    loggamma1 = sp.special.gammaln((n_replicates-1)/2)
    loggamma2 = sp.special.gammaln(n_replicates/2)
    k_n = np.sqrt((n_replicates-1)/2) * np.exp(loggamma1 - loggamma2)

    return k_n * bessel_std


def _bias(mean_c, asymptotic_free_energy=None):
    """Compute the bias from the mean."""
    asymptotic_free_energy = _estimate_asymptotic_free_energy(
        mean_c, asymptotic_free_energy)
    return mean_c - asymptotic_free_energy


def _estimate_asymptotic_free_energy(data, asymptotic_free_energy=None):
    """Shortcut to estimate the asymptotic free energy as mean_DG[-1] if asymptotic_free_energy is None.

    Parameters
    ----------
    data : numpy.ndarray
        If 2D this is handled as the free energy trajectories. If 1D
        this is considered the mean free energy.

    """
    if asymptotic_free_energy is None:
        if len(data.shape) == 2:
            data = data.mean(axis=0)
        return data[-1]
    return asymptotic_free_energy


# =============================================================================
# Relative efficiency definitions
# =============================================================================

def _relative_efficiency(stats_A, stats_B):
    """Encapsulate the definition of relative efficiency.

    Parameters
    ----------
    stats_A : np.ndarray
    stats_B : np.ndarray
        If 1D, stats_A[c] is the statistic at the c-th computational cost.
        If 2D, stats_A[r][c] is the r-th bootstrap sample of the statistic
        at the c-th computational cost.

    """
    assert stats_A.shape == stats_B.shape
    # Check if this is a bootstrapped version or not.
    if len(stats_A.shape) == 1:
        axis = None  # sum
        axis = -1  # trapz
    else:
        axis = 1

    # sum_A = np.sum(stats_A, axis=axis)
    # sum_B = np.sum(stats_B, axis=axis)
    sum_A = sp.integrate.trapz(stats_A, axis=axis)
    sum_B = sp.integrate.trapz(stats_B, axis=axis)
    return np.log10(sum_A / sum_B)


def _std_relative_efficiency(free_energy_A, free_energy_B, model=None):
    """Shortcut to compute the standard deviation relative efficiency from the free energy trajectories."""
    if model == 'normal':
        std_c_A = _normal_unbiased_std(free_energy_A)
        std_c_B = _normal_unbiased_std(free_energy_B)
    else:
        assert model is None
        std_c_A = np.std(free_energy_A, axis=0, ddof=1)
        std_c_B = np.std(free_energy_B, axis=0, ddof=1)
    return _relative_efficiency(std_c_A, std_c_B)


def _abs_bias_relative_efficiency_from_params(
        mean_c_A, mean_c_B,
        asymptotic_free_energy_A=None,
        asymptotic_free_energy_B=None
):
    """Shortcut to compute the absolute bias relative efficiency from mean and asymptotic free energy."""
    bias_c_A = _bias(mean_c_A, asymptotic_free_energy_A)
    bias_c_B = _bias(mean_c_B, asymptotic_free_energy_B)
    return _relative_efficiency(np.abs(bias_c_A), np.abs(bias_c_B))


def _abs_bias_relative_efficiency(
        free_energy_A, free_energy_B,
        asymptotic_free_energy_A=None,
        asymptotic_free_energy_B=None
):
    """Shortcut to compute the absolute bias relative efficiency from the free energy trajectories."""
    mean_c_A = free_energy_A.mean(axis=0)
    mean_c_B = free_energy_B.mean(axis=0)
    return _abs_bias_relative_efficiency_from_params(
        mean_c_A, mean_c_B, asymptotic_free_energy_A, asymptotic_free_energy_B)


def _rmse_relative_efficiency_from_params(
        mean_c_A, mean_c_B, var_c_A, var_c_B,
        asymptotic_free_energy_A=None,
        asymptotic_free_energy_B=None
):
    """Shortcut to compute the RMSE relative efficiency from mean, var, and asymptotic free energy."""
    bias_c_A = _bias(mean_c_A, asymptotic_free_energy_A)
    bias_c_B = _bias(mean_c_B, asymptotic_free_energy_B)
    rmse_c_A = np.sqrt(var_c_A + bias_c_A**2)
    rmse_c_B = np.sqrt(var_c_B + bias_c_B**2)
    return _relative_efficiency(rmse_c_A, rmse_c_B)


def _rmse_relative_efficiency(
        free_energy_A, free_energy_B,
        asymptotic_free_energy_A=None,
        asymptotic_free_energy_B=None
):
    """Shortcut to compute the RMSE relative efficiency from the free energy trajectories."""
    asymptotic_free_energy_A = _estimate_asymptotic_free_energy(
        free_energy_A, asymptotic_free_energy_A)
    asymptotic_free_energy_B = _estimate_asymptotic_free_energy(
        free_energy_B, asymptotic_free_energy_B)
    rmse_c_A = np.sqrt(np.mean((free_energy_A - asymptotic_free_energy_A)**2, axis=0))
    rmse_c_B = np.sqrt(np.mean((free_energy_B - asymptotic_free_energy_B)**2, axis=0))
    return _relative_efficiency(rmse_c_A, rmse_c_B)


# =============================================================================
# Parametric bootstrap sampling with normal statistics
# =============================================================================

def _generate_normal_std_sample(
        std_c, n_replicates,
        n_bootstrap_samples=None
):
    """Generate a bootstrap sample for the standard deviation assuming the
    data is normally-distributed.

    In this function, free_energy[:, c] at a fixed computational cost is
    assumed to be normally distributed. Under this assumption the standard
    deviation is chi distributed with n_replicates-1 degrees of freedom,
    which allows us to generate a bootstrap sample for the standard
    deviation as a function of the computational cost from a parametric
    distribution.

    Parameters
    ----------
    std_c : numpy.ndarray
        std_c[c] is the standard deviation of the free energy estimate
        at the c-th computational cost.
    n_replicates : int
        The number of replicates used to compute the standard deviations.
    n_bootstrap_samples : int or None, optional
        If not None, multiple bootstrap samples are generated.

    Returns
    -------
    bootstrap_std : np.ndarray
        If n_bootstrap_samples is None, bootstrap_std[c] is the bootstrap
        sample of the standard deviation at the c-th computational cost.
        Otherwise, bootstrap_std[r][c] is the r-th bootstrap sample of the
        standard deviation at the c-th computational cost.

    """
    # Sample from the chi distribution with the correct scale.
    df = n_replicates - 1
    n_costs = len(std_c)
    if n_bootstrap_samples is None:
        size = n_costs
    else:
        size = (n_bootstrap_samples, n_costs)
    return sp.stats.chi.rvs(df=df, scale=std_c/np.sqrt(df), size=size)


def _generate_normal_std_rel_eff_sample(
        std_c_A, std_c_B, n_replicates,
        n_bootstrap_samples=None
):
    """Generate a bootstrap sample for the standard deviation relative efficiency
    assuming the data is normally-distributed.

    See Also
    --------
    _generate_normal_std_sample

    Parameters
    ----------
    std_c_A : numpy.ndarray
        std_c_A[c] is the standard deviation of the free energy estimate
        computed by method A at the c-th computational cost.
    std_c_A : numpy.ndarray
        std_c_B[c] is the standard deviation of the free energy estimate
        computed by method B at the c-th computational cost.
    n_replicates : int
        The number of replicates used to compute the standard deviations.
    n_bootstrap_samples : int or None, optional
        If not None, multiple bootstrap samples are generated.

    Returns
    -------
    bootstrap_std_rel_eff : np.ndarray
        If n_bootstrap_samples is None, bootstrap_std_rel_eff[c] is the
        bootstrap sample of the standard deviation relative efficiency at
        the c-th computational cost. Otherwise, bootstrap_std_rel_eff[r][c]
        is the r-th bootstrap sample of the standard deviation relative
        efficiency at the c-th computational cost.

    """
    bootstrap_std_c_A = _generate_normal_std_sample(std_c_A, n_replicates, n_bootstrap_samples)
    bootstrap_std_c_B = _generate_normal_std_sample(std_c_B, n_replicates, n_bootstrap_samples)
    # Compute the relative efficiency.
    return _relative_efficiency(bootstrap_std_c_A, bootstrap_std_c_B)


def _generate_normal_abs_bias_sample(
        mean_c, std_c, n_replicates,
        asymptotic_free_energy=None,
        n_bootstrap_samples=None
):
    """Generate a bootstrap sample for the absolute bias assuming the
    data to be normally-distributed.

    In this function, free_energy[:, c] at a fixed computational cost is
    assumed to be normally distributed. Under this assumption the sample
    mean is t-distributed with n_replicates-1 degrees of freedom, which
    allows us to generate a bootstrap sample for the absolute bias as a
    function of the computational cost from a parametric distribution.

    Parameters
    ----------
    mean_c : numpy.ndarray
        mean_c[c] is the mean of the free energy estimates computed
        at the c-th computational cost.
    std_c : numpy.ndarray
        std_c[c] is the standard deviation of the free energy estimate
        computed at the c-th computational cost.
    n_replicates : int
        The number of replicates used to compute the means and standard
        deviations.
    asymptotic_free_energy : float, optional
        If given, this will be used as the asymptotic free energy of
        to compute the bias. Otherwise, this is estimated as mean_c[-1].
    n_bootstrap_samples : int or None, optional
        If not None, multiple bootstrap samples are generated.

    Returns
    -------
    bootstrap_abs_bias : np.ndarray
        If n_bootstrap_samples is None, bootstrap_abs_bias[c] is the bootstrap
        sample of the absolute bias at the c-th computational cost. Otherwise,
        bootstrap_abs_bias[r][c] is the r-th bootstrap sample of the absolute
        bias at the c-th computational cost.

    """
    # Sample mean the t distribution with the correct scale and mean.
    df = n_replicates - 1
    bias_c = _bias(mean_c, asymptotic_free_energy)
    n_costs = len(std_c)
    if n_bootstrap_samples is None:
        size = n_costs
    else:
        size = (n_bootstrap_samples, n_costs)
    bootstrap_abs_bias_c = sp.stats.t.rvs(df=df, scale=std_c/np.sqrt(n_replicates), size=size) + bias_c
    return np.abs(bootstrap_abs_bias_c)


def _generate_normal_abs_bias_rel_eff_sample(
        mean_c_A, mean_c_B, std_c_A, std_c_B, n_replicates,
        asymptotic_free_energy_A=None, asymptotic_free_energy_B=None,
        n_bootstrap_samples=None
):
    """Generate a bootstrap sample for the absolute bias relative efficiency
    assuming the data to be normally distributed.

    See Also
    --------
    _generate_normal_abs_bias_sample

    Parameters
    ----------
    mean_c_A : numpy.ndarray
        mean_c_A[c] is the mean of the free energy estimates computed
        by method A at the c-th computational cost.
    mean_c_B : numpy.ndarray
        mean_c_B[c] is the mean of the free energy estimates computed
        by method B at the c-th computational cost.
    std_c_A : numpy.ndarray
        std_c_A[c] is the standard deviation of the free energy estimate
        computed by method A at the c-th computational cost.
    std_c_B : numpy.ndarray
        std_c_B[c] is the standard deviation of the free energy estimate
        computed by method B at the c-th computational cost.
    n_replicates : int
        The number of replicates used to compute the means and standard
        deviations.
    asymptotic_free_energy_A : float, optional
        If given, this will be used as the asymptotic free energy of
        method A to compute the bias. Otherwise, this is estimated
        as mean_c_A[-1].
    asymptotic_free_energy_B : float, optional
        If given, this will be used as the asymptotic free energy of
        method B to compute the bias. Otherwise, this is estimated
        as mean_c_B[-1].
    n_bootstrap_samples : int or None, optional
        If not None, multiple bootstrap samples are generated.

    Returns
    -------
    bootstrap_abs_bias_rel_eff : np.ndarray
        If n_bootstrap_samples is None, bootstrap_abs_bias_rel_eff[c] is
        the bootstrap sample of the absolute bias relative efficiency at the
        c-th computational cost. Otherwise, bootstrap_abs_bias_rel_eff[r][c]
        is the r-th bootstrap sample of the absolute bias relative efficiency
        at the c-th computational cost.

    """
    bootstrap_bias_c_A = _generate_normal_abs_bias_sample(
        mean_c_A, std_c_A, n_replicates,
        asymptotic_free_energy=asymptotic_free_energy_A,
        n_bootstrap_samples=n_bootstrap_samples
    )
    bootstrap_bias_c_B = _generate_normal_abs_bias_sample(
        mean_c_B, std_c_B, n_replicates,
        asymptotic_free_energy=asymptotic_free_energy_B,
        n_bootstrap_samples=n_bootstrap_samples
    )
    # Compute the relative efficiency.
    return _relative_efficiency(bootstrap_bias_c_A, bootstrap_bias_c_B)


def _generate_normal_rmse_sample(
        mean_c, std_c, n_replicates,
        asymptotic_free_energy=None,
        n_bootstrap_samples=None
):
    """Generate a bootstrap sample for the RMSE assuming the data to be
    normally distributed.

    In this function, free_energy[:, c] at a fixed computational cost is
    assumed to be normally distributed. Under this assumption the sample
    mean is t-distributed and the standard deviation is chi-distributed
    with n_replicates-1 degrees of freedom, which allows us to generate
    a bootstrap sample for the absolute bias as a function of the computational
    cost from a parametric distribution.

    Parameters
    ----------
    mean_c : numpy.ndarray
        mean_c[c] is the mean of the free energy estimates computed
        at the c-th computational cost.
    std_c : numpy.ndarray
        std_c[c] is the standard deviation of the free energy estimate
        computed at the c-th computational cost.
    n_replicates : int
        The number of replicates used to compute the means and standard
        deviations.
    asymptotic_free_energy : float, optional
        If given, this will be used as the asymptotic free energy of
        to compute the bias. Otherwise, this is estimated as mean_c[-1].
    n_bootstrap_samples : int or None, optional
        If not None, multiple bootstrap samples are generated.

    Returns
    -------
    bootstrap_rmse : np.ndarray
        If n_bootstrap_samples is None, bootstrap_rmse[c] is the bootstrap
        sample of the RMSE at the c-th computational cost. Otherwise,
        bootstrap_rmse[r][c] is the r-th bootstrap sample of the RMSE at
        the c-th computational cost.

    """
    bootstrap_std_c = _generate_normal_std_sample(
        std_c, n_replicates, n_bootstrap_samples)
    bootstrap_abs_bias_c = _generate_normal_abs_bias_sample(
        mean_c, std_c, n_replicates, asymptotic_free_energy, n_bootstrap_samples)
    return np.sqrt(bootstrap_std_c**2 + bootstrap_abs_bias_c**2)


def _generate_normal_rmse_rel_eff_sample(
        mean_c_A, mean_c_B, std_c_A, std_c_B, n_replicates,
        asymptotic_free_energy_A=None, asymptotic_free_energy_B=None,
        n_bootstrap_samples=None
):
    """Generate a bootstrap sample for the RMSE relative efficiency
    assuming the data to be normally distributed.

    See Also
    --------
    _generate_normal_rmse_sample

    Parameters
    ----------
    mean_c_A : numpy.ndarray
        mean_c_A[c] is the mean of the free energy estimates computed
        by method A at the c-th computational cost.
    mean_c_B : numpy.ndarray
        mean_c_B[c] is the mean of the free energy estimates computed
        by method B at the c-th computational cost.
    std_c_A : numpy.ndarray
        std_c_A[c] is the standard deviation of the free energy estimate
        computed by method A at the c-th computational cost.
    std_c_A : numpy.ndarray
        std_c_B[c] is the standard deviation of the free energy estimate
        computed by method B at the c-th computational cost.
    n_replicates : int
        The number of replicates used to compute the means and standard
        deviations.
    asymptotic_free_energy_A : float, optional
        If given, this will be used as the asymptotic free energy of
        method A to compute the bias. Otherwise, this is estimated
        as mean_c_A[-1].
    asymptotic_free_energy_B : float, optional
        If given, this will be used as the asymptotic free energy of
        method B to compute the bias. Otherwise, this is estimated
        as mean_c_B[-1].
    n_bootstrap_samples : int or None, optional
        If not None, multiple bootstrap samples are generated.

    Returns
    -------
    bootstrap_rmse_rel_eff : np.ndarray
        If n_bootstrap_samples is None, bootstrap_rmse_rel_eff[c] is the
        bootstrap sample of the RMSE relative efficiency at the c-th
        computational cost. Otherwise, bootstrap_rmse_rel_eff[r][c] is
        the r-th bootstrap sample of the RMSE relative efficiency at
        the c-th computational cost.

    """
    # Generate samples for the std and bias.
    bootstrap_rmse_c_A = _generate_normal_rmse_sample(
        mean_c_A, std_c_A, n_replicates, asymptotic_free_energy_A, n_bootstrap_samples)
    bootstrap_rmse_c_B = _generate_normal_rmse_sample(
        mean_c_B, std_c_B, n_replicates, asymptotic_free_energy_B, n_bootstrap_samples)
    return _relative_efficiency(bootstrap_rmse_c_A, bootstrap_rmse_c_B)


# =============================================================================
# Wrappers of bootstrap sampling functions for arch.bootstrap.IIDBootstrap
# =============================================================================

class _IIDBootstrapNotEqual(arch.bootstrap.IIDBootstrap):
    """A bootstrap facility class that avoid generating bootstrap sampling
    concentrating all the distribution on a single data point."""

    def bootstrap(self, reps):
        for _ in range(reps):
            indices = np.asarray(self.update_indices())
            # Regenerate indices until there is at least one that is different
            # to avoid generating trajectories with std == 0.0.
            while np.allclose(indices, indices[1]):
                indices = np.asarray(self.update_indices())
            self._index = indices
            yield self._resample()


def _cache_params_arch(
        free_energy_A, free_energy_B, params_cache, compute_mean=True,
):
    """Utility function to handle the cache of the free energy mean and std."""
    if params_cache is None:
        # Generate parameters just for this call.
        params_cache = {}

    if isinstance(params_cache, dict) and len(params_cache) == 0:
        # Cache number of replicates.
        assert free_energy_A.shape == free_energy_B.shape
        params_cache['n_replicates'] = free_energy_A.shape[0]

        for i, free_energy in enumerate([free_energy_A, free_energy_B]):
            suffix = 'A' if i == 0 else 'B'

            # Cache standard deviation.
            params_cache['std_c_' + suffix] = _normal_unbiased_std(free_energy)
            # Cache mean bias if requested.
            if compute_mean:
                params_cache['mean_c_' + suffix] = free_energy.mean(axis=0)

    return params_cache


def _generate_normal_std_rel_eff_sample_arch(
        free_energy_A, free_energy_B,
        params=None, state=None,
        params_cache=None
):
    """Wraps around _generate_normal_std_rel_eff_sample for use with arch.bootstrap."""
    params_cache = _cache_params_arch(free_energy_A, free_energy_B,
                                      params_cache, compute_mean=False)
    if params is None:
        return _relative_efficiency(params_cache['std_c_A'], params_cache['std_c_B'])
    return _generate_normal_std_rel_eff_sample(**params_cache)


def _generate_normal_abs_bias_rel_eff_sample_arch(
        free_energy_A, free_energy_B, params=None, state=None,
        asymptotic_free_energy_A=None, asymptotic_free_energy_B=None,
        params_cache=None
):
    """Wraps around _generate_normal_abs_bias_rel_eff_sample for use with arch.bootstrap."""
    params_cache = _cache_params_arch(free_energy_A, free_energy_B,
                                      params_cache, compute_mean=True)

    if params is None:
        return _abs_bias_relative_efficiency_from_params(
            params_cache['mean_c_A'], params_cache['mean_c_B'],
            asymptotic_free_energy_A, asymptotic_free_energy_B
        )
    return _generate_normal_abs_bias_rel_eff_sample(
        asymptotic_free_energy_A=asymptotic_free_energy_A,
        asymptotic_free_energy_B=asymptotic_free_energy_B,
        **params_cache
    )


def _generate_normal_rmse_rel_eff_sample_arch(
        free_energy_A, free_energy_B, params=None, state=None,
        asymptotic_free_energy_A=None, asymptotic_free_energy_B=None,
        params_cache=None
):
    """Wraps around _generate_normal_abs_bias_rel_eff_sample for use with arch.bootstrap."""
    params_cache = _cache_params_arch(free_energy_A, free_energy_B,
                                      params_cache, compute_mean=True)
    if params is None:
        return _rmse_relative_efficiency_from_params(
            params_cache['mean_c_A'], params_cache['mean_c_B'],
            params_cache['std_c_A']**2, params_cache['std_c_B']**2,
            asymptotic_free_energy_A, asymptotic_free_energy_B,
        )

    return _generate_normal_rmse_rel_eff_sample(
        asymptotic_free_energy_A=asymptotic_free_energy_A,
        asymptotic_free_energy_B=asymptotic_free_energy_B,
        **params_cache
    )
