#!/usr/bin/env python

"""Functions for the efficiency statistical analysis."""


# =============================================================================
# GLOBAL IMPORTS
# =============================================================================

import copy

import numpy as np
import scipy as sp


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

class RelativeEfficiency:
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

    """

    def __init__(
            self, free_energy_A, free_energy_B,
            asymptotic_free_energy_A=None,
            asymptotic_free_energy_B=None
    ):
        # Check that all means and stds have the same number of data points.
        shapes = {x.shape for x in [free_energy_A, free_energy_B]}
        if len(shapes) != 1:
            raise ValueError('free_energy_A and free_energy_B must have the same shape')

        self._free_energy_A = free_energy_A
        self._free_energy_B = free_energy_B
        self.asymptotic_free_energy_A = asymptotic_free_energy_A
        self.asymptotic_free_energy_B = asymptotic_free_energy_B

    @property
    def n_replicates(self):
        """The number of replicate free energy trajectories."""
        return self._free_energy_A.shape[0]

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
        return _compute_normal_unbiased_std(self._free_energy_A)

    @property
    def std_c_B(self):
        """std_c_B[c] is the standard deviation of the free energy estimate
        computed by method B at the c-th computational cost."""
        return _compute_normal_unbiased_std(self._free_energy_B)

    @property
    def bias_c_A(self):
        """bias_c_A[c] is the bias of the free energy estimates computed
        by method A at the c-th computational cost."""
        mean_c_A = self.mean_c_A
        if self.asymptotic_free_energy_A is None:
            return mean_c_A - mean_c_A[-1]
        return mean_c_A - self.asymptotic_free_energy_A

    @property
    def bias_c_B(self):
        """bias_c_B[c] is the bias of the free energy estimates computed
        by method A at the c-th computational cost."""
        mean_c_B = self.mean_c_B
        if self.asymptotic_free_energy_B is None:
            return mean_c_B - mean_c_B[-1]
        return mean_c_B - self.asymptotic_free_energy_B

    def compute_mean_relative_efficiencies(
            self,
            weighted=True,
            n_bootstrap_samples=10000,
            confidence_interval=None,
    ):
        """Compute the arithmetic mean relative std, absolute bias, and RMSD efficiency for the data.

        Parameters
        ----------
        weighted : bool, optional
            If True, the mean is weighted by the inverse variance of the
            efficiency.
        confidence_interval : float or None
            The inter-percentile range used for the confidence interval. If
            None, no confidence interval is returned. Default is None.
        n_bootstrap_samples : int, optional
            The number of bootstrap samples per computational cost used to
            compute the mean weights and the mean relative efficiencies confidence
            intervals. Default is 10000.

        Returns
        -------
        std_mean_efficiency : float
            The mean std relative efficiency.
        abs_bias_mean_efficiency : float
            The mean absolute bias relative efficiency.
        rmse_mean_efficiency : float
            The mean RMSE relative efficiency.
        std_mean_efficiency_CI : List[float], optional
            If confidence_interval is not None, this is the pair (low_bound, high_bound)
            defining the bootstrap confidence interval of the std mean relative
            efficiency computed with bootstrapping.
        abs_bias_mean_efficiency_CI : List[float], optional
            If confidence_interval is not None, this is the pair (low_bound, high_bound)
            defining the bootstrap confidence interval of the absolute bias mean relative
            efficiency computed with bootstrapping.
        rmse_mean_efficiency_CI : List[float], optional
            If confidence_interval is not None, this is the pair (low_bound, high_bound)
            defining the bootstrap confidence interval of the RMSE mean relative
            efficiency computed with bootstrapping.

        """
        # Compute relative efficiencies.
        relative_efficiencies = self.compute_relative_efficiencies()

        # Generate mean weights with bootstrapping if required.
        if weighted:
            statistics_weights = self.compute_relative_efficiencies_mean_weights(
                n_bootstrap_samples)
        else:
            statistics_weights = [None for _ in relative_efficiencies]

        # Compute mean relative efficiencies.
        mean_efficiencies = []
        for rel_eff_c, weights in zip(relative_efficiencies, statistics_weights):
            mean_eff = np.average(rel_eff_c, weights=weights)
            mean_efficiencies.append(mean_eff)

        # Compute confidence intervals.
        confidence_intervals = []
        if confidence_interval is not None:
            from arch.bootstrap import IIDBootstrap
            bs = IIDBootstrap(self._free_energy_A, self._free_energy_B)

            for statistic_idx, func in enumerate([
                _generate_normal_std_mean_rel_eff_sample_arch,
                _generate_normal_abs_bias_mean_rel_eff_sample_arch,
                _generate_normal_rmse_mean_rel_eff_sample_arch,
            ]):
                if statistic_idx == 0:
                    extra_kwargs = {}
                else:
                    extra_kwargs = {
                        'asymptotic_free_energy_A': self.asymptotic_free_energy_A,
                        'asymptotic_free_energy_B': self.asymptotic_free_energy_B,
                    }

                ci = bs.conf_int(
                    func, reps=n_bootstrap_samples, method='basic',
                    size=confidence_interval, sampling='parametric',
                    extra_kwargs={
                        'params_cache': {},
                        'weights': statistics_weights[statistic_idx],
                        **extra_kwargs
                    },
                )
                confidence_intervals.append(ci)

        return mean_efficiencies + confidence_intervals

    def compute_relative_efficiencies(
            self,
            std_A=None, bias_A=None,
            std_B=None, bias_B=None,
            confidence_interval=None,
            n_bootstrap_samples=10000,
    ):
        """Compute std, absolute bias, and RMSE relative effiencies for all computational costs.

        Parameters
        ----------
        std_A : np.ndarray, optional
            This can be a 1D or 2D array of standard deviations for method A.
            If 1D, std_c[c] could be the standard deviation of the free energy
            estimate at the c-th computational cost. If 2D, std_c[i][c] is
            the i-th bootstrap sample for the c-th computational cost. If not
            given self.std_c_A is used.
        std_B : np.ndarray, optional
            1D or 2D array of standard deviation for method B. If not given
            self.std_c_B is used.
        bias_A : np.ndarray, optional
            1D or 2D array of biases for method A. If not given self.bias_c_A
            is used.
        bias_B : np.ndarray, optional
            1D or 2D array of biases for method B. If not given self.bias_c_B
            is used.
        confidence_interval : float or None
            The inter-percentile range used for the confidence interval. If
            None, no confidence interval is returned. Default is None.
        n_bootstrap_samples : int, optional
            The number of bootstrap samples per computational cost used to
            compute the mean weights and the mean relative efficiencies confidence
            intervals. Default is 10000.

        Returns
        -------
        std_efficiencies : np.ndarray
            The array of std relative efficiencies. The shape is the same as std.
        abs_bias_efficiencies : np.ndarray
            The array of absolute bias relative efficiencies. The shape is the
            same as std.
        rmse_efficiencies : np.ndarray
            The array of RMSE relative efficiencies. The shape is the same as std.
        std_efficiency_CI : List[float], optional
            If confidence_interval is not None, this is the pair (low_bound, high_bound)
            defining the bootstrap confidence interval of the std mean relative
            efficiency computed with bootstrapping.
        abs_bias_efficiency_CI : List[float], optional
            If confidence_interval is not None, this is the pair (low_bound, high_bound)
            defining the bootstrap confidence interval of the absolute bias mean relative
            efficiency computed with bootstrapping.
        rmse_efficiency_CI : List[float], optional
            If confidence_interval is not None, this is the pair (low_bound, high_bound)
            defining the bootstrap confidence interval of the RMSE mean relative
            efficiency computed with bootstrapping.

        """
        if std_A is None:
            std_A = self.std_c_A
        if bias_A is None:
            bias_A = self.bias_c_A
        if std_B is None:
            std_B = self.std_c_B
        if bias_B is None:
            bias_B = self.bias_c_B

        rmse_A = np.sqrt(std_A**2 + bias_A**2)
        rmse_B = np.sqrt(std_B**2 + bias_B**2)

        # Compute the relative efficiencies.
        relative_efficiencies = [
            _compute_relative_efficiency(std_A, std_B),
            _compute_relative_efficiency(np.abs(bias_A), np.abs(bias_B)),
            _compute_relative_efficiency(rmse_A, rmse_B)
        ]

        # Compute confidence intervals.
        confidence_intervals = []
        if confidence_interval is not None:
            from arch.bootstrap import IIDBootstrap
            bs = IIDBootstrap(self._free_energy_A, self._free_energy_B)

            for statistic_idx, func in enumerate([
                _generate_normal_std_rel_eff_sample_arch,
                _generate_normal_abs_bias_rel_eff_sample_arch,
                _generate_normal_rmse_rel_eff_sample_arch,
            ]):
                if statistic_idx == 0:
                    extra_kwargs = {}
                else:
                    extra_kwargs = {
                        'asymptotic_free_energy_A': self.asymptotic_free_energy_A,
                        'asymptotic_free_energy_B': self.asymptotic_free_energy_B,
                    }

                ci = bs.conf_int(
                    func, reps=n_bootstrap_samples, method='basic',
                    size=confidence_interval, sampling='parametric',
                    extra_kwargs={'params_cache': {}, **extra_kwargs},
                )
                confidence_intervals.append(ci)

        return relative_efficiencies + confidence_intervals

    def compute_relative_efficiencies_mean_weights(self, n_bootstrap_samples=10000):
        """
        Compute the weights to compute the mean relative efficiencies.

        Parameters
        ----------
        n_bootstrap_samples : int, optional
            The number of bootstrap samples per computational cost used to
            compute the mean weights and the mean relative efficiencies confidence
            intervals. Default is 10000.

        Returns
        -------
        std_efficiency_weights_c : numpy.ndarray
        abs_bias_efficiency_weights_c : numpy.ndarray
        rmse_efficiency_weights_c : numpy.ndarray
            X_efficiency_weights_c[c] is the mean weight for the c-th
            computational cost for the standard deviation, absolute bias
            and RMSE relative efficiencies.

        """
        # Generate relative efficiency bootstrap distributions.
        bootstrap_distributions = []
        for method in ['A', 'B']:
            std_c = getattr(self, 'std_c_' + method)
            mean_c = getattr(self, 'mean_c_' + method)
            asymptotic_free_energy = getattr(self, 'asymptotic_free_energy_' + method)

            # Generate std and abs bias bootstrap distributions.
            bootstrap_std_c = _generate_normal_std_sample(
                std_c, self.n_replicates, n_bootstrap_samples)
            bootstrap_abs_bias_c = _generate_normal_abs_bias_sample(
                mean_c, std_c, self.n_replicates, asymptotic_free_energy, n_bootstrap_samples)

            # Update args to pass to compute_relative_efficiencies().
            bootstrap_distributions.extend((bootstrap_std_c, bootstrap_abs_bias_c))

        # Compute relative efficiency bootstrap distributions.
        bootstrap_rel_effs = self.compute_relative_efficiencies(*bootstrap_distributions)

        # Compute weights.
        statistics_weights = []
        for bootstrap_rel_eff in bootstrap_rel_effs:
            statistics_weights.append(_compute_inverse_variance_weights(bootstrap_rel_eff))

        return statistics_weights

# =============================================================================
# Internal-usage utility functions
# =============================================================================

def _compute_normal_unbiased_std(data):
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


def _compute_bias(mean_c, asymptotic_free_energy=None):
    """Compute the bias from the mean."""
    if asymptotic_free_energy is None:
        asymptotic_free_energy = mean_c[-1]
    return mean_c - asymptotic_free_energy


def _compute_relative_efficiency(stats_c_A, stats_c_B):#, discard_zeroes=True):
    """Encapsulate the definition of relative efficiency."""
    # if discard_zeroes:
    #     non_zero_indices = ~np.isclose(stats_c_A, 0.0)
    #     stats_c_A = stats_c_A[non_zero_indices]
    #     stats_c_B = stats_c_B[non_zero_indices]
    #     non_zero_indices = ~np.isclose(stats_c_B, 0.0)
    #     stats_c_A = stats_c_A[non_zero_indices]
    #     stats_c_B = stats_c_B[non_zero_indices]
    # return stats_c_A / stats_c_B
    return stats_c_A - stats_c_B


def _compute_abs_bias_relative_efficiency(
        mean_c_A, mean_c_B,
        asymptotic_free_energy_A=None,
        asymptotic_free_energy_B=None
):
    bias_c_A = _compute_bias(mean_c_A, asymptotic_free_energy_A)
    bias_c_B = _compute_bias(mean_c_B, asymptotic_free_energy_B)
    return _compute_relative_efficiency(np.abs(bias_c_A), np.abs(bias_c_B))


def _compute_rmse_relative_efficiency(
        mean_c_A, mean_c_B, std_c_A, std_c_B,
        asymptotic_free_energy_A=None,
        asymptotic_free_energy_B=None
):
    bias_c_A = _compute_bias(mean_c_A, asymptotic_free_energy_A)
    bias_c_B = _compute_bias(mean_c_B, asymptotic_free_energy_B)
    rmse_c_A = np.sqrt(std_c_A**2 + bias_c_A**2)
    rmse_c_B = np.sqrt(std_c_B**2 + bias_c_B**2)
    return _compute_relative_efficiency(rmse_c_A, rmse_c_B)


def _compute_inverse_variance_weights(efficiency_distribution):
    """Compute the mean weights from the bootstrap distribution of the efficiency.

    Parameters
    ----------
    efficiency_distribution : numpy.ndarray
        efficiency_distribution[i][c] is the i-th bootstrap sample
        relative efficiency generated at the c-th computational cost.

    Returns
    -------
    efficiency_weights_c : numpy.ndarray
        efficiency_weights_c[c] is the mean weight for the c-th
        computational cost.

    """

    n_costs = efficiency_distribution.shape[1]
    efficiency_vars = np.array([np.var(efficiency_distribution[:,c], ddof=1) for c in range(n_costs)])

    # Compute the normalized weights.
    efficiency_weights = 1 / efficiency_vars
    efficiency_weights /= np.sum(efficiency_weights)
    return efficiency_weights


# =============================================================================
# Parametric bootstrap sampling
# =============================================================================

def _generate_normal_std_sample(std_c, n_replicates, n_bootstrap_samples=None):
    """Generate a bootstrap sample for the standard deviation.

    In this function, free_energy[c] at a fixed computational cost is
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


def _generate_normal_std_rel_eff_sample(std_c_A, std_c_B, n_replicates, n_bootstrap_samples=None):
    """Generate a bootstrap sample for the standard deviation relative efficiency.

    In this function, free_energy[c] at a fixed computational cost is
    assumed to be normally distributed. Under this assumption the standard
    deviation is chi distributed with n_replicates-1 degrees of freedom,
    which allows us to generate a bootstrap sample for the standard
    deviation as a function of the computational cost from a parametric
    distribution.

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
        If n_bootstrap_samples is None, bootstrap_std_rel_eff[c] is
        the bootstrap sample of the standard deviation relative efficiency
        at the c-th computational cost. Otherwise, bootstrap_std_rel_eff[r][c]
        is the r-th bootstrap sample of the standard deviation relative
        efficiency at the c-th computational cost.

    """
    bootstrap_std_c_A = _generate_normal_std_sample(std_c_A, n_replicates, n_bootstrap_samples)
    bootstrap_std_c_B = _generate_normal_std_sample(std_c_B, n_replicates, n_bootstrap_samples)
    # Compute the relative efficiency.
    return _compute_relative_efficiency(bootstrap_std_c_A, bootstrap_std_c_B)


def _generate_normal_std_mean_rel_eff_sample(*args, weights=None, **kwargs):
    """Generate a bootstrap sample for the standard deviation relative efficiency.

    Parameters
    ----------
    *args
        Positional arguments to forward to _generate_normal_std_rel_eff_sample.
    weights : numpy.ndarray or None, optional
        If given, weights[c] is the weight to be used in the weighted
        mean of the relative efficiency. If None, the average won't be
        weighted. Default is None.
    **kwargs
        Keyword arguments to forward to _generate_normal_std_rel_eff_sample.

    Returns
    -------
    bootstrap_std_mean_rel_eff : float
        bootstrap_std_mean_rel_eff is the bootstrap sample of the standard
        deviation mean relative efficiency.

    See Also
    --------
    _generate_normal_std_rel_eff_sample

    """
    bootstrap_std_rel_eff = _generate_normal_std_rel_eff_sample(*args, **kwargs)
    return np.average(bootstrap_std_rel_eff, weights=weights, axis=0)


def _generate_normal_abs_bias_sample(
        mean_c, std_c, n_replicates,
        asymptotic_free_energy=None,
        n_bootstrap_samples=None
):
    """Generate a bootstrap sample for the absolute bias.

    In this function, free_energy[c] at a fixed computational cost is
    assumed to be normally distributed. Under this assumption the sample
    mean is t-distributed with n_replicates-1 degrees of freedom,
    which allows us to generate a bootstrap sample for the absolute
    bias as a function of the computational cost from a parametric
    distribution.

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
    bias_c = _compute_bias(mean_c, asymptotic_free_energy)
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
    """Generate a bootstrap sample for the absolute bias relative efficiency.

    In this function, free_energy[c] at a fixed computational cost is
    assumed to be normally distributed. Under this assumption the sample
    mean is t-distributed with n_replicates-1 degrees of freedom,
    which allows us to generate a bootstrap sample for the absolute
    bias as a function of the computational cost from a parametric
    distribution.

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
        the bootstrap sample of the absolute bias relative efficiency at
        the c-th computational cost. Otherwise, bootstrap_abs_bias_rel_eff[r][c]
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
    return _compute_relative_efficiency(bootstrap_bias_c_A, bootstrap_bias_c_B)


def _generate_normal_abs_bias_mean_rel_eff_sample(
        *args, weights=None, discard_asymptotic=False,
        asymptotic_free_energy_A=None, asymptotic_free_energy_B=None,
        **kwargs
):
    """Generate a bootstrap sample for the absolute bias mean relative efficiency.

    Parameters
    ----------
    *args
        Positional arguments to forward to _generate_normal_abs_bias_rel_eff_sample.
    weights : numpy.ndarray or None, optional
        If given, weights[c] is the weight to be used in the weighted
        mean of the relative efficiency. If None, the average won't be
        weighted. Default is None.
    discard_asymptotic : bool, optional
        If True and at least one of the two biases is forced to 0 in the last
        data point, this is ignored in the average.
    asymptotic_free_energy_A : float, optional
        If given, this will be used as the asymptotic free energy of
        method A to compute the bias. Otherwise, this is estimated
        as mean_c_A[-1].
    asymptotic_free_energy_B : float, optional
        If given, this will be used as the asymptotic free energy of
        method B to compute the bias. Otherwise, this is estimated
        as mean_c_B[-1].
    **kwargs
        Keyword arguments to forward to _generate_normal_abs_bias_rel_eff_sample.

    Returns
    -------
    bootstrap_abs_bias_mean_rel_eff : float
        bootstrap_abs_bias_mean_rel_eff is the bootstrap sample of the
        absolute bias mean relative efficiency.

    See Also
    --------
    _generate_normal_abs_bias_rel_eff_sample

    """
    # Compute the mean relative efficiency.
    bootstrap_abs_bias_rel_eff = _generate_normal_abs_bias_rel_eff_sample(
        *args, asymptotic_free_energy_A=asymptotic_free_energy_A,
        asymptotic_free_energy_B=asymptotic_free_energy_B, **kwargs
    )
    if discard_asymptotic and (asymptotic_free_energy_A is None or asymptotic_free_energy_B is None):
        bootstrap_abs_bias_rel_eff = bootstrap_abs_bias_rel_eff[:-1]
        if weights is not None:
            weights = weights[:-1]
    return np.average(bootstrap_abs_bias_rel_eff, weights=weights, axis=0)


def _generate_normal_rmse_sample(
        mean_c, std_c, n_replicates,
        asymptotic_free_energy=None,
        n_bootstrap_samples=None
):
    """Generate a bootstrap sample for the RMSE.

    In this function, free_energy[c] at a fixed computational cost is
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
    """Generate a bootstrap sample for the RMSE relative efficiency.

    In this function, free_energy[c] at a fixed computational cost is
    assumed to be normally distributed. Under this assumption the sample
    mean is t-distributed and the standard deviation is chi-distributed
    with n_replicates-1 degrees of freedom, which allows us to generate
    a bootstrap sample for the absolute bias as a function of the computational
    cost from a parametric distribution.

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
    return _compute_relative_efficiency(bootstrap_rmse_c_A, bootstrap_rmse_c_B)


def _generate_normal_rmse_mean_rel_eff_sample(*args, weights=None, **kwargs):
    """Generate a bootstrap sample for the absolute bias mean relative efficiency.

    Parameters
    ----------
    *args
        Positional arguments to forward to _generate_normal_rmse_rel_eff_sample.
    weights : numpy.ndarray or None, optional
        If given, weights[c] is the weight to be used in the weighted
        mean of the relative efficiency. If None, the average won't be
        weighted. Default is None.
    **kwargs
        Keyword arguments to forward to _generate_normal_rmse_rel_eff_sample.

    Returns
    -------
    bootstrap_abs_bias_mean_rel_eff : float
        bootstrap_abs_bias_mean_rel_eff is the bootstrap sample of the
        absolute bias mean relative efficiency.

    See Also
    --------
    _generate_normal_abs_bias_rel_eff_sample

    """
    # Compute the mean relative efficiency.
    bootstrap_rmse_rel_eff = _generate_normal_rmse_rel_eff_sample(*args, **kwargs)
    return np.average(bootstrap_rmse_rel_eff, weights=weights, axis=0)


# =============================================================================
# Wrappers for arch.bootstrap.IIDBootstrap
# =============================================================================

def _cache_params_arch(
        free_energy_A, free_energy_B, params_cache, compute_mean=True,
):
    """Utility function to handle the cache of the free energy mean, std, and asymptotic DG."""
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
            params_cache['std_c_' + suffix] = _compute_normal_unbiased_std(free_energy)
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
        return _compute_relative_efficiency(params_cache['std_c_A'], params_cache['std_c_B'])
    return _generate_normal_std_rel_eff_sample(**params_cache)


def _generate_normal_std_mean_rel_eff_sample_arch(
        free_energy_A, free_energy_B, params=None, state=None,
        params_cache=None, weights=None
):
    """Wraps around _generate_normal_std_mean_rel_eff_sample for use with arch.bootstrap."""
    if params is None:
        std_rel_eff = _generate_normal_std_rel_eff_sample_arch(
            free_energy_A, free_energy_B,
            params_cache=params_cache
        )
        return np.average(std_rel_eff, weights=weights)
    return _generate_normal_std_mean_rel_eff_sample(weights=weights, **params_cache)


def _generate_normal_abs_bias_rel_eff_sample_arch(
        free_energy_A, free_energy_B, params=None, state=None,
        asymptotic_free_energy_A=None, asymptotic_free_energy_B=None,
        params_cache=None
):
    """Wraps around _generate_normal_abs_bias_rel_eff_sample for use with arch.bootstrap."""
    params_cache = _cache_params_arch(free_energy_A, free_energy_B,
                                      params_cache, compute_mean=True)

    if params is None:
        abs_bias_rel_eff = _compute_abs_bias_relative_efficiency(
            params_cache['mean_c_A'], params_cache['mean_c_B'],
            asymptotic_free_energy_A, asymptotic_free_energy_B
        )
    else:
        abs_bias_rel_eff = _generate_normal_abs_bias_rel_eff_sample(
            asymptotic_free_energy_A=asymptotic_free_energy_A,
            asymptotic_free_energy_B=asymptotic_free_energy_B,
            **params_cache
        )

    # We discard the last data point if we force it to have zero
    # bias to avoid numerical instabilities in the BCa method.
    if asymptotic_free_energy_A is None or asymptotic_free_energy_B is None:
        abs_bias_rel_eff = abs_bias_rel_eff[:-1]

    return abs_bias_rel_eff


def _generate_normal_abs_bias_mean_rel_eff_sample_arch(
        free_energy_A, free_energy_B, params=None, state=None, weights=None,
        params_cache=None, asymptotic_free_energy_A=None, asymptotic_free_energy_B=None
):
    """Wraps around _generate_normal_abs_bias_mean_rel_eff_sample for use with arch.bootstrap."""
    if params is None:
        abs_bias_rel_eff = _generate_normal_abs_bias_rel_eff_sample_arch(
            free_energy_A, free_energy_B, params_cache=params_cache,
            asymptotic_free_energy_A=asymptotic_free_energy_A,
            asymptotic_free_energy_B=asymptotic_free_energy_B
        )

        # _generate_normal_abs_bias_rel_eff_sample_arch already discard the last
        # data point so we need to discard only the last parameter and weight.
        if (asymptotic_free_energy_A is None or asymptotic_free_energy_B is None) and weights is not None:
            weights = weights[:-1]
        return np.average(abs_bias_rel_eff, weights=weights)

    return _generate_normal_abs_bias_mean_rel_eff_sample(
        weights=weights, discard_asymptotic=True,
        asymptotic_free_energy_A=asymptotic_free_energy_A,
        asymptotic_free_energy_B=asymptotic_free_energy_B,
        **params_cache
    )


def _generate_normal_rmse_rel_eff_sample_arch(
        free_energy_A, free_energy_B, params=None, state=None, n_replicates=None,
        asymptotic_free_energy_A=None, asymptotic_free_energy_B=None,
        params_cache=None
):
    """Wraps around _generate_normal_abs_bias_rel_eff_sample for use with arch.bootstrap."""
    params_cache = _cache_params_arch(free_energy_A, free_energy_B,
                                      params_cache, compute_mean=True)
    if params is None:
        return _compute_rmse_relative_efficiency(
            params_cache['mean_c_A'], params_cache['mean_c_B'],
            params_cache['std_c_A'], params_cache['std_c_B'],
            asymptotic_free_energy_A, asymptotic_free_energy_B,
        )

    return _generate_normal_rmse_rel_eff_sample(
        asymptotic_free_energy_A=asymptotic_free_energy_A,
        asymptotic_free_energy_B=asymptotic_free_energy_B,
        **params_cache
    )


def _generate_normal_rmse_mean_rel_eff_sample_arch(
        free_energy_A, free_energy_B, params=None, state=None, weights=None,
        params_cache=None, asymptotic_free_energy_A=None, asymptotic_free_energy_B=None
):
    """Wraps around _generate_normal_abs_bias_mean_rel_eff_sample for use with arch.bootstrap."""
    if params is None:
        rmse_rel_eff = _generate_normal_rmse_rel_eff_sample_arch(
            free_energy_A, free_energy_B, params_cache=params_cache,
            asymptotic_free_energy_A=asymptotic_free_energy_A,
            asymptotic_free_energy_B=asymptotic_free_energy_B
        )
        return np.average(rmse_rel_eff, weights=weights)

    return _generate_normal_rmse_mean_rel_eff_sample(
        weights=weights,
        asymptotic_free_energy_A=asymptotic_free_energy_A,
        asymptotic_free_energy_B=asymptotic_free_energy_B,
        **params_cache
    )
