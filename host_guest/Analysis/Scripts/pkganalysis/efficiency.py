#!/usr/bin/env python

"""Functions for the efficiency statistical analysis."""


# =============================================================================
# GLOBAL IMPORTS
# =============================================================================

import copy

import numpy as np
import scipy as sp


# =============================================================================
# UTILITIES
# =============================================================================

def discard_initial_zeros(mean_c_A, mean_c_B, std_c_A, std_c_B):
    """Discard the initial data points for which there's no data.

    Some submissions don't have estimates for the first part of
    the trajectory. This function removes the initial data points
    with mean free energy 0.0.

    Parameters
    ----------
    mean_c_A : numpy.ndarray
        mean_c_A[c] is the mean free energy for method A at the c-th
        computational cost.
    std_c_A : numpy.ndarray
        std_c_A[c] is the standard deviation of the free energy estimate
        for method A at the c-th computational cost.
    mean_c_B : numpy.ndarray
        mean_c_B[c] is the mean free energy for method B at the c-th
        computational cost.
    std_c_B : numpy.ndarray
        std_c_B[c] is the standard deviation of the free energy estimate
        for method A at the c-th computational cost.

    Returns
    -------
    mean_c_A : numpy.ndarray
    std_c_A : numpy.ndarray
    mean_c_B : numpy.ndarray
    std_c_B : numpy.ndarray
        The new data after removing the initial zeros.
    """
    first_nonzero_idx_A = np.nonzero(mean_c_A)[0][0]
    first_nonzero_idx_B = np.nonzero(mean_c_B)[0][0]
    first_nonzero_idx = max(first_nonzero_idx_A, first_nonzero_idx_B)

    mean_c_A = copy.deepcopy(mean_c_A[first_nonzero_idx:])
    mean_c_B = copy.deepcopy(mean_c_B[first_nonzero_idx:])
    std_c_A = copy.deepcopy(std_c_A[first_nonzero_idx:])
    std_c_B = copy.deepcopy(std_c_B[first_nonzero_idx:])
    return mean_c_A, std_c_A, mean_c_B, std_c_B


# =============================================================================
# MAIN FUNCTIONS
# =============================================================================

def compute_mean_relative_efficiencies(
        mean_c_A, std_c_A, asymptotic_free_energy_A,
        mean_c_B, std_c_B, asymptotic_free_energy_B,
        n_replicates, n_bootstrap_samples=10000,
        confidence_interval=0.95, weighted=True,
        discard_zero_bias=True
):
    """Compute the arithmetic mean relative std, absolute bias, and RMSD efficiency for the data.

    Parameters
    ----------
    mean_c_A : np.ndarray
        mean_c_A[c] is the mean free energy of method A at the c-th
        computational cost.
    std_c_A : np.ndarray
        std_c_A[c] is the standard deviation of the free energy estimate
        of method A at the c-th computational cost.
    asymptotic_free_energy_A : float
        The asymptotic free energy of method A used to compute the bias.
    mean_c_B : np.ndarray
        mean_c_B[c] is the mean free energy of method B at the c-th
        computational cost.
    std_c_B : np.ndarray
        std_c_B[c] is the standard deviation of the free energy estimate
        of method B at the c-th computational cost.
    asymptotic_free_energy_B : float
        The asymptotic free energy of method B used to compute the bias.
    n_replicates : int
        The number of independent replicate calculations used to estimate
        the standard deviation.
    n_bootstrap_samples : int, optional
        The number of bootstrap samples per computational cost used to
        compute the weights and the mean relative efficiencies confidence
        intervals. Default is 10000.
    confidence_interval : float or None
        The inter-percentile range used for the confidence interval. If
        None, no confidence interval is returned.
    weighted : bool, optional
        If True, the mean is weighted by the inverse variance of the
        efficiency.
    discard_zero_bias : bool, optional
        If True, data points that have bias that is exactly 0.0 are discarded
        for the calculation of the absolute bias relative efficiency to avoid
        having an undefined relative efficiency.

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
        defining the percentile confidence interval of the std mean relative
        efficiency computed with bootstrapping.
    abs_bias_mean_efficiency_CI : List[float], optional
        If confidence_interval is not None, this is the pair (low_bound, high_bound)
        defining the percentile confidence interval of the absolute bias mean relative
        efficiency computed with bootstrapping.
    rmse_mean_efficiency_CI : List[float], optional
        If confidence_interval is not None, this is the pair (low_bound, high_bound)
        defining the percentile confidence interval of the RMSE mean relative
        efficiency computed with bootstrapping.

    """
    # Check that all means and stds have the same number of data points.
    lengths = {len(x) for x in [mean_c_A, mean_c_B, std_c_A, std_c_B]}
    if len(lengths) != 1:
        raise ValueError('All means and stds must have the same length')

    # Compute bias.
    bias_c_A = mean_c_A - asymptotic_free_energy_A
    bias_c_B = mean_c_B - asymptotic_free_energy_B

    # Compute the relative efficiencies for each computational cost.
    # relative_efficiencies is a triple with the relative efficiencies
    # of std, absolute bias, and RMSE.
    relative_efficiencies = compute_relative_efficiencies(
        std_c_A, bias_c_A, std_c_B, bias_c_B)
    n_statistics = len(relative_efficiencies)

    # Generate relative efficiency bootstrap distributions if we need
    # to compute the weighted mean and or confidence intervals.
    if weighted or (confidence_interval is not None):
        # efficiency_distributions is a triple with the bootstrap distributions
        # of std, absolute bias, and RMSE relative efficiencies.
        efficiency_distributions = generate_relative_efficiency_bootstrap_distribution(
            mean_c_A, std_c_A, asymptotic_free_energy_A,
            mean_c_B, std_c_B, asymptotic_free_energy_B,
            n_replicates, n_bootstrap_samples
        )

    # Determine weights if requested.
    if weighted:
        weights = np.empty(shape=(n_statistics, len(mean_c_A)))
        for i, efficiency_distribution in enumerate(efficiency_distributions):
            weights[i] = compute_efficiencies_mean_weights(efficiency_distribution)
    else:
        weights = [None for _ in range(n_statistics)]

    # Take the mean of the relative efficiencies.
    mean_efficiencies = [None for _ in range(n_statistics)]
    for i, rel_eff_c in enumerate(relative_efficiencies):
        mean_efficiencies[i] = _compute_weighted_mean(rel_eff_c, weights[i], remove_nans=True)

    # Determine confidence intervals.
    if confidence_interval is not None:
        confidence_intervals = [[None, None] for _ in range(n_statistics)]
        for i, efficiency_distribution in enumerate(efficiency_distributions):
            confidence_intervals[i] = compute_mean_efficiency_confidence_interval(
                efficiency_distribution, weights_c=weights[i], percentile=confidence_interval)
    else:
        # This list is appended to the return value.
        confidence_intervals = []

    return mean_efficiencies + confidence_intervals


def _compute_weighted_mean(x, weights, remove_nans=True):
    """Compute weighted average of the data discarding NaN data points."""
    not_nan_indices = ~np.isnan(x)
    x = x[not_nan_indices]
    if weights is not None:
        weights = weights[not_nan_indices]
    return float(np.average(x, weights=weights))


# TODO: REMOVE THIS WHEN YOU FINALIZE DECISION ABOUT RELATIVE EFFICIENCY DEFINITION
def compute_relative_efficiencies(*args):
    return compute_relative_efficiencies_diff(*args)


def compute_relative_efficiencies_ratio(
        std_A, bias_A, std_B, bias_B
):
    """Compute std, absolute bias, and RMSE relative effiencies for all computational costs.

    Absolute bias relative efficiencies for computational costs that are
    exactly zeros are set to np.nan.

    Parameters
    ----------
    std_A : np.ndarray
        This can be a 1D or 2D array of standard deviations for method A.
        If 1D, std_c[c] could be the standard deviation of the free energy
        estimate at the c-th computational cost. If 2D, std_c[i][c] is
        the i-th bootstrap sample for the c-th computational cost.
    bias_A : np.ndarray
        1D or 2D array of biases.
    std_B : np.ndarray
        1D or 2D array of standard deviation for the reference calculation.
    bias_B : np.ndarray
        1D or 2D array of biases for the reference calculation.

    Returns
    -------
    std_efficiencies : np.ndarray
        The array of std relative efficiencies. The shape is the same as std.
    abs_bias_efficiencies : np.ndarray
        The array of absolute bias relative efficiencies. The shape is the
        same as std.
    rmse_efficiencies : np.ndarray
        The array of RMSE relative efficiencies. The shape is the same as std.

    """
    rmse_c_A = np.sqrt(std_A**2 + bias_A**2)
    rmse_c_B = np.sqrt(std_B**2 + bias_B**2)

    std_efficiencies =  std_A / std_B
    rmse_efficiencies = rmse_c_A / rmse_c_B
    abs_bias_efficiencies = np.abs(bias_A) / np.abs(bias_B)

    # Substitute 0.0 or np.inf absolute bias with np.nan.
    inf_indices = np.isinf(abs_bias_efficiencies)
    zero_indices = np.isclose(abs_bias_efficiencies, 0.0)
    for indices in [inf_indices, zero_indices]:
        abs_bias_efficiencies[indices] = np.nan

    return std_efficiencies, abs_bias_efficiencies, rmse_efficiencies


def compute_relative_efficiencies_diff(
        std_A, bias_A, std_B, bias_B
):
    """Compute std, absolute bias, and RMSE relative effiencies for all computational costs.

    Parameters
    ----------
    std_A : np.ndarray
        This can be a 1D or 2D array of standard deviations for method A.
        If 1D, std_c[c] could be the standard deviation of the free energy
        estimate at the c-th computational cost. If 2D, std_c[i][c] is
        the i-th bootstrap sample for the c-th computational cost.
    bias_A : np.ndarray
        1D or 2D array of biases.
    std_B : np.ndarray
        1D or 2D array of standard deviation for the reference calculation.
    bias_B : np.ndarray
        1D or 2D array of biases for the reference calculation.

    Returns
    -------
    std_efficiencies : np.ndarray
        The array of std relative efficiencies. The shape is the same as std.
    abs_bias_efficiencies : np.ndarray
        The array of absolute bias relative efficiencies. The shape is the
        same as std.
    rmse_efficiencies : np.ndarray
        The array of RMSE relative efficiencies. The shape is the same as std.

    """
    rmse_c_A = np.sqrt(std_A**2 + bias_A**2)
    rmse_c_B = np.sqrt(std_B**2 + bias_B**2)

    std_efficiencies =  std_A - std_B
    rmse_efficiencies = rmse_c_A - rmse_c_B
    abs_bias_efficiencies = np.abs(bias_A) - np.abs(bias_B)

    return std_efficiencies, abs_bias_efficiencies, rmse_efficiencies


def compute_relative_efficiencies_old(
        std_A, bias_A, std_B, bias_B,
        discard_zero_bias=True
):
    """Compute std, absolute bias, and RMSE relative effiencies for all computational costs.

    Parameters
    ----------
    std_A : np.ndarray
        This can be a 1D or 2D array of standard deviations for method A.
        If 1D, std_c[c] could be the standard deviation of the free energy
        estimate at the c-th computational cost. If 2D, std_c[i][c] is
        the i-th bootstrap sample for the c-th computational cost.
    bias_A : np.ndarray
        1D or 2D array of biases.
    std_B : np.ndarray
        1D or 2D array of standard deviation for the reference calculation.
    bias_B : np.ndarray
        1D or 2D array of biases for the reference calculation.
    discard_zero_bias : bool, optional
        If True, zeros in the bias are discarded when computing the
        absolute bias efficiency to avoid having an undefined relative
        efficiency.

    Returns
    -------
    std_efficiencies : np.ndarray
        The array of std relative efficiencies. The shape is the same as std.
    abs_bias_efficiencies : np.ndarray
        The array of absolute bias relative efficiencies. The shape is the
        same as std.
    rmse_efficiencies : np.ndarray
        The array of RMSE relative efficiencies. The shape is the same as std.

    """
    rmse_c_A = np.sqrt(std_A**2 + bias_A**2)
    rmse_c_B = np.sqrt(std_B**2 + bias_B**2)

    std_efficiencies =  std_A / std_B
    rmse_efficiencies = rmse_c_A / rmse_c_B

    def _discard_zero_bias(a, b):
        for searched in [a, b]:
            nonzero_indices = ~np.isclose(searched, 0.0)
            a = a[nonzero_indices]
            b = b[nonzero_indices]
        return a, b

    # For the bias, discard all elements that have exactly 0.0 bias,
    # which can cause the the relative efficiency to be infinity or 0.
    # In RMSE, the contribution from std prevents it to be undefined.
    if not discard_zero_bias:
        abs_bias_efficiencies = np.abs(bias_A) / np.abs(bias_B)
    elif len(std_efficiencies.shape) == 1:
        # This is a 1D array.
        bias_A, bias_B = _discard_zero_bias(bias_A, bias_B)
        abs_bias_efficiencies = np.abs(bias_A) / np.abs(bias_B)
    else:
        # This is a 2D array.
        abs_bias_efficiencies = [None for _ in std_efficiencies]
        for i in range(len(abs_bias_efficiencies)):
            bias_A_i, bias_B_i = _discard_zero_bias(bias_A[i], bias_B[i])
            abs_bias_efficiencies[i] = np.abs(bias_A_i) / np.abs(bias_B_i)
        abs_bias_efficiencies = np

    return std_efficiencies, abs_bias_efficiencies, rmse_efficiencies


def _unbiased_std_Kn(n_replicates):
    """Compute the factor used to unbias the Bessel-corrected sample std"""
    # Use the log gamma function for numerical stability.
    loggamma1 = sp.special.gammaln((n_replicates-1)/2)
    loggamma2 = sp.special.gammaln(n_replicates/2)
    k_n = np.sqrt((n_replicates-1)/2) * np.exp(loggamma1 - loggamma2)
    return k_n


def generate_std_bootstrap_distribution(
        std_c, n_replicates, n_bootstrap_samples=10000
):
    """Generate a bootstrap distribution for the standard deviations.

    The underlying data is assumed to be normally distributed so that
    the std is supposed to be chi-distributed with n_replicates-1 degrees
    of freedom.

    Parameters
    ----------
    std_c : array-like
        std_c[c] is the standard deviation of the free energy estimate
        at the c-th computational cost.
    n_replicates : int
        The number of independent replicate calculations used to estimate
        the standard deviation.
    n_bootstrap_samples : int, optional
        The number of bootstrap samples per computational cost. Default
        is 10000.

    Returns
    -------
    std_distribution : array-like
        std_distribution[i][c] is the i-th bootstrap standard deviation
        at the c-th computational cost.

    """
    n_stds = len(std_c)
    df = n_replicates - 1
    size = (n_bootstrap_samples, n_stds)
    kn = _unbiased_std_Kn(n_replicates)
    scale = kn * std_c / np.sqrt(df)
    # return sp.stats.chi.rvs(df=df, scale=std_c/np.sqrt(df), size=size)
    return sp.stats.chi.rvs(df=df, scale=scale, size=size)


def generate_bias_bootstrap_distribution(
        mean_c, std_c, asymptotic_free_energy,
        n_replicates, n_bootstrap_samples=10000
):
    """Generate a bootstrap distribution for the means.

    The underlying data is assumed to be normally distributed so that
    the mean is supposed to be t-distributed with N-1 degrees of freedom.

    Parameters
    ----------
    mean_c : float or array-like
        mean_c[c] is the mean of the free energy at the c-th computational
        cost.
    std_c : array-like
        std_c[c] is the standard deviation of the free energy estimate
        at the c-th computational cost.
    asymptotic_free_energy : float
        The asymptotic free energy used to compute the bias.
    n_replicates : int
        The number of independent replicate calculations used to estimate
        the standard deviation.
    n_bootstrap_samples : int, optional
        The number of bootstrap samples per computational cost. Default
        is 10000.

    Returns
    -------
    bias_distribution : array-like
        bias_distribution[i][c] is the i-th bootstrap bias at the c-th
        computational cost.
    """
    n_means = len(mean_c)
    df = n_replicates - 1
    size = (n_bootstrap_samples, n_means)
    bias_c = mean_c - asymptotic_free_energy
    return sp.stats.t.rvs(df=df, scale=std_c/np.sqrt(n_replicates), size=size) + bias_c


def generate_relative_efficiency_bootstrap_distribution(
        mean_c_A, std_c_A, asymptotic_free_energy_A,
        mean_c_B, std_c_B, asymptotic_free_energy_B,
        n_replicates, n_bootstrap_samples=10000
):
    """Generate bootstrap distribution for std, absolute bias, and RMSE efficiencies.

    Parameters
    ----------
    mean_c_A : np.ndarray
        mean_c_A[c] is the mean free energy of method A at the c-th
        computational cost.
    std_c_A : np.ndarray
        std_c_A[c] is the standard deviation of the free energy estimate
        of method A at the c-th computational cost.
    asymptotic_free_energy_A : float
        The asymptotic free energy of method A used to compute the bias.
    mean_c_B : np.ndarray
        mean_c_B[c] is the mean free energy of method B at the c-th
        computational cost.
    std_c_B : np.ndarray
        std_c_B[c] is the standard deviation of the free energy estimate
        of method B at the c-th computational cost.
    asymptotic_free_energy_B : float
        The asymptotic free energy of method B used to compute the bias.
    n_replicates : int
        The number of independent replicate calculations used to estimate
        the standard deviation.
    n_bootstrap_samples : int, optional
        The number of bootstrap samples per computational cost. Default
        is 10000.

    Returns
    -------
    std_efficiency_distribution : np.ndarray
        std_efficiency_distribution[i][c] is the i-th bootstrap sample
        std relative efficiency generated at the c-th computational cost.
    abs_bias_efficiency_distribution : np.ndarray
        abs_bias_efficiency_distribution[i][c] is the i-th bootstrap sample
        absolute bias relative efficiency generated at the c-th computational
        cost.
    rmse_efficiency_distribution : np.ndarray
        rmse_efficiency_distribution[i][c] is the i-th bootstrap sample
        RMSE relative efficiency generated at the c-th computational cost.
    """
    kwargs = {'n_replicates': n_replicates, 'n_bootstrap_samples': n_bootstrap_samples}

    # Create std and bias bootstrap distribution.
    std_distribution_A = generate_std_bootstrap_distribution(std_c_A, **kwargs)
    std_distribution_B = generate_std_bootstrap_distribution(std_c_B, **kwargs)
    bias_distribution_A = generate_bias_bootstrap_distribution(
        mean_c_A, std_c_A, asymptotic_free_energy_A, **kwargs)
    bias_distribution_B = generate_bias_bootstrap_distribution(
        mean_c_B, std_c_B, asymptotic_free_energy_B, **kwargs)

    # Propagate distributions.
    efficiency_distributions = compute_relative_efficiencies(
        std_distribution_A, bias_distribution_A, std_distribution_B, bias_distribution_B)
    return efficiency_distributions


def compute_efficiencies_mean_weights(efficiency_distribution):
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


def compute_relative_efficiencies_confidence_interval(efficiency_distribution, percentile=0.95):
    """Compute the mean weights from the bootstrap distribution of the efficiency.

    Parameters
    ----------
    efficiency_distribution : numpy.ndarray
        efficiency_distribution[i][c] is the i-th bootstrap sample
        relative efficiency generated at the c-th computational cost.
    percentile : float
        The significance of the percentile confidence interval.

    Returns
    -------
    low_bound_c : numpy.ndarray
        low_bound_c[c] is the low bound of the confidence interval for
        the c-th computational cost.
    high_bound_c : numpy.ndarray
        high_bound_c[c] is the high bound of the confidence interval for
        the c-th computational cost.

    """
    low_percentile = 100 * (1 - percentile) / 2
    high_percentile = 100* percentile + low_percentile

    n_costs = efficiency_distribution.shape[1]
    low_bound_c = np.empty(shape=n_costs)
    high_bound_c = np.empty(shape=n_costs)
    for c in range(n_costs):
        low_bound_c[c], high_bound_c[c] = np.percentile(efficiency_distribution[:,c],
                                                        q=[low_percentile, high_percentile])
    return low_bound_c, high_bound_c


def compute_mean_efficiency_confidence_interval(
        efficiency_distribution, weights_c=None, percentile=0.95
):
    """Compute a percentile confidence interval for the mean relative efficiency.

    Parameters
    ----------
    efficiency_distribution : numpy.ndarray
        efficiency_distribution[i][c] is the i-th bootstrap sample
        relative efficiency generated at the c-th computational cost.
    weights_c : numpy.ndarray
        weights_c[c] is the weight of the relative efficiency at the
        c-th computational cost.
    percentile : float
        The significance of the percentile confidence interval.

    Returns
    -------
    low_bound : float
        The low bound of the percentile confidence interval.
    high_bound : float
        The high bound of the percentile confidence interval.

    """
    # Use the relative efficiency bootstrap distribution to generate
    # a bootstrap distribution of the mean relative efficiency.
    mean_efficiency_distribution = np.empty(len(efficiency_distribution))
    for i, rel_eff_c in enumerate(efficiency_distribution):
        mean_efficiency_distribution[i] = np.average(rel_eff_c, weights=weights_c)

    # Compute percentile confidence intervals.
    low_percentile = (1 - percentile) / 2 * 100
    high_percentile = 100* percentile + low_percentile
    return np.percentile(mean_efficiency_distribution, q=[low_percentile, high_percentile]).tolist()

