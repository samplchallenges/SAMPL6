#!/usr/bin/env python

"""Utility functions for statistical analysis."""


# =============================================================================
# GLOBAL IMPORTS
# =============================================================================

import numpy as np
import scipy.stats
import scipy.special


# =============================================================================
# STATISTIC ESTIMATORS
# =============================================================================

def kendall_tau(data):
    x, y = data.T
    correlation, p_value = scipy.stats.kendalltau(x, y)
    return correlation


def r2(data):
    x, y = data.T
    slope, intercept, r_value, p_value, stderr = scipy.stats.linregress(x, y)
    return r_value**2


def slope(data):
    x, y = data.T
    slope, intercept, r_value, p_value, stderr = scipy.stats.linregress(x, y)
    return slope


def me(data):
    x, y = data.T
    error = np.array(x) - np.array(y)
    return error.mean()


def mae(data):
    x, y = data.T
    error = np.abs(np.array(x) - np.array(y))
    return error.mean()


def rmse(data):
    x, y = data.T
    error = np.array(x) - np.array(y)
    rmse = np.sqrt((error**2).mean())
    return rmse


# =============================================================================
# STANDARD ERROR OF THE MEAN, VARIANCE, AND STANDARD DEVIATION
# =============================================================================

def _unbiased_std_Kn(data):
    """Compute the factor used to unbias the Bessel-corrected sample std"""
    n = len(data)
    # Use the log gamma function for numerical stability.
    loggamma1 = scipy.special.gammaln((n-1)/2)
    loggamma2 = scipy.special.gammaln(n/2)
    k_n = np.sqrt((n-1)/2) * np.exp(loggamma1 - loggamma2)
    return k_n


def unbiased_std(data):
    """For n ~< 10, the sqrt of the Bessel-corrected variance is a biased estimate of the std.

    This function implement the unbiased estimate.

    See http://web.eecs.umich.edu/~fessler/papers/files/tr/stderr.pdf.

    """
    bessel_std = np.std(data, ddof=1)
    k_n = _unbiased_std_Kn(data)
    return k_n * bessel_std


def ci_unbiased_std(data, confidence=0.95):
    u_std = unbiased_std(data)
    k_n = _unbiased_std_Kn(data)
    # The error is chi-distributed
    n = len(data)
    scale = k_n * u_std / np.sqrt(n-1)
    chi_interval = scipy.stats.chi.interval(alpha=confidence, df=n-1, scale=scale)
    return chi_interval


def unbiased_sem(data):
    """Compute the standard error of the mean with the unbiased estimate of the std."""
    n = len(data)
    return unbiased_std(data) / np.sqrt(n)


def unbiased_mean_confidence_interval(data, confidence=0.95):
    norm_statistics = scipy.stats.norm.interval(alpha=confidence)[1]
    sem = unbiased_sem(data)
    mean = np.mean(data)
    return mean, norm_statistics * sem


def mean_confidence_interval(data, confidence=0.95):
    """Compute mean and t-based confidence interval for the data.

    Parameters
    ----------
    data : array like
    confidence : float

    Returns
    -------
    mean : numpy.ndarray
        The data mean.
    ci : float
        The confidence interval around the mean.
    """
    t_statistics = scipy.stats.t.interval(alpha=confidence, df=len(data)-1)[1]
    sem = scipy.stats.sem(data)
    mean = np.mean(data)
    return mean, t_statistics * sem


# =============================================================================
# BOOTSTRAP
# =============================================================================

def compute_bootstrap_statistics(samples, stats_funcs, percentile=0.95,
                                 n_bootstrap_samples=10000, sems=None):
    """Compute bootstrap confidence interval for the given statistics functions.

    Parameters
    ----------
    samples : np.ndarray
        The series of samples that will be bootstrapped. Each bootstrap
        sample is passed to the estimator as stat_func(samples[bootstrap_indices]).
    stats_funcs : list of callables
        The statistics estimators functions with signature stat_func(samples).
    percentile : float, optional
        The bootstrap percentile. Default is 0.95.
    n_bootstrap_samples : int, optional
        The number of bootstrap samples to sample. Default is 10000.
    sems : np.ndarray, optional
        The standard error of the means of the samples. If passed, each
        bootstrap cycle will re-sample each data point from a normal
        distribution with its standard deviation. This must have the same
        shape of samples. Data points for which SEMs is 0.0 won't be
        re-sampled.

    Returns
    -------
    bootstrap_statistics : list of tuples
        bootstrap_statistics[i] is (statistics, confidence_interval, bootstrap_samples)
        of stats_funcs[i]. confidence_interval is a pair (lower_bound, upper_bound),
        and bootstrap_samples are the (ordered) bootstrap statistics used to compute
        the confidence interval.

    """
    # Handle case where only a single function is passed.
    try:
        len(stats_funcs)
    except TypeError:
        stats_funcs = [stats_funcs]

    # Compute mean statistics.
    statistics = [stats_func(samples) for stats_func in stats_funcs]

    # Generate bootstrap statistics.
    bootstrap_samples_statistics = np.zeros((len(statistics), n_bootstrap_samples))
    for bootstrap_sample_idx in range(n_bootstrap_samples):
        # Re-sample from a normal distribution if the SEMs are given.
        if sems is None:
            bootstrap_samples = samples
        else:
            bootstrap_samples = resample_from_normal(samples, stds=sems)

        # Re-sample with replacement.
        samples_indices = np.random.randint(low=0, high=len(samples), size=len(samples))
        bootstrap_samples = bootstrap_samples[samples_indices]

        # Compute statistics for the bootstrap sample.
        for stats_func_idx, stats_func in enumerate(stats_funcs):
            bootstrap_samples_statistics[stats_func_idx][bootstrap_sample_idx] = stats_func(bootstrap_samples)

    # Compute confidence intervals.
    percentile_index = int(np.floor(n_bootstrap_samples * (1 - percentile) / 2)) - 1
    bootstrap_statistics = []
    for stats_func_idx, samples_statistics in enumerate(bootstrap_samples_statistics):
        samples_statistics.sort()
        stat_lower_percentile = samples_statistics[percentile_index]
        stat_higher_percentile = samples_statistics[-percentile_index+1]
        confidence_interval = (stat_lower_percentile, stat_higher_percentile)
        bootstrap_statistics.append([statistics[stats_func_idx], confidence_interval, samples_statistics])

    return bootstrap_statistics


def resample_from_normal(samples, stds):
    """Resample from a normal distribution.

    This handle SEMs == 0.0 that make numpy.random.normal throw an exception.

    Parameters
    ----------
    samples : 2D np.ndarray
        The series of samples that will be re-sampled..
    sems : 2D np.ndarray, optional
        The standard deviations of each sample. It must have the same shape
        of samples. Samples for which the standard deviations is 0.0 are not
        re-sampled.

    Returns
    -------
    new_samples : np.ndarray
        The samples resampled from a normal distribution with mean old_sample
        and standard deviation std.
    """
    new_samples = np.empty(samples.shape, dtype=samples.dtype)
    for sample_idx, (sample, std) in enumerate(zip(samples, stds)):
        # Handle case where individual samples are tuples and only a
        # subset of them has std 0.0
        for sample_subidx, (sample_i, std_i) in enumerate(zip(sample, std)):
            if std_i == 0.0:
                new_samples[sample_idx][sample_subidx] = sample_i
            else:
                new_samples[sample_idx][sample_subidx] = np.random.normal(loc=sample_i, scale=std_i)
    return new_samples
