#!/usr/bin/env python

"""Utility functions for statistical analysis."""


# =============================================================================
# GLOBAL IMPORTS
# =============================================================================

import numpy as np
import scipy.stats


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
# STANDARD ERROR OF THE MEAN
# =============================================================================

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
    t_statistics = scipy.stats.t.interval(alpha=0.95, df=len(data)-1)[1]
    sem = scipy.stats.sem(data)
    mean = np.mean(data)
    return mean, t_statistics * sem


# =============================================================================
# BOOTSTRAP
# =============================================================================

def compute_bootstrap_statistics(samples, stats_funcs, percentile=0.95, n_bootstrap_samples=10000):
    """Compute bootstrap confidence interval for the given statistics functions.

    Parameters
    ----------
    samples : array-like
        The series of samples that will be bootstrapped. Each bootstrap
        sample is passed to the estimator as stat_func(samples[bootstrap_indices]).
    stats_funcs : list of callables
        The statistics estimators functions with signature stat_func(samples).
    percentile : float, optional
        The bootstrap percentile.
    n_bootstrap_samples : int
        The number of bootstrap samples to sample.

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
        samples_indices = np.random.randint(low=0, high=len(samples), size=len(samples))
        for stats_func_idx, stats_func in enumerate(stats_funcs):
            bootstrap_samples_statistics[stats_func_idx][bootstrap_sample_idx] = stats_func(samples[samples_indices])

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
