import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import binom


def log_likelihood(p: float, k_reads, n_reads, epsilon: float) -> float:
    f_i = p * (1 - epsilon) + (1 - p) * epsilon
    f_i = np.clip(f_i, 1e-10, 1 - 1e-10)
    log_L_i = binom.logpmf(k_reads, n_reads, f_i)
    return float(np.sum(log_L_i))


def bayesian_estimation(
    k_reads,
    n_reads,
    epsilon: float,
    p_max: float = 0.05,
    p_grid_size: int = 4000,
    ci: float = 0.95,
):
    k_reads = np.asarray(k_reads)
    n_reads = np.asarray(n_reads)

    p_grid = np.linspace(0, p_max, p_grid_size)
    log_L_values = np.array([log_likelihood(p, k_reads, n_reads, epsilon) for p in p_grid])

    L_values = np.exp(log_L_values - np.max(log_L_values))
    posterior = L_values / np.sum(L_values)

    p_est_mean = float(np.sum(p_grid * posterior))
    cumulative_posterior = np.cumsum(posterior)

    interp_cdf = interp1d(
        cumulative_posterior,
        p_grid,
        bounds_error=False,
        fill_value=(p_grid[0], p_grid[-1]),
        assume_sorted=True,
    )

    p_est_median = float(interp_cdf(0.5))
    p_est_mode = float(p_grid[np.argmax(posterior)])

    alpha = (1 - ci) / 2
    lower_bound = float(interp_cdf(alpha))
    upper_bound = float(interp_cdf(1 - alpha))

    return {
        "p_grid": p_grid,
        "posterior": posterior,
        "p_est_mean": p_est_mean,
        "p_est_median": p_est_median,
        "p_est_mode": p_est_mode,
        "lower_bound": lower_bound,
        "upper_bound": upper_bound,
        "ci": ci,
        "epsilon": epsilon,
        "p_max": p_max,
    }
