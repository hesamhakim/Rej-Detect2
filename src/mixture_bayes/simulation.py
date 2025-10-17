import numpy as np

from .inference import bayesian_estimation


def generate_informative_positions(genome_length: int, n_positions: int, rng: np.random.Generator | None = None):
    if rng is None:
        rng = np.random.default_rng()
    return rng.choice(np.arange(1, genome_length + 1), n_positions, replace=False)


def simulate_reads(
    p_true: float,
    n_positions: int,
    mean_depth: float,
    epsilon: float,
    positions=None,
    rng: np.random.Generator | None = None,
):
    if rng is None:
        rng = np.random.default_rng()

    if positions is None:
        positions = np.arange(1, n_positions + 1)

    n_reads = rng.poisson(lam=mean_depth, size=n_positions)

    f_i = p_true * (1 - epsilon) + (1 - p_true) * epsilon
    k_reads = rng.binomial(n_reads, f_i)

    return k_reads, n_reads, np.asarray(positions)


def simulate_and_estimate(
    p_true: float,
    n_positions: int,
    mean_depth: float,
    epsilon: float,
    positions=None,
    p_max: float = 0.05,
    p_grid_size: int = 4000,
    rng: np.random.Generator | None = None,
):
    k_reads, n_reads, positions = simulate_reads(
        p_true=p_true,
        n_positions=n_positions,
        mean_depth=mean_depth,
        epsilon=epsilon,
        positions=positions,
        rng=rng,
    )

    results = bayesian_estimation(
        k_reads=k_reads,
        n_reads=n_reads,
        epsilon=epsilon,
        p_max=p_max,
        p_grid_size=p_grid_size,
    )
    results["p_true"] = p_true
    results["positions"] = positions
    results["k_reads"] = k_reads
    results["n_reads"] = n_reads
    return results
