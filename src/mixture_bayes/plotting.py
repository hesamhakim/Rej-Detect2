import matplotlib.pyplot as plt
import numpy as np


def plot_posterior_distributions(results_list, title: str = "Posterior Distributions"):
    plt.figure(figsize=(12, 8))

    for results in results_list:
        plt.plot(
            results["p_grid"] * 100,
            results["posterior"],
            label=f"True p = {results.get('p_true', float('nan')) * 100:.2f}%",
            alpha=0.7,
        )
        if "p_true" in results:
            plt.axvline(results["p_true"] * 100, color="black", linestyle="--", alpha=0.35)

    plt.xlabel("Mixture Proportion p (%)")
    plt.ylabel("Posterior Probability")
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_true_vs_estimated(results_list, title: str = "True vs Estimated Mixture Proportions"):
    true_p_values = np.array([res["p_true"] for res in results_list])
    estimated_p_means = np.array([res["p_est_mean"] for res in results_list])
    estimated_p_medians = np.array([res["p_est_median"] for res in results_list])
    estimated_p_modes = np.array([res["p_est_mode"] for res in results_list])

    plt.figure(figsize=(10, 6))
    plt.plot(true_p_values * 100, estimated_p_means * 100, "o-", label="Estimated p (Mean)")
    plt.plot(true_p_values * 100, estimated_p_medians * 100, "s--", label="Estimated p (Median)")
    plt.plot(true_p_values * 100, estimated_p_modes * 100, "d:", label="Estimated p (Mode)")

    x_min = 0
    x_max = max(float(np.max(true_p_values) * 100), 2.0)
    plt.plot([x_min, x_max], [x_min, x_max], "k--", label="Ideal")

    plt.xlabel("True Mixture Proportion p (%)")
    plt.ylabel("Estimated Mixture Proportion p (%)")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.show()
