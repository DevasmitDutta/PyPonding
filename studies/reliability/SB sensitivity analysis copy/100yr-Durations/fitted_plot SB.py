import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm, binom
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os

# --- Filtering function with tolerance ---
def filter_by_median_im(IM, num_gms, num_collapse, tol=1e-6):
    """
    Keep only one representative (median IM) per probability level,
    grouping probabilities within a tolerance.
    """
    prob = np.array(num_collapse) / np.array(num_gms)
    df = pd.DataFrame({
        'IM': IM,
        'num_gms': num_gms,
        'num_collapse': num_collapse,
        'prob': prob
    }).sort_values('prob').reset_index(drop=True)

    filtered_rows = []
    used = np.zeros(len(df), dtype=bool)

    for i in range(len(df)):
        if used[i]:
            continue
        p_ref = df.loc[i, 'prob']
        mask = np.isclose(df['prob'], p_ref, atol=tol)
        group = df[mask]
        used[mask] = True
        median_im = group['IM'].median()
        closest_row = group.iloc[(group['IM'] - median_im).abs().argmin()]
        filtered_rows.append(closest_row)

    filtered_df = pd.DataFrame(filtered_rows)
    return (filtered_df['IM'].values,
            filtered_df['num_gms'].values,
            filtered_df['num_collapse'].values)


# --- MLE fit ---
def mlefit(params, num_gms, num_collapse, IM):
    if params[1] < 0:
        return 1e10
    p = norm.cdf(np.log(IM), loc=params[0], scale=params[1])
    p = np.clip(p, 1e-10, 1 - 1e-10)
    likelihood = binom.pmf(num_collapse, num_gms, p)
    likelihood[likelihood == 0] = np.finfo(float).tiny
    return -np.sum(np.log(likelihood))


def fn_mle_pc(IM, num_gms, num_collapse, filter_data=True):
    if filter_data:
        IM, num_gms, num_collapse = filter_by_median_im(IM, num_gms, num_collapse)

    x0 = [np.mean(np.log(IM)), np.std(np.log(IM))]
    result = minimize(mlefit, x0, args=(num_gms, num_collapse, IM),
                      options={'maxiter': 1000}, method='Nelder-Mead')
    log_mean, log_stdev = result.x
    return log_mean, log_stdev, (IM, num_gms, num_collapse)


# --- Main execution ---
durations = ['x0.25hr','x1hr','x2hr','x3hr']
# shape_name = ['W16X26','W12X22','W8X24','W10X15']
dark_colors = list(mcolors.TABLEAU_COLORS.values())[:len(durations)]

for i, duration in enumerate(durations):
    # for j, shape in enumerate(shape_name):
        plt.figure(figsize=(20, 12))
        data_1 = pd.read_csv(
            os.path.join(os.getcwd(),
                         f'studies/reliability/SB sensitivity analysis copy/100yr-Durations/Fragility_trials_{duration}/fragility_data_points1976.csv')
        )

        IM = np.linspace(1, 16, 50)
        num_collapse = data_1[duration][~np.isnan(data_1[duration])].values * 200
        num_gms = np.full(len(num_collapse), 200)

        # Fit fragility with filtering
        mu_opt, sigma_opt, (IM_filt, num_gms_filt, num_collapse_filt) = fn_mle_pc(
            IM, num_gms, num_collapse, filter_data=True
        )

        # Scatter: filtered points only
        plt.scatter(IM_filt, num_collapse_filt / num_gms_filt,
                    color=dark_colors[i], marker='^', s=50,
                    label=f'Filtered Data ({duration})')

        # Plot fitted curve
        IM_plot = np.linspace(0.1, 20, 200)
        p_fit = norm.cdf(np.log(IM_plot), loc=mu_opt, scale=sigma_opt)
        plt.plot(IM_plot, p_fit, color=dark_colors[i], linestyle='--',
                 label=f'Fitted Lognormal CDF ({duration})\nlog_mean={mu_opt:.3f}, log_stdev={sigma_opt:.3f}')

        plt.xlabel('IM')
        plt.ylabel('Probability of Limit State of Exceedance')
        # Ensure output directory exists
        out_dir = f'studies/reliability/SB sensitivity analysis copy/100yr-Durations/Fragility_fit_median_{duration}'
        os.makedirs(out_dir, exist_ok=True)

        plt.title(f'Fragility Function Fitting (100yr return period) - Santa Barbara ({duration} duration)')
        plt.legend(loc='best', fontsize=8)
        plt.grid(True)
        plt.savefig(f'studies/reliability/SB sensitivity analysis copy/100yr-Durations/Fragility_fit_median_{duration}/fitted_plot.png')
        plt.show()

        print(f"Duration: {duration}")
        print(f"Log-Mean (mu): {mu_opt:.4f}")
        print(f"Log-StdDev (sigma): {sigma_opt:.4f}")
