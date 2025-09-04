import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import norm
import os
import matplotlib.colors as mcolors
from mle_fit import fn_mle_pc
from scipy.stats import norm, binom

durations = ['x100yr']

shape_name = ['W16X26','W12X22','W8X24','W10X15']


num_durations = len(durations)

# Define a dark color palette
dark_colors = list(mcolors.TABLEAU_COLORS.values())[:num_durations]

# data_2 = pd.read_csv(os.path.join(os.getcwd(), 'studies/reliability/SB/final_precipitation_intensity_data.xlsx - Sheet1.csv'))

for i, duration in enumerate(durations):
    for j, shape in enumerate(shape_name):
        # Create a new figure for each duration and shape
        plt.figure(figsize=(20, 12))

        # Load the data from the CSV file
        data_1 = pd.read_csv(os.path.join(os.getcwd(), f'studies/reliability/LAX_Design_analysis/15 min/Fragility_trials_{duration}_{shape_name[j]}/fragility_data_points1976.csv'))

        # Extract the durations and intensities
        IM = np.linspace(1,16,50)
        num_collapse = data_1[duration][~np.isnan(data_1[duration])].values * 200
        num_gms = np.full(len(data_1[duration][~np.isnan(data_1[duration])].values), 200)

        # Use the fn_mle_pc function to get the optimized parameters
        mu_opt, sigma_opt = fn_mle_pc(np.transpose(IM), np.transpose(num_gms), np.transpose(num_collapse))

        # Scatter plot of observed data points with smaller triangle markers
        plt.scatter(IM, num_collapse / num_gms, color=dark_colors[i], marker='^', s=50, label=f'Observed Data ({duration})')

        # Generate IM values for plotting the fitted curve
        IM_plot = np.linspace(0, 20, 100)
        p_fit = norm.cdf(np.log(IM_plot), loc=mu_opt, scale=sigma_opt)

        # Plot the fitted fragility curve
        plt.plot(IM_plot, p_fit, color=dark_colors[i], linestyle='--', label=f'Fitted Lognormal CDF ({duration})\nmu={mu_opt:.3f}, sigma={sigma_opt:.3f}')

        # Add labels and title
        plt.xlabel('IM')
        plt.ylabel('Probability of Limit State of Exceedance')
        # Ensure output directory exists
        out_dir = f'studies/reliability/LAX_Design_analysis/15 min/Fragility_fit_1000samp_{duration}_{shape_name[j]}'
        os.makedirs(out_dir, exist_ok=True)
        plt.title(f'Fragility Function Fitting (100yr return period)- LAX ({duration} duration)')
        plt.legend(loc='best', fontsize=8)
        plt.grid(True)
        plt.savefig(f'studies/reliability/LAX_Design_analysis/15 min/Fragility_fit_1000samp_{duration}_{shape_name[j]}/fitted_plot.png')
        plt.show()

