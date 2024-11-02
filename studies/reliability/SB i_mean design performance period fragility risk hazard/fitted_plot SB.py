import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import norm
import os
import matplotlib.colors as mcolors
from mle_fit import fn_mle_pc
from scipy.stats import norm, binom

def fit_fragility_curve():    
        # Define a dark color palette
        results_df = pd.DataFrame(columns=['mu_opt', 'sigma_opt'])

        plt.figure(figsize=(20, 12))

        # Load the data from the CSV file
        data_1 = pd.read_csv(os.path.join(os.getcwd(), f'studies/reliability/SB i_mean design performance period fragility risk hazard/Fragility_trials/fragility_data_points1976.csv'))
        # Assume the CSV file has columns named 'IM' and 'num_collapse' and 'num_gms'
        IM = np.linspace(1,4,100)
        num_collapse = data_1['xnohr'].values * 200
        num_gms = np.full(len(data_1['xnohr'].values), 200)
        
        # Use the fn_mle_pc function to get the optimized parameters
        mu_opt, sigma_opt = fn_mle_pc(np.transpose(IM), np.transpose(num_gms), np.transpose(num_collapse))
        results_df = results_df._append({'mu_opt': mu_opt, 'sigma_opt': sigma_opt}, ignore_index=True)

        # Scatter plot of observed data points with smaller triangle markers
        plt.scatter(IM, num_collapse / num_gms, marker='^', s=50, label=f'Observed Data')
        print(num_collapse / num_gms)

        # Generate IM values for plotting the fitted curve
        IM_plot = np.linspace(0, 7, 500)
        p_fit = norm.cdf(np.log(IM_plot), loc=mu_opt, scale=sigma_opt)

        # Plot the fitted fragility curve
        plt.plot(IM_plot, p_fit, linestyle='--', label=f'Fitted Lognormal CDF ()\nmu={mu_opt:.3f}, sigma={sigma_opt:.3f}')
        results_df.to_csv(f'studies/reliability/SB i_mean design performance period fragility risk hazard/Fragility_fit/mu_sigma_opt_durationwise.csv', index=False)

        # Add labels and title
        plt.xlabel('IM')
        plt.ylabel('Probability of Limit State of Exceedance')
        plt.title(f'Fragility Function Fitting - Santa Barbara')
        plt.legend(loc='best', fontsize=8)
        plt.grid(True)
        plt.savefig(f'studies/reliability/SB i_mean design performance period fragility risk hazard/Fragility_fit/fitted_plot.png')

        # Show the plot
        plt.show()  

fit_fragility_curve()
