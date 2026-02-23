import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import norm
import os
import matplotlib.colors as mcolors
from mle_fit import fn_mle_pc
from scipy.stats import norm, binom

periods = ['x2yr','x5yr','x10yr','x25yr','x50yr','x100yr']

num_periods = len(periods)

# Define a dark color palette
dark_colors = list(mcolors.TABLEAU_COLORS.values())[:num_periods]

# data_2 = pd.read_csv(os.path.join(os.getcwd(), 'studies/reliability/SB/final_precipitation_intensity_data.xlsx - Sheet1.csv'))

for i, period in enumerate(periods):

    if period == 'x100yr':
        plt.figure(figsize=(20, 12))

        # Load the data from the CSV file
        data_1 = pd.read_csv(os.path.join(os.getcwd(), f'studies/reliability/LAX sensitivity analysis/15min-return periods/Fragility_trials_{period}/fragility_data_points1976.csv'))

        # Extract the periods and intensities
        # intensities = data_2[period][~np.isnan(data_2[period])].sort_values(ascending=True).values

        # Assume the CSV file has columns named 'IM' and 'num_collapse' and 'num_gms'
        # IM = intensities
        IM = np.linspace(1,5,50)
        num_collapse = data_1[period][~np.isnan(data_1[period])].values * 200
        # num_gms = np.ones(len(data_1[period].values),1) * 200
        num_gms = np.full(len(data_1[period][~np.isnan(data_1[period])].values), 200)
        
        # Use the fn_mle_pc function to get the optimized parameters
        mu_opt, sigma_opt = fn_mle_pc(np.transpose(IM), np.transpose(num_gms), np.transpose(num_collapse))

        # Scatter plot of observed data points with smaller triangle markers
        plt.scatter(IM, num_collapse / num_gms, color=dark_colors[i], marker='^', s=50, label=f'Observed Data ({period})')

        # Generate IM values for plotting the fitted curve
        IM_plot = np.linspace(0, 7, 500)
        p_fit = norm.cdf(np.log(IM_plot), loc=mu_opt, scale=sigma_opt)

        # Plot the fitted fragility curve
        plt.plot(IM_plot, p_fit, color=dark_colors[i], linestyle='--', label=f'Fitted Lognormal CDF ({period})\nmu={mu_opt:.3f}, sigma={sigma_opt:.3f}')

        # Add labels and title
        plt.xlabel('IM')
        plt.ylabel('Probability of Limit State of Exceedance')
        plt.title(f'Fragility Function Fitting (15min duration)- LAX ({period} return period)')
        plt.legend(loc='best', fontsize=8)
        plt.grid(True)
        plt.savefig(f'studies/reliability/LAX sensitivity analysis/15min-return periods/Fragility_fit/fitted_plot_{period}.png')

        # Show the plot
        plt.show()  

plt.figure(figsize=(20, 12))

results_df = pd.DataFrame(columns=['period', 'mu_opt', 'sigma_opt'])

for i, period in enumerate(periods):
    # Load the data from the CSV file
    data_1 = pd.read_csv(os.path.join(os.getcwd(), f'studies/reliability/LAX sensitivity analysis/15min-return periods/Fragility_trials_{period}/fragility_data_points1976.csv'))
    
    # Assume the CSV file has columns named 'IM' and 'num_collapse' and 'num_gms'
    IM = np.linspace(1, 5, 50)
    num_collapse = data_1[period][~np.isnan(data_1[period])].values * 200
    num_gms = np.full(len(data_1[period][~np.isnan(data_1[period])].values), 200)

    # Use the fn_mle_pc function to get the optimized parameters
    mu_opt, sigma_opt = fn_mle_pc(np.transpose(IM), np.transpose(num_gms), np.transpose(num_collapse))
    results_df = results_df._append({'period': period, 'mu_opt': mu_opt, 'sigma_opt': sigma_opt}, ignore_index=True)

    # Plot the fitted fragility curve
    IM_plot = np.linspace(0, 7, 100)
    p_fit = norm.cdf(np.log(IM_plot), loc=mu_opt, scale=sigma_opt)
    plt.plot(IM_plot, p_fit, color=dark_colors[i], linestyle='--', label=f'Fitted Lognormal CDF ({period})\nmu={mu_opt:.3f}, sigma={sigma_opt:.3f}')

# Save the results DataFrame to a CSV file
results_df.to_csv(f'studies/reliability/LAX sensitivity analysis/15min-return periods/Fragility_fit/mu_sigma_opt_periodwise.csv', index=False)

# Add labels and title
plt.xlabel('IM')
plt.ylabel('Probability of Limit State of Exceedance') 
plt.title('Theoritical Fragility Curves - LAX (15min duration)')
plt.legend(loc='best', fontsize=8)
plt.grid(True)
plt.savefig(f'studies/reliability/LAX sensitivity analysis/15min-return periods/Fragility_fit/period_comparision_plot.png')

# Show the plot
plt.show()