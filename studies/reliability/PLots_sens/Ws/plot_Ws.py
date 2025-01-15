import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import norm
import os
import matplotlib.colors as mcolors
from mle_fit import fn_mle_pc 
from scipy.stats import norm, binom

durations = ['x0.25hr']
num_durations = len(durations)
ws_vector = [6, 8, 12]

# Define a dark color palette
dark_colors = list(mcolors.TABLEAU_COLORS.values())[:len(ws_vector)] 
print(dark_colors)    

plt.figure(figsize=(20, 12))

for duration in durations:
    # Load the data from the CSV file
    for i, ws in enumerate(ws_vector):
        data_1 = pd.read_csv(os.path.join(os.getcwd(), f'studies/reliability/LAX analysis sens Ws_{ws}/100yr-Durations/Fragility_trials_{duration}/fragility_data_points1976.csv'))
        # Assume the CSV file has columns named 'IM' and 'num_collapse' and 'num_gms'
        IM = np.linspace(1,16,50)
        num_collapse = data_1[duration][~np.isnan(data_1[duration])].values * 200
        num_gms = np.full(len(data_1[duration][~np.isnan(data_1[duration])].values), 200)

        # Use the fn_mle_pc function to get the optimized parameters
        mu_opt, sigma_opt = fn_mle_pc(np.transpose(IM), np.transpose(num_gms), np.transpose(num_collapse))
        # Scatter plot of observed data points with smaller triangle markers
        # plt.scatter(IM, num_collapse / num_gms, color=dark_colors[i], marker='^', s=50, label=f'Observed Data ({duration})')

        # Plot the fitted fragility curve
        IM_plot = np.linspace(0, 20, 100)
        p_fit = norm.cdf(np.log(IM_plot), loc=mu_opt, scale=sigma_opt)
        plt.plot(IM_plot, p_fit, color=dark_colors[i], linestyle='--', label=f'Fitted Lognormal wsF ({duration})\nmu={mu_opt:.3f}, sigma={sigma_opt:.3f}, ws ={ws:.3f}')

# Add labels and title
plt.xlabel('IM')
plt.ylabel('Probability of Limit State of Exceedance') 
plt.title('Theoretical Fragility Curves - LAX (100 yr period)')
plt.legend(loc='best', fontsize=8)
plt.grid(True)
plt.savefig(f'studies/reliability/PLots_sens/Ws/ws_comparision_plot.png')

# Show the plot
plt.show()