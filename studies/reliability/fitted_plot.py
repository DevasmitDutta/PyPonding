import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import norm
import os
import matplotlib.colors as mcolors
from mle_fit import fn_mle_pc
from scipy.stats import norm, binom

durations = ['x1hr', 'x2hr', 'x3hr', 'x6hr', 'x12hr', 'x24hr']
num_durations = len(durations)

# Define a dark color palette
dark_colors = list(mcolors.TABLEAU_COLORS.values())[:num_durations]

data_2 = pd.read_csv(os.path.join(os.getcwd(), 'studies/reliability/final_precipitation_intensity_data.xlsx - Sheet1.csv'))

plt.figure(figsize=(20, 12))

for i, duration in enumerate(durations):
    # Load the data from the CSV file
    data_1 = pd.read_csv(os.path.join(os.getcwd(), f'studies/reliability/Fragility_trials_{duration}/fragility_data_points1976.csv'))

    # Extract the durations and intensities
    intensities = data_2[duration][~np.isnan(data_2[duration])].sort_values(ascending=True).values

    # Assume the CSV file has columns named 'IM' and 'num_collapse' and 'num_gms'
    IM = intensities
    num_collapse = data_1[duration][~np.isnan(data_1[duration])].values * 200
    # num_gms = np.ones(len(data_1[duration].values),1) * 200
    num_gms = np.full(len(data_1[duration][~np.isnan(data_1[duration])].values), 200)

    # Define the log-likelihood function as per Equation (1) in Baker (2015)
    def log_likelihood(params, IM, num_collapse, num_gms):
        # ** Penalize any negative beta with a very large loglik value **
        if params[1] < 0:
            return 1e10
    
        mu, sigma = params
        p = norm.cdf(np.log(IM), loc=mu, scale=sigma)
        
        # # Replace zeros in p with a small positive number to avoid log(0)
        # p[p == 0] = np.finfo(float).tiny
        
        # log_lik = -np.sum(num_collapse * np.log(p) + (num_gms - num_collapse) * np.log(1 - p))

        # Likelihood of observing num_collapse collapses, given num_gms observations
        likelihood = binom.pmf(num_collapse, num_gms, p)

        # ** Cannot have zero likelihood value, so replace every zero likelihood 
        # value with the smallest positive normalized fixed-point value **
        likelihood[likelihood == 0] = np.finfo(float).tiny
        
        # Sum negative log likelihood (we take the negative value because we want
        # the maximum log likelihood, and the function is searching for a minimum)
        log_lik = -np.sum(np.log(likelihood))

        return log_lik

    # Initial guess for parameters [mu, sigma]
    initial_guess = [np.mean(np.log(IM)), np.std(np.log(IM))]

    # Minimize the negative log-likelihood
    result = minimize(log_likelihood, initial_guess, args=(IM, num_collapse, num_gms), method='Nelder-Mead')
    mu_opt, sigma_opt = result.x

    # Use the fn_mle_pc function to get the optimized parameters
    # mu_opt, sigma_opt = fn_mle_pc(np.transpose(IM), np.transpose(num_gms), np.transpose(num_collapse))

    # Scatter plot of observed data points with smaller triangle markers
    plt.scatter(IM, num_collapse / num_gms, color=dark_colors[i], marker='^', s=50, label=f'Observed Data ({duration})')

    # Generate IM values for plotting the fitted curve
    IM_plot = np.linspace(min(IM), 3, 500)
    p_fit = norm.cdf(np.log(IM_plot), loc=mu_opt, scale=sigma_opt)

    # Plot the fitted fragility curve
    plt.plot(IM_plot, p_fit, color=dark_colors[i], linestyle='--', label=f'Fitted Lognormal CDF ({duration})\nmu={mu_opt:.3f}, sigma={sigma_opt:.3f}')

# Add labels and title
plt.xlabel('IM')
plt.ylabel('Probability of Collapse')
plt.title('Fragility Function Fitting')
plt.legend(loc='best', fontsize=8)
plt.grid(True)
plt.savefig('studies/reliability/Fragility_fit/fitted_plot.png')

# Show the plot
plt.show()