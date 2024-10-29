import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import genextreme as gev
from scipy.stats import gumbel_r as gumbel
import os
import scipy.stats as stats
from scipy.stats import norm
from scipy.stats import ecdf

np.random.seed(0)
# exec(open('/home/devasmit/Desktop/PyPonding/studies/reliability/LAX/hazard/fragility_assessment LAX-2.py').read())
intensity_data = pd.read_csv(os.path.join(os.getcwd(),'studies/reliability/SB 2/final_precipitation_intensity_data.xlsx - Sheet1.csv'));

# Extract valid data, removing NaNs
x_point_25hr = intensity_data['x0.25hr'].dropna()
x_point_50hr = intensity_data['x0.50hr'].dropna()
# x_point_75hr = intensity_data['x0.75hr'].dropna()
x1hr = intensity_data['x1hr'].dropna()
x2hr = intensity_data['x2hr'].dropna()
x3hr = intensity_data['x3hr'].dropna()

label = 1976

f1, sorted_data_1 = ecdf(x_point_25hr).cdf.probabilities, np.unique(x_point_25hr)
f2, sorted_data_2 = ecdf(x_point_50hr).cdf.probabilities, np.unique(x_point_50hr)
# f3, sorted_data_3 = ecdf(x_point75hr).cdf.probabilities, np.unique(x_point75hr)
f4, sorted_data_4 = ecdf(x1hr).cdf.probabilities, np.unique(x1hr)
f5, sorted_data_5 = ecdf(x2hr).cdf.probabilities, np.unique(x2hr)
f6, sorted_data_6 = ecdf(x3hr).cdf.probabilities, np.unique(x3hr)  

# Define the list of window sizes (moving window sizes)
window_sizes = [2, 5, 10, 25, 50] #np.arange(5,41,1) #[5, 8, 10, 15, 20, 30, 40, 50, 60]  # Generates values from 1 to 10 with a step size of 1
line_styles = ['-', '--', ':', '-.', '.-.']  # Different line styles for each window size

columns = ['x0.25hr', 'x0.50hr', 'x1hr', 'x2hr', 'x3hr']
# columns = ['x1hr']    

# Create an empty DataFrame with the specified rows and columns
df_collect_pf_exceedance = pd.DataFrame(index=window_sizes, columns=columns)

# Define colors for each duration
colors = ['r', 'g', 'b', 'c', 'm', 'k']  # Red, Green, Blue, Cyan, Magenta, Black

# Define the function to fit and plot GEV
def duration_window_wise(window_size, duration):

    num_samples = 1000
    mean_data = pd.read_csv(f'studies/reliability/SB {window_size}/mean_intensity_conf_interval.csv')
    _95thpercentile_data = pd.read_csv(f'studies/reliability/SB {window_size}/95_percent_intensity_conf_interval.csv')

    mean_data.set_index('Dur (hr)', inplace=True)
    mean = mean_data.loc[duration, str(window_size) + 'yr']

    _95thpercentile_data.set_index('Dur (hr)', inplace=True)
    percentile_95 = _95thpercentile_data.loc[duration, str(window_size) + 'yr']
    
    gamma = 0.57721
    beta = (mean - percentile_95) / (np.log(-np.log(0.95)) + gamma)
    mu = mean - beta * gamma
    random_samples = np.sort(gumbel.rvs(loc=mu, scale=beta, size=num_samples))

    plt.hist(random_samples, edgecolor='black', bins = 100, alpha=0.7, label=f'{window_size} yrs, shape: {0:.2f}, loc: {mu:.2f}, scale: {beta:.2f}')

# Loop through each window size and call the plotting function
for (i, duration) in enumerate(columns):
    plt.figure(figsize=(10, 6))

    for (j, window_size) in enumerate(window_sizes):
        duration_window_wise(window_size, duration)

    plt.title(f'Generated random samples from GEV for different performance periods for {duration}')
    plt.xlabel('Intensity Measure')
    plt.ylabel('Frequency')
    plt.legend()
    plt.savefig(f'studies/reliability/SB risk conf interval/Random_samples_from_GEV_distribution_{duration}.png')    

for (i, duration) in enumerate(columns):

    for (j, window_size) in enumerate(window_sizes):
        plt.figure(figsize=(10, 6))
        duration_window_wise(window_size, duration)

        plt.title(f'Generated random samples from GEV for different {window_size} yr period {duration}')
        plt.xlabel('Intensity Measure')
        plt.ylabel('Frequency')
        plt.legend()
        plt.savefig(f'studies/reliability/SB risk conf interval/Random_samples_from_GEV_distribution_{window_size}_yr_{duration}.png')    