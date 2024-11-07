import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gumbel_r as gumbel
import os

np.random.seed(0)

intensity_data = pd.read_csv(os.path.join(os.getcwd(),'studies/reliability/SB performance period fragility risk hazard/final_precipitation_intensity_data.xlsx - Sheet1.csv'))

# Extract valid data, removing NaNs
x_point_25hr = intensity_data['x0.25hr'].dropna()
x_point_50hr = intensity_data['x0.50hr'].dropna()
x1hr = intensity_data['x1hr'].dropna()
x2hr = intensity_data['x2hr'].dropna()
x3hr = intensity_data['x3hr'].dropna()

# Define the list of window sizes (moving window sizes)
window_sizes = [2, 25, 50, 100]
line_styles = ['-', '--', ':', '-.']

columns = ['x0.25hr', 'x0.50hr', 'x1hr', 'x2hr', 'x3hr']

# Create an empty DataFrame with the specified rows and columns
df_collect_pf_exceedance = pd.DataFrame(index=window_sizes, columns=columns)

# Define colors for each duration
colors = ['r', 'g', 'b', 'c', 'm', 'k']

# Define the function to fit and plot GEV
def duration_window_wise(window_size, duration, color, linestyle):
    mean_data = pd.read_csv(f'studies/reliability/SB i_mean design performance period fragility risk hazard/mean_intensity_conf_interval.csv')
    _95thpercentile_data = pd.read_csv(f'studies/reliability/SB i_mean design performance period fragility risk hazard/95_percent_intensity_conf_interval.csv')

    mean_data.set_index('Dur (hr)', inplace=True)
    mean = mean_data.loc[duration, str(window_size) + 'yr']

    _95thpercentile_data.set_index('Dur (hr)', inplace=True)
    percentile_95 = _95thpercentile_data.loc[duration, str(window_size) + 'yr']
    
    gamma = 0.57721
    beta = (mean - percentile_95) / (np.log(-np.log(0.95)) + gamma)
    mu = mean - beta * gamma

    x = np.linspace(0, 6, 1000)
    one_minus_cdf = 1 - gumbel.cdf(x, loc=mu, scale=beta)
    plt.plot(x, one_minus_cdf, label=f'{duration} - {window_size} yrs', color=color, linestyle=linestyle)

# Create a single figure for all plots
plt.figure(figsize=(12, 8))

# Loop through each window size and call the plotting function
for (i, duration) in enumerate(columns):
    for (j, window_size) in enumerate(window_sizes):
        color = colors[i % len(colors)]
        linestyle = line_styles[j % len(line_styles)]
        duration_window_wise(window_size, duration, color, linestyle)

plt.title('Intensity Exceedance Plot for All Durations and Return Periods')
plt.xlabel('Intensity Measure')
plt.ylabel('P(IM > i)')
plt.legend()
# plt.grid(True)  # Remove or comment out this line to remove the grid
plt.savefig('studies/reliability/SB i_mean design performance period fragility risk hazard/SB hazard return period/imexceedance.png')
plt.show()