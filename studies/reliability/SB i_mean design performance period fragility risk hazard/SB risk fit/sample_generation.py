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
intensity_data = pd.read_csv(os.path.join(os.getcwd(),'studies/reliability/SB i_mean design performance period fragility risk hazard/final_precipitation_intensity_data.xlsx - Sheet1.csv'));

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
window_sizes = [10, 15, 20, 30, 40] #np.arange(5,41,1) #[5, 8, 10, 15, 20, 30, 40, 50, 60]  # Generates values from 1 to 10 with a step size of 1
line_styles = ['-', '--', ':', '-.', '.-.']  # Different line styles for each window size

columns = ['x0.25hr', 'x0.50hr', 'x1hr', 'x2hr', 'x3hr']
# columns = ['x1hr']    

# Create an empty DataFrame with the specified rows and columns
df_collect_pf_exceedance = pd.DataFrame(index=window_sizes, columns=columns)

# Define colors for each duration
colors = ['r', 'g', 'b', 'c', 'm', 'k']  # Red, Green, Blue, Cyan, Magenta, Black

# Define the function to fit and plot GEV
def duration_window_wise(window_size, duration):
    if duration == 'x0.50hr':
      max_index = len(x_point_50hr)
      gev_params = []
    elif duration == 'x0.25hr' or duration == 'x0.75hr':
      max_index = len(x_point_25hr)
      gev_params = []
    else:
      max_index = len(x3hr)
      gev_params = [] 
    
    data = intensity_data[duration].dropna()

    for i in range(max_index - window_size + 1):
        window_data = data.iloc[i:i + window_size]

        params_window = gev.fit(window_data)

        gev_params.append(params_window)

    mean_gev_params = np.mean(gev_params, axis=0)
        # Assuming mean_gev_params_1 contains the shape, location, and scale parameters
    shape, loc, scale = mean_gev_params

    num_samples = 1000
    # Generate 100 random samples from the GEV distribution
    random_samples = np.sort(stats.genextreme.rvs(c=shape, loc=loc, scale=scale, size=num_samples))

    plt.hist(random_samples, edgecolor='black', bins = 20, alpha=0.7, label=f'{window_size} yrs, shape: {mean_gev_params[0]:.2f}, loc: {mean_gev_params[1]:.2f}, scale: {mean_gev_params[2]:.2f}')

# Loop through each window size and call the plotting function
for (i, duration) in enumerate(columns):
    plt.figure(figsize=(10, 6))

    for (j, window_size) in enumerate(window_sizes):
        duration_window_wise(window_size, duration)

    plt.title(f'Generated random samples from GEV for different performance periods for {duration}')
    plt.xlabel('Intensity Measure')
    plt.ylabel('Frequency')
    plt.legend()
    plt.savefig(f'studies/reliability/SB i_mean design performance period fragility risk hazard/SB risk fit/Random_samples_from_GEV_distribution_{duration}.png')    

for (i, duration) in enumerate(columns):

    for (j, window_size) in enumerate(window_sizes):
        plt.figure(figsize=(10, 6))
        duration_window_wise(window_size, duration)

        plt.title(f'Generated random samples from GEV for different {window_size} yr period {duration}')
        plt.xlabel('Intensity Measure')
        plt.ylabel('Frequency')
        plt.legend()
        plt.savefig(f'studies/reliability/SB i_mean design performance period fragility risk hazard/SB risk fit/Random_samples_from_GEV_distribution_{window_size}_yr_{duration}.png')    