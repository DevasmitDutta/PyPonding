import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import genextreme as gev
from scipy.stats import ecdf
import os

intensity_data = pd.read_csv(os.path.join(os.getcwd(),'studies/reliability/SB 2/final_precipitation_intensity_data.xlsx - Sheet1.csv'));

# Extract valid data, removing NaNs
x_point25hr = intensity_data['x0.25hr'].dropna()
x_point50hr = intensity_data['x0.50hr'].dropna()
x_point75hr = intensity_data['x0.75hr'].dropna()
x1hr = intensity_data['x1hr'].dropna()
x2hr = intensity_data['x2hr'].dropna()
x3hr = intensity_data['x3hr'].dropna()

f1, sorted_data_1 = ecdf(x_point25hr).cdf.probabilities, np.unique(x_point25hr)
f2, sorted_data_2 = ecdf(x_point50hr).cdf.probabilities, np.unique(x_point50hr)
f3, sorted_data_3 = ecdf(x_point75hr).cdf.probabilities, np.unique(x_point75hr)
f4, sorted_data_4 = ecdf(x1hr).cdf.probabilities, np.unique(x1hr)
f5, sorted_data_5 = ecdf(x2hr).cdf.probabilities, np.unique(x2hr)
f6, sorted_data_6 = ecdf(x3hr).cdf.probabilities, np.unique(x3hr)      

# Define the list of window sizes (moving window sizes)
window_sizes = [20, 30, 40]  # Example window sizes
line_styles = ['-', '--', ':', '-.', '.-.']  # Different line styles for each window size

# Define colors for each duration
colors = ['r', 'g', 'b', 'c', 'm', 'k']  # Red, Green, Blue, Cyan, Magenta, Black

# Create figure to plot everything
plt.figure(figsize=(15, 10))
plt.scatter(sorted_data_1, 1 - f1, marker='^', color='r', label='Observed - 0.25 hr')
plt.scatter(sorted_data_2, 1 - f2, marker='^', color='g', label='Observed - 0.50 hr')
plt.scatter(sorted_data_3, 1 - f3, marker='^', color='b', label='Observed - 0.75 hr')
plt.scatter(sorted_data_4, 1 - f4, marker='^', color='c', label='Observed - 1 hr')
plt.scatter(sorted_data_5, 1 - f5, marker='^', color='m', label='Observed - 2 hr')
plt.scatter(sorted_data_6, 1 - f6, marker='^', color='k', label='Observed - 3 hr')

# Define the function to fit and plot GEV
def plot_window_wise(window_size, x1hr, x2hr, x3hr, x6hr, x12hr, x24hr, line_style, colors):
    max_index_others = len(x1hr)
    max_index_24hr = len(x24hr)
    
    gev_params_1, gev_params_2, gev_params_3 = [], [], []
    gev_params_4, gev_params_5, gev_params_6 = [], [], []
    
    for i in range(max_index_others - window_size + 1):
        window_data_1 = x1hr.iloc[i:i + window_size]
        window_data_2 = x2hr.iloc[i:i + window_size]
        window_data_3 = x3hr.iloc[i:i + window_size]
        window_data_4 = x6hr.iloc[i:i + window_size]
        window_data_5 = x12hr.iloc[i:i + window_size]

        params_window_1 = gev.fit(window_data_1)
        params_window_2 = gev.fit(window_data_2)
        params_window_3 = gev.fit(window_data_3)
        params_window_4 = gev.fit(window_data_4)
        params_window_5 = gev.fit(window_data_5)

        gev_params_1.append(params_window_1)
        gev_params_2.append(params_window_2)
        gev_params_3.append(params_window_3)
        gev_params_4.append(params_window_4)
        gev_params_5.append(params_window_5)

    for i in range(max_index_24hr - window_size + 1):
        window_data_6 = x24hr.iloc[i:i + window_size]

        params_window_6 = gev.fit(window_data_6)

        gev_params_6.append(params_window_6)

    x = np.linspace(0, max(np.max(sorted_data_1), np.max(sorted_data_2), np.max(sorted_data_3),
                           np.max(sorted_data_4), np.max(sorted_data_5), np.max(sorted_data_6)), 1000)

    # Plot fitted GEV for each window size
    plt.plot(x, 1 - gev.cdf(x, *np.mean(gev_params_1, axis=0)), line_style, color=colors[0], label=f'0.25 hr - perform period {window_size} Mean')
    plt.plot(x, 1 - gev.cdf(x, *np.mean(gev_params_2, axis=0)), line_style, color=colors[1], label=f'0.50 hr - perform period {window_size} Mean')
    plt.plot(x, 1 - gev.cdf(x, *np.mean(gev_params_3, axis=0)), line_style, color=colors[2], label=f'0.75 hr - perform period {window_size} Mean')
    plt.plot(x, 1 - gev.cdf(x, *np.mean(gev_params_4, axis=0)), line_style, color=colors[3], label=f'1 hr - perform period {window_size} Mean')
    plt.plot(x, 1 - gev.cdf(x, *np.mean(gev_params_5, axis=0)), line_style, color=colors[4], label=f'2 hr - perform period {window_size} Mean')
    plt.plot(x, 1 - gev.cdf(x, *np.mean(gev_params_6, axis=0)), line_style, color=colors[5], label=f'3 hr - perform period {window_size} Mean')

# Loop through each window size and call the plotting function
for i, window_size in enumerate(window_sizes):
    line_style = line_styles[i]
    plot_window_wise(window_size, x_point25hr, x_point50hr, x_point75hr, x1hr, x2hr, x3hr, line_style, colors)

# Finalize figure with labels and title
plt.xlabel('Intensities, i (inc/hr)')
plt.ylabel('P(IM > i)')
plt.title('Probability of Intensity Exceedance for Different Durations and Performance Periods (Santa Barbara)')
plt.legend()
plt.savefig('studies/reliability/SB hazard/movingwindowcombined.png')    
