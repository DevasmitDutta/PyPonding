import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import genextreme as gev
import os
import scipy.stats as stats
from scipy.stats import norm

np.random.seed(0)
# exec(open('/home/devasmit/Desktop/PyPonding/studies/reliability/LAX/hazard/fragility_assessment LAX-2.py').read())
intensity_data = pd.read_csv(os.path.join(os.getcwd(),'studies/reliability/SB copy/final_precipitation_intensity_data.xlsx - Sheet1.csv'));

# Extract valid data, removing NaNs
x_point_25hr = intensity_data['x0.25hr'].dropna()
x_point_50hr = intensity_data['x0.50hr'].dropna()
x_point_75hr = intensity_data['x0.75hr'].dropna()
x1hr = intensity_data['x1hr'].dropna()
x2hr = intensity_data['x2hr'].dropna()
x3hr = intensity_data['x3hr'].dropna()

label = 1976

# Calculate ECDF for each duration
def ecdf(data):
    """Compute empirical cumulative distribution function (ECDF)."""
    x = np.sort(data)
    y = np.arange(1, len(x) + 1) / len(x)
    return y, x

f1, sorted_data_1 = ecdf(x_point_25hr)
f2, sorted_data_2 = ecdf(x_point_50hr)
f3, sorted_data_3 = ecdf(x_point_75hr)
f4, sorted_data_4 = ecdf(x1hr)
f5, sorted_data_5 = ecdf(x2hr)
f6, sorted_data_6 = ecdf(x3hr)

# Define the list of window sizes (moving window sizes)
window_sizes = [100] #np.arange(5,41,1) #[5, 8, 10, 15, 20, 30, 40, 50, 60]  # Generates values from 1 to 10 with a step size of 1
line_styles = ['-', '--', ':', '-.', '.-.']  # Different line styles for each window size

columns = ['x0.25hr', 'x0.50hr', 'x1hr', 'x2hr', 'x3hr']

# Create an empty DataFrame with the specified rows and columns
df_collect_pf_exceedance = pd.DataFrame(index=window_sizes, columns=columns)

# Define colors for each duration
colors = ['r', 'g', 'b', 'c', 'm', 'k']  # Red, Green, Blue, Cyan, Magenta, Black

# Create figure to plot everything
plt.figure(figsize=(10, 6))
# plt.scatter(sorted_data_1, f1, marker='^', color='r', label='Observed - 1 hr')
# plt.scatter(sorted_data_2, f2, marker='^', color='g', label='Observed - 2 hr')
# plt.scatter(sorted_data_3, f3, marker='^', color='b', label='Observed - 3 hr')
# plt.scatter(sorted_data_4, f4, marker='^', color='c', label='Observed - 6 hr')
# plt.scatter(sorted_data_5, f5, marker='^', color='m', label='Observed - 12 hr')
# plt.scatter(sorted_data_6, f6, marker='^', color='k', label='Observed - 24 hr')

# Define the function to fit and plot GEV
def duration_window_wise(window_size, duration):
    # if duration == 'x0.50hr':
    #   max_index = len(x_point_50hr)
    #   gev_params = []
    # elif duration == 'x0.25hr' or duration == 'x0.75hr':
    #   max_index = len(x_point_25hr)
    #   gev_params = []
    # else:
    #   max_index = len(x3hr)
    #   gev_params = [] 
    
    # data = intensity_data[duration].dropna()

    # for i in range(max_index - window_size + 1):
    #     window_data = data.iloc[i:i + window_size]

    #     params_window = gev.fit(window_data)

    #     gev_params.append(params_window)

    # x = np.linspace(0, max(np.max(sorted_data_1), np.max(sorted_data_2), np.max(sorted_data_3),
    #                        np.max(sorted_data_4), np.max(sorted_data_5), np.max(sorted_data_6)), 1000)

    # mean_gev_params = np.mean(gev_params, axis=0)
        # Assuming mean_gev_params_1 contains the shape, location, and scale parameters
    # shape, loc, scale = mean_gev_params

    num_samples = 1000
    # Generate 100 random samples from the GEV distribution
    # random_samples = np.sort(stats.genextreme.rvs(c=shape, loc=loc, scale=scale, size=num_samples))
    # print(f'Random samples: {random_samples}')

    mean_data = pd.read_csv('studies/reliability/SB copy/mean_intensity_conf_interval.csv')
    _95thpercentile_data = pd.read_csv('studies/reliability/LAX/95_percent_intensity_conf_interval.csv')

    mean_data.set_index('Dur (hr)', inplace=True)
    mean = mean_data.loc[duration, window_size]

    _95thpercentile_data.set_index('Dur (hr)', inplace=True)
    percentile_95 = _95thpercentile_data.loc[duration, window_size]
    
    gamma = 0.57721
    beta = (mean - percentile_95) / (np.log(-np.log(0.95)) + gamma)
    mu = mean - beta * gamma
    random_samples = np.sort(gumbel.rvs(loc=mu, scale=beta, size=num_samples))
    
    Pf = np.zeros(1000)
    threshold_num = np.random.uniform(0, 1, num_samples)
    cnt = 0
    for (i, IM) in enumerate(random_samples):
        print(f'Entering IM count: {i} with IM: {IM}')
        # Pf[i] = run_analysis('W14X22',360.0,0.020833333333333332,6.944444444444444e-05, label, IM, duration, window_size)

        # Read the CSV file into a DataFrame
        df = pd.read_csv('studies/reliability/SB copy/Fragility_fit/mu_sigma_opt_durationwise.csv')

        # Extract the optimized mu and sigma values for the current duration
        mu_opt = df.loc[df['duration'].str.strip() == duration, 'mu_opt'].values[0]
        sigma_opt = df.loc[df['duration'].str.strip() == duration, 'sigma_opt'].values[0]

        Pf[i] = norm.cdf(np.log(IM), loc=mu_opt, scale=sigma_opt)
        # num = np.random.uniform(0, 1)
        print(f'Probability of exceedance: {Pf[i]} < Random number: {threshold_num[i]}')

        if Pf[i] < threshold_num[i]:
            print('True')
            cnt += 1  
        else:
            print('False')

    #  print(f'shape: {mean_gev_params[0]}, loc: {mean_gev_params[1]}, scale: {mean_gev_params[2]} for {duration} and {window_size} years')

    # plt.hist(threshold_num, bins = 100, edgecolor='black', alpha=0.7)
    # plt.title('Histogram of Uniform Samples')
    # plt.xlabel('Value')
    # plt.ylabel('Frequency')
    # plt.show()

    # plt.hist(random_samples, edgecolor='black', bins = 100, alpha=0.7, label=f'{window_size} yrs, shape: {mean_gev_params[0]:.2f}, loc: {mean_gev_params[1]:.2f}, scale: {mean_gev_params[2]:.2f}')
    # plt.hist(random_samples, edgecolor='black', bins = 100, alpha=0.7, label=f'{window_size} yrs, shape: {0:.2f}, loc: {mu:.2f}, scale: {beta:.2f}')

     
    Prob_of_exceedance = cnt/num_samples
    print(f'Probability of exceedance: {Prob_of_exceedance} in {window_size} years over {duration} duration') 
    df_collect_pf_exceedance.loc[window_size, duration] = Prob_of_exceedance 


# Loop through each window size and call the plotting function
for (i, window_size) in enumerate(window_sizes):
    for (j, duration) in enumerate(columns):
        duration_window_wise(window_size, duration)

# plt.title('Generated random samples from GEV distribution for x1hr')
# plt.xlabel('Intensity Measure')
# plt.ylabel('Frequency')
# plt.legend()
# plt.savefig('studies/reliability/SB copy/Fragility_fit/Random_samples_from_GEV_distribution_x1hr.png')

plt.figure(figsize=(10, 6))
# Finalize figure with labels and title
for (i, duration) in enumerate(columns):
    plt.plot(window_sizes, df_collect_pf_exceedance[duration], label = f'{duration}', color = colors[i])

plt.xlabel('Return Period, years')
plt.ylabel('Probability of Limit State of Exceedance')
plt.title('Probability of Limit State of Exceedance in X years over different durations')
plt.legend()
plt.savefig('studies/reliability/SB copy/Fragility_fit/Probability_of_limit_state_exceedance_X_years.png')