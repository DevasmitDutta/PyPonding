import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import genextreme as gev
import os
import scipy.stats as stats
from scipy.stats import norm

np.random.seed(0)
exec(open('/home/devasmit/Desktop/PyPonding/studies/reliability/LAX/hazard/fragility_assessment LAX-2.py').read())
intensity_data = pd.read_csv(os.path.join(os.getcwd(),'studies/reliability/LAX/final_precipitation_intensity_data.xlsx - Sheet1.csv'));

# Extract valid data, removing NaNs
x1hr = intensity_data['x1hr'].dropna()
x2hr = intensity_data['x2hr'].dropna()
x3hr = intensity_data['x3hr'].dropna()
x6hr = intensity_data['x6hr'].dropna()
x12hr = intensity_data['x12hr'].dropna()
x24hr = intensity_data['x24hr'].dropna()

label = 1976

# Calculate ECDF for each duration
def ecdf(data):
    """Compute empirical cumulative distribution function (ECDF)."""
    x = np.sort(data)
    y = np.arange(1, len(x) + 1) / len(x)
    return y, x

f1, sorted_data_1 = ecdf(x1hr)
f2, sorted_data_2 = ecdf(x2hr)
f3, sorted_data_3 = ecdf(x3hr)
f4, sorted_data_4 = ecdf(x6hr)
f5, sorted_data_5 = ecdf(x12hr)
f6, sorted_data_6 = ecdf(x24hr)

# Define the list of window sizes (moving window sizes)
window_sizes = np.arange(5,61,1) #[5, 8, 10, 15, 20, 30, 40, 50, 60]  # Generates values from 1 to 10 with a step size of 1
line_styles = ['-', '--', ':', '-.', '.-.']  # Different line styles for each window size

columns = ['x1hr', 'x2hr', 'x3hr', 'x6hr', 'x12hr', 'x24hr']

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
    if duration != 'x24hr':
      max_index = len(x1hr)
      gev_params = []
    else:
      max_index = len(x24hr)
      gev_params = []
    
    data = intensity_data[duration].dropna()

    for i in range(max_index - window_size + 1):
        window_data = data.iloc[i:i + window_size]

        params_window = gev.fit(window_data)

        gev_params.append(params_window)

    x = np.linspace(0, max(np.max(sorted_data_1), np.max(sorted_data_2), np.max(sorted_data_3),
                           np.max(sorted_data_4), np.max(sorted_data_5), np.max(sorted_data_6)), 1000)

    mean_gev_params = np.mean(gev_params, axis=0)
        # Assuming mean_gev_params_1 contains the shape, location, and scale parameters
    shape, loc, scale = mean_gev_params

    # Generate 100 random samples from the GEV distribution
    random_samples = np.sort(stats.genextreme.rvs(c=shape, loc=loc, scale=scale, size=1000))
    # print(f'Random samples: {random_samples}')
    
    Pf = np.zeros(1000)
    cnt = 0
    for (i, IM) in enumerate(random_samples):
        print(f'Entering IM count: {i} with IM: {IM}')
        # Pf[i] = run_analysis('W14X22',360.0,0.020833333333333332,6.944444444444444e-05, label, IM, duration, window_size)

        # Read the CSV file into a DataFrame
        df = pd.read_csv('studies/reliability/LAX/Fragility_fit/mu_sigma_opt_durationwise.csv')

        # Extract the optimized mu and sigma values for the current duration
        mu_opt = df.loc[df['duration'].str.strip() == duration, 'mu_opt'].values[0]
        sigma_opt = df.loc[df['duration'].str.strip() == duration, 'sigma_opt'].values[0]

        Pf[i] = norm.cdf(np.log(IM), loc=mu_opt, scale=sigma_opt)
        num = np.random.uniform(0, 1)
        if Pf[i] < num:
            cnt += 1  

    Prob_of_exceedance = cnt/1000
    print(f'Probability of exceedance: {Prob_of_exceedance} in {window_size} years over {duration} duration') 
    df_collect_pf_exceedance.loc[window_size, duration] = Prob_of_exceedance 


# Loop through each window size and call the plotting function
for (i, window_size) in enumerate(window_sizes):
    for (j, duration) in enumerate(columns):
        duration_window_wise(window_size, duration)

# Finalize figure with labels and title
for (i, duration) in enumerate(columns):
    plt.plot(window_sizes, df_collect_pf_exceedance[duration], label = f'{duration}', color = colors[i])

plt.xlabel('Return Period, years')
plt.ylabel('Probability of Limit State of Exceedance')
plt.title('Probability of Limit State of Exceedance in X years over different durations')
plt.legend()
plt.savefig('studies/reliability/LAX/Fragility_fit/Probability_of_limit_state_exceedance_X_years.png')