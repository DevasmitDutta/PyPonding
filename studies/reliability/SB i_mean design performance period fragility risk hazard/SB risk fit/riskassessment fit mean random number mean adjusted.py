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
window_sizes = [5, 10, 15, 20, 25, 30, 35, 40] # [5, 10, 15, 20, 25, 30, 35, 40] #np.arange(5,41,1)  # Generates values from 1 to 10 with a step size of 1
line_styles = ['-', '--', ':', '-.', '.-.']  # Different line styles for each window size

columns = ['x0.25hr', 'x0.50hr', 'x1hr', 'x2hr', 'x3hr']
# columns = ['x0.25hr']

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

    num_samples = 100
    # Generate 100 random samples from the GEV distribution
    random_samples = np.sort(stats.genextreme.rvs(c=shape, loc=loc, scale=scale, size=num_samples))
    
    Pf = np.zeros(num_samples)
    cnt = 0
    for (i, IM) in enumerate(random_samples):
        print(f'Entering IM count: {i} with IM: {IM}')
        # Pf[i] = run_analysis('W14X22',360.0,0.020833333333333332,6.944444444444444e-05, label, IM, duration, window_size)

        # Read the CSV file into a DataFrame
        df = pd.read_csv(f'studies/reliability/SB i_mean design performance period fragility risk hazard/Fragility_fit/mu_sigma_opt_durationwise.csv')
        # print(df.values[0])
        # Extract the optimized mu and sigma values for the current duration
        mu_opt = df.values[0][0] #df.loc[df['duration'].str.strip() == duration, 'mu_opt'].values[0]
        sigma_opt = df.values[0][1] #df.loc[df['duration'].str.strip() == duration, 'sigma_opt'].values[0]


        Pf[i] = norm.cdf(np.log(IM), loc=mu_opt, scale=sigma_opt)
        threshold_num = np.random.uniform(0, norm.cdf(np.log(np.mean(random_samples)), loc=mu_opt, scale=sigma_opt))
        # num = np.random.uniform(0, 1)
        print(f'For duration {duration} period {window_size} Optimized mu: {mu_opt}, Optimized sigma: {sigma_opt} at sample_mean {np.mean(random_samples)}')
        print(f'Probability of exceedance: {Pf[i]} < Random number: {threshold_num}')

        # if Pf[i] < threshold_num:
        #     print('True')
        #     cnt += 1  
        # else:
        #     print('False')

        if Pf[i] > threshold_num:
            print('True')
            cnt += 1  
        else:
            print('False')
     
    Prob_of_exceedance = cnt/num_samples
    print(f'Probability of exceedance: {Prob_of_exceedance} in {window_size} years over {duration} duration') 
    df_collect_pf_exceedance.loc[window_size, duration] = Prob_of_exceedance 


# Loop through each window size and call the plotting function
for (i, duration) in enumerate(columns):
    for (j, window_size) in enumerate(window_sizes):
        duration_window_wise(window_size, duration)

plt.figure(figsize=(10, 6))
# Finalize figure with labels and title
for (i, duration) in enumerate(columns):
    plt.plot(window_sizes, df_collect_pf_exceedance[duration], label = f'{duration}', color = colors[i])

plt.xlabel('Performance Period, years')
plt.ylabel('Probability of Limit State of Exceedance')
plt.title('Probability of Limit State of Exceedance in X years over different durations')
plt.legend()
plt.savefig('studies/reliability/SB i_mean design performance period fragility risk hazard/SB risk fit/Probability_of_limit_state_exceedance_X_years_mean.png')