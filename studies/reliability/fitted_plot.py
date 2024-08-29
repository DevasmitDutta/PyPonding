import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import norm
import os
import matplotlib.colors as mcolors
from mle_fit import fn_mle_pc
from scipy.stats import norm, binom

def fragility_fit(choice):
    if choice == 'LAX':
       
        durations = ['x1hr','x2hr','x3hr','x6hr','x12hr','x24hr']
        num_durations = len(durations)

        # Define a dark color palette
        dark_colors = list(mcolors.TABLEAU_COLORS.values())[:num_durations]

        data_2 = pd.read_csv(os.path.join(os.getcwd(), 'studies/reliability/LAX/final_precipitation_intensity_data.xlsx - Sheet1.csv'))

        for i, duration in enumerate(durations):

            if duration == 'x1hr':

                plt.figure(figsize=(20, 12))

                # Load the data from the CSV file
                data_1 = pd.read_csv(os.path.join(os.getcwd(), f'studies/reliability/LAX/Fragility_trials_{duration}/fragility_data_points1976.csv'))

                # Extract the durations and intensities
                intensities = data_2[duration][~np.isnan(data_2[duration])].sort_values(ascending=True).values

                # Assume the CSV file has columns named 'IM' and 'num_collapse' and 'num_gms'
                IM = intensities
                num_collapse = data_1[duration][~np.isnan(data_1[duration])].values * 200
                # num_gms = np.ones(len(data_1[duration].values),1) * 200
                num_gms = np.full(len(data_1[duration][~np.isnan(data_1[duration])].values), 200)

                # Use the fn_mle_pc function to get the optimized parameters
                mu_opt, sigma_opt = fn_mle_pc(np.transpose(IM), np.transpose(num_gms), np.transpose(num_collapse))

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
                plt.title('Theoritical Fragility Function Fitting (1hr duration)- LAX Airport')
                plt.legend(loc='best', fontsize=8)
                plt.grid(True)
                plt.savefig('studies/reliability/LAX/Fragility_fit/fitted_plot_1hr.png')

                # Show the plot
                plt.show()  

        plt.figure(figsize=(20, 12))

        for i, duration in enumerate(durations):
                
            # Load the data from the CSV file
            data_1 = pd.read_csv(os.path.join(os.getcwd(), f'studies/reliability/LAX/Fragility_trials_{duration}/fragility_data_points1976.csv'))

            # Extract the durations and intensities
            intensities = data_2[duration][~np.isnan(data_2[duration])].sort_values(ascending=True).values

            # Assume the CSV file has columns named 'IM' and 'num_collapse' and 'num_gms'
            IM = intensities
            num_collapse = data_1[duration][~np.isnan(data_1[duration])].values * 200
            # num_gms = np.ones(len(data_1[duration].values),1) * 200
            num_gms = np.full(len(data_1[duration][~np.isnan(data_1[duration])].values), 200)

            # Use the fn_mle_pc function to get the optimized parameters
            mu_opt, sigma_opt = fn_mle_pc(np.transpose(IM), np.transpose(num_gms), np.transpose(num_collapse))

            # Scatter plot of observed data points with smaller triangle markers
            # plt.scatter(IM, num_collapse / num_gms, color=dark_colors[i], marker='^', s=50, label=f'Observed Data ({duration})')

            # Generate IM values for plotting the fitted curve
            IM_plot = np.linspace(min(IM), 3, 500)
            p_fit = norm.cdf(np.log(IM_plot), loc=mu_opt, scale=sigma_opt)

            # Plot the fitted fragility curve
            plt.plot(IM_plot, p_fit, color=dark_colors[i], linestyle='--', label=f'Fitted Lognormal CDF ({duration})\nmu={mu_opt:.3f}, sigma={sigma_opt:.3f}')

        # Add labels and title
        plt.xlabel('IM')
        plt.ylabel('Probability of Collapse')
        plt.title('Theoritical Fragility Curves - LAX Airport')
        plt.legend(loc='best', fontsize=8)
        plt.grid(True)
        plt.savefig('studies/reliability/LAX/Fragility_fit/duration_comparision_plot.png')

        # Show the plot
        plt.show()

    elif choice == 'SB':
        durations = ['x0.25hr','x0.50hr','x0.75hr','x1hr','x2hr','x3hr']
        num_durations = len(durations)

        # Define a dark color palette
        dark_colors = list(mcolors.TABLEAU_COLORS.values())[:num_durations]

        data_2 = pd.read_csv(os.path.join(os.getcwd(), 'studies/reliability/SB/final_precipitation_intensity_data.xlsx - Sheet1.csv'))

        for i, duration in enumerate(durations):

            if duration == 'x1hr':

                plt.figure(figsize=(20, 12))

                # Load the data from the CSV file
                data_1 = pd.read_csv(os.path.join(os.getcwd(), f'studies/reliability/SB/Fragility_trials_{duration}/fragility_data_points1976.csv'))

                # Extract the durations and intensities
                intensities = data_2[duration][~np.isnan(data_2[duration])].sort_values(ascending=True).values

                # Assume the CSV file has columns named 'IM' and 'num_collapse' and 'num_gms'
                IM = intensities
                num_collapse = data_1[duration][~np.isnan(data_1[duration])].values * 200
                # num_gms = np.ones(len(data_1[duration].values),1) * 200
                num_gms = np.full(len(data_1[duration][~np.isnan(data_1[duration])].values), 200)

                # Use the fn_mle_pc function to get the optimized parameters
                mu_opt, sigma_opt = fn_mle_pc(np.transpose(IM), np.transpose(num_gms), np.transpose(num_collapse))

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
                plt.title('Fragility Function Fitting (1 hr duration)- Santa Barbara')
                plt.legend(loc='best', fontsize=8)
                plt.grid(True)
                plt.savefig('studies/reliability/SB/Fragility_fit/fitted_plot_1hr.png')

                # Show the plot
                plt.show()  

        plt.figure(figsize=(20, 12))

        for i, duration in enumerate(durations):
                
                # Load the data from the CSV file
            data_1 = pd.read_csv(os.path.join(os.getcwd(), f'studies/reliability/SB/Fragility_trials_{duration}/fragility_data_points1976.csv'))

            # Extract the durations and intensities
            intensities = data_2[duration][~np.isnan(data_2[duration])].sort_values(ascending=True).values

            # Assume the CSV file has columns named 'IM' and 'num_collapse' and 'num_gms'
            IM = intensities
            num_collapse = data_1[duration][~np.isnan(data_1[duration])].values * 200
            # num_gms = np.ones(len(data_1[duration].values),1) * 200
            num_gms = np.full(len(data_1[duration][~np.isnan(data_1[duration])].values), 200)

            # Use the fn_mle_pc function to get the optimized parameters
            mu_opt, sigma_opt = fn_mle_pc(np.transpose(IM), np.transpose(num_gms), np.transpose(num_collapse))

            # # Scatter plot of observed data points with smaller triangle markers
            # plt.scatter(IM, num_collapse / num_gms, color=dark_colors[i], marker='^', s=50, label=f'Observed Data ({duration})')

            # Generate IM values for plotting the fitted curve
            IM_plot = np.linspace(min(IM), 3, 500)
            p_fit = norm.cdf(np.log(IM_plot), loc=mu_opt, scale=sigma_opt)

            # Plot the fitted fragility curve
            plt.plot(IM_plot, p_fit, color=dark_colors[i], linestyle='--', label=f'Fitted Lognormal CDF ({duration})\nmu={mu_opt:.3f}, sigma={sigma_opt:.3f}')

        # Add labels and title
        plt.xlabel('IM')
        plt.ylabel('Probability of Collapse')
        plt.title('Theoritical Fragility Curves - Santa Barbara')
        plt.legend(loc='best', fontsize=8)
        plt.grid(True)
        plt.savefig('studies/reliability/SB/Fragility_fit/duration_comparision_plot.png')

        # Show the plot
        plt.show()

fragility_fit('SB')