import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# Fragility parameters: log-mean (mu) and log-std-dev (sigma)
fragility_data = {
    '100 yr': {
        'mu': [3.364, 2.934, 3.072],
        'sigma': [0.345, 0.281, 0.263]
    },
    '500 yr': {
        'mu': [3.168, 2.942, 3.942],
        'sigma': [0.287, 0.242, 0.242]
    },
    '1000 yr': {
        'mu': [3.503, 2.79, 2.868],
        'sigma': [0.33, 0.262, 0.252]
    }
}

durations = ['0.25 hr', '1 hr', '2 hr']
intensity_range = np.linspace(5, 10, 500)  # Avoid log(0) by starting > 0

# Plot fragility curves separately for each return period
for year in fragility_data:
    plt.figure(figsize=(8, 5))
    for i, duration in enumerate(durations):
        mu = fragility_data[year]['mu'][i]
        sigma = fragility_data[year]['sigma'][i]
        prob_exceed = norm.cdf(np.log(intensity_range), loc=mu, scale=sigma)
        plt.plot(intensity_range, prob_exceed, label=f'{duration}')
    
    plt.title(f'Fragility Curves – {year}')
    plt.xlabel('Rainfall Intensity (in/hr)')
    plt.ylabel('Probability of Exceedance')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(title='Duration')
    plt.tight_layout()

plt.show()
