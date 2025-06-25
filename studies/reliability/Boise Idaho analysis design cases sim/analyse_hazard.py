import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import genextreme

# GEV parameters for each duration (shape = c, scale = σ, location = μ)
gev_params = {
    '0.25 hr': {'shape': -3.13695870563135, 'scale': 6.62381785498274e-05, 'loc': 0.400020760542928},
    '1 hr':    {'shape': -0.850094508367361, 'scale': 0.0773127978946366, 'loc': 0.231436631889813},
    '2 hr':    {'shape': -0.671148088734158, 'scale': 0.0432541766522506, 'loc': 0.177422510786878},
    '3 hr':    {'shape': -0.362127548423937, 'scale': 0.0391903509566718, 'loc': 0.147856188799176}
}

# Intensity range for plotting
x = np.linspace(0, 2, 1000)  # Adjust range based on intensity domain

# Plot CCDF (1 - CDF)
plt.figure(figsize=(8, 5))
for duration, params in gev_params.items():
    c, loc, scale = params['shape'], params['loc'], params['scale']
    ccdf = 1 - genextreme.cdf(x, c=c, loc=loc, scale=scale)
    plt.plot(x, ccdf, label=duration)

plt.title('GEV Hazard Curves (1 - CDF)')
plt.xlabel('Rainfall Intensity (in/hr)')
plt.ylabel('Exceedance Probability')
# plt.yscale('log')  # Optional: log scale to visualize small probabilities
# plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.show()
