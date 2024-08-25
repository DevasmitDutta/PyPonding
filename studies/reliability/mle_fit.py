import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm, binom

def fn_mle_pc(IM, num_gms, num_collapse):
    """
    This function fits a lognormal CDF to observed probability of collapse 
    data using optimization on the likelihood function for the data.
    
    Inputs:
    IM            1xn           IM levels of interest
    num_gms       1x1 or 1xn    number of ground motions used at each IM level
    num_collapse 	1xn           number of collapses observed at each IM level
    
    Outputs:
    theta         1x1           median of fragility function
    beta          1x1           lognormal standard deviation of fragility function
    """

    # Initial guess for the fragility function parameters theta and beta
    # ** Use method of moments **
    x0 = [np.mean(np.log(IM)), np.std(np.log(IM))]

    # Run optimization
    result = minimize(mlefit, x0, args=(num_gms, num_collapse, IM), 
                      options={'maxiter': 1000}, method='Nelder-Mead')

    # Extract optimized parameters
    theta = np.exp(result.x[0])  # return theta in linear space
    beta = result.x[1]

    return theta, beta

def mlefit(params, num_gms, num_collapse, IM):
    """
    Objective function to be optimized in the minimization process.
    
    Inputs:
    params        list          current guess of the parameters [log(theta), beta]
    num_gms       1x1 or 1xn    number of ground motions used at each IM level
    num_collapse 	1xn           number of collapses observed at each IM level
    IM            1xn           IM levels of interest
    
    Outputs:
    loglik        1x1           negative log-likelihood to be minimized
    """

    # ** Penalize any negative beta with a very large loglik value **
    if params[1] < 0:
        return 1e10
    
    # Estimated probabilities of collapse, given the current fragility function
    # parameter estimates
    p = norm.cdf(np.log(IM), loc=params[0], scale=params[1])

    # Likelihood of observing num_collapse(i) collapses, given num_gms
    # observations, using the current fragility function parameter estimates
    likelihood = binom.pmf(num_collapse, num_gms, p)

    # ** Cannot have zero likelihood value, so replace every zero likelihood 
    # value with the smallest positive normalized fixed-point value **
    likelihood[likelihood == 0] = np.finfo(float).tiny
    
    # Sum negative log likelihood (we take the negative value because we want
    # the maximum log likelihood, and the function is searching for a minimum)
    loglik = -np.sum(np.log(likelihood))

    return loglik
