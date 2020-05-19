from scipy.stats import norm
from numpy import std, sqrt, log


def vuong_closeness(LA, LB, N, pdiff):
    """Performs Vuong's closeness test

    Args:
        LA (np.array): pointwise likelihood of model A
        LB (np.array): pointwise likelihood of model B
        N (int): Number of observations
        pdiff (int): Difference in number of model paramaters (A - B)
    Returns:
        float: z-score
        float: pvalue
    """

    LR = LA.sum() - LB.sum() - (pdiff/2)*log(N)
    omega = std(LA-LB)
    Z = LR/(omega*sqrt(N))
    pval = 1 - norm.cdf(Z)
    return(Z, pval)
