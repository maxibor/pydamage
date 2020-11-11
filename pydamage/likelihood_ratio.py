import numpy as np
from scipy.stats import chi2


def LR(L0, L1, df):
    """Performs likelihood ratio test

    Args:
        L0 (np array): Log likelihood of null model
        L1 (np array): Log likelihood of alternative model
        df (int): chi2 degrees of freedom: number of constrained params
    Returns:
        float: LR test statistic lambda
        float: pvalue
    """

    LR_lambda = -2 * (L0.sum() - L1.sum())
    pval = 1 - chi2.cdf(LR_lambda, df)
    return(LR_lambda, pval)
