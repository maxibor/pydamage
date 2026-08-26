#!/usr/bin/env python3

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import binom, chi2
from pydamage.utils import RMSE


def optim(function, parameters: list, xdata: np.ndarray, ydata: np.ndarray, bounds: tuple, loss: str = "huber") -> tuple:
    """Find optimal parameters given data using curve_fit.

    Args:
        function (callable): Function to optimize
        parameters (list): Parameter names
        xdata (np.ndarray): Independent variable values
        ydata (np.ndarray): Dependent variable values
        bounds (tuple): ((min_params...), (max_params...))
        loss (str): Loss function for robust estimation

    Returns:
        tuple: (popt_dict, perr_dict)
    """
    popt, pcov = curve_fit(function, xdata=xdata, ydata=ydata, bounds=bounds, loss=loss)
    perr = np.sqrt(np.diag(pcov))
    popt_dict = {k: float(v) for (k, v) in zip(parameters, popt)}
    perr_dict = {f"{k}_stdev": float(v) for (k, v) in zip(parameters, perr)}
    return popt_dict, perr_dict


def LR(L0: np.ndarray, L1: np.ndarray, df: int) -> tuple:
    """Performs likelihood ratio test between null model L0 and alternative model L1.

    Args:
        L0 (np.ndarray): Log-likelihood of null model
        L1 (np.ndarray): Log-likelihood of alternative model
        df (int): Degrees of freedom (difference in number of parameters)

    Returns:
        tuple: (LR_lambda, pvalue)
    """
    LR_lambda = float(-2 * (np.sum(L0) - np.sum(L1)))
    LR_lambda = max(0.0, LR_lambda)
    pval = float(1.0 - chi2.cdf(LR_lambda, df))
    return LR_lambda, pval


def fit_models(
    ref: str,
    model_A,
    model_B,
    damage: np.ndarray,
    mut_count: np.ndarray,
    conserved_count: np.ndarray,
    verbose: bool = False,
) -> dict:
    """Performs model fitting and runs Likelihood Ratio test.

    Args:
        ref (str): Name of reference in alignment file
        model_A (pydamage.models.damage_model): PyDamage H1 Damage model
        model_B (pydamage.models.null_model): PyDamage H0 Null model
        damage (np.ndarray): Amount of damage at each position
        mut_count (np.ndarray): Number of mutations per position
        conserved_count (np.ndarray): Number of conserved bases per position
        verbose (bool): Verbose output

    Returns:
        dict: Fitted parameters and test statistics
    """
    xdata = np.arange(damage.size)

    # 1. Fit Damage Model (H1)
    optim_A, stdev_A = optim(
        function=model_A.fit,
        parameters=model_A.kwds,
        xdata=xdata,
        ydata=damage,
        bounds=model_A.bounds,
    )
    if optim_A["pmax"] < optim_A["pmin"]:
        optim_A["pmax"] = optim_A["pmin"]

    # 2. Fit Null Model (H0)
    optim_B, stdev_B = optim(
        function=model_B.fit,
        parameters=model_B.kwds,
        xdata=xdata,
        ydata=damage,
        bounds=model_B.bounds,
    )

    # 3. Binomial Likelihoods
    c_sites = xdata
    c2c_count_per_site = conserved_count
    binom_k = mut_count
    binom_n = c2c_count_per_site + binom_k

    binom_p_damage = np.clip(model_A.fit(c_sites, **optim_A), 1e-15, 1.0 - 1e-15)
    LA = binom.logpmf(k=binom_k, n=binom_n, p=binom_p_damage)

    binom_p_null = np.clip(model_B.fit(c_sites, **optim_B), 1e-15, 1.0 - 1e-15)
    LB = binom.logpmf(k=binom_k, n=binom_n, p=binom_p_null)

    # 4. Likelihood Ratio Test
    pdiff = len(model_A.kwds) - len(model_B.kwds)
    LR_lambda, pval = LR(L0=LB, L1=LA, df=pdiff)

    res = {}
    res.update(optim_A)
    res.update(stdev_A)
    res.update(optim_B)
    res.update(stdev_B)
    res.update({"pvalue": pval})
    res.update(
        {
            "model_params": list(optim_A.values())
            + list(optim_B.values())
            + list(stdev_A.values())
            + list(stdev_B.values())
        }
    )
    fitted_vals = model_A.fit(x=xdata, **optim_A)
    res["residuals"] = damage - fitted_vals
    res["RMSE"] = RMSE(res["residuals"])
    return res
