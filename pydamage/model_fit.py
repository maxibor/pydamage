#!/usr/bin/env python3

import numpy as np
from pydamage.optim import optim
from scipy.stats import binom
from pydamage.utils import sort_dict_by_keys, RMSE, create_damage_dict
from pydamage.likelihood_ratio import LR


def fit_models(
    ref,
    model_A,
    model_B,
    damage,
    mut_count,
    conserved_count,
    verbose,
):
    """Performs model fitting and runs Likelihood ratio test

    Args:
        ref (str): name of referene in alignment file
        model_A (pydamage.models): Pydamage H1 Damage model
        model_B (pydamage.models): Pydamage H0 Null model
        damage (np.array of float): Amount of damage at each position where CtoT or GtoA (if reversed) transitions were observed
        mut_count(np.array of int): Number of CtoT and GtoA per position
        conserved_count(np.array of int): Number of conserved C and G per position
        verbose (bool): verbose mode
    """
    #################
    # MODEL FITTING #
    #################

    # Only fitting model to interval [0,wlen]
    xdata = np.arange(damage.size)

    res = {}
    optim_A, stdev_A = optim(
        function=model_A.fit,  # damage model
        parameters=model_A.kwds,
        xdata=xdata,
        ydata=damage,
        bounds=model_A.bounds,
    )
    if optim_A["pmax"] < optim_A["pmin"]:  # making sure that fitting makes sense
        optim_A["pmax"] = optim_A["pmin"]

    optim_B, stdev_B = optim(
        function=model_B.fit,  # null model
        parameters=model_B.kwds,
        xdata=xdata,
        ydata=damage,
        bounds=model_B.bounds,
    )

    ##########################
    # LIKELIHOOD CALCULATION #
    ##########################

    # position of sites where C in reference
    c_sites = xdata
    # counts of conserved positions at each site
    c2c_count_per_site = conserved_count
    # counts of CtoT and GtoA mutations at each site
    binom_k = mut_count
    # all C per site
    binom_n = c2c_count_per_site + binom_k

    # Likelihood for model A - Damage Model
    binom_p_damage = model_A.fit(c_sites, **optim_A)
    LA = binom.logpmf(k=binom_k, n=binom_n, p=binom_p_damage)

    # Likelihood for model B - Null Model
    binom_p_null = model_B.fit(c_sites, **optim_B)
    LB = binom.logpmf(k=binom_k, n=binom_n, p=binom_p_null)

    # Difference of number of paramters between model A and model B
    pdiff = len(model_A.kwds) - len(model_B.kwds)

    ################
    # LR TEST #
    ################

    LR_lambda, pval = LR(L0=LB, L1=LA, df=pdiff)

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
    res["residuals"] = damage - model_A.fit(x=xdata, **optim_A)
    res["RMSE"] = RMSE(res["residuals"])
    return res
