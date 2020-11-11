#!/usr/bin/env python3

import numpy as np
from pydamage.optim import optim
from scipy.stats import binom
import collections
from pydamage.utils import sort_dict_by_keys, RMSE, create_ct_cc_dict
from pydamage.likelihood_ratio import LR


def fit_models(
    ref, model_A, model_B, ct_data, cc_data, ga_data, all_bases, wlen, verbose
):
    """Performs model fitting and runs Likelihood ratio test

    Args:
        ref (str): name of referene in alignment file
        model_A (pydamage.models): Pydamage H1 Damage model
        model_B (pydamage.models): Pydamage H0 Null model
        ct_data (list of int): List of positions where CtoT transitions were observed
        ga_data (list of int): List of positions where GtoA transitions were observed
        cc_data (list of int): List of positions where C in ref and query
        all_bases (list of int): List of positions where a base is aligned
        wlen (int): window length
        verbose (bool): verbose mode
    """
    all_bases_pos, all_bases_counts = np.unique(np.sort(all_bases), return_counts=True)
    c2t_pos, c2t_counts = np.unique(np.sort(ct_data), return_counts=True)
    g2a_pos, g2a_counts = np.unique(np.sort(ga_data), return_counts=True)
    c2t = dict(zip(c2t_pos, c2t_counts))
    g2a = dict(zip(g2a_pos, g2a_counts))

    # Getting Positions and counts of C2C and C2T
    c2t_dict, c2c_dict = create_ct_cc_dict(ct_data=ct_data, cc_data=cc_data, wlen=wlen)

    # Adding zeros at positions where no damage is observed
    for i in all_bases_pos:
        if i < len(all_bases_counts):
            if all_bases_counts[i] > 0:
                if i not in c2t:
                    c2t[i] = 0
                if i not in g2a:
                    g2a[i] = 0
    c2t = sort_dict_by_keys(c2t)
    g2a = sort_dict_by_keys(g2a)

    xdata = np.array(list(c2t.keys()))
    counts = np.array(list(c2t.values()))

    ydata = list(counts / counts.sum())
    qlen = len(ydata)
    ydata_counts = {i: c for i, c in enumerate(ydata)}
    ctot_out = {f"CtoT-{k}": v for k, v in enumerate(ydata)}

    g2a_counts = np.array(list(g2a.values()))
    y_ga = list(g2a_counts / g2a_counts.sum())
    gtoa_out = {f"GtoA-{k}": v for k, v in enumerate(y_ga)}

    for i in range(qlen):
        if i not in ydata_counts:
            ydata_counts[i] = np.nan
        if f"CtoT-{i}" not in ctot_out:
            ctot_out[f"CtoT-{i}"] = np.nan
        if f"GtoA-{i}" not in gtoa_out:
            gtoa_out[f"GtoA-{i}"] = np.nan

    #################
    # MODEL FITTING #
    #################

    # Only fitting model to interval [0,wlen]
    xdata = xdata[:wlen]
    ydata = ydata[:wlen]

    res = {}
    optim_A, stdev_A = optim(
        function=model_A.fit,  # damage model
        parameters=model_A.kwds,
        xdata=xdata,
        ydata=ydata,
        bounds=model_A.bounds,
    )
    if optim_A["pmax"] < optim_A["pmin"]:  # making sure that fitting makes sense
        optim_A["pmax"] = optim_A["pmin"]

    optim_B, stdev_B = optim(
        function=model_B.fit,  # null model
        parameters=model_B.kwds,
        xdata=xdata,
        ydata=ydata,
        bounds=model_B.bounds,
    )

    ##########################
    # LIKELIHOOD CALCULATION #
    ##########################

    # position of sites where C in reference
    c_sites = np.array(list(c2t_dict.keys()))
    # counts of C2C at each site
    c2c_count_per_site = np.array(list(c2c_dict.values()))
    # counts of C2T at each site
    binom_k = np.array(list(c2t_dict.values()))
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

    res.update(ydata_counts)
    res.update(ctot_out)
    res.update(gtoa_out)
    res.update(optim_A)
    res.update(stdev_A)
    res.update(optim_B)
    res.update(stdev_B)
    res.update({"pvalue": pval})
    res.update({"base_cov": all_bases_counts})
    res.update(
        {
            "model_params": list(optim_A.values())
            + list(optim_B.values())
            + list(stdev_A.values())
            + list(stdev_B.values())
        }
    )
    res["qlen"] = qlen
    res["residuals"] = ydata - model_A.fit(x=xdata, **optim_A)
    res["RMSE"] = RMSE(res["residuals"])
    res["wlen"] = wlen
    return res
