#!/usr/bin/env python3

import numpy as np
from scipy.stats import norm
from pydamage.optim import optim
import collections
from pydamage.utils import sort_dict_by_keys

def vuong_closeness(ref, model_A, model_B, ct_data, ga_data, all_bases, wlen, verbose):
    """Performs model fitting and runs Vuong's closeness test

    Args:
        ref (str): name of referene in alignment file
        model_A (pydamage.models): Pydamage H1 model
        model_B (pydamage.models): Pydamage H0 model
        ct_data (list of int): List of positions where CtoT transitions were observerd
        ga_data (list of int): List of positions where GtoA transitions were observerd
        all_bases (list of int): List of positions where a base is aligned
        wlen (int): window length
        verbose (bool): verbose mode
    """
    all_bases_pos, all_bases_counts = np.unique(np.sort(all_bases), return_counts=True)
    c2t_pos, c2t_counts = np.unique(np.sort(ct_data), return_counts=True)
    g2a_pos, g2a_counts = np.unique(np.sort(ga_data), return_counts=True)
    c2t = dict(zip(c2t_pos, c2t_counts))
    g2a = dict(zip(g2a_pos, g2a_counts))

    # Adding zeros at positions where no damage is observed
    for i in all_bases_pos:
        if all_bases_counts[i] > 0:
            if i not in c2t:
                c2t[i] = 0
            if i not in g2a:
                g2a[i] = 0
    c2t = sort_dict_by_keys(c2t)
    g2a = sort_dict_by_keys(g2a)

    xdata = np.array(list(c2t.keys()))
    counts = np.array(list(c2t.values()))

    ydata = list(counts/counts.sum())
    qlen = len(ydata)
    ydata_counts = {i:c for i,c in enumerate(ydata)}
    ctot_out = {f"CtoT-{k}":v for k,v in enumerate(ydata)}
    
    g2a_counts = np.array(list(g2a.values()))
    y_ga = list(g2a_counts/g2a_counts.sum())
    y_g2a_counts = {i:c for i,c in enumerate(y_ga)}
    gtoa_out = {f"GtoA-{k}":v for k,v in enumerate(y_ga)}

   
    # for i in all_bases_pos:
    #     print(f"{ref},{i},{all_bases_counts[i]}")
    for i in range(qlen):
        if i not in ydata_counts :
            ydata_counts[i] = np.nan
        if f"CtoT-{i}" not in ctot_out:
            ctot_out[f"CtoT-{i}"] = np.nan
        if f"GtoA-{i}" not in gtoa_out:
            gtoa_out[f"GtoA-{i}"] = np.nan


    xdata = xdata[:wlen]
    ydata = ydata[:wlen]
    ct_data = np.array(ct_data)
    ct_data = ct_data[ct_data<wlen]

    res = {}
    optim_A, stdev_A = optim(function=model_A.pmf,
                    parameters=model_A.kwds,
                    xdata=xdata,
                    ydata=ydata,
                    bounds=model_A.bounds)
    optim_B, stdev_B = optim(function=model_B.pmf,
                    parameters=model_B.kwds,
                    xdata=xdata,
                    ydata=ydata,
                    bounds=model_B.bounds)

    LA = model_A.log_pmf(x=ct_data, **optim_A)
    LB = model_B.log_pmf(x=ct_data, **optim_B)
    pdiff = len(model_A.kwds) - len(model_B.kwds)
    LR = LA.sum() - LB.sum() - (pdiff/2)*np.log(len(ct_data))
    omega = np.std(LA-LB)
    Z = LR/(np.sqrt(len(ct_data))*omega)
    pval = norm.cdf(Z)
    # if verbose:
    #     print(f"\nReference: {ref}")
    #     print(f"Vuong closeness test Z-score for {ref}: {round(Z, 4)}")
    #     print(f"Vuong closeness test p-value for {ref}: {round(pval, 4)}")
    #     print(f"Model A parameters for {ref}: {optim_A}")
    #     print(f"Model B parameters for {ref}: {optim_B}")
    res.update(ydata_counts)
    res.update(ctot_out)
    res.update(gtoa_out)
    res.update(optim_A)
    res.update(stdev_A)
    res.update(optim_B)
    res.update(stdev_B)
    res.update({'pvalue': pval})
    res.update({'base_cov': all_bases_counts})
    res.update({'model_params': list(optim_A.values()) + list(optim_B.values())+list(stdev_A.values())+list(stdev_B.values())})
    res['qlen'] = qlen
    return(res)
