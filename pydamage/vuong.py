#!/usr/bin/env python3

import numpy as np
from scipy.stats import norm
from pydamage.optim import optim


def vuong_closeness(ref, model_A, model_B, ct_data, ga_data, wlen, verbose):
    xdata, counts = np.unique(np.sort(ct_data), return_counts=True)
    ydata = list(counts/counts.sum())
    qlen = len(ydata)
    ydata_counts = {i:c for i,c in enumerate(ydata)}
    # print(ydata_counts)
    ctot_out = {f"CtoT-{k}":v for k,v in enumerate(ydata)}
    
    ga_pos, ga_counts = np.unique(np.sort(ga_data), return_counts=True)
    y_ga = list(ga_counts/ga_counts.sum())
    y_ga_counts = {i:c for i,c in enumerate(y_ga)}
    gtoa_out = {f"GtoA-{k}":v for k,v in enumerate(y_ga)}

    for i in range(qlen):
        if i not in ydata_counts:
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
    if verbose:
        print(f"\nReference: {ref}")
        print(f"Vuong closeness test Z-score for {ref}: {round(Z, 4)}")
        print(f"Vuong closeness test p-value for {ref}: {round(pval, 4)}")
        print(f"Model A parameters for {ref}: {optim_A}")
        print(f"Model B parameters for {ref}: {optim_B}")
    res.update(ydata_counts)
    res.update(ctot_out)
    res.update(gtoa_out)
    res.update(optim_A)
    res.update(stdev_A)
    res.update(optim_B)
    res.update(stdev_B)
    res.update({'pvalue': pval})
    res['qlen'] = qlen
    return(res)
