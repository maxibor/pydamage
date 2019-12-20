#!/usr/bin/env python3

import numpy as np
from scipy.stats import norm
from optim import optim


def vuong_closeness(model_A, model_B, data):
    xdata, counts = np.unique(np.sort(data), return_counts=True)
    ydata = counts/counts.sum()
    res = {}
    optim_A = optim(function=model_A.pmf,
                    parameters=model_A.kwds,
                    xdata=xdata,
                    ydata=ydata,
                    bounds=model_A.bounds)
    optim_B = optim(function=model_B.pmf,
                    parameters=model_B.kwds,
                    xdata=xdata,
                    ydata=ydata,
                    bounds=model_B.bounds)

    LA = model_A.log_pmf(x=data, **optim_A)
    LB = model_B.log_pmf(x=data, **optim_B)
    pdiff = len(model_A.kwds) - len(model_B.kwds)
    LR = LA.sum() - LB.sum() - (pdiff/2)*np.log(len(data))
    omega = np.std(LA-LB)
    Z = LR/(np.sqrt(len(data))*omega)
    pval = norm.cdf(Z)
    print(f"Vuong closeness test Z-score: {round(Z, 4)}")
    print(f"Vuong closeness test p-value: {round(pval, 4)}")
    print(f"Model A parameters: {optim_A}")
    print(f"Model B parameters: {optim_B}")
    res.update(optim_A)
    res.update(optim_B)
    res.update({'pvalue': pval})
    return(res)
