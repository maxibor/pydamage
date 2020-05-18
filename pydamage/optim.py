#!/usr/bin/env python3


from scipy.optimize import curve_fit
from numpy import sqrt, diag


def optim(function, parameters, xdata, ydata, bounds, loss='huber'):
    """Find optimal parameters given data

    Args:
        function (function): function to optimize
        xdata (np array): x values
        ydata (np array): y values
        bounds (tuple of tuple): optimization bounds
                    ((par1_min, par2_min), (par1_max, par2_max))
        loss (str): loss function. See https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html#scipy.optimize.least_squares
    Returns:
        dict: popt_dict - 'parameter_name':'parameter_value'
        dict: perr_dict - 'parameter_name':'standard_deviation'
    """
    popt, pcov = curve_fit(function, xdata=xdata,
                           ydata=ydata, bounds=bounds, loss=loss)
    perr = sqrt(diag(pcov))
    popt_dict = {k: v for (k, v) in zip(parameters, popt)}
    perr_dict = {f"{k}_stdev": v for (k, v) in zip(parameters, perr)}
    return(popt_dict, perr_dict)
