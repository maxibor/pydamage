#!/usr/bin/env python3


from scipy.optimize import curve_fit


def optim(function, parameters, xdata, ydata):
    """Find optimal parameters given data

    Args:
        function (function): function to optimize
        xdata (np array): x values
        ydata (np array): y values
    Returns:
        (dict): 'parameter_name':'parameter_value'
    """
    op_par = curve_fit(function, xdata=xdata, ydata=ydata)
    op_par_dict = {k: v for (k, v) in zip(parameters, op_par)}
    return(op_par_dict)
