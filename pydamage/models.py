#!/usr/bin/env python3

import numpy as np


class geom_mod():
    def __init__(self):
        self.kwds = ('p', 'pmin', 'pmax')

    def __repr__(self):
        return(
            f"""
            A modified geometric function of formula:
            y = ((((1-p)**x)*p) - xmin)/(xmax - xmin))*(pmax - pmin) + pmin
            With parameters:
            - p
            - pmin
            - pmax
            - xmin
            - xmax
            """
        )

    def pmf(self, x, p, pmin, pmax):
        """Probability mass function of the discrete geometric function

        Args:
            x (numpy array) data
            p (float): shape parameter
            pmin (float): min y value
            pmax (float): max y value
        Returns:
            (numpy array): PMF(x)
        """
        base_geom = ((1-p)**x)*p
        xmin = min(base_geom)
        xmax = max(base_geom)
        scaled_geom = ((base_geom - xmin)/(xmax - xmin)) * (pmax - pmin) + pmin
        return(scaled_geom)

    def log_pmf(self, x, p, pmin, pmax):
        """Log Probability mass function of the discrete geometric function

        Args:
            x (numpy array) data
            p (float): shape parameter
            pmin (float): min y value
            pmax (float): max y value
        Returns:
            (numpy array): LogPMF(x)
        """
        return(np.log(self.pmf(x=x, p=p, pmin=pmin, pmax=pmax)))


class unif_mod():
    def __init__(self):
        self.kwds = ('pmin')

    def __repr__(self):
        return(
            f"""
            A modified uniform function of formula
            y = pmin
            With parameters:
            - pmin
            """
        )

    def pmf(self, x, pmin):
        """Probability mass function of the discrete geometric function

        Args:
            x (numpy array) data
            pmin (float): target y value
        Returns:
            (numpy array): PMF(x)
        """
        return(np.array([pmin]*len(x)))

    def log_pmf(self, x, pmin):
        """Log probability mass function of the discrete geometric function

        Args:
            x (numpy array) data
            pmin (float): target y value
        Returns:
            (numpy array): LogPMF(x)
        """
        return(np.log(self.pmf(x=x, pmin=pmin)))


def vuong_closeness(model_A, model_B, data):
    xdata, counts = np.unique(np.sort(data), return_counts=True)
    ydata = counts/counts.sum()

    opt_par_A = model_A.optim(xdata=xdata, ydata=ydata)
    opt_par_B = model_B.optim(xdata=xdata, ydata=ydata)

    log_pmf_A = model_A.log_pmf(x=data, **opt_par_A)
    log_pmf_B = model_B.log_pmf(x=data, **opt_par_B)
