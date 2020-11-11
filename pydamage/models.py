#!/usr/bin/env python3

import numpy as np


class damage_model():
    def __init__(self):
        self.kwds = ['p', 'pmin', 'pmax']
        self.bounds = ((1e-15, 1e-15, 1e-15), (0.99, 0.2, 0.99))

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

    def _geom_pmf(self, x, p):
        """Probability mass function of the discrete geometric distribution

        Args:
            x (int): position
            p (float): parameter of distribition
        Returns:
            float: PMF(x)
        """
        return(((1-p)**x)*p)

    def fit(self, x, p, pmin, pmax, wlen=35):
        """Damage model function 

        Args:
            x (numpy ndarray or int) data
            p (float): shape parameter
            pmin (float): min y value
            pmax (float): max y value
            wlen(int):  window length
        Returns:
            np.array: PMF(x)
        """
        vec_base_geom = np.vectorize(self._geom_pmf)

        if type(x) == int:
            x = np.array(x)

        base_geom = vec_base_geom(x, p)

        xmax = self._geom_pmf(0, p)
        xmin = self._geom_pmf(wlen-1, p)
        scaled_geom = ((base_geom - xmin)/(xmax - xmin)) * \
            (pmax - pmin) + pmin
        return(scaled_geom)


class null_model():
    def __init__(self):
        self.kwds = ('p0',)
        self.bounds = ((1e-15,), (0.2,))

    def __repr__(self):
        return(
            f"""
            A modified uniform function of formula
            y = p0
            With parameters:
            - p0
            """
        )

    def fit(self, x, p0):
        """Null model function

        Args:
            x (numpy array) data
            p0 (float): target y value
        Returns:
            np.array: PMF(x)
        """
        return(np.array([p0]*len(x)))
