#!/usr/bin/env python3

import numpy as np


class geom_mod():
    def __init__(self):
        self.kwds = ['geom_p', 'geom_pmin', 'geom_pmax']
        self.bounds = ((1e-15, 1e-15, 1e-15), (0.99, 0.2, 0.99))

    def __repr__(self):
        return(
            f"""
            A modified geometric function of formula:
            y = ((((1-geom_p)**x)*geom_p) - xmin)/(xmax - xmin))*(geom_pmax - geom_pmin) + geom_pmin
            With parameters:
            - geom_p
            - geom_pmin
            - geom_pmax
            - xmin
            - xmax
            """
        )

    def _geom_pmf(self, x, geom_p):
        """Probability mass function of the discrete geometric distribution

        Args:
            x (int): position
            geom_p (float): parameter of distribition
        Returns:
            float: PMF(x)
        """
        return(((1-geom_p)**x)*geom_p)

    def pmf(self, x, geom_p, geom_pmin, geom_pmax, wlen=35):
        """Probability mass function of the damage model 

        Args:
            x (numpy array) data
            geom_p (float): shape parameter
            geom_pmin (float): min y value
            geom_pmax (float): max y value
            wlen(int):  window length
        Returns:
            np.array: PMF(x)
        """
        vec_base_geom = np.vectorize(self._geom_pmf)

        base_geom = vec_base_geom(x, geom_p)

        xmax = self._geom_pmf(0, geom_p)
        xmin = self._geom_pmf(wlen-1, geom_p)
        scaled_geom = ((base_geom - xmin)/(xmax - xmin)) * \
            (geom_pmax - geom_pmin) + geom_pmin
        return(scaled_geom)

    def log_pmf(self, x, geom_p, geom_pmin, geom_pmax, wlen):
        """Log Probability mass function of the discrete geometric function

        Args:
            x (numpy array) data
            geom_p (float): shape parameter
            geom_pmin (float): min y value
            geom_pmax (float): max y value
            wlen(int):  window length
        Returns:
            np.array: LogPMF(x)
        """
        return(np.log(self.pmf(x, geom_p, geom_pmin, geom_pmax, wlen)))


class unif_mod():
    def __init__(self):
        self.kwds = ('unif_pmin',)
        self.bounds = ((1e-15,), (0.2,))

    def __repr__(self):
        return(
            f"""
            A modified uniform function of formula
            y = unif_pmin
            With parameters:
            - unif_pmin
            """
        )

    def pmf(self, x, unif_pmin):
        """Probability mass function of the discrete geometric function

        Args:
            x (numpy array) data
            unif_pmin (float): target y value
        Returns:
            np.array: PMF(x)
        """
        return(np.array([unif_pmin]*len(x)))

    def log_pmf(self, x, unif_pmin):
        """Log probability mass function of the discrete geometric function

        Args:
            x (numpy array) data
            unif_pmin (float): target y value
        Returns:
            np.array: LogPMF(x)
        """
        return(np.log(self.pmf(x, unif_pmin)))
