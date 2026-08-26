#!/usr/bin/env python3

import numpy as np
import pandas as pd


class damage_model:
    """Modified geometric damage model for ancient DNA nucleotide substitutions."""

    def __init__(self):
        self.kwds = ["p", "pmin", "pmax"]
        self.bounds = ((1e-15, 1e-15, 1e-15), (0.99, 0.2, 0.99))

    def __repr__(self):
        return """
            A modified geometric function of formula:
            y = ((((1-p)**x)*p) - xmin)/(xmax - xmin))*(pmax - pmin) + pmin
            With parameters:
            - p
            - pmin
            - pmax
            - xmin
            - xmax
            """

    def _geom_pmf(self, x, p):
        """Probability mass function of the discrete geometric distribution."""
        return ((1.0 - p) ** x) * p

    def fit(self, x, p, pmin, pmax, wlen=35, **kwargs):
        """Damage model function evaluated at positions x.

        Args:
            x (numpy.ndarray or int or float): Data positions
            p (float): Shape parameter
            pmin (float): Minimum y value
            pmax (float): Maximum y value
            wlen (int): Window length

        Returns:
            numpy.ndarray: Model values at x
        """
        x_arr = np.asarray(x, dtype=np.float64)
        base_geom = ((1.0 - p) ** x_arr) * p

        xmax = p
        xmin = ((1.0 - p) ** (wlen - 1)) * p
        denom = xmax - xmin

        if denom == 0:
            scaled_geom = np.full_like(base_geom, pmin)
        else:
            scaled_geom = ((base_geom - xmin) / denom) * (pmax - pmin) + pmin

        if x_arr.ndim > 0 and len(scaled_geom) > wlen:
            scaled_geom[wlen:] = pmin

        return scaled_geom


class null_model:
    """Null (uniform) background model for nucleotide substitutions."""

    def __init__(self):
        self.kwds = ("p0",)
        self.bounds = ((1e-15,), (0.2,))

    def __repr__(self):
        return """
            A modified uniform function of formula
            y = p0
            With parameters:
            - p0
            """

    def fit(self, x, p0, **kwargs):
        """Null model function.

        Args:
            x (numpy.ndarray or sequence): Data positions
            p0 (float): Target baseline value

        Returns:
            numpy.ndarray: Uniform array of p0 values
        """
        x_arr = np.asarray(x)
        if x_arr.ndim == 0:
            return np.array(p0, dtype=np.float64)
        return np.full(len(x_arr), p0, dtype=np.float64)


# Parameters for the predicted accuracy GLM model with a logit link function
glm_model_params = {
    "intercept": -1.78842427736877,
    "actual_cov": 0.0271448079667983,
    "damage": 20.4193418075295,
    "contiglength": 4.19608849503652e-05,
}


def prepare_data(pd_df: pd.DataFrame) -> pd.DataFrame:
    """Prepare PyDamage result data for accuracy modeling."""
    df = pd_df[["coverage", "reflen", "damage_model_pmax"]].rename(
        columns={
            "damage_model_pmax": "damage",
            "reflen": "contiglength",
            "coverage": "actual_cov",
        }
    )
    return df.astype({"actual_cov": float, "contiglength": int, "damage": float})


def glm_predict(df: pd.DataFrame, model_params: dict) -> pd.DataFrame:
    """Predict accuracy using logistic regression model."""
    linear_pred = (
        model_params["intercept"]
        + model_params["actual_cov"] * df["actual_cov"]
        + model_params["damage"] * df["damage"]
        + model_params["contiglength"] * df["contiglength"]
    )
    pred_accuracy = 1.0 / (1.0 + np.exp(-linear_pred))
    return pred_accuracy.to_frame(name="predicted_accuracy")
