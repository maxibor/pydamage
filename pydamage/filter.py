#!/usr/bin/env python3

import logging
import os
from kneed import KneeLocator
import numpy as np
import pandas as pd

from pydamage.utils import df_to_csv


def define_threshold(pydam_df: pd.DataFrame, min_knee: float = 0.5, alpha: float = 0.05) -> float:
    """Find kneedle point in PyDamage results for the predicted accuracy threshold.

    Args:
        pydam_df (pd.DataFrame): PyDamage results
        min_knee (float): Minimum predicted accuracy threshold
        alpha (float): Maximum q-value threshold

    Returns:
        float: Optimal threshold value
    """
    thresholds = [round(float(i), 2) for i in np.arange(min_knee, 1.0, 0.01)]
    nb_contigs = [
        pydam_df.query(f"predicted_accuracy >= {t} & qvalue <= {alpha}").shape[0]
        for t in thresholds
    ]

    kneedle = KneeLocator(
        thresholds,
        nb_contigs,
        S=1.0,
        curve="convex",
        direction="decreasing",
        online=True,
    )
    return float(kneedle.knee) if kneedle.knee is not None else min_knee


def filter_pydamage_results(pydam_df: pd.DataFrame, acc_thresh: float, alpha: float = 0.05) -> pd.DataFrame:
    """Filter PyDamage results by predicted accuracy and q-value thresholds."""
    return pydam_df.query(f"predicted_accuracy >= {acc_thresh} & qvalue <= {alpha}")


def apply_filter(csv: str, threshold: float, outdir: str, alpha: float = 0.05) -> pd.DataFrame:
    """Apply PyDamage filtering to result CSV file and write output."""
    df = pd.read_csv(csv)
    outfile = "pydamage_filtered_results.csv"

    if threshold == 0 or threshold is None:
        threshold = define_threshold(df, alpha=alpha)
        logging.info(f"Optimal prediction accuracy threshold found to be: {threshold}")

    filt_df = filter_pydamage_results(df, acc_thresh=threshold, alpha=alpha)
    logging.info(f"Filtering PyDamage results with qvalue <= {alpha} and predicted_accuracy >= {threshold}")

    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    df_to_csv(filt_df, outdir, outfile)
    logging.info(f"Filtered PyDamage results written to {outdir}/{outfile}")
    return filt_df
