from kneed import KneeLocator
from numpy import arange
from pydamage.utils import df_to_csv
from pandas import read_csv
import os


def define_threshold(pydam_df, min_knee=0.5, alpha=0.05):
    """Find kneedle point in PyDamage results

    Finding the kneedle point to get the optimal
    tradeoff between FP and FN, for the predicted
    accurary threshold

    Args:
        pydam_df (pandas df): pydamage results
        min_knee (float, optional): Min predicted_accuracy threshold.
        alpha(float, optional): Alpha q-value threshold
    """
    thresholds = [i.round(2) for i in arange(min_knee, 1, 0.01)]
    nb_contigs = list()
    nb_contigs = []
    for i in thresholds:
        nb_contigs.append(
            pydam_df.query(f"predicted_accuracy >= {i} & qvalue <= {alpha}").shape[0]
        )
    kneedle = KneeLocator(
        thresholds,
        nb_contigs,
        S=1.0,
        curve="convex",
        direction="decreasing",
        online=True,
    )
    print(thresholds)
    print(nb_contigs)
    return kneedle.knee


def filter_pydamage_results(pydam_df, acc_thresh, alpha=0.05):
    """Filter pydamage results on pred_accuracy and qvalue

    Args:
        pydam_df (pandas df): pydamage results
        acc_thresh (float): predicted accuracy threshold
        alpha (float, optional): Alpha q-value threshold. Defaults to 0.05.
    """

    return pydam_df.query(f"predicted_accuracy >= {acc_thresh} & qvalue <= {alpha}")


def apply_filter(csv, threshold, outdir, alpha=0.05):
    """Apply pydamage filtering

    Args:
        csv (str): path to pydamage result file
        threshold(float): Treshold value. 0 is for finding threshold with kneed method
        outdir (str): Path to output directory
        alpha (float, optional): Alpha q-value threshold. Defaults to 0.05.
    """

    df = read_csv(csv)
    outfile = "pydamage_filtered_results.csv"
    if threshold == 0:
        threshold = define_threshold(df)
        print(f"Optimal prediction accuracy threshold found to be: {threshold}")
    filt_df = filter_pydamage_results(df, acc_thresh=threshold)
    print(
        f"Filtering PyDamage results with qvalue <= {alpha} and predicted_accuracy >= {threshold}"
    )
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    df_to_csv(filt_df, outdir, outfile)
    print(f"Filtered PyDamage results written to {outdir}/{outfile}")
    return filt_df
