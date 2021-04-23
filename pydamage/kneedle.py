from kneed import KneeLocator
from numpy import arange
from pydamage.utils import df_to_csv
from pandas import read_csv


def find_knee(pydam_df, min_knee=0.5, alpha=0.05):
    """Find kneedle point in PyDamage results

    Finding the kneedle point to get the optimal
    tradeoff between FP and FN, for the predicted
    accurary threshold

    Args:
        pydam_df (pandas df): pydamage results
        min_knee (float, optional): Min pred_accuracy threshold. Defaults to 0.5
        alpha(float, optional): Alpha q-value threshold
    """
    thresholds = [i.round(2) for i in arange(min_knee, 1, 0.01)]
    nb_contigs = list()
    nb_contigs = []
    for i in thresholds:
        nb_contigs.append(
            pydam_df.query(f"pred_accuracy >= {i} & qvalue <= {alpha}").shape[0]
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
        acc_thresh (float): predictiona accuracy threshold
        alpha (float, optional): Alpha q-value threshold. Defaults to 0.05.
    """

    return pydam_df.query(f"pred_accuracy >= {acc_thresh} & qvalue <= {alpha}")


def apply_filter(csv, outdir, alpha=0.05):
    """Apply pydamage filtering

    Args:
        csv (str): path to pydamage result file
        outdir (str): Path to output directory
        alpha (float, optional): Alpha q-value threshold. Defaults to 0.05.
    """

    df = read_csv(csv)
    outfile = "pydamage_filtered_results.csv"
    knee = find_knee(df)
    print(f"Optimal prediction accuracy threshold found to be: {knee}")
    filt_df = filter_pydamage_results(df, acc_thresh=knee)
    print(
        f"Filtering PyDamage results with qvalue <={alpha} and pred_accuracy >= {knee}"
    )
    df_to_csv(filt_df, outdir, outfile)
    print(f"Filtered PyDamage results written to {outdir}/{outfile}")
    return filt_df
