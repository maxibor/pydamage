import pkg_resources
import pandas as pd
import numpy as np
import pickle


def load_model():
    """Returns the pmml model"""
    # This is a stream,like object. If you want the actual info, call
    # stream.read()
    stream = pkg_resources.resource_stream(__name__, "models/glm_accuracy_model.pickle")
    return pickle.load(stream)


def prepare_data(pd_df):
    """Prepare pydamage result data for accuracy modelling

    Args:
        pd_df (pandas DataFrame):pydamage df result
    """
    coverage_bins = pd.IntervalIndex.from_tuples(
        [
            (0, 2),
            (2, 3),
            (3, 5),
            (5, 10),
            (10, 20),
            (20, 50),
            (50, 100),
            (100, 200),
            (200, np.inf),
        ]
    )
    coverage_bins_labels = [
        "1-2",
        "2-3",
        "3-5",
        "5-10",
        "10-20",
        "20-50",
        "50-100",
        "100-200",
        "200-500",
    ]

    reflen_bins = pd.IntervalIndex.from_tuples(
        [
            (0, 1000),
            (1000, 2000),
            (2000, 5000),
            (5000, 10000),
            (10000, 20000),
            (20000, 50000),
            (50000, 100000),
            (100000, 200000),
            (200000, np.inf),
        ]
    )

    reflen_bins_labels = [
        "500-1000",
        "1000-2000",
        "2000-5000",
        "5000-10000",
        "10000-20000",
        "20000-50000",
        "50000-100000",
        "100000-200000",
        "200000-500000",
    ]
    simu_cov = pd.cut(pd_df["coverage"], coverage_bins)
    simu_cov.cat.rename_categories(coverage_bins_labels, inplace=True)
    pd_df["simuCov"] = simu_cov

    simu_contig_length = pd.cut(pd_df["reflen"], reflen_bins)
    simu_contig_length.cat.rename_categories(reflen_bins_labels, inplace=True)

    pd_df["simuContigLength"] = simu_contig_length
    pd_df = pd_df[
        ["simuCov", "simuContigLength", "damage_model_pmax", "gc_content"]
    ].rename(columns={"damage_model_pmax": "damage", "gc_content": "GCcontent"})

    return pd_df


def fit_model(df, model):
    """Fit GLM model to data

    Args:
        df (pandas DataFrame): prepared pydamage results
        model (pypmml model): GLM accuracy model
    """
    return model.predict(df).to_frame(name="pred_accuracy")
