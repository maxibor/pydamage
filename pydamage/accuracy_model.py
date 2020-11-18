import pkg_resources
import pandas as pd
import numpy as np
import pickle
import gzip


def load_model():
    """Returns the pmml model"""
    model_path = pkg_resources.resource_stream(
        __name__, "models/accuracy_model_v2_python.pickle.gz"
    )
    with gzip.open(model_path, 'rb') as mod:
        return pickle.load(mod)


def prepare_data(pd_df):
    """Prepare pydamage result data for accuracy modelling

    Args:
        pd_df (pandas DataFrame):pydamage df result
    """
    pd_df = pd_df[["coverage", "reflen", "damage_model_pmax"]].rename(
        columns={"damage_model_pmax": "damage", "reflen": "contiglength"}
    )
    pd_df = pd_df.astype({"coverage": float, "contiglength": int, "damage": float})
    return pd_df


def fit_model(df, model):
    """Fit GLM model to data

    Args:
        df (pandas DataFrame): prepared pydamage results
        model (pypmml model): GLM accuracy model
    """
    return model.predict(df).to_frame(name="pred_accuracy")
