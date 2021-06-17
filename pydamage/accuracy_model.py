import pkg_resources
from pypmml import Model


def load_model():
    """Returns the pmml model"""
    model_path = pkg_resources.resource_stream(
        __name__, "models/pydamage_glm_model.pmml"
    )
    model = Model.load(model_path)
    return model


def prepare_data(pd_df):
    """Prepare pydamage result data for accuracy modelling

    Args:
        pd_df (pandas DataFrame):pydamage df result
    """
    pd_df = pd_df[["coverage", "reflen", "damage_model_pmax"]].rename(
        columns={
            "damage_model_pmax": "damage",
            "reflen": "contiglength",
            "coverage": "actual_cov",
        }
    )
    pd_df = pd_df.astype({"actual_cov": float, "contiglength": int, "damage": float})
    return pd_df


def fit_model(df, model):
    """Fit GLM model to data

    Args:
        df (pandas DataFrame): prepared pydamage results
        model (pypmml model): GLM accuracy model
    """
    prediction = list(model.predict(df)["Predicted_sig"])
    df["predicted_accuracy"] = prediction
    return df["predicted_accuracy"].to_frame()
