from numpy import exp


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


def glm_predict(df, model_params):
    """Predict accuracy

    Args:
        df (pandas DataFrame): prepared pydamage results
        model_params (dict): model parameters
    """
    df["predicted_accuracy"] = (
        model_params["intercept"]
        + model_params["actual_cov"] * df["actual_cov"]
        + model_params["damage"] * df["damage"]
        + model_params["contiglength"] * df["contiglength"]
    )
    df["predicted_accuracy"] = 1 / (1 + exp(-df["predicted_accuracy"]))
    return df["predicted_accuracy"].to_frame()
