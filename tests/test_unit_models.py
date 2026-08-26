import pytest
import numpy as np
import pandas as pd
from pydamage import models
from pydamage.model_fit import optim, LR
from pydamage.models import prepare_data, glm_predict
from pydamage.utils import RMSE


def test_damage_model_properties():
    dmg = models.damage_model()
    assert "p" in dmg.kwds
    assert "pmin" in dmg.kwds
    assert "pmax" in dmg.kwds
    assert len(dmg.bounds[0]) == 3
    assert len(dmg.bounds[1]) == 3


def test_damage_model_fit_array():
    dmg = models.damage_model()
    x = np.arange(10)
    y = dmg.fit(x=x, p=0.5, pmin=0.01, pmax=0.5, wlen=10)
    assert len(y) == 10
    assert y[0] == pytest.approx(0.5, abs=1e-4)
    assert y[-1] == pytest.approx(0.01, abs=1e-4)
    assert np.all(np.diff(y) <= 0)  # monotonically decreasing


def test_null_model_properties():
    null = models.null_model()
    assert "p0" in null.kwds
    assert len(null.bounds[0]) == 1
    assert len(null.bounds[1]) == 1


def test_null_model_fit():
    null = models.null_model()
    x = np.arange(15)
    y = null.fit(x=x, p0=0.05)
    assert len(y) == 15
    assert np.all(y == 0.05)


def test_likelihood_ratio_test():
    # When L1 is much higher than L0 (damage model explains data much better)
    L0 = np.array([-10.0, -10.0, -10.0])
    L1 = np.array([-2.0, -2.0, -2.0])
    lr_lambda, pval = LR(L0=L0, L1=L1, df=2)
    assert lr_lambda > 0
    assert pval < 0.001

    # When L1 equals L0
    lr_lambda_zero, pval_one = LR(L0=L0, L1=L0, df=2)
    assert lr_lambda_zero == 0.0
    assert pval_one == pytest.approx(1.0)


def test_rmse():
    residuals = np.array([1.0, -1.0, 1.0, -1.0])
    assert RMSE(residuals) == pytest.approx(1.0)
    assert RMSE(np.zeros(5)) == 0.0


def test_glm_accuracy_prediction():
    df = pd.DataFrame(
        {
            "coverage": [0.5, 5.0, 50.0],
            "reflen": [1000, 10000, 100000],
            "damage_model_pmax": [0.01, 0.15, 0.30],
        }
    )
    prep = prepare_data(df)
    assert "actual_cov" in prep.columns
    assert "contiglength" in prep.columns
    assert "damage" in prep.columns
    
    pred = glm_predict(prep, models.glm_model_params)
    assert "predicted_accuracy" in pred.columns
    # Probabilities must be strictly bounded in (0, 1)
    assert (pred["predicted_accuracy"] >= 0.0).all()
    assert (pred["predicted_accuracy"] <= 1.0).all()
    # Higher damage and coverage should result in higher predicted accuracy
    assert pred["predicted_accuracy"].iloc[2] > pred["predicted_accuracy"].iloc[0]
