#!/usr/bin/env python3

import pytest
import numpy as np

from pydamage import models
from pydamage.optim import optim


@pytest.fixture(autouse=True)
def generate_data():
    x = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
                  1, 1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9])
    xdata, counts = np.unique(x, return_counts=True)
    ydata = counts/counts.sum()
    return(x, xdata, ydata)


def test_damage_model_fit(generate_data):
    g = models.damage_model()
    assert g.fit(x=generate_data[0], p=0.5, pmin=0.01, pmax=0.5).all() == np.array([0.3, 0.3, 0.3, 0.3, 0.3,
                                                                                    0.3, 0.3, 0.3, 0.3, 0.3,
                                                                                    0.15471624, 0.15471624, 0.15471624, 0.15471624, 0.15471624,
                                                                                    0.08207436, 0.08207436, 0.04575342, 0.02759295, 0.01851272,
                                                                                    0.0139726, 0.01170254, 0.01056751, 0.01]).all()


def test_damage_model_optim(generate_data):
    g = models.damage_model()
    o, e = optim(function=g.fit,
                 parameters=g.kwds,
                 xdata=generate_data[1],
                 ydata=generate_data[2],
                 bounds=g.bounds,
                 loss='linear')

    target = {'p': 0.6039535547658853,
              'pmin': 0.03637474931290336,
              'pmax': 0.4211432052501663}
    for k in o:
        assert round(o[k], 3) == round(target[k], 3)


def test_null_model_fit(generate_data):
    u = models.null_model()
    assert u.fit(x=generate_data[0], p0=0.1).all() == np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                                                                0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]).all()


def test_null_model_optim(generate_data):
    u = models.null_model()
    o, e = optim(function=u.fit,
                 parameters=u.kwds,
                 xdata=generate_data[1],
                 ydata=generate_data[2],
                 bounds=u.bounds,
                 loss='linear')
    assert o == {'p0': 0.1000000000000005}


if __name__ == "__main__":
    data, xdata, ydata = generate_data()
    g = models.damage_model()
    o = optim(function=g.fit, parameters=g.kwds, xdata=xdata, ydata=ydata)
