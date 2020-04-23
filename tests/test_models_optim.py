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


def test_geom_pmf(generate_data):
    g = models.geom_mod()
    assert g.pmf(x=generate_data[0], geom_p=0.5, geom_pmin=0.01, geom_pmax=0.5).all() == np.array([0.3, 0.3, 0.3, 0.3, 0.3,
                                                                                                   0.3, 0.3, 0.3, 0.3, 0.3,
                                                                                                   0.15471624, 0.15471624, 0.15471624, 0.15471624, 0.15471624,
                                                                                                   0.08207436, 0.08207436, 0.04575342, 0.02759295, 0.01851272,
                                                                                                   0.0139726, 0.01170254, 0.01056751, 0.01]).all()


def test_geom_log_pmf(generate_data):
    g = models.geom_mod()
    assert g.log_pmf(x=generate_data[0], geom_p=0.5, geom_pmin=0.01, geom_pmax=0.5).all() == np.array([-1.2039728, -1.2039728, -1.2039728, -1.2039728, -1.2039728,
                                                                                                       -1.2039728, -1.2039728, -1.2039728, -1.2039728, -1.2039728,
                                                                                                       -1.86616253, -1.86616253, -1.86616253, -1.86616253, -1.86616253,
                                                                                                       -2.50012956, -2.50012956, -3.08448863, -3.59019479, -3.98929721,
                                                                                                       -4.27065681, -4.44794902, -4.54997064, -4.60517019]).all()


def test_geom_optim(generate_data):
    g = models.geom_mod()
    o, e = optim(function=g.pmf,
              parameters=g.kwds,
              xdata=generate_data[1],
              ydata=generate_data[2],
              bounds=g.bounds,
              loss='linear')

    target = {'geom_p': 0.6039535547658853,
              'geom_pmin': 0.03637474931290336,
              'geom_pmax': 0.4211432052501663}
    for k in o:
        assert round(o[k], 3) == round(target[k], 3)


def test_unif_pmf(generate_data):
    u = models.unif_mod()
    assert u.pmf(x=generate_data[0], unif_pmin=0.1).all() == np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                                                                       0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]).all()


def test_unif_log_pmf(generate_data):
    u = models.unif_mod()
    assert u.log_pmf(x=generate_data[0], unif_pmin=0.1).all() == np.array([-2.30258509, -2.30258509, -2.30258509, -2.30258509, -2.30258509,
                                                                           -2.30258509, -2.30258509, -2.30258509, -2.30258509, -2.30258509,
                                                                           -2.30258509, -2.30258509, -2.30258509, -2.30258509, -2.30258509,
                                                                           -2.30258509, -2.30258509, -2.30258509, -2.30258509, -2.30258509,
                                                                           -2.30258509, -2.30258509, -2.30258509, -2.30258509]).all()


def test_unif_optim(generate_data):
    u = models.unif_mod()
    o, e = optim(function=u.pmf,
              parameters=u.kwds,
              xdata=generate_data[1],
              ydata=generate_data[2],
              bounds=u.bounds,
              loss='linear')
    assert o == {'unif_pmin': 0.1}


if __name__ == "__main__":
    data, xdata, ydata = generate_data()
    g = models.geom_mod()
    o = optim(function=g.pmf, parameters=g.kwds, xdata=xdata, ydata=ydata)
