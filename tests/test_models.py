#!/usr/bin/env python3

import pytest
import numpy as np

from pydamage import models
from pydamage.optim import optim


@pytest.fixture(autouse=True)
def generate_data():
    x = np.array([3, 1, 3, 1, 1, 1, 1, 2, 1, 3])
    xdata, counts = np.unique(x, return_counts=True)
    ydata = counts/counts.sum()
    return(x, xdata, ydata)


def test_geom_pmf(generate_data):
    g = models.geom_mod()
    assert g.pmf(x=generate_data[0], p=0.5, pmin=0, pmax=0.5).all() == np.array(
        [0, 0.5, 0, 0.5, 0.5, 0.5, 0.5, 0.16666667, 0.5, 0]).all()


def test_geom_log_pmf(generate_data):
    g = models.geom_mod()
    assert g.log_pmf(x=generate_data[0], p=0.5, pmin=0, pmax=0.5).all() == np.array(
        [-np.inf, -0.69314718, -np.inf, -0.69314718, -0.69314718, -0.69314718, -0.69314718, -1.79175947, -0.69314718, -np.inf]).all()


def test_geom_optim(generate_data):
    g = models.geom_mod()
    optim(function=g.pmf,
          parameters=g.kwds,
          xdata=generate_data[1],
          ydata=generate_data[2])


def test_unif_pmf(generate_data):
    u = models.unif_mod()
    assert u.pmf(x=generate_data[0], pmin=0.1).all() == np.array(
        [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]).all()


def test_unif_log_pmf(generate_data):
    u = models.unif_mod()
    assert u.log_pmf(x=generate_data[0], pmin=0.1).all() == np.array([-2.30258509, -2.30258509, -2.30258509, -2.30258509, -2.30258509,
                                                                      -2.30258509, -2.30258509, -2.30258509, -2.30258509, -2.30258509]).all()
