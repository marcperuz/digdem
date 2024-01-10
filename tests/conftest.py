# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 16:19:18 2023

@author: peruzzetto
"""

import pytest

import numpy as np

import digdem.section


@pytest.fixture
def synthetic_case():
    nx = 301
    ny = 201
    x = np.linspace(0, 1.5, 301)
    y = np.linspace(-0.5, 0.5, 201)
    xmesh, ymesh = np.meshgrid(x, y, indexing="ij")
    # Initial topography
    slope_angle = 45
    zold = -np.tan(np.deg2rad(slope_angle)) * xmesh + 0.5 * ymesh**2
    zold = zold - np.min(zold)

    # Mask
    xc = 0.6
    yc = 0.0
    xwidth = 0.30
    ywidth = 0.20
    tmp = 1 - (xmesh - xc) ** 2 / xwidth**2 - (ymesh - yc) ** 2 / ywidth**2
    mask = tmp >= 0

    xtest = x[nx // 2]
    ytest = y[ny // 2]
    ztest = zold[nx // 2, ny // 2]

    return dict(
        x=x, y=y, zold=zold, mask=mask, xtest=xtest, ytest=ytest, ztest=ztest
    )


@pytest.fixture
def synthetic_profiles():
    # Longitudinal cross-section
    point1 = (0.1, 0.0)
    point2 = (1.2, 0.0)
    sec_long = digdem.section.Section(
        [point1, point2], "A-B", normalized=False, dinterp=0.05
    )

    # Transversal cross-section
    point1 = (0.6, -0.3)
    point2 = (0.6, 0.3)
    sec_trans = digdem.section.Section(
        [point1, point2], "C-D", normalized=False, dinterp=0.05
    )
    return dict(sec_long=sec_long, sec_trans=sec_trans)
