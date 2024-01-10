# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 13:11:09 2024

@author: peruzzetto
"""

import numpy as np

import digdem.surfmod
import digdem.section


def test_Section_creation():
    point1 = (0.1, 0.0)
    point2 = (1.2, 0.0)
    sec_long = digdem.section.Section(
        [point1, point2],
        "A-B",
        normalized=False,
    )
    assert np.abs(sec_long.line.length - 1.1) < 1e-6


def test_SurfMod_creation(synthetic_case, synthetic_profiles):
    new_dem = digdem.surfmod.SurfMod(
        synthetic_case["zold"],
        synthetic_case["mask"],
        xaxis=synthetic_case["x"],
        yaxis=synthetic_case["y"],
        sections=[
            synthetic_profiles["sec_long"],
            synthetic_profiles["sec_trans"],
        ],
        dinterp=0.05,
        smooth=0,
        sparse_interp=1,
    )
    ztest = new_dem.fz_old(
        synthetic_case["xtest"], synthetic_case["ytest"], grid=False
    )
    # print(ztest, synthetic_case["ztest"])
    # assert ztest == synthetic_case["ztest"]
    assert np.abs(ztest - synthetic_case["ztest"]) < 1e-6
    # assert ztest == synthetic_case["ztest"]


def test_SurfMod_interpolation(synthetic_case, synthetic_profiles):
    new_dem = digdem.surfmod.SurfMod(
        synthetic_case["zold"],
        synthetic_case["mask"],
        xaxis=synthetic_case["x"],
        yaxis=synthetic_case["y"],
        sections=[
            synthetic_profiles["sec_long"],
            synthetic_profiles["sec_trans"],
        ],
        dinterp=0.05,
        smooth=0,
        sparse_interp=1,
    )

    pt = synthetic_profiles["sec_long"].line.intersection(
        synthetic_profiles["sec_trans"].line
    )
    xint = pt.xy[0][0]
    yint = pt.xy[1][0]
    zint_old = new_dem.fz_old(xint, yint, grid=False)

    new_dem.set_intersections(["A-B"], ["C-D"], [zint_old - 0.2])
    new_dem.update_interpolation()

    znew = new_dem.fz_new(xint, yint)

    assert np.abs(zint_old - 0.2 - znew) < 1e-4
