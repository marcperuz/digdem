#!/localstorage/peruzzetto/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 11:10:35 2020

@author: peruzzetto
"""

import numpy as np
import matplotlib.pyplot as plt

import digdem.surfmod
import digdem.section

xmesh, ymesh = np.meshgrid(
    np.linspace(0, 1, 300), np.linspace(0, 1, 300), indexing="ij"
)

# Initial topography
slope_angle = 45
zold = -np.tan(np.deg2rad(slope_angle)) * xmesh
zold = zold - np.min(zold)

# Mask where topography must be modified
xc = 0.5
yc = 0.5
xwidth = 0.25
ywidth = 0.15
tmp = 1 - (xmesh - xc) ** 2 / xwidth**2 - (ymesh - yc) ** 2 / ywidth**2
mask = tmp >= 0

# Create controling sections
point1 = (0.1, 0.5)
point2 = (0.9, 0.5)
sec_long = digdem.section.Section(
    [point1, point2], "long", normalized=False, dinterp=0.05
)
# sec_long.extreme_points_names=['a','b']
point1 = (0.5, 0.2)
point2 = (0.5, 0.8)
sec_trans = digdem.section.Section(
    [point1, point2], "trans", normalized=False, dinterp=0.05
)

# Create SurfMod

new_dem = digdem.surfmod.SurfMod(
    zold,
    mask,
    xaxis=xmesh[:, 0],
    yaxis=ymesh[0, :],
    sections=[sec_long, sec_trans],
    dinterp=0.05,
    smooth=0,
    sparse_interp=1,
)

# %%Ue only intersection for topography interpolation
new_dem.set_intersections(["long"], ["trans"], [0.4])

new_dem.update_interpolation()

new_dem.plot(
    figsize=None,
    axe=None,
    contour_step=0.05,
    azdeg=315,
    altdeg=45,
    contour_properties={},
    section_properties={"lw": 3},
    text_properties={},
    control_points_properties=None,
    point_name_dist=None,
    npts=2,
)

new_dem["long"].plot()
new_dem["trans"].plot()

# %% Change intersection point height

new_dem.set_intersections(["long"], ["trans"], [0.45])

new_dem.update_interpolation()

new_dem.plot(
    figsize=None,
    axe=None,
    contour_step=0.05,
    azdeg=315,
    altdeg=45,
    contour_properties={},
    section_properties={"lw": 3},
    text_properties={},
    control_points_properties=None,
    point_name_dist=None,
    npts=2,
)

new_dem["long"].plot()
new_dem["trans"].plot()

# %% Add control points along sections

sec_long.set_control_points([0.3, 0.5], [0.48, 0.3])
sec_trans.set_control_points([0.2, 0.4], [0.42, 0.42])
# new_dem.set_intersections(['long'],['trans'],[0.45])

new_dem.update_interpolation()

new_dem.plot(
    figsize=None,
    axe=None,
    contour_step=0.05,
    azdeg=315,
    altdeg=45,
    contour_properties={},
    section_properties={"lw": 3},
    text_properties={},
    control_points_properties=None,
    point_name_dist=None,
    npts=2,
)

new_dem["long"].plot()
new_dem["trans"].plot()

plt.show()
