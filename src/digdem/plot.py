#!/localstorage/peruzzetto/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 11:49:44 2020

@author: peruzzetto
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

BOLD_CONTOURS_INTV = [
    0.1,
    0.2,
    0.5,
    1,
    2.0,
    5,
    10,
    20,
    50,
    100,
    200,
    500,
    1000,
]
NB_THIN_CONTOURS = 10
NB_BOLD_CONTOURS = 3


def get_contour_intervals(
    zmin, zmax, nb_bold_contours=None, nb_thin_contours=None
):
    if nb_thin_contours is None:
        nb_thin_contours = NB_THIN_CONTOURS
    if nb_bold_contours is None:
        nb_bold_contours = NB_BOLD_CONTOURS

    intv = (zmax - zmin) / nb_bold_contours
    i = np.argmin(np.abs(np.array(BOLD_CONTOURS_INTV) - intv))

    bold_intv = BOLD_CONTOURS_INTV[i]
    if BOLD_CONTOURS_INTV[i] != BOLD_CONTOURS_INTV[0]:
        if bold_intv - intv > 0:
            bold_intv = BOLD_CONTOURS_INTV[i - 1]

    if nb_thin_contours is None:
        thin_intv = bold_intv / NB_THIN_CONTOURS
        if (zmax - zmin) / bold_intv > 5:
            thin_intv = thin_intv * 2
    else:
        thin_intv = bold_intv / nb_thin_contours

    return bold_intv, thin_intv


def plot_topo(
    z,
    x,
    y,
    indexing="ij",
    contour_step=None,
    nlevels=None,
    level_min=None,
    step_contour_bold="auto",
    contour_labels_properties=None,
    label_contour=True,
    contour_label_effect=None,
    axe=None,
    vert_exag=1,
    fraction=1,
    ndv=-9999,
    uniform_grey=None,
    contours_prop=None,
    contours_bold_prop=None,
    figsize=None,
    interpolation=None,
    sea_level=0,
    sea_color=None,
    alpha=1,
    azdeg=315,
    altdeg=45,
    zmin=None,
    zmax=None,
):
    """
    Plot topography with hillshading.

    Parameters
    ----------
    z : TYPE
        DESCRIPTION.
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    contour_step : TYPE, optional
        DESCRIPTION. The default is None.
    nlevels : TYPE, optional
        DESCRIPTION. The default is 25.
    contour_labels_properties : TYPE, optional
        DESCRIPTION. The default is None.
    axe : TYPE, optional
        DESCRIPTION. The default is None.
    vert_exag : TYPE, optional
        DESCRIPTION. The default is 1.
    contour_color : TYPE, optional
        DESCRIPTION. The default is 'k'.
    contour_alpha : TYPE, optional
        DESCRIPTION. The default is 1.
    figsize : TYPE, optional
        DESCRIPTION. The default is (10,10).
    interpolation : TYPE, optional
        DESCRIPTION. The default is None.
    sea_level : TYPE, optional
        DESCRIPTION. The default is 0.
    sea_color : TYPE, optional
        DESCRIPTION. The default is None.
    alpha : TYPE, optional
        DESCRIPTION. The default is 1.
    azdeg : TYPE, optional
        DESCRIPTION. The default is 315.
    altdeg : TYPE, optional
        DESCRIPTION. The default is 45.
    fraction : TYPE, optional
        DESCRIPTION. The default is 1.
    ndv : TYPE, optional
        DESCRIPTION. The default is -9999.

    Returns
    -------
    None.

    """
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    im_extent = [x[0] - dx / 2, x[-1] + dx / 2, y[0] - dy / 2, y[-1] + dy / 2]
    ls = mcolors.LightSource(azdeg=azdeg, altdeg=altdeg)

    if indexing == "ij":
        z = np.flip(z.T, axis=0)

    z[z == ndv] = np.nan
    auto_bold_intv = None

    if nlevels is None and contour_step is None:
        auto_bold_intv, contour_step = get_contour_intervals(
            np.nanmin(z), np.nanmax(z)
        )

    if level_min is None:
        if contour_step is not None:
            level_min = np.ceil(np.nanmin(z) / contour_step) * contour_step
        else:
            level_min = np.nanmin(z)
    if contour_step is not None:
        levels = np.arange(level_min, np.nanmax(z), contour_step)
    else:
        levels = np.linspace(level_min, np.nanmax(z), nlevels)

    if axe is None:
        fig = plt.figure(figsize=figsize)
        axe = fig.gca()
    else:
        fig = axe.figure
    axe.set_ylabel("Y (m)")
    axe.set_xlabel("X (m)")
    axe.set_aspect("equal")

    if uniform_grey is None:
        shaded_topo = ls.hillshade(
            z, vert_exag=vert_exag, dx=dx, dy=dy, fraction=1
        )
    else:
        shaded_topo = np.ones(z.shape) * uniform_grey
    shaded_topo[z == ndv] = np.nan
    axe.imshow(
        shaded_topo,
        cmap="gray",
        extent=im_extent,
        interpolation=interpolation,
        alpha=alpha,
        vmin=0,
        vmax=1,
    )

    if contours_prop is None:
        contours_prop = dict(alpha=0.5, colors="k", linewidths=0.5)
    axe.contour(
        x,
        y,
        np.flip(z, axis=0),
        extent=im_extent,
        levels=levels,
        **contours_prop,
    )

    if contours_bold_prop is None:
        contours_bold_prop = dict(alpha=0.8, colors="k", linewidths=0.8)

    if step_contour_bold == "auto":
        if auto_bold_intv is None:
            auto_bold_intv, _ = get_contour_intervals(
                np.nanmin(z), np.nanmax(z)
            )
        step_contour_bold = auto_bold_intv

    if step_contour_bold > 0:
        lmin = np.ceil(np.nanmin(z) / step_contour_bold) * step_contour_bold
        if lmin < level_min:
            lmin = lmin + step_contour_bold
        levels = np.arange(lmin, np.nanmax(z), step_contour_bold)
        cs = axe.contour(
            x,
            y,
            np.flip(z, axis=0),
            extent=im_extent,
            levels=levels,
            **contours_bold_prop,
        )
        if label_contour:
            if contour_labels_properties is None:
                contour_labels_properties = {}
            clbls = axe.clabel(cs, **contour_labels_properties)
            if contour_label_effect is not None:
                plt.setp(clbls, path_effects=contour_label_effect)

    if sea_color is not None:
        cmap_sea = mcolors.ListedColormap([sea_color])
        cmap_sea.set_under(color="w", alpha=0)
        mask_sea = (z <= sea_level) * 1
        if mask_sea.any():
            axe.imshow(
                mask_sea,
                extent=im_extent,
                cmap=cmap_sea,
                vmin=0.5,
                origin="lower",
                interpolation="none",
            )
