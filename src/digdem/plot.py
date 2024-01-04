#!/localstorage/peruzzetto/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 11:49:44 2020

@author: peruzzetto
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LightSource

def hillshade(array, azimuth, angle_altitude,xaxis=None,yaxis=None,maxlight=[0,1]):
    
    if (xaxis is None) or (yaxis is None):    
        x, y = np.gradient(array)
    else:
        x,y=np.gradient(array,xaxis,yaxis,edge_order=2)
    slope = np.pi/2. - np.arctan(np.sqrt(x*x + y*y))
    aspect = np.arctan2(-x, y)
    azimuthrad = azimuth*np.pi / 180.
    altituderad = angle_altitude*np.pi / 180.
     
 
    shaded = np.sin(altituderad) * np.sin(slope)\
     + np.cos(altituderad) * np.cos(slope)\
     * np.cos(azimuthrad - aspect)
     
    out=255*(shaded + 1)/2
    out=(maxlight[1]-maxlight[0])*out+maxlight[0]
    return out

def plot_topo(z, x, y, contour_step=None, nlevels=25, level_min=None,
              step_contour_bold=0, contour_labels_properties=None,
              label_contour=True, contour_label_effect=None,
              axe=None, vert_exag=1,
              contours_prop=None, contours_bold_prop=None,
              figsize=(10, 10),
              interpolation=None,
              sea_level=0, sea_color=None, alpha=1, azdeg=315, altdeg=45,
              zmin=None, zmax=None,
              plot_terrain=False, cmap_terrain='gist_earth',
              blend_mode='overlay', fraction=1, ndv=0,
              plot_colorbar=False, colorbar_kwargs=None):
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
    plot_terrain : TYPE, optional
        DESCRIPTION. The default is False.
    cmap_terrain : TYPE, optional
        DESCRIPTION. The default is 'gist_earth'.
    blend_mode : TYPE, optional
        DESCRIPTION. The default is 'overlay'.
    fraction : TYPE, optional
        DESCRIPTION. The default is 1.
    ndv : TYPE, optional
        DESCRIPTION. The default is 0.
    plot_colorbar : TYPE, optional
        DESCRIPTION. The default is False.
    colorbar_kwargs : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    None.

    """
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    im_extent = [x[0]-dx/2,
                 x[-1]+dx/2,
                 y[0]-dy/2,
                 y[-1]+dy/2]
    ls = LightSource(azdeg=azdeg, altdeg=altdeg)

    if level_min is None:
        if contour_step is not None:
            level_min = np.ceil(np.nanmin(z)/contour_step)*contour_step
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
    axe.set_ylabel('Y (m)')
    axe.set_xlabel('X (m)')
    axe.set_aspect('equal')

    if plot_terrain:
        z2 = np.copy(z)
        if zmin is None:
            zmin = np.nanmin(z2)
        if zmax is None:
            zmax = np.nanmax(z2)
        z2[np.isnan(z2)] = ndv
        z2 = np.flip(z2.T, axis=0)
        rgb = ls.shade(z2, plt.get_cmap(cmap_terrain),
                       blend_mode=blend_mode, fraction=fraction,
                       vert_exag=vert_exag, dx=dx, dy=dy,
                       vmin=zmin, vmax=zmax)
        ind = z2 == ndv
        ind = np.tile(ind[:, :, np.newaxis], (1, 1, 4))
        rgb[ind] = np.nan
        alpha=1
        plt.imshow(rgb, extent=im_extent,
                   interpolation=interpolation, alpha=alpha)
        if plot_colorbar:
            if colorbar_kwargs is None:
                colorbar_kwargs = {}
            norm = Normalize(vmin=zmin, vmax=zmax)
            cc = colorbar(cm.ScalarMappable(norm=norm, cmap=cmap_terrain),
                          ax=axe, **colorbar_kwargs)
            if 'position' in colorbar_kwargs:
                if colorbar_kwargs['position'] in ['bottom', 'top']:
                    cc.ax.set_xlabel('Altitude (m)')
                else:
                    cc.ax.set_ylabel('Altitude (m)')
            else:
                cc.ax.set_ylabel('Altitude (m)')
    else:
        shaded_topo = ls.hillshade(np.flip(z.T, axis=0),
                                   vert_exag=vert_exag, dx=dx, dy=dy,
                                   fraction=fraction)
        axe.imshow(shaded_topo, cmap='gray', extent=im_extent,
                   interpolation=interpolation, alpha=alpha)

    if contours_prop is None:
        contours_prop = dict(alpha=0.5, colors='k',
                             linewidths=0.4)
    axe.contour(z.T, extent=im_extent,
                levels=levels,
                **contours_prop)

    if contours_bold_prop is None:
        contours_bold_prop = dict(alpha=0.5, colors='k',
                                  linewidths=0.6)

    if step_contour_bold > 0:
        lmin = np.ceil(np.nanmin(z)/step_contour_bold)*step_contour_bold
        levels = np.arange(lmin, np.nanmax(z), step_contour_bold)
        cs = axe.contour(z.T, extent=im_extent,
                         levels=levels,
                         **contours_bold_prop)
        if label_contour:
            if contour_labels_properties is None:
                contour_labels_properties = {}
            clbls = axe.clabel(cs, **contour_labels_properties)
            if contour_label_effect is not None:
                plt.setp(clbls, path_effects=contour_label_effect)

    if sea_color is not None:
        cmap_sea = mcolors.ListedColormap([sea_color])
        cmap_sea.set_under(color='w', alpha=0)
        mask_sea = (z.T <= sea_level)*1
        if mask_sea.any():
            axe.imshow(mask_sea, extent=im_extent, cmap=cmap_sea,
                       vmin=0.5, origin='lower', interpolation='none')