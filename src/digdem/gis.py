# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 09:28:44 2024

@author: peruzzetto
"""

import importlib

import numpy as np
import geopandas as gpd

from pandas.api.types import is_numeric_dtype


def read_raster(file):
    if file.endswith(".asc") or file.endswith(".txt"):
        return read_ascii(file)
    elif file.endswith(".tif") or file.endswith(".tif"):
        return read_tiff(file)


def write_raster(x, y, z, file_out, fmt=None, **kwargs):
    # File format read from file_out overrides fmt
    fmt_from_file = file_out.split(".")
    if len(fmt_from_file) > 1:
        fmt = fmt_from_file[-1]
    else:
        file_out = file_out + "." + fmt

    if fmt not in ["asc", "ascii", "txt", "tif", "tiff"]:
        raise ValueError("File format not implemented in write_raster")

    if fmt.startswith("tif"):
        if importlib.util.find_spec("rasterio") is None:
            print(
                (
                    "rasterio is required to write tif files.",
                    " Switching to asc format",
                )
            )
            fmt = "asc"

    if fmt in ["asc", "ascii", "txt"]:
        write_ascii(x, y, z, file_out)
    elif fmt in ["tif", "tiff"]:
        write_tiff(x, y, z, file_out, **kwargs)
    else:
        raise NotImplementedError()


def read_tiff(file):
    import rasterio

    with rasterio.open(file, "r") as src:
        dem = src.read(1)
        ny, nx = dem.shape
        x = np.linspace(src.bounds.left, src.bounds.right, nx)
        y = np.linspace(src.bounds.bottom, src.bounds.top, ny)
    return x, y, dem


def read_ascii(file):
    """
    Read ascii grid file to numpy ndarray.

    Parameters
    ----------
    file : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    dem = np.loadtxt(file, skiprows=6)
    grid = {}
    with open(file, "r") as fid:
        for i in range(6):
            tmp = fid.readline().split()
            grid[tmp[0]] = float(tmp[1])
    try:
        x0 = grid["xllcenter"]
        y0 = grid["yllcenter"]
    except KeyError:
        x0 = grid["xllcorner"]
        y0 = grid["yllcorner"]
    nx = int(grid["ncols"])
    ny = int(grid["nrows"])
    dx = dy = grid["cellsize"]
    x = np.linspace(x0, x0 + (nx - 1) * dx, nx)
    y = np.linspace(y0, y0 + (ny - 1) * dy, ny)

    return x, y, dem


def write_tiff(x, y, z, file_out, raster_ref=None, **kwargs):
    """
    Write raster as tif file

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    z : TYPE
        DESCRIPTION.
    file_out : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    import rasterio
    from rasterio.transform import Affine

    if "driver" not in kwargs:
        kwargs["driver"] = "GTiff"

    if raster_ref is None:
        res = (x[-1] - x[0]) / (len(x) - 1)
        transform = Affine.translation(
            x[0] - res / 2, y[-1] - res / 2
        ) * Affine.scale(res, -res)
        profile = dict(
            height=z.shape[0],
            width=z.shape[1],
            count=1,
            dtype=z.dtype,
            transform=transform,
            **kwargs
        )
    else:
        profile = rasterio.open(raster_ref).profile

    with rasterio.open(file_out, "w", **profile) as dst:
        dst.write(z, 1)


def write_ascii(x, y, z, file_out):
    """
    Write raster as ascii file.

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    z : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    nx = z.shape[1]
    ny = z.shape[0]
    cellsize = x[1] - x[0]
    header_txt = (
        "ncols {:.0f}\nnrows {:.0f}\nxllcorner {:.5f}\nyllcorner {:.5f}\n"
    )
    header_txt += "cellsize {:.4f}\nnodata_value -99999"
    header_txt = header_txt.format(nx, ny, x[0], y[0], cellsize)
    np.savetxt(file_out, z, header=header_txt, comments="")


def rasterize_shapefile(
    shapefile,
    file_out,
    field=None,
    lookup_table=None,
    raster_ref=None,
    ndv=-99999,
    resolution=None,
):
    import rasterio
    import rasterio.features

    if isinstance(shapefile, str):
        shapefile = gpd.read_file(shapefile)

    if field is None:
        field = shapefile.columns[0]

    if field == "id":
        shapefile["id"] = np.arange(1, len(shapefile) + 1)

    # Remove null values
    shapefile = shapefile[~shapefile[field].isnull()]

    if not is_numeric_dtype(shapefile[field]):
        if lookup_table is None:
            # if no lookup table is given, all values are considered and
            # numbered by ascending order
            vals = shapefile[field].unique()
            vals.sort()
            lookup_table = dict()
            for i in range(len(vals)):
                lookup_table[vals[i]] = i + 1
        else:
            # Otherwise, shapes with field values not in the
            # lookup table are removed
            shapefile = shapefile[shapefile[field].isin(lookup_table.keys())]
        field_lookup = "lookup"
        # Add field with corresponding ids
        shapefile[field_lookup] = ndv
        for key in lookup_table:
            shapefile.loc[
                shapefile[field] == key, field_lookup
            ] = lookup_table[key]
    else:
        field_lookup = field

    if raster_ref is not None:
        # If raster is given as reference, use it to get the characteristics
        # of the new raster
        rst = rasterio.open(raster_ref)
        meta = rst.meta.copy()
        meta.update(compress="lzw")
        with rasterio.open(file_out, "w+", **meta) as out:
            out_arr = out.read(1)

            # this is where we create a generator of geom,
            # value pairs to use in rasterizing
            shapes = (
                (geom, value)
                for geom, value in zip(
                    shapefile.geometry, shapefile[field_lookup]
                )
            )

            burned = rasterio.features.rasterize(
                shapes=shapes, fill=ndv, out=out_arr, transform=out.transform
            )
            if "nodata" in meta:
                burned[burned == meta["nodata"]] = np.nan
            out.write_band(1, burned)

    else:
        if resolution is None:
            raise TypeError("resolution must be a strictly positive number")
        from geocube.api.core import make_geocube

        cube = make_geocube(
            shapefile,
            measurements=[field_lookup],
            resolution=(resolution, -resolution),
        )
        cube.lookup.rio.to_raster("file_out.tif")

    return lookup_table
