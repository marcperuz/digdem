# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 13:55:55 2019

@author: peruzzetto
"""

import itertools
import warnings
import logging
import scipy.interpolate

import shapely.geometry as geom
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd

from collections import defaultdict
from digdem.section import Section
from scipy.signal import convolve2d

import digdem.plot
import digdem.gis


class Tree(defaultdict):
    def __init__(self, value=None):
        super(Tree, self).__init__(Tree)
        self.value = value


class InterpPoints:
    def __init__(self, *args):
        if len(args) == 1:
            # If only one argument, it is assumed to be a text file containing
            # points coordinates separated by spaces and three columns
            # for x, y and z
            data = np.loadtxt(args[0])
            if data.ndim == 1:
                data = data[np.newaxis, :]
            args = [list(data[:, 0]), list(data[:, 1]), list(data[:, 2])]
        if len(args) == 2:
            self.points = args[0]
            self.z = args[1]
        elif len(args) == 3:
            self.z = args[2]
            assert len(args[0]) == len(args[1])
            points = [[args[0][i], args[1][i]] for i in range(len(args[0]))]
            self.points = geom.MultiPoint(points)
        else:
            self.points = geom.MultiPoint()
            self.z = []
        assert len(self.points.geoms) == len(self.z)

    @property
    def xy(self):
        return np.array(geom.LineString(self.points.geoms).coords)

    def add(self, *args):
        args = list(args)
        for i in range(len(args)):
            if not isinstance(args[i], list):
                args[i] = [args[i]]
        if len(args) == 3:
            x = args[0]
            y = args[1]
            z = args[2]
            points = [geom.Point([x[i], y[i]]) for i in range(len(x))]
        elif len[args] == 2:
            z = args[1]
            points = args[0]
        self.z = self.z + z
        self.points = geom.MultiPoint(list(self.points.geoms) + points)

    def unique(self):
        xy = self.xy
        xy, ind = np.unique(xy, axis=0, return_index=True)
        self.z = [self.z[i] for i in ind]
        self.points = geom.MultiPoint(xy)

    def __len__(self):
        return len(self.points.geoms)


class SurfMod(dict):
    @staticmethod
    def read_shapefile(
        shapefile,
        direction_field_name="dir",
        orientation_field_name="ori",
        normalized=False,
    ):
        sections = []
        try:
            assert shapefile[-4:] == ".shp"
            shp_reader = gpd.read_file(shapefile)
            # print(shapefile)
            for row in shp_reader.index:
                line = shp_reader.loc[row, "geometry"]
                direction = None
                orientation = None
                if direction_field_name in shp_reader.columns:
                    direction = shp_reader.loc[row, direction_field_name]
                if orientation_field_name in shp_reader.columns:
                    orientation = shp_reader.loc[row, orientation_field_name]
                name = shp_reader.loc[row, "name"]
                sections.append(
                    Section(
                        line,
                        name,
                        direction=direction,
                        orientation=orientation,
                        normalized=normalized,
                    )
                )
        except AssertionError:
            print("Format of " + shapefile + " file not recognized")
            exit(1)

        return sections

    @staticmethod
    def read_control_points_from_textfile(textfile, sep=" "):
        pos = {}
        z = {}
        z_intersections = Tree()
        with open(textfile, "r") as file:
            lines = [line.strip() for line in file if line.strip]
        assert (
            lines[0][0] == "#"
        ), "File must begin with a section name, identified with #"
        for line in lines:
            if line[-1] == "\n":
                line = line[:-1]
            if line[0] == "#":
                key = line[1:]
                pos[key] = []
                z[key] = []
            else:
                vals = line.split(sep)
                try:
                    if vals[0][0].isdigit():
                        pos[key].append(float(vals[0]))
                        z[key].append(float(vals[1]))
                    else:
                        z_intersections[key][vals[0]] = float(vals[1])
                except ValueError:
                    pass
        return pos, z, z_intersections

    def __init__(
        self,
        surf_old,
        mask,
        xaxis=None,
        yaxis=None,
        indexing="ij",
        sections=[],
        contour=None,
        normalized=False,
        sparse_interp=1,
        dinterp=10,
        ndv=-9999.9,
        smooth=0.5,
        surf_max=None,
        surf_min=None,
    ):
        self.logger = logging.getLogger(__name__)

        try:
            assert indexing in ["raster", "ij"], (
                "indexing must be " "raster" " or " "ij" ""
            )
            self.indexing = indexing
            if indexing == "raster":
                surf_old = np.flip(surf_old, axis=0).T
                mask = np.flip(mask, axis=0).T
                if surf_max is not None:
                    surf_max = np.flip(surf_max, axis=0).T
                if surf_min is not None:
                    surf_min = np.flip(surf_min, axis=0).T
            elif indexing == "ij":
                pass

            if xaxis is None:
                xaxis = np.arange(0, surf_old.shape[0], 1)
            else:
                assert isinstance(xaxis, np.ndarray), "xaxis is not np.ndarray"
            if yaxis is None:
                yaxis = np.arange(0, surf_old.shape[1], 1)
            else:
                assert isinstance(yaxis, np.ndarray), "yaxis is not np.ndarray"

            sizex_coherent = (
                mask.shape[0] == surf_old.shape[0] == xaxis.shape[0]
            )
            sizey_coherent = (
                mask.shape[1] == surf_old.shape[1] == yaxis.shape[0]
            )
            assert sizex_coherent & sizey_coherent, "Uncoherent shapes"

            self.surf_old = np.copy(surf_old)
            self.surf_old[np.isnan(self.surf_old)] = ndv
            self.surf_new = np.copy(self.surf_old)
            self.xaxis = xaxis
            self.yaxis = yaxis
            self.dx = np.abs(xaxis[1] - xaxis[0])
            self.dy = np.abs(yaxis[1] - yaxis[0])
            self.mask = np.copy(mask)
            self.mask[np.isnan(self.mask)] = 0
            self.normalized = normalized
            self.sparse_interp = sparse_interp
            self.intersections = Tree()
            self.ndv = ndv
            self.interp_points = InterpPoints()
            self.sup_interp_points = InterpPoints()
            self.dinterp = dinterp
            self.fz_old = None
            self.fz_new = None
            self.smooth = smooth
            self.surf_max = surf_max
            self.surf_min = surf_min

            self.contour = contour
            self.sections = sections

            if self.contour is None:
                self.set_contour()

            if len(self.sections) > 0:
                self.set_intersections()

            self.set_fz_old()

        except AssertionError as error:
            self.logger.error(error)
            raise

    @property
    def surf_min(self):
        return self._surf_min

    @surf_min.setter
    def surf_min(self, array):
        if array is None:
            self._surf_min = array
        else:
            if self.indexing == "raster":
                self._surf_min = np.flip(array, axis=1).T
            elif self.indexing == "ij":
                self._surf_min = array

    @surf_min.getter
    def surf_min(self):
        if self._surf_min is None:
            return self._surf_min
        else:
            if self.indexing == "raster":
                return np.flip(self._surf_min.T, axis=1)
            elif self.indexing == "ij":
                return self._surf_min

    @property
    def surf_max(self):
        return self._surf_max

    @surf_max.setter
    def surf_max(self, array):
        if array is None:
            self._surf_max = array
        else:
            if self.indexing == "raster":
                self._surf_max = np.flip(array, axis=1).T
            elif self.indexing == "ij":
                self._surf_max = array

    @surf_max.getter
    def surf_max(self):
        if self._surf_max is None:
            return self._surf_max
        else:
            if self.indexing == "raster":
                return np.flip(self._surf_max.T, axis=1)
            elif self.indexing == "ij":
                return self._surf_max

    @property
    def sections(self):
        return self._sections

    @sections.setter
    def sections(self, sections, **kwargs):
        if isinstance(sections, list):
            prof = []
            for section in sections:
                if isinstance(section, Section):
                    prof.append(section)
        elif isinstance(sections, str):
            if sections[-4:] == ".shp":
                kwargs["normalized"] = self.normalized
                prof = SurfMod.read_shapefile(sections, **kwargs)
            else:
                warnings.warn(
                    "sections initiation : file format not recognized"
                )
        self._sections = prof

        self.update_sections_extents()

        for section in self._sections:
            self[section.name] = section
            section.dinterp = self.dinterp

        for name1, name2 in itertools.product(self, self):
            for ff in ["points", "pos", "z"]:
                self.intersections[name1][name2][ff] = []

    @property
    def sections_extent(self):
        if len(self.sections) > 0:
            sec = self.sections[0]
            xx, yy = sec.line.xy
            xmin = min(xx)
            xmax = max(xx)
            ymin = min(yy)
            ymax = max(yy)
            for sec in self.sections[1:]:
                xx, yy = sec.line.xy
                xmin = min(xmin, min(xx))
                xmax = max(xmax, max(xx))
                ymin = min(ymin, min(yy))
                ymax = max(ymax, max(yy))
        else:
            xmin = None
            ymin = None
            xmax = None
            ymax = None
        return xmin, xmax, ymin, ymax

    def set_fz_old(self, finterp=None, outside_mask=False):
        if finterp is None:
            xmin, xmax, ymin, ymax = self.sections_extent
            indx = np.argwhere((self.xaxis >= xmin) & (self.xaxis <= xmax))[
                :, 0
            ]
            indy = np.argwhere((self.yaxis >= ymin) & (self.yaxis <= ymax))[
                :, 0
            ]
            xinterp = self.xaxis[indx[0] : indx[-1] + 1 : self.sparse_interp]
            yinterp = self.yaxis[indy[0] : indy[-1] + 1 : self.sparse_interp]
            zinterp = self.surf_old[
                indx[0] : indx[-1] + 1 : self.sparse_interp,
                indy[0] : indy[-1] + 1 : self.sparse_interp,
            ]
            if outside_mask:
                xinterp2, yinterp2 = np.meshgrid(
                    xinterp, yinterp, indexing="ij"
                )
                # Indexs inside the mask
                mask2 = self.mask[
                    indx[0] : indx[-1] + 1 : self.sparse_interp,
                    indy[0] : indy[-1] + 1 : self.sparse_interp,
                ]
                ind_mask = mask2 > 0
                # Indexes of mask contour
                mask3 = np.copy(mask2)
                mask3[ind_mask] = 1
                mask3 = convolve2d(
                    mask3,
                    np.ones((self.sparse_interp, self.sparse_interp)),
                    mode="same",
                )
                ind_cont = (mask3 > 0) & (mask2 == 0)
                # Interpoling function using contour
                fz_tmp = scipy.interpolate.Rbf(
                    xinterp2[ind_cont],
                    yinterp2[ind_cont],
                    zinterp[ind_cont],
                    function="multiquadric",
                    smooth=self.smooth,
                )
                # Modify altitudes within the mask for fz_old definition
                zinterp[ind_mask] = fz_tmp(
                    xinterp2[ind_mask], yinterp2[ind_mask]
                )
            self.fz_old = scipy.interpolate.RectBivariateSpline(
                xinterp, yinterp, zinterp
            )
        else:
            self.fz_old = finterp
        for section in self.sections:
            section.fz_old = self.fz_old
            if len(section.extent.pos) > 0:
                x0, y0 = section.extent.points.geoms[0].xy
                x1, y1 = section.extent.points.geoms[1].xy
                section.extent.z = [self.ndv, self.ndv]
                section.extent.z[0] = self.fz_old(x0[0], y0[0], grid=False)
                section.extent.z[-1] = self.fz_old(x1[0], y1[0], grid=False)

    def set_fz_new(self):
        if len(self.interp_points.z) > 0:
            self.fz_new = scipy.interpolate.Rbf(
                self.interp_points.xy[:, 0],
                self.interp_points.xy[:, 1],
                self.interp_points.z,
                function="multiquadric",
                smooth=self.smooth,
            )

    def set_surf_new(self):
        if self.fz_new is None:
            self.set_fz_new()
        if self.fz_new is not None:
            ind = self.mask > 0
            xmesh, ymesh = np.meshgrid(self.xaxis, self.yaxis, indexing="ij")
            self.surf_new[ind] = self.fz_new(xmesh[ind], ymesh[ind])
            if self.surf_min is not None:
                self.surf_new[ind] = np.maximum(
                    self.surf_new[ind], self.surf_min[ind]
                )
            if self.surf_max is not None:
                self.surf_new[ind] = np.minimum(
                    self.surf_new[ind], self.surf_max[ind]
                )

    def add_section(self, section):
        try:
            assert isinstance(
                section, Section
            ), "section is not a Section instance"
            self.sections.append(section)
            self[section.name] = section
            if self.contour is not None:
                self[section.name].set_extent(self.contour)
        except AssertionError:
            raise

    def update_sections_extents(self):
        if self.contour is not None:
            for section in self.sections:
                section.set_extent(self.contour)

    def set_contour(self, mask=None, x=None, y=None):
        plt.ioff()
        if mask is None:
            mask = self.mask
        if x is None:
            x = self.xaxis
        if y is None:
            y = self.yaxis

        try:
            backend = plt.get_backend()
            plt.switch_backend("Agg")
            fig, axe = plt.subplots(1, 1)
            cs = axe.contour(x, y, mask.T, [0.5])
            try:
                # Correct formulation for Matplotlib >=3.8
                p = cs.get_paths()
            except AttributeError:
                # Ensure backwards compatibility for older versions of matplotlib
                p = cs.collections[0].get_paths()
            assert len(p) > 0, "mask is either all 1 or all 0"
            p = p[0]
            self.contour = geom.Polygon(p.vertices)
            self.update_sections_extents()
            plt.close("all")
            plt.switch_backend(backend)

        except ValueError:
            raise
        except TypeError:
            raise
        except AssertionError:
            raise
        plt.ion()

    def set_intersections(self, *args):
        names = [name for name in self]
        nsections = len(names)

        if len(args) == 0:
            for i in range(nsections):
                for j in range(i + 1, nsections):
                    intersect = self[names[i]].line.intersection(
                        self[names[j]].line
                    )
                    # print(intersect)
                    if isinstance(intersect, geom.Point):
                        intersect = geom.MultiPoint([intersect])
                    if isinstance(intersect, geom.LineString):
                        continue
                    self.intersections[names[i]][names[j]][
                        "points"
                    ] = intersect
                    self.intersections[names[j]][names[i]][
                        "points"
                    ] = intersect
                    posi = [
                        self[names[i]].line.project(
                            point, normalized=self.normalized
                        )
                        for point in intersect.geoms
                    ]
                    posj = [
                        self[names[j]].line.project(
                            point, normalized=self.normalized
                        )
                        for point in intersect.geoms
                    ]
                    self.intersections[names[i]][names[j]]["pos"] = posi
                    self.intersections[names[j]][names[i]]["pos"] = posj
                    self.intersections[names[i]][names[j]]["z"] = [
                        np.nan for point in intersect.geoms
                    ]
                    self.intersections[names[j]][names[i]]["z"] = [
                        np.nan for point in intersect.geoms
                    ]
        #                    print(names[i],names[j],posi)
        #                    print(names[j],names[i],posj)
        else:
            names1 = args[0]
            names2 = args[1]
            zs = args[2]
            assert len(names1) == len(names2)
            assert len(names1) == len(zs)
            for key1, key2, z in zip(names1, names2, zs):
                if (key1 in names) and (key2 in names):
                    self.intersections[key1][key2]["z"] = [z]
                    self.intersections[key2][key1]["z"] = [z]

                    self[key1].set_intersection_points(
                        key2,
                        self.intersections[key1][key2]["pos"],
                        self.intersections[key1][key2]["z"],
                    )
                    self[key2].set_intersection_points(
                        key1,
                        self.intersections[key2][key1]["pos"],
                        self.intersections[key2][key1]["z"],
                    )

    def set_control_points(self, *args, sep=" "):
        if len(args) == 2:
            pos = args[0]
            z = args[1]
        elif len(args) == 1:
            (
                pos,
                z,
                z_intersections,
            ) = SurfMod.read_control_points_from_textfile(args[0], sep=sep)
            for name1 in z_intersections:
                for name2 in z_intersections[name1]:
                    zz = z_intersections[name1][name2]
                    self.intersections[name1][name2]["z"] = [zz]
                    self.intersections[name2][name1]["z"] = [zz]
        else:
            pos = []
            z = []

        if (isinstance(pos, dict)) & (isinstance(z, dict)):
            for key in self:
                if key in pos:
                    self[key].set_control_points(pos[key], z[key])
                for key2 in self.intersections[key]:
                    self[key].set_intersection_points(
                        key2,
                        self.intersections[key][key2]["pos"],
                        self.intersections[key][key2]["z"],
                    )
                #                        print(key,key2,self.intersections[key][key2]['pos'],self.intersections[key][key2]['z'])
                self[key].filter_control_points()

        else:
            raise ValueError("set_control_points recieved incorrect arguments")

    def set_interp_points(self):
        interp_points = []
        interp_z = []
        if self.contour is None:
            self.set_contour()

        # Add points on the contour
        if self.contour is not None:
            exterior = self.contour.exterior
            npts = int(np.floor(exterior.length / self.dinterp) + 1)
            ss = np.linspace(0, exterior.length, npts)
            new_pts = [exterior.interpolate(s) for s in ss]
            interp_points = interp_points + new_pts
            coords = [pt.coords[:][0] for pt in new_pts]
            interp_z = interp_z + [
                self.fz_old(coord[0], coord[1], grid=False) for coord in coords
            ]

        # Add points from the sections
        for section in self.sections:
            section.set_interp_spline()
            section.set_interp_points()
            interp_points = interp_points + [
                pp for pp in section.interp_points.points.geoms
            ]
            interp_z = interp_z + section.interp_points.z

        # Add supplementary control points
        interp_points = interp_points + [
            pp for pp in self.sup_interp_points.points.geoms
        ]
        interp_z = interp_z + self.sup_interp_points.z

        self.interp_points.points = geom.MultiPoint(interp_points)
        self.interp_points.z = interp_z
        self.interp_points.unique()

    def update_interpolation(self):
        self.set_interp_points()
        self.set_fz_new()
        self.set_surf_new()

    def set_sup_interp_points(self, *args):
        self.sup_interp_points = InterpPoints(*args)

    def add_sup_interp_points(self, *args):
        self.sup_interp_points.add(*args)

    def save_to_raster(self, file_out, which="new", fmt=None, **kwargs):
        if which == "old":
            z = self.surf_old.copy()
        elif which == "new":
            z = self.surf_new.copy()
        z = np.flip(z.T, axis=0)
        digdem.gis.write_raster(
            self.xaxis, self.yaxis, z, file_out, fmt=fmt, **kwargs
        )

    def plot(
        self,
        figsize=None,
        axe=None,
        contour_properties=None,
        section_properties=None,
        text_properties=None,
        control_points_properties=None,
        plot_interp_points=False,
        interp_points_properties=None,
        point_name_dist=None,
        npts=2,
        **kwargs,
    ):
        if axe is None:
            fig, axe = plt.subplots(1, 1, figsize=figsize)
        else:
            fig = axe.figure

        # Plot topography

        digdem.plot.plot_topo(
            self.surf_new,
            self.xaxis,
            self.yaxis,
            ndv=self.ndv,
            # indexing=self.indexing,
            axe=axe,
            **kwargs,
        )

        # Plot contour
        if contour_properties is None:
            contour_properties = dict()
        xc, yc = self.contour.exterior.coords.xy
        axe.plot(xc, yc, **contour_properties)

        # Plot sections
        if control_points_properties is None:
            control_points_properties = {
                "marker": "o",
                "ms": 5,
                "mec": "white",
                "mfc": "k",
            }
        for section in self.sections:
            section.plot_on_map(
                axe=axe,
                section_properties=section_properties,
                text_properties=text_properties,
                control_points_properties=control_points_properties,
                npts=npts,
                point_name_dist=point_name_dist,
            )

        for pt in self.sup_interp_points.points.geoms:
            axe.plot(
                pt.coords[0][0], pt.coords[0][1], **control_points_properties
            )

        if plot_interp_points:
            if interp_points_properties is None:
                interp_points_properties = {
                    "linestyle": "None",
                    "marker": "o",
                    "ms": 3,
                    "mfc": "w",
                    "mec": "k",
                }
            for pt in self.interp_points.points.geoms:
                axe.plot(
                    pt.coords[0][0],
                    pt.coords[0][1],
                    **interp_points_properties,
                )

        return fig, axe
