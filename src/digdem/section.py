#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 21:16:48 2019

@author: peruzzetto
"""

import shapely
import warnings
import copy
import scipy.interpolate

import shapely.geometry as geom
import numpy as np
import matplotlib.pyplot as plt


# class Extent():
#    def __init__(self,coords=None,pos=None):
#        self.coords=coords
#        self.pos=pos


class SectionPoints:
    def __init__(self, *args, normalized=False, z=[]):
        self.normalized = normalized
        self.pos = []
        self.points = []
        self.z = []
        try:
            if len(args) == 2:
                linear_geom = args[0]
                self.line = linear_geom
                points = args[1]
                if isinstance(points, geom.LineString):
                    points = geom.MultiPoint(points.coords[:])

                if isinstance(points, geom.MultiPoint):
                    if not points.is_empty:
                        self.points = points
                        self.pos = [
                            linear_geom.project(point, normalized=normalized)
                            for point in points.geoms
                        ]
                        self.z = z
                elif isinstance(points, list):
                    self.pos = points
                    self.points = geom.MultiPoint(
                        [
                            linear_geom.interpolate(
                                s, normalized=self.normalized
                            )
                            for s in points
                        ]
                    )
                    self.z = z
                else:
                    warnings.warn(
                        "Second argument must be either shapely MultiPoint, LinearString or float list. \
                                  SectionPoints is created as empty"
                    )
                    self.set_empty()

                self.order_from_pos()

            else:
                try:
                    assert (
                        len(args) == 0
                    ), "0 or 2 arguments must be given for SectionPoints creation"
                    self.set_empty()
                except AssertionError:
                    raise
        except AssertionError:
            raise

    @property
    def z(self):
        return self._z

    @z.setter
    def z(self, z):
        assert isinstance(z, list), "z must be a float list"
        assert all(
            type(zz) in [float, int, np.float32, np.float64] for zz in z
        ), "z contains non int or float elements"
        if isinstance(self.pos, list):
            #            print(self.pos,z)
            assert (len(z) == 0) | (
                len(z) == len(self.pos)
            ), "length of z does not match points input"
        self._z = z

    def set_empty(self):
        self.points = geom.MultiPoint()
        self.pos = []
        self.z = []

    def order_from_pos(self):
        if len(self.pos) > 1:
            ind = sorted(range(len(self.pos)), key=lambda k: self.pos[k])
            self.points = geom.MultiPoint([self.points.geoms[i] for i in ind])
            self.pos = [self.pos[i] for i in ind]
            if len(self.z) > 1:
                self.z = [self.z[i] for i in ind]

    def filter_by_pos(self, posmin, posmax):
        tmp = np.array(self.pos)
        ind = np.argwhere((tmp > posmin) & (tmp < posmax))
        if len(ind) > 0:
            self.points = geom.MultiPoint(
                [self.points.geoms[i[0]] for i in ind]
            )
            self.pos = [self.pos[i[0]] for i in ind]
            if len(self.z) > 1:
                self.z = [self.z[i[0]] for i in ind]
        else:
            self.points = geom.MultiPoint()
            self.pos = []
            self.z = []

    def add_points(self, points, z=[]):
        if not isinstance(points, list):
            points = [points]

        points2 = []
        pos2 = []
        z2 = []

        for point in points:
            if type(point) in [float, int, np.float32, np.float64]:
                pos2.append(point)
                points2.append(
                    geom.Point(
                        self.line.interpolate(
                            point, normalized=self.normalized
                        )
                    )
                )
            elif point.isinstance(SectionPoints):
                try:
                    assert (
                        point.line == self.line
                    ), "Combining SectionPoints along different lines"
                    pos2 = pos2 + point.pos
                    points2 = points2 + [pp for pp in point.points]
                    z2 = z2 + point.z
                except AssertionError:
                    raise
            else:
                if isinstance(point, geom.Point):
                    point = [point]
                points2 = points2 + [pp for pp in point]
                pos2 = pos2 + [
                    self.line.project(pp, normalized=self.normalized)
                    for pp in point
                ]
        self.pos = self.pos + pos2
        self.points = geom.MultiPoint([pp for pp in self.points] + points2)
        if len(z2) > 0:
            self.z = self.z + z2
        else:
            self.z = self.z + z
        self.order_from_pos()


class Section:
    def __init__(
        self,
        points,
        name,
        orientation=None,
        direction="",
        normalized=False,
        topography=None,
        name_sep="-",
        fz_old=None,
        control_points=None,
        dinterp=10,
        interp_spline_options={},
        line_points_as_control_points=True,
        point_name_dist=None,
    ):
        if isinstance(points, geom.linestring.LineString):
            self.line = points
        elif isinstance(points, list):
            self.line = geom.LineString(points)
        elif isinstance(points, np.ndarray):
            tmp = [
                (points[i, 0], points[i, 1]) for i in range(points.shape[0])
            ]
            self.line = geom.LineString(tmp)

        self.name = name
        self.direction = direction
        self.normalized = normalized
        self.topography = topography
        self.fz_old = fz_old
        self.interp_spline = None
        self.interp_spline_options = interp_spline_options
        self.control_points = control_points
        self.intersection_points = {}
        self.orientation = orientation
        self.extent = None
        self.extreme_points_names = []
        self.name_sep = name_sep
        self.point_name_dist = point_name_dist
        self.dinterp = dinterp
        self.interp_points = SectionPoints()

        self.set_extreme_points_names()

        if (control_points is None) & line_points_as_control_points:
            self.control_points = SectionPoints(
                self.line, self.line, normalized=self.normalized
            )

    @property
    def orientation(self):
        return self._orientation

    @orientation.setter
    def orientation(self, orientation):
        self._orientation = orientation
        if isinstance(orientation, str):
            if orientation == "S-N":
                reverse = self.line.coords[0][1] > self.line.coords[-1][1]
            elif orientation == "N-S":
                reverse = self.line.coords[0][1] < self.line.coords[-1][1]
            elif orientation == "W-E":
                reverse = self.line.coords[0][0] > self.line.coords[-1][0]
            elif orientation == "E-W":
                reverse = self.line.coords[0][0] < self.line.coords[-1][0]
            else:
                warnings.warn(
                    "orientation not set, must be either 'S-N', 'N-S', 'W-E' or 'E-W'"
                )
            if reverse:
                self.line = shapely.ops.substring(
                    self.line, self.line.length, 0
                )

    def set_extreme_points_names(self, sep=None):
        if sep is None:
            sep = self.name_sep
        if sep != "":
            names = self.name.split(self.name_sep)
            if len(names) == 2:
                self.extreme_points_names = names
            else:
                warnings.warn(
                    "Extreme points names could not be set, invalid string separator"
                )

    def set_extent(self, other):
        if other is not None:
            intersections = self.line.intersection(other)
            self.extent = SectionPoints(
                self.line, intersections, normalized=self.normalized
            )
            self.filter_control_points()
            if self.fz_old is not None and len(self.extent.pos) > 0:
                coords = np.array(
                    geom.LineString(self.extent.points.geoms).coords
                )
                self.extent.z = list(
                    self.fz_old(coords[:, 0], coords[:, 1], grid=False)
                )

    def filter_control_points(self):
        if self.extent is not None:
            if (self.control_points is not None) & (len(self.extent.pos) > 1):
                self.control_points.filter_by_pos(
                    self.extent.pos[0], self.extent.pos[-1]
                )
            if len(self.extent.pos) == 0:
                self.control_points.set_empty()

    def filter_intersection_points(self):
        if self.extent is not None:
            for key in self.intersection_points:
                if (self.intersection_points[key] is not None) & (
                    len(self.extent.pos) > 1
                ):
                    self.intersection_points[key].filter_by_pos(
                        self.extent.pos[0], self.extent.pos[-1]
                    )
                if len(self.extent.pos) == 0:
                    self.intersection_points[key].set_empty()

    def set_control_points(self, pos, z):
        try:
            assert len(pos) == len(z), "pos and z must have same length"
            self.control_points = SectionPoints(
                self.line, pos, z=z, normalized=self.normalized
            )
            self.filter_control_points()
        except AssertionError:
            raise

    def set_intersection_points(self, name2, pos, z=[]):
        try:
            assert len(pos) == len(z), "pos and z must have same length"
            self.intersection_points[name2] = SectionPoints(
                self.line, pos, z=z, normalized=self.normalized
            )
            self.filter_intersection_points()
        except AssertionError:
            raise

    def set_interp_spline(self):
        ss = copy.copy(self.control_points.pos)
        zz = copy.copy(self.control_points.z)
        for key in self.intersection_points:
            ss = ss + copy.copy(self.intersection_points[key].pos)
            zz = zz + copy.copy(self.intersection_points[key].z)
        ind = sorted(range(len(ss)), key=lambda k: ss[k])
        ss.sort()
        zz = [zz[i] for i in ind]
        if len(ss) > 0:
            if self.extent is not None:
                if (len(self.extent.pos) > 0) & (len(self.extent.z) > 0):
                    ss.append(self.extent.pos[-1])
                    zz.append(self.extent.z[-1])
                    ss.insert(0, self.extent.pos[0])
                    zz.insert(0, self.extent.z[0])
            npts = len(ss)
            self.interp_spline = scipy.interpolate.splrep(
                ss, zz, s=0, k=min(3, npts - 1), **self.interp_spline_options
            )

    def set_interp_points(self):
        if len(self.interp_points.pos) > 0 or len(self.extent.pos) > 1:
            if self.extent is not None and len(self.extent.pos) > 1:
                #                print(self.name,len(self.extent.pos))
                pos1 = self.extent.pos[0]
                pos2 = self.extent.pos[-1]
            else:
                pos1 = 0
                pos2 = self.line.length
            pos = list(np.arange(pos1, pos2 + self.dinterp / 2, self.dinterp))
            if pos[-1] != pos2:
                pos.append(pos2)
            if self.interp_spline is None:
                z = []
            else:
                z = scipy.interpolate.splev(pos, self.interp_spline, der=0)
            self.interp_points = SectionPoints(self.line, pos, z=list(z))
        else:
            self.interp_points = SectionPoints()

    def plot(
        self,
        axe=None,
        figsize=None,
        npos=200,
        old_topo_properties=None,
        new_topo_properties=None,
        extent_properties=None,
        control_points_properties=None,
        fz_new=None,
        comp_surf=None,
        plot_old_topo=True,
        fz_old=None,
    ):
        if axe is None:
            fig, axe = plt.subplots(1, 1, figsize=figsize)
        else:
            fig = axe.figure

        if old_topo_properties is None:
            old_topo_properties = {"lw": 2, "color": "k", "linestyle": "-"}
        if new_topo_properties is None:
            new_topo_properties = {
                "lw": 2,
                "color": [0.5, 0.5, 0.5],
                "linestyle": "-",
            }

        if extent_properties is None:
            extent_properties = {
                "linestyle": "None",
                "marker": "o",
                "ms": 6,
                "mec": "black",
                "mfc": "white",
            }

        if control_points_properties is None:
            control_points_properties = {
                "linestyle": "None",
                "marker": "o",
                "ms": 6,
                "mec": "black",
                "mfc": [0.5, 0.5, 0.5],
            }

        pos = np.linspace(0, self.line.length, npos)
        pts = geom.LineString([self.line.interpolate(p) for p in pos])
        coords = np.array(pts.coords)

        if self.extent is not None and len(self.extent.pos) > 0:
            ind = np.argwhere(
                (pos >= self.extent.pos[0]) & (pos <= self.extent.pos[1])
            )
            ind = np.squeeze(ind)
        elif len(self.control_points.pos) > 1:
            ind = np.argwhere(
                (pos >= self.control_points.pos[0])
                & (pos <= self.control_points.pos[1])
            )
            ind = np.squeeze(ind)
        else:
            ind = []

        if fz_old is None:
            fz_old = self.fz_old

        if fz_old is not None:
            old_topo = fz_old(coords[:, 0], coords[:, 1], grid=False)
        else:
            old_topo = None

        if len(ind) > 0:

            def evaluate_new_topo(*args):
                if fz_new is None:
                    out = scipy.interpolate.splev(
                        pos[ind], self.interp_spline, der=0
                    )
                else:
                    out = fz_new(
                        coords[ind[0] : ind[-1] + 1, 0],
                        coords[ind[0] : ind[-1] + 1, 1],
                    )
                    if old_topo is not None:
                        if comp_surf == "min":
                            out = np.minimum(out, old_topo[ind])
                        elif comp_surf == "max":
                            out = np.maximum(out, old_topo[ind])
                return out

            if plot_old_topo:
                axe.plot(pos, old_topo, **old_topo_properties)
                new_topo = np.copy(old_topo)
                new_topo[ind] = evaluate_new_topo(old_topo)
                pos_new = pos
            else:
                new_topo = evaluate_new_topo()
                pos_new = pos[ind]

            axe.plot(pos_new, new_topo, **new_topo_properties)

            # Plot extent and control points
            if self.extent is not None:
                if len(self.extent.pos) > 1:
                    for i in [0, -1]:
                        axe.plot(
                            self.extent.pos[i],
                            self.extent.z[i],
                            **extent_properties,
                        )

            for i in range(len(self.control_points.points.geoms)):
                axe.plot(
                    self.control_points.pos[i],
                    self.control_points.z[i],
                    **control_points_properties,
                )

            for key in self.intersection_points:
                for i in range(
                    len(self.intersection_points[key].points.geoms)
                ):
                    axe.plot(
                        self.intersection_points[key].pos[i],
                        self.intersection_points[key].z[i],
                        **control_points_properties,
                    )

            axe.set_aspect("equal")
            axe.set_xlabel("Distance (m)")
            axe.set_ylabel("Altitude (m)")
            axe.set_title(self.name)
            axe.grid("both")

            fig.tight_layout()

        return fig, axe

    def plot_on_map(
        self,
        axe=None,
        figsize=None,
        point_name_dist=None,
        npts=2,
        section_properties=None,
        text_properties=None,
        control_points_properties=None,
    ):
        xx, yy = self.line.coords.xy
        if control_points_properties is None:
            control_points_properties = {
                "linestyle": "None",
                "marker": "o",
                "ms": 5,
                "color": "k",
                "mfc": "k",
            }

        if axe is None:
            fig, axe = plt.subplots(1, 1, figsize=figsize)
        # else:
        #     fig = axe.figure

        if section_properties is None:
            section_properties = dict()
        axe.plot(xx, yy, **section_properties)
        if len(self.control_points.points.geoms) > 0:
            for p in self.control_points.points.geoms:
                axe.plot(
                    p.coords[0][0], p.coords[0][1], **control_points_properties
                )

        for key in self.intersection_points:
            if len(self.intersection_points[key].points.geoms) > 0:
                for p in self.intersection_points[key].points.geoms:
                    axe.plot(
                        p.coords[0][0],
                        p.coords[0][1],
                        **control_points_properties,
                    )

        if len(self.extreme_points_names) > 0:
            if point_name_dist is None:
                point_name_dist = self.point_name_dist

            if point_name_dist is None:
                point_name_dist = self.line.length / 15

            if text_properties is None:
                text_properties = dict()

            for i in [0, 1]:
                xtmp = [xx[int((1 - 2 * npts) * i + npts - 1)], xx[-i]]
                ytmp = [yy[int((1 - 2 * npts) * i + npts - 1)], yy[-i]]

                vecx = xtmp[1] - xtmp[0]
                vecy = ytmp[1] - ytmp[0]
                nn = np.sqrt(vecx**2 + vecy**2)
                vecx = vecx / nn
                vecy = vecy / nn

                axe.text(
                    xx[-i] + point_name_dist * vecx,
                    yy[-i] + point_name_dist * vecy,
                    self.extreme_points_names[i],
                    horizontalalignment="center",
                    verticalalignment="center",
                    **text_properties,
                )
