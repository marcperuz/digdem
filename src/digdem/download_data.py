# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 10:39:05 2023

@author: peruzzetto
"""

import requests

import os


def import_rainier_data(folder_out=None):
    if folder_out is None:
        folder_out = "."

    file_save = os.path.join(folder_out, "raster_rainier.tif")
    url_base = (
        "https://github.com/marcperuz/digdem/raw/"
        + "master/data/rainier_example/"
    )
    r = requests.get(url_base + "raster_rainier.tif")
    with open(file_save, "wb") as f:
        f.write(r.content)

    for fmt in ["shp", "dbf", "shx", "prj"]:
        for file in ["sections", "scar"]:
            file_save = os.path.join(folder_out, "{}.{}".format(file, fmt))
            r = requests.get(url_base + "{}.{}".format(file, fmt))
            with open(file_save, "wb") as f:
                f.write(r.content)

    for text_file in ["sections_control_points.txt", "sup_interp_points.txt"]:
        file_save = os.path.join(folder_out, text_file)
        r = requests.get(url_base + text_file)
        with open(file_save, "w") as f:
            f.write(r.text)
