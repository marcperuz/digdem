# digdem

- [Description](#description)
- [Installation](#installation)
- [Quick start](#quick-start)
  
## Description <a name="description"></a>

`digdem` package can be used to modify semi-automatically Digital Elevation Models (or any raster data) in specific regions, using a limited number of
control points and control profiles. It follows the following steps:

1. Create instance of `SurfMod` with initial DEM, and an array specifying where the DEM will be modified;
2. Generate controlling sections (typically longitudinal and transverse sections), and specify the new altitude of their intersections if any;
3. Add additional control points along the sections by specifying their altitude, and if needed additional control points located within the mask but not on the sections;
4. `digdem` then interpolates the new DEM within the mask by:
   - Interpolating splines along each section
   - Interpolating the new DEM with Radial Basis Functions, using points along the splines, points on the contour of the mask, and if applicable additional control points within the mask.
5. Plot the new topography

Note that `digdem` is still under development, thus only minimal documentation is available at the moment, and testing is underway.
Contributions are feedback are most welcome. 

## Installation <a name="installation"></a>

To install `digdem` from GitHub or PyPi, you'll need to have `pip` installed on your computer. `digdem` is not yet available on `conda-forge`. 

It is strongly recommended to install `digdem` in a virtual environnement dedicated to this package. This can be done with `virtualenv`
(see the documentation e.g. [here](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/)).
Create the environnement with :

```bash
python -m venv /path/to/myenv
```

and activate it, on Linux:

```bash
source /path/to/myenv/bin/activate
```

and on Windows:

```cmd.exe
\path\to\myenv\Scripts\activate
```

Alternatively, if you are more used to Anaconda :

```bash
conda create -n digdem pip
conda activate digdem
```

or equivalently with Mamba :

```bash
mamba create -n digdem pip
mamba activate digdem
```

Before installing with `pip`, make sure `pip`, `steuptools` and `wheel` are up to date

```
python -m pip install --upgrade pip setuptools wheel
```

### Latest stable realease from PyPi <a name="pypi-install"></a>

```
python -m pip install tilupy
```

### Development version on from GitHub <a name="source-install"></a>

Download the GithHub repository [here](https://github.com/marcperuz/digdem), or clone it with

```
git clone https://github.com/marcperuz/tilupy.git
```

Open a terminal in the created folder and type:

```
python -m pip install .
```

## Quick start <a name="quick-start"></a>

See jupyter notebooks in examples.
