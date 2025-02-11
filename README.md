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

`digdem` can be installed from GitHub, [PyPi](https://pypi.org/project/digdem/) or [Anaconda](https://anaconda.org/conda-forge/digdem). Supposedly stable releases are distributed to PyPi and Anaconda. Distribution to Anaconda is made manually from the PyPi package, thus the last available version on conda-forge may be an older version than the one distributed on PyPi. 

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
### Latest stable release from Anaconda

With Anaconda:

```bash
conda install digdem
```

or equivalently with Mamba: 

```bash
mamba install digdem
```

### Latest stable realease from PyPi <a name="pypi-install"></a>

Before installing with `pip`, make sure `pip`, `steuptools` and `wheel` are up to date

```
python -m pip install --upgrade pip setuptools wheel
```

Then run

```
python -m pip install digdem
```

### Development version on from GitHub <a name="source-install"></a>

Download the GithHub repository [here](https://github.com/marcperuz/digdem), or clone it with

```
git clone https://github.com/marcperuz/digdem.git
```

Open a terminal in the created folder and type:

```
python -m pip install .
```

If you want to developp and test `digdem` and have your modification directly taken into account when running your code, use:

```
python -m pip install -e .
```

Alternatively, if you don't want to install the package with pip, you can also add the `src` folder to your Python path.

## Quick start <a name="quick-start"></a>

See jupyter notebooks in examples.
