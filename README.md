# digdem

- [Description](#description)
- [Installation](#installation)
- [Quick start](#quick-start)
  
## Description <a name="description"></a>

`digdem` package can be used to modify semi-automatically Digital Elevation Models (or any raster data) in specific regions, using a limited number of
control points and control profiles.

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

To be completed
