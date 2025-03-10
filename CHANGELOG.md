# CHANGELOG

## [v0.2.5](https://github.com/marcperuz/digdem/releases/tag/v0.2.5) - 2025-03-06 08:53:45+00:00

Bug fix in file reading

### Bug Fixes

- surfmod:
  - Correct control points reading ([43feb15](https://github.com/marcperuz/digdem/commit/43feb151681fbd18cf2afb32cc89c08add6c378a))

## [v0.2.4](https://github.com/marcperuz/digdem/releases/tag/v0.2.4) - 2025-02-19 07:45:08+00:00

Minor bugs correction

### Bug Fixes

- surfmod:
  - deal with ending/trailing spaces in read_control_points ([199af61](https://github.com/marcperuz/digdem/commit/199af61e1f20f106f99a681b3393b2a7767df59e))
  - change surf_min and surf_max orientation depending on indexing ([1aef6ad](https://github.com/marcperuz/digdem/commit/1aef6ad67b0239e067c1dcd428004568bebf8c36))

- gis:
  - change ndv by np.nan in rasterize_shapefile ([01727aa](https://github.com/marcperuz/digdem/commit/01727aab664f33d371f2a5041f6a95a08df25439))

## [v0.2.3](https://github.com/marcperuz/digdem/releases/tag/v0.2.3) - 2024-01-16 09:53:41

Documentation corrected/updated with support for Anaconda installation

### Documentation

- README:
  - update installation information ([3e7f16d](https://github.com/marcperuz/digdem/commit/3e7f16db2ab9fd1161a6ed4b30eed904110606cb))
  - fix installation ([7201394](https://github.com/marcperuz/digdem/commit/7201394feed9db2371859732581a80282333a5be))

## [v0.2.2](https://github.com/marcperuz/digdem/releases/tag/v0.2.2) - 2024-01-10 17:32:58

Minor bugs fix

### Bug Fixes

- plot:
  - fix issue with ndv value leading to incorrect contour plot ([3bbc9da](https://github.com/marcperuz/digdem/commit/3bbc9dad880ae8171ed5d5e150a76a092ac906c2))

### Refactor

- surfmod:
  - adapt to future releases of matplotlib ([f6094d6](https://github.com/marcperuz/digdem/commit/f6094d6d3444aca1d55a4d56a40ea0024dc7ec12))

- tests:
  - add instance of SurfMod as pytest fixture ([628bc72](https://github.com/marcperuz/digdem/commit/628bc727962a3126af5ea75bf288eee08be1ee03))

- examples:
  - delete examples/synthetic_example.py ([73f3beb](https://github.com/marcperuz/digdem/commit/73f3bebe2d3511d0bd41e7e4b7027d27563333e8))

## [v0.2.1](https://github.com/marcperuz/digdem/releases/tag/v0.2.1) - 2024-01-09 16:32:11

Fix dependencies requirements

### Bug Fixes

- pyproject:
  - add rasterio dependency ([af6938e](https://github.com/marcperuz/digdem/commit/af6938e9dbf9d0a8ce4554b200dbd1ab2dce40cd))

## [v0.2.0](https://github.com/marcperuz/digdem/releases/tag/v0.2.0) - 2024-01-09 16:02:06

New features and examples to work with raster and input text files

## [v0.1.4.rc1](https://github.com/marcperuz/digdem/releases/tag/v0.1.4.rc1) - 2024-01-09 15:57:50

New features with examples with input rasters and text files

### Feature

- examples:
  - add example scripts for Mount Rainier ([ba07db7](https://github.com/marcperuz/digdem/commit/ba07db7fcef5ef6062f64939f64d98e33e54f03b))

- surfmod:
  - change default behaviour of save_to_raster ([70d4388](https://github.com/marcperuz/digdem/commit/70d43884bdb71f6b699dc1e6d0587c2d093ef868))
  - add feature to read InterpPoints from text file ([07c73ce](https://github.com/marcperuz/digdem/commit/07c73ce3c2ee0d88fff028d9d501f2ee279504a3))
  - add method to save surf to raster ([f286a53](https://github.com/marcperuz/digdem/commit/f286a538819fc2468510b3fb7ebfda2258bc7e23))

- data:
  - update unput files for rainier example ([f530174](https://github.com/marcperuz/digdem/commit/f5301744702f98433d5e2d2ea7fdc34c199b2fcc))
  - update function to download rainier data ([d4b1700](https://github.com/marcperuz/digdem/commit/d4b17000ff5ee13bbeb893372421e97a308e9944))
  - update data for rainier example ([a630c2c](https://github.com/marcperuz/digdem/commit/a630c2cbb37f7e7445d1312a25553d798b152fc8))
  - add function to dowload sample data ([1e19e4a](https://github.com/marcperuz/digdem/commit/1e19e4a380f996ca42de92fb9ab8312b52ecd16f))
  - add data for example with real topo ([99f824b](https://github.com/marcperuz/digdem/commit/99f824b4b0f60f097cd0e40fc092c322ef10d331))

- gis:
  - add feature for reference raster in write_tiff ([5838628](https://github.com/marcperuz/digdem/commit/58386286ad0d6cdb90ef3ccdb9f81706e4912e61))
  - add gis module ([1467ef6](https://github.com/marcperuz/digdem/commit/1467ef685695d891c5fbe41746d9542eab0ff8ea))

### Bug Fixes

- section:
  - set normalized attribute when reading sections from shapefile ([fc1d2c4](https://github.com/marcperuz/digdem/commit/fc1d2c4e3a6b662e9be29d7e5b81abd75072ede7))
  - change attribute name for distance between interp_points ([b6c0616](https://github.com/marcperuz/digdem/commit/b6c06163ea7dc6e12539e23066c9bd4cdca1441b))

- plot:
  - fix bug in kwargs for interp_points plotting ([54fbd8e](https://github.com/marcperuz/digdem/commit/54fbd8e9c625c923a01390ef78360014ce604d0b))

- topo:
  - fix raster transformation ([275e28f](https://github.com/marcperuz/digdem/commit/275e28ff97ef5e49721eafb0e8a188e9d63e6c15))

- surfmod:
  - fix bug in section creation with shapefile ([87b7260](https://github.com/marcperuz/digdem/commit/87b7260ae7f589602adca2bde53a8ca6f925e853))

- gis:
  - correct import ([30a62ab](https://github.com/marcperuz/digdem/commit/30a62ab4090dc19b7259990c14536b6479235562))

- data:
  - correct projection of shapefiles for rainier example ([9d26c85](https://github.com/marcperuz/digdem/commit/9d26c854b69a6f7e9904e4fd525e2912b6630b1e))

- sup_interp_point:
  - fix bug in sup coontrol point addition ([06d6102](https://github.com/marcperuz/digdem/commit/06d6102adce8e92393a51d3c38875f216aaec679))

### Documentation

- data:
  - add README with citation for Mount rainier topographic data ([e4f94f9](https://github.com/marcperuz/digdem/commit/e4f94f910344b5fbfe969bbc5d005d648433c6e7))

- examples:
  - move examples to a dedicated folder ([13ca8b4](https://github.com/marcperuz/digdem/commit/13ca8b4e715da3b9bee991b30ab28d02f7632c22))

## [v0.1.4](https://github.com/marcperuz/digdem/releases/tag/v0.1.4) - 2024-01-08 15:23:55

Stable release to github

## [v0.1.1b](https://github.com/marcperuz/digdem/releases/tag/v0.1.1b) - 2024-01-05 14:59:19

Test for release to pypi

### Documentation

- versioning:
  - add CHANGELOG file ([401b665](https://github.com/marcperuz/digdem/commit/401b665c7ac157c430a6bb5b78e8886d257afa19))

## [v0.1.1](https://github.com/marcperuz/digdem/releases/tag/v0.1.1) - 2024-01-05 14:45:11

First initial release with minimum working example

### Feature

- general:
  - improve various aspects to get minimal working example ([96af94b](https://github.com/marcperuz/digdem/commit/96af94b2576bd1eed3c381bf8476934a553f18e9))

- example:
  - add jupyter notebook example ([ff63564](https://github.com/marcperuz/digdem/commit/ff6356426c3fc9553e56bdc83b6dd5610b70621d))

- plot:
  - improve plot functions ([3f42277](https://github.com/marcperuz/digdem/commit/3f42277bd1ff89872dc32baeff73bcf0073af589))

- initiation:
  - add initial files ([45a2926](https://github.com/marcperuz/digdem/commit/45a2926788b13e31e13f7c88bac0a18b6f8f7ec7))

### Bug Fixes

- imports:
  - change gdal import to geopandas import ([dce561f](https://github.com/marcperuz/digdem/commit/dce561f0f362c6a7460c41fbd2a140e90f21c0ed))

- toml:
  - fix project name declaration ([aafcf40](https://github.com/marcperuz/digdem/commit/aafcf408d3a5c8d6fe5c5acdac9eedb6b7ee45b4))

### Documentation

- readme:
  - improve description in readme ([278521a](https://github.com/marcperuz/digdem/commit/278521a0db40d663ddedea9e7ac45dd271cecd70))

\* *This CHANGELOG was automatically generated by [auto-generate-changelog](https://github.com/BobAnkh/auto-generate-changelog)*
