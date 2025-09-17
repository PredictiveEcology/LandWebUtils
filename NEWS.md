# LandWebUtils (development version)

* add package dependencies: `purrr`, `scales`, `sf`, `stats`, `terra`, `withr`;
* remove package dependencies `sp`, `raster`;
* use native R pipe and remove `magrittr` dependency;
* use `achubaty/SDMTools`, a fork of the orphaned `SDMTools` package which allows installation on R 4.5;
* add plot functions from `LandWeb_preamble` module;
* add study area creation / extraction functions from `LandWeb_preamble` module;
* use `usethis` approach to managing package imports;
* improve documentation;

# LandWebUtils 1.0.3 (2025-02-11)

* add `fs` to Imports;
* new helper `findOldSimFiles()` to find old simulation files;

# LandWebUtils 1.0.2 (2024-10-31)

* update/fix `landmine_burn1()` for "burny" scenarios;

# LandWebUtils 1.0.1 (2024-09-06)

* after further testing and discussion, will use TSF instead of SAM by default;

# LandWebUtils 1.0.0 (2023-11-06)

* use current condition TSF instead of SAM by default;

# LandWebUtils 0.1.7 (2023-09-14)

* add `ggplot2` and `rasterVis` to Imports;
* add `LandMine` diagnostic plots from module;
* improved cluster setup and use;

# LandWebUtils 0.1.6 (2023-09-11)

* require R >= 4.2 due to changes in dependency packages;
* add `LandMine` module functions;
* add helper `findSimFile()` to check previously- and currently-used simulation file names;
* fix CC TSF filename issue in `LeadingVegTypeByAgeClass()` and `LargePatches()`;
* add LandWeb-specific  functions form `LandR` package;

# LandWebUtils 0.1.5 (2022-11-08)

* require R >= 4.0 due to changes in dependency packages;
* add `future.apply` and `purrr` to Imports;
* add `rasterListByPoly()`;
* add `analysesOutputTimes()`;
* use `.ageClassCutOffs` as integer;
* improve documentation for `LeadingVegTypeByAgeClass()` and `LargePatches()`;
* add argument `crop2poly` to `LargePatches()`;
* update `LargePatches()` to handle VTM extents smaller than TSF;

# LandWebUtils 0.1.4 (2021-05-06)

* fix `LargePatches()` for `LBstatus`;

# LandWebUtils 0.1.3 (2021-03-10)

* move `simFile()` to package `SpaDES.core`;

# LandWebUtils 0.1.2 (2020-08-28)

* require R >= 3.6 due to changes in dependency packages;
* add `dplyr` to Imports;
* use `PredictiveEcology/map@development`;
* add `covr` to Suggests;
* add total area of large patches to csv outputs;
* calculate the 12.5% and 87.5% quantiles in addition to the quartiles;

# LandWebUtils 0.1.1 (2020-06-16)

* add `cleanAreaName()`;

# LandWebUtils 0.1.0 (2020-05-07)

* moved `updateSpeciesTable()` to `LandR` package;
* remove LandWeb-specific species parameter values;

# LandWebUtils 0.0.3 (2020-05-07)

* use `PredictiveEcology/map@LandWeb` and `PredictiveEcology/LandR@LandWeb`;
* restructure and fix issue in `updateSpeciesTable()`;
* use unique polygon names for leading `allCombos` in `LeadingVegTypeByAgeClass()`;
* order boxplot age classes by age instead of alphabetically;
* add argument `ext` to `simFile()` to allow files other than `.rds`;
* ensure destination path for `.csv` files exists when creating plots;
* decrease default shade tolerances for LandWeb species;
* ensure consistent species names in summary outputs;
* add `.ageClassCutOffs`;

# LandWebUtils 0.0.2 (2019-08-27)

* require R >= 3.5 due to changes in dependency packages;
* add `data.table` and `SpaDES.core` to Imports;
* add `LandR`, `map` and `pemisc` Remotes;
* add histograms for leading vegetation cover;
* add helper `simFile()` to find simulation output files;
* use more consistent species names in outputs;
* use default species growth curve parameter values;
* tweak resprouting-related parameter values;
* allow updating dispersal parameter values;
* increase `highDispersal` scenario dispersal distances;
* adjust aspen growth curve parameter value;
* add `updateSpeciesTable()` from `Biomass_borealDataPrep` module;
* write leading boxplot data to `.csv`;
* improve figures and plots;
* fix calculation and plotting bugs in `LargePatches()`;
* fix CC 'red line' calculation;

# LandWebUtils 0.0.1 (2019-03-20)

* initial version based on LandWeb-specific components from `map` package;
