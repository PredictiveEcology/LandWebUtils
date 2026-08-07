# LandWebUtils (development version)

* **breaking:** reporting-polygon layers are now keyed by a curated `NAME_SHORT` rather than their long descriptive name, and output keys come from `refCodeFor()` (= `.slug(NAME_SHORT)`) rather than `abbreviate(name, minlength = 8)`. The old keying gave no cross-unit uniqueness guarantee -- only 327 of the 399 crossed FMA x ANSR names abbreviated uniquely, so 68 units silently overwrote each other's parquet aggregates and figures. Existing `_aggregates/` directories and figure folders are named under the old scheme and will be rebuilt;
* `tenureShortNames()` supplies the curated, partner-facing short name for each member tenure of the v10 boundary layer, and `validateShortNames()` enforces the invariants the keying depends on (present, filesystem-safe, <= 16 characters, unique) when the tables load; `build_studyarea_crosswalk()` gains a `name_short` column and errors on an uncurated member rather than falling back to a collision-prone long name;
* `buildCrossedReportingPolygons()` restores the v2 **tenure x sub-region** reporting units that v3 was not producing: v2 clipped each sub-region layer to a single FMA, whereas v3 clips to `studyArea`, which is now an ecoregion group unioning several tenures. Each tenure is crossed with the layers marked `cross` in `reportingPolygonLayers()` (Caribou ranges, Alberta Natural Subregions, BC biogeoclimatic zones, NWT/national ecoregions, national ecozones), with its own active/passive landbase coverage(s), and -- where both exist -- with the two together (the `tenure x landbase x Caribou` triple). Crossing happens **before** the name-grouping, so two tenures' shares of the same sub-region are no longer pooled;
* reporting-layer labels come from `labelCols`, a **comma-separated priority list** coalesced left to right, replacing the single `labelCol`. The v10 layers were assembled by merging per-jurisdiction sources, and the Caribou ranges populate `RANGE_NAME` for only 21 of 74 features (Ontario and Manitoba) -- British Columbia uses `HERD_NAME`, Alberta `LOCALRANGE`, Saskatchewan `CONUNIT`/`RNGEUNIT`. Labelling on `RANGE_NAME` alone left 53 ranges unnamed, and unnamed features are dropped, so **no tenure x Caribou reporting unit existed outside ON/MB**; placeholder values (`"Not Applicable"` and similar) are also treated as missing;
* a reporting layer may now be assembled from **several sources**: rows of `reportingPolygonLayers()` sharing a `NAME_SHORT` are built independently and merged (reduced to the common `Name`), and the new `where` column keeps only the features a source is meant to contribute. This carries the **interim NWT caribou supplement**: the v10 Caribou layer has no Northwest Territories ranges (0 of its 74 features intersect either NWT tenure), so the `Northwest Territories (NWT)` polygon from the layer LandWeb used before v10 is merged in, restoring the `FMANWT FP Caribou` / `FMANWT FR Caribou` units v2 produced. A candidate replacement (GNWT *Boreal Caribou Range Planning Regions*, 6 named regions) is documented alongside the table;
* `buildCrossedReportingPolygons()` gains `members`, restricting the crossings to the study area's own tenures. v10 tenure polygons overlap, so clipping to a study-area group boundary also catches neighbouring tenures the crosswalk assigned to a *different* group (e.g. `NorthernAlbertaUplands`, whose sole member is `FortNelson_TSA`, clips in `Tolko_Norbord_LC`); crossing those split one partner's numbers across two study-area reports;
* a crossing that fails now warns immediately as `REPORTING UNIT MISSING: crossing '<unit>' failed`, naming the unit that was lost rather than reporting a bare geometry error;
* `joinReportingPolygons()` is ported off the retired `sp` (returned a `Spatial*` object) to `sf`/`SpatVector`, returning the same class it was given; it repairs geometries with `spatialutils::repair_geoms()` instead of `reproducible::fixErrors()`, since `sf::st_make_valid()` collapses the reversed ring-winding-order polygons in the v10 FMA layer. Its `Name.1`/`Name.2` re-concatenation and `ACTIVE|PASSIVE` de-duplication are preserved (both are required by the triple crossings). The intersection now runs through `terra::intersect()` rather than `sf::st_intersection()`, which threw GEOS `TopologyException: side location conflict` on some already-repaired tenure x Caribou pairs (Dawson Creek TSA lost its Caribou units to this, despite 17,405 km^2 of genuine overlap across 9 ranges); one engine also keeps the result deterministic rather than dependent on which pairs happen to trip GEOS;
* the tenure reporting layer is now labelled from `.fmaMemberIdentity()` (the coalesced per-jurisdiction identity) rather than `FMA_NAME` alone, which is `NA` for every non-Alberta member and so left the tenure layer empty for the BC/SK/ON/MB study-area groups;
* `buildCrossedReportingPolygons()` skips tenures contributing less than `min_area_km2` (default 1 km^2) to the study area: a neighbouring tenure that merely grazes the study-area boundary leaves a sliver in the clipped layer, and crossing it would mint a full set of reporting units -- each with its own `refCode`, parquet aggregate directory, and figures -- over a handful of pixels;
* `buildLandbasePolygons()` gates on the **tenures present in the study area** (`tenures=`) instead of a regex matched against the study-area name (`applies_to`/`studyAreaName=`). The old gating silently stopped matching anything once study areas became ecoregion groups (`"WesternAlbertaUpland"` matches no company name), so no landbase was ever fetched;
* generated download directory and file names now use underscores instead of dots (e.g. `National_Ecozones`, not `National.Ecozones`), which is friendlier on Windows; previously cached downloads under the old dotted names will be re-fetched once, or can be renamed in place;
* `buildReportingPolygons()` now locates extracted vector files recursively, fixing layers whose shapefile unpacks into a subdirectory (the national ecoregion and ecozone archives) being dropped as though they did not intersect the study area; a layer whose file cannot be found now warns instead of being dropped silently;
* `prepEcoregionLayer()` and `prepEcoprovinceLayer()` likewise search recursively, and now fail with an informative error rather than a subscript-out-of-bounds when no shapefile is found;
* drop support for R 4.2 and 4.3 due to changes in dependency packages;
* add package dependencies: `purrr`, `scales`, `sf`, `stats`, `terra`, `withr`;
* remove package dependencies `sp`, `raster`;
* use native R pipe and remove `magrittr` dependency;
* replace use of orphaned `SDMTools` package with `landscapemetrics`;
* use `parallelly` to setup cluster in `landmine_optim_clusterSetup()`;
* add plot functions from `LandWeb_preamble` module;
* add study area creation / extraction functions from `LandWeb_preamble` module;
* use `usethis` approach to managing package imports;
* improve documentation;
* `buildReportingPolygons()` and `reportingPolygonLayers()` assemble the reporting-polygon list using `spatialutils::prep_vector()` and the `workflowtools` `*_once` download helpers, replacing the per-study-area hardcoding;
* `build_studyarea_crosswalk()` (with `studyAreaCrosswalk()`) groups the LandWeb FMA/TSA/FML boundary polygons into ecologically-coherent study areas by assigning each tenure to its dominant ecological unit (an ecoregion by default; supply any eco\* layer for a different granularity) via verified spatial overlap, emitting a crosswalk with `company` and `province` columns so licensee-driven reruns still map cleanly; reversed-winding polygons are repaired via `spatialutils::repair_geoms()` rather than `sf::st_make_valid()` (which collapses them);
* `extractFMA()` falls back to matching `LandWebStudyAreas$Description` against an `FMA_NAME` column when the canonical short names are absent (the updated FMA boundaries layer);
* `prepEcoregionLayer()` and `prepEcoprovinceLayer()` download and prepare the AAFC ecostrat ecoregion/ecoprovince layers for use as the study-area grouping layer (named with a `Layer` suffix so they do not shadow `LandR::prepEcoregions()`, which `Biomass_borealDataPrep` calls);
* `prepStudyArea()` now resolves an ecoregion study-area group to the union of its member FMAs via the grouping crosswalk, while still handling Alberta FMU names, the `random` test area, and legacy `LandWebStudyAreas` names;
* drop the `LandR` dependency by vendoring the `equivalentName()` name-lookup helpers;
* add package dependencies: `googledrive`, `spatialutils`, `workflowtools`;

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
