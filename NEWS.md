# LandWebUtils 1.0.3.9031

* `landmine_fri_summary()` gains `studyArea`, so the fix for the denominator lives at the source rather than only in `landmine_fri_metrics()`. Callers that do not pass it keep their previous (wrong) numbers until updated deliberately;
* `landmine_plot_fri_zones()` now builds its two panels with `patchwork` and gives each its own fill scale. Faceting forced a single shared scale, on which log-ratios spanning about one unit collapsed entirely against target intervals spanning 170 years -- panel B rendered as a flat block of colour. Panel A is also now a sequential single-hue ramp (it encodes a magnitude) rather than a diverging one, and the caption wraps instead of running off the page;
* `landmine_plot_fri_drivers()` uses a reserved status palette -- neutral for on-target, amber for outside tolerance, red for severe -- rather than an arbitrary hue cycle that coloured "ok" red;

# LandWebUtils 1.0.3.9030

* new `landmine_fri_metrics()`, `landmine_fri_verdict()`, `landmine_plot_fri_zones()` and `landmine_plot_fri_drivers()`: per-FRI-zone attainment diagnostics, so a user does not have to read a simulation log to learn that a zone missed its target, nor mount a bespoke investigation to find out why;
* **these immediately found that there was no attainment problem at all.** `landmine_fri_summary()` divides burns by every flammable pixel carrying a fire return interval, but LandMine masks both its rate-of-spread map and its spread probability to `sim$studyArea`, and the interval raster is not similarly masked. Pixels outside the polygon cannot be ignited or spread into, yet each one inflates its zone's achieved interval. On WesternAlbertaUpland 988,437 of 3,361,240 flammable zone pixels (29.4%) were outside, burning at a mean rate of 0.51 against 16.59 inside; zone 170 had **73 burnable pixels out of 146,627**. Uncorrected, two zones appeared to under-burn by 3.5x and 12.2x. With the denominator restricted to the burnable area, **every zone falls between 0.90 and 1.07** -- the fire model has been hitting its targets throughout;
* `landmine_fri_metrics()` therefore takes a `studyArea` and reports both figures, `ratio` over the burnable area and `ratioUnmasked` the way the module's own summary computes it, plus `pctInStudyArea`, so the discrepancy is visible rather than mistaken for a fire-model failure. `landmine_fri_verdict()` names it explicitly in its one-line summary;
* the structural columns (`pctFlam`, `nPatches`, `largestPatchPct`, `pctNoCohorts`, `pctTreedLCC`, `cohortGap`, `pctCohortsLost`) exist to tell the remaining failure modes apart: porosity from fragmentation, genuinely non-forest land from forest that was never populated with cohorts, and land that never had cohorts from land that lost them;
* a zone that is entirely non-flammable is now **reported** rather than silently dropped, as `landmine_fri_summary()` does by masking before it enumerates zones;
* every raster argument is checked against one grid. Reading rasters of different lengths into one table recycles the shorter one and yields a complete set of plausible, wrong numbers -- and this is a live hazard here, not a hypothetical one: the dataPrep land-cover and stand-age rasters are written on the larger `biomassParam` extent while the rest of the stage's outputs are on the simulation grid;

# LandWebUtils 1.0.3.9029

* new `landmine_fire_ros()` and `landmine_known_species()`, promoted from `LandMine`'s `fireROS()` -- the last and largest of the module's untested internals, and the surviving hypothesis for why FRI zones 55 and 170 under-burn. It was previously reachable only through a simulation, because it read eight values off `sim`/`mod`/`P(sim)`; it is now a pure function of rasters and tables, so per-zone rate-of-spread distributions can be compared in seconds rather than hours;
* correctness is established by **equivalence testing against the module's inline implementation** over 200 randomised landscapes (varying species mixes, stand ages including `NA`, flammability including `NA`, and both `ROStype` values), plus a sweep in which each fuel type in turn is absent from the landscape;
* documents, and pins by test, two behaviours that are almost certainly bugs but are **deliberately not fixed** here, because either would change the fire regime and both sit in the code path currently under investigation:
    - **mixedwood stands never receive the `ROSTable`'s mixed rates.** `mixed` is the one fuel type invented for attribute-table entries matching no species, rather than derived from `landmine_known_species()`, so the species join leaves its `leading` value `NA` and the subsequent join onto `ROSTable` drops it. Mixedwood therefore falls through to `ROSother` -- 30 rather than 12 or 17 with the default table, i.e. the mature-spruce rate at every age -- and both `mixed` rows of the table are dead entries;
    - **the young age class's species filter is applied under the opposite condition to its two siblings.** The git history shows this is an unfinished copy-paste (introduced as a copy of the mature guard that still named `"mature"` in its condition while assigning to the young class; the follow-up commit repaired the subset key but not the dropped `!`). The new `youngGuard` argument selects `"legacy"` or `"symmetric"` so the difference can be measured rather than argued about;
* the stand-age class boundaries are now the `ageCutoffs` argument instead of `40` and `120` written three times. They are the only place the `ROSTable`'s age *labels* acquire a numeric meaning, since the table itself carries no ages;
* the runtime self-test is kept (as `assertions`, still defaulting to `getOption("LandR.assertions", TRUE)`) so development-mode runs are still checked mid-simulation, and the two classes of pixel it silently skipped -- `tsf == 0` and `tsf` past its 999-year bin, both dropped by a right-closed `cut()` -- are covered by tests instead;
* inputs are now validated rather than producing a plausible-looking map: a vegetation map whose attribute table was dropped (`terra::writeRaster()` stores it in an `.aux.xml` sidecar), an attribute table that is not `1..n` (the lookup uses table row *positions* as raster *values*), rasters on different grids, and a `ROSTable` missing the rows the `ROSother` and `"burny"` paths need;
* `ROSother`'s second assertion no longer compares a value against 5% of itself;

# LandWebUtils 1.0.3.9028

* new `landmine_reburn_budget()`, promoted from both branches of `LandMine`'s reburn loop -- the densest concentration of positional contracts in that module, and the one its own comments warn is order-sensitive. It rebuilds the per-zone fire counts and target sizes for the next round, for phase 1 (each too-small fire keeps its full original target) and phase 2 (new fires sized to the remaining shortfall);
* three contracts are now tested, none of which errored when broken -- they simply produced a plausible but wrong fire regime: the right join is driven by `friByPolygon`, so `numFiresThisPeriod` comes out zones-ascending with the `NA` row last, matching the `[.GRP]` indexing at the call site; `na.omit()` on the whole table is load-bearing, dropping both zero-shortfall zones and the `NA`-FRI row so `fireSizesInPixels` stays positionally paired with the start cells; and `remainingSize` is assigned by position onto the too-small rows;
* correctness is established by **equivalence testing against the original inline logic** over 200 randomised multi-zone inputs (interleaved zones, both phases), rather than against hand-computed expectations -- the failure mode is subtle enough that hand-derived cases would not be convincing;
* the input table is no longer modified by reference, which the inline `:=` version did;

# LandWebUtils 1.0.3.9027

* new `landmine_ignition_budget()`, promoted from `LandMine`'s `Init()`. It masks zero-FRI and non-flammable pixels out of the fire-return-interval raster, tabulates pixels per FRI zone, and converts that to expected fires per year (`(area / meanFireSize) / FRI`). Verified behaviour-neutral against the previous inline block on the real WesternAlbertaUpland rasters: the masked raster, the per-zone FRI vector and `numFiresPerYear` (including names and order) are all identical;
* the promotion exists to pin an **ordering contract that is consumed positionally downstream**. The reburn loop indexes per-zone fire counts with `numFiresThisPeriod[.GRP]`, which is valid only because (1) `terra::freq()` returns rows ascending by value, (2) the `value = NA` row is appended last, and (3) dividing by the `NA` return interval makes that entry's *value* `NA` so a later `na.omit()` drops it -- an `NA` *name* alone would not. Break any one and zones silently receive each other's fire counts, with no error. All three are now tested at the point where the contract is created, including with a deliberately unsorted input raster;
* also dropped dead code in the same block: a `value = seq_len(NROW(.))` column that was then `order()`ed by, a trivially-identity permutation since the column *is* the row index;

# LandWebUtils 1.0.3.9026

* new `landmine_fri_summary()`, `landmine_area_burned_by_zone()`, `landmine_draw_num_fires()`, `landmine_sizes_to_pixels()` and `landmine_reburn_ceiling()`, promoted from the `LandMine` module so they can be unit-tested. The module has no working test suite -- its `tests/` are the untouched SpaDES stub, with hardcoded `c:/Eliot/...` paths and a call to a function that does not exist -- so this logic was previously exercised only by running a 1000-year simulation, where each failure mode below produces plausible numbers rather than an error;
* `landmine_area_burned_by_zone()` **fixes a bug**: burned pixels were counted with `unname(table(currBurn[ids])[2])`, positional indexing into a `table()` that assumes the levels are always `{0, 1}`. A zone that burns *completely* has only one level, so `[2]` is `NA`, and the module's following `npix[is.na(npix)] <- 0` -- whose comment says it handles the *no-burn* case -- silently turned a total burn into **zero hectares burned**. Now counted with `sum(vals == 1, na.rm = TRUE)`, with tests for the fully-burned, partly-burned and unburned cases;
* `landmine_fri_summary()` **removes a verbatim duplication**: the same 20 lines appeared in both the single- and multi-mode summary functions, so the two could silently diverge. Two edge cases are deliberately preserved rather than "fixed", and pinned by tests so any future change is deliberate: a zone that never burns yields `Inf`, and the un-`na.rm`'d `sum()` makes a zone `NA` if a burn pixel is `NA` where `lthfc` is not (an implicit contract that the two `NA` masks agree);
* `landmine_draw_num_fires()` restores the zone names that `rnbinom()` drops -- the counts are indexed by zone downstream, so losing them would assign zones each other's fire counts;
* `landmine_sizes_to_pixels()` is tested for the property that justifies it existing at all: probabilistic rounding preserves `E[pixels]`, so sub-pixel fires are not systematically lost;

# LandWebUtils (development version)

* New `landmine_plot_calibConvergence()` plots DEoptim's best value against iteration and annotates where improvement stopped, so it is visible whether a run converged or was still improving when it hit `itermax`. On the first 240 m/120 m calibrations the last improvement came at iterations 149 and 128 of 200, with no gain in the final 50 -- i.e. iteration count was not the binding constraint, objective noise was;

* `landmine_optim_fitAndison()` pins the RNG kind when `crnSeed` is given, and restores the caller's generator on exit. `set.seed()` reseeds whichever generator is current, and cluster workers run L'Ecuyer-CMRG (from `clusterSetRNGStream()`) while a plain session runs Mersenne-Twister -- so the same parameters and the same `crnSeed` gave objective values differing by ~8x (0.260 vs 2.117) depending on where they were evaluated, and a calibration could not be reproduced outside the cluster that produced it;

* New `landmine_optim_islandWeight()`, used by default to scale the remnant-island term by resolution. Andison rasterised at **50 m** (quarter-hectare) pixels, and his Table 3.6 puts 65% of island area in the 2-5 ha class for fires under 1,000 ha. At 240 m a single pixel is already 5.76 ha, so that entire class is unreachable by construction and no parameter choice can meet the target. The weight is `min(1, 50 / pixelSize)`: 1.0 at 50 m, 0.42 at 120 m, 0.21 at 240 m. The weights used are recorded on the result's `"weights"` attribute;

* `landmine_optim_fireSizes()` converts fire sizes from **hectares** to pixels for a given `pixelSize`, and `landmine_optim_calibrate()` now takes `fireSizesHa`. The calibration targets are defined in hectares, but the model takes pixels, so a fixed pixel-count vector meant a different physical fire at every resolution: `c(10, ..., 1e5)` px is 10-100,000 ha at 100 m but 57.6-576,000 ha at 240 m. The v3 move to 240 m therefore put nearly the whole calibration outside the 10-3,000 ha range Andison fitted, and made calibrating at several pixel sizes uninterpretable -- the runs were fitting different physical fires, so agreement could not demonstrate that parameters transfer across resolutions. Warns for areas outside the fitted range, or under 5 pixels;

* `landmine_optim_clusterSetup()` no longer requires `nnodes >= 70` and a multiple of 7. Those are constraints on DEoptim's population size, not the worker count; conflating them made the calibration unrunnable on a 64-core host and tied both `NP` and the per-worker RNG streams to `availableCores()`, i.e. to the machine and its load. It also gains a `seed` argument, since it previously drew its RNG stream seed from ambient state. `landmine_optim_calibrate()` defaults `NP` to 70 independently of the worker count;

* `landmine_optim_params_append()` records the `objective` a fit was made against, because parameter values are only comparable within one objective;

* **BREAKING (fire behaviour):** `landmine_burn1()`'s second spawn rule now tests `numActive >= nActiveCells1[2]` rather than `>`. With `>`, the `numActive == nActiveCells1[2]` case matched none of the three rules and fell through to `spawnNewActive[1]` -- the value intended for a barely-active fire -- so the *highest* spawn probability was applied to one of the most active states. A further gap is documented but deliberately not patched: when `numActive >= nActiveCells1[2]` and `size >= sizeCutoffs[1]` no rule applies, and choosing a value there needs a re-fit. These thresholds are LandWeb-derived, not Andison's: the thesis has no firelet-parameter table and states that its calibration was manual;

* **BREAKING (calibration):** the "did the fire reach its target size" penalty in `landmine_optim_fitSN()` and `landmine_optim_fitSN2()` was dead code. `burnedMap` stores `initialPixels` (the ignition cell index) at each burned cell, so `sum(burnedMap[], na.rm = TRUE)` summed pixel IDs (~5e9) rather than counting cells, and the comparison against a target of at most 1e5 was never true. It now uses the new `nBurned` column. Parameter rows in `LandMine_DEoptim_params.csv` fitted before this are **not comparable** to rows fitted after;

* New `landmine_optim_fitAndison()` fits against the statistics Andison (1996) actually reported, replacing the single constant perimeter-to-area target: the shape-vs-area relationship `SHAPE = 1.770 + 0.041*log10(AREA)^4` (§3.3.1), a continuous size shortfall rather than a step penalty, and the remnant-island fractions of Table 3.5. Supporting functions: `landmine_optim_shapeIndex()`, `landmine_optim_shapeTarget()`, `landmine_optim_islands()`, `landmine_optim_islandTarget()`. It accepts `crnSeed` for common random numbers. Note the defaults use the **empirical** column of each Andison table (shape intercept 1.770, island fractions 4.0/6.0/9.0); the other column in each is that model's own output, and targeting it would be circular;

* `landmine_optim_calibrate()` gains `objective` (`"andison"`, the new default, or `"sn"` to reproduce historical fits), `crnSeed`, and `replicates`;

* `landmine_optim_burnFun()` returns `nBurned` and `shapeIndex` in its `LM` table, and passes `spreadProbRel` to `spread2()` as a numeric vector -- bit-identical, and about twice as fast per objective evaluation;

* New `landmine_optim_unpack()`, `landmine_optim_params_read()`, `landmine_optim_params_append()`, `landmine_optim_landscape()` and `landmine_optim_calibrate()` move the LandMine fire-spread calibration out of `LandMine.Rmd`, which is knitted with `eval = FALSE` and had drifted: it referenced the removed `raster` package API, called `terra::cellFromRowCol()` with `raster`'s argument names, and used an invalid `read.csv(file, row.names = FALSE)`. The `10^` parameter convention had been repeated at eight call sites and the CSV schema at four, two of which disagreed (one appended a row, the other overwrote the file). Behaviour is unchanged -- the existing calibration rows remain valid and comparable;

* `landmine_burn1()` accepts numeric vectors for `spreadProb` and `spreadProbRel` and forwards them to `spread2()`. Prefer these to rasters: `spread2()` re-materialises a raster in full on every iteration, so a raster costs one O(ncell) read per spread step, and `landmine_burn1()` steps one iteration at a time. Numeric `spreadProbRel` requires SpaDES.tools >= 2.1.2.9000 (PredictiveEcology/SpaDES.tools#106), now declared in `Imports`;

* `reportingPolygonLayers()` no longer includes **`Parks`**. Reporting units are `studyAreaReporting` intersected by tenure, and protected areas are carved *out* of forest tenures, so a parks layer has almost nothing to report against: for `WesternAlbertaUpland`, 81 of 732 parks merely abut the FMAs and the entire intersection is 0.05 km^2 (against 6,538 km^2 for the buffered `studyArea`). v2 did not report parks either. Meaningful parks reporting needs a wider domain than one study area -- and wider than one province, since provincial parks sit within a province but national parks cross provincial boundaries -- so the layer is excluded pending that review;

* **fix:** `.caribouNameFixes()` had the `Snake-Sahtaneh` entry **backwards**, rewriting British Columbia's *correct* spelling into ECCC's typo. CGNDB has **Sahtaneh River** (BC) and no match for "Sahtahneh", so the federal layer is the one in error; the entry is removed and BC's name kept. Alberta's `East Side Athabasca` / `West Side Athabasca` are likewise left alone: both are **official enumerated domain values** in Alberta's published metadata ("Local Population Ranges" / "Population Subunit", authority Wildlife Management, Alberta Environment and Parks), so dropping "River" is Alberta nomenclature rather than a truncation — and not a field-width artifact either, since `LOCALRANGE` is 25 characters wide, its longest value is 19, and "West Side Athabasca River" is exactly 25. `Bischto` → `Bistcho` stands (CGNDB: Bistcho Lake; no such place as "Bischto");

* correct the Manitoba caribou provenance note: the layer was provided by the Government of Manitoba in **2018** and is **current** (the `_2015` in the filename is the delineation year, not the delivery date). It was previously described as a 2015 copy worth re-requesting, which understated it -- Manitoba has published nothing since, so it is the authoritative version;

* **breaking:** caribou reporting ranges are now **assembled from the six jurisdictional sources** (`caribouRangeLayers()` + `buildCaribouRanges()`) instead of fetched as one pre-combined file. `reportingPolygonLayers()` carries a single `Caribou` row with the new `source = "assembly"`, dispatched to a named assembler. This replaces both the hand-combined v10 layer — whose manual combine had **dropped the Northwest Territories ranges entirely** (0 of its 74 features intersect either NWT tenure) and left range names scattered across five per-jurisdiction columns — and the interim GNWT-only supplement that patched the NWT hole. Reporting goes to *jurisdictional* partners, who each want their own management-unit names, so the jurisdictional sources are authoritative; the ECCC national layer is retained only as a documented comparison in `scripts/make_caribou_reference.R`, **not** as a reporting layer. **Caribou reporting-unit names therefore change from v2**, which labelled its boreal units from ECCC (`Bistcho (BIS)`, `Saskatchewan Boreal Plains (BPL)`, …); v2↔v3 caribou comparisons need a crosswalk rather than matching on name;
* locally **extirpated herds are excluded** from the caribou ranges. Two sources flag it in different vocabulary — AB `STATUS == "Extirp"` (the Banff range, lost to an avalanche in 2009, and one of AB's *mountain* ranges) and BC `HERD_STATUS == "Extirpated"`. `statusCols` is a candidate list because BC's field is `HERD_STATUS` over its WFS but truncated to `HERD_STAT` in a shapefile export;
* British Columbia is now fetched **live from the BC Data Catalogue WFS** rather than a cached order-based BCGW export, whose URL was not reproducible and whose content was stale: the WFS release has 55 features against the old export's 72, and flags **5** extirpated herds where the old one flagged 1;
* caribou range names are **corrected where a source is demonstrably wrong** (`Bischto` → `Bistcho`). This matters beyond cosmetics: the summaries group by *name* rather than by feature, so a range two provinces spell differently appears as two reporting units — the same mechanism that stitches `Chinchaga`, `Narraway` and `Redrock-Prairie Creek` back into single units across the AB/BC border. The bar is **evidence of an error, not disagreement between sources**, each entry verified against NRCan's Canadian Geographical Names Database (CGNDB);
* Alberta is labelled on `SUBUNIT` before `LOCALRANGE` — **subunit granularity throughout**, decided deliberately rather than to match v2, since no single setting does. v2 reported Alberta's *mountain* ranges at subunit level (`A La Peche`, `Narraway` and `Redrock-Prairie Creek` were three separate v2 units, all subunits of the one `LOCALRANGE` `"West Central"`) but its *boreal* ranges whole, because it sourced AB boreal from ECCC at range level. `SUBUNIT` therefore matches v2 for the mountain ranges and reports finer than v2 for exactly one boreal range: `East Side Athabasca River` resolves into its 7 subunits (`Agnes`, `Algar`, `Bohn`, `Christina`, `Egg-Pony`, `Wandering`, `Wiau`). Every other Alberta range is 1:1 with both v2 and ECCC, and the finer breakdown is kept as an improvement;
* Ontario is included. It is easy to omit and its loss is silent: `RANGE_NAME` was populated for exactly the ON and MB features of the old combined layer, so leaving Ontario out cost every ON tenure × Caribou reporting unit;
* `.fetchVectorFile()` factors the three transports out of `buildReportingPolygons()` so the assembler and the layer builder fetch identically;

* **breaking:** reporting-polygon layers are now keyed by a curated `NAME_SHORT` rather than their long descriptive name, and output keys come from `refCodeFor()` (= `.slug(NAME_SHORT)`) rather than `abbreviate(name, minlength = 8)`. The old keying gave no cross-unit uniqueness guarantee -- only 327 of the 399 crossed FMA x ANSR names abbreviated uniquely, so 68 units silently overwrote each other's parquet aggregates and figures. Existing `_aggregates/` directories and figure folders are named under the old scheme and will be rebuilt;
* `tenureShortNames()` supplies the curated, partner-facing short name for each member tenure of the v10 boundary layer, and `validateShortNames()` enforces the invariants the keying depends on (present, filesystem-safe, <= 16 characters, unique) when the tables load; `build_studyarea_crosswalk()` gains a `name_short` column and errors on an uncurated member rather than falling back to a collision-prone long name;
* `buildCrossedReportingPolygons()` restores the v2 **tenure x sub-region** reporting units that v3 was not producing: v2 clipped each sub-region layer to a single FMA, whereas v3 clips to `studyArea`, which is now an ecoregion group unioning several tenures. Each tenure is crossed with the layers marked `cross` in `reportingPolygonLayers()` (Caribou ranges, Alberta Natural Subregions, BC biogeoclimatic zones, NWT/national ecoregions, national ecozones), with its own active/passive landbase coverage(s), and -- where both exist -- with the two together (the `tenure x landbase x Caribou` triple). Crossing happens **before** the name-grouping, so two tenures' shares of the same sub-region are no longer pooled;
* reporting-layer labels come from `labelCols`, a **comma-separated priority list** coalesced left to right, replacing the single `labelCol`. The v10 layers were assembled by merging per-jurisdiction sources, and the Caribou ranges populate `RANGE_NAME` for only 21 of 74 features (Ontario and Manitoba) -- British Columbia uses `HERD_NAME`, Alberta `LOCALRANGE`, Saskatchewan `CONUNIT`/`RNGEUNIT`. Labelling on `RANGE_NAME` alone left 53 ranges unnamed, and unnamed features are dropped, so **no tenure x Caribou reporting unit existed outside ON/MB**; placeholder values (`"Not Applicable"` and similar) are also treated as missing;
* a reporting layer may now be assembled from **several sources**: rows of `reportingPolygonLayers()` sharing a `NAME_SHORT` are built independently and merged (reduced to the common `Name`), and the new `where` column keeps only the features a source is meant to contribute. This carries the **NWT caribou supplement**: the v10 Caribou layer has no Northwest Territories ranges (0 of its 74 features intersect either NWT tenure), which had cost the `FMANWT FP Caribou` / `FMANWT FR Caribou` units v2 produced. The gap is filled from GNWT *Boreal Caribou Range Planning Regions* (MapServer layer 97), which subdivides NT1 into 6 named planning regions the way `LOCALRANGE`/`HERD_NAME` subdivide AB/BC -- so an NWT tenure gets a real sub-region breakdown rather than the single whole-territory polygon an earlier stopgap supplied. Labelled on `REGION` rather than `NAME`: both exist, but `NAME` carries diacritics and apostrophes (`Sahtú`, `Wek'èezhìi`, `Gwich'in`) and this label is slugged into the `refCode` that keys parquet aggregates and figure filenames. No `where` filter is needed -- the layer is NWT-only by construction, and its `Yukon` region simply does not intersect a study area that stops short of Yukon;
* `reportingPolygonLayers()` gains a third `source` type, **`"geojson"`**, for a URL that returns GeoJSON directly (an ArcGIS REST `query` endpoint) rather than a zipped shapefile: it is fetched with no archive-extraction step, and named from the layer key since `basename()` of a query URL is not a usable filename. Unknown `source` values now error when the table loads -- previously anything that was not `"drive"` fell through to the zipped-URL path, so a typo surfaced later as an opaque "no vector file found";
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
