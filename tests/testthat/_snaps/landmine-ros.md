# ageCutoffs must be length 2 and increasing

    Code
      .ros_call(f, ageCutoffs = 40)
    Condition
      Error in `landmine_fire_ros()`:
      ! length(ageCutoffs) == 2L is not TRUE

---

    Code
      .ros_call(f, ageCutoffs = c(120, 40))
    Condition
      Error in `landmine_fire_ros()`:
      ! ageCutoffs[[1]] < ageCutoffs[[2]] is not TRUE

# a vegetation map with no attribute table is a clear error, not a silent NA map

    Code
      .ros_call(f)
    Condition
      Error in `landmine_fire_ros()`:
      ! `vegTypeMap` has no raster attribute table, so its values cannot be mapped to species. `terra::writeRaster()` stores the attribute table in an `.aux.xml` sidecar -- copying the `.tif` alone silently drops it.

# an attribute table with duplicate or missing values is refused

    Code
      LandWebUtils:::.ros_check_rat(data.table::data.table(ID = c(1L, 1L, 3L),
      Species = c("a", "b", "c")))
    Condition
      Error in `LandWebUtils:::.ros_check_rat()`:
      ! `vegTypeMap`'s attribute table has duplicate values, so a raster value cannot be mapped to one vegetation type.

---

    Code
      LandWebUtils:::.ros_check_rat(data.table::data.table(ID = c(1L, NA, 3L),
      Species = c("a", "b", "c")))
    Condition
      Error in `LandWebUtils:::.ros_check_rat()`:
      ! `vegTypeMap`'s attribute table has missing values, so a raster value cannot be mapped to one vegetation type.

# rasters on a different grid are refused

    Code
      .ros_call(f)
    Condition
      Error in `.as_cell_vector()`:
      ! `rstTimeSinceFire` has 19 values but `vegTypeMap` has 20 cells; they must be on the same grid.

# ROSother is validated against the table's range and against mature spruce

    Code
      .ros_call(f, ROSother = 500L)
    Condition
      Error in `landmine_fire_ros()`:
      ! `ROSother` (500) is outside the range of `ROSTable$ros` (6-30).

---

    Code
      .ros_call(f, ROSother = 20L)
    Condition
      Error in `landmine_fire_ros()`:
      ! `ROSother` (20) is not within 5% of the mature spruce rate of spread (30).

# a table with no mature spruce row cannot validate ROSother

    Code
      .ros_call(f)
    Condition
      Error in `landmine_fire_ros()`:
      ! `ROSTable` must have exactly one mature-spruce row to validate `ROSother` against; found 0.

# 'burny' needs an immature_young deciduous row

    Code
      .ros_call(f, ROStype = "burny")
    Condition
      Error in `landmine_fire_ros()`:
      ! `ROStype = "burny"` needs exactly one immature_young deciduous row in `ROSTable` to use for non-flammable pixels; found 0.

