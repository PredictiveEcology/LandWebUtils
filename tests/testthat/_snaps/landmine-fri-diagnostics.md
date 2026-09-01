# A RASTER ON A DIFFERENT GRID IS REFUSED, not silently recycled

    Code
      .fri_call(f, lccMap = bigger, treeClasses = 210L)
    Condition
      Error in `.fri_check_geom()`:
      ! `lccMap` is not on the same grid as `lthfc` (600 vs 400 cells). Align it first -- reading rasters of different lengths into one table recycles the shorter one and yields a complete set of plausible, wrong numbers.

---

    Code
      .fri_call(f, vegTypeMapInit = bigger)
    Condition
      Error in `.fri_check_geom()`:
      ! `vegTypeMapInit` is not on the same grid as `lthfc` (600 vs 400 cells). Align it first -- reading rasters of different lengths into one table recycles the shorter one and yields a complete set of plausible, wrong numbers.

