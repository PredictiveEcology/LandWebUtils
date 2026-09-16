# .cleanLandWebStudyArea() requires a fire-return-interval column

    Code
      .cleanLandWebStudyArea(lthfc_sf(1, col = "other"))
    Condition
      Error in `.cleanLandWebStudyArea()`:
      ! any(c("LTHFC", "LTHRC") %in% names(poly)) is not TRUE

# polygonClean() dispatches on type and rejects what it does not know

    Code
      polygonClean(x)
    Condition
      Error in `polygonClean()`:
      ! Either fn or type must be specified

---

    Code
      polygonClean(x, type = "other")
    Condition
      Error in `polygonClean()`:
      ! Unknown type

