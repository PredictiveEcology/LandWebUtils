# build_studyarea_crosswalk errors on a missing eco_field

    Code
      build_studyarea_crosswalk(fmas, eco)
    Condition
      Error in `build_studyarea_crosswalk()`:
      ! `eco_field` 'REGION_NAM' is not a column of `eco`.

# prepStudyArea rejects an unknown study area

    Code
      prepStudyArea("ZZZNonExistent", destinationPath = tempdir(), crosswalk = cw)
    Condition
      Error in `prepStudyArea()`:
      ! study area 'ZZZNonExistent' is not a grouping (see `studyAreaCrosswalk()$group`), an Alberta FMU, or a legacy study area (`LandWebStudyAreas$Name`).

