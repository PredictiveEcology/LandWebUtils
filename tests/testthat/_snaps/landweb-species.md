# landweb_sppEquiv() rejects tables it cannot map

    Code
      landweb_sppEquiv(as.data.frame(lr_sppEquiv()))
    Condition
      Error:
      ! `sppEquiv` must be a data.table.
    Code
      landweb_sppEquiv(lr_sppEquiv()[, !"SCANFI"])
    Condition
      Error:
      ! `sppEquiv` has no `SCANFI` column to map species from.
    Code
      landweb_sppEquiv(lr_sppEquiv()[!SCANFI %in% names(landweb_species_map())])
    Condition
      Error:
      ! `sppEquiv` is empty: no tree species reached this stage. LandWeb study areas are forested, so this indicates a broken species mapping (e.g. `sppEquiv`/`sppEquivCol`), which the upstream modules would otherwise treat as a valid no-species run.

# landweb_require_species() stops on an empty species input

    Code
      landweb_require_species(data.table::data.table(species = character(0)),
      "cohortData")
    Condition
      Error:
      ! `cohortData` is empty: no tree species reached this stage. LandWeb study areas are forested, so this indicates a broken species mapping (e.g. `sppEquiv`/`sppEquivCol`), which the upstream modules would otherwise treat as a valid no-species run.
    Code
      landweb_require_species(NULL, "sppEquiv")
    Condition
      Error:
      ! `sppEquiv` is empty: no tree species reached this stage. LandWeb study areas are forested, so this indicates a broken species mapping (e.g. `sppEquiv`/`sppEquivCol`), which the upstream modules would otherwise treat as a valid no-species run.
    Code
      landweb_require_species(list(speciesLayers = NULL), "speciesLayers")
    Condition
      Error:
      ! `speciesLayers` is empty: no tree species reached this stage. LandWeb study areas are forested, so this indicates a broken species mapping (e.g. `sppEquiv`/`sppEquivCol`), which the upstream modules would otherwise treat as a valid no-species run.
    Code
      landweb_require_species(list(other = 1), "speciesLayers")
    Condition
      Error:
      ! `speciesLayers` is empty: no tree species reached this stage. LandWeb study areas are forested, so this indicates a broken species mapping (e.g. `sppEquiv`/`sppEquivCol`), which the upstream modules would otherwise treat as a valid no-species run.

