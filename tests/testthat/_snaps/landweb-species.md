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

