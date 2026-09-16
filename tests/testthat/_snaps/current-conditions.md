# cc_age_composite() rejects an unusable fill year

    Code
      cc_age_composite(cc_grid(1), cc_poly(), cc_grid(1))
    Condition
      Error:
      ! `fillYear` is required when `fill` is supplied.

---

    Code
      cc_age_composite(cc_grid(1), cc_poly(), cc_grid(1), fillYear = 2030)
    Condition
      Error:
      ! `fillYear` (2030) is after `epoch` (2025).

# cc_age_pct_missing() errors when there is no forest to measure

    Code
      cc_age_pct_missing(cc_grid(80), cc_poly(), cc_grid(18L))
    Condition
      Error:
      ! No forest (classes 1, 2, 5, 6) in the study area, so the share of forest missing an age is undefined.

