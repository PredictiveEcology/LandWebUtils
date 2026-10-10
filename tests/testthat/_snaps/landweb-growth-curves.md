# growth traits keep the five traits and their source, and refuse missing ones

    Code
      landweb_growth_traits(bad)
    Condition
      Error:
      ! growth-curve traits are missing for: Pinu_con

---

    Code
      landweb_growth_traits(growth_species_table()[, !"longevity"])
    Condition
      Error:
      ! the species table has no longevity

# a unit without a dominant species, or whose dominant has no traits, stops

    Code
      landweb_unit_growth_traits(sp, tr, data.table::data.table(unit = "Pice_gla",
        dominant = "Pice_gla"))
    Condition
      Error:
      ! no dominant species for unit(s): Abie_spp

---

    Code
      landweb_unit_growth_traits(sp, tr, units)
    Condition
      Error:
      ! no growth-curve traits for: Abie_bal

# damage agent codes: bark beetles and defoliators, by plot source

    Code
      landweb_damage_codes("fire")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "barkBeetles", "defoliators"

# a bootstrap stops when no plot lies in the fitting area

    Code
      landweb_resample_psp(psp_fixture(), far, seed = 1L)
    Condition
      Error:
      ! no plot lies inside the fitting area

# trait frequencies count each species' trait sets over the refits

    Code
      landweb_growth_trait_frequency(refits[, !"resample"], tr)
    Condition
      Error:
      ! the refits' traits have no resample

