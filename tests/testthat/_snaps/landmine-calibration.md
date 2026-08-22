# landmine_optim_unpack rejects malformed input

    Code
      landmine_optim_unpack(c(1, 2, 3))
    Condition
      Error:
      ! `par` must have 7 elements (4 spawnNewActive, 2 sizeCutoffs, 1 spreadProb), not 3.

---

    Code
      landmine_optim_unpack(c(-0.22, NA, -0.65, -1.48, 3.21, 4.72, 0.856))
    Condition
      Error:
      ! `par` must not contain NA.

# landmine_optim_params_read rejects a bad rowID

    Code
      landmine_optim_params_read(f, rowID = 99L)
    Condition
      Error:
      ! `rowID` must be a single row number between 1 and 1.

---

    Code
      landmine_optim_params_read("no-such-file.csv")
    Condition
      Error:
      ! calibration parameter file not found: no-such-file.csv

