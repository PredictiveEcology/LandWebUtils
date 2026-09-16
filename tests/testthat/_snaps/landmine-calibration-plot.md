# landmine_plot_calibDiscrimination requires the expected columns

    Code
      landmine_plot_calibDiscrimination(data.frame(a = 1))
    Condition
      Error:
      ! `scores` needs column(s): group, objective. Got: a.

# landmine_plot_calibConvergence errors informatively on the wrong input

    Code
      landmine_plot_calibConvergence(list(a = 1))
    Condition
      Error:
      ! could not find `member$bestvalit`; pass a DEoptim result or the trace vector.

