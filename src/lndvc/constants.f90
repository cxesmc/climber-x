module const_m
    ! Physical constants for the lndvc framework.
    !
    ! These are re-exported from src/main/constants.f90 rather than declared
    ! again. The rule was always that the values must match the reference
    ! model, and a second table cannot enforce it: Lf stayed at 334.e3 while
    ! the reference moved to the standard 333.5e3, and z_sfl stayed at 100 m
    ! after the reference went to 10 m. Naming the framework's constants here
    ! keeps every `use const_m` working unchanged while there is one
    ! definition site -- the same move as delegating lndvc_grid to lnd_grid.
    !
    ! Add a name here when the framework needs it. If a quantity genuinely
    ! belongs to lndvc alone, declare it in the module that uses it.

    use constants, only : pi, T0, frac_vu, rho_i, cap_i, cap_a, lambda_i, &
                          Lf, Le, Ls, Rd, Rv, sigma, karman, g, z_sfl

    implicit none

end module const_m
