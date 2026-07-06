module const_m
    ! Physical constants for the lndvc framework.
    !
    ! Brought in from src/main/constants.f90 so the framework is self-contained.
    ! Values must match the reference model. Expanded as more physics is ported.

    use precision, only : wp

    implicit none

    real(wp), parameter :: pi       = 3.141592653589793_wp
    real(wp), parameter :: T0       = 273.15_wp          !! freezing point [K]
    real(wp), parameter :: frac_vu  = 0.45_wp            !! fraction of solar spectrum in visible+UV
    real(wp), parameter :: rho_i    = 910._wp            !! kg/m3, density of ice
    real(wp), parameter :: cap_i    = 2110._wp           !! J/kg/K, specific heat capacity of ice
    real(wp), parameter :: cap_a    = 1000._wp           !! J/kg/K, specific heat capacity of air
    real(wp), parameter :: lambda_i = 2.2_wp             !! W/m/K, thermal conductivity of ice
    real(wp), parameter :: Lf       = 334.e3_wp          !! J/kg, latent heat of fusion
    real(wp), parameter :: Le       = 2501.e3_wp         !! J/kg, latent heat of evaporation
    real(wp), parameter :: Ls       = Le + Lf            !! J/kg, latent heat of sublimation
    real(wp), parameter :: Rd       = 287.058_wp         !! J/kg/K, gas constant of dry air
    real(wp), parameter :: Rv       = 461.5_wp           !! J/kg/K, gas constant of water vapor
    real(wp), parameter :: sigma    = 5.670373e-8_wp     !! W/m2/K4, Stefan-Boltzmann constant
    real(wp), parameter :: karman   = 0.4_wp             !! von Karman constant
    real(wp), parameter :: g        = 9.81_wp            !! m/s2, gravitational acceleration
    real(wp), parameter :: z_sfl    = 100._wp            !! m, surface layer height

end module const_m
