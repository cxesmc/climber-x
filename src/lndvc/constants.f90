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
    real(wp), parameter :: lambda_i = 2.2_wp             !! W/m/K, thermal conductivity of ice
    real(wp), parameter :: Lf       = 334.e3_wp          !! J/kg, latent heat of fusion

end module const_m
