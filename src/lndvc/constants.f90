module const_m
    ! Physical constants for the lndvc framework.
    !
    ! Brought in from src/main/constants.f90 so the framework is self-contained.
    ! Values must match the reference model. Expanded as more physics is ported.

    use precision, only : wp

    implicit none

    real(wp), parameter :: pi      = 3.141592653589793_wp
    real(wp), parameter :: T0      = 273.15_wp          !! freezing point [K]
    real(wp), parameter :: frac_vu = 0.45_wp            !! fraction of solar spectrum in visible+UV

end module const_m
