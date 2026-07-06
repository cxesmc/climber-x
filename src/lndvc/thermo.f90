module thermo_m
    ! Thermodynamic helper functions for the lndvc framework.
    !
    ! Brought in from src/main/constants.f90 (bodies unchanged) so the framework
    ! is self-contained: saturation specific humidity over ice and its
    ! temperature derivative, and air density.

    use precision, only : sp, dp, wp
    use const_m,   only : T0, Rd, Rv

    implicit none

    private
    public :: q_sat_i, dqsat_dT_i, rho_a, fqsat

    interface q_sat_i
        module procedure q_sat_i_sp
        module procedure q_sat_i_dp
    end interface q_sat_i

    interface fqsat
        module procedure fqsat_sp
        module procedure fqsat_dp
    end interface fqsat

contains

    pure function e_sat_i(temp) ! temp in K
        implicit none
        real(wp), intent(in) :: temp
        real(wp) :: e_sat_i
        e_sat_i = 6.1121d2 * exp( 22.587_wp * (temp-T0) / (273.86_wp + (temp-T0)))   ! Pa, ice
    end function e_sat_i

    pure function q_sat_i_sp(temp,p) ! temp in K, p in Pa
        implicit none
        real(sp), intent(in) :: temp, p
        real(sp) :: q_sat_i_sp
        q_sat_i_sp = 380.1726_wp * exp( 22.587_wp * (temp-T0) / (temp+0.71_wp)) / p  ! Pa, ice
    end function q_sat_i_sp

    pure function q_sat_i_dp(temp,p) ! temp in K, p in Pa
        implicit none
        real(dp), intent(in) :: temp, p
        real(dp) :: q_sat_i_dp
        q_sat_i_dp = 380.1726_wp * exp( 22.587_wp * (temp-T0) / (temp+0.71_wp)) / p  ! Pa, ice
    end function q_sat_i_dp

    pure function dqsat_dT_i(temp,p)
        implicit none
        real(wp), intent(in) :: temp, p
        real(wp) :: dqsat_dT_i
        real(wp) :: t, desat_dT_i
        real(wp), parameter :: Lv = 2834.d3 ! J/kg, ice
        t = temp - T0  ! °C
        desat_dT_i = Lv * e_sat_i(temp) / (Rv * temp**2)   ! Clausius - Clapeyron
        dqsat_dT_i = 0.622_wp * desat_dT_i / p  ! approximation
    end function dqsat_dT_i

    pure function rho_a(temp,p)
        implicit none
        real(wp), intent(in) :: temp, p
        real(wp) :: rho_a
        rho_a = p / (Rd*temp)
    end function rho_a

    pure function fqsat_sp(T,p)
        implicit none
        real(sp), intent(in) :: T, p
        real(sp) :: fqsat_sp
        real(sp), parameter :: Ti=248.
        real(sp) :: r_w, qsatw, qsati
        if (T.ge.T0)  then
          qsatw = 380.0047_wp * exp( 17.625_wp * (T-T0) / (T-30.11_wp)) / p  ! Pa, water
          fqsat_sp=qsatw
        elseif((T.gt.Ti).and.(T.lt.T0)) then
          r_w = 1.-((T0-T)/(T0-Ti))
          qsatw = 380.0047_wp * exp( 17.625_wp * (T-T0) / (T-30.11_wp)) / p  ! Pa, water
          qsati = 380.1726_wp * exp( 22.587_wp * (T-T0) / (T+0.71_wp)) / p  ! Pa, ice
          fqsat_sp=r_w*qsatw+(1.-r_w)*qsati
        else
          qsati = 380.1726_wp * exp( 22.587_wp * (T-T0) / (T+0.71_wp)) / p  ! Pa, ice
          fqsat_sp=qsati
        endif
    end function fqsat_sp

    pure function fqsat_dp(T,p)
        implicit none
        real(dp), intent(in) :: T, p
        real(dp) :: fqsat_dp
        real(dp), parameter :: Ti=248.
        real(dp) :: r_w, qsatw, qsati
        if (T.ge.T0)  then
          qsatw = 380.0047_wp * exp( 17.625_wp * (T-T0) / (T-30.11_wp)) / p  ! Pa, water
          fqsat_dp=qsatw
        elseif((T.gt.Ti).and.(T.lt.T0)) then
          r_w = 1.-((T0-T)/(T0-Ti))
          qsatw = 380.0047_wp * exp( 17.625_wp * (T-T0) / (T-30.11_wp)) / p  ! Pa, water
          qsati = 380.1726_wp * exp( 22.587_wp * (T-T0) / (T+0.71_wp)) / p  ! Pa, ice
          fqsat_dp=r_w*qsatw+(1.-r_w)*qsati
        else
          qsati = 380.1726_wp * exp( 22.587_wp * (T-T0) / (T+0.71_wp)) / p  ! Pa, ice
          fqsat_dp=qsati
        endif
    end function fqsat_dp

end module thermo_m
