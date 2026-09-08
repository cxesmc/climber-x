!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c o n s t a n t s
!
!  Purpose : definition of constants
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit
!
! This file is part of CLIMBER-X.
!
! CLIMBER-X is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! CLIMBER-X is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with CLIMBER-X.  If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module constants

  use precision, only : wp, sp, dp
  implicit none

  character(len=*), parameter :: map_gen = "coords"   !! mapping weight generator: "coords" (in-package) or "cdo" (external)

  real(wp), parameter :: pi = 3.141592653589793_wp

  real(wp), parameter :: rho_i   = 910._wp      !! kg/m3, density of ice
  real(wp), parameter :: rho_w   = 1000._wp     !! kg/m3, density of pure water
  real(wp), parameter :: rho_sw  = 1028._wp     !! kg/m3, density of seawater
  real(wp), parameter :: rho_as  = 3.3e3_wp     !! kg/m3, density of the asthenosphere

  real(wp), parameter :: cap_a   = 1000._wp     !! J/kg/K, specific heat capacity of air
  real(wp), parameter :: cap_i   = 2110._wp     !! J/kg/K, specific heat capacity of ice
  real(wp), parameter :: cap_w   = 4187._wp     !! J/kg/K, specific heat capacity of water

  real(wp), parameter :: lambda_a = 0.023_wp    !! W/m/K, thermal conductivity of air
  real(wp), parameter :: lambda_i = 2.2_wp      !! W/m/K, thermal conductivity of ice
  real(wp), parameter :: lambda_w = 0.6_wp      !! W/m/K, thermal conductivity of water

  real(wp), parameter :: emis_snow = 0.99_wp    !! longwave emissivity of snow
  real(wp), parameter :: emis_w = 0.98_wp       !! longwave emissivity of water

  real(wp), parameter :: Le      = 2501.e3_wp   !! J/kg, latent heat of evaporation
  real(wp), parameter :: Lf      = 334.e3_wp    !! J/kg, latent heat of fusion
  real(wp), parameter :: Ls      = Le + Lf      !! J/kg, latent heat of sublimation

  real(wp), parameter :: Rd      = 287.058_wp   !! J/kg/K, specific gas constant of dry air 
  real(wp), parameter :: Rv      = 461.5_wp     !! J/kg/K, specific gas constant of water vapor

  real(wp), parameter :: R_earth = 6371000._wp  !! m, radius of the Earth
  real(wp), parameter :: a_earth = 6378137._wp, b_earth = 6356752.3142_wp
  real(wp), parameter :: omega = 7.2921e-5_wp  !! earth angular velocity [1/s]
  real(wp), parameter :: fcoriolis = 2._wp*omega  !! coriolis factor, 2*omega [1/s]

  real(wp), parameter :: sigma = 5.670373e-8_wp   !! W/m2/K4, Stefan-Boltzmann constant
  real(wp), parameter :: karman = 0.4_wp        !! von Karman constant

  real(wp), parameter :: T0 = 273.15_wp
  real(wp), parameter :: g = 9.81_wp

  real(wp), parameter :: ppm_to_PgC = 2.123_wp    !! PgC/ppm, conversion factor

  real(wp), parameter :: c13_c12_std = 0.0112372_wp !! standard VPDB 13C/12C ratio
  real(wp), parameter :: c14_c_std = 1.170e-12_wp !! background preidustrial atmospheric value, Orr 2017 (OMIP)
  real(wp), parameter :: c14_tdec=1._wp/8267._wp/31556926._wp !! radiocarbon decay rate (1/s)
 
  real(wp), parameter :: k_boltz = 8.62e-5_wp     !! Boltzmann constant for Arrhenius function (eV/K) 

  real(wp), parameter :: z_sfl = 10._wp        !! surface layer height (m) 

  real(wp), parameter :: frac_vu = 0.45_wp      !! fraction of solar spectrum in visible and ultraviolet

  interface fqsat
    module procedure fqsat_sp
    module procedure fqsat_dp
  end interface fqsat

  interface q_sat_i
    module procedure q_sat_i_sp
    module procedure q_sat_i_dp
  end interface q_sat_i

  interface q_sat_w
    module procedure q_sat_w_sp
    module procedure q_sat_w_dp
  end interface q_sat_w


contains


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  e _ s a t _ w
  !   Purpose    :  compute saturation vapor pressure over water
  !                 from Alduchov, 1996
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pure function e_sat_w(temp) ! temp in K

    implicit none

    real(wp), intent(in) :: temp
    real(wp) :: e_sat_w

    e_sat_w = 6.1094d2 * exp( 17.625_wp * (temp-T0) / (243.04_wp + (temp-T0)))   ! Pa, water

  end function e_sat_w


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  e _ s a t _ i
  !   Purpose    :  compute saturation vapor pressure over ice
  !                 from Alduchov, 1996
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pure function e_sat_i(temp) ! temp in K

    implicit none

    real(wp), intent(in) :: temp
    real(wp) :: e_sat_i

    e_sat_i = 6.1121d2 * exp( 22.587_wp * (temp-T0) / (273.86_wp + (temp-T0)))   ! Pa, ice

  end function e_sat_i


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  q _ s a t _ w _ s p
  !   Purpose    :  compute saturation mixing ratio over water
  !                 from Alduchov, 1996
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pure function q_sat_w_sp(temp,p) ! temp in K, p in Pa

    implicit none

    real(sp), intent(in) :: temp, p
    real(sp) :: q_sat_w_sp

!    e_s = 6.1094d2 * exp( 17.625_wp * (temp-T0) / (243.04_wp + (temp-T0)))   ! Pa, water
!    q_sat_w_sp = 0.622_wp * e_s / (p - 0.378_wp * e_s)   ! kg/kg
    ! approximate
    q_sat_w_sp = 380.0047_wp * exp( 17.625_wp * (temp-T0) / (temp-30.11_wp)) / p  ! Pa, water

  end function q_sat_w_sp


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  q _ s a t _ w _ d p
  !   Purpose    :  compute saturation mixing ratio over water
  !                 from Alduchov, 1996
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pure function q_sat_w_dp(temp,p) ! temp in K, p in Pa

    implicit none

    real(dp), intent(in) :: temp, p
    real(dp) :: q_sat_w_dp

!    e_s = 6.1094d2 * exp( 17.625_wp * (temp-T0) / (243.04_wp + (temp-T0)))   ! Pa, water
!    q_sat_w_dp = 0.622_wp * e_s / (p - 0.378_wp * e_s)   ! kg/kg
    ! approximate
    q_sat_w_dp = 380.0047_wp * exp( 17.625_wp * (temp-T0) / (temp-30.11_wp)) / p  ! Pa, water

  end function q_sat_w_dp


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  q _ s a t _ i _ s p
  !   Purpose    :  compute saturation mixing ratio over ice
  !                 from Alduchov, 1996
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pure function q_sat_i_sp(temp,p) ! temp in K, p in Pa

    implicit none

    real(sp), intent(in) :: temp, p
    real(sp) :: q_sat_i_sp

!    e_s = 6.1121d2 * exp( 22.587_wp * (temp-T0) / (273.86_wp + (temp-T0)))   ! Pa, ice
!    q_sat_i_sp = 0.622_wp * e_s / (p - 0.378_wp * e_s)   ! kg/kg
    ! approximate
    q_sat_i_sp = 380.1726_wp * exp( 22.587_wp * (temp-T0) / (temp+0.71_wp)) / p  ! Pa, ice

  end function q_sat_i_sp


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  q _ s a t _ i _ d p
  !   Purpose    :  compute saturation mixing ratio over ice
  !                 from Alduchov, 1996
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pure function q_sat_i_dp(temp,p) ! temp in K, p in Pa

    implicit none

    real(dp), intent(in) :: temp, p
    real(dp) :: q_sat_i_dp

!    e_s = 6.1121d2 * exp( 22.587_wp * (temp-T0) / (273.86_wp + (temp-T0)))   ! Pa, ice
!    q_sat_i_dp = 0.622_wp * e_s / (p - 0.378_wp * e_s)   ! kg/kg
    ! approximate
    q_sat_i_dp = 380.1726_wp * exp( 22.587_wp * (temp-T0) / (temp+0.71_wp)) / p  ! Pa, ice

  end function q_sat_i_dp


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  q _ t o _ e
  !   Purpose    :  convert air specific humidity (kg/kg) to vapor pressure (Pa)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pure function q_to_e(q,p)

    implicit none

    real(wp), intent(in) :: q, p
    real(wp) :: q_to_e

     q_to_e = q * p / (0.622_wp + 0.378_wp *q) ! Pa

  end function q_to_e


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  d q s a t _ d T _ w
  !   Purpose    :  derivative of qsat with temperature
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pure function dqsat_dT_w(temp,p)

    implicit none

    real(wp), intent(in) :: temp, p
    real(wp) :: dqsat_dT_w

    real(wp) :: t, desat_dT_w, Lv

     t = temp - T0  ! °C

!     Lv = ( 2500.8_wp - 2.36_wp * t + 0.0016_wp * t**2 - 0.00006_wp * t**3 ) * 1000._wp ! J/kg, water
     Lv = ( 2500.8_wp - 2.36_wp * t ) * 1000._wp ! J/kg, water

     desat_dT_w = Lv * e_sat_w(temp) / (Rv * temp**2)   ! Clausius - Clapeyron
!     dqsat_dT_w = p * 0.622_wp *  desat_dT_w / (p - 0.387_wp * e_sat_w(temp))**2
     dqsat_dT_w = 0.622_wp * desat_dT_w / p  ! approximation


  end function dqsat_dT_w


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  d q s a t _ d T _ i
  !   Purpose    :  derivative of qsat with temperature
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pure function dqsat_dT_i(temp,p)

    implicit none

    real(wp), intent(in) :: temp, p
    real(wp) :: dqsat_dT_i

    real(wp) :: t, desat_dT_i

    real(wp), parameter :: Lv = 2834.d3 ! J/kg, ice

     t = temp - T0  ! °C

     desat_dT_i = Lv * e_sat_i(temp) / (Rv * temp**2)   ! Clausius - Clapeyron
!     dqsat_dT_i = p * 0.622_wp *  desat_dT_i / (p - 0.387_wp * e_sat_i(temp))**2
     dqsat_dT_i = 0.622_wp * desat_dT_i / p  ! approximation


  end function dqsat_dT_i

 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  r h o _ a
  !   Purpose    :  compute air density (assuming dry air)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pure function rho_a(temp,p)

    implicit none

    real(wp), intent(in) :: temp, p
    real(wp) :: rho_a

    rho_a = p / (Rd*temp)

  end function rho_a


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  f q s a t _ s p
  !   Purpose    :  Saturated specific humidity for water, ice and intermediate case
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pure function FQSAT_sp(T,p)

    implicit none
      
    real(sp), intent(in) :: T, p
    real(sp) :: fqsat_sp

    real(sp), parameter :: Ti=248.
    real(sp) :: r_w, qsatw, qsati

    if (T.ge.T0)  then 
      r_w = 1.
    elseif((T.gt.Ti).and.(T.lt.T0)) then
      r_w = 1.-((T0-T)/(T0-Ti)) 
    else
      r_w= 0.
    endif            

    qsatw = 380.0047_wp * exp( 17.625_wp * (T-T0) / (T-30.11_wp)) / p  ! Pa, water
    qsati = 380.1726_wp * exp( 22.587_wp * (T-T0) / (T+0.71_wp)) / p  ! Pa, ice
    FQSAT_sp=r_w*QSATW+(1.-r_w)*QSATI

  end function fqsat_sp

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  f q s a t _ d p
  !   Purpose    :  Saturated specific humidity for water, ice and intermediate case
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pure function fqsat_dp(T,p)

    implicit none
      
    real(dp), intent(in) :: T, p
    real(dp) :: fqsat_dp

    real(dp), parameter :: Ti=248.
    real(dp) :: r_w, qsatw, qsati

    if (T.ge.T0)  then 
      qsatw = 380.0047_wp * exp( 17.625_wp * (T-T0) / (T-30.11_wp)) / p  ! Pa, water
      FQSAT_dp=qsatw
    elseif((T.gt.Ti).and.(T.lt.T0)) then
      r_w = 1.-((T0-T)/(T0-Ti)) 
      qsatw = 380.0047_wp * exp( 17.625_wp * (T-T0) / (T-30.11_wp)) / p  ! Pa, water
      qsati = 380.1726_wp * exp( 22.587_wp * (T-T0) / (T+0.71_wp)) / p  ! Pa, ice
      FQSAT_dp=r_w*qsatw+(1.-r_w)*qsati
    else
      qsati = 380.1726_wp * exp( 22.587_wp * (T-T0) / (T+0.71_wp)) / p  ! Pa, ice
      FQSAT_dp=qsati
    endif            

  end function fqsat_dp

end module constants
