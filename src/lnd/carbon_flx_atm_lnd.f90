!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c a r b o n _ f l u x _ a t m _ l n d _ m o d
!
!  Purpose : compute net atmosphere-land carbon flux
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit and Andrey Ganopolski
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
module carbon_flux_atm_lnd_mod

  use precision, only : wp
  use timer, only : sec_day, day_mon
  use lnd_grid, only : ic_min, ic_peat, ic_shelf, ic_ice, ic_lake

  implicit none

  private
  public :: carbon_flux_atm_lnd

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  c a r b o n _ f l u x _ a t m _ l n d
  !   Purpose    :  compute net atmosphere-land carbon flux
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine carbon_flux_atm_lnd(f_veg, f_peat, f_ice_grd, f_shelf, f_lake,area, &
      npp_real, npp13_real, npp14_real, soil_resp, soil_resp13, soil_resp14, &
      fire_c_flux,fire_c13_flux,fire_c14_flux, &      
      Cflx_atm_lnd_2d, C13flx_atm_lnd_2d, C14flx_atm_lnd_2d, &
      Cflx_atm_lnd, C13flx_atm_lnd, C14flx_atm_lnd)

    implicit none

    real(wp), intent(in) :: f_veg, f_peat, f_ice_grd, f_shelf, f_lake, area
    real(wp), intent(in) :: npp_real, npp13_real, npp14_real
    real(wp), dimension(:), intent(in) :: soil_resp, soil_resp13, soil_resp14
    real(wp), intent(in) :: fire_c_flux, fire_c13_flux, fire_c14_flux    
    real(wp), intent(inout) :: Cflx_atm_lnd_2d, C13flx_atm_lnd_2d, C14flx_atm_lnd_2d
    real(wp), intent(inout) :: Cflx_atm_lnd, C13flx_atm_lnd, C14flx_atm_lnd

    real(wp) :: npp_ij, npp13_ij, npp14_ij, sresp_ij, sresp13_ij, sresp14_ij
    real(wp) :: fire_ij, fire13_ij, fire14_ij


    ! NPP
    npp_ij = npp_real * f_veg * sec_day*day_mon * area ! kgC/mon
    npp13_ij = npp13_real * f_veg * sec_day*day_mon * area 
    npp14_ij = npp14_real * f_veg * sec_day*day_mon * area 
    ! soil respiration
    sresp_ij = ( soil_resp(ic_min)*(f_veg-f_peat) &
      + soil_resp(ic_peat)*f_peat &
      + soil_resp(ic_shelf)*f_shelf &
      + soil_resp(ic_lake)*f_lake &
      + soil_resp(ic_ice)*f_ice_grd) &
      * sec_day*day_mon * area ! kgC/mon
    sresp13_ij = ( soil_resp13(ic_min)*(f_veg-f_peat) &
      + soil_resp13(ic_peat)*f_peat &
      + soil_resp13(ic_shelf)*f_shelf &
      + soil_resp13(ic_lake)*f_lake &
      + soil_resp13(ic_ice)*f_ice_grd) &
      * sec_day*day_mon * area ! kgC/mon
    sresp14_ij = ( soil_resp14(ic_min)*(f_veg-f_peat) &
      + soil_resp14(ic_peat)*f_peat &
      + soil_resp14(ic_shelf)*f_shelf &
      + soil_resp14(ic_lake)*f_lake &
      + soil_resp14(ic_ice)*f_ice_grd) &
      * sec_day*day_mon * area ! kgC/mon

    ! fire emissions
    fire_ij = fire_c_flux * f_veg * sec_day*day_mon * area ! kgC/mon
    fire13_ij = fire_c13_flux * f_veg * sec_day*day_mon * area
    fire14_ij = fire_c14_flux * f_veg * sec_day*day_mon * area

    Cflx_atm_lnd_2d = (npp_ij-sresp_ij-fire_ij) / (sec_day*day_mon)     ! kgC/s
    C13flx_atm_lnd_2d = (npp13_ij-sresp13_ij-fire13_ij) / (sec_day*day_mon)
    C14flx_atm_lnd_2d = (npp14_ij-sresp14_ij-fire14_ij) / (sec_day*day_mon)

    Cflx_atm_lnd = Cflx_atm_lnd + (npp_ij-sresp_ij-fire_ij)
    C13flx_atm_lnd = C13flx_atm_lnd + (npp13_ij-sresp13_ij-fire13_ij)
    C14flx_atm_lnd = C14flx_atm_lnd + (npp14_ij-sresp14_ij-fire14_ij)

    return

  end subroutine carbon_flux_atm_lnd

end module carbon_flux_atm_lnd_mod
