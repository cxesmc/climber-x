!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : n 2 o _ e m i s _ m o d
!
!  Purpose : N2O emissions
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
module n2o_emis_mod

  use precision, only : wp
  use constants, only : T0
  use timer, only : sec_day
  use lnd_params, only : n2o_par

  implicit none

  private
  public :: n2o_emission

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  n 2 o _ e m i s s i o n 
  !   Purpose    :  n2o emissions 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine n2o_emission(t_soil, theta_w, theta_sat, &
                      n2o_emis)

  implicit none

  real(wp), intent(in) :: t_soil        !! soil temperature of the top layer [K]
  real(wp), intent(in) :: theta_w       !! liquid soil water content of the top layer [m3/m3]
  real(wp), intent(in) :: theta_sat     !! soil porosity [m3/m3]

  real(wp), intent(out) :: n2o_emis     !! n2o emissions [kg N2O-N/m2/s]

  real(wp) :: fac_t, fac_sm
  real(wp), parameter :: t_soil_ref = T0+10._wp  ! K


  ! temperature factor
  fac_t = n2o_par%q10_n2o**((t_soil-t_soil_ref)/10._wp) 

  ! soil moisture factor 
  fac_sm = (theta_w/theta_sat)**n2o_par%theta_exp

  ! N2O emissions
  n2o_emis = n2o_par%k_n2o/sec_day * fac_t * fac_sm     ! kg N2O-N/m2/s

  return

  end subroutine n2o_emission

end module n2o_emis_mod
