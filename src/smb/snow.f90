!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s n o w _ m o d
!
!  Purpose : snow layer update for SEMIX model
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2024 Potsdam Institute for Climate Impact Research,
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
module snow_mod

  use precision, only : wp
  use constants, only : Lf, T0, g
  use timer, only : sec_day
  use control, only : check_water
  use smb_params, only : dt, rdt
  use smb_params, only : snow_par

  implicit none

  private
  public :: snow_update

contains


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s n o w _ u p d a t e 
  !   Purpose    :  update snow layer
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine snow_update(mask_snow,evp, &
                              w_snow,w_snow_old,w_snow_max, &
                              h_snow)


    implicit none

    integer, intent(inout) :: mask_snow
    real(wp), intent(in) :: evp
    real(wp), intent(inout) :: w_snow, w_snow_old, w_snow_max
    real(wp), intent(out) :: h_snow


    ! snow water equivalent evolution
    ! remove sublimation, snowfall has been added already and snowmelt already removed and refreezing added during soil temperature update
    if (mask_snow.eq.1) then
      w_snow = w_snow - evp*dt  ! kg/m2
    endif

    if (check_water .and. w_snow.lt.0._wp) print *,'WARNING w_snow < 0',w_snow
    w_snow = max(0._wp,w_snow)  ! not strictly conserving water here!

    ! limit w_snow 
    if (w_snow.gt.snow_par%w_snow_max) w_snow = snow_par%w_snow_max

    ! update snow thickness
    h_snow = w_snow / snow_par%rho  ! m

    ! update snow mask
    if( w_snow .gt. snow_par%w_snow_crit ) then
      mask_snow = 1
    else
      mask_snow = 0
    endif

    ! save seasonal maximum snow swe
    if (w_snow.gt.w_snow_old .or. w_snow.eq.snow_par%w_snow_max) w_snow_max = w_snow

    return

  end subroutine snow_update

end module snow_mod


