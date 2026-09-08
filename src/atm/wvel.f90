!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : w v e l _ m o d
!
!  Purpose : vertical velocities at particular levels
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Andrey Ganopolski and Matteo Willeit
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
module wvel_mod

  use atm_params, only : wp
  use atm_params, only : c_woro, c_weff, nsmooth_weff
  use atm_grid, only : im, jm, kweff
  use smooth_atm_mod, only : smooth2
  !$ use omp_lib

  implicit none

  private
  public :: wvel

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  w v e l
  !   Purpose    :  vertical velocities for parameterisations
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine wvel(w3_nt, wsyn, wind, sigoro, &
      wcld, woro, weff)

    implicit none

    real(wp), intent(in ) :: w3_nt(:,:,:)
    real(wp), intent(in ) :: wsyn(:,:)
    real(wp), intent(in ) :: wind(:,:)
    real(wp), intent(in ) :: sigoro(:,:)

    real(wp), intent(out) :: wcld(:,:)
    real(wp), intent(out) :: woro(:,:)
    real(wp), intent(out) :: weff(:,:)

    integer :: i, j


    do j=1,jm
      do i=1,im

        ! large scale mean vertical velocity at cloud level, with or without contribution by topography 
        wcld(i,j) = w3_nt(i,j,kweff(i,j))

        ! vertical velocity due to subgrid scale orography
        woro(i,j) = c_woro*wind(i,j)*sigoro(i,j)

        ! effective vertical velocity at cloud level (mean + synoptic + orographic)
        weff(i,j) = wcld(i,j) + c_weff*(wsyn(i,j)+woro(i,j))   

      enddo 
    enddo

    ! smooth
    call smooth2(wcld,nsmooth_weff) 
    call smooth2(woro,nsmooth_weff) 
    call smooth2(weff,nsmooth_weff) 


    return 

  end subroutine wvel

end module wvel_mod
