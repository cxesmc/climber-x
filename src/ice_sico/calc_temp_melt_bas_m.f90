!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ t e m p _ m e l t _ b a s _ m
!
!> @file
!!
!! Computation of the melting and basal temperatures.
!!
!! @section Copyright
!!
!! Copyright 2009-2017 Ralf Greve
!!
!! @section License
!!
!! This file is part of SICOPOLIS.
!!
!! SICOPOLIS is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! SICOPOLIS is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with SICOPOLIS.  If not, see <http://www.gnu.org/licenses/>.
!<
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------------------------------------------------------------
!> Computation of the melting and basal temperatures.
!<------------------------------------------------------------------------------
module calc_temp_melt_bas_m

  use sico_types_m
  use sico_state
  use sico_grid_mod

  implicit none

  private
  public :: calc_temp_melt, calc_temp_bas

contains

!-------------------------------------------------------------------------------
!> Computation of the melting temperatures.
!<------------------------------------------------------------------------------
  subroutine calc_temp_melt(st,grd)

  implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class), intent(in) :: grd

  integer :: i, j, kc, kt

!-------- Compute the melting temperatures --------

  do i=0, grd%IMAX
  do j=0, grd%JMAX

     do kt=0, grd%KTMAX
        st%temp_t_m(kt,j,i) = -(BETA*st%H_c(j,i)+grd%atm2(kt)*st%H_t(j,i))
     end do

     do kc=0, grd%KCMAX
        st%temp_c_m(kc,j,i) = -grd%atm1(kc)*st%H_c(j,i)
     end do

  end do
  end do

  end subroutine calc_temp_melt

!-------------------------------------------------------------------------------
!> Computation of the basal temperatures.
!<------------------------------------------------------------------------------
  subroutine calc_temp_bas(st,grd)

  implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class), intent(in) :: grd

  integer :: i, j

!-------- Computation of the basal temperatures --------

  do i=0, grd%IMAX
  do j=0, grd%JMAX

     if ( (st%maske(j,i) == 0).or.(st%maske(j,i) == 3) ) then
                                   ! glaciated land or floating ice

        if (st%n_cts(j,i) == -1) then   ! cold ice base

           st%temp_b(j,i)  = st%temp_c(0,j,i)
           st%temph_b(j,i) = st%temp_c(0,j,i) - st%temp_c_m(0,j,i)
                          ! relative to the pressure melting point

        else   ! n_cts(j,i) == 0 or 1, temperate ice base

           st%temp_b(j,i)  = st%temp_t_m(0,j,i)
           st%temph_b(j,i) = 0.0_wp
                          ! relative to the pressure melting point

        end if

     else   ! st%maske(j,i) == 1 or 2, ice-free land or sea

        st%temp_b(j,i)  = st%temp_c(0,j,i)
        st%temph_b(j,i) = st%temp_c(0,j,i) - st%temp_c_m(0,j,i)
                       ! relative to the pressure melting point

     end if

  end do
  end do

  end subroutine calc_temp_bas

!-------------------------------------------------------------------------------

end module calc_temp_melt_bas_m
!
