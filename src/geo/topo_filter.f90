!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : t o p o _ f i l t e r _ mo d
!
!  Purpose : filter topography
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
module topo_filter_mod

  use precision, only : wp
  use constants, only : pi, R_earth
  use coords, only : grid_class

  implicit none
    
  real(wp), dimension(:,:,:), allocatable :: weight_save

  private
  public :: topo_filter

contains

  subroutine topo_filter(grid, z_topo, sigma_filter, &
      z_topo_fil)

    implicit none

    type(grid_class), intent(in) :: grid
    real(wp), intent(in)  :: z_topo(:,:)
    real(wp), intent(in)  :: sigma_filter
    real(wp), intent(out) :: z_topo_fil(:,:)

    integer :: i, j, ii, jj, i_f, j_f, n_filter
    real(wp) :: dx, dy, dist, weight, sum_weight, ztopofil

    logical, save :: firstcall = .true.


    n_filter     = ceiling(2.0_wp*sigma_filter/(grid%G%dx*pi*R_earth/180._wp*1e-3_wp))

    if (firstcall) then
      ! compute and save weights for filter
      allocate(weight_save(grid%G%ny,-n_filter:n_filter,-n_filter:n_filter))
      !$omp parallel do private(j,ii,jj,j_f,weight,sum_weight,dx,dy,dist)
      do j=1,grid%G%ny
        sum_weight = 0.0_wp
        weight_save(j,:,:) = 0._wp
        do ii=-n_filter, n_filter
          do jj=-n_filter, n_filter
            j_f = j+jj
            if (j_f.ge.1 .and. j_f.le.grid%G%ny) then
              dx = abs(ii)*grid%G%dx*pi*R_earth/180._wp*cos(grid%y(1,j_f)*pi/180._wp)*1.e-3_wp ! km 
              dy = abs(jj)*grid%G%dy*pi*R_earth/180._wp*1.e-3_wp ! km 
              dist       = sqrt(dx**2+dy**2)
              weight = exp(-(dist/sigma_filter)**2)
              sum_weight = sum_weight + weight
              weight_save(j,ii,jj) = weight
            endif
          end do
        end do
        weight_save(j,:,:) = weight_save(j,:,:)/sum_weight
      end do
      !$omp end parallel do
    endif

    !$omp parallel do collapse(2) private(i,j,ii,jj,i_f,j_f,ztopofil)
    do j=1,grid%G%ny
      do i=1,grid%G%nx
        ztopofil = 0._wp
        do jj=-n_filter, n_filter
          do ii=-n_filter, n_filter
            i_f = i+ii
            j_f = j+jj
            if (i_f <  1) i_f = grid%G%nx + i_f 
            if (i_f > grid%G%nx) i_f = i_f - grid%G%nx
            if (j_f.ge.1 .and. j_f.le.grid%G%ny) then
              ztopofil = ztopofil + weight_save(j,ii,jj)*z_topo(i_f,j_f)
            endif
          end do
          z_topo_fil(i,j) = ztopofil
        end do
      end do
    end do
    !$omp end parallel do

    if (firstcall) firstcall=.false.

    return

  end subroutine topo_filter

end module topo_filter_mod
