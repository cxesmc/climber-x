!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : t o p o _ f i l l _ m o d
!
!  Purpose : fill topography
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
module topo_fill_mod

  use precision, only : wp 
  !$ use omp_lib

  implicit none

  private
  public :: topo_fill  

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  t o p o _ f i l l 
  !   Purpose    :  fill topography
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine topo_fill(z_topo, mask, &
      z_topo_fill)

    implicit none

    real(wp), intent(in) :: z_topo(:,:)
    integer, intent(in) :: mask(:,:)
    real(wp), intent(out) :: z_topo_fill(:,:)

    integer :: i, j, ii, jj, ni, nj
    logical :: something_done
    real(wp) :: eps
    real(wp), dimension(:,:), allocatable :: z_topo_ext
    real(wp), dimension(:,:), allocatable :: z_topo_fill_ext


    ! get grid size
    ni = size(z_topo,1)
    nj = size(z_topo,2)

    ! allocate
    allocate(z_topo_ext(0:ni+1,0:nj+1))
    allocate(z_topo_fill_ext(0:ni+1,0:nj+1))

    ! extend z_topography
    z_topo_ext(1:ni,1:nj) = z_topo
    z_topo_ext(:,0) = 9999._wp ! south pole
    z_topo_ext(:,nj+1) = 9999._wp ! north pole
    z_topo_ext(0,:) = z_topo_ext(ni,:)
    z_topo_ext(ni+1,:) = z_topo_ext(1,:)

    !--------------------------------------------------------------
    ! fill depressions using algorithm from Planchon & Darboux (2001)
    !--------------------------------------------------------------

    ! initalisation of the surface to large altitudes
    where (mask.ne.2 .and. mask.ne.3)
      z_topo_fill_ext(1:ni,1:nj) = 9999._wp
    elsewhere
      z_topo_fill_ext(1:ni,1:nj) = z_topo_ext(1:ni,1:nj)
    endwhere
    z_topo_fill_ext(:,0) = 9999._wp ! south pole
    z_topo_fill_ext(:,nj+1) = 9999._wp ! north pole
    z_topo_fill_ext(0,:) = z_topo_fill_ext(ni,:)
    z_topo_fill_ext(ni+1,:) = z_topo_fill_ext(1,:)

    eps = 1.e-4_wp   ! small slope to ensure all points are routed to the ocean
    ! fill algorithm, Table 5 in Planchon & Darboux (2001)
    something_done = .true.
    do while (something_done)
      something_done = .false.
      !$omp parallel do private(i,j,ii,jj)
      do j=1,nj
        do i=1,ni
          if (mask(i,j).ne.2 .and. mask(i,j).ne.3) then ! only internal points
            if (z_topo_fill_ext(i,j)>z_topo_ext(i,j)) then
              do jj=j-1,j+1
                do ii=i-1,i+1
                  if (ii.eq.i .and. jj.eq.j) then
                  else
                    if (z_topo_ext(i,j).ge.(z_topo_fill_ext(ii,jj)+eps)) then
                      z_topo_fill_ext(i,j) = z_topo_ext(i,j)
                      something_done = .true.
                    else if (z_topo_fill_ext(i,j).gt.(z_topo_fill_ext(ii,jj)+eps)) then
                      z_topo_fill_ext(i,j) = z_topo_fill_ext(ii,jj)+eps
                      something_done = .true.
                    endif
                  endif
                enddo
              enddo
            endif
          endif
        enddo
      enddo
      !$omp end parallel do
    enddo

    z_topo_fill = z_topo_fill_ext(1:ni,1:nj)    ! for output

    deallocate(z_topo_ext)
    deallocate(z_topo_fill_ext)

   return

  end subroutine topo_fill


end module topo_fill_mod


