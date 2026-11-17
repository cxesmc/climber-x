!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : o c n _ m a s k _ m o d
!
!  Purpose : coarse ocean mask
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
module ocn_mask_mod

  use precision, only : wp
  use geo_params, only : f_crit, l_ocn_topo
  use geo_params, only : l_close_panama, l_close_bering, l_close_iberia 
  use geo_grid, only : ni, nj, n_topo_sur, i_topo_sur, j_topo_sur, i0_topo, i1_topo, j0_topo, j1_topo
  use fill_ocean_mod, only : fill_ocean_lowres
  use ocn_topology_mod, only : ocn_topology

  implicit none

  private
  public :: ocn_mask, f_crit_cell

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  o c n _ m a s k
  ! Purpose    :  derive the coarse ocean mask (the only place where it is decided): threshold the ocean fraction,
  !               apply explicit strait closures, optionally repair the connectivity (l_ocn_topo) and remove ocean
  !               cells not connected to the ocean origin
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine ocn_mask(hires_mask, hires_z_bed, lon, lat, lon_ocn_origin, lat_ocn_origin, &
      mask_ocn)

  implicit none

  integer, intent(inout) :: hires_mask(:,:)
  real(wp), intent(in) :: hires_z_bed(:,:)   !! high resolution bedrock elevation [m]
  real(wp), intent(in) :: lat(:)
  real(wp), intent(in) :: lon(:)
  real(wp), intent(in) :: lat_ocn_origin  
  real(wp), intent(in) :: lon_ocn_origin

  integer, intent(out) :: mask_ocn(:,:)      !! coarse ocean mask (1 = ocean, 0 = land)

  integer :: i, j, nocn
  real(wp), dimension(ni,nj) :: f_ocn
  real(wp), dimension(ni,nj) :: f_ocn_max
  logical, dimension(ni,nj) :: mask_fix

  integer, dimension(:), allocatable :: mask_cell
  integer, dimension(:,:), allocatable :: mask_ocn_connect


  allocate( mask_cell(n_topo_sur) )

  !$omp parallel do private(i,j,nocn,mask_cell)
  do j=1,nj
    do i=1,ni

      mask_cell = pack(hires_mask(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j)),.true.)

      ! ocean fraction: count number of ocean points in grid cell
      nocn = count(mask_cell.eq.2)
      f_ocn(i,j) = real(nocn,wp)/real(n_topo_sur,wp)

      ! largest ocean fraction the cell can be widened to without removing ice sheet: ocean, ice-free land and lakes
      f_ocn_max(i,j) = real(count(mask_cell.eq.1 .or. mask_cell.eq.2 .or. mask_cell.eq.4),wp)/real(n_topo_sur,wp)

      mask_ocn(i,j) = 0
      if (f_ocn(i,j).gt.0._wp) mask_ocn(i,j) = 1

      ! apply critical threshold for ocean/land
      if (f_ocn(i,j).lt.f_crit_cell(j)) then
        mask_ocn(i,j) = 0
      endif

    enddo
  enddo
  !$omp end parallel do

  ! explicit strait closures, these cells remain land also in the topology repair
  mask_fix = .false.
  ! fix for Panama
  if (l_close_panama) then
    ! close Panama
    mask_fix(16,23) = .true.
    mask_fix(17:18,22) = .true.
    mask_fix(19:20,21) = .true.
    mask_fix(21,20) = .true.
  endif
  ! option to close Bering 
  if (l_close_bering) then
    ! close Bering
    mask_fix(1:4,32) = .true.
  endif
  ! option to close the spurious strait through Iberia: the Alboran cell (-2.5E,37.5N) and the Biscay cell
  ! (-2.5E,42.5N) are both fractional ocean cells, so their common v-face is open and ~4 Sv circulate
  ! Atlantic -> Gibraltar -> Alboran -> Biscay -> Atlantic; removing the Biscay cell closes both this face
  ! and the u-face to the Gulf of Lion cell (2.5E,42.5N)
  if (l_close_iberia) then
    mask_fix(36,27) = .true.
  endif
  where (mask_fix) mask_ocn = 0

  ! topology-preserving repair of the coarse mask; cells re-opened below the threshold are widened to f_crit in
  ! hires_to_lowres, so cells that cannot reach f_crit without removing ice sheet can not be re-opened
  if (l_ocn_topo) then
    do j=1,nj
      do i=1,ni
        if (f_ocn_max(i,j).lt.f_crit_cell(j) .and. mask_ocn(i,j).eq.0) mask_fix(i,j) = .true.
      enddo
    enddo
    call ocn_topology(hires_mask, hires_z_bed, lon, lat, lon_ocn_origin, lat_ocn_origin, f_ocn, mask_fix, &
      mask_ocn)
  endif

  ! make sure that all ocean cells are connected on the low resolution grid
  ! where not, set ocean fraction to zero to avoid having 'ocean lakes'
  allocate( mask_ocn_connect(ni,nj))
  call fill_ocean_lowres(real(mask_ocn,wp),lon,lat,lon_ocn_origin,lat_ocn_origin, &
    mask_ocn_connect)
  where (mask_ocn_connect.eq.0)
    mask_ocn = 0
  endwhere

  !$omp parallel do private(i,j,nocn,mask_cell)
  do j=1,nj
    do i=1,ni

      mask_cell = pack(hires_mask(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j)),.true.)

      ! count number of ocean points in grid cell
      nocn = count(mask_cell.eq.2)
      ! if number of ocean cells is larger than 0 but the ocean fraction on the
      ! low resolution grid is zero, then set high resolution ocean cells to land and floating ice cells to grounded ice
      if (nocn.gt.0 .and. mask_ocn_connect(i,j).eq.0) then
        where (hires_mask(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j))==2)  ! ocean
          hires_mask(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j)) = 1   ! land
        endwhere
        where (hires_mask(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j))==3)  ! floating ice
          hires_mask(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j)) = 0   ! grounded ice
        endwhere
      endif

    enddo
  enddo
  !$omp end parallel do

  deallocate( mask_ocn_connect)

  deallocate( mask_cell )


  end subroutine ocn_mask 


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  f _ c r i t _ c e l l
  ! Purpose  :  minimum ocean fraction of an ocean cell in row j
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(wp) function f_crit_cell(j)

  implicit none

  integer, intent(in) :: j

  if (j==1 .or. j.eq.nj) then
    f_crit_cell = f_crit*2._wp  ! fixme?!
  else
    f_crit_cell = f_crit
  endif

  end function f_crit_cell

end module ocn_mask_mod
