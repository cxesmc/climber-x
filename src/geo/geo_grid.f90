!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : g e o _ g r i d
!
!  Purpose : geography model grid
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
module geo_grid

  use precision, only : wp, dp
  use ncio
  use constants, only: pi, R_earth 
  use control, only : out_dir
  use climber_grid, only : ni, nj, lon, lat, dlon, dlat
  use coords, only : grid_class, grid_init
  use geo_params, only : geo_ref_file

  implicit none

  integer :: ni_topo, nj_topo
  integer :: n_topo_sur
  integer :: i_topo_sur
  integer :: j_topo_sur
  real(wp), dimension(:,:), allocatable :: area
  real(dp), dimension(:,:), allocatable :: area_dp
  real(wp), dimension(:), allocatable :: lon_topo, lat_topo
  integer, dimension(:), allocatable :: i0_topo, i1_topo
  integer, dimension(:), allocatable :: j0_topo, j1_topo
  integer, dimension(:), allocatable :: i_lowres, j_lowres

  !private
  !public :: geo_grid_init

contains

    subroutine geo_grid_init(grid, grid_hires)

    implicit none

    type(grid_class), intent(out) :: grid
    type(grid_class), intent(out) :: grid_hires

    integer :: i, j, ii, jj
    integer :: ppos
    real(wp) :: lon1, lon2, lat1, lat2
    real(wp) :: dlon_topo, dlat_topo


    ! generate grid object
    call grid_init(grid,name="GEO-5x5",mtype="latlon",units="degrees", x=real(lon,dp),y=real(lat,dp))

    ! read topography/bathymetry on high resolution
    ni_topo = nc_size(trim(geo_ref_file),"lon")
    nj_topo = nc_size(trim(geo_ref_file),"lat")
    allocate( lon_topo(ni_topo) )
    allocate( lat_topo(nj_topo) )
    allocate( area(ni_topo,nj_topo))
    allocate( area_dp(ni_topo,nj_topo))
    allocate( i_lowres(ni_topo))
    allocate( j_lowres(nj_topo))

    call nc_read(trim(geo_ref_file),"lon",lon_topo)
    call nc_read(trim(geo_ref_file),"lat",lat_topo)
    dlon_topo = lon_topo(2)-lon_topo(1)
    dlat_topo = lat_topo(2)-lat_topo(1)

    ! generate high resolution grid object
    ppos = scan(trim(geo_ref_file),".", BACK= .true.)-1
    call grid_init(grid_hires,name=trim(geo_ref_file(7:ppos)),mtype="latlon",units="degrees", &
      x0=real(lon_topo(1),dp),dx=real(dlon_topo,dp),nx=ni_topo,y0=real(lat_topo(1),dp),dy=real(dlat_topo,dp),ny=nj_topo)

    ! number of topo points in one model surface grid cell
    n_topo_sur = (ni_topo*nj_topo)/(ni*nj)
    i_topo_sur = ni_topo/ni
    j_topo_sur = nj_topo/nj

    allocate( i0_topo(ni) )
    allocate( i1_topo(ni) )
    allocate( j0_topo(nj) )
    allocate( j1_topo(nj) )

    ! derive indexes of high resolution bounding box of coarse grid
    i0_topo = 0
    i1_topo = 0
    do i=1,ni
      lon1 = lon(i)-0.5_wp*dlon
      lon2 = lon(i)+0.5_wp*dlon
      do ii=1,ni_topo-1
        if (lon_topo(ii+1).gt.lon1 .and. lon_topo(ii).le.lon1) then
          i0_topo(i) = ii+1
        endif
        if (lon_topo(ii+1).gt.lon2 .and. lon_topo(ii).le.lon2) then
          i1_topo(i) = ii
        endif
      enddo
      if (i0_topo(i).eq.0) i0_topo(i) = 1
      if (i1_topo(i).eq.0) i1_topo(i) = ni_topo
    enddo

    j0_topo = 0
    j1_topo = 0
    do j=1,nj
      lat1 = lat(j)-0.5_wp*dlat
      lat2 = lat(j)+0.5_wp*dlat
      do jj=1,nj_topo-1
        if (lat_topo(jj+1).gt.lat1 .and. lat_topo(jj).le.lat1) then
          j0_topo(j) = jj+1
        endif
        if (lat_topo(jj+1).gt.lat2 .and. lat_topo(jj).le.lat2) then
          j1_topo(j) = jj
        endif
      enddo
      if (j0_topo(j).eq.0) j0_topo(j) = 1
      if (j1_topo(j).eq.0) j1_topo(j) = nj_topo
    enddo

    ! derive correspondence between indexes on high and low resolution grids
    do i=1,ni
      lon1 = lon(i)-0.5_wp*dlon
      if (i.eq.1) lon1 = lon1-1.e-3_wp
      lon2 = lon(i)+0.5_wp*dlon
      do ii=1,ni_topo
        if (lon_topo(ii).gt.lon1 .and. lon_topo(ii).le.lon2) then
          i_lowres(ii) = i
        endif
      enddo
    enddo
    do j=1,nj
      lat1 = lat(j)-0.5_wp*dlat
      if (j.eq.1) lat1 = lat1-1.e-3_wp
      lat2 = lat(j)+0.5_wp*dlat
      do jj=1,nj_topo
        if (lat_topo(jj).gt.lat1 .and. lat_topo(jj).le.lat2) then
          j_lowres(jj) = j
        endif
      enddo
    enddo

    ! grid cell area 
    do i=1,ni_topo
      do j=1,nj_topo
        area(i,j) = 2._wp*pi*cos(pi*lat_topo(j)/180._wp)*R_earth*dlon_topo/360._wp*R_earth*pi*dlat_topo/180._wp
        area_dp(i,j) = 2._wp*pi*cos(pi*lat_topo(j)/180._wp)*R_earth*dlon_topo/360._wp*R_earth*pi*dlat_topo/180._wp
      enddo
    enddo


    return

    end subroutine geo_grid_init

end module geo_grid
