!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : i c e _ i d _ m o d 
!
!  Purpose : read mask of ice ids
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
module ice_id_mod

  use precision, only : wp, dp

  use ncio
  use coord, only : grid_class, grid_init
  use coord, only : map_scrip_class, map_scrip_init, map_scrip_field

  implicit none

  character (len=256) :: ice_id_file = "input/mask_id_ice.nc"
  integer :: ice_id_eu = 1
  integer :: ice_id_grl = 2
  integer :: ice_id_es = 3
  integer :: ice_id_tib = 4
  integer :: ice_id_na = 5
  integer :: ice_id_rm = 6
  integer :: ice_id_al = 7

  private
  public :: set_ice_id
  public :: ice_id_eu, ice_id_grl, ice_id_es, ice_id_tib, ice_id_na, ice_id_rm, ice_id_al

contains


  subroutine set_ice_id(ice_grid, ice_id_mask)

    implicit none

    type(grid_class), intent(in) :: ice_grid
    integer, intent(out) :: ice_id_mask(:,:)

    integer :: ppos, spos
    integer :: ni_ice_id, nj_ice_id
    real(wp), dimension(:), allocatable :: lon_ice_id, lat_ice_id
    integer, dimension(:,:), allocatable :: ice_id_in
    type(grid_class) :: ice_id_grid
    type(map_scrip_class) :: maps_ice_id_to_ice


    ! read from file
    ni_ice_id = nc_size(trim(ice_id_file),"x")
    nj_ice_id = nc_size(trim(ice_id_file),"y")
    allocate( ice_id_in(ni_ice_id,nj_ice_id) )
    allocate( lon_ice_id(ni_ice_id) )
    allocate( lat_ice_id(nj_ice_id) )
    call nc_read(trim(ice_id_file),"y",lat_ice_id)
    call nc_read(trim(ice_id_file),"x",lon_ice_id)
    call nc_read(trim(ice_id_file),"ice_id",ice_id_in)

    ! grid definition
    spos = scan(trim(ice_id_file),"/", BACK= .true.)+1
    ppos = scan(trim(ice_id_file),".", BACK= .true.)-1
    call grid_init(ice_id_grid,name=trim(ice_id_file(spos:ppos)),mtype="latlon",units="degrees",x=real(lon_ice_id,dp),y=real(lat_ice_id,dp))
    ! map to ice grid
    call map_scrip_init(maps_ice_id_to_ice,ice_id_grid,ice_grid,method="nn",fldr="maps",load=.TRUE.,clean=.FALSE.)
    call map_scrip_field(maps_ice_id_to_ice,"ice_id",ice_id_in,ice_id_mask,method="mean",missing_value=-9999._dp)

    deallocate(ice_id_in, lon_ice_id, lat_ice_id)


  end subroutine set_ice_id

end module ice_id_mod
