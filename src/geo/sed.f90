!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s e d _ m o d 
!
!  Purpose : read sediment thickness
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
module sed_mod

  use precision, only : wp, dp
  use geo_params, only : sed_file

  use ncio
  use coord, only : grid_class, grid_init
  use coord, only : map_scrip_class, map_scrip_init, map_scrip_field

  implicit none

  private
  public :: sed

contains


  subroutine sed(geo_grid, h_sed)

    implicit none

    type(grid_class), intent(in) :: geo_grid
    real(wp), intent(out) :: h_sed(:,:)

    integer :: ppos, spos
    integer :: ni_sed, nj_sed
    real(wp), dimension(:), allocatable :: lon_sed, lat_sed
    real(wp), dimension(:,:), allocatable :: h_sed_in
    type(grid_class) :: sed_grid
    type(map_scrip_class) :: maps_sed_to_geo


    ! read from file
    ni_sed = nc_size(trim(sed_file),"lon")
    nj_sed = nc_size(trim(sed_file),"lat")
    allocate( h_sed_in(ni_sed,nj_sed) )
    allocate( lon_sed(ni_sed) )
    allocate( lat_sed(nj_sed) )
    call nc_read(trim(sed_file),"lat",lat_sed)
    call nc_read(trim(sed_file),"lon",lon_sed)
    call nc_read(trim(sed_file),"h_sed",h_sed_in)

    ! grid definition
    spos = scan(trim(sed_file),"/", BACK= .true.)+1
    ppos = scan(trim(sed_file),".", BACK= .true.)-1
    call grid_init(sed_grid,name=trim(sed_file(spos:ppos)),mtype="latlon",units="degrees",x=real(lon_sed,dp),y=real(lat_sed,dp))
    ! map to geo grid
    call map_scrip_init(maps_sed_to_geo,sed_grid,geo_grid,method="bil",fldr="maps",load=.TRUE.,clean=.FALSE.)
    call map_scrip_field(maps_sed_to_geo,"h_sed",h_sed_in,h_sed,method="mean",missing_value=-9999._dp)
      !filt_method="gaussian",filt_par=[1._wp,geo_grid%G%dx])

    deallocate(h_sed_in, lon_sed, lat_sed)


  end subroutine sed

end module sed_mod
