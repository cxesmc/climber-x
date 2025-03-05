!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : q _ g e o _ m o d
!
!  Purpose : geothermal heat flux
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
module q_geo_mod

  use precision, only : wp, dp
  use geo_params, only : i_q_geo, q_geo_const, q_geo_file

  use ncio
  use coord, only : grid_class, grid_init
  use coord, only : map_class, map_init, map_field
  use coord, only : map_scrip_class, map_scrip_init, map_scrip_field

  implicit none

  private
  public :: geo_heat

contains


  subroutine geo_heat(geo_grid, q_geo)

    implicit none

    type(grid_class), intent(in) :: geo_grid
    real(wp), intent(out) :: q_geo(:,:)

    integer :: ppos, spos
    integer :: ni_qgeo, nj_qgeo
    real(wp), dimension(:), allocatable :: lon_qgeo, lat_qgeo
    real(wp), dimension(:,:), allocatable :: q_geo_in
    type(grid_class) :: qgeo_grid
    type(map_class) :: map_qgeo_to_geo
    type(map_scrip_class) :: maps_qgeo_to_geo


    if (i_q_geo.eq.1) then

      ! uniform value
      q_geo = q_geo_const 

    else if (i_q_geo.eq.2) then

      ! read from file
      ni_qgeo = nc_size(trim(q_geo_file),"lon")
      nj_qgeo = nc_size(trim(q_geo_file),"lat")
      allocate( q_geo_in(ni_qgeo,nj_qgeo) )
      allocate( lon_qgeo(ni_qgeo) )
      allocate( lat_qgeo(nj_qgeo) )
      call nc_read(trim(q_geo_file),"lat",lat_qgeo)
      call nc_read(trim(q_geo_file),"lon",lon_qgeo)
      call nc_read(trim(q_geo_file),"q_geo",q_geo_in)

      ! grid definition
      spos = scan(trim(q_geo_file),"/", BACK= .true.)+1
      ppos = scan(trim(q_geo_file),".", BACK= .true.)-1
      call grid_init(qgeo_grid,name=trim(q_geo_file(spos:ppos)),mtype="latlon",units="degrees",x=real(lon_qgeo,dp),y=real(lat_qgeo,dp))
      ! map to geo grid
      call map_scrip_init(maps_qgeo_to_geo,qgeo_grid,geo_grid,method="bil",fldr="maps",load=.TRUE.,clean=.FALSE.)
      call map_scrip_field(maps_qgeo_to_geo,"q_geo",q_geo_in,q_geo,method="mean",missing_value=-9999._dp)
        !filt_method="gaussian",filt_par=[1._wp,geo_grid%G%dx])

      deallocate(q_geo_in, lon_qgeo, lat_qgeo)

    endif

  end subroutine geo_heat

end module q_geo_mod
