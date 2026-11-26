!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module : g r i d _ n c
!
!  Purpose : write CF grid-mapping information of a (projected) grid
!            to a NetCDF file
!
!  This file is part of CLIMBER-X.
!
!  CLIMBER-X is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  CLIMBER-X is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with CLIMBER-X.  If not, see <https://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module grid_nc

  use ncio
  use coords, only : grid_class

  implicit none

  ! name of the CF grid-mapping variable written by nc_write_map
  character(len=*), parameter :: grid_mapping_name = "crs"

  private
  public :: grid_mapping_name, grid_nc_write_map

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  g r i d _ n c _ w r i t e _ m a p
  ! Purpose    :  write the grid-mapping variable (crs) and the 2D
  !               lon/lat fields of grid to file fnm, as yelmo_grid_write
  !               does. The file must exist, be closed, and have the x/y
  !               dimensions dim_x/dim_y. x_e/y_n are written in the units
  !               of the x/y axes, given by xy_fac (axis units per metre,
  !               1e-3 for km axes).
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine grid_nc_write_map(fnm, grid, dim_x, dim_y, xy_fac)

    implicit none

    character(len=*), intent(in) :: fnm
    type(grid_class), intent(in) :: grid
    character(len=*), intent(in) :: dim_x, dim_y
    double precision, intent(in) :: xy_fac


    if (.not. grid%cs%is_projection) return

    call nc_write_map(fnm, grid%cs%mtype, grid%cs%proj%lambda, phi=grid%cs%proj%phi, &
      alpha=grid%cs%proj%alpha, x_e=grid%cs%proj%x_e*xy_fac, y_n=grid%cs%proj%y_n*xy_fac, &
      k0=grid%cs%proj%k0, &
      is_sphere=grid%cs%planet%is_sphere, semi_major_axis=grid%cs%planet%a, &
      inverse_flattening=grid%cs%planet%inverse_flattening)

    call nc_write_attr(fnm, dim_x, "standard_name", "projection_x_coordinate")
    call nc_write_attr(fnm, dim_y, "standard_name", "projection_y_coordinate")

    call nc_write(fnm, "lon2D", grid%lon, dim1=dim_x, dim2=dim_y, &
      long_name="longitude", grid_mapping=grid_mapping_name)
    call nc_write_attr(fnm, "lon2D", "units", "degrees_east")
    call nc_write_attr(fnm, "lon2D", "standard_name", "longitude")
    call nc_write(fnm, "lat2D", grid%lat, dim1=dim_x, dim2=dim_y, &
      long_name="latitude", grid_mapping=grid_mapping_name)
    call nc_write_attr(fnm, "lat2D", "units", "degrees_north")
    call nc_write_attr(fnm, "lat2D", "standard_name", "latitude")

    return

  end subroutine grid_nc_write_map

end module grid_nc
