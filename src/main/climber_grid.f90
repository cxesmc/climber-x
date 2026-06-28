!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c l i m b e r _ g r i d
!
!  Purpose : climber grid definition
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
module climber_grid

  use precision, only : wp, dp
  use constants, only : pi, r_earth
  use control, only : in_dir, out_dir
  use coords, only : grid_class, grid_init
  use ncio
  use nml

  implicit none

  integer, parameter :: ni = 72
  integer, parameter :: nj = 36
  real(wp), parameter :: lon0 = -180._wp
  real(wp), parameter :: lat0 = -90._wp 
  real(wp), parameter :: dlon = 360._wp/ni 
  real(wp), parameter :: dlat = 180._wp/nj
  real(wp) :: lon(ni)
  real(wp) :: lonu(ni+1)
  real(wp) :: lat(nj)
  real(wp) :: latv(nj+1)
  real(wp) :: area(ni,nj)
  real(wp) :: area_tot
  integer, dimension(ni,nj) :: basin_mask
  integer, dimension(ni,nj) :: basin_mask2
  integer, parameter :: i_atlantic = 1
  integer, parameter :: i_pacific = 2
  integer, parameter :: i_indian = 3
  integer, parameter :: i_southern = 4
  integer, parameter :: i_medi = 5

contains

  subroutine climber_grid_init

    implicit none

    integer :: i, j

    ! 
    do i=1,ni
     lon(i) = lon0+dlon/2._wp + (i-1)*dlon
    end do
    do i=1,ni+1
     lonu(i) = lon0 + (i-1)*dlon
    end do
    do j=1,nj
     lat(j) = lat0+dlat/2._wp + (j-1)*dlat
    end do
    do j=1,nj+1
     latv(j) = lat0 + (j-1)*dlat
    end do

    ! grid cell area
    do i=1,ni
     do j=1,nj
      area(i,j) = 2._wp*pi*cos(pi*lat(j)/180._wp)*dlon/360._wp*pi*dlat/180._wp*R_earth**2;
     enddo
    enddo

    ! total surface area of the Earth
    area_tot = sum(area)

    ! read ocean basin mask
    call nc_read(trim(in_dir)//"basin_mask_5x5.nc","basin_mask",basin_mask)
    call nc_read(trim(in_dir)//"basin_mask_5x5.nc","basin_mask2",basin_mask2)

  return

  end subroutine climber_grid_init


  subroutine ice_grid_init(domain, ice_grid)

    implicit none

    character(len=*), intent(in)  :: domain
    type(grid_class), intent(out) :: ice_grid

    integer             :: nx, ny
    character(len=256)  :: fnm
    character (len=256) :: type
    real(wp)            :: dx, dy
    real(wp)            :: x0, y0
    real(wp)            :: lambda, phi
    real(wp)            :: alpha

    ! define grid based on domain name
    fnm = trim(out_dir)//"/ice_grids.nml"
    call nml_read(fnm,domain,"grid_type"  ,type  )
    call nml_read(fnm,domain,"grid_dx"    ,dx    )
    call nml_read(fnm,domain,"grid_nx"    ,nx    )
    call nml_read(fnm,domain,"grid_ny"    ,ny    )
    call nml_read(fnm,domain,"grid_x0"    ,x0    )
    call nml_read(fnm,domain,"grid_y0"    ,y0    )
    call nml_read(fnm,domain,"grid_lambda",lambda)
    call nml_read(fnm,domain,"grid_phi"   ,phi   )
    if (trim(type)=="stereographic") then
      ! alpha needed for oblique stereo
      call nml_read(fnm,domain,"grid_alpha",alpha)
    else
      ! dummy value
      alpha = 1.e-3_wp
    endif
    dy = dx ! km
    ! define grid object
    call grid_init(ice_grid,name=domain,mtype=type,units="kilometers", &
      x0=real(x0,dp),y0=real(y0,dp),dx=real(dx,dp),nx=nx,dy=real(dy,dp),ny=ny, &
      lambda=real(lambda,dp),phi=real(phi,dp),alpha=real(alpha,dp),lon180=.true.)               

  return

  end subroutine ice_grid_init

end module climber_grid
