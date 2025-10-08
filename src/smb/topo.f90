!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : t o p o _ m o d
!
!  Purpose : topography-related functions for SMB model
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2024 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit and Reinhard Calov
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
module topo_mod

  use precision, only : wp
  use constants, only : pi, r_earth
  use coord, only : grid_class, grid_allocate, map_class, map_field, map_scrip_class, map_scrip_field
  use smb_params, only : h_atm, p0, prc_par, surf_par, gamma

  implicit none

  private
  public :: topo_filter, topo_grad_map, topo_factors

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  t o p o _ f i l t e r
  !   Purpose    :  filter topography
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine topo_filter(grid, z_sur, z_sur_fil)

  implicit none

  type(grid_class), intent(in) :: grid
  real(wp), dimension(:,:), intent(in) :: z_sur
  real(wp), dimension(:,:), intent(out) :: z_sur_fil

  integer :: i, j, ii, jj, i_f, j_f, nx, ny, n_filter
  real(wp) :: dx, dy
  real(wp) :: sigma_filter
  real(wp) :: dist, weigh, sum_weigh


  nx = grid%G%nx
  ny = grid%G%ny

  dx = grid%G%dx ! km
  dy = grid%G%dy ! km

  sigma_filter = prc_par%topo_filter_width/grid%G%dx   ! half span of filtered area, in grid points
  n_filter     = ceiling(2.0_wp*sigma_filter)

  do i=1,nx 
    do j=1,ny

      sum_weigh = 0.0_wp
      z_sur_fil(i,j) = 0._wp

      do ii=-n_filter, n_filter
        do jj=-n_filter, n_filter

          i_f = i+ii
          j_f = j+jj

          if (i_f <  1) i_f = 1
          if (i_f > nx) i_f = nx

          if (j_f <  1) j_f = 1
          if (j_f > ny) j_f = ny

          dist      = sqrt(real(ii,wp)**2+real(jj,wp)**2)
          weigh     = exp(-(dist/sigma_filter)**2)
          sum_weigh = sum_weigh + weigh

          z_sur_fil(i,j) = z_sur_fil(i,j) + weigh*z_sur(i_f,j_f)

        end do
      end do

      z_sur_fil(i,j) = z_sur_fil(i,j)/sum_weigh

    end do
  end do

  return

  end subroutine topo_filter


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  t o p o _ g r a d _ m a p 
  !   Purpose    :  derive topography gradients
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine topo_grad_map(grid_latlon, maps_to_latlon, maps_from_latlon, z_sur, &
                      dz_sur, dz_dx_sur, dz_dy_sur)

  implicit none

  type(grid_class), intent(in) :: grid_latlon
  type(map_scrip_class), intent(in) :: maps_to_latlon
  type(map_scrip_class), intent(in) :: maps_from_latlon
  real(wp), dimension(:,:), intent(in) :: z_sur
  real(wp), dimension(:,:), intent(out) :: dz_sur, dz_dx_sur, dz_dy_sur

  real(wp), dimension(:,:), allocatable :: z_sur_latlon
  real(wp), dimension(:,:), allocatable :: dz_dx_sur_latlon
  real(wp), dimension(:,:), allocatable :: dz_dy_sur_latlon

  integer :: i, j, im, ip, jm, jp, nx, ny


  ! ++++++++++++++++++++++++
  ! topography gradients

  ! map surface elevation to lat-lon
  call grid_allocate(grid_latlon, z_sur_latlon)
  call grid_allocate(grid_latlon, dz_dx_sur_latlon)
  call grid_allocate(grid_latlon, dz_dy_sur_latlon)
  call map_scrip_field(maps_to_latlon,"z_sur",z_sur,z_sur_latlon,method="mean")

  ! compute topography gradients in spherical coordinates
  nx = grid_latlon%G%nx
  ny = grid_latlon%G%ny

  do j=1,ny
    do i=1,nx

      im=i-1
      if (im.lt.1) im=1
      ip=i+1
      if (ip.gt.nx) ip=nx
      jm=j-1
      if (jm.lt.1) jm=1
      jp=j+1
      if (jp.gt.ny) jp=ny

      dz_dx_sur_latlon(i,j) = (1._wp/(r_earth*cos(pi*grid_latlon%lat(i,j)/180._wp))) &
        * (z_sur_latlon(ip,j)-z_sur_latlon(im,j))/(pi/180._wp*(grid_latlon%lon(ip,j)-grid_latlon%lon(im,j)))
      dz_dy_sur_latlon(i,j) = 1._wp/r_earth &
        * (z_sur_latlon(i,jp)-z_sur_latlon(i,jm))/(pi/180._wp*(grid_latlon%lat(i,jp)-grid_latlon%lat(i,jm)))

    enddo
  enddo

  ! map topography gradients back to stereographic projection
  call map_scrip_field(maps_from_latlon,"dz_dx_sur",dz_dx_sur_latlon,dz_dx_sur,method="mean")
  call map_scrip_field(maps_from_latlon,"dz_dy_sur",dz_dy_sur_latlon,dz_dy_sur,method="mean")

  dz_sur = sqrt(dz_dx_sur**2+dz_dy_sur**2)

  deallocate(z_sur_latlon)
  deallocate(dz_dx_sur_latlon)
  deallocate(dz_dy_sur_latlon)

  return

  end subroutine topo_grad_map


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  t o p o _ f a c t o r s
  !   Purpose    :  derive topography factors
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine topo_factors(grid, z_sur, z_sur_i, &
                      f_ele, pressure)

  implicit none

  type(grid_class), intent(in) :: grid
  real(wp), dimension(:,:), intent(in) :: z_sur, z_sur_i
  real(wp), dimension(:,:), intent(out) :: f_ele, pressure

  integer :: i, j, nx, ny
  real(wp) :: dT, dz


  nx = grid%G%nx
  ny = grid%G%ny

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! elevation correction factor for precipitation using Clausius-Clapeyron

  do j=1,ny
    do i=1,nx
      if (prc_par%l_elevation_corr .and. z_sur(i,j).ge.prc_par%z_sur_crit_fele) then
        dz = z_sur(i,j)-z_sur_i(i,j)
        dT = -gamma*dz 
        f_ele(i,j) = exp(prc_par%dP_dT*dT)
      else 
        f_ele(i,j) = 1._wp
      endif
      ! additionally reduce precip above some high elevation
      if (z_sur(i,j).ge.prc_par%z_sur_high_fele) then
        f_ele(i,j) = 0.1_wp*f_ele(i,j)
      endif
    enddo
  enddo

  ! ++++++++++++++++++++++++
  ! surface pressure (elevation dependence only)

  pressure = p0*exp(-z_sur/h_atm)


  return

end subroutine topo_factors


end module topo_mod
