!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : f l u x _ a d j _ m o d
!
!  Purpose : flux adjustment
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
module flux_adj_mod

  use precision, only : wp
  use climber_grid, only : lat, lon, basin_mask, i_atlantic, i_pacific, i_southern
  use ocn_grid, only : maxi, maxj, dx, dy
  use constants, only: rho_w
  use ocn_params, only: flux_adj_atl, lat_min_flux_adj_atl, lat_max_flux_adj_atl
  use ocn_params, only: flux_adj_ant, nj_flux_adj_ant
  use ocn_params, only: flux_adj_pac, lat_min_flux_adj_pac, lat_max_flux_adj_pac

  implicit none

  integer , allocatable :: j_flux_adj_atl(:,:)
  integer , allocatable :: j_flux_adj_ant(:,:)
  integer , allocatable :: j_flux_adj_pac(:,:)
  real(wp), allocatable :: rflux_adj_atl(:,:)
  real(wp), allocatable :: rflux_adj_ant(:,:)
  real(wp), allocatable :: rflux_adj_pac(:,:)

  private
  public :: flux_adj_init, flux_adj_update

contains


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  f l u x _ a d j _ i n i t
  !   Purpose    :  Atlantic-Pacific freshwater flux adjustment
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine flux_adj_init(f_ocn)

    implicit none

    real(wp), intent(in) :: f_ocn(:,:)

    integer :: i, j, cnt


    allocate(rflux_adj_atl(maxi,maxj))
    allocate(rflux_adj_ant(maxi,maxj))
    allocate(rflux_adj_pac(maxi,maxj))
    allocate(j_flux_adj_atl(maxi,2))
    allocate(j_flux_adj_ant(maxi,2))
    allocate(j_flux_adj_pac(maxi,2))

    ! find index of northern and southern boundary of Atlantic box where flux_adj_atl should be applied
    do j=2,maxj
     ! southern boundary
     if ((lat(j-1).le.lat_min_flux_adj_atl) .and. (lat(j).ge.lat_min_flux_adj_atl)) then
        j_flux_adj_atl(:,1) = j
     endif
     ! northern boundary
     if ((lat(j-1).le.lat_max_flux_adj_atl) .and. (lat(j).ge.lat_max_flux_adj_atl)) then
        j_flux_adj_atl(:,2) = j-1
     endif
    enddo

    ! find index of northern and southern boundary of Antarctic ocean belt where flux_adj_ant should be applied
    do i=1,maxi
      cnt = 0
      do j=2,maxj
        if (cnt.lt.nj_flux_adj_ant .and. f_ocn(i,j).gt.0._wp) then  ! first ocean cells from south
          cnt = cnt+1
          if (cnt.eq.1)               j_flux_adj_ant(i,1) = j
          if (cnt.eq.nj_flux_adj_ant) j_flux_adj_ant(i,2) = j
        endif
      enddo
    enddo

    ! find index of northern and southern boundary of Pacific box where flux_adj_pac should be applied
    do j=2,maxj
     ! southern boundary
     if ((lat(j-1).le.lat_min_flux_adj_pac) .and. (lat(j).ge.lat_min_flux_adj_pac)) then
        j_flux_adj_pac(:,1) = j
     endif
     ! northern boundary
     if ((lat(j-1).le.lat_max_flux_adj_pac) .and. (lat(j).ge.lat_max_flux_adj_pac)) then
        j_flux_adj_pac(:,2) = j-1
     endif
    enddo

   return

  end subroutine flux_adj_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  f l u x _ a d j _ u p d a t e
  !   Purpose    :  freshwater flux adjustment
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine flux_adj_update(f_ocn, fw_flux_adj)

    implicit none

    real(wp), intent(in) :: f_ocn(:,:)
    real(wp), intent(out) :: fw_flux_adj(:,:)

    integer :: i, j
    real(wp) :: area_flux_adj_atl
    real(wp) :: area_flux_adj_ant
    real(wp) :: area_flux_adj_pac
    real(wp) :: area_flux_comp


    area_flux_adj_atl = 0._wp
    ! compute total flux_adj Atlantic area
    do i=1,maxi
      do j=j_flux_adj_atl(i,1),j_flux_adj_atl(i,2)
        ! Atlantic
        if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_atlantic) area_flux_adj_atl = area_flux_adj_atl + dx(j)*dy*f_ocn(i,j)
      enddo
    enddo

    area_flux_adj_ant = 0._wp
    ! compute total flux_adj Antarctic ocean area
    do i=1,maxi
      do j=j_flux_adj_ant(i,1),j_flux_adj_ant(i,2)
        ! Atlantic
        if (f_ocn(i,j).gt.0._wp) area_flux_adj_ant = area_flux_adj_ant + dx(j)*dy*f_ocn(i,j)
      enddo
    enddo

    area_flux_adj_pac = 0._wp
    ! compute total flux_adj Pacific area
    do i=1,maxi
      do j=j_flux_adj_pac(i,1),j_flux_adj_pac(i,2)
        ! Atlantic
        if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_pacific) area_flux_adj_pac = area_flux_adj_pac + dx(j)*dy*f_ocn(i,j)
      enddo
    enddo

    area_flux_comp = 0._wp
    ! compute total ocean area
    do i=1,maxi
      do j=1,maxj
        if (f_ocn(i,j).gt.0._wp) area_flux_comp = area_flux_comp + dx(j)*dy*f_ocn(i,j)
      enddo
    enddo

    ! derive flux_adj weight for each cell
    rflux_adj_atl = 0._wp
    do i=1,maxi
      do j=j_flux_adj_atl(i,1),j_flux_adj_atl(i,2)
        ! Atlantic
        if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_atlantic) rflux_adj_atl(i,j) = 1.d6/area_flux_adj_atl  ! m3/s/Sv / m2 = m/s/Sv
      enddo
    enddo
    rflux_adj_ant = 0._wp
    do i=1,maxi
      do j=j_flux_adj_ant(i,1),j_flux_adj_ant(i,2)
        ! Antarctic ocean
        if (f_ocn(i,j).gt.0._wp) rflux_adj_ant(i,j) = 1.d6/area_flux_adj_ant  ! m3/s/Sv / m2 = m/s/Sv
      enddo
    enddo
    rflux_adj_pac = 0._wp
    do i=1,maxi
      do j=j_flux_adj_pac(i,1),j_flux_adj_pac(i,2)
        ! Pacific
        if (f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_pacific) rflux_adj_pac(i,j) = 1.d6/area_flux_adj_pac  ! m3/s/Sv / m2 = m/s/Sv
      enddo
    enddo

    ! compensate over the global ocean
    do i=1,maxi
      do j=1,maxj
        if (f_ocn(i,j).gt.0._wp) rflux_adj_atl(i,j) = rflux_adj_atl(i,j) - 1.d6/area_flux_comp  ! m3/s/Sv / m2 = m/s/Sv
        if (f_ocn(i,j).gt.0._wp) rflux_adj_ant(i,j) = rflux_adj_ant(i,j) - 1.d6/area_flux_comp  ! m3/s/Sv / m2 = m/s/Sv
        if (f_ocn(i,j).gt.0._wp) rflux_adj_pac(i,j) = rflux_adj_pac(i,j) - 1.d6/area_flux_comp  ! m3/s/Sv / m2 = m/s/Sv
      enddo
    enddo

    ! flux adjustment freshwater forcing
    fw_flux_adj = (flux_adj_atl*rflux_adj_atl + flux_adj_pac*rflux_adj_pac + flux_adj_ant*rflux_adj_ant) * rho_w  ! Sv * m/s/Sv * kg/m3 = kg/m2/s (freshwater density)


   return

  end subroutine flux_adj_update

end module flux_adj_mod
