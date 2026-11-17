!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : f i l l _ o c e a n _ m o d 
!
!  Purpose : fill ocean domain starting from deep ocean origin point
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
module fill_ocean_mod

  use precision, only : wp

  implicit none

  integer :: ni, nj
  logical, dimension(:,:), allocatable :: mask_sea_level
  integer, dimension(:,:), allocatable :: mask_ocean

  private
  public :: fill_ocean, fill_ocean_lowres

  contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  f i l l _ o c e a n
  ! Purpose  :  fill whole ocean domain mask
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine fill_ocean(z_topo,lon,lat,lon_ocn_origin,lat_ocn_origin,mask) 

    implicit none

    real(wp), intent(in) :: z_topo(:,:)
    real(wp), intent(in) :: lon(:), lat(:)
    real(wp), intent(in) :: lon_ocn_origin, lat_ocn_origin
    integer, intent(inout) :: mask(:,:)

    integer :: i_ocn_origin, j_ocn_origin


    ni = size(z_topo,1)
    nj = size(z_topo,2)

    allocate(mask_sea_level(ni,nj))
    allocate(mask_ocean(ni,nj))

    ! potential ocean mask based on sea level (current sea level is 0 by definition)
    where (z_topo.le.0._wp) 
      mask_sea_level = .true.
    elsewhere
      mask_sea_level = .false.
    endwhere

    ! initialize
    mask_ocean = 0

    ! determine index of origin point
    i_ocn_origin = minloc(abs(lon-lon_ocn_origin),1)
    j_ocn_origin = minloc(abs(lat-lat_ocn_origin),1)

    ! recursively fill ocean domain starting from origin point
    ! all points below sea level directly connected to the origin point are ocean points
    call fill_ocean_sub(i_ocn_origin,j_ocn_origin)   

    ! set ocean 
    where (mask_ocean==1 .and. mask.ne.3)
      mask = 2  ! ocean point
    endwhere

    ! set disconnected floating ice points to grounded ice
    where (mask_ocean==0 .and. mask.eq.3)
      mask = 0  ! grounded ice
    endwhere

    deallocate(mask_sea_level)
    deallocate(mask_ocean)

    return

  end subroutine fill_ocean


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  f i l l _ o c e a n _ l o w r e s
  ! Purpose  :  determine ocean cells that are connected to origin point
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine fill_ocean_lowres(f_ocn,lon,lat,lon_ocn_origin,lat_ocn_origin, &
      mask_ocn_connect) 

    implicit none

    real(wp), intent(in) :: f_ocn(:,:)
    real(wp), intent(in) :: lon(:), lat(:)
    real(wp), intent(in) :: lon_ocn_origin, lat_ocn_origin
    integer, intent(out) :: mask_ocn_connect(:,:)

    integer :: i_ocn_origin, j_ocn_origin


    ni = size(f_ocn,1)
    nj = size(f_ocn,2)

    allocate(mask_sea_level(ni,nj))
    allocate(mask_ocean(ni,nj))

    ! ocean mask 
    where (f_ocn.gt.0._wp) 
      mask_sea_level = .true.
    elsewhere
      mask_sea_level = .false.
    endwhere

    ! initialize
    mask_ocean = 0

    ! determine index of origin point
    i_ocn_origin = minloc(abs(lon-lon_ocn_origin),1)
    j_ocn_origin = minloc(abs(lat-lat_ocn_origin),1)

    ! recursively fill ocean domain starting from origin point
    ! all ocean points directly connected to the origin point are connected ocean points
    call fill_ocean_sub(i_ocn_origin,j_ocn_origin)   

    ! set ocean 
    mask_ocn_connect = mask_ocean

    deallocate(mask_sea_level)
    deallocate(mask_ocean)

    return

  end subroutine fill_ocean_lowres

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  f i l l _ o c e a n _ s u b
  ! Purpose  :  recursive subroutine to fill whole ocean domain mask
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  recursive subroutine fill_ocean_sub(i,j) 

    implicit none

    integer, intent(in) :: i, j

    integer :: im1, ip1, jm1, jp1


    if (mask_sea_level(i,j) .and. mask_ocean(i,j)==0) then

      ! This is a contiguous cell - index it
      mask_ocean(i,j) = 1

      ! Now recursively search neighbouring cells
      ! West
      im1 = i-1
      if (im1 > 0) then
        call fill_ocean_sub(im1,j)
      endif
      ! East
      ip1 = i+1
      if (ip1 < (ni+1)) then
        call fill_ocean_sub(ip1,j)
      endif
      ! South
      jm1 = j-1
      if (jm1 > 0) then
        call fill_ocean_sub(i,jm1)
      endif
      ! North
      jp1 = j+1
      if (jp1 < (nj+1)) then
        call fill_ocean_sub(i,jp1)
      endif
!      ! North-west
!      im1 = i-1
!      jp1 = j+1
!      if (im1 > 0 .and. jp1 < (nj+1)) then
!        call fill_ocean_sub(im1,jp1)
!      endif
!      ! South-west
!      im1 = i-1
!      jm1 = j-1
!      if (im1 > 0 .and. jm1 > 0) then
!        call fill_ocean_sub(im1,jm1)
!      endif
!      ! North-east
!      ip1 = i+1
!      jp1 = j+1
!      if (ip1 < (ni+1) .and. jp1 < (nj+1)) then
!        call fill_ocean_sub(ip1,jp1)
!      endif
!      ! South-east
!      ip1 = i+1
!      jm1 = j-1
!      if (ip1 < (ni+1) .and. jm1 > 0) then
!        call fill_ocean_sub(ip1,jm1)
!      endif

      ! wrap east-west
      ! Jump west
      if (i == 1) then
        call fill_ocean_sub(ni,j)
      endif
      ! Jump east
      if (i == ni) then
        call fill_ocean_sub(1,j)
      endif

    endif

  end subroutine fill_ocean_sub

end module fill_ocean_mod
