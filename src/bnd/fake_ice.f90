!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : f a k e _ i c e _ m o d
!
!  Purpose : fake ice
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
module fake_ice_mod

  use precision, only : wp, dp
  use ncio
  use control, only : ifake_ice, fake_ice_const_file, fake_ice_var_file
  use timer, only : sec_year, n_year_geo
  use coords, only : grid_class, grid_init

  implicit none

  type grid_ice_to_cmn_type
    integer, dimension(:,:), allocatable :: i_lowres 
    integer, dimension(:,:), allocatable :: j_lowres
    integer, dimension(:,:), allocatable :: ncells
  end type

  type fake_ice_type
     type(grid_class) :: grid
     ! fake ice <-> cmn grid correspondence
     type(grid_ice_to_cmn_type) :: grid_ice_to_cmn
     real(wp), dimension(:,:), allocatable :: h_ice
     real(wp), dimension(:,:), allocatable :: dh_ice_dt
   end type fake_ice_type

   real(wp), dimension(:,:), allocatable :: h_ice_0, h_ice_1
   real(wp), dimension(:), allocatable :: time_ice
   real(wp), dimension(:), allocatable :: lon_ice, lat_ice
   integer :: ntime, ni_ice, nj_ice, i0, i1, i1_old

   private
   public :: fake_ice_init, fake_ice_update, fake_ice_type


contains

  subroutine fake_ice_init(time_now,cmn_grid,fake_ice)


  implicit none

  real(wp), intent(in) :: time_now
  type(grid_class), intent(in) :: cmn_grid
  type(fake_ice_type) :: fake_ice

  integer :: i, j, imin
  integer :: ii, jj
  real(wp) :: lon1, lon2, lat1, lat2, dlon_sur, dlat_sur
  integer :: ppos, spos
  real(wp) :: w0, w1
  real(wp) :: dlon_ice, dlat_ice
  character(len=256) :: fnm


  if (ifake_ice.eq.0) then
    fnm = fake_ice_const_file
  else if (ifake_ice.eq.1 .or. ifake_ice.eq.2) then
    fnm = fake_ice_var_file
  endif

  ! read lat/lon 
  ni_ice = nc_size(trim(fnm),"lon")
  nj_ice = nc_size(trim(fnm),"lat")
  allocate( lon_ice(ni_ice) )
  allocate( lat_ice(nj_ice) )

  allocate(fake_ice%h_ice(ni_ice,nj_ice))
  allocate(fake_ice%dh_ice_dt(ni_ice,nj_ice))

  call nc_read(trim(fnm),"lon",lon_ice)
  call nc_read(trim(fnm),"lat",lat_ice)
  dlon_ice = lon_ice(2)-lon_ice(1)
  dlat_ice = lat_ice(2)-lat_ice(1)

  if (ifake_ice.eq.0) then
    ! constant ice thickness

    ! generate grid object
    spos = scan(trim(fnm),"/", BACK= .true.)+1
    ppos = scan(trim(fnm),".", BACK= .true.)-1
    call grid_init(fake_ice%grid,name=trim(fnm(spos:ppos)),mtype="latlon",units="degrees", &
      x0=real(lon_ice(1),dp),dx=real(dlon_ice,dp),nx=ni_ice,y0=real(lat_ice(1),dp),dy=real(dlat_ice,dp),ny=nj_ice)

    call nc_read(trim(fake_ice_const_file),"ice_thickness",fake_ice%h_ice)

  else if (ifake_ice.eq.1) then
    ! variable ice thickness 

    ! generate grid object
    spos = scan(trim(fnm),"/", BACK= .true.)+1
    ppos = scan(trim(fnm),".", BACK= .true.)-1
    call grid_init(fake_ice%grid,name=trim(fnm(spos:ppos)),mtype="latlon",units="degrees", &
      x0=real(lon_ice(1),dp),dx=real(dlon_ice,dp),nx=ni_ice,y0=real(lat_ice(1),dp),dy=real(dlat_ice,dp),ny=nj_ice)

    ! read time dimension
    ntime = nc_size(trim(fnm),"time")
    allocate( time_ice(ntime) )
    call nc_read(trim(fnm),"time",time_ice)

    ! look for time slice closest to current time
    if (time_now.le.time_ice(lbound(time_ice,1))) then
      i0 = 1
      i1 = 1
      w0 = 1._wp
      w1 = 0._wp
      !stop 'ERROR, fake_ice not defined for current time'
    else if (time_now.gt.time_ice(ubound(time_ice,1))) then
      i0 = ubound(time_ice,1)
      i1 = ubound(time_ice,1)
      w0 = 0._wp
      w1 = 1._wp
      !stop 'ERROR, fake_ice not defined for current time'
    else
      imin = minloc(abs(time_ice-time_now),1) 
      if (time_ice(imin).lt.time_now) then
        i0 = imin
        i1 = imin+1
      else
        i0 = imin-1
        i1 = imin
      endif
      w0 = 1._wp - dble(abs(time_ice(i0)-time_now))/dble((time_ice(i1)-time_ice(i0)))
      w1 = 1._wp - w0
    endif

    allocate(h_ice_0(ni_ice,nj_ice))
    allocate(h_ice_1(ni_ice,nj_ice))
    h_ice_0 = 0._wp
    h_ice_1 = 0._wp

    call nc_read(trim(fake_ice_var_file),"ice_thickness",h_ice_0, start=[1,1,i0], count=[ni_ice,nj_ice,1] ) 
    call nc_read(trim(fake_ice_var_file),"ice_thickness",h_ice_1, start=[1,1,i1], count=[ni_ice,nj_ice,1] ) 

    fake_ice%h_ice = w0*h_ice_0 + w1*h_ice_1

    fake_ice%dh_ice_dt = 0._wp 

  else if (ifake_ice.eq.2) then
    ! constant ice thickness closest to initial time read from variable file

    ! generate grid object
    spos = scan(trim(fnm),"/", BACK= .true.)+1
    ppos = scan(trim(fnm),".", BACK= .true.)-1
    call grid_init(fake_ice%grid,name=trim(fnm(spos:ppos)),mtype="latlon",units="degrees", &
      x0=real(lon_ice(1),dp),dx=real(dlon_ice,dp),nx=ni_ice,y0=real(lat_ice(1),dp),dy=real(dlat_ice,dp),ny=nj_ice)

    ! read time dimension
    ntime = nc_size(trim(fnm),"time")
    allocate( time_ice(ntime) )
    call nc_read(trim(fnm),"time",time_ice)

    ! look for time slice closest to current time
    if (time_now.le.time_ice(lbound(time_ice,1))) then
      i0 = 1
      !stop 'ERROR, fake_ice not defined for current time'
    else if (time_now.gt.time_ice(ubound(time_ice,1))) then
      i0 = ubound(time_ice,1)
      !stop 'ERROR, fake_ice not defined for current time'
    else
      imin = minloc(abs(time_ice-time_now),1) 
      if (time_ice(imin).lt.time_now) then
        i0 = imin
      else
        i0 = imin-1
      endif
    endif

    call nc_read(trim(fake_ice_var_file),"ice_thickness",fake_ice%h_ice, start=[1,1,i0], count=[ni_ice,nj_ice,1] ) 

  endif

  ! derive correspondence between indexes on fake geo and coupler grids
  allocate(fake_ice%grid_ice_to_cmn%i_lowres(fake_ice%grid%G%nx,fake_ice%grid%G%ny))
  allocate(fake_ice%grid_ice_to_cmn%j_lowres(fake_ice%grid%G%nx,fake_ice%grid%G%ny))
  allocate(fake_ice%grid_ice_to_cmn%ncells(cmn_grid%G%nx,cmn_grid%G%ny))
  dlon_sur = cmn_grid%lon(2,1)-cmn_grid%lon(1,1)
  dlat_sur = cmn_grid%lat(1,2)-cmn_grid%lat(1,1)
  fake_ice%grid_ice_to_cmn%ncells = 0
  do i=1,cmn_grid%G%nx
    do j=1,cmn_grid%G%ny
      lon1 = cmn_grid%lon(i,j)-0.5_wp*dlon_sur
      if (i.eq.1) lon1 = lon1-1.e-3_wp
      lon2 = cmn_grid%lon(i,j)+0.5_wp*dlon_sur
      lat1 = cmn_grid%lat(i,j)-0.5_wp*dlat_sur
      if (j.eq.1) lat1 = lat1-1.e-3_wp
      lat2 = cmn_grid%lat(i,j)+0.5_wp*dlat_sur
      do ii=1,fake_ice%grid%G%nx
        do jj=1,fake_ice%grid%G%ny 
          if (fake_ice%grid%lon(ii,jj).gt.lon1 .and. fake_ice%grid%lon(ii,jj).le.lon2 &
            .and. fake_ice%grid%lat(ii,jj).gt.lat1 .and. fake_ice%grid%lat(ii,jj).le.lat2) then
            fake_ice%grid_ice_to_cmn%i_lowres(ii,jj) = i
            fake_ice%grid_ice_to_cmn%j_lowres(ii,jj) = j
            fake_ice%grid_ice_to_cmn%ncells(i,j) = fake_ice%grid_ice_to_cmn%ncells(i,j) + 1
          endif
        enddo
      enddo
    enddo
  enddo

!  where (fake_ice%grid%lat.lt.-60._wp .or. (fake_ice%grid%lat.gt.60._wp .and. fake_ice%grid%lon.gt.-60._wp .and.fake_ice%grid%lon.lt.0._wp))
!  elsewhere
!    fake_ice%h_ice = min(100._wp,fake_ice%h_ice)
!  endwhere

  return

  end subroutine fake_ice_init


  subroutine fake_ice_update(time_now,fake_ice)

  implicit none

  real(wp), intent(in) :: time_now
  type(fake_ice_type), intent(inout) :: fake_ice

  integer :: imin
  real(wp) :: w0, w1


  if (ifake_ice.eq.1) then

    ! look for time slice closest to current time
    i1_old = i1
    if (time_now.lt.time_ice(lbound(time_ice,1))) then
      i0 = 1
      i1 = 1
      w0 = 1._wp
      w1 = 0._wp
      !stop 'ERROR, ice mask not defined for current time'
    else if (time_now.gt.time_ice(ubound(time_ice,1))) then
      w0 = 0._wp
      w1 = 1._wp
      !stop 'ERROR, ice mask not defined for current time'
    else
      imin = minloc(abs(time_ice-time_now),1) 
      if (time_ice(imin).lt.time_now) then
        i0 = imin
        i1 = imin+1
      else
        i0 = imin-1
        i1 = imin
      endif
      w0 = 1._wp - dble(abs(time_ice(i0)-time_now))/dble(time_ice(i1)-time_ice(i0))
      w1 = 1._wp - w0
    endif

    ! read new time slice if necessary
    if (i0.eq.i1_old) then
      h_ice_0 = h_ice_1
      call nc_read(trim(fake_ice_var_file),"ice_thickness",h_ice_1,start=[1,1,i1],count=[ni_ice,nj_ice,1] )
    endif

    ! rate of change of ice sheet thickness
    fake_ice%dh_ice_dt = (w0*h_ice_0 + w1*h_ice_1 - fake_ice%h_ice) / (n_year_geo*sec_year)   ! m/s

    fake_ice%h_ice = w0*h_ice_0 + w1*h_ice_1

  endif

  return

  end subroutine fake_ice_update


end module fake_ice_mod
