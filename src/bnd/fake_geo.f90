!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : f a k e _ g e o _ m o d
!
!  Purpose : fake geo
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
module fake_geo_mod

  use precision, only : wp, dp
  use ncio
  use control, only : ifake_geo, fake_geo_const_file, fake_geo_var_file, fake_geo_ref_file
  use coords, only : grid_class, grid_init

  implicit none

  type grid_geo_to_cmn_type
    integer, dimension(:,:), allocatable :: i_lowres 
    integer, dimension(:,:), allocatable :: j_lowres
    integer, dimension(:,:), allocatable :: ncells
  end type

  type fake_geo_type
     type(grid_class) :: grid
     ! fake geo <-> cmn grid correspondence
     type(grid_geo_to_cmn_type) :: grid_geo_to_cmn
     real(wp), dimension(:,:), allocatable :: z_bed
     real(wp), dimension(:,:), allocatable :: z_bed_ref
   end type fake_geo_type

   real(wp), dimension(:,:), allocatable :: z_bed_0, z_bed_1
   real(wp), dimension(:), allocatable :: time_geo
   real(wp), dimension(:), allocatable :: lon_geo, lat_geo
   integer :: ntime, ni_geo, nj_geo, i0, i1, i1_old

   private
   public :: fake_geo_init, fake_geo_update, fake_geo_type


contains

  subroutine fake_geo_init(time_now,cmn_grid,sea_level,fake_geo)


  implicit none

  real(wp), intent(in) :: time_now
  type(grid_class), intent(in) :: cmn_grid
  real(wp), intent(in) :: sea_level
  type(fake_geo_type) :: fake_geo

  integer :: imin
  integer :: ppos, spos
  integer :: i, j, ii, jj
  integer :: ni_geo_ref, nj_geo_ref
  real(wp) :: lon1, lon2, lat1, lat2, dlon_sur, dlat_sur
  real(wp) :: w0, w1
  real(wp) :: dlon_geo, dlat_geo
  character(len=256) :: fnm


  if (ifake_geo.eq.0) then
    fnm = fake_geo_const_file
  else if (ifake_geo.eq.1 .or. ifake_geo.eq.2) then
    fnm = fake_geo_var_file
  endif

  ! read lat/lon 
  ni_geo = nc_size(trim(fnm),"lon")
  nj_geo = nc_size(trim(fnm),"lat")
  allocate( lon_geo(ni_geo) )
  allocate( lat_geo(nj_geo) )

  allocate(fake_geo%z_bed(ni_geo,nj_geo))

  call nc_read(trim(fnm),"lon",lon_geo)
  call nc_read(trim(fnm),"lat",lat_geo)
  dlon_geo = lon_geo(2)-lon_geo(1)
  dlat_geo = lat_geo(2)-lat_geo(1)

  ! read reference bedrock topography file 
  ! used to compute anomalies to be added on top of the high-resolution reference geo bedrock topography
  ni_geo_ref = nc_size(trim(fake_geo_ref_file),"lon")
  nj_geo_ref = nc_size(trim(fake_geo_ref_file),"lat")
  if (ni_geo_ref.ne.ni_geo .or. nj_geo_ref.ne.nj_geo) then
    print *,'WARNING: dimensions of fake_geo_ref_file does not match those in fake_geo_const/var_file'
    print *,'This warning can be safely ignored if fake_geo_ref_file==geo_ref_file'
  endif
  allocate(fake_geo%z_bed_ref(ni_geo_ref,nj_geo_ref))
  call nc_read(trim(fake_geo_ref_file),"bedrock_topography",fake_geo%z_bed_ref)

  if (ifake_geo.eq.0) then
    ! constant topography

    ! generate grid object
    spos = scan(trim(fnm),"/", BACK= .true.)+1
    ppos = scan(trim(fnm),".", BACK= .true.)-1
    call grid_init(fake_geo%grid,name=trim(fnm(spos:ppos)),mtype="latlon",units="degrees", &
      x0=real(lon_geo(1),dp),dx=real(dlon_geo,dp),nx=ni_geo,y0=real(lat_geo(1),dp),dy=real(dlat_geo,dp),ny=nj_geo)

    call nc_read(trim(fake_geo_const_file),"bedrock_topography",fake_geo%z_bed)

    ! account for prescribed sea level
    fake_geo%z_bed = fake_geo%z_bed - sea_level

  else if (ifake_geo.eq.1) then
    ! variable bedrock elevation from file

    ! generate grid object
    spos = scan(trim(fnm),"/", BACK= .true.)+1
    ppos = scan(trim(fnm),".", BACK= .true.)-1
    call grid_init(fake_geo%grid,name=trim(fnm(spos:ppos)),mtype="latlon",units="degrees", &
      x0=real(lon_geo(1),dp),dx=real(dlon_geo,dp),nx=ni_geo,y0=real(lat_geo(1),dp),dy=real(dlat_geo,dp),ny=nj_geo)

    ! read time dimension
    ntime = nc_size(trim(fnm),"time")
    allocate( time_geo(ntime) )
    call nc_read(trim(fnm),"time",time_geo)

    ! look for time slice closest to current time
    if (time_now.le.time_geo(lbound(time_geo,1))) then
      i0 = 1
      i1 = 1
      w0 = 1._wp
      w1 = 0._wp
      !stop 'ERROR, fake_geo not defined for current time'
    else if (time_now.gt.time_geo(ubound(time_geo,1))) then
      i0 = ubound(time_geo,1)
      i1 = ubound(time_geo,1)
      w0 = 0._wp
      w1 = 1._wp
      !stop 'ERROR, fake_geo not defined for current time'
    else
      imin = minloc(abs(time_geo-time_now),1) 
      if (time_geo(imin).lt.time_now) then
        i0 = imin
        i1 = imin+1
      else
        i0 = imin-1
        i1 = imin
      endif
      w0 = 1._wp - dble(abs(time_geo(i0)-time_now))/dble((time_geo(i1)-time_geo(i0)))
      w1 = 1._wp - w0
    endif

    ! bedrock elevation
    allocate(z_bed_0(ni_geo,nj_geo))
    allocate(z_bed_1(ni_geo,nj_geo))
    z_bed_0 = 0._wp
    z_bed_1 = 0._wp

    call nc_read(trim(fake_geo_var_file),"bedrock_topography",z_bed_0, start=[1,1,i0], count=[ni_geo,nj_geo,1] ) ! bedrock elevation
    call nc_read(trim(fake_geo_var_file),"bedrock_topography",z_bed_1, start=[1,1,i1], count=[ni_geo,nj_geo,1] ) ! bedrock elevation

    fake_geo%z_bed = w0*z_bed_0 + w1*z_bed_1

  else if (ifake_geo.eq.2) then
    ! constant bedrock topography closest to initial time read from variable file

    ! generate grid object
    spos = scan(trim(fnm),"/", BACK= .true.)+1
    ppos = scan(trim(fnm),".", BACK= .true.)-1
    call grid_init(fake_geo%grid,name=trim(fnm(spos:ppos)),mtype="latlon",units="degrees", &
      x0=real(lon_geo(1),dp),dx=real(dlon_geo,dp),nx=ni_geo,y0=real(lat_geo(1),dp),dy=real(dlat_geo,dp),ny=nj_geo)

    ! read time dimension
    ntime = nc_size(trim(fnm),"time")
    allocate( time_geo(ntime) )
    call nc_read(trim(fnm),"time",time_geo)

    ! look for time geo closest to current time
    if (time_now.le.time_geo(lbound(time_geo,1))) then
      i0 = 1
      !stop 'ERROR, fake_geo not defined for current time'
    else if (time_now.gt.time_geo(ubound(time_geo,1))) then
      i0 = ubound(time_geo,1)
      !stop 'ERROR, fake_geo not defined for current time'
    else
      imin = minloc(abs(time_geo-time_now),1) 
      if (time_geo(imin).lt.time_now) then
        i0 = imin
      else
        i0 = imin-1
      endif
    endif

    ! bedrock elevation
    call nc_read(trim(fake_geo_var_file),"bedrock_topography",fake_geo%z_bed, start=[1,1,i0], count=[ni_geo,nj_geo,1] ) ! bedrock elevation

  endif

  ! derive correspondence between indexes on fake geo and coupler grids
  allocate(fake_geo%grid_geo_to_cmn%i_lowres(fake_geo%grid%G%nx,fake_geo%grid%G%ny))
  allocate(fake_geo%grid_geo_to_cmn%j_lowres(fake_geo%grid%G%nx,fake_geo%grid%G%ny))
  allocate(fake_geo%grid_geo_to_cmn%ncells(cmn_grid%G%nx,cmn_grid%G%ny))
  dlon_sur = cmn_grid%lon(2,1)-cmn_grid%lon(1,1)
  dlat_sur = cmn_grid%lat(1,2)-cmn_grid%lat(1,1)
  fake_geo%grid_geo_to_cmn%ncells = 0
  do i=1,cmn_grid%G%nx
    do j=1,cmn_grid%G%ny
      lon1 = cmn_grid%lon(i,j)-0.5_wp*dlon_sur
      if (i.eq.1) lon1 = lon1-1.e-3_wp
      lon2 = cmn_grid%lon(i,j)+0.5_wp*dlon_sur
      lat1 = cmn_grid%lat(i,j)-0.5_wp*dlat_sur
      if (j.eq.1) lat1 = lat1-1.e-3_wp
      lat2 = cmn_grid%lat(i,j)+0.5_wp*dlat_sur
      do ii=1,fake_geo%grid%G%nx
        do jj=1,fake_geo%grid%G%ny 
          if (fake_geo%grid%lon(ii,jj).gt.lon1 .and. fake_geo%grid%lon(ii,jj).le.lon2 &
            .and. fake_geo%grid%lat(ii,jj).gt.lat1 .and. fake_geo%grid%lat(ii,jj).le.lat2) then
            fake_geo%grid_geo_to_cmn%i_lowres(ii,jj) = i
            fake_geo%grid_geo_to_cmn%j_lowres(ii,jj) = j
            fake_geo%grid_geo_to_cmn%ncells(i,j) = fake_geo%grid_geo_to_cmn%ncells(i,j) + 1
          endif
        enddo
      enddo
    enddo
  enddo


  return

  end subroutine fake_geo_init


  subroutine fake_geo_update(time_now,fake_geo)

  implicit none

  real(wp), intent(in) :: time_now
  type(fake_geo_type), intent(inout) :: fake_geo

  integer :: imin
  real(wp) :: w0, w1


  if (ifake_geo.eq.1) then

    ! look for time slice closest to current time
    i1_old = i1
    if (time_now.lt.time_geo(lbound(time_geo,1))) then
      i0 = 1
      i1 = 1
      w0 = 1._wp
      w1 = 0._wp
      !stop 'ERROR, geo mask not defined for current time'
    else if (time_now.gt.time_geo(ubound(time_geo,1))) then
      w0 = 0._wp
      w1 = 1._wp
      !stop 'ERROR, geo mask not defined for current time'
    else
      imin = minloc(abs(time_geo-time_now),1) 
      if (time_geo(imin).lt.time_now) then
        i0 = imin
        i1 = imin+1
      else
        i0 = imin-1
        i1 = imin
      endif
      w0 = 1._wp - dble(abs(time_geo(i0)-time_now))/dble(time_geo(i1)-time_geo(i0))
      w1 = 1._wp - w0
    endif

    ! read new time slice if necessary
    if (i0.eq.i1_old) then
      z_bed_0 = z_bed_1
      call nc_read(trim(fake_geo_var_file),"bedrock_topography",z_bed_1,start=[1,1,i1],count=[ni_geo,nj_geo,1] ) ! bedrock elevation
    endif

    fake_geo%z_bed  = w0*z_bed_0  + w1*z_bed_1

  endif


  return

  end subroutine fake_geo_update


end module fake_geo_mod
