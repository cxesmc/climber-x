!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : f a k e _ d u s t _ m o d
!
!  Purpose : fake dust
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
module fake_dust_mod

  use precision, only : wp
  use timer, only: doy, nday_year, time_soy_bnd, monthly2daily
  use climber_grid, only :  ni, nj
  use control, only : ifake_dust, fake_dust_const_file, fake_dust_var_file
  use ncio

  implicit none

  type fake_dust_type
    real(wp), dimension(:,:), allocatable :: dust_dep
  end type
    
  real(wp), dimension(:,:,:), allocatable :: dust_dep
  real(wp), dimension(:,:,:), allocatable :: dust_dep_0
  real(wp), dimension(:,:,:), allocatable :: dust_dep_1
  real(wp), dimension(:), allocatable :: time_dust

  integer :: i0, i1, i1_old
  integer, dimension(nday_year) :: m0, m1
  real(wp), dimension(nday_year) :: wtm0, wtm1
 
  private
  public :: fake_dust_init, fake_dust_update, fake_dust_type

contains

  subroutine fake_dust_init(time_now,dust)

  implicit none

  real(wp), intent(in) :: time_now
  type(fake_dust_type) :: dust

  integer :: imin, ntime
  real(wp) :: w0, w1


  ! allocate dust type variables
  allocate(dust%dust_dep(ni,nj))

  ! allocate local variables
  allocate(dust_dep(ni,nj,0:13))
  if (ifake_dust.eq.1) then
    allocate(dust_dep_0(ni,nj,12))
    allocate(dust_dep_1(ni,nj,12))
  endif

  ! get weights for interpolation from monthly to daily
  call monthly2daily(m0,m1,wtm0,wtm1)

  ! read constant climate forcing
  if (ifake_dust.eq.0) then

    call nc_read(fake_dust_const_file,"dust_dep",dust_dep(:,:,1:12),start=[1,1,1],count=[ni,nj,12])
    where (dust_dep.lt.0._wp)
      dust_dep = 0._wp
    endwhere

  ! initilaize variable climate forcing
  else if (ifake_dust.eq.1) then

    ntime = nc_size(trim(fake_dust_var_file),"time")
    allocate( time_dust(ntime) )
    call nc_read(trim(fake_dust_var_file),"time",time_dust)
    ! look for time slice closest to current time
    if (time_now.lt.time_dust(lbound(time_dust,1))) then
      i0 = 1
      i1 = 1
      w0 = 1._wp
      w1 = 0._wp
      !stop 'ERROR, climate forcing not defined for current time'
    else if (time_now.gt.time_dust(ubound(time_dust,1))) then
      stop 'ERROR, climate forcing not defined for current time'
    else
      imin = minloc(abs(time_dust-time_now),1) 
      if (time_dust(imin).lt.time_now) then
        i0 = imin
        i1 = imin+1
      else
        i0 = imin-1
        i1 = imin
      endif
      w0 = 1._wp - abs(time_dust(i0)-time_now)/(time_dust(i1)-time_dust(i0))
      w1 = 1._wp - w0
    endif
    call nc_read(fake_dust_var_file,"dust_dep",dust_dep_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_dust_var_file,"dust_dep",dust_dep_1,start=[1,1,1,i1],count=[ni,nj,12,1] )

    dust_dep(:,:,1:12) = w0*dust_dep_0 + w1*dust_dep_1

    where (dust_dep.lt.0._wp)
      dust_dep = 0._wp
    endwhere

  endif

  dust_dep(:,:,0) = dust_dep(:,:,12)
  dust_dep(:,:,13) = dust_dep(:,:,1)

  ! initialize dust
  dust%dust_dep  = dust_dep(:,:,1) 


  return

  end subroutine fake_dust_init


  subroutine fake_dust_update(time_now,dust)

  implicit none

  real(wp), intent(in) :: time_now
  type(fake_dust_type)  :: dust

  integer :: imin
  real(wp) :: w0, w1


  ! constant climate forcing
  if (ifake_dust.eq.0) then

    dust%dust_dep  = wtm0(doy)*dust_dep(:,:,m0(doy))    + wtm1(doy)*dust_dep(:,:,m1(doy))

    ! variable climate forcing
  else if (ifake_dust.eq.1) then
    if (time_soy_bnd) then
      ! look for time slice closest to current time
      i1_old = i1
      if (time_now.lt.time_dust(lbound(time_dust,1))) then
        i0 = 1
        i1 = 1
        w0 = 1._wp
        w1 = 0._wp
        !stop 'ERROR, fake_atm not defined for current time'
      else if (time_now.gt.time_dust(ubound(time_dust,1))) then
        w0 = 0._wp
        w1 = 1._wp
        !stop 'ERROR, fake_atm not defined for current time'
      else
        imin = minloc(abs(time_dust-time_now),1) 
        if (time_dust(imin).lt.time_now) then
          i0 = imin
          i1 = imin+1
        else
          i0 = imin-1
          i1 = imin
        endif
        w0 = 1._wp - abs(time_dust(i0)-time_now)/(time_dust(i1)-time_dust(i0))
        w1 = 1._wp - w0
      endif
      ! read new time slice if necessary
      if (i0.eq.i1_old) then
        dust_dep_0 = dust_dep_1
        call nc_read(trim(fake_dust_var_file),"dust_dep",dust_dep_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
      endif
      dust_dep(:,:,1:12) = w0*dust_dep_0 + w1*dust_dep_1

      where (dust_dep.lt.0._wp)
        dust_dep = 0._wp
      endwhere

      ! extended monthly grid
      dust_dep(:,:,0) = dust_dep(:,:,12)
      dust_dep(:,:,13) = dust_dep(:,:,1)

    endif

    dust%dust_dep  = wtm0(doy)*dust_dep(:,:,m0(doy))    + wtm1(doy)*dust_dep(:,:,m1(doy))

  endif


  return

  end subroutine fake_dust_update


end module fake_dust_mod
