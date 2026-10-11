!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : f a k e _ s i c _ m o d
!
!  Purpose : fake sea ice
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
module fake_sic_mod

  use precision, only : wp
  use timer, only: doy, nday_year, time_soy_bnd, monthly2daily
  use climber_grid, only :  ni, nj
  use control, only : ifake_sic, fake_sic_const_file, fake_sic_var_file, prc_forcing
  use ncio

  implicit none

  type fake_sic_type
    real(wp), dimension(:,:), allocatable :: f_sic
  end type
    
  real(wp), dimension(:,:,:), allocatable :: f_sic, f_sic_0, f_sic_1
  real(wp), dimension(:), allocatable :: time_sic

  integer :: i0, i1, i1_old
  integer, dimension(nday_year) :: m0, m1
  real(wp), dimension(nday_year) :: wtm0, wtm1
 
  private
  public :: fake_sic_init, fake_sic_update, fake_sic_type

contains

  subroutine fake_sic_init(time_now,sic)

  implicit none

  real(wp), intent(in) :: time_now
  type(fake_sic_type) :: sic

  integer :: imin, ntime
  real(wp) :: w0, w1


  ! allocate sic type variables
  allocate(sic%f_sic(ni,nj))

  ! allocate local variables
  allocate(f_sic(ni,nj,0:13))
  if (ifake_sic.eq.1) then
    allocate(f_sic_0(ni,nj,12))
    allocate(f_sic_1(ni,nj,12))
  endif

  ! get weights for interpolation from monthly to daily
  call monthly2daily(m0,m1,wtm0,wtm1)

  ! read constant climate forcing
  if (ifake_sic.eq.0) then

    call nc_read(fake_sic_const_file,"f_sic",f_sic(:,:,1:12),start=[1,1,1],count=[ni,nj,12])
    where (f_sic(:,:,1:12).gt.0.95_wp) 
      f_sic(:,:,1:12) = 1._wp
    endwhere

  ! initilaize variable climate forcing
  else if (ifake_sic.eq.1) then

    ntime = nc_size(trim(fake_sic_var_file),"time")
    allocate( time_sic(ntime) )
    call nc_read(trim(fake_sic_var_file),"time",time_sic)
    ! look for time slice closest to current time
    if (time_now.lt.time_sic(lbound(time_sic,1))) then
      i0 = 1
      i1 = 1
      w0 = 1._wp
      w1 = 0._wp
      !stop 'ERROR, sea ice forcing not defined for current time'
    else if (time_now.gt.time_sic(ubound(time_sic,1))) then
      stop 'ERROR, sea ice forcing not defined for current time'
    else
      imin = minloc(abs(time_sic-time_now),1) 
      if (time_sic(imin).lt.time_now) then
        i0 = imin
        i1 = imin+1
      else
        i0 = imin-1
        i1 = imin
      endif
      w0 = 1._wp - abs(time_sic(i0)-time_now)/(time_sic(i1)-time_sic(i0))
      w1 = 1._wp - w0
    endif
    call nc_read(fake_sic_var_file,"f_sic",f_sic_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_sic_var_file,"f_sic",f_sic_1,start=[1,1,1,i1],count=[ni,nj,12,1] )

    f_sic(:,:,1:12) = w0*f_sic_0 + w1*f_sic_1

  endif

  f_sic(:,:,0)  = f_sic(:,:,12)
  f_sic(:,:,13) = f_sic(:,:,1)


  return

  end subroutine fake_sic_init


  subroutine fake_sic_update(time_now,sic)

  implicit none

  real(wp), intent(in) :: time_now
  type(fake_sic_type)  :: sic

  integer :: imin
  real(wp) :: w0, w1


  ! constant climate forcing
  if (ifake_sic.eq.0) then

    sic%f_sic    = wtm0(doy)*f_sic(:,:,m0(doy)) + wtm1(doy)*f_sic(:,:,m1(doy)) ! N/m2

    ! variable climate forcing
  else if (ifake_sic.eq.1) then
    if (time_soy_bnd) then
      ! look for time slice closest to current time
      i1_old = i1
      if (time_now.lt.time_sic(lbound(time_sic,1))) then
        i0 = 1
        i1 = 1
        w0 = 1._wp
        w1 = 0._wp
        !stop 'ERROR, fake_sic not defined for current time'
      else if (time_now.gt.time_sic(ubound(time_sic,1))) then
        stop 'ERROR, fake_sic not defined for current time'
      else
        imin = minloc(abs(time_sic-time_now),1) 
        if (time_sic(imin).lt.time_now) then
          i0 = imin
          i1 = imin+1
        else
          i0 = imin-1
          i1 = imin
        endif
        w0 = 1._wp - abs(time_sic(i0)-time_now)/(time_sic(i1)-time_sic(i0))
        w1 = 1._wp - w0
      endif
      ! read new time slice if necessary
      if (i0.eq.i1_old) then
        f_sic_0 = f_sic_1
        call nc_read(trim(fake_sic_var_file),"f_sic",f_sic_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
      endif
      f_sic(:,:,1:12) = w0*f_sic_0 + w1*f_sic_1

      ! extended monthly grid
      f_sic(:,:,0)  = f_sic(:,:,12)
      f_sic(:,:,13) = f_sic(:,:,1)

    endif

    sic%f_sic = wtm0(doy)*f_sic(:,:,m0(doy)) + wtm1(doy)*f_sic(:,:,m1(doy))

  endif


  return

  end subroutine fake_sic_update


end module fake_sic_mod
