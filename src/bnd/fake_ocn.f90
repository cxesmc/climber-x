!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : f a k e _ o c n _ m o d
!
!  Purpose : fake ocean
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
module fake_ocn_mod

  use precision, only : wp
  use constants, only : T0
  use timer, only : doy, nday_year, time_soy_bnd, monthly2daily
  use climber_grid, only : ni, nj
  use control, only : ifake_ocn, fake_ocn_const_file, fake_ocn_var_file
  use ncio

  implicit none

  type fake_ocn_type
    real(wp), dimension(:), allocatable :: depth
    real(wp), dimension(:,:), allocatable :: sst, sss, t1l, s1l, t_shelf, uo1, vo1
    real(wp), dimension(:,:,:), allocatable :: t, s
  end type fake_ocn_type

  real(wp), dimension(:), allocatable :: time_ocn
  real(wp), dimension(:,:,:), allocatable :: sst, sss, t1l, s1l, t_shelf, uo1, vo1
  real(wp), dimension(:,:,:), allocatable :: sst_0, sss_0, t1l_0, s1l_0, t_shelf_0, uo1_0, vo1_0
  real(wp), dimension(:,:,:), allocatable :: sst_1, sss_1, t1l_1, s1l_1, t_shelf_1, uo1_1, vo1_1
  real(wp), dimension(:,:,:,:), allocatable :: t, s
  real(wp), dimension(:,:,:,:), allocatable :: t_0, s_0
  real(wp), dimension(:,:,:,:), allocatable :: t_1, s_1

  integer :: nk, k_s, k_1l
  integer :: i0, i1, i1_old
  integer, dimension(nday_year) :: m0, m1
  real(wp), dimension(nday_year) :: wtm0, wtm1
  
  private
  public :: fake_ocn_init, fake_ocn_update, fake_ocn_type

contains

  subroutine fake_ocn_init(time_now,ocn)

  implicit none

  real(wp), intent(in) :: time_now
  type(fake_ocn_type) :: ocn

  integer :: imin, ntime
  real(wp) :: w0, w1


  ! allocate ocn type variables
  allocate(ocn%sst(ni,nj))
  allocate(ocn%sss(ni,nj))
  allocate(ocn%t1l(ni,nj))
  allocate(ocn%s1l(ni,nj))
  allocate(ocn%t_shelf(ni,nj))
  allocate(ocn%uo1(ni,nj))
  allocate(ocn%vo1(ni,nj))

  ! allocate local variables
  allocate(sst(ni,nj,0:13))
  if (ifake_ocn.eq.1) then
    allocate(sst_0(ni,nj,12))
    allocate(sst_1(ni,nj,12))
  endif
  allocate(sss(ni,nj,0:13))
  if (ifake_ocn.eq.1) then
    allocate(sss_0(ni,nj,12))
    allocate(sss_1(ni,nj,12))
  endif
  allocate(t1l(ni,nj,0:13))
  if (ifake_ocn.eq.1) then
    allocate(t1l_0(ni,nj,12))
    allocate(t1l_1(ni,nj,12))
  endif
  allocate(s1l(ni,nj,0:13))
  if (ifake_ocn.eq.1) then
    allocate(s1l_0(ni,nj,12))
    allocate(s1l_1(ni,nj,12))
  endif
  allocate(t_shelf(ni,nj,0:13))
  if (ifake_ocn.eq.1) then
    allocate(t_shelf_0(ni,nj,12))
    allocate(t_shelf_1(ni,nj,12))
  endif
  allocate(uo1(ni,nj,0:13))
  if (ifake_ocn.eq.1) then
    allocate(uo1_0(ni,nj,12))
    allocate(uo1_1(ni,nj,12))
  endif
  allocate(vo1(ni,nj,0:13))
  if (ifake_ocn.eq.1) then
    allocate(vo1_0(ni,nj,12))
    allocate(vo1_1(ni,nj,12))
  endif

  ! get weights for interpolation from monthly to daily
  call monthly2daily(m0,m1,wtm0,wtm1)

  k_s  = 1  ! index of surface layer
  k_1l = 1  ! index of first layer

  if (ifake_ocn.eq.0) then
    ! climatological ocean 

    nk = nc_size(trim(fake_ocn_const_file),"depth")
    allocate(ocn%t(ni,nj,nk))
    allocate(ocn%s(ni,nj,nk))
    allocate(ocn%depth(nk))
    allocate(t(ni,nj,nk,0:13))
    if (ifake_ocn.eq.1) then
      allocate(t_0(ni,nj,nk,12))
      allocate(t_1(ni,nj,nk,12))
    endif
    allocate(s(ni,nj,nk,0:13))
    if (ifake_ocn.eq.1) then
      allocate(s_0(ni,nj,nk,12))
      allocate(s_1(ni,nj,nk,12))
    endif
    call nc_read(fake_ocn_const_file,"depth",ocn%depth(:) )

    call nc_read(fake_ocn_const_file,"t",sst(:,:,1:12), start=[1,1,k_s,1], count=[ni,nj,1,12] )
    call nc_read(fake_ocn_const_file,"s",sss(:,:,1:12), start=[1,1,k_s,1], count=[ni,nj,1,12] )
    call nc_read(fake_ocn_const_file,"t",t1l(:,:,1:12), start=[1,1,k_1l,1], count=[ni,nj,1,12] )
    call nc_read(fake_ocn_const_file,"s",s1l(:,:,1:12), start=[1,1,k_1l,1], count=[ni,nj,1,12] )
    call nc_read(fake_ocn_const_file,"t",t_shelf(:,:,1:12), start=[1,1,k_1l,1], count=[ni,nj,1,12] )
    call nc_read(fake_ocn_const_file,"t",t(:,:,:,1:12), start=[1,1,1,1], count=[ni,nj,nk,12] )
    call nc_read(fake_ocn_const_file,"s",s(:,:,:,1:12), start=[1,1,1,1], count=[ni,nj,nk,12] )
    call nc_read(fake_ocn_const_file,"u_sic",uo1(:,:,1:12),start=[1,1,1],count=[ni,nj,12])
    call nc_read(fake_ocn_const_file,"v_sic",vo1(:,:,1:12),start=[1,1,1],count=[ni,nj,12])

  else if (ifake_ocn.eq.1) then
    ! time dependent ocean 

    ntime = nc_size(trim(fake_ocn_var_file),"time")
    allocate( time_ocn(ntime) )
    call nc_read(trim(fake_ocn_var_file),"time",time_ocn)
    ! look for time slice closest to current time
    if (time_now.lt.time_ocn(lbound(time_ocn,1))) then
      i0 = 1
      i1 = 1
      w0 = 1._wp
      w1 = 0._wp
      !stop 'ERROR, ocean forcing not defined for current time'
    else if (time_now.gt.time_ocn(ubound(time_ocn,1))) then
      stop 'ERROR, ocean forcing not defined for current time'
    else
      imin = minloc(abs(time_ocn-time_now),1) 
      if (time_ocn(imin).lt.time_now) then
        i0 = imin
        i1 = imin+1
      else
        i0 = imin-1
        i1 = imin
      endif
      w0 = 1._wp - abs(time_ocn(i0)-time_now)/(time_ocn(i1)-time_ocn(i0))
      w1 = 1._wp - w0
    endif

    nk = nc_size(trim(fake_ocn_var_file),"depth")
    allocate(ocn%t(ni,nj,nk))
    allocate(ocn%s(ni,nj,nk))
    allocate(ocn%depth(nk))
    if (ifake_ocn.eq.1) then
      allocate(t_0(ni,nj,nk,12))
      allocate(t_1(ni,nj,nk,12))
    endif
    allocate(s(ni,nj,nk,0:13))
    if (ifake_ocn.eq.1) then
      allocate(s_0(ni,nj,nk,12))
      allocate(s_1(ni,nj,nk,12))
    endif
    call nc_read(fake_ocn_var_file,"depth",ocn%depth(:) )

    call nc_read(fake_ocn_var_file,"t",sst_0,start=[1,1,k_s,1,i0],count=[ni,nj,1,12,1] )
    call nc_read(fake_ocn_var_file,"t",sst_1,start=[1,1,k_s,1,i1],count=[ni,nj,1,12,1] )
    sst(:,:,1:12) = w0*sst_0 + w1*sst_1

    call nc_read(fake_ocn_var_file,"s",sst_0,start=[1,1,k_s,1,i0],count=[ni,nj,1,12,1] )
    call nc_read(fake_ocn_var_file,"s",sst_1,start=[1,1,k_s,1,i1],count=[ni,nj,1,12,1] )
    sss(:,:,1:12) = w0*sss_0 + w1*sss_1

    call nc_read(fake_ocn_var_file,"t",t1l_0,start=[1,1,k_1l,1,i0],count=[ni,nj,1,12,1] )
    call nc_read(fake_ocn_var_file,"t",t1l_1,start=[1,1,k_1l,1,i1],count=[ni,nj,1,12,1] )
    t1l(:,:,1:12) = w0*t1l_0 + w1*t1l_1

    call nc_read(fake_ocn_var_file,"s",s1l_0,start=[1,1,k_1l,1,i0],count=[ni,nj,1,12,1] )
    call nc_read(fake_ocn_var_file,"s",s1l_1,start=[1,1,k_1l,1,i1],count=[ni,nj,1,12,1] )
    s1l(:,:,1:12) = w0*s1l_0 + w1*s1l_1

    call nc_read(fake_ocn_var_file,"t",t_shelf_0,start=[1,1,k_1l,1,i0],count=[ni,nj,1,12,1] )
    call nc_read(fake_ocn_var_file,"t",t_shelf_1,start=[1,1,k_1l,1,i1],count=[ni,nj,1,12,1] )
    t_shelf(:,:,1:12) = w0*t_shelf_0 + w1*t_shelf_1

    call nc_read(fake_ocn_var_file,"t",t_0,start=[1,1,1,1,i0],count=[ni,nj,nk,12,1] )
    call nc_read(fake_ocn_var_file,"t",t_1,start=[1,1,1,1,i1],count=[ni,nj,nk,12,1] )
    t(:,:,:,1:12) = w0*t_0 + w1*t_1

    call nc_read(fake_ocn_var_file,"s",s_0,start=[1,1,1,1,i0],count=[ni,nj,nk,12,1] )
    call nc_read(fake_ocn_var_file,"s",s_1,start=[1,1,1,1,i1],count=[ni,nj,nk,12,1] )
    s(:,:,:,1:12) = w0*s_0 + w1*s_1

    call nc_read(fake_ocn_var_file,"u_sic",uo1_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_ocn_var_file,"u_sic",uo1_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    uo1(:,:,1:12) = w0*uo1_0 + w1*uo1_1

    call nc_read(fake_ocn_var_file,"v_sic",vo1_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_ocn_var_file,"v_sic",vo1_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    vo1(:,:,1:12) = w0*vo1_0 + w1*vo1_1

  endif

  ! periodic boundary conditions

  sst(:,:,0)  = sst(:,:,12)
  sst(:,:,13) = sst(:,:,1)
  ocn%sst     = sst(:,:,1) 

  sss(:,:,0)  = sss(:,:,12)
  sss(:,:,13) = sss(:,:,1)
  ocn%sss     = sss(:,:,1) 

  t1l(:,:,0)  = t1l(:,:,12)
  t1l(:,:,13) = t1l(:,:,1)
  ocn%t1l     = t1l(:,:,1) 

  s1l(:,:,0)  = s1l(:,:,12)
  s1l(:,:,13) = s1l(:,:,1)
  ocn%s1l     = s1l(:,:,1) 

  t_shelf(:,:,0)  = t_shelf(:,:,12)
  t_shelf(:,:,13) = t_shelf(:,:,1)
  ocn%t_shelf     = t_shelf(:,:,1) 

  t(:,:,:,0)  = t(:,:,:,12)
  t(:,:,:,13) = t(:,:,:,1)
  ocn%t       = t(:,:,:,1) 

  s(:,:,:,0)  = s(:,:,:,12)
  s(:,:,:,13) = s(:,:,:,1)
  ocn%s       = s(:,:,:,1) 

  uo1(:,:,0)  = uo1(:,:,12)
  uo1(:,:,13) = uo1(:,:,1)
  ocn%uo1     = uo1(:,:,1) 

  vo1(:,:,0)  = vo1(:,:,12)
  vo1(:,:,13) = vo1(:,:,1)
  ocn%vo1     = vo1(:,:,1) 
 


  return

  end subroutine fake_ocn_init


  subroutine fake_ocn_update(time_now,ocn)

  implicit none

  real(wp), intent(in) :: time_now
  type(fake_ocn_type), intent(inout) :: ocn

  integer :: imin
  real(wp) :: w0, w1


  if (ifake_ocn.eq.0) then

    ocn%sst = wtm0(doy)*sst(:,:,m0(doy)) + wtm1(doy)*sst(:,:,m1(doy))
    ocn%sss = wtm0(doy)*sss(:,:,m0(doy)) + wtm1(doy)*sss(:,:,m1(doy))
    ocn%t1l = wtm0(doy)*t1l(:,:,m0(doy)) + wtm1(doy)*t1l(:,:,m1(doy))
    ocn%s1l = wtm0(doy)*s1l(:,:,m0(doy)) + wtm1(doy)*s1l(:,:,m1(doy))
    ocn%t_shelf = wtm0(doy)*t_shelf(:,:,m0(doy)) + wtm1(doy)*t_shelf(:,:,m1(doy)) + T0
    ocn%t = wtm0(doy)*t(:,:,:,m0(doy)) + wtm1(doy)*t(:,:,:,m1(doy)) + T0
    ocn%s = wtm0(doy)*s(:,:,:,m0(doy)) + wtm1(doy)*s(:,:,:,m1(doy)) 
    ocn%uo1 = wtm0(doy)*uo1(:,:,m0(doy)) + wtm1(doy)*uo1(:,:,m1(doy))
    ocn%vo1 = wtm0(doy)*vo1(:,:,m0(doy)) + wtm1(doy)*vo1(:,:,m1(doy))

    ! variable climate forcing
  else if (ifake_ocn.eq.1) then

    if (time_soy_bnd) then
      ! look for time slice closest to current time
      i1_old = i1
      if (time_now.lt.time_ocn(lbound(time_ocn,1))) then
        i0 = 1
        i1 = 1
        w0 = 1._wp
        w1 = 0._wp
        !stop 'ERROR, ocean forcing not defined for current time'
      else if (time_now.gt.time_ocn(ubound(time_ocn,1))) then
        w0 = 0._wp
        w1 = 1._wp
        !stop 'ERROR, ocean forcing not defined for current time'
      else
        imin = minloc(abs(time_ocn-time_now),1) 
        if (time_ocn(imin).lt.time_now) then
          i0 = imin
          i1 = imin+1
        else
          i0 = imin-1
          i1 = imin
        endif
        w0 = 1._wp - abs(time_ocn(i0)-time_now)/(time_ocn(i1)-time_ocn(i0))
        w1 = 1._wp - w0
      endif
      ! read new time slice if necessary
      if (i0.eq.i1_old) then
        sst_0 = sst_1
        call nc_read(trim(fake_ocn_var_file),"t",sst_1,start=[1,1,k_s,1,i1],count=[ni,nj,1,12,1] )
        sss_0 = sss_1
        call nc_read(trim(fake_ocn_var_file),"s",sss_1,start=[1,1,k_s,1,i1],count=[ni,nj,1,12,1] )
        t1l_0 = t1l_1
        call nc_read(trim(fake_ocn_var_file),"t",t1l_1,start=[1,1,k_1l,1,i1],count=[ni,nj,1,12,1] )
        s1l_0 = s1l_1
        call nc_read(trim(fake_ocn_var_file),"s",s1l_1,start=[1,1,k_1l,1,i1],count=[ni,nj,1,12,1] )
        t_shelf_0 = t_shelf_1
        call nc_read(trim(fake_ocn_var_file),"t",t_shelf_1,start=[1,1,k_1l,1,i1],count=[ni,nj,1,12,1] )
        t_0 = t_1
        call nc_read(trim(fake_ocn_var_file),"t",t_1,start=[1,1,1,1,i1],count=[ni,nj,nk,12,1] )
        s_0 = s_1
        call nc_read(trim(fake_ocn_var_file),"s",s_1,start=[1,1,1,1,i1],count=[ni,nj,nk,12,1] )
        uo1_0 = uo1_1
        call nc_read(trim(fake_ocn_var_file),"u_sic",uo1_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        vo1_0 = vo1_1
        call nc_read(trim(fake_ocn_var_file),"v_sic",vo1_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
      endif
      sst(:,:,1:12) = w0*sst_0 + w1*sst_1
      sst(:,:,0)    = sst(:,:,12)
      sst(:,:,13)   = sst(:,:,1)
      sss(:,:,1:12) = w0*sss_0 + w1*sss_1
      sss(:,:,0)    = sss(:,:,12)
      sss(:,:,13)   = sss(:,:,1)
      t1l(:,:,1:12) = w0*t1l_0 + w1*t1l_1
      t1l(:,:,0)    = t1l(:,:,12)
      t1l(:,:,13)   = t1l(:,:,1)
      s1l(:,:,1:12) = w0*s1l_0 + w1*s1l_1
      s1l(:,:,0)    = s1l(:,:,12)
      s1l(:,:,13)   = s1l(:,:,1)
      t_shelf(:,:,1:12) = w0*t_shelf_0 + w1*t_shelf_1
      t_shelf(:,:,0)    = t_shelf(:,:,12)
      t_shelf(:,:,13)   = t_shelf(:,:,1)
      t(:,:,:,1:12) = w0*t_0 + w1*t_1
      t(:,:,:,0)    = t(:,:,:,12)
      t(:,:,:,13)   = t(:,:,:,1)
      s(:,:,:,1:12) = w0*s_0 + w1*s_1
      s(:,:,:,0)    = s(:,:,:,12)
      s(:,:,:,13)   = s(:,:,:,1)
      uo1(:,:,1:12) = w0*uo1_0 + w1*uo1_1
      uo1(:,:,0)    = uo1(:,:,12)
      uo1(:,:,13)   = uo1(:,:,1)
      vo1(:,:,1:12) = w0*vo1_0 + w1*vo1_1
      vo1(:,:,0)    = vo1(:,:,12)
      vo1(:,:,13)   = vo1(:,:,1)
    endif

    ocn%sst = wtm0(doy)*sst(:,:,m0(doy)) + wtm1(doy)*sst(:,:,m1(doy))
    ocn%sss = wtm0(doy)*sss(:,:,m0(doy)) + wtm1(doy)*sss(:,:,m1(doy))
    ocn%t1l = wtm0(doy)*t1l(:,:,m0(doy)) + wtm1(doy)*t1l(:,:,m1(doy))
    ocn%s1l = wtm0(doy)*s1l(:,:,m0(doy)) + wtm1(doy)*s1l(:,:,m1(doy))
    ocn%t_shelf = wtm0(doy)*t_shelf(:,:,m0(doy)) + wtm1(doy)*t_shelf(:,:,m1(doy)) + T0
    ocn%t= wtm0(doy)*t(:,:,:,m0(doy)) + wtm1(doy)*t(:,:,:,m1(doy)) + T0
    ocn%s= wtm0(doy)*s(:,:,:,m0(doy)) + wtm1(doy)*s(:,:,:,m1(doy)) 
    ocn%uo1 = wtm0(doy)*uo1(:,:,m0(doy)) + wtm1(doy)*uo1(:,:,m1(doy))
    ocn%vo1 = wtm0(doy)*vo1(:,:,m0(doy)) + wtm1(doy)*vo1(:,:,m1(doy))

  endif


  return

  end subroutine fake_ocn_update


end module fake_ocn_mod
