!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : f a k e _ a t m _ m o d
!
!  Purpose : fake atmosphere
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
module fake_atm_mod

  use precision, only : wp
  use timer, only: doy, mon, nday_year, time_soy_bnd, monthly2daily
  use climber_grid, only :  ni, nj
  use control, only : ifake_atm, fake_atm_const_file, fake_atm_var_file, flag_smb, prc_forcing
  use ncio

  implicit none

  type fake_atm_type
    real(wp), dimension(:,:), allocatable :: tair, qair, rain, snow, prc, cld, cld_day, pressure, slp, tsl, htrop, rbatm, swtoa, swdown, swdown_c, swdown_d, lwdown, usur, vsur, wind, taux, tauy
    real(wp), dimension(:,:), allocatable :: tair_min_mon
  end type
    
  real(wp), dimension(:,:,:), allocatable :: tair, qair, rain, snow, prc, cld, cld_day, pressure, slp, tsl, htrop, rbatm, swtoa, swdown, swdown_c, swdown_d, lwdown, usur, vsur, wind, taux, tauy
  real(wp), dimension(:,:), allocatable :: tair_min_mon
  real(wp), dimension(:,:,:), allocatable :: tair_0, qair_0, rain_0, snow_0, prc_0, cld_0, cld_day_0, pressure_0, slp_0, tsl_0, htrop_0, rbatm_0, swtoa_0, swdown_0, swdown_c_0, swdown_d_0, lwdown_0, usur_0, vsur_0, wind_0, taux_0, tauy_0
  real(wp), dimension(:,:,:), allocatable :: tair_1, qair_1, rain_1, snow_1, prc_1, cld_1, cld_day_1, pressure_1, slp_1, tsl_1, htrop_1, rbatm_1, swtoa_1, swdown_1, swdown_c_1, swdown_d_1, lwdown_1, usur_1, vsur_1, wind_1, taux_1, tauy_1
  real(wp), dimension(:), allocatable :: time_atm

  integer :: i0, i1, i1_old
  integer, dimension(nday_year) :: m0, m1
  real(wp), dimension(nday_year) :: wtm0, wtm1
 
  private
  public :: fake_atm_init, fake_atm_update, fake_atm_type

contains

  subroutine fake_atm_init(time_now,atm)

  implicit none

  real(wp), intent(in) :: time_now
  type(fake_atm_type) :: atm

  integer :: imin, ntime
  real(wp) :: w0, w1


  ! allocate atm type variables
  allocate(atm%taux(ni,nj))
  allocate(atm%tauy(ni,nj))
  allocate(atm%tair(ni,nj))
  allocate(atm%tair_min_mon(ni,nj))
  allocate(atm%qair(ni,nj))
  allocate(atm%prc(ni,nj))
  allocate(atm%cld(ni,nj))
  allocate(atm%cld_day(ni,nj))
  allocate(atm%rain(ni,nj))
  allocate(atm%snow(ni,nj))
  allocate(atm%pressure(ni,nj))
  allocate(atm%slp(ni,nj))
  allocate(atm%tsl(ni,nj))
  allocate(atm%htrop(ni,nj))
  allocate(atm%rbatm(ni,nj))
  allocate(atm%swtoa(ni,nj))
  allocate(atm%swdown(ni,nj))
  allocate(atm%swdown_c(ni,nj))
  allocate(atm%swdown_d(ni,nj))
  allocate(atm%lwdown(ni,nj))
  allocate(atm%usur(ni,nj))
  allocate(atm%vsur(ni,nj))
  allocate(atm%wind(ni,nj))

  ! allocate local variables
  allocate(tair(ni,nj,0:13))
  allocate(tair_min_mon(ni,nj))  
  allocate(qair(ni,nj,0:13))
  allocate(rain(ni,nj,0:13))
  allocate(snow(ni,nj,0:13))
  allocate(prc(ni,nj,0:13))
  allocate(cld(ni,nj,0:13))
  allocate(cld_day(ni,nj,0:13))
  allocate(pressure(ni,nj,0:13))
  allocate(slp(ni,nj,0:13))
  allocate(tsl(ni,nj,0:13))
  allocate(htrop(ni,nj,0:13))
  allocate(rbatm(ni,nj,0:13))
  allocate(swtoa(ni,nj,0:13))
  allocate(swdown(ni,nj,0:13))
  allocate(swdown_c(ni,nj,0:13))
  allocate(swdown_d(ni,nj,0:13))
  allocate(lwdown(ni,nj,0:13))
  allocate(usur(ni,nj,0:13))
  allocate(vsur(ni,nj,0:13))
  allocate(wind(ni,nj,0:13))
  allocate(taux(ni,nj,0:13))
  allocate(tauy(ni,nj,0:13))
  if (ifake_atm.eq.1) then
    allocate(tair_0(ni,nj,12))
    allocate(qair_0(ni,nj,12))
    allocate(rain_0(ni,nj,12))
    allocate(snow_0(ni,nj,12))
    allocate(prc_0(ni,nj,12))
    allocate(pressure_0(ni,nj,12))
    allocate(slp_0(ni,nj,12))
    allocate(tsl_0(ni,nj,12))
    allocate(htrop_0(ni,nj,12))
    allocate(cld_0(ni,nj,12))
    allocate(swdown_0(ni,nj,12))
    allocate(swdown_c_0(ni,nj,12))
    allocate(swdown_d_0(ni,nj,12))
    allocate(lwdown_0(ni,nj,12))
    allocate(wind_0(ni,nj,12))
    allocate(usur_0(ni,nj,12))
    allocate(vsur_0(ni,nj,12))
    allocate(taux_0(ni,nj,12))
    allocate(tauy_0(ni,nj,12))
    allocate(tair_1(ni,nj,12))
    allocate(qair_1(ni,nj,12))
    allocate(rain_1(ni,nj,12))
    allocate(snow_1(ni,nj,12))
    allocate(prc_1(ni,nj,12))
    allocate(pressure_1(ni,nj,12))
    allocate(slp_1(ni,nj,12))
    allocate(tsl_1(ni,nj,12))
    allocate(htrop_1(ni,nj,12))
    allocate(cld_1(ni,nj,12))
    allocate(swdown_1(ni,nj,12))
    allocate(swdown_c_1(ni,nj,12))
    allocate(swdown_d_1(ni,nj,12))
    allocate(lwdown_1(ni,nj,12))
    allocate(usur_1(ni,nj,12))
    allocate(vsur_1(ni,nj,12))
    allocate(wind_1(ni,nj,12))
    allocate(taux_1(ni,nj,12))
    allocate(tauy_1(ni,nj,12))
  endif

  ! get weights for interpolation from monthly to daily
  call monthly2daily(m0,m1,wtm0,wtm1)

  ! read constant climate forcing
  if (ifake_atm.eq.0) then

    call nc_read(fake_atm_const_file,"tair",tair(:,:,1:12) )
    call nc_read(fake_atm_const_file,"qair",qair(:,:,1:12) )
    if (prc_forcing.eq.0) then
      ! read rain and snow separately
      call nc_read(fake_atm_const_file,"rain",rain(:,:,1:12) )
      call nc_read(fake_atm_const_file,"snow",snow(:,:,1:12) )
    else if (prc_forcing.eq.1) then
      call nc_read(fake_atm_const_file,"prc",prc(:,:,1:12) )
    endif
    call nc_read(fake_atm_const_file,"cld",cld(:,:,1:12) )
    call nc_read(fake_atm_const_file,"cld_day",cld_day(:,:,1:12) )
    call nc_read(fake_atm_const_file,"pressure",pressure(:,:,1:12) )
    call nc_read(fake_atm_const_file,"slp",slp(:,:,1:12) )
    call nc_read(fake_atm_const_file,"tsl",tsl(:,:,1:12) )
    call nc_read(fake_atm_const_file,"htrop",htrop(:,:,1:12) )
    call nc_read(fake_atm_const_file,"rbatm",rbatm(:,:,1:12) )
    call nc_read(fake_atm_const_file,"swdown",swdown(:,:,1:12) )
    if (flag_smb) then
      call nc_read(fake_atm_const_file,"swdown_c",swdown_c(:,:,1:12) )
    else
      swdown_c(:,:,1:12) = swdown(:,:,1:12)
    endif
    call nc_read(fake_atm_const_file,"lwdown",lwdown(:,:,1:12) )
    call nc_read(fake_atm_const_file,"u10",usur(:,:,1:12) )
    call nc_read(fake_atm_const_file,"v10",vsur(:,:,1:12) )
    call nc_read(fake_atm_const_file,"wind",wind(:,:,1:12) )
    call nc_read(fake_atm_const_file,"taux",taux(:,:,1:12) )
    call nc_read(fake_atm_const_file,"tauy",tauy(:,:,1:12) )

    tair_min_mon = minval(tair(:,:,1:12),3)

  ! initilaize variable climate forcing
  else if (ifake_atm.eq.1) then

    ntime = nc_size(trim(fake_atm_var_file),"time")
    allocate( time_atm(ntime) )
    call nc_read(trim(fake_atm_var_file),"time",time_atm)
    ! look for time slice closest to current time
    if (time_now.lt.time_atm(lbound(time_atm,1))) then
      i0 = 1
      i1 = 1
      w0 = 1._wp
      w1 = 0._wp
      !stop 'ERROR, climate forcing not defined for current time'
    else if (time_now.gt.time_atm(ubound(time_atm,1))) then
      stop 'ERROR, climate forcing not defined for current time'
    else
      imin = minloc(abs(time_atm-time_now),1) 
      if (time_atm(imin).lt.time_now) then
        i0 = imin
        i1 = imin+1
      else
        i0 = imin-1
        i1 = imin
      endif
      w0 = 1._wp - abs(time_atm(i0)-time_now)/(time_atm(i1)-time_atm(i0))
      w1 = 1._wp - w0
    endif
    call nc_read(fake_atm_var_file,"tair",tair_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"qair",qair_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    if (prc_forcing.eq.0) then
      ! read rain and snow separately
      call nc_read(fake_atm_var_file,"rain",rain_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
      call nc_read(fake_atm_var_file,"snow",snow_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    else if (prc_forcing.eq.1) then
      call nc_read(fake_atm_var_file,"prc",prc_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    endif
    call nc_read(fake_atm_var_file,"pressure",pressure_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"slp",slp_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"tsl",tsl_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"htrop",htrop_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"cld",cld_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"swdown",swdown_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    if (flag_smb) then
      call nc_read(fake_atm_var_file,"swdown_c",swdown_c_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    else
      swdown_c_0 = swdown_0
    endif
    call nc_read(fake_atm_var_file,"lwdown",lwdown_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"u10",usur_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"v10",vsur_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"wind",wind_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"taux",taux_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"tauy",tauy_0,start=[1,1,1,i0],count=[ni,nj,12,1] )

    call nc_read(fake_atm_var_file,"tair",tair_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"qair",qair_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    if (prc_forcing.eq.0) then
      ! read rain and snow separately
      call nc_read(fake_atm_var_file,"rain",rain_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
      call nc_read(fake_atm_var_file,"snow",snow_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    else if (prc_forcing.eq.1) then
      call nc_read(fake_atm_var_file,"prc",prc_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    endif
    call nc_read(fake_atm_var_file,"pressure",pressure_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"slp",slp_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"tsl",tsl_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"htrop",htrop_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"cld",cld_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"swdown",swdown_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    if (flag_smb) then
      call nc_read(fake_atm_var_file,"swdown_c",swdown_c_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    else
      swdown_c_1 = swdown_1
    endif
    call nc_read(fake_atm_var_file,"lwdown",lwdown_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"u10",usur_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"v10",vsur_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"wind",wind_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"taux",taux_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_var_file,"tauy",tauy_1,start=[1,1,1,i1],count=[ni,nj,12,1] )

    tair(:,:,1:12) = w0*tair_0 + w1*tair_1
    qair(:,:,1:12) = w0*qair_0 + w1*qair_1
    if (prc_forcing.eq.0) then
      rain(:,:,1:12) = w0*rain_0 + w1*rain_1
      snow(:,:,1:12) = w0*snow_0 + w1*snow_1
    else if (prc_forcing.eq.1) then
      prc(:,:,1:12) = w0*prc_0 + w1*prc_1
    endif
    usur(:,:,1:12) = w0*usur_0 + w1*usur_1
    vsur(:,:,1:12) = w0*vsur_0 + w1*vsur_1
    wind(:,:,1:12) = w0*wind_0 + w1*wind_1
    taux(:,:,1:12) = w0*taux_0 + w1*taux_1
    tauy(:,:,1:12) = w0*tauy_0 + w1*tauy_1
    pressure(:,:,1:12) = w0*pressure_0 + w1*pressure_1
    slp(:,:,1:12) = w0*slp_0 + w1*slp_1
    tsl(:,:,1:12) = w0*tsl_0 + w1*tsl_1
    htrop(:,:,1:12) = w0*htrop_0 + w1*htrop_1
    cld(:,:,1:12) = w0*cld_0 + w1*cld_1
    lwdown(:,:,1:12) = w0*lwdown_0 + w1*lwdown_1
    swdown(:,:,1:12) = w0*swdown_0 + w1*swdown_1
    if (flag_smb) then
      swdown_c(:,:,1:12) = w0*swdown_c_0 + w1*swdown_c_1
    else
      swdown_c(:,:,1:12) = swdown(:,:,1:12)
    endif

    tair_min_mon = minval(tair(:,:,1:12),3)

  endif

  tair(:,:,0)  = tair(:,:,12)
  qair(:,:,0) = qair(:,:,12)
  cld(:,:,0)  = cld(:,:,12)
  cld_day(:,:,0)  = cld_day(:,:,12)
  if (prc_forcing.eq.0) then
    rain(:,:,0)  = rain(:,:,12)
    snow(:,:,0) = snow(:,:,12)
  else if (prc_forcing.eq.1) then
    prc(:,:,0)  = prc(:,:,12)
  endif
  pressure(:,:,0)  = pressure(:,:,12)
  slp(:,:,0)  = slp(:,:,12)
  tsl(:,:,0)  = tsl(:,:,12)
  htrop(:,:,0)  = htrop(:,:,12)
  rbatm(:,:,0) = rbatm(:,:,12)
  swdown(:,:,0) = swdown(:,:,12)
  swdown_c(:,:,0) = swdown_c(:,:,12)
  lwdown(:,:,0)  = lwdown(:,:,12)
  usur(:,:,0) = usur(:,:,12)
  vsur(:,:,0) = vsur(:,:,12)
  wind(:,:,0) = wind(:,:,12)
  taux(:,:,0)  = taux(:,:,12)
  tauy(:,:,0) = tauy(:,:,12)

  tair(:,:,13)  = tair(:,:,1)
  qair(:,:,13) = qair(:,:,1)
  cld(:,:,13)  = cld(:,:,1)
  cld_day(:,:,13)  = cld_day(:,:,1)
  if (prc_forcing.eq.0) then
    rain(:,:,13)  = rain(:,:,1)
    snow(:,:,13) = snow(:,:,1)
  else if (prc_forcing.eq.1) then
    prc(:,:,13)  = prc(:,:,1)
  endif
  pressure(:,:,13)  = pressure(:,:,1)
  slp(:,:,13)  = slp(:,:,1)
  tsl(:,:,13)  = tsl(:,:,1)
  htrop(:,:,13)  = htrop(:,:,1)
  rbatm(:,:,13) = rbatm(:,:,1)
  swdown(:,:,13) = swdown(:,:,1)
  swdown_c(:,:,13) = swdown_c(:,:,1)
  lwdown(:,:,13)  = lwdown(:,:,1)
  usur(:,:,13) = usur(:,:,1)
  vsur(:,:,13) = vsur(:,:,1)
  wind(:,:,13) = wind(:,:,1)
  taux(:,:,13)  = taux(:,:,1)
  tauy(:,:,13) = tauy(:,:,1)

  ! derive diffuse downward surface radiation
  where (cld.gt.0.1_wp)
    swdown_d = (swdown-(1._wp-cld)*swdown_c) / cld
  elsewhere
    swdown_d = swdown
  endwhere

  where (wind .le. 0._wp)
    wind = 1.e-3_wp
  endwhere
  where (taux.lt.-900._wp)
    taux = 0._wp
  endwhere
  where (tauy.lt.-900._wp)
    tauy = 0._wp
  endwhere

  ! initialize atm
  atm%taux     = taux(:,:,1) 
  atm%tauy     = tauy(:,:,1) 
  if (prc_forcing.eq.0) then
    atm%rain     = rain(:,:,1) 
    atm%snow     = snow(:,:,1) 
  else if (prc_forcing.eq.1) then
    atm%prc      = prc(:,:,1)  
  endif
  atm%cld      = cld(:,:,1)   
  atm%cld_day  = cld_day(:,:,1)   
  atm%pressure = pressure(:,:,1) 
  atm%slp = slp(:,:,1) 
  atm%tsl = tsl(:,:,1) 
  atm%htrop = htrop(:,:,1) 
  atm%rbatm = rbatm(:,:,1)
  atm%lwdown  = lwdown(:,:,1) 
  atm%swdown  = swdown(:,:,1)
  atm%swdown_c  = swdown_c(:,:,1)
  atm%swdown_d  = swdown_d(:,:,1)
  atm%usur     = usur(:,:,1)    
  atm%vsur     = vsur(:,:,1)   
  atm%wind     = wind(:,:,1) 
  atm%tair  = tair(:,:,1)   
  atm%qair  = qair(:,:,1)   
  atm%tair_min_mon = tair_min_mon


  return

  end subroutine fake_atm_init


  subroutine fake_atm_update(time_now,atm)

  implicit none

  real(wp), intent(in) :: time_now
  type(fake_atm_type)  :: atm

  integer :: imin
  real(wp) :: w0, w1


  ! constant climate forcing
  if (ifake_atm.eq.0) then

    atm%taux     = wtm0(doy)*taux(:,:,m0(doy))     + wtm1(doy)*taux(:,:,m1(doy)) ! N/m2
    atm%tauy     = wtm0(doy)*tauy(:,:,m0(doy))     + wtm1(doy)*tauy(:,:,m1(doy)) ! N/m2
    if (prc_forcing.eq.0) then
      atm%rain     = wtm0(doy)*rain(:,:,m0(doy))     + wtm1(doy)*rain(:,:,m1(doy))
      atm%snow     = wtm0(doy)*snow(:,:,m0(doy))     + wtm1(doy)*snow(:,:,m1(doy))
    else if (prc_forcing.eq.1) then
      atm%prc      = wtm0(doy)*prc(:,:,m0(doy))     + wtm1(doy)*prc(:,:,m1(doy))
    endif
    atm%cld      = wtm0(doy)*cld(:,:,m0(doy))     + wtm1(doy)*cld(:,:,m1(doy))
    atm%cld_day  = wtm0(doy)*cld_day(:,:,m0(doy)) + wtm1(doy)*cld_day(:,:,m1(doy))
    atm%pressure = wtm0(doy)*pressure(:,:,m0(doy)) + wtm1(doy)*pressure(:,:,m1(doy))
    atm%slp = wtm0(doy)*slp(:,:,m0(doy)) + wtm1(doy)*slp(:,:,m1(doy))
    atm%tsl = wtm0(doy)*tsl(:,:,m0(doy)) + wtm1(doy)*tsl(:,:,m1(doy))
    atm%htrop = wtm0(doy)*htrop(:,:,m0(doy)) + wtm1(doy)*htrop(:,:,m1(doy))
    atm%rbatm   = wtm0(doy)*rbatm(:,:,m0(doy))    + wtm1(doy)*rbatm(:,:,m1(doy))
    atm%lwdown  = wtm0(doy)*lwdown(:,:,m0(doy))   + wtm1(doy)*lwdown(:,:,m1(doy))
    atm%swdown  = wtm0(doy)*swdown(:,:,m0(doy))   + wtm1(doy)*swdown(:,:,m1(doy))
    atm%swdown_c  = wtm0(doy)*swdown_c(:,:,m0(doy))   + wtm1(doy)*swdown_c(:,:,m1(doy))
    atm%swdown_d  = wtm0(doy)*swdown_d(:,:,m0(doy))   + wtm1(doy)*swdown_d(:,:,m1(doy))
    atm%usur     = wtm0(doy)*usur(:,:,m0(doy))     + wtm1(doy)*usur(:,:,m1(doy))
    atm%vsur     = wtm0(doy)*vsur(:,:,m0(doy))     + wtm1(doy)*vsur(:,:,m1(doy))
    atm%wind     = wtm0(doy)*wind(:,:,m0(doy))     + wtm1(doy)*wind(:,:,m1(doy))
    !atm%tair  = wtm0(doy)*tair(:,:,m0(doy))    + wtm1(doy)*tair(:,:,m1(doy))
    atm%tair  = tair(:,:,mon)
    atm%qair  = wtm0(doy)*qair(:,:,m0(doy))    + wtm1(doy)*qair(:,:,m1(doy))
    atm%tair_min_mon = tair_min_mon

    ! variable climate forcing
  else if (ifake_atm.eq.1) then

    if (time_soy_bnd) then
      ! look for time slice closest to current time
      i1_old = i1
      if (time_now.lt.time_atm(lbound(time_atm,1))) then
        i0 = 1
        i1 = 1
        w0 = 1._wp
        w1 = 0._wp
        !stop 'ERROR, fake_atm not defined for current time'
      else if (time_now.gt.time_atm(ubound(time_atm,1))) then
        w0 = 0._wp
        w1 = 1._wp
        !stop 'ERROR, fake_atm not defined for current time'
      else
        imin = minloc(abs(time_atm-time_now),1) 
        if (time_atm(imin).lt.time_now) then
          i0 = imin
          i1 = imin+1
        else
          i0 = imin-1
          i1 = imin
        endif
        w0 = 1._wp - abs(time_atm(i0)-time_now)/(time_atm(i1)-time_atm(i0))
        w1 = 1._wp - w0
      endif
      ! read new time slice if necessary
      if (i0.eq.i1_old) then
        tair_0 = tair_1
        qair_0 = qair_1
        if (prc_forcing.eq.0) then
          rain_0 = rain_1
          snow_0 = snow_1
        else if (prc_forcing.eq.1) then
          prc_0 = prc_1
        endif
        wind_0 = wind_1
        usur_0 = usur_1
        vsur_0 = vsur_1
        pressure_0 = pressure_1
        slp_0 = slp_1
        tsl_0 = tsl_1
        htrop_0 = htrop_1
        cld_0 = cld_1
        lwdown_0 = lwdown_1
        swdown_0 = swdown_1
        swdown_c_0 = swdown_c_1
        call nc_read(trim(fake_atm_var_file),"tair",tair_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_var_file),"qair",qair_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        if (prc_forcing.eq.0) then
          call nc_read(trim(fake_atm_var_file),"rain",rain_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
          call nc_read(trim(fake_atm_var_file),"snow",snow_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
          prc_1(:,:,1:12) = rain_1(:,:,1:12) + snow_1(:,:,1:12)
        else if (prc_forcing.eq.1) then
          call nc_read(trim(fake_atm_var_file),"prc",prc_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        endif
        call nc_read(trim(fake_atm_var_file),"u10",wind_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_var_file),"v10",wind_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_var_file),"wind",wind_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_var_file),"taux",taux_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_var_file),"tauy",tauy_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_var_file),"pressure",pressure_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_var_file),"slp",slp_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_var_file),"tsl",tsl_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_var_file),"htrop",htrop_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_var_file),"cld",cld_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_var_file),"lwdown",lwdown_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_var_file),"swdown",swdown_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        if (flag_smb) then
          call nc_read(fake_atm_var_file,"swdown_c",swdown_c_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        else
          swdown_c_1 = swdown_1
        endif
      endif
      tair(:,:,1:12) = w0*tair_0 + w1*tair_1
      qair(:,:,1:12) = w0*qair_0 + w1*qair_1
      if (prc_forcing.eq.0) then
        rain(:,:,1:12) = w0*rain_0 + w1*rain_1
        snow(:,:,1:12) = w0*snow_0 + w1*snow_1
      else if (prc_forcing.eq.1) then
        prc(:,:,1:12) = w0*prc_0 + w1*prc_1
      endif
      usur(:,:,1:12)    = w0*usur_0 + w1*usur_1
      vsur(:,:,1:12)    = w0*vsur_0 + w1*vsur_1
      wind(:,:,1:12) = w0*wind_0 + w1*wind_1
      taux(:,:,1:12) = w0*taux_0 + w1*taux_1
      tauy(:,:,1:12) = w0*tauy_0 + w1*tauy_1
      pressure(:,:,1:12) = w0*pressure_0 + w1*pressure_1
      slp(:,:,1:12) = w0*slp_0 + w1*slp_1
      tsl(:,:,1:12) = w0*tsl_0 + w1*tsl_1
      htrop(:,:,1:12) = w0*htrop_0 + w1*htrop_1
      cld(:,:,1:12) = w0*cld_0 + w1*cld_1
      lwdown(:,:,1:12) = w0*lwdown_0 + w1*lwdown_1
      swdown(:,:,1:12) = w0*swdown_0 + w1*swdown_1
      swdown_c(:,:,1:12) = w0*swdown_c_0 + w1*swdown_c_1

      atm%tair_min_mon = tair_min_mon

      ! extended monthly grid
      tair(:,:,0)  = tair(:,:,12)
      qair(:,:,0) = qair(:,:,12)
      cld(:,:,0)  = cld(:,:,12)
      if (prc_forcing.eq.0) then
        rain(:,:,0)  = rain(:,:,12)
        snow(:,:,0) = snow(:,:,12)
      else if (prc_forcing.eq.1) then
        prc(:,:,0)  = prc(:,:,12)
      endif
      pressure(:,:,0)  = pressure(:,:,12)
      slp(:,:,0)  = slp(:,:,12)
      tsl(:,:,0)  = tsl(:,:,12)
      htrop(:,:,0)  = htrop(:,:,12)
      swdown(:,:,0) = swdown(:,:,12)
      swdown_c(:,:,0) = swdown_c(:,:,12)
      lwdown(:,:,0)  = lwdown(:,:,12)
      usur(:,:,0) = usur(:,:,12)
      vsur(:,:,0) = vsur(:,:,12)
      wind(:,:,0) = wind(:,:,12)
      taux(:,:,0)  = taux(:,:,12)
      tauy(:,:,0) = tauy(:,:,12)

      tair(:,:,13)  = tair(:,:,1)
      qair(:,:,13) = qair(:,:,1)
      cld(:,:,13)  = cld(:,:,1)
      if (prc_forcing.eq.0) then
        rain(:,:,13)  = rain(:,:,1)
        snow(:,:,13) = snow(:,:,1)
      else if (prc_forcing.eq.1) then
        prc(:,:,13)  = prc(:,:,1)
      endif
      pressure(:,:,13)  = pressure(:,:,1)
      slp(:,:,13)  = slp(:,:,1)
      tsl(:,:,13)  = tsl(:,:,1)
      htrop(:,:,13)  = htrop(:,:,1)
      swdown(:,:,13) = swdown(:,:,1)
      swdown_c(:,:,13) = swdown_c(:,:,1)
      lwdown(:,:,13)  = lwdown(:,:,1)
      usur(:,:,13) = usur(:,:,1)
      vsur(:,:,13) = vsur(:,:,1)
      wind(:,:,13) = wind(:,:,1)
      taux(:,:,13)  = taux(:,:,1)
      tauy(:,:,13) = tauy(:,:,1)

    endif

    atm%tair  = wtm0(doy)*tair(:,:,m0(doy))    + wtm1(doy)*tair(:,:,m1(doy))
    atm%qair  = wtm0(doy)*qair(:,:,m0(doy))    + wtm1(doy)*qair(:,:,m1(doy))
    if (prc_forcing.eq.0) then
      atm%rain  = wtm0(doy)*rain(:,:,m0(doy))    + wtm1(doy)*rain(:,:,m1(doy))
      atm%snow  = wtm0(doy)*snow(:,:,m0(doy))    + wtm1(doy)*snow(:,:,m1(doy))
    else if (prc_forcing.eq.1) then
      atm%prc   = wtm0(doy)*prc(:,:,m0(doy))     + wtm1(doy)*prc(:,:,m1(doy))
    endif
    atm%usur  = wtm0(doy)*usur(:,:,m0(doy))    + wtm1(doy)*usur(:,:,m1(doy))
    atm%vsur  = wtm0(doy)*vsur(:,:,m0(doy))    + wtm1(doy)*vsur(:,:,m1(doy))
    atm%wind  = wtm0(doy)*wind(:,:,m0(doy))    + wtm1(doy)*wind(:,:,m1(doy))
    atm%taux  = wtm0(doy)*taux(:,:,m0(doy))    + wtm1(doy)*taux(:,:,m1(doy))
    atm%tauy  = wtm0(doy)*tauy(:,:,m0(doy))    + wtm1(doy)*tauy(:,:,m1(doy))
    atm%pressure  = wtm0(doy)*pressure(:,:,m0(doy))    + wtm1(doy)*pressure(:,:,m1(doy))
    atm%slp   = wtm0(doy)*slp(:,:,m0(doy))     + wtm1(doy)*slp(:,:,m1(doy))
    atm%tsl   = wtm0(doy)*tsl(:,:,m0(doy))     + wtm1(doy)*tsl(:,:,m1(doy))
    atm%htrop   = wtm0(doy)*htrop(:,:,m0(doy))     + wtm1(doy)*htrop(:,:,m1(doy))
    atm%cld   = wtm0(doy)*cld(:,:,m0(doy))     + wtm1(doy)*cld(:,:,m1(doy))
    atm%lwdown  = wtm0(doy)*lwdown(:,:,m0(doy))    + wtm1(doy)*lwdown(:,:,m1(doy))
    atm%swdown  = wtm0(doy)*swdown(:,:,m0(doy))    + wtm1(doy)*swdown(:,:,m1(doy))
    atm%swdown_c  = wtm0(doy)*swdown_c(:,:,m0(doy))    + wtm1(doy)*swdown_c(:,:,m1(doy))

    ! derive diffuse downward surface radiation
    where (atm%cld.gt.0.1_wp)
      atm%swdown_d = (atm%swdown-(1._wp-atm%cld)*atm%swdown_c) / atm%cld
    elsewhere
      atm%swdown_d = atm%swdown
    endwhere

    where (atm%qair.lt.1e-6_wp)
      atm%qair = 1.e-6_wp
    endwhere
    where (atm%wind .le. 0._wp)
      atm%wind = 1.e-3_wp
    endwhere

  endif


  return

  end subroutine fake_atm_update


end module fake_atm_mod
