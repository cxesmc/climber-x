!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : fake_atm_hires_mod 
!
!  Purpose : prescribed high-res atmospheric variables
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2024 Potsdam Institute for Climate Impact Research,
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
module fake_atm_hires_mod

  use precision, only : wp, dp
  use timer, only: doy, nday_year, time_soy_bnd, monthly2daily
  use smb_params, only : ifake_atm_hires, fake_atm_hires_const_file, fake_atm_hires_var_file, prc_forcing, wind_forcing
  use coord, only : grid_init, grid_class
  use coord, only : map_scrip_init, map_scrip_class
  use ncio

  implicit none

  type fake_atm_hires_type
    type(grid_class) :: grid
    type(map_scrip_class) :: maps_atm_to_smb
    real(wp), dimension(:,:), allocatable :: tair, tstd, qair, rain, snow, prc, pressure, swdown, lwdown, usur, vsur, wind, cld, cod, alb
  end type
    
  real(wp), dimension(:,:,:), allocatable :: tair, tstd, qair, rain, snow, prc, pressure, swdown, lwdown, usur, vsur, wind, cld, cod, alb
  real(wp), dimension(:,:,:), allocatable :: tair_0, tstd_0, qair_0, rain_0, snow_0, prc_0, pressure_0, swdown_0, lwdown_0, usur_0, vsur_0, wind_0, cld_0, cod_0, alb_0
  real(wp), dimension(:,:,:), allocatable :: tair_1, tstd_1, qair_1, rain_1, snow_1, prc_1, pressure_1, swdown_1, lwdown_1, usur_1, vsur_1, wind_1, cld_1, cod_1, alb_1
  real(wp), dimension(:), allocatable :: time_atm

  integer :: ni, nj
  integer :: i0, i1, i1_old
  integer, dimension(nday_year) :: m0, m1
  real(wp), dimension(nday_year) :: wtm0, wtm1
 
  private
  public :: fake_atm_hires_init, fake_atm_hires_update, fake_atm_hires_type

contains

  subroutine fake_atm_hires_init(smb_grid,time_now,atm)

  implicit none

  type(grid_class), intent(in) :: smb_grid
  real(wp), intent(in) :: time_now
  type(fake_atm_hires_type), intent(inout) :: atm

  integer :: imin, ntime
  integer :: ppos, spos
  real(wp) :: w0, w1
  real(wp) :: dlon, dlat
  real(wp), dimension(:), allocatable :: lon, lat
  character(len=256) :: fnm


  if (ifake_atm_hires.eq.0) then
    fnm = fake_atm_hires_const_file
  else if (ifake_atm_hires.eq.1) then
    fnm = fake_atm_hires_var_file
  endif

  ! read dimensions
  ni = nc_size(trim(fnm),"lon")
  nj = nc_size(trim(fnm),"lat")
  ! read lat lon
  allocate(lon(ni))
  allocate(lat(nj))
  call nc_read(trim(fnm),"lon",lon)
  call nc_read(trim(fnm),"lat",lat)
  dlon = lon(2)-lon(1)
  dlat = lat(2)-lat(1)

  ! create atm grid object
  spos = scan(trim(fnm),"/", BACK= .true.)+1
  ppos = scan(trim(fnm),".", BACK= .true.)-1
  call grid_init(atm%grid,name=trim(fnm(spos:ppos)),mtype="latlon",units="degrees", &
    x0=real(lon(1),dp),dx=real(dlon,dp),nx=ni,y0=real(lat(1),dp),dy=real(dlat,dp),ny=nj)

  ! initialize map
  call map_scrip_init(atm%maps_atm_to_smb,atm%grid,smb_grid,method="con",fldr="maps",load=.TRUE.,clean=.FALSE.)

  ! allocate atm type variables
  allocate(atm%tair(ni,nj))
  allocate(atm%tstd(ni,nj))
  allocate(atm%qair(ni,nj))
  allocate(atm%prc(ni,nj))
  allocate(atm%rain(ni,nj))
  allocate(atm%snow(ni,nj))
  allocate(atm%pressure(ni,nj))
  allocate(atm%swdown(ni,nj))
  allocate(atm%lwdown(ni,nj))
  allocate(atm%usur(ni,nj))
  allocate(atm%vsur(ni,nj))
  allocate(atm%wind(ni,nj))
  allocate(atm%cld(ni,nj))
  allocate(atm%cod(ni,nj))
  allocate(atm%alb(ni,nj))

  ! allocate local variables
  allocate(tair(ni,nj,0:13))
  allocate(tstd(ni,nj,0:13))
  allocate(qair(ni,nj,0:13))
  allocate(rain(ni,nj,0:13))
  allocate(snow(ni,nj,0:13))
  allocate(prc(ni,nj,0:13))
  allocate(pressure(ni,nj,0:13))
  allocate(swdown(ni,nj,0:13))
  allocate(lwdown(ni,nj,0:13))
  allocate(usur(ni,nj,0:13))
  allocate(vsur(ni,nj,0:13))
  allocate(wind(ni,nj,0:13))
  allocate(cld(ni,nj,0:13))
  allocate(cod(ni,nj,0:13))
  allocate(alb(ni,nj,0:13))
  if (ifake_atm_hires.eq.1) then
    allocate(tair_0(ni,nj,12))
    allocate(tstd_0(ni,nj,12))
    allocate(qair_0(ni,nj,12))
    allocate(rain_0(ni,nj,12))
    allocate(snow_0(ni,nj,12))
    allocate(prc_0(ni,nj,12))
    allocate(pressure_0(ni,nj,12))
    allocate(swdown_0(ni,nj,12))
    allocate(lwdown_0(ni,nj,12))
    allocate(wind_0(ni,nj,12))
    allocate(cld_0(ni,nj,12))
    allocate(cod_0(ni,nj,12))
    allocate(alb_0(ni,nj,12))
    allocate(usur_0(ni,nj,12))
    allocate(vsur_0(ni,nj,12))
    allocate(tair_1(ni,nj,12))
    allocate(tstd_1(ni,nj,12))
    allocate(qair_1(ni,nj,12))
    allocate(rain_1(ni,nj,12))
    allocate(snow_1(ni,nj,12))
    allocate(prc_1(ni,nj,12))
    allocate(pressure_1(ni,nj,12))
    allocate(swdown_1(ni,nj,12))
    allocate(lwdown_1(ni,nj,12))
    allocate(usur_1(ni,nj,12))
    allocate(vsur_1(ni,nj,12))
    allocate(wind_1(ni,nj,12))
    allocate(cld_1(ni,nj,12))
    allocate(cod_1(ni,nj,12))
    allocate(alb_1(ni,nj,12))
  endif

  ! get weights for interpolation from monthly to daily
  call monthly2daily(m0,m1,wtm0,wtm1)

  ! read constant climate forcing
  if (ifake_atm_hires.eq.0) then

    call nc_read(fake_atm_hires_const_file,"t_air",tair(:,:,1:12) )
    !call nc_read(fake_atm_hires_const_file,"tstd",tstd(:,:,1:12) )
    call nc_read(fake_atm_hires_const_file,"tmaxstd",tstd(:,:,1:12) )
    call nc_read(fake_atm_hires_const_file,"qair",qair(:,:,1:12) )
    if (prc_forcing.eq.0) then
      ! read rain and snow separately
      call nc_read(fake_atm_hires_const_file,"rain",rain(:,:,1:12) )
      call nc_read(fake_atm_hires_const_file,"snow",snow(:,:,1:12) )
    else if (prc_forcing.eq.1) then
      call nc_read(fake_atm_hires_const_file,"prc",prc(:,:,1:12) )
    endif
    call nc_read(fake_atm_hires_const_file,"pressure",pressure(:,:,1:12) )
    call nc_read(fake_atm_hires_const_file,"swdown",swdown(:,:,1:12) )
    call nc_read(fake_atm_hires_const_file,"lwdown",lwdown(:,:,1:12) )
    if (wind_forcing.eq.0) then
      ! read u and v wind separately
      call nc_read(fake_atm_hires_const_file,"usur",usur(:,:,1:12) )
      call nc_read(fake_atm_hires_const_file,"vsur",vsur(:,:,1:12) )
      ! derive wind speed, not correct but the best that can be done
      wind(:,:,1:12) = sqrt(usur(:,:,1:12)**2+vsur(:,:,1:12)**2)
    else
      call nc_read(fake_atm_hires_const_file,"wind",wind(:,:,1:12) )
    endif
    call nc_read(fake_atm_hires_const_file,"cld",cld(:,:,1:12) )
    call nc_read(fake_atm_hires_const_file,"cod",cod(:,:,1:12) )
    call nc_read(fake_atm_hires_const_file,"alb",alb(:,:,1:12) )

  ! initilaize variable climate forcing
  else if (ifake_atm_hires.eq.1) then

    ntime = nc_size(trim(fake_atm_hires_var_file),"time")
    allocate( time_atm(ntime) )
    call nc_read(trim(fake_atm_hires_var_file),"time",time_atm)
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
    call nc_read(fake_atm_hires_var_file,"tair",tair_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"tstd",tstd_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"qair",qair_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"rain",rain_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"snow",snow_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    prc_0(:,:,1:12) = rain_0(:,:,1:12) + snow_0(:,:,1:12)
    call nc_read(fake_atm_hires_var_file,"pressure",pressure_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"swdown",swdown_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"lwdown",lwdown_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"usur",usur_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"vsur",vsur_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"wind",wind_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"cld",cld_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"cod",cod_0,start=[1,1,1,i0],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"alb",alb_0,start=[1,1,1,i0],count=[ni,nj,12,1] )

    call nc_read(fake_atm_hires_var_file,"tair",tair_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"tstd",tstd_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"qair",qair_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"rain",rain_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"snow",snow_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    prc_1(:,:,1:12) = rain_1(:,:,1:12) + snow_1(:,:,1:12)
    call nc_read(fake_atm_hires_var_file,"pressure",pressure_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"swdown",swdown_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"lwdown",lwdown_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"usur",usur_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"vsur",vsur_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"wind",wind_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"cld",cld_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"cod",cod_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
    call nc_read(fake_atm_hires_var_file,"alb",alb_1,start=[1,1,1,i1],count=[ni,nj,12,1] )

    tair(:,:,1:12) = w0*tair_0 + w1*tair_1
    tstd(:,:,1:12) = w0*tstd_0 + w1*tstd_1
    qair(:,:,1:12) = w0*qair_0 + w1*qair_1
    rain(:,:,1:12) = w0*rain_0 + w1*rain_1
    snow(:,:,1:12) = w0*snow_0 + w1*snow_1
    prc(:,:,1:12) = w0*prc_0 + w1*prc_1
    usur(:,:,1:12) = w0*usur_0 + w1*usur_1
    vsur(:,:,1:12) = w0*vsur_0 + w1*vsur_1
    wind(:,:,1:12) = w0*wind_0 + w1*wind_1
    cld(:,:,1:12) = w0*cld_0 + w1*cld_1
    cod(:,:,1:12) = w0*cod_0 + w1*cod_1
    alb(:,:,1:12) = w0*alb_0 + w1*alb_1
    pressure(:,:,1:12) = w0*pressure_0 + w1*pressure_1
    lwdown(:,:,1:12) = w0*lwdown_0 + w1*lwdown_1
    swdown(:,:,1:12) = w0*swdown_0 + w1*swdown_1

  endif

  tair(:,:,0)  = tair(:,:,12)
  tstd(:,:,0)  = tstd(:,:,12)
  qair(:,:,0) = qair(:,:,12)
  prc(:,:,0)  = prc(:,:,12)
  rain(:,:,0)  = rain(:,:,12)
  snow(:,:,0) = snow(:,:,12)
  pressure(:,:,0)  = pressure(:,:,12)
  swdown(:,:,0) = swdown(:,:,12)
  lwdown(:,:,0)  = lwdown(:,:,12)
  usur(:,:,0) = usur(:,:,12)
  vsur(:,:,0) = vsur(:,:,12)
  wind(:,:,0) = wind(:,:,12)
  cld(:,:,0) = cld(:,:,12)
  cod(:,:,0) = cod(:,:,12)
  alb(:,:,0) = alb(:,:,12)

  tair(:,:,13)  = tair(:,:,1)
  tstd(:,:,13)  = tstd(:,:,1)
  qair(:,:,13) = qair(:,:,1)
  prc(:,:,13)  = prc(:,:,1)
  rain(:,:,13)  = rain(:,:,1)
  snow(:,:,13) = snow(:,:,1)
  pressure(:,:,13)  = pressure(:,:,1)
  swdown(:,:,13) = swdown(:,:,1)
  lwdown(:,:,13)  = lwdown(:,:,1)
  usur(:,:,13) = usur(:,:,1)
  vsur(:,:,13) = vsur(:,:,1)
  wind(:,:,13) = wind(:,:,1)
  cld(:,:,13) = cld(:,:,1)
  cod(:,:,13) = cod(:,:,1)
  alb(:,:,13) = alb(:,:,1)

  where (wind .le. 0._wp)
    wind = 1.e-3_wp
  endwhere

  ! initialize atm
  atm%tair  = tair(:,:,1)   
  atm%tstd  = tstd(:,:,1)   
  atm%qair  = qair(:,:,1)   
  if (prc_forcing.eq.0) then
    atm%rain     = rain(:,:,1) 
    atm%snow     = snow(:,:,1) 
  else if (prc_forcing.eq.1) then
    atm%prc      = prc(:,:,1)  
  endif
  atm%pressure = pressure(:,:,1) 
  atm%lwdown  = lwdown(:,:,1) 
  atm%swdown  = swdown(:,:,1)
  atm%usur     = usur(:,:,1)    
  atm%vsur     = vsur(:,:,1)   
  atm%wind     = wind(:,:,1) 
  atm%cld     = cld(:,:,1) 
  atm%cod     = cod(:,:,1) 
  atm%alb     = alb(:,:,1) 


  return

  end subroutine fake_atm_hires_init


  subroutine fake_atm_hires_update(time_now,atm)

  implicit none

  real(wp), intent(in) :: time_now
  type(fake_atm_hires_type)  :: atm

  integer :: imin
  real(wp) :: w0, w1


  ! constant climate forcing
  if (ifake_atm_hires.eq.0) then

    if (prc_forcing.eq.0) then
      atm%rain     = wtm0(doy)*rain(:,:,m0(doy))     + wtm1(doy)*rain(:,:,m1(doy))
      atm%snow     = wtm0(doy)*snow(:,:,m0(doy))     + wtm1(doy)*snow(:,:,m1(doy))
    else if (prc_forcing.eq.1) then
      atm%prc      = wtm0(doy)*prc(:,:,m0(doy))     + wtm1(doy)*prc(:,:,m1(doy))
    endif
    atm%pressure = wtm0(doy)*pressure(:,:,m0(doy)) + wtm1(doy)*pressure(:,:,m1(doy))
    atm%lwdown  = wtm0(doy)*lwdown(:,:,m0(doy))   + wtm1(doy)*lwdown(:,:,m1(doy))
    atm%swdown  = wtm0(doy)*swdown(:,:,m0(doy))   + wtm1(doy)*swdown(:,:,m1(doy))
    atm%wind     = wtm0(doy)*wind(:,:,m0(doy))     + wtm1(doy)*wind(:,:,m1(doy))
    atm%cld     = wtm0(doy)*cld(:,:,m0(doy))     + wtm1(doy)*cld(:,:,m1(doy))
    atm%cod     = wtm0(doy)*cod(:,:,m0(doy))     + wtm1(doy)*cod(:,:,m1(doy))
    atm%alb     = wtm0(doy)*alb(:,:,m0(doy))     + wtm1(doy)*alb(:,:,m1(doy))
    atm%tair  = wtm0(doy)*tair(:,:,m0(doy))    + wtm1(doy)*tair(:,:,m1(doy))
    atm%tstd  = wtm0(doy)*tstd(:,:,m0(doy))    + wtm1(doy)*tstd(:,:,m1(doy))
    atm%qair  = wtm0(doy)*qair(:,:,m0(doy))    + wtm1(doy)*qair(:,:,m1(doy))

    where (atm%qair.lt.1e-4_wp)
      atm%qair = 1.e-4_wp
    endwhere
    where (atm%swdown.lt.0._wp)
      atm%swdown = 0._wp
    endwhere
    where (atm%rain.lt.0._wp)
      atm%rain = 0._wp
    endwhere
    where (atm%snow.lt.0._wp)
      atm%snow = 0._wp
    endwhere

    ! variable climate forcing
  else if (ifake_atm_hires.eq.1) then
    if (time_soy_bnd) then
      ! look for time slice closest to current time
      i1_old = i1
      if (time_now.lt.time_atm(lbound(time_atm,1))) then
        i0 = 1
        i1 = 1
        w0 = 1._wp
        w1 = 0._wp
        !stop 'ERROR, fake_atm_hires not defined for current time'
      else if (time_now.gt.time_atm(ubound(time_atm,1))) then
        stop 'ERROR, fake_atm_hires mask not defined for current time'
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
        tstd_0 = tstd_1
        qair_0 = qair_1
        rain_0 = rain_1
        prc_0 = prc_1
        snow_0 = snow_1
        wind_0 = wind_1
        cld_0 = cld_1
        cod_0 = cod_1
        alb_0 = alb_1
        usur_0 = usur_1
        vsur_0 = vsur_1
        pressure_0 = pressure_1
        lwdown_0 = lwdown_1
        swdown_0 = swdown_1
        call nc_read(trim(fake_atm_hires_var_file),"tair",tair_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_hires_var_file),"tstd",tstd_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_hires_var_file),"qair",qair_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_hires_var_file),"rain",rain_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_hires_var_file),"snow",snow_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        prc_1(:,:,1:12) = rain_1(:,:,1:12) + snow_1(:,:,1:12)
        call nc_read(trim(fake_atm_hires_var_file),"usur",usur_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_hires_var_file),"vsur",vsur_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_hires_var_file),"wind",wind_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_hires_var_file),"cld",cld_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_hires_var_file),"cod",cod_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_hires_var_file),"alb",alb_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_hires_var_file),"pressure",pressure_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_hires_var_file),"lwdown",lwdown_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
        call nc_read(trim(fake_atm_hires_var_file),"swdown",swdown_1,start=[1,1,1,i1],count=[ni,nj,12,1] )
      endif
      tair(:,:,1:12) = w0*tair_0 + w1*tair_1
      tstd(:,:,1:12) = w0*tstd_0 + w1*tstd_1
      qair(:,:,1:12) = w0*qair_0 + w1*qair_1
      rain(:,:,1:12) = w0*rain_0 + w1*rain_1
      prc(:,:,1:12) = w0*prc_0 + w1*prc_1
      snow(:,:,1:12) = w0*snow_0 + w1*snow_1
      usur(:,:,1:12)    = w0*usur_0 + w1*usur_1
      vsur(:,:,1:12)    = w0*vsur_0 + w1*vsur_1
      wind(:,:,1:12) = w0*wind_0 + w1*wind_1
      cld(:,:,1:12) = w0*cld_0 + w1*cld_1
      cod(:,:,1:12) = w0*cod_0 + w1*cod_1
      alb(:,:,1:12) = w0*alb_0 + w1*alb_1
      pressure(:,:,1:12) = w0*pressure_0 + w1*pressure_1
      lwdown(:,:,1:12) = w0*lwdown_0 + w1*lwdown_1
      swdown(:,:,1:12) = w0*swdown_0 + w1*swdown_1

      ! extended monthly grid
      tair(:,:,0)  = tair(:,:,12)
      tstd(:,:,0)  = tstd(:,:,12)
      qair(:,:,0) = qair(:,:,12)
      prc(:,:,0)  = prc(:,:,12)
      rain(:,:,0)  = rain(:,:,12)
      snow(:,:,0) = snow(:,:,12)
      pressure(:,:,0)  = pressure(:,:,12)
      swdown(:,:,0) = swdown(:,:,12)
      lwdown(:,:,0)  = lwdown(:,:,12)
      usur(:,:,0) = usur(:,:,12)
      vsur(:,:,0) = vsur(:,:,12)
      wind(:,:,0) = wind(:,:,12)
      cld(:,:,0) = cld(:,:,12)
      cod(:,:,0) = cod(:,:,12)
      alb(:,:,0) = alb(:,:,12)

      tair(:,:,13)  = tair(:,:,1)
      tstd(:,:,13)  = tstd(:,:,1)
      qair(:,:,13) = qair(:,:,1)
      prc(:,:,13)  = prc(:,:,1)
      rain(:,:,13)  = rain(:,:,1)
      snow(:,:,13) = snow(:,:,1)
      pressure(:,:,13)  = pressure(:,:,1)
      swdown(:,:,13) = swdown(:,:,1)
      lwdown(:,:,13)  = lwdown(:,:,1)
      usur(:,:,13) = usur(:,:,1)
      vsur(:,:,13) = vsur(:,:,1)
      wind(:,:,13) = wind(:,:,1)
      cld(:,:,13) = cld(:,:,1)
      cod(:,:,13) = cod(:,:,1)
      alb(:,:,13) = alb(:,:,1)

    endif

    atm%tair  = wtm0(doy)*tair(:,:,m0(doy))    + wtm1(doy)*tair(:,:,m1(doy))
    atm%tstd  = wtm0(doy)*tstd(:,:,m0(doy))    + wtm1(doy)*tstd(:,:,m1(doy))
    atm%qair  = wtm0(doy)*qair(:,:,m0(doy))    + wtm1(doy)*qair(:,:,m1(doy))
    atm%rain  = wtm0(doy)*rain(:,:,m0(doy))    + wtm1(doy)*rain(:,:,m1(doy))
    atm%snow  = wtm0(doy)*snow(:,:,m0(doy))    + wtm1(doy)*snow(:,:,m1(doy))
    atm%prc   = wtm0(doy)*prc(:,:,m0(doy))     + wtm1(doy)*prc(:,:,m1(doy))
    atm%usur  = wtm0(doy)*usur(:,:,m0(doy))    + wtm1(doy)*usur(:,:,m1(doy))
    atm%vsur  = wtm0(doy)*vsur(:,:,m0(doy))    + wtm1(doy)*vsur(:,:,m1(doy))
    atm%wind  = wtm0(doy)*wind(:,:,m0(doy))    + wtm1(doy)*wind(:,:,m1(doy))
    atm%cld   = wtm0(doy)*cld(:,:,m0(doy))     + wtm1(doy)*cld(:,:,m1(doy))
    atm%cod   = wtm0(doy)*cod(:,:,m0(doy))     + wtm1(doy)*cod(:,:,m1(doy))
    atm%alb   = wtm0(doy)*alb(:,:,m0(doy))     + wtm1(doy)*alb(:,:,m1(doy))
    atm%pressure  = wtm0(doy)*pressure(:,:,m0(doy))    + wtm1(doy)*pressure(:,:,m1(doy))
    atm%lwdown  = wtm0(doy)*lwdown(:,:,m0(doy))    + wtm1(doy)*lwdown(:,:,m1(doy))
    atm%swdown  = wtm0(doy)*swdown(:,:,m0(doy))    + wtm1(doy)*swdown(:,:,m1(doy))

    where (atm%qair.lt.1e-4_wp)
      atm%qair = 1.e-4_wp
    endwhere

  endif


  return

  end subroutine fake_atm_hires_update


end module fake_atm_hires_mod
