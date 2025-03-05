!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s i c _ o u t
!
!  Purpose : sea ice model diagnostics and output
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
module sic_out

  use precision, only : wp
  use constants, only : pi, g, cap_w, Lf
  use dim_name, only: dim_lon, dim_lat, dim_depth, dim_time, dim_month, dim_day
  use timer, only : n_accel, nyears, year, year_clim, year_now, sec_day, mon, nmon_year, doy, nday_year
  use timer, only : nstep_year_sic, nstep_mon_sic, nyout_sic, ny_out_ts, y_out_ts_clim, time_out_ts_clim
  use timer, only : time_soy_sic, time_eoy_sic, time_out_sic
  use control, only : out_dir
  use climber_grid, only : basin_mask, basin_mask2, i_atlantic, i_pacific, i_indian, i_southern
  use climber_grid, only : lon, lat, lonu, latv
  use sic_grid, only : maxi, maxj, phi0, dphi, sv, dxv, area
  use sic_params, only : dt, l_daily_output, l_diag_dyn, rho_sic
  use sic_def, only : sic_class
  use ncio

  implicit none

  private
  public :: sic_diag, sic_diag_init

  integer :: nout

  integer :: j_fram, i_fram(2)
  integer :: j_denmark, i_denmark(2)
  integer, parameter :: nlatv_buoy = 6
  real(wp), dimension(nlatv_buoy) :: latv_buoy = (/40._wp,45._wp,50._wp,55._wp,60._wp,65._wp/)
  integer, parameter :: ilatv_buoy_sel = 4

  type ts_out
     integer :: ncells
     real(wp) :: area
     real(wp) :: a_nh, a_sh, e_nh, e_sh, v_nh, v_sh
     real(wp) :: a_nh_min, a_nh_max, a_sh_min, a_sh_max
     real(wp) :: e_nh_min, e_nh_max, e_sh_min, e_sh_max
     real(wp) :: v_nh_min, v_nh_max, v_sh_min, v_sh_max
     real(wp) :: fram_exp, denmark_exp, buoy_sic_NA(nlatv_buoy)
  end type

  type s_out
     real(wp), dimension(:,:), allocatable :: hsic, hsnow, focn, fsic, tsic, tocn
     real(wp), dimension(:,:), allocatable :: dh_sic_dt_therm, dh_sic_dt_dyn, dh_sic_dt_therm_ocn, dh_sic_dt_therm_sic
     real(wp), dimension(:,:), allocatable :: lh, sh, lw, fx, fw, fw_brines, e, p_e
     real(wp), dimension(:,:), allocatable :: lh_sic, sh_sic, lw_sic, fx_sic, fw_sic
     real(wp), dimension(:,:), allocatable :: lh_ocn, sh_ocn, lw_ocn, fx_ocn, fw_ocn
     real(wp), dimension(:,:), allocatable :: flx_melt_top, flx_melt_bot
     real(wp), dimension(:,:), allocatable :: e_sic, e_ocn
     real(wp), dimension(:,:), allocatable :: cde_ocn
     real(wp), dimension(:,:), allocatable :: cde_sic
     real(wp), dimension(:,:), allocatable :: cdh_ocn
     real(wp), dimension(:,:), allocatable :: cdh_sic
     real(wp), dimension(:,:), allocatable :: snow_grain
     real(wp), dimension(:,:), allocatable :: dust_con
     real(wp), dimension(:,:), allocatable :: alb_sic, alb_ocn
     real(wp), dimension(:,:), allocatable :: rain, snow
     real(wp), dimension(:,:), allocatable :: usic, vsic
     real(wp), dimension(:,:), allocatable :: tauxa, tauya
     real(wp), dimension(:,:), allocatable :: tauxo, tauyo
     real(wp), dimension(:,:), allocatable :: str_d ! The divergence stress tensor component [Pa m].
     real(wp), dimension(:,:), allocatable :: str_t ! The tension stress tensor component [Pa m].
     real(wp), dimension(:,:), allocatable :: str_s ! The shearing stress tensor component [Pa m].
     real(wp), dimension(:,:), allocatable :: fxic   ! Zonal force due to internal stresses [Pa].
     real(wp), dimension(:,:), allocatable :: fxic_d ! Zonal force due to divergence internal stress [Pa].
     real(wp), dimension(:,:), allocatable :: fxic_t ! Zonal force due to tension internal stress [Pa].
     real(wp), dimension(:,:), allocatable :: fxic_s ! Zonal force due to shearing internal stress [Pa].
     real(wp), dimension(:,:), allocatable :: Cor_u  ! Zonal Coriolis acceleration [m s-2].
     real(wp), dimension(:,:), allocatable :: PFu    ! Zonal hydrostatic pressure driven acceleration [m s-2].
     real(wp), dimension(:,:), allocatable :: fyic   ! Meridional force due to internal stresses [Pa].
     real(wp), dimension(:,:), allocatable :: fyic_d ! Meridional force due to divergence internal stress [Pa].
     real(wp), dimension(:,:), allocatable :: fyic_t ! Meridional force due to tension internal stress [Pa].
     real(wp), dimension(:,:), allocatable :: fyic_s ! Meridional force due to shearing internal stress [Pa].
     real(wp), dimension(:,:), allocatable :: Cor_v  ! Meridional Coriolis acceleration [m s-2].
     real(wp), dimension(:,:), allocatable :: PFv    ! Meridional hydrostatic pressure driven acceleration [m s-2].
     real(wp), dimension(:,:), allocatable :: fwt, fwp, fwa
  end type

  type(ts_out), allocatable :: mon_ts(:,:), ann_ts(:)
  type(s_out) :: day_s(nday_year), mon_s(nmon_year), ann_s


contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s i c _ d i a g _ i n i t
  ! Purpose  :  Initialize netcdf output for sea ice
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_diag_init

    implicit none

    integer :: i, j, k
    real(wp) :: lat_fram, lon_fram_1, lon_fram_2, dist_fram
    real(wp) :: lat_denmark, lon_denmark_1, lon_denmark_2, dist_denmark


    nout = 0

    ! initialise output
    call ts_nc(trim(out_dir)//"/sic_ts.nc")
    call sic_nc(trim(out_dir)//"/sic.nc")
    if (l_daily_output) call sic_daily_nc(trim(out_dir)//"/sic_daily.nc")

    ! locate Fram Strait
    lon_fram_1 = -30.
    i_fram(1) = 0
    dist_fram = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_fram_1/180.0)),2.0*pi).lt.dist_fram) then
          i_fram(1) = i
          dist_fram = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_fram_1/180.0)),2.0*pi)
       endif
    enddo
    lon_fram_2 = 15.
    i_fram(2) = 0
    dist_fram = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_fram_2/180.0)),2.0*pi).lt.dist_fram) then
          i_fram(2) = i
          dist_fram = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_fram_2/180.0)),2.0*pi)
       endif
    enddo
    lat_fram = 80.
    j_fram = 0
    dist_fram = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_fram/180.0)).lt.dist_fram) then
          j_fram = j
          dist_fram = abs(sv(j)-sin(pi*lat_fram/180.0))
       endif
    enddo

    ! locate denmark Strait
    lon_denmark_1 = -47.5
    i_denmark(1) = 0
    dist_denmark = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_denmark_1/180.0)),2.0*pi).lt.dist_denmark) then
          i_denmark(1) = i
          dist_denmark = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_denmark_1/180.0)),2.0*pi)
       endif
    enddo
    lon_denmark_2 = -17.5
    i_denmark(2) = 0
    dist_denmark = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_denmark_2/180.0)),2.0*pi).lt.dist_denmark) then
          i_denmark(2) = i
          dist_denmark = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_denmark_2/180.0)),2.0*pi)
       endif
    enddo
    lat_denmark = 67.5
    j_denmark = 0
    dist_denmark = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_denmark/180.0)).lt.dist_denmark) then
          j_denmark = j
          dist_denmark = abs(sv(j)-sin(pi*lat_denmark/180.0))
       endif
    enddo

    ! allocate
    allocate(ann_ts(ny_out_ts))
    allocate(mon_ts(ny_out_ts,nmon_year))

    allocate(ann_s%hsic(maxi,maxj))
    allocate(ann_s%hsnow(maxi,maxj))
    allocate(ann_s%focn(maxi,maxj))
    allocate(ann_s%fsic(maxi,maxj))
    allocate(ann_s%tsic(maxi,maxj))
    allocate(ann_s%dh_sic_dt_therm(maxi,maxj))
    allocate(ann_s%dh_sic_dt_therm_sic(maxi,maxj))
    allocate(ann_s%dh_sic_dt_therm_ocn(maxi,maxj))
    allocate(ann_s%dh_sic_dt_dyn(maxi,maxj))
    allocate(ann_s%flx_melt_top(maxi,maxj))
    allocate(ann_s%flx_melt_bot(maxi,maxj))
    allocate(ann_s%tocn(maxi,maxj))
    allocate(ann_s%lh(maxi,maxj))
    allocate(ann_s%sh(maxi,maxj))
    allocate(ann_s%lw(maxi,maxj))
    allocate(ann_s%fx(maxi,maxj))
    allocate(ann_s%fw(maxi,maxj))
    allocate(ann_s%fw_brines(maxi,maxj))
    allocate(ann_s%e(maxi,maxj))
    allocate(ann_s%p_e(maxi,maxj))
    allocate(ann_s%lh_sic(maxi,maxj))
    allocate(ann_s%sh_sic(maxi,maxj))
    allocate(ann_s%lw_sic(maxi,maxj))
    allocate(ann_s%fx_sic(maxi,maxj))
    allocate(ann_s%fw_sic(maxi,maxj))
    allocate(ann_s%e_sic(maxi,maxj))
    allocate(ann_s%lh_ocn(maxi,maxj))
    allocate(ann_s%sh_ocn(maxi,maxj))
    allocate(ann_s%lw_ocn(maxi,maxj))
    allocate(ann_s%fx_ocn(maxi,maxj))
    allocate(ann_s%fw_ocn(maxi,maxj))
    allocate(ann_s%e_ocn(maxi,maxj))
    allocate(ann_s%cde_ocn(maxi,maxj))
    allocate(ann_s%cde_sic(maxi,maxj))
    allocate(ann_s%cdh_ocn(maxi,maxj))
    allocate(ann_s%cdh_sic(maxi,maxj))
    allocate(ann_s%snow_grain(maxi,maxj))
    allocate(ann_s%dust_con(maxi,maxj))
    allocate(ann_s%alb_ocn(maxi,maxj))
    allocate(ann_s%alb_sic(maxi,maxj))
    allocate(ann_s%rain(maxi,maxj))
    allocate(ann_s%snow(maxi,maxj))
    allocate(ann_s%usic(maxi,maxj))
    allocate(ann_s%vsic(maxi,maxj))
    allocate(ann_s%tauxa(maxi,maxj))
    allocate(ann_s%tauya(maxi,maxj))
    allocate(ann_s%tauxo(maxi,maxj))
    allocate(ann_s%tauyo(maxi,maxj))
    allocate(ann_s%str_d(maxi,maxj))
    allocate(ann_s%str_t(maxi,maxj))
    allocate(ann_s%str_s(maxi,maxj))
    allocate(ann_s%fxic  (maxi,maxj)) 
    allocate(ann_s%fxic_d(maxi,maxj)) 
    allocate(ann_s%fxic_t(maxi,maxj)) 
    allocate(ann_s%fxic_s(maxi,maxj)) 
    allocate(ann_s%Cor_u (maxi,maxj)) 
    allocate(ann_s%PFu   (maxi,maxj)) 
    allocate(ann_s%fyic  (maxi,maxj)) 
    allocate(ann_s%fyic_d(maxi,maxj)) 
    allocate(ann_s%fyic_t(maxi,maxj)) 
    allocate(ann_s%fyic_s(maxi,maxj)) 
    allocate(ann_s%Cor_v (maxi,maxj)) 
    allocate(ann_s%PFv   (maxi,maxj)) 
    allocate(ann_s%fwt(3,maxj))
    allocate(ann_s%fwp(3,maxj))
    allocate(ann_s%fwa(3,maxj))

    do k=1,nmon_year
     allocate(mon_s(k)%hsic(maxi,maxj))
     allocate(mon_s(k)%hsnow(maxi,maxj))
     allocate(mon_s(k)%fsic(maxi,maxj))
     allocate(mon_s(k)%tsic(maxi,maxj))
     allocate(mon_s(k)%tocn(maxi,maxj))
     allocate(mon_s(k)%dh_sic_dt_therm(maxi,maxj))
     allocate(mon_s(k)%dh_sic_dt_therm_sic(maxi,maxj))
     allocate(mon_s(k)%dh_sic_dt_therm_ocn(maxi,maxj))
     allocate(mon_s(k)%dh_sic_dt_dyn(maxi,maxj))
     allocate(mon_s(k)%flx_melt_top(maxi,maxj))
     allocate(mon_s(k)%flx_melt_bot(maxi,maxj))
     allocate(mon_s(k)%lh(maxi,maxj))
     allocate(mon_s(k)%sh(maxi,maxj))
     allocate(mon_s(k)%lw(maxi,maxj))
     allocate(mon_s(k)%fx(maxi,maxj))
     allocate(mon_s(k)%fw(maxi,maxj))
     allocate(mon_s(k)%fw_brines(maxi,maxj))
     allocate(mon_s(k)%e(maxi,maxj))
     allocate(mon_s(k)%p_e(maxi,maxj))
     allocate(mon_s(k)%lh_sic(maxi,maxj))
     allocate(mon_s(k)%sh_sic(maxi,maxj))
     allocate(mon_s(k)%lw_sic(maxi,maxj))
     allocate(mon_s(k)%fx_sic(maxi,maxj))
     allocate(mon_s(k)%fw_sic(maxi,maxj))
     allocate(mon_s(k)%e_sic(maxi,maxj))
     allocate(mon_s(k)%lh_ocn(maxi,maxj))
     allocate(mon_s(k)%sh_ocn(maxi,maxj))
     allocate(mon_s(k)%lw_ocn(maxi,maxj))
     allocate(mon_s(k)%fx_ocn(maxi,maxj))
     allocate(mon_s(k)%fw_ocn(maxi,maxj))
     allocate(mon_s(k)%e_ocn(maxi,maxj))
     allocate(mon_s(k)%cde_ocn(maxi,maxj))
     allocate(mon_s(k)%cde_sic(maxi,maxj))
     allocate(mon_s(k)%cdh_ocn(maxi,maxj))
     allocate(mon_s(k)%cdh_sic(maxi,maxj))
     allocate(mon_s(k)%snow_grain(maxi,maxj))
     allocate(mon_s(k)%dust_con(maxi,maxj))
     allocate(mon_s(k)%alb_ocn(maxi,maxj))
     allocate(mon_s(k)%alb_sic(maxi,maxj))
     allocate(mon_s(k)%rain(maxi,maxj))
     allocate(mon_s(k)%snow(maxi,maxj))
     allocate(mon_s(k)%usic(maxi,maxj))
     allocate(mon_s(k)%vsic(maxi,maxj))
     allocate(mon_s(k)%tauxa(maxi,maxj))
     allocate(mon_s(k)%tauya(maxi,maxj))
     allocate(mon_s(k)%tauxo(maxi,maxj))
     allocate(mon_s(k)%tauyo(maxi,maxj))
     allocate(mon_s(k)%str_d(maxi,maxj))
     allocate(mon_s(k)%str_t(maxi,maxj))
     allocate(mon_s(k)%str_s(maxi,maxj))
     allocate(mon_s(k)%fxic  (maxi,maxj)) 
     allocate(mon_s(k)%fxic_d(maxi,maxj)) 
     allocate(mon_s(k)%fxic_t(maxi,maxj)) 
     allocate(mon_s(k)%fxic_s(maxi,maxj)) 
     allocate(mon_s(k)%Cor_u (maxi,maxj)) 
     allocate(mon_s(k)%PFu   (maxi,maxj)) 
     allocate(mon_s(k)%fyic  (maxi,maxj)) 
     allocate(mon_s(k)%fyic_d(maxi,maxj)) 
     allocate(mon_s(k)%fyic_t(maxi,maxj)) 
     allocate(mon_s(k)%fyic_s(maxi,maxj)) 
     allocate(mon_s(k)%Cor_v (maxi,maxj)) 
     allocate(mon_s(k)%PFv   (maxi,maxj)) 
     allocate(mon_s(k)%fwt(3,maxj))
     allocate(mon_s(k)%fwp(3,maxj))
     allocate(mon_s(k)%fwa(3,maxj))
    enddo


   return

  end subroutine sic_diag_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s i c _ d i a g
  !   Purpose    :  sea ice diagnostics
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_diag(sic)

    implicit none

    type(sic_class), intent(inout) :: sic

    integer :: i, j, n, m, k, y, JNS
    real(wp) :: sumv(2), suma(2), sume(2)
    character (len=256) :: fnm
    real(wp) :: mon_avg, ann_avg
    real(wp) :: fram_exp
    real(wp) :: denmark_exp
    real(wp) :: buoy_sic_NA(nlatv_buoy), S0, S_sic, S_star, beta
    real(wp) :: tv2, tv3
    real(wp) :: fwt(3,maxj), fwp(3,maxj), fwa(3,maxj)


    ! current index
    y = y_out_ts_clim

    mon_avg = 1._wp/nstep_mon_sic
    ann_avg = 1._wp/nstep_year_sic

    if (time_soy_sic) then
       mon_ts(y,:)%a_nh = 0._wp
       mon_ts(y,:)%a_sh = 0._wp
       mon_ts(y,:)%e_nh = 0._wp
       mon_ts(y,:)%e_sh = 0._wp
       mon_ts(y,:)%v_nh = 0._wp
       mon_ts(y,:)%v_sh = 0._wp
       ann_ts(y)%fram_exp = 0._wp
       ann_ts(y)%denmark_exp = 0._wp
       ann_ts(y)%buoy_sic_NA(:) = 0._wp
    endif
   
    sumv = 0._wp
    suma = 0._wp
    sume = 0._wp
    do j=1,maxj
       do i=1,maxi
          if (sic%f_ocn(i,j).gt.0._wp) then
             ! Northern/Southern hemisphere totals
             if (j.gt.(maxj/2)) then
                ! Northern hemisphere total volume, area and extent
                sumv(1) = sumv(1) + sic%h_sic_mean(i,j)*area(i,j) ! m3
                suma(1) = suma(1) + sic%f_sic(i,j)*area(i,j) ! m2
                if (sic%f_sic(i,j).gt.0.15_wp) sume(1) = sume(1) + area(i,j) ! m2
             else
                ! Southern hemisphere total volume, area and extent
                sumv(2) = sumv(2) + sic%h_sic_mean(i,j)*area(i,j) ! m3
                suma(2) = suma(2) + sic%f_sic(i,j)*area(i,j) ! m2
                if (sic%f_sic(i,j).gt.0.15_wp) sume(2) = sume(2) + area(i,j) ! m2
             endif
          endif
       enddo
    enddo

    mon_ts(y,mon)%a_nh = mon_ts(y,mon)%a_nh  +   suma(1)*1.e-12 * mon_avg  ! mln km2
    mon_ts(y,mon)%a_sh = mon_ts(y,mon)%a_sh  +   suma(2)*1.e-12 * mon_avg  ! mln km2
    mon_ts(y,mon)%e_nh = mon_ts(y,mon)%e_nh  +   sume(1)*1.e-12 * mon_avg  ! mln km2
    mon_ts(y,mon)%e_sh = mon_ts(y,mon)%e_sh  +   sume(2)*1.e-12 * mon_avg  ! mln km2
    mon_ts(y,mon)%v_nh = mon_ts(y,mon)%v_nh  +   sumv(1)*1.e-12 * mon_avg  ! 1e4 km3
    mon_ts(y,mon)%v_sh = mon_ts(y,mon)%v_sh  +   sumv(2)*1.e-12 * mon_avg  ! 1e4 km3

    ! sea ice export through Fram Strait (Sv) 
    fram_exp = 0.0
    do i=i_fram(1),i_fram(2)
      fram_exp = fram_exp + (sic%fay_sic(i,j_fram)+sic%fdy_sic(i,j_fram))/dt *1.e-6 ! m3/s * Sv/(m3/s) = Sv
    enddo
    ann_ts(y)%fram_exp = ann_ts(y)%fram_exp + fram_exp*ann_avg 

    ! sea ice export through denmark Strait (Sv) 
    denmark_exp = 0.0
    do i=i_denmark(1),i_denmark(2)
      denmark_exp = denmark_exp + (sic%fay_sic(i,j_denmark)+sic%fdy_sic(i,j_denmark))/dt *1.e-6 ! m/s*m*m * Sv/(m3/s) = Sv
    enddo
    ann_ts(y)%denmark_exp = ann_ts(y)%denmark_exp + denmark_exp*ann_avg 

    ! bouyancy flux from sea ice export through the southern boundary of the North Atlantic
    S0 = 34.7_wp    ! psu, reference salinity
    S_sic = 0._wp   ! psu, sea ice salinity
    S_star = S0-S_sic 
    beta = 0.8_wp   ! kg/m3/psu, haline contraction coefficient
    do n=1,nlatv_buoy
      JNS = minloc(abs(sv(:)-sin(pi*latv_buoy(n)/180.0)),1) + lbound(sv,1) - 1    ! -1 because index of sv array starts from 0!
      buoy_sic_NA(n) = 0.0
      do i=1,maxi
        if (lon(i).gt.-70. .and. basin_mask(i,JNS).eq.i_atlantic) then
          buoy_sic_NA(n) = buoy_sic_NA(n) &
            + g*beta*S_star*(sic%fay_sic(i,JNS)+sic%fdy_sic(i,JNS))  ! freshwater export, m/s2*kg/m3/psu*psu*m3 = kg*m/s2 = N
        endif
      enddo
    enddo
    ann_ts(y)%buoy_sic_NA(:) = ann_ts(y)%buoy_sic_NA(:) + buoy_sic_NA(:)
    sic%buoy_sic_NA(:) = buoy_sic_NA(:)


    if( time_eoy_sic ) then

     call ts_ave( mon_ts(y,:), ann_ts(y) )

     ! count number of active sea ice grid cells
     ann_ts(y)%ncells = 0
     ann_ts(y)%area = 0_wp
     do j=1,maxj
       do i=1,maxi
         if (sic%f_ocn(i,j).gt.0._wp) then
           ann_ts(y)%ncells = ann_ts(y)%ncells + 1
           ann_ts(y)%area   = ann_ts(y)%area + area(i,j) * 1.e-12_wp ! mln km2
         endif
       enddo
     enddo

     ! write to file
     if (time_out_ts_clim) then
       fnm = trim(out_dir)//"/sic_ts.nc"
       do k = 1, nmon_year
         call ts_nc_write(fnm,mon_ts(1:y,k),k,year_clim-y+1,y)
       end do
       call ts_nc_write(fnm,ann_ts(1:y),nmon_year+1,year_clim-y+1,y)
     endif

     ! print header
     if (mod(year,10).eq.1) then
        print '(a7,a9,6a7)','sic','year','NHmax','NHmin','SHmax','SHmin','FramExp','DenExp'
     endif

     ! print values
     print '(a7,i9,4F7.1,2F7.2)', &
     'sic',year_now,maxval(mon_ts(y,:)%a_nh),minval(mon_ts(y,:)%a_nh),maxval(mon_ts(y,:)%a_sh),minval(mon_ts(y,:)%a_sh), &
     -ann_ts(y)%fram_exp, -ann_ts(y)%denmark_exp

    endif


    ! spatially explicit output
    if ( time_out_sic ) then

     if( time_soy_sic ) then
      do m=1,nmon_year
       mon_s(m)%hsic  = 0._wp
       mon_s(m)%hsnow = 0._wp
       mon_s(m)%fsic  = 0._wp
       mon_s(m)%tsic  = 0._wp
       mon_s(m)%tocn  = 0._wp
       mon_s(m)%dh_sic_dt_therm = 0._wp
       mon_s(m)%dh_sic_dt_therm_sic = 0._wp
       mon_s(m)%dh_sic_dt_therm_ocn = 0._wp
       mon_s(m)%dh_sic_dt_dyn   = 0._wp
       mon_s(m)%flx_melt_top = 0._wp
       mon_s(m)%flx_melt_bot = 0._wp
       mon_s(m)%lh    = 0._wp
       mon_s(m)%sh    = 0._wp
       mon_s(m)%lw    = 0._wp
       mon_s(m)%fx    = 0._wp
       mon_s(m)%fw    = 0._wp
       mon_s(m)%fw_brines = 0._wp
       mon_s(m)%e     = 0._wp
       mon_s(m)%p_e   = 0._wp
       mon_s(m)%lh_sic    = 0._wp
       mon_s(m)%sh_sic    = 0._wp
       mon_s(m)%lw_sic    = 0._wp
       mon_s(m)%fx_sic    = 0._wp
       mon_s(m)%fw_sic    = 0._wp
       mon_s(m)%e_sic     = 0._wp
       mon_s(m)%lh_ocn    = 0._wp
       mon_s(m)%sh_ocn    = 0._wp
       mon_s(m)%lw_ocn    = 0._wp
       mon_s(m)%fx_ocn    = 0._wp
       mon_s(m)%fw_ocn    = 0._wp
       mon_s(m)%e_ocn     = 0._wp
       mon_s(m)%cde_ocn    = 0._wp
       mon_s(m)%cde_sic    = 0._wp
       mon_s(m)%cdh_ocn    = 0._wp
       mon_s(m)%cdh_sic    = 0._wp
       mon_s(m)%snow_grain = 0._wp
       mon_s(m)%dust_con   = 0._wp
       mon_s(m)%alb_ocn   = 0._wp
       mon_s(m)%alb_sic   = 0._wp
       mon_s(m)%rain     = 0._wp
       mon_s(m)%snow     = 0._wp
       mon_s(m)%usic     = 0._wp
       mon_s(m)%vsic     = 0._wp
       mon_s(m)%tauxa    = 0._wp
       mon_s(m)%tauya    = 0._wp
       mon_s(m)%tauxo    = 0._wp
       mon_s(m)%tauyo    = 0._wp
       mon_s(m)%str_d  = 0._wp  
       mon_s(m)%str_t  = 0._wp
       mon_s(m)%str_s  = 0._wp
       mon_s(m)%fxic   = 0._wp 
       mon_s(m)%fxic_d = 0._wp 
       mon_s(m)%fxic_t = 0._wp 
       mon_s(m)%fxic_s = 0._wp 
       mon_s(m)%Cor_u  = 0._wp 
       mon_s(m)%PFu    = 0._wp 
       mon_s(m)%fyic   = 0._wp 
       mon_s(m)%fyic_d = 0._wp 
       mon_s(m)%fyic_t = 0._wp 
       mon_s(m)%fyic_s = 0._wp 
       mon_s(m)%Cor_v  = 0._wp 
       mon_s(m)%PFv    = 0._wp 
       mon_s(m)%fwt    = 0._wp 
       mon_s(m)%fwa    = 0._wp 
       mon_s(m)%fwp    = 0._wp 
      enddo
     endif

     if (l_daily_output) then
       day_s(doy)%usic  = sic%u*100.  ! cm/s
       day_s(doy)%vsic  = sic%v*100.  ! cm/s
       day_s(doy)%hsic  = sic%h_sic_mean
       day_s(doy)%hsnow = sic%h_snow_mean   
       day_s(doy)%fsic  = sic%f_sic    
     endif

     where (sic%f_ocn.gt.0._wp) 
       mon_s(mon)%hsic  = mon_s(mon)%hsic   + sic%h_sic    * mon_avg
       mon_s(mon)%hsnow = mon_s(mon)%hsnow  + sic%h_snow   * mon_avg
       mon_s(mon)%fsic  = mon_s(mon)%fsic   + sic%f_sic    * mon_avg
       mon_s(mon)%tsic  = mon_s(mon)%tsic   + sic%t_skin_sic * mon_avg
       mon_s(mon)%tocn  = mon_s(mon)%tocn   + sic%t_skin_ocn * mon_avg
       mon_s(mon)%dh_sic_dt_therm= mon_s(mon)%dh_sic_dt_therm   + sic%dh_sic_dt_therm*sec_day    * mon_avg
       mon_s(mon)%dh_sic_dt_therm_sic= mon_s(mon)%dh_sic_dt_therm_sic   + sic%dh_sic_sic/dt*sec_day    * mon_avg
       mon_s(mon)%dh_sic_dt_therm_ocn= mon_s(mon)%dh_sic_dt_therm_ocn   + sic%dh_sic_ocn/dt*sec_day    * mon_avg
       mon_s(mon)%dh_sic_dt_dyn  = mon_s(mon)%dh_sic_dt_dyn     + sic%dh_sic_dt_dyn*sec_day      * mon_avg
       mon_s(mon)%flx_melt_top = mon_s(mon)%flx_melt_top  + sic%flx_melt_top   * mon_avg
       mon_s(mon)%flx_melt_bot = mon_s(mon)%flx_melt_bot  + sic%flx_melt_bot   * mon_avg
       mon_s(mon)%lh    = mon_s(mon)%lh     + (sic%f_sic*sic%flx_lh_sic + (1._wp-sic%f_sic)*sic%flx_lh_ocn)     * mon_avg
       mon_s(mon)%sh    = mon_s(mon)%sh     + (sic%f_sic*sic%flx_sh_sic + (1._wp-sic%f_sic)*sic%flx_sh_ocn)     * mon_avg
       mon_s(mon)%lw    = mon_s(mon)%lw     + (sic%f_sic*sic%flx_lwu_sic + (1._wp-sic%f_sic)*sic%flx_lwu_ocn)     * mon_avg
       mon_s(mon)%fx    = mon_s(mon)%fx     + sic%flx_ocn    * mon_avg
       mon_s(mon)%fw    = mon_s(mon)%fw     + sic%fw_ocn*sec_day     * mon_avg  ! kg/m2/day
       mon_s(mon)%fw_brines = mon_s(mon)%fw_brines + sic%fw_brines*sec_day     * mon_avg  ! kg/m2/day
       mon_s(mon)%e     = mon_s(mon)%e      + (sic%f_sic*sic%evp_sic + (1._wp-sic%f_sic)*sic%evp_ocn)*sec_day        * mon_avg ! kg/m2/day
       mon_s(mon)%p_e     = mon_s(mon)%p_e      + (sic%rain + sic%snow -(sic%f_sic*sic%evp_sic + (1._wp-sic%f_sic)*sic%evp_ocn))*sec_day        * mon_avg ! kg/m2/day
       mon_s(mon)%lh_sic    = mon_s(mon)%lh_sic     + sic%flx_lh_sic      * mon_avg
       mon_s(mon)%sh_sic    = mon_s(mon)%sh_sic     + sic%flx_sh_sic     * mon_avg
       mon_s(mon)%lw_sic    = mon_s(mon)%lw_sic     + sic%flx_lwu_sic    * mon_avg
       mon_s(mon)%fx_sic    = mon_s(mon)%fx_sic     + sic%flx_ocn_sic    * mon_avg
       mon_s(mon)%fw_sic    = mon_s(mon)%fw_sic     + sic%fw_ocn_sic*sec_day     * mon_avg  ! kg/m2/day
       mon_s(mon)%e_sic     = mon_s(mon)%e_sic      + sic%evp_sic*sec_day         * mon_avg ! kg/m2/day
       mon_s(mon)%lh_ocn    = mon_s(mon)%lh_ocn     + sic%flx_lh_ocn     * mon_avg
       mon_s(mon)%sh_ocn    = mon_s(mon)%sh_ocn     + sic%flx_sh_ocn     * mon_avg
       mon_s(mon)%lw_ocn    = mon_s(mon)%lw_ocn     + sic%flx_lwu_ocn     * mon_avg
       mon_s(mon)%fx_ocn    = mon_s(mon)%fx_ocn     + sic%flx_ocn_ocn    * mon_avg
       mon_s(mon)%fw_ocn    = mon_s(mon)%fw_ocn     + sic%fw_ocn_ocn*sec_day     * mon_avg  ! kg/m2/day
       mon_s(mon)%e_ocn     = mon_s(mon)%e_ocn      + sic%evp_ocn*sec_day        * mon_avg ! kg/m2/day
       mon_s(mon)%cde_ocn   = mon_s(mon)%cde_ocn    + sic%cde_ocn        * mon_avg 
       mon_s(mon)%cde_sic   = mon_s(mon)%cde_sic    + sic%cde_sic        * mon_avg 
       mon_s(mon)%cdh_ocn   = mon_s(mon)%cdh_ocn    + sic%cdh_ocn        * mon_avg 
       mon_s(mon)%cdh_sic   = mon_s(mon)%cdh_sic    + sic%cdh_sic        * mon_avg 
       mon_s(mon)%snow_grain= mon_s(mon)%snow_grain + sic%snow_grain        * mon_avg 
       mon_s(mon)%dust_con  = mon_s(mon)%dust_con   + sic%dust_con*1.e6     * mon_avg   ! mg/kg
       mon_s(mon)%alb_ocn   = mon_s(mon)%alb_ocn    + sic%albedo_ocn        * mon_avg 
       mon_s(mon)%alb_sic   = mon_s(mon)%alb_sic    + sic%albedo_sic        * mon_avg 
       mon_s(mon)%rain      = mon_s(mon)%rain       + sic%rain*sec_day        * mon_avg ! kg/m2/day
       mon_s(mon)%snow      = mon_s(mon)%snow       + sic%snow*sec_day        * mon_avg ! kg/m2/day
     endwhere
     mon_s(mon)%usic     = mon_s(mon)%usic      + sic%u*100.       * mon_avg ! cm/s
     mon_s(mon)%vsic     = mon_s(mon)%vsic      + sic%v*100.       * mon_avg ! cm/s
     mon_s(mon)%tauxa    = mon_s(mon)%tauxa     + sic%tauxa        * mon_avg ! N/m2
     mon_s(mon)%tauya    = mon_s(mon)%tauya     + sic%tauya        * mon_avg ! N/m2
     mon_s(mon)%tauxo    = mon_s(mon)%tauxo     + sic%tauxo        * mon_avg ! N/m2
     mon_s(mon)%tauyo    = mon_s(mon)%tauyo     + sic%tauyo        * mon_avg ! N/m2
     mon_s(mon)%str_d    = mon_s(mon)%str_d     + sic%str_d        * mon_avg 
     mon_s(mon)%str_t    = mon_s(mon)%str_t     + sic%str_t        * mon_avg
     mon_s(mon)%str_s    = mon_s(mon)%str_s     + sic%str_s        * mon_avg
     mon_s(mon)%fxic     = mon_s(mon)%fxic      + sic%fxic         * mon_avg
     mon_s(mon)%fxic_d   = mon_s(mon)%fxic_d    + sic%fxic_d       * mon_avg
     mon_s(mon)%fxic_t   = mon_s(mon)%fxic_t    + sic%fxic_t       * mon_avg
     mon_s(mon)%fxic_s   = mon_s(mon)%fxic_s    + sic%fxic_s       * mon_avg
     mon_s(mon)%Cor_u    = mon_s(mon)%Cor_u     + sic%Cor_u        * mon_avg
     mon_s(mon)%PFu      = mon_s(mon)%PFu       + sic%PFu          * mon_avg
     mon_s(mon)%fyic     = mon_s(mon)%fyic      + sic%fyic         * mon_avg
     mon_s(mon)%fyic_d   = mon_s(mon)%fyic_d    + sic%fyic_d       * mon_avg
     mon_s(mon)%fyic_t   = mon_s(mon)%fyic_t    + sic%fyic_t       * mon_avg
     mon_s(mon)%fyic_s   = mon_s(mon)%fyic_s    + sic%fyic_s       * mon_avg
     mon_s(mon)%Cor_v    = mon_s(mon)%Cor_v     + sic%Cor_v        * mon_avg
     mon_s(mon)%PFv      = mon_s(mon)%PFv       + sic%PFv          * mon_avg

     fwt = 0._wp
     fwa = 0._wp
     fwp = 0._wp
     do j=1,maxj
       do i=1,maxi
         ! advective transport
         tv2 = sic%fay_sic(i,j)/dt  ! m3/s
         ! diffusive transport
         tv3 = sic%fdy_sic(i,j)/dt  ! m3/s
         fwt(1,j) = fwt(1,j) + tv2 + tv3
         fwt(2,j) = fwt(2,j) + tv2
         fwt(3,j) = fwt(3,j) + tv3
         if (basin_mask2(i,j).eq.i_pacific .or. basin_mask2(i,j).eq.i_indian) then
           fwp(1,j) = fwp(1,j) + tv2 + tv3
           fwp(2,j) = fwp(2,j) + tv2
           fwp(3,j) = fwp(3,j) + tv3
         else if (basin_mask2(i,j).eq.i_atlantic) then
           fwa(1,j) = fwa(1,j) + tv2 + tv3
           fwa(2,j) = fwa(2,j) + tv2
           fwa(3,j) = fwa(3,j) + tv3
         endif
       enddo
     enddo

     mon_s(mon)%fwt = mon_s(mon)%fwt + fwt*1.e-6_wp * mon_avg ! Sv
     mon_s(mon)%fwa = mon_s(mon)%fwa + fwa*1.e-6_wp * mon_avg ! Sv
     mon_s(mon)%fwp = mon_s(mon)%fwp + fwp*1.e-6_wp * mon_avg ! Sv

   endif

    if (time_out_sic .and. time_eoy_sic) then
       ann_s%focn  = sic%f_ocn
       call sic_diag_out
    endif


   return

  end subroutine sic_diag


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s i c _ d i a g _ o u t
  ! Purpose  :  write sea ice netcdf output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_diag_out

    implicit none

    integer :: k, ncid
    character (len=256) :: fnm


    nout = nout +1 

    ! Get annual values
    call sic_ave( mon_s,ann_s )

    ! write to file
    fnm = trim(out_dir)//"/sic.nc"
    call nc_open(fnm,ncid)
    call nc_write(fnm,dim_time,real(year_now,wp), dim1=dim_time, start=[nout], count=[1],ncid=ncid)    
    do k = 1, nmon_year
       call sic_nc_write(fnm,ncid,mon_s(k),k,nout)
    end do
    call sic_nc_write(fnm,ncid,ann_s,nmon_year+1,nout)
    call nc_close(ncid)

    if (l_daily_output) then
      ! write to file
      fnm = trim(out_dir)//"/sic_daily.nc"
      call nc_open(fnm,ncid)
      call nc_write(fnm,dim_time,real(year_now,wp), dim1=dim_time, start=[nout], count=[1],ncid=ncid)    
      do k = 1, nday_year
        call sic_daily_nc_write(fnm,ncid,day_s(k),k,nout)
      end do
      call nc_close(ncid)
    endif


   return

  end subroutine sic_diag_out
  

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t s _ n c
  ! Purpose  :  initialize netcdf file for time series output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ts_nc(fnm)

    implicit none

    character (len=*) :: fnm
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_write_dim(fnm, dim_time, x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE.)
    call nc_write_dim(fnm, dim_month, x=1._wp, dx=1._wp, nx=13, axis="e", &
    units="months")
    call nc_write_dim(fnm, dim_lat, x=latv_buoy, axis="y", units="degN")
    call nc_write_dim(fnm, dim_lon, x=1, axis="x", units="1")

    return

  end subroutine ts_nc
 

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t s _ n c _ w r i t e
  ! Purpose  :  write time series to netcdf
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ts_nc_write(fnm,vars,ndat,nout,y)

    implicit none

    type(ts_out) :: vars(:)

    character (len=*) :: fnm
    integer :: ndat, nout, n, y, ncid, i

    call nc_open(fnm,ncid)
    call nc_write(fnm,"time", real([(i,i=(year_now-(y-1)*n_accel),(year_now),(n_accel))],wp), &
    dim1=dim_time,start=[nout],count=[y],ncid=ncid)    
    if (ndat.eq.13) then
       call nc_write(fnm,"ncells", vars%ncells, dims=[dim_time],start=[nout],count=[y],long_name="number of active sea ice grid cells",units="1",ncid=ncid)
       call nc_write(fnm,"area", vars%area, dims=[dim_time],start=[nout],count=[y],long_name="area of sea ice domain",units="mln km2",ncid=ncid)
       call nc_write(fnm,"a_nh_min", vars%a_nh_min, dims=[dim_time],start=[nout],count=[y],long_name="minimum NH sea ice area (4.4 NSIDC)",units="mln km2",ncid=ncid)
       call nc_write(fnm,"a_nh_max", vars%a_nh_max, dims=[dim_time],start=[nout],count=[y],long_name="maximum NH sea ice area (13 NSIDC)",units="mln km2",ncid=ncid)
       call nc_write(fnm,"a_sh_min", vars%a_sh_min, dims=[dim_time],start=[nout],count=[y],long_name="minimum SH sea ice area (2 NSIDC)",units="mln km2",ncid=ncid)
       call nc_write(fnm,"a_sh_max", vars%a_sh_max, dims=[dim_time],start=[nout],count=[y],long_name="maximum SH sea ice area (14.5 NSIDC)",units="mln km2",ncid=ncid)
       call nc_write(fnm,"e_nh_min", vars%e_nh_min, dims=[dim_time],start=[nout],count=[y],long_name="minimum NH sea ice extent ",units="mln km2",ncid=ncid)
       call nc_write(fnm,"e_nh_max", vars%e_nh_max, dims=[dim_time],start=[nout],count=[y],long_name="maximum NH sea ice extent ",units="mln km2",ncid=ncid)
       call nc_write(fnm,"e_sh_min", vars%e_sh_min, dims=[dim_time],start=[nout],count=[y],long_name="minimum SH sea ice extent ",units="mln km2",ncid=ncid)
       call nc_write(fnm,"e_sh_max", vars%e_sh_max, dims=[dim_time],start=[nout],count=[y],long_name="maximum SH sea ice extent ",units="mln km2",ncid=ncid)
       call nc_write(fnm,"v_nh_min", vars%v_nh_min, dims=[dim_time],start=[nout],count=[y],long_name="minimum NH sea ice volume ",units="10^4 km3",ncid=ncid)
       call nc_write(fnm,"v_nh_max", vars%v_nh_max, dims=[dim_time],start=[nout],count=[y],long_name="maximum NH sea ice volume ",units="10^4 km3",ncid=ncid)
       call nc_write(fnm,"v_sh_min", vars%v_sh_min, dims=[dim_time],start=[nout],count=[y],long_name="minimum SH sea ice volume ",units="10^4 km3",ncid=ncid)
       call nc_write(fnm,"v_sh_max", vars%v_sh_max, dims=[dim_time],start=[nout],count=[y],long_name="maximum SH sea ice volume ",units="10^4 km3",ncid=ncid)
       call nc_write(fnm,"fram_exp", vars%fram_exp, dims=[dim_time],start=[nout],count=[y], &
         long_name="sea ice export through the Fram strait, positive northward (-0.1 Sv)",units="Sv",ncid=ncid)
       call nc_write(fnm,"denmark_exp", vars%denmark_exp, dims=[dim_time],start=[nout],count=[y], &
         long_name="sea ice export through the denmark strait, positive northward (-0.1 Sv)",units="Sv",ncid=ncid)
       call nc_write(fnm,"buoy_sic_NA55", vars%buoy_sic_NA(ilatv_buoy_sel), dims=[dim_time],start=[nout],count=[y], &
         long_name="Bouyancy flux from export of sea ice through southern border of North Atlantic",units="N",ncid=ncid)
       do n=1,nlatv_buoy
         call nc_write(fnm,"buoy_sic_NA",    vars%buoy_sic_NA(n),dim1=dim_lat,dim2=dim_time,start=[n,nout],count=[1,y],long_name="Bouyancy flux from export of sea ice through southern border of North Atlantic",units="N",ncid=ncid)
       enddo
    endif
    call nc_write(fnm,"a_nh", vars%a_nh, dims=[dim_month,dim_time],start=[ndat,nout],count=[1,y],long_name="monthly NH sea ice area",units="mln km2",ncid=ncid)
    call nc_write(fnm,"a_sh", vars%a_sh, dims=[dim_month,dim_time],start=[ndat,nout],count=[1,y],long_name="monthly SH sea ice area",units="mln km2",ncid=ncid)
    call nc_write(fnm,"e_nh", vars%e_nh, dims=[dim_month,dim_time],start=[ndat,nout],count=[1,y],long_name="monthly NH sea ice extent (>15 percent)",units="mln km2",ncid=ncid)
    call nc_write(fnm,"e_sh", vars%e_sh, dims=[dim_month,dim_time],start=[ndat,nout],count=[1,y],long_name="monthly SH sea ice extent (>15 percent)",units="mln km2",ncid=ncid)
    call nc_write(fnm,"v_nh", vars%v_nh, dims=[dim_month,dim_time],start=[ndat,nout],count=[1,y],long_name="monthly NH sea ice volume",units="10^4 km3",ncid=ncid)
    call nc_write(fnm,"v_sh", vars%v_sh, dims=[dim_month,dim_time],start=[ndat,nout],count=[1,y],long_name="monthly SH sea ice volume",units="10^4 km3",ncid=ncid)
    call nc_close(ncid)


   return

  end subroutine ts_nc_write


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t s _ a v e
  ! Purpose  :  Average (or sum) the time series as needed
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ts_ave(d,ave)

    implicit none

    type(ts_out) :: d(:), ave

    integer :: k, n
    real(wp) :: div

    n = size(d)
    div = real(n,wp)

    ! Set all values to zero
    ave%a_nh = 0._wp
    ave%a_sh = 0._wp
    ave%e_nh = 0._wp
    ave%e_sh = 0._wp
    ave%v_nh = 0._wp
    ave%v_sh = 0._wp

    do k = 1, n
     ave%a_nh = ave%a_nh + d(k)%a_nh / div
     ave%a_sh = ave%a_sh + d(k)%a_sh / div
     ave%e_nh = ave%e_nh + d(k)%e_nh / div
     ave%e_sh = ave%e_sh + d(k)%e_sh / div
     ave%v_nh = ave%v_nh + d(k)%v_nh / div
     ave%v_sh = ave%v_sh + d(k)%v_sh / div
    end do

    ave%a_nh_min = minval(d(:)%a_nh)
    ave%a_nh_max = maxval(d(:)%a_nh)
    ave%a_sh_min = minval(d(:)%a_sh)
    ave%a_sh_max = maxval(d(:)%a_sh)
    ave%e_nh_min = minval(d(:)%e_nh)
    ave%e_nh_max = maxval(d(:)%e_nh)
    ave%e_sh_min = minval(d(:)%e_sh)
    ave%e_sh_max = maxval(d(:)%e_sh)
    ave%v_nh_min = minval(d(:)%v_nh)
    ave%v_nh_max = maxval(d(:)%v_nh)
    ave%v_sh_min = minval(d(:)%v_sh)
    ave%v_sh_max = maxval(d(:)%v_sh)


   return

  end subroutine ts_ave


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s i c _ n c
  ! Purpose  :  Initialize sea ice netcdf output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_nc(fnm)

    implicit none

    character (len=*) :: fnm
    integer :: ncid
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,dim_time, x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE., ncid=ncid)
    call nc_write_dim(fnm, dim_month, x=1._wp, dx=1._wp, nx=13, axis="e", &
    units="months", ncid=ncid)
    call nc_write_dim(fnm, dim_lon, x=lon, axis="x", ncid=ncid)
    call nc_write_dim(fnm, dim_lat, x=lat, axis="y", ncid=ncid)
    call nc_write_dim(fnm,"lonu",x=lonu(2:maxi+1),axis="y",ncid=ncid)
    call nc_write_dim(fnm,"latv",x=latv(2:maxj+1),axis="y",ncid=ncid)
    call nc_write_dim(fnm,"type",x=1._wp,dx=1._wp,nx=3,units="[tot,adv,diff]",ncid=ncid)
    call nc_close(ncid)

   return

  end subroutine sic_nc


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s i c _ n c _ w r i t e
  ! Purpose  :  Output of sea ice netcdf
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_nc_write(fnm,ncid,vars,ndat,nout)

    implicit none

    type(s_out) :: vars

    character (len=*) :: fnm
    integer :: ndat, nout, ncid

    character(len=4), parameter :: dim_lonu = "lonu"
    character(len=4), parameter :: dim_latv = "latv"
    
    if (ndat.eq.13) then
    call nc_write(fnm,"focn",  sngl(vars%focn),  dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[maxi,maxj,1],long_name="ocean fraction",units="/",ncid=ncid)
    endif

    call nc_write(fnm,"hsic",  sngl(vars%hsic),  dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="sea ice thickness of sea ice fraction",units="m",ncid=ncid)
    call nc_write(fnm,"hsnow", sngl(vars%hsnow), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="snow thickness of sea ice fraction",units="m",ncid=ncid)
    call nc_write(fnm,"fsic",  sngl(vars%fsic),  dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="sea ice fraction",units="/",ncid=ncid)
    call nc_write(fnm,"tsic",  sngl(vars%tsic),  dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="sea ice skin temperature",units="�C",ncid=ncid)
    call nc_write(fnm,"tocn",  sngl(vars%tocn),  dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="ocean skin temperature",units="�C",ncid=ncid)
    call nc_write(fnm,"dh_sic_dt_therm",  sngl(vars%dh_sic_dt_therm),  dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="sea ice thickness change from thermodynamics",units="m/day",ncid=ncid)
    call nc_write(fnm,"dh_sic_dt_therm_sic",  sngl(vars%dh_sic_dt_therm_sic),  dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="sea ice thickness change from thermodynamics over sea ice area",units="m/day",ncid=ncid)
    call nc_write(fnm,"dh_sic_dt_therm_ocn",  sngl(vars%dh_sic_dt_therm_ocn),  dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="sea ice thickness change from thermodynamics in leads",units="m/day",ncid=ncid)
    call nc_write(fnm,"dh_sic_dt_dyn",  sngl(vars%dh_sic_dt_dyn),  dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="sea ice thickness change from transport",units="m/day",ncid=ncid)
    call nc_write(fnm,"flx_melt_top", sngl(vars%flx_melt_top), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="snow/ice melt flux from top",units="W/m2",ncid=ncid)
    call nc_write(fnm,"flx_melt_bot", sngl(vars%flx_melt_bot), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="snow/ice melt flux from bottom",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lh",    sngl(vars%lh),    dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="latent heat flux",units="W/m2",ncid=ncid)
    call nc_write(fnm,"sh",    sngl(vars%sh),    dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="sensible heat flux",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lwu",   sngl(vars%lw),    dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="upwelling longwave radiation",units="W/m2",ncid=ncid)
    call nc_write(fnm,"fx",    sngl(vars%fx),    dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="heat flux to the ocean",units="W/m2",ncid=ncid)
    call nc_write(fnm,"fw",    sngl(vars%fw),    dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="sea ice freshwater flux to the ocean",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"fw_brines",    sngl(vars%fw_brines),    dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="freshwater flux to the ocean from brine rejection",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"evp",   sngl(vars%e),     dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="evaporation/sublimation",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"p-e",   sngl(vars%p_e),   dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="precipitation - evaporation",units="kg/m2/day",ncid=ncid)

    call nc_write(fnm,"lh_sic",    sngl(vars%lh_sic),    dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="latent heat flux over sea ice",units="W/m2",ncid=ncid)
    call nc_write(fnm,"sh_sic",    sngl(vars%sh_sic),    dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="sensible heat flux over sea ice",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lwu_sic",   sngl(vars%lw_sic),    dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="upwelling longwave radiation over sea ice",units="W/m2",ncid=ncid)
    call nc_write(fnm,"fx_sic",    sngl(vars%fx_sic),    dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="heat flux to the ocean from sea ice fraction",units="W/m2",ncid=ncid)
    call nc_write(fnm,"fw_sic",    sngl(vars%fw_sic),    dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="sea ice freshwater flux to the ocean from sea ice fraction",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"evp_sic",   sngl(vars%e_sic),     dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="sublimation from sea ice",units="kg/m2/day",ncid=ncid)

    call nc_write(fnm,"lh_ocn",    sngl(vars%lh_ocn),    dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="latent heat flux from ocean",units="W/m2",ncid=ncid)
    call nc_write(fnm,"sh_ocn",    sngl(vars%sh_ocn),    dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="sensible heat flux from ocean",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lwu_ocn",   sngl(vars%lw_ocn),    dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="upwelling longwave radiation from ocean",units="W/m2",ncid=ncid)
    call nc_write(fnm,"fx_ocn",    sngl(vars%fx_ocn),    dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="heat flux to the ocean from ice free fraction",units="W/m2",ncid=ncid)
    call nc_write(fnm,"fw_ocn",    sngl(vars%fw_ocn),    dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="sea ice freshwater flux to the ocean from ice free fraction",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"evp_ocn",   sngl(vars%e_ocn),     dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="evaporation from ocean",units="kg/m2/day",ncid=ncid)

    call nc_write(fnm,"cde_ocn",   sngl(vars%cde_ocn),   dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="drag coefficient for moisture over ocean water",units="/",ncid=ncid)
    call nc_write(fnm,"cde_sic",   sngl(vars%cde_sic),   dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="drag coefficient for moisture over sea ice",units="/",ncid=ncid)
    call nc_write(fnm,"cdh_ocn",   sngl(vars%cdh_ocn),   dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="drag coefficient for heat over ocean water",units="/",ncid=ncid)
    call nc_write(fnm,"cdh_sic",   sngl(vars%cdh_sic),   dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="drag coefficient for heat over sea ice",units="/",ncid=ncid)

    call nc_write(fnm,"snow_grain",sngl(vars%snow_grain),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="snow grain size",units="micro m",ncid=ncid)
    call nc_write(fnm,"dust_con",  sngl(vars%dust_con),  dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="dust concentration in snow",units="mg/kg",ncid=ncid)
    call nc_write(fnm,"alb_ocn",   sngl(vars%alb_ocn),   dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="ocean albedo",units="/",ncid=ncid)
    call nc_write(fnm,"alb_sic",   sngl(vars%alb_sic),   dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="sea ice albedo",units="/",ncid=ncid)

    call nc_write(fnm,"rain",   sngl(vars%rain),   dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="rainfall",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"snow",   sngl(vars%snow),   dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="snowfall",units="kg/m2/day",ncid=ncid)

    call nc_write(fnm,"usic",    sngl(vars%usic),   dims=[dim_lonu,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="zonal sea ice drift on u-grid",units="cm/s",ncid=ncid)
    call nc_write(fnm,"vsic",    sngl(vars%vsic),   dims=[dim_lon,dim_latv,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="meridional sea ice drift on v-grid",units="cm/s",ncid=ncid)
    call nc_write(fnm,"tauxa",   sngl(vars%tauxa),  dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="zonal wind stress on sea ice",units="N/m2",ncid=ncid)
    call nc_write(fnm,"tauya",   sngl(vars%tauya),  dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="meridional wind stress on sea ice",units="N/m2",ncid=ncid)
    call nc_write(fnm,"tauxo",   sngl(vars%tauxo),  dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="zonal ocean stress on sea ice",units="N/m2",ncid=ncid)
    call nc_write(fnm,"tauyo",   sngl(vars%tauyo),  dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="meridional ocean stress on sea ice",units="N/m2",ncid=ncid)
    call nc_write(fnm,"str_d ", sngl(vars%str_d  ), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="The divergence stress tensor component",units="Pa m",ncid=ncid) 
    call nc_write(fnm,"str_t ", sngl(vars%str_t  ), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="The tension stress tensor component",units="Pa m",ncid=ncid)
    call nc_write(fnm,"str_s ", sngl(vars%str_s  ), dims=[dim_lonu,dim_latv,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="The shearing stress tensor component",units="Pa m",ncid=ncid)

    call nc_write(fnm,"fwt   ", sngl(vars%fwt    ), dims=["type",dim_latv,dim_month,dim_time],start=[1,1,ndat,nout],count=[3,maxj,1,1],long_name="Meridional freshwater transport by sea ice",units="Sv",ncid=ncid)
    call nc_write(fnm,"fwa   ", sngl(vars%fwa    ), dims=["type",dim_latv,dim_month,dim_time],start=[1,1,ndat,nout],count=[3,maxj,1,1],long_name="Meridional freshwater transport by sea ice in the Atlantic",units="Sv",ncid=ncid)
    call nc_write(fnm,"fwp   ", sngl(vars%fwp    ), dims=["type",dim_latv,dim_month,dim_time],start=[1,1,ndat,nout],count=[3,maxj,1,1],long_name="Meridional freshwater transport by sea ice in the Indo-Pacific",units="Sv",ncid=ncid)

    if (l_diag_dyn) then
    call nc_write(fnm,"fxic  ", sngl(vars%fxic   ), dims=[dim_lonu,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="Zonal force due to internal stresses",units="Pa",ncid=ncid)
    call nc_write(fnm,"fxic_d", sngl(vars%fxic_d ), dims=[dim_lonu,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="Zonal force due to divergence internal stress",units="Pa",ncid=ncid)
    call nc_write(fnm,"fxic_t", sngl(vars%fxic_t ), dims=[dim_lonu,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="Zonal force due to tension internal stress",units="Pa",ncid=ncid)
    call nc_write(fnm,"fxic_s", sngl(vars%fxic_s ), dims=[dim_lonu,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="Zonal force due to shearing internal stress",units="Pa",ncid=ncid)
    call nc_write(fnm,"Cor_u ", sngl(vars%Cor_u  ), dims=[dim_lonu,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="Zonal Coriolis acceleration",units="m/s2",ncid=ncid)
    call nc_write(fnm,"PFu   ", sngl(vars%PFu    ), dims=[dim_lonu,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="Zonal hydrostatic pressure driven acceleration",units="m/s2",ncid=ncid)
    call nc_write(fnm,"fyic  ", sngl(vars%fyic   ), dims=[dim_lon,dim_latv,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="Meridional force due to internal stresses",units="Pa",ncid=ncid)
    call nc_write(fnm,"fyic_d", sngl(vars%fyic_d ), dims=[dim_lon,dim_latv,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="Meridional force due to divergence internal stress",units="Pa",ncid=ncid)
    call nc_write(fnm,"fyic_t", sngl(vars%fyic_t ), dims=[dim_lon,dim_latv,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="Meridional force due to tension internal stress",units="Pa",ncid=ncid)
    call nc_write(fnm,"fyic_s", sngl(vars%fyic_s ), dims=[dim_lon,dim_latv,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="Meridional force due to shearing internal stress",units="Pa",ncid=ncid)
    call nc_write(fnm,"Cor_v ", sngl(vars%Cor_v  ), dims=[dim_lon,dim_latv,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="Meridional Coriolis acceleration",units="m/s2",ncid=ncid)
    call nc_write(fnm,"PFv   ", sngl(vars%PFv    ), dims=[dim_lon,dim_latv,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="Meridional hydrostatic pressure driven acceleration",units="m/s2",ncid=ncid)
    endif

   return

  end subroutine sic_nc_write


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s i c _ d a i l y _ n c
  ! Purpose  :  Initialize sea ice daily netcdf output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_daily_nc(fnm)

    implicit none

    character (len=*) :: fnm
    integer :: ncid
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,dim_time, x=empty_time, axis="t", units="years BP", unlimited=.TRUE.,ncid=ncid)
    call nc_write_dim(fnm,dim_lon,x=lon,axis="x",ncid=ncid)
    call nc_write_dim(fnm,dim_lat,x=lat,axis="y",ncid=ncid)
    call nc_write_dim(fnm, dim_day, x=1._wp, dx=1._wp, nx=nday_year, units="doy", axis="e", ncid=ncid)
    call nc_write_dim(fnm,"lonu",x=lonu(2:maxi+1),axis="y",ncid=ncid)
    call nc_write_dim(fnm,"latv",x=latv(2:maxj+1),axis="y",ncid=ncid)
    call nc_close(ncid)

   return

  end subroutine sic_daily_nc


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s i c _ d a i l y _ n c _ w r i t e
  ! Purpose  :  Output of daily sea ice netcdf
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_daily_nc_write(fnm,ncid,vars,ndat,nout)

    implicit none

    type(s_out) :: vars

    character (len=*) :: fnm
    integer :: ndat, nout, ncid

    character(len=4), parameter :: dim_lonu = "lonu"
    character(len=4), parameter :: dim_latv = "latv"
    
    call nc_write(fnm,"hsic",  sngl(vars%hsic),  dims=[dim_lon,dim_lat,dim_day,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="mean gridcell ice thickness",units="m",ncid=ncid)
    call nc_write(fnm,"hsnow", sngl(vars%hsnow), dims=[dim_lon,dim_lat,dim_day,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="mean gridcell snow thickness",units="m",ncid=ncid)
    call nc_write(fnm,"fsic",  sngl(vars%fsic),  dims=[dim_lon,dim_lat,dim_day,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="sea ice fraction",units="/",ncid=ncid)
    call nc_write(fnm,"usic",   sngl(vars%usic), dims=[dim_lonu,dim_lat,dim_day,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="zonal sea ice drift",units="cm/s",ncid=ncid)
    call nc_write(fnm,"vsic",   sngl(vars%vsic), dims=[dim_lon,dim_latv,dim_day,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="meridional sea ice drift",units="cm/s",ncid=ncid)

   return

  end subroutine sic_daily_nc_write


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s i c _ a v e
  ! Purpose  :  Average (or sum) the sea ice fields
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_ave(d,ave)

    implicit none

    type(s_out) :: d(:), ave

    integer :: k, n
    real(wp) :: div

    n = size(d)
    div = real(n,wp)

    ! Set all values to zero
    ave%hsic  = 0._wp
    ave%hsnow = 0._wp
    ave%fsic  = 0._wp
    ave%tsic  = 0._wp
    ave%tocn  = 0._wp
    ave%dh_sic_dt_therm= 0._wp
    ave%dh_sic_dt_therm_sic= 0._wp
    ave%dh_sic_dt_therm_ocn= 0._wp
    ave%dh_sic_dt_dyn  = 0._wp
    ave%flx_melt_top = 0._wp
    ave%flx_melt_bot = 0._wp
    ave%lh    = 0._wp
    ave%sh    = 0._wp
    ave%lw    = 0._wp
    ave%fx    = 0._wp
    ave%fw    = 0._wp
    ave%fw_brines = 0._wp
    ave%e     = 0._wp
    ave%p_e   = 0._wp
    ave%lh_sic    = 0._wp
    ave%sh_sic    = 0._wp
    ave%lw_sic    = 0._wp
    ave%fx_sic    = 0._wp
    ave%fw_sic    = 0._wp
    ave%e_sic     = 0._wp
    ave%lh_ocn    = 0._wp
    ave%sh_ocn    = 0._wp
    ave%lw_ocn    = 0._wp
    ave%fx_ocn    = 0._wp
    ave%fw_ocn    = 0._wp
    ave%e_ocn     = 0._wp
    ave%cde_ocn    = 0._wp
    ave%cde_sic    = 0._wp
    ave%cdh_ocn    = 0._wp
    ave%cdh_sic    = 0._wp
    ave%snow_grain = 0._wp
    ave%dust_con   = 0._wp
    ave%alb_ocn    = 0._wp
    ave%alb_sic    = 0._wp
    ave%rain     = 0._wp
    ave%snow     = 0._wp
    ave%usic     = 0._wp
    ave%vsic     = 0._wp
    ave%tauxa    = 0._wp
    ave%tauya    = 0._wp
    ave%tauxo    = 0._wp
    ave%tauyo    = 0._wp
    ave%str_d    = 0._wp 
    ave%str_t    = 0._wp
    ave%str_s    = 0._wp
    ave%fxic     = 0._wp
    ave%fxic_d   = 0._wp
    ave%fxic_t   = 0._wp
    ave%fxic_s   = 0._wp
    ave%Cor_u    = 0._wp
    ave%PFu      = 0._wp
    ave%fyic     = 0._wp
    ave%fyic_d   = 0._wp
    ave%fyic_t   = 0._wp
    ave%fyic_s   = 0._wp
    ave%Cor_v    = 0._wp
    ave%PFv      = 0._wp
    ave%fwt      = 0._wp
    ave%fwa      = 0._wp
    ave%fwp      = 0._wp

    ! Loop over the time indices to sum up and average (if necessary)
    do k = 1, n
       ave%hsic  = ave%hsic   + d(k)%hsic   / div
       ave%hsnow = ave%hsnow  + d(k)%hsnow  / div
       ave%fsic  = ave%fsic   + d(k)%fsic   / div
       ave%tsic  = ave%tsic   + d(k)%tsic   / div
       ave%tocn  = ave%tocn   + d(k)%tocn   / div
       ave%dh_sic_dt_therm= ave%dh_sic_dt_therm   + d(k)%dh_sic_dt_therm   / div
       ave%dh_sic_dt_therm_sic= ave%dh_sic_dt_therm_sic   + d(k)%dh_sic_dt_therm_sic   / div
       ave%dh_sic_dt_therm_ocn= ave%dh_sic_dt_therm_ocn   + d(k)%dh_sic_dt_therm_ocn   / div
       ave%dh_sic_dt_dyn  = ave%dh_sic_dt_dyn     + d(k)%dh_sic_dt_dyn     / div
       ave%flx_melt_top = ave%flx_melt_top + d(k)%flx_melt_top  / div
       ave%flx_melt_bot = ave%flx_melt_bot + d(k)%flx_melt_bot  / div
       ave%lh    = ave%lh     + d(k)%lh     / div
       ave%sh    = ave%sh     + d(k)%sh     / div
       ave%lw    = ave%lw     + d(k)%lw     / div
       ave%fx    = ave%fx     + d(k)%fx     / div
       ave%fw    = ave%fw     + d(k)%fw     / div
       ave%fw_brines = ave%fw_brines + d(k)%fw_brines     / div
       ave%e     = ave%e      + d(k)%e      / div
       ave%p_e   = ave%p_e    + d(k)%p_e    / div
       ave%lh_sic    = ave%lh_sic     + d(k)%lh_sic     / div
       ave%sh_sic    = ave%sh_sic     + d(k)%sh_sic     / div
       ave%lw_sic    = ave%lw_sic     + d(k)%lw_sic     / div
       ave%fx_sic    = ave%fx_sic     + d(k)%fx_sic     / div
       ave%fw_sic    = ave%fw_sic     + d(k)%fw_sic     / div
       ave%e_sic     = ave%e_sic      + d(k)%e_sic      / div
       ave%lh_ocn    = ave%lh_ocn     + d(k)%lh_ocn     / div
       ave%sh_ocn    = ave%sh_ocn     + d(k)%sh_ocn     / div
       ave%lw_ocn    = ave%lw_ocn     + d(k)%lw_ocn     / div
       ave%fx_ocn    = ave%fx_ocn     + d(k)%fx_ocn     / div
       ave%fw_ocn    = ave%fw_ocn     + d(k)%fw_ocn     / div
       ave%e_ocn     = ave%e_ocn      + d(k)%e_ocn      / div
       ave%cde_ocn    = ave%cde_ocn   + d(k)%cde_ocn    / div
       ave%cde_sic    = ave%cde_sic   + d(k)%cde_sic    / div
       ave%cdh_ocn    = ave%cdh_ocn   + d(k)%cdh_ocn    / div
       ave%cdh_sic    = ave%cdh_sic   + d(k)%cdh_sic    / div
       ave%snow_grain = ave%snow_grain+ d(k)%snow_grain / div
       ave%dust_con   = ave%dust_con  + d(k)%dust_con   / div
       ave%alb_ocn    = ave%alb_ocn   + d(k)%alb_ocn    / div
       ave%alb_sic    = ave%alb_sic   + d(k)%alb_sic    / div
       ave%rain     = ave%rain      + d(k)%rain      / div
       ave%snow     = ave%snow      + d(k)%snow      / div
       ave%usic     = ave%usic      + d(k)%usic      / div
       ave%vsic     = ave%vsic      + d(k)%vsic      / div
       ave%tauxa    = ave%tauxa     + d(k)%tauxa     / div
       ave%tauya    = ave%tauya     + d(k)%tauya     / div
       ave%tauxo    = ave%tauxo     + d(k)%tauxo     / div
       ave%tauyo    = ave%tauyo     + d(k)%tauyo     / div
       ave%str_d    = ave%str_d   + d(k)%str_d       / div 
       ave%str_t    = ave%str_t   + d(k)%str_t       / div
       ave%str_s    = ave%str_s   + d(k)%str_s       / div
       ave%fxic     = ave%fxic    + d(k)%fxic        / div
       ave%fxic_d   = ave%fxic_d  + d(k)%fxic_d      / div
       ave%fxic_t   = ave%fxic_t  + d(k)%fxic_t      / div
       ave%fxic_s   = ave%fxic_s  + d(k)%fxic_s      / div
       ave%Cor_u    = ave%Cor_u   + d(k)%Cor_u       / div
       ave%PFu      = ave%PFu     + d(k)%PFu         / div
       ave%fyic     = ave%fyic    + d(k)%fyic        / div
       ave%fyic_d   = ave%fyic_d  + d(k)%fyic_d      / div
       ave%fyic_t   = ave%fyic_t  + d(k)%fyic_t      / div
       ave%fyic_s   = ave%fyic_s  + d(k)%fyic_s      / div
       ave%Cor_v    = ave%Cor_v   + d(k)%Cor_v       / div
       ave%PFv      = ave%PFv     + d(k)%PFv         / div
       ave%fwt      = ave%fwt     + d(k)%fwt         / div
       ave%fwa      = ave%fwa     + d(k)%fwa         / div
       ave%fwp      = ave%fwp     + d(k)%fwp         / div
    end do

   return

  end subroutine sic_ave

end module sic_out
