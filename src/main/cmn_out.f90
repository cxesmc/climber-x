!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c m n _ o u t
!
!  Purpose : common grid diagnostics and output
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
module cmn_out
  
  use precision, only : wp
  use dim_name, only: dim_time, dim_lon, dim_lat, dim_month, dim_nsurf
  use timer, only : n_accel, doy, mon, year, year_clim, sec_day, sec_mon, sec_year, day_mon, &
  nday_year, nmon_year, year_now, year_out_start
  use timer, only : time_out_cmn, time_soy, time_eoy, time_eom, nstep_mon, nyout_cmn, &
  ny_out_ts, y_out_ts_clim, time_out_ts_clim
  use control, only : out_dir
  use constants, only : T0, Lf
  use climber_grid, only : ni, nj, area, lon, lat
  use climber_grid, only : basin_mask, i_atlantic, i_pacific, i_indian, i_southern
  use coupler, only : nsurf 

  use ncio
  
  implicit none

  integer :: nout

  type ts_out
    real(wp) :: t2m
    real(wp) :: lh, sh, lwd, lwu, swn, lwn, lsnow, ebal
    real(wp) :: lh_l, sh_l, lwd_l, lwu_l, swn_l, lwn_l, lsnow_l, ebal_l
    real(wp) :: lh_o, sh_o, lwd_o, lwu_o, swn_o, lwn_o, lsnow_o, ebal_o
    real(wp) :: prc, prc_atl, prc_pac, prc_ind, prc_so
    real(wp) :: evp, evp_atl, evp_pac, evp_ind, evp_so
    real(wp) :: runoff, runoff_veg, runoff_ice, runoff_atl, runoff_pac, runoff_ind, runoff_so
    real(wp) :: calving, calving_veg, calving_ice, calving_atl, calving_pac, calving_ind, calving_so
    real(wp) :: bmelt, bmelt_grd, bmelt_flt, bmelt_atl, bmelt_pac, bmelt_ind, bmelt_so
    real(wp) :: fw_atl, fw_natl, fw_pac, fw_ind, fw_so
  end type

  type sur_out
    integer, dimension(:,:), allocatable :: mask_ice
    integer, dimension(:,:), allocatable :: mask_smb
    real(wp), dimension(:,:,:), allocatable :: f_stp
    real(wp), dimension(:,:), allocatable :: f_ocn, f_lnd, f_ice
    real(wp), dimension(:,:), allocatable :: prc, rain, snow, wind
    real(wp), dimension(:,:), allocatable :: lh, evp, sh, lwd, lwu, swn, lwn, lsnow, ebal
!-----------PET/CWD/MCWD ----------------------------------------------------
    real(wp), dimension(:,:), allocatable :: pet     !! potential ET [mm/day]
    real(wp), dimension(:,:), allocatable :: CWD     !! water deficit [mm]
    real(wp), dimension(:,:), allocatable :: MCWD    !! max deficit [mm]
    real(wp), dimension(:,:), allocatable :: MCWD_yearly !! annual maximum CWD [mm]
!----------------------------------------------------------------------------
    real(wp), dimension(:,:), allocatable :: lh_l, sh_l, swn_l, lwn_l, ebal_l
    real(wp), dimension(:,:), allocatable :: lh_o, sh_o, swn_o, lwn_o, ebal_o
    real(wp), dimension(:,:), allocatable :: tskin, t2m, q2m
    real(wp), dimension(:,:), allocatable :: runoff, runoff_veg, runoff_ice
    real(wp), dimension(:,:), allocatable :: calving, calving_veg, calving_ice
    real(wp), dimension(:,:), allocatable :: bmelt, bmelt_grd, bmelt_flt
    real(wp), dimension(:), allocatable :: fw_atl
    real(wp), dimension(:), allocatable :: fw_natl
    real(wp), dimension(:,:), allocatable :: alb, alb_dif, alb_dir
  end type
 
  type(ts_out), allocatable :: ann_ts(:)
  type(ts_out) :: mon_ts(nmon_year)
  type(sur_out) :: mon_sur(nmon_year), ann_sur

  private
  public :: cmn_diag_init, cmn_diag


contains 
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c m n _ d i a g _ i n i t
  ! Purpose  :  Initialize netcdf output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine cmn_diag_init

    implicit none

    integer :: k


    nout = 0

    ! time series output
    call ts_nc(trim(out_dir)//"/cmn_ts.nc")
    call surf_nc(trim(out_dir)//"/cmn.nc")

    ! allocate
    allocate(ann_ts(ny_out_ts))

    do k=1,nmon_year
     allocate(mon_sur(k)%prc(ni,nj))
     allocate(mon_sur(k)%rain(ni,nj))
     allocate(mon_sur(k)%snow(ni,nj))
     allocate(mon_sur(k)%wind(ni,nj))
     allocate(mon_sur(k)%lh(ni,nj))
     allocate(mon_sur(k)%evp(ni,nj))
     allocate(mon_sur(k)%sh(ni,nj))
!-----------PET/CWD/MCWD ----------------------------------------------------
     allocate(mon_sur(k)%pet(ni,nj))
     allocate(mon_sur(k)%CWD(ni,nj))
     allocate(mon_sur(k)%MCWD(ni,nj))
     allocate(mon_sur(k)%MCWD_yearly(ni,nj))
!----------------------------------------------------------------------------     
     allocate(mon_sur(k)%lwd(ni,nj))
     allocate(mon_sur(k)%lwu(ni,nj))
     allocate(mon_sur(k)%swn(ni,nj))
     allocate(mon_sur(k)%lwn(ni,nj))
     allocate(mon_sur(k)%lsnow(ni,nj))
     allocate(mon_sur(k)%ebal(ni,nj))
     allocate(mon_sur(k)%lh_l(ni,nj))
     allocate(mon_sur(k)%sh_l(ni,nj))
     allocate(mon_sur(k)%swn_l(ni,nj))
     allocate(mon_sur(k)%lwn_l(ni,nj))
     allocate(mon_sur(k)%ebal_l(ni,nj))
     allocate(mon_sur(k)%lh_o(ni,nj))
     allocate(mon_sur(k)%sh_o(ni,nj))
     allocate(mon_sur(k)%swn_o(ni,nj))
     allocate(mon_sur(k)%lwn_o(ni,nj))
     allocate(mon_sur(k)%ebal_o(ni,nj))
     allocate(mon_sur(k)%tskin(ni,nj))
     allocate(mon_sur(k)%t2m(ni,nj))
     allocate(mon_sur(k)%q2m(ni,nj))
     allocate(mon_sur(k)%runoff(ni,nj))
     allocate(mon_sur(k)%runoff_veg(ni,nj))
     allocate(mon_sur(k)%runoff_ice(ni,nj))
     allocate(mon_sur(k)%calving(ni,nj))
     allocate(mon_sur(k)%calving_veg(ni,nj))
     allocate(mon_sur(k)%calving_ice(ni,nj))
     allocate(mon_sur(k)%bmelt(ni,nj))
     allocate(mon_sur(k)%bmelt_grd(ni,nj))
     allocate(mon_sur(k)%bmelt_flt(ni,nj))
     allocate(mon_sur(k)%alb(ni,nj))
     allocate(mon_sur(k)%alb_dif(ni,nj))
     allocate(mon_sur(k)%alb_dir(ni,nj))
    enddo

     allocate(ann_sur%mask_ice(ni,nj))
     allocate(ann_sur%mask_smb(ni,nj))
     allocate(ann_sur%f_stp(ni,nj,nsurf))
     allocate(ann_sur%f_ocn(ni,nj))
     allocate(ann_sur%f_lnd(ni,nj))
     allocate(ann_sur%f_ice(ni,nj))
     allocate(ann_sur%prc(ni,nj))
     allocate(ann_sur%rain(ni,nj))
     allocate(ann_sur%snow(ni,nj))
     allocate(ann_sur%wind(ni,nj))
     allocate(ann_sur%lh(ni,nj))
     allocate(ann_sur%evp(ni,nj))
     allocate(ann_sur%sh(ni,nj))
!-----------PET/CWD/MCWD ----------------------------------------------------
     allocate(ann_sur%pet(ni,nj))
     allocate(ann_sur%CWD(ni,nj))
     allocate(ann_sur%MCWD(ni,nj))
     allocate(ann_sur%MCWD_yearly(ni,nj))
!----------------------------------------------------------------------------
     allocate(ann_sur%lwd(ni,nj))
     allocate(ann_sur%lwu(ni,nj))
     allocate(ann_sur%swn(ni,nj))
     allocate(ann_sur%lwn(ni,nj))
     allocate(ann_sur%lsnow(ni,nj))
     allocate(ann_sur%ebal(ni,nj))
     allocate(ann_sur%lh_l(ni,nj))
     allocate(ann_sur%sh_l(ni,nj))
     allocate(ann_sur%swn_l(ni,nj))
     allocate(ann_sur%lwn_l(ni,nj))
     allocate(ann_sur%ebal_l(ni,nj))
     allocate(ann_sur%lh_o(ni,nj))
     allocate(ann_sur%sh_o(ni,nj))
     allocate(ann_sur%swn_o(ni,nj))
     allocate(ann_sur%lwn_o(ni,nj))
     allocate(ann_sur%ebal_O(ni,nj))
     allocate(ann_sur%tskin(ni,nj))
     allocate(ann_sur%t2m(ni,nj))
     allocate(ann_sur%q2m(ni,nj))
     allocate(ann_sur%runoff(ni,nj))
     allocate(ann_sur%runoff_veg(ni,nj))
     allocate(ann_sur%runoff_ice(ni,nj))
     allocate(ann_sur%calving(ni,nj))
     allocate(ann_sur%calving_veg(ni,nj))
     allocate(ann_sur%calving_ice(ni,nj))
     allocate(ann_sur%bmelt(ni,nj))
     allocate(ann_sur%bmelt_grd(ni,nj))
     allocate(ann_sur%bmelt_flt(ni,nj))
     allocate(ann_sur%fw_atl(nj))
     allocate(ann_sur%fw_natl(nj))
     allocate(ann_sur%alb(ni,nj))
     allocate(ann_sur%alb_dif(ni,nj))
     allocate(ann_sur%alb_dir(ni,nj))


     return

  end subroutine cmn_diag_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c m n _ d i a g
  ! Purpose  :  compute diagnostis
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine cmn_diag(cmn)

    use coupler, only : cmn_class

    implicit none

    type(cmn_class), intent(in) :: cmn

    integer :: i, j, k, m, y
    real(wp) :: tmp
    real(wp) :: mon_avg


    ! current index
    y = y_out_ts_clim

    mon_avg = 1._wp/nstep_mon

    if( time_soy ) then
      do k=1,nmon_year
        mon_ts(k)%t2m = 0._wp
        mon_ts(k)%lh = 0._wp
        mon_ts(k)%sh = 0._wp
        mon_ts(k)%lwd = 0._wp
        mon_ts(k)%lwu = 0._wp
        mon_ts(k)%swn = 0._wp
        mon_ts(k)%lwn = 0._wp
        mon_ts(k)%lsnow = 0._wp
        mon_ts(k)%lh_l = 0._wp
        mon_ts(k)%sh_l = 0._wp
        mon_ts(k)%lwd_l = 0._wp
        mon_ts(k)%lwu_l = 0._wp
        mon_ts(k)%swn_l = 0._wp
        mon_ts(k)%lwn_l = 0._wp
        mon_ts(k)%lsnow_l = 0._wp
        mon_ts(k)%lh_o = 0._wp
        mon_ts(k)%sh_o = 0._wp
        mon_ts(k)%lwd_o = 0._wp
        mon_ts(k)%lwu_o = 0._wp
        mon_ts(k)%swn_o = 0._wp
        mon_ts(k)%lwn_o = 0._wp
        mon_ts(k)%lsnow_o = 0._wp
        mon_ts(k)%runoff = 0._wp
        mon_ts(k)%runoff_veg = 0._wp
        mon_ts(k)%runoff_ice = 0._wp
        mon_ts(k)%calving = 0._wp
        mon_ts(k)%calving_veg = 0._wp
        mon_ts(k)%calving_ice = 0._wp
        mon_ts(k)%bmelt = 0._wp
        mon_ts(k)%bmelt_grd = 0._wp
        mon_ts(k)%bmelt_flt = 0._wp
        mon_ts(k)%evp = 0._wp
        mon_ts(k)%prc = 0._wp
        mon_ts(k)%prc_atl = 0._wp
        mon_ts(k)%prc_pac = 0._wp
        mon_ts(k)%prc_ind = 0._wp
        mon_ts(k)%prc_so  = 0._wp
        mon_ts(k)%evp_atl = 0._wp
        mon_ts(k)%evp_pac = 0._wp
        mon_ts(k)%evp_ind = 0._wp
        mon_ts(k)%evp_so  = 0._wp
        mon_ts(k)%runoff_atl = 0._wp
        mon_ts(k)%runoff_pac = 0._wp
        mon_ts(k)%runoff_ind = 0._wp
        mon_ts(k)%runoff_so  = 0._wp
        mon_ts(k)%calving_atl = 0._wp
        mon_ts(k)%calving_pac = 0._wp
        mon_ts(k)%calving_ind = 0._wp
        mon_ts(k)%calving_so  = 0._wp
        mon_ts(k)%bmelt_atl = 0._wp
        mon_ts(k)%bmelt_pac = 0._wp
        mon_ts(k)%bmelt_ind = 0._wp
        mon_ts(k)%bmelt_so  = 0._wp
        mon_ts(k)%fw_atl = 0._wp
        mon_ts(k)%fw_natl = 0._wp
        mon_ts(k)%fw_pac = 0._wp
        mon_ts(k)%fw_ind = 0._wp
        mon_ts(k)%fw_so  = 0._wp
      enddo
    endif


    ! globally integrated values
    mon_ts(mon)%runoff = mon_ts(mon)%runoff + sum(cmn%runoff_o * cmn%f_ocn*area) ! kg/s
    mon_ts(mon)%runoff_veg = mon_ts(mon)%runoff_veg + sum(cmn%runoff_veg_o * cmn%f_ocn*area) ! kg/s
    mon_ts(mon)%runoff_ice = mon_ts(mon)%runoff_ice + sum(cmn%runoff_ice_o * cmn%f_ocn*area) ! kg/s
    mon_ts(mon)%calving = mon_ts(mon)%calving + sum(cmn%calving_o * cmn%f_ocn*area) ! kg/s
    mon_ts(mon)%calving_veg = mon_ts(mon)%calving_veg + sum(cmn%calving_veg_o * cmn%f_ocn*area) ! kg/s
    mon_ts(mon)%calving_ice = mon_ts(mon)%calving_ice + sum(cmn%calving_ice_o * cmn%f_ocn*area) ! kg/s
    mon_ts(mon)%bmelt = mon_ts(mon)%bmelt + sum(cmn%bmelt_o * cmn%f_ocn*area) ! kg/s
    mon_ts(mon)%bmelt_grd = mon_ts(mon)%bmelt_grd + sum(cmn%bmelt_grd_o * cmn%f_ocn*area) ! kg/s
    mon_ts(mon)%bmelt_flt = mon_ts(mon)%bmelt_flt + sum(cmn%bmelt_flt_o * cmn%f_ocn*area) ! kg/s
    mon_ts(mon)%evp = mon_ts(mon)%evp + sum(sum(cmn%evp*cmn%f_stp,3) * area) ! kg/s
    mon_ts(mon)%prc = mon_ts(mon)%prc + sum(sum((cmn%rain+cmn%snow)*cmn%f_stp,3) * area) ! kg/s

    mon_ts(mon)%prc_atl = mon_ts(mon)%prc_atl + sum(cmn%prc * cmn%f_ocn*area, basin_mask.eq.i_atlantic) ! kg/s
    mon_ts(mon)%prc_pac = mon_ts(mon)%prc_pac + sum(cmn%prc * cmn%f_ocn*area, basin_mask.eq.i_pacific) ! kg/s
    mon_ts(mon)%prc_ind = mon_ts(mon)%prc_ind + sum(cmn%prc * cmn%f_ocn*area, basin_mask.eq.i_indian) ! kg/s
    mon_ts(mon)%prc_so  = mon_ts(mon)%prc_so  + sum(cmn%prc * cmn%f_ocn*area, basin_mask.eq.i_southern) ! kg/s

    mon_ts(mon)%evp_atl = mon_ts(mon)%evp_atl + sum(((1.-cmn%f_sic)*cmn%evp(:,:,1)+cmn%f_sic*cmn%evp(:,:,2)) * cmn%f_ocn*area, basin_mask.eq.i_atlantic) ! kg/s
    mon_ts(mon)%evp_pac = mon_ts(mon)%evp_pac + sum(((1.-cmn%f_sic)*cmn%evp(:,:,1)+cmn%f_sic*cmn%evp(:,:,2)) * cmn%f_ocn*area, basin_mask.eq.i_pacific) ! kg/s
    mon_ts(mon)%evp_ind = mon_ts(mon)%evp_ind + sum(((1.-cmn%f_sic)*cmn%evp(:,:,1)+cmn%f_sic*cmn%evp(:,:,2)) * cmn%f_ocn*area, basin_mask.eq.i_indian) ! kg/s
    mon_ts(mon)%evp_so  = mon_ts(mon)%evp_so  + sum(((1.-cmn%f_sic)*cmn%evp(:,:,1)+cmn%f_sic*cmn%evp(:,:,2)) * cmn%f_ocn*area, basin_mask.eq.i_southern) ! kg/s

    mon_ts(mon)%runoff_atl = mon_ts(mon)%runoff_atl + sum(cmn%runoff_o * cmn%f_ocn*area, basin_mask.eq.i_atlantic) ! kg/s
    mon_ts(mon)%runoff_pac = mon_ts(mon)%runoff_pac + sum(cmn%runoff_o * cmn%f_ocn*area, basin_mask.eq.i_pacific) ! kg/s
    mon_ts(mon)%runoff_ind = mon_ts(mon)%runoff_ind + sum(cmn%runoff_o * cmn%f_ocn*area, basin_mask.eq.i_indian) ! kg/s
    mon_ts(mon)%runoff_so  = mon_ts(mon)%runoff_so  + sum(cmn%runoff_o * cmn%f_ocn*area, basin_mask.eq.i_southern) ! kg/s

    mon_ts(mon)%calving_atl = mon_ts(mon)%calving_atl + sum(cmn%calving_o * cmn%f_ocn*area, basin_mask.eq.i_atlantic) ! kg/s
    mon_ts(mon)%calving_pac = mon_ts(mon)%calving_pac + sum(cmn%calving_o * cmn%f_ocn*area, basin_mask.eq.i_pacific) ! kg/s
    mon_ts(mon)%calving_ind = mon_ts(mon)%calving_ind + sum(cmn%calving_o * cmn%f_ocn*area, basin_mask.eq.i_indian) ! kg/s
    mon_ts(mon)%calving_so  = mon_ts(mon)%calving_so  + sum(cmn%calving_o * cmn%f_ocn*area, basin_mask.eq.i_southern) ! kg/s

    mon_ts(mon)%bmelt_atl = mon_ts(mon)%bmelt_atl + sum(cmn%bmelt_o * cmn%f_ocn*area, basin_mask.eq.i_atlantic) ! kg/s
    mon_ts(mon)%bmelt_pac = mon_ts(mon)%bmelt_pac + sum(cmn%bmelt_o * cmn%f_ocn*area, basin_mask.eq.i_pacific) ! kg/s
    mon_ts(mon)%bmelt_ind = mon_ts(mon)%bmelt_ind + sum(cmn%bmelt_o * cmn%f_ocn*area, basin_mask.eq.i_indian) ! kg/s
    mon_ts(mon)%bmelt_so  = mon_ts(mon)%bmelt_so  + sum(cmn%bmelt_o * cmn%f_ocn*area, basin_mask.eq.i_southern) ! kg/s

    mon_ts(mon)%fw_atl = mon_ts(mon)%fw_atl &
      + sum((cmn%prc-((1.-cmn%f_sic)*cmn%evp(:,:,1)+cmn%f_sic*cmn%evp(:,:,2))+cmn%runoff_o+cmn%calving_o+cmn%bmelt_o) &
      * cmn%f_ocn*area, basin_mask.eq.i_atlantic) ! kg/s
    mon_ts(mon)%fw_pac = mon_ts(mon)%fw_pac &
      + sum((cmn%prc-((1.-cmn%f_sic)*cmn%evp(:,:,1)+cmn%f_sic*cmn%evp(:,:,2))+cmn%runoff_o+cmn%calving_o+cmn%bmelt_o) &
      * cmn%f_ocn*area, basin_mask.eq.i_pacific) ! kg/s
    mon_ts(mon)%fw_ind = mon_ts(mon)%fw_ind &
      + sum((cmn%prc-((1.-cmn%f_sic)*cmn%evp(:,:,1)+cmn%f_sic*cmn%evp(:,:,2))+cmn%runoff_o+cmn%calving_o+cmn%bmelt_o) &
      * cmn%f_ocn*area, basin_mask.eq.i_indian) ! kg/s
    mon_ts(mon)%fw_so  = mon_ts(mon)%fw_so  &
      + sum((cmn%prc-((1.-cmn%f_sic)*cmn%evp(:,:,1)+cmn%f_sic*cmn%evp(:,:,2))+cmn%runoff_o+cmn%calving_o+cmn%bmelt_o) &
      * cmn%f_ocn*area, basin_mask.eq.i_southern) ! kg/s
    do i=1,ni
      mon_ts(mon)%fw_natl= mon_ts(mon)%fw_natl &
        + sum((cmn%prc(i,:)-((1.-cmn%f_sic(i,:))*cmn%evp(i,:,1)+cmn%f_sic(i,:)*cmn%evp(i,:,2))+cmn%runoff_o(i,:)+cmn%calving_o(i,:)+cmn%bmelt_o(i,:)) &
        * cmn%f_ocn(i,:)*area(i,:), basin_mask(i,:).eq.i_atlantic .and. lat>0.) ! kg/s
    enddo

    ! global mean values
    mon_ts(mon)%t2m = mon_ts(mon)%t2m + sum(sum(cmn%t2m*cmn%f_stp,3) * area) / sum(area) ! K
    mon_ts(mon)%sh  = mon_ts(mon)%sh  + sum(sum(cmn%sh*cmn%f_stp,3) * area) / sum(area) ! W/m2
    mon_ts(mon)%lh  = mon_ts(mon)%lh  + sum(sum(cmn%lh*cmn%f_stp,3) * area) / sum(area) ! W/m2
    mon_ts(mon)%swn = mon_ts(mon)%swn + sum(sum(cmn%swnet*cmn%f_stp,3) * area) / sum(area) ! W/m2
    mon_ts(mon)%lwd = mon_ts(mon)%lwd + sum(sum(cmn%lwd*cmn%f_stp,3) * area) / sum(area) ! W/m2
    mon_ts(mon)%lwu = mon_ts(mon)%lwu + sum(sum(cmn%lwu*cmn%f_stp,3) * area) / sum(area) ! W/m2
    mon_ts(mon)%lwn = mon_ts(mon)%lwn + sum((sum(cmn%lwd*cmn%f_stp,3)-sum(cmn%lwu*cmn%f_stp,3)) * area) / sum(area) ! W/m2
    mon_ts(mon)%lsnow= mon_ts(mon)%lsnow + sum(sum(cmn%snow*cmn%f_stp,3) * area) / sum(area) * Lf ! W/m2
    if (sum(cmn%f_ocn).gt.0._wp) then
      mon_ts(mon)%sh_o  = mon_ts(mon)%sh_o  + sum(sum(cmn%sh(:,:,1:2)*cmn%f_stp(:,:,1:2),3)    * area) / sum(cmn%f_ocn*area) ! W/m2
      mon_ts(mon)%lh_o  = mon_ts(mon)%lh_o  + sum(sum(cmn%lh(:,:,1:2)*cmn%f_stp(:,:,1:2),3)    * area) / sum(cmn%f_ocn*area) ! W/m2
      mon_ts(mon)%swn_o = mon_ts(mon)%swn_o + sum(sum(cmn%swnet(:,:,1:2)*cmn%f_stp(:,:,1:2),3) * area) / sum(cmn%f_ocn*area) ! W/m2
      mon_ts(mon)%lwd_o = mon_ts(mon)%lwd_o + sum(sum(cmn%lwd(:,:,1:2)*cmn%f_stp(:,:,1:2),3)   * area) / sum(cmn%f_ocn*area) ! W/m2
      mon_ts(mon)%lwu_o = mon_ts(mon)%lwu_o + sum(sum(cmn%lwu(:,:,1:2)*cmn%f_stp(:,:,1:2),3)   * area) / sum(cmn%f_ocn*area) ! W/m2
      mon_ts(mon)%lwn_o = mon_ts(mon)%lwn_o + sum(sum((cmn%lwd(:,:,1:2)-cmn%lwu(:,:,1:2))*cmn%f_stp(:,:,1:2),3)   * area) / sum(cmn%f_ocn*area) ! W/m2
      mon_ts(mon)%lsnow_o= mon_ts(mon)%lsnow_o + sum(sum(cmn%snow(:,:,1:2)*cmn%f_stp(:,:,1:2),3) * area) / sum(cmn%f_ocn*area) * Lf ! W/m2
    endif
    if (sum(1._wp-cmn%f_ocn).gt.0._wp) then
      mon_ts(mon)%sh_l  = mon_ts(mon)%sh_l  + sum(sum(cmn%sh(:,:,3:nsurf)*cmn%f_stp(:,:,3:nsurf),3)    * area) / sum((1._wp-cmn%f_ocn)*area) ! W/m2
      mon_ts(mon)%lh_l  = mon_ts(mon)%lh_l  + sum(sum(cmn%lh(:,:,3:nsurf)*cmn%f_stp(:,:,3:nsurf),3)    * area) / sum((1._wp-cmn%f_ocn)*area) ! W/m2
      mon_ts(mon)%swn_l = mon_ts(mon)%swn_l + sum(sum(cmn%swnet(:,:,3:nsurf)*cmn%f_stp(:,:,3:nsurf),3) * area) / sum((1._wp-cmn%f_ocn)*area) ! W/m2
      mon_ts(mon)%lwd_l = mon_ts(mon)%lwd_l + sum(sum(cmn%lwd(:,:,3:nsurf)*cmn%f_stp(:,:,3:nsurf),3)   * area) / sum((1._wp-cmn%f_ocn)*area) ! W/m2
      mon_ts(mon)%lwu_l = mon_ts(mon)%lwu_l + sum(sum(cmn%lwu(:,:,3:nsurf)*cmn%f_stp(:,:,3:nsurf),3)   * area) / sum((1._wp-cmn%f_ocn)*area) ! W/m2
      mon_ts(mon)%lwn_l = mon_ts(mon)%lwn_l + sum(sum((cmn%lwd(:,:,3:nsurf)-cmn%lwu(:,:,3:nsurf))*cmn%f_stp(:,:,3:nsurf),3)   * area) / sum((1._wp-cmn%f_ocn)*area) ! W/m2
      mon_ts(mon)%lsnow_l= mon_ts(mon)%lsnow_l + sum(sum(cmn%snow(:,:,3:nsurf)*cmn%f_stp(:,:,3:nsurf),3) * area) / sum((1._wp-cmn%f_ocn)*area) * Lf ! W/m2
    endif

    if (time_eom) then
      mon_ts(mon)%t2m = mon_ts(mon)%t2m * mon_avg - T0
      mon_ts(mon)%sh = mon_ts(mon)%sh * mon_avg 
      mon_ts(mon)%lh = mon_ts(mon)%lh * mon_avg 
      mon_ts(mon)%lwd = mon_ts(mon)%lwd * mon_avg 
      mon_ts(mon)%lwu = mon_ts(mon)%lwu * mon_avg 
      mon_ts(mon)%swn = mon_ts(mon)%swn * mon_avg
      mon_ts(mon)%lwn = mon_ts(mon)%lwn * mon_avg
      mon_ts(mon)%lsnow = mon_ts(mon)%lsnow * mon_avg
      mon_ts(mon)%ebal= mon_ts(mon)%swn + mon_ts(mon)%lwn - mon_ts(mon)%sh - mon_ts(mon)%lh - mon_ts(mon)%lsnow
      mon_ts(mon)%sh_o = mon_ts(mon)%sh_o * mon_avg
      mon_ts(mon)%lh_o = mon_ts(mon)%lh_o * mon_avg 
      mon_ts(mon)%lwd_o = mon_ts(mon)%lwd_o * mon_avg 
      mon_ts(mon)%lwu_o = mon_ts(mon)%lwu_o * mon_avg 
      mon_ts(mon)%swn_o = mon_ts(mon)%swn_o * mon_avg
      mon_ts(mon)%lwn_o = mon_ts(mon)%lwn_o * mon_avg 
      mon_ts(mon)%lsnow_o = mon_ts(mon)%lsnow_o * mon_avg
      mon_ts(mon)%ebal_o= mon_ts(mon)%swn_o + mon_ts(mon)%lwn_o - mon_ts(mon)%sh_o - mon_ts(mon)%lh_o - mon_ts(mon)%lsnow_o
      mon_ts(mon)%sh_l = mon_ts(mon)%sh_l * mon_avg 
      mon_ts(mon)%lh_l = mon_ts(mon)%lh_l * mon_avg 
      mon_ts(mon)%lwd_l = mon_ts(mon)%lwd_l * mon_avg 
      mon_ts(mon)%lwu_l = mon_ts(mon)%lwu_l * mon_avg 
      mon_ts(mon)%swn_l = mon_ts(mon)%swn_l * mon_avg
      mon_ts(mon)%lwn_l = mon_ts(mon)%lwn_l * mon_avg
      mon_ts(mon)%lsnow_l = mon_ts(mon)%lsnow_l * mon_avg
      mon_ts(mon)%ebal_l= mon_ts(mon)%swn_l + mon_ts(mon)%lwn_l - mon_ts(mon)%sh_l - mon_ts(mon)%lh_l - mon_ts(mon)%lsnow_l
      mon_ts(mon)%runoff = mon_ts(mon)%runoff * 1.d-15 * sec_mon * mon_avg ! 10^15 kg/mon
      mon_ts(mon)%runoff_veg = mon_ts(mon)%runoff_veg * 1.d-15 * sec_mon * mon_avg ! 10^15 kg/mon
      mon_ts(mon)%runoff_ice = mon_ts(mon)%runoff_ice * 1.d-15 * sec_mon * mon_avg ! 10^15 kg/mon
      mon_ts(mon)%runoff_atl = mon_ts(mon)%runoff_atl * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%runoff_pac = mon_ts(mon)%runoff_pac * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%runoff_ind = mon_ts(mon)%runoff_ind * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%runoff_so  = mon_ts(mon)%runoff_so  * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%calving = mon_ts(mon)%calving * 1.d-15 * sec_mon * mon_avg ! 10^15 kg/mon
      mon_ts(mon)%calving_veg = mon_ts(mon)%calving_veg * 1.d-15 * sec_mon * mon_avg ! 10^15 kg/mon
      mon_ts(mon)%calving_ice = mon_ts(mon)%calving_ice * 1.d-15 * sec_mon * mon_avg ! 10^15 kg/mon
      mon_ts(mon)%calving_atl = mon_ts(mon)%calving_atl * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%calving_pac = mon_ts(mon)%calving_pac * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%calving_ind = mon_ts(mon)%calving_ind * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%calving_so  = mon_ts(mon)%calving_so  * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%bmelt = mon_ts(mon)%calving * 1.d-15 * sec_mon * mon_avg ! 10^15 kg/mon
      mon_ts(mon)%bmelt_grd = mon_ts(mon)%bmelt_grd * 1.d-15 * sec_mon * mon_avg ! 10^15 kg/mon
      mon_ts(mon)%bmelt_flt = mon_ts(mon)%bmelt_flt * 1.d-15 * sec_mon * mon_avg ! 10^15 kg/mon
      mon_ts(mon)%bmelt_atl = mon_ts(mon)%bmelt_atl * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%bmelt_pac = mon_ts(mon)%bmelt_pac * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%bmelt_ind = mon_ts(mon)%bmelt_ind * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%bmelt_so  = mon_ts(mon)%bmelt_so  * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%evp = mon_ts(mon)%evp * 1.d-15 * sec_mon * mon_avg ! 10^15 kg/mon
      mon_ts(mon)%evp_atl = mon_ts(mon)%evp_atl * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%evp_pac = mon_ts(mon)%evp_pac * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%evp_ind = mon_ts(mon)%evp_ind * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%evp_so  = mon_ts(mon)%evp_so  * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%prc = mon_ts(mon)%prc * 1.d-15 * sec_mon * mon_avg ! 10^15 kg/mon
      mon_ts(mon)%prc_atl = mon_ts(mon)%prc_atl * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%prc_pac = mon_ts(mon)%prc_pac * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%prc_ind = mon_ts(mon)%prc_ind * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%prc_so  = mon_ts(mon)%prc_so  * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%fw_atl = mon_ts(mon)%fw_atl * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%fw_natl= mon_ts(mon)%fw_natl * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%fw_pac = mon_ts(mon)%fw_pac * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%fw_ind = mon_ts(mon)%fw_ind * 1.d-6 * 1d-3 * mon_avg ! Sv
      mon_ts(mon)%fw_so  = mon_ts(mon)%fw_so  * 1.d-6 * 1d-3 * mon_avg ! Sv
    endif

    if( time_eoy ) then

     call ts_ave(mon_ts,ann_ts(y))

     ! write to netcdf file 
     if (time_out_ts_clim) then
       call ts_nc_write(trim(out_dir)//"/cmn_ts.nc",ann_ts(1:y),year_clim-y+1,y)
     endif

     ! write to standard output
     if (mod(year,10).eq.1) then
       print '(a7,a9,13a7)','cmn','year','CO2','CO2r','t2m','prc','evp','run','sh','lh','lwn','swn','ebal','ebal_o','ebal_l'
     endif

     print '(a7,i9,10F7.1,3F7.2)', &
     'cmn',year_now,cmn%co2,cmn%co2_rad,ann_ts(y)%t2m,ann_ts(y)%prc,ann_ts(y)%evp,ann_ts(y)%runoff,ann_ts(y)%sh,ann_ts(y)%lh,ann_ts(y)%lwn,ann_ts(y)%swn,ann_ts(y)%ebal,ann_ts(y)%ebal_o,ann_ts(y)%ebal_l
    endif


    ! 2D-3D output 
    if( time_out_cmn ) then

      if( time_soy ) then
        do m=1,nmon_year
          mon_sur(m)%prc      = 0._wp
          mon_sur(m)%rain     = 0._wp
          mon_sur(m)%snow     = 0._wp
          mon_sur(m)%evp      = 0._wp
!-----------PET/CWD/MCWD ----------------------------------------------------
          mon_sur(m)%pet       = 0._wp
          mon_sur(m)%CWD       = 0._wp
          mon_sur(m)%MCWD      = 0._wp
!--------------------------------------------------------------------------
          mon_sur(m)%runoff   = 0._wp
          mon_sur(m)%runoff_veg   = 0._wp
          mon_sur(m)%runoff_ice   = 0._wp
          mon_sur(m)%calving  = 0._wp
          mon_sur(m)%calving_veg  = 0._wp
          mon_sur(m)%calving_ice  = 0._wp
          mon_sur(m)%bmelt    = 0._wp
          mon_sur(m)%bmelt_grd= 0._wp
          mon_sur(m)%bmelt_flt= 0._wp
          mon_sur(m)%wind     = 0._wp
          mon_sur(m)%lh       = 0._wp
          mon_sur(m)%sh       = 0._wp
          mon_sur(m)%lwd      = 0._wp
          mon_sur(m)%lwu      = 0._wp
          mon_sur(m)%swn      = 0._wp
          mon_sur(m)%lwn      = 0._wp
          mon_sur(m)%lsnow    = 0._wp
          mon_sur(m)%lh_l     = 0._wp
          mon_sur(m)%sh_l     = 0._wp
          mon_sur(m)%swn_l    = 0._wp
          mon_sur(m)%lwn_l    = 0._wp
          mon_sur(m)%lh_o     = 0._wp
          mon_sur(m)%sh_o     = 0._wp
          mon_sur(m)%swn_o    = 0._wp
          mon_sur(m)%lwn_o    = 0._wp
          mon_sur(m)%tskin    = 0._wp
          mon_sur(m)%t2m      = 0._wp
          mon_sur(m)%q2m      = 0._wp
          mon_sur(m)%alb      = 0._wp
          mon_sur(m)%alb_dir  = 0._wp
          mon_sur(m)%alb_dif  = 0._wp
        enddo
      endif

      mon_sur(mon)%rain = mon_sur(mon)%rain + sum(cmn%rain*cmn%f_stp,3) ! kg/m2/s
      mon_sur(mon)%snow = mon_sur(mon)%snow + sum(cmn%snow*cmn%f_stp,3) ! kg/m2/s
      mon_sur(mon)%prc = mon_sur(mon)%rain + mon_sur(mon)%snow ! kg/m2/s
      mon_sur(mon)%evp = mon_sur(mon)%evp + sum(cmn%evp*cmn%f_stp,3) ! kg/m2/s
!-----------PET/CWD/MCWD ----------------------------------------------------
      ! PET: accumulate for monthly average
      mon_sur(mon)%pet = mon_sur(mon)%pet + sum(cmn%pet*cmn%f_stp,3) ! mm/day
!----------------------------------------------------------------------------
      mon_sur(mon)%runoff = mon_sur(mon)%runoff + cmn%runoff_o ! kg/m2/s
      mon_sur(mon)%runoff_veg = mon_sur(mon)%runoff_veg + cmn%runoff_veg_o ! kg/m2/s
      mon_sur(mon)%runoff_ice = mon_sur(mon)%runoff_ice + cmn%runoff_ice_o ! kg/m2/s
      mon_sur(mon)%calving = mon_sur(mon)%calving + cmn%calving_o ! kg/m2/s
      mon_sur(mon)%calving_veg = mon_sur(mon)%calving_veg + cmn%calving_veg_o ! kg/m2/s
      mon_sur(mon)%calving_ice = mon_sur(mon)%calving_ice + cmn%calving_ice_o ! kg/m2/s
      mon_sur(mon)%bmelt = mon_sur(mon)%bmelt + cmn%bmelt_o ! kg/m2/s
      mon_sur(mon)%bmelt_grd = mon_sur(mon)%bmelt_grd + cmn%bmelt_grd_o ! kg/m2/s
      mon_sur(mon)%bmelt_flt = mon_sur(mon)%bmelt_flt + cmn%bmelt_flt_o ! kg/m2/s
      mon_sur(mon)%wind = mon_sur(mon)%wind + sum(cmn%wind*cmn%f_stp,3) ! m/s
      mon_sur(mon)%lh = mon_sur(mon)%lh + sum(cmn%lh*cmn%f_stp,3) ! W/m2
      mon_sur(mon)%sh = mon_sur(mon)%sh + sum(cmn%sh*cmn%f_stp,3) ! W/m2
      mon_sur(mon)%lwd = mon_sur(mon)%lwd + sum(cmn%lwd*cmn%f_stp,3) ! W/m2
      mon_sur(mon)%lwu = mon_sur(mon)%lwu + sum(cmn%lwu*cmn%f_stp,3) ! W/m2
      mon_sur(mon)%swn = mon_sur(mon)%swn + sum(cmn%swnet*cmn%f_stp,3) ! W/m2
      mon_sur(mon)%lwn = mon_sur(mon)%lwn + (sum(cmn%lwd*cmn%f_stp,3)-sum(cmn%lwu*cmn%f_stp,3)) ! W/m2
      mon_sur(mon)%lsnow = mon_sur(mon)%lsnow + sum(cmn%snow*cmn%f_stp,3) * Lf ! W/m2
      mon_sur(mon)%lh_l = mon_sur(mon)%lh_l + sum(cmn%lh(:,:,3:nsurf)*cmn%f_stp(:,:,3:nsurf),3) ! W/m2
      mon_sur(mon)%sh_l = mon_sur(mon)%sh_l + sum(cmn%sh(:,:,3:nsurf)*cmn%f_stp(:,:,3:nsurf),3) ! W/m2
      mon_sur(mon)%swn_l = mon_sur(mon)%swn_l + sum(cmn%swnet(:,:,3:nsurf)*cmn%f_stp(:,:,3:nsurf),3) ! W/m2
      mon_sur(mon)%lwn_l = mon_sur(mon)%lwn_l + (sum(cmn%lwd(:,:,3:nsurf)*cmn%f_stp(:,:,3:nsurf),3)-sum(cmn%lwu(:,:,3:nsurf)*cmn%f_stp(:,:,3:nsurf),3)) ! W/m2
      mon_sur(mon)%lh_o = mon_sur(mon)%lh_o + sum(cmn%lh(:,:,1:2)*cmn%f_stp(:,:,1:2),3) ! W/m2
      mon_sur(mon)%sh_o = mon_sur(mon)%sh_o + sum(cmn%sh(:,:,1:2)*cmn%f_stp(:,:,1:2),3) ! W/m2
      mon_sur(mon)%swn_o = mon_sur(mon)%swn_o + sum(cmn%swnet(:,:,1:2)*cmn%f_stp(:,:,1:2),3) ! W/m2
      mon_sur(mon)%lwn_o = mon_sur(mon)%lwn_o + (sum(cmn%lwd(:,:,1:2)*cmn%f_stp(:,:,1:2),3)-sum(cmn%lwu(:,:,1:2)*cmn%f_stp(:,:,1:2),3)) ! W/m2
      mon_sur(mon)%tskin = mon_sur(mon)%tskin + sum(cmn%t_skin*cmn%f_stp,3) ! K
      mon_sur(mon)%t2m   = mon_sur(mon)%t2m + sum(cmn%t2m*cmn%f_stp,3) ! K
      mon_sur(mon)%q2m   = mon_sur(mon)%q2m + sum(cmn%q2m*cmn%f_stp,3) ! kg/kg
      mon_sur(mon)%alb_dir = mon_sur(mon)%alb_dir + (0.6_wp*sum(cmn%alb_vis_dir*cmn%f_stp,3) + 0.4_wp*sum(cmn%alb_nir_dir*cmn%f_stp,3))
      mon_sur(mon)%alb_dif = mon_sur(mon)%alb_dif + (0.6_wp*sum(cmn%alb_vis_dif*cmn%f_stp,3) + 0.4_wp*sum(cmn%alb_nir_dif*cmn%f_stp,3))

      if (time_eom) then
!-----------PET/CWD/MCWD ----------------------------------------------------
        ! CWD and MCWD: copy from coupler (already monthly values)
        mon_sur(mon)%CWD = cmn%CWD_monthly
        mon_sur(mon)%MCWD = cmn%MCWD_monthly
!----------------------------------------------------------------------------
        mon_sur(mon)%prc = mon_sur(mon)%prc * mon_avg * sec_day ! kg/m2/day
        mon_sur(mon)%rain = mon_sur(mon)%rain * mon_avg * sec_day ! kg/m2/day
        mon_sur(mon)%snow = mon_sur(mon)%snow * mon_avg * sec_day ! kg/m2/day
        mon_sur(mon)%evp = mon_sur(mon)%evp * mon_avg * sec_day ! kg/m2/day
!-----------PET/CWD/MCWD ----------------------------------------------------
        ! Average PET (no sec_day conversion, already in mm/day)
        mon_sur(mon)%pet = mon_sur(mon)%pet * mon_avg ! mm/day
!----------------------------------------------------------------------------
        mon_sur(mon)%runoff = mon_sur(mon)%runoff * mon_avg * sec_day ! kg/m2/day
        mon_sur(mon)%runoff_veg = mon_sur(mon)%runoff_veg * mon_avg * sec_day ! kg/m2/day
        mon_sur(mon)%runoff_ice = mon_sur(mon)%runoff_ice * mon_avg * sec_day ! kg/m2/day
        mon_sur(mon)%calving = mon_sur(mon)%calving * mon_avg * sec_day ! kg/m2/day
        mon_sur(mon)%calving_veg = mon_sur(mon)%calving_veg * mon_avg * sec_day ! kg/m2/day
        mon_sur(mon)%calving_ice = mon_sur(mon)%calving_ice * mon_avg * sec_day ! kg/m2/day
        mon_sur(mon)%bmelt = mon_sur(mon)%bmelt * mon_avg * sec_day ! kg/m2/day
        mon_sur(mon)%bmelt_grd = mon_sur(mon)%bmelt_grd * mon_avg * sec_day ! kg/m2/day
        mon_sur(mon)%bmelt_flt = mon_sur(mon)%bmelt_flt * mon_avg * sec_day ! kg/m2/day
        mon_sur(mon)%wind = mon_sur(mon)%wind * mon_avg
        mon_sur(mon)%lh = mon_sur(mon)%lh * mon_avg
        mon_sur(mon)%sh = mon_sur(mon)%sh * mon_avg
        mon_sur(mon)%lwd = mon_sur(mon)%lwd * mon_avg
        mon_sur(mon)%lwu = mon_sur(mon)%lwu * mon_avg
        mon_sur(mon)%swn = mon_sur(mon)%swn * mon_avg
        mon_sur(mon)%lwn = mon_sur(mon)%lwn * mon_avg
        mon_sur(mon)%lsnow = mon_sur(mon)%lsnow * mon_avg
        mon_sur(mon)%ebal= mon_sur(mon)%swn + mon_sur(mon)%lwn - mon_sur(mon)%sh - mon_sur(mon)%lh - mon_sur(mon)%lsnow
        mon_sur(mon)%lh_l = mon_sur(mon)%lh_l * mon_avg
        mon_sur(mon)%sh_l = mon_sur(mon)%sh_l * mon_avg
        mon_sur(mon)%swn_l = mon_sur(mon)%swn_l * mon_avg
        mon_sur(mon)%lwn_l = mon_sur(mon)%lwn_l * mon_avg
        where (cmn%f_lnd.gt.0._wp)
          mon_sur(mon)%ebal_l= mon_sur(mon)%swn_l + mon_sur(mon)%lwn_l - mon_sur(mon)%sh_l - mon_sur(mon)%lh_l - mon_sur(mon)%lsnow
        elsewhere
          mon_sur(mon)%ebal_l= 0._wp
        endwhere
        mon_sur(mon)%lh_o = mon_sur(mon)%lh_o * mon_avg
        mon_sur(mon)%sh_o = mon_sur(mon)%sh_o * mon_avg
        mon_sur(mon)%swn_o = mon_sur(mon)%swn_o * mon_avg
        mon_sur(mon)%lwn_o = mon_sur(mon)%lwn_o * mon_avg
        where (cmn%f_ocn.gt.0._wp)
          mon_sur(mon)%ebal_o= mon_sur(mon)%swn_o + mon_sur(mon)%lwn_o - mon_sur(mon)%sh_o - mon_sur(mon)%lh_o - mon_sur(mon)%lsnow
        elsewhere
          mon_sur(mon)%ebal_o= 0._wp
        endwhere
        mon_sur(mon)%tskin = mon_sur(mon)%tskin * mon_avg
        mon_sur(mon)%t2m   = mon_sur(mon)%t2m   * mon_avg
        mon_sur(mon)%q2m   = mon_sur(mon)%q2m   * mon_avg
        mon_sur(mon)%alb_dir = mon_sur(mon)%alb_dir * mon_avg
        mon_sur(mon)%alb_dif = mon_sur(mon)%alb_dif * mon_avg
      endif

      if (time_eoy) then

        ann_sur%mask_ice  = cmn%mask_ice
        ann_sur%mask_smb = cmn%mask_smb

        ann_sur%f_stp = cmn%f_stp
        ann_sur%f_ocn = cmn%f_ocn
        ann_sur%f_lnd = cmn%f_lnd
        ann_sur%f_ice = cmn%f_ice

        ! Get annual values
        call surf_ave(mon_sur,ann_sur)
!-----------PET/CWD/MCWD ----------------------------------------------------
        ! Copy annual MCWD from coupler
        ann_sur%MCWD_yearly = cmn%MCWD_yearly      
!--------------------------------------------------------------------------
        ! net Atlantic freshwater flux integrated from the North Pole
        do j=nj,1,-1
          tmp = 0.
          do i=1,ni
            if (basin_mask(i,j).eq.i_atlantic) then
              tmp = tmp + cmn%f_ocn(i,j)*area(i,j)*(ann_sur%prc(i,j)-ann_sur%evp(i,j)+ann_sur%runoff(i,j)+ann_sur%calving(i,j)+ann_sur%bmelt(i,j))/sec_day*1.e-6*1.e-3  ! kg/m2/day * m2 * day/s * Sv/(m3/s) * m3/kg = Sv
            endif
          enddo
          if (j.eq.nj) then
            ann_sur%fw_atl(j) = tmp
          else
            ann_sur%fw_atl(j) = ann_sur%fw_atl(j+1) + tmp
          endif
        enddo

        ! write output
        call cmn_diag_out

      endif

    endif


    return

  end subroutine cmn_diag


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c m n _ d i a g _ o u t
  ! Purpose  :  Initialize netcdf output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine cmn_diag_out

    implicit none

    integer :: k, ncid
    character (len=256) :: fnm


    nout = nout+1


    fnm = trim(out_dir)//"/cmn.nc"
    call nc_open(fnm,ncid)
    call nc_write(fnm,dim_time,dble(year_now),dim1=dim_time,start=[nout],count=[1],ncid=ncid)
    do k = 1, nmon_year
      call surf_nc_write(fnm,ncid,mon_sur(k),k,nout)
    end do
    call surf_nc_write(fnm,ncid,ann_sur,nmon_year+1,nout)
    call nc_close(ncid)


    return
 
   end subroutine cmn_diag_out


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t s _ n c
  ! Purpose  :  Initialize netcdf output 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ts_nc(fnm)

    implicit none

    character (len=*) :: fnm
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_write_dim(fnm,dim_time, x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE.)
    call nc_write_dim(fnm, dim_lat, x=1, axis="y", units="1")
    call nc_write_dim(fnm, dim_lon, x=1, axis="x", units="1")

    return

  end subroutine ts_nc


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t s _ n c _ w r i t e
  ! Purpose  :  Output timestep of netcdf output 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ts_nc_write(fnm,vars,ndat,y)

    implicit none

    type(ts_out) :: vars(:)

    character (len=*) :: fnm
    integer :: ndat, y, ncid, i

    call nc_open(fnm,ncid)
    call nc_write(fnm,"time", dble([(i,i=(year_now-(y-1)*n_accel),(year_now),(n_accel))]), &
    dim1=dim_time,start=[ndat],count=[y],ncid=ncid)    
    call nc_write(fnm,"t2m",     vars%t2m,    dim1=dim_time,start=[ndat],count=[y],long_name="global mean annual surface air temperaure",units="°C",ncid=ncid)
    call nc_write(fnm,"sh",      vars%sh,     dim1=dim_time,start=[ndat],count=[y],long_name="global mean annual surface sensible heat flux",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lh",      vars%lh,     dim1=dim_time,start=[ndat],count=[y],long_name="global mean annual surface latent heat flux",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lwd",     vars%lwd,    dim1=dim_time,start=[ndat],count=[y],long_name="global mean annual downward longwave radiation at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lwu",     vars%lwu,    dim1=dim_time,start=[ndat],count=[y],long_name="global mean annual upward longwave radiation at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"swnet",   vars%swn,    dim1=dim_time,start=[ndat],count=[y],long_name="global mean annual net surface shortwave radiation",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lwnet",   vars%lwn,    dim1=dim_time,start=[ndat],count=[y],long_name="global mean annual net surface longwave radiation",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lsnow",   vars%lsnow,  dim1=dim_time,start=[ndat],count=[y],long_name="global mean annual latent energy of snowfall",units="W/m2",ncid=ncid)
    call nc_write(fnm,"ebal_sur",vars%ebal,   dim1=dim_time,start=[ndat],count=[y],long_name="global mean annual net surface energy balance",units="W/m2",ncid=ncid)
    call nc_write(fnm,"prc",     vars%prc,    dim1=dim_time,start=[ndat],count=[y],long_name="globally integrated annual precipitation",units="10^15 kg/yr",ncid=ncid)
    call nc_write(fnm,"evp",     vars%evp,    dim1=dim_time,start=[ndat],count=[y],long_name="globally integrated annual evapotranspiration",units="10^15 kg/yr",ncid=ncid)
    call nc_write(fnm,"runoff",  vars%runoff, dim1=dim_time,start=[ndat],count=[y],long_name="globally integrated annual runoff into the ocean",units="10^15 kg/yr",ncid=ncid)
    call nc_write(fnm,"runoff_veg",  vars%runoff_veg, dim1=dim_time,start=[ndat],count=[y],long_name="globally integrated annual runoff into the ocean from vegetated part",units="10^15 kg/yr",ncid=ncid)
    call nc_write(fnm,"runoff_ice",  vars%runoff_ice, dim1=dim_time,start=[ndat],count=[y],long_name="globally integrated annual runoff into the ocean from ice sheets",units="10^15 kg/yr",ncid=ncid)
    call nc_write(fnm,"calving", vars%calving,dim1=dim_time,start=[ndat],count=[y],long_name="globally integrated annual calving into the ocean",units="10^15 kg/yr",ncid=ncid)
    call nc_write(fnm,"calving_veg", vars%calving_veg,dim1=dim_time,start=[ndat],count=[y],long_name="globally integrated annual calving into the ocean from vegetated part",units="10^15 kg/yr",ncid=ncid)
    call nc_write(fnm,"calving_ice", vars%calving_ice,dim1=dim_time,start=[ndat],count=[y],long_name="globally integrated annual calving into the ocean from ice sheets",units="10^15 kg/yr",ncid=ncid)
    call nc_write(fnm,"bmelt",  vars%bmelt,   dim1=dim_time,start=[ndat],count=[y],long_name="globally integrated annual basal melt into the ocean",units="10^15 kg/yr",ncid=ncid)
    call nc_write(fnm,"bmelt_grd",  vars%bmelt_grd,   dim1=dim_time,start=[ndat],count=[y],long_name="globally integrated annual basal melt from grouned ice into the ocean",units="10^15 kg/yr",ncid=ncid)
    call nc_write(fnm,"bmelt_flt",  vars%bmelt_flt,   dim1=dim_time,start=[ndat],count=[y],long_name="globally integrated annual basal melt from floating ice into the ocean",units="10^15 kg/yr",ncid=ncid)
    call nc_write(fnm,"prc_atl", vars%prc_atl,dim1=dim_time,start=[ndat],count=[y],long_name="annual precipitation into the Atlantic ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"prc_pac", vars%prc_pac,dim1=dim_time,start=[ndat],count=[y],long_name="annual precipitation into the Pacific Ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"prc_ind", vars%prc_ind,dim1=dim_time,start=[ndat],count=[y],long_name="annual precipitation into the Indian Ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"prc_so",  vars%prc_so, dim1=dim_time,start=[ndat],count=[y],long_name="annual precipitation into the Southern Ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"evp_atl", vars%evp_atl,dim1=dim_time,start=[ndat],count=[y],long_name="annual evaporation from the Atlantic ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"evp_pac", vars%evp_pac,dim1=dim_time,start=[ndat],count=[y],long_name="annual evaporation from the Pacific Ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"evp_ind", vars%evp_ind,dim1=dim_time,start=[ndat],count=[y],long_name="annual evaporation from the Indian Ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"evp_so",  vars%evp_so, dim1=dim_time,start=[ndat],count=[y],long_name="annual evaporation from the Southern Ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_atl", vars%runoff_atl,dim1=dim_time,start=[ndat],count=[y],long_name="annual runoff into the Atlantic ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_pac", vars%runoff_pac,dim1=dim_time,start=[ndat],count=[y],long_name="annual runoff into the Pacific Ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_ind", vars%runoff_ind,dim1=dim_time,start=[ndat],count=[y],long_name="annual runoff into the Indian Ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_so",  vars%runoff_so, dim1=dim_time,start=[ndat],count=[y],long_name="annual runoff into the Southern Ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"calving_atl", vars%calving_atl,dim1=dim_time,start=[ndat],count=[y],long_name="annual calving into the Atlantic ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"calving_pac", vars%calving_pac,dim1=dim_time,start=[ndat],count=[y],long_name="annual calving into the Pacific Ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"calving_ind", vars%calving_ind,dim1=dim_time,start=[ndat],count=[y],long_name="annual calving into the Indian Ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"calving_so",  vars%calving_so, dim1=dim_time,start=[ndat],count=[y],long_name="annual calving into the Southern Ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"bmelt_atl", vars%bmelt_atl,dim1=dim_time,start=[ndat],count=[y],long_name="annual ice sheet basal melt into the Atlantic ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"bmelt_pac", vars%bmelt_pac,dim1=dim_time,start=[ndat],count=[y],long_name="annual ice sheet basal melt into the Pacific Ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"bmelt_ind", vars%bmelt_ind,dim1=dim_time,start=[ndat],count=[y],long_name="annual ice sheet basal melt into the Indian Ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"bmelt_so",  vars%bmelt_so, dim1=dim_time,start=[ndat],count=[y],long_name="annual ice sheet basal melt into the Southern Ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"fw_atl", vars%fw_atl,dim1=dim_time,start=[ndat],count=[y],long_name="annual freshwater flux into the Atlantic ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"fw_natl",vars%fw_natl,dim1=dim_time,start=[ndat],count=[y],long_name="annual freshwater flux into the North Atlantic ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"fw_pac", vars%fw_pac,dim1=dim_time,start=[ndat],count=[y],long_name="annual freshwater flux into the Pacific Ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"fw_ind", vars%fw_ind,dim1=dim_time,start=[ndat],count=[y],long_name="annual freshwater flux into the Indian Ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"fw_so",  vars%fw_so, dim1=dim_time,start=[ndat],count=[y],long_name="annual freshwater flux into the Southern Ocean",units="Sv",ncid=ncid)
    call nc_close(ncid)

   return

  end subroutine ts_nc_write


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t s _ a v e
  ! Purpose  :  Average (or sum) the global time series output as needed
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ts_ave(d,ave)
    
    implicit none
    
    type(ts_out) :: d(:), ave

    integer :: k, n
    real(wp) :: div
    
    n = size(d)
    div = dble(n)
   
    ! Set all values to zero
    ave%t2m = 0._wp
    ave%sh = 0._wp
    ave%lh = 0._wp
    ave%lwd = 0._wp
    ave%lwu = 0._wp
    ave%swn = 0._wp
    ave%lwn = 0._wp
    ave%lsnow = 0._wp
    ave%ebal= 0._wp
    ave%sh_o = 0._wp
    ave%lh_o = 0._wp
    ave%lwd_o = 0._wp
    ave%lwu_o = 0._wp
    ave%swn_o = 0._wp
    ave%lwn_o = 0._wp
    ave%ebal_o= 0._wp
    ave%sh_l = 0._wp
    ave%lh_l = 0._wp
    ave%lwd_l = 0._wp
    ave%lwu_l = 0._wp
    ave%swn_l = 0._wp
    ave%lwn_l = 0._wp
    ave%ebal_l= 0._wp
    ave%runoff = 0._wp
    ave%runoff_veg = 0._wp
    ave%runoff_ice = 0._wp
    ave%calving = 0._wp
    ave%calving_veg = 0._wp
    ave%calving_ice = 0._wp
    ave%bmelt = 0._wp
    ave%bmelt_grd = 0._wp
    ave%bmelt_flt = 0._wp
    ave%evp = 0._wp
    ave%prc = 0._wp
    ave%prc_atl = 0._wp
    ave%prc_pac = 0._wp
    ave%prc_ind = 0._wp
    ave%prc_so  = 0._wp
    ave%evp_atl = 0._wp
    ave%evp_pac = 0._wp
    ave%evp_ind = 0._wp
    ave%evp_so  = 0._wp
    ave%runoff_atl = 0._wp
    ave%runoff_pac = 0._wp
    ave%runoff_ind = 0._wp
    ave%runoff_so  = 0._wp
    ave%calving_atl = 0._wp
    ave%calving_pac = 0._wp
    ave%calving_ind = 0._wp
    ave%calving_so  = 0._wp
    ave%bmelt_atl = 0._wp
    ave%bmelt_pac = 0._wp
    ave%bmelt_ind = 0._wp
    ave%bmelt_so  = 0._wp
    ave%fw_atl = 0._wp
    ave%fw_natl= 0._wp
    ave%fw_pac = 0._wp
    ave%fw_ind = 0._wp
    ave%fw_so  = 0._wp
    ! Loop over the time indices to sum up
    do k=1,n
      ave%t2m = ave%t2m + d(k)%t2m / div
      ave%sh = ave%sh + d(k)%sh / div
      ave%lh = ave%lh + d(k)%lh / div
      ave%lwd = ave%lwd + d(k)%lwd / div
      ave%lwu = ave%lwu + d(k)%lwu / div
      ave%swn = ave%swn + d(k)%swn / div
      ave%lwn = ave%lwn + d(k)%lwn / div
      ave%lsnow = ave%lsnow + d(k)%lsnow / div
      ave%ebal= ave%ebal+ d(k)%ebal / div
      ave%sh_o = ave%sh_o + d(k)%sh_o / div
      ave%lh_o = ave%lh_o + d(k)%lh_o / div
      ave%lwd_o = ave%lwd_o + d(k)%lwd_o / div
      ave%lwu_o = ave%lwu_o + d(k)%lwu_o / div
      ave%swn_o = ave%swn_o + d(k)%swn_o / div
      ave%lwn_o = ave%lwn_o + d(k)%lwn_o / div
      ave%ebal_o= ave%ebal_o+ d(k)%ebal_o / div
      ave%sh_l = ave%sh_l + d(k)%sh_l / div
      ave%lh_l = ave%lh_l + d(k)%lh_l / div
      ave%lwd_l = ave%lwd_l + d(k)%lwd_l / div
      ave%lwu_l = ave%lwu_l + d(k)%lwu_l / div
      ave%swn_l = ave%swn_l + d(k)%swn_l / div
      ave%lwn_l = ave%lwn_l + d(k)%lwn_l / div
      ave%ebal_l= ave%ebal_l+ d(k)%ebal_l / div
      ave%runoff = ave%runoff + d(k)%runoff
      ave%runoff_veg = ave%runoff_veg + d(k)%runoff_veg
      ave%runoff_ice = ave%runoff_ice + d(k)%runoff_ice
      ave%calving = ave%calving + d(k)%calving
      ave%calving_veg = ave%calving_veg + d(k)%calving_veg
      ave%calving_ice = ave%calving_ice + d(k)%calving_ice
      ave%bmelt = ave%bmelt + d(k)%bmelt
      ave%bmelt_grd = ave%bmelt_grd + d(k)%bmelt_grd
      ave%bmelt_flt = ave%bmelt_flt + d(k)%bmelt_flt
      ave%evp = ave%evp + d(k)%evp
      ave%prc = ave%prc + d(k)%prc
      ave%prc_atl = ave%prc_atl + d(k)%prc_atl / div
      ave%prc_pac = ave%prc_pac + d(k)%prc_pac / div
      ave%prc_ind = ave%prc_ind + d(k)%prc_ind / div
      ave%prc_so  = ave%prc_so  + d(k)%prc_so / div
      ave%evp_atl = ave%evp_atl + d(k)%evp_atl / div
      ave%evp_pac = ave%evp_pac + d(k)%evp_pac / div
      ave%evp_ind = ave%evp_ind + d(k)%evp_ind / div
      ave%evp_so  = ave%evp_so  + d(k)%evp_so / div
      ave%runoff_atl = ave%runoff_atl + d(k)%runoff_atl / div
      ave%runoff_pac = ave%runoff_pac + d(k)%runoff_pac / div
      ave%runoff_ind = ave%runoff_ind + d(k)%runoff_ind / div
      ave%runoff_so  = ave%runoff_so  + d(k)%runoff_so / div
      ave%calving_atl = ave%calving_atl + d(k)%calving_atl / div
      ave%calving_pac = ave%calving_pac + d(k)%calving_pac / div
      ave%calving_ind = ave%calving_ind + d(k)%calving_ind / div
      ave%calving_so  = ave%calving_so  + d(k)%calving_so / div
      ave%bmelt_atl = ave%bmelt_atl + d(k)%bmelt_atl / div
      ave%bmelt_pac = ave%bmelt_pac + d(k)%bmelt_pac / div
      ave%bmelt_ind = ave%bmelt_ind + d(k)%bmelt_ind / div
      ave%bmelt_so  = ave%bmelt_so  + d(k)%bmelt_so / div
      ave%fw_atl = ave%fw_atl + d(k)%fw_atl / div
      ave%fw_natl= ave%fw_natl+ d(k)%fw_natl/ div
      ave%fw_pac = ave%fw_pac + d(k)%fw_pac / div
      ave%fw_ind = ave%fw_ind + d(k)%fw_ind / div
      ave%fw_so  = ave%fw_so  + d(k)%fw_so / div
    end do
      
   return
    
  end subroutine ts_ave
  

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s u r f _ n c
  ! Purpose  :  Initialize netcdf output 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine surf_nc(fnm)
    
    implicit none
    
    character (len=*) :: fnm
    integer :: ncid
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,dim_time, x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE.,ncid=ncid)    
    call nc_write_dim(fnm, dim_month, x=1._wp, dx=1._wp, nx=13, axis="e", &
    units="months", ncid=ncid)
    call nc_write_dim(fnm,dim_nsurf, x=1._wp,dx=1._wp, nx=nsurf, units="n/a", &
    axis="z", ncid=ncid)
    call nc_write_dim(fnm, dim_lat, x=lat, axis="y", ncid=ncid)
    call nc_write_dim(fnm, dim_lon, x=lon, axis="x", ncid=ncid)
    call nc_close(ncid)

    return
  
  end subroutine surf_nc
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s u r f _ n c _ w r i t e
  ! Purpose  :  Output timestep of netcdf output 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine surf_nc_write(fnm,ncid,vars,ndat,nout)
    
    implicit none
    
    type(sur_out) :: vars
 
    character (len=*) :: fnm
    integer :: ndat, nout, ncid


    if (ndat.eq.13) then
      if (nout.eq.1) then
        call nc_write(fnm,"mask_ice", vars%mask_ice,dims=[dim_lon,dim_lat],start=[1,1],count=[ni,nj],long_name="ice sheet model mask on coupler grid",units="/",ncid=ncid)
        call nc_write(fnm,"mask_smb", vars%mask_smb,dims=[dim_lon,dim_lat],start=[1,1],count=[ni,nj],long_name="smb model mask on coupler grid",units="/",ncid=ncid)
      endif
      call nc_write(fnm,"f_stp", vars%f_stp,dims=[dim_lon,dim_lat,dim_nsurf,dim_time],start=[1,1,1,nout],count=[ni,nj,nsurf,1],long_name="surface type fractions",units="/",ncid=ncid)
      call nc_write(fnm,"f_ocn", vars%f_ocn,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="ocean fraction",units="/",ncid=ncid)
      call nc_write(fnm,"f_lnd", vars%f_lnd,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="land fraction",units="/",ncid=ncid)
      call nc_write(fnm,"f_ice", vars%f_ice,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="ice fraction",units="/",ncid=ncid)
      call nc_write(fnm,"fw_atl",vars%fw_atl,dims=[dim_lat,dim_time],start=[1,nout],count=[nj,1],long_name="net Atlantic freshwater (P-E+R) integrated from the North Pole",units="Sv",ncid=ncid)
!-----------MCWD -------------------------------------------------
      call nc_write(fnm,"MCWD_yearly", vars%MCWD_yearly,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout], count=[ni,nj,1],long_name="annual maximum climatic water deficit",units="mm", ncid=ncid)
!-----------------------------------------------------------------
    endif

    call nc_write(fnm,"prc",     vars%prc, dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="precipitation",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"rain",    vars%rain,   dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="rainfall",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"snow",    vars%snow,   dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="snowfall",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"evp",     vars%evp,dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="evaporation",units="kg/m2/day",ncid=ncid)
!-----------PET/CWD/MCWD ----------------------------------------------------
    call nc_write(fnm,"pet",     vars%pet, dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="potential evapotranspiration",units="mm/day", ncid=ncid)
    call nc_write(fnm,"CWD",     vars%CWD, dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="cumulative climatic water deficit",units="mm", ncid=ncid)
    call nc_write(fnm,"MCWD",    vars%MCWD,dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="maximum climatic water deficit", units="mm", ncid=ncid)
!----------------------------------------------------------------------------
    call nc_write(fnm,"runoff",  vars%runoff,dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="runoff into the ocean",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"runoff_veg",  vars%runoff_veg,dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="runoff into the ocean from vegetated part",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"runoff_ice",  vars%runoff_ice,dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="runoff into the ocean from ice sheets",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"calving", vars%calving,dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="calving into the ocean",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"calving_veg", vars%calving_veg,dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="calving into the ocean from vegetated part",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"calving_ice", vars%calving_ice,dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="calving into the ocean from ice sheets",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"bmelt",   vars%bmelt,dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="basal melt flux into the ocean",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"bmelt_grd",   vars%bmelt_grd,dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="basal melt flux from grounded ice into the ocean",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"bmelt_flt",   vars%bmelt_flt,dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="basal melt flux from floating ice into the ocean",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"wind",    vars%wind,  dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="surface wind speed",units="m/s",ncid=ncid)
    call nc_write(fnm,"lh",      vars%lh,dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="latent heat flux",units="W/m2",ncid=ncid)
    call nc_write(fnm,"sh",      vars%sh,dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="sensible heat flux",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lwd",     vars%lwd, dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="downward longwave radiation at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lwu",     vars%lwu, dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="upward longwave radiation at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"swnet",   vars%swn, dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="net surface shortwave radiation",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lwnet",   vars%lwn, dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="net surface longwave radiation",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lsnow",   vars%lsnow, dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="snowfall latent energy",units="W/m2",ncid=ncid)
    call nc_write(fnm,"ebal_sur",vars%ebal, dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="net surface energy balance",units="W/m2",ncid=ncid)
    call nc_write(fnm,"ebal_surl",vars%ebal_l, dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="net surface energy balance",units="W/m2",ncid=ncid)
    call nc_write(fnm,"ebal_suro",vars%ebal_o, dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="net surface energy balance",units="W/m2",ncid=ncid)
    call nc_write(fnm,"tskin",   vars%tskin, dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="skin temperature",units="K",ncid=ncid)
    call nc_write(fnm,"t2m",     vars%t2m, dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="surface air temperature",units="K",ncid=ncid)
    call nc_write(fnm,"q2m",     vars%q2m, dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="surface air specific humidity",units="kg/kg",ncid=ncid)
    call nc_write(fnm,"alb_dir", vars%alb_dir, dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="direct beam surface albedo",units="/",ncid=ncid)
    call nc_write(fnm,"alb_dif", vars%alb_dif, dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[ni,nj,1,1],long_name="diffuse radiation surface albedo",units="/",ncid=ncid)


   return              
                       
  end subroutine surf_nc_write
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s u r f _ a v e
  ! Purpose  :  Average (or sum) as needed
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine surf_ave(d,ave)
    
    implicit none
    
    type(sur_out) :: d(:), ave

    integer :: k, n
    real(wp) :: div


    n = size(d)
    div = dble(n)

    ave%prc      = 0._wp
    ave%rain     = 0._wp
    ave%snow     = 0._wp
    ave%evp      = 0._wp
    ave%runoff   = 0._wp
    ave%runoff_veg   = 0._wp
    ave%runoff_ice   = 0._wp
    ave%calving  = 0._wp
    ave%calving_veg  = 0._wp
    ave%calving_ice  = 0._wp
    ave%bmelt    = 0._wp
    ave%bmelt_grd= 0._wp
    ave%bmelt_flt= 0._wp
    ave%wind     = 0._wp
    ave%lh       = 0._wp
    ave%sh       = 0._wp
    ave%lwd      = 0._wp
    ave%lwu      = 0._wp
    ave%swn      = 0._wp
    ave%lwn      = 0._wp
    ave%lsnow    = 0._wp
    ave%ebal     = 0._wp
    ave%ebal_l   = 0._wp
    ave%ebal_o   = 0._wp
    ave%tskin    = 0._wp
    ave%t2m      = 0._wp
    ave%q2m      = 0._wp
    ave%alb_dir  = 0._wp
    ave%alb_dif  = 0._wp
    ! Loop over the time indices to sum up and average (if necessary)
    do k = 1, n
     ave%prc      = ave%prc     + d(k)%prc     / div    
     ave%rain     = ave%rain    + d(k)%rain    / div
     ave%snow     = ave%snow    + d(k)%snow    / div
     ave%evp      = ave%evp     + d(k)%evp     / div
     ave%runoff   = ave%runoff  + d(k)%runoff  / div
     ave%runoff_veg   = ave%runoff_veg  + d(k)%runoff_veg  / div
     ave%runoff_ice   = ave%runoff_ice  + d(k)%runoff_ice  / div
     ave%calving  = ave%calving + d(k)%calving / div
     ave%calving_veg  = ave%calving_veg + d(k)%calving_veg / div
     ave%calving_ice  = ave%calving_ice + d(k)%calving_ice / div
     ave%bmelt    = ave%bmelt   + d(k)%bmelt   / div
     ave%bmelt_grd= ave%bmelt_grd+ d(k)%bmelt_grd   / div
     ave%bmelt_flt= ave%bmelt_flt+ d(k)%bmelt_flt   / div
     ave%wind     = ave%wind    + d(k)%wind    / div
     ave%lh       = ave%lh      + d(k)%lh      / div
     ave%sh       = ave%sh      + d(k)%sh      / div
     ave%lwd      = ave%lwd     + d(k)%lwd     / div
     ave%lwu      = ave%lwu     + d(k)%lwu     / div
     ave%swn      = ave%swn     + d(k)%swn     / div
     ave%lwn      = ave%lwn     + d(k)%lwn     / div
     ave%lsnow    = ave%lsnow   + d(k)%lsnow   / div
     ave%ebal     = ave%ebal    + d(k)%ebal    / div
     ave%ebal_l   = ave%ebal_l  + d(k)%ebal_l  / div
     ave%ebal_o   = ave%ebal_o  + d(k)%ebal_o  / div
     ave%tskin    = ave%tskin   + d(k)%tskin   / div
     ave%t2m      = ave%t2m     + d(k)%t2m     / div
     ave%q2m      = ave%q2m     + d(k)%q2m     / div
     ave%alb_dir  = ave%alb_dir + d(k)%alb_dir / div
     ave%alb_dif  = ave%alb_dif + d(k)%alb_dif / div
   end do
      

   return
    
  end subroutine surf_ave
  

end module cmn_out



