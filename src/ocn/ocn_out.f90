!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : o c n _ out
!
!  Purpose : ocean model diagnostics and output
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Neil R. Edwards and Matteo Willeit
!
! This file is part of CLIMBER-X.
!
! This file was partly ported from the original c-GOLDSTEIN model,
! see Edwards and Marsh (2005)
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
module ocn_out

  use precision, only : wp
  use dim_name, only: len, dim_lon, dim_lat, dim_time, dim_month, dim_day, dim_type, dim_isles

  use timer, only : nmon_year, nstep_mon_ocn, nstep_year_ocn, time_soy_ocn, time_eoy_ocn, &
  time_out_ocn, nyout_ocn, ny_out_ts, y_out_ts_clim, time_out_ts_clim
  use timer, only : n_accel, year, year_clim, year_now, mon, doy, nyears, sec_day, sec_year, nday_year
  use constants, only : pi, cap_w, g
  use control, only: out_dir
  use climber_grid, only : lon, lat, lonu, latv, basin_mask, basin_mask2, i_atlantic, i_pacific, i_indian, i_southern
  use ocn_grid, only: mask_c, mask_v, k1, k1_shelf, k1_1000, k1_3000, topo, bathy ,maxi, maxj, maxk
  use ocn_grid, only : zro, zw, dx, dy, dz, dza, dxv, ocn_area, ocn_area_tot, ocn_vol, ocn_vol_tot, rdy, dphi, phi0, s, sv
  use ocn_grid, only: maxisles, n_isles
  use ocn_grid, only : map_isles, map_edge
  use ocn_params, only : l_daily_output, l_output_extended
  use ocn_params, only: dt, age_tracer, i_age, dye_tracer, i_dye, cons_tracer, i_cons, l_cfc, i_cfc11, i_cfc12, rho0, n_tracers_ocn
  use ocn_params, only: diff_dia, drag, drag_bcl, l_hosing, l_flux_adj_atl, l_flux_adj_ant, l_flux_adj_pac, i_fwf_buoy, i_alphabeta
  use ocn_params, only : l_noise_fw, l_noise_flx
  use ocn_params, only : depth_buoy
  use ocn_params, only : l_diff_dia_strat
  use ocn_params, only : l_bering_flow
  use transport_ocn_mod, only: drho_dx, drho_dy, drho_dz, Ri
  use eos_mod, only : eos
  use ocn_def, only : ocn_class
  use ncio

  !$use omp_lib

  implicit none

  integer :: nout

  real, parameter :: missing_value = -9999.
  real(wp) :: mon_avg, ann_avg
  integer :: k_over
  integer :: k_drho
  integer :: js_drho1, jn_drho1
  integer :: js_drho2, jn1_drho2, jn2_drho2
  integer :: jsf    !! j index where Southern ocean finishes and Atlantic/Pacific/Indian oceans start
  integer :: J26N, jas, jan1, jan2, jps, JNS, k1_buoy
  integer, dimension(:), allocatable :: jan
  integer, dimension(2) :: loc_drake, loc_bering, loc_davis, loc_medi, loc_indo, loc_agulhas
  integer :: j_fram, i_fram(2)
  integer :: j_denmark, i_denmark(2)
  integer :: j_atlN(2), i_atlN(2)
  integer :: j_atlN50(2), i_atlN50(2)
  integer :: j_lab(2), i_lab(2)
  integer :: j_irm(2), i_irm(2)
  integer :: j_gin(2), i_gin(2)
  integer :: j_bkn(2), i_bkn(2)
  integer :: j_wedd(2), i_wedd(2)
  integer :: j_ross(2), i_ross(2)
  integer :: j_so(2), i_so(2)
  integer :: j_soS60(2), i_soS60(2)
  integer :: j_ibe, i_ibe
  integer, parameter :: nlatv_buoy = 6
  real(wp), dimension(nlatv_buoy) :: latv_buoy = (/40._wp,45._wp,50._wp,55._wp,60._wp,65._wp/)
  integer, parameter :: ilatv_buoy_sel = 4
  integer, parameter :: ilatv_buoy_sel2 = 6

  real(wp), dimension(:,:), allocatable :: t_atl, t_pac, t_ind, t_so
  real(wp), dimension(:,:), allocatable :: s_atl, s_pac, s_ind, s_so
  real(wp), dimension(:,:), allocatable :: rho_atl, rho_pac, rho_ind, rho_so

  type ts_out
     integer :: ncells
     real(wp) :: area, vol
     real(wp) :: sst, sss, t, s, cons, age, dye
     real(wp) :: tdocn, tdocn_atl, tdocn_pac, tdocn_ind, tdocn_so
     real(wp), dimension(6) :: ohc, ohc700, ohc2000
     real(wp) :: amoc26N, omaxa, omina, oinfa, omaxp, ominp, omaxs, omins
     real(wp) :: hmaxa, h55a, hmaxs, h60s, hmaxp, fmaxa
     real(wp) :: fov, fovs, fovn
     real(wp) :: faz, fazs, fazn
     real(wp) :: fw_lab
     real(wp) :: buoy, buoyT, buoyS
     real(wp) :: buoy_N, buoyT_N, buoyS_N
     real(wp) :: buoy_tr, buoyT_tr, buoyS_tr
     real(wp) :: buoy_S, buoyT_S, buoyS_S
     real(wp) :: buoyT_NA(nlatv_buoy)
     real(wp) :: buoyS_NA(nlatv_buoy)
     real(wp) :: buoySw_NA(nlatv_buoy)
     real(wp) :: buoySw_lab
     real(wp) :: buoySi_NA(nlatv_buoy)
     real(wp) :: buoy_NA(nlatv_buoy)
     real(wp) :: hosing
     real(wp) :: noise_fw, noise_flx, fw_noise
     real(wp) :: cfc11, cfc12
     real(wp) :: drhoatl1
     real(wp) :: drhoatl2, drhoTatl2, drhoSatl2
     real(wp) :: drhoatl3, drhoTatl3, drhoSatl3
     real(wp) :: saln0, dvsf
     real(wp), dimension(10) :: flx, fw, fw_corr, vsf, p_e_sic, runoff, runoff_veg, runoff_ice, runoff_lake, calving, bmelt, icemelt
     real(wp) :: drake, bering, davis, fram, denmark, medi, indo, agulhas
     real(wp) :: fw_bering, fw_davis, fw_fram, fw_denmark
     real(wp) :: shelf
     real(wp) :: sl_steric
     real(wp) :: mld_atlN50, mld_lab, mld_irm, mld_gin, mld_bkn, mld_wedd, mld_ross, mld_so
     real(wp) :: mldst_atlN50, mldst_lab, mldst_irm, mldst_gin, mldst_bkn, mldst_wedd, mldst_ross, mldst_so
     real(wp) :: pe_atlN, pe_atlN50, pe_lab, pe_irm, pe_gin, pe_bkn, pe_wedd, pe_ross, pe_so
     real(wp) :: buoy_lab, buoy_irm, buoy_gin, buoy_bkn, buoy_wedd, buoy_ross, buoy_so, buoy_soS60, buoyT_soS60, buoyS_soS60
     real(wp) :: t_atlN50, t_lab, t_irm, t_gin, t_bkn, t_wedd, t_ross, t_so
     real(wp) :: s_atlN50, s_lab, s_irm, s_gin, s_bkn, s_wedd, s_ross, s_so
     real(wp) :: t_ibe
  end type

  type o_out
     real(wp), dimension(:,:), allocatable :: f_ocn
     integer, dimension(:,:), allocatable :: mask_ocn
     real(wp), dimension(:,:), allocatable :: area
     integer, dimension(:,:), allocatable :: k1
     real(wp), dimension(:,:), allocatable :: topo
     real(wp), dimension(:,:), allocatable :: bathy
     real(wp), dimension(:,:,:), allocatable :: vol
     integer, dimension(:,:), allocatable :: map_isles
     integer, dimension(:,:,:), allocatable :: map_edge
     real(wp), dimension(:,:), allocatable :: drag
     real(wp), dimension(:,:), allocatable :: drag_bcl
     real(wp), dimension(:,:,:), allocatable :: t, s, age, dye, cons, rho, rho2, u, v, w
     real(wp), dimension(:,:,:), allocatable :: dt_dt_adv, dt_dt_diff
     real(wp), dimension(:,:,:), allocatable :: ds_dt_adv, ds_dt_diff
     real(wp), dimension(:,:), allocatable :: dt_dt_flxsur, ds_dt_flxsur
     real(wp), dimension(:,:), allocatable :: sst, tbot
     real(wp), dimension(:,:), allocatable :: sss, sbot
     real(wp), dimension(:,:,:), allocatable :: fdx, fdy, fdz
     real(wp), dimension(:,:), allocatable :: taux, tauy
     real(wp), dimension(:,:,:), allocatable :: cfc11, cfc12
     real(wp), dimension(:,:,:), allocatable :: diffdia, drho_dx, drho_dy, drho_dz, Ri
     real(wp), dimension(:,:), allocatable :: t_atl, t_pac, t_ind, t_so
     real(wp), dimension(:,:), allocatable :: s_atl, s_pac, s_ind, s_so
     real(wp), dimension(:,:), allocatable :: a_atl, a_pac, a_ind, a_so
     real(wp), dimension(:,:), allocatable :: d_atl, d_pac, d_ind, d_so
     real(wp), dimension(:,:), allocatable :: rho_atl, rho_pac, rho_ind, rho_so
     real(wp), dimension(:,:), allocatable :: drho_0_1000, drho_0_3000
     real(wp), dimension(:,:), allocatable :: ub, vb, psi
     real(wp), dimension(:,:,:), allocatable :: ubisl, vbisl
     real(wp), dimension(:,:), allocatable :: opsi, opsia, opsip, opsii
     real(wp), dimension(:,:), allocatable :: flx, fw, vsf, p_e_sic, runoff, runoffSv, runoffSv_ice, calving, bmelt, fw_hosing, fw_flux_adj
     real(wp), dimension(:,:), allocatable :: fw_noise, flx_noise
     real(wp), dimension(:,:), allocatable :: buoy, buoyS, buoyT
     real(wp), dimension(:,:), allocatable :: hft, hfp, hfa
     real(wp), dimension(:,:), allocatable :: hftz, hfpz, hfaz
     real(wp), dimension(:,:), allocatable :: fwt, fwp, fwa
     real(wp), dimension(:,:), allocatable :: fwtz, fwpz, fwaz
     real(wp), dimension(:,:), allocatable :: fayti, fdyti
     real(wp), dimension(:,:), allocatable :: faysi, fdysi 
     real(wp), dimension(:,:), allocatable :: mld, mldmax, mldst, ke_tau
     real(wp), dimension(:,:), allocatable :: dconv, dven, nconv, kven
     real(wp), dimension(:,:), allocatable :: ssh
     real(wp), dimension(:,:), allocatable :: q_geo
  end type

  type(ts_out), allocatable :: ann_ts(:)
  type(o_out) :: day_o(nday_year), mon_o(nmon_year), ann_o


  private
  public :: ocn_diag_init, ocn_diag

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  o c n _ d i a g _ i n i t
  ! Purpose  :  Initialize netcdf output for ocean
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocn_diag_init

    implicit none

    integer :: i, j, k
    real(wp) :: tv1
    real(wp) :: dist_drake, dist_bering, dist_davis, dist_fram, dist_denmark, dist_medi, dist_indo, dist_agulhas
    real(wp) :: dist_atlN, dist_atlN50, dist_lab, dist_irm, dist_gin, dist_bkn, dist_wedd, dist_ross, dist_so, dist_soS60, dist_ibe
    real(wp), parameter :: lon_drake=-67.5_wp, lat_drake=-60._wp
    real(wp), parameter :: lon_bering=-167.5_wp, lat_bering=67.5_wp
    real(wp), parameter :: lon_davis=-62.5_wp, lat_davis=67.5_wp
    real(wp), parameter :: lon_indo=117.5_wp, lat_indo=-7.5_wp
    real(wp), parameter :: lon_agulhas=22._wp, lat_agulhas=-37.5_wp
    real(wp), parameter :: lon_medi=-2.5_wp, lat_medi=37.5_wp
    real(wp), parameter :: lon_fram_1 = -30., lon_fram_2 = 15., lat_fram = 80.
    real(wp), parameter :: lon_denmark_1 = -47.5, lon_denmark_2 = -17.5, lat_denmark = 67.5
    real(wp), parameter :: lon_atlN_1=-100_wp, lon_atlN_2=100._wp, lat_atlN_1=50._wp , lat_atlN_2=90._wp
    real(wp), parameter :: lon_atlN50_1=-100_wp, lon_atlN50_2=100._wp, lat_atlN50_1=50._wp , lat_atlN50_2=90._wp
    real(wp), parameter :: lon_lab_1=-65_wp, lon_lab_2=-40._wp, lat_lab_1=55._wp , lat_lab_2=70._wp
    real(wp), parameter :: lon_irm_1=-40_wp, lon_irm_2=-20._wp, lat_irm_1=55._wp , lat_irm_2=70._wp
    real(wp), parameter :: lon_gin_1=-20_wp, lon_gin_2=20._wp, lat_gin_1=65._wp , lat_gin_2=80._wp
    real(wp), parameter :: lon_bkn_1=20_wp, lon_bkn_2=80._wp, lat_bkn_1=70._wp , lat_bkn_2=85._wp
    real(wp), parameter :: lon_wedd_1=-80_wp, lon_wedd_2=-10._wp, lat_wedd_1=-85._wp , lat_wedd_2=-65._wp
    real(wp), parameter :: lon_ross_1=-180_wp, lon_ross_2=-140._wp, lat_ross_1=-85._wp , lat_ross_2=-65._wp
    real(wp), parameter :: lon_so_1=-180_wp, lon_so_2=180._wp, lat_so_1=-85._wp , lat_so_2=-55._wp
    real(wp), parameter :: lon_soS60_1=-180_wp, lon_soS60_2=180._wp, lat_soS60_1=-85._wp , lat_soS60_2=-60._wp
    real(wp), parameter :: lon_ibe=-7.5_wp, lat_ibe=37.5_wp


    nout = 0

    ! initialize netcdf files
    call ts_nc(trim(out_dir)//"/ocn_ts.nc")
    call ocn_nc(trim(out_dir)//"/ocn.nc")
    if (l_daily_output) call ocn_daily_nc(trim(out_dir)//"/ocn_daily.nc")

    ! allocate
    allocate(ann_ts(ny_out_ts))

    allocate(ann_o%f_ocn(maxi,maxj))
    allocate(ann_o%mask_ocn(maxi,maxj))
    allocate(ann_o%area(maxi,maxj))
    allocate(ann_o%k1(maxi,maxj))
    allocate(ann_o%map_isles(maxi,maxj))
    allocate(ann_o%map_edge(maxi,maxj,maxisles))
    allocate(ann_o%vol(maxi,maxj,maxk))
    allocate(ann_o%t(maxi,maxj,maxk))
    allocate(ann_o%dt_dt_adv(maxi,maxj,maxk))
    allocate(ann_o%dt_dt_diff(maxi,maxj,maxk))
    allocate(ann_o%dt_dt_flxsur(maxi,maxj))
    allocate(ann_o%sst(maxi,maxj))
    allocate(ann_o%tbot(maxi,maxj))
    allocate(ann_o%t_atl(maxj,maxk))
    allocate(ann_o%t_pac(maxj,maxk))
    allocate(ann_o%t_ind(maxj,maxk))
    allocate(ann_o%t_so(maxj,maxk))
    allocate(ann_o%s(maxi,maxj,maxk))
    allocate(ann_o%ds_dt_adv(maxi,maxj,maxk))
    allocate(ann_o%ds_dt_diff(maxi,maxj,maxk))
    allocate(ann_o%ds_dt_flxsur(maxi,maxj))
    allocate(ann_o%s_atl(maxj,maxk))
    allocate(ann_o%s_pac(maxj,maxk))
    allocate(ann_o%s_ind(maxj,maxk))
    allocate(ann_o%s_so(maxj,maxk))
    allocate(ann_o%sss(maxi,maxj))
    allocate(ann_o%sbot(maxi,maxj))
    allocate(ann_o%a_atl(maxj,maxk))
    allocate(ann_o%a_pac(maxj,maxk))
    allocate(ann_o%a_ind(maxj,maxk))
    allocate(ann_o%a_so(maxj,maxk))
    allocate(ann_o%d_atl(maxj,maxk))
    allocate(ann_o%d_pac(maxj,maxk))
    allocate(ann_o%d_ind(maxj,maxk))
    allocate(ann_o%d_so(maxj,maxk))
    allocate(ann_o%diffdia(maxi,maxj,maxk))
    allocate(ann_o%drho_dx(maxi,maxj,maxk))
    allocate(ann_o%drho_dy(maxi,maxj,maxk))
    allocate(ann_o%drho_dz(maxi,maxj,maxk))
    allocate(ann_o%drho_0_1000(maxi,maxj))
    allocate(ann_o%drho_0_3000(maxi,maxj))
    allocate(ann_o%Ri(maxi,maxj,maxk))
    allocate(ann_o%age(maxi,maxj,maxk))
    allocate(ann_o%dye(maxi,maxj,maxk))
    allocate(ann_o%cons(maxi,maxj,maxk))
    allocate(ann_o%cfc11(maxi,maxj,maxk))
    allocate(ann_o%cfc12(maxi,maxj,maxk))
    allocate(ann_o%rho(maxi,maxj,maxk))
    allocate(ann_o%rho2(maxi,maxj,maxk))
    allocate(ann_o%rho_atl(maxj,maxk))
    allocate(ann_o%rho_pac(maxj,maxk))
    allocate(ann_o%rho_ind(maxj,maxk))
    allocate(ann_o%rho_so(maxj,maxk))
    allocate(ann_o%u(maxi,maxj,maxk))
    allocate(ann_o%v(maxi,maxj,maxk))
    allocate(ann_o%w(maxi,maxj,maxk))
    allocate(ann_o%fdx(maxi,maxj,maxk))
    allocate(ann_o%fdy(maxi,maxj,maxk))
    allocate(ann_o%fdz(maxi,maxj,maxk))
    allocate(ann_o%ub(maxi,maxj))
    allocate(ann_o%vb(maxi,maxj))
    allocate(ann_o%ubisl(maxi,maxj,maxisles))
    allocate(ann_o%vbisl(maxi,maxj,maxisles))
    allocate(ann_o%taux(maxi,maxj))
    allocate(ann_o%tauy(maxi,maxj))
    allocate(ann_o%psi(maxi,maxj))
    allocate(ann_o%opsi(0:maxj,maxk))
    allocate(ann_o%opsia(0:maxj,maxk))
    allocate(ann_o%opsip(0:maxj,maxk))
    allocate(ann_o%opsii(0:maxj,maxk))
    allocate(ann_o%flx(maxi,maxj))
    allocate(ann_o%fw(maxi,maxj))
    allocate(ann_o%vsf(maxi,maxj))
    allocate(ann_o%p_e_sic(maxi,maxj))
    allocate(ann_o%runoff(maxi,maxj))
    allocate(ann_o%runoffSv(maxi,maxj))
    allocate(ann_o%runoffSv_ice(maxi,maxj))
    allocate(ann_o%calving(maxi,maxj))
    allocate(ann_o%bmelt(maxi,maxj))
    allocate(ann_o%fw_hosing(maxi,maxj))
    allocate(ann_o%fw_flux_adj(maxi,maxj))
    allocate(ann_o%fw_noise(maxi,maxj))
    allocate(ann_o%flx_noise(maxi,maxj))
    allocate(ann_o%buoy(maxi,maxj))
    allocate(ann_o%buoyS(maxi,maxj))
    allocate(ann_o%buoyT(maxi,maxj))
    allocate(ann_o%hft(5,maxj))
    allocate(ann_o%hfp(5,maxj))
    allocate(ann_o%hftz(maxj,maxk))
    allocate(ann_o%hfpz(maxj,maxk))
    allocate(ann_o%hfaz(maxj,maxk))
    allocate(ann_o%hfa(5,maxj))
    allocate(ann_o%fwt(5,maxj))
    allocate(ann_o%fwp(5,maxj))
    allocate(ann_o%fwa(5,maxj))
    allocate(ann_o%fwtz(maxj,maxk))
    allocate(ann_o%fwpz(maxj,maxk))
    allocate(ann_o%fwaz(maxj,maxk))
    allocate(ann_o%fayti(maxi,maxj))
    allocate(ann_o%fdyti(maxi,maxj))
    allocate(ann_o%faysi(maxi,maxj))
    allocate(ann_o%fdysi(maxi,maxj))
    allocate(ann_o%nconv(maxi,maxj))
    allocate(ann_o%kven(maxi,maxj))
    allocate(ann_o%dconv(maxi,maxj))
    allocate(ann_o%dven(maxi,maxj))
    allocate(ann_o%mld(maxi,maxj))
    allocate(ann_o%mldmax(maxi,maxj))
    allocate(ann_o%mldst(maxi,maxj))
    allocate(ann_o%ke_tau(maxi,maxj))
    allocate(ann_o%ssh(maxi,maxj))
    allocate(ann_o%q_geo(maxi,maxj))

    do k=1,nmon_year
     allocate(mon_o(k)%t(maxi,maxj,maxk))
     allocate(mon_o(k)%dt_dt_adv(maxi,maxj,maxk))
     allocate(mon_o(k)%dt_dt_diff(maxi,maxj,maxk))
     allocate(mon_o(k)%dt_dt_flxsur(maxi,maxj))
     allocate(mon_o(k)%sst(maxi,maxj))
     allocate(mon_o(k)%tbot(maxi,maxj))
     allocate(mon_o(k)%t_atl(maxj,maxk))
     allocate(mon_o(k)%t_pac(maxj,maxk))
     allocate(mon_o(k)%t_ind(maxj,maxk))
     allocate(mon_o(k)%t_so(maxj,maxk))
     allocate(mon_o(k)%s(maxi,maxj,maxk))
     allocate(mon_o(k)%ds_dt_adv(maxi,maxj,maxk))
     allocate(mon_o(k)%ds_dt_diff(maxi,maxj,maxk))
     allocate(mon_o(k)%ds_dt_flxsur(maxi,maxj))
     allocate(mon_o(k)%sss(maxi,maxj))
     allocate(mon_o(k)%sbot(maxi,maxj))
     allocate(mon_o(k)%s_atl(maxj,maxk))
     allocate(mon_o(k)%s_pac(maxj,maxk))
     allocate(mon_o(k)%s_ind(maxj,maxk))
     allocate(mon_o(k)%s_so(maxj,maxk))
     allocate(mon_o(k)%diffdia(maxi,maxj,maxk))
     allocate(mon_o(k)%drho_dx(maxi,maxj,maxk))
     allocate(mon_o(k)%drho_dy(maxi,maxj,maxk))
     allocate(mon_o(k)%drho_dz(maxi,maxj,maxk))
     allocate(mon_o(k)%drho_0_1000(maxi,maxj))
     allocate(mon_o(k)%drho_0_3000(maxi,maxj))
     allocate(mon_o(k)%Ri(maxi,maxj,maxk))
     allocate(mon_o(k)%rho(maxi,maxj,maxk))
     allocate(mon_o(k)%rho2(maxi,maxj,maxk))
     allocate(mon_o(k)%rho_atl(maxj,maxk))
     allocate(mon_o(k)%rho_pac(maxj,maxk))
     allocate(mon_o(k)%rho_ind(maxj,maxk))
     allocate(mon_o(k)%rho_so(maxj,maxk))
     allocate(mon_o(k)%u(maxi,maxj,maxk))
     allocate(mon_o(k)%v(maxi,maxj,maxk))
     allocate(mon_o(k)%w(maxi,maxj,maxk))
     allocate(mon_o(k)%fdx(maxi,maxj,maxk))
     allocate(mon_o(k)%fdy(maxi,maxj,maxk))
     allocate(mon_o(k)%fdz(maxi,maxj,maxk))
     allocate(mon_o(k)%ub(maxi,maxj))
     allocate(mon_o(k)%vb(maxi,maxj))
     allocate(mon_o(k)%ubisl(maxi,maxj,maxisles))
     allocate(mon_o(k)%vbisl(maxi,maxj,maxisles))
     allocate(mon_o(k)%taux(maxi,maxj))
     allocate(mon_o(k)%tauy(maxi,maxj))
     allocate(mon_o(k)%psi(maxi,maxj))
     allocate(mon_o(k)%opsi(0:maxj,maxk))
     allocate(mon_o(k)%opsia(0:maxj,maxk))
     allocate(mon_o(k)%opsip(0:maxj,maxk))
     allocate(mon_o(k)%opsii(0:maxj,maxk))
     allocate(mon_o(k)%flx(maxi,maxj))
     allocate(mon_o(k)%fw(maxi,maxj))
     allocate(mon_o(k)%vsf(maxi,maxj))
     allocate(mon_o(k)%p_e_sic(maxi,maxj))
     allocate(mon_o(k)%runoff(maxi,maxj))
     allocate(mon_o(k)%runoffSv(maxi,maxj))
     allocate(mon_o(k)%runoffSv_ice(maxi,maxj))
     allocate(mon_o(k)%calving(maxi,maxj))
     allocate(mon_o(k)%bmelt(maxi,maxj))
     allocate(mon_o(k)%fw_hosing(maxi,maxj))
     allocate(mon_o(k)%fw_flux_adj(maxi,maxj))
     allocate(mon_o(k)%fw_noise(maxi,maxj))
     allocate(mon_o(k)%flx_noise(maxi,maxj))
     allocate(mon_o(k)%buoy(maxi,maxj))
     allocate(mon_o(k)%buoyS(maxi,maxj))
     allocate(mon_o(k)%buoyT(maxi,maxj))
     allocate(mon_o(k)%hftz(maxj,maxk))
     allocate(mon_o(k)%hft(5,maxj))
     allocate(mon_o(k)%hfpz(maxj,maxk))
     allocate(mon_o(k)%hfp(5,maxj))
     allocate(mon_o(k)%hfaz(maxj,maxk))
     allocate(mon_o(k)%hfa(5,maxj))
     allocate(mon_o(k)%fwt(5,maxj))
     allocate(mon_o(k)%fwp(5,maxj))
     allocate(mon_o(k)%fwa(5,maxj))
     allocate(mon_o(k)%fwtz(maxj,maxk))
     allocate(mon_o(k)%fwpz(maxj,maxk))
     allocate(mon_o(k)%fwaz(maxj,maxk))
     allocate(mon_o(k)%fayti(maxi,maxj))
     allocate(mon_o(k)%fdyti(maxi,maxj))
     allocate(mon_o(k)%faysi(maxi,maxj))
     allocate(mon_o(k)%fdysi(maxi,maxj))
     allocate(mon_o(k)%nconv(maxi,maxj))
     allocate(mon_o(k)%kven(maxi,maxj))
     allocate(mon_o(k)%dconv(maxi,maxj))
     allocate(mon_o(k)%dven(maxi,maxj))
     allocate(mon_o(k)%mld(maxi,maxj))
     allocate(mon_o(k)%mldst(maxi,maxj))
     allocate(mon_o(k)%ke_tau(maxi,maxj))
     allocate(mon_o(k)%ssh(maxi,maxj))
    enddo

    allocate(t_atl(maxj,maxk))
    allocate(t_pac(maxj,maxk))
    allocate(t_ind(maxj,maxk))
    allocate(t_so(maxj,maxk))
    allocate(s_atl(maxj,maxk))
    allocate(s_pac(maxj,maxk))
    allocate(s_ind(maxj,maxk))
    allocate(s_so(maxj,maxk))
    allocate(rho_atl(maxj,maxk))
    allocate(rho_pac(maxj,maxk))
    allocate(rho_ind(maxj,maxk))
    allocate(rho_so(maxj,maxk))
 
    mon_avg = 1._wp/nstep_mon_ocn
    ann_avg = 1._wp/nstep_year_ocn

    ! find j index separating Southern ocean from the other basins
    jsf = 1
    do while (basin_mask(1,jsf).eq.0 .or. basin_mask(1,jsf).eq.i_southern)
      jsf = jsf+1
    enddo
    ! find j start index of Atlantic ocean 
    do i=1,maxi
      do j=2,maxj
        if (basin_mask(i,j).eq.i_atlantic .and. basin_mask(i,j-1).eq.i_southern) jas = j
      enddo
    enddo
    ! find j start index of Pacific ocean 
    do i=1,maxi
      do j=2,maxj
        if (basin_mask(i,j).eq.i_pacific .and. basin_mask(i,j-1).eq.i_southern) jps = j
      enddo
    enddo
    
    J26N = minloc(abs(sv(:)-sin(pi*26._wp/180.0)),1) + lbound(sv,1) - 1    ! -1 because index of sv array starts from 0!

    jan1 = minloc(abs(sv(:)-sin(pi*70._wp/180.0)),1) + lbound(sv,1) - 1    ! -1 because index of sv array starts from 0!
    jan2 = minloc(abs(sv(:)-sin(pi*75._wp/180.0)),1) + lbound(sv,1) - 1 
    allocate(jan(maxi))
    do i=1,maxi
      if (lon(i).lt.-40._wp) then
        ! Canadian Arctic Archipelago
        jan(i) = jan1
      else 
        ! Fram Strait + Barents Sea
        jan(i) = jan2
      endif
    enddo

    ! find strait locations
    loc_drake(1) = 0
    dist_drake = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_drake/180.0)),2.0*pi).lt.dist_drake) then
          loc_drake(1) = i
          dist_drake = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_drake/180.0)),2.0*pi)
       endif
    enddo
    loc_drake(2) = 0
    dist_drake = 999.
    do j=1,maxj
       if (abs(s(j)-sin(pi*lat_drake/180.0)).lt.dist_drake) then
          loc_drake(2) = j
          dist_drake = abs(s(j)-sin(pi*lat_drake/180.0))
       endif
    enddo

    loc_bering(1) = 0
    dist_bering = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_bering/180.0)),2.0*pi).lt.dist_bering) then
          loc_bering(1) = i
          dist_bering = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_bering/180.0)),2.0*pi)
       endif
    enddo
    loc_bering(2) = 0
    dist_bering = 999.
    do j=1,maxj
       if (abs(s(j)-sin(pi*lat_bering/180.0)).lt.dist_bering) then
          loc_bering(2) = j
          dist_bering = abs(sv(j)-sin(pi*lat_bering/180.0))
       endif
    enddo

    loc_davis(1) = 0
    dist_davis = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_davis/180.0)),2.0*pi).lt.dist_davis) then
          loc_davis(1) = i
          dist_davis = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_davis/180.0)),2.0*pi)
       endif
    enddo
    loc_davis(2) = 0
    dist_davis = 999.
    do j=1,maxj
       if (abs(s(j)-sin(pi*lat_davis/180.0)).lt.dist_davis) then
          loc_davis(2) = j
          dist_davis = abs(sv(j)-sin(pi*lat_davis/180.0))
       endif
    enddo

    i_fram(1) = 0
    dist_fram = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_fram_1/180.0)),2.0*pi).lt.dist_fram) then
          i_fram(1) = i
          dist_fram = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_fram_1/180.0)),2.0*pi)
       endif
    enddo
    i_fram(2) = 0
    dist_fram = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_fram_2/180.0)),2.0*pi).lt.dist_fram) then
          i_fram(2) = i
          dist_fram = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_fram_2/180.0)),2.0*pi)
       endif
    enddo
    j_fram = 0
    dist_fram = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_fram/180.0)).lt.dist_fram) then
          j_fram = j
          dist_fram = abs(sv(j)-sin(pi*lat_fram/180.0))
       endif
    enddo

    i_denmark(1) = 0
    dist_denmark = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_denmark_1/180.0)),2.0*pi).lt.dist_denmark) then
          i_denmark(1) = i
          dist_denmark = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_denmark_1/180.0)),2.0*pi)
       endif
    enddo
    i_denmark(2) = 0
    dist_denmark = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_denmark_2/180.0)),2.0*pi).lt.dist_denmark) then
          i_denmark(2) = i
          dist_denmark = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_denmark_2/180.0)),2.0*pi)
       endif
    enddo
    j_denmark = 0
    dist_denmark = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_denmark/180.0)).lt.dist_denmark) then
          j_denmark = j
          dist_denmark = abs(sv(j)-sin(pi*lat_denmark/180.0))
       endif
    enddo

    i_atlN(1) = 0
    dist_atlN = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_atlN_1/180.0)),2.0*pi).lt.dist_atlN) then
          i_atlN(1) = i
          dist_atlN = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_atlN_1/180.0)),2.0*pi)
       endif
    enddo
    i_atlN(2) = 0
    dist_atlN = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_atlN_2/180.0)),2.0*pi).lt.dist_atlN) then
          i_atlN(2) = i
          dist_atlN = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_atlN_2/180.0)),2.0*pi)
       endif
    enddo
    j_atlN(1) = 0
    dist_atlN = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_atlN_1/180.0)).lt.dist_atlN) then
          j_atlN(1) = j
          dist_atlN = abs(sv(j)-sin(pi*lat_atlN_1/180.0))
       endif
    enddo
    j_atlN(2) = 0
    dist_atlN = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_atlN_2/180.0)).lt.dist_atlN) then
          j_atlN(2) = j
          dist_atlN = abs(sv(j)-sin(pi*lat_atlN_2/180.0))
       endif
    enddo

    i_atlN50(1) = 0
    dist_atlN50 = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_atlN50_1/180.0)),2.0*pi).lt.dist_atlN50) then
          i_atlN50(1) = i
          dist_atlN50 = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_atlN50_1/180.0)),2.0*pi)
       endif
    enddo
    i_atlN50(2) = 0
    dist_atlN50 = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_atlN50_2/180.0)),2.0*pi).lt.dist_atlN50) then
          i_atlN50(2) = i
          dist_atlN50 = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_atlN50_2/180.0)),2.0*pi)
       endif
    enddo
    j_atlN50(1) = 0
    dist_atlN50 = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_atlN50_1/180.0)).lt.dist_atlN50) then
          j_atlN50(1) = j
          dist_atlN50 = abs(sv(j)-sin(pi*lat_atlN50_1/180.0))
       endif
    enddo
    j_atlN50(2) = 0
    dist_atlN50 = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_atlN50_2/180.0)).lt.dist_atlN50) then
          j_atlN50(2) = j
          dist_atlN50 = abs(sv(j)-sin(pi*lat_atlN50_2/180.0))
       endif
    enddo

    i_lab(1) = 0
    dist_lab = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_lab_1/180.0)),2.0*pi).lt.dist_lab) then
          i_lab(1) = i
          dist_lab = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_lab_1/180.0)),2.0*pi)
       endif
    enddo
    i_lab(2) = 0
    dist_lab = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_lab_2/180.0)),2.0*pi).lt.dist_lab) then
          i_lab(2) = i
          dist_lab = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_lab_2/180.0)),2.0*pi)
       endif
    enddo
    j_lab(1) = 0
    dist_lab = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_lab_1/180.0)).lt.dist_lab) then
          j_lab(1) = j
          dist_lab = abs(sv(j)-sin(pi*lat_lab_1/180.0))
       endif
    enddo
    j_lab(2) = 0
    dist_lab = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_lab_2/180.0)).lt.dist_lab) then
          j_lab(2) = j
          dist_lab = abs(sv(j)-sin(pi*lat_lab_2/180.0))
       endif
    enddo

    i_irm(1) = 0
    dist_irm = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_irm_1/180.0)),2.0*pi).lt.dist_irm) then
          i_irm(1) = i
          dist_irm = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_irm_1/180.0)),2.0*pi)
       endif
    enddo
    i_irm(2) = 0
    dist_irm = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_irm_2/180.0)),2.0*pi).lt.dist_irm) then
          i_irm(2) = i
          dist_irm = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_irm_2/180.0)),2.0*pi)
       endif
    enddo
    j_irm(1) = 0
    dist_irm = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_irm_1/180.0)).lt.dist_irm) then
          j_irm(1) = j
          dist_irm = abs(sv(j)-sin(pi*lat_irm_1/180.0))
       endif
    enddo
    j_irm(2) = 0
    dist_irm = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_irm_2/180.0)).lt.dist_irm) then
          j_irm(2) = j
          dist_irm = abs(sv(j)-sin(pi*lat_irm_2/180.0))
       endif
    enddo

    i_gin(1) = 0
    dist_gin = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_gin_1/180.0)),2.0*pi).lt.dist_gin) then
          i_gin(1) = i
          dist_gin = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_gin_1/180.0)),2.0*pi)
       endif
    enddo
    i_gin(2) = 0
    dist_gin = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_gin_2/180.0)),2.0*pi).lt.dist_gin) then
          i_gin(2) = i
          dist_gin = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_gin_2/180.0)),2.0*pi)
       endif
    enddo
    j_gin(1) = 0
    dist_gin = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_gin_1/180.0)).lt.dist_gin) then
          j_gin(1) = j
          dist_gin = abs(sv(j)-sin(pi*lat_gin_1/180.0))
       endif
    enddo
    j_gin(2) = 0
    dist_gin = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_gin_2/180.0)).lt.dist_gin) then
          j_gin(2) = j
          dist_gin = abs(sv(j)-sin(pi*lat_gin_2/180.0))
       endif
    enddo

    i_bkn(1) = 0
    dist_bkn = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_bkn_1/180.0)),2.0*pi).lt.dist_bkn) then
          i_bkn(1) = i
          dist_bkn = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_bkn_1/180.0)),2.0*pi)
       endif
    enddo
    i_bkn(2) = 0
    dist_bkn = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_bkn_2/180.0)),2.0*pi).lt.dist_bkn) then
          i_bkn(2) = i
          dist_bkn = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_bkn_2/180.0)),2.0*pi)
       endif
    enddo
    j_bkn(1) = 0
    dist_bkn = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_bkn_1/180.0)).lt.dist_bkn) then
          j_bkn(1) = j
          dist_bkn = abs(sv(j)-sin(pi*lat_bkn_1/180.0))
       endif
    enddo
    j_bkn(2) = 0
    dist_bkn = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_bkn_2/180.0)).lt.dist_bkn) then
          j_bkn(2) = j
          dist_bkn = abs(sv(j)-sin(pi*lat_bkn_2/180.0))
       endif
    enddo

    i_wedd(1) = 0
    dist_wedd = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_wedd_1/180.0)),2.0*pi).lt.dist_wedd) then
          i_wedd(1) = i
          dist_wedd = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_wedd_1/180.0)),2.0*pi)
       endif
    enddo
    i_wedd(2) = 0
    dist_wedd = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_wedd_2/180.0)),2.0*pi).lt.dist_wedd) then
          i_wedd(2) = i
          dist_wedd = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_wedd_2/180.0)),2.0*pi)
       endif
    enddo
    j_wedd(1) = 0
    dist_wedd = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_wedd_1/180.0)).lt.dist_wedd) then
          j_wedd(1) = j
          dist_wedd = abs(sv(j)-sin(pi*lat_wedd_1/180.0))
       endif
    enddo
    j_wedd(2) = 0
    dist_wedd = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_wedd_2/180.0)).lt.dist_wedd) then
          j_wedd(2) = j
          dist_wedd = abs(sv(j)-sin(pi*lat_wedd_2/180.0))
       endif
    enddo

    i_ross(1) = 0
    dist_ross = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_ross_1/180.0)),2.0*pi).lt.dist_ross) then
          i_ross(1) = i
          dist_ross = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_ross_1/180.0)),2.0*pi)
       endif
    enddo
    i_ross(2) = 0
    dist_ross = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_ross_2/180.0)),2.0*pi).lt.dist_ross) then
          i_ross(2) = i
          dist_ross = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_ross_2/180.0)),2.0*pi)
       endif
    enddo
    j_ross(1) = 0
    dist_ross = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_ross_1/180.0)).lt.dist_ross) then
          j_ross(1) = j
          dist_ross = abs(sv(j)-sin(pi*lat_ross_1/180.0))
       endif
    enddo
    j_ross(2) = 0
    dist_ross = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_ross_2/180.0)).lt.dist_ross) then
          j_ross(2) = j
          dist_ross = abs(sv(j)-sin(pi*lat_ross_2/180.0))
       endif
    enddo

    i_so(1) = 1
    i_so(2) = maxi
    j_so(1) = 0
    dist_so = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_so_1/180.0)).lt.dist_so) then
          j_so(1) = j
          dist_so = abs(sv(j)-sin(pi*lat_so_1/180.0))
       endif
    enddo
    j_so(2) = 0
    dist_so = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_so_2/180.0)).lt.dist_so) then
          j_so(2) = j
          dist_so = abs(sv(j)-sin(pi*lat_so_2/180.0))
       endif
    enddo

    i_soS60(1) = 1
    i_soS60(2) = maxi
    j_soS60(1) = 0
    dist_soS60 = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_soS60_1/180.0)).lt.dist_soS60) then
          j_soS60(1) = j
          dist_soS60 = abs(sv(j)-sin(pi*lat_soS60_1/180.0))
       endif
    enddo
    j_soS60(2) = 0
    dist_soS60 = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_soS60_2/180.0)).lt.dist_soS60) then
          j_soS60(2) = j
          dist_soS60 = abs(sv(j)-sin(pi*lat_soS60_2/180.0))
       endif
    enddo

    i_ibe = 0
    dist_ibe = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_ibe/180.0)),2.0*pi).lt.dist_ibe) then
          i_ibe = i
          dist_ibe = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_ibe/180.0)),2.0*pi)
       endif
    enddo
    j_ibe = 0
    dist_ibe = 999.
    do j=1,maxj
       if (abs(sv(j)-sin(pi*lat_ibe/180.0)).lt.dist_ibe) then
          j_ibe = j
          dist_ibe = abs(sv(j)-sin(pi*lat_ibe/180.0))
       endif
    enddo

    loc_medi(1) = 0
    dist_medi = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_medi/180.0)),2.0*pi).lt.dist_medi) then
          loc_medi(1) = i
          dist_medi = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_medi/180.0)),2.0*pi)
       endif
    enddo
    loc_medi(2) = 0
    dist_medi = 999.
    do j=1,maxj
       if (abs(s(j)-sin(pi*lat_medi/180.0)).lt.dist_medi) then
          loc_medi(2) = j
          dist_medi = abs(s(j)-sin(pi*lat_medi/180.0))
       endif
    enddo

    loc_agulhas(1) = 0
    dist_agulhas = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_agulhas/180.0)),2.0*pi).lt.dist_agulhas) then
          loc_agulhas(1) = i
          dist_agulhas = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_agulhas/180.0)),2.0*pi)
       endif
    enddo
    loc_agulhas(2) = 0
    dist_agulhas = 999.
    do j=1,maxj
       if (abs(s(j)-sin(pi*lat_agulhas/180.0)).lt.dist_agulhas) then
          loc_agulhas(2) = j
          dist_agulhas = abs(s(j)-sin(pi*lat_agulhas/180.0))
       endif
    enddo

    loc_indo(1) = 0
    dist_indo = 999.
    do i=1,maxi
       if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_indo/180.0)),2.0*pi).lt.dist_indo) then
          loc_indo(1) = i
          dist_indo = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_indo/180.0)),2.0*pi)
       endif
    enddo
    loc_indo(2) = 0
    dist_indo = 999.
    do j=1,maxj
       if (abs(s(j)-sin(pi*lat_indo/180.0)).lt.dist_indo) then
          loc_indo(2) = j
          dist_indo = abs(s(j)-sin(pi*lat_indo/180.0))
       endif
    enddo


    ! variable to control the level below which max/min overturning is calculated 
    ! currently it is set for overturning values on levels > 700 m deep
    k_over = maxk
    do k=maxk,1,-1
     tv1 = zw(k)
     if (tv1.gt.-700._wp) then
      k_over = k - 1
     endif
    enddo
    print*,'depth levels for OPSIT max/min overturning : 1 to', k_over

    ! latitudinal bands for atlantic meridional density gradient
    js_drho1 = jas
    jn_drho1 = 0
    tv1 = 999.
    do j=1,maxj
       if (abs(s(j)-sin(pi*(52.5)/180.0)).lt.tv1) then
          jn_drho1 = j
          tv1 = abs(s(j)-sin(pi*(52.5)/180.0))
       endif
    enddo
    js_drho2 = js_drho1
    jn1_drho2 = 0
    tv1 = 999.
    do j=1,maxj
       if (abs(s(j)-sin(pi*(42.5)/180.0)).lt.tv1) then
          jn1_drho2 = j
          tv1 = abs(s(j)-sin(pi*(42.5)/180.0))
       endif
    enddo
    jn2_drho2 = 0
    tv1 = 999.
    do j=1,maxj
       if (abs(s(j)-sin(pi*(57.5)/180.0)).lt.tv1) then
          jn2_drho2 = j
          tv1 = abs(s(j)-sin(pi*(57.5)/180.0))
       endif
    enddo
    print *,'js_drho1,jn1_drho2,jn2_drho2',js_drho1,jn1_drho2,jn2_drho2

    ! depth level for atlantic meridional density gradient 
    k_drho = maxk
    do k=maxk,1,-1
     tv1 = zw(k)
     if (tv1.gt.-700._wp) then
      k_drho = k - 1
     endif
    enddo

    ! define depth of 'mixed layer' for North Atlantic buoyancy diagnostics
    k1_buoy = maxk
    tv1 = 1000._wp
    do k=maxk,1,-1
      if (abs(zw(k-1)+depth_buoy).lt.tv1) then
        k1_buoy = k
        tv1 = abs(zw(k-1)+depth_buoy)
      endif
    enddo
    print *,'k1_buoy',k1_buoy
    print *,'depth buoy',zw(k1_buoy-1),zro(k1_buoy)

   return

  end subroutine ocn_diag_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  o c n _ d i a g
  !   Purpose    :  ocean diagnostics
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocn_diag(ocn)

    implicit none

    type(ocn_class), intent(inout) :: ocn

    integer :: i, j, k, kr, l, m, n, y, ntot, nx, nxa, nxp
    real(wp) :: tv2, tv3
    real(wp) :: sum_2d(n_tracers_ocn), sum_3d(n_tracers_ocn)
    real(wp) :: opsi(0:maxj,0:maxk), ou(maxj,maxk)
    real(wp) :: opsia(0:maxj,0:maxk)
    real(wp) :: opsip(0:maxj,0:maxk)
    real(wp) :: opsii(0:maxj,0:maxk)
    integer :: iposa(2)
    real(wp), dimension(10) :: fw, fw_corr, p_e_sic, runoff, runoff_veg, runoff_ice, runoff_lake, calving, bmelt, icemelt, flx, vsf
    real(wp) :: fw_noise
    real(wp) :: hft(3,maxj), hfp(3,maxj), hfa(3,maxj)
    real(wp) :: hftz(maxj,maxk), hfpz(maxj,maxk), hfaz(maxj,maxk)
    real(wp) :: fwt(3,maxj), fwp(3,maxj), fwa(3,maxj)
    real(wp) :: fwtz(maxj,maxk), fwpz(maxj,maxk), fwaz(maxj,maxk)
    real(wp) :: vz(maxj,maxk), vza(maxj,maxk), vzp(maxj,maxk) 
    real(wp) :: tz(maxj,maxk), tza(maxj,maxk), tzp(maxj,maxk) 
    real(wp) :: sz(maxj,maxk), sza(maxj,maxk), szp(maxj,maxk) 
    real(wp) :: dxz(maxj,maxk), dxza(maxj,maxk), dxzp(maxj,maxk) 
    real(wp) :: hfto(maxj), hfao(maxj), hfpo(maxj), hftg(maxj), hfag(maxj), hfpg(maxj)
    real(wp) :: fwto(maxj), fwao(maxj), fwpo(maxj), fwtg(maxj), fwag(maxj), fwpg(maxj)
    real(wp) :: global_tocn, global_socn, global_sst, global_sss, global_cons, global_age, global_dye
    real(wp), dimension(6) :: ohc, ohc700, ohc2000
    real(wp) :: tdocn, tdocn_atl, tdocn_pac, tdocn_ind, tdocn_so
    real(wp) :: ocnvol_tot, ocnvol_atl, ocnvol_pac, ocnvol_ind, ocnvol_so
    real(wp) :: global_cfc11, global_cfc12
    real(wp), save :: ocn_area_tot0
    real(wp), save :: sl_steric0
    real(wp) :: sl_steric
    logical :: int_drake, int_bering, int_davis, int_medi, int_indo, int_agulhas
    real(wp) :: tf_drake, tf_bering, tf_davis, tf_fram, tf_denmark, tf_medi, tf_indo, tf_agulhas
    real(wp) :: fw_bering, fw_davis, fw_fram, fw_denmark
    real(wp) :: area_atlN50, area_lab, area_irm, area_gin, area_bkn, area_wedd, area_ross, area_so
    real(wp) :: mldst, rho_k, rho_maxk
    real(wp) :: mld_atlN50, mld_lab, mld_irm, mld_gin, mld_bkn, mld_wedd, mld_ross, mld_so
    real(wp) :: mldst_atlN50, mldst_lab, mldst_irm, mldst_gin, mldst_bkn, mldst_wedd, mldst_ross, mldst_so
    real(wp) :: pe_atlN, pe_atlN50, pe_lab, pe_irm, pe_gin, pe_bkn, pe_wedd, pe_ross, pe_so
    real(wp) :: buoy_lab, buoy_irm, buoy_gin, buoy_bkn, buoy_wedd, buoy_ross, buoy_so, buoy_soS60, buoyT_soS60, buoyS_soS60
    real(wp) :: t_atlN50, t_lab, t_irm, t_gin, t_bkn, t_wedd, t_ross, t_sos
    real(wp) :: s_atlN50, s_lab, s_irm, s_gin, s_bkn, s_wedd, s_ross, s_sos
    real(wp) :: rho_s1, rho_n1
    real(wp) :: rho_b2, rho_n2, rho_n3, rhoT_b2, rhoT_n2, rhoT_n3, rhoS_b2, rhoS_n2, rhoS_n3, area_b, area_n, area_n3, area_ij
    real(wp) :: ocnvol
    real(wp) :: shelf
    integer :: bmask, bmask2
    real(wp) :: fov, fovs, fovn, faz, fazn, fazs, vbar, dxdz, vzab, vzam, sza1, sza2, fazz
    real(wp) :: buoy, buoyT, buoyS
    real(wp) :: buoy_N, buoyT_N, buoyS_N
    real(wp) :: buoy_tr, buoyT_tr, buoyS_tr
    real(wp) :: buoy_S, buoyT_S, buoyS_S
    real(wp) :: buoy_NA(nlatv_buoy), buoyT_NA(nlatv_buoy), buoyS_NA(nlatv_buoy)
    real(wp) :: buoySw_NA(nlatv_buoy), buoySi_NA(nlatv_buoy), buoySw_lab, fw_lab
    real(wp) :: alpha, beta, rho1, rho2
    real(wp), allocatable, dimension(:,:,:) :: rho


    ! current index
    y = y_out_ts_clim

    sum_2d = 0.0
    sum_3d = 0._wp
    ohc = 0._wp
    ohc700 = 0._wp
    ohc2000 = 0._wp
    tdocn = 0._wp
    tdocn_atl = 0._wp
    tdocn_pac = 0._wp
    tdocn_ind = 0._wp
    tdocn_so = 0._wp
    ocnvol_tot = 0._wp
    ocnvol_atl = 0._wp
    ocnvol_pac = 0._wp
    ocnvol_ind = 0._wp
    ocnvol_so = 0._wp
    global_cons = 0._wp
    global_age = 0._wp
    global_dye = 0._wp
    global_cfc11 = 0.0_wp
    global_cfc12 = 0.0_wp
    sl_steric = 0._wp
    fw = 0._wp
    fw_corr = 0._wp
    fw_noise = 0._wp
    p_e_sic = 0._wp
    runoff = 0._wp
    runoff_veg = 0._wp
    runoff_ice = 0._wp
    runoff_lake = 0._wp
    icemelt = 0._wp
    calving = 0._wp
    bmelt = 0._wp
    vsf = 0._wp
    flx = 0._wp
    hft = 0._wp
    hfp = 0._wp
    hfa = 0._wp    
    hftz = 0._wp
    hfpz = 0._wp
    hfaz = 0._wp    
    fwt = 0._wp
    fwp = 0._wp
    fwa = 0._wp
    fwtz = 0._wp
    fwpz = 0._wp
    fwaz = 0._wp

    if (year.eq.1) then
      ocn_area_tot0 = ocn_area_tot
    endif

    !$omp parallel do collapse(2) private(i,j,k,l,ocnvol,tv2,tv3,bmask,bmask2) &
    !$omp reduction(+:sum_2d,sum_3d,ohc,ohc700,ohc2000,tdocn,tdocn_atl,tdocn_pac,tdocn_ind,tdocn_so,ocnvol_tot) &
    !$omp reduction(+:ocnvol_atl,ocnvol_pac,ocnvol_ind,ocnvol_so,global_cons,global_age,global_dye,global_cfc11,global_cfc12) &
    !$omp reduction(+:sl_steric,fw,fw_corr,fw_noise,p_e_sic,runoff,runoff_veg,runoff_ice,runoff_lake,icemelt,calving,bmelt,vsf,flx) &
    !$omp reduction(+:hft,hfp,hfa,hftz,hfpz,hfaz,fwt,fwp,fwa,fwtz,fwpz,fwaz)
    do i=1,maxi
      do j=1,maxj
        bmask = basin_mask(i,j)
        bmask2 = basin_mask2(i,j)
        if (ocn%f_ocn(i,j).gt.0._wp) then
          ! mean SST and SSS
          do l=1,n_tracers_ocn
            sum_2d(l) = sum_2d(l) + ocn%ts(i,j,maxk,l)*dx(j)*dy*ocn%f_ocn(i,j)
          enddo
          ! freshwater and heat fluxes (1=global, 2=Atlantic, 3=Pacific, 4=Indian, 5=Southern)
          ! Global ocean
          fw(1) = fw(1) + ((ocn%p_e_sic(i,j)+ocn%runoff(i,j)+ocn%calving(i,j)+ocn%bmelt(i,j)+ocn%fw_hosing(i,j)+ocn%fw_flux_adj(i,j))*ocn%grid%ocn_area(i,j))  ! kg/m2/s * m2 = kg/s
          fw_corr(1) = fw_corr(1) + (ocn%fw_corr(i,j)*ocn%grid%ocn_area(i,j))  ! kg/m2/s * m2 = kg/s
          fw_noise = fw_noise + (ocn%fw_noise(i,j)*ocn%grid%ocn_area(i,j))  ! kg/m2/s * m2 = kg/s
          p_e_sic(1) = p_e_sic(1) + ocn%p_e_sic(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
          runoff(1) = runoff(1) + ocn%runoff(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
          runoff_veg(1) = runoff_veg(1) + ocn%runoff_veg(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
          runoff_ice(1) = runoff_ice(1) + ocn%runoff_ice(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
          runoff_lake(1) = runoff_lake(1) + ocn%runoff_lake(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
          calving(1) = calving(1) + ocn%calving(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
          bmelt(1) = bmelt(1) + ocn%bmelt(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
          icemelt(1) = icemelt(1) + ocn%melt_ice(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
          vsf(1) = vsf(1) + ocn%flx_sur(i,j,2)*rho0/ocn%saln0*ocn%grid%ocn_area(i,j)  ! m/s*psu * kg/m3 / psu * m2 = kg/s
          flx(1) = flx(1) + (ocn%flx(i,j)*ocn%grid%ocn_area(i,j))   ! W/m2 * m2 = W
          ! Atlantic basin cells
          if (bmask.eq.i_atlantic) then
            fw(2) = fw(2) + ((ocn%p_e_sic(i,j)+ocn%runoff(i,j)+ocn%calving(i,j)+ocn%bmelt(i,j)+ocn%fw_hosing(i,j)+ocn%fw_flux_adj(i,j))*ocn%grid%ocn_area(i,j))
            p_e_sic(2) = p_e_sic(2) + ocn%p_e_sic(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            runoff(2) = runoff(2) + ocn%runoff(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            runoff_veg(2) = runoff_veg(2) + ocn%runoff_veg(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            runoff_ice(2) = runoff_ice(2) + ocn%runoff_ice(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            runoff_lake(2) = runoff_lake(2) + ocn%runoff_lake(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            calving(2) = calving(2) + ocn%calving(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            bmelt(2) = bmelt(2) + ocn%bmelt(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            icemelt(2) = icemelt(2) + ocn%melt_ice(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            vsf(2) = vsf(2) + ocn%flx_sur(i,j,2)*rho0/ocn%saln0*ocn%grid%ocn_area(i,j)  ! m/s*psu * kg/m3 / psu * m2 = kg/s
            flx(2) = flx(2) + (ocn%flx(i,j)*ocn%grid%ocn_area(i,j))
            if (lat(j).gt.0._wp) then
              ! North Atlantic north of 30N
              fw(8) = fw(8) + ((ocn%p_e_sic(i,j)+ocn%runoff(i,j)+ocn%calving(i,j)+ocn%bmelt(i,j)+ocn%fw_hosing(i,j)+ocn%fw_flux_adj(i,j))*ocn%grid%ocn_area(i,j))
              p_e_sic(8) = p_e_sic(8) + ocn%p_e_sic(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              runoff(8) = runoff(8) + ocn%runoff(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              runoff_veg(8) = runoff_veg(8) + ocn%runoff_veg(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              runoff_ice(8) = runoff_ice(8) + ocn%runoff_ice(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              runoff_lake(8) = runoff_lake(8) + ocn%runoff_lake(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              calving(8) = calving(8) + ocn%calving(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              bmelt(8) = bmelt(8) + ocn%bmelt(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              icemelt(8) = icemelt(8) + ocn%melt_ice(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              vsf(8) = vsf(8) + ocn%flx_sur(i,j,2)*rho0/ocn%saln0*ocn%grid%ocn_area(i,j)  ! m/s*psu * kg/m3 / psu * m2 = kg/s
              flx(8) = flx(8) + (ocn%flx(i,j)*ocn%grid%ocn_area(i,j))
            endif
            if (lat(j).gt.30._wp) then
              ! North Atlantic north of 30N
              fw(6) = fw(6) + ((ocn%p_e_sic(i,j)+ocn%runoff(i,j)+ocn%calving(i,j)+ocn%bmelt(i,j)+ocn%fw_hosing(i,j)+ocn%fw_flux_adj(i,j))*ocn%grid%ocn_area(i,j))
              p_e_sic(6) = p_e_sic(6) + ocn%p_e_sic(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              runoff(6) = runoff(6) + ocn%runoff(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              runoff_veg(6) = runoff_veg(6) + ocn%runoff_veg(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              runoff_ice(6) = runoff_ice(6) + ocn%runoff_ice(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              runoff_lake(6) = runoff_lake(6) + ocn%runoff_lake(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              calving(6) = calving(6) + ocn%calving(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              bmelt(6) = bmelt(6) + ocn%bmelt(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              icemelt(6) = icemelt(6) + ocn%melt_ice(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              vsf(6) = vsf(6) + ocn%flx_sur(i,j,2)*rho0/ocn%saln0*ocn%grid%ocn_area(i,j)  ! m/s*psu * kg/m3 / psu * m2 = kg/s
              flx(6) = flx(6) + (ocn%flx(i,j)*ocn%grid%ocn_area(i,j))
            endif
            if (lat(j).gt.50._wp) then
              ! North Atlantic north of 50N
              fw(7) = fw(7) + ((ocn%p_e_sic(i,j)+ocn%runoff(i,j)+ocn%calving(i,j)+ocn%bmelt(i,j)+ocn%fw_hosing(i,j)+ocn%fw_flux_adj(i,j))*ocn%grid%ocn_area(i,j))
              p_e_sic(7) = p_e_sic(7) + ocn%p_e_sic(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              runoff(7) = runoff(7) + ocn%runoff(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              runoff_veg(7) = runoff_veg(7) + ocn%runoff_veg(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              runoff_ice(7) = runoff_ice(7) + ocn%runoff_ice(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              runoff_lake(7) = runoff_lake(7) + ocn%runoff_lake(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              calving(7) = calving(7) + ocn%calving(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              bmelt(7) = bmelt(7) + ocn%bmelt(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              icemelt(7) = icemelt(7) + ocn%melt_ice(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              vsf(7) = vsf(7) + ocn%flx_sur(i,j,2)*rho0/ocn%saln0*ocn%grid%ocn_area(i,j)  ! m/s*psu * kg/m3 / psu * m2 = kg/s
              flx(7) = flx(7) + (ocn%flx(i,j)*ocn%grid%ocn_area(i,j))
            endif
            if (lat(j).gt.55._wp) then
              ! North Atlantic north of 50N
              fw(9) = fw(9) + ((ocn%p_e_sic(i,j)+ocn%runoff(i,j)+ocn%calving(i,j)+ocn%bmelt(i,j)+ocn%fw_hosing(i,j)+ocn%fw_flux_adj(i,j))*ocn%grid%ocn_area(i,j))
              p_e_sic(9) = p_e_sic(9) + ocn%p_e_sic(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              runoff(9) = runoff(9) + ocn%runoff(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              runoff_veg(9) = runoff_veg(9) + ocn%runoff_veg(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              runoff_ice(9) = runoff_ice(9) + ocn%runoff_ice(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              runoff_lake(9) = runoff_lake(9) + ocn%runoff_lake(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              calving(9) = calving(9) + ocn%calving(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              bmelt(9) = bmelt(9) + ocn%bmelt(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              icemelt(9) = icemelt(9) + ocn%melt_ice(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              vsf(9) = vsf(9) + ocn%flx_sur(i,j,2)*rho0/ocn%saln0*ocn%grid%ocn_area(i,j)  ! m/s*psu * kg/m3 / psu * m2 = kg/s
              flx(9) = flx(9) + (ocn%flx(i,j)*ocn%grid%ocn_area(i,j))
            endif
            ! Pacific basin cells
          else if (bmask.eq.i_pacific) then
            fw(3) = fw(3) + ((ocn%p_e_sic(i,j)+ocn%runoff(i,j)+ocn%calving(i,j)+ocn%bmelt(i,j)+ocn%fw_hosing(i,j)+ocn%fw_flux_adj(i,j))*ocn%grid%ocn_area(i,j))
            p_e_sic(3) = p_e_sic(3) + ocn%p_e_sic(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            runoff(3) = runoff(3) + ocn%runoff(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            runoff_veg(3) = runoff_veg(3) + ocn%runoff_veg(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            runoff_ice(3) = runoff_ice(3) + ocn%runoff_ice(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            runoff_lake(3) = runoff_lake(3) + ocn%runoff_lake(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            calving(3) = calving(3) + ocn%calving(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            bmelt(3) = bmelt(3) + ocn%bmelt(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            icemelt(3) = icemelt(3) + ocn%melt_ice(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            vsf(3) = vsf(3) + ocn%flx_sur(i,j,2)*rho0/ocn%saln0*ocn%grid%ocn_area(i,j)  ! m/s*psu * kg/m3 / psu * m2 = kg/s
            flx(3) = flx(3) + (ocn%flx(i,j)*ocn%grid%ocn_area(i,j))
            ! Indian basin cells
          else if (bmask.eq.i_indian) then
            fw(4) = fw(4) + ((ocn%p_e_sic(i,j)+ocn%runoff(i,j)+ocn%calving(i,j)+ocn%bmelt(i,j)+ocn%fw_hosing(i,j)+ocn%fw_flux_adj(i,j))*ocn%grid%ocn_area(i,j))
            p_e_sic(4) = p_e_sic(4) + ocn%p_e_sic(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            runoff(4) = runoff(4) + ocn%runoff(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            runoff_veg(4) = runoff_veg(4) + ocn%runoff_veg(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            runoff_ice(4) = runoff_ice(4) + ocn%runoff_ice(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            runoff_lake(4) = runoff_lake(4) + ocn%runoff_lake(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            calving(4) = calving(4) + ocn%calving(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            bmelt(4) = bmelt(4) + ocn%bmelt(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            icemelt(4) = icemelt(4) + ocn%melt_ice(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            vsf(4) = vsf(4) + ocn%flx_sur(i,j,2)*rho0/ocn%saln0*ocn%grid%ocn_area(i,j)  ! m/s*psu * kg/m3 / psu * m2 = kg/s
            flx(4) = flx(4) + (ocn%flx(i,j)*ocn%grid%ocn_area(i,j))
            ! Southern basin cells
          else if (bmask.eq.i_southern) then
            fw(5) = fw(5) + ((ocn%p_e_sic(i,j)+ocn%runoff(i,j)+ocn%calving(i,j)+ocn%bmelt(i,j)+ocn%fw_hosing(i,j)+ocn%fw_flux_adj(i,j))*ocn%grid%ocn_area(i,j))
            p_e_sic(5) = p_e_sic(5) + ocn%p_e_sic(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            runoff(5) = runoff(5) + ocn%runoff(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            runoff_veg(5) = runoff_veg(5) + ocn%runoff_veg(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            runoff_ice(5) = runoff_ice(5) + ocn%runoff_ice(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            runoff_lake(5) = runoff_lake(5) + ocn%runoff_lake(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            calving(5) = calving(5) + ocn%calving(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            bmelt(5) = bmelt(5) + ocn%bmelt(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            icemelt(5) = icemelt(5) + ocn%melt_ice(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
            vsf(5) = vsf(5) + ocn%flx_sur(i,j,2)*rho0/ocn%saln0*ocn%grid%ocn_area(i,j)  ! m/s*psu * kg/m3 / psu * m2 = kg/s
            flx(5) = flx(5) + (ocn%flx(i,j)*ocn%grid%ocn_area(i,j))
            if (lat(j).lt.-60._wp) then
              ! Southern Ocean south of 60S
              fw(10) = fw(10) + ((ocn%p_e_sic(i,j)+ocn%runoff(i,j)+ocn%calving(i,j)+ocn%bmelt(i,j)+ocn%fw_hosing(i,j)+ocn%fw_flux_adj(i,j))*ocn%grid%ocn_area(i,j))
              p_e_sic(10) = p_e_sic(10) + ocn%p_e_sic(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              runoff(10) = runoff(10) + ocn%runoff(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              runoff_veg(10) = runoff_veg(10) + ocn%runoff_veg(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              runoff_ice(10) = runoff_ice(10) + ocn%runoff_ice(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              runoff_lake(10) = runoff_lake(10) + ocn%runoff_lake(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              calving(10) = calving(10) + ocn%calving(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              bmelt(10) = bmelt(10) + ocn%bmelt(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              icemelt(10) = icemelt(10) + ocn%melt_ice(i,j)*ocn%grid%ocn_area(i,j)  ! kg/m2/s * m2 = kg/s
              vsf(10) = vsf(10) + ocn%flx_sur(i,j,2)*rho0/ocn%saln0*ocn%grid%ocn_area(i,j)  ! m/s*psu * kg/m3 / psu * m2 = kg/s
              flx(10) = flx(10) + (ocn%flx(i,j)*ocn%grid%ocn_area(i,j))
            endif
          endif
          do k=k1(i,j),maxk
            ocnvol = ocn_vol(i,j,k)
            ohc(1) = ohc(1) + ocn%ts(i,j,k,1)*ocnvol*cap_w*rho0     ! K*m3 * J/kg/K * kg/m3 = J
            if (zro(k).ge.-700._wp) ohc700(1) = ohc700(1) + ocn%ts(i,j,k,1)*ocnvol*cap_w*rho0     ! K*m3 * J/kg/K * kg/m3 = J
            if (zro(k).ge.-2000._wp) ohc2000(1) = ohc2000(1) + ocn%ts(i,j,k,1)*ocnvol*cap_w*rho0  ! K*m3 * J/kg/K * kg/m3 = J
            if (bmask.eq.i_atlantic) then
              ohc(2) = ohc(2) + ocn%ts(i,j,k,1)*ocnvol*cap_w*rho0     ! K*m3 * J/kg/K * kg/m3 = J
              if (zro(k).ge.-700._wp) ohc700(2) = ohc700(2) + ocn%ts(i,j,k,1)*ocnvol*cap_w*rho0     ! K*m3 * J/kg/K * kg/m3 = J
              if (zro(k).ge.-2000._wp) ohc2000(2) = ohc2000(2) + ocn%ts(i,j,k,1)*ocnvol*cap_w*rho0  ! K*m3 * J/kg/K * kg/m3 = J
              if (lat(j).gt.30._wp) then
                ohc(6) = ohc(6) + ocn%ts(i,j,k,1)*ocnvol*cap_w*rho0     ! K*m3 * J/kg/K * kg/m3 = J
                if (zro(k).ge.-700._wp) ohc700(6) = ohc700(6) + ocn%ts(i,j,k,1)*ocnvol*cap_w*rho0     ! K*m3 * J/kg/K * kg/m3 = J
                if (zro(k).ge.-2000._wp) ohc2000(6) = ohc2000(6) + ocn%ts(i,j,k,1)*ocnvol*cap_w*rho0  ! K*m3 * J/kg/K * kg/m3 = J
              endif
            else if (bmask.eq.i_pacific) then
              ohc(3) = ohc(3) + ocn%ts(i,j,k,1)*ocnvol*cap_w*rho0     ! K*m3 * J/kg/K * kg/m3 = J
              if (zro(k).ge.-700._wp) ohc700(3) = ohc700(3) + ocn%ts(i,j,k,1)*ocnvol*cap_w*rho0     ! K*m3 * J/kg/K * kg/m3 = J
              if (zro(k).ge.-2000._wp) ohc2000(3) = ohc2000(3) + ocn%ts(i,j,k,1)*ocnvol*cap_w*rho0  ! K*m3 * J/kg/K * kg/m3 = J
            else if (bmask.eq.i_indian) then
              ohc(4) = ohc(4) + ocn%ts(i,j,k,1)*ocnvol*cap_w*rho0     ! K*m3 * J/kg/K * kg/m3 = J
              if (zro(k).ge.-700._wp) ohc700(4) = ohc700(4) + ocn%ts(i,j,k,1)*ocnvol*cap_w*rho0     ! K*m3 * J/kg/K * kg/m3 = J
              if (zro(k).ge.-2000._wp) ohc2000(4) = ohc2000(4) + ocn%ts(i,j,k,1)*ocnvol*cap_w*rho0  ! K*m3 * J/kg/K * kg/m3 = J
            else if (bmask.eq.i_southern) then
              ohc(5) = ohc(5) + ocn%ts(i,j,k,1)*ocnvol*cap_w*rho0     ! K*m3 * J/kg/K * kg/m3 = J
              if (zro(k).ge.-700._wp) ohc700(5) = ohc700(5) + ocn%ts(i,j,k,1)*ocnvol*cap_w*rho0     ! K*m3 * J/kg/K * kg/m3 = J
              if (zro(k).ge.-2000._wp) ohc2000(5) = ohc2000(5) + ocn%ts(i,j,k,1)*ocnvol*cap_w*rho0  ! K*m3 * J/kg/K * kg/m3 = J
            endif
          enddo
          do k=k1(i,j),maxk
            ! deep ocean temperatures
            if (zro(k).le.-2500._wp) then
              ocnvol = ocn_vol(i,j,k)
              ocnvol_tot = ocnvol_tot + ocnvol
              tdocn = tdocn + ocn%ts(i,j,k,1)*ocnvol
              if (bmask.eq.i_atlantic) then
                ocnvol_atl = ocnvol_atl + ocnvol
                tdocn_atl = tdocn_atl + ocn%ts(i,j,k,1)*ocnvol
              else if (bmask.eq.i_pacific) then
                ocnvol_pac = ocnvol_pac + ocnvol
                tdocn_pac = tdocn_pac + ocn%ts(i,j,k,1)*ocnvol
              else if (bmask.eq.i_indian) then
                ocnvol_ind = ocnvol_ind + ocnvol
                tdocn_ind = tdocn_ind + ocn%ts(i,j,k,1)*ocnvol
              else if (bmask.eq.i_southern) then
                ocnvol_so = ocnvol_so + ocnvol
                tdocn_so = tdocn_so + ocn%ts(i,j,k,1)*ocnvol
              endif
            endif
          enddo
        endif
        ! northward heat transport
        if (j.lt.maxj .and. k1(i,j).le.maxk .and. k1(i,j+1).le.maxk) then
          do k=k1(i,j),maxk
            kr = maxk-k+1
            ! advective transport
            tv2 = ocn%fay(i,j,k,1)*cap_w*ocn%rho(i,j,k)/dt  ! W
            ! diffusive transport
            tv3 = ocn%fdy(i,j,k,1)*cap_w*ocn%rho(i,j,k)/dt  ! W
            hftz(j,kr) = hftz(j,kr) + (tv2 + tv3)/dz(k)
            hft(1,j) = hft(1,j) + tv2 + tv3
            hft(2,j) = hft(2,j) + tv2
            hft(3,j) = hft(3,j) + tv3
            if (bmask2.eq.i_pacific .or. bmask2.eq.i_indian) then
              hfpz(j,kr) = hfpz(j,kr) + (tv2 + tv3)/dz(k)
              hfp(1,j) = hfp(1,j) + tv2 + tv3
              hfp(2,j) = hfp(2,j) + tv2
              hfp(3,j) = hfp(3,j) + tv3
            else if (bmask2.eq.i_atlantic) then
              hfaz(j,kr) = hfaz(j,kr) + (tv2 + tv3)/dz(k)
              hfa(1,j) = hfa(1,j) + tv2 + tv3
              hfa(2,j) = hfa(2,j) + tv2
              hfa(3,j) = hfa(3,j) + tv3
            endif
          enddo
        endif
        ! northward freshwater transport
        if (j.lt.maxj .and. k1(i,j).le.maxk .and. k1(i,j+1).le.maxk) then
          do k=k1(i,j),maxk
            kr = maxk-k+1
            ! advective transport
            tv2 = - 1._wp/ocn%saln0*ocn%fay(i,j,k,2)/dt  ! m3/s
            ! diffusive transport
            tv3 = - 1._wp/ocn%saln0*ocn%fdy(i,j,k,2)/dt  ! m3/s
            fwtz(j,kr) = fwtz(j,kr) + (tv2 + tv3)/dz(k)    
            fwt(1,j) = fwt(1,j) + tv2 + tv3
            fwt(2,j) = fwt(2,j) + tv2
            fwt(3,j) = fwt(3,j) + tv3
            if (bmask2.eq.i_pacific .or. bmask2.eq.i_indian) then
              fwpz(j,kr) = fwpz(j,kr) + (tv2 + tv3)/dz(k)
              fwp(1,j) = fwp(1,j) + tv2 + tv3
              fwp(2,j) = fwp(2,j) + tv2
              fwp(3,j) = fwp(3,j) + tv3
            else if (bmask2.eq.i_atlantic) then
              fwaz(j,kr) = fwaz(j,kr) + (tv2 + tv3)/dz(k)
              fwa(1,j) = fwa(1,j) + tv2 + tv3
              fwa(2,j) = fwa(2,j) + tv2
              fwa(3,j) = fwa(3,j) + tv3
            endif
          enddo
        endif

        do k=1,maxk
          if (mask_c(i,j,k).eq.1) then
            ocnvol = ocn_vol(i,j,k)
            ! temperature and salinity integrals
            do l=1,2
              sum_3d(l) = sum_3d(l) + ocn%ts(i,j,k,l)*ocnvol
            enddo
            ! CFC integrals
            if (l_cfc) then
              global_cfc11 = global_cfc11 + ocn%ts(i,j,k,i_cfc11) * ocnvol ! mol
              global_cfc12 = global_cfc12 + ocn%ts(i,j,k,i_cfc12) * ocnvol ! mol
            endif
            ! conservative tracer
            if (cons_tracer) global_cons = global_cons + ocn%ts(i,j,k,i_cons)*ocnvol
            ! age tracer
            if (age_tracer) global_age = global_age + ocn%ts(i,j,k,i_age)*ocnvol/ocn_vol_tot
            ! dye tracer
            if (dye_tracer) global_dye = global_dye + ocn%ts(i,j,k,i_dye)*ocnvol/ocn_vol_tot
            ! steric sea level
            sl_steric = sl_steric + ocnvol*rho0/ocn%rho(i,j,k) / ocn_area_tot0     ! m
          endif
        enddo
      enddo
    enddo
    !$omp end parallel do

    if (year.eq.1) then
      sl_steric0 = sl_steric
    endif

    global_sst = sum_2d(1)/ocn_area_tot
    global_sss = sum_2d(2)/ocn_area_tot
    global_tocn = sum_3d(1)/ocn_vol_tot
    global_socn = sum_3d(2)/ocn_vol_tot

    tdocn     = tdocn/ocnvol_tot
    tdocn_atl = tdocn_atl/ocnvol_atl
    tdocn_pac = tdocn_pac/ocnvol_pac
    tdocn_ind = tdocn_ind/ocnvol_ind
    tdocn_so  = tdocn_so /ocnvol_so 


    ! calculate meridional overturning streamfunctions

    ! initialize
    opsi  = 0._wp
    opsia = 0._wp
    opsip = 0._wp
    opsii = 0._wp

    ! global
    do j=1,maxj-1
       do k=1,maxk
          ou(j,k) = 0._wp
          do i=1,maxi
             ou(j,k) = ou(j,k) + ocn%u(2,i,j,k)
          enddo
          opsi(j,k) = opsi(j,k-1) - ou(j,k)*dxv(j)*dz(k)  ! m3/s
       enddo
    enddo

    ! Pacific
    do j=1,maxj-1
       do k=1,maxk
          ou(j,k) = 0
             do i=1,maxi
               if (basin_mask(i,j).eq.i_pacific) then
                ou(j,k) = ou(j,k) + ocn%u(2,i,j,k)
              endif
             enddo
          opsip(j,k) = opsip(j,k-1) - ou(j,k)*dxv(j)*dz(k)
       enddo
    enddo

    ! Indian
    do j=1,maxj-1
       do k=1,maxk
          ou(j,k) = 0
             do i=1,maxi
               if (basin_mask(i,j).eq.i_indian) then
                ou(j,k) = ou(j,k) + ocn%u(2,i,j,k)
              endif
             enddo
          opsii(j,k) = opsii(j,k-1) - ou(j,k)*dxv(j)*dz(k)
       enddo
    enddo

    ! Atlantic
    do j=1,maxj-1
       do k=1,maxk
          ou(j,k) = 0
          do i=1,maxi
            if (basin_mask(i,j).eq.i_atlantic) then
             ou(j,k) = ou(j,k) + ocn%u(2,i,j,k)
           endif
          enddo
          opsia(j,k) = opsia(j,k-1) - ou(j,k)*dxv(j)*dz(k)
       enddo
    enddo

    ! freshwater fluxes by the MOC and gyre into the Atlantic at its southern border (Liu 2017, eq 1)
    j = jas    ! index of south atlantic border

    vbar = 0._wp
    dxdz = 0._wp
    do k=1,maxk
      do i=1,maxi
        if (basin_mask(i,j).eq.i_atlantic) then
          if (mask_v(i,j,k).eq.1) then
            ! barotropic
            vbar = vbar + ocn%u(2,i,j,k)*dxv(j)*dz(k)
            dxdz = dxdz + dxv(j)*dz(k)
          endif
        endif
      enddo
    enddo
    vbar = vbar/dxdz

    fovs = 0._wp
    fazs = 0._wp
    do k=1,maxk
      vzab = 0._wp
      vzam = 0._wp
      sza1 = 0._wp
      sza2 = 0._wp
      nxa  = 0
      do i=1,maxi
        if (basin_mask(i,j).eq.i_atlantic) then
          if (mask_v(i,j,k).eq.1) then
            ! zonal integral of baroclinic meridional velocity
            vzab = vzab + (ocn%u(2,i,j,k)-ocn%ub(2,i,j))*dxv(j) ! m2/s
            ! zonal mean meridional velocity
            vzam = vzam + ocn%u(2,i,j,k) ! m/s
            ! zonal mean salinity
            sza1 = sza1 + ocn%ts(i,j,k,2)
            sza2 = sza2 + ocn%ts(i,j+1,k,2)
            nxa = nxa+1
          endif
        endif
      enddo
      if (nxa>0) sza1 = sza1/nxa
      if (nxa>0) sza2 = sza2/nxa
      if (nxa>0) vzam = vzam/nxa
      !vzab = (vzam-vbar)*nxa*dxv(j)
      !print *,vzab/(nxa*dxv(j)),vzam-vbar
      fazz = 0._wp
      do i=1,maxi
        if (basin_mask(i,j).eq.i_atlantic) then
          if (mask_v(i,j,k).eq.1) then
            ! zonal integral of v'S'
            fazz = fazz + (ocn%u(2,i,j,k)-vzam) * 0.5_wp*(ocn%ts(i,j,k,2)-sza1 + ocn%ts(i,j+1,k,2)-sza2) *dxv(j) ! m2/s
          endif
        endif
      enddo
      ! integrate vertically
      fovs = fovs - vzab*0.5_wp*(sza1+sza2)*dz(k)    ! m3/s * psu
      fazs = fazs - fazz*dz(k)    ! m3/s * psu
    enddo
    fovs = fovs/ocn%saln0    ! m3/s
    fazs = fazs/ocn%saln0    ! m3/s

    ! freshwater fluxes by the MOC and gyre out of the Atlantic into the Arctic (Liu 2017, eq 1)
    fovn = 0._wp
    fazn = 0._wp
    do k=1,maxk
      vzab = 0._wp
      vzam = 0._wp
      sza1 = 0._wp
      sza2 = 0._wp
      nxa  = 0
      do i=1,maxi
        j = jan(i)
        if (lon(i).gt.-125._wp .and. lon(i).lt.100._wp .and. basin_mask(i,j).eq.i_atlantic) then
          if (mask_v(i,j,k).eq.1) then
            ! zonal integral of baroclinic meridional velocity
            vzab = vzab + (ocn%u(2,i,j,k)-ocn%ub(2,i,j))*dxv(j) ! m2/s
            ! zonal mean meridional velocity
            vzam = vzam + ocn%u(2,i,j,k) ! m/s
            ! zonal mean salinity
            sza1 = sza1 + ocn%ts(i,j,k,2)
            sza2 = sza2 + ocn%ts(i,j+1,k,2)
            nxa = nxa+1
          endif
        endif
      enddo
      if (nxa>0) sza1 = sza1/nxa
      if (nxa>0) sza2 = sza2/nxa
      if (nxa>0) vzam = vzam/nxa
      fazz = 0._wp
      do i=1,maxi
        j = jan(i)
        if (lon(i).gt.-125._wp .and. lon(i).lt.100._wp .and. basin_mask(i,j).eq.i_atlantic) then
          if (mask_v(i,j,k).eq.1) then
            ! zonal integral of v'S'
            fazz = fazz + (ocn%u(2,i,j,k)-vzam) * 0.5_wp*(ocn%ts(i,j,k,2)-sza1 + ocn%ts(i,j+1,k,2)-sza2) *dxv(j) ! m2/s
          endif
        endif
      enddo
      ! integrate vertically
      fovn = fovn - vzab*0.5_wp*(sza1+sza2)*dz(k)    ! m3/s * psu
      fazn = fazn - fazz*dz(k)    ! m3/s * psu
    enddo
    fovn = fovn/ocn%saln0    ! m3/s
    fazn = fazn/ocn%saln0    ! m3/s

    fov = fovs-fovn
    faz = fazs-fazn

    ! global buoyancy flux
    buoyT = 0._wp
    buoyS = 0._wp
    buoyT_N = 0._wp
    buoyS_N = 0._wp
    buoyT_tr = 0._wp
    buoyS_tr = 0._wp
    buoyT_S = 0._wp
    buoyS_S = 0._wp
    do j=1,maxj
      do i=1,maxi
        if (ocn%f_ocn(i,j).gt.0._wp) then
          if (i_alphabeta.eq.1) then
            rho1 = eos(ocn%ts(i,j,maxk,1)+0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
            rho2 = eos(ocn%ts(i,j,maxk,1)-0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
            alpha = rho2-rho1
          else if (i_alphabeta.eq.2) then
            alpha = (519._wp+122._wp*ocn%ts(i,j,maxk,1))*1.e-4_wp    ! kg/m3/K
          endif
          buoyT = buoyT + g/cap_w*ocn%flx(i,j)*alpha/rho0*ocn%grid%ocn_area(i,j)*dt  ! m/s2 * kg*K/J * J/s/m2 * kg/m3/K * m3/kg * m2 *s = kg*m/s2 = N
          if (lat(j).gt.30._wp) then
            buoyT_N = buoyT_N + g/cap_w*ocn%flx(i,j)*alpha/rho0*ocn%grid%ocn_area(i,j)*dt  ! m/s2 * kg*K/J * J/s/m2 * kg/m3/K * m3/kg * m2 *s = kg*m/s2 = N
          endif
          if (lat(j).lt.-30._wp) then
            buoyT_S = buoyT_S + g/cap_w*ocn%flx(i,j)*alpha/rho0*ocn%grid%ocn_area(i,j)*dt  ! m/s2 * kg*K/J * J/s/m2 * kg/m3/K * m3/kg * m2 *s = kg*m/s2 = N
          endif
          if (lat(j).gt.-30._wp .and. lat(j).lt.30._wp) then
            buoyT_tr = buoyT_tr + g/cap_w*ocn%flx(i,j)*alpha/rho0*ocn%grid%ocn_area(i,j)*dt  ! m/s2 * kg*K/J * J/s/m2 * kg/m3/K * m3/kg * m2 *s = kg*m/s2 = N
          endif
          if (i_alphabeta.eq.1) then
            rho1 = eos(ocn%ts(i,j,maxk,1),ocn%ts(i,j,maxk,2)+0.5_wp,0._wp)    
            rho2 = eos(ocn%ts(i,j,maxk,1),ocn%ts(i,j,maxk,2)-0.5_wp,0._wp)    
            beta = rho1-rho2
          else if (i_alphabeta.eq.2) then
            beta = 0.8_wp   ! kg/m3/psu
          endif
          if (i_fwf_buoy.eq.1) then 
            buoyS = buoyS + g*ocn%fw_corr(i,j)*beta*ocn%saln0/rho0*ocn%grid%ocn_area(i,j)*dt  ! m/s2 * kg/m2/s kg/m3/psu * m3/kg * psu * m2 *s = kg*m/s2 = N
            if (lat(j).gt.30._wp) then
              buoyS_N = buoyS_N + g*ocn%fw_corr(i,j)*beta*ocn%saln0/rho0*ocn%grid%ocn_area(i,j)*dt  ! m/s2 * kg/m2/s kg/m3/psu * m3/kg * psu * m2 *s = kg*m/s2 = N
            endif
            if (lat(j).lt.-30._wp) then
              buoyS_S = buoyS_S + g*ocn%fw_corr(i,j)*beta*ocn%saln0/rho0*ocn%grid%ocn_area(i,j)*dt  ! m/s2 * kg/m2/s kg/m3/psu * m3/kg * psu * m2 *s = kg*m/s2 = N
            endif
            if (lat(j).gt.-30._wp .and. lat(j).lt.30._wp) then
              buoyS_tr = buoyS_tr + g*ocn%fw_corr(i,j)*beta*ocn%saln0/rho0*ocn%grid%ocn_area(i,j)*dt  ! m/s2 * kg/m2/s kg/m3/psu * m3/kg * psu * m2 *s = kg*m/s2 = N
            endif
          else if (i_fwf_buoy.eq.2) then
            buoyS = buoyS + g*ocn%flx_sur(i,j,2)*beta*ocn%grid%ocn_area(i,j)*dt  ! m/s2 * m/s*psu * kg/m3/psu * m2 *s = kg*m/s2 = N
            if (lat(j).gt.30._wp) then
              buoyS_N = buoyS_N + g*ocn%flx_sur(i,j,2)*beta*ocn%grid%ocn_area(i,j)*dt  ! m/s2 * m/s*psu * kg/m3/psu * m2 *s = kg*m/s2 = N
            endif
            if (lat(j).lt.-30._wp) then
              buoyS_S = buoyS_S + g*ocn%flx_sur(i,j,2)*beta*ocn%grid%ocn_area(i,j)*dt  ! m/s2 * m/s*psu * kg/m3/psu * m2 *s = kg*m/s2 = N
            endif
            if (lat(j).gt.-30._wp .and. lat(j).lt.30._wp) then
              buoyS_tr = buoyS_tr + g*ocn%flx_sur(i,j,2)*beta*ocn%grid%ocn_area(i,j)*dt  ! m/s2 * m/s*psu * kg/m3/psu * m2 *s = kg*m/s2 = N
            endif
          endif
        endif
      enddo
    enddo
    buoy = buoyT + buoyS
    buoy_N = buoyT_N + buoyS_N
    buoy_tr = buoyT_tr + buoyS_tr
    buoy_S = buoyT_S + buoyS_S

    ! buoyancy balance of the North Atlantic at different latitudes
    do n=1,nlatv_buoy
      JNS = minloc(abs(sv(:)-sin(pi*latv_buoy(n)/180.0)),1) + lbound(sv,1) - 1    ! -1 because index of sv array starts from 0!
      ! thermal effect
      buoyT_NA(n) = 0._wp
      do j=JNS+1,maxj
        do i=1,maxi
          if (ocn%f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_atlantic) then
            if (i_alphabeta.eq.1) then
              rho1 = eos(ocn%ts(i,j,maxk,1)+0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
              rho2 = eos(ocn%ts(i,j,maxk,1)-0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
              alpha = rho2-rho1
            else if (i_alphabeta.eq.2) then
              alpha = (519._wp+122._wp*ocn%ts(i,j,maxk,1))*1.e-4_wp    ! kg/m3/K
            endif
            !print *,'T,S,alpha,alpha1',ocn%ts(i,j,maxk,1),ocn%ts(i,j,maxk,2),alpha,rho2-rho1
            buoyT_NA(n) = buoyT_NA(n) + g/cap_w*ocn%flx(i,j)*alpha/rho0*ocn%grid%ocn_area(i,j)*dt  ! m/s2 * kg*K/J * J/s/m2 * kg/m3/K * m3/kg * m2 *s = kg*m/s2 = N
          endif
        enddo
      enddo
      ! haline effect
      buoyS_NA(n) = 0._wp
      do j=JNS+1,maxj
        do i=1,maxi
          if (ocn%f_ocn(i,j).gt.0._wp .and. basin_mask(i,j).eq.i_atlantic) then
            if (i_alphabeta.eq.1) then
              rho1 = eos(ocn%ts(i,j,maxk,1),ocn%ts(i,j,maxk,2)+0.5_wp,0._wp)    
              rho2 = eos(ocn%ts(i,j,maxk,1),ocn%ts(i,j,maxk,2)-0.5_wp,0._wp)    
              beta = rho1-rho2
            else if (i_alphabeta.eq.2) then
              beta = 0.8_wp   ! kg/m3/psu
            endif
            !print *,'beta,beta1',beta,rho2-rho1
            if (i_fwf_buoy.eq.1) then 
              buoyS_NA(n) = buoyS_NA(n) + g*ocn%fw_corr(i,j)*beta*ocn%saln0/rho0*ocn%grid%ocn_area(i,j)*dt  ! m/s2 * kg/m2/s kg/m3/psu * m3/kg * psu * m2 *s = kg*m/s2 = N
            else if (i_fwf_buoy.eq.2) then
              buoyS_NA(n) = buoyS_NA(n) + g*ocn%flx_sur(i,j,2)*beta*ocn%grid%ocn_area(i,j)*dt  ! m/s2 * m/s*psu * kg/m3/psu * m2 *s = kg*m/s2 = N
            endif
          endif
        enddo
      enddo
      ! buoyancy from meridional salinity flux relative to saln0 in the mixed layer
      buoySw_NA(n) = 0._wp
      beta = 0.8_wp   ! kg/m3/psu
      J = JNS
      do k=maxk,k1_buoy,-1
        do i=1,maxi
          if (lon(i).gt.-70. .and. basin_mask(i,j).eq.i_atlantic .and. mask_v(i,J,k).eq.1) then
            buoySw_NA(n) = buoySw_NA(n) - g*beta*(0.5*(ocn%ts(i,j,k,2)+ocn%ts(i,j+1,k,2))-ocn%saln0)*ocn%u(2,i,J,k)*dxv(J)*dz(k)*dt    
            ! m/s2 * kg/m3/psu * psu * m/s*m2*s = kg*m/s2 = N
          endif
        enddo
      enddo
      ! buoyancy from sea ice export
      buoySi_NA(n) = ocn%buoy_sic_NA(n)
      ! sum
      buoy_NA(n) = buoyT_NA(n) + buoyS_NA(n) 
    enddo
    
    ! Labrador current only
    fw_lab = 0._wp
    buoySw_lab = 0._wp
    J = minloc(abs(sv(:)-sin(pi*latv_buoy(ilatv_buoy_sel)/180.0)),1) + lbound(sv,1) - 1    ! -1 because index of sv array starts from 0!
    do k=maxk,k1_buoy,-1
      do i=1,maxi
        if (lon(i).gt.-70. .and. lon(i).lt.-30. .and. basin_mask(i,j).eq.i_atlantic .and. mask_v(i,J,k).eq.1) then
          fw_lab = fw_lab + (0.5*(ocn%ts(i,j,k,2)+ocn%ts(i,j+1,k,2))-ocn%saln0)/ocn%saln0*ocn%u(2,i,J,k)*dxv(J)*dz(k)    ! psu/psu * m/s * m2 = m3/s
          buoySw_lab = buoySw_lab - g*beta*(0.5*(ocn%ts(i,j,k,2)+ocn%ts(i,j+1,k,2))-ocn%saln0)*ocn%u(2,i,J,k)*dxv(J)*dz(k)*dt    
            ! m/s2 * kg/m3/psu * psu * m/s*m2*s = kg*m/s2 = N
        endif
      enddo
    enddo
    fw_lab = fw_lab*1.e-6_wp  ! Sv


    ! 1. Drake Passage through flow (Sv)
    tf_drake = 0.0
    int_drake = .true.
    do j=loc_drake(2),maxj
       if (int_drake.and.k1(loc_drake(1),j).le.maxk) then
          do k=k1(loc_drake(1),j),maxk
             tf_drake = tf_drake + ocn%u(1,loc_drake(1),j,k)*dy*dz(k) ! m3/s
!	print *,'u',ocn%u(1,loc_drake(1),j,k),loc_drake(1),j,k
          enddo
       else
          int_drake = .false.
       endif
    enddo
    do j=loc_drake(2)-1,1,-1
       if (k1(loc_drake(1),j).le.maxk) then
          do k=k1(loc_drake(1),j),maxk
             tf_drake = tf_drake + ocn%u(1,loc_drake(1),j,k)*dy*dz(k) ! m3/s
!        print *,'u',ocn%u(1,loc_drake(1),j,k),loc_drake(1),j,k
          enddo
       endif
    enddo
    tf_drake = tf_drake*1.e-6_wp  ! Sv
!        print *,'loc_drake',loc_drake
!        print *,'tf_drake',tf_drake
!        print *,''


    ! 2. Bering strait through flow (Sv)
    tf_bering = 0.0
    fw_bering = 0.0
    int_bering = .true.
    do i=loc_bering(1),maxi
       if (int_bering.and.k1(i,loc_bering(2)).le.maxk) then
          do k=k1(i,loc_bering(2)),maxk
             tf_bering = tf_bering + ocn%u(2,i,loc_bering(2),k)*dxv(loc_bering(2))*dz(k) ! m3/s
             fw_bering = fw_bering + ocn%u(2,i,loc_bering(2),k)*(1.0-0.5*(ocn%ts(i,loc_bering(2)+1,k,2) + ocn%ts(i,loc_bering(2),k,2))/ocn%saln0) &
                       * dxv(loc_bering(2))*dz(k) ! m3/s
!	print *,'u',ocn%u(2,i,loc_bering(2),k),i,loc_bering(2),k
          enddo
       else
          int_bering = .false.
       endif
    enddo
    int_bering = .true.
    do i=loc_bering(1)-1,1,-1
       if (int_bering.and.k1(i,loc_bering(2)).le.maxk) then
          do k=k1(i,loc_bering(2)),maxk
             tf_bering = tf_bering + ocn%u(2,i,loc_bering(2),k)*dxv(loc_bering(2))*dz(k) ! m3/s
             fw_bering = fw_bering + ocn%u(2,i,loc_bering(2),k)*(1.0-0.5*(ocn%ts(i,loc_bering(2)+1,k,2) + ocn%ts(i,loc_bering(2),k,2))/ocn%saln0) &
                       * dxv(loc_bering(2))*dz(k) ! m3/s
!        print *,'u',ocn%u(2,i,loc_bering(2),k),i,loc_bering(2),k
          enddo
       else
          int_bering = .false.
       endif
    enddo
    tf_bering = tf_bering*1.e-6_wp  ! Sv
    fw_bering = fw_bering*1.e-6_wp  ! Sv
    if (l_bering_flow) then ! if parameterisation for throughflow enabled
      tf_bering = ocn%bering_tf
      fw_bering = ocn%bering_fw
    endif

!        print *,'loc_drake',loc_drake
!        print *,'tf_drake',tf_drake
!        print *,''

    ! 3. Davis strait through flow (Sv)
    tf_davis = 0.0
    fw_davis = 0.0
    int_davis = .true.
    do i=loc_davis(1),maxi
       if (int_davis.and.k1(i,loc_davis(2)).le.maxk) then
          do k=k1(i,loc_davis(2)),maxk
             tf_davis = tf_davis + ocn%u(2,i,loc_davis(2),k)*dxv(loc_davis(2))*dz(k) ! m3/s
             fw_davis = fw_davis + ocn%u(2,i,loc_davis(2),k)*(1.0-0.5*(ocn%ts(i,loc_davis(2)+1,k,2) + ocn%ts(i,loc_davis(2),k,2))/ocn%saln0) &
                       * dxv(loc_davis(2))*dz(k) ! m3/s
!	print *,'u',ocn%u(2,i,loc_davis(2),k),i,loc_davis(2),k
          enddo
       else
          int_davis = .false.
       endif
    enddo
    int_davis = .true.
    do i=loc_davis(1)-1,1,-1
       if (int_davis.and.k1(i,loc_davis(2)).le.maxk) then
          do k=k1(i,loc_davis(2)),maxk
             tf_davis = tf_davis + ocn%u(2,i,loc_davis(2),k)*dxv(loc_davis(2))*dz(k) ! m3/s
             fw_davis = fw_davis + ocn%u(2,i,loc_davis(2),k)*(1.0-0.5*(ocn%ts(i,loc_davis(2)+1,k,2) + ocn%ts(i,loc_davis(2),k,2))/ocn%saln0) &
                       * dxv(loc_davis(2))*dz(k) ! m3/s
!        print *,'u',ocn%u(2,i,loc_davis(2),k),i,loc_davis(2),k
          enddo
       else
          int_davis = .false.
       endif
    enddo
    tf_davis = tf_davis*1.e-6_wp  ! Sv
    fw_davis = fw_davis*1.e-6_wp  ! Sv
!        print *,'loc_davis',loc_davis
!        print *,'tf_davis',tf_davis
!        print *,''

    ! 4. Fram strait through flow (Sv)
    tf_fram = 0.0
    fw_fram = 0.0
    do i=i_fram(1),i_fram(2)
      do k=k1(i,j_fram),maxk
        tf_fram = tf_fram + ocn%u(2,i,j_fram,k)*dxv(j_fram)*dz(k) ! m3/s
        fw_fram = fw_fram + ocn%u(2,i,j_fram,k)*(1.0-0.5*(ocn%ts(i,j_fram+1,k,2) + ocn%ts(i,j_fram,k,2))/ocn%saln0) &
          * dxv(j_fram)*dz(k) ! m3/s
      enddo
    enddo
    tf_fram = tf_fram*1.e-6_wp  ! Sv
    fw_fram = fw_fram*1.e-6_wp  ! Sv

    ! 4. Denmark strait through flow (Sv)
    tf_denmark = 0.0
    fw_denmark = 0.0
    do i=i_denmark(1),i_denmark(2)
      do k=k1(i,j_denmark),maxk
        tf_denmark = tf_denmark + ocn%u(2,i,j_denmark,k)*dxv(j_denmark)*dz(k) ! m3/s
        fw_denmark = fw_denmark + ocn%u(2,i,j_denmark,k)*(1.0-0.5*(ocn%ts(i,j_denmark+1,k,2) + ocn%ts(i,j_denmark,k,2))/ocn%saln0) &
          * dxv(j_denmark)*dz(k) ! m3/s
      enddo
    enddo
    tf_denmark = tf_denmark*1.e-6_wp  ! Sv
    fw_denmark = fw_denmark*1.e-6_wp  ! Sv

    ! 5. Gibraltar Strait through flow (Sv)
    tf_medi = 0.0
    int_medi = .true.
    do j=loc_medi(2),maxj
       if (int_medi.and.k1(loc_medi(1),j).le.maxk) then
          do k=k1(loc_medi(1),j),maxk
             tf_medi = tf_medi + ocn%u(1,loc_medi(1),j,k)*dy*dz(k) ! m3/s
!	print *,'u',ocn%u(1,loc_medi(1),j,k),loc_medi(1),j,k
          enddo
       else
          int_medi = .false.
       endif
    enddo
    int_medi = .true.
    do j=loc_medi(2)-1,1,-1
       if (int_medi.and.k1(loc_medi(1),j).le.maxk) then
          do k=k1(loc_medi(1),j),maxk
             tf_medi = tf_medi + ocn%u(1,loc_medi(1),j,k)*dy*dz(k) ! m3/s
!        print *,'u',ocn%u(1,loc_medi(1),j,k),loc_medi(1),j,k
          enddo
       else
          int_medi = .false.
       endif
    enddo
    tf_medi = tf_medi*1.e-6_wp  ! Sv
!        print *,'loc_drake',loc_drake
!        print *,'tf_drake',tf_drake
!        print *,''

    ! 6. Indonesian Passage through flow (Sv)
    tf_indo = 0.0
    int_indo = .true.
    do j=loc_indo(2),maxj
       if (int_indo.and.k1(loc_indo(1),j).le.maxk) then
          do k=k1(loc_indo(1),j),maxk
             tf_indo = tf_indo + ocn%u(1,loc_indo(1),j,k)*dy*dz(k) ! m3/s
!	print *,'u',ocn%u(1,loc_indo(1),j,k),loc_indo(1),j,k
          enddo
       else
          int_indo = .false.
       endif
    enddo
    int_indo = .true.
    do j=loc_indo(2)-1,1,-1
       if (int_indo.and.k1(loc_indo(1),j).le.maxk) then
          do k=k1(loc_indo(1),j),maxk
             tf_indo = tf_indo + ocn%u(1,loc_indo(1),j,k)*dy*dz(k) ! m3/s
!        print *,'u',ocn%u(1,loc_indo(1),j,k),loc_indo(1),j,k
          enddo
       else
          int_indo = .false.
       endif
    enddo
    tf_indo = tf_indo*1.e-6_wp  ! Sv
!        print *,'loc_drake',loc_drake
!        print *,'tf_drake',tf_drake
!        print *,''


    ! 7. Angulhas Passage through flow (Sv)
    tf_agulhas = 0.0
    int_agulhas = .true.
    do j=loc_agulhas(2),maxj
       if (int_agulhas.and.k1(loc_agulhas(1),j).le.maxk) then
          do k=k1(loc_agulhas(1),j),maxk
             ! only if westward flow
             if (ocn%u(1,loc_agulhas(1),j,k).lt.0._wp) then
                tf_agulhas = tf_agulhas + ocn%u(1,loc_agulhas(1),j,k)*dy*dz(k) ! m3/s
                !print *,'u',ocn%u(1,loc_agulhas(1),j,k),loc_agulhas(1),j,k
             endif
          enddo
       else
          int_agulhas = .false.
       endif
    enddo
    int_agulhas = .true.
    do j=loc_agulhas(2)-1,1,-1
       if (int_agulhas.and.k1(loc_agulhas(1),j).le.maxk .and. lat(j).lt.-50._wp) then
          do k=k1(loc_agulhas(1),j),maxk
             ! only if westward flow
             if (ocn%u(1,loc_agulhas(1),j,k).lt.0._wp) then
                tf_agulhas = tf_agulhas + ocn%u(1,loc_agulhas(1),j,k)*dy*dz(k) ! m3/s
                !print *,'u',ocn%u(1,loc_agulhas(1),j,k),loc_agulhas(1),j,k
             endif
          enddo
       else
          int_agulhas = .false.
       endif
    enddo
    tf_agulhas = tf_agulhas*1.e-6_wp  ! Sv

    ! Atlantic meridional density gradient at 750 m depth between 52.5N and 32.5S
    ntot = count(basin_mask(:,js_drho1).eq.i_atlantic .and. k_drho.ge.k1(1:maxi,js_drho1))
    rho_s1 = sum(ocn%rho(:,js_drho1,k_drho), basin_mask(:,js_drho1).eq.i_atlantic .and. k_drho.ge.k1(1:maxi,js_drho1)) / ntot
    ntot = count(basin_mask(:,jn_drho1).eq.i_atlantic .and. k_drho.ge.k1(1:maxi,jn_drho1))
    rho_n1 = sum(ocn%rho(:,jn_drho1,k_drho), basin_mask(:,jn_drho1).eq.i_atlantic .and. k_drho.ge.k1(1:maxi,jn_drho1)) / ntot
    ! Atlantic meridional density gradient following Bonan 2022, difference
    ! between average 40-60N and the whole Atl basin south of 60N
    area_n = 0._wp
    area_n3 = 0._wp
    area_b = 0._wp
    rho_b2 = 0._wp
    rho_n2 = 0._wp
    rho_n3 = 0._wp
    rhoT_b2 = 0._wp
    rhoT_n2 = 0._wp
    rhoT_n3 = 0._wp
    rhoS_b2 = 0._wp
    rhoS_n2 = 0._wp
    rhoS_n3 = 0._wp
    alpha = (519._wp+122._wp*5._wp)*1.e-4_wp    ! kg/m3/K, assume value at 5 degC
    beta = 0.8_wp   ! kg/m3/psu
    do j=1,maxj
      do i=1,maxi
        if (basin_mask(i,j).eq.i_atlantic .and. k_drho.ge.k1(i,j)) then
          if (j.ge.js_drho2 .and. j.le.jn2_drho2) then
            area_ij = ocn%f_ocn(i,j)*ocn%grid%ocn_area(i,j)
            area_b = area_b + area_ij 
            rho_b2 = rho_b2 + ocn%rho(i,j,k_drho)*area_ij
            rhoT_b2 = rhoT_b2 - alpha*ocn%ts(i,j,k_drho,1)*area_ij
            rhoS_b2 = rhoS_b2 + beta*ocn%ts(i,j,k_drho,2)*area_ij
          endif
          if (j.ge.jn1_drho2 .and. j.le.jn2_drho2) then
            area_ij = ocn%f_ocn(i,j)*ocn%grid%ocn_area(i,j)
            area_n = area_n + area_ij 
            rho_n2 = rho_n2 + ocn%rho(i,j,k_drho)*area_ij
            rhoT_n2 = rhoT_n2 - alpha*ocn%ts(i,j,k_drho,1)*area_ij
            rhoS_n2 = rhoS_n2 + beta*ocn%ts(i,j,k_drho,2)*area_ij
          endif
          if (j.ge.(jn1_drho2+2) .and. j.le.(jn2_drho2+2)) then
            area_ij = ocn%f_ocn(i,j)*ocn%grid%ocn_area(i,j)
            area_n3 = area_n3 + area_ij 
            rho_n3 = rho_n3 + ocn%rho(i,j,k_drho)*area_ij
            rhoT_n3 = rhoT_n3 - alpha*ocn%ts(i,j,k_drho,1)*area_ij
            rhoS_n3 = rhoS_n3 + beta*ocn%ts(i,j,k_drho,2)*area_ij
          endif
        endif
      enddo
    enddo
    rho_b2 = rho_b2/area_b
    rho_n2 = rho_n2/area_n
    rho_n3 = rho_n3/area_n3
    rhoT_b2 = rhoT_b2/area_b
    rhoT_n2 = rhoT_n2/area_n
    rhoT_n3 = rhoT_n3/area_n3
    rhoS_b2 = rhoS_b2/area_b
    rhoS_n2 = rhoS_n2/area_n
    rhoS_n3 = rhoS_n3/area_n3

    ! set to zero at start of the year
    if( time_soy_ocn ) then
      ann_ts(y)%sst = 0._wp
      ann_ts(y)%sss = 0._wp
      ann_ts(y)%t   = 0._wp
      ann_ts(y)%s   = 0._wp
      ann_ts(y)%tdocn = 0._wp
      ann_ts(y)%tdocn_atl = 0._wp
      ann_ts(y)%tdocn_pac = 0._wp
      ann_ts(y)%tdocn_ind = 0._wp
      ann_ts(y)%tdocn_so = 0._wp
      ann_ts(y)%cons= 0._wp
      ann_ts(y)%age= 0._wp
      ann_ts(y)%dye= 0._wp
      ann_ts(y)%ohc   = 0._wp
      ann_ts(y)%ohc700   = 0._wp
      ann_ts(y)%ohc2000  = 0._wp
      ann_ts(y)%fw = 0._wp
      ann_ts(y)%fw_corr = 0._wp
      ann_ts(y)%saln0 = 0._wp
      ann_ts(y)%dvsf = 0._wp
      ann_ts(y)%fw_noise = 0._wp
      ann_ts(y)%p_e_sic = 0._wp
      ann_ts(y)%runoff = 0._wp
      ann_ts(y)%runoff_veg = 0._wp
      ann_ts(y)%runoff_ice = 0._wp
      ann_ts(y)%runoff_lake = 0._wp
      ann_ts(y)%calving = 0._wp
      ann_ts(y)%bmelt = 0._wp
      ann_ts(y)%icemelt = 0._wp
      ann_ts(y)%vsf = 0._wp
      ann_ts(y)%flx = 0._wp
      ann_ts(y)%drake = 0._wp
      ann_ts(y)%bering = 0._wp
      ann_ts(y)%davis = 0._wp
      ann_ts(y)%fram = 0._wp
      ann_ts(y)%denmark = 0._wp
      ann_ts(y)%medi = 0._wp
      ann_ts(y)%indo = 0._wp
      ann_ts(y)%agulhas = 0._wp
      ann_ts(y)%fov    = 0._wp
      ann_ts(y)%fovs   = 0._wp
      ann_ts(y)%fovn   = 0._wp
      ann_ts(y)%faz    = 0._wp
      ann_ts(y)%fazs   = 0._wp
      ann_ts(y)%fazn   = 0._wp
      ann_ts(y)%fw_lab   = 0._wp
      ann_ts(y)%buoy = 0._wp
      ann_ts(y)%buoyT = 0._wp
      ann_ts(y)%buoyS = 0._wp
      ann_ts(y)%buoy_N = 0._wp
      ann_ts(y)%buoyT_N = 0._wp
      ann_ts(y)%buoyS_N = 0._wp
      ann_ts(y)%buoy_tr = 0._wp
      ann_ts(y)%buoyT_tr = 0._wp
      ann_ts(y)%buoyS_tr = 0._wp
      ann_ts(y)%buoy_S = 0._wp
      ann_ts(y)%buoyT_S = 0._wp
      ann_ts(y)%buoyS_S = 0._wp
      ann_ts(y)%buoyT_NA(:) = 0._wp
      ann_ts(y)%buoyS_NA(:) = 0._wp
      ann_ts(y)%buoySw_NA(:) = 0._wp
      ann_ts(y)%buoySw_lab   = 0._wp
      ann_ts(y)%buoySi_NA(:) = 0._wp
      ann_ts(y)%buoy_NA(:) = 0._wp
      ann_ts(y)%fw_bering = 0._wp
      ann_ts(y)%fw_davis = 0._wp
      ann_ts(y)%fw_fram = 0._wp
      ann_ts(y)%fw_denmark = 0._wp
      ann_ts(y)%drhoatl1 = 0._wp
      ann_ts(y)%drhoatl2 = 0._wp
      ann_ts(y)%drhoTatl2 = 0._wp
      ann_ts(y)%drhoSatl2 = 0._wp
      ann_ts(y)%drhoatl3 = 0._wp
      ann_ts(y)%drhoTatl3 = 0._wp
      ann_ts(y)%drhoSatl3 = 0._wp
      ann_ts(y)%cfc11   = 0._wp
      ann_ts(y)%cfc12   = 0._wp
      ann_ts(y)%sl_steric = 0._wp
      ann_ts(y)%mld_atlN50 = 0._wp
      ann_ts(y)%mld_lab = 0._wp
      ann_ts(y)%mld_irm = 0._wp
      ann_ts(y)%mld_gin = 0._wp
      ann_ts(y)%mld_bkn = 0._wp
      ann_ts(y)%mld_wedd = 0._wp
      ann_ts(y)%mld_ross = 0._wp
      ann_ts(y)%mld_so = 0._wp
      ann_ts(y)%mldst_atlN50 = 0._wp
      ann_ts(y)%mldst_lab = 0._wp
      ann_ts(y)%mldst_irm = 0._wp
      ann_ts(y)%mldst_gin = 0._wp
      ann_ts(y)%mldst_bkn = 0._wp
      ann_ts(y)%mldst_wedd = 0._wp
      ann_ts(y)%mldst_ross = 0._wp
      ann_ts(y)%mldst_so = 0._wp
      ann_ts(y)%pe_atlN = 0._wp
      ann_ts(y)%pe_atlN50 = 0._wp
      ann_ts(y)%pe_lab = 0._wp
      ann_ts(y)%pe_irm = 0._wp
      ann_ts(y)%pe_gin = 0._wp
      ann_ts(y)%pe_bkn = 0._wp
      ann_ts(y)%pe_wedd = 0._wp
      ann_ts(y)%pe_ross = 0._wp
      ann_ts(y)%pe_so = 0._wp
      ann_ts(y)%buoy_lab = 0._wp
      ann_ts(y)%buoy_irm = 0._wp
      ann_ts(y)%buoy_gin = 0._wp
      ann_ts(y)%buoy_bkn = 0._wp
      ann_ts(y)%buoy_wedd = 0._wp
      ann_ts(y)%buoy_ross = 0._wp
      ann_ts(y)%buoy_so = 0._wp
      ann_ts(y)%buoy_soS60 = 0._wp
      ann_ts(y)%buoyT_soS60 = 0._wp
      ann_ts(y)%buoyS_soS60 = 0._wp
      ann_ts(y)%t_atlN50 = 0._wp
      ann_ts(y)%t_lab = 0._wp
      ann_ts(y)%t_irm = 0._wp
      ann_ts(y)%t_gin = 0._wp
      ann_ts(y)%t_bkn = 0._wp
      ann_ts(y)%t_wedd = 0._wp
      ann_ts(y)%t_ross = 0._wp
      ann_ts(y)%t_so = 0._wp
      ann_ts(y)%t_ibe = 0._wp
      ann_ts(y)%s_atlN50 = 0._wp
      ann_ts(y)%s_lab = 0._wp
      ann_ts(y)%s_irm = 0._wp
      ann_ts(y)%s_gin = 0._wp
      ann_ts(y)%s_bkn = 0._wp
      ann_ts(y)%s_wedd = 0._wp
      ann_ts(y)%s_ross = 0._wp
      ann_ts(y)%s_so = 0._wp

      ann_o%opsi  = 0._wp
      ann_o%opsia = 0._wp
      ann_o%opsip = 0._wp
      ann_o%opsii = 0._wp

      ann_o%hft   = 0._wp
      ann_o%hfp   = 0._wp
      ann_o%hfa   = 0._wp
      ann_o%hftz  = 0._wp
      ann_o%hfpz  = 0._wp
      ann_o%hfaz  = 0._wp
      ann_o%fwa   = 0._wp

    endif

    ! sum up and average over the year
    ann_ts(y)%sst    = ann_ts(y)%sst    + global_sst         * ann_avg ! °C
    ann_ts(y)%sss    = ann_ts(y)%sss    + global_sss * ann_avg ! psu
    ann_ts(y)%t      = ann_ts(y)%t      + global_tocn        * ann_avg ! °C
    ann_ts(y)%s      = ann_ts(y)%s      + global_socn* ann_avg ! psu
    ann_ts(y)%tdocn  = ann_ts(y)%tdocn  + tdocn        * ann_avg ! °C
    ann_ts(y)%tdocn_atl  = ann_ts(y)%tdocn_atl  + tdocn_atl        * ann_avg ! °C
    ann_ts(y)%tdocn_pac  = ann_ts(y)%tdocn_pac  + tdocn_pac        * ann_avg ! °C
    ann_ts(y)%tdocn_ind  = ann_ts(y)%tdocn_ind  + tdocn_ind        * ann_avg ! °C
    ann_ts(y)%tdocn_so   = ann_ts(y)%tdocn_so   + tdocn_so         * ann_avg ! °C
    if (cons_tracer) ann_ts(y)%cons = ann_ts(y)%cons + global_cons  * ann_avg
    if (age_tracer) ann_ts(y)%age = ann_ts(y)%age + global_age  * ann_avg
    if (dye_tracer) ann_ts(y)%dye = ann_ts(y)%dye + global_dye  * ann_avg
    ann_ts(y)%ohc = ann_ts(y)%ohc + ohc        * ann_avg ! J
    ann_ts(y)%ohc700 = ann_ts(y)%ohc700 + ohc700        * ann_avg ! J
    ann_ts(y)%ohc2000 = ann_ts(y)%ohc2000 + ohc2000        * ann_avg ! J
    ann_ts(y)%fw = ann_ts(y)%fw + fw/rho0*1.e-6_wp       * ann_avg ! Sv
    ann_ts(y)%fw_corr = ann_ts(y)%fw_corr + fw_corr/rho0*1.e-6_wp       * ann_avg ! Sv
    ann_ts(y)%saln0 = ann_ts(y)%saln0 + ocn%saln0 * ann_avg ! psu
    ann_ts(y)%dvsf = ann_ts(y)%dvsf + ocn%dvsf/ocn%saln0*ocn_area_tot*1.e-6_wp * ann_avg ! Sv
    ann_ts(y)%fw_noise = ann_ts(y)%fw_noise + fw_noise/rho0*1.e-6_wp       * ann_avg ! Sv
    ann_ts(y)%p_e_sic = ann_ts(y)%p_e_sic + p_e_sic/rho0*1.e-6_wp       * ann_avg ! Sv
    ann_ts(y)%runoff = ann_ts(y)%runoff + runoff/rho0*1.e-6_wp       * ann_avg ! Sv
    ann_ts(y)%runoff_veg = ann_ts(y)%runoff_veg + runoff_veg/rho0*1.e-6_wp       * ann_avg ! Sv
    ann_ts(y)%runoff_ice = ann_ts(y)%runoff_ice + runoff_ice/rho0*1.e-6_wp       * ann_avg ! Sv
    ann_ts(y)%runoff_lake = ann_ts(y)%runoff_lake + runoff_lake/rho0*1.e-6_wp       * ann_avg ! Sv
    ann_ts(y)%calving = ann_ts(y)%calving + calving/rho0*1.e-6_wp       * ann_avg ! Sv
    ann_ts(y)%bmelt = ann_ts(y)%bmelt + bmelt/rho0*1.e-6_wp       * ann_avg ! Sv
    ann_ts(y)%icemelt = ann_ts(y)%icemelt + icemelt/rho0*1.e-6_wp       * ann_avg ! Sv
    ann_ts(y)%vsf = ann_ts(y)%vsf + vsf/rho0*1.e-6_wp       * ann_avg ! Sv
    ann_ts(y)%flx = ann_ts(y)%flx + flx*1.e-15_wp   * ann_avg ! PW
    ann_ts(y)%drake = ann_ts(y)%drake + tf_drake    * ann_avg ! Sv
    ann_ts(y)%bering = ann_ts(y)%bering + tf_bering    * ann_avg ! Sv
    ann_ts(y)%davis = ann_ts(y)%davis + tf_davis    * ann_avg ! Sv
    ann_ts(y)%fram = ann_ts(y)%fram + tf_fram    * ann_avg ! Sv
    ann_ts(y)%denmark = ann_ts(y)%denmark + tf_denmark    * ann_avg ! Sv
    ann_ts(y)%medi = ann_ts(y)%medi + tf_medi    * ann_avg ! Sv
    ann_ts(y)%indo = ann_ts(y)%indo + tf_indo    * ann_avg ! Sv
    ann_ts(y)%agulhas = ann_ts(y)%agulhas + tf_agulhas    * ann_avg ! Sv
    ann_ts(y)%fov  = ann_ts(y)%fov + fov*1e-6 * ann_avg  ! Sv
    ann_ts(y)%fovs  = ann_ts(y)%fovs + fovs*1e-6 * ann_avg  ! Sv
    ann_ts(y)%fovn  = ann_ts(y)%fovn + fovn*1e-6 * ann_avg  ! Sv
    ann_ts(y)%faz  = ann_ts(y)%faz + faz*1e-6 * ann_avg  ! Sv
    ann_ts(y)%fazs  = ann_ts(y)%fazs + fazs*1e-6 * ann_avg  ! Sv
    ann_ts(y)%fazn  = ann_ts(y)%fazn + fazn*1e-6 * ann_avg  ! Sv
    ann_ts(y)%fw_lab    = ann_ts(y)%fw_lab + fw_lab * ann_avg ! Sv
    ann_ts(y)%buoy  = ann_ts(y)%buoy + buoy
    ann_ts(y)%buoyT  = ann_ts(y)%buoyT + buoyT
    ann_ts(y)%buoyS  = ann_ts(y)%buoyS + buoyS
    ann_ts(y)%buoy_N  = ann_ts(y)%buoy_N + buoy_N
    ann_ts(y)%buoyT_N  = ann_ts(y)%buoyT_N + buoyT_N
    ann_ts(y)%buoyS_N  = ann_ts(y)%buoyS_N + buoyS_N
    ann_ts(y)%buoy_tr  = ann_ts(y)%buoy_tr + buoy_tr
    ann_ts(y)%buoyT_tr  = ann_ts(y)%buoyT_tr + buoyT_tr
    ann_ts(y)%buoyS_tr  = ann_ts(y)%buoyS_tr + buoyS_tr
    ann_ts(y)%buoy_S  = ann_ts(y)%buoy_S + buoy_S
    ann_ts(y)%buoyT_S  = ann_ts(y)%buoyT_S + buoyT_S
    ann_ts(y)%buoyS_S  = ann_ts(y)%buoyS_S + buoyS_S
    ann_ts(y)%buoyT_NA(:)  = ann_ts(y)%buoyT_NA(:) + buoyT_NA(:)
    ann_ts(y)%buoyS_NA(:)  = ann_ts(y)%buoyS_NA(:) + buoyS_NA(:)
    ann_ts(y)%buoySw_NA(:)  = ann_ts(y)%buoySw_NA(:) + buoySw_NA(:)
    ann_ts(y)%buoySw_lab    = ann_ts(y)%buoySw_lab + buoySw_lab
    ann_ts(y)%buoySi_NA(:)  = ann_ts(y)%buoySi_NA(:) + buoySi_NA(:)
    ann_ts(y)%buoy_NA(:)  = ann_ts(y)%buoy_NA(:) + buoy_NA(:)
    ann_ts(y)%fw_bering = ann_ts(y)%fw_bering + fw_bering    * ann_avg ! Sv
    ann_ts(y)%fw_davis = ann_ts(y)%fw_davis + fw_davis    * ann_avg ! Sv
    ann_ts(y)%fw_fram = ann_ts(y)%fw_fram + fw_fram    * ann_avg ! Sv
    ann_ts(y)%fw_denmark = ann_ts(y)%fw_denmark + fw_denmark    * ann_avg ! Sv
    ann_ts(y)%drhoatl1 = ann_ts(y)%drhoatl1 + (rho_n1-rho_s1)    * ann_avg ! kg/m3
    ann_ts(y)%drhoatl2 = ann_ts(y)%drhoatl2 + (rho_n2-rho_b2)    * ann_avg ! kg/m3
    ann_ts(y)%drhoTatl2 = ann_ts(y)%drhoTatl2 + (rhoT_n2-rhoT_b2)    * ann_avg ! kg/m3
    ann_ts(y)%drhoSatl2 = ann_ts(y)%drhoSatl2 + (rhoS_n2-rhoS_b2)    * ann_avg ! kg/m3
    ann_ts(y)%drhoatl3 = ann_ts(y)%drhoatl3 + (rho_n3-rho_b2)    * ann_avg ! kg/m3
    ann_ts(y)%drhoTatl3 = ann_ts(y)%drhoTatl3 + (rhoT_n3-rhoT_b2)    * ann_avg ! kg/m3
    ann_ts(y)%drhoSatl3 = ann_ts(y)%drhoSatl3 + (rhoS_n3-rhoS_b2)    * ann_avg ! kg/m3
    if (l_cfc) then
      ann_ts(y)%cfc11      = ann_ts(y)%cfc11      + global_cfc11*1.e-6_wp        * ann_avg ! Mmol
      ann_ts(y)%cfc12      = ann_ts(y)%cfc12      + global_cfc12*1.e-6_wp        * ann_avg ! Mmol
    endif
    ann_ts(y)%sl_steric = ann_ts(y)%sl_steric + (sl_steric-sl_steric0)    * ann_avg ! m

    ann_o%opsi  = ann_o%opsi  + opsi(0:maxj,1:maxk)*1e-6  * ann_avg ! Sv
    ann_o%opsia = ann_o%opsia + opsia(0:maxj,1:maxk)*1e-6 * ann_avg ! Sv
    ann_o%opsip = ann_o%opsip + opsip(0:maxj,1:maxk)*1e-6 * ann_avg ! Sv
    ann_o%opsii = ann_o%opsii + opsii(0:maxj,1:maxk)*1e-6 * ann_avg ! Sv

    ann_o%hft(1:3,:) = ann_o%hft(1:3,:) + hft*1e-15 * ann_avg ! PW
    ann_o%hfp(1:3,:) = ann_o%hfp(1:3,:) + hfp*1e-15 * ann_avg ! PW
    ann_o%hfa(1:3,:) = ann_o%hfa(1:3,:) + hfa*1e-15 * ann_avg ! PW
    ann_o%hftz(:,:) = ann_o%hftz(:,:) + hftz*1e-15 * ann_avg ! PW/m
    ann_o%hfpz(:,:) = ann_o%hfpz(:,:) + hfpz*1e-15 * ann_avg ! PW/m
    ann_o%hfaz(:,:) = ann_o%hfaz(:,:) + hfaz*1e-15 * ann_avg ! PW/m
    ann_o%fwa(1:3,:) = ann_o%fwa(1:3,:) + fwa*1e-6 * ann_avg  ! Sv

    pe_atlN = 0._wp
    do i=i_atlN(1),i_atlN(2)
      do j=j_atlN(1),j_atlN(2)
        if (ocn%f_ocn(i,j).gt.0._wp) then
          pe_atlN = pe_atlN + ocn%conv_pe(i,j)*ocn%grid%ocn_area(i,j)   ! J/m2*m2=J
        endif
      enddo
    enddo

    mld_atlN50 = 0._wp
    mldst_atlN50 = 0._wp
    pe_atlN50 = 0._wp
    t_atlN50 = 0._wp
    s_atlN50 = 0._wp
    area_atlN50 = 0._wp
    do i=i_atlN50(1),i_atlN50(2)
      do j=j_atlN50(1),j_atlN50(2)
        if (ocn%f_ocn(i,j).gt.0._wp) then
          mldst = -zw(k1(i,j))
          rho_maxk = eos(ocn%ts(i,j,maxk,1),ocn%ts(i,j,maxk,2),0._wp)    ! in-situ density at the surface
          do k=maxk,k1(i,j),-1
            rho_k = eos(ocn%ts(i,j,k,1),ocn%ts(i,j,k,2),0._wp)   ! in-situ density
            if ((rho_k-rho_maxk).gt.0.125_wp) then
              mldst = -zw(k)
              exit
            endif
          enddo
          mldst_atlN50 = max(mldst_atlN50,mldst)
          mld_atlN50 = max(mld_atlN50,-ocn%mld(i,j))
          pe_atlN50 = pe_atlN50 + ocn%conv_pe(i,j)*ocn%grid%ocn_area(i,j)   ! J/m2*m2=J
          t_atlN50 = t_atlN50 + ocn%ts(i,j,maxk,1)*ocn%grid%ocn_area(i,j)
          s_atlN50 = s_atlN50 + ocn%ts(i,j,maxk,2)*ocn%grid%ocn_area(i,j)
          area_atlN50 = area_atlN50 + ocn%grid%ocn_area(i,j)
        endif
      enddo
    enddo
    t_atlN50 = t_atlN50/area_atlN50
    s_atlN50 = s_atlN50/area_atlN50


    mld_lab = 0._wp
    mldst_lab = 0._wp
    pe_lab = 0._wp
    buoy_lab = 0._wp
    t_lab = 0._wp
    s_lab = 0._wp
    area_lab = 0._wp
    do i=i_lab(1),i_lab(2)
      do j=j_lab(1),j_lab(2)
        if (ocn%f_ocn(i,j).gt.0._wp) then
          mldst = -zw(k1(i,j))
          rho_maxk = eos(ocn%ts(i,j,maxk,1),ocn%ts(i,j,maxk,2),0._wp)    ! in-situ density at the surface
          do k=maxk,k1(i,j),-1
            rho_k = eos(ocn%ts(i,j,k,1),ocn%ts(i,j,k,2),0._wp)   ! in-situ density
            if ((rho_k-rho_maxk).gt.0.125_wp) then
              mldst = -zw(k)
              exit
            endif
          enddo
          mldst_lab = max(mldst_lab,mldst)
          mld_lab = max(mld_lab,-ocn%mld(i,j))
          pe_lab = pe_lab + ocn%conv_pe(i,j)*ocn%grid%ocn_area(i,j)   ! J/m2*m2=J
          if (i_alphabeta.eq.1) then
            rho1 = eos(ocn%ts(i,j,maxk,1)+0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
            rho2 = eos(ocn%ts(i,j,maxk,1)-0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
            alpha = rho2-rho1
          else if (i_alphabeta.eq.2) then
            alpha = (519._wp+122._wp*ocn%ts(i,j,maxk,1))*1.e-4_wp    ! kg/m3/K
          endif
          beta = 0.8_wp   ! kg/m3/psu
          buoy_lab = buoy_lab + (ocn%flx(i,j)*alpha/cap_w+ocn%fw_corr(i,j)*beta*ocn%saln0) * g/rho0*ocn%grid%ocn_area(i,j)*dt  ! N 
          t_lab = t_lab + ocn%ts(i,j,maxk,1)*ocn%grid%ocn_area(i,j)
          s_lab = s_lab + ocn%ts(i,j,maxk,2)*ocn%grid%ocn_area(i,j)
          area_lab = area_lab + ocn%grid%ocn_area(i,j)
        endif
      enddo
    enddo
    t_lab = t_lab/area_lab
    s_lab = s_lab/area_lab

    mld_irm = 0._wp
    mldst_irm = 0._wp
    pe_irm = 0._wp
    buoy_irm = 0._wp
    t_irm = 0._wp
    s_irm = 0._wp
    area_irm = 0._wp
    do i=i_irm(1),i_irm(2)
      do j=j_irm(1),j_irm(2)
        if (ocn%f_ocn(i,j).gt.0._wp) then
          mldst = -zw(k1(i,j))
          rho_maxk = eos(ocn%ts(i,j,maxk,1),ocn%ts(i,j,maxk,2),0._wp)    ! in-situ density at the surface
          do k=maxk,k1(i,j),-1
            rho_k = eos(ocn%ts(i,j,k,1),ocn%ts(i,j,k,2),0._wp)   ! in-situ density
            if ((rho_k-rho_maxk).gt.0.125_wp) then
              mldst = -zw(k)
              exit
            endif
          enddo
          mldst_irm = max(mldst_irm,mldst)
          mld_irm = max(mld_irm,-ocn%mld(i,j))
          pe_irm = pe_irm + ocn%conv_pe(i,j)*ocn%grid%ocn_area(i,j)   ! J/m2*m2=J
          if (i_alphabeta.eq.1) then
            rho1 = eos(ocn%ts(i,j,maxk,1)+0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
            rho2 = eos(ocn%ts(i,j,maxk,1)-0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
            alpha = rho2-rho1
          else if (i_alphabeta.eq.2) then
            alpha = (519._wp+122._wp*ocn%ts(i,j,maxk,1))*1.e-4_wp    ! kg/m3/K
          endif
          beta = 0.8_wp   ! kg/m3/psu
          buoy_irm = buoy_irm + (ocn%flx(i,j)*alpha/cap_w+ocn%fw_corr(i,j)*beta*ocn%saln0) * g/rho0*ocn%grid%ocn_area(i,j)*dt  ! N 
          t_irm = t_irm + ocn%ts(i,j,maxk,1)*ocn%grid%ocn_area(i,j)
          s_irm = s_irm + ocn%ts(i,j,maxk,2)*ocn%grid%ocn_area(i,j)
          area_irm = area_irm + ocn%grid%ocn_area(i,j)
        endif
      enddo
    enddo
    t_irm = t_irm/area_irm
    s_irm = s_irm/area_irm

    mld_gin = 0._wp
    mldst_gin = 0._wp
    pe_gin = 0._wp
    buoy_gin = 0._wp
    t_gin = 0._wp
    s_gin = 0._wp
    area_gin = 0._wp
    do i=i_gin(1),i_gin(2)
      do j=j_gin(1),j_gin(2)
        if (ocn%f_ocn(i,j).gt.0._wp) then
          mldst = -zw(k1(i,j))
          rho_maxk = eos(ocn%ts(i,j,maxk,1),ocn%ts(i,j,maxk,2),0._wp)    ! in-situ density at the surface
          do k=maxk,k1(i,j),-1
            rho_k = eos(ocn%ts(i,j,k,1),ocn%ts(i,j,k,2),0._wp)   ! in-situ density
            if ((rho_k-rho_maxk).gt.0.125_wp) then
              mldst = -zw(k)
              exit
            endif
          enddo
          mldst_gin = max(mldst_gin,mldst)
          mld_gin = max(mld_gin,-ocn%mld(i,j))
          pe_gin = pe_gin + ocn%conv_pe(i,j)*ocn%grid%ocn_area(i,j)
          if (i_alphabeta.eq.1) then
            rho1 = eos(ocn%ts(i,j,maxk,1)+0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
            rho2 = eos(ocn%ts(i,j,maxk,1)-0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
            alpha = rho2-rho1
          else if (i_alphabeta.eq.2) then
            alpha = (519._wp+122._wp*ocn%ts(i,j,maxk,1))*1.e-4_wp    ! kg/m3/K
          endif
          beta = 0.8_wp   ! kg/m3/psu
          buoy_gin = buoy_gin + (ocn%flx(i,j)*alpha/cap_w+ocn%fw_corr(i,j)*beta*ocn%saln0) * g/rho0*ocn%grid%ocn_area(i,j)*dt  ! N 
          t_gin = t_gin + ocn%ts(i,j,maxk,1)*ocn%grid%ocn_area(i,j)
          s_gin = s_gin + ocn%ts(i,j,maxk,2)*ocn%grid%ocn_area(i,j)
          area_gin = area_gin + ocn%grid%ocn_area(i,j)
        endif
      enddo
    enddo
    t_gin = t_gin/area_gin
    s_gin = s_gin/area_gin

    mld_bkn = 0._wp
    mldst_bkn = 0._wp
    pe_bkn = 0._wp
    buoy_bkn = 0._wp
    t_bkn = 0._wp
    s_bkn = 0._wp
    area_bkn = 0._wp
    do i=i_bkn(1),i_bkn(2)
      do j=j_bkn(1),j_bkn(2)
        if (ocn%f_ocn(i,j).gt.0._wp) then
          mldst = -zw(k1(i,j))
          rho_maxk = eos(ocn%ts(i,j,maxk,1),ocn%ts(i,j,maxk,2),0._wp)    ! in-situ density at the surface
          do k=maxk,k1(i,j),-1
            rho_k = eos(ocn%ts(i,j,k,1),ocn%ts(i,j,k,2),0._wp)   ! in-situ density
            if ((rho_k-rho_maxk).gt.0.125_wp) then
              mldst = -zw(k)
              exit
            endif
          enddo
          mldst_bkn = max(mldst_bkn,mldst)
          mld_bkn = max(mld_bkn,-ocn%mld(i,j))
          pe_bkn = pe_bkn + ocn%conv_pe(i,j)*ocn%grid%ocn_area(i,j)
          if (i_alphabeta.eq.1) then
            rho1 = eos(ocn%ts(i,j,maxk,1)+0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
            rho2 = eos(ocn%ts(i,j,maxk,1)-0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
            alpha = rho2-rho1
          else if (i_alphabeta.eq.2) then
            alpha = (519._wp+122._wp*ocn%ts(i,j,maxk,1))*1.e-4_wp    ! kg/m3/K
          endif
          beta = 0.8_wp   ! kg/m3/psu
          buoy_bkn = buoy_bkn + (ocn%flx(i,j)*alpha/cap_w+ocn%fw_corr(i,j)*beta*ocn%saln0) * g/rho0*ocn%grid%ocn_area(i,j)*dt  ! N 
          t_bkn = t_bkn + ocn%ts(i,j,maxk,1)*ocn%grid%ocn_area(i,j)
          s_bkn = s_bkn + ocn%ts(i,j,maxk,2)*ocn%grid%ocn_area(i,j)
          area_bkn = area_bkn + ocn%grid%ocn_area(i,j)
        endif
      enddo
    enddo
    t_bkn = t_bkn/area_bkn
    s_bkn = s_bkn/area_bkn

    mld_wedd = 0._wp
    mldst_wedd = 0._wp
    pe_wedd = 0._wp
    buoy_wedd = 0._wp
    t_wedd = 0._wp
    s_wedd = 0._wp
    area_wedd = 0._wp
    do i=i_wedd(1),i_wedd(2)
      do j=j_wedd(1),j_wedd(2)
        if (ocn%f_ocn(i,j).gt.0._wp) then
          mldst = -zw(k1(i,j))
          rho_maxk = eos(ocn%ts(i,j,maxk,1),ocn%ts(i,j,maxk,2),0._wp)    ! in-situ density at the surface
          do k=maxk,k1(i,j),-1
            rho_k = eos(ocn%ts(i,j,k,1),ocn%ts(i,j,k,2),0._wp)   ! in-situ density
            if ((rho_k-rho_maxk).gt.0.125_wp) then
              mldst = -zw(k)
              exit
            endif
          enddo
          mldst_wedd = max(mldst_wedd,mldst)
          mld_wedd = max(mld_wedd,-ocn%mld(i,j))
          pe_wedd = pe_wedd + ocn%conv_pe(i,j)*ocn%grid%ocn_area(i,j)
          if (i_alphabeta.eq.1) then
            rho1 = eos(ocn%ts(i,j,maxk,1)+0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
            rho2 = eos(ocn%ts(i,j,maxk,1)-0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
            alpha = rho2-rho1
          else if (i_alphabeta.eq.2) then
            alpha = (519._wp+122._wp*ocn%ts(i,j,maxk,1))*1.e-4_wp    ! kg/m3/K
          endif
          beta = 0.8_wp   ! kg/m3/psu
          buoy_wedd = buoy_wedd + (ocn%flx(i,j)*alpha/cap_w+ocn%fw_corr(i,j)*beta*ocn%saln0) * g/rho0*ocn%grid%ocn_area(i,j)*dt  ! N 
          t_wedd = t_wedd + ocn%ts(i,j,maxk,1)*ocn%grid%ocn_area(i,j)
          s_wedd = s_wedd + ocn%ts(i,j,maxk,2)*ocn%grid%ocn_area(i,j)
          area_wedd = area_wedd + ocn%grid%ocn_area(i,j)
        endif
      enddo
    enddo
    t_wedd = t_wedd/area_wedd
    s_wedd = s_wedd/area_wedd

    mld_ross = 0._wp
    mldst_ross = 0._wp
    pe_ross = 0._wp
    buoy_ross = 0._wp
    t_ross = 0._wp
    s_ross = 0._wp
    area_ross = 0._wp
    do i=i_ross(1),i_ross(2)
      do j=j_ross(1),j_ross(2)
        if (ocn%f_ocn(i,j).gt.0._wp) then
          mldst = -zw(k1(i,j))
          rho_maxk = eos(ocn%ts(i,j,maxk,1),ocn%ts(i,j,maxk,2),0._wp)    ! in-situ density at the surface
          do k=maxk,k1(i,j),-1
            rho_k = eos(ocn%ts(i,j,k,1),ocn%ts(i,j,k,2),0._wp)   ! in-situ density
            if ((rho_k-rho_maxk).gt.0.125_wp) then
              mldst = -zw(k)
              exit
            endif
          enddo
          mldst_ross = max(mldst_ross,mldst)
          mld_ross = max(mld_ross,-ocn%mld(i,j))
          pe_ross = pe_ross + ocn%conv_pe(i,j)*ocn%grid%ocn_area(i,j)
          if (i_alphabeta.eq.1) then
            rho1 = eos(ocn%ts(i,j,maxk,1)+0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
            rho2 = eos(ocn%ts(i,j,maxk,1)-0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
            alpha = rho2-rho1
          else if (i_alphabeta.eq.2) then
            alpha = (519._wp+122._wp*ocn%ts(i,j,maxk,1))*1.e-4_wp    ! kg/m3/K
          endif
          beta = 0.8_wp   ! kg/m3/psu
          buoy_ross = buoy_ross + (ocn%flx(i,j)*alpha/cap_w+ocn%fw_corr(i,j)*beta*ocn%saln0) * g/rho0*ocn%grid%ocn_area(i,j)*dt  ! N 
          t_ross = t_ross + ocn%ts(i,j,maxk,1)*ocn%grid%ocn_area(i,j)
          s_ross = s_ross + ocn%ts(i,j,maxk,2)*ocn%grid%ocn_area(i,j)
          area_ross = area_ross + ocn%grid%ocn_area(i,j)
        endif
      enddo
    enddo
    t_ross = t_ross/area_ross
    s_ross = s_ross/area_ross

    mld_so = 0._wp
    mldst_so = 0._wp
    pe_so = 0._wp
    buoy_so = 0._wp
    t_sos = 0._wp
    s_sos = 0._wp
    area_so = 0._wp
    do i=i_so(1),i_so(2)
      do j=j_so(1),j_so(2)
        if (ocn%f_ocn(i,j).gt.0._wp) then
          mldst = -zw(k1(i,j))
          rho_maxk = eos(ocn%ts(i,j,maxk,1),ocn%ts(i,j,maxk,2),0._wp)    ! in-situ density at the surface
          do k=maxk,k1(i,j),-1
            rho_k = eos(ocn%ts(i,j,k,1),ocn%ts(i,j,k,2),0._wp)   ! in-situ density
            if ((rho_k-rho_maxk).gt.0.125_wp) then
              mldst = -zw(k)
              exit
            endif
          enddo
          mldst_so = max(mldst_so,mldst)
          mld_so = max(mld_so,-ocn%mld(i,j))
          pe_so = pe_so + ocn%conv_pe(i,j)*ocn%grid%ocn_area(i,j)
          if (i_alphabeta.eq.1) then
            rho1 = eos(ocn%ts(i,j,maxk,1)+0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
            rho2 = eos(ocn%ts(i,j,maxk,1)-0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
            alpha = rho2-rho1
          else if (i_alphabeta.eq.2) then
            alpha = (519._wp+122._wp*ocn%ts(i,j,maxk,1))*1.e-4_wp    ! kg/m3/K
          endif
          beta = 0.8_wp   ! kg/m3/psu
          buoy_so = buoy_so + (ocn%flx(i,j)*alpha/cap_w+ocn%fw_corr(i,j)*beta*ocn%saln0) * g/rho0*ocn%grid%ocn_area(i,j)*dt  ! N 

          t_sos = t_sos + ocn%ts(i,j,maxk,1)*ocn%grid%ocn_area(i,j)
          s_sos = s_sos + ocn%ts(i,j,maxk,2)*ocn%grid%ocn_area(i,j)
          area_so = area_so + ocn%grid%ocn_area(i,j)
        endif
      enddo
    enddo
    t_sos = t_sos/area_so
    s_sos = s_sos/area_so

    buoy_soS60 = 0._wp
    buoyT_soS60 = 0._wp
    buoyS_soS60 = 0._wp
    do i=i_soS60(1),i_soS60(2)
      do j=j_soS60(1),j_soS60(2)
        if (ocn%f_ocn(i,j).gt.0._wp) then
          if (i_alphabeta.eq.1) then
            rho1 = eos(ocn%ts(i,j,maxk,1)+0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
            rho2 = eos(ocn%ts(i,j,maxk,1)-0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
            alpha = rho2-rho1
          else if (i_alphabeta.eq.2) then
            alpha = (519._wp+122._wp*ocn%ts(i,j,maxk,1))*1.e-4_wp    ! kg/m3/K
          endif
          beta = 0.8_wp   ! kg/m3/psu
          buoy_soS60  = buoy_soS60  + (ocn%flx(i,j)*alpha/cap_w+ocn%fw_corr(i,j)*beta*ocn%saln0) * g/rho0*ocn%grid%ocn_area(i,j)*dt  ! N 
          buoyT_soS60 = buoyT_soS60 + ocn%flx(i,j)*alpha/cap_w * g/rho0*ocn%grid%ocn_area(i,j)*dt  ! N 
          buoyS_soS60 = buoyS_soS60 + ocn%fw_corr(i,j)*beta*ocn%saln0 * g/rho0*ocn%grid%ocn_area(i,j)*dt  ! N 
        endif
      enddo
    enddo

    ann_ts(y)%mld_atlN50= max(ann_ts(y)%mld_atlN50,mld_atlN50)
    ann_ts(y)%mld_lab   = max(ann_ts(y)%mld_lab,mld_lab)
    ann_ts(y)%mld_irm   = max(ann_ts(y)%mld_irm,mld_irm)
    ann_ts(y)%mld_gin  = max(ann_ts(y)%mld_gin,mld_gin)
    ann_ts(y)%mld_bkn  = max(ann_ts(y)%mld_bkn,mld_bkn)
    ann_ts(y)%mld_wedd = max(ann_ts(y)%mld_wedd,mld_wedd)
    ann_ts(y)%mld_ross = max(ann_ts(y)%mld_ross,mld_ross)
    ann_ts(y)%mld_so   = max(ann_ts(y)%mld_so,mld_so)
    ann_ts(y)%mldst_atlN50   = max(ann_ts(y)%mldst_atlN50,mldst_atlN50)
    ann_ts(y)%mldst_lab   = max(ann_ts(y)%mldst_lab,mldst_lab)
    ann_ts(y)%mldst_irm   = max(ann_ts(y)%mldst_irm,mldst_irm)
    ann_ts(y)%mldst_gin  = max(ann_ts(y)%mldst_gin,mldst_gin)
    ann_ts(y)%mldst_bkn  = max(ann_ts(y)%mldst_bkn,mldst_bkn)
    ann_ts(y)%mldst_wedd = max(ann_ts(y)%mldst_wedd,mldst_wedd)
    ann_ts(y)%mldst_ross = max(ann_ts(y)%mldst_ross,mldst_ross)
    ann_ts(y)%mldst_so   = max(ann_ts(y)%mldst_so,mldst_so)
    ann_ts(y)%pe_atlN   = ann_ts(y)%pe_atlN   + pe_atlN    
    ann_ts(y)%pe_atlN50   = ann_ts(y)%pe_atlN50   + pe_atlN50    
    ann_ts(y)%pe_lab   = ann_ts(y)%pe_lab   + pe_lab    
    ann_ts(y)%pe_irm   = ann_ts(y)%pe_irm   + pe_irm
    ann_ts(y)%pe_gin  = ann_ts(y)%pe_gin  + pe_gin  
    ann_ts(y)%pe_bkn  = ann_ts(y)%pe_bkn  + pe_bkn  
    ann_ts(y)%pe_wedd = ann_ts(y)%pe_wedd + pe_wedd 
    ann_ts(y)%pe_ross = ann_ts(y)%pe_ross + pe_ross 
    ann_ts(y)%pe_so   = ann_ts(y)%pe_so   + pe_so
    ann_ts(y)%buoy_lab   = ann_ts(y)%buoy_lab   + buoy_lab    
    ann_ts(y)%buoy_irm   = ann_ts(y)%buoy_irm   + buoy_irm
    ann_ts(y)%buoy_gin  = ann_ts(y)%buoy_gin  + buoy_gin  
    ann_ts(y)%buoy_bkn  = ann_ts(y)%buoy_bkn  + buoy_bkn  
    ann_ts(y)%buoy_wedd = ann_ts(y)%buoy_wedd + buoy_wedd 
    ann_ts(y)%buoy_ross = ann_ts(y)%buoy_ross + buoy_ross 
    ann_ts(y)%buoy_so   = ann_ts(y)%buoy_so   + buoy_so
    ann_ts(y)%buoy_soS60   = ann_ts(y)%buoy_soS60   + buoy_soS60
    ann_ts(y)%buoyT_soS60   = ann_ts(y)%buoyT_soS60   + buoyT_soS60
    ann_ts(y)%buoyS_soS60   = ann_ts(y)%buoyS_soS60   + buoyS_soS60
    ann_ts(y)%t_atlN50   = ann_ts(y)%t_atlN50   + t_atlN50   * ann_avg 
    ann_ts(y)%t_lab   = ann_ts(y)%t_lab   + t_lab   * ann_avg 
    ann_ts(y)%t_irm   = ann_ts(y)%t_irm   + t_irm   * ann_avg 
    ann_ts(y)%t_gin  = ann_ts(y)%t_gin  + t_gin  * ann_avg
    ann_ts(y)%t_bkn  = ann_ts(y)%t_bkn  + t_bkn  * ann_avg
    ann_ts(y)%t_wedd = ann_ts(y)%t_wedd + t_wedd * ann_avg
    ann_ts(y)%t_ross = ann_ts(y)%t_ross + t_ross * ann_avg
    ann_ts(y)%t_so   = ann_ts(y)%t_so   + t_sos  * ann_avg
    ann_ts(y)%t_ibe  = ann_ts(y)%t_ibe  + 0.5_wp*(ocn%ts(i_ibe,j_ibe,maxk,1)+ocn%ts(i_ibe-1,j_ibe,maxk,1))  * ann_avg
    ann_ts(y)%s_atlN50   = ann_ts(y)%s_atlN50   + s_atlN50   * ann_avg 
    ann_ts(y)%s_lab   = ann_ts(y)%s_lab   + s_lab   * ann_avg 
    ann_ts(y)%s_irm   = ann_ts(y)%s_irm   + s_irm   * ann_avg 
    ann_ts(y)%s_gin  = ann_ts(y)%s_gin  + s_gin  * ann_avg
    ann_ts(y)%s_bkn  = ann_ts(y)%s_bkn  + s_bkn  * ann_avg
    ann_ts(y)%s_wedd = ann_ts(y)%s_wedd + s_wedd * ann_avg
    ann_ts(y)%s_ross = ann_ts(y)%s_ross + s_ross * ann_avg
    ann_ts(y)%s_so   = ann_ts(y)%s_so   + s_sos  * ann_avg


    if( time_eoy_ocn ) then

       ! Pacific
       ann_ts(y)%ominp = 0
       ann_ts(y)%omaxp = 0
       do j=jps,maxj-1
          do k=1,maxk-1
             if (ann_o%opsip(j,k).lt.ann_ts(y)%ominp.and.k.le.k_over) ann_ts(y)%ominp = ann_o%opsip(j,k)
             if (ann_o%opsip(j,k).gt.ann_ts(y)%omaxp.and.k.le.k_over) ann_ts(y)%omaxp = ann_o%opsip(j,k)
          enddo
       enddo
   
       ! Atlantic
       ann_ts(y)%omina = 0.0
       ann_ts(y)%omaxa = 0.0
       ann_ts(y)%oinfa = 0.0
       do j=jas,maxj-1
          do k=1,maxk-1
             if (ann_o%opsia(j,k).lt.ann_ts(y)%omina.and.k.le.k_over) ann_ts(y)%omina = ann_o%opsia(j,k)
             if (ann_o%opsia(j,k).gt.ann_ts(y)%omaxa.and.k.le.k_over) then
                ann_ts(y)%omaxa = ann_o%opsia(j,k)
                ! position of Atlantic maximum overturning now recorded
                iposa(1) = j
                iposa(2) = k
             endif
             ! southern Atlantic inflow
             if (j.eq.jas.and.ann_o%opsia(j,k).gt.ann_ts(y)%oinfa.and.k.le.k_over) ann_ts(y)%oinfa = ann_o%opsia(j,k)
          enddo
       enddo
       ann_ts(y)%amoc26N = 0._wp
       do k=1,maxk-1
         if (k.le.k_over) then
           ann_ts(y)%amoc26N = max(ann_ts(y)%amoc26N,ann_o%opsia(J26N,k))
         endif
       enddo
       ! save in ocn variable
       ocn%amoc = ann_ts(y)%omaxa
   
       ! Southern ocean
       ann_ts(y)%omins = 0.0
       ann_ts(y)%omaxs = 0.0
       !do j=1,jsf
       do j=1,8
          do k=1,maxk-1
             if (ann_o%opsi(j,k).lt.ann_ts(y)%omins) ann_ts(y)%omins = ann_o%opsi(j,k)
             if (ann_o%opsi(j,k).gt.ann_ts(y)%omaxs) ann_ts(y)%omaxs = ann_o%opsi(j,k)
          enddo
       enddo

       ! maximum northward Atlantic heat transport
       ann_ts(y)%hmaxa = maxval(ann_o%hfa(1,:))
       ann_ts(y)%h55a  = ann_o%hfa(1,29)    ! at 55N
       ! maximum northward Pacific heat transport
       ann_ts(y)%hmaxp = maxval(ann_o%hfp(1,:))
       ! maximum southward total heat transport
       ann_ts(y)%hmaxs = -minval(ann_o%hft(1,:))
       ann_ts(y)%h60s  = -ann_o%hft(1,6)

       ! maximum northward Atlantic freshwater transport
       ann_ts(y)%fmaxa = maxval(ann_o%fwa)

       ! freshwater hosing 
       ann_ts(y)%hosing = ocn%hosing

       ! noise
       ann_ts(y)%noise_fw = ocn%noise_fw*sec_day    ! kg/m2/day
       ann_ts(y)%noise_flx = ocn%noise_flx          ! W/m2

       ! ocean surface area and volume
       ann_ts(y)%area = ocn_area_tot * 1.e-12_wp ! mln km2
       ann_ts(y)%vol  = ocn_vol_tot * 1.e-15_wp ! mln km3

       ! number of active ocean grid cells
       ann_ts(y)%ncells = 0
       do j=1,maxj
         do i=1,maxi
           if (ocn%f_ocn(i,j).gt.0._wp) then
             ann_ts(y)%ncells = ann_ts(y)%ncells + 1
           endif
         enddo
       enddo

       ! total continental shelf ocean area
       shelf = 0._wp
       do i=1,maxi
         do j=1,maxj
           if (k1(i,j).ge.k1_shelf .and. k1(i,j).le.maxk) then
             shelf = shelf + ocn%f_ocn(i,j)*ocn%grid%ocn_area(i,j) ! m2
           endif
         enddo
       enddo
       ann_ts(y)%shelf = shelf * 1.e-12_wp ! mln km2


       ! write to standard output
       if (mod(year,10).eq.1) then
         print '(a7,a9,14a7)','ocn','year','SST','SSS','TVOL','SVOL','AMOC','PMOC','AABW','hAMOC','FWFnet','FLXnet','Drake','Indo','Bering','Ashelf'
       endif

       print '(a7,i9,8F7.1,2F7.2,4F7.1)', &
         'ocn',year_now,ann_ts(y)%sst,ann_ts(y)%sss,ann_ts(y)%t,ann_ts(y)%s,ann_ts(y)%omaxa,ann_ts(y)%omaxp,-ann_ts(y)%omins, &
         ann_ts(y)%hmaxa,ann_ts(y)%fw(1),ann_ts(y)%flx(1)*1.e15_wp/ocn%grid%ocn_area_tot,ann_ts(y)%drake,ann_ts(y)%indo,ann_ts(y)%bering,ann_ts(y)%shelf

       ! write to netcdf file 
       if (time_out_ts_clim) then
         call ts_nc_write(trim(out_dir)//"/ocn_ts.nc",ann_ts(1:y),year_clim-y+1,y)
       endif

     endif


     ! spatially explicit output
     if ( time_out_ocn ) then

       if( time_soy_ocn ) then
         do m=1,nmon_year
           do i=1,maxi
             do j=1,maxj
               do k=1,maxk
                 kr = maxk-k+1
                 if (mask_c(i,j,kr).eq.1) then
                   mon_o(m)%t   (i,j,k)  = 0._wp
                   mon_o(m)%s   (i,j,k)  = 0._wp
                   mon_o(m)%dt_dt_adv(i,j,k)  = 0._wp
                   mon_o(m)%dt_dt_diff(i,j,k) = 0._wp
                   mon_o(m)%ds_dt_adv(i,j,k)  = 0._wp
                   mon_o(m)%ds_dt_diff(i,j,k) = 0._wp
                   mon_o(m)%rho (i,j,k)  = 0._wp
                   mon_o(m)%rho2(i,j,k)  = 0._wp
                   mon_o(m)%diffdia(i,j,k) = 0._wp
                   mon_o(m)%drho_dx(i,j,k) = 0._wp
                   mon_o(m)%drho_dy(i,j,k) = 0._wp
                   mon_o(m)%drho_dz(i,j,k) = 0._wp
                   mon_o(m)%Ri(i,j,k) = 0._wp
                 else
                   mon_o(m)%t   (i,j,k)  = missing_value 
                   mon_o(m)%s   (i,j,k)  = missing_value 
                   mon_o(m)%dt_dt_adv(i,j,k)  = missing_value 
                   mon_o(m)%dt_dt_diff(i,j,k) = missing_value 
                   mon_o(m)%ds_dt_adv(i,j,k)  = missing_value 
                   mon_o(m)%ds_dt_diff(i,j,k) = missing_value 
                   mon_o(m)%rho (i,j,k)  = missing_value 
                   mon_o(m)%rho2(i,j,k)  = missing_value 
                   mon_o(m)%diffdia(i,j,k) = missing_value 
                   mon_o(m)%drho_dx(i,j,k) = missing_value 
                   mon_o(m)%drho_dy(i,j,k) = missing_value 
                   mon_o(m)%drho_dz(i,j,k) = missing_value 
                   mon_o(m)%Ri(i,j,k) = missing_value 
                 endif
              enddo
            enddo
          enddo
          mon_o(m)%ub    = 0._wp
          mon_o(m)%vb    = 0._wp
          mon_o(m)%u     = 0._wp
          mon_o(m)%v     = 0._wp
          mon_o(m)%w     = 0._wp
          mon_o(m)%fdx   = 0._wp
          mon_o(m)%fdy   = 0._wp
          mon_o(m)%fdz   = 0._wp
          mon_o(m)%sst   = 0._wp
          mon_o(m)%tbot  = 0._wp
          mon_o(m)%t_atl = 0._wp
          mon_o(m)%t_pac = 0._wp
          mon_o(m)%t_ind = 0._wp
          mon_o(m)%t_so  = 0._wp
          mon_o(m)%sss   = 0._wp
          mon_o(m)%sbot  = 0._wp
          mon_o(m)%s_atl = 0._wp
          mon_o(m)%s_pac = 0._wp
          mon_o(m)%s_ind = 0._wp
          mon_o(m)%s_so  = 0._wp
          mon_o(m)%rho_atl = 0._wp
          mon_o(m)%rho_pac = 0._wp
          mon_o(m)%rho_ind = 0._wp
          mon_o(m)%rho_so  = 0._wp
          mon_o(m)%drho_0_1000 = 0._wp
          mon_o(m)%drho_0_3000 = 0._wp
          mon_o(m)%psi   = 0._wp
          mon_o(m)%opsi  = 0._wp
          mon_o(m)%opsia = 0._wp
          mon_o(m)%opsip = 0._wp
          mon_o(m)%opsii = 0._wp
          mon_o(m)%hft   = 0._wp
          mon_o(m)%hfp   = 0._wp
          mon_o(m)%hfa   = 0._wp
          mon_o(m)%hftz  = 0._wp
          mon_o(m)%hfpz  = 0._wp
          mon_o(m)%hfaz  = 0._wp
          mon_o(m)%fwt   = 0._wp
          mon_o(m)%fwp   = 0._wp
          mon_o(m)%fwa   = 0._wp
          mon_o(m)%fwtz  = 0._wp
          mon_o(m)%fwpz  = 0._wp
          mon_o(m)%fwaz  = 0._wp
          mon_o(m)%fayti = 0._wp
          mon_o(m)%fdyti = 0._wp
          mon_o(m)%faysi = 0._wp
          mon_o(m)%fdysi = 0._wp
          mon_o(m)%taux  = 0._wp
          mon_o(m)%tauy  = 0._wp
          mon_o(m)%ssh = 0._wp
          where (ocn%f_ocn.gt.0._wp) 
            mon_o(m)%dt_dt_flxsur = 0._wp
            mon_o(m)%ds_dt_flxsur = 0._wp
            mon_o(m)%flx   = 0._wp
            mon_o(m)%fw    = 0._wp
            mon_o(m)%vsf   = 0._wp
            mon_o(m)%p_e_sic   = 0._wp
            mon_o(m)%runoff= 0._wp
            mon_o(m)%runoffSv= 0._wp
            mon_o(m)%runoffSv_ice= 0._wp
            mon_o(m)%calving= 0._wp
            mon_o(m)%bmelt= 0._wp
            mon_o(m)%fw_hosing = 0._wp
            mon_o(m)%fw_flux_adj = 0._wp
            mon_o(m)%fw_noise = 0._wp
            mon_o(m)%flx_noise = 0._wp
            mon_o(m)%buoy  = 0._wp
            mon_o(m)%buoyS = 0._wp
            mon_o(m)%buoyT = 0._wp
            mon_o(m)%nconv = 0._wp
            mon_o(m)%kven  = 0._wp
            mon_o(m)%dconv = 0._wp
            mon_o(m)%dven  = 0._wp
            mon_o(m)%mld   = 0._wp
            mon_o(m)%mldst = 0._wp
            mon_o(m)%ke_tau = 0._wp
          elsewhere
            mon_o(m)%dt_dt_flxsur = missing_value 
            mon_o(m)%ds_dt_flxsur = missing_value 
            mon_o(m)%flx   = missing_value
            mon_o(m)%fw    = missing_value
            mon_o(m)%vsf   = missing_value
            mon_o(m)%p_e_sic   = missing_value
            mon_o(m)%runoff= missing_value
            mon_o(m)%runoffSv= missing_value
            mon_o(m)%runoffSv_ice= missing_value
            mon_o(m)%calving= missing_value
            mon_o(m)%bmelt= missing_value
            mon_o(m)%fw_hosing = missing_value
            mon_o(m)%fw_flux_adj = missing_value
            mon_o(m)%fw_noise = missing_value
            mon_o(m)%flx_noise = missing_value
            mon_o(m)%buoy  = missing_value
            mon_o(m)%buoyS = missing_value
            mon_o(m)%buoyT = missing_value
            mon_o(m)%nconv = missing_value
            mon_o(m)%kven  = missing_value
            mon_o(m)%dconv = missing_value
            mon_o(m)%dven  = missing_value
            mon_o(m)%mld   = missing_value
            mon_o(m)%mldst = missing_value
            mon_o(m)%ke_tau   = missing_value
          endwhere
          do n=1,n_isles
            mon_o(m)%ubisl(:,:,n) = 0._wp
            mon_o(m)%vbisl(:,:,n) = 0._wp
          enddo
          do n=n_isles+1,maxisles
            mon_o(m)%ubisl(:,:,n) = missing_value   
            mon_o(m)%vbisl(:,:,n) = missing_value
          enddo

        enddo
      endif

       ! compute in-situ density
       allocate(rho(maxi,maxj,maxk))
       do i=1,maxi
         do j=1,maxj
           do k=1,maxk
             if (mask_c(i,j,k).eq.1) then
              rho(i,j,k) = eos(ocn%ts(i,j,k,1),ocn%ts(i,j,k,2),0._wp)
            endif
           enddo
         enddo
       enddo

    ! density difference between surface and different depths
    do j=1,maxj
      do i=1,maxi
        if (mask_c(i,j,k1_1000).eq.1) then
          mon_o(mon)%drho_0_1000(i,j) = mon_o(mon)%drho_0_1000(i,j) + (rho(i,j,maxk)-rho(i,j,k1_1000)) * mon_avg
        else
          mon_o(mon)%drho_0_1000(i,j) = missing_value
        endif
        if (mask_c(i,j,k1_3000).eq.1) then
          mon_o(mon)%drho_0_3000(i,j) = mon_o(mon)%drho_0_3000(i,j) + (rho(i,j,maxk)-rho(i,j,k1_3000)) * mon_avg
        else
          mon_o(mon)%drho_0_3000(i,j) = missing_value
        endif
      enddo
    enddo

    ! zonal basin averages
    do j=1,maxj
       do k=1,maxk
          ntot = count(basin_mask2(:,j).eq.i_atlantic .and. k.ge.k1(1:maxi,j))
          if (ntot.gt.0) then
            t_atl(j,k) = sum(ocn%ts(:,j,k,1), basin_mask2(:,j).eq.i_atlantic .and. k.ge.k1(1:maxi,j)) / ntot
            s_atl(j,k) = sum(ocn%ts(:,j,k,2), basin_mask2(:,j).eq.i_atlantic .and. k.ge.k1(1:maxi,j)) / ntot
            rho_atl(j,k) = sum(rho(:,j,k), basin_mask2(:,j).eq.i_atlantic .and. k.ge.k1(1:maxi,j)) / ntot
          else
            t_atl(j,k) = missing_value
            s_atl(j,k) = missing_value
            rho_atl(j,k) = missing_value
          endif
          ntot = count(basin_mask2(:,j).eq.i_pacific .and. k.ge.k1(1:maxi,j))
          if (ntot.gt.0) then
            t_pac(j,k) = sum(ocn%ts(:,j,k,1), basin_mask2(:,j).eq.i_pacific .and. k.ge.k1(1:maxi,j)) / ntot
            s_pac(j,k) = sum(ocn%ts(:,j,k,2), basin_mask2(:,j).eq.i_pacific .and. k.ge.k1(1:maxi,j)) / ntot
            rho_pac(j,k) = sum(rho(:,j,k), basin_mask2(:,j).eq.i_pacific .and. k.ge.k1(1:maxi,j)) / ntot
          else
            t_pac(j,k) = missing_value
            s_pac(j,k) = missing_value
            rho_pac(j,k) = missing_value
          endif
          ntot = count(basin_mask2(:,j).eq.i_indian .and. k.ge.k1(1:maxi,j))
          if (ntot.gt.0) then
            t_ind(j,k) = sum(ocn%ts(:,j,k,1), basin_mask2(:,j).eq.i_indian .and. k.ge.k1(1:maxi,j)) / ntot
            s_ind(j,k) = sum(ocn%ts(:,j,k,2), basin_mask2(:,j).eq.i_indian .and. k.ge.k1(1:maxi,j)) / ntot
            rho_ind(j,k) = sum(rho(:,j,k), basin_mask2(:,j).eq.i_indian .and. k.ge.k1(1:maxi,j)) / ntot
          else
            t_ind(j,k) = missing_value
            s_ind(j,k) = missing_value
            rho_ind(j,k) = missing_value
          endif
          ntot = count(basin_mask(:,j).eq.i_southern .and. k.ge.k1(1:maxi,j))
          if (ntot.gt.0) then
            t_so(j,k)  = sum(ocn%ts(:,j,k,1), basin_mask(:,j).eq.i_southern .and. k.ge.k1(1:maxi,j)) / ntot
            s_so(j,k)  = sum(ocn%ts(:,j,k,2), basin_mask(:,j).eq.i_southern .and. k.ge.k1(1:maxi,j)) / ntot
            rho_so(j,k)  = sum(rho(:,j,k), basin_mask(:,j).eq.i_southern .and. k.ge.k1(1:maxi,j)) / ntot
          else
            t_so(j,k)  = missing_value
            s_so(j,k)  = missing_value
            rho_so(j,k)  = missing_value
          endif
       enddo
    enddo

    do j=1,maxj
      do k=1,maxk
        kr = maxk-k+1
        do i=1,maxi
          if (mask_c(i,j,kr).eq.1) then
            mon_o(mon)%t      (i,j,k) = mon_o(mon)%t      (i,j,k)   + ocn%ts  (i,j,kr,1) * mon_avg ! °C
            mon_o(mon)%s      (i,j,k) = mon_o(mon)%s      (i,j,k)   + ocn%ts  (i,j,kr,2) * mon_avg ! psu
            if (kr.eq.maxk) then
              mon_o(mon)%dt_dt_adv(i,j,k) = mon_o(mon)%dt_dt_adv(i,j,k) + (ocn%dts_dt_adv(i,j,kr,1)+ocn%flx_sur(i,j,1)/dz(kr)) * sec_year * mon_avg ! °C/year
              mon_o(mon)%ds_dt_adv(i,j,k) = mon_o(mon)%ds_dt_adv(i,j,k) + (ocn%dts_dt_adv(i,j,kr,2)+ocn%flx_sur(i,j,2)/dz(kr)) * sec_year * mon_avg ! psu/year
            else
              mon_o(mon)%dt_dt_adv(i,j,k)  = mon_o(mon)%dt_dt_adv(i,j,k)   + ocn%dts_dt_adv(i,j,kr,1)  * sec_year * mon_avg ! °C/year
              mon_o(mon)%ds_dt_adv(i,j,k)  = mon_o(mon)%ds_dt_adv(i,j,k)   + ocn%dts_dt_adv(i,j,kr,2)  * sec_year * mon_avg ! psu/year
            endif
            mon_o(mon)%dt_dt_diff(i,j,k) = mon_o(mon)%dt_dt_diff(i,j,k)  + ocn%dts_dt_diff(i,j,kr,1) * sec_year * mon_avg ! °C/year
            mon_o(mon)%ds_dt_diff(i,j,k) = mon_o(mon)%ds_dt_diff(i,j,k)  + ocn%dts_dt_diff(i,j,kr,2) * sec_year * mon_avg ! psu/year
            mon_o(mon)%diffdia(i,j,k) = mon_o(mon)%diffdia(i,j,k)   + diff_dia(i,j,kr)   * mon_avg ! m2/s
            mon_o(mon)%drho_dx(i,j,k) = mon_o(mon)%drho_dx(i,j,k)   + drho_dx (i,j,kr)   * mon_avg ! m2/s
            mon_o(mon)%drho_dy(i,j,k) = mon_o(mon)%drho_dy(i,j,k)   + drho_dy (i,j,kr)   * mon_avg ! m2/s
            mon_o(mon)%drho_dz(i,j,k) = mon_o(mon)%drho_dz(i,j,k)   + drho_dz (i,j,kr)   * mon_avg ! m2/s
            mon_o(mon)%Ri(i,j,k)      = mon_o(mon)%Ri(i,j,k)        + Ri      (i,j,kr)   * mon_avg 
            mon_o(mon)%rho    (i,j,k) = mon_o(mon)%rho    (i,j,k)   + rho     (i,j,kr)   * mon_avg ! kg/m3
            mon_o(mon)%rho2   (i,j,k) = mon_o(mon)%rho2   (i,j,k)   + ocn%rho (i,j,kr)   * mon_avg ! kg/m3
          endif
          mon_o(mon)%u(i,j,k) = mon_o(mon)%u(i,j,k) + ocn%u(1,i,j,kr) * mon_avg ! m/s
          mon_o(mon)%v(i,j,k) = mon_o(mon)%v(i,j,k) + ocn%u(2,i,j,kr) * mon_avg ! m/s
          mon_o(mon)%w(i,j,k) = mon_o(mon)%w(i,j,k) + ocn%u(3,i,j,kr) * mon_avg ! m/s
          mon_o(mon)%fdx(i,j,k) = mon_o(mon)%fdx(i,j,k) + ocn%fdx(i,j,kr,1) * mon_avg ! m3*K
          mon_o(mon)%fdy(i,j,k) = mon_o(mon)%fdy(i,j,k) + ocn%fdy(i,j,kr,1) * mon_avg ! m3*K
          mon_o(mon)%fdz(i,j,k) = mon_o(mon)%fdz(i,j,k) + ocn%fdz(i,j,kr,1) * mon_avg ! m3*K
        enddo
        mon_o(mon)%t_atl  (j,k) = mon_o(mon)%t_atl  (j,k) + t_atl  (j,kr) * mon_avg
        mon_o(mon)%t_pac  (j,k) = mon_o(mon)%t_pac  (j,k) + t_pac  (j,kr) * mon_avg
        mon_o(mon)%t_ind  (j,k) = mon_o(mon)%t_ind  (j,k) + t_ind  (j,kr) * mon_avg
        mon_o(mon)%t_so   (j,k) = mon_o(mon)%t_so   (j,k) + t_so   (j,kr) * mon_avg
        mon_o(mon)%s_atl  (j,k) = mon_o(mon)%s_atl  (j,k) + s_atl  (j,kr) * mon_avg
        mon_o(mon)%s_pac  (j,k) = mon_o(mon)%s_pac  (j,k) + s_pac  (j,kr) * mon_avg
        mon_o(mon)%s_ind  (j,k) = mon_o(mon)%s_ind  (j,k) + s_ind  (j,kr) * mon_avg
        mon_o(mon)%s_so   (j,k) = mon_o(mon)%s_so   (j,k) + s_so   (j,kr) * mon_avg
        mon_o(mon)%rho_atl(j,k) = mon_o(mon)%rho_atl(j,k) + rho_atl(j,kr) * mon_avg
        mon_o(mon)%rho_pac(j,k) = mon_o(mon)%rho_pac(j,k) + rho_pac(j,kr) * mon_avg
        mon_o(mon)%rho_ind(j,k) = mon_o(mon)%rho_ind(j,k) + rho_ind(j,kr) * mon_avg
        mon_o(mon)%rho_so (j,k) = mon_o(mon)%rho_so (j,k) + rho_so (j,kr) * mon_avg
        mon_o(mon)%opsi   (j,k) = mon_o(mon)%opsi   (j,k) + opsi   (j,kr)*1.e-6_wp * mon_avg ! Sv
        mon_o(mon)%opsia  (j,k) = mon_o(mon)%opsia  (j,k) + opsia  (j,kr)*1.e-6_wp * mon_avg ! Sv
        mon_o(mon)%opsip  (j,k) = mon_o(mon)%opsip  (j,k) + opsip  (j,kr)*1.e-6_wp * mon_avg ! Sv
        mon_o(mon)%opsii  (j,k) = mon_o(mon)%opsii  (j,k) + opsii  (j,kr)*1.e-6_wp * mon_avg ! Sv
      enddo
    enddo

    mon_o(mon)%hft(1:3,:)   = mon_o(mon)%hft(1:3,:)   + hft*1.e-15_wp * mon_avg ! PW
    mon_o(mon)%hfp(1:3,:)   = mon_o(mon)%hfp(1:3,:)   + hfp*1.e-15_wp * mon_avg ! PW
    mon_o(mon)%hfa(1:3,:)   = mon_o(mon)%hfa(1:3,:)   + hfa*1.e-15_wp * mon_avg ! PW
    mon_o(mon)%hftz(:,:)    = mon_o(mon)%hftz(:,:)    + hftz*1.e-15_wp * mon_avg ! PW/m
    mon_o(mon)%hfpz(:,:)    = mon_o(mon)%hfpz(:,:)    + hfpz*1.e-15_wp * mon_avg ! PW/m
    mon_o(mon)%hfaz(:,:)    = mon_o(mon)%hfaz(:,:)    + hfaz*1.e-15_wp * mon_avg ! PW/m
    mon_o(mon)%fwt(1:3,:)   = mon_o(mon)%fwt(1:3,:)   + fwt*1.e-6_wp  * mon_avg ! Sv
    mon_o(mon)%fwp(1:3,:)   = mon_o(mon)%fwp(1:3,:)   + fwp*1.e-6_wp  * mon_avg ! Sv
    mon_o(mon)%fwa(1:3,:)   = mon_o(mon)%fwa(1:3,:)   + fwa*1.e-6_wp  * mon_avg ! Sv
    mon_o(mon)%fwtz(:,:)    = mon_o(mon)%fwtz(:,:)    + fwtz*1.e-6_wp  * mon_avg ! Sv/m
    mon_o(mon)%fwpz(:,:)    = mon_o(mon)%fwpz(:,:)    + fwpz*1.e-6_wp  * mon_avg ! Sv/m
    mon_o(mon)%fwaz(:,:)    = mon_o(mon)%fwaz(:,:)    + fwaz*1.e-6_wp  * mon_avg ! Sv/m


    ! compute zonal mean meridional velocity, global and for each basin
    do j=1,maxj-1
      do k=1,maxk
        vz(j,k)  = 0._wp
        vza(j,k) = 0._wp
        vzp(j,k) = 0._wp
        dxz(j,k)  = 0._wp
        dxza(j,k) = 0._wp
        dxzp(j,k) = 0._wp
        nx  = 0
        nxa = 0
        nxp = 0
        do i=1,maxi
          if (mask_v(i,j,k).eq.1) then
            nx = nx+1
            dxz(j,k) = dxz(j,k) + dxv(j) 
            vz(j,k) = vz(j,k) + ocn%u(2,i,j,k) 
            if (basin_mask2(i,j).eq.i_atlantic) then
              nxa = nxa+1
              dxza(j,k) = dxza(j,k) + dxv(j) 
              vza(j,k) = vza(j,k) + ocn%u(2,i,j,k) 
            endif
            if (basin_mask2(i,j).eq.i_pacific .or. basin_mask2(i,j).eq.i_indian) then
              nxp = nxp+1
              dxzp(j,k) = dxzp(j,k) + dxv(j) 
              vzp(j,k) = vzp(j,k) + ocn%u(2,i,j,k) 
            endif
          endif
        enddo
        if (nx>0)  vz(j,k)  = vz(j,k)/nx
        if (nxa>0) vza(j,k) = vza(j,k)/nxa
        if (nxp>0) vzp(j,k) = vzp(j,k)/nxp
      enddo
    enddo

    ! compute zonal mean potential temperature and salinity, global and for each basin
    do j=1,maxj
      do k=1,maxk
        tz(j,k)  = 0._wp
        tza(j,k) = 0._wp
        tzp(j,k) = 0._wp
        sz(j,k)  = 0._wp
        sza(j,k) = 0._wp
        szp(j,k) = 0._wp
        nx  = 0
        nxa = 0
        nxp = 0
        do i=1,maxi
          if (mask_c(i,j,k).eq.1) then
            nx = nx+1
            tz(j,k) = tz(j,k) + ocn%ts(i,j,k,1) 
            sz(j,k) = sz(j,k) + ocn%ts(i,j,k,2) 
            if (basin_mask2(i,j).eq.i_atlantic) then
              nxa = nxa+1
              tza(j,k) = tza(j,k) + ocn%ts(i,j,k,1) 
              sza(j,k) = sza(j,k) + ocn%ts(i,j,k,2) 
            endif
            if (basin_mask2(i,j).eq.i_pacific .or. basin_mask2(i,j).eq.i_indian) then
              nxp = nxp+1
              tzp(j,k) = tzp(j,k) + ocn%ts(i,j,k,1) 
              szp(j,k) = szp(j,k) + ocn%ts(i,j,k,2) 
            endif
          endif
        enddo
        if (nx>0)  tz(j,k)  = tz(j,k)/nx
        if (nxa>0) tza(j,k) = tza(j,k)/nxa
        if (nxp>0) tzp(j,k) = tzp(j,k)/nxp
        if (nx>0)  sz(j,k)  = sz(j,k)/nx
        if (nxa>0) sza(j,k) = sza(j,k)/nxa
        if (nxp>0) szp(j,k) = szp(j,k)/nxp
      enddo
    enddo

    ! meridional heat and freshwater transport by overturning
    hfto(:) = 0._wp
    hfao(:) = 0._wp
    hfpo(:) = 0._wp
    fwto(:) = 0._wp
    fwao(:) = 0._wp
    fwpo(:) = 0._wp
    do j=1,maxj-1
      do k=1,maxk
        hfto(j) = hfto(j) + vz(j,k)  * 0.5_wp*(tz(j,k) +tz(j+1,k))  * dxz(j,k)*dz(k)  * rho0*cap_w   ! m/s * K * m2 * kg/m3 * J/kg/K = W
        hfao(j) = hfao(j) + vza(j,k) * 0.5_wp*(tza(j,k)+tza(j+1,k)) * dxza(j,k)*dz(k) * rho0*cap_w   ! m/s * K * m2 * kg/m3 * J/kg/K = W
        hfpo(j) = hfpo(j) + vzp(j,k) * 0.5_wp*(tzp(j,k)+tzp(j+1,k)) * dxzp(j,k)*dz(k) * rho0*cap_w   ! m/s * K * m2 * kg/m3 * J/kg/K = W
        fwto(j) = fwto(j) + vz(j,k)  * 0.5_wp*(sz(j,k) +sz(j+1,k))  * dxz(j,k)*dz(k)  * (-1._wp/ocn%saln0)   ! m/s * psu * m2 / psu = m3/s 
        fwao(j) = fwao(j) + vza(j,k) * 0.5_wp*(sza(j,k)+sza(j+1,k)) * dxza(j,k)*dz(k) * (-1._wp/ocn%saln0)   ! m/s * psu * m2 / psu = m3/s 
        fwpo(j) = fwpo(j) + vzp(j,k) * 0.5_wp*(szp(j,k)+szp(j+1,k)) * dxzp(j,k)*dz(k) * (-1._wp/ocn%saln0)   ! m/s * psu * m2 / psu = m3/s 
      enddo
    enddo

    ! meridional heat and freshwater transport by gyres
    hftg(:) = 0._wp
    hfag(:) = 0._wp
    hfpg(:) = 0._wp
    fwtg(:) = 0._wp
    fwag(:) = 0._wp
    fwpg(:) = 0._wp
    do j=1,maxj-1
      do k=1,maxk
        do i=1,maxi
          if (mask_v(i,j,k).eq.1) then
            hftg(j) = hftg(j) + (ocn%u(2,i,j,k)-vz(j,k))  * 0.5_wp*(ocn%ts(i,j,k,1)+ocn%ts(i,j+1,k,1) - tz(j,k)-tz(j+1,k))  * dxv(j)*dz(k) * rho0*cap_w   ! m/s * K * m2 * kg/m3 * J/kg/K = W
            fwtg(j) = fwtg(j) + (ocn%u(2,i,j,k)-vz(j,k))  * 0.5_wp*(ocn%ts(i,j,k,2)+ocn%ts(i,j+1,k,2) - sz(j,k)-sz(j+1,k))  * dxv(j)*dz(k) * (-1._wp/ocn%saln0)   ! m/s * psu * m2 / psu = m3/s
            if (basin_mask2(i,j).eq.i_atlantic) then
              hfag(j) = hfag(j) + (ocn%u(2,i,j,k)-vza(j,k))  * 0.5_wp*(ocn%ts(i,j,k,1)+ocn%ts(i,j+1,k,1) - tza(j,k)-tza(j+1,k))  * dxv(j)*dz(k) * rho0*cap_w   ! m/s * K * m2 * kg/m3 * J/kg/K = W
              fwag(j) = fwag(j) + (ocn%u(2,i,j,k)-vza(j,k))  * 0.5_wp*(ocn%ts(i,j,k,2)+ocn%ts(i,j+1,k,2) - sza(j,k)-sza(j+1,k))  * dxv(j)*dz(k) * (-1._wp/ocn%saln0)   ! m/s * psu * m2 / psu = m3/s
            endif
            if (basin_mask2(i,j).eq.i_pacific .or. basin_mask2(i,j).eq.i_indian) then
              hfpg(j) = hfpg(j) + (ocn%u(2,i,j,k)-vzp(j,k))  * 0.5_wp*(ocn%ts(i,j,k,1)+ocn%ts(i,j+1,k,1) - tzp(j,k)-tzp(j+1,k))  * dxv(j)*dz(k) * rho0*cap_w   ! m/s * K * m2 * kg/m3 * J/kg/K = W
              fwpg(j) = fwpg(j) + (ocn%u(2,i,j,k)-vzp(j,k))  * 0.5_wp*(ocn%ts(i,j,k,2)+ocn%ts(i,j+1,k,2) - szp(j,k)-szp(j+1,k))  * dxv(j)*dz(k) * (-1._wp/ocn%saln0)   ! m/s * psu * m2 / psu = m3/s
            endif
          endif
        enddo
      enddo
    enddo

    mon_o(mon)%hft(4,:)   = mon_o(mon)%hft(4,:)   + hfto*1.e-15_wp * mon_avg ! PW
    mon_o(mon)%hft(5,:)   = mon_o(mon)%hft(5,:)   + hftg*1.e-15_wp * mon_avg ! PW
    mon_o(mon)%hfa(4,:)   = mon_o(mon)%hfa(4,:)   + hfao*1.e-15_wp * mon_avg ! PW
    mon_o(mon)%hfa(5,:)   = mon_o(mon)%hfa(5,:)   + hfag*1.e-15_wp * mon_avg ! PW
    mon_o(mon)%hfp(4,:)   = mon_o(mon)%hfp(4,:)   + hfpo*1.e-15_wp * mon_avg ! PW
    mon_o(mon)%hfp(5,:)   = mon_o(mon)%hfp(5,:)   + hfpg*1.e-15_wp * mon_avg ! PW
    mon_o(mon)%fwt(4,:)   = mon_o(mon)%fwt(4,:)   + fwto*1.e-6_wp * mon_avg ! Sv 
    mon_o(mon)%fwt(5,:)   = mon_o(mon)%fwt(5,:)   + fwtg*1.e-6_wp * mon_avg ! Sv 
    mon_o(mon)%fwa(4,:)   = mon_o(mon)%fwa(4,:)   + fwao*1.e-6_wp * mon_avg ! Sv 
    mon_o(mon)%fwa(5,:)   = mon_o(mon)%fwa(5,:)   + fwag*1.e-6_wp * mon_avg ! Sv 
    mon_o(mon)%fwp(4,:)   = mon_o(mon)%fwp(4,:)   + fwpo*1.e-6_wp * mon_avg ! Sv 
    mon_o(mon)%fwp(5,:)   = mon_o(mon)%fwp(5,:)   + fwpg*1.e-6_wp * mon_avg ! Sv 

    do j=1,maxj
      do i=1,maxi
        mon_o(mon)%fayti(i,j) = mon_o(mon)%fayti(i,j) + sum(ocn%fay(i,j,k1(i,j):maxk,1)*ocn%rho(i,j,k1(i,j):maxk))*cap_w/dt * mon_avg ! W
        mon_o(mon)%fdyti(i,j) = mon_o(mon)%fdyti(i,j) + sum(ocn%fdy(i,j,k1(i,j):maxk,1)*ocn%rho(i,j,k1(i,j):maxk))*cap_w/dt * mon_avg ! W
      enddo
    enddo

    mon_o(mon)%taux  = mon_o(mon)%taux  + ocn%tau(1,:,:) * mon_avg
    mon_o(mon)%tauy  = mon_o(mon)%tauy  + ocn%tau(2,:,:) * mon_avg
    mon_o(mon)%ssh = mon_o(mon)%ssh  + ocn%ssh * mon_avg
    mon_o(mon)%ub    = mon_o(mon)%ub    + ocn%ub(1,1:maxi,1:maxj) * mon_avg ! m/s
    mon_o(mon)%vb    = mon_o(mon)%vb    + ocn%ub(2,1:maxi,1:maxj) * mon_avg ! m/s
    where (ocn%f_ocn.gt.0._wp) 
      mon_o(mon)%sst   = mon_o(mon)%sst   + ocn%ts(:,:,maxk,1) * mon_avg ! °C
      mon_o(mon)%sss   = mon_o(mon)%sss   + ocn%ts(:,:,maxk,2) * mon_avg ! psu
      mon_o(mon)%psi   = mon_o(mon)%psi   + ocn%psi(1:maxi,1:maxj)/rho0*1.e-6_wp * mon_avg ! Sv
      mon_o(mon)%dt_dt_flxsur = mon_o(mon)%dt_dt_flxsur - ocn%flx_sur(:,:,1)/dz(maxk) * sec_year * mon_avg ! °C/year
      mon_o(mon)%ds_dt_flxsur = mon_o(mon)%ds_dt_flxsur - ocn%flx_sur(:,:,2)/dz(maxk) * sec_year * mon_avg ! psu/year
      mon_o(mon)%flx   = mon_o(mon)%flx   + ocn%flx                              * mon_avg ! W/m2
      mon_o(mon)%fw    = mon_o(mon)%fw    + (ocn%p_e_sic+ocn%runoff+ocn%calving+ocn%bmelt+ocn%fw_hosing+ocn%fw_flux_adj)*sec_day         * mon_avg ! kg/m2/day
      mon_o(mon)%vsf   = mon_o(mon)%vsf   + ocn%flx_sur(:,:,2)*rho0/ocn%saln0*sec_day                   * mon_avg ! kg/m2/day
      mon_o(mon)%p_e_sic   = mon_o(mon)%p_e_sic   + ocn%p_e_sic*sec_day         * mon_avg ! kg/m2/day
      mon_o(mon)%runoff= mon_o(mon)%runoff+ ocn%runoff*sec_day         * mon_avg ! kg/m2/day
      mon_o(mon)%runoffSv= mon_o(mon)%runoffSv+ ocn%runoff*ocn%grid%ocn_area/rho0*1.e-6_wp         * mon_avg ! Sv
      mon_o(mon)%runoffSv_ice= mon_o(mon)%runoffSv_ice+ ocn%runoff_ice*ocn%grid%ocn_area/rho0*1.e-6_wp         * mon_avg ! Sv
      mon_o(mon)%calving= mon_o(mon)%calving+ ocn%calving*sec_day         * mon_avg ! kg/m2/day
      mon_o(mon)%bmelt = mon_o(mon)%bmelt + ocn%bmelt*sec_day         * mon_avg ! kg/m2/day
      mon_o(mon)%fw_hosing = mon_o(mon)%fw_hosing + ocn%fw_hosing*sec_day        * mon_avg ! kg/m2/day
      mon_o(mon)%fw_flux_adj = mon_o(mon)%fw_flux_adj + ocn%fw_flux_adj*sec_day  * mon_avg ! kg/m2/day
      mon_o(mon)%fw_noise = mon_o(mon)%fw_noise + ocn%fw_noise*sec_day       * mon_avg ! kg/m2/day
      mon_o(mon)%flx_noise = mon_o(mon)%flx_noise + ocn%flx_noise        * mon_avg ! W/m2
      mon_o(mon)%nconv = mon_o(mon)%nconv + real(ocn%nconv,wp)                   * mon_avg 
      mon_o(mon)%kven  = mon_o(mon)%kven  + real(ocn%kven,wp)                    * mon_avg
      mon_o(mon)%dconv = mon_o(mon)%dconv + ocn%dconv                            * mon_avg
      mon_o(mon)%dven  = mon_o(mon)%dven  + ocn%dven                             * mon_avg
      mon_o(mon)%mld   = mon_o(mon)%mld   + (-ocn%mld)                           * mon_avg
      mon_o(mon)%ke_tau = mon_o(mon)%ke_tau   + ocn%ke_tau/dt*1000._wp     * mon_avg ! mW/m2
    elsewhere
      mon_o(mon)%sst         = missing_value 
      mon_o(mon)%sss         = missing_value 
      mon_o(mon)%psi         = missing_value 
      mon_o(mon)%dt_dt_flxsur = missing_value 
      mon_o(mon)%ds_dt_flxsur = missing_value 
      mon_o(mon)%flx         = missing_value 
      mon_o(mon)%fw          = missing_value 
      mon_o(mon)%vsf         = missing_value 
      mon_o(mon)%p_e_sic         = missing_value 
      mon_o(mon)%runoff      = missing_value 
      mon_o(mon)%runoffSv    = missing_value 
      mon_o(mon)%runoffSv_ice= missing_value 
      mon_o(mon)%calving     = missing_value 
      mon_o(mon)%bmelt       = missing_value 
      mon_o(mon)%fw_hosing   = missing_value 
      mon_o(mon)%fw_flux_adj = missing_value 
      mon_o(mon)%fw_noise    = missing_value 
      mon_o(mon)%flx_noise   = missing_value 
      mon_o(mon)%nconv       = missing_value 
      mon_o(mon)%kven        = missing_value 
      mon_o(mon)%dconv       = missing_value 
      mon_o(mon)%dven        = missing_value 
      mon_o(mon)%mld         = missing_value 
      mon_o(mon)%ke_tau      = missing_value 
    endwhere
    do n=1,n_isles
      mon_o(mon)%ubisl(:,:,n) = mon_o(mon)%ubisl(:,:,n) + ocn%ub_isl(1,1:maxi,1:maxj,n)         * mon_avg ! m/s
      mon_o(mon)%vbisl(:,:,n) = mon_o(mon)%vbisl(:,:,n) + ocn%ub_isl(2,1:maxi,1:maxj,n)         * mon_avg ! m/s
    enddo

    ! bottom ocean temperature and salinity
    do j=1,maxj
      do i=1,maxi
        if (ocn%f_ocn(i,j).gt.0._wp) then
          mon_o(mon)%tbot(i,j) = mon_o(mon)%tbot(i,j) + ocn%ts(i,j,k1(i,j),1) * mon_avg ! °C
          mon_o(mon)%sbot(i,j) = mon_o(mon)%sbot(i,j) + ocn%ts(i,j,k1(i,j),2) * mon_avg ! °C
        else
          mon_o(mon)%tbot(i,j) = missing_value
          mon_o(mon)%sbot(i,j) = missing_value 
        endif
      enddo
    enddo

    ! buoyancy
    do j=1,maxj
      do i=1,maxi
        if (ocn%f_ocn(i,j).gt.0._wp) then
          if (i_alphabeta.eq.1) then
            rho1 = eos(ocn%ts(i,j,maxk,1)+0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
            rho2 = eos(ocn%ts(i,j,maxk,1)-0.5_wp,ocn%ts(i,j,maxk,2),0._wp)    
            alpha = rho2-rho1
          else if (i_alphabeta.eq.2) then
            alpha = (519._wp+122._wp*ocn%ts(i,j,maxk,1))*1.e-4_wp    ! kg/m3/K
          endif
          mon_o(mon)%buoyT(i,j) = mon_o(mon)%buoyT(i,j) + g/cap_w*ocn%flx(i,j)*alpha/rho0 * mon_avg ! m/s2 * kg*K/J * J/s/m2 * kg/m3/K * m3/kg = kg/m/s3 = N/m2/s
          if (i_alphabeta.eq.1) then
            rho1 = eos(ocn%ts(i,j,maxk,1),ocn%ts(i,j,maxk,2)+0.5_wp,0._wp)    
            rho2 = eos(ocn%ts(i,j,maxk,1),ocn%ts(i,j,maxk,2)-0.5_wp,0._wp)    
            beta = rho1-rho2
          else if (i_alphabeta.eq.2) then
            beta = 0.8_wp   ! kg/m3/psu
          endif
          if (i_fwf_buoy.eq.1) then 
            mon_o(mon)%buoyS(i,j) = mon_o(mon)%buoyS(i,j) + g*ocn%fw_corr(i,j)*beta*ocn%saln0/rho0 * mon_avg  ! m/s2 * kg/m2/s kg/m3/psu * m3/kg * psu = kg*m/s2 = N/m2/s
          else if (i_fwf_buoy.eq.2) then 
            mon_o(mon)%buoyS(i,j) = mon_o(mon)%buoyS(i,j) + g*ocn%flx_sur(i,j,2)*beta * mon_avg  ! m/s2 * kg/m2/s kg/m3/psu * m3/kg * psu = kg*m/s2 = N/m2/s
          endif
          mon_o(mon)%buoy(i,j) = mon_o(mon)%buoyT(i,j) + mon_o(mon)%buoyS(i,j)
        endif
      enddo
    enddo

    ! mixed layer depth from sigma t criterion
    do j=1,maxj
      do i=1,maxi
        if (ocn%f_ocn(i,j).gt.0._wp) then
          mldst = -zw(k1(i,j))
          do k=maxk,k1(i,j),-1
            if (rho(i,j,k)-(rho(i,j,maxk)).gt.0.125_wp) then
              mldst = -zw(k)
              exit
            endif
          enddo
          mon_o(mon)%mldst(i,j) = mon_o(mon)%mldst(i,j) + mldst * mon_avg
        endif
      enddo
    enddo

     if (l_daily_output) then
       day_o(doy)%t = ocn%ts(:,:,:,1)  ! °C
       day_o(doy)%s = ocn%ts(:,:,:,2)  ! psu
       day_o(doy)%rho = rho(:,:,:)  ! kg/m3
       day_o(doy)%mld = -ocn%mld(:,:)  ! m
       if (cons_tracer) day_o(doy)%cons = ocn%ts(:,:,:,i_cons)
       day_o(doy)%u = ocn%u(1,1:maxi,1:maxj,1:maxk)  ! m/s
       day_o(doy)%v = ocn%u(2,1:maxi,1:maxj,1:maxk)  ! m/s
       day_o(doy)%w = ocn%u(3,1:maxi,1:maxj,1:maxk)  ! m/s
       day_o(doy)%ub= ocn%ub(1,1:maxi,1:maxj) * mon_avg ! m/s
       day_o(doy)%vb= ocn%ub(2,1:maxi,1:maxj) * mon_avg ! m/s
     endif

       deallocate(rho)
   endif

   if (time_out_ocn .and. time_eoy_ocn) then

     ann_o%f_ocn = ocn%f_ocn
     ann_o%mask_ocn = ocn%grid%mask_ocn
     ann_o%area = ocn%f_ocn*ocn%grid%ocn_area
     ann_o%topo = topo
     ann_o%bathy = bathy
     ann_o%k1 = ocn%grid%k1
     ann_o%drag  = drag(1,1:maxi,:)
     ann_o%drag_bcl  = drag_bcl(1,1:maxi,:)
     ann_o%map_isles = map_isles
     do n=1,n_isles+1
       ann_o%map_edge(:,:,n) = map_edge(:,:,n)
     enddo
     do n=n_isles+2,maxisles
       ann_o%map_edge(:,:,n) = missing_value
     enddo

     ann_o%q_geo = ocn%q_geo

     do j=1,maxj
       do i=1,maxi
         ann_o%mldmax(i,j) = 0._wp
         do m=1,nmon_year
           ann_o%mldmax(i,j) = max(ann_o%mldmax(i,j),mon_o(m)%mld(i,j))
         enddo
       enddo
     enddo

     do j=1,maxj
       do k=1,maxk
         kr = maxk-k+1
         do i=1,maxi
           if (mask_c(i,j,kr).eq.1) then
             ann_o%vol(i,j,k) = ocn_vol(i,j,kr)
             if (age_tracer)  ann_o%age(i,j,k)  = ocn%ts(i,j,kr,i_age)
             if (dye_tracer)  ann_o%dye(i,j,k)  = ocn%ts(i,j,kr,i_dye)
             if (cons_tracer) ann_o%cons(i,j,k) = ocn%ts(i,j,kr,i_cons)
             if (l_cfc) then
               ann_o%cfc11(i,j,k)  = ocn%ts(i,j,kr,i_cfc11) * 1.e9_wp ! mol/m3 -> pmol/l
               ann_o%cfc12(i,j,k)  = ocn%ts(i,j,kr,i_cfc12) * 1.e9_wp ! mol/m3 -> pmol/l
             endif
           else
             ann_o%vol(i,j,k)  = missing_value 
             if (age_tracer)  ann_o%age(i,j,k)  = missing_value 
             if (dye_tracer)  ann_o%dye(i,j,k)  = missing_value 
             if (cons_tracer) ann_o%cons(i,j,k) = missing_value 
             if (l_cfc) then
               ann_o%cfc11(i,j,k)  = missing_value 
               ann_o%cfc12(i,j,k)  = missing_value 
             endif
           endif
         enddo
         ntot = count(basin_mask2(:,j).eq.i_atlantic .and. kr.ge.k1(1:maxi,j))
         if (ntot.gt.0) then
           if (age_tracer) ann_o%a_atl(j,k) = sum(ocn%ts(:,j,kr,i_age), basin_mask2(:,j).eq.i_atlantic .and. kr.ge.k1(1:maxi,j)) / ntot
           if (dye_tracer) ann_o%d_atl(j,k) = sum(ocn%ts(:,j,kr,i_dye), basin_mask2(:,j).eq.i_atlantic .and. kr.ge.k1(1:maxi,j)) / ntot
         else
           if (age_tracer) ann_o%a_atl(j,k) = missing_value
           if (dye_tracer) ann_o%d_atl(j,k) = missing_value
         endif
         ntot = count(basin_mask2(:,j).eq.i_pacific .and. kr.ge.k1(1:maxi,j))
         if (ntot.gt.0) then
           if (age_tracer) ann_o%a_pac(j,k) = sum(ocn%ts(:,j,kr,i_age), basin_mask2(:,j).eq.i_pacific .and. kr.ge.k1(1:maxi,j)) / ntot
           if (dye_tracer) ann_o%d_pac(j,k) = sum(ocn%ts(:,j,kr,i_dye), basin_mask2(:,j).eq.i_pacific .and. kr.ge.k1(1:maxi,j)) / ntot
         else
           if (age_tracer) ann_o%a_pac(j,k) = missing_value
           if (dye_tracer) ann_o%d_pac(j,k) = missing_value
         endif
         ntot = count(basin_mask2(:,j).eq.i_indian .and. kr.ge.k1(1:maxi,j))
         if (ntot.gt.0) then
           if (age_tracer) ann_o%a_ind(j,k) = sum(ocn%ts(:,j,kr,i_age), basin_mask2(:,j).eq.i_indian .and. kr.ge.k1(1:maxi,j)) / ntot
           if (dye_tracer) ann_o%d_ind(j,k) = sum(ocn%ts(:,j,kr,i_dye), basin_mask2(:,j).eq.i_indian .and. kr.ge.k1(1:maxi,j)) / ntot
         else
           if (age_tracer) ann_o%a_ind(j,k) = missing_value
           if (dye_tracer) ann_o%d_ind(j,k) = missing_value
         endif
         ntot = count(basin_mask(:,j).eq.i_southern .and. kr.ge.k1(1:maxi,j))
         if (ntot.gt.0) then
           if (age_tracer) ann_o%a_so(j,k)  = sum(ocn%ts(:,j,kr,i_age), basin_mask(:,j).eq.i_southern .and. kr.ge.k1(1:maxi,j)) / ntot
           if (dye_tracer) ann_o%d_so(j,k)  = sum(ocn%ts(:,j,kr,i_dye), basin_mask(:,j).eq.i_southern .and. kr.ge.k1(1:maxi,j)) / ntot
         else
           if (age_tracer) ann_o%a_so(j,k)  = missing_value
           if (dye_tracer) ann_o%d_so(j,k)  = missing_value
         endif
       enddo
     enddo

     call ocn_diag_out

    endif


   return

  end subroutine ocn_diag


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  o c n _ d i a g _ o u t
  ! Purpose  :  write ocean netcdf output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocn_diag_out

    implicit none

    integer :: k, ncid
    character (len=256) :: fnm


    nout = nout + 1

    ! Get annual values
    call ocn_ave( mon_o,ann_o )

    ! write to file
    fnm = trim(out_dir)//"/ocn.nc"
    call nc_open(fnm,ncid)
    call nc_write(fnm,dim_time,real(year_now,wp), dim1=dim_time, start=[nout], count=[1],ncid=ncid)
    do k = 1, nmon_year
       call ocn_nc_write(fnm,ncid,mon_o(k),k,nout)
    end do
    call ocn_nc_write(fnm,ncid,ann_o,nmon_year+1,nout)
    call nc_close(ncid)

    if (l_daily_output) then
      ! write to file
      fnm = trim(out_dir)//"/ocn_daily.nc"
      call nc_open(fnm,ncid)
      call nc_write(fnm,dim_time,real(year_now,wp), dim1=dim_time, start=[nout], count=[1],ncid=ncid)    
      do k = 1, nday_year
        call ocn_daily_nc_write(fnm,ncid,day_o(k),k,nout)
      end do
      call nc_close(ncid)
    endif


   return

  end subroutine ocn_diag_out
  

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t s _ n c
  ! Purpose  :  initialize netcdf file for time series output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ts_nc(fnm)

    implicit none

    character (len=*) :: fnm
    real(8) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_write_dim(fnm,dim_time, x=empty_time, units="years BP", unlimited=.TRUE.)
    call nc_write_dim(fnm,dim_lat, x=latv_buoy, axis="y", units="degN")

    return

  end subroutine ts_nc
 

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t s _ n c _ w r i t e
  ! Purpose  :  write time series to netcdf
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ts_nc_write(fnm,vars,ndat,y)

    implicit none

    type(ts_out) :: vars(:)

    character (len=*) :: fnm
    integer :: ndat, y, ncid, n, i


    call nc_open(fnm,ncid)
    call nc_write(fnm,"time", real([(i,i=(year_now-(y-1)*n_accel),(year_now),(n_accel))],wp), &
  dim1=dim_time,start=[ndat],count=[y],ncid=ncid)    
    call nc_write(fnm,"sst",     vars%sst,    dim1=dim_time,start=[ndat],count=[y],long_name="sea surface temperature",units="C",ncid=ncid)
    call nc_write(fnm,"sss",     vars%sss,    dim1=dim_time,start=[ndat],count=[y],long_name="sea surface salinity",units="psu",ncid=ncid)
    call nc_write(fnm,"tvol",    vars%t,      dim1=dim_time,start=[ndat],count=[y],long_name="volume averaged ocean potential temperature",units="C",ncid=ncid)
    call nc_write(fnm,"svol",    vars%s,      dim1=dim_time,start=[ndat],count=[y],long_name="volume averaged ocean salinity",units="psu",ncid=ncid)
    call nc_write(fnm,"t_deep",  vars%tdocn,  dim1=dim_time,start=[ndat],count=[y],long_name="volume averaged deep ocean potential temperature (below 2500 m)",units="C",ncid=ncid)
    call nc_write(fnm,"t_deep_atl",  vars%tdocn_atl,  dim1=dim_time,start=[ndat],count=[y],long_name="volume averaged deep Atlantic ocean potential temperature (below 2500 m)",units="C",ncid=ncid)
    call nc_write(fnm,"t_deep_pac",  vars%tdocn_pac,  dim1=dim_time,start=[ndat],count=[y],long_name="volume averaged deep Pacific ocean potential temperature (below 2500 m)",units="C",ncid=ncid)
    call nc_write(fnm,"t_deep_ind",  vars%tdocn_ind,  dim1=dim_time,start=[ndat],count=[y],long_name="volume averaged deep Indian ocean potential temperature (below 2500 m)",units="C",ncid=ncid)
    call nc_write(fnm,"t_deep_so ",  vars%tdocn_so ,  dim1=dim_time,start=[ndat],count=[y],long_name="volume averaged deep Southern ocean potential temperature (below 2500 m)",units="C",ncid=ncid)
    if (cons_tracer) call nc_write(fnm,"cons",    vars%cons,   dim1=dim_time,start=[ndat],count=[y],long_name="conservative tracer",units="/",ncid=ncid)
    if (age_tracer) call nc_write(fnm,"age",    vars%age,   dim1=dim_time,start=[ndat],count=[y],long_name="volume averaged ideal age tracer",units="years",ncid=ncid)
    if (dye_tracer) call nc_write(fnm,"dye",    vars%dye,   dim1=dim_time,start=[ndat],count=[y],long_name="volume averaged dye tracer",units="/",ncid=ncid)
    call nc_write(fnm,"ohc",vars%ohc(1),dim1=dim_time,start=[ndat],count=[y],long_name="global ocean heat content ",units="J",ncid=ncid)
    call nc_write(fnm,"ohc_atl", vars%ohc(2),dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic ocean heat content ",units="J",ncid=ncid)
    call nc_write(fnm,"ohc_natl",vars%ohc(6),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic ocean heat content ",units="J",ncid=ncid)
    call nc_write(fnm,"ohc_pac", vars%ohc(3),dim1=dim_time,start=[ndat],count=[y],long_name="Pacific ocean heat content ",units="J",ncid=ncid)
    call nc_write(fnm,"ohc_ind", vars%ohc(4),dim1=dim_time,start=[ndat],count=[y],long_name="Indian ocean ocean heat content ",units="J",ncid=ncid)
    call nc_write(fnm,"ohc_so",  vars%ohc(5),dim1=dim_time,start=[ndat],count=[y],long_name="Southern ocean ocean heat content ",units="J",ncid=ncid)
    call nc_write(fnm,"ohc700",vars%ohc700(1),dim1=dim_time,start=[ndat],count=[y],long_name="global ocean heat content (top 700 m)",units="J",ncid=ncid)
    call nc_write(fnm,"ohc700_atl", vars%ohc700(2),dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic ocean heat content (top 700 m)",units="J",ncid=ncid)
    call nc_write(fnm,"ohc700_natl",vars%ohc700(6),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic ocean heat content (top 700 m)",units="J",ncid=ncid)
    call nc_write(fnm,"ohc700_pac", vars%ohc700(3),dim1=dim_time,start=[ndat],count=[y],long_name="Pacific ocean heat content (top 700 m)",units="J",ncid=ncid)
    call nc_write(fnm,"ohc700_ind", vars%ohc700(4),dim1=dim_time,start=[ndat],count=[y],long_name="Indian ocean ocean heat content (top 700 m)",units="J",ncid=ncid)
    call nc_write(fnm,"ohc700_so",  vars%ohc700(5),dim1=dim_time,start=[ndat],count=[y],long_name="Southern ocean ocean heat content (top 700 m)",units="J",ncid=ncid)
    call nc_write(fnm,"ohc2000",vars%ohc2000(1),dim1=dim_time,start=[ndat],count=[y],long_name="global ocean heat content (top 2000 m)",units="J",ncid=ncid)
    call nc_write(fnm,"ohc2000_atl", vars%ohc2000(2),dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic ocean heat content (top 2000 m)",units="J",ncid=ncid)
    call nc_write(fnm,"ohc2000_natl",vars%ohc2000(6),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic ocean heat content (top 2000 m)",units="J",ncid=ncid)
    call nc_write(fnm,"ohc2000_pac", vars%ohc2000(3),dim1=dim_time,start=[ndat],count=[y],long_name="Pacific ocean heat content (top 2000 m)",units="J",ncid=ncid)
    call nc_write(fnm,"ohc2000_ind", vars%ohc2000(4),dim1=dim_time,start=[ndat],count=[y],long_name="Indian ocean ocean heat content (top 2000 m)",units="J",ncid=ncid)
    call nc_write(fnm,"ohc2000_so",  vars%ohc2000(5),dim1=dim_time,start=[ndat],count=[y],long_name="Southern ocean ocean heat content (top 2000 m)",units="J",ncid=ncid)
    call nc_write(fnm,"amoc26N", vars%amoc26N,dim1=dim_time,start=[ndat],count=[y],long_name="maximum Atlantic overturning at 26N (RAPID)",units="Sv",ncid=ncid)
    call nc_write(fnm,"omaxa",   vars%omaxa,  dim1=dim_time,start=[ndat],count=[y],long_name="maximum Atlantic overturning",units="Sv",ncid=ncid)
    call nc_write(fnm,"omina",   vars%omina,  dim1=dim_time,start=[ndat],count=[y],long_name="minimum Atlantic overturning",units="Sv",ncid=ncid)
    call nc_write(fnm,"oinfa",   vars%oinfa,  dim1=dim_time,start=[ndat],count=[y],long_name="maximum southern Atlantic inflow",units="Sv",ncid=ncid)
    call nc_write(fnm,"omaxp",   vars%omaxp,  dim1=dim_time,start=[ndat],count=[y],long_name="maximum Pacific overturning",units="Sv",ncid=ncid)
    call nc_write(fnm,"pmoc",    vars%omaxp,  dim1=dim_time,start=[ndat],count=[y],long_name="maximum Pacific overturning",units="Sv",ncid=ncid)
    call nc_write(fnm,"ominp",   vars%ominp,  dim1=dim_time,start=[ndat],count=[y],long_name="minimum Pacific overturning",units="Sv",ncid=ncid)
    call nc_write(fnm,"omaxs",   vars%omaxs,  dim1=dim_time,start=[ndat],count=[y],long_name="maximum Southern ocean overturning",units="Sv",ncid=ncid)
    call nc_write(fnm,"omins",   vars%omins,  dim1=dim_time,start=[ndat],count=[y],long_name="minimum Southern Ocean overturning",units="Sv",ncid=ncid)
    call nc_write(fnm,"aabw",    -vars%omins,  dim1=dim_time,start=[ndat],count=[y],long_name="Antarctic bottom water formation rate",units="Sv",ncid=ncid)
    call nc_write(fnm,"hmaxa",   vars%hmaxa,  dim1=dim_time,start=[ndat],count=[y],long_name="maximum Atlantic northward meridional heat transport",units="PW",ncid=ncid)
    call nc_write(fnm,"hN55a",   vars%h55a,  dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic northward meridional heat transport at 55N",units="PW",ncid=ncid)
    call nc_write(fnm,"hmaxp",   vars%hmaxp,  dim1=dim_time,start=[ndat],count=[y],long_name="maximum Pacific northward meridional heat transport",units="PW",ncid=ncid)
    call nc_write(fnm,"hmaxs",   vars%hmaxs,  dim1=dim_time,start=[ndat],count=[y],long_name="maximum southward meridional ocean heat transport",units="PW",ncid=ncid)
    call nc_write(fnm,"h60s",    vars%h60s,   dim1=dim_time,start=[ndat],count=[y],long_name="southward meridional ocean heat transport at 60S",units="PW",ncid=ncid)
    call nc_write(fnm,"fmaxa",   vars%fmaxa,  dim1=dim_time,start=[ndat],count=[y],long_name="maximum Atlantic meridional freshwater transport",units="Sv",ncid=ncid)
    call nc_write(fnm,"fw_glob", vars%fw(1),dim1=dim_time,start=[ndat],count=[y],long_name="global net freshwater flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"fw_atl",  vars%fw(2),dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic net freshwater flux to ocean (-0.28 +/- 0.04 Talley2008)",units="Sv",ncid=ncid)
    call nc_write(fnm,"fw_atlN", vars%fw(8),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic net freshwater flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"fw_atlN30", vars%fw(6),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>30N) net freshwater flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"fw_atlN50", vars%fw(7),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>50N) net freshwater flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"fw_atlN55", vars%fw(9),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>55N) net freshwater flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"fw_pac",  vars%fw(3),dim1=dim_time,start=[ndat],count=[y],long_name="Pacific net freshwater flux to ocean (0.04 +/- 0.09 Talley2008)",units="Sv",ncid=ncid)
    call nc_write(fnm,"fw_ind",  vars%fw(4),dim1=dim_time,start=[ndat],count=[y],long_name="Indian ocean net freshwater flux to ocean (-0.37 +/- 0-10 Talley2008)",units="Sv",ncid=ncid)
    call nc_write(fnm,"fw_so",   vars%fw(5),dim1=dim_time,start=[ndat],count=[y],long_name="Southern ocean net freshwater flux to ocean (0.61 +/- 0.13 Talley2008)",units="Sv",ncid=ncid)
    call nc_write(fnm,"fw_soS60", vars%fw(10),dim1=dim_time,start=[ndat],count=[y],long_name="Southern ocean (<60S) net freshwater flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"fw_corr_glob", vars%fw_corr(1),dim1=dim_time,start=[ndat],count=[y],long_name="global corrected net freshwater flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"saln0", vars%saln0,dim1=dim_time,start=[ndat],count=[y],long_name="reference salinity",units="psu",ncid=ncid)
    call nc_write(fnm,"dvsf", vars%dvsf,dim1=dim_time,start=[ndat],count=[y],long_name="global virtual salinity flux correction",units="Sv",ncid=ncid)
    call nc_write(fnm,"hosing",  vars%hosing,dim1=dim_time,start=[ndat],count=[y],long_name="Freshwater hosing flux",units="Sv",ncid=ncid)
    call nc_write(fnm,"fw_noise",  vars%fw_noise,dim1=dim_time,start=[ndat],count=[y],long_name="Global freshwater flux from noise",units="Sv",ncid=ncid)
    call nc_write(fnm,"noise_fw",  vars%noise_fw,dim1=dim_time,start=[ndat],count=[y],long_name="Noise applied to freshwater flux",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"noise_flx",  vars%noise_flx,dim1=dim_time,start=[ndat],count=[y],long_name="Noise applied to heat flux",units="W/m2",ncid=ncid)
    call nc_write(fnm,"vsf_glob", vars%vsf(1),dim1=dim_time,start=[ndat],count=[y],long_name="global virtual salinity flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"vsf_atl",  vars%vsf(2),dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic virtual salinity flux",units="Sv",ncid=ncid)
    call nc_write(fnm,"vsf_atlN", vars%vsf(8),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic virtual salinity flux",units="Sv",ncid=ncid)
    call nc_write(fnm,"vsf_atlN30", vars%vsf(6),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>30N) virtual salinity flux",units="Sv",ncid=ncid)
    call nc_write(fnm,"vsf_atlN50", vars%vsf(7),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>50N) virtual salinity flux",units="Sv",ncid=ncid)
    call nc_write(fnm,"vsf_atlN55", vars%vsf(9),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>55N) virtual salinity flux",units="Sv",ncid=ncid)
    call nc_write(fnm,"vsf_pac",  vars%vsf(3),dim1=dim_time,start=[ndat],count=[y],long_name="Pacific virtual salinity flux",units="Sv",ncid=ncid)
    call nc_write(fnm,"vsf_ind",  vars%vsf(4),dim1=dim_time,start=[ndat],count=[y],long_name="Indian ocean virtual salinity flux",units="Sv",ncid=ncid)
    call nc_write(fnm,"vsf_so",   vars%vsf(5),dim1=dim_time,start=[ndat],count=[y],long_name="Southern ocean virtual salinity flux",units="Sv",ncid=ncid)
    call nc_write(fnm,"vsf_soS60", vars%vsf(10),dim1=dim_time,start=[ndat],count=[y],long_name="Southern ocean (<60S) virtual salinity flux",units="Sv",ncid=ncid)
    call nc_write(fnm,"p_e_sic_glob", vars%p_e_sic(1),dim1=dim_time,start=[ndat],count=[y],long_name="global P-E+sea ice freshwater flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"p_e_sic_atl",  vars%p_e_sic(2),dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic P-E+sea ice freshwater flux",units="Sv",ncid=ncid)
    call nc_write(fnm,"p_e_sic_atlN", vars%p_e_sic(8),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic P-E+sea ice freshwater flux",units="Sv",ncid=ncid)
    call nc_write(fnm,"p_e_sic_atlN30", vars%p_e_sic(6),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>30N) P-E+sea ice freshwater flux",units="Sv",ncid=ncid)
    call nc_write(fnm,"p_e_sic_atlN50", vars%p_e_sic(7),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>50N) P-E+sea ice freshwater flux",units="Sv",ncid=ncid)
    call nc_write(fnm,"p_e_sic_atlN55", vars%p_e_sic(9),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>55N) P-E+sea ice freshwater flux",units="Sv",ncid=ncid)
    call nc_write(fnm,"p_e_sic_pac",  vars%p_e_sic(3),dim1=dim_time,start=[ndat],count=[y],long_name="Pacific P-E+sea ice freshwater flux",units="Sv",ncid=ncid)
    call nc_write(fnm,"p_e_sic_ind",  vars%p_e_sic(4),dim1=dim_time,start=[ndat],count=[y],long_name="Indian ocean P-E+sea ice freshwater flux",units="Sv",ncid=ncid)
    call nc_write(fnm,"p_e_sic_so",   vars%p_e_sic(5),dim1=dim_time,start=[ndat],count=[y],long_name="Southern ocean P-E+sea ice freshwater flux",units="Sv",ncid=ncid)
    call nc_write(fnm,"p_e_sic_soS60",   vars%p_e_sic(10),dim1=dim_time,start=[ndat],count=[y],long_name="Southern ocean (<60S) P-E+sea ice freshwater flux",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_glob", vars%runoff(1),dim1=dim_time,start=[ndat],count=[y],long_name="global runoff flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_atl",  vars%runoff(2),dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic runoff flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_atlN", vars%runoff(8),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic runoff flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_atlN30", vars%runoff(6),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>30N) runoff flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_atlN50", vars%runoff(7),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>50N) runoff flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_atlN55", vars%runoff(9),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>55N) runoff flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_pac",  vars%runoff(3),dim1=dim_time,start=[ndat],count=[y],long_name="Pacific runoff flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_ind",  vars%runoff(4),dim1=dim_time,start=[ndat],count=[y],long_name="Indian ocean runoff flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_so",   vars%runoff(5),dim1=dim_time,start=[ndat],count=[y],long_name="Southern ocean runoff flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_veg_glob", vars%runoff_veg(1),dim1=dim_time,start=[ndat],count=[y],long_name="global runoff flux from ice- and lake-free land to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_veg_atl",  vars%runoff_veg(2),dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic runoff flux from ice- and lake-free land to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_veg_atlN30", vars%runoff_veg(6),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>30N) runoff flux from ice- and lake-free land to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_veg_atlN50", vars%runoff_veg(7),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>50N) runoff flux from ice- and lake-free land to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_veg_pac",  vars%runoff_veg(3),dim1=dim_time,start=[ndat],count=[y],long_name="Pacific runoff flux from ice- and lake-free land to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_veg_ind",  vars%runoff_veg(4),dim1=dim_time,start=[ndat],count=[y],long_name="Indian ocean runoff flux from ice- and lake-free land to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_veg_so",   vars%runoff_veg(5),dim1=dim_time,start=[ndat],count=[y],long_name="Southern ocean runoff flux from ice- and lake-free land to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_ice_glob", vars%runoff_ice(1),dim1=dim_time,start=[ndat],count=[y],long_name="global runoff flux from ice sheets to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_ice_atl",  vars%runoff_ice(2),dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic runoff flux from ice sheets to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_ice_atlN30", vars%runoff_ice(6),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>30N) runoff flux from ice sheets to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_ice_atlN50", vars%runoff_ice(7),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>50N) runoff flux from ice sheets to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_ice_pac",  vars%runoff_ice(3),dim1=dim_time,start=[ndat],count=[y],long_name="Pacific runoff flux from ice sheets to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_ice_ind",  vars%runoff_ice(4),dim1=dim_time,start=[ndat],count=[y],long_name="Indian ocean runoff flux from ice sheets to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_ice_so",   vars%runoff_ice(5),dim1=dim_time,start=[ndat],count=[y],long_name="Southern ocean runoff flux from ice sheets to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_lake_glob", vars%runoff_lake(1),dim1=dim_time,start=[ndat],count=[y],long_name="global runoff flux from lakes to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_lake_atl",  vars%runoff_lake(2),dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic runoff flux from lakes to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_lake_atlN30", vars%runoff_lake(6),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>30N) runoff flux from lakes to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_lake_atlN50", vars%runoff_lake(7),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>50N) runoff flux from lakes to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_lake_pac",  vars%runoff_lake(3),dim1=dim_time,start=[ndat],count=[y],long_name="Pacific runoff flux from lakes to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_lake_ind",  vars%runoff_lake(4),dim1=dim_time,start=[ndat],count=[y],long_name="Indian ocean runoff flux from lakes to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"runoff_lake_so",   vars%runoff_lake(5),dim1=dim_time,start=[ndat],count=[y],long_name="Southern ocean runoff flux from lakes to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"calving_glob", vars%calving(1),dim1=dim_time,start=[ndat],count=[y],long_name="global calving flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"calving_glob", vars%calving(1),dim1=dim_time,start=[ndat],count=[y],long_name="global calving flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"calving_atl",  vars%calving(2),dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic calving flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"calving_atlN30", vars%calving(6),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>30N) calving flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"calving_atlN50", vars%calving(7),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>50N) calving flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"calving_pac",  vars%calving(3),dim1=dim_time,start=[ndat],count=[y],long_name="Pacific calving flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"calving_ind",  vars%calving(4),dim1=dim_time,start=[ndat],count=[y],long_name="Indian ocean calving flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"calving_so",   vars%calving(5),dim1=dim_time,start=[ndat],count=[y],long_name="Southern ocean calving flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"bmelt_glob", vars%bmelt(1),dim1=dim_time,start=[ndat],count=[y],long_name="global ice basal melt flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"bmelt_atl",  vars%bmelt(2),dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic ice basal melt flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"bmelt_atlN30", vars%bmelt(6),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>30N) ice basal melt flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"bmelt_atlN50", vars%bmelt(7),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>50N) ice basal melt flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"bmelt_pac",  vars%bmelt(3),dim1=dim_time,start=[ndat],count=[y],long_name="Pacific ice basal flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"bmelt_ind",  vars%bmelt(4),dim1=dim_time,start=[ndat],count=[y],long_name="Indian ocean ice basal melt flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"bmelt_so",   vars%bmelt(5),dim1=dim_time,start=[ndat],count=[y],long_name="Southern ocean ice basal melt flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"icemelt_atl",    vars%icemelt(2),dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic ice-sheet melt flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"icemelt_atlN30", vars%icemelt(6),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>30) ice-sheet melt flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"icemelt_atlN50", vars%icemelt(7),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>50N) ice-sheet melt flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"icemelt_pac",    vars%icemelt(3),dim1=dim_time,start=[ndat],count=[y],long_name="Pacific ice-sheet melt flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"icemelt_ind",    vars%icemelt(4),dim1=dim_time,start=[ndat],count=[y],long_name="Indian ocean ice-sheet melt flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"icemelt_so",     vars%icemelt(5),dim1=dim_time,start=[ndat],count=[y],long_name="Southern ocean ice-sheet melt flux to ocean",units="Sv",ncid=ncid)
    call nc_write(fnm,"flx_glob_avg",vars%flx(1)*1.e15_wp/ocn_area_tot,dim1=dim_time,start=[ndat],count=[y],long_name="global average net heat flux to ocean",units="W/m2",ncid=ncid)
    call nc_write(fnm,"flx_glob",vars%flx(1),dim1=dim_time,start=[ndat],count=[y],long_name="global net heat flux to ocean",units="PW",ncid=ncid)
    call nc_write(fnm,"flx_atl", vars%flx(2),dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic heat flux to ocean",units="PW",ncid=ncid)
    call nc_write(fnm,"flx_atlN30", vars%flx(6),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>30N) heat flux to ocean",units="PW",ncid=ncid)
    call nc_write(fnm,"flx_atlN50", vars%flx(7),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>50N) heat flux to ocean",units="PW",ncid=ncid)
    call nc_write(fnm,"flx_atlN55", vars%flx(9),dim1=dim_time,start=[ndat],count=[y],long_name="North Atlantic (>55N) heat flux to ocean",units="PW",ncid=ncid)
    call nc_write(fnm,"flx_pac", vars%flx(3),dim1=dim_time,start=[ndat],count=[y],long_name="Pacific heat flux to ocean",units="PW",ncid=ncid)
    call nc_write(fnm,"flx_ind", vars%flx(4),dim1=dim_time,start=[ndat],count=[y],long_name="Indian ocean heat flux to ocean",units="PW",ncid=ncid)
    call nc_write(fnm,"flx_so",  vars%flx(5),dim1=dim_time,start=[ndat],count=[y],long_name="Southern ocean heat flux to ocean",units="PW",ncid=ncid)
    call nc_write(fnm,"flx_soS60", vars%flx(10),dim1=dim_time,start=[ndat],count=[y],long_name="Southern ocean (<60S) heat flux to ocean",units="PW",ncid=ncid)
    call nc_write(fnm,"drake",   vars%drake,dim1=dim_time,start=[ndat],count=[y],long_name="Drake passage throughflow (positive=east) (128 +/- 8)",units="SV",ncid=ncid)
    call nc_write(fnm,"bering",   vars%bering,dim1=dim_time,start=[ndat],count=[y],long_name="Bering Strait throughflow (positive=north) (0.8 +/- 0.2)",units="SV",ncid=ncid)
    call nc_write(fnm,"davis",   vars%davis,dim1=dim_time,start=[ndat],count=[y],long_name="Davis Strait throughflow (positive=north) (~-2)",units="SV",ncid=ncid)
    call nc_write(fnm,"fram",   vars%fram,dim1=dim_time,start=[ndat],count=[y],long_name="Fram Strait throughflow (positive=north)",units="SV",ncid=ncid)
    call nc_write(fnm,"denmark",   vars%denmark,dim1=dim_time,start=[ndat],count=[y],long_name="Denmark Strait throughflow (positive=north)",units="SV",ncid=ncid)
    call nc_write(fnm,"medi",   vars%medi,dim1=dim_time,start=[ndat],count=[y],long_name="Gibraltar Strait throughflow",units="SV",ncid=ncid)
    call nc_write(fnm,"indo",   vars%indo,dim1=dim_time,start=[ndat],count=[y],long_name="Indonesian Passage throughflow (-15 +/- 4)",units="SV",ncid=ncid)
    call nc_write(fnm,"agulhas",   vars%agulhas,dim1=dim_time,start=[ndat],count=[y],long_name="Flow from Indian into Atlantic around the tip of South Africa",units="SV",ncid=ncid)
    call nc_write(fnm,"fov",    vars%fov,dim1=dim_time,start=[ndat],count=[y],long_name="Net freshwater flux by the MOC into the Atlantic",units="Sv",ncid=ncid)
    call nc_write(fnm,"fovs",    vars%fovs,dim1=dim_time,start=[ndat],count=[y],long_name="Freshwater flux by the MOC into the Atlantic at its southern border",units="Sv",ncid=ncid)
    call nc_write(fnm,"fovn",    vars%fovn,dim1=dim_time,start=[ndat],count=[y],long_name="Freshwater flux by the MOC out of the Atlantic at its northern border",units="Sv",ncid=ncid)
    call nc_write(fnm,"faz",    vars%faz,dim1=dim_time,start=[ndat],count=[y],long_name="Net freshwater flux by the azonal circulation into the Atlantic",units="Sv",ncid=ncid)
    call nc_write(fnm,"fazs",    vars%fazs,dim1=dim_time,start=[ndat],count=[y],long_name="Freshwater flux by the azonal circulation into the Atlantic at its southern border",units="Sv",ncid=ncid)
    call nc_write(fnm,"fazn",    vars%fazn,dim1=dim_time,start=[ndat],count=[y],long_name="Freshwater flux by the azonal circulation out of the Atlantic at its northern border",units="Sv",ncid=ncid)

    call nc_write(fnm,"fw_lab",    vars%fw_lab,dim1=dim_time,start=[ndat],count=[y],long_name="Freshwater export through Labrador current",units="Sv",ncid=ncid)
    call nc_write(fnm,"buoy",    vars%buoy,dim1=dim_time,start=[ndat],count=[y],long_name="Global surface bouyancy flux",units="N",ncid=ncid)
    call nc_write(fnm,"buoyT",    vars%buoyT,dim1=dim_time,start=[ndat],count=[y],long_name="Global thermal surface bouyancy flux",units="N",ncid=ncid)
    call nc_write(fnm,"buoyS",    vars%buoyS,dim1=dim_time,start=[ndat],count=[y],long_name="Global haline surface bouyancy flux",units="N",ncid=ncid)
    call nc_write(fnm,"buoy_N",    vars%buoy_N,dim1=dim_time,start=[ndat],count=[y],long_name="Surface bouyancy flux north of 30N",units="N",ncid=ncid)
    call nc_write(fnm,"buoyT_N",    vars%buoyT_N,dim1=dim_time,start=[ndat],count=[y],long_name="Thermal surface bouyancy flux north of 30N ",units="N",ncid=ncid)
    call nc_write(fnm,"buoyS_N",    vars%buoyS_N,dim1=dim_time,start=[ndat],count=[y],long_name="Haline surface bouyancy flux north of 30N",units="N",ncid=ncid)
    call nc_write(fnm,"buoy_tr",    vars%buoy_tr,dim1=dim_time,start=[ndat],count=[y],long_name="Surface bouyancy flux between 30S and 30N",units="N",ncid=ncid)
    call nc_write(fnm,"buoyT_tr",    vars%buoyT_tr,dim1=dim_time,start=[ndat],count=[y],long_name="Thermal surface bouyancy flux between 30S and 30N",units="N",ncid=ncid)
    call nc_write(fnm,"buoyS_tr",    vars%buoyS_tr,dim1=dim_time,start=[ndat],count=[y],long_name="Haline surface bouyancy flux between 30S and 30N",units="N",ncid=ncid)
    call nc_write(fnm,"buoy_S",    vars%buoy_S,dim1=dim_time,start=[ndat],count=[y],long_name="Surface bouyancy flux south of 30S",units="N",ncid=ncid)
    call nc_write(fnm,"buoyT_S",    vars%buoyT_S,dim1=dim_time,start=[ndat],count=[y],long_name="Thermal surface bouyancy flux south of 30S",units="N",ncid=ncid)
    call nc_write(fnm,"buoyS_S",    vars%buoyS_S,dim1=dim_time,start=[ndat],count=[y],long_name="Haline surface bouyancy flux south of 30S",units="N",ncid=ncid)
    call nc_write(fnm,"buoyT_NA55",    vars%buoyT_NA(ilatv_buoy_sel),dim1=dim_time,start=[ndat],count=[y],long_name="Thermal component of surface bouyancy flux over North Atlantic",units="N",ncid=ncid)
    call nc_write(fnm,"buoyS_NA55",    vars%buoyS_NA(ilatv_buoy_sel),dim1=dim_time,start=[ndat],count=[y],long_name="Haline component of surface bouyancy flux over North Atlantic",units="N",ncid=ncid)
    call nc_write(fnm,"buoySw_NA55",    vars%buoySw_NA(ilatv_buoy_sel),dim1=dim_time,start=[ndat],count=[y],long_name="Buoyancy flux from net salinity transport into the North Atlantic",units="N",ncid=ncid)
    call nc_write(fnm,"buoySw_lab",    vars%buoySw_lab,dim1=dim_time,start=[ndat],count=[y],long_name="Buoyancy flux from net freshwater export through Labrador current",units="N",ncid=ncid)
    call nc_write(fnm,"buoySi_NA55",    vars%buoySi_NA(ilatv_buoy_sel),dim1=dim_time,start=[ndat],count=[y],long_name="Sea ice export component of surface bouyancy flux over North Atlantic",units="N",ncid=ncid)
    call nc_write(fnm,"buoy_NA55",    vars%buoy_NA(ilatv_buoy_sel),dim1=dim_time,start=[ndat],count=[y],long_name="Bouyancy flux over North Atlantic",units="N",ncid=ncid)
    do n=1,nlatv_buoy
      call nc_write(fnm,"buoyT_NA",    vars%buoyT_NA(n),dim1=dim_lat,dim2=dim_time,start=[n,ndat],count=[1,y],long_name="Thermal component of surface bouyancy flux over North Atlantic",units="N",ncid=ncid)
      call nc_write(fnm,"buoyS_NA",    vars%buoyS_NA(n),dim1=dim_lat,dim2=dim_time,start=[n,ndat],count=[1,y],long_name="Haline component of surface bouyancy flux over North Atlantic",units="N",ncid=ncid)
      call nc_write(fnm,"buoySw_NA",    vars%buoySw_NA(n),dim1=dim_lat,dim2=dim_time,start=[n,ndat],count=[1,y],long_name="Buoyancy flux from net salinity transport into the North Atlantic",units="N",ncid=ncid)
      call nc_write(fnm,"buoySi_NA",    vars%buoySi_NA(n),dim1=dim_lat,dim2=dim_time,start=[n,ndat],count=[1,y],long_name="Sea ice export component of surface bouyancy flux over North Atlantic",units="N",ncid=ncid)
      call nc_write(fnm,"buoy_NA",    vars%buoy_NA(n),dim1=dim_lat,dim2=dim_time,start=[n,ndat],count=[1,y],long_name="Bouyancy flux over North Atlantic",units="N",ncid=ncid)
    enddo

    call nc_write(fnm,"fw_bering",   vars%fw_bering,dim1=dim_time,start=[ndat],count=[y],long_name="Bering Strait freshwater transport relative to saln0",units="SV",ncid=ncid)
    call nc_write(fnm,"fw_davis",   vars%fw_davis,dim1=dim_time,start=[ndat],count=[y],long_name="Davis Strait freshwater transport relative to saln0",units="SV",ncid=ncid)
    call nc_write(fnm,"fw_fram",   vars%fw_fram,dim1=dim_time,start=[ndat],count=[y],long_name="Fram Strait freshwater transport relative to saln0",units="SV",ncid=ncid)
    call nc_write(fnm,"fw_denmark",   vars%fw_denmark,dim1=dim_time,start=[ndat],count=[y],long_name="Denmark Strait freshwater transport relative to saln0",units="SV",ncid=ncid)
    call nc_write(fnm,"drhoatl1",   vars%drhoatl1,dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic meridional density gradient at 750 m depth between 52.5N and 32.5S",units="kg/m3",ncid=ncid)
    call nc_write(fnm,"drhoatl2",   vars%drhoatl2,dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic meridional density gradient at 750 m depth between 40N-60N and 30S-60N",units="kg/m3",ncid=ncid)
    call nc_write(fnm,"drhoTatl2",   vars%drhoTatl2,dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic meridional density gradient due to temperature at 750 m depth between 40N-60N and 30S-60N",units="kg/m3",ncid=ncid)
    call nc_write(fnm,"drhoSatl2",   vars%drhoSatl2,dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic meridional density gradient due to salinity at 750 m depth between 40N-60N and 30S-60N",units="kg/m3",ncid=ncid)
    call nc_write(fnm,"drhoatl3",   vars%drhoatl3,dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic meridional density gradient at 750 m depth between 50N-70N and 30S-60N",units="kg/m3",ncid=ncid)
    call nc_write(fnm,"drhoTatl3",   vars%drhoTatl3,dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic meridional density gradient due to temperature at 750 m depth between 50N-70N and 30S-60N",units="kg/m3",ncid=ncid)
    call nc_write(fnm,"drhoSatl3",   vars%drhoSatl3,dim1=dim_time,start=[ndat],count=[y],long_name="Atlantic meridional density gradient due to salinity at 750 m depth between 50N-70N and 30S-60N",units="kg/m3",ncid=ncid)

    if (l_cfc) then
      call nc_write(fnm,"cfc11",     vars%cfc11,    dim1=dim_time,start=[ndat],count=[y],long_name="CFC11 integral",units="Mmol",ncid=ncid)
      call nc_write(fnm,"cfc12",     vars%cfc12,    dim1=dim_time,start=[ndat],count=[y],long_name="CFC12 integral",units="Mmol",ncid=ncid)
    endif

    call nc_write(fnm,"ncells",  vars%ncells,dim1=dim_time,start=[ndat],count=[y],long_name="Number of active ocean grid cells",units="1",ncid=ncid)
    call nc_write(fnm,"area",    vars%area,dim1=dim_time,start=[ndat],count=[y],long_name="Total surface ocean area",units="mln km2",ncid=ncid)
    call nc_write(fnm,"vol",     vars%vol,dim1=dim_time,start=[ndat],count=[y],long_name="Total ocean volume",units="mln km3",ncid=ncid)
    call nc_write(fnm,"shelf",   vars%shelf,dim1=dim_time,start=[ndat],count=[y],long_name="Total area of ocean shelf",units="mln km2",ncid=ncid)
    call nc_write(fnm,"sl_steric",   vars%sl_steric,dim1=dim_time,start=[ndat],count=[y],long_name="Steric sea level change relative to first simulation year (includes also contribution form land ice volume changes if ice sheets are interactive)",units="m",ncid=ncid)

    call nc_write(fnm,"mld_atlN50",   vars%mld_atlN50,dim1=dim_time,start=[ndat],count=[y],long_name="Maximum mixed layer depth in Atlantic >50N",units="m",ncid=ncid)
    call nc_write(fnm,"mld_lab",   vars%mld_lab,dim1=dim_time,start=[ndat],count=[y],long_name="Maximum mixed layer depth in Labrador Sea",units="m",ncid=ncid)
    call nc_write(fnm,"mld_irm",   vars%mld_irm,dim1=dim_time,start=[ndat],count=[y],long_name="Maximum mixed layer depth in Irminger Sea",units="m",ncid=ncid)
    call nc_write(fnm,"mld_gin",  vars%mld_gin,dim1=dim_time,start=[ndat],count=[y],long_name="Maximum mixed layer depth in the GIN seas",units="m",ncid=ncid)
    call nc_write(fnm,"mld_bkn",  vars%mld_bkn,dim1=dim_time,start=[ndat],count=[y],long_name="Maximum mixed layer depth in the Barents-Kara-Nansen seas",units="m",ncid=ncid)
    call nc_write(fnm,"mld_wedd", vars%mld_wedd,dim1=dim_time,start=[ndat],count=[y],long_name="Maximum mixed layer depth in the Weddel Sea",units="m",ncid=ncid)
    call nc_write(fnm,"mld_ross", vars%mld_ross,dim1=dim_time,start=[ndat],count=[y],long_name="Maximum mixed layer depth in the Ross Sea",units="m",ncid=ncid)
    call nc_write(fnm,"mld_so", vars%mld_so,dim1=dim_time,start=[ndat],count=[y],long_name="Maximum mixed layer depth in the Southern Ocean around Antarctica",units="m",ncid=ncid)

    call nc_write(fnm,"mldst_atlN50",   vars%mldst_atlN50,dim1=dim_time,start=[ndat],count=[y],long_name="Maximum mixed layer depth based on sigma-t in Atlantic <50N",units="m",ncid=ncid)
    call nc_write(fnm,"mldst_lab",   vars%mldst_lab,dim1=dim_time,start=[ndat],count=[y],long_name="Maximum mixed layer depth based on sigma-t in Labrador Sea",units="m",ncid=ncid)
    call nc_write(fnm,"mldst_irm",   vars%mldst_irm,dim1=dim_time,start=[ndat],count=[y],long_name="Maximum mixed layer depth based on sigma-t in Irminger Sea",units="m",ncid=ncid)
    call nc_write(fnm,"mldst_gin",  vars%mldst_gin,dim1=dim_time,start=[ndat],count=[y],long_name="Maximum mixed layer depth based on sigma-t in the GIN seas",units="m",ncid=ncid)
    call nc_write(fnm,"mldst_bkn",  vars%mldst_bkn,dim1=dim_time,start=[ndat],count=[y],long_name="Maximum mixed layer depth based on sigma-t in the Barents-Kara-Nansen seas",units="m",ncid=ncid)
    call nc_write(fnm,"mldst_wedd", vars%mldst_wedd,dim1=dim_time,start=[ndat],count=[y],long_name="Maximum mixed layer depth based on sigma-t in the Weddel Sea",units="m",ncid=ncid)
    call nc_write(fnm,"mldst_ross", vars%mldst_ross,dim1=dim_time,start=[ndat],count=[y],long_name="Maximum mixed layer depth based on sigma-t in the Ross Sea",units="m",ncid=ncid)
    call nc_write(fnm,"mldst_so", vars%mldst_so,dim1=dim_time,start=[ndat],count=[y],long_name="Maximum mixed layer depth based on sigma-t in the Southern Ocean around Antarctica",units="m",ncid=ncid)

    call nc_write(fnm,"pe_atlN",   vars%pe_atlN,dim1=dim_time,start=[ndat],count=[y],long_name="Potential energy released by convection in North Atlantic",units="J",ncid=ncid)
    call nc_write(fnm,"pe_atlN50",   vars%pe_atlN50,dim1=dim_time,start=[ndat],count=[y],long_name="Potential energy released by convection in Atlantic >50N",units="J",ncid=ncid)
    call nc_write(fnm,"pe_lab",   vars%pe_lab,dim1=dim_time,start=[ndat],count=[y],long_name="Potential energy released by convection in Labrador Sea",units="J",ncid=ncid)
    call nc_write(fnm,"pe_irm",   vars%pe_irm,dim1=dim_time,start=[ndat],count=[y],long_name="Potential energy released by convection in Irminger Sea",units="J",ncid=ncid)
    call nc_write(fnm,"pe_gin",  vars%pe_gin,dim1=dim_time,start=[ndat],count=[y],long_name="Potential energy released by convection in the GIN seas",units="J",ncid=ncid)
    call nc_write(fnm,"pe_bkn",  vars%pe_bkn,dim1=dim_time,start=[ndat],count=[y],long_name="Potential energy released by convection in the Barents-Kara-Nansen seas",units="J",ncid=ncid)
    call nc_write(fnm,"pe_wedd", vars%pe_wedd,dim1=dim_time,start=[ndat],count=[y],long_name="Potential energy released by convection in the Weddel Sea",units="J",ncid=ncid)
    call nc_write(fnm,"pe_ross", vars%pe_ross,dim1=dim_time,start=[ndat],count=[y],long_name="Potential energy released by convection in the Ross Sea",units="J",ncid=ncid)
    call nc_write(fnm,"pe_so",   vars%pe_so,dim1=dim_time,start=[ndat],count=[y],long_name="Potential energy released by convection in the Southern Ocean around Antarctica",units="J",ncid=ncid)

    call nc_write(fnm,"buoy_lab",   vars%buoy_lab,dim1=dim_time,start=[ndat],count=[y],long_name="Buoyancy flux over Labrador Sea",units="N",ncid=ncid)
    call nc_write(fnm,"buoy_irm",   vars%buoy_irm,dim1=dim_time,start=[ndat],count=[y],long_name="Buoyancy flux over Irminger Sea",units="N",ncid=ncid)
    call nc_write(fnm,"buoy_gin",  vars%buoy_gin,dim1=dim_time,start=[ndat],count=[y],long_name="Buoyancy flux over the GIN seas",units="N",ncid=ncid)
    call nc_write(fnm,"buoy_bkn",  vars%buoy_bkn,dim1=dim_time,start=[ndat],count=[y],long_name="Buoyancy flux over the Barents-Kara-Nansen seas",units="N",ncid=ncid)
    call nc_write(fnm,"buoy_wedd", vars%buoy_wedd,dim1=dim_time,start=[ndat],count=[y],long_name="Buoyancy flux over the Weddel Sea",units="N",ncid=ncid)
    call nc_write(fnm,"buoy_ross", vars%buoy_ross,dim1=dim_time,start=[ndat],count=[y],long_name="Buoyancy flux over the Ross Sea",units="N",ncid=ncid)
    call nc_write(fnm,"buoy_so",   vars%buoy_so,dim1=dim_time,start=[ndat],count=[y],long_name="Buoyancy flux over the Southern Ocean around Antarctica",units="N",ncid=ncid)
    call nc_write(fnm,"buoy_soS60",   vars%buoy_soS60,dim1=dim_time,start=[ndat],count=[y],long_name="Buoyancy flux over the Southern Ocean south of 60S",units="N",ncid=ncid)
    call nc_write(fnm,"buoyT_soS60",   vars%buoyT_soS60,dim1=dim_time,start=[ndat],count=[y],long_name="Thermal component of buoyancy flux over the Southern Ocean south of 60S",units="N",ncid=ncid)
    call nc_write(fnm,"buoyS_soS60",   vars%buoyS_soS60,dim1=dim_time,start=[ndat],count=[y],long_name="Haline component of buoyancy flux over the Southern Ocean south of 60S",units="N",ncid=ncid)

    call nc_write(fnm,"t_atlN50",   vars%t_atlN50,dim1=dim_time,start=[ndat],count=[y],  long_name="Average sea surface temperature in Atlantic >50N",units="degC",ncid=ncid)
    call nc_write(fnm,"t_lab",   vars%t_lab,dim1=dim_time,start=[ndat],count=[y],  long_name="Average sea surface temperature in Labrador Sea",units="degC",ncid=ncid)
    call nc_write(fnm,"t_irm",   vars%t_irm,dim1=dim_time,start=[ndat],count=[y],  long_name="Average sea surface temperature in Irminger Sea",units="degC",ncid=ncid)
    call nc_write(fnm,"t_gin",  vars%t_gin,dim1=dim_time,start=[ndat],count=[y], long_name="Average sea surface temperature in the GIN seas",units="degC",ncid=ncid)
    call nc_write(fnm,"t_bkn",  vars%t_bkn,dim1=dim_time,start=[ndat],count=[y], long_name="Average sea surface temperature in the Barents-Kara-Nansen seas",units="degC",ncid=ncid)
    call nc_write(fnm,"t_wedd", vars%t_wedd,dim1=dim_time,start=[ndat],count=[y],long_name="Average sea surface temperature in the Weddel Sea",units="degC",ncid=ncid)
    call nc_write(fnm,"t_ross", vars%t_ross,dim1=dim_time,start=[ndat],count=[y],long_name="Average sea surface temperature in the Ross Sea",units="degC",ncid=ncid)
    call nc_write(fnm,"t_so",   vars%t_so,dim1=dim_time,start=[ndat],count=[y],  long_name="Average sea surface temperature in the Southern Ocean around Antarctica",units="degC",ncid=ncid)
    call nc_write(fnm,"t_ibe",  vars%t_ibe,dim1=dim_time,start=[ndat],count=[y], long_name="Annual mean sea surface temperature at the Iberian margin",units="degC",ncid=ncid)

    call nc_write(fnm,"s_atlN50",   vars%s_atlN50,dim1=dim_time,start=[ndat],count=[y],  long_name="Average sea surface salinity in Atlantic >50N",units="degC",ncid=ncid)
    call nc_write(fnm,"s_lab",   vars%s_lab,dim1=dim_time,start=[ndat],count=[y],  long_name="Average sea surface salinity in Labrador Sea",units="degC",ncid=ncid)
    call nc_write(fnm,"s_irm",   vars%s_irm,dim1=dim_time,start=[ndat],count=[y],  long_name="Average sea surface salinity in Irminger Sea",units="degC",ncid=ncid)
    call nc_write(fnm,"s_gin",  vars%s_gin,dim1=dim_time,start=[ndat],count=[y], long_name="Average sea surface salinity in the GIN seas",units="degC",ncid=ncid)
    call nc_write(fnm,"s_bkn",  vars%s_bkn,dim1=dim_time,start=[ndat],count=[y], long_name="Average sea surface salinity in the Barents-Kara-Nansen seas",units="degC",ncid=ncid)
    call nc_write(fnm,"s_wedd", vars%s_wedd,dim1=dim_time,start=[ndat],count=[y],long_name="Average sea surface salinity in the Weddel Sea",units="degC",ncid=ncid)
    call nc_write(fnm,"s_ross", vars%s_ross,dim1=dim_time,start=[ndat],count=[y],long_name="Average sea surface salinity in the Ross Sea",units="degC",ncid=ncid)
    call nc_write(fnm,"s_so",   vars%s_so,dim1=dim_time,start=[ndat],count=[y],  long_name="Average sea surface salinity in the Southern Ocean around Antarctica",units="degC",ncid=ncid)

    call nc_close(ncid)


   return

  end subroutine ts_nc_write


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  o c n _ n c
  ! Purpose  :  Initialize ocean netcdf output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocn_nc(fnm)

    implicit none

    character (len=*) :: fnm
    integer :: ncid
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,dim_time, x=empty_time, units="years BP", unlimited=.TRUE.,ncid=ncid)
    call nc_write_dim(fnm,dim_lon,x=lon,axis="x",ncid=ncid)
    call nc_write_dim(fnm,dim_lat,x=lat,axis="y",ncid=ncid)
    call nc_write_dim(fnm,"lonu",x=lonu(2:maxi+1),axis="y",ncid=ncid)
    call nc_write_dim(fnm,"latv",x=latv(2:maxj+1),axis="y",ncid=ncid)
    call nc_write_dim(fnm,"latv1",x=latv(1:maxj+1),axis="y",ncid=ncid)
    call nc_write_dim(fnm,"lev",x=-zro(maxk:1:-1),units="m",axis="z",ncid=ncid)
    call nc_write_dim(fnm,"levw",x=-zw(maxk:1:-1),units="m",axis="z",ncid=ncid)
    call nc_write_dim(fnm,dim_isles,x=1._wp,dx=1._wp,nx=maxisles,ncid=ncid)
    call nc_write_dim(fnm,dim_type,x=1._wp,dx=1._wp,nx=5,units="[tot,adv,diff,over,gyre]",ncid=ncid)
    call nc_write_dim(fnm,dim_month,x=1._wp,dx=1._wp,nx=13,units="months",ncid=ncid)
    call nc_write(fnm,"dx", dx, dims=[dim_lat],start=[1],count=[maxj],long_name="longitude grid distance at cell centers",units="m",ncid=ncid)
    call nc_write(fnm,"dxv", dxv, dims=["latv1"],start=[1],count=[maxj+1],long_name="longitude grid distance at cell edges on v-grid",units="m",ncid=ncid)
    call nc_write(fnm,"dy", dy, long_name="latitudinal grid distance at cell centers and cell edges",units="m",ncid=ncid)
    call nc_write(fnm,"dz", dz(maxk:1:-1), dims=["lev"],start=[1],count=[maxk],long_name="thickness of grid cells",units="m",ncid=ncid)
    call nc_write(fnm,"dzw", dza(maxk:1:-1), dims=["levw"],start=[1],count=[maxk],long_name="thickness between layer centers",units="m",ncid=ncid)
    call nc_close(ncid)

   return

  end subroutine ocn_nc


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  o c n _ n c _ w r i t e
  ! Purpose  :  Output of ocean netcdf
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocn_nc_write(fnm,ncid,vars,ndat,nout)

    implicit none

    type(o_out) :: vars

    character (len=*) :: fnm
    integer :: ndat, nout, ncid

    if (ndat.eq.13) then
    call nc_write(fnm,"mask_ocn",    vars%mask_ocn,  dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[maxi,maxj,1],long_name="surface ocean mask",units="/",missing_value=int(missing_value),ncid=ncid)
    call nc_write(fnm,"f_ocn",    sngl(vars%f_ocn),  dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[maxi,maxj,1],long_name="surface ocean fraction",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"area",    sngl(vars%area), dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[maxi,maxj,1],long_name="surface ocean area",units="m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"k1",    vars%k1,  dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[maxi,maxj,1],long_name="index of first (bottom) layer",units="/",missing_value=int(missing_value),ncid=ncid)
    call nc_write(fnm,"topo",    sngl(vars%topo),  dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[maxi,maxj,1],long_name="topography",units="m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"bathy",    sngl(vars%bathy),  dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[maxi,maxj,1],long_name="bathymetry",units="m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"vol",    sngl(vars%vol),  dims=[dim_lon,dim_lat,"lev",dim_time],start=[1,1,1,nout],count=[maxi,maxj,maxk,1],long_name="ocean volume",units="m3",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"drag",    sngl(vars%drag),  dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[maxi,maxj,1],long_name="drag",units="s^-1",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"drag_bcl",    sngl(vars%drag_bcl),  dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[maxi,maxj,1],long_name="drag",units="s^-1",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"map_isles",    vars%map_isles,  dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[maxi,maxj,1],long_name="island map",units="/",missing_value=int(missing_value),ncid=ncid)
    call nc_write(fnm,"map_edge",    vars%map_edge,  dims=[dim_lon,dim_lat,dim_isles,dim_time],start=[1,1,1,nout],count=[maxi,maxj,maxisles,1],long_name="island edges",units="/",missing_value=int(missing_value),ncid=ncid)
    call nc_write(fnm,"mldmax",    sngl(vars%mldmax),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[maxi,maxj,1],long_name="maximum mixed layer depth from mixed layer scheme",units="m",missing_value=missing_value,ncid=ncid)
    endif
    call nc_write(fnm,"t",        sngl(vars%t),  dims=[dim_lon,dim_lat,"lev",dim_month,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="potential temperature",units="C",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"sst", sngl(vars%sst),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="sea surface temperature",units="C",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"tbot", sngl(vars%tbot),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="bottom ocean potential temperature",units="C",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"t_atl",    sngl(vars%t_atl),  dims=[dim_lat,"lev",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj,maxk,1,1],long_name="Zonal mean Atlantic potential temperature",units="C",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"t_pac",    sngl(vars%t_pac),  dims=[dim_lat,"lev",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj,maxk,1,1],long_name="Zonal mean Pacific potential temperature",units="C",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"t_ind",    sngl(vars%t_ind),  dims=[dim_lat,"lev",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj,maxk,1,1],long_name="Zonal mean Indian Ocean potential temperature",units="C",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"t_so",    sngl(vars%t_so),  dims=[dim_lat,"lev",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj,maxk,1,1],long_name="Zonal mean Southern Ocean potential temperature",units="C",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"s",        sngl(vars%s),  dims=[dim_lon,dim_lat,"lev",dim_month,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="salinity",units="psu",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"sss", sngl(vars%sss),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="sea surface salinity",units="psu",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"sbot", sngl(vars%sbot),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="bottom ocean salinity",units="psu",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"s_atl",    sngl(vars%s_atl),  dims=[dim_lat,"lev",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj,maxk,1,1],long_name="Zonal mean Atlantic salinity",units="psu",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"s_pac",    sngl(vars%s_pac),  dims=[dim_lat,"lev",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj,maxk,1,1],long_name="Zonal mean Pacific salinity",units="psu",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"s_ind",    sngl(vars%s_ind),  dims=[dim_lat,"lev",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj,maxk,1,1],long_name="Zonal mean Indian Ocean salinity",units="psu",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"s_so",    sngl(vars%s_so),  dims=[dim_lat,"lev",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj,maxk,1,1],long_name="Zonal mean Southern Ocean salinity",units="psu",missing_value=missing_value,ncid=ncid)
    if (age_tracer .and. ndat.eq.13) then
    call nc_write(fnm,"age",      sngl(vars%age),dims=[dim_lon,dim_lat,"lev",dim_time],start=[1,1,1,nout],count=[maxi,maxj,maxk,1],long_name="age",units="years",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"a_atl",    sngl(vars%a_atl),  dims=[dim_lat,"lev",dim_time],start=[1,1,nout],count=[maxj,maxk,1],long_name="Zonal mean Atlantic age",units="years",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"a_pac",    sngl(vars%a_pac),  dims=[dim_lat,"lev",dim_time],start=[1,1,nout],count=[maxj,maxk,1],long_name="Zonal mean Pacific age",units="years",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"a_ind",    sngl(vars%a_ind),  dims=[dim_lat,"lev",dim_time],start=[1,1,nout],count=[maxj,maxk,1],long_name="Zonal mean Indian Ocean age",units="years",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"a_so",     sngl(vars%a_so),   dims=[dim_lat,"lev",dim_time],start=[1,1,nout],count=[maxj,maxk,1],long_name="Zonal mean Southern Ocean age",units="years",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"dT_dt_adv", sngl(vars%dt_dt_adv),dims=[dim_lon,dim_lat,"lev",dim_time],start=[1,1,1,nout],count=[maxi,maxj,maxk,1],long_name="Temperature tendency due to advection",units="degC/year",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"dT_dt_diff", sngl(vars%dt_dt_diff),dims=[dim_lon,dim_lat,"lev",dim_time],start=[1,1,1,nout],count=[maxi,maxj,maxk,1],long_name="Temperature tendency due to diffusion",units="degC/year",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"dT_dt", sngl(vars%dt_dt_adv+vars%dt_dt_diff),dims=[dim_lon,dim_lat,"lev",dim_time],start=[1,1,1,nout],count=[maxi,maxj,maxk,1],long_name="Temperature tendency due to advection and diffusion",units="degC/year",missing_value=missing_value*2,ncid=ncid)
    call nc_write(fnm,"dT_dt_flxsur", sngl(vars%dt_dt_flxsur),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[maxi,maxj,1],long_name="Sea surface temperature tendency due to surface heat flux",units="degC/year",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"dS_dt_adv", sngl(vars%ds_dt_adv),dims=[dim_lon,dim_lat,"lev",dim_time],start=[1,1,1,nout],count=[maxi,maxj,maxk,1],long_name="Salinity tendency due to advection",units="psu/year",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"dS_dt_diff", sngl(vars%ds_dt_diff),dims=[dim_lon,dim_lat,"lev",dim_time],start=[1,1,1,nout],count=[maxi,maxj,maxk,1],long_name="Salinity tendency due to diffusion",units="psu/year",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"dS_dt", sngl(vars%ds_dt_adv+vars%ds_dt_diff),dims=[dim_lon,dim_lat,"lev",dim_time],start=[1,1,1,nout],count=[maxi,maxj,maxk,1],long_name="Salinity tendency due to advection and diffusion",units="psu/year",missing_value=missing_value*2,ncid=ncid)
    call nc_write(fnm,"dS_dt_flxsur", sngl(vars%ds_dt_flxsur),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[maxi,maxj,1],long_name="Sea surface salinity tendency due to surface freshwater flux",units="psu/year",missing_value=missing_value,ncid=ncid)
    endif
    if (dye_tracer .and. ndat.eq.13) then
    call nc_write(fnm,"dye",      sngl(vars%dye),dims=[dim_lon,dim_lat,"lev",dim_time],start=[1,1,1,nout],count=[maxi,maxj,maxk,1],long_name="dye tracer",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"d_atl",    sngl(vars%d_atl),  dims=[dim_lat,"lev",dim_time],start=[1,1,nout],count=[maxj,maxk,1],long_name="Zonal mean Atlantic dye",units="years",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"d_pac",    sngl(vars%d_pac),  dims=[dim_lat,"lev",dim_time],start=[1,1,nout],count=[maxj,maxk,1],long_name="Zonal mean Pacific dye",units="years",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"d_ind",    sngl(vars%d_ind),  dims=[dim_lat,"lev",dim_time],start=[1,1,nout],count=[maxj,maxk,1],long_name="Zonal mean Indian Ocean dye",units="years",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"d_so",     sngl(vars%d_so),   dims=[dim_lat,"lev",dim_time],start=[1,1,nout],count=[maxj,maxk,1],long_name="Zonal mean Southern Ocean dye",units="years",missing_value=missing_value,ncid=ncid)
    endif
    if (cons_tracer .and. ndat.eq.13) then
    call nc_write(fnm,"cons",      sngl(vars%cons),dims=[dim_lon,dim_lat,"lev",dim_time],start=[1,1,1,nout],count=[maxi,maxj,maxk,1],long_name="conservative tracer",units="/",missing_value=missing_value,ncid=ncid)
    endif
    if (l_cfc .and. ndat.eq.13) then
    call nc_write(fnm,"cfc11", sngl(vars%cfc11),dims=[dim_lon,dim_lat,"lev",dim_time],start=[1,1,1,nout],count=[maxi,maxj,maxk,1],long_name="CFC11",units="pmol/l",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"cfc12", sngl(vars%cfc12),dims=[dim_lon,dim_lat,"lev",dim_time],start=[1,1,1,nout],count=[maxi,maxj,maxk,1],long_name="CFC12",units="pmol/l",missing_value=missing_value,ncid=ncid)
    endif
    call nc_write(fnm,"rho",      sngl(vars%rho),dims=[dim_lon,dim_lat,"lev",dim_month,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="in-situ density",units="kg/m3",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"drho_0_1000", sngl(vars%drho_0_1000),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="density difference between the surface and 1000 m depth",units="kg/m3",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"drho_0_3000", sngl(vars%drho_0_3000),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="density difference between the surface and 3000 m depth",units="kg/m3",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"rho2",      sngl(vars%rho2),dims=[dim_lon,dim_lat,"lev",dim_month,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="density",units="kg/m3",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"rho_atl",    sngl(vars%rho_atl),  dims=[dim_lat,"lev",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj,maxk,1,1],long_name="Zonal mean Atlantic in-situ density",units="kg/m3",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"rho_pac",    sngl(vars%rho_pac),  dims=[dim_lat,"lev",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj,maxk,1,1],long_name="Zonal mean Pacific in-situ density",units="kg/m3",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"rho_ind",    sngl(vars%rho_ind),  dims=[dim_lat,"lev",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj,maxk,1,1],long_name="Zonal mean Indian Ocean in-situ density",units="kg/m3",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"rho_so",    sngl(vars%rho_so),  dims=[dim_lat,"lev",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj,maxk,1,1],long_name="Zonal mean Southern Ocean in-situ density",units="kg/m3",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"taux",       sngl(vars%taux), dims=["lonu      ",dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="zonal wind/sea ice stress",units="N/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"tauy",       sngl(vars%tauy), dims=[dim_lon,"latv",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="meridional wind/sea ice stress",units="N/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"u",        sngl(vars%u),  dims=["lonu      ",dim_lat,"lev       ",dim_month,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="zonal velocity",units="m/s",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"v",        sngl(vars%v),  dims=[dim_lon,"latv","lev",dim_month,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="meridional velocity",units="m/s",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"w",        sngl(vars%w),  dims=[dim_lon,dim_lat,"levw",dim_month,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="vertical velocity",units="m/s",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ub",       sngl(vars%ub), dims=["lonu      ",dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="barotropic zonal velocity",units="m/s",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"vb",       sngl(vars%vb), dims=[dim_lon,"latv      ",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="barotropic meridional velocity",units="m/s",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ubisl",       sngl(vars%ubisl), dims=["lonu      ",dim_lat,dim_isles,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxisles,1,1],long_name="islands barotropic zonal velocity",units="m/s",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"vbisl",       sngl(vars%vbisl), dims=[dim_lon,"latv      ",dim_isles,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxisles,1,1],long_name="islands barotropic meridional velocity",units="m/s",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"psi",      sngl(vars%psi),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="barotropic streamfunction",units="Sv",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"opsi",     sngl(vars%opsi),dims=["latv1     ","levw      ",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj+1,maxk,1,1],long_name="global overturning circulation",units="Sv",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"opsi_a",   sngl(vars%opsia),dims=["latv1     ","levw      ",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj+1,maxk,1,1],long_name="Atlantic overturning circulation",units="Sv",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"opsi_p",   sngl(vars%opsip),dims=["latv1     ","levw      ",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj+1,maxk,1,1],long_name="Pacific overturning circulation",units="Sv",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"opsi_i",   sngl(vars%opsii),dims=["latv1     ","levw      ",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj+1,maxk,1,1],long_name="Indian overturning circulation",units="Sv",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"flx",      sngl(vars%flx),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="net ocean heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fw",       sngl(vars%fw),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="net ocean freshwater flux",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"vsf",      sngl(vars%vsf),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="virtual salinity flux",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"p_e_sic", sngl(vars%p_e_sic),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="net ocean P-E+sea ice fluxes",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"runoff",   sngl(vars%runoff),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="runoff to ocean",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"runoffSv", sngl(vars%runoffSv),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="runoff to ocean in Sv",units="Sv",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"runoffSv_ice", sngl(vars%runoffSv_ice),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="runoff from ice sheets to ocean in Sv",units="Sv",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"calving",  sngl(vars%calving),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="calving to ocean",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"bmelt",    sngl(vars%bmelt),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="basal melt freshwater flux to ocean",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    if (l_hosing) then
    call nc_write(fnm,"fw_hosing",       sngl(vars%fw_hosing),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="freshwater hosing flux",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    endif
    if (l_flux_adj_atl .or. l_flux_adj_ant .or. l_flux_adj_pac) then
    call nc_write(fnm,"fw_flux_adj",       sngl(vars%fw_flux_adj),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="freshwater flux adjustment",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    endif
    if (l_noise_fw) then
    call nc_write(fnm,"fw_noise",       sngl(vars%fw_noise),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="freshwater flux noise",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    endif
    if (l_noise_flx) then
    call nc_write(fnm,"flx_noise",       sngl(vars%flx_noise),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="heat flux noise",units="W/m2",missing_value=missing_value,ncid=ncid)
    endif
    call nc_write(fnm,"buoy",       sngl(vars%buoy),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="surface buoyancy flux",units="N/m2/s",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"buoyT",       sngl(vars%buoyT),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="thermal component of surface buoyancy flux",units="N/m2/s",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"buoyS",       sngl(vars%buoyS),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="haline component of surface buoyancy flux",units="N/m2/s",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"hft",      sngl(vars%hft),dims=[dim_type,"latv",dim_month,dim_time],start=[1,1,ndat,nout],count=[5,maxj,1,1],long_name="Global ocean northward heat transport",units="PW",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"hfp",      sngl(vars%hfp),dims=[dim_type,"latv",dim_month,dim_time],start=[1,1,ndat,nout],count=[5,maxj,1,1],long_name="Indo-Pacific ocean northward heat transport",units="PW",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"hfa",      sngl(vars%hfa),dims=[dim_type,"latv",dim_month,dim_time],start=[1,1,ndat,nout],count=[5,maxj,1,1],long_name="Atlantic ocean northward heat transport",units="PW",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"hftz", sngl(vars%hftz),dims=["latv","lev",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj,maxk,1,1],long_name="Global ocean northward heat transport",units="PW/m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"hfpz", sngl(vars%hfpz),dims=["latv","lev",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj,maxk,1,1],long_name="Indo-Pacific ocean northward heat transport",units="PW/m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"hfaz", sngl(vars%hfaz),dims=["latv","lev",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj,maxk,1,1],long_name="Atlantic ocean northward heat transport",units="PW/m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fwt",      sngl(vars%fwt),dims=[dim_type,"latv",dim_month,dim_time],start=[1,1,ndat,nout],count=[5,maxj,1,1],long_name="Global ocean northward freshwater transport",units="Sv",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fwp",      sngl(vars%fwp),dims=[dim_type,"latv",dim_month,dim_time],start=[1,1,ndat,nout],count=[5,maxj,1,1],long_name="Indo-Pacific ocean northward freshwater transport",units="Sv",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fwa",      sngl(vars%fwa),dims=[dim_type,"latv",dim_month,dim_time],start=[1,1,ndat,nout],count=[5,maxj,1,1],long_name="Atlantic ocean northward freshwater transport",units="Sv",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fwtz", sngl(vars%fwtz),dims=["latv","lev",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj,maxk,1,1],long_name="Global ocean northward freshwater transport",units="Sv/m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fwpz", sngl(vars%fwpz),dims=["latv","lev",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj,maxk,1,1],long_name="Indo-Pacific ocean northward freshwater transport",units="Sv/m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fwaz", sngl(vars%fwaz),dims=["latv","lev",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxj,maxk,1,1],long_name="Atlantic ocean northward freshwater transport",units="Sv/m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fayti",    sngl(vars%fayti),dims=[dim_lon,"latv",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="vertically integrated northward advective heat transport",units="W",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fdyti",    sngl(vars%fdyti),dims=[dim_lon,"latv",dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="vertically integrated northward diffusive heat transport",units="W",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"nconv",     sngl(vars%nconv),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="number of mixed layers",units="\",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"kven",      sngl(vars%kven),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="maximum integer surface water ventilation",units="PW",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"dconv",     sngl(vars%dconv),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="maximum depth of convection",units="m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"dven",      sngl(vars%dven),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="surface water ventilation depth",units="m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"mld",       sngl(vars%mld),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="mixed layer depth from mixed layer scheme",units="m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"mldst",     sngl(vars%mldst),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="mixed layer depth from sigma-t criterion",units="m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ke_tau",  sngl(vars%ke_tau),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="kinetic energy input into the ocean by wind stress",units="mW/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ssh",      sngl(vars%ssh),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="elevation of the free surface",units="m",missing_value=missing_value,ncid=ncid)

    if (l_output_extended) then
    call nc_write(fnm,"fdx",        sngl(vars%fdx),  dims=["lonu      ",dim_lat,"lev       ",dim_month,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="zonal diffusive tracer volume flux",units="m3*K",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fdy",        sngl(vars%fdy),  dims=[dim_lon,"latv      ","lev       ",dim_month,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="meridional diffusive tracer volume flux",units="m3*K",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fdz",        sngl(vars%fdz),  dims=[dim_lon,dim_lat,"levw      ",dim_month,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="vertical diffusive tracer volume flux",units="m3*K",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"drho_dx",      sngl(vars%drho_dx),dims=[dim_lon,dim_lat,"lev",dim_month,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="x density gradient",units="kg/m4",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"drho_dy",      sngl(vars%drho_dy),dims=[dim_lon,dim_lat,"lev",dim_month,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="y density gradient",units="kg/m4",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"drho_dz",      sngl(vars%drho_dz),dims=[dim_lon,dim_lat,"lev",dim_month,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="z density gradient",units="kg/m4",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Ri",           sngl(vars%Ri),dims=[dim_lon,dim_lat,"lev",dim_month,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="Richardson number",units="1",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"slope_x",      sngl(abs(vars%drho_dx/min(vars%drho_dz,-1.e-6))),dims=[dim_lon,dim_lat,"lev",dim_month,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="x isopycnal slope",units="m/m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"slope_y",      sngl(abs(vars%drho_dy/min(vars%drho_dz,-1.e-6))),dims=[dim_lon,dim_lat,"lev",dim_month,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="y isopycnal slope",units="m/m",missing_value=missing_value,ncid=ncid)
    endif
    if (l_output_extended .or. l_diff_dia_strat) then
      call nc_write(fnm,"diffdia",      sngl(vars%diffdia),dims=[dim_lon,dim_lat,"lev",dim_month,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="diapycnal diffusivity",units="m2/s",missing_value=missing_value,ncid=ncid)
    endif

   return

  end subroutine ocn_nc_write


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  o c n _ d a i l y _ n c
  ! Purpose  :  Initialize ocean netcdf output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocn_daily_nc(fnm)

    implicit none

    character (len=*) :: fnm
    integer :: ncid
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,dim_time, x=empty_time, units="years BP", unlimited=.TRUE.,ncid=ncid)
    call nc_write_dim(fnm,dim_lon,x=lon,axis="x",ncid=ncid)
    call nc_write_dim(fnm,dim_lat,x=lat,axis="y",ncid=ncid)
    call nc_write_dim(fnm,"lev",x=zro,units="m",axis="z",ncid=ncid)
    call nc_write_dim(fnm,dim_day,x=1._wp,dx=1._wp,nx=nday_year,units="doy",ncid=ncid)
    call nc_close(ncid)

   return

  end subroutine ocn_daily_nc


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  o c n _ d a i l y _ n c _ w r i t e
  ! Purpose  :  Output of ocean netcdf
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocn_daily_nc_write(fnm,ncid,vars,ndat,nout)

    implicit none

    type(o_out) :: vars

    character (len=*) :: fnm
    integer :: ndat, nout, ncid

    call nc_write(fnm,"mld",      sngl(vars%mld),dims=[dim_lon,dim_lat,dim_day,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="mixed layer depth",units="m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"t",        sngl(vars%t),  dims=[dim_lon,dim_lat,"lev",dim_day,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="temperature",units="C",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"s",        sngl(vars%s),  dims=[dim_lon,dim_lat,"lev",dim_day,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="salinity",units="psu",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"rho",      sngl(vars%rho),  dims=[dim_lon,dim_lat,"lev",dim_day,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="density",units="kg/m3",missing_value=missing_value,ncid=ncid)
    if (cons_tracer) call nc_write(fnm,"cons",        sngl(vars%cons),  dims=[dim_lon,dim_lat,"lev",dim_day,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="conservative tracer",units="1",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"u",        sngl(vars%u),  dims=[dim_lon,dim_lat,"lev",dim_day,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="zonal velocity",units="m/s",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"v",        sngl(vars%v),  dims=[dim_lon,dim_lat,"lev",dim_day,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="meridional velocity",units="m/s",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"w",        sngl(vars%w),  dims=[dim_lon,dim_lat,"lev",dim_day,dim_time],start=[1,1,1,ndat,nout],count=[maxi,maxj,maxk,1,1],long_name="vertical velocity",units="m/s",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ub",       sngl(vars%ub),  dims=[dim_lon,dim_lat,dim_day,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="zonal barotropic velocity",units="m/s",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"vb",       sngl(vars%vb),  dims=[dim_lon,dim_lat,dim_day,dim_time],start=[1,1,ndat,nout],count=[maxi,maxj,1,1],long_name="meridional barotropic velocity",units="m/s",missing_value=missing_value,ncid=ncid)

   return

  end subroutine ocn_daily_nc_write


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  o c n _ a v e
  ! Purpose  :  Average (or sum) the ocn fields
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocn_ave(d,ave)

    implicit none

    type(o_out) :: d(:), ave

    integer :: k, n
    real(wp) :: div

    n = size(d)
    div = real(n,wp)

    ! Set all values to zero
    ave%t     = 0._wp
    ave%sst   = 0._wp
    ave%tbot  = 0._wp
    ave%t_atl = 0._wp
    ave%t_pac = 0._wp
    ave%t_ind = 0._wp
    ave%t_so  = 0._wp
    ave%s     = 0._wp
    ave%sss   = 0._wp
    ave%sbot  = 0._wp
    ave%s_atl = 0._wp
    ave%s_pac = 0._wp
    ave%s_ind = 0._wp
    ave%s_so  = 0._wp
    ave%dt_dt_adv  = 0._wp
    ave%dt_dt_diff = 0._wp
    ave%dt_dt_flxsur = 0._wp
    ave%ds_dt_adv  = 0._wp
    ave%ds_dt_diff = 0._wp
    ave%ds_dt_flxsur = 0._wp
    ave%diffdia  = 0._wp
    ave%rho   = 0._wp
    ave%drho_0_1000 = 0._wp
    ave%drho_0_3000 = 0._wp
    ave%rho2  = 0._wp
    ave%rho_atl = 0._wp
    ave%rho_pac = 0._wp
    ave%rho_ind = 0._wp
    ave%rho_so  = 0._wp
    ave%drho_dx = 0._wp
    ave%drho_dy = 0._wp
    ave%drho_dz = 0._wp
    ave%Ri = 0._wp
    ave%taux  = 0._wp
    ave%tauy  = 0._wp
    ave%u     = 0._wp
    ave%v     = 0._wp
    ave%w     = 0._wp
    ave%ub    = 0._wp
    ave%vb    = 0._wp
    ave%ubisl = 0._wp
    ave%vbisl = 0._wp
    ave%psi   = 0._wp
    ave%opsi  = 0._wp
    ave%opsia = 0._wp
    ave%opsip = 0._wp
    ave%opsii = 0._wp
    ave%fdx   = 0._wp
    ave%fdy   = 0._wp
    ave%fdz   = 0._wp
    ave%flx   = 0._wp
    ave%fw    = 0._wp
    ave%vsf   = 0._wp
    ave%p_e_sic   = 0._wp
    ave%runoff= 0._wp
    ave%runoffSv= 0._wp
    ave%runoffSv_ice= 0._wp
    ave%calving= 0._wp
    ave%bmelt = 0._wp
    ave%fw_hosing = 0._wp
    ave%fw_flux_adj = 0._wp
    ave%fw_noise = 0._wp
    ave%flx_noise = 0._wp
    ave%buoy  = 0._wp
    ave%buoyS = 0._wp
    ave%buoyT = 0._wp
    ave%hft   = 0._wp
    ave%hfp   = 0._wp
    ave%hfa   = 0._wp
    ave%hftz   = 0._wp
    ave%hfpz   = 0._wp
    ave%hfaz   = 0._wp
    ave%fwt   = 0._wp
    ave%fwp   = 0._wp
    ave%fwa   = 0._wp
    ave%fwtz   = 0._wp
    ave%fwpz   = 0._wp
    ave%fwaz   = 0._wp
    ave%fayti = 0._wp
    ave%fdyti = 0._wp
    ave%nconv = 0._wp
    ave%kven  = 0._wp
    ave%dconv = 0._wp
    ave%dven  = 0._wp
    ave%mld   = 0._wp
    ave%mldst = 0._wp
    ave%ke_tau   = 0._wp
    ave%ssh = 0._wp

    ! Loop over the time indices to sum up and average (if necessary)
    do k = 1, n
       ave%t       = ave%t        + d(k)%t        / div
       ave%sst     = ave%sst      + d(k)%sst      / div
       ave%tbot    = ave%tbot     + d(k)%tbot     / div
       ave%t_atl   = ave%t_atl    + d(k)%t_atl    / div
       ave%t_pac   = ave%t_pac    + d(k)%t_pac    / div
       ave%t_ind   = ave%t_ind    + d(k)%t_ind    / div
       ave%t_so    = ave%t_so     + d(k)%t_so     / div
       ave%s       = ave%s        + d(k)%s        / div
       ave%sss     = ave%sss      + d(k)%sss      / div
       ave%sbot    = ave%sbot     + d(k)%sbot     / div
       ave%s_atl   = ave%s_atl    + d(k)%s_atl    / div
       ave%s_pac   = ave%s_pac    + d(k)%s_pac    / div
       ave%s_ind   = ave%s_ind    + d(k)%s_ind    / div
       ave%s_so    = ave%s_so     + d(k)%s_so     / div
       ave%dt_dt_adv  = ave%dt_dt_adv  + d(k)%dt_dt_adv  / div
       ave%dt_dt_diff = ave%dt_dt_diff + d(k)%dt_dt_diff / div
       ave%dt_dt_flxsur = ave%dt_dt_flxsur + d(k)%dt_dt_flxsur / div
       ave%ds_dt_adv  = ave%ds_dt_adv  + d(k)%ds_dt_adv  / div
       ave%ds_dt_diff = ave%ds_dt_diff + d(k)%ds_dt_diff / div
       ave%ds_dt_flxsur = ave%ds_dt_flxsur + d(k)%ds_dt_flxsur / div
       ave%diffdia = ave%diffdia  + d(k)%diffdia  / div
       ave%drho_dx = ave%drho_dx  + d(k)%drho_dx  / div
       ave%drho_dy = ave%drho_dy  + d(k)%drho_dy  / div
       ave%drho_dz = ave%drho_dz  + d(k)%drho_dz  / div
       ave%Ri      = ave%Ri       + d(k)%Ri       / div
       ave%rho     = ave%rho      + d(k)%rho      / div
       ave%drho_0_1000 = ave%drho_0_1000 + d(k)%drho_0_1000 / div
       ave%drho_0_3000 = ave%drho_0_3000 + d(k)%drho_0_3000 / div
       ave%rho2    = ave%rho2     + d(k)%rho2     / div
       ave%rho_atl = ave%rho_atl  + d(k)%rho_atl  / div
       ave%rho_pac = ave%rho_pac  + d(k)%rho_pac  / div
       ave%rho_ind = ave%rho_ind  + d(k)%rho_ind  / div
       ave%rho_so  = ave%rho_so   + d(k)%rho_so   / div
       ave%taux    = ave%taux     + d(k)%taux     / div
       ave%tauy    = ave%tauy     + d(k)%tauy     / div
       ave%u       = ave%u        + d(k)%u        / div
       ave%v       = ave%v        + d(k)%v        / div
       ave%w       = ave%w        + d(k)%w        / div
       ave%ub      = ave%ub       + d(k)%ub       / div
       ave%vb      = ave%vb       + d(k)%vb       / div
       ave%ubisl   = ave%ubisl    + d(k)%ubisl    / div
       ave%vbisl   = ave%vbisl    + d(k)%vbisl    / div
       ave%psi     = ave%psi      + d(k)%psi      / div
       ave%opsi    = ave%opsi     + d(k)%opsi     / div
       ave%opsia   = ave%opsia    + d(k)%opsia    / div
       ave%opsip   = ave%opsip    + d(k)%opsip    / div
       ave%opsii   = ave%opsii    + d(k)%opsii    / div
       ave%fdx     = ave%fdx      + d(k)%fdx      / div
       ave%fdy     = ave%fdy      + d(k)%fdy      / div
       ave%fdz     = ave%fdz      + d(k)%fdz      / div
       ave%flx     = ave%flx      + d(k)%flx      / div
       ave%fw      = ave%fw       + d(k)%fw       / div
       ave%vsf     = ave%vsf      + d(k)%vsf      / div
       ave%p_e_sic     = ave%p_e_sic      + d(k)%p_e_sic      / div
       ave%runoff  = ave%runoff   + d(k)%runoff   / div
       ave%runoffSv  = ave%runoffSv  + d(k)%runoffSv   / div
       ave%runoffSv_ice  = ave%runoffSv_ice   + d(k)%runoffSv_ice   / div
       ave%calving = ave%calving  + d(k)%calving  / div
       ave%bmelt   = ave%bmelt    + d(k)%bmelt    / div
       ave%fw_hosing = ave%fw_hosing + d(k)%fw_hosing / div
       ave%fw_flux_adj = ave%fw_flux_adj + d(k)%fw_flux_adj / div
       ave%fw_noise = ave%fw_noise + d(k)%fw_noise / div
       ave%flx_noise = ave%flx_noise + d(k)%flx_noise / div
       ave%buoy    = ave%buoy     + d(k)%buoy     / div
       ave%buoyS   = ave%buoyS    + d(k)%buoyS    / div
       ave%buoyT   = ave%buoyT    + d(k)%buoyT    / div
       ave%hft     = ave%hft      + d(k)%hft      / div
       ave%hfp     = ave%hfp      + d(k)%hfp      / div
       ave%hfa     = ave%hfa      + d(k)%hfa      / div
       ave%hftz    = ave%hftz     + d(k)%hftz     / div
       ave%hfpz    = ave%hfpz     + d(k)%hfpz     / div
       ave%hfaz    = ave%hfaz     + d(k)%hfaz     / div
       ave%fwt     = ave%fwt      + d(k)%fwt      / div
       ave%fwp     = ave%fwp      + d(k)%fwp      / div
       ave%fwa     = ave%fwa      + d(k)%fwa      / div
       ave%fwtz    = ave%fwtz     + d(k)%fwtz     / div
       ave%fwpz    = ave%fwpz     + d(k)%fwpz     / div
       ave%fwaz    = ave%fwaz     + d(k)%fwaz     / div
       ave%fayti   = ave%fayti    + d(k)%fayti    / div
       ave%fdyti   = ave%fdyti    + d(k)%fdyti    / div
       ave%nconv   = ave%nconv    + d(k)%nconv    / div
       ave%kven    = ave%kven     + d(k)%kven     / div
       ave%dconv   = ave%dconv    + d(k)%dconv    / div
       ave%dven    = ave%dven     + d(k)%dven     / div
       ave%mld     = ave%mld      + d(k)%mld      / div
       ave%mldst   = ave%mldst    + d(k)%mldst    / div
       ave%ke_tau= ave%ke_tau + d(k)%ke_tau / div
       ave%ssh = ave%ssh  + d(k)%ssh  / div
    end do

   return

  end subroutine ocn_ave


end module ocn_out
