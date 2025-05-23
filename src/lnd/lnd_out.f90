!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : l n d _ o u t
!
!  Purpose : land model diagnostic and output
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit and Andrey Ganopolski
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
module lnd_out
  
  use precision, only : wp, sp
  use dim_name, only: dim_time, dim_lon, dim_lat, dim_depth, dim_depth1, &
  & dim_day, dim_month, dim_nsurf, dim_ncarb, dim_nlit, dim_npft, dim_nsoil
  use timer, only : n_accel, doy, mon, year, year_clim, nyears, sec_day, sec_year, day_mon, &
  nday_year, nmon_year, year_ini, year_now
  use timer, only : time_out_lnd, time_soy_lnd, time_eoy_lnd, time_eom_lnd, nstep_mon_lnd, nyout_lnd, ny_out_ts, y_out_ts_clim, time_out_ts_clim
  use control, only : out_dir
  use constants, only : T0, c13_c12_std, c14_c_std
  use constants, only : e_sat_w, q_to_e, frac_vu
  use climber_grid, only : lat, lon, area
  use lnd_grid, only : nx, ny, npft, nsurf, nsoil, ncarb, nl, nl_l, nlc, dz, dz_c, z, z_c, rdz, z_int, z_l
  use lnd_grid, only : flag_tree, flag_grass, flag_shrub, flag_veg, flag_pft
  use lnd_grid, only : ic_min, ic_peat, ic_shelf, ic_ice, ic_lake
  use lnd_grid, only : is_veg, is_ice, is_lake
  use lnd_grid, only : i_bare, i_pft, i_surf, i_soil
  use lnd_params, only : dt, dt_day, peat_par, soilc_par
  use lnd_params, only : i_weathering, weath_gemco2_par, weath_uhh_par
  use lnd_params, only : write_surf, write_surf_n, write_carbon, write_soil, write_soil_par, write_lake, write_cons, l_daily_output
  use lnd_def, only : lnd_2d_class, lnd_0d_class

  use ncio
  !$use omp_lib

  implicit none

  private
  public :: lnd_diag_init, lnd_diag

  real(sp), parameter :: missing_value = -9999.

  integer :: nout

  integer :: nlit

  integer, dimension(:,:), allocatable :: pi_perm_mask

  type ts_out
    real(wp) :: co2, land, veg, ice, snow, shelf, lake, perm
    real(wp) :: landc, Cflx_atm_lnd, C13flx_atm_lnd, C14flx_atm_lnd, d13Cflx_atm_lnd, D14Cflx_atm_lnd
    real(wp) :: Cflx_burial, d13Cflx_burial
    real(wp) :: gpp, npp, npp_pimask, sresp, npp14, sresp14
    real(wp) :: runoff, runsur, calving, drain, evp, esur, trans, prc, t, wet, wettrop, wetextrop, inund, peat, peatpot
    real(wp) :: forest, grass, shrub, desert, crop, pasture, vegc, vegc_pimask
    real(wp) :: soilc, soilc60N, litterc, fastc, slowc, minc, peatc, icec, shelfc, lakec, permc, soilc_pimask, inertc
    real(wp) :: soilc1m, soilc60N1m, minc1m, peatc1m, icec1m, shelfc1m, lakec1m, permc1m, soilc1m_pimask
    real(wp) :: soilc3m, minc3m, peatc3m, icec3m, shelfc3m, lakec3m, permc3m, soilc3m_pimask
    real(wp) :: ch4, ch4trop, ch4shelf, ch4lake,  ch4extrop, d13ch4
    real(wp) :: dust_e
    real(wp) :: weath_carb, weath_sil, weath_loess
    real(wp) :: doc_export, poc_export
    real(wp), dimension(npft) :: pfts
  end type

  type surf_out
    real(wp), allocatable, dimension(:,:) :: f_land
    real(wp), allocatable, dimension(:,:) :: f_veg
    real(wp), allocatable, dimension(:,:) :: gdd5
    real(wp), allocatable, dimension(:,:) :: t2m_min_mon
    real(wp), allocatable, dimension(:,:) :: tlake_sur
    real(wp), allocatable, dimension(:,:,:) :: frac_surf
    real(wp), allocatable, dimension(:,:,:) :: swe
    real(wp), allocatable, dimension(:,:,:) :: swe_max
    integer, allocatable, dimension(:,:,:) :: msnow
    real(wp), allocatable, dimension(:,:,:) :: hsnow
    real(wp), allocatable, dimension(:,:,:) :: tsnow
    real(wp), allocatable, dimension(:,:,:) :: albsnw
    real(wp), allocatable, dimension(:,:) :: cti_lim
    real(wp), allocatable, dimension(:,:) :: fwet
    real(wp), allocatable, dimension(:,:) :: fwetmax
    real(wp), allocatable, dimension(:,:) :: fwetmaxpot
    real(wp), allocatable, dimension(:,:) :: wtab
    real(wp), allocatable, dimension(:,:) :: inf
    real(wp), allocatable, dimension(:,:) :: alt
    real(wp), allocatable, dimension(:,:) :: vpd
    real(wp), allocatable, dimension(:,:,:) :: alb
    real(wp), allocatable, dimension(:,:,:) :: alb_dif
    real(wp), allocatable, dimension(:,:,:) :: alb_dir
    real(wp), allocatable, dimension(:,:,:) :: wind
    real(wp), allocatable, dimension(:,:,:) :: ra
    real(wp), allocatable, dimension(:,:,:) :: rag
    real(wp), allocatable, dimension(:,:,:) :: Ri
    real(wp), allocatable, dimension(:,:,:) :: z0m
    real(wp), allocatable, dimension(:,:,:) :: z0h
    real(wp), allocatable, dimension(:,:,:) :: Ch
    real(wp), allocatable, dimension(:,:,:) :: rs
    real(wp), allocatable, dimension(:,:,:) :: betas
    real(wp), allocatable, dimension(:,:,:) :: le
    real(wp), allocatable, dimension(:,:,:) :: etot
    real(wp), allocatable, dimension(:,:,:) :: ecan
    real(wp), allocatable, dimension(:,:,:) :: esur
    real(wp), allocatable, dimension(:,:,:) :: trans
    real(wp), allocatable, dimension(:,:,:) :: sh
    real(wp), allocatable, dimension(:,:,:) :: slh
    real(wp), allocatable, dimension(:,:,:) :: ef
    real(wp), allocatable, dimension(:,:,:) :: lwnet
    real(wp), allocatable, dimension(:,:,:) :: swnet
    real(wp), allocatable, dimension(:,:,:) :: g
    real(wp), allocatable, dimension(:,:,:) :: flx_melt
    real(wp), allocatable, dimension(:,:,:) :: tskin
    real(wp), allocatable, dimension(:,:,:) :: tskin_amp
    real(wp), allocatable, dimension(:,:,:) :: scf
    real(wp), allocatable, dimension(:,:,:) :: fcanwat
    real(wp), allocatable, dimension(:,:,:) :: fcansn
    real(wp), allocatable, dimension(:,:,:) :: runoff
    real(wp), allocatable, dimension(:,:,:) :: runsur
    real(wp), allocatable, dimension(:,:,:) :: calving
    real(wp), allocatable, dimension(:,:,:) :: drain
    real(wp), allocatable, dimension(:,:,:) :: rsursub
    real(wp), allocatable, dimension(:,:,:) :: snowmelt
    real(wp), allocatable, dimension(:,:,:) :: icemelt
    real(wp), allocatable, dimension(:,:) :: dust_e_d
    real(wp), allocatable, dimension(:,:) :: dust_e_g
    real(wp), allocatable, dimension(:,:) :: dust_e_s
    real(wp), allocatable, dimension(:,:) :: dust_e
    real(wp), allocatable, dimension(:,:) :: dust_dep
    real(wp), allocatable, dimension(:,:,:) :: snow_grain
    real(wp), allocatable, dimension(:,:,:) :: dust_con
  end type
  type surf_g_out
    real(wp), dimension(nx,ny) :: le, etot, ecan, esur, trans, sh, slh, ef, swnet, lwnet, g, flx_melt, tskin, tskin_amp, alb, alb_dif, alb_dir, ra, rag, Ri, Ch, rs
  end type
 
  type carbon_out
    real(wp), dimension(nx,ny) :: fpeat, fpeatpot, dCpeatdt, acroh, catoh
    real(wp), dimension(nx,ny,npft) :: lai, sai, xi, wue, gcan, gpp, npp, npp13, npp14, aresp, veg_h, pfts, seeds, vegc, lambda
    real(wp), dimension(nx,ny,npft) :: gamma_dist, gamma_luc, gamma_ice
    real(wp), dimension(nx,ny) :: bare
    real(wp), dimension(nx,ny,ncarb) :: sresp, sresp13, sresp14, soilc, litter
    real(wp), allocatable, dimension(:,:,:,:) :: litter_prof, litter13_prof, litter14_prof
    real(wp), allocatable, dimension(:,:,:,:) :: litterc_prof, fastc_prof, slowc_prof, soilc_prof
    real(wp), allocatable, dimension(:,:,:,:) :: litterc13_prof, fastc13_prof, slowc13_prof, soilc13_prof
    real(wp), allocatable, dimension(:,:,:,:) :: litterc14_prof, fastc14_prof, slowc14_prof, soilc14_prof
    real(wp), dimension(nx,ny) :: ch4, ch4wet, ch4shelf, ch4lake, ch4peat
    real(wp), dimension(nx,ny,npft) :: disc, d13c_veg, D14c_veg
    real(wp), dimension(nx,ny,ncarb) :: d13c_soil, D14c_soil
    real(wp), allocatable, dimension(:,:,:,:) :: d13c_litter_prof, D14c_litter_prof
    real(wp), allocatable, dimension(:,:,:,:) :: d13c_prof, D14c_prof
    real(wp), dimension(nx,ny) :: doc_export, poc_export 
    real(wp), dimension(nx,ny) :: weath_carb, weath_sil
    real(wp), allocatable, dimension(:,:,:) :: lithology
  end type
  type carbon_g_out
    real(wp), dimension(nx,ny) :: lai, xi, wue, gcan, gpp, npp, npp13, npp14, aresp, sresp, sresp13, sresp14, veg_h, litter, vegc, soilc
    real(wp), dimension(nx,ny) :: Cflx_atm_lnd, C13flx_atm_lnd, C14flx_atm_lnd, d13Cflx_atm_lnd, D14Cflx_atm_lnd 
    real(wp), dimension(nx,ny) :: disc, d13c_veg, d13c_soil, D14c_veg, D14c_soil
  end type

  type soil_out
    real(wp), allocatable, dimension(:,:,:) :: tsoil, tice, tsublake, thetaw, thetai, theta, thetas, fthetas
  end type

  type soil_par_out
    real(wp), allocatable, dimension(:,:,:) :: lambda_i, cap_i, lambda_if, cap_if, kappa
    real(wp), allocatable, dimension(:,:,:) :: ksat, psisat, bi
    real(wp), allocatable, dimension(:,:,:) :: ftemp, fmoist, fdepth
    real(wp), allocatable, dimension(:,:,:) :: fsoc, theta_sat, theta_field, theta_wilt, psi
  end type

  type lake_out
    real(wp), allocatable, dimension(:,:,:) :: t_lake, f_i_lake, lambda_lake
    real(wp), allocatable, dimension(:,:) :: lakebal
    real(wp), allocatable, dimension(:,:) :: f_lake_ice
    real(wp), allocatable, dimension(:,:) :: h_lake
    real(wp), allocatable, dimension(:,:) :: h_lake_conv
    real(wp), allocatable, dimension(:,:) :: h_lake_mix
  end type

  type cons_out
    real(wp), dimension(nx,ny,nsurf) :: econs_su1, econs_su2
    real(wp), dimension(nx,ny) :: econs_so
    real(wp), dimension(nx,ny,nsoil) :: wcons
    real(wp), dimension(nx,ny,ncarb) :: ccons_s, ccons13_s, ccons14_s
    real(wp), dimension(nx,ny) :: ccons_v, ccons13_v, ccons14_v
  end type


  type(ts_out), allocatable :: mon_ts(:,:), ann_ts(:)
  type(surf_out) :: day_su(nday_year), mon_su(nmon_year), ann_su
  type(surf_g_out) :: mon_su_g(nmon_year), ann_su_g
  type(carbon_out) :: mon_c(nmon_year), ann_c
  type(carbon_g_out) :: mon_c_g(nmon_year), ann_c_g
  type(soil_out) :: mon_s(nmon_year), ann_s
  type(soil_par_out) :: mon_sp(nmon_year), ann_sp
  type(lake_out) :: mon_l(nmon_year), ann_l
  type(cons_out) :: mon_b(nmon_year), ann_b



contains 
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  l a n d _ d i a g _ i n i t
  ! Purpose  :  Initialize netcdf output for paladyn
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lnd_diag_init

    implicit none

    integer :: k


    nout = 0

    if (i_weathering.eq.1) then
      nlit = weath_gemco2_par%nlit
    else if (i_weathering.eq.2) then
      nlit = weath_uhh_par%nlit
    endif

    ! initialize output files
    call ts_nc(trim(out_dir)//"/lnd_ts.nc")
    if (write_surf) call surf_nc(trim(out_dir)//"/lnd_surf.nc")
    if (write_surf .and. l_daily_output) call surf_daily_nc(trim(out_dir)//"/lnd_surf_daily.nc")
    if (write_soil) call soil_nc(trim(out_dir)//"/lnd_soil.nc")
    if (write_carbon) call carbon_nc(trim(out_dir)//"/lnd_carb.nc")
    if (write_soil_par) call soil_par_nc(trim(out_dir)//"/lnd_soil_par.nc")
    if (write_lake) call lake_nc(trim(out_dir)//"/lnd_lake.nc")
    if (write_cons) call cons_nc(trim(out_dir)//"/lnd_conservation.nc")

    ! allocate
    allocate(mon_ts(ny_out_ts,nmon_year))
    allocate(ann_ts(ny_out_ts))

    do k=1,nmon_year
      allocate(mon_su(k)%f_land(nx,ny))
      allocate(mon_su(k)%f_veg(nx,ny))
      allocate(mon_su(k)%tlake_sur(nx,ny))
      allocate(mon_su(k)%frac_surf(nx,ny,nsurf))
      allocate(mon_su(k)%swe(nx,ny,nsoil))
      allocate(mon_su(k)%swe_max(nx,ny,nsoil))
      allocate(mon_su(k)%hsnow(nx,ny,nsoil))
      allocate(mon_su(k)%tsnow(nx,ny,nsoil))
      allocate(mon_su(k)%albsnw(nx,ny,nsoil))
      allocate(mon_su(k)%cti_lim(nx,ny))
      allocate(mon_su(k)%fwet(nx,ny))
      allocate(mon_su(k)%fwetmax(nx,ny))
      allocate(mon_su(k)%wtab(nx,ny))
      allocate(mon_su(k)%inf(nx,ny))
      allocate(mon_su(k)%alt(nx,ny))
      allocate(mon_su(k)%vpd(nx,ny))
      allocate(mon_su(k)%alb(nx,ny,nsurf))
      allocate(mon_su(k)%alb_dif(nx,ny,nsurf))
      allocate(mon_su(k)%alb_dir(nx,ny,nsurf))
      allocate(mon_su(k)%wind(nx,ny,nsurf))
      allocate(mon_su(k)%ra(nx,ny,nsurf))
      allocate(mon_su(k)%rag(nx,ny,npft))
      allocate(mon_su(k)%Ri(nx,ny,nsurf))
      allocate(mon_su(k)%z0m(nx,ny,nsurf))
      allocate(mon_su(k)%z0h(nx,ny,nsurf))
      allocate(mon_su(k)%Ch(nx,ny,nsurf))
      allocate(mon_su(k)%rs(nx,ny,nsurf))
      allocate(mon_su(k)%betas(nx,ny,nsurf))
      allocate(mon_su(k)%le(nx,ny,nsurf))
      allocate(mon_su(k)%etot(nx,ny,nsurf))
      allocate(mon_su(k)%ecan(nx,ny,nsurf))
      allocate(mon_su(k)%esur(nx,ny,nsurf))
      allocate(mon_su(k)%trans(nx,ny,nsurf))
      allocate(mon_su(k)%sh(nx,ny,nsurf))
      allocate(mon_su(k)%slh(nx,ny,nsurf))
      allocate(mon_su(k)%ef(nx,ny,nsurf))
      allocate(mon_su(k)%lwnet(nx,ny,nsurf))
      allocate(mon_su(k)%swnet(nx,ny,nsurf))
      allocate(mon_su(k)%g(nx,ny,nsurf))
      allocate(mon_su(k)%flx_melt(nx,ny,nsurf))
      allocate(mon_su(k)%tskin(nx,ny,nsurf))
      allocate(mon_su(k)%tskin_amp(nx,ny,nsurf))
      allocate(mon_su(k)%scf(nx,ny,nsurf))
      allocate(mon_su(k)%fcanwat(nx,ny,nsurf))
      allocate(mon_su(k)%fcansn(nx,ny,nsurf))
      allocate(mon_su(k)%runoff(nx,ny,nsoil))
      allocate(mon_su(k)%runsur(nx,ny,nsoil))
      allocate(mon_su(k)%calving(nx,ny,nsoil))
      allocate(mon_su(k)%drain(nx,ny,nsoil))
      allocate(mon_su(k)%rsursub(nx,ny,nsoil))
      allocate(mon_su(k)%snowmelt(nx,ny,nsoil))
      allocate(mon_su(k)%icemelt(nx,ny,nsoil))
      allocate(mon_su(k)%dust_e_d(nx,ny))
      allocate(mon_su(k)%dust_e_g(nx,ny))
      allocate(mon_su(k)%dust_e_s(nx,ny))
      allocate(mon_su(k)%dust_e(nx,ny))
      allocate(mon_su(k)%dust_dep(nx,ny))
      allocate(mon_su(k)%snow_grain(nx,ny,nsoil))
      allocate(mon_su(k)%dust_con(nx,ny,nsoil))

      allocate(mon_s(k)%tsoil(nx,ny,nl))
      allocate(mon_s(k)%tice(nx,ny,nl))
      allocate(mon_s(k)%tsublake(nx,ny,nl))
      allocate(mon_s(k)%thetaw(nx,ny,nl))
      allocate(mon_s(k)%thetai(nx,ny,nl))
      allocate(mon_s(k)%theta(nx,ny,nl))
      allocate(mon_s(k)%thetas(nx,ny,nl))
      allocate(mon_s(k)%fthetas(nx,ny,nl))

      allocate(mon_sp(k)%lambda_if(nx,ny,0:nl))
      allocate(mon_sp(k)%cap_if(nx,ny,0:nl))
      allocate(mon_sp(k)%lambda_i(nx,ny,0:nl))
      allocate(mon_sp(k)%cap_i(nx,ny,0:nl))
      allocate(mon_sp(k)%kappa(nx,ny,nl))
      allocate(mon_sp(k)%ksat(nx,ny,nl))
      allocate(mon_sp(k)%psisat(nx,ny,nl))
      allocate(mon_sp(k)%bi(nx,ny,nl))
      allocate(mon_sp(k)%theta_sat(nx,ny,nl))
      allocate(mon_sp(k)%theta_field(nx,ny,nl))
      allocate(mon_sp(k)%theta_wilt(nx,ny,nl))
      allocate(mon_sp(k)%fsoc(nx,ny,nl))
      allocate(mon_sp(k)%psi(nx,ny,nl))
      allocate(mon_sp(k)%ftemp(nx,ny,nl))
      allocate(mon_sp(k)%fmoist(nx,ny,nl))
      allocate(mon_sp(k)%fdepth(nx,ny,nl))

      allocate(mon_l(k)%t_lake(nx,ny,nl_l))
      allocate(mon_l(k)%f_i_lake(nx,ny,nl_l))
      allocate(mon_l(k)%f_lake_ice(nx,ny))
      allocate(mon_l(k)%lambda_lake(nx,ny,nl_l))
      allocate(mon_l(k)%lakebal(nx,ny))
      allocate(mon_l(k)%h_lake_conv(nx,ny))
      allocate(mon_l(k)%h_lake_mix(nx,ny))

      allocate(mon_c(k)%litter_prof(nx,ny,nlc,ncarb))
      allocate(mon_c(k)%litter13_prof(nx,ny,nlc,ncarb))
      allocate(mon_c(k)%litter14_prof(nx,ny,nlc,ncarb))
    enddo

    do k=1,nday_year
      allocate(day_su(k)%frac_surf(nx,ny,nsurf))
      allocate(day_su(k)%le(nx,ny,nsurf))
      allocate(day_su(k)%sh(nx,ny,nsurf))
      allocate(day_su(k)%lwnet(nx,ny,nsurf))
      allocate(day_su(k)%swnet(nx,ny,nsurf))
      allocate(day_su(k)%g(nx,ny,nsurf))
      allocate(day_su(k)%flx_melt(nx,ny,nsurf))
      allocate(day_su(k)%tskin(nx,ny,nsurf))
      allocate(day_su(k)%ra(nx,ny,nsurf))
      allocate(day_su(k)%rag(nx,ny,npft))
      allocate(day_su(k)%Ri(nx,ny,nsurf))
      allocate(day_su(k)%rs(nx,ny,nsurf))
      allocate(day_su(k)%betas(nx,ny,nsurf))
      allocate(day_su(k)%ecan(nx,ny,nsurf))
      allocate(day_su(k)%msnow(nx,ny,nsoil))
      allocate(day_su(k)%hsnow(nx,ny,nsoil))
      allocate(day_su(k)%swe(nx,ny,nsoil))
      allocate(day_su(k)%swe_max(nx,ny,nsoil))
      allocate(day_su(k)%tsnow(nx,ny,nsoil))
      allocate(day_su(k)%albsnw(nx,ny,nsoil))
    enddo

    allocate(ann_su%f_land(nx,ny))
    allocate(ann_su%f_veg(nx,ny))
    allocate(ann_su%gdd5(nx,ny))
    allocate(ann_su%t2m_min_mon(nx,ny))
    allocate(ann_su%tlake_sur(nx,ny))
    allocate(ann_su%frac_surf(nx,ny,nsurf))
    allocate(ann_su%swe(nx,ny,nsoil))
    allocate(ann_su%swe_max(nx,ny,nsoil))
    allocate(ann_su%hsnow(nx,ny,nsoil))
    allocate(ann_su%tsnow(nx,ny,nsoil))
    allocate(ann_su%albsnw(nx,ny,nsoil))
    allocate(ann_su%cti_lim(nx,ny))
    allocate(ann_su%fwet(nx,ny))
    allocate(ann_su%fwetmax(nx,ny))
    allocate(ann_su%fwetmaxpot(nx,ny))
    allocate(ann_su%wtab(nx,ny))
    allocate(ann_su%inf(nx,ny))
    allocate(ann_su%alt(nx,ny))
    allocate(ann_su%vpd(nx,ny))
    allocate(ann_su%alb(nx,ny,nsurf))
    allocate(ann_su%alb_dif(nx,ny,nsurf))
    allocate(ann_su%alb_dir(nx,ny,nsurf))
    allocate(ann_su%wind(nx,ny,nsurf))
    allocate(ann_su%ra(nx,ny,nsurf))
    allocate(ann_su%rag(nx,ny,npft))
    allocate(ann_su%Ri(nx,ny,nsurf))
    allocate(ann_su%z0m(nx,ny,nsurf))
    allocate(ann_su%z0h(nx,ny,nsurf))
    allocate(ann_su%Ch(nx,ny,nsurf))
    allocate(ann_su%rs(nx,ny,nsurf))
    allocate(ann_su%betas(nx,ny,nsurf))
    allocate(ann_su%le(nx,ny,nsurf))
    allocate(ann_su%etot(nx,ny,nsurf))
    allocate(ann_su%ecan(nx,ny,nsurf))
    allocate(ann_su%esur(nx,ny,nsurf))
    allocate(ann_su%trans(nx,ny,nsurf))
    allocate(ann_su%sh(nx,ny,nsurf))
    allocate(ann_su%slh(nx,ny,nsurf))
    allocate(ann_su%ef(nx,ny,nsurf))
    allocate(ann_su%lwnet(nx,ny,nsurf))
    allocate(ann_su%swnet(nx,ny,nsurf))
    allocate(ann_su%g(nx,ny,nsurf))
    allocate(ann_su%flx_melt(nx,ny,nsurf))
    allocate(ann_su%tskin(nx,ny,nsurf))
    allocate(ann_su%tskin_amp(nx,ny,nsurf))
    allocate(ann_su%scf(nx,ny,nsurf))
    allocate(ann_su%fcanwat(nx,ny,nsurf))
    allocate(ann_su%fcansn(nx,ny,nsurf))
    allocate(ann_su%runoff(nx,ny,nsoil))
    allocate(ann_su%runsur(nx,ny,nsoil))
    allocate(ann_su%calving(nx,ny,nsoil))
    allocate(ann_su%drain(nx,ny,nsoil))
    allocate(ann_su%rsursub(nx,ny,nsoil))
    allocate(ann_su%snowmelt(nx,ny,nsoil))
    allocate(ann_su%icemelt(nx,ny,nsoil))
    allocate(ann_su%dust_e_d(nx,ny))
    allocate(ann_su%dust_e_g(nx,ny))
    allocate(ann_su%dust_e_s(nx,ny))
    allocate(ann_su%dust_e(nx,ny))
    allocate(ann_su%dust_dep(nx,ny))
    allocate(ann_su%snow_grain(nx,ny,nsoil))
    allocate(ann_su%dust_con(nx,ny,nsoil))

    allocate(ann_s%tsoil(nx,ny,nl))
    allocate(ann_s%tice(nx,ny,nl))
    allocate(ann_s%tsublake(nx,ny,nl))
    allocate(ann_s%thetaw(nx,ny,nl))
    allocate(ann_s%thetai(nx,ny,nl))
    allocate(ann_s%theta(nx,ny,nl))
    allocate(ann_s%thetas(nx,ny,nl))
    allocate(ann_s%fthetas(nx,ny,nl))

    allocate(ann_c%litter_prof(nx,ny,nlc,ncarb))
    allocate(ann_c%litter13_prof(nx,ny,nlc,ncarb))
    allocate(ann_c%litter14_prof(nx,ny,nlc,ncarb))
    allocate(ann_c%litterc_prof(nx,ny,nlc,ncarb))
    allocate(ann_c%fastc_prof(nx,ny,nlc,ncarb))
    allocate(ann_c%slowc_prof(nx,ny,nlc,ncarb))
    allocate(ann_c%soilc_prof(nx,ny,nlc,ncarb))
    allocate(ann_c%litterc13_prof(nx,ny,nlc,ncarb))
    allocate(ann_c%fastc13_prof(nx,ny,nlc,ncarb))
    allocate(ann_c%slowc13_prof(nx,ny,nlc,ncarb))
    allocate(ann_c%soilc13_prof(nx,ny,nlc,ncarb))
    allocate(ann_c%litterc14_prof(nx,ny,nlc,ncarb))
    allocate(ann_c%fastc14_prof(nx,ny,nlc,ncarb))
    allocate(ann_c%slowc14_prof(nx,ny,nlc,ncarb))
    allocate(ann_c%soilc14_prof(nx,ny,nlc,ncarb))
    allocate(ann_c%d13c_litter_prof(nx,ny,nlc,ncarb))
    allocate(ann_c%D14c_litter_prof(nx,ny,nlc,ncarb))
    allocate(ann_c%d13c_prof(nx,ny,nlc,ncarb))
    allocate(ann_c%D14c_prof(nx,ny,nlc,ncarb))
    allocate(ann_c%lithology(nx,ny,nlit))

    allocate(ann_sp%lambda_if(nx,ny,0:nl))
    allocate(ann_sp%cap_if(nx,ny,0:nl))
    allocate(ann_sp%lambda_i(nx,ny,0:nl))
    allocate(ann_sp%cap_i(nx,ny,0:nl))
    allocate(ann_sp%kappa(nx,ny,nl))
    allocate(ann_sp%ksat(nx,ny,nl))
    allocate(ann_sp%psisat(nx,ny,nl))
    allocate(ann_sp%bi(nx,ny,nl))
    allocate(ann_sp%theta_sat(nx,ny,nl))
    allocate(ann_sp%theta_field(nx,ny,nl))
    allocate(ann_sp%theta_wilt(nx,ny,nl))
    allocate(ann_sp%fsoc(nx,ny,nl))
    allocate(ann_sp%psi(nx,ny,nl))
    allocate(ann_sp%ftemp(nx,ny,nl))
    allocate(ann_sp%fmoist(nx,ny,nl))
    allocate(ann_sp%fdepth(nx,ny,nl))

    allocate(ann_l%t_lake(nx,ny,nl_l))
    allocate(ann_l%lambda_lake(nx,ny,nl_l))
    allocate(ann_l%f_i_lake(nx,ny,nl_l))
    allocate(ann_l%f_lake_ice(nx,ny))
    allocate(ann_l%lakebal(nx,ny))
    allocate(ann_l%h_lake(nx,ny))
    allocate(ann_l%h_lake_conv(nx,ny))
    allocate(ann_l%h_lake_mix(nx,ny))

    ! read PI permafrost area mask
    allocate(pi_perm_mask(nx,ny))
    call nc_read('input/pi_perm_mask.nc', "pi_perm_mask",pi_perm_mask)


  end subroutine lnd_diag_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  l n d _ d i a g
  ! Purpose  :  netcdf output for paladyn
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lnd_diag(lnd,lndp)

    implicit none

    type(lnd_2d_class), intent(in) :: lnd(:,:)
    type(lnd_0d_class), intent(in) :: lndp

    integer :: k, i, j, m, n, y
    real(wp) :: mon_avg
    real(wp), dimension(npft) :: frac

    ! global values
    real(wp) :: global_gpp
    real(wp) :: global_npp, global_npp_pimask, global_npp13, global_npp14, global_sresp, global_sresp13, global_sresp14
    real(wp) :: global_t, global_snow
    real(wp) :: global_run, global_runsur, global_calving, global_drain, global_evp, global_esur, global_trans, global_prc
    real(wp) :: global_wet, global_wettrop, global_wetextrop, global_inund, global_peat, global_peatpot, tmp
    real(wp) :: global_forest, global_grass, global_shrub, global_desert, global_crop, global_pasture
    real(wp), dimension(npft) :: global_pfts
    real(wp) :: global_vegc, global_vegc_pimask
    real(wp) :: litterc, fastc, slowc
    real(wp) :: minc, peatc, icec, shelfc, lakec
    real(wp) :: global_litterc, global_fastc, global_slowc
    real(wp) :: global_soilc60N, global_soilc60N1m
    real(wp) :: global_minc, global_peatc, global_permc, global_soilc_pimask, global_icec, global_shelfc, global_lakec, global_inertc
    real(wp) :: global_minc1m, global_peatc1m, global_permc1m, global_soilc1m_pimask, global_icec1m, global_shelfc1m, global_lakec1m
    real(wp) :: global_minc3m, global_peatc3m, global_permc3m, global_soilc3m_pimask, global_icec3m, global_shelfc3m, global_lakec3m
    real(wp) :: global_ch4, global_ch4trop, global_ch4shelf, global_ch4lake, global_ch4extrop, global_c13h4
    real(wp) :: global_dust_e
    real(wp) :: npp_ij, npp13_ij, npp14_ij, sresp_ij, sresp13_ij, sresp14_ij
    real(wp) :: sum_frac


    ! current index
    y = y_out_ts_clim

    mon_avg = 1._wp/nstep_mon_lnd

    ! initialize globally integrated values
    global_gpp = 0._wp
    global_run = 0._wp
    global_runsur = 0._wp
    global_calving= 0._wp
    global_drain = 0._wp
    global_evp   = 0._wp
    global_esur  = 0._wp
    global_trans = 0._wp
    global_t   = 0._wp
    global_prc   = 0._wp
    global_snow = 0._wp
    global_dust_e = 0._wp

    !$omp parallel do collapse(2) private(i,j,k) &
    !$omp reduction(+:global_gpp,global_evp,global_esur,global_trans,global_t,global_prc,global_runsur,global_calving,global_drain,global_run,global_snow,global_dust_e)
    do j = 1, ny
      do i = 1, nx
        if( lnd(i,j)%f_land.gt.0._wp ) then
          ! GPP
          do k = 1, npft
            if (lnd(i,j)%frac_surf(k).gt.0._wp) then
              global_gpp = global_gpp + lnd(i,j)%gpp(k) * lnd(i,j)%frac_surf(k) * area(i,j) * 1.d-12 * dt_day ! PgC
            endif
          enddo
          do k = 1, nsurf
            if (lnd(i,j)%frac_surf(k).gt.0._wp) then
              ! evapotranspiration
              global_evp  = global_evp   + lnd(i,j)%et(k) * lnd(i,j)%frac_surf(k) * area(i,j) * 1.d-15 * dt ! 10^15 kg/m2/sum_dt
              ! surface evaporation
              global_esur  = global_esur   + lnd(i,j)%evap_surface(k) * lnd(i,j)%frac_surf(k) * area(i,j) * 1.d-15 * dt ! 10^15 kg/m2/sum_dt
              ! evapotranspiration
              global_trans = global_trans  + lnd(i,j)%transpiration(k) * lnd(i,j)%frac_surf(k) * area(i,j) * 1.d-15 * dt ! 10^15 kg/m2/sum_dt
            endif
            ! land temperature
            if(flag_veg(k).eq.1 .and. lnd(i,j)%frac_surf(k).gt.0._wp) global_t  = global_t   + lnd(i,j)%t2m(k) * lnd(i,j)%frac_surf(k) * area(i,j) 
          enddo
          ! precipitation
          if (lnd(i,j)%f_shelf.lt.1._wp) then
            global_prc = global_prc + sum((lnd(i,j)%rain(:)+lnd(i,j)%snow(:))*lnd(i,j)%frac_surf(:)) * area(i,j) * 1.d-15 * dt  ! 10^15 kg/m2/sum_dt
          endif
          ! surface runoff
          global_runsur = global_runsur &
          + (lnd(i,j)%f_veg*lnd(i,j)%runoff_sur(is_veg) + lnd(i,j)%f_ice*lnd(i,j)%runoff_sur(is_ice) &
          + lnd(i,j)%f_lake*lnd(i,j)%runoff_sur(is_lake)) &
          * area(i,j) * 1.d-15 * dt! 10^15 kg/m2/sum_dt
          ! 'calving' 
          global_calving = global_calving &
          + (lnd(i,j)%f_veg*lnd(i,j)%calving(is_veg) + lnd(i,j)%f_ice*lnd(i,j)%calving(is_ice) &
          + lnd(i,j)%f_lake*lnd(i,j)%calving(is_lake)) &
          * area(i,j) * 1.d-15 * dt! 10^15 kg/m2/sum_dt
          ! drainage
          global_drain = global_drain &
          + (lnd(i,j)%f_veg*lnd(i,j)%drainage(is_veg) + lnd(i,j)%f_ice*lnd(i,j)%drainage(is_ice) &
          + lnd(i,j)%f_lake*lnd(i,j)%drainage(is_lake)) &
          * area(i,j) * 1.d-15 * dt! 10^15 kg/m2/sum_dt
          ! total runoff (surface runoff + drainage)
          global_run = global_run &
          + (lnd(i,j)%f_veg*lnd(i,j)%runoff(is_veg) + lnd(i,j)%f_ice*lnd(i,j)%runoff(is_ice) + lnd(i,j)%f_lake*lnd(i,j)%runoff(is_lake)) &
          * area(i,j) * 1.d-15 * dt! 10^15 kg/m2/sum_dt
          ! snow volume
          global_snow = global_snow &
            + (lnd(i,j)%w_snow(is_veg)/1000._wp * lnd(i,j)%f_veg &
            + lnd(i,j)%w_snow(is_lake)/1000._wp * lnd(i,j)%f_lake) &
            * area(i,j) * 1.d-12 ! 10^12 m3
          ! dust emissions
          global_dust_e = global_dust_e + lnd(i,j)%dust_emis * area(i,j) * 1.e-9 * dt ! kg/m2/s*s*m2*Tg/kg = Tg
        endif
      enddo
    enddo
    !$omp end parallel do

    ! reset at first call of year
    if( time_soy_lnd ) then
      do k=1,nmon_year
        mon_ts(y,k)%gpp = 0._wp
        mon_ts(y,k)%runoff = 0._wp
        mon_ts(y,k)%runsur = 0._wp
        mon_ts(y,k)%calving = 0._wp
        mon_ts(y,k)%drain = 0._wp
        mon_ts(y,k)%evp   = 0._wp
        mon_ts(y,k)%esur   = 0._wp
        mon_ts(y,k)%trans  = 0._wp
        mon_ts(y,k)%prc   = 0._wp
        mon_ts(y,k)%t   = 0._wp
        mon_ts(y,k)%snow= 0._wp
        mon_ts(y,k)%dust_e= 0._wp
      enddo
    endif

    ! cumulate/average from daily to monthly
    mon_ts(y,mon)%gpp    = mon_ts(y,mon)%gpp + global_gpp 
    mon_ts(y,mon)%runoff = mon_ts(y,mon)%runoff + global_run 
    mon_ts(y,mon)%runsur = mon_ts(y,mon)%runsur + global_runsur
    mon_ts(y,mon)%calving= mon_ts(y,mon)%calving+ global_calving
    mon_ts(y,mon)%drain  = mon_ts(y,mon)%drain + global_drain
    mon_ts(y,mon)%evp    = mon_ts(y,mon)%evp   + global_evp
    mon_ts(y,mon)%esur   = mon_ts(y,mon)%esur   + global_esur
    mon_ts(y,mon)%trans  = mon_ts(y,mon)%trans  + global_trans
    mon_ts(y,mon)%prc    = mon_ts(y,mon)%prc   + global_prc
    mon_ts(y,mon)%t      = mon_ts(y,mon)%t   + (global_t / sum(lnd%f_veg * area))  * mon_avg  ! K
    mon_ts(y,mon)%snow   = mon_ts(y,mon)%snow+ global_snow * mon_avg
    mon_ts(y,mon)%dust_e = mon_ts(y,mon)%dust_e+ global_dust_e


    ! at the end of month
    if( time_eom_lnd ) then

      ! initialize global values
      global_npp = 0._wp
      global_npp_pimask = 0._wp
      global_npp13 = 0._wp
      global_npp14 = 0._wp
      global_sresp  = 0._wp
      global_sresp13= 0._wp
      global_sresp14= 0._wp
      global_peat   = 0._wp
      global_peatpot   = 0._wp
      global_wet = 0._wp
      global_wettrop = 0._wp
      global_wetextrop = 0._wp
      global_inund  = 0._wp
      global_forest = 0._wp
      global_pfts   = 0._wp
      global_grass  = 0._wp
      global_shrub  = 0._wp
      global_desert = 0._wp
      global_crop   = 0._wp
      global_pasture= 0._wp
      global_vegc   = 0._wp
      global_vegc_pimask   = 0._wp
      global_litterc= 0._wp
      global_fastc  = 0._wp
      global_slowc  = 0._wp
      global_soilc60N   = 0._wp
      global_soilc60N1m   = 0._wp
      global_minc   = 0._wp
      global_peatc  = 0._wp
      global_permc  = 0._wp
      global_soilc_pimask  = 0._wp
      global_icec   = 0._wp
      global_shelfc = 0._wp
      global_lakec  = 0._wp
      global_inertc = 0._wp
      global_minc1m   = 0._wp
      global_peatc1m  = 0._wp
      global_permc1m  = 0._wp
      global_soilc1m_pimask  = 0._wp
      global_icec1m   = 0._wp
      global_shelfc1m = 0._wp
      global_lakec1m  = 0._wp
      global_minc3m   = 0._wp
      global_peatc3m  = 0._wp
      global_permc3m  = 0._wp
      global_soilc3m_pimask  = 0._wp
      global_icec3m   = 0._wp
      global_shelfc3m = 0._wp
      global_lakec3m  = 0._wp
      global_ch4    = 0._wp
      global_ch4trop= 0._wp
      global_ch4shelf = 0._wp
      global_ch4lake  = 0._wp
      global_ch4extrop= 0._wp
      global_c13h4  = 0._wp

      !$omp parallel do collapse(2) private(i,j,k,n,npp_ij,npp13_ij,npp14_ij,sresp_ij,sresp13_ij,sresp14_ij,tmp,litterc,fastc, slowc,minc,peatc,icec,shelfc,lakec) &
      !$omp reduction(+:global_npp,global_npp_pimask,global_npp13,global_npp14,global_sresp,global_sresp13,global_sresp14,global_inund) &
      !$omp reduction(+:global_wettrop,global_wet,global_wetextrop,global_peat,global_peatpot) &
      !$omp reduction(+:global_ch4,global_ch4trop,global_ch4extrop,global_ch4shelf,global_ch4lake) &
      !$omp reduction(+:global_c13h4,global_pfts,global_vegc,global_vegc_pimask,global_desert,global_forest,global_grass,global_shrub,global_crop,global_pasture) &
      !$omp reduction(+:global_litterc,global_fastc,global_slowc,global_soilc60N,global_soilc60N1m) &
      !$omp reduction(+:global_minc,global_peatc,global_icec,global_shelfc,global_lakec,global_permc,global_soilc_pimask,global_inertc) &
      !$omp reduction(+:global_minc1m,global_peatc1m,global_icec1m,global_shelfc1m,global_lakec1m,global_permc1m,global_soilc1m_pimask) &
      !$omp reduction(+:global_minc3m,global_peatc3m,global_icec3m,global_shelfc3m,global_lakec3m,global_permc3m,global_soilc3m_pimask) 
      do i = 1, nx
        do j = 1, ny
          if( lnd(i,j)%f_land.gt.0._wp ) then
            ! NPP
            npp_ij = lnd(i,j)%npp_real * lnd(i,j)%f_veg * sec_day*day_mon * area(i,j) * 1.d-12 ! PgC/mon
            npp13_ij = lnd(i,j)%npp13_real * lnd(i,j)%f_veg * sec_day*day_mon * area(i,j) * 1.d-12
            npp14_ij = lnd(i,j)%npp14_real * lnd(i,j)%f_veg * sec_day*day_mon * area(i,j) * 1.d-12
            global_npp = global_npp + npp_ij
            if( pi_perm_mask(i,j).eq.1 ) then            
              global_npp_pimask = global_npp_pimask + npp_ij
            endif
            global_npp13 = global_npp13 + npp13_ij
            global_npp14 = global_npp14 + npp14_ij 
            ! soil respiration
            sresp_ij = ( lnd(i,j)%soil_resp(ic_min)*(lnd(i,j)%f_veg-lnd(i,j)%f_peat) &
            + lnd(i,j)%soil_resp(ic_peat)*lnd(i,j)%f_peat &
            + lnd(i,j)%soil_resp(ic_shelf)*lnd(i,j)%f_shelf &
            + lnd(i,j)%soil_resp(ic_lake)*lnd(i,j)%f_lake &
            + lnd(i,j)%soil_resp(ic_ice)*lnd(i,j)%f_ice_grd) &
            * sec_day*day_mon * area(i,j) * 1.d-12 ! PgC/mon
            sresp13_ij = ( lnd(i,j)%soil_resp13(ic_min)*(lnd(i,j)%f_veg-lnd(i,j)%f_peat) &
            + lnd(i,j)%soil_resp13(ic_peat)*lnd(i,j)%f_peat &
            + lnd(i,j)%soil_resp13(ic_shelf)*lnd(i,j)%f_shelf &
            + lnd(i,j)%soil_resp13(ic_lake)*lnd(i,j)%f_lake &
            + lnd(i,j)%soil_resp13(ic_ice)*lnd(i,j)%f_ice_grd) &
            * sec_day*day_mon * area(i,j) * 1.d-12 ! PgC/mon
            sresp14_ij = ( lnd(i,j)%soil_resp14(ic_min)*(lnd(i,j)%f_veg-lnd(i,j)%f_peat) &
            + lnd(i,j)%soil_resp14(ic_peat)*lnd(i,j)%f_peat &
            + lnd(i,j)%soil_resp14(ic_shelf)*lnd(i,j)%f_shelf &
            + lnd(i,j)%soil_resp14(ic_lake)*lnd(i,j)%f_lake &
            + lnd(i,j)%soil_resp14(ic_ice)*lnd(i,j)%f_ice_grd) &
            * sec_day*day_mon * area(i,j) * 1.d-12 ! PgC/mon
            global_sresp = global_sresp + sresp_ij
            global_sresp13 = global_sresp13 + sresp13_ij
            global_sresp14 = global_sresp14 + sresp14_ij
            ! inundated area excluding peatlands
            global_inund = global_inund + max(0._wp,(lnd(i,j)%f_wetland*lnd(i,j)%f_veg-lnd(i,j)%f_peat)) * area(i,j) * 1.d-12   ! mln km2
            ! wetland area
            tmp = lnd(i,j)%f_wetland*lnd(i,j)%f_veg * area(i,j) * 1.d-12   ! mln km2
            global_wet = global_wet + tmp
            ! tropical wetland area
            if (lat(j).ge.-30._wp.and.lat(j).le.30._wp) then
              global_wettrop = global_wettrop + tmp
            endif
            ! extratropical wetland area
            if (abs(lat(j)).gt.30._wp) then
              global_wetextrop = global_wetextrop + tmp
            endif
            ! peatland area
            global_peat = global_peat + lnd(i,j)%f_peat * area(i,j) * 1.d-12   ! mln km2
            global_peatpot = global_peatpot + lnd(i,j)%f_peat_pot * area(i,j) * 1.d-12   ! mln km2
            ! CH4
            tmp = (lnd(i,j)%ch4_emis_wetland*max(0._wp,lnd(i,j)%f_wetland*lnd(i,j)%f_veg-lnd(i,j)%f_peat) &  ! wetland
            + lnd(i,j)%ch4_emis_shelf*lnd(i,j)%f_shelf & ! shelf
            + lnd(i,j)%ch4_emis_lake*lnd(i,j)%f_lake & ! lakes
            + lnd(i,j)%ch4_emis_peat*lnd(i,j)%f_peat) & ! peatland
            * sec_year * area(i,j) * 1.d-9   ! TgCH4/yr
            global_ch4 = global_ch4 + tmp
            ! tropical CH4
            if (lat(j).ge.-30.and.lat(j).le.30) then
              global_ch4trop = global_ch4trop + tmp
            endif
            ! extratropical CH4
            if (abs(lat(j)).gt.30) then
              global_ch4extrop = global_ch4extrop + tmp
            endif
            ! shelf CH4
            global_ch4shelf = global_ch4shelf + lnd(i,j)%ch4_emis_shelf*lnd(i,j)%f_shelf &
            * sec_year * area(i,j) * 1.d-9   ! TgCH4/yr
            ! lake CH4
            global_ch4lake = global_ch4lake + lnd(i,j)%ch4_emis_lake*lnd(i,j)%f_lake &
            * sec_year * area(i,j) * 1.d-9   ! TgCH4/yr
            ! C13 of CH4
            global_c13h4 = global_c13h4 &
            + (lnd(i,j)%c13h4_emis_wetland*max(0._wp,lnd(i,j)%f_wetland*lnd(i,j)%f_veg-lnd(i,j)%f_peat) &  ! wetland
            + lnd(i,j)%c13h4_emis_shelf*lnd(i,j)%f_shelf & ! shelf
            + lnd(i,j)%c13h4_emis_peat*lnd(i,j)%f_peat) & ! peatland
            * sec_year * area(i,j) * 1.d-9   ! TgC13H4/yr
            do k = 1, npft
              ! PFTs area
              global_pfts(k) = global_pfts(k) + lnd(i,j)%frac_surf(k) * area(i,j) * 1.d-12  ! mln km2
              ! forest area
              if( flag_tree(k) .eq. 1 )   global_forest = global_forest + lnd(i,j)%frac_surf(k) * area(i,j) * 1.d-12  ! mln km2
              ! grass area
              if( flag_grass(k) .eq. 1 )  global_grass  = global_grass  + lnd(i,j)%frac_surf(k) * area(i,j) * 1.d-12  ! mln km2
              ! shrub area
              if( flag_shrub(k) .eq. 1 )  global_shrub  = global_shrub  + lnd(i,j)%frac_surf(k) * area(i,j) * 1.d-12  ! mln km2
              ! vegetation carbon
              global_vegc = global_vegc + lnd(i,j)%veg_c(k) * lnd(i,j)%frac_surf(k) * area(i,j) * 1.d-12  ! mln km2
              if ( pi_perm_mask(i,j).eq.1) global_vegc_pimask = global_vegc_pimask + lnd(i,j)%veg_c(k) * lnd(i,j)%frac_surf(k) * area(i,j) * 1.d-12  ! mln km2
            enddo
            ! desert area
            global_desert = global_desert + lnd(i,j)%frac_surf(i_bare) * area(i,j) * 1.d-12  ! mln km2
            ! cropland area
            global_crop = global_crop + lnd(i,j)%f_crop*lnd(i,j)%f_veg * area(i,j) * 1.d-12  ! mln km2
            ! pasture area
            global_pasture = global_pasture + lnd(i,j)%f_pasture*lnd(i,j)%f_veg * area(i,j) * 1.d-12  ! mln km2
            ! soil carbon
            do n=1,nl
              if (n.eq.1) then
                litterc = (lnd(i,j)%litter_c(n)*(lnd(i,j)%f_veg-lnd(i,j)%f_peat) &
                  + lnd(i,j)%litter_c_peat/dz(n)*lnd(i,j)%f_peat &
                  + lnd(i,j)%litter_c_ice(n)*lnd(i,j)%f_ice_grd &
                  + lnd(i,j)%litter_c_shelf(n)*lnd(i,j)%f_shelf &
                  + lnd(i,j)%litter_c_lake(n)*lnd(i,j)%f_lake) &
                  * area(i,j) ! kgC/m
                fastc = (lnd(i,j)%fast_c(n)*(lnd(i,j)%f_veg-lnd(i,j)%f_peat) &
                  + lnd(i,j)%acro_c/dz(n)*lnd(i,j)%f_peat &
                  + lnd(i,j)%fast_c_ice(n)*lnd(i,j)%f_ice_grd &
                  + lnd(i,j)%fast_c_shelf(n)*lnd(i,j)%f_shelf &
                  + lnd(i,j)%fast_c_lake(n)*lnd(i,j)%f_lake) &
                  * area(i,j) ! kgC/m
                slowc = (lnd(i,j)%slow_c(n)*(lnd(i,j)%f_veg-lnd(i,j)%f_peat) &
                  + lnd(i,j)%slow_c_ice(n)*lnd(i,j)%f_ice_grd &
                  + lnd(i,j)%slow_c_shelf(n)*lnd(i,j)%f_shelf &
                  + lnd(i,j)%slow_c_lake(n)*lnd(i,j)%f_lake) &
                  * area(i,j) ! kgC/m
              else
                litterc = (lnd(i,j)%litter_c(n)*(lnd(i,j)%f_veg-lnd(i,j)%f_peat) &
                  + lnd(i,j)%litter_c_ice(n)*lnd(i,j)%f_ice_grd &
                  + lnd(i,j)%litter_c_shelf(n)* lnd(i,j)%f_shelf &
                  + lnd(i,j)%litter_c_lake(n)*lnd(i,j)%f_lake) &
                  * area(i,j) ! kgC/m
                fastc = (lnd(i,j)%fast_c(n)*(lnd(i,j)%f_veg-lnd(i,j)%f_peat) &
                  + lnd(i,j)%fast_c_ice(n)*lnd(i,j)%f_ice_grd &
                  + lnd(i,j)%fast_c_shelf(n)*lnd(i,j)%f_shelf &
                  + lnd(i,j)%fast_c_lake(n)*lnd(i,j)%f_lake) &
                  * area(i,j) ! kgC/m
                slowc = (lnd(i,j)%slow_c(n)*(lnd(i,j)%f_veg-lnd(i,j)%f_peat) &
                  + lnd(i,j)%cato_c(n)*lnd(i,j)%f_peat &
                  + lnd(i,j)%slow_c_ice(n)*lnd(i,j)%f_ice_grd &
                  + lnd(i,j)%slow_c_shelf(n)*lnd(i,j)%f_shelf &
                  + lnd(i,j)%slow_c_lake(n)*lnd(i,j)%f_lake) &
                  * area(i,j) ! kgC/m
              endif
              minc  = (lnd(i,j)%litter_c(n) + lnd(i,j)%fast_c(n) + lnd(i,j)%slow_c(n)) * (lnd(i,j)%f_veg-lnd(i,j)%f_peat) * area(i,j) ! kgC/m
              if (n.eq.1) then
                peatc = (lnd(i,j)%litter_c_peat/dz(n) + lnd(i,j)%acro_c/dz(n) + lnd(i,j)%cato_c(n)) * lnd(i,j)%f_peat * area(i,j) ! kgC/m
              else
                peatc = lnd(i,j)%cato_c(n) * lnd(i,j)%f_peat * area(i,j) ! kgC/m 
              endif
              icec   = (lnd(i,j)%litter_c_ice(n) + lnd(i,j)%fast_c_ice(n) + lnd(i,j)%slow_c_ice(n)) * lnd(i,j)%f_ice_grd * area(i,j) ! kgC/m
              shelfc = (lnd(i,j)%litter_c_shelf(n) + lnd(i,j)%fast_c_shelf(n) + lnd(i,j)%slow_c_shelf(n)) * lnd(i,j)%f_shelf * area(i,j) ! kgC/m
              lakec = (lnd(i,j)%litter_c_lake(n) + lnd(i,j)%fast_c_lake(n) + lnd(i,j)%slow_c_lake(n)) * lnd(i,j)%f_lake * area(i,j) ! kgC/m
              ! depth integrated soil carbon
              global_litterc= global_litterc + litterc*dz(n)  ! kgC, litter carbon
              global_fastc  = global_fastc + fastc*dz(n)  ! kgC, fast soil carbon
              global_slowc  = global_slowc + slowc*dz(n)  ! kgC, slow soil carbon
              global_minc   = global_minc + minc*dz(n)  ! kgC, mineral soil carbon
              global_peatc  = global_peatc + peatc*dz(n)  ! kgC, peat carbon
              global_icec   = global_icec + icec*dz(n)  ! kgC, carbon below ice sheets
              global_shelfc = global_shelfc + shelfc*dz(n)  ! kgC, carbon below water on ocean shelf
              global_lakec  = global_lakec + lakec*dz(n)  ! kgC, carbon below lake water
              if (lat(j).gt.60._wp) global_soilc60N = global_soilc60N + minc*dz(n) + peatc*dz(n) + icec*dz(n) + shelfc*dz(n) + lakec*dz(n)
              if( lnd(i,j)%alt .gt. -1._wp ) global_permc = global_permc + (minc+peatc)*dz(n)  ! kgC, carbon in permafrost
              if (pi_perm_mask(i,j).eq.1) global_soilc_pimask = global_soilc_pimask + (minc+peatc)*dz(n)  ! kgC, carbon in permafrost
              ! top meter soil carbon
              if (z_int(n).le.1._wp) then
                global_minc1m   = global_minc1m + minc*dz(n)  ! kgC, mineral soil carbon
                global_peatc1m  = global_peatc1m + peatc*dz(n)  ! kgC, peat carbon
                global_icec1m   = global_icec1m + icec*dz(n)  ! kgC, carbon below ice sheets
                global_shelfc1m = global_shelfc1m + shelfc*dz(n)  ! kgC, carbon below water on ocean shelf
                global_lakec1m  = global_lakec1m + lakec*dz(n)  ! kgC, carbon below lake water
                if( lnd(i,j)%alt .gt. -1._wp ) global_permc1m = global_permc1m + (minc+peatc)*dz(n)  ! kgC
                if (pi_perm_mask(i,j).eq.1) global_soilc1m_pimask = global_soilc1m_pimask + (minc+peatc)*dz(n)  ! kgC, carbon in permafrost
                if (lat(j).gt.60._wp) global_soilc60N1m = global_soilc60N1m + minc*dz(n) + peatc*dz(n) + icec*dz(n) + shelfc*dz(n) + lakec*dz(n)
              else if (z_int(n).gt.1._wp .and. z_int(n-1).le.1._wp) then
                global_minc1m   = global_minc1m + minc*(dz(n)-(z_int(n)-1._wp))  ! kgC, mineral soil carbon
                global_peatc1m  = global_peatc1m + peatc*(dz(n)-(z_int(n)-1._wp))  ! kgC, peat carbon
                global_icec1m   = global_icec1m + icec*(dz(n)-(z_int(n)-1._wp))  ! kgC, carbon below ice sheets
                global_shelfc1m = global_shelfc1m + shelfc*(dz(n)-(z_int(n)-1._wp))  ! kgC, carbon below water on ocean shelf
                global_lakec1m  = global_lakec1m + lakec*(dz(n)-(z_int(n)-1._wp))  ! kgC, carbon below lake water
                if( lnd(i,j)%alt .gt. -1._wp ) global_permc1m = global_permc1m + (minc+peatc)*(dz(n)-(z_int(n)-1._wp))  ! kgC
                if (pi_perm_mask(i,j).eq.1) global_soilc1m_pimask = global_soilc1m_pimask + (minc+peatc)*(dz(n)-(z_int(n)-1._wp)) ! kgC, in permafrost
                if (lat(j).gt.60._wp) global_soilc60N1m = global_soilc60N1m + minc*(dz(n)-(z_int(n)-1._wp)) + peatc*(dz(n)-(z_int(n)-1._wp)) &
                  + icec*(dz(n)-(z_int(n)-1._wp)) + shelfc*(dz(n)-(z_int(n)-1._wp)) + lakec*(dz(n)-(z_int(n)-1._wp))
              endif
              ! top 3 meters soil carbon
              if (z_int(n).le.3._wp) then
                global_minc3m   = global_minc3m + minc*dz(n)  ! kgC, mineral soil carbon
                global_peatc3m  = global_peatc3m + peatc*dz(n)  ! kgC, peat carbon
                global_icec3m   = global_icec3m + icec*dz(n)  ! kgC, carbon below ice sheets
                global_shelfc3m = global_shelfc3m + shelfc*dz(n)  ! kgC, carbon below water on ocean shelf
                global_lakec3m  = global_lakec3m + lakec*dz(n)  ! kgC, carbon below lake water
                if( lnd(i,j)%alt .gt. -1._wp ) global_permc3m = global_permc3m + (minc+peatc)*dz(n)  ! kgC
                if (pi_perm_mask(i,j).eq.1) global_soilc3m_pimask = global_soilc3m_pimask + (minc+peatc)*dz(n)  ! kgC, carbon in permafrost
              else if (z_int(n).gt.3._wp .and. z_int(n-1).le.3._wp) then
                global_minc3m   = global_minc3m + minc*(dz(n)-(z_int(n)-3._wp))  ! kgC, mineral soil carbon
                global_peatc3m  = global_peatc3m + peatc*(dz(n)-(z_int(n)-3._wp))  ! kgC, peat carbon
                global_icec3m   = global_icec3m + icec*(dz(n)-(z_int(n)-3._wp))  ! kgC, carbon below ice sheets
                global_shelfc3m = global_shelfc3m + shelfc*(dz(n)-(z_int(n)-3._wp))  ! kgC, carbon below water on ocean shelf
                global_lakec3m  = global_lakec3m + lakec*(dz(n)-(z_int(n)-3._wp))  ! kgC, carbon below lake water
                if( lnd(i,j)%alt .gt. -1._wp ) global_permc3m = global_permc3m + (minc+peatc)*(dz(n)-(z_int(n)-3._wp))  ! kgC
                if (pi_perm_mask(i,j).eq.1) global_soilc3m_pimask = global_soilc3m_pimask + (minc+peatc)*(dz(n)-(z_int(n)-3._wp)) ! kgC, in permafrost
              endif
              ! inert carbon
              if (lnd(i,j)%theta_i(n)/max(1e-20_wp,lnd(i,j)%theta_i(n)+lnd(i,j)%theta_w(n)) .gt. 0.5_wp) then
                global_inertc = global_inertc + minc*dz(n) ! kgC, inert carbon
              endif
              if (lnd(i,j)%theta_i_shelf(n)/max(1e-20_wp,lnd(i,j)%theta_i_shelf(n)+lnd(i,j)%theta_w_shelf(n)) .gt. 0.5_wp) then
                global_inertc = global_inertc + shelfc*dz(n) ! kgC, inert carbon
              endif
              if (lnd(i,j)%theta_i_sublake(n)/max(1e-20_wp,lnd(i,j)%theta_i_sublake(n)+lnd(i,j)%theta_w_sublake(n)) .gt. 0.5_wp) then
                global_inertc = global_inertc + lakec*dz(n) ! kgC, inert carbon
              endif
              global_inertc = global_inertc + icec*dz(n) ! kgC, inert carbon
            enddo
          endif
        enddo
      enddo
      !$omp end parallel do

      ! convert from kgC to PgC
      global_litterc = global_litterc * 1.e-12
      global_fastc = global_fastc * 1.e-12
      global_slowc = global_slowc * 1.e-12
      global_soilc60N = global_soilc60N * 1.e-12
      global_soilc60N1m = global_soilc60N1m * 1.e-12
      global_minc = global_minc * 1.e-12
      global_peatc = global_peatc * 1.e-12
      global_icec = global_icec * 1.e-12
      global_shelfc = global_shelfc * 1.e-12
      global_lakec = global_lakec * 1.e-12
      global_permc = global_permc * 1.e-12
      global_soilc_pimask = global_soilc_pimask * 1.e-12
      global_inertc = global_inertc * 1.e-12
      global_minc1m = global_minc1m * 1.e-12
      global_peatc1m = global_peatc1m * 1.e-12
      global_icec1m = global_icec1m * 1.e-12
      global_shelfc1m = global_shelfc1m * 1.e-12
      global_lakec1m = global_lakec1m * 1.e-12
      global_permc1m = global_permc1m * 1.e-12
      global_soilc1m_pimask = global_soilc1m_pimask * 1.e-12
      global_minc3m = global_minc3m * 1.e-12
      global_peatc3m = global_peatc3m * 1.e-12
      global_icec3m = global_icec3m * 1.e-12
      global_shelfc3m = global_shelfc3m * 1.e-12
      global_lakec3m = global_lakec3m * 1.e-12
      global_permc3m = global_permc3m * 1.e-12
      global_soilc3m_pimask = global_soilc3m_pimask * 1.e-12

      ! assign to monthly time series type
      mon_ts(y,mon)%npp      = global_npp
      mon_ts(y,mon)%npp_pimask      = global_npp_pimask
      mon_ts(y,mon)%sresp    = global_sresp
      mon_ts(y,mon)%npp14      = global_npp14
      mon_ts(y,mon)%sresp14    = global_sresp14
      mon_ts(y,mon)%forest   = global_forest
      mon_ts(y,mon)%grass    = global_grass
      mon_ts(y,mon)%shrub    = global_shrub
      mon_ts(y,mon)%desert   = global_desert
      mon_ts(y,mon)%crop     = global_crop
      mon_ts(y,mon)%pasture  = global_pasture
      mon_ts(y,mon)%pfts     = global_pfts
      mon_ts(y,mon)%vegc     = global_vegc
      mon_ts(y,mon)%vegc_pimask     = global_vegc_pimask
      mon_ts(y,mon)%peat     = global_peat
      mon_ts(y,mon)%peatpot  = global_peatpot
      mon_ts(y,mon)%wet    = global_wet
      mon_ts(y,mon)%wettrop    = global_wettrop
      mon_ts(y,mon)%wetextrop    = global_wetextrop
      mon_ts(y,mon)%inund    = global_inund
      mon_ts(y,mon)%soilc    = global_minc + global_peatc + global_icec + global_shelfc + global_lakec
      mon_ts(y,mon)%soilc60N = global_soilc60N
      mon_ts(y,mon)%litterc  = global_litterc
      mon_ts(y,mon)%fastc    = global_fastc
      mon_ts(y,mon)%slowc    = global_slowc
      mon_ts(y,mon)%minc     = global_minc
      mon_ts(y,mon)%peatc    = global_peatc
      mon_ts(y,mon)%permc    = global_permc
      mon_ts(y,mon)%soilc_pimask    = global_soilc_pimask
      mon_ts(y,mon)%icec     = global_icec
      mon_ts(y,mon)%shelfc   = global_shelfc
      mon_ts(y,mon)%lakec    = global_lakec
      mon_ts(y,mon)%inertc   = global_inertc
      mon_ts(y,mon)%soilc1m    = global_minc1m + global_peatc1m + global_icec1m + global_shelfc1m + global_lakec1m
      mon_ts(y,mon)%soilc60N1m = global_soilc60N1m 
      mon_ts(y,mon)%minc1m     = global_minc1m
      mon_ts(y,mon)%peatc1m    = global_peatc1m
      mon_ts(y,mon)%permc1m    = global_permc1m
      mon_ts(y,mon)%soilc1m_pimask    = global_soilc1m_pimask
      mon_ts(y,mon)%icec1m     = global_icec1m
      mon_ts(y,mon)%shelfc1m   = global_shelfc1m
      mon_ts(y,mon)%lakec1m    = global_lakec1m
      mon_ts(y,mon)%soilc3m    = global_minc3m + global_peatc3m + global_icec3m + global_shelfc3m + global_lakec3m
      mon_ts(y,mon)%minc3m     = global_minc3m
      mon_ts(y,mon)%peatc3m    = global_peatc3m
      mon_ts(y,mon)%permc3m    = global_permc3m
      mon_ts(y,mon)%soilc3m_pimask    = global_soilc3m_pimask
      mon_ts(y,mon)%icec3m     = global_icec3m
      mon_ts(y,mon)%shelfc3m   = global_shelfc3m
      mon_ts(y,mon)%lakec3m    = global_lakec3m
      mon_ts(y,mon)%ch4      = global_ch4
      mon_ts(y,mon)%ch4trop   = global_ch4trop
      mon_ts(y,mon)%ch4shelf   = global_ch4shelf
      mon_ts(y,mon)%ch4lake    = global_ch4lake
      mon_ts(y,mon)%ch4extrop  = global_ch4extrop
      if( global_ch4 .gt. 0._wp ) then
        mon_ts(y,mon)%d13ch4   = (global_c13h4/global_ch4 / c13_c12_std - 1._wp ) * 1000._wp  ! permil
      else
        mon_ts(y,mon)%d13ch4   = 0._wp
      endif

    endif ! time_eoy_lnd

    ! at end of year
    if( time_eoy_lnd ) then

      ! atmospheric CO2
      ann_ts(y)%co2 = lndp%co2
      ! land area
      ann_ts(y)%land = sum(lnd%f_land * area) * 1.d-12 ! mln km2
      ! vegetated area (including desert)
      ann_ts(y)%veg = sum(lnd%f_veg * area) * 1.d-12 ! mln km2
      ! ice sheet covered area
      ann_ts(y)%ice = sum(lnd%f_ice * area) * 1.d-12 ! mln km2
      ! flooded shelf area
      ann_ts(y)%shelf = sum(lnd%f_shelf * area) * 1.d-12 ! mln km2
      ! lake area
      ann_ts(y)%lake = sum(lnd%f_lake * area) * 1.d-12 ! mln km2
      ! permafrost area
      ann_ts(y)%perm = sum(lnd%f_veg * area,mask=lnd%alt.gt.-1._wp) * 1.d-12   ! mln km2
      ! organic carbon weathering
      ann_ts(y)%doc_export = sum(lnd%doc_export * lnd%f_veg*area) * 1.d-12    ! PgC/yr
      ann_ts(y)%poc_export = sum(lnd%poc_export * lnd%f_veg*area) * 1.d-12    ! PgC/yr
      ! carbonate weathering
      ann_ts(y)%weath_carb = sum(lnd%weath_carb * lnd%f_veg*area) * 1.d-12    ! Tmol C/yr
      ! silicate weathering
      ann_ts(y)%weath_sil = sum(lnd%weath_sil * lnd%f_veg*area) * 1.d-12    ! Tmol C/yr
      ! loess carbonate weathering
      ann_ts(y)%weath_loess = sum(lnd%weath_loess * lnd%f_veg*area) * 1.d-12    ! Tmol C/yr

      ! cumulate/average monthly to annual
      call ts_ave( mon_ts(y,:), ann_ts(y) )

      ! total land carbon
      ann_ts(y)%landc = ann_ts(y)%soilc + ann_ts(y)%vegc 
      ann_ts(y)%Cflx_atm_lnd   = lndp%Cflx_atm_lnd*1.e-12 ! PgC/yr
      ann_ts(y)%C13flx_atm_lnd = lndp%C13flx_atm_lnd*1.e-12 ! PgC/yr 
      ann_ts(y)%C14flx_atm_lnd = lndp%C14flx_atm_lnd*1.e-12 ! PgC/yr

      ! air-land carbon fluxes
      if (abs(ann_ts(y)%Cflx_atm_lnd).gt.1.e-10) then
        ann_ts(y)%d13Cflx_atm_lnd = (ann_ts(y)%C13flx_atm_lnd/ann_ts(y)%Cflx_atm_lnd / c13_c12_std - 1._wp ) * 1000._wp
        ann_ts(y)%D14Cflx_atm_lnd = (ann_ts(y)%C14flx_atm_lnd/ann_ts(y)%Cflx_atm_lnd*(0.975_wp/(1._wp+ann_ts(y)%d13Cflx_atm_lnd/1000._wp))**2 &
          / c14_c_std - 1._wp) * 1000._wp 
      else
        ann_ts(y)%d13Cflx_atm_lnd = missing_value
        ann_ts(y)%D14Cflx_atm_lnd = missing_value
      endif

      ann_ts(y)%Cflx_burial  = lndp%Cflx_burial*1.e-12 ! PgC/yr
      if (abs(lndp%Cflx_burial).gt.1.e-10_wp) then
        ann_ts(y)%d13Cflx_burial = ((lndp%C13flx_burial/lndp%Cflx_burial)/c13_c12_std-1._wp) * 1000._wp
      else
        ann_ts(y)%d13Cflx_burial = missing_value 
      endif

      ! write time series of global annual variables to file
      if (time_out_ts_clim) then
        call ts_nc_write(trim(out_dir)//"/lnd_ts.nc",ann_ts(1:y),year_clim-y+1,y)
      endif

      ! write to standard output
      if (mod(year,10).eq.1) then
        print '(a7,a9,17a7)','lnd','year','Cflx','NPP','CH4e','Aperm','prc','evp','run','Awet','Apeat', &
        'Cveg','Csoil','Cperm','Cpeat','Cice','Clake','Cshelf','dust_e'
      endif
      print '(a7,i9,F7.3,10F7.1,6F7.0)', &
      'lnd',year_now,ann_ts(y)%Cflx_atm_lnd,ann_ts(y)%npp,ann_ts(y)%ch4,ann_ts(y)%perm, &
      ann_ts(y)%prc,ann_ts(y)%evp,ann_ts(y)%runoff,ann_ts(y)%wet,ann_ts(y)%peat,ann_ts(y)%vegc,ann_ts(y)%soilc, &
      ann_ts(y)%permc,ann_ts(y)%peatc,ann_ts(y)%icec,ann_ts(y)%lakec,ann_ts(y)%shelfc,ann_ts(y)%dust_e

    endif ! time_eoy_lnd


    !-------------------------------------------------------
    ! 2D-3D output
    !-------------------------------------------------------

    if( time_out_lnd ) then

      if( time_soy_lnd ) then
        do m=1,nmon_year

          if( write_surf ) then
            do i=1,nx
              do j=1,ny
                if (lnd(i,j)%f_land.gt.0._wp) then
                  mon_su(m)%le(i,j,:)       = 0._wp
                  mon_su(m)%etot(i,j,:)     = 0._wp
                  mon_su(m)%ecan(i,j,:)     = 0._wp
                  mon_su(m)%esur(i,j,:)     = 0._wp
                  mon_su(m)%trans(i,j,:)    = 0._wp
                  mon_su(m)%sh(i,j,:)       = 0._wp
                  mon_su(m)%slh(i,j,:)      = 0._wp
                  mon_su(m)%ef(i,j,:)       = 0._wp
                  mon_su(m)%swnet(i,j,:)    = 0._wp
                  mon_su(m)%lwnet(i,j,:)    = 0._wp
                  mon_su(m)%g(i,j,:)        = 0._wp
                  mon_su(m)%flx_melt(i,j,:) = 0._wp
                  mon_su(m)%tskin(i,j,:)    = 0._wp
                  mon_su(m)%tskin_amp(i,j,:)= 0._wp
                  mon_su(m)%tlake_sur(i,j)  = 0._wp
                  mon_su(m)%alb(i,j,:)      = 0._wp
                  mon_su(m)%alb_dir(i,j,:)  = 0._wp
                  mon_su(m)%alb_dif(i,j,:)  = 0._wp
                  mon_su(m)%albsnw(i,j,:)   = 0._wp
                  mon_su(m)%scf(i,j,:)      = 0._wp
                  mon_su(m)%fcanwat(i,j,:)  = 0._wp
                  mon_su(m)%fcansn(i,j,:)   = 0._wp
                  mon_su(m)%swe(i,j,:)      = 0._wp
                  mon_su(m)%swe_max(i,j,:)  = 0._wp
                  mon_su(m)%hsnow(i,j,:)    = 0._wp
                  mon_su(m)%tsnow(i,j,:)    = 0._wp
                  mon_su(m)%wind(i,j,:)     = 0._wp
                  mon_su(m)%ra(i,j,:)       = 0._wp
                  mon_su(m)%rag(i,j,:)      = 0._wp
                  mon_su(m)%Ri(i,j,:)       = 0._wp
                  mon_su(m)%z0m(i,j,:)      = 0._wp
                  mon_su(m)%z0h(i,j,:)      = 0._wp
                  mon_su(m)%Ch(i,j,:)       = 0._wp
                  mon_su(m)%rs(i,j,:)       = 0._wp
                  mon_su(m)%betas(i,j,:)    = 0._wp
                  mon_su(m)%cti_lim(i,j)    = 0._wp
                  mon_su(m)%fwet(i,j)       = 0._wp
                  mon_su(m)%wtab(i,j)       = 0._wp
                  mon_su(m)%inf(i,j)        = 0._wp
                  mon_su(m)%runoff(i,j,:)   = 0._wp
                  mon_su(m)%runsur(i,j,:)   = 0._wp
                  mon_su(m)%calving(i,j,:)  = 0._wp
                  mon_su(m)%drain(i,j,:)    = 0._wp
                  mon_su(m)%snowmelt(i,j,:) = 0._wp
                  mon_su(m)%icemelt(i,j,:)  = 0._wp
                  mon_su(m)%vpd(i,j)        = 0._wp
                  mon_su(m)%dust_e_d(i,j)   = 0._wp
                  mon_su(m)%dust_e_g(i,j)   = 0._wp
                  mon_su(m)%dust_e_s(i,j)   = 0._wp
                  mon_su(m)%dust_e(i,j)     = 0._wp
                  mon_su(m)%dust_dep(i,j)   = 0._wp
                  mon_su(m)%snow_grain(i,j,:) = 0._wp
                  mon_su(m)%dust_con(i,j,:) = 0._wp
                else
                  mon_su(m)%le(i,j,:)       = missing_value 
                  mon_su(m)%etot(i,j,:)     = missing_value 
                  mon_su(m)%ecan(i,j,:)     = missing_value 
                  mon_su(m)%esur(i,j,:)     = missing_value 
                  mon_su(m)%trans(i,j,:)    = missing_value 
                  mon_su(m)%sh(i,j,:)       = missing_value 
                  mon_su(m)%slh(i,j,:)      = missing_value 
                  mon_su(m)%ef(i,j,:)       = missing_value 
                  mon_su(m)%swnet(i,j,:)    = missing_value 
                  mon_su(m)%lwnet(i,j,:)    = missing_value 
                  mon_su(m)%g(i,j,:)        = missing_value 
                  mon_su(m)%flx_melt(i,j,:) = missing_value 
                  mon_su(m)%tskin(i,j,:)    = missing_value 
                  mon_su(m)%tskin_amp(i,j,:)= missing_value 
                  mon_su(m)%tlake_sur(i,j)  = missing_value 
                  mon_su(m)%alb(i,j,:)      = missing_value 
                  mon_su(m)%alb_dir(i,j,:)  = missing_value 
                  mon_su(m)%alb_dif(i,j,:)  = missing_value 
                  mon_su(m)%albsnw(i,j,:)   = missing_value 
                  mon_su(m)%scf(i,j,:)      = missing_value 
                  mon_su(m)%fcanwat(i,j,:)  = missing_value 
                  mon_su(m)%fcansn(i,j,:)   = missing_value 
                  mon_su(m)%swe(i,j,:)      = missing_value 
                  mon_su(m)%swe_max(i,j,:)  = missing_value 
                  mon_su(m)%hsnow(i,j,:)    = missing_value 
                  mon_su(m)%tsnow(i,j,:)    = missing_value 
                  mon_su(m)%wind(i,j,:)     = missing_value 
                  mon_su(m)%ra(i,j,:)       = missing_value 
                  mon_su(m)%rag(i,j,:)      = missing_value 
                  mon_su(m)%Ri(i,j,:)       = missing_value 
                  mon_su(m)%z0m(i,j,:)      = missing_value 
                  mon_su(m)%z0h(i,j,:)      = missing_value 
                  mon_su(m)%Ch(i,j,:)       = missing_value 
                  mon_su(m)%rs(i,j,:)       = missing_value 
                  mon_su(m)%betas(i,j,:)    = missing_value 
                  mon_su(m)%cti_lim(i,j)    = missing_value 
                  mon_su(m)%fwet(i,j)       = missing_value 
                  mon_su(m)%wtab(i,j)       = missing_value 
                  mon_su(m)%inf(i,j)        = missing_value 
                  mon_su(m)%runoff(i,j,:)   = missing_value 
                  mon_su(m)%runsur(i,j,:)   = missing_value 
                  mon_su(m)%calving(i,j,:)  = missing_value 
                  mon_su(m)%drain(i,j,:)    = missing_value 
                  mon_su(m)%snowmelt(i,j,:) = missing_value 
                  mon_su(m)%icemelt(i,j,:)  = missing_value 
                  mon_su(m)%vpd(i,j)        = missing_value 
                  mon_su(m)%dust_e_d(i,j)   = missing_value 
                  mon_su(m)%dust_e_g(i,j)   = missing_value 
                  mon_su(m)%dust_e_s(i,j)   = missing_value 
                  mon_su(m)%dust_e(i,j)     = missing_value 
                  mon_su(m)%dust_dep(i,j)   = missing_value 
                  mon_su(m)%snow_grain(i,j,:) = missing_value
                  mon_su(m)%dust_con(i,j,:) = missing_value 
                endif
              enddo
            enddo
          endif

          if( write_carbon ) then
            mon_c(m)%lai   = 0._wp
            mon_c(m)%xi    = 0._wp
            mon_c(m)%wue   = 0._wp
            mon_c(m)%gcan  = 0._wp
            mon_c(m)%gpp   = 0._wp
            mon_c(m)%npp   = 0._wp
            mon_c(m)%npp13 = 0._wp
            mon_c(m)%npp14 = 0._wp
            mon_c(m)%aresp   = 0._wp
            mon_c(m)%disc    = 0._wp
          endif

          if( write_soil ) then
            mon_s(m)%tsoil  = 0._wp
            mon_s(m)%tice   = 0._wp
            mon_s(m)%tsublake  = 0._wp
            mon_s(m)%thetaw = 0._wp
            mon_s(m)%thetai = 0._wp
            mon_s(m)%theta  = 0._wp
            mon_s(m)%thetas = 0._wp
            mon_s(m)%fthetas = 0._wp
          endif

          if( write_soil_par ) then
            mon_sp(m)%lambda_if = 0._wp
            mon_sp(m)%cap_if = 0._wp
            mon_sp(m)%lambda_i = 0._wp
            mon_sp(m)%cap_i  = 0._wp
            mon_sp(m)%kappa  = 0._wp
            mon_sp(m)%ksat   = 0._wp
            mon_sp(m)%psisat = 0._wp
            mon_sp(m)%bi     = 0._wp
            mon_sp(m)%theta_sat = 0._wp
            mon_sp(m)%theta_field = 0._wp
            mon_sp(m)%theta_wilt = 0._wp
            mon_sp(m)%ftemp  = 0._wp
            mon_sp(m)%fmoist = 0._wp
            mon_sp(m)%fdepth = 0._wp
          endif

          if ( write_lake ) then
            mon_l(m)%lakebal  = 0._wp
            mon_l(m)%t_lake   = 0._wp
            mon_l(m)%lambda_lake   = 0._wp
            mon_l(m)%f_i_lake = 0._wp
            mon_l(m)%f_lake_ice = 0._wp
            mon_l(m)%h_lake_conv  = 0._wp
            mon_l(m)%h_lake_mix   = 0._wp
          endif

          if( write_cons ) then
            mon_b(m)%econs_su1 = 0._wp
            mon_b(m)%econs_su2 = 0._wp
            mon_b(m)%econs_so  = 0._wp
            mon_b(m)%wcons     = 0._wp
          endif

        enddo
      endif

      if( write_surf ) then
        if (l_daily_output) then
          !!$omp parallel do private(i,j,n)
          do i=1,nx
            do j=1,ny
                day_su(doy)%frac_surf(i,j,:)= lnd(i,j)%frac_surf    
                day_su(doy)%le(i,j,:)       = lnd(i,j)%flx_lh    
                day_su(doy)%sh(i,j,:)       = lnd(i,j)%flx_sh
                day_su(doy)%lwnet(i,j,:)    = lnd(i,j)%lwnet 
                day_su(doy)%swnet(i,j,:)    = lnd(i,j)%swnet  
                day_su(doy)%g(i,j,:)        = lnd(i,j)%flx_g   
                day_su(doy)%flx_melt(i,j,:) = lnd(i,j)%flx_melt
                day_su(doy)%tskin(i,j,:)    = lnd(i,j)%t_skin   
                day_su(doy)%ra(i,j,:)       = lnd(i,j)%r_a   
                day_su(doy)%rs(i,j,:)       = lnd(i,j)%r_s   
                day_su(doy)%betas(i,j,:)    = lnd(i,j)%beta_s   
                day_su(doy)%rag(i,j,:)      = lnd(i,j)%r_a_can   
                day_su(doy)%Ri(i,j,:)       = lnd(i,j)%Ri
                day_su(doy)%ecan(i,j,:)     = (lnd(i,j)%evap_can+lnd(i,j)%subl_can)*sec_day ! mm/day
                day_su(doy)%msnow(i,j,:)    = lnd(i,j)%mask_snow
                day_su(doy)%hsnow(i,j,:)    = lnd(i,j)%h_snow
                day_su(doy)%swe(i,j,:)      = lnd(i,j)%w_snow
                day_su(doy)%swe_max(i,j,:)  = lnd(i,j)%w_snow_max
                day_su(doy)%tsnow(i,j,1)    = lnd(i,j)%t_soil(0) 
                day_su(doy)%tsnow(i,j,2)    = lnd(i,j)%t_ice(0) 
                day_su(doy)%tsnow(i,j,3)    = lnd(i,j)%t_lake(0) 
                day_su(doy)%albsnw(i,j,:)   = & 
                  + frac_vu         * 0.5_wp*(lnd(i,j)%alb_snow_vis_dir+lnd(i,j)%alb_snow_vis_dif)  &
                  + (1._wp-frac_vu) * 0.5_wp*(lnd(i,j)%alb_snow_nir_dir+lnd(i,j)%alb_snow_nir_dif) 
            enddo
          enddo
          !!$omp end parallel do
        endif
        !!$omp parallel do private(i,j,n)
        do i=1,nx
          do j=1,ny
            where (lnd(i,j)%frac_surf.gt.0._wp)
              mon_su(mon)%le(i,j,:)       = mon_su(mon)%le(i,j,:)       + lnd(i,j)%flx_lh     * mon_avg
              mon_su(mon)%etot(i,j,:)     = mon_su(mon)%etot(i,j,:)     + lnd(i,j)%et*sec_day           * mon_avg ! mm/day
              mon_su(mon)%ecan(i,j,:)     = mon_su(mon)%ecan(i,j,:)     + (lnd(i,j)%evap_can+lnd(i,j)%subl_can)*sec_day  * mon_avg ! mm/day
              mon_su(mon)%esur(i,j,:)     = mon_su(mon)%esur(i,j,:)     + lnd(i,j)%evap_surface*sec_day * mon_avg ! mm/day
              mon_su(mon)%trans(i,j,:)    = mon_su(mon)%trans(i,j,:)    + lnd(i,j)%transpiration*sec_day* mon_avg ! mm/day
              mon_su(mon)%sh(i,j,:)       = mon_su(mon)%sh(i,j,:)       + lnd(i,j)%flx_sh   * mon_avg
              mon_su(mon)%ef(i,j,:)       = mon_su(mon)%ef(i,j,:)       + min(1._wp,max(0._wp,lnd(i,j)%flx_lh/(lnd(i,j)%flx_lh+lnd(i,j)%flx_sh))) * mon_avg
              mon_su(mon)%slh(i,j,:)      = mon_su(mon)%slh(i,j,:)      + (lnd(i,j)%flx_lh+lnd(i,j)%flx_sh)   * mon_avg
              mon_su(mon)%lwnet(i,j,:)    = mon_su(mon)%lwnet(i,j,:)    + lnd(i,j)%lwnet                * mon_avg
              mon_su(mon)%swnet(i,j,:)    = mon_su(mon)%swnet(i,j,:)    + lnd(i,j)%swnet                * mon_avg
              mon_su(mon)%g(i,j,:)        = mon_su(mon)%g(i,j,:)        + lnd(i,j)%flx_g     * mon_avg
              mon_su(mon)%flx_melt(i,j,:) = mon_su(mon)%flx_melt(i,j,:) + lnd(i,j)%flx_melt  * mon_avg
              mon_su(mon)%tskin(i,j,:)    = mon_su(mon)%tskin(i,j,:)    + lnd(i,j)%t_skin               * mon_avg
              mon_su(mon)%tskin_amp(i,j,:)= mon_su(mon)%tskin_amp(i,j,:)+ lnd(i,j)%t_skin_amp           * mon_avg
              mon_su(mon)%scf(i,j,:)      = mon_su(mon)%scf(i,j,:)      + lnd(i,j)%f_snow               * mon_avg
              mon_su(mon)%fcanwat(i,j,:)  = mon_su(mon)%fcanwat(i,j,:)  + lnd(i,j)%f_wat_can            * mon_avg
              mon_su(mon)%fcansn(i,j,:)   = mon_su(mon)%fcansn(i,j,:)   + lnd(i,j)%f_snow_can           * mon_avg
              mon_su(mon)%alb(i,j,:)      = mon_su(mon)%alb(i,j,:)      + lnd(i,j)%albedo               * mon_avg
              mon_su(mon)%alb_dir(i,j,:)  = mon_su(mon)%alb_dir(i,j,:)  &
                + (frac_vu*lnd(i,j)%alb_vis_dir+(1._wp-frac_vu)*lnd(i,j)%alb_nir_dir)* mon_avg
              mon_su(mon)%alb_dif(i,j,:)  = mon_su(mon)%alb_dif(i,j,:)  &
                + (frac_vu*lnd(i,j)%alb_vis_dif+(1._wp-frac_vu)*lnd(i,j)%alb_nir_dif)* mon_avg
              mon_su(mon)%wind(i,j,:)     = mon_su(mon)%wind(i,j,:)     + lnd(i,j)%wind                 * mon_avg
              mon_su(mon)%ra(i,j,:)       = mon_su(mon)%ra(i,j,:)       + lnd(i,j)%r_a                  * mon_avg
              mon_su(mon)%Ri(i,j,:)       = mon_su(mon)%Ri(i,j,:)       + lnd(i,j)%Ri                   * mon_avg
              mon_su(mon)%z0m(i,j,:)      = mon_su(mon)%z0m(i,j,:)      + lnd(i,j)%rough_m              * mon_avg
              mon_su(mon)%z0h(i,j,:)      = mon_su(mon)%z0h(i,j,:)      + lnd(i,j)%rough_h              * mon_avg
              mon_su(mon)%Ch(i,j,:)       = mon_su(mon)%Ch(i,j,:)       + lnd(i,j)%Ch                   * mon_avg
              mon_su(mon)%rs(i,j,:)       = mon_su(mon)%rs(i,j,:)       + lnd(i,j)%r_s                  * mon_avg
              mon_su(mon)%betas(i,j,:)    = mon_su(mon)%betas(i,j,:)    + lnd(i,j)%beta_s               * mon_avg
            endwhere
            do n=1,npft
              if (lnd(i,j)%frac_surf(n).gt.0._wp) then
                mon_su(mon)%rag(i,j,n)      = mon_su(mon)%rag(i,j,n)      + lnd(i,j)%r_a_can(n)              * mon_avg
              endif
            enddo
            if (lnd(i,j)%f_land.gt.0._wp) then
              mon_su(mon)%albsnw(i,j,:)   = mon_su(mon)%albsnw(i,j,:)   &
                + frac_vu         * 0.5_wp*(lnd(i,j)%alb_snow_vis_dir+lnd(i,j)%alb_snow_vis_dif)       * mon_avg &
                + (1._wp-frac_vu) * 0.5_wp*(lnd(i,j)%alb_snow_nir_dir+lnd(i,j)%alb_snow_nir_dif)       * mon_avg
              mon_su(mon)%swe(i,j,:)      = mon_su(mon)%swe(i,j,:)      + lnd(i,j)%w_snow                  * mon_avg
              mon_su(mon)%swe_max(i,j,:)  = mon_su(mon)%swe_max(i,j,:)  + lnd(i,j)%w_snow_max              * mon_avg
              mon_su(mon)%hsnow(i,j,:)    = mon_su(mon)%hsnow(i,j,:)    + lnd(i,j)%h_snow               * mon_avg
              mon_su(mon)%tsnow(i,j,1)    = mon_su(mon)%tsnow(i,j,1)    + lnd(i,j)%t_soil(0)            * mon_avg
              mon_su(mon)%tsnow(i,j,2)    = mon_su(mon)%tsnow(i,j,2)    + lnd(i,j)%t_ice(0)            * mon_avg
              mon_su(mon)%tsnow(i,j,3)    = mon_su(mon)%tsnow(i,j,3)    + lnd(i,j)%t_lake(0)            * mon_avg
              mon_su(mon)%tlake_sur(i,j)  = mon_su(mon)%tlake_sur(i,j)  + lnd(i,j)%t_lake(1) * mon_avg
              mon_su(mon)%wtab(i,j)       = mon_su(mon)%wtab(i,j)       + lnd(i,j)%w_table              * mon_avg
              mon_su(mon)%inf(i,j)        = mon_su(mon)%inf(i,j)        + lnd(i,j)%infiltration*sec_day * mon_avg
              mon_su(mon)%runoff(i,j,:)   = mon_su(mon)%runoff(i,j,:)   + (lnd(i,j)%runoff_sur+lnd(i,j)%drainage)*sec_day     * mon_avg
              mon_su(mon)%runsur(i,j,:)   = mon_su(mon)%runsur(i,j,:)   + lnd(i,j)%runoff_sur*sec_day     * mon_avg
              mon_su(mon)%calving(i,j,:)  = mon_su(mon)%calving(i,j,:)  + lnd(i,j)%calving*sec_day     * mon_avg
              mon_su(mon)%drain(i,j,:)    = mon_su(mon)%drain(i,j,:)    + lnd(i,j)%drainage*sec_day     * mon_avg
              mon_su(mon)%snowmelt(i,j,:) = mon_su(mon)%snowmelt(i,j,:) + lnd(i,j)%snowmelt*sec_day     * mon_avg
              mon_su(mon)%icemelt(i,j,:)  = mon_su(mon)%icemelt(i,j,:)  + lnd(i,j)%icemelt*sec_day      * mon_avg
              mon_su(mon)%vpd(i,j)        = mon_su(mon)%vpd(i,j)        + (e_sat_w(lnd(i,j)%t2m(1))-q_to_e(lnd(i,j)%q2m(1),lnd(i,j)%pressure(1))) * mon_avg
              mon_su(mon)%snow_grain(i,j,:) = mon_su(mon)%snow_grain(i,j,:) + lnd(i,j)%snow_grain * mon_avg
              mon_su(mon)%dust_con(i,j,:) = mon_su(mon)%dust_con(i,j,:) + lnd(i,j)%dust_con*1.e6 * mon_avg  ! mg/kg
              mon_su(mon)%dust_e_d(i,j)   = mon_su(mon)%dust_e_d(i,j)   + lnd(i,j)%dust_emis_d                * 1000._wp*dt ! g/m2/mon
              mon_su(mon)%dust_e_g(i,j)   = mon_su(mon)%dust_e_g(i,j)   + lnd(i,j)%dust_emis_g * 1000._wp*dt
              mon_su(mon)%dust_e_s(i,j)   = mon_su(mon)%dust_e_s(i,j)   + lnd(i,j)%dust_emis_s * 1000._wp*dt
              mon_su(mon)%dust_e(i,j)     = mon_su(mon)%dust_e(i,j)     + lnd(i,j)%dust_emis * 1000._wp*dt
              mon_su(mon)%dust_dep(i,j)   = mon_su(mon)%dust_dep(i,j)   + lnd(i,j)%dust_dep * 1000._wp*dt
            endif
          enddo
        enddo
        !!$omp end parallel do
        if( time_eom_lnd ) then
          where (lnd%f_land.gt.0._wp) 
            mon_su(mon)%fwet    = lnd%f_wetland
            mon_su(mon)%cti_lim = lnd%cti_lim
          endwhere
          if(time_eoy_lnd) then
            where (lnd%f_land.gt.0._wp) 
              ann_su%fwetmaxpot = lnd%f_wet_max
            elsewhere 
              ann_su%fwetmaxpot = 0._wp
            endwhere
            do i=1,nx
              do j=1,ny 
                ann_su%fwetmax(i,j) = mon_su(1)%fwet(i,j)
                do m=1,nmon_year
                  ann_su%fwetmax(i,j) = max(ann_su%fwetmax(i,j),mon_su(m)%fwet(i,j))
                enddo
              enddo
            enddo
          endif
          where ( (mon_su(mon)%drain+mon_su(mon)%runsur) .gt. 0 ) 
            mon_su(mon)%rsursub = mon_su(mon)%runsur/(mon_su(mon)%drain+mon_su(mon)%runsur)
          elsewhere
            mon_su(mon)%rsursub = 0._wp
          endwhere
        endif

      endif

      if( write_carbon ) then
        !!$omp parallel do private(i,j)
        do i=1,nx
          do j=1,ny
            if (lnd(i,j)%f_land.gt.0._wp) then
              mon_c(mon)%lai(i,j,:)   = mon_c(mon)%lai(i,j,:)   + lnd(i,j)%lai            * mon_avg
              mon_c(mon)%xi(i,j,:)    = mon_c(mon)%xi(i,j,:)    + lnd(i,j)%ci*1.e6_wp/lndp%co2    * mon_avg 
              mon_c(mon)%wue(i,j,:)   = mon_c(mon)%wue(i,j,:)   + (lndp%co2-lnd(i,j)%ci*1.e6_wp)/1.6_wp    * mon_avg 
              mon_c(mon)%gcan(i,j,:)  = mon_c(mon)%gcan(i,j,:)  + lnd(i,j)%g_can          * mon_avg
              mon_c(mon)%gpp(i,j,:)   = mon_c(mon)%gpp(i,j,:)   + lnd(i,j)%gpp*real(nday_year,wp)   * mon_avg  ! kgC/m2/yr
              mon_c(mon)%npp(i,j,:)   = mon_c(mon)%npp(i,j,:)   + lnd(i,j)%npp*real(nday_year,wp)   * mon_avg  ! kgC/m2/yr
              mon_c(mon)%npp13(i,j,:) = mon_c(mon)%npp13(i,j,:) + lnd(i,j)%npp13*real(nday_year,wp)   * mon_avg  ! kgC/m2/yr
              mon_c(mon)%npp14(i,j,:) = mon_c(mon)%npp14(i,j,:) + lnd(i,j)%npp14*real(nday_year,wp)   * mon_avg  ! kgC/m2/yr
              mon_c(mon)%aresp(i,j,:) = mon_c(mon)%aresp(i,j,:) + lnd(i,j)%aresp*real(nday_year,wp)   * mon_avg  ! kgC/m2/yr
              mon_c(mon)%disc(i,j,:)  = mon_c(mon)%disc(i,j,:)  + lnd(i,j)%discrimination*lnd(i,j)%gpp*real(nday_year,wp) * mon_avg  ! permil*kgC/m2/s
            endif
          enddo
        enddo
        !!$omp end parallel do
      endif

      if( write_soil_par ) then
        !!$omp parallel do private(i,j)
        do i=1,nx
          do j=1,ny
            if (lnd(i,j)%f_land.gt.0._wp) then
              mon_sp(mon)%lambda_if(i,j,:) = mon_sp(mon)%lambda_if(i,j,:) + lnd(i,j)%lambda_soil      * mon_avg
              mon_sp(mon)%cap_if(i,j,:)    = mon_sp(mon)%cap_if(i,j,:)    + lnd(i,j)%cap_soil           * mon_avg
              mon_sp(mon)%lambda_i(i,j,:)  = mon_sp(mon)%lambda_i(i,j,:) + lnd(i,j)%lambda_ice      * mon_avg
              mon_sp(mon)%cap_i(i,j,:)     = mon_sp(mon)%cap_i(i,j,:)    + lnd(i,j)%cap_ice           * mon_avg
              mon_sp(mon)%kappa(i,j,:)  = mon_sp(mon)%kappa(i,j,:)  + lnd(i,j)%kappa_int*sec_day       * mon_avg
              mon_sp(mon)%ksat(i,j,:)  = mon_sp(mon)%ksat(i,j,:)  + lnd(i,j)%k_sat*sec_day       * mon_avg
              mon_sp(mon)%psisat(i,j,:)  = mon_sp(mon)%psisat(i,j,:)  + lnd(i,j)%psi_sat       * mon_avg
              mon_sp(mon)%bi(i,j,:)  = mon_sp(mon)%bi(i,j,:)  - lnd(i,j)%psi_exp       * mon_avg
              mon_sp(mon)%theta_sat(i,j,:) = mon_sp(mon)%theta_sat(i,j,:) + lnd(i,j)%theta_sat       * mon_avg
              mon_sp(mon)%theta_field(i,j,:) = mon_sp(mon)%theta_field(i,j,:) + lnd(i,j)%theta_field       * mon_avg
              mon_sp(mon)%theta_wilt(i,j,:) = mon_sp(mon)%theta_wilt(i,j,:) + lnd(i,j)%theta_wilt       * mon_avg
            endif
          enddo
        enddo
        !!omp end parallel do
      endif

      if( write_soil ) then
        !!$omp parallel do private(i,j)
        do i=1,nx
          do j=1,ny
            if (lnd(i,j)%f_land.gt.0._wp) then
              mon_s(mon)%tsoil(i,j,:)  = mon_s(mon)%tsoil(i,j,:)  + lnd(i,j)%t_soil(1:nl)      * mon_avg
              mon_s(mon)%tice(i,j,:)  = mon_s(mon)%tice(i,j,:)  + lnd(i,j)%t_ice(1:nl)      * mon_avg
              mon_s(mon)%tsublake(i,j,:)  = mon_s(mon)%tsublake(i,j,:)  + lnd(i,j)%t_sublake(1:nl)      * mon_avg
              mon_s(mon)%thetaw(i,j,:) = mon_s(mon)%thetaw(i,j,:) + lnd(i,j)%theta_w                 * mon_avg
              mon_s(mon)%thetai(i,j,:) = mon_s(mon)%thetai(i,j,:) + lnd(i,j)%theta_i                 * mon_avg
              mon_s(mon)%theta(i,j,:)  = mon_s(mon)%theta(i,j,:)  + lnd(i,j)%theta                   * mon_avg
              mon_s(mon)%thetas(i,j,:) = mon_s(mon)%thetas(i,j,:) + lnd(i,j)%theta_sat       * mon_avg
            endif
          enddo
        enddo
        !!$omp end parallel do
        if( time_eom_lnd ) then 
          do k=1,nl
            where (lnd%f_land.gt.0._wp) 
              mon_s(mon)%fthetas(:,:,k) = mon_s(mon)%theta(:,:,k)/mon_s(mon)%thetas(:,:,k)
            endwhere
          enddo
        endif
      endif

      if (write_lake) then
        do i=1,nx
          do j=1,ny
            if (lnd(i,j)%f_land.gt.0._wp) then
              mon_l(mon)%lakebal(i,j)    = mon_l(mon)%lakebal(i,j)    + lnd(i,j)%lake_water_tendency*sec_day              * mon_avg
              mon_l(mon)%t_lake(i,j,:)   = mon_l(mon)%t_lake(i,j,:)   + lnd(i,j)%t_lake(1:nl_l)      * mon_avg
              mon_l(mon)%lambda_lake(i,j,:)   = mon_l(mon)%lambda_lake(i,j,:)   + lnd(i,j)%lambda_lake(1:nl_l)      * mon_avg
              mon_l(mon)%f_i_lake(i,j,:) = mon_l(mon)%f_i_lake(i,j,:) + lnd(i,j)%f_i_lake(1:nl_l)      * mon_avg
              mon_l(mon)%f_lake_ice(i,j) = mon_l(mon)%f_lake_ice(i,j) + lnd(i,j)%f_lake_ice      * mon_avg
              mon_l(mon)%h_lake_conv(i,j)  = mon_l(mon)%h_lake_conv(i,j)  + lnd(i,j)%h_lake_conv     * mon_avg
              mon_l(mon)%h_lake_mix(i,j)   = mon_l(mon)%h_lake_mix(i,j)   + lnd(i,j)%h_lake_mix      * mon_avg
            endif
          enddo
        enddo
        if (time_eoy_lnd) then
          ann_l%h_lake(:,:) = lnd(:,:)%h_lake
        endif
      endif

      if( write_cons ) then
        do i=1,nx
          do j=1,ny
            if (lnd(i,j)%f_land.gt.0._wp) then
              mon_b(mon)%econs_su1(i,j,:) = mon_b(mon)%econs_su1(i,j,:) + lnd(i,j)%energy_cons_surf1     * mon_avg
              mon_b(mon)%econs_su2(i,j,:) = mon_b(mon)%econs_su2(i,j,:) + lnd(i,j)%energy_cons_surf2     * mon_avg
              mon_b(mon)%wcons(i,j,:)     = mon_b(mon)%wcons(i,j,:)     + lnd(i,j)%water_cons    * mon_avg
            endif
          enddo
        enddo
        mon_b(mon)%econs_so  = mon_b(mon)%econs_so  + lnd%energy_cons_soil      * mon_avg
      endif

      ! monthly values 
      if( time_eom_lnd ) then

        if( write_surf ) then
          mon_su(mon)%f_land = lnd%f_land
          mon_su(mon)%f_veg = lnd%f_veg
          ann_su%f_land = lnd%f_land
          ann_su%f_veg = lnd%f_veg
          ann_su%alt= lnd%alt
          ann_su%gdd5= lnd%gdd5
          ann_su%t2m_min_mon= lnd%t2m_min_mon - T0  ! degC
          do i=1,nx
            do j=1,ny
              sum_frac = sum(lnd(i,j)%frac_surf)
              if (sum_frac.gt.0._wp) then
                mon_su(mon)%frac_surf(i,j,:) = lnd(i,j)%frac_surf/sum_frac
              else
                mon_su(mon)%frac_surf(i,j,:) = 0._wp
              endif
            enddo
          enddo
        endif

        if( write_soil_par ) then
          do i=1,nx
            do j=1,ny
              mon_sp(mon)%fsoc(i,j,:)     = lnd(i,j)%frac_soc
              mon_sp(mon)%psi(i,j,:)      = lnd(i,j)%psi
              mon_sp(mon)%ftemp(i,j,:)    = lnd(i,j)%ftemp
              mon_sp(mon)%fmoist(i,j,:)   = lnd(i,j)%fmoist
              mon_sp(mon)%fdepth(i,j,:)   = lnd(i,j)%fdepth
            enddo
          enddo
        endif

        if( write_carbon ) then

          where (mon_c(mon)%gpp.gt.0._wp)
            mon_c(mon)%disc = mon_c(mon)%disc / mon_c(mon)%gpp
          elsewhere
            mon_c(mon)%disc = 0._wp
          endwhere

          do i=1,nx
            do j=1,ny
              mon_c(mon)%sresp(i,j,:)    = lnd(i,j)%soil_resp*sec_year  ! kgC/m2/yr
              mon_c(mon)%sresp13(i,j,:)  = lnd(i,j)%soil_resp13*sec_year  ! kgC/m2/yr
              mon_c(mon)%sresp14(i,j,:)  = lnd(i,j)%soil_resp14*sec_year  ! kgC/m2/yr
              do n=1,ncarb
                mon_c(mon)%litter(i,j,n) = sum(lnd(i,j)%litterfall(:,n))*sec_year ! kgC/m2/yr
                mon_c(mon)%litter_prof(i,j,:,n) = lnd(i,j)%litterfall(:,n)*sec_year ! kgC/m2/yr
                mon_c(mon)%litter13_prof(i,j,:,n) = lnd(i,j)%litterfall13(:,n)*sec_year ! kgC/m2/yr
                mon_c(mon)%litter14_prof(i,j,:,n) = lnd(i,j)%litterfall14(:,n)*sec_year ! kgC/m2/yr
              enddo
            enddo
          enddo
          mon_c(mon)%ch4    = (lnd%ch4_emis_wetland*max(0._wp,lnd%f_wetland*lnd%f_veg-lnd%f_peat) &
          + lnd%ch4_emis_shelf*lnd%f_shelf &
          + lnd%ch4_emis_lake*lnd%f_lake &
          + lnd%ch4_emis_peat*lnd%f_peat) &
          * sec_year * 1.d3 ! gC/m2 gridcell/yr
          mon_c(mon)%ch4wet = lnd%ch4_emis_wetland*max(0._wp,lnd%f_wetland*lnd%f_veg-lnd%f_peat)*sec_year * 1.d3 ! gC/m2 gridcell/yr
          mon_c(mon)%ch4shelf = lnd%ch4_emis_shelf*lnd%f_shelf*sec_year * 1.d3 ! gC/m2 gridcell/yr
          mon_c(mon)%ch4lake  = lnd%ch4_emis_lake*lnd%f_lake*sec_year * 1.d3 ! gC/m2 gridcell/yr
          mon_c(mon)%ch4peat= lnd%ch4_emis_peat*lnd%f_peat*sec_year * 1.d3 ! gC/m2 gridcell/yr

          ! grid cell averages
          do i=1,nx
            do j=1,ny
              frac(:) = lnd(i,j)%frac_surf(1:npft)
              mon_c_g(mon)%lai(i,j)   = sum(mon_c(mon)%lai(i,j,1:npft)*frac)
              mon_c_g(mon)%gcan(i,j)  = sum(mon_c(mon)%gcan(i,j,1:npft)*frac)
              mon_c_g(mon)%gpp(i,j)   = sum(mon_c(mon)%gpp(i,j,1:npft)*frac)
              mon_c_g(mon)%aresp(i,j) = sum(mon_c(mon)%aresp(i,j,1:npft)*frac)
              if (lnd(i,j)%f_veg.gt.0._wp) then
                mon_c_g(mon)%disc(i,j)  = sum(mon_c(mon)%disc(i,j,1:npft)*frac,mask=mon_c(mon)%gpp(i,j,1:npft)>0._wp)/lnd(i,j)%f_veg
              else
                mon_c_g(mon)%disc(i,j)  = 0._wp
              endif
            enddo
          enddo

          mon_c_g(mon)%npp   = lnd%npp_real*sec_year  ! kgC/m2/yr
          mon_c_g(mon)%sresp = mon_c(mon)%sresp(:,:,ic_min)*(lnd%f_veg-lnd%f_peat) &
          + mon_c(mon)%sresp(:,:,ic_peat)*lnd%f_peat &
          + mon_c(mon)%sresp(:,:,ic_shelf)*lnd%f_shelf &
          + mon_c(mon)%sresp(:,:,ic_lake)*lnd%f_lake &
          + mon_c(mon)%sresp(:,:,ic_ice)*lnd%f_ice_grd ! kgC/m2/yr
          mon_c_g(mon)%npp13 = lnd%npp13_real*sec_year  ! kgC/m2/yr
          mon_c_g(mon)%sresp13 = mon_c(mon)%sresp13(:,:,ic_min)*(lnd%f_veg-lnd%f_peat) &
          + mon_c(mon)%sresp13(:,:,ic_peat)*lnd%f_peat &
          + mon_c(mon)%sresp13(:,:,ic_shelf)*lnd%f_shelf &
          + mon_c(mon)%sresp13(:,:,ic_lake)*lnd%f_lake &
          + mon_c(mon)%sresp13(:,:,ic_ice)*lnd%f_ice_grd ! kgC/m2/yr
          mon_c_g(mon)%npp14 = lnd%npp14_real*sec_year  ! kgC/m2/yr
          mon_c_g(mon)%sresp14 = mon_c(mon)%sresp14(:,:,ic_min)*(lnd%f_veg-lnd%f_peat) &
          + mon_c(mon)%sresp14(:,:,ic_peat)*lnd%f_peat &
          + mon_c(mon)%sresp14(:,:,ic_shelf)*lnd%f_shelf &
          + mon_c(mon)%sresp14(:,:,ic_lake)*lnd%f_lake &
          + mon_c(mon)%sresp14(:,:,ic_ice)*lnd%f_ice_grd ! kgC/m2/yr
          mon_c_g(mon)%Cflx_atm_lnd  = mon_c_g(mon)%npp*lnd%f_veg - mon_c_g(mon)%sresp ! kgC/m2/yr
          mon_c_g(mon)%C13flx_atm_lnd = mon_c_g(mon)%npp13*lnd%f_veg - mon_c_g(mon)%sresp13 ! kgC/m2/yr 
          mon_c_g(mon)%C14flx_atm_lnd = mon_c_g(mon)%npp14*lnd%f_veg - mon_c_g(mon)%sresp14 ! kgC/m2/yr

          ! at end of year
          if( time_eoy_lnd ) then

            do i=1,nx
              do j=1,ny
                ann_c%veg_h(i,j,:) = lnd(i,j)%veg_h
                ann_c%pfts(i,j,:)  = lnd(i,j)%pft_frac
                ann_c%seeds(i,j,:)  = lnd(i,j)%seed_frac
                ann_c%sai(i,j,:)   = lnd(i,j)%sai
                ann_c%lambda(i,j,:)   = lnd(i,j)%lambda
                ann_c%gamma_dist(i,j,:)  = 1._wp/max(1.e-4_wp,lnd(i,j)%gamma_dist*sec_year)    ! years
                ann_c%gamma_luc(i,j,:)   = 1._wp/max(1.e-4_wp,lnd(i,j)%gamma_luc*sec_year)     ! years
                ann_c%gamma_ice(i,j,:)   = 1._wp/max(1.e-4_wp,lnd(i,j)%gamma_ice*sec_year)     ! years
                ann_c%vegc(i,j,:)   = lnd(i,j)%veg_c
                ann_c%soilc(i,j,:)  = lnd(i,j)%soil_c_tot
              enddo
            enddo
            ann_c%bare   = 1._wp - sum(ann_c%pfts,3)
            ann_c%fpeat= lnd%f_peat
            ann_c%fpeatpot= lnd%f_peat_pot
            ann_c%dCpeatdt= lnd%dCpeat_dt*sec_year*1.d3 ! gC/m2/yr
            ann_c%acroh= lnd%acro_h
            ann_c%catoh= lnd%cato_h
            do i=1,nx
              do j=1,ny
                ann_c%litterc_prof(i,j,:,ic_min)   = lnd(i,j)%litter_c
                ann_c%fastc_prof(i,j,:,ic_min)  = lnd(i,j)%fast_c
                ann_c%slowc_prof(i,j,:,ic_min)  = lnd(i,j)%slow_c
                ann_c%slowc_prof(i,j,:,ic_peat)   = lnd(i,j)%cato_c
                ann_c%litterc_prof(i,j,:,ic_shelf)   = lnd(i,j)%litter_c_shelf
                ann_c%fastc_prof(i,j,:,ic_shelf)  = lnd(i,j)%fast_c_shelf
                ann_c%slowc_prof(i,j,:,ic_shelf)  = lnd(i,j)%slow_c_shelf
                ann_c%litterc_prof(i,j,:,ic_lake)   = lnd(i,j)%litter_c_lake
                ann_c%fastc_prof(i,j,:,ic_lake)  = lnd(i,j)%fast_c_lake
                ann_c%slowc_prof(i,j,:,ic_lake)  = lnd(i,j)%slow_c_lake
                ann_c%litterc_prof(i,j,:,ic_ice)   = lnd(i,j)%litter_c_ice
                ann_c%fastc_prof(i,j,:,ic_ice)  = lnd(i,j)%fast_c_ice
                ann_c%slowc_prof(i,j,:,ic_ice)  = lnd(i,j)%slow_c_ice
                ann_c%litterc13_prof(i,j,:,ic_min)   = lnd(i,j)%litter_c13
                ann_c%fastc13_prof(i,j,:,ic_min)  = lnd(i,j)%fast_c13
                ann_c%slowc13_prof(i,j,:,ic_min)  = lnd(i,j)%slow_c13
                ann_c%slowc13_prof(i,j,:,ic_peat)   = lnd(i,j)%cato_c13
                ann_c%litterc13_prof(i,j,:,ic_shelf)   = lnd(i,j)%litter_c13_shelf
                ann_c%fastc13_prof(i,j,:,ic_shelf)  = lnd(i,j)%fast_c13_shelf
                ann_c%slowc13_prof(i,j,:,ic_shelf)  = lnd(i,j)%slow_c13_shelf
                ann_c%litterc13_prof(i,j,:,ic_lake)   = lnd(i,j)%litter_c13_lake
                ann_c%fastc13_prof(i,j,:,ic_lake)  = lnd(i,j)%fast_c13_lake
                ann_c%slowc13_prof(i,j,:,ic_lake)  = lnd(i,j)%slow_c13_lake
                ann_c%litterc13_prof(i,j,:,ic_ice)   = lnd(i,j)%litter_c13_ice
                ann_c%fastc13_prof(i,j,:,ic_ice)  = lnd(i,j)%fast_c13_ice
                ann_c%slowc13_prof(i,j,:,ic_ice)  = lnd(i,j)%slow_c13_ice
                ann_c%litterc14_prof(i,j,:,ic_min)   = lnd(i,j)%litter_c14
                ann_c%fastc14_prof(i,j,:,ic_min)  = lnd(i,j)%fast_c14
                ann_c%slowc14_prof(i,j,:,ic_min)  = lnd(i,j)%slow_c14
                ann_c%slowc14_prof(i,j,:,ic_peat)   = lnd(i,j)%cato_c14
                ann_c%litterc14_prof(i,j,:,ic_shelf)   = lnd(i,j)%litter_c14_shelf
                ann_c%fastc14_prof(i,j,:,ic_shelf)  = lnd(i,j)%fast_c14_shelf
                ann_c%slowc14_prof(i,j,:,ic_shelf)  = lnd(i,j)%slow_c14_shelf
                ann_c%litterc14_prof(i,j,:,ic_lake)   = lnd(i,j)%litter_c14_lake
                ann_c%fastc14_prof(i,j,:,ic_lake)  = lnd(i,j)%fast_c14_lake
                ann_c%slowc14_prof(i,j,:,ic_lake)  = lnd(i,j)%slow_c14_lake
                ann_c%litterc14_prof(i,j,:,ic_ice)   = lnd(i,j)%litter_c14_ice
                ann_c%fastc14_prof(i,j,:,ic_ice)  = lnd(i,j)%fast_c14_ice
                ann_c%slowc14_prof(i,j,:,ic_ice)  = lnd(i,j)%slow_c14_ice
                do n=1,ncarb
                  do k=1,nlc
                    ann_c%litter_prof(i,j,k,n)   = 0._wp
                    ann_c%litter13_prof(i,j,k,n) = 0._wp 
                    ann_c%litter14_prof(i,j,k,n) = 0._wp
                    do m=1,nmon_year
                      ann_c%litter_prof(i,j,k,n)   = ann_c%litter_prof(i,j,k,n)   + mon_c(m)%litter_prof(i,j,k,n)/nmon_year  ! kgC/m2/yr 
                      ann_c%litter13_prof(i,j,k,n) = ann_c%litter13_prof(i,j,k,n) + mon_c(m)%litter13_prof(i,j,k,n)/nmon_year
                      ann_c%litter14_prof(i,j,k,n) = ann_c%litter14_prof(i,j,k,n) + mon_c(m)%litter14_prof(i,j,k,n)/nmon_year
                    enddo
                  enddo
                enddo
              enddo
            enddo
            ann_c%litterc_prof(:,:,:,ic_peat) = 0._wp
            ann_c%fastc_prof(:,:,:,ic_peat) = 0._wp
            ann_c%litterc_prof(:,:,1,ic_peat)  = lnd%litter_c_peat*rdz(1) ! kgC/m3
            ann_c%fastc_prof(:,:,1,ic_peat) = lnd%acro_c*rdz(1) ! kgC/m3
            ann_c%soilc_prof = ann_c%litterc_prof + ann_c%fastc_prof + ann_c%slowc_prof
            ann_c%litterc13_prof(:,:,:,ic_peat) = 0._wp
            ann_c%fastc13_prof(:,:,:,ic_peat) = 0._wp
            ann_c%litterc13_prof(:,:,1,ic_peat)  = lnd%litter_c13_peat*rdz(1) ! kgC/m3
            ann_c%fastc13_prof(:,:,1,ic_peat) = lnd%acro_c13*rdz(1) ! kgC/m3
            ann_c%soilc13_prof = ann_c%litterc13_prof + ann_c%fastc13_prof + ann_c%slowc13_prof
            ann_c%litterc14_prof(:,:,:,ic_peat) = 0._wp
            ann_c%fastc14_prof(:,:,:,ic_peat) = 0._wp
            ann_c%litterc14_prof(:,:,1,ic_peat)  = lnd%litter_c14_peat*rdz(1) ! kgC/m3
            ann_c%fastc14_prof(:,:,1,ic_peat) = lnd%acro_c14*rdz(1) ! kgC/m3
            ann_c%soilc14_prof = ann_c%litterc14_prof + ann_c%fastc14_prof + ann_c%slowc14_prof

            where (ann_c%litter_prof.gt.epsilon(1._wp))
              ann_c%d13c_litter_prof = (ann_c%litter13_prof/ann_c%litter_prof / c13_c12_std - 1._wp ) * 1000._wp
            elsewhere
              ann_c%d13c_litter_prof = missing_value 
            endwhere
            where (ann_c%litter_prof.gt.epsilon(1._wp) .and. ann_c%litter13_prof.gt.epsilon(1._wp))
              ann_c%D14c_litter_prof = (ann_c%litter14_prof/ann_c%litter_prof*(0.975_wp/(1._wp+ann_c%d13c_litter_prof/1000._wp))**2 / c14_c_std &
              - 1._wp) * 1000._wp
            elsewhere
              ann_c%D14c_litter_prof = missing_value
            endwhere

            where (ann_c%soilc_prof.gt.epsilon(1._wp))
              ann_c%d13c_prof = (ann_c%soilc13_prof/ann_c%soilc_prof / c13_c12_std - 1._wp ) * 1000._wp
            elsewhere
              ann_c%d13c_prof = missing_value 
            endwhere
            where (ann_c%soilc_prof.gt.epsilon(1._wp) .and. ann_c%soilc13_prof.gt.epsilon(1._wp))
              ann_c%D14c_prof = (ann_c%soilc14_prof/ann_c%soilc_prof*(0.975_wp/(1._wp+ann_c%d13c_prof/1000._wp))**2 / c14_c_std &
              - 1._wp) * 1000._wp
            elsewhere
              ann_c%D14c_prof = missing_value 
            endwhere

            do i=1,nx
              do j=1,ny
                where( lnd(i,j)%veg_c.gt.epsilon(1._wp) ) 
                  ann_c%d13c_veg(i,j,:) = (lnd(i,j)%veg_c13/lnd(i,j)%veg_c / c13_c12_std - 1._wp ) * 1000._wp ! permil
                elsewhere
                  ann_c%d13c_veg(i,j,:) = missing_value 
                endwhere
                where( lnd(i,j)%veg_c.gt.epsilon(1._wp) .and. lnd(i,j)%veg_c13.gt.epsilon(1._wp) ) 
                  ann_c%D14c_veg(i,j,:) = (lnd(i,j)%veg_c14/lnd(i,j)%veg_c*(0.975_wp/(1._wp+ann_c%d13c_veg(i,j,:)/1000._wp))**2 / c14_c_std &
                    - 1._wp) * 1000._wp
                elsewhere
                  ann_c%D14c_veg(i,j,:) = missing_value
                endwhere
                where( lnd(i,j)%soil_c_tot.gt.epsilon(1._wp) ) 
                  ann_c%d13c_soil(i,j,:) = (lnd(i,j)%soil_c13_tot/lnd(i,j)%soil_c_tot / c13_c12_std - 1._wp ) * 1000._wp 
                elsewhere
                  ann_c%d13c_soil(i,j,:) = missing_value 
                endwhere
                ! carbon weighted radiocarbon
                do n=1,ncarb
                  if (lnd(i,j)%soil_c_tot(n).gt.epsilon(1._wp) .and. lnd(i,j)%soil_c13_tot(n).gt.epsilon(1._wp)) then
                    ann_c%D14c_soil(i,j,n) = 0._wp
                    do k=1,nl
                      if (ann_c%soilc_prof(i,j,k,n).gt.0._wp) then
                        ann_c%D14c_soil(i,j,n) = ann_c%D14c_soil(i,j,n) &
                          + (ann_c%soilc14_prof(i,j,k,n)/ann_c%soilc_prof(i,j,k,n)*(0.975_wp/(1._wp+ann_c%d13c_soil(i,j,n)/1000._wp))**2/c14_c_std &
                          - 1._wp) * 1000._wp * ann_c%soilc_prof(i,j,k,n)*dz(k)/lnd(i,j)%soil_c_tot(n)
                      endif
                    enddo
                  else
                    ann_c%D14c_soil(i,j,n) = missing_value 
                  endif
                enddo
              enddo
            enddo

            ann_c%doc_export = lnd%doc_export
            ann_c%poc_export = lnd%poc_export
            ann_c%weath_carb = lnd%weath_carb
            ann_c%weath_sil  = lnd%weath_sil
            do i=1,nx
              do j=1,ny
                if (i_weathering.eq.1) then
                  ann_c%lithology(i,j,:)  = lnd(i,j)%lithology_gemco2(:)
                else if (i_weathering.eq.2) then
                  ann_c%lithology(i,j,:)  = lnd(i,j)%lithology_uhh(:)
                endif
              enddo
            enddo

            ! grid cell averages
            do i=1,nx
              do j=1,ny
                frac(:) = lnd(i,j)%frac_surf(1:npft)
                ann_c_g%vegc(i,j)  = sum(ann_c%vegc(i,j,1:npft)*frac)
                if (lnd(i,j)%f_veg.gt.0._wp) then 
                  ann_c_g%d13c_veg(i,j) = sum(ann_c%d13c_veg(i,j,1:npft)*frac)/lnd(i,j)%f_veg
                  ann_c_g%D14c_veg(i,j) = sum(ann_c%D14c_veg(i,j,1:npft)*frac)/lnd(i,j)%f_veg
                else
                  ann_c_g%d13c_veg(i,j) = missing_value
                  ann_c_g%D14c_veg(i,j) = missing_value 
                endif
                ann_c_g%soilc(i,j) = lnd(i,j)%soil_c_tot(ic_min) * (lnd(i,j)%f_veg-lnd(i,j)%f_peat) &
                  + lnd(i,j)%soil_c_tot(ic_peat) * lnd(i,j)%f_peat &
                  + lnd(i,j)%soil_c_tot(ic_shelf) * lnd(i,j)%f_shelf &
                  + lnd(i,j)%soil_c_tot(ic_lake) * lnd(i,j)%f_lake &
                  + lnd(i,j)%soil_c_tot(ic_ice) * lnd(i,j)%f_ice_grd
              enddo
            enddo

          endif

        endif

        if( write_cons ) then
          do i=1,nx
            do j=1,ny
              mon_b(mon)%ccons_s(i,j,:) = lnd(i,j)%carbon_cons_soil
              mon_b(mon)%ccons13_s(i,j,:) = lnd(i,j)%carbon13_cons_soil
              mon_b(mon)%ccons14_s(i,j,:) = lnd(i,j)%carbon14_cons_soil
            enddo
          enddo
          mon_b(mon)%ccons_v = lnd%carbon_cons_veg
          mon_b(mon)%ccons13_v = lnd%carbon13_cons_veg
          mon_b(mon)%ccons14_v = lnd%carbon14_cons_veg
        endif

      endif

      if (time_out_lnd .and. time_eoy_lnd) then
        call lnd_diag_out
      endif


    endif
    
  end subroutine lnd_diag


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  l n d _ d i a g _ o u t
  ! Purpose  :  Initialize netcdf output for paladyn
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lnd_diag_out

    implicit none

    integer :: k, ncid
    character (len=256) :: fnm


    nout = nout + 1

!    !$omp parallel sections

!    !$omp section
    ! surface
    if( write_surf ) then
      call surf_ave( mon_su,ann_su,mon_su_g,ann_su_g )
      fnm = trim(out_dir)//"/lnd_surf.nc"
      call nc_open(fnm,ncid)
      call nc_write(fnm,dim_time,real(year_now,wp),dim1=dim_time,start=[nout],count=[1],ncid=ncid)
      do k = 1, nmon_year
        call surf_nc_write(fnm,ncid,mon_su(k),mon_su_g(k),k,nout)
      end do
      call surf_nc_write(fnm,ncid,ann_su,ann_su_g,nmon_year+1,nout)
      call nc_close(ncid)
    endif

    if ( write_surf .and. l_daily_output) then
      ! write to file
      fnm = trim(out_dir)//"/lnd_surf_daily.nc"
      call nc_open(fnm,ncid)
      call nc_write(fnm,dim_time,real(year_now,wp), dim1=dim_time, start=[nout], count=[1],ncid=ncid)    
      do k = 1, nday_year
        call surf_daily_nc_write(fnm,ncid,day_su(k),k,nout)
      end do
      call nc_close(ncid)
    endif

!    !$omp section
    ! carbon
    if( write_carbon ) then
      call carbon_ave( mon_c,ann_c,mon_c_g,ann_c_g )
      fnm = trim(out_dir)//"/lnd_carb.nc"
      call nc_open(fnm,ncid)
      call nc_write(fnm,dim_time,real(year_now,wp),dim1=dim_time,start=[nout],count=[1],ncid=ncid)
      do k = 1, nmon_year
        call carbon_nc_write(fnm,ncid,mon_c(k),mon_c_g(k),k,nout)
      end do
      call carbon_nc_write(fnm,ncid,ann_c,ann_c_g,nmon_year+1,nout)
      call nc_close(ncid)
    endif

!    !$omp section
    ! soil
    if( write_soil ) then
      call soil_ave( mon_s,ann_s )
      fnm = trim(out_dir)//"/lnd_soil.nc"
      call nc_open(fnm,ncid)
      call nc_write(fnm,dim_time,real(year_now,wp),dim1=dim_time,start=[nout],count=[1],ncid=ncid)
      do k = 1, nmon_year
        call soil_nc_write(fnm,ncid,mon_s(k),k,nout)
      end do
      call soil_nc_write(fnm,ncid,ann_s,nmon_year+1,nout)
      call nc_close(ncid)
    endif

!    !$omp section
    ! soil parameters
    if( write_soil_par ) then
      call soil_par_ave( mon_sp,ann_sp )
      fnm = trim(out_dir)//"/lnd_soil_par.nc"
      call nc_open(fnm,ncid)
      call nc_write(fnm,dim_time,real(year_now,wp),dim1=dim_time,start=[nout],count=[1],ncid=ncid)
      do k = 1, nmon_year
        call soil_par_nc_write(fnm,ncid,mon_sp(k),k,nout)
      end do
      call soil_par_nc_write(fnm,ncid,ann_sp,nmon_year+1,nout)
      call nc_close(ncid)
    endif

    ! lake
    if( write_lake ) then
      call lake_ave( mon_l,ann_l )
      fnm = trim(out_dir)//"/lnd_lake.nc"
      call nc_open(fnm,ncid)
      call nc_write(fnm,dim_time,real(year_now,wp),dim1=dim_time,start=[nout],count=[1],ncid=ncid)
      do k = 1, nmon_year
        call lake_nc_write(fnm,ncid,mon_l(k),k,nout)
      end do
      call lake_nc_write(fnm,ncid,ann_l,nmon_year+1,nout)
      call nc_close(ncid)
    endif

!    !$omp section
    ! balance
    if( write_cons ) then
      call cons_ave( mon_b,ann_b )
      fnm = trim(out_dir)//"/lnd_conservation.nc"
      call nc_open(fnm,ncid)
      call nc_write(fnm,dim_time,real(year_now,wp),dim1=dim_time,start=[nout],count=[1],ncid=ncid)
      do k = 1, nmon_year
        call cons_nc_write(fnm,ncid,mon_b(k),k,nout)
      end do
      call cons_nc_write(fnm,ncid,ann_b,nmon_year+1,nout)
      call nc_close(ncid)
    endif

!    !$omp end parallel sections


   end subroutine lnd_diag_out


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t s _ n c
  ! Purpose  :  Initialize netcdf output for paladyn
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ts_nc(fnm)

    implicit none

    character (len=*) :: fnm
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_write_dim(fnm,dim_time, x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE.)
    !call nc_write_dim(fnm,dim_time, x=[year_ini:year_ini+nyears], axis="t", units="years BP", unlimited=.TRUE.)
    call nc_write_dim(fnm, dim_npft, x=i_pft, axis="z", units="n/a")
    call nc_write_dim(fnm, dim_lat, x=1, axis="y", units="1")
    call nc_write_dim(fnm, dim_lon, x=1, axis="x", units="1")

    return

  end subroutine ts_nc


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t s _ n c _ w r i t e
  ! Purpose  :  Output timestep of netcdf output for paladyn
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ts_nc_write(fnm,vars,ndat,y)

    implicit none

    type(ts_out) :: vars(:)

    character (len=*) :: fnm
    integer :: ndat, y, ncid, n, i

    call nc_open(fnm,ncid)
    call nc_write(fnm,"time", real([(i,i=(year_now-(y-1)*n_accel),(year_now),(n_accel))],wp), &
    dim1=dim_time,start=[ndat],count=[y],ncid=ncid)   
    call nc_write(fnm,"co2", sngl(vars%co2),    dim1=dim_time,start=[ndat],count=[y],long_name="atmospheric CO2",units="ppm",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"temp",    sngl(vars%t-T0),   dim1=dim_time,start=[ndat],count=[y],long_name="land temperaure",units="degC",missing_value=missing_value,ncid=ncid)

    call nc_write(fnm,"Aland",    sngl(vars%land),   dim1=dim_time,start=[ndat],count=[y],long_name="global land area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Aveg",     sngl(vars%veg),    dim1=dim_time,start=[ndat],count=[y],long_name="global vegetated area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Aice",     sngl(vars%ice),    dim1=dim_time,start=[ndat],count=[y],long_name="global ice area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Ashelf",   sngl(vars%shelf),  dim1=dim_time,start=[ndat],count=[y],long_name="global shelf area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Alake",    sngl(vars%lake),   dim1=dim_time,start=[ndat],count=[y],long_name="global lake area",units="mln km^2",missing_value=missing_value,ncid=ncid)

    call nc_write(fnm,"Cland",   sngl(vars%landc),  dim1=dim_time,start=[ndat],count=[y],long_name="total land carbon",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cflx_atm_lnd",    sngl(vars%Cflx_atm_lnd),   dim1=dim_time,start=[ndat],count=[y],long_name="net atmosphere-land carbon flux",units="PgC/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"C13flx_atm_lnd",    sngl(vars%C13flx_atm_lnd),   dim1=dim_time,start=[ndat],count=[y],long_name="net atmosphere-land carbon 13 flux",units="PgC/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"C14flx_atm_lnd",    sngl(vars%C14flx_atm_lnd),   dim1=dim_time,start=[ndat],count=[y],long_name="net atmosphere-land carbon 14 flux",units="PgC/yr",missing_value=missing_value,ncid=ncid)
!    call nc_write(fnm,"d13Cflx_atm_lnd",    sngl(vars%d13Cflx_atm_lnd),   dim1=dim_time,start=[ndat],count=[y],long_name="d13C of net atmosphere-land carbon flux",units="permil",missing_value=missing_value,ncid=ncid)
!    call nc_write(fnm,"D14Cflx_atm_lnd",    sngl(vars%D14Cflx_atm_lnd),   dim1=dim_time,start=[ndat],count=[y],long_name="D14C of net atmosphere-land carbon flux",units="permil",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cflx_burial", sngl(vars%Cflx_burial),  dim1=dim_time,start=[ndat],count=[y],long_name="carbon burial flux (lost from the system)",units="PgC/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"d13Cflx_burial", sngl(vars%d13Cflx_burial),  dim1=dim_time,start=[ndat],count=[y],long_name="d13C of carbon burial flux",units="permil",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"gpp",     sngl(vars%gpp),    dim1=dim_time,start=[ndat],count=[y],long_name="gross primary productivity",units="PgC/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"npp",     sngl(vars%npp),    dim1=dim_time,start=[ndat],count=[y],long_name="net primary productivity",units="PgC/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"npp_pimask",     sngl(vars%npp_pimask),    dim1=dim_time,start=[ndat],count=[y],long_name="net primary productivity in pre-industrial permafrost area",units="PgC/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"sresp",   sngl(vars%sresp),  dim1=dim_time,start=[ndat],count=[y],long_name="soil respiration",units="PgC/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"runoff",  sngl(vars%runoff), dim1=dim_time,start=[ndat],count=[y],long_name="total runoff (surface runoff + drainage)",units="10^15 kg/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"runsur",  sngl(vars%runsur), dim1=dim_time,start=[ndat],count=[y],long_name="surface runoff",units="10^15 kg/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"calving", sngl(vars%calving),dim1=dim_time,start=[ndat],count=[y],long_name="calving",units="10^15 kg/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"drain",   sngl(vars%drain),  dim1=dim_time,start=[ndat],count=[y],long_name="drainage",units="10^15 kg/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"evp",     sngl(vars%evp),    dim1=dim_time,start=[ndat],count=[y],long_name="evapotranspiration",units="10^15 kg/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"esur",    sngl(vars%esur),   dim1=dim_time,start=[ndat],count=[y],long_name="surface evaporation",units="10^15 kg/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"trans",   sngl(vars%trans),  dim1=dim_time,start=[ndat],count=[y],long_name="transpiration",units="10^15 kg/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"prc",     sngl(vars%prc),    dim1=dim_time,start=[ndat],count=[y],long_name="precipitation",units="10^15 kg/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Awet",     sngl(vars%wet),    dim1=dim_time,start=[ndat],count=[y],long_name="wetland area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Awettrop",     sngl(vars%wettrop),    dim1=dim_time,start=[ndat],count=[y],long_name="tropical wetland area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Awetextrop",     sngl(vars%wetextrop),    dim1=dim_time,start=[ndat],count=[y],long_name="extratropical wetland area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Ainund",   sngl(vars%inund),  dim1=dim_time,start=[ndat],count=[y],long_name="inundated area excluding peatlands",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Vsnow",    sngl(vars%snow),   dim1=dim_time,start=[ndat],count=[y],long_name="total snow volume",units="10^12 m^3",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Aforest",  sngl(vars%forest), dim1=dim_time,start=[ndat],count=[y],long_name="forest area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Agrass",   sngl(vars%grass),  dim1=dim_time,start=[ndat],count=[y],long_name="grass area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Ashrub",   sngl(vars%shrub),  dim1=dim_time,start=[ndat],count=[y],long_name="shrub area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Adesert",  sngl(vars%desert), dim1=dim_time,start=[ndat],count=[y],long_name="desert area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Acrop",    sngl(vars%crop),   dim1=dim_time,start=[ndat],count=[y],long_name="cropland area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Apasture", sngl(vars%pasture),dim1=dim_time,start=[ndat],count=[y],long_name="pasture area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"A_BL",    sngl(vars%pfts(1)), dim1=dim_time,start=[ndat],count=[y],long_name="Broadleaf forest area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"A_NL",    sngl(vars%pfts(2)), dim1=dim_time,start=[ndat],count=[y],long_name="Needleleaf forest area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"A_C3",    sngl(vars%pfts(3)), dim1=dim_time,start=[ndat],count=[y],long_name="C3 grassland area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"A_C4",    sngl(vars%pfts(4)), dim1=dim_time,start=[ndat],count=[y],long_name="C4 grassland area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"A_SH",    sngl(vars%pfts(5)), dim1=dim_time,start=[ndat],count=[y],long_name="Shrubland area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cveg",    sngl(vars%vegc),   dim1=dim_time,start=[ndat],count=[y],long_name="vegetation carbon",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cveg_pimask",    sngl(vars%vegc_pimask),   dim1=dim_time,start=[ndat],count=[y],long_name="vegetation carbon in pre-industrial permafrost area",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Csoil",   sngl(vars%soilc),  dim1=dim_time,start=[ndat],count=[y],long_name="total soil carbon",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Csoil60N",   sngl(vars%soilc60N),  dim1=dim_time,start=[ndat],count=[y],long_name="total soil carbon north of 60N",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Clitter", sngl(vars%litterc),   dim1=dim_time,start=[ndat],count=[y],long_name="litter carbon",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cfast",   sngl(vars%fastc),   dim1=dim_time,start=[ndat],count=[y],long_name="fast soil carbon",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cslow",   sngl(vars%slowc),   dim1=dim_time,start=[ndat],count=[y],long_name="slow soil carbon",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cmin",    sngl(vars%minc),   dim1=dim_time,start=[ndat],count=[y],long_name="mineral soil carbon",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cpeat",   sngl(vars%peatc),  dim1=dim_time,start=[ndat],count=[y],long_name="carbon in peat",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cice",    sngl(vars%icec),   dim1=dim_time,start=[ndat],count=[y],long_name="carbon below ice sheets",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cshelf",  sngl(vars%shelfc), dim1=dim_time,start=[ndat],count=[y],long_name="carbon below ocean shelf water",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Clake",   sngl(vars%lakec),  dim1=dim_time,start=[ndat],count=[y],long_name="carbon below lake water",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cperm",   sngl(vars%permc),  dim1=dim_time,start=[ndat],count=[y],long_name="carbon in permafrost area",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cperm_pimask",   sngl(vars%soilc_pimask),  dim1=dim_time,start=[ndat],count=[y],long_name="carbon in pre-industrial permafrost area",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cinert",   sngl(vars%inertc),  dim1=dim_time,start=[ndat],count=[y],long_name="inert land carbon (frozen soils and below ice sheets)",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Csoil1m",   sngl(vars%soilc1m),  dim1=dim_time,start=[ndat],count=[y],long_name="top meter soil carbon",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Csoil60N1m",   sngl(vars%soilc60N1m),  dim1=dim_time,start=[ndat],count=[y],long_name="top meter soil carbon north of 60N",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cmin1m",    sngl(vars%minc1m),   dim1=dim_time,start=[ndat],count=[y],long_name="top meter mineral soil carbon",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cpeat1m",   sngl(vars%peatc1m),  dim1=dim_time,start=[ndat],count=[y],long_name="top meter carbon in peat",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cice1m",    sngl(vars%icec1m),   dim1=dim_time,start=[ndat],count=[y],long_name="top meter carbon below ice sheets",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cshelf1m",  sngl(vars%shelfc1m), dim1=dim_time,start=[ndat],count=[y],long_name="top meter carbon below ocean shelf water",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Clake1m",   sngl(vars%lakec1m),  dim1=dim_time,start=[ndat],count=[y],long_name="top meter carbon below lake water",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cperm1m",   sngl(vars%permc1m),  dim1=dim_time,start=[ndat],count=[y],long_name="top meter carbon in permafrost area",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cperm1m_pimask",   sngl(vars%soilc1m_pimask),  dim1=dim_time,start=[ndat],count=[y],long_name="top meter carbon in pre-industrial permafrost area",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Csoil3m",   sngl(vars%soilc3m),  dim1=dim_time,start=[ndat],count=[y],long_name="top 3 meters soil carbon",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cmin3m",    sngl(vars%minc3m),   dim1=dim_time,start=[ndat],count=[y],long_name="top 3 meters mineral soil carbon",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cpeat3m",   sngl(vars%peatc3m),  dim1=dim_time,start=[ndat],count=[y],long_name="top 3 meters carbon in peat",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cice3m",    sngl(vars%icec3m),   dim1=dim_time,start=[ndat],count=[y],long_name="top 3 meters carbon below ice sheets",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cshelf3m",  sngl(vars%shelfc3m), dim1=dim_time,start=[ndat],count=[y],long_name="top 3 meters carbon below ocean shelf water",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Clake3m",   sngl(vars%lakec3m),  dim1=dim_time,start=[ndat],count=[y],long_name="top 3 meters carbon below lake water",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cperm3m",   sngl(vars%permc3m),  dim1=dim_time,start=[ndat],count=[y],long_name="top 3 meters carbon in permafrost area",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Cperm3m_pimask",   sngl(vars%soilc3m_pimask),  dim1=dim_time,start=[ndat],count=[y],long_name="top 3 meters carbon in pre-industrial permafrost area",units="PgC",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Apeat",    sngl(vars%peat),   dim1=dim_time,start=[ndat],count=[y],long_name="peatland area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Apeatpot", sngl(vars%peatpot),   dim1=dim_time,start=[ndat],count=[y],long_name="potential peatland area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Aperm",    sngl(vars%perm),   dim1=dim_time,start=[ndat],count=[y],long_name="permafrost area",units="mln km^2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ch4",     sngl(vars%ch4),    dim1=dim_time,start=[ndat],count=[y],long_name="methane emissions",units="TgCH4/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ch4trop",  sngl(vars%ch4trop), dim1=dim_time,start=[ndat],count=[y],long_name="tropical methane emissions",units="TgCH4/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ch4extrop", sngl(vars%ch4extrop),dim1=dim_time,start=[ndat],count=[y],long_name="extratropical methane emissions",units="TgCH4/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ch4shelf",sngl(vars%ch4shelf), dim1=dim_time,start=[ndat],count=[y],long_name="methane emissions from ocean shelf",units="TgCH4/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ch4lake",sngl(vars%ch4lake), dim1=dim_time,start=[ndat],count=[y],long_name="methane emissions from lakes",units="TgCH4/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"d13ch4",  sngl(vars%d13ch4), dim1=dim_time,start=[ndat],count=[y],long_name="d13C of methane emissions",units="permil",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"dust_e",  sngl(vars%dust_e), dim1=dim_time,start=[ndat],count=[y],long_name="total dust emission",units="Tg/yr",missing_value=missing_value,ncid=ncid)

    call nc_write(fnm,"doc_export",  sngl(vars%doc_export), dim1=dim_time,start=[ndat],count=[y],long_name="dissolved organic carbon export through rivers",units="PgC/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"poc_export",  sngl(vars%poc_export), dim1=dim_time,start=[ndat],count=[y],long_name="particulate organic carbon export through rivers",units="PgC/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"weath_carb",  sngl(vars%weath_carb), dim1=dim_time,start=[ndat],count=[y],long_name="carbonate weathering",units="Tmol C/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"weath_sil",  sngl(vars%weath_sil), dim1=dim_time,start=[ndat],count=[y],long_name="silicate weathering",units="Tmol C/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"weath_loess", sngl(vars%weath_loess), dim1=dim_time,start=[ndat],count=[y],long_name="loess carbonate weathering",units="Tmol C/yr",missing_value=missing_value,ncid=ncid)

    call nc_close(ncid)

   return

  end subroutine ts_nc_write


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t s _ a v e
  ! Purpose  :  Average (or sum) the paladyn fields as needed
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ts_ave(d,ave)
    
    implicit none
    
    type(ts_out) :: d(:), ave

    integer :: k, n
    real(wp) :: div
    
    n = size(d)
    div = real(n,wp)
   
    ! Set all values to zero
    ave%gpp     = 0._wp
    ave%npp     = 0._wp
    ave%npp_pimask     = 0._wp
    ave%sresp   = 0._wp
    ave%npp14   = 0._wp
    ave%sresp14 = 0._wp
    ave%runoff  = 0._wp
    ave%runsur  = 0._wp
    ave%calving = 0._wp
    ave%drain   = 0._wp
    ave%evp     = 0._wp
    ave%esur    = 0._wp
    ave%trans   = 0._wp
    ave%prc     = 0._wp
    ave%t       = 0._wp
    ave%wet     = 0._wp
    ave%wettrop     = 0._wp
    ave%wetextrop     = 0._wp
    ave%inund   = 0._wp
    ave%snow    = 0._wp
    ave%peat    = 0._wp
    ave%peatpot = 0._wp
    ave%forest  = 0._wp
    ave%desert  = 0._wp
    ave%grass   = 0._wp
    ave%shrub   = 0._wp
    ave%crop    = 0._wp
    ave%pasture = 0._wp
    ave%pfts    = 0._wp
    ave%vegc    = 0._wp
    ave%vegc_pimask    = 0._wp
    ave%soilc   = 0._wp
    ave%soilc60N   = 0._wp
    ave%litterc    = 0._wp
    ave%fastc    = 0._wp
    ave%slowc    = 0._wp
    ave%minc    = 0._wp
    ave%peatc   = 0._wp
    ave%permc   = 0._wp
    ave%soilc_pimask   = 0._wp
    ave%inertc  = 0._wp
    ave%icec    = 0._wp
    ave%shelfc  = 0._wp
    ave%lakec   = 0._wp
    ave%soilc1m   = 0._wp
    ave%soilc60N1m   = 0._wp
    ave%minc1m    = 0._wp
    ave%peatc1m   = 0._wp
    ave%permc1m   = 0._wp
    ave%soilc1m_pimask   = 0._wp
    ave%icec1m    = 0._wp
    ave%shelfc1m  = 0._wp
    ave%lakec1m   = 0._wp
    ave%soilc3m   = 0._wp
    ave%minc3m    = 0._wp
    ave%peatc3m   = 0._wp
    ave%permc3m   = 0._wp
    ave%soilc3m_pimask   = 0._wp
    ave%icec3m    = 0._wp
    ave%shelfc3m  = 0._wp
    ave%lakec3m   = 0._wp
    ave%ch4     = 0._wp
    ave%ch4trop  = 0._wp
    ave%ch4shelf= 0._wp
    ave%ch4lake = 0._wp
    ave%ch4extrop = 0._wp
    ave%d13ch4  = 0._wp
    ave%dust_e = 0._wp

    ! Loop over the time indices to sum up
    do k = 1, n
     ave%gpp     = ave%gpp        + d(k)%gpp
     ave%npp     = ave%npp        + d(k)%npp
     ave%npp_pimask     = ave%npp_pimask        + d(k)%npp_pimask
     ave%sresp   = ave%sresp      + d(k)%sresp
     ave%npp14   = ave%npp14      + d(k)%npp14
     ave%sresp14 = ave%sresp14    + d(k)%sresp14
     ave%runoff  = ave%runoff     + d(k)%runoff
     ave%runsur  = ave%runsur     + d(k)%runsur
     ave%calving = ave%calving    + d(k)%calving
     ave%drain   = ave%drain      + d(k)%drain
     ave%evp     = ave%evp        + d(k)%evp
     ave%esur    = ave%esur       + d(k)%esur
     ave%trans   = ave%trans      + d(k)%trans
     ave%prc     = ave%prc        + d(k)%prc
     ave%t       = ave%t          + d(k)%t      / div
     ave%wet     = ave%wet        + d(k)%wet    / div
     ave%wettrop     = ave%wettrop        + d(k)%wettrop    / div
     ave%wetextrop     = ave%wetextrop        + d(k)%wetextrop    / div
     ave%inund   = ave%inund      + d(k)%inund  / div
     ave%snow    = ave%snow       + d(k)%snow   / div
     ave%peat    = ave%peat       + d(k)%peat   / div
     ave%peatpot = ave%peatpot    + d(k)%peatpot   / div
     ave%forest  = ave%forest     + d(k)%forest / div
     ave%grass   = ave%grass      + d(k)%grass  / div
     ave%shrub   = ave%shrub      + d(k)%shrub  / div
     ave%desert  = ave%desert     + d(k)%desert / div
     ave%crop    = ave%crop       + d(k)%crop   / div
     ave%pasture = ave%pasture    + d(k)%pasture/ div
     ave%pfts    = ave%pfts       + d(k)%pfts   / div
     ave%vegc    = ave%vegc       + d(k)%vegc   / div
     ave%vegc_pimask    = ave%vegc_pimask       + d(k)%vegc_pimask   / div
     ave%soilc   = ave%soilc      + d(k)%soilc  / div
     ave%soilc60N   = ave%soilc60N      + d(k)%soilc60N  / div
     ave%litterc    = ave%litterc       + d(k)%litterc   / div
     ave%fastc    = ave%fastc       + d(k)%fastc   / div
     ave%slowc    = ave%slowc       + d(k)%slowc   / div
     ave%minc    = ave%minc       + d(k)%minc   / div
     ave%peatc   = ave%peatc      + d(k)%peatc  / div
     ave%permc   = ave%permc      + d(k)%permc  / div
     ave%soilc_pimask   = ave%soilc_pimask      + d(k)%soilc_pimask  / div
     ave%inertc  = ave%inertc     + d(k)%inertc / div
     ave%icec    = ave%icec       + d(k)%icec   / div
     ave%shelfc  = ave%shelfc     + d(k)%shelfc / div
     ave%lakec   = ave%lakec      + d(k)%lakec  / div
     ave%soilc1m   = ave%soilc1m      + d(k)%soilc1m  / div
     ave%soilc60N1m   = ave%soilc60N1m      + d(k)%soilc60N1m  / div
     ave%minc1m    = ave%minc1m       + d(k)%minc1m   / div
     ave%peatc1m   = ave%peatc1m      + d(k)%peatc1m  / div
     ave%permc1m   = ave%permc1m      + d(k)%permc1m  / div
     ave%soilc1m_pimask   = ave%soilc1m_pimask      + d(k)%soilc1m_pimask  / div
     ave%icec1m    = ave%icec1m       + d(k)%icec1m   / div
     ave%shelfc1m  = ave%shelfc1m     + d(k)%shelfc1m / div
     ave%lakec1m   = ave%lakec1m      + d(k)%lakec1m  / div
     ave%soilc3m   = ave%soilc3m      + d(k)%soilc3m  / div
     ave%minc3m    = ave%minc3m       + d(k)%minc3m   / div
     ave%peatc3m   = ave%peatc3m      + d(k)%peatc3m  / div
     ave%permc3m   = ave%permc3m      + d(k)%permc3m  / div
     ave%soilc3m_pimask   = ave%soilc3m_pimask      + d(k)%soilc3m_pimask  / div
     ave%icec3m    = ave%icec3m       + d(k)%icec3m   / div
     ave%shelfc3m  = ave%shelfc3m     + d(k)%shelfc3m / div
     ave%lakec3m   = ave%lakec3m      + d(k)%lakec3m  / div
     ave%ch4     = ave%ch4        + d(k)%ch4    / div
     ave%ch4trop = ave%ch4trop     + d(k)%ch4trop / div
     ave%ch4shelf= ave%ch4shelf   + d(k)%ch4shelf / div
     ave%ch4lake = ave%ch4lake    + d(k)%ch4lake  / div
     ave%ch4extrop = ave%ch4extrop    + d(k)%ch4extrop/ div
     ave%d13ch4  = ave%d13ch4     + d(k)%d13ch4/ div
     ave%dust_e     = ave%dust_e        + d(k)%dust_e
    end do
      
   return
    
  end subroutine ts_ave
  

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s u r f _ n c
  ! Purpose  :  Initialize netcdf output for paladyn
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine surf_nc(fnm)
    
    implicit none
    
    character (len=*) :: fnm
    integer :: ncid
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm, dim_time, x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE.,ncid=ncid)    
    call nc_write_dim(fnm, dim_month, x=1._wp, dx=1._wp, nx=13, axis="e", &
    units="months", ncid=ncid)
    call nc_write_dim(fnm, dim_nsurf, x=i_surf, units="n/a", axis="z", ncid=ncid)
    call nc_write_dim(fnm, dim_lat, x=lat, axis="y", units="degrees_north", ncid=ncid)
    call nc_write_dim(fnm, dim_lon, x=lon, axis="x", units="degrees_east", ncid=ncid)
    call nc_write_dim(fnm, dim_npft,x=i_pft,units="n/a" ,ncid=ncid)
    call nc_write_dim(fnm,dim_nsoil,x=i_soil,units="n/a", ncid=ncid)
    call nc_close(ncid)

    return
  
  end subroutine surf_nc
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s u r f _ n c _ w r i t e
  ! Purpose  :  Output timestep of netcdf output for paladyn
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine surf_nc_write(fnm,ncid,vars,vars_g,ndat,nout)
    
    implicit none
    
    type(surf_out) :: vars
    type(surf_g_out) :: vars_g   
 
    character (len=*) :: fnm
    integer :: ndat, nout, ncid

    
    call nc_write(fnm,"fland",    sngl(vars%f_land),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="land fraction",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fveg",     sngl(vars%f_veg),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="vegetation fraction",units="/",missing_value=missing_value,ncid=ncid)

    if( ndat .eq. nmon_year+1 ) then
      call nc_write(fnm,"fsurf",    sngl(vars%frac_surf),dims=[dim_lon,dim_lat,dim_nsurf,dim_time],start=[1,1,1,nout],count=[nx,ny,nsurf,1],long_name="surface type fractions",units="/",missing_value=missing_value,ncid=ncid)
    endif

    if (write_surf_n) then
      ! write output for each surface type

    call nc_write(fnm,"etot",     sngl(vars%etot), dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1], &
      long_name="evapotranspiration",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ecan",     sngl(vars%ecan), dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1], &
      long_name="evaporation from canopy",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"esur",     sngl(vars%esur), dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1], & 
      long_name="evaporation from soil",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"trans",    sngl(vars%trans),dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1], & 
      long_name="transpiration from vegetation",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"le",       sngl(vars%le),   dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1], & 
      long_name="latent heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"sh",       sngl(vars%sh),   dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1], & 
      long_name="sensible heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"slh",      sngl(vars%slh),   dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1], & 
      long_name="sum of sensible and latent heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ef",       sngl(vars%ef),   dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1], & 
      long_name="evaporative fraction",units="1",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"swnet",    sngl(vars%swnet),dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1], & 
      long_name="net shortwave radiation at the surface",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"lwnet",    sngl(vars%lwnet),dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1], & 
      long_name="net longwave radiation at the surface",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"g",        sngl(vars%g),    dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1], & 
      long_name="ground heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"flx_melt", sngl(vars%flx_melt),dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1], & 
      long_name="heat flux going into snowmelt",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"tskin",    sngl(vars%tskin),dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1], & 
      long_name="skin temperature",units="K",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"tskin_amp",sngl(vars%tskin_amp),dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="amplitude of diurnal cycle of skin temperature",units="K",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"alb",      sngl(vars%alb),  dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="surface albedo",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"alb_dir",  sngl(vars%alb_dir),dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="surface direct beam albedo",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"alb_dif",  sngl(vars%alb_dif),dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="surface diffuse albedo",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"scf",      sngl(vars%scf),dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="snow cover fraction",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fcanwat", sngl(vars%fcanwat),dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="fraction of canopy covered by water",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fcansn",   sngl(vars%fcansn),dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="fraction of canopy covered by snow",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"wind",     sngl(vars%wind),   dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="wind speed",units="m/s",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ra",       sngl(vars%ra),   dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="aerodynamic resistance",units="s/m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"rag",      sngl(vars%rag),  dims=[dim_lon,dim_lat,dim_npft,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,npft,1,1],long_name="aerodynamic resistance below canopy",units="s/m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Ri",     sngl(vars%Ri), dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="bulk richardson number",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"z0m",    sngl(vars%z0m), dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="surface roughness length for momentum",units="m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"z0h",    sngl(vars%z0h), dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="surface roughness length for scalars",units="m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Ch",     sngl(vars%Ch), dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="exchange coefficient for scalars",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"rs",       sngl(vars%rs),   dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="surface resistance",units="s/m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"betas",    sngl(vars%betas), dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="surface evaporation coefficient",units="1",missing_value=missing_value,ncid=ncid)

    endif


    call nc_write(fnm,"albsnw",   sngl(vars%albsnw),  dims=[dim_lon,dim_lat,dim_nsoil,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="snow albedo",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"swe",      sngl(vars%swe),dims=[dim_lon,dim_lat,dim_nsoil,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="snow water equivalent",units="kg/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"swe_max",  sngl(vars%swe_max),dims=[dim_lon,dim_lat,dim_nsoil,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="seasonal maximum snow water equivalent",units="kg/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"hsnow",    sngl(vars%hsnow),dims=[dim_lon,dim_lat,dim_nsoil,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="snow thickness",units="m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"tsnow",    sngl(vars%tsnow),dims=[dim_lon,dim_lat,dim_nsoil,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="snow layer temperature",units="K",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"tlake",     sngl(vars%tlake_sur), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="top layer lake temperature",units="K",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"cti_lim",  sngl(vars%cti_lim), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="critical CTI index",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fwet",     sngl(vars%fwet), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="saturated fraction",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"wtab",     sngl(vars%wtab), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="water table depth",units="m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"inf",      sngl(vars%inf),  dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="infiltration",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"runoff",     sngl(vars%runoff), dims=[dim_lon,dim_lat,dim_nsoil,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="total runoff (surface runoff + drainage)",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"runsur",     sngl(vars%runsur), dims=[dim_lon,dim_lat,dim_nsoil,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="surface unoff",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"calving",     sngl(vars%calving), dims=[dim_lon,dim_lat,dim_nsoil,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="calving",units="kg/m2/day weq",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"drain",    sngl(vars%drain),dims=[dim_lon,dim_lat,dim_nsoil,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="drainage",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"rsursub",  sngl(vars%rsursub),dims=[dim_lon,dim_lat,dim_nsoil,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="ratio surface runoff to drainage",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"snowmelt",  sngl(vars%snowmelt),dims=[dim_lon,dim_lat,dim_nsoil,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="snowmelt",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"icemelt",  sngl(vars%icemelt),dims=[dim_lon,dim_lat,dim_nsoil,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="icemelt",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"vpd",         sngl(vars%vpd), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="vapor pressure deficit",units="Pa",missing_value=missing_value,ncid=ncid)

    call nc_write(fnm,"dust_e_d",     sngl(vars%dust_e_d), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="dust emission from desert",units="g/m2/mon",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"dust_e_g",     sngl(vars%dust_e_g), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="dust emission from grassland",units="g/m2/mon",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"dust_e_s",     sngl(vars%dust_e_s), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="dust emission from shrubland",units="g/m2/mon",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"dust_e",     sngl(vars%dust_e), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="total dust emission",units="g/m2/mon",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"dust_dep",     sngl(vars%dust_dep), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="total dust deposition",units="g/m2/mon",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"dust_con", sngl(vars%dust_con),dims=[dim_lon,dim_lat,dim_nsoil,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="dust_concentration in snow",units="mg/kg",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"snow_grain", sngl(vars%snow_grain),dims=[dim_lon,dim_lat,dim_nsoil,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="snow grain size",units="micro m",missing_value=missing_value,ncid=ncid)

    ! grid cell averages

    call nc_write(fnm,"etot_g",     sngl(vars_g%etot), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="evapotranspiration",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ecan_g",     sngl(vars_g%ecan), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="canopy evaporation",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"esur_g",     sngl(vars_g%esur), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="surface evaporation",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"trans_g",    sngl(vars_g%trans),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="transpiration",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"le_g",       sngl(vars_g%le),   dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="latent heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"sh_g",       sngl(vars_g%sh),   dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="sensible heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"slh_g",      sngl(vars_g%slh),   dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="sum of sensible and latent heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ef_g",       sngl(vars_g%ef),   dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="evaporative fraction",units="1",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"lwnet_g",    sngl(vars_g%lwnet),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="net longwave radiation at the surface",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"swnet_g",    sngl(vars_g%swnet),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="net shortwave radiation at the surface",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"g_g",        sngl(vars_g%g),    dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="ground heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"flx_melt_g", sngl(vars_g%flx_melt),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="heat flux going into snowmelt",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"tskin_g",    sngl(vars_g%tskin),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="skin temperature",units="K",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"tskin_amp_g",    sngl(vars_g%tskin_amp),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="amplitude of diurnal cycle of skin temperature",units="K",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"alb_g",      sngl(vars_g%alb),  dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="surface albedo",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"alb_dir_g",  sngl(vars_g%alb_dir),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="surface albedo",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"alb_dif_g",  sngl(vars_g%alb_dif),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="surface diffuse albedo",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ra_g",       sngl(vars_g%ra),   dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="aerodynamic resistance",units="s/m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"rag_g",      sngl(vars_g%rag),  dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="aerodynamic resistance below canopy",units="s/m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Ri_g",     sngl(vars_g%Ri), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="bulk Richardson number",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Ch_g",     sngl(vars_g%Ch), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="exchange coefficient for scalars",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"rs_g",       sngl(vars_g%rs),   dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="surface resistance",units="s/m",missing_value=missing_value,ncid=ncid)


    ! seasonally invariant variables
    if( ndat .eq. nmon_year+1 ) then
     call nc_write(fnm,"gdd5",       sngl(vars%gdd5),  dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="growing degree days above 5 degC",units="K",missing_value=missing_value,ncid=ncid)
     call nc_write(fnm,"t2m_min_mon",sngl(vars%t2m_min_mon),  dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="near-surface air temperature of the coldest month",units="degC",missing_value=missing_value,ncid=ncid)
     call nc_write(fnm,"alt",        sngl(vars%alt),  dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="active layer thickness",units="m",missing_value=missing_value,ncid=ncid)
     call nc_write(fnm,"fwetmax",    sngl(vars%fwetmax),  dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="maxmiumum monthly wetland extent",units="m",missing_value=missing_value,ncid=ncid)
     call nc_write(fnm,"fwetmaxpot",    sngl(vars%fwetmaxpot), dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="potential maximum wetland extent",units="m",missing_value=missing_value,ncid=ncid)
    endif


   return
    
  end subroutine surf_nc_write
   
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s u r f _ a v e
  ! Purpose  :  Average (or sum) the paladyn fields as needed
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine surf_ave(d,ave,d_g,ave_g)
    
    implicit none
    
    type(surf_out) :: d(:), ave
    type(surf_g_out) :: d_g(:), ave_g

    integer :: k, kk, n
    real(wp) :: div
    
    n = size(d)
    div = real(n,wp)
   
    ! Set all values to zero
    ave%frac_surf = 0._wp
    ave%etot    = 0._wp
    ave%ecan    = 0._wp
    ave%esur    = 0._wp
    ave%trans   = 0._wp
    ave%le      = 0._wp
    ave%sh      = 0._wp
    ave%slh     = 0._wp
    ave%ef      = 0._wp
    ave%lwnet   = 0._wp
    ave%swnet   = 0._wp
    ave%g       = 0._wp
    ave%flx_melt= 0._wp
    ave%tskin   = 0._wp
    ave%tskin_amp= 0._wp
    ave%alb     = 0._wp
    ave%alb_dir = 0._wp
    ave%alb_dif = 0._wp
    ave%albsnw  = 0._wp
    ave%scf     = 0._wp
    ave%fcanwat = 0._wp
    ave%fcansn  = 0._wp
    ave%swe     = 0._wp
    ave%swe_max = 0._wp
    ave%hsnow   = 0._wp
    ave%tsnow   = 0._wp
    ave%tlake_sur = 0._wp
    ave%wind    = 0._wp      
    ave%ra      = 0._wp      
    ave%rag     = 0._wp
    ave%Ri      = 0._wp
    ave%z0m     = 0._wp
    ave%z0h     = 0._wp
    ave%Ch      = 0._wp
    ave%rs      = 0._wp
    ave%betas   = 0._wp
    ave%cti_lim = 0._wp
    ave%fwet    = 0._wp
    ave%wtab    = 0._wp
    ave%inf     = 0._wp
    ave%runoff  = 0._wp
    ave%runsur  = 0._wp
    ave%calving = 0._wp
    ave%drain   = 0._wp
    ave%rsursub = 0._wp
    ave%snowmelt= 0._wp
    ave%icemelt = 0._wp
    ave%vpd     = 0._wp
    ave%dust_e_d = 0._wp
    ave%dust_e_g = 0._wp
    ave%dust_e_s = 0._wp
    ave%dust_e = 0._wp
    ave%dust_dep = 0._wp
    ave%dust_con = 0._wp
    ave%snow_grain = 0._wp

    do k = 1, n
    d_g(k)%etot    = 0._wp
    d_g(k)%ecan    = 0._wp
    d_g(k)%esur    = 0._wp
    d_g(k)%trans   = 0._wp
    d_g(k)%le      = 0._wp
    d_g(k)%sh      = 0._wp
    d_g(k)%slh     = 0._wp
    d_g(k)%ef      = 0._wp
    d_g(k)%lwnet   = 0._wp
    d_g(k)%swnet   = 0._wp
    d_g(k)%g       = 0._wp
    d_g(k)%flx_melt= 0._wp
    d_g(k)%tskin   = 0._wp
    d_g(k)%tskin_amp= 0._wp
    d_g(k)%alb     = 0._wp
    d_g(k)%alb_dir = 0._wp
    d_g(k)%alb_dif = 0._wp
    d_g(k)%ra      = 0._wp
    d_g(k)%rag     = 0._wp
    d_g(k)%Ri      = 0._wp
    d_g(k)%Ch      = 0._wp
    d_g(k)%rs      = 0._wp
    enddo

    ave_g%etot    = 0._wp
    ave_g%ecan    = 0._wp
    ave_g%esur    = 0._wp
    ave_g%trans   = 0._wp
    ave_g%le      = 0._wp
    ave_g%sh      = 0._wp
    ave_g%slh     = 0._wp
    ave_g%ef      = 0._wp
    ave_g%lwnet   = 0._wp
    ave_g%swnet   = 0._wp
    ave_g%g       = 0._wp
    ave_g%flx_melt= 0._wp
    ave_g%tskin   = 0._wp
    ave_g%tskin_amp= 0._wp
    ave_g%alb     = 0._wp
    ave_g%alb_dir = 0._wp
    ave_g%alb_dif = 0._wp
    ave_g%ra      = 0._wp
    ave_g%rag     = 0._wp
    ave_g%Ri      = 0._wp
    ave_g%Ch      = 0._wp
    ave_g%rs      = 0._wp

    ! Loop over the time indices to sum up and average (if necessary)
    do k = 1, n
      ave%frac_surf= ave%frac_surf   + d(k)%frac_surf    / div
      ave%etot    = ave%etot       + d(k)%etot      / div
      ave%ecan    = ave%ecan       + d(k)%ecan      / div
      ave%esur    = ave%esur       + d(k)%esur      / div
      ave%trans   = ave%trans      + d(k)%trans     / div
      ave%le      = ave%le         + d(k)%le        / div
      ave%sh      = ave%sh         + d(k)%sh        / div
      ave%slh     = ave%slh        + d(k)%slh       / div
      ave%ef      = ave%ef         + d(k)%ef        / div
      ave%lwnet   = ave%lwnet      + d(k)%lwnet     / div
      ave%swnet   = ave%swnet      + d(k)%swnet     / div
      ave%g       = ave%g          + d(k)%g         / div
      ave%flx_melt= ave%flx_melt   + d(k)%flx_melt  / div
      ave%tskin   = ave%tskin      + d(k)%tskin     / div
      ave%tskin_amp= ave%tskin_amp + d(k)%tskin_amp / div
      ave%alb     = ave%alb        + d(k)%alb       / div
      ave%alb_dir = ave%alb_dir    + d(k)%alb_dir   / div
      ave%alb_dif = ave%alb_dif    + d(k)%alb_dif   / div
      ave%albsnw  = ave%albsnw     + d(k)%albsnw    / div
      ave%scf     = ave%scf        + d(k)%scf       / div
      ave%fcanwat = ave%fcanwat    + d(k)%fcanwat   / div
      ave%fcansn  = ave%fcansn     + d(k)%fcansn    / div
      ave%swe     = ave%swe        + d(k)%swe       / div
      ave%swe_max = ave%swe_max    + d(k)%swe_max   / div
      ave%hsnow   = ave%hsnow      + d(k)%hsnow     / div
      ave%tsnow   = ave%tsnow      + d(k)%tsnow     / div
      ave%tlake_sur   = ave%tlake_sur      + d(k)%tlake_sur     / div
      ave%wind    = ave%wind       + d(k)%wind      / div
      ave%ra      = ave%ra         + d(k)%ra        / div
      ave%rag     = ave%rag        + d(k)%rag       / div
      ave%Ri      = ave%Ri         + d(k)%Ri        / div
      ave%z0m     = ave%z0m        + d(k)%z0m       / div
      ave%z0h     = ave%z0h        + d(k)%z0h       / div
      ave%Ch      = ave%Ch         + d(k)%Ch        / div
      ave%rs      = ave%rs         + d(k)%rs        / div
      ave%betas   = ave%betas      + d(k)%betas     / div
      ave%cti_lim = ave%cti_lim    + d(k)%cti_lim   / div
      ave%fwet    = ave%fwet       + d(k)%fwet      / div
      ave%wtab    = ave%wtab       + d(k)%wtab      / div
      ave%inf     = ave%inf        + d(k)%inf       / div
      ave%runoff  = ave%runoff     + d(k)%runoff    / div
      ave%runsur  = ave%runsur     + d(k)%runsur    / div
      ave%calving = ave%calving    + d(k)%calving   / div
      ave%drain   = ave%drain      + d(k)%drain     / div
      ave%rsursub = ave%rsursub    + d(k)%rsursub   / div
      ave%snowmelt= ave%snowmelt   + d(k)%snowmelt  / div
      ave%icemelt = ave%icemelt    + d(k)%icemelt   / div
      ave%vpd     = ave%vpd        + d(k)%vpd      / div
      ave%dust_e_d    = ave%dust_e_d       + d(k)%dust_e_d      / div
      ave%dust_e_g    = ave%dust_e_g       + d(k)%dust_e_g      / div
      ave%dust_e_s    = ave%dust_e_s       + d(k)%dust_e_s      / div
      ave%dust_e    = ave%dust_e       + d(k)%dust_e      / div
      ave%dust_dep    = ave%dust_dep       + d(k)%dust_dep      / div
      ave%dust_con    = ave%dust_con       + d(k)%dust_con      / div
      ave%snow_grain  = ave%snow_grain     + d(k)%snow_grain    / div
    end do
      
    ! grid cell averages
    do k=1,n
     do kk=1,nsurf
      d_g(k)%etot   = d_g(k)%etot   + d(k)%etot(:,:,kk)  * d(k)%frac_surf(:,:,kk)
      d_g(k)%ecan   = d_g(k)%ecan   + d(k)%ecan(:,:,kk)  * d(k)%frac_surf(:,:,kk)
      d_g(k)%esur   = d_g(k)%esur   + d(k)%esur(:,:,kk)  * d(k)%frac_surf(:,:,kk)
      d_g(k)%trans  = d_g(k)%trans  + d(k)%trans(:,:,kk) * d(k)%frac_surf(:,:,kk)
      d_g(k)%le     = d_g(k)%le     + d(k)%le(:,:,kk)    * d(k)%frac_surf(:,:,kk)
      d_g(k)%sh     = d_g(k)%sh     + d(k)%sh(:,:,kk)    * d(k)%frac_surf(:,:,kk)
      d_g(k)%slh    = d_g(k)%slh    + d(k)%slh(:,:,kk)   * d(k)%frac_surf(:,:,kk)
      d_g(k)%ef     = d_g(k)%ef     + d(k)%ef(:,:,kk)    * d(k)%frac_surf(:,:,kk)
      d_g(k)%lwnet  = d_g(k)%lwnet  + d(k)%lwnet(:,:,kk) * d(k)%frac_surf(:,:,kk)
      d_g(k)%swnet  = d_g(k)%swnet  + d(k)%swnet(:,:,kk) * d(k)%frac_surf(:,:,kk)
      d_g(k)%g      = d_g(k)%g      + d(k)%g(:,:,kk)     * d(k)%frac_surf(:,:,kk)
      d_g(k)%flx_melt= d_g(k)%flx_melt+ d(k)%flx_melt(:,:,kk)     * d(k)%frac_surf(:,:,kk)
      d_g(k)%tskin  = d_g(k)%tskin  + d(k)%tskin(:,:,kk) * d(k)%frac_surf(:,:,kk)
      d_g(k)%tskin_amp = d_g(k)%tskin_amp + d(k)%tskin_amp(:,:,kk) * d(k)%frac_surf(:,:,kk)
      d_g(k)%alb    = d_g(k)%alb    + d(k)%alb(:,:,kk)   * d(k)%frac_surf(:,:,kk)
      d_g(k)%alb_dir= d_g(k)%alb_dir+ d(k)%alb_dir(:,:,kk)* d(k)%frac_surf(:,:,kk)
      d_g(k)%alb_dif= d_g(k)%alb_dif+ d(k)%alb_dif(:,:,kk)* d(k)%frac_surf(:,:,kk)
      d_g(k)%ra     = d_g(k)%ra     + d(k)%ra(:,:,kk)    * d(k)%frac_surf(:,:,kk)
      d_g(k)%Ri     = d_g(k)%Ri     + d(k)%Ri(:,:,kk)    * d(k)%frac_surf(:,:,kk)
      d_g(k)%Ch     = d_g(k)%Ch     + d(k)%Ch(:,:,kk)    * d(k)%frac_surf(:,:,kk)
      d_g(k)%rs     = d_g(k)%rs     + d(k)%rs(:,:,kk)    * d(k)%frac_surf(:,:,kk)
     enddo
     do kk=1,npft
      d_g(k)%rag    = d_g(k)%rag    + d(k)%rag(:,:,kk)   * d(k)%frac_surf(:,:,kk)
     enddo
    enddo

    do k=1,nsurf
      ave_g%etot   = ave_g%etot   + ave%etot(:,:,k)  * ave%frac_surf(:,:,k)
      ave_g%ecan   = ave_g%ecan   + ave%ecan(:,:,k)  * ave%frac_surf(:,:,k)
      ave_g%esur   = ave_g%esur   + ave%esur(:,:,k)  * ave%frac_surf(:,:,k)
      ave_g%trans  = ave_g%trans  + ave%trans(:,:,k) * ave%frac_surf(:,:,k)
      ave_g%le     = ave_g%le     + ave%le(:,:,k)    * ave%frac_surf(:,:,k)
      ave_g%sh     = ave_g%sh     + ave%sh(:,:,k)    * ave%frac_surf(:,:,k)
      ave_g%slh    = ave_g%slh    + ave%slh(:,:,k)   * ave%frac_surf(:,:,k)
      ave_g%ef     = ave_g%ef     + ave%ef(:,:,k)    * ave%frac_surf(:,:,k)
      ave_g%lwnet  = ave_g%lwnet  + ave%lwnet(:,:,k) * ave%frac_surf(:,:,k)
      ave_g%swnet  = ave_g%swnet  + ave%swnet(:,:,k) * ave%frac_surf(:,:,k)
      ave_g%g      = ave_g%g      + ave%g(:,:,k)     * ave%frac_surf(:,:,k)
      ave_g%flx_melt= ave_g%flx_melt + ave%flx_melt(:,:,k)     * ave%frac_surf(:,:,k)
      ave_g%tskin  = ave_g%tskin  + ave%tskin(:,:,k) * ave%frac_surf(:,:,k)
      ave_g%tskin_amp = ave_g%tskin_amp + ave%tskin_amp(:,:,k) * ave%frac_surf(:,:,k)
      ave_g%alb    = ave_g%alb    + ave%alb(:,:,k)   * ave%frac_surf(:,:,k)
      ave_g%alb_dir= ave_g%alb_dir+ ave%alb_dir(:,:,k)* ave%frac_surf(:,:,k)
      ave_g%alb_dif= ave_g%alb_dif+ ave%alb_dif(:,:,k)* ave%frac_surf(:,:,k)
      ave_g%ra     = ave_g%ra     + ave%ra(:,:,k)    * ave%frac_surf(:,:,k)
      ave_g%Ri     = ave_g%Ri     + ave%Ri(:,:,k)    * ave%frac_surf(:,:,k)
      ave_g%Ch     = ave_g%Ch     + ave%Ch(:,:,k)    * ave%frac_surf(:,:,k)
      ave_g%rs     = ave_g%rs     + ave%rs(:,:,k)    * ave%frac_surf(:,:,k)
    enddo
    do k=1,npft
      ave_g%rag    = ave_g%rag    + ave%rag(:,:,k)   * ave%frac_surf(:,:,k)
    enddo

   return
    
  end subroutine surf_ave
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s u r f _ d a i l y _ n c
  ! Purpose  :  Initialize netcdf output for paladyn
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine surf_daily_nc(fnm)
    
    implicit none
    
    character (len=*) :: fnm
    integer :: ncid
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,dim_time, x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE.,ncid=ncid)    
    call nc_write_dim(fnm, dim_day, x=1._wp, dx=1._wp, nx=nday_year, axis="e", &
    units="days", ncid=ncid)
    call nc_write_dim(fnm, dim_nsurf, x=i_surf, units="n/a", axis="z", ncid=ncid)
    call nc_write_dim(fnm, dim_lat, x=lat, axis="y", units="degrees_north", ncid=ncid)
    call nc_write_dim(fnm, dim_lon, x=lon, axis="x", units="degrees_east", ncid=ncid)
    call nc_write_dim(fnm, dim_npft, x=i_pft, units="n/a", ncid=ncid)
    call nc_write_dim(fnm, dim_nsoil, x=i_soil, units="n/a", ncid=ncid)
    call nc_close(ncid)

    return
  
  end subroutine surf_daily_nc
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s u r f _ d a i l y _ n c _ w r i t e
  ! Purpose  :  Output timestep of netcdf output for paladyn
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine surf_daily_nc_write(fnm,ncid,vars,ndat,nout)
    
    implicit none
    
    type(surf_out) :: vars
 
    character (len=*) :: fnm
    integer :: ndat, nout, ncid

    
    call nc_write(fnm,"frac_surf",       sngl(vars%frac_surf),   dims=[dim_lon,dim_lat,dim_nsurf,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="latent heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"le",       sngl(vars%le),   dims=[dim_lon,dim_lat,dim_nsurf,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="latent heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"sh",       sngl(vars%sh),   dims=[dim_lon,dim_lat,dim_nsurf,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="sensible heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"swnet",    sngl(vars%swnet),dims=[dim_lon,dim_lat,dim_nsurf,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="net shortwave radiation at the surface",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"lwnet",    sngl(vars%lwnet),dims=[dim_lon,dim_lat,dim_nsurf,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="net longwave radiation at the surface",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"g",        sngl(vars%g),    dims=[dim_lon,dim_lat,dim_nsurf,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="ground heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"flx_melt",sngl(vars%flx_melt),dims=[dim_lon,dim_lat,dim_nsurf,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="ground heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"tskin",    sngl(vars%tskin),dims=[dim_lon,dim_lat,dim_nsurf,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="skin temperature",units="K",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ra",        sngl(vars%ra),    dims=[dim_lon,dim_lat,dim_nsurf,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="ground heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"rs",        sngl(vars%rs),    dims=[dim_lon,dim_lat,dim_nsurf,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="ground heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"betas",     sngl(vars%betas),    dims=[dim_lon,dim_lat,dim_nsurf,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="ground heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"rag",        sngl(vars%rag),    dims=[dim_lon,dim_lat,dim_npft,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,npft,1,1],long_name="ground heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"Ri",      sngl(vars%Ri),    dims=[dim_lon,dim_lat,dim_nsurf,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="ground heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ecan",     sngl(vars%ecan),    dims=[dim_lon,dim_lat,dim_nsurf,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="ground heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"msnow",    float(vars%msnow),dims=[dim_lon,dim_lat,dim_nsoil,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="snow thickness",units="m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"hsnow",    sngl(vars%hsnow),dims=[dim_lon,dim_lat,dim_nsoil,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="snow thickness",units="m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"w_snow",    sngl(vars%swe),dims=[dim_lon,dim_lat,dim_nsoil,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="snow water equivalent",units="kg/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"w_snow_max",    sngl(vars%swe_max),dims=[dim_lon,dim_lat,dim_nsoil,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="max snow water equivalent",units="kg/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"albsnw",    sngl(vars%albsnw),dims=[dim_lon,dim_lat,dim_nsoil,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="snow thickness",units="m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"tsnow",    sngl(vars%tsnow),dims=[dim_lon,dim_lat,dim_nsoil,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="snow layer temperature",units="K",missing_value=missing_value,ncid=ncid)

    call nc_write(fnm,"le_g",       sngl(sum(vars%le*vars%frac_surf,3)),   dims=[dim_lon,dim_lat,dim_day,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="latent heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"sh_g",       sngl(sum(vars%sh*vars%frac_surf,3)),   dims=[dim_lon,dim_lat,dim_day,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="latent heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"lwnet_g",       sngl(sum(vars%lwnet*vars%frac_surf,3)),   dims=[dim_lon,dim_lat,dim_day,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="latent heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"g_g",       sngl(sum(vars%g*vars%frac_surf,3)),   dims=[dim_lon,dim_lat,dim_day,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="latent heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"flx_melt_g",       sngl(sum(vars%flx_melt*vars%frac_surf,3)),   dims=[dim_lon,dim_lat,dim_day,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="latent heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"tskin_g",       sngl(sum(vars%tskin*vars%frac_surf,3)),   dims=[dim_lon,dim_lat,dim_day,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="latent heat flux",units="W/m2",missing_value=missing_value,ncid=ncid)

   return
    
  end subroutine surf_daily_nc_write
   

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c o n s _ n c
  ! Purpose  :  Initialize netcdf output for paladyn
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine cons_nc(fnm)
    
    implicit none
    
    character (len=*) :: fnm
    integer :: ncid
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,dim_time, x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE.,ncid=ncid)
    call nc_write_dim(fnm, dim_month, x=1._wp, dx=1._wp, nx=13, units="months", ncid=ncid)
    call nc_write_dim(fnm, dim_nsurf, x=i_surf, units="n/a", axis="z", ncid=ncid)
    call nc_write_dim(fnm, dim_lat, x=lat, axis="y", units="degrees_north", ncid=ncid)
    call nc_write_dim(fnm, dim_lon, x=lon, axis="x", units="degrees_east", ncid=ncid)
    call nc_write_dim(fnm, dim_nsoil, x=i_soil, units="n/a", ncid=ncid)
    call nc_write_dim(fnm, dim_ncarb, x=1._wp, dx=1._wp, nx=ncarb, units="n/a", ncid=ncid)
    call nc_close(ncid)
    
    return
  
  end subroutine cons_nc
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c o n s _ n c _ w r i t e
  ! Purpose  :  Output timestep of netcdf output for paladyn
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine cons_nc_write(fnm,ncid,vars,ndat,nout)
    
    implicit none
    
    type(cons_out) :: vars
 
    character (len=*) :: fnm
    integer :: ndat, nout, ncid


    call nc_write(fnm,"econs_su1", sngl(vars%econs_su1), dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="surface energy balance",units="w/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"econs_su2", sngl(vars%econs_su2), dims=[dim_lon,dim_lat,dim_nsurf,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsurf,1,1],long_name="surface energy balance",units="w/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"econs_so",  sngl(vars%econs_so),  dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="soil energy balance",units="w/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"wcons",     sngl(vars%wcons),     dims=[dim_lon,dim_lat,dim_nsoil,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nsoil,1,1],long_name="water balance",units="kg/m2",missing_value=missing_value,ncid=ncid)

    call nc_write(fnm,"ccons_s",     sngl(vars%ccons_s),   dims=[dim_lon,dim_lat,dim_ncarb,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,ncarb,1,1],long_name="soil carbon balance",units="kgC/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ccons13_s",   sngl(vars%ccons13_s), dims=[dim_lon,dim_lat,dim_ncarb,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,ncarb,1,1],long_name="soil carbon 13 balance",units="kgC/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ccons14_s",   sngl(vars%ccons14_s), dims=[dim_lon,dim_lat,dim_ncarb,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,ncarb,1,1],long_name="soil carbon 14 balance",units="kgC/m2",missing_value=missing_value,ncid=ncid)


    call nc_write(fnm,"ccons_v",     sngl(vars%ccons_v),   dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="vegetation carbon balance",units="kgC/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ccons13_v",   sngl(vars%ccons13_v), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="vegetation carbon 13 balance",units="kgC/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ccons14_v",   sngl(vars%ccons14_v), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="vegetation carbon 14 balance",units="kgC/m2",missing_value=missing_value,ncid=ncid)


   return
    
  end subroutine cons_nc_write
   
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c o n s _ a v e
  ! Purpose  :  Average (or sum) the paladyn fields as needed
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine cons_ave(d,ave)
    
    implicit none
    
    type(cons_out) :: d(:), ave

    integer :: k, n
    real(wp) :: div
    
    n = size(d)
    div = real(n,wp)
   
    ! Set all values to zero
    ave%econs_su1= 0._wp
    ave%econs_su2= 0._wp
    ave%econs_so = 0._wp
    ave%wcons    = 0._wp
    ave%ccons_s  = 0._wp
    ave%ccons_v  = 0._wp
    ave%ccons13_s  = 0._wp
    ave%ccons13_v  = 0._wp
    ave%ccons14_s  = 0._wp
    ave%ccons14_v  = 0._wp

    ! Loop over the time indices to sum up and average (if necessary)
    do k = 1, n
      ave%econs_su1= ave%econs_su1   + d(k)%econs_su1  / div
      ave%econs_su2= ave%econs_su2   + d(k)%econs_su2  / div
      ave%econs_so = ave%econs_so    + d(k)%econs_so   / div
      ave%wcons    = ave%wcons       + d(k)%wcons      / div
      ave%ccons_s  = ave%ccons_s     + d(k)%ccons_s    / div
      ave%ccons_v  = ave%ccons_v     + d(k)%ccons_v    / div
      ave%ccons13_s  = ave%ccons13_s     + d(k)%ccons13_s    / div
      ave%ccons13_v  = ave%ccons13_v     + d(k)%ccons13_v    / div
      ave%ccons14_s  = ave%ccons14_s     + d(k)%ccons14_s    / div
      ave%ccons14_v  = ave%ccons14_v     + d(k)%ccons14_v    / div
    end do
      
   return
    
  end subroutine cons_ave
  

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c a r b o n _ n c
  ! Purpose  :  Initialize netcdf output for paladyn
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine carbon_nc(fnm)
    
    implicit none
    
    character (len=*) :: fnm
    integer :: ncid
    real(wp) :: empty_time(0)
 
    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm, dim_time, x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE.,ncid=ncid)
    call nc_write_dim(fnm, dim_month, x=1._wp, dx=1._wp, nx=13, axis="e", &
    units="months", ncid=ncid)
    call nc_write_dim(fnm, dim_npft, x=i_pft, axis="z", units="n/a", ncid=ncid)
    call nc_write_dim(fnm, dim_lon, x=lon, axis="x", ncid=ncid)
    call nc_write_dim(fnm, dim_lat, x=lat, axis="y", ncid=ncid)
    call nc_write_dim(fnm, dim_depth, x=z_c(1:nlc), units="m", axis="z", ncid=ncid)
    call nc_write_dim(fnm, dim_ncarb, x=1._wp, dx=1._wp, nx=ncarb, units="n/a", ncid=ncid)
    call nc_write_dim(fnm,dim_nlit,x=1._wp,dx=1._wp,nx=nlit,units="n/a",ncid=ncid)
    call nc_close(ncid)
 
    return
  
  end subroutine carbon_nc
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c a r b o n _ n c _ w r i t e
  ! Purpose  :  Output timestep of netcdf output for paladyn
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine carbon_nc_write(fnm,ncid,vars,vars_g,ndat,nout)
    
    implicit none
    
    type(carbon_out) :: vars
    type(carbon_g_out) :: vars_g
 
    character (len=*) :: fnm
    integer :: ndat, nout, ncid
    

    call nc_write(fnm,"Cflx_atm_lnd",      sngl(vars_g%Cflx_atm_lnd),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="net land carbon flux",units="kgC/m2/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"C13flx_atm_lnd",      sngl(vars_g%C13flx_atm_lnd),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="net land carbon 13 flux",units="kgC/m2/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"C14flx_atm_lnd",      sngl(vars_g%C14flx_atm_lnd),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="net land carbon 14 flux",units="kgC/m2/yr",missing_value=missing_value,ncid=ncid)
!    call nc_write(fnm,"d13Cflx_atm_lnd",    sngl(vars_g%d13Cflx_atm_lnd),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="d13C of atmosphere-land carbon flux",units="permil",missing_value=missing_value,ncid=ncid)
!    call nc_write(fnm,"D14Cflx_atm_lnd",    sngl(vars_g%D14Cflx_atm_lnd),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="D14C of atmosphere-land carbon flux",units="permil",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"gpp",       sngl(vars%gpp),   dims=[dim_lon,dim_lat,dim_npft,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,npft,1,1],long_name="gross primary productivity",units="kgC/m2/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"npp",       sngl(vars%npp),   dims=[dim_lon,dim_lat,dim_npft,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,npft,1,1],long_name="net primary productivity",units="kgC/m2/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"aresp",     sngl(vars%aresp), dims=[dim_lon,dim_lat,dim_npft,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,npft,1,1],long_name="autotrophic respiration",units="kgC/m2/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"lai",       sngl(vars%lai),   dims=[dim_lon,dim_lat,dim_npft,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,npft,1,1],long_name="leaf area index",units="m2/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"xi",        sngl(vars%xi),    dims=[dim_lon,dim_lat,dim_npft,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,npft,1,1],long_name="ratio of leaf internal to ambient partial pressure of CO2",units="1",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"wue",       sngl(vars%wue),   dims=[dim_lon,dim_lat,dim_npft,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,npft,1,1],long_name="intrinsic water use efficiency",units="micro mol/mol",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"gcan",      sngl(vars%gcan),  dims=[dim_lon,dim_lat,dim_npft,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,npft,1,1],long_name="canopy conductance",units="m/s",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"litter",    sngl(vars%litter), dims=[dim_lon,dim_lat,dim_ncarb,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,ncarb,1,1],long_name="litterfall",units="kgC/m2/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"sresp",     sngl(vars%sresp),  dims=[dim_lon,dim_lat,dim_ncarb,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,ncarb,1,1],long_name="soil respiration",units="kgC/m2/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ch4",       sngl(vars%ch4),    dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="methane emissions",units="gCH4/m2 gridcell/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ch4wet",    sngl(vars%ch4wet), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="methane emissions from wetlands",units="gCH4/m2 gridcell/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ch4shelf",  sngl(vars%ch4shelf), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="methane emissions from ocean shelf",units="gCH4/m2 gridcell/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ch4lake",   sngl(vars%ch4lake), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="methane emissions from lakes",units="gCH4/m2 gridcell/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ch4peat",   sngl(vars%ch4peat),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="methane emissions from peatlands",units="gCH4/m2 gridcell/yr",missing_value=missing_value,ncid=ncid)

    ! grid cell averages

    call nc_write(fnm,"gpp_g",       sngl(vars_g%gpp),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="gross primary productivity",units="kgC/m2/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"npp_g",       sngl(vars_g%npp),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="net primary productivity",units="kgC/m2/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"aresp_g",     sngl(vars_g%aresp),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="autotrophic respiration",units="kgC/m2/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"lai_g",       sngl(vars_g%lai),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="leaf area index",units="m2/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"gcan_g",      sngl(vars_g%gcan),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="canopy conductance",units="m/s",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"sresp_g",     sngl(vars_g%sresp),dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="soil respiration",units="kgC/m2/yr",missing_value=missing_value,ncid=ncid)

    if (ndat.eq.13) then

    call nc_write(fnm,"seeds",     sngl(vars%seeds),  dims=[dim_lon,dim_lat,dim_npft,dim_time],start=[1,1,1,nout],count=[nx,ny,npft,1],long_name="seed fraction",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"pfts",      sngl(vars%pfts),  dims=[dim_lon,dim_lat,dim_npft,dim_time],start=[1,1,1,nout],count=[nx,ny,npft,1],long_name="pft fraction of icefree area",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"bare",      sngl(vars%bare),  dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="bare soil fraction of icefree area",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"lambda",sngl(vars%lambda),  dims=[dim_lon,dim_lat,dim_npft,dim_time],start=[1,1,1,nout],count=[nx,ny,npft,1],long_name="NPP partitioning factor",units="1",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"gamma_dist",sngl(vars%gamma_dist),  dims=[dim_lon,dim_lat,dim_npft,dim_time],start=[1,1,1,nout],count=[nx,ny,npft,1],long_name="disturbance rate",units="years",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"gamma_luc",sngl(vars%gamma_luc),  dims=[dim_lon,dim_lat,dim_npft,dim_time],start=[1,1,1,nout],count=[nx,ny,npft,1],long_name="disturbance rate from land use change",units="years",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"gamma_ice",sngl(vars%gamma_ice),  dims=[dim_lon,dim_lat,dim_npft,dim_time],start=[1,1,1,nout],count=[nx,ny,npft,1],long_name="disturbance rate from ice sheets",units="years",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"vegc",      sngl(vars%vegc),  dims=[dim_lon,dim_lat,dim_npft,dim_time],start=[1,1,1,nout],count=[nx,ny,npft,1],long_name="vegetation carbon",units="kgC/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"soilc",     sngl(vars%soilc),  dims=[dim_lon,dim_lat,dim_ncarb,dim_time],start=[1,1,1,nout],count=[nx,ny,ncarb,1],long_name="soil carbon",units="kgC/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"sai",       sngl(vars%sai),   dims=[dim_lon,dim_lat,dim_npft,dim_time],start=[1,1,1,nout],count=[nx,ny,npft,1],long_name="stem area index",units="m2/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"veg_h",     sngl(vars%veg_h), dims=[dim_lon,dim_lat,dim_npft,dim_time],start=[1,1,1,nout],count=[nx,ny,npft,1],long_name="vegetation height",units="m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fpeat",       sngl(vars%fpeat),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="peatland fraction",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fpeatpot",    sngl(vars%fpeatpot),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="potential peatland fraction",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"dCpeatdt",    sngl(vars%dCpeatdt),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="peat carbon accumulation rate",units="gC/m2/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"acroh",       sngl(vars%acroh),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="acrotelm thickness",units="m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"catoh",       sngl(vars%catoh),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="catotelm thickness",units="m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"litterc_prof",     sngl(vars%litterc_prof),dims=[dim_lon,dim_lat,dim_depth,dim_ncarb,dim_time],start=[1,1,1,1,nout],count=[nx,ny,nlc,ncarb,1],long_name="litter carbon",units="kgC/m3",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fastc_prof",       sngl(vars%fastc_prof),  dims=[dim_lon,dim_lat,dim_depth,dim_ncarb,dim_time],start=[1,1,1,1,nout],count=[nx,ny,nlc,ncarb,1],long_name="fast soil carbon",units="kgC/m3",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"slowc_prof",       sngl(vars%slowc_prof),  dims=[dim_lon,dim_lat,dim_depth,dim_ncarb,dim_time],start=[1,1,1,1,nout],count=[nx,ny,nlc,ncarb,1],long_name="slow soil carbon",units="kgC/m3",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"soilc_prof",       sngl(vars%soilc_prof),  dims=[dim_lon,dim_lat,dim_depth,dim_ncarb,dim_time],start=[1,1,1,1,nout],count=[nx,ny,nlc,ncarb,1],long_name="total soil carbon",units="kgC/m3",missing_value=missing_value,ncid=ncid)

    call nc_write(fnm,"d13C_veg",     sngl(vars%d13c_veg),   dims=[dim_lon,dim_lat,dim_npft,dim_time],start=[1,1,1,nout],count=[nx,ny,npft,1],long_name="d13C vegetation carbon",units="permil",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"D14C_veg",     sngl(vars%D14c_veg),   dims=[dim_lon,dim_lat,dim_npft,dim_time],start=[1,1,1,nout],count=[nx,ny,npft,1],long_name="D14C vegetation carbon",units="permil",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"d13C_soilc",     sngl(vars%d13c_soil), dims=[dim_lon,dim_lat,dim_ncarb,dim_time],start=[1,1,1,nout],count=[nx,ny,ncarb,1],long_name="d13C soil carbon",units="permil",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"D14C_soilc",     sngl(vars%D14c_soil), dims=[dim_lon,dim_lat,dim_ncarb,dim_time],start=[1,1,1,nout],count=[nx,ny,ncarb,1],long_name="D14C soil carbon",units="permil",missing_value=missing_value,ncid=ncid)

    call nc_write(fnm,"litter_prof",     sngl(vars%litter_prof),dims=[dim_lon,dim_lat,dim_depth,dim_ncarb,dim_time],start=[1,1,1,1,nout],count=[nx,ny,nlc,ncarb,1],long_name="profile of litter input into the soil",units="kgC/m2/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"d13C_litter_prof",     sngl(vars%d13c_litter_prof),dims=[dim_lon,dim_lat,dim_depth,dim_ncarb,dim_time],start=[1,1,1,1,nout],count=[nx,ny,nlc,ncarb,1],long_name="d13C profile of litter input into the soil",units="permil",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"D14C_litter_prof",     sngl(vars%D14c_litter_prof),dims=[dim_lon,dim_lat,dim_depth,dim_ncarb,dim_time],start=[1,1,1,1,nout],count=[nx,ny,nlc,ncarb,1],long_name="D14C profile of litter input into the soil",units="permil",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"d13C_soilc_prof",     sngl(vars%d13c_prof),dims=[dim_lon,dim_lat,dim_depth,dim_ncarb,dim_time],start=[1,1,1,1,nout],count=[nx,ny,nlc,ncarb,1],long_name="d13C soil carbon profile",units="permil",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"D14C_soilc_prof",     sngl(vars%D14c_prof),dims=[dim_lon,dim_lat,dim_depth,dim_ncarb,dim_time],start=[1,1,1,1,nout],count=[nx,ny,nlc,ncarb,1],long_name="D14C soil carbon profile",units="permil",missing_value=missing_value,ncid=ncid)

    call nc_write(fnm,"doc_export",       sngl(vars%doc_export),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="dissolved organic output export through rivers",units="kgC/m2/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"poc_export",       sngl(vars%poc_export),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="particulate organic output export through rivers",units="kgC/m2/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"weath_carb",       sngl(vars%weath_carb),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="carbonate weathering",units="mol C/m2/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"weath_sil",       sngl(vars%weath_sil),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="silicate weathering",units="mol C/m2/yr",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"lithology",     sngl(vars%lithology),dims=[dim_lon,dim_lat,dim_nlit,dim_time],start=[1,1,1,nout],count=[nx,ny,nlit,1],long_name="lithology",units="/",missing_value=missing_value,ncid=ncid)

    ! grid cell averages
    call nc_write(fnm,"vegc_g",      sngl(vars_g%vegc),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="vegetation carbon",units="kgC/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"soilc_g",     sngl(vars_g%soilc),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="soil carbon",units="kgC/m2",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"d13C_veg_g",   sngl(vars_g%d13c_veg),   dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="d13C vegetation carbon",units="permil",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"D14C_veg_g",   sngl(vars_g%D14c_veg),   dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="D14C vegetation carbon",units="permil",missing_value=missing_value,ncid=ncid)
  endif

    call nc_write(fnm,"disc",         sngl(vars%disc),       dims=[dim_lon,dim_lat,dim_npft,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,npft,1,1],long_name="discrimination during photosynthesis",units="permil",missing_value=missing_value,ncid=ncid)

    ! grid cell averages

    call nc_write(fnm,"disc_g",       sngl(vars_g%disc),       dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="discrimination",units="permil",missing_value=missing_value,ncid=ncid)

    return
    
  end subroutine carbon_nc_write
   

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c a r b o n _ a v e
  ! Purpose  :  Average (or sum) the paladyn fields as needed
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine carbon_ave(d,ave,d_g,ave_g)
    
    implicit none
    
    type(carbon_out) :: d(:), ave
    type(carbon_g_out) :: d_g(:), ave_g

    integer :: k, n
    real(wp) :: div
    
    n = size(d)
    div = real(n,wp)
   
    ! Set all values to zero
    ave%lai      = 0._wp
    ave%xi       = 0._wp
    ave%wue      = 0._wp
    ave%gcan     = 0._wp
    ave%gpp      = 0._wp
    ave%npp      = 0._wp
    ave%aresp    = 0._wp
    ave%litter   = 0._wp
    ave%sresp    = 0._wp
    ave%ch4      = 0._wp
    ave%ch4wet   = 0._wp
    ave%ch4shelf = 0._wp
    ave%ch4lake  = 0._wp
    ave%ch4peat  = 0._wp
    ave%disc     = 0._wp

    ave_g%lai      = 0._wp
    ave_g%gcan     = 0._wp
    ave_g%gpp      = 0._wp
    ave_g%npp      = 0._wp
    ave_g%Cflx_atm_lnd = 0._wp
    ave_g%C13flx_atm_lnd = 0._wp
    ave_g%C14flx_atm_lnd = 0._wp
    ave_g%d13Cflx_atm_lnd = 0._wp
    ave_g%D14Cflx_atm_lnd = 0._wp
    ave_g%aresp    = 0._wp
    ave_g%sresp    = 0._wp
    ave_g%disc     = 0._wp

    ! Loop over the time indices to sum up and average (if necessary)
    do k = 1, n
      ave%lai      = ave%lai         + d(k)%lai    / div
      ave%xi       = ave%xi          + d(k)%xi     / div
      ave%wue      = ave%wue         + d(k)%wue    / div
      ave%gcan     = ave%gcan        + d(k)%gcan   / div
      ave%gpp      = ave%gpp         + d(k)%gpp    / div
      ave%npp      = ave%npp         + d(k)%npp    / div
      ave%aresp    = ave%aresp       + d(k)%aresp  / div
      ave%litter   = ave%litter      + d(k)%litter / div
      ave%sresp    = ave%sresp       + d(k)%sresp  / div
      ave%ch4   = ave%ch4      + d(k)%ch4 / div
      ave%ch4wet   = ave%ch4wet      + d(k)%ch4wet / div
      ave%ch4shelf  = ave%ch4shelf     + d(k)%ch4shelf / div
      ave%ch4lake   = ave%ch4lake      + d(k)%ch4lake / div
      ave%ch4peat   = ave%ch4peat      + d(k)%ch4peat / div
      ave%disc     = ave%disc        + d(k)%disc*d(k)%gpp / div
    end do
    where (ave%gpp>0._wp) 
      ave%disc = ave%disc/ave%gpp
    endwhere
      
    ! grid cell averages
    do k=1,n
     ave_g%npp    = ave_g%npp      + d_g(k)%npp    / div
     ave_g%Cflx_atm_lnd = ave_g%Cflx_atm_lnd + d_g(k)%Cflx_atm_lnd / div
     ave_g%C13flx_atm_lnd = ave_g%C13flx_atm_lnd + d_g(k)%C13flx_atm_lnd / div
     ave_g%C14flx_atm_lnd = ave_g%C14flx_atm_lnd + d_g(k)%C14flx_atm_lnd / div
     ave_g%d13Cflx_atm_lnd = ave_g%d13Cflx_atm_lnd + d_g(k)%d13Cflx_atm_lnd / div
     ave_g%D14Cflx_atm_lnd = ave_g%D14Cflx_atm_lnd + d_g(k)%D14Cflx_atm_lnd / div
     ave_g%sresp   = ave_g%sresp     + d_g(k)%sresp   / div
     ave_g%lai    = ave_g%lai    + d_g(k)%lai   /div 
     ave_g%gcan   = ave_g%gcan   + d_g(k)%gcan  /div
     ave_g%gpp    = ave_g%gpp    + d_g(k)%gpp   /div
     ave_g%aresp  = ave_g%aresp  + d_g(k)%aresp /div
     ave_g%disc   = ave_g%disc   + d_g(k)%disc*d_g(k)%gpp  /div
    enddo
    where (ave_g%gpp>0._wp) 
      ave_g%disc = ave_g%disc/ave_g%gpp
    endwhere


   return
    
  end subroutine carbon_ave
  

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s o i l _ n c
  ! Purpose  :  Initialize netcdf output for paladyn
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine soil_nc(fnm)
    
    implicit none
    
    character (len=*) :: fnm
    integer :: ncid
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm, dim_time, x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE., ncid=ncid)
    call nc_write_dim(fnm, dim_month, x=1._wp, dx=1._wp, nx=13, axis="e", &
    units="months", ncid=ncid)
    call nc_write_dim(fnm, dim_depth, x=z(1:nl), units="m", axis="z", ncid=ncid)
    call nc_write_dim(fnm, dim_lat, x=lat, axis="y", ncid=ncid)
    call nc_write_dim(fnm, dim_lon, x=lon, axis="x", ncid=ncid)
    call nc_close(ncid)

    return
  
  end subroutine soil_nc
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s o i l _ n c _ w r i t e
  ! Purpose  :  Output timestep of netcdf output for paladyn
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine soil_nc_write(fnm,ncid,vars,ndat,nout)
    
    implicit none
    
    type(soil_out) :: vars
 
    character (len=*) :: fnm
    integer :: ndat, nout, ncid
 

    call nc_write(fnm,"tsoil",     sngl(vars%tsoil),  dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nl,1,1],long_name="soil temperature",units="K",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"tice",     sngl(vars%tice),  dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nl,1,1],long_name="ice temperature",units="K",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"tsublake",     sngl(vars%tsublake),  dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nl,1,1],long_name="soil temperature below lake",units="K",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"thetaw",    sngl(vars%thetaw), dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nl,1,1],long_name="soil liquid water content",units="m3/m3",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"thetai",    sngl(vars%thetai), dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nl,1,1],long_name="soil frozen water content",units="m3/m3",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"theta",     sngl(vars%theta),  dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nl,1,1],long_name="total soil water content",units="m3/m3",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"thetas",    sngl(vars%thetas),  dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nl,1,1],long_name="soil porosity",units="m3/m3",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fthetas",   sngl(vars%fthetas),  dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nl,1,1],long_name="fraction of porosity filled with water",units="/",missing_value=missing_value,ncid=ncid)


    return
    
  end subroutine soil_nc_write
   
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s o i l _ a v e
  ! Purpose  :  Average (or sum) the paladyn fields as needed
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine soil_ave(d,ave)
    
    implicit none
    
    type(soil_out) :: d(:), ave

    
    integer :: k, n
    real(wp) :: div
    
    n = size(d)
    div = real(n,wp)

    ! Set all values to zero
    ave%tsoil      = 0._wp
    ave%tice       = 0._wp
    ave%tsublake      = 0._wp
    ave%thetaw     = 0._wp
    ave%thetai     = 0._wp
    ave%theta      = 0._wp
    ave%thetas     = 0._wp
    ave%fthetas    = 0._wp


    ! Loop over the time indices to sum up and average (if necessary)
    do k = 1, n
      ave%tsoil      = ave%tsoil         + d(k)%tsoil     / div
      ave%tice       = ave%tice         + d(k)%tice     / div
      ave%tsublake   = ave%tsublake      + d(k)%tsublake     / div
      ave%thetaw     = ave%thetaw        + d(k)%thetaw    / div
      ave%thetai     = ave%thetai        + d(k)%thetai    / div
      ave%theta      = ave%theta         + d(k)%theta     / div
      ave%thetas     = ave%thetas        + d(k)%thetas    / div
      ave%fthetas    = ave%fthetas       + d(k)%fthetas   / div
    end do
    
   return
    
  end subroutine soil_ave
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s o i l _ p a r _ n c
  ! Purpose  :  Initialize netcdf output for paladyn
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine soil_par_nc(fnm)
    
    implicit none
    
    character (len=*) :: fnm
    integer :: ncid
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm, dim_time, x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE., ncid=ncid)
    call nc_write_dim(fnm, dim_month, x=1._wp, dx=1._wp, nx=13, axis="e", &
    units="months", ncid=ncid)
    call nc_write_dim(fnm, dim_depth,  x=z, units="m", axis="z", ncid=ncid)
    call nc_write_dim(fnm, dim_lat,x=lat,axis="y",ncid=ncid)
    call nc_write_dim(fnm, dim_lon,x=lon,axis="x",ncid=ncid)
    call nc_write_dim(fnm, dim_ncarb,x=1._wp,dx=1._wp,nx=ncarb,units="n/a",ncid=ncid)
    call nc_close(ncid)

    return
  
  end subroutine soil_par_nc
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s o i l _ p a r _ n c _ w r i t e
  ! Purpose  :  Output timestep of netcdf output for paladyn
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine soil_par_nc_write(fnm,ncid,vars,ndat,nout)
    
    implicit none
    
    type(soil_par_out) :: vars
 
    character (len=*) :: fnm
    integer :: ndat, nout, ncid


    call nc_write(fnm,"lambda_if",   sngl(vars%lambda_if),dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nl+1,1,1],long_name="heat conductivity",units="W/m/K",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"cap_if",      sngl(vars%cap_if),   dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nl+1,1,1],long_name="heat capacity",units="J/m3/K",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"lambda_i",   sngl(vars%lambda_i),dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nl+1,1,1],long_name="heat conductivity of ice",units="W/m/K",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"cap_i",      sngl(vars%cap_i),   dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nl+1,1,1],long_name="heat capacity of ice",units="J/m3/K",missing_value=missing_value,ncid=ncid)

    call nc_write(fnm,"kappa",    sngl(vars%kappa),dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,2,ndat,nout],count=[nx,ny,nl,1,1],long_name="hydraulic conductivity",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"ksat",    sngl(vars%ksat),dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,2,ndat,nout],count=[nx,ny,nl,1,1],long_name="hydraulic conductivity",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"psisat",    sngl(vars%psisat),dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,2,ndat,nout],count=[nx,ny,nl,1,1],long_name="hydraulic conductivity",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"b",    sngl(vars%bi),dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,2,ndat,nout],count=[nx,ny,nl,1,1],long_name="hydraulic conductivity",units="kg/m2/day",missing_value=missing_value,ncid=ncid)

    call nc_write(fnm,"theta_sat",   sngl(vars%theta_sat),dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,2,ndat,nout],count=[nx,ny,nl,1,1],long_name="soil porosity",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"theta_field",   sngl(vars%theta_field),dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,2,ndat,nout],count=[nx,ny,nl,1,1],long_name="soil moisture at field capacity",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"theta_wilt",   sngl(vars%theta_wilt),dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,2,ndat,nout],count=[nx,ny,nl,1,1],long_name="soil moisture at wilting point",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fsoc",   sngl(vars%fsoc),  dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,2,ndat,nout],count=[nx,ny,nl,1,1],long_name="soil organic carbon fraction",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"psi",    sngl(vars%psi),  dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,2,ndat,nout],count=[nx,ny,nl,1,1],long_name="soil matric potential",units="m",missing_value=missing_value,ncid=ncid)

    call nc_write(fnm,"ftemp",  sngl(vars%ftemp),  dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nl,1,1],long_name="temperature limitation factor for soil respiration",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fmoist",  sngl(vars%fmoist), dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nl,1,1],long_name="moisture limitation factor for soil respiration",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"fdepth",  sngl(vars%fdepth),  dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nl,1,1],long_name="depth factor for soil respiration",units="/",missing_value=missing_value,ncid=ncid)

    return
    
  end subroutine soil_par_nc_write
   
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s o i l _ p a r _ a v e
  ! Purpose  :  Average (or sum) the paladyn fields as needed
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine soil_par_ave(d,ave)
    
    implicit none
    
    type(soil_par_out) :: d(:), ave

    
    integer :: k, n
    real(wp) :: div
    
    n = size(d)
    div = real(n,wp)

    ! Set all values to zero
    ave%lambda_if  = 0._wp
    ave%cap_if     = 0._wp
    ave%lambda_i  = 0._wp
    ave%cap_i     = 0._wp
    ave%kappa   = 0._wp
    ave%ksat   = 0._wp
    ave%psisat   = 0._wp
    ave%bi   = 0._wp
    ave%theta_sat  = 0._wp
    ave%theta_field  = 0._wp
    ave%theta_wilt  = 0._wp
    ave%fsoc    = 0._wp
    ave%psi     = 0._wp
    ave%ftemp   = 0._wp
    ave%fmoist  = 0._wp
    ave%fdepth  = 0._wp

    ! Loop over the time indices to sum up and average (if necessary)
    do k = 1, n
      ave%lambda_if = ave%lambda_if  + d(k)%lambda_if  / div
      ave%cap_if    = ave%cap_if     + d(k)%cap_if     / div
      ave%lambda_i = ave%lambda_i  + d(k)%lambda_i  / div
      ave%cap_i    = ave%cap_i     + d(k)%cap_i     / div
      ave%kappa  = ave%kappa   + d(k)%kappa   / div
      ave%ksat  = ave%ksat   + d(k)%ksat   / div
      ave%psisat  = ave%psisat   + d(k)%psisat   / div
      ave%bi  = ave%bi   + d(k)%bi   / div
      ave%theta_sat = ave%theta_sat  + d(k)%theta_sat  / div
      ave%theta_field = ave%theta_field  + d(k)%theta_field  / div
      ave%theta_wilt = ave%theta_wilt  + d(k)%theta_wilt  / div
      ave%ftemp  = ave%ftemp   + d(k)%ftemp   / div
      ave%fmoist = ave%fmoist  + d(k)%fmoist  / div
      ave%fdepth = ave%fdepth  + d(k)%fdepth  / div
      ave%fsoc   = ave%fsoc    + d(k)%fsoc    / div
      ave%psi    = ave%psi     + d(k)%psi     / div
    end do
    
   return
    
  end subroutine soil_par_ave
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  l a k e _ n c
  ! Purpose  :  Initialize netcdf output for paladyn
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lake_nc(fnm)
    
    implicit none
    
    character (len=*) :: fnm
    integer :: ncid
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm, dim_time, x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE., ncid=ncid)
    call nc_write_dim(fnm, dim_month, x=1._wp, dx=1._wp, nx=13, axis="e", &
    units="months", ncid=ncid)
    call nc_write_dim(fnm, dim_depth, x=z_l(1:nl_l), units="m", axis="z", ncid=ncid)
    call nc_write_dim(fnm, dim_lat, x=lat, axis="y", ncid=ncid)
    call nc_write_dim(fnm, dim_lon, x=lon, axis="x", ncid=ncid)
    call nc_close(ncid)

    return
  
  end subroutine lake_nc
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  l a k e _ n c _ w r i t e
  ! Purpose  :  Output timestep of netcdf output for paladyn
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lake_nc_write(fnm,ncid,vars,ndat,nout)
    
    implicit none
    
    type(lake_out) :: vars
 
    character (len=*) :: fnm
    integer :: ndat, nout, ncid
 

    call nc_write(fnm,"t_lake",    sngl(vars%t_lake), dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nl_l,1,1],long_name="lake temperature",units="K",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"lambda_lake",    sngl(vars%lambda_lake), dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nl_l,1,1],long_name="lake thermal conductivity",units="W/m/K",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"f_i_lake",     sngl(vars%f_i_lake),  dims=[dim_lon,dim_lat,dim_depth,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[nx,ny,nl_l,1,1],long_name="frozen lake water fraction",units="/",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"f_lake_ice",    sngl(vars%f_lake_ice), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="ice fraction over lake",units="m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"lake_bal", sngl(vars%lakebal), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="lake surface water balance",units="kg/m2/day",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"h_lake_conv",    sngl(vars%h_lake_conv), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="convection depth in lake",units="m",missing_value=missing_value,ncid=ncid)
    call nc_write(fnm,"h_lake_mix",    sngl(vars%h_lake_mix), dims=[dim_lon,dim_lat,dim_month,dim_time],start=[1,1,ndat,nout],count=[nx,ny,1,1],long_name="mixed layer depth in lake (without convection)",units="m",missing_value=missing_value,ncid=ncid)
    if( ndat .eq. nmon_year+1 ) then
     call nc_write(fnm,"h_lake",     sngl(vars%h_lake),  dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[nx,ny,1],long_name="lake depth",units="m",missing_value=missing_value,ncid=ncid)
    endif


    return
    
  end subroutine lake_nc_write
   
  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  l a k e _ a v e
  ! Purpose  :  Average (or sum) the paladyn fields as needed
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lake_ave(d,ave)
    
    implicit none
    
    type(lake_out) :: d(:), ave

    
    integer :: k, n
    real(wp) :: div
    
    n = size(d)
    div = real(n,wp)

    ! Set all values to zero
    ave%t_lake     = 0._wp
    ave%lambda_lake     = 0._wp
    ave%f_i_lake   = 0._wp
    ave%f_lake_ice   = 0._wp
    ave%lakebal = 0._wp
    ave%h_lake_conv    = 0._wp
    ave%h_lake_mix     = 0._wp

    ! Loop over the time indices to sum up and average (if necessary)
    do k = 1, n
      ave%t_lake      = ave%t_lake         + d(k)%t_lake    / div
      ave%lambda_lake      = ave%lambda_lake         + d(k)%lambda_lake    / div
      ave%f_i_lake    = ave%f_i_lake       + d(k)%f_i_lake  / div
      ave%f_lake_ice    = ave%f_lake_ice       + d(k)%f_lake_ice  / div
      ave%lakebal = ave%lakebal    + d(k)%lakebal      / div
      ave%h_lake_conv     = ave%h_lake_conv        + d(k)%h_lake_conv   / div
      ave%h_lake_mix      = ave%h_lake_mix         + d(k)%h_lake_mix    / div
    end do
    
   return
    
  end subroutine lake_ave
  

end module lnd_out


