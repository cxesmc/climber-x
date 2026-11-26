!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : a t m _ p a r a m s
!
!  Purpose : definition of atmospheric model parameters
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Andrey Ganopolski and Matteo Willeit
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
module atm_params

  use precision, only : wp, sp, dp
  use nml, only : nml_read
  use constants, only : Le, Ls, g, Rd, T0
  use control, only : out_dir, flag_dust, flag_co2
  use timer, only : dt_atm, sec_day

  implicit none

  real(wp) :: tstep
  integer :: nstep_fast

  real(wp) :: fcormin
  integer :: i_fcorg

  real(wp) :: f_ice_pow
  real(wp) :: r_scat
  real(wp) :: a1_w 
  real(wp) :: a2_w 
  real(wp) :: b1_w
  real(wp) :: b2_w
  real(wp) :: c_itf_c
  real(wp) :: c_itf_cc
  integer :: i_swr_ray
  real(wp) :: tau_ray
  integer :: i_swr_gas
  real(wp) :: a_gas_vu
  real(wp) :: b_gas_vu
  real(wp) :: a_gas_ir
  real(wp) :: b_gas_ir
  real(wp) :: n_wp
  real(wp) :: p_cosz_sw

  logical :: l_so4_de
  logical :: l_so4_ie
  real(wp) :: beta_so4
  real(wp) :: sigma_so4
  real(wp) :: r_so4
  real(wp) :: N_so4_nat

  integer :: i_mmc_fer
  integer :: i_mmc_fer_shape    !! shape of the Ferrel lower branch: 0 = half-sine, 2 = local psi_F = aferz(j) (requires i_mmc_fer=2)
  real(wp) :: c_mmc_fer_q      !! exponent of the edge taper |sin(ff)|^q of the Ferrel branch, used only if i_mmc_fer_shape=2
  integer :: i_mmc_had
  real(wp) :: c_mmc_had
  real(wp) :: c_mmc_q
  real(wp) :: q_mmc_ref
  real(wp) :: c_mmc_fer
  real(wp) :: c_mmc_pol
  real(wp) :: c_mmc_1
  integer  :: i_mmc_geo
  real(wp) :: c_mmc_edge
  real(wp) :: c_mmc_2
  integer  :: i_mmc_itcz
  real(wp) :: c_mmc_efe
  real(wp) :: c_mmc_dt0
  real(wp) :: c_mmc_dt1
  real(wp) :: hdv_mmc_ref
  real(wp) :: c_mmc_mt        !! weight of the mountain torque in the lower-branch momentum budget of zslp (0 = off)
  integer  :: nsmooth_mtz     !! passes of the meridional 1-2-1 smoothing of the zonal mean mountain torque (0 = off)
  real(wp) :: pbl             !! relative pressure of the PBL top (p/ps), uniform in latitude

  real(wp) :: c_uter_pol
  real(wp) :: c_uter_eq
  integer :: i_mass_com
  real(wp) :: dps_com_topo
  real(wp) :: dp_com_topo

  real(wp) :: c_slp_1
  real(wp) :: c_slp_2
  real(wp) :: c_aslp_eq_hb
  real(wp) :: c_aslp_eq_lat0
  real(wp) :: c_slp_clim_f
  real(wp) :: c_slp_clim_tau
  integer  :: n_slp_clim_harm
  logical :: l_aslp_temp_adv
  real(wp) :: c_aslp_temp_tau
  real(wp) :: c_aslp_temp_dT
  real(wp) :: c_aslp_temp_wtau
  logical :: l_aslp_topo
  real(wp) :: c_aslp_sw_amp
  real(wp) :: c_aslp_sw_tau
  real(wp) :: c_aslp_sw_h
  real(wp) :: c_aslp_sw_umin
  real(wp) :: c_aslp_sw_fr
  real(wp) :: c_aslp_sw_nbv
  real(wp) :: c_aslp_sw_lateq
  integer :: n_aslp_sw

  real(wp) :: zmax
  real(wp) :: dpc
  integer  :: i_dpc
  real(wp) :: dp_com
  real(wp) :: dps_com
  real(wp) :: ptopdyn

  real(wp) :: hcld_base
  integer :: i_lw_cld
  real(wp) :: c_lw_clot
  integer :: i_lw_cld_ice
  real(wp) :: c_lw_clot_ice
  real(wp) :: T_lw_ice
  real(wp) :: ak_co2
  real(wp) :: beta_co2
  real(wp) :: ak_wv
  real(wp) :: a_vap
  real(wp) :: beta_vap
  real(wp) :: a2_vap
  real(wp) :: beta2_vap
  real(wp) :: a3_vap
  integer :: i_lwr_band
  real(wp) :: nu_win1, nu_win2, nu_co21, nu_co22
  real(wp) :: k_cont, c_cont_T, k_wvw, ak_wvw, f_wv_co2
  real(wp) :: ak_co2_3b, beta_co2_3b, a0_co2_3b, a1_co2_3b
  real(wp) :: ak_wv_3b, a_vap_3b, beta_vap_3b, a2_vap_3b, beta2_vap_3b, a3_vap_3b
  real(wp) :: ak_o3_3b, a_o3_3b, beta_o3_3b
  real(wp) :: k_co2w, k_co2r, b_co2w, c_co2_hi
  real(wp) :: z_top_lw
  logical :: l_o3
  real(wp) :: ecs_scale
  real(wp) :: ecs_scale_dT

  real(wp) :: c_gam_1   
  real(wp) :: c_gam_2
  real(wp) :: c_gam_3
  real(wp) :: gams_min
  real(wp) :: dgams         !! K/m, excess of the unstable surface-layer lapse rate over the free-tropospheric gamb, over ice-free ocean and land
  real(wp) :: sh_gams_pos
  real(wp) :: sh_gams_0
  real(wp) :: sh_gams_neg
  real(wp) :: hgams
  real(wp) :: hgamt
  integer :: nsmooth_gam
  real(wp) :: c_gam_tau
  real(wp) :: c_hrm_tau


  real(wp) :: c_hrs_1   
  real(wp) :: c_hrs_2
  real(wp) :: c_hrs_3
  real(wp) :: c_hrs_4
  real(wp) :: c_hrs_5
  real(wp) :: c_hrs_6

  real(wp) :: c_clot_1   
  real(wp) :: c_clot_2
  real(wp) :: c_clot_3   
  real(wp) :: c_clot_4
  real(wp) :: c_clot_5

  real(wp) :: c_syn_1
  real(wp) :: c_syn_2
  real(wp) :: c_syn_3 
  real(wp) :: c_syn_4
  integer  :: i_syn_diss   !! synoptic energy dissipation: 1 = (c_syn_3+c_syn_4*cda)*EKE^1.5, 2 = linear, (c_syn_3+c_syn_4*cda)*sqrt(syn_e0)*EKE
  real(wp) :: syn_e0       !! m2/s2, reference synoptic energy of the linear dissipation (i_syn_diss=2)
  real(wp) :: dz_syn_min   !! m, minimum depth of the layer over which the EKE production shear and stability are taken
  real(wp) :: c_eke_slope  !! sloping-boundary reduction of EKE production: x max(0,1-c_eke_slope*delta), delta = terrain slope / isentrope slope (0 = off)
  real(wp) :: L_eke_slope  !! m, e-folding scale of the Gaussian terrain smoothing for the EKE slope term (eddy scale; 0 = local grid slope)
  real(wp) :: c_syn_5      !! s, macro-diffusivity factor for the synoptic energy, K = c_syn_5*EKE
  real(wp) :: c_syn_6
  real(wp) :: c_syn_7
  real(wp) :: windmin
  real(wp) :: synsurmin
  real(wp) :: tau_fac
  real(wp) :: c_tau_syn
  logical  :: l_tau_stab
  real(wp) :: f_tau_Ri_unstab
  real(wp) :: f_tau_Ri_stab
  real(wp) :: c_diff       !! s, macro-diffusivity factor, K = c_diff*EKE
  real(wp) :: h_diff       !! m, e-folding height of the macro-diffusivity, K(z) = K*exp(-z/h_diff) (Branscome 1983); <=0 = vertically uniform
  integer  :: i_hdiff_ref  !! reference height for h_diff: 1 = local (face) surface, 2 = sea level
  real(wp) :: f_diff_wtr   !! factor on the EKE macro-diffusivity for water vapour (not dust/CO2), K_q = f_diff_wtr*K; 1 = same K as DSE
  integer  :: i_diff_prof  !! vertical profile of the macro-diffusivity: 0 = exp(-z/h_diff) up to 0.2*p0, 1 = linear from 1 at the surface to 0 at the tropopause, 2 = as 1 but linear from 1 at sea level, 3 = uniform up to the tropopause
  real(wp) :: diff_stab   !! m2/s, uniform additional ZONAL macro-diffusivity (all tracers) for numerical stability of i_adv=2; 0 = off
  integer :: i_adv         !! advection scheme: 1 = first order upstream (donor cell), 2 = flux corrected transport (Zalesak 1979)
  real(wp) :: c_filt_conv  !! strength [0..1] of the conservative high-order 2dx (Shapiro) filter applied to the convergences (convdse/convwtr) to suppress the grid-scale checkerboard; 0 = off
  integer :: nord_filt_conv !! order (selectivity) of the convergence 2dx filter; higher = more selective (spares resolved scales)
  integer :: i_filt_conv_polar !! zonal Fourier filter of the convergences poleward of lat_filt_conv_polar (Arakawa-Lamb taper on the physical wavelength): 0 = off, 1 = on
  real(wp) :: lat_filt_conv_polar !! deg, cutoff latitude of the polar convergence filter

  logical  :: l_kata_wind
  real(wp) :: h_kata
  real(wp) :: c_kata
  logical  :: l_kata_upw

  real(wp) :: cd0_ocn
  real(wp) :: cd0_sic
  integer :: i_zoro
  integer :: i_acbar
  real(wp) :: acbar_max
  real(wp) :: acbar_scale
  real(wp) :: c_acbar_0
  real(wp) :: c_acbar_f
  real(wp) :: c_acbar_p
  real(wp) :: c_acbar_wind
  integer  :: i_usur_eps
  real(wp) :: h_usur_eps

  real(wp) :: c_cld_1
  real(wp) :: c_cld_2
  real(wp) :: c_cld_3
  real(wp) :: c_cld_4
  real(wp) :: c_cld_5
  real(wp) :: c_cld_6
  real(wp) :: c_cld_60_ocn
  real(wp) :: c_cld_60_lnd
  real(wp) :: c_cld_7
  real(wp) :: c_cld_8
  real(wp) :: c_cld_9
  real(wp) :: c_cld_10
  real(wp) :: c_cld_5_sc
  real(wp) :: cld_max
  integer :: i_cld_low
  real(wp) :: q_cld_low
  integer :: nsmooth_cld

  real(wp) :: c_hcld_low
  real(wp) :: c_hcld_dz_min
  real(wp) :: c_hcld_1
  real(wp) :: c_hcld_2
  real(wp) :: c_hcld_3
  real(wp) :: c_hcld_4

  real(wp) :: c_weff
  real(wp) :: c_woro

  integer :: i_rbstr
  real(wp) :: c_trop_1
  real(wp) :: c_trop_2
  real(wp) :: c_trop_3


  real(wp) :: hpbl

  real(wp) :: rh_max_ocn
  real(wp) :: rh_max_lnd
  real(wp) :: rh_strat

  real(wp) :: z_q2
  real(wp) :: z_q2_ocn

  logical :: l_dust
  logical :: l_dust_rad
  real(wp) :: c_dhs_1
  real(wp) :: c_dhs_2
  real(wp) :: c_dust_dry
  real(wp) :: c_dust_wet
  real(wp) :: c_dust_mec

  logical :: l_co2d

  real(wp) :: ars_ot
  real(wp) :: ars_im
  real(wp) :: h_ars

  integer :: nsmooth_cda
  integer :: nsmooth_weff
  integer :: n_aslp_shapiro
  real(wp) :: aslp_polar_lat
  integer :: nsmooth_acbar

  logical :: l_write_timer

  logical :: l_daily_output
  logical :: l_output_flx3d
  logical :: l_output_extended

  real(wp) :: tam_init = 30.0_wp   !! initial atm peak temperature at equator

  real(wp), parameter :: cle = Le         !! J/kg, latent heat of evaporation
  real(wp), parameter :: cls = Ls         !! J/kg, latent heat of melting + evaporation (sublimation)
  real(wp), parameter :: cp = 1000._wp       !! J/kg/K, specific heat of air at constant pressure
  real(wp), parameter :: cv = 715._wp       !! J/kg/K, specific heat of air at constant volume
  real(wp), parameter :: gad = g/cp       !! K/m, adiabatic lapse rate
  real(wp), parameter :: eps = 0.001_wp      !! small number
  real(wp), parameter :: atm_mass = 5.12e18_wp         !! kg, total mass of atmosphere
  real(wp), parameter :: hatm = Rd*T0/g      !! m, atmospheric height scale 
  logical :: l_p0_var   !! variable mean sea level pressure?
  real(wp) :: p0      !! Pa, average sea level pressure 
  real(wp) :: ps0     !! Pa, average surface pressure 
  real(wp) :: amas    !! kg/m2, average mass of atmospheric column
  real(wp) :: ra      !! kg/m3, air density at p0

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  a t m _ p a r a m s _ i n i t
  !   Purpose    :  reading/initialisation of atmosphere parameters
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine atm_params_init

    implicit none

    ! read namelist parameters
    call atm_par_load(trim(out_dir)//"/atm_par.nml")

    ! time step
    tstep = dt_atm/nstep_fast   ! s

    return

  end subroutine atm_params_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  a t m _ p a r _ l o a d 
  !   Purpose    :  read atmosphere parameters namelist
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine atm_par_load(filename)

    implicit none

    character(len=*), intent(in) :: filename


    ! Read parameters from file
    write(*,*) "atmosphere parameters ==========="
    call nml_read(filename,"atm_par","nstep_fast",nstep_fast)
    call nml_read(filename,"atm_par","fcormin",fcormin)
    call nml_read(filename,"atm_par","i_fcorg",i_fcorg)
    call nml_read(filename,"atm_par","f_ice_pow",f_ice_pow)
    call nml_read(filename,"atm_par","r_scat",r_scat)
    call nml_read(filename,"atm_par","a1_w",a1_w)
    a2_w = 1._wp-a1_w
    call nml_read(filename,"atm_par","b1_w",b1_w)
    call nml_read(filename,"atm_par","b2_w",b2_w)
    call nml_read(filename,"atm_par","c_itf_c",c_itf_c)
    call nml_read(filename,"atm_par","c_itf_cc",c_itf_cc)
    call nml_read(filename,"atm_par","i_swr_ray",i_swr_ray)
    call nml_read(filename,"atm_par","tau_ray",tau_ray)
    call nml_read(filename,"atm_par","i_swr_gas",i_swr_gas)
    call nml_read(filename,"atm_par","a_gas_vu",a_gas_vu)
    call nml_read(filename,"atm_par","b_gas_vu",b_gas_vu)
    call nml_read(filename,"atm_par","a_gas_ir",a_gas_ir)
    call nml_read(filename,"atm_par","b_gas_ir",b_gas_ir)
    call nml_read(filename,"atm_par","n_wp",n_wp)
    call nml_read(filename,"atm_par","p_cosz_sw",p_cosz_sw)
    call nml_read(filename,"atm_par","i_mmc_fer",i_mmc_fer)
    call nml_read(filename,"atm_par","i_mmc_fer_shape",i_mmc_fer_shape)
    call nml_read(filename,"atm_par","c_mmc_fer_q",c_mmc_fer_q)
    if (i_mmc_fer_shape.eq.2 .and. i_mmc_fer.ne.2) stop 'ERROR atm_params: i_mmc_fer_shape=2 requires i_mmc_fer=2'
    call nml_read(filename,"atm_par","i_mmc_had",i_mmc_had)
    call nml_read(filename,"atm_par","c_mmc_had",c_mmc_had)
    call nml_read(filename,"atm_par","c_mmc_q",c_mmc_q)
    call nml_read(filename,"atm_par","q_mmc_ref",q_mmc_ref)
    call nml_read(filename,"atm_par","c_mmc_fer",c_mmc_fer)
    call nml_read(filename,"atm_par","c_mmc_pol",c_mmc_pol)
    call nml_read(filename,"atm_par","c_mmc_1",c_mmc_1)
    call nml_read(filename,"atm_par","i_mmc_geo",i_mmc_geo)
    call nml_read(filename,"atm_par","c_mmc_edge",c_mmc_edge)
    call nml_read(filename,"atm_par","c_mmc_2",c_mmc_2)
    call nml_read(filename,"atm_par","i_mmc_itcz",i_mmc_itcz)
    call nml_read(filename,"atm_par","c_mmc_efe",c_mmc_efe)
    call nml_read(filename,"atm_par","c_mmc_dt0",c_mmc_dt0)
    call nml_read(filename,"atm_par","c_mmc_dt1",c_mmc_dt1)
    call nml_read(filename,"atm_par","hdv_mmc_ref",hdv_mmc_ref)
    call nml_read(filename,"atm_par","c_mmc_mt",c_mmc_mt)
    nsmooth_mtz = 1
    call nml_read(filename,"atm_par","nsmooth_mtz",nsmooth_mtz)
    call nml_read(filename,"atm_par","c_uter_pol",c_uter_pol)
    call nml_read(filename,"atm_par","c_uter_eq",c_uter_eq)
    call nml_read(filename,"atm_par","i_mass_com",i_mass_com)
    call nml_read(filename,"atm_par","dps_com_topo",dps_com_topo)
    call nml_read(filename,"atm_par","dp_com_topo",dp_com_topo)
    call nml_read(filename,"atm_par","l_p0_var",l_p0_var)
    call nml_read(filename,"atm_par","p0",p0)
    call nml_read(filename,"atm_par","c_slp_1",c_slp_1)
    call nml_read(filename,"atm_par","c_slp_2",c_slp_2)
    call nml_read(filename,"atm_par","c_aslp_eq_hb",c_aslp_eq_hb)
    call nml_read(filename,"atm_par","c_aslp_eq_lat0",c_aslp_eq_lat0)
    c_slp_clim_f = 1._wp
    call nml_read(filename,"atm_par","c_slp_clim_f",c_slp_clim_f)
    c_slp_clim_tau = 5._wp
    call nml_read(filename,"atm_par","c_slp_clim_tau",c_slp_clim_tau)
    n_slp_clim_harm = 4
    call nml_read(filename,"atm_par","n_slp_clim_harm",n_slp_clim_harm)
    call nml_read(filename,"atm_par","l_aslp_temp_adv",l_aslp_temp_adv)
    call nml_read(filename,"atm_par","c_aslp_temp_tau",c_aslp_temp_tau)
    c_aslp_temp_tau = c_aslp_temp_tau * sec_day
    c_aslp_temp_dT = 0._wp
    call nml_read(filename,"atm_par","c_aslp_temp_dT",c_aslp_temp_dT)
    c_aslp_temp_wtau = 20._wp
    call nml_read(filename,"atm_par","c_aslp_temp_wtau",c_aslp_temp_wtau)
    call nml_read(filename,"atm_par","l_aslp_topo",l_aslp_topo)
    call nml_read(filename,"atm_par","c_aslp_sw_amp",c_aslp_sw_amp)
    call nml_read(filename,"atm_par","c_aslp_sw_tau",c_aslp_sw_tau)
    c_aslp_sw_tau = c_aslp_sw_tau * sec_day
    call nml_read(filename,"atm_par","c_aslp_sw_h",c_aslp_sw_h)
    call nml_read(filename,"atm_par","c_aslp_sw_umin",c_aslp_sw_umin)
    call nml_read(filename,"atm_par","c_aslp_sw_fr",c_aslp_sw_fr)
    call nml_read(filename,"atm_par","c_aslp_sw_nbv",c_aslp_sw_nbv)
    call nml_read(filename,"atm_par","c_aslp_sw_lateq",c_aslp_sw_lateq)
    call nml_read(filename,"atm_par","n_aslp_sw",n_aslp_sw)
    call nml_read(filename,"atm_par","zmax",zmax)
    call nml_read(filename,"atm_par","dpc",dpc)
    call nml_read(filename,"atm_par","i_dpc",i_dpc)
    call nml_read(filename,"atm_par","dp_com",dp_com)
    call nml_read(filename,"atm_par","dps_com",dps_com)
    call nml_read(filename,"atm_par","ptopdyn",ptopdyn)
    call nml_read(filename,"atm_par","pbl",pbl)
    call nml_read(filename,"atm_par","c_gam_1",c_gam_1)
    call nml_read(filename,"atm_par","c_gam_2",c_gam_2)
    call nml_read(filename,"atm_par","c_gam_3",c_gam_3)
    call nml_read(filename,"atm_par","gams_min",gams_min)
    call nml_read(filename,"atm_par","dgams",dgams)
    call nml_read(filename,"atm_par","sh_gams_pos",sh_gams_pos)
    call nml_read(filename,"atm_par","sh_gams_0",sh_gams_0)
    call nml_read(filename,"atm_par","sh_gams_neg",sh_gams_neg)
    if (sh_gams_neg.le.0._wp) sh_gams_neg = sh_gams_pos
    call nml_read(filename,"atm_par","hgams",hgams)
    call nml_read(filename,"atm_par","hgamt",hgamt)
    call nml_read(filename,"atm_par","c_gam_tau",c_gam_tau)
    call nml_read(filename,"atm_par","c_hrm_tau",c_hrm_tau)
    call nml_read(filename,"atm_par","nsmooth_gam",nsmooth_gam)
    call nml_read(filename,"atm_par","nsmooth_cld",nsmooth_cld)
    call nml_read(filename,"atm_par","hcld_base",hcld_base)
    call nml_read(filename,"atm_par","i_lw_cld",i_lw_cld)
    call nml_read(filename,"atm_par","c_lw_clot",c_lw_clot)
    call nml_read(filename,"atm_par","i_lw_cld_ice",i_lw_cld_ice)
    call nml_read(filename,"atm_par","c_lw_clot_ice",c_lw_clot_ice)
    call nml_read(filename,"atm_par","T_lw_ice",T_lw_ice)
    call nml_read(filename,"atm_par","ak_co2",ak_co2)
    call nml_read(filename,"atm_par","beta_co2",beta_co2)
    call nml_read(filename,"atm_par","ak_wv",ak_wv)
    call nml_read(filename,"atm_par","a_vap",a_vap)
    call nml_read(filename,"atm_par","beta_vap",beta_vap)
    call nml_read(filename,"atm_par","a2_vap",a2_vap)
    call nml_read(filename,"atm_par","beta2_vap",beta2_vap)
    call nml_read(filename,"atm_par","a3_vap",a3_vap)
    call nml_read(filename,"atm_par","i_lwr_band",i_lwr_band)
    call nml_read(filename,"atm_par","nu_win1",nu_win1)
    call nml_read(filename,"atm_par","nu_win2",nu_win2)
    call nml_read(filename,"atm_par","nu_co21",nu_co21)
    call nml_read(filename,"atm_par","nu_co22",nu_co22)
    call nml_read(filename,"atm_par","k_cont",k_cont)
    call nml_read(filename,"atm_par","c_cont_T",c_cont_T)
    call nml_read(filename,"atm_par","k_wvw",k_wvw)
    call nml_read(filename,"atm_par","ak_wvw",ak_wvw)
    call nml_read(filename,"atm_par","f_wv_co2",f_wv_co2)
    call nml_read(filename,"atm_par","ak_co2_3b",ak_co2_3b)
    call nml_read(filename,"atm_par","beta_co2_3b",beta_co2_3b)
    call nml_read(filename,"atm_par","a0_co2_3b",a0_co2_3b)
    call nml_read(filename,"atm_par","a1_co2_3b",a1_co2_3b)
    call nml_read(filename,"atm_par","ak_wv_3b",ak_wv_3b)
    call nml_read(filename,"atm_par","a_vap_3b",a_vap_3b)
    call nml_read(filename,"atm_par","beta_vap_3b",beta_vap_3b)
    call nml_read(filename,"atm_par","a2_vap_3b",a2_vap_3b)
    call nml_read(filename,"atm_par","beta2_vap_3b",beta2_vap_3b)
    call nml_read(filename,"atm_par","a3_vap_3b",a3_vap_3b)
    call nml_read(filename,"atm_par","ak_o3_3b",ak_o3_3b)
    call nml_read(filename,"atm_par","a_o3_3b",a_o3_3b)
    call nml_read(filename,"atm_par","beta_o3_3b",beta_o3_3b)
    call nml_read(filename,"atm_par","k_co2w",k_co2w)
    call nml_read(filename,"atm_par","k_co2r",k_co2r)
    call nml_read(filename,"atm_par","b_co2w",b_co2w)
    call nml_read(filename,"atm_par","c_co2_hi",c_co2_hi)
    call nml_read(filename,"atm_par","z_top_lw",z_top_lw)
    call nml_read(filename,"atm_par","l_o3",l_o3)
    call nml_read(filename,"atm_par","ecs_scale",ecs_scale)
    call nml_read(filename,"atm_par","ecs_scale_dT",ecs_scale_dT)
    call nml_read(filename,"atm_par","c_hrs_1",c_hrs_1)
    call nml_read(filename,"atm_par","c_hrs_2",c_hrs_2)
    call nml_read(filename,"atm_par","c_hrs_3",c_hrs_3)
    call nml_read(filename,"atm_par","c_hrs_4",c_hrs_4)
    call nml_read(filename,"atm_par","c_hrs_5",c_hrs_5)
    call nml_read(filename,"atm_par","c_hrs_6",c_hrs_6)
    call nml_read(filename,"atm_par","c_clot_1",c_clot_1)
    call nml_read(filename,"atm_par","c_clot_2",c_clot_2)
    call nml_read(filename,"atm_par","c_clot_3",c_clot_3)
    call nml_read(filename,"atm_par","c_clot_4",c_clot_4)
    call nml_read(filename,"atm_par","c_clot_5",c_clot_5)
    call nml_read(filename,"atm_par","c_syn_1",c_syn_1)
    call nml_read(filename,"atm_par","c_syn_2",c_syn_2)
    call nml_read(filename,"atm_par","c_syn_3",c_syn_3)
    call nml_read(filename,"atm_par","c_syn_4",c_syn_4)
    call nml_read(filename,"atm_par","i_syn_diss",i_syn_diss)
    call nml_read(filename,"atm_par","c_eke_slope",c_eke_slope)
    call nml_read(filename,"atm_par","L_eke_slope",L_eke_slope)
    call nml_read(filename,"atm_par","dz_syn_min",dz_syn_min)
    call nml_read(filename,"atm_par","syn_e0",syn_e0)
    call nml_read(filename,"atm_par","c_syn_5",c_syn_5)
    call nml_read(filename,"atm_par","c_syn_6",c_syn_6)
    call nml_read(filename,"atm_par","c_syn_7",c_syn_7)
    call nml_read(filename,"atm_par","windmin",windmin)
    call nml_read(filename,"atm_par","synsurmin",synsurmin)
    call nml_read(filename,"atm_par","tau_fac",tau_fac)
    call nml_read(filename,"atm_par","c_tau_syn",c_tau_syn)
    call nml_read(filename,"atm_par","l_tau_stab",l_tau_stab)
    call nml_read(filename,"atm_par","f_tau_Ri_unstab",f_tau_Ri_unstab)
    call nml_read(filename,"atm_par","f_tau_Ri_stab",f_tau_Ri_stab)
    call nml_read(filename,"atm_par","c_diff",c_diff)
    call nml_read(filename,"atm_par","h_diff",h_diff)
    call nml_read(filename,"atm_par","i_hdiff_ref",i_hdiff_ref)
    call nml_read(filename,"atm_par","i_diff_prof",i_diff_prof)
    call nml_read(filename,"atm_par","f_diff_wtr",f_diff_wtr)
    call nml_read(filename,"atm_par","diff_stab",diff_stab)
    call nml_read(filename,"atm_par","i_adv",i_adv)
    call nml_read(filename,"atm_par","c_filt_conv",c_filt_conv)
    call nml_read(filename,"atm_par","nord_filt_conv",nord_filt_conv)
    call nml_read(filename,"atm_par","i_filt_conv_polar",i_filt_conv_polar)
    call nml_read(filename,"atm_par","lat_filt_conv_polar",lat_filt_conv_polar)
    call nml_read(filename,"atm_par","l_kata_wind",l_kata_wind)
    call nml_read(filename,"atm_par","h_kata",h_kata)
    call nml_read(filename,"atm_par","c_kata",c_kata)
    call nml_read(filename,"atm_par","l_kata_upw",l_kata_upw)
    call nml_read(filename,"atm_par","cd0_ocn",cd0_ocn)
    call nml_read(filename,"atm_par","cd0_sic",cd0_sic)
    call nml_read(filename,"atm_par","i_zoro",i_zoro)
    call nml_read(filename,"atm_par","i_acbar",i_acbar)
    call nml_read(filename,"atm_par","acbar_max",acbar_max)
    call nml_read(filename,"atm_par","acbar_scale",acbar_scale)
    call nml_read(filename,"atm_par","c_acbar_0",c_acbar_0)
    call nml_read(filename,"atm_par","c_acbar_f",c_acbar_f)
    call nml_read(filename,"atm_par","c_acbar_p",c_acbar_p)
    call nml_read(filename,"atm_par","c_acbar_wind",c_acbar_wind)
    call nml_read(filename,"atm_par","i_usur_eps",i_usur_eps)
    call nml_read(filename,"atm_par","h_usur_eps",h_usur_eps)
    call nml_read(filename,"atm_par","c_cld_1",c_cld_1)
    call nml_read(filename,"atm_par","c_cld_2",c_cld_2)
    call nml_read(filename,"atm_par","c_cld_3",c_cld_3)
    call nml_read(filename,"atm_par","c_cld_4",c_cld_4)
    call nml_read(filename,"atm_par","c_cld_5",c_cld_5)
    call nml_read(filename,"atm_par","c_cld_6",c_cld_6)
    call nml_read(filename,"atm_par","c_cld_60_ocn",c_cld_60_ocn)
    call nml_read(filename,"atm_par","c_cld_60_lnd",c_cld_60_lnd)
    call nml_read(filename,"atm_par","c_cld_7",c_cld_7)
    call nml_read(filename,"atm_par","c_cld_8",c_cld_8)
    call nml_read(filename,"atm_par","c_cld_9",c_cld_9)
    call nml_read(filename,"atm_par","c_cld_10",c_cld_10)
    call nml_read(filename,"atm_par","c_cld_5_sc",c_cld_5_sc)
    call nml_read(filename,"atm_par","i_cld_low",i_cld_low)
    call nml_read(filename,"atm_par","q_cld_low",q_cld_low)
    call nml_read(filename,"atm_par","c_hcld_low",c_hcld_low)
    call nml_read(filename,"atm_par","c_hcld_dz_min",c_hcld_dz_min)
    call nml_read(filename,"atm_par","c_hcld_1",c_hcld_1)
    call nml_read(filename,"atm_par","c_hcld_2",c_hcld_2)
    call nml_read(filename,"atm_par","c_hcld_3",c_hcld_3)
    call nml_read(filename,"atm_par","c_hcld_4",c_hcld_4)
    call nml_read(filename,"atm_par","cld_max",cld_max)
    call nml_read(filename,"atm_par","c_woro",c_woro)
    call nml_read(filename,"atm_par","c_weff",c_weff)
    call nml_read(filename,"atm_par","i_rbstr",i_rbstr)
    call nml_read(filename,"atm_par","hpbl",hpbl)
    call nml_read(filename,"atm_par","rh_max_ocn",rh_max_ocn)
    call nml_read(filename,"atm_par","rh_max_lnd",rh_max_lnd)
    call nml_read(filename,"atm_par","rh_strat",rh_strat)
    call nml_read(filename,"atm_par","z_q2",z_q2)
    call nml_read(filename,"atm_par","z_q2_ocn",z_q2_ocn)
    call nml_read(filename,"atm_par","c_trop_1",c_trop_1)
    call nml_read(filename,"atm_par","c_trop_2",c_trop_2)
    call nml_read(filename,"atm_par","c_trop_3",c_trop_3)
    call nml_read(filename,"atm_par","l_dust",l_dust)
    call nml_read(filename,"atm_par","l_dust_rad",l_dust_rad)
    if (.not.flag_dust) then
      l_dust = .false.
      l_dust_rad = .false.
    endif
    call nml_read(filename,"atm_par","l_co2d",l_co2d)
    if (.not.flag_co2) then
      l_co2d = .false.
    endif
    call nml_read(filename,"atm_par","c_dhs_1",c_dhs_1)
    call nml_read(filename,"atm_par","c_dhs_2",c_dhs_2)
    call nml_read(filename,"atm_par","c_dust_dry",c_dust_dry)
    call nml_read(filename,"atm_par","c_dust_wet",c_dust_wet)
    call nml_read(filename,"atm_par","c_dust_mec",c_dust_mec)
    call nml_read(filename,"atm_par","l_so4_de",l_so4_de)
    call nml_read(filename,"atm_par","l_so4_ie",l_so4_ie)
    call nml_read(filename,"atm_par","beta_so4",beta_so4)
    call nml_read(filename,"atm_par","sigma_so4",sigma_so4)
    call nml_read(filename,"atm_par","r_so4",r_so4)
    call nml_read(filename,"atm_par","N_so4_nat",N_so4_nat)
    call nml_read(filename,"atm_par","ars_ot",ars_ot)
    call nml_read(filename,"atm_par","ars_im",ars_im)
    call nml_read(filename,"atm_par","h_ars",h_ars)
    call nml_read(filename,"atm_par","nsmooth_cda",nsmooth_cda)
    call nml_read(filename,"atm_par","nsmooth_weff",nsmooth_weff)
    call nml_read(filename,"atm_par","n_aslp_shapiro",n_aslp_shapiro)
    call nml_read(filename,"atm_par","aslp_polar_lat",aslp_polar_lat)
    call nml_read(filename,"atm_par","nsmooth_acbar",nsmooth_acbar)

    call nml_read(filename,"atm_par","l_write_timer",l_write_timer)

    call nml_read(filename,"atm_par","l_daily_output",l_daily_output)
    call nml_read(filename,"atm_par","l_output_flx3d",l_output_flx3d)
    call nml_read(filename,"atm_par","l_output_extended",l_output_extended)

    call nml_read(filename,"atm_par","tam_init",tam_init)
 
    return

  end subroutine atm_par_load

end module atm_params
