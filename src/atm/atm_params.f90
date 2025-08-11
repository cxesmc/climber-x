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
  use timer, only : dt_atm

  implicit none

  real(wp) :: tstep
  integer :: nstep_fast

  real(wp) :: fcormin
  real(wp) :: fcoramin
  real(wp) :: fcorumin

  logical :: l_sct_0
  logical :: l_alb_0
  real(wp) :: f_ice_pow
  real(wp) :: r_scat
  real(wp) :: a1_w 
  real(wp) :: a2_w 
  real(wp) :: b1_w
  real(wp) :: b2_w
  real(wp) :: c_itf_c
  real(wp) :: c_itf_cc

  logical :: l_so4_de
  logical :: l_so4_ie
  real(wp) :: beta_so4
  real(wp) :: sigma_so4
  real(wp) :: r_so4
  real(wp) :: N_so4_nat

  integer :: i_mmc
  real(wp) :: c_mmc_had
  real(wp) :: c_mmc_fer
  real(wp) :: c_mmc_pol
  real(wp) :: c_mmc_z
  real(wp) :: c_mmc_1
  real(wp) :: c_mmc_2
  real(wp) :: c_mmc_3
  real(wp) :: c_mmc_4
  real(wp) :: pblp 
  real(wp) :: pble

  real(wp) :: c_uter_pol
  real(wp) :: c_uter_eq
  logical :: l_mass_com_topo

  real(wp) :: c_slp_1
  real(wp) :: c_slp_2
  real(wp) :: c_slp_3
  real(wp) :: c_slp_4
  real(wp) :: c_slp_5
  logical :: l_aslp_topo
  real(wp) :: c_aslp_topo_1
  real(wp) :: c_aslp_topo_2
  real(wp) :: c_aslp_topo_3
  real(wp) :: c_aslp_topo_4

  real(wp) :: zmax
  real(wp) :: dpc
  real(wp) :: pcmin
  real(wp) :: pcmax
  real(wp) :: ptopdyn

  real(wp) :: hcld_base
  integer :: i_lw_cld
  real(wp) :: c_lw_clot
  real(wp) :: ak_co2
  real(wp) :: beta_co2
  real(wp) :: ak_wv
  real(wp) :: a_vap
  real(wp) :: beta_vap
  real(wp) :: a2_vap
  real(wp) :: beta2_vap
  real(wp) :: a3_vap
  logical :: l_o3
  real(wp) :: ecs_scale
  real(wp) :: ecs_scale_dT

  real(wp) :: c_gam_1   
  real(wp) :: c_gam_2
  real(wp) :: c_gam_3
  real(wp) :: c_gam_4
  real(wp) :: c_gam_5
  real(wp) :: c_gam_6
  real(wp) :: gams_max_lnd
  real(wp) :: gams_min_ocn
  real(wp) :: gams_max_ocn
  real(wp) :: hgams
  real(wp) :: hgamt
  integer :: nsmooth_gam
  real(wp) :: c_gam_rel
  integer :: i_tsl
  integer :: i_tslz
  real(wp) :: c_tsl_gam
  real(wp) :: c_tsl_gam_ice
  real(wp) :: tsl_gams_min_lnd
  real(wp) :: tsl_gams_min_ice

  real(wp) :: c_wrt_1
  real(wp) :: c_wrt_2
  real(wp) :: c_wrt_3   
  real(wp) :: c_wrt_4

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

  real(wp) :: c_syn_1
  real(wp) :: c_syn_2
  real(wp) :: c_syn_3 
  real(wp) :: c_syn_4
  real(wp) :: c_syn_5    
  real(wp) :: c_syn_6
  real(wp) :: c_syn_7
  real(wp) :: c_syn_8
  real(wp) :: windmin 
  real(wp) :: synsurmin
  real(wp) :: c_wind_ele 
  real(wp) :: c_diffx_dse
  real(wp) :: c_diffx_wtr
  real(wp) :: c_diff_dse
  integer :: i_diff_wtr
  real(wp) :: c_diff_wtr
  integer :: i_diff_dst

  integer :: i_kata_wind
  real(wp) :: h_kata

  real(wp) :: cd0_ocn
  real(wp) :: cd0_sic
  integer :: i_zoro
  integer :: i_acbar
  real(wp) :: acbar_max
  real(wp) :: acbar_scale

  real(wp) :: c_cld_1
  real(wp) :: c_cld_2
  real(wp) :: c_cld_3
  real(wp) :: c_cld_4
  real(wp) :: c_cld_5
  real(wp) :: c_cld_55
  real(wp) :: c_cld_6
  real(wp) :: c_cld_7
  real(wp) :: c_cld_8
  logical :: l_cld_low_ice
  real(wp) :: cld_max
  integer :: nsmooth_cld

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

  integer :: i_pbl

  real(wp) :: hpbl

  real(wp) :: rh_max
  real(wp) :: rskin_ocn_min
  real(wp) :: rh_strat

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

  integer :: nsmooth_cda
  integer :: nsmooth_weff
  integer :: nsmooth_aslp
  integer :: nsmooth_aslp_eq
  integer :: nj_eq
  integer :: nsmooth_aslp_topo
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
    call nml_read(filename,"atm_par","fcoramin",fcoramin)
    call nml_read(filename,"atm_par","fcorumin",fcorumin)
    call nml_read(filename,"atm_par","f_ice_pow",f_ice_pow)
    call nml_read(filename,"atm_par","r_scat",r_scat)
    call nml_read(filename,"atm_par","l_sct_0",l_sct_0)
    call nml_read(filename,"atm_par","l_alb_0",l_alb_0)
    call nml_read(filename,"atm_par","a1_w",a1_w)
    a2_w = 1._wp-a1_w
    call nml_read(filename,"atm_par","b1_w",b1_w)
    call nml_read(filename,"atm_par","b2_w",b2_w)
    call nml_read(filename,"atm_par","c_itf_c",c_itf_c)
    call nml_read(filename,"atm_par","c_itf_cc",c_itf_cc)
    call nml_read(filename,"atm_par","i_mmc",i_mmc)
    call nml_read(filename,"atm_par","c_mmc_had",c_mmc_had)
    call nml_read(filename,"atm_par","c_mmc_fer",c_mmc_fer)
    call nml_read(filename,"atm_par","c_mmc_pol",c_mmc_pol)
    call nml_read(filename,"atm_par","c_mmc_z",c_mmc_z)
    call nml_read(filename,"atm_par","c_mmc_1",c_mmc_1)
    call nml_read(filename,"atm_par","c_mmc_2",c_mmc_2)
    call nml_read(filename,"atm_par","c_mmc_3",c_mmc_3)
    call nml_read(filename,"atm_par","c_mmc_4",c_mmc_4)
    call nml_read(filename,"atm_par","c_uter_pol",c_uter_pol)
    call nml_read(filename,"atm_par","c_uter_eq",c_uter_eq)
    call nml_read(filename,"atm_par","l_mass_com_topo",l_mass_com_topo)
    call nml_read(filename,"atm_par","l_p0_var",l_p0_var)
    call nml_read(filename,"atm_par","p0",p0)
    call nml_read(filename,"atm_par","c_slp_1",c_slp_1)
    call nml_read(filename,"atm_par","c_slp_2",c_slp_2)
    call nml_read(filename,"atm_par","c_slp_3",c_slp_3)
    call nml_read(filename,"atm_par","c_slp_4",c_slp_4)
    call nml_read(filename,"atm_par","c_slp_5",c_slp_5)
    call nml_read(filename,"atm_par","l_aslp_topo",l_aslp_topo)
    call nml_read(filename,"atm_par","c_aslp_topo_1",c_aslp_topo_1)
    call nml_read(filename,"atm_par","c_aslp_topo_2",c_aslp_topo_2)
    call nml_read(filename,"atm_par","c_aslp_topo_3",c_aslp_topo_3)
    call nml_read(filename,"atm_par","c_aslp_topo_4",c_aslp_topo_4)
    call nml_read(filename,"atm_par","zmax",zmax)
    call nml_read(filename,"atm_par","dpc",dpc)
    call nml_read(filename,"atm_par","pcmin",pcmin)
    call nml_read(filename,"atm_par","pcmax",pcmax)
    call nml_read(filename,"atm_par","ptopdyn",ptopdyn)
    call nml_read(filename,"atm_par","pblp",pblp)
    call nml_read(filename,"atm_par","pble",pble)
    call nml_read(filename,"atm_par","c_gam_1",c_gam_1)
    call nml_read(filename,"atm_par","c_gam_2",c_gam_2)
    call nml_read(filename,"atm_par","c_gam_3",c_gam_3)
    call nml_read(filename,"atm_par","c_gam_4",c_gam_4)
    call nml_read(filename,"atm_par","c_gam_5",c_gam_5)
    call nml_read(filename,"atm_par","c_gam_6",c_gam_6)
    call nml_read(filename,"atm_par","gams_max_lnd",gams_max_lnd)
    call nml_read(filename,"atm_par","gams_min_ocn",gams_min_ocn)
    call nml_read(filename,"atm_par","gams_max_ocn",gams_max_ocn)
    call nml_read(filename,"atm_par","hgams",hgams)
    call nml_read(filename,"atm_par","hgamt",hgamt)
    call nml_read(filename,"atm_par","c_gam_rel",c_gam_rel)
    call nml_read(filename,"atm_par","nsmooth_gam",nsmooth_gam)
    call nml_read(filename,"atm_par","nsmooth_cld",nsmooth_cld)
    call nml_read(filename,"atm_par","i_tsl",i_tsl)
    call nml_read(filename,"atm_par","i_tslz",i_tslz)
    call nml_read(filename,"atm_par","c_tsl_gam",c_tsl_gam)
    call nml_read(filename,"atm_par","c_tsl_gam_ice",c_tsl_gam_ice)
    call nml_read(filename,"atm_par","tsl_gams_min_lnd",tsl_gams_min_lnd)
    call nml_read(filename,"atm_par","tsl_gams_min_ice",tsl_gams_min_ice)
    call nml_read(filename,"atm_par","c_wrt_1",c_wrt_1)
    call nml_read(filename,"atm_par","c_wrt_2",c_wrt_2)
    call nml_read(filename,"atm_par","c_wrt_3",c_wrt_3)
    call nml_read(filename,"atm_par","c_wrt_4",c_wrt_4)
    call nml_read(filename,"atm_par","hcld_base",hcld_base)
    call nml_read(filename,"atm_par","i_lw_cld",i_lw_cld)
    call nml_read(filename,"atm_par","c_lw_clot",c_lw_clot)
    call nml_read(filename,"atm_par","ak_co2",ak_co2)
    call nml_read(filename,"atm_par","beta_co2",beta_co2)
    call nml_read(filename,"atm_par","ak_wv",ak_wv)
    call nml_read(filename,"atm_par","a_vap",a_vap)
    call nml_read(filename,"atm_par","beta_vap",beta_vap)
    call nml_read(filename,"atm_par","a2_vap",a2_vap)
    call nml_read(filename,"atm_par","beta2_vap",beta2_vap)
    call nml_read(filename,"atm_par","a3_vap",a3_vap)
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
    call nml_read(filename,"atm_par","c_syn_1",c_syn_1)
    call nml_read(filename,"atm_par","c_syn_2",c_syn_2)
    call nml_read(filename,"atm_par","c_syn_3",c_syn_3)
    call nml_read(filename,"atm_par","c_syn_4",c_syn_4)
    call nml_read(filename,"atm_par","c_syn_5",c_syn_5)
    call nml_read(filename,"atm_par","c_syn_6",c_syn_6)
    call nml_read(filename,"atm_par","c_syn_7",c_syn_7)
    call nml_read(filename,"atm_par","c_syn_8",c_syn_8)
    call nml_read(filename,"atm_par","windmin",windmin)
    call nml_read(filename,"atm_par","synsurmin",synsurmin)
    call nml_read(filename,"atm_par","c_wind_ele",c_wind_ele)
    call nml_read(filename,"atm_par","c_diffx_dse",c_diffx_dse)
    call nml_read(filename,"atm_par","c_diffx_wtr",c_diffx_wtr)
    call nml_read(filename,"atm_par","c_diff_dse",c_diff_dse)
    call nml_read(filename,"atm_par","i_diff_wtr",i_diff_wtr)
    call nml_read(filename,"atm_par","c_diff_wtr",c_diff_wtr)
    call nml_read(filename,"atm_par","i_diff_dst",i_diff_dst)
    call nml_read(filename,"atm_par","i_kata_wind",i_kata_wind)
    call nml_read(filename,"atm_par","h_kata",h_kata)
    call nml_read(filename,"atm_par","cd0_ocn",cd0_ocn)
    call nml_read(filename,"atm_par","cd0_sic",cd0_sic)
    call nml_read(filename,"atm_par","i_zoro",i_zoro)
    call nml_read(filename,"atm_par","i_acbar",i_acbar)
    call nml_read(filename,"atm_par","acbar_max",acbar_max)
    call nml_read(filename,"atm_par","acbar_scale",acbar_scale)
    call nml_read(filename,"atm_par","c_cld_1",c_cld_1)
    call nml_read(filename,"atm_par","c_cld_2",c_cld_2)
    call nml_read(filename,"atm_par","c_cld_3",c_cld_3)
    call nml_read(filename,"atm_par","c_cld_4",c_cld_4)
    call nml_read(filename,"atm_par","c_cld_5",c_cld_5)
    call nml_read(filename,"atm_par","c_cld_55",c_cld_55)
    call nml_read(filename,"atm_par","c_cld_6",c_cld_6)
    call nml_read(filename,"atm_par","c_cld_7",c_cld_7)
    call nml_read(filename,"atm_par","c_cld_8",c_cld_8)
    call nml_read(filename,"atm_par","l_cld_low_ice",l_cld_low_ice)
    call nml_read(filename,"atm_par","c_hcld_1",c_hcld_1)
    call nml_read(filename,"atm_par","c_hcld_2",c_hcld_2)
    call nml_read(filename,"atm_par","c_hcld_3",c_hcld_3)
    call nml_read(filename,"atm_par","c_hcld_4",c_hcld_4)
    call nml_read(filename,"atm_par","cld_max",cld_max)
    call nml_read(filename,"atm_par","c_woro",c_woro)
    call nml_read(filename,"atm_par","c_weff",c_weff)
    call nml_read(filename,"atm_par","i_rbstr",i_rbstr)
    call nml_read(filename,"atm_par","i_pbl",i_pbl)
    call nml_read(filename,"atm_par","hpbl",hpbl)
    call nml_read(filename,"atm_par","rh_max",rh_max)
    call nml_read(filename,"atm_par","rskin_ocn_min",rskin_ocn_min)
    call nml_read(filename,"atm_par","rh_strat",rh_strat)
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
    call nml_read(filename,"atm_par","nsmooth_cda",nsmooth_cda)
    call nml_read(filename,"atm_par","nsmooth_weff",nsmooth_weff)
    call nml_read(filename,"atm_par","nsmooth_aslp",nsmooth_aslp)
    call nml_read(filename,"atm_par","nsmooth_aslp_eq",nsmooth_aslp_eq)
    call nml_read(filename,"atm_par","nj_eq",nj_eq)
    call nml_read(filename,"atm_par","nsmooth_aslp_topo",nsmooth_aslp_topo)
    call nml_read(filename,"atm_par","nsmooth_acbar",nsmooth_acbar)

    call nml_read(filename,"atm_par","l_write_timer",l_write_timer)

    call nml_read(filename,"atm_par","l_daily_output",l_daily_output)
    call nml_read(filename,"atm_par","l_output_flx3d",l_output_flx3d)
    call nml_read(filename,"atm_par","l_output_extended",l_output_extended)

    call nml_read(filename,"atm_par","tam_init",tam_init)
 
    return

  end subroutine atm_par_load

end module atm_params
