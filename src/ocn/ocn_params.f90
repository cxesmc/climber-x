!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : o c n _ p a r a m s
!
!  Purpose : ocean parameter definition
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Neil R. Edwards and Matteo Willeit
!
! This file is part of CLIMBER-X.
!
! This file was ported from the original c-GOLDSTEIN model,
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
module ocn_params

  use precision, only : wp
  use nml
  use timer, only : dt_ocn, sec_day
  use control, only : out_dir

  implicit none

     real(wp) :: dt !! time step [s]

     integer :: n_tracers_tot    !! total number of tracers (ocn+bgc) []
     integer :: n_tracers_ocn    !! number of ocean model tracers []
     integer :: n_tracers_bgc    !! number of bgc tracers []
     integer :: n_tracers_trans  !! number of transported tracers []
     integer, allocatable :: idx_tracers_trans(:)  !! index of transported tracers
     integer :: i_age   !! index of age tracer
     integer :: i_dye   !! index of dye tracer
     integer :: i_cons  !! index of conservative tracer
     integer :: i_cfc11 !! index of CFC11 tracer
     integer :: i_cfc12 !! index of CFC12 tracer

     ! Namelist parameters and flags
     integer :: i_init
     integer :: nlayers
     integer :: i_smooth
     real(wp) :: smooth_fac
     real(wp), dimension(:), allocatable :: zw_in
     integer :: i_isl
     logical :: l_isl_ant
     logical :: l_isl_aus
     logical :: l_isl_grl
     logical :: l_isl_ame
     logical :: l_bering_flow
     real(wp) :: c_bering
     integer :: n_isl
     real(wp) :: lon_isl(20)
     real(wp) :: lat_isl(20)
     real(wp) :: dbl

     real(wp) :: rho0 !! reference density ocean water [kg/m3]

     real(wp) :: urelax
     real(wp) :: tau_scale

     type drag_par_type
       real(wp) :: adrag
       integer :: drag_topo_n
       real(wp) :: drag_topo_fac
       real(wp) :: drag_frac_fac
       real(wp) :: drag_topo_scale_eq
       real(wp) :: z_drag_shallow
     end type
     type(drag_par_type) :: drag_par

     integer :: i_advection
     integer :: i_conv_shuffle
     logical :: l_conv_shuffle_passive
     logical :: l_mix_bgc_all
     logical :: l_mld
     real(wp) :: pe_buoy_coeff, ke_tau_coeff, ke_wind_dec
     integer :: i_eos
     real(wp) :: drhcor_max
     integer :: i_diff
     integer :: i_diff_dia
     logical :: l_diff33_impl
     real(wp) :: diff_iso, diff_dia_zref
     real(wp) :: diff_dia_min, diff_dia_bgc_min, diff_dia_max, diff_dia_ref
     logical :: l_diff_dia_strat
     real(wp) :: alpha_strat
     real(wp) :: brunt_vaisala_ref
     real(wp) :: slope_max 
     real(wp) :: diff_gm
     real(wp) :: tau_sst, tau_sss
     integer :: n_cells_dist_runoff
     integer :: n_cells_dist_calving
     integer :: n_cells_dist_weath
     logical :: l_fw_corr
     logical :: l_fw_melt_ice_sep
     integer :: i_fw
     integer :: i_brines
     integer :: i_brines_z
     real(wp) :: frac_brines
     real(wp) :: z_mix_brines
     real(wp) :: relax_run
     real(wp) :: relax_calv
     real(wp) :: relax_bmelt
     real(wp) :: scale_runoff_ice
     ! peak and background temperatures for i_init option 3
     real(wp) :: init3_peak = 25.0_wp
     real(wp) :: init3_bg   = 15.0_wp
     integer :: i_saln0
     real(wp) :: saln0_const
     logical :: l_salinity_restore
     real(wp) :: shelf_depth
     logical :: age_tracer
     logical :: dye_tracer
     logical :: cons_tracer
     logical :: l_cfc
     logical :: l_hosing
     integer :: i_hosing
     integer :: hosing_basin
     integer :: year_hosing_ini
     integer :: year_hosing_end
     integer :: year_hosing_ramp
     real(wp) :: hosing_ini, hosing_trend, hosing_sigma
     real(wp) :: lat_min_hosing, lat_max_hosing
     real(wp) :: lon_min_hosing, lon_max_hosing
     integer :: hosing_comp_basin
     real(wp) :: lat_min_hosing_comp, lat_max_hosing_comp
     real(wp) :: lon_min_hosing_comp, lon_max_hosing_comp
     character (len=256) :: hosing_file
     logical :: l_noise_fw 
     logical :: is_fw_noise_Sv
     logical :: l_noise_flx 
     integer :: i_noise 
     logical :: l_noise_only_coast
     integer :: noise_basin 
     real(wp) :: noise_amp_fw 
     real(wp) :: noise_amp_flx 
     real(wp) :: noise_autocorr
     real(wp) :: noise_period 
     real(wp) :: lat_min_noise 
     real(wp) :: lat_max_noise 
     real(wp) :: lon_min_noise 
     real(wp) :: lon_max_noise 

     logical :: l_flux_adj_atl, l_flux_adj_ant, l_flux_adj_pac
     real(wp) :: flux_adj_atl, flux_adj_ant, flux_adj_pac
     integer :: nj_flux_adj_ant
     real(wp) :: lat_min_flux_adj_atl, lat_max_flux_adj_atl
     real(wp) :: lat_min_flux_adj_pac, lat_max_flux_adj_pac

     integer :: i_fwf_buoy 
     integer :: i_alphabeta
     real(wp) :: depth_buoy

     logical :: l_ocn_input_fix
     logical :: l_ocn_fix_wind           ! used prescribed wind stress from average of first 30 years of simulation?
     logical :: l_ocn_fix_fw             ! used prescribed surface freshwater flux from average of first 30 years of simulation?
     logical :: l_ocn_fix_flx            ! used prescribed surface heat flux from average of first 30 years of simulation?
     integer :: i_ocn_input_fix          ! what fixed ocean input to use, only if any l_ocn_fix_*==T : 
     logical :: l_ocn_input_fix_write    ! write average daily ocean input variables from last 30 years of simulation to file?
     character(len=256) :: ocn_input_fix_file  ! file with daily climatology of ocean input fields to be used for simulations with any l_ocn_fix_*==T

     logical :: l_q_geo

     logical :: l_daily_output
     logical :: l_output_extended

     real(wp) :: fcormin
     real(wp) :: fcormin_ref

     real(wp), dimension(:), allocatable :: fcor, fcorv
     real(wp), dimension(:,:,:), allocatable :: drag
     real(wp), dimension(:,:,:), allocatable :: drag_bcl

     real(wp), dimension(:), allocatable :: mlddec, mlddecd

     real(wp), dimension(:,:,:), allocatable :: diff_dia
     real(wp), dimension(:,:,:), allocatable :: diff_dia_bgc
     real(wp), dimension(:,:), allocatable :: slope_crit
     real(wp), dimension(:), allocatable :: diffx_max
     real(wp), dimension(:), allocatable :: diffy_max
    
     real(wp), dimension(:,:), allocatable :: rtv, rtv3

contains

    subroutine ocn_params_init

    implicit none

     ! time step
     dt = dt_ocn

    ! read model settings and parameters from namelist
    call ocn_par_load(trim(out_dir)//"/ocn_par.nml")

 
    return

    end subroutine ocn_params_init

subroutine ocn_par_load(filename)

    implicit none

    character (len=*) :: filename

    ! Read parameters from file
    write(*,*) "ocean parameters ==========="
    call nml_read(filename,"ocn_par","i_init",i_init)
    call nml_read(filename,"ocn_par","nlayers",nlayers)
    allocate(zw_in(nlayers+1))
    call nml_read(filename,"ocn_par","levels",zw_in)
    call nml_read(filename,"ocn_par","i_smooth",i_smooth)
    call nml_read(filename,"ocn_par","smooth_fac",smooth_fac)
    call nml_read(filename,"ocn_par","i_isl",i_isl)
    call nml_read(filename,"ocn_par","l_isl_ant",l_isl_ant)
    call nml_read(filename,"ocn_par","l_isl_aus",l_isl_aus)
    call nml_read(filename,"ocn_par","l_isl_grl",l_isl_grl)
    call nml_read(filename,"ocn_par","l_isl_ame",l_isl_ame)
    call nml_read(filename,"ocn_par","l_bering_flow",l_bering_flow)
    call nml_read(filename,"ocn_par","c_bering",c_bering)
    call nml_read(filename,"ocn_par","n_isl",n_isl)
    call nml_read(filename,"ocn_par","lon_isl",lon_isl)
    call nml_read(filename,"ocn_par","lat_isl",lat_isl)
    call nml_read(filename,"ocn_par","dbl",dbl)

    call nml_read(filename,"ocn_par","rho0",rho0)
    call nml_read(filename,"ocn_par","urelax",urelax)
    call nml_read(filename,"ocn_par","tau_scale",tau_scale)
    call nml_read(filename,"ocn_par","adrag",drag_par%adrag)
    call nml_read(filename,"ocn_par","drag_topo_n",drag_par%drag_topo_n)
    call nml_read(filename,"ocn_par","drag_topo_fac",drag_par%drag_topo_fac)
    call nml_read(filename,"ocn_par","drag_topo_scale_eq",drag_par%drag_topo_scale_eq)
    call nml_read(filename,"ocn_par","drag_frac_fac",drag_par%drag_frac_fac)
    call nml_read(filename,"ocn_par","z_drag_shallow",drag_par%z_drag_shallow)

    call nml_read(filename,"ocn_par","i_eos",i_eos)

    call nml_read(filename,"ocn_par","drhcor_max",drhcor_max)

    call nml_read(filename,"ocn_par","fcormin",fcormin)
    fcormin_ref = fcormin

    call nml_read(filename,"ocn_par","i_advection",i_advection)

    call nml_read(filename,"ocn_par","i_diff",i_diff)
    call nml_read(filename,"ocn_par","i_diff_dia",i_diff_dia)
    call nml_read(filename,"ocn_par","l_diff33_impl",l_diff33_impl)
    call nml_read(filename,"ocn_par","diff_iso",diff_iso)
    call nml_read(filename,"ocn_par","diff_dia_zref",diff_dia_zref)
    call nml_read(filename,"ocn_par","diff_dia_min",diff_dia_min)
    call nml_read(filename,"ocn_par","diff_dia_max",diff_dia_max)
    call nml_read(filename,"ocn_par","diff_dia_bgc_min",diff_dia_bgc_min)
    call nml_read(filename,"ocn_par","diff_dia_ref",diff_dia_ref)
    call nml_read(filename,"ocn_par","l_diff_dia_strat",l_diff_dia_strat)
    call nml_read(filename,"ocn_par","alpha_strat",alpha_strat)
    call nml_read(filename,"ocn_par","brunt_vaisala_ref",brunt_vaisala_ref)
    call nml_read(filename,"ocn_par","slope_max",slope_max)
    call nml_read(filename,"ocn_par","diff_gm",diff_gm)

    call nml_read(filename,"ocn_par","i_conv_shuffle",i_conv_shuffle)
    call nml_read(filename,"ocn_par","l_conv_shuffle_passive",l_conv_shuffle_passive)
    call nml_read(filename,"ocn_par","l_mix_bgc_all",l_mix_bgc_all)

    call nml_read(filename,"ocn_par","l_mld",l_mld)
    call nml_read(filename,"ocn_par","pe_buoy_coeff",pe_buoy_coeff)
    call nml_read(filename,"ocn_par","ke_tau_coeff",ke_tau_coeff)
    call nml_read(filename,"ocn_par","ke_wind_dec",ke_wind_dec)

    call nml_read(filename,"ocn_par","tau_sst",tau_sst)
    call nml_read(filename,"ocn_par","tau_sss",tau_sss)

    call nml_read(filename,"ocn_par","shelf_depth",shelf_depth)

    call nml_read(filename,"ocn_par","age_tracer",age_tracer)
    call nml_read(filename,"ocn_par","dye_tracer",dye_tracer)
    call nml_read(filename,"ocn_par","cons_tracer",cons_tracer)
    call nml_read(filename,"ocn_par","l_cfc",l_cfc)

    call nml_read(filename,"ocn_par","n_cells_dist_runoff",n_cells_dist_runoff)
    call nml_read(filename,"ocn_par","n_cells_dist_calving",n_cells_dist_calving)
    call nml_read(filename,"ocn_par","n_cells_dist_weath",n_cells_dist_weath)
    call nml_read(filename,"ocn_par","l_fw_corr",l_fw_corr)
    call nml_read(filename,"ocn_par","l_fw_melt_ice_sep",l_fw_melt_ice_sep)
    call nml_read(filename,"ocn_par","i_fw",i_fw)
    call nml_read(filename,"ocn_par","init3_peak",init3_peak)
    call nml_read(filename,"ocn_par","init3_bg",init3_bg)
    call nml_read(filename,"ocn_par","i_saln0",i_saln0)
    call nml_read(filename,"ocn_par","saln0_const",saln0_const)
    call nml_read(filename,"ocn_par","l_salinity_restore",l_salinity_restore)
    call nml_read(filename,"ocn_par","i_brines",i_brines)
    call nml_read(filename,"ocn_par","i_brines_z",i_brines_z)
    call nml_read(filename,"ocn_par","frac_brines",frac_brines)
    call nml_read(filename,"ocn_par","z_mix_brines",z_mix_brines)
    call nml_read(filename,"ocn_par","relax_run",relax_run)
    call nml_read(filename,"ocn_par","relax_calv",relax_calv)
    call nml_read(filename,"ocn_par","relax_bmelt",relax_bmelt)
    call nml_read(filename,"ocn_par","scale_runoff_ice",scale_runoff_ice)

    call nml_read(filename,"ocn_par","l_hosing",l_hosing)
    call nml_read(filename,"ocn_par","i_hosing",i_hosing)
    call nml_read(filename,"ocn_par","hosing_basin",hosing_basin)
    call nml_read(filename,"ocn_par","lat_min_hosing",lat_min_hosing)
    call nml_read(filename,"ocn_par","lat_max_hosing",lat_max_hosing)
    call nml_read(filename,"ocn_par","lon_min_hosing",lon_min_hosing)
    call nml_read(filename,"ocn_par","lon_max_hosing",lon_max_hosing)
    call nml_read(filename,"ocn_par","hosing_comp_basin",hosing_comp_basin)
    call nml_read(filename,"ocn_par","lat_min_hosing_comp",lat_min_hosing_comp)
    call nml_read(filename,"ocn_par","lat_max_hosing_comp",lat_max_hosing_comp)
    call nml_read(filename,"ocn_par","lon_min_hosing_comp",lon_min_hosing_comp)
    call nml_read(filename,"ocn_par","lon_max_hosing_comp",lon_max_hosing_comp)
    call nml_read(filename,"ocn_par","hosing_ini",hosing_ini)
    call nml_read(filename,"ocn_par","hosing_trend",hosing_trend)
    call nml_read(filename,"ocn_par","hosing_sigma",hosing_sigma)
    call nml_read(filename,"ocn_par","year_hosing_ini",year_hosing_ini)
    call nml_read(filename,"ocn_par","year_hosing_end",year_hosing_end)
    call nml_read(filename,"ocn_par","year_hosing_ramp",year_hosing_ramp)
    call nml_read(filename,"ocn_par","hosing_file",hosing_file)

    call nml_read(filename,"ocn_par","l_noise_fw   ",l_noise_fw   )
    call nml_read(filename,"ocn_par","is_fw_noise_Sv",is_fw_noise_Sv)
    call nml_read(filename,"ocn_par","l_noise_flx  ",l_noise_flx  )
    call nml_read(filename,"ocn_par","i_noise      ",i_noise      )
    call nml_read(filename,"ocn_par","l_noise_only_coast",l_noise_only_coast)
    call nml_read(filename,"ocn_par","noise_basin  ",noise_basin  )
    call nml_read(filename,"ocn_par","noise_amp_fw ",noise_amp_fw )
    if (.not.is_fw_noise_Sv) then
      noise_amp_fw = noise_amp_fw/sec_day     ! kg/m2/day -> kg/m2/s
    endif
    call nml_read(filename,"ocn_par","noise_amp_flx",noise_amp_flx)
    call nml_read(filename,"ocn_par","noise_autocorr",noise_autocorr )
    call nml_read(filename,"ocn_par","noise_period ",noise_period )
    call nml_read(filename,"ocn_par","lat_min_noise",lat_min_noise)
    call nml_read(filename,"ocn_par","lat_max_noise",lat_max_noise)
    call nml_read(filename,"ocn_par","lon_min_noise",lon_min_noise)
    call nml_read(filename,"ocn_par","lon_max_noise",lon_max_noise)

    call nml_read(filename,"ocn_par","l_flux_adj_atl",l_flux_adj_atl)
    call nml_read(filename,"ocn_par","flux_adj_atl",flux_adj_atl)
    call nml_read(filename,"ocn_par","l_flux_adj_ant",l_flux_adj_ant)
    call nml_read(filename,"ocn_par","flux_adj_ant",flux_adj_ant)
    call nml_read(filename,"ocn_par","nj_flux_adj_ant",nj_flux_adj_ant)
    call nml_read(filename,"ocn_par","l_flux_adj_pac",l_flux_adj_pac)
    call nml_read(filename,"ocn_par","flux_adj_pac",flux_adj_pac)
    call nml_read(filename,"ocn_par","lat_min_flux_adj_atl",lat_min_flux_adj_atl)
    call nml_read(filename,"ocn_par","lat_max_flux_adj_atl",lat_max_flux_adj_atl)
    call nml_read(filename,"ocn_par","lat_min_flux_adj_pac",lat_min_flux_adj_pac)
    call nml_read(filename,"ocn_par","lat_max_flux_adj_pac",lat_max_flux_adj_pac)
    if (.not.l_flux_adj_atl) flux_adj_atl = 0._wp
    if (.not.l_flux_adj_pac) flux_adj_pac = 0._wp
    if (.not.l_flux_adj_ant) flux_adj_ant = 0._wp

    call nml_read(filename,"ocn_par","i_fwf_buoy",i_fwf_buoy)
    call nml_read(filename,"ocn_par","i_alphabeta",i_alphabeta)
    call nml_read(filename,"ocn_par","depth_buoy",depth_buoy)

    call nml_read(filename,"ocn_par","l_ocn_fix_wind",l_ocn_fix_wind)
    call nml_read(filename,"ocn_par","l_ocn_fix_fw",  l_ocn_fix_fw )
    call nml_read(filename,"ocn_par","l_ocn_fix_flx", l_ocn_fix_flx)
    if (l_ocn_fix_wind .or. l_ocn_fix_fw .or. l_ocn_fix_flx) then
      l_ocn_input_fix = .true.
    else
      l_ocn_input_fix = .false.
    endif
    call nml_read(filename,"ocn_par","i_ocn_input_fix",i_ocn_input_fix)
    call nml_read(filename,"ocn_par","l_ocn_input_fix_write",l_ocn_input_fix_write)
    call nml_read(filename,"ocn_par","ocn_input_fix_file",ocn_input_fix_file)

    call nml_read(filename,"ocn_par","l_q_geo",l_q_geo)

    call nml_read(filename,"ocn_par","l_daily_output",l_daily_output)
    call nml_read(filename,"ocn_par","l_output_extended",l_output_extended)

   return

end subroutine ocn_par_load

end module ocn_params
