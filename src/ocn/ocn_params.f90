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
  use timer, only : dt_ocn, sec_day
  use control, only : out_dir
  use nml
  use ncio

  implicit none

     real(wp) :: dt !! time step [s]

     ! ramped ocean sub-cycling for cold-start spin-up
     logical :: l_ocn_dt_ramp       !! enable ramped sub-cycling of the ocean time step on cold start (ocn_restart=F)?
     integer :: ocn_dt_ramp_nsub0   !! initial number of sub-steps N0 (power of two), effective initial step = dt/N0 []
     integer :: ocn_dt_ramp_years   !! number of years over which N halves from N0 down to 1 [years]
     integer :: n_sub               !! current number of sub-steps per ocean call []

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
     logical :: l_limit_ratio   !! limit the neighbour ocean-depth ratio of the bathymetry by deepening shallow cells (despike for barotropic stability)?
     real(wp) :: depth_ratio_max  !! maximum allowed neighbour depth ratio (threshold for l_limit_ratio) []
     real(wp), dimension(:), allocatable :: zw_in
     integer :: i_isl
     real(wp) :: isl_area_min
     logical :: l_isl_ant
     logical :: l_isl_aus
     logical :: l_isl_grl
     logical :: l_isl_ame
     logical :: l_bering_flow
     real(wp) :: c_bering
     real(wp) :: bering_tf_fix   !! Sv, prescribed barotropic Bering Strait throughflow (positive north), < -900 = off
     real(wp) :: drake_tf_fix    !! Sv, prescribed barotropic Drake Passage throughflow (positive east), < -900 = off
     real(wp) :: indo_tf_fix     !! Sv, prescribed barotropic Indonesian throughflow (positive east), < -900 = off
     real(wp) :: davis_tf_fix    !! Sv, prescribed barotropic Davis Strait throughflow (positive north), < -900 = off
     integer :: n_isl
     real(wp) :: lon_isl(20)
     real(wp) :: lat_isl(20)
     real(wp) :: dbl

     real(wp) :: rho0 !! reference density ocean water [kg/m3]

     real(wp) :: urelax
     real(wp) :: tau_scale

     type drag_par_type
       real(wp) :: adrag           !! background drag of the barotropic and baroclinic velocities (namelist in days)
       real(wp) :: drag_frac       !! partial-cell drag rate at f_ocn=0 (namelist in days, 0 = off), times max(0,1-f_face/f_frac_crit)
       real(wp) :: f_frac_crit     !! ocean fraction of a face above which there is no partial-cell drag
       real(wp) :: cd_bot          !! bottom drag coefficient, 0 = no bottom friction
       real(wp) :: u_bot           !! m/s, rms near-bottom speed (tides, eddies) linearising the quadratic bottom drag
       real(wp) :: drag_cb         !! checkerboard friction rate (namelist in days, 0 = off), see stencil_coef in invert.f90
       real(wp) :: adrag_bcl       !! baroclinic drag away from the equator (namelist in days, 0 = adrag)
       real(wp) :: lat_drag_bcl    !! deg, e-folding latitude of the transition from adrag (equator) to adrag_bcl
       integer  :: i_fric_ret      !! return flow of the barotropic friction: 0 = uniform over the column, 1 = in the bottom level
       real(wp) :: lat_fric_ret    !! deg, equatorial regularisation of 1/f in the bottom return, f/(f^2+f0^2) with f0 = f(lat_fric_ret)
     end type
     type(drag_par_type) :: drag_par

     integer :: i_conv_shuffle
     logical :: l_conv_shuffle_passive
     logical :: l_downslope           !! downsloping flow of dense shelf water (Campin and Goosse 1999)?
     real(wp) :: c_downslope         !! coefficient of the gravity-current speed sqrt(g'H) []
     real(wp) :: ent_downslope       !! entrainment of ambient water by the plume, fraction of its volume per km of descent []
     real(wp) :: drho_downslope_min  !! minimum shelf-minus-ambient density contrast for a downslope flow [kg/m3]
     real(wp) :: f_downslope_max     !! maximum fraction of a cell volume exchanged per time step []
     real(wp) :: z_downslope_max     !! maximum bottom depth of a source (shelf) cell [m]
     real(wp) :: pe_buoy_coeff, ke_tau_coeff, ke_wind_dec
     integer :: i_eos
     real(wp) :: drhcor_max
     integer :: i_cor_form      
     integer :: nbw             !! half band width of the streamfunction matrix (maxi+1 for the 5-point stencil, 2*maxi with drag_cb > 0)
     real(wp), allocatable :: qcor(:,:)  !! f/H at tracer points, used by the J3 Coriolis term (i_cor_form=1)
     integer :: i_diff_dia
     real(wp) :: diff_iso, diff_dia_zref
     real(wp) :: diff_dia_min, diff_dia_max, diff_dia_ref
     logical :: l_diff_dia_strat
     real(wp) :: alpha_strat
     real(wp) :: brunt_vaisala_ref
     real(wp) :: slope_max 
     integer :: i_taper         !! isoneutral slope taper: 1 = (Sc/S)^2 (Gerdes 1991), 2 = tanh (Danabasoglu and McWilliams 1995)
     real(wp) :: slope_width    !! width of the tanh taper, used only if i_taper==2 []
     real(wp) :: drho_hbl       !! potential density threshold relative to 10 m for the boundary layer depth [kg/m3]
     real(wp) :: diff_gm
     real(wp) :: tau_sst, tau_sss
     integer :: n_cells_dist_runoff
     integer :: n_cells_dist_calving
     integer :: n_cells_dist_weath
     logical :: l_fw_corr
     logical :: l_fw_melt_ice_sep
     integer :: i_fw
     logical :: l_brines
     real(wp) :: frac_brines
     real(wp) :: drho_brines_coast
     real(wp) :: relax_run
     real(wp) :: relax_calv
     real(wp) :: relax_bmelt
     real(wp) :: scale_runoff_ice
     real(wp) :: scale_calving_ice
     integer :: i_scale_dhdt_ice
     real(wp) :: scale_dhdt_ice
     character (len=256) :: scale_dhdt_ice_file
     ! peak and background temperatures for i_init option 3
     real(wp) :: init3_peak = 25.0_wp
     real(wp) :: init3_bg   = 15.0_wp
     integer :: i_saln0
     real(wp) :: saln0_const
     logical :: l_salinity_restore
     real(wp) :: ocn_depth_min
     logical :: age_tracer
     logical :: dye_tracer
     logical :: cons_tracer
     logical :: l_cfc
     logical :: l_hosing
     integer, parameter :: n_hosing_domain_max=10  !! maximum number of hosing domains
     integer :: n_hosing_domain  !! number of hosing domains
     character (len=256) :: hosing_domain_name(n_hosing_domain_max)
     integer :: i_hosing_comp
     integer :: hosing_comp_basin
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
     real(wp), dimension(:,:), allocatable :: slope_crit
     real(wp), dimension(:), allocatable :: diffx_max
     real(wp), dimension(:), allocatable :: diffy_max
    
     real(wp), dimension(:,:), allocatable :: rtv, rtv3

     real(wp), allocatable :: scale_dhdt_ice_time(:)
     real(wp), allocatable :: scale_dhdt_ice_data(:)

contains

    subroutine ocn_params_init

    implicit none

     ! time step
     dt = dt_ocn

     ! default to a single sub-step (overridden each year in ocn_update when l_ocn_dt_ramp)
     n_sub = 1

    ! read model settings and parameters from namelist
    call ocn_par_load(trim(out_dir)//"/ocn_par.nml")

    i_age = 0
    i_dye = 0
 
    return

    end subroutine ocn_params_init

subroutine ocn_par_load(filename)

    implicit none

    character (len=*) :: filename
    integer :: n
    integer :: ntime


    ! Read parameters from file
    write(*,*) "ocean parameters ==========="
    call nml_read(filename,"ocn_par","l_ocn_dt_ramp",l_ocn_dt_ramp)
    call nml_read(filename,"ocn_par","ocn_dt_ramp_nsub0",ocn_dt_ramp_nsub0)
    call nml_read(filename,"ocn_par","ocn_dt_ramp_years",ocn_dt_ramp_years)
    call nml_read(filename,"ocn_par","i_init",i_init)
    call nml_read(filename,"ocn_par","nlayers",nlayers)
    allocate(zw_in(nlayers+1))
    call nml_read(filename,"ocn_par","levels",zw_in)
    call nml_read(filename,"ocn_par","i_smooth",i_smooth)
    call nml_read(filename,"ocn_par","smooth_fac",smooth_fac)
    call nml_read(filename,"ocn_par","l_limit_ratio",l_limit_ratio)
    call nml_read(filename,"ocn_par","depth_ratio_max",depth_ratio_max)
    call nml_read(filename,"ocn_par","i_isl",i_isl)
    call nml_read(filename,"ocn_par","isl_area_min",isl_area_min)
    call nml_read(filename,"ocn_par","l_isl_ant",l_isl_ant)
    call nml_read(filename,"ocn_par","l_isl_aus",l_isl_aus)
    call nml_read(filename,"ocn_par","l_isl_grl",l_isl_grl)
    call nml_read(filename,"ocn_par","l_isl_ame",l_isl_ame)
    call nml_read(filename,"ocn_par","l_bering_flow",l_bering_flow)
    call nml_read(filename,"ocn_par","c_bering",c_bering)
    bering_tf_fix = -999._wp
    call nml_read(filename,"ocn_par","bering_tf_fix",bering_tf_fix)
    drake_tf_fix = -999._wp
    call nml_read(filename,"ocn_par","drake_tf_fix",drake_tf_fix)
    indo_tf_fix = -999._wp
    call nml_read(filename,"ocn_par","indo_tf_fix",indo_tf_fix)
    davis_tf_fix = -999._wp
    call nml_read(filename,"ocn_par","davis_tf_fix",davis_tf_fix)
    call nml_read(filename,"ocn_par","n_isl",n_isl)
    call nml_read(filename,"ocn_par","lon_isl",lon_isl)
    call nml_read(filename,"ocn_par","lat_isl",lat_isl)
    call nml_read(filename,"ocn_par","dbl",dbl)

    call nml_read(filename,"ocn_par","rho0",rho0)
    call nml_read(filename,"ocn_par","urelax",urelax)
    call nml_read(filename,"ocn_par","tau_scale",tau_scale)
    call nml_read(filename,"ocn_par","adrag",drag_par%adrag)
    call nml_read(filename,"ocn_par","drag_frac",drag_par%drag_frac)
    call nml_read(filename,"ocn_par","f_frac_crit",drag_par%f_frac_crit)
    call nml_read(filename,"ocn_par","cd_bot",drag_par%cd_bot)
    call nml_read(filename,"ocn_par","u_bot",drag_par%u_bot)
    call nml_read(filename,"ocn_par","drag_cb",drag_par%drag_cb)
    call nml_read(filename,"ocn_par","adrag_bcl",drag_par%adrag_bcl)
    call nml_read(filename,"ocn_par","lat_drag_bcl",drag_par%lat_drag_bcl)
    call nml_read(filename,"ocn_par","i_fric_ret",drag_par%i_fric_ret)
    call nml_read(filename,"ocn_par","lat_fric_ret",drag_par%lat_fric_ret)
    if (drag_par%lat_drag_bcl.le.0._wp .or. drag_par%lat_fric_ret.le.0._wp) then
      print *,'ERROR: lat_drag_bcl and lat_fric_ret must be > 0',drag_par%lat_drag_bcl,drag_par%lat_fric_ret
      stop
    endif
    if (drag_par%i_fric_ret.lt.0 .or. drag_par%i_fric_ret.gt.1) then
      print *,'ERROR: i_fric_ret must be 0 or 1, i_fric_ret =',drag_par%i_fric_ret
      stop
    endif

    call nml_read(filename,"ocn_par","i_eos",i_eos)
    if (i_eos<0 .or. i_eos>4) then
      print *,'ERROR: i_eos must be 0-4, i_eos =',i_eos
      stop
    endif

    call nml_read(filename,"ocn_par","drhcor_max",drhcor_max)
    call nml_read(filename,"ocn_par","i_cor_form",i_cor_form)

    call nml_read(filename,"ocn_par","fcormin",fcormin)
    fcormin_ref = fcormin


    call nml_read(filename,"ocn_par","i_diff_dia",i_diff_dia)
    call nml_read(filename,"ocn_par","diff_iso",diff_iso)
    call nml_read(filename,"ocn_par","diff_dia_zref",diff_dia_zref)
    call nml_read(filename,"ocn_par","diff_dia_min",diff_dia_min)
    call nml_read(filename,"ocn_par","diff_dia_max",diff_dia_max)
    call nml_read(filename,"ocn_par","diff_dia_ref",diff_dia_ref)
    call nml_read(filename,"ocn_par","l_diff_dia_strat",l_diff_dia_strat)
    call nml_read(filename,"ocn_par","alpha_strat",alpha_strat)
    call nml_read(filename,"ocn_par","brunt_vaisala_ref",brunt_vaisala_ref)
    call nml_read(filename,"ocn_par","slope_max",slope_max)
    call nml_read(filename,"ocn_par","i_taper",i_taper)
    call nml_read(filename,"ocn_par","slope_width",slope_width)
    call nml_read(filename,"ocn_par","drho_hbl",drho_hbl)
    call nml_read(filename,"ocn_par","diff_gm",diff_gm)

    call nml_read(filename,"ocn_par","i_conv_shuffle",i_conv_shuffle)
    call nml_read(filename,"ocn_par","l_conv_shuffle_passive",l_conv_shuffle_passive)
    call nml_read(filename,"ocn_par","l_downslope",l_downslope)
    call nml_read(filename,"ocn_par","c_downslope",c_downslope)
    call nml_read(filename,"ocn_par","ent_downslope",ent_downslope)
    call nml_read(filename,"ocn_par","drho_downslope_min",drho_downslope_min)
    call nml_read(filename,"ocn_par","f_downslope_max",f_downslope_max)
    call nml_read(filename,"ocn_par","z_downslope_max",z_downslope_max)

    call nml_read(filename,"ocn_par","pe_buoy_coeff",pe_buoy_coeff)
    call nml_read(filename,"ocn_par","ke_tau_coeff",ke_tau_coeff)
    call nml_read(filename,"ocn_par","ke_wind_dec",ke_wind_dec)

    call nml_read(filename,"ocn_par","tau_sst",tau_sst)
    call nml_read(filename,"ocn_par","tau_sss",tau_sss)

    call nml_read(filename,"ocn_par","ocn_depth_min",ocn_depth_min)

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
    call nml_read(filename,"ocn_par","l_brines",l_brines)
    call nml_read(filename,"ocn_par","frac_brines",frac_brines)
    call nml_read(filename,"ocn_par","drho_brines_coast",drho_brines_coast)
    call nml_read(filename,"ocn_par","relax_run",relax_run)
    call nml_read(filename,"ocn_par","relax_calv",relax_calv)
    call nml_read(filename,"ocn_par","relax_bmelt",relax_bmelt)
    call nml_read(filename,"ocn_par","scale_runoff_ice",scale_runoff_ice)
    call nml_read(filename,"ocn_par","scale_calving_ice",scale_calving_ice)
    call nml_read(filename,"ocn_par","i_scale_dhdt_ice",i_scale_dhdt_ice)
    call nml_read(filename,"ocn_par","scale_dhdt_ice",scale_dhdt_ice)
    call nml_read(filename,"ocn_par","scale_dhdt_ice_file",scale_dhdt_ice_file)
    if (i_scale_dhdt_ice.eq.1) then
      ntime = nc_size(trim(scale_dhdt_ice_file),"time")
      allocate( scale_dhdt_ice_time(ntime) )
      allocate( scale_dhdt_ice_data(ntime) )
      call nc_read(trim(scale_dhdt_ice_file),"time",scale_dhdt_ice_time)    
      call nc_read(trim(scale_dhdt_ice_file),"scale_dhdt_ice",scale_dhdt_ice_data) 
    endif

    call nml_read(filename,"ocn_par","l_hosing",l_hosing)
    call nml_read(filename,"ocn_par","hosing_domain_name",hosing_domain_name)
    n_hosing_domain = 0
    do n=1,10
      if (len_trim(hosing_domain_name(n)).ne.0) n_hosing_domain=n_hosing_domain+1
    enddo
    if (n_hosing_domain>n_hosing_domain_max) then
      print *,'ERROR: too many hosing domains, n_hosing_domain>n_hosing_domain_max',n_hosing_domain,n_hosing_domain_max
      stop
    endif
    call nml_read(filename,"ocn_par","i_hosing_comp",i_hosing_comp)
    call nml_read(filename,"ocn_par","hosing_comp_basin",hosing_comp_basin)

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
