!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c o n t r o l
!
!  Purpose : reading of control namelist 
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
module control

  use precision, only : wp
  use nml

  implicit none

  logical :: atm_restart, co2_restart, ch4_restart, lnd_restart, ocn_restart, sic_restart, bgc_restart, geo_restart, ice_restart, smb_restart, bmb_restart
  character (len=256) :: restart_in_dir 
  integer :: i_write_restart 
  integer :: n_year_write_restart 
  integer :: years_write_restart(10)

  logical :: flag_atm, flag_co2, flag_ch4, flag_lnd, flag_dust, flag_lakes, flag_ocn, flag_sic, flag_bgc, flag_ice, flag_smb, flag_bmb, flag_geo
  logical :: ocn_restore_sal, ocn_restore_temp
  logical :: atm_fix_tau
  character (len=256) :: ice_model_name
  !character (len=65), allocatable :: ice_domain_name(:)
  integer, parameter :: n_ice_domain_max=10  !! maximum number of ice sheet model domains
  integer :: n_ice_domain  !! number of ice sheet model domains
  character (len=256) :: ice_domain_name(n_ice_domain_max)
  logical :: l_aquaplanet
  logical :: l_aqua_slab

  logical :: l_spinup_cc
  logical :: l_daily_input_save_ocn
  logical :: l_daily_input_save_bgc
  integer :: nyears_spinup_bgc
  integer :: year_start_offline
  integer :: nyear_avg_offline

  logical :: l_feedbacks

  integer :: iorbit
  real(wp) :: ecc_const, obl_const, per_const
  character (len=256) :: orbit_file

  integer :: isol
  real(wp) :: sol_const
  character (len=256) :: sol_file

  integer :: ivolc
  real(wp) :: volc_const
  character (len=256) :: volc_file
  real(wp) :: volc_scale

  integer :: isea_level
  real(wp) :: sea_level_const
  real(wp) :: sea_level_init
  character (len=256) :: sea_level_file

  integer :: ico2, id13c, iD14c  
  real(wp) :: co2_const, d13c_atm_const, D14c_atm_const 
  character (len=256) :: co2_file, d13c_atm_file, D14c_atm_file
  real(wp) :: dco2_dt
  real(wp) :: co2_max

  integer :: ico2_rad
  real(wp) :: co2_ref
  real(wp) :: co2_rad_const
  character (len=256) :: co2_rad_file

  integer :: iC14_production
  real(wp) :: C14_production_const
  character (len=256) :: C14_production_file

  integer :: ico2_degas
  real(wp) :: co2_degas_const, d13c_degas
  character (len=256) :: co2_degas_file

  logical :: l_weathering
  real(wp) :: d13c_weath

  integer :: ico2_emis
  real(wp) :: co2_emis_const
  character (len=256) :: co2_emis_file
  real(wp) :: co2_pulse
  real(wp) :: k_emis_fb
  real(wp) :: C_emis_fb
  real(wp) :: co2_emis_min

  integer :: id13C_emis 
  real(wp) :: d13C_emis_const 
  character (len=256) :: d13C_emis_file

  logical :: l_c13
  logical :: l_c14
  logical :: l_ocn_co2
  logical :: l_lnd_co2

  integer :: ich4
  real(wp) :: ch4_const
  integer :: i_ch4_tau
  real(wp) :: ch4_tau_const
  character (len=256) :: ch4_file
  character (len=256) :: ch4_tau_file
  character (len=256) :: CO_NOx_VOC_file

  integer :: ich4_rad
  real(wp) :: ch4_ref
  real(wp) :: ch4_rad_const
  character (len=256) :: ch4_rad_file

  integer :: ich4_emis
  real(wp) :: ch4_emis_const
  character (len=256) :: ch4_emis_file
  real(wp) :: ch4_emis_other_const

  integer :: in2o
  real(wp) :: n2o_ref
  real(wp) :: n2o_const
  character (len=256) :: n2o_file

  integer :: iso4
  real(wp) :: so4_const
  character (len=256) :: so4_file

  integer :: io3
  real(wp) :: o3_const
  character (len=256) :: o3_file_const
  character (len=256) :: o3_file_var

  integer :: icfc
  real(wp) :: cfc11_const
  real(wp) :: cfc12_const
  character (len=256) :: cfc_file

  integer :: iluc
  character (len=256) :: luc_file

  integer :: idist
  character (len=256) :: dist_file

  character (len=256) :: in_dir, out_dir

  integer :: ifake_atm
  character (len=256) :: fake_atm_const_file, fake_atm_var_file
  integer :: prc_forcing
  integer :: ifake_dust
  character (len=256) :: fake_dust_const_file, fake_dust_var_file
  integer :: ifake_ocn
  character (len=256) :: fake_ocn_const_file, fake_ocn_var_file
  integer :: ifake_sic
  character (len=256) :: fake_sic_const_file, fake_sic_var_file
  character (len=256) :: fake_lnd_const_file
  integer :: ifake_ice
  character (len=256) :: fake_ice_const_file, fake_ice_var_file
  integer :: ifake_geo
  character (len=256) :: fake_geo_const_file, fake_geo_var_file, fake_geo_ref_file

  logical :: check_energy, check_water, check_carbon

  logical :: l_debug_main_loop
  logical :: l_write_timer

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  c o n t r o l _ l o a d
  !   Purpose    :  read climber control parameters
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine control_load(filename)

    implicit none

    character (len=*), intent(in) :: filename

    integer :: n

    ! Read parameters from file
    write(*,*) "control parameters ==========="
    call nml_read(filename,"control","flag_co2",flag_co2)
    call nml_read(filename,"control","flag_ch4",flag_ch4)
    call nml_read(filename,"control","flag_atm",flag_atm)
    call nml_read(filename,"control","flag_lnd",flag_lnd)
    call nml_read(filename,"control","flag_dust",flag_dust)
    call nml_read(filename,"control","flag_lakes",flag_lakes)
    call nml_read(filename,"control","flag_ocn",flag_ocn)
    call nml_read(filename,"control","flag_sic",flag_sic)
    call nml_read(filename,"control","flag_bgc",flag_bgc)
    call nml_read(filename,"control","flag_ice",flag_ice)
    call nml_read(filename,"control","ice_model_name",ice_model_name)
    call nml_read(filename,"control","ice_domain_name",ice_domain_name)
    n_ice_domain = 0
    do n=1,10
      if (len_trim(ice_domain_name(n)).ne.0) n_ice_domain=n_ice_domain+1
    enddo
    if (n_ice_domain>n_ice_domain_max) then
      print *,'ERROR: too many ice domains, n_ice_domain>n_ice_domain_max',n_ice_domain,n_ice_domain_max
      stop
    endif
    print *,'number of ice domains',n_ice_domain
    call nml_read(filename,"control","flag_geo",flag_geo)
    call nml_read(filename,"control","l_aquaplanet",l_aquaplanet)
    call nml_read(filename,"control","l_aqua_slab",l_aqua_slab)
    call nml_read(filename,"control","l_feedbacks",l_feedbacks)
    call nml_read(filename,"control","l_spinup_cc",l_spinup_cc)
    if (l_spinup_cc) then
      l_daily_input_save_ocn = .true.
      l_daily_input_save_bgc = .true.
    else
      l_daily_input_save_ocn = .false.
      l_daily_input_save_bgc = .false.
    endif
    call nml_read(filename,"control","nyears_spinup_bgc",nyears_spinup_bgc)
    call nml_read(filename,"control","year_start_offline",year_start_offline)
    call nml_read(filename,"control","nyear_avg_offline",nyear_avg_offline)
    call nml_read(filename,"control","flag_smb",flag_smb)
    call nml_read(filename,"control","flag_bmb",flag_bmb)
    call nml_read(filename,"control","ocn_restore_sal",ocn_restore_sal)
    call nml_read(filename,"control","ocn_restore_temp",ocn_restore_temp)
    call nml_read(filename,"control","atm_fix_tau",atm_fix_tau)

    call nml_read(filename,"control","iorbit",iorbit)
    call nml_read(filename,"control","ecc_const",ecc_const)
    call nml_read(filename,"control","obl_const",obl_const)
    call nml_read(filename,"control","per_const",per_const)
    call nml_read(filename,"control","orbit_file",orbit_file)
    call nml_read(filename,"control","isol",isol)
    call nml_read(filename,"control","sol_const",sol_const)
    call nml_read(filename,"control","sol_file",sol_file)
    call nml_read(filename,"control","ivolc",ivolc)
    call nml_read(filename,"control","volc_const",volc_const)
    call nml_read(filename,"control","volc_file",volc_file)
    call nml_read(filename,"control","volc_scale",volc_scale)
    call nml_read(filename,"control","isea_level",isea_level)
    call nml_read(filename,"control","sea_level_const",sea_level_const)
    call nml_read(filename,"control","sea_level_init",sea_level_init)
    call nml_read(filename,"control","ico2",ico2)
    call nml_read(filename,"control","co2_ref",co2_ref)
    call nml_read(filename,"control","co2_const",co2_const)
    call nml_read(filename,"control","co2_file",co2_file)
    call nml_read(filename,"control","dco2_dt",dco2_dt)
    call nml_read(filename,"control","co2_max",co2_max)
    call nml_read(filename,"control","ico2_rad",ico2_rad)
    call nml_read(filename,"control","co2_rad_const",co2_rad_const)
    call nml_read(filename,"control","co2_rad_file",co2_rad_file)
    call nml_read(filename,"control","id13c",id13c)
    call nml_read(filename,"control","d13c_atm_const",d13c_atm_const)
    call nml_read(filename,"control","d13c_atm_file",d13c_atm_file)
    call nml_read(filename,"control","iD14c",iD14c)
    call nml_read(filename,"control","D14c_atm_const",D14c_atm_const)    
    call nml_read(filename,"control","D14c_atm_file",D14c_atm_file)
    call nml_read(filename,"control","iC14_production",iC14_production)
    call nml_read(filename,"control","C14_production_const",C14_production_const)
    call nml_read(filename,"control","C14_production_file",C14_production_file)
    call nml_read(filename,"control","ico2_degas",ico2_degas)
    call nml_read(filename,"control","co2_degas_const",co2_degas_const)
    call nml_read(filename,"control","co2_degas_file",co2_degas_file)
    call nml_read(filename,"control","d13c_degas",d13c_degas)
    call nml_read(filename,"control","l_weathering",l_weathering)
    call nml_read(filename,"control","d13c_weath",d13c_weath)
    call nml_read(filename,"control","ico2_emis",ico2_emis)
    call nml_read(filename,"control","co2_emis_const",co2_emis_const)
    call nml_read(filename,"control","co2_emis_file",co2_emis_file)
    call nml_read(filename,"control","co2_pulse",co2_pulse)
    call nml_read(filename,"control","k_emis_fb",k_emis_fb)
    call nml_read(filename,"control","C_emis_fb",C_emis_fb)
    call nml_read(filename,"control","co2_emis_min",co2_emis_min)
    call nml_read(filename,"control","id13C_emis",id13C_emis)
    call nml_read(filename,"control","d13C_emis_const",d13C_emis_const)
    call nml_read(filename,"control","d13C_emis_file",d13C_emis_file)
    call nml_read(filename,"control","l_c13",l_c13)
    call nml_read(filename,"control","l_c14",l_c14)
    call nml_read(filename,"control","l_ocn_co2",l_ocn_co2)
    call nml_read(filename,"control","l_lnd_co2",l_lnd_co2)
    call nml_read(filename,"control","ich4",ich4)
    call nml_read(filename,"control","ch4_const",ch4_const)
    call nml_read(filename,"control","i_ch4_tau",i_ch4_tau)
    call nml_read(filename,"control","ch4_tau_const",ch4_tau_const)
    call nml_read(filename,"control","ch4_tau_file",ch4_tau_file)
    call nml_read(filename,"control","CO_NOx_VOC_file",CO_NOx_VOC_file)
    call nml_read(filename,"control","ch4_file",ch4_file)
    call nml_read(filename,"control","ich4_rad",ich4_rad)
    call nml_read(filename,"control","ch4_ref",ch4_ref)
    call nml_read(filename,"control","ch4_rad_const",ch4_rad_const)
    call nml_read(filename,"control","ch4_rad_file",ch4_rad_file)
    call nml_read(filename,"control","ich4_emis",ich4_emis)
    call nml_read(filename,"control","ch4_emis_const",ch4_emis_const)
    call nml_read(filename,"control","ch4_emis_file",ch4_emis_file)
    call nml_read(filename,"control","ch4_emis_other_const",ch4_emis_other_const)
    call nml_read(filename,"control","in2o",in2o)
    call nml_read(filename,"control","n2o_ref",n2o_ref)
    call nml_read(filename,"control","n2o_const",n2o_const)
    call nml_read(filename,"control","n2o_file",n2o_file)
    call nml_read(filename,"control","iso4",iso4)
    call nml_read(filename,"control","so4_const",so4_const)
    call nml_read(filename,"control","so4_file",so4_file)
    call nml_read(filename,"control","io3",io3)
    call nml_read(filename,"control","o3_const",o3_const)
    call nml_read(filename,"control","o3_file_const",o3_file_const)
    call nml_read(filename,"control","o3_file_var",o3_file_var)
    call nml_read(filename,"control","icfc",icfc)
    call nml_read(filename,"control","cfc11_const",cfc11_const)
    call nml_read(filename,"control","cfc12_const",cfc12_const)
    call nml_read(filename,"control","cfc_file",cfc_file)
    call nml_read(filename,"control","iluc",iluc)
    call nml_read(filename,"control","luc_file",luc_file)
    call nml_read(filename,"control","idist",idist)
    call nml_read(filename,"control","dist_file",dist_file)
    call nml_read(filename,"control","atm_restart",atm_restart)
    call nml_read(filename,"control","co2_restart",co2_restart)
    call nml_read(filename,"control","ch4_restart",ch4_restart)
    call nml_read(filename,"control","lnd_restart",lnd_restart)
    call nml_read(filename,"control","ocn_restart",ocn_restart)
    call nml_read(filename,"control","bgc_restart",bgc_restart)
    call nml_read(filename,"control","sic_restart",sic_restart)
    call nml_read(filename,"control","geo_restart",geo_restart)
    call nml_read(filename,"control","ice_restart",ice_restart)
    call nml_read(filename,"control","smb_restart",smb_restart)
    call nml_read(filename,"control","bmb_restart",bmb_restart)
    call nml_read(filename,"control","restart_in_dir",restart_in_dir)
    call nml_read(filename,"control","i_write_restart",i_write_restart)
    call nml_read(filename,"control","n_year_write_restart",n_year_write_restart)
    years_write_restart(:) = -999999
    call nml_read(filename,"control","years_write_restart",years_write_restart)
    call nml_read(filename,"control","ifake_atm",ifake_atm)
    call nml_read(filename,"control","fake_atm_const_file",fake_atm_const_file)
    call nml_read(filename,"control","fake_atm_var_file",fake_atm_var_file)
    call nml_read(filename,"control","ifake_dust",ifake_dust)
    call nml_read(filename,"control","fake_dust_const_file",fake_dust_const_file)
    call nml_read(filename,"control","fake_dust_var_file",fake_dust_var_file)
    call nml_read(filename,"control","ifake_ocn",ifake_ocn)
    call nml_read(filename,"control","fake_ocn_const_file",fake_ocn_const_file)
    call nml_read(filename,"control","fake_ocn_var_file",fake_ocn_var_file)
    call nml_read(filename,"control","fake_lnd_const_file",fake_lnd_const_file)
    call nml_read(filename,"control","ifake_sic",ifake_sic)
    call nml_read(filename,"control","fake_sic_const_file",fake_sic_const_file)
    call nml_read(filename,"control","fake_sic_var_file",fake_sic_var_file)
    call nml_read(filename,"control","ifake_ice",ifake_ice)
    call nml_read(filename,"control","fake_ice_const_file",fake_ice_const_file)
    call nml_read(filename,"control","fake_ice_var_file",fake_ice_var_file)
    call nml_read(filename,"control","ifake_geo",ifake_geo)
    call nml_read(filename,"control","fake_geo_const_file",fake_geo_const_file)
    call nml_read(filename,"control","fake_geo_var_file",fake_geo_var_file)
    call nml_read(filename,"control","fake_geo_ref_file",fake_geo_ref_file)
 
    call nml_read(filename,"control","prc_forcing",prc_forcing)

    call nml_read(filename,"control","sea_level_file",sea_level_file)
    
    call nml_read(filename,"control","check_water",check_water)
    call nml_read(filename,"control","check_energy",check_energy)
    call nml_read(filename,"control","check_carbon",check_carbon)

    call nml_read(filename,"control","l_debug_main_loop",l_debug_main_loop)
    call nml_read(filename,"control","l_write_timer",l_write_timer)

   return

  end subroutine control_load

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  a r g s
  ! Purpose    :  Get command line arguments from program call
  !               *only obtains character arguments
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  subroutine args(folder)

    implicit none

    integer :: i, narg

    character (len=256)  :: string
    character (len=*)    :: folder

    ! Get the number of arguments total
    narg = command_argument_count()

    ! Set default values (in case no arguments provided)
    ! Default value assumes program is running within
    ! current folder.
    folder = "./"

    if (narg .gt. 0) then
      ! Get the last argument, the rest are ignored
      ! note: this should be arg #1 (after the executable)
      call get_command_argument(narg,string)
      folder = trim(adjustl(string))
    end if

    i = len(trim(folder))
    if ( scan(trim(folder),"/",back=.TRUE.) .ne. i ) folder = trim(folder)//"/"

    write(*,"(a1,5x,a,a)") "e","output folder: ", trim(folder)

   return

  end subroutine args


end module control
