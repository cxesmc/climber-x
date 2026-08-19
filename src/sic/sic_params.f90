!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s i c _ p a r a m s
!
!  Purpose : sea ice parameters
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
module sic_params

  use precision, only : wp
  use nml, only : nml_read
  use timer, only : sec_day, dt_sic
  use control, only : out_dir

  implicit none

  real(wp) :: dt
  real(wp) :: f_sic_min
  real(wp) :: h_sic_min
  real(wp) :: h_snow_min
  real(wp) :: h_sic_max
  real(wp) :: h_snow_max
  real(wp) :: diffsic
  real(wp) :: i_ice_albedo
  real(wp) :: alb_sic_vis
  real(wp) :: alb_sic_nir
  real(wp) :: c_alb_sic
  integer :: i_adv
  integer :: i_fsic
  real(wp) :: f_sic_max
  real(wp) :: f_sic_pow
  real(wp) :: h0
  ! regional sea-ice prescription (sensitivity experiment)
  logical :: l_sic_prescribe
  real(wp) :: sic_prescribe_lon_w
  real(wp) :: sic_prescribe_lon_e
  real(wp) :: sic_prescribe_lat_s
  real(wp) :: sic_prescribe_lat_n
  real(wp) :: sic_prescribe_fsic
  real(wp) :: sic_prescribe_h_sic
  real(wp) :: sic_prescribe_h_snow
  real(wp) :: u_star
  integer :: i_cd_ocn
  logical :: l_neutral_ocn
  logical :: l_neutral_sic
  real(wp) :: z0m_ocn   !! ocean water roughness length for momentum (m)
  real(wp) :: z0m_sic   !! sea ice roughness length for momentum (m)
  real(wp), parameter :: z0m_snow = 0.0024_wp !! snow roughness_length for momentum (m)
  real(wp) :: Cde0
  real(wp) :: Cdh0
  real(wp) :: f_Ri_stab
  real(wp) :: f_Ri_unstab
  logical :: l_diag_dyn
  logical :: l_daily_output 

  ! parameters for sea ice dynamics
  type :: dyn_par_type
    real(wp) :: p0 = 2.75e4             ! A constant in the expression for the ice strength, P* in Hunke & Dukowicz 1997.[Pa]
    real(wp) :: p0_rho                  ! The pressure constant divided by ice density [N m kg-1].
    real(wp) :: c0 = 20.                ! A constant in the exponent of the expression for the ice strength, c* in Hunke & Dukowicz 1997.
    real(wp) :: cdw = 3.24e-3           ! The drag coefficient between the sea ice and water. [nondim]
    real(wp) :: EC = 2.0                ! The ellipticity coefficient for the plastic yield curve in the sea-ice rheology.  
                                        ! For an infinite ellipticity (i.e., a cavitating fluid rheology), use 0. units="Nondim"
    real(wp) :: Rho_ocean = 1030.       ! The nominal density of sea water [kg m-3].
    real(wp) :: Rho_ice = 905.          ! The nominal density of sea ice [kg m-3].
    real(wp) :: drag_bg_vel2 = 0.0      ! A background (subgridscale) velocity for drag with the ocean squared [m2 s-2].  This is always 0 for now.
    real(wp) :: Tdamp = -0.2            ! The damping timescale associated with the elastic terms in the sea-ice dynamics equations (if positive) 
                                        ! or the fraction of DT_ICE_DYNAMICS (if negative). units="s or nondim"
    real(wp) :: del_sh_min_scale = 2.0  ! A scaling factor for the minimum permitted value of minimum shears 
                                        ! used in the denominator of the stress equations [nondim]. Probably needs to be greater than 1.
    logical :: project_drag_vel = .true.! If true, project forward the ice velocity used in the drag calculation to avoid an instability that can occur 
                                        ! when an finite stress is applied to thin ice moving with the velocity of the ocean.
    logical :: project_f_sic = .true.      ! If true, project the ice concentration and related ice strength changes due to the convergent or divergent ice flow.
    integer :: evp_sub_steps = 432      ! The number of iterations in the EVP dynamics for each slow time step.
    real(wp) :: cfl_fac = 0.5_wp 
  end type dyn_par_type
  type(dyn_par_type) :: dyn_par

  ! snow parameters
  type snow_par_type
    integer :: i_snow_albedo
    logical :: l_snow_aging
    logical :: l_snow_dust
    real(wp) :: alb_snow_vis_dif_new = 0.99_wp
    real(wp) :: alb_snow_nir_dif_new = 0.65_wp
    real(wp) :: snow_grain_fresh
    real(wp) :: snow_grain_old
    real(wp) :: f_age_t 
    real(wp) :: snow_0
    real(wp) :: snow_1
  end type
  type(snow_par_type) :: snow_par

  real(wp), parameter :: rho_sic  = 910._wp      ! kg/m3, Density of sea ice
  real(wp), parameter :: rho_snow = 250._wp      ! kg/m3, Density of snow on sea ice
  real(wp), parameter :: emis_sic = 0.99_wp
  real(wp) :: lambda_sic
  real(wp) :: lambda_snow
  real(wp) :: h_k_crit
  

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s i c _ p a r a m s _ i n i t
  !   Purpose    :  reading/initialisation of sea ice parameters
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_params_init

    implicit none


    ! time step
    dt = dt_sic

    call sic_par_load(trim(out_dir)//"/sic_par.nml")


    return

  end subroutine sic_params_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s i c _ p a r _ l o a d 
  !   Purpose    :  reading of sea ice parameters namelist
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_par_load(filename)

    implicit none

    character (len=*) :: filename


    ! Read parameters from file
    write(*,*) "sea ice parameters ==========="
    call nml_read(filename,"sic_par","diffsic",diffsic)
    call nml_read(filename,"sic_par","f_sic_min",f_sic_min)
    call nml_read(filename,"sic_par","h_sic_min",h_sic_min)
    call nml_read(filename,"sic_par","h_snow_min",h_snow_min)
    call nml_read(filename,"sic_par","h_sic_max",h_sic_max)
    call nml_read(filename,"sic_par","h_snow_max",h_snow_max)
    call nml_read(filename,"sic_par","i_ice_albedo",i_ice_albedo)
    call nml_read(filename,"sic_par","alb_sic_vis",alb_sic_vis)
    call nml_read(filename,"sic_par","alb_sic_nir",alb_sic_nir)
    call nml_read(filename,"sic_par","c_alb_sic",c_alb_sic)
    call nml_read(filename,"sic_par","i_fsic",i_fsic)
    call nml_read(filename,"sic_par","i_adv",i_adv)
    call nml_read(filename,"sic_par","f_sic_max",f_sic_max)
    call nml_read(filename,"sic_par","f_sic_pow",f_sic_pow)
    call nml_read(filename,"sic_par","h0",h0)
    call nml_read(filename,"sic_par","l_sic_prescribe",l_sic_prescribe)
    call nml_read(filename,"sic_par","sic_prescribe_lon_w",sic_prescribe_lon_w)
    call nml_read(filename,"sic_par","sic_prescribe_lon_e",sic_prescribe_lon_e)
    call nml_read(filename,"sic_par","sic_prescribe_lat_s",sic_prescribe_lat_s)
    call nml_read(filename,"sic_par","sic_prescribe_lat_n",sic_prescribe_lat_n)
    call nml_read(filename,"sic_par","sic_prescribe_fsic",sic_prescribe_fsic)
    call nml_read(filename,"sic_par","sic_prescribe_h_sic",sic_prescribe_h_sic)
    call nml_read(filename,"sic_par","sic_prescribe_h_snow",sic_prescribe_h_snow)
    call nml_read(filename,"sic_par","lambda_snow",lambda_snow)
    call nml_read(filename,"sic_par","lambda_sic",lambda_sic)
    call nml_read(filename,"sic_par","h_k_crit",h_k_crit)
    call nml_read(filename,"sic_par","u_star",u_star)
    call nml_read(filename,"sic_par","i_cd_ocn",i_cd_ocn)
    call nml_read(filename,"sic_par","l_neutral_ocn",l_neutral_ocn)
    call nml_read(filename,"sic_par","l_neutral_sic",l_neutral_sic)
    call nml_read(filename,"sic_par","z0m_ocn",z0m_ocn)
    call nml_read(filename,"sic_par","z0m_sic",z0m_sic)
    call nml_read(filename,"sic_par","Cde0",Cde0)
    call nml_read(filename,"sic_par","Cdh0",Cdh0)
    call nml_read(filename,"sic_par","f_Ri_stab",f_Ri_stab)
    call nml_read(filename,"sic_par","f_Ri_unstab",f_Ri_unstab)
    call nml_read(filename,"sic_par","l_daily_output",l_daily_output)

    call nml_read(filename,"sic_par","i_snow_albedo",snow_par%i_snow_albedo)
    call nml_read(filename,"sic_par","l_snow_aging",snow_par%l_snow_aging)
    call nml_read(filename,"sic_par","l_snow_dust",snow_par%l_snow_dust)
    call nml_read(filename,"sic_par","snow_grain_fresh" ,snow_par%snow_grain_fresh )
    call nml_read(filename,"sic_par","snow_grain_old" ,snow_par%snow_grain_old )
    call nml_read(filename,"sic_par","f_age_t" ,snow_par%f_age_t)
    call nml_read(filename,"sic_par","snow_0" ,snow_par%snow_0 )
    snow_par%snow_0 = snow_par%snow_0/sec_day
    call nml_read(filename,"sic_par","snow_1" ,snow_par%snow_1 )

    call nml_read(filename,"sic_par","p0              ", dyn_par%p0)               
    call nml_read(filename,"sic_par","c0              ", dyn_par%c0) 
    call nml_read(filename,"sic_par","cdw             ", dyn_par%cdw) 
    call nml_read(filename,"sic_par","EC              ", dyn_par%EC) 
    call nml_read(filename,"sic_par","Tdamp           ", dyn_par%Tdamp) 
    call nml_read(filename,"sic_par","del_sh_min_scale", dyn_par%del_sh_min_scale)
    call nml_read(filename,"sic_par","project_drag_vel", dyn_par%project_drag_vel)
    call nml_read(filename,"sic_par","project_f_sic   ", dyn_par%project_f_sic) 
    call nml_read(filename,"sic_par","evp_sub_steps   ", dyn_par%evp_sub_steps)
    call nml_read(filename,"sic_par","cfl_fac         ", dyn_par%cfl_fac)
    dyn_par%p0_rho = dyn_par%p0 / dyn_par%Rho_ice

    return

  end subroutine sic_par_load

end module sic_params
