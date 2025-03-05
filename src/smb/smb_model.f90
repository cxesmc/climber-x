!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s m b _ m o d e l 
!
!  Purpose : Main surface mass balance model 
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2024 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit, Reihard Calov and Andrey Ganopolski
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
module smb_model

  use nml
  use ncio
  use precision, only : wp, dp
  use timer, only : year, year_ini, year_now, doy, mon, nmon_year, nday_year, day_year, nstep_mon_smb, nstep_year_smb, sec_year, sec_mon
  use timer, only : time_soy_smb, time_eoy_smb, time_eom_smb
  use constants, only : pi, rho_i, fcoriolis, T0, r_earth
  use control, only : out_dir, restart_in_dir, smb_restart
  use smb_grid, only : smb_grid_init, nl
  use smb_params, only : i_smb, gamma, i_z_sur_eff, alpha_zstd, l_regional_climate_forcing, l_diurnal_cycle, prc_par
  use smb_params, only : smb_params_init, map_method, filt_sigma, dt, snow_par, surf_par, smb_crit_mask, n_smb_mask_ext, nday_update_climate
  use smb_params, only : l_smb_bias_corr, year_ini_smb_ref, year_end_smb_ref
  use smb_params, only : i_t2m_bias_corr, i_prc_bias_corr,  bias_corr_file, t2m_bias_scale_fac, t2m_bias_corr_uniform
  use smb_params, only : l_maxice, mask_maxice_file, smb_maxice
  use smb_params, only : l_Tvar_ann, Tvar_ann_amp, Tvar_ann_period
  use smb_params, only : l_Tvar_day, Tvar_day_amp, Tvar_day_period
  use smb_params, only : l_write_timer
  use smb_def, only : smb_in_class, smb_class

  use fake_atm_hires_mod, only : fake_atm_hires_type, fake_atm_hires_init, fake_atm_hires_update
  use coord, only : grid_allocate, grid_init, grid_class
  use coord, only : map_scrip_class, map_scrip_init, map_scrip_field
  use filter, only : filter1d

  use topo_mod, only : topo_filter, topo_grad_map1, topo_grad_map2, topo_factors
  use ice_mod, only : frac_ice, albedo_ice, margin_ice
  use smb_simple_m, only : smb_simple
  use smb_pdd_m, only : smb_pdd
  use semi_m, only : semi
  use bias_corr_mod, only : smb_ref_write, smb_bias_corr, t2m_bias_corr, prc_bias_corr

  implicit none

  type(fake_atm_hires_type) :: atm        ! high resolution atmosphere forcing type (from regional climate models)

  private
  public :: smb_init, smb_update, smb_end, smb_write_restart

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s m b _ u p d a t e
  !   Purpose    :  update surface energy and mass balance interface
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_update(smb_in,smb)

  !$  use omp_lib

  implicit none

  type(smb_in_class), intent(in) :: smb_in
  type(smb_class) :: smb

  integer :: i, j, ii, jj, n, nn, k
  character(len=256) :: fnm
  integer :: ncid
  !$ real(wp) :: time1,time2


  !-------------------------------------
  ! initialize 
  !-------------------------------------

  if (time_soy_smb) then
    ! initialize simple smb cumulated values at start of year
    smb%simple%pdd      = 0.0_wp
    smb%simple%snow_cum = 0.0_wp
    smb%simple%rain_cum = 0.0_wp
    smb%simple%t2m_cum  = 0.0_wp

    ! initalize annual values
    smb%ann_prc       = 0._wp 
    smb%ann_snow      = 0._wp 
    smb%ann_smb       = 0._wp 
    smb%ann_melt      = 0._wp 
    smb%ann_icemelt   = 0._wp 
    smb%ann_ablation  = 0._wp 
    smb%ann_evp       = 0._wp 
    smb%ann_runoff    = 0._wp 
    smb%ann_refreezing= 0._wp 
    smb%dt_snowfree   = 0._wp
  endif  

  !-------------------------------------
  ! topography related quantities
  !-------------------------------------

  ! deal with topography related quantities at start of year
  if (time_soy_smb) then

    ! interpolate low resolution elevation
    call map_scrip_field(smb%maps_cmn_to_ice,"z_sur",smb_in%z_sur, smb%z_sur_i,method="mean",missing_value=-9999._dp)

    ! filter topography for precipitation downscaling (orographic enhancement of precipitation follows a smoothed topography (Pedgley 1970))
    call topo_filter(smb%grid, smb%z_sur, smb%z_sur_fil)

    ! compute topography gradients of filtered topography
    if (i_smb.eq.3 .or. prc_par%l_slope_effect) then
      !call topo_grad_map2(smb%grid, smb%grid_latlon, smb%maps_to_latlon, smb%maps_from_latlon, smb%z_sur-smb%z_sur_i, &
      call topo_grad_map2(smb%grid, smb%grid_latlon, smb%maps_to_latlon, smb%maps_from_latlon, smb%z_sur_fil, &
        smb%dz_sur, smb%dz_dx_sur, smb%dz_dy_sur)
    else
      smb%dz_sur    = 0._wp
      smb%dz_dx_sur = 0._wp
      smb%dz_dy_sur = 0._wp
    endif

    ! compute topography related factors
    call topo_factors(smb%grid, smb%z_sur, smb%z_sur_i, &
      smb%f_ele, smb%pressure)

    ! derive ice margin mask
    call margin_ice(smb%mask_ice, &
      smb%mask_margin)

    ! sub-grid standard deviation of topography
    smb%z_sur_std = max(0._wp,smb%z_bed_std-smb%h_ice)*min(smb%z_sur,1._wp)

    ! effective surface elevation 
    if (i_z_sur_eff.eq.0) then
      ! equal to actual elevation
      smb%z_sur_eff = smb%z_sur
    else if (i_z_sur_eff.eq.1) then
      ! accounting for sub-grid standard deviation of topography
      ! suppress elevation correction when ice is growing thick relative to the std and over ocean
      smb%z_sur_eff = smb%z_sur + alpha_zstd*smb%z_sur_std
    endif

  endif

  !-------------------------------------
  ! ice fraction and albedo 
  !-------------------------------------

  if (time_soy_smb) then

    ! save old ice cover fraction
    smb%f_ice_old = smb%f_ice
    ! ice sheet fraction 
    call frac_ice(smb%mask_ice, smb%mask_margin, smb%h_ice, smb%z_bed_std*min(smb%z_sur,1._wp), &   ! in
      smb%f_ice)    ! out

  endif

  if (time_eoy_smb) then

    call albedo_ice(smb%mask_ice, smb%mask_margin, smb%f_ice, smb%f_ice_old, smb%ann_icemelt, smb%ann_smb, smb%dt_snowfree, &  ! in
      smb%alb_ice)  ! inout

  endif

  !-------------------------------------
  ! generate artificial temperature variability
  !-------------------------------------
  if (time_soy_smb .and. l_Tvar_ann) then
    if (mod(year,Tvar_ann_period).eq.0) then
      ! change sign every Tvar_ann_period years
      smb%dTvar_ann = -smb%dTvar_ann
    endif
  endif
  if (l_Tvar_day) then
    if (mod(doy,Tvar_day_period).eq.0) then
      ! change sign every Tvar_day_period days
      smb%dTvar_day = -smb%dTvar_day
    endif
  endif
  smb%dTvar = smb%dTvar_ann + smb%dTvar_day

  !-------------------------------------
  ! smb PDD model
  !-------------------------------------

  ! call smb PDD at the end of month
  ! smb PDD is called in any case, either it is used as model for the
  ! surface mass balance or it is just used to define the mask where the
  ! physically-based surface energy and mass balance model (SEMI) is applied
  if (i_smb.ne.3) then

    if (time_eom_smb) then

      !-------------------------------------
      ! map required fields from low to high resolution grid

      !$omp parallel sections
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"t2m",smb_in%t2m, smb%t2m_i,method="mean",missing_value=-9999._dp)
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"prc",smb_in%prc, smb%prc_i,method="mean",missing_value=-9999._dp)
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"u700",smb_in%u700, smb%u700_i,method="mean",missing_value=-9999._dp)
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"v700",smb_in%v700, smb%v700_i,method="mean",missing_value=-9999._dp)
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"wind",smb_in%wind, smb%wind_i,method="mean",missing_value=-9999._dp)
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"cld",smb_in%cld, smb%cld_i,method="mean",missing_value=-9999._dp)
      !$omp end parallel sections

      !$omp parallel do private (i,j) 
      do i = 1,smb%grid%G%nx
        do j = 1,smb%grid%G%ny

          !-------------------------------------
          ! call smb PDD
          call smb_pdd(smb%z_sur_eff(i,j), smb%z_sur_i(i,j), smb%dz_dx_sur(i,j), smb%dz_dy_sur(i,j), smb%dz_sur(i,j), smb%f_ele(i,j), &   ! in
            smb%t2m_i(i,j), smb%t2m_bias_i(i,j,doy), smb%u700_i(i,j), smb%v700_i(i,j), smb%wind_i(i,j), smb%prc_i(i,j), smb%prc_bias_i(i,j,doy), &    ! in
            smb%t2m(i,j), smb%u700(i,j), smb%v700(i,j), smb%wind(i,j), smb%snow(i,j), smb%rain(i,j), smb%prc(i,j), smb%f_wind(i,j), &   ! out
            smb%simple%pdd(i,j), smb%simple%t2m_cum(i,j), smb%simple%snow_cum(i,j), smb%simple%rain_cum(i,j),&    ! inout
            smb%simple%smb(i,j), smb%simple%melt(i,j), smb%simple%melt_star(i,j), smb%simple%runoff(i,j), smb%simple%t_ice(i,j))  ! out at end of year

          ! if smb PDD is used for surface mass balance, assign relevant variables
          if (i_smb==2) then
            ! cumulate over the year
            smb%ann_prc(i,j)   = smb%ann_prc(i,j)  + smb%prc(i,j) *sec_mon    ! kg/m2
            smb%ann_snow(i,j)  = smb%ann_snow(i,j) + smb%snow(i,j)*sec_mon    ! kg/m2
            if (time_eoy_smb) then
              ! annually integrated values
              if (smb%mask_maxice(i,j).eq.1) then
                smb%ann_smb(i,j)       = smb%simple%smb(i,j)*sec_year    ! kg/m2
              else
                ! add large negative value outside of allowed ice mask
                smb%ann_smb(i,j)       = smb%simple%smb(i,j)*sec_year + smb_maxice
              endif
              smb%ann_melt(i,j)      = smb%simple%melt(i,j)*sec_year   ! kg/m2
              smb%ann_ablation(i,j)  = smb%simple%melt(i,j)*sec_year   ! kg/m2
              smb%ann_evp(i,j)       = 0._wp    ! kg/m2
              smb%ann_runoff(i,j)    = smb%simple%runoff(i,j)*sec_year   ! kg/m2
              smb%ann_refreezing(i,j)= smb%simple%melt_star(i,j)*sec_year   ! kg/m2 ! fixme?
              ! monthly runoff, same for all months
              smb%mon_runoff(i,j,:)      = smb%simple%runoff(i,j)   ! kg/m2/s
              ! monthly icemelt, same for all months, simply taken as half of the runoff
              smb%mon_icemelt(i,j,:)      = 0.5_wp*smb%simple%runoff(i,j)   ! kg/m2/s
              ! 10m firn temperature
              smb%t_ice(i,j)         = smb%simple%t_ice(i,j) ! degC
            endif
          endif

        enddo
      enddo
      !$omp end parallel do 

    endif

  endif


  !-------------------------------------
  ! SEMI 
  !-------------------------------------

  ! call actual semi model
  if (i_smb==1) then

    if (time_soy_smb) then

      ! update active cells
      n = 0  ! counter
      smb%ncells = 0
      smb%grid_smb_to_cmn%ncells_ice(:,:) = 0
      do i = 1,smb%grid%G%nx
        do j = 1,smb%grid%G%ny
          n = n + 1
          if (smb%mask_smb(i,j).eq.1) then
            smb%ncells = smb%ncells + 1 ! sum to get total number of active cells
            smb%idx_cell_active(smb%ncells) = n
          else
            ! reset variables outside of smb mask
            smb%t_skin(i,j)           = T0
            smb%t_prof(i,j,:)         = T0
            smb%snowmelt(i,j)         = 0._wp 
            smb%icemelt(i,j)          = 0._wp 
            smb%refreezing(i,j)       = 0._wp 
            smb%refreezing_sum(i,j)   = 0._wp 
            smb%runoff(i,j)           = 0._wp 
            smb%evp(i,j)              = 0._wp 
            smb%mask_snow(i,j)        = 0
            smb%w_snow(i,j)           = 0._wp
            smb%w_snow_max(i,j)       = 0._wp
            smb%h_snow(i,j)           = 0._wp 
            smb%snow_grain(i,j)       = snow_par%snow_grain_fresh
            smb%dust_con(i,j)         = 0._wp 
          endif
          ! count number of ice sheet cells in each climate model coarse grid cell
          if (smb%mask_ice(i,j).eq.1) then
            ii = smb%grid_smb_to_cmn%i_lowres(i,j) 
            jj = smb%grid_smb_to_cmn%j_lowres(i,j) 
            smb%grid_smb_to_cmn%ncells_ice(ii,jj) = smb%grid_smb_to_cmn%ncells_ice(ii,jj)+1
          endif
        enddo
      enddo

      ! initialize
      smb%mon_runoff = 0._wp
      smb%mon_icemelt = 0._wp
      smb%t_ice = 0._wp

    endif  

    !--------------------------------------------------------
    ! map fields from low (climate) to high (smb) resolution

    if (doy.eq.1 .or. mod(doy,nday_update_climate)==0) then

      !$ time1 = omp_get_wtime()

      !$omp parallel sections
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"tam",smb_in%tam, smb%tam_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"ram",smb_in%ram, smb%ram_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"gam",smb_in%gam, smb%gam_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"tstd",smb_in%tstd, smb%tstd_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"prc",smb_in%prc, smb%prc_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"u700",smb_in%u700, smb%u700_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"v700",smb_in%v700, smb%v700_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"wind",smb_in%wind, smb%wind_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"cld",smb_in%cld, smb%cld_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"dust",smb_in%dust, smb%dust_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      if (time_soy_smb) then
        call map_scrip_field(smb%maps_cmn_to_ice,"t_ground",smb_in%t_ground, smb%t_ground_i,method="mean",missing_value=-9999._dp)
      endif
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"swd_toa",smb_in%swd_toa, smb%swd_toa_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      if (l_diurnal_cycle) then
        call map_scrip_field(smb%maps_cmn_to_ice,"swd_toa_min",smb_in%swd_toa_min, smb%swd_toa_min_i,method="mean",missing_value=-9999._dp, &
          mask_pack=(smb%mask_smb==1))
      endif
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"swd_sur_vis_dir",smb_in%swd_sur_vis_dir, smb%swd_sur_vis_dir_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"swd_sur_nir_dir",smb_in%swd_sur_nir_dir, smb%swd_sur_nir_dir_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"swd_sur_vis_dif",smb_in%swd_sur_vis_dif, smb%swd_sur_vis_dif_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"swd_sur_nir_dif",smb_in%swd_sur_nir_dif, smb%swd_sur_nir_dif_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"dswd_dalb_vis_dir",smb_in%dswd_dalb_vis_dir, smb%dswd_dalb_vis_dir_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"dswd_dalb_nir_dir",smb_in%dswd_dalb_nir_dir, smb%dswd_dalb_nir_dir_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"dswd_dalb_vis_dif",smb_in%dswd_dalb_vis_dif, smb%dswd_dalb_vis_dif_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"dswd_dalb_nir_dif",smb_in%dswd_dalb_nir_dif, smb%dswd_dalb_nir_dif_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"dswd_dz_nir_dir",smb_in%dswd_dz_nir_dir, smb%dswd_dz_nir_dir_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"dswd_dz_nir_dif",smb_in%dswd_dz_nir_dif, smb%dswd_dz_nir_dif_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"alb_vis_dir",smb_in%alb_vis_dir, smb%alb_vis_dir_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"alb_nir_dir",smb_in%alb_nir_dir, smb%alb_nir_dir_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"alb_vis_dif",smb_in%alb_vis_dif, smb%alb_vis_dif_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"alb_nir_dif",smb_in%alb_nir_dif, smb%alb_nir_dif_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"coszm",smb_in%coszm, smb%coszm_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"lwdown",smb_in%lwdown, smb%lwdown_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp section
      !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()
      call map_scrip_field(smb%maps_cmn_to_ice,"gam_lw",smb_in%gam_lw, smb%gam_lw_i,method="mean",missing_value=-9999._dp, &
        mask_pack=(smb%mask_smb==1))
      !$omp end parallel sections

      ! total precipitation, todo
      !prc_i_total = sum(smb%prc_i)

      !$ time2 = omp_get_wtime()
      !$ if(l_write_timer) print *,'map',time2-time1

    endif

    if (l_regional_climate_forcing) then
      !---------------------------------------------------------------
      ! get high-resolution climate forcing from regional climate model

      call fake_atm_hires_update(real(year_now,wp),atm)
      call map_scrip_field(atm%maps_atm_to_smb,"t2m",atm%tair, smb%t2m,method="mean") 
      call map_scrip_field(atm%maps_atm_to_smb,"tstd",atm%tstd, smb%tstd_i,method="mean") 
      call map_scrip_field(atm%maps_atm_to_smb,"q2m",atm%qair, smb%q2m,method="mean") 
      call map_scrip_field(atm%maps_atm_to_smb,"rain",atm%rain, smb%rain,method="mean") 
      call map_scrip_field(atm%maps_atm_to_smb,"snow",atm%snow, smb%snow,method="mean") 
      call map_scrip_field(atm%maps_atm_to_smb,"swdown",atm%swdown, smb%swdown,method="mean") 
      call map_scrip_field(atm%maps_atm_to_smb,"lwdown",atm%lwdown, smb%lwdown,method="mean") 
      call map_scrip_field(atm%maps_atm_to_smb,"wind",atm%wind, smb%wind,method="mean") 
      call map_scrip_field(atm%maps_atm_to_smb,"cod",atm%cod, smb%cod,method="mean") 
      call map_scrip_field(atm%maps_atm_to_smb,"albedo",atm%alb, smb%albedo,method="mean") 
      call map_scrip_field(atm%maps_atm_to_smb,"pressure",atm%pressure, smb%pressure,method="mean") 

      !    fnm = "test_smb.nc"
      !    call nc_create(fnm)
      !    call nc_open(fnm,ncid)
      !    call nc_write_dim(fnm,"nx",x=1._wp,dx=1._wp,nx=smb%grid%G%nx,ncid=ncid)
      !    call nc_write_dim(fnm,"ny",x=1._wp,dx=1._wp,nx=smb%grid%G%ny,ncid=ncid)
      !    call nc_write(fnm,"t2m",      smb%t2m,dims=["nx","ny"],start=[1,1],count=[smb%grid%G%nx,smb%grid%G%ny],ncid=ncid)
      !    call nc_write(fnm,"q2m",     smb%q2m,dims=["nx","ny"],start=[1,1],count=[smb%grid%G%nx,smb%grid%G%ny],ncid=ncid)
      !    call nc_write(fnm,"rain",     smb%rain,dims=["nx","ny"],start=[1,1],count=[smb%grid%G%nx,smb%grid%G%ny],ncid=ncid)
      !    call nc_write(fnm,"snow",     smb%snow,dims=["nx","ny"],start=[1,1],count=[smb%grid%G%nx,smb%grid%G%ny],ncid=ncid)
      !    call nc_write(fnm,"lwdown",     smb%lwdown,dims=["nx","ny"],start=[1,1],count=[smb%grid%G%nx,smb%grid%G%ny],ncid=ncid)
      !    call nc_write(fnm,"swdown",     smb%swdown,dims=["nx","ny"],start=[1,1],count=[smb%grid%G%nx,smb%grid%G%ny],ncid=ncid)
      !    call nc_write(fnm,"wind",     smb%wind,dims=["nx","ny"],start=[1,1],count=[smb%grid%G%nx,smb%grid%G%ny],ncid=ncid)
      !    call nc_write(fnm,"pressure",     smb%pressure,dims=["nx","ny"],start=[1,1],count=[smb%grid%G%nx,smb%grid%G%ny],ncid=ncid)
      !    call nc_close(ncid)
    endif


    ! loop over grid points where surface energy and mass balance has to be computed (1D)

    !$ time1 = omp_get_wtime()
    !$omp parallel do  &
    !$omp private (i,j,n,nn) 
    do nn=1,smb%ncells

      n = smb%idx_cell_active(nn)
      i = smb%ij_1d(1, n)
      j = smb%ij_1d(2, n)
      !!$ print *,i,j,n,smb%ncells,smb%grid%G%nx*smb%grid%G%ny,omp_get_thread_num(),'/',omp_get_num_threads()

      !-------------------------------------
      ! call SEMI 
      call semi(i ,j ,smb%mask_ice(i,j), smb%mask_ice_old(i,j), smb%f_ice(i,j), smb%alb_ice(i,j), & ! in
        smb%z_sur_eff(i,j), smb%z_sur_i(i,j), smb%z_sur_std(i,j), &     ! in
        smb%dz_dx_sur(i,j), smb%dz_dy_sur(i,j), smb%dz_sur(i,j), smb%f_ele(i,j), & ! in
        smb%tam_i(i,j), smb%t2m_bias_i(i,j,doy), smb%dTvar, smb%gam_i(i,j), smb%tstd_i(i,j), smb%ram_i(i,j), &
        smb%pressure(i,j), smb%u700_i(i,j), smb%v700_i(i,j), smb%wind_i(i,j), smb%prc_i(i,j), smb%prc_bias_i(i,j,doy), &    ! in
        smb%alb_vis_dir_i(i,j), smb%alb_nir_dir_i(i,j), smb%alb_vis_dif_i(i,j), smb%alb_nir_dif_i(i,j), &     ! in
        smb%swd_sur_vis_dir_i(i,j), smb%swd_sur_nir_dir_i(i,j), smb%swd_sur_vis_dif_i(i,j), smb%swd_sur_nir_dif_i(i,j), &     ! in
        smb%dswd_dalb_vis_dir_i(i,j), smb%dswd_dalb_nir_dir_i(i,j), smb%dswd_dalb_vis_dif_i(i,j), smb%dswd_dalb_nir_dif_i(i,j), & ! in
        smb%dswd_dz_nir_dir_i(i,j), smb%dswd_dz_nir_dif_i(i,j), smb%dust_i(i,j), smb%coszm_i(i,j), &
        smb%swd_toa_i(i,j), smb%swd_toa_min_i(i,j), smb%cld_i(i,j), smb%lwdown_i(i,j), smb%gam_lw_i(i,j), &
        smb%tam(i,j), smb%t2m(i,j), smb%t_skin(i,j), smb%t_skin_old(i,j), smb%t_skin_amp(i,j), &  ! out
        smb%t_prof(i,j,:), smb%t_prof_old(i,j,:), & ! out
        smb%q2m(i,j), smb%qsat(i,j), smb%dqsatdT(i,j), smb%r_a(i,j), &     ! out
        smb%u700(i,j), smb%v700(i,j), smb%wind(i,j), smb%snow(i,j), smb%rain(i,j), smb%prc(i,j), smb%f_wind(i,j), &   ! out
        smb%mask_snow(i,j), smb%f_snow(i,j), smb%h_snow(i,j), smb%w_snow(i,j), smb%w_snow_old(i,j), smb%w_snow_max(i,j), &    ! out
        smb%snow_grain(i,j), smb%dust_con(i,j), & ! out
        smb%alb_vis_dir(i,j), smb%alb_nir_dir(i,j), smb%alb_vis_dif(i,j), smb%alb_nir_dif(i,j), &   ! out
        smb%alb_snow_vis_dir(i,j), smb%alb_snow_nir_dir(i,j), smb%alb_snow_vis_dif(i,j), smb%alb_snow_nir_dif(i,j), &   ! out
        smb%dt_snowfree(i,j), smb%alb_bg(i,j), smb%albedo(i,j), smb%cld(i,j), smb%swnet(i,j), smb%swnet_min(i,j), smb%swdown(i,j), smb%lwdown(i,j), &  ! out
        smb%flx_g(i,j), smb%dflxg_dT(i,j), smb%flx_melt(i,j), smb%flx_sh(i,j), smb%flx_lwu(i,j), smb%flx_lh(i,j), smb%evp(i,j), &   ! out
        smb%num_lh(i,j), smb%num_sh(i,j), smb%num_sw(i,j), smb%num_lw(i,j), smb%denom_lh(i,j), smb%denom_sh(i,j), smb%denom_lw(i,j), &  ! out
        smb%f_sh(i,j), smb%f_e(i,j), smb%f_lh(i,j), smb%f_lw(i,j), & ! out
        smb%snowmelt(i,j), smb%icemelt(i,j), smb%refreezing(i,j), smb%refreezing_sum(i,j), smb%f_rfz_to_snow(i,j))   ! out

      ! derive runoff
      smb%runoff(i,j) = smb%snowmelt(i,j) + smb%icemelt(i,j) + smb%rain(i,j) - smb%refreezing(i,j)
      ! cumulate runoff values for each month
      smb%mon_runoff(i,j,mon) = smb%mon_runoff(i,j,mon) + smb%runoff(i,j)
      ! cumulate icemelt values for each month
      smb%mon_icemelt(i,j,mon) = smb%mon_icemelt(i,j,mon) + smb%icemelt(i,j)
      ! annual mean ice temperature
      smb%t_ice(i,j) = smb%t_ice(i,j) + smb%t_prof(i,j,nl) ! bottom layer

      smb%ann_prc(i,j)       = smb%ann_prc(i,j)        + smb%prc(i,j)  * dt
      smb%ann_snow(i,j)      = smb%ann_snow(i,j)       + smb%snow(i,j) * dt
      if (smb%mask_maxice(i,j).eq.1) then
        smb%ann_smb(i,j)       = smb%ann_smb(i,j)        + (smb%snow(i,j)-smb%evp(i,j)-smb%snowmelt(i,j)-smb%icemelt(i,j)+smb%refreezing(i,j))*dt
      else
        ! set large negative value outside of allowed ice mask
        smb%ann_smb(i,j)       = smb%ann_smb(i,j)        + (smb%snow(i,j)-smb%evp(i,j)-smb%snowmelt(i,j)-smb%icemelt(i,j)+smb%refreezing(i,j))*dt &
                                                         + smb_maxice*dt/sec_year  ! kg/m2/yr * s /s*yr = kg/m2
      endif
      smb%ann_melt(i,j)      = smb%ann_melt(i,j)       + (smb%snowmelt(i,j) + smb%icemelt(i,j)) * dt
      smb%ann_icemelt(i,j)   = smb%ann_icemelt(i,j)    + smb%icemelt(i,j) * dt
      smb%ann_ablation(i,j)  = smb%ann_ablation(i,j)   + (smb%snowmelt(i,j) + smb%icemelt(i,j) + smb%evp(i,j)) * dt
      smb%ann_evp(i,j)       = smb%ann_evp(i,j)        + smb%evp(i,j) * dt
      smb%ann_runoff(i,j)    = smb%ann_runoff(i,j)     + (smb%snowmelt(i,j) + smb%icemelt(i,j) + smb%rain(i,j) - smb%refreezing(i,j)) * dt
      smb%ann_refreezing(i,j)= smb%ann_refreezing(i,j) + smb%refreezing(i,j) * dt

    enddo
    !$omp end parallel do 

    !$ time2 = omp_get_wtime()
    !$ if(l_write_timer) print *,'smb',time2-time1

    ! Monthly mean runoff and icemelt
    if (time_eom_smb) then
      smb%mon_runoff(:,:,mon)  = smb%mon_runoff(:,:,mon)/real(nstep_mon_smb,wp) ! kg/m2/s
      smb%mon_icemelt(:,:,mon) = smb%mon_icemelt(:,:,mon)/real(nstep_mon_smb,wp) ! kg/m2/s
      where (smb%mon_runoff.lt.0._wp)  smb%mon_runoff  = 0._wp
      where (smb%mon_icemelt.lt.0._wp) smb%mon_icemelt = 0._wp
    endif

    if (time_eoy_smb) then

      ! ground temperature at ice sheet model elevation, using fixed (atmospheric) lapse rate of 6.5 K/km (supported by data)
      smb%t_ground = smb%t_ground_i-6.5e-3_wp*(smb%z_sur_eff-smb%z_sur_i)
      ! bias correction
      smb%t_ground = smb%t_ground - sum(smb%t2m_bias_i(:,:,:),3)/day_year

      ! annual mean ice (10 m firn) temperature
      where (smb%mask_smb.eq.1)
        smb%t_ice = smb%t_ice/real(nstep_year_smb,wp)-T0   ! degC
      elsewhere
        smb%t_ice = min(-1e-3_wp,smb%t_ground)  ! degC
      endwhere

      ! correct if SMB anomaly approach is enabled
      if (l_smb_bias_corr) then
        smb%ann_smb = smb%ann_smb_ref-smb%ann_smb_cx_ref + smb%ann_smb
        smb%ann_prc = smb%ann_prc_ref-smb%ann_prc_cx_ref + smb%ann_prc
        smb%ann_evp = smb%ann_evp_ref-smb%ann_evp_cx_ref + smb%ann_evp
      endif
      ! compute average SMB to be written to output file to be used as reference
      ! in future simulations with l_smb_bias_corr==T
      if (year_now.ge.year_ini_smb_ref .and. year_now.le.year_end_smb_ref) then
        smb%ann_smb_avg_ref = smb%ann_smb_avg_ref + smb%ann_smb/real(year_end_smb_ref-year_ini_smb_ref+1,wp)
        smb%ann_prc_avg_ref = smb%ann_prc_avg_ref + smb%ann_prc/real(year_end_smb_ref-year_ini_smb_ref+1,wp)
        smb%ann_evp_avg_ref = smb%ann_evp_avg_ref + smb%ann_evp/real(year_end_smb_ref-year_ini_smb_ref+1,wp)
      endif
      ! write to file at the end of the reference period
      if (year_now.eq.year_end_smb_ref) then
        call smb_ref_write(smb%grid, year_ini_smb_ref, year_end_smb_ref, &
          smb%ann_smb_avg_ref, smb%ann_prc_avg_ref, smb%ann_evp_avg_ref)
      endif

      ! update mask where mass balance has to be computed with SEMI
      do i = 1,smb%grid%G%nx
        do j = 1,smb%grid%G%ny
          smb%mask_smb(i,j) = 0
          if (smb%simple%smb(i,j)*sec_year.gt.smb_crit_mask) then
            ! use information from smb PDD 
            smb%mask_smb(i,j) = 1
          else if (smb%mask_ice(i,j)==1) then
            ! include all ice points 
            smb%mask_smb(i,j) = 1
          else if (i>1 .and. i<smb%grid%G%nx .and. j>1 .and. j<smb%grid%G%ny) then
            ! additionally include all ice neighbors
            if (smb%mask_ice(i+1,j).eq.1 .or. smb%mask_ice(i-1,j).eq.1 .or. smb%mask_ice(i,j+1).eq.1 .or. smb%mask_ice(i,j-1).eq.1) then
              smb%mask_smb(i,j) = 1
            endif
          endif
        enddo
      enddo
      ! enlarge mask along borders
      do n=1,n_smb_mask_ext
        smb%mask_smb_tmp = smb%mask_smb
        do i = 2,smb%grid%G%nx-1
          do j = 2,smb%grid%G%ny-1
            if (smb%mask_smb(i,j)==0) then
              if (smb%mask_smb(i+1,j).eq.1 .or. smb%mask_smb(i-1,j).eq.1 &
                .or. smb%mask_smb(i,j+1).eq.1 .or. smb%mask_smb(i,j-1).eq.1 &
                .or. smb%mask_smb(i+1,j+1).eq.1 .or. smb%mask_smb(i+1,j-1).eq.1 &
                .or. smb%mask_smb(i-1,j+1).eq.1 .or. smb%mask_smb(i-1,j-1).eq.1 &
                ) then
                smb%mask_smb_tmp(i,j) = 1
              endif
            endif
          enddo
        enddo
        smb%mask_smb = smb%mask_smb_tmp
      enddo

    endif

  endif


  !-------------------------------------
  ! simple SMB scheme
  !-------------------------------------

  if (i_smb==3) then

    if (time_soy_smb) smb%t_ice(:,:) = 0._wp

    if (time_eom_smb) then

      !-------------------------------------
      ! map temperature from low to high resolution grid
      call map_scrip_field(smb%maps_cmn_to_ice,"t2m",smb_in%t2m, smb%t2m_i,method="mean",missing_value=-9999._dp)

      ! 2m temperature at ice sheet elevation
      smb%t2m(:,:) = smb%t2m_i(:,:) + gamma*(smb%z_sur_i(:,:)-smb%z_sur_eff(:,:))

      ! annual average temperature
      smb%t_ice(:,:) = smb%t_ice(:,:) + (smb%t2m(:,:)-T0)/real(nmon_year,wp)    ! degC

      if (time_eoy_smb) then
        !-------------------------------------
        ! surface mass balance 
        call smb_simple(smb%grid, smb_in%co2, smb_in%Smax65N, smb%mask_maxice, smb%z_sur_eff, smb%dz_dx_sur, smb%dz_dy_sur, &   ! in
          smb%ann_smb)      ! out

        ! set to zero, not available for this SMB scheme
        smb%mon_runoff(:,:,:) = 0._wp
        smb%mon_icemelt(:,:,:) = 0._wp
        smb%ann_prc(:,:)  = 0._wp
        smb%ann_runoff(:,:) = 0._wp
      endif

    endif
  endif


  return

  end subroutine smb_update
      

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s m b _ i n i t
  !   Purpose    :  initialize surface energy and mass balance interface
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_init(smb_in,smb,i_domain,grid,cmn_grid,z_bed_1min,lon_1min,lat_1min)

    implicit none

    type(smb_in_class) :: smb_in
    type(smb_class) :: smb
    integer, intent(in) :: i_domain
    type(grid_class), intent(in) :: grid 
    type(grid_class), intent(in) :: cmn_grid
    real(wp), intent(in) :: z_bed_1min(:,:)
    real(wp), intent(in) :: lon_1min(:)
    real(wp), intent(in) :: lat_1min(:)

    integer :: i, j, n, ii, jj, jj1, jj2, d
    real(wp) :: lon1, lon2, lat1, lat2, dlon_sur, dlat_sur
    integer :: nlon, nlat, ncells
    integer :: n_lon_1min, n_lat_1min, n_lon_1min_ext
    real(wp) :: lon_min, lon_max, dlon, lat_min, lat_max, dlat
    real(wp) :: z_bed_cell_mean
    real(wp) :: lon_arr_min, lon_arr_max
    real(wp) :: lon_arr(4), lat_arr(4)
    real(wp), allocatable, dimension(:) :: lon_1min_ext
    real(wp), allocatable, dimension(:,:) :: z_bed_1min_ext
    real(wp), allocatable, dimension(:) :: z_bed_cell
    integer, parameter :: ncells_max = 10000
    integer :: ni, nj
    integer :: ppos, spos
    type(grid_class) :: mask_maxice_grid
    type(map_scrip_class) :: maps_maxice_to_ice
    integer, dimension(:,:), allocatable :: maxi
    real(wp), dimension(:), allocatable :: lon_maxi, lat_maxi


    call smb_params_init
    call smb_grid_init

    ! store grid information
    smb%grid = grid

    ! define smb_in grid = common grid
    smb_in%grid = cmn_grid

    ! derive correspondence between indexes on smb and coupler grids
    allocate(smb%grid_smb_to_cmn%i_lowres(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%grid_smb_to_cmn%j_lowres(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%grid_smb_to_cmn%ncells(cmn_grid%G%nx,cmn_grid%G%ny))
    allocate(smb%grid_smb_to_cmn%ncells_ice(cmn_grid%G%nx,cmn_grid%G%ny))
    dlon_sur = cmn_grid%lon(2,1)-cmn_grid%lon(1,1)
    dlat_sur = cmn_grid%lat(1,2)-cmn_grid%lat(1,1)
    smb%grid_smb_to_cmn%ncells = 0
    do i=1,cmn_grid%G%nx
      do j=1,cmn_grid%G%ny
        lon1 = cmn_grid%lon(i,j)-0.5_wp*dlon_sur
        if (i.eq.1) lon1 = lon1-1.e-3_wp
        lon2 = cmn_grid%lon(i,j)+0.5_wp*dlon_sur
        lat1 = cmn_grid%lat(i,j)-0.5_wp*dlat_sur
        if (j.eq.1) lat1 = lat1-1.e-3_wp
        lat2 = cmn_grid%lat(i,j)+0.5_wp*dlat_sur
        do ii=1,smb%grid%G%nx
          do jj=1,smb%grid%G%ny 
            if (smb%grid%lon(ii,jj).gt.lon1 .and. smb%grid%lon(ii,jj).le.lon2 .and. smb%grid%lat(ii,jj).gt.lat1 .and. smb%grid%lat(ii,jj).le.lat2) then
              smb%grid_smb_to_cmn%i_lowres(ii,jj) = i
              smb%grid_smb_to_cmn%j_lowres(ii,jj) = j
              smb%grid_smb_to_cmn%ncells(i,j) = smb%grid_smb_to_cmn%ncells(i,j) + 1
            endif
          enddo
        enddo
      enddo
    enddo

    ! Generate mapping from cmn to smb/ice
    call map_scrip_init(smb%maps_cmn_to_ice,smb_in%grid,smb%grid,method=map_method,fldr="maps",load=.TRUE.,clean=.FALSE.)

    ! compute standard deviation of sub-grid topography from 1 min resolution bedrock topography
    allocate(z_bed_cell(ncells_max))
    allocate(smb%z_bed_std(smb%grid%G%nx,smb%grid%G%ny))
    lat_min = minval(smb%grid%lat(:,:))
    lat_max = maxval(smb%grid%lat(:,:))
    jj1 = minloc(abs(lat_1min-lat_min),1)
    jj2 = minloc(abs(lat_1min-lat_max),1)
    n_lon_1min = size(lon_1min)
    n_lat_1min = size(lat_1min)
    allocate(lon_1min_ext(3*n_lon_1min))
    lon_1min_ext(1:n_lon_1min)                = lon_1min - 360._wp
    lon_1min_ext(n_lon_1min+1:n_lon_1min*2)   = lon_1min
    lon_1min_ext(n_lon_1min*2+1:n_lon_1min*3) = lon_1min + 360._wp
    n_lon_1min_ext = size(lon_1min_ext)
    allocate(z_bed_1min_ext(3*n_lon_1min,n_lat_1min))
    z_bed_1min_ext(1:n_lon_1min,:)                = z_bed_1min
    z_bed_1min_ext(n_lon_1min+1:n_lon_1min*2,:)   = z_bed_1min
    z_bed_1min_ext(n_lon_1min*2+1:n_lon_1min*3,:) = z_bed_1min
    !$omp parallel do collapse(2) private(i,j,ii,jj,lon_arr,lon_arr_min,lon_arr_max,lat_arr,lon1,lon2,lat1,lat2,ncells,z_bed_cell,z_bed_cell_mean)
    do j=1,smb%grid%G%ny-1
      do i=1,smb%grid%G%nx-1
        lon_arr = (/ smb%grid%lon(i,j),smb%grid%lon(i+1,j),smb%grid%lon(i,j+1),smb%grid%lon(i+1,j+1) /)
        lon_arr_min = minval(lon_arr)
        lon_arr_max = maxval(lon_arr)
        if (lon_arr_max-lon_arr_min .gt. 180._wp) then
          lon1 = minval(pack(lon_arr,lon_arr>0._wp)) - 360._wp
          lon2 = maxval(pack(lon_arr,lon_arr<0._wp))
        else
          lon1 = lon_arr_min
          lon2 = lon_arr_max
        endif
        lat_arr = (/ smb%grid%lat(i,j),smb%grid%lat(i+1,j),smb%grid%lat(i,j+1),smb%grid%lat(i+1,j+1) /)
        lat1 = minval(lat_arr)
        lat2 = maxval(lat_arr)
        ncells = 0
        outer: do jj=jj1,jj2
          if (lat_1min(jj).gt.lat1 .and. lat_1min(jj).le.lat2) then
            do ii=1,n_lon_1min_ext
              if (lon_1min_ext(ii).gt.lon1 .and. lon_1min_ext(ii).le.lon2) then
                ncells = ncells + 1
                z_bed_cell(ncells) = z_bed_1min_ext(ii,jj)
                if (ncells.eq.ncells_max) exit outer
              endif
            enddo
          endif
        enddo outer
        !if (ncells.eq.ncells_max) then
        !  print *
        !  print *,'ncells too large',ncells
        !  print *,'i,j',i,j
        !  print *,'lon12,lat12',lon1,lon2,lat1,lat2
        !  print *,'lon_arr',lon_arr
        !endif
        if (ncells.gt.0) then
          z_bed_cell_mean = sum(z_bed_cell(1:ncells))/ncells
          smb%z_bed_std(i,j) = sqrt(sum((z_bed_cell(1:ncells)-z_bed_cell_mean)**2) / real(ncells-1,wp))
        else
          smb%z_bed_std(i,j) = 0._wp
        endif
      enddo
    enddo
    !$omp end parallel do
    smb%z_bed_std(smb%grid%G%nx,:) = smb%z_bed_std(smb%grid%G%nx-1,:)
    smb%z_bed_std(:,smb%grid%G%ny) = smb%z_bed_std(:,smb%grid%G%ny-1)
    deallocate(z_bed_cell)
    deallocate(lon_1min_ext)
    deallocate(z_bed_1min_ext)

    ! generate regular lat-lon grid covering smb domain 
    if (i_smb.eq.3 .or. prc_par%l_slope_effect) then

      lon_min = floor(minval(smb%grid%lon))  ! minimum lon
      lon_max = ceiling(maxval(smb%grid%lon))  ! maximum lon
      lat_min = floor(minval(smb%grid%lat))  ! minimum lat
      lat_max = ceiling(maxval(smb%grid%lat))  ! maximum lat
      dlon = (smb%grid%G%dx*1.e3_wp * 360._wp) / (2._wp*pi*r_earth*cos(0.5_wp*(lat_max+lat_min)*pi/180._wp)) ! zonal resolution in deg
      nlon = ceiling((lon_max-lon_min)/dlon)
      dlat = (smb%grid%G%dy*1.e3_wp * 360._wp) / (2._wp*pi*r_earth*cos(0.5_wp*(lat_max+lat_min)*pi/180._wp)) ! meridional resolution in deg
      nlat = ceiling((lat_max-lat_min)/dlat)
      call grid_init(smb%grid_latlon,name=trim(smb%grid%name)//trim("_latlon"),mtype="latlon",units="degrees", &
        x0=real(lon_min,dp),dx=real(dlon,dp),nx=nlon,y0=real(lat_min,dp),dy=real(dlat,dp),ny=nlat)
      call map_scrip_init(smb%maps_to_latlon,smb%grid,smb%grid_latlon,fldr="maps",load=.TRUE.,clean=.FALSE.)
      call map_scrip_init(smb%maps_from_latlon,smb%grid_latlon,smb%grid,fldr="maps",load=.TRUE.,clean=.FALSE.)

    endif


    call grid_allocate(smb%grid, smb%mask_maxice      )

    if (l_maxice) then
      ! read maximum ice extent mask

      ni = nc_size(trim(mask_maxice_file),"lon")
      nj = nc_size(trim(mask_maxice_file),"lat")
      allocate( maxi(ni,nj) )
      allocate( lon_maxi(ni) )
      allocate( lat_maxi(nj) )
      call nc_read(trim(mask_maxice_file),"lat",lat_maxi)
      call nc_read(trim(mask_maxice_file),"lon",lon_maxi)
      call nc_read(trim(mask_maxice_file),"mask",maxi)

      ! grid definition
      spos = scan(trim(mask_maxice_file),"/", BACK= .true.)+1
      ppos = scan(trim(mask_maxice_file),".", BACK= .true.)-1
      call grid_init(mask_maxice_grid,name=trim(mask_maxice_file(spos:ppos)),mtype="latlon",units="degrees",x=lon_maxi,y=lat_maxi)
      ! map to ice grid
      call map_scrip_init(maps_maxice_to_ice,mask_maxice_grid,smb%grid,method="nn",fldr="maps",load=.TRUE.,clean=.FALSE.)
      call map_scrip_field(maps_maxice_to_ice,"mask",maxi,smb%mask_maxice,method="mean",missing_value=-9999._wp, &
        filt_method="none",filt_par=[5.0*smb%grid%G%dx,smb%grid%G%dx])

      deallocate(maxi, lon_maxi, lat_maxi)

    else

      smb%mask_maxice(:,:) = 1

    endif


    allocate(smb%idx_cell_active(smb%grid%G%nx*smb%grid%G%ny))
    allocate(smb%ij_1d(2,smb%grid%G%nx*smb%grid%G%ny))
    allocate(smb%id_map(smb%grid%G%nx,smb%grid%G%ny))
    ! mapping between 1d <-> 2d
    n = 0  ! counter
    do i = 1,smb%grid%G%nx
      do j = 1,smb%grid%G%ny
        n = n + 1
        smb%id_map(i, j) = n
        smb%ij_1d(1, n) = i
        smb%ij_1d(2, n) = j
      enddo
    enddo

    if (l_regional_climate_forcing) then
      ! initialize regional climate forcing
      call fake_atm_hires_init(smb%grid,real(year_ini,wp),atm)
    endif

    ! Allocate arrays 
    if (i_domain==1) then
    call grid_allocate(smb_in%grid, smb_in%z_sur          )
    call grid_allocate(smb_in%grid, smb_in%t2m            )
    call grid_allocate(smb_in%grid, smb_in%tam            )
    call grid_allocate(smb_in%grid, smb_in%ram            )
    call grid_allocate(smb_in%grid, smb_in%gam            )
    call grid_allocate(smb_in%grid, smb_in%tstd           )
    call grid_allocate(smb_in%grid, smb_in%prc            )
    call grid_allocate(smb_in%grid, smb_in%u700           )
    call grid_allocate(smb_in%grid, smb_in%v700           )
    call grid_allocate(smb_in%grid, smb_in%wind           )
    call grid_allocate(smb_in%grid, smb_in%cld            )
    call grid_allocate(smb_in%grid, smb_in%dust           )
    call grid_allocate(smb_in%grid, smb_in%swd_toa     )
    call grid_allocate(smb_in%grid, smb_in%swd_toa_min )
    call grid_allocate(smb_in%grid, smb_in%swd_sur_vis_dir )
    call grid_allocate(smb_in%grid, smb_in%swd_sur_nir_dir )
    call grid_allocate(smb_in%grid, smb_in%swd_sur_vis_dif )
    call grid_allocate(smb_in%grid, smb_in%swd_sur_nir_dif )
    call grid_allocate(smb_in%grid, smb_in%dswd_dalb_vis_dir)  
    call grid_allocate(smb_in%grid, smb_in%dswd_dalb_nir_dir)
    call grid_allocate(smb_in%grid, smb_in%dswd_dalb_vis_dif)
    call grid_allocate(smb_in%grid, smb_in%dswd_dalb_nir_dif)
    call grid_allocate(smb_in%grid, smb_in%dswd_dz_nir_dir)
    call grid_allocate(smb_in%grid, smb_in%dswd_dz_nir_dif)
    call grid_allocate(smb_in%grid, smb_in%alb_vis_dir )
    call grid_allocate(smb_in%grid, smb_in%alb_nir_dir )
    call grid_allocate(smb_in%grid, smb_in%alb_vis_dif )
    call grid_allocate(smb_in%grid, smb_in%alb_nir_dif )
    call grid_allocate(smb_in%grid, smb_in%coszm          )
    call grid_allocate(smb_in%grid, smb_in%t_ground       )
    call grid_allocate(smb_in%grid, smb_in%lwdown         )
    call grid_allocate(smb_in%grid, smb_in%gam_lw         )
    allocate(smb_in%t2m_bias(smb_in%grid%G%nx,smb_in%grid%G%ny,nday_year))
    allocate(smb_in%prc_bias(smb_in%grid%G%nx,smb_in%grid%G%ny,nday_year))
    endif

    call grid_allocate(smb%grid, smb%mask_smb         )
    call grid_allocate(smb%grid, smb%mask_smb_tmp     )
    call grid_allocate(smb%grid, smb%mask_ice         )
    call grid_allocate(smb%grid, smb%mask_ice_old     )
    call grid_allocate(smb%grid, smb%mask_margin      )
    call grid_allocate(smb%grid, smb%z_sur            )
    call grid_allocate(smb%grid, smb%z_sur_eff        )
    call grid_allocate(smb%grid, smb%z_sur_fil        )
    call grid_allocate(smb%grid, smb%h_ice            )
    call grid_allocate(smb%grid, smb%dz_dx_sur        )
    call grid_allocate(smb%grid, smb%dz_dy_sur        )
    call grid_allocate(smb%grid, smb%dz_sur           )
    call grid_allocate(smb%grid, smb%coszm_i          )
    call grid_allocate(smb%grid, smb%z_sur_i          )
    call grid_allocate(smb%grid, smb%t2m_i            )
    call grid_allocate(smb%grid, smb%ram_i            )
    call grid_allocate(smb%grid, smb%tam_i            )
    call grid_allocate(smb%grid, smb%gam_i            )
    call grid_allocate(smb%grid, smb%tstd_i           )
    call grid_allocate(smb%grid, smb%prc_i            )
    call grid_allocate(smb%grid, smb%u700_i           )
    call grid_allocate(smb%grid, smb%v700_i           )
    call grid_allocate(smb%grid, smb%wind_i           )
    call grid_allocate(smb%grid, smb%cld_i            )
    call grid_allocate(smb%grid, smb%dust_i           )
    call grid_allocate(smb%grid, smb%swd_toa_i     )
    call grid_allocate(smb%grid, smb%swd_toa_min_i )
    call grid_allocate(smb%grid, smb%swd_sur_vis_dir_i )
    call grid_allocate(smb%grid, smb%swd_sur_nir_dir_i )
    call grid_allocate(smb%grid, smb%swd_sur_vis_dif_i )
    call grid_allocate(smb%grid, smb%swd_sur_nir_dif_i )
    call grid_allocate(smb%grid, smb%dswd_dalb_vis_dir_i)  
    call grid_allocate(smb%grid, smb%dswd_dalb_nir_dir_i)
    call grid_allocate(smb%grid, smb%dswd_dalb_vis_dif_i)
    call grid_allocate(smb%grid, smb%dswd_dalb_nir_dif_i)
    call grid_allocate(smb%grid, smb%dswd_dz_nir_dir_i)
    call grid_allocate(smb%grid, smb%dswd_dz_nir_dif_i)
    call grid_allocate(smb%grid, smb%alb_vis_dir_i )
    call grid_allocate(smb%grid, smb%alb_nir_dir_i )
    call grid_allocate(smb%grid, smb%alb_vis_dif_i )
    call grid_allocate(smb%grid, smb%alb_nir_dif_i )
    call grid_allocate(smb%grid, smb%lwdown_i         )
    call grid_allocate(smb%grid, smb%gam_lw_i         )
    call grid_allocate(smb%grid, smb%t_ground_i       )
    call grid_allocate(smb%grid, smb%tam              )
    call grid_allocate(smb%grid, smb%t2m              )
    call grid_allocate(smb%grid, smb%t_skin           )
    call grid_allocate(smb%grid, smb%t_skin_amp       )
    call grid_allocate(smb%grid, smb%t_skin_old       )
    call grid_allocate(smb%grid, smb%t_ice            )
    call grid_allocate(smb%grid, smb%q2m             )
    call grid_allocate(smb%grid, smb%pressure         )
    call grid_allocate(smb%grid, smb%mask_snow        )
    call grid_allocate(smb%grid, smb%f_snow           )
    call grid_allocate(smb%grid, smb%h_snow           )
    call grid_allocate(smb%grid, smb%w_snow           )
    call grid_allocate(smb%grid, smb%w_snow_old       )
    call grid_allocate(smb%grid, smb%w_snow_max       )
    call grid_allocate(smb%grid, smb%snowmelt         )
    call grid_allocate(smb%grid, smb%icemelt          )
    call grid_allocate(smb%grid, smb%refreezing       )
    call grid_allocate(smb%grid, smb%refreezing_sum   )
    call grid_allocate(smb%grid, smb%f_rfz_to_snow)
    call grid_allocate(smb%grid, smb%runoff           )
    call grid_allocate(smb%grid, smb%cod              )
    call grid_allocate(smb%grid, smb%albedo           )
    call grid_allocate(smb%grid, smb%f_ice            )
    call grid_allocate(smb%grid, smb%f_ice_old        )
    call grid_allocate(smb%grid, smb%alb_bg           )
    call grid_allocate(smb%grid, smb%dt_snowfree      )
    call grid_allocate(smb%grid, smb%alb_ice          )
    call grid_allocate(smb%grid, smb%alb_vis_dir )
    call grid_allocate(smb%grid, smb%alb_nir_dir )
    call grid_allocate(smb%grid, smb%alb_vis_dif )
    call grid_allocate(smb%grid, smb%alb_nir_dif )
    call grid_allocate(smb%grid, smb%alb_snow_vis_dir )
    call grid_allocate(smb%grid, smb%alb_snow_nir_dir )
    call grid_allocate(smb%grid, smb%alb_snow_vis_dif )
    call grid_allocate(smb%grid, smb%alb_snow_nir_dif )
    call grid_allocate(smb%grid, smb%snow_grain       )
    call grid_allocate(smb%grid, smb%dust_con         )
    call grid_allocate(smb%grid, smb%cld              )
    call grid_allocate(smb%grid, smb%swnet            )
    call grid_allocate(smb%grid, smb%swnet_min        )
    call grid_allocate(smb%grid, smb%swdown           )
    call grid_allocate(smb%grid, smb%lwdown           )
    call grid_allocate(smb%grid, smb%r_a              )
    call grid_allocate(smb%grid, smb%flx_g            )
    call grid_allocate(smb%grid, smb%dflxg_dT         )
    call grid_allocate(smb%grid, smb%flx_melt         )
    call grid_allocate(smb%grid, smb%flx_sh           )
    call grid_allocate(smb%grid, smb%flx_lwu          )
    call grid_allocate(smb%grid, smb%flx_lh           )
    call grid_allocate(smb%grid, smb%evp              )
    call grid_allocate(smb%grid, smb%prc              )
    call grid_allocate(smb%grid, smb%f_ele            )
    call grid_allocate(smb%grid, smb%f_wind           )
    call grid_allocate(smb%grid, smb%rain             )
    call grid_allocate(smb%grid, smb%snow             )
    call grid_allocate(smb%grid, smb%u700             )
    call grid_allocate(smb%grid, smb%v700             )
    call grid_allocate(smb%grid, smb%wind             )
    call grid_allocate(smb%grid, smb%t_ground         )
    call grid_allocate(smb%grid, smb%ann_smb          ) 
    call grid_allocate(smb%grid, smb%ann_smb_ref      ) 
    call grid_allocate(smb%grid, smb%ann_smb_cx_ref   ) 
    call grid_allocate(smb%grid, smb%ann_smb_avg_ref  ) 
    call grid_allocate(smb%grid, smb%ann_prc          ) 
    call grid_allocate(smb%grid, smb%ann_prc_ref      ) 
    call grid_allocate(smb%grid, smb%ann_prc_cx_ref   ) 
    call grid_allocate(smb%grid, smb%ann_prc_avg_ref  ) 
    call grid_allocate(smb%grid, smb%ann_snow         ) 
    call grid_allocate(smb%grid, smb%ann_ablation     ) 
    call grid_allocate(smb%grid, smb%ann_melt         ) 
    call grid_allocate(smb%grid, smb%ann_icemelt      ) 
    call grid_allocate(smb%grid, smb%ann_evp          ) 
    call grid_allocate(smb%grid, smb%ann_evp_ref      ) 
    call grid_allocate(smb%grid, smb%ann_evp_cx_ref   ) 
    call grid_allocate(smb%grid, smb%ann_evp_avg_ref  ) 
    call grid_allocate(smb%grid, smb%ann_runoff       ) 
    call grid_allocate(smb%grid, smb%ann_refreezing   ) 
    call grid_allocate(smb%grid, smb%num_lh           )
    call grid_allocate(smb%grid, smb%num_sh           )
    call grid_allocate(smb%grid, smb%num_sw           )
    call grid_allocate(smb%grid, smb%num_lw           )
    call grid_allocate(smb%grid, smb%denom_lh         )
    call grid_allocate(smb%grid, smb%denom_sh         )
    call grid_allocate(smb%grid, smb%denom_lw         )
    call grid_allocate(smb%grid, smb%f_sh             )
    call grid_allocate(smb%grid, smb%f_e              )
    call grid_allocate(smb%grid, smb%f_lh             )
    call grid_allocate(smb%grid, smb%f_lw             )
    call grid_allocate(smb%grid, smb%qsat             )
    call grid_allocate(smb%grid, smb%dqsatdT          )

    allocate(smb%t2m_bias_i(smb%grid%G%nx,smb%grid%G%ny,nday_year))
    allocate(smb%prc_bias_i(smb%grid%G%nx,smb%grid%G%ny,nday_year))
    allocate(smb%t_prof(smb%grid%G%nx,smb%grid%G%ny,0:nl))
    allocate(smb%t_prof_old(smb%grid%G%nx,smb%grid%G%ny,0:nl))
    allocate(smb%mon_runoff(smb%grid%G%nx,smb%grid%G%ny,nmon_year))
    allocate(smb%mon_icemelt(smb%grid%G%nx,smb%grid%G%ny,nmon_year))

    call grid_allocate(smb%grid, smb%simple%pdd             )
    call grid_allocate(smb%grid, smb%simple%t2m_cum         )
    call grid_allocate(smb%grid, smb%simple%snow_cum        )
    call grid_allocate(smb%grid, smb%simple%rain_cum        )
    call grid_allocate(smb%grid, smb%simple%melt            )
    call grid_allocate(smb%grid, smb%simple%melt_star       )
    call grid_allocate(smb%grid, smb%simple%runoff          )
    call grid_allocate(smb%grid, smb%simple%smb             )
    call grid_allocate(smb%grid, smb%simple%t_ice           )

    ! initialize smb_in variables
    smb_in%z_sur = 0._wp
    smb_in%t2m = 0._wp
    smb_in%tam = 0._wp
    smb_in%ram = 0._wp
    smb_in%gam = 0._wp
    smb_in%tstd = 0._wp
    smb_in%prc = 0._wp
    smb_in%u700 = 0._wp
    smb_in%v700 = 0._wp
    smb_in%wind = 0._wp
    smb_in%cld = 0._wp
    smb_in%dust = 0._wp
    smb_in%t_ground = 0._wp
    smb_in%swd_toa = 0._wp
    smb_in%swd_toa_min = 0._wp
    smb_in%swd_sur_vis_dir = 0._wp
    smb_in%swd_sur_nir_dir = 0._wp
    smb_in%swd_sur_vis_dif = 0._wp
    smb_in%swd_sur_nir_dif = 0._wp
    smb_in%dswd_dalb_vis_dir = 0._wp
    smb_in%dswd_dalb_nir_dir = 0._wp
    smb_in%dswd_dalb_vis_dif = 0._wp
    smb_in%dswd_dalb_nir_dif = 0._wp
    smb_in%dswd_dz_nir_dir = 0._wp
    smb_in%dswd_dz_nir_dif = 0._wp
    smb_in%alb_vis_dir = 0._wp
    smb_in%alb_nir_dir = 0._wp
    smb_in%alb_vis_dif = 0._wp
    smb_in%alb_nir_dif = 0._wp
    smb_in%coszm = 0._wp

    ! initialize interpolated variables
    smb%z_sur            = 0._wp 
    smb%z_sur_eff        = 0._wp 
    smb%z_sur_fil        = 0._wp 
    smb%h_ice            = 0._wp 
    smb%dz_dx_sur        = 0._wp 
    smb%dz_dy_sur        = 0._wp 
    smb%dz_sur           = 0._wp 
    smb%coszm_i          = 0._wp 
    smb%z_sur_i          = 0._wp 
    smb%t2m_i            = 0._wp 
    smb%tam_i            = 0._wp 
    smb%gam_i            = 0._wp 
    smb%ram_i            = 0._wp 
    smb%prc_i            = 0._wp 
    smb%u700_i           = 0._wp 
    smb%v700_i           = 0._wp 
    smb%wind_i           = 0._wp 
    smb%cld_i            = 0._wp 
    smb%dust_i           = 0._wp 
    smb%swd_toa_i     = 0._wp 
    smb%swd_toa_min_i = 0._wp 
    smb%swd_sur_vis_dir_i = 0._wp 
    smb%swd_sur_nir_dir_i = 0._wp 
    smb%swd_sur_vis_dif_i = 0._wp 
    smb%swd_sur_nir_dif_i = 0._wp 
    smb%dswd_dalb_vis_dir_i = 0._wp
    smb%dswd_dalb_nir_dir_i = 0._wp
    smb%dswd_dalb_vis_dif_i = 0._wp
    smb%dswd_dalb_nir_dif_i = 0._wp
    smb%dswd_dz_nir_dir_i = 0._wp
    smb%dswd_dz_nir_dif_i = 0._wp
    smb%alb_vis_dir_i = 0._wp 
    smb%alb_nir_dir_i = 0._wp 
    smb%alb_vis_dif_i = 0._wp 
    smb%alb_nir_dif_i = 0._wp 
    smb%t_ground_i    = 0._wp

    ! downscaled variables
    smb%t2m              = T0+5._wp 
    smb%t_ground         = 0._wp
    smb%q2m              = 0._wp 
    smb%pressure         = 0._wp 
    smb%prc              = 0._wp 
    smb%f_ele            = 0._wp 
    smb%f_wind           = 0._wp 
    smb%rain             = 0._wp 
    smb%snow             = 0._wp 
    smb%u700             = 0._wp 
    smb%v700             = 0._wp 
    smb%wind             = 0._wp 
    smb%cod              = 0._wp 
    smb%albedo           = 0._wp 
    smb%alb_snow_vis_dir = 0._wp 
    smb%alb_snow_nir_dir = 0._wp 
    smb%alb_snow_vis_dif = 0._wp 
    smb%alb_snow_nir_dif = 0._wp 
    smb%cld              = 0._wp 
    smb%swnet            = 0._wp 
    smb%swnet_min        = 0._wp 
    smb%swdown           = 0._wp 
    smb%lwdown           = 0._wp 
    smb%r_a              = 0._wp 
    smb%flx_g            = 0._wp 
    smb%dflxg_dT         = 0._wp 
    smb%flx_melt         = 0._wp 
    smb%flx_sh           = 0._wp 
    smb%flx_lwu          = 0._wp 
    smb%flx_lh           = 0._wp 

    ! old variables
    smb%w_snow_old = 0._wp
    smb%t_prof_old = T0
    smb%t_skin_old = T0

    smb%ann_smb_avg_ref = 0._wp 
    smb%ann_prc_avg_ref = 0._wp 
    smb%ann_evp_avg_ref = 0._wp 

    ! read files for SMB anomaly approach 
    if (l_smb_bias_corr) then
      call smb_bias_corr(smb%grid, &
        smb%ann_smb_ref, smb%ann_smb_cx_ref, &
        smb%ann_prc_ref, smb%ann_prc_cx_ref, &
        smb%ann_evp_ref, smb%ann_evp_cx_ref)
    endif

    ! initialization for bias correction 
    if (i_domain==1) then
      if (i_t2m_bias_corr.eq.1 .or. i_t2m_bias_corr.eq.2) then
        ! read bias correction file and initialize
        call t2m_bias_corr(i_t2m_bias_corr, bias_corr_file, smb_in%t2m_bias)
      else
        smb_in%t2m_bias= 0._wp
      endif
      if (i_prc_bias_corr.eq.1 .or. i_prc_bias_corr.eq.2) then
        ! read bias correction file and initialize
        call prc_bias_corr(i_prc_bias_corr, bias_corr_file, smb_in%prc_bias)
      else
        smb_in%prc_bias= 1._wp  ! precipitation bias is multiplicative
      endif
    endif
    ! interpolate to ice sheet grid
    do d=1,nday_year
      ! filter close to poles
      call filter_smb(smb_in%t2m_bias(:,:,d),real(smb_in%grid%lat(1,:),wp)) 
      call map_scrip_field(smb%maps_cmn_to_ice,"t2m_bias",smb_in%t2m_bias(:,:,d), smb%t2m_bias_i(:,:,d),method="mean",missing_value=-9999._dp)
      call map_scrip_field(smb%maps_cmn_to_ice,"prc_bias",smb_in%prc_bias(:,:,d), smb%prc_bias_i(:,:,d),method="mean",missing_value=-9999._dp)
    enddo
    ! apply scaling factor and add additional uniform bias correction
    smb%t2m_bias_i = t2m_bias_scale_fac*smb%t2m_bias_i - t2m_bias_corr_uniform

    ! initialize artificial temperature variability
    if (l_Tvar_ann) then
      smb%dTvar_ann = Tvar_ann_amp
    else
      smb%dTvar_ann = 0._wp
    endif
    if (l_Tvar_day) then
      smb%dTvar_day = Tvar_day_amp
    else
      smb%dTvar_day = 0._wp
    endif
    smb%dTvar = smb%dTvar_ann + smb%dTvar_day


    if (smb_restart) then

      ! read restart file 
      call smb_read_restart(trim(restart_in_dir)//"/smb_"//trim(smb%grid%name)//"_restart.nc",smb)
      print *,'read restart file ',trim(restart_in_dir)//"/smb_"//trim(smb%grid%name)//"_restart.nc"

    else

      ! initialize for new run
      smb%t2m            = T0+5._wp

      smb%t_ice            = T0 
      smb%mon_runoff       = 0._wp
      smb%mon_icemelt      = 0._wp

      smb%snowmelt         = 0._wp 
      smb%icemelt          = 0._wp 
      smb%evp              = 0._wp 
      smb%refreezing       = 0._wp 
      smb%refreezing_sum   = 0._wp 
      smb%runoff           = 0._wp 

      smb%mask_smb = 1
      smb%mask_ice = 1
      smb%mask_ice_old = 1
      smb%mask_snow = 0
      smb%f_ice = 1._wp
      smb%w_snow = 0._wp
      smb%w_snow_max = 0._wp
      smb%h_snow = 0._wp 
      smb%t_skin = T0
      smb%t_prof = T0

      smb%snow_grain       = snow_par%snow_grain_fresh
      smb%dust_con         = 0._wp 

      smb%alb_ice = surf_par%alb_firn
      smb%alb_vis_dir = 0._wp 
      smb%alb_nir_dir = 0._wp 
      smb%alb_vis_dif = 0._wp 
      smb%alb_nir_dif = 0._wp 

    endif

    ! count number of ice sheet cells in each climate model coarse grid cell
    smb%grid_smb_to_cmn%ncells_ice(:,:) = 0
    do i = 1,smb%grid%G%nx
      do j = 1,smb%grid%G%ny
        if (smb%mask_ice(i,j).eq.1) then
          ii = smb%grid_smb_to_cmn%i_lowres(i,j) 
          jj = smb%grid_smb_to_cmn%j_lowres(i,j) 
          smb%grid_smb_to_cmn%ncells_ice(ii,jj) = smb%grid_smb_to_cmn%ncells_ice(ii,jj)+1
        endif
      enddo
    enddo

    print*
    print*,'======================================================='
    print*,' Initialisation of smb ', trim(smb%grid%name),' complete'
    print*,'======================================================='
    print*

  return

  end subroutine smb_init
      
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  f i l t e r _ s m b
  ! Purpose  :  filter smb_in variables
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine filter_smb(var2d, lat)

    implicit none

    real(wp), dimension(:,:), intent(inout) :: var2d
    real(wp), dimension(:), intent(in) :: lat

    integer :: j, nj, nfil

    nj = size(var2d,2)

    do j=1,nj
      if (j.eq.1 .or. j.eq.nj) then
        nfil = -1
      else
        nfil = nint(1._wp/cos(lat(j)*pi/180._wp))-1
      endif
      call filter1d(var2d(:,j), nfil) 
    enddo

    return

  end subroutine filter_smb

     
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s m b _ e n d 
  ! Purpose  :  end smb
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_end(smb_in, smb)

    implicit none

    type(smb_in_class) :: smb_in
    type(smb_class) :: smb


    ! Deallocate all state variables to free memory
    call smb_dealloc(smb_in, smb)


    return

  end subroutine smb_end


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s m b _ d e a l l o c 
  ! Purpose  :  deallocate smb variables
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_dealloc(smb_in, smb)

    implicit none 

    type(smb_in_class) :: smb_in
    type(smb_class) :: smb


    ! Deallocate all state variables to free memory
    deallocate(smb%grid_smb_to_cmn%i_lowres)
    deallocate(smb%grid_smb_to_cmn%j_lowres)
    deallocate(smb%grid_smb_to_cmn%ncells)
    deallocate(smb%grid_smb_to_cmn%ncells_ice)
    deallocate(smb%idx_cell_active)
    deallocate(smb%ij_1d)
    deallocate(smb%id_map)

    if (allocated(smb_in%z_sur)) then
    deallocate(smb_in%z_sur          )
    deallocate(smb_in%t2m            )
    deallocate(smb_in%t2m_bias       )
    deallocate(smb_in%tam            )
    deallocate(smb_in%gam            )
    deallocate(smb_in%ram            )
    deallocate(smb_in%tstd           )
    deallocate(smb_in%prc            )
    deallocate(smb_in%prc_bias       )
    deallocate(smb_in%wind           )
    deallocate(smb_in%cld            )
    deallocate(smb_in%dust           )
    deallocate(smb_in%swd_toa     )
    deallocate(smb_in%swd_toa_min )
    deallocate(smb_in%swd_sur_vis_dir )
    deallocate(smb_in%swd_sur_nir_dir )
    deallocate(smb_in%swd_sur_vis_dif )
    deallocate(smb_in%swd_sur_nir_dif )
    deallocate(smb_in%dswd_dalb_vis_dir)  
    deallocate(smb_in%dswd_dalb_nir_dir)
    deallocate(smb_in%dswd_dalb_vis_dif)
    deallocate(smb_in%dswd_dalb_nir_dif)
    deallocate(smb_in%dswd_dz_nir_dir)
    deallocate(smb_in%dswd_dz_nir_dif)
    deallocate(smb_in%alb_vis_dir )
    deallocate(smb_in%alb_nir_dir )
    deallocate(smb_in%alb_vis_dif )
    deallocate(smb_in%alb_nir_dif )
    deallocate(smb_in%coszm          )
    deallocate(smb_in%t_ground       )
    deallocate(smb_in%lwdown         )
    deallocate(smb_in%gam_lw         )
    deallocate(smb_in%u700)
    deallocate(smb_in%v700) 
    endif

    deallocate(smb%mask_smb         )
    deallocate(smb%mask_smb_tmp     )
    deallocate(smb%mask_ice         )
    deallocate(smb%mask_ice_old     )
    deallocate(smb%mask_maxice      )
    deallocate(smb%mask_margin      )
    deallocate(smb%z_sur            )
    deallocate(smb%z_sur_eff        )
    deallocate(smb%z_sur_fil        )
    deallocate(smb%h_ice            )
    deallocate(smb%dz_dx_sur        )
    deallocate(smb%dz_dy_sur        )
    deallocate(smb%dz_sur           )
    deallocate(smb%coszm_i          )
    deallocate(smb%z_sur_i          )
    deallocate(smb%t2m_i            )
    deallocate(smb%t2m_bias_i       )
    deallocate(smb%tam_i            )
    deallocate(smb%gam_i            )
    deallocate(smb%ram_i            )
    deallocate(smb%tstd_i           )
    deallocate(smb%prc_i            )
    deallocate(smb%prc_bias_i       )
    deallocate(smb%u700_i           )
    deallocate(smb%v700_i           )
    deallocate(smb%wind_i           )
    deallocate(smb%cld_i            )
    deallocate(smb%dust_i           )
    deallocate(smb%swd_toa_i     )
    deallocate(smb%swd_toa_min_i )
    deallocate(smb%swd_sur_vis_dir_i )
    deallocate(smb%swd_sur_nir_dir_i )
    deallocate(smb%swd_sur_vis_dif_i )
    deallocate(smb%swd_sur_nir_dif_i )
    deallocate(smb%dswd_dalb_vis_dir_i)  
    deallocate(smb%dswd_dalb_nir_dir_i)
    deallocate(smb%dswd_dalb_vis_dif_i)
    deallocate(smb%dswd_dalb_nir_dif_i)
    deallocate(smb%dswd_dz_nir_dir_i)
    deallocate(smb%dswd_dz_nir_dif_i)
    deallocate(smb%alb_vis_dir_i )
    deallocate(smb%alb_nir_dir_i )
    deallocate(smb%alb_vis_dif_i )
    deallocate(smb%alb_nir_dif_i )
    deallocate(smb%lwdown_i         )
    deallocate(smb%gam_lw_i         )
    deallocate(smb%t_ground_i       )
    deallocate(smb%tam              )
    deallocate(smb%t2m              )
    deallocate(smb%t_skin           )
    deallocate(smb%t_skin_amp       )
    deallocate(smb%t_skin_old       )
    deallocate(smb%t_ice            )
    deallocate(smb%q2m             )
    deallocate(smb%pressure         )
    deallocate(smb%mask_snow        )
    deallocate(smb%f_snow           )
    deallocate(smb%h_snow           )
    deallocate(smb%w_snow           )
    deallocate(smb%w_snow_old       )
    deallocate(smb%w_snow_max       )
    deallocate(smb%snowmelt         )
    deallocate(smb%icemelt          )
    deallocate(smb%refreezing       )
    deallocate(smb%refreezing_sum   )
    deallocate(smb%f_rfz_to_snow)
    deallocate(smb%runoff           )
    deallocate(smb%cod              )
    deallocate(smb%albedo           )
    deallocate(smb%f_ice            )
    deallocate(smb%f_ice_old        )
    deallocate(smb%alb_bg           )
    deallocate(smb%alb_ice          )
    deallocate(smb%dt_snowfree      )
    deallocate(smb%alb_vis_dir )
    deallocate(smb%alb_nir_dir )
    deallocate(smb%alb_vis_dif )
    deallocate(smb%alb_nir_dif )
    deallocate(smb%alb_snow_vis_dir )
    deallocate(smb%alb_snow_nir_dir )
    deallocate(smb%alb_snow_vis_dif )
    deallocate(smb%alb_snow_nir_dif )
    deallocate(smb%snow_grain       )
    deallocate(smb%dust_con         )
    deallocate(smb%cld              )
    deallocate(smb%swnet            )
    deallocate(smb%swnet_min        )
    deallocate(smb%swdown           )
    deallocate(smb%lwdown           )
    deallocate(smb%r_a              )
    deallocate(smb%flx_g            )
    deallocate(smb%dflxg_dT         )
    deallocate(smb%flx_melt         )
    deallocate(smb%flx_sh           )
    deallocate(smb%flx_lwu          )
    deallocate(smb%flx_lh           )
    deallocate(smb%evp              )
    deallocate(smb%prc              )
    deallocate(smb%f_ele            )
    deallocate(smb%f_wind           )
    deallocate(smb%rain             )
    deallocate(smb%snow             )
    deallocate(smb%u700             )
    deallocate(smb%v700             )
    deallocate(smb%wind             )
    deallocate(smb%t_ground         )
    deallocate(smb%ann_smb          ) 
    deallocate(smb%ann_smb_ref      ) 
    deallocate(smb%ann_smb_cx_ref   ) 
    deallocate(smb%ann_smb_avg_ref  ) 
    deallocate(smb%ann_prc          ) 
    deallocate(smb%ann_prc_ref      ) 
    deallocate(smb%ann_prc_cx_ref   ) 
    deallocate(smb%ann_prc_avg_ref  ) 
    deallocate(smb%ann_snow         ) 
    deallocate(smb%ann_ablation     ) 
    deallocate(smb%ann_melt         ) 
    deallocate(smb%ann_icemelt      ) 
    deallocate(smb%ann_evp          ) 
    deallocate(smb%ann_evp_ref      ) 
    deallocate(smb%ann_evp_cx_ref   ) 
    deallocate(smb%ann_evp_avg_ref  ) 
    deallocate(smb%ann_runoff       ) 
    deallocate(smb%ann_refreezing   ) 
    deallocate(smb%num_lh           )
    deallocate(smb%num_sh           )
    deallocate(smb%num_sw           )
    deallocate(smb%num_lw           )
    deallocate(smb%denom_lh         )
    deallocate(smb%denom_sh         )
    deallocate(smb%denom_lw         )
    deallocate(smb%f_sh             )
    deallocate(smb%f_e              )
    deallocate(smb%f_lh             )
    deallocate(smb%f_lw             )
    deallocate(smb%qsat             )
    deallocate(smb%dqsatdT          )
    deallocate(smb%t_prof)
    deallocate(smb%t_prof_old)
    deallocate(smb%mon_runoff)
    deallocate(smb%mon_icemelt)

    deallocate(smb%simple%pdd             )
    deallocate(smb%simple%t2m_cum         )
    deallocate(smb%simple%snow_cum        )
    deallocate(smb%simple%rain_cum        )
    deallocate(smb%simple%melt            )
    deallocate(smb%simple%melt_star       )
    deallocate(smb%simple%runoff          )
    deallocate(smb%simple%smb             )
    deallocate(smb%simple%t_ice           )

    return

  end subroutine smb_dealloc 


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s m b _ w r i t e _ r e s t a r t
  ! Purpose  :  Write restart file
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_write_restart(fnm,smb)

    use dim_name, only: dim_x, dim_y, dim_month, dim_depth

    implicit none

    type(smb_class) :: smb

    character (len=*) :: fnm
    integer :: ncid
    integer :: nx, ny


    nx = smb%grid%G%nx
    ny = smb%grid%G%ny

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,dim_x,x=smb%grid%G%x0,dx=smb%grid%G%dx,nx=nx,axis="x",units="m",ncid=ncid)
    call nc_write_dim(fnm,dim_y,x=smb%grid%G%y0,dx=smb%grid%G%dy,nx=ny,axis="y",units="m",ncid=ncid)
    call nc_write_dim(fnm,dim_month,x=1,dx=1,nx=nmon_year,units="mon",ncid=ncid)
    call nc_write_dim(fnm,dim_depth,x=0,dx=1,nx=nl+1,axis="z",units="#",ncid=ncid)

    call nc_write(fnm,"mask_smb", smb%mask_smb, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="mask where to compute surface mass balance",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

    call nc_write(fnm,"mask_ice", smb%mask_ice, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="ice mask",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

    call nc_write(fnm,"mask_ice_old", smb%mask_ice_old, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="old ice mask",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

    call nc_write(fnm,"f_ice", smb%f_ice, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="ice cover fraction",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

    call nc_write(fnm,"mask_snow", smb%mask_snow, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="snow mask",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

    call nc_write(fnm,"t_skin", smb%t_skin, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="skin temperature",grid_mapping="polar_stereographic",units="K",ncid=ncid)    

    call nc_write(fnm,"t2m", smb%t2m, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="surface air 2m temperature",grid_mapping="polar_stereographic",units="K",ncid=ncid)    

    call nc_write(fnm,"t_prof", smb%t_prof, dims=[dim_x,dim_y,dim_depth],start=[1,1,1],count=[nx,ny,nl+1], &
      long_name="snow+ice/soil temperature profile",grid_mapping="polar_stereographic",units="K",ncid=ncid)    

    call nc_write(fnm,"t_ice", smb%t_ice, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="annual mean surface ice temperature",grid_mapping="polar_stereographic",units="degC",ncid=ncid)    

    call nc_write(fnm,"runoff", smb%mon_runoff, dims=[dim_x,dim_y,dim_month],start=[1,1,1],count=[nx,ny,nmon_year], &
      long_name="monthly mean runoff",grid_mapping="polar_stereographic",units="kg/m2/s",ncid=ncid)    

    call nc_write(fnm,"mon_icemelt", smb%mon_icemelt, dims=[dim_x,dim_y,dim_month],start=[1,1,1],count=[nx,ny,nmon_year], &
      long_name="monthly mean icemelt",grid_mapping="polar_stereographic",units="kg/m2/s",ncid=ncid)    

    call nc_write(fnm,"snowmelt", smb%snowmelt, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="snowmelt",grid_mapping="polar_stereographic",units="kg/m2/s",ncid=ncid)    

    call nc_write(fnm,"icemelt", smb%icemelt, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="icemelt",grid_mapping="polar_stereographic",units="kg/m2/s",ncid=ncid)    

    call nc_write(fnm,"evp", smb%evp, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="evaporation",grid_mapping="polar_stereographic",units="kg/m2/s",ncid=ncid)    

    call nc_write(fnm,"refreezing", smb%refreezing, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="refreezing",grid_mapping="polar_stereographic",units="kg/m2/s",ncid=ncid)    

    call nc_write(fnm,"refreezing_sum", smb%refreezing_sum, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="refreezing summed up over the year",grid_mapping="polar_stereographic",units="kg/m2",ncid=ncid)    

    call nc_write(fnm,"w_snow", smb%w_snow, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="snow water equivalent",grid_mapping="polar_stereographic",units="kg/m2",ncid=ncid)    

    call nc_write(fnm,"w_snow_max", smb%w_snow_max, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="seasonal max snow water equivalent",grid_mapping="polar_stereographic",units="kg/m2",ncid=ncid)    

    call nc_write(fnm,"h_snow", smb%h_snow, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="snow thickness",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

    call nc_write(fnm,"snow_grain", smb%snow_grain, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="snow grain size",grid_mapping="polar_stereographic",units="um",ncid=ncid)    

    call nc_write(fnm,"dust_con", smb%dust_con, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="dust concentration in snow",grid_mapping="polar_stereographic",units="kg/kg",ncid=ncid)    

    call nc_write(fnm,"alb_ice", smb%alb_ice, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="ice albedo",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

    call nc_write(fnm,"alb_vis_dir", smb%alb_vis_dir, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="clear-sky visible albedo",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

    call nc_write(fnm,"alb_nir_dir", smb%alb_nir_dir, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="clear-sky near-infrared albedo",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

    call nc_write(fnm,"alb_vis_dif", smb%alb_vis_dif, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="diffuse visible albedo",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

    call nc_write(fnm,"alb_nir_dif", smb%alb_nir_dif, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="diffuse near-infrared albedo",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

    call nc_close(ncid)


  end subroutine smb_write_restart


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s m b _ r e a d _ r e s t a r t
  ! Purpose  :  read smb restart file
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_read_restart(fnm,smb)

    implicit none

    type(smb_class), intent(inout) :: smb

    character (len=*) :: fnm


    call nc_read(fnm,"mask_smb", smb%mask_smb)

    call nc_read(fnm,"mask_ice", smb%mask_ice)

    call nc_read(fnm,"mask_ice_old", smb%mask_ice_old)

    call nc_read(fnm,"f_ice", smb%f_ice)

    call nc_read(fnm,"mask_snow", smb%mask_snow)

    call nc_read(fnm,"t_skin", smb%t_skin)

    call nc_read(fnm,"t2m", smb%t2m)

    call nc_read(fnm,"t_prof", smb%t_prof)

    call nc_read(fnm,"t_ice", smb%t_ice)

    call nc_read(fnm,"runoff", smb%mon_runoff)

    call nc_read(fnm,"mon_icemelt", smb%mon_icemelt)

    call nc_read(fnm,"snowmelt", smb%snowmelt)

    call nc_read(fnm,"icemelt", smb%icemelt)

    call nc_read(fnm,"evp", smb%evp)

    call nc_read(fnm,"refreezing", smb%refreezing)

    call nc_read(fnm,"refreezing_sum", smb%refreezing_sum)

    call nc_read(fnm,"w_snow", smb%w_snow)

    call nc_read(fnm,"w_snow_max", smb%w_snow_max)

    call nc_read(fnm,"h_snow", smb%h_snow)

    call nc_read(fnm,"snow_grain", smb%snow_grain)

    call nc_read(fnm,"dust_con", smb%dust_con)

    call nc_read(fnm,"alb_ice", smb%alb_ice)

    call nc_read(fnm,"alb_vis_dir", smb%alb_vis_dir)
    call nc_read(fnm,"alb_vis_dif", smb%alb_vis_dif)
    call nc_read(fnm,"alb_nir_dir", smb%alb_nir_dir)
    call nc_read(fnm,"alb_nir_dif", smb%alb_nir_dif)


    return

  end subroutine smb_read_restart

end module smb_model
