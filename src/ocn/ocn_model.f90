!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : o c n _ m o d e l
!
!  Purpose : main ocean model
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
module ocn_model

    use precision, only : wp, dp
    use nml
    use ncio

    use timer, only: time_soy_ocn, time_eoy_ocn, sec_year, year, nyears, nday_year, nstep_year_ocn, doy
    use control, only : ocn_restart, restart_in_dir, out_dir
    use constants, only : pi, cap_w, Lf, omega
    use climber_grid, only : lon, lat
    use ocn_grid, only : grid_class, ocn_grid_init, ocn_grid_update
    use ocn_grid, only : maxi, maxj, maxk, maxisles, c, dzz, zw, zro, dz, dza, mask_ocn, mask_c, k1, k1_pot, k_mix_brines, ocn_area, ocn_area_tot, ocn_vol
    use ocn_params, only : ocn_params_init, i_init, dbl
    use ocn_params, only : dt, rho0, init3_peak, init3_bg, i_saln0, saln0_const, i_fw, l_fw_corr, l_fw_melt_ice_sep, i_brines, i_brines_z, frac_brines
    use ocn_params, only : n_tracers_tot, n_tracers_ocn, n_tracers_bgc, idx_tracers_trans, age_tracer, dye_tracer, cons_tracer, l_cfc
    use ocn_params, only : i_age, i_dye, i_cons, i_cfc11, i_cfc12
    use ocn_params, only : l_mld, l_hosing, hosing_ini, l_flux_adj_atl, l_flux_adj_ant, l_flux_adj_pac, l_salinity_restore, l_q_geo
    use ocn_params, only : l_ocn_input_fix, i_ocn_input_fix, l_ocn_input_fix_write, ocn_input_fix_file
    use ocn_params, only : l_noise_fw, l_noise_flx
    use ocn_params, only : tau_scale, ke_tau_coeff
    use ocn_params, only : l_bering_flow
    use ocn_def, only : ocn_class
    use ocn_check, only : check_vel

    use eos_mod, only : eos
    use momentum_mod, only : momentum_init, momentum, fcormin_increase, fcormin_reset
    use transport_ocn_mod, only : transport_init, transport
    use hosing_mod, only : hosing_init, hosing_update
    use flux_adj_mod, only : flux_adj_init, flux_adj_update
    use noise_mod, only : noise_init, noise_fw_update, noise_flx_update
    use restore_salinity_mod, only : restore_salinity 
    use ocn_grid_update_state_mod, only : ocn_grid_update_state
    use cfc_flux_mod, only : cfc_flux
    use free_surface_mod, only : free_surface
    use bering_mod, only : bering

    !$  use omp_lib

    implicit none

    real(wp), allocatable :: sst_min_act(:,:)  !! minimum annual surface layer temperature for corals [degC]
    real(wp), allocatable :: sst_max_act(:,:)  !! maximum annual surface layer temperature for corals [degC]

    private
    public :: ocn_init, ocn_update, ocn_end
    public :: ocn_read_restart, ocn_write_restart


contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  o c n _ u p d a t e
  !   Purpose    :  update ocean
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocn_update(ocn)

    implicit none

    type(ocn_class) :: ocn

    integer :: i, j, k, n, k_mix, n_mix, k1_max
    logical :: flag_brines
    real(wp) :: tau
    real(wp) :: sal_before, sal_after
    real(wp) :: vsf_saln0, vsf_saloc
    logical :: error, error_eq, error_noneq
    real(wp) :: avg, tv1

    !$ logical, parameter :: print_omp = .false.
    !$ real(wp) :: time1,time2


    !$ if(print_omp) print *

    if (time_soy_ocn) then

      !------------------------------------------------------------------------
      ! update ocean area and volume accounting for actual ocean fraction
      !------------------------------------------------------------------------
      call ocn_grid_update(ocn%f_ocn,ocn%grid)  

      !------------------------------------------------------------------------
      ! adjust velocities and tracer concentrations after changes in ocean mask/volume
      !------------------------------------------------------------------------
      call ocn_grid_update_state(ocn%grid,ocn%ts,ocn%u)

      ! update flux adjustment to account for changes in ocean fraction
      if (l_flux_adj_atl .or. l_flux_adj_ant .or. l_flux_adj_pac) then
        call flux_adj_update(ocn%f_ocn, ocn%fw_flux_adj)
      endif

      if (i_saln0.eq.2) then
        ! compute average salinity to be used as reference for virtual salinity flux
        ocn%saln0 = sum(ocn%ts(:,:,:,2)*ocn%grid%ocn_vol(:,:,:))/ocn%grid%ocn_vol_tot ! psu
      endif

    endif

    !------------------------------------------------------------------------
    ! compute climatology of ocn input fields over last 30 years of simulation
    !------------------------------------------------------------------------

    if (l_ocn_input_fix_write .and. year.gt.(nyears-30)) then
      ! average over 30 years to remove possible noise
      avg = real(nday_year,wp)/real(nstep_year_ocn,wp) / 30._wp
      ocn%daily_input_save%stressxu(:,:,doy)    = ocn%daily_input_save%stressxu(:,:,doy)    + ocn%stressxu(:,:)     * avg 
      ocn%daily_input_save%stressyv(:,:,doy)    = ocn%daily_input_save%stressyv(:,:,doy)    + ocn%stressyv(:,:)     * avg 
      ocn%daily_input_save%stressxv(:,:,doy)    = ocn%daily_input_save%stressxv(:,:,doy)    + ocn%stressxv(:,:)     * avg 
      ocn%daily_input_save%stressyu(:,:,doy)    = ocn%daily_input_save%stressyu(:,:,doy)    + ocn%stressyu(:,:)     * avg 
      ocn%daily_input_save%wind(:,:,doy)        = ocn%daily_input_save%wind(:,:,doy)        + ocn%wind(:,:)         * avg
      ocn%daily_input_save%f_sic(:,:,doy)       = ocn%daily_input_save%f_sic(:,:,doy)       + ocn%f_sic(:,:)        * avg
      ocn%daily_input_save%slp(:,:,doy)         = ocn%daily_input_save%slp(:,:,doy)         + ocn%slp(:,:)          * avg
      ocn%daily_input_save%p_e_sic(:,:,doy)     = ocn%daily_input_save%p_e_sic(:,:,doy)     + ocn%p_e_sic(:,:)      * avg  
      ocn%daily_input_save%fw_brines(:,:,doy)   = ocn%daily_input_save%fw_brines(:,:,doy)   + ocn%fw_brines(:,:)    * avg
      ocn%daily_input_save%runoff(:,:,doy)      = ocn%daily_input_save%runoff(:,:,doy)      + ocn%runoff(:,:)       * avg
      ocn%daily_input_save%runoff_veg(:,:,doy)  = ocn%daily_input_save%runoff_veg(:,:,doy)  + ocn%runoff_veg(:,:)   * avg 
      ocn%daily_input_save%runoff_ice(:,:,doy)  = ocn%daily_input_save%runoff_ice(:,:,doy)  + ocn%runoff_ice(:,:)   * avg 
      ocn%daily_input_save%runoff_lake(:,:,doy) = ocn%daily_input_save%runoff_lake(:,:,doy) + ocn%runoff_lake(:,:)  * avg 
      ocn%daily_input_save%calving(:,:,doy)     = ocn%daily_input_save%calving(:,:,doy)     + ocn%calving(:,:)      * avg
      ocn%daily_input_save%bmelt_grd(:,:,doy)   = ocn%daily_input_save%bmelt_grd(:,:,doy)   + ocn%bmelt_grd(:,:)    * avg
      ocn%daily_input_save%bmelt_flt(:,:,doy)   = ocn%daily_input_save%bmelt_flt(:,:,doy)   + ocn%bmelt_flt(:,:)    * avg
      ocn%daily_input_save%bmelt(:,:,doy)       = ocn%daily_input_save%bmelt(:,:,doy)       + ocn%bmelt(:,:)        * avg
      ocn%daily_input_save%flx(:,:,doy)         = ocn%daily_input_save%flx(:,:,doy)         + ocn%flx(:,:)          * avg      
    endif

    !------------------------------------------------------------------------
    ! calculate wind stresses
    !------------------------------------------------------------------------

    !$ time1 = omp_get_wtime()
    do j=1,maxj
      do i=1,maxi
        if (mask_ocn(i,j).eq.1) then
          ocn%tau(1,i,j) = tau_scale*ocn%stressxu(i,j)  ! N/m2
          ocn%tau(2,i,j) = tau_scale*ocn%stressyv(i,j)  ! N/m2
          ocn%dtau_dz2(1,i,j) = 2._wp*ocn%stressxu(i,j)/dbl**2 !dzz  ! N/m4
          ocn%dtau_dz2(2,i,j) = 2._wp*ocn%stressyu(i,j)/dbl**2 !dzz  ! N/m4
          ocn%dtav_dz2(1,i,j) = 2._wp*ocn%stressxv(i,j)/dbl**2 !dzz  ! N/m4
          ocn%dtav_dz2(2,i,j) = 2._wp*ocn%stressyv(i,j)/dbl**2 !dzz  ! N/m4
        else
          ocn%tau(1,i,j) = 0._wp
          ocn%tau(2,i,j) = 0._wp
        endif
      enddo
    enddo
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'tau',time2-time1

    !------------------------------------------------------------------------
    ! wind calculations needed for mixed layer
    !------------------------------------------------------------------------

    !$ time1 = omp_get_wtime()
    if (l_mld) then
      do j=1,maxj
        do i=1,maxi
          if (mask_ocn(i,j).eq.1) then
            ! compute wind stress from magnitude of surface wind: tau = Cd*rho_a*V^2
            tau = 1.3e-3_wp*1.3_wp*ocn%wind(i,j)**2
            ! kinetic energy input into the ocean by wind stress, proportional to tau^(3/2), only over ice-free fraction
            ocn%ke_tau(i,j) = ocn%f_ocn2(i,j)/ocn%f_ocn(i,j)*(1._wp-ocn%f_sic(i,j))*ke_tau_coeff*tau**1.5_wp / sqrt(rho0) * dt  ! J/m2 or kg/s2
          else
            ocn%ke_tau(i,j) = 0._wp
          endif
        enddo
      enddo
    endif
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'mld',time2-time1

    !------------------------------------------------------------------------
    ! assign heat and freshwater fluxes 
    !------------------------------------------------------------------------

    ! surface freshwater flux, P-E+sea ice fluxes
    ocn%fw = ocn%p_e_sic

    ! update freshwater noise and add it to freshwater flux, if needed
    if (l_noise_fw .and. time_soy_ocn) then
      call noise_fw_update(real(year,wp),ocn%f_ocn,ocn%mask_coast,ocn%noise_fw,ocn%fw_noise)
    endif
    ocn%fw = ocn%fw + ocn%fw_noise

    ! freshwater Atlantic-Pacific flux adjustment and freshwater flux correction around Antarctica
    if (l_flux_adj_atl .or. l_flux_adj_ant .or. l_flux_adj_pac) then
      ocn%fw = ocn%fw + ocn%fw_flux_adj
    endif

    !$ time1 = omp_get_wtime()

    if (l_fw_melt_ice_sep) then
      ! add runoff (with contribution by melting of ice sheets removed), calving and basal melt freshwater fluxes
      ocn%fw = ocn%fw + (ocn%runoff-ocn%melt_ice) + ocn%calving + ocn%bmelt_grd + ocn%bmelt_flt
    else
      ! add runoff, calving and basal melt freshwater fluxes
      ocn%fw = ocn%fw + ocn%runoff + ocn%calving + ocn%bmelt_grd + ocn%bmelt_flt
    endif

    if (l_fw_corr) then
      ! compute annual net freshwater flux to be used for correction
      if (time_soy_ocn) ocn%fw_glob_tmp = 0._wp
      ocn%fw_glob_tmp = ocn%fw_glob_tmp + sum(ocn%fw*ocn_area)*dt   ! kg/m2/s * m2 * s = kg
      if (time_eoy_ocn) ocn%fw_glob = ocn%fw_glob_tmp/ocn_area_tot/sec_year    ! kg / m2 / s 
      ! correct freshwater flux to ensure net zero freshwater flux into ocean
      ocn%fw_corr = ocn%fw - ocn%fw_glob 
    else
      ocn%fw_corr = ocn%fw
    endif

    if (l_fw_melt_ice_sep) then
      ! apply freshwater flux from ice sheet melt and update salinity
      ! total salinity before
      sal_before = sum(ocn%ts(:,:,:,2)*ocn_vol)
      where (mask_ocn.eq.1)
        ! change in surface salinity due to applied freshwater flux from ice sheet melt
        ocn%ts(:,:,maxk,2) = ocn%ts(:,:,maxk,2) - ocn%melt_ice(:,:)*ocn%ts(:,:,maxk,2)/rho0/dz(maxk)*dt   ! kg/m2/s * psu * m3/kg / m * s = psu
      endwhere
      ! total salinity after
      sal_after = sum(ocn%ts(:,:,:,2)*ocn_vol)
      ! make sure net effect of ice sheet meltwater on total ocean salinity is zero
      where (mask_c.eq.1) 
        ocn%ts(:,:,:,2) = ocn%ts(:,:,:,2) - (sal_after-sal_before)/ocn%grid%ocn_vol_tot
      endwhere
      !print *,'sal_before,sal_mid,sal_after',sal_before/ocn%grid%ocn_vol_tot,sal_after/ocn%grid%ocn_vol_tot,sum(ocn%ts(:,:,:,2)*ocn_vol)/ocn%grid%ocn_vol_tot
    endif

    ! update freshwater hosing and add it to freshwater flux, if needed
    if (l_hosing .and. time_soy_ocn) then
      call hosing_update(real(year,wp),ocn%f_ocn,ocn%amoc,ocn%hosing,ocn%fw_hosing)
    endif
    ocn%fw_corr = ocn%fw_corr + ocn%fw_hosing 

    ! remove latent heat needed to melt ice reaching the ocean through calving and basal melt below ice shelfs
    ocn%flx = ocn%flx - ocn%calving*Lf - ocn%bmelt_flt*Lf    ! kg/m2/s * J/kg = W/m2 

    ! update heat flux noise and add it to heat flux, if needed
    if (l_noise_flx .and. time_soy_ocn) then
      call noise_flx_update(real(year,wp),ocn%f_ocn,ocn%mask_coast,ocn%noise_flx,ocn%flx_noise)
    endif
    ocn%flx = ocn%flx + ocn%flx_noise

    ! input tracer fluxes

    ocn%flx_sur = 0._wp  ! by default no input flux for tracers
    ocn%flx_bot = 0._wp  ! by default no input flux for tracers

    if (i_fw.eq.2) then
      vsf_saln0 = sum(ocn%fw_corr(:,:)*ocn_area(:,:))*ocn%saln0/rho0 / ocn_area_tot  ! kg/m2/s * m2 * psu * m3/kg /m2
      vsf_saloc = sum(ocn%fw_corr(:,:)*ocn%ts(:,:,maxk,2)*ocn_area(:,:))/rho0 / ocn_area_tot
      ocn%dvsf = vsf_saln0-vsf_saloc        ! m/s*psu
    else
      ocn%dvsf = 0._wp
    endif

    do j=1,maxj
      do i=1,maxi
        if (mask_ocn(i,j).eq.1) then
          ! convert fluxes

          ! heat flux
          ocn%flx_sur(i,j,1) = -ocn%flx(i,j)/cap_w/rho0  ! W/m2 -> m/s*K, positive upward

          ! virtual salinity flux

          if (i_brines.eq.0) then
            ! no brines
            flag_brines = .false.
          else if (i_brines.eq.1) then
            ! brines everywhere
            flag_brines = .true.
          else if (i_brines.eq.2) then
            ! brines only along coast
            if (ocn%mask_coast(i,j).eq.1) then
              flag_brines = .true.
            else
              flag_brines = .false.
            endif
          else if (i_brines.eq.3) then
            ! brines only in SH
            if (j.le.maxj/2) then
              flag_brines = .true.
            else
              flag_brines = .false.
            endif
          else if (i_brines.eq.4) then
            ! brines only in SH and along coast
            if (ocn%mask_coast(i,j).eq.1 .and. j.le.maxj/2) then
              flag_brines = .true.
            else
              flag_brines = .false.
            endif
          endif

          if (flag_brines) then
            ! freshwater flux excluding brines

            ! flux without brines
            if (i_fw.eq.1) then
              ! virtual salinity flux using reference saln0, could be bad (e.g. Yin 2010)
              ocn%flx_sur(i,j,2) = (ocn%fw_corr(i,j)-frac_brines*ocn%fw_brines(i,j))*ocn%saln0/rho0   ! kg/m2/s -> m/s*psu
            else if (i_fw.eq.2) then
              ! virtual salinity flux using local salinity, compensate over the whole surface ocean to conserve salinity
              ocn%flx_sur(i,j,2) = (ocn%fw_corr(i,j)*ocn%ts(i,j,maxk,2)-frac_brines*ocn%fw_brines(i,j)*ocn%saln0)/rho0 + ocn%dvsf  ! kg/m2/s -> m/s*psu 
            else if (i_fw.eq.3) then
              ! virtual salinity flux using local salinity
              ocn%flx_sur(i,j,2) = (ocn%fw_corr(i,j)*ocn%ts(i,j,maxk,2)-frac_brines*ocn%fw_brines(i,j)*ocn%saln0)/rho0  ! kg/m2/s -> m/s*psu 
            endif

            ! put brines at ocean bottom
            if (i_brines_z.eq.1) then
              k = k1(i,j)
              ocn%ts(i,j,k,2) = ocn%ts(i,j,k,2) - frac_brines*ocn%fw_brines(i,j)*ocn%saln0/rho0/dz(k)*dt  ! kg/m2/s * psu  * m3/kg / m * s -> psu 
            else if (i_brines_z.eq.2) then
              k1_max = k1(i,j)
              tv1 = 5000._wp
              do k=maxk,1,-1
                if (abs(ocn%z_ocn_max(i,j)-zw(k-1)).lt.tv1) then
                  k1_max = k
                  tv1 = abs(ocn%z_ocn_max(i,j)-zw(k-1))
                endif
              enddo
              k = max(k1(i,j),k1_max)
              ocn%ts(i,j,k,2) = ocn%ts(i,j,k,2) - frac_brines*ocn%fw_brines(i,j)*ocn%saln0/rho0/dz(k)*dt  ! kg/m2/s * psu  * m3/kg / m * s -> psu 
            else if (i_brines_z.eq.3) then
              k1_max = k1(i,j)
              tv1 = 5000._wp
              do k=maxk,1,-1
                if (abs(ocn%z_ocn_max(i,j)-zw(k-1)).lt.tv1) then
                  k1_max = k
                  tv1 = abs(ocn%z_ocn_max(i,j)-zw(k-1))
                endif
              enddo
              ! distribute freshwater flux from brine rejection over several layers and update salinity
              k1_max = max(k1(i,j),k1_max)
              do k=k1_max,maxk
                ocn%ts(i,j,k,2) = ocn%ts(i,j,k,2) - frac_brines*ocn%fw_brines(i,j)*dz(k)/sum(dz(k1_max:maxk))*ocn%saln0/rho0/dz(k)*dt  ! kg/m2/s * psu  * m3/kg / m * s -> psu 
              enddo
            else if (i_brines_z.eq.4) then
              do k=k1(i,j),maxk
                ocn%ts(i,j,k,2) = ocn%ts(i,j,k,2) - frac_brines*ocn%fw_brines(i,j)*dz(k)/sum(dz(k1(i,j):maxk))*ocn%saln0/rho0/dz(k)*dt  ! kg/m2/s * psu  * m3/kg / m * s -> psu 
              enddo
            endif

          else
            ! no special treatment of brines

            if (i_fw.eq.1) then
              ! virtual salinity flux using reference saln0, could be bad (e.g. Yin 2010)
              ocn%flx_sur(i,j,2) = ocn%fw_corr(i,j)*ocn%saln0/rho0   ! kg/m2/s -> m/s*psu
            else if (i_fw.eq.2) then
              ! virtual salinity flux using local salinity, compensate over the whole surface ocean to conserve salinity
              ocn%flx_sur(i,j,2) = ocn%fw_corr(i,j)*ocn%ts(i,j,maxk,2)/rho0 + ocn%dvsf  ! kg/m2/s -> m/s*psu 
            else if (i_fw.eq.3) then
              ! virtual salinity flux using local salinity
              ocn%flx_sur(i,j,2) = ocn%fw_corr(i,j)*ocn%ts(i,j,maxk,2)/rho0   ! kg/m2/s -> m/s*psu 
            endif

          endif

          ! geothermal heat flux
          if (l_q_geo) then
            ocn%flx_bot(i,j,1) = ocn%q_geo(i,j)/cap_w/rho0  ! W/m2 -> m/s*K
          endif

          ! age tracer
          if (age_tracer) then
            ocn%ts(i,j,maxk,i_age) = 0._wp ! set age to zero at surface 
            do k=1,maxk-1
              ocn%ts(i,j,k,i_age) = ocn%ts(i,j,k,i_age) + dt/sec_year ! increment age 
            enddo
          endif

          ! dye tracer
          if (dye_tracer) then
            if (j.gt.maxj/2.) then
              ocn%ts(i,j,maxk,i_dye) = 1._wp ! set to 1 in NH
            else
              ocn%ts(i,j,maxk,i_dye) = -1._wp ! set to -1 in SH
            endif
          endif

        endif
      enddo
    enddo
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'fluxes',time2-time1

    !------------------------------------------------------------------------
    ! CFC air-sea flux
    !------------------------------------------------------------------------

    if (l_cfc) then
      call cfc_flux(ocn%f_sic,ocn%wind,ocn%slp,ocn%ts(:,:,maxk,1),ocn%ts(:,:,maxk,2), &
                    ocn%cfc11_atm,ocn%cfc12_atm,ocn%ts(:,:,maxk,i_cfc11),ocn%ts(:,:,maxk,i_cfc12), &
                    ocn%flx_sur(:,:,i_cfc11), ocn%flx_sur(:,:,i_cfc12))
    endif

    !------------------------------------------------------------------------
    ! solve frictional-geostrophic balance equation
    !------------------------------------------------------------------------

    error = .true.

    n = 0
    do while (error)

      !$ time1 = omp_get_wtime()
      call momentum(ocn%f_ocn,ocn%tau,ocn%dtau_dz2,ocn%dtav_dz2,ocn%rho, &
        ocn%ub,ocn%ub_isl,ocn%u,ocn%psi)
      !$ time2 = omp_get_wtime()
      !$ if(print_omp) print *,'momentum',time2-time1

      ! check velocity range
      call check_vel(ocn%u, ocn%error, error_eq, error_noneq)

      if (error_eq) then
        n = n+1
        print *
        print *,'WARNING: CFL stability criterium not met around the equator'
        print *,'min Coriolis parameter increased temporarily by a factor of ',2**n,' for stability'
        print *,'velocities are re-diagnosed...'
        call fcormin_increase(2._wp)
        if (n.eq.2) exit
      endif

      error = error_eq !.or. error_noneq

    enddo

    ! reset min Coriolis parameter
    if (n.ge.1) then 
      call fcormin_reset
    endif

    !------------------------------------------------------------------------
    ! tracer transport (advection+diffusion+convection)
    !------------------------------------------------------------------------

    ! if CFL criterium not fullfilled, temporary reduce timestep
    error = ocn%error
    if (error) then
      dt = 0.5_wp*dt
      ocn%error = .false.
    endif

    !print *
    !print *,'before'
    !do i=1,n_tracers_tot
    !  print *,i,sum(ocn%ts(:,:,:,i)*ocn_vol(:,:,:))
    !enddo
    !$ time1 = omp_get_wtime()
    call transport(ocn%l_tracers_trans,ocn%l_tracer_dic,ocn%l_tracers_isodiff,ocn%grid%l_large_vol_change, &
                  ocn%u,ocn%ke_tau,ocn%flx_sur,ocn%flx_bot,ocn%f_ocn,ocn%mask_coast,ocn%z_ocn_max, &
                  ocn%ts,ocn%rho,ocn%nconv,ocn%dconv,ocn%kven,ocn%dven,ocn%conv_pe, &
                  ocn%mld,ocn%fdx,ocn%fdy,ocn%fdz,ocn%fax,ocn%fay,ocn%faz,ocn%dts_dt_adv,ocn%dts_dt_diff, ocn%error)
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'transport',time2-time1
    !print *,'after'
    !do i=1,n_tracers_tot
    !  print *,i,sum(ocn%ts(:,:,:,i)*ocn_vol(:,:,:))
    !enddo

    ! restore timestep 
    if (error) dt = 2._wp*dt

    !------------------------------------------------------------------------
    ! restore global salinity to reference value
    !------------------------------------------------------------------------

    if (l_salinity_restore) then
      call restore_salinity(ocn%ts(:,:,:,2))
    endif

    !------------------------------------------------------------------------
    ! diagnose free surface elevation, needed by sea ice 
    !------------------------------------------------------------------------
    !$ time1 = omp_get_wtime()
    call free_surface(ocn%rho, ocn%ssh)
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'free surface',time2-time1

    !------------------------------------------------------------------------
    ! parameterisation for Bering Strait throughflow 
    !------------------------------------------------------------------------
    if (l_bering_flow) then
      call bering(ocn%A_bering, ocn%f_ocn, ocn%ssh, ocn%saln0, &
                  ocn%ts(:,:,:,2), ocn%bering_tf, ocn%bering_fw)
    endif

    !------------------------------------------------------------------------
    ! annual min and max surface layer temperature (for corals)
    !------------------------------------------------------------------------

    !$ time1 = omp_get_wtime()
    if (time_soy_ocn) sst_min_act(:,:) = ocn%ts(:,:,maxk,1)
    if (time_soy_ocn) sst_max_act(:,:) = ocn%ts(:,:,maxk,1)
    where (ocn%ts(:,:,maxk,1).lt.ocn%sst_min) sst_min_act = ocn%ts(:,:,maxk,1)
    where (ocn%ts(:,:,maxk,1).gt.ocn%sst_max) sst_max_act = ocn%ts(:,:,maxk,1)
    if (time_eoy_ocn) ocn%sst_min = sst_min_act 
    if (time_eoy_ocn) ocn%sst_max = sst_max_act
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'sst',time2-time1
 

   return

  end subroutine ocn_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  o c n _ i n i t
  !   Purpose    :  initialize ocean
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocn_init(ocn,f_ocn,z_ocn,z_ocn_max,mask_coast,ocn_vol_tot_real,A_bering,l_daily_input_save)

    use ncio

    implicit none

    type(ocn_class), intent(out) :: ocn
    real(wp), dimension(:,:), intent(in) :: f_ocn   !! ocean fraction including floating ice
    real(wp), dimension(:,:), intent(in) :: z_ocn   !! 'real' ocean bathymetry [m]
    real(wp), dimension(:,:), intent(in) :: z_ocn_max   !! 'real' ocean bathymetry [m]
    integer, dimension(:,:), intent(in) :: mask_coast   !! 'real' ocean bathymetry [m]
    real(dp), intent(in) :: ocn_vol_tot_real
    real(dp), intent(in) :: A_bering
    logical, intent(in) :: l_daily_input_save

    integer :: i, j, k
    integer :: ni, nj, nk
    real(wp), dimension(:,:,:), allocatable :: tmp
    real(wp), dimension(:), allocatable :: tmp_depth
    real(wp) :: tmp_sum, tmp_cnt
    integer :: ncid
    character (len=256) :: fnm

    !------------------------------------------------------------------------
    ! initialize model parameters
    !------------------------------------------------------------------------
    call ocn_params_init

    !------------------------------------------------------------------------
    ! setup grid 
    !------------------------------------------------------------------------
    call ocn_grid_init(f_ocn,z_ocn,z_ocn_max,mask_coast,ocn_vol_tot_real,ocn%grid)

    !------------------------------------------------------------------------
    ! number and index of tracers
    !------------------------------------------------------------------------

    ! index 1 is reserved for temperature and 2 for salinity
    n_tracers_ocn = 2
    ! age tracer
    if (age_tracer) then
      n_tracers_ocn = n_tracers_ocn + 1
      i_age = n_tracers_ocn
    endif
    ! dye tracer
    if (dye_tracer) then
      n_tracers_ocn = n_tracers_ocn + 1
      i_dye = n_tracers_ocn
    endif
    ! conservative tracer
    if (cons_tracer) then
      n_tracers_ocn = n_tracers_ocn + 1
      i_cons = n_tracers_ocn
    endif
    ! CFC11 and CFC12
    if (l_cfc) then
      n_tracers_ocn = n_tracers_ocn + 1
      i_cfc11 = n_tracers_ocn
      n_tracers_ocn = n_tracers_ocn + 1
      i_cfc12 = n_tracers_ocn
    endif

    ! total number of tracers, including tracers for biogeochemistry
    n_tracers_tot = n_tracers_ocn + n_tracers_bgc

    !------------------------------------------------------------------------
    ! allocate 
    !------------------------------------------------------------------------
    call ocn_alloc(ocn, l_daily_input_save)

    ! default is to transport all tracers, can be changed by bgc
    ocn%l_tracers_trans = .true.
    allocate(idx_tracers_trans(n_tracers_tot))
    ! DIC tracer flag
    ocn%l_tracer_dic = .false.
    ! default is to apply isopycnal diffusion to all tracers, can be changed by bgc
    ocn%l_tracers_isodiff = .true.
    

    !------------------------------------------------------------------------
    ! read daily climatology of input fields, if required
    !------------------------------------------------------------------------
    if (l_ocn_input_fix .and. i_ocn_input_fix.eq.2) then
    
      fnm = trim(ocn_input_fix_file)
      call nc_open(fnm,ncid)
      call nc_read(fnm,"stressxu   ",ocn%daily_input_save%stressxu   ,ncid=ncid)
      call nc_read(fnm,"stressyv   ",ocn%daily_input_save%stressyv   ,ncid=ncid)
      call nc_read(fnm,"stressxv   ",ocn%daily_input_save%stressxv   ,ncid=ncid)
      call nc_read(fnm,"stressyu   ",ocn%daily_input_save%stressyu   ,ncid=ncid)
      call nc_read(fnm,"wind       ",ocn%daily_input_save%wind       ,ncid=ncid)
      call nc_read(fnm,"f_sic      ",ocn%daily_input_save%f_sic      ,ncid=ncid)
      call nc_read(fnm,"slp        ",ocn%daily_input_save%slp        ,ncid=ncid)
      call nc_read(fnm,"p_e_sic    ",ocn%daily_input_save%p_e_sic    ,ncid=ncid)
      call nc_read(fnm,"fw_brines  ",ocn%daily_input_save%fw_brines  ,ncid=ncid)
      call nc_read(fnm,"runoff     ",ocn%daily_input_save%runoff     ,ncid=ncid)
      call nc_read(fnm,"runoff_veg ",ocn%daily_input_save%runoff_veg ,ncid=ncid)
      call nc_read(fnm,"runoff_ice ",ocn%daily_input_save%runoff_ice ,ncid=ncid)
      call nc_read(fnm,"runoff_lake",ocn%daily_input_save%runoff_lake,ncid=ncid)
      call nc_read(fnm,"calving    ",ocn%daily_input_save%calving    ,ncid=ncid)
      call nc_read(fnm,"bmelt_grd  ",ocn%daily_input_save%bmelt_grd  ,ncid=ncid)
      call nc_read(fnm,"bmelt_flt  ",ocn%daily_input_save%bmelt_flt  ,ncid=ncid)
      call nc_read(fnm,"bmelt      ",ocn%daily_input_save%bmelt      ,ncid=ncid)
      call nc_read(fnm,"flx        ",ocn%daily_input_save%flx        ,ncid=ncid)
      call nc_close(ncid)

    endif


    !------------------------------------------------------------------------
    ! restart 
    !------------------------------------------------------------------------
    if (ocn_restart) then

      !------------------------------------------------------------------------
      ! read restart file
      !------------------------------------------------------------------------
      call ocn_read_restart(trim(restart_in_dir)//"/ocn_restart.nc",ocn)

      if (cons_tracer) then
        ocn%ts(:,:,:,i_cons) = 1._wp
      endif

      if (l_cfc) then
        ocn%ts(:,:,:,i_cfc11) = 0._wp
        ocn%ts(:,:,:,i_cfc12) = 0._wp
      endif

    else  

      ! set initial conditions   

      if (i_init.eq.1) then
        ! set latitude dependence for initial temperature and uniform initial salinity values

        !do k=1,maxk
          do j=1,maxj
            !do i=1,maxi
              ocn%ts(:,j,:,1) = 30._wp*c(j)
            !enddo
          enddo
        !enddo
        ocn%ts(:,:,:,2) = saln0_const

      else if (i_init.eq.2) then
        ! initialize using present day observations

        ocn%ts(:,:,:,1) = 0._wp
        ocn%ts(:,:,:,2) = saln0_const

        ! initialize 3D temperature and salinity from World Ocean Atlas 2013 data, annual mean
        fnm = "input/WOA13_5x5.nc"
        ni = nc_size(fnm,"lon")
        nj = nc_size(fnm,"lat")
        nk = nc_size(fnm,"depth")
        allocate(tmp_depth(nk))
        call nc_read(fnm,"depth",tmp_depth,start=[1],count=[nk])
        allocate(tmp(ni,nj,nk))
        ! temperature
        call nc_read(fnm,"t",tmp,start=[1,1,1,13],count=[ni,nj,nk,1])
        do k=1,maxk
          do j=1,maxj
            do i=1,maxi
              if (k.ge.k1(i,j)) then
                tmp_sum = sum(tmp(i,j,:), tmp(i,j,:)>-990. .and. tmp_depth<-zw(k-1) .and. tmp_depth>=-zw(k))
                tmp_cnt = count(tmp(i,j,:)>-990. .and. tmp_depth<-zw(k-1) .and. tmp_depth>=-zw(k)) 
                if (tmp_cnt.gt.0._wp) then
                  ocn%ts(i,j,k,1) = tmp_sum/tmp_cnt + 0.1_wp/1000._wp*zro(k) ! convert to potential temperature using 'lapse rate' of 0.1 K/km
                else
                  ocn%ts(i,j,k,1) = 0._wp
                endif
              else ! outside of domain
                ocn%ts(i,j,k,1) = 0._wp
              endif
            enddo
          enddo
        enddo
        ! salinity
        call nc_read(fnm,"s",tmp,start=[1,1,1,13],count=[ni,nj,nk,1])
        do k=1,maxk
          do j=1,maxj
            do i=1,maxi
              if (k.ge.k1(i,j)) then
                tmp_sum = sum(tmp(i,j,:), tmp(i,j,:)>-990. .and. tmp_depth<-zw(k-1) .and. tmp_depth>=-zw(k))
                tmp_cnt = count(tmp(i,j,:)>-990. .and. tmp_depth<-zw(k-1) .and. tmp_depth>=-zw(k))
                if (tmp_cnt.gt.0._wp) then
                  ocn%ts(i,j,k,2) = tmp_sum/tmp_cnt
                else
                  ocn%ts(i,j,k,2) = saln0_const
                endif
              else ! outside of domain
                ocn%ts(i,j,k,2) = saln0_const
              endif
            enddo
          enddo
        enddo

      else if (i_init.eq.3) then
         ! initial conditions for Eocene following Lunt et al, GMD, 2017, sec 4.2.6
         ! The ocean should be initialized as stationary,
         ! with no initial sea ice, and a zonally symmetric
         ! temperature (T, \degC) and globally constant
         ! salinity (S, psu) distribution given by the following:
         ! T[\degC] = { (\frac{5000-z}{5000} 25\cos(\phi)) + 15 if z <= 5000m
         !            { 15                                      if z >  5000m
         ! S[psu] = 34.7
         ! where \phi is latitude, and z is depth of the ocean (metres below surface)
         ! NOTE: zro(k) is negative with depth!
         do k=1,maxk
            do j=1,maxj
               !do i=1,maxi
                  if (zro(k).ge.-5000) then
                     !ocn%ts(:,j,k,1) = ((5000.0_wp+zro(k))/5000.0_wp)*25.0_wp*c(j) + 15.0_wp
                     ocn%ts(:,j,k,1) = ((5000.0_wp+zro(k))/5000.0_wp)*init3_peak*c(j) + init3_bg
                  else
                     !ocn%ts(:,j,k,1) = 15.0_wp
                     ocn%ts(:,j,k,1) = init3_bg
                  endif
               !enddo
            enddo
         enddo
         ocn%ts(:,:,:,2) = saln0_const

      else
         print *,'Error: unknow value for &ocn_par/i_init', i_init
      endif

      do k=1,maxk
         print *,'k, zro(k)', k, zro(k), ((5000.0_wp+zro(k))/5000.0_wp)
         do j=1,maxj
            !do i=1,maxi
                print *, j, c(j), ocn%ts(1,j,k,1)
            !enddo
         enddo
      enddo

      ! initialize velocity to zero
      ocn%u = 0._wp

      ! initialize psi and barotropic velocity
      ocn%psi = 0._wp
      ocn%ub  = 0._wp

      if (age_tracer) then
        ocn%ts(:,:,:,i_age) = 0._wp
      endif

      if (dye_tracer) then
        ocn%ts(:,:,:,i_dye) = 0._wp
      endif

      if (cons_tracer) then
        ocn%ts(:,:,:,i_cons) = 1._wp
      endif

      if (l_cfc) then
        ocn%ts(:,:,:,i_cfc11) = 0._wp
        ocn%ts(:,:,:,i_cfc12) = 0._wp
      endif

      ocn%fw_glob = 0._wp

    endif

    if (i_saln0.eq.1) then
      ocn%saln0 = saln0_const
    else if (i_saln0.eq.2) then
      ocn%saln0 = sum(ocn%ts(:,:,:,2)*ocn%grid%ocn_vol(:,:,:))/ocn%grid%ocn_vol_tot ! psu
    endif
    
    ocn%sst_min = ocn%ts(:,:,maxk,1)-1._wp
    ocn%sst_max = ocn%ts(:,:,maxk,1)+1._wp

    ! compute density
    do i=1,maxi
      do j=1,maxj
        do k=1,maxk
          ocn%rho(i,j,k) = eos(ocn%ts(i,j,k,1),ocn%ts(i,j,k,2),zro(k))
        enddo
      enddo
    enddo

    ! convection arrays
    do i=1,maxi
      do j=1,maxj
        ocn%nconv(i,j) = 0
        ocn%kven(i,j) = 0
        ocn%dconv(i,j) = 0._wp
        ocn%dven(i,j) = 0._wp
        ocn%conv_pe(i,j) = 0._wp
      enddo
    enddo

    ! initialize variables for frictional geostrophic balance equation
    call momentum_init

    ! initialize transport variables
    call transport_init

    ! initialize freshwater hosing
    if (l_hosing) then
      call hosing_init(f_ocn,real(year,wp))
      ocn%hosing = hosing_ini
    else
      ocn%fw_hosing = 0._wp
      ocn%hosing = 0._wp
    endif

    ! initialize freshwater flux correction
    if (l_flux_adj_atl .or. l_flux_adj_ant .or. l_flux_adj_pac) then
      call flux_adj_init(f_ocn)
    else
      ocn%fw_flux_adj = 0._wp
    endif

    ! initialize noise
    if (l_noise_fw .or. l_noise_flx) then
      call noise_init(f_ocn)
    endif
    if (.not.l_noise_fw) then
      ocn%fw_noise = 0._wp
      ocn%noise_fw = 0._wp
    endif
    if (.not.l_noise_flx) then
      ocn%flx_noise = 0._wp
      ocn%noise_flx = 0._wp
    endif

    ocn%melt_ice(:,:) = 0._wp

    ocn%amoc = 0._wp

    ocn%A_bering = A_bering
    ocn%bering_tf = 0._wp
    ocn%bering_fw = 0._wp

    ocn%error = .false.

    print*
    print*,'======================================================='
    print*,' Initialisation of OCEAN complete'
    print*,'======================================================='
    print*


  return

  end subroutine ocn_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  o c n _ a l l o c
  !   Purpose    :  allocate ocean
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocn_alloc(ocn, l_daily_input_save)

    implicit none 

    type(ocn_class) :: ocn
    logical, intent(in) :: l_daily_input_save


    allocate(ocn%l_tracers_trans(n_tracers_tot))
    allocate(ocn%l_tracer_dic(n_tracers_tot))
    allocate(ocn%l_tracers_isodiff(n_tracers_tot))

    allocate(ocn%z_ocn_max(maxi,maxj))
    allocate(ocn%f_ocn(maxi,maxj))
    allocate(ocn%f_ocn2(maxi,maxj))
    allocate(ocn%f_sic(maxi,maxj))
    allocate(ocn%mask_coast(maxi,maxj))

    allocate(ocn%wind(maxi,maxj))
    allocate(ocn%slp(maxi,maxj))
    allocate(ocn%stressxu(maxi,maxj))
    allocate(ocn%stressxv(maxi,maxj)) 
    allocate(ocn%stressyu(maxi,maxj)) 
    allocate(ocn%stressyv(maxi,maxj)) 
    allocate(ocn%flx(maxi,maxj))
    allocate(ocn%fw(maxi,maxj))
    allocate(ocn%p_e_sic(maxi,maxj))
    allocate(ocn%fw_brines(maxi,maxj))
    allocate(ocn%runoff(maxi,maxj))
    allocate(ocn%runoff_veg(maxi,maxj))
    allocate(ocn%runoff_ice(maxi,maxj))
    allocate(ocn%runoff_lake(maxi,maxj))
    allocate(ocn%melt_ice(maxi,maxj))
    allocate(ocn%calving(maxi,maxj))
    allocate(ocn%bmelt(maxi,maxj))
    allocate(ocn%bmelt_grd(maxi,maxj))
    allocate(ocn%bmelt_flt(maxi,maxj))

    allocate(ocn%u(3,0:maxi,0:maxj,maxk))
    allocate(ocn%ub(2,0:maxi+1,0:maxj))
    allocate(ocn%ub_isl(2,maxi,maxj,maxisles))
    allocate(ocn%ts(maxi,maxj,maxk,n_tracers_tot))
    allocate(ocn%fdx(0:maxi,0:maxj,0:maxk,n_tracers_tot))
    allocate(ocn%fdy(0:maxi,0:maxj,0:maxk,n_tracers_tot))
    allocate(ocn%fdz(0:maxi,0:maxj,0:maxk,n_tracers_tot))
    allocate(ocn%fax(0:maxi,0:maxj,0:maxk,n_tracers_tot))
    allocate(ocn%fay(0:maxi,0:maxj,0:maxk,n_tracers_tot))
    allocate(ocn%faz(0:maxi,0:maxj,0:maxk,n_tracers_tot))
    allocate(ocn%dts_dt_adv(maxi,maxj,maxk,n_tracers_tot))
    allocate(ocn%dts_dt_diff(maxi,maxj,maxk,n_tracers_tot))
    allocate(ocn%sst_min(maxi,maxj))
    allocate(ocn%sst_max(maxi,maxj))
    allocate(ocn%flx_sur(maxi,maxj,n_tracers_tot))
    allocate(ocn%flx_bot(maxi,maxj,n_tracers_tot))
    allocate(ocn%tau(2,maxi,maxj))
    allocate(ocn%dtau_dz2(2,maxi,maxj))
    allocate(ocn%dtav_dz2(2,maxi,maxj))
    allocate(ocn%rho(maxi,maxj,maxk))
    allocate(ocn%fw_hosing(maxi,maxj))
    allocate(ocn%fw_flux_adj(maxi,maxj))
    allocate(ocn%fw_noise(maxi,maxj))
    allocate(ocn%flx_noise(maxi,maxj))
    allocate(ocn%psi(0:maxi,0:maxj))
    allocate(ocn%ke_tau(maxi,maxj))
    allocate(ocn%mld(maxi,maxj))
    allocate(ocn%nconv(maxi,maxj))
    allocate(ocn%kven(maxi,maxj))
    allocate(ocn%dconv(maxi,maxj))
    allocate(ocn%dven(maxi,maxj))
    allocate(ocn%conv_pe(maxi,maxj))
    allocate(ocn%ssh(maxi,maxj))
    allocate(ocn%q_geo(maxi,maxj))

    allocate(sst_min_act(maxi,maxj))
    allocate(sst_max_act(maxi,maxj))

    if (l_daily_input_save .or. l_ocn_input_fix .or. l_ocn_input_fix_write) then
      allocate(ocn%daily_input_save%stressxu(maxi,maxj,nday_year))
      allocate(ocn%daily_input_save%stressyv(maxi,maxj,nday_year))
      allocate(ocn%daily_input_save%stressxv(maxi,maxj,nday_year))
      allocate(ocn%daily_input_save%stressyu(maxi,maxj,nday_year))
      allocate(ocn%daily_input_save%wind(maxi,maxj,nday_year))
      allocate(ocn%daily_input_save%f_sic(maxi,maxj,nday_year))
      allocate(ocn%daily_input_save%slp(maxi,maxj,nday_year))
      allocate(ocn%daily_input_save%p_e_sic(maxi,maxj,nday_year))
      allocate(ocn%daily_input_save%fw_brines(maxi,maxj,nday_year))
      allocate(ocn%daily_input_save%runoff(maxi,maxj,nday_year))
      allocate(ocn%daily_input_save%runoff_veg(maxi,maxj,nday_year))
      allocate(ocn%daily_input_save%runoff_ice(maxi,maxj,nday_year))
      allocate(ocn%daily_input_save%runoff_lake(maxi,maxj,nday_year))
      allocate(ocn%daily_input_save%calving(maxi,maxj,nday_year))
      allocate(ocn%daily_input_save%bmelt_grd(maxi,maxj,nday_year))
      allocate(ocn%daily_input_save%bmelt_flt(maxi,maxj,nday_year))
      allocate(ocn%daily_input_save%bmelt(maxi,maxj,nday_year))
      allocate(ocn%daily_input_save%flx(maxi,maxj,nday_year))
      ocn%daily_input_save%stressxu = 0._wp
      ocn%daily_input_save%stressyv = 0._wp
      ocn%daily_input_save%stressxv = 0._wp
      ocn%daily_input_save%stressyu = 0._wp
      ocn%daily_input_save%wind = 0._wp
      ocn%daily_input_save%f_sic = 0._wp
      ocn%daily_input_save%slp = 0._wp
      ocn%daily_input_save%p_e_sic = 0._wp
      ocn%daily_input_save%fw_brines = 0._wp
      ocn%daily_input_save%runoff = 0._wp
      ocn%daily_input_save%runoff_veg = 0._wp
      ocn%daily_input_save%runoff_ice = 0._wp
      ocn%daily_input_save%runoff_lake = 0._wp
      ocn%daily_input_save%calving = 0._wp
      ocn%daily_input_save%bmelt_grd = 0._wp
      ocn%daily_input_save%bmelt_flt = 0._wp
      ocn%daily_input_save%bmelt = 0._wp
      ocn%daily_input_save%flx = 0._wp
    endif

   return

  end subroutine ocn_alloc 


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  o c n _ e n d
  !   Purpose    :  end ocean
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocn_end(ocn)

    implicit none

    type(ocn_class) :: ocn

    integer :: ncid
    character(len=256) :: fnm
    real(wp), parameter :: missing_value = -9999._wp

    if (l_ocn_input_fix_write) then
    
      fnm = trim(out_dir)//"/ocn_input_fix.nc"
      call nc_create(fnm)
      call nc_open(fnm,ncid)
      call nc_write_dim(fnm,"lon",x=lon,axis="x",ncid=ncid)
      call nc_write_dim(fnm,"lat",x=lat,axis="y",ncid=ncid)
      call nc_write_dim(fnm,"day",x=1._wp,dx=1._wp,nx=nday_year,units="doy",ncid=ncid)
      call nc_write(fnm,"stressxu   ",ocn%daily_input_save%stressxu    ,dims=["lon","lat","day"],start=[1,1,1],count=[maxi,maxj,nday_year],long_name="",units="",missing_value=missing_value,ncid=ncid)
      call nc_write(fnm,"stressyv   ",ocn%daily_input_save%stressyv    ,dims=["lon","lat","day"],start=[1,1,1],count=[maxi,maxj,nday_year],long_name="",units="",missing_value=missing_value,ncid=ncid)
      call nc_write(fnm,"stressxv   ",ocn%daily_input_save%stressxv    ,dims=["lon","lat","day"],start=[1,1,1],count=[maxi,maxj,nday_year],long_name="",units="",missing_value=missing_value,ncid=ncid)
      call nc_write(fnm,"stressyu   ",ocn%daily_input_save%stressyu    ,dims=["lon","lat","day"],start=[1,1,1],count=[maxi,maxj,nday_year],long_name="",units="",missing_value=missing_value,ncid=ncid)
      call nc_write(fnm,"wind       ",ocn%daily_input_save%wind        ,dims=["lon","lat","day"],start=[1,1,1],count=[maxi,maxj,nday_year],long_name="",units="",missing_value=missing_value,ncid=ncid)
      call nc_write(fnm,"f_sic      ",ocn%daily_input_save%f_sic       ,dims=["lon","lat","day"],start=[1,1,1],count=[maxi,maxj,nday_year],long_name="",units="",missing_value=missing_value,ncid=ncid)
      call nc_write(fnm,"slp        ",ocn%daily_input_save%slp         ,dims=["lon","lat","day"],start=[1,1,1],count=[maxi,maxj,nday_year],long_name="",units="",missing_value=missing_value,ncid=ncid)
      call nc_write(fnm,"p_e_sic    ",ocn%daily_input_save%p_e_sic     ,dims=["lon","lat","day"],start=[1,1,1],count=[maxi,maxj,nday_year],long_name="",units="",missing_value=missing_value,ncid=ncid)
      call nc_write(fnm,"fw_brines  ",ocn%daily_input_save%fw_brines   ,dims=["lon","lat","day"],start=[1,1,1],count=[maxi,maxj,nday_year],long_name="",units="",missing_value=missing_value,ncid=ncid)
      call nc_write(fnm,"runoff     ",ocn%daily_input_save%runoff      ,dims=["lon","lat","day"],start=[1,1,1],count=[maxi,maxj,nday_year],long_name="",units="",missing_value=missing_value,ncid=ncid)
      call nc_write(fnm,"runoff_veg ",ocn%daily_input_save%runoff_veg  ,dims=["lon","lat","day"],start=[1,1,1],count=[maxi,maxj,nday_year],long_name="",units="",missing_value=missing_value,ncid=ncid)
      call nc_write(fnm,"runoff_ice ",ocn%daily_input_save%runoff_ice  ,dims=["lon","lat","day"],start=[1,1,1],count=[maxi,maxj,nday_year],long_name="",units="",missing_value=missing_value,ncid=ncid)
      call nc_write(fnm,"runoff_lake",ocn%daily_input_save%runoff_lake ,dims=["lon","lat","day"],start=[1,1,1],count=[maxi,maxj,nday_year],long_name="",units="",missing_value=missing_value,ncid=ncid)
      call nc_write(fnm,"calving    ",ocn%daily_input_save%calving     ,dims=["lon","lat","day"],start=[1,1,1],count=[maxi,maxj,nday_year],long_name="",units="",missing_value=missing_value,ncid=ncid)
      call nc_write(fnm,"bmelt_grd  ",ocn%daily_input_save%bmelt_grd   ,dims=["lon","lat","day"],start=[1,1,1],count=[maxi,maxj,nday_year],long_name="",units="",missing_value=missing_value,ncid=ncid)
      call nc_write(fnm,"bmelt_flt  ",ocn%daily_input_save%bmelt_flt   ,dims=["lon","lat","day"],start=[1,1,1],count=[maxi,maxj,nday_year],long_name="",units="",missing_value=missing_value,ncid=ncid)
      call nc_write(fnm,"bmelt      ",ocn%daily_input_save%bmelt       ,dims=["lon","lat","day"],start=[1,1,1],count=[maxi,maxj,nday_year],long_name="",units="",missing_value=missing_value,ncid=ncid)
      call nc_write(fnm,"flx        ",ocn%daily_input_save%flx         ,dims=["lon","lat","day"],start=[1,1,1],count=[maxi,maxj,nday_year],long_name="",units="",missing_value=missing_value,ncid=ncid)
      call nc_close(ncid)

    endif

    call ocn_dealloc(ocn)
 

   return

  end subroutine ocn_end


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  o c n _ d e a l l o c
  !   Purpose    :  deallocate ocean
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocn_dealloc(ocn)

    implicit none 

    type(ocn_class) :: ocn


    deallocate(ocn%z_ocn_max)
    deallocate(ocn%f_ocn)
    deallocate(ocn%f_ocn2)
    deallocate(ocn%f_sic)
    deallocate(ocn%mask_coast)

    deallocate(ocn%wind)
    deallocate(ocn%slp)
    deallocate(ocn%stressxu)
    deallocate(ocn%stressxv)
    deallocate(ocn%stressyu)
    deallocate(ocn%stressyv)
    deallocate(ocn%flx)
    deallocate(ocn%fw)
    deallocate(ocn%p_e_sic)
    deallocate(ocn%fw_brines)

    deallocate(ocn%u)
    deallocate(ocn%ub)
    deallocate(ocn%ub_isl)
    deallocate(ocn%ts)
    deallocate(ocn%fdx)
    deallocate(ocn%fdy)
    deallocate(ocn%fdz)
    deallocate(ocn%fax)
    deallocate(ocn%fay)
    deallocate(ocn%faz)
    deallocate(ocn%dts_dt_adv)
    deallocate(ocn%dts_dt_diff)
    deallocate(ocn%flx_sur)
    deallocate(ocn%flx_bot)
    deallocate(ocn%tau)
    deallocate(ocn%dtau_dz2)
    deallocate(ocn%dtav_dz2)
    deallocate(ocn%rho)
    deallocate(ocn%fw_hosing)
    deallocate(ocn%fw_flux_adj)
    deallocate(ocn%fw_noise)
    deallocate(ocn%flx_noise)
    deallocate(ocn%psi)
    deallocate(ocn%ke_tau)
    deallocate(ocn%mld)
    deallocate(ocn%nconv)
    deallocate(ocn%kven)
    deallocate(ocn%dconv)
    deallocate(ocn%dven)
    deallocate(ocn%conv_pe)
    deallocate(ocn%ssh)
    deallocate(ocn%q_geo)

    deallocate(sst_min_act)
    deallocate(sst_max_act)

   return

  end subroutine ocn_dealloc 


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  o c n _ w r i t e _ r e s t a r t
  ! Purpose  :  Write restart netcdf file 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocn_write_restart(fnm,ocn)

    implicit none

    character (len=*) :: fnm
    type(ocn_class) :: ocn

    integer :: i 

    call nc_create(fnm)
    call nc_write_dim(fnm,"lon",x=ocn%grid%lon,axis="x")
    call nc_write_dim(fnm,"lat",x=ocn%grid%lat,axis="y")
    call nc_write_dim(fnm,"lon1",x=[2.*ocn%grid%lon(1)-ocn%grid%lon(2),ocn%grid%lon],axis="x")
    call nc_write_dim(fnm,"lat1",x=[2.*ocn%grid%lat(1)-ocn%grid%lat(2),ocn%grid%lat],axis="y")
    call nc_write_dim(fnm,"zro",x=ocn%grid%zro,units="m",axis="z")
    call nc_write_dim(fnm,"zw",x=ocn%grid%zw,units="m",axis="z")
    call nc_write_dim(fnm,"dz",x=ocn%grid%dz,units="m",axis="z")
    call nc_write_dim(fnm,"dza",x=ocn%grid%dza,units="m",axis="z")
    call nc_write_dim(fnm,"dir2",x=[1,2])
    call nc_write_dim(fnm,"dir3",x=[1,2,3])
    call nc_write_dim(fnm,"ntrc",x=[(i, i=1,n_tracers_ocn)],units="#tracer")
    call nc_write_dim(fnm,"c",x=1)

    call nc_write(fnm,"f_ocn", ocn%f_ocn,   dims=["lon  ","lat  "],long_name="ocean fraction",units="/")
    call nc_write(fnm,"k1_pot",ocn%grid%k1_pot,  dims=["lon  ","lat  "],long_name="index of bottom layer",units="/")
    call nc_write(fnm,"ocn_area", ocn%grid%ocn_area,  dims=["lon  ","lat  "],long_name="ocean area",units="/")
    call nc_write(fnm,"ocn_vol",  ocn%grid%ocn_vol,   dims=["lon  ","lat  ", "zro  "],long_name="ocean volume",units="/")

    call nc_write(fnm,"u",     ocn%u,    dims=["dir3","lon1","lat1","zro "],long_name="3D velocity field",units="m/s")
    call nc_write(fnm,"ub",    ocn%ub(:,0:maxi,0:maxj),   dims=["dir2 ","lon1 ","lat1 "],long_name="barotropic velocity",units="m/s")
    call nc_write(fnm,"ts",    ocn%ts(:,:,:,1:2),  dims=["lon ","lat ","zro ","ntrc"],long_name="3D tracer fields",units="/")
    if (age_tracer) then
      call nc_write(fnm,"age",    ocn%ts(:,:,:,i_age),  dims=["lon ","lat ","zro "],long_name="age tracer fields",units="/")
    endif
    if (age_tracer) then
      call nc_write(fnm,"dye",    ocn%ts(:,:,:,i_dye),  dims=["lon ","lat ","zro "],long_name="dye tracer fields",units="/")
    endif

    call nc_write(fnm,"fw_glob", ocn%fw_glob,   dims=["c"],long_name="net ocean freshwater flux",units="kg/m2/s")

   return

  end subroutine ocn_write_restart


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  o c n _ r e a d _ r e s t a r t
  ! Purpose  :  read restart netcdf file 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocn_read_restart(fnm,ocn)

    implicit none

    character (len=*) :: fnm
    type(ocn_class) :: ocn

    call nc_read(fnm,"zro",  ocn%grid%zro)
    call nc_read(fnm,"zw",   ocn%grid%zw)
    call nc_read(fnm,"dz",   ocn%grid%dz)
    call nc_read(fnm,"dza",  ocn%grid%dza)

    !call nc_read(fnm,"f_ocn",  ocn%f_ocn)
    call nc_read(fnm,"k1_pot", ocn%grid%k1_pot)
    call nc_read(fnm,"ocn_area", ocn%grid%ocn_area)
    call nc_read(fnm,"ocn_vol",  ocn%grid%ocn_vol)

    k1_pot = ocn%grid%k1_pot
    k1(1:maxi,1:maxj) = k1_pot
    ! periodic boundary conditions
    k1(0,1:maxj) = k1(maxi,1:maxj)
    k1(maxi+1,1:maxj) = k1(1,1:maxj)
    ! North Pole and South Pole 'islands'
    k1(:,0) = 99
    k1(:,maxj+1) = 99

    call nc_read(fnm,"u",     ocn%u(:,0:maxi,0:maxj,1:maxk))
    call nc_read(fnm,"ub",    ocn%ub(:,0:maxi,0:maxj))
    ocn%ub(:,maxi+1,0:maxj) = ocn%ub(:,1,0:maxj)
    call nc_read(fnm,"ts",    ocn%ts(:,:,:,1:2))    ! temperature and salinity

    if (age_tracer) then
      if (nc_exists_var(fnm,"age")) then
        call nc_read(fnm,"age",    ocn%ts(:,:,:,i_age)) 
      else
        ocn%ts(:,:,:,i_age) = 0._wp
      endif
    endif

    if (dye_tracer) then
      if (nc_exists_var(fnm,"dye")) then
        call nc_read(fnm,"dye",    ocn%ts(:,:,:,i_dye)) 
      else
        ocn%ts(:,:,:,i_dye) = 0._wp
      endif
    endif

    call nc_read(fnm,"fw_glob",   ocn%fw_glob)

    print *,'read restart file ',fnm


   return

  end subroutine ocn_read_restart

end module ocn_model

