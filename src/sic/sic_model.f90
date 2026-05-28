!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s i c _ m o d e l
!
!  Purpose : main sea ice model
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
module sic_model

    use precision, only : wp
    use constants, only : T0, Lf, sigma
    use timer, only: sec_day, nday_year, doy, time_soy_sic
    use climber_grid, only: lon, lat, dlon, dlat
    use control, only: out_dir, restart_in_dir, sic_restart, check_water

    use sic_grid, only: maxi, maxj, area, area_old, dx, dy
    use sic_grid, only : sic_grid_init, sic_grid_update
    use sic_params, only : sic_params_init
    use sic_params, only: f_sic_min, h_sic_max, h_sic_min, h_snow_max, h_snow_min
    use sic_params, only : i_fsic, f_sic_max, f_sic_pow, h0, rho_sic, rho_snow, dt
    use sic_def, only : sic_class

    use ebal_sic_mod, only : ebal_sic
    use ebal_ocn_mod, only : ebal_ocn
    use surface_par_sic
    use sic_dyn_mod, only : sic_dyn
    use transport_sic_mod, only : transport
    use ncio

    !$use omp_lib

    implicit none

    private
    public :: sic_init, sic_update, sic_end
    public :: sic_read_restart, sic_write_restart

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s i c _ u p d a t e
  ! Purpose  :  update sea ice
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_update(sic)

    implicit none

    type(sic_class) :: sic

    integer :: i, j
    real(wp) :: t_freeze
    real(wp), dimension(maxi,maxj) :: flx_sic_cut, fw_sic_cut

    ! per-cell + global water-budget check (only used when check_water is on)
    real(wp), dimension(maxi,maxj) :: water_mass_old, fw_ocn_pre_redist
    real(wp) :: water_mass_new, flux_per_area, water_res, fw_cut_per_area
    real(wp) :: water_res_global
    real(wp) :: water_mass_old_tot, water_mass_new_tot, water_flux_tot, water_res_tot
    integer  :: n_cell_fail
    real(wp), parameter :: water_check_tol_sic = 1.e-6_wp  ! kg/m2 per step
    real(wp), parameter :: water_check_tol_sic_glob = 1.e4_wp  ! kg per step


    if (time_soy_sic) then
      ! update ocean area accounting for actual ocean fraction
      call sic_grid_update(sic%f_ocn)
      ! account for ocean area changes, kill dead cells and init new cells
      call sic_update_cells(sic)
    endif

    ! snapshot per-cell sea ice + snow water mass for the water-budget check
    if (check_water) then
      do j=1,maxj
        do i=1,maxi
          if (sic%f_ocn(i,j).gt.0._wp) then
            water_mass_old(i,j) = sic%h_sic_mean(i,j)*rho_sic + sic%h_snow_mean(i,j)*rho_snow  ! kg/m2
          else
            water_mass_old(i,j) = 0._wp
          endif
        enddo
      enddo
    endif

    !$omp parallel do collapse(2) private(i,j,t_freeze)
    do j=1,maxj
      do i=1,maxi

        ! compute albedo of snow on sea ice
        call snow_albedo(sic%f_sic(i,j), sic%w_snow(i,j), sic%w_snow_max(i,j), sic%t_skin_sic(i,j), &   ! in
          sic%snow(i,j), sic%dust_dep(i,j), sic%coszm(doy,j), &   ! in
          sic%alb_snow_vis_dir(i,j), sic%alb_snow_vis_dif(i,j), sic%alb_snow_nir_dir(i,j), sic%alb_snow_nir_dif(i,j), & ! out
          sic%snow_grain(i,j), sic%dust_con(i,j)) ! out

        ! compute sea ice and ocean albedoes
        call ocn_sic_albedo(sic%f_ocn(i,j), sic%f_sic(i,j), sic%h_snow(i,j), sic%t_skin_sic(i,j), sic%coszm(doy,j), &   ! in
          sic%alb_snow_vis_dir(i,j), sic%alb_snow_vis_dif(i,j), sic%alb_snow_nir_dir(i,j), sic%alb_snow_nir_dif(i,j), & ! in
          sic%alb_ocn_vis_dir(i,j), sic%alb_ocn_vis_dif(i,j), sic%alb_ocn_nir_dir(i,j), sic%alb_ocn_nir_dif(i,j), &    ! out
          sic%alb_sic_vis_dir(i,j), sic%alb_sic_vis_dif(i,j), sic%alb_sic_nir_dir(i,j), sic%alb_sic_nir_dif(i,j), &    ! out
          sic%albedo_sic(i,j), sic%albedo_ocn(i,j))    ! out

        ! compute drag coefficients over sea ice and ocean
        call cdrag(sic%f_ocn(i,j), sic%t_skin_ocn(i,j), sic%t_air_ocn(i,j), sic%t_skin_sic(i,j), sic%t_air_sic(i,j), sic%wind(i,j), &
          sic%rough_m_ocn(i,j), sic%rough_h_ocn(i,j), sic%Cde_ocn(i,j), sic%Cdh_ocn(i,j), &
          sic%rough_m_sic(i,j), sic%rough_h_sic(i,j), sic%Cde_sic(i,j), sic%Cdh_sic(i,j))

        if (sic%f_ocn(i,j).gt.0._wp) then

          ! freezing temperature below sea ice as function of surface ocean salinity, from Marsh et al 2011, from Millero 1978
          t_freeze = -0.0575_wp*sic%sss(i,j)+0.0017_wp*sic%sss(i,j)**1.5_wp-0.0002_wp*sic%sss(i,j)**2 + T0  ! K

          if (sic%f_sic(i,j) .gt. 0._wp) then
            ! solve energy balance equation over sea ice and compute skin temperature
            ! diagnose surface energy fluxes
            ! compute changes in sea ice and snow thickness (thermodynamic only) over sea ice part
            call ebal_sic(sic%sst(i,j), t_freeze, sic%t_air_sic(i,j), &
            sic%q_air_sic(i,j), sic%pressure(i,j), &
            sic%flx_swnet_sic(i,j), sic%flx_lwd_sic(i,j), sic%snow(i,j), sic%rain(i,j), sic%wind(i,j), sic%Cde_sic(i,j), sic%Cdh_sic(i,j), &
            sic%t_skin_sic(i,j), sic%h_snow(i,j), sic%h_sic(i,j), &
            sic%flx_ocn_sic(i,j), sic%flx_sh_sic(i,j), sic%flx_lwu_sic(i,j), sic%flx_lh_sic(i,j), &
            sic%fw_ocn_sic(i,j), sic%evp_sic(i,j), &
            sic%flx_melt_top(i,j), sic%flx_melt_bot(i,j), &
            sic%dh_snow(i,j), sic%dh_sic_sic(i,j), i, j)
          else
            sic%flx_ocn_sic(i,j) = 0._wp
            sic%flx_sh_sic(i,j) = 0._wp
            sic%flx_lwu_sic(i,j) = 0._wp
            sic%flx_lh_sic(i,j) = 0._wp
            sic%fw_ocn_sic(i,j) = 0._wp
            sic%evp_sic(i,j) = 0._wp
            sic%flx_melt_top(i,j) = 0._wp
            sic%flx_melt_bot(i,j) = 0._wp
            sic%dh_snow(i,j) = 0._wp
            sic%dh_sic_sic(i,j) = 0._wp
          endif

          if (sic%f_sic(i,j) .lt. 1._wp) then
            ! solve energy balance equation over ocean and compute skin temperature
            ! diagnose surface energy fluxes
            ! compute changes in sea ice thickness (thermodynamic only) over ocean part
            call ebal_ocn(sic%sst(i,j), t_freeze, sic%t_air_ocn(i,j), sic%q_air_ocn(i,j), sic%pressure(i,j), &
              sic%flx_swnet_ocn(i,j), sic%flx_lwd_ocn(i,j), sic%snow(i,j), sic%rain(i,j), sic%wind(i,j), sic%Cde_ocn(i,j), sic%Cdh_ocn(i,j), &
              sic%t_skin_ocn(i,j), &
              sic%flx_ocn_ocn(i,j), sic%flx_sh_ocn(i,j), sic%flx_lwu_ocn(i,j), sic%flx_lh_ocn(i,j), &
              sic%fw_ocn_ocn(i,j), sic%evp_ocn(i,j), &
              sic%dh_sic_ocn(i,j), i, j)
          else
            sic%flx_ocn_ocn(i,j) = 0._wp 
            sic%flx_sh_ocn(i,j) = 0._wp 
            sic%flx_lwu_ocn(i,j) = 0._wp
            sic%flx_lh_ocn(i,j) = 0._wp
            sic%fw_ocn_ocn(i,j) = 0._wp
            sic%evp_ocn(i,j) = 0._wp
            sic%dh_sic_ocn(i,j) = 0._wp 
          endif

        endif

        if (sic%f_ocn(i,j).gt.0._wp) then
          ! grid cell mean snow and ice thickness changes from thermodynamics 
          sic%dh_sic_dt_therm(i,j) = (sic%f_sic(i,j)*sic%dh_sic_sic(i,j) + (1._wp-sic%f_sic(i,j))*sic%dh_sic_ocn(i,j))/dt  ! m/s
          sic%dh_snow_dt_therm(i,j) = (sic%f_sic(i,j)*sic%dh_snow(i,j))/dt  ! m/s
          sic%h_sic_mean(i,j) = sic%h_sic_mean(i,j) + sic%dh_sic_dt_therm(i,j)*dt
          sic%h_snow_mean(i,j) = sic%h_snow_mean(i,j) + sic%dh_snow_dt_therm(i,j)*dt
        else
          sic%dh_sic_dt_therm(i,j)  = 0._wp 
          sic%dh_snow_dt_therm(i,j) = 0._wp 
          sic%h_sic_mean(i,j)       = 0._wp 
          sic%h_snow_mean(i,j)      = 0._wp 
        endif

        ! average fluxes (must use the f_sic that drove the dh_*_dt_therm / h_*_mean update above
        sic%flx_sh(i,j) = (sic%f_sic(i,j)*sic%flx_sh_sic(i,j) + (1._wp-sic%f_sic(i,j))*sic%flx_sh_ocn(i,j))  ! W/m2
        sic%flx_lh(i,j) = (sic%f_sic(i,j)*sic%flx_lh_sic(i,j) + (1._wp-sic%f_sic(i,j))*sic%flx_lh_ocn(i,j))  ! W/m2
        sic%flx_lwu(i,j) = (sic%f_sic(i,j)*sic%flx_lwu_sic(i,j) + (1._wp-sic%f_sic(i,j))*sic%flx_lwu_ocn(i,j))  ! W/m2
        sic%evp(i,j) = (sic%f_sic(i,j)*sic%evp_sic(i,j) + (1._wp-sic%f_sic(i,j))*sic%evp_ocn(i,j))  ! W/m2
        sic%flx_ocn(i,j) = (sic%f_sic(i,j)*sic%flx_ocn_sic(i,j) + (1._wp-sic%f_sic(i,j))*sic%flx_ocn_ocn(i,j))  ! W/m2
        sic%fw_ocn(i,j)  = (sic%f_sic(i,j)*sic%fw_ocn_sic(i,j)  + (1._wp-sic%f_sic(i,j))*sic%fw_ocn_ocn(i,j))   ! kg/m2/s

        ! update sea ice fraction, thermodynamic contribution only
        if (i_fsic.eq.2 .or. i_fsic.eq.3) then
          ! prognostic sea ice fraction following Fichefet et al, 1997, and MPIOMv1.1 manual
          if (sic%f_ocn(i,j).gt.0._wp) then
            if (sic%dh_sic_ocn(i,j).gt.0._wp) then
              ! growing sea ice
              if (i_fsic.eq.2) then
                sic%f_sic(i,j) = sic%f_sic(i,j) + sqrt(1._wp-sic%f_sic(i,j)**2)*(1._wp-sic%f_sic(i,j)) * sic%dh_sic_ocn(i,j)/h0
              else if (i_fsic.eq.3) then
                sic%f_sic(i,j) = sic%f_sic(i,j) + (1._wp-sic%f_sic(i,j)) * sic%dh_sic_ocn(i,j)/h0
              endif
            endif
            if (sic%h_sic(i,j).gt.0._wp .and. sic%dh_sic_sic(i,j).lt.0._wp) then
              ! melting sea ice, fraction reduction assuming that ice thickness is
              ! distributed (linearly) between 0 and 2*h_sic
              sic%f_sic(i,j) = sic%f_sic(i,j) + sic%f_sic(i,j)*(sic%dh_sic_sic(i,j))/(2._wp*sic%h_sic(i,j))
            endif
            sic%f_sic(i,j) = max(0._wp,sic%f_sic(i,j))
            sic%f_sic(i,j) = min(f_sic_max,sic%f_sic(i,j))
          else
            sic%f_sic(i,j) = 0._wp
          endif
        endif

      enddo
    enddo
    !$omp end parallel do

    ! sea ice dynamics on c-grid, Hunke & Dukowics 1997 as implemented in SIS2 at GFDL 
    call sic_dyn(dt, sic%f_sic, sic%h_sic_mean, sic%h_snow_mean, sic%uo, sic%vo, sic%tauxa, sic%tauya, sic%ssh, &   ! in
      sic%u, sic%v, sic%str_d, sic%str_t, sic%str_s, &        ! inout
      sic%tauxo, sic%tauyo, sic%PFu, sic%PFv, sic%Cor_u, sic%Cor_v, & ! out
      sic%fxic, sic%fxic_d, sic%fxic_t, sic%fxic_s, sic%fyic, sic%fyic_d, sic%fyic_t, sic%fyic_s) ! out

    ! integrate advection-diffusion equation for sea ice thickness, snow thickness and sea ice concentration

    !$omp parallel sections
    !$omp section
    ! transport sea ice thickness
    call transport(sic%u,sic%v,sic%f_ocn, &
      sic%h_sic_mean,sic%dh_sic_dt_dyn, &
      sic%fax_sic,sic%fay_sic,sic%fdx_sic,sic%fdy_sic)
    !$omp section
    ! transport snow thickness
    call transport(sic%u,sic%v,sic%f_ocn, &
      sic%h_snow_mean,sic%dh_snow_dt_dyn, &
      sic%fax_snow,sic%fay_snow,sic%fdx_snow,sic%fdy_snow)
    !$omp section
    if (i_fsic.eq.2 .or. i_fsic.eq.3) then
      ! transport sea ice concentration
      call transport(sic%u,sic%v,sic%f_ocn, &
        sic%f_sic)
    endif
    !$omp end parallel sections

    !$omp parallel do collapse(2) private(i,j) 
    do j=1,maxj
      do i=1,maxi

        flx_sic_cut(i,j) = 0._wp
        fw_sic_cut(i,j) = 0._wp

        ! limit sea ice concentration
        sic%f_sic(i,j) = max(0._wp,sic%f_sic(i,j))
        sic%f_sic(i,j) = min(f_sic_max,sic%f_sic(i,j))

        ! set sea ice concentration to 1 when threshold thickness crossed
        if (sic%h_sic_mean(i,j)>h_sic_max/2._wp) sic%f_sic(i,j) = 1._wp 

        ! snow to ice conversion: when snow on ice exceeds h_snow_max, convert
        ! the excess (per ice area) into sea ice. 
        if (sic%f_sic(i,j) .gt. 0._wp .and. sic%h_snow_mean(i,j) .gt. sic%f_sic(i,j)*h_snow_max) then
          sic%h_sic_mean(i,j)  = sic%h_sic_mean(i,j) &
            + (sic%h_snow_mean(i,j) - sic%f_sic(i,j)*h_snow_max)*rho_snow/rho_sic
          sic%h_snow_mean(i,j) = sic%f_sic(i,j)*h_snow_max
        endif

        ! warn for unphysically large negative grid-mean thicknesses
        if (sic%h_sic_mean(i,j).lt.-1e-1_wp) then
          print *
          print *,'warning, sea ice thickness < 0'
          print *,'h_sic_mean,fsic',sic%h_sic_mean(i,j),sic%f_sic(i,j),i,j,sic%f_ocn(i,j)
          print *,'dhis,dhis_sic,dhis_ocn',sic%dh_sic_dt_therm(i,j)*dt,sic%dh_sic_sic(i,j),sic%dh_sic_ocn(i,j)
        endif
        if (sic%h_snow_mean(i,j).lt.-1e-1_wp) then
          print *
          print *,'warning, snow thickness < 0'
          print *,'h_snow_mean,fsic',sic%h_snow_mean(i,j),sic%f_sic(i,j),i,j,sic%f_ocn(i,j)
          print *,'dhsnow',sic%dh_snow_dt_therm(i,j)*dt
        endif

        if (sic%f_ocn(i,j).gt.0._wp) then
          ! grid average freshwater flux related to brine rejection (negative)
          sic%fw_brines(i,j) = -max(0._wp,sic%dh_sic_dt_therm(i,j)*rho_sic)  ! kg/m2/s
        else
          sic%flx_ocn(i,j) = 0._wp
          sic%fw_ocn(i,j)  = 0._wp
          sic%fw_brines(i,j)  = 0._wp
        endif

        ! limit snow and sea ice thickness
        if (sic%f_ocn(i,j).gt.0._wp) then
          ! if sea ice thickness < h_sic_min then remove sea ice and snow and update fluxes
          if (sic%h_sic_mean(i,j) .lt. h_sic_min) then
            sic%flx_ocn(i,j) = sic%flx_ocn(i,j) - sic%h_sic_mean(i,j)*rho_sic*Lf/dt - sic%h_snow_mean(i,j)*rho_snow*Lf/dt  ! W/m2
            sic%fw_ocn(i,j)  = sic%fw_ocn(i,j)  + sic%h_sic_mean(i,j)*rho_sic/dt    + sic%h_snow_mean(i,j)*rho_snow/dt  ! kg/m2/s
            sic%h_sic_mean(i,j) = 0._wp
            sic%h_snow_mean(i,j) = 0._wp
          endif
          ! if sea ice thickness > h_sic_max then remove excess sea ice and update fluxes
          if (sic%h_sic_mean(i,j) .gt. h_sic_max) then
            flx_sic_cut(i,j) = - (sic%h_sic_mean(i,j)-h_sic_max)*rho_sic*Lf/dt*area(i,j)  ! W
            fw_sic_cut(i,j)  =   (sic%h_sic_mean(i,j)-h_sic_max)*rho_sic/dt*area(i,j)   ! kg/s
            sic%h_sic_mean(i,j) = h_sic_max
          endif
          ! if snow thickness < h_snow_min then remove snow and update fluxes
          if (sic%h_snow_mean(i,j) .lt. h_snow_min) then
            sic%flx_ocn(i,j) = sic%flx_ocn(i,j) - sic%h_snow_mean(i,j)*rho_snow*Lf/dt  ! W
            sic%fw_ocn(i,j)  = sic%fw_ocn(i,j)  + sic%h_snow_mean(i,j)*rho_snow/dt  ! kg/s
            sic%h_snow_mean(i,j) = 0._wp
          endif
          ! if snow thickness > h_snow_max then remove excess snow and update fluxes
          if (sic%h_snow_mean(i,j) .gt. h_snow_max) then
            sic%flx_ocn(i,j) = sic%flx_ocn(i,j) - (sic%h_snow_mean(i,j)-h_snow_max)*rho_snow*Lf/dt  ! W
            sic%fw_ocn(i,j)  = sic%fw_ocn(i,j)  + (sic%h_snow_mean(i,j)-h_snow_max)*rho_snow/dt  ! kg/s
            sic%h_snow_mean(i,j) = h_snow_max
          endif
        endif

        ! compute fractional sea ice area
        if (i_fsic.eq.1) then
          ! thickness-area relation from CLIMBER-2
          if (sic%f_ocn(i,j).gt.0._wp) then
            if (sic%h_sic_mean(i,j) .ge. h_sic_min) then
              sic%f_sic(i,j) = (sic%h_sic_mean(i,j)/2._wp)**f_sic_pow
            else
              sic%f_sic(i,j) = 0._wp
            endif
          else
            sic%f_sic(i,j) = 0._wp
          endif
          sic%f_sic(i,j) = min(f_sic_max,sic%f_sic(i,j))
        endif

        ! cut small sea ice concentrations 
        if (sic%f_ocn(i,j).gt.0._wp) then
          ! if sea ice fraction < f_sic_min then remove sea ice and snow and update fluxes
          if (sic%f_sic(i,j) .lt. f_sic_min) then
            sic%flx_ocn(i,j) = sic%flx_ocn(i,j) - sic%h_sic_mean(i,j)*rho_sic*Lf/dt - sic%h_snow_mean(i,j)*rho_snow*Lf/dt  ! W/m2
            sic%fw_ocn(i,j)  = sic%fw_ocn(i,j)  + sic%h_sic_mean(i,j)*rho_sic/dt    + sic%h_snow_mean(i,j)*rho_snow/dt  ! kg/m2/s
            sic%h_sic_mean(i,j) = 0._wp
            sic%h_snow_mean(i,j) = 0._wp
            sic%f_sic(i,j) = 0._wp
          endif
        endif

        ! save seasonal maximum snow water equivalent
        if (sic%w_snow(i,j).gt.sic%h_snow(i,j)*rho_snow) sic%w_snow_max(i,j) = sic%w_snow(i,j)

        ! compute sea ice and snow thickness of sea ice covered area
        if (sic%f_ocn(i,j).gt.0._wp .and. sic%f_sic(i,j).gt.0._wp) then
          sic%h_sic(i,j)  = sic%h_sic_mean(i,j)/sic%f_sic(i,j)
          sic%h_snow(i,j) = sic%h_snow_mean(i,j)/sic%f_sic(i,j)
          sic%w_snow(i,j) = sic%h_snow(i,j)*rho_snow    ! kg/m2
        else
          sic%h_sic(i,j)  = 0._wp
          sic%h_snow(i,j) = 0._wp
          sic%w_snow(i,j) = 0._wp
        endif

        if (sic%h_sic(i,j).eq.0._wp) then
          sic%f_sic(i,j) = 0._wp
        endif

        ! check for NaNs
        if (sic%u(i,j).ne.sic%u(i,j)) sic%error = .true.
        if (sic%v(i,j).ne.sic%v(i,j)) sic%error = .true.
        if (sic%f_sic(i,j).ne.sic%f_sic(i,j)) sic%error = .true.
        if (sic%h_sic_mean(i,j).ne.sic%h_sic_mean(i,j)) sic%error = .true.
        if (sic%h_snow_mean(i,j).ne.sic%h_snow_mean(i,j)) sic%error = .true.
      enddo
    enddo
    !$omp end parallel do

    ! save fw_ocn before global redistribution of the h_sic_max cut, so the
    ! per-cell water-budget check below sees the actual local fw_ocn 
    if (check_water) fw_ocn_pre_redist = sic%fw_ocn

    ! distribute heat flux from cutting sea ice/snow over the whole ocean surface
    sic%flx_ocn(:,:) = sic%flx_ocn(:,:) + sum(flx_sic_cut)/sum(area) ! W/m2
    ! distribute freshwater flux from cutting sea ice/snow over the whole ocean surface
    sic%fw_ocn(:,:)  = sic%fw_ocn(:,:)  + sum(fw_sic_cut)/sum(area) ! kg/m2/s

    ! per-cell water-budget check
    ! the budget per unit ocean area, over a single sic_update step, is
    !   ( m_new - m_old ) = ( rain + snow - evp - fw_ocn_local
    !                       + dh_sic_dt_dyn*rho_sic + dh_snow_dt_dyn*rho_snow
    !                       - fw_sic_cut/area ) * dt
    ! where m = h_sic_mean*rho_sic + h_snow_mean*rho_snow [kg/m2].
    if (check_water) then
      water_res_global = 0._wp
      water_mass_old_tot = 0._wp
      water_mass_new_tot = 0._wp
      water_flux_tot     = 0._wp
      n_cell_fail = 0
      do j=1,maxj
        do i=1,maxi
          if (sic%f_ocn(i,j).gt.0._wp) then
            water_mass_new = sic%h_sic_mean(i,j)*rho_sic + sic%h_snow_mean(i,j)*rho_snow  ! kg/m2
            fw_cut_per_area = fw_sic_cut(i,j)/area(i,j)  ! kg/m2/s; nonzero only where h_sic_max cut
            flux_per_area = sic%snow(i,j) + sic%rain(i,j) - sic%evp(i,j) - fw_ocn_pre_redist(i,j) &
              + sic%dh_sic_dt_dyn(i,j)*rho_sic + sic%dh_snow_dt_dyn(i,j)*rho_snow &
              - fw_cut_per_area  ! kg/m2/s
            water_res = (water_mass_new - water_mass_old(i,j)) - flux_per_area*dt  ! kg/m2
            water_res_global = water_res_global + water_res*area(i,j)
            ! independent global accumulators
            ! d(total water)/dt = sum( (rain + snow - evp - fw_ocn) * area )
            water_mass_old_tot = water_mass_old_tot + water_mass_old(i,j)*area(i,j)
            water_mass_new_tot = water_mass_new_tot + water_mass_new*area(i,j)
            water_flux_tot     = water_flux_tot &
              + (sic%snow(i,j) + sic%rain(i,j) - sic%evp(i,j) - sic%fw_ocn(i,j))*area(i,j)*dt
            if (abs(water_res).gt.water_check_tol_sic) then
              n_cell_fail = n_cell_fail + 1
              if (n_cell_fail.le.5) then
                print *,''
                print *,'WARNING: sic per-cell water budget residual ',water_res,' kg/m2 at (i,j) ',i,j
                print *,'  m_old, m_new                  [kg/m2]    = ',water_mass_old(i,j), water_mass_new
                print *,'  (rain+snow)*dt                [kg/m2]    = ',(sic%snow(i,j)+sic%rain(i,j))*dt
                print *,'  evp*dt                        [kg/m2]    = ',sic%evp(i,j)*dt
                print *,'  fw_ocn_local*dt               [kg/m2]    = ',fw_ocn_pre_redist(i,j)*dt
                print *,'  fw_sic_cut/area*dt            [kg/m2]    = ',fw_cut_per_area*dt
                print *,'  transport_div*dt              [kg/m2]    = ', &
                  (sic%dh_sic_dt_dyn(i,j)*rho_sic + sic%dh_snow_dt_dyn(i,j)*rho_snow)*dt
                print *,'  f_ocn, f_sic                             = ',sic%f_ocn(i,j), sic%f_sic(i,j)
              endif
            endif
          endif
        enddo
      enddo
      if (n_cell_fail.gt.0) then
        print *,'WARNING: sic per-cell water budget failed in ',n_cell_fail,' cell(s); ', &
          'domain-integrated residual ',water_res_global,' kg per step'
      endif
      ! independent global check: storage change should equal integrated
      ! atmospheric + ocean fluxes (transport is domain-conservative). 
      water_res_tot = (water_mass_new_tot - water_mass_old_tot) - water_flux_tot
      if (abs(water_res_tot).gt.water_check_tol_sic_glob) then
        print *,''
        print *,'WARNING: sic global water budget residual ',water_res_tot,' kg per step'
        print *,'  mass_old_total [kg]      = ',water_mass_old_tot
        print *,'  mass_new_total [kg]      = ',water_mass_new_tot
        print *,'  flux_total*dt  [kg]      = ',water_flux_tot
      endif
    endif


   return

  end subroutine sic_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s i c _ u p d a t e _ c e l l s
  ! Purpose  :  update sea ice cells
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_update_cells(sic)

    implicit none

    type(sic_class) :: sic

    integer :: i, j, ii, jj, iii, jjj
    real(wp) :: ice_excess, snow_excess, area_neighbour


    ! update sea ice and snow thickness to conserve ice/snow accounting for ocean area changes
    do j=1,maxj
      do i=1,maxi
        if (area_old(i,j).gt.area(i,j) .and. area(i,j).gt.0._wp) then
          ! cell area is shrinking
          sic%h_sic_mean(i,j) = sic%h_sic_mean(i,j)*area_old(i,j)/area(i,j)
          sic%h_snow_mean(i,j) = sic%h_snow_mean(i,j)*area_old(i,j)/area(i,j)
        endif
        if (area_old(i,j).lt.area(i,j) .and. area_old(i,j).gt.0._wp) then
          ! cell area is increasing
          sic%h_sic_mean(i,j) = sic%h_sic_mean(i,j)*area_old(i,j)/area(i,j)
          sic%h_snow_mean(i,j) = sic%h_snow_mean(i,j)*area_old(i,j)/area(i,j)
        endif
      enddo
    enddo
    ! initialize newly formed cells, partly using neighbouring cells? , todo
    do j=1,maxj
      do i=1,maxi
        sic%mask_new_cell(i,j) = 0
        if (area_old(i,j).eq.0._wp .and. area(i,j).gt.0._wp) then
          sic%mask_new_cell(i,j) = 1
          sic%f_sic(i,j)      = 0._wp
          sic%h_sic_mean(i,j) = 0._wp
          sic%h_sic(i,j)      = 0._wp
          sic%h_snow_mean(i,j) = 0._wp
          sic%h_snow(i,j)     = 0._wp
          sic%w_snow(i,j)     = 0._wp
          sic%w_snow_max(i,j) = 0._wp
          sic%t_skin_sic(i,j) = T0
          sic%t_skin_ocn(i,j) = T0
          ! set also ocean variables which are not defined yet because
          ! sic_model is called before ocn_model
          sic%sst(i,j) = T0
          sic%t_skin_sic(i,j) = T0
          sic%t_skin_ocn(i,j) = T0
          ! set also ocean variables which are not defined yet because
          ! sic_model is called before ocn_model
          sic%sst(i,j) = T0
        endif
      enddo
    enddo
    ! end cell and move ice/snow to neighbours to conserve mass
    do j=1,maxj
      do i=1,maxi
        if (area_old(i,j).gt.0._wp .and. area(i,j).eq.0._wp) then
          ! compute ice and snow volume of dead grid cell
          ice_excess = sic%h_sic_mean(i,j)*area_old(i,j) ! m3
          snow_excess = sic%h_snow_mean(i,j)*area_old(i,j) ! m3
          ! remove ice and snow from dead cell
          sic%h_sic_mean(i,j) = 0._wp
          sic%h_snow_mean(i,j) = 0._wp
          ! move excess ice to neighbouring cells
          area_neighbour = 0._wp
          do ii=i-1,i+1
            do jj=j-1,j+1
              iii = ii
              if (iii.eq.0) iii = maxi
              if (iii.eq.maxi+1) iii = 1
              jjj = jj
              jjj = max(1,jjj)
              jjj = min(maxj,jjj)
              if (area(iii,jjj).gt.0._wp) then
                area_neighbour = area_neighbour + area(iii,jjj)
              endif
            enddo
          enddo
          if (area_neighbour.gt.0._wp) then
            do ii=i-1,i+1
              do jj=j-1,j+1
                iii = ii
                if (iii.eq.0) iii = maxi
                if (iii.eq.maxi+1) iii = 1
                jjj = jj
                jjj = max(1,jjj)
                jjj = min(maxj,jjj)
                if (area(iii,jjj).gt.0._wp) then
                  sic%h_sic_mean(iii,jjj) = sic%h_sic_mean(iii,jjj) + ice_excess/area_neighbour
                  sic%h_snow_mean(iii,jjj) = sic%h_snow_mean(iii,jjj) + snow_excess/area_neighbour
                endif
              enddo
            enddo
          else
            ! no neighbours, distribute ice and snow to all cells with sea ice
            area_neighbour = sum(area,mask=sic%f_sic.gt.0._wp)
            where (sic%f_sic.gt.0._wp .and. area.gt.0._wp)
              sic%h_sic_mean = sic%h_sic_mean + ice_excess/area_neighbour
              sic%h_snow_mean = sic%h_snow_mean + snow_excess/area_neighbour
            endwhere
          endif
        endif
      enddo
    enddo

    return

  end subroutine sic_update_cells


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s i c _ i n i t
  ! Purpose  :  initialize sea ice
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_init(sic, f_ocn, dz1)

    implicit none

    type(sic_class) :: sic
    real(wp), dimension(:,:), intent(in) :: f_ocn
    real(wp), intent(in) :: dz1


    ! initialize sea ice model grid
    call sic_grid_init(f_ocn, dz1)

    ! initialize/read sea ice model parameters
    call sic_params_init

    ! allocate 
    call sic_alloc(sic)

    if (sic_restart) then
      ! read restart file
      call sic_read_restart(trim(restart_in_dir)//"/sic_restart.nc",sic)
    else
      ! initialise prognostic variables
      sic%u          = 0._wp
      sic%v          = 0._wp
      sic%h_sic_mean = 0._wp
      sic%h_sic      = 0._wp
      sic%f_sic      = 0._wp
      sic%h_snow_mean = 0._wp
      sic%h_snow     = 0._wp
      sic%w_snow     = 0._wp
      sic%w_snow_max = 0._wp
      sic%str_d      = 0._wp
      sic%str_t      = 0._wp
      sic%str_s      = 0._wp
      sic%t_skin_sic = T0
      sic%t_skin_ocn = T0
      sic%sst        = T0
      sic%albedo_sic = 0.1_wp
      sic%albedo_ocn = 0.1_wp
      sic%alb_sic_vis_dir = 0.1_wp  
      sic%alb_sic_vis_dif = 0.1_wp
      sic%alb_sic_nir_dir = 0.1_wp
      sic%alb_sic_nir_dif = 0.1_wp
      sic%alb_ocn_vis_dir = 0.1_wp
      sic%alb_ocn_vis_dif = 0.1_wp
      sic%alb_ocn_nir_dir = 0.1_wp
      sic%alb_ocn_nir_dif = 0.1_wp
      sic%rough_m_ocn = 1.e-3_wp  
      sic%rough_m_sic = 1.e-3_wp
      sic%rough_h_ocn = 1.e-3_wp
      sic%rough_h_sic = 1.e-3_wp
      sic%Cdh_ocn     = 1.e-3_wp
      sic%Cdh_sic     = 1.e-3_wp 
      sic%flx_sh  = 0._wp
      sic%flx_lwu = sigma*T0**4
      sic%flx_lh  = 0._wp
      sic%evp     = 0._wp
      sic%flx_sh_sic = 0._wp
      sic%flx_lwu_sic = 0._wp
      sic%flx_lh_sic = 0._wp
      sic%evp_sic = 0._wp
      sic%flx_sh_ocn = 0._wp
      sic%flx_lwu_ocn = 0._wp
      sic%flx_lh_ocn = 0._wp
      sic%evp_ocn = 0._wp
    endif

    ! other initialisations
    sic%dh_sic_dt_dyn = 0._wp
    sic%dh_sic_dt_therm = 0._wp
    sic%dh_sic_sic = 0._wp
    sic%dh_sic_ocn = 0._wp
    sic%dh_snow_dt_dyn = 0._wp
    sic%dh_snow_dt_therm = 0._wp
    sic%dh_snow = 0._wp

    sic%flx_ocn_sic = 0._wp
    sic%fw_ocn_sic = 0._wp
    sic%flx_ocn_ocn = 0._wp
    sic%fw_ocn_ocn = 0._wp

    sic%error = .false.

    print*
    print*,'======================================================='
    print*,' Initialisation of SEA ICE complete'
    print*,'======================================================='
    print*


    return

  end subroutine sic_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s i c _ end
  ! Purpose  :  end sea ice
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_end(sic)

    implicit none

    type(sic_class) :: sic


    ! Deallocate all state variables to free memory
    call sic_dealloc(sic)


    return

  end subroutine sic_end


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s i c _ a l l o c
  ! Purpose  :  allocate sea ice
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_alloc(sic)

    implicit none 

    type(sic_class) :: sic

    ! Allocate all state variables with dimensions

    allocate(sic%mask_new_cell(maxi,maxj))
    allocate(sic%f_ocn(maxi,maxj))
    allocate(sic%u(maxi,maxj))
    allocate(sic%v(maxi,maxj))
    allocate(sic%uo(maxi,maxj))
    allocate(sic%vo(maxi,maxj))
    allocate(sic%tauxa(maxi,maxj))
    allocate(sic%tauya(maxi,maxj))
    allocate(sic%tauxo(maxi,maxj))
    allocate(sic%tauyo(maxi,maxj))
    allocate(sic%ssh(maxi,maxj))
    allocate(sic%str_d(maxi,maxj))
    allocate(sic%str_t(maxi,maxj))
    allocate(sic%str_s(maxi,maxj))
    allocate(sic%fxic  (maxi,maxj)) 
    allocate(sic%fxic_d(maxi,maxj)) 
    allocate(sic%fxic_t(maxi,maxj)) 
    allocate(sic%fxic_s(maxi,maxj)) 
    allocate(sic%Cor_u (maxi,maxj)) 
    allocate(sic%PFu   (maxi,maxj)) 
    allocate(sic%fyic  (maxi,maxj)) 
    allocate(sic%fyic_d(maxi,maxj)) 
    allocate(sic%fyic_t(maxi,maxj)) 
    allocate(sic%fyic_s(maxi,maxj)) 
    allocate(sic%Cor_v (maxi,maxj)) 
    allocate(sic%PFv   (maxi,maxj)) 
    allocate(sic%h_sic_mean(maxi,maxj))
    allocate(sic%h_sic(maxi,maxj))
    allocate(sic%f_sic(maxi,maxj))
    allocate(sic%dh_sic_dt_dyn(maxi,maxj))
    allocate(sic%dh_sic_dt_therm(maxi,maxj))
    allocate(sic%dh_sic_sic(maxi,maxj))
    allocate(sic%dh_sic_ocn(maxi,maxj))
    allocate(sic%fax_sic(0:maxi,maxj))
    allocate(sic%fay_sic(maxi,0:maxj))
    allocate(sic%fdx_sic(0:maxi,maxj))
    allocate(sic%fdy_sic(maxi,0:maxj))
    allocate(sic%fax_snow(0:maxi,maxj))
    allocate(sic%fay_snow(maxi,0:maxj))
    allocate(sic%fdx_snow(0:maxi,maxj))
    allocate(sic%fdy_snow(maxi,0:maxj))
    allocate(sic%h_snow_mean(maxi,maxj))    
    allocate(sic%h_snow(maxi,maxj))
    allocate(sic%w_snow(maxi,maxj))
    allocate(sic%w_snow_max(maxi,maxj))
    allocate(sic%dh_snow_dt_dyn(maxi,maxj))
    allocate(sic%dh_snow_dt_therm(maxi,maxj))
    allocate(sic%dh_snow(maxi,maxj))
    allocate(sic%flx_melt_top(maxi,maxj))
    allocate(sic%flx_melt_bot(maxi,maxj))
    allocate(sic%t_skin_sic(maxi,maxj))
    allocate(sic%t_skin_ocn(maxi,maxj))
    allocate(sic%sst(maxi,maxj))
    allocate(sic%sss(maxi,maxj))
    allocate(sic%flx_lwd_sic(maxi,maxj))
    allocate(sic%flx_lwd_ocn(maxi,maxj))
    allocate(sic%flx_swnet_sic(maxi,maxj))
    allocate(sic%flx_swnet_ocn(maxi,maxj))
    allocate(sic%wind(maxi,maxj))
    allocate(sic%pressure(maxi,maxj))
    allocate(sic%t_air_sic(maxi,maxj))
    allocate(sic%q_air_sic(maxi,maxj))
    allocate(sic%t_air_ocn(maxi,maxj))
    allocate(sic%q_air_ocn(maxi,maxj))
    allocate(sic%snow(maxi,maxj))
    allocate(sic%rain(maxi,maxj))
    allocate(sic%flx_ocn(maxi,maxj))
    allocate(sic%fw_ocn(maxi,maxj))
    allocate(sic%fw_brines(maxi,maxj))
    allocate(sic%flx_sh(maxi,maxj))
    allocate(sic%flx_lh(maxi,maxj))
    allocate(sic%flx_lwu(maxi,maxj))
    allocate(sic%evp(maxi,maxj))
    allocate(sic%flx_sh_sic(maxi,maxj))
    allocate(sic%flx_lh_sic(maxi,maxj))
    allocate(sic%flx_lwu_sic(maxi,maxj))
    allocate(sic%evp_sic(maxi,maxj))
    allocate(sic%flx_ocn_sic(maxi,maxj))
    allocate(sic%fw_ocn_sic(maxi,maxj))
    allocate(sic%flx_sh_ocn(maxi,maxj))
    allocate(sic%flx_lh_ocn(maxi,maxj))
    allocate(sic%flx_lwu_ocn(maxi,maxj))
    allocate(sic%evp_ocn(maxi,maxj))
    allocate(sic%flx_ocn_ocn(maxi,maxj))
    allocate(sic%fw_ocn_ocn(maxi,maxj))
    allocate(sic%snow_grain(maxi,maxj))
    allocate(sic%dust_dep(maxi,maxj))
    allocate(sic%dust_con(maxi,maxj))
    allocate(sic%alb_snow_vis_dir(maxi,maxj))
    allocate(sic%alb_snow_vis_dif(maxi,maxj))
    allocate(sic%alb_snow_nir_dir(maxi,maxj))
    allocate(sic%alb_snow_nir_dif(maxi,maxj))
    allocate(sic%alb_ocn_vis_dir(maxi,maxj))
    allocate(sic%alb_ocn_vis_dif(maxi,maxj))
    allocate(sic%alb_ocn_nir_dir(maxi,maxj))
    allocate(sic%alb_ocn_nir_dif(maxi,maxj))
    allocate(sic%alb_sic_vis_dir(maxi,maxj))
    allocate(sic%alb_sic_vis_dif(maxi,maxj))
    allocate(sic%alb_sic_nir_dir(maxi,maxj))
    allocate(sic%alb_sic_nir_dif(maxi,maxj))
    allocate(sic%albedo_ocn(maxi,maxj))
    allocate(sic%albedo_sic(maxi,maxj))
    allocate(sic%rough_m_ocn(maxi,maxj))
    allocate(sic%rough_h_ocn(maxi,maxj))
    allocate(sic%rough_m_sic(maxi,maxj))
    allocate(sic%rough_h_sic(maxi,maxj))
    allocate(sic%Cde_ocn(maxi,maxj))
    allocate(sic%Cdh_ocn(maxi,maxj))
    allocate(sic%Cde_sic(maxi,maxj))
    allocate(sic%Cdh_sic(maxi,maxj))
    allocate(sic%coszm(nday_year,maxj))


    return

  end subroutine sic_alloc 


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s i c _ d e a l l o c
  ! Purpose  :  deallocate sea ice
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_dealloc(sic)

    implicit none 

    type(sic_class) :: sic


    ! Deallocate all state variables to free memory
    deallocate(sic%mask_new_cell)
    deallocate(sic%f_ocn)
    deallocate(sic%u)
    deallocate(sic%v)
    deallocate(sic%uo)
    deallocate(sic%vo)
    deallocate(sic%tauxa)
    deallocate(sic%tauya)
    deallocate(sic%tauxo)
    deallocate(sic%tauyo)
    deallocate(sic%ssh)
    deallocate(sic%str_d)
    deallocate(sic%str_t)
    deallocate(sic%str_s)
    deallocate(sic%fxic  ) 
    deallocate(sic%fxic_d) 
    deallocate(sic%fxic_t) 
    deallocate(sic%fxic_s) 
    deallocate(sic%Cor_u ) 
    deallocate(sic%PFu   ) 
    deallocate(sic%fyic  ) 
    deallocate(sic%fyic_d) 
    deallocate(sic%fyic_t) 
    deallocate(sic%fyic_s) 
    deallocate(sic%Cor_v ) 
    deallocate(sic%PFv   ) 
    deallocate(sic%h_sic_mean)
    deallocate(sic%h_sic)
    deallocate(sic%f_sic)
    deallocate(sic%dh_sic_dt_dyn)
    deallocate(sic%dh_sic_dt_therm)
    deallocate(sic%dh_sic_sic)
    deallocate(sic%dh_sic_ocn)
    deallocate(sic%fax_sic)
    deallocate(sic%fay_sic)
    deallocate(sic%fdx_sic)
    deallocate(sic%fdy_sic)
    deallocate(sic%fax_snow)
    deallocate(sic%fay_snow)
    deallocate(sic%fdx_snow)
    deallocate(sic%fdy_snow)
    deallocate(sic%h_snow_mean)
    deallocate(sic%h_snow)
    deallocate(sic%w_snow)
    deallocate(sic%w_snow_max)
    deallocate(sic%dh_snow_dt_dyn)
    deallocate(sic%dh_snow_dt_therm)
    deallocate(sic%dh_snow)
    deallocate(sic%flx_melt_top)
    deallocate(sic%flx_melt_bot)
    deallocate(sic%t_skin_sic)
    deallocate(sic%t_skin_ocn)
    deallocate(sic%sst)
    deallocate(sic%sss)
    deallocate(sic%flx_lwd_sic)
    deallocate(sic%flx_lwd_ocn)
    deallocate(sic%flx_swnet_sic)
    deallocate(sic%flx_swnet_ocn)
    deallocate(sic%wind)
    deallocate(sic%pressure)
    deallocate(sic%t_air_sic)
    deallocate(sic%q_air_sic)
    deallocate(sic%t_air_ocn)
    deallocate(sic%q_air_ocn)
    deallocate(sic%snow)
    deallocate(sic%rain)
    deallocate(sic%flx_ocn)
    deallocate(sic%fw_ocn)
    deallocate(sic%fw_brines)
    deallocate(sic%flx_sh)
    deallocate(sic%flx_lh)
    deallocate(sic%flx_lwu)
    deallocate(sic%evp)
    deallocate(sic%flx_sh_sic)
    deallocate(sic%flx_lh_sic)
    deallocate(sic%flx_lwu_sic)
    deallocate(sic%evp_sic)
    deallocate(sic%flx_ocn_sic)
    deallocate(sic%fw_ocn_sic)
    deallocate(sic%flx_sh_ocn)
    deallocate(sic%flx_lh_ocn)
    deallocate(sic%flx_lwu_ocn)
    deallocate(sic%evp_ocn)
    deallocate(sic%flx_ocn_ocn)
    deallocate(sic%fw_ocn_ocn)
    deallocate(sic%snow_grain)
    deallocate(sic%dust_dep)
    deallocate(sic%dust_con)
    deallocate(sic%alb_snow_vis_dir)
    deallocate(sic%alb_snow_vis_dif)
    deallocate(sic%alb_snow_nir_dir)
    deallocate(sic%alb_snow_nir_dif)
    deallocate(sic%alb_ocn_vis_dir)
    deallocate(sic%alb_ocn_vis_dif)
    deallocate(sic%alb_ocn_nir_dir)
    deallocate(sic%alb_ocn_nir_dif)
    deallocate(sic%alb_sic_vis_dir)
    deallocate(sic%alb_sic_vis_dif)
    deallocate(sic%alb_sic_nir_dir)
    deallocate(sic%alb_sic_nir_dif)
    deallocate(sic%albedo_ocn)
    deallocate(sic%albedo_sic)
    deallocate(sic%rough_m_ocn)
    deallocate(sic%rough_h_ocn)
    deallocate(sic%rough_m_sic)
    deallocate(sic%rough_h_sic)
    deallocate(sic%Cde_ocn)
    deallocate(sic%Cdh_ocn)
    deallocate(sic%Cde_sic)
    deallocate(sic%Cdh_sic)
    deallocate(sic%coszm)


    return

  end subroutine sic_dealloc 


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s i c _ w r i t e _ r e s t a r t
  ! Purpose  :  Write restart netcdf file 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_write_restart(fnm,sic)

    implicit none

    character (len=*) :: fnm
    type(sic_class) :: sic


    call nc_create(fnm)
    call nc_write_dim(fnm,"lon",x=lon,axis="x")
    call nc_write_dim(fnm,"lat",x=lat,axis="y")

    call nc_write(fnm,"u",     sic%u,     dims=["lon","lat"],long_name="zonal sea ice drift velocity",units="m/s")
    call nc_write(fnm,"v",     sic%v,     dims=["lon","lat"],long_name="meridional sea ice drift velocity",units="m/s")
    call nc_write(fnm,"h_sic_mean",   sic%h_sic_mean,   dims=["lon","lat"],long_name="grid cell mean sea ice thickness",units="m")
    call nc_write(fnm,"f_sic", sic%f_sic, dims=["lon","lat"],long_name="sea ice fraction",units="/")
    call nc_write(fnm,"h_sic", sic%h_sic, dims=["lon","lat"],long_name="sea ice thickness",units="m")
    call nc_write(fnm,"h_snow_mean",sic%h_snow_mean,dims=["lon","lat"],long_name="grid cell mean snow thickness",units="m")
    call nc_write(fnm,"h_snow",sic%h_snow,dims=["lon","lat"],long_name="snow thickness",units="m")
    call nc_write(fnm,"w_snow",sic%w_snow,dims=["lon","lat"],long_name="snow water equivalent",units="kg/m2")
    call nc_write(fnm,"w_snow_max",sic%w_snow_max,dims=["lon","lat"],long_name="seasonal maximum snow water equivalent",units="kg/m2")
    call nc_write(fnm,"str_s",sic%str_s,dims=["lon","lat"],long_name="The divergence stress tensor component ",units="Pa m")
    call nc_write(fnm,"str_d",sic%str_d,dims=["lon","lat"],long_name="The tension stress tensor component    ",units="Pa m")
    call nc_write(fnm,"str_t",sic%str_t,dims=["lon","lat"],long_name="The shearing stress tensor component   ",units="Pa m")
    call nc_write(fnm,"t_skin_sic", sic%t_skin_sic, dims=["lon","lat"],long_name="sea ice skin temperature",units="K")
    call nc_write(fnm,"t_skin_ocn", sic%t_skin_ocn, dims=["lon","lat"],long_name="ocean skin temperature",units="K")
    call nc_write(fnm,"sst", sic%sst, dims=["lon","lat"],long_name="top ocean temperature",units="K")
    call nc_write(fnm,"albedo_sic", sic%albedo_sic, dims=["lon","lat"],long_name="sea ice albedo",units="/")
    call nc_write(fnm,"alb_sic_vis_dir", sic%alb_sic_vis_dir, dims=["lon","lat"],long_name="sea ice albedo",units="/")
    call nc_write(fnm,"alb_sic_vis_dif", sic%alb_sic_vis_dif, dims=["lon","lat"],long_name="sea ice albedo",units="/")
    call nc_write(fnm,"alb_sic_nir_dir", sic%alb_sic_nir_dir, dims=["lon","lat"],long_name="sea ice albedo",units="/")
    call nc_write(fnm,"alb_sic_nir_dif", sic%alb_sic_nir_dif, dims=["lon","lat"],long_name="sea ice albedo",units="/")
    call nc_write(fnm,"alb_ocn", sic%albedo_ocn, dims=["lon","lat"],long_name="ocean albedo",units="/")
    call nc_write(fnm,"alb_ocn_vis_dir", sic%alb_ocn_vis_dir, dims=["lon","lat"],long_name="sea ice albedo",units="/")
    call nc_write(fnm,"alb_ocn_vis_dif", sic%alb_ocn_vis_dif, dims=["lon","lat"],long_name="sea ice albedo",units="/")
    call nc_write(fnm,"alb_ocn_nir_dir", sic%alb_ocn_nir_dir, dims=["lon","lat"],long_name="sea ice albedo",units="/")
    call nc_write(fnm,"alb_ocn_nir_dif", sic%alb_ocn_nir_dif, dims=["lon","lat"],long_name="sea ice albedo",units="/")
    call nc_write(fnm,"rough_m_ocn", sic%rough_m_ocn, dims=["lon","lat"],long_name="roughness length for momentum over ocean",units="m")
    call nc_write(fnm,"rough_m_sic", sic%rough_m_sic, dims=["lon","lat"],long_name="roughness length for momentum over sea ice",units="m")
    call nc_write(fnm,"rough_h_ocn", sic%rough_h_ocn, dims=["lon","lat"],long_name="roughness length for heat/moisture over ocean",units="m")
    call nc_write(fnm,"rough_h_sic", sic%rough_h_sic, dims=["lon","lat"],long_name="roughness length for heat/moisture over sea ice",units="m")
    call nc_write(fnm,"Cdh_ocn", sic%Cdh_ocn, dims=["lon","lat"],long_name="drag coefficient for heat/moisture over ocean",units="m")
    call nc_write(fnm,"Cdh_sic", sic%Cdh_sic, dims=["lon","lat"],long_name="drag coefficient for heat/moisture over sea ice",units="m")
    call nc_write(fnm,"flx_sh",  sic%flx_sh,  dims=["lon","lat"],long_name="sensible heat flux",units="W/m2")
    call nc_write(fnm,"flx_sh_ocn",  sic%flx_sh_ocn,  dims=["lon","lat"],long_name="sensible heat flux",units="W/m2")
    call nc_write(fnm,"flx_sh_sic",  sic%flx_sh_sic,  dims=["lon","lat"],long_name="sensible heat flux",units="W/m2")
    call nc_write(fnm,"flx_lh",  sic%flx_lh,  dims=["lon","lat"],long_name="latent heat flux",units="W/m2")
    call nc_write(fnm,"flx_lh_ocn",  sic%flx_lh_ocn,  dims=["lon","lat"],long_name="latent heat flux",units="W/m2")
    call nc_write(fnm,"flx_lh_sic",  sic%flx_lh_sic,  dims=["lon","lat"],long_name="latent heat flux",units="W/m2")
    call nc_write(fnm,"flx_lwu", sic%flx_lwu, dims=["lon","lat"],long_name="upward longwave at the surface",units="W/m2")
    call nc_write(fnm,"flx_lwu_ocn", sic%flx_lwu_ocn, dims=["lon","lat"],long_name="upward longwave at the surface",units="W/m2")
    call nc_write(fnm,"flx_lwu_sic", sic%flx_lwu_sic, dims=["lon","lat"],long_name="upward longwave at the surface",units="W/m2")
    call nc_write(fnm,"evp",     sic%evp,     dims=["lon","lat"],long_name="evaporation",units="kg/m2/s")
    call nc_write(fnm,"evp_ocn",     sic%evp_ocn,     dims=["lon","lat"],long_name="evaporation",units="kg/m2/s")
    call nc_write(fnm,"evp_sic",     sic%evp_sic,     dims=["lon","lat"],long_name="evaporation",units="kg/m2/s")


   return

  end subroutine sic_write_restart


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s i c _ r e a d _ r e s t a r t
  ! Purpose  :  read restart netcdf file 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sic_read_restart(fnm,sic)

    implicit none

    character (len=*) :: fnm
    type(sic_class) :: sic


    call nc_read(fnm,"u",     sic%u)
    call nc_read(fnm,"v",     sic%v)
    call nc_read(fnm,"h_sic_mean",   sic%h_sic_mean)
    call nc_read(fnm,"f_sic", sic%f_sic)
    call nc_read(fnm,"h_sic", sic%h_sic)
    call nc_read(fnm,"h_snow_mean",sic%h_snow_mean)
    call nc_read(fnm,"h_snow",sic%h_snow)
    call nc_read(fnm,"w_snow",sic%w_snow)
    call nc_read(fnm,"w_snow_max",sic%w_snow_max)
    call nc_read(fnm,"str_s",sic%str_s)
    call nc_read(fnm,"str_t",sic%str_t)
    call nc_read(fnm,"str_d",sic%str_d)
    call nc_read(fnm,"t_skin_sic", sic%t_skin_sic)
    call nc_read(fnm,"t_skin_ocn", sic%t_skin_ocn)
    call nc_read(fnm,"sst", sic%sst)
    call nc_read(fnm,"albedo_sic", sic%albedo_sic)
    call nc_read(fnm,"alb_sic_vis_dir", sic%alb_sic_vis_dir)
    call nc_read(fnm,"alb_sic_vis_dif", sic%alb_sic_vis_dif)
    call nc_read(fnm,"alb_sic_nir_dir", sic%alb_sic_nir_dir)
    call nc_read(fnm,"alb_sic_nir_dif", sic%alb_sic_nir_dif)
    call nc_read(fnm,"alb_ocn", sic%albedo_ocn)
    call nc_read(fnm,"alb_ocn_vis_dir", sic%alb_ocn_vis_dir)
    call nc_read(fnm,"alb_ocn_vis_dif", sic%alb_ocn_vis_dif)
    call nc_read(fnm,"alb_ocn_nir_dir", sic%alb_ocn_nir_dir)
    call nc_read(fnm,"alb_ocn_nir_dif", sic%alb_ocn_nir_dif)
    call nc_read(fnm,"rough_m_ocn", sic%rough_m_ocn)
    call nc_read(fnm,"rough_m_sic", sic%rough_m_sic)
    call nc_read(fnm,"rough_h_ocn", sic%rough_h_ocn)
    call nc_read(fnm,"rough_h_sic", sic%rough_h_sic)
    call nc_read(fnm,"Cdh_ocn", sic%Cdh_ocn)
    call nc_read(fnm,"Cdh_sic", sic%Cdh_sic)
    call nc_read(fnm,"flx_sh",      sic%flx_sh)
    call nc_read(fnm,"flx_sh_ocn",  sic%flx_sh_ocn)
    call nc_read(fnm,"flx_sh_sic",  sic%flx_sh_sic)
    call nc_read(fnm,"flx_lh",      sic%flx_lh)
    call nc_read(fnm,"flx_lh_ocn",  sic%flx_lh_ocn)
    call nc_read(fnm,"flx_lh_sic",  sic%flx_lh_sic)
    call nc_read(fnm,"flx_lwu",     sic%flx_lwu)
    call nc_read(fnm,"flx_lwu_ocn", sic%flx_lwu_ocn)
    call nc_read(fnm,"flx_lwu_sic", sic%flx_lwu_sic)
    call nc_read(fnm,"evp",         sic%evp)
    call nc_read(fnm,"evp_sic",     sic%evp_sic)
    call nc_read(fnm,"evp_ocn",     sic%evp_ocn)

    print *,'read restart file ',fnm

   return

  end subroutine sic_read_restart


end module sic_model

