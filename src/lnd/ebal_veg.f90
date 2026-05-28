!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : e b a l _ v e g _ m o d
!
!  Purpose : energy balance over vegetated and bare land
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
module ebal_veg_mod

   use precision, only : wp
   use timer, only : doy
   use constants, only : Le, Ls, Lf, q_sat_w, q_sat_i, dqsat_dt_w, dqsat_dt_i, T0, pi
   use constants, only : rho_a, cap_a, sigma
   use control, only : check_energy
   use lnd_grid, only : npft, nsurf, nveg, nl, z, flag_veg, i_bare, flag_pft
   use lnd_params, only : dt, rdt
   use lnd_params, only : surf_par, l_diurnal_cycle
   use lnd_params, only : pft_par
   use wiso_params, only : l_wiso, nwiso, i_o18, Rstd

   implicit none

   private
   public :: ebal_veg, update_tskin_veg

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  e b a l _ v e g
  !   Purpose    :  compute skin temperature and diagnose surface energy fluxes
  !              :  by solving the surface energy balance equation
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ebal_veg(frac_surf, mask_snow, h_snow, w_snow, lambda_veg, evap_canopy, subl_canopy, &
                     t_skin, t_skin_old, t_soil, tatm, qatm, pressure, swnet, swnet_min, lwdown, &
                     beta_s, r_s, beta_s_can, r_s_can, r_a, r_a_can, &
                     flx_g, dflxg_dT, flx_melt, flx_g_mean, dflxg_dT_mean, flx_melt_mean, t_skin_amp, &
                     num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw, &
                     f_sh, f_e, f_t, f_le, f_lt, f_lw, lh_ecan, qsat_e, dqsatdT_e, qsat_t, dqsatdT_t, &
                     energy_cons_surf1, i, j)

    implicit none

    integer, intent(in) :: mask_snow
    real(wp), intent(in) :: h_snow
    real(wp), intent(in) :: w_snow
    real(wp), dimension(:), intent(in) :: frac_surf, evap_canopy, subl_canopy
    real(wp), dimension(0:), intent(in) :: t_soil, lambda_veg
    real(wp), dimension(:), intent(in) :: tatm, qatm, lwdown, swnet, swnet_min
    real(wp), dimension(:), intent(in) :: pressure
    real(wp), dimension(:), intent(in) :: beta_s, r_s, r_a
    real(wp), dimension(:), intent(in) :: r_a_can, beta_s_can, r_s_can
    real(wp), dimension(:), intent(inout) :: t_skin, t_skin_old
    real(wp), dimension(:), intent(inout) :: flx_g, dflxg_dT, flx_melt
    real(wp), intent(out) :: flx_g_mean, dflxg_dT_mean, flx_melt_mean
    real(wp), dimension(:), intent(inout) :: t_skin_amp
    real(wp), dimension(:) :: num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw
    real(wp), dimension(:) :: f_sh, f_e, f_t, f_le, f_lt, f_lw, lh_ecan, qsat_e, dqsatdT_e, qsat_t, dqsatdT_t
    real(wp), dimension(:), intent(inout) :: energy_cons_surf1
 
    integer :: n
    real(wp) :: tsoil, p, rhoa
    real(wp) :: sh, lh, g, lw, sh_0, lh_0, g_0, lw_0
    real(wp) :: num, denom, num_g, denom_g
    real(wp) :: emiss, emiss_soil, f_veg
    real(wp) :: L_v_skin, L_v_sur
    real(wp) :: sh_pos, lh_pos, g_pos, lw_pos
    real(wp) :: t_skin_max, t_skin_min, t_skin_pos
    real(wp) :: t1, dt_pos, Ts, acos_fac, sqrt_fac


    integer :: i,j

    f_veg = sum(frac_surf,mask=flag_veg.eq.1)

    flx_g_mean    = 0._wp
    dflxg_dT_mean = 0._wp
    flx_melt      = 0._wp
    flx_melt_mean = 0._wp

    do n=1,nveg

      p = pressure(n)

      if (frac_surf(n).gt.0._wp) then

        if (n.eq.i_bare) then ! bare soil
        !---------------------------------------------------------
        ! bare soil
        !---------------------------------------------------------

          t_skin_old(n) = t_skin(n)
          if (mask_snow.eq.0) then
            tsoil = t_soil(1) ! first soil layer temperature
            emiss = surf_par%emissivity(n)
            L_v_skin = Le
            qsat_e(n) = q_sat_w(t_skin_old(n),p)
            dqsatdT_e(n) = dqsat_dT_w(t_skin_old(n),p)
          else
            tsoil = t_soil(0) ! snow temperature
            emiss = surf_par%emissivity_snow
            L_v_skin = Ls
            qsat_e(n) = q_sat_i(t_skin_old(n),p)
            dqsatdT_e(n) = dqsat_dT_i(t_skin_old(n),p)
          endif

          rhoa = rho_a(tatm(n),p)
          f_sh(n) = rhoa*cap_a/r_a(n)
          f_le(n) = L_v_skin * rhoa * beta_s(n) / (r_a(n) + r_s(n))
          f_e(n)  = f_le(n) / L_v_skin
          f_lw(n) = emiss*sigma

          lh_ecan(n)  = 0._wp

          ! ground heat flux
          if (mask_snow.eq.1) then
            num_g   = lambda_veg(0)/(0.5_wp*h_snow) * tsoil
            denom_g = lambda_veg(0)/(0.5_wp*h_snow)
          else
            num_g   = lambda_veg(1)/z(1) * tsoil
            denom_g = lambda_veg(1)/z(1)
          endif

          num_lh(n)   = - f_le(n) * (qsat_e(n) - dqsatdT_e(n)*t_skin_old(n) - qatm(n))  ! bare soil evaporation
          denom_lh(n) = f_le(n) * dqsatdT_e(n)

          num_sh(n)   = f_sh(n) * tatm(n)
          denom_sh(n) = f_sh(n)

          num_sw(n)  = swnet(n)

          num_lw(n)  = emiss*lwdown(n) + 3._wp*emiss*sigma*t_skin_old(n)**4
          denom_lw(n)= 4._wp * emiss*sigma*t_skin_old(n)**3

          num   = num_sw(n) + num_lw(n) + num_sh(n) + num_lh(n) + num_g

          denom = denom_lw(n) + denom_sh(n) + denom_lh(n) + denom_g

          ! new skin temperature        
          t_skin(n) = num/denom

          ! compute ground heat flux and derivative to be used as input to the soil temperature solver
          if (mask_snow.eq.1) then
            flx_g(n) = lambda_veg(0)/(0.5_wp*h_snow) * (t_skin(n) - tsoil)
            dflxg_dT(n) = - lambda_veg(0)/(0.5_wp*h_snow)
          else
            flx_g(n) = lambda_veg(1)/z(1) * (t_skin(n) - tsoil)
            dflxg_dT(n) = - lambda_veg(1)/z(1)
          endif

          ! diagnose fluxes
          sh = f_sh(n) * (t_skin(n) - tatm(n))
          lh = f_le(n) * (qsat_e(n) + dqsatdT_e(n)*(t_skin(n)-t_skin_old(n)) - qatm(n))
          lw = (1._wp-emiss)*lwdown(n) + emiss*sigma * (t_skin_old(n)**4 + 4._wp*t_skin_old(n)**3*(t_skin(n)-t_skin_old(n)))
          g = flx_g(n)

          if (.not.l_diurnal_cycle) then

            ! limit t_skin to <= 0°C if snow is on the ground
            if (h_snow.gt.0._wp .and. t_skin(n).gt.T0) then

              ! re-diagnose fluxes with t_skin = T0
              sh_0 = f_sh(n) * (T0 - tatm(n))
              lh_0 = f_le(n) * (qsat_e(n) + dqsatdT_e(n)*(T0-t_skin_old(n)) - qatm(n))
              lw_0 = (1._wp-emiss)*lwdown(n) + emiss*sigma * (t_skin_old(n)**4 + 4._wp*t_skin_old(n)**3*(T0-t_skin_old(n)))
              g_0 = lambda_veg(0)/(0.5_wp*max(h_snow,0.04_wp)) * (T0 - tsoil)

              ! diagnose remaining energy flux which can be used to melt snow
              flx_melt(n) = sh + lh + g + lw &
                - (sh_0+lh_0+g_0+lw_0)

              if (flx_melt(n).le.w_snow*Lf*rdt) then

                ! reset t_skin to melt temperature
                t_skin(n) = T0
                ! update fluxes
                lw = lw_0
                sh = sh_0
                lh = lh_0
                g  = g_0
                flx_g(n) = g

              else

                ! limit flx_melt to the amount of energy needed to melt the whole snow layer
                flx_melt(n) = w_snow*Lf*rdt

                ! rediagnose skin temperature substracting snowmelt energy
                t_skin(n)   = (num_sw(n) + num_lw(n) + num_sh(n) + num_lh(n) + num_g - flx_melt(n)) / denom

                ! diagnose fluxes with new skin temperature
                sh = f_sh(n) * (t_skin(n) - tatm(n))
                lh = f_lt(n) * (qsat_t(n) + dqsatdT_t(n)*(t_skin(n)-t_skin_old(n)) - qatm(n)) &
                  + f_le(n) * (qsat_e(n) + dqsatdT_e(n)*(t_skin(n)-t_skin_old(n)) - qatm(n)) &
                  + lh_ecan(n)
                lw = (1._wp-emiss)*lwdown(n) + emiss*sigma * (t_skin_old(n)**4 + 4._wp*t_skin_old(n)**3*(t_skin(n)-t_skin_old(n)))
                if(mask_snow.eq.1) then
                  flx_g(n) = lambda_veg(0)/(0.5_wp*max(h_snow,0.04_wp)) * (t_skin(n) - tsoil)
                else
                  flx_g(n) = lambda_veg(1)/z(1) * (t_skin(n) - tsoil)
                endif
                g = flx_g(n)

              endif

            endif

            t_skin_amp(n) = 0._wp

          else ! parameterisation of diurnal cycle of skin temperature

            ! diurnal maximum skin temperature derived by conserving daily shortwave
            ! radiation integral, assuming sinusoidal diurnal cycle of SW
            t_skin_max = (2._wp*swnet(n)-swnet_min(n) + num_lw(n) + num_sh(n) + num_lh(n) + num_g) / denom

            ! amplitude of diurnal cycle of skin temperature
            t_skin_amp(n) = max(0._wp, (t_skin_max-t_skin(n)))
            ! t_skin_max after scaling applied
            t_skin_max = t_skin(n)+t_skin_amp(n)
            ! diurnal minimum skin temperature (symmetric)
            t_skin_min = t_skin(n)-t_skin_amp(n)

            ! compute energy for melt from diurnal cycle, following (partly) Krapp et al. 2016 (SEMIC)
            if (h_snow.gt.0._wp) then

              if (t_skin_max.gt.T0) then
                ! re-diagnose fluxes with t_skin = T0
                sh_0 = f_sh(n) * (T0 - tatm(n))
                lh_0 = f_le(n) * (qsat_e(n) + dqsatdT_e(n)*(T0-t_skin_old(n)) - qatm(n))
                lw_0 = (1._wp-emiss)*lwdown(n) + emiss*sigma * (t_skin_old(n)**4 + 4._wp*t_skin_old(n)**3*(T0-t_skin_old(n)))
                if (mask_snow.eq.1) then
                  g_0 = lambda_veg(0)/(0.5_wp*h_snow) * (T0 - tsoil)
                else
                  g_0 = lambda_veg(1)/z(1) * (T0 - tsoil)
                endif
              endif

              if (t_skin_min.gt.T0) then ! no need to resolve diurnal cycle, just use daily mean

                ! diagnose remaining energy flux which can be used to melt snow/ice
                flx_melt(n) = sh + lh + g + lw &
                  - (sh_0+lh_0+g_0+lw_0)

                if (flx_melt(n).le.w_snow*Lf*rdt) then

                  ! reset t_skin to melt temperature
                  t_skin(n) = T0
                  ! update fluxes
                  lw = lw_0
                  sh = sh_0
                  lh = lh_0
                  g  = g_0
                  flx_g(n) = g

                else

                  ! limit flx_melt to the amount of energy needed to melt the whole snow layer
                  flx_melt(n) = w_snow*Lf*rdt

                  ! rediagnose skin temperature substracting snowmelt energy
                  t_skin(n)   = (num_sw(n) + num_lw(n) + num_sh(n) + num_lh(n) + num_g - flx_melt(n)) / denom

                  ! diagnose fluxes with new skin temperature
                  sh = f_sh(n) * (t_skin(n) - tatm(n))
                  lh = f_lt(n) * (qsat_t(n) + dqsatdT_t(n)*(t_skin(n)-t_skin_old(n)) - qatm(n)) &
                    + f_le(n) * (qsat_e(n) + dqsatdT_e(n)*(t_skin(n)-t_skin_old(n)) - qatm(n)) &
                    + lh_ecan(n)
                  lw = (1._wp-emiss)*lwdown(n) + emiss*sigma * (t_skin_old(n)**4 + 4._wp*t_skin_old(n)**3*(t_skin(n)-t_skin_old(n)))
                  if(mask_snow.eq.1) then
                    flx_g(n) = lambda_veg(0)/(0.5_wp*h_snow) * (t_skin(n) - tsoil)
                  else
                    flx_g(n) = lambda_veg(1)/z(1) * (t_skin(n) - tsoil)
                  endif
                  g = flx_g(n)

                endif

              else if (t_skin_max.gt.T0 .and. t_skin_amp(n)>0._wp) then

                Ts = t_skin(n)-T0
                acos_fac = acos(Ts/t_skin_amp(n))
                sqrt_fac = sqrt(1._wp-Ts**2/t_skin_amp(n)**2)
                t1 = dt/(2._wp*pi) * acos_fac

                ! compute flux used to melt snow
                if (t_skin_max.gt.T0) then   
                  ! compute mean skin temperature above T0
                  dt_pos = dt-2._wp*t1 ! s, time that temperature is above freezing point
                  t_skin_pos = T0 + dt/(pi*dt_pos)*(-Ts*acos_fac+t_skin_amp(n)*sqrt_fac+pi*Ts)
                  ! diagnose fluxes with t_skin = t_skin_pos
                  sh_pos = f_sh(n) * (t_skin_pos - tatm(n))
                  if (mask_snow.eq.1) then
                    g_pos  = lambda_veg(0)/(0.5_wp*h_snow) * (t_skin_pos - tsoil)
                  else
                    g_pos  = lambda_veg(1)/z(1) * (t_skin_pos - tsoil)
                  endif
                  lh_pos = f_le(n) * (qsat_e(n) + dqsatdT_e(n)*(t_skin_pos-t_skin_old(n)) - qatm(n))
                  lw_pos = (1._wp-emiss)*lwdown(n) + emiss*sigma * (t_skin_old(n)**4 + 4._wp*t_skin_old(n)**3*(t_skin_pos-t_skin_old(n)))
                  ! diagnose energy flux which can be used to melt snow/ice
                  flx_melt(n) = sh_pos + lh_pos + g_pos + lw_pos &
                    - (sh_0+lh_0+g_0+lw_0)
                  ! scale with dt_pos
                  flx_melt(n) = flx_melt(n)*dt_pos/dt
                endif

                flx_melt(n) = min(flx_melt(n),w_snow*Lf*rdt)

                ! rediagnose skin temperature substracting snowmelt energy
                t_skin(n)   = (num_sw(n) + num_lw(n) + num_sh(n) + num_lh(n) + num_g - flx_melt(n)) / denom

                ! diagnose fluxes with new skin temperature
                sh = f_sh(n) * (t_skin(n) - tatm(n))
                lh = f_le(n) * (qsat_e(n) + dqsatdT_e(n)*(t_skin(n)-t_skin_old(n)) - qatm(n))
                lw = (1._wp-emiss)*lwdown(n) + emiss*sigma * (t_skin_old(n)**4 + 4._wp*t_skin_old(n)**3*(t_skin(n)-t_skin_old(n)))
                if(mask_snow.eq.1) then
                  flx_g(n) = lambda_veg(0)/(0.5_wp*h_snow) * (t_skin(n) - tsoil)
                else
                  flx_g(n) = lambda_veg(1)/z(1) * (t_skin(n) - tsoil)
                endif
                g = flx_g(n)

              endif

            endif

          endif


        else
        !---------------------------------------------------------
        ! vegetation
        !---------------------------------------------------------

          t_skin_old(n) = t_skin(n)

          if (mask_snow.eq.0) then
            tsoil = t_soil(1) ! first soil layer temperature
            L_v_sur = Le
            qsat_e(n) = q_sat_w(t_skin_old(n),p)
            dqsatdT_e(n) = dqsat_dT_w(t_skin_old(n),p)
            emiss_soil = surf_par%emissivity(i_bare)
          else
            tsoil = t_soil(0) ! snow temperature
            L_v_sur = Ls
            qsat_e(n) = q_sat_i(t_skin_old(n),p)
            dqsatdT_e(n) = dqsat_dT_i(t_skin_old(n),p)
            emiss_soil = surf_par%emissivity_snow
          endif

          emiss = surf_par%emissivity(n)
          qsat_t(n) = q_sat_w(t_skin_old(n),p)
          dqsatdT_t(n) = dqsat_dT_w(t_skin_old(n),p)

          rhoa = rho_a(tatm(n),p)
          f_sh(n) = rhoa*cap_a/r_a(n)
          f_lt(n) = Le * rhoa * beta_s(n) / (r_a(n) + r_s(n))   ! transpiration latent heat factor
          f_le(n) = L_v_sur * rhoa * beta_s_can(n) / (r_a_can(n) + r_a(n) + r_s_can(n)) ! surface evaporation latent heat factor
          f_t(n)  = f_lt(n) / Le    ! transpiration factor
          f_e(n)  = f_le(n) / L_v_sur   ! surface evaporation factor
          f_lw(n) = emiss*sigma

          ! ground heat flux
          if (mask_snow.eq.1) then
            num_g   = lambda_veg(0)/(0.5_wp*h_snow) * tsoil
            denom_g = lambda_veg(0)/(0.5_wp*h_snow)
          else
            num_g   = lambda_veg(1)/z(1) * tsoil
            denom_g = lambda_veg(1)/z(1)
          endif

          lh_ecan(n)  = Le * evap_canopy(n) + Ls * subl_canopy(n) ! latent heat flux of wet canopy evaporation/sublimation
          num_lh(n)   = - f_lt(n) * (qsat_t(n) - dqsatdT_t(n)*t_skin_old(n) - qatm(n)) &  ! transpiration
            - f_le(n) * (qsat_e(n) - dqsatdT_e(n)*t_skin_old(n) - qatm(n)) &  ! evaporation from soil below the canopy
            - lh_ecan(n)  ! wet canopy evaporation
          denom_lh(n) = f_lt(n) * dqsatdT_t(n) + f_le(n) * dqsatdT_e(n) 

          num_sh(n)   = f_sh(n) * tatm(n)
          denom_sh(n) = f_sh(n)

          num_sw(n)  = swnet(n)

          num_lw(n)  = emiss*lwdown(n) + 3._wp*emiss*sigma*t_skin_old(n)**4
          denom_lw(n)= 4._wp * emiss*sigma*t_skin_old(n)**3

          num   = num_sw(n) + num_lw(n) + num_sh(n) + num_lh(n) + num_g

          denom = denom_lw(n) + denom_sh(n) + denom_lh(n) + denom_g

          ! new skin temperature        
          t_skin(n) = num/denom

          ! compute ground heat flux and derivative to be used as input to the soil temperature solver
          if (mask_snow.eq.1) then
            flx_g(n) = lambda_veg(0)/(0.5_wp*h_snow) * (t_skin(n) - tsoil)
            dflxg_dT(n) = - lambda_veg(0)/(0.5_wp*h_snow)
          else
            flx_g(n) = lambda_veg(1)/z(1) * (t_skin(n) - tsoil)
            dflxg_dT(n) = - lambda_veg(1)/z(1)
          endif

          ! diagnose fluxes
          sh = f_sh(n) * (t_skin(n) - tatm(n))
          lh = f_lt(n) * (qsat_t(n) + dqsatdT_t(n)*(t_skin(n)-t_skin_old(n)) - qatm(n)) &
            + f_le(n) * (qsat_e(n) + dqsatdT_e(n)*(t_skin(n)-t_skin_old(n)) - qatm(n)) &
            + lh_ecan(n)
          lw = (1._wp-emiss)*lwdown(n) + emiss*sigma * (t_skin_old(n)**4 + 4._wp*t_skin_old(n)**3*(t_skin(n)-t_skin_old(n)))
          g  = flx_g(n)

          if (.not.l_diurnal_cycle) then

            ! limit t_skin to <= 0°C if snow is on the ground
            if (mask_snow.eq.1 .and. t_skin(n).gt.T0) then

              ! re-diagnose fluxes with t_skin = T0
              sh_0 = f_sh(n) * (T0 - tatm(n))
              lh_0 = f_lt(n) * (qsat_t(n) + dqsatdT_t(n)*(T0-t_skin_old(n)) - qatm(n)) &
                + f_le(n) * (qsat_e(n) + dqsatdT_e(n)*(T0-t_skin_old(n)) - qatm(n)) &
                + lh_ecan(n)
              lw_0 = (1._wp-emiss)*lwdown(n) + emiss*sigma * (t_skin_old(n)**4 + 4._wp*t_skin_old(n)**3*(T0-t_skin_old(n)))
              g_0 = lambda_veg(0)/(0.5_wp*h_snow) * (T0 - tsoil)

              ! diagnose remaining energy flux which can be used to melt snow
              flx_melt(n) = sh + lh + g + lw &
                - (sh_0+lh_0+g_0+lw_0)

              if (flx_melt(n).le.w_snow*Lf*rdt) then

                ! reset t_skin to melt temperature
                t_skin(n) = T0
                ! update fluxes
                lw = lw_0
                sh = sh_0
                lh = lh_0
                g  = g_0
                flx_g(n) = g

              else

                ! limit flx_melt to the amount of energy needed to melt the whole snow layer
                flx_melt(n) = w_snow*Lf*rdt

                ! rediagnose skin temperature substracting snowmelt energy
                t_skin(n)   = (num_sw(n) + num_lw(n) + num_sh(n) + num_lh(n) + num_g - flx_melt(n)) / denom

                ! diagnose fluxes with new skin temperature
                sh = f_sh(n) * (t_skin(n) - tatm(n))
                lh = f_lt(n) * (qsat_t(n) + dqsatdT_t(n)*(t_skin(n)-t_skin_old(n)) - qatm(n)) &
                  + f_le(n) * (qsat_e(n) + dqsatdT_e(n)*(t_skin(n)-t_skin_old(n)) - qatm(n)) &
                  + lh_ecan(n)
                lw = (1._wp-emiss)*lwdown(n) + emiss*sigma * (t_skin_old(n)**4 + 4._wp*t_skin_old(n)**3*(t_skin(n)-t_skin_old(n)))
                if(mask_snow.eq.1) then
                  flx_g(n) = lambda_veg(0)/(0.5_wp*h_snow) * (t_skin(n) - tsoil)
                else
                  flx_g(n) = lambda_veg(1)/z(1) * (t_skin(n) - tsoil)
                endif
                g = flx_g(n)

              endif

            endif

            t_skin_amp(n) = 0._wp 

          else ! parameterisation of diurnal cycle of skin temperature

            ! diurnal maximum skin temperature derived by conserving daily shortwave
            ! radiation integral, assuming sinusoidal diurnal cycle of SW
            t_skin_max = (2._wp*swnet(n)-swnet_min(n) + num_lw(n) + num_sh(n) + num_lh(n) + num_g) / denom

            ! amplitude of diurnal cycle of skin temperature
            t_skin_amp(n) = max(0._wp, (t_skin_max-t_skin(n)))
            ! t_skin_max after scaling applied
            t_skin_max = t_skin(n)+t_skin_amp(n)
            ! diurnal minimum skin temperature (symmetric)
            t_skin_min = t_skin(n)-t_skin_amp(n)

            ! compute energy for melt from diurnal cycle, following (partly) Krapp et al. 2016 (SEMIC)
            if (h_snow.gt.0._wp) then

              if (t_skin_max.gt.T0) then
                ! re-diagnose fluxes with t_skin = T0
                sh_0 = f_sh(n) * (T0 - tatm(n))
                lh_0 = f_lt(n) * (qsat_t(n) + dqsatdT_t(n)*(T0-t_skin_old(n)) - qatm(n)) &
                  + f_le(n) * (qsat_e(n) + dqsatdT_e(n)*(T0-t_skin_old(n)) - qatm(n)) &
                  + lh_ecan(n)
                lw_0 = (1._wp-emiss)*lwdown(n) + emiss*sigma * (t_skin_old(n)**4 + 4._wp*t_skin_old(n)**3*(T0-t_skin_old(n)))
                if (mask_snow.eq.1) then
                  g_0 = lambda_veg(0)/(0.5_wp*h_snow) * (T0 - tsoil)
                else
                  g_0 = lambda_veg(1)/z(1) * (T0 - tsoil)
                endif
              endif

              if (t_skin_min.gt.T0) then ! no need to resolve diurnal cycle, just use daily mean

                ! diagnose remaining energy flux which can be used to melt snow/ice
                flx_melt(n) = sh + lh + g + lw &
                  - (sh_0+lh_0+g_0+lw_0)

                if (flx_melt(n).le.w_snow*Lf*rdt) then

                  ! reset t_skin to melt temperature
                  t_skin(n) = T0
                  ! update fluxes
                  lw = lw_0
                  sh = sh_0
                  lh = lh_0
                  g  = g_0
                  flx_g(n) = g

                else

                  ! limit flx_melt to the amount of energy needed to melt the whole snow layer
                  flx_melt(n) = w_snow*Lf*rdt

                  ! rediagnose skin temperature substracting snowmelt energy
                  t_skin(n)   = (num_sw(n) + num_lw(n) + num_sh(n) + num_lh(n) + num_g - flx_melt(n)) / denom

                  ! diagnose fluxes with new skin temperature
                  sh = f_sh(n) * (t_skin(n) - tatm(n))
                  lh = f_lt(n) * (qsat_t(n) + dqsatdT_t(n)*(t_skin(n)-t_skin_old(n)) - qatm(n)) &
                    + f_le(n) * (qsat_e(n) + dqsatdT_e(n)*(t_skin(n)-t_skin_old(n)) - qatm(n)) &
                    + lh_ecan(n)
                  lw = (1._wp-emiss)*lwdown(n) + emiss*sigma * (t_skin_old(n)**4 + 4._wp*t_skin_old(n)**3*(t_skin(n)-t_skin_old(n)))
                  if(mask_snow.eq.1) then
                    flx_g(n) = lambda_veg(0)/(0.5_wp*h_snow) * (t_skin(n) - tsoil)
                  else
                    flx_g(n) = lambda_veg(1)/z(1) * (t_skin(n) - tsoil)
                  endif
                  g = flx_g(n)

                endif

              else if (t_skin_max.gt.T0 .and. t_skin_amp(n)>0._wp) then

                Ts = t_skin(n)-T0
                acos_fac = acos(Ts/t_skin_amp(n))
                sqrt_fac = sqrt(1._wp-Ts**2/t_skin_amp(n)**2)
                t1 = dt/(2._wp*pi) * acos_fac

                ! compute flux used to melt snow
                if (t_skin_max.gt.T0) then   
                  ! compute mean skin temperature above T0
                  dt_pos = dt-2._wp*t1 ! s, time that temperature is above freezing point
                  t_skin_pos = T0 + dt/(pi*dt_pos)*(-Ts*acos_fac+t_skin_amp(n)*sqrt_fac+pi*Ts)
                  ! diagnose fluxes with t_skin = t_skin_pos
                  sh_pos = f_sh(n) * (t_skin_pos - tatm(n))
                  if (mask_snow.eq.1) then
                    g_pos  = lambda_veg(0)/(0.5_wp*h_snow) * (t_skin_pos - tsoil)
                  else
                    g_pos  = lambda_veg(1)/z(1) * (t_skin_pos - tsoil)
                  endif
                  lh_pos = f_lt(n) * (qsat_t(n) + dqsatdT_t(n)*(t_skin_pos-t_skin_old(n)) - qatm(n)) &
                    + f_le(n) * (qsat_e(n) + dqsatdT_e(n)*(t_skin_pos-t_skin_old(n)) - qatm(n)) &
                    + lh_ecan(n)
                  lw_pos = (1._wp-emiss)*lwdown(n) + emiss*sigma * (t_skin_old(n)**4 + 4._wp*t_skin_old(n)**3*(t_skin_pos-t_skin_old(n)))
                  ! diagnose energy flux which can be used to melt snow/ice
                  flx_melt(n) = sh_pos + lh_pos + g_pos + lw_pos &
                    - (sh_0+lh_0+g_0+lw_0)
                  ! scale with dt_pos
                  flx_melt(n) = flx_melt(n)*dt_pos/dt
                endif

                flx_melt(n) = min(flx_melt(n),w_snow*Lf*rdt)

                ! rediagnose skin temperature substracting snowmelt energy
                t_skin(n)   = (num_sw(n) + num_lw(n) + num_sh(n) + num_lh(n) + num_g - flx_melt(n)) / denom

                ! diagnose fluxes with new skin temperature
                sh = f_sh(n) * (t_skin(n) - tatm(n))
                lh = f_lt(n) * (qsat_t(n) + dqsatdT_t(n)*(t_skin(n)-t_skin_old(n)) - qatm(n)) &
                  + f_le(n) * (qsat_e(n) + dqsatdT_e(n)*(t_skin(n)-t_skin_old(n)) - qatm(n)) &
                  + lh_ecan(n)
                lw = (1._wp-emiss)*lwdown(n) + emiss*sigma * (t_skin_old(n)**4 + 4._wp*t_skin_old(n)**3*(t_skin(n)-t_skin_old(n)))
                if(mask_snow.eq.1) then
                  flx_g(n) = lambda_veg(0)/(0.5_wp*h_snow) * (t_skin(n) - tsoil)
                else
                  flx_g(n) = lambda_veg(1)/z(1) * (t_skin(n) - tsoil)
                endif
                g = flx_g(n)

              endif

            endif

          endif

        endif


        if (check_energy) then
          ! energy conservation check
          energy_cons_surf1(n) = swnet(n) + lwdown(n) - lw - sh - lh - g - flx_melt(n)
          if (abs(energy_cons_surf1(n)).gt.1.d-9) then
            print *,''
            print *,'energy balance',n,energy_cons_surf1(n),mask_snow
            print *,'swnet,lwd,lwu,sh,lh,g,melt',swnet(n),lwdown(n),lw,sh,lh,g,flx_melt(n)
            print *,'h_snow',h_snow
            print *,'t_skin',t_skin(n)
            print *,'t_skin_old',t_skin_old(n)
            print *,'t_skin_amp',t_skin_amp(n)
            print *,'t_skin_min',t_skin_min
            print *,'t_skin_max',t_skin_max
            print *,'t_soil',tsoil
          endif
        endif

        if (t_skin(n).gt.350._wp .or. t_skin(n).lt.150._wp) then
        !if (i.eq.19.and.j.eq.31) then
          print *,''
          print *,'T_SKIN out of range!!!',t_skin(n)
          print *,'i,j,n',i,j,n
          print *,'doy',doy
          print *,'t_skin, t_skin_old',t_skin(n),t_skin_old(n)
          print *,'t_atm',tatm(n)
          print *,'t_soil',t_soil
          print *,'swnet',swnet(n)
          print *,'lwdown',lwdown(n)
          print *,'lw',lw
          print *,'sh',sh
          print *,'lh',lh
          print *,'g',g
          print *,'melt',flx_melt(n)
          print *,'mask_snow',mask_snow
          print *,'h_snow',h_snow
          print *,'lambda_veg',lambda_veg
          print *,'frac_s',frac_surf
          print *,'evap_canopy',evap_canopy(n)
          print *,'beta_s',beta_s(n)
          print *,'r_s',r_s(n)
          print *,'r_a',r_a(n)
          if (flag_pft(n).eq.1) then
            print *,'r_a_can',r_a_can(n)
            print *,'r_s_can',r_s_can(n)
            print *,'beta_s_can',beta_s_can(n)
          endif
          print *,''
        endif

        ! compute mean ground heat flux and derivative over vegetated area
        flx_g_mean = flx_g_mean + flx_g(n) * frac_surf(n)/f_veg
        dflxg_dT_mean = dflxg_dT_mean + dflxg_dT(n) * frac_surf(n)/f_veg
        ! compute mean energy flux going to melt snow
        flx_melt_mean = flx_melt_mean + flx_melt(n) * frac_surf(n)/f_veg

      endif
    enddo

    if (maxval(t_skin(1:npft+1)).gt.350._wp .or. minval(t_skin(1:npft+1)).lt.150._wp) then
      print *,''
      print *,'t_skin out of range!!!',t_skin(1:npft+1)
    endif

    return

  end subroutine ebal_veg


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  u p d a t e _ t s k i n _ v e g
  !   Purpose    :  update skin temperature using new ground heat flux computed from new topsoil temperature 
  !              :  and re-diagnose surface energy fluxes
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine update_tskin_veg(frac_surf, mask_snow, t_skin_old, dflxg_dT, tatm, qatm, swnet, lwdown, &
                         t_soil, t_soil_old, evap_canopy, subl_canopy, &
                         flx_g, flx_melt, t_skin, t_skin_veg, &
                         flx_sh, flx_lwu, lwnet, flx_lh, evap_surface, transpiration, et, &
                         num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw, &
                         f_sh, f_e, f_t, f_le, f_lt, f_lw, lh_ecan, qsat_e, dqsatdT_e, qsat_t, dqsatdT_t, &
                         energy_cons_surf2, &
                         w_w,w_w_iso,wilt,evap_canopy_iso,subl_canopy_iso, &
                         evap_surface_iso,transpiration_iso,et_iso)

    implicit none

    integer, intent(in) :: mask_snow
    real(wp), dimension(:), intent(in) :: frac_surf, t_skin_old, dflxg_dT
    real(wp), dimension(:), intent(in) :: tatm, qatm, lwdown, swnet
    real(wp), dimension(0:), intent(in) :: t_soil, t_soil_old
    real(wp), dimension(:), intent(in) :: evap_canopy, subl_canopy
    real(wp), dimension(:), intent(inout) :: flx_g, flx_melt, t_skin
    real(wp), dimension(:), intent(inout) :: flx_sh, flx_lwu, lwnet
    real(wp), dimension(:), intent(inout) :: flx_lh, evap_surface, transpiration, et
    real(wp), dimension(:), intent(inout) ::   num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw
    real(wp), dimension(:), intent(inout) :: f_sh, f_e, f_t, f_le, f_lt, f_lw, lh_ecan, qsat_e, dqsatdT_e, qsat_t, dqsatdT_t
    real(wp), dimension(:), intent(inout) :: energy_cons_surf2
    real(wp), intent(out) :: t_skin_veg
    ! water-isotope siblings (always passed; values are 0 unless l_wiso=.true.)
    real(wp), dimension(:),   intent(in)    :: w_w
    real(wp), dimension(:,:), intent(in)    :: w_w_iso, wilt
    real(wp), dimension(:,:), intent(in)    :: evap_canopy_iso, subl_canopy_iso
    real(wp), dimension(:,:), intent(inout) :: evap_surface_iso, transpiration_iso, et_iso

    integer :: n, k, iso
    real(wp) :: num, denom, f_veg, wilt_tot, r_top, r_xylem(nwiso)


    t_skin_veg = 0._wp

    ! initialize
    et(:) = 0._wp
    evap_surface(:) = 0._wp
    transpiration(:) = 0._wp
    if (l_wiso) then
      et_iso(:,:)            = 0._wp
      evap_surface_iso(:,:)  = 0._wp
      transpiration_iso(:,:) = 0._wp
    endif

    f_veg = sum(frac_surf,mask=flag_veg.eq.1)

    ! top-soil iso ratio for surface evaporation (donor = layer 1)
    if (l_wiso .and. w_w(1).gt.0._wp) then
      ! r_top is just a scalar at the top-layer ratio for each iso species; used inline below
    endif

    do n=1,nveg

      if( frac_surf(n) .gt. 0._wp ) then

        ! update soil heat flux
        if( mask_snow .eq. 1 ) then
          flx_g(n)  = flx_g(n) + dflxg_dT(n) * (t_soil(0) - t_soil_old(0))
        else
          flx_g(n)  = flx_g(n) + dflxg_dT(n) * (t_soil(1) - t_soil_old(1))
        endif

        num   = num_sw(n) + num_lw(n) + num_sh(n) + num_lh(n)  &
          - flx_g(n) - flx_melt(n)
        denom = denom_lw(n) + denom_sh(n) + denom_lh(n)

        ! new skin temperature
        t_skin(n) = num/denom

        ! diagnose surface energy fluxes from skin temperature
        flx_sh(n)  = f_sh(n) * (t_skin(n) - tatm(n))
        flx_lwu(n) = (1._wp-f_lw(n)/sigma)*lwdown(n) + f_lw(n) * (t_skin_old(n)**4 + 4._wp*t_skin_old(n)**3*(t_skin(n)-t_skin_old(n)))
        lwnet(n)   = lwdown(n) - flx_lwu(n)
        flx_lh(n)  = f_lt(n) * (qsat_t(n) + dqsatdT_t(n)*(t_skin(n)-t_skin_old(n)) - qatm(n)) &
          + f_le(n) * (qsat_e(n) + dqsatdT_e(n)*(t_skin(n)-t_skin_old(n)) - qatm(n)) &
          + lh_ecan(n)
        if(flag_pft(n) .eq. 1) then
          ! evaporation from surface (soil/snow)
          evap_surface(n)  = f_e(n) * (qsat_e(n) + dqsatdT_e(n)*(t_skin(n)-t_skin_old(n)) - qatm(n))
          ! transpiration
          transpiration(n) = f_t(n) * (qsat_t(n) + dqsatdT_t(n)*(t_skin(n)-t_skin_old(n)) - qatm(n))
        else
          ! evaporation from surface (soil/snow)
          evap_surface(n)  = f_e(n) * (qsat_e(n) + dqsatdT_e(n)*(t_skin(n)-t_skin_old(n)) - qatm(n))
          transpiration(n) = 0._wp
        endif
        ! total evapotranspiration
        et(n) = transpiration(n) + evap_surface(n) + evap_canopy(n) + subl_canopy(n)

        ! ---------- water-isotope ET assembly (passive, donor-pool ratio) ----------
        if (l_wiso) then
          do iso=1,nwiso
            ! surface evaporation: at top-soil ratio if outflow, else at VSMOW (dew from atm placeholder)
            if (evap_surface(n).ge.0._wp .and. w_w(1).gt.0._wp) then
              evap_surface_iso(n,iso) = evap_surface(n) * w_w_iso(1,iso)/w_w(1)
            else
              evap_surface_iso(n,iso) = evap_surface(n) * Rstd(iso)
            endif
            ! transpiration: outflow at root-weighted soil ratio, inflow (condensation) at Rstd
            if (flag_pft(n) .eq. 1) then
              if (transpiration(n) .ge. 0._wp) then
                wilt_tot = sum(wilt(:,n) * pft_par%root_frac(:,n))
                r_xylem(iso) = 0._wp
                if (wilt_tot.gt.1.e-20_wp) then
                  do k=1,nl
                    if (w_w(k).gt.0._wp) then
                      r_xylem(iso) = r_xylem(iso) &
                                    + pft_par%root_frac(k,n) * wilt(k,n) * w_w_iso(k,iso)/w_w(k) / wilt_tot
                    endif
                  enddo
                endif
                transpiration_iso(n,iso) = transpiration(n) * r_xylem(iso)
              else
                ! negative transpiration = vapor flux from atm onto canopy/soil at standard ratio
                transpiration_iso(n,iso) = transpiration(n) * Rstd(iso)
              endif
            else
              transpiration_iso(n,iso) = 0._wp
            endif
            et_iso(n,iso) = transpiration_iso(n,iso) + evap_surface_iso(n,iso) &
                          + evap_canopy_iso(n,iso) + subl_canopy_iso(n,iso)
          enddo
        endif

        if( check_energy ) then
          ! energy conservation check
          energy_cons_surf2(n) = swnet(n) &
          + lwdown(n)        &
          - flx_lwu(n) &
          - flx_sh(n) &
          - flx_lh(n) &
          - flx_g(n)  &
          - flx_melt(n)
          if( abs(energy_cons_surf2(n)) .gt. 1.e-9_wp ) then
            print *,''
            print *,'energy balance',n,energy_cons_surf2(n),mask_snow
          endif

        endif

        ! compute mean skin temperature over vegetation
        t_skin_veg = t_skin_veg + t_skin(n) * frac_surf(n)/f_veg

      endif
    enddo

    ! assign mean skin temperature to non-existing surface types
    do n=1,nveg
      if (frac_surf(n).eq.0._wp) then
        t_skin(n) = t_skin_veg
      endif
    enddo


    if( maxval(t_skin).gt.350. .or. minval(t_skin).lt.150. ) then
      print *,''
      print *,'t_skin out of range update_tskin!!!',t_skin
      print *,'t_skin_old',t_skin_old
      print *,'t_soil',t_soil
      !stop
    endif

    return

  end subroutine update_tskin_veg

end module ebal_veg_mod

