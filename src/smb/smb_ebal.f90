!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s m b _ e b a l
!
!  Purpose : Surface energy balance for SEMIX
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
module smb_ebal_mod

   use precision, only : wp
   use constants, only : Ls, q_sat_i, dqsat_dt_i
   use constants, only : rho_a, cap_a, sigma, lambda_i
   use constants, only : T0, pi 
   use control, only : check_energy
   use smb_grid, only : z, nl
   use smb_params, only : surf_par, snow_par, l_diurnal_cycle, tstd_scale, dt, l_dew

   implicit none

   private
   public :: ebal, update_tskin

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  e b a l 
  !   Purpose    :  compute skin temperature and diagnose surface energy fluxes
  !              :  by solving the surface energy balance equation
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ebal(mask_snow, h_snow, &
                     t_prof, t2m, tstd, q2m, pressure, swnet, swnet_min, lwdown, r_a, &
                     t_skin, &
                     t_skin_old, t_skin_amp, flx_g, dflxg_dT, flx_melt, &
                     num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw, &
                     f_sh, f_e, f_lh, f_lw, qsat, dqsatdT, i,j)

    implicit none

    integer, intent(in) :: mask_snow
    real(wp), intent(in) :: h_snow
    real(wp), dimension(0:), intent(in) :: t_prof
    real(wp), intent(in) :: t2m, tstd, q2m, swnet, swnet_min
    real(wp), intent(in) :: pressure, lwdown
    real(wp), intent(in) :: r_a

    real(wp), intent(inout) :: t_skin

    real(wp), intent(out) :: t_skin_old, t_skin_amp
    real(wp), intent(out) :: flx_g, dflxg_dT, flx_melt
    real(wp), intent(out) :: num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw
    real(wp), intent(out) :: f_sh, f_e, f_lh, f_lw, qsat, dqsatdT

    real(wp) :: energy_cons_surf1

    real(wp) :: tsoil, rhoa
    real(wp) :: sh, lh, g, lw, sh_0, lh_0, g_0, lw_0 
    real(wp) :: sh_pos, lh_pos, g_pos, lw_pos
    real(wp) :: num, denom, num_g, denom_g, emiss
    real(wp) :: t_skin_max, t_skin_min, t_skin_pos, t_skin_neg
    real(wp) :: t1, dt_pos, dt_neg, Ts, acos_fac, sqrt_fac
    real(wp) :: flx_melt_0, flx_melt_diurnal

    integer :: i,j


    flx_melt = 0._wp
    t_skin_amp = 0._wp

    t_skin_old = t_skin ! old skin temperature

    if(mask_snow.eq.1) then
      ! snow
      tsoil = t_prof(0) ! snow temperature
      emiss = surf_par%emissivity_snow
    else
      tsoil = t_prof(1) ! top layer temperature
      emiss = surf_par%emissivity_ice
    endif
    qsat = q_sat_i(t_skin,pressure)
    dqsatdT = dqsat_dT_i(t_skin,pressure)

    rhoa = rho_a(t2m,pressure)
    f_sh = rhoa*cap_a/r_a
    if (.not.l_dew .and. q2m.gt.qsat) then
      ! inhibit dew/frost deposition (negative latent heat flux)
      f_lh = 0._wp
      f_e = 0._wp
    else
      f_lh = Ls / r_a * rhoa
      f_e  = 1._wp /(r_a * rhoa)
    endif
    f_lw = emiss*sigma

    ! ground heat flux
    if(mask_snow.eq.1) then
      num_g   = snow_par%lambda/(0.5_wp*h_snow) * tsoil
      denom_g = snow_par%lambda/(0.5_wp*h_snow)
    else
      num_g   = lambda_i/z(1) * tsoil
      denom_g = lambda_i/z(1)
    endif

    num_lh   = - f_lh * (qsat - dqsatdT*t_skin - q2m)  ! surface sublimation
    denom_lh = f_lh * dqsatdT 

    num_sh   = f_sh * t2m
    denom_sh = f_sh

    num_sw   = swnet

    num_lw   = emiss*lwdown + 3._wp*emiss*sigma*t_skin**4
    denom_lw = 4._wp * emiss*sigma*t_skin**3

    num   = num_sw + num_lw + num_sh + num_lh + num_g
    denom = denom_lw + denom_sh + denom_lh + denom_g

    ! new skin temperature        
    t_skin   = num/denom

    ! compute ground heat flux and derivative to be used as input to the temperature profile solver
    if(mask_snow.eq.1) then
      flx_g = snow_par%lambda/(0.5_wp*h_snow) * (t_skin - tsoil)
      dflxg_dT = - snow_par%lambda/(0.5_wp*h_snow)
    else
      flx_g = lambda_i/z(1) * (t_skin - tsoil)
      dflxg_dT = - lambda_i/z(1)
    endif

    ! diagnose fluxes 
    sh = f_sh * (t_skin - t2m)
    lh = f_lh * (qsat + dqsatdT*(t_skin-t_skin_old) - q2m)
    lw = (1._wp-emiss)*lwdown + emiss*sigma * (t_skin_old**4 + 4._wp*t_skin_old**3*(t_skin-t_skin_old))
    g  = flx_g

    if (l_diurnal_cycle) then   
      ! parameterisation of diurnal cycle of skin temperature and melt

      ! diurnal maximum skin temperature derived by conserving daily shortwave
      ! radiation integral, assuming sinusoidal diurnal cycle of SW
      t_skin_max = (2._wp*swnet-swnet_min + num_lw + num_sh + num_lh + num_g) / denom

      ! amplitude of diurnal cycle of skin temperature
      ! add standard deviation of air/skin temperature to diurnal cycle amplitude to represent weather
      t_skin_amp = max(0._wp, (t_skin_max-t_skin) + tstd*tstd_scale)
      ! t_skin_max after scaling applied
      t_skin_max = t_skin+t_skin_amp
      ! diurnal minimum skin temperature (symmetric)
      t_skin_min = t_skin-t_skin_amp

      ! diagnose fluxes with t_skin = T0
      sh_0 = f_sh * (T0 - t2m)
      if (mask_snow.eq.1) then
        g_0  = snow_par%lambda/(0.5_wp*h_snow) * (T0 - tsoil)
      else
        g_0  = lambda_i/z(1) * (T0 - tsoil)
      endif
      lh_0 = f_lh * (qsat + dqsatdT*(T0-t_skin_old) - q2m)
      lw_0 = (1._wp-emiss)*lwdown + emiss*sigma * (t_skin_old**4 + 4._wp*t_skin_old**3*(T0-t_skin_old))

      ! diagnose melt flux without diurnal cycle parameterisation
      flx_melt_0 = sh + lh + g + lw &
        - (sh_0+lh_0+g_0+lw_0)
      flx_melt_0 = max(0._wp,flx_melt_0)

      if (t_skin_min.gt.T0) then
        ! skin temperature above freezing all day, no need to resolve diurnal cycle 

        ! melt flux is equal to melt flux without diurnal cycle parameterisation already computed above
        flx_melt = flx_melt_0 

        ! reset t_skin to freezing point 
        t_skin   = T0
        ! update fluxes
        lw = lw_0
        sh = sh_0
        lh = lh_0
        g  = g_0
        flx_g = g

      else if (t_skin_max.gt.T0 .and. t_skin_amp>0._wp) then
        ! daily maximum skin temperature is above melting point, compute diurnal
        ! cycle and diagnose flux available for melt

        ! compute diurnal cycle of temperature following Krapp et al. 2016 (SEMIC)
        Ts = t_skin-T0
        acos_fac = acos(Ts/t_skin_amp)
        sqrt_fac = sqrt(1._wp-Ts**2/t_skin_amp**2)
        t1 = dt/(2._wp*pi) * acos_fac
        dt_pos = dt-2._wp*t1 ! s, time that temperature is above freezing point
        ! compute mean skin temperature above T0
        t_skin_pos = T0 + dt/(pi*dt_pos)*(-Ts*acos_fac+t_skin_amp*sqrt_fac+pi*Ts)
        ! diagnose fluxes with t_skin = t_skin_pos
        sh_pos = f_sh * (t_skin_pos - t2m)
        if (mask_snow.eq.1) then
          g_pos  = snow_par%lambda/(0.5_wp*h_snow) * (t_skin_pos - tsoil)
        else
          g_pos  = lambda_i/z(1) * (t_skin_pos - tsoil)
        endif
        lh_pos = f_lh * (qsat + dqsatdT*(t_skin_pos-t_skin_old) - q2m)
        lw_pos = (1._wp-emiss)*lwdown + emiss*sigma * (t_skin_old**4 + 4._wp*t_skin_old**3*(t_skin_pos-t_skin_old))
        ! compute energy flux which can be used to melt snow/ice
        flx_melt = sh_pos + lh_pos + g_pos + lw_pos &
          - (sh_0+lh_0+g_0+lw_0)
        ! scale with dt_pos
        flx_melt = flx_melt*dt_pos/dt

        ! reset skin temperature to T0 if above
        if (t_skin>T0) then
          t_skin = T0
          ! update fluxes
          lw = lw_0
          sh = sh_0
          lh = lh_0
          g  = g_0
          flx_g = g
        endif

        ! extract contribution to melt flux by diurnal cycle only
        flx_melt_diurnal = flx_melt-flx_melt_0

        ! add refreezing flux from diurnal cycle to snow/ice heat flux
        g = g - flx_melt_diurnal
        flx_g = g

      endif

    else
      ! no diurnal cycle

      ! limit t_skin to <= T0 

      ! re-diagnose fluxes with t_skin = T0
      sh_0 = f_sh * (T0 - t2m)
      if (mask_snow.eq.1) then
        g_0  = snow_par%lambda/(0.5_wp*h_snow) * (T0 - tsoil)
      else
        g_0  = lambda_i/z(1) * (T0 - tsoil)
      endif
      lh_0 = f_lh * (qsat + dqsatdT*(T0-t_skin_old) - q2m)
      lw_0 = (1._wp-emiss)*lwdown + emiss*sigma * (t_skin_old**4 + 4._wp*t_skin_old**3*(T0-t_skin_old))

      ! diagnose remaining energy flux which can be used to melt snow/ice
      flx_melt = sh + lh + g + lw &
        - (sh_0+lh_0+g_0+lw_0)

      if (flx_melt.le.0._wp) then
        flx_melt = 0._wp
      else
        ! reset t_skin to melt temperature
        t_skin   = T0
        ! update fluxes
        lw = lw_0
        sh = sh_0
        lh = lh_0
        g  = g_0
        flx_g = g
      endif

    endif


    ! energy conservation check
    if( check_energy ) then
      energy_cons_surf1 = swnet &
      + lwdown        &
      - lw &
      - sh &
      - lh &
      - flx_g &
      - flx_melt 
      if( abs(energy_cons_surf1) .gt. 1.e-10_wp) then
        print *,''
        print *,'ebal'
        print *,'energy balance',energy_cons_surf1
        print *,'mask_snow',mask_snow
        print *,'h_snow',h_snow
        print *,'t_skin,t_skin_old',t_skin,t_skin_old
        print *,'t_skin_max,t_skin_min,t_skin_amp',t_skin_max,t_skin_min,t_skin_amp
        print *,'t_skin_pos,t_skin_neg',t_skin_pos,t_skin_neg
        print *,'dt_pos,dt_neg',dt_pos,dt_neg
        print *,'sw,lw_d,lw_u,sh,lh,g',swnet,emiss*lwdown,lw,sh,lh,g
        print *,'flx_melt',flx_melt
        print *,''
        !stop
      endif
    endif

    if( t_skin.gt.350._wp .or. t_skin.lt.150._wp ) then
      print *,''
      print *,i,j
      print *,'t_skin over ice out of range!!!',t_skin
      print *,'mask_snow',mask_snow
      print *,'h_snow',h_snow
      print *,'sw,lw_d,lw_u,sh,lh,g',swnet,emiss*lwdown,lw,sh,lh,g
      print *,'t_prof',t_prof
      print *,'t2m',t2m
      print *,'q2m',q2m
      print *,'r_a',r_a
      print *,''
      !stop
    endif


    return

  end subroutine ebal


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  u p d a t e _ t s k i n 
  !   Purpose    :  update skin temperature using new ground heat flux computed from new topsoil temperature 
  !              :  and re-diagnose surface energy fluxes
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine update_tskin(mask_snow,t_skin_old,dflxg_dT,t2m,q2m,swnet,lwdown, &
                         t_prof,t_prof_old,flx_melt, &
                         num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw, &
                         f_sh, f_e, f_lh, f_lw, qsat, dqsatdT, &
                         t_skin,flx_g,flx_sh,flx_lwu,flx_lh,evp,i,j)

    implicit none

    integer, intent(in) :: mask_snow
    real(wp), intent(in) :: t_skin_old, dflxg_dT
    real(wp), intent(in) :: t2m, q2m, swnet, lwdown
    real(wp), intent(in) :: flx_melt
    real(wp), dimension(0:), intent(in) :: t_prof, t_prof_old
    real(wp), intent(in) :: num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw
    real(wp), intent(in) :: f_sh, f_e, f_lh, f_lw, qsat, dqsatdT
    real(wp), intent(inout) :: t_skin, flx_g, flx_sh, flx_lwu, flx_lh, evp

    real(wp) :: num, denom
    real(wp) :: energy_cons_surf2

    integer :: i,j


    if (mask_snow.eq.1) then
      flx_g  = flx_g + dflxg_dT * (t_prof(0) - t_prof_old(0))
    else
      flx_g  = flx_g + dflxg_dT * (t_prof(1) - t_prof_old(1))
    endif
    num   = num_sw + num_lw + num_sh + num_lh  &
          - flx_g - flx_melt 

    denom = denom_lw + denom_sh + denom_lh 

    ! new skin temperature 
    t_skin   = num/denom

    ! diagnose surface energy fluxes from skin temperature 
    flx_sh  = f_sh * (t_skin - t2m)
    flx_lwu = f_lw * (t_skin_old**4 + 4._wp*t_skin_old**3*(t_skin-t_skin_old))
    flx_lh  = f_lh * (qsat + dqsatdT*(t_skin-t_skin_old) - q2m)

    ! evaporation/sublimation
    evp = f_e * (qsat + dqsatdT*(t_skin-t_skin_old) - q2m)

    if( check_energy ) then
      ! energy conservation check
      energy_cons_surf2 = swnet &
      + f_lw/sigma*lwdown        &
      - flx_lwu &
      - flx_sh &
      - flx_lh &
      - flx_g &
      - flx_melt 

          if( abs(energy_cons_surf2) .gt. 1.e-10_wp ) then
      !if (i.eq.9 .and. j.eq.105) then
            print *,''
            print *,'update_tskin'
            print *,'energy balance',energy_cons_surf2
        print *,'mask_snow',mask_snow
        print *,'t_skin,t_skin_old',t_skin,t_skin_old
      if( mask_snow.eq.1) then
        print *,'t_snow,t_snow_old',t_prof(0),t_prof_old(0)
      endif
      print *,'t_prof,t_prof_old',t_prof(1),t_prof_old(1)
        print *,'sw,lw_d,lw_u,sh,lh,g,flx_melt',swnet,f_lw/sigma*lwdown,flx_lwu,flx_sh,flx_lh,flx_g,flx_melt
        print *,''
        !stop
          endif

    endif

    if( t_skin.gt.350. .or. t_skin.lt.150. ) then
      print *,''
      print *,i,j
      print *,'t_skin over ice out of range update_tskin!!!',t_skin
      print *,'t_skin_old',t_skin_old
      print *,'mask_snow',mask_snow
      print *,'t2m',t2m
      if( mask_snow.eq.1) then
        print *,'t_snow,t_snow_old',t_prof(0),t_prof_old(0)
      endif
      print *,'t_prof,t_prof_old',t_prof(1),t_prof_old(1)
      print *,'sw',swnet
      print *,'lw',flx_lwu
      print *,'sh',flx_sh
      print *,'lh',flx_lh
      print *,'g',flx_g
      print *,'dflxg_dT',dflxg_dT
      print *,''
      !stop
    endif

    return

  end subroutine update_tskin


end module smb_ebal_mod

