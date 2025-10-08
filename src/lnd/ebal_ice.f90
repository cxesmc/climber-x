!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : e b a l _ i c e _ m o d
!
!  Purpose : energy balance over ice sheets
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
module ebal_ice_mod

   use precision, only : wp
   use constants, only : Ls, Lf, q_sat_i, dqsat_dt_i, T0, pi
   use constants, only : rho_a, cap_a, sigma
   use control, only : check_energy
   use lnd_grid, only : i_ice, z, nl
   use lnd_params, only : surf_par, hydro_par, l_diurnal_cycle, dt, rdt

   implicit none

   private
   public :: ebal_ice, update_tskin_ice

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  e b a l _ i c e
  !   Purpose    :  compute skin temperature and diagnose surface energy fluxes
  !              :  by solving the surface energy balance equation
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ebal_ice(mask_snow, h_snow, lambda_ice, &
                     t_skin, t_skin_old, t_ice, tatm, qatm, pressure, swnet, swnet_min, lwdown, &
                     beta_s, r_s, r_a, &
                     flx_g, dflxg_dT, flx_melt, t_skin_amp, &
                     num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw, &
                     f_sh, f_e, f_lh, f_lw, qsat, dqsatdT, &
                     energy_cons_surf1, i, j)

    implicit none

    integer, intent(in) :: i,j
    integer, intent(in) :: mask_snow
    real(wp), intent(in) :: h_snow
    real(wp), dimension(0:), intent(in) :: t_ice, lambda_ice
    real(wp), intent(in) :: tatm, qatm, swnet, swnet_min
    real(wp), intent(in) :: pressure, lwdown
    real(wp), intent(in) :: beta_s, r_s, r_a
    real(wp), intent(inout) :: t_skin
    real(wp), intent(inout) :: t_skin_old
    real(wp), intent(out) :: flx_g, dflxg_dT, flx_melt
    real(wp), intent(inout) :: t_skin_amp
    real(wp), intent(out) :: energy_cons_surf1
    real(wp), intent(inout) :: num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw
    real(wp), intent(inout) :: f_sh, f_e, f_lh, f_lw, qsat, dqsatdT

    real(wp) :: tice, p, rhoa
    real(wp) :: sh, lh, g, lw, sh_0, lh_0, g_0, lw_0
    real(wp) :: num, denom, num_g, denom_g, emiss
    real(wp) :: sh_pos, lh_pos, g_pos, lw_pos
    real(wp) :: t_skin_max, t_skin_min, t_skin_pos
    real(wp) :: t1, dt_pos, Ts, acos_fac, sqrt_fac


    flx_melt = 0._wp

    t_skin_old = t_skin
    if( mask_snow .eq. 0 ) then
      tice = t_ice(1) ! first ice layer temperature
    else
      tice = t_ice(0) ! snow temperature
    endif
    p = pressure

    if( mask_snow .eq. 1) then
      emiss = surf_par%emissivity_snow
    else
      emiss = surf_par%emissivity(i_ice)
    endif
    qsat = q_sat_i(t_skin_old,p)
    dqsatdT = dqsat_dT_i(t_skin_old,p) 

    rhoa = rho_a(tatm,p)
    f_sh = rhoa*cap_a/r_a
    f_lh = Ls * beta_s / (r_a + r_s) * rhoa
    f_e  = f_lh / Ls
    f_lw = emiss*sigma

    ! soil heat flux
    if( mask_snow .eq. 1) then
      num_g   = lambda_ice(0)/(0.5_wp*h_snow) * tice
      denom_g = lambda_ice(0)/(0.5_wp*h_snow)
    else
      num_g   = lambda_ice(1)/z(1) * tice
      denom_g = lambda_ice(1)/z(1)
    endif

    num_lh   = - f_lh * (qsat - dqsatdT*t_skin_old - qatm)  ! surface sublimation
    denom_lh = f_lh * dqsatdT 

    num_sh   = f_sh * tatm
    denom_sh = f_sh

    num_sw   = swnet

    num_lw   = emiss*lwdown + 3._wp*emiss*sigma*t_skin_old**4
    denom_lw = 4._wp * emiss*sigma*t_skin_old**3

    num   = num_sw + num_lw + num_sh+ num_lh + num_g

    denom = denom_lw + denom_sh + denom_lh + denom_g

    ! new skin temperature        
    t_skin   = num/denom

    ! diagnose fluxes 

    ! compute ice heat flux and derivative to be used as input to the ice temperature solver
    if( mask_snow .eq. 1) then
      flx_g = lambda_ice(0)/(0.5_wp*h_snow) * (t_skin - tice)
      dflxg_dT = - lambda_ice(0)/(0.5_wp*h_snow)
    else
      flx_g = lambda_ice(1)/z(1) * (t_skin - tice)
      dflxg_dT = - lambda_ice(1)/z(1)
    endif

    sh = f_sh * (t_skin - tatm)
    lh = f_lh * (qsat + dqsatdT*(t_skin-t_skin_old) - qatm)
    lw = (1._wp-emiss)*lwdown + emiss*sigma * (t_skin_old**4 + 4._wp*t_skin_old**3*(t_skin-t_skin_old))
    g  = flx_g

    if (.not.l_diurnal_cycle) then

      ! limit t_skin to <= 0°C 
      if (t_skin.gt.T0) then

        ! re-diagnose fluxes with t_skin = T0
        sh_0 = f_sh * (T0 - tatm)
        lh_0 = f_lh * (qsat + dqsatdT*(T0-t_skin_old) - qatm)
        if( mask_snow .eq. 1) then
          g_0 = lambda_ice(0)/(0.5_wp*h_snow) * (T0 - tice)
        else
          g_0 = lambda_ice(1)/z(1) * (T0 - tice)
        endif
        lw_0 = (1._wp-emiss)*lwdown + emiss*sigma * (t_skin_old**4 + 4._wp*t_skin_old**3*(T0-t_skin_old))

        ! diagnose remaining energy flux which can be used to melt snow/ice 
        flx_melt = sh + lh + g + lw &
          - (sh_0+lh_0+g_0+lw_0)

        if (flx_melt.le.0._wp) then
          flx_melt = 0._wp
        else 
          ! reset t_skin to 0°
          t_skin = T0
          ! update fluxes
          lw = lw_0
          sh = sh_0
          lh = lh_0
          g  = g_0
          flx_g = g
        endif

      endif

      t_skin_amp = 0._wp

    else ! parameterisation of diurnal cycle of skin temperature

      ! diurnal maximum skin temperature derived by conserving daily shortwave
      ! radiation integral, assuming sinusoidal diurnal cycle of SW
      t_skin_max = (2._wp*swnet-swnet_min + num_lw + num_sh + num_lh + num_g) / denom

      ! amplitude of diurnal cycle of skin temperature
      t_skin_amp = max(0._wp, (t_skin_max-t_skin))
      ! t_skin_max after scaling applied
      t_skin_max = t_skin+t_skin_amp
      ! diurnal minimum skin temperature (symmetric)
      t_skin_min = t_skin-t_skin_amp

      ! compute energy for melt from diurnal cycle, following (partly) Krapp et al. 2016 (SEMIC)
      if (t_skin_max.gt.T0) then
        ! re-diagnose fluxes with t_skin = T0
        sh_0 = f_sh * (T0 - tatm)
        lh_0 = f_lh * (qsat + dqsatdT*(T0-t_skin_old) - qatm)
        lw_0 = (1._wp-emiss)*lwdown + emiss*sigma * (t_skin_old**4 + 4._wp*t_skin_old**3*(T0-t_skin_old))
        if( mask_snow .eq. 1) then
          g_0 = lambda_ice(0)/(0.5_wp*h_snow) * (T0 - tice)
        else
          g_0 = lambda_ice(1)/z(1) * (T0 - tice)
        endif
      endif

      if (t_skin_min.gt.T0) then ! no need to resolve diurnal cycle, just use daily mean

        ! diagnose remaining energy flux which can be used to melt snow/ice
        flx_melt = sh + lh + g + lw &
          - (sh_0+lh_0+g_0+lw_0)

        if (flx_melt.le.0._wp) then
          flx_melt = 0._wp
        else 
          ! reset t_skin to 0°
          t_skin = T0
          ! update fluxes
          lw = lw_0
          sh = sh_0
          lh = lh_0
          g  = g_0
          flx_g = g
        endif

      else if (t_skin_max.gt.T0 .and. t_skin_amp>0._wp) then

        Ts = t_skin-T0
        acos_fac = acos(Ts/t_skin_amp)
        sqrt_fac = sqrt(1._wp-Ts**2/t_skin_amp**2)
        t1 = dt/(2._wp*pi) * acos_fac

        ! compute flux used to melt snow
        if (t_skin_max.gt.T0) then   
          ! compute mean skin temperature above T0
          dt_pos = dt-2._wp*t1 ! s, time that temperature is above freezing point
          t_skin_pos = T0 + dt/(pi*dt_pos)*(-Ts*acos_fac+t_skin_amp*sqrt_fac+pi*Ts)
          ! diagnose fluxes with t_skin = t_skin_pos
          sh_pos = f_sh * (t_skin_pos - tatm)
          if (mask_snow.eq.1) then
            g_pos  = lambda_ice(0)/(0.5_wp*h_snow) * (t_skin_pos - tice)
          else
            g_pos  = lambda_ice(1)/z(1) * (t_skin_pos - tice)
          endif
          lh_pos = f_lh * (qsat + dqsatdT*(t_skin_pos-t_skin_old) - qatm)
          lw_pos = (1._wp-emiss)*lwdown + emiss*sigma * (t_skin_old**4 + 4._wp*t_skin_old**3*(t_skin_pos-t_skin_old))
          ! diagnose energy flux which can be used to melt snow/ice
          flx_melt = sh_pos + lh_pos + g_pos + lw_pos &
            - (sh_0+lh_0+g_0+lw_0)
          ! scale with dt_pos
          flx_melt = flx_melt*dt_pos/dt
        endif

        ! rediagnose skin temperature substracting snowmelt energy
        t_skin   = (num_sw + num_lw + num_sh + num_lh + num_g - flx_melt) / denom

        ! diagnose fluxes with new skin temperature
        sh = f_sh * (t_skin - tatm)
        lh = f_lh * (qsat + dqsatdT*(t_skin-t_skin_old) - qatm)
        lw = (1._wp-emiss)*lwdown + emiss*sigma * (t_skin_old**4 + 4._wp*t_skin_old**3*(t_skin-t_skin_old))
        if(mask_snow.eq.1) then
          flx_g = lambda_ice(0)/(0.5_wp*h_snow) * (t_skin - tice)
        else
          flx_g = lambda_ice(1)/z(1) * (t_skin - tice)
        endif
        g = flx_g

      endif

    endif


    if( check_energy ) then
      ! energy conservation check
      energy_cons_surf1 = swnet &
      + lwdown        &
      - lw &
      - sh &
      - lh &
      - g &
      - flx_melt
      if( abs(energy_cons_surf1) .gt. 1.d-10 ) then
        print *,''
        print *,'energy balance',energy_cons_surf1,mask_snow
        print *,'sw,lw_d,lw_u,sh,lh,g',swnet,emiss*lwdown,lw,sh,lh,g,flx_melt
        print *,'t_skin,t_skin_old,tice',t_skin,t_skin_old,tice
        print *,''
      endif
    endif

    if( t_skin.gt.350._wp .or. t_skin.lt.100._wp ) then
      print *,''
      print *,'t_skin over ice out of range!!!',t_skin,i,j
      print *,'tskinold',t_skin_old
      print *,'mask_snow',mask_snow
      print *,'h_snow',h_snow
      print *,'sw,lw_d,lw_u,sh,lh,g',swnet,emiss*lwdown,lw,sh,lh,g
      print *,'lambda_ice',lambda_ice
      print *,'t_ice',t_ice
      print *,'tatm',tatm
      print *,'qatm',qatm
      print *,'beta_s',beta_s
      print *,'r_s',r_s
      print *,'r_a',r_a
      print *,''
      !stop
    endif


    return

  end subroutine ebal_ice


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  u p d a t e _ t s k i n 
  !   Purpose    :  update skin temperature using new ground heat flux computed from new topsoil temperature 
  !              :  and re-diagnose surface energy fluxes
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine update_tskin_ice(mask_snow, t_skin_old, dflxg_dT, tatm, qatm, swnet, lwdown, &
                         t_ice, t_ice_old, flx_g, flx_melt, & 
                         t_skin, flx_sh, flx_lwu, flx_lh, evap_surface, et, &
                         num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw, &
                         f_sh, f_e, f_lh, f_lw, qsat, dqsatdT, &
                         energy_cons_surf2, i, j)

    implicit none

    integer, intent(in) :: i, j
    integer, intent(in) :: mask_snow
    real(wp), intent(in) :: t_skin_old, dflxg_dT
    real(wp), intent(in) :: tatm, qatm, swnet
    real(wp), intent(in) :: lwdown
    real(wp), dimension(0:), intent(in) :: t_ice, t_ice_old
    real(wp), intent(inout) :: flx_g, flx_melt, t_skin
    real(wp), intent(inout) :: flx_sh, flx_lwu
    real(wp), intent(inout) :: flx_lh,evap_surface,et
    real(wp), intent(inout) :: num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw
    real(wp), intent(inout) :: f_sh, f_e, f_lh, f_lw, qsat, dqsatdT
    real(wp), intent(inout) :: energy_cons_surf2

    real(wp) :: num, denom


    if( mask_snow .eq. 1 ) then
      flx_g  = flx_g + dflxg_dT * (t_ice(0) - t_ice_old(0))
    else
      flx_g  = flx_g + dflxg_dT * (t_ice(1) - t_ice_old(1))
    endif
    num   = num_sw + num_lw + num_sh + num_lh  &
          - flx_g - flx_melt

    denom = denom_lw + denom_sh + denom_lh 

    ! new skin temperature 
    t_skin   = num/denom

    ! diagnose surface energy fluxes from skin temperature 
    flx_sh  = f_sh * (t_skin - tatm)
    flx_lwu = (1._wp-f_lw/sigma)*lwdown + f_lw * (t_skin_old**4 + 4._wp*t_skin_old**3*(t_skin-t_skin_old))
    flx_lh  = f_lh * (qsat + dqsatdT*(t_skin-t_skin_old) - qatm)

    evap_surface = f_e * (qsat + dqsatdT*(t_skin-t_skin_old) - qatm)

    et = evap_surface

    if( check_energy ) then
      ! energy conservation check
      energy_cons_surf2 = swnet &
      + lwdown        &
      - flx_lwu &
      - flx_sh &
      - flx_lh &
      - flx_g &
      - flx_melt
          if( abs(energy_cons_surf2) .gt. 1.e-10_wp ) then
            print *,''
            print *,'energy balance',energy_cons_surf2,mask_snow
          endif
    endif

    if( t_skin.gt.350. .or. t_skin.lt.100. ) then
      print *,''
      print *,'t_skin over ice out of range update_tskin!!!',t_skin,i,j
      print *,'t_skin_old',t_skin_old
      print *,'mask_snow',mask_snow
      print *,'tatm',tatm
      if( mask_snow.eq.1) then
        print *,'t_snow,t_snow_old',t_ice(0),t_ice_old(0)
      endif
      print *,'t_ice,t_ice_old',t_ice(1),t_ice_old(1)
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

  end subroutine update_tskin_ice

end module ebal_ice_mod

