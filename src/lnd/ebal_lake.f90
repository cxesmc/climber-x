!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : e b a l _ l a k e _ m o d
!
!  Purpose : energy balance over lakes
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
module ebal_lake_mod

   use precision, only : wp
   use constants, only : Le, Ls, Lf, q_sat_i, dqsat_dt_i, q_sat_w, dqsat_dt_w
   use constants, only : T0, rho_a, cap_a, sigma, pi
   use control, only : check_energy
   use lnd_grid, only : i_lake, i_ice, z_l
   use lnd_params, only : dt, rdt
   use lnd_params, only : surf_par, l_diurnal_cycle

   implicit none

   private
   public :: ebal_lake, update_tskin_lake

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  e b a l _ l a k e
  !   Purpose    :  compute surface energy fluxes
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ebal_lake(mask_snow, h_snow, lambda_lake, &
                     t_skin, t_skin_old, t_lake, tatm, qatm, pressure, swnet, lwdown, &
                     beta_s, r_s, r_a, &
                     flx_g, dflxg_dT, flx_melt, &
                     t_skin_amp, &
                     num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw, &
                     f_sh, f_e, f_le, f_lw, qsat_e, dqsatdT_e, &
                     energy_cons_surf1, i, j)

    implicit none

    integer, intent(in) :: i,j
    integer, intent(in) :: mask_snow
    real(wp), intent(in) :: h_snow
    real(wp), dimension(0:), intent(in) :: t_lake, lambda_lake
    real(wp), intent(in) :: tatm, qatm, swnet
    real(wp), intent(in) :: pressure, lwdown
    real(wp), intent(in) :: beta_s, r_s, r_a
    real(wp), intent(inout) :: t_skin
    real(wp), intent(inout) :: t_skin_old
    real(wp), intent(out) :: flx_g, dflxg_dT, flx_melt
    real(wp), intent(out) :: t_skin_amp
    real(wp), intent(out) :: energy_cons_surf1
    real(wp), intent(inout) :: num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw
    real(wp), intent(inout) :: f_sh, f_e, f_le, f_lw, qsat_e, dqsatdT_e
 

    integer :: n
    real(wp) :: p, rhoa, t_snow
    real(wp) :: sh, lh, g, lw, sh_0, lh_0, g_0, lw_0
    real(wp) :: num, denom, num_g, denom_g, emiss
    real(wp) :: L_v_skin


    n = i_lake

    t_skin_old = t_skin

    t_skin_amp = 0._wp
    flx_melt = 0._wp

    p = pressure

    if (mask_snow.eq.0) then

      ! top lake layer temperature is used to compute the lake-atmosphere fluxes

      t_skin = t_lake(1)

      if (t_skin.gt.0._wp) then
        emiss = surf_par%emissivity(i_lake)
        L_v_skin = Le
        qsat_e = q_sat_w(t_skin,p)
      else
        emiss = surf_par%emissivity_snow
        L_v_skin = Ls
        qsat_e = q_sat_i(t_skin,p)
      endif

      rhoa = rho_a(tatm,p)
      f_sh = rhoa*cap_a/r_a
      f_le = L_v_skin * rhoa * beta_s / (r_a + r_s)
      f_e  = f_le / L_v_skin
      f_lw = emiss*sigma

      sh = f_sh * (t_skin - tatm)
      lh = f_le * (qsat_e  - qatm)
      lw = (1._wp-emiss)*lwdown + emiss*sigma * t_skin**4 

      flx_g = swnet + lwdown - lw - sh - lh 
      dflxg_dT = 0._wp


    else ! mask_snow==1

      t_snow = t_lake(0) ! snow temperature

      emiss = surf_par%emissivity_snow
      qsat_e = q_sat_i(t_skin_old,p)
      dqsatdT_e = dqsat_dT_i(t_skin_old,p) 

      rhoa = rho_a(tatm,p)
      f_sh = rhoa*cap_a/r_a
      f_le = Ls * beta_s / (r_a + r_s) * rhoa
      f_e  = f_le / Ls
      f_lw = emiss*sigma

      num_g   = lambda_lake(0)/(0.5_wp*h_snow) * t_snow
      denom_g = lambda_lake(0)/(0.5_wp*h_snow)

      num_lh   = - f_le * (qsat_e - dqsatdT_e*t_skin_old - qatm)  ! surface sublimation
      denom_lh = f_le * dqsatdT_e 

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
      flx_g = lambda_lake(0)/(0.5_wp*h_snow) * (t_skin - t_snow)
      dflxg_dT = - lambda_lake(0)/(0.5_wp*h_snow)

      sh = f_sh * (t_skin - tatm)
      lh = f_le * (qsat_e + dqsatdT_e*(t_skin-t_skin_old) - qatm)
      lw = (1._wp-emiss)*lwdown + emiss*sigma * (t_skin_old**4 + 4._wp*t_skin_old**3*(t_skin-t_skin_old))
      g  = flx_g

      ! limit t_skin to <= 0°C if snow is on the lake
      if (t_skin.gt.T0) then

        ! re-diagnose fluxes with t_skin = T0
        sh_0 = f_sh * (T0 - tatm)
        lh_0 = f_le * (qsat_e + dqsatdT_e*(T0-t_skin_old) - qatm)
        g_0 = lambda_lake(0)/(0.5_wp*h_snow) * (T0 - t_snow)
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

    endif

    if( check_energy ) then
      ! energy conservation check
      energy_cons_surf1 = swnet &
      + lwdown        &
      - lw &
      - sh &
      - lh &
      - flx_g &
      - flx_melt
      if( abs(energy_cons_surf1) .gt. 1.d-10 ) then
        print *,''
        print *,'energy balance',n,energy_cons_surf1,mask_snow
        print *,'sw,lw_d,lw_u,sh,lh,flx_g',swnet,emiss*lwdown,lw,sh,lh,flx_g
        print *,'t_skin,t_skin_old,t_lake',t_skin,t_skin_old,t_lake(1)
        print *,''
      endif
    endif

    if( t_skin_old.gt.350._wp .or. t_skin_old.lt.190._wp ) then
      print *,'t_skin_old over lake out of range!!!',t_skin_old,i,j
      print *,'t_skin',t_skin
      print *,'lambda_lake',lambda_lake
      print *,'h_snow',h_snow
      print *,'sw,lw_d,lw_u,sh,lh,flx_g',swnet,emiss*lwdown,lw,sh,lh,flx_g

    endif

    if( t_skin.gt.350._wp .or. t_skin.lt.190._wp ) then
    !if( i.eq.21 .and. j.eq.33 ) then
      print *,''
      print *,'t_skin over lake out of range!!!',t_skin,i,j
      print *,'mask_snow',mask_snow
      print *,'h_snow',h_snow
      print *,'sw,lw_d,lw_u,sh,lh,flx_g',swnet,emiss*lwdown,lw,sh,lh,flx_g
      print *,'lambda_lake',lambda_lake
      print *,'tskinold',t_skin_old
      print *,'t_lake',t_lake
      print *,'tatm',tatm
      print *,'qatm',qatm
      print *,'beat_s',beta_s
      print *,'r_s',r_s
      print *,'r_a',r_a
      print *,''
      !stop
    endif


    return

  end subroutine ebal_lake


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  u p d a t e _ t s k i n 
  !   Purpose    :  update skin temperature using new ground heat flux computed from new topsoil temperature 
  !              :  and re-diagnose surface energy fluxes
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine update_tskin_lake(mask_snow, t_skin_old, dflxg_dT, tatm, qatm, swnet, lwdown, &
                         t_lake, t_lake_old, flx_g, flx_melt, & 
                         t_skin, flx_sh, flx_lwu, flx_lh, evap_surface, et, &
                         num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw, &
                         f_sh, f_e, f_le, f_lw, qsat_e, dqsatdT_e, &
                         energy_cons_surf2, i, j)

    implicit none

    integer, intent(in) :: i, j
    integer, intent(in) :: mask_snow
    real(wp), intent(in) :: t_skin_old, dflxg_dT
    real(wp), intent(in) :: tatm, qatm, swnet
    real(wp), intent(in) :: lwdown
    real(wp), dimension(0:), intent(in) :: t_lake, t_lake_old
    real(wp), intent(inout) :: flx_g, flx_melt, t_skin
    real(wp), intent(inout) :: flx_sh, flx_lwu
    real(wp), intent(inout) :: flx_lh,evap_surface,et
    real(wp), intent(inout) :: num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw
    real(wp), intent(inout) :: f_sh, f_e, f_le, f_lw, qsat_e, dqsatdT_e
    real(wp), intent(inout) :: energy_cons_surf2

    real(wp) :: num, denom


    if( mask_snow .eq. 1 ) then

      flx_g  = flx_g + dflxg_dT * (t_lake(0) - t_lake_old(0))
      num   = num_sw + num_lw + num_sh + num_lh  &
        - flx_g - flx_melt

      denom = denom_lw + denom_sh + denom_lh 

      ! new skin temperature 
      t_skin   = num/denom

      ! diagnose surface energy fluxes from skin temperature 
      flx_sh  = f_sh * (t_skin - tatm)
      flx_lwu = (1._wp-f_lw/sigma)*lwdown + f_lw * (t_skin_old**4 + 4._wp*t_skin_old**3*(t_skin-t_skin_old))
      flx_lh  = f_le * (qsat_e + dqsatdT_e*(t_skin-t_skin_old) - qatm)

      evap_surface = f_e * (qsat_e + dqsatdT_e*(t_skin-t_skin_old) - qatm)

    else

      flx_sh = f_sh * (t_skin - tatm)
      flx_lh = f_le * (qsat_e  - qatm)
      flx_lwu = (1._wp-f_lw/sigma)*lwdown + f_lw * t_skin**4 

      evap_surface = f_e * (qsat_e - qatm)

    endif

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
      print *,'t_skin over lake out of range update_tskin!!!',t_skin,i,j
      print *,'t_skin_old',t_skin_old
      print *,'mask_snow',mask_snow
      print *,'tatm',tatm
      if( mask_snow.eq.1) then
        print *,'t_snow,t_snow_old',t_lake(0),t_lake_old(0)
      endif
      print *,'t_lake,t_lake_old',t_lake(1),t_lake_old(1)
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

  end subroutine update_tskin_lake

end module ebal_lake_mod

