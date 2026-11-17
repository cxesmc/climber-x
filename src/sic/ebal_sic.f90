!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : e b a l _ s i c _ m o d
!
!  Purpose : energy balance over sea ice
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
module ebal_sic_mod

   use precision, only : wp
   use control, only: check_energy, check_water
   use constants, only : Le, Lf, Ls, q_sat_i, dqsat_dt_i
   use constants, only : cap_a, cap_w, rho_a, rho_w, emis_snow, sigma, T0
   use sic_params, only : rho_sic, rho_snow, dt, emis_sic, lambda_sic, lambda_snow, h_k_crit, u_star

   implicit none

   private
   public :: ebal_sic

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  e b a l _ s i c
  !   Purpose    :  compute skin temperature and diagnose surface energy fluxes
  !              :  over sea ice by solving the surface energy balance equation
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ebal_sic(t_ocn,t_freeze,t_air,qair,pressure,swnet,lwdown,snow,rain,wind,Cde,Cdh, &
                      t_skin,h_snow,h_sic, &
                      flx_ocn,flx_sh,flx_lwu,flx_lh,fw_ocn,evp, &
                      flx_melt_top, flx_melt_bot, &
                      dh_snow,dh_sic,i,j)

    implicit none

    real(wp), intent(in   ) :: t_ocn, t_freeze, t_air, qair
    real(wp), intent(in   ) :: pressure, lwdown, swnet, snow, rain, wind, Cde, Cdh
    real(wp), intent(in   ) :: h_snow, h_sic
    real(wp), intent(inout) :: t_skin
    real(wp), intent(out  ) :: flx_ocn, flx_sh, flx_lwu, flx_lh, fw_ocn, evp
    real(wp), intent(out  ) :: flx_melt_top, flx_melt_bot
    real(wp), intent(out  ) :: dh_snow, dh_sic
 
    integer :: i, j
    real(wp) :: sh, lh, g, lw, sh_new, lh_new, g_new, lw_new, t_skin_old
    real(wp) :: flux_melt, flux_snowmelt_top, flux_sicmelt_top, flux_excess
    real(wp) :: flux_bot_net, flux_snowmelt_bot, flux_sicmelt_bot
    real(wp) :: num, denom, num_g, denom_g, emiss, rhoa, k_scale, h_e
    real(wp) :: num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw
    real(wp) :: f_sh, f_e, f_lh, f_lw, f_g, qsat, dqsatdT, h_snow_tmp, h_sic_tmp
    real(wp) :: energy_cons_top, energy_cons_bot, water_cons

    real(wp), parameter :: Ch = 0.0058  ! exchange coefficient, McPhee 1992 
    real(wp), parameter :: e = exp(1._wp)
    real(wp), parameter :: eps = 0.1  ! m


    t_skin_old = t_skin

    ! heat conductivity scaling to represent sub-grid thickness distribution following Fichefet 1997, eq. 2 and 3
    ! scaling suppressed for large ice thicknesses
    h_e = lambda_snow*lambda_sic/(lambda_snow+lambda_sic)*(h_snow/lambda_snow+h_sic/lambda_sic) ! effective thickness
    if (h_e.gt.(e*eps/2._wp)) then
      k_scale = 1._wp+0.5_wp*log(2._wp*h_e/(e*eps))*exp(-max(0._wp,h_sic-h_k_crit)/5._wp)  
    else 
      k_scale = 1._wp
    endif

    if (h_snow .gt. 0._wp) then
      emiss = emis_snow
    else
      emiss = emis_sic
    endif
    qsat = q_sat_i(t_skin,pressure)
    dqsatdT = dqsat_dT_i(t_skin,pressure) 

    rhoa = rho_a(t_air,pressure)
    f_sh = Cdh*wind*rhoa*cap_a
    f_lh = Ls*Cde*wind*rhoa
    f_e  = Cde*wind*rhoa
    f_lw = emiss*sigma
    f_g  = lambda_snow*k_scale/(h_snow+(h_sic*lambda_snow/lambda_sic))

    ! sea ice heat flux
    num_g    = f_g * t_freeze
    denom_g  = f_g

    num_lh   = - f_lh * (qsat - dqsatdT*t_skin - qair)  ! surface sublimation
    denom_lh = f_lh * dqsatdT 

    num_sh   = f_sh * t_air
    denom_sh = f_sh

    num_sw   = swnet

    num_lw   = emiss*lwdown + 3._wp*emiss*sigma*t_skin**4
    denom_lw = 4._wp * emiss*sigma*t_skin**3

    num   = num_sw + num_lw + num_sh+ num_lh + num_g

    denom = denom_lw + denom_sh + denom_lh + denom_g

    ! new skin temperature        
    t_skin   = num/denom

    ! diagnose fluxes with the updated skin temperature
    sh = f_sh * (t_skin - t_air)
    lh = f_lh * (qsat + dqsatdT*(t_skin-t_skin_old) - qair)
    lw = (1._wp-emiss)*lwdown + emiss*sigma * (t_skin_old**4 + 4._wp*t_skin_old**3*(t_skin-t_skin_old))
    g  = f_g * (t_skin - t_freeze) 

    dh_snow = 0._wp
    dh_sic = 0._wp
    flx_ocn = 0._wp
    fw_ocn = 0._wp
    flux_melt = 0._wp
    flux_snowmelt_top = 0._wp
    flux_sicmelt_top  = 0._wp
    flux_snowmelt_bot = 0._wp
    flux_sicmelt_bot  = 0._wp

    ! add snowfall
    dh_snow = dh_snow + snow/rho_snow*dt

    if (t_skin .le. T0) then

      ! remove sublimation
      dh_snow = dh_snow - lh/(rho_snow*Ls)*dt
      if (h_snow+dh_snow .lt. 0._wp) then
        ! not enough snow to sublimate, sublimation from sea ice
        dh_sic = (h_snow + dh_snow)*rho_snow/rho_sic
        dh_snow = -h_snow  ! m
        if (h_sic+dh_sic .lt. 0._wp) then
          ! not enough ice to sublimate, evaporate from ocean and correct heat flux for (Ls-Le) difference
          flx_ocn = flx_ocn - (h_sic + dh_sic)*rho_sic/dt * (Ls-Le)  ! W/m2
          dh_sic = -h_sic  ! m
        endif
      endif

    else ! limit t_skin to <= 0°C

      ! reset t_skin to 0°C
      t_skin = T0

      ! re-diagnose fluxes with t_skin = T0
      sh_new = f_sh * (t_skin - t_air)
      lh_new = f_lh * (qsat + dqsatdT*(t_skin-t_skin_old) - qair)
      lw_new = (1._wp-emiss)*lwdown + emiss*sigma * (t_skin_old**4 + 4._wp*t_skin_old**3*(t_skin-t_skin_old))
      g_new  = f_g * (t_skin - t_freeze)

      ! diagnose remaining energy flux which can be used to melt snow or ice
      flux_melt = sh + lh + g + lw &
      - (sh_new+lh_new+g_new+lw_new)

      lw = lw_new
      sh = sh_new
      lh = lh_new
      g  = g_new

      ! remove sublimation
      dh_snow = dh_snow - lh/(rho_snow*Ls)*dt
      if (h_snow+dh_snow .lt. 0._wp) then
        ! not enough snow to sublimate, sublimation from sea ice
        dh_sic = (h_snow + dh_snow)*rho_snow/rho_sic
        dh_snow = -h_snow  ! m
        if (h_sic+dh_sic .lt. 0._wp) then
          ! not enough ice to sublimate, evaporate from ocean and correct heat flux for (Ls-Le) difference
          flx_ocn = flx_ocn - (h_sic + dh_sic)*rho_sic/dt * (Ls-Le)  ! W/m2
          dh_sic = -h_sic  ! m
        endif
      endif

      ! virtual snow and sea ice thickness before melting
      h_snow_tmp = h_snow + dh_snow
      h_sic_tmp  = h_sic + dh_sic
      ! melt snow and/or ice
      ! first melt snow and then ice if necessary
      dh_snow = dh_snow - flux_melt / (rho_snow*Lf) * dt  ! m
      flux_snowmelt_top = flux_melt
      if (h_snow+dh_snow .lt. 0._wp) then
        ! not enough snow to melt, melt all snow and use excess energy to melt ice below
        flux_snowmelt_top = h_snow_tmp * rho_snow*Lf / dt
        flux_excess = flux_melt - flux_snowmelt_top  ! W/m2
        dh_snow = -h_snow  ! m
        dh_sic = dh_sic - flux_excess / (rho_sic*Lf) * dt  ! m
        flux_sicmelt_top = flux_excess
        if (h_sic+dh_sic .lt. 0._wp) then
          ! not enough ice to melt, melt all ice and add excess energy to the sea ice conductive heat flux
          flux_sicmelt_top = h_sic_tmp * rho_sic*Lf / dt
          flux_excess = flux_melt - flux_snowmelt_top - flux_sicmelt_top  ! W/m2
          dh_sic = -h_sic  ! m
          g = g + flux_excess  ! W/m2
        endif
      endif

    endif

    ! ice ablation/accretion from bottom
    ! ocean heat flux, formulation from McPhee 1992, Weaver 2001 
    flx_ocn = flx_ocn + Ch*u_star * rho_w*cap_w * (t_freeze - t_ocn)
    ! net heat flux at the lower sea ice boundary
    ! used for ablation or accretion (if flux positive or negative, respectively)
    flux_bot_net = g - flx_ocn

    ! temporary snow and sea ice thickness before accretion/ablation
    h_snow_tmp = h_snow + dh_snow
    h_sic_tmp  = h_sic + dh_sic

    ! ablation/accretion of sea ice from the bottom
    dh_sic = dh_sic - flux_bot_net / (rho_sic*Lf) * dt  ! m

    flux_sicmelt_bot = flux_bot_net
    if (h_sic+dh_sic .lt. 0._wp) then
      ! all sea ice melting
      flux_sicmelt_bot = h_sic_tmp * rho_sic*Lf / dt
      flux_excess = flux_bot_net - flux_sicmelt_bot
      dh_sic = -h_sic
      dh_snow = dh_snow - flux_excess / (rho_snow*Lf) * dt  ! m
      flux_snowmelt_bot = flux_excess
      if (h_snow+dh_snow .lt. 0._wp) then
        ! not enough snow to melt, melt all snow and add excess energy to the ocean heat flux
        flux_snowmelt_bot = h_snow_tmp * rho_snow*Lf / dt
        flux_excess = flux_bot_net - flux_sicmelt_bot - flux_snowmelt_bot  ! W/m2
        dh_snow = -h_snow  ! m
        flx_ocn = flx_ocn + flux_excess  ! W/m2
      endif
    endif

    flx_sh  = sh
    flx_lwu = lw
    flx_lh  = lh
    evp     = lh/Ls 

    ! freshwater flux to the ocean
    fw_ocn  = rain + snow - evp - dh_snow*rho_snow/dt - dh_sic*rho_sic/dt  ! kg/m2/s

    ! melt fluxes for diagnostics
    flx_melt_top = flux_sicmelt_top + flux_snowmelt_top
    flx_melt_bot = flux_sicmelt_bot + flux_snowmelt_bot


    if( check_energy ) then
      ! energy conservation check
      energy_cons_top = swnet &
      + lwdown        &
      - lw &
      - sh &
      - lh &
      - g &
      - flux_snowmelt_top - flux_sicmelt_top
      if( abs(energy_cons_top) .gt. 1.d-9 ) then
        print *,''
        print *,'energy balance sic TOP',energy_cons_top
        print *,'tskin',t_skin
        print *,'t_air,t_ocn,t_freeze',t_air,t_ocn,t_freeze
        print *,'sw,lw_d,lw_u,sh,lh,g',swnet,emiss*lwdown,lw,sh,lh,g
        print *,'flux_melt',flux_melt
        print *,'flux_snowmelt_top,flux_sicmelt_top',flux_snowmelt_top,flux_sicmelt_top
        print *,'flux_bot_net',flux_bot_net
        print *,'flx_ocn',flx_ocn,Ch*u_star * rho_w*cap_w * (t_freeze - t_ocn)
        print *,'flux_snowmelt_bot,flux_sicmelt_bot',flux_snowmelt_bot,flux_sicmelt_bot
        print *,'h_sic,h_snow',h_sic,h_snow
        print *,'dh_sic,dh_snow',dh_sic,dh_snow
        print *,''
      endif

      energy_cons_bot = g &
      - flx_ocn      &
      - flux_snowmelt_bot - flux_sicmelt_bot 

      if( abs(energy_cons_bot) .gt. 1.d-10 ) then
        print *,''
        print *,'energy balance sic BOT',energy_cons_bot
        print *,'flx_ocn,g',flx_ocn,g
        print *,'flux_snowmelt,flux_sicmelt',flux_snowmelt_bot,flux_sicmelt_bot
        print *,'flux_bot_net',flux_bot_net
        print *,''
      endif

    endif

    if( check_water ) then
      ! water conservation check
      water_cons = dh_snow*rho_snow/dt &
      + dh_sic*rho_sic/dt &
      - snow &
      + evp &
      + (flux_snowmelt_top+flux_snowmelt_bot)/Lf &
      + (flux_sicmelt_top+flux_sicmelt_bot)/Lf
      if( abs(water_cons) .gt. 1.d-10 .and. (h_sic+dh_sic).gt.0._wp) then
        print *,''
        print *,'water balance sic',water_cons*dt
        print *,i,j
        print *,'tskin',t_skin
        print *,'snow',snow*dt
        print *,'evp',evp*dt
        print *,'dh_snow',dh_snow*rho_snow
        print *,'dh_sic',dh_sic*rho_sic
        print *,'h_sic,h_snow',h_sic*rho_sic,h_snow*rho_snow
        print *,'snowmelt',(flux_snowmelt_top+flux_snowmelt_bot)/Lf*dt
        print *,'sicmelt',(flux_sicmelt_top+flux_sicmelt_bot)/Lf*dt
        print *,'sicmelt_top',(flux_sicmelt_top)/Lf*dt
        print *,'sicmelt_bot',(flux_sicmelt_bot)/Lf*dt
      endif
    endif


    if( t_skin_old.gt.T0+40. .or. t_skin_old.lt.150._wp ) then
      print *,'t_skin_old over sic out of range!!!',t_skin_old,i,j
      print *,'t_skin',t_skin
      print *,'h_snow',h_snow
      print *,'sw,lw_d,lw_u,sh,lh,g,dflx_sic_dtsic',swnet,emiss*lwdown,lw,sh,lh,g
    endif

    if( t_skin.gt.T0 .or. t_skin.lt.150._wp ) then
      print *,''
      print *,'t_skin over sic out of range!!!',t_skin,i,j
      print *,'h_snow',h_snow
      print *,'sw,lw_d,lw_u,sh,lh,g',swnet,emiss*lwdown,lw,sh,lh,g
      print *,'tskinold',t_skin_old
      print *,'t_air',t_air
      print *,'t_ocn',t_ocn
      print *,'flx_ocn',flx_ocn
      print *,'qair',qair
      print *,'Cde,Cdh',Cde,Cdh
      print *,''
    endif

    return

  end subroutine ebal_sic

end module ebal_sic_mod

