!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s o i l _ c a r b o n _ p a r _ m o d
!
!  Purpose : soil carbon parameters
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
module soil_carbon_par_mod

  use precision, only : wp
  use timer, only : doy, mon, year, nday_mon, nmon_year, nstep_mon_lnd, time_eoy_lnd
  use constants, only : T0, k_boltz
  use control, only : lnd_restart
  use lnd_grid, only : z_int, z, dz, nl
  use lnd_grid, only : z_int_c, z_c, dz_c, nlc
  use lnd_params, only : soilc_par, peat_par, ch4_par, dust_par, nmonwet

  implicit none

  private
  public :: soil_carbon_par

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s o i l _ c a r b o n _ p a r
  !   Purpose    :  update soil carbon decomposition rate, diffusivity and advection
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine soil_carbon_par(f_veg,theta_field,theta_sat,litter_c_peat,acro_c,cato_c,dust_dep, &
                            t_soil_cum,theta_w_cum,theta_i_cum,psi,f_wet_cum,w_table_cum,f_wet_mon,f_wet_long,w_table_mon,t_soil_max, &
                            ftemp,fmoist,fdepth, &
                            k_litter,k_fast,k_slow,diff_soilc,adv_soilc,k_litter_wet,k_fast_wet,k_slow_wet, &
                            k_litter_peat,k_acro,k_cato,k_litter_peat_anox,k_acro_anox,ch4_frac_wet,ch4_frac_peat, &
                            f_peat_pot,f_oxic_peat,f_wetland,w_table_min,w_table_peat,alt,acro_h,cato_h,peat_c_ini_year)

    implicit none

    real(wp), intent(in) :: f_veg
    real(wp), dimension(:), intent(in) :: theta_field, theta_sat, psi
    real(wp), intent(in) :: litter_c_peat, acro_c
    real(wp), dimension(:), intent(in) :: cato_c
    real(wp), intent(in) :: dust_dep
    real(wp), intent(inout) :: f_wet_cum, w_table_cum
    real(wp), dimension(:), intent(inout) :: f_wet_mon, w_table_mon
    real(wp), dimension(:), intent(inout) :: f_wet_long
    real(wp), dimension(:), intent(inout) :: t_soil_cum, theta_w_cum, theta_i_cum, t_soil_max, ch4_frac_wet, ch4_frac_peat
    real(wp), dimension(:), intent(inout) :: ftemp, fmoist, fdepth
    real(wp), dimension(:), intent(inout) :: k_litter, k_fast, k_slow, k_litter_wet, k_fast_wet, k_slow_wet, k_cato, diff_soilc, adv_soilc
    real(wp), intent(inout) :: k_litter_peat, k_acro, k_litter_peat_anox, k_acro_anox
    real(wp), intent(inout) :: f_peat_pot, f_oxic_peat, f_wetland, w_table_min, w_table_peat, alt, acro_h, cato_h, peat_c_ini_year

    integer :: n, k, m, lsup, bubble
    real(wp), dimension(nl) :: t_soil_mean, theta_w_mean, theta_i_mean
    real(wp), dimension(nlc) :: adv, advsoilc, diff, diffsoilc
    real(wp), dimension(nl) :: fmoist_wet, ftemp_ch4
    real(wp), dimension(nl) :: klitter, kfast, kslow, klitter_wet, kfast_wet, kslow_wet
    real(wp) :: klitter_peat, kacro, klitter_peat_anox, kacro_anox
    real(wp), dimension(nl) :: kcato
    real(wp) :: temp, f_peat_pot_old
    real(wp), dimension(nmonwet) :: f_wet_ordered


    ! average cumulated monthly variables
    t_soil_mean      = t_soil_cum  / nstep_mon_lnd
    theta_w_mean     = theta_w_cum / nstep_mon_lnd
    theta_i_mean     = theta_i_cum / nstep_mon_lnd
    f_wet_mon(mon)   = f_wet_cum   / nstep_mon_lnd
    f_wetland        = f_wet_cum   / nstep_mon_lnd
    w_table_mon(mon) = w_table_cum / nstep_mon_lnd

    ! compute potential peatland fraction at end of year, based on monthly values of wetland extent
    if( time_eoy_lnd ) then

      ! shift to make space for current year, discard oldest year data
      f_wet_long(1:nmonwet-nmon_year) = f_wet_long(nmon_year+1:nmonwet)
      ! add monthly wetland fraction of the current year
      f_wet_long(nmonwet-nmon_year+1:nmonwet) = f_wet_mon

      if (lnd_restart .and. year.eq.1) then
        ! fill with current year
        do n=1,peat_par%nyearwet_peat
          f_wet_long((n-1)*nmon_year+1:nmon_year*n) = f_wet_mon
        enddo
      endif

      ! bubblesort f_wet_mon in ascending order
      f_wet_ordered = f_wet_long
      lsup = nmonwet
      do while (lsup > 1)
        bubble = 0 !bubble in the greatest element out of order
        do m = 1, (lsup-1)
          if (f_wet_ordered(m) > f_wet_ordered(m+1)) then
            temp = f_wet_ordered(m)
            f_wet_ordered(m) = f_wet_ordered(m+1)
            f_wet_ordered(m+1) = temp
            bubble = m
          endif
        enddo
        lsup = bubble
      enddo
      if (peat_par%peat_area) then
        f_peat_pot_old = f_peat_pot
        f_peat_pot = f_wet_ordered(nmonwet-peat_par%nmonwet_peat*peat_par%nyearwet_peat) ! in relative grid cell fraction
        f_peat_pot = max(peat_par%f_peat_min*f_veg,f_peat_pot)
      else
        f_peat_pot = 0._wp
      endif

      ! annual minimum water table depth
      w_table_min = minval(w_table_mon)
    endif

    ! min and max temperature for active layer thickness
    ! reset at first call of the year
    if( doy .eq. nday_mon ) then
      t_soil_max = 200._wp
    endif
    ! update values
    do k=1,nl
      t_soil_max(k) = max(t_soil_max(k),t_soil_mean(k))
    enddo
    ! compute active layer thickness at end of year
    if( time_eoy_lnd ) then
      alt = 0._wp
      if( t_soil_max(nl) .gt. T0 ) then
        alt = -1._wp
      else
        do k=1,nl
          if (t_soil_max(k) .eq. T0) alt = z(k)
          if (t_soil_max(k) .gt. T0) alt = z_int(k)
        enddo
      endif
    endif

    ! temperature dependence of decomposition rate
    do k=1,nl
      if (soilc_par%iresp_temp.eq.1) then ! Lloyd & Taylor 1994
        if( t_soil_mean(k) .gt. 240._wp ) then
          ftemp(k) = exp(308.56_wp * (1._wp/56.02_wp - 1._wp/(46.02_wp+t_soil_mean(k)-T0)) ) &
                   * 1._wp/(1._wp+exp(0.3_wp*(t_soil_mean(k)-(T0+soilc_par%temp_max_crit_resp))))  ! logistic suppression at high temperatures, reflecting enzyme denaturation
        else
          ftemp(k) = 0._wp
        endif
      elseif (soilc_par%iresp_temp.eq.2) then ! Arrhenius
        if( t_soil_mean(k) .gt. 260._wp ) then
          ftemp(k) = exp(soilc_par%Ea*(1._wp/(k_boltz*283.15_wp)-1._wp/(k_boltz*t_soil_mean(k)))) &
                   * 1._wp/(1._wp+exp(0.3_wp*(t_soil_mean(k)-(T0+soilc_par%temp_max_crit_resp))))  ! logistic suppression at high temperatures, reflecting enzyme denaturation
        else
          ftemp(k) = 0._wp
        endif
      elseif (soilc_par%iresp_temp.eq.3) then ! Q10
        ftemp(k) = soilc_par%q10_c**((t_soil_mean(k)-283.15_wp)/10._wp) &
                 * 1._wp/(1._wp+exp(0.3_wp*(t_soil_mean(k)-(T0+soilc_par%temp_max_crit_resp))))  ! logistic suppression at high temperatures, reflecting enzyme denaturation
      endif
    enddo

    ! temperature dependence of methane emission fraction of soil respiration
    if (ch4_par%ich4_ftemp.eq.1) then
      do k=1,nl
        ftemp_ch4(k) = exp(ch4_par%Ea_ch4*(1._wp/(k_boltz*303.15_wp)-1._wp/(k_boltz*t_soil_mean(k))))
      enddo
    else if (ch4_par%ich4_ftemp.eq.2) then
      ! Q10: Riley 2011, Kleinen 2019
      do k=1,nl
        ftemp_ch4(k) = ch4_par%q10_ch4**((t_soil_mean(k)-295._wp)/10._wp) 
      enddo
    endif
    ch4_frac_wet = ch4_par%ch4_frac_wet * ftemp_ch4
    ch4_frac_peat= ch4_par%ch4_frac_peat* ftemp_ch4

    ! soil moisture dependence of decomposition rate
    do k=1,nl
      if (soilc_par%iresp_moist.eq.0) then 
        fmoist(k) = 1._wp
      else if (soilc_par%iresp_moist.eq.1) then ! Porporato 2003
        if( theta_w_mean(k) .le. theta_field(k) ) then
          fmoist(k) = max(0._wp,theta_w_mean(k)-soilc_par%theta_crit_fmoist)/(theta_field(k)-soilc_par%theta_crit_fmoist) ! linear increase below field capacity
        else
          fmoist(k) = theta_field(k)/theta_w_mean(k) ! hyperbolic decrease above field capacity
        endif
      else if (soilc_par%iresp_moist.eq.2) then ! Koven 2013
        fmoist(k) = log(soilc_par%psi_min/psi(k))/log(soilc_par%psi_min/soilc_par%psi_max) 
        fmoist(k) = max(0._wp,fmoist(k))
        fmoist(k) = min(1._wp,fmoist(k))
      else if (soilc_par%iresp_moist.eq.3) then ! Porporato 2003, but using liquid + frozen water
        if( (theta_w_mean(k)+theta_i_mean(k)) .le. theta_field(k) ) then
          fmoist(k) = max(0._wp,(theta_w_mean(k)+theta_i_mean(k))-soilc_par%theta_crit_fmoist)/(theta_field(k)-soilc_par%theta_crit_fmoist) ! linear increase below field capacity
        else
          fmoist(k) = theta_field(k)/(theta_w_mean(k)+theta_i_mean(k)) ! hyperbolic decrease above field capacity
        endif
      endif
      fmoist_wet(k) = theta_field(k)/theta_sat(k)
    enddo

    ! soil depth dependence of decomposition rate, Koven 2013
    do k=1,nl
      fdepth(k) = exp(-z(k)/soilc_par%z_tau)
    enddo

    ! decomposition rate, 1/s
    klitter     = max(soilc_par%k_min, soilc_par%k10_litter * ftemp*fmoist*fdepth    )  
    kfast       = max(soilc_par%k_min, soilc_par%k10_fast   * ftemp*fmoist*fdepth    ) 
    kslow       = max(soilc_par%k_min, soilc_par%k10_slow   * ftemp*fmoist*fdepth    )
    klitter_wet = max(soilc_par%k_min, soilc_par%k10_litter * ftemp*fmoist_wet*fdepth) 
    kfast_wet   = max(soilc_par%k_min, soilc_par%k10_fast   * ftemp*fmoist_wet*fdepth) 
    kslow_wet   = max(soilc_par%k_min, soilc_par%k10_slow   * ftemp*fmoist_wet*fdepth) 

    if( peat_par%peat_carb ) then
      ! acrotelm thickness based on acrotelm carbon content and acrotelm carbon density, Kleinen 2012
      acro_h = acro_c / peat_par%rho_acro  ! m, kgC/m2 / (kgC/m3)
      ! catotelm thickness based on catotelm carbon content and catotelm carbon density, Kleinen 2012
      cato_h = sum(cato_c(1:nl)*dz(1:nl)) / peat_par%rho_cato  ! m, kgC/m3*m / (kgC/m3)
      ! oxic fraction of acrotelm decomposition in top layer based on water table depth 
      !f_oxic = min(max(0._wp, w_table_mon(mon,i,j)-w_table_min(i,j)-z_int(k-1)), dz(k)) / dz(k)
      !f_oxic = min(max(0._wp,w_table_mon(mon,i,j)-w_table_min(i,j)), dz(1)) / dz(1)
      !f_oxic_peat = min(max(0._wp,w_table_mon(mon)-w_table_min), 0.3_wp) / 0.3_wp ! acrotelm in the top 30 cm
      w_table_peat = min(max(0._wp,w_table_mon(mon)-w_table_min), acro_h)
      f_oxic_peat = w_table_peat / max(1.d-10,acro_h)

      klitter_peat = soilc_par%k10_litter * ftemp(1) &
        * (f_oxic_peat + (1._wp-f_oxic_peat)*peat_par%fmoist_peat)  ! litter, oxic and anoxic decomposition
      kacro        = peat_par%k10_acro   * ftemp(1) &
        * (f_oxic_peat + (1._wp-f_oxic_peat)*peat_par%fmoist_peat)  ! acrotelm, oxic and anoxic decomposition
      kcato        = peat_par%k10_cato   * ftemp*peat_par%fmoist_peat ! *fdepth    ! catotelm
      ! anoxic decomposition rate for methane emission
      klitter_peat_anox = soilc_par%k10_litter * ftemp(1) * peat_par%fmoist_peat
      kacro_anox        = peat_par%k10_acro   * ftemp(1) * peat_par%fmoist_peat

    endif

    ! diffusivity for soil carbon, m2/s 
    if( alt .gt. -1._wp ) then
      ! cryoturbation, Koven 2009, 2013, dependence on active layer thickness
      do k=1,nlc
        if( z_c(k) .lt. alt ) then
          diff(k) = soilc_par%diff_cryo * exp(-z(k)/soilc_par%z_diff)
        else if( z_c(k).ge.alt .and. z_c(k).lt.soilc_par%n_alt*alt ) then
          diff(k) = soilc_par%diff_cryo * (1._wp - (z_c(k)-alt)/((soilc_par%n_alt-1._wp)*alt)) * exp(-z_c(k)/soilc_par%z_diff)
        else
          diff(k) = soilc_par%diff_min
        endif
      enddo
    else
      ! bioturbation, Brakhekke 2011 based on mixing length theory
      diff = soilc_par%diff_bio * exp(-z_c(1:nlc)/soilc_par%z_diff)  ! m2/s
    endif
    ! carbon diffusivity at the soil levels (interfaces), m2/s
    do k=1,nl-1
      if( diff(k) .eq. 0._wp .and. diff(k+1) .eq. 0._wp ) then
        diffsoilc(k) = 0._wp
      else
        diffsoilc(k) = diff(k)*diff(k+1)*(z_c(k+1)-z_c(k))  &
          / (diff(k)*(z_c(k+1)-z_int_c(k)) + diff(k+1)*(z_int_c(k)-z_c(k)))
      endif
    enddo
    diffsoilc(nl:nlc) = 0._wp

    ! advection (sedimentation) for soil carbon, m/s 
    if (soilc_par%iadv_soilc.eq.0) then
      ! no advection
      adv = 0._wp
    else if (soilc_par%iadv_soilc.eq.1) then
      ! advection using dust deposition 
      adv = dust_dep/soilc_par%rho_loess
    else if (soilc_par%iadv_soilc.eq.2) then
      ! constant rate
      adv = soilc_par%adv_soil
    else if (soilc_par%iadv_soilc.eq.3) then
      ! advection using dust deposition + constant rate
      adv = soilc_par%adv_soil + dust_dep/soilc_par%rho_loess
    endif
    ! carbon advection at the soil levels (interfaces), m/s
    do k=1,nlc-1
      advsoilc(k) = 0.5_wp*(adv(k)+adv(k+1)) 
    enddo
    advsoilc(nlc) = 0._wp

    if (mon.eq.1) then
      diff_soilc = 0._wp
      adv_soilc = 0._wp
      k_litter = 0._wp 
      k_fast = 0._wp 
      k_slow = 0._wp
      k_litter_wet = 0._wp
      k_fast_wet = 0._wp
      k_slow_wet = 0._wp
      k_litter_peat = 0._wp
      k_acro = 0._wp
      k_cato = 0._wp
      k_litter_peat_anox = 0._wp
      k_acro_anox = 0._wp
    endif

    ! cumulate for soil carbon call
    diff_soilc = diff_soilc + diffsoilc ! diffusion
    adv_soilc  = adv_soilc + advsoilc   ! advection

    ! non-wetland
    k_litter(1:nl)   = k_litter(1:nl) + klitter
    k_fast(1:nl)     = k_fast(1:nl) + kfast
    k_slow(1:nl)     = k_slow(1:nl) + kslow

    ! wetland (exluding peatland)
    k_litter_wet(1:nl) = k_litter_wet(1:nl) + klitter_wet
    k_fast_wet(1:nl)   = k_fast_wet(1:nl) + kfast_wet
    k_slow_wet(1:nl)   = k_slow_wet(1:nl) + kslow_wet

    ! peatlands
    if( peat_par%peat_carb ) then
      k_litter_peat = k_litter_peat + klitter_peat
      k_acro     = k_acro + kacro
      k_cato(1:nl)     = k_cato(1:nl) + kcato
      k_litter_peat_anox = k_litter_peat_anox + klitter_peat_anox
      k_acro_anox        = k_acro_anox + kacro_anox

      ! save peat carbon at beginning of year
      if( mon .eq. 1) peat_c_ini_year = litter_c_peat+acro_c+sum(cato_c(1:nlc)*dz_c(1:nlc))
    endif

    ! no decomposition in the burial layer
    k_litter(nlc) = 0._wp
    k_fast(nlc)   = 0._wp
    k_slow(nlc)   = 0._wp
    k_litter_wet(nlc) = 0._wp
    k_fast_wet(nlc)   = 0._wp
    k_slow_wet(nlc)   = 0._wp
    k_cato(nlc) = 0._wp

    ! reset cumulated variables
    t_soil_cum  = 0._wp
    theta_w_cum = 0._wp
    theta_i_cum = 0._wp
    f_wet_cum   = 0._wp
    w_table_cum = 0._wp


    return

  end subroutine soil_carbon_par


end module soil_carbon_par_mod
