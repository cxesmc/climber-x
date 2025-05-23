!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : h y d r o l o g y _ m o d
!
!  Purpose : surface hydrology
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
module hydrology_mod

  use precision, only : wp
  use timer, only : sec_day
  use constants, only : q_sat_w, q_sat_i, Lf, rho_a, g, T0
  use control, only : check_water
  use lnd_grid, only : nsurf, npft, nveg, nsoil, i_ice, i_lake, is_veg, is_ice, is_lake
  use lnd_grid, only : flag_veg, flag_pft, flag_tree
  use lnd_grid, only : z_int, dz, nl
  use lnd_params, only : dt, rdt
  use lnd_params, only : snow_par, hydro_par

  implicit none

  private
  public :: canopy_water, surface_hydrology

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  c a n o p y _ w a t e r
  !   Purpose    :  solve analytical model of canopy water 
  !              :  including interception and throughfall
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine canopy_water(frac_surf,lai,sai,r_a,t_skin,pressure,qair,rain,snow, &
                         w_can,w_can_old,s_can,s_can_old, &
                         rain_ground,snow_ground,evap_can,subl_can,f_wat_can,f_snow_can)

    implicit none

    real(wp), dimension(:), intent(in) :: frac_surf, r_a, t_skin, qair
    real(wp), dimension(:), intent(in) :: lai, sai
    real(wp), dimension(:), intent(in) :: pressure
    real(wp), dimension(:), intent(in) :: rain, snow
    real(wp), dimension(:), intent(inout) :: w_can, w_can_old, s_can, s_can_old
    real(wp), dimension(:), intent(inout) :: rain_ground, snow_ground, evap_can, subl_can
    real(wp), dimension(:), intent(inout) :: f_wat_can, f_snow_can

    integer :: n
    real(wp) :: fac_e_w, fac_e_s, fac_i_w, fac_i_s, w_can_max, s_can_max
    real(wp) :: dw_can, rhoa
    real(wp) :: tau_s, fac_lai
    logical :: flag_w, flag_s


    do n=1,npft

      if (frac_surf(n).gt.0._wp) then

        if (hydro_par%l_prc_intercept) then
          flag_w = (rain(n).gt.0._wp .or. w_can(n).gt.0._wp)
          flag_s = (snow(n).gt.0._wp .or. s_can(n).gt.0._wp)
        else
          flag_w = .false.
          flag_s = .false.
        endif

        ! check if fac_lai needed and avoid computing it twice
        if( flag_w .or. flag_s ) then
          fac_lai = (1._wp - exp(-0.5_wp*(lai(n)+sai(n))))
        endif

        w_can_old(n) = w_can(n)
        s_can_old(n) = s_can(n)

        ! canopy liquid water interception only for trees
        if( flag_w ) then

          w_can_max = hydro_par%can_max_w * ( lai(n) + sai(n) ) ! maximum canopy water, kg/m2

          rhoa = rho_a(t_skin(n),pressure(n))
          if( .not. hydro_par%l_dew ) then ! exclude dew deposition (negative evaporation/sublimation)
            fac_e_w = rhoa/r_a(n) * max(0._wp, (q_sat_w(t_skin(n),pressure(n)) - qair(n))) ! evaporation factor
          else ! allow dew deposition
            fac_e_w = rhoa/r_a(n) * (q_sat_w(t_skin(n),pressure(n)) - qair(n)) ! evaporation factor
          endif

          fac_i_w = hydro_par%alpha_int_w(n) * fac_lai ! interception factor for water

          ! update canopy water
          w_can(n) = (fac_i_w*rain(n) + w_can(n)*rdt) / (1._wp*rdt + fac_e_w/w_can_max + 1._wp/hydro_par%tau_w)
          if( w_can(n) .lt. 0._wp ) then
            dw_can = -w_can(n)
            w_can(n) = 0._wp
          endif
          if( w_can(n) .lt. 1.e-30_wp ) w_can(n) = 0._wp
          if( w_can(n) .gt. w_can_max ) w_can(n) = w_can_max

          evap_can(n) = fac_e_w * w_can(n)/w_can_max

          rain_ground(n) = rain(n) - evap_can(n) - (w_can(n) - w_can_old(n)) * rdt

          ! fraction of water-covered canopy
          f_wat_can(n) = w_can(n)/w_can_max

        else

          rain_ground(n) = rain(n)
          evap_can(n) = 0._wp
          f_wat_can(n) = 0._wp

        endif


        ! snow interception
        if( flag_s ) then

          s_can_max = hydro_par%can_max_s * ( lai(n) + sai(n) ) ! maximum canopy snow, kg/m2

          rhoa = rho_a(t_skin(n),pressure(n))
          if( .not. hydro_par%l_dew ) then ! exclude dew deposition (negative evaporation/sublimation)
            fac_e_s = rhoa/r_a(n) * max(0._wp, (q_sat_i(t_skin(n),pressure(n)) - qair(n))) ! sublimation factor
          else ! allow dew deposition
            fac_e_s = rhoa/r_a(n) * (q_sat_i(t_skin(n),pressure(n)) - qair(n)) ! sublimation factor
          endif

          fac_i_s = hydro_par%alpha_int_s * fac_lai ! interception factor for snow

          ! update canopy snow
          if( t_skin(n) .lt. T0 ) then
            tau_s = 10._wp*hydro_par%tau_s
          else
            tau_s = 1._wp*hydro_par%tau_s
          endif
          s_can(n) = (fac_i_s*snow(n) + s_can(n)*rdt) / (1._wp*rdt + fac_e_s/s_can_max + 1._wp/tau_s)
          if( s_can(n) .lt. 1.e-20_wp ) s_can(n) = 0._wp
          if( s_can(n) .gt. s_can_max ) s_can(n) = s_can_max  

          subl_can(n) = fac_e_s * s_can(n)/s_can_max

          snow_ground(n) = snow(n) - subl_can(n) - (s_can(n) - s_can_old(n)) * rdt

          ! fraction of snow-covered canopy
          f_snow_can(n) = s_can(n)/s_can_max

        else

          snow_ground(n) = snow(n)
          subl_can(n) = 0._wp
          f_snow_can(n) = 0._wp

        endif


        ! reset negative rain
        if( rain_ground(n) .lt. 0._wp ) then
          if( check_water .and. rain_ground(n).lt.-1.e-10_wp) print *,'rain < 0 ',n,rain_ground(n)*dt,rain(n)*dt
          rain_ground(n) = 0._wp
        endif

        ! reset negative snow
        if( snow_ground(n) .lt. 0._wp ) then
          if( check_water .and. snow_ground(n).lt.-1.e-10_wp) print *,'snow < 0 ',snow_ground(n)
          snow_ground(n) = 0._wp
        endif


      else ! PFT not present

        evap_can(n) = 0._wp
        subl_can(n) = 0._wp
        w_can(n)    = 0._wp
        s_can(n)    = 0._wp
        f_wat_can(n) = 0._wp
        f_snow_can(n) = 0._wp

      endif

    enddo

    where (flag_pft.eq.0) 
        rain_ground = rain
        snow_ground = snow
        evap_can = 0._wp
        subl_can = 0._wp
        f_wat_can = 0._wp
        f_snow_can = 0._wp
    endwhere

    return

  end subroutine canopy_water


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s u r f a c e _ h y d r o l o g y
  !   Purpose    :  surface hydrology
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine surface_hydrology(frac_surf,mask_snow,evap_surface,rain_ground,snow_ground,snowmelt,icemelt, &
                              theta,theta_sat,theta_field,k_sat,cap_soil,cap_ice,cap_lake,f_peat,w_table_peat, &
                              cti_mean, cti_cdf, &
                              dyptop_k,dyptop_v,dyptop_xm,dyptop_fmax, &
                              w_snow_old,w_snow,w_snow_max,w_w,w_i,w_table_cum,f_wet_cum,tatm,t_soil,t_ice,t_lake, &
                              h_snow,calving,runoff_sur,infiltration,w_table,f_wet,f_wet_max,cti_lim,lake_water_tendency)


    implicit none

    integer, dimension(:), intent(inout) :: mask_snow
    real(wp), dimension(:), intent(in) :: frac_surf, evap_surface, snow_ground, tatm
    real(wp), dimension(:), intent(inout) :: rain_ground
    real(wp), intent(in) :: f_peat, w_table_peat
    real(wp), dimension(:), intent(in) :: theta, theta_sat, theta_field, k_sat
    real(wp), intent(in) :: cap_soil, cap_ice, cap_lake
    real(wp), intent(in) :: cti_mean, cti_cdf(:)
    real(wp), intent(in) :: dyptop_k,dyptop_v,dyptop_xm,dyptop_fmax
    real(wp), dimension(:), intent(inout) :: w_snow_old, w_snow, w_snow_max, snowmelt, icemelt
    real(wp), intent(inout) :: w_table_cum, f_wet_cum
    real(wp), dimension(0:), intent(inout) :: t_soil, t_ice, t_lake
    real(wp), dimension(:), intent(inout) :: w_w, w_i
    real(wp), dimension(:), intent(out) :: h_snow, calving, runoff_sur
    real(wp), intent(out) :: lake_water_tendency
    real(wp), intent(out) :: infiltration, w_table, f_wet, f_wet_max, cti_lim

    integer :: n, i1, i2
    real(wp) :: f_sat, w1, w2, phi, rain_g, snow_g, H, H_m, H_star, wsnowold_m
    real(wp) :: z, d_z, deficit, integral, psi
    real(wp) :: f_ice_grd, f_veg, f_lake
    real(wp) :: infiltration_max, q_liq, sublimation, evaporation, dws


    runoff_sur = 0._wp
    calving = 0._wp
    infiltration = 0._wp
    w_table = 0._wp
    f_wet = 0._wp
    f_wet_max = 0._wp
    cti_lim = 0._wp
    lake_water_tendency = 0._wp

    f_ice_grd = frac_surf(i_ice)
    f_veg = sum(frac_surf,mask=flag_veg.eq.1)
    f_lake = frac_surf(i_lake)

    !************************
    ! vegetated grid part
    !************************
    if( f_veg .gt. 0._wp ) then

      sublimation = 0._wp
      if( mask_snow(is_veg) .eq. 1 ) then
        do n=1,nveg
          ! mean sublimation from snow
          if (frac_surf(n).gt.0._wp) sublimation = sublimation + evap_surface(n)*frac_surf(n)/f_veg  ! kg/m2/s 
        enddo
      endif

      ! mean snow on the ground over vegetated part
      snow_g = 0._wp
      do n=1,nsurf
        if (flag_veg(n).eq.1 .and. frac_surf(n).gt.0._wp) snow_g = snow_g + snow_ground(n) * frac_surf(n)/f_veg
      enddo

      ! add snowfall to snow layer and remove sublimation, snowmelt already removed during soil temperature update
      w_snow(is_veg) = w_snow(is_veg) + snow_g * dt - sublimation * dt  ! kg/m2

      ! if snowmass below critical snow mass for explicit snow layer, use possible first soil layer excess energy to melt snow
      ! and update snowmelt
      if (w_snow(is_veg).gt.0._wp .and. w_snow(is_veg).lt.snow_par%w_snow_crit .and. t_soil(1).gt.T0) then
        H = cap_soil*rdt*(t_soil(1) - T0)    ! W/m2, energy available to melt snow
        wsnowold_m = w_snow(is_veg)
        H_m = H*dt/Lf ! kg/m2, snow that can be melted
        ! update w_snow
        w_snow(is_veg) = max( 0._wp, w_snow(is_veg)-H_m )  ! kg/m2
        H_star = H - Lf*rdt * (wsnowold_m - w_snow(is_veg)) ! heat not used to melt snow
        t_soil(1) = T0 + dt/cap_soil * H_star
        ! update snowmelt
        snowmelt(is_veg) = snowmelt(is_veg) + (wsnowold_m - w_snow(is_veg)) * rdt  ! kg/m2/s
      endif

      ! if w_snow negative, reset to 0 and remove required ice from the top soil layer (sublimation)
      if( w_snow(is_veg) .lt. 0._wp ) then
        ! if not enough ice and water in first layer, remove also ice from second layer
        if(-w_snow(is_veg) .gt. (w_i(1)+w_w(1))) then
          if (check_water) print *,'WARNING: not enough ice or water to sublimate in top layer!', w_snow(is_veg),w_i(1),w_w(1)
          dws = -(w_snow(is_veg)+w_i(1)+w_w(1))
          w_i(2) = w_i(2) - min(w_i(2),dws)
          w_i(1) = 0._wp
          w_w(1) = 0._wp
          ! if not enough ice in first layer, remove liquid water instead to keep water balance
        elseif(-w_snow(is_veg) .gt. w_i(1)) then
          if (check_water) print *,'WARNING: not enough ice to sublimate in top layer, sublimate also liquid water!', w_snow(is_veg),w_i(1)
          dws = -(w_snow(is_veg)+w_i(1))
          w_w(1) = w_w(1) - min(w_w(1),dws)
          w_i(1) = 0._wp
        else
          ! enough ice to sublimate in top layer
          w_i(1) = w_i(1) + w_snow(is_veg)
        endif
        w_snow(is_veg) = 0._wp
      endif

      ! limit w_snow and add to 'calving'
      if( w_snow(is_veg) .gt. snow_par%w_snow_off ) then
        calving(is_veg) = (w_snow(is_veg) - snow_par%w_snow_off) * rdt ! kg/m2/s
        w_snow(is_veg) = snow_par%w_snow_off
      endif
      ! update snow height
      h_snow(is_veg) = w_snow(is_veg) / snow_par%rho_snow

      ! save seasonal maximum snow swe
      if (w_snow(is_veg).gt.w_snow_old(is_veg)) w_snow_max(is_veg) = w_snow(is_veg)

      ! water table depth, from soil water content
      if (hydro_par%i_wtab.eq.1) then

        w_table = z_int(nl) - sum(theta / theta_sat * dz(1:nl))

      else if (hydro_par%i_wtab.eq.2) then

        ! water table after Niu et al 2005, section 2.4
        ! column soil water deficit
        deficit = sum((theta_sat-theta) * dz(1:nl))
        d_z = 0.1_wp
        w_table = 0._wp
        do while (integral<deficit)
          integral = 0._wp      
          z = 0._wp
          do while (z<w_table)
            integral = integral + (theta_sat(1)-theta_sat(1)*((-0.2-(w_table-z))/-0.2)**(-1._wp/6.)) * d_z
            z = z+d_z
          enddo
          w_table = w_table + d_z
        enddo

      else if (hydro_par%i_wtab.eq.3) then

        w_table = hydro_par%wtab_scale * (1._wp-theta(1)/theta_sat(1))

      elseif (hydro_par%i_wtab.eq.4) then

        w_table = z_int(nl) - sum(theta(1:nl-1) / theta_sat(1:nl-1) * dz(1:nl-1)) * sum(dz(1:nl))/sum(dz(1:nl-1))

      else if (hydro_par%i_wtab.eq.5) then

        w_table = z_int(nl) - 0.5_wp * (sum(theta/theta_sat*dz(1:nl)) + theta(1)/theta_sat(1)*sum(dz(1:nl)))

      else if (hydro_par%i_wtab.eq.6) then

        w_table = z_int(nl) - 0.5_wp * (sum(theta/theta_field*dz(1:nl)) + theta(1)/theta_field(1)*sum(dz(1:nl)))

      else if (hydro_par%i_wtab.eq.7) then

        w_table = z_int(nl) - theta(1)/theta_sat(1)*sum(dz(1:nl))

      else if (hydro_par%i_wtab.eq.8) then

        ! Kleinen 2020
        !k = nl
        !do while (psi.lt.0.95_wp .and. k.le.1) 
        !  psi = theta(k)/theta_field(k)
        !  k = k-1
        !enddo
        !if (k.eq.0) then
        !  w_table = 0._wp
        !else
        !  w_table = z_int(k) - psi*dz(k+1)
        !endif

      endif

      w_table_cum = w_table_cum + w_table
        
      ! max possible wetland extent (w_table=0)
      if (cti_mean.gt.14._wp) then
        f_wet_max = 0._wp
      else
        i1 = max(cti_mean,hydro_par%cti_min)
        i2 = i1+1
        w2 = (max(cti_mean,hydro_par%cti_min)-i1)/(i2-i1)
        w1 = 1._wp-w2
        f_wet_max = 1._wp-(w1*cti_cdf(i1)+w2*cti_cdf(i2))
      endif
      ! no inundation if CTI lower than critical value (5.5 in Kleinen 2020)
      if (cti_mean.le.hydro_par%cti_mean_crit) f_wet_max = 0._wp 

      ! saturated grid cell fraction
      if( hydro_par%i_fwet .eq. 1 ) then

        ! fraction of icefree surface at saturation, SIMTOP, Niu 2005
        f_sat = f_wet_max * exp(-hydro_par%f_wtab * w_table)
        ! wetland area, excluding snow
        if( mask_snow(is_veg) .eq. 1 ) then
          f_wet = 0._wp ! no wetland where snow
        else
          f_wet = f_sat
        endif

      elseif( hydro_par%i_fwet .eq. 2 ) then

        ! DYPTOP, Stocker 2014
        phi = (1._wp+dyptop_v*exp(-dyptop_k*(-w_table-dyptop_xm)))**(-1._wp/dyptop_v)
        f_sat = min(dyptop_fmax,phi)
        ! wetland area, excluding snow
        if( mask_snow(is_veg) .eq. 1 ) then
          f_wet = 0._wp ! no wetland where snow
        else
          if( dyptop_fmax .gt. hydro_par%fmax_crit ) then
            f_wet = f_sat
          else
            f_wet = 0._wp
          endif
        endif

      elseif( hydro_par%i_fwet .eq. 3 ) then

        ! TOPMODEL following Kleinen et al 2020
        cti_lim = cti_mean + hydro_par%f_wtab*w_table   ! NOTE: w_table is positive!
        cti_lim = max(cti_lim,hydro_par%cti_min)
        cti_lim = max(1._wp,cti_lim)
        if (cti_lim.gt.14._wp) then
          f_sat = 0._wp
        else
          i1 = cti_lim
          i2 = i1+1
          w2 = (cti_lim-i1)/(i2-i1)
          w1 = 1._wp-w2
          f_sat = 1._wp-(w1*cti_cdf(i1)+w2*cti_cdf(i2))
        endif
        ! no inundation if CTI lower than critical value (5.5 in Kleinen 2020)
        if (cti_mean.le.hydro_par%cti_mean_crit) f_sat = 0._wp 
        if( mask_snow(is_veg) .eq. 1 ) then
          f_wet = 0._wp ! no wetland where snow
        else
          f_wet = f_sat
        endif

      elseif( hydro_par%i_fwet .eq. 4 ) then

        ! TOPMODEL following Kleinen et al 2020
        cti_lim = hydro_par%cti_min + hydro_par%f_wtab*w_table   ! NOTE: w_table is positive!
        cti_lim = max(1._wp,cti_lim)
        if (cti_lim.gt.14._wp) then
          f_sat = 0._wp
        else
          i1 = cti_lim
          i2 = i1+1
          w2 = (cti_lim-i1)/(i2-i1)
          w1 = 1._wp-w2
          f_sat = 1._wp-(w1*cti_cdf(i1)+w2*cti_cdf(i2))
        endif
        ! no inundation if CTI lower than critical value (5.5 in Kleinen 2020)
        if (cti_mean.le.hydro_par%cti_mean_crit) f_sat = 0._wp 
        if( mask_snow(is_veg) .eq. 1 ) then
          f_wet = 0._wp ! no wetland where snow
        else
          f_wet = f_sat
        endif
        
      endif

      f_wet_cum = f_wet_cum + f_wet

      infiltration_max = k_sat(1) * (1._wp-f_sat)

      ! mean rain on the ground over vegetated part
      rain_g = 0._wp
      do n=1,nsurf
        if (flag_veg(n).eq.1 .and. frac_surf(n).gt.0._wp) rain_g = rain_g + rain_ground(n) * frac_surf(n)/f_veg
      enddo
      q_liq = rain_g + snowmelt(is_veg) ! kg/m2/s

      ! surface runoff, kg/m2/s
      runoff_sur(is_veg) = f_sat * q_liq &  ! all into runoff over saturated or frozen surface
        + (1._wp - f_sat) * max(0._wp, q_liq - infiltration_max)    ! this condition is never met!
      ! soil liquid water infiltration, kg/m2/s
      infiltration = rain_g + snowmelt(is_veg) - runoff_sur(is_veg)  ! kg/m2/s

    endif

    !************************
    ! ice
    !************************
    if( f_ice_grd .gt. 0._wp ) then

      sublimation = evap_surface(i_ice)

      ! add snowfall to snow layer and remove sublimation, snowmelt already removed during soil temperature update
      w_snow(is_ice) = w_snow(is_ice) + snow_ground(i_ice) * dt - sublimation * dt  ! kg/m2
      if( check_water .and. w_snow(is_ice).lt.0._wp ) print *,'WARNING w_snow < 0 over ice',w_snow(is_ice)
      w_snow(is_ice) = max(0._wp,w_snow(is_ice))  ! not strictly conserving water here!

      ! limit w_snow and add to 'calving'
      if( w_snow(is_ice) .gt. snow_par%w_snow_off ) then
        calving(is_ice) = (w_snow(is_ice) - snow_par%w_snow_off) * rdt ! kg/m2/s
        w_snow(is_ice) = snow_par%w_snow_off
      endif
      ! update snow height
      h_snow(is_ice) = w_snow(is_ice) / snow_par%rho_snow

      ! save seasonal maximum snow swe
      if (w_snow(is_ice).gt.w_snow_old(is_ice)) w_snow_max(is_ice) = w_snow(is_ice)

      ! total liquid water runoff 
      if (hydro_par%l_runoff_icemelt) then
        runoff_sur(is_ice) = rain_ground(i_ice) + snowmelt(is_ice) + icemelt(is_ice) ! kg/m2/s 
      else
        runoff_sur(is_ice) = rain_ground(i_ice) + snowmelt(is_ice) ! kg/m2/s 
      endif

      ! update snow mask
      if( w_snow(is_ice) .gt. snow_par%w_snow_crit ) then
        mask_snow(is_ice) = 1
      else
        mask_snow(is_ice) = 0
      endif

    endif

    !************************
    ! lake
    !************************
    if( f_lake .gt. 0._wp ) then

      if( mask_snow(is_lake) .eq. 1 ) then
        sublimation = evap_surface(i_lake)
        evaporation = 0._wp
      else
        sublimation = 0._wp
        evaporation = evap_surface(i_lake)
      endif

      ! add snowfall to snow layer and remove snowmelt and sublimation
      w_snow(is_lake) = w_snow(is_lake) + snow_ground(i_lake) * dt - snowmelt(is_lake)*dt - sublimation * dt  ! kg/m2

      ! if snowmass below critical snow mass for explicit snow layer, use possible first lake layer excess energy to melt snow
      ! and update snowmelt
      if (w_snow(is_lake).gt.0._wp .and. w_snow(is_lake).lt.snow_par%w_snow_crit .and. t_lake(1).gt.T0) then
        H = cap_lake*rdt*(t_lake(1) - T0)    ! W/m2, energy available to melt snow
        wsnowold_m = w_snow(is_lake)
        H_m = H*dt/Lf ! kg/m2, snow that can be melted
        ! update w_snow
        w_snow(is_lake) = max( 0._wp, w_snow(is_lake)-H_m )  ! kg/m2
        H_star = H - Lf*rdt * (wsnowold_m - w_snow(is_lake)) ! heat not used to melt snow
        t_lake(1) = T0 + dt/cap_lake * H_star
        ! update snowmelt
        snowmelt(is_lake) = snowmelt(is_lake) + (wsnowold_m - w_snow(is_lake)) * rdt  ! kg/m2/s
      endif

      if( check_water .and. w_snow(is_lake).lt.0._wp ) print *,'WARNING w_snow < 0 over lake',w_snow(is_lake)
      w_snow(is_lake) = max(0._wp,w_snow(is_lake))  ! not strictly conserving water here!

      ! limit w_snow and add to 'calving' 
      if( w_snow(is_lake) .gt. snow_par%w_snow_off ) then
        calving(is_lake) = (w_snow(is_lake) - snow_par%w_snow_off) * rdt ! kg/m2/s
        w_snow(is_lake) = snow_par%w_snow_off
      endif
      ! save seasonal maximum snow swe
      if (w_snow(is_lake).gt.w_snow_old(is_lake)) w_snow_max(is_lake) = w_snow(is_lake)

      ! update snow height
      h_snow(is_lake) = w_snow(is_lake) / snow_par%rho_snow

      ! net lake surface water balance
      lake_water_tendency = rain_ground(i_lake) + snowmelt(is_lake) - evaporation ! kg/m2/s

      ! runoff from lakes is computed in the coupler!
      runoff_sur(is_lake) = 0._wp 

      ! update snow mask
      if( w_snow(is_lake) .gt. snow_par%w_snow_crit ) then
        mask_snow(is_lake) = 1
      else
        mask_snow(is_lake) = 0
      endif

    endif


    return

  end subroutine surface_hydrology

end module hydrology_mod


