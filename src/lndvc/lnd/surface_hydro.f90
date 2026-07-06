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
module lndvc_hydrology_mod

  use precision, only : wp
  use timer, only : sec_day
  use constants, only : q_sat_w, q_sat_i, Lf, rho_a, g, T0
  use control, only : check_water
  use lnd_grid, only : nsurf, npft, nveg, nsoil, i_ice, i_lake, is_veg, is_ice, is_lake
  use lnd_grid, only : flag_veg, flag_pft, flag_tree
  use lnd_grid, only : z_int, dz, nl
  use lnd_params, only : dt, rdt
  use lnd_params, only : snow_par, hydro_par
  use wiso_params, only : l_wiso, nwiso, i_o18, Rstd

  implicit none

  private
  public :: canopy_water, surface_hydrology_lake

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  c a n o p y _ w a t e r
  !   Purpose    :  solve analytical model of canopy water 
  !              :  including interception and throughfall
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine canopy_water(frac_surf,lai,sai,r_a,t_skin,pressure,qair,rain,snow, &
                         w_can,w_can_old,s_can,s_can_old, &
                         rain_ground,snow_ground,evap_can,subl_can,f_wat_can,f_snow_can, &
                         rain_iso,snow_iso, &
                         w_can_iso,w_can_iso_old,s_can_iso,s_can_iso_old, &
                         rain_ground_iso,snow_ground_iso,evap_can_iso,subl_can_iso)

    implicit none

    real(wp), dimension(:), intent(in) :: frac_surf, r_a, t_skin, qair
    real(wp), dimension(:), intent(in) :: lai, sai
    real(wp), dimension(:), intent(in) :: pressure
    real(wp), dimension(:), intent(in) :: rain, snow
    real(wp), dimension(:), intent(inout) :: w_can, w_can_old, s_can, s_can_old
    real(wp), dimension(:), intent(inout) :: rain_ground, snow_ground, evap_can, subl_can
    real(wp), dimension(:), intent(inout) :: f_wat_can, f_snow_can
    ! water-isotope arguments (always passed; values are 0 unless l_wiso=.true.)
    real(wp), dimension(:,:), intent(in)    :: rain_iso, snow_iso
    real(wp), dimension(:,:), intent(inout) :: w_can_iso, w_can_iso_old
    real(wp), dimension(:,:), intent(inout) :: s_can_iso, s_can_iso_old
    real(wp), dimension(:,:), intent(inout) :: rain_ground_iso, snow_ground_iso
    real(wp), dimension(:,:), intent(inout) :: evap_can_iso, subl_can_iso

    integer :: n, iso
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
        if (l_wiso) then
          do iso=1,nwiso
            w_can_iso_old(n,iso) = w_can_iso(n,iso)
            s_can_iso_old(n,iso) = s_can_iso(n,iso)
          enddo
        endif

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

          ! water isotopes: same implicit equation, same denominator
          if (l_wiso) then
            do iso=1,nwiso
              w_can_iso(n,iso) = (fac_i_w*rain_iso(n,iso) + w_can_iso(n,iso)*rdt) &
                                / (1._wp*rdt + fac_e_w/w_can_max + 1._wp/hydro_par%tau_w)
              if (w_can_iso(n,iso) .lt. 0._wp)              w_can_iso(n,iso) = 0._wp
              if (w_can_iso(n,iso) .lt. 1.e-30_wp*Rstd(iso)) w_can_iso(n,iso) = 0._wp
              if (w_can_iso(n,iso) .gt. w_can_max*Rstd(iso)) w_can_iso(n,iso) = w_can_max*Rstd(iso)
              evap_can_iso(n,iso)    = fac_e_w * w_can_iso(n,iso)/w_can_max
              rain_ground_iso(n,iso) = rain_iso(n,iso) - evap_can_iso(n,iso) &
                                      - (w_can_iso(n,iso) - w_can_iso_old(n,iso)) * rdt
            enddo
          endif

        else

          rain_ground(n) = rain(n)
          evap_can(n) = 0._wp
          f_wat_can(n) = 0._wp
          if (l_wiso) then
            do iso=1,nwiso
              rain_ground_iso(n,iso) = rain_iso(n,iso)
              evap_can_iso(n,iso)    = 0._wp
            enddo
          endif

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

          if (l_wiso) then
            do iso=1,nwiso
              s_can_iso(n,iso) = (fac_i_s*snow_iso(n,iso) + s_can_iso(n,iso)*rdt) &
                                / (1._wp*rdt + fac_e_s/s_can_max + 1._wp/tau_s)
              if (s_can_iso(n,iso) .lt. 1.e-20_wp*Rstd(iso)) s_can_iso(n,iso) = 0._wp
              if (s_can_iso(n,iso) .gt. s_can_max*Rstd(iso)) s_can_iso(n,iso) = s_can_max*Rstd(iso)
              subl_can_iso(n,iso)    = fac_e_s * s_can_iso(n,iso)/s_can_max
              snow_ground_iso(n,iso) = snow_iso(n,iso) - subl_can_iso(n,iso) &
                                      - (s_can_iso(n,iso) - s_can_iso_old(n,iso)) * rdt
            enddo
          endif

        else

          snow_ground(n) = snow(n)
          subl_can(n) = 0._wp
          f_snow_can(n) = 0._wp
          if (l_wiso) then
            do iso=1,nwiso
              snow_ground_iso(n,iso) = snow_iso(n,iso)
              subl_can_iso(n,iso)    = 0._wp
            enddo
          endif

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

        if (l_wiso) then
          do iso=1,nwiso
            if (rain_ground_iso(n,iso) .lt. 0._wp) rain_ground_iso(n,iso) = 0._wp
            if (snow_ground_iso(n,iso) .lt. 0._wp) snow_ground_iso(n,iso) = 0._wp
          enddo
        endif


      else ! PFT not present

        evap_can(n) = 0._wp
        subl_can(n) = 0._wp
        w_can(n)    = 0._wp
        s_can(n)    = 0._wp
        f_wat_can(n) = 0._wp
        f_snow_can(n) = 0._wp
        if (l_wiso) then
          do iso=1,nwiso
            evap_can_iso(n,iso) = 0._wp
            subl_can_iso(n,iso) = 0._wp
            w_can_iso(n,iso)    = 0._wp
            s_can_iso(n,iso)    = 0._wp
          enddo
        endif

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
    if (l_wiso) then
      do iso=1,nwiso
        where (flag_pft.eq.0)
          rain_ground_iso(:,iso) = rain_iso(:,iso)
          snow_ground_iso(:,iso) = snow_iso(:,iso)
          evap_can_iso(:,iso) = 0._wp
          subl_can_iso(:,iso) = 0._wp
        endwhere
      enddo
    endif

    return

  end subroutine canopy_water



  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s u r f a c e _ h y d r o l o g y _ l a k e
  !   Purpose    :  surface hydrology for a single lake tile
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine surface_hydrology_lake(mask_snow,evap_surface,snow_ground,rain_ground,snowmelt, &
                                   cap_lake,t_lake1,w_snow_old,w_snow,w_snow_max,h_snow, &
                                   calving,runoff_sur,lake_water_tendency, &
                                   evap_surface_iso,snow_ground_iso,rain_ground_iso, &
                                   snowmelt_iso,w_snow_iso,calving_iso,runoff_sur_iso)

    implicit none

    integer,  intent(inout) :: mask_snow
    real(wp), intent(in)    :: evap_surface, snow_ground, rain_ground
    real(wp), intent(inout) :: snowmelt
    real(wp), intent(in)    :: cap_lake
    real(wp), intent(inout) :: t_lake1
    real(wp), intent(in)    :: w_snow_old
    real(wp), intent(inout) :: w_snow, w_snow_max
    real(wp), intent(out)   :: h_snow
    real(wp), intent(out)   :: calving, runoff_sur
    real(wp), intent(out)   :: lake_water_tendency
    ! water-isotope siblings (always passed; values are 0 unless l_wiso=.true.)
    real(wp), dimension(:), intent(in)    :: evap_surface_iso, snow_ground_iso, rain_ground_iso
    real(wp), dimension(:), intent(inout) :: snowmelt_iso, w_snow_iso
    real(wp), dimension(:), intent(out)   :: calving_iso, runoff_sur_iso

    integer :: iso
    real(wp) :: sublimation, evaporation, w_snow_pre, wsnowold_m
    real(wp) :: H, H_m, H_star
    real(wp) :: sublimation_iso(nwiso), evaporation_iso(nwiso), w_snow_iso_pre(nwiso)
    real(wp) :: wsnowold_iso(nwiso), ratio_snow(nwiso)


    calving = 0._wp
    runoff_sur = 0._wp
    lake_water_tendency = 0._wp
    if (l_wiso) then
      calving_iso    = 0._wp
      runoff_sur_iso = 0._wp
    endif

    if( mask_snow .eq. 1 ) then
      sublimation = evap_surface
      evaporation = 0._wp
      if (l_wiso) then
        sublimation_iso(:) = evap_surface_iso(:)
        evaporation_iso(:) = 0._wp
      endif
    else
      sublimation = 0._wp
      evaporation = evap_surface
      if (l_wiso) then
        sublimation_iso(:) = 0._wp
        evaporation_iso(:) = evap_surface_iso(:)
      endif
    endif

    ! add snowfall to snow layer and remove sublimation, snowmelt already removed during lake temperature update
    w_snow_pre = w_snow + snow_ground * dt - sublimation * dt  ! kg/m2
    if (l_wiso) then
      do iso=1,nwiso
        w_snow_iso_pre(iso) = w_snow_iso(iso) &
                            + snow_ground_iso(iso) * dt &
                            - sublimation_iso(iso) * dt
      enddo
    endif

    ! if sublimation depleted the snowpack, redirect the deficit to lake-water evaporation
    ! (lake water is treated as inexhaustible in this model)
    if (w_snow_pre .lt. 0._wp) then
      evaporation = evaporation + (-w_snow_pre * rdt)   ! kg/m2/s
      if (l_wiso) then
        do iso=1,nwiso
          evaporation_iso(iso) = evaporation_iso(iso) + (-w_snow_iso_pre(iso) * rdt)
        enddo
      endif
    endif

    w_snow = max(0._wp, w_snow_pre)
    if (l_wiso) then
      do iso=1,nwiso
        w_snow_iso(iso) = max(0._wp, w_snow_iso_pre(iso))
      enddo
    endif

    ! if snowmass below critical snow mass for explicit snow layer, use possible first lake layer excess energy to melt snow
    ! and update snowmelt
    if (w_snow.gt.0._wp .and. w_snow.lt.snow_par%w_snow_crit .and. t_lake1.gt.T0) then
      H = cap_lake*rdt*(t_lake1 - T0)    ! W/m2, energy available to melt snow
      wsnowold_m = w_snow
      if (l_wiso) wsnowold_iso = w_snow_iso(:)
      H_m = H*dt/Lf ! kg/m2, snow that can be melted
      ! update w_snow
      w_snow = max( 0._wp, w_snow-H_m )  ! kg/m2
      H_star = H - Lf*rdt * (wsnowold_m - w_snow) ! heat not used to melt snow
      t_lake1 = T0 + dt/cap_lake * H_star
      ! update snowmelt
      snowmelt = snowmelt + (wsnowold_m - w_snow) * rdt  ! kg/m2/s
      if (l_wiso .and. wsnowold_m.gt.0._wp) then
        do iso=1,nwiso
          ratio_snow(iso) = wsnowold_iso(iso) / wsnowold_m
          w_snow_iso(iso)   = ratio_snow(iso) * w_snow
          snowmelt_iso(iso) = snowmelt_iso(iso) &
                                     + ratio_snow(iso) * (wsnowold_m - w_snow) * rdt
        enddo
      endif
    endif

    ! limit w_snow and add to 'calving'
    if( w_snow .gt. snow_par%w_snow_off ) then
      calving = (w_snow - snow_par%w_snow_off) * rdt ! kg/m2/s
      if (l_wiso .and. w_snow.gt.0._wp) then
        do iso=1,nwiso
          ratio_snow(iso) = w_snow_iso(iso) / w_snow
          calving_iso(iso) = ratio_snow(iso) * (w_snow - snow_par%w_snow_off) * rdt
          w_snow_iso(iso)  = ratio_snow(iso) * snow_par%w_snow_off
        enddo
      endif
      w_snow = snow_par%w_snow_off
    endif
    ! save seasonal maximum snow swe
    if (w_snow.gt.w_snow_old) w_snow_max = w_snow

    ! update snow height
    h_snow = w_snow / snow_par%rho_snow

    ! net lake surface water balance
    lake_water_tendency = rain_ground + snowmelt - evaporation ! kg/m2/s

    ! route the lake surface water balance (P + M - E) to runoff_sur;
    ! the coupler transfers it to the ocean (when interactive lakes are off)
    ! evaporation is included as negative runoff in order to close the lake water budget
    runoff_sur = lake_water_tendency
    if (l_wiso) then
      do iso=1,nwiso
        runoff_sur_iso(iso) = rain_ground_iso(iso) + snowmelt_iso(iso) - evaporation_iso(iso)
      enddo
    endif

    ! update snow mask
    if( w_snow .gt. snow_par%w_snow_crit ) then
      mask_snow = 1
    else
      mask_snow = 0
    endif

    return

  end subroutine surface_hydrology_lake

end module lndvc_hydrology_mod


