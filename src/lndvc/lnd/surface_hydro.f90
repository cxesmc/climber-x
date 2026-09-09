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
  use lnd_params, only : dt, rdt
  use lnd_params, only : snow_par, hydro_par
  use wiso_params, only : l_wiso, nwiso, i_o18, Rstd

  implicit none

  private
  public :: canopy_water, surface_hydrology_lake, surface_hydrology_veg

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


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s u r f a c e _ h y d r o l o g y _ v e g
  !   Purpose    :  vegetated-land surface hydrology (port C.2a):
  !              :  snow-layer update, water-table depth, wetland (TOPMODEL/
  !              :  DYPTOP) fraction, surface runoff and infiltration.
  !   Scope      :  the veg branch of the reference all-class surface_hydrology,
  !              :  scalarised to the vc's single shared soil/snow column. Its
  !              :  bare+PFT sub-tiles remain multi-tile arrays (the frac-
  !              :  weighted sublimation / snow / rain sums over the vc's own
  !              :  tiles). The ice branch is redundant (SEMI owns ice-tile snow)
  !              :  and the lake branch lives in surface_hydrology_lake, so both
  !              :  are dropped along with their exclusive arguments
  !              :  (icemelt/icesub/et, cap_lake/t_lake, lake_water_tendency).
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine surface_hydrology_veg(frac_surf,mask_snow,evap_surface,rain_ground,snow_ground,snowmelt, &
                                  k_sat,cap_soil, &
                                  cti_mean, cti_cdf, &
                                  w_snow_old,w_snow,w_snow_max,w_w,w_i,w_table_cum,f_wet_cum,t_soil, &
                                  h_snow,calving,runoff_sur,infiltration,w_table_eff,fz_eff,f_wet,f_wet_max,cti_lim, &
                                  evap_surface_iso,rain_ground_iso,snow_ground_iso,snowmelt_iso, &
                                  w_snow_iso,w_w_iso,w_i_iso, &
                                  calving_iso,runoff_sur_iso,infiltration_iso)

    implicit none

    integer, intent(in) :: mask_snow
    real(wp), dimension(:), intent(in) :: frac_surf, snow_ground, evap_surface, rain_ground
    real(wp), dimension(:), intent(in) :: k_sat
    real(wp), intent(in) :: cap_soil
    real(wp), intent(in) :: cti_mean, cti_cdf(:)
    real(wp), intent(in) :: w_snow_old
    real(wp), intent(inout) :: snowmelt
    real(wp), intent(inout) :: w_snow, w_snow_max
    real(wp), intent(inout) :: w_table_cum, f_wet_cum
    real(wp), dimension(0:), intent(inout) :: t_soil
    real(wp), dimension(:), intent(inout) :: w_w, w_i
    real(wp), intent(out) :: h_snow, calving, runoff_sur
    real(wp), intent(out) :: infiltration, f_wet, f_wet_max, cti_lim
    real(wp), intent(in) :: w_table_eff  ! m, effective water table: the shallower of the aquifer
                                         ! table and the table perched on the frost table
    real(wp), intent(in) :: fz_eff       ! -, f*z of whichever of the two regimes is the shallower,
                                         ! i.e. the TOPMODEL shift of the CTI threshold
    ! water-isotope siblings (always passed; values are 0 unless l_wiso=.true.)
    real(wp), dimension(:,:), intent(in)    :: evap_surface_iso, snow_ground_iso, rain_ground_iso
    real(wp), dimension(:),   intent(inout) :: snowmelt_iso, w_snow_iso
    real(wp), dimension(:,:), intent(inout) :: w_w_iso, w_i_iso
    real(wp), dimension(:),   intent(out)   :: calving_iso, runoff_sur_iso, infiltration_iso

    integer :: n, iso
    real(wp) :: f_sat, rain_g, snow_g, H, H_m, H_star, wsnowold_m
    real(wp) :: f_veg
    real(wp) :: infiltration_max, f_run_sur, q_liq, sublimation, dws
    real(wp) :: snow_g_iso(nwiso), rain_g_iso(nwiso), sublimation_iso(nwiso)
    real(wp) :: q_liq_iso(nwiso), ratio_snow(nwiso), wsnowold_iso(nwiso)


    runoff_sur = 0._wp
    calving = 0._wp
    infiltration = 0._wp
    f_wet = 0._wp
    f_wet_max = 0._wp
    cti_lim = 0._wp
    if (l_wiso) then
      runoff_sur_iso   = 0._wp
      calving_iso      = 0._wp
      infiltration_iso = 0._wp
    endif

    f_veg = sum(frac_surf,mask=flag_veg.eq.1)

    !************************
    ! vegetated grid part
    !************************
    if( f_veg .gt. 0._wp ) then

      sublimation = 0._wp
      if (l_wiso) sublimation_iso = 0._wp
      if( mask_snow .eq. 1 ) then
        do n=1,nveg
          ! mean sublimation from snow
          if (frac_surf(n).gt.0._wp) then
            sublimation = sublimation + evap_surface(n)*frac_surf(n)/f_veg  ! kg/m2/s
            if (l_wiso) then
              do iso=1,nwiso
                sublimation_iso(iso) = sublimation_iso(iso) + evap_surface_iso(n,iso)*frac_surf(n)/f_veg
              enddo
            endif
          endif
        enddo
      endif

      ! mean snow on the ground over vegetated part
      snow_g = 0._wp
      if (l_wiso) snow_g_iso = 0._wp
      do n=1,nsurf
        if (flag_veg(n).eq.1 .and. frac_surf(n).gt.0._wp) then
          snow_g = snow_g + snow_ground(n) * frac_surf(n)/f_veg
          if (l_wiso) then
            do iso=1,nwiso
              snow_g_iso(iso) = snow_g_iso(iso) + snow_ground_iso(n,iso) * frac_surf(n)/f_veg
            enddo
          endif
        endif
      enddo

      ! add snowfall to snow layer and remove sublimation, snowmelt already removed during soil temperature update
      w_snow = w_snow + snow_g * dt - sublimation * dt  ! kg/m2
      if (l_wiso) then
        do iso=1,nwiso
          w_snow_iso(iso) = w_snow_iso(iso) + snow_g_iso(iso) * dt - sublimation_iso(iso) * dt
        enddo
      endif

      ! if snowmass below critical snow mass for explicit snow layer, use possible first soil layer excess energy to melt snow
      ! and update snowmelt
      if (w_snow.gt.0._wp .and. w_snow.lt.snow_par%w_snow_crit .and. t_soil(1).gt.T0) then
        H = cap_soil*rdt*(t_soil(1) - T0)    ! W/m2, energy available to melt snow
        wsnowold_m = w_snow
        if (l_wiso) wsnowold_iso = w_snow_iso(:)
        H_m = H*dt/Lf ! kg/m2, snow that can be melted
        ! update w_snow
        w_snow = max( 0._wp, w_snow-H_m )  ! kg/m2
        H_star = H - Lf*rdt * (wsnowold_m - w_snow) ! heat not used to melt snow
        t_soil(1) = T0 + dt/cap_soil * H_star
        ! update snowmelt
        snowmelt = snowmelt + (wsnowold_m - w_snow) * rdt  ! kg/m2/s
        if (l_wiso .and. wsnowold_m.gt.0._wp) then
          do iso=1,nwiso
            ! melted mass at current snowpack ratio
            ratio_snow(iso) = wsnowold_iso(iso) / wsnowold_m
            w_snow_iso(iso) = ratio_snow(iso) * w_snow
            snowmelt_iso(iso) = snowmelt_iso(iso) &
                                      + ratio_snow(iso) * (wsnowold_m - w_snow) * rdt
          enddo
        endif
      endif

      ! if w_snow negative, reset to 0 and remove required ice from the top soil layer (sublimation)
      if( w_snow .lt. 0._wp ) then
        ! if not enough ice and water in first layer, remove also ice from second layer
        if(-w_snow .gt. (w_i(1)+w_w(1))) then
          if (check_water) print *,'WARNING: not enough ice or water to sublimate in top layer!', w_snow,w_i(1),w_w(1)
          dws = -(w_snow+w_i(1)+w_w(1))
          if (l_wiso) then
            do iso=1,nwiso
              if (w_i(2).gt.0._wp) then
                w_i_iso(2,iso) = w_i_iso(2,iso) - (w_i_iso(2,iso)/w_i(2)) * min(w_i(2),dws)
              endif
              w_i_iso(1,iso) = 0._wp
              w_w_iso(1,iso) = 0._wp
            enddo
          endif
          w_i(2) = w_i(2) - min(w_i(2),dws)
          w_i(1) = 0._wp
          w_w(1) = 0._wp
          ! if not enough ice in first layer, remove liquid water instead to keep water balance
        elseif(-w_snow .gt. w_i(1)) then
          if (check_water) print *,'WARNING: not enough ice to sublimate in top layer, sublimate also liquid water!', w_snow,w_i(1)
          dws = -(w_snow+w_i(1))
          if (l_wiso) then
            do iso=1,nwiso
              if (w_w(1).gt.0._wp) then
                w_w_iso(1,iso) = w_w_iso(1,iso) - (w_w_iso(1,iso)/w_w(1)) * min(w_w(1),dws)
              endif
              w_i_iso(1,iso) = 0._wp
            enddo
          endif
          w_w(1) = w_w(1) - min(w_w(1),dws)
          w_i(1) = 0._wp
        else
          ! enough ice to sublimate in top layer
          if (l_wiso .and. w_i(1).gt.0._wp) then
            do iso=1,nwiso
              w_i_iso(1,iso) = w_i_iso(1,iso) + (w_i_iso(1,iso)/w_i(1)) * w_snow
            enddo
          endif
          w_i(1) = w_i(1) + w_snow
        endif
        w_snow = 0._wp
        if (l_wiso) w_snow_iso(:) = 0._wp
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
      ! update snow height
      h_snow = w_snow / snow_par%rho_snow

      ! save seasonal maximum snow swe
      if (w_snow.gt.w_snow_old) w_snow_max = w_snow

      ! the water table is prognostic and updated in groundwater() after the soil hydrology,
      ! so it is used here with a one time step lag, exactly as in soil_hydro.
      ! w_table_eff, not the aquifer table, is what wets the surface: where the ground freezes the
      ! saturated zone that matters sits on the frost table, not tens of metres down in the aquifer.
      ! w_table_cum feeds w_table_mon and hence the peat acrotelm oxic fraction, so that sees it too.

      w_table_cum = w_table_cum + w_table_eff

      ! max possible wetland extent (w_table=0)
      f_wet_max = f_cti_exceed(max(cti_mean,hydro_par%cti_min), cti_cdf)
      ! no inundation if CTI lower than critical value (5.5 in Kleinen 2020)
      if (cti_mean.le.hydro_par%cti_mean_crit) f_wet_max = 0._wp

      ! saturated grid cell fraction, TOPMODEL following Kleinen et al 2020.
      ! f_drain is the TOPMODEL transmissivity decay factor, the same parameter that sets the
      ! aquifer baseflow recession in groundwater().
      ! The SIMTOP (Niu 2005) and DYPTOP (Stocker 2014) alternatives were removed: both were
      ! calibrated for a water table within a metre or two of the surface, and with the
      ! prognostic aquifer, which puts it near 13 m, both collapse to essentially no wetland
      ! (0.8 and 0.008 mln km2 against ~3.1 observed), DYPTOP additionally overflowing exp().
      ! Neither has a free parameter left to recalibrate with, since f_wtab was merged into
      ! f_drain and the DYPTOP shape parameters are read from a file fitted elsewhere.
      cti_lim = cti_mean + fz_eff   ! f_drain*w_table, or f_drain_perch*w_table_perch
      cti_lim = max(cti_lim,hydro_par%cti_min)
      f_sat = f_cti_exceed(cti_lim, cti_cdf)
      ! no inundation if CTI lower than critical value (5.5 in Kleinen 2020)
      if (cti_mean.le.hydro_par%cti_mean_crit) f_sat = 0._wp
      if( mask_snow .eq. 1 ) then
        f_wet = 0._wp ! no wetland where snow
      else
        f_wet = f_sat
      endif

      f_wet_cum = f_wet_cum + f_wet

      ! Infiltration capacity, kg/m2/s. The relevant conductivity is the one at SATURATION of
      ! the pore space open to liquid water, not the one at the antecedent water content: during
      ! an event the top layer wets up towards saturation and the Green-Ampt capacity decays to
      ! k_sat from above, never below it.
      infiltration_max = k_sat(1) * (1._wp-f_sat)

      ! mean rain on the ground over vegetated part
      rain_g = 0._wp
      if (l_wiso) rain_g_iso = 0._wp
      do n=1,nsurf
        if (flag_veg(n).eq.1 .and. frac_surf(n).gt.0._wp) then
          rain_g = rain_g + rain_ground(n) * frac_surf(n)/f_veg
          if (l_wiso) then
            do iso=1,nwiso
              rain_g_iso(iso) = rain_g_iso(iso) + rain_ground_iso(n,iso) * frac_surf(n)/f_veg
            enddo
          endif
        endif
      enddo
      q_liq = rain_g + snowmelt ! kg/m2/s

      ! surface runoff, kg/m2/s
      runoff_sur = f_sat * q_liq &  ! all into runoff over the saturated fraction
        + (1._wp - f_sat) * max(0._wp, q_liq - infiltration_max)    ! infiltration excess
      ! NOTE the infiltration-excess term is unreachable: k_sat(1) is ~100 mm/day while q_liq is
      ! a DAILY MEAN at 5 degrees and never approaches that. Activating it needs a sub-daily/
      ! sub-grid rainfall intensity distribution, not a different conductivity.
      ! soil liquid water infiltration, kg/m2/s
      infiltration = rain_g + snowmelt - runoff_sur  ! kg/m2/s

      ! iso: surface runoff removes water at the isotopic ratio of the incoming flux, so the
      ! bulk runoff fraction applies unchanged to the iso fluxes (no fractionation).
      ! NOTE the split must not be re-derived from the iso fluxes themselves: applying the
      ! infiltration_max threshold to q_liq_iso separately does not preserve the incoming ratio
      if (l_wiso) then
        if (q_liq .gt. 0._wp) then
          f_run_sur = runoff_sur / q_liq
        else
          f_run_sur = 0._wp
        endif
        do iso=1,nwiso
          q_liq_iso(iso) = rain_g_iso(iso) + snowmelt_iso(iso)
          runoff_sur_iso(iso) = f_run_sur * q_liq_iso(iso)
          infiltration_iso(iso) = q_liq_iso(iso) - runoff_sur_iso(iso)
        enddo
      endif

    endif

    return

  end subroutine surface_hydrology_veg


  ! ------------------------------------------------------------------------------------------
  ! Fraction of the grid cell with a compound topographic index above cti, i.e. the exceedance
  ! probability 1-cdf(cti), from the tabulated CTI distribution (Marthews et al. 2015, one
  ! entry per integer CTI bin).
  !
  ! The exceedance probability decays close to geometrically with CTI (area-weighted global
  ! means: 0.105, 0.067, 0.042, 0.025, 0.0128, 0.0050 for bins 9..14, a ratio of about 0.6 per
  ! unit), so it is interpolated log-linearly. Interpolating the cdf linearly instead, as was
  ! done before, overshoots systematically between the nodes because the tail is convex, and
  ! inflated the global wetland area by about 4% (+0.2 mln km2). The last bin is exactly zero
  ! because the distribution ends there, which log-linear interpolation cannot represent, so
  ! that one interval falls back to linear.
  ! ------------------------------------------------------------------------------------------
  function f_cti_exceed(cti, cti_cdf) result(f)

    implicit none

    real(wp), intent(in) :: cti
    real(wp), intent(in) :: cti_cdf(:)
    real(wp) :: f

    integer :: i1, i2, ncti
    real(wp) :: w2, e1, e2

    ncti = size(cti_cdf)

    if (cti .ge. real(ncti,wp)) then
      ! above the top of the tabulated distribution, which the cdf reaches with a value of 1
      f = 0._wp
    else
      i1 = int(max(1._wp,cti))
      i2 = i1+1
      w2 = max(1._wp,cti)-real(i1,wp)
      e1 = 1._wp-cti_cdf(i1)
      e2 = 1._wp-cti_cdf(i2)
      if (e1.gt.0._wp .and. e2.gt.0._wp) then
        f = e1**(1._wp-w2) * e2**w2
      else
        f = (1._wp-w2)*e1 + w2*e2
      endif
    endif

  end function f_cti_exceed

end module lndvc_hydrology_mod


