!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module : w a t e r _ c h e c k _ m o d
!
!  Purpose : water (and water-isotope) conservation checks
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
module water_check_mod

   use precision, only : wp
   use lnd_grid, only : nveg, is_veg, is_ice, i_ice
   use lnd_params, only : dt, hydro_par
   use wiso_params, only : l_wiso, nwiso

   implicit none

   ! per-timestep closure tolerance, kg/m2; same value used for bulk and each iso species
   real(wp), parameter :: tol_water  = 1.d-10
   ! residual magnitude (kg/m2) above which the run is aborted — likely a real bug
   real(wp), parameter :: stop_water = 0.1_wp

   private
   public :: water_balance_check

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  w a t e r _ b a l a n c e _ c h e c k
  !   Purpose    :  per-timestep water and water-isotope conservation
  !                 check, per surface type. Lake budget is intentionally
  !                 not enforced (lake water is treated as inexhaustible).
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine water_balance_check(i, j, frac_surf, f_veg, &
                                 rain, snow, et, runoff_sur, calving, drainage, icemelt, icesub, &
                                 w_w, w_i, w_snow, w_can, s_can, &
                                 w_w_old, w_i_old, w_snow_old, w_can_old, s_can_old, &
                                 lake_water_tendency, &
                                 rain_iso, snow_iso, et_iso, &
                                 runoff_sur_iso, calving_iso, drainage_iso, icemelt_iso, icesub_iso, &
                                 w_w_iso, w_i_iso, w_snow_iso, w_can_iso, s_can_iso, &
                                 w_w_iso_old, w_i_iso_old, w_snow_iso_old, w_can_iso_old, s_can_iso_old, &
                                 water_cons, water_iso_cons)

    implicit none

    integer,  intent(in) :: i, j
    real(wp), intent(in) :: f_veg
    real(wp), intent(in) :: lake_water_tendency
    real(wp), dimension(:),   intent(in) :: frac_surf
    real(wp), dimension(:),   intent(in) :: rain, snow, et
    real(wp), dimension(:),   intent(in) :: runoff_sur, calving, drainage, icemelt, icesub
    real(wp), dimension(:),   intent(in) :: w_can, s_can, w_can_old, s_can_old
    real(wp), dimension(:),   intent(in) :: w_w, w_i, w_w_old, w_i_old
    real(wp), dimension(:),   intent(in) :: w_snow, w_snow_old
    real(wp), dimension(:,:), intent(in) :: rain_iso, snow_iso, et_iso
    real(wp), dimension(:,:), intent(in) :: runoff_sur_iso, calving_iso, drainage_iso, icemelt_iso, icesub_iso
    real(wp), dimension(:,:), intent(in) :: w_can_iso, s_can_iso, w_can_iso_old, s_can_iso_old
    real(wp), dimension(:,:), intent(in) :: w_w_iso, w_i_iso, w_w_iso_old, w_i_iso_old
    real(wp), dimension(:,:), intent(in) :: w_snow_iso, w_snow_iso_old
    real(wp), dimension(:),   intent(out) :: water_cons
    real(wp), dimension(:,:), intent(out) :: water_iso_cons

    integer  :: n, iso
    real(wp) :: rain_veg, snow_veg, et_veg, dw_can_veg
    real(wp) :: in_b, out_b, store_b
    real(wp) :: rain_veg_iso(nwiso), snow_veg_iso(nwiso), et_veg_iso(nwiso), dw_can_veg_iso(nwiso)
    real(wp) :: in_i, out_i, store_i

    water_cons     = 0._wp
    water_iso_cons = 0._wp

    !----------------------------------------------------------------
    ! vegetated grid part
    !----------------------------------------------------------------
    if (f_veg .gt. 0._wp) then

      ! bulk aggregation over veg PFTs
      rain_veg = 0._wp; snow_veg = 0._wp; et_veg = 0._wp; dw_can_veg = 0._wp
      do n = 1, nveg
        rain_veg   = rain_veg   + rain(n) * frac_surf(n)/f_veg
        snow_veg   = snow_veg   + snow(n) * frac_surf(n)/f_veg
        et_veg     = et_veg     + et(n)   * frac_surf(n)/f_veg
        dw_can_veg = dw_can_veg + (w_can(n) - w_can_old(n) + s_can(n) - s_can_old(n)) * frac_surf(n)/f_veg
      enddo

      in_b    = (rain_veg + snow_veg) * dt
      out_b   = (runoff_sur(is_veg) + calving(is_veg) + drainage(is_veg) + et_veg) * dt
      store_b = sum(w_w - w_w_old + w_i - w_i_old) + (w_snow(is_veg) - w_snow_old(is_veg)) + dw_can_veg
      water_cons(is_veg) = in_b - out_b - store_b

      if (abs(water_cons(is_veg)) .gt. tol_water) then
        call print_bulk_failure('is_veg', i, j, water_cons(is_veg), &
                                rain_veg*dt, snow_veg*dt, 0._wp, 0._wp, &
                                et_veg*dt, runoff_sur(is_veg)*dt, drainage(is_veg)*dt, calving(is_veg)*dt, &
                                sum(w_w-w_w_old), sum(w_i-w_i_old), w_snow(is_veg)-w_snow_old(is_veg), dw_can_veg)
        if (abs(water_cons(is_veg)) .gt. stop_water) stop 'water_balance_check: bulk imbalance over SOIL exceeds stop_water'
      endif

      ! iso, per species
      if (l_wiso) then
        rain_veg_iso(:)   = 0._wp
        snow_veg_iso(:)   = 0._wp
        et_veg_iso(:)     = 0._wp
        dw_can_veg_iso(:) = 0._wp
        do iso = 1, nwiso
          do n = 1, nveg
            rain_veg_iso(iso)   = rain_veg_iso(iso)   + rain_iso(n,iso) * frac_surf(n)/f_veg
            snow_veg_iso(iso)   = snow_veg_iso(iso)   + snow_iso(n,iso) * frac_surf(n)/f_veg
            et_veg_iso(iso)     = et_veg_iso(iso)     + et_iso(n,iso)   * frac_surf(n)/f_veg
            dw_can_veg_iso(iso) = dw_can_veg_iso(iso) &
                                + (w_can_iso(n,iso) - w_can_iso_old(n,iso) &
                                 + s_can_iso(n,iso) - s_can_iso_old(n,iso)) * frac_surf(n)/f_veg
          enddo

          in_i    = (rain_veg_iso(iso) + snow_veg_iso(iso)) * dt
          out_i   = (runoff_sur_iso(is_veg,iso) + calving_iso(is_veg,iso) &
                   + drainage_iso(is_veg,iso) + et_veg_iso(iso)) * dt
          store_i = sum(w_w_iso(:,iso) - w_w_iso_old(:,iso) + w_i_iso(:,iso) - w_i_iso_old(:,iso)) &
                  + (w_snow_iso(is_veg,iso) - w_snow_iso_old(is_veg,iso)) + dw_can_veg_iso(iso)
          water_iso_cons(is_veg,iso) = in_i - out_i - store_i

          if (abs(water_iso_cons(is_veg,iso)) .gt. tol_water) then
            call print_iso_failure('is_veg', i, j, iso, water_iso_cons(is_veg,iso), &
                                   rain_veg_iso(iso)*dt, snow_veg_iso(iso)*dt, 0._wp, 0._wp, &
                                   et_veg_iso(iso)*dt, runoff_sur_iso(is_veg,iso)*dt, &
                                   drainage_iso(is_veg,iso)*dt, calving_iso(is_veg,iso)*dt, &
                                   sum(w_w_iso(:,iso)-w_w_iso_old(:,iso)), &
                                   sum(w_i_iso(:,iso)-w_i_iso_old(:,iso)), &
                                   w_snow_iso(is_veg,iso)-w_snow_iso_old(is_veg,iso), dw_can_veg_iso(iso))
            if (abs(water_iso_cons(is_veg,iso)) .gt. stop_water) stop 'water_balance_check: iso imbalance over SOIL exceeds stop_water'
          endif
        enddo
      endif

    endif

    !----------------------------------------------------------------
    ! ice grid part
    !----------------------------------------------------------------
    if (frac_surf(i_ice) .gt. 0._wp) then

      ! icemelt only enters the lnd budget when routed into runoff
      in_b    = (rain(i_ice) + snow(i_ice) + icesub(is_ice)) * dt &
              + merge(icemelt(is_ice)*dt, 0._wp, hydro_par%l_runoff_icemelt)
      out_b   = (runoff_sur(is_ice) + calving(is_ice) + drainage(is_ice) + et(i_ice)) * dt
      store_b = w_snow(is_ice) - w_snow_old(is_ice)
      water_cons(is_ice) = in_b - out_b - store_b

      if (abs(water_cons(is_ice)) .gt. tol_water) then
        call print_bulk_failure('is_ice', i, j, water_cons(is_ice), &
                                rain(i_ice)*dt, snow(i_ice)*dt, &
                                merge(icemelt(is_ice)*dt, 0._wp, hydro_par%l_runoff_icemelt), icesub(is_ice)*dt, &
                                et(i_ice)*dt, runoff_sur(is_ice)*dt, drainage(is_ice)*dt, calving(is_ice)*dt, &
                                0._wp, 0._wp, w_snow(is_ice)-w_snow_old(is_ice), 0._wp)
        if (abs(water_cons(is_ice)) .gt. stop_water) stop 'water_balance_check: bulk imbalance over ICE exceeds stop_water'
      endif

      if (l_wiso) then
        do iso = 1, nwiso
          in_i    = (rain_iso(i_ice,iso) + snow_iso(i_ice,iso) + icesub_iso(is_ice,iso)) * dt &
                  + merge(icemelt_iso(is_ice,iso)*dt, 0._wp, hydro_par%l_runoff_icemelt)
          out_i   = (runoff_sur_iso(is_ice,iso) + calving_iso(is_ice,iso) &
                   + drainage_iso(is_ice,iso) + et_iso(i_ice,iso)) * dt
          store_i = w_snow_iso(is_ice,iso) - w_snow_iso_old(is_ice,iso)
          water_iso_cons(is_ice,iso) = in_i - out_i - store_i

          if (abs(water_iso_cons(is_ice,iso)) .gt. tol_water) then
            call print_iso_failure('is_ice', i, j, iso, water_iso_cons(is_ice,iso), &
                                   rain_iso(i_ice,iso)*dt, snow_iso(i_ice,iso)*dt, &
                                   merge(icemelt_iso(is_ice,iso)*dt, 0._wp, hydro_par%l_runoff_icemelt), &
                                   icesub_iso(is_ice,iso)*dt, &
                                   et_iso(i_ice,iso)*dt, runoff_sur_iso(is_ice,iso)*dt, &
                                   drainage_iso(is_ice,iso)*dt, calving_iso(is_ice,iso)*dt, &
                                   0._wp, 0._wp, w_snow_iso(is_ice,iso)-w_snow_iso_old(is_ice,iso), 0._wp)
            if (abs(water_iso_cons(is_ice,iso)) .gt. stop_water) stop 'water_balance_check: iso imbalance over ICE exceeds stop_water'
          endif
        enddo
      endif

    endif

    ! Lake water is treated as inexhaustible in this model; no budget enforcement here.
    ! lake_water_tendency is accepted as an argument for diagnostic purposes only.

    return

  end subroutine water_balance_check


  ! ---- internal pretty-printers ----------------------------------------

  subroutine print_bulk_failure(tag, i, j, residual, &
                                rain_dt, snow_dt, icemelt_dt, icesub_dt, &
                                et_dt, runoff_dt, drainage_dt, calving_dt, &
                                dw_w, dw_i, dw_snow, dw_can)
    character(len=*), intent(in) :: tag
    integer,  intent(in) :: i, j
    real(wp), intent(in) :: residual
    real(wp), intent(in) :: rain_dt, snow_dt, icemelt_dt, icesub_dt
    real(wp), intent(in) :: et_dt, runoff_dt, drainage_dt, calving_dt
    real(wp), intent(in) :: dw_w, dw_i, dw_snow, dw_can
    real(wp) :: in_tot, out_tot, store_tot

    in_tot    = rain_dt + snow_dt + icemelt_dt + icesub_dt
    out_tot   = et_dt + runoff_dt + drainage_dt + calving_dt
    store_tot = dw_w + dw_i + dw_snow + dw_can

    print '(A,A,A,I0,A,I0,A,ES12.4,A,ES10.2,A)', &
      'WATER IMBALANCE ',tag,' at (i=',i,',j=',j,') = ',residual,' kg/m2  (tol=',tol_water,')'
    print '(A,4(A,ES10.3),A,ES10.3)', &
      '  IN:   rain=',rain_dt,'  snow=',snow_dt,'  icemelt=',icemelt_dt,'  icesub=',icesub_dt,'  -> ',in_tot
    print '(A,4(A,ES10.3),A,ES10.3)', &
      '  OUT:  et=',et_dt,'  runoff=',runoff_dt,'  drainage=',drainage_dt,'  calving=',calving_dt,'  -> ',out_tot
    print '(A,4(A,ES10.3),A,ES10.3)', &
      '  STORE:dw_w=',dw_w,'  dw_i=',dw_i,'  dw_snow=',dw_snow,'  dw_can=',dw_can,'  -> ',store_tot
  end subroutine print_bulk_failure

  subroutine print_iso_failure(tag, i, j, iso, residual, &
                               rain_dt, snow_dt, icemelt_dt, icesub_dt, &
                               et_dt, runoff_dt, drainage_dt, calving_dt, &
                               dw_w, dw_i, dw_snow, dw_can)
    character(len=*), intent(in) :: tag
    integer,  intent(in) :: i, j, iso
    real(wp), intent(in) :: residual
    real(wp), intent(in) :: rain_dt, snow_dt, icemelt_dt, icesub_dt
    real(wp), intent(in) :: et_dt, runoff_dt, drainage_dt, calving_dt
    real(wp), intent(in) :: dw_w, dw_i, dw_snow, dw_can
    real(wp) :: in_tot, out_tot, store_tot

    in_tot    = rain_dt + snow_dt + icemelt_dt + icesub_dt
    out_tot   = et_dt + runoff_dt + drainage_dt + calving_dt
    store_tot = dw_w + dw_i + dw_snow + dw_can

    print '(A,A,A,I0,A,I0,A,I0,A,ES12.4,A,ES10.2,A)', &
      'WATER-ISO IMBALANCE ',tag,' at (i=',i,',j=',j,') iso=',iso,' : ',residual,' kg/m2  (tol=',tol_water,')'
    print '(A,4(A,ES10.3),A,ES10.3)', &
      '  IN:   rain=',rain_dt,'  snow=',snow_dt,'  icemelt=',icemelt_dt,'  icesub=',icesub_dt,'  -> ',in_tot
    print '(A,4(A,ES10.3),A,ES10.3)', &
      '  OUT:  et=',et_dt,'  runoff=',runoff_dt,'  drainage=',drainage_dt,'  calving=',calving_dt,'  -> ',out_tot
    print '(A,4(A,ES10.3),A,ES10.3)', &
      '  STORE:dw_w=',dw_w,'  dw_i=',dw_i,'  dw_snow=',dw_snow,'  dw_can=',dw_can,'  -> ',store_tot
  end subroutine print_iso_failure

end module water_check_mod
