!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module : g r o u n d w a t e r _ m o d
!
!  Purpose : unconfined aquifer below the soil column and prognostic water table depth
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
module lndvc_groundwater_mod

  ! Simple groundwater model following Niu et al. 2007 (SIMGM), as used in CLM4 and Noah-MP.
  ! An unconfined aquifer of thickness z_wtab_max and specific yield sy_aqf is placed below the
  ! resolved soil column. It is recharged by the drainage out of the soil column and discharged
  ! by a TOPMODEL baseflow that decays exponentially with the water table depth,
  !
  !   rho_w * sy_aqf * dz/dt = kappa_max * exp(-f_drain*z) - recharge  ,   z positive downward
  !
  ! This replaced a set of diagnostic parameterisations that mapped the column water content
  ! onto a depth and therefore scaled with the total soil depth z_int(nl). The equilibrium depth
  ! here is z = 1/f_drain * log(kappa_max/recharge), i.e. it is set by the climate and by two
  ! hydraulic parameters and is independent of the vertical grid.
  !
  ! The aquifer is in series below the soil column: the Richards bottom boundary condition in
  ! soil_hydro stays free drainage and that drainage is the recharge here.
  ! NOTE a variant that replaced the bottom boundary condition by the Darcy flux towards the
  ! water table, q(nl) = kappa_int(nl)*(1+psi(nl)/(w_table-z(nl))), was implemented and tested
  ! in a 5000 yr spinup and then removed. It makes no difference where the soil is wet, because
  ! there the recharge/baseflow balance sets the water table anyway. Where the soil is dry both
  ! fluxes vanish (recharge as kappa ~ theta^9, baseflow as exp(-f_drain*z)) and the water table
  ! drifts towards the hydrostatic level z = z(nl)+|psi(nl)| at a rate of ~0.2 m per 1000 years,
  ! so it never equilibrates and the answer depends on the length of the run.

  use precision, only : wp
  use timer, only : sec_day
  use constants, only : rho_w
  use lnd_grid, only : dz, z_int, nl
  use lnd_params, only : dt, rdt, hydro_par
  use wiso_params, only : l_wiso, nwiso

   implicit none

   private
   public :: groundwater, w_table_ini, perched_water_table

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  p e r c h e d _ w a t e r _ t a b l e
  !   Purpose    :  diagnose a water table perched on the frost table, and the effective
  !              :  water table seen by the wetland and peat schemes
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Where the soil freezes, the saturated zone that matters for surface wetness sits in the
  ! active layer on top of the frost table, not in the aquifer tens of metres below. SIMGM has
  ! no such zone, so without this the model has no boreal wetlands at all once liquid water is
  ! correctly prevented from passing through frozen ground into the aquifer (see soil_hydro).
  !
  ! Purely diagnostic: everything is recomputed each step from w_w, theta_i and theta_sat, all
  ! of which are already prognostic and already in the restart. No new restart field.
  subroutine perched_water_table(theta_sat,theta_i,w_w,w_table, w_table_perch,w_table_eff,fz_eff)

    implicit none

    real(wp), dimension(:), intent(in) :: theta_sat, theta_i, w_w
    real(wp), intent(in) :: w_table                  ! m, aquifer water table, positive down
    real(wp), intent(out) :: w_table_perch           ! m, perched water table, z_wtab_max if no frost table
    real(wp), intent(out) :: w_table_eff             ! m, what the peat scheme should use
    real(wp), intent(out) :: fz_eff                  ! -, TOPMODEL shift of cti_lim, f*z of the winning regime

    integer :: k, k_frz
    real(wp) :: deficit, w_w_max_k

    ! frost table = shallowest layer whose ice fraction exceeds f_ice_perch
    k_frz = 0
    do k=1,nl
      if (theta_i(k)/max(1.e-10_wp,theta_sat(k)) .gt. hydro_par%f_ice_perch) then
        k_frz = k
        exit
      endif
    enddo

    if (k_frz .gt. 1) then
      ! walk up from the frost table while the layers are saturated with liquid water
      w_table_perch = z_int(k_frz-1)
      do k=k_frz-1,1,-1
        w_w_max_k = (theta_sat(k)-theta_i(k))*dz(k)*rho_w      ! ice-reduced liquid capacity
        deficit = w_w_max_k - w_w(k)                           ! kg/m2
        if (deficit .le. 0._wp) then
          w_table_perch = z_int(k-1)                           ! saturated, table at or above the top
        else
          ! the table lies inside layer k: place it by the deficit over the effective porosity,
          ! so that it stays continuous within the layer instead of snapping to layer interfaces
          w_table_perch = z_int(k-1) &
                        + min(dz(k), deficit/(rho_w*max(1.e-6_wp,theta_sat(k)-theta_i(k))))
          exit
        endif
      enddo
    else
      ! no frost table above the base of the column, or the top layer itself is frozen.
      ! z_wtab_max means "no perched zone": it is finite so the diagnostic stays usable, and
      ! since the aquifer table never exceeds it the min() below is then a no-op
      w_table_perch = hydro_par%z_wtab_max
    endif

    ! The shallower of the two controls surface saturation: a perched table above frozen ground
    ! wets the surface whatever the aquifer is doing. The two regimes carry DIFFERENT TOPMODEL
    ! decay factors, so the shift of the CTI threshold has to be taken from whichever one wins,
    ! it cannot be a single f times a single depth. f_drain is calibrated for the deep mineral
    ! aquifer, which varies over metres; the perched zone lives in a sub-metre active layer in
    ! organic soil, where the transmissivity decay is far steeper (f of order 1-10 1/m is
    ! standard in peatland hydrology). Using f_drain for both made the perched zone saturate
    ! about 40% of every permafrost cell, three times the observed wetland fraction.
    if (w_table_perch .lt. w_table) then
      w_table_eff = w_table_perch
      fz_eff      = hydro_par%f_drain_perch * w_table_perch
    else
      w_table_eff = w_table
      fz_eff      = hydro_par%f_drain * w_table
    endif

   return

  end subroutine perched_water_table

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  w _ t a b l e _ i n i
  !   Purpose    :  initial water table depth for a cold start
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function w_table_ini()

    implicit none

    real(wp) :: w_table_ini

    ! equilibrium depth for a recharge of 5% of the maximum baseflow, a typical humid-mid-latitude
    ! value, so that the initial state scales with the parameters instead of being a fixed depth
    w_table_ini = min(hydro_par%z_wtab_max, log(20._wp)/hydro_par%f_drain)

  end function w_table_ini


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  g r o u n d w a t e r
  !   Purpose    :  update the water table depth and compute baseflow
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine groundwater(theta,theta_w,drainage,w_table,runoff_gw, &
                        drainage_iso,w_gw_iso,runoff_gw_iso)

    implicit none

    real(wp), dimension(:), intent(in) :: theta, theta_w
    real(wp), intent(in) :: drainage      ! kg/m2/s, recharge from the bottom of the soil column
    real(wp), intent(inout) :: w_table    ! m, water table depth, positive downward
    real(wp), intent(out) :: runoff_gw    ! kg/m2/s, baseflow
    ! water-isotope siblings (always passed; values are 0 unless l_wiso=.true.)
    real(wp), dimension(:), intent(in)    :: drainage_iso
    real(wp), dimension(:), intent(inout) :: w_gw_iso
    real(wp), dimension(:), intent(out)   :: runoff_gw_iso

    integer :: iso
    real(wp) :: z, kappa, rsy, f_liq, w_gw_new
    real(wp) :: w_gw_iso_new(nwiso)


    ! no baseflow from a frozen column, same shutdown as in the SIMTOP subsurface runoff
    f_liq = sum(theta_w*dz(1:nl)) / max(1.e-10_wp, sum(theta*dz(1:nl)))
    kappa = hydro_par%kappa_max/sec_day * f_liq   ! kg/m2/s

    ! water table rise per unit of water added to the aquifer, m/(kg/m2)
    rsy = 1._wp/(rho_w*hydro_par%sy_aqf)

    ! 1) recharge from the soil column, raises the water table
    z = max(0._wp, w_table - drainage*dt*rsy)

    ! 2) baseflow, exact solution of rho_w*sy*dz/dt = kappa*exp(-f_drain*z) over one time step.
    !    Written with the exponential factored out so that exp(-f_drain*z) underflows harmlessly
    !    for a deep water table instead of exp(+f_drain*z) overflowing.
    z = z + 1._wp/hydro_par%f_drain &
          * log(1._wp + hydro_par%f_drain*kappa*dt*rsy*exp(-hydro_par%f_drain*z))
    z = min(z, hydro_par%z_wtab_max)

    ! 3) diagnose the baseflow from the change in aquifer storage, so that water is conserved
    !    exactly whatever the clipping in 1) and 2) did. Note that z can only have increased
    !    relative to step 1), so runoff_gw is always positive.
    runoff_gw = drainage + rho_w*hydro_par%sy_aqf*(z - w_table)*rdt   ! kg/m2/s

    if (l_wiso) then
      ! well mixed aquifer, recharge is mixed in before the baseflow is withdrawn
      w_gw_new     = rho_w*hydro_par%sy_aqf*(hydro_par%z_wtab_max - w_table) + drainage*dt  ! kg/m2
      w_gw_iso_new = w_gw_iso + drainage_iso*dt
      if (w_gw_new .gt. 0._wp) then
        do iso=1,nwiso
          runoff_gw_iso(iso) = runoff_gw * w_gw_iso_new(iso)/w_gw_new
        enddo
      else
        runoff_gw_iso = 0._wp
      endif
      w_gw_iso = w_gw_iso_new - runoff_gw_iso*dt
    else
      runoff_gw_iso = 0._wp
    endif

    w_table = z


   return

  end subroutine groundwater

end module lndvc_groundwater_mod
