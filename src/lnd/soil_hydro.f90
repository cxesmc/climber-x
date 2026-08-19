!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s o i l _ h y d r o _ m o d
!
!  Purpose : soil hydrology
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
module soil_hydro_mod

  use precision, only : wp
  use timer, only : time_soy_lnd
  use constants, only : rho_w, rho_i
  use control, only : check_water
  use lnd_grid, only : dz, rdz_neg, rdz_pos, nl
  use lnd_grid, only : nsurf, nveg, npft, nsoil, flag_veg, is_veg, flag_pft
  use lnd_params, only : dt, rdt
  use lnd_params, only : pft_par, hydro_par, snow_par
  use tridiag, only : tridiag_solve
  use wiso_params, only : l_wiso, nwiso, i_o18, Rstd

   implicit none

   private
   public :: soil_hydro

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s o i l _ h y d r o
  !   Purpose    :  update soil liquid water content
  !              :  by solving the tridiagonal system
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine soil_hydro(frac_surf,mask_snow,theta_sat,k_sat,k_exp,psi_exp,kappa_int,psi, &
                       transpiration,evap_surface,infiltration,wilt, &
                       w_snow,w_w,w_i,theta_w,theta_i,theta,theta_w_cum,theta_i_cum,theta_fire_cum, &
                       drainage,runoff_exc, &
                       transpiration_iso,evap_surface_iso,infiltration_iso, &
                       w_w_iso,drainage_iso,runoff_exc_iso)

    implicit none

    integer, intent(inout) :: mask_snow
    real(wp), dimension(:), intent(in) :: frac_surf
    real(wp), intent(in) :: infiltration
    real(wp), dimension(:), intent(in) :: theta_sat, k_sat, kappa_int, psi
    integer, dimension(:), intent(in) :: psi_exp, k_exp
    real(wp), dimension(:), intent(in) :: transpiration, evap_surface
    real(wp), dimension(:,:), intent(in) :: wilt
    real(wp), intent(inout) :: w_snow
    real(wp), dimension(:), intent(inout) :: w_w, w_i, theta_w, theta_i, theta, theta_w_cum, theta_i_cum
    real(wp), intent(inout) :: theta_fire_cum
    real(wp), dimension(:), intent(out) :: drainage
    real(wp), intent(out) :: runoff_exc      ! saturation excess leaving at the surface, kg/m2/s
    ! water-isotope siblings (always passed; values are 0 unless l_wiso=.true.)
    real(wp), dimension(:,:), intent(in)    :: transpiration_iso, evap_surface_iso
    real(wp), dimension(:),   intent(in)    :: infiltration_iso
    real(wp), dimension(:,:), intent(inout) :: w_w_iso
    real(wp), dimension(:,:), intent(out)   :: drainage_iso
    real(wp), dimension(:),   intent(out)   :: runoff_exc_iso

    integer :: i,j
    integer :: n, k, kk, i_here, i_print, j_print, n_print, iso
    real(wp) :: wilt_tot, dw, ddw, f_veg, water_bal
    real(wp), dimension(nl) :: inf
    real(wp), dimension(nl) :: q, e, dk_dtheta, dpsi_dtheta, dq_1_dtheta_1, dq_1_dtheta, dq_dtheta, dq_dtheta1
    real(wp), dimension(nl) :: a, b, c, r, x, ww_old
    real(wp), dimension(nl) :: w_w_max
    real(wp), dimension(nl) :: q_new
    real(wp), dimension(nl,nwiso) :: inf_iso_l, e_iso, q_iso, ww_old_iso
    real(wp) :: ratio(nwiso), excess_iso(nwiso), deficit_iso(nwiso)
    real(wp) :: excess, acc, dw_iso
    real(wp), parameter :: w_w_min = 0.01_wp ! mm, kg/m2


    i_print = 23
    j_print = 10
    n_print = 3
 

    ww_old = w_w
    runoff_exc = 0._wp
    if (l_wiso) runoff_exc_iso = 0._wp
    if (l_wiso) ww_old_iso = w_w_iso

    inf(:) = 0._wp
    if (l_wiso) inf_iso_l(:,:) = 0._wp
    ! all infiltration into the top layer
    inf(1) = infiltration
    if (l_wiso) then
      do iso=1,nwiso
        inf_iso_l(1,iso) = infiltration_iso(iso)
      enddo
    endif

    e(:) = 0._wp
    if (l_wiso) e_iso(:,:) = 0._wp

    drainage = 0._wp
    if (l_wiso) drainage_iso = 0._wp

    f_veg = sum(frac_surf,mask=flag_veg.eq.1)

    ! evaporation/transpiration from layers
    do n=1,nveg
     if( frac_surf(n) .gt. 0._wp ) then
      ! evaporation from or dew deposition to first soil layer
      if(mask_snow .eq. 0) then ! no snow
       e(1) = e(1) + evap_surface(n) * frac_surf(n)/f_veg
       if (l_wiso) then
         do iso=1,nwiso
           e_iso(1,iso) = e_iso(1,iso) + evap_surface_iso(n,iso) * frac_surf(n)/f_veg
         enddo
       endif
      endif
      ! transpiration from vegetation
      if( flag_pft(n) .eq. 1 ) then
       wilt_tot = sum(wilt(:,n) * pft_par%root_frac(:,n))
       do k=1,nl
        e(k) = e(k) &
             + (transpiration(n) * pft_par%root_frac(k,n) * wilt(k,n) / max(1.d-20,wilt_tot)) & ! kg/m2/s
             * frac_surf(n)/f_veg
        if (l_wiso) then
          do iso=1,nwiso
            e_iso(k,iso) = e_iso(k,iso) &
                         + (transpiration_iso(n,iso) * pft_par%root_frac(k,n) * wilt(k,n) / max(1.d-20,wilt_tot)) &
                         * frac_surf(n)/f_veg
          enddo
        endif
       enddo
      endif
     endif
    enddo


    do k=1,nl-1
     q(k) = -kappa_int(k) * ( (psi(k+1) - psi(k)) * rdz_pos(k) - 1._wp )
     ! liquid water only, consistent with kappa_int in soil_par_hydro()
     dk_dtheta(k) = k_exp(k) * k_sat(k) &
                  * ( max(hydro_par%theta_min, (theta_w(k)+theta_w(k+1))/(theta_sat(k)+theta_sat(k+1))) ) **(k_exp(k)-1) &
                  / (theta_sat(k) + theta_sat(k+1))
    enddo
    ! bottom boundary condition
    k = nl
    q(k) = kappa_int(k) ! free drainage
    dk_dtheta(k) = k_exp(k) * k_sat(k) &
                 * (max(hydro_par%theta_min, theta_w(k)/theta_sat(k)))**(k_exp(k)-1) / theta_sat(k)

    do k=1,nl
     dpsi_dtheta(k) = psi_exp(k) * psi(k) / max(hydro_par%theta_min*theta_sat(k), theta_w(k))
    enddo

    do k=2,nl
     dq_1_dtheta_1(k) = kappa_int(k-1) * rdz_neg(k) * dpsi_dtheta(k-1) & 
                      - dk_dtheta(k-1)*((psi(k) - psi(k-1)) * rdz_neg(k) - 1._wp)
     dq_1_dtheta(k)   = - (kappa_int(k-1) * rdz_neg(k) * dpsi_dtheta(k)) & 
                      - dk_dtheta(k-1)*((psi(k) - psi(k-1)) * rdz_neg(k) - 1._wp)
    enddo
    do k=1,nl-1
     dq_dtheta(k)     = kappa_int(k)  * rdz_pos(k) * dpsi_dtheta(k) &    
                      - dk_dtheta(k)  *((psi(k+1)   - psi(k)) * rdz_pos(k) - 1._wp)
     dq_dtheta1(k)    = - (kappa_int(k) * rdz_pos(k) * dpsi_dtheta(k+1)) &    
                      - dk_dtheta(k)  *((psi(k+1)   - psi(k)) * rdz_pos(k) - 1._wp)
    enddo

    ! top layer, k=1
    k = 1
    a(k) = 0._wp
    b(k) = - dq_dtheta(k) - rho_w*dz(k)*rdt
    c(k) = - dq_dtheta1(k)
    r(k) = - inf(k) + q(k) + e(k)

    ! intermediate layers
    do k=2,nl-1
     a(k) = dq_1_dtheta_1(k)
     b(k) = - dq_dtheta(k) + dq_1_dtheta(k) - rho_w*dz(k)*rdt
     c(k) = - dq_dtheta1(k)
     r(k) = - inf(k) - q(k-1) + q(k) + e(k)
    enddo

    ! bottom layer, k=nl, free drainage boundary condition: q_N = -k_N
    k = nl
    a(k) = dq_1_dtheta_1(k)
    b(k) = - dk_dtheta(k) + dq_1_dtheta(k) - rho_w*dz(k)*rdt
    c(k) = 0._wp
    r(k) = - q(k-1) + q(k) + e(k)

    ! solve tridiagonal system for liquid volumetric water content change, x [m3/m3]
    call tridiag_solve(a,b,c,r,x,nl)

    ! update soil water content, [kg/m2 or mm]
    w_w = w_w + x*dz(1:nl)*rho_w


    ! drainage is equal to -q(nl) at the new time step, expanded in Taylor series:
    drainage(is_veg) = kappa_int(nl) + dk_dtheta(nl)*x(nl)

    ! ---------- water isotope advection ----------
    ! Reconstruct new-time inter-layer bulk fluxes from the solver (same equations as bulk),
    ! advect iso mass using upwind donor ratio (OLD-time), and update w_w_iso by layer mass balance.
    if (l_wiso) then
      ! new-time bulk fluxes at each interface (positive = downward)
      do k=1,nl-1
        q_new(k) = q(k) + dq_dtheta(k)*x(k) + dq_dtheta1(k)*x(k+1)
      enddo
      q_new(nl) = q(nl) + dk_dtheta(nl)*x(nl)   ! free drainage
      ! iso fluxes at each interface (upwind from donor's OLD ratio)
      do k=1,nl-1
        if (q_new(k) .ge. 0._wp) then
          if (ww_old(k).gt.0._wp) then
            do iso=1,nwiso
              q_iso(k,iso) = q_new(k) * ww_old_iso(k,iso)/ww_old(k)
            enddo
          else
            q_iso(k,:) = 0._wp
          endif
        else
          if (ww_old(k+1).gt.0._wp) then
            do iso=1,nwiso
              q_iso(k,iso) = q_new(k) * ww_old_iso(k+1,iso)/ww_old(k+1)
            enddo
          else
            q_iso(k,:) = 0._wp
          endif
        endif
      enddo
      ! bottom flux iso
      if (q_new(nl).gt.0._wp .and. ww_old(nl).gt.0._wp) then
        do iso=1,nwiso
          q_iso(nl,iso) = q_new(nl) * ww_old_iso(nl,iso)/ww_old(nl)
        enddo
      else
        q_iso(nl,:) = 0._wp
      endif
      ! per-layer iso mass update (q_iso(0) := 0; influx from above = q_iso(k-1), outflux below = q_iso(k))
      do iso=1,nwiso
        w_w_iso(1,iso) = ww_old_iso(1,iso) &
                        + dt * (inf_iso_l(1,iso) - q_iso(1,iso) - e_iso(1,iso))
        do k=2,nl
          w_w_iso(k,iso) = ww_old_iso(k,iso) &
                          + dt * (inf_iso_l(k,iso) + q_iso(k-1,iso) - q_iso(k,iso) &
                                 - e_iso(k,iso))
        enddo
      enddo
      ! drainage iso, consistent with the bulk: the flux at the bottom of the column
      drainage_iso(is_veg,:) = q_iso(nl,:)
    endif

    ! check wether w_w_min < w_w < (theta_sat-theta_i)*dz, w_w_min = 0.01 [kg/m2]
    ! first check for excess liquid water
    !
    ! Redistribute liquid water in excess of the ice-reduced pore space. Every transfer is
    ! limited by the FREE CAPACITY of the receiving layer, so a frozen layer, whose w_w_max is
    ! small because the pores are filled with ice, actually blocks the flow. Water that cannot
    ! move down is passed back up and leaves at the surface as saturation excess.
    w_w_max = (theta_sat-theta_i)*dz(1:nl)*rho_w

    ! downward, never more than the layer below can accept
    do k=1,nl-1
     excess = w_w(k) - w_w_max(k)
     if (excess .gt. 0._wp) then
      acc = min(excess, max(0._wp, w_w_max(k+1)-w_w(k+1)))
      if (acc .gt. 0._wp) then
       if (l_wiso .and. w_w(k).gt.0._wp) then
         do iso=1,nwiso
           dw_iso = w_w_iso(k,iso)/w_w(k) * acc
           w_w_iso(k+1,iso) = w_w_iso(k+1,iso) + dw_iso
           w_w_iso(k,iso)   = w_w_iso(k,iso)   - dw_iso
         enddo
       endif
       w_w(k+1) = w_w(k+1) + acc
       w_w(k)   = w_w(k)   - acc
      endif
     endif
    enddo

    ! genuine saturation excess out of the base of the column, this feeds the aquifer
    if(w_w(nl) .gt. w_w_max(nl)) then
     excess = w_w(nl) - w_w_max(nl)
     if (l_wiso .and. w_w(nl).gt.0._wp) then
       do iso=1,nwiso
         dw_iso = w_w_iso(nl,iso)/w_w(nl) * excess
         drainage_iso(is_veg,iso) = drainage_iso(is_veg,iso) + dw_iso * rdt
         w_w_iso(nl,iso) = w_w_iso(nl,iso) - dw_iso
       enddo
     endif
     drainage(is_veg) = drainage(is_veg) + excess * rdt ! add excess to drainage, kg/m2/s
     w_w(nl) = w_w_max(nl)
    endif

    ! upward, for water stalled below a frozen barrier. Unlike the downward pass this one is NOT
    ! capacity limited: the destination is the surface, which is an unlimited sink, so there is
    ! no barrier to respect and the water simply leaves the column. 
    do k=nl,2,-1
     excess = w_w(k) - w_w_max(k)
     if (excess .gt. 0._wp) then
      if (l_wiso .and. w_w(k).gt.0._wp) then
        do iso=1,nwiso
          dw_iso = w_w_iso(k,iso)/w_w(k) * excess
          w_w_iso(k-1,iso) = w_w_iso(k-1,iso) + dw_iso
          w_w_iso(k,iso)   = w_w_iso(k,iso)   - dw_iso
        enddo
      endif
      w_w(k-1) = w_w(k-1) + excess
      w_w(k)   = w_w_max(k)
     endif
    enddo

    ! whatever the column cannot hold at all leaves at the surface
    if(w_w(1) .gt. w_w_max(1)) then
     excess = w_w(1) - w_w_max(1)
     if (l_wiso .and. w_w(1).gt.0._wp) then
       do iso=1,nwiso
         dw_iso = w_w_iso(1,iso)/w_w(1) * excess
         runoff_exc_iso(iso) = runoff_exc_iso(iso) + dw_iso * rdt
         w_w_iso(1,iso) = w_w_iso(1,iso) - dw_iso
       enddo
     endif
     runoff_exc = runoff_exc + excess * rdt
     w_w(1) = w_w_max(1)
    endif


    ! then check for negative liquid water
    if( minval(w_w) .lt. w_w_min) then
     do k=1,nl-1
      if(w_w(k) .lt. w_w_min) then
       ! pull water from layer below; iso at donor (k+1) ratio, Rstd fallback when donor empty
       if (l_wiso) then
         do iso=1,nwiso
           if (w_w(k+1) .gt. 0._wp) then
             deficit_iso(iso) = (w_w_iso(k+1,iso)/w_w(k+1)) * (w_w_min - w_w(k))
           else
             deficit_iso(iso) = Rstd(iso) * (w_w_min - w_w(k))   ! donor empty: use standard ratio for mass balance
           endif
           w_w_iso(k+1,iso) = w_w_iso(k+1,iso) - deficit_iso(iso)
           w_w_iso(k,iso)   = w_w_min * Rstd(iso)
         enddo
       endif
       w_w(k+1) = w_w(k+1) - (w_w_min - w_w(k)) ! take necessary water from layer below
       w_w(k) = w_w_min    ! set to w_w_min
      endif
     enddo
     if(w_w(nl) .lt. w_w_min) then ! bottom layer negative liquid water
      dw = 0._wp
      i_here = 0
lp1:  do k=nl-1,1,-1 ! search for the required amount of water in layers above
       dw = dw + w_w(k) - w_w_min
       if(dw .ge. w_w_min-w_w(nl)) then ! found water enough to fill bottom layer to w_w_min
        i_here = 1
        ddw = 0._wp
        do kk=nl-1,k+1,-1 ! extract all water from the layers nl-1,k+1
         if (l_wiso .and. w_w(kk).gt.0._wp) then
           do iso=1,nwiso
             w_w_iso(kk,iso) = w_w_iso(kk,iso) - (w_w_iso(kk,iso)/w_w(kk)) * (w_w(kk) - w_w_min)
           enddo
         endif
         ddw = ddw + (w_w(kk) - w_w_min)
         w_w(kk) = w_w_min ! reset to w_w_min
        enddo
        if (l_wiso .and. w_w(k).gt.0._wp) then
          do iso=1,nwiso
            w_w_iso(k,iso) = w_w_iso(k,iso) - (w_w_iso(k,iso)/w_w(k)) * (w_w_min-w_w(nl)-ddw)
            w_w_iso(nl,iso) = w_w_min * Rstd(iso)
          enddo
        endif
        w_w(k) = w_w(k) - (w_w_min-w_w(nl)-ddw) ! extract only the required amount from layer k
        w_w(nl) = w_w_min ! set bottom layer to w_w_min
        exit lp1
       endif
      enddo lp1
      if(i_here .eq. 0) then ! not enough water in the layers, remove additional necessary water from drainage
       if (l_wiso) then
         do iso=1,nwiso
           ! drainage absorbs the (negative) deficit at VSMOW (no donor available)
           drainage_iso(is_veg,iso) = drainage_iso(is_veg,iso) - (w_w_min - w_w(nl) - dw) * Rstd(iso) * rdt
           w_w_iso(:,iso) = w_w_min * Rstd(iso)
         enddo
       endif
       drainage(is_veg)  = drainage(is_veg) - (w_w_min - w_w(nl) - dw) * rdt
       w_w = w_w_min ! set w_w to w_w_min in all layers
      endif

      if(check_water .and. drainage(is_veg)*dt .lt. -0.1_wp) then
       print *,' '
       print *,' '
       print *,'WARNING, negative drainage!! ',drainage(is_veg)*dt,' mm/day'
       print *,'infiltration',infiltration*dt
       print *,' '
       print *,' '
       !stop
      endif

     endif
    endif

    ! update volumetric water content, m3/m3
    do k=1,nl
     theta(k) = w_w(k)/(dz(k)*rho_w) + w_i(k)/(dz(k)*rho_i) ! total (liquid + frozen) volumetric water content
     theta_w(k) = w_w(k)/(dz(k)*rho_w) ! liquid water content
     theta_i(k) = w_i(k)/(dz(k)*rho_i) ! frozen water content
    enddo


    if( check_water ) then
     ! water balance
     water_bal =  sum(inf)*dt - sum(e)*dt - drainage(is_veg)*dt - runoff_exc*dt &
               - sum(w_w-ww_old)

     if(abs(water_bal).gt.1.d-7) then
      print *,'mask_snow,i,j',mask_snow,i,j
      print *,'water balance',water_bal
      print *,'sum(frac_surf)',sum(frac_surf)
      print *,'frac_sur',frac_surf
      print *,'dw_soil',sum(w_w-ww_old)
      print *,'sum(inf)',sum(inf)*dt,infiltration*dt
      print *,'sum(e)',sum(e)*dt
      print *,'drain',drainage(is_veg)*dt
      print *,'runoff_exc',runoff_exc*dt
      stop
     endif
    endif

    ! update snow mask over icefree grid cell
    if( w_snow .gt. snow_par%w_snow_crit ) then
     mask_snow = 1
    else
     mask_snow = 0
    endif

    if (time_soy_lnd) then
      theta_w_cum = 0._wp
      theta_i_cum = 0._wp
      theta_fire_cum = 0._wp
    endif
    ! cumulate soil moisture for soil carbon decomposition
    theta_w_cum = theta_w_cum + theta_w
    theta_i_cum = theta_i_cum + theta_i
    ! cumulate top soil moisture for fire disturbance rate
    theta_fire_cum = theta_fire_cum + theta_w(1) + theta_i(1)

    return

  end subroutine soil_hydro

end module soil_hydro_mod
