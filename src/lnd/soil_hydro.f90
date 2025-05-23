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
  use timer, only : sec_day, time_soy_lnd
  use constants, only : rho_w, rho_i
  use control, only : check_water
  use lnd_grid, only : dz, rdz_neg, rdz_pos, nl
  use lnd_grid, only : nsurf, nveg, npft, nsoil, flag_veg, is_veg, flag_pft
  use lnd_params, only : dt, rdt
  use lnd_params, only : pft_par, hydro_par, snow_par
  use tridiag, only : tridiag_solve

   implicit none

   private
   public :: soil_hydro

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s o i l _ h y d r o
  !   Purpose    :  update soil liquid water content
  !              :  by solving the tridiagonal system
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine soil_hydro(frac_surf,mask_snow,theta_sat,k_sat,k_exp,psi_exp,kappa_int,psi,w_table, &
                       transpiration,evap_surface,infiltration,wilt, &
                       w_snow,w_w,w_i,theta_w,theta_i,theta,theta_w_cum,theta_i_cum,theta_fire_cum, &
                       drainage)

    implicit none

    integer, intent(inout) :: mask_snow
    real(wp), dimension(:), intent(in) :: frac_surf
    real(wp), intent(in) :: infiltration, w_table
    real(wp), dimension(:), intent(in) :: theta_sat, k_sat, kappa_int, psi 
    integer, dimension(:), intent(in) :: psi_exp, k_exp
    real(wp), dimension(:), intent(in) :: transpiration, evap_surface
    real(wp), dimension(:,:), intent(in) :: wilt
    real(wp), intent(inout) :: w_snow
    real(wp), dimension(:), intent(inout) :: w_w, w_i, theta_w, theta_i, theta, theta_w_cum, theta_i_cum
    real(wp), intent(inout) :: theta_fire_cum
    real(wp), dimension(:), intent(out) :: drainage

    integer :: i,j
    integer :: n, k, kk, i_here, i_print, j_print, n_print
    real(wp) :: frac_inf_2, wilt_tot, dw, ddw, run_sub, f_veg, water_bal
    real(wp), dimension(nl) :: inf
    real(wp), dimension(nl) :: r_sub
    real(wp), dimension(nl) :: q, e, dk_dtheta, dpsi_dtheta, dq_1_dtheta_1, dq_1_dtheta, dq_dtheta, dq_dtheta1
    real(wp), dimension(nl) :: a, b, c, r, x, ww_old
    real(wp), dimension(nl) :: w_w_max
    real(wp), parameter :: w_w_min = 0.01_wp ! mm, kg/m2


    i_print = 23
    j_print = 10
    n_print = 3
 

    ww_old = w_w

    inf(:) = 0._wp
    ! all infiltration into the top layer or spread over top two layers
    frac_inf_2 = hydro_par%frac_inf_2 * (1._wp - theta_i(2)/(theta_w(2)+theta_i(2)+0.001_wp))
    inf(1) = (1._wp-frac_inf_2)*infiltration
    inf(2) = frac_inf_2*infiltration

    e(:) = 0._wp
    
    drainage = 0._wp

    f_veg = sum(frac_surf,mask=flag_veg.eq.1)

    ! evaporation/transpiration from layers
    do n=1,nveg
     if( frac_surf(n) .gt. 0._wp ) then
      ! evaporation from or dew deposition to first soil layer
      if(mask_snow .eq. 0) then ! no snow
       e(1) = e(1) + evap_surface(n) * frac_surf(n)/f_veg
      endif
      ! transpiration from vegetation
      if( flag_pft(n) .eq. 1 ) then
       wilt_tot = sum(wilt(:,n) * pft_par%root_frac(:,n))
       do k=1,nl
        e(k) = e(k) &
             + (transpiration(n) * pft_par%root_frac(k,n) * wilt(k,n) / max(1.d-20,wilt_tot)) & ! kg/m2/s
             * frac_surf(n)/f_veg
       enddo
      endif
     endif
    enddo

    if( hydro_par%i_runoff .eq. 1 ) then
     r_sub(:) = 0._wp
    elseif( hydro_par%i_runoff .eq. 2 ) then
     ! subsurface runoff
     run_sub = hydro_par%kappa_max/sec_day * exp(-2._wp*hydro_par%f_wtab*w_table) * sum(theta_w*dz(1:nl))/sum(theta*dz(1:nl)) ! kg/m2/s
     do k=1,nl
      r_sub(k) = run_sub * kappa_int(k)*dz(k) / sum(kappa_int*dz(1:nl))
     enddo
    endif


    do k=1,nl-1
     q(k) = -kappa_int(k) * ( (psi(k+1) - psi(k)) * rdz_pos(k) - 1._wp )
     if( hydro_par%i_cond_theta .eq. 1 ) then
      ! use liquid AND frozen water
      dk_dtheta(k) = k_exp(k) * k_sat(k) &
                   * ( max(hydro_par%theta_min, (theta(k)+theta(k+1))/(theta_sat(k)+theta_sat(k+1))) )**(k_exp(k)-1) &
                   / (theta_sat(k) + theta_sat(k+1))
     elseif( hydro_par%i_cond_theta .eq. 2 ) then
      ! use only liquid water
      dk_dtheta(k) = k_exp(k) * k_sat(k) &
                   * ( max(hydro_par%theta_min, (theta_w(k)+theta_w(k+1))/(theta_sat(k)+theta_sat(k+1))) ) **(k_exp(k)-1) &
                   / (theta_sat(k) + theta_sat(k+1))
     elseif( hydro_par%i_cond_theta .eq. 3 ) then
      ! use only liquid water
      dk_dtheta(k) = k_exp(k) * k_sat(k) &
                   * ( max(hydro_par%theta_min, 0.5_wp*(theta(k)+theta(k+1)+theta_w(k)+theta_w(k+1))/(theta_sat(k)+theta_sat(k+1))) ) **(k_exp(k)-1) &
                   / (theta_sat(k) + theta_sat(k+1))
     endif
    enddo
    ! bottom boundary condition
    k = nl
    if( hydro_par%i_runoff .eq. 1 ) then
     q(k) = kappa_int(k) ! free drainage
     if( hydro_par%i_cond_theta .eq. 1 ) then
      ! use liquid AND frozen water
      dk_dtheta(k) = k_exp(k) * k_sat(k) &
                   * (max(hydro_par%theta_min, theta(k)/theta_sat(k)))**(k_exp(k)-1) / theta_sat(k)
     elseif( hydro_par%i_cond_theta .eq. 2 ) then
      ! use only liquid water
      dk_dtheta(k) = k_exp(k) * k_sat(k) &
                   * (max(hydro_par%theta_min, theta_w(k)/theta_sat(k)))**(k_exp(k)-1) / theta_sat(k)
     elseif( hydro_par%i_cond_theta .eq. 3 ) then
      ! use only liquid water
      dk_dtheta(k) = k_exp(k) * k_sat(k) &
                   * (max(hydro_par%theta_min, 0.5_wp*(theta(k)+theta_w(k))/theta_sat(k)))**(k_exp(k)-1) / theta_sat(k)
     endif
    elseif( hydro_par%i_runoff .eq. 2 ) then
     q(k) = 0._wp                ! no drainage
     dk_dtheta(k) = 0._wp
    endif

    do k=1,nl
    if( hydro_par%i_cond_theta .eq. 1 ) then
      ! use liquid AND frozen water
      dpsi_dtheta(k) = psi_exp(k) * psi(k) / max(hydro_par%theta_min*theta_sat(k), theta(k))
     elseif( hydro_par%i_cond_theta .eq. 2 ) then
      ! use only liquid water
      dpsi_dtheta(k) = psi_exp(k) * psi(k) / max(hydro_par%theta_min*theta_sat(k), theta_w(k))
     elseif( hydro_par%i_cond_theta .eq. 3 ) then
      ! use only liquid water
      dpsi_dtheta(k) = psi_exp(k) * psi(k) / max(hydro_par%theta_min*theta_sat(k), 0.5_wp*(theta(k)+theta_w(k)))
     endif
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
    r(k) = - inf(k) + q(k) + e(k) + r_sub(k) 

    ! intermediate layers
    do k=2,nl-1
     a(k) = dq_1_dtheta_1(k)
     b(k) = - dq_dtheta(k) + dq_1_dtheta(k) - rho_w*dz(k)*rdt
     c(k) = - dq_dtheta1(k)
     r(k) = - inf(k) - q(k-1) + q(k) + e(k) + r_sub(k)
    enddo

    ! bottom layer, k=nl, free drainage boundary condition: q_N = -k_N
    k = nl
    a(k) = dq_1_dtheta_1(k)
    b(k) = - dk_dtheta(k) + dq_1_dtheta(k) - rho_w*dz(k)*rdt
    c(k) = 0._wp
    r(k) = - q(k-1) + q(k) + e(k) + r_sub(k)

    ! solve tridiagonal system for liquid volumetric water content change, x [m3/m3]
    call tridiag_solve(a,b,c,r,x,nl)

    ! update soil water content, [kg/m2 or mm]
    w_w = w_w + x*dz(1:nl)*rho_w


    if( hydro_par%i_runoff .eq. 1 ) then
     ! drainage is equal to -q(nl) at the new time step, expanded in Taylor series:
     drainage(is_veg) = kappa_int(nl) + dk_dtheta(nl)*x(nl)
    elseif( hydro_par%i_runoff .eq. 2 ) then
     drainage(is_veg) = run_sub        
    endif

    ! check wether w_w_min < w_w < (theta_sat-theta_i)*dz, w_w_min = 0.01 [kg/m2]
    ! first check for excess liquid water 
    w_w_max = (theta_sat-theta_i)*dz(1:nl)*rho_w
    do k=1,nl-1 ! start from top layer
     if(w_w(k) .gt. w_w_max(k)) then
      if(k.eq.1) then
       ! half of excess water to drainage and half to layer below
       w_w(k+1) = w_w(k+1) + 0.5_wp*(w_w(k) - w_w_max(k)) ! move excess water one layer down
       drainage(is_veg) = drainage(is_veg) + 0.5_wp*(w_w(k) - w_w_max(k)) * rdt 
      else     
       w_w(k+1) = w_w(k+1) + w_w(k) - w_w_max(k) ! move excess water one layer down
      endif 
      w_w(k) = w_w_max(k) ! reset water to saturation value
     endif
    enddo
    if(w_w(nl) .gt. w_w_max(nl)) then ! bottom layer oversaturated
     drainage(is_veg) = drainage(is_veg) &
              + (w_w(nl) - w_w_max(nl)) * rdt ! add excess to drainage, kg/m2/s
     w_w(nl) = w_w_max(nl)
    endif

 
    ! then check for negative liquid water
    if( minval(w_w) .lt. w_w_min) then
     do k=1,nl-1
      if(w_w(k) .lt. w_w_min) then
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
         ddw = ddw + (w_w(kk) - w_w_min)
         w_w(kk) = w_w_min ! reset to w_w_min
        enddo
        w_w(k) = w_w(k) - (w_w_min-w_w(nl)-ddw) ! extract only the required amount from layer k
        w_w(nl) = w_w_min ! set bottom layer to w_w_min
        exit lp1
       endif
      enddo lp1
      if(i_here .eq. 0) then ! not enough water in the layers, remove additional necessary water from drainage
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
     water_bal =  sum(inf)*dt - sum(e)*dt - drainage(is_veg)*dt &
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
