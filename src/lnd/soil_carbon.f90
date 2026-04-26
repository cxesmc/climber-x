!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s o i l _ c a r b o n _ m o d
!
!  Purpose : soil carbon dynamics
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
module soil_carbon_mod

  use precision, only : wp
  use timer, only : day_mon
  use constants, only : c14_tdec
  use control, only : check_carbon
  use lnd_grid, only : z, dz_c, rdz_c, rdz_pos_c, rdz_neg_c
  use lnd_grid, only : nl, nlc, ncarb, ic_min
  use lnd_params, only : dt_c, dt_day_c
  use lnd_params, only : soilc_par, ch4_par
  use tridiag, only : tridiag_solve

  implicit none

  private
  public :: soil_carbon

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s o i l _ c a r b o n
  !   Purpose    :  solve soil carbon prognostic equations
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine soil_carbon(f_veg,f_wetland,f_peat,f_crop,f_pasture,litterfall,litterfall13,litterfall14, &
                        litter_c,fast_c,slow_c,litter_c13,fast_c13,slow_c13,litter_c14,fast_c14,slow_c14, &
                        k_litter,k_fast,k_slow,k_litter_wet,k_fast_wet,k_slow_wet,diff_soilc,adv_soilc,ch4_frac_wet, &
                        soil_resp,soil_resp13,soil_resp14,soil_resp_l,soil_c_tot,soil_c13_tot,soil_c14_tot, &
                        ch4_emis_wetland,c13h4_emis_wetland,carbon_bal_soil,carbon13_bal_soil,carbon14_bal_soil)

    implicit none

    real(wp), intent(in) :: f_veg, f_wetland, f_peat, f_crop, f_pasture
    real(wp), dimension(:), intent(in) :: litterfall, litterfall13, litterfall14
    real(wp), dimension(:), intent(inout) :: litter_c, fast_c, slow_c, litter_c13, fast_c13, slow_c13, litter_c14, fast_c14, slow_c14
    real(wp), dimension(:), intent(inout) :: k_litter, k_fast, k_slow, k_litter_wet, k_fast_wet, k_slow_wet, diff_soilc, adv_soilc
    real(wp), dimension(:), intent(inout) :: ch4_frac_wet
    real(wp), intent(inout) :: soil_resp, soil_resp13, soil_resp14
    real(wp), dimension(:), intent(inout) :: soil_resp_l
    real(wp), intent(out) :: soil_c_tot, soil_c13_tot, soil_c14_tot
    real(wp), intent(out) :: ch4_emis_wetland, c13h4_emis_wetland
    real(wp), intent(out) :: carbon_bal_soil, carbon13_bal_soil, carbon14_bal_soil

    integer :: k
    real(wp), dimension(1:nlc) :: a, b, c, r, x
    real(wp), dimension(1:nlc) :: litter_c_old, fast_c_old, slow_c_old
    real(wp), dimension(1:nl) :: litter_resp, fast_resp, slow_resp
    real(wp), dimension(1:nlc) :: litter_c13_old, fast_c13_old, slow_c13_old
    real(wp), dimension(1:nl) :: litter_resp13, fast_resp13, slow_resp13
    real(wp), dimension(1:nlc) :: litter_c14_old, fast_c14_old, slow_c14_old
    real(wp), dimension(1:nl) :: litter_resp14, fast_resp14, slow_resp14
    real(wp), dimension(1:nlc) :: klitter, kfast, kslow, diff, adv
    real(wp) :: f_inund


    ! save old prognostic variables
    litter_c_old = litter_c
    fast_c_old = fast_c
    slow_c_old = slow_c
    litter_c13_old = litter_c13
    fast_c13_old = fast_c13
    slow_c13_old = slow_c13
    litter_c14_old = litter_c14
    fast_c14_old = fast_c14
    slow_c14_old = slow_c14

    ! decomposition rate, 1/s
    ! inudated fraction of vegetated gridcell, excluding peatlands
    ! f_wetland is relative to the vegetated gridcell part, f_peat is relative to the whole gridcell
    ! f_inund is the inundated fraction of the vegetated gridcell part
    f_inund = max(0._wp,f_wetland-f_peat/f_veg)
    ! grid cell mean decomposition rates assuming NO respiration (k=0) of carbon buried under the ice sheets
    klitter = ((1._wp-f_inund)*k_litter + f_inund*k_litter_wet) &  ! for vegetated land fraction
            / real(dt_day_c,wp)*day_mon
    kfast   = ((1._wp-f_inund)*k_fast   + f_inund*k_fast_wet) &
            / real(dt_day_c,wp)*day_mon
    kslow   = ((1._wp-f_inund)*k_slow   + f_inund*k_slow_wet) &
            / real(dt_day_c,wp)*day_mon

    ! add additional contribution due to land use change
    do k=1,nl
      ! litter affected by both cropland and pasture expansion
      klitter(k) = klitter(k) + soilc_par%gamma_luc_Clitter_crop*f_crop + soilc_par%gamma_luc_Clitter_pasture*f_pasture 
      if (z(k).lt.0.5_wp) then
        ! fast soil carbon pool affected only by cropland and only in the top 0.5 meter of soil
        kfast(k) = kfast(k) + soilc_par%gamma_luc_Csoil_crop*f_crop
      endif
    enddo

    ! vertical carbon diffusivity
    diff = diff_soilc / real(dt_day_c,wp)*day_mon
    ! vertical carbon advection
    adv = adv_soilc / real(dt_day_c,wp)*day_mon

    ! reset cumulated values
    k_litter = 0._wp
    k_fast = 0._wp
    k_slow = 0._wp
    k_litter_wet = 0._wp
    k_fast_wet = 0._wp
    k_slow_wet = 0._wp
    diff_soilc = 0._wp
    adv_soilc  = 0._wp

    ! LITTER

    ! upstream scheme for advection
    ! top soil layer
    k = 1
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - klitter(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - adv(k)*rdz_c(k)
    c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k)
    r(k) = -litter_c(k)/dt_c - litterfall(k)*rdz_c(k)  ! kgC/m3/s

    ! other layers
    do k=2,nlc
     a(k) = diff(k-1)*rdz_neg_c(k)*rdz_c(k) + adv(k-1)*rdz_c(k)
     b(k) = -1._wp/dt_c - klitter(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - diff(k-1)*rdz_neg_c(k)*rdz_c(k) - adv(k)*rdz_c(k) 
     c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k)
     r(k) = -litter_c(k)/dt_c - litterfall(k)*rdz_c(k)
    enddo

    ! solve tridiagonal system
    call tridiag_solve(a,b,c,r,x,nlc)
    ! assign new litter carbon
    litter_c = x

    ! litter respiration to atmosphere
    litter_resp = soilc_par%f_resp_litter*klitter(1:nl)*litter_c(1:nl)*dz_c(1:nl)   ! kgC/m2/s

    ! FAST SOIL CARBON

    ! top soil layer
    k = 1
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - kfast(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - adv(k)*rdz_c(k)
    c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
    r(k) = -fast_c(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_fast*klitter(k)*litter_c(k)  ! kgC/m3/s

    ! intermediate layers
    do k=2,nlc
     a(k) = diff(k-1)*rdz_neg_c(k)*rdz_c(k) + adv(k-1)*rdz_c(k)
     b(k) = -1._wp/dt_c - kfast(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - diff(k-1)*rdz_neg_c(k)*rdz_c(k) - adv(k)*rdz_c(k) 
     c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
     r(k) = -fast_c(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_fast*klitter(k)*litter_c(k)
    enddo

    ! solve tridiagonal system
    call tridiag_solve(a,b,c,r,x,nlc)
    ! assign new litter carbon
    fast_c = x

    fast_resp = kfast(1:nl)*fast_c(1:nl)*dz_c(1:nl)  ! kgC/m2/s

    ! SLOW SOIL CARBON

    ! top soil layer
    k = 1
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - kslow(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - adv(k)*rdz_c(k)
    c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
    r(k) = -slow_c(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_slow*klitter(k)*litter_c(k)  ! kgC/m3/s

    ! intermediate layers
    do k=2,nlc
     a(k) = diff(k-1)*rdz_neg_c(k)*rdz_c(k) + adv(k-1)*rdz_c(k)
     b(k) = -1._wp/dt_c - kslow(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - diff(k-1)*rdz_neg_c(k)*rdz_c(k) - adv(k)*rdz_c(k) 
     c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
     r(k) = -slow_c(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_slow*klitter(k)*litter_c(k)
    enddo

    ! solve tridiagonal system
    call tridiag_solve(a,b,c,r,x,nlc)
    ! assign new litter carbon
    slow_c = x

    slow_resp = kslow(1:nl)*slow_c(1:nl)*dz_c(1:nl)   ! kgC/m2/s

    ! total soil respiration
    soil_resp = sum(litter_resp) + sum(fast_resp) + sum(slow_resp) ! kgC/m2/s
    
    ! soil respiration for each layer to be used for heat generated by decomposition
    soil_resp_l = (litter_resp + fast_resp + slow_resp) *rdz_c(1:nl) ! kgC/m3/s

    ! total soil carbon, kgC/m2, excluding burial layer
    soil_c_tot = sum(litter_c(1:nl)*dz_c(1:nl)) + sum(fast_c(1:nl)*dz_c(1:nl)) + sum(slow_c(1:nl)*dz_c(1:nl))

    ! methane emissions from wetland
    ch4_emis_wetland = sum(ch4_frac_wet * soil_resp_l*dz_c(1:nl)) * ch4_par%c_ch4_conv ! kgCH4/m2/s

    ! carbon conservation check
    if (check_carbon) then
      carbon_bal_soil = sum((litter_c-litter_c_old)*dz_c(1:nlc)) &
        + sum((fast_c-fast_c_old)*dz_c(1:nlc)) &
        + sum((slow_c-slow_c_old)*dz_c(1:nlc)) &
        - sum(litterfall)*dt_c &
        + soil_resp*dt_c
      if( abs(carbon_bal_soil) .gt. 1.d-10) then
        print *,'soil carbon balance',carbon_bal_soil
        stop
      endif
    endif

!......C13...............

    ! LITTER

    ! top soil layer
    k = 1
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - klitter(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - adv(k)*rdz_c(k)
    c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
    r(k) = -litter_c13(k)/dt_c - litterfall13(k)*rdz_c(k)  ! kgC/m3/s

    ! intermediate layers
    do k=2,nlc
     a(k) = diff(k-1)*rdz_neg_c(k)*rdz_c(k) + adv(k-1)*rdz_c(k)
     b(k) = -1._wp/dt_c - klitter(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - diff(k-1)*rdz_neg_c(k)*rdz_c(k) - adv(k)*rdz_c(k) 
     c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
     r(k) = -litter_c13(k)/dt_c - litterfall13(k)*rdz_c(k)
    enddo

    ! solve tridiagonal system
    call tridiag_solve(a,b,c,r,x,nlc)
    ! assign new litter carbon
    litter_c13 = x

    ! litter respiration to atmosphere
    litter_resp13 = soilc_par%f_resp_litter*klitter(1:nl)*litter_c13(1:nl)*dz_c(1:nl)   ! kgC/m2/s

    ! FAST SOIL CARBON

    ! top soil layer
    k = 1
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - kfast(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - adv(k)*rdz_c(k)
    c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
    r(k) = -fast_c13(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_fast*klitter(k)*litter_c13(k)  ! kgC/m3/s

    ! intermediate layers
    do k=2,nlc
     a(k) = diff(k-1)*rdz_neg_c(k)*rdz_c(k) + adv(k-1)*rdz_c(k)
     b(k) = -1._wp/dt_c - kfast(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - diff(k-1)*rdz_neg_c(k)*rdz_c(k) - adv(k)*rdz_c(k) 
     c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
     r(k) = -fast_c13(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_fast*klitter(k)*litter_c13(k)
    enddo

    ! solve tridiagonal system
    call tridiag_solve(a,b,c,r,x,nlc)
    ! assign new litter carbon
    fast_c13 = x

    fast_resp13 = kfast(1:nl)*fast_c13(1:nl)*dz_c(1:nl)  ! kgC/m2/s

    ! SLOW SOIL CARBON

    ! top soil layer
    k = 1
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - kslow(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - adv(k)*rdz_c(k)
    c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
    r(k) = -slow_c13(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_slow*klitter(k)*litter_c13(k)  ! kgC/m3/s

    ! intermediate layers
    do k=2,nlc
     a(k) = diff(k-1)*rdz_neg_c(k)*rdz_c(k) + adv(k-1)*rdz_c(k)
     b(k) = -1._wp/dt_c - kslow(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - diff(k-1)*rdz_neg_c(k)*rdz_c(k) - adv(k)*rdz_c(k) 
     c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k)
     r(k) = -slow_c13(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_slow*klitter(k)*litter_c13(k)
    enddo

    ! solve tridiagonal system
    call tridiag_solve(a,b,c,r,x,nlc)
    ! assign new litter carbon
    slow_c13 = x

    slow_resp13 = kslow(1:nl)*slow_c13(1:nl)*dz_c(1:nl)   ! kgC/m2/s

    ! total soil respiration
    soil_resp13 = sum(litter_resp13) + sum(fast_resp13) + sum(slow_resp13) ! kgC/m2/s

    ! total soil carbon, kgC/m2
    soil_c13_tot = sum(litter_c13(1:nl)*dz_c(1:nl)) + sum(fast_c13(1:nl)*dz_c(1:nl)) + sum(slow_c13(1:nl)*dz_c(1:nl))

    ! methane emissions from wetland
    c13h4_emis_wetland = sum(ch4_frac_wet * (litter_resp13+fast_resp13+slow_resp13)) * ch4_par%c_ch4_conv ! kgCH4/m2/s

    ! carbon conservation check
    if (check_carbon) then
      carbon13_bal_soil = sum((litter_c13-litter_c13_old)*dz_c(1:nlc)) &
        + sum((fast_c13-fast_c13_old)*dz_c(1:nlc)) &
        + sum((slow_c13-slow_c13_old)*dz_c(1:nlc)) &
        - sum(litterfall13)*dt_c &
        + soil_resp13*dt_c
      if( abs(carbon13_bal_soil) .gt. 1.d-10) then
        print *,'soil carbon 13 balance',carbon13_bal_soil
        stop
      endif
    endif

      
!.....C14........

    ! LITTER

    ! top soil layer
    k = 1
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - klitter(k) - c14_tdec - diff(k)*rdz_pos_c(k)*rdz_c(k) - adv(k)*rdz_c(k)
    c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
    r(k) = -litter_c14(k)/dt_c - litterfall14(k)*rdz_c(k)  ! kgC/m3/s

    ! intermediate layers
    do k=2,nlc
     a(k) = diff(k-1)*rdz_neg_c(k)*rdz_c(k) + adv(k-1)*rdz_c(k)
     b(k) = -1._wp/dt_c - klitter(k) - c14_tdec - diff(k)*rdz_pos_c(k)*rdz_c(k) - diff(k-1)*rdz_neg_c(k)*rdz_c(k) - adv(k)*rdz_c(k) 
     c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k)
     r(k) = -litter_c14(k)/dt_c - litterfall14(k)*rdz_c(k)
    enddo

    ! solve tridiagonal system
    call tridiag_solve(a,b,c,r,x,nlc)
    ! assign new litter carbon
    litter_c14 = x

    ! litter respiration to atmosphere
    litter_resp14 = soilc_par%f_resp_litter*klitter(1:nl)*litter_c14(1:nl)*dz_c(1:nl)   ! kgC/m2/s


    ! FAST SOIL CARBON

    ! top soil layer
    k = 1
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - kfast(k) - c14_tdec - diff(k)*rdz_pos_c(k)*rdz_c(k) - adv(k)*rdz_c(k)
    c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
    r(k) = -fast_c14(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_fast*klitter(k)*litter_c14(k)  ! kgC/m3/s

    ! intermediate layers
    do k=2,nlc
     a(k) = diff(k-1)*rdz_neg_c(k)*rdz_c(k) + adv(k-1)*rdz_c(k)
     b(k) = -1._wp/dt_c - kfast(k) - c14_tdec - diff(k)*rdz_pos_c(k)*rdz_c(k) - diff(k-1)*rdz_neg_c(k)*rdz_c(k) - adv(k)*rdz_c(k) 
     c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
     r(k) = -fast_c14(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_fast*klitter(k)*litter_c14(k)
    enddo

    ! solve tridiagonal system
    call tridiag_solve(a,b,c,r,x,nlc)
    ! assign new litter carbon
    fast_c14 = x

    fast_resp14 = kfast(1:nl)*fast_c14(1:nl)*dz_c(1:nl)  ! kgC/m2/s


    ! SLOW SOIL CARBON

    ! top soil layer
    k = 1
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - kslow(k) - c14_tdec - diff(k)*rdz_pos_c(k)*rdz_c(k) - adv(k)*rdz_c(k)
    c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
    r(k) = -slow_c14(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_slow*klitter(k)*litter_c14(k)  ! kgC/m3/s

    ! intermediate layers
    do k=2,nlc
     a(k) = diff(k-1)*rdz_neg_c(k)*rdz_c(k) + adv(k-1)*rdz_c(k)
     b(k) = -1._wp/dt_c - kslow(k) - c14_tdec - diff(k)*rdz_pos_c(k)*rdz_c(k) - diff(k-1)*rdz_neg_c(k)*rdz_c(k) - adv(k)*rdz_c(k) 
     c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
     r(k) = -slow_c14(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_slow*klitter(k)*litter_c14(k)
    enddo

    ! solve tridiagonal system
    call tridiag_solve(a,b,c,r,x,nlc)
    ! assign new litter carbon
    slow_c14 = x

    slow_resp14 = kslow(1:nl)*slow_c14(1:nl)*dz_c(1:nl)   ! kgC/m2/s


    ! total soil respiration
    soil_resp14 = sum(litter_resp14) + sum(fast_resp14) + sum(slow_resp14) ! kgC/m2/s

    ! total soil carbon, kgC/m2
    soil_c14_tot = sum(litter_c14(1:nl)*dz_c(1:nl)) + sum(fast_c14(1:nl)*dz_c(1:nl)) + sum(slow_c14(1:nl)*dz_c(1:nl))

    ! carbon conservation check
    if (check_carbon) then
      carbon14_bal_soil = sum((litter_c14-litter_c14_old)*dz_c(1:nlc)) &
        + sum((fast_c14-fast_c14_old)*dz_c(1:nlc)) &
        + sum((slow_c14-slow_c14_old)*dz_c(1:nlc)) &
        - sum(litterfall14)*dt_c &
        + soil_resp14*dt_c &
        + c14_tdec*soil_c14_tot*dt_c
      if( carbon14_bal_soil .gt. 1.d-20) then
        print *,'soil carbon 14 balance',carbon14_bal_soil
        stop
      endif
    endif


  end subroutine soil_carbon

end module soil_carbon_mod

