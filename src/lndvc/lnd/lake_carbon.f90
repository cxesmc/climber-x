!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : l a k e _ c a r b o n _ m o d
!
!  Purpose : carbon in soil below lakes
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
module lndvc_lake_carbon_mod

  use precision, only : wp
  use timer, only : day_mon
  use constants, only : c14_tdec
  use control, only : check_carbon
  use lnd_grid, only : dz_c, rdz_c, rdz_pos_c, rdz_neg_c
  use lnd_grid, only : nl, nlc, ncarb, ic_lake
  use lnd_params, only : dt_c, dt_day_c
  use lnd_params, only : soilc_par, ch4_par
  use tridiag, only : tridiag_solve

  implicit none

  private
  public :: lake_carbon

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l a k e _ c a r b o n
  !   Purpose    :  compute carbon evolution below lake water
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lake_carbon(litterfall, litterfall13, litterfall14, ch4_frac_lake, &
                        litter_c_lake,fast_c_lake,slow_c_lake, &
                        litter_c13_lake,fast_c13_lake,slow_c13_lake,litter_c14_lake,fast_c14_lake,slow_c14_lake, &
                        k_litter_lake,k_fast_lake,k_slow_lake,diff_lakec,adv_lakec, &
                        soil_resp,soil_resp13,soil_resp14,soil_resp_l,soil_c_tot,soil_c13_tot,soil_c14_tot, &
                        ch4_emis_lake,c13h4_emis_lake,carbon_bal_soil,carbon13_bal_soil,carbon14_bal_soil)

    implicit none

    real(wp), dimension(:), intent(in) :: litterfall, litterfall13, litterfall14
    real(wp), dimension(:), intent(in) :: ch4_frac_lake
    real(wp), dimension(:), intent(inout) :: litter_c_lake, fast_c_lake, slow_c_lake
    real(wp), dimension(:), intent(inout) :: litter_c13_lake, fast_c13_lake, slow_c13_lake, litter_c14_lake, fast_c14_lake, slow_c14_lake
    real(wp), dimension(:), intent(inout) :: k_litter_lake, k_fast_lake, k_slow_lake, diff_lakec, adv_lakec
    real(wp), intent(inout) :: soil_resp, soil_resp13, soil_resp14
    real(wp), dimension(:), intent(inout) :: soil_resp_l
    real(wp), intent(out) :: soil_c_tot, soil_c13_tot, soil_c14_tot
    real(wp), intent(out) :: ch4_emis_lake, c13h4_emis_lake
    real(wp), intent(out) :: carbon_bal_soil, carbon13_bal_soil, carbon14_bal_soil

    integer :: k
    real(wp), dimension(nlc) :: a, b, c, r, x
    real(wp), dimension(nlc) :: litter_c_lake_old, fast_c_lake_old, slow_c_lake_old
    real(wp), dimension(nl) :: litter_resp, fast_resp, slow_resp
    real(wp), dimension(nlc) :: litter_c13_lake_old, fast_c13_lake_old, slow_c13_lake_old
    real(wp), dimension(nl) :: litter_resp13, fast_resp13, slow_resp13
    real(wp), dimension(nlc) :: litter_c14_lake_old, fast_c14_lake_old, slow_c14_lake_old
    real(wp), dimension(nl) :: litter_resp14, fast_resp14, slow_resp14
    real(wp), dimension(nlc) :: klitter, kfast, kslow, diff, adv


    ! save old prognostic variables
    litter_c_lake_old = litter_c_lake
    fast_c_lake_old = fast_c_lake
    slow_c_lake_old = slow_c_lake
    litter_c13_lake_old = litter_c13_lake
    fast_c13_lake_old = fast_c13_lake
    slow_c13_lake_old = slow_c13_lake
    litter_c14_lake_old = litter_c14_lake
    fast_c14_lake_old = fast_c14_lake
    slow_c14_lake_old = slow_c14_lake

    ! decomposition rate, 1/s
    klitter = k_litter_lake / real(dt_day_c,wp)*day_mon
    kfast   = k_fast_lake / real(dt_day_c,wp)*day_mon
    kslow   = k_slow_lake / real(dt_day_c,wp)*day_mon
    ! vertical carbon diffusivity
    diff = diff_lakec / real(dt_day_c,wp)*day_mon
    ! vertical carbon advection
    adv = adv_lakec / real(dt_day_c,wp)*day_mon

    ! reset cumulated values
    k_litter_lake = 0._wp
    k_fast_lake = 0._wp
    k_slow_lake = 0._wp
    diff_lakec = 0._wp
    adv_lakec = 0._wp

    ! LITTER

    ! top soil layer
    k = 1
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - klitter(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - adv(k)*rdz_c(k)
    c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
    r(k) = -litter_c_lake(k)/dt_c - litterfall(k)*rdz_c(k)  ! kgC/m3/s

    do k=2,nlc
     a(k) = diff(k-1)*rdz_neg_c(k)*rdz_c(k) + adv(k-1)*rdz_c(k)
     b(k) = -1._wp/dt_c - klitter(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - diff(k-1)*rdz_neg_c(k)*rdz_c(k) - adv(k)*rdz_c(k) 
     c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k)
     r(k) = -litter_c_lake(k)/dt_c - litterfall(k)*rdz_c(k)  ! kgC/m3/s
    enddo

    ! solve tridiagonal system
    call tridiag_solve(a,b,c,r,x,nlc)
    ! assign new litter carbon
    litter_c_lake = x

    ! litter respiration to atmosphere
    litter_resp = soilc_par%f_resp_litter*klitter(1:nl)*litter_c_lake(1:nl)*dz_c(1:nl)   ! kgC/m2/s

    ! FAST SOIL CARBON

    ! top soil layer
    k = 1
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - kfast(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - adv(k)*rdz_c(k)
    c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
    r(k) = -fast_c_lake(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_fast*klitter(k)*litter_c_lake(k)  ! kgC/m3/s

    do k=2,nlc
     a(k) = diff(k-1)*rdz_neg_c(k)*rdz_c(k) + adv(k-1)*rdz_c(k)
     b(k) = -1._wp/dt_c - kfast(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - diff(k-1)*rdz_neg_c(k)*rdz_c(k) - adv(k)*rdz_c(k) 
     c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
     r(k) = -fast_c_lake(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_fast*klitter(k)*litter_c_lake(k)
    enddo

    ! solve tridiagonal system
    call tridiag_solve(a,b,c,r,x,nlc)
    ! assign new litter carbon
    fast_c_lake = x

    fast_resp = kfast(1:nl)*fast_c_lake(1:nl)*dz_c(1:nl)  ! kgC/m2/s

    ! SLOW SOIL CARBON

    ! top soil layer
    k = 1
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - kslow(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - adv(k)*rdz_c(k)
    c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
    r(k) = -slow_c_lake(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_slow*klitter(k)*litter_c_lake(k)  ! kgC/m3/s

    ! intermediate layers
    do k=2,nlc
     a(k) = diff(k-1)*rdz_neg_c(k)*rdz_c(k) + adv(k-1)*rdz_c(k)
     b(k) = -1._wp/dt_c - kslow(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - diff(k-1)*rdz_neg_c(k)*rdz_c(k) - adv(k)*rdz_c(k) 
     c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
     r(k) = -slow_c_lake(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_slow*klitter(k)*litter_c_lake(k)
    enddo

    ! solve tridiagonal system
    call tridiag_solve(a,b,c,r,x,nlc)
    ! assign new litter carbon
    slow_c_lake = x

    slow_resp = kslow(1:nl)*slow_c_lake(1:nl)*dz_c(1:nl)   ! kgC/m2/s


    ! total soil respiration
    soil_resp = sum(litter_resp) + sum(fast_resp) + sum(slow_resp) ! kgC/m2/s

    ! soil respiration for each layer to be used for heat generated by decomposition
    soil_resp_l = (litter_resp + fast_resp + slow_resp) *rdz_c(1:nl) ! kgC/m3/s

    ! total soil carbon, kgC/m2
    soil_c_tot = sum(litter_c_lake(1:nl)*dz_c(1:nl)) + sum(fast_c_lake(1:nl)*dz_c(1:nl)) + sum(slow_c_lake(1:nl)*dz_c(1:nl))

    ! methane emissions from wetland
    ch4_emis_lake = sum(ch4_frac_lake * soil_resp_l*dz_c(1:nl)) * ch4_par%c_ch4_conv ! kgCH4/m2/s

    ! carbon conservation check
    if (check_carbon) then
      carbon_bal_soil = sum((litter_c_lake-litter_c_lake_old)*dz_c(1:nlc)) &
        + sum((fast_c_lake-fast_c_lake_old)*dz_c(1:nlc)) &
        + sum((slow_c_lake-slow_c_lake_old)*dz_c(1:nlc)) &
        - sum(litterfall)*dt_c &
        + soil_resp*dt_c
      if( abs(carbon_bal_soil) .gt. 1.d-10) then
        print *,'lake carbon balance',carbon_bal_soil
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
    r(k) = -litter_c13_lake(k)/dt_c - litterfall13(k)*rdz_c(k)  ! kgC/m3/s

    ! intermediate layers
    do k=2,nlc
     a(k) = diff(k-1)*rdz_neg_c(k)*rdz_c(k) + adv(k-1)*rdz_c(k)
     b(k) = -1._wp/dt_c - klitter(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - diff(k-1)*rdz_neg_c(k)*rdz_c(k) - adv(k)*rdz_c(k) 
     c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
     r(k) = -litter_c13_lake(k)/dt_c - litterfall13(k)*rdz_c(k)  ! kgC/m3/s
    enddo

    ! solve tridiagonal system
    call tridiag_solve(a,b,c,r,x,nlc)
    ! assign new litter carbon
    litter_c13_lake = x

    ! litter respiration to atmosphere
    litter_resp13 = soilc_par%f_resp_litter*klitter(1:nl)*litter_c13_lake(1:nl)*dz_c(1:nl)   ! kgC/m2/s

    ! FAST SOIL CARBON

    ! top soil layer
    k = 1
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - kfast(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - adv(k)*rdz_c(k)
    c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
    r(k) = -fast_c13_lake(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_fast*klitter(k)*litter_c13_lake(k)  ! kgC/m3/s

    do k=2,nlc
     a(k) = diff(k-1)*rdz_neg_c(k)*rdz_c(k) + adv(k-1)*rdz_c(k)
     b(k) = -1._wp/dt_c - kfast(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - diff(k-1)*rdz_neg_c(k)*rdz_c(k) - adv(k)*rdz_c(k) 
     c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k)
     r(k) = -fast_c13_lake(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_fast*klitter(k)*litter_c13_lake(k)
    enddo

    ! solve tridiagonal system
    call tridiag_solve(a,b,c,r,x,nlc)
    ! assign new litter carbon
    fast_c13_lake = x

    fast_resp13 = kfast(1:nl)*fast_c13_lake(1:nl)*dz_c(1:nl)  ! kgC/m2/s

    ! SLOW SOIL CARBON

    ! top soil layer
    k = 1
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - kslow(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - adv(k)*rdz_c(k)
    c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
    r(k) = -slow_c13_lake(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_slow*klitter(k)*litter_c13_lake(k)  ! kgC/m3/s

    do k=2,nlc
     a(k) = diff(k-1)*rdz_neg_c(k)*rdz_c(k) + adv(k-1)*rdz_c(k)
     b(k) = -1._wp/dt_c - kslow(k) - diff(k)*rdz_pos_c(k)*rdz_c(k) - diff(k-1)*rdz_neg_c(k)*rdz_c(k) - adv(k)*rdz_c(k) 
     c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
     r(k) = -slow_c13_lake(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_slow*klitter(k)*litter_c13_lake(k)
    enddo

    ! solve tridiagonal system
    call tridiag_solve(a,b,c,r,x,nlc)
    ! assign new litter carbon
    slow_c13_lake = x

    slow_resp13 = kslow(1:nl)*slow_c13_lake(1:nl)*dz_c(1:nl)   ! kgC/m2/s

    ! total soil respiration
    soil_resp13 = sum(litter_resp13) + sum(fast_resp13) + sum(slow_resp13) ! kgC/m2/s

    ! total soil carbon, kgC/m2, excluding burial layer
    soil_c13_tot = sum(litter_c13_lake(1:nl)*dz_c(1:nl)) + sum(fast_c13_lake(1:nl)*dz_c(1:nl)) + sum(slow_c13_lake(1:nl)*dz_c(1:nl))

    ! methane emissions from wetland
    c13h4_emis_lake = sum(ch4_frac_lake * (litter_resp13+fast_resp13+slow_resp13)) * ch4_par%c_ch4_conv ! kgCH4/m2/s

    ! carbon conservation check
    if (check_carbon) then
      carbon13_bal_soil = sum((litter_c13_lake-litter_c13_lake_old)*dz_c(1:nlc)) &
        + sum((fast_c13_lake-fast_c13_lake_old)*dz_c(1:nlc)) &
        + sum((slow_c13_lake-slow_c13_lake_old)*dz_c(1:nlc)) &
        - sum(litterfall13)*dt_c &
        + soil_resp13*dt_c
      if( abs(carbon13_bal_soil) .gt. 1.d-10) then
        print *,'lake carbon 13 balance',carbon13_bal_soil
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
    r(k) = -litter_c14_lake(k)/dt_c - litterfall14(k)*rdz_c(k)  ! kgC/m3/s

    do k=2,nlc
     a(k) = diff(k-1)*rdz_neg_c(k)*rdz_c(k) + adv(k-1)*rdz_c(k)
     b(k) = -1._wp/dt_c - klitter(k) - c14_tdec - diff(k)*rdz_pos_c(k)*rdz_c(k) - diff(k-1)*rdz_neg_c(k)*rdz_c(k) - adv(k)*rdz_c(k) 
     c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
     r(k) = -litter_c14_lake(k)/dt_c - litterfall14(k)*rdz_c(k)  ! kgC/m3/s
    enddo

    ! solve tridiagonal system
    call tridiag_solve(a,b,c,r,x,nlc)
    ! assign new litter carbon
    litter_c14_lake = x

    ! litter respiration to atmosphere
    litter_resp14 = soilc_par%f_resp_litter*klitter(1:nl)*litter_c14_lake(1:nl)*dz_c(1:nl)   ! kgC/m2/s

    ! FAST SOIL CARBON

    ! top soil layer
    k = 1
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - kfast(k) - c14_tdec - diff(k)*rdz_pos_c(k)*rdz_c(k) - adv(k)*rdz_c(k)
    c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
    r(k) = -fast_c14_lake(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_fast*klitter(k)*litter_c14_lake(k)  ! kgC/m3/s

    do k=2,nlc
     a(k) = diff(k-1)*rdz_neg_c(k)*rdz_c(k) + adv(k-1)*rdz_c(k)
     b(k) = -1._wp/dt_c - kfast(k) - c14_tdec - diff(k)*rdz_pos_c(k)*rdz_c(k) - diff(k-1)*rdz_neg_c(k)*rdz_c(k) - adv(k)*rdz_c(k) 
     c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
     r(k) = -fast_c14_lake(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_fast*klitter(k)*litter_c14_lake(k)
    enddo

    ! solve tridiagonal system
    call tridiag_solve(a,b,c,r,x,nlc)
    ! assign new litter carbon
    fast_c14_lake = x

    fast_resp14 = kfast(1:nl)*fast_c14_lake(1:nl)*dz_c(1:nl)  ! kgC/m2/s

    ! SLOW SOIL CARBON

    ! top soil layer
    k = 1
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - kslow(k) - c14_tdec - diff(k)*rdz_pos_c(k)*rdz_c(k) - adv(k)*rdz_c(k)
    c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
    r(k) = -slow_c14_lake(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_slow*klitter(k)*litter_c14_lake(k)  ! kgC/m3/s

    do k=2,nlc
     a(k) = diff(k-1)*rdz_neg_c(k)*rdz_c(k) + adv(k-1)*rdz_c(k)
     b(k) = -1._wp/dt_c - kslow(k) - c14_tdec - diff(k)*rdz_pos_c(k)*rdz_c(k) - diff(k-1)*rdz_neg_c(k)*rdz_c(k) - adv(k)*rdz_c(k) 
     c(k) = diff(k)*rdz_pos_c(k)*rdz_c(k) 
     r(k) = -slow_c14_lake(k)/dt_c - (1._wp-soilc_par%f_resp_litter)*soilc_par%f_litter_to_slow*klitter(k)*litter_c14_lake(k)
    enddo

    ! solve tridiagonal system
    call tridiag_solve(a,b,c,r,x,nlc)
    ! assign new litter carbon
    slow_c14_lake = x

    slow_resp14 = kslow(1:nl)*slow_c14_lake(1:nl)*dz_c(1:nl)   ! kgC/m2/s


    ! total soil respiration
    soil_resp14 = sum(litter_resp14) + sum(fast_resp14) + sum(slow_resp14) ! kgC/m2/s

    ! total soil carbon, kgC/m2
    soil_c14_tot = sum(litter_c14_lake(1:nl)*dz_c(1:nl)) + sum(fast_c14_lake(1:nl)*dz_c(1:nl)) + sum(slow_c14_lake(1:nl)*dz_c(1:nl))

    ! carbon conservation check
    if (check_carbon) then
      carbon14_bal_soil = sum((litter_c14_lake-litter_c14_lake_old)*dz_c(1:nlc)) &
        + sum((fast_c14_lake-fast_c14_lake_old)*dz_c(1:nlc)) &
        + sum((slow_c14_lake-slow_c14_lake_old)*dz_c(1:nlc)) &
        - sum(litterfall14)*dt_c &
        + soil_resp14*dt_c &
        + c14_tdec*soil_c14_tot*dt_c
      if( carbon14_bal_soil .gt. 1.d-10) then
        print *,'lake carbon 14 balance',carbon14_bal_soil
        stop
      endif
    endif


    return

  end subroutine lake_carbon

end module lndvc_lake_carbon_mod
