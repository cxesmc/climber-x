!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : p e a t _ c a r b o n _ m o d
!
!  Purpose : peat carbon dynamics
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
module lndvc_peat_carbon_mod

  use precision, only : wp
  use timer, only : day_mon, day_year, time_eoy_lnd
  use constants, only : c14_tdec
  use control, only : check_carbon
  use lnd_grid, only : dz_c, rdz_c, nl, nlc, ncarb, ic_peat
  use lnd_params, only : dt_c, dt_day_c
  use lnd_params, only : soilc_par, peat_par, ch4_par 
  use tridiag, only : tridiag_solve

  implicit none

  private
  public :: peat_carbon

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  p e a t _ c a r b o n
  !   Purpose    :  compute soil carbon evolution
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine peat_carbon(f_oxic_peat,litterfall,litterfall13,litterfall14, &
                        litter_c_peat,acro_c,cato_c,litter_c13_peat,acro_c13,cato_c13,litter_c14_peat,acro_c14,cato_c14, &
                        k_litter_peat,k_acro,k_cato,k_litter_peat_anox,k_acro_anox,ch4_frac_peat, &
                        soil_resp,soil_resp13,soil_resp14,soil_resp_l,soil_c_tot,soil_c13_tot,soil_c14_tot, &
                        ch4_emis_peat,c13h4_emis_peat, &
                        carbon_bal_soil,carbon13_bal_soil,carbon14_bal_soil,peat_c_ini_year,dCpeat_dt)

    implicit none

    real(wp), intent(in) :: f_oxic_peat
    real(wp), dimension(:), intent(in) :: litterfall, litterfall13, litterfall14
    real(wp), dimension(:), intent(inout) :: cato_c, cato_c13, cato_c14
    real(wp), intent(inout) :: litter_c_peat, acro_c, litter_c13_peat, acro_c13, litter_c14_peat, acro_c14
    real(wp), intent(inout) :: k_litter_peat, k_acro, k_litter_peat_anox, k_acro_anox
    real(wp), dimension(:), intent(inout) :: k_cato
    real(wp), dimension(:), intent(inout) :: ch4_frac_peat
    real(wp), intent(inout) :: soil_resp, soil_resp13, soil_resp14
    real(wp), intent(out) :: soil_c_tot, soil_c13_tot, soil_c14_tot
    real(wp), dimension(:), intent(out) :: soil_resp_l
    real(wp), intent(out) :: ch4_emis_peat, c13h4_emis_peat
    real(wp), intent(out) :: carbon_bal_soil, carbon13_bal_soil, carbon14_bal_soil
    real(wp), intent(inout) :: peat_c_ini_year, dCpeat_dt

    integer :: k
    real(wp)                  :: peat_c, dt_p
    real(wp), dimension(1:nlc) :: a, b, c, r, x
    real(wp)                  :: litter_c_peat_old, acro_c_old
    real(wp)                  :: litter_c13_peat_old, acro_c13_old
    real(wp)                  :: litter_c14_peat_old, acro_c14_old
    real(wp), dimension(1:nlc) :: cato_c_old, cato_c13_old, cato_c14_old
    real(wp)                   :: litter_resp, acro_resp, litter_resp13, acro_resp13, litter_resp14, acro_resp14
    real(wp), dimension(1:nl) :: cato_resp, cato_resp13, cato_resp14
    real(wp)                  :: klitter, kacro, klitter_anox, kacro_anox, k_a_to_c
    real(wp), dimension(1:nlc) :: kcato
    real(wp), dimension(1:nlc) :: c_old


    ! save old prognostic variables
    litter_c_peat_old = litter_c_peat
    acro_c_old = acro_c
    cato_c_old = cato_c
    litter_c13_peat_old = litter_c13_peat
    acro_c13_old = acro_c13
    cato_c13_old = cato_c13
    litter_c14_peat_old = litter_c14_peat
    acro_c14_old = acro_c14
    cato_c14_old = cato_c14

    ! decomposition rate, 1/s
    klitter = k_litter_peat / real(dt_day_c,wp)*day_mon
    kacro   = k_acro / real(dt_day_c,wp)*day_mon
    kcato   = k_cato / real(dt_day_c,wp)*day_mon
    klitter_anox = k_litter_peat_anox / real(dt_day_c,wp)*day_mon
    kacro_anox   = k_acro_anox / real(dt_day_c,wp)*day_mon

    ! reset cumulated values
    k_litter_peat = 0._wp
    k_acro = 0._wp
    k_cato = 0._wp
    k_litter_peat_anox = 0._wp
    k_acro_anox = 0._wp


    ! LITTER

    ! new litter carbon, implicitly solved, kgC/m2, Kleinen 2012
    litter_c_peat = (litter_c_peat + sum(litterfall)*dt_c) / (1._wp+klitter*dt_c)
    litter_c13_peat = (litter_c13_peat + sum(litterfall13)*dt_c) / (1._wp+klitter*dt_c)
    litter_c14_peat = (litter_c14_peat + sum(litterfall14)*dt_c) / (1._wp+(klitter+c14_tdec)*dt_c)

    ! litter respiration to atmosphere
    litter_resp = soilc_par%f_resp_litter*klitter*litter_c_peat   ! kgC/m2/s
    litter_resp13 = soilc_par%f_resp_litter*klitter*litter_c13_peat   ! kgC/m2/s
    litter_resp14 = soilc_par%f_resp_litter*klitter*litter_c14_peat   ! kgC/m2/s


    ! ACROTELM, FAST SOIL CARBON

    ! acrotelm to catotelm transfer rate
    if( (acro_c+litter_c_peat) .gt. peat_par%acroc_crit ) then ! transfer carbon to catotelm only when acrotelm carbon > acroc_crit
     k_a_to_c = peat_par%k_acro_to_cato
    else
     k_a_to_c = 0._wp
    endif

    ! new acrotelm carbon, implicitly solved, kgC/m2, Kleinen 2012
    acro_c = (acro_c+(1._wp-soilc_par%f_resp_litter)*klitter*litter_c_peat*dt_c) / (1._wp+(k_a_to_c+kacro)*dt_c)
    acro_c13 = (acro_c13+(1._wp-soilc_par%f_resp_litter)*klitter*litter_c13_peat*dt_c) / (1._wp+(k_a_to_c+kacro)*dt_c)
    acro_c14 = (acro_c14+(1._wp-soilc_par%f_resp_litter)*klitter*litter_c14_peat*dt_c) / (1._wp+(k_a_to_c+kacro+c14_tdec)*dt_c)

    acro_resp = kacro*acro_c  ! kgC/m2/s
    acro_resp13 = kacro*acro_c13  ! kgC/m2/s
    acro_resp14 = kacro*acro_c14  ! kgC/m2/s

    ! CATOTELM, SLOW SOIL CARBON
       
    ! start from 2nd soil layer, acrotelm above
    k = 2
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - kcato(k)
    c(k) = 0._wp
    r(k) = -cato_c(k)/dt_c - k_a_to_c*acro_c*rdz_c(k)  ! kgC/m3/s

    ! intermediate layers
    do k=3,nlc-1
     a(k) = 0._wp
     b(k) = -1._wp/dt_c - kcato(k) 
     c(k) = 0._wp
     r(k) = -cato_c(k)/dt_c
    enddo

    ! bottom layer, k=nlc
    k = nlc
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - kcato(k) 
    c(k) = 0._wp
    r(k) = -cato_c(k)/dt_c

    ! solve tridiagonal system
    call tridiag_solve(a(2:nlc),b(2:nlc),c(2:nlc),r(2:nlc),x(2:nlc),nlc-1)
    ! assign new litter carbon
    cato_c(2:nlc) = x(2:nlc)

    cato_resp = kcato(1:nl)*cato_c(1:nl)*dz_c(1:nl)   ! kgC/m2/s

!......C13.....
    ! start from 2nd soil layer, acrotelm above
    k = 2
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - kcato(k)
    c(k) = 0._wp
    r(k) = -cato_c13(k)/dt_c - k_a_to_c*acro_c13*rdz_c(k)  ! kgC/m3/s

    ! intermediate layers
    do k=3,nlc-1
     a(k) = 0._wp
     b(k) = -1._wp/dt_c - kcato(k) 
     c(k) = 0._wp
     r(k) = -cato_c13(k)/dt_c
    enddo

    ! bottom layer, k=nl
    k = nlc
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - kcato(k) 
    c(k) = 0._wp
    r(k) = -cato_c13(k)/dt_c

    ! solve tridiagonal system
    call tridiag_solve(a(2:nlc),b(2:nlc),c(2:nlc),r(2:nlc),x(2:nlc),nlc-1)
    ! assign new litter carbon
    cato_c13(2:nlc) = x(2:nlc)

    cato_resp13 = kcato(1:nl)*cato_c13(1:nl)*dz_c(1:nl)   ! kgC/m2/s


!......C14.....
    ! start from 2nd soil layer, acrotelm above
    k = 2
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - kcato(k) - c14_tdec 
    c(k) = 0._wp
    r(k) = -cato_c14(k)/dt_c - k_a_to_c*acro_c14*rdz_c(k)  ! kgC/m3/s

    ! intermediate layers
    do k=3,nlc-1
     a(k) = 0._wp
     b(k) = -1._wp/dt_c - kcato(k) - c14_tdec 
     c(k) = 0._wp
     r(k) = -cato_c14(k)/dt_c
    enddo

    ! bottom layer, k=nl
    k = nlc
    a(k) = 0._wp
    b(k) = -1._wp/dt_c - kcato(k) - c14_tdec 
    c(k) = 0._wp
    r(k) = -cato_c14(k)/dt_c

    ! solve tridiagonal system
    call tridiag_solve(a(2:nlc),b(2:nlc),c(2:nlc),r(2:nlc),x(2:nlc),nlc-1)
    ! assign new litter carbon
    cato_c14(2:nlc) = x(2:nlc)

    cato_resp14 = kcato(1:nl)*cato_c14(1:nl)*dz_c(1:nl)   ! kgC/m2/s

    ! increase catotelm thickness if catotelm carbon density larger than threshhold
    do k=2,nlc-1
     if( cato_c(k) .gt. peat_par%rho_cato ) then
      c_old(k) = cato_c(k)
      cato_c(k+1)   = cato_c(k+1) + (cato_c(k) - peat_par%rho_cato) * dz_c(k) *rdz_c(k+1)
      cato_c(k) = peat_par%rho_cato
      cato_c13(k+1) = cato_c13(k+1) + cato_c13(k)*(1._wp-cato_c(k)/c_old(k)) * dz_c(k) *rdz_c(k+1)
      cato_c13(k) = cato_c13(k)*cato_c(k)/c_old(k)
      cato_c14(k+1) = cato_c14(k+1) + cato_c14(k)*(1._wp-cato_c(k)/c_old(k)) * dz_c(k) *rdz_c(k+1)
      cato_c14(k) = cato_c14(k)*cato_c(k)/c_old(k)
     endif
    enddo

    ! total soil respiration
    soil_resp   = litter_resp + acro_resp + sum(cato_resp) ! kgC/m2/s
    soil_resp13 = litter_resp13 + acro_resp13 + sum(cato_resp13) ! kgC/m2/s
    soil_resp14 = litter_resp14 + acro_resp14 + sum(cato_resp14) ! kgC/m2/s

    ! soil respiration for each layer to be used for heat generated by decomposition
    soil_resp_l(1) = (litter_resp + acro_resp) *rdz_c(1)  ! kgC/m3/s
    soil_resp_l(2:nl) = cato_resp(2:nl) *rdz_c(2:nl) ! kgC/m3/s

    ! total soil carbon, kgC/m2
    soil_c_tot = litter_c_peat + acro_c + sum(cato_c(1:nl)*dz_c(1:nl))
    soil_c13_tot = litter_c13_peat + acro_c13 + sum(cato_c13(1:nl)*dz_c(1:nl))
    soil_c14_tot = litter_c14_peat + acro_c14 + sum(cato_c14(1:nl)*dz_c(1:nl))

    ! methane emissions from peatland
    ch4_emis_peat = (ch4_frac_peat(1)*(1._wp-f_oxic_peat)*(litter_resp+acro_resp) + sum(ch4_frac_peat*cato_resp)) * ch4_par%c_ch4_conv ! kgCH4/m2/s
    c13h4_emis_peat = (ch4_frac_peat(1)*(1._wp-f_oxic_peat)*(litter_resp13+acro_resp13) + sum(ch4_frac_peat*cato_resp13)) * ch4_par%c_ch4_conv ! kgC13H4/m2/s


    ! check carbon conservation
    if (check_carbon) then
      carbon_bal_soil = (litter_c_peat-litter_c_peat_old) &
        + (acro_c-acro_c_old) &
        + sum((cato_c-cato_c_old)*dz_c(1:nlc)) &
        - sum(litterfall)*dt_c &
        + soil_resp*dt_c

      carbon13_bal_soil = (litter_c13_peat-litter_c13_peat_old) &
        + (acro_c13-acro_c13_old) &
        + sum((cato_c13-cato_c13_old)*dz_c(1:nlc)) &
        - sum(litterfall13)*dt_c &
        + soil_resp13*dt_c

      carbon14_bal_soil = (litter_c14_peat-litter_c14_peat_old) &
        + (acro_c14-acro_c14_old) &
        + sum((cato_c14-cato_c14_old)*dz_c(1:nlc)) &
        - sum(litterfall14)*dt_c &
        + soil_resp14*dt_c &
        + c14_tdec*soil_c14_tot*dt_c

      if( abs(carbon_bal_soil) .gt. 1.d-10) then
        print *,'peat carbon balance',ic_peat,carbon_bal_soil
        stop
      endif
      if( abs(carbon13_bal_soil) .gt. 1.d-10) then
        print *,'peat carbon 13 balance',ic_peat,carbon13_bal_soil
        stop
      endif
      if( carbon14_bal_soil .gt. 1.d-10) then
        print *,'peat carbon 14 balance',ic_peat,carbon14_bal_soil
        stop
      endif
    endif


    ! compute peat carbon yearly accumulation rate
    if( time_eoy_lnd ) then

     dt_p = dt_c*day_year/real(dt_day_c,wp)

     peat_c = litter_c_peat+acro_c+sum(cato_c(1:nlc)*dz_c(1:nlc)) ! kg/m2

     dCpeat_dt = (peat_c-peat_c_ini_year)/dt_p ! kg/m2/s

    endif


    return

  end subroutine peat_carbon

end module lndvc_peat_carbon_mod
