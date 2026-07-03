!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : e o s _ m o d
!
!  Purpose : equation of state of seawater
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Neil R. Edwards and Matteo Willeit
!
! This file is part of CLIMBER-X.
!
! This file was ported from the original c-GOLDSTEIN model,
! see Edwards and Marsh (2005)
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
module eos_mod

  use precision, only : wp
  use ocn_grid
  use ocn_params, only : i_eos
  use constants, only : rho_w

  implicit none

  real(wp), dimension(5) :: ec = [-0.0559_wp,0.7968_wp,-0.0063_wp,3.7315e-5_wp,2.5e-5_wp]

  real(wp), dimension(3) :: ai = [3.65e-4_wp,8.9e-7_wp,8.3e-5_wp]
  real(wp), dimension(8) :: bi = [999.84_wp,6.794e-2_wp,9.095e-3_wp,0.824_wp,4.09e-3_wp,5.72e-3_wp,1.02e-4_wp,4.83e-4_wp]
  real(wp), dimension(7) :: ci = [19652.2_wp,148.4_wp,54.67_wp,0.6_wp,3.24_wp,1.4e-3_wp,2.28e-3_wp]

  real(wp), dimension(6) :: a0 = [999.842594,6.793952e-2,-9.095290e-3,1.001685e-4,-1.120083e-6,6.536336e-9]
  real(wp), dimension(5) :: a1 = [8.24493e-1,-4.0899e-3,7.6438e-5,-8.2467e-7,5.3875e-9]
  real(wp), dimension(3) :: b1 = [-5.72466e-3,1.0227e-4,-1.6546e-6]
  real(wp) :: c1 = 4.8314e-4
  real(wp), dimension(7) :: c2 = [19659.33_wp,144.4304_wp,52.848_wp,0.3101089_wp,3.186519_wp,2.212276e-2_wp,6.704388e-3_wp]
  real(wp), dimension(26) :: a2 = [1.965933e4_wp,1.444304e2_wp,-1.706103_wp,9.648704e-3_wp,-4.190254e-5_wp,52.84855_wp,-0.3101089_wp,6.283263e-3_wp,-5.084188e-5_wp,3.886640e-1_wp,9.085835e-3_wp,-4.619924e-4_wp,3.186519_wp,2.212276e-2_wp,-2.984642e-4_wp,1.956415e-6_wp,6.704388e-3_wp,-1.847318e-4_wp,2.059331e-7_wp,1.480266e-4_wp,2.102898e-4_wp,-1.202016e-5_wp,1.394680e-7_wp,-2.040237e-6_wp,6.128773e-8_wp,6.207323e-10_wp]

  private
  public :: eos, eos_tb

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  e o s
  !   Purpose    :  equation of state, calculates density
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pure function eos(t,s,z) result (rho)

    implicit none

    real(wp), intent(in) :: t, s, z
    real(wp) :: rho
    real(wp) :: tc, p, pk, r00, r0, A, B, C
    real(wp) :: t2, t3, t4, s32, p2


    if (i_eos.eq.0) then

      ! state equation following WS 1993	
      rho = rho_w + ec(1)*t + ec(2)*s + ec(3)*t**2 + ec(4)*t**3  ! t in degC and s in psu; rho_w = pure-water baseline (NOT Boussinesq rho0)

    elseif (i_eos.eq.1) then

      ! Thermobaricity term (T*z) added as option. Optimised for -1<deep T<6, S=34.9
      ! but is an order of magnitude improvement even in other parts of parameter space. It does not
      ! change surface densities which are fine anyway.
      rho = rho_w + ec(1)*t + ec(2)*s + ec(3)*t**2 + ec(4)*t**3 + ec(5)*t*z  ! rho_w = pure-water baseline (NOT Boussinesq rho0)

    elseif (i_eos.eq.2) then

      ! EOS80 equation of state, Millero and Poisson 1981, but only limited number of terms! As in CLIMBER-2
      p = 0.1_wp*(-z)  ! pressure in bar
      ! is this term a correction from potential to in-situ temperature?
      tc = (t+ai(1)*p+ai(2)*p**2)/(1._wp-ai(3)*p)
      ! potential density 
      r0 = bi(1)+bi(2)*tc-bi(3)*tc*tc+s*(bi(4)-bi(5)*tc) &
        + s*sqrt(s)*(-bi(6)+bi(7)*tc)+bi(8)*s**2

      pk = ci(1)+ci(2)*tc+s*(ci(3)-ci(4)*tc)+p*(ci(5)+ci(6)*tc) &
        + p*s*ci(7)

      rho = r0/(1._wp-p/pk)

    elseif (i_eos.eq.3) then

      ! EOS80 equation of state, Millero and Poisson 1981, but only limited number of terms!
      ! parameters for potential instead of in-situ temperature (after Jackett and McDougall, 1995)

      p = 0.1_wp*(-z)  ! pressure in bar

      t2 = t**2
      t3 = t**3
      t4 = t**4

      ! potential density -
      r00 = a0(1) + a0(2)*t + a0(3)*t2 + a0(4)*t3 + a0(5)*t4 + a0(6)*t**5
      A = a1(1) + a1(2)*t + a1(3)*t2 + a1(4)*t3 + a1(5)*t4 
      B = b1(1) + b1(2)*t + b1(3)*t2
      C = c1
      r0 = r00 + A*s + B*s*sqrt(s) + C*s**2

      ! bulk secant modulus
      pk = c2(1)+c2(2)*t+s*(c2(3)-c2(4)*t)+p*(c2(5)+c2(6)*t) &
        + p*s*c2(7)

      rho = r0/(1._wp-p/pk)

    elseif (i_eos.eq.4) then

      ! EOS80 equation of state, Millero and Poisson 1981, with all terms in the bulk secant modulus!
      ! parameters for potential instead of in-situ temperature (after Jackett and McDougall, 1995)

      p = 0.1_wp*(-z)  ! pressure in bar

      t2 = t**2
      t3 = t**3
      t4 = t**4
      s32 = s*sqrt(s)
      p2 = p**2

      ! potential density -
      r00 = a0(1) + a0(2)*t + a0(3)*t2 + a0(4)*t3 + a0(5)*t4 + a0(6)*t**5
      A = a1(1) + a1(2)*t + a1(3)*t2 + a1(4)*t3 + a1(5)*t4 
      B = b1(1) + b1(2)*t + b1(3)*t2
      C = c1
      r0 = r00 + A*s + B*s32 + C*s**2

      ! bulk secant modulus
      pk = a2(1) &
        + a2(2)*t + a2(3)*t2 + a2(4)*t3 + a2(5)*t4 &
        + a2(6)*s + a2(7)*s*t + a2(8)*s*t2 + a2(9)*s*t3 &
        + a2(10)*s32 + a2(11)*s32*t + a2(12)*s32*t2 &
        + a2(13)*p + a2(14)*p*t + a2(15)*p*t2 + a2(16)*p*t3 &
        + a2(17)*p*s + a2(18)*p*s*t + a2(19)*p*s*t2 + a2(20)*p*s32 &
        + a2(21)*p2 + a2(22)*p2*t + a2(23)*p2*t2 &
        + a2(24)*p2*s + a2(25)*p2*s*t + a2(26)*p2*s*t2

      rho = r0/(1._wp-p/pk)

    elseif (i_eos.eq.5) then

      ! “stiffened” EOS derived from the compressibility of sea water and the UNESCO EOS80, following Ma 2020
      ! but only limited number of terms!
      ! parameters for theta instead of in-situ temperature (after Jackett and McDougall, 1995)

      p = 0.059808*(exp(-0.025*(-z)) - 1) + 0.100766*(-z) + 2.28405e-7*(-z)**2 !  pressure in bar

      t2 = t**2
      t3 = t**3
      t4 = t**4

      ! potential density
      r00 = a0(1) + a0(2)*t + a0(3)*t2 + a0(4)*t3 + a0(5)*t4 + a0(6)*t**5
      A = a1(1) + a1(2)*t + a1(3)*t2 + a1(4)*t3 + a1(5)*t4 
      B = b1(1) + b1(2)*t + b1(3)*t2
      C = c1
      r0 = r00 + A*s + B*s*sqrt(s) + C*s**2

      ! bulk secant modulus
      pk = c2(1)+c2(2)*t+s*(c2(3)-c2(4)*t)+p*(c2(5)+c2(6)*t) &
        + p*s*c2(7)

      rho = r0/(1._wp-p/pk)


    endif

   return

  end function eos


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  e o s _ t b
  !   Purpose    :  calculates thermobaric density from in-situ density
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine eos_tb(rho,z, &
                 rho_tb)

    implicit none

    real(wp), intent(in) :: rho, z
    real(wp), intent(out) :: rho_tb

    real(wp) :: p, rp


    p = 0.059808*(exp(-0.025*(-z)) - 1) + 0.100766*(-z) + 2.28405e-7*(-z)**2 !  pressure in bar
    rp = 1.02819 - 2.93161e-4*exp(-0.05*p) + 4.4004e-5*p 

    ! thermobaric density
    rho_tb = rho/rp


  end subroutine eos_tb

end module eos_mod

