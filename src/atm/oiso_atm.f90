!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : o i s o _ a t m _ m o d
!
!  Purpose : compute fractionation of oxygen isotopes in the atmosphere
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2026 Potsdam Institute for Climate Impact Research,
!                         Andrey Ganopolski, Matteo Willeit and
!                         Julius Eberhard
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
module oiso_atm_mod

  use atm_params, only : wp
  use atm_params, only : coef_fro18_sv_kin
  use constants, only : T0   ! 0 degC in Kelvin
  !$ use omp_lib

  implicit none

  private
  public :: oiso_fract_prcw, oiso_fract_prcs, d18o_atm
  public :: Rstd_o16, Rstd_o18

  ! reference O16/O18 water-mass ratios relative to total water (VSMOW).
  real(wp), parameter :: Rstd_o16 = 1._wp           ! O16 water / total water
  real(wp), parameter :: Rstd_o18 = 2.0052e-3_wp    ! O18 water / total water (VSMOW, 2005.20 ppm)

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  o i s o _ f r a c t _ p r c w
  !   Purpose    :  compute fractions of O16 and O18 in rain
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine oiso_fract_prcw(teff, &                  ! in
                             fro16_prcw, fro18_prcw)  ! out

    ! Note: teff needs to be in K.

    implicit none

    real(wp), intent(in ) :: teff

    real(wp), intent(out) :: fro16_prcw, fro18_prcw

    !-------------------------------
    ! compute equilibrium fractionation
    ! factor at temperature teff
    ! following Roche (2013) and Majoube (1971b)
    !-------------------------------

    ! TODO(2026-05-28,eberhard): Decide whether to compute actual O16 fractionation
    !                            or to leave it at alpha16 = 1 (O16 = total water).
    fro18_prcw = exp(1137._wp/teff**2 - 0.4156_wp/teff - 0.0020667_wp)  ! dimensionless alpha18
    !fro16_prcw = 1._wp - fro18_prcw
    fro16_prcw = 1._wp                                                  ! dimensionless alpha16

    return

  end subroutine oiso_fract_prcw


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  o i s o _ f r a c t _ p r c s
  !   Purpose    :  compute fractions of O16 and O18 in snow
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine oiso_fract_prcs(teff, &                  ! in
                             fro16_prcs, fro18_prcs)  ! out

    ! Note: teff needs to be in degC.

    implicit none

    real(wp), intent(in ) :: teff

    real(wp), intent(out) :: fro16_prcs, fro18_prcs

    real(wp) :: fro18_sv, D18_D

    !-------------------------------
    ! compute kinetic fractionation
    ! factor at temperature teff
    ! following Roche (2013), Majoube (1971a)
    ! and Merlivat & Jouzel (1979)
    !-------------------------------

    ! equilibrium fractionation coefficient between water vapor and solid water
    ! following Majoube (1971a)
    fro18_sv = exp(11.839_wp/teff - 0.028244_wp)  ! in units of 1

    ! ratio of molecular diffusivities of O18 water and average water
    ! following Merlivat (1978)
    D18_D = 0.9723_wp

    ! TODO(2026-05-28,eberhard): Decide whether to compute actual O16 fractionation
    !                            or to leave it at alpha16 = 1 (O16 = total water).
    fro18_prcs = fro18_sv*(1._wp-coef_fro18_sv_kin*(teff-T0)) &
      / (1._wp+fro18_sv*(-coef_fro18_sv_kin*(teff-T0))/D18_D)  ! dimensionless alpha18
    !fro16_prcs = 1._wp - fro18_prcs
    fro16_prcs = 1._wp                                    ! dimensionless alpha16

    return

  end subroutine oiso_fract_prcs


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  d 1 8 o _ a t m
  !   Purpose    :  compute deviation of O18/O16 compared to
  !                 Vienna Standard Ocean Mean Water (VSMOW)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine d18o_atm(o16, o18, &  ! in
                      d18o)        ! out

    implicit none

    real(wp), intent(in ) :: o16, o18

    real(wp), intent(out) :: d18o

    d18o = 1000.*o18/o16 / Rstd_o18 - 1000._wp  ! permil; VSMOW mean value = 2005.20 ppm

    return

  end subroutine d18o_atm

end module oiso_atm_mod

