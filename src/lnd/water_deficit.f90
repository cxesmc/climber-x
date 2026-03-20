!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : w a t e r _ d e f i c i t _ m o d
!
!  Purpose : compute cumulative water deficit 
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Da Nian and Matteo Willeit
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
module water_deficit_mod

  use precision, only : wp
  use timer, only : mon, sec_day
  use constants, only : T0, Le, Rd, cap_a, rho_a
  use constants, only : e_sat_w, dqsat_dT_w, q_to_e
  use lnd_grid, only : nsurf

  implicit none

  private
  public :: calculate_pet, calculate_cwd

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  !  Subroutine : c a l c u l a t e _ P E T
  !
  !  Purpose : Calculate potential evapotranspiration using FAO-56
  !            Penman-Monteith method
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine calculate_pet(frac_surf, t2m, q2m, ps, swnet, lwd, lwu, ra, pet)

    implicit none

    real(wp), intent(in) :: frac_surf(:)      ! surface type fractions [1]
    real(wp), intent(in) :: t2m(:)      ! 2m temperature [K]
    real(wp), intent(in) :: q2m(:)      ! 2m specific humidity [kg/kg]
    real(wp), intent(in) :: ps(:)       ! Surface pressure [Pa]
    real(wp), intent(in) :: swnet(:)    ! Net SW radiation [W/m²]
    real(wp), intent(in) :: lwd(:)      ! Downward LW radiation [W/m²]
    real(wp), intent(in) :: lwu(:)      ! Upward LW radiation [W/m²]
    real(wp), intent(in) :: ra(:)       ! aerodynamic resistance [s/m]
    real(wp), intent(out) :: pet(:)     ! Potential ET [kg/m2/day]

    ! Local variables
    integer :: n
    real(wp) :: T_K, P, es, ea, VPD, dqsat_dT, Delta, gamma
    real(wp) :: Rn, G, num, denom

    ! Constants
    real(wp), parameter :: eps = 0.622_wp         ! Ratio molecular weights

    do n = 1, nsurf

      if (frac_surf(n).gt.0._wp) then

        T_K = t2m(n)

        ! Skip calculation if very cold (frozen conditions)
        if (T_K < (T0 - 10.0_wp)) then
          pet(n) = 0.0_wp
          cycle
        end if

        P = ps(n)

        ! ================================================================
        ! Vapor pressures (same as photosynthesis)
        ! ================================================================
        es = e_sat_w(T_K)                    ! Saturation VP [Pa]
        ea = q_to_e(q2m(n), P)               ! Actual VP [Pa]
        VPD = max(es - ea, 0.1_wp)           ! Vapor pressure deficit [Pa]

        ! ================================================================
        ! Slope of saturation curve 
        ! ================================================================
        dqsat_dT = dqsat_dT_w(T_K, P)        ! d(qsat)/dT [kg/kg/K]
        Delta = (P / eps) * dqsat_dT         ! d(es)/dT [Pa/K]

        ! ================================================================
        ! Psychrometric constant
        ! ================================================================
        gamma = (cap_a * P) / (eps * Le)    ! Psychrometric constant [Pa/K]

        ! ================================================================
        ! Net radiation
        ! ================================================================
        Rn = swnet(n) + (lwd(n) - lwu(n))  ! [W/m²]

        ! ═══════════════════════════════════════════════════════
        ! Soil heat flux [W/m²] 
        ! ═══════════════════════════════════════════════════════
        G = 0.0_wp  ! Daily: G ≈ 0

        ! ================================================================
        ! FAO-56 Penman-Monteith equation
        ! ================================================================
        ! λE = [Δ(Rn−G) + ρcₚ·VPD/rₐ] / [Δ + γ(1+rₛ/rₐ)]. PET when rₛ=0, no surface resistance
        num = Delta * (Rn - G) + rho_a(T_K,P)*cap_a*VPD/ra(n)
        denom = Delta + gamma

        if (denom > 0.0_wp) then
          pet(n) = (num / denom) / Le   ! kg/m2/s
        else
          pet(n) = 0.0_wp
        end if

        ! Ensure non-negative
        pet(n) = max(0.0_wp, pet(n))

      else 

        pet(n) = 0._wp

      endif

    end do

    return

  end subroutine calculate_pet


  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  !  Subroutine : c a l c u l a t e _ C W D
  !
  !  Purpose : Calculate cumulative water deficit 
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine calculate_cwd(frac_surf, rain, snow, pet, cwd_mon)

    implicit none

    real(wp), intent(in) :: frac_surf(:)      ! surface type fractions [1]
    real(wp), intent(in) :: rain(:)           ! rainfall rate [kg/m2/s]
    real(wp), intent(in) :: snow(:)           ! snowfall rate [kg/m2/s]
    real(wp), intent(in) :: pet(:)            ! potential evapotranspiration [kg/m2/s]
    real(wp), intent(inout) :: cwd_mon(:)        ! cumulative water deficit [kg/m2]

    integer :: n
    real(wp) :: prc_daily, pet_daily

    prc_daily = 0.0_wp
    pet_daily = 0.0_wp

    ! Calculate grid-cell daily values
    do n=1,nsurf
      if (frac_surf(n) > 0.0_wp) then
        ! Precipitation [kg/m2/day]
        prc_daily = prc_daily + (rain(n) + snow(n)) * sec_day * frac_surf(n)
        ! PET already in mm/day
        pet_daily = pet_daily + pet(n) * sec_day * frac_surf(n)
      end if
    end do

    ! Accumulate monthly totals [kg/m2 or mm]
    cwd_mon(mon) = max(0._wp, cwd_mon(mon) + pet_daily - prc_daily)

    return

  end subroutine calculate_cwd

end module water_deficit_mod
