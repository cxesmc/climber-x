!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module : l n d v c _ i n i t _ c e l l _ m o d
!
!  Purpose : one-time physical initialization of a vegetated land virtual cell
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
module lndvc_init_cell_mod

  use precision, only : wp
  use constants, only : rho_w, rho_i, T0
  use lnd_grid, only : dz, nl, npft, i_bare
  use lnd_params, only : i_init_veg, veg_par, pft_par, peat_par, surf_par

  implicit none

  private
  public :: lndvc_init_cell_veg

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l n d v c _ i n i t _ c e l l _ v e g
  !   Purpose    :  initialize the prognostic veg/soil state of a vegetated
  !                 land virtual cell. Faithful de-scoped port of the reference
  !                 init_cell_veg (src/lnd/init_cell.f90): the shelf-restart and
  !                 cross-class skin-temperature-mean branches do not apply to a
  !                 single-class veg vc, so this uses the reference cold-start
  !                 path (t_skin_mean = T0). Soil parameters are assumed already
  !                 seeded (theta_sat intent(in)); the soil carbon pools are C.4
  !                 and left to their neutral init. Surface arrays are sized over
  !                 the nveg (bare+PFT) sub-tiles, veg arrays over npft.
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lndvc_init_cell_veg(c13_c12_atm, c14_c_atm, f_veg, &
                          pft_frac, t_skin, t_soil, w_can, s_can, &
                          w_snow, w_snow_max, h_snow, mask_snow, &
                          theta, theta_w, theta_i, w_w, w_i, theta_sat, &
                          alt, gdd5, gdd, phen, phen_acc, lai_bal, lai, sai, root_frac, litter_in_frac, &
                          gamma_dist, gamma_fire, npp_ann, npp13_ann, npp14_ann, &
                          leaf_c, root_c, stem_c, veg_h, veg_c, veg_c13, veg_c14, &
                          f_peat, f_peat_pot, w_table_min, w_table_peat, dCpeat_dt, z0m)

    implicit none

    real(wp), intent(in) :: c13_c12_atm, c14_c_atm, f_veg
    real(wp), dimension(:), intent(inout) :: pft_frac         ! npft
    real(wp), dimension(:), intent(inout) :: t_skin           ! nveg
    real(wp), dimension(0:), intent(inout) :: t_soil          ! 0:nl
    real(wp), dimension(:), intent(inout) :: w_can, s_can     ! nveg
    real(wp), intent(inout) :: w_snow, w_snow_max, h_snow
    integer,  intent(inout) :: mask_snow
    real(wp), dimension(:), intent(inout) :: theta, theta_w, theta_i, w_w, w_i   ! nl
    real(wp), dimension(:), intent(in)    :: theta_sat        ! nl
    real(wp), intent(inout) :: alt, gdd5
    real(wp), dimension(:), intent(inout) :: gdd, phen, phen_acc, lai_bal, lai, sai, gamma_dist, gamma_fire  ! npft
    real(wp), dimension(:,:), intent(inout) :: root_frac      ! nl,npft
    real(wp), dimension(:), intent(inout) :: litter_in_frac   ! nl
    real(wp), dimension(:), intent(inout) :: npp_ann, npp13_ann, npp14_ann        ! npft
    real(wp), dimension(:), intent(inout) :: leaf_c, root_c, stem_c, veg_h, veg_c, veg_c13, veg_c14  ! npft
    real(wp), intent(inout) :: f_peat, f_peat_pot, w_table_min, w_table_peat, dCpeat_dt
    real(wp), dimension(:), intent(inout) :: z0m              ! nveg

    integer :: n
    real(wp) :: t_skin_mean


    ! initialise prognostic variables
    alt = -1._wp

    ! vegetation properties
    gdd5     = 5000._wp
    gdd      = 0._wp
    phen     = 0._wp
    phen_acc = 0._wp
    if (i_init_veg.eq.1) then
      ! desert
      pft_frac = veg_par%seed_fraction
      lai_bal  = pft_par%lai_min
    else if (i_init_veg.eq.2) then
      ! forest
      pft_frac = veg_par%seed_fraction
      lai_bal  = pft_par%lai_min
      pft_frac(1:2) = 0.5_wp
      lai_bal(1:2)  = pft_par%lai_max(1:2)
    endif
    lai       = lai_bal
    sai       = lai_bal * veg_par%sai_scale
    npp_ann   = 0._wp
    npp13_ann = 0._wp
    npp14_ann = 0._wp
    leaf_c    = lai_bal/pft_par%sla
    root_c    = leaf_c
    stem_c    = pft_par%awl*lai_bal**pft_par%bwl
    veg_h     = pft_par%awh * lai_bal
    veg_c     = leaf_c + root_c + stem_c
    veg_c13   = veg_c * c13_c12_atm
    veg_c14   = veg_c * c14_c_atm
    gamma_dist = 0.001_wp
    gamma_fire = 0._wp
    root_frac  = pft_par%root_frac
    litter_in_frac = pft_par%litter_in_frac

    w_can = 0._wp
    s_can = 0._wp

    ! cold-start skin temperature: a veg vc holds no ice/lake sibling tiles, so
    ! the reference cross-class frac-weighted mean reduces to the cold-start T0
    t_skin_mean = T0
    t_skin(:) = t_skin_mean

    w_snow     = 0._wp
    w_snow_max = w_snow
    h_snow     = 0._wp
    mask_snow  = 0

    ! initialize soil temperature to skin temperature; soil saturated, phase by T
    t_soil(0)    = T0
    t_soil(1:nl) = t_skin_mean
    if (t_skin_mean .lt. T0) then
      theta_i = theta_sat
      theta_w = 0._wp
    else
      theta_i = 0._wp
      theta_w = theta_sat
    endif
    w_w   = theta_w * dz(1:nl) * rho_w
    w_i   = theta_i * dz(1:nl) * rho_i
    theta = w_w/(rho_w*dz(1:nl)) + w_i/(rho_i*dz(1:nl))

    ! initialize peatlands
    if (peat_par%peat_area .or. peat_par%peat_carb) then
      f_peat     = peat_par%f_peat_min*f_veg
      f_peat_pot = peat_par%f_peat_min*f_veg
    else
      f_peat     = 0._wp
      f_peat_pot = 0._wp
    endif
    w_table_min  = 0._wp
    w_table_peat = 0._wp
    dCpeat_dt    = 0._wp

    do n=1,npft
      z0m(n) = 0.1_wp * veg_h(n)
    enddo
    z0m(i_bare) = surf_par%z0m_bare

    return

  end subroutine lndvc_init_cell_veg

end module lndvc_init_cell_mod
