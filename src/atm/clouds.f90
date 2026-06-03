!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module :  c l o u d s _ m o d
!
!  Purpose : computation of cloud properties
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Andrey Ganopolski and Matteo Willeit
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
module clouds_mod

  use atm_params, only : wp
  use constants, only : T0, pi
  use atm_params, only : c_cld_1, c_cld_2, c_cld_3, c_cld_4, c_cld_5, c_cld_55, c_cld_6, c_cld_7, c_cld_8, l_cld_low_ice, cld_max, nsmooth_cld
  use atm_params, only : c_hcld_1, c_hcld_2, c_hcld_3, c_hcld_4
  use atm_params, only : c_clot_1, c_clot_2, c_clot_3, c_clot_4
  use atm_params, only : l_so4_ie, r_so4, N_so4_nat
  use atm_grid, only : im, jm, i_ice
  use smooth_atm_mod, only : smooth2
  !$ use omp_lib

  implicit none

  private
  public :: clouds 

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  c l o u d s 
  !   Purpose    :  compute cloud fraction, top height and optical depth
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine clouds(frst, weff, wcld, zsa, t2a, ram, qam, rskina, wcon, htrop, so4, sam2, &
      fweff, cld_rh, cld_low, cld, hcld, clot)

    implicit none

    real(wp), intent(in   ) :: frst(:,:,:)
    real(wp), intent(in   ) :: zsa(:,:)
    real(wp), intent(in   ) :: t2a(:,:)
    real(wp), intent(in   ) :: ram(:,:)
    real(wp), intent(in   ) :: qam(:,:,:)
    real(wp), intent(in   ) :: rskina(:,:)
    real(wp), intent(in   ) :: wcon(:,:)
    real(wp), intent(in   ) :: htrop(:,:)
    real(wp), intent(in   ) :: weff(:,:)
    real(wp), intent(in   ) :: wcld(:,:)
    real(wp), intent(in   ) :: so4(:,:)
    real(wp), intent(in   ) :: sam2(:,:)

    real(wp), intent(out  ) :: fweff(:,:)
    real(wp), intent(out  ) :: cld_rh(:,:)
    real(wp), intent(out  ) :: cld_low(:,:)

    real(wp), intent(inout) :: cld(:,:)
    real(wp), intent(inout) :: hcld(:,:)
    real(wp), intent(inout) :: clot(:,:)

    integer :: i, j
    real(wp) :: dr, fr, f_freezedry
    real(wp) :: cldn, hcldl, clotl, tcldm, ftemp
    real(wp) :: L_so4_ant, N_so4_ant_0, N_so4, f_mod

    real(wp), parameter :: cld_min=0.1_wp          ! minimum cloud fraction
    real(wp), parameter :: alpha_c = 2.5e-9_wp     ! m^3
    real(wp), parameter :: rho_so4 = 1.769e3_wp    ! kg/m3, sulfate aerosol density
    real(wp), parameter :: H_so4 = 1500._wp        ! m, sulfate aerosol height scale


    !$omp parallel do collapse(2) private(i, j, dr, fr, f_freezedry, hcldl, clotl, tcldm, ftemp) &
    !$omp private(L_so4_ant, N_so4_ant_0, N_so4, f_mod)
    do j=1,jm
      do i=1,im

        ! effective vertical velocity factor for cloud parameterization, scaled between -1 and 1
        fweff(i,j) = tanh(c_cld_3*weff(i,j))

        !--------------------------------------------
        ! cloud fraction
        !--------------------------------------------

        ! near-surface relative humidity gradient, a measure of surface inversion
        dr = rskina(i,j)-ram(i,j)
        dr = min(dr,c_cld_6)
        dr = max(dr,-c_cld_6)

        ! low clouds related to surface inversion
        ! 'freezedry' reduction of cloud cover, Vavrus & Walliser (2008)
        f_freezedry = 0.1_wp+0.9_wp*qam(i,j,1)/(c_cld_7+1.e-20_wp)
        f_freezedry = min(1._wp,f_freezedry)
        ! relative weight of low clouds 
        fr = f_freezedry*(dr+c_cld_6)/(2._wp*c_cld_6+1.e-20_wp) 
        if (l_cld_low_ice) then
          cld_low(i,j) = c_cld_5*fr*ram(i,j)**c_cld_55
        else
          cld_low(i,j) = (1._wp-frst(i,j,i_ice))*c_cld_5*fr*ram(i,j)**c_cld_55
        endif

        ! clouds related to large scale atmospheric relative humidity
        cld_rh(i,j) = (c_cld_1+c_cld_2*fweff(i,j))*ram(i,j)**c_cld_4 + c_cld_8*max(0._wp,sam2(i,j)-20._wp)

        !--------------------------------------------
        ! cloud height
        !--------------------------------------------

        hcldl = c_hcld_1 + c_hcld_2*htrop(i,j) * (1._wp+c_hcld_3*(wcld(i,j)-c_hcld_4))
        hcldl = min(hcldl,htrop(i,j)-1.e3_wp)
        hcldl = max(hcldl,zsa(i,j)+2.5e3_wp)

        ! smooth in time
        hcld(i,j) = 0.1_wp*hcldl + 0.9_wp*hcld(i,j)

        !--------------------------------------------
        ! cloud optical thickness
        !--------------------------------------------

        tcldm = t2a(i,j)-T0-c_clot_1
        ftemp = 1._wp+tanh(-tcldm/c_clot_2)
        ftemp = min(1._wp,ftemp)
        clotl = c_clot_3*ftemp*(cld(i,j)*wcon(i,j))**c_clot_4 
        clotl = min(10._wp, clotl)           

        ! indirect effect of sulfate aerosols on cloud optical thickness
        if (l_so4_ie) then
          ! column integrated anthropogenic aerosol number burden (m^-2)
          L_so4_ant = so4(i,j)/(rho_so4*4._wp/3._wp*pi*r_so4**3)
          ! anthropogenic sulphate aerosol number concentration at the surface (m^-3)
          N_so4_ant_0 = L_so4_ant/H_so4
          ! anthropogenic and natural sulphate aerosol number concentration at the cloud base (m^-3)
          N_so4 = N_so4_ant_0*exp(-1._wp) + N_so4_nat

          f_mod = (1._wp-exp(-alpha_c*N_so4))/(1._wp-exp(-alpha_c*N_so4_nat))

          clotl = clotl*f_mod**0.33_wp
        endif

        ! smooth in time
        clot(i,j) = 0.1_wp*clotl + 0.9_wp*clot(i,j)

      enddo
    enddo
    !$omp end parallel do

    call smooth2(cld_low,nsmooth_cld)

    do j=1,jm
      do i=1,im
        ! total cloud fraction
        cldn = 1._wp-(1._wp-cld_rh(i,j))*(1._wp-cld_low(i,j))
        cldn = max(cldn,cld_min)       
        cldn = min(cldn,cld_max)
        cld(i,j) = 0.1_wp*cldn + 0.9_wp*cld(i,j)
      enddo
    enddo

    return 

  end subroutine clouds

end module clouds_mod
