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
  use constants, only : T0, pi, fqsat
  use atm_params, only : p0, c_cld_1, c_cld_2, c_cld_3, c_cld_4, c_cld_5, c_cld_6, c_cld_60_ocn, c_cld_60_lnd, c_cld_7, i_cld_low, q_cld_low
  use atm_params, only : c_cld_8, c_cld_9, c_cld_10, c_cld_5_sc
  use atm_params, only : c_hcld_low
  use atm_params, only : cld_max, nsmooth_cld
  use atm_params, only : c_hcld_1, c_hcld_2, c_hcld_3, c_hcld_4
  use atm_params, only : c_clot_1, c_clot_2, c_clot_3, c_clot_4, c_clot_5
  use atm_params, only : l_so4_ie, r_so4, N_so4_nat, N_so4_nat_lnd
  use atm_grid, only : im, jm, i_ice, i_ocn, i_sic, i_lnd, i_lake
  use smooth_atm_mod, only : smooth2
  use vesta_mod, only : t_prof
  !$ use omp_lib

  implicit none

  private
  public :: clouds 

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  c l o u d s 
  !   Purpose    :  compute cloud fraction, top height and optical depth
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine clouds(frst, weff, wcld, zsa, t2a, ram, qam, r2a, q2a, wcon, htrop, so4, &
      tam, gams, gamb, gamt, &
      fweff, cld_rh, cld_low, cld, hcld, hcld_rh, hcld_low, clot)

    implicit none

    real(wp), intent(in   ) :: frst(:,:,:)
    real(wp), intent(in   ) :: zsa(:,:)
    real(wp), intent(in   ) :: t2a(:,:)
    real(wp), intent(in   ) :: ram(:,:)
    real(wp), intent(in   ) :: qam(:,:)
    real(wp), intent(in   ) :: r2a(:,:)
    real(wp), intent(in   ) :: q2a(:,:)
    real(wp), intent(in   ) :: wcon(:,:)
    real(wp), intent(in   ) :: htrop(:,:)
    real(wp), intent(in   ) :: weff(:,:)
    real(wp), intent(in   ) :: wcld(:,:)
    real(wp), intent(in   ) :: so4(:,:)
    real(wp), intent(in   ) :: tam(:,:)
    real(wp), intent(in   ) :: gams(:,:)
    real(wp), intent(in   ) :: gamb(:,:)
    real(wp), intent(in   ) :: gamt(:,:)

    real(wp), intent(out  ) :: fweff(:,:)
    real(wp), intent(out  ) :: cld_rh(:,:)
    real(wp), intent(out  ) :: cld_low(:,:)
    real(wp), intent(out  ) :: hcld_rh(:,:)
    real(wp), intent(out  ) :: hcld_low(:,:)

    real(wp), intent(inout) :: cld(:,:)
    real(wp), intent(inout) :: hcld(:,:)
    real(wp), intent(inout) :: clot(:,:)

    integer :: i, j
    real(wp) :: dr, fr, f_freezedry, f1, c_cld_60, f_warm, f_sc, c5_eff
    real(wp) :: cldn, hcldl, clotl, tcldm, ftemp, fdyn
    real(wp) :: f_rh, a_rh, a_low, t_rh, t_low, t_eff
    real(wp) :: L_so4_ant, N_so4_ant_0, N_so4, f_mod

    real(wp), parameter :: cld_min=0.1_wp          ! minimum cloud fraction
    real(wp), parameter :: alpha_c = 2.5e-9_wp     ! m^3
    real(wp), parameter :: rho_so4 = 1.769e3_wp    ! kg/m3, sulfate aerosol density
    real(wp), parameter :: H_so4 = 1500._wp        ! m, sulfate aerosol height scale


    !$omp parallel do collapse(2) private(i, j, c_cld_60, dr, fr, f_freezedry, f1, f_warm, f_sc, c5_eff, hcldl, clotl, tcldm, ftemp, fdyn) &
    !$omp private(L_so4_ant, N_so4_ant_0, N_so4, f_mod)
    do j=1,jm
      do i=1,im

        ! effective vertical velocity factor for cloud parameterization, scaled between -1 and 1
        fweff(i,j) = tanh(c_cld_3*weff(i,j))

        !--------------------------------------------
        ! cloud fraction
        !--------------------------------------------

        ! relative humidity gradient, a measure of surface inversion
        c_cld_60 = (frst(i,j,i_ocn)+frst(i,j,i_sic)+frst(i,j,i_lake))*c_cld_60_ocn + frst(i,j,i_lnd)*c_cld_60_lnd
        dr = (r2a(i,j)-ram(i,j)-c_cld_60)/(c_cld_6+1.e-20_wp)
        dr = min(dr, 1._wp)
        dr = max(dr,-1._wp)

        ! low clouds related to surface inversion
        ! 'freezedry' reduction of cloud cover, Vavrus & Walliser (2008)
        f_freezedry = min(1._wp, 0.5_wp+0.5_wp*qam(i,j)/(c_cld_7+1.e-20_wp))

        ! relative weight of low clouds, linear in the normalised dr
        fr = f_freezedry*0.5_wp*(1._wp+dr)

        if (i_cld_low.eq.1) then
          
          cld_low(i,j) = (1._wp-frst(i,j,i_ice))*c_cld_5*fr

        else if (i_cld_low.eq.2) then
          
          cld_low(i,j) = (1._wp-frst(i,j,i_ice))*c_cld_5*fr*ram(i,j)

        else if (i_cld_low.eq.3) then

          cld_low(i,j) = (1._wp-frst(i,j,i_ice))*c_cld_5*fr * min(1._wp,q_cld_low/fqsat(t2a(i,j),p0))

        else if (i_cld_low.eq.4) then

          ! Multiplicative suppression of low clouds in warm, moist columns.
          ! fr measures the surface inversion through the humidity contrast r2a-ram, but every
          ! such contrast measure increases with warming, so on its own it gives low clouds a
          ! negative amount feedback. f_warm adds the missing absolute-warmth dependence, the
          ! analogue of the negative SST term in the observational two-predictor low cloud
          ! models: a deeper, more decoupled boundary layer entrains more dry air and sustains less stratocumulus. 
          f_warm = exp(-max(0._wp,qam(i,j)-c_cld_8)/(c_cld_9+1.e-20_wp))

          ! f_sc selects the warm subtropical Sc columns and
          ! gives them their own amplitude c_cld_5_sc. The ramp c_cld_10 is deliberately much
          ! narrower than c_cld_9, so it has saturated before it reaches the decks and the
          ! amplitude carries no sensitivity of its own; c_cld_5_sc then sets the present-day
          ! deck cloud and c_cld_9 the response to warming, with little cross-talk.
          f_sc = 0.5_wp*(1._wp+tanh((qam(i,j)-c_cld_8)/(c_cld_10+1.e-20_wp)))
          c5_eff = c_cld_5 + (c_cld_5_sc-c_cld_5)*f_sc

          cld_low(i,j) = (1._wp-frst(i,j,i_ice))*c5_eff*fr * f_warm
          cld_low(i,j) = min(1._wp,cld_low(i,j))
          cld_low(i,j) = max(0._wp,cld_low(i,j))

        endif

        ! clouds related to large scale atmospheric relative humidity
        cld_rh(i,j) = (c_cld_1+c_cld_2*fweff(i,j))*ram(i,j)**c_cld_4 

        !--------------------------------------------
        ! cloud height
        !--------------------------------------------

        ! top of the clouds associated with the large scale relative humidity
        hcldl = c_hcld_1 + c_hcld_2*htrop(i,j) * (1._wp+c_hcld_3*(wcld(i,j)-c_hcld_4))
        hcldl = min(hcldl,htrop(i,j)-1.e3_wp)
        hcldl = max(hcldl,zsa(i,j)+2.5e3_wp)
        hcld_rh(i,j) = hcldl

        ! top of the low clouds, a fixed height above the surface so that it warms with the
        ! surface and its longwave effect stays close to zero under warming
        hcld_low(i,j) = min(zsa(i,j)+c_hcld_low, htrop(i,j)-1.e3_wp)

        ! the merging of the two cloud tops and the time smoothing of hcld need cld_low after
        ! the horizontal smoothing below, so they are done in the second loop

        !--------------------------------------------
        ! cloud optical thickness
        !--------------------------------------------

        tcldm = t2a(i,j)-T0-c_clot_1
        ftemp = 1._wp+tanh(-tcldm/c_clot_2)
        ftemp = min(1._wp,ftemp)
        fdyn = max(0.1_wp, 1._wp+c_clot_5*fweff(i,j))
        clotl = c_clot_3*ftemp*(cld(i,j)*wcon(i,j))**c_clot_4 * fdyn

        ! cloud droplet number concentration at the cloud base (m^-3). The natural background is
        ! higher over land than over ocean because continental air carries far more CCN (observed
        ! ~250-600 cm-3 against ~60-150 cm-3 over the ocean); at a given water path this is what
        ! makes continental cloud optically thicker. 
        N_so4 = N_so4_nat + (N_so4_nat_lnd-N_so4_nat)*frst(i,j,i_lnd)

        ! indirect effect of anthropogenic sulfate aerosols
        if (l_so4_ie) then
          ! column integrated anthropogenic aerosol number burden (m^-2)
          L_so4_ant = so4(i,j)/(rho_so4*4._wp/3._wp*pi*r_so4**3)
          ! anthropogenic sulphate aerosol number concentration at the surface (m^-3)
          N_so4_ant_0 = L_so4_ant/H_so4
          ! add it at the cloud base
          N_so4 = N_so4 + N_so4_ant_0*exp(-1._wp)
        endif

        ! Twomey: at a fixed water path the optical thickness scales as N^(1/3), with the
        ! Boucher & Lohmann saturation at high N. Normalised to the marine background, so
        ! N_so4_nat_lnd = N_so4_nat reproduces the previous behaviour exactly.
        f_mod = (1._wp-exp(-alpha_c*N_so4))/(1._wp-exp(-alpha_c*N_so4_nat))
        clotl = clotl*f_mod**0.33_wp

        clotl = min(10._wp, clotl)           

        ! smooth in time
        clot(i,j) = 0.1_wp*clotl + 0.9_wp*clot(i,j)

      enddo
    enddo
    !$omp end parallel do

    call smooth2(cld_low,nsmooth_cld)

    !$omp parallel do collapse(2) private(i, j, cldn, hcldl, f_rh, a_rh, a_low, t_rh, t_low, t_eff)
    do j=1,jm
      do i=1,im
        ! total cloud fraction
        cldn = 1._wp-(1._wp-cld_rh(i,j))*(1._wp-cld_low(i,j))
        cldn = max(cldn,cld_min)
        cldn = min(cldn,cld_max)
        cld(i,j) = 0.1_wp*cldn + 0.9_wp*cld(i,j)

        !--------------------------------------------
        ! effective radiating cloud top
        !--------------------------------------------

        ! The two cloud types share one optical thickness and therefore one emissivity, so the
        ! (1-eps) clear-sky terms of the two cloudy skies recombine with their areas exactly and
        ! only the emission has to be matched:
        !   cld*T_eff^4 = a_rh*T_rh^4 + a_low*T_low^4
        ! The areas are the random overlap that cld = 1-(1-cld_rh)*(1-cld_low) already implies,
        ! with the low cloud radiatively visible only where there is no high cloud above it.
        f_rh  = min(1._wp,max(0._wp,cld_rh(i,j)))
        a_rh  = f_rh
        a_low = min(1._wp,max(0._wp,cld_low(i,j)))*(1._wp-f_rh)
        if ((a_rh+a_low).gt.1.e-10_wp) then
          t_rh  = t_prof(zsa(i,j), hcld_rh(i,j),  tam(i,j), gams(i,j), gamb(i,j), gamt(i,j), htrop(i,j), 1)
          t_low = t_prof(zsa(i,j), hcld_low(i,j), tam(i,j), gams(i,j), gamb(i,j), gamt(i,j), htrop(i,j), 1)
          t_eff = ((a_rh*t_rh**4 + a_low*t_low**4)/(a_rh+a_low))**0.25_wp
          hcldl = z_of_tprof(t_eff, zsa(i,j), tam(i,j), gams(i,j), gamb(i,j), gamt(i,j), htrop(i,j))
        else
          hcldl = hcld_rh(i,j)
        endif
        hcldl = min(hcldl,htrop(i,j)-1.e3_wp)
        hcldl = max(hcldl,zsa(i,j))

        ! smooth in time
        hcld(i,j) = 0.1_wp*hcldl + 0.9_wp*hcld(i,j)

      enddo
    enddo
    !$omp end parallel do

    return

  end subroutine clouds


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  z _ o f _ t p r o f
  !   Purpose    :  height at which the temperature profile equals a given temperature
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Inverse of t_prof, needed to turn the effective radiating temperature of the merged
  ! cloud back into the cloud top height the radiation expects. t_prof is piecewise in the
  ! surface layer, so it is inverted by bisection rather than analytically; the profile
  ! decreases monotonically with height, and 30 bisections resolve the whole troposphere to
  ! well below a metre. Temperatures outside the profile range return the corresponding end.
  pure function z_of_tprof(t, zs, tam, gams, gamb, gamt, htrop)

    implicit none

    real(wp), intent(in) :: t
    real(wp), intent(in) :: zs
    real(wp), intent(in) :: tam
    real(wp), intent(in) :: gams
    real(wp), intent(in) :: gamb
    real(wp), intent(in) :: gamt
    real(wp), intent(in) :: htrop

    real(wp) :: z_of_tprof

    integer :: n
    real(wp) :: za, zb, zm, ta, tb

    integer, parameter :: n_iter = 30


    za = zs
    zb = max(htrop,zs)
    ta = t_prof(zs, za, tam, gams, gamb, gamt, htrop, 1)
    tb = t_prof(zs, zb, tam, gams, gamb, gamt, htrop, 1)

    if (t.ge.ta .or. (ta-tb).lt.1.e-10_wp) then
      z_of_tprof = za
    else if (t.le.tb) then
      z_of_tprof = zb
    else
      do n=1,n_iter
        zm = 0.5_wp*(za+zb)
        if (t_prof(zs, zm, tam, gams, gamb, gamt, htrop, 1).gt.t) then
          za = zm
        else
          zb = zm
        endif
      enddo
      z_of_tprof = 0.5_wp*(za+zb)
    endif

    return

  end function z_of_tprof

end module clouds_mod
