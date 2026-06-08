!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s i c _ d y n _ m o d
!
!  Purpose : update sea-ice dynamics using elastic-viscous-plastic rheology 
!            with a C-grid discretization
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit
!
! This file is part of CLIMBER-X.
!
! This file was a part of SIS2.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                                                                              !
! C-grid SEA ICE DYNAMICS using ELASTIC-VISCOUS-PLASTIC RHEOLOGY adapted from  !
! Hunke and Dukowicz (JPO 1997, H&D hereafter) with some derivation from SIS1  !
! and with guidance from the C-grid implementations of sea-ice in MITgcm as    !
! documented in MITgcm user notes by Martin Losch and in LIM3 by S. Bouillon   !
! et al. (Ocean Modelling, 2009 & 2013). This code initially written by        !
! Robert Hallberg in 2013.                                                     !
!                                                                              !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
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
module sic_dyn_mod

  use precision, only : wp
  use constants, only : g
  use sic_params, only : dyn_par, rho_sic, rho_snow, l_diag_dyn
  use sic_grid, only : maxi, maxj, dx, dy, dxv, area, area_full, mask_u, mask_v, mask_q, fcorv


implicit none 

private

public :: sic_dyn

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> sic_dyn takes a single dynamics timestep with EVP subcycles
subroutine sic_dyn(dt, f_sic, h_sic_mean, h_snow_mean, uo, vo, tauxa, tauya, ssh, &
    ui, vi, str_d, str_t, str_s, &
    tauxo, tauyo, PFu, PFv, Cor_u, Cor_v, fxic, fxic_d, fxic_t, fxic_s, fyic, fyic_d, fyic_t, fyic_s)

  real(wp),                 intent(in   ) :: dt ! The amount of time over which the ice dynamics are to be advanced [s].
  real(wp), dimension(:,:), intent(in   ) :: f_sic ! Sea ice concentration [nondim]
  real(wp), dimension(:,:), intent(in   ) :: h_sic_mean    ! grid-cell mean sea ice thickness [m]
  real(wp), dimension(:,:), intent(in   ) :: h_snow_mean   ! grid-cell mean show thickness on sea ice [m] 
  real(wp), dimension(:,:), intent(in   ) :: uo    ! Zonal ocean velocity on u-grid [m s-1]
  real(wp), dimension(:,:), intent(in   ) :: vo    ! Meridional ocean velocity on v-grid [m s-1]
  real(wp), dimension(:,:), intent(in   ) :: tauxa ! Zonal air stress on ice on u-grid [Pa]
  real(wp), dimension(:,:), intent(in   ) :: tauya ! Meridional air stress on ice on v-grid [Pa]
  real(wp), dimension(:,:), intent(in   ) :: ssh ! The height of the sea surface [m].

  real(wp), dimension(:,:), intent(inout) :: ui    ! Zonal ice velocity on u-grid [m s-1]
  real(wp), dimension(:,:), intent(inout) :: vi    ! Meridional ice velocity on v-grid [m s-1]
  real(wp), dimension(:,:), intent(inout) :: str_d ! The divergence stress tensor component [Pa m].
  real(wp), dimension(:,:), intent(inout) :: str_t ! The tension stress tensor component [Pa m].
  real(wp), dimension(:,:), intent(inout) :: str_s ! The shearing stress tensor component [Pa m].

  real(wp), dimension(:,:), intent(  out) :: tauxo ! Zonal ice stress on ocean on u-grid [Pa]
  real(wp), dimension(:,:), intent(  out) :: tauyo ! Meridional ice stress on ocean on v-grid [Pa]

  real(wp), dimension(:,:) :: fxic   ! Zonal force due to internal stresses [Pa].
  real(wp), dimension(:,:) :: fxic_d ! Zonal force due to divergence internal stress [Pa].
  real(wp), dimension(:,:) :: fxic_t ! Zonal force due to tension internal stress [Pa].
  real(wp), dimension(:,:) :: fxic_s ! Zonal force due to shearing internal stress [Pa].
  real(wp), dimension(:,:) :: Cor_u  ! Zonal Coriolis acceleration [m s-2].
  real(wp), dimension(:,:) :: PFu    ! Zonal hydrostatic pressure driven acceleration [m s-2].
  real(wp), dimension(:,:) :: fyic   ! Meridional force due to internal stresses [Pa].
  real(wp), dimension(:,:) :: fyic_d ! Meridional force due to divergence internal stress [Pa].
  real(wp), dimension(:,:) :: fyic_t ! Meridional force due to tension internal stress [Pa].
  real(wp), dimension(:,:) :: fyic_s ! Meridional force due to shearing internal stress [Pa].
  real(wp), dimension(:,:) :: Cor_v  ! Meridional Coriolis acceleration [m s-2].
  real(wp), dimension(:,:) :: PFv    ! Meridional hydrostatic pressure driven acceleration [m s-2].

  real(wp), dimension(maxi,maxj) :: mis   ! Mass per unit ocean area of sea ice, snow and melt pond water [kg m-2]
  real(wp), dimension(maxi,maxj) :: mice  ! Mass per unit ocean area of sea ice [kg m-2]

  integer, dimension(maxi,maxj) :: mask_dyn

  real(wp), dimension(maxi,maxj) :: sh_Dt    ! sh_Dt is the horizontal tension (du/dx - dv/dy) including all metric terms [s-1].
  real(wp), dimension(maxi,maxj) :: sh_Dd    ! sh_Dd is the flow divergence (du/dx + dv/dy) including all metric terms [s-1].
  real(wp), dimension(maxi,maxj) :: sh_Ds  ! sh_Ds is the horizontal shearing strain (du/dy + dv/dx) including all metric terms [s-1].

  real(wp), dimension(maxi,maxj) :: pres_mice ! The ice internal pressure per unit column mass [N m kg-1].
  real(wp), dimension(maxi,maxj) :: f_sic_proj  ! The projected ice concentration [nondim].
  real(wp), dimension(maxi,maxj) :: zeta     ! The ice bulk viscosity [Pa m s] (i.e., [N s m-1]).
  real(wp), dimension(maxi,maxj) :: del_sh   ! The magnitude of the shear rates [s-1].
  real(wp), dimension(maxi,maxj) :: del_sh_min_pr  ! When multiplied by pres_mice, this gives the minimum value of del_sh that is used in the calculation of zeta [s-1].
                ! This is set based on considerations of numerical stability, and varies with the grid spacing.

  real(wp), dimension(maxi,maxj) :: u_tmp ! A temporary copy of the old values of ui [m s-1].
  real(wp), dimension(maxi,maxj) :: mi_u  ! The total ice and snow mass interpolated to u points [kg m-2].
  real(wp), dimension(maxi,maxj) :: f2dt_u  ! The squared effective Coriolis parameter at u-points times a time step [s-1].
  real(wp), dimension(maxi,maxj) :: I1_f2dt2_u  ! 1 / ( 1 + f^2 dt_evp^2) at u-points [nondim].

  real(wp), dimension(maxi,0:maxj) :: mi_v  ! The total ice and snow mass interpolated to v points [kg m-2].
  real(wp), dimension(maxi,maxj) :: f2dt_v ! The squared effective Coriolis parameter at v-points times a time step [s-1].
  real(wp), dimension(maxi,maxj) :: I1_f2dt2_v  ! 1 / ( 1 + f^2 dt_evp^2) at v-points [nondim].

  real(wp), dimension(maxi,maxj) :: mi_ratio_A_q ! A ratio of the masses interpolated to the faces around a vorticity point 
                                                         ! that ranges between (4 mi_min/mi_max) and 1, divided by the sum of the ocean areas 
                                                         ! around a point [m-2].
  real(wp), dimension(maxi,0:maxj) :: q ! A potential-vorticity-like field for the ice, the Coriolis parameter 
                                              ! divided by a spatially averaged mass per unit area [s-1 m2 kg-1].

  real(wp), dimension(maxi,maxj) :: &
    azon, bzon, & !  _zon & _mer are the values of the Coriolis force which
    czon, dzon, & ! are applied to the neighboring values of vi & ui,
    amer, bmer, & ! respectively to get the barotropic inertial rotation,
    cmer, dmer    ! in units of [s-1].  azon and amer couple the same pair of
                  ! velocities, but with the influence going in opposite
                  ! directions.

  real(wp), dimension(maxj) :: ui_min_cfl
  real(wp), dimension(maxj) :: ui_max_cfl
  real(wp), dimension(maxj) :: vi_min_cfl
  real(wp), dimension(maxj) :: vi_max_cfl

  real(wp) :: Cor       ! A Coriolis accleration [m s-2].
  real(wp) :: fxic_now  ! Zonal ice internal stress convergence [kg m-1 s-2].
  real(wp) :: fyic_now  ! Meridional ice internal stress convergence [kg m-1 s-2].
  real(wp) :: drag_u, drag_v ! Drag rates with the ocean at u & v points [kg m-2 s-1].
  real(wp) :: dxharm    ! The harmonic mean of the x- and y- grid spacings [m].
  real(wp) :: muq2, mvq2  ! The product of the u- and v-face masses per unit cell area surrounding a vorticity point [kg2 m-4].
  real(wp) :: muq, mvq    ! The u- and v-face masses per unit cell area extrapolated to a vorticity point on the coast [kg m-2].
  real(wp) :: I_1pdt_T    ! 1.0 / (1.0 + dt_2Tdamp) [nondim].
  real(wp) :: I_1pE2dt_T  ! 1.0 / (1.0 + EC^2 * dt_2Tdamp) [nondim].

  real(wp) :: v2_at_u     ! The squared v-velocity interpolated to u points [m2 s-2].
  real(wp) :: u2_at_v     ! The squared u-velocity interpolated to v points [m2 s-2].
  real(wp) :: uio_init    ! Ice-ocean velocity differences [m s-1]
  real(wp) :: vio_init    ! Ice-ocean velocity differences [m s-1]
  real(wp) :: m_uio_explicit ! Ice-ocean x-velocity differences times the ice mass [kg m-1 s-1]
  real(wp) :: m_vio_explicit ! Ice-ocean y-velocity differences times the ice mass [kg m-1 s-1]
  real(wp) :: uio_pred    ! Ice-ocean x-velocity differences [m s-1]
  real(wp) :: vio_pred    ! Ice-ocean y-velocity differences [m s-1]
  real(wp) :: I_cdRhoDt   ! The inverse of the product of the drag coefficient, ocean density and timestep [m3 kg-1 s-1].
  real(wp) :: cdRho       ! The ice density times the drag coefficient and rescaling factors [kg m-3]
  real(wp) :: b_vel0      ! The initial difference between the velocity magnitude and the absolute value of the u- or v- component, plus
                      ! the ice thickness divided by the time step and the drag coefficient [m s-1].
  real(wp) :: uio_C   ! A u-velocity difference between the ocean and ice [m s-1].
  real(wp) :: vio_C   ! A v-velocity difference between the ocean and ice [m s-1].

  real(wp) :: Tdamp   ! The damping timescale of the stress tensor components toward their equilibrium solution due to the elastic terms [s].
  real(wp) :: dt_evp  ! The short timestep associated with the EVP dynamics [s].
  real(wp) :: dt_2Tdamp ! The ratio of the timestep to the elastic damping timescale [nondim].
  real(wp) :: dt_cumulative ! The elapsed time within this call to EVP dynamics [s].
  real(wp) :: I_sub_steps  ! The number inverse of the number of EVP time steps per slow time step.
  real(wp) :: EC2     ! EC^2, where EC is the yield curve axis ratio.
  real(wp) :: I_EC2   ! 1/EC^2, where EC is the yield curve axis ratio.
  real(wp) :: I_EC    ! 1/EC, where EC is the yield curve axis ratio.
  real(wp) :: I_2EC   ! 1/(2*EC), where EC is the yield curve axis ratio.
  real(wp), parameter :: H_subroundoff = 1e-10_wp ! A negligible ice thickness [m].
  real(wp) :: m_neglect  ! A tiny mass per unit area [kg m-2].
  real(wp) :: m_neglect2 ! A tiny mass per unit area squared [kg2 m-4].
  real(wp) :: m_neglect4 ! A tiny mass per unit area to the 4th power [kg4 m-8].
  real(wp) :: sum_area   ! The sum of ocean areas around a vorticity point [m2].

  integer :: i, j, n, im1, ip1, jm1, jp1


  dt_evp = dt/dyn_par%evp_sub_steps

  cdRho = dyn_par%cdw * dyn_par%Rho_ocean

  I_cdRhoDt = 1._wp / (dyn_par%cdw * dyn_par%Rho_ocean * dt_evp)

  EC2 = dyn_par%EC**2
  I_EC = 0._wp ; if (dyn_par%EC > 0._wp) I_EC = 1._wp / dyn_par%EC
  I_2EC = 0._wp ; if (dyn_par%EC > 0._wp) I_2EC = 0.5_wp / dyn_par%EC
  I_EC2 = 0._wp ; if (EC2 > 0._wp) I_EC2 = 1._wp / EC2

  Tdamp = dyn_par%Tdamp
  if (dyn_par%Tdamp == 0._wp) then
    ! Hunke (2001) chooses a specified multiple (0.36) of dt for Tdamp, and shows that
    ! stability requires Tdamp > 2*dt_evp.  Here 0.2 is used instead for greater stability.
    Tdamp = max(0.2_wp*dt, 3._wp*dt_evp)
  elseif (dyn_par%Tdamp < 0._wp) then
    Tdamp = max(-dyn_par%Tdamp*dt, 3._wp*dt_evp)
  endif
  dt_2Tdamp = dt_evp / (2._wp * Tdamp)

  I_1pdt_T = 1._wp / (1._wp + dt_2Tdamp)
  I_1pE2dt_T = 1._wp / (1._wp + EC2*dt_2Tdamp)

  m_neglect = H_subroundoff*dyn_par%Rho_ice
  m_neglect2 = m_neglect**2 
  m_neglect4 = m_neglect**4

  ! velocity limits for CFL
  ! the teoretical limit for 2D advection is 1./sqrt(1./dx(j)**2+1./dy**2)/dt
  ! so use dx for both u and v components, as it is always the limiting factor
  do j=1,maxj
    ui_min_cfl(j) = -dyn_par%cfl_fac * dx(j) / dt
    ui_max_cfl(j) =  dyn_par%cfl_fac * dx(j) / dt
    vi_min_cfl(j) = -dyn_par%cfl_fac * dxv(J) / dt
    vi_max_cfl(j) =  dyn_par%cfl_fac * dxv(J) / dt
  enddo

  !$omp parallel 

  !$omp do collapse(2) private(i,j)
  do j=1,maxj
    do i=1,maxi
      ! Mass per unit ocean area of sea ice [kg m-2]
      mice(i,j) = h_sic_mean(i,j)*rho_sic
      ! Mass per unit ocean area of sea ice, snow and melt pond water [kg m-2]
      mis(i,j)  = h_sic_mean(i,j)*rho_sic + h_snow_mean(i,j)*rho_snow

      ! Zero these arrays to accumulate sums.
      tauxo(I,j) = 0._wp ;  tauyo(i,J) = 0._wp
      fxic(I,j) = 0._wp ;   fyic(i,J) = 0._wp
      Cor_u(I,j) = 0._wp ;  Cor_v(i,J) = 0._wp
      fxic_d(I,j) = 0._wp ; fyic_d(i,J) = 0._wp
      fxic_t(I,j) = 0._wp ; fyic_t(i,J) = 0._wp
      fxic_s(I,j) = 0._wp ; fyic_s(i,J) = 0._wp

      ! initialize because values outside mask_dyn area accessed
      sh_Ds(I,J) = 0._wp
      zeta(i,j) = 0._wp
      amer(I,j) = 0._wp
      cmer(i,J) = 0._wp
      u_tmp(I,j) = 0._wp
    enddo
  enddo
  !$omp end do

  !$omp do collapse(2) private(i,j,dxharm)
  do j=1,maxj
    do i=1,maxi
      f_sic_proj(i,j) = f_sic(i,j)

      ! Precompute pres_mice and the minimum value of del_sh for stability.
      pres_mice(i,j) = dyn_par%p0_rho*exp(-dyn_par%c0*max(1._wp-f_sic(i,j),0._wp))

      dxharm = 2._wp*dx(j)*dy / (dx(j) + dy)
      !   Setting a minimum value of del_sh is sufficient to guarantee numerical
      ! stability of the overall time-stepping for the velocities and stresses.
      ! Setting a minimum value of the shear magnitudes is equivalent to setting
      ! a maximum value of the effective lateral viscosities.
      ! I think that this is stable when dyn_par%del_sh_min_scale >= 1.  -RWH
      if (dxharm > 0._wp) then
        del_sh_min_pr(i,j) = (2._wp*dyn_par%del_sh_min_scale * dt_evp**2) / (Tdamp * dxharm**2)
      else
        del_sh_min_pr(i,j) = 0._wp
      endif
    enddo 
  enddo
  !$omp end do

  ! Ensure that the input stresses are not larger than could be justified by
  ! the ice pressure now, as the ice might have melted or been advected away
  ! during the thermodynamic and transport phases.
  call limit_stresses(pres_mice, mice(:,:), str_d, str_t, str_s)

  !$omp do collapse(2) private(i,j,im1,ip1,jm1,jp1,sum_area,muq2,mvq2,muq,mvq)
  do j=1,maxj
    do i=1,maxi

      im1 = i-1 
      if (im1.eq.0) im1=maxi
      ip1 = i+1 
      if (ip1.eq.maxi+1) ip1=1
      jm1 = j-1 
      if (jm1.eq.0) jm1=1
      jp1 = j+1 
      if (jp1.eq.maxj+1) jp1=maxj

      ! Zero out ice velocities with no mass.
      if (mask_u(I,j) * (mis(i,j)+mis(ip1,j)) == 0._wp) ui(I,j) = 0._wp
      if (mask_v(i,J) * (mis(i,j)+mis(i,jp1)) == 0._wp) vi(i,J) = 0._wp

      ! define mask where to compute sea ice velocities
      if (( mis(im1,jm1)+mis(i,jm1)+mis(ip1,jm1) + mis(im1,j)+mis(i,j)+mis(ip1,j) + mis(im1,jp1)+mis(i,jp1)+mis(ip1,jp1)) == 0._wp) then
        mask_dyn(i,j) = 0
      else
        mask_dyn(i,j) = 1
      endif

      !
      if (j.lt.maxj) then
        sum_area = (area(i,j) + area(ip1,j+1)) + (area(i,j+1) + area(ip1,j))
      else
        sum_area = (area(i,j) + area(ip1,j))
      endif

      if (sum_area <= 0._wp) then
        ! This is a land point.
        mi_ratio_A_q(I,J) = 0._wp
      elseif (mask_q(I,J)>0.5_wp) then
        ! This is an interior ocean point.
        !   Determine an appropriately averaged mass on q-points. The following
        ! expression for mi_q is mi when the masses are all equal, and goes to 4
        ! times the smallest mass averaged onto the 4 adjacent velocity points.  It
        ! comes from taking the harmonic means of the harmonic means of the
        ! arithmetic mean masses at the velocity points.  mi_ratio goes from 4 times
        ! the ratio of the smallest mass at velocity points over the largest mass
        ! at velocity points up to 1.
        muq2 = 0.25_wp * (mis(i,j) + mis(ip1,j)) * (mis(i,j+1) + mis(ip1,j+1))
        mvq2 = 0.25_wp * (mis(i,j) + mis(i,j+1)) * (mis(ip1,j) + mis(ip1,j+1))
        mi_ratio_A_q(I,J) = 32._wp * muq2 * mvq2 / ((m_neglect4 + (muq2 + mvq2) * &
          ((mis(i,j) + mis(ip1,j+1)) + (mis(i,j+1) + mis(ip1,j)))**2) * sum_area)
      elseif (j.lt.maxj .and. (mask_u(I,j) + mask_u(I,jp1)) + (mask_v(i,J) + mask_v(ip1,J)) > 1.5_wp) then
        !   This is a corner point, and there are 1 unmasked u-point and 1 v-point.
        ! The ratio below goes from 4 times the ratio of the smaller of the two
        ! masses at velocity points over the larger up to 1.
        muq = 0.5_wp * (mask_u(I,j) * (mis(i,j) + mis(ip1,j)) + &
          mask_u(I,j+1) * (mis(i,j+1) + mis(ip1,j+1)) )
        mvq = 0.5_wp * (mask_v(i,J) * (mis(i,j) + mis(i,j+1)) + &
          mask_v(ip1,J) * (mis(ip1,j) + mis(ip1,j+1)) )
        mi_ratio_A_q(I,J) = 4._wp * muq * mvq / ((m_neglect2 + (muq + mvq)**2) * sum_area)
      else
        ! This is a straight coastline or all neighboring velocity points are
        ! masked out.  In any case, with just 1 point, the ratio is always 1.
        mi_ratio_A_q(I,J) = 1._wp / sum_area
      endif

    enddo
  enddo
  !$omp end do

  !$omp do collapse(2) private(i,j,ip1,sum_area)
  do j=0,maxj
    do i=1,maxi

      ip1 = i+1 
      if (ip1.eq.maxi+1) ip1=1

      ! ice mass on velocity points
      if (j.gt.0) then
        mi_u(I,j) = 0.5_wp*(mis(ip1,j) + mis(i,j))
      endif
      if (j.eq.0) then
        mi_v(i,J) = 0.5_wp*mis(i,j+1)
      else if (j.eq.maxj) then
        mi_v(i,J) = 0.5_wp*mis(i,j)
      else
        mi_v(i,J) = 0.5_wp*(mis(i,j+1) + mis(i,j))
      endif

      if (j.eq.0) then
        sum_area = (area_full(ip1,j+1)) + (area_full(i,j+1))
        q(I,J) = fcorv(J) * sum_area / &
          (((area_full(i,j+1) * mis(i,j+1)) + &
          (area_full(ip1,j+1) * mis(ip1,j+1))) + sum_area * m_neglect)
      else if (j.eq.maxj) then
        sum_area = (area_full(i,j)) + (area_full(ip1,j))
        q(I,J) = fcorv(J) * sum_area / &
          (((area_full(i,j) * mis(i,j)) + &
          (area_full(ip1,j) * mis(ip1,j))) + sum_area * m_neglect)
      else
        sum_area = (area_full(i,j) + area_full(ip1,j+1)) + (area_full(i,j+1) + area_full(ip1,j))
        q(I,J) = fcorv(J) * sum_area / &
          (((area_full(i,j) * mis(i,j) + area_full(ip1,j+1) * mis(ip1,j+1)) + &
          (area_full(ip1,j) * mis(ip1,j) + area_full(i,j+1) * mis(i,j+1))) + sum_area * m_neglect)
      endif

    enddo 
  enddo
  !$omp end do

  !$omp do collapse(2) private(i,j,im1,ip1)
  do j=1,maxj
    do i=1,maxi

      im1 = i-1 
      if (im1.eq.0) im1=maxi
      ip1 = i+1 
      if (ip1.eq.maxi+1) ip1=1

      ! Calculate terms related to the Coriolis force on the zonal velocity.
      azon(I,j) = 0.25_wp * mi_v(ip1,J) * q(I,J)
      bzon(I,j) = 0.25_wp * mi_v(i,J) * q(I,J)
      czon(I,j) = 0.25_wp * mi_v(i,J-1) * q(I,J-1)
      dzon(I,j) = 0.25_wp * mi_v(ip1,J-1) * q(I,J-1)

      f2dt_u(I,j) = dt_evp * 4._wp * ((azon(I,j)**2 + czon(I,j)**2) + &
        (bzon(I,j)**2 + dzon(I,j)**2))
      I1_f2dt2_u(I,j) = 1._wp / ( 1._wp + dt_evp * f2dt_u(I,j) )

      ! Calculate the zonal acceleration due to the sea level slope.
      PFu(I,j) = -g*(ssh(ip1,j)-ssh(i,j)) / dx(j) 

      if (j.lt.maxj) then
        ! Calculate terms related to the Coriolis force on the meridional velocity.
        amer(Im1,j) = 0.25_wp * mi_u(Im1,j) * q(Im1,J)
        bmer(I,j) = 0.25_wp * mi_u(I,j) * q(I,J)
        cmer(I,j+1) = 0.25_wp * mi_u(I,j+1) * q(I,J)
        dmer(Im1,j+1) = 0.25_wp * mi_u(Im1,j+1) * q(Im1,J)

        f2dt_v(i,J) = dt_evp * 4._wp * ((amer(Im1,j)**2 + cmer(I,j+1)**2) + &
          (bmer(I,j)**2 + dmer(Im1,j+1)**2))
        I1_f2dt2_v(i,J) = 1._wp / ( 1._wp + dt_evp * f2dt_v(i,J) )

        ! Calculate the meridional acceleration due to the sea level slope.
        PFv(i,J) = -g*(ssh(i,j+1)-ssh(i,j)) / dy 
      endif
    enddo 
  enddo
  !$omp end do

  !$omp end parallel

  dt_cumulative = 0._wp

  ! Do the iterative time steps.
  ! A single parallel region spans the whole EVP subcycle
  !$omp parallel private(n)
  do n=1,dyn_par%evp_sub_steps

    ! single thread accumulates the elapsed time (kept exact for reproducibility)
    !$omp single
    dt_cumulative = dt_cumulative + dt_evp
    !$omp end single

    !    Calculate the strain tensor for viscosities and forcing elastic eqn.
    !  The following are the forms of the horizontal tension and hori-
    !  shearing strain advocated by Smagorinsky (1993) and discussed in
    !  Griffies and Hallberg (MWR, 2000).  Similar forms are used in the sea
    !  ice model of Bouillon et al. (Ocean Modelling, 2009).

    !$omp do collapse(2) private(i,j,ip1,im1)
    do j=1,maxj
      do i=1,maxi

        if (mask_dyn(i,j).eq.1) then

        ip1 = i+1 
        if (ip1.eq.maxi+1) ip1=1
        im1 = i-1 
        if (im1.eq.0) im1=maxi

        ! This uses a no-slip boundary condition.
        if (j.lt.maxj) then
          sh_Ds(I,J) = (2._wp-mask_q(I,J)) * &
            (dxv(J)/dy*(ui(I,j+1)/dx(j+1) - ui(I,j)/dx(j)) + &
            dy/dxv(J)*(vi(ip1,J)/dy      - vi(i,J)/dy))
        else
          sh_Ds(I,J) = (2._wp-mask_q(I,J)) * &
            (dy/dxv(J)*(vi(ip1,J)/dy      - vi(i,J)/dy))
        endif
        if (j.gt.1) then
          sh_Dt(i,j) = (dy/dx(j)   *(1._wp/dy *ui(I,j)      - 1._wp/dy *ui(Im1,j)) - &
            dx(j)/dy   *(1._wp/dxv(J) * vi(i,J) - 1._wp/dxv(J-1) *vi(i,J-1)))
          sh_Dd(i,j) = (1._wp/(dx(j)*dy) *(dy*ui(I,j)       - dy*ui(Im1,j)) + &
            1._wp/(dx(j)*dy) *(dxv(J) * vi(i,J) - dxv(J-1) *vi(i,J-1)))
        else
          sh_Dt(i,j) = (dy/dx(j)   *(1._wp/dy *ui(I,j)      - 1./dy *ui(Im1,j)) - &
            dx(j)/dy   *(1._wp/dxv(J) * vi(i,J)))
          sh_Dd(i,j) = (1._wp/(dx(j)*dy) *(dy*ui(I,j)       - dy*ui(Im1,j)) + &
            1._wp/(dx(j)*dy) *(dxv(J) * vi(i,J)))
        endif

        if (dyn_par%project_f_sic) then
          ! Estimate future ice concentrations from the approximate expression
          !   d f_sic / dt_evp = - f_sic * sh_Dt
          ! The choice to base this on the final velocity, the initial concentration
          ! and the elapsed time is because it is that final velocity that will drive
          ! ice convergence.
          f_sic_proj(i,j) = f_sic(i,j) * exp(-dt_cumulative*sh_Dd(i,j))
          ! Recompute pres_mice.
          pres_mice(i,j) = dyn_par%p0_rho*exp(-dyn_par%c0*max(1._wp-f_sic_proj(i,j),0._wp))
        endif

      endif

      enddo 
    enddo
    !$omp end do

    ! calculate viscosities - how often should we do this ?
    !$omp do collapse(2) private(i,j,im1,jm1)
    do j=1,maxj
      do i=1,maxi

        if (mask_dyn(i,j).eq.1) then

        im1 = i-1 
        if (im1.eq.0) im1=maxi
        jm1 = j-1 
        if (jm1.eq.0) jm1=1

        ! Averaging the squared shearing strain is larger than squaring
        ! the averaged strain.  I don't know what is better. -RWH

        del_sh(i,j) = sqrt(sh_Dd(i,j)**2 + I_EC2 * (sh_Dt(i,j)**2 + &
          (0.25_wp * ((sh_Ds(Im1,Jm1) + sh_Ds(I,J)) + &
          (sh_Ds(Im1,J) + sh_Ds(I,Jm1))))**2 ) ) ! H&D eqn 9

        if (max(del_sh(i,j), del_sh_min_pr(i,j)*pres_mice(i,j)) /=0._wp) then
          zeta(i,j) = 0.5_wp*pres_mice(i,j)*mice(i,j) / &
            max(del_sh(i,j), del_sh_min_pr(i,j)*pres_mice(i,j))
        else
          zeta(i,j) = 0._wp
        endif

        endif

      enddo 
    enddo
    !$omp end do

    ! Step the stress component equations semi-implicitly.
    !$omp do collapse(2) private(i,j,ip1)
    do j=1,maxj
      do i=1,maxi

        if (mask_dyn(i,j).eq.1) then

        ip1 = i+1 
        if (ip1.eq.maxi+1) ip1=1

        ! This expression uses that Pres=2*del_sh*zeta with an elliptic yield curve.
        str_d(i,j) = I_1pdt_T * ( str_d(i,j) + dt_2Tdamp * &
          ( zeta(i,j) * sh_Dd(i,j) - 0.5_wp*pres_mice(i,j)*mice(i,j) ) )
        str_t(i,j) = I_1pdt_T * ( str_t(i,j) + (I_EC2 * dt_2Tdamp) * &
          ( zeta(i,j) * sh_Dt(i,j) ) )
        ! zeta is already set to 0 over land.
        if (j.lt.maxj) then
          str_s(I,J) = I_1pdt_T * ( str_s(I,J) + (I_EC2 * dt_2Tdamp) * &
            ( ((area_full(i,j)*zeta(i,j) + area_full(ip1,j+1)*zeta(ip1,j+1)) + &
            (area_full(ip1,j)*zeta(ip1,j) + area_full(i,j+1)*zeta(i,j+1))) * &
            mi_ratio_A_q(I,J) * sh_Ds(I,J) ) )
        else
          str_s(I,J) = I_1pdt_T * ( str_s(I,J) + (I_EC2 * dt_2Tdamp) * &
            ( ((area_full(i,j)*zeta(i,j) ) + &
            (area_full(ip1,j)*zeta(ip1,j) )) * &
            mi_ratio_A_q(I,J) * sh_Ds(I,J) ) )
        endif

        endif

      enddo 
    enddo
    !$omp end do

    !$omp do collapse(2) private(i,j,ip1,Cor,fxic_now,v2_at_u,uio_init,drag_u,m_uio_explicit,b_vel0,uio_pred,uio_C)
    do j=1,maxj
      do i=1,maxi

        if (mask_dyn(i,j)==1) then

        ip1 = i+1 
        if (ip1.eq.maxi+1) ip1=1

        ! Save the current values of u for later use in updating v.
        u_tmp(I,j) = ui(I,j)

        if (j.gt.1) then
          Cor = ((azon(I,j) * vi(ip1,J) + czon(I,j) * vi(i,J-1)) + &
            (bzon(I,j) * vi(i,J) + dzon(I,j) * vi(ip1,J-1))) ! - Cor_ref_u(I,j)
        else
          Cor = ((azon(I,j) * vi(ip1,J)) + &
            (bzon(I,j) * vi(i,J))) ! - Cor_ref_u(I,j)
        endif

        !  Evaluate 1/m x.Div(m strain).  This expressions include all metric terms
        !  for an orthogonal grid.  The str_d term integrates out to no curl, while
        !  str_s & str_t terms impose no divergence and do not act on solid body rotation.
        if (j.gt.1) then
          fxic_now = 1._wp/dx(j) * (str_d(ip1,j) - str_d(i,j)) + &
            (1._wp/dy *(dy**2*str_t(ip1,j) - dy**2*str_t(i,j)) + &
            1._wp/dx(j) *(dxv(J)**2  *str_s(I,J) - dxv(J-1)**2*str_s(I,J-1)) ) / (dx(j)*dy) 
          v2_at_u =  dyn_par%drag_bg_vel2 + 0.25_wp * &
            (((vi(i,J)-vo(i,J))**2 + (vi(ip1,J-1)-vo(ip1,J-1))**2) + &
            ((vi(ip1,J)-vo(ip1,J))**2 + (vi(i,J-1)-vo(i,J-1))**2))
        else
          fxic_now = 1._wp/dx(j) * (str_d(ip1,j) - str_d(i,j)) + &
            (1._wp/dy *(dy**2*str_t(ip1,j) - dy**2*str_t(i,j)) + &
            1._wp/dx(j) *(dxv(j)**2  *str_s(I,J)) ) / (dx(j)*dy) 
          v2_at_u =  dyn_par%drag_bg_vel2 + 0.25_wp * &
            (((vi(i,J)-vo(i,J))**2) + &
            ((vi(ip1,J)-vo(ip1,J))**2))
        endif

        uio_init = (ui(I,j)-uo(I,j))

        if (dyn_par%project_drag_vel) then
          ! Project the new u-velocity using a quasi-analytic implicit treatment for
          ! drag, but explicit treatments for everything else, to estimate the drag
          ! coefficient, then take the larger of the two estimates of
          ! the ice-ocean drag.
          drag_u = 0._wp
          if (mask_u(I,j) > 0._wp) then
            m_uio_explicit = uio_init*mi_u(I,j) + dt_evp * &
              ((Cor + PFu(I,j))*mi_u(I,j) + (fxic_now + tauxa(I,j)))
            b_vel0 = mi_u(I,j) * I_cdRhoDt + &
              ( sqrt(uio_init**2 + v2_at_u) - abs(uio_init) )
            if (b_vel0**2 > 1.e8_wp*I_cdRhoDt*abs(m_uio_explicit)) then
              uio_pred = m_uio_explicit * I_cdRhoDt / b_vel0
            else
              uio_pred = 0.5_wp * (sqrt(b_vel0**2 + 4._wp*I_cdRhoDt*abs(m_uio_explicit)) - b_vel0)
            endif
            drag_u = cdRho * sqrt(max(uio_init**2, uio_pred**2) + v2_at_u )
          endif
        else
          drag_u = cdRho * sqrt(uio_init**2 + v2_at_u )
        endif

        if (l_diag_dyn) then
          ! Determine the Coriolis acceleration and sum for averages...
          Cor_u(I,j) = Cor_u(I,j) + (Cor - f2dt_u(I,j) * ui(I,j)) * I1_f2dt2_u(I,j)
          ! sum accelerations to take averages.
          fxic(I,j) = fxic(I,j) + fxic_now
          fxic_d(I,j) = fxic_d(I,j) + mask_u(I,j) * 1._wp/dx(j) * (str_d(ip1,j) - str_d(i,j))
          fxic_t(I,j) = fxic_t(I,j) + mask_u(I,j) * 1._wp/dy *(dy**2* str_t(ip1,j) - dy**2* str_t(i,j) ) / (dx(j)*dy) 
        endif

        if (j.gt.1) fxic_s(I,j) = fxic_s(I,j) + mask_u(I,j) * 1._wp/dx(j) *(dxv(j)**2*str_s(I,J) - dxv(j)**2*str_s(I,J-1)) / (dx(j)*dy) 
        !   This is a quasi-implicit timestep of Coriolis, followed by an explicit
        ! update of the other terms and an implicit bottom drag calculation.
        uio_C =  mask_u(I,j) * ( mi_u(I,j) * &
          ((ui(I,j) + dt_evp * Cor) * I1_f2dt2_u(I,j) - uo(I,j)) + &
          dt_evp * (mi_u(I,j) * PFu(I,j) + (fxic_now + tauxa(I,j))) ) / &
          (mi_u(I,j) + m_neglect + dt_evp * drag_u)

        ui(I,j) = (uio_C + uo(I,j)) * mask_u(I,j)

        if (ui(I,j) < ui_min_cfl(j)) then
          ui(I,j) = ui_min_cfl(j)
        elseif (ui(I,j) > ui_max_cfl(j)) then
          ui(I,j) = ui_max_cfl(j)
        endif

        ! Note that tauxo is the stress felt by the ocean.
        tauxo(I,j) = tauxo(I,j) + drag_u*uio_C

      endif

      enddo
    enddo
    !$omp end do

    !$omp do collapse(2) private(i,j,im1,Cor,fyic_now,u2_at_v,vio_init,drag_v,m_vio_explicit,b_vel0,vio_pred,vio_C)
    do j=1,maxj-1
      do i=1,maxi

        if (mask_dyn(i,j)==1) then

        im1 = i-1 
        if (im1.eq.0) im1=maxi

        Cor = -1._wp*((amer(Im1,j) * u_tmp(Im1,j) + cmer(I,j+1) * u_tmp(I,j+1)) + &
          (bmer(I,j) * u_tmp(I,j) + dmer(Im1,j+1) * u_tmp(Im1,j+1)))
        !  Evaluate 1/m y.Div(m strain).  This expressions include all metric terms
        !  for an orthogonal grid.  The str_d term integrates out to no curl, while
        !  str_s & str_t terms impose no divergence and do not act on solid body rotation.
        fyic_now = 1._wp/dy * (str_d(i,j+1)-str_d(i,j)) + &
          (-1._wp/dxv(J) *(dx(j+1)**2*str_t(i,j+1) - &
          dx(j)**2  *str_t(i,j)) + &
          1._wp/dy *(dy**2  *str_s(I,J) - &
          dy**2  *str_s(Im1,J)) ) / (dxv(J)*dy) 
        u2_at_v = dyn_par%drag_bg_vel2 + 0.25_wp * &
          (((u_tmp(I,j)-uo(I,j))**2 + (u_tmp(Im1,j+1)-uo(Im1,j+1))**2) + &
          ((u_tmp(I,j+1)-uo(I,j+1))**2 + (u_tmp(Im1,j)-uo(Im1,j))**2))

        vio_init = (vi(i,J)-vo(i,J))

        if (dyn_par%project_drag_vel) then
          ! Project the new v-velocity using a quasi-analytic implicit treatment for
          ! drag, but explicit treatments for everything else, to estimate the drag
          ! coefficient, then take the larger of the two estimates of
          ! the ice-ocean drag.

          drag_v = 0._wp
          if (mask_v(i,J) > 0._wp) then
            m_vio_explicit = vio_init*mi_v(i,J) + dt_evp * &
              ((Cor + PFv(i,J))*mi_v(i,J) + (fyic_now + tauya(i,J)))
            b_vel0 = mi_v(i,J) * I_cdRhoDt + (sqrt(vio_init**2 + u2_at_v) - abs(vio_init))
            if (b_vel0**2 > 1.e8_wp*I_cdRhoDt*abs(m_vio_explicit)) then
              vio_pred = m_vio_explicit * I_cdRhoDt / b_vel0
            else
              vio_pred = 0.5_wp * (sqrt(b_vel0**2 + 4._wp*I_cdRhoDt*abs(m_vio_explicit)) - b_vel0)
            endif
            drag_v = cdRho * sqrt(max(vio_init**2, vio_pred**2) + u2_at_v )
          endif
        else
          drag_v = cdRho * sqrt(vio_init**2 + u2_at_v )
        endif

        if (l_diag_dyn) then
          ! Determine the Coriolis acceleration and sum for averages...
          Cor_v(I,J) = Cor_v(I,J) + (Cor - f2dt_v(i,J) * vi(i,J)) * I1_f2dt2_v(i,J)
          ! sum accelerations to take averages.
          fyic(i,J) = fyic(i,J) + fyic_now
          fyic_d(i,J) = fyic_d(i,J) + mask_v(i,J) * 1._wp/dy * (str_d(i,j+1)-str_d(i,j))
          fyic_t(i,J) = fyic_t(i,J) + mask_v(i,J) * (1._wp/dxv(J)*(dx(j+1)**2*(-str_t(i,j+1)) - dx(j)**2  *(-str_t(i,j))) ) / (dxv(J)*dy) 
          fyic_s(i,J) = fyic_s(i,J) + mask_v(i,J) * (1._wp/dy *(dy**2*str_s(I,J) - dy**2*str_s(Im1,J)) ) / (dxv(J)*dy)
        endif

        !   This is a quasi-implicit timestep of Coriolis, followed by an explicit
        ! update of the other terms and an implicit bottom drag calculation.
        vio_C =  mask_v(i,J) * ( mi_v(i,J) * &
          ((vi(i,J) + dt_evp * Cor) * I1_f2dt2_v(i,J) - vo(i,J)) + &
          dt_evp * (mi_v(i,J) * PFv(i,J) + (fyic_now + tauya(i,J))) ) / &
          (mi_v(i,J) + m_neglect + dt_evp * drag_v)

        vi(i,J) = (vio_C + vo(i,J)) * mask_v(i,J)

        if (vi(i,J) < vi_min_cfl(J)) then
          vi(i,J) = vi_min_cfl(J)
        elseif (vi(i,J) > vi_max_cfl(J)) then
          vi(i,J) = vi_max_cfl(J)
        endif

        ! Note that tauyo is the stress felt by the ocean.
        tauyo(i,J) = tauyo(i,J) + drag_v*vio_C

      endif

      enddo
    enddo
    !$omp end do

  enddo ! l=1,evp_sub_steps
  !$omp end parallel

  ! make averages
  I_sub_steps = 1._wp/dyn_par%evp_sub_steps
  do j=1,maxj
    do i=1,maxi
      tauxo(I,j) = tauxo(I,j) * (mask_u(I,j) * I_sub_steps)
      tauyo(i,J) = tauyo(i,J) * (mask_v(i,J) * I_sub_steps)
      if (l_diag_dyn) then
        fxic(I,j) = fxic(I,j) * (mask_u(I,j) * I_sub_steps)
        Cor_u(I,j) = Cor_u(I,j) * (mask_u(I,j) * I_sub_steps)
        fxic_d(I,j) = fxic_d(I,j) * (mask_u(I,j) * I_sub_steps)
        fxic_t(I,j) = fxic_t(I,j) * (mask_u(I,j) * I_sub_steps)
        fxic_s(I,j) = fxic_s(I,j) * (mask_u(I,j) * I_sub_steps)
        fyic(i,J) = fyic(i,J) * (mask_v(i,J) * I_sub_steps)
        Cor_v(i,J) = Cor_v(i,J) * (mask_v(i,J) * I_sub_steps)
        fyic_d(i,J) = fyic_d(i,J) * (mask_v(i,J) * I_sub_steps)
        fyic_t(i,J) = fyic_t(i,J) * (mask_v(i,J) * I_sub_steps)
        fyic_s(i,J) = fyic_s(i,J) * (mask_v(i,J) * I_sub_steps)
      endif
    enddo 
  enddo

end subroutine sic_dyn

!> limit_stresses ensures that the input stresses are not larger than could be justified by the ice
!! pressure now, as the ice might have melted or been advected away during the thermodynamic and
!! transport phases, or the ice flow convergence or divergence may have altered the ice concentration.
subroutine limit_stresses(pres_mice, mice, str_d, str_t, str_s, limit)
  real(wp), dimension(:,:), intent(in)    :: pres_mice ! The ice internal pressure per unit column mass [N m kg-1].
  real(wp), dimension(:,:), intent(in)    :: mice  ! The mass per unit total area (ice covered and ice free) of the ice [kg m-2].
  real(wp), dimension(:,:), intent(inout) :: str_d ! The divergence stress tensor component [Pa m].
  real(wp), dimension(:,:), intent(inout) :: str_t ! The tension stress tensor component [Pa m].
  real(wp), dimension(:,:), intent(inout) :: str_s ! The shearing stress tensor component [Pa m].
  real(wp), optional,       intent(in)    :: limit ! A factor by which the strength limits are changed.

!   This subroutine ensures that the input stresses are not larger than could
! be justified by the ice pressure now, as the ice might have melted or been
! advected away during the thermodynamic and transport phases, or the
! ice flow convergence or divergence may have altered the ice concentration.

  ! Local variables
  real(wp) :: pressure  ! The integrated internal ice pressure at a point [Pa m].
  real(wp) :: pres_avg  ! The average of the internal ice pressures around a point [Pa m].
  real(wp) :: sum_area  ! The sum of ocean areas around a vorticity point [m2].
  real(wp) :: I_2EC     ! 1/(2*EC), where EC is the yield curve axis ratio.
  real(wp) :: lim       ! A local copy of the factor by which the limits are changed.
  real(wp) :: lim_2     ! The limit divided by 2.

  integer :: i, j, ip1

  lim = 1._wp ; if (present(limit)) lim = limit
  I_2EC = 0._wp ; if (dyn_par%EC > 0._wp) I_2EC = (0.5_wp*lim) / dyn_par%EC
  lim_2 = 0.5_wp * lim

  ! The rescaling here is done separately for each component.
  !$omp do collapse(2) private(i,j,ip1,pressure,sum_area,pres_avg)
  do j=1,maxj
    do i=1,maxi

      ip1 = i+1 
      if (ip1.eq.maxi+1) ip1=1

      pressure = pres_mice(i,j)*mice(i,j)
      if (str_d(i,j) < -pressure) str_d(i,j) = -pressure
      if (dyn_par%EC*str_t(i,j) > lim_2*pressure) str_t(i,j) = I_2EC*pressure
      if (dyn_par%EC*str_t(i,j) < -lim_2*pressure) str_t(i,j) = -I_2EC*pressure
      if (j.lt.maxj) then
        sum_area = (area(i,j) + area(ip1,j+1)) + (area(i,j+1) + area(ip1,j))
      else
        sum_area = (area(i,j)) + (area(ip1,j))
      endif
      pres_avg = 0._wp
      if (j.lt.maxj) then
        if (sum_area > 0._wp) &
          pres_avg = ((area(i,j)   * (pres_mice(i,j)*mice(i,j)) + &
            area(ip1,j+1)*(pres_mice(ip1,j+1)*mice(ip1,j+1))) + &
            (area(ip1,j) * (pres_mice(ip1,j)*mice(ip1,j)) + &
            area(i,j+1) * (pres_mice(i,j+1)*mice(i,j+1)))) / sum_area
      else
        if (sum_area > 0._wp) &
          pres_avg = ((area(i,j)   * (pres_mice(i,j)*mice(i,j))) + &
            (area(ip1,j) * (pres_mice(ip1,j)*mice(ip1,j)))) / sum_area
      endif

      if (dyn_par%EC*str_s(I,J) > lim_2*pres_avg) str_s(I,J) = I_2EC*pres_avg
      if (dyn_par%EC*str_s(I,J) < -lim_2*pres_avg) str_s(I,J) = -I_2EC*pres_avg
    enddo 
  enddo
  !$omp end do

end subroutine limit_stresses

end module sic_dyn_mod
