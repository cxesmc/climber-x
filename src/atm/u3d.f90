!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : u 3 d _ m o d
!
!  Purpose : 3-D wind field
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
module u3d_mod

  use atm_params, only : wp
  use constants, only : g, T0, r_earth, pi
  use atm_params, only : amas, ra, i_pbl, dpc, dp_com, ptop_com, ptopdyn
  use atm_params, only : i_mass_com_topo, dps_com_topo, i_mass_com_vert
  use atm_params, only : c_uter_pol, c_uter_eq
  use atm_grid, only : im, imc, jm, jmc, km, kmc, k500, k700, dxt, dxu, dy, zl, sqr, aim
  use atm_grid, only : fcort, cost, sint
  use atm_grid, only : pl, dplx, dply, dplxo, dplyo, plx, ply, pblt, pblu

  implicit none
  
  private
  public :: u3d

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s o l v e _ p o i s s o n
  !   Purpose    :  solve laplacian(psi) = rhs on the sphere
  !
  !   The discrete operator is the one whose convergence matches the mass flux
  !   stencil used for the columns,
  !       conv(i,j) = f(i,j)-f(i+1,j) + f(i,j+1)-f(i,j)  ,
  !   with the fluxes formed as  fx = dy/dxt * (psi(i-1)-psi(i))  and
  !   fy = dxu/dy * (psi(j)-psi(j-1)). The correction therefore cancels the
  !   column convergence exactly rather than approximately.
  !
  !   Solved by Fourier transform in longitude (the coefficients depend on
  !   latitude only) and a tridiagonal solve in latitude for each zonal
  !   wavenumber. dxu vanishes at both poles, which imposes the no-flux boundary
  !   condition. For wavenumber zero the operator is singular - psi is defined up
  !   to a constant - so psi is pinned at the first row; this requires the global
  !   sum of rhs to vanish, which is the global mass constraint.
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine solve_poisson(rhs, psi)

    use, intrinsic :: iso_c_binding
    use atm_params, only : dp

    implicit none

    include 'fftw3.f03'

    real(wp), intent(in ) :: rhs(:,:)
    real(wp), intent(out) :: psi(:,:)

    integer :: j, n, nn
    real(wp) :: lamb
    real(wp) :: wx(jm), wy(jmc)
    real(wp), dimension(jm) :: aa, bb, cc
    complex(dp), dimension(jm) :: rr, gam, xx
    complex(dp) :: bet
    complex(dp), dimension(im/2+1,jm) :: rhs_fft
    real(dp), dimension(im), save :: row
    complex(dp), dimension(im/2+1), save :: row_fft
    type(C_PTR), save :: plan_r2c, plan_c2r
    logical, save :: plans_ready = .false.


    ! metric weights of the flux-form gradient
    do j=1,jm
      wx(j) = dy/dxt(j)
    enddo
    do j=1,jmc
      wy(j) = dxu(j)/dy
    enddo

    if (.not.plans_ready) then
      plan_r2c = fftw_plan_dft_r2c_1d(im, row, row_fft, FFTW_ESTIMATE)
      plan_c2r = fftw_plan_dft_c2r_1d(im, row_fft, row, FFTW_ESTIMATE)
      plans_ready = .true.
    endif

    do j=1,jm
      row(:) = rhs(:,j)
      call fftw_execute_dft_r2c(plan_r2c, row, row_fft)
      rhs_fft(:,j) = row_fft(:)
    enddo

    do n=1,im/2+1
      nn = n-1
      ! eigenvalue of the zonal second difference for wavenumber nn
      lamb = -4._wp*sin(pi*real(nn,wp)/real(im,wp))**2
      do j=1,jm
        aa(j) = wy(j)
        cc(j) = wy(j+1)
        bb(j) = -(wy(j)+wy(j+1)) + lamb*wx(j)
        rr(j) = rhs_fft(n,j)
      enddo
      if (nn.eq.0) then
        bb(1) = 1._wp
        cc(1) = 0._wp
        rr(1) = (0._dp,0._dp)
      endif
      ! Thomas algorithm
      bet = bb(1)
      xx(1) = rr(1)/bet
      do j=2,jm
        gam(j) = cc(j-1)/bet
        bet = bb(j)-aa(j)*gam(j)
        xx(j) = (rr(j)-aa(j)*xx(j-1))/bet
      enddo
      do j=jm-1,1,-1
        xx(j) = xx(j)-gam(j+1)*xx(j+1)
      enddo
      rhs_fft(n,:) = xx(:)
    enddo

    do j=1,jm
      row_fft(:) = rhs_fft(:,j)
      call fftw_execute_dft_c2r(plan_c2r, row_fft, row)
      psi(:,j) = row(:) * aim   ! aim accounts for the FFTW normalisation
    enddo

    return

  end subroutine solve_poisson

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  u 3 d
  !   Purpose    :  computation of 3D wind field and advective mass fluxes
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine u3d(niter, pzsa, ptrop, ugb, vgb, uab, vab, t3, &
        ua, va, uter, vter, uterf, vterf, u3, v3, w3, uz500, &
        fax, faxo, fay, fayo, fac, fac_topo, psi, psi_topo, &
        fax_psi, fay_psi, fax_psi_topo, fay_psi_topo)

    implicit none

    integer,  intent(in   ) :: niter
    real(wp), intent(in   ) :: pzsa(:,:)
    real(wp), intent(in   ) :: ptrop(:)
    real(wp), intent(in   ) :: ugb(:,:)
    real(wp), intent(in   ) :: vgb(:,:)
    real(wp), intent(in   ) :: uab(:,:)
    real(wp), intent(in   ) :: vab(:,:)
    real(wp), intent(in   ) :: t3(:,:,:)

    real(wp), intent(out  ) :: ua(:,:,:)
    real(wp), intent(out  ) :: va(:,:,:)
    real(wp), intent(out  ) :: uter(:,:,:)
    real(wp), intent(out  ) :: vter(:,:,:)
    real(wp), intent(out  ) :: uterf(:,:,:)
    real(wp), intent(out  ) :: vterf(:,:,:)
    real(wp), intent(inout) :: u3(:,:,:)
    real(wp), intent(inout) :: v3(:,:,:)
    real(wp), intent(inout) :: w3(:,:,:)
    real(wp), intent(inout) :: uz500(:)
    real(wp), intent(out  ) :: fax(:,:,:)
    real(wp), intent(out  ) :: faxo(:,:,:)
    real(wp), intent(out  ) :: fay(:,:,:)
    real(wp), intent(out  ) :: fayo(:,:,:)
    real(wp), intent(out  ) :: fac(:,:)
    real(wp), intent(out  ) :: fac_topo(:,:)
    ! flux potential of the column mass correction, kg/s
    real(wp), intent(out  ) :: psi(:,:)
    ! flux potential of its topographic part, zero unless i_mass_com_topo==1
    real(wp), intent(out  ) :: psi_topo(:,:)
    ! the corrective mass flux itself, kept per level so that adifa can form the
    ! dry static energy convergence it carries with the same upstream values the
    ! advection uses; _topo is the grad(psi_topo) part, zero unless i_mass_com_topo==1
    real(wp), intent(out  ) :: fax_psi(:,:,:)
    real(wp), intent(out  ) :: fay_psi(:,:,:)
    real(wp), intent(out  ) :: fax_psi_topo(:,:,:)
    real(wp), intent(out  ) :: fay_psi_topo(:,:,:)

    integer :: i, j, k, n, ipl, imi, jmi, kpl
    real(wp) :: pzx, pzy, dp_c, pbl_t, pbl_u, dp, pc1, pc2, dp_l, pl1, pl2, fxpbl, fypbl, uabc, vabc, ctv
    real(wp) :: faz
    real(wp) :: fcx, fcy, fcxt, fcyt, fmean, fmean_topo, divbar, ptrx, ptry
    real(wp) :: dfx, dfy, ftrop, fsum
    real(wp) :: fcx_geo, fcx_ter, fcy_geo, fcy_ter, fmean_geo, fmean_ter
    ! per-component flux potentials and convergences, only used when i_mass_com_vert==1
    real(wp), dimension(im,jm) :: psi_geo, psi_ter, fac_geo, fac_topo_geo
    ! per-component mass fluxes; the fourth index selects the component
    ! (1 = barotropic geostrophic, 2 = thermal wind, 3 = ageostrophic PBL)
    real(wp), allocatable, dimension(:,:,:,:), save :: fax_c, fay_c
    real(wp), dimension(imc,jm) :: fxcol
    real(wp), dimension(im,jmc) :: fycol
    real(wp) :: u_g, v_g, ug_b, vg_b
    real(wp) :: u3k, u3kp1, v3k, v3kp1
    real(wp) :: dplxdy, dplydxu
    real(wp) :: c_damp_eq, c_damp_pol

    real(wp), dimension(kmc) :: uteri
    real(wp), dimension(kmc) :: vteri
    real(wp), dimension(km) :: rdpl   ! 1/(pl(k)-pl(k+1)), column-invariant


    ! per-component mass fluxes, allocated on first use; needed by i_mass_com_vert==1
    if (.not.allocated(fax_c)) then
      allocate(fax_c(imc,jm,km,3))
      allocate(fay_c(im,jmc,km,3))
      fax_c = 0._wp; fay_c = 0._wp
    endif

    ! precompute reciprocal layer-pressure thickness (depends on k only, not on i,j)
    do k=1,km
      rdpl(k) = 1._wp/(pl(k)-pl(k+1))
    enddo

    !$omp parallel do collapse(2) private(i,j,k,ipl,imi,jmi,pzx,pzy,dp,pbl_t,pbl_u,pc1,pc2,dp_l,pl1,pl2,fxpbl,fypbl,uabc,vabc,uteri,vteri,ctv,c_damp_eq,c_damp_pol)
    do j=1,jm
      do i=1,im

        ipl = modulo(i,im) + 1
        imi = modulo(i - 2, im) + 1
        jmi = max(1,j-1)

        !-------------------------------------------------------
        ! Vertical profile of ageostrophic wind

        ! x-component on u-points 

        ! velocity in the PBL
        if (i_pbl.eq.1) then
          pzx = 0.5_wp*(pzsa(i,j)+pzsa(imi,j))
        else if (i_pbl.eq.2) then
          pzx = 1._wp
        endif
        pbl_t = pzx+pblt(j)-1._wp
        do k=1,km
          dp = (pl(k)-pbl_t)*rdpl(k)
          dp = min(1._wp,dp)
          dp = max(0._wp,dp)
          ua(i,j,k) = uab(i,j) * dp
        enddo

        ! compensatory velocity in the upper troposphere
        if (i_pbl.eq.1) then
          fxpbl = uab(i,j)*(1._wp-pblt(j))        
        else if (i_pbl.eq.2) then
          fxpbl = uab(i,j)*max(0._wp,pzx-pblt(j))        
        endif
        uabc = -fxpbl/dpc
        pc2 = ptrop(j)
        pc1 = pc2+dpc
        do k=1,km
          pl1 = max(min(pl(k),pc1),pl(k+1))
          pl2 = min(max(pl(k+1),pc2),pl(k))
          dp_l = pl1-pl2
          ua(i,j,k) = ua(i,j,k)+uabc*dp_l*rdpl(k)
        enddo     

        ! y-component on v-points 

        ! velocity in the PBL
        if (i_pbl.eq.1) then
          pzy = 0.5_wp*(pzsa(i,j)+pzsa(i,jmi))
        else if (i_pbl.eq.2) then
          pzy = 1._wp
        endif
        pbl_u = pzy+pblu(j)-1._wp
        do k=1,km
          dp = (pl(k)-pbl_u)*rdpl(k)
          dp = min(1._wp,dp)
          dp = max(0._wp,dp)
          va(i,j,k) = vab(i,j) * dp
        enddo

        ! compensatory velocity in the upper troposphere
        if (i_pbl.eq.1) then
          fypbl = vab(i,j)*(1._wp-pblu(j))        
        else if (i_pbl.eq.2) then
          fypbl = vab(i,j)*max(0._wp,pzy-pblu(j))        
        endif
        vabc = -fypbl/dpc
        pc2 = 0.5_wp*(ptrop(j)+ptrop(jmi))  
        pc1 = pc2+dpc
        do k=1,km
          pl1 = max(min(pl(k),pc1),pl(k+1))
          pl2 = min(max(pl(k+1),pc2),pl(k))
          dp_l = pl1-pl2
          va(i,j,k) = va(i,j,k)+vabc*dp_l*rdpl(k)
        enddo  

        !-------------------------------------------------------
        ! Thermal wind on T-points

        if (j.eq.1 .or. j.eq.jm) then

          uter(i,j,:) = 0._wp
          vter(i,j,:) = 0._wp
          uterf(i,j,:) = 0._wp
          vterf(i,j,:) = 0._wp

        else

          ! thermal wind on levels
          uteri(1) = 0._wp
          vteri(1) = 0._wp
          uteri(km+1) = 0._wp
          vteri(km+1) = 0._wp
          do k=1,km-1
            ctv = (zl(k+1)-zl(k))*g/(T0*fcort(j)) 
            uteri(k+1) = uteri(k)-ctv*(t3(i,jmi,k)-t3(i,j+1,k))/(2._wp*dy)   
            vteri(k+1) = vteri(k)+ctv*(t3(ipl,j,k)-t3(imi,j,k))/(2._wp*dxt(j)) 
          enddo

          ! thermal wind in layers, dampened at equator and poles
          c_damp_pol = min(1._wp,c_uter_pol*cost(j)**2)
          c_damp_eq  = min(1._wp, c_uter_eq*sint(j)**2)
          do k=1,km
            uter(i,j,k) = 0.5_wp*(uteri(k)+uteri(k+1)) * c_damp_eq * c_damp_pol 
            vter(i,j,k) = 0.5_wp*(vteri(k)+vteri(k+1)) * c_damp_eq * c_damp_pol 
            ! thermal wind without polar damping needed for EKE production
            uterf(i,j,k) = 0.5_wp*(uteri(k)+uteri(k+1)) * c_damp_eq 
            vterf(i,j,k) = 0.5_wp*(vteri(k)+vteri(k+1)) * c_damp_eq
          enddo

        endif

      enddo
    enddo
    !$omp end parallel do

    !-------------------------------------------------------
    ! advective mass transport
    !-------------------------------------------------------

    !$omp parallel do collapse(2) private(i,j,k,imi,u_g,v_g,dplxdy,dplydxu,ftrop)
    do k=1,km
      do j=1,jm

        ! fraction of the layer that lies below the dynamical top; the geostrophic
        ! wind is limited to the troposphere, see below
        if (pl(k+1).ge.ptopdyn) then
          ftrop = 1._wp
        elseif (pl(k+1).lt.ptopdyn .and. pl(k).ge.ptopdyn) then
          ftrop = (pl(k)-ptopdyn)*rdpl(k)
        else
          ftrop = 0._wp
        endif

        ! x-components

        do i=1,im
          imi = modulo(i - 2, im) + 1
          ! geostrophic zonal wind on u-points, limited to troposphere
          if (pl(k+1).ge.ptopdyn) then
            u_g  = 0.5_wp*(ugb(imi,j)+ugb(i,j)) + 0.5_wp*(uter(imi,j,k)+uter(i,j,k))
          elseif (pl(k+1).lt.ptopdyn.and.pl(k).ge.ptopdyn) then
            u_g  = (0.5_wp*(ugb(imi,j)+ugb(i,j)) + 0.5_wp*(uter(imi,j,k)+uter(i,j,k))) *(pl(k)-ptopdyn)*rdpl(k)
          else
            u_g  = 0._wp
          endif 
          dplxdy = dplx(i,j,k)*dy
          ! mass flux
          fax(i,j,k) = (u_g+ua(i,j,k))*dplxdy ! m/s * kg/m2 * m = kg/s
          ! orographic component of mass flux 
          faxo(i,j,k) = u_g*(dplxo(i,j,k)*dy-dplxdy) 
          ! the same mass flux, split by wind component; the three add up to fax
          ! exactly, since ftrop reproduces the branch taken above
          fax_c(i,j,k,1) = 0.5_wp*(ugb(imi,j)+ugb(i,j))*ftrop * dplxdy
          fax_c(i,j,k,2) = 0.5_wp*(uter(imi,j,k)+uter(i,j,k))*ftrop * dplxdy
          fax_c(i,j,k,3) = ua(i,j,k) * dplxdy
        enddo
        ! periodic boundary conditions
        fax(imc,j,k) = fax(1,j,k)
        faxo(imc,j,k) = faxo(1,j,k)
        fax_c(imc,j,k,:) = fax_c(1,j,k,:)

        ! y-components

        do i=1,im
          if (j.eq.1) then
            ! N and S boundary conditions, no flux
            fay(i,1,k)  = 0._wp
            fayo(i,1,k) = 0._wp
            fay(i,jmc,k)  = 0._wp
            fayo(i,jmc,k) = 0._wp
            fay_c(i,1,k,:)   = 0._wp
            fay_c(i,jmc,k,:) = 0._wp
          else
            ! geostrophic meridional wind on u-points, limited to troposphere
            if (pl(k+1).ge.ptopdyn) then
              v_g  = 0.5_wp*(vgb(i,j-1)+vgb(i,j)) + 0.5_wp*(vter(i,j-1,k)+vter(i,j,k))
            elseif (pl(k+1).lt.ptopdyn.and.pl(k).ge.ptopdyn) then
              v_g  = (0.5_wp*(vgb(i,j-1)+vgb(i,j)) + 0.5_wp*(vter(i,j-1,k)+vter(i,j,k))) *(pl(k)-ptopdyn)*rdpl(k)
            else
              v_g  = 0._wp
            endif 
            dplydxu = dply(i,j,k)*dxu(j)
            ! mass flux
            fay(i,j,k) = (v_g+va(i,j,k))*dplydxu  ! m/s * kg/m2 * m = kg/s
            ! orographic component of mass flux for temperature
            fayo(i,j,k) = v_g*(dplyo(i,j,k)*dxu(j)-dplydxu) 
            ! the same mass flux, split by wind component
            fay_c(i,j,k,1) = 0.5_wp*(vgb(i,j-1)+vgb(i,j))*ftrop * dplydxu
            fay_c(i,j,k,2) = 0.5_wp*(vter(i,j-1,k)+vter(i,j,k))*ftrop * dplydxu
            fay_c(i,j,k,3) = va(i,j,k) * dplydxu
          endif
        enddo

      enddo
    enddo 
    !$omp end parallel do

    !-------------------------------------------------------
    ! per-column mass conservation
    !-------------------------------------------------------
    ! The column mass is diagnosed from topography and does not evolve, so the
    ! column-integrated mass flux divergence has to vanish: the same fluxes advect
    ! dry static energy and moisture, and advecting a uniform tracer must leave it
    ! unchanged. The ageostrophic wind closes by construction (see above), so what
    ! remains is the geostrophic and thermal-wind divergence, a large part of which
    ! is the topographic term V.grad(p_s) - a wind crossing an elevation gradient
    ! carries a different column mass in than out. That term is physically real and
    ! is what drives dp_s/dt in a model with a prognostic surface pressure; here it
    ! has nowhere to go and has to be removed.
    !
    ! The residual is removed with a corrective mass flux F_c = grad(psi)
    ! obtained from  laplacian(psi) = -(column convergence),  discretised so that its
    ! convergence matches the stencil used above. This is the smallest correction
    ! that closes the budget and, unlike adding the whole residual to one face, it
    ! has no preferred direction. Enforcing it per column also makes the zonal mean
    ! meridional mass flux vanish at every latitude, so no separate step is needed:
    ! summing the column constraint zonally telescopes the zonal fluxes away and,
    ! with no flux through the poles, leaves sum_i fay(i,j) = 0 for all j.
    !
    ! The weighting within the band is not a degree of freedom: dplx is the pressure
    ! thickness of the air-filled part of the layer, so dp_l*rdpl*dplx reduces to
    ! dp_l for any layer clear of the ground and every mode spreads its correction
    ! uniformly in pressure across its band. That reduction fails only where a layer
    ! is cut by the surface and still reaches into the band, which cannot happen for
    ! the topographic band below (it would need a layer thicker than dps_com_topo)
    ! and needs Antarctic or Tibetan terrain for the dp_com one.
    !
    ! i_mass_com_vert selects the vertical placement of what is left after the
    ! topographic term:
    !   0  the whole residual in a slab of depth dp_com hanging below the tropopause,
    !      so the band follows the tropopause down at high latitudes.
    !   1  split by wind component - one flux potential for the barotropic geostrophic
    !      part and one for the thermal wind - and each returned through the same
    !      layers, in the same proportion, that the component's own flux occupies at
    !      that face. The correction then carries the same dry static energy and
    !      moisture per unit mass as the flux that created the imbalance, so it is
    !      neutral in both by construction. The barotropic wind is height independent
    !      through the troposphere and so makes its divergence distributed by layer
    !      mass, while the thermal wind is a shear anchored at zero and makes its
    !      divergence aloft; dp_com returns both at the same place.
    !   2  as 0, but the slab hangs below the FIXED level ptop_com, so the correction
    !      is returned at the same pressure everywhere. Its depth is still dp_com, so
    !      unlike a fixed pressure PAIR it cannot collapse when the tropopause drops
    !      below its lower edge.
    !
    ! With i_mass_com_topo==1 the topographic term is separated from the rest and
    ! given its own potential, because the two belong at different heights. The
    ! mass that upslope flow piles into a column is in reality lifted over the
    ! barrier and leaves just above it, whereas the divergence of the balanced
    ! flow is compensated by the secondary circulation high up. The level matters:
    ! the return flow removes the moisture that the upslope flow just brought in
    ! if it is placed near the surface, and removes far too much dry static energy
    ! (which grows with height) if it is placed at the tropopause.

    ! column convergence before compensation, kg/s (also the diagnostic fac)
    !$omp parallel do collapse(2) private(i,j,k)
    do j=1,jm
      do i=1,im
        fac(i,j) = 0._wp
        do k=1,km
          fac(i,j) = fac(i,j) + fax(i,j,k)-fax(i+1,j,k) + fay(i,j+1,k)-fay(i,j,k)
        enddo
      enddo
    enddo
    !$omp end parallel do


    ! The global sum of the column convergence vanishes identically; removing the
    ! mean only guards against round-off, and is the solvability condition of the
    ! Poisson problem.
    fmean = 0._wp
    do j=1,jm
      do i=1,im
        fmean = fmean + fac(i,j)
      enddo
    enddo
    fmean = fmean*aim/real(jm,wp)

    ! Column convergence of the barotropic geostrophic component alone. The
    ! ageostrophic wind closes column-wise by construction, so its column flux
    ! vanishes at every face and fac = fac_geo + fac_ter exactly; the thermal wind
    ! part is therefore taken as the difference and the budget closes to round-off
    ! whatever the discretisation does.
    if (i_mass_com_vert.eq.1) then
      !$omp parallel do collapse(2) private(i,j,k)
      do j=1,jm
        do i=1,im
          fac_geo(i,j) = 0._wp
          do k=1,km
            fac_geo(i,j) = fac_geo(i,j) + fax_c(i,j,k,1)-fax_c(i+1,j,k,1) &
                                        + fay_c(i,j+1,k,1)-fay_c(i,j,k,1)
          enddo
        enddo
      enddo
      !$omp end parallel do
    else
      fac_geo(:,:) = 0._wp
    endif

    if (i_mass_com_topo.eq.1) then

      !-------------------------------------------------------
      ! split off the topographic part of the column convergence
      !-------------------------------------------------------
      ! Writing div(F) = div(p_s*V) = p_s*div(V) + V.grad(p_s) with V the mass
      ! weighted mean wind, the first term is the divergence of the flow itself and
      ! the second the purely topographic one. Discretely the mean wind is taken
      ! from the column integrated flux at the faces, and the topographic part is
      ! the residual, so that the two add up to fac to round-off.

      !$omp parallel do collapse(2) private(i,j,k)
      do j=1,jm
        do i=1,imc
          fxcol(i,j) = 0._wp
          do k=1,km
            fxcol(i,j) = fxcol(i,j) + fax(i,j,k)
          enddo
        enddo
      enddo
      !$omp end parallel do
      ! no flux through the poles
      fycol(:,1)   = 0._wp
      fycol(:,jmc) = 0._wp
      !$omp parallel do collapse(2) private(i,j,k)
      do j=2,jm
        do i=1,im
          fycol(i,j) = 0._wp
          do k=1,km
            fycol(i,j) = fycol(i,j) + fay(i,j,k)
          enddo
        enddo
      enddo
      !$omp end parallel do

      !$omp parallel do collapse(2) private(i,j,divbar)
      do j=1,jm
        do i=1,im
          ! convergence of the mass weighted mean wind, m2/s
          divbar = 0._wp
          if (plx(i,j).gt.0._wp)   divbar = divbar + fxcol(i,j)/plx(i,j)
          if (plx(i+1,j).gt.0._wp) divbar = divbar - fxcol(i+1,j)/plx(i+1,j)
          if (ply(i,j+1).gt.0._wp) divbar = divbar + fycol(i,j+1)/ply(i,j+1)
          if (ply(i,j).gt.0._wp)   divbar = divbar - fycol(i,j)/ply(i,j)
          ! kg/s, the column mass is the same combination used to build dplx/dply
          fac_topo(i,j) = fac(i,j) - (pzsa(i,j)-pl(kmc))*amas*divbar
        enddo
      enddo
      !$omp end parallel do

      ! Unlike the total, the topographic part has no reason to integrate to zero
      ! (it is the mountain-torque-like term), so its own mean is removed here and
      ! left to the dynamic part, which then carries the round-off residual fmean.
      fmean_topo = 0._wp
      do j=1,jm
        do i=1,im
          fmean_topo = fmean_topo + fac_topo(i,j)
        enddo
      enddo
      fmean_topo = fmean_topo*aim/real(jm,wp)

      if (i_mass_com_vert.eq.1) then
        ! the same construction applied to the barotropic component alone, so that
        ! fac_topo_geo + fac_topo_ter = fac_topo to round-off
        !$omp parallel do collapse(2) private(i,j,k)
        do j=1,jm
          do i=1,imc
            fxcol(i,j) = 0._wp
            do k=1,km
              fxcol(i,j) = fxcol(i,j) + fax_c(i,j,k,1)
            enddo
          enddo
        enddo
        !$omp end parallel do
        fycol(:,1)   = 0._wp
        fycol(:,jmc) = 0._wp
        !$omp parallel do collapse(2) private(i,j,k)
        do j=2,jm
          do i=1,im
            fycol(i,j) = 0._wp
            do k=1,km
              fycol(i,j) = fycol(i,j) + fay_c(i,j,k,1)
            enddo
          enddo
        enddo
        !$omp end parallel do
        !$omp parallel do collapse(2) private(i,j,divbar)
        do j=1,jm
          do i=1,im
            divbar = 0._wp
            if (plx(i,j).gt.0._wp)   divbar = divbar + fxcol(i,j)/plx(i,j)
            if (plx(i+1,j).gt.0._wp) divbar = divbar - fxcol(i+1,j)/plx(i+1,j)
            if (ply(i,j+1).gt.0._wp) divbar = divbar + fycol(i,j+1)/ply(i,j+1)
            if (ply(i,j).gt.0._wp)   divbar = divbar - fycol(i,j)/ply(i,j)
            fac_topo_geo(i,j) = fac_geo(i,j) - (pzsa(i,j)-pl(kmc))*amas*divbar
          enddo
        enddo
        !$omp end parallel do
      else
        fac_topo_geo(:,:) = 0._wp
      endif

      call solve_poisson(-(fac_topo-fmean_topo), psi_topo)

      if (i_mass_com_vert.ne.1) then
        call solve_poisson(-(fac-fac_topo+fmean_topo-fmean), psi)
      endif

    else

      fac_topo(:,:) = 0._wp
      fac_topo_geo(:,:) = 0._wp
      psi_topo(:,:) = 0._wp
      if (i_mass_com_vert.ne.1) then
        call solve_poisson(-(fac-fmean), psi)
      endif

    endif

    if (i_mass_com_vert.eq.1) then

      ! The two remaining parts, (fac_geo-fac_topo_geo) and its complement, add up to
      ! (fac-fac_topo) identically, so removing each one's own mean - the solvability
      ! condition of the Poisson problem - leaves the three corrections closing the
      ! column budget to the same round-off the unsplit form leaves in fmean.
      fmean_geo = 0._wp
      fmean_ter = 0._wp
      do j=1,jm
        do i=1,im
          fmean_geo = fmean_geo + (fac_geo(i,j)-fac_topo_geo(i,j))
          fmean_ter = fmean_ter + (fac(i,j)-fac_geo(i,j))-(fac_topo(i,j)-fac_topo_geo(i,j))
        enddo
      enddo
      fmean_geo = fmean_geo*aim/real(jm,wp)
      fmean_ter = fmean_ter*aim/real(jm,wp)

      call solve_poisson(-(fac_geo-fac_topo_geo-fmean_geo), psi_geo)
      call solve_poisson(-((fac-fac_geo)-(fac_topo-fac_topo_geo)-fmean_ter), psi_ter)
      psi(:,:) = 0._wp

    else

      psi_geo(:,:) = 0._wp
      psi_ter(:,:) = 0._wp

    endif

    fax_psi(:,:,:)      = 0._wp
    fay_psi(:,:,:)      = 0._wp
    fax_psi_topo(:,:,:) = 0._wp
    fay_psi_topo(:,:,:) = 0._wp

    ! zonal component of the corrective flux, on u-points
    !$omp parallel do private(i,j,k,n,imi,fcx,fcxt,fcx_geo,fcx_ter,fsum,dfx,dp_c,pc1,pc2,pl1,pl2,dp_l,ptrx)
    do j=1,jm
      do i=1,im
        imi = modulo(i - 2, im) + 1
        fcx = dy/dxt(j)*(psi(imi,j)-psi(i,j))
        fcx_geo = 0._wp
        fcx_ter = 0._wp
        if (i_mass_com_vert.eq.1) then
          fcx_geo = dy/dxt(j)*(psi_geo(imi,j)-psi_geo(i,j))
          fcx_ter = dy/dxt(j)*(psi_ter(imi,j)-psi_ter(i,j))
        endif
        if (i_mass_com_topo.eq.1) then
          fcxt = dy/dxt(j)*(psi_topo(imi,j)-psi_topo(i,j))
          ! topographic part: terrain following band, from dps_com_topo above the
          ! local surface up to the tropopause. The correction is uniform in
          ! pressure across the band, so its centroid is the band midpoint: it
          ! follows the barrier down only where the surface is high and the band is
          ! correspondingly short (~0.34 for a Tibetan column), and sits in the mid
          ! troposphere over the ocean (~0.52-0.56), well below the dp_com slab.
          ptrx = 0.5_wp*(pzsa(i,j)+pzsa(imi,j))
          pc2 = ptrop(j)
          pc1 = max(ptrx-dps_com_topo,pc2+0.05_wp)
          dp_c = 0._wp
          do k=1,km
            pl1 = max(min(pl(k),pc1),pl(k+1))
            pl2 = min(max(pl(k+1),pc2),pl(k))
            dp_c = dp_c + (pl1-pl2)*rdpl(k)*dplx(i,j,k)
          enddo
          if (dp_c.gt.0._wp) then
            do k=1,km
              pl1 = max(min(pl(k),pc1),pl(k+1))
              pl2 = min(max(pl(k+1),pc2),pl(k))
              dp_l = pl1-pl2
              dfx = fcxt*dp_l*rdpl(k)*dplx(i,j,k)/dp_c
              fax(i,j,k) = fax(i,j,k) + dfx
              fax_psi_topo(i,j,k) = dfx
            enddo
          elseif (plx(i,j).gt.0._wp) then
            ! degenerate column, fall back to the column mean
            do k=1,km
              dfx = fcxt*dplx(i,j,k)/plx(i,j)
              fax(i,j,k) = fax(i,j,k) + dfx
              fax_psi_topo(i,j,k) = dfx
            enddo
          endif
        endif
        if (i_mass_com_vert.eq.1) then
          ! Each component's residual is returned through the same layers, and in
          ! the same proportion, that the component's own flux occupies at this
          ! face. The correction then carries the same dry static energy and
          ! moisture per unit mass as the flux that created the imbalance, so it is
          ! neutral in both by construction. For the barotropic geostrophic wind
          ! |fax_c(..,1)| is proportional to ftrop*dplx, i.e. this is the layer mass
          ! weighting; for the thermal wind it follows its own top-heavy shear.
          do n=1,2
            if (n.eq.1) then
              fcx = fcx_geo
            else
              fcx = fcx_ter
            endif
            fsum = 0._wp
            do k=1,km
              fsum = fsum + abs(fax_c(i,j,k,n))
            enddo
            if (fsum.gt.0._wp) then
              do k=1,km
                dfx = fcx*abs(fax_c(i,j,k,n))/fsum
                fax(i,j,k) = fax(i,j,k) + dfx
                fax_psi(i,j,k) = fax_psi(i,j,k) + dfx
              enddo
            elseif (plx(i,j).gt.0._wp) then
              ! no wind of this component in the column, fall back to the column mean
              do k=1,km
                dfx = fcx*dplx(i,j,k)/plx(i,j)
                fax(i,j,k) = fax(i,j,k) + dfx
                fax_psi(i,j,k) = fax_psi(i,j,k) + dfx
              enddo
            endif
          enddo
        else
          ! Confine the correction to a slab of thickness dp_com, hanging below the
          ! tropopause with i_mass_com_vert==0 and below the fixed level ptop_com with
          ! i_mass_com_vert==2. Either way the depth is dp_com, so the band cannot
          ! collapse the way a fixed pressure PAIR does once the tropopause drops
          ! below its lower edge.
          if (i_mass_com_vert.eq.2) then
            pc2 = ptop_com
          else
            pc2 = ptrop(j)
          endif
          pc1 = pc2+dp_com
          dp_c = 0._wp
          do k=1,km
            pl1 = max(min(pl(k),pc1),pl(k+1))
            pl2 = min(max(pl(k+1),pc2),pl(k))
            dp_c = dp_c + (pl1-pl2)*rdpl(k)*dplx(i,j,k)
          enddo
          if (dp_c.gt.0._wp) then
            do k=1,km
              pl1 = max(min(pl(k),pc1),pl(k+1))
              pl2 = min(max(pl(k+1),pc2),pl(k))
              dp_l = pl1-pl2
              dfx = fcx*dp_l*rdpl(k)*dplx(i,j,k)/dp_c
              fax(i,j,k) = fax(i,j,k) + dfx
              fax_psi(i,j,k) = fax_psi(i,j,k) + dfx
            enddo
          endif
        endif
      enddo
      ! periodic boundary condition
      do k=1,km
        fax(imc,j,k) = fax(1,j,k)
        fax_psi(imc,j,k) = fax_psi(1,j,k)
        fax_psi_topo(imc,j,k) = fax_psi_topo(1,j,k)
      enddo
    enddo
    !$omp end parallel do

    ! meridional component, on v-points; no flux through the poles
    !$omp parallel do collapse(2) private(i,j,k,n,fcy,fcyt,fcy_geo,fcy_ter,fsum,dfy,dp_c,pc1,pc2,pl1,pl2,dp_l,ptry)
    do j=2,jm
      do i=1,im
        fcy = dxu(j)/dy*(psi(i,j)-psi(i,j-1))
        fcy_geo = 0._wp
        fcy_ter = 0._wp
        if (i_mass_com_vert.eq.1) then
          fcy_geo = dxu(j)/dy*(psi_geo(i,j)-psi_geo(i,j-1))
          fcy_ter = dxu(j)/dy*(psi_ter(i,j)-psi_ter(i,j-1))
        endif
        if (i_mass_com_topo.eq.1) then
          fcyt = dxu(j)/dy*(psi_topo(i,j)-psi_topo(i,j-1))
          ptry = 0.5_wp*(pzsa(i,j)+pzsa(i,j-1))
          pc2 = 0.5_wp*(ptrop(j)+ptrop(j-1))
          pc1 = max(ptry-dps_com_topo,pc2+0.05_wp)
          dp_c = 0._wp
          do k=1,km
            pl1 = max(min(pl(k),pc1),pl(k+1))
            pl2 = min(max(pl(k+1),pc2),pl(k))
            dp_c = dp_c + (pl1-pl2)*rdpl(k)*dply(i,j,k)
          enddo
          if (dp_c.gt.0._wp) then
            do k=1,km
              pl1 = max(min(pl(k),pc1),pl(k+1))
              pl2 = min(max(pl(k+1),pc2),pl(k))
              dp_l = pl1-pl2
              dfy = fcyt*dp_l*rdpl(k)*dply(i,j,k)/dp_c
              fay(i,j,k) = fay(i,j,k) + dfy
              fay_psi_topo(i,j,k) = dfy
            enddo
          elseif (ply(i,j).gt.0._wp) then
            do k=1,km
              dfy = fcyt*dply(i,j,k)/ply(i,j)
              fay(i,j,k) = fay(i,j,k) + dfy
              fay_psi_topo(i,j,k) = dfy
            enddo
          endif
        endif
        if (i_mass_com_vert.eq.1) then
          ! as for the zonal component, each part returned through its own flux
          do n=1,2
            if (n.eq.1) then
              fcy = fcy_geo
            else
              fcy = fcy_ter
            endif
            fsum = 0._wp
            do k=1,km
              fsum = fsum + abs(fay_c(i,j,k,n))
            enddo
            if (fsum.gt.0._wp) then
              do k=1,km
                dfy = fcy*abs(fay_c(i,j,k,n))/fsum
                fay(i,j,k) = fay(i,j,k) + dfy
                fay_psi(i,j,k) = fay_psi(i,j,k) + dfy
              enddo
            elseif (ply(i,j).gt.0._wp) then
              do k=1,km
                dfy = fcy*dply(i,j,k)/ply(i,j)
                fay(i,j,k) = fay(i,j,k) + dfy
                fay_psi(i,j,k) = fay_psi(i,j,k) + dfy
              enddo
            endif
          enddo
        else
          if (i_mass_com_vert.eq.2) then
            pc2 = ptop_com
          else
            pc2 = 0.5_wp*(ptrop(j)+ptrop(j-1))
          endif
          pc1 = pc2+dp_com
          dp_c = 0._wp
          do k=1,km
            pl1 = max(min(pl(k),pc1),pl(k+1))
            pl2 = min(max(pl(k+1),pc2),pl(k))
            dp_c = dp_c + (pl1-pl2)*rdpl(k)*dply(i,j,k)
          enddo
          if (dp_c.gt.0._wp) then
            do k=1,km
              pl1 = max(min(pl(k),pc1),pl(k+1))
              pl2 = min(max(pl(k+1),pc2),pl(k))
              dp_l = pl1-pl2
              dfy = fcy*dp_l*rdpl(k)*dply(i,j,k)/dp_c
              fay(i,j,k) = fay(i,j,k) + dfy
              fay_psi(i,j,k) = fay_psi(i,j,k) + dfy
            enddo
          endif
        endif
      enddo
    enddo
    !$omp end parallel do



    if (niter.eq.1) then

      ! total 3-D wind on T-points, interpolate to levels

      !$omp parallel do collapse(2) private(i,j,ipl,k,kpl,ug_b,vg_b,u3k,v3k,u3kp1,v3kp1,faz)
      do j=1,jm
        do i=1,im

          ipl = modulo(i,im) + 1

          ug_b = ugb(i,j)    
          vg_b = vgb(i,j)
          u3(i,j,1) = 0.5_wp*(ua(i,j,1)+ua(ipl,j,1)) + ug_b
          v3(i,j,1) = 0.5_wp*(va(i,j,1)+va(i,min(j+1,jm),1)) + vg_b
          do k=1,km-1     
            u3k   = 0.5_wp*(ua(i,j,k)+ua(ipl,j,k)) + ug_b + uter(i,j,k)
            v3k   = 0.5_wp*(va(i,j,k)+va(i,min(j+1,jm),k)) + vg_b + vter(i,j,k)
            u3kp1 = 0.5_wp*(ua(i,j,k+1)+ua(ipl,j,k+1)) + ug_b + uter(i,j,k+1)
            v3kp1 = 0.5_wp*(va(i,j,k+1)+va(i,min(j+1,jm),k+1)) + vg_b + vter(i,j,k+1)
            u3(i,j,k+1) = 0.5_wp*(u3k+u3kp1)
            v3(i,j,k+1) = 0.5_wp*(v3k+v3kp1) 
          enddo

          ! z-component and vertical velocity        

          faz = 0._wp
          w3(i,j,1) = 0._wp

          do k=1,km 
            kpl = min(k+1,km)
            faz = faz + fax(i,j,k)-fax(i+1,j,k) + fay(i,j+1,k)-fay(i,j,k)
            ! vertical velocity
            w3(i,j,k+1) = faz/(sqr(i,j)*ra*pl(kpl))    ! kg/s /m2 *m3/kg = m/s
          enddo

        enddo
      enddo
      !$omp end parallel do

      ! zonal mean 500 hPa zonal wind
      do j=1,jm
        uz500(j) = 0._wp
        do i=1,im
          uz500(j) = uz500(j)+u3(i,j,k500)/im
        enddo
      enddo

    endif


    return

  end subroutine u3d


end module u3d_mod
