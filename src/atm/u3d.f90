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
  use constants, only : g, T0, r_earth, pi, omega
  use atm_params, only : amas, ra, dpc, i_dpc, dp_com, ptopdyn
  use atm_params, only : i_mass_com, dps_com_topo, dp_com_topo, dps_com
  use atm_params, only : c_uter_pol, c_uter_eq, c_uz500_pol
  use atm_grid, only : im, imc, jm, jmc, km, kmc, k500, k700, dxt, dxu, dy, zl, sqr, aim
  use atm_grid, only : fcort, cost, sint, cosu, sinu
  use atm_grid, only : pl, rdpl, dplt, dplx, dply, dplxo, dplyo, plx, ply, pblt, pblu

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
    complex(dp), dimension(jm) :: rr, xx
    complex(dp) :: bet
    complex(dp), dimension(im/2+1,jm) :: rhs_fft
    real(dp), dimension(im), save :: row
    complex(dp), dimension(im/2+1), save :: row_fft
    type(C_PTR), save :: plan_r2c, plan_c2r
    logical, save :: plans_ready = .false.
    ! The operator depends on the grid only, so the Thomas elimination factors of each
    ! wavenumber are computed once: gam_t, bet_t and the sub-diagonal aa_t
    real(wp), dimension(jm), save :: aa_t
    complex(dp), dimension(jm,im/2+1), save :: gam_t, bet_t


    if (.not.plans_ready) then
      plan_r2c = fftw_plan_dft_r2c_1d(im, row, row_fft, FFTW_ESTIMATE)
      plan_c2r = fftw_plan_dft_c2r_1d(im, row_fft, row, FFTW_ESTIMATE)

      ! metric weights of the flux-form gradient
      do j=1,jm
        wx(j) = dy/dxt(j)
      enddo
      do j=1,jmc
        wy(j) = dxu(j)/dy
      enddo

      do n=1,im/2+1
        nn = n-1
        ! eigenvalue of the zonal second difference for wavenumber nn
        lamb = -4._wp*sin(pi*real(nn,wp)/real(im,wp))**2
        do j=1,jm
          aa(j) = wy(j)
          cc(j) = wy(j+1)
          bb(j) = -(wy(j)+wy(j+1)) + lamb*wx(j)
        enddo
        if (nn.eq.0) then
          bb(1) = 1._wp
          cc(1) = 0._wp
        endif
        ! Thomas algorithm, elimination factors
        bet = bb(1)
        bet_t(1,n) = bet
        gam_t(1,n) = (0._dp,0._dp)
        do j=2,jm
          gam_t(j,n) = cc(j-1)/bet
          bet = bb(j)-aa(j)*gam_t(j,n)
          bet_t(j,n) = bet
        enddo
      enddo
      aa_t(:) = wy(1:jm)

      plans_ready = .true.
    endif

    do j=1,jm
      row(:) = rhs(:,j)
      call fftw_execute_dft_r2c(plan_r2c, row, row_fft)
      rhs_fft(:,j) = row_fft(:)
    enddo

    ! the wavenumbers are independent
    !$omp parallel do private(n,j,rr,xx)
    do n=1,im/2+1
      do j=1,jm
        rr(j) = rhs_fft(n,j)
      enddo
      if (n.eq.1) rr(1) = (0._dp,0._dp)   ! psi pinned at the first row for wavenumber zero
      ! Thomas algorithm, substitution
      xx(1) = rr(1)/bet_t(1,n)
      do j=2,jm
        xx(j) = (rr(j)-aa_t(j)*xx(j-1))/bet_t(j,n)
      enddo
      do j=jm-1,1,-1
        xx(j) = xx(j)-gam_t(j+1,n)*xx(j+1)
      enddo
      rhs_fft(n,:) = xx(:)
    enddo
    !$omp end parallel do

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
        ua, va, uter, vter, uterf, vterf, u3, u3f, v3, w3, w3_nt, uz500, &
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
    real(wp), intent(inout) :: u3f(:,:,:)   !! as u3, but without the polar damping of the thermal wind
    real(wp), intent(inout) :: v3(:,:,:)
    real(wp), intent(inout) :: w3(:,:,:)
    real(wp), intent(inout) :: w3_nt(:,:,:)
    real(wp), intent(inout) :: uz500(:)
    real(wp), intent(out  ) :: fax(:,:,:)
    real(wp), intent(out  ) :: faxo(:,:,:)
    real(wp), intent(out  ) :: fay(:,:,:)
    real(wp), intent(out  ) :: fayo(:,:,:)
    real(wp), intent(out  ) :: fac(:,:)
    real(wp), intent(out  ) :: fac_topo(:,:)
    ! flux potential of the column mass correction, kg/s
    real(wp), intent(out  ) :: psi(:,:)
    ! flux potential of its topographic part, zero unless i_mass_com==1
    real(wp), intent(out  ) :: psi_topo(:,:)
    ! the corrective mass flux itself, kept per level so that adifa can form the
    ! dry static energy convergence it carries with the same upstream values the
    ! advection uses; _topo is the grad(psi_topo) part, zero unless i_mass_com==1
    real(wp), intent(out  ) :: fax_psi(:,:,:)
    real(wp), intent(out  ) :: fay_psi(:,:,:)
    real(wp), intent(out  ) :: fax_psi_topo(:,:,:)
    real(wp), intent(out  ) :: fay_psi_topo(:,:,:)

    integer :: i, j, k, ipl, imi, jmi, kpl
    real(wp) :: pzx, pzy, dp_c, pbl_t, pbl_u, dp, pc1, pc2, dp_l, pl1, pl2, fxpbl, fypbl, uabc, vabc, ctv
    real(wp) :: faz, faz_nt, fxw, fxe, fys, fyn
    real(wp) :: fcx, fcy, fcxt, fcyt, fmean, fmean_topo, divbar, ptrx, ptry, ptrj
    real(wp) :: dfx, dfy
    real(wp), dimension(imc,jm) :: fxcol
    real(wp), dimension(im,jmc) :: fycol
    real(wp) :: u_g, v_g, ug_b, vg_b
    real(wp), dimension(km) :: ftop   ! fraction of each layer below ptopdyn
    real(wp) :: u3k, u3kp1, v3k, v3kp1
    real(wp) :: dplxdy, dplydxu
    real(wp) :: c_damp_eq, c_damp_pol
    real(wp), dimension(im,jm,km) :: uteru   ! thermal wind in layers without polar damping, for u3f

    real(wp), dimension(kmc) :: uteri
    real(wp), dimension(kmc) :: vteri
    real(wp), dimension(kmc) :: uterif, vterif   ! thermal wind on levels for the EKE production
    real(wp), dimension(kmc) :: uterid           ! zonal thermal wind on levels from the polar damped azonal temperature
    real(wp) :: cdamp_uz
    real(wp), dimension(jm)  :: cdamp            ! polar damping factor, applied to the azonal temperature
    real(wp), dimension(jm,kmc) :: t3z           ! zonal mean temperature on levels
    real(wp) :: pst



    ! polar damping factor of the thermal wind and zonal mean temperature
    do j=1,jm
      cdamp(j) = min(1._wp,c_uter_pol*cost(j)**2)
      do k=1,kmc
        t3z(j,k) = sum(t3(:,j,k))*aim
      enddo
    enddo

    !$omp parallel do collapse(2) private(i,j,k,ipl,imi,jmi,pzx,pzy,dp,pbl_t,pbl_u,pc1,pc2,dp_l,pl1,pl2,fxpbl,fypbl,uabc,vabc,uteri,vteri,uterif,vterif,uterid,pst,ctv,c_damp_eq,c_damp_pol)
    do j=1,jm
      do i=1,im

        ipl = modulo(i,im) + 1
        imi = modulo(i - 2, im) + 1
        jmi = max(1,j-1)

        !-------------------------------------------------------
        ! Vertical profile of ageostrophic wind

        ! x-component on u-points

        ! velocity in the PBL
        pzx = 0.5_wp*(pzsa(i,j)+pzsa(imi,j))
        pbl_t = pzx+pblt(j)-1._wp
        do k=1,km
          dp = (pl(k)-pbl_t)*rdpl(k)
          dp = min(1._wp,dp)
          dp = max(0._wp,dp)
          ua(i,j,k) = uab(i,j) * dp
        enddo

        ! compensatory velocity in the upper troposphere
        fxpbl = uab(i,j)*(1._wp-pblt(j))
        pc2 = ptrop(j)
        if (i_dpc.eq.1) then
          ! fraction dpc of the troposphere in pressure, from sea level to the zonal mean tropopause
          pc1 = pc2+dpc*(1._wp-pc2)
          uabc = -fxpbl/(pc1-pc2)
        else
          pc1 = pc2+dpc
          uabc = -fxpbl/dpc
        endif
        do k=1,km
          pl1 = max(min(pl(k),pc1),pl(k+1))
          pl2 = min(max(pl(k+1),pc2),pl(k))
          dp_l = pl1-pl2
          ua(i,j,k) = ua(i,j,k)+uabc*dp_l*rdpl(k)
        enddo

        ! y-component on v-points

        ! velocity in the PBL
        pzy = 0.5_wp*(pzsa(i,j)+pzsa(i,jmi))
        pbl_u = pzy+pblu(j)-1._wp
        do k=1,km
          dp = (pl(k)-pbl_u)*rdpl(k)
          dp = min(1._wp,dp)
          dp = max(0._wp,dp)
          va(i,j,k) = vab(i,j) * dp
        enddo

        ! compensatory velocity in the upper troposphere
        fypbl = vab(i,j)*(1._wp-pblu(j))
        pc2 = 0.5_wp*(ptrop(j)+ptrop(jmi))
        if (i_dpc.eq.1) then
          pc1 = pc2+dpc*(1._wp-pc2)
          vabc = -fypbl/(pc1-pc2)
        else
          pc1 = pc2+dpc
          vabc = -fypbl/dpc
        endif
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
          uteru(i,j,:) = 0._wp

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

          ! thermal wind for the EKE production: a layer adds to the shear only if its centre,
          ! sqrt(pl(k)*pl(k+1)) = exp(-zc/hatm), is above the highest surface of the gradient stencil
          pst = min(pzsa(i,j),pzsa(i,jmi),pzsa(i,j+1),pzsa(ipl,j),pzsa(imi,j))
          uterif(1) = 0._wp
          vterif(1) = 0._wp
          uterif(km+1) = 0._wp
          vterif(km+1) = 0._wp
          do k=1,km-1
            if (sqrt(pl(k)*pl(k+1)).gt.pst) then
              uterif(k+1) = uterif(k)
              vterif(k+1) = vterif(k)
            else
              uterif(k+1) = uterif(k) + (uteri(k+1)-uteri(k))
              vterif(k+1) = vterif(k) + (vteri(k+1)-vteri(k))
            endif
          enddo

          ! thermal wind in layers, dampened at equator and poles
          c_damp_pol = min(1._wp,c_uter_pol*cost(j)**2)
          c_damp_eq  = min(1._wp, c_uter_eq*sint(j)**2)
          ! Polar damping. Scaling the wind by a latitude dependent factor c is the
          ! geostrophic wind of f/c, whose beta, beta*(1+2tan^2) at high latitudes, makes the
          ! flux divergent, and the column mass correction would then return that mass in its
          ! band, damping the azonal temperature. Instead, the azonal temperature of each row
          ! is scaled by c before the gradient, which keeps the azonal thermal wind
          ! non-divergent apart from the true beta, at the cost of an extra zonal wind
          ! -K*T_az*dc/dy. The zonal mean part carries no meridional flux and is scaled by
          ! c(j); v is unchanged, as c(j) factors out of the zonal gradient.
          ! The equatorial damping scales the wind: inside c, K*dc/dy ~ phi/f would stay
          ! finite at the equator and add a tropical zonal wind of several m/s.
          uterid(1) = 0._wp
          uterid(km+1) = 0._wp
          do k=1,km-1
            ctv = (zl(k+1)-zl(k))*g/(T0*fcort(j))
            uterid(k+1) = uterid(k)-ctv*( cdamp(jmi)*(t3(i,jmi,k)-t3z(jmi,k)) - cdamp(j+1)*(t3(i,j+1,k)-t3z(j+1,k)) &
                                        + cdamp(j)*(t3z(jmi,k)-t3z(j+1,k)) )/(2._wp*dy)
          enddo
          do k=1,km
            uter(i,j,k) = 0.5_wp*(uterid(k)+uterid(k+1)) * c_damp_eq
            vter(i,j,k) = 0.5_wp*(vteri(k)+vteri(k+1)) * c_damp_eq * c_damp_pol 
            ! thermal wind without polar damping needed for EKE production
            uterf(i,j,k) = 0.5_wp*(uterif(k)+uterif(k+1)) * c_damp_eq 
            vterf(i,j,k) = 0.5_wp*(vterif(k)+vterif(k+1)) * c_damp_eq
            ! thermal wind without polar damping for u3f
            uteru(i,j,k) = 0.5_wp*(uteri(k)+uteri(k+1)) * c_damp_eq
          enddo

        endif

      enddo
    enddo
    !$omp end parallel do

    !-------------------------------------------------------
    ! advective mass transport
    !-------------------------------------------------------

    ! The geostrophic (ugb,vgb) and the thermal wind (uter,vter) components are carried
    ! up to the top pressure ptopdyn; the layer holding the top carries the part of it
    ! below the top.
    do k=1,km
      ftop(k) = max(0._wp, min(1._wp, (pl(k)-ptopdyn)*rdpl(k)))
    enddo

    !$omp parallel do collapse(2) private(i,j,k,imi,u_g,v_g,dplxdy,dplydxu)
    do k=1,km
      do j=1,jm

        ! x-components

        do i=1,im
          imi = modulo(i - 2, im) + 1
          ! geostrophic zonal wind on u-points, limited to troposphere
          dplxdy = dplx(i,j,k)*dy
          ! geostrophic and thermal wind, averaged from the T-points to the face
          u_g = 0.5_wp*(ugb(imi,j)+ugb(i,j))*ftop(k) + 0.5_wp*(uter(imi,j,k)+uter(i,j,k))*ftop(k)
          ! mass flux
          fax(i,j,k) = (u_g+ua(i,j,k))*dplxdy ! m/s * kg/m2 * m = kg/s
          ! orographic component of mass flux 
          faxo(i,j,k) = u_g*(dplxo(i,j,k)*dy-dplxdy) 
        enddo
        ! periodic boundary conditions
        fax(imc,j,k) = fax(1,j,k)
        faxo(imc,j,k) = faxo(1,j,k)

        ! y-components

        do i=1,im
          if (j.eq.1) then
            ! N and S boundary conditions, no flux
            fay(i,1,k)  = 0._wp
            fayo(i,1,k) = 0._wp
            fay(i,jmc,k)  = 0._wp
            fayo(i,jmc,k) = 0._wp
          else
            ! geostrophic meridional wind on v-points, limited to troposphere
            dplydxu = dply(i,j,k)*dxu(j)
            v_g = 0.5_wp*(vgb(i,j-1)+vgb(i,j))*ftop(k) + 0.5_wp*(vter(i,j-1,k)+vter(i,j,k))*ftop(k)
            ! mass flux
            fay(i,j,k) = (v_g+va(i,j,k))*dplydxu  ! m/s * kg/m2 * m = kg/s
            ! orographic component of mass flux for temperature
            fayo(i,j,k) = v_g*(dplyo(i,j,k)*dxu(j)-dplydxu) 
          endif
        enddo

      enddo
    enddo
    !$omp end parallel do

    !-------------------------------------------------------
    ! per-column mass conservation
    !-------------------------------------------------------
    ! The column mass is diagnosed from topography and does not evolve, so the column-integrated mass flux divergence has to vanish: 
    ! the same fluxes advect dry static energy and moisture, and advecting a uniform tracer must leave it unchanged. 
    ! The ageostrophic wind closes by construction (see above), so what remains is the geostrophic and thermal-wind divergence, 
    ! a large part of which is the topographic term V.grad(p_s).
    !
    ! The residual is removed with a corrective mass flux F_c = grad(psi) obtained from  laplacian(psi) = -(column convergence)
    !
    ! i_mass_com selects the vertical placement of the correction:
    !   0  the whole residual in a slab of depth dp_com hanging below the tropopause, so the band follows the tropopause down at high latitudes.
    !   1  as 0, but the topographic term is separated from the rest and given its own potential, returned in a terrain following band from 
    !      dps_com_topo above the local surface, dp_com_topo thick and capped at the tropopause. 
    !   2  the whole residual, topographic part included, as one correction in a terrain following band from dps_com above the local surface,
    !      dp_com thick and capped at the tropopause, uniform in pressure across the band (the same wind at every height of the band).
    !      With dps_com=0 and dp_com=1 the band is the whole column up to the tropopause, the smallest 3-D correction that closes every column:
    !      minimising sum_k |dF_k|^2/dp_k subject to sum_k dF_k = grad(psi) gives dF_k ~ dp_k
    !      It is also what a prognostic surface pressure would do in the limit of fast barotropic adjustment: 
    !      a p_s anomaly exerts the same pressure gradient force at every height, so the flow that removes it is vertically uniform.

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

    if (i_mass_com.eq.1) then

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

      call solve_poisson(-(fac_topo-fmean_topo), psi_topo)
      call solve_poisson(-(fac-fac_topo+fmean_topo-fmean), psi)

    else

      fac_topo(:,:) = 0._wp
      psi_topo(:,:) = 0._wp
      call solve_poisson(-(fac-fmean), psi)

    endif

    !$omp parallel do collapse(2) private(j,k)
    do k=1,km
      do j=1,jmc
        if (j.le.jm) then
          fax_psi(:,j,k)      = 0._wp
          fax_psi_topo(:,j,k) = 0._wp
        endif
        fay_psi(:,j,k)      = 0._wp
        fay_psi_topo(:,j,k) = 0._wp
      enddo
    enddo
    !$omp end parallel do

    ! zonal component of the corrective flux, on u-points
    !$omp parallel do private(i,j,k,imi,fcx,fcxt,dfx,dp_c,pc1,pc2,pl1,pl2,dp_l,ptrx)
    do j=1,jm
      do i=1,im
        imi = modulo(i - 2, im) + 1
        fcx = dy/dxt(j)*(psi(imi,j)-psi(i,j))
        if (i_mass_com.eq.1) then
          fcxt = dy/dxt(j)*(psi_topo(imi,j)-psi_topo(i,j))
          ! topographic part: terrain following band of thickness dp_com_topo, from
          ! dps_com_topo above the local surface upwards, its top capped at the
          ! tropopause (dp_com_topo>=1 gives the whole band up to the tropopause)
          ptrx = 0.5_wp*(pzsa(i,j)+pzsa(imi,j))
          pc1 = max(ptrx-dps_com_topo,ptrop(j)+0.05_wp)
          pc2 = max(pc1-dp_com_topo,ptrop(j))
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
        if (i_mass_com.eq.2) then
          ! terrain following band of thickness dp_com, from dps_com above the local
          ! surface upwards, its top capped at the tropopause, as for the topographic part
          ptrx = 0.5_wp*(pzsa(i,j)+pzsa(imi,j))
          pc1 = max(ptrx-dps_com,ptrop(j)+0.05_wp)
          pc2 = max(pc1-dp_com,ptrop(j))
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
          elseif (plx(i,j).gt.0._wp) then
            do k=1,km
              dfx = fcx*dplx(i,j,k)/plx(i,j)
              fax(i,j,k) = fax(i,j,k) + dfx
              fax_psi(i,j,k) = fax_psi(i,j,k) + dfx
            enddo
          endif
        else
          ! slab of thickness dp_com hanging below the tropopause
          pc2 = ptrop(j)
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
    !$omp parallel do collapse(2) private(i,j,k,fcy,fcyt,dfy,dp_c,pc1,pc2,pl1,pl2,dp_l,ptry,ptrj)
    do j=2,jm
      do i=1,im
        fcy = dxu(j)/dy*(psi(i,j)-psi(i,j-1))
        if (i_mass_com.eq.1) then
          fcyt = dxu(j)/dy*(psi_topo(i,j)-psi_topo(i,j-1))
          ptry = 0.5_wp*(pzsa(i,j)+pzsa(i,j-1))
          ptrj = 0.5_wp*(ptrop(j)+ptrop(j-1))
          pc1 = max(ptry-dps_com_topo,ptrj+0.05_wp)
          pc2 = max(pc1-dp_com_topo,ptrj)
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
        if (i_mass_com.eq.2) then
          ptry = 0.5_wp*(pzsa(i,j)+pzsa(i,j-1))
          ptrj = 0.5_wp*(ptrop(j)+ptrop(j-1))
          pc1 = max(ptry-dps_com,ptrj+0.05_wp)
          pc2 = max(pc1-dp_com,ptrj)
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
          elseif (ply(i,j).gt.0._wp) then
            do k=1,km
              dfy = fcy*dply(i,j,k)/ply(i,j)
              fay(i,j,k) = fay(i,j,k) + dfy
              fay_psi(i,j,k) = fay_psi(i,j,k) + dfy
            enddo
          endif
        else
          pc2 = 0.5_wp*(ptrop(j)+ptrop(j-1))
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

      !$omp parallel do collapse(2) private(i,j,ipl,k,kpl,ug_b,vg_b,u3k,v3k,u3kp1,v3kp1,faz,faz_nt,fxw,fxe,fys,fyn)
      do j=1,jm
        do i=1,im

          ipl = modulo(i,im) + 1

          ug_b = ugb(i,j)    
          ! on the polar rows v at the T-point is the v of the one interior face (v does not vanish at the pole)
          vg_b = vgb(i,j)
          u3(i,j,1) = 0.5_wp*(ua(i,j,1)+ua(ipl,j,1)) + ug_b
          u3f(i,j,1) = u3(i,j,1)
          v3(i,j,1) = 0.5_wp*(va(i,max(j,2),1)+va(i,min(j+1,jm),1)) + vg_b
          do k=1,km-1     
            u3k   = 0.5_wp*(ua(i,j,k)+ua(ipl,j,k)) + ug_b + uter(i,j,k)
            v3k   = 0.5_wp*(va(i,max(j,2),k)+va(i,min(j+1,jm),k)) + vg_b + vter(i,j,k)
            u3kp1 = 0.5_wp*(ua(i,j,k+1)+ua(ipl,j,k+1)) + ug_b + uter(i,j,k+1)
            v3kp1 = 0.5_wp*(va(i,max(j,2),k+1)+va(i,min(j+1,jm),k+1)) + vg_b + vter(i,j,k+1)
            u3(i,j,k+1) = 0.5_wp*(u3k+u3kp1)
            ! the same without the polar damping of the thermal wind
            u3f(i,j,k+1) = u3(i,j,k+1) + 0.5_wp*(uteru(i,j,k)-uter(i,j,k) + uteru(i,j,k+1)-uter(i,j,k+1))
            v3(i,j,k+1) = 0.5_wp*(v3k+v3kp1) 
          enddo

          ! z-component and vertical velocity        

          ! w3 is the mass flux convergence accumulated from the surface upward, so at
          ! any level it carries the convergence of every layer below it, including the
          ! layers the terrain cuts through. Writing the convergence of a single layer as
          ! div(U*dp) = dp*div(U) + U*grad(dp), the second term is the level-wise form of
          ! the topographic term V.grad(p_s) that is split off as fac_topo above: a wind
          ! crossing an elevation gradient carries a different layer mass in than out. It
          ! is non-zero only where the surface cuts the layer, which is always below the
          ! cloud level kweff, and it is the dominant signal in w3 just above sloping
          ! terrain. w3_nt keeps the dp*div(U) term alone

          faz = 0._wp
          faz_nt = 0._wp
          w3(i,j,1) = 0._wp
          w3_nt(i,j,1) = 0._wp

          do k=1,km 
            kpl = min(k+1,km)
            faz = faz + fax(i,j,k)-fax(i+1,j,k) + fay(i,j+1,k)-fay(i,j,k)
            ! vertical velocity
            w3(i,j,k+1) = faz/(sqr(i,j)*ra*pl(kpl))    ! kg/s /m2 *m3/kg = m/s
            ! face winds, m2/s; a face with no air carries no flux either. ipl is the
            ! index of face i+1 in the face arrays, which are dimensioned im and closed
            ! periodically, so it also covers the face imc = 1
            fxw = 0._wp
            fxe = 0._wp
            fys = 0._wp
            fyn = 0._wp
            if (dplx(i,j,k).gt.0._wp)   fxw = (fax(i,j,k)  -fax_psi_topo(i,j,k))  /dplx(i,j,k)
            if (dplx(ipl,j,k).gt.0._wp) fxe = (fax(i+1,j,k)-fax_psi_topo(i+1,j,k))/dplx(ipl,j,k)
            ! no flux through the poles, and dply is not defined there
            if (j.gt.1) then
              if (dply(i,j,k).gt.0._wp) fys = (fay(i,j,k)-fay_psi_topo(i,j,k))/dply(i,j,k)
            endif
            if (j.lt.jm) then
              if (dply(i,j+1,k).gt.0._wp) fyn = (fay(i,j+1,k)-fay_psi_topo(i,j+1,k))/dply(i,j+1,k)
            endif
            faz_nt = faz_nt + dplt(i,j,k)*(fxw-fxe+fyn-fys)   ! kg/m2 *m2/s = kg/s
            ! vertical velocity without the topographically induced part
            w3_nt(i,j,k+1) = faz_nt/(sqr(i,j)*ra*pl(kpl))
          enddo

        enddo
      enddo
      !$omp end parallel do

      ! zonal mean 500 hPa zonal wind
      do j=1,jm
        ! polar damping of the thermal wind in the wind for the sea level pressure, independent of c_uter_pol
        cdamp_uz = 1._wp
        if (c_uz500_pol.gt.0._wp) cdamp_uz = min(1._wp,c_uz500_pol*cost(j)**2)
        uz500(j) = 0._wp
        do i=1,im
          uz500(j) = uz500(j) + (u3f(i,j,k500) &
            - (1._wp-cdamp_uz)*0.5_wp*(uteru(i,j,k500-1)+uteru(i,j,k500)))/im
        enddo
      enddo

    endif


    return

  end subroutine u3d


end module u3d_mod
