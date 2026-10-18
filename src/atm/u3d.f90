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
  use atm_params, only : amas, ra, i_pbl, dpc, dp_com, ptopdyn
  use atm_params, only : i_pbl_az, dp_pbl_az
  use atm_params, only : hatm, h_pbl_az
  use atm_params, only : aslp_pbl_az_1, aslp_pbl_az_2
  use atm_params, only : i_mass_com, dps_com_topo, dp_com_topo, dps_com, gad
  use atm_params, only : i_dp_com_topo, aslp_com_topo_1, aslp_com_topo_2, p0
  use atm_params, only : c_uter_pol, c_uter_eq, c_uter_eq_f, c_uter_pol_f
  use atm_grid, only : im, imc, jm, jmc, km, kmc, k500, k700, dxt, dxu, dy, zl, sqr, aim
  use atm_grid, only : fcort, cost, sint
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
  subroutine u3d(niter, pzsa, ptrop, ugb, vgb, uab, vab, t3, gams, aslp, &
        ua, va, uter, vter, uterf, vterf, u3, u3f, v3, w3, w3_nt, uz500, &
        fax, faxo, fay, fayo, fac, fac_topo, psi, psi_topo, &
        fax_psi, fay_psi, fax_psi_topo, fay_psi_topo, froude, &
        dp_ret_topo, p_ret_topo, dp_ret_ageo, p_ret_ageo)

    implicit none

    integer,  intent(in   ) :: niter
    real(wp), intent(in   ) :: pzsa(:,:)
    real(wp), intent(in   ) :: ptrop(:)
    real(wp), intent(in   ) :: ugb(:,:)
    real(wp), intent(in   ) :: vgb(:,:)
    real(wp), intent(in   ) :: uab(:,:)
    real(wp), intent(in   ) :: vab(:,:)
    real(wp), intent(in   ) :: t3(:,:,:)
    real(wp), intent(in   ) :: gams(:,:)
    real(wp), intent(in   ) :: aslp(:,:)   ! Pa, azonal sea level pressure

    real(wp), intent(out  ) :: ua(:,:,:)
    real(wp), intent(out  ) :: va(:,:,:)
    real(wp), intent(out  ) :: uter(:,:,:)
    real(wp), intent(out  ) :: vter(:,:,:)
    real(wp), intent(out  ) :: uterf(:,:,:)
    real(wp), intent(out  ) :: vterf(:,:,:)
    real(wp), intent(inout) :: u3(:,:,:)
    real(wp), intent(inout) :: u3f(:,:,:)   !! as u3, but with its own damping of the thermal wind towards the poles (c_uter_pol_f) and the equator (c_uter_eq_f)
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
    ! flux potential of its topographic part, zero unless i_mass_com==1 or 3
    real(wp), intent(out  ) :: psi_topo(:,:)
    ! the corrective mass flux itself, kept per level so that adifa can form the
    ! dry static energy convergence it carries with the same upstream values the
    ! advection uses; _topo is the grad(psi_topo) part, zero unless i_mass_com==1 or 3
    real(wp), intent(out  ) :: fax_psi(:,:,:)
    real(wp), intent(out  ) :: fay_psi(:,:,:)
    real(wp), intent(out  ) :: fax_psi_topo(:,:,:)
    real(wp), intent(out  ) :: fay_psi_topo(:,:,:)
    ! Froude number U/(N*H) of the low-level flow against the local crest (diagnostic)
    real(wp), intent(out  ) :: froude(:,:)
    ! hPa, thickness and mid pressure of the bands in which the topographic correction and the
    ! azonal ageostrophic PBL flux are returned (diagnostics)
    real(wp), intent(out  ) :: dp_ret_topo(:,:)
    real(wp), intent(out  ) :: p_ret_topo(:,:)
    real(wp), intent(out  ) :: dp_ret_ageo(:,:)
    real(wp), intent(out  ) :: p_ret_ageo(:,:)

    integer :: i, j, k, ipl, imi, jmi, kpl
    real(wp) :: pzx, pzy, dp_c, pbl_t, pbl_u, dp, pc1, pc2, dp_l, pl1, pl2, fxpbl, fypbl, uabc, vabc, ctv
    real(wp) :: faz, faz_nt, fxw, fxe, fys, fyn
    real(wp) :: fcx, fcy, fcxt, fcyt, fmean, fmean_topo, divbar, ptrx, ptry, ptrj
    real(wp) :: dfx, dfy
    real(wp), dimension(imc,jm) :: fxcol
    real(wp), dimension(im,jmc) :: fycol
    real(wp) :: u_g, v_g, ug_b, vg_b
    real(wp) :: u3k, u3kp1, v3k, v3kp1
    real(wp) :: dplxdy, dplydxu
    real(wp) :: c_damp_eq, c_damp_pol
    real(wp), dimension(im,jm,km) :: uteru   ! thermal wind in layers with its own (weaker) polar and equatorial damping, for u3f

    real(wp), dimension(kmc) :: uteri
    real(wp), dimension(kmc) :: vteri
    real(wp), dimension(kmc) :: uterif, vterif   ! thermal wind on levels for the EKE production
    real(wp) :: pst

    real(wp), dimension(imc,jm) :: uabm
    real(wp), dimension(im,jmc) :: vabm
    real(wp) :: ab_mean, az_mean, ub, vb
    real(wp), dimension(im,jm) :: fdpt, hpblt
    real(wp), dimension(imc,jm) :: dpzx, fzx
    real(wp), dimension(im,jmc) :: dpzy, fzy
    real(wp), dimension(imc,jm) :: dphx, uazf
    real(wp), dimension(im,jmc) :: dphy, vazf
    real(wp) :: pzf
    ! Froude number diagnostic
    real(wp), dimension(im,jm) :: pcr
    real(wp) :: n_bv, u_k, h_cr
    real(wp) :: fdct(im,jm), dp_t, f_t
    real(wp), parameter :: froude_max = 10._wp


    !-------------------------------------------------------
    ! Azonal part of the PBL Ekman mass flux 
    !-------------------------------------------------------
    ! The zonal mean of the PBL ageostrophic flux is the model's mean meridional circulation and
    ! keeps the fixed PBL (pblt/pblu) and the return below the tropopause. Its deviation from the
    ! zonal mean (the stationary overturning cells under the azonal highs and lows) is carried by a
    ! PBL of height h, so its mass flux is (uab-zonal mean)*(PBL depth), and it is returned in a band
    ! of depth dp_pbl_az right above that PBL: a shallow cold-core high is compensated in the lower
    ! troposphere, not by subsidence through the whole column. The response is set by the height of
    ! the return band (PBL depth + dp_pbl_az/2)
    ! h is the constant h_pbl_az
    ! The azonal flux is re-centred on every row, so the zonal mean PBL flux and the mean meridional circulation, are unchanged
    ! With i_pbl_az=2 the depth of the return band follows the azonal sea level pressure: dp_pbl_az
    ! under highs (aslp above aslp_pbl_az_1), up to the tropopause under lows (aslp below
    ! aslp_pbl_az_2), linear in between, so a shallow cold-core high is compensated right above
    ! the PBL and the inflow into a heat or convective low through the whole troposphere.
    ! With i_pbl_az=3 the band keeps the depth dp_pbl_az and is lifted instead: right above the
    ! PBL under highs, its top at the tropopause under lows, linear in between, so the return
    ! moves from the lower to the upper troposphere without spreading out.
    ! fdpt is the deep fraction on T-points; faces take the mean of their two cells. dpzx/dpzy is
    ! the depth of the return band on the faces, fzx/fzy the fraction of the possible lift.
    if (i_pbl_az.ge.1) then

      do j=1,jm
        do i=1,im
          ! PBL height on T-points, m
          hpblt(i,j) = h_pbl_az
          ! deep fraction of the return band on T-points, from the azonal sea level pressure
          if (i_pbl_az.eq.2 .or. i_pbl_az.eq.3) then
            fdpt(i,j) = (aslp_pbl_az_1-aslp(i,j))/(aslp_pbl_az_1-aslp_pbl_az_2)
            fdpt(i,j) = min(1._wp,max(0._wp,fdpt(i,j)))
          else
            fdpt(i,j) = 0._wp
          endif
        enddo
      enddo

      ! x-component on u-points
      do j=1,jm
        ab_mean = sum(uab(1:im,j))/real(im,wp)
        do i=1,im
          imi = modulo(i - 2, im) + 1
          pzf = 0.5_wp*(pzsa(i,j)+pzsa(imi,j))
          dphx(i,j) = pzf*(1._wp-exp(-0.5_wp*(hpblt(i,j)+hpblt(imi,j))/hatm))
          if (i_pbl_az.eq.2) then
            dpzx(i,j) = dp_pbl_az + 0.5_wp*(fdpt(i,j)+fdpt(imi,j))*max(0._wp,1._wp-dp_pbl_az)
            fzx(i,j)  = 0._wp
          else
            dpzx(i,j) = dp_pbl_az
            fzx(i,j)  = 0.5_wp*(fdpt(i,j)+fdpt(imi,j))
          endif
          uabm(i,j) = ab_mean
          uazf(i,j) = (uab(i,j)-uabm(i,j))*dphx(i,j)
        enddo
        az_mean = sum(uazf(1:im,j))/real(im,wp)
        do i=1,im
          uazf(i,j) = uazf(i,j)-az_mean
        enddo
        uabm(imc,j) = uabm(1,j)
        dphx(imc,j) = dphx(1,j)
        dpzx(imc,j) = dpzx(1,j)
        fzx(imc,j)  = fzx(1,j)
        uazf(imc,j) = uazf(1,j)
      enddo

      ! y-component on v-points; no flux through the poles
      vabm(:,1)   = vab(:,1)
      vabm(:,jmc) = vab(:,jmc)
      dphy(:,1)   = 0._wp
      dphy(:,jmc) = 0._wp
      dpzy(:,1)   = dp_pbl_az
      dpzy(:,jmc) = dp_pbl_az
      fzy(:,1)    = 0._wp
      fzy(:,jmc)  = 0._wp
      vazf(:,1)   = 0._wp
      vazf(:,jmc) = 0._wp
      do j=2,jm
        ab_mean = sum(vab(1:im,j))/real(im,wp)
        do i=1,im
          pzf = 0.5_wp*(pzsa(i,j)+pzsa(i,j-1))
          dphy(i,j) = pzf*(1._wp-exp(-0.5_wp*(hpblt(i,j)+hpblt(i,j-1))/hatm))
          if (i_pbl_az.eq.2) then
            dpzy(i,j) = dp_pbl_az + 0.5_wp*(fdpt(i,j)+fdpt(i,j-1))*max(0._wp,1._wp-dp_pbl_az)
            fzy(i,j)  = 0._wp
          else
            dpzy(i,j) = dp_pbl_az
            fzy(i,j)  = 0.5_wp*(fdpt(i,j)+fdpt(i,j-1))
          endif
          vabm(i,j) = ab_mean
          vazf(i,j) = (vab(i,j)-vabm(i,j))*dphy(i,j)
        enddo
        az_mean = sum(vazf(1:im,j))/real(im,wp)
        do i=1,im
          vazf(i,j) = vazf(i,j)-az_mean
        enddo
      enddo

    endif


    !$omp parallel do collapse(2) private(i,j,k,ipl,imi,jmi,pzx,pzy,dp,pbl_t,pbl_u,pc1,pc2,dp_l,pl1,pl2,fxpbl,fypbl,uabc,vabc,uteri,vteri,uterif,vterif,pst,ctv,c_damp_eq,c_damp_pol) &
    !$omp private(ub,vb,pzf)
    do j=1,jm
      do i=1,im

        ipl = modulo(i,im) + 1
        imi = modulo(i - 2, im) + 1
        jmi = max(1,j-1)

        if (i_pbl_az.ge.1) then
          ub = uabm(i,j)
          vb = vabm(i,j)
        else
          ub = uab(i,j)
          vb = vab(i,j)
        endif

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
          ua(i,j,k) = ub * dp
        enddo

        ! compensatory velocity in the upper troposphere
        if (i_pbl.eq.1) then
          fxpbl = ub*(1._wp-pblt(j))
        else if (i_pbl.eq.2) then
          fxpbl = ub*max(0._wp,pzx-pblt(j))
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

        ! azonal part in the diagnostic PBL and its return right above it; velocities are
        ! referred to the face layer mass, so outflow and return balance exactly
        if (i_pbl_az.ge.1) then
        if (dphx(i,j).gt.0._wp) then
          pzf = 0.5_wp*(pzsa(i,j)+pzsa(imi,j))
          pc1 = pzf
          pc2 = pzf-dphx(i,j)
          do k=1,km
            pl1 = max(min(pl(k),pc1),pl(k+1))
            pl2 = min(max(pl(k+1),pc2),pl(k))
            dp_l = pl1-pl2
            if (dp_l.gt.0._wp .and. dplx(i,j,k).gt.0._wp) ua(i,j,k) = ua(i,j,k)+uazf(i,j)/dphx(i,j)*dp_l*amas/dplx(i,j,k)
          enddo
          pc1 = pzf-dphx(i,j)
          ! lift of the band (i_pbl_az=3), at most until its top reaches the tropopause
          pc1 = pc1 - fzx(i,j)*max(0._wp,pc1-dpzx(i,j)-ptrop(j))
          pc2 = max(pc1-dpzx(i,j),ptrop(j))
          if (pc1-pc2.gt.0._wp) then
            do k=1,km
              pl1 = max(min(pl(k),pc1),pl(k+1))
              pl2 = min(max(pl(k+1),pc2),pl(k))
              dp_l = pl1-pl2
              if (dp_l.gt.0._wp .and. dplx(i,j,k).gt.0._wp) ua(i,j,k) = ua(i,j,k)-uazf(i,j)/(pc1-pc2)*dp_l*amas/dplx(i,j,k)
            enddo
          endif
        endif
        endif

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
          va(i,j,k) = vb * dp
        enddo

        ! compensatory velocity in the upper troposphere
        if (i_pbl.eq.1) then
          fypbl = vb*(1._wp-pblu(j))
        else if (i_pbl.eq.2) then
          fypbl = vb*max(0._wp,pzy-pblu(j))
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

        ! azonal part in the diagnostic PBL and its return right above it
        if (i_pbl_az.ge.1) then
        if (dphy(i,j).gt.0._wp) then
          pzf = 0.5_wp*(pzsa(i,j)+pzsa(i,jmi))
          pc1 = pzf
          pc2 = pzf-dphy(i,j)
          do k=1,km
            pl1 = max(min(pl(k),pc1),pl(k+1))
            pl2 = min(max(pl(k+1),pc2),pl(k))
            dp_l = pl1-pl2
            if (dp_l.gt.0._wp .and. dply(i,j,k).gt.0._wp) va(i,j,k) = va(i,j,k)+vazf(i,j)/dphy(i,j)*dp_l*amas/dply(i,j,k)
          enddo
          pc1 = pzf-dphy(i,j)
          pc1 = pc1 - fzy(i,j)*max(0._wp,pc1-dpzy(i,j)-0.5_wp*(ptrop(j)+ptrop(jmi)))
          pc2 = max(pc1-dpzy(i,j),0.5_wp*(ptrop(j)+ptrop(jmi)))
          if (pc1-pc2.gt.0._wp) then
            do k=1,km
              pl1 = max(min(pl(k),pc1),pl(k+1))
              pl2 = min(max(pl(k+1),pc2),pl(k))
              dp_l = pl1-pl2
              if (dp_l.gt.0._wp .and. dply(i,j,k).gt.0._wp) va(i,j,k) = va(i,j,k)-vazf(i,j)/(pc1-pc2)*dp_l*amas/dply(i,j,k)
            enddo
          endif
        endif
        endif

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
          do k=1,km
            uter(i,j,k) = 0.5_wp*(uteri(k)+uteri(k+1)) * c_damp_eq * c_damp_pol 
            vter(i,j,k) = 0.5_wp*(vteri(k)+vteri(k+1)) * c_damp_eq * c_damp_pol 
            ! thermal wind without polar damping needed for EKE production
            uterf(i,j,k) = 0.5_wp*(uterif(k)+uterif(k+1)) * c_damp_eq 
            vterf(i,j,k) = 0.5_wp*(vterif(k)+vterif(k+1)) * c_damp_eq
            ! thermal wind for u3f: separate polar (c_uter_pol_f) and equatorial (c_uter_eq_f) damping; the
            ! weak equatorial damping keeps the subtropical jet but removes the fcormin-limited thermal wind
            ! of the deep tropics
            uteru(i,j,k) = 0.5_wp*(uteri(k)+uteri(k+1)) * min(1._wp, c_uter_eq_f*sint(j)**2) * min(1._wp, c_uter_pol_f*cost(j)**2)
          enddo

        endif

      enddo
    enddo
    !$omp end parallel do

    !-------------------------------------------------------
    ! advective mass transport
    !-------------------------------------------------------

    !$omp parallel do collapse(2) private(i,j,k,imi,u_g,v_g,dplxdy,dplydxu)
    do k=1,km
      do j=1,jm

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
          endif
        enddo

      enddo
    enddo 
    !$omp end parallel do

    !-------------------------------------------------------
    ! crest of the 3x3 neighbourhood and Froude number of the low-level flow
    !-------------------------------------------------------
    ! pcr is the surface pressure of the highest surface of the 3x3 neighbourhood,
    ! the crest the low-level air of the column runs into. The Froude number
    ! Fr = U/(N*H), H = crest height above the local surface, U = geostrophic plus
    ! thermal wind of the layer holding the surface, N from gams, is a diagnostic of
    ! the flow regime: below Fr~1 the low-level flow is blocked, and the blocked depth
    ! above the local surface is H*max(0,1-Fr). Capped at froude_max where the column
    ! is its own crest.
    do j=1,jm
      do i=1,im
        pcr(i,j) = pzsa(i,j)
        do jmi=max(1,j-1),min(jm,j+1)
          do ipl=-1,1
            pcr(i,j) = min(pcr(i,j),pzsa(modulo(i+ipl-1,im)+1,jmi))
          enddo
        enddo
        h_cr = hatm*log(pzsa(i,j)/pcr(i,j))
        n_bv = sqrt(max(1.e-6_wp, g/T0*(gad-gams(i,j))))
        k = 1
        do while (k.lt.km .and. pl(k+1).ge.pzsa(i,j))
          k = k+1
        enddo
        u_k = max(0.5_wp, sqrt((ugb(i,j)+uter(i,j,k))**2 + (vgb(i,j)+vter(i,j,k))**2))
        if (n_bv*h_cr*froude_max.gt.u_k) then
          froude(i,j) = u_k/(n_bv*h_cr)
        else
          froude(i,j) = froude_max
        endif
      enddo
    enddo

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
    !   3  as 2 for the non-topographic part, with the topographic term split off and returned in its own terrain following band as in 1.

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

    if (i_mass_com.eq.1 .or. i_mass_com.eq.3) then

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

    fax_psi(:,:,:)      = 0._wp
    fay_psi(:,:,:)      = 0._wp
    fax_psi_topo(:,:,:) = 0._wp
    fay_psi_topo(:,:,:) = 0._wp

    ! Depth of the band of the topographic correction. With i_dp_com_topo=1 it follows the
    ! azonal sea level pressure, as the return of the azonal PBL flux does for i_pbl_az=2:
    ! dp_com_topo under highs (aslp above aslp_com_topo_1), the whole column up to the
    ! tropopause under lows (aslp below aslp_com_topo_2), linear in between. Air lifted over
    ! terrain in a stable cold high stays low, while in a heat or convective low the
    ! forced ascent reaches through the troposphere. A shallow band everywhere dries
    ! summer Central Eurasia (the compensation removes moist, low-DSE air at 1-2 km), a
    ! deep band everywhere undoes the LGM gain on the cold SE flank of the Laurentide.
    ! With i_dp_com_topo=2 the band keeps the depth dp_com_topo and is lifted instead: from
    ! dps_com_topo above the surface under highs to its top at the tropopause under lows,
    ! linear in between.
    ! fdct is the deep fraction on T-points; faces take the mean of their two cells.
    if (i_dp_com_topo.eq.1 .or. i_dp_com_topo.eq.2) then
      do j=1,jm
        do i=1,im
          fdct(i,j) = (aslp_com_topo_1-aslp(i,j))/(aslp_com_topo_1-aslp_com_topo_2)
          fdct(i,j) = min(1._wp,max(0._wp,fdct(i,j)))
        enddo
      enddo
    else
      fdct(:,:) = 0._wp
    endif

    !-------------------------------------------------------
    ! diagnostics: bands in which the mass is returned, on T-points, hPa
    !-------------------------------------------------------
    ! Topographic correction: the terrain following band from dps_com_topo above the
    ! surface, of depth dp_com_topo (deepened under lows with i_dp_com_topo=1, lifted with
    ! i_dp_com_topo=2), capped at the tropopause; zero unless i_mass_com=1 or 3.
    ! Azonal ageostrophic PBL flux: with i_pbl_az>=1 the band right above the PBL, of depth
    ! dp_pbl_az (deepened under lows with i_pbl_az=2, lifted with i_pbl_az=3), capped at the tropopause; with
    ! i_pbl_az=0 the slab of depth dpc below the tropopause, with the zonal mean PBL flux.
    do j=1,jm
      do i=1,im
        if (i_mass_com.eq.1 .or. i_mass_com.eq.3) then
          pc1 = max(pzsa(i,j)-dps_com_topo,ptrop(j)+0.05_wp)
          if (i_dp_com_topo.eq.2) then
            dp_t = dp_com_topo
            pc1 = pc1 - fdct(i,j)*max(0._wp,pc1-dp_t-ptrop(j))
          else
            dp_t = dp_com_topo + fdct(i,j)*max(0._wp,1._wp-dp_com_topo)
          endif
          pc2 = max(pc1-dp_t,ptrop(j))
          dp_ret_topo(i,j) = (pc1-pc2)*p0/100._wp
          p_ret_topo(i,j)  = 0.5_wp*(pc1+pc2)*p0/100._wp
        else
          dp_ret_topo(i,j) = 0._wp
          p_ret_topo(i,j)  = 0._wp
        endif
        if (i_pbl_az.ge.1) then
          pc1 = pzsa(i,j)*exp(-hpblt(i,j)/hatm)
          if (i_pbl_az.eq.3) then
            pc1 = pc1 - fdpt(i,j)*max(0._wp,pc1-dp_pbl_az-ptrop(j))
            pc2 = max(pc1-dp_pbl_az,ptrop(j))
          else
            pc2 = max(pc1-(dp_pbl_az+fdpt(i,j)*max(0._wp,1._wp-dp_pbl_az)),ptrop(j))
          endif
        else
          pc2 = ptrop(j)
          pc1 = pc2+dpc
        endif
        dp_ret_ageo(i,j) = (pc1-pc2)*p0/100._wp
        p_ret_ageo(i,j)  = 0.5_wp*(pc1+pc2)*p0/100._wp
      enddo
    enddo

    ! zonal component of the corrective flux, on u-points
    !$omp parallel do private(i,j,k,imi,fcx,fcxt,dfx,dp_c,pc1,pc2,pl1,pl2,dp_l,ptrx,dp_t,f_t)
    do j=1,jm
      do i=1,im
        imi = modulo(i - 2, im) + 1
        fcx = dy/dxt(j)*(psi(imi,j)-psi(i,j))
        if (i_mass_com.eq.1 .or. i_mass_com.eq.3) then
          fcxt = dy/dxt(j)*(psi_topo(imi,j)-psi_topo(i,j))
          ! topographic part: terrain following band of thickness dp_com_topo, from
          ! dps_com_topo above the local surface upwards, its top capped at the
          ! tropopause (dp_com_topo>=1 gives the whole band up to the tropopause);
          ! with i_dp_com_topo=1 deepened to the tropopause under lows, with
          ! i_dp_com_topo=2 lifted under lows at constant depth
          ptrx = 0.5_wp*(pzsa(i,j)+pzsa(imi,j))
          pc1 = max(ptrx-dps_com_topo,ptrop(j)+0.05_wp)
          f_t = 0.5_wp*(fdct(i,j)+fdct(imi,j))
          if (i_dp_com_topo.eq.2) then
            dp_t = dp_com_topo
            pc1 = pc1 - f_t*max(0._wp,pc1-dp_t-ptrop(j))
          else
            dp_t = dp_com_topo + f_t*max(0._wp,1._wp-dp_com_topo)
          endif
          pc2 = max(pc1-dp_t,ptrop(j))
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
        if (i_mass_com.eq.2 .or. i_mass_com.eq.3) then
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
    !$omp parallel do collapse(2) private(i,j,k,fcy,fcyt,dfy,dp_c,pc1,pc2,pl1,pl2,dp_l,ptry,ptrj,dp_t,f_t)
    do j=2,jm
      do i=1,im
        fcy = dxu(j)/dy*(psi(i,j)-psi(i,j-1))
        if (i_mass_com.eq.1 .or. i_mass_com.eq.3) then
          fcyt = dxu(j)/dy*(psi_topo(i,j)-psi_topo(i,j-1))
          ptry = 0.5_wp*(pzsa(i,j)+pzsa(i,j-1))
          ptrj = 0.5_wp*(ptrop(j)+ptrop(j-1))
          pc1 = max(ptry-dps_com_topo,ptrj+0.05_wp)
          f_t = 0.5_wp*(fdct(i,j)+fdct(i,j-1))
          if (i_dp_com_topo.eq.2) then
            dp_t = dp_com_topo
            pc1 = pc1 - f_t*max(0._wp,pc1-dp_t-ptrj)
          else
            dp_t = dp_com_topo + f_t*max(0._wp,1._wp-dp_com_topo)
          endif
          pc2 = max(pc1-dp_t,ptrj)
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
        if (i_mass_com.eq.2 .or. i_mass_com.eq.3) then
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
          vg_b = vgb(i,j)
          u3(i,j,1) = 0.5_wp*(ua(i,j,1)+ua(ipl,j,1)) + ug_b
          u3f(i,j,1) = u3(i,j,1)
          v3(i,j,1) = 0.5_wp*(va(i,j,1)+va(i,min(j+1,jm),1)) + vg_b
          do k=1,km-1     
            u3k   = 0.5_wp*(ua(i,j,k)+ua(ipl,j,k)) + ug_b + uter(i,j,k)
            v3k   = 0.5_wp*(va(i,j,k)+va(i,min(j+1,jm),k)) + vg_b + vter(i,j,k)
            u3kp1 = 0.5_wp*(ua(i,j,k+1)+ua(ipl,j,k+1)) + ug_b + uter(i,j,k+1)
            v3kp1 = 0.5_wp*(va(i,j,k+1)+va(i,min(j+1,jm),k+1)) + vg_b + vter(i,j,k+1)
            u3(i,j,k+1) = 0.5_wp*(u3k+u3kp1)
            ! the same with the undamped thermal wind
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
        uz500(j) = 0._wp
        do i=1,im
          uz500(j) = uz500(j)+u3(i,j,k500)/im
        enddo
      enddo

    endif


    return

  end subroutine u3d


end module u3d_mod
