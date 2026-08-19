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

  use atm_params, only : wp, dp
  use constants, only : g, T0, r_earth, pi
  use atm_params, only : amas, ra, i_pbl, dp_com, ptopdyn, i_ptopdyn
  use atm_params, only : vprof_stl, vprof_a, vprof_b
  use atm_params, only : i_mass_com_topo, dps_com_topo, i_mass_com_vert
  use atm_params, only : c_uter_eq, i_uter_damp
  use atm_params, only : l_output_flx3d, l_diag_wcomp
  use atm_grid, only : im, imc, jm, jmc, km, kmc, k500, k700, dxt, dxu, dy, zl, sqr, aim
  use atm_grid, only : fcort, cost, sint, cdamp_pol
  use atm_grid, only : pl, dplx, dply, plx, ply, ptopt, ptopu

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
  !   Function   :  v p r o f _ m a s s
  !   Purpose    :  mass carried by the lower branch below a level
  !
  !   Cumulative mass, in units of the sigma coordinate and per unit surface wind, that the
  !   observed profile carries between the level sig and the surface sig_s,
  !       M/M_max = 1 - (1-zeta)**vprof_a ,   zeta = (sig_s-sig)/(sig_s-sig_stl)
  !   saturating at M_max = (sig_s-sig_stl)/vprof_a at and above the steering level. 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pure function vprof_mass(sig, sig_s, sig_stl) result(m)

    implicit none

    real(wp), intent(in) :: sig
    real(wp), intent(in) :: sig_s
    real(wp), intent(in) :: sig_stl
    real(wp) :: m

    real(wp) :: d, z

    d = sig_s-sig_stl
    if (d.le.0._wp) then
      m = 0._wp
    else
      z = (sig_s-sig)/d
      z = min(1._wp,max(0._wp,z))
      m = d/vprof_a*(1._wp-(1._wp-z)**vprof_a)
    endif

  end function vprof_mass


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  v p r o f _ w r e t
  !   Purpose    :  cumulative shape of the return branch
  !
  !   Runs from 0 at the steering level to 1 at the cell top, following the observed
  !       1 - M/M_max = xi**vprof_b ,   xi = (sig_stl-sig)/(sig_stl-sig_top)
  !   so the return flow itself grows as xi**(vprof_b-1): zero at the steering level and
  !   largest at the cell top. Only the shape is used; u3d scales it to close the column
  !   mass budget exactly on the discrete layer thicknesses.
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pure function vprof_wret(sig, sig_stl, sig_top) result(w)

    implicit none

    real(wp), intent(in) :: sig
    real(wp), intent(in) :: sig_stl
    real(wp), intent(in) :: sig_top
    real(wp) :: w

    real(wp) :: x

    if (sig_stl.le.sig_top) then
      w = 0._wp
    else
      x = (sig_stl-sig)/(sig_stl-sig_top)
      x = min(1._wp,max(0._wp,x))
      w = x**vprof_b
    endif

  end function vprof_wret


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  u 3 d
  !   Purpose    :  computation of 3D wind field and advective mass fluxes
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine u3d(niter, pzsa, ptrop, ugb, vgb, ugbu, vgbv, uab, vab, t3, &
        ua, va, uter, vter, uterf, vterf, u3, v3, w3, uz500, &
        fax, fay, fac, fac_topo, fac_geo, ucor, vcor, &
        w3_geo, w3_ter, w3_ageo)

    implicit none

    integer,  intent(in   ) :: niter
    real(wp), intent(in   ) :: pzsa(:,:)
    real(wp), intent(in   ) :: ptrop(:)
    real(wp), intent(in   ) :: ugb(:,:)
    real(wp), intent(in   ) :: vgb(:,:)
    real(wp), intent(in   ) :: ugbu(:,:)   ! barotropic geostrophic zonal wind on u-points
    real(wp), intent(in   ) :: vgbv(:,:)   ! barotropic geostrophic meridional wind on v-points
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
    real(wp), intent(out  ) :: fay(:,:,:)
    real(wp), intent(out  ) :: fac(:,:)
    real(wp), intent(out  ) :: fac_topo(:,:)
    real(wp), intent(out  ) :: fac_geo(:,:)
    real(wp), intent(out  ) :: ucor(:,:)
    real(wp), intent(out  ) :: vcor(:,:)
    real(wp), intent(inout) :: w3_geo(:,:,:)
    real(wp), intent(inout) :: w3_ter(:,:,:)
    real(wp), intent(inout) :: w3_ageo(:,:,:)

    integer :: i, j, k, n, ipl, imi, jmi, kpl
    real(wp) :: pzx, pzy, dp_c, pc1, pc2, dp_l, pl1, pl2, fxpbl, fypbl, uabc, vabc, ctv
    real(wp) :: p_stl, p_top
    ! vertical shape of the compensating return flow, normalised only in the sense that its
    ! amplitude is set below by the requirement that the column mass flux vanish
    real(wp), dimension(km) :: wcom
    real(wp) :: fai, dfa, faz
    real(wp) :: u_g, v_g, ug_b, vg_b
    real(wp) :: ftrop, faz_c(3)
    real(wp) :: fcx, fcy, fmean
    real(wp), dimension(im,jm) :: psi
    ! flux potential of the topographic part, only used when i_mass_com_topo==1
    real(wp), dimension(im,jm) :: psi_topo
    ! per-component flux potentials, only used when i_mass_com_vert==1
    real(wp), dimension(im,jm) :: psi_geo, psi_ter
    ! topographic part of the barotropic geostrophic convergence, i_mass_com_vert==1
    real(wp), dimension(im,jm) :: fac_topo_geo
    real(wp), dimension(imc,jm) :: fxcol
    real(wp), dimension(im,jmc) :: fycol
    real(wp) :: fcxt, fcyt, divbar, ptrx, ptry, fmean_topo
    real(wp) :: fcx_geo, fcx_ter, fcy_geo, fcy_ter, fsum
    real(wp) :: fmean_geo, fmean_ter
    ! per-component mass fluxes; the third index selects the component
    ! (1 = barotropic geostrophic, 2 = thermal wind, 3 = ageostrophic PBL)
    real(wp), allocatable, dimension(:,:,:,:), save :: fax_c, fay_c
    real(wp) :: u3k, u3kp1, v3k, v3kp1
    real(wp) :: dplxdy, dplydxu
    real(wp) :: c_damp_eq, c_damp_pol
    ! pressure above which the geostrophic wind is disregarded, per latitude
    real(wp) :: ptdyn(jm)

    real(wp), dimension(kmc) :: uteri
    real(wp), dimension(kmc) :: vteri
    ! thermal wind of the zonal mean temperature, kept separate when i_uter_damp==2
    real(wp), dimension(kmc) :: uterz
    ! azonal thermal wind without the polar damping, for uterf/vterf
    real(wp), dimension(kmc) :: uterif, vterif
    ! azonal temperature conditioned in potential space, and its zonal mean.
    ! t3az carries the polar damping, t3azf does not, so that
    ! equatorial one, so that uterf/vterf keep their meaning of a thermal wind
    ! without the polar cap.
    real(wp), allocatable, dimension(:,:,:), save :: t3az, t3azf
    real(wp), allocatable, dimension(:,:), save :: t3zm
    real(wp), dimension(km) :: rdpl   ! 1/(pl(k)-pl(k+1)), column-invariant


    ! Allocate the per-component mass fluxes on first use. These are needed by the
    ! diagnostic w3_geo/w3_ter/w3_ageo and, with i_mass_com_vert==1, by the column
    ! mass correction itself, so they are no longer gated on l_diag_wcomp.
    if (.not.allocated(fax_c)) then
      allocate(fax_c(imc,jm,km,3))
      allocate(fay_c(im,jmc,km,3))
      fax_c = 0._wp; fay_c = 0._wp
    endif

    ! precompute reciprocal layer-pressure thickness (depends on k only, not on i,j)
    do k=1,km
      rdpl(k) = 1._wp/(pl(k)-pl(k+1))
    enddo

    !-------------------------------------------------------
    ! conditioning of the temperature that generates the thermal wind
    !-------------------------------------------------------
    ! Multiplying the thermal wind by a latitude dependent factor c(lat) leaves a
    ! field that is not the geostrophic wind of any scalar: it carries an extra
    ! divergence v*dc/dy with no physical basis, which poleward of ~55 deg exceeds
    ! the real beta term by an order of magnitude. Conditioning the temperature
    ! instead, as azslp does for the azonal sea level pressure, keeps the thermal
    ! wind exactly balanced with respect to a real temperature field, so that its
    ! divergence is the beta term and nothing else.
    !
    ! Only the azonal part can be treated this way: the full field carries a large
    ! mean and scaling that would produce a huge spurious gradient. The zonal mean
    ! part needs no such care - it has v=0 and no zonal structure, so its divergence
    ! vanishes identically and the cap on the jet can stay in wind space.
    !
    ! Replacing the polar damping by a zonal filter of the azonal temperature was
    ! tried and does not work: the planetary waves that survive any amount of
    ! filtering are enough to destabilise the model, so it is the amplitude
    ! reduction that is load bearing, not the removal of short zonal scales.
    if (i_uter_damp.eq.2) then

      if (.not.allocated(t3az)) allocate(t3az(im,jm,km), t3azf(im,jm,km), t3zm(jm,km))

      do k=1,km
        do j=1,jm
          t3zm(j,k) = 0._wp
          do i=1,im
            t3zm(j,k) = t3zm(j,k)+t3(i,j,k)
          enddo
          t3zm(j,k) = t3zm(j,k)*aim
        enddo
        ! the same polar damping as in mode 1, but applied to the potential rather
        ! than to the wind. The amplitude reduction is therefore identical, while
        ! the compensating zonal wind -Psi*dc/dy that keeps the field balanced now
        ! appears by itself when the gradients are taken.
        !
        ! The equatorial damping deliberately stays on the wind. Measured against
        ! the real beta term, the spurious divergence of a wind space factor c is
        ! dlnc/dy / (beta/f), which is 2*tan(lat)^2 for the polar factor 3cos^2 -
        ! 6 at 60 deg, 28 at 75, 114 at 82.5 - but exactly 2 for the equatorial
        ! factor 5sin^2, at every latitude. So there is little to gain at the
        ! equator, and a lot to lose: the compensating wind scales as cot(lat) and
        ! so diverges there, reaching several times the thermal wind itself inside
        ! the tropics. Enforcing exact balance where f is only kept finite by
        ! fcormin, and where geostrophy does not hold anyway, injects a large
        ! spurious zonal jet.
        !
        ! The compensating wind is -K*T_az*dc/dy: it scales with the anomaly
        ! itself, not with its gradient, so over a broad warm anomaly it is a
        ! monopole where the physical thermal wind is a dipole. Everything
        ! therefore depends on where dc/dy is put, which is what i_uter_pol sets;
        ! cdamp_pol is built once in atm_grid_init.
        do j=1,jm
          do i=1,im
            t3azf(i,j,k) = t3(i,j,k)-t3zm(j,k)
            t3az(i,j,k)  = t3azf(i,j,k) * cdamp_pol(j)
          enddo
        enddo
      enddo

    endif

    !$omp parallel do collapse(2) private(i,j,k,ipl,imi,jmi,pzx,pzy,dp_c,pc1,pc2,dp_l,pl1,pl2,fxpbl,fypbl,uabc,vabc,uteri,vteri,uterz,uterif,vterif,ctv,c_damp_eq,c_damp_pol,p_stl,p_top,wcom)
    do j=1,jm
      do i=1,im

        ipl = modulo(i,im) + 1
        imi = modulo(i - 2, im) + 1
        jmi = max(1,j-1)

        !-------------------------------------------------------
        ! Vertical profile of ageostrophic wind
        !
        ! The observed shape of the zonal mean meridional circulation, derived in
        !   benchmark/out/ps_eval/diag_shape.jl from 1324 cell columns of ERA5, five CMIP6
        !   piControl, six CMIP5/PMIP3 LGM and five CMIP5 aquaControl members, three seasons,
        !   and the Hadley, Ferrel and polar cells.  Writing
        !     zeta = (ps-p)/(ps-p_stl)         height above the ground over the depth of the
        !                                      lower branch, with p_stl = vprof_stl*ps the
        !                                      steering level where the wind changes sign
        !     xi   = (p_stl-p)/(p_stl-p_top)   the same for the return branch
        !   the mass flux carried below a level collapses onto
        !     M(p)/M_max = 1 - (1-zeta)**vprof_a       in the lower branch
        !     M(p)/M_max = 1 - xi**vprof_b             in the return branch
        !   with a scatter of 0.08 and 0.10 respectively across all of that.  So the wind
        !   decays as (1-zeta)**(vprof_a-1) from its surface value to zero at the steering
        !   level, and the return flow grows as xi**(vprof_b-1) from zero there to a maximum
        !   at the cell top.  This replaced a PBL slab of fixed sigma depth with the return
        !   flow spread uniformly below the tropopause, which fits the same data at rms 0.174
        !   against 0.023.
        !
        !   Near the ground M tends to ps-p either way, so the profile agrees with the slab
        !   there and the cross-isobar relation that produces uab/vab in slp.f90 is untouched.
        !   Unlike the slab this profile carries no tuned depth: the mass it transports follows
        !   from the steering level, M_max = (ps-p_stl)/vprof_a.
        !
        !   The cell top is a fixed function of latitude, ptopt/ptopu, and deliberately not the
        !   tropopause: across PI, LGM and aquaplanet the observed top moves by only 15-33 hPa,
        !   while p_top/p_trop moves by 26 per cent and is besides 1.42 in the Hadley cell
        !   against 0.73 in the Ferrel cell.

        ! x-component on u-points 

        ! surface pressure the profile is hung from
        if (i_pbl.eq.1) then
          pzx = 0.5_wp*(pzsa(i,j)+pzsa(imi,j))
        else if (i_pbl.eq.2) then
          pzx = 1._wp
        endif

        ! the layer mean of the wind is the mass the profile carries across the layer
        ! divided by the layer thickness
        p_stl = vprof_stl*pzx
        p_top = min(ptopt(j),0.5_wp*p_stl)
        do k=1,km
          ua(i,j,k) = uab(i,j) * (vprof_mass(pl(k+1),pzx,p_stl) &
                                 -vprof_mass(pl(k)  ,pzx,p_stl))*rdpl(k)
          wcom(k) =              (vprof_wret(pl(k+1),p_stl,p_top) &
                                 -vprof_wret(pl(k)  ,p_stl,p_top))*rdpl(k)
        enddo

        ! Compensatory velocity in the upper troposphere.
        ! The ageostrophic wind closes the mean meridional circulation, so its
        ! column-integrated mass flux has to vanish. The compensation is therefore
        ! derived from the mass flux the lower branch actually carries, using the
        ! layer thickness dplx that the mass flux itself uses. Deriving it instead
        ! from the nominal pressure depth leaves a residual wherever the profile
        ! intersects the topography, growing with surface elevation. wcom carries
        ! the shape of the return flow only; its amplitude is set here.
        fxpbl = 0._wp
        dp_c = 0._wp
        do k=1,km
          fxpbl = fxpbl + ua(i,j,k)*dplx(i,j,k)
          dp_c  = dp_c  + wcom(k)*dplx(i,j,k)
        enddo
        if (dp_c.gt.0._wp) then
          uabc = -fxpbl/dp_c
          do k=1,km
            ua(i,j,k) = ua(i,j,k)+uabc*wcom(k)
          enddo
        endif

        ! y-component on v-points 

        if (i_pbl.eq.1) then
          pzy = 0.5_wp*(pzsa(i,j)+pzsa(i,jmi))
        else if (i_pbl.eq.2) then
          pzy = 1._wp
        endif

        p_stl = vprof_stl*pzy
        p_top = min(ptopu(j),0.5_wp*p_stl)
        do k=1,km
          va(i,j,k) = vab(i,j) * (vprof_mass(pl(k+1),pzy,p_stl) &
                                 -vprof_mass(pl(k)  ,pzy,p_stl))*rdpl(k)
          wcom(k) =              (vprof_wret(pl(k+1),p_stl,p_top) &
                                 -vprof_wret(pl(k)  ,p_stl,p_top))*rdpl(k)
        enddo

        ! compensatory velocity in the upper troposphere, as for the x-component
        fypbl = 0._wp
        dp_c = 0._wp
        do k=1,km
          fypbl = fypbl + va(i,j,k)*dply(i,j,k)
          dp_c  = dp_c  + wcom(k)*dply(i,j,k)
        enddo
        if (dp_c.gt.0._wp) then
          vabc = -fypbl/dp_c
          do k=1,km
            va(i,j,k) = va(i,j,k)+vabc*wcom(k)
          enddo
        endif

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
          uterz(1) = 0._wp
          uterif(1) = 0._wp
          vterif(1) = 0._wp
          uteri(km+1) = 0._wp
          vteri(km+1) = 0._wp
          uterz(km+1) = 0._wp
          uterif(km+1) = 0._wp
          vterif(km+1) = 0._wp
          if (i_uter_damp.eq.1) then
            do k=1,km-1
              ctv = (zl(k+1)-zl(k))*g/(T0*fcort(j)) 
              uteri(k+1) = uteri(k)-ctv*(t3(i,jmi,k)-t3(i,j+1,k))/(2._wp*dy)   
              vteri(k+1) = vteri(k)+ctv*(t3(ipl,j,k)-t3(imi,j,k))/(2._wp*dxt(j)) 
            enddo
          else
            ! azonal and zonal mean parts separately, the first from the already
            ! conditioned temperature, the second untouched
            do k=1,km-1
              ctv = (zl(k+1)-zl(k))*g/(T0*fcort(j)) 
              uteri(k+1) = uteri(k)-ctv*(t3az(i,jmi,k)-t3az(i,j+1,k))/(2._wp*dy)
              vteri(k+1) = vteri(k)+ctv*(t3az(ipl,j,k)-t3az(imi,j,k))/(2._wp*dxt(j))
              uterif(k+1) = uterif(k)-ctv*(t3azf(i,jmi,k)-t3azf(i,j+1,k))/(2._wp*dy)
              vterif(k+1) = vterif(k)+ctv*(t3azf(ipl,j,k)-t3azf(imi,j,k))/(2._wp*dxt(j))
              uterz(k+1) = uterz(k)-ctv*(t3zm(jmi,k)-t3zm(j+1,k))/(2._wp*dy)
            enddo
          endif

          ! thermal wind in layers, dampened at equator and poles
          c_damp_pol = cdamp_pol(j)
          c_damp_eq  = min(1._wp, c_uter_eq*sint(j)**2)
          if (i_uter_damp.eq.1) then
            do k=1,km
              uter(i,j,k) = 0.5_wp*(uteri(k)+uteri(k+1)) * c_damp_eq * c_damp_pol
              vter(i,j,k) = 0.5_wp*(vteri(k)+vteri(k+1)) * c_damp_eq * c_damp_pol
              ! thermal wind without polar damping, for EKE production
              uterf(i,j,k) = 0.5_wp*(uteri(k)+uteri(k+1)) * c_damp_eq
              vterf(i,j,k) = 0.5_wp*(vteri(k)+vteri(k+1)) * c_damp_eq
            enddo
          else
            do k=1,km
              uter(i,j,k) = (0.5_wp*(uteri(k)+uteri(k+1)) &
                          + 0.5_wp*(uterz(k)+uterz(k+1)) * c_damp_pol) * c_damp_eq
              vter(i,j,k) = 0.5_wp*(vteri(k)+vteri(k+1)) * c_damp_eq
              ! as above but without the polar cap anywhere, for EKE production
              uterf(i,j,k) = (0.5_wp*(uterif(k)+uterif(k+1)) &
                           + 0.5_wp*(uterz(k)+uterz(k+1))) * c_damp_eq
              vterf(i,j,k) = 0.5_wp*(vterif(k)+vterif(k+1)) * c_damp_eq
            enddo
          endif

        endif

      enddo
    enddo
    !$omp end parallel do


    !-------------------------------------------------------
    ! advective mass transport
    !-------------------------------------------------------

    ! Top of the layer that carries the geostrophic wind. A fixed pressure cuts
    ! the wind off well below the tropopause in the tropics and well above it at
    ! high latitudes, so it can optionally follow the tropopause instead.
    do j=1,jm
      if (i_ptopdyn.eq.1) then
        ptdyn(j) = ptopdyn
      else
        ptdyn(j) = ptrop(j)
      endif
    enddo

    !$omp parallel do collapse(2) private(i,j,k,imi,u_g,v_g,dplxdy,dplydxu,ftrop)
    do k=1,km
      do j=1,jm

        ! fraction of the layer that lies below the dynamical top; the geostrophic
        ! wind is limited to the troposphere, see below
        if (pl(k+1).ge.ptdyn(j)) then
          ftrop = 1._wp
        elseif (pl(k+1).lt.ptdyn(j) .and. pl(k).ge.ptdyn(j)) then
          ftrop = (pl(k)-ptdyn(j))*rdpl(k)
        else
          ftrop = 0._wp
        endif

        ! x-components

        do i=1,im
          imi = modulo(i - 2, im) + 1
          ! geostrophic zonal wind on u-points, limited to troposphere.
          ! ugbu is built on the faces in u2d.f90; with i_ugb_psi=0 it is exactly the
          ! average 0.5*(ugb(imi,j)+ugb(i,j)) this line used to form here, and with
          ! i_ugb_psi=1 it is the curl of the corner streamfunction, which must NOT be
          ! routed through T-points or the non-divergence is lost.
          if (pl(k+1).ge.ptdyn(j)) then
            u_g  = ugbu(i,j) + 0.5_wp*(uter(imi,j,k)+uter(i,j,k))
          elseif (pl(k+1).lt.ptdyn(j).and.pl(k).ge.ptdyn(j)) then
            u_g  = (ugbu(i,j) + 0.5_wp*(uter(imi,j,k)+uter(i,j,k))) *(pl(k)-ptdyn(j))*rdpl(k)
          else
            u_g  = 0._wp
          endif 
          dplxdy = dplx(i,j,k)*dy
          ! mass flux
          fax(i,j,k) = (u_g+ua(i,j,k))*dplxdy ! m/s * kg/m2 * m = kg/s
          ! the same mass flux, split by wind component; the three add up to fax
          ! exactly, since ftrop reproduces the branch taken above
          fax_c(i,j,k,1) = ugbu(i,j)*ftrop * dplxdy
          fax_c(i,j,k,2) = 0.5_wp*(uter(imi,j,k)+uter(i,j,k))*ftrop * dplxdy
          fax_c(i,j,k,3) = ua(i,j,k) * dplxdy
        enddo
        ! periodic boundary conditions
        fax(imc,j,k) = fax(1,j,k)
        fax_c(imc,j,k,:) = fax_c(1,j,k,:)

        ! y-components

        do i=1,im
          if (j.eq.1) then
            ! N and S boundary conditions, no flux
            fay(i,1,k)  = 0._wp
            fay(i,jmc,k)  = 0._wp
            fay_c(i,1,k,:)   = 0._wp
            fay_c(i,jmc,k,:) = 0._wp
          else
            ! geostrophic meridional wind on v-points, limited to troposphere; see the
            ! note on ugbu above
            if (pl(k+1).ge.ptdyn(j)) then
              v_g  = vgbv(i,j) + 0.5_wp*(vter(i,j-1,k)+vter(i,j,k))
            elseif (pl(k+1).lt.ptdyn(j).and.pl(k).ge.ptdyn(j)) then
              v_g  = (vgbv(i,j) + 0.5_wp*(vter(i,j-1,k)+vter(i,j,k))) *(pl(k)-ptdyn(j))*rdpl(k)
            else
              v_g  = 0._wp
            endif 
            dplydxu = dply(i,j,k)*dxu(j)
            ! mass flux
            fay(i,j,k) = (v_g+va(i,j,k))*dplydxu  ! m/s * kg/m2 * m = kg/s
            ! the same mass flux, split by wind component
            fay_c(i,j,k,1) = vgbv(i,j)*ftrop * dplydxu
            fay_c(i,j,k,2) = 0.5_wp*(vter(i,j-1,k)+vter(i,j,k))*ftrop * dplydxu
            fay_c(i,j,k,3) = va(i,j,k) * dplydxu
          endif
        enddo

      enddo
    enddo 
    !$omp end parallel do


    if (l_output_flx3d .and. niter.eq.1) then
      ! column convergence before compensation
      do j=1,jm
        do i=1,im
          fac(i,j) = 0._wp
          do k=1,km 
            fac(i,j) = fac(i,j) + fax(i,j,k) -fax(i+1,j,k) +fay(i,j+1,k) -fay(i,j,k)
          enddo
        enddo
      enddo
    endif

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
    ! The residual is removed with a corrective mass flux F_c = grad(psi) obtained
    ! from  laplacian(psi) = -(column convergence),  discretised so that its
    ! convergence matches the stencil used above. This is the smallest correction
    ! that closes the budget and, unlike adding the whole residual to one face, it
    ! has no preferred direction. Enforcing it per column also makes the zonal mean
    ! meridional mass flux vanish at every latitude, so no separate step is needed:
    ! summing the column constraint zonally telescopes the zonal fluxes away and,
    ! with no flux through the poles, leaves sum_i fay(i,j) = 0 for all j.
    !
    ! With i_mass_com_topo==1 the topographic term is separated from the rest and
    ! given its own potential, because the two belong at different heights. The
    ! mass that upslope flow piles into a column is in reality lifted over the
    ! barrier and leaves just above it, whereas the divergence of the balanced
    ! flow is compensated by the secondary circulation high up. The level matters:
    ! the return flow removes the moisture that the upslope flow just brought in
    ! if it is placed near the surface, and removes far too much dry static energy
    ! (which grows with height) if it is placed at the tropopause.
    !
    ! i_mass_com_vert==1 carries the same argument through to what is left. The
    ! remaining residual is not one thing either: the barotropic geostrophic wind is
    ! height independent through the troposphere, so its divergence arrives spread by
    ! layer mass, while the thermal wind is a shear anchored at zero and makes its
    ! divergence aloft. Measured over 55-140E / 50-70N in JJA, the centroid of the
    ! divergence generation is 3.3 km for the barotropic part and 5.3-6.6 km for the
    ! thermal wind, against a dp_com band at 6.5-12 km. So each component gets its own
    ! potential and is returned through the same layers, in the same proportion, that
    ! its own flux occupies at that face - which makes the correction carry exactly the
    ! dry static energy and moisture per unit mass that the erroneous flux carried, and
    ! therefore neutral in both. For the barotropic component |fax_c(..,1)| is
    ! proportional to ftrop*dplx, so that rule reduces to layer mass weighting.
    !
    ! The ageostrophic wind needs no potential: it closes column-wise by construction,
    ! so its column flux vanishes at every face and fac = fac_geo + fac_ter exactly.

    ! column convergence, kg/s (also reported as the diagnostic fac)
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
      endif

      call solve_poisson(-(fac_topo-fmean_topo), psi_topo)

      if (i_mass_com_vert.eq.0) then
        call solve_poisson(-(fac-fac_topo+fmean_topo-fmean), psi)
      endif

    else

      fac_topo(:,:) = 0._wp
      fac_topo_geo(:,:) = 0._wp
      psi_topo(:,:) = 0._wp
      if (i_mass_com_vert.eq.0) then
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

    ! zonal component of the corrective flux, on u-points
    !$omp parallel do collapse(2) private(i,j,k,n,imi,fcx,fcxt,fcx_geo,fcx_ter,fsum,dp_c,pc1,pc2,pl1,pl2,dp_l,ptrx)
    do j=1,jm
      do i=1,im
        imi = modulo(i - 2, im) + 1
        fcx = dy/dxt(j)*(psi(imi,j)-psi(i,j))
        fcxt = 0._wp
        if (i_mass_com_topo.eq.1) fcxt = dy/dxt(j)*(psi_topo(imi,j)-psi_topo(i,j))
        fcx_geo = 0._wp
        fcx_ter = 0._wp
        if (i_mass_com_vert.eq.1) then
          fcx_geo = dy/dxt(j)*(psi_geo(imi,j)-psi_geo(i,j))
          fcx_ter = dy/dxt(j)*(psi_ter(imi,j)-psi_ter(i,j))
        endif
        ucor(i,j) = 0._wp
        if (plx(i,j).gt.0._wp) ucor(i,j) = (fcx+fcxt+fcx_geo+fcx_ter)/(plx(i,j)*dy)
        if (i_mass_com_topo.eq.1) then
          ! topographic part: terrain following band, from dps_com_topo above the
          ! local surface up to the tropopause, weighted by the layer mass so that
          ! most of it sits just above the barrier
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
              fax(i,j,k) = fax(i,j,k) + fcxt*dp_l*rdpl(k)*dplx(i,j,k)/dp_c
            enddo
          elseif (plx(i,j).gt.0._wp) then
            ! degenerate column, fall back to the column mean
            do k=1,km
              fax(i,j,k) = fax(i,j,k) + fcxt*dplx(i,j,k)/plx(i,j)
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
                fax(i,j,k) = fax(i,j,k) + fcx*abs(fax_c(i,j,k,n))/fsum
              enddo
            elseif (plx(i,j).gt.0._wp) then
              ! no wind of this component in the column, fall back to the column mean
              do k=1,km
                fax(i,j,k) = fax(i,j,k) + fcx*dplx(i,j,k)/plx(i,j)
              enddo
            endif
          enddo
        else
          ! confine the correction to a slab of thickness dp_com hanging below the
          ! tropopause, the same construction the ageostrophic closure uses above. A
          ! fixed pressure pair instead collapsed the band to a 10 hPa sliver wherever
          ! the tropopause fell below it, which is everywhere poleward of ~45 deg.
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
              fax(i,j,k) = fax(i,j,k) + fcx*dp_l*rdpl(k)*dplx(i,j,k)/dp_c
            enddo
          endif
        endif
      enddo
      ! periodic boundary condition
      do k=1,km
        fax(imc,j,k) = fax(1,j,k)
      enddo
    enddo
    !$omp end parallel do

    ! meridional component, on v-points; no flux through the poles
    vcor(:,1) = 0._wp
    !$omp parallel do collapse(2) private(i,j,k,n,fcy,fcyt,fcy_geo,fcy_ter,fsum,dp_c,pc1,pc2,pl1,pl2,dp_l,ptry)
    do j=2,jm
      do i=1,im
        fcy = dxu(j)/dy*(psi(i,j)-psi(i,j-1))
        fcyt = 0._wp
        if (i_mass_com_topo.eq.1) fcyt = dxu(j)/dy*(psi_topo(i,j)-psi_topo(i,j-1))
        fcy_geo = 0._wp
        fcy_ter = 0._wp
        if (i_mass_com_vert.eq.1) then
          fcy_geo = dxu(j)/dy*(psi_geo(i,j)-psi_geo(i,j-1))
          fcy_ter = dxu(j)/dy*(psi_ter(i,j)-psi_ter(i,j-1))
        endif
        vcor(i,j) = 0._wp
        if (ply(i,j).gt.0._wp .and. dxu(j).gt.0._wp) &
          vcor(i,j) = (fcy+fcyt+fcy_geo+fcy_ter)/(ply(i,j)*dxu(j))
        if (i_mass_com_topo.eq.1) then
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
              fay(i,j,k) = fay(i,j,k) + fcyt*dp_l*rdpl(k)*dply(i,j,k)/dp_c
            enddo
          elseif (ply(i,j).gt.0._wp) then
            do k=1,km
              fay(i,j,k) = fay(i,j,k) + fcyt*dply(i,j,k)/ply(i,j)
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
                fay(i,j,k) = fay(i,j,k) + fcy*abs(fay_c(i,j,k,n))/fsum
              enddo
            elseif (ply(i,j).gt.0._wp) then
              do k=1,km
                fay(i,j,k) = fay(i,j,k) + fcy*dply(i,j,k)/ply(i,j)
              enddo
            endif
          enddo
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
              fay(i,j,k) = fay(i,j,k) + fcy*dp_l*rdpl(k)*dply(i,j,k)/dp_c
            enddo
          endif
        endif
      enddo
    enddo
    !$omp end parallel do


    if (niter.eq.1) then

      ! total 3-D wind on T-points, interpolate to levels

      !$omp parallel do collapse(2) private(i,j,ipl,k,kpl,ug_b,vg_b,u3k,v3k,u3kp1,v3kp1,faz,faz_c)
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

          ! diagnostic: the same vertical integral, one wind component at a time.
          ! The three add up to w3 (up to the mass-restoring applied elsewhere), so
          ! their relative size shows which component supplies the divergence.
          if (l_diag_wcomp) then
            faz_c(:) = 0._wp
            w3_geo(i,j,1)  = 0._wp
            w3_ter(i,j,1)  = 0._wp
            w3_ageo(i,j,1) = 0._wp
            do k=1,km
              kpl = min(k+1,km)
              faz_c(:) = faz_c(:) + fax_c(i,j,k,:)-fax_c(i+1,j,k,:) + fay_c(i,j+1,k,:)-fay_c(i,j,k,:)
              w3_geo(i,j,k+1)  = faz_c(1)/(sqr(i,j)*ra*pl(kpl))
              w3_ter(i,j,k+1)  = faz_c(2)/(sqr(i,j)*ra*pl(kpl))
              w3_ageo(i,j,k+1) = faz_c(3)/(sqr(i,j)*ra*pl(kpl))
            enddo
          endif

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
