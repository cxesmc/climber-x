!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s l p _ m o d
!
!  Purpose : sea level pressure
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
module slp_mod

  use atm_params, only : wp, dp
  use constants, only : pi, r_earth, omega, g, Rd, T0, Le
  use atm_params, only : ra, p0, hatm
  use atm_params, only : c_slp_1, c_slp_2, c_slp_4, c_slp_5
  use atm_params, only : c_slp_trop_f, c_slp_trop_tau, c_slp_trop_lat0, c_slp_trop_lat1
  use atm_params, only : l_aslp_temp_adv, c_aslp_temp_tau
  use atm_params, only : l_aslp_topo, c_aslp_sw_amp, c_aslp_sw_tau, c_aslp_sw_h, c_aslp_sw_umin, c_aslp_sw_lateq, n_aslp_sw
  use atm_params, only : c_aslp_sw_fr, c_aslp_sw_nbv
  use atm_params, only : c_aslp_sw_eddy, i_aslp_sw_eddy, n_aslp_sw_eddy
  use atm_params, only : l_aslp_gill, c_aslp_gill, c_aslp_gill_tau, c_aslp_gill_c, c_aslp_gill_lat0, c_aslp_gill_lat1, n_aslp_gill
  use atm_params, only : i_aslp_gill_lat, c_aslp_gill_trel
  use timer, only : dt_atm, sec_day
  use atm_params, only : cp
  use atm_params, only : c_mmc_had, c_mmc_fer, c_mmc_pol, c_mmc_1, c_mmc_2, i_mmc_geo, c_mmc_edge
  use atm_params, only : i_mmc_itcz, c_mmc_efe
  use atm_params, only : i_mmc_fer
  use atm_params, only : i_mmc_had, c_mmc_q, q_mmc_ref
  use atm_params, only : c_mmc_dt0, c_mmc_dt1, hdv_mmc_ref, c_mmc_mt, amas
  use atm_params, only : nsmooth_aslp, aslp_polar_lat
  use atm_grid, only : im, jm, jmc, aim, jeq, jts, jtn, jps, jpn, dy, pl, k500, k850, k300, i_ice
  use atm_grid, only : fcorua, sint, cost, fiu, fit, dxt, pblu
  use atm_grid, only : zc, k_dse_lo, k_dse_up
  use smooth_atm_mod, only : smooth2_m, zona
  !$ use omp_lib

  implicit none

  private
  public :: zslp, azslp, efe_diag, col_prev_unset

  real(wp), parameter :: col_prev_unset = -9.e9_wp   !! sentinel of tcol_prev/wcon_prev: no previous step (efe_diag)

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  azslp
  !   Purpose    :  compute azonal component of sea level pressure
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine azslp(tsksl, htrop, zsa, uz500, u3f, sam, eke, prc, ficz, &
      aslp, &
      aslp_temp, aslp_topo, aslp_eddy, aslp_gill, aslp_temp_slow, dz500, atsl)

    use, intrinsic :: iso_c_binding 

    implicit none

    ! include FFTW3 library for Fourier Transform (https://www.fftw.org/)
    include 'fftw3.f03'

    real(wp), intent(in   ) :: tsksl(:,:)
    real(wp), intent(in   ) :: htrop(:,:)
    real(wp), intent(in   ) :: zsa(:,:)
    real(wp), intent(in   ) :: uz500(:)       !! zonal mean 500 hPa zonal wind for the upstream displacement of aslp_temp, thermal wind damped with c_uz500_pol
    real(wp), intent(in   ) :: u3f(:,:,:)     !! 3D zonal wind with undamped thermal wind, previous step, used if l_aslp_topo
    real(wp), intent(in   ) :: sam(:,:)       !! synoptic eddy kinetic energy, used only if l_aslp_topo and c_aslp_sw_eddy>0
    real(wp), intent(in   ) :: eke(:,:)       !! ERA-Interim 2-6 day eddy kinetic energy climatology, used only if i_aslp_sw_eddy==2
    real(wp), intent(in   ) :: prc(:,:)       !! precipitation (kg/m2/s), previous step, used only if l_aslp_gill
    real(wp), intent(in   ) :: ficz        !! ITCZ latitude of the meridional cell geometry from the previous step (radians)

    real(wp), intent(inout) :: aslp(:,:)
    real(wp), intent(inout) :: aslp_temp(:,:)
    real(wp), intent(inout) :: aslp_topo(:,:)
    real(wp), intent(inout) :: aslp_eddy(:,:)
    real(wp), intent(inout) :: aslp_gill(:,:)
    real(wp), intent(inout) :: aslp_temp_slow(:,:)   !! slowly relaxed aslp_temp, used only if c_slp_trop_f<1

    real(wp), intent(out  ) :: atsl(:,:)
    real(wp), intent(out  ) :: dz500(:,:)

    integer :: i, j
    real(wp) :: k, uz
    real(wp) :: tslz(jm)
    real(wp) :: htropz(jm)
    real(wp) :: uz500s(jm)
    type(C_PTR) :: plan_r2c_temp, plan_c2r_temp
    real(wp) :: aslp_frc(im,jm)
    real(wp) :: aslp_gill_new(im,jm), cgill
    real(wp) :: cslow, ftrop, wtrop, alat
    real(dp), dimension(im) :: aslp_temp_dp
    complex(dp), dimension(im/2+1) :: aslp_temp_fft


    ! zonal mean zonal wind for the upstream displacement of aslp_temp:
    ! smooth zonal mean 500 hPa zonal wind
    do j=2,jm-1
      uz500s(j) = 0.4_wp*uz500(j)+0.3_wp*(uz500(j+1)+uz500(j-1))
    enddo
    uz500s(1)  = 0.5_wp*(uz500(1)+uz500(2))
    uz500s(jm) = 0.5_wp*(uz500(jm)+uz500(jm-1))

    ! zonal mean sea level temperature and tropopause height
    do j=1,jm
      tslz(j) = 0._wp      
      htropz(j) = 0._wp      
      do i=1,im
        tslz(j) = tslz(j) + tsksl(i,j)*aim 
        htropz(j) = htropz(j) + htrop(i,j)*aim 
      enddo
    enddo

    !------------------------------------------------
    ! temperature related azonal sea level pressure 
    !------------------------------------------------

    !$omp parallel do private(i,j)
    do j=1,jm

      do i=1,im
        ! azonal sea level temperature
        atsl(i,j) = tsksl(i,j)-tslz(j)
      enddo

      do i=1,im
        if (j.eq.1 .or. j.eq.jm) then
          ! azonal SLP vanishes at the Poles
          aslp_temp(i,j) = 0._wp
        else
          ! as in CLIMBER-2, Petoukhov 2000, eq. (17)
          aslp_temp(i,j) = -c_slp_1*g*p0*10000._wp/(2._wp*Rd*T0**2)*atsl(i,j)
          !aslp_temp(i,j) = -c_slp_1*g*p0*htropz(j)/(2._wp*Rd*T0**2)*atsl(i,j)
        endif
      enddo

      ! ensure zero zonal mean
      aslp_temp(:,j) = aslp_temp(:,j) - sum(aslp_temp(:,j))*aim

    enddo
    !$omp end parallel do


    !------------------------------------------------
    ! upstream displacement of temperature related azonal sea level pressure
    !------------------------------------------------

    ! The hydrostatic surface pressure anomaly is set by the depth-integrated cooling of
    ! the air column, not by the local sea level temperature. Air crossing a cold continent
    ! in the mean westerly flow keeps cooling as it goes, so the temperature minimum
    ! accumulates at the downstream edge while the cooling that drives the surface high is
    ! centred upstream of it. atsl therefore lags the pressure anomaly: in ERA-Interim the
    ! DJF azonal sea level pressure correlates best with the azonal sea level temperature
    ! displaced 15-20 deg to the east (r=-0.83 over 40-75N, r=-0.93 over North America,
    ! against only -0.70 and -0.54 at zero lag).
    ! The equilibrium thermal pressure anomaly is therefore displaced upstream, as the
    ! steady state of
    !   -uz*d(aslp_temp)/dx = (aslp_temp_eq-aslp_temp)/c_aslp_temp_tau ,
    ! solved in Fourier space, where each zonal wavenumber n is shifted upstream by the
    ! phase angle atan(k*n*uz*c_aslp_temp_tau) and damped by 1/sqrt(1+(k*n*uz*tau)**2).

    if (l_aslp_temp_adv) then

      ! Make forward and backward plans for the FFT
      plan_r2c_temp = fftw_plan_dft_r2c_1d(im, aslp_temp_dp, aslp_temp_fft, FFTW_ESTIMATE)
      plan_c2r_temp = fftw_plan_dft_c2r_1d(im, aslp_temp_fft, aslp_temp_dp, FFTW_ESTIMATE)

      ! azonal SLP vanishes at the Poles, nothing to displace there
      do j=2,jm-1

        k = 2._wp*pi/(2._wp*pi*r_earth*cost(j))      ! lowest zonal wavenumber
        uz = max(0.1_wp,uz500s(j))                   ! advecting wind, westerly only

        aslp_temp_dp(:) = aslp_temp(:,j)

        !  forward transform the data
        call fftw_execute_dft_r2c(plan_r2c_temp, aslp_temp_dp, aslp_temp_fft)

        ! negative imaginary part = displacement upstream (to the west) for westerly uz
        do i=1,im/2+1
          aslp_temp_fft(i) = aslp_temp_fft(i) &
            / cmplx(1._dp, -real(k*(i-1)*uz*c_aslp_temp_tau,dp), dp)
        enddo

        ! backward transform the data
        call fftw_execute_dft_c2r(plan_c2r_temp, aslp_temp_fft, aslp_temp_dp)

        aslp_temp(:,j) = aslp_temp_dp(:) * aim   ! aim accounts for the FFTW normalisation

      enddo

      ! destroy FFT plans
      call fftw_destroy_plan(plan_r2c_temp)
      call fftw_destroy_plan(plan_c2r_temp)

    endif

    !------------------------------------------------
    ! topographic stationary planetary waves
    !------------------------------------------------

    if (l_aslp_topo) then

      ! linear barotropic stationary waves, coupled in latitude
      if (i_aslp_sw_eddy.eq.2) then
        call stationary_waves(zsa, u3f, eke, aslp_topo, aslp_eddy)
      else
        call stationary_waves(zsa, u3f, sam, aslp_topo, aslp_eddy)
      endif

    else

      aslp_topo(:,:) = 0._wp
      aslp_eddy(:,:) = 0._wp

    endif


    !------------------------------------------------
    ! baroclinic (Gill) response to the azonal latent heating
    !------------------------------------------------

    if (l_aslp_gill) then
      if (c_aslp_gill_trel.gt.0._wp) then
        ! The Gill response relaxes towards the response to the current heating with its own time scale.
        ! Heating -> Kelvin low east of it -> convergence -> heating further east closes an eastward
        ! propagating loop (10-40 d tropical waves); a relaxation time longer than the loop period breaks it
        ! while keeping the seasonal monsoon/ITCZ response, which only lags by atan(2pi*trel/1yr).
        call gill_response(prc, ficz, aslp_gill_new)
        cgill = min(1._wp, dt_atm/(c_aslp_gill_trel*sec_day))
        aslp_gill = aslp_gill + cgill*(aslp_gill_new-aslp_gill)
      else
        call gill_response(prc, ficz, aslp_gill)
      endif
    else
      aslp_gill(:,:) = 0._wp
    endif


    !------------------------------------------------
    ! damping of the fast part of the thermal azonal sea level pressure
    !------------------------------------------------

    ! Over land the thermal low closes a delayed negative feedback (warm skin -> low ->
    ! ageostrophic convergence -> ascent cooling) that oscillates with a ~10 d period in JAS over
    ! West Africa, and with 20-40 d over Tibet/Central Asia in MJJASO (both grow with CO2).
    ! Lowering c_slp_2 kills the tropical one only by halving the equilibrium tropical aslp (the
    ! equatorial damping accumulates as c S/(1-(1-c) S)), which changes the mean climate. Here only
    ! the departure of aslp_temp from its slowly relaxed state is damped, by f, with a relaxation
    ! time tau: the time mean is unchanged and the seasonal cycle lags by ~(1-f)*tau.
    ! f is c_slp_trop_f equatorward of c_slp_trop_lat0 and 1 (no damping) poleward of c_slp_trop_lat1,
    ! with a cosine transition in between. Extratropical damping (tested down to tau 5 d) creates a
    ! new ~10 d mode over Tibet at the LGM.
    if (c_slp_trop_f.lt.1._wp) then
      if (any(aslp_temp_slow.gt.1.e19_wp)) aslp_temp_slow = aslp_temp
      cslow = min(1._wp, dt_atm/(c_slp_trop_tau*sec_day))
      do j=1,jm
        alat = abs(asin(sint(j)))*180._wp/pi
        if (alat.le.c_slp_trop_lat0) then
          wtrop = 1._wp
        else if (alat.ge.c_slp_trop_lat1) then
          wtrop = 0._wp
        else
          wtrop = 0.5_wp*(1._wp+cos(pi*(alat-c_slp_trop_lat0)/(c_slp_trop_lat1-c_slp_trop_lat0)))
        endif
        ftrop = wtrop*c_slp_trop_f + (1._wp-wtrop)
        do i=1,im
          aslp_temp_slow(i,j) = aslp_temp_slow(i,j) + cslow*(aslp_temp(i,j)-aslp_temp_slow(i,j))
          aslp_temp(i,j) = aslp_temp_slow(i,j) + ftrop*(aslp_temp(i,j)-aslp_temp_slow(i,j))
        enddo
      enddo
    endif


    !-------------------------------------------------------------------------------
    ! geopotential height zonal anomalies at 500 hPa from planetary stationary waves
    !-------------------------------------------------------------------------------

    ! height anomaly of the equivalent barotropic (~300 hPa) level, scaled as the sea level pressure
    dz500(:,:) = aslp_topo(:,:)/(ra*g)  ! m


    !------------------------------------------------
    ! azonal component of sea level pressure 
    !------------------------------------------------

    ! the forcing is tapered once per call through its vorticity (does not accumulate), while smoothing and
    ! equatorial damping act on the time-relaxed state at every call, so they accumulate: the equilibrium
    ! is c S X/(1-(1-c) S) instead of S X (analysis/stationary_wave/aslp_smooth_loop.jl)
    do j=1,jm
      do i=1,im
        aslp_frc(i,j) = aslp_temp(i,j)+aslp_topo(i,j)+aslp_gill(i,j)
      enddo
    enddo
    call vorticity_taper(aslp_frc)
    do j=1,jm
      do i=1,im
        ! relax in time
        aslp(i,j) = (1._wp-c_slp_2)*aslp(i,j)+c_slp_2*aslp_frc(i,j)
      enddo
    enddo

    ! smooth in space
    call smooth2_m(aslp,nsmooth_aslp)

    ! equatorial damping
    do j=1,jm
      do i=1,im
        aslp(i,j) = aslp(i,j) * min(1._wp,c_slp_5+c_slp_4*(sint(j)**2-sint(jm/2)**2))
      enddo
    enddo


    return

  end subroutine azslp


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  v o r t i c i t y _ t a p e r
  !   Purpose    :  polar taper of an azonal field through its vorticity
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Multiplying the field by a latitude taper T adds the curvature 2 grad(T).grad(p) + p lap(T), which
  ! for a deep subpolar low is a cyclonic band equatorward and an anticyclonic band poleward of the
  ! taper edge (analysis/amoc_stab/polar_taper_offline.jl). Instead, poleward of aslp_polar_lat the
  ! Laplacian of the field is multiplied by T (raised cosine, 1 at aslp_polar_lat, 0 at the pole) and
  ! inverted per zonal wavenumber n,
  !   lap_n(p_n) = T * lap_n(x_n) ,   p_n = x_n at the last row equatorward of aslp_polar_lat ,
  ! with the same discrete Laplacian (5 point, face cosines, zero flux through the pole face). The
  ! vorticity then decays monotonically to the pole and never changes sign, and each wavenumber decays
  ! harmonically (regularly) towards the pole. The zonal mean (n=0) is left unchanged.
  subroutine vorticity_taper(x)

    use, intrinsic :: iso_c_binding

    implicit none

    include 'fftw3.f03'

    real(wp), intent(inout) :: x(:,:)

    integer :: j, k, n, m, jb, hemi, nrow
    integer :: rows(jm)
    real(wp) :: dphi, dlon, al, tap, den
    real(wp) :: a(jm), b(jm), c(jm)
    complex(dp) :: rhs(jm), cp(jm), dpv(jm), sol(jm), pb
    type(C_PTR) :: plan_r2c, plan_c2r
    real(dp), dimension(im) :: row
    complex(dp), dimension(im/2+1) :: rowc
    complex(dp) :: xh(im/2+1,jm)


    dphi = pi/real(jm,wp)
    dlon = 2._wp*pi/real(im,wp)

    plan_r2c = fftw_plan_dft_r2c_1d(im, row, rowc, FFTW_ESTIMATE)
    plan_c2r = fftw_plan_dft_c2r_1d(im, rowc, row, FFTW_ESTIMATE)

    do j=1,jm
      row(:) = x(:,j)
      call fftw_execute_dft_r2c(plan_r2c, row, rowc)
      xh(:,j) = rowc(:)
    enddo

    do hemi=1,2

      ! rows of the cap ordered from the edge to the pole (internal j=1 is the North Pole row);
      ! jb is the Dirichlet row, the last one equatorward of aslp_polar_lat
      nrow = 0
      if (hemi.eq.1) then
        do j=jm/2,1,-1
          if (fit(j)*180._wp/pi.gt.aslp_polar_lat) then
            nrow = nrow+1
            rows(nrow) = j
          endif
        enddo
        if (nrow.eq.0) cycle
        jb = rows(1)+1
      else
        do j=jm/2+1,jm
          if (-fit(j)*180._wp/pi.gt.aslp_polar_lat) then
            nrow = nrow+1
            rows(nrow) = j
          endif
        enddo
        if (nrow.eq.0) cycle
        jb = rows(1)-1
      endif
      m = nrow

      do n=1,im/2
        do k=1,m
          al = abs(fit(rows(k)))
          a(k) = cos(al-0.5_wp*dphi)/(r_earth**2*cos(al)*dphi**2)
          c(k) = max(0._wp,cos(al+0.5_wp*dphi))/(r_earth**2*cos(al)*dphi**2)
          if (k.eq.m) c(k) = 0._wp   ! no flux through the pole face
          b(k) = -(a(k)+c(k)) - (2._wp-2._wp*cos(real(n,wp)*dlon))/(r_earth*cos(al)*dlon)**2
        enddo
        ! tapered Laplacian of the field
        pb = xh(n+1,jb)
        do k=1,m
          al = abs(fit(rows(k)))*180._wp/pi
          tap = 0.5_wp*(1._wp+cos(pi*(al-aslp_polar_lat)/(90._wp-aslp_polar_lat)))
          if (k.eq.1) then
            rhs(k) = a(k)*pb
          else
            rhs(k) = a(k)*xh(n+1,rows(k-1))
          endif
          rhs(k) = rhs(k) + b(k)*xh(n+1,rows(k))
          if (k.lt.m) rhs(k) = rhs(k) + c(k)*xh(n+1,rows(k+1))
          rhs(k) = tap*rhs(k)
        enddo
        rhs(1) = rhs(1) - a(1)*pb
        ! Thomas algorithm
        cp(1) = c(1)/b(1)
        dpv(1) = rhs(1)/b(1)
        do k=2,m
          den = b(k) - a(k)*real(cp(k-1),wp)
          cp(k) = c(k)/den
          dpv(k) = (rhs(k) - a(k)*dpv(k-1))/den
        enddo
        sol(m) = dpv(m)
        do k=m-1,1,-1
          sol(k) = dpv(k) - cp(k)*sol(k+1)
        enddo
        do k=1,m
          xh(n+1,rows(k)) = sol(k)
        enddo
      enddo

      ! back to grid space for the cap rows; the zonal mean (n=0) is kept
      do k=1,m
        rowc(:) = xh(:,rows(k))
        call fftw_execute_dft_c2r(plan_c2r, rowc, row)
        x(:,rows(k)) = row(:)*aim
      enddo

    enddo

    call fftw_destroy_plan(plan_r2c)
    call fftw_destroy_plan(plan_c2r)

    return

  end subroutine vorticity_taper


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  stationary_waves
  !   Purpose    :  sea level pressure of the topographically forced stationary
  !                 planetary waves from a linear barotropic model
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine stationary_waves(zsa, u3f, sam, aslp_sw, aslp_eddy)

    ! Steady linear barotropic vorticity equation at the equivalent barotropic level (~300 hPa),
    ! linearised about the zonal mean wind U(phi) of that level,
    !   U/(a cos) d(zeta')/dlambda + q_y v' = -r zeta' + F ,
    !   q_y = 2 Omega cos/a - (1/a) d/dphi [ (1/(a cos)) d(U cos)/dphi ] ,
    ! forced by the low-level (~850 hPa) zonal mean flow over the topography,
    !   F = -(U_850 f/H) (1/(a cos)) d(zsa)/dlambda .
    ! In Fourier space each zonal wavenumber n is a second order ODE in latitude,
    !   (i n U/(a cos) + r) lap_n(psi_n) + i n q_y/(a cos) psi_n = F_n ,
    ! solved as one tridiagonal system with psi = 0 at the poles. Unlike a row-by-row
    ! Charney-Eliassen solution the meridional structure is not prescribed, so
    ! the waves can propagate across latitudes (Rockies -> Labrador trough -> European ridge)
    ! and the solution is not tied to the K^2 = beta/U resonance of a fixed meridional wavenumber.
    ! The height anomaly f*psi/g is converted to sea level pressure assuming an equivalent
    ! barotropic structure, aslp_sw = c_aslp_sw_amp * ra*g * f*psi/g.
    !
    ! Optional transient-eddy forcing (c_aslp_sw_eddy>0). For non-divergent eddies the vorticity
    ! forcing depends only on the eddy momentum flux u'v' and the anisotropy M = (u'^2-v'^2)/2,
    !   F_e = (d_yy - d_xx) u'v' + 2 d_xy M ,
    ! the eddy kinetic energy itself drops out. The azonal parts of M and u'v' are closed with the
    ! synoptic kinetic energy E = sam (i_aslp_sw_eddy=1) or the ERA-Interim 2-6 day eke_2_6_int
    ! climatology (i_aslp_sw_eddy=2) and its gradients,
    !   (M, u'v') = c_aslp_sw_eddy * [ c1 E' + c2 a dE/dx + c3 a dE/dy ] ,
    ! E' = E - [E]: eddies are zonally elongated where they grow and meridionally elongated where
    ! they decay along the storm track. c1..c3 are fitted to 300 hPa high-pass (<9 d) NCEP R1
    ! eddy statistics, DJF 1981-2010. Offline with the observed EKE the forcing adds the
    ! storm-track exit troughs (Aleutian, Icelandic) and the European ridge.
    !
    ! The response to the azonal latent heating is the baroclinic Gill term (gill_response), not a
    ! barotropic Rossby wave source here: the equivalent barotropic conversion maps the upper-level
    ! anticyclone over the heating to the surface, i.e. a high where the observed surface low is.

    use, intrinsic :: iso_c_binding

    implicit none

    include 'fftw3.f03'

    real(wp), intent(in ) :: zsa(:,:)
    real(wp), intent(in ) :: u3f(:,:,:)
    real(wp), intent(in ) :: sam(:,:)
    real(wp), intent(out) :: aslp_sw(:,:)
    real(wp), intent(out) :: aslp_eddy(:,:)

    integer :: i, j, n, nmax
    real(wp) :: dphi, ri
    real(wp) :: phi(jm), uup(jm), ulo(jm), ue(jm), zbar(jm), qy(jm)
    real(wp) :: cph(0:jm)            ! cos at the half latitudes between j and j+1 (0 at the poles)
    real(dp) :: frc(im,jm)
    real(dp), dimension(im) :: rowr
    complex(dp), dimension(im/2+1) :: rowc
    complex(dp) :: frc_hat(im/2+1,jm), psi_hat(im/2+1,jm)
    complex(dp) :: frce_hat(im/2+1,jm), psie_hat(im/2+1,jm)
    real(dp) :: frce(im,jm)
    integer :: ip, imn
    real(wp) :: samz, dsdx, dsdy, tap, sgn
    real(wp) :: hfr
    real(wp) :: zeff(im,jm)
    integer :: nmaxe
    real(wp), dimension(im,jm) :: eddm, edduv, mfx, mfy
    ! closure coefficients for (M, u'v') on (E', a dE/dx, a dE/dy), NCEP R1 300 hPa high-pass DJF fit (NH);
    ! in the SH the terms odd under y -> -y (c_edd_m(3), c_edd_uv(1:2)) change sign
    real(wp), parameter :: c_edd_m(3)  = [-0.29_wp, 0.039_wp, 0.015_wp]
    real(wp), parameter :: c_edd_uv(3) = [ 0.16_wp,-0.056_wp, 0.015_wp]
    complex(dp) :: al(jm), ad(jm), au(jm), rhs(jm), x(jm)
    type(C_PTR) :: plan_r2c, plan_c2r

    nmax = min(n_aslp_sw, im/2-1)
    ! the eddy forcing is ~ the third derivative of E, which amplifies the small-scale noise of the model EKE:
    ! keep only the planetary-scale part (equivalent to truncating E, the closure being linear in E)
    nmaxe = min(n_aslp_sw_eddy, nmax)
    ri = 1._wp/c_aslp_sw_tau
    dphi = dy/r_earth

    do j=1,jm
      phi(j) = asin(sint(j))
      ! u3f, not u3: the thermal wind in u3 is damped towards the poles and the equator, which removes
      ! the subtropical jet waveguide and weakens the subpolar westerlies the waves propagate on
      uup(j) = sum(u3f(:,j,k300))*aim
      ulo(j) = sum(u3f(:,j,k850))*aim
      ue(j)  = max(c_aslp_sw_umin, uup(j))
    enddo
    cph(0)  = 0._wp
    cph(jm) = 0._wp
    do j=1,jm-1
      cph(j) = cos(0.5_wp*(phi(j)+phi(j+1)))
    enddo

    ! gradient of zonal mean absolute vorticity, orientation independent (centred differences in phi)
    zbar(:) = 0._wp
    do j=2,jm-1
      zbar(j) = -(uup(j+1)*cost(j+1)-uup(j-1)*cost(j-1))/(phi(j+1)-phi(j-1)) / (r_earth*cost(j))
    enddo
    do j=1,jm
      qy(j) = 2._wp*omega*cost(j)/r_earth
    enddo
    do j=2,jm-1
      qy(j) = qy(j) + (zbar(j+1)-zbar(j-1))/(phi(j+1)-phi(j-1)) / r_earth
    enddo

    ! effective height of the topography seen by the forcing flow (c_aslp_sw_fr>0, Froude blocking).
    ! Where N h/U >> 1 the low-level flow is blocked and goes around the obstacle instead of over it, so the
    ! vortex stretching saturates at a height of order U/N (e.g. Valdes & Hoskins 1991):
    !   h_eff = h0 tanh(h/h0) ,  h0 = c_aslp_sw_fr * max(|U_850|, c_aslp_sw_umin) / N .
    ! Offline (analysis/aslp_topo/froude.jl, c = 2, tau = 4 d) this raises the PI correlation with the ERA dynamic SLP
    ! part (DJF/MAM/SON 0.52/0.26/0.31 -> 0.70/0.54/0.59) and brings the ratio of the PI to the LGM-PI (CMIP5) amplitude
    ! deficit from 3.9/1.4/1.2 to 1.9/1.4/0.9: the linear forcing over the 3 km ice sheets is otherwise far too strong.
    do j=1,jm
      if (c_aslp_sw_fr.gt.0._wp) then
        hfr = c_aslp_sw_fr*max(abs(ulo(j)),c_aslp_sw_umin)/c_aslp_sw_nbv
        do i=1,im
          zeff(i,j) = hfr*tanh(max(zsa(i,j),0._wp)/hfr)
        enddo
      else
        zeff(:,j) = zsa(:,j)
      endif
    enddo

    ! topographic forcing
    frc(:,:) = 0._dp
    do j=2,jm-1
      do i=1,im
        frc(i,j) = -ulo(j)*2._wp*omega*sint(j)/c_aslp_sw_h &
          * (zeff(modulo(i,im)+1,j)-zeff(modulo(i-2,im)+1,j))/(2._wp*dxt(j))
      enddo
    enddo

    ! transient-eddy forcing
    frce(:,:) = 0._dp
    if (c_aslp_sw_eddy.gt.0._wp) then
      ! azonal anisotropy and momentum flux from the synoptic kinetic energy
      eddm(:,:) = 0._wp
      edduv(:,:) = 0._wp
      do j=2,jm-1
        samz = sum(sam(:,j))*aim
        sgn = sign(1._wp, sint(j))
        do i=1,im
          dsdx = r_earth*(sam(modulo(i,im)+1,j)-sam(modulo(i-2,im)+1,j))/(2._wp*dxt(j))
          dsdy = (sam(i,j+1)-sam(i,j-1))/(phi(j+1)-phi(j-1))
          eddm(i,j)  = c_aslp_sw_eddy*(c_edd_m(1)*(sam(i,j)-samz)  + c_edd_m(2)*dsdx  + sgn*c_edd_m(3)*dsdy)
          edduv(i,j) = c_aslp_sw_eddy*(sgn*c_edd_uv(1)*(sam(i,j)-samz) + sgn*c_edd_uv(2)*dsdx + c_edd_uv(3)*dsdy)
        enddo
      enddo
      ! eddy momentum flux convergence, u'u' = M, v'v' = -M (azonal parts only)
      mfx(:,:) = 0._wp
      mfy(:,:) = 0._wp
      do j=2,jm-1
        do i=1,im
          ip  = modulo(i,im)+1
          imn = modulo(i-2,im)+1
          mfx(i,j) = -( (eddm(ip,j)-eddm(imn,j))/(2._wp*dxt(j)) &
            + (edduv(i,j+1)*cost(j+1)**2-edduv(i,j-1)*cost(j-1)**2)/((phi(j+1)-phi(j-1))*r_earth*cost(j)**2) )
          mfy(i,j) = -( (edduv(ip,j)-edduv(imn,j))/(2._wp*dxt(j)) &
            - (eddm(i,j+1)*cost(j+1)-eddm(i,j-1)*cost(j-1))/((phi(j+1)-phi(j-1))*r_earth*cost(j)) &
            + eddm(i,j)*sint(j)/(cost(j)*r_earth) )
        enddo
      enddo
      ! curl of the eddy momentum flux convergence, tapered to zero equatorward of 20 deg and, where the
      ! 1/cos metric terms amplify the closure error (2x ERA5 at 77.5N), poleward of 80 deg
      do j=3,jm-2
        tap = min(1._wp, max(0._wp, (abs(phi(j))*180._wp/pi-20._wp)/10._wp), &
                  max(0._wp, (80._wp-abs(phi(j))*180._wp/pi)/10._wp))
        do i=1,im
          ip  = modulo(i,im)+1
          imn = modulo(i-2,im)+1
          frce(i,j) = tap * ( (mfy(ip,j)-mfy(imn,j))/(2._wp*dxt(j)) &
            - (mfx(i,j+1)*cost(j+1)-mfx(i,j-1)*cost(j-1))/((phi(j+1)-phi(j-1))*r_earth*cost(j)) )
        enddo
      enddo
    endif

    plan_r2c = fftw_plan_dft_r2c_1d(im, rowr, rowc, FFTW_ESTIMATE)
    plan_c2r = fftw_plan_dft_c2r_1d(im, rowc, rowr, FFTW_ESTIMATE)

    do j=1,jm
      rowr(:) = frc(:,j)
      call fftw_execute_dft_r2c(plan_r2c, rowr, rowc)
      frc_hat(:,j) = rowc(:)
      rowr(:) = frce(:,j)
      call fftw_execute_dft_r2c(plan_r2c, rowr, rowc)
      frce_hat(:,j) = rowc(:)
    enddo

    psi_hat(:,:) = (0._dp,0._dp)
    psie_hat(:,:) = (0._dp,0._dp)
    do n=1,nmax
      do j=1,jm
        ! Laplacian of wavenumber n, tridiagonal in latitude
        al(j) = cph(j-1)/(r_earth**2*dphi**2*cost(j))
        au(j) = cph(j)  /(r_earth**2*dphi**2*cost(j))
        ad(j) = -(cph(j-1)+cph(j))/(r_earth**2*dphi**2*cost(j)) - real(n,wp)**2/(r_earth*cost(j))**2
        ! (i n U/(a cos) + r) lap + i n q_y/(a cos)
        al(j) = cmplx(ri, real(n,wp)*ue(j)/(r_earth*cost(j)), dp) * al(j)
        au(j) = cmplx(ri, real(n,wp)*ue(j)/(r_earth*cost(j)), dp) * au(j)
        ad(j) = cmplx(ri, real(n,wp)*ue(j)/(r_earth*cost(j)), dp) * ad(j) &
          + cmplx(0._dp, real(n,wp)*qy(j)/(r_earth*cost(j)), dp)
        rhs(j) = frc_hat(n+1,j)
      enddo
      call tridiag_c(al, ad, au, rhs, x)
      psi_hat(n+1,:) = x(:)
      if (c_aslp_sw_eddy.gt.0._wp .and. n.le.nmaxe) then
        call tridiag_c(al, ad, au, frce_hat(n+1,:), x)
        psie_hat(n+1,:) = x(:)
      endif
    enddo

    do j=1,jm
      rowc(:) = psi_hat(:,j)
      call fftw_execute_dft_c2r(plan_c2r, rowc, rowr)
      if (abs(phi(j)).lt.c_aslp_sw_lateq*pi/180._wp) then
        aslp_sw(:,j) = 0._wp
      else
        ! f*psi/g is the height anomaly, aim the FFTW normalisation
        aslp_sw(:,j) = c_aslp_sw_amp*ra*g * 2._wp*omega*sint(j)/g * rowr(:)*aim
      endif
    enddo

    ! transient-eddy part, included in aslp_sw
    aslp_eddy(:,:) = 0._wp
    if (c_aslp_sw_eddy.gt.0._wp) then
      do j=1,jm
        rowc(:) = psie_hat(:,j)
        call fftw_execute_dft_c2r(plan_c2r, rowc, rowr)
        if (abs(phi(j)).ge.c_aslp_sw_lateq*pi/180._wp) then
          aslp_eddy(:,j) = c_aslp_sw_amp*ra*g * 2._wp*omega*sint(j)/g * rowr(:)*aim
        endif
      enddo
      aslp_sw(:,:) = aslp_sw(:,:) + aslp_eddy(:,:)
    endif

    call fftw_destroy_plan(plan_r2c)
    call fftw_destroy_plan(plan_c2r)

    return

  end subroutine stationary_waves


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  gill_response
  !   Purpose    :  sea level pressure of the first baroclinic mode response
  !                 to the azonal latent heating (Gill 1980) on the sphere
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine gill_response(prc, ficz, aslp_gill)

    ! Deep (latent) heating warms the column above the surface layer, so the surface low it
    ! produces is invisible to aslp_temp, which only sees the surface-anchored tsksl. In ERA-Interim
    ! the JJA low over E China and the W Pacific monsoon flank is carried by the warm free troposphere
    ! (azonal T +0.2 K at 850, +2.5 K at 500, +3.9 K at 300 hPa over E China) and by a(Z500), not by the
    ! surface temperature. The steady, damped, linear shallow-water equations of the first baroclinic
    ! mode (Gill 1980; Matsuno 1966), with phi the lower-level geopotential anomaly,
    !   eps u - f v = -(1/(a cos)) d(phi)/dlambda
    !   eps v + f u = -(1/a) d(phi)/dlat
    !   eps phi + c^2 div(u,v) = -Q ,
    ! give the surface low beneath and west of the heating, the monsoon trough and the subtropical
    ! high to the east of it (Rodwell & Hoskins 2001). The heating forcing is the thickness tendency
    ! of the 1000-500 hPa layer heated by the azonal latent heat release,
    !   Q = Rd ln2 g Le P' / (cp dp) ,  dp = 500 hPa ,
    ! tapered to zero beyond c_aslp_gill_lat0..lat1, where precipitation is frontal and the heating is
    ! balanced by horizontal advection instead (Hoskins & Karoly 1981). With i_aslp_gill_lat=2 the
    ! distance is taken from the ITCZ (ficz) instead of the equator, so the window follows the tropical
    ! belt into the summer hemisphere (Asian monsoon, Meiyu-Baiu) and excludes the winter-hemisphere
    ! subtropical storm-track entrances, which a fixed window counts as deep heating. Offline (CX prc,
    ! ERA surface-anchored thermal part, Asia-Pacific 0-50N r) 25/35 deg from ficz vs 30/40 from the
    ! equator: DJF 0.71 vs 0.70 (no Gill 0.72), MAM 0.64 vs 0.66, JJA 0.95 vs 0.91, SON 0.61 vs 0.54.
    ! In Fourier space each zonal wavenumber n is one banded system in latitude (see below).
    ! aslp_gill = c_aslp_gill * ra * phi. Offline against ERA-Interim, JJA, 0-50N 40E-120W, adding it to the
    ! surface-anchored thermal part raises the correlation with the azonal msl from 0.77 to 0.86-0.90
    ! for c_aslp_gill ~ 1 (eps^-1 = 3-10 d, c = 30-50 m/s); the DJF correlation is unchanged.

    use, intrinsic :: iso_c_binding

    implicit none

    include 'fftw3.f03'

    real(wp), intent(in ) :: prc(:,:)
    real(wp), intent(in ) :: ficz          !! ITCZ latitude (radians)
    real(wp), intent(out) :: aslp_gill(:,:)

    integer :: i, j, n, nmax, f, jj, k, r
    real(wp) :: eps, c2, dphis, alat, tap, prcz, rn
    real(wp) :: phi(jm), cph(jm-1), fcj(jm), fcf(jm-1)
    real(dp), dimension(im) :: rowr
    complex(dp), dimension(im/2+1) :: rowc
    complex(dp) :: q_hat(im/2+1,jm), phi_hat(im/2+1,jm)
    complex(dp) :: amat(2*jm-1,2*jm-1), rhs(2*jm-1), x(2*jm-1)
    real(dp) :: frc(im,jm)
    type(C_PTR) :: plan_r2c, plan_c2r

    nmax  = min(n_aslp_gill, im/2-1)
    eps   = 1._wp/(c_aslp_gill_tau*86400._wp)
    c2    = c_aslp_gill_c**2

    do j=1,jm
      phi(j) = asin(sint(j))
    enddo
    dphis = phi(2)-phi(1)     ! signed: v and the latitude derivatives are taken northward whatever the grid orientation
    do j=1,jm
      fcj(j) = 2._wp*omega*sint(j)
    enddo
    do f=1,jm-1
      cph(f) = cos(0.5_wp*(phi(f)+phi(f+1)))
      fcf(f) = 2._wp*omega*sin(0.5_wp*(phi(f)+phi(f+1)))
    enddo

    ! heating forcing, tapered with distance from the equator (i_aslp_gill_lat=1) or from the ITCZ (2)
    frc(:,:) = 0._dp
    do j=1,jm
      if (i_aslp_gill_lat.eq.2) then
        alat = abs(phi(j)-ficz)*180._wp/pi
      else
        alat = abs(phi(j))*180._wp/pi
      endif
      if (alat.le.c_aslp_gill_lat0) then
        tap = 1._wp
      else if (alat.ge.c_aslp_gill_lat1) then
        tap = 0._wp
      else
        tap = 0.5_wp*(1._wp+cos(pi*(alat-c_aslp_gill_lat0)/(c_aslp_gill_lat1-c_aslp_gill_lat0)))
      endif
      prcz = sum(prc(:,j))*aim
      frc(:,j) = tap*Rd*log(2._wp)*g*Le*(prc(:,j)-prcz)/(cp*5.e4_wp)
    enddo

    plan_r2c = fftw_plan_dft_r2c_1d(im, rowr, rowc, FFTW_ESTIMATE)
    plan_c2r = fftw_plan_dft_c2r_1d(im, rowc, rowr, FFTW_ESTIMATE)

    do j=1,jm
      rowr(:) = frc(:,j)
      call fftw_execute_dft_r2c(plan_r2c, rowr, rowc)
      q_hat(:,j) = rowc(:)
    enddo

    ! C grid in latitude, as the 2D C-grid solution it reproduces (converged at 5 deg: r 0.999 against 0.5 deg):
    ! phi and u on the cell latitudes, v on the half latitude f between rows f and f+1 (v = 0 at the poles).
    ! u is local and is eliminated, u_j = [ f_j (v_{j-1/2}+v_{j+1/2})/2 - i n phi_j/(a cos_j) ] / eps;
    ! the remaining unknowns are ordered x = (phi_1, v_1, phi_2, v_2, ..., phi_jm), a banded system
    ! solved with Gaussian elimination and partial pivoting.
    phi_hat(:,:) = (0._dp,0._dp)
    ! the wavenumbers are independent; each dense solve is O(jm^3)
    !$omp parallel do schedule(dynamic,1) private(n,rn,amat,rhs,x,j,r,k,f,jj)
    do n=1,nmax
      rn = real(n,wp)
      amat(:,:) = (0._dp,0._dp)
      ! continuity on row j: eps phi + c^2 [ i n u/(a cos) + (v_{j+1/2} cph - v_{j-1/2} cph)/(a cos dphis) ] = -Q
      do j=1,jm
        r = 2*j-1
        amat(r,r) = amat(r,r) + eps + c2*rn**2/(eps*(r_earth*cost(j))**2)
        do k=0,1
          f = j-1+k     ! k=0: half latitude below row j (towards row j-1), k=1: above (towards row j+1)
          if (f.ge.1 .and. f.le.jm-1) then
            amat(r,2*f) = amat(r,2*f) + c2*cmplx((2*k-1)*cph(f)/(r_earth*cost(j)*dphis), &
              rn*fcj(j)/(2._wp*eps*r_earth*cost(j)), dp)
          endif
        enddo
        rhs(r) = -q_hat(n+1,j)
      enddo
      ! meridional momentum on half latitude f: eps v + f (u_f + u_{f+1})/2 = -(phi_{f+1}-phi_f)/(a dphis)
      do f=1,jm-1
        r = 2*f
        amat(r,r) = amat(r,r) + eps
        do jj=f,f+1
          amat(r,2*jj-1) = amat(r,2*jj-1) + cmplx(0._wp, -0.5_wp*fcf(f)*rn/(eps*r_earth*cost(jj)), dp)
          do k=jj-1,jj
            if (k.ge.1 .and. k.le.jm-1) amat(r,2*k) = amat(r,2*k) + 0.25_wp*fcf(f)*fcj(jj)/eps
          enddo
        enddo
        amat(r,2*f+1) = amat(r,2*f+1) + 1._wp/(r_earth*dphis)
        amat(r,2*f-1) = amat(r,2*f-1) - 1._wp/(r_earth*dphis)
        rhs(r) = (0._dp,0._dp)
      enddo
      call solve_c_dense(amat, rhs, x)
      do j=1,jm
        phi_hat(n+1,j) = x(2*j-1)
      enddo
    enddo
    !$omp end parallel do

    do j=1,jm
      rowc(:) = phi_hat(:,j)
      call fftw_execute_dft_c2r(plan_c2r, rowc, rowr)
      aslp_gill(:,j) = c_aslp_gill*ra*rowr(:)*aim    ! aim accounts for the FFTW normalisation
    enddo

    call fftw_destroy_plan(plan_r2c)
    call fftw_destroy_plan(plan_c2r)

    return

  end subroutine gill_response


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  solve_c_dense
  !   Purpose    :  complex linear system, Gaussian elimination with partial pivoting
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine solve_c_dense(a, b, x)

    implicit none

    complex(dp), intent(inout) :: a(:,:)   !! matrix, overwritten
    complex(dp), intent(inout) :: b(:)     !! right hand side, overwritten
    complex(dp), intent(out  ) :: x(:)

    integer :: k, i, p, nk
    complex(dp) :: tmp(size(b)), t, fac

    nk = size(b)
    do k=1,nk-1
      p = k-1+maxloc(abs(a(k:nk,k)),1)
      if (p.ne.k) then
        tmp(:) = a(k,:); a(k,:) = a(p,:); a(p,:) = tmp(:)
        t = b(k); b(k) = b(p); b(p) = t
      endif
      do i=k+1,nk
        if (a(i,k).ne.(0._dp,0._dp)) then
          fac = a(i,k)/a(k,k)
          a(i,k:nk) = a(i,k:nk) - fac*a(k,k:nk)
          b(i) = b(i) - fac*b(k)
        endif
      enddo
    enddo
    do k=nk,1,-1
      x(k) = (b(k) - sum(a(k,k+1:nk)*x(k+1:nk)))/a(k,k)
    enddo

    return

  end subroutine solve_c_dense


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  tridiag_c
  !   Purpose    :  complex tridiagonal solver (Thomas algorithm)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine tridiag_c(a, b, c, d, x)

    implicit none

    complex(dp), intent(in ) :: a(:), b(:), c(:), d(:)   !! sub-, main and super-diagonal, right hand side
    complex(dp), intent(out) :: x(:)

    integer :: k, nk
    complex(dp) :: cpr(size(b)), dpr(size(b)), den

    nk = size(b)
    cpr(1) = c(1)/b(1)
    dpr(1) = d(1)/b(1)
    do k=2,nk
      den    = b(k)-a(k)*cpr(k-1)
      cpr(k) = c(k)/den
      dpr(k) = (d(k)-a(k)*dpr(k-1))/den
    enddo
    x(nk) = dpr(nk)
    do k=nk-1,1,-1
      x(k) = dpr(k)-cpr(k)*x(k+1)
    enddo

    return

  end subroutine tridiag_c


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  zslp
  !   Purpose    :  compute zonally averaged sea level pressure component
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine zslp(sin_cos_acbar, tsl, aslp, pzsa, fdydse, t3, q3, htrop, ttrop, efe, &
          ficz, mtz, slp, had_fi, had_width)

    implicit none

    real(wp), intent(in ) :: sin_cos_acbar(:,:)
    real(wp), intent(in ) :: tsl(:,:)
    real(wp), intent(in ) :: aslp(:,:)
    real(wp), intent(in ) :: pzsa(:,:)     !! relative surface pressure of the (smoothed) topography, exp(-zsa_smooth/hatm)
    real(wp), intent(in ) :: fdydse(:,:)
    real(wp), intent(in ) :: t3(:,:,:)
    real(wp), intent(in ) :: q3(:,:,:)
    real(wp), intent(in ) :: htrop(:,:)
    real(wp), intent(in ) :: ttrop(:,:)
    real(wp), intent(in ) :: efe          !! energy flux equator from the previous step (deg), see efe_diag

    real(wp), intent(out) :: ficz         !! ITCZ latitude of the cell geometry (radians)
    real(wp), intent(out) :: mtz(:)       !! zonal mean mountain torque [p_s dz_s/dx] on T-points (N/m2)
    real(wp), intent(out) :: slp(:,:)
    real(wp), intent(out) :: had_fi
    real(wp), intent(out) :: had_width

    real(wp), parameter :: p6=pi/6._wp
    real(wp), parameter :: dse_min=1.e3_wp   !! J/kg, floor on the branch DSE contrast

    integer :: i, j, k
    integer :: jc, jtn1, jtn2, jpn1, jpn2, jts1, jts2, jps1, jps2
    real(wp) :: wtn1, wtn2, wpn1, wpn2, wts1, wts2, wps1, wps2
    real(wp) :: tnh, tsh, scosh
    real(wp) :: dtpn, dtps, dtfn, dtfs, ttrpmx, tbhn, tbhs, dthn, dths, cocfn, cocfs
    real(wp) :: hadwidsc
    real(wp) :: hadwid
    real(wp) :: ff, coc, psum, csum
    real(wp) :: dzdse, aferzn, aferzs
    real(wp) :: qtrp, wtrp, fqhad
    real(wp) :: htrpz, sbotz, stopz, tsfcz, wwid, hdv

    real(wp) :: fisu(jmc), fist(jm), vsz(jmc), tslz(jm), psz(jmc), acbarz(jm)
    real(wp) :: dsez(jm), fedz(jmc), aferz(jm)
    real(wp) :: mtu, sgnf


    ! zonal mean sea level temperature
    tslz(:) = 0._wp
    do j=1,jm 
      do i=1,im
        tslz(j) = tslz(j) + tsl(i,j)*aim 
      enddo        
    enddo

    ! NH and SH mean sea level temperatures
    tnh = 0._wp
    tsh = 0._wp
    scosh = 0._wp
    do j=1,jm 
      do i=1,im
        if (j.le.jeq) then 
          tnh = tnh + tsl(i,j)*cost(j)
          scosh = scosh + cost(j)
        else
          tsh = tsh + tsl(i,j)*cost(j)
        endif        
      enddo
    enddo 
    tnh = tnh/scosh     
    tsh = tsh/scosh

    ! zonally averaged cross-isobar angle 
    acbarz(:) = 0._wp
    do j=2,jm      
      do i=1,im
        acbarz(j) = acbarz(j) + 0.5_wp*(sin_cos_acbar(i,j)+sin_cos_acbar(i,j-1)) * aim
      enddo
    enddo

    ! maximum tropical sea level temperature, drives the Hadley branches
    ttrpmx = 0._wp
    do j=jtn,jts+1
      ttrpmx = max(ttrpmx,tslz(j))
    enddo 

    ! ITCZ position, depends on temperature difference between the two hemispheres 

    ! ITCZ latitude, selected by i_mmc_itcz.  The zero of the cell phase ff below, i.e. the boundary
    ! between the two Hadley cells, sits exactly at fi = ficz.
    !
    ! i_mmc_itcz = 0  Proportional to the interhemispheric sea level temperature contrast,
    !                   ficz = c_mmc_2*(T_NH - T_SH).
    !                 One coefficient serves the seasonal cycle and forced changes alike, so the ITCZ moves
    !                 4-10 times too far under LGM and 4xCO2 forcing, whose contrast sits at high latitudes.
    !
    ! i_mmc_itcz = 1  Energy flux equator closure: ficz = c_mmc_efe*efe.  efe is the zero of the northward
    !                 atmospheric energy transport implied by the boundary fluxes and the column storage
    !                 (efe_diag), so it does not depend on the circulation this routine prescribes.  It is
    !                 from the previous step (from the restart file on the first step; 0 after a cold start).
    if (i_mmc_itcz.eq.1) then
      ficz = c_mmc_efe*efe*pi/180._wp
    else
      ficz = c_mmc_2*(tnh-tsh)
    endif

    ! Width of the Hadley cells: c_mmc_dt0 degrees in total, scaled by the Held (2000) supercriticality
    ! criterion.  The Hadley cell ends where the angular momentum conserving flow first becomes
    ! baroclinically unstable, at
    !   fi_H ~ [ g H Dv / (Omega**2 a**2) ]**(1/4),
    ! with H the tropopause height and Dv = (theta(H)-theta(sfc))/theta_sfc the bulk dry stability of the
    ! tropical troposphere.  Only H*Dv varies here, so
    !   hadwid = c_mmc_dt0 * (H*Dv/hdv_mmc_ref)**c_mmc_dt1,
    ! and hadwidsc follows from width = 180/(3*hadwidsc) deg.
    ! Both factors grow under warming - the tropopause rises and the troposphere gets more stable - which is
    ! the accepted explanation of the observed expansion.  Note that Dv is the DRY stability up to the tropopause
    ! and does rise here, unlike the gross moist stability between fixed levels.
    htrpz = 0._wp
    stopz = 0._wp
    sbotz = 0._wp
    tsfcz = 0._wp
    wwid  = 0._wp
    do j=jtn,jts+1
      do i=1,im
        ! ttrop is t_prof evaluated at htrop: t_prof clamps its argument at the tropopause,
        ! so vesta's profile loop leaves it holding exactly the tropopause temperature.
        htrpz = htrpz + htrop(i,j)*cost(j)
        stopz = stopz + (ttrop(i,j)+g*htrop(i,j)/cp)*cost(j)
        sbotz = sbotz + (t3(i,j,1)+g*zc(1)/cp)*cost(j)
        tsfcz = tsfcz + t3(i,j,1)*cost(j)
        wwid  = wwid  + cost(j)
      enddo
    enddo
    htrpz = htrpz/wwid
    stopz = stopz/wwid
    sbotz = sbotz/wwid
    tsfcz = tsfcz/wwid
    hdv   = htrpz*(stopz-sbotz)/tsfcz
    hadwid = c_mmc_dt0*(max(hdv,1.e-20_wp)/hdv_mmc_ref)**c_mmc_dt1

    hadwidsc = 60._wp/max(hadwid,1._wp)

    hadwidsc = max(hadwidsc,0.5_wp)
    hadwidsc = min(hadwidsc,1.5_wp)       

    ! position of borders of cells

    do j=1,jm
      fisu(j) = cell_phase(fiu(j),ficz,hadwidsc)
      fist(j) = cell_phase(fit(j),ficz,hadwidsc)
    enddo
    do j=1,jm
      ! N boundary between Hadley and Ferrel cells
      if (fisu(j).ge.pi .and. fisu(j+1).lt.pi) then
        jc= j
        if (fist(jc).lt.pi) then
          jtn1 = jc-1
          jtn2 = jc
        else
          jtn1 = jc
          jtn2 = jc+1
        endif
        wtn2 = 1._wp-(pi-fist(jtn2))/(fist(jtn1)-fist(jtn2))
        wtn1 = 1._wp-wtn2
      endif
      ! N boundary between Ferrel and Polar cells
      if (fisu(j).ge.2._wp*pi .and. fisu(j+1).lt.2._wp*pi) then
        jc= j
        if (fist(jc).lt.2.*pi) then
          jpn1 = jc-1
          jpn2 = jc
        else
          jpn1 = jc
          jpn2 = jc+1
        endif
        wpn2 = 1._wp-(2.*pi-fist(jpn2))/(fist(jpn1)-fist(jpn2))
        wpn1 = 1._wp-wpn2
      endif
      ! S boundary between Hadley and Ferrel cells
      if (fisu(j).ge.-pi .and. fisu(j+1).lt.-pi) then
        jc= j
        if (fist(jc).lt.-pi) then
          jts1 = jc-1
          jts2 = jc
        else
          jts1 = jc
          jts2 = jc+1
        endif
        wts2 = 1._wp-(-pi-fist(jts2))/(fist(jts1)-fist(jts2))
        wts1 = 1._wp-wts2
      endif
      ! S boundary between Ferrel and Polar cells
      if (fisu(j).ge.-2.*pi .and. fisu(j+1).lt.-2.*pi) then
        jc= j
        if (fist(jc).lt.-2.*pi) then
          jps1 = jc-1
          jps2 = jc
        else
          jps1 = jc
          jps2 = jc+1
        endif
        wps2 = 1._wp-(-2.*pi-fist(jps2))/(fist(jps1)-fist(jps2))
        wps1 = 1._wp-wps2
      endif
    enddo

    ! Hadley cell width and ITCZ position 
    had_fi = 0.5_wp*((wtn1*fit(jtn1)+wtn2*fit(jtn2)) + (wts1*fit(jts1)+wts2*fit(jts2)))
    had_width = (wtn1*fit(jtn1)+wtn2*fit(jtn2)) - (wts1*fit(jts1)+wts2*fit(jts2)) 

    ! Sea level temperature contrasts across the cells

    ! temperature gradients in the polar cells
    dtpn = (0.5_wp*tslz(jpn)+0.5_wp*tslz(jpn+1)) - tslz(1)
    dtps = (0.5_wp*tslz(jps)+0.5_wp*tslz(jps+1)) - tslz(jm)

    ! temperature gradients in the ferrel cells
    dtfn = (0.5*tslz(jtn)+0.5*tslz(jtn+1)) - (0.5*tslz(jpn)+0.5*tslz(jpn+1))
    dtfs = (0.5*tslz(jts)+0.5*tslz(jts+1)) - (0.5*tslz(jps)+0.5*tslz(jps+1))

    ! temperature gradients in Hadley cells    
    tbhn = 0.5*tslz(jtn)+0.5*tslz(jtn+1)
    tbhs = 0.5*tslz(jts)+0.5*tslz(jts+1)

    dthn = max(0._wp,ttrpmx-tbhn)
    dths = max(0._wp,ttrpmx-tbhs)

    ! Thermodynamic scaling of the Hadley amplitude, selected by i_mmc_had
    !
    ! i_mmc_had = 0  The amplitude is the sea level temperature contrast across the cell alone,
    !                  psi_H = c_mmc_had * dth.
    !                This has no thermodynamic damping, so the cells follow whatever the near surface
    !                temperature contrast does and strengthen under warming, opposite to CMIP.
    !
    ! i_mmc_had = 1  The same contrast, scaled by the moisture mass budget of the cell.  The lower branch
    !                imports water vapour at a rate psi_H*q and this has to balance the net precipitation
    !                of the ascending branch, which is set by the radiative cooling and is close to
    !                invariant, so
    !                  psi_H ~ (P-E)/q  ->  psi_H = c_mmc_had * dth * (q_mmc_ref/qtrp)**c_mmc_q,
    !                with qtrp the tropical mean boundary layer specific humidity.  This is the
    !                Held & Soden (2006) weakening of the tropical overturning: q rises at ~5 %/K while
    !                the net precipitation is nearly fixed, so the mass flux has to fall.
    if (i_mmc_had.eq.1) then
      qtrp = 0._wp
      wtrp = 0._wp
      do j=jtn,jts+1
        do i=1,im
          qtrp = qtrp + q3(i,j,1)*cost(j)
          wtrp = wtrp + cost(j)
        enddo
      enddo
      qtrp = qtrp/wtrp
      fqhad = (q_mmc_ref/max(qtrp,1.e-6_wp))**c_mmc_q
      dthn = dthn*fqhad
      dths = dths*fqhad
    endif

    ! Amplitude of the Ferrel branches, selected by i_mmc_fer
    !
    ! i_mmc_fer = 1  The amplitude is proportional to the meridional sea level temperature contrast ACROSS the cell,
    !                  psi_F = c_mmc_fer * dtf,   dtf = tslz(Hadley edge) - tslz(polar edge)
    !
    ! i_mmc_fer = 2  Transformed Eulerian mean form, resolved in latitude.  The Ferrel cell is the Eulerian residue of the baroclinic eddy fluxes,
    !                  psi_F ~ -(2 pi a cos(fi)/g) [v'th'] / (dth/dp),
    !                i.e. the mass circulation that returns, against the mean dry static energy gradient, the DSE the synoptic eddies transport poleward.  
    !                Written as a bulk balance over the depth of the cell that is
    !                  psi_F(j) = c_mmc_fer * |F_eddy(j)| / Ds(j),
    !                with F_eddy the zonally integrated synoptic eddy DSE flux (W) and Ds the DSE contrast between the two branches (J/kg).

    if (i_mmc_fer.eq.1) then

      cocfn = c_mmc_fer*dtfn
      cocfs = c_mmc_fer*dtfs

    else if (i_mmc_fer.eq.2) then

      ! Bulk dry static energy contrast between the lower and the upper branch of the cells.
      ! It is the denominator of the Ferrel amplitude below. 
      dzdse = zc(k_dse_up)-zc(k_dse_lo)
      do j=1,jm
        dsez(j) = 0._wp
        do i=1,im
          dsez(j) = dsez(j) + (cp*(t3(i,j,k_dse_up)-t3(i,j,k_dse_lo)) + g*dzdse)*aim
        enddo
        dsez(j) = max(dsez(j),dse_min)
      enddo

      ! zonally integrated transient eddy DSE flux, W (fdydse carries kg/s*K)
      fedz(:) = 0._wp
      do j=1,jmc
        do i=1,im
          fedz(j) = fedz(j) + fdydse(i,j)*cp
        enddo
      enddo

      do j=2,jm
        aferz(j) = c_mmc_fer * abs(fedz(j))/(0.5_wp*(dsez(j-1)+dsez(j))) * 1.e-10_wp
      enddo

      ! average over latitudinal belt
      aferzn = 0._wp
      do j=jpn+1,jtn
        aferzn = aferzn + aferz(j)/real(jtn-jpn,wp)
      enddo
      aferzs = 0._wp
      do j=jts+1,jps
        aferzs = aferzs + aferz(j)/real(jps-jts,wp)
      enddo
      do j=2,jm
        cocfn = aferzn
        cocfs = aferzs
      enddo

    endif

    ! Surface meridional ageostrofic wind 

    vsz(:) = 0._wp

    do j=2,jm

      ff = cell_phase(fiu(j),ficz,hadwidsc)

      coc = 0._wp

      if (ff.ge.0._wp .and. ff.lt.pi)              coc = c_mmc_had*dthn
      if (ff.ge.pi .and. ff.lt.2._wp*pi)           coc = cocfn
      if (ff.ge.2._wp*pi .and. ff.lt.3._wp*pi)     coc = c_mmc_pol*dtpn

      if (-ff.gt.0._wp .and. (-ff).le.pi)          coc = c_mmc_had*dths
      if (-ff.gt.pi .and. (-ff).le.2._wp*pi)       coc = cocfs
      if (-ff.gt.2._wp*pi .and. (-ff).le.3._wp*pi) coc = c_mmc_pol*dtps

      vsz(j) = -coc*sin(ff) 

    enddo

    ! Mountain torque
    !
    ! vsz prescribes the lower branch of the cells, i.e. its Coriolis torque f*M.  In the zonal mean that
    ! torque is balanced by the surface friction AND by the form drag of the orography,
    !   f*M = [tau] + [p_s dz_s/dx],
    ! and only the friction part comes with zonal mean surface westerlies and so with a meridional gradient
    ! of the zonal mean SLP.  Where the azonal SLP puts high pressure on the upstream flank of a mountain
    ! range (the LGM ice sheets, the Rockies), the lower branch is carried by the geostrophic flow against
    ! that pressure dipole instead: u3d's vgb over the layer mass above the terrain, whose zonal integral is
    ! exactly [p_s dz_s/dx]/f.  The zonal part of the SLP drops out of the zonal integral, and with
    ! p_s = slp*pzsa, z_s = -hatm*ln(pzsa),
    !   [p_s dz_s/dx] = hatm*[pzsa d(aslp)/dx].
    ! The friction (Ekman) part of the branch is vsz - MT/(f*M_pbl), with M_pbl = (1-pblu)*amas the mass of
    ! the PBL that carries the zonal mean ageostrophic flux in u3d.  In the SLP increment below f cancels,
    ! so the torque term needs no regularisation at the equator.
    do j=1,jm
      mtz(j) = 0._wp
      if (j.gt.1 .and. j.lt.jm) then
        do i=1,im
          mtz(j) = mtz(j) + pzsa(i,j)*(aslp(modulo(i,im)+1,j)-aslp(modulo(i-2,im)+1,j))/(2._wp*dxt(j))
        enddo
        mtz(j) = hatm*mtz(j)*aim
      endif
    enddo

    ! Zonally averaged SLP   
    ! Integration from NP to SP     

    psz(:) = 0._wp
    do j=2,jm
      mtu  = 0.5_wp*(mtz(j-1)+mtz(j))
      if (fiu(j).gt.0._wp) then
        sgnf = 1._wp
      else if (fiu(j).lt.0._wp) then
        sgnf = -1._wp
      else
        sgnf = 0._wp
      endif
      psz(j) = psz(j-1) + (vsz(j)*fcorua(j) - c_mmc_mt*sgnf*mtu/((1._wp-pblu(j))*amas))*ra/acbarz(j) * dy
    enddo

    ! SLP = zonal+azonal

    do i=1,im
      do j=1,jm
        slp(i,j)  = psz(j) + aslp(i,j)
      enddo
    enddo

    ! Restoring of atmospheric mass  

    psum = 0._wp
    csum = 0._wp
    do i=1,im
      do j=1,jm
        psum = psum + slp(i,j)*cost(j)
        csum = csum + cost(j)
      enddo    
    enddo 
    psum = psum/csum

    do i=1,im
      do j=1,jm
        slp(i,j)  = slp(i,j)  + p0-psum
      enddo 
    enddo 

    return

  end subroutine zslp


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  e f e _ d i a g
  !   Purpose    :  diagnostic energy flux equator (EFE), from boundary fluxes only
  !
  !   The net energy input into the atmospheric column is NEI = rb_atm + sha + lha, the moist static
  !   energy source of the column budget in time_step (latent heat enters through evaporation).
  !   Its storage is dE/dt with E = cheat*(tam+dtcol) + Le*wcon, the column enthalpy plus latent
  !   energy, differenced between successive calls.  The difference is taken in the column state,
  !   cheat*d(tam+dtcol) + Le*d(wcon) with the current cheat, so that a change of the column mass
  !   with the topography (atm_grid_update, or a restart into a new geometry) is not read as storage.  NEI - dE/dt, with its global mean removed so
  !   that the implied transport closes at both poles, integrated from the south pole gives the
  !   northward atmospheric energy transport F at every latitude boundary.  Nothing here uses the
  !   model's own transports, so the EFE is independent of the circulation that the MMC
  !   parameterisation prescribes.
  !
  !   aht_eq  : F across the equator (PW), positive northward
  !   efe     : latitude (deg) of the zero crossing of F closest to the equator within +-25 deg,
  !             linearly interpolated between v-points (the full EFE)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine efe_diag(rb_atm, sha, lha, cheat, tam, dtcol, wcon, sqr, dt, &
      tcol_prev, wcon_prev, aht_eq, efe)

    implicit none

    real(wp), intent(in   ) :: rb_atm(:,:)
    real(wp), intent(in   ) :: sha(:,:)
    real(wp), intent(in   ) :: lha(:,:)
    real(wp), intent(in   ) :: cheat(:,:)
    real(wp), intent(in   ) :: tam(:,:)
    real(wp), intent(in   ) :: dtcol(:,:)
    real(wp), intent(in   ) :: wcon(:,:)
    real(wp), intent(in   ) :: sqr(:,:)
    real(wp), intent(in   ) :: dt
    real(wp), intent(inout) :: tcol_prev(:,:)   !! tam+dtcol at the previous call (K), col_prev_unset before the first
    real(wp), intent(inout) :: wcon_prev(:,:)   !! wcon at the previous call (kg m-2)
    real(wp), intent(out  ) :: aht_eq
    real(wp), intent(out  ) :: efe

    integer :: i, j
    real(wp), parameter :: fi_efe_max = 25._wp   !! deg, search window for the zero crossing
    real(dp) :: atot, nglob
    real(dp) :: nrow(jm), arow(jm), f(jmc)
    real(wp) :: src(im,jm)
    real(wp) :: la, lb, x

    ! column energy storage since the previous call; zero storage on the first call (cold start, or a
    ! restart file without tcol_prev)
    where (tcol_prev.le.col_prev_unset)
      src = rb_atm + sha + lha
    elsewhere
      src = rb_atm + sha + lha - (cheat*(tam+dtcol-tcol_prev) + Le*(wcon-wcon_prev))/dt
    endwhere
    tcol_prev = tam+dtcol
    wcon_prev = wcon

    ! zonal integrals per row (W) and the global mean (W/m2)
    do j=1,jm
      nrow(j) = 0._dp
      arow(j) = 0._dp
      do i=1,im
        nrow(j) = nrow(j) + real(src(i,j),dp)*real(sqr(i,j),dp)
        arow(j) = arow(j) + real(sqr(i,j),dp)
      enddo
    enddo
    atot  = sum(arow)
    nglob = sum(nrow)/atot

    ! northward transport across the northern edge of each row, fiu(j), integrated from the south
    ! pole (j=1 is the northern-most row); f(jmc) = 0 at the south pole
    f(jmc) = 0._dp
    do j=jm,1,-1
      f(j) = f(j+1) + (nrow(j) - nglob*arow(j))
    enddo

    ! transport across the equator, fiu(jeq+1) = 0
    aht_eq = real(f(jeq+1),wp)*1.e-15_wp

    ! zero crossing of f closest to the equator
    efe = -999._wp
    do j=1,jm
      la = fiu(j)*180._wp/pi
      lb = fiu(j+1)*180._wp/pi
      if (abs(la).gt.fi_efe_max+1.e-6_wp .or. abs(lb).gt.fi_efe_max+1.e-6_wp) cycle
      if ((f(j).gt.0._dp .and. f(j+1).le.0._dp) .or. (f(j).lt.0._dp .and. f(j+1).ge.0._dp)) then
        x = lb + (la-lb)*real(-f(j+1)/(f(j)-f(j+1)),wp)
        if (efe.lt.-900._wp .or. abs(x).lt.abs(efe)) efe = x
      endif
    enddo
    ! no crossing in the window: the EFE lies beyond it, on the side opposite to the transport
    if (efe.lt.-900._wp) efe = -sign(fi_efe_max, aht_eq)

    return

  end subroutine efe_diag

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function :  c e l l _ p h a s e
  !   Purpose  :  phase of the meridional cell pattern at latitude fi;
  !               0 at the ITCZ, +-pi at the Hadley edges, +-2pi at the
  !               Ferrel/polar boundaries
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function cell_phase(fi, ficz, hadwidsc) result(ff)

    implicit none

    real(wp), intent(in) :: fi        !! latitude (radians)
    real(wp), intent(in) :: ficz      !! ITCZ latitude (radians)
    real(wp), intent(in) :: hadwidsc  !! width scale, one cell spans pi/(6*hadwidsc) radians
    real(wp) :: ff

    real(wp), parameter :: dfi_min = 5._wp*pi/180._wp   !! minimum latitudinal extent of a Hadley cell
    real(wp) :: hw, fen, fes

    ! i_mmc_geo = 0  The whole pattern is displaced towards the ITCZ by ficz/(1+c_mmc_1*(fi-ficz)**2).
    !                The displacement decays with the distance from the ITCZ, so the summer Hadley edge,
    !                which lies close to the ITCZ, moves by ~0.6-0.75*ficz and the winter edge by ~0.2-0.3*ficz.
    !                The summer subtropical high and trades end up ~10 deg too far poleward in JJA, while the
    !                observed edges move by only ~0.25-0.35 times the ITCZ excursion in both hemispheres.
    !
    ! i_mmc_geo = 1  Both Hadley edges sit at +-hw + c_mmc_edge*ficz, with hw = pi/(6*hadwidsc) half the
    !                Hadley width, and the phase is linear in latitude within each Hadley cell, so the cell
    !                boundary stays at ficz and the summer cell narrows.  The Ferrel and polar cells keep the
    !                width hw poleward of the edges.
    if (i_mmc_geo.eq.1) then
      hw  = pi/(6._wp*hadwidsc)
      fen = max(hw + c_mmc_edge*ficz, ficz + dfi_min)
      fes = min(-hw + c_mmc_edge*ficz, ficz - dfi_min)
      if (fi.ge.fen) then
        ff = pi + 6._wp*hadwidsc*(fi-fen)
      else if (fi.ge.ficz) then
        ff = pi*(fi-ficz)/(fen-ficz)
      else if (fi.ge.fes) then
        ff = pi*(fi-ficz)/(ficz-fes)
      else
        ff = -pi + 6._wp*hadwidsc*(fi-fes)
      endif
    else
      ff = 6._wp*hadwidsc*(fi-ficz/(c_mmc_1*(fi-ficz)**2+1._wp))
    endif

  end function cell_phase

end module slp_mod
