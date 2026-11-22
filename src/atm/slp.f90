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
  use atm_params, only : ra, p0, hatm, pbl
  use atm_params, only : c_slp_1, c_slp_2, c_aslp_eq_hb, c_aslp_eq_lat0
  use atm_params, only : c_slp_clim_f, c_slp_clim_tau, n_slp_clim_harm
  use atm_params, only : l_aslp_temp_adv, c_aslp_temp_tau, c_aslp_temp_dT, c_aslp_temp_wtau
  use atm_params, only : l_aslp_topo, c_aslp_sw_amp, c_aslp_sw_tau, c_aslp_sw_h, c_aslp_sw_umin, c_aslp_sw_lateq, n_aslp_sw
  use atm_params, only : c_aslp_sw_fr, c_aslp_sw_nbv
  use timer, only : nday_year
  use atm_params, only : cp
  use atm_params, only : c_mmc_had, c_mmc_fer, c_mmc_pol, c_mmc_1, c_mmc_2, i_mmc_geo, c_mmc_edge
  use atm_params, only : i_mmc_itcz, c_mmc_efe
  use atm_params, only : i_mmc_fer, i_mmc_fer_shape, c_mmc_fer_q
  use atm_params, only : i_mmc_had, c_mmc_q, q_mmc_ref
  use atm_params, only : c_mmc_dt0, c_mmc_dt1, hdv_mmc_ref, c_mmc_mt, nsmooth_mtz, amas
  use atm_params, only : aslp_polar_lat, n_aslp_shapiro
  use atm_grid, only : im, jm, jmc, aim, jeq, jts, jtn, jps, jpn, dy, pl, k500, k850, k300, i_ice
  use atm_grid, only : fcorua, sint, cost, fiu, fit, dxt
  use atm_grid, only : z3, k_dse_lo, k_dse_up
  use smooth_atm_mod, only : zona
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
  subroutine azslp(doy, tlay, htrop, zsa, uz500, u3f, frlnd, tsl, &
      aslp, &
      aslp_temp, aslp_topo, &
      aslp_temp_hcoef, aslp_temp_hacc, aslp_temp_hnday, aslp_temp_hnyr, aslp_temp_lc)

    use, intrinsic :: iso_c_binding 

    implicit none

    ! include FFTW3 library for Fourier Transform (https://www.fftw.org/)
    include 'fftw3.f03'

    integer,  intent(in   ) :: doy            !! day of the year, 1..nday_year (azslp is called once per day)
    real(wp), intent(in   ) :: tlay(:,:)
    real(wp), intent(in   ) :: htrop(:,:)
    real(wp), intent(in   ) :: zsa(:,:)
    real(wp), intent(in   ) :: uz500(:)       !! zonal mean 500 hPa zonal wind for the upstream displacement of aslp_temp, thermal wind damped with c_uz500_pol
    real(wp), intent(in   ) :: u3f(:,:,:)     !! 3D zonal wind with undamped thermal wind, previous step, used if l_aslp_topo
    real(wp), intent(in   ) :: frlnd(:,:)     !! land fraction, used only if c_aslp_temp_dT>0
    real(wp), intent(in   ) :: tsl(:,:)       !! air temperature profile reduced to sea level (K), previous step, used only if c_aslp_temp_dT>0

    real(wp), intent(inout) :: aslp(:,:)
    real(wp), intent(inout) :: aslp_temp(:,:)
    real(wp), intent(inout) :: aslp_topo(:,:)
    real(wp), intent(inout) :: aslp_temp_hcoef(:,:,:) !! harmonic coefficients of the aslp_temp seasonal climatology, used only if c_slp_clim_f<1
    real(wp), intent(inout) :: aslp_temp_hacc(:,:,:)  !! projection of the current year's aslp_temp on the harmonics
    real(wp), intent(inout) :: aslp_temp_hnday        !! days accumulated in aslp_temp_hacc
    real(wp), intent(inout) :: aslp_temp_hnyr         !! complete years entered into aslp_temp_hcoef
    real(wp), intent(inout) :: aslp_temp_lc(:)        !! K, time-relaxed land-mean azonal sea level temperature per row (9999 = not initialised)


    integer :: i, j
    real(wp) :: k, uz
    real(wp) :: tlayz(jm)
    real(wp) :: rho_sw(jm)
    real(wp) :: htropz(jm)
    real(wp) :: uz500s(jm)
    real(wp) :: tau_j(jm), clnd(jm), sl, sw, wr, cw
    real(wp) :: atsl(im,jm)          !! K, azonal air temperature reduced to sea level (tsl), for the land contrast
    real(wp) :: atlay(im,jm)         !! K, azonal sea level-500 hPa layer-mean temperature, drives aslp_temp
    integer :: jj
    type(C_PTR) :: plan_r2c_temp, plan_c2r_temp
    real(wp) :: aslp_frc(im,jm)
    real(wp) :: weq
    integer :: h, nharm
    real(wp) :: phase, wyr
    real(wp) :: basis(2*n_slp_clim_harm+1), wproj(2*n_slp_clim_harm+1)
    real(wp) :: aslp_temp_clim(im,jm)
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
      tlayz(j) = 0._wp      
      htropz(j) = 0._wp      
      do i=1,im
        tlayz(j) = tlayz(j) + tlay(i,j)*aim 
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
        atlay(i,j) = tlay(i,j)-tlayz(j)
      enddo

      do i=1,im
        if (j.eq.1 .or. j.eq.jm) then
          ! azonal SLP vanishes at the Poles
          aslp_temp(i,j) = 0._wp
        else
          ! hypsometric gain of the sea level-500 hPa layer, linearised about the zonal mean layer temperature:
          ! with z500 = Rd*Tm/g*ln(slp/p500) held fixed (no azonal pressure at 500 hPa),
          ! d slp/d Tm = -slp*ln(slp/p500)/Tm = -p0*ln(2)/Tm (p500 = p0/2, as in the tlay layer; 256 Pa/K at 273 K).
          ! tlayz is the zonal mean of tlay, i.e. the zonal mean layer temperature. (The fixed CLIMBER-2 gain,
          ! Petoukhov 2000 eq. 17, g*p0*5000 m/(Rd*T0**2) = 232 Pa/K, was used before.)
          aslp_temp(i,j) = -c_slp_1*p0*log(2._wp)/tlayz(j)*atlay(i,j)
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
    ! centred upstream of it. The temperature anomaly therefore lags the pressure anomaly: in ERA-Interim the
    ! DJF azonal sea level pressure correlates best with the azonal sea level temperature
    ! displaced 15-20 deg to the east (r=-0.83 over 40-75N, r=-0.93 over North America,
    ! against only -0.70 and -0.54 at zero lag).
    ! The equilibrium thermal pressure anomaly is therefore displaced upstream, as the
    ! steady state of
    !   -uz*d(aslp_temp)/dx = (aslp_temp_eq-aslp_temp)/c_aslp_temp_tau ,
    ! solved in Fourier space, where each zonal wavenumber n is shifted upstream by the
    ! phase angle atan(k*n*uz*c_aslp_temp_tau) and damped by 1/sqrt(1+(k*n*uz*tau)**2).

    ! Weighting of the shift by the land temperature contrast (c_aslp_temp_dT>0). The upstream offset is a
    ! cold-continent feature: in ERA-Interim (40-75N, analysis/stationary_wave/era_shift_contrast.jl) the SLP is
    ! 10-20 deg upstream of the azonal temperature from October to March, when the continents are colder
    ! than the zonal mean, and in phase (or slightly downstream) from May to August, when they are warmer and
    ! the thermal lows sit on their heating. Per row, tau = c_aslp_temp_tau * min(1, max(0, -C/c_aslp_temp_dT)),
    ! with C the land-fraction weighted mean of the azonal sea level temperature atsl over the row and its two
    ! neighbours (weights 1/4,1/2,1/4),
    ! relaxed in time with c_aslp_temp_wtau (days, azslp is called once per day) so that the weight does not
    ! follow synoptic atsl and adds no fast feedback. Rows (and neighbours) with less than 5% land get no shift.
    ! With c_aslp_temp_dT=0 the shift uses the constant c_aslp_temp_tau everywhere.
    ! The contrast is taken from atsl, the azonal air temperature profile reduced to sea level (40-75N PI land
    ! DJF -3.5, JJA +2.7 K, close to the ERA t2m behaviour), not from the azonal layer mean atlay that is shifted:
    ! its land contrast is weak in the model (DJF -0.8, JJA +0.5 K) and has no seasonal cycle at 42.5S (tried
    ! before as i_aslp_temp_lc=1; both score the same at PI and LGM, output/lcshift2).
    if (l_aslp_temp_adv) then
      do j=1,jm
        atsl(:,j) = tsl(:,j) - sum(tsl(:,j))*aim
      enddo
      do j=1,jm
        sl = 0._wp
        sw = 0._wp
        do jj=max(1,j-1),min(jm,j+1)
          wr = 0.25_wp
          if (jj.eq.j) wr = 0.5_wp
          sl = sl + wr*sum(frlnd(:,jj)*atsl(:,jj))
          sw = sw + wr*sum(frlnd(:,jj))
        enddo
        if (sw.gt.0.05_wp*real(im,wp)) then
          clnd(j) = sl/sw
        else
          clnd(j) = 0._wp
        endif
      enddo
      cw = min(1._wp, 1._wp/max(c_aslp_temp_wtau,1._wp))
      do j=1,jm
        if (aslp_temp_lc(j).gt.9000._wp) then
          aslp_temp_lc(j) = clnd(j)
        else
          aslp_temp_lc(j) = aslp_temp_lc(j) + cw*(clnd(j)-aslp_temp_lc(j))
        endif
        if (c_aslp_temp_dT.gt.0._wp) then
          tau_j(j) = c_aslp_temp_tau*min(1._wp, max(0._wp, -aslp_temp_lc(j)/c_aslp_temp_dT))
        else
          tau_j(j) = c_aslp_temp_tau
        endif
      enddo
    endif

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
            / cmplx(1._dp, -real(k*(i-1)*uz*tau_j(j),dp), dp)
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

      ! density converting the stationary-wave height anomaly to sea level pressure. The azonal SLP is the linearised hypsometric split of the sea level-500 hPa layer,
      !   a(slp) = p0*g/(Rd*Tm)*a(Z500) - p0*ln(2)/Tm*a(Tm) ,
      ! aslp_temp being the second term; the stationary-wave height is then read as a(Z500) and converted with
      ! the same zonal mean layer temperature (tlayz = Tm) instead of the fixed density ra at T0
      ! (ERA-Interim: p0/(Rd*Tm) = 0.98-1.10 ra, analysis/itsl0_rossby/aslp_topo_hyps.jl).
      ! The propagation wind stays at 300 hPa (u3f at k300), which gives the best pattern
      ! (analysis/aslp_topo/level_test.jl: DJF NH r 0.64 vs 0.40 with the 500 hPa wind).
      do j=1,jm
        rho_sw(j) = p0/(Rd*tlayz(j))
      enddo

      ! linear barotropic stationary waves, coupled in latitude
      call stationary_waves(zsa, u3f, rho_sw, aslp_topo)

    else

      aslp_topo(:,:) = 0._wp

    endif


    !------------------------------------------------
    ! damping of the departure of the thermal azonal SLP from its own seasonal climatology
    !------------------------------------------------

    ! The thermal-low feedback (aslp_temp -> geostrophic/ageostrophic wind -> advection, convergence
    ! -> atlay) produces spurious travelling modes at every sub-seasonal period: the ~10 d W African
    ! oscillator, the 16 d k=2 and 40-50 d k=1 extratropical waves, the 20-40 d Central Asian and the
    ! LGM ~10 d Tibet modes. A relaxed slow state (aslp = slow + f*(aslp-slow), tried before) separates
    ! them only for periods << tau, lags the seasonal cycle by (1-f)*tau and makes the 40-50 d wave worse.
    ! Here the slow part is the model's own seasonal climatology of aslp_temp, a mean plus
    ! n_slp_clim_harm annual harmonics, fitted to each complete year and averaged over the years
    ! (cumulative mean for the first c_slp_clim_tau years, then e-folding c_slp_clim_tau years).
    ! Only the departure from it is damped, aslp_temp = clim + c_slp_clim_f*(aslp_temp-clim): no
    ! seasonal lag at any f, the mean and seasonal cycle are kept, and every sub-seasonal and
    ! interannual anomaly loses loop gain by the same factor, at all latitudes. The cost is that
    ! interannual anomalies and the response to a forcing change are damped/lag by ~c_slp_clim_tau.
    ! The filter starts after the first complete year (hnyr>=1).
    if (c_slp_clim_f.lt.1._wp) then
      nharm = 2*n_slp_clim_harm+1
      phase = 2._wp*pi*(real(doy,wp)-0.5_wp)/real(nday_year,wp)
      basis(1) = 1._wp
      wproj(1) = 1._wp/real(nday_year,wp)
      do h=1,n_slp_clim_harm
        basis(2*h)   = cos(real(h,wp)*phase)
        basis(2*h+1) = sin(real(h,wp)*phase)
        wproj(2*h)   = 2._wp/real(nday_year,wp)
        wproj(2*h+1) = 2._wp/real(nday_year,wp)
      enddo
      ! accumulate the unfiltered aslp_temp of this year
      if (doy.eq.1) then
        aslp_temp_hacc = 0._wp
        aslp_temp_hnday = 0._wp
      endif
      do h=1,nharm
        aslp_temp_hacc(:,:,h) = aslp_temp_hacc(:,:,h) + wproj(h)*basis(h)*aslp_temp(:,:)
      enddo
      aslp_temp_hnday = aslp_temp_hnday + 1._wp
      ! damp the departure from the climatology of the previous years
      if (aslp_temp_hnyr.ge.1._wp) then
        aslp_temp_clim = 0._wp
        do h=1,nharm
          aslp_temp_clim(:,:) = aslp_temp_clim(:,:) + basis(h)*aslp_temp_hcoef(:,:,h)
        enddo
        aslp_temp(:,:) = aslp_temp_clim(:,:) + c_slp_clim_f*(aslp_temp(:,:)-aslp_temp_clim(:,:))
      endif
      ! end of year: enter the fit of a complete year into the climatology
      if (doy.eq.nday_year) then
        if (nint(aslp_temp_hnday).eq.nday_year) then
          wyr = 1._wp/min(aslp_temp_hnyr+1._wp, max(1._wp,c_slp_clim_tau))
          aslp_temp_hcoef = aslp_temp_hcoef + wyr*(aslp_temp_hacc-aslp_temp_hcoef)
          aslp_temp_hnyr = aslp_temp_hnyr + 1._wp
        endif
        aslp_temp_hacc = 0._wp
        aslp_temp_hnday = 0._wp
      endif
    endif


    !------------------------------------------------
    ! azonal component of sea level pressure 
    !------------------------------------------------

    ! the forcing is tapered once per call through its vorticity (does not accumulate), while smoothing and
    ! equatorial damping act on the time-relaxed state at every call, so they accumulate: the equilibrium
    ! is c S X/(1-(1-c) S) instead of S X (analysis/stationary_wave/aslp_smooth_loop.jl)
    do j=1,jm
      do i=1,im
        aslp_frc(i,j) = aslp_temp(i,j)+aslp_topo(i,j)
      enddo
    enddo
    ! smooth equatorial taper of the TOTAL forcing, 1-(1-c_aslp_eq_hb)*exp(-(lat/c_aslp_eq_lat0)^2) (c_aslp_eq_hb=1: none).
    ! It is applied before the relaxation, so it does not accumulate and leaves the relaxation time latitude-independent.
    ! A damping of the relaxed state at every call (tried before, min(1,c5+c4*(sin2-sin2eq))) accumulates to
    ! c_slp_2*S/(1-(1-c_slp_2)*S), a factor-6 cliff in one row near 15 deg, shortens the equatorial response time to ~1 d
    ! and crashed (equatorial aslp/wind blow-up in year 1) for c_slp_2*c5 >= 0.1.
    do j=1,jm
      weq = 1._wp - (1._wp-c_aslp_eq_hb)*exp(-(asin(sint(j))*180._wp/pi/c_aslp_eq_lat0)**2)
      aslp_frc(:,j) = weq*aslp_frc(:,j)
    enddo
    ! Shapiro filter of the forcing, once per call (does not accumulate, independent of c_slp_2)
    call shapiro_filter(aslp_frc, n_aslp_shapiro)
    call vorticity_taper(aslp_frc)
    do j=1,jm
      do i=1,im
        ! relax in time
        aslp(i,j) = (1._wp-c_slp_2)*aslp(i,j)+c_slp_2*aslp_frc(i,j)
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
  !   Subroutine :  s h a p i r o _ f i l t e r
  !   Purpose    :  Shapiro low-pass filter of order n of an azonal field
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! In each direction x <- x - (-delta2/4)**n x, with delta2 the 3-point second difference (periodic in
  ! longitude, mirror rows at the poles), applied first in longitude and then in latitude (Shapiro 1970).
  ! The response per direction is 1 - sin(k*dx/2)**(2n): exactly zero at 2 dx, monotonic, never negative,
  ! and flat at large scales (n=2: 0.94 at 6 dx, 0.98 at 8 dx; n=4: 0.996 at 6 dx). It is applied to the
  ! forcing once per call, so it does not accumulate through c_slp_2 (a 5-point smoothing of the relaxed
  ! state, used before, did: equilibrium c S/(1-(1-c) S), ~40% weaker and with a -0.2 checkerboard response).
  ! Both directions are linear and zonally invariant, so a zero zonal mean stays zero.
  subroutine shapiro_filter(x, n)

    implicit none

    real(wp), intent(inout) :: x(:,:)
    integer, intent(in) :: n

    integer :: i, j, k, ni, nj, im1, ip1, jm1, jp1
    real(wp), allocatable :: d(:,:), e(:,:)

    if (n.le.0) return
    ni = size(x,1)
    nj = size(x,2)
    allocate(d(ni,nj), e(ni,nj))

    ! longitude
    d = x
    do k=1,n
      do j=1,nj
        do i=1,ni
          im1 = i-1
          if (im1.lt.1) im1 = ni
          ip1 = i+1
          if (ip1.gt.ni) ip1 = 1
          e(i,j) = -0.25_wp*(d(ip1,j)-2._wp*d(i,j)+d(im1,j))
        enddo
      enddo
      d = e
    enddo
    x = x - d

    ! latitude
    d = x
    do k=1,n
      do j=1,nj
        jm1 = max(1,j-1)
        jp1 = min(nj,j+1)
        do i=1,ni
          e(i,j) = -0.25_wp*(d(i,jp1)-2._wp*d(i,j)+d(i,jm1))
        enddo
      enddo
      d = e
    enddo
    x = x - d

    deallocate(d, e)

    return

  end subroutine shapiro_filter


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  stationary_waves
  !   Purpose    :  sea level pressure of the topographically forced stationary
  !                 planetary waves from a linear barotropic model
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine stationary_waves(zsa, u3f, rho_sw, aslp_sw)

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
    ! barotropic structure, aslp_sw = c_aslp_sw_amp * rho_sw*g * f*psi/g, with rho_sw = ra, or
    ! p0/(Rd*Tm) with the zonal mean layer temperature (psi read as Z500).
    !
    ! The azonal latent heating is not used as a barotropic Rossby wave source: the equivalent barotropic
    ! conversion maps the upper-level anticyclone over the heating to the surface, i.e. a high where the
    ! observed surface low is. A Gill-type baroclinic response to it (tried before) moved the equatorial low
    ! to the central Pacific and removed the E Pacific high, collapsing the coupled Pacific SST gradient.

    use, intrinsic :: iso_c_binding

    implicit none

    include 'fftw3.f03'

    real(wp), intent(in ) :: zsa(:,:)
    real(wp), intent(in ) :: u3f(:,:,:)
    real(wp), intent(in ) :: rho_sw(:)      !! kg/m3, zonal mean density converting the height anomaly to sea level pressure
    real(wp), intent(out) :: aslp_sw(:,:)

    integer :: i, j, n, nmax
    real(wp) :: dphi, ri
    real(wp) :: phi(jm), uup(jm), ulo(jm), ue(jm), zbar(jm), qy(jm)
    real(wp) :: cph(0:jm)            ! cos at the half latitudes between j and j+1 (0 at the poles)
    real(dp) :: frc(im,jm)
    real(dp), dimension(im) :: rowr
    complex(dp), dimension(im/2+1) :: rowc
    complex(dp) :: frc_hat(im/2+1,jm), psi_hat(im/2+1,jm)
    real(wp) :: hfr
    real(wp) :: zeff(im,jm)
    complex(dp) :: al(jm), ad(jm), au(jm), rhs(jm), x(jm)
    type(C_PTR) :: plan_r2c, plan_c2r

    nmax = min(n_aslp_sw, im/2-1)
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

    plan_r2c = fftw_plan_dft_r2c_1d(im, rowr, rowc, FFTW_ESTIMATE)
    plan_c2r = fftw_plan_dft_c2r_1d(im, rowc, rowr, FFTW_ESTIMATE)

    do j=1,jm
      rowr(:) = frc(:,j)
      call fftw_execute_dft_r2c(plan_r2c, rowr, rowc)
      frc_hat(:,j) = rowc(:)
    enddo

    psi_hat(:,:) = (0._dp,0._dp)
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
    enddo

    do j=1,jm
      rowc(:) = psi_hat(:,j)
      call fftw_execute_dft_c2r(plan_c2r, rowc, rowr)
      if (abs(phi(j)).lt.c_aslp_sw_lateq*pi/180._wp) then
        aslp_sw(:,j) = 0._wp
      else
        ! f*psi/g is the height anomaly, aim the FFTW normalisation
        aslp_sw(:,j) = c_aslp_sw_amp*rho_sw(j)*g * 2._wp*omega*sint(j)/g * rowr(:)*aim
      endif
    enddo

    call fftw_destroy_plan(plan_r2c)
    call fftw_destroy_plan(plan_c2r)

    return

  end subroutine stationary_waves


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
    real(wp) :: aferzn, aferzs
    real(wp) :: qtrp, wtrp, fqhad
    real(wp) :: htrpz, sbotz, stopz, tsfcz, wwid, hdv

    real(wp) :: fisu(jmc), fist(jm), vsz(jmc), tslz(jm), psz(jmc), acbarz(jm)
    real(wp) :: dsez(jm), fedz(jmc), aferz(jm)
    real(wp) :: mtu, sgnf
    real(wp) :: mtz0(jm)
    integer :: ks


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
        sbotz = sbotz + (t3(i,j,1)+g*z3(i,j,1)/cp)*cost(j)
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
      ! with the geometric distance of the two layer centres of each column (vesta)
      do j=1,jm
        dsez(j) = 0._wp
        do i=1,im
          dsez(j) = dsez(j) + (cp*(t3(i,j,k_dse_up)-t3(i,j,k_dse_lo)) + g*(z3(i,j,k_dse_up)-z3(i,j,k_dse_lo)))*aim
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

      ! i_mmc_fer_shape = 2: Ferrel branches from the latitude-resolved TEM closure instead of the half-sine with the
      ! belt-mean amplitude. aferz(j) = c_mmc_fer*|F_eddy|/Ds is the cell's mass flux at this v-point, the velocity of the
      ! branch scales with 1/cos(lat) (reference 45 deg, absorbed in c_mmc_fer), and |sin(ff)|^c_mmc_fer_q tapers it to
      ! zero at the cell edges only. F/Ds is much flatter across the cell than sin(ff), which moves pressure gradient
      ! onto the poleward flank of the westerlies (offline emulator: analysis/ferrel_fds_shape_emulator.jl).
      if (i_mmc_fer_shape.eq.2) then
        if ((ff.ge.pi .and. ff.lt.2._wp*pi) .or. (-ff.gt.pi .and. (-ff).le.2._wp*pi)) then
          vsz(j) = -sign(1._wp,sin(ff)) * aferz(j) * cos(0.25_wp*pi)/max(cos(fiu(j)),0.1_wp) * abs(sin(ff))**c_mmc_fer_q
        endif
      endif

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
    ! The friction (Ekman) part of the branch is vsz - MT/(f*M_pbl), with M_pbl = (1-pbl)*amas the mass of
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
    ! Meridional 1-2-1 smoothing (nsmooth_mtz passes; the pole rows stay zero). On the 5 deg grid single rows pick
    ! up individual ranges (Tibet, Rockies, Greenland), and the torque, a pressure x slope covariance, is a small
    ! residual that is very sensitive to the phase of aslp: its row-to-row scatter (rms 2nd difference / rms over
    ! 27.5-72.5N 0.57 DJF, 0.39 JJA, CMIP 0.25-0.67) goes straight into the zonal SLP gradient. One pass halves it
    ! and leaves the 10 deg band means nearly unchanged (analysis/mtz_robust.jl).
    do ks=1,nsmooth_mtz
      mtz0 = mtz
      do j=2,jm-1
        mtz(j) = 0.25_wp*mtz0(j-1) + 0.5_wp*mtz0(j) + 0.25_wp*mtz0(j+1)
      enddo
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
      psz(j) = psz(j-1) + (vsz(j)*fcorua(j) - c_mmc_mt*sgnf*mtu/((1._wp-pbl)*amas))*ra/acbarz(j) * dy
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
    real(wp), parameter :: dfer_min = 10._wp*pi/180._wp !! minimum latitudinal extent of a Ferrel cell (i_mmc_geo=2)
    real(wp) :: hw, fen, fes, fpn, fps, fmn, fms

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
    !
    ! i_mmc_geo = 2  Hadley edges as in i_mmc_geo = 1.  The Ferrel/polar boundaries sit at fpn/fps = +-2*hw,
    !                independent of ficz: they follow FORCED changes of the Hadley width (Held scaling) but not the
    !                seasonal displacement of the Hadley edges.  In CMIP6/PMIP the polar edge stays at ~61-65 deg all
    !                year while the Hadley edge moves 7-9 deg (the winter Ferrel cell is the wide one), and the
    !                annual-mean polar edge is ~2x the Hadley edge in every state (|P| = 3.8 + 1.94|H|, rms 3.6 deg over
    !                aquaplanet, PI, 4xCO2 and LGM; ratio 1.96-2.24).  Within each Ferrel cell the phase is piecewise
    !                linear: 1.5 pi (maximum of the branch) is reached at hw/2 from the Hadley edge (the mid-point of
    !                the cell for ficz = 0), kept within the central half of the cell; the polar cells keep the
    !                angular rate of i_mmc_geo = 1 (pi per hw).
    if (i_mmc_geo.eq.2) then
      hw  = pi/(6._wp*hadwidsc)
      fen = max(hw + c_mmc_edge*ficz, ficz + dfi_min)
      fes = min(-hw + c_mmc_edge*ficz, ficz - dfi_min)
      fpn = max( 2._wp*hw, fen + dfer_min)
      fps = min(-2._wp*hw, fes - dfer_min)
      ! mid-point of the Ferrel cell (phase 1.5 pi, maximum of the branch) at hw/2 from the Hadley edge,
      ! kept within the central half of the cell
      fmn = fen + min(max(0.5_wp*hw, 0.25_wp*(fpn-fen)), 0.75_wp*(fpn-fen))
      fms = fes - min(max(0.5_wp*hw, 0.25_wp*(fes-fps)), 0.75_wp*(fes-fps))
      if (fi.ge.fpn) then
        ff = 2._wp*pi + pi*(fi-fpn)/hw
      else if (fi.ge.fmn) then
        ff = 1.5_wp*pi + 0.5_wp*pi*(fi-fmn)/(fpn-fmn)
      else if (fi.ge.fen) then
        ff = pi + 0.5_wp*pi*(fi-fen)/(fmn-fen)
      else if (fi.ge.ficz) then
        ff = pi*(fi-ficz)/(fen-ficz)
      else if (fi.ge.fes) then
        ff = pi*(fi-ficz)/(ficz-fes)
      else if (fi.ge.fms) then
        ff = -pi - 0.5_wp*pi*(fes-fi)/(fes-fms)
      else if (fi.ge.fps) then
        ff = -1.5_wp*pi - 0.5_wp*pi*(fms-fi)/(fms-fps)
      else
        ff = -2._wp*pi + pi*(fi-fps)/hw
      endif
    else if (i_mmc_geo.eq.1) then
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
