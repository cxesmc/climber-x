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
  use atm_params, only : c_slp_1, c_slp_2, c_slp_3, c_slp_4, c_slp_5
  use atm_params, only : l_aslp_temp_adv, c_aslp_temp_tau
  use atm_params, only : l_aslp_topo, c_aslp_topo_1, c_aslp_topo_2, c_aslp_topo_3, c_aslp_topo_4
  use atm_params, only : cp
  use atm_params, only : c_mmc_had, c_mmc_fer, c_mmc_pol, c_mmc_1, c_mmc_2
  use atm_params, only : i_mmc_itcz, c_mmc_efe
  use atm_params, only : i_mmc_fer
  use atm_params, only : i_mmc_had, c_mmc_q, q_mmc_ref
  use atm_params, only : c_mmc_dt0, c_mmc_dt1, hdv_mmc_ref
  use atm_params, only : nsmooth_aslp, nsmooth_aslp_topo
  use atm_grid, only : im, jm, jmc, aim, jeq, jts, jtn, jps, jpn, dy, pl, k500, i_ice
  use atm_grid, only : fcorua, sint, cost, fiu, fit
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
  subroutine azslp(tsksl, htrop, zsa, uz500, &
      aslp, &
      aslp_temp, aslp_topo, dz500, atsl)

    use, intrinsic :: iso_c_binding 

    implicit none

    ! include FFTW3 library for Fourier Transform (https://www.fftw.org/)
    include 'fftw3.f03'

    real(wp), intent(in   ) :: tsksl(:,:)
    real(wp), intent(in   ) :: htrop(:,:)
    real(wp), intent(in   ) :: zsa(:,:)
    real(wp), intent(in   ) :: uz500(:)

    real(wp), intent(inout) :: aslp(:,:)
    real(wp), intent(inout) :: aslp_temp(:,:)
    real(wp), intent(inout) :: aslp_topo(:,:)

    real(wp), intent(out  ) :: atsl(:,:)
    real(wp), intent(out  ) :: dz500(:,:)

    integer :: i, j
    real(wp) :: cor, beta, k
    real(wp) :: m, r, uz
    real(wp) :: tslz(jm)
    real(wp) :: htropz(jm)
    real(wp) :: uz500s(jm)
    real(wp) :: u500(jm)
    real(wp) :: dz500o(im)
    type(C_PTR) :: plan_r2c, plan_c2r
    type(C_PTR) :: plan_r2c_temp, plan_c2r_temp
    real(wp), dimension(im) :: eps, Kn2
    real(wp) :: zsa_smooth(im,jm)
    real(dp) :: zsa_smooth_dp(im,jm)
    real(dp), dimension(im) :: psi
    complex(dp), dimension(im/2+1) :: zsa_fft
    complex(dp), dimension(im/2+1) :: psi_fft
    real(dp), dimension(im) :: aslp_temp_dp
    complex(dp), dimension(im/2+1) :: aslp_temp_fft


    ! smooth zonal mean 500 hPa zonal wind
    do j=2,jm-1
      uz500s(j) = 0.4_wp*uz500(j)+0.3_wp*(uz500(j+1)+uz500(j-1))
    enddo
    uz500s(1)  = 0.5_wp*(uz500(1)+uz500(2))
    uz500s(jm) = 0.5_wp*(uz500(jm)+uz500(jm-1))

    do j=1,jm
      u500(j) = max(0.1_wp,uz500s(j)) * c_aslp_topo_3
    enddo

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
          !aslp_temp(i,j) = -c_slp_1*g*p0*12000._wp/(2._wp*Rd*T0**2)*atsl(i,j)
          aslp_temp(i,j) = -c_slp_1*g*p0*htropz(j)/(2._wp*Rd*T0**2)*atsl(i,j)
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

      ! smooth topography
      do j=1,jm
        zsa_smooth(:,j) = zsa(:,j) * c_aslp_topo_4      ! optional scaling to mimic difference between surface wind and 500 hPa wind
      enddo
      call smooth2_m(zsa_smooth,nsmooth_aslp_topo)

      zsa_smooth_dp = zsa_smooth

      ! Make forward and backward plans for the FFT
      plan_r2c = fftw_plan_dft_r2c_1d(im, zsa_smooth_dp, zsa_fft, FFTW_ESTIMATE)
      plan_c2r = fftw_plan_dft_c2r_1d(im, psi_fft, psi, FFTW_ESTIMATE)

      do j=1,jm

        cor  = 2._wp*omega*sint(j)  
        beta = 2._wp*omega*cost(j)/r_earth  
        k = 2._wp*pi/(2._wp*pi*r_earth*cost(j))      ! lowest zonal wavenumber 

        uz = u500(j)
        r = c_aslp_topo_1
        m = pi/c_aslp_topo_2

        do i=1,im
          Kn2(i) = (k*(i-1))**2+m**2
        enddo
        do i=2,im
          eps(i) = r*Kn2(i)/(k*(i-1)*uz)
        enddo
        eps(1) = 0._wp

        !  forward transform the data
        call fftw_execute_dft_r2c(plan_r2c, zsa_smooth_dp(:,j), zsa_fft)

        do i=1,im/2+1
          psi_fft(i) = cor*zsa_fft(i)/(htropz(j)*(Kn2(i)-beta/uz-cmplx(0._wp,eps(i)) ))
        enddo

        ! backward transform the data
        call fftw_execute_dft_c2r(plan_c2r, psi_fft, psi )

        ! zonal anomalies
        psi = psi - sum(psi)*aim
        ! convert to geopotential height perturbation
        dz500o = psi * abs(cor)/g * aim  ! m

        ! convert from geopotential height anomalies at ~500 hPa (dz500o) to sea level pressure anomalies 
        ! assuming dz is constant throughout the lower troposphere, see PIK-report eq. 7.23, 7.24, 4.73  
        aslp_topo(:,j) = dz500o(:) * ra*pl(k500)*g    ! Pa

      enddo

      ! destrox FFT plans
      call fftw_destroy_plan(plan_r2c)
      call fftw_destroy_plan(plan_c2r)

    else

      aslp_topo(:,:) = 0._wp

    endif


    !-------------------------------------------------------------------------------
    ! geopotential height zonal anomalies at 500 hPa from planetary stationary waves
    !-------------------------------------------------------------------------------

    dz500(:,:) = aslp_topo(:,:)/(ra*pl(k500)*g)  ! m


    !------------------------------------------------
    ! azonal component of sea level pressure 
    !------------------------------------------------

    do j=1,jm
      do i=1,im
        ! relax in time
        aslp(i,j) = (1._wp-c_slp_2)*aslp(i,j)+c_slp_2*(aslp_temp(i,j)+aslp_topo(i,j)) 
      enddo
    enddo

    ! smooth in space
    call smooth2_m(aslp,nsmooth_aslp)

    ! polar and equatorial damping
    do j=1,jm
      do i=1,im
        aslp(i,j)=aslp(i,j) &
          * min(1._wp,c_slp_3*(cost(j)-cost(1))**2) * min(1._wp,c_slp_5+c_slp_4*(sint(j)**2-sint(jm/2)**2))
      enddo
    enddo


    return

  end subroutine azslp


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  zslp
  !   Purpose    :  compute zonally averaged sea level pressure component
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine zslp(sin_cos_acbar, tsl, aslp, fdydse, t3, q3, htrop, ttrop, efe, &
          ficz, slp, had_fi, had_width)

    implicit none

    real(wp), intent(in ) :: sin_cos_acbar(:,:)
    real(wp), intent(in ) :: tsl(:,:)
    real(wp), intent(in ) :: aslp(:,:)
    real(wp), intent(in ) :: fdydse(:,:)
    real(wp), intent(in ) :: t3(:,:,:)
    real(wp), intent(in ) :: q3(:,:,:)
    real(wp), intent(in ) :: htrop(:,:)
    real(wp), intent(in ) :: ttrop(:,:)
    real(wp), intent(in ) :: efe          !! energy flux equator from the previous step (deg), see efe_diag

    real(wp), intent(out) :: ficz         !! ITCZ latitude of the cell geometry (radians)
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
      fisu(j) = 6._wp*hadwidsc*(fiu(j)-ficz/(c_mmc_1*(fiu(j)-ficz)**2+1._wp))
      fist(j) = 6._wp*hadwidsc*(fit(j)-ficz/(c_mmc_1*(fit(j)-ficz)**2+1._wp))
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

      ff = 6._wp*hadwidsc*(fiu(j)-ficz/(c_mmc_1*(fiu(j)-ficz)**2+1._wp))

      coc = 0._wp

      if (ff.ge.0._wp .and. ff.lt.pi)              coc = c_mmc_had*dthn
      if (ff.ge.pi .and. ff.lt.2._wp*pi)           coc = cocfn
      if (ff.ge.2._wp*pi .and. ff.lt.3._wp*pi)     coc = c_mmc_pol*dtpn

      if (-ff.gt.0._wp .and. (-ff).le.pi)          coc = c_mmc_had*dths
      if (-ff.gt.pi .and. (-ff).le.2._wp*pi)       coc = cocfs
      if (-ff.gt.2._wp*pi .and. (-ff).le.3._wp*pi) coc = c_mmc_pol*dtps

      vsz(j) = -coc*sin(ff) 

    enddo

    ! Zonally averaged SLP   
    ! Integration from NP to SP     

    psz(:) = 0._wp
    do j=2,jm
      psz(j) = psz(j-1) + vsz(j)*fcorua(j)*ra/acbarz(j) * dy
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
      if (abs(la).gt.fi_efe_max .or. abs(lb).gt.fi_efe_max) cycle
      if ((f(j).gt.0._dp .and. f(j+1).le.0._dp) .or. (f(j).lt.0._dp .and. f(j+1).ge.0._dp)) then
        x = lb + (la-lb)*real(-f(j+1)/(f(j)-f(j+1)),wp)
        if (efe.lt.-900._wp .or. abs(x).lt.abs(efe)) efe = x
      endif
    enddo
    ! no crossing in the window: the EFE lies beyond it, on the side opposite to the transport
    if (efe.lt.-900._wp) efe = -sign(fi_efe_max, aht_eq)

    return

  end subroutine efe_diag

end module slp_mod
