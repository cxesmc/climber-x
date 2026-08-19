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
  use constants, only : pi, r_earth, omega, g, Rd, T0
  use atm_params, only : ra, p0, hatm
  use atm_params, only : c_slp_1, c_slp_2, c_slp_3, c_slp_4, c_slp_5
  use atm_params, only : l_aslp_temp_adv, c_aslp_temp_tau
  use atm_params, only : l_aslp_topo, c_aslp_topo_1, c_aslp_topo_2, c_aslp_topo_3, c_aslp_topo_4
  use atm_params, only : cp
  use atm_params, only : i_mmc_had, c_mmc_had, c_mmc_te, n_mmc_te
  use atm_params, only : i_mmc_fer, c_mmc_fer, c_mmc_fer_e, c_mmc_pol, c_mmc_1, c_mmc_2
  use atm_params, only : c_mmc_dt0, c_mmc_dt1, i_fzsa, c_mmc_z
  use atm_params, only : nsmooth_aslp, nsmooth_aslp_topo
  use atm_grid, only : im, jm, jmc, aim, jeq, jts, jtn, jps, jpn, dy, pl, k500
  use atm_grid, only : fcorua, sint, cost, fiu, fit
  use atm_grid, only : zc, k_dse_lo, k_dse_up
  use smooth_atm_mod, only : smooth2_m, smooth2eq
  !$ use omp_lib

  implicit none

  private
  public :: zslp, azslp 

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
          aslp_temp(i,j) = -c_slp_1*g*p0*10000._wp/(2._wp*Rd*T0**2)*atsl(i,j) 
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
  subroutine zslp(sin_cos_acbar, tsl, aslp, zsa, fdydse, fsydseg, t3, &
          slp, had_fi, had_width)

    implicit none

    real(wp), intent(in ) :: sin_cos_acbar(:,:)
    real(wp), intent(in ) :: tsl(:,:)
    real(wp), intent(in ) :: aslp(:,:)
    real(wp), intent(in ) :: zsa(:,:)
    real(wp), intent(in ) :: fdydse(:,:)
    real(wp), intent(in ) :: fsydseg(:)
    real(wp), intent(in ) :: t3(:,:,:)

    real(wp), intent(out) :: slp(:,:)
    real(wp), intent(out) :: had_fi
    real(wp), intent(out) :: had_width

    real(wp), parameter :: p6=pi/6._wp
    real(wp), parameter :: dse_min=1.e3_wp   !! J/kg, floor on the branch DSE contrast
    real(wp), parameter :: fi_dt_eq =10._wp*pi/180._wp  !! equatorial band, |fi| below this
    real(wp), parameter :: fi_dt_s1 =25._wp*pi/180._wp  !! subtropical band, equatorward limit
    real(wp), parameter :: fi_dt_s2 =35._wp*pi/180._wp  !! subtropical band, poleward limit

    integer :: i, j
    integer :: jc, jtn1, jtn2, jts1, jts2
    real(wp) :: wtn1, wtn2, wts1, wts2
    real(wp) :: tnh, tsh, scosh
    real(wp) :: dtpn, dtps, dtfn, dtfs, ficz, ttrpmx, tbhn, tbhs, dthn, dths
    real(wp) :: ttrpmn, wte, wtesum, fi_te
    real(wp) :: hadwidsc
    real(wp) :: dttrp, hadwid, teqsum, tsubsum, weqsum, wsubsum
    real(wp) :: ff, coc, psum, csum
    real(wp) :: dzdse

    real(wp) :: fisu(jmc), fist(jm), vsz(jmc), tslz(jm), psz(jmc), acbarz(jm)
    real(wp) :: dsez(jm), fedz(jmc), aferz(jm), cocfn(jm), cocfs(jm)
    real(wp) :: zsaz(jm), fzsat(jm), fzsa(jm)


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

    ! maximum tropical temperature, drives the Hadley branches
    ttrpmx = 0._wp
    do j=jtn,jts+1
      ttrpmx = max(ttrpmx,tslz(j))
    enddo

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

    ! ITCZ position, depends on temperature difference between the two hemispheres

    ficz = c_mmc_2*(tnh-tsh)

    ! Width of the Hadley cells.  The cells span ff in [-pi,pi], i.e. a total width of pi/(3*hadwidsc) radians
    ! The predictor is the tropical temperature CONTRAST. The total width is linear in
    !   dttrp = <tsl>(|fi|<10 deg) - <tsl>(25-35 deg),
    ! a stronger equator-to-subtropics contrast gives a narrower cell

    teqsum  = 0._wp; weqsum  = 0._wp
    tsubsum = 0._wp; wsubsum = 0._wp
    do j=1,jm
      if (abs(fit(j)).lt.fi_dt_eq) then
        teqsum  = teqsum  + tslz(j)*cost(j)
        weqsum  = weqsum  + cost(j)
      else if (abs(fit(j)).gt.fi_dt_s1 .and. abs(fit(j)).lt.fi_dt_s2) then
        tsubsum = tsubsum + tslz(j)*cost(j)
        wsubsum = wsubsum + cost(j)
      endif
    enddo
    dttrp = teqsum/max(weqsum,1.e-20_wp) - tsubsum/max(wsubsum,1.e-20_wp)

    ! total width in degrees, then hadwidsc from width = 180/(3*hadwidsc) deg
    hadwid   = c_mmc_dt0 - c_mmc_dt1*dttrp
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
    enddo

    ! Hadley cell width and ITCZ position 
    had_fi = 0.5_wp*((wtn1*fit(jtn1)+wtn2*fit(jtn2)) + (wts1*fit(jts1)+wts2*fit(jts2)))
    had_width = (wtn1*fit(jtn1)+wtn2*fit(jtn2)) - (wts1*fit(jts1)+wts2*fit(jts2)) 

    ! temperature gradients in the polar cells
    dtpn = (0.5_wp*tslz(jpn)+0.5_wp*tslz(jpn+1)) - tslz(1)
    dtps = (0.5_wp*tslz(jps)+0.5_wp*tslz(jps+1)) - tslz(jm)

    ! Amplitude of the Hadley branches, selected by i_mmc_had.
    !
    ! i_mmc_had = 1  The original closure.  Each branch is driven by the sea level temperature
    !                contrast between the warmest tropical latitude and the fixed 30 degree
    !                border of its own hemisphere,
    !                  coc_H = c_mmc_had * (ttrpmx - tslz(30 deg)).
    !                The two branches share ttrpmx and differ ONLY through tslz(30N) versus
    !                tslz(30S), a subtropical and largely continental contrast.  That has the
    !                right sign but the wrong seasonality: it splits the two cells about
    !                equally in DJF and in JJA, whereas the observed split is far more extreme
    !                in boreal summer.  Against ERA5 and the CMIP median the scheme gives a
    !                winter-to-summer cell ratio of 3.8 in DJF against 5.7 observed and 6.2 in
    !                JJA against 20.1, the southern JJA cell coming out at 16.7 against 19.8 at
    !                the pre-industrial and 15.2 against 21.4 at the LGM.
    !
    ! i_mmc_had = 2  The same contrast, scaled by the displacement of the THERMAL EQUATOR,
    !                  coc_H = c_mmc_had * (ttrpmx - tslz(30 deg)) * (1 -+ c_mmc_te*fi_te*|fi_te|),
    !                the minus taken in the northern hemisphere.  fi_te in degrees is the
    !                latitude of the tropical temperature maximum, so a thermal equator in the
    !                north weakens the northern cell and strengthens the southern
    !                cross-equatorial one, which is what the boreal summer does: the ascending
    !                branch is driven to 19 N by the northern land masses and the winter cell
    !                that reaches across the equator to it is correspondingly deep, while the
    !                austral summer only pulls it to 14 S over an almost entirely oceanic
    !                southern tropics.  This is the asymmetry the incumbent has no way of
    !                expressing, and it is what the c_mmc_2 ITCZ excursion alone cannot supply,
    !                that entering the cell GEOMETRY and being limited to a few degrees before
    !                the shape function distorts.
    !
    !                The response is QUADRATIC in the displacement, not linear, because that is
    !                what the cells ask for.  The austral summer thermal equator reaches 6-7 deg
    !                and the winter cell then needs no correction at all; the boreal summer one
    !                reaches 11-13 deg and the winter cell needs 27 per cent at the
    !                pre-industrial and 40 at the LGM.  A linear factor tuned on JJA therefore
    !                over-corrects DJF and drives the southern summer cell from 4.5 to 2.9
    !                against a reference of 4.1, where the quadratic one leaves it at 3.7.
    !
    ! c_mmc_had needs retuning between the two: the displacement factor raises the winter
    ! branches, which carry nearly all of the circulation, so the coefficient comes down by
    ! about a tenth at c_mmc_te = 0.002.
    tbhn = 0.5*tslz(jtn)+0.5*tslz(jtn+1)
    tbhs = 0.5*tslz(jts)+0.5*tslz(jts+1)
    dthn = max(0._wp,ttrpmx-tbhn)
    dths = max(0._wp,ttrpmx-tbhs)

    if (i_mmc_had.eq.2) then

      ! Latitude of the thermal equator: the cos-weighted centroid of the tropical sea level
      ! temperature measured above the coldest point of the same band, raised to n_mmc_te.  The
      ! exponent sharpens the centroid towards the true maximum; it is not the maximum itself,
      ! because the argmax on this grid moves in 5 degree steps and a branch amplitude that
      ! jumps by that much from one time step to the next is not usable.  n_mmc_te = 4 tracks
      ! the maximum to about a degree in the annual mean and stays smooth.
      ttrpmn = ttrpmx
      do j=jtn,jts+1
        ttrpmn = min(ttrpmn,tslz(j))
      enddo
      wtesum = 0._wp
      fi_te  = 0._wp
      do j=jtn,jts+1
        wte    = max(0._wp,tslz(j)-ttrpmn)**n_mmc_te * cost(j)
        wtesum = wtesum + wte
        fi_te  = fi_te + wte*fit(j)
      enddo
      if (wtesum.gt.0._wp) then
        fi_te = fi_te/wtesum * 180._wp/pi     ! deg
      else
        fi_te = 0._wp
      endif

      dthn = dthn * max(0._wp,1._wp-c_mmc_te*fi_te*abs(fi_te))
      dths = dths * max(0._wp,1._wp+c_mmc_te*fi_te*abs(fi_te))

    else if (i_mmc_had.ne.1) then

      stop 'i_mmc_had'

    endif

    ! Amplitude of the Ferrel branches, selected by i_mmc_fer.
    !
    ! i_mmc_fer = 1  The original closure.  The amplitude is proportional to the meridional
    !                sea level temperature contrast ACROSS the cell,
    !                  psi_F = c_mmc_fer * dtf,   dtf = tslz(Hadley edge) - tslz(polar edge),
    !                one scalar per hemisphere, so the amplitude is uniform along the branch.
    !
    ! i_mmc_fer = 2  Transformed Eulerian mean form, resolved in latitude.  The Ferrel cell is
    !                the Eulerian residue of the baroclinic eddy fluxes,
    !                  psi_F ~ -(2 pi a cos(fi)/g) [v'th'] / (dth/dp),
    !                i.e. the mass circulation that returns, against the mean dry static energy
    !                gradient, the DSE the synoptic eddies transport poleward.  Written as a
    !                bulk balance over the depth of the cell that is
    !                  psi_F(j) = c_mmc_fer_e * |F_eddy(j)| / Ds(j),
    !                with F_eddy the zonally integrated synoptic eddy DSE flux (W) and Ds the
    !                DSE contrast between the two branches (J/kg).
    !
    ! i_mmc_fer = 3  As 2, but F_eddy is the TOTAL eddy DSE flux, transient plus stationary,
    !                instead of the transient (diffusive) part alone.  Transformed Eulerian
    !                mean theory has the residual circulation respond to the whole eddy flux,
    !                and closure 2 sees only half of it: the standing waves carry a large
    !                poleward DSE flux in the northern hemisphere and almost none in the
    !                southern, so leaving them out builds in the wrong hemispheric asymmetry.
    !                With the diffusive flux alone the model gives a winter northern-to-
    !                southern Ferrel ratio of 0.66 against 0.86-1.10 in ERA5 and the CMIP6
    !                piControl ensemble.  fsydseg is the stationary part with the mean
    !                meridional circulation already removed in adifa, which it must be: that
    !                part IS the cell being parameterised.
    !
    ! NOTE c_mmc_fer and c_mmc_fer_e multiply completely different quantities (K versus
    ! 10^10 kg/s) and are NOT interchangeable; each carries its own tuned value.

    if (i_mmc_fer.eq.1) then

      ! meridional sea level temperature contrast across each Ferrel cell.  The cell edges are
      ! the fixed grid boundaries jtn/jpn, as for the Hadley and polar contrasts above.
      dtfn = (0.5_wp*tslz(jtn)+0.5_wp*tslz(jtn+1)) - (0.5_wp*tslz(jpn)+0.5_wp*tslz(jpn+1))
      dtfs = (0.5_wp*tslz(jts)+0.5_wp*tslz(jts+1)) - (0.5_wp*tslz(jps)+0.5_wp*tslz(jps+1))
      do j=2,jm
        cocfn(j) = c_mmc_fer*dtfn
        cocfs(j) = c_mmc_fer*dtfs
      enddo

    else if (i_mmc_fer.eq.2 .or. i_mmc_fer.eq.3) then

      ! zonally integrated transient eddy DSE flux, W (fdydse carries kg/s*K)
      fedz(:) = 0._wp
      do j=1,jmc
        do i=1,im
          fedz(j) = fedz(j) + fdydse(i,j)*cp
        enddo
      enddo

      ! add the standing wave flux, already zonally integrated and already free of the mean
      ! meridional circulation.  The sum is taken before the absolute value so that the two
      ! can oppose each other where they physically do.
      if (i_mmc_fer.eq.3) then
        do j=1,jmc
          fedz(j) = fedz(j) + fsydseg(j)*cp
        enddo
      endif

      do j=2,jm
        aferz(j) = c_mmc_fer_e * abs(fedz(j))/(0.5_wp*(dsez(j-1)+dsez(j))) * 1.e-10_wp
      enddo

      do j=2,jm
        cocfn(j) = aferz(j)
        cocfs(j) = aferz(j)
      enddo

    else

      stop 'i_mmc_fer'

    endif

    cocfn(1) = 0._wp
    cocfs(1) = 0._wp

    ! Topographic mask on the amplitude of all cell branches.  Elevated topography blocks the
    ! near-surface meridional flow that closes the cells, so the amplitude is reduced where the
    ! surface stands high above sea level.  fzsa is on the u grid, where the amplitudes are used.

    if (i_fzsa.eq.0) then

      ! no mask
      fzsa(:) = 1._wp

    else if (i_fzsa.eq.1) then

      ! zonal mean elevation first, so that a partly elevated latitude circle is masked as a whole
      zsaz(:) = 0._wp
      do j=1,jm
        do i=1,im
          zsaz(j) = zsaz(j) + zsa(i,j)*aim
        enddo
      enddo
      fzsa(1) = 1._wp - min(1._wp,max(0._wp,zsaz(1)/c_mmc_z))
      do j=2,jm
        fzsa(j) = 1._wp - min(1._wp,max(0._wp,0.5_wp*(zsaz(j-1)+zsaz(j))/c_mmc_z))
      enddo

    else if (i_fzsa.eq.2) then

      ! per grid cell first, then zonal mean, so that a latitude circle which is only partly
      ! covered by high topography keeps the flow over the low part
      fzsat(:) = 0._wp
      do j=1,jm
        do i=1,im
          fzsat(j) = fzsat(j) + max(0._wp,1._wp-zsa(i,j)/c_mmc_z)*aim
        enddo
      enddo
      fzsa(1) = fzsat(1)
      do j=2,jm
        fzsa(j) = 0.5_wp*(fzsat(j-1)+fzsat(j))
      enddo

    else

      stop 'i_fzsa'

    endif

    ! Surface meridional ageostrofic wind

    vsz(:) = 0._wp

    do j=2,jm

      ff = 6._wp*hadwidsc*(fiu(j)-ficz/(c_mmc_1*(fiu(j)-ficz)**2+1._wp))

      coc = 0._wp

      if (ff.ge.0._wp .and. ff.lt.pi)              coc = c_mmc_had*dthn
      if (ff.ge.pi .and. ff.lt.2._wp*pi)           coc = cocfn(j)
      if (ff.ge.2._wp*pi .and. ff.lt.3._wp*pi)     coc = c_mmc_pol*dtpn

      if (-ff.gt.0._wp .and. (-ff).le.pi)          coc = c_mmc_had*dths
      if (-ff.gt.pi .and. (-ff).le.2._wp*pi)       coc = cocfs(j)
      if (-ff.gt.2._wp*pi .and. (-ff).le.3._wp*pi) coc = c_mmc_pol*dtps

      vsz(j) = -coc*fzsa(j)*sin(ff)

    enddo

    ! Zonally averaged SLP
    ! Integration from NP to SP

    psz(:) = 0._wp
    do j=2,jm
      psz(j) = psz(j-1) + vsz(j)*fcorua(j)*ra*dy/acbarz(j)
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

end module slp_mod
