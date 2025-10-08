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
  use atm_params, only : l_aslp_topo, c_aslp_topo_1, c_aslp_topo_2, c_aslp_topo_3, c_aslp_topo_4
  use atm_params, only : i_mmc, c_mmc_had, c_mmc_fer, c_mmc_pol, c_mmc_z, c_mmc_1, c_mmc_2, c_mmc_3, c_mmc_4
  use atm_params, only : nsmooth_aslp, nsmooth_aslp_eq, nj_eq, nsmooth_aslp_topo
  use atm_grid, only : im, jm, jmc, aim, jeq, jts, jtn, jps, jpn, dy, pl, k500, i_ice
  use atm_grid, only : fcorua, sint, cost, fiu, fit
  use smooth_atm_mod, only : smooth2_m, smooth2eq, zona
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
    real(wp), dimension(im) :: eps, Kn2
    real(wp) :: zsa_smooth(im,jm)
    real(dp) :: zsa_smooth_dp(im,jm)
    real(dp), dimension(im) :: psi
    complex(dp), dimension(im/2+1) :: zsa_fft
    complex(dp), dimension(im/2+1) :: psi_fft


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
    call smooth2eq(aslp,nj_eq,nsmooth_aslp_eq)
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
  subroutine zslp(zsa, sin_cos_acbar, tsl, aslp, &
          slp, had_fi, had_width)

    implicit none

    real(wp), intent(in ) :: zsa(:,:)
    real(wp), intent(in ) :: sin_cos_acbar(:,:)
    real(wp), intent(in ) :: tsl(:,:)
    real(wp), intent(in ) :: aslp(:,:)

    real(wp), intent(out) :: slp(:,:)
    real(wp), intent(out) :: had_fi
    real(wp), intent(out) :: had_width

    real(wp), parameter :: p6=pi/6._wp

    integer :: i, j
    integer :: jc, jtn1, jtn2, jpn1, jpn2, jts1, jts2, jps1, jps2
    real(wp) :: wtn1, wtn2, wpn1, wpn2, wts1, wts2, wps1, wps2
    real(wp) :: tnh, tsh, scosh
    real(wp) :: ttrp, dtpn, dtps, dtfn, dtfs, ficz, ttrpmx, tbhn, tbhs, dthn, dths
    real(wp) :: hadwidsc
    real(wp) :: ff, coc, psum, csum
   
    real(wp) :: fisu(jmc), fist(jm), vsz(jmc), tslz(jm), psz(jmc), acbarz(jm), zsaz(jm), fzsa(jm)


    ! zonal mean sea level temperature and topography 
    tslz(:) = 0._wp
    zsaz(:) = 0._wp
    do j=1,jm 
      do i=1,im
        tslz(j) = tslz(j) + tsl(i,j)*aim 
        zsaz(j) = zsaz(j) + zsa(i,j)*aim 
      enddo        
    enddo

    ! factor to account for zonal mean topography in the PBL
    do j=1,jm 
      fzsa(j) = 1._wp-min(1._wp,max(0._wp,0.5_wp*(zsaz(max(1,j-1))+zsaz(j))/c_mmc_z))
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

    ! tropical temperature  
    ttrp = zona(tslz,cost,jtn+1,jts)
    ttrp = max(ttrp,c_mmc_4+50._wp)
    ttrpmx = 0._wp
    do j=jtn,jts+1
      ttrpmx = max(ttrpmx,tslz(j))
    enddo 

    ! ITCZ position, depends on temperature difference between the two hemispheres 

    ficz = c_mmc_2*(tnh-tsh)

    ! width of Hadley cell, increases with tropical temperature 

    ! scaling of Hadley cell width
    hadwidsc = c_mmc_3/(ttrp-c_mmc_4)
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

    if (i_mmc.eq.1) then

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

    else if (i_mmc.eq.2) then

      ! temperature gradients in the polar cells
      dtpn = (wpn1*tslz(jpn1)+wpn2*tslz(jpn2)) - tslz(1)
      dtps = (wps1*tslz(jps1)+wps2*tslz(jps2)) - tslz(jm)

      ! temperature gradients in the ferrel cells
      dtfn = (wtn1*tslz(jtn1)+wtn2*tslz(jtn2)) - (wpn1*tslz(jpn1)+wpn2*tslz(jpn2))
      dtfs = (wts1*tslz(jts1)+wts2*tslz(jts2)) - (wps1*tslz(jps1)+wps2*tslz(jps2))

      ! temperature gradients in Hadley cells    
      tbhn = wtn1*tslz(jtn1)+wtn2*tslz(jtn2)
      tbhs = wts1*tslz(jts1)+wts2*tslz(jts2)
      dthn = max(0._wp,ttrpmx-tbhn)
      dths = max(0._wp,ttrpmx-tbhs)

    endif

    ! Surface meridional ageostrofic wind 

    vsz(:) = 0._wp

    do j=2,jm

      ff = 6._wp*hadwidsc*(fiu(j)-ficz/(c_mmc_1*(fiu(j)-ficz)**2+1._wp))

      coc = 0._wp

      if (ff.ge.0._wp .and. ff.lt.pi)              coc = c_mmc_had*dthn * fzsa(j) 
      if (ff.ge.pi .and. ff.lt.2._wp*pi)           coc = c_mmc_fer*dtfn * fzsa(j) 
      if (ff.ge.2._wp*pi .and. ff.lt.3._wp*pi)     coc = c_mmc_pol*dtpn * fzsa(j) 

      if (-ff.gt.0._wp .and. (-ff).le.pi)          coc = c_mmc_had*dths * fzsa(j) 
      if (-ff.gt.pi .and. (-ff).le.2._wp*pi)       coc = c_mmc_fer*dtfs * fzsa(j) 
      if (-ff.gt.2._wp*pi .and. (-ff).le.3._wp*pi) coc = c_mmc_pol*dtps * fzsa(j) 

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

end module slp_mod
