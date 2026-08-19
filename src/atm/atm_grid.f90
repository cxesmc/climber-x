!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : a t m _ g r i d
!
!  Purpose : definition of atmospheric grid
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
module atm_grid

  use atm_params, only : wp
  use nml, only : nml_read
  use constants, only : pi, r_earth, omega, g, Rd, T0
  use climber_grid, only: ni, nj, dlat
  use control, only : out_dir
  use atm_params, only : atm_mass, hatm, amas, ra, hcld_base, fcormin, i_fcorg
  use atm_params, only : l_p0_var, p0, ps0, vprof_tope, vprof_topp
  use atm_params, only : i_uter_pol, c_uter_pol, c_uter_pol_min, lat_uter_pol_1, lat_uter_pol_2
  use smooth_atm_mod, only : smooth2

  implicit none

  integer, parameter :: im = ni
  integer, parameter :: jm = nj
  integer, parameter :: imc = im+1
  integer, parameter :: jmc = jm+1
  integer, parameter :: jeq=jm/2
  integer, parameter :: jpn=jm/6
  integer, parameter :: jtn=jm/3
  integer, parameter :: jts=jm*2/3
  integer, parameter :: jps=jm*5/6
  integer :: km
  integer :: kmc
  integer, parameter :: nm = 5  !! number of macro surface types
  !! index of surface types
  integer, parameter :: i_ocn = 1
  integer, parameter :: i_sic = 2
  integer, parameter :: i_lnd = 3
  integer, parameter :: i_ice = 4  
  integer, parameter :: i_lake = 5
  real(wp) :: zsa_scale
  real(wp) :: zsa_scale_dyn
  integer :: nsmooth_zsa
  integer :: llwr
  integer :: nlwr1
  integer :: nlwr2
  integer :: nlwr3
  integer :: nlwr4
  integer :: llwr1
  integer :: llwr2
  integer :: llwr3
  integer :: llwr4
      
  real(wp) :: fit(jm)
  real(wp) :: fiu(jmc)
  real(wp) :: thetat(jm)
  real(wp) :: cost(jm)
  real(wp) :: sint(jm)
  real(wp) :: costhetat(jm)
  real(wp) :: sinthetat(jm)
  real(wp) :: signf(jm)
  real(wp) :: cosu(jmc)
  real(wp) :: sinu(jmc)      
  real(wp) :: dxt(jm)
  real(wp) :: dxu(jmc)
  real(wp) :: sqr(im,jm)
  real(wp) :: esqr      !! Earth surface area
  real(wp) :: dy
  real(wp) :: aim
 
  real(wp) :: fcort(jm)
  real(wp) :: fcorg(jm)      !! s, the reciprocal-Coriolis factor of the geostrophic PBL wind, see i_fcorg
  real(wp) :: fcorgu(jmc)    !! s, the same factor on the v-face latitudes, for the streamfunction form (i_ugb_psi=1)
  real(wp) :: fcorta_sqrt(jm)
  real(wp) :: fcorta(jm)
  real(wp) :: fcorua(jmc)
  real(wp) :: cdamp_pol(jm)  !! 1, the polar amplitude reduction of the thermal wind, see i_uter_pol
       
  real(wp) :: plx(imc,jm)
  real(wp) :: plx_trop(imc,jm)   !! tropospheric (k<=km-2) zonal column mass, for implicit zonal diffusion
  real(wp) :: ply(im,jmc)
  real(wp) :: ply_trop(im,jmc)   !! tropospheric (k<=km-2) meridional column mass, for implicit meridional diffusion
  real(wp) :: ptopt(jm)     !! top of the mean meridional cell in t-points
  real(wp) :: ptopu(jmc)    !! top of the mean meridional cell in u-points
  integer :: k1(im,jm)  
  integer :: kweff(im,jm)  !! k-index for effective vertical velocity leve (w-grid)
  integer :: k1000
  integer :: k900
  integer :: k850
  integer :: k700
  integer :: k500
  integer :: k300
  ! layer centres bracketing the two branches of the Ferrel cell, used by slp_mod::zslp to form
  ! the dry static energy contrast of the eddy closure
  integer :: k_dse_lo
  integer :: k_dse_up

  real(wp), allocatable :: pl(:)
  real(wp), allocatable :: dpl(:)
  real(wp), allocatable :: zl(:)
  real(wp), allocatable :: zc(:)
  real(wp), allocatable :: dplx(:,:,:)
  real(wp), allocatable :: dply(:,:,:)
  real(wp), allocatable :: exp_zc(:)

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  a t m _ g r i d _ i n i t 
  !   Purpose    :  initialize atmospheric grid
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine atm_grid_init

    implicit none

    integer :: i, j, k
    real(wp) :: fcortp, fcorup
    real(wp) :: xpol
    character(len=256) :: fnm


    fnm = trim(out_dir)//"/atm_par.nml"
    call nml_read(fnm,"atm_par","llwr",llwr)
    call nml_read(fnm,"atm_par","km",km)
    call nml_read(fnm,"atm_par","zsa_scale",zsa_scale)
    call nml_read(fnm,"atm_par","zsa_scale_dyn",zsa_scale_dyn)
    call nml_read(fnm,"atm_par","nsmooth_zsa",nsmooth_zsa)
    kmc = km+1

    ! allocate
    allocate(pl(kmc))
    allocate(dpl(km))
    allocate(zl(kmc))
    allocate(zc(km))
    allocate(dplx(im,jm,km))
    allocate(dply(im,jm,km))
    allocate(exp_zc(km))

    ! read pressure levels from namelist
    call nml_read(fnm,"atm_par","pl",pl)

    ! layers for longwave radiation
    if (mod(llwr,5).ne.0) then
      print *,'abort llwr is not a multiple of 5'
      stop
    endif      
    nlwr1 = llwr/5*2 
    nlwr2 = llwr/5 
    nlwr3 = llwr/5 
    nlwr4 = llwr/5 
    llwr1 = nlwr1
    llwr2 = llwr1+nlwr2
    llwr3 = llwr2+nlwr3
    llwr4 = llwr3+nlwr4 
    if (llwr4.ne.llwr) then
      print *,'abort llwr4 ne llwr'
      stop
    endif       

    ! grid size and area   

    aim = 1._wp/im   

    dy = pi*r_earth/jm

    do j=1,jm
      fit(j) = pi*(90._wp+dlat/2._wp-dlat*j)/180._wp
      thetat(j) = pi/2.-fit(j)
      cost(j) = cos(fit(j))
      sint(j) = sin(fit(j))
      costhetat(j) = cos(thetat(j))
      sinthetat(j) = sin(thetat(j))
      dxt(j) = dy*cost(j)
    enddo

    do j=1,jmc
      fiu(j) = pi*(90._wp+dlat-dlat*j)/180._wp
      cosu(j) = cos(fiu(j))
      sinu(j) = sin(fiu(j))       
      dxu(j) = dy*cosu(j)
    enddo

    do j=1,jm
      do i=1,im
        sqr(i,j) = dxt(j)*dy
      enddo

      if (j.le.jm/2) then
        signf(j) = 1._wp
      else
        signf(j) = -1._wp
      endif

    enddo 

    esqr = 0._wp
    do i=1,im
      do j=1,jm
        esqr = esqr+sqr(i,j)
      enddo
    enddo

    ! Coriolis parameter

    ! fcort is floored in MAGNITUDE at fcormin but keeps its sign, so any quantity built as
    ! 1/fcort reverses sign discontinuously across the equator at the full floored magnitude.
    ! With fcormin = 1e-5 the floor binds for |lat| < 3.93 deg, i.e. on the two rows either
    ! side of the equator, and that is where the geostrophic PBL wind ugb jumps.
    !
    ! fcorg is the factor u2d.f90 actually divides the SLP gradient by, and i_fcorg chooses how
    ! the equatorial singularity in it is handled:
    !
    !   i_fcorg = 0   fcorg = 1/fcort, the original.  Floored in magnitude, so |fcorg| is at its
    !                 LARGEST on the two rows either side of the equator and changes sign between
    !                 them: ugb jumps by 2*|dpdy|/(fcormin*ra), which is about 9.5 m/s in the
    !                 zonal mean at the solstices.
    !   i_fcorg = 1   fcorg = f/(f**2 + fcormin**2), the Rayleigh-damped geostrophic balance, in
    !                 which fcormin doubles as the linear drag rate r.  This is the same
    !                 regularisation but with the opposite behaviour at the equator: |fcorg| peaks
    !                 at 1/(2*fcormin) where |f| = fcormin and goes to ZERO at f = 0, so ugb
    !                 crosses zero continuously instead of reversing at full amplitude.  It is
    !                 also the factor consistent with uab/vab below, which already carry the exact
    !                 down-gradient companion -dpdx*<sin a cos a>/(|f|*ra) of the same balance.
    !
    ! Using fcormin as the drag rate is dimensionally the right scale: cd*|V|/h_pbl with
    ! cd = 1.3e-3, |V| = 7 m/s and h_pbl = 1 km gives 9e-6, essentially the 1e-5 the floor
    ! already carries.  |f| = fcormin at |lat| = 3.93 deg.
    !
    ! Rebuilding us from the archived ugb, vgb and acbar of output/cacbar/pi (the reconstruction
    ! reproduces the archived us to 4%) gives, for the PI zonal mean:
    !
    !     drag rate     jump at eq DJF     u rms |lat|<25   DJF / JJA / ANN
    !     floor 1e-5          9.58 m/s                 2.23 / 2.12 / 0.59
    !     1.0e-5              4.34                     1.36 / 1.35 / 0.88
    !     1.5e-5              2.30                     1.29 / 1.31 / 1.08
    !     2.0e-5              1.38                     1.44 / 1.43 / 1.26
    !
    ! so the solstice seasons want 1.0-1.5e-5.  The annual mean gets WORSE because it was living
    ! off the cancellation of the two seasonal jumps, which have opposite sign; once that
    ! cancellation is gone the model's equatorial easterlies show up as too weak (-0.6 against
    ! -1.6 m/s in CMIP5 at 2.5N).  That bias is pre-existing and was hidden, not created here,
    ! and the missing easterly belongs in the down-gradient term uab, not in geostrophy.
    !
    ! The other branch of the wind that divides by fcort, the thermal wind of u3d.f90:496, is
    ! already regular at the equator: uter is multiplied by c_damp_eq = min(1,c_uter_eq*sin**2),
    ! so its 1/f is cancelled to leading order and it goes to zero at f = 0 - qualitatively what
    ! i_fcorg=1 now does for ugb.  It was only the barotropic branch that was left on the floor.
    !
    ! NOTE fcormin now sets both the magnitude floor of fcort/fcorta/fcorua and, with i_fcorg=1,
    ! the drag rate of fcorg, so raising it widens the floored band as well: at 4.2e-5 the floor
    ! would bind out to |lat| = 16.7 deg and flatten the whole trade-wind belt.
    do j=1,jm
      fcortp = 2._wp*omega*sint(j)
      fcort(j) = signf(j)*max(ABS(fcortp),fcormin)
      fcorta(j) = max(ABS(fcortp),fcormin)
      fcorta_sqrt(j) = sqrt(abs(fcorta(j)))
      if (i_fcorg.eq.0) then
        fcorg(j) = 1._wp/fcort(j)
      else if (i_fcorg.eq.1) then
        fcorg(j) = fcortp/(fcortp**2 + fcormin**2)
      else
        stop 'i_fcorg'
      endif
    enddo

    do j=1,jm
      fcorup = 2._wp*omega*sinu(j)
      fcorua(j) = max(ABS(fcorup),fcormin)
    enddo

    ! The same reciprocal-Coriolis factor as fcorg, but evaluated on the v-face latitudes.
    ! The streamfunction form of the barotropic geostrophic wind (i_ugb_psi=1) carries the
    ! streamfunction on the cell corners, which sit at fiu, so it needs the factor there.
    ! Both poles are included: fiu(1) and fiu(jmc) are +-90 deg, where dxu vanishes and the
    ! streamfunction is held zonally constant so that no mass crosses the pole.
    do j=1,jmc
      fcorup = 2._wp*omega*sinu(j)
      if (i_fcorg.eq.0) then
        fcorgu(j) = 1._wp/(sign(1._wp,fcorup)*max(ABS(fcorup),fcormin))
      else if (i_fcorg.eq.1) then
        fcorgu(j) = fcorup/(fcorup**2 + fcormin**2)
      else
        stop 'i_fcorg'
      endif
    enddo

    ! Polar amplitude reduction of the thermal wind.
    !
    ! With i_uter_damp=2 this factor multiplies the azonal TEMPERATURE, so the wind picks up
    ! the extra term -K*T_az*dc/dy on top of the damped thermal wind -K*c*dT_az/dy.  That term
    ! is proportional to the anomaly itself rather than to its gradient, so over a broad warm
    ! anomaly it is a monopole where the physical thermal wind is a dipole, and it does not
    ! weaken the wave, it adds a different one.  Everything therefore depends on WHERE dc/dy is
    ! put, not just on how much total damping there is.
    !
    !   i_uter_pol = 0   c = min(1, c_uter_pol*cos(lat)**2), the original.  The min() makes
    !                    dc/dy vanish equatorward of 54.7 deg (for c_uter_pol=3) and then switch
    !                    on discontinuously, and c keeps falling all the way to the pole, so
    !                    dc/dy is non-zero over the whole 55-90 deg band - which in JJA is
    !                    exactly where the Siberian azonal temperature maximum sits.  Measured
    !                    from the t3 of output/iuterdamp/itrdmp.2, the resulting spurious zonal
    !                    wind is 1.96 m/s rms over 55-140E / 50-75N, 0.98 over the NH 30-87N.
    !
    !   i_uter_pol = 1   c = 1 equatorward of lat_uter_pol_1, c = c_uter_pol_min poleward of
    !                    lat_uter_pol_2, joined by a quintic smoothstep.  dc/dy is then exactly
    !                    zero inside the polar cap and in mid-latitudes, and confined to the
    !                    transition band, where it can be spread as thinly as wanted by widening
    !                    the band.  c_uter_pol is not used.  With the namelist defaults the same
    !                    diagnostic gives 0.88 m/s rms over Siberia and 0.52 over the NH, i.e.
    !                    -55% and -47%, at the same area-weighted total damping poleward of
    !                    55 deg (0.513 against 0.505).
    !
    ! A constant c also repairs the defect of i_uter_damp=1, for the same reason.  Scaling the
    ! wind by c is identical to computing the thermal wind with f_eff = f/c, so mode 1 is not
    ! unbalanced at all - it is the geostrophic wind of a planet whose Coriolis parameter is
    ! f/c, and its divergence is exactly that planet's beta term.  What is wrong is the size of
    ! that beta: beta_eff/beta = 1 + 2*tan(lat)^2, i.e. 7 at 60 deg, 29 at 75, 116 at 82.5.
    ! Where c is constant f_eff = f/c_uter_pol_min, so beta_eff/f_eff = beta/f EXACTLY; the
    ! discrete maximum over 40-85N falls from 140 to 4.  Inside the cap both modes are therefore
    ! clean and the choice of i_uter_damp stops mattering there.
    !
    ! BUT a constant c is not admissible in the polar cap on this grid, which is why i_uter_pol
    ! defaults to 0.  The damping exists to keep the zonal Courant number u*tstep/dxt finite
    ! (dxt = 24 km at 87.5 deg against 556 at the equator) and to kill the 1/cos(lat) singularity
    ! of v_ter, which is built from dT/dx = dT/(2*dxt).  c ~ cos(lat)**2 does both: Cx ~ 3cos(lat),
    ! bounded and decreasing poleward, and v_ter ~ cos(lat).  A constant 0.45 gives Cx = 0.98 at
    ! 82.5 deg and max|v_ter| there of 12 m/s against 1.4 now; c ~ cos gives Cx = 1.3-1.45 at
    ! 67.5-77.5 deg, outright unstable.  So c must fall at least as fast as cos(lat)**2 near the
    ! pole, and any such c has |dln(c)/dlat| >= 2*tan(lat).  Polar CFL, a small dc/dy and a
    ! multiplicative latitude factor cannot all hold at once; the way out is a limiter or a polar
    ! zonal filter that acts only where the Courant number is actually violated, not a smooth
    ! amplitude factor spread over the whole cap.
    !
    ! The smoothstep 10x^3-15x^4+6x^5 has zero first AND second derivative at both ends, so the
    ! join carries no kink of the kind min() introduces.
    do j=1,jm
      if (i_uter_pol.eq.0) then
        cdamp_pol(j) = min(1._wp, c_uter_pol*cost(j)**2)
      else if (i_uter_pol.eq.1) then
        xpol = (abs(fit(j))*180._wp/pi - lat_uter_pol_1) / (lat_uter_pol_2-lat_uter_pol_1)
        xpol = min(1._wp, max(0._wp, xpol))
        cdamp_pol(j) = 1._wp - (1._wp-c_uter_pol_min) * xpol**3*(10._wp+xpol*(6._wp*xpol-15._wp))
      else
        stop 'i_uter_pol'
      endif
    enddo

    ! Top of the mean meridional cells, used by the ageostrophic wind profile in u3d.
    ! The observed 5 % level of the streamfunction sits at 121 hPa in the Hadley cell, 155 in
    ! the Ferrel cell and 177 in the polar cell, and moves by only 15-33 hPa between PI, LGM
    ! and aquaplanet, so this is a fixed function of latitude with no climate dependence.
    do j=1,jm
      ptopt(j) = vprof_topp-(vprof_topp-vprof_tope)*cost(j)**2
      ptopu(j) = vprof_topp-(vprof_topp-vprof_tope)*cosu(j)**2
    enddo
    ! u3d never reaches j=jmc, but this slot should not be left undefined
    ptopu(jmc) = vprof_topp

    ! pressure levels

    ra = p0/(Rd*T0)     ! kg/m3, air density at pressure p0 and temperature T0

    ps0 = atm_mass/esqr*g               !! Pa, average surface pressure
   
    do k=1,km
      dpl(k) = pl(k)-pl(k+1)
    enddo

    ! model z-levels

    do k=1,km
      zl(k) = -hatm*log(pl(k))
    enddo
    zl(kmc) = 30.e3_wp

    ! layer centers
    do k=1,km
      zc(k) = 0.5_wp*(zl(k+1)+zl(k))
      exp_zc(k) = exp(-zc(k)/hatm)
    enddo

    ! k-index of selected pressure levels
    k1000 = 1
    k900  = minloc(abs(pl-0.9_wp),1)
    k850  = minloc(abs(pl-0.85_wp),1)
    k700  = minloc(abs(pl-0.7_wp),1)
    k500  = minloc(abs(pl-0.5_wp),1)
    k300  = minloc(abs(pl-0.3_wp),1)

    ! layer centres nearest the lower and upper branch of the Ferrel cell (~1300 and ~7300 m,
    ! i.e. roughly 850 and 400 hPa); zc, not zl, because t3 is carried at layer centres
    k_dse_lo = minloc(abs(zc-1300._wp),1)
    k_dse_up = minloc(abs(zc-7300._wp),1)

    print *,'k 1000 hPa',k1000, ', z 1000 hPa',zl(k1000)
    print *,'k 900  hPa',k900,  ', z 900  hPa',zl(k900)
    print *,'k 850  hPa',k850,  ', z 850  hPa',zl(k850)
    print *,'k 700  hPa',k700,  ', z 700  hPa',zl(k700)
    print *,'k 500  hPa',k500,  ', z 500  hPa',zl(k500)
    print *,'k 300  hPa',k300,  ', z 500  hPa',zl(k300)
    print *,'k dse lo  ',k_dse_lo, ', zc      ',zc(k_dse_lo)
    print *,'k dse up  ',k_dse_up, ', zc      ',zc(k_dse_up)

    ! initialize, needed by vesta
    kweff(:,:) = 4

    return

  end subroutine atm_grid_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  a t m _ g r i d _ u p d a t e
  !   Purpose    :  update atmospheric grid
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine atm_grid_update(zs, frst, zsa, &
      zsa_smooth, slope, slope_x, slope_y, ra2, ra2a, pzsa0, pzsa, psa, ps)

    implicit none

    real(wp), intent(inout) :: zs(:,:,:)
    real(wp), intent(in   ) :: frst(:,:,:)
    real(wp), intent(inout) :: zsa(:,:)

    real(wp), intent(out  ) :: zsa_smooth(:,:)
    real(wp), intent(out  ) :: slope(:,:)
    real(wp), intent(out  ) :: slope_x(:,:)
    real(wp), intent(out  ) :: slope_y(:,:)
    real(wp), intent(out  ) :: ra2(:,:,:)
    real(wp), intent(out  ) :: ra2a(:,:)
    real(wp), intent(out  ) :: pzsa0(:,:)
    real(wp), intent(out  ) :: pzsa(:,:)
    real(wp), intent(out  ) :: psa(:,:)
    real(wp), intent(out  ) :: ps(:,:,:)

    integer :: i, j, k, n, imi, ipl
    real(wp) :: dp, px, py


    zs = zsa_scale * zs
    zsa = zsa_scale * zsa

    ! scale/smooth topography for dynamics
    zsa_smooth = zsa
    zsa_smooth = zsa_scale_dyn*zsa_smooth
    call smooth2(zsa_smooth,nsmooth_zsa)

    ! topography slopes
    do i=1,im
      imi=i-1
      if (imi.lt.1) imi=im  
      ipl=i+1
      if (ipl.gt.im) ipl=1
      do j=1,jm
        slope_x(i,j) = (min(3000._wp,zsa(i,j))-min(3000._wp,zsa(imi,j)))/dxt(j)    ! on u-grid
      enddo
      slope_y(i,1) = 0._wp
      do j=2,jm
        slope_y(i,j) = (min(3000._wp,zsa(i,j-1))-min(3000._wp,zsa(i,j)))/dy        ! on v-grid
      enddo
    enddo

    do i=1,im
      imi=i-1
      if (imi.lt.1) imi=im  
      ipl=i+1
      if (ipl.gt.im) ipl=1
      do j=2,jm-1
        slope(i,j) = sqrt((0.5_wp*(slope_x(i,j)+slope_x(ipl,j)))**2+(0.5_wp*(slope_y(i,j)+slope_y(i,j+1)))**2)
      enddo
      slope(i,1) = 0._wp
      slope(i,jm) = 0._wp
    enddo
    call smooth2(slope,1)

    ! pressure at surface of smooth topography
    pzsa0 = exp(-zsa/hatm)

    ! pressure at surface of smooth topography
    pzsa = exp(-zsa_smooth/hatm)

    if (l_p0_var) then

      ! compute mean sea level pressure from topography and mean surface pressure (conserved)
      dp = 0._wp
      do i=1,im
        do j=1,jm
          dp = dp + exp(-zsa(i,j)/hatm) * sqr(i,j)/esqr
        enddo
      enddo
      p0 = ps0/dp ! Pa, average sea level pressure

    endif

    amas = p0/g ! kg/m2, average mass of atmospheric column

    ra = p0/(Rd*T0)     ! kg/m3, air density at pressure p0 and temperature T0

    ! k-index of first layer above topography 
    do i=1,im
      do j=1,jm
        do k=1,km
          if (zc(k).ge.zsa(i,j)) then
            k1(i,j) = k
            exit
          endif
        enddo
      enddo
    enddo

    ! k-index for vertical velocity for clouds
    do i=1,im
      do j=1,jm
        do k=1,kmc
          if (zl(k).gt.max(zsa(i,j),hcld_base)) then
            kweff(i,j) = k
            exit
          endif
        enddo
      enddo
    enddo

    ! surface air density and pressure, function only of elevation
    do i=1,im
      do j=1,jm
        psa(i,j) = p0*exp(-zsa(i,j)/hatm)
        ra2a(i,j) = 0._wp
        do n=1,nm       
          ra2(i,j,n) = ra*exp(-zs(i,j,n)/hatm)
          ra2a(i,j) = ra2a(i,j)+frst(i,j,n)*ra2(i,j,n)        
          ps(i,j,n) = p0*exp(-zs(i,j,n)/hatm)
        enddo
      enddo
    enddo

    do i=1,im
      do j=1,jm

        imi=i-1
        if (imi.lt.1) imi=im 

        px = (pzsa(i,j)+pzsa(imi,j))*0.5_wp

        plx(i,j) = 0._wp
        plx_trop(i,j) = 0._wp

        do k=1,km
          if (px.le.pl(k+1))then
            dplx(i,j,k) = 0._wp
          elseif (px.lt.pl(k)) then
            dplx(i,j,k) = (px-pl(k+1))*amas
          else
            dplx(i,j,k) = (pl(k)-pl(k+1))*amas
          endif
          plx(i,j) = plx(i,j)+dplx(i,j,k)
          ! tropospheric column mass (diffusion is limited to k<=km-2, see adifa)
          if (k.le.km-2) plx_trop(i,j) = plx_trop(i,j)+dplx(i,j,k)
        enddo

      enddo
    enddo    
    ! periodic wrap point, so that plx can be addressed at i+1 up to imc
    plx(imc,:)      = plx(1,:)
    plx_trop(imc,:) = plx_trop(1,:)


    do i=1,im
      do j=2,jm

        py = (pzsa(i,j)+pzsa(i,j-1))*0.5_wp

        ply(i,j) = 0._wp
        ply_trop(i,j) = 0._wp

        do k=1,km
          if (py.le.pl(k+1))then
            dply(i,j,k) = 0._wp
          elseif (py.lt.pl(k)) then
            dply(i,j,k) = (py-pl(k+1))*amas
          else
            dply(i,j,k) = (pl(k)-pl(k+1))*amas
          endif
          ply(i,j) = ply(i,j)+dply(i,j,k)
          ! tropospheric column mass (diffusion is limited to k<=km-2, see adifa)
          if (k.le.km-2) ply_trop(i,j) = ply_trop(i,j)+dply(i,j,k)
        enddo

      enddo
    enddo
    dply(:,1,:) = 0._wp
    ! no flux through the poles, so these carry no mass
    ply(:,1)        = 0._wp
    ply(:,jmc)      = 0._wp
    ply_trop(:,1)   = 0._wp
    ply_trop(:,jmc) = 0._wp

    return

  end subroutine atm_grid_update

end module atm_grid
