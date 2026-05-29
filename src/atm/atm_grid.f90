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
  use atm_params, only : atm_mass, hatm, amas, ra, hcld_base, fcormin, fcoramin, fcorumin
  use atm_params, only : l_p0_var, p0, ps0, pble, pblp
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
  integer, parameter :: lm = 3  !! number of water tracers, incl. oxygen isotopes
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
  real(wp) :: fcortf(jm)
  real(wp) :: fcorta_sqrt(jm)
  real(wp) :: fcoru(jmc)
  real(wp) :: fcorta(jm)
  real(wp) :: fcorua(jmc)
       
  real(wp) :: plx(imc,jm)
  real(wp) :: ply(im,jmc)
  real(wp) :: pblt(jm)      !! planetary boundary height in t-points
  real(wp) :: pblu(jmc)     !! planetary boundary height in u-points
  integer :: k1(im,jm)  
  integer :: kweff(im,jm)  !! k-index for effective vertical velocity leve (w-grid)
  integer :: k1000
  integer :: k900
  integer :: k850
  integer :: k700
  integer :: k500
  integer :: k300

  real(wp), allocatable :: pl(:)
  real(wp), allocatable :: dpl(:)
  real(wp), allocatable :: zl(:)
  real(wp), allocatable :: zc(:)
  real(wp), allocatable :: dplx(:,:,:)
  real(wp), allocatable :: dply(:,:,:)
  real(wp), allocatable :: dplxo(:,:,:)
  real(wp), allocatable :: dplyo(:,:,:)
  integer, allocatable :: kplxo(:,:)
  integer, allocatable :: kplyo(:,:)
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
    allocate(dplxo(im,jm,km))
    allocate(dplyo(im,jm,km))
    allocate(kplxo(im,jm))
    allocate(kplyo(im,jm))
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

    do j=1,jm
      fcortp = 2._wp*omega*sint(j) 
      fcort(j) = signf(j)*max(ABS(fcortp),fcormin)
      fcortf(j) = signf(j)*max(ABS(fcortp),fcorumin)
      fcorta(j) = max(ABS(fcortp),fcoramin)        
      fcorta_sqrt(j) = sqrt(abs(fcorta(j)))
    enddo

    do j=1,jm
      fcorup = 2._wp*omega*sinu(j) 
      fcoru(j) = signf(j)*max(ABS(fcorup),fcormin)
      fcorua(j) = max(ABS(fcorup),fcoramin)
    enddo

    ! PBL thickness
    do j=1,jm
      pblt(j) = pblp-(pblp-pble)*cost(j)**2      
      pblu(j) = pblp-(pblp-pble)*cosu(j)**2 
    enddo       

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

    print *,'k 1000 hPa',k1000, ', z 1000 hPa',zl(k1000)
    print *,'k 900  hPa',k900,  ', z 900  hPa',zl(k900)
    print *,'k 850  hPa',k850,  ', z 850  hPa',zl(k850)
    print *,'k 700  hPa',k700,  ', z 700  hPa',zl(k700)
    print *,'k 500  hPa',k500,  ', z 500  hPa',zl(k500)
    print *,'k 300  hPa',k300,  ', z 500  hPa',zl(k300)

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

        do k=1,km         
          if (px.le.pl(k+1))then
            dplx(i,j,k) = 0._wp
          elseif (px.lt.pl(k)) then
            dplx(i,j,k) = (px-pl(k+1))*amas
          else
            dplx(i,j,k) = (pl(k)-pl(k+1))*amas
          endif
          plx(i,j) = plx(i,j)+dplx(i,j,k)         
          ! orographic component
          dplxo(i,j,k) = max(0._wp,min(pl(k),max(pzsa(imi,j),pzsa(i,j)))-pl(k+1))*amas
        enddo

      enddo
    enddo    


    do i=1,im
      do j=2,jm

        py = (pzsa(i,j)+pzsa(i,j-1))*0.5_wp

        ply(i,j) = 0._wp

        do k=1,km         
          if (py.le.pl(k+1))then
            dply(i,j,k) = 0._wp
          elseif (py.lt.pl(k)) then
            dply(i,j,k) = (py-pl(k+1))*amas
          else
            dply(i,j,k) = (pl(k)-pl(k+1))*amas
          endif
          ply(i,j) = ply(i,j)+dply(i,j,k)         
          ! orographic component
          dplyo(i,j,k) = max(0._wp,min(pl(k),max(pzsa(i,j-1),pzsa(i,j)))-pl(k+1))*amas
        enddo

      enddo
    enddo
    dply(:,1,:) = 0._wp
    dplyo(:,1,:) = 0._wp

    ! k-index of first layer above topography
    do i=1,im
      do j=1,jm
        do k=1,km
          if (dplx(i,j,k).gt.0._wp) then  ! first layer above topography
            kplxo(i,j) = k
            exit
          endif
        enddo
      enddo
    enddo
    do i=1,im
      do j=2,jm
        do k=1,km
          if (dply(i,j,k).gt.0._wp) then  ! first layer above topography
            kplyo(i,j) = k
            exit
          endif
        enddo
      enddo
    enddo
    kplyo(:,1) = 1

    return

  end subroutine atm_grid_update

end module atm_grid
