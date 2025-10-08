!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : i n s o l a t i o n
!
!  Purpose : insolation at the top of the atmosphere
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Andrey Ganopolsky, Alexander Robinson and 
!                         Matteo Willeit
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
module insolation

   use precision, only : wp
   use ncio
   use constants, only: pi
   use timer, only: day_year, nday_year
   use control, only : iorbit, ecc_const, obl_const, per_const, l_aquaplanet

   implicit none

    real(wp) :: ecc, per, xob
    real(wp), dimension(:), allocatable :: time_dat, ecc_dat, per_dat, obl_dat

    ! present day orbital parameters 
!    real(wp), parameter :: btime = 0._wp
!    real(wp), parameter :: ecc = 0.1670236225492288d-1
!    real(wp), parameter :: xob = 0.4090928042223415_wp * 180._wp/pi ! deg
!    real(wp), parameter :: per = 0.1796256991128036d1 * 180._wp/pi ! deg


contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine : i n i _ s i n s o l
  ! Purpose    : read orbital parameters from Laskar
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ini_sinsol(orbit_file)

    implicit none

   character (len=*), intent(in) :: orbit_file

   integer :: ntime


   ntime = nc_size(trim(orbit_file),"time")
   allocate( time_dat(ntime) )
   allocate( ecc_dat(ntime) )
   allocate( per_dat(ntime) )
   allocate( obl_dat(ntime) )
   call nc_read(trim(orbit_file),"time",time_dat)    ! time in years BP
   call nc_read(trim(orbit_file),"ecc",ecc_dat) 
   call nc_read(trim(orbit_file),"per",per_dat) 
   call nc_read(trim(orbit_file),"obl",obl_dat) 

      return

  end subroutine ini_sinsol


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  s i n s o l
  ! Purpose    :  compute solar insolation from orbital parameters
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sinsol(s0, btime, lat, &
                   eccentricity, precession, obliquity, solar, solarm, solarmin, solarmax, cosz, coszm)

    implicit none

    real(wp), intent(in) :: s0      ! solar 'constant'
    real(wp), intent(in) :: btime   ! astronomical time, present day 1950 = 0.0
    real(wp), intent(in) :: lat(:)  ! array of latitude values, in radians!
    real(wp), intent(out) :: eccentricity
    real(wp), intent(out) :: precession
    real(wp), intent(out) :: obliquity
    real(wp), dimension(:,:,:), intent(out) :: solar   ! solar radiation at the top of the atmosphere (with diurnal cycle)
    real(wp), dimension(:,:), intent(out) :: solarm    ! daily mean solar radiation at the top of the atmosphere 
    real(wp), dimension(:,:), intent(out) :: solarmin  ! daily minimum solar radiation at the top of the atmosphere 
    real(wp), dimension(:,:), intent(out) :: solarmax  ! daily maximum solar radiation at the top of the atmosphere
    real(wp), dimension(:,:,:), intent(out) :: cosz    ! solar zenith angle (with diurnal cycle)
    real(wp), dimension(:,:), intent(out) :: coszm     ! daytime averaged solar zenith angle
    
    integer :: j, k, h, h1, ny
    real(wp) :: s, cosn, cosp, acosp, don, fi
    real(wp) :: tperi, zavexpe
    real(wp) :: pclock, pytime, pdisse, pzen1, pzen2, pzen3, prae
    integer, parameter :: nh = 24*10
    real(wp), allocatable :: sdom(:,:)


    ny = size(lat)

    allocate(sdom(nday_year,ny))

    ! get current orbital parameters
    if (iorbit.eq.0 .or. iorbit.eq.2) then
      call orbit_par(btime)
    else
      ecc = ecc_const
      per = per_const
      xob = obl_const
    endif

    ! settings for aquaplanet
    if (l_aquaplanet) then
      ecc = 0._wp
      xob = 0._wp
    endif

    eccentricity = ecc
    precession = sin(per*pi/180._wp)
    obliquity = xob

!...1) Berger program calculates orbital parameters
!   ===========================================   
    call bergor(tperi,zavexpe)
!   =========================================== 
    
!...2) Daily insolation is calculated by 1/2 hourly integration for each day
     
    do j = 1, ny
      fi = lat(j) * pi/180._wp
      do k = 1, nday_year
        pytime=k*2._wp*pi/day_year
        solarm(k,j)=0._wp
        solarmin(k,j)=1.e6_wp
        solarmax(k,j)=0._wp
        coszm(k,j) =0._wp        
        sdom(k,j) =0._wp    
        h1 = 0
        do h = 1, nh
          pclock=h*2._wp*pi/nh  
          !   =================================================================          
          call orbit(tperi,zavexpe,pclock,pytime,pdisse,pzen1,pzen2,pzen3,prae)
          !   =================================================================                    
          cosp=pzen1*sin(fi)+pzen2*cos(fi)
          cosn=max(cosp,0.0_wp)
          cosp=max(cosp,0.05_wp)
          acosp=1._wp/cosp
          acosp=min(acosp,10._wp)
          if (cosn.eq.0._wp) then
            don=0._wp
          else
            don=1._wp
          endif

          s=s0*cosn*pdisse
          if (modulo(h,nh/24/2+h1*nh/24).eq.0) then
            h1 = h1+1
            solar(k,h1,j) = s
            cosz(k,h1,j) = cosn
          endif
          solarm(k,j)=solarm(k,j)+s
          if (s<solarmin(k,j)) solarmin(k,j) = s
          if (s>solarmax(k,j)) solarmax(k,j) = s          
          coszm(k,j)=coszm(k,j)+cosn
          sdom(k,j)=sdom(k,j)+don
        enddo
        if (solarmin(k,j)>solarmax(k,j)) solarmin(k,j) = 0._wp
      enddo

      ! Daily insolation and zenite angle
      do k = 1, nday_year
        solarm(k,j)=solarm(k,j)/nh
        if (sdom(k,j) .gt. 0._wp) then
          coszm(k,j)=coszm(k,j)/sdom(k,j)
        else
          coszm(k,j)=0.0_wp
        endif
      enddo

    enddo ! end j-loop

    deallocate(sdom)

    return
  
  end subroutine sinsol 

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine : o r b i t
  ! Purpose    : 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine orbit(tperi,zavexpe,pclock,pytime, &
                  pdisse,pzen1,pzen2,pzen3,prae)


  implicit none

      real(wp), intent(in) :: tperi, zavexpe, pclock, pytime
      real(wp), intent(out) :: pdisse, pzen1, pzen2, pzen3, prae

      integer :: niter
      real(wp) :: zclock, zytime, ttrop, time, eold, enew, zeps, cose
      real(wp) :: E, zdisse, zsqecc, ztgean, znu, zlambda, xobche
      real(wp) :: zsinde, zdecli, zzen1, zzen2, zzen3
      real(wp), parameter :: eps = 1.d-6
      real(wp), parameter :: zrae = 0.1277d-2


      zclock=pclock
      zytime=pytime

!     TROPICAL YEAR

      ttrop=360._wp
!
!     EXACT CALCULATION WITH KEPLER'S LAW (ELLIPSE)
!     (MONIN, 1986 : AN INTRODUCTION TO THE THEORY OF CLIMATE)
!     USE FIRST DERIVATIVE (NEWTON) TO CALCULATE SOLUTION OF
!     EQUATION OF eccENTRIC ANOMALY *E*
!
      time=zytime-2._wp*pi*tperi/ttrop
      eold=time/(1._wp-ecc)
      enew=time
      niter=0
  250 continue
      zeps=eold-enew
      if (niter.ge.30) go to 270
      if (abs(zeps).lt.eps) go to 280
      niter=niter+1
      eold=enew
      cose=cos(enew)
      enew=(time+ecc*(sin(enew)-enew*cose))/(1._wp-ecc*cose)
      go to 250
  270 print*,' SUBROUTINE *orbit* - eccentric anomaly not found!'
      print*,' ERROR IN   *orbit* -- STOP'
      stop
  280 continue
      E=enew
      zdisse=(1._wp/(1._wp-ecc*COS(E)))**2
!
!     Change in the calculation of the declination.
!     Used are not the formulas from Paltridge/Platt as in the ECHAM model
!     with fixed constants for contemporanious conditions
!     but the exact equations from Monin (s.a. - Keplers law)
!     are solved. Day of vernal equinox is fixed for a 360 day year on the
!     21. Maerz noon, with start of year at 01/01/00 GMT (s.a.).
!
      zsqecc=sqrt((1._wp+ecc)/(1._wp-ecc))
      ztgean=tan(E/2._wp)
!
!     znu: true anomaly (actual angle of Earth's position from perihelion)
!     zlambda: true longitude of the Earth (actual angle from vernal equinox)
!
      znu=2._wp*atan(zsqecc*ztgean)
      zlambda=znu+zavexpe
      !xobche=xobch/180.*pi
      xobche=xob/180._wp*pi
      zsinde=sin(xobche)*sin(zlambda)
      zdecli=asin(zsinde)
 
!
      zzen1=SIN(zdecli)
      zzen2=COS(zdecli)*COS(zclock)
      zzen3=COS(zdecli)*SIN(zclock)
!
!
      pdisse=zdisse
      pzen1=zzen1
      pzen2=zzen2
      pzen3=zzen3
      prae=zrae
 
      return

   end subroutine orbit


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine : b e r g o r
  ! Purpose    : calculate time of perihelion in year:
  !              An Introdiuction to the Theory of Climate (Reidel Publ.)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   subroutine bergor(tper,zan)

   implicit none

     real(wp), intent(out) :: tper
     real(wp), intent(out) :: zan

     real(wp) :: ang, tgnu, sqecc, e, tian, days


      ! - per: perihelion from aut. equinox in degrees (BERGER)
      ! - ecc: eccentricity of Earth's orbit (BERGER)

      ang=per-180._wp
      if (ang.lt.0._wp) ang=ang+360._wp
      zan=ang/180._wp*pi                !   =zavexpe
      tgnu=tan(zan/2._wp)
      sqecc=sqrt((1._wp-ecc)/(1._wp+ecc))
      e=2._wp*atan(sqecc*tgnu)          !   Eccentric Anomaly

!   time angle in radians of perihelion from vernal equinox

      tian=e-ecc*sin(e)
      days=tian*180._wp/pi                !   360 day year only
      if (days.lt.0._wp) days=days+360._wp !   days from ver.eq. to perh.

!   time in days from begin of year: vernal equinox fixed at 3/21, 12 GMT
!    = 80.5 days in 360 day year

      tper=days+80.5_wp
      if (tper.gt.360._wp) tper=tper-360._wp

            
      return

    end subroutine bergor


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine : o r b i t _ p a r
  ! Purpose    : read orbital parameters from Laskar
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! input : btime- Time in years (negative for the past, reference year 1950 a. D.)
  ! output: Earths orbital parameters
  !         per  - Longitude of Perihelion (measured from vernal equinox,
  !                relative to observation from the earth, i.e. angle from
  !                autumnal equinox to perihelion)
  !         ecc  - Eccentricity
  !         xob  - Obliquity
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine orbit_par(btime)

    implicit none

    real(wp), intent(in) :: btime
    real(wp) :: per_dat1, per_dat2
    integer :: imin, i0, i1
    real(wp) :: w0, w1


      ! interpolate to current year
      if (btime.le.minval(time_dat)) then
        ecc = ecc_dat(lbound(ecc_dat,1))
        xob = obl_dat(lbound(obl_dat,1))
        per = per_dat(lbound(per_dat,1))
      elseif (btime.ge.maxval(time_dat)) then
        ecc = ecc_dat(ubound(ecc_dat,1))
        xob = obl_dat(ubound(obl_dat,1))
        per = per_dat(ubound(per_dat,1))
      else
        imin = minloc(abs(time_dat-btime),1) 
        if (time_dat(imin).lt.btime) then
          i0 = imin
          i1 = imin+1
        else
          i0 = imin-1
          i1 = imin
        endif
        w0 = 1._wp - abs(time_dat(i0)-btime)/(time_dat(i1)-time_dat(i0))
        w1 = 1._wp - w0
        ecc = w0*ecc_dat(i0) + w1*ecc_dat(i1)
        xob = w0*obl_dat(i0) + w1*obl_dat(i1)
        per_dat1=per_dat(i0)
        per_dat2=per_dat(i1)
        if (per_dat2.lt.per_dat1) per_dat2=per_dat2+360._wp
        per = w0*per_dat1 + w1*per_dat2
      endif
      if (per.gt.360._wp) per=per-360._wp

    return

   end subroutine orbit_par


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  d a y _ l e n g t h
  !   Purpose    :  computes daylength for each latitude and day of year
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine day_length(lat, &
                       daylength)

    implicit none

    real(wp), intent(in) :: lat(:)
    real(wp), dimension(:,:), intent(out) :: daylength

    integer :: j, k, ny
    real(wp) :: delta, u, v, hh


    ny = size(lat)

    do j=1,ny
      do k=1,nday_year

        delta=pi/180._wp*(-xob*cos(2._wp*pi*(k+10.0_wp)/day_year));
        u=sin(pi/180._wp*(lat(j)))*sin(delta);
        v=cos(pi/180._wp*(lat(j)))*cos(delta);

        if(u .ge. v) then
          daylength(k,j)=24._wp
          elseif(u .le. -v) then
          daylength(k,j)=0._wp
        else
          hh=acos(-u/v)
          daylength(k,j)=24._wp*hh/pi
        endif

      enddo
    enddo

    return

  end subroutine day_length

end module insolation
