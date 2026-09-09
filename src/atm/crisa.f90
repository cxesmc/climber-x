!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c r i s a _ m o d
!
!  Purpose : computation of cross-isobar angle
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
module crisa_mod

  use atm_params, only : wp
  use constants, only : pi, karman, omega
  use atm_params, only : cd0_ocn, cd0_sic, acbar_max, acbar_scale, nsmooth_acbar
  use atm_params, only : i_acbar, c_acbar_0, c_acbar_f, c_acbar_wind
  use atm_grid, only : im, jm, nm, fcorta_sqrt, fcorta, i_ocn, i_sic, i_lake
  use smooth_atm_mod, only : smooth2
  !$ use omp_lib

  implicit none

  ! |f| at 45 degrees, the reference Coriolis parameter of the i_acbar=2 closure.  Written with
  ! the numeric value of sin(pi/4) because sin() is not allowed in a constant expression.
  real(wp), parameter :: fcor45 = 2._wp*omega*0.70710678118654752_wp

  private
  public :: crisa

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  c r i s a
  !   Purpose    :  computation of cross-isobar angle
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine crisa(frst, z0m, zoro, &
        cd, cda, cd0, cd0a, acbar, sin_cos_acbar, cos_acbar, sin_acbar, epsa)

    implicit none

    real(wp), intent(in )  :: frst(:,:,:)
    real(wp), intent(in )  :: z0m(:,:,:)
    real(wp), intent(in )  :: zoro(:,:)

    real(wp), intent(out) :: cd(:,:,:)
    real(wp), intent(out) :: cda(:,:)
    real(wp), intent(out) :: cd0(:,:,:)
    real(wp), intent(out) :: cd0a(:,:)
    real(wp), intent(out) :: acbar(:,:)
    real(wp), intent(out) :: sin_cos_acbar(:,:)
    real(wp), intent(out) :: cos_acbar(:,:,:)
    real(wp), intent(out) :: sin_acbar(:,:,:)
    real(wp), intent(out) :: epsa(:,:,:)

    integer :: i, j, n
    real(wp) :: alfa, acbarn, acbarw

    real(wp), parameter :: z_ref = 100._wp      ! m, reference height

    !
    ! TWO ANGLES, because the one angle has two jobs and the data says they want different
    ! numbers.  The closure angle a solves sin(a)/sqrt(1-sin(2a)) = R and is what slp.f90 needs:
    ! it enters there only through sin(a)cos(a), which converts the ageostrophic mass flux of the
    ! mean meridional circulation into a zonal mean SLP gradient.  usur uses the same a twice, to
    ! ROTATE the geostrophic wind and to SCALE it by epsa = cos(a)-sin(a).  The rotation is what
    ! c_acbar_wind acts on; epsa keeps the closure angle, so that
    !     sqrt(us**2+vs**2) = epsa*|Vg|
    ! is EXACTLY unchanged and the scalar wind, all turbulent fluxes and the magnitude of the
    ! surface stress do not move at all.  Only the direction does.
    !
    ! The two jobs are measured separately and they disagree.  Because epsa cancels in the ratio,
    !     atan(signf*<vs>/<us>) = a
    ! is a measurement of the turning angle free of any amplitude assumption, and on the CMIP5
    ! aquaControl ensemble it gives 22-29 deg through the trades and 6-9 deg in the midlatitudes
    ! where the same ensemble's SLP-implied angle - which i_acbar=2 reproduces to 0.5 deg - is
    ! 12-15 and 5-7.  The factor between them is 1.6, not an offset: the demanded increment is
    ! +0.9 deg where a = 7 deg and +8.6 deg where a = 13 deg.  Fitting one factor to the zonal
    ! mean 10 m wind gives 1.70 (aquaplanet), 1.55 (PI) and 1.65 (LGM) independently, with a flat
    ! optimum over 1.4-1.9.
    !
    ! There is a second, independent reason to expect a factor above one, which also says where
    ! the missing wind went.  usur and slp.f90 share the angle but not the amplitude, so
    !     <vs>/vsz = epsa/cos(a) = 1 - tan(a),
    ! i.e. the surface meridional wind is 32 per cent weaker in the trades than the ageostrophic
    ! wind the SLP field was built from.  The model's own vabz is within 20 per cent of the
    ! observed v10 there while <vs> is 30-35 per cent too weak: the meridional flow is already
    ! right, usur discards a third of it.  A slab Ekman balance - the same one that produces
    ! sin(a)cos(a) in slp.f90 - gives the surface wind as cos(a)*Rot(a)*Vg, i.e. epsa = cos(a)
    ! rather than cos(a)-sin(a), and then vs would equal vsz identically.  That is NOT done here,
    ! because it would also multiply us by 1/(1-tan(a)) = 1.48 in the trades, which us cannot
    ! afford; the rotation delivers the same gain on vs (sin(1.6a)/sin(a) = 1.56 at a = 18 deg)
    ! while leaving us nearly alone (cos(1.6a)/cos(a) = 0.92).
    !
    ! What this does NOT fix is the zonal surface wind.  Rotation leaves its rms error flat
    ! (1.18 -> 1.19 m/s in PI); the aquaplanet has us too strong at 42.5-47.5 and too weak at
    ! 32.5 and poleward of 57.5, i.e. a surface westerly belt that is too narrow and too peaked.
    ! That is an error in the zonal mean SLP gradient feeding ugb, not in the partition.
    !
    ! See benchmark/out/ps_eval/diag_acbar_offset.jl, diag_acbar_eps.jl and diag_acbar_rot.jl.
    !

    if (i_acbar.ne.1 .and. i_acbar.ne.2) stop 'i_acbar'

    !$omp parallel do collapse(2) private(i, j, n, alfa, acbarn, acbarw)
    do j=1,jm
      do  i=1,im

        !------------------------------------------------
        ! drag coefficient 
        !------------------------------------------------

        do n=1,nm

          !------------------------------------------------
          ! neutral drag coefficient

          if (n.eq.i_ocn) then
            ! ocean

            cd(i,j,n) = cd0_ocn
            cd0(i,j,n) = cd(i,j,n)

          else if (n.eq.i_sic) then
            ! sea ice

            cd(i,j,n) = cd0_sic
            cd0(i,j,n) = cd(i,j,n)

          else if (n.eq.i_lake) then
            ! lake

            cd(i,j,n) = cd0_ocn
            cd0(i,j,n) = cd(i,j,n)

          else
            ! land or ice sheet

            if (frst(i,j,n).gt.0._wp) then
              cd(i,j,n) = (karman/log(z_ref/(z0m(i,j,n)+zoro(i,j))))**2
              cd0(i,j,n)   = (karman/log(z_ref/z0m(i,j,n)))**2
            else
              cd(i,j,n) = 0.01_wp
              cd0(i,j,n) = 0.01_wp
            endif

          endif

        enddo

        ! grid cell average
        cda(i,j) = sum(cd(i,j,:)*frst(i,j,:)) 

        ! drag coefficient without orographic component
        cd0a(i,j) = sum(cd0(i,j,:)*frst(i,j,:)) 


        !------------------------------------------------
        ! solve of the equation for cross-isobar angle   
        !------------------------------------------------

        alfa = acbar_cd(cda(i,j), j)
        acbar(i,j) = alfa
        sin_cos_acbar(i,j) = sin(alfa)*cos(alfa)

        ! Solving of the equation for cross-isobar angle for each surface type
        do n=1,nm
          acbarn = acbar_cd(cd(i,j,n), j)
          epsa(i,j,n) = sqrt(1._wp-sin(2._wp*acbarn))
          ! the angle usur rotates the geostrophic wind by, see c_acbar_wind below
          acbarw = c_acbar_wind*acbarn
          cos_acbar(i,j,n) = cos(acbarw)
          sin_acbar(i,j,n) = sin(acbarw)

        enddo

      enddo
    enddo
    !$omp end parallel do

    ! smooth in space
    call smooth2(acbar,nsmooth_acbar)
    call smooth2(sin_cos_acbar,nsmooth_acbar)

    return

  end subroutine crisa


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  a c b a r _ c d
  !   Purpose    :  cross-isobar angle belonging to one drag coefficient
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  ! Both options solve the same relation,
  !   sin(a)/sqrt(1-sin(2a)) = R,
  ! and differ only in R.  Since 1-sin(2a) = (cos(a)-sin(a))**2 the relation is
  ! tan(a) = R/(1+R) in closed form.
  !
  ! i_acbar = 1  The original,  R = cd/sqrt(|f|),  solved by bisection on [0,pi/4].  The
  !              bisection is therefore solving something that has a closed form; it is kept
  !              exactly as it was, iteration included, so that earlier runs are reproduced bit
  !              for bit.  R is NOT dimensionally homogeneous: cd is dimensionless while f is
  !              1/s, so R carries s**(1/2) and a constant of 1 s**(-1/2) is hidden in the
  !              calibration.
  !
  ! i_acbar = 2  R = cd/(c_acbar_0 + c_acbar_f*phi**2),  phi = |f|/f_45,  in closed form.
  !
  !              Fitted to the angle IMPLIED by the relation the angle is actually used in.
  !              slp.f90 defines it through
  !                psz(j) = psz(j-1) + vsz*fcorua*ra*dy/sin(a)cos(a),
  !              which inverts to sin(a)cos(a) = vsz/(signf*ug)
  !              with ug from the zonal mean SLP gradient.  Measured that way on the CMIP5
  !              aquaControl ensemble, the angle is
  !                             5 deg   10    15    20    30    35    40    45    60
  !                measured      15.5   14.4  12.5  11.5   7.1   6.6   5.4   4.8   4.6
  !                i_acbar=1     15.0   11.6   9.9   8.8   7.5   7.1   6.8   6.5   5.9
  !              i.e. the real angle has considerably MORE latitudinal contrast than cd/sqrt(f)
  !              carries.  This form reproduces it to 0.46 deg rms and 0.92 deg at worst,
  !              against 1.86 and 2.85 for i_acbar=1, over the 16 rows where the inversion is
  !              well conditioned.  The hemispheric asymmetry of the measurement itself, on a
  !              symmetric aquaplanet, is 0.08 deg, so the residual is well above the noise but
  !              the improvement is far larger than it.
  !
  !              A power law in f cannot do this: the log-log slope of R against |f| steepens
  !              from -0.28 in the tropics to -1.45 in the subtropics, and the best power law
  !              lands at 1.54 deg rms, barely better than the original.  What the measured R
  !              does follow is 1/R linear in f**2, which is this form.
  !
  !              Two structural gains over i_acbar=1.  It is dimensionally homogeneous, both
  !              constants carrying the units of cd, i.e. none.  And it is bounded at the
  !              equator by construction, at cd/c_acbar_0, so the fcormin floor no longer has to
  !              hold the angle up there - with the shipped constants the ocean angle runs from
  !              15.9 deg at the equator (16.2 in the phi -> 0 limit, the floor costing 0.3) to
  !              3.0 deg at the pole.
  !
  !              acbar_max still binds over rough orography, as it does for i_acbar=1: cd = 0.02
  !              at 45 deg gives 28.7 deg here against 33.5 for the original, both on the cap.
  !
  !              The cd dependence is LINEAR, the same as i_acbar=1 carries.  That is not a free
  !              choice made here: regressing the demanded R on ln(cd) and on the azonal SLP
  !              amplitude over the PI and LGM ensembles gives an exponent of 0.89-0.96 with the
  !              two predictors only 0.16 correlated.  It is, however, only that - an exponent
  !              near one - and the PI/LGM data cannot pin the amplitude, because the northern
  !              extratropical inversion is not a measurement of a turning angle there (see the
  !              note on c_acbar_0 in atm_par.nml).
  !
  function acbar_cd(cd, j) result(alfa)

    implicit none

    real(wp), intent(in) :: cd
    integer,  intent(in) :: j
    real(wp) :: alfa

    integer :: iter
    real(wp) :: rhs, rhsn, alfa0, alfa1


    if (i_acbar.eq.1) then

      rhs = cd/fcorta_sqrt(j)
      alfa0 = 0._wp
      alfa1 = pi/4._wp
      ! Iteration loop
      do iter=1,10
        alfa = 0.5_wp*(alfa0+alfa1)
        rhsn = sin(alfa)/sqrt(1._wp-sin(2._wp*alfa))
        if (rhsn.gt.rhs) then
          alfa1 = alfa
        else
          alfa0 = alfa
        endif
      enddo

    else

      ! fcorta is the floored |f|, which is what the fit used; the floor changes phi**2 by less
      ! than 1 per cent of c_acbar_0 anyway, since this form does not need it
      rhs = cd/(c_acbar_0 + c_acbar_f*(fcorta(j)/fcor45)**2)
      alfa = atan(rhs/(1._wp+rhs))

    endif

    alfa = alfa*acbar_scale
    alfa = max(alfa,0.05_wp)
    alfa = min(alfa,acbar_max)

    return

  end function acbar_cd

end module crisa_mod
