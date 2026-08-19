!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : v e s t a _ m o d
!
!  Purpose : vertical structure of atmosphere
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
module vesta_mod

  use atm_params, only : wp
  use constants, only : fqsat, pi
  use atm_params, only : gad, hatm, p0, ra, zmax
  use atm_params, only : c_gam_1, c_gam_2, c_gam_3, hgams, hgamt, c_gam_rel, nsmooth_gam
  use atm_params, only : gams_min, gams_max, sh_gams
  use atm_params, only : c_hrs_1, c_hrs_2, c_hrs_3, c_hrs_4, c_hrs_5, c_hrs_6, rh_strat
  use atm_params, only : i_zpbl, h_pbl_min
  use atm_params, only : i_rh_free, rh_free, c_rhf_1, c_rhf_2, c_rhf_3, c_rhf_4, rhf_min, rhf_max
  use atm_params, only : c_dhs_1, c_dhs_2
  use atm_params, only : c_trop_1, c_trop_2, c_trop_3
  use atm_params, only : l_dust
  use atm_grid, only : im, jm, km, aim, zl, fit, exp_zc
  use smooth_atm_mod, only : smooth2
  !$ use omp_lib

  implicit none

  private
  public :: hscales, vesta, t_prof, rh_prof, tropoheight

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  h s c a l e s
  !   Purpose    :  computation of lapse rate and height scales of moisture and dust
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine hscales(ra2a, sha, qam, wcon, wcld, had_fi, had_width, &
      gams, gamb, gamt, hrm, rhfree, &
      hqeff, hdust)

    implicit none

    real(wp), intent(in) :: ra2a(:,:)
    real(wp), intent(in) :: sha(:,:)   !! grid mean sensible heat flux, positive upward (W/m2)
    real(wp), intent(in) :: qam(:,:)
    real(wp), intent(in) :: wcon(:,:)
    real(wp), intent(in) :: wcld(:,:)
    real(wp), intent(in) :: had_fi
    real(wp), intent(in) :: had_width

    real(wp), intent(inout) :: gams(:,:)
    real(wp), intent(inout) :: gamb(:,:)
    real(wp), intent(inout) :: gamt(:,:)
    real(wp), intent(inout) :: hrm(:,:)
    real(wp), intent(inout) :: rhfree(:,:)

    real(wp), intent(out) :: hqeff(:,:)
    real(wp), intent(out) :: hdust(:,:)

    integer :: i, j
    real(wp) :: hrs, fi, f_trop, fi_rhf, f_rhf, rhfs
    real(wp), dimension(im,jm) :: gam_s, gam_b, gam_t

    ! The floor on hrm was 1000 m, which is above the values the observed profiles ask for
    ! once rh_prof relaxes towards a background free tropospheric humidity (i_rh_free): the
    ! best fit to ERA-Interim and CMIP5 is 300-2000 m by band, geometric mean about 730 m.
    ! With the floor at 1000 m the residual (ram-rh_free)*exp(-2000/hrm) term still overshoots
    ! the observed free troposphere by 0.04-0.06, which would have to be absorbed by tuning
    ! rh_free downwards. At 500 m that residual is 0.003 and rh_free can be set to the
    ! observed free tropospheric relative humidity directly.
    real(wp), parameter :: hrs_min = 500._wp
    real(wp), parameter :: hrs_max = 10.e3_wp


    !$omp parallel do collapse(2) private(i,j,hrs,fi,f_trop,fi_rhf,f_rhf,rhfs)
    do j=1,jm
      do i=1,im

        !----------------------------------------------
        ! lapse rate

        ! lapse rate in the boundary layer.
        gam_s(i,j) = gams_min + (gams_max-gams_min)*0.5_wp*(1._wp+tanh(sha(i,j)/sh_gams))

        ! bottom
        gam_b(i,j) = c_gam_1 - c_gam_2*qam(i,j) 

        ! top
        gam_t(i,j) = gam_b(i,j) + c_gam_3

        !----------------------------------------------
        ! height scale for relative humidity       

        fi = c_hrs_6*(fit(j)-had_fi)/(0.5_wp*had_width)
        fi = min(fi,pi/2._wp)
        fi = max(fi,-pi/2._wp)
        f_trop = 1._wp-sin(fi)**8
        hrs = f_trop * c_hrs_1*exp(c_hrs_2*wcld(i,j)) + (1._wp-f_trop) * c_hrs_1*c_hrs_3 
        hrs = max(hrs,hrs_min)   
        hrs = min(hrs,hrs_max)   
        hrm(i,j) = 0.9_wp*hrm(i,j) + 0.1_wp*hrs

        !----------------------------------------------
        ! background free tropospheric relative humidity
        ! The profile in rh_prof relaxes towards this value above the boundary layer instead of
        ! decaying to zero. Observations put it at 0.68 poleward of 60 deg, 0.56 at 40-60 deg,
        ! 0.37 at 20-40 deg and about 0.44 over the ocean in the deep tropics, i.e. a minimum in
        ! the subsidence belts and maxima in the ascending tropics and at the poles. That is the
        ! same shape the Hadley weight and the cloud level vertical velocity wcld already
        ! describe, so the same construction as for hrs is used, with its own width parameter
        ! c_rhf_4 because the free tropospheric humidity leaves the tropical regime much further
        ! equatorward than hrs does, and it is relaxed in time in the same way.
        if (i_rh_free.eq.0) then
          ! uniform, taken straight from the namelist and not relaxed or clamped, so that
          ! rh_free = 0 reproduces the original purely multiplicative profile exactly
          rhfree(i,j) = rh_free
        else
          fi_rhf = c_rhf_4*(fit(j)-had_fi)/(0.5_wp*had_width)
          fi_rhf = min(fi_rhf,pi/2._wp)
          fi_rhf = max(fi_rhf,-pi/2._wp)
          f_rhf  = 1._wp-sin(fi_rhf)**8
          rhfs = f_rhf * (c_rhf_1 + c_rhf_2*wcld(i,j)) + (1._wp-f_rhf) * c_rhf_3
          rhfs = max(rhfs,rhf_min)
          rhfs = min(rhfs,rhf_max)
          rhfree(i,j) = 0.9_wp*rhfree(i,j) + 0.1_wp*rhfs
        endif

        !----------------------------------------------
        ! effective moisture height scale

        hqeff(i,j) = wcon(i,j)/(ra2a(i,j)*qam(i,j))

        !----------------------------------------------
        ! dust height scale 

        hdust(i,j) = c_dhs_1+c_dhs_2*wcld(i,j)

      enddo
    enddo
    !$omp end parallel do

    !----------------------------------------------
    ! smoothing and time relaxation of lapse rate

    !$omp parallel sections
    !$omp section
    call smooth2(gam_s,nsmooth_gam)
    do j=1,jm
      do i=1,im
        gams(i,j) = c_gam_rel*gams(i,j) + (1._wp-c_gam_rel)*gam_s(i,j)
      enddo
    enddo
    !$omp section
    call smooth2(gam_b,nsmooth_gam)
    do j=1,jm
      do i=1,im
        gamb(i,j) = c_gam_rel*gamb(i,j) + (1._wp-c_gam_rel)*gam_b(i,j)
      enddo
    enddo
    !$omp section
    call smooth2(gam_t,nsmooth_gam)
    do j=1,jm
      do i=1,im
        gamt(i,j) = c_gam_rel*gamt(i,j) + (1._wp-c_gam_rel)*gam_t(i,j)
      enddo
    enddo
    !$omp end parallel sections

    return

  end subroutine hscales


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  v e s t a
  !   Purpose    :  vertical structure of atmosphere 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine vesta(zsa, tam, gams, gamb, gamt, htrop, ram, hrm, rhfree, dam, hdust, &
      A_trop, W_strat, t3, q3, tp, d3, ttrop)

    implicit none

    real(wp), intent(in ) :: zsa(:,:)
    real(wp), intent(in ) :: tam(:,:)
    real(wp), intent(in ) :: gams(:,:)
    real(wp), intent(in ) :: gamb(:,:)
    real(wp), intent(in ) :: gamt(:,:)
    real(wp), intent(in ) :: htrop(:,:)
    real(wp), intent(in ) :: ram(:,:)
    real(wp), intent(in ) :: hrm(:,:)
    real(wp), intent(in ) :: rhfree(:,:)
    real(wp), intent(in ) :: dam(:,:)
    real(wp), intent(in ) :: hdust(:,:)

    real(wp), intent(out) :: A_trop(:,:)   
    real(wp), intent(out) :: W_strat(:,:)  
    real(wp), intent(out) :: t3(:,:,:)
    real(wp), intent(out) :: q3(:,:,:)
    real(wp), intent(out) :: tp(:,:,:)
    real(wp), intent(out) :: d3(:,:,:)
    real(wp), intent(out) :: ttrop(:,:)

    integer :: i, j, k
    logical :: flag_strat
    real(wp) :: z_sur, taml, htropl
    real(wp) :: gamsl, gambl, gamtl, z
    real(wp) :: t, rsur, hrml, rhfl, rh, rh_unit, rh_zero, q, q_unit, q_zero, qsat, A_l, W_l
    real(wp) :: dvol


    !$omp parallel do collapse(2) private(i,j,k,z_sur,taml,htropl,gamsl,gambl,gamtl,z,t,rsur,hrml,rhfl,rh,rh_unit,rh_zero,q,q_unit,q_zero,qsat,A_l,W_l,dvol,flag_strat)
    do j=1,jm
      do i=1,im

        ! 2D fields

        z_sur  = zsa(i,j)
        taml   = tam(i,j)
        gamsl  = gams(i,j)
        gambl  = gamb(i,j)
        gamtl  = gamt(i,j)
        htropl = htrop(i,j)
        rsur   = ram(i,j)
        hrml   = hrm(i,j)
        rhfl   = rhfree(i,j)

        ! 3D fields of temperature and humidity

        A_l   = 0._wp   ! ram-dependent coefficient ∫ E(z)·qsat·ρ dz over the troposphere
        W_l   = 0._wp   ! ram-independent intercept: free tropospheric background + stratosphere
        flag_strat = .false.
        do k=1,km

          z = 0.5_wp*(zl(k)+zl(k+1))

          if (.not.flag_strat) then
            ! construct vertical temperature profile
            t = t_prof(z_sur, z, taml, gamsl, gambl, gamtl, htropl, 1)
            ! derive specific humidity profile from temperature and relative humidity profiles
            rh = rh_prof(z_sur, z, rsur, hrml, rhfl, htropl)
            ! rh_prof is AFFINE in ram: rh = ram·E(z) + rh_free·(1-E(z)) in the troposphere,
            ! with E(z) the exponential decay factor. Split it into the ram-dependent slope
            ! E(z) and the ram-independent intercept, so that the column water stays linear
            ! in ram and time_step can still invert wcon = ram·A_trop + W_strat.
            ! With rh_free = 0 this reduces to rh_zero = 0 and rh_unit = E(z), i.e. the
            ! original formulation.
            rh_zero = rh_prof(z_sur, z, 0._wp, hrml, rhfl, htropl)
            rh_unit = rh_prof(z_sur, z, 1._wp, hrml, rhfl, htropl) - rh_zero
            qsat = fqsat(t,p0*exp_zc(k))
            q = rh*qsat
            q_unit = rh_unit*qsat
            q_zero = rh_zero*qsat
          endif

          ! volume element (mass per unit area for this layer above surface)
          if (zl(k).ge.z_sur) then
            dvol = ra*exp_zc(k)*(zl(k+1)-zl(k))
          else if (zl(k).lt.z_sur .and. zl(k+1).gt.z_sur) then
            dvol = ra*exp_zc(k)*(zl(k+1)-z_sur)
          else
            dvol = 0._wp
          endif

          ! vertical integral of water content split into the ram-dependent slope (A) and the
          ! ram-independent intercept (W); the total column water wcon = ram·A_trop + W_strat
          ! is reconstructed in time_step.
          if (z.le.htropl+1._wp) then
            ! tropospheric: q = ram·E(z)·qsat + rh_free·(1-E(z))·qsat
            A_l = A_l + q_unit*dvol
            W_l = W_l + q_zero*dvol   ! zero unless rh_free > 0
          else
            ! stratospheric: q = rh_strat·qsat ⇒ dW = q·dvol (independent of ram)
            W_l = W_l + q*dvol
          endif

          t3(i,j,k) = t
          q3(i,j,k) = q

          ! potential temperature
          tp(i,j,k) = t + gad*min(z,zmax)

          ! dust profile
          if (l_dust) d3(i,j,k) = dam(i,j)*min(1._wp,exp(-(z-z_sur)/hdust(i,j))) ! kg/kg

          if (z.gt.htropl) flag_strat = .true.

        enddo

        A_trop(i,j)  = A_l
        W_strat(i,j) = W_l

        ! tropopause temperature

        ttrop(i,j) = t

      enddo
    enddo
    !$omp end parallel do

    return 

  end subroutine vesta


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  t _ p r o f
  !   Purpose    :  compute vertical temperature profile
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pure function t_prof(zs, z, tam, gams, gamb, gamt, htrop, iflag)

    implicit none

    real(wp), intent(in) :: zs
    real(wp), intent(in) :: z
    real(wp), intent(in) :: tam
    real(wp), intent(in) :: gams
    real(wp), intent(in) :: gamb
    real(wp), intent(in) :: gamt
    real(wp), intent(in) :: htrop
    integer, intent(in) :: iflag

    real(wp) :: t_prof

    real(wp) :: zk


    zk = min(z,htrop)

    if (iflag.eq.0) then
      ! temperature profile ignoring surface layer

      t_prof = tam - gamb*(zk-zs) - (gamt-gamb)*(zk**2-zs**2)/(2._wp*hgamt) 

    else
      ! temperature profile with surface layer

      if (zk.lt.zs) then
        ! virtual temperature profile below surface
        t_prof = tam - gamb*(zk-zs) - (gamt-gamb)*(zk**2-zs**2)/(2._wp*hgamt) 
      else if (zk.gt.(zs+hgams)) then
        ! temperature profile above boundary layer
        t_prof = tam - gams*hgams - gamb*(zk-(zs+hgams)) - (gamt-gamb)*(zk**2-zs**2)/(2._wp*hgamt) 
      else
        ! temperature profile in boundary layer
        t_prof = tam - gams*(zk-zs) - (gamt-gamb)*(zk**2-zs**2)/(2._wp*hgamt)
      endif

    endif

    return

  end function t_prof


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  r h _ p r o f
  !   Purpose    :  compute vertical relative humidity profile
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pure function rh_prof(zs, z, ram, h_rh, rh_fr, htrop)

    implicit none

    real(wp), intent(in) :: zs
    real(wp), intent(in) :: z
    real(wp), intent(in) :: ram
    real(wp), intent(in) :: h_rh
    real(wp), intent(in) :: rh_fr
    real(wp), intent(in) :: htrop

    real(wp) :: rh_prof

    real(wp) :: z_pbl

    ! Above the boundary layer the relative humidity relaxes exponentially from the
    ! boundary layer value ram towards the background free tropospheric value rh_fr, which
    ! is either the namelist scalar rh_free or the 2d field computed in hscales (i_rh_free).
    ! rh_fr = 0 gives back the purely multiplicative profile rh = ram*exp(-(z-z_pbl)/h_rh),
    ! for which rh(z)/ram depends on h_rh alone and is therefore identical over land and
    ! ocean. Observations show a free troposphere that is nearly the same over land and
    ! ocean while the boundary layer differs strongly, which no value of h_rh can reproduce
    ! without a background term.
    ! Height at which the boundary layer value gives way to the free troposphere.
    ! With i_zpbl=0 this follows the terrain, so over elevated land the relative humidity is
    ! held at ram up to c_hrs_5 above the local surface - 2.1 km above sea level at the mean
    ! NH midlatitude land elevation, against 1.0 km over the ocean. The free troposphere is
    ! quasi-horizontal, so once the profile above z_pbl is sharp this makes the elevated land
    ! boundary layer moister than the ocean through 1.3-1.8 km, where observations have the
    ! ocean moister. i_zpbl=1 references the transition to sea level instead, keeping only a
    ! minimum depth h_pbl_min above the local surface.
    if (i_zpbl.eq.0) then
      z_pbl = zs+c_hrs_5
    else
      z_pbl = max(zs+h_pbl_min, c_hrs_5)
    endif
    if (z.le.z_pbl) then
      rh_prof = ram
    else if (z.gt.z_pbl.and.z.le.(zs+c_hrs_4)) then
      rh_prof = rh_fr + (ram-rh_fr)*exp(-(z-z_pbl)/h_rh)
    else if (z.gt.(zs+c_hrs_4).and.z.le.(htrop+1.)) then
      rh_prof = rh_fr + (ram-rh_fr)*exp(-(zs+c_hrs_4-z_pbl)/h_rh)
    else
      rh_prof = rh_strat
    endif

    return

  end function rh_prof


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  t r o p o h e i g h t
  !   Purpose    :  compute height of tropopause
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine tropoheight(had_fi, had_width, rb_str, hcld, &
      htrop, ptrop)

    implicit none
    
    real(wp), intent(in   ) :: had_fi
    real(wp), intent(in   ) :: had_width
    real(wp), intent(in   ) :: rb_str(:,:)
    real(wp), intent(in   ) :: hcld(:,:)
    real(wp), intent(inout) :: htrop(:,:)
    real(wp), intent(out  ) :: ptrop(:)

    integer :: i, j
    real(wp) :: fi, fic, x, sheat, rbstr, dhtrop, htropp
    real(wp), parameter :: x1 = asin(0.1_wp**(1._wp/8._wp))  
    real(wp), parameter :: h_trop_min = 6.e3_wp
    real(wp), parameter :: h_trop_max = 25.e3_wp


    do j=1,jm
      fic = had_width/2._wp
      x = x1/fic
      fi = x*(fit(j)-had_fi)  
      if (fi.gt.pi/2._wp)  fi = pi/2._wp
      if (fi.lt.-pi/2._wp) fi = -pi/2._wp
      sheat = c_trop_2*(1._wp-c_trop_3*(1._wp-sin(fi)**8))
      ptrop(j) = 0._wp
      do i=1,im
        rbstr = rb_str(i,j) + sheat 
        dhtrop = -c_trop_1*rbstr
        htropp = htrop(i,j)+dhtrop
        htropp = max(htropp,h_trop_min) 
        htropp = max(htropp,hcld(i,j)+1000._wp)      
        htropp = min(htropp,h_trop_max)
        ! height of tropopause
        htrop(i,j) = htropp
        ! pressure at tropopause
        ptrop(j) = ptrop(j) + exp(-htrop(i,j)/hatm)*aim
      enddo
    enddo

    return

  end subroutine tropoheight

end module vesta_mod
