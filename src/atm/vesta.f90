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
  use constants, only : fqsat, pi, Rd, g
  use atm_params, only : cp, gad, hatm, p0, ra, zmax
  use atm_params, only : c_gam_1, c_gam_2, c_gam_3, hgams, hgamt, c_gam_tau, c_hrm_tau, nsmooth_gam, tstep
  use atm_params, only : gams_min, dgams, sh_gams_pos, sh_gams_0, sh_gams_neg
  use atm_params, only : c_hrs_1, c_hrs_2, c_hrs_3, c_hrs_4, c_hrs_5, c_hrs_6, rh_strat
  use atm_params, only : c_dhs_1, c_dhs_2
  use atm_params, only : c_trop_1, c_trop_2, c_trop_3
  use atm_params, only : l_dust
  use timer, only : sec_day
  use atm_grid, only : im, jm, km, kmc, aim, zl, zc, z3, z3l, fit, cost, exp_zc, i_lnd, i_ocn, dplt
  use smooth_atm_mod, only : smooth2
  !$ use omp_lib

  implicit none

  private
  public :: hscales, vesta, t_prof, rh_prof, tropoheight, tropopause_pressure

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  h s c a l e s
  !   Purpose    :  computation of lapse rate and height scales of moisture and dust
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine hscales(frst, ra2a, sha, tskina, tam, qam, ram, tsl, wcon, wcld, had_fi, had_width, &
      gams, gamb, gamt, hrm, &
      hqeff, hdust)

    implicit none

    real(wp), intent(in) :: frst(:,:,:)
    real(wp), intent(in) :: ra2a(:,:)
    real(wp), intent(in) :: sha(:,:)   !! grid mean sensible heat flux, positive upward (W/m2)
    real(wp), intent(in) :: tskina(:,:)
    real(wp), intent(in) :: tam(:,:)
    real(wp), intent(in) :: qam(:,:)
    real(wp), intent(in) :: ram(:,:)   !! surface relative humidity
    real(wp), intent(in) :: tsl(:,:)   !! sea level temperature (K), extrapolated with the previous gamb/gamt
    real(wp), intent(in) :: wcon(:,:)
    real(wp), intent(in) :: wcld(:,:)
    real(wp), intent(in) :: had_fi
    real(wp), intent(in) :: had_width

    real(wp), intent(inout) :: gams(:,:)
    real(wp), intent(inout) :: gamb(:,:)
    real(wp), intent(inout) :: gamt(:,:)
    real(wp), intent(inout) :: hrm(:,:)

    real(wp), intent(out) :: hqeff(:,:)
    real(wp), intent(out) :: hdust(:,:)

    integer :: i, j
    real(wp) :: hrs, fi, f_trop, cgam, chrm
    real(wp) :: sh_w, q_gam, gmax
    real(wp), dimension(im,jm) :: gam_s, gam_b, gam_t

    real(wp), parameter :: hrs_min = 1.e3_wp
    real(wp), parameter :: hrs_max = 10.e3_wp


    ! time relaxation weights per fast step, hscales is called every fast step
    cgam = min(1._wp, tstep/(c_gam_tau*sec_day))
    chrm = min(1._wp, tstep/(c_hrm_tau*sec_day))

    !$omp parallel do collapse(2) private(i,j,hrs,fi,f_trop,sh_w,q_gam,gmax)
    do j=1,jm
      do i=1,im

        !----------------------------------------------
        ! lapse rate

        ! humidity controlling the free-tropospheric lapse rate.
        ! gamb is the lapse rate at sea level (t_prof is linear in absolute z), so it uses the humidity
        ! at sea level, ram*qsat(tsl,p0), not the surface value qam, which over high terrain is the much
        ! drier air at zsa. Over ocean and lowlands the two are nearly identical.
        ! tsl is the previous step's value, extrapolated below zsa with gamb itself. The loop is a
        ! negative feedback (larger gamb -> warmer tsl -> moister -> smaller gamb) and is damped by the
        ! c_gam_tau relaxation of gamb.
        q_gam = ram(i,j)*fqsat(tsl(i,j),p0)

        ! bottom
        gam_b(i,j) = c_gam_1 - c_gam_2*q_gam

        ! top
        gam_t(i,j) = gam_b(i,j) - c_gam_2*q_gam + c_gam_3

        ! lapse rate in the surface layer (zs..zs+hgams), on top of which t_prof adds the same curvature
        ! term as in the free troposphere.
        ! One tanh in the grid-mean sensible heat flux, between a stable and an unstable limit:
        ! - sha > sh_gams_0 (unstable): gam_s -> gam_b + dgams over ice-free ocean and land (gam_b elsewhere).
        !   Over the ocean the subcloud/cloud layer is ~2 K/km steeper than the free troposphere (ERA 2.2,
        !   CMIP median 2.0). Over land the daily-mean excess is smaller (ERA +1.2, CMIP median +1.0 K/km,
        !   SH > 5 W/m2), because the daytime superadiabatic layer cancels against the nocturnal inversion;
        !   the same dgams is used anyway, the steeper land surface layer partly compensating a too cold
        !   land column (2 m temperature). The tanh is not fully saturated at typical SH, so the realised
        !   excess is somewhat below dgams.
        ! - sha < sh_gams_0 (stable, surface-based inversion): gam_s -> gams_min.
        ! Width sh_gams_pos on the unstable side and sh_gams_neg on the stable side. A narrow stable
        ! width saturates calm, weakly cooled (decoupled) surfaces at gams_min, while the floor
        ! gams_min caps the inversion of strongly cooled, wind-mixed surfaces (katabatic ice sheets).
        if (sha(i,j).ge.sh_gams_0) then
          sh_w = sh_gams_pos
        else
          sh_w = sh_gams_neg
        endif
        gmax = gam_b(i,j) + dgams*(frst(i,j,i_ocn)+frst(i,j,i_lnd))
        gam_s(i,j) = gams_min + (gmax-gams_min)*0.5_wp*(1._wp+tanh((sha(i,j)-sh_gams_0)/sh_w))

        !----------------------------------------------
        ! height scale for relative humidity       

        fi = c_hrs_6*(fit(j)-had_fi)/(0.5_wp*had_width)
        fi = min(fi,pi/2._wp)
        fi = max(fi,-pi/2._wp)
        f_trop = 1._wp-sin(fi)**8
        hrs = f_trop * c_hrs_1*exp(c_hrs_2*wcld(i,j)) + (1._wp-f_trop) * c_hrs_1*c_hrs_3 
        hrs = max(hrs,hrs_min)   
        hrs = min(hrs,hrs_max)   
        hrm(i,j) = (1._wp-chrm)*hrm(i,j) + chrm*hrs 

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
        gams(i,j) = (1._wp-cgam)*gams(i,j) + cgam*gam_s(i,j)
      enddo
    enddo
    !$omp section
    call smooth2(gam_b,nsmooth_gam)
    do j=1,jm
      do i=1,im
        gamb(i,j) = (1._wp-cgam)*gamb(i,j) + cgam*gam_b(i,j)
      enddo
    enddo
    !$omp section
    call smooth2(gam_t,nsmooth_gam)
    do j=1,jm
      do i=1,im
        gamt(i,j) = (1._wp-cgam)*gamt(i,j) + cgam*gam_t(i,j)
      enddo
    enddo
    !$omp end parallel sections

    return

  end subroutine hscales


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  v e s t a
  !   Purpose    :  vertical structure of atmosphere 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! t_prof and rh_prof are formulated in geometric height (lapse rates in K/m, zsa, hgams,
  ! htrop, hrm), while the dynamics treat t3 as the temperature on the pressure levels pl*p0
  ! (thermal wind, layer masses). Each level is therefore placed at the geometric height that
  ! the column's own temperature profile gives it, from the hypsometric equation integrated up-
  ! and downward from the surface (hypso_column); z3/z3l (atm_grid) keep these heights for
  ! synop, zslp and the output. The level pressure p0*exp_zc in qsat and the layer mass
  ! dvol (= dp/g) do not depend on it, so A_trop and W_strat are mass integrals.
  subroutine vesta(zsa, tam, gams, gamb, gamt, htrop, ram, hrm, dam, hdust, &
      A_trop, W_strat, dtcol, t3, q3, tp, d3, ttrop)

    implicit none

    real(wp), intent(in ) :: zsa(:,:)
    real(wp), intent(in ) :: tam(:,:)
    real(wp), intent(in ) :: gams(:,:)
    real(wp), intent(in ) :: gamb(:,:)
    real(wp), intent(in ) :: gamt(:,:)
    real(wp), intent(in ) :: htrop(:,:)
    real(wp), intent(in ) :: ram(:,:)
    real(wp), intent(in ) :: hrm(:,:)
    real(wp), intent(in ) :: dam(:,:)
    real(wp), intent(in ) :: hdust(:,:)

    real(wp), intent(out) :: A_trop(:,:)   
    real(wp), intent(out) :: W_strat(:,:)  
    real(wp), intent(out) :: dtcol(:,:)    
    real(wp), intent(out) :: t3(:,:,:)
    real(wp), intent(out) :: q3(:,:,:)
    real(wp), intent(out) :: tp(:,:,:)
    real(wp), intent(out) :: d3(:,:,:)
    real(wp), intent(out) :: ttrop(:,:)

    integer :: i, j, k
    logical :: flag_strat
    real(wp) :: z_sur, taml, htropl
    real(wp) :: gamsl, gambl, gamtl, z
    real(wp) :: t, rsur, hrml, rh, rh_unit, q, q_unit, qsat, A_l, W_l
    real(wp) :: dvol, S_num, S_den


    !$omp parallel do collapse(2) private(i,j,k,z_sur,taml,htropl,gamsl,gambl,gamtl,z,t,rsur,hrml,rh,rh_unit,q,q_unit,qsat,A_l,W_l,dvol,S_num,S_den,flag_strat)
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

        ! geometric heights of the layer centres and interfaces of this column
        call hypso_column(z_sur, taml, gamsl, gambl, gamtl, htropl, z3(i,j,:), z3l(i,j,:))

        ! 3D fields of temperature and humidity

        A_l   = 0._wp   ! tropospheric coefficient ∫ g(z)·qsat·ρ dz
        W_l   = 0._wp   ! stratospheric intercept ∫ rh_strat·qsat·ρ dz
        S_num = 0._wp   ! ∫ (t_prof - tam) dp, the shape offset of the column
        S_den = 0._wp   ! ∫ dp
        flag_strat = .false.
        do k=1,km

          z = z3(i,j,k)

          if (.not.flag_strat) then
            ! construct vertical temperature profile
            t = t_prof(z_sur, z, taml, gamsl, gambl, gamtl, htropl, 1)
            ! derive specific humidity profile from temperature and relative humidity profiles
            rh = rh_prof(z_sur, z, rsur, hrml, htropl)
            ! rh_prof is proportional to ram in the troposphere: rh = ram*E(z).
            ! Evaluating it at ram = 1 recovers E(z) without duplicating the profile.
            rh_unit = rh_prof(z_sur, z, 1._wp, hrml, htropl)
            qsat = fqsat(t,p0*exp_zc(k))
            q = rh*qsat
            q_unit = rh_unit*qsat
          endif

          ! volume element (mass per unit area for this layer above surface)
          if (zl(k).ge.z_sur) then
            dvol = ra*exp_zc(k)*(zl(k+1)-zl(k))
          else if (zl(k).lt.z_sur .and. zl(k+1).gt.z_sur) then
            dvol = ra*exp_zc(k)*(zl(k+1)-z_sur)
          else
            dvol = 0._wp
          endif

          ! vertical integral of water content split into tropospheric slope (A) and stratospheric intercept (W);
          ! the total column water wcon = ram·A_trop + W_strat is reconstructed in time_step.
          if (z.le.htropl+1._wp) then
            ! tropospheric: q = ram·E(z)·qsat ⇒ dA = q_unit·dvol
            A_l = A_l + q_unit*dvol
          else
            ! stratospheric: q = rh_strat·qsat ⇒ dW = q·dvol (independent of ram)
            W_l = W_l + q*dvol
          endif

          t3(i,j,k) = t
          q3(i,j,k) = q

          ! Offset between the mass-weighted column temperature and the profile anchor tam.
          ! The temperature budget advances tam (time_step.f90) while t_prof(zs) = tam exactly,
          ! so the column heat content is cheat*(tam + dtcol) and NOT cheat*tam, cheat = pzsa*amas*cp
          ! being the column enthalpy per K (atm_grid) and sum(dplt) = pzsa*amas the same mass.  Whenever
          ! the profile shape moves - gams/gamb/gamt with qam and sha, htrop with the state -
          ! the column's heat content changes with no term in deba.  
          S_num = S_num + (t-taml)*dplt(i,j,k)
          S_den = S_den + dplt(i,j,k)

          ! potential temperature
          tp(i,j,k) = t + gad*min(z,zmax)

          ! dust profile
          if (l_dust) d3(i,j,k) = dam(i,j)*min(1._wp,exp(-(z-z_sur)/hdust(i,j))) ! kg/kg

          if (z.gt.htropl) flag_strat = .true.

        enddo

        A_trop(i,j)  = A_l
        W_strat(i,j) = W_l
        dtcol(i,j)   = S_num/max(S_den,1.e-20_wp)

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
  pure function rh_prof(zs, z, ram, h_rh, htrop)

    implicit none

    real(wp), intent(in) :: zs
    real(wp), intent(in) :: z
    real(wp), intent(in) :: ram
    real(wp), intent(in) :: h_rh   !! hrm
    real(wp), intent(in) :: htrop

    real(wp) :: rh_prof

    real(wp) :: z_pbl

    z_pbl = max(zs+500._wp, c_hrs_5)
    if (z.le.z_pbl) then
      rh_prof = ram
    else if (z.gt.z_pbl.and.z.le.(zs+c_hrs_4)) then
      rh_prof = ram*exp(-(z-z_pbl)/h_rh)
    else if (z.gt.(zs+c_hrs_4).and.z.le.(htrop+1.)) then
      rh_prof = ram*exp(-(zs+c_hrs_4-z_pbl)/h_rh)
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
      htrop, rb_str_ref)

    implicit none

    real(wp), intent(in   ) :: had_fi
    real(wp), intent(in   ) :: had_width
    real(wp), intent(in   ) :: rb_str(:,:)
    real(wp), intent(in   ) :: hcld(:,:)
    real(wp), intent(inout) :: htrop(:,:)
    ! optional reference stratospheric radiative balance, used for the feedback analysis.
    ! If present, the tropopause is relaxed towards the state in which the stratospheric
    ! radiative balance matches rb_str_ref instead of the balance with the prescribed
    ! stratospheric heating sheat. The dynamical heating is then held fixed and cancels.
    real(wp), intent(in   ), optional :: rb_str_ref(:,:)

    integer :: i, j
    real(wp) :: fi, fic, x, sheat, rbstr, dhtrop, htropp
    logical :: fixed_heat
    real(wp), parameter :: x1 = asin(0.1_wp**(1._wp/8._wp))  
    real(wp), parameter :: h_trop_min = 6.e3_wp
    real(wp), parameter :: h_trop_max = 25.e3_wp


    fixed_heat = present(rb_str_ref)

    do j=1,jm
      fic = had_width/2._wp
      x = x1/fic
      fi = x*(fit(j)-had_fi)  
      if (fi.gt.pi/2._wp)  fi = pi/2._wp
      if (fi.lt.-pi/2._wp) fi = -pi/2._wp
      sheat = c_trop_2*(1._wp-c_trop_3*(1._wp-sin(fi)**8))
      do i=1,im
        if (fixed_heat) then
          rbstr = rb_str(i,j) - rb_str_ref(i,j)
        else
          rbstr = rb_str(i,j) + sheat
        endif
        dhtrop = -c_trop_1*rbstr
        htropp = htrop(i,j)+dhtrop
        htropp = max(htropp,h_trop_min) 
        htropp = max(htropp,hcld(i,j)+1000._wp)      
        htropp = min(htropp,h_trop_max)
        ! height of tropopause
        htrop(i,j) = htropp
      enddo
    enddo

    return

  end subroutine tropoheight


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  h y p s o _ c o l u m n
  !   Purpose    :  geometric height of the pressure levels of one column
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pure subroutine hypso_column(zs, tam, gams, gamb, gamt, htrop, zc_geo, zl_geo)

    ! Solves dz/dlnp = -Rd*T(z)/g with T = t_prof(z) exactly, from the surface, at pressure
    ! p0*exp(-zs/hatm) (the surface of the layer masses), upward to the levels above it and
    ! downward, along the sub-surface extrapolation of t_prof, to those below it. Targets are
    ! the interfaces zl and the layer centres zc, i.e. ln(p/p0) = -zl/hatm and -zc/hatm.
    ! Within each segment of t_prof (0: below zs, 1: surface layer zs..zs+hgams, 2: free
    ! troposphere up to htrop) T is quadratic in z, T(z0+u) = T0 - (gam+2*qq*z0)*u - qq*u**2 with
    ! qq = (gamt-gamb)/(2*hgamt), so in l = ln(p/p0) the equation is a constant-coefficient Riccati
    ! equation du/dl = a + b*u + c*u**2, a = -kz*T0, b = kz*(gam+2*qq*z0), c = kz*qq, kz = Rd/g,
    ! with the solution u(0) = 0
    !   u = a*S/(C - b/2*S),  C = cosh(r*dl), S = sinh(r*dl)/r,  r**2 = b**2/4 - a*c
    ! (cos/sin for r**2 < 0; both are the same power series in x = r**2*dl**2, summed here to x**9,
    ! exact to roundoff for |x| < 2, while |x| < 0.4 for any tropospheric segment). The pressure
    ! at the top of a segment follows from the inverse, dl = rho*atanh(r*rho)/(r*rho) with
    ! rho = u/(a + b*u/2) (atan for r**2 < 0). Above htrop T = t_prof(htrop) and
    ! z = htrop - (l - l_trop)*kz*ttrop. Every level is solved from the base of its segment, so the
    ! level loop has no dependency and vectorises. Agrees with an RK4 reference (dlnp 1e-5) to
    ! 5e-8 m, against ~1 m for the former RK2 with dlnp_max = 0.05, at 1/4.6 of the cost
    ! (analysis/hypso_accuracy/hypso_exact_test.f90).

    implicit none

    real(wp), intent(in ) :: zs, tam, gams, gamb, gamt, htrop
    real(wp), intent(out) :: zc_geo(:)
    real(wp), intent(out) :: zl_geo(:)

    integer :: k, n, nt
    real(wp) :: lnp_t(2*km+1), z_t(2*km+1), lnps, ttrop, qq, kz, zbl, ztr, lnp_bl, lnp_tr
    real(wp) :: a0, b0, a1, b1, a2, b2, c, t_bl
    real(wp) :: a, b, z0, l0, dl, x, cc, sl

    kz = Rd/g
    qq = (gamt-gamb)/(2._wp*hgamt)
    c  = kz*qq
    ztr = max(htrop,zs)
    zbl = min(zs+hgams, ztr)
    lnps = -zs/hatm
    ttrop = t_prof(zs, htrop, tam, gams, gamb, gamt, htrop, 1)

    ! segments 0 (below ground) and 1 (surface layer), from the surface
    a0 = -kz*tam
    b0 = kz*(gamb + 2._wp*qq*zs)
    a1 = a0
    b1 = kz*(gams + 2._wp*qq*zs)
    lnp_bl = lnps + ric_dl(a1, b1, c, zbl-zs)
    ! segment 2 (free troposphere), from the top of the surface layer
    t_bl = tam - gams*(zbl-zs) - qq*(zbl**2-zs**2)
    a2 = -kz*t_bl
    b2 = kz*(gamb + 2._wp*qq*zbl)
    lnp_tr = lnp_bl + ric_dl(a2, b2, c, ztr-zbl)

    ! targets from the bottom up: interface 1, centre 1, interface 2, ..., interface kmc
    nt = 2*km+1
    do k=1,km
      lnp_t(2*k-1) = -zl(k)/hatm
      lnp_t(2*k)   = -zc(k)/hatm
    enddo
    lnp_t(nt) = -zl(kmc)/hatm

    !$omp simd private(a,b,z0,l0,dl,x,cc,sl)
    do n=1,nt
      if (lnp_t(n).ge.lnps) then
        a = a0 ; b = b0 ; z0 = zs  ; l0 = lnps
      else if (lnp_t(n).ge.lnp_bl) then
        a = a1 ; b = b1 ; z0 = zs  ; l0 = lnps
      else
        a = a2 ; b = b2 ; z0 = zbl ; l0 = lnp_bl
      endif
      dl = lnp_t(n)-l0
      x = (0.25_wp*b*b - a*c)*dl*dl
      ! C = sum x^m/(2m)!, S/dl = sum x^m/(2m+1)!
      cc = 1._wp+x*(1._wp/2._wp+x*(1._wp/24._wp+x*(1._wp/720._wp+x*(1._wp/40320._wp+x*(1._wp/3628800._wp &
         + x*(1._wp/479001600._wp+x*(1._wp/87178291200._wp+x*(1._wp/20922789888000._wp &
         + x*(1._wp/6402373705728000._wp)))))))))
      sl = 1._wp+x*(1._wp/6._wp+x*(1._wp/120._wp+x*(1._wp/5040._wp+x*(1._wp/362880._wp+x*(1._wp/39916800._wp &
         + x*(1._wp/6227020800._wp+x*(1._wp/1307674368000._wp+x*(1._wp/355687428096000._wp &
         + x*(1._wp/121645100408832000._wp)))))))))
      z_t(n) = z0 + a*sl*dl/(cc - 0.5_wp*b*sl*dl)
      ! isothermal stratosphere
      if (lnp_t(n).lt.lnp_tr) z_t(n) = ztr - (lnp_t(n)-lnp_tr)*kz*ttrop
    enddo

    do k=1,km
      zl_geo(k) = z_t(2*k-1)
      zc_geo(k) = z_t(2*k)
    enddo
    zl_geo(kmc) = z_t(nt)

  contains

    ! dl at which du/dl = a + b*u + c*u**2, u(0) = 0, reaches u
    pure function ric_dl(a, b, c, u) result(dl)
      real(wp), intent(in) :: a, b, c, u
      real(wp) :: dl, rho, y, f, rr
      rho = u/(a + 0.5_wp*b*u)
      y = (0.25_wp*b*b - a*c)*rho*rho
      if (abs(y).lt.1.e-8_wp) then
        f = 1._wp + y/3._wp
      else if (y.gt.0._wp) then
        rr = sqrt(y)
        f = atanh(rr)/rr
      else
        rr = sqrt(-y)
        f = atan(rr)/rr
      endif
      dl = rho*f
    end function ric_dl

  end subroutine hypso_column


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  t r o p o p a u s e _ p r e s s u r e
  !   Purpose    :  zonal mean relative tropopause pressure from the hypsometric equation
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine tropopause_pressure(zsa, tam, gams, gamb, gamt, htrop, ptrop)

    ! p(htrop)/p0 = exp(-zsa/hatm)*exp(-int_zsa^htrop g/(Rd*T) dz), T = t_prof, consistent with
    ! the geometric heights of the pressure levels in vesta.

    implicit none

    real(wp), intent(in ) :: zsa(:,:)
    real(wp), intent(in ) :: tam(:,:)
    real(wp), intent(in ) :: gams(:,:)
    real(wp), intent(in ) :: gamb(:,:)
    real(wp), intent(in ) :: gamt(:,:)
    real(wp), intent(in ) :: htrop(:,:)
    real(wp), intent(out) :: ptrop(:)

    real(wp), parameter :: dz_max = 250._wp
    integer :: i, j, n, ns
    real(wp) :: h, s

    !$omp parallel do private(i,j,n,ns,h,s)
    do j=1,jm
      ptrop(j) = 0._wp
      do i=1,im
        ns = max(1,ceiling(max(0._wp,htrop(i,j)-zsa(i,j))/dz_max))
        h = max(0._wp,htrop(i,j)-zsa(i,j))/real(ns,wp)
        s = 0._wp
        do n=1,ns
          s = s + g*h/(Rd*t_prof(zsa(i,j), zsa(i,j)+(n-0.5_wp)*h, tam(i,j), gams(i,j), gamb(i,j), gamt(i,j), htrop(i,j), 1))
        enddo
        ptrop(j) = ptrop(j) + exp(-zsa(i,j)/hatm-s)*aim
      enddo
    enddo

    return

  end subroutine tropopause_pressure

end module vesta_mod
