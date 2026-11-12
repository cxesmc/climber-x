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
  use constants, only : fqsat, pi, Le, Rd, Rv, g
  use atm_params, only : cp, gad, hatm, p0, ra, zmax
  use atm_params, only : i_gamb_q, q_gam_ref, c_gam_w, c_gam_twid
  use atm_params, only : c_gam_1, c_gam_2, c_gam_3, hgams, hgamt, c_gam_tau, c_hrm_tau, nsmooth_gam, tstep
  use atm_params, only : gams_min, gams_max, sh_gams, gams_max_conv, sh_gams_conv, dsh_gams_conv
  use atm_params, only : c_hrs_1, c_hrs_2, c_hrs_3, c_hrs_4, c_hrs_5, c_hrs_6, rh_strat
  use atm_params, only : c_dhs_1, c_dhs_2
  use atm_params, only : c_trop_1, c_trop_2, c_trop_3
  use atm_params, only : l_dust
  use timer, only : sec_day
  use atm_grid, only : im, jm, km, kmc, aim, zl, zc, z3, z3l, fit, cost, exp_zc, i_lnd, dplt
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
  subroutine hscales(frst, ra2a, sha, tskina, tam, qam, wcon, wcld, had_fi, had_width, &
      gams, gamb, gamt, hrm, &
      hqeff, hdust)

    implicit none

    real(wp), intent(in) :: frst(:,:,:)
    real(wp), intent(in) :: ra2a(:,:)
    real(wp), intent(in) :: sha(:,:)   !! grid mean sensible heat flux, positive upward (W/m2)
    real(wp), intent(in) :: tskina(:,:)
    real(wp), intent(in) :: tam(:,:)
    real(wp), intent(in) :: qam(:,:)
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
    real(wp) :: qtrp, wtrp, dq_gam, x_gam
    real(wp) :: t1had, p1, qs1, gm_had, w_had, dlat, f_gam
    real(wp), dimension(im,jm) :: gam_s, gam_b, gam_t

    real(wp), parameter :: hrs_min = 1.e3_wp
    real(wp), parameter :: hrs_max = 10.e3_wp
    real(wp), parameter :: lat_gam_trop = 30._wp*pi/180._wp
    real(wp), parameter :: x1_gam = asin(0.1_wp**(1._wp/8._wp))


    ! time relaxation weights per fast step, hscales is called every fast step
    cgam = min(1._wp, tstep/(c_gam_tau*sec_day))
    chrm = min(1._wp, tstep/(c_hrm_tau*sec_day))

    ! tropical mean (30S-30N) qam, which drives the free tropospheric lapse rate for i_gamb_q = 1
    qtrp = 0._wp
    wtrp = 0._wp
    if (i_gamb_q.eq.1) then
      do j=1,jm
        if (abs(fit(j)).lt.lat_gam_trop) then
          do i=1,im
            qtrp = qtrp + qam(i,j)*cost(j)
            wtrp = wtrp + cost(j)
          enddo
        endif
      enddo
      qtrp = qtrp/wtrp
    endif

    ! i_gamb_q = 2: moist adiabatic lapse rate of the mean state over the Hadley cells.
    ! T1 is the mean temperature at the top of the surface layer (tam - gams*hgams), averaged
    ! between the two Hadley cell edges (had_fi -/+ had_width/2, edge rows by their fractional
    ! overlap), so the averaging region follows the circulation through the season and across
    ! climates. Gamma_m is evaluated at T1 and at p(z*=hgams), in K per geometric metre like t_prof.
    ! No reference state enters.
    gm_had = 0._wp
    if (i_gamb_q.eq.2) then
      dlat = abs(fit(2)-fit(1))
      t1had = 0._wp
      wtrp = 0._wp
      do j=1,jm
        w_had = cost(j)*max(0._wp,min(1._wp,(0.5_wp*had_width-abs(fit(j)-had_fi))/dlat+0.5_wp))
        if (w_had.gt.0._wp) then
          do i=1,im
            t1had = t1had + (tam(i,j)-gams(i,j)*hgams)*w_had
            wtrp = wtrp + w_had
          enddo
        endif
      enddo
      t1had = t1had/wtrp
      p1 = p0*exp(-hgams/hatm)
      qs1 = fqsat(t1had,p1)
      gm_had = g/cp*(1._wp+Le*qs1/(Rd*t1had)) / (1._wp+Le**2*qs1*(Rd/Rv)/(cp*Rd*t1had**2))
    endif

    !$omp parallel do collapse(2) private(i,j,hrs,fi,f_trop,dq_gam,x_gam,f_gam)
    do j=1,jm
      do i=1,im

        !----------------------------------------------
        ! lapse rate

        ! lapse rate in the boundary layer.
        ! One tanh in the grid-mean sensible heat flux replaces the per-surface-type branching
        ! (9 constants and a switch over ocean / land / lake / ice). sha > 0 is an unstable
        ! surface layer and drives gam_s towards the well mixed limit gams_max; sha < 0 is a
        ! stable, decoupled surface layer and drives it towards gams_min. sh_gams sets the flux
        ! scale over which the surface layer changes regime.
        ! A second tanh, centred at the much larger flux sh_gams_conv, raises the ceiling from
        ! gams_max to gams_max_conv only in very unstable (strongly convective) surface layers,
        ! leaving the curve unchanged at moderate fluxes.
        gam_s(i,j) = gams_min + (gams_max-gams_min)*0.5_wp*(1._wp+tanh(sha(i,j)/sh_gams)) &
                   + (gams_max_conv-gams_max)*0.5_wp*(1._wp+tanh((sha(i,j)-sh_gams_conv)/dsh_gams_conv))

        if (i_gamb_q.eq.1) then
          ! The tropical free troposphere follows the moist adiabat of the convecting regions and
          ! weak temperature gradients spread it, so the lapse rate responds to the tropical mean
          ! moisture anomaly, not to the local humidity. The anomaly vanishes in the reference
          ! climate, so c_gam_1 and c_gam_3 set the mean state and c_gam_2 the feedback alone.
          ! Same latitude taper as the stratospheric heating in tropoheight.
          x_gam = x1_gam/(0.5_wp*had_width)*(fit(j)-had_fi)
          x_gam = max(-pi/2._wp,min(pi/2._wp,x_gam))
          dq_gam = (1._wp-sin(x_gam)**8)*(qtrp-q_gam_ref)
        else
          dq_gam = qam(i,j)
        endif

        if (i_gamb_q.eq.2) then
          ! bottom: relax c_gam_1 towards the moist adiabat of the Hadley cells, by a fraction
          ! c_gam_w in the tropics, tapered to 0 outside them (taper as tropoheight, half-width
          ! scaled by c_gam_twid). top: same relation as the moisture closure below, with
          ! c_gam_2*q replaced by c_gam_1-gam_b.
          x_gam = x1_gam/(c_gam_twid*0.5_wp*had_width)*(fit(j)-had_fi)
          x_gam = max(-pi/2._wp,min(pi/2._wp,x_gam))
          f_gam = 1._wp-sin(x_gam)**8
          gam_b(i,j) = c_gam_1 + c_gam_w*f_gam*(gm_had-c_gam_1)
          gam_t(i,j) = 2._wp*gam_b(i,j) - c_gam_1 + c_gam_3
        else
          ! bottom
          gam_b(i,j) = c_gam_1 - c_gam_2*dq_gam

          ! top
          gam_t(i,j) = gam_b(i,j) - c_gam_2*dq_gam + c_gam_3
        endif

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

    ! Integrates dz/dlnp = -Rd*T(z)/g with T = t_prof(z) from the surface, at pressure
    ! p0*exp(-zs/hatm) (the surface of the layer masses), upward to the levels above it and
    ! downward, along the sub-surface extrapolation of t_prof, to those below it. Targets are
    ! the interfaces zl and the layer centres zc, i.e. ln(p/p0) = -zl/hatm and -zc/hatm.
    ! Midpoint (RK2) steps of at most dlnp_max. Above the tropopause T = t_prof(htrop) is
    ! constant and the upward integration there is exact, z1 = z0 - (lnp1-lnp0)*Rd*ttrop/g.
    ! dlnp_max = 0.05 keeps the error within ~1 m of an RK4 reference (offline test,
    ! analysis/hypso_accuracy); 0.1/0.2 give 5/9 m for only ~25% fewer t_prof calls.

    implicit none

    real(wp), intent(in ) :: zs, tam, gams, gamb, gamt, htrop
    real(wp), intent(out) :: zc_geo(:)
    real(wp), intent(out) :: zl_geo(:)

    real(wp), parameter :: dlnp_max = 0.05_wp
    integer :: k, n, n0, nt
    real(wp) :: lnp_t(2*km+1), z_t(2*km+1), z, lnp, lnps, ttrop

    ttrop = t_prof(zs, htrop, tam, gams, gamb, gamt, htrop, 1)

    ! targets from the bottom up: interface 1, centre 1, interface 2, ..., interface kmc
    nt = 2*km+1
    do k=1,km
      lnp_t(2*k-1) = -zl(k)/hatm
      lnp_t(2*k)   = -zc(k)/hatm
    enddo
    lnp_t(nt) = -zl(kmc)/hatm
    lnps = -zs/hatm

    ! first target at or above the surface
    n0 = nt+1
    do n=1,nt
      if (lnp_t(n).le.lnps) then
        n0 = n
        exit
      endif
    enddo

    ! upward from the surface
    z = zs
    lnp = lnps
    do n=n0,nt
      z = hypso_dz(z, lnp, lnp_t(n))
      lnp = lnp_t(n)
      z_t(n) = z
    enddo

    ! downward from the surface, below ground
    z = zs
    lnp = lnps
    do n=n0-1,1,-1
      z = hypso_dz(z, lnp, lnp_t(n))
      lnp = lnp_t(n)
      z_t(n) = z
    enddo

    do k=1,km
      zl_geo(k) = z_t(2*k-1)
      zc_geo(k) = z_t(2*k)
    enddo
    zl_geo(kmc) = z_t(nt)

  contains

    pure function hypso_dz(z0, lnp0, lnp1) result(z1)
      real(wp), intent(in) :: z0, lnp0, lnp1
      real(wp) :: z1, h, zm
      integer :: ns, s
      ns = max(1,ceiling(abs(lnp1-lnp0)/dlnp_max))
      h = (lnp1-lnp0)/real(ns,wp)
      z1 = z0
      do s=1,ns
        if (h.lt.0._wp .and. z1.ge.htrop) then
          ! upward in the isothermal stratosphere: exact for the rest of the interval
          z1 = z1 - (lnp1-(lnp0+real(s-1,wp)*h))*Rd*ttrop/g
          return
        endif
        zm = z1 - 0.5_wp*h*Rd*t_prof(zs, z1, tam, gams, gamb, gamt, htrop, 1)/g
        z1 = z1 - h*Rd*t_prof(zs, zm, tam, gams, gamb, gamt, htrop, 1)/g
      enddo
    end function hypso_dz

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
