!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : t i m e _ s t e p _ m o d
!
!  Purpose : time integration of equations for temperature, humidity and dust
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
module time_step_mod

  use atm_params, only : wp
  use precision, only : dp
  use constants, only : T0, fqsat, q_sat_w, q_sat_i
  use timer, only : sec_day, year, doy
  use atm_params, only : tstep, amas, hatm, ra, cv, cle, cls, l_dust, rh_max, rskin_ocn_min, gams_max_ocn, tsl_gams_min_lnd, tsl_gams_min_ice, i_tsl, i_tslz, c_tsl_gam, c_tsl_gam_ice, hgams
  use atm_params, only : c_wrt_1, c_wrt_2, c_wrt_3, c_wrt_4
  use control, only : check_water, check_energy
  use atm_grid, only : im, jm, nm, i_ocn, i_sic, i_lake, i_ice, i_lnd, sqr
  use vesta_mod, only : t_prof
  !$ use omp_lib

  implicit none

  private
  public :: time_step
  
contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  t i m e _ s t e p
  !   Purpose    :  time integration of equations for temperature, humidity and dust
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine time_step(frst, zs, zsa, ps, psa, ra2a, slope, evpa, convwtr, wcon, A_trop, W_strat, sam, eke, sam2, &
      tskin, convdse, rb_atm, rb_sur, sha, gams, gamb, gamt, &
      convdst, dust_emis, dust_dep, hdust, &
      convco2, co2flx, &
      tam, qam, dam, cam, prc, prcw, prcs, prc_conv, prc_wcon, prc_over, &
      q2, q2a, ram, r2, r2a, rskina, tsl, tsksl, t2, t2a, tskina, error)

    implicit none

    real(wp), intent(in   ) :: frst(:,:,:)
    real(wp), intent(in   ) :: zs(:,:,:)
    real(wp), intent(in   ) :: zsa(:,:)
    real(wp), intent(in   ) :: ps(:,:,:)
    real(wp), intent(in   ) :: psa(:,:)
    real(wp), intent(in   ) :: ra2a(:,:)
    real(wp), intent(in   ) :: slope(:,:)
    real(wp), intent(in   ) :: evpa(:,:)
    real(wp), intent(in   ) :: convwtr(:,:)
    real(wp), intent(inout) :: wcon(:,:)
    real(wp), intent(in   ) :: A_trop(:,:)
    real(wp), intent(in   ) :: W_strat(:,:)
    real(wp), intent(in   ) :: sam(:,:)
    real(wp), intent(in   ) :: eke(:,:)
    real(wp), intent(in   ) :: sam2(:,:)
    real(wp), intent(in   ) :: tskin(:,:,:)
    real(wp), intent(in   ) :: convdse(:,:)
    real(wp), intent(in   ) :: rb_atm(:,:)
    real(wp), intent(in   ) :: rb_sur(:,:)
    real(wp), intent(in   ) :: sha(:,:)
    real(wp), intent(in   ) :: gams(:,:)
    real(wp), intent(in   ) :: gamb(:,:)
    real(wp), intent(in   ) :: gamt(:,:)
    real(wp), intent(in   ) :: convdst(:,:)
    real(wp), intent(in   ) :: dust_emis(:,:)
    real(wp), intent(in   ) :: dust_dep(:,:)
    real(wp), intent(in   ) :: hdust(:,:)
    real(wp), intent(in   ) :: convco2(:,:)
    real(wp), intent(in   ) :: co2flx(:,:)

    real(wp), intent(inout) :: tam(:,:)
    real(wp), intent(inout) :: qam(:,:)
    real(wp), intent(inout) :: dam(:,:)
    real(wp), intent(inout) :: cam(:,:)
    real(wp), intent(inout) :: prc(:,:)
    real(wp), intent(inout) :: prcw(:,:,:)
    real(wp), intent(inout) :: prcs(:,:,:)
    real(wp), intent(inout) :: prc_conv(:,:)
    real(wp), intent(inout) :: prc_wcon(:,:)
    real(wp), intent(inout) :: prc_over(:,:)

    real(wp), intent(out  ) :: q2(:,:,:)
    real(wp), intent(out  ) :: q2a(:,:)
    real(wp), intent(inout) :: ram(:,:)
    real(wp), intent(out  ) :: r2(:,:,:)
    real(wp), intent(out  ) :: r2a(:,:)
    real(wp), intent(out  ) :: rskina(:,:)
    real(wp), intent(out  ) :: tsl(:,:)
    real(wp), intent(out  ) :: tsksl(:,:)
    real(wp), intent(inout) :: t2(:,:,:)
    real(wp), intent(inout) :: t2a(:,:)
    real(wp), intent(out  ) :: tskina(:,:)

    logical, intent(inout) :: error

    integer :: i, j, n
    real(wp) :: dwdt, q2sat, qsat, rr, deba, dtdt, dddt, dcdt, frsnw, heff, rh, tam_zs
    real(wp) :: convwtr_slope
    real(wp) :: prc_tmp, prcw_tmp, prcs_tmp
    real(wp) :: prc_ocn, prc_ocn_conv, prc_ocn_wcon
    real(wp) :: prc_lnd, prc_lnd_conv, prc_lnd_wcon
    real(wp) :: qold, A_loc, wcon_budget, wcon_in, water_res
    real(wp) :: frocn
    real(wp), dimension(nm) :: rskin

    real(wp), parameter :: A_trop_min = 1.e-3_wp   ! safeguard against division by ~0 in extreme cold profile
    real(wp), parameter :: water_check_tol = 1.e-9_wp  ! kg/m2 per fast step; anything above is real, not round-off
    real(wp), parameter :: atm_heat_tol = 1.e-3_wp     ! W/m2, max allowed atmospheric energy imbalance
    real(wp) :: tam_old_loc
    real(dp) :: e_dh, e_conv, e_rad, e_sha, e_lat   ! global atm energy budget: heat-content change [J] and source terms [W*m2]
    real(dp) :: area_tot
    real(wp) :: res_trans, res_tot, e_net


    ! global atmospheric energy budget accumulators (for check_energy)
    e_dh = 0._dp; e_conv = 0._dp; e_rad = 0._dp; e_sha = 0._dp; e_lat = 0._dp

    !$omp parallel do private(i,j,n,rr,convwtr_slope,prc_tmp,prcw_tmp,prcs_tmp,prc_ocn,prc_ocn_conv,prc_ocn_wcon,prc_lnd,prc_lnd_conv,prc_lnd_wcon) &
    !$omp private(frocn,heff,qold,A_loc,wcon_budget,wcon_in,water_res,dwdt,q2sat,qsat,frsnw,deba,dtdt,dddt,dcdt,tam_zs,rh,rskin,tam_old_loc) &
    !$omp reduction(+:e_dh,e_conv,e_rad,e_sha,e_lat)
    do j=1,jm
      do i=1,im

        !---------------------------
        ! precipitation
        !---------------------------

        rr = (ram(i,j)/rh_max)

        ! moisture convergence due to synoptic activity on slope
        convwtr_slope = c_wrt_3*sqrt(sam(i,j))*slope(i,j)*ra*qam(i,j) + c_wrt_4*max(0._wp,sam2(i,j)-20._wp)*ra*qam(i,j)  ! m/s * kg/m3 * kg/kg = kg/m2/s
          !+ c_wrt_4*max(0._wp,eke(i,j)-20._wp)*ra*qam(i,j)  ! m/s * kg/m3 * kg/kg = kg/m2/s

        ! precipitation from moisture convergence and evaporation
        prc_ocn_conv = max(0._wp,convwtr(i,j)+convwtr_slope+evpa(i,j))*rr
        prc_lnd_conv = prc_ocn_conv

        ! additional precipitation from water content and residence time of water in the atmosphere
        ! ocean
        prc_ocn_wcon = wcon(i,j)*rr/(c_wrt_1*sec_day)
        prc_ocn = prc_ocn_conv+prc_ocn_wcon 
        ! land
        prc_lnd_wcon = wcon(i,j)*rr/(c_wrt_2*sec_day)
        prc_lnd = prc_lnd_conv+prc_lnd_wcon

        ! weighted mean of land and ocean
        frocn = frst(i,j,i_ocn)+frst(i,j,i_sic)
        prc_tmp = frocn*prc_ocn + (1._wp-frocn)*prc_lnd 

        ! diagnostics
        prc_conv(i,j) = frocn*prc_ocn_conv + (1._wp-frocn)*prc_lnd_conv
        prc_wcon(i,j) = frocn*prc_ocn_wcon + (1._wp-frocn)*prc_lnd_wcon

        !-------------------------------------
        ! update prognostic atmospheric humidity
        !-------------------------------------

        ! column water content tendency
        dwdt = evpa(i,j)-prc_tmp+convwtr(i,j)  ! kg/m2/s

        ! snapshot column water before the budget update (for optional per-column check)
        if (check_water) wcon_in = wcon(i,j)

        ! prognostic update of column water content
        wcon(i,j) = wcon(i,j) + dwdt*tstep
        ! save post-budget wcon (pre-cap) for the mass-conservative prc correction
        wcon_budget = wcon(i,j)

        ! invert vesta's linear relation wcon = ram·A_trop + W_strat
        ! to derive surface relative humidity ram from the prognostic wcon
        A_loc = max(A_trop(i,j), A_trop_min)
        ram(i,j) = (wcon(i,j) - W_strat(i,j)) / A_loc

        ! saturated specific humidity at the surface 
        qsat = fqsat(tam(i,j),psa(i,j))
        qam(i,j) = ram(i,j)*qsat

        ! overflow: cap ram at rh_max 
        prc_over(i,j) = 0._wp
        if (ram(i,j).gt.rh_max) then
          ram(i,j) = rh_max
          qam(i,j) = ram(i,j)*qsat
          wcon(i,j) = ram(i,j)*A_loc + W_strat(i,j)
        endif

        ! underflow: cap qam at 1e-6
        if (qam(i,j).lt.1.e-6_wp) then
          qam(i,j) = 1.e-6_wp
          if (qsat.gt.0._wp) ram(i,j) = qam(i,j)/qsat
          wcon(i,j) = ram(i,j)*A_loc + W_strat(i,j)
        endif

        ! mass-consistent precipitation correction: water added/removed by capping is reflected in the precipitation 
        prc_over(i,j) = (wcon_budget - wcon(i,j))/tstep
        prc_tmp = prc_tmp + prc_over(i,j)

        if (prc_tmp.lt.0._wp) then
          print *,'WARNING: prc<0 ',prc_tmp,' in (i,j) ',i,j
          if (prc_tmp.lt.-1.) error = .true.
          prc_tmp = 0._wp
        endif

        ! per-column water-budget consistency check:
        ! (wcon_new - wcon_old) should equal (E - P + CW)·dt 
        if (check_water) then
          water_res = (wcon(i,j) - wcon_in) - (evpa(i,j) - prc_tmp + convwtr(i,j))*tstep
          if (abs(water_res) .gt. water_check_tol) then
            print *,'WARNING: atm water budget residual ',water_res,' kg/m2 in (i,j) ',i,j
          endif
        endif

        ! separate precipitation into rain and snow using 2m temperature,
        ! separately for each macro surface type
        prcw_tmp = 0._wp
        prcs_tmp = 0._wp
        do n=1,nm
          if (frst(i,j,n).gt.0._wp) then
            if (t2(i,j,n).lt.T0-5._wp) then
              frsnw = 1._wp
            else if (t2(i,j,n).gt.T0+5._wp) then
              frsnw = 0._wp
            else 
              frsnw = 0.1_wp*(T0+5._wp-t2(i,j,n))
            endif
            prcw(i,j,n) = prcw(i,j,n) + prc_tmp*(1._wp-frsnw) 
            prcs(i,j,n) = prcs(i,j,n) + prc_tmp*frsnw
            prcw_tmp = prcw_tmp + prc_tmp*(1._wp-frsnw) * frst(i,j,n)
            prcs_tmp = prcs_tmp + prc_tmp*frsnw         * frst(i,j,n)
          endif
        enddo

        ! cumulate precipitation over fast time steps
        prc(i,j) = prc(i,j) + prc_tmp

        !-----------------------
        ! update dust
        !-----------------------

        if (l_dust) then
          ! effective dust height scale
          heff = hdust(i,j)*hatm/(hdust(i,j)+hatm)    ! m
          ! tendency of dust content
          dddt = convdst(i,j) + dust_emis(i,j) - dust_dep(i,j) ! kg/m2/s
          ! new column dust content
          dam(i,j) = dam(i,j) + dddt/(heff*ra)*tstep     ! kg/m2/s / m * m3/kg * s = kg/m2

          if (dam(i,j).lt.0._wp) then
            if (dam(i,j).lt.-1.e-7_wp) print *,'WARNING: dust<0 ',dam(i,j),' in (i,j) ',i,j
            if (dam(i,j).lt.-1._wp) error = .true.
            dam(i,j) = 0._wp
          endif
        endif

        !-----------------------
        ! update CO2 mass mixing ratio
        !-----------------------

        ! tendency of CO2
        dcdt = convco2(i,j) + co2flx(i,j) ! kgCO2/m2/s
        ! new CO2 
        cam(i,j) = cam(i,j) + dcdt/(hatm*ra2a(i,j))*tstep     ! kgCO2/m2/s / m * m3/kg * s = kgCO2/kg

        !-------------------------------------
        ! update atmospheric prognostic temperature
        !-------------------------------------

        ! net energy balance of atmospheric column
        deba = convdse(i,j) + rb_atm(i,j) + sha(i,j) + (cle*prcw_tmp+cls*prcs_tmp)   ! W/m2
        ! temperature tendency
        dtdt=deba/(amas*cv)    ! J/m2/s * m2/kg * kg*K/J = K/s
        ! atmospheric column temperature before the update (for the energy conservation check)
        if (check_energy) tam_old_loc = tam(i,j)
        ! new atmospheric temperature
        tam(i,j)=tam(i,j)+dtdt*tstep

        ! accumulate the global atmospheric energy budget [J] and source terms [W*m2]
        ! convdse is the dry-static-energy transport convergence; rb_atm/sha/latent are the column sources/sinks
        if (check_energy) then
          e_dh   = e_dh   + real(amas*cv*(tam(i,j)-tam_old_loc),dp)*real(sqr(i,j),dp)
          e_conv = e_conv + real(convdse(i,j),dp)*real(sqr(i,j),dp)
          e_rad  = e_rad  + real(rb_atm(i,j),dp)*real(sqr(i,j),dp)
          e_sha  = e_sha  + real(sha(i,j),dp)*real(sqr(i,j),dp)
          e_lat  = e_lat  + real(cle*prcw_tmp+cls*prcs_tmp,dp)*real(sqr(i,j),dp)
        endif

        !-------------------------------------
        ! 2m temperature and humidity 
        !-------------------------------------

        rskin(:) = 0._wp
        do n=1,nm
          if (frst(i,j,n).gt.0._wp) then

            ! atmospheric temperature at the elevation of each surface type
            tam_zs = t_prof(zsa(i,j), zs(i,j,n), tam(i,j), gams(i,j), gamb(i,j), gamt(i,j), 30.e3_wp, 0)

            ! 2m temperature
            t2(i,j,n) = 0.5_wp*(tskin(i,j,n)+tam_zs)

            if (n.eq.i_ocn .or. tskin(i,j,n).gt.T0) then
              ! saturation over water
              q2sat = q_sat_w(t2(i,j,n),ps(i,j,n))
              qsat  = q_sat_w(tskin(i,j,n),ps(i,j,n))
            else
              ! saturation over ice
              q2sat = q_sat_i(t2(i,j,n),ps(i,j,n))
              qsat  = q_sat_i(tskin(i,j,n),ps(i,j,n))
            endif

            ! skin relative humidity keeping constant specific humidity (qam)
            rskin(n) = min(qsat,qam(i,j))/qsat
            if (n.eq.i_ocn) then
              rskin(n) = max(rskin(n),rskin_ocn_min)
            endif

            ! rh as average between ram and rskina
!            if (n.eq.i_ocn) then
!              rh = 0.5_wp*(rh_max+rskin(n)) 
!            else
              rh = 0.5_wp*(ram(i,j)+rskin(n)) 
!            endif
            q2(i,j,n) = rh*q2sat

            r2(i,j,n) = q2(i,j,n)/q2sat

          endif
        enddo

        ! grid cell averages
        t2a(i,j) = sum(t2(i,j,:)*frst(i,j,:))
        q2a(i,j) = sum(q2(i,j,:)*frst(i,j,:))
        r2a(i,j) = sum(r2(i,j,:)*frst(i,j,:))
        tskina(i,j) = sum(tskin(i,j,:)*frst(i,j,:)) 
        rskina(i,j) = sum(rskin(:)*frst(i,j,:)) 

        !-------------------------------------
        ! sea-level temperature 
        !-------------------------------------

        ! reduce temperature to sea level
        if (i_tslz.eq.1) then
          tsl(i,j)=t_prof(zsa(i,j), 0._wp, tam(i,j), gams(i,j), gamb(i,j), gamt(i,j), 30.e3_wp, 0)
        else if (i_tslz.eq.2) then
          tsl(i,j)=t_prof(zsa(i,j), 0._wp, tam(i,j), gams(i,j), gamb(i,j), gamt(i,j), 30.e3_wp, 0) &
            + (gamb(i,j)-gams(i,j))*(1500._wp*(frst(i,j,i_ocn)+frst(i,j,i_sic)))
        endif

        ! sea level temperature for azonal sea level pressure, using skin temperature
        if (i_tsl.eq.1) then
          tsksl(i,j)=max(tam(i,j),tskina(i,j))+c_tsl_gam*zsa(i,j) 
          if (rb_sur(i,j).gt.0._wp) then
            tsksl(i,j)=tsksl(i,j) + (frst(i,j,i_ocn)+frst(i,j,i_sic)+frst(i,j,i_lake))*(gams_max_ocn-gams(i,j))*hgams
          endif
        else if (i_tsl.eq.2) then
          tsksl(i,j) = max(tam(i,j),tskina(i,j)) &
            - max(frst(i,j,i_lnd)*tsl_gams_min_lnd+frst(i,j,i_ice)*tsl_gams_min_ice,gams(i,j))*hgams &
            + ((1._wp-frst(i,j,i_ice))*c_tsl_gam + frst(i,j,i_ice)*c_tsl_gam_ice)*zsa(i,j) 
        else 
          stop 'itsl'
        endif

        !-------------------------------------
        ! check for NaNs

        if (t2a(i,j).ne.t2a(i,j)) error=.true.
        if (q2a(i,j).ne.q2a(i,j)) error=.true.

      enddo
    enddo
    !$omp end parallel do

    !-------------------------------------
    ! global atmospheric energy conservation check
    !-------------------------------------
    ! The column temperature is updated explicitly by deba = convdse + rb_atm + sha + latent,
    ! so the per-column budget closes by construction. The meaningful global tests are:
    !  - the dry-static-energy transport (convdse) must integrate to zero (transport conserves energy)
    !  - the total heat-content change must equal the net column energy input (update consistency)
    if (check_energy) then
      area_tot = sum(real(sqr,dp))
      ! global-mean heat-transport convergence [W/m2]; should be ~0 if the transport conserves energy
      res_trans = real(e_conv/area_tot,wp)
      ! global-mean consistency residual [W/m2]: heat-content change minus net energy input
      res_tot = real((e_dh/tstep-(e_conv+e_rad+e_sha+e_lat))/area_tot,wp)
      ! global-mean net energy input to the atmospheric column [W/m2]
      e_net = real((e_conv+e_rad+e_sha+e_lat)/area_tot,wp)
      if (abs(res_trans).gt.atm_heat_tol .or. abs(res_tot).gt.atm_heat_tol) then
        print *
        print *,'WARNING: atmosphere energy not conserved! year, doy ',year, doy
        print *,'heat transport convergence [W/m2] ',res_trans
        print *,'update consistency residual[W/m2] ',res_tot
        print *,'net energy input to column [W/m2] ',e_net
        print *,'  radiative                [W/m2] ',real(e_rad/area_tot,wp)
        print *,'  sensible heat            [W/m2] ',real(e_sha/area_tot,wp)
        print *,'  latent heat (precip)     [W/m2] ',real(e_lat/area_tot,wp)
        print *,'  transport (convdse)      [W/m2] ',res_trans
      endif
    endif

    return

  end subroutine time_step

end module time_step_mod

