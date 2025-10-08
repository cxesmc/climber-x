!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s m b _ t e m p _ m o d
!
!  Purpose : soil/ice temperature for SEMIX
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2024 Potsdam Institute for Climate Impact Research,
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
module smb_temp_mod

  use precision, only : wp
  use constants, only : Lf, T0, cap_i, lambda_i, rho_i
  use control, only : check_energy
  use smb_grid, only : z, z_int, dz, rdz_pos, rdz_neg, nl
  use smb_params, only : dt, rdt
  use smb_params, only : snow_par, surf_par
  use tridiag, only : tridiag_solve

  implicit none

  private
  public :: smb_temp

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s m b _ t e m p
  !   Purpose    :  compute soil/ice temperature
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_temp(mask_snow, h_snow, rain, flx_g, dflxg_dT, flx_melt, refreezing_sum, &
                     t_prof, w_snow, snowmelt, icemelt, refreezing, f_rfz_to_snow, t_prof_old, i,j)

    implicit none

    integer, intent(in) :: mask_snow
    real(wp), intent(in) :: h_snow
    real(wp), intent(in) :: rain
    real(wp), intent(in) :: flx_g, dflxg_dT, flx_melt
    real(wp), intent(in) :: refreezing_sum
    real(wp), intent(inout) :: w_snow
    real(wp), dimension(0:), intent(inout) :: t_prof
    real(wp), intent(out) :: snowmelt
    real(wp), intent(out) :: icemelt
    real(wp), intent(out) :: refreezing
    real(wp), intent(out) :: f_rfz_to_snow
    real(wp), dimension(0:), intent(out) :: t_prof_old

    integer :: i, j, k
    real(wp), dimension(0:nl) :: a, b, c, r, x
    real(wp), dimension(0:nl) :: lambda, lambda_int, cap, cap_dz, rcap_dz
    real(wp), dimension(0:nl) :: dz_loc, z_loc, rdz_loc_pos, rdz_loc_neg
    real(wp) :: H, w_melt, H_star, w_snow_old, dw_snow, dw_ice, energy_warm_ice, flx_excess
    real(wp) :: melt_rain, flx_cold_content, refreezing_avail, refreezing_pot
    real(wp) :: refreezing_tot, refreezing_to_snow, refreezing_to_ice
    real(wp) :: refreezing_max, f_refreezing
    real(wp) :: energy_cons
    real(wp) :: i_p, j_p
    logical :: i_print

    i_print = .false.
    i_p = 99
    j_p = 1

    ! initialize to zero
    snowmelt = 0._wp
    icemelt = 0._wp
    refreezing = 0._wp
    f_rfz_to_snow = 0._wp
    energy_warm_ice = 0._wp
    flx_excess = 0._wp

    dz_loc(0) = h_snow
    z_loc(0) = -0.5_wp * dz_loc(0)
    dz_loc(1:nl) = dz(1:nl)
    z_loc(1:nl) = z(1:nl)

    rdz_loc_pos(1:nl-1) = rdz_pos(1:nl-1)
    rdz_loc_neg(2:nl) = rdz_neg(2:nl)
    rdz_loc_pos(0) = 1._wp / (z_loc(1) - z_loc(0))
    rdz_loc_neg(1) = 1._wp / (z_loc(1) - z_loc(0))

    ! heat capacity, J/m3/K and thermal conductivity, W/m/K
    cap(0) = cap_i*snow_par%rho    ! snow
    lambda(0) = snow_par%lambda 
    cap(1:nl) = rho_i * cap_i
    lambda(1:nl) = lambda_i

    ! assign some arrays
    cap_dz = cap*dz_loc ! J/m2/K
    rcap_dz(1:nl) = 1._wp/cap_dz(1:nl)
    if( dz_loc(0) .gt. 0._wp ) rcap_dz(0) = 1._wp/cap_dz(0)

    ! thermal conductivity at the layer interfaces, W/m/K
    do k=0,nl-1
      lambda_int(k) = lambda(k) * lambda(k+1) * (z_loc(k+1) - z_loc(k))  &
        / ( lambda(k) * (z_loc(k+1) - z_int(k)) + lambda(k+1) * (z_int(k) - z_loc(k)) )
    enddo
    lambda_int(nl) = 0._wp  ! no heat flux at the bottom

    ! save old prognostic variables
    t_prof_old = t_prof
    w_snow_old = w_snow

    if (mask_snow.eq.1) then

      ! top layer, snow, k=0
      k = 0
      a(k) = 0._wp
      b(k) = 1._wp + dt * rcap_dz(k) * (lambda_int(k) * rdz_loc_pos(k) - dflxg_dT)
      c(k) = - dt * rcap_dz(k) * lambda_int(k) * rdz_loc_pos(k)
      r(k) = t_prof(k) + dt * rcap_dz(k) * (flx_g - dflxg_dT*t_prof(k))

      ! intermediate layers
      do k=1,nl-1
        a(k) = - dt * rcap_dz(k) * lambda_int(k-1) * rdz_loc_neg(k)
        b(k) = 1._wp + dt * rcap_dz(k) * ( lambda_int(k-1) * rdz_loc_neg(k) + lambda_int(k) * rdz_loc_pos(k) )
        c(k) = - dt * rcap_dz(k) * lambda_int(k) * rdz_loc_pos(k)
        r(k) = t_prof(k)
      enddo

      ! bottom layer, k=nl
      k = nl
      a(k) = - dt * rcap_dz(k) * lambda_int(k-1) * rdz_loc_neg(k)
      b(k) = 1._wp + dt * rcap_dz(k) * lambda_int(k-1) * rdz_loc_neg(k)
      c(k) = 0._wp
      r(k) = t_prof(k)

      ! solve tridiagonal system
      call tridiag_solve(a,b,c,r,x,nl+1)
      ! assign new soil temperature
      t_prof = x

      ! snow melt if snow layer temperature above freezing
      if (t_prof(0).gt.T0) then
        ! energy available for melting
        H = dflxg_dT * (T0-t_prof(0)) - cap_dz(0)*rdt*(T0-t_prof(0))    ! W/m2
        ! snow water equivalent that can be melted
        w_melt = H*dt/Lf    ! kg/m2
        ! actual melt, limited by snow availability
        snowmelt = min(w_snow_old, w_melt) * rdt  ! kg/m2/s
        ! new snow water equivalent
        w_snow = max(0._wp, w_snow_old-w_melt)
        ! determine energy not used during melt and adjust temperature if necessary
        t_prof(0) = T0
        H_star = H - Lf * (w_snow_old - w_snow) * rdt   ! W/m2
        ! add excess snow energy to top ice/soil layer
        t_prof(1) = t_prof(1) + dt*rcap_dz(1) * H_star
      endif

      ! use flx_melt to directly melt part of the snow layer 
      ! without the need for the whole snow layer to reach melting point
      if (flx_melt.gt.0._wp) then
        ! snow that can be melted
        dw_snow = flx_melt / (rdt*((T0-t_prof(0))*cap_i+Lf))    ! kg/m2
        if (dw_snow.le.w_snow) then
          ! enough snow to melt
          snowmelt = snowmelt + dw_snow*rdt  ! kg/m2/s
          w_snow = w_snow-dw_snow
        else
          ! not enough snow to melt
          ! melt all
          snowmelt = snowmelt + w_snow*rdt  ! kg/m2/s
          ! determine energy not used during melt 
          flx_excess = flx_melt - w_snow*rdt*Lf - rdt*(T0-t_prof(0))*cap_i*w_snow  ! W/m2
          w_snow = 0._wp
          ! add excess snow energy to top ice/soil layer
          t_prof(1) = t_prof(1) + dt*rcap_dz(1) * flx_excess
        endif
      endif

      ! refreeze (part of) snow meltwater and/or rain using snow 'cold content'
      ! 'cold content' of snow layer
      flx_cold_content = (cap_i*w_snow)*rdt*(T0-t_prof(0)) !&    ! snow layer, W/m2
                      ! + cap_dz(1)*rdt*(T0-t_prof(1))   ! first ice layer

      ! liquid water available for refreezing (snowmelt + rain)
      melt_rain = (snowmelt + rain)*dt ! kg/m2

      if (flx_cold_content>0._wp .and. melt_rain>0._wp) then

        ! maximum water that can be refrozen, linear dependence on snow thickness
        if (snow_par%i_rfz.eq.1) then
          refreezing_max = snow_par%porosity*w_snow  ! kg/m2
        else if (snow_par%i_rfz.eq.2) then
          refreezing_max = snow_par%porosity*w_snow - refreezing_sum ! kg/m2
          refreezing_max = max(0._wp,refreezing_max)
        endif

        ! fraction of snow meltwater + rain that can be refrozen
        f_refreezing = min(1._wp,refreezing_max/melt_rain)

        ! total liquid water actually available for refreezing 
        refreezing_avail = melt_rain * f_refreezing ! kg/m2

        ! potential water that can be refrozen using the 'cold content'
        refreezing_pot = flx_cold_content/Lf*dt   ! kg/m2

        ! actual water that is refrozen is the minimum between available and potential
        refreezing_tot = snow_par%f_rfz_max * min(refreezing_avail,refreezing_pot)
        refreezing = refreezing + refreezing_tot*rdt  ! kg/m2/s

        ! fraction of refreezing occuring in snow
        f_rfz_to_snow = snow_par%f_rfz_to_snow_max * min(1._wp, w_snow/snow_par%wsnow_crit_rfz)

        ! refreezing to snow
        refreezing_to_snow = f_rfz_to_snow * refreezing_tot ! kg/m2
        ! remaining refreezing water forms superimposed ice
        refreezing_to_ice = refreezing_tot-refreezing_to_snow   ! kg/m2

        !print *
        !print *,'T snow/ice before',t_prof(0:1)
        !print *,'melt+rain',melt_rain
        !print *,'w_snow',w_snow
        !print *,'refreezing_max',refreezing_max
        !print *,'f_refreezing',f_refreezing
        !print *,'refreezing_tot',refreezing_tot
        !print *,'refreezing_avail',refreezing_avail
        !print *,'refreezing_pot',refreezing_pot
        !print *,'f_rfz_to_snow',f_rfz_to_snow
        !print *,'refreezing_to_snow',refreezing_to_snow
        !print *,'refreezing_to_ice',refreezing_to_ice

        ! compute new snow layer temperature
        ! release latent heat to warm snow layer 
        ! additionally warm snow because refreezing water has a temperature of 0 degc
        ! also correct for derivative of ground heat flux
        ! solve implicitely
        t_prof(0) = (refreezing_to_snow*Lf/(cap_i*w_snow*dt) + T0*refreezing_to_snow/(w_snow*dt) - t_prof(0)*dflxg_dT/(cap_i*w_snow) + t_prof(0)*rdt) &
          / (rdt+refreezing_to_snow/(w_snow*dt)-dflxg_dT/(cap_i*w_snow))
        ! account for refreezing to ice and warm the top ice layer with latent heat
        t_prof(1) = t_prof(1) + rcap_dz(1) * refreezing_to_ice*Lf   ! kg/m2/s * J/kg * s * K*m2/J = K

        ! add snow refreezing water to snow water equivalent
        w_snow = w_snow + refreezing_to_snow

        !print *,'T snow/ice after',t_prof(0:1)
        !print *,'w_snow after',w_snow

      endif

      ! ice melt if temperature above freezing
      H = 0._wp
      do k=1,nl
        if (t_prof(k).gt.T0) then
          H = H - cap_dz(k)*rdt*(T0-t_prof(k))    ! W/m2
          t_prof(k) = T0
        endif
      enddo
      icemelt = H/Lf ! kg/m2/s

    else ! w_snow lower than w_snow_crit, no explicit snow layer

      ! top layer, k=1
      k = 1
      a(k) = 0._wp
      b(k) = 1._wp + dt * rcap_dz(k) * (lambda_int(k) * rdz_loc_pos(k) - dflxg_dT)
      c(k) = - dt * rcap_dz(k) * lambda_int(k) * rdz_loc_pos(k)
      r(k) = t_prof(k) + dt * rcap_dz(k) * (flx_g - dflxg_dT*t_prof(k))

      ! intermediate layers
      do k=2,nl-1
        a(k) = - dt * rcap_dz(k) * lambda_int(k-1) * rdz_loc_neg(k)
        b(k) = 1._wp + dt * rcap_dz(k) * ( lambda_int(k-1) * rdz_loc_neg(k) + lambda_int(k) * rdz_loc_pos(k) )
        c(k) = - dt * rcap_dz(k) * lambda_int(k) * rdz_loc_pos(k)
        r(k) = t_prof(k) 
      enddo

      ! bottom layer, k=nl
      k = nl
      a(k) = - dt * rcap_dz(k) * lambda_int(k-1) * rdz_loc_neg(k)
      b(k) = 1._wp + dt * rcap_dz(k) * lambda_int(k-1) * rdz_loc_neg(k)
      c(k) = 0._wp
      r(k) = t_prof(k) 

      ! solve tridiagonal system
      call tridiag_solve(a(1:nl),b(1:nl),c(1:nl),r(1:nl),x(1:nl),nl)
      ! assign new soil temperature
      t_prof(1:nl) = x(1:nl)

      ! use flx_melt to directly melt (not-explicit) snow (if present)
      if (flx_melt.gt.0._wp) then
        dw_snow = flx_melt / (rdt*Lf) ! snow that can be melted
        if (dw_snow.le.w_snow) then
          ! enough snow to melt
          snowmelt = snowmelt + dw_snow*rdt  ! kg/m2/s
          w_snow = w_snow-dw_snow
          flx_excess = 0._wp
        else
          ! not enough snow to melt
          ! melt all
          snowmelt = snowmelt + w_snow*rdt  ! kg/m2/s
          ! determine energy not used during melt
          flx_excess = flx_melt - w_snow*rdt*Lf ! W/m2
          w_snow = 0._wp
        endif
      else
        flx_excess = 0._wp
      endif

      ! ice melt

      ! top layer ice melt
      k = 1
      if (t_prof(k).gt.T0) then
        H = dflxg_dT * (T0-t_prof(k)) - cap_dz(k)*rdt*(T0-t_prof(k)) 
        icemelt = H/Lf
        t_prof(k) = T0
      endif

      ! use flx_excess to directly melt ice 
      ! without the need for the whole ice layer to reach melting point
      if (flx_excess.gt.0._wp) then
        dw_ice = flx_excess / (rdt*(max(0._wp,T0-t_prof(1))*cap_i+Lf)) ! kg/m2, ice that can be melted
        icemelt = icemelt + dw_ice*rdt  ! kg/m2/s
        energy_warm_ice = rdt*max(0._wp,T0-t_prof(1))*cap_i*dw_ice ! W/m2
      endif

      ! additionally melt ice if temperature of any layer above freezing point
      H = 0._wp
      do k=1,nl
        if (t_prof(k).gt.T0) then
          H = H - cap_dz(k)*rdt*(T0-t_prof(k))
          t_prof(k) = T0
        endif
      enddo
      icemelt = icemelt + H/Lf ! kg/m2/s

      ! set temperature of non existing snow layer equal to top ice/soil layer
      t_prof(0) = min(T0, t_prof(1))

    endif


    ! check energy conservation 
    if( check_energy ) then

      ! energy balance
      if (mask_snow.eq.1) then
        energy_cons =  flx_g + dflxg_dT*(t_prof(0)-t_prof_old(0)) + flx_melt  &
          - Lf*(snowmelt+icemelt-refreezing) &
          + (w_snow-w_snow_old)*cap_i*T0*rdt &
          - (cap_i*w_snow*rdt*t_prof(0) - cap_i*w_snow_old*rdt*t_prof_old(0)) &
          - sum(cap_dz(1:nl)*rdt * (t_prof(1:nl)-t_prof_old(1:nl)))
      else
        energy_cons = flx_g + dflxg_dT*(t_prof(1)-t_prof_old(1)) + flx_melt &
          - Lf*(snowmelt+icemelt-refreezing) - energy_warm_ice &
          - sum(cap_dz(1:nl)*rdt * (t_prof(1:nl)-t_prof_old(1:nl)))
      endif

      if(abs(energy_cons).gt.1.e-10_wp ) then
      !if (i.eq.40 .and. j.eq.105) then
        print *,' '
        print *,'energy conservation soil/ice',energy_cons
        print *,'mask_snow',mask_snow
        print *,'t_soil_old',t_prof_old
        print *,'t_soil_after_cond',x
        print *,'t_soil_new',t_prof
        print *,'dg_dT',dflxg_dT
        if(mask_snow.eq.1) then
          print *,'g,dg',flx_g,dflxg_dT*(t_prof(0)-t_prof_old(0))
        else
          print *,'g,dg',flx_g,dflxg_dT*(t_prof(1)-t_prof_old(1))
        endif
        print *,'flx_melt [W/m2]',flx_melt
        print *,'flx_excess',flx_excess
        print *,'flx_cold_content [W/m2]',flx_cold_content
        print *,'snowmelt + rain [kg/m2]',(rain+snowmelt)*dt
        print *,'rain [kg/m2]',rain*dt
        print *,'snowmelt [kg/m2]',snowmelt*dt
        print *,'icemelt [kg/m2]',icemelt*dt
        print *,'refreezing [kg/m2]',refreezing*dt
        print *,'f_refreezing',f_refreezing
        print *,'refreezing_avail,refreezing_pot [kg/m2]',refreezing_avail, refreezing_pot
        print *,'energy_internal [W/m2]',rdt * (cap_dz*t_prof-cap_dz*t_prof_old)
        print *,'snowmelt energy [W/m2]',Lf*snowmelt
        print *,'icemelt energy [W/m2]',Lf*icemelt
        print *,'refreezing energy [W/m2]',Lf*refreezing
        print *,'energy_warm_ice [W/m2]',energy_warm_ice
        print *,'snow_new,snow_old [kg/m2]',w_snow,w_snow_old
        print *,'cap_dz',cap_dz
        print *,'lambda_int',lambda_int
        print *,' '
      if(abs(energy_cons).gt.1._wp ) stop
      endif
    endif

    ! check temperature range
    if(maxval(t_prof) .gt. 350._wp .or. minval(t_prof) .lt. 150._wp) then
      print *,''
      print *,'WARNING t_prof out of range',t_prof
      print *,i,j
      print *,'mask_snow',mask_snow
      print *,'t_soil_old',t_prof_old
      print *,'t_soil_after_cond',x
      print *,'t_soil_new',t_prof
      print *,'dg_dT',dflxg_dT
      if(mask_snow.eq.1) then
        print *,'g,dg',flx_g,dflxg_dT*(t_prof(0)-t_prof_old(0))
      else
        print *,'g,dg',flx_g,dflxg_dT*(t_prof(1)-t_prof_old(1))
      endif
      print *,'flx_melt [W/m2]',flx_melt
      print *,'flx_cold_content [W/m2]',flx_cold_content
      print *,'snowmelt + rain [kg/m2]',(rain+snowmelt)*dt
      print *,'rain [kg/m2]',rain*dt
      print *,'snowmelt [kg/m2]',snowmelt*dt
      print *,'icemelt [kg/m2]',icemelt*dt
      print *,'refreezing [kg/m2]',refreezing*dt
      print *,'f_refreezing',f_refreezing
      print *,'refreezing_avail,refreezing_pot [kg/m2]',refreezing_avail, refreezing_pot
      print *,'energy_internal [W/m2]',rdt * (cap_dz*t_prof-cap_dz*t_prof_old)
      print *,'snowmelt energy [W/m2]',Lf*snowmelt
      print *,'icemelt energy [W/m2]',Lf*icemelt
      print *,'refreezing energy [W/m2]',Lf*refreezing
      print *,'energy_warm_ice [W/m2]',energy_warm_ice
      print *,'snow_new,snow_old [kg/m2]',w_snow,w_snow_old
      print *,'cap_dz',cap_dz
      print *,'lambda_int',lambda_int
      print *,' '
      !stop
    endif

    return

  end subroutine smb_temp


end module smb_temp_mod

