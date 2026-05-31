!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : l a k e _ t e m p _ m o d
!
!  Purpose : lake temperature
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
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
module lake_temp_mod

  use precision, only : wp
  use constants, only : rho_w, Lf, g, T0, cap_w, cap_i
  use control, only : check_energy, check_water
  use lnd_grid, only : z_l, z_int_l, dz_l, nl_l, rdz_pos_l, rdz_neg_l
  use lnd_params, only : dt, rdt, dz_lake_ice
  use lake_convection_mod, only : lake_convection
  use tridiag, only : tridiag_solve

  implicit none

  private
  public :: lake_temp

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l a k e _ t e m p
  !   Purpose    :  compute lake temperature and phase changes 
  !              :  by solving the tridiagonal system
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lake_temp(mask_snow, h_snow, h_lake, cap_lake, lambda_int, &
                      flx_g, dflxg_dT, flx_melt, &
                      t_lake, w_snow, w_w, w_i, f_i_lake, f_lake_ice, snowmelt, &
                      t_lake_old, w_snow_old, &
                      h_lake_conv, h_lake_mix, &
                      energy_cons_lake,i,j)

    implicit none

    integer, intent(in) :: mask_snow
    real(wp), intent(in) :: h_snow
    real(wp), intent(in) :: h_lake
    real(wp), intent(in) :: flx_g, dflxg_dT, flx_melt
    real(wp), dimension(0:), intent(in) :: cap_lake, lambda_int
    real(wp), intent(inout) :: w_snow
    real(wp), dimension(:), intent(inout) :: w_w, w_i, f_i_lake
    real(wp), intent(inout) :: f_lake_ice
    real(wp), dimension(0:), intent(inout) :: t_lake
    real(wp), intent(out) :: snowmelt
    real(wp), intent(inout) :: w_snow_old
    real(wp), dimension(0:), intent(out) :: t_lake_old
    real(wp), intent(out) :: h_lake_conv
    real(wp), intent(out) :: h_lake_mix
    real(wp), intent(out) :: energy_cons_lake

    integer :: i, j, k
    real(wp), dimension(0:nl_l) :: a, b, c, r, x
    real(wp), dimension(0:nl_l) :: cap, rcap
    real(wp), dimension(0:nl_l) :: dz_loc, z_loc, rdz_loc_pos, rdz_loc_neg
    real(wp), dimension(nl_l) :: w, w_w_old, w_i_old
    real(wp) :: t_freeze, t_freeze_lake
    real(wp) :: w_melt, w_ice, w_liq, w_liq_max, H, H_star, w_snow_tmp, dw_snow, energy_warm_snow, flx_excess
    real(wp) :: dz_l_sum
    real(wp) :: i_p, j_p, n_p, sum_water
    logical :: melt, i_print

    i_print = .false.
    n_p = 1
    i_p = 21
    j_p = 33
    i_print = i_print .and. i.eq.i_p .and. j.eq.j_p

    ! set snowmelt to zero
    snowmelt = 0._wp
    energy_warm_snow = 0._wp

    dz_loc(0) = h_snow
    z_loc(0) = -0.5_wp * dz_loc(0)
    z_loc(1:nl_l-1) = z_l(1:nl_l-1)
    z_loc(nl_l) = (h_lake+0.5_wp*z_loc(nl_l-1))/1.5_wp 

    dz_loc(1) = 0.5_wp * ( z_loc(1) + z_loc(2) )
    do k=2,nl_l-1
     dz_loc(k) = 0.5_wp * ( z_loc(k+1) - z_loc(k-1) )
    enddo
    dz_loc(nl_l) = z_loc(nl_l) - z_loc(nl_l-1)

    rdz_loc_pos(1:nl_l-2) = rdz_pos_l(1:nl_l-2)
    rdz_loc_neg(2:nl_l-1) = rdz_neg_l(2:nl_l-1)
    rdz_loc_pos(0) = 1._wp / (z_loc(1) - z_loc(0))
    rdz_loc_neg(1) = 1._wp / (z_loc(1) - z_loc(0))
    rdz_loc_pos(nl_l-1) = 1._wp / (z_loc(nl_l) - z_loc(nl_l-1))
    rdz_loc_neg(nl_l) = 1._wp / (z_loc(nl_l) - z_loc(nl_l-1))

    ! assign some arrays
    cap = cap_lake*dz_loc  ! J/m2/K
    rcap(1:nl_l) = 1._wp/cap(1:nl_l)
    if( cap(0).gt.0._wp ) rcap(0) = 1._wp/cap(0)

    ! total water equivalent in each layer 
    w(1:nl_l) = dz_loc(1:nl_l) * rho_w  ! kg/m2

    w_i(1:nl_l) = w(1:nl_l)*f_i_lake(1:nl_l)
    w_w(1:nl_l) = w(1:nl_l)-w_i(1:nl_l)
    
    ! save old prognostic variables
    t_lake_old = t_lake
    w_w_old = w_w
    w_i_old = w_i
    w_snow_old = w_snow

    ! freezing temperature of lake water as function of lake salinity, from Marsh et al 2011, from Millero 1978
    t_freeze_lake = T0
    !t_freeze_lake = -0.0575_wp*sal+0.0017_wp*sal**1.5_wp-0.0002_wp*sal**2 + T0  ! K TODO

    if( mask_snow .eq. 1 ) then 

      ! top layer, snow, k=0
      k = 0
      a(k) = 0._wp
      b(k) = 1._wp + dt * rcap(k) * ( lambda_int(k) * rdz_loc_pos(k) - dflxg_dT )
      c(k) = - dt * rcap(k) * lambda_int(k) * rdz_loc_pos(k)
      r(k) = t_lake(k) + dt * rcap(k) * ( flx_g - dflxg_dT*t_lake(k) )

      ! intermediate layers
      do k=1,nl_l-1
        a(k) = - dt * rcap(k) * lambda_int(k-1) * rdz_loc_neg(k)
        b(k) = 1._wp + dt * rcap(k) * ( lambda_int(k-1) * rdz_loc_neg(k) + lambda_int(k) * rdz_loc_pos(k) )
        c(k) = - dt * rcap(k) * lambda_int(k) * rdz_loc_pos(k)
        r(k) = t_lake(k) 
      enddo

      ! bottom layer, k=nl_l
      k = nl_l
      a(k) = - dt * rcap(k) * lambda_int(k-1) * rdz_loc_neg(k)
      b(k) = 1._wp + dt * rcap(k) * lambda_int(k-1) * rdz_loc_neg(k)
      c(k) = 0._wp
      r(k) = t_lake(k) 

      ! solve tridiagonal system
      call tridiag_solve(a,b,c,r,x,nl_l+1)
      ! assign new lake temperature
      t_lake = x


      ! PHASE CHANGE

      do k=0,nl_l

        if( k .eq. 0 ) then 
          t_freeze = T0
          w_ice = w_snow
          w_liq = 0._wp
          w_liq_max = -1._wp
        else
          t_freeze = t_freeze_lake
          w_ice = w_i(k)
          w_liq = w_w(k)
          if( t_lake(k) .lt. t_freeze ) then
            w_liq_max = 0._wp 
          else
            w_liq_max = dz_loc(k) * rho_w 
          endif
        endif

        ! freezing or melting occuring
        if( (t_lake(k).gt.t_freeze .and. w_ice.gt.0._wp) .or. (k.ne.0 .and. t_lake(k).lt.t_freeze .and. w_liq.gt.w_liq_max) ) then
          if( (t_lake(k).gt.t_freeze .and. w_ice.gt.0._wp) ) then
            melt = .true.
          else
            melt = .false.
          endif

          if( k .eq. 0 ) then ! snow layer
            H = (dflxg_dT -cap(k)*rdt) * (t_freeze - t_lake(k)) ! W/m2
          else ! lake layers
            H = -cap(k)*rdt*(t_freeze - t_lake(k))
          endif

          w_melt = H*dt/Lf  ! kg/m2

          ! melting
          if( w_melt .gt. 0._wp .and. melt ) then
            if( k .eq. 0 ) then
              t_lake(k) = t_freeze
              if (w_melt.le.w_ice) then
                snowmelt = w_melt*rdt  ! kg/m2/s
                w_snow = w_ice-w_melt
              else ! not enough snow to melt
                snowmelt = w_ice*rdt  ! kg/m2/s
                w_snow = 0._wp
                H_star = Lf*(w_melt-w_ice)*rdt
                ! add excess energy to first lake layer
                t_lake(1) = t_lake(1) + dt*rcap(1)*H_star
              endif
            else
              if (w_melt.le.w_ice) then
                w_i(k) = w_ice-w_melt
                t_lake(k) = t_freeze
              else ! not enough ice to melt
                w_i(k) = 0._wp
                H_star = Lf*(w_melt-w_ice)*rdt
                t_lake(k) = t_freeze + dt*rcap(k)*H_star
              endif
            endif

          ! freezing
          else if( w_melt.lt.0._wp .and. (.not.melt) ) then

            if( w_liq+w_ice .ge. w_liq_max ) then
              w_i(k) = min(w_liq+w_ice-w_liq_max, w_ice-w_melt) 
            else
              w_i(k) = 0._wp
            endif

            ! determine energy not used/released during melt/freezing processes and adjust temperature if necessary
            H_star = H - Lf * (w_ice-w_i(k)) * rdt
            t_lake(k) = t_freeze + dt*rcap(k) * H_star

          endif

          ! update liquid water content
          if( k .ne. 0 ) w_w(k) = w_liq+w_ice - w_i(k)

        endif

        ! use flx_melt to directly melt part of the snow layer 
        ! without the need for the whole snow layer to reach melting point
        if (k.eq.0 .and. flx_melt.gt.0._wp) then
          dw_snow = flx_melt / (rdt*((t_freeze-t_lake(0))*cap_i+Lf))
          if (dw_snow.le.w_snow) then
            ! enough snow to melt
            snowmelt = snowmelt + dw_snow*rdt  ! kg/m2/s
            w_snow = w_snow-dw_snow
            energy_warm_snow = rdt*(t_freeze-t_lake(0))*cap_i*dw_snow ! W/m2
          else
            ! not enough snow to melt
            ! melt all
            snowmelt = snowmelt + w_snow*rdt  ! kg/m2/s
            ! determine energy not used during melt and adjust temperature
            flx_excess = flx_melt + dflxg_dT*(t_freeze-t_lake(0)) - w_snow*rdt*Lf - rdt*(t_freeze-t_lake(0))*cap_i*w_snow  ! W/m2
            w_snow = 0._wp
            t_lake(0) = t_freeze
            ! add excess energy to first lake layer
            t_lake(1) = t_lake(1) + dt*rcap(1) * flx_excess
          endif
        endif

      enddo


    else ! w_snow lower than w_snow_crit, no explicit snow layer

      ! top layer, k=1
      k = 1
      a(k) = 0._wp
      b(k) = 1._wp + dt * rcap(k) * ( lambda_int(k) * rdz_loc_pos(k) - dflxg_dT )
      c(k) = - dt * rcap(k) * lambda_int(k) * rdz_loc_pos(k)
      r(k) = t_lake(k) + dt * rcap(k) * ( flx_g - dflxg_dT*t_lake(k) ) 

      ! intermediate layers
      do k=2,nl_l-1
        a(k) = - dt * rcap(k) * lambda_int(k-1) * rdz_loc_neg(k)
        b(k) = 1._wp + dt * rcap(k) * ( lambda_int(k-1) * rdz_loc_neg(k) + lambda_int(k) * rdz_loc_pos(k) )
        c(k) = - dt * rcap(k) * lambda_int(k) * rdz_loc_pos(k)
        r(k) = t_lake(k) 
      enddo

      ! bottom layer, k=nl_l
      k = nl_l
      a(k) = - dt * rcap(k) * lambda_int(k-1) * rdz_loc_neg(k)
      b(k) = 1._wp + dt * rcap(k) * lambda_int(k-1) * rdz_loc_neg(k)
      c(k) = 0._wp
      r(k) = t_lake(k) 

      ! solve tridiagonal system
      call tridiag_solve(a(1:nl_l),b(1:nl_l),c(1:nl_l),r(1:nl_l),x(1:nl_l),nl_l)
      ! assign new lake temperature
      t_lake(1:nl_l) = x(1:nl_l)


      ! PHASE CHANGE

      t_freeze = t_freeze_lake

      do k=1,nl_l
        w_ice = w_i(k)
        w_liq = w_w(k)
        ! maximum liquid water content, kg/m2, Niu 2006 (WRONG in CLM technical report)
        if( t_lake(k) .lt. t_freeze ) then
          w_liq_max = 0._wp
        else
          w_liq_max = dz_loc(k) * rho_w
        endif

        ! freezing or melting occuring
        if( (t_lake(k).gt.t_freeze .and. w_ice.gt.0._wp) .or. (t_lake(k).lt.t_freeze .and. w_liq.gt.0._wp) ) then
          if( (t_lake(k).gt.t_freeze .and. w_ice.gt.0._wp) ) then
            melt = .true.
          else
            melt = .false.
          endif

          if( k .eq. 1 ) then
            H = dflxg_dT * (t_freeze - t_lake(k)) - cap(k)*rdt*(t_freeze -t_lake(k))  ! W/m2
          else
            ! intermediate lake layers
            H = -cap(k)*rdt*(t_freeze - t_lake(k))
          endif

          ! update ice content
          w_melt = H*dt/Lf  ! kg/m2

          ! melting
          if( w_melt .gt. 0._wp .and. melt ) then

            ! instead of melting ice in the first lake layer, first use it to melt 'invisible' snow
            if( k.eq.1 ) then
              w_snow = max(0._wp, w_snow-w_melt)
              snowmelt = (w_snow_old - w_snow) * rdt  ! kg/m2/s
              H = H - Lf * snowmelt  ! W/m2
              w_melt = H*dt/Lf  ! kg/m2
            endif

            w_i(k) = max(0._wp, w_ice-w_melt)
            ! determine energy not used during melt and adjust temperature if necessary
            H_star = H - Lf * (w_ice - w_i(k)) * rdt
            if( k .eq. 1 ) then
              t_lake(k) = t_freeze + dt*rcap(k) * H_star / (1._wp - dt*rcap(k)*dflxg_dT)
            else
              t_lake(k) = t_freeze + dt*rcap(k) * H_star
            endif

          ! freezing
          else if( w_melt .lt. 0._wp .and. (.not. melt) ) then
            if( (w_liq+w_ice) .ge. w_liq_max ) then
              w_i(k) = min(w_liq+w_ice-w_liq_max, w_ice-w_melt) 
            else
              w_i(k) = 0._wp
            endif
            ! determine energy not released during freezing and adjust temperature if necessary
            H_star = H - Lf * (w_ice - w_i(k)) * rdt
            if( k .eq. 1 ) then
              t_lake(k) = t_freeze + dt*rcap(k) * H_star / (1._wp - dt*rcap(k)*dflxg_dT)
            else
              t_lake(k) = t_freeze + dt*rcap(k) * H_star
            endif

          endif

          ! update liquid water content
          w_w(k) = w_liq+w_ice - w_i(k)

        endif

      enddo

      ! use flx_melt to directly melt (not-explicit) snow (if present)
      if (flx_melt.gt.0._wp) then
        dw_snow = flx_melt / (rdt*Lf) ! snow that can be melted
        if (dw_snow.le.w_snow) then
          ! enough snow to melt
          snowmelt = snowmelt + dw_snow*rdt  ! kg/m2/s
          w_snow = w_snow-dw_snow
        else
          ! not enough snow to melt
          ! melt all
          snowmelt = snowmelt + w_snow*rdt  ! kg/m2/s
          ! determine energy not used during melt
          flx_excess = flx_melt - w_snow*rdt*Lf ! W/m2
          w_snow = 0._wp
          ! use flx_excess to warm top lake layer
          t_lake(1) = (t_lake(1) + dt*rcap(1) * (flx_excess-dflxg_dT*t_lake(1))) / (1._wp-dflxg_dT*dt*rcap(1))
        endif
      endif

      ! if snow present but no explicit snow layer, use possible first lake layer excess energy to melt snow
      if( w_snow.gt.0._wp ) then
        H = dflxg_dT * (t_freeze - t_lake(1)) - cap(1)*rdt*(t_freeze - t_lake(1))

        if( H .gt. 0._wp ) then

          w_melt = H*dt/Lf ! kg/m2
          w_snow_tmp = w_snow
          w_snow = max( 0._wp, w_snow-w_melt )  ! kg/m2
          H_star = H - Lf*rdt * (w_snow_tmp - w_snow)

          w_ice = w_i(1)
          w_i(1) = max(0._wp, w_ice-w_melt)
          w_w(1) = w_w(1)+w_ice - w_i(1)

          t_lake(1) = t_freeze + dt*rcap(1) * H_star / (1._wp - dt*rcap(1) * dflxg_dT)

          snowmelt = snowmelt + (w_snow_tmp - w_snow) * rdt  ! kg/m2/s

        endif
      endif

      ! set temperature of non existing snow layer equal to first lake layer
      t_lake(0) = min(T0, t_lake(1))

    endif

    ! diagnose fraction of ice in layers
    do k=1,nl_l
      f_i_lake(k) = w_i(k)/(w_i(k)+w_w(k))
    enddo

    ! average ice fraction over given layer thickness
    f_lake_ice = 0._wp
    dz_l_sum = 0._wp
    do k=1,nl_l
      if (z_int_l(k).le.dz_lake_ice) then
        dz_l_sum = dz_l_sum + dz_l(k)
        f_lake_ice = f_lake_ice + f_i_lake(k)*dz_l(k)
      !else if (z_int_l(k).gt.dz_lake_ice .and. z_l(k).le.dz_lake_ice) then
      !  dz_l_sum = dz_l_sum + (dz_l(k)-(z_int_l(k)-dz_lake_ice))
      !  f_lake_ice = f_lake_ice + f_i_lake(k)*(dz_l(k)-(z_int_l(k)-dz_lake_ice))
      else if (z_int_l(k).gt.dz_lake_ice .and. z_int_l(k-1).le.dz_lake_ice) then
        dz_l_sum = dz_l_sum + (dz_lake_ice-z_int_l(k-1))
        f_lake_ice = f_lake_ice + f_i_lake(k)*(dz_lake_ice-z_int_l(k-1))
      endif
    enddo
    if (dz_l_sum.gt.0._wp) then
      f_lake_ice = f_lake_ice/dz_l_sum
    else
      f_lake_ice = f_i_lake(1)
    endif
    !print *
    !print *,z_l(1:4)
    !print *,z_int_l(1:4)
    !print *,dz_l(1:4)
    !print *,f_i_lake(1:4)
    !print *,f_lake_ice
    !print *,dz_l_sum

    !----------------------------------------
    ! convective adjustment
    call lake_convection(t_lake(1:nl_l), f_i_lake, dz_loc(1:nl_l), cap(1:nl_l), t_freeze_lake, &
                         h_lake_conv, h_lake_mix)


    ! update ice and water content
    do k=1,nl_l
      w_i(k) = w(k)*f_i_lake(k)
      w_w(k) = w(k)-w_i(k)
    enddo


    if( check_energy ) then
      ! lake energy balance
      if( mask_snow .eq. 1 ) then
        energy_cons_lake = flx_g + dflxg_dT*(t_lake(0)-t_lake_old(0)) + flx_melt  &
        - sum(Lf*rdt * (w_i_old-w_i)) &
        - Lf*snowmelt - energy_warm_snow &
        - sum(cap*rdt * (t_lake-t_lake_old))
      else
        energy_cons_lake =  flx_g + dflxg_dT*(t_lake(1)-t_lake_old(1)) + flx_melt  &
        - sum(Lf*rdt * (w_i_old-w_i)) &
        - Lf*snowmelt &
        - sum(cap(1:nl_l)*rdt * (t_lake(1:nl_l)-t_lake_old(1:nl_l)))
      endif


      if(abs(energy_cons_lake).gt.1.e-8_wp ) then
        print *,' '
        print *,'energy conservation lake ',energy_cons_lake
        print *,'i,j,mask_snow',i,j,mask_snow
        print *,'t_lake_old',t_lake_old
        print *,'t_lake_after_cond',x
        print *,'t_lake_new',t_lake
        print *,'w_i_old',w_i_old
        print *,'w_i_new',w_i
        print *,'f_i_lake',f_i_lake
        if(mask_snow.eq.1) then
          print *,'g,dg_dT',flx_g,dflxg_dT*(t_lake(0)-t_lake_old(0)),flx_g+dflxg_dT*(t_lake(0)-t_lake_old(0))
          print *,'dsnow dT',energy_warm_snow
        else
          print *,'g,dg_dT',flx_g,dflxg_dT*(t_lake(1)-t_lake_old(1)),flx_g+dflxg_dT*(t_lake(1)-t_lake_old(1))
        endif
        print *,'flx_melt',flx_melt
        print *,'energy_internal',cap*rdt * (t_lake-t_lake_old)
        print *,'phase change',Lf*rdt * (w_i_old-w_i)
        print *,'snowmelt energy',Lf*snowmelt
        print *,'snow_new,snow_old,snowmelt',w_snow,w_snow_old,snowmelt*dt
        print *,'lambda_int',lambda_int
        print *,' '
        !       if(abs(energy_cons_lake).gt.10.) stop
      endif
    endif

    if( check_water ) then
      ! lake water balance
      sum_water = sum(w_w - w_w_old + w_i - w_i_old)
      if( abs(sum_water+w_snow - w_snow_old+snowmelt*dt) .gt. 1.d-10 ) then
        print *,'water conservation ',sum_water+snowmelt*dt+w_snow - w_snow_old 
        print *,w_snow,w_snow_old,snowmelt
        print *,w_w(k),w_w_old(k)
        print *,w_i(k),w_i_old(k)
        stop
      endif
    endif


    ! check temperature range
    if((mask_snow.eq.1 .and. (maxval(t_lake).gt.350._wp .or. minval(t_lake).lt.150._wp)) &
      .or. (mask_snow.eq.0 .and. (maxval(t_lake(1:)).gt.350._wp .or. minval(t_lake(1:)).lt.150._wp))) then
      print *,''
      print *,'WARNING t_lake out of range',t_lake
      print *,'i,j',i,j
      print *,'masksnow',mask_snow
      k = 0
      print *,'SNOW layer'
      print *,'mask_snow,snowmelt',mask_snow, snowmelt*dt
      print *,'t_snow, before, after',t_lake_old(k),t_lake(k)
      print *,'w_snow, before, after',w_snow_old,w_snow

      k=1
      print *,'top lake layer'
      print *,'t_lake, before, after',t_lake_old(k),t_lake(k)
      print *,'w_w, before, after',w_w_old(k),w_w(k)
      print *,'w_i, before, after',w_i_old(k),w_i(k)
      print *,'g,dg_dT*dT',flx_g,dflxg_dT,dflxg_dT*(t_lake(1)-t_lake_old(1))
      print *,'flx_melt',flx_melt
      print *,'cap,lambda',cap, lambda_int
      print *,''
      !if(i.eq.54.and.j.eq.25) stop
    endif

    if(i_print) then
      print *,''
      print *,'WARNING t_lake out of range',t_lake
      print *,'i,j',i,j
      print *,'masksnow',mask_snow
      k = 0
      print *,'SNOW layer'
      print *,'mask_snow,snowmelt',mask_snow, snowmelt*dt
      print *,'t_snow, before, after',t_lake_old(k),t_lake(k)
      print *,'w_snow, before, after',w_snow_old,w_snow
      k=1
      print *,'top lake layer'
      print *,'t_lake, before, after',t_lake_old(k),t_lake(k)
      print *,'w_w, before, after',w_w_old(k),w_w(k)
      print *,'w_i, before, after',w_i_old(k),w_i(k)
      print *,'g,dg_dT',flx_g,dflxg_dT*(t_lake(1)-t_lake_old(1))
      print *,'cap,lambda',cap, lambda_int
      print *,''
      print *,'t_lake_old',t_lake_old
      print *,'t_lake_after_cond',x
      print *,'t_lake_new',t_lake
      print *,'w_i_old',w_i_old
      print *,'w_i_new',w_i
      print *,'f_i_lake',f_i_lake
    endif


    return

  end subroutine lake_temp

end module lake_temp_mod

