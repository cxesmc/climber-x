!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : d i f f u s i o n _ m o d
!
!  Purpose : ocean tracer diffusion
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Neil R. Edwards and Matteo Willeit
!
! This file is part of CLIMBER-X.
!
! This file was ported from the original c-GOLDSTEIN model,
! see Edwards and Marsh (2005)
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
module diffusion_mod

  use precision, only : wp
  use constants, only : r_earth
  use tridiag, only : tridiag_solve
  use ocn_params, only : dt, diff_gm, diffx_max, diffy_max, l_hdiff_bl, i_taper, slope_width, slope_max
  use ocn_grid, only : k1, maxi, maxj, maxk, mask_ocn, mask_c, mask_u, mask_v, mask_w, dx, dxv, dy, dz, dza, rdx, rdy, rdza, zw, zro, ocn_area

  implicit none

  real(wp), parameter :: eps = 1.e-15_wp

  private
  public :: diffusion, slope_taper

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  d i f f u s i o n
  !   Purpose    :  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine diffusion(i_diff, f_ocn,tracer,diff_iso,diff_dia,drho_dx,drho_dy,drho_dz,slope_crit,mld, &
                      fdx, fdy, fdz, tracer_tendency)

    implicit none

    integer, intent(in) :: i_diff
    real(wp), dimension(:,:), intent(in) :: f_ocn
    real(wp), dimension(:,:,:), intent(in) :: tracer
    real(wp), intent(in) :: diff_iso
    real(wp), dimension(:,:,:), intent(in) :: diff_dia
    real(wp), dimension(:,:,:), intent(in) :: drho_dx, drho_dy, drho_dz
    real(wp), dimension(:,:), intent(in) :: slope_crit
    real(wp), dimension(:,:), intent(in) :: mld   ! mixed layer depth of the previous time step (negative) [m]

    real(wp), dimension(0:,0:,0:), intent(out) :: fdx
    real(wp), dimension(0:,0:,0:), intent(out) :: fdy
    real(wp), dimension(0:,0:,0:), intent(out) :: fdz
    real(wp), dimension(:,:,:), intent(out) :: tracer_tendency

    integer :: i, j, k, i0, j0, i1, i2, ip1, nnp, knp, k1, k2
    integer :: jp, jm, j0c
    real(wp) :: dts_dx, dts_dy, dts_dz, tv
    real(wp) :: slope2, slope_x, slope_y, slope_x_avg, slope_y_avg, drhodx, drhody, drhodz, tracer_ijk
    real(wp) :: taper, mc_x, yfac
    real(wp) :: taper_h, tap_sum, tap_n   ! face-mean slope taper of the horizontal diagonal terms (A11, A22)
    real(wp) :: a33                      ! A33 coefficient on the w-face, K*S^2 averaged over the 4 corners [m2/s]
    real(wp), allocatable, dimension(:,:,:) :: k33
    real(wp), dimension(maxk) :: ta, tb, tc, tr, tx   ! tridiagonal system of the implicit A33 solve
    integer :: kb, n, kk
    real(wp) :: dx_dz, dy_dz
    real(wp), parameter :: f_min = 1.e-30_wp   ! floor for f_ocn in branchless denominator (land f_ocn=0)


    if (i_diff.eq.1) allocate(k33(maxi,maxj,maxk), source=0._wp)


    do k=1,maxk
      do j=1,maxj
        do i=1,maxi

          tracer_ijk = tracer(i,j,k)

          ! flux to east
          if (mask_u(I,j,k).eq.1) then  
            ip1 = modulo(i,maxi) + 1
            dy_dz = dy*dz(k)
            taper_h = 1._wp
            if (i_diff.eq.1) then
              tap_sum = 0._wp
              tap_n = 0._wp
              do knp=0,1
                do nnp=0,1
                  i0 = modulo(i+nnp-1, maxi) + 1
                  k1 = max(1, k+knp-1)
                  if (mask_w(i0,j,k1).eq.1) then
                    slope_x = drho_dx(I,j,k)/min(-eps,drho_dz(i0,j,k1))
                    tap_sum = tap_sum + slope_taper(slope_x**2,slope_crit(j,k))
                    tap_n = tap_n + 1._wp
                  endif
                enddo
              enddo
              if (tap_n.gt.0._wp) taper_h = tap_sum/tap_n
              ! boundary layer: no taper inside the mixed layer of both columns
              if (l_hdiff_bl) then
                if (zro(k).ge.max(mld(i,j),mld(ip1,j))) taper_h = 1._wp
              endif
            endif
            fdx(I,j,k) = - rdx(j) * taper_h * min(diffx_max(j),diff_iso) * (tracer(ip1,j,k)-tracer_ijk) *dy_dz*dt  ! volume flux
            if (i_diff.eq.1 .and. diff_iso.ne.diff_gm) then
              drhodx = drho_dx(I,j,k)
              tv = 0._wp
              do knp=0,1
                do nnp=0,1
                  i0 = modulo(i+nnp-1, maxi) + 1   ! wraps i+nnp = maxi+1 back to 1
                  k1 = max(1, k+knp-1)             ! when k+knp-1=0, k1==k2 below => dts_dz=0 (corner vanishes)
                  k2 = k+knp
                  dts_dz = (tracer(i0,j,k2)-tracer(i0,j,k1)) * rdza(k1)  ! K/m or psu/m
                  slope_x = drhodx/min(-eps,drho_dz(i0,j,k1))
                  taper = slope_taper(slope_x**2,slope_crit(j,k))
                  ! mask_w as a multiply (0 on closed face)
                  tv = tv + real(mask_w(i0,j,k1),wp) * taper*(min(diffx_max(j),diff_iso)-min(diffx_max(j),diff_gm))*slope_x*dts_dz
                enddo
              enddo
              tv = 0.25_wp*tv *dy_dz*dt  ! m/s*K *m2*s = m3 * K
              fdx(I,j,k) = fdx(I,j,k) + tv
            endif
          else
            fdx(I,j,k) = 0._wp
          endif

          ! flux to north
          if (j.lt.maxj .and. mask_v(i,J,k).eq.1) then 
            dx_dz = dxv(J)*dz(k)
            ! tapered horizontal diagonal term (A22), see above
            taper_h = 1._wp
            if (i_diff.eq.1) then
              tap_sum = 0._wp
              tap_n = 0._wp
              do knp=0,1
                do nnp=0,1
                  j0 = j+nnp
                  k1 = max(1, k+knp-1)
                  if (mask_w(i,j0,k1).eq.1) then
                    slope_y = drho_dy(i,J,k)/min(-eps,drho_dz(i,j0,k1))
                    tap_sum = tap_sum + slope_taper(slope_y**2,slope_crit(j0,k))
                    tap_n = tap_n + 1._wp
                  endif
                enddo
              enddo
              if (tap_n.gt.0._wp) taper_h = tap_sum/tap_n
              ! boundary layer: no taper inside the mixed layer of both columns
              if (l_hdiff_bl) then
                if (zro(k).ge.max(mld(i,j),mld(i,j+1))) taper_h = 1._wp
              endif
            endif
            fdy(i,J,k) = - rdy * taper_h * min(diffy_max(J),diff_iso) * (tracer(i,j+1,k)-tracer_ijk) *dx_dz*dt
            if (i_diff.eq.1 .and. diff_iso.ne.diff_gm) then
              ! add isoneutral diffusion
              drhody = drho_dy(i,J,k)
              tv = 0._wp
              do knp=0,1
                do nnp=0,1
                  j0 = j+nnp
                  k1 = max(1, k+knp-1)             ! when k+knp-1=0, k1==k2 below => dts_dz=0 (corner vanishes)
                  k2 = k+knp
                  dts_dz = (tracer(i,j0,k2)-tracer(i,j0,k1)) * rdza(k1)  ! K/m or psu/m
                  slope_y = drhody/min(-eps,drho_dz(i,j0,k1))
                  taper = slope_taper(slope_y**2,slope_crit(j0,k))
                  ! mask_w as a multiply (0 on closed face)
                  tv = tv + real(mask_w(i,j0,k1),wp) * taper*(min(diffy_max(j),diff_iso)-min(diffy_max(j),diff_gm))*slope_y*dts_dz
                enddo
              enddo
              tv = 0.25_wp*tv *dx_dz*dt  ! m/s*K *m2*s = m3 * K
              fdy(i,J,k) = fdy(i,J,k) + tv
            endif
          else
            fdy(i,J,k) = 0._wp
          endif

          ! flux up
          if (mask_w(i,j,K).eq.1) then
            ! z-derivative of tracer field
            dts_dz = (tracer(i,j,k+1) - tracer_ijk)*rdza(k)  ! K/m or psu/m            
            fdz(i,j,K) = - diff_dia(i,j,K) * dts_dz *ocn_area(i,j)*dt
            if (i_diff.eq.1) then
              ! add isoneutral diffusion
              drhodz = min(-eps,drho_dz(i,j,K)) ! negative
              ! compute horizontal derivatives of tracer
              tv = 0._wp
              a33 = 0._wp
              slope_x_avg = 0._wp
              slope_y_avg = 0._wp
              do knp=0,1
                do nnp=0,1
                  ! phi derivative of tracer field (branchless: mask_c as a multiply, i1/i2 wrapped via modulo)
                  i0 = modulo(i-2+2*nnp, maxi) + 1   ! phi mask point, wraps 0 -> maxi and maxi+1 -> 1
                  i1 = modulo(i+nnp-2, maxi) + 1     ! = i+nnp-1, wrapped 0 -> maxi
                  i2 = modulo(i+nnp-1, maxi) + 1     ! = i+nnp,   wrapped maxi+1 -> 1
                  mc_x = real(mask_c(i0,j,k+knp),wp)
                  dts_dx = mc_x * (tracer(i2,j,k+knp)-tracer(i1,j,k+knp)) * rdx(j)
                  slope_x = mc_x * drho_dx(i1,j,k+knp)/drhodz
                  ! theta-derivative (branchless: clamp j indices, zero the contribution outside [1,maxj])
                  j0 = j-1+2*nnp
                  j0c = min(max(j0,1),maxj)          ! clamped theta mask index (only used when in-bounds)
                  jp  = min(j+nnp, maxj)
                  jm  = max(j+nnp-1, 1)
                  yfac = merge(real(mask_c(i,j0c,k+knp),wp), 0._wp, j0.ge.1 .and. j0.le.maxj)
                  dts_dy = yfac * (tracer(i,jp,k+knp)-tracer(i,jm,k+knp)) * rdy
                  slope_y = yfac * drho_dy(i,jm,k+knp)/drhodz
                  ! sum up over 4 corners
                  slope2 = slope_x**2+slope_y**2+eps
                  taper = slope_taper(slope2,slope_crit(j,K))
                  tv = tv + (diff_iso+diff_gm)*taper * (slope_x*dts_dx + slope_y*dts_dy)
                  ! A33 term, solved implicitly below
                  a33 = a33 + diff_iso*taper*slope2
                  ! average slopes
                  slope_x_avg = slope_x_avg + slope_x
                  slope_y_avg = slope_y_avg + slope_y
                enddo
              enddo
              fdz(i,j,K) = fdz(i,j,K) + 0.25_wp*tv *ocn_area(i,j)*dt  ! m/s*K *m2*s = m3 * K
              k33(i,j,K) = 0.25_wp*a33
            endif
          else
            fdz(i,j,K) = 0._wp
          endif

        enddo
      enddo
    enddo


    ! western boundary fluxes, periodic boundary condition
    fdx(0,:,:) = fdx(maxi,:,:)

    ! southern boundary fluxes
    fdy(:,0,:) = 0._wp

    ! bottom boundary fluxes
    fdz(:,:,0) = 0._wp


    do k=1,maxk
      do j=1,maxj
        do i=1,maxi
          ! branchless mask_c; guard f_ocn (=0 on land) in the denominator
          tracer_tendency(i,j,k) = - real(mask_c(i,j,k),wp) &
                                 * (fdx(I,j,k)-fdx(I-1,j,k) + fdy(i,J,k)-fdy(i,J-1,k) + fdz(i,j,K)-fdz(i,j,K-1)) &   ! m3 * K
                                 / (dx(j)*dy*dz(k)*max(f_ocn(i,j),f_min) * dt)
        enddo
      enddo
    enddo

    ! implicit A33 (backward Euler, per column): dz dT/dt = d/dz(K33 dT/dz), no flux at surface and bottom;
    ! the change is added to the tendency, so it is booked as diffusion like the explicit terms
    if (i_diff.eq.1) then
      do j=1,maxj
        do i=1,maxi
          if (mask_c(i,j,maxk).eq.1) then
            kb = maxk
            do kk=maxk,1,-1
              if (mask_c(i,j,kk).eq.1) kb = kk
              if (mask_c(i,j,kk).eq.0) exit
            enddo
            n = maxk-kb+1
            if (n.ge.2) then
              do kk=kb,maxk
                k = kk-kb+1
                ta(k) = 0._wp
                tc(k) = 0._wp
                if (kk.gt.kb) ta(k) = -dt*k33(i,j,kk-1)/(dz(kk)*dza(kk-1))   ! face below cell kk
                if (kk.lt.maxk) tc(k) = -dt*k33(i,j,kk)/(dz(kk)*dza(kk))     ! face above cell kk
                tb(k) = 1._wp - ta(k) - tc(k)
                tr(k) = tracer(i,j,kk)
              enddo
              call tridiag_solve(ta(1:n),tb(1:n),tc(1:n),tr(1:n),tx(1:n),n)
              do kk=kb,maxk
                tracer_tendency(i,j,kk) = tracer_tendency(i,j,kk) + (tx(kk-kb+1)-tracer(i,j,kk))/dt
              enddo
            endif
          endif
        enddo
      enddo
      deallocate(k33)
    endif

   return

  end subroutine diffusion


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function :  s l o p e _ t a p e r
  !   Purpose  :  taper of the isoneutral diffusion and GM tensor for the squared slope slope2
  !               and the critical slope slope_c.
  !               i_taper==1: (Sc/S)^2 for S>Sc (Gerdes et al 1991); A33 = K*taper*S^2 stays at K*Sc^2
  !                           however steep the slope, i.e. a vertical diffusivity K*Sc^2 in near-neutral columns
  !               i_taper==2: 0.5*(1+tanh((Sc-S)/Sw)) (Danabasoglu and McWilliams 1995); all terms vanish for S >> Sc
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  elemental function slope_taper(slope2, slope_c) result(taper)

    implicit none

    real(wp), intent(in) :: slope2, slope_c
    real(wp) :: taper

    if (i_taper.eq.2) then
      taper = 0.5_wp*(1._wp+tanh((slope_c-sqrt(slope2))/slope_width))
    else
      taper = min(1._wp,slope_c**2/(slope2+eps))
    endif

  end function slope_taper

end module diffusion_mod
