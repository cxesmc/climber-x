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
  use ocn_params, only : dt, l_diff33_impl, diff_gm, diffx_max, diffy_max
  use ocn_grid, only : k1, maxi, maxj, maxk, mask_ocn, mask_c, mask_u, mask_v, mask_w, dx, dxv, dy, dz, dza, rdx, rdy, rdza, zw, ocn_area

  implicit none

  real(wp), parameter :: eps = 1.e-15_wp

  private
  public :: diffusion, diffusion_33

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  d i f f u s i o n
  !   Purpose    :  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine diffusion(i_diff, f_ocn,tracer,diff_iso,diff_dia,drho_dx,drho_dy,drho_dz,slope_crit, &
                      fdx, fdy, fdz, slope2_w, tracer_tendency)

    implicit none

    integer, intent(in) :: i_diff
    real(wp), dimension(:,:), intent(in) :: f_ocn
    real(wp), dimension(:,:,:), intent(in) :: tracer
    real(wp), intent(in) :: diff_iso
    real(wp), dimension(:,:,:), intent(in) :: diff_dia
    real(wp), dimension(:,:,:), intent(in) :: drho_dx, drho_dy, drho_dz
    real(wp), dimension(:,:), intent(in) :: slope_crit

    real(wp), dimension(0:,0:,0:), intent(out) :: fdx
    real(wp), dimension(0:,0:,0:), intent(out) :: fdy
    real(wp), dimension(0:,0:,0:), intent(out) :: fdz
    real(wp), dimension(:,:,:), intent(out) :: slope2_w
    real(wp), dimension(:,:,:), intent(out) :: tracer_tendency

    integer :: i, j, k, i0, j0, i1, i2, ip1, nnp, knp, k1, k2
    real(wp) :: dts_dx, dts_dy, dts_dz, tv
    real(wp) :: slope2, slope_x, slope_y, slope_x_avg, slope_y_avg, drhodx, drhody, drhodz, tracer_ijk
    real(wp) :: taper
    real(wp) :: dx_dz, dy_dz


    do k=1,maxk
      do j=1,maxj
        do i=1,maxi

          tracer_ijk = tracer(i,j,k)

          ! flux to east
          if (mask_u(I,j,k).eq.1) then  
            ip1 = i+1
            if (ip1.eq.maxi+1) ip1=1
            dy_dz = dy*dz(k)
            fdx(I,j,k) = - rdx(j) * min(diffx_max(j),diff_iso) * (tracer(ip1,j,k)-tracer_ijk) *dy_dz*dt  ! volume flux
            if (i_diff.eq.1 .and. diff_iso.ne.diff_gm) then
              drhodx = drho_dx(I,j,k)
              tv = 0._wp
              do knp=0,1
                do nnp=0,1
                  i0 = i+nnp
                  if (i0.eq.0) i0=maxi
                  if (i0.eq.maxi+1) i0=1
                  k1 = k+knp-1
                  k2 = k+knp
                  if (k1.gt.0) then
                    if (mask_w(i0,j,k1).eq.1) then
                      dts_dz = (tracer(i0,j,k2)-tracer(i0,j,k1)) * rdza(k1)  ! K/m or psu/m            
                      slope_x = drhodx/min(-eps,drho_dz(i0,j,k1))
                      ! slope tapering following Gerdes 1991
                      taper = min(1._wp,(slope_crit(j,k)/slope_x)**2)
                    else
                      dts_dz = 0._wp
                      slope_x = 0._wp
                      taper = 1._wp
                    endif
                    ! sum up over 4 corners
                    tv = tv + taper * min(diffx_max(j),(diff_iso-diff_gm))*slope_x*dts_dz 
                  endif
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
            fdy(i,J,k) = - rdy * min(diffy_max(J),diff_iso) * (tracer(i,j+1,k)-tracer_ijk) *dx_dz*dt
            if (i_diff.eq.1 .and. diff_iso.ne.diff_gm) then
              ! add isoneutral diffusion
              drhody = drho_dy(i,J,k)
              tv = 0._wp
              do knp=0,1
                do nnp=0,1
                  j0 = j+nnp
                  k1 = k+knp-1
                  k2 = k+knp
                  if (k1.gt.0) then
                    if (mask_w(i,j0,k1).eq.1) then
                      dts_dz = (tracer(i,j0,k2)-tracer(i,j0,k1)) * rdza(k1)  ! K/m or psu/m            
                      slope_y = drhody/min(-eps,drho_dz(i,j0,k1))
                      ! slope tapering following Gerdes 1991
                      taper = min(1._wp,(slope_crit(j0,k)/slope_y)**2)
                    else
                      dts_dz = 0._wp
                      slope_y = 0._wp
                      taper = 1._wp
                    endif
                    ! sum up over 4 corners
                    tv = tv + taper * min(diffy_max(j),(diff_iso-diff_gm))*slope_y*dts_dz 
                  endif
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
            if (.not.l_diff33_impl) then
              fdz(i,j,K) = - diff_dia(i,j,K) * dts_dz *ocn_area(i,j)*dt
            else
              fdz(i,j,K) = 0._wp
            endif
            if (i_diff.eq.1) then
              ! add isoneutral diffusion
              drhodz = min(-eps,drho_dz(i,j,K)) ! negative
              ! compute horizontal derivatives of tracer
              tv = 0._wp
              slope_x_avg = 0._wp
              slope_y_avg = 0._wp
              do knp=0,1
                do nnp=0,1
                  ! phi derivative of tracer field
                  i0 = i-1+2*nnp
                  if (i0.eq.0) i0=maxi
                  if (i0.eq.maxi+1) i0=1
                  if (mask_c(i0,j,k+knp).eq.1) then
                    i1 = i+nnp-1
                    if (i1.eq.0) i1=maxi
                    i2 = i+nnp
                    if (i2.eq.maxi+1) i2=1
                    dts_dx = (tracer(i2,j,k+knp)-tracer(i1,j,k+knp)) * rdx(j)
                    slope_x = drho_dx(I1,j,k+knp)/drhodz
                  else
                    dts_dx = 0._wp
                    slope_x = 0._wp
                  endif
                  ! theta-derivative of tracer field
                  j0 = j-1+2*nnp
                  if (j0.eq.0 .or. j0.eq.maxj+1) then
                    dts_dy = 0._wp
                    slope_y = 0._wp
                  else if (mask_c(i,j0,k+knp).eq.0) then
                    dts_dy = 0._wp
                    slope_y = 0._wp
                  else
                    dts_dy = (tracer(i,j+nnp,k+knp)-tracer(i,j+nnp-1,k+knp)) * rdy
                    slope_y = drho_dy(i,J+nnp-1,k+knp)/drhodz
                  endif
                  ! sum up over 4 corners
                  slope2 = slope_x**2+slope_y**2+eps
                  taper = min(1._wp,slope_crit(j,K)**2/slope2)
                  tv = tv + (diff_iso+diff_gm) *taper* (slope_x*dts_dx + slope_y*dts_dy)
                  if (.not.l_diff33_impl) then
                    ! explicitely include A33 term 
                    tv = tv - diff_iso*taper*slope2*dts_dz
                  endif
                  ! average slopes
                  slope_x_avg = slope_x_avg + slope_x
                  slope_y_avg = slope_y_avg + slope_y
                enddo
              enddo
              slope2 = (slope_x_avg/4._wp)**2+(slope_y_avg/4._wp)**2+eps
              slope2_w(i,j,K) = slope2      ! save for diffusion_33
              ! tapering for large slopes, e.g. http://mitgcm.org/public/r2_manual/latest/online_documents/node240.html
              ! slope tapering following Gerdes 1991
              fdz(i,j,K) = fdz(i,j,K) + 0.25_wp*tv *ocn_area(i,j)*dt  ! m/s*K *m2*s = m3 * K
            endif
          else
            slope2_w(i,j,K) = eps
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
          if (mask_c(i,j,k).eq.1) then
            tracer_tendency(i,j,k) = - (fdx(I,j,k)-fdx(I-1,j,k) + fdy(i,J,k)-fdy(i,J-1,k) + fdz(i,j,K)-fdz(i,j,K-1)) &   ! m3 * K
                                   / (dx(j)*dy*dz(k)*f_ocn(i,j) * dt) 
          endif
        enddo
      enddo
    enddo

   return

  end subroutine diffusion


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  d i f f u s i o n _ 3 3
  !   Purpose    :  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine diffusion_33(i_diff,tracer,diff_iso,diff_dia,slope2_w,slope_max)

    implicit none

    integer, intent(in) :: i_diff
    real(wp), dimension(:,:,:), intent(inout) :: tracer
    real(wp), intent(in) :: diff_iso
    real(wp), dimension(:,:,:), intent(in) :: diff_dia
    real(wp), dimension(:,:,:), intent(in) :: slope2_w 
    real(wp), intent(in) :: slope_max

    integer :: i, j, k, kk, k1_ij
    real(wp) :: slope2, taper
    real(wp), dimension(maxk) :: a, b, c, r
    real(wp), dimension(maxk) :: diff_int


    do j=1,maxj
      do i=1,maxi
        if (mask_ocn(i,j).eq.1) then

          ! index of bottom layer
          k1_ij = k1(i,j)

          ! diffusivity profile on w-grid
          do k=k1_ij,maxk
            slope2 = slope2_w(i,j,K)
            ! tapering for large slopes, e.g. http://mitgcm.org/public/r2_manual/latest/online_documents/node240.html
            ! slope tapering following Gerdes 1991
            taper = min(1._wp,slope_max**2/slope2)
            if (i_diff.eq.1) then
              diff_int(K) = diff_dia(i,j,K) + diff_iso*taper*slope2
            else
              diff_int(K) = diff_dia(i,j,K) 
            endif
          enddo

          ! bottom layer
          kk = 1
          k = k1_ij
          a(kk) = 0._wp
          b(kk) = 1._wp + dt/dz(k) * diff_int(K)/dza(K)
          c(kk) = - dt/dz(k) * diff_int(K)/dza(K)
          r(kk) = tracer(i,j,k) 

          ! intermediate layers
          do k=k1_ij+1,maxk-1
            kk = kk+1
            a(kk) = - dt/dz(k) * diff_int(K-1)/dza(K-1)
            b(kk) = 1._wp + dt/dz(k) * ( diff_int(K-1)/dza(K-1) + diff_int(K)/dza(K) )
            c(kk) = - dt/dz(k) * diff_int(K)/dza(K)
            r(kk) = tracer(i,j,k)
          enddo

          ! top layer
          kk = kk+1
          k = maxk
          a(kk) = - dt/dz(k) * diff_int(K-1)/dza(K-1)
          b(kk) = 1._wp + dt/dz(k) * diff_int(K-1)/dza(K-1)
          c(kk) = 0._wp
          r(kk) = tracer(i,j,k) 

          ! solve tridiagonal system
          call tridiag_solve(a(1:kk),b(1:kk),c(1:kk),r(1:kk),tracer(i,j,k1_ij:maxk),kk)

        endif
      enddo
    enddo

    return

  end subroutine diffusion_33

end module diffusion_mod
