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
  use ocn_params, only : f_face_dif, drag_par
  use ocn_params, only : dt, diff_gm, diffx_max, diffy_max, i_taper, slope_width, slope_max
  use ocn_grid, only : k1, maxi, maxj, maxk, mask_c, mask_u, mask_v, mask_w, dx, dxv, dy, dz, dza, rdx, rdy, rdza, zw, ocn_area

  implicit none

  real(wp), parameter :: eps = 1.e-15_wp

  ! tracer-independent coefficients of the diffusive fluxes, computed once per time step by diffusion_coef
  ! and applied to each tracer by diffusion; fluxes are volume fluxes times tracer [m3 * tracer unit]
  real(wp), allocatable, dimension(:,:,:) :: cfx, cfy, cfz   ! diagonal terms: horizontal (tapered A11, A22) and diapycnal [m3]
  real(wp), allocatable, dimension(:,:,:,:) :: cex, cey      ! horizontal faces, 4 corners: skew/Redi cross terms (diff_iso/=diff_gm) [m3*m]
  real(wp), allocatable, dimension(:,:,:,:) :: czx, czy      ! vertical faces, 4 corners: isoneutral + GM cross terms (A31, A32) [m3]
  real(wp), allocatable, dimension(:,:,:) :: ta33, cp33, rb33   ! implicit A33: sub-diagonal, Thomas c' and 1/pivot

  private
  public :: diffusion_coef, diffusion, slope_taper

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  d i f f u s i o n _ c o e f
  !   Purpose    :  tracer-independent coefficients of the isoneutral, GM, horizontal and
  !                 diapycnal diffusion (slopes, tapers, boundary layer, implicit A33 matrix)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine diffusion_coef(diff_iso,diff_dia,drho_dx,drho_dy,drho_dz,slope_crit,hbl,sx_bl,sy_bl)

    implicit none

    real(wp), intent(in) :: diff_iso
    real(wp), dimension(:,:,:), intent(in) :: diff_dia
    real(wp), dimension(:,:,:), intent(in) :: drho_dx, drho_dy, drho_dz
    real(wp), dimension(:,:), intent(in) :: slope_crit
    real(wp), dimension(:,:), intent(in) :: hbl   ! surface boundary layer depth (positive) [m]
    real(wp), dimension(:,:), intent(in) :: sx_bl, sy_bl   ! tapered slopes at the base of the boundary layer []

    integer :: i, j, k, i0, j0, i1, nnp, knp, kc, kk, n
    integer :: jm, j0c
    real(wp) :: slope2, slope_x, slope_y, drhodz
    real(wp) :: taper, mc_x, yfac
    real(wp) :: tap_sum, tap_n   ! face-mean slope taper of the horizontal diagonal terms (A11, A22)
    real(wp) :: a33              ! A33 coefficient on the w-face, K*S^2 averaged over the 4 corners [m2/s]
    real(wp) :: dkx, dky
    real(wp) :: tb, tc, rm
    logical :: l_bl        ! w-point inside the boundary layer, linear slope profile
    real(wp) :: f_bl       ! depth/hbl of the w-point inside the boundary layer
    real(wp), dimension(maxk) :: k33


    if (.not.allocated(cfx)) then
      allocate(cfx(maxi,maxj,maxk), cfy(maxi,maxj,maxk), cfz(maxi,maxj,maxk), source=0._wp)
      allocate(cex(4,maxi,maxj,maxk), cey(4,maxi,maxj,maxk), source=0._wp)
      allocate(czx(4,maxi,maxj,maxk), czy(4,maxi,maxj,maxk), source=0._wp)
      allocate(ta33(maxi,maxj,maxk), cp33(maxi,maxj,maxk), rb33(maxi,maxj,maxk), source=0._wp)
    endif

    !$omp parallel do collapse(2) private(i,j,k,i0,j0,i1,nnp,knp,kc,jm,j0c,slope2,slope_x,slope_y,drhodz, &
    !$omp taper,mc_x,yfac,tap_sum,tap_n,a33,dkx,dky,l_bl,f_bl)
    do k=1,maxk
      do j=1,maxj
        do i=1,maxi

          ! flux to east
          cfx(i,j,k) = 0._wp
          cex(:,i,j,k) = 0._wp
          if (mask_u(i,j,k).eq.1) then
            tap_sum = 0._wp
            tap_n = 0._wp
            do knp=0,1
              do nnp=0,1
                i0 = modulo(i+nnp-1, maxi) + 1
                kc = max(1, k+knp-1)
                if (mask_w(i0,j,kc).eq.1) then
                  if (zw(kc).ge.-hbl(i0,j)) then
                    tap_sum = tap_sum + 1._wp   ! boundary layer: no taper
                  else
                    slope_x = drho_dx(i,j,k)/min(-eps,drho_dz(i0,j,kc))
                    tap_sum = tap_sum + slope_taper(slope_x**2,slope_crit(j,k))
                  endif
                  tap_n = tap_n + 1._wp
                endif
              enddo
            enddo
            taper = 1._wp
            if (tap_n.gt.0._wp) taper = tap_sum/tap_n
            cfx(i,j,k) = rdx(j) * taper * min(diffx_max(j),diff_iso) * dy*dz(k)*dt * f_face_dif(1,i,j)
            if (diff_iso.ne.diff_gm) then
              dkx = min(diffx_max(j),diff_iso)-min(diffx_max(j),diff_gm)
              do knp=0,1
                do nnp=0,1
                  i0 = modulo(i+nnp-1, maxi) + 1   ! wraps i+nnp = maxi+1 back to 1
                  kc = max(1, k+knp-1)             ! when k+knp-1=0 the corner has dts_dz=0 (vanishes)
                  if (zw(kc).ge.-hbl(i0,j)) then
                    ! boundary layer: slope linear in depth from its value at the base, no taper
                    slope_x = -zw(kc)/hbl(i0,j)*sx_bl(i0,j)
                    taper = 1._wp
                  else
                    slope_x = drho_dx(i,j,k)/min(-eps,drho_dz(i0,j,kc))
                    taper = slope_taper(slope_x**2,slope_crit(j,k))
                  endif
                  ! mask_w as a multiply (0 on closed face)
                  cex(1+nnp+2*knp,i,j,k) = 0.25_wp*real(mask_w(i0,j,kc),wp) * taper*dkx*slope_x*rdza(kc) *dy*dz(k)*dt * f_face_dif(1,i,j)
                enddo
              enddo
            endif
          endif

          ! flux to north
          cfy(i,j,k) = 0._wp
          cey(:,i,j,k) = 0._wp
          if (j.lt.maxj .and. mask_v(i,j,k).eq.1) then
            tap_sum = 0._wp
            tap_n = 0._wp
            do knp=0,1
              do nnp=0,1
                j0 = j+nnp
                kc = max(1, k+knp-1)
                if (mask_w(i,j0,kc).eq.1) then
                  if (zw(kc).ge.-hbl(i,j0)) then
                    tap_sum = tap_sum + 1._wp   ! boundary layer: no taper
                  else
                    slope_y = drho_dy(i,j,k)/min(-eps,drho_dz(i,j0,kc))
                    tap_sum = tap_sum + slope_taper(slope_y**2,slope_crit(j0,k))
                  endif
                  tap_n = tap_n + 1._wp
                endif
              enddo
            enddo
            taper = 1._wp
            if (tap_n.gt.0._wp) taper = tap_sum/tap_n
            cfy(i,j,k) = rdy * taper * min(diffy_max(j),diff_iso) * dxv(j)*dz(k)*dt * f_face_dif(2,i,j)
            if (diff_iso.ne.diff_gm) then
              dky = min(diffy_max(j),diff_iso)-min(diffy_max(j),diff_gm)
              do knp=0,1
                do nnp=0,1
                  j0 = j+nnp
                  kc = max(1, k+knp-1)             ! when k+knp-1=0 the corner has dts_dz=0 (vanishes)
                  if (zw(kc).ge.-hbl(i,j0)) then
                    ! boundary layer: slope linear in depth from its value at the base, no taper
                    slope_y = -zw(kc)/hbl(i,j0)*sy_bl(i,j0)
                    taper = 1._wp
                  else
                    slope_y = drho_dy(i,j,k)/min(-eps,drho_dz(i,j0,kc))
                    taper = slope_taper(slope_y**2,slope_crit(j0,k))
                  endif
                  ! mask_w as a multiply (0 on closed face)
                  cey(1+nnp+2*knp,i,j,k) = 0.25_wp*real(mask_w(i,j0,kc),wp) * taper*dky*slope_y*rdza(kc) *dxv(j)*dz(k)*dt * f_face_dif(2,i,j)
                enddo
              enddo
            endif
          endif

          ! flux up
          cfz(i,j,k) = 0._wp
          czx(:,i,j,k) = 0._wp
          czy(:,i,j,k) = 0._wp
          ta33(i,j,k) = 0._wp   ! used as temporary storage of K33 on the w-face
          if (mask_w(i,j,k).eq.1) then
            cfz(i,j,k) = diff_dia(i,j,k)*rdza(k) *ocn_area(i,j)*dt
            drhodz = min(-eps,drho_dz(i,j,k)) ! negative
            l_bl = zw(k).ge.-hbl(i,j)
            f_bl = merge(-zw(k)/hbl(i,j), 0._wp, l_bl)   ! linear slope profile inside the boundary layer
            a33 = 0._wp
            do knp=0,1
              do nnp=0,1
                ! phi corner (mask_c as a multiply, i1 wrapped via modulo)
                i0 = modulo(i-2+2*nnp, maxi) + 1   ! phi mask point, wraps 0 -> maxi and maxi+1 -> 1
                i1 = modulo(i+nnp-2, maxi) + 1     ! = i+nnp-1, wrapped 0 -> maxi
                mc_x = real(mask_c(i0,j,k+knp),wp)
                if (l_bl) then
                  slope_x = mc_x * f_bl*sx_bl(i,j)
                else
                  slope_x = mc_x * drho_dx(i1,j,k+knp)/drhodz
                endif
                ! theta corner (clamp j indices, zero the contribution outside [1,maxj])
                j0 = j-1+2*nnp
                j0c = min(max(j0,1),maxj)          ! clamped theta mask index (only used when in-bounds)
                jm  = max(j+nnp-1, 1)
                yfac = merge(real(mask_c(i,j0c,k+knp),wp), 0._wp, j0.ge.1 .and. j0.le.maxj)
                if (l_bl) then
                  slope_y = yfac * f_bl*sy_bl(i,j)
                else
                  slope_y = yfac * drho_dy(i,jm,k+knp)/drhodz
                endif
                slope2 = slope_x**2+slope_y**2+eps
                if (l_bl) then
                  taper = 1._wp
                else
                  taper = slope_taper(slope2,slope_crit(j,k))
                endif
                ! cross terms, applied to the corner tracer differences in diffusion (mask factors folded in)
                ! with i_frac_dyn=2 the area of the vertical cross flux is set by the face it pairs with, so that the
                ! tensor keeps the symmetry between the lateral and the vertical cross terms
                if (drag_par%i_frac_dyn.eq.2 .or. drag_par%i_frac_dyn.ge.4) then
                  czx(1+nnp+2*knp,i,j,k) = 0.25_wp*(diff_iso+diff_gm)*taper*slope_x*mc_x*rdx(j) *dx(j)*dy*f_face_dif(1,i1,j)*dt
                  czy(1+nnp+2*knp,i,j,k) = 0.25_wp*(diff_iso+diff_gm)*taper*slope_y*yfac*rdy *dx(j)*dy*f_face_dif(2,i,jm)*dt
                else
                  czx(1+nnp+2*knp,i,j,k) = 0.25_wp*(diff_iso+diff_gm)*taper*slope_x*mc_x*rdx(j) *ocn_area(i,j)*dt
                  czy(1+nnp+2*knp,i,j,k) = 0.25_wp*(diff_iso+diff_gm)*taper*slope_y*yfac*rdy *ocn_area(i,j)*dt
                endif
                ! A33 term, solved implicitly
                a33 = a33 + diff_iso*taper*slope2
              enddo
            enddo
            ta33(i,j,k) = 0.25_wp*a33
          endif

        enddo
      enddo
    enddo
    !$omp end parallel do

    ! implicit A33 (backward Euler, per column): dz dT/dt = d/dz(K33 dT/dz), no flux at surface and bottom.
    ! LU (Thomas) factorization of the tridiagonal matrix, shared by all tracers
    !$omp parallel do collapse(2) private(i,j,k,kk,n,k33,tb,tc,rm)
    do j=1,maxj
      do i=1,maxi
        k33(:) = ta33(i,j,:)
        ta33(i,j,:) = 0._wp
        cp33(i,j,:) = 0._wp
        rb33(i,j,:) = 0._wp
        if (mask_c(i,j,maxk).eq.1) then
          n = maxk-k1(i,j)+1
          if (n.ge.2) then
            do kk=k1(i,j),maxk
              ta33(i,j,kk) = 0._wp
              tc = 0._wp
              if (kk.gt.k1(i,j)) ta33(i,j,kk) = -dt*k33(kk-1)/(dz(kk)*dza(kk-1))   ! face below cell kk
              if (kk.lt.maxk) tc = -dt*k33(kk)/(dz(kk)*dza(kk))                      ! face above cell kk
              tb = 1._wp - ta33(i,j,kk) - tc
              if (kk.eq.k1(i,j)) then
                rm = 1._wp/tb
              else
                rm = 1._wp/(tb-cp33(i,j,kk-1)*ta33(i,j,kk))
              endif
              rb33(i,j,kk) = rm
              cp33(i,j,kk) = tc*rm
            enddo
          endif
        endif
      enddo
    enddo
    !$omp end parallel do

   return

  end subroutine diffusion_coef


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  d i f f u s i o n
  !   Purpose    :  diffusive fluxes and tendency of one tracer, with the coefficients of diffusion_coef
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine diffusion(f_ocn,tracer,diff_iso,fdx, fdy, fdz, tracer_tendency)

    implicit none

    real(wp), dimension(:,:), intent(in) :: f_ocn
    real(wp), dimension(:,:,:), intent(in) :: tracer
    real(wp), intent(in) :: diff_iso

    real(wp), dimension(0:,0:,0:), intent(out) :: fdx
    real(wp), dimension(0:,0:,0:), intent(out) :: fdy
    real(wp), dimension(0:,0:,0:), intent(out) :: fdz
    real(wp), dimension(:,:,:), intent(out) :: tracer_tendency

    integer :: i, j, k, i0, j0, i1, i2, ip1, nnp, knp, kc, kk, n
    integer :: jp, jm
    real(wp) :: tv, tracer_ijk
    real(wp), dimension(maxk) :: dp, tx
    real(wp), parameter :: f_min = 1.e-30_wp   ! floor for f_ocn in branchless denominator (land f_ocn=0)


    !$omp parallel do collapse(2) private(i,j,k,i0,j0,i1,i2,ip1,nnp,knp,kc,jp,jm,tv,tracer_ijk)
    do k=1,maxk
      do j=1,maxj
        do i=1,maxi

          tracer_ijk = tracer(i,j,k)

          ! flux to east
          ip1 = modulo(i,maxi) + 1
          tv = - cfx(i,j,k) * (tracer(ip1,j,k)-tracer_ijk)
          if (diff_iso.ne.diff_gm) then
            do knp=0,1
              do nnp=0,1
                i0 = modulo(i+nnp-1, maxi) + 1
                kc = max(1, k+knp-1)
                tv = tv + cex(1+nnp+2*knp,i,j,k)*(tracer(i0,j,min(k+knp,maxk))-tracer(i0,j,kc))
              enddo
            enddo
          endif
          fdx(i,j,k) = tv

          ! flux to north
          if (j.lt.maxj) then
            tv = - cfy(i,j,k) * (tracer(i,j+1,k)-tracer_ijk)
            if (diff_iso.ne.diff_gm) then
              do knp=0,1
                do nnp=0,1
                  j0 = j+nnp
                  kc = max(1, k+knp-1)
                  tv = tv + cey(1+nnp+2*knp,i,j,k)*(tracer(i,j0,min(k+knp,maxk))-tracer(i,j0,kc))
                enddo
              enddo
            endif
            fdy(i,j,k) = tv
          else
            fdy(i,j,k) = 0._wp
          endif

          ! flux up
          if (mask_w(i,j,k).eq.1) then
            tv = - cfz(i,j,k) * (tracer(i,j,k+1) - tracer_ijk)
            do knp=0,1
              do nnp=0,1
                i1 = modulo(i+nnp-2, maxi) + 1     ! = i+nnp-1, wrapped 0 -> maxi
                i2 = modulo(i+nnp-1, maxi) + 1     ! = i+nnp,   wrapped maxi+1 -> 1
                jp = min(j+nnp, maxj)
                jm = max(j+nnp-1, 1)
                tv = tv + czx(1+nnp+2*knp,i,j,k)*(tracer(i2,j,k+knp)-tracer(i1,j,k+knp)) &
                        + czy(1+nnp+2*knp,i,j,k)*(tracer(i,jp,k+knp)-tracer(i,jm,k+knp))
              enddo
            enddo
            fdz(i,j,k) = tv
          else
            fdz(i,j,k) = 0._wp
          endif

        enddo
      enddo
    enddo
    !$omp end parallel do

    ! western boundary fluxes, periodic boundary condition
    fdx(0,:,:) = fdx(maxi,:,:)

    ! southern boundary fluxes
    fdy(:,0,:) = 0._wp

    ! bottom boundary fluxes
    fdz(:,:,0) = 0._wp

    !$omp parallel do collapse(2) private(i,j,k,kk,n,dp,tx)
    do j=1,maxj
      do i=1,maxi
        do k=1,maxk
          ! branchless mask_c; guard f_ocn (=0 on land) in the denominator
          tracer_tendency(i,j,k) = - real(mask_c(i,j,k),wp) &
                                 * (fdx(i,j,k)-fdx(i-1,j,k) + fdy(i,j,k)-fdy(i,j-1,k) + fdz(i,j,k)-fdz(i,j,k-1)) &   ! m3 * K
                                 / (dx(j)*dy*dz(k)*max(f_ocn(i,j),f_min) * dt)
        enddo
        ! implicit A33 with the factorization of diffusion_coef: forward and back substitution;
        ! the change is added to the tendency, so it is booked as diffusion like the explicit terms
        if (mask_c(i,j,maxk).eq.1) then
          n = maxk-k1(i,j)+1
          if (n.ge.2) then
            kk = k1(i,j)
            dp(kk) = tracer(i,j,kk)*rb33(i,j,kk)
            do kk=k1(i,j)+1,maxk
              dp(kk) = (tracer(i,j,kk)-dp(kk-1)*ta33(i,j,kk))*rb33(i,j,kk)
            enddo
            tx(maxk) = dp(maxk)
            do kk=maxk-1,k1(i,j),-1
              tx(kk) = dp(kk)-cp33(i,j,kk)*tx(kk+1)
            enddo
            do kk=k1(i,j),maxk
              tracer_tendency(i,j,kk) = tracer_tendency(i,j,kk) + (tx(kk)-tracer(i,j,kk))/dt
            enddo
          endif
        endif
      enddo
    enddo
    !$omp end parallel do

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
