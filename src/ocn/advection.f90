!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : a d v e c t i o n _ m o d
!
!  Purpose : ocean tracer advection
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
module advection_mod

  use precision, only : wp
  use ocn_grid, only : maxi, maxj, maxk, k1, dx, dxv, dy, dz, dza, mask_ocn, mask_c, mask_u, mask_v, mask_w
  use ocn_params , only : dt, diff_iso, diff_dia
  use constants, only : r_earth

  implicit none

  private
  public :: advection_upstream, advection_fct

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  a d v e c t i o n _ u p s t r e a m
  !   Purpose    :  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine advection_upstream(f_ocn,u,tracer,flx_sur,flx_bot, &
                               fax,fay,faz,tracer_tendency)

    implicit none

    real(wp), dimension(:,:), intent(in) :: f_ocn
    real(wp), dimension(:,:,:,:), intent(in) :: u
    real(wp), dimension(:,:,:), intent(in) :: tracer
    real(wp), dimension(:,:), intent(in) :: flx_sur
    real(wp), dimension(:,:), intent(in) :: flx_bot

    real(wp), dimension(0:,0:,0:), intent(out) :: fax
    real(wp), dimension(0:,0:,0:), intent(out) :: fay
    real(wp), dimension(0:,0:,0:), intent(out) :: faz
    real(wp), dimension(:,:,:), intent(out) :: tracer_tendency

    integer :: i, j, k, ip1
    real(wp) :: ups(3), pec(3)
    

    do i=1,maxi
      ip1 = i+1
      if (ip1.eq.maxi+1) ip1=1
      do j=1,maxj
        !if (mask_ocn(i,j).eq.1) then
          do k=1,maxk
            pec(1) = u(1,i,j,k)*dx(j)/diff_iso
            ups(1) = pec(1) / (2.0 + abs(pec(1)))
            pec(2) = u(2,i,j,k)*dy/diff_iso
            ups(2) = pec(2) / (2.0 + abs(pec(2)))
            pec(3) = u(3,i,j,k)*dza(k)/diff_dia(i,j,k)
            ups(3) = pec(3) / (2.0 + abs(pec(3)))

            ! flux to east
            if (mask_u(i,j,k).eq.1) then
              fax(i,j,k) = u(1,i,j,k)*0.5_wp*((1.-ups(1))*tracer(ip1,j,k) + (1.+ups(1))*tracer(i,j,k)) * dy*dz(k)*dt  ! m/s*K * m2*s = m3 * K
            else
              fax(i,j,k) = 0
            endif
            ! flux to north
            if (j.lt.maxj .and. mask_v(i,j,k).eq.1) then
              fay(i,j,k) = u(2,i,j,k)*0.5_wp*((1.-ups(2))*tracer(i,j+1,k) + (1.+ups(2))*tracer(i,j,k)) * dxv(j)*dz(k)*dt  ! m/s*K * m2*s = m3 * K
            else
              fay(i,j,k) = 0._wp
            endif
            ! flux up
            if (mask_w(i,j,k).eq.1) then
              faz(i,j,k) = u(3,i,j,k)*0.5_wp*((1.-ups(3))*tracer(i,j,k+1) + (1.+ups(3))*tracer(i,j,k)) * dx(j)*dy*dt ! m/s*K * m2*s = m3 * K
            else if (k.eq.maxk) then
              faz(i,j,k) = flx_sur(i,j) * dx(j)*dy*f_ocn(i,j)*dt ! m/s*K * m2*s = m3 * K, atmosphere-ocean flux
            else if (k.eq.(k1(i,j)-1)) then
              faz(i,j,k) = flx_bot(i,j) * dx(j)*dy*f_ocn(i,j)*dt ! m/s*K * m2*s = m3 * K, bottom ocean flux
            else
              faz(i,j,k) = 0._wp
            endif
          enddo
        !endif
      enddo
    enddo

    ! western boundary fluxes, periodic boundary condition
    fax(0,:,:) = fax(maxi,:,:)

    ! southern boundary fluxes
    fay(:,0,:) = 0._wp

    ! bottom boundary fluxes
    faz(:,:,0) = 0._wp

    ! tracer tendency due to advection, K/s
    do i=1,maxi
      do j=1,maxj
        do k=1,maxk
          if (mask_c(i,j,k).eq.1) then
            tracer_tendency(i,j,k) = -(fax(i,j,k)-fax(i-1,j,k) + fay(i,j,k)-fay(i,j-1,k) + faz(i,j,k)-faz(i,j,k-1)) &
                            / (dx(j)*dy*dz(k)*f_ocn(i,j)*dt)
          endif
        enddo
      enddo
    enddo

   return

  end subroutine advection_upstream


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  a d v e c t i o n _ f c t
  !   Purpose    :  Flux-Correction Transport scheme (Zalesak 1979)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine advection_fct(f_ocn,u,tracer,flx_sur,flx_bot, &
                          fax,fay,faz,tracer_tendency)

    implicit none

    real(wp), dimension(:,:), intent(in) :: f_ocn
    real(wp), dimension(:,:,:,:), intent(in) :: u
    real(wp), dimension(:,:,:), intent(in) :: tracer
    real(wp), dimension(:,:), intent(in) :: flx_sur
    real(wp), dimension(:,:), intent(in) :: flx_bot
    real(wp), dimension(0:,0:,0:), intent(out) :: fax
    real(wp), dimension(0:,0:,0:), intent(out) :: fay
    real(wp), dimension(0:,0:,0:), intent(out) :: faz
    real(wp), dimension(:,:,:), intent(out) :: tracer_tendency

    real(wp), allocatable, dimension(:,:,:) :: tracer_tmp
    real(wp) :: tracer_ijk
    real(wp), dimension(2) :: tracer_i, tracer_j, tracer_k
    real(wp), dimension(:,:,:), allocatable :: fxl, fyl, fzl, afx, afy, afz
    real(wp), dimension(:,:,:), allocatable :: xa, xb, rp, rn
    real(wp) :: fxh, fyh, fzh
    integer :: i, j, k
    integer :: im1, ip1, ip2, jm1, jp1, jp2, km1, kp1, kp2
    real(wp) :: uvel, vvel, wvel
    real(wp) :: dy_dz, dx_dz
    real(wp) :: pp, pm, qp, qm, qdp
    real(wp) :: xmax, xmin
    real(wp) :: tracer_tmp1

    allocate(tracer_tmp(maxi,maxj,maxk))

    allocate(fxl(0:maxi,0:maxj,0:maxk), source=0._wp)
    allocate(fyl(0:maxi,0:maxj,0:maxk), source=0._wp)
    allocate(fzl(0:maxi,0:maxj,0:maxk), source=0._wp)
    allocate(afx(0:maxi,0:maxj,0:maxk), source=0._wp)
    allocate(afy(0:maxi,0:maxj,0:maxk), source=0._wp)
    allocate(afz(0:maxi,0:maxj,0:maxk), source=0._wp)

    allocate(xa(maxi,maxj,maxk))
    allocate(xb(maxi,maxj,maxk))
    allocate(rp(maxi,maxj,maxk))
    allocate(rn(maxi,maxj,maxk))


    ! volume flux calculation 
    ! for low and high oder solutions

    ! X-direction
    !!$omp parallel do collapse(3) private(i,j,k,uvel,ip1,tracer_i,fxh)
    do k=1,maxk
      do j=1,maxj
        do i=1,maxi
          if (mask_u(i,j,k).eq.1) then
            ip1 = modulo(i,maxi) + 1
            uvel = u(1,i,j,k)
            tracer_i = tracer((/i,ip1/),j,k)
            dy_dz = dy*dz(k)
            ! low order, ustream
            fxl(i,j,k) = uvel * merge(tracer_i(1), tracer_i(2), uvel > 0._wp) * dy_dz * dt  ! m/s*K * m2*s = m3 * K
            ! high order, centered differences
            fxh=0.5_wp*(tracer_i(1)+tracer_i(2))*uvel*dy_dz*dt
            afx(i,j,k)=fxh-fxl(i,j,k)
            !if (i.eq.maxi) then
            !  fxl(0,j,k) = fxl(maxi,j,k)
            !  afx(0,j,k) = afx(maxi,j,k)
            !endif
          else
            fxl(i,j,k)=0._wp
            afx(i,j,k)=0._wp
          endif
        enddo
      enddo
    enddo
    !!$omp end parallel do

    ! periodic boundary conditions
    fxl(0,:,:) = fxl(maxi,:,:)
    afx(0,:,:) = afx(maxi,:,:)

    ! Y-direction
    !!$omp parallel do collapse(3) private(i,j,k,vvel,tracer_j,fyh)
    do k=1,maxk
      do j=1,maxj-1
        do i=1,maxi
          if (mask_v(i,j,k).eq.1) then
            vvel = u(2,i,j,k)
            tracer_j = tracer(i,j:j+1,k)
            dx_dz = dxv(j)*dz(k)
            ! low order
            fyl(i,j,k) = vvel * merge(tracer_j(1), tracer_j(2), vvel > 0._wp) * dx_dz * dt  ! m/s*K * m2*s = m3 * K
            ! high order
            fyh=0.5_wp*(tracer_j(1)+tracer_j(2))*vvel*dx_dz*dt
            afy(i,j,k)=fyh-fyl(i,j,k)
          else
            fyl(i,j,k)=0._wp
            afy(i,j,k)=0._wp
          endif
        enddo
      enddo
    enddo
    !!$omp end parallel do

    ! no meridional flux across South Pole and North Pole
    fyl(:,0,:) = 0._wp
    afy(:,0,:) = 0._wp
    fyl(:,maxj,:) = 0._wp
    afy(:,maxj,:) = 0._wp

    ! Z-direction
    !!$omp parallel do collapse(3) private(i,j,k,wvel,tracer_k,fzh)
    do k=1,maxk
      do j=1,maxj
        do i=1,maxi
          if (mask_w(i,j,k).eq.1) then
            wvel = u(3,i,j,k)
            tracer_k = tracer(i,j,k:k+1)
            ! low order
            fzl(i,j,k) = wvel * merge(tracer_k(1), tracer_k(2), wvel > 0._wp) * dx(j)*dy * dt  ! m/s*K * m2*s = m3 * K
            ! high order
            fzh=0.5_wp*(tracer_k(1)+tracer_k(2))*wvel*dx(j)*dy*dt
            afz(i,j,k)=fzh-fzl(i,j,k)
          else if (k.eq.maxk .and. mask_ocn(i,j).eq.1) then
            ! atmosphere-ocean flux
            fzl(i,j,k)=flx_sur(i,j)*dx(j)*dy*f_ocn(i,j)*dt  ! m/s * K * m2*s = m3*K
            afz(i,j,k)=0._wp
          else if (k.eq.(k1(i,j)-1)) then
            ! bottom ocean flux
            fzl(i,j,k)=flx_bot(i,j)*dx(j)*dy*f_ocn(i,j)*dt  ! m/s * K * m2*s = m3*K
            afz(i,j,k)=0._wp
          else
            fzl(i,j,k)=0._wp
            afz(i,j,k)=0._wp
          endif
        enddo
      enddo
    enddo
    !!$omp end parallel do
    fzl(:,:,0) = 0._wp
    afz(:,:,0) = 0._wp

    ! STEP I: LOWER ORDER SOLUTION
    !!$omp parallel do collapse(3) private(i,j,k)
    do k=1,maxk
      do j=1,maxj
        do i=1,maxi
          tracer_ijk = tracer(i,j,k)
          if (mask_c(i,j,k).eq.1) then
            ! low order solution
            tracer_tmp(i,j,k) = tracer_ijk-(fxl(i,j,k)-fxl(i-1,j,k)+fyl(i,j,k)-fyl(i,j-1,k)+fzl(i,j,k)-fzl(i,j,k-1)) &
                              / (dx(j)*dy*dz(k)*f_ocn(i,j))
          else
            tracer_tmp(i,j,k) = 0._wp
          endif 
          ! for fluxes limits 
          xa(i,j,k)=max(tracer_ijk,tracer_tmp(i,j,k))
          xb(i,j,k)=min(tracer_ijk,tracer_tmp(i,j,k))
        enddo
      enddo      
    enddo
    !!$omp end parallel do

    ! flux limiters

    ! apply eq. 14' of Zalesak 1979
    do k=1,maxk
      do j=1,maxj
        do i=1,maxi
          ip1 = modulo(i,maxi) + 1
          im1 = modulo(i - 2, maxi) + 1
          ip2 = modulo(i+1, maxi) + 1
          jm1=max(1,j-1)
          jp1=min(maxj,j+1)
          jp2=min(maxj,j+2)
          km1=max(1,k-1)
          kp1=min(maxk,k+1)
          kp2=min(maxk,k+2)
          if (afx(I,j,k)*(tracer_tmp(ip1,j,k)-tracer_tmp(i,j,k)).lt.0._wp &
            .and. (afx(I,j,k)*(tracer_tmp(ip2,j,k)-tracer_tmp(ip1,j,k)).lt.0._wp .or. afx(I,j,k)*(tracer_tmp(i,j,k)-tracer_tmp(im1,j,k)).lt.0._wp)) then
            !print *,'afx set to zero'
            afx(I,j,k) = 0._wp
          endif
          if (afy(i,J,k)*(tracer_tmp(i,jp1,k)-tracer_tmp(i,j,k)).lt.0._wp &
            .and. (afy(i,J,k)*(tracer_tmp(i,jp2,k)-tracer_tmp(i,jp1,k)).lt.0._wp .or. afy(i,J,k)*(tracer_tmp(i,j,k)-tracer_tmp(i,jm1,k)).lt.0._wp)) then
            !print *,'afy set to zero'
            afy(i,J,k) = 0._wp
          endif
          if (afz(i,j,K)*(tracer_tmp(i,j,kp1)-tracer_tmp(i,j,k)).lt.0._wp &
            .and. (afz(i,j,K)*(tracer_tmp(i,j,kp2)-tracer_tmp(i,j,kp1)).lt.0._wp .or. afz(i,j,K)*(tracer_tmp(i,j,k)-tracer_tmp(i,j,km1)).lt.0._wp)) then
            !print *,'afz set to zero'
            afz(i,j,K) = 0._wp
          endif
        enddo
      enddo
    enddo

    !!$omp parallel do collapse(3) private(i,j,k,im1,ip1,jm1,jp1,km1,kp1,xmin,xmax,pp,pm,qp,qm,qdp)
    do k=1,maxk    
      do j=1,maxj
        do i=1,maxi
          im1 = modulo(i-2, maxi) + 1
          ip1 = modulo(i,maxi) + 1
          jm1=max(1,j-1)
          jp1=min(maxj,j+1)
          km1=max(1,k-1)
          kp1=min(maxk,k+1)
          xmax = max(xa(im1,j,k),xa(ip1,j,k),xa(i,jm1,k),xa(i,jp1,k),xa(i,j,km1),xa(i,j,kp1),xa(i,j,k))
          xmin = min(xb(im1,j,k),xb(ip1,j,k),xb(i,jm1,k),xb(i,jp1,k),xb(i,j,km1),xb(i,j,kp1),xb(i,j,k))
          pp = max(0._wp,afx(im1,j,k))-min(0._wp,afx(i,j,k))+max(0._wp,afy(i,jm1,k))-min(0._wp,afy(i,j,k))+max(0._wp,afz(i,j,km1))-min(0._wp,afz(i,j,k))
          pm = max(0._wp,afx(i,j,k))-min(0._wp,afx(im1,j,k))+max(0._wp,afy(i,j,k))-min(0._wp,afy(i,jm1,k))+max(0._wp,afz(i,j,k))-min(0._wp,afz(i,j,km1))
          qp = (xmax-tracer_tmp(i,j,k))*(dx(j)*dy*dz(k)*f_ocn(i,j))
          qm = (tracer_tmp(i,j,k)-xmin)*(dx(j)*dy*dz(k)*f_ocn(i,j))
          qdp = qp / max(pp, 1.e-3_wp)
          rp(i,j,k) = merge(min(1._wp, qdp), 0._wp, pp >= 1.e-3_wp)
          qdp = qm / max(pm, 1.e-3_wp)
          rn(i,j,k) = merge(min(1._wp, qdp), 0._wp, pm >= 1.e-3_wp)
        enddo
      enddo
    enddo
    !!$omp end parallel do

    !!$omp parallel do collapse(3) private(i,j,k,ip1)
    do k=1,maxk
      do j=1,maxj
        do i=1,maxi
          if (mask_u(i,j,k).eq.1) then
            ip1 = modulo(i,maxi) + 1
            afx(i,j,k) = afx(i,j,k) * merge(min(rp(ip1,j,k),rn(i,j,k)), min(rp(i,j,k),rn(ip1,j,k)), afx(i,j,k) >= 0._wp)
          endif
          if (j.lt.maxj .and. mask_v(i,j,k).eq.1) then
            afy(i,j,k) = afy(i,j,k) * merge(min(rp(i,j+1,k),rn(i,j,k)), min(rp(i,j,k),rn(i,j+1,k)), afy(i,j,k) >= 0._wp)
          endif
          if (mask_w(i,j,k).eq.1) then
            afz(i,j,k) = afz(i,j,k) * merge(min(rp(i,j,k+1),rn(i,j,k)), min(rp(i,j,k),rn(i,j,k+1)), afz(i,j,k) >= 0._wp)
          endif
        enddo
      enddo
    enddo
    !!$omp end parallel do

    !periodic b.c.
    afx(0,:,:) = afx(maxi,:,:)


    ! STEP II: HIGH ORDER CORRECTION   

    !!$omp parallel do collapse(3) private(i,j,k,tracer_tmp1)
    do k=1,maxk
       do j=1,maxj
          do i=1,maxi
             if (mask_c(i,j,k).eq.1) then
                tracer_tmp1 = tracer_tmp(i,j,k)-(afx(i,j,k)-afx(i-1,j,k)+afy(i,j,k)-afy(i,j-1,k)+afz(i,j,k)-afz(i,j,k-1)) &
                            / (dx(j)*dy*dz(k)*f_ocn(i,j))
                ! calculate the overall tracer tendency due to advection
                tracer_tendency(i,j,k) = (tracer_tmp1 - tracer(i,j,k)) / dt
             endif
          enddo
       enddo 
    enddo
    !!$omp end parallel do

    fax(0:,:,:) = fxl(0:,:,:) + afx(0:,:,:)
    fay(:,0:,:) = fyl(:,0:,:) + afy(:,0:,:)
    faz(:,:,0:) = fzl(:,:,0:) + afz(:,:,0:)

    deallocate(tracer_tmp)
    deallocate(fxl,fyl,fzl,afx,afy,afz)
    deallocate(xa,xb,rp,rn)

   return

  end subroutine advection_fct


end module advection_mod
