!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : t r a n s p o r t _ s i c _ m o d
!
!  Purpose : sea ice transport
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
module transport_sic_mod

   use precision, only : wp
   use sic_grid, only : dx, dxv, dy, rdx, rdy, maxi, maxj
   use sic_params, only : dt, i_adv
   use sic_params, only : diffsic

   implicit none

   private
   public :: transport

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  t r a n s p o r t
  !   Purpose    :  sea ice and snow on sea ice transport (advection+diffusion)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine transport(u, v, f_ocn, &
                      his, dh_dt_dyn, fax, fay, fdx, fdy)

    implicit none

    real(wp), dimension(:,:), intent(in) :: u   !! zonal sea ice drift velocity [m/s]
    real(wp), dimension(:,:), intent(in) :: v   !! meridional sea ice drift velocity [m/s]
    real(wp), dimension(:,:), intent(in) :: f_ocn   !! ocean fraction []

    real(wp), dimension(:,:), intent(inout) :: his  !! sea ice or snow thickness [m]
    real(wp), dimension(:,:), intent(out), optional :: dh_dt_dyn    !! sea ice or snow thickness change due to transport [m]
    real(wp), dimension(0:,:), intent(out), optional :: fax    !! zonal volume flux from advection [m3]
    real(wp), dimension(:,0:), intent(out), optional :: fay    !! meridional volume flux from advection [m3]
    real(wp), dimension(0:,:), intent(out), optional :: fdx    !! zonal volume flux from diffusion [m3]
    real(wp), dimension(:,0:), intent(out), optional :: fdy    !! meridional volume flux from diffusion [m3]

    logical :: is_sic_snow
    real(wp), dimension(maxi,maxj) :: his_old, his_tmp
    real(wp), dimension(0:maxi,maxj) :: fxl, afx
    real(wp), dimension(maxi,0:maxj) :: fyl, afy
    real(wp), dimension(0:maxi,maxj) :: fdx1
    real(wp), dimension(maxi,0:maxj) :: fdy1
    real(wp), dimension(maxi,maxj) :: adv_tendency, diff_tendency
    real(wp) :: fxh, fyh
    real(wp), dimension(maxi,maxj) :: xa, xb, rp, rn
    integer :: i, j
    integer :: im1, ip1, ip2, jm1, jp1, jp2
    real(wp) :: uvel, vvel
    real(wp) :: pp, pm, qp, qm, qdp
    real(wp) :: xmax, xmin
    real(wp) :: cx, cy
    real(wp) :: diffx_max
    real(wp) :: diffy_max


    if (present(dh_dt_dyn)) then
      is_sic_snow = .true.
    else
      is_sic_snow = .false.
    endif

    his_old = his

    adv_tendency = 0._wp
    diff_tendency = 0._wp


    if (i_adv.eq.1) then

      !--------------------------------------------------
      ! advection (upstream)
      !--------------------------------------------------

      ! volume flux calculation 

      ! X-direction
      do j=1,maxj
        do i=1,maxi
          if (i.lt.maxi) then
            ip1 = i+1
          else if (i.eq.maxi) then
            ip1 = 1
          endif
          if (f_ocn(i,j).eq.0._wp .or. f_ocn(ip1,j).eq.0._wp) then
            fxl(I,j)=0._wp
          else
            uvel = u(I,j)
            ! low order
            if (uvel.gt.0._wp) then
              fxl(I,j)=uvel*his(i,j)*dy*dt  ! m/s*m*m*s = m3
            else
              fxl(I,j)=uvel*his(ip1,j)*dy*dt
            endif
          endif
        enddo
      enddo
      ! periodic boundary conditions
      fxl(0,:) = fxl(maxI,:)

      ! Y-direction
      do j=1,maxj-1
        do i=1,maxi
          if (f_ocn(i,j).eq.0._wp .or. f_ocn(i,j+1).eq.0._wp) then
            fyl(i,J)=0._wp
          else
            vvel = v(i,J)
            ! low order
            if (vvel.gt.0._wp) then
              fyl(i,J)=vvel*his(i,j)*dxv(J)*dt  ! m/s*m*m*s = m3 
            else
              fyl(i,J)=vvel*his(i,j+1)*dxv(J)*dt  
            endif
          endif
        enddo
      enddo
      ! no meridional flux across South Pole and North Pole
      fyl(:,maxJ) = 0._wp
      fyl(:,0) = 0._wp

      do j=1,maxj
        do i=1,maxi
          if (f_ocn(i,j).gt.0._wp) then
            adv_tendency(i,j) = - (fxl(I,j)-fxl(I-1,j) + fyl(i,J)-fyl(i,J-1)) &
              / (dx(j)*dy*f_ocn(i,j)) & ! area of ocean in grid cell
              / dt
          endif
        enddo
      enddo

      if (is_sic_snow) then
        fax(0:,:) = fxl(0:,:)
        fay(:,0:) = fyl(:,0:)
      endif

    else if (i_adv.eq.2) then

      !--------------------------------------------------
      ! advection (FCT)
      !--------------------------------------------------

      ! volume flux calculation 
      ! for lower and high oder solutions

      ! X-direction
      do j=1,maxj
        do i=1,maxi
          if (i.lt.maxi) then
            ip1 = i+1
          else if (i.eq.maxi) then
            ip1 = 1
          endif
          if (f_ocn(i,j).eq.0._wp .or. f_ocn(ip1,j).eq.0._wp) then
            fxl(I,j)=0._wp
            afx(I,j)=0._wp
          else
            uvel = u(I,j)
            ! low order
            if (uvel.gt.0._wp) then
              fxl(I,j)=uvel*his(i,j)*dy*dt  ! m/s*m*m*s = m3
            else
              fxl(I,j)=uvel*his(ip1,j)*dy*dt
            endif
            ! high order
            fxh=0.5_wp*(his(i,j)+his(ip1,j))*uvel*dy*dt  ! m/s*m*m*s = m3
            ! high-low order
            afx(I,j)=fxh-fxl(I,j)
          endif
        enddo
      enddo
      ! periodic boundary conditions
      fxl(0,:) = fxl(maxI,:)
      afx(0,:) = afx(maxI,:)

      ! Y-direction
      do j=1,maxj-1
        do i=1,maxi
          if (f_ocn(i,j).eq.0._wp .or. f_ocn(i,j+1).eq.0._wp) then
            fyl(i,J)=0._wp
            afy(i,J)=0._wp
          else
            vvel = v(i,J)
            ! low order
            if (vvel.gt.0._wp) then
              fyl(i,J)=vvel*his(i,j)*dxv(J)*dt  ! m/s*m*m*s = m3 
            else
              fyl(i,J)=vvel*his(i,j+1)*dxv(J)*dt  
            endif
            ! high order
            fyh=0.5_wp*(his(i,j)+his(i,j+1))*vvel*dxv(J)*dt  ! m/s*m*m*s = m3
            ! high-low order
            afy(i,J)=fyh-fyl(i,J)
          endif
        enddo
      enddo
      ! no meridional flux across South Pole and North Pole
      fyl(:,maxJ) = 0._wp
      afy(:,maxJ) = 0._wp
      fyl(:,0) = 0._wp
      afy(:,0) = 0._wp

      ! STEP I: LOWER ORDER SOLUTION
      do j=1,maxj
        do i=1,maxi
          if (f_ocn(i,j).gt.0._wp) then
            his_tmp(i,j) = his(i,j)-(fxl(I,j)-fxl(I-1,j)+fyl(i,J)-fyl(i,J-1)) &
              / (dx(j)*dy*f_ocn(i,j)) ! ocean area in grid cell
          else
            his_tmp(i,j) = 0._wp
          endif
        enddo
      enddo

      !------------------------------------------------------
      ! 2D flux limiters
      !------------------------------------------------------

      ! apply eq. 14' of Zalesak 1979
      do j=1,maxj
        do i=1,maxi
          ip1=i+1
          if (ip1.eq.maxi+1) ip1=1
          im1=i-1
          if (im1.eq.0) im1=maxi
          ip2=i+2
          if (ip2.eq.maxi+1) ip2=1
          if (ip2.eq.maxi+2) ip2=2
          jm1=max(1,j-1)
          jp1=min(maxj,j+1)
          jp2=min(maxj,j+2)
          if (afx(I,j)*(his_tmp(ip1,j)-his_tmp(i,j)).lt.0._wp &
            .and. (afx(I,j)*(his_tmp(ip2,j)-his_tmp(ip1,j)).lt.0._wp .or. afx(I,j)*(his_tmp(i,j)-his_tmp(im1,j)).lt.0._wp)) then
            afx(I,j) = 0._wp
          endif
          if (afy(i,J)*(his_tmp(i,jp1)-his_tmp(i,j)).lt.0._wp &
            .and. (afy(i,J)*(his_tmp(i,jp2)-his_tmp(i,jp1)).lt.0._wp .or. afy(i,J)*(his_tmp(i,j)-his_tmp(i,jm1)).lt.0._wp)) then
            afy(i,J) = 0._wp
          endif
        enddo
      enddo
      afx(0,:) = afx(maxI,:)
    
      do j=1,maxj
        do i=1,maxi
          if (f_ocn(i,j).gt.0._wp) then
            xa(i,j)=max(his(i,j),his_tmp(i,j))
            xb(i,j)=min(his(i,j),his_tmp(i,j))
          else
            ! land cells must not bound the neighbouring ocean cells
            xa(i,j)=-1.e30_wp
            xb(i,j)=1.e30_wp
          endif
        enddo
      enddo
    
      do j=1,maxj
        do i=1,maxi
          im1=i-1
          if (im1.eq.0) im1=maxi
          ip1=i+1
          if (ip1.eq.maxi+1) ip1=1
          jm1=max(1,j-1)
          jp1=min(maxj,j+1)
          xmax = max(xa(im1,j),xa(ip1,j),xa(i,jm1),xa(i,jp1),xa(i,j))
          xmin = min(xb(im1,j),xb(ip1,j),xb(i,jm1),xb(i,jp1),xb(i,j))
          pp = max(0._wp,afx(Im1,j))-min(0._wp,afx(I,j))+max(0._wp,afy(i,Jm1))-min(0._wp,afy(i,J))
          pm = max(0._wp,afx(I,j))-min(0._wp,afx(Im1,j))+max(0._wp,afy(i,J))-min(0._wp,afy(i,Jm1))
          qp = (xmax-his_tmp(i,j))*(dx(j)*dy*f_ocn(i,j))
          qm = (his_tmp(i,j)-xmin)*(dx(j)*dy*f_ocn(i,j))
          if (pp.lt.1.e-3_wp) then
            rp(i,j) = 0._wp
          else
            qdp = qp/pp
            rp(i,j) = min(1._wp,qdp)
          endif
          if (pm.lt.1.e-3_wp) then
            rn(i,j) = 0._wp
          else
            qdp = qm/pm
            rn(i,j) = min(1._wp,qdp)
          endif
        enddo
      enddo
    
      do j=1,maxj
        do i=1,maxi
          ip1 = i+1
          if (ip1.eq.maxi+1) ip1=1
          if (f_ocn(i,j).eq.0._wp .or. f_ocn(ip1,j).eq.0._wp) then
            afx(I,j) = 0._wp
          else
            if (afx(I,j).ge.0._wp) then
              cx = min(rp(ip1,j),rn(i,j))
            else
              cx = min(rp(i,j),rn(ip1,j))
            endif
            afx(I,j) = cx*afx(I,j)
          endif
        enddo
      enddo
      afx(0,:) = afx(maxI,:)
      do j=1,maxj-1 ! afy(:,0) and afy(:,maxj) are already set =0 above
        do i=1,maxi
          if (f_ocn(i,j).eq.0._wp .or. f_ocn(i,j+1).eq.0._wp) then
            afy(i,J) = 0._wp
          else
            if (afy(i,J).ge.0._wp) then
              cy = min(rp(i,j+1),rn(i,j))
            else
              cy = min(rp(i,j),rn(i,j+1))
            endif
            afy(i,J) = cy*afy(i,J)
          endif
        enddo
      enddo
    
      ! STEP II: HIGH ORDER CORRECTION   
      do j=1,maxj
        do i=1,maxi
          if (f_ocn(i,j).gt.0._wp) then
            his_tmp(i,j) = his_tmp(i,j)-(afx(I,j)-afx(I-1,j)+afy(i,J)-afy(i,J-1)) &
              / (dx(j)*dy*f_ocn(i,j))  ! ocean area in grid cell
          endif
        enddo
      enddo
    
      ! calculate the overall advective fluxes
      do j=1,maxj
        do i=1,maxi
          if (f_ocn(i,j).gt.0._wp) then
            adv_tendency(i,j) = (his_tmp(i,j) - his(i,j)) / dt
          endif
        enddo
      enddo
    
      if (is_sic_snow) then
        fax(0:,:) = fxl(0:,:) + afx(0:,:)
        fay(:,0:) = fyl(:,0:) + afy(:,0:)
      endif

    endif

    !--------------------------------------------------
    ! diffusion
    !--------------------------------------------------

    ! flux to east
    do j=1,maxj
      do i=1,maxi
        ip1 = i+1
        if (ip1.eq.maxi+1) ip1=1
        if (f_ocn(i,j).gt.0._wp .and. f_ocn(ip1,j).gt.0._wp) then
          diffx_max = 0.2_wp * 0.5_wp*dx(j)**2/dt    ! maximum diffusivity for numerical stability
          fdx1(I,j) = - (his(ip1,j) - his(i,j)) * rdx(j) * min(diffx_max,diffsic) * dy*dt  ! m/m*m2/s*m*s = m3
        else
          fdx1(I,j) = 0._wp
        endif
      enddo
    enddo
    fdx1(0,:) = fdx1(maxI,:)

    ! flux to north
    do j=1,maxj-1
      do i=1,maxi
        if (f_ocn(i,j).gt.0._wp .and. f_ocn(i,j+1).gt.0._wp) then
          diffy_max = 0.2_wp * 0.5_wp*dy**2/dt    ! maximum diffusivity for numerical stability
          fdy1(i,J) = - (his(i,j+1) - his(i,j)) * rdy * min(diffy_max,diffsic) * dxv(J)*dt  ! m/m*m2/s*m*s = m3
        else
          fdy1(i,J) = 0._wp
        endif
      enddo
    enddo
    fdy1(:,0) = 0._wp
    fdy1(:,maxJ) = 0._wp

    do j=1,maxj
      do i=1,maxi
        if (f_ocn(i,j).gt.0._wp) then
          diff_tendency(i,j) = - (fdx1(I,j)-fdx1(I-1,j) + fdy1(i,J)-fdy1(i,J-1)) &
            / (dx(j)*dy*f_ocn(i,j)) & ! area of ocean in grid cell
            / dt
        endif
      enddo
    enddo

    !--------------------------------------------------
    ! time step
    !--------------------------------------------------

    do j=1,maxj
      do i=1,maxi
        if (f_ocn(i,j).gt.0._wp) then
          his(i,j) = his(i,j) + dt*(adv_tendency(i,j) + diff_tendency(i,j))
          if(his(i,j).lt.-1.e-1_wp) then
            print *,'WARNING! prognostic sea ice variable < 0'
            print *,'his,his_old,i,j',his(i,j),his_old(i,j),i,j
            if (i_adv.eq.2) print *,'his_tmp',his_tmp(i,j)
            print *,'adv,diff',dt*adv_tendency(i,j),dt*diff_tendency(i,j)
            print *,'u,v',u(i,j),v(i,j)
            print *,'fdx,fdy',-(fdx1(I,j)-fdx1(I-1,j)),-(fdy1(i,J)-fdy1(i,J-1))
            if (i_adv.eq.2) print *,'afx,afy',-(afx(I,j)-afx(I-1,j)),-(afy(i,J)-afy(i,J-1))
          endif
        endif
      enddo
    enddo

    if (is_sic_snow) then
      ! calculate dynamical component 
      ! this is purely diagnostic
      where (f_ocn.gt.0._wp)
        dh_dt_dyn = (his - his_old)/dt 
      elsewhere
        dh_dt_dyn = 0._wp
      endwhere
      ! fluxes for output
      fdx(0:,:) = fdx1(0:,:)
      fdy(:,0:) = fdy1(:,0:)
    endif


  return

  end subroutine transport

end module transport_sic_mod
