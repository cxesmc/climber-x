!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : u 3 d _ m o d
!
!  Purpose : 3-D wind field
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
module u3d_mod

  use atm_params, only : wp
  use constants, only : g, T0, r_earth
  use atm_params, only : amas, ra, i_pbl, dpc, pcmin, pcmax, ptopdyn
  use atm_params, only : l_mass_com_topo
  use atm_params, only : c_uter_pol, c_uter_eq
  use atm_params, only : l_output_flx3d
  use atm_grid, only : im, imc, jm, jmc, km, kmc, k500, k700, dxt, dxu, dy, zl, sqr
  use atm_grid, only : fcort, cost, sint
  use atm_grid, only : pl, dplx, dply, dplxo, dplyo, kplxo, kplyo, pblt, pblu

  implicit none
  
  private
  public :: u3d

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  u 3 d
  !   Purpose    :  computation of 3D wind field and advective mass fluxes
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine u3d(niter, pzsa, ptrop, ugb, vgb, uab, vab, t3, &
        ua, va, uter, vter, uterf, vterf, u3, v3, w3, uz500, &
        fax, faxo, fay, fayo, fac)

    implicit none

    integer,  intent(in   ) :: niter
    real(wp), intent(in   ) :: pzsa(:,:)
    real(wp), intent(in   ) :: ptrop(:)
    real(wp), intent(in   ) :: ugb(:,:)
    real(wp), intent(in   ) :: vgb(:,:)
    real(wp), intent(in   ) :: uab(:,:)
    real(wp), intent(in   ) :: vab(:,:)
    real(wp), intent(in   ) :: t3(:,:,:)

    real(wp), intent(out  ) :: ua(:,:,:)
    real(wp), intent(out  ) :: va(:,:,:)
    real(wp), intent(out  ) :: uter(:,:,:)
    real(wp), intent(out  ) :: vter(:,:,:)
    real(wp), intent(out  ) :: uterf(:,:,:)
    real(wp), intent(out  ) :: vterf(:,:,:)
    real(wp), intent(inout) :: u3(:,:,:)
    real(wp), intent(inout) :: v3(:,:,:)
    real(wp), intent(inout) :: w3(:,:,:)
    real(wp), intent(inout) :: uz500(:)
    real(wp), intent(out  ) :: fax(:,:,:)
    real(wp), intent(out  ) :: faxo(:,:,:)
    real(wp), intent(out  ) :: fay(:,:,:)
    real(wp), intent(out  ) :: fayo(:,:,:)
    real(wp), intent(out  ) :: fac(:,:)

    integer :: i, j, k, n, ipl, imi, jmi, kpl
    real(wp) :: pzx, pzy, dp_c, pbl_t, pbl_u, dp, pc1, pc2, dp_l, pl1, pl2, fxpbl, fypbl, uabc, vabc, ctv
    real(wp) :: fai, dfa, faz
    real(wp) :: u_g, v_g, ug_b, vg_b
    real(wp) :: u3k, u3kp1, v3k, v3kp1
    real(wp) :: dplxdy, dplydxu
    real(wp) :: c_damp_eq, c_damp_pol

    real(wp), dimension(kmc) :: uteri
    real(wp), dimension(kmc) :: vteri 


    !$omp parallel do collapse(2) private(i,j,k,ipl,imi,jmi,pzx,pzy,dp,pbl_t,pbl_u,pc1,pc2,dp_l,pl1,pl2,fxpbl,fypbl,uabc,vabc,uteri,vteri,ctv,c_damp_eq,c_damp_pol)
    do j=1,jm
      do i=1,im

        ipl = modulo(i,im) + 1
        imi = modulo(i - 2, im) + 1
        jmi = max(1,j-1)

        !-------------------------------------------------------
        ! Vertical profile of ageostrophic wind

        ! x-component on u-points 

        ! velocity in the PBL
        if (i_pbl.eq.1) then
          pzx = 0.5_wp*(pzsa(i,j)+pzsa(imi,j))
        else if (i_pbl.eq.2) then
          pzx = 1._wp
        endif
        pbl_t = pzx+pblt(j)-1._wp
        do k=1,km
          dp = (pl(k)-pbl_t)/(pl(k)-pl(k+1))
          dp = min(1._wp,dp)
          dp = max(0._wp,dp)
          ua(i,j,k) = uab(i,j) * dp
        enddo

        ! compensatory velocity in the upper troposphere
        if (i_pbl.eq.1) then
          fxpbl = uab(i,j)*(1._wp-pblt(j))        
        else if (i_pbl.eq.2) then
          fxpbl = uab(i,j)*max(0._wp,pzx-pblt(j))        
        endif
        uabc = -fxpbl/dpc
        pc2 = ptrop(j)
        pc1 = pc2+dpc
        do k=1,km
          pl1 = max(min(pl(k),pc1),pl(k+1))
          pl2 = min(max(pl(k+1),pc2),pl(k))
          dp_l = pl1-pl2
          ua(i,j,k) = ua(i,j,k)+uabc*dp_l/(pl(k)-pl(k+1))
        enddo     

        ! y-component on v-points 

        ! velocity in the PBL
        if (i_pbl.eq.1) then
          pzy = 0.5_wp*(pzsa(i,j)+pzsa(i,jmi))
        else if (i_pbl.eq.2) then
          pzy = 1._wp
        endif
        pbl_u = pzy+pblu(j)-1._wp
        do k=1,km
          dp = (pl(k)-pbl_u)/(pl(k)-pl(k+1))
          dp = min(1._wp,dp)
          dp = max(0._wp,dp)
          va(i,j,k) = vab(i,j) * dp
        enddo

        ! compensatory velocity in the upper troposphere
        if (i_pbl.eq.1) then
          fypbl = vab(i,j)*(1._wp-pblu(j))        
        else if (i_pbl.eq.2) then
          fypbl = vab(i,j)*max(0._wp,pzy-pblu(j))        
        endif
        vabc = -fypbl/dpc
        pc2 = 0.5_wp*(ptrop(j)+ptrop(jmi))  
        pc1 = pc2+dpc
        do k=1,km
          pl1 = max(min(pl(k),pc1),pl(k+1))
          pl2 = min(max(pl(k+1),pc2),pl(k))
          dp_l = pl1-pl2
          va(i,j,k) = va(i,j,k)+vabc*dp_l/(pl(k)-pl(k+1))
        enddo  

        !-------------------------------------------------------
        ! Thermal wind on T-points

        if (j.eq.1 .or. j.eq.jm) then

          uter(i,j,:) = 0._wp
          vter(i,j,:) = 0._wp
          uterf(i,j,:) = 0._wp
          vterf(i,j,:) = 0._wp

        else

          ! thermal wind on levels
          uteri(1) = 0._wp
          vteri(1) = 0._wp
          uteri(km+1) = 0._wp
          vteri(km+1) = 0._wp
          do k=1,km-1
            ctv = (zl(k+1)-zl(k))*g/(T0*fcort(j)) 
            uteri(k+1) = uteri(k)-ctv*(t3(i,jmi,k)-t3(i,j+1,k))/(2._wp*dy)   
            vteri(k+1) = vteri(k)+ctv*(t3(ipl,j,k)-t3(imi,j,k))/(2._wp*dxt(j)) 
          enddo

          ! thermal wind in layers, dampened at equator and poles
          c_damp_pol = min(1._wp,c_uter_pol*cost(j)**2)
          c_damp_eq  = min(1._wp, c_uter_eq*sint(j)**2)
          do k=1,km
            uter(i,j,k) = 0.5_wp*(uteri(k)+uteri(k+1)) * c_damp_eq * c_damp_pol 
            vter(i,j,k) = 0.5_wp*(vteri(k)+vteri(k+1)) * c_damp_eq * c_damp_pol 
            ! thermal wind without polar damping needed for EKE production
            uterf(i,j,k) = 0.5_wp*(uteri(k)+uteri(k+1)) * c_damp_eq 
            vterf(i,j,k) = 0.5_wp*(vteri(k)+vteri(k+1)) * c_damp_eq
          enddo

        endif

      enddo
    enddo      
    !$omp end parallel do


    !-------------------------------------------------------
    ! advective mass transport 
    !-------------------------------------------------------

    !$omp parallel do collapse(2) private(i,j,k,imi,u_g,v_g,dplxdy,dplydxu)
    do k=1,km
      do j=1,jm

        ! x-components

        do i=1,im
          imi = modulo(i - 2, im) + 1
          ! geostrophic zonal wind on u-points, limited to troposphere
          if (pl(k+1).ge.ptopdyn) then
            u_g  = 0.5_wp*(ugb(imi,j)+ugb(i,j)) + 0.5_wp*(uter(imi,j,k)+uter(i,j,k))
          elseif (pl(k+1).lt.ptopdyn.and.pl(k).ge.ptopdyn) then
            u_g  = (0.5_wp*(ugb(imi,j)+ugb(i,j)) + 0.5_wp*(uter(imi,j,k)+uter(i,j,k))) *(pl(k)-ptopdyn)/(pl(k)-pl(k+1))
          else
            u_g  = 0._wp
          endif 
          dplxdy = dplx(i,j,k)*dy
          ! mass flux
          fax(i,j,k) = (u_g+ua(i,j,k))*dplxdy ! m/s * kg/m2 * m = kg/s
          ! orographic component of mass flux 
          faxo(i,j,k) = u_g*(dplxo(i,j,k)*dy-dplxdy) 
        enddo
        ! periodic boundary conditions
        fax(imc,j,k) = fax(1,j,k)
        faxo(imc,j,k) = faxo(1,j,k)

        ! y-components

        do i=1,im
          if (j.eq.1) then
            ! N and S boundary conditions, no flux
            fay(i,1,k)  = 0._wp
            fayo(i,1,k) = 0._wp
            fay(i,jmc,k)  = 0._wp
            fayo(i,jmc,k) = 0._wp
          else
            ! geostrophic meridional wind on u-points, limited to troposphere
            if (pl(k+1).ge.ptopdyn) then
              v_g  = 0.5_wp*(vgb(i,j-1)+vgb(i,j)) + 0.5_wp*(vter(i,j-1,k)+vter(i,j,k))
            elseif (pl(k+1).lt.ptopdyn.and.pl(k).ge.ptopdyn) then
              v_g  = (0.5_wp*(vgb(i,j-1)+vgb(i,j)) + 0.5_wp*(vter(i,j-1,k)+vter(i,j,k))) *(pl(k)-ptopdyn)/(pl(k)-pl(k+1))
            else
              v_g  = 0._wp
            endif 
            dplydxu = dply(i,j,k)*dxu(j)
            ! mass flux
            fay(i,j,k) = (v_g+va(i,j,k))*dplydxu  ! m/s * kg/m2 * m = kg/s
            ! orographic component of mass flux for temperature
            fayo(i,j,k) = v_g*(dplyo(i,j,k)*dxu(j)-dplydxu) 
          endif
        enddo

      enddo
    enddo 
    !$omp end parallel do

    !-------------------------------------------------------
    ! mass transport compensation
    !-------------------------------------------------------

    if (l_mass_com_topo) then

      ! deal with vertical mass flux generated by orography
      ! make mass flux follow orography
      !$omp parallel do collapse(2) private(i,j,k)
      do j=1,jm
        do i=1,im
          do k=1,k500
            ! add mass flux in first layer above topography
            fax(i,j,kplxo(i,j)) = fax(i,j,kplxo(i,j)) + faxo(i,j,k)
            fay(i,j,kplyo(i,j)) = fay(i,j,kplyo(i,j)) + fayo(i,j,k)
          enddo
        enddo
      enddo
      !$omp end parallel do

    endif


    if (l_output_flx3d .and. niter.eq.1) then
      ! column convergence before compensation
      do j=1,jm
        do i=1,im
          fac(i,j) = 0._wp
          do k=1,km 
            fac(i,j) = fac(i,j) + fax(i,j,k) -fax(i+1,j,k) +fay(i,j+1,k) -fay(i,j,k)
          enddo
        enddo
      enddo
    endif

    ! enforce that zonal mean meridional mass flux is zero
    !$omp parallel do private(i,j,k,n,dfa)
    do j=1,jm
      ! zonal and vertical mean 
      dfa = 0._wp
      n = 0
      do k=1,km
        do i=1,im
          if (dply(i,j,k).gt.0._wp) then
            dfa = dfa + fay(i,j,k)
            n = n+1
          endif
        enddo
      enddo
      if (n.gt.0) then
        dfa = dfa/n   
        do k=1,km
          do i=1,im
            if (dply(i,j,k).gt.0._wp) then
              fay(i,j,k) = fay(i,j,k) - dfa
            endif
          enddo
        enddo
      endif
    enddo
    !$omp end parallel do

    ! apply mass compensation to make sure that each column is non-divergent
    !$omp parallel do private(i,j,k,ipl,fai,pzx,dp_c,pc1,pc2,dp_l,pl1,pl2)
    do j=1,jm
      do i=1,im

        ipl = modulo(i,im) + 1

        ! column integrated convergence
        fai  = 0._wp
        do k=1,km 
          fai = fai + fax(i,j,k)-fax(i+1,j,k) + fay(i,j+1,k)-fay(i,j,k)
        enddo
        ! compensate in zonal component somewhere in the upper troposphere
        pzx = 0.5_wp*(pzsa(ipl,j)+pzsa(i,j))
        do k=1,km
          if (pl(k).lt.pzx) exit
        enddo
        pc2 = max(pcmin,ptrop(j))
        pc1 = min(pcmax,pl(k))
        pc1 = max(pc2+0.01_wp,pc1)
        dp_c = pc1-pc2 
        do k=1,km
          pl1 = max(min(pl(k),pc1),pl(k+1))
          pl2 = min(max(pl(k+1),pc2),pl(k))
          dp_l = pl1-pl2
          fax(i+1,j,k)=fax(i+1,j,k)+fai*dp_l/dp_c
        enddo     

      enddo
    enddo
    !$omp end parallel do


    if (niter.eq.1) then

      ! total 3-D wind on T-points, interpolate to levels

      !$omp parallel do collapse(2) private(i,j,ipl,k,kpl,ug_b,vg_b,u3k,v3k,u3kp1,v3kp1,faz)
      do j=1,jm
        do i=1,im

          ipl = modulo(i,im) + 1

          ug_b = ugb(i,j)    
          vg_b = vgb(i,j)
          u3(i,j,1) = 0.5_wp*(ua(i,j,1)+ua(ipl,j,1)) + ug_b
          v3(i,j,1) = 0.5_wp*(va(i,j,1)+va(i,min(j+1,jm),1)) + vg_b
          do k=1,km-1     
            u3k   = 0.5_wp*(ua(i,j,k)+ua(ipl,j,k)) + ug_b + uter(i,j,k)
            v3k   = 0.5_wp*(va(i,j,k)+va(i,min(j+1,jm),k)) + vg_b + vter(i,j,k)
            u3kp1 = 0.5_wp*(ua(i,j,k+1)+ua(ipl,j,k+1)) + ug_b + uter(i,j,k+1)
            v3kp1 = 0.5_wp*(va(i,j,k+1)+va(i,min(j+1,jm),k+1)) + vg_b + vter(i,j,k+1)
            u3(i,j,k+1) = 0.5_wp*(u3k+u3kp1)
            v3(i,j,k+1) = 0.5_wp*(v3k+v3kp1) 
          enddo

          ! z-component and vertical velocity        

          faz = 0._wp
          w3(i,j,1) = 0._wp

          do k=1,km 
            kpl = min(k+1,km)
            faz = faz + fax(i,j,k)-fax(i+1,j,k) + fay(i,j+1,k)-fay(i,j,k)
            ! vertical velocity
            w3(i,j,k+1) = faz/(sqr(i,j)*ra*pl(kpl))    ! kg/s /m2 *m3/kg = m/s
          enddo

        enddo
      enddo
      !$omp end parallel do

      ! zonal mean 500 hPa zonal wind
      do j=1,jm
        uz500(j) = 0._wp
        do i=1,im
          uz500(j) = uz500(j)+u3(i,j,k500)/im
        enddo
      enddo

    endif


    return

  end subroutine u3d

end module u3d_mod
