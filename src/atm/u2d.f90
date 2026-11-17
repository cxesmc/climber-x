!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : u 2 d _ m o d
!
!  Purpose : PBL and surface wind
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
module u2d_mod

  use atm_params, only : wp
  use constants, only : g
  use atm_params, only : ra, l_kata_wind, h_kata, c_kata, l_kata_upw, hgams
  use atm_grid, only : im, imc, jm, jmc, nm, dxt, dy
  use atm_grid, only : fcorg, fcorta, fcorua, signf
  !$use omp_lib

  implicit none
  
  private
  public :: u2d, usur

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  u 2 d
  !   Purpose    :  computation of geostrophic and ageostrophic wind in 
  !                 planetary boundary layer
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine u2d(slp, sin_cos_acbar, &
      ugb, vgb, uab, vab)

    implicit none

    real(wp), intent(in ) :: slp(:,:)
    real(wp), intent(in ) :: sin_cos_acbar(:,:)

    real(wp), intent(out) :: ugb(:,:)
    real(wp), intent(out) :: vgb(:,:)
    real(wp), intent(out) :: uab(:,:)
    real(wp), intent(out) :: vab(:,:)

    integer :: i, j, ipl, imi
    real(wp) :: dpdx, dpdy, dpdxa, dpdya, acbarb


    !$omp parallel do collapse(2) private(i,j,ipl,imi,dpdx,dpdy,dpdxa,dpdya,acbarb)
    do j=1,jm
      do i=1,im 

        ipl = modulo(i,im) + 1
        imi = modulo(i - 2, im) + 1

        !------------------------------------------------------
        ! PBL geostrophic wind components in T-points

        ! Horizontal SLP gradient in T-points      
        dpdx = 0.5_wp*(slp(ipl,j)-slp(imi,j))/dxt(j) 
        if (j.eq.1) then
          ! centred through the pole, using the cell at the opposite longitude
          dpdy = 0.5_wp*(slp(modulo(i-1+im/2,im)+1,1)-slp(i,2))/dy
        else if (j.eq.jm) then
          dpdy = 0.5_wp*(slp(i,jm-1)-slp(modulo(i-1+im/2,im)+1,jm))/dy
        else
          dpdy = 0.5_wp*(slp(i,j-1)-slp(i,j+1))/dy
        endif
        ! fcorg is 1/fcort for i_fcorg=0 and the regular f/(f**2+fcormin**2) for i_fcorg=1.
        ! See the fcorg block in atm_grid.f90.
        ugb(i,j) = -dpdy*fcorg(j)/ra
        vgb(i,j) =  dpdx*fcorg(j)/ra

        !------------------------------------------------------
        ! Ageostrophic wind components in PBL (on U-points)

        dpdxa = (slp(i,j)-slp(imi,j))/dxt(j)
        acbarb = 0.5_wp*(sin_cos_acbar(imi,j)+sin_cos_acbar(i,j))
        uab(i,j) = -1._wp/(fcorta(j)*ra)*(dpdxa*acbarb)

        if (j.gt.1) then
          dpdya = (slp(i,j-1)-slp(i,j))/dy
          acbarb = 0.5_wp*(sin_cos_acbar(i,j-1)+sin_cos_acbar(i,j))
          vab(i,j) = -1._wp/(fcorua(j)*ra)*(dpdya*acbarb)
        endif

        if (j.eq.jm) then
          vab(i,1)   = 0._wp
          vab(i,jmc) = 0._wp
        endif

        ! periodic boundary conditions
        if (i.eq.1) then
          uab(imc,j) = uab(1,j)
        endif

      enddo
    enddo 
    !$omp end parallel do

    return

  end subroutine u2d


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  u s u r
  !   Purpose    :  computation of near-surface wind components
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine usur(ugb, vgb, epsa, cos_acbar, sin_acbar, tskina, cd0a, &
      frst, frlnd, zsa, gams, gamb, synsur, &
      usk, vsk, &
      us, vs)

    implicit none

    real(wp), intent(in   ) :: ugb(:,:)
    real(wp), intent(in   ) :: vgb(:,:)
    real(wp), intent(in   ) :: epsa(:,:,:)
    real(wp), intent(in   ) :: cos_acbar(:,:,:)
    real(wp), intent(in   ) :: sin_acbar(:,:,:)
    real(wp), intent(in   ) :: tskina(:,:)
    real(wp), intent(in   ) :: cd0a(:,:)
    real(wp), intent(in   ) :: frst(:,:,:)
    real(wp), intent(in   ) :: frlnd(:,:)
    real(wp), intent(in   ) :: zsa(:,:)
    real(wp), intent(in   ) :: gams(:,:)
    real(wp), intent(in   ) :: gamb(:,:)
    real(wp), intent(in   ) :: synsur(:,:,:)

    real(wp), intent(inout) :: usk(:,:)
    real(wp), intent(inout) :: vsk(:,:)

    real(wp), intent(out  ) :: us(:,:,:)
    real(wp), intent(out  ) :: vs(:,:,:)

    integer :: i, j, n, ipl, imi, iter
    real(wp) :: dtheta
    real(wp) :: fx, fy, fw, fe, fn, fs, fcor, ul, vl, syn2, ukk, vkk, w, r, den

    real(wp), dimension(im,jm)  :: uk, vk
    real(wp), dimension(im,jm)  :: bk

    ! buoyancy of the inversion layer, g*dtheta/theta0, for the upslope forcing
    if (l_kata_wind .and. l_kata_upw) then
      !$omp parallel do collapse(2) private(i,j)
      do j=1,jm
        do i=1,im
          bk(i,j) = g*c_kata*max(gamb(i,j)-gams(i,j),0._wp)*hgams*(1._wp-h_kata/(2._wp*hgams))/tskina(i,j)
        enddo
      enddo
      !$omp end parallel do
    endif

    !$omp parallel do collapse(2) private(i,j,n,ipl,imi,iter,dtheta,fx,fy,fw,fe,fn,fs,fcor,ul,vl,syn2,ukk,vkk,w,r,den)
    do j=1,jm
      do i=1,im 

        ipl = modulo(i,im) + 1
        imi = modulo(i - 2, im) + 1

        if (j.gt.1 .and. j.lt.jm) then

          !------------------------------------------------------
          ! Near surface wind components in T-points
          do n=1,nm
            us(i,j,n) = epsa(i,j,n) * (ugb(i,j)*cos_acbar(i,j,n) - signf(j)*vgb(i,j)*sin_acbar(i,j,n))
            vs(i,j,n) = epsa(i,j,n) * (vgb(i,j)*cos_acbar(i,j,n) + signf(j)*ugb(i,j)*sin_acbar(i,j,n))
          enddo

        endif

        !------------------------------------------------------
        ! katabatic surface wind 

        uk(i,j) = 0._wp
        vk(i,j) = 0._wp

        if (l_kata_wind) then

          ! inversion-layer slab of depth h_kata with Coriolis (Ball 1960, Mahrt 1982), complex notation V = u + i v:
          !   0 = F - i f V_k - r V_k  ->  V_k = F / (r + i f),  r = cd0a*W/h_kata
          !   F = g * dtheta/theta0 * (-grad zsa), downslope buoyancy force, T-points, unclipped topography
          !   dtheta = c_kata * layer-mean deficit vs the free atmosphere extrapolated to the ground (vesta profile):
          !            surface deficit (gamb-gams)*hgams, linear to 0 at hgams -> mean over h_kata: *(1-h_kata/(2 hgams))
          !   W = total wind speed incl. synoptic part, sqrt(|V_ls+V_k|^2 + synsur^2); fixed-point iteration from the previous V_k
          if (j.gt.1 .and. j.lt.jm .and. frlnd(i,j).gt.0._wp) then
            if (l_kata_upw) then
              ! slope on the four faces, each with the buoyancy of its upslope cell (j=1 is north), averaged to the T-point:
              ! a mixed coastal cell feels the air draining from the ice sheet, not its own deficit lowered by the open water
              fw = -(zsa(i,j)-zsa(imi,j))/dxt(j)
              fe = -(zsa(ipl,j)-zsa(i,j))/dxt(j)
              fn = -(zsa(i,j-1)-zsa(i,j))/dy
              fs = -(zsa(i,j)-zsa(i,j+1))/dy
              fx = 0.5_wp*(fw*merge(bk(imi,j),bk(i,j),fw.gt.0._wp) + fe*merge(bk(i,j),bk(ipl,j),fe.gt.0._wp))
              fy = 0.5_wp*(fn*merge(bk(i,j),bk(i,j-1),fn.gt.0._wp) + fs*merge(bk(i,j+1),bk(i,j),fs.gt.0._wp))
            else
              dtheta = c_kata*max(gamb(i,j)-gams(i,j),0._wp)*hgams*(1._wp-h_kata/(2._wp*hgams))
              fx = -(zsa(ipl,j)-zsa(imi,j))/(2._wp*dxt(j))     ! eastward component of -grad zsa
              fy = -(zsa(i,j-1)-zsa(i,j+1))/(2._wp*dy)        ! northward component of -grad zsa (j=1 is north)
              fx = g*dtheta/tskina(i,j)*fx
              fy = g*dtheta/tskina(i,j)*fy
            endif
            if (fx.ne.0._wp .or. fy.ne.0._wp) then
              fcor = signf(j)*fcorta(j)
              ! large-scale wind (grid-cell mean) and synoptic variance
              ul = sum(us(i,j,:)*frst(i,j,:))
              vl = sum(vs(i,j,:)*frst(i,j,:))
              syn2 = sum(synsur(i,j,:)**2*frst(i,j,:))
              ukk = usk(i,j)
              vkk = vsk(i,j)
              do iter=1,5
                w = max(sqrt((ul+ukk)**2+(vl+vkk)**2+syn2),1._wp)
                r = cd0a(i,j)*w/h_kata
                den = r**2+fcor**2
                ukk = (fx*r+fy*fcor)/den
                vkk = (fy*r-fx*fcor)/den
              enddo
              uk(i,j) = ukk
              vk(i,j) = vkk
            endif
          endif

        endif

      enddo
    enddo 
    !$omp end parallel do

    ! Wind near the poles
    us(:,1,:)  = 0.5_wp*us(:,2,:)       
    us(:,jm,:) = 0.5_wp*us(:,jm-1,:) 
    vs(:,1,:)  = 0.5_wp*vs(:,2,:)       
    vs(:,jm,:) = 0.5_wp*vs(:,jm-1,:)

    !$omp parallel do collapse(2) private(i,j,n)
    do j=1,jm
      do i=1,im 

        ! relax in time
        usk(i,j) = 0.1_wp*uk(i,j) + 0.9_wp*usk(i,j)
        vsk(i,j) = 0.1_wp*vk(i,j) + 0.9_wp*vsk(i,j)

        ! add katabatic wind to surface wind
        do n=1,nm
          us(i,j,n) = us(i,j,n) + usk(i,j)
          vs(i,j,n) = vs(i,j,n) + vsk(i,j)
        enddo

      enddo
    enddo
    !$omp end parallel do

    return

  end subroutine usur

end module u2d_mod
