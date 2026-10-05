!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module : d o w n s l o p e _ m o d
!
!  Purpose : density-driven downsloping flow of dense shelf water
!            (Campin and Goosse 1999, Tellus 51A)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2026 Potsdam Institute for Climate Impact Research,
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
module downslope_mod

  ! Dense water formed by convection on the continental shelves (Antarctica, Nordic Seas) cannot
  ! reach the abyss in a z-coordinate model at this resolution: it leaves the shelf column
  ! laterally at the depth of the shelf bottom (500-1000 m), and from there the only way down is
  ! convective adjustment of the neighbouring deep column, which mixes it into the whole column,
  ! i.e. the maximum-entrainment limit. This routine provides the missing pathway following
  ! Campin and Goosse (1999):
  !
  !  - for every ocean cell whose bottom is shallower than that of a horizontal neighbour, the water
  !    in its bottom layer is compared with the neighbour's water at the same level. If it is
  !    denser (by more than drho_downslope_min), a downslope volume flux
  !        Q = c_downslope * sqrt(g' H) * H * L,   g' = g*drho/rho0
  !    is taken from the shelf bottom layer (H its thickness, L the length of the shared face;
  !    a gravity-current speed sqrt(g'H) times a Froude-like coefficient that also stands for the
  !    fraction of the 5 deg face occupied by the plume),
  !  - the plume descends the deep column level by level. At each level it entrains a fraction
  !    ent_downslope*dz/1000 of its own volume from the ambient water, and it continues as long as
  !    it is still denser than the ambient water AT THE PRESSURE OF THAT LEVEL (the thermobaric
  !    effect is what makes cold shelf water sink to the bottom: a comparison at the shelf
  !    pressure would stop the plume at 2 km),
  !  - the plume is deposited at the deepest level where it is still denser (its neutral
  !    buoyancy level, or the bottom); the compensating return flow rises through the deep
  !    column, feeds the entrainment on the way, and re-enters the shelf cell at the shelf bottom
  !    level. All exchanges are volume-conserving, so tracers are conserved to round-off.
  !
  ! The exchanged volume per time step (including entrainment) is capped at f_downslope_max of the
  ! shelf bottom cell and of every deep cell it passes through, which keeps the explicit update bounded
  ! with the one-day ocean time step (a 5x5 deg, 200 m cell holds 1.6e13 m3, i.e. 1 Sv moves 0.5% of it per day).

  use precision, only : wp
  use constants, only : g
  use ocn_grid, only : maxi, maxj, maxk, k1, mask_ocn, zro, dz, dy, dxv, ocn_vol
  use ocn_params, only : dt, rho0, n_tracers_tot
  use ocn_params, only : c_downslope, ent_downslope, drho_downslope_min, f_downslope_max
  use eos_mod

  implicit none

  private
  public :: downslope

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  d o w n s l o p e
  !   Purpose    :  downsloping flow of dense shelf water into the neighbouring deep columns
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine downslope(l_trans_tracers,f_ocn,ts,fds,zds)

    implicit none

    logical, dimension(:), intent(in) :: l_trans_tracers
    real(wp), dimension(:,:), intent(in) :: f_ocn
    real(wp), dimension(:,:,:,:), intent(inout) :: ts
    real(wp), dimension(:,:), intent(out) :: fds    !! downslope volume flux out of the shelf cell [m3/s]
    real(wp), dimension(:,:), intent(out) :: zds    !! flux-weighted depth at which the plume is deposited [m]

    integer :: i, j, n, ii, jj, ks, kd, k, kt, l
    real(wp) :: rho_s, rho_d, drho, gp, H, flen, q, vex, vp, vent, w, vcell
    real(wp), dimension(n_tracers_tot) :: cp, cs, cp_prev
    real(wp) :: vp_prev
    real(wp), dimension(maxk) :: vin, w_up
    real(wp), dimension(maxk,n_tracers_tot) :: ts_old
    integer, dimension(4) :: di = (/1,-1,0,0/), dj = (/0,0,1,-1/)

    fds(:,:) = 0._wp
    zds(:,:) = 0._wp

    do j=1,maxj
      do i=1,maxi
        if (mask_ocn(i,j).ne.1) cycle
        ks = k1(i,j)    ! bottom level of the (potential) shelf cell

        do n=1,4
          ii = i+di(n)
          jj = j+dj(n)
          if (jj.lt.1 .or. jj.gt.maxj) cycle
          ii = modulo(ii-1,maxi)+1
          if (mask_ocn(ii,jj).ne.1) cycle
          kd = k1(ii,jj)
          if (kd.ge.ks) cycle   ! neighbour is not deeper

          ! density of the shelf bottom water and of the neighbour at the same level
          rho_s = eos(ts(i,j,ks,1),ts(i,j,ks,2),zro(ks))
          rho_d = eos(ts(ii,jj,ks,1),ts(ii,jj,ks,2),zro(ks))
          drho = rho_s-rho_d
          if (drho.le.drho_downslope_min) cycle

          ! downslope volume flux across the shared face
          gp = g*drho/rho0
          H  = dz(ks)
          if (n.le.2) then
            flen = dy               ! zonal neighbour: meridional face
          else
            flen = dxv(min(j,jj))   ! meridional neighbour: zonal face at the cell edge
          endif
          flen = flen*min(f_ocn(i,j),f_ocn(ii,jj))
          q = c_downslope*sqrt(gp*H)*H*flen      ! m3/s
          vex = min(q*dt, f_downslope_max*ocn_vol(i,j,ks))   ! m3 exchanged this step

          ! old tracer values of the deep column (all fluxes use them, so the update is conservative)
          ts_old(kd:maxk,:) = ts(ii,jj,kd:maxk,:)
          cs(:) = ts(i,j,ks,:)

          ! descend the deep column: plume properties cp, volume vp; vin(k) = volume entrained at level k.
          ! Entrainment at a level is only kept once the plume is found to continue below it, so the
          ! deposition level itself entrains nothing (its water is where the plume ends up anyway).
          cp(:) = cs(:)
          vp = vex
          vin(:) = 0._wp
          kt = ks
          do k=ks-1,kd,-1
            cp_prev(:) = cp(:)
            vp_prev = vp
            if (kt.lt.ks) then
              ! tentative entrainment at the level above (kt = k+1), an intermediate level if the plume goes on
              vent = ent_downslope*dz(kt)/1000._wp*vp
              do l=1,n_tracers_tot
                cp(l) = (cp(l)*vp + ts_old(kt,l)*vent)/(vp+vent)
              enddo
              vp = vp+vent
              vin(kt) = vent
            endif
            rho_s = eos(cp(1),cp(2),zro(k))
            rho_d = eos(ts_old(k,1),ts_old(k,2),zro(k))
            if (rho_s.le.rho_d) then
              ! the plume stops at kt: undo the tentative entrainment there
              cp(:) = cp_prev(:)
              vp = vp_prev
              if (kt.lt.ks) vin(kt) = 0._wp
              exit
            endif
            kt = k
          enddo
          if (kt.eq.ks) cycle   ! not denser than the neighbour one level down: no downslope flow

          ! keep every exchange below f_downslope_max of the cells involved: the plume volume vp is the largest
          ! flow (it passes the deposition level), so scale the whole circulation down if needed. A uniform
          ! scaling leaves the plume's mixing ratios, hence cp and the deposition level, unchanged.
          vcell = ocn_vol(i,j,ks)
          do k=kt,ks
            vcell = min(vcell,ocn_vol(ii,jj,k))
          enddo
          if (vp.gt.f_downslope_max*vcell) then
            w = f_downslope_max*vcell/vp
            vex = vex*w
            vp  = vp*w
            vin(:) = vin(:)*w
          endif

          ! return flow up the deep column from the deposition level to the shelf bottom level:
          ! w_up(k) = volume leaving level k upward = volume arriving from below minus what the plume entrained at k
          w = vp
          do k=kt,ks-1
            w_up(k) = w
            w = w-vin(k+1)   ! vin(ks) = 0 by construction
          enddo
          ! the level below ks is the last with an upward flow; what arrives at ks equals vex

          ! tracer updates (explicit, old values on the right-hand side)
          do l=1,n_tracers_tot
            if (.not.l_trans_tracers(l)) cycle
            ! deposition level: gains the plume, loses w_up(kt) upward
            vcell = ocn_vol(ii,jj,kt)
            ts(ii,jj,kt,l) = ts(ii,jj,kt,l) + (cp(l)*vp - ts_old(kt,l)*w_up(kt))/vcell
            ! intermediate levels: gain from below, lose entrainment to the plume and the upward flow
            do k=kt+1,ks-1
              vcell = ocn_vol(ii,jj,k)
              ts(ii,jj,k,l) = ts(ii,jj,k,l) + (ts_old(k-1,l)*w_up(k-1) - ts_old(k,l)*(vin(k)+w_up(k)))/vcell
            enddo
            ! shelf bottom level of the deep column: gains from below, loses vex to the shelf cell
            vcell = ocn_vol(ii,jj,ks)
            ts(ii,jj,ks,l) = ts(ii,jj,ks,l) + (ts_old(ks-1,l)*w_up(ks-1) - ts_old(ks,l)*vex)/vcell
            ! shelf cell: loses vex of its own water, gains vex from the deep column at the same level
            vcell = ocn_vol(i,j,ks)
            ts(i,j,ks,l) = ts(i,j,ks,l) + (ts_old(ks,l) - cs(l))*vex/vcell
          enddo

          ! diagnostics
          fds(i,j) = fds(i,j) + vex/dt
          zds(i,j) = zds(i,j) + vex/dt*zro(kt)

        enddo   ! neighbours

        if (fds(i,j).gt.0._wp) zds(i,j) = zds(i,j)/fds(i,j)

      enddo
    enddo

    return

  end subroutine downslope

end module downslope_mod
