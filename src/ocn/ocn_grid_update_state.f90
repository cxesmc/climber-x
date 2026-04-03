!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : o c n _ g r i d _ u p d a t e _ s t a t e _ m o d
!
!  Purpose : update ocean state following grid update
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
module ocn_grid_update_state_mod

  use precision, only : wp
  use ocn_params, only : n_tracers_tot, i_age, i_dye
  use ocn_grid, only : grid_class, maxi, maxj, maxk, mask_u, mask_v, mask_w

  implicit none

  private
  public :: ocn_grid_update_state

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  o c n _ g r i d _ u p d a t e _ s t a t e
  ! Purpose  :  update ocean state after changes in land/sea mask
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocn_grid_update_state(grid,tracers,u)

    implicit none

    type(grid_class), intent(in) :: grid
    real(wp), dimension(:,:,:,:), intent(inout) :: tracers
    real(wp), dimension(:,0:,0:,:), intent(inout) :: u

    integer :: i, j, k, l, n, ii, jj, kk, iii, jjj, kkk
    real(wp), dimension(n_tracers_tot) :: tr_before, tr_after, dtr, tr_sum
    real(wp), allocatable, dimension(:,:,:) :: ocn_vol_tmp


    !-------------------------------------------------------
    ! adjust tracer concentrations for new ocean volume (conserve tracer inventories)
    !-------------------------------------------------------

    allocate(ocn_vol_tmp(maxi,maxj,maxk))
    ocn_vol_tmp = grid%ocn_vol_old

    ! total tracer amount before grid update
    do l=2,n_tracers_tot    ! exclude temperature (index==1)
      tr_before(l) = sum(tracers(:,:,:,l) * grid%ocn_vol_old)
    enddo

    ! initialize tracer concentration in newly formed grid cells 
    ! using concentration in neighbouring cells 
    do i=1,maxi
      do j=1,maxj
        do k=maxk,grid%k1(i,j),-1
          if (grid%ocn_vol_old(i,j,k).eq.0._wp .and. grid%ocn_vol(i,j,k).gt.0._wp) then
            ! new cell, initialize 
            n = 0
            tr_sum = 0._wp
            do ii=i-1,i+1
              do jj=j-1,j+1
                do kk=k-1,k+1
                  iii = ii
                  if (iii.eq.0) iii = maxi
                  if (iii.eq.maxi+1) iii = 1
                  jjj = jj
                  jjj = max(1,jjj)
                  jjj = min(maxj,jjj)
                  kkk = kk
                  kkk = max(1,kkk)
                  kkk = min(maxk,kkk)
                  if (ocn_vol_tmp(iii,jjj,kkk).gt.0._wp) then
                    n = n+1
                    do l=1,n_tracers_tot
                      tr_sum(l) = tr_sum(l) + tracers(iii,jjj,kkk,l)
                    enddo
                  endif
                enddo
              enddo
            enddo
            if (n.eq.0) then
              print *,'no ocean neighbours to initialize tracer concentrations in cell (i,j,k)',i,j,k,', whole ocean level average used instead!'
              n = count(grid%ocn_vol(:,:,k).gt.0._wp)
              do l=1,n_tracers_tot
                tr_sum(l) = sum(tracers(:,:,k,l),mask=grid%ocn_vol(:,:,k).gt.0._wp)
                print *,l,tr_sum(l)/real(n,wp),n
              enddo
            endif
            do l=1,n_tracers_tot
              tracers(i,j,k,l) = tr_sum(l)/real(n,wp)
            enddo
            ! update temporary volume with filled cells
            ocn_vol_tmp(i,j,k) = grid%ocn_vol(i,j,k)
          endif
        enddo
      enddo
    enddo

    ! set tracer concentration to zero outside domain, i.e. in dead cells
    do l=1,n_tracers_tot
      where (grid%mask_c.eq.0)
        tracers(:,:,:,l) = 0._wp
      endwhere
    enddo

    ! total tracer amount after grid update
    do l=2,n_tracers_tot
      tr_after(l) = sum(tracers(:,:,:,l) * grid%ocn_vol)
    enddo

    ! change in global ocean concentration needed to conserve tracers
    do l=2,n_tracers_tot
      dtr(l) = (tr_after(l)-tr_before(l))/grid%ocn_vol_tot
      if (l.eq.i_age .or. l.eq.i_dye) then
        where (grid%mask_c.eq.1)
          tracers(:,:,:,l) = tracers(:,:,:,l) - dtr(l)
        endwhere
      else
        where (grid%mask_c.eq.1)
          tracers(:,:,:,l) = tracers(:,:,:,l) * (1._wp - dtr(l)/max(1.e-30_wp,tr_before(l)))
        endwhere
      endif
    enddo


    ! new total tracer amount after grid update, just to check conservation
    do l=2,n_tracers_tot
      tr_after(l) = sum(tracers(:,:,:,l) * grid%ocn_vol)
    enddo

    !print *,'ocn_vol_tot',grid%ocn_vol_tot
    !do l=2,n_tracers_tot
    !  print *,l,tr_before(l),tr_after(l),tr_after(l)-tr_before(l)
    !enddo

    deallocate(ocn_vol_tmp)

    !-------------------------------------------------------
    ! set velocities to zero outside of domain 
    !-------------------------------------------------------

    do i=1,maxi
      do j=1,maxj
        do k=1,maxk
          if (mask_u(i,j,k).eq.0) then
            u(1,i,j,k) = 0._wp
          endif
          if (mask_v(i,j,k).eq.0) then
            u(2,i,j,k) = 0._wp
          endif
          if (mask_w(i,j,k).eq.0) then
            u(3,i,j,k) = 0._wp
          endif
        enddo
      enddo
    enddo

    return
    
  end subroutine ocn_grid_update_state

end module ocn_grid_update_state_mod

