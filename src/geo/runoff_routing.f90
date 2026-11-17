!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : r u n o f f _ r o u t i n g _ m o d
!
!  Purpose : runoff routing to ocean
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
module runoff_routing_mod

  use precision, only : wp 
  use lake_mod, only : lake_type
  !$ use omp_lib

  implicit none

  private
  public :: runoff_routing, runoff_routing_lakes

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  r u n o f f _ r o u t i n g
  !   Purpose    :  runoff routing of filled topography
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine runoff_routing(z_topo_fill, mask, area, &
      i_runoff, j_runoff, flow_acc)

    implicit none

    real(wp), intent(in) :: z_topo_fill(:,:)
    integer, intent(in) :: mask(:,:)
    real(wp), intent(in) :: area(:,:)
    integer, intent(out) :: i_runoff(:,:), j_runoff(:,:)
    real(wp), intent(out) :: flow_acc(:,:)

    integer :: i, j, n, ir, jr, loop, ni, nj
    integer :: im1, ip1, jm1, jp1
    integer, parameter :: ndir = 8 !! number of runoff directions
    integer, parameter :: irn=1, irne=2, ire=3, irse=4, irs=5, irsw=6, irw=7, irnw=8
    real(wp) :: z_topo_ij
    real(wp), dimension(:,:), allocatable :: z_topo_fill_ext
    real(wp), dimension(ndir) :: slope
    integer, dimension(:,:), allocatable :: rundir


    ! get grid size
    ni = size(z_topo_fill,1)
    nj = size(z_topo_fill,2)

    ! allocate
    allocate(z_topo_fill_ext(0:ni+1,0:nj+1))
    allocate(rundir(ni,nj))

    ! extend z_topo_fill

    z_topo_fill_ext(1:ni,1:nj) = z_topo_fill(1:ni,1:nj)
    z_topo_fill_ext(:,0) = 9999._wp ! south pole
    z_topo_fill_ext(:,nj+1) = 9999._wp ! north pole
    z_topo_fill_ext(0,:) = z_topo_fill_ext(ni,:)
    z_topo_fill_ext(ni+1,:) = z_topo_fill_ext(1,:)

    !--------------------------------------------------------------
    ! determine runoff directions of filled z_topography
    !--------------------------------------------------------------

    !$omp parallel do private(i,j,n,im1,ip1,jm1,jp1,z_topo_ij,slope)
    do j=1,nj
       do i=1,ni
         if (mask(i,j).ne.2 .and. mask(i,j).ne.3) then
           z_topo_ij = z_topo_fill_ext(i,j)
           im1 = i-1
           ip1 = i+1
           jm1 = j-1
           jp1 = j+1
           ! compute topography slopes in all directions
           do n=1,ndir
             if (n==1) then
               ! north
               slope(n) = z_topo_fill_ext(i,jp1)-z_topo_ij
             else if (n==2) then
               ! north-east
               slope(n) = z_topo_fill_ext(ip1,jp1)-z_topo_ij
             else if (n==3) then
               ! east
               slope(n) = z_topo_fill_ext(ip1,j)-z_topo_ij
             else if (n==4) then
               ! south-east
               slope(n) = z_topo_fill_ext(ip1,jm1)-z_topo_ij
             else if (n==5) then
               ! south
               slope(n) = z_topo_fill_ext(i,jm1)-z_topo_ij
             else if (n==6) then
               ! south-west
               slope(n) = z_topo_fill_ext(im1,jm1)-z_topo_ij
             else if (n==7) then
               ! west
               slope(n) = z_topo_fill_ext(im1,j)-z_topo_ij
             else if (n==8) then
               ! north-west
               slope(n) = z_topo_fill_ext(im1,jp1)-z_topo_ij
             endif
           enddo
           ! find steepest slope direction
           rundir(i,j) = minloc(slope,1)
         else ! ocean
           rundir(i,j) = 0
         endif
       enddo
     enddo
     !$omp end parallel do

    ! river routing
    ! matrix (i_runoff(i,j),j_runoff(i,j)) defines where to put the runoff from point (i,j) 
    flow_acc = 0._wp
    !$omp parallel do private(i,j,ir,jr,loop)
    do j=1,nj
       do i=1,ni
         if (mask(i,j).ne.2) then
          ir = i
          jr = j
          ! follow slope down to the first ocean point
          loop = 0
          do while(rundir(ir,jr) .gt. 0)
             if (rundir(ir,jr) .eq. irn) then
                jr = jr + 1
             else if (rundir(ir,jr) .eq. irne) then
                ir = ir + 1
                jr = jr + 1
             else if (rundir(ir,jr) .eq. ire) then
                ir = ir + 1
             else if (rundir(ir,jr) .eq. irse) then
                ir = ir + 1
                jr = jr - 1
             else if (rundir(ir,jr) .eq. irs) then
                jr = jr - 1
             else if (rundir(ir,jr) .eq. irsw) then
                ir = ir - 1
                jr = jr - 1
             else if (rundir(ir,jr) .eq. irw) then
                ir = ir - 1
             else if (rundir(ir,jr) .eq. irnw) then
                ir = ir - 1
                jr = jr + 1
             endif

             ! periodic b.c.
             if (ir.eq.ni+1) then
                ir = 1
             else if (ir.eq.0)then
                ir = ni
             endif

             ! avoid inf. loops
             loop = loop + 1
             if (loop.gt.100000) then
                print*,'There is a problem calculating runoff'
                print*,'Located at (',i,',',j,')'
                print*,'ir = ',ir
                print*,'jr = ',jr
                print*,'rundir(ir,jr)',rundir(ir,jr)
                print *,'z_topo 3x3'
                print *,z_topo_fill_ext(ir-1:ir+1,jr+1)
                print *,z_topo_fill_ext(ir-1:ir+1,jr)
                print *,z_topo_fill_ext(ir-1:ir+1,jr-1)
                print *,'rundir 3x3'
                print *,rundir(ir-1:ir+1,jr+1)
                print *,rundir(ir-1:ir+1,jr)
                print *,rundir(ir-1:ir+1,jr-1)
                stop 'problem calculating runoff'
             endif

             ! flow accumulation, upstream drainage area: add the area of the source cell (i,j)
             ! to every cell downstream of it; atomic because paths from different j rows
             ! (threads) can pass through the same cell
             !$omp atomic
             flow_acc(ir,jr) = flow_acc(ir,jr) + area(i,j)

          enddo
          ! save the indexes of the runoff destination point (ocean point)
          i_runoff(i,j) = ir
          j_runoff(i,j) = jr
        else
          ! set to 0 for ocean points
          i_runoff(i,j) = 0
          j_runoff(i,j) = 0
        endif
       enddo
    enddo
    !$omp end parallel do

    deallocate(z_topo_fill_ext)
    deallocate(rundir)


   return

  end subroutine runoff_routing


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  r u n o f f _ r o u t i n g _ l a k e s
  !   Purpose    :  mapping of runoff routing to lakes
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine runoff_routing_lakes(z_topo, z_topo_fill, mask, mask_lake_pot, i_runoff, j_runoff, &
      map_runoff, lake)

    implicit none

    real(wp), intent(in) :: z_topo(:,:)
    real(wp), intent(in) :: z_topo_fill(:,:)
    integer, intent(in) :: mask(:,:)
    integer, intent(in) :: mask_lake_pot(:,:)
    integer, intent(in) :: i_runoff(:,:), j_runoff(:,:)
    integer, intent(out) :: map_runoff(:,:)
    type(lake_type), intent(inout) :: lake(:)

    integer :: i, j, n, ir, jr, loop, ni, nj
    integer :: im1, ip1, jm1, jp1
    integer, parameter :: ndir = 8 !! number of runoff directions
    integer, parameter :: irn=1, irne=2, ire=3, irse=4, irs=5, irsw=6, irw=7, irnw=8
    real(wp) :: z_topo_ij
    real(wp), dimension(:,:), allocatable :: z_topo_fill_ext
    real(wp), dimension(ndir) :: slope
    integer, dimension(:,:), allocatable :: rundir


    !--------------------------------------------------------------
    ! runoff routing from lakes to ocean 
    !--------------------------------------------------------------

    do n=1,size(lake)
      lake(n)%hires%i_runoff = i_runoff(lake(n)%hires%i,lake(n)%hires%j)
      lake(n)%hires%j_runoff = j_runoff(lake(n)%hires%i,lake(n)%hires%j)
    enddo

    !--------------------------------------------------------------
    ! map runoff to lakes
    !--------------------------------------------------------------

    ! get grid size
    ni = size(z_topo_fill,1)
    nj = size(z_topo_fill,2)

    allocate(z_topo_fill_ext(0:ni+1,0:nj+1))
    allocate(rundir(ni,nj))

    ! fill z_topo where no lakes
    where (mask_lake_pot==0)
      z_topo_fill_ext(1:ni,1:nj) = z_topo_fill
    elsewhere
      z_topo_fill_ext(1:ni,1:nj) = z_topo
    endwhere
    z_topo_fill_ext(:,0) = 9999._wp ! south pole
    z_topo_fill_ext(:,nj+1) = 9999._wp ! north pole
    z_topo_fill_ext(0,:) = z_topo_fill_ext(ni,:)
    z_topo_fill_ext(ni+1,:) = z_topo_fill_ext(1,:)

    !$omp parallel do private(i,j,n,im1,ip1,jm1,jp1,z_topo_ij,slope)
    do j=1,nj
       do i=1,ni
         if (mask(i,j).ne.2 .and. mask(i,j).ne.3) then
         z_topo_ij = z_topo_fill_ext(i,j)
         im1 = i-1
         ip1 = i+1
         jm1 = j-1
         jp1 = j+1
         do n=1,ndir
           if (n==1) then
             ! north
             slope(n) = z_topo_fill_ext(i,jp1)-z_topo_ij
           else if (n==2) then
             ! north-east
             slope(n) = z_topo_fill_ext(ip1,jp1)-z_topo_ij
           else if (n==3) then
             ! east
             slope(n) = z_topo_fill_ext(ip1,j)-z_topo_ij
           else if (n==4) then
             ! south-east
             slope(n) = z_topo_fill_ext(ip1,jm1)-z_topo_ij
           else if (n==5) then
             ! south
             slope(n) = z_topo_fill_ext(i,jm1)-z_topo_ij
           else if (n==6) then
             ! south-west
             slope(n) = z_topo_fill_ext(im1,jm1)-z_topo_ij
           else if (n==7) then
             ! west
             slope(n) = z_topo_fill_ext(im1,j)-z_topo_ij
           else if (n==8) then
             ! north-west
             slope(n) = z_topo_fill_ext(im1,jp1)-z_topo_ij
           endif
         enddo
         rundir(i,j) = minloc(slope,1)
         if (minval(slope).ge.0._wp) then
           ! local depression
           rundir(i,j) = -1
         endif
       else ! ocean
         rundir(i,j) = 0
       endif
       enddo
     enddo
     !$omp end parallel do

    ! river routing
    !$omp parallel do private(i,j,ir,jr,loop)
    do j=1,nj
       do i=1,ni
         if (mask(i,j).ne.2) then
          ir = i
          jr = j
          loop = 0
          do while(rundir(ir,jr) .gt. 0)
             if (rundir(ir,jr) .eq. irn) then
                jr = jr + 1
             else if (rundir(ir,jr) .eq. irne) then
                ir = ir + 1
                jr = jr + 1
             else if (rundir(ir,jr) .eq. ire) then
                ir = ir + 1
             else if (rundir(ir,jr) .eq. irse) then
                ir = ir + 1
                jr = jr - 1
             else if (rundir(ir,jr) .eq. irs) then
                jr = jr - 1
             else if (rundir(ir,jr) .eq. irsw) then
                ir = ir - 1
                jr = jr - 1
             else if (rundir(ir,jr) .eq. irw) then
                ir = ir - 1
             else if (rundir(ir,jr) .eq. irnw) then
                ir = ir - 1
                jr = jr + 1
             endif

             ! periodic b.c.
             if (ir.eq.ni+1) then
                ir = 1
             else if (ir.eq.0) then
                ir = ni
             endif

             ! avoid inf. loops
             loop = loop + 1
             if (loop.gt.100000) then
                print*,'There is a problem calculating runoff'
                print*,'Located at (',i,',',j,')'
                print*,'ir = ',ir
                print*,'jr = ',jr
                print*,'rundir(ir,jr',rundir(ir,jr),')'
                print *,rundir(ir-1:ir+1,jr+1)
                print *,rundir(ir-1:ir+1,jr)
                print *,rundir(ir-1:ir+1,jr-1)
                print *,z_topo_fill_ext(ir-1:ir+1,jr+1)
                print *,z_topo_fill_ext(ir-1:ir+1,jr)
                print *,z_topo_fill_ext(ir-1:ir+1,jr-1)
                stop 'problem calculating runoff'
             endif

          enddo

          ! force runoff from lake catchment into a single (deepest) point
          map_runoff(i,j) = mask_lake_pot(ir,jr)

        else

          map_runoff(i,j) = 0

        endif
       enddo
    enddo
    !$omp end parallel do

    deallocate(z_topo_fill_ext)
    deallocate(rundir)


   return

  end subroutine runoff_routing_lakes


end module runoff_routing_mod

