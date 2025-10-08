!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : m o m e n t u m _ m o d
!
!  Purpose : momentum equation 
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
module momentum_mod

  use ncio
  use dim_name

  use precision, only : wp
  use constants, only : omega
  use timer, only : time_soy_ocn, doy, year
  use control, only : out_dir
  use climber_grid, only : lon, latv
  use ocn_params, only : drag, drag_bcl, rtv, rtv3, fcor, fcorv, fcormin, fcormin_ref, dt
  use ocn_params, only : drag_par
  use ocn_params, only : i_eos
  use ocn_grid, only : maxi, maxj, maxk, mpxi, mpxj 
  use ocn_grid, only : maxisles, n_isles, psiles
  use ocn_grid, only : dx, dy, dz, zw, zro, s, sv, cv, k1, h

  use eos_mod, only : eos_tb
  use invert_mod, only : invert
  use ubarsolv_mod, only : ubarsolv
  use jbar_mod, only : jbar
  use island_mod, only : island
  use matinv_mod, only : matinv, matmult
  use wind_mod, only : wind
  use velc_mod, only : velc

  !$  use omp_lib

  implicit none

  real(wp), allocatable :: ubar_wind(:,:)
  real(wp), allocatable :: ubar_jbar(:,:)
  real(wp), allocatable :: ratm(:,:)
  real(wp), allocatable :: gb(:)
  real(wp), allocatable :: bp(:,:,:)
  real(wp), allocatable :: gap(:,:)
  real(wp), allocatable :: psisl(:,:,:)
  real(wp), allocatable :: ubisl(:,:,:,:)
  real(wp), allocatable :: erisl(:,:)
  real(wp), allocatable :: psibc(:)
  real(wp), allocatable :: tmpdrg(:,:)
  real(wp), allocatable :: rho_tb(:,:,:)

  private
  public :: momentum_init, momentum
  public :: bp, ubar_wind, ubar_jbar   ! just for output

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  m o m e n t u m
  !   Purpose    :  solve frictional-geostrophic balance equation
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine momentum(f_ocn,tau,dtau_dz2,dtav_dz2,rho, &
                        ub,ub_isl,u,psi)

    implicit none

    real(wp), dimension(:,:), intent(in) :: f_ocn
    real(wp), dimension(:,:,:), intent(in) :: tau
    real(wp), dimension(:,:,:), intent(in) :: dtau_dz2
    real(wp), dimension(:,:,:), intent(in) :: dtav_dz2
    real(wp), dimension(:,:,:), intent(in) :: rho

    real(wp), dimension(:,0:,0:,:), intent(inout) :: u
    real(wp), dimension(:,0:,0:), intent(inout) :: ub
    real(wp), dimension(:,:,:,:), intent(out) :: ub_isl
    real(wp), dimension(0:,0:), intent(inout) :: psi

    integer :: i, j, k, l, isl

    !$ logical, parameter :: print_omp = .false.
    !$ real(wp) :: time1,time2


    if (time_soy_ocn) then
      ! update to new ocean grid
      call momentum_update_grid(f_ocn,tau)
    endif

    if (i_eos.eq.5) then
      ! compute thermobaric density (Ma 2020, Dukowicz 2001)
      do j=1,maxj
        do i=1,maxi
          do k=k1(i,j),maxk
            call eos_tb(rho(i,j,k),zro(k), &
                        rho_tb(i,j,k))
          enddo
        enddo
      enddo
    endif

    ! wind stress term for barotropic streamfunction
    !$ time1 = omp_get_wtime()
    call wind(tau, &    ! in 
              ubar_wind)       ! out
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'momentum: wind',time2-time1

    ! jbar term for barotropic streamfunction
    !$ time1 = omp_get_wtime()
    if (i_eos.ne.5) then
      call jbar(rho, &    ! in
                bp, ubar_jbar)     ! out
    else if (i_eos.eq.5) then
      call jbar(rho_tb, &    ! in
                bp, ubar_jbar)     ! out
    endif
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'momentum: jbar',time2-time1

    ! convert wind and jbar terms to 1d-vector
    gb(:) = 0._wp
    do j=1,maxj
      do i=1,maxi
        l = i + j*maxi
        gb(l) = ubar_wind(i,j) + ubar_jbar(i,j)
      enddo
    enddo

    !$ time1 = omp_get_wtime()
    call ubarsolv(ratm,gap, &
                  gb, &
                  ub,psi)
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'momentum: ubarsolv',time2-time1

    ! find island path integral due to wind and jbar terms 
    !$ time1 = omp_get_wtime()
    do isl=1,n_isles
       call island(ub,tau,bp, &
                   isl,1, &
                   erisl(isl,n_isles+1))
    enddo
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'momentum: islands',time2-time1

    ! solve system of simultaneous equations
    !$ time1 = omp_get_wtime()
    if (n_isles > 1) then
       call matmult(n_isles,erisl,erisl(:,n_isles+1))
       do isl=1,n_isles
          psibc(isl) = - erisl(isl,n_isles+1)
       enddo
     else if (n_isles==1) then
       psibc(1) = - erisl(1,2)/erisl(1,1)
    endif
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'momentum: matmult',time2-time1

    !$ time1 = omp_get_wtime()
    ! sum up barotropic velocity contribution from all islands
    do j=1,maxj
       do i=0,maxi+1
          do isl=1,n_isles
             ub(1,i,j) = ub(1,i,j) + ubisl(1,i,j,isl)*psibc(isl)
             ub(2,i,j) = ub(2,i,j) + ubisl(2,i,j,isl)*psibc(isl)
          enddo
       enddo
    enddo

    ! for output
    ub_isl = 0._wp
    do isl=1,n_isles
      ub_isl(:,:,:,isl) = ubisl(:,1:maxi,1:maxj,isl)*psibc(isl)
    enddo

    ! update diagnostic barotropic streamfunction
    do j=0,maxj
       do i=0,maxi
          do isl=1,n_isles
             psi(i,j) = psi(i,j) + psisl(i,j,isl)*psibc(isl)
          enddo
       enddo
    enddo
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'momentum: psi',time2-time1

    ! update velocities
    !$ time1 = omp_get_wtime()
    if (i_eos.ne.5) then
      call velc(rho,dtau_dz2,dtav_dz2,ub, &
                u)
    else if (i_eos.eq.5) then
      call velc(rho_tb,dtau_dz2,dtav_dz2,ub, &
                u)
    endif
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'momentum: velc',time2-time1


    return

  end subroutine momentum


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  m o m e n t u m _ u p d a t e _ g r i d
  !   Purpose    :  update grid for momentum 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine momentum_update_grid(f_ocn,tau)

    implicit none

    real(wp), dimension(:,:), intent(in) :: f_ocn
    real(wp), dimension(:,:,:), intent(in) :: tau

    integer :: i, j, l, isl, isol
    integer :: i1, i1p, j1, ii, ip1
    real(wp) :: tmp, min_dep, dep_fac, topo_fac, min_frac


    ! calculate drag at psi points (grid cell corners!)   

    do j=0,maxj
      do i=0,maxi
        tmpdrg(i,j) = drag_par%adrag
      enddo
    enddo

    ! increase drag in shallow water regions depeding on water depth
    do j=0,maxj
      do i=0,maxi
        if (drag_par%drag_topo_n.eq.1) then
          min_dep = min(-zw(min(maxk,k1(i,j))),-zw(min(maxk,k1(i+1,j))),-zw(min(maxk,k1(i,j+1))),-zw(min(maxk,k1(i+1,j+1))))
        else if (drag_par%drag_topo_n.eq.2) then
          min_dep = -zw(min(maxk,k1(i,j)))
          do j1=max(0,j-1),min(maxj+1,j+2)
            do i1=i-1,i+2
              i1p = 1 + mod(maxi + i1-1,maxi)
              min_dep = min(min_dep,-zw(min(maxk,k1(i1p,j1))))
            enddo
          enddo 
        endif
        dep_fac = max(0._wp,(drag_par%z_drag_shallow-min_dep)/drag_par%z_drag_shallow)

        !topo_fac = drag_par%drag_topo_fac*(1._wp+drag_par%drag_topo_scale_eq*cv(j))
        if (latv(j+1).ge.-20._wp .and. latv(j+1).le.20._wp) then
          topo_fac = drag_par%drag_topo_fac*(1._wp+drag_par%drag_topo_scale_eq)
        else
          topo_fac = drag_par%drag_topo_fac
        endif

        ii = i
        if (ii.eq.0) ii=maxi
        ip1 = i+1
        if (ip1.eq.maxi+1) ip1=1
        min_frac = min(f_ocn(ii,max(1,j)),f_ocn(ip1,max(1,j)),f_ocn(ii,min(maxj,j+1)),f_ocn(ip1,min(maxj,j+1)))

        !tmp = drag_par%adrag*(1._wp + min(drag_par%drag_topo_fac,topo_fac*dep_fac) + drag_par%drag_frac_fac*(1._wp-min_frac))
        tmp = drag_par%adrag*(1._wp + topo_fac*dep_fac + drag_par%drag_frac_fac*(1._wp-min_frac))

        if (tmp.gt.tmpdrg(i,j)) then
          tmpdrg(i,j) = tmp
        endif

      enddo
    enddo

    ! interpolate to velocity points
    do j=1,maxj
       do i=1,maxi
          drag(1,i,j) = 0.5_wp*(tmpdrg(i,j) + tmpdrg(i,j-1))    ! u-grid
          drag(2,i,j) = 0.5_wp*(tmpdrg(i,j) + tmpdrg(i-1,j))    ! v-grid
       enddo
    enddo

    ! periodic boundary conditions
    do j=1,maxj
       drag(2,maxi+1,j) = drag(2,1,j)
    enddo


    ! for baroclinic velocity

    do j=0,maxj
      do i=0,maxi
        tmpdrg(i,j) = drag_par%adrag
      enddo
    enddo

    ! interpolate to velocity points
    do j=1,maxj
       do i=1,maxi
          drag_bcl(1,i,j) = 0.5_wp*(tmpdrg(i,j) + tmpdrg(i,j-1))    ! u-grid
          drag_bcl(2,i,j) = 0.5_wp*(tmpdrg(i,j) + tmpdrg(i-1,j))    ! v-grid
       enddo
    enddo

    ! boundary conditions, assuming no flow out of N or S bdy
    do j=1,maxj
       drag_bcl(2,maxi+1,j) = drag_bcl(2,1,j)
    enddo


    ! arrays for efficiency
    do j=1,maxj
      do i=1,maxi
        rtv(i,j) = 1._wp/(fcor(j)**2 + drag_bcl(1,i,j)*drag_bcl(1,i,j))
        rtv3(i,j) = 1._wp/(fcorv(j)**2 + drag_bcl(2,i,j)*drag_bcl(2,i,j))
      enddo
    enddo


    ! for barotropic velocity

    call invert(gap,ratm)

    do isol=1,n_isles

     ! set source term to 1 on the ith island (i+1th landmass) only

       do j=0,maxj
          do i=1,maxi
             l=i + j*maxi
             if(psiles(l).eq.isol + 1)then
                gb(l) = 1._wp
             else
                gb(l) = 0._wp
             endif
          enddo
       enddo
       call ubarsolv(ratm,gap, &
                     gb,ubisl(:,:,:,isol),psisl(:,:,isol))

       ! find island path integral due to unit source on boundary
       do isl=1,n_isles
          call island(ubisl(:,:,:,isol),tau,bp, &  ! tau and bp not used in this call (indj==0)!
                      isl,0, &
                      erisl(isl,isol))  ! erisl in m/s2
       enddo
    enddo

    !print*,'island path integrals due to unit sources'
    !do isl=1,n_isles
    !  do isol=1,n_isles
    !    print *,isl,isol,erisl(isl,isol)
    !  enddo
    !enddo

    ! partially invert inland integral error matrix for psi bc calc
    call matinv(n_isles,erisl)  ! inverted erisl in s2/m?

    !print*,'island path integrals due to unit sources after matinv'
    !do isl=1,n_isles
    !  do isol=1,n_isles
    !    print *,isl,isol,erisl(isl,isol)
    !  enddo
    !enddo


    return

  end subroutine momentum_update_grid


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  m o m e n t u m _ i n i t
  !   Purpose    :  initialize momentum equation variables
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine momentum_init

    implicit none

    integer :: j


    ! allocate variables
    allocate(ubar_wind(maxi,maxj))
    allocate(ubar_jbar(maxi,maxj))
    allocate(ratm(mpxi*mpxj,mpxi+1))
    allocate(gb(mpxi*mpxj))
    allocate(gap(mpxi*mpxj,2*mpxi+3))
    allocate(bp(maxi+1,maxj,maxk))
    allocate(psisl(0:maxi,0:maxj,maxisles))
    allocate(ubisl(2,0:maxi+1,0:maxj,maxisles))
    allocate(erisl(maxisles,maxisles+1))
    allocate(psibc(maxisles))
    allocate(drag(2,maxi+1,maxj))
    allocate(drag_bcl(2,maxi+1,maxj))
    allocate(tmpdrg(0:maxi,0:maxj))
    allocate(fcor(0:maxj+1))
    allocate(fcorv(0:maxj+1))
    allocate(rtv(maxi,maxj))
    allocate(rtv3(maxi,maxj))
    allocate(rho_tb(maxi,maxj,maxk))


    drag_par%adrag = 1._wp/(drag_par%adrag*86400._wp) ! 1/s

    ! Coriolis parameter
    do j=1,maxj
      fcor(j) = sign(max(abs(2._wp*omega*s(j)),fcormin),s(j))
    enddo
    do j=0,maxj
      fcorv(j) = sign(max(abs(2._wp*omega*sv(j)),fcormin),sv(j))
    enddo

    return

  end subroutine momentum_init

end module momentum_mod
 
