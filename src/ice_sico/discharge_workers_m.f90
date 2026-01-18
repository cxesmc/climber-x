!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  d i s c h a r g e _ w o r k e r s _ m
!
!> @file
!!
!! Ice discharge parameterization for the Greenland ice sheet.
!!
!! References:
!! @li Calov, R. and A. Robinson,  and M. Perrette, \n
!!     and A. Ganopolski 2015.\n
!!     Simulating the Greenland ice sheet under present-day and 
!!     palaeo constraints including a new discharge parameterization.\n
!!     The Cryosphere 9: 179-196.
!!
!! @section Copyright
!!
!! Copyright 2011-2018 Reinhard Calov, Andrey Ganopolski
!! Potsdam Institute for Climate Impact Research  
!!
!! @section License
!!
!! This file is part of SICOPOLIS.
!!
!! SICOPOLIS is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! SICOPOLIS is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with SICOPOLIS.  If not, see <http://www.gnu.org/licenses/>.
!<
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------------------------------------------------------------
!> Ice discharge parameterization for the Greenland ice sheet.
!<------------------------------------------------------------------------------
module discharge_workers_m

  use constants, only : pi
  use sico_types_m
  use sico_grid_mod, only : sico_grid_class
  use sico_state, only : sico_state_class
  use sico_params, only : sico_par_class
  use sico_timer, only : sico_timer_class
  
  use compare_float_m

  implicit none

  private

  public :: discharge

contains


!-------------------------------------------------------------------------------
!> Ice discharge parameterization main formula, controler (general).
!! [Compute ice discharge via a parameterization using distance of ice 
!! margin to coast and ice thickness as parameters.]
!<------------------------------------------------------------------------------
  subroutine discharge(st, grd, tmr, par)

  ! Authors: Reinhard Calov, Andrey Ganopolski
  ! Institution: Potsdam Institute for Climate Impact Research  
  ! Date: 11.6.16

  implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class), intent(in) :: grd
  type(sico_timer_class), intent(inout) :: tmr
  type(sico_par_class), intent(in) :: par

  real(wp) :: cos_tc

  integer :: i, j, ic, jc
  real(wp) :: t_freeze
  real(wp) :: f_sd, f_l, f_r
  real(wp), parameter :: dl = 32._wp  ! km


  cos_tc = cos(par%alpha_tc*pi/180._wp)

  tmr%n_discharge_call = tmr%n_discharge_call + 1

  if ((mod(tmr%n_discharge_call, tmr%iter_mar_coa)==0).or.(tmr%n_discharge_call==1)) then
     ! define marginal ring
     call marginal_ring(st, grd, tmr, par)
     ! compute distance to coast and save indexes of nearest coastal point
     call coastal_distance(st, grd, tmr, par)
  end if

  do j=0, grd%JMAX
    do i=0, grd%IMAX
      if (st%mask_mar(j,i).eq.1.and.st%cos_grad_tc(j,i).ge.cos_tc) then
        ! indexes of nearest coastal point
        ic = st%i_cst(j,i)
        jc = st%j_cst(j,i)
        ! freezing temperature after Beckmann & Goose 2003 (eq. 2), no depth dependence
        t_freeze = 0.0939_wp - 0.057_wp*st%s_ocn(jc,ic)   ! degC
        if (par%ii_disc.eq.1) then
          st%dis_perp(j,i) = tanh(max(0._wp,st%zl_std(j,i))/par%zl_std_crit_disc) * grd%DX**par%m_R/st%cst_dist(j,i)**par%m_D &
            * (par%c_dis_bm*max(0._wp,st%t_ocn(jc,ic)-t_freeze)**par%m_T + par%c_dis_clv*(st%H_c(j,i)+st%H_t(j,i))**par%m_H)
        else if (par%ii_disc.eq.2) then
          st%dis_perp(j,i) = tanh(max(0._wp,st%zl_std(j,i))/par%zl_std_crit_disc) * (grd%DX/32._wp)**par%m_R &
            * (32._wp*1000._wp/st%cst_dist(j,i))**par%m_D * (st%H_c(j,i)+st%H_t(j,i))**par%m_H &
            * (par%c_dis_clv + par%c_dis_bm*max(0._wp,st%t_ocn(jc,ic)-t_freeze)**par%m_T) / (86400._wp*365._wp)
        else if (par%ii_disc.eq.3) then
          ! Calculate scaling factors (roughness, distance to coast and resolution)
          f_sd = tanh(st%zl_std(j,i)/par%zl_std_crit_disc)
          f_l  = ( dl / (dl + max(dl,st%cst_dist(j,i)/1000._wp)) )**par%m_D
          f_r  = ( grd%DX/dl )**par%m_R 
          ! Calculate subgrid discharge mass balance rate 
          st%dis_perp(j,i) =  f_sd * f_l * f_r * ((st%H_c(j,i)+st%H_t(j,i)) / par%tau_dmb) 
        endif
      else
        st%dis_perp(j,i) = 0.0_wp
      endif
    enddo
  enddo
  st%dis_perp = max(0.0_wp, st%dis_perp) ! ensure positive values


  end subroutine discharge


!-------------------------------------------------------------------------------
!> Distance to the coast (general).
!! [Compute distance to the coast for every land point.]
!<------------------------------------------------------------------------------
  subroutine coastal_distance(st, grd, tmr, par)

  ! Author: Reinhard Calov
  ! Institution: Potsdam Institute for Climate Impact Research  
  ! Date: 2.11.11

  implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class), intent(in) :: grd
  type(sico_timer_class), intent(in) :: tmr
  type(sico_par_class), intent(in) :: par

  ! --------------------------------------------------
  real(wp), parameter :: c_smooth=0.2_wp
  ! --------------------------------------------------
  integer :: i, j, i_pos, j_pos, l, l_e, d_l
  real(wp) :: cst_dist_tmp, qrt_1, qrt_2
  logical :: leave_loop

  real(wp), allocatable, dimension(:,:) :: zs_tmp, grad_zs_x, grad_zs_y, &
  grad_dist_x, grad_dist_y

  allocate(zs_tmp(0:grd%JMAX,0:grd%IMAX), grad_zs_x(0:grd%JMAX,0:grd%IMAX), grad_zs_y(0:grd%JMAX,0:grd%IMAX), &
  grad_dist_x(0:grd%JMAX,0:grd%IMAX), grad_dist_y(0:grd%JMAX,0:grd%IMAX))

    select case (par%i_dist_coast)

    case(1)
    ! brute force computation of minimal distance to coast for every land point

      do i_pos=0, grd%IMAX
      do j_pos=0, grd%JMAX
        if(st%maske(j_pos,i_pos).ne.2) then
          st%cst_dist(j_pos,i_pos)=1.d+20
          do i=0, grd%IMAX
          do j=0, grd%JMAX
            if(st%maske(j,i).eq.2) then
              cst_dist_tmp=sqrt((grd%xi(i_pos)-grd%xi(i))**2+(grd%eta(j_pos)-grd%eta(j))**2)
              if(cst_dist_tmp.le.st%cst_dist(j_pos,i_pos)) then
                st%cst_dist(j_pos,i_pos)=cst_dist_tmp
                st%i_cst(j_pos,i_pos) = i
                st%j_cst(j_pos,i_pos) = j
              end if
            end if
          end do
          end do
        else
          st%cst_dist(j_pos,i_pos)=par%dist_coast_min
          st%i_cst(j_pos,i_pos) = i_pos
          st%j_cst(j_pos,i_pos) = j_pos
        end if
      end do
      end do

    case(2)
      ! 2 : Defining a square with two grid distances side length around the 
      ! actual grid point and enlanging the square successively by
      ! one grid. Because this should have been a circle, the
      ! square is enlarged again if the coast line was found; this time 
      ! enlargement is by a rectangles a the respective four side of 
      ! the square. Smaller side lenght equals the rectangle root of 2 of the
      ! larger side. 

      do i_pos=0, grd%IMAX
      do j_pos=0, grd%JMAX
        if(st%maske(j_pos,i_pos).ne.2) then
          leave_loop=.false.
          st%cst_dist(j_pos,i_pos)=1.e20_wp

          do l=1, max(grd%IMAX,grd%JMAX)

            do i=max(i_pos-l,0), min(i_pos+l,grd%IMAX)
              j=max(j_pos-l,0); j=min(j,grd%JMAX)
              if(st%maske(j,i).eq.2) then
                leave_loop=.true.
                cst_dist_tmp=sqrt((grd%xi(i_pos)-grd%xi(i))**2+(grd%eta(j_pos)-grd%eta(j))**2)
                if(cst_dist_tmp.le.st%cst_dist(j_pos,i_pos)) then
                  st%cst_dist(j_pos,i_pos)=cst_dist_tmp
                  st%i_cst(j_pos,i_pos) = i
                  st%j_cst(j_pos,i_pos) = j
                end if
              end if
              j=min(j_pos+l,grd%JMAX)
              if(st%maske(j,i).eq.2) then
                leave_loop=.true.
                cst_dist_tmp=sqrt((grd%xi(i_pos)-grd%xi(i))**2+(grd%eta(j_pos)-grd%eta(j))**2)
                if(cst_dist_tmp.le.st%cst_dist(j_pos,i_pos)) then
                  st%cst_dist(j_pos,i_pos)=cst_dist_tmp
                  st%i_cst(j_pos,i_pos) = i
                  st%j_cst(j_pos,i_pos) = j
                end if
              end if
            end do
            do j=max(j_pos-l+1,0), min(j_pos+l-1,grd%JMAX)
              i=max(i_pos-l,0); i=min(i,grd%IMAX)
              if(st%maske(j,i).eq.2) then
                leave_loop=.true.
                cst_dist_tmp=sqrt((grd%xi(i_pos)-grd%xi(i))**2+(grd%eta(j_pos)-grd%eta(j))**2)
                if(cst_dist_tmp.le.st%cst_dist(j_pos,i_pos)) then
                  st%cst_dist(j_pos,i_pos)=cst_dist_tmp
                  st%i_cst(j_pos,i_pos) = i
                  st%j_cst(j_pos,i_pos) = j
                end if
              end if
              i=min(i_pos+l,grd%IMAX)
              if(st%maske(j,i).eq.2) then
                leave_loop=.true.
                cst_dist_tmp=sqrt((grd%xi(i_pos)-grd%xi(i))**2+(grd%eta(j_pos)-grd%eta(j))**2)
                if(cst_dist_tmp.le.st%cst_dist(j_pos,i_pos)) then
                  st%cst_dist(j_pos,i_pos)=cst_dist_tmp
                  st%i_cst(j_pos,i_pos) = i
                  st%j_cst(j_pos,i_pos) = j
                end if
              end if
            end do

            if(leave_loop) then
              l_e=l
              d_l=ceiling(l_e*sqrt(2.0_wp)+0.001_wp)
              exit   ! leave loop in l
            end if

          end do
          ! left
          do i=max(i_pos-(l_e+d_l),0),min(i_pos-(l_e-1),grd%IMAX)
          do j=max(j_pos-l_e,0),min(j_pos+l_e,grd%JMAX)
            if(st%maske(j,i).eq.2) then
              cst_dist_tmp=sqrt((grd%xi(i_pos)-grd%xi(i))**2+(grd%eta(j_pos)-grd%eta(j))**2)
              if(cst_dist_tmp.le.st%cst_dist(j_pos,i_pos)) then
                st%cst_dist(j_pos,i_pos)=cst_dist_tmp
                st%i_cst(j_pos,i_pos) = i
                st%j_cst(j_pos,i_pos) = j
              end if
            end if
          end do
          end do
          ! right
          do i=max(i_pos+l_e+1,0),min(i_pos+l_e+d_l,grd%IMAX)
          do j=max(j_pos-l_e,0),min(j_pos+l_e,grd%JMAX)
            if(st%maske(j,i).eq.2) then
              cst_dist_tmp=sqrt((grd%xi(i_pos)-grd%xi(i))**2+(grd%eta(j_pos)-grd%eta(j))**2)
              if(cst_dist_tmp.le.st%cst_dist(j_pos,i_pos)) then
                st%cst_dist(j_pos,i_pos)=cst_dist_tmp
                st%i_cst(j_pos,i_pos) = i
                st%j_cst(j_pos,i_pos) = j
              end if
            end if
          end do
          end do
          ! lower
          do i=max(i_pos-l_e,0),min(i_pos+l_e,grd%IMAX)
          do j=max(j_pos-(l_e+d_l),0),min(j_pos-(l_e-1),grd%JMAX)
            if(st%maske(j,i).eq.2) then
              cst_dist_tmp=sqrt((grd%xi(i_pos)-grd%xi(i))**2+(grd%eta(j_pos)-grd%eta(j))**2)
              if(cst_dist_tmp.le.st%cst_dist(j_pos,i_pos)) then
                st%cst_dist(j_pos,i_pos)=cst_dist_tmp
                st%i_cst(j_pos,i_pos) = i
                st%j_cst(j_pos,i_pos) = j
              end if
            end if
          end do
          end do
          ! upper
          do i=max(i_pos-l_e,0),min(i_pos+l_e,grd%IMAX)
          do j=max(j_pos+l_e+1,0),min(j_pos+l_e+d_l,grd%JMAX)
            if(st%maske(j,i).eq.2) then
              cst_dist_tmp=sqrt((grd%xi(i_pos)-grd%xi(i))**2+(grd%eta(j_pos)-grd%eta(j))**2)
              if(cst_dist_tmp.le.st%cst_dist(j_pos,i_pos)) then
                st%cst_dist(j_pos,i_pos)=cst_dist_tmp
                st%i_cst(j_pos,i_pos) = i
                st%j_cst(j_pos,i_pos) = j
              end if
            end if
          end do
          end do
        else
          st%cst_dist(j_pos,i_pos)=par%dist_coast_min
          st%i_cst(j_pos,i_pos) = i_pos
          st%j_cst(j_pos,i_pos) = j_pos
        end if
      end do
      end do

    case(3)
      ! 2 : Defining a square with two grid distances side length around the 
      ! actual grid point and enlanging the square successively by
      ! one grid. Because this should have been a circle, the
      ! square is enlarged again if the coast line was found; this time 
      ! enlargement is by a rectangles a the respective four side of 
      ! the square. Smaller side lenght equals the rectangle root of 2 of the
      ! larger side. 

      do i_pos=0, grd%IMAX
      do j_pos=0, grd%JMAX
        if(st%maske(j_pos,i_pos).ne.2 .and. st%maske(j_pos,i_pos).ne.3) then
          leave_loop=.false.
          st%cst_dist(j_pos,i_pos)=1.e20_wp

          do l=1, max(grd%IMAX,grd%JMAX)

            do i=max(i_pos-l,0), min(i_pos+l,grd%IMAX)
              j=max(j_pos-l,0); j=min(j,grd%JMAX)
              if(st%maske(j,i).eq.2 .or. st%maske(j,i).eq.3) then
                leave_loop=.true.
                cst_dist_tmp=sqrt((grd%xi(i_pos)-grd%xi(i))**2+(grd%eta(j_pos)-grd%eta(j))**2)
                if(cst_dist_tmp.le.st%cst_dist(j_pos,i_pos)) then
                  st%cst_dist(j_pos,i_pos)=cst_dist_tmp
                  st%i_cst(j_pos,i_pos) = i
                  st%j_cst(j_pos,i_pos) = j
                end if
              end if
              j=min(j_pos+l,grd%JMAX)
              if(st%maske(j,i).eq.2 .or. st%maske(j,i).eq.3) then
                leave_loop=.true.
                cst_dist_tmp=sqrt((grd%xi(i_pos)-grd%xi(i))**2+(grd%eta(j_pos)-grd%eta(j))**2)
                if(cst_dist_tmp.le.st%cst_dist(j_pos,i_pos)) then
                  st%cst_dist(j_pos,i_pos)=cst_dist_tmp
                  st%i_cst(j_pos,i_pos) = i
                  st%j_cst(j_pos,i_pos) = j
                end if
              end if
            end do
            do j=max(j_pos-l+1,0), min(j_pos+l-1,grd%JMAX)
              i=max(i_pos-l,0); i=min(i,grd%IMAX)
              if(st%maske(j,i).eq.2 .or. st%maske(j,i).eq.3) then
                leave_loop=.true.
                cst_dist_tmp=sqrt((grd%xi(i_pos)-grd%xi(i))**2+(grd%eta(j_pos)-grd%eta(j))**2)
                if(cst_dist_tmp.le.st%cst_dist(j_pos,i_pos)) then
                  st%cst_dist(j_pos,i_pos)=cst_dist_tmp
                  st%i_cst(j_pos,i_pos) = i
                  st%j_cst(j_pos,i_pos) = j
                end if
              end if
              i=min(i_pos+l,grd%IMAX)
              if(st%maske(j,i).eq.2 .or. st%maske(j,i).eq.3) then
                leave_loop=.true.
                cst_dist_tmp=sqrt((grd%xi(i_pos)-grd%xi(i))**2+(grd%eta(j_pos)-grd%eta(j))**2)
                if(cst_dist_tmp.le.st%cst_dist(j_pos,i_pos)) then
                  st%cst_dist(j_pos,i_pos)=cst_dist_tmp
                  st%i_cst(j_pos,i_pos) = i
                  st%j_cst(j_pos,i_pos) = j
                end if
              end if
            end do

            if(leave_loop) then
              l_e=l
              d_l=ceiling(l_e*sqrt(2.0_wp)+0.001_wp)
              exit   ! leave loop in l
            end if

          end do
          ! left
          do i=max(i_pos-(l_e+d_l),0),min(i_pos-(l_e-1),grd%IMAX)
          do j=max(j_pos-l_e,0),min(j_pos+l_e,grd%JMAX)
            if(st%maske(j,i).eq.2 .or. st%maske(j,i).eq.3) then
              cst_dist_tmp=sqrt((grd%xi(i_pos)-grd%xi(i))**2+(grd%eta(j_pos)-grd%eta(j))**2)
              if(cst_dist_tmp.le.st%cst_dist(j_pos,i_pos)) then
                st%cst_dist(j_pos,i_pos)=cst_dist_tmp
                st%i_cst(j_pos,i_pos) = i
                st%j_cst(j_pos,i_pos) = j
              end if
            end if
          end do
          end do
          ! right
          do i=max(i_pos+l_e+1,0),min(i_pos+l_e+d_l,grd%IMAX)
          do j=max(j_pos-l_e,0),min(j_pos+l_e,grd%JMAX)
            if(st%maske(j,i).eq.2 .or. st%maske(j,i).eq.3) then
              cst_dist_tmp=sqrt((grd%xi(i_pos)-grd%xi(i))**2+(grd%eta(j_pos)-grd%eta(j))**2)
              if(cst_dist_tmp.le.st%cst_dist(j_pos,i_pos)) then
                st%cst_dist(j_pos,i_pos)=cst_dist_tmp
                st%i_cst(j_pos,i_pos) = i
                st%j_cst(j_pos,i_pos) = j
              end if
            end if
          end do
          end do
          ! lower
          do i=max(i_pos-l_e,0),min(i_pos+l_e,grd%IMAX)
          do j=max(j_pos-(l_e+d_l),0),min(j_pos-(l_e-1),grd%JMAX)
            if(st%maske(j,i).eq.2 .or. st%maske(j,i).eq.3) then
              cst_dist_tmp=sqrt((grd%xi(i_pos)-grd%xi(i))**2+(grd%eta(j_pos)-grd%eta(j))**2)
              if(cst_dist_tmp.le.st%cst_dist(j_pos,i_pos)) then
                st%cst_dist(j_pos,i_pos)=cst_dist_tmp
                st%i_cst(j_pos,i_pos) = i
                st%j_cst(j_pos,i_pos) = j
              end if
            end if
          end do
          end do
          ! upper
          do i=max(i_pos-l_e,0),min(i_pos+l_e,grd%IMAX)
          do j=max(j_pos+l_e+1,0),min(j_pos+l_e+d_l,grd%JMAX)
            if(st%maske(j,i).eq.2 .or. st%maske(j,i).eq.3) then
              cst_dist_tmp=sqrt((grd%xi(i_pos)-grd%xi(i))**2+(grd%eta(j_pos)-grd%eta(j))**2)
              if(cst_dist_tmp.le.st%cst_dist(j_pos,i_pos)) then
                st%cst_dist(j_pos,i_pos)=cst_dist_tmp
                st%i_cst(j_pos,i_pos) = i
                st%j_cst(j_pos,i_pos) = j
              end if
            end if
          end do
          end do
        else
          st%cst_dist(j_pos,i_pos)=par%dist_coast_min
          st%i_cst(j_pos,i_pos) = i_pos
          st%j_cst(j_pos,i_pos) = j_pos
        end if
      end do
      end do

      ! make sure distance is 0 for grounded ice points at the grounding line,
      ! to be consistent with Yelmo
      st%cst_dist = st%cst_dist-grd%DX*1000._wp

    end select

  ! set minimum distance
    st%cst_dist = max(par%dist_coast_min, st%cst_dist)

  ! forbit too large values
    st%cst_dist = min(1.0e10_wp, st%cst_dist)

    zs_tmp = st%zs

  ! smoothing the surface topography
    do i=1, grd%IMAX-1
    do j=1, grd%JMAX-1
       zs_tmp(j,i) = c_smooth*(st%zs(j,i-1)+st%zs(j,i+1)+st%zs(j-1,i)+st%zs(j+1,i)) &
                     +(1.0_wp-4.0_wp*c_smooth)*st%zs(j,i)
    end do
    end do

  ! gradients (2. order for now ...)

!  ------ x-derivatives

    do i=1, grd%IMAX-1
    do j=0, grd%JMAX
      grad_zs_x(j,i)    = (zs_tmp(j,i+1)-zs_tmp(j,i-1))*0.5_wp*grd%dxi_inv*grd%insq_g11_g(j,i)
      grad_dist_x(j,i)  = (st%cst_dist(j,i+1)-st%cst_dist(j,i-1))*0.5_wp*grd%dxi_inv*grd%insq_g11_g(j,i)
    end do
    end do

    do j=0, grd%JMAX
      grad_zs_x(j,0)       = (zs_tmp(j,1)-zs_tmp(j,0))*grd%dxi_inv*grd%insq_g11_g(j,0)
      grad_zs_x(j,grd%IMAX)    = (zs_tmp(j,grd%IMAX)-zs_tmp(j,grd%IMAX-1))*grd%dxi_inv*grd%insq_g11_g(j,grd%IMAX)
      grad_dist_x(j,0)     = (st%cst_dist(j,1)-st%cst_dist(j,0))*grd%dxi_inv*grd%insq_g11_g(j,0)
      grad_dist_x(j,grd%IMAX)  = (st%cst_dist(j,grd%IMAX)-st%cst_dist(j,grd%IMAX-1))*grd%dxi_inv*grd%insq_g11_g(j,grd%IMAX)
    end do

!  ------ y-derivatives

    do i=0, grd%IMAX
    do j=1, grd%JMAX-1
      grad_zs_y(j,i)    = (zs_tmp(j+1,i)-zs_tmp(j-1,i))*0.5_wp*grd%deta_inv*grd%insq_g22_g(j,i)
      grad_dist_y(j,i)  = (st%cst_dist(j+1,i)-st%cst_dist(j-1,i))*0.5_wp*grd%deta_inv*grd%insq_g22_g(j,i)
    end do
    end do

    do i=0, grd%IMAX
      grad_zs_y(0,i)       = (zs_tmp(1,i)-zs_tmp(0,i))*grd%deta_inv*grd%insq_g22_g(0,i)
      grad_zs_y(grd%JMAX,i)    = (zs_tmp(grd%JMAX,i)-zs_tmp(grd%JMAX-1,i))*grd%deta_inv*grd%insq_g22_g(grd%JMAX,i)
      grad_dist_y(0,i)     = (st%cst_dist(1,i)-st%cst_dist(0,i))*grd%deta_inv*grd%insq_g22_g(0,i)
      grad_dist_y(grd%JMAX,i)  = (st%cst_dist(grd%JMAX,i)-st%cst_dist(grd%JMAX-1,i))*grd%deta_inv*grd%insq_g22_g(grd%JMAX,i)
    end do

  ! angle between topographical gradient and coastal distance gradient

    do i=0, grd%IMAX
    do j=0, grd%JMAX
      qrt_1=grad_zs_x(j,i)*grad_zs_x(j,i)+grad_zs_y(j,i)*grad_zs_y(j,i)
      qrt_2=grad_dist_x(j,i)*grad_dist_x(j,i)+grad_dist_y(j,i)*grad_dist_y(j,i)
      if(qrt_1.ne.0.0_wp.and.qrt_2.ne.0.0_wp) then
  !!! top_cst_alpha(j,i)=dacos( (grad_zs_x(j,i)*grad_dist_x(j,i)+grad_zs_y(j,i)*grad_dist_y(j,i))/ &
  !!! (sqrt(qrt_1)*sqrt(qrt_2)) )

        st%cos_grad_tc(j,i)= (grad_zs_x(j,i)*grad_dist_x(j,i)+grad_zs_y(j,i)*grad_dist_y(j,i))/ &
        (sqrt(qrt_1)*sqrt(qrt_2))
      else
        st%cos_grad_tc(j,i)=-1.0_wp
      end if
    end do
    end do

    deallocate(zs_tmp, grad_zs_x, grad_zs_y, grad_dist_x, grad_dist_y)

  end subroutine coastal_distance

!-------------------------------------------------------------------------------
!> Ring along an ice sheet margin (general).
!! [Compute marginal ring.]
!<------------------------------------------------------------------------------
  subroutine marginal_ring(st, grd, tmr, par)

  ! Author: Reinhard Calov
  ! Institution: Potsdam Institute for Climate Impact Research  
  ! Date: 28.10.11

  ! Purpose: Determine an r_max_eff wide ring along ice margins 
  ! towards the interior of the ice sheet. Small stripes of land are
  ! not considered as land. For the logics De Morgan is used often.

  ! Methode: Two staggered loops in i,j. The inner i,j loop acts inside
  ! a rectangle defined by r_mar_eff. r_mar_eff should be small compared
  ! to the domain size. 

    implicit none

    type(sico_state_class), intent(inout) :: st
    type(sico_grid_class), intent(in) :: grd
    type(sico_timer_class), intent(in) :: tmr
    type(sico_par_class), intent(in) :: par

    integer :: i, j, i_pos, j_pos, di_eff, dj_eff
    integer :: count_tmp
    real(wp) :: r_p

    if(par%r_mar_eff.le.1.0e6_wp) then

      st%mask_mar=0
      do i_pos=1, grd%IMAX-1
      do j_pos=1, grd%JMAX-1
        if((st%maske(j_pos,i_pos).eq.1.or.st%maske(j_pos,i_pos).eq.2).and.     &
           .not.(st%maske(j_pos,i_pos).eq.1.and.st%maske(j_pos,i_pos-1).eq.0   &
                 .and.st%maske(j_pos,i_pos+1).eq.0.or. &
                  st%maske(j_pos,i_pos).eq.1.and.st%maske(j_pos-1,i_pos).eq.0  &
                  .and.st%maske(j_pos+1,i_pos).eq.0.or. &
                   st%maske(j_pos,i_pos).eq.1.and.st%maske(j_pos,i_pos-1).eq.3 &
                   .and.st%maske(j_pos,i_pos+1).eq.3.or. &
                   st%maske(j_pos,i_pos).eq.1.and.st%maske(j_pos-1,i_pos).eq.3 &
                   .and.st%maske(j_pos+1,i_pos).eq.3.or. &
                   st%maske(j_pos,i_pos).eq.1.and.st%maske(j_pos,i_pos-1).eq.3 &
                   .and.st%maske(j_pos,i_pos+1).eq.0.or. &
                   st%maske(j_pos,i_pos).eq.1.and.st%maske(j_pos-1,i_pos).eq.0 &
                   .and.st%maske(j_pos+1,i_pos).eq.3)) then ! outside ice sheet, exclude isolated land stripes

          di_eff=int(par%r_mar_eff/grd%dxi)+1; dj_eff=int(par%r_mar_eff/grd%deta)+1 ! only for grid=0, 1 yet!

          do i=max(i_pos-di_eff,0), min(i_pos+di_eff,grd%IMAX)
          do j=max(j_pos-dj_eff,0), min(j_pos+dj_eff,grd%JMAX)
            r_p=sqrt((grd%xi(i_pos)-grd%xi(i))**2+(grd%eta(j_pos)-grd%eta(j))**2)

            if(r_p.le.par%r_mar_eff.and.(st%maske(j,i).eq.0.or.st%maske(j,i).eq.3)) then

              st%mask_mar(j,i)=1
            end if
          end do
          end do
        end if
      end do
      end do

      where(st%maske.ge.1)
        st%mask_mar=1
      end where

    else ! the ring encompassed entire Greenland 
      st%mask_mar=1
    end if

  end subroutine marginal_ring

end module discharge_workers_m
!
