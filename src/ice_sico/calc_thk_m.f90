!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ t h k _ m
!
!> @file
!!
!! Computation of the ice thickness.
!!
!! @section Copyright
!!
!! Copyright 2009-2017 Ralf Greve, Reinhard Calov, Tatsuru Sato
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
!> Computation of the ice thickness.
!<------------------------------------------------------------------------------
module calc_thk_m

! Include header for lis solver fortran interface
#include "lisf.h"

  use timer, only : sec_year

  use sico_types_m
  use sico_state
  use sico_grid_mod
  use sico_timer
  use sico_params, only : sico_par_class, eps_dp, eps_H, RHO_SW, pi_180
  use sico_maths_m, only : sor_sprs
  use topograd_m

  implicit none

  logical,                            save :: flag_solver_explicit

  private
  public :: calc_thk_init
  public :: calc_thk_sia_expl, calc_thk_sia_impl_sor, calc_thk_sia_impl_lis
  public :: calc_thk_expl, calc_thk_impl_sor, calc_thk_impl_lis
  public :: calc_thk_mask_update
  public :: account_mb_source

contains

!-------------------------------------------------------------------------------
!> Initialisations for the ice thickness computation.
!<------------------------------------------------------------------------------
subroutine calc_thk_init(st,grd,par)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_par_class), intent(inout) :: par

integer :: i, j


!-------- Computation/initialisation of the ice base topography
!                                       and its time derivative --------

if (par%margin==1 .or. par%margin==2) then ! /* only grounded ice */

  st%zb_neu   = st%zl_neu
  st%dzb_dtau = st%dzl_dtau

else if (par%margin==3) then !  /* grounded and floating ice */

  where (st%maske <= 1)   ! grounded ice or ice-free land
    st%zb_neu   = st%zl_neu
    st%dzb_dtau = st%dzl_dtau
  elsewhere   ! (maske >= 2; ocean or floating ice)
    st%zb_neu   = st%zb       ! initialisation,
    st%dzb_dtau = 0.0_wp   ! will be overwritten later
  end where

endif

!-------- Initialisation of the ice thickness
!                               and surface topography --------

st%H = st%H_c + st%H_t

st%zs_neu = st%zs   ! initialisation,
st%H_neu  = st%H    ! will be overwritten later

!-------- Solver type --------

if (par%calcthk==1 .or. par%calcthk==4) then ! /* explicit solver */

  flag_solver_explicit = .true.

else if (par%calcthk==2 .or. par%calcthk==3 .or. par%calcthk==5 .or. par%calcthk==6) then ! /* implicit solver */

  flag_solver_explicit = .false.

endif

if(par%i_advance==1) then
! rc no positive smb over ocean for ice advance scheme
!  where(st%maske==2) st%as_perp=min(st%as_perp,-0.001_wp/sec_year)

   do i=1, grd%IMAX-1
   do j=1, grd%JMAX-1

      if ( (st%maske(j,i)==2) &   ! ocean
           .and. &
             (     (st%maske(j,i+1)==2)   &   ! with
              .and.(st%maske(j,i-1)==2)   &   ! all ocean
              .and.(st%maske(j+1,i)==2)   &   ! neighbouring
              .and.(st%maske(j-1,i)==2) ) &   ! points
         ) &
         st%as_perp(j,i) = min(st%as_perp(j,i),-0.001_wp/sec_year)

       enddo
     enddo

end if

if(par%i_tibet==1) then
  where(67.0_wp*pi_180.le.grd%lambda.and.grd%lambda.le.106.0_wp*pi_180 &
   .and.23.0_wp*pi_180.le.grd%phi.and.grd%phi.le.47.0_wp*pi_180) &
                                        st%as_perp=min(st%as_perp,-0.1_wp/sec_year)
end if

!-------- Source term for the ice thickness equation --------

  st%mb_source = st%as_perp - st%Q_b_tot - st%calving

end subroutine calc_thk_init

!-------------------------------------------------------------------------------
!> Explicit solver for the diffusive SIA ice surface equation.
!<------------------------------------------------------------------------------
subroutine calc_thk_sia_expl(st,grd,tmr,par)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_timer_class), intent(in) :: tmr
type(sico_par_class), intent(in) :: par

integer                       :: i, j

!-------- Abbreviations --------

st%czs2 = 0.0_wp
st%czs3 = 0.0_wp

do i=0, grd%IMAX-1
do j=0, grd%JMAX
   st%czs2(j,i) = grd%azs2*0.5_wp*(st%h_diff(j,i)+st%h_diff(j,i+1)) &
               *grd%sq_g22_sgx(j,i)*grd%insq_g11_sgx(j,i)
end do
end do

do i=0, grd%IMAX
do j=0, grd%JMAX-1
   st%czs3(j,i) = grd%azs3*0.5_wp*(st%h_diff(j,i)+st%h_diff(j+1,i)) &
               *grd%sq_g11_sgy(j,i)*grd%insq_g22_sgy(j,i)
end do
end do

!-------- Solution of the explicit scheme --------

do i=0, grd%IMAX
  do j=0, grd%JMAX

    if ( grd%flag_inner_point(j,i) )  then

      st%zs_neu(j,i) = st%zs(j,i) &
        + tmr%dtime*(st%mb_source(j,i)+st%dzb_dtau(j,i)) &
        + ( st%czs2(j,i)  *(st%zs(j,i+1)-st%zs(j,i)  ) &
        -st%czs2(j,i-1)*(st%zs(j,i)  -st%zs(j,i-1)) &
        +st%czs3(j,i)  *(st%zs(j+1,i)-st%zs(j,i)  ) &
        -st%czs3(j-1,i)*(st%zs(j,i)  -st%zs(j-1,i)) ) &
        *grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i)

    else
      st%zs_neu(j,i) = st%zb_neu(j,i)   ! zero-thickness boundary condition
    end if

  end do
end do


!-------- Ice thickness --------

st%H_neu = st%zs_neu - st%zb_neu

!-------- Applying the source term --------

call apply_mb_source(st,grd,tmr,par)

!-------- Adjusting the ice thickness, if needed --------

call thk_adjust(st,tmr,par)

end subroutine calc_thk_sia_expl

!-------------------------------------------------------------------------------
!> Over-implicit solver for the diffusive SIA ice surface equation.
!<------------------------------------------------------------------------------
subroutine calc_thk_sia_impl_sor(st,grd,tmr,par)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_timer_class), intent(in) :: tmr
type(sico_par_class), intent(in) :: par

integer                       :: i, j
integer                       :: k, nnz

integer                            :: ierr
integer                            :: nc, nr
integer                            :: nmax 
integer                            :: n_sprs 
integer, allocatable, dimension(:) :: lgs_a_ptr, lgs_a_index
integer, allocatable, dimension(:) :: lgs_a_diag_index
integer, allocatable, dimension(:) :: lgs_a_index_trim
real(wp),     allocatable, dimension(:) :: lgs_a_value, lgs_b_value, lgs_x_value
real(wp),     allocatable, dimension(:) :: lgs_a_value_trim
real(wp)                                :: eps_sor

nmax = (grd%IMAX+1)*(grd%JMAX+1)
n_sprs = 10*(grd%IMAX+1)*(grd%JMAX+1)

!-------- Abbreviations --------

st%czs2 = 0.0_wp
st%czs3 = 0.0_wp

do i=0, grd%IMAX-1
do j=0, grd%JMAX
   st%czs2(j,i) = grd%azs2*0.5_wp*(st%h_diff(j,i)+st%h_diff(j,i+1)) &
               *grd%sq_g22_sgx(j,i)*grd%insq_g11_sgx(j,i)
end do
end do

do i=0, grd%IMAX
do j=0, grd%JMAX-1
   st%czs3(j,i) = grd%azs3*0.5_wp*(st%h_diff(j,i)+st%h_diff(j+1,i)) &
               *grd%sq_g11_sgy(j,i)*grd%insq_g22_sgy(j,i)
end do
end do

!-------- Assembly of the system of linear equations
!                     (matrix storage: compressed sparse row CSR) --------

allocate(lgs_a_value(n_sprs), lgs_a_index(n_sprs), lgs_a_ptr(nmax+1))
allocate(lgs_a_diag_index(nmax), lgs_b_value(nmax), lgs_x_value(nmax))

lgs_a_value = 0.0_wp
lgs_a_index = 0
lgs_a_ptr   = 0

lgs_b_value = 0.0_wp
lgs_x_value = 0.0_wp

lgs_a_ptr(1) = 1

k = 0

do nr=1, nmax   ! loop over rows

   i = grd%n2i(nr)
   j = grd%n2j(nr)

   if ( grd%flag_inner_point(j,i) ) then

      nc = grd%ij2n(j,i-1)   ! smallest nc (column counter), for zs(j,i-1)
      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k) = -st%czs2(j,i-1)*par%ovi_weight &
                        *grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i)
      lgs_a_index(k) = nc

      nc = grd%ij2n(j-1,i)   ! next nc (column counter), for zs(j-1,i)
      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k) = -st%czs3(j-1,i)*par%ovi_weight &
                        *grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i)
      lgs_a_index(k) = nc

      nc = grd%ij2n(j,i)     ! next nc (column counter), for zs(j,i)
      ! if (nc /= nr) &                     (diagonal element)
      !    stop ' >>> calc_thk_sia_impl: Check for diagonal element failed!'
      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k) = 1.0_wp &
                       + (st%czs2(j,i)+st%czs2(j,i-1)+st%czs3(j,i)+st%czs3(j-1,i)) &
                         *par%ovi_weight &
                         *grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i)
      lgs_a_diag_index(nr) = k
      lgs_a_index(k) = nc

      nc = grd%ij2n(j+1,i)   ! next nc (column counter), for zs(j+1,i)
      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k) = -st%czs3(j,i)*par%ovi_weight &
                        *grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i)
      lgs_a_index(k) = nc

      nc = grd%ij2n(j,i+1)   ! largest nc (column counter), for zs(j,i+1)
      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k) = -st%czs2(j,i)*par%ovi_weight &
                        *grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i)
      lgs_a_index(k) = nc

      lgs_b_value(nr) = st%zs(j,i) &
                          + tmr%dtime*(st%mb_source(j,i)+st%dzb_dtau(j,i)) &
                          + ( st%czs2(j,i)*(st%zs(j,i+1)-st%zs(j,i)) &
                             -st%czs2(j,i-1)*(st%zs(j,i)-st%zs(j,i-1)) &
                             +st%czs3(j,i)*(st%zs(j+1,i)-st%zs(j,i)) &
                             -st%czs3(j-1,i)*(st%zs(j,i)-st%zs(j-1,i)) ) &
                            *(1.0_wp-par%ovi_weight) &
                            *grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i)
                                                          ! right-hand side

   else   !  zero-thickness boundary condition

      k = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k)       = 1.0_wp   ! diagonal element only
      lgs_a_diag_index(nr) = k
      lgs_a_index(k)       = nr
      lgs_b_value(nr)      = st%zb_neu(j,i)

   end if

   lgs_x_value(nr) = st%zs(j,i)   ! old surface topography,
                               ! initial guess for solution vector

   lgs_a_ptr(nr+1) = k+1   ! row is completed, store index to next row

end do

nnz = k   ! number of non-zero elements of the matrix

!-------- Solution of the system of linear equations --------

!  ------ Solution with the built-in SOR solver

allocate(lgs_a_value_trim(nnz), lgs_a_index_trim(nnz))

do k=1, nnz   ! relocate matrix to trimmed arrays
   lgs_a_value_trim(k) = lgs_a_value(k)
   lgs_a_index_trim(k) = lgs_a_index(k)
end do

deallocate(lgs_a_value, lgs_a_index)

eps_sor = 1.0e-05_wp*st%mean_accum*tmr%dtime   ! convergence parameter

call sor_sprs(lgs_a_value_trim, &
              lgs_a_index_trim, lgs_a_diag_index, lgs_a_ptr, &
              lgs_b_value, &
              nnz, nmax, par%omega_sor, eps_sor, lgs_x_value, ierr)

do nr=1, nmax
   i = grd%n2i(nr)
   j = grd%n2j(nr)
   st%zs_neu(j,i) = lgs_x_value(nr)
end do

deallocate(lgs_a_value_trim, lgs_a_index_trim, lgs_a_ptr)
deallocate(lgs_a_diag_index, lgs_b_value, lgs_x_value)

!-------- Ice thickness --------

st%H_neu = st%zs_neu - st%zb_neu

!-------- Applying the source term --------

call apply_mb_source(st,grd,tmr,par)

!-------- Adjusting the ice thickness, if needed --------

call thk_adjust(st,tmr,par)

end subroutine calc_thk_sia_impl_sor

!-------------------------------------------------------------------------------
!> Over-implicit solver for the diffusive SIA ice surface equation.
!<------------------------------------------------------------------------------
subroutine calc_thk_sia_impl_lis(st,grd,tmr,par)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_timer_class), intent(in) :: tmr
type(sico_par_class), intent(in) :: par

integer                       :: i, j
integer                       :: k, nnz

LIS_INTEGER                            :: ierr
LIS_INTEGER                            :: iter
LIS_INTEGER                            :: nc, nr
LIS_INTEGER                            :: nmax   
LIS_INTEGER                            :: n_sprs 
LIS_INTEGER, allocatable, dimension(:) :: lgs_a_ptr, lgs_a_index
LIS_INTEGER, allocatable, dimension(:) :: lgs_a_diag_index
LIS_MATRIX                             :: lgs_a
LIS_VECTOR                             :: lgs_b, lgs_x
LIS_SCALAR,  allocatable, dimension(:) :: lgs_a_value, lgs_b_value, lgs_x_value
LIS_SOLVER                             :: solver
character(len=256)                     :: ch_solver_set_option

nmax = (grd%IMAX+1)*(grd%JMAX+1)
n_sprs = 10*(grd%IMAX+1)*(grd%JMAX+1)

!-------- Abbreviations --------

st%czs2 = 0.0_wp
st%czs3 = 0.0_wp

do i=0, grd%IMAX-1
do j=0, grd%JMAX
   st%czs2(j,i) = grd%azs2*0.5_wp*(st%h_diff(j,i)+st%h_diff(j,i+1)) &
               *grd%sq_g22_sgx(j,i)*grd%insq_g11_sgx(j,i)
end do
end do

do i=0, grd%IMAX
do j=0, grd%JMAX-1
   st%czs3(j,i) = grd%azs3*0.5_wp*(st%h_diff(j,i)+st%h_diff(j+1,i)) &
               *grd%sq_g11_sgy(j,i)*grd%insq_g22_sgy(j,i)
end do
end do

!-------- Assembly of the system of linear equations
!                     (matrix storage: compressed sparse row CSR) --------

allocate(lgs_a_value(n_sprs), lgs_a_index(n_sprs), lgs_a_ptr(nmax+1))
allocate(lgs_a_diag_index(nmax), lgs_b_value(nmax), lgs_x_value(nmax))

lgs_a_value = 0.0_wp
lgs_a_index = 0
lgs_a_ptr   = 0

lgs_b_value = 0.0_wp
lgs_x_value = 0.0_wp

lgs_a_ptr(1) = 1

k = 0

do nr=1, nmax   ! loop over rows

   i = grd%n2i(nr)
   j = grd%n2j(nr)

   if ( grd%flag_inner_point(j,i) ) then

      nc = grd%ij2n(j,i-1)   ! smallest nc (column counter), for zs(j,i-1)
      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k) = -st%czs2(j,i-1)*par%ovi_weight &
                        *grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i)
      lgs_a_index(k) = nc

      nc = grd%ij2n(j-1,i)   ! next nc (column counter), for zs(j-1,i)
      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k) = -st%czs3(j-1,i)*par%ovi_weight &
                        *grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i)
      lgs_a_index(k) = nc

      nc = grd%ij2n(j,i)     ! next nc (column counter), for zs(j,i)
      ! if (nc /= nr) &                     (diagonal element)
      !    stop ' >>> calc_thk_sia_impl: Check for diagonal element failed!'
      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k) = 1.0_wp &
                       + (st%czs2(j,i)+st%czs2(j,i-1)+st%czs3(j,i)+st%czs3(j-1,i)) &
                         *par%ovi_weight &
                         *grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i)
      lgs_a_diag_index(nr) = k
      lgs_a_index(k) = nc

      nc = grd%ij2n(j+1,i)   ! next nc (column counter), for zs(j+1,i)
      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k) = -st%czs3(j,i)*par%ovi_weight &
                        *grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i)
      lgs_a_index(k) = nc

      nc = grd%ij2n(j,i+1)   ! largest nc (column counter), for zs(j,i+1)
      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k) = -st%czs2(j,i)*par%ovi_weight &
                        *grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i)
      lgs_a_index(k) = nc

      lgs_b_value(nr) = st%zs(j,i) &
                          + tmr%dtime*(st%mb_source(j,i)+st%dzb_dtau(j,i)) &
                          + ( st%czs2(j,i)*(st%zs(j,i+1)-st%zs(j,i)) &
                             -st%czs2(j,i-1)*(st%zs(j,i)-st%zs(j,i-1)) &
                             +st%czs3(j,i)*(st%zs(j+1,i)-st%zs(j,i)) &
                             -st%czs3(j-1,i)*(st%zs(j,i)-st%zs(j-1,i)) ) &
                            *(1.0_wp-par%ovi_weight) &
                            *grd%insq_g11_g(j,i)*grd%insq_g22_g(j,i)
                                                          ! right-hand side

   else   !  zero-thickness boundary condition

      k = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k)       = 1.0_wp   ! diagonal element only
      lgs_a_diag_index(nr) = k
      lgs_a_index(k)       = nr
      lgs_b_value(nr)      = st%zb_neu(j,i)

   end if

   lgs_x_value(nr) = st%zs(j,i)   ! old surface topography,
                               ! initial guess for solution vector

   lgs_a_ptr(nr+1) = k+1   ! row is completed, store index to next row

end do

nnz = k   ! number of non-zero elements of the matrix

!-------- Solution of the system of linear equations --------

!  ------ Settings for Lis

call lis_initialize(ierr)
call CHKERR(ierr)

call lis_matrix_create(LIS_COMM_WORLD, lgs_a, ierr)
call lis_vector_create(LIS_COMM_WORLD, lgs_b, ierr)
call lis_vector_create(LIS_COMM_WORLD, lgs_x, ierr)

call lis_matrix_set_size(lgs_a, 0, nmax, ierr)
call lis_vector_set_size(lgs_b, 0, nmax, ierr)
call lis_vector_set_size(lgs_x, 0, nmax, ierr)

do nr=1, nmax

   do nc=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
      call lis_matrix_set_value(LIS_INS_VALUE, nr, lgs_a_index(nc), &
                                               lgs_a_value(nc), lgs_a, ierr)
   end do

   call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_b_value(nr), lgs_b, ierr)
   call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_x_value(nr), lgs_x, ierr)

end do

call lis_matrix_set_type(lgs_a, LIS_MATRIX_CSR, ierr)
call lis_matrix_assemble(lgs_a, ierr)

!  ------ Solution with Lis

call lis_solver_create(solver, ierr)

ch_solver_set_option = '-i bicg -p ilu '// &
                       '-maxiter 1000 -tol 1.0e-12 -initx_zeros false'

call lis_solver_set_option(trim(ch_solver_set_option), solver, ierr)
call CHKERR(ierr)

call lis_solve(lgs_a, lgs_b, lgs_x, solver, ierr)
call CHKERR(ierr)

call lis_solver_get_iter(solver, iter, ierr)
write(6,'(10x,a,i0)') 'calc_thk_sia_impl: iter = ', iter

lgs_x_value = 0.0_wp
call lis_vector_gather(lgs_x, lgs_x_value, ierr)
call lis_matrix_destroy(lgs_a, ierr)
call lis_vector_destroy(lgs_b, ierr)
call lis_vector_destroy(lgs_x, ierr)
call lis_solver_destroy(solver, ierr)

call lis_finalize(ierr)      
call CHKERR(ierr)

do nr=1, nmax
   i = grd%n2i(nr)
   j = grd%n2j(nr)
   st%zs_neu(j,i) = lgs_x_value(nr)
end do

deallocate(lgs_a_value, lgs_a_index, lgs_a_ptr)
deallocate(lgs_a_diag_index, lgs_b_value, lgs_x_value)

!-------- Ice thickness --------

st%H_neu = st%zs_neu - st%zb_neu

!-------- Applying the source term --------

call apply_mb_source(st,grd,tmr,par)

!-------- Adjusting the ice thickness, if needed --------

call thk_adjust(st,tmr,par)

end subroutine calc_thk_sia_impl_lis

!-------------------------------------------------------------------------------
!> Explicit solver for the general ice thickness equation.
!<------------------------------------------------------------------------------
subroutine calc_thk_expl(st,grd,tmr,par)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_timer_class), intent(in) :: tmr
type(sico_par_class), intent(in) :: par

integer                       :: i, j
real(wp)                      :: H_eff_tmp


! test: ensure zero velocities in ocean
!--do i=0, grd%IMAX
!--do j=0, grd%JMAX
!--   if (grd%flag_inner_point(j,i)) then
!--     if(st%maske(j,i).eq.2.and.st%maske(j,i+1).eq.2) st%vx_m(j,i) = 0.0_wp
!--     if(st%maske(j,i).eq.2.and.st%maske(j+1,i).eq.2) st%vy_m(j,i) = 0.0_wp
!--   end if
!--end do
!--end do

!-------- Special control of advance of shelf ice front  --------

if(par%i_advance==1) then

! is this really necessary?
! st%flag_calving_front_1(j,i) upto date?
   st%flag_calving_front_1 = .false.
   do i=1, grd%IMAX-1
   do j=1, grd%JMAX-1

      if ( (st%maske(j,i)==3) &   ! floating ice
           .and. &
             (    (st%maske(j,i+1)==2)   &   ! with
              .or.(st%maske(j,i-1)==2)   &   ! one
              .or.(st%maske(j+1,i)==2)   &   ! neighbouring
              .or.(st%maske(j-1,i)==2) ) &   ! sea point
         ) &
         st%flag_calving_front_1(j,i) = .true.   ! detection
                                                 ! of the calving front
   end do
   end do


   st%H_eff=0
   do i=0, grd%IMAX
   do j=0, grd%JMAX
      if ( st%flag_calving_front_1(j,i) .and. grd%flag_inner_point(j,i) ) then
         H_eff_tmp = 1.e30_wp
        ! eventually exclude maske=0 query
        if( (st%maske(j,i-1).eq.0 .or. st%maske(j,i-1).eq.3).and. &
              .not.st%flag_calving_front_1(j,i-1) .and. st%H(j,i-1)<H_eff_tmp) &
            H_eff_tmp=st%H(j,i-1)
         if( (st%maske(j,i+1).eq.0 .or. st%maske(j,i+1).eq.3).and. &
              .not.st%flag_calving_front_1(j,i+1) .and. st%H(j,i+1)<H_eff_tmp) &
            H_eff_tmp=st%H(j,i+1)
         if( (st%maske(j-1,i).eq.0 .or. st%maske(j-1,i).eq.3).and. &
              .not.st%flag_calving_front_1(j-1,i) .and. st%H(j-1,i)<H_eff_tmp) &
            H_eff_tmp=st%H(j-1,i)
         if( (st%maske(j+1,i).eq.0 .or. st%maske(j+1,i).eq.3).and. &
              .not.st%flag_calving_front_1(j+1,i) .and. st%H(j+1,i)<H_eff_tmp) &
            H_eff_tmp=st%H(j+1,i)
 
         st%H_eff(j,i) = H_eff_tmp

! test: no upstream flow into shelf ice ocean front 
         if (st%vx_m(j,i-1) >= 0.0_wp .and. st%maske(j,i-1).eq.2) st%vx_m(j,i-1) = 0.0_wp
         if (st%vx_m(j,i)   <  0.0_wp .and. st%maske(j,i+1).eq.2) st%vx_m(j,i)   = 0.0_wp 
         if (st%vy_m(j-1,i) >= 0.0_wp .and. st%maske(j-1,i).eq.2) st%vy_m(j-1,i) = 0.0_wp
         if (st%vy_m(j,i)   <  0.0_wp .and. st%maske(j+1,i).eq.2) st%vy_m(j,i)   = 0.0_wp

         if(st%H(j,i)<H_eff_tmp) then
!           determine downstream velocity and set zero (st%maske(j,i).eq.3)
            if (st%maske(j,i-1).eq.2) st%vx_m(j,i-1) = 0.0_wp
            if (st%maske(j,i+1).eq.2) st%vx_m(j,i)   = 0.0_wp 
            if (st%maske(j-1,i).eq.2) st%vy_m(j-1,i) = 0.0_wp
            if (st%maske(j+1,i).eq.2) st%vy_m(j,i)   = 0.0_wp
         end if
      end if
   end do
   end do

end if

!-------- Abbreviations for upstream scheme --------

do i=0, grd%IMAX
do j=0, grd%JMAX

   if (grd%flag_inner_point(j,i)) then

      st%vx_m_1(j,i) = st%vx_m(j,i-1)
      st%vx_m_2(j,i) = st%vx_m(j,i)
      st%vy_m_1(j,i) = st%vy_m(j-1,i)
      st%vy_m_2(j,i) = st%vy_m(j,i)

      if (st%vx_m_1(j,i) >= 0.0_wp) then    ! up
         st%upH_x_1(j,i) = st%H(j,i-1)
      else                                  ! down
         st%upH_x_1(j,i) = st%H(j,i)
      end if

      if (st%vx_m_2(j,i) >= 0.0_wp) then    ! down
         st%upH_x_2(j,i) = st%H(j,i)
      else                                  ! up
         st%upH_x_2(j,i) = st%H(j,i+1)
      end if

      if (st%vy_m_1(j,i) >= 0.0_wp) then     ! up
         st%upH_y_1(j,i) = st%H(j-1,i)
      else                                   ! down
         st%upH_y_1(j,i) = st%H(j,i)
      end if

      if (st%vy_m_2(j,i) >= 0.0_wp) then     ! down
         st%upH_y_2(j,i) = st%H(j,i)
      else                                   ! up
         st%upH_y_2(j,i) = st%H(j+1,i)
      end if

   else   ! .not.(flag_inner_point(j,i))

      st%vx_m_1(j,i) = 0.0_wp
      st%vx_m_2(j,i) = 0.0_wp
      st%vy_m_1(j,i) = 0.0_wp
      st%vy_m_2(j,i) = 0.0_wp

      st%upH_x_1(j,i) = 0.0_wp
      st%upH_x_2(j,i) = 0.0_wp
      st%upH_y_1(j,i) = 0.0_wp
      st%upH_y_2(j,i) = 0.0_wp

   end if

end do
end do

!-------- Solution of the explicit scheme --------

  where ( grd%flag_inner_point )   ! inner point

    st%H_neu = st%H + tmr%dtime*st%mb_source &
      - grd%dt_darea &
      * (  ( st%vx_m_2*st%upH_x_2*grd%sq_g22_x_2*grd%deta   &
      -st%vx_m_1*st%upH_x_1*grd%sq_g22_x_1*grd%deta ) &
      + ( st%vy_m_2*st%upH_y_2*grd%sq_g11_y_2*grd%dxi    &
      -st%vy_m_1*st%upH_y_1*grd%sq_g11_y_1*grd%dxi  ) )

  elsewhere
    st%H_neu = 0.0_wp   ! zero-thickness boundary condition
  end where


!-------- Applying the source term --------

call apply_mb_source(st,grd,tmr,par)

!-------- Adjusting the ice thickness, if needed --------

call thk_adjust(st,tmr,par)

end subroutine calc_thk_expl

!-------------------------------------------------------------------------------
!> Over-implicit solver for the general ice thickness equation.
!<------------------------------------------------------------------------------
subroutine calc_thk_impl_sor(st,grd,tmr,par)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_timer_class), intent(in) :: tmr
type(sico_par_class), intent(in) :: par

integer                       :: i, j
integer                       :: k, nnz

integer                            :: ierr
integer                            :: nc, nr
integer                            :: nmax   
integer                            :: n_sprs 
integer, allocatable, dimension(:) :: lgs_a_ptr, lgs_a_index
integer, allocatable, dimension(:) :: lgs_a_diag_index
integer, allocatable, dimension(:) :: lgs_a_index_trim
real(wp),     allocatable, dimension(:) :: lgs_a_value, lgs_b_value, lgs_x_value
real(wp),     allocatable, dimension(:) :: lgs_a_value_trim
real(wp)                                :: eps_sor

nmax = (grd%IMAX+1)*(grd%JMAX+1)
n_sprs = 10*(grd%IMAX+1)*(grd%JMAX+1)

!-------- Abbreviations --------

do i=0, grd%IMAX
do j=0, grd%JMAX

   if (grd%flag_inner_point(j,i)) then

      st%vx_m_1(j,i) = st%vx_m(j,i-1)
      st%vx_m_2(j,i) = st%vx_m(j,i)
      st%vy_m_1(j,i) = st%vy_m(j-1,i)
      st%vy_m_2(j,i) = st%vy_m(j,i)

      if (st%vx_m_1(j,i) >= 0.0_wp) then
         st%upH_x_1(j,i) = st%H(j,i-1)
      else
         st%upH_x_1(j,i) = st%H(j,i)
      end if

      if (st%vx_m_2(j,i) >= 0.0_wp) then
         st%upH_x_2(j,i) = st%H(j,i)
      else
         st%upH_x_2(j,i) = st%H(j,i+1)
      end if

      if (st%vy_m_1(j,i) >= 0.0_wp) then
         st%upH_y_1(j,i) = st%H(j-1,i)
      else
         st%upH_y_1(j,i) = st%H(j,i)
      end if

      if (st%vy_m_2(j,i) >= 0.0_wp) then
         st%upH_y_2(j,i) = st%H(j,i)
      else
         st%upH_y_2(j,i) = st%H(j+1,i)
      end if

   else   ! .not.(flag_inner_point(j,i))

      st%vx_m_1(j,i) = 0.0_wp
      st%vx_m_2(j,i) = 0.0_wp
      st%vy_m_1(j,i) = 0.0_wp
      st%vy_m_2(j,i) = 0.0_wp

      st%upH_x_1(j,i) = 0.0_wp
      st%upH_x_2(j,i) = 0.0_wp
      st%upH_y_1(j,i) = 0.0_wp
      st%upH_y_2(j,i) = 0.0_wp

   end if

end do
end do

!-------- Assembly of the system of linear equations
!                     (matrix storage: compressed sparse row CSR) --------

allocate(lgs_a_value(n_sprs), lgs_a_index(n_sprs), lgs_a_ptr(nmax+1))
allocate(lgs_a_diag_index(nmax), lgs_b_value(nmax), lgs_x_value(nmax))

lgs_a_value = 0.0_wp
lgs_a_index = 0
lgs_a_ptr   = 0

lgs_b_value = 0.0_wp
lgs_x_value = 0.0_wp

lgs_a_ptr(1) = 1

k = 0

do nr=1, nmax   ! loop over rows

   i = grd%n2i(nr)
   j = grd%n2j(nr)

   if ( grd%flag_inner_point(j,i) )  then

      k=k+1 ; nc=grd%ij2n(j,i-1) ; lgs_a_index(k)=nc   ! for H(j,i-1)
      if (st%vx_m_1(j,i) > 0.0_wp) &
         lgs_a_value(k) = -grd%dt_darea(j,i)*st%vx_m_1(j,i) &
                                        *grd%sq_g22_x_1(j,i)*grd%deta*par%ovi_weight

      k=k+1 ; nc=grd%ij2n(j-1,i) ; lgs_a_index(k)=nc   ! for H(j-1,i)
      if (st%vy_m_1(j,i) > 0.0_wp) &
         lgs_a_value(k) = -grd%dt_darea(j,i)*st%vy_m_1(j,i) &
                                        *grd%sq_g11_y_1(j,i)*grd%dxi*par%ovi_weight

      k=k+1 ; lgs_a_index(k)=nr ; lgs_a_diag_index(nr)=k  ! for H(j,i)
      lgs_a_value(k) = 1.0_wp                             ! (diagonal element)
      if (st%vy_m_1(j,i) < 0.0_wp) &
         lgs_a_value(k) = lgs_a_value(k) &
                          - grd%dt_darea(j,i)*st%vy_m_1(j,i) &
                                         *grd%sq_g11_y_1(j,i)*grd%dxi*par%ovi_weight
      if (st%vx_m_1(j,i) < 0.0_wp) &
         lgs_a_value(k) = lgs_a_value(k) &
                          - grd%dt_darea(j,i)*st%vx_m_1(j,i) &
                                         *grd%sq_g22_x_1(j,i)*grd%deta*par%ovi_weight
      if (st%vx_m_2(j,i) > 0.0_wp) &
         lgs_a_value(k) = lgs_a_value(k) &
                          + grd%dt_darea(j,i)*st%vx_m_2(j,i) &
                                         *grd%sq_g22_x_2(j,i)*grd%deta*par%ovi_weight
      if (st%vy_m_2(j,i) > 0.0_wp) &
         lgs_a_value(k) = lgs_a_value(k) &
                          + grd%dt_darea(j,i)*st%vy_m_2(j,i) &
                                         *grd%sq_g11_y_2(j,i)*grd%dxi*par%ovi_weight

      k=k+1 ; nc=grd%ij2n(j+1,i) ; lgs_a_index(k)=nc   ! for H(j+1,i)
      if (st%vy_m_2(j,i) < 0.0_wp) &
         lgs_a_value(k) = grd%dt_darea(j,i)*st%vy_m_2(j,i) &
                                       *grd%sq_g11_y_2(j,i)*grd%dxi*par%ovi_weight

      k=k+1 ; nc=grd%ij2n(j,i+1) ; lgs_a_index(k)=nc   ! for H(j,i+1)
      if (st%vx_m_2(j,i) < 0.0_wp) &
         lgs_a_value(k) = grd%dt_darea(j,i)*st%vx_m_2(j,i) &
                                       *grd%sq_g22_x_2(j,i)*grd%deta*par%ovi_weight

      lgs_b_value(nr) = st%H(j,i) &
                        +tmr%dtime*st%mb_source(j,i) &
                        -(1.0_wp-par%ovi_weight) &
                           * grd%dt_darea(j,i) &
                             * (  ( st%vx_m_2(j,i)*st%upH_x_2(j,i) &
                                               *grd%sq_g22_x_2(j,i)*grd%deta   &
                                   -st%vx_m_1(j,i)*st%upH_x_1(j,i) &
                                               *grd%sq_g22_x_1(j,i)*grd%deta ) &
                                + ( st%vy_m_2(j,i)*st%upH_y_2(j,i) &
                                               *grd%sq_g11_y_2(j,i)*grd%dxi    &
                                   -st%vy_m_1(j,i)*st%upH_y_1(j,i) &
                                               *grd%sq_g11_y_1(j,i)*grd%dxi  ) )
                                                          ! right-hand side

   else   ! zero-thickness boundary condition

      k = k+1
      lgs_a_value(k)       = 1.0_wp   ! diagonal element only
      lgs_a_diag_index(nr) = k
      lgs_a_index(k)       = nr
      lgs_b_value(nr)      = 0.0_wp

   end if

   lgs_x_value(nr) = st%H(j,i)   ! old ice thickness,
                              ! initial guess for solution vector

   lgs_a_ptr(nr+1) = k+1   ! row is completed, store index to next row

end do

nnz = k   ! number of non-zero elements of the matrix

!-------- Solution of the system of linear equations --------

!  ------ Solution with the built-in SOR solver

allocate(lgs_a_value_trim(nnz), lgs_a_index_trim(nnz))

do k=1, nnz   ! relocate matrix to trimmed arrays
   lgs_a_value_trim(k) = lgs_a_value(k)
   lgs_a_index_trim(k) = lgs_a_index(k)
end do

deallocate(lgs_a_value, lgs_a_index)

eps_sor = 1.0e-05_wp*st%mean_accum*tmr%dtime   ! convergence parameter

call sor_sprs(lgs_a_value_trim, &
              lgs_a_index_trim, lgs_a_diag_index, lgs_a_ptr, &
              lgs_b_value, &
              nnz, nmax, par%omega_sor, eps_sor, lgs_x_value, ierr)

do nr=1, nmax
   i = grd%n2i(nr)
   j = grd%n2j(nr)
   st%H_neu(j,i) = lgs_x_value(nr)
end do

deallocate(lgs_a_value_trim, lgs_a_index_trim, lgs_a_ptr)
deallocate(lgs_a_diag_index, lgs_b_value, lgs_x_value)

!-------- Applying the source term --------

call apply_mb_source(st,grd,tmr,par)

!-------- Adjusting the ice thickness, if needed --------

call thk_adjust(st,tmr,par)

end subroutine calc_thk_impl_sor

!-------------------------------------------------------------------------------
!> Over-implicit solver for the general ice thickness equation.
!<------------------------------------------------------------------------------
subroutine calc_thk_impl_lis(st,grd,tmr,par)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_timer_class), intent(in) :: tmr
type(sico_par_class), intent(in) :: par

integer                       :: i, j
integer                       :: k, nnz

LIS_INTEGER                            :: ierr
LIS_INTEGER                            :: iter
LIS_INTEGER                            :: nc, nr
LIS_INTEGER                            :: nmax  
LIS_INTEGER                            :: n_sprs
LIS_INTEGER, allocatable, dimension(:) :: lgs_a_ptr, lgs_a_index
LIS_INTEGER, allocatable, dimension(:) :: lgs_a_diag_index
LIS_MATRIX                             :: lgs_a
LIS_VECTOR                             :: lgs_b, lgs_x
LIS_SCALAR,  allocatable, dimension(:) :: lgs_a_value, lgs_b_value, lgs_x_value
LIS_SOLVER                             :: solver
character(len=256)                     :: ch_solver_set_option

nmax = (grd%IMAX+1)*(grd%JMAX+1)
n_sprs = 10*(grd%IMAX+1)*(grd%JMAX+1)

!-------- Abbreviations --------

do i=0, grd%IMAX
do j=0, grd%JMAX

   if (grd%flag_inner_point(j,i)) then

      st%vx_m_1(j,i) = st%vx_m(j,i-1)
      st%vx_m_2(j,i) = st%vx_m(j,i)
      st%vy_m_1(j,i) = st%vy_m(j-1,i)
      st%vy_m_2(j,i) = st%vy_m(j,i)

      if (st%vx_m_1(j,i) >= 0.0_wp) then
         st%upH_x_1(j,i) = st%H(j,i-1)
      else
         st%upH_x_1(j,i) = st%H(j,i)
      end if

      if (st%vx_m_2(j,i) >= 0.0_wp) then
         st%upH_x_2(j,i) = st%H(j,i)
      else
         st%upH_x_2(j,i) = st%H(j,i+1)
      end if

      if (st%vy_m_1(j,i) >= 0.0_wp) then
         st%upH_y_1(j,i) = st%H(j-1,i)
      else
         st%upH_y_1(j,i) = st%H(j,i)
      end if

      if (st%vy_m_2(j,i) >= 0.0_wp) then
         st%upH_y_2(j,i) = st%H(j,i)
      else
         st%upH_y_2(j,i) = st%H(j+1,i)
      end if

   else   ! .not.(flag_inner_point(j,i))

      st%vx_m_1(j,i) = 0.0_wp
      st%vx_m_2(j,i) = 0.0_wp
      st%vy_m_1(j,i) = 0.0_wp
      st%vy_m_2(j,i) = 0.0_wp

      st%upH_x_1(j,i) = 0.0_wp
      st%upH_x_2(j,i) = 0.0_wp
      st%upH_y_1(j,i) = 0.0_wp
      st%upH_y_2(j,i) = 0.0_wp

   end if

end do
end do

!-------- Assembly of the system of linear equations
!                     (matrix storage: compressed sparse row CSR) --------

allocate(lgs_a_value(n_sprs), lgs_a_index(n_sprs), lgs_a_ptr(nmax+1))
allocate(lgs_a_diag_index(nmax), lgs_b_value(nmax), lgs_x_value(nmax))

lgs_a_value = 0.0_wp
lgs_a_index = 0
lgs_a_ptr   = 0

lgs_b_value = 0.0_wp
lgs_x_value = 0.0_wp

lgs_a_ptr(1) = 1

k = 0

do nr=1, nmax   ! loop over rows

   i = grd%n2i(nr)
   j = grd%n2j(nr)

   if ( grd%flag_inner_point(j,i) )  then

      k=k+1 ; nc=grd%ij2n(j,i-1) ; lgs_a_index(k)=nc   ! for H(j,i-1)
      if (st%vx_m_1(j,i) > 0.0_wp) &
         lgs_a_value(k) = -grd%dt_darea(j,i)*st%vx_m_1(j,i) &
                                        *grd%sq_g22_x_1(j,i)*grd%deta*par%ovi_weight

      k=k+1 ; nc=grd%ij2n(j-1,i) ; lgs_a_index(k)=nc   ! for H(j-2,i)
      if (st%vy_m_1(j,i) > 0.0_wp) &
         lgs_a_value(k) = -grd%dt_darea(j,i)*st%vy_m_1(j,i) &
                                        *grd%sq_g11_y_1(j,i)*grd%dxi*par%ovi_weight

      k=k+1 ; lgs_a_index(k)=nr ; lgs_a_diag_index(nr)=k  ! for H(j,i)
      lgs_a_value(k) = 1.0_wp                             ! (diagonal element)
      if (st%vy_m_1(j,i) < 0.0_wp) &
         lgs_a_value(k) = lgs_a_value(k) &
                          - grd%dt_darea(j,i)*st%vy_m_1(j,i) &
                                         *grd%sq_g11_y_1(j,i)*grd%dxi*par%ovi_weight
      if (st%vx_m_1(j,i) < 0.0_wp) &
         lgs_a_value(k) = lgs_a_value(k) &
                          - grd%dt_darea(j,i)*st%vx_m_1(j,i) &
                                         *grd%sq_g22_x_1(j,i)*grd%deta*par%ovi_weight
      if (st%vx_m_2(j,i) > 0.0_wp) &
         lgs_a_value(k) = lgs_a_value(k) &
                          + grd%dt_darea(j,i)*st%vx_m_2(j,i) &
                                         *grd%sq_g22_x_2(j,i)*grd%deta*par%ovi_weight
      if (st%vy_m_2(j,i) > 0.0_wp) &
         lgs_a_value(k) = lgs_a_value(k) &
                          + grd%dt_darea(j,i)*st%vy_m_2(j,i) &
                                         *grd%sq_g11_y_2(j,i)*grd%dxi*par%ovi_weight

      k=k+1 ; nc=grd%ij2n(j+1,i) ; lgs_a_index(k)=nc   ! for H(j+1,i)
      if (st%vy_m_2(j,i) < 0.0_wp) &
         lgs_a_value(k) = grd%dt_darea(j,i)*st%vy_m_2(j,i) &
                                       *grd%sq_g11_y_2(j,i)*grd%dxi*par%ovi_weight

      k=k+1 ; nc=grd%ij2n(j,i+1) ; lgs_a_index(k)=nc   ! for H(j,i+1)
      if (st%vx_m_2(j,i) < 0.0_wp) &
         lgs_a_value(k) = grd%dt_darea(j,i)*st%vx_m_2(j,i) &
                                       *grd%sq_g22_x_2(j,i)*grd%deta*par%ovi_weight

      lgs_b_value(nr) = st%H(j,i) &
                        +tmr%dtime*st%mb_source(j,i) &
                        -(1.0_wp-par%ovi_weight) &
                           * grd%dt_darea(j,i) &
                             * (  ( st%vx_m_2(j,i)*st%upH_x_2(j,i) &
                                               *grd%sq_g22_x_2(j,i)*grd%deta   &
                                   -st%vx_m_1(j,i)*st%upH_x_1(j,i) &
                                               *grd%sq_g22_x_1(j,i)*grd%deta ) &
                                + ( st%vy_m_2(j,i)*st%upH_y_2(j,i) &
                                               *grd%sq_g11_y_2(j,i)*grd%dxi    &
                                   -st%vy_m_1(j,i)*st%upH_y_1(j,i) &
                                               *grd%sq_g11_y_1(j,i)*grd%dxi  ) )
                                                          ! right-hand side

   else   ! zero-thickness boundary condition

      k = k+1
      lgs_a_value(k)       = 1.0_wp   ! diagonal element only
      lgs_a_diag_index(nr) = k
      lgs_a_index(k)       = nr
      lgs_b_value(nr)      = 0.0_wp

   end if

   lgs_x_value(nr) = st%H(j,i)   ! old ice thickness,
                              ! initial guess for solution vector

   lgs_a_ptr(nr+1) = k+1   ! row is completed, store index to next row

end do

nnz = k   ! number of non-zero elements of the matrix

!-------- Solution of the system of linear equations --------

!  ------ Settings for Lis

call lis_initialize(ierr)
call CHKERR(ierr)

call lis_matrix_create(LIS_COMM_WORLD, lgs_a, ierr)
call lis_vector_create(LIS_COMM_WORLD, lgs_b, ierr)
call lis_vector_create(LIS_COMM_WORLD, lgs_x, ierr)

call lis_matrix_set_size(lgs_a, 0, nmax, ierr)
call lis_vector_set_size(lgs_b, 0, nmax, ierr)
call lis_vector_set_size(lgs_x, 0, nmax, ierr)

do nr=1, nmax

   do nc=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
      call lis_matrix_set_value(LIS_INS_VALUE, nr, lgs_a_index(nc), &
                                               lgs_a_value(nc), lgs_a, ierr)
   end do

   call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_b_value(nr), lgs_b, ierr)
   call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_x_value(nr), lgs_x, ierr)

end do

call lis_matrix_set_type(lgs_a, LIS_MATRIX_CSR, ierr)
call lis_matrix_assemble(lgs_a, ierr)

!  ------ Solution with Lis

call lis_solver_create(solver, ierr)

ch_solver_set_option = '-i bicg -p ilu '// &
                       '-maxiter 1000 -tol 1.0e-12 -initx_zeros false'

call lis_solver_set_option(trim(ch_solver_set_option), solver, ierr)
call CHKERR(ierr)

call lis_solve(lgs_a, lgs_b, lgs_x, solver, ierr)
call CHKERR(ierr)

call lis_solver_get_iter(solver, iter, ierr)
write(6,'(10x,a,i0)') 'calc_thk_impl: iter = ', iter

lgs_x_value = 0.0_wp
call lis_vector_gather(lgs_x, lgs_x_value, ierr)
call lis_matrix_destroy(lgs_a, ierr)
call lis_vector_destroy(lgs_b, ierr)
call lis_vector_destroy(lgs_x, ierr)
call lis_solver_destroy(solver, ierr)

call lis_finalize(ierr)      
call CHKERR(ierr)

do nr=1, nmax
   i = grd%n2i(nr)
   j = grd%n2j(nr)
   st%H_neu(j,i) = lgs_x_value(nr)
end do

deallocate(lgs_a_value, lgs_a_index, lgs_a_ptr)
deallocate(lgs_a_diag_index, lgs_b_value, lgs_x_value)

!-------- Applying the source term --------

call apply_mb_source(st,grd,tmr,par)

!-------- Adjusting the ice thickness, if needed --------

call thk_adjust(st,tmr,par)

end subroutine calc_thk_impl_lis

!-------------------------------------------------------------------------------
!> Ice thickness evolution due to the source term (surface mass balance,
!! basal mass balance and calving).
!<------------------------------------------------------------------------------
subroutine apply_mb_source(st,grd,tmr,par)

  implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_timer_class), intent(in) :: tmr
type(sico_par_class), intent(in) :: par

integer :: i, j

!-------- Compute new ice thickness H_neu_flow due to glacial flow only
!         (no source term considered) --------

  st%H_neu_flow = 0.0_wp

  do i=0, grd%IMAX
  do j=0, grd%JMAX
     if (grd%flag_inner_point(j,i)) then
        st%H_neu_flow(j,i) = st%H_neu(j,i) - tmr%dtime*st%mb_source(j,i)
                          ! new ice thickness due to glacial flow only
     end if
  end do
  end do

!-------- Apply source term --------

  st%H_neu = 0.0_wp

  if (par%mb_account==0) then

    do i=0, grd%IMAX
      do j=0, grd%JMAX
        if (grd%flag_inner_point(j,i)) then
          st%H_neu(j,i) = max(st%H_neu_flow(j,i) + tmr%dtime*st%mb_source(j,i), 0.0_wp)
        end if
      end do
    end do

  else if (par%mb_account==1) then

    do i=2, grd%IMAX-2   ! outermost two grid points excepted,
      do j=2, grd%JMAX-2   ! required for accurate mass balance accounting
        st%H_neu(j,i) = max(st%H_neu_flow(j,i) + tmr%dtime*st%mb_source(j,i), 0.0_wp)
      end do
    end do

  endif

!-------- Reset new ice thickness if needed --------

if (par%margin==1) then

  do i=0, grd%IMAX
    do j=0, grd%JMAX
      if (st%maske(j,i) >= 2 .and. st%H_neu(j,i) > eps_H) then
        st%H_neu(j,i) = 0.0_wp
      end if
    end do
  end do

else   !/* MARGIN==2 or 3  */

if (par%l_abyss==1) then
  do i=0, grd%IMAX
    do j=0, grd%JMAX
      if (st%zl0(j,i) <= par%z_abyss) then
        st%H_neu(j,i) = 0.0_wp
      end if
    end do
  end do
else if (par%l_abyss==2) then
  do i=0, grd%IMAX
    do j=0, grd%JMAX
      if (st%zl(j,i) <= par%z_abyss) then
        st%H_neu(j,i) = 0.0_wp
      end if
    end do
  end do
end if

endif

end subroutine apply_mb_source

!-------------------------------------------------------------------------------
!> Adjustment of the newly computed ice thickness distribution.
!<------------------------------------------------------------------------------
  subroutine thk_adjust(st,tmr,par)

  implicit none

type(sico_state_class), intent(inout) :: st
type(sico_timer_class), intent(in) :: tmr
type(sico_par_class), intent(in) :: par

real(wp) :: target_topo_tau

!-------- Saving computed H_neu before any adjustments --------

  st%H_neu_tmp = st%H_neu

!-------- Correct negative thickness values --------

  where (st%H_neu < 0.0_wp) st%H_neu = 0.0_wp

  if (par%l_exclude_grl) then
    ! exclude update of ice thickness over Greenland
    where(st%id_mask.eq.2) 
      st%H_neu = st%H
    endwhere
  endif

!-------- Further adjustments --------

if (par%thk_evol==0) then

!  ------ No evolution of the ice thickness

  st%H_neu  = st%H   ! newly computed ice thickness is discarded

else if (par%thk_evol==1) then

!  ------ No adjustment, ice thickness evolves freely;
!         thus nothing to be done

else if (par%thk_evol==2) then

!  ------ Adjustment due to prescribed target topography

stop 'thk_evol==2 not implemented yet'

  if (tmr%time >= par%time_target_topo_final) then

     st%H_neu = st%H_target

  else if (tmr%time >= par%time_target_topo_init) then

     target_topo_tau = (par%time_target_topo_final-tmr%time)/3.0_wp

     st%H_neu =   ( target_topo_tau*st%H_neu + tmr%dtime*st%H_target ) &
             / ( target_topo_tau       + tmr%dtime )

  end if

else if (par%thk_evol==3) then

stop 'thk_evol==3 not implemented yet'

!  ------ Adjustment due to prescribed target topography with relaxation time

  st%H_neu =   ( target_topo_tau*st%H_neu + tmr%dtime*st%H_target ) &
          / ( par%target_topo_tau       + tmr%dtime )

else if (par%thk_evol==4) then

!  ------ Maximum ice extent constrained by prescribed mask

  where (st%mask_maxextent == 0) &   ! not allowed to glaciate
     st%H_neu = 0.0_wp

endif

!-------- Computation of the mass balance adjustment --------

  st%smb_corr = (st%H_neu-st%H_neu_tmp)*tmr%dtime_inv

  st%as_perp = st%as_perp + st%smb_corr
  st%accum   = st%accum   + max(st%smb_corr, 0.0_wp)
  st%runoff  = st%runoff  - min(st%smb_corr, 0.0_wp)
                    ! runoff is counted as positive for mass loss
    
end subroutine thk_adjust

!-------------------------------------------------------------------------------
!> Update of the ice-land-ocean mask etc.
!<------------------------------------------------------------------------------
subroutine calc_thk_mask_update(st,grd,tmr,par,n_calc_thk_mask_update_aux)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_timer_class), intent(in) :: tmr
type(sico_par_class), intent(in) :: par

integer, intent(in) :: n_calc_thk_mask_update_aux

integer                       :: i, j
real(wp)                      :: year_sec_inv

!-------- Term abbreviations --------

year_sec_inv = 1.0_wp/sec_year

!-------- Saving computed H_neu before any modifications --------

st%H_neu_tmp = st%H_neu

!-------- Update of the mask --------

if (n_calc_thk_mask_update_aux == 1) then
   call calc_thk_mask_update_aux1(st,grd,tmr,par)
else if (n_calc_thk_mask_update_aux == 2) then
   call calc_thk_mask_update_aux2(st,grd,tmr,par)
else if (n_calc_thk_mask_update_aux == 3) then
   call calc_thk_mask_update_aux3(st,grd,tmr,par)
else
   write(6, fmt='(a)') ' >>> calc_thk_mask_update:'
   write(6, fmt='(a)') '        n_calc_thk_mask_update_aux has no valid value!'
   stop
end if

!  ------ Enforce connectivity of the ocean

!call ocean_connect(st,grd)

!-------- Time derivatives --------

st%dzs_dtau  = (st%zs_neu-st%zs)*tmr%dtime_inv
st%dzb_dtau  = (st%zb_neu-st%zb)*tmr%dtime_inv
st%dzm_dtau  = st%dH_t_dtau+st%dzb_dtau
st%dH_c_dtau = st%dzs_dtau-st%dzm_dtau

if (par%thk_evol==2) then
if ( abs((tmr%time-par%time_target_topo_final)*year_sec_inv) < eps ) then
   st%dzs_dtau  = 0.0_wp   ! Introduced
   st%dzb_dtau  = 0.0_wp   ! by
   st%dzm_dtau  = 0.0_wp   ! Tatsuru Sato
   st%dH_c_dtau = 0.0_wp   ! for
   st%dH_t_dtau = 0.0_wp   ! stability reasons
end if
endif

!-------- New gradients --------

if (par%topograd==0) then
  call topograd_1(st,grd, 2)
else if (par%topograd==1) then
  call topograd_2(st,grd, 2)
endif

!-------- Compute the volume flux across the CTS, am_perp --------

if (par%calcmod==1) then

  do i=1, grd%IMAX-1
  do j=1, grd%JMAX-1
  
     if ( ((st%maske(j,i) == 0).or.(st%maske(j,i) == 3)) &
          .and.(st%n_cts(j,i) == 1)) then
  
        st%am_perp_st(j,i) = &
                  0.5_wp*(st%vx_c(0,j,i)+st%vx_c(0,j,i-1))*st%dzm_dxi_g(j,i) &
                + 0.5_wp*(st%vy_c(0,j,i)+st%vy_c(0,j-1,i))*st%dzm_deta_g(j,i) &
                - 0.5_wp*(st%vz_c(0,j,i)+st%vz_t(grd%KTMAX-1,j,i))
        st%am_perp(j,i) = st%am_perp_st(j,i) + st%dzm_dtau(j,i)
  
     else
        st%am_perp_st(j,i) = 0.0_wp
        st%am_perp(j,i)    = 0.0_wp
     end if
  
  end do
  end do

endif

end subroutine calc_thk_mask_update

!-------------------------------------------------------------------------------
!> Update of the ice-land-ocean mask for SIA-only dynamics of ice sheets
!! without ice shelves.
!<------------------------------------------------------------------------------
subroutine calc_thk_mask_update_aux1(st,grd,tmr,par)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_timer_class), intent(in) :: tmr
type(sico_par_class), intent(in) :: par

integer :: i, j
real(wp)     :: rhosw_rho_ratio
real(wp)     :: H_inv


!-------- Term abbreviations --------

rhosw_rho_ratio = RHO_SW/RHO

! ocean without contact to grounded ice should stay ocean
if (par%ice_control==1) then
if (par%thk_evol>=1) then
do i=1, grd%IMAX-1
do j=1, grd%JMAX-1
! ocean without contact to ice (grounded of shelf ice) should stay ocean
  if( .not.(st%maske(j,i)==1 .or. st%maske(j,i)==0 .or.   &
            (st%maske(j,i)==2 .and. & ! ocean before
             (st%maske(j,i-1)==0 .or. st%maske(j,i+1)==0 .or. &
              st%maske(j-1,i)==0 .or. st%maske(j+1,i)==0))) ) then
    st%maske(j,i) = 2
    st%H_neu(j,i) = 0.0_wp
  end if
end do
end do
end if
end if

!-------- Update of the mask --------

st%zs_neu    = st%zb_neu + st%H_neu   ! ice surface topography
st%H_sea_neu = st%z_sl - st%zl_neu    ! sea depth
st%H_balance = st%H_sea_neu*rhosw_rho_ratio

if (par%thk_evol>=1) then

do i=0, grd%IMAX
do j=0, grd%JMAX

   if (st%maske(j,i) <= 1) then   ! grounded ice or ice-free land

      if (st%H_neu(j,i) >= eps_H) then
         st%maske(j,i) = 0   ! grounded ice
      else
        ! MW, land or sea depending on sea level. Was this a bug? fixme
        if (st%zl_neu(j,i)>st%z_sl(j,i)) then
          st%maske(j,i) = 1   ! ice-free land
        else
          st%maske(j,i) = 2   ! sea
        endif
      end if

   else   ! (maske(j,i) == 2); sea

     if (par%margin==2) then

      if (st%H_neu(j,i) >= eps_H) then

         if ( st%H_neu(j,i) < (rhosw_rho_ratio*st%H_sea_neu(j,i)) ) then
           if (par%marine_ice_formation==1) then
             st%maske(j,i) = 2   ! floating ice cut off -> sea
           else if (par%marine_ice_formation==2) then
             st%maske(j,i) = 0   ! "underwater ice"
           endif
         else
            st%maske(j,i) = 0   ! grounded ice
         end if

         if (par%marine_ice_calving==2 .or. par%marine_ice_calving==4 .or. par%marine_ice_calving==6) then
           if (st%zl0(j,i) < st%z_mar) st%maske(j,i) = 2   ! sea
         else if (par%marine_ice_calving==3 .or. par%marine_ice_calving==5 .or. par%marine_ice_calving==7) then
           if (st%zl_neu(j,i) < st%z_mar) st%maske(j,i) = 2   ! sea
         endif

      else
         st%maske(j,i) = 2   ! sea
      end if

     endif

   end if

end do
end do

endif

!  ------ Adjustment due to prescribed target topography

if (par%thk_evol==2) then
if (tmr%time >= par%time_target_topo_final) st%maske = st%maske_target
endif

!-------- Correction of zs_neu and zb_neu for ice-free land and sea --------

do i=0, grd%IMAX
do j=0, grd%JMAX

   if (st%maske(j,i) == 1) then   ! ice-free land

      st%zs_neu(j,i) = st%zb_neu(j,i)   ! this prevents zs_neu(j,i)
      st%H_neu(j,i)  = 0.0_wp        ! from being below zb_neu(j,i)

   else if (st%maske(j,i) == 2) then   ! sea

      st%zs_neu(j,i) = st%zb_neu(j,i)   ! this prevents zs_neu(j,i)
      st%H_neu(j,i)  = 0.0_wp        ! from being below zb_neu(j,i)

   else if (st%maske(j,i) == 3) then   ! floating ice

      write(6, fmt='(a)') ' >>> calc_thk_mask_update_aux1:'
      write(6, fmt='(a)') '          maske(j,i)==3 not allowed for'
      write(6, fmt='(a)') '          SIA-only dynamics!'
      stop

   end if

end do
end do

!-------- Special scheme for alpine glaciation --------

if (par%l_alpine) then
  call alpine_ice(st,grd,tmr,par)
endif

!-------- Limit thickness of isolated ice points --------

call limit_thickness_isolated_ice(st,grd,par)

!-------- Computation of further quantities --------

do i=0, grd%IMAX
do j=0, grd%JMAX

   if (st%maske(j,i) == 0) then   ! grounded ice

      if (st%n_cts(j,i) == 1) then
         if (st%H(j,i) > 0.0_wp) then
            H_inv        = 1.0_wp/st%H(j,i)
            st%H_c_neu(j,i) = st%H_c(j,i) * st%H_neu(j,i)*H_inv
            st%H_t_neu(j,i) = st%H_t(j,i) * st%H_neu(j,i)*H_inv
         else
            st%H_c_neu(j,i) = 0.99_wp * st%H_neu(j,i)   ! this case should not occur,
            st%H_t_neu(j,i) = 0.01_wp * st%H_neu(j,i)   ! just for safety
         end if
         st%zm_neu(j,i) = st%zb_neu(j,i)+st%H_t_neu(j,i)
      else
         st%H_c_neu(j,i) = st%H_neu(j,i)
         st%H_t_neu(j,i) = 0.0_wp
         st%zm_neu(j,i)  = st%zb_neu(j,i)
      end if

   else   ! maske(j,i) == 1 or 2, ice-free land or sea

      st%H_c_neu(j,i) = 0.0_wp
      st%H_t_neu(j,i) = 0.0_wp
      st%H_neu(j,i)   = 0.0_wp
      st%zs_neu(j,i)  = st%zb_neu(j,i)
      st%zm_neu(j,i)  = st%zb_neu(j,i)

   end if

end do
end do

end subroutine calc_thk_mask_update_aux1

!-------------------------------------------------------------------------------
!> Update of the ice-land-ocean mask for hybrid SIA/SStA dynamics of ice sheets
!! without ice shelves.
!<------------------------------------------------------------------------------
subroutine calc_thk_mask_update_aux2(st,grd,tmr,par)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_timer_class), intent(in) :: tmr
type(sico_par_class), intent(in) :: par

integer :: i, j
real(wp)     :: rhosw_rho_ratio
real(wp)     :: H_inv


!-------- Term abbreviations --------

rhosw_rho_ratio = RHO_SW/RHO

! ocean without contact to grounded ice should stay ocean
if (par%ice_control==1) then
if (par%thk_evol>=1) then
do i=1, grd%IMAX-1
do j=1, grd%JMAX-1
! ocean without contact to ice (grounded of shelf ice) should stay ocean
  if( .not.(st%maske(j,i)==1 .or. st%maske(j,i)==0 .or.   &
            (st%maske(j,i)==2 .and. & ! ocean before
             (st%maske(j,i-1)==0 .or. st%maske(j,i+1)==0 .or. &
              st%maske(j-1,i)==0 .or. st%maske(j+1,i)==0))) ) then
    st%maske(j,i) = 2
    st%H_neu(j,i) = 0.0_wp
  end if
end do
end do
end if
end if

!-------- Update of the mask --------

st%zs_neu    = st%zb_neu + st%H_neu   ! ice surface topography
st%H_sea_neu = st%z_sl - st%zl_neu    ! sea depth
st%H_balance = st%H_sea_neu*rhosw_rho_ratio

if (par%thk_evol>=1) then

do i=0, grd%IMAX
do j=0, grd%JMAX

   if (st%maske(j,i) <= 1) then   ! grounded ice or ice-free land

      if (st%H_neu(j,i) >= eps_H) then
         st%maske(j,i) = 0   ! grounded ice
      else
        ! MW, land or sea depending on sea level. Was this a bug? fixme
        if (st%zl_neu(j,i)>st%z_sl(j,i)) then
          st%maske(j,i) = 1   ! ice-free land
        else
          st%maske(j,i) = 2   ! sea
        endif
      end if

   else   ! (st%maske(j,i) == 2); sea

     if (par%margin==2) then

      if (st%H_neu(j,i) >= eps_H) then

         if ( st%H_neu(j,i) < (rhosw_rho_ratio*st%H_sea_neu(j,i)) ) then
           if (par%marine_ice_formation==1) then
             st%maske(j,i) = 2   ! floating ice cut off -> sea
           else if (par%marine_ice_formation==2) then
             st%maske(j,i) = 0   ! "underwater ice"
           endif
         else
            st%maske(j,i) = 0   ! grounded ice
         end if

         if (par%marine_ice_calving==2 .or. par%marine_ice_calving==4 .or. par%marine_ice_calving==6) then
           if (st%zl0(j,i) < st%z_mar) st%maske(j,i) = 2   ! sea
         else if (par%marine_ice_calving==3 .or. par%marine_ice_calving==5 .or. par%marine_ice_calving==7) then
           if (st%zl_neu(j,i) < st%z_mar) st%maske(j,i) = 2   ! sea
         endif

      else
         st%maske(j,i) = 2   ! sea
      end if

    endif

   end if

end do
end do

endif

!  ------ Adjustment due to prescribed target topography

if (par%thk_evol==2) then
if (tmr%time >= par%time_target_topo_final) st%maske = st%maske_target
endif

!-------- Correction of zs_neu and zb_neu for ice-free land and sea --------

do i=0, grd%IMAX
do j=0, grd%JMAX

   if (st%maske(j,i) == 1) then   ! ice-free land

      st%zs_neu(j,i) = st%zb_neu(j,i)   ! this prevents zs_neu(j,i)
      st%H_neu(j,i)  = 0.0_wp        ! from being below zb_neu(j,i)

   else if (st%maske(j,i) == 2) then   ! sea

      st%zs_neu(j,i) = st%zb_neu(j,i)   ! this prevents zs_neu(j,i)
      st%H_neu(j,i)  = 0.0_wp        ! from being below zb_neu(j,i)

   else if (st%maske(j,i) == 3) then   ! floating ice

      write(6, fmt='(a)') ' >>> calc_thk_mask_update_aux2:'
      write(6, fmt='(a)') '          maske(j,i)==3 not allowed for'
      write(6, fmt='(a)') '          hybrid SIA/SStA dynamics of ice sheets'
      write(6, fmt='(a)') '          without ice shelves!'
      stop

   end if

end do
end do

!-------- Special scheme for alpine glaciation --------

if (par%l_alpine) then
  call alpine_ice(st,grd,tmr,par)
endif

!-------- Limit thickness of isolated ice points --------

call limit_thickness_isolated_ice(st,grd,par)

!-------- Computation of further quantities --------

do i=0, grd%IMAX
do j=0, grd%JMAX

   if (st%maske(j,i) == 0) then   ! grounded ice

      if (st%n_cts(j,i) == 1) then
         if (st%H(j,i) > 0.0_wp) then
            H_inv        = 1.0_wp/st%H(j,i)
            st%H_c_neu(j,i) = st%H_c(j,i) * st%H_neu(j,i)*H_inv
            st%H_t_neu(j,i) = st%H_t(j,i) * st%H_neu(j,i)*H_inv
         else
            st%H_c_neu(j,i) = 0.99_wp * st%H_neu(j,i)   ! this case should not occur,
            st%H_t_neu(j,i) = 0.01_wp * st%H_neu(j,i)   ! just for safety
         end if
         st%zm_neu(j,i) = st%zb_neu(j,i)+st%H_t_neu(j,i)
      else
         st%H_c_neu(j,i) = st%H_neu(j,i)
         st%H_t_neu(j,i) = 0.0_wp
         st%zm_neu(j,i)  = st%zb_neu(j,i)
      end if

   else   ! maske(j,i) == 1 or 2, ice-free land or sea

      st%H_c_neu(j,i) = 0.0_wp
      st%H_t_neu(j,i) = 0.0_wp
      st%H_neu(j,i)   = 0.0_wp
      st%zs_neu(j,i) = st%zb_neu(j,i)
      st%zm_neu(j,i) = st%zb_neu(j,i)

   end if

end do
end do

end subroutine calc_thk_mask_update_aux2

!-------------------------------------------------------------------------------
!> Update of the ice-land-ocean mask for coupled SIA/SSA or
!! SIA/SStA/SSA dynamics of ice sheets with ice shelves.
!<------------------------------------------------------------------------------
subroutine calc_thk_mask_update_aux3(st,grd,tmr,par)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_timer_class), intent(in) :: tmr
type(sico_par_class), intent(in) :: par

integer :: i, j
real(wp)     :: dtime_tau_calv, rhosw_rho_ratio, rho_rhosw_ratio
real(wp)     :: H_inv
logical      :: flag_calving_event


!-------- Term abbreviations --------

dtime_tau_calv  = tmr%dtime/max(par%tau_calv,eps)
rhosw_rho_ratio = RHO_SW/RHO
rho_rhosw_ratio = RHO/RHO_SW



! ocean without contact to ice (grounded or shelf ice) should stay ocean
if (par%ice_control==1) then
if (par%thk_evol>=1) then
do i=1, grd%IMAX-1
do j=1, grd%JMAX-1
! ocean without contact to ice (grounded of shelf ice) should stay ocean
  if( .not.(st%maske(j,i)==1 .or. st%maske(j,i)==0 .or. st%maske(j,i)==3 .or.   &
            (st%maske(j,i)==2 .and. & ! ocean before
             (st%maske(j,i-1)==0 .or. st%maske(j,i-1)==3 .or. &
              st%maske(j,i+1)==0 .or. st%maske(j,i+1)==3 .or. &
              st%maske(j-1,i)==0 .or. st%maske(j-1,i)==3 .or. &
              st%maske(j+1,i)==0 .or. st%maske(j+1,i)==3))) ) then
    st%maske(j,i) = 2
    st%H_neu(j,i) = 0.0_wp
  end if
end do
end do
end if
end if

!-------- Update of the mask --------

st%zs_neu    = st%zb_neu + st%H_neu   ! ice surface topography
st%H_sea_neu = st%z_sl - st%zl_neu    ! sea depth
st%H_balance = st%H_sea_neu*rhosw_rho_ratio

if (par%thk_evol>=1) then

do i=1, grd%IMAX-1
do j=1, grd%JMAX-1

   ! grounding_line migration check

   if ( ( st%maske(j,i) <= 1 &
         .and.    (st%maske(j,i+1)>1.or.st%maske(j,i-1)>1 &
               .or.st%maske(j+1,i)>1.or.st%maske(j-1,i)>1) ) &
        .or. &
        ( st%maske(j,i)>=2 &
         .and.    (st%maske(j,i+1)==0.or.st%maske(j,i-1)==0  &
               .or.st%maske(j+1,i)==0.or.st%maske(j-1,i)==0) ) ) then

      if (st%H_neu(j,i) >= eps_H) then

         if (st%H_neu(j,i)<st%H_balance(j,i).and.st%zl_neu(j,i)<st%z_sl(j,i)) then
            st%maske(j,i)    = 3
            st%zb_neu(j,i)   = st%z_sl(j,i)-rho_rhosw_ratio*st%H_neu(j,i)
            st%dzb_dtau(j,i) = tmr%dtime_inv*(st%zb_neu(j,i)-st%zb(j,i))
            st%zs_neu(j,i)   = st%zb_neu(j,i)+st%H_neu(j,i)
         else
            st%maske(j,i)    = 0
            st%zb_neu(j,i)   = st%zl_neu(j,i)
            st%dzb_dtau(j,i) = st%dzl_dtau(j,i)
            st%zs_neu(j,i)   = st%zb_neu(j,i)+st%H_neu(j,i)
         end if

      else   ! if (st%H_neu(j,i) <= eps_H then

         if (st%zl_neu(j,i)>st%z_sl(j,i)) then
            st%maske(j,i)    = 1
            st%zb_neu(j,i)   = st%zl_neu(j,i)
            st%dzb_dtau(j,i) = tmr%dtime_inv*(st%zb_neu(j,i)-st%zb(j,i))
            st%zs_neu(j,i)   = st%zb_neu(j,i)
         else
            st%maske(j,i)    = 2
            st%zb_neu(j,i)   = st%z_sl(j,i)
            st%dzb_dtau(j,i) = 0.0_wp
            st%zs_neu(j,i)   = st%z_sl(j,i)
         end if

      end if

   else if (st%maske(j,i) <= 1) then   ! grounded ice or ice-free land

      if (st%H_neu(j,i) >= eps_H) then   ! can change maske 0 or 1
         st%maske(j,i)    = 0
         st%zb_neu(j,i)   = st%zl_neu(j,i)
         st%dzb_dtau(j,i) = st%dzl_dtau(j,i)
         st%zs_neu(j,i)   = st%zb_neu(j,i)+st%H_neu(j,i)
      else
         ! MW bug? glaciated areas don't turn to ocean when ice retreating! fixme
         !st%maske(j,i)    = 1
         !st%zb_neu(j,i)   = st%zl_neu(j,i)
         !st%dzb_dtau(j,i) = st%dzl_dtau(j,i)
         !st%zs_neu(j,i)   = st%zl_neu(j,i)
         if (st%zl_neu(j,i)>st%z_sl(j,i)) then
           st%maske(j,i)    = 1
           st%zb_neu(j,i)   = st%zl_neu(j,i)
           st%dzb_dtau(j,i) = st%dzl_dtau(j,i)
           st%zs_neu(j,i)   = st%zl_neu(j,i)
         else
           st%maske(j,i)    = 2
           st%zb_neu(j,i)   = st%z_sl(j,i)
           st%dzb_dtau(j,i) = 0.0_wp
           st%zs_neu(j,i)   = st%z_sl(j,i)
         end if
      end if

   else   ! if (maske(j,i)==2.or.maske(j,i)==3) then

      if (st%H_neu(j,i) > eps_H) then

         if (st%H_neu(j,i)<st%H_balance(j,i).and.st%zl_neu(j,i)<st%z_sl(j,i)) then
            st%maske(j,i)    = 3
            st%zb_neu(j,i)   = st%z_sl(j,i)-rho_rhosw_ratio*st%H_neu(j,i)
            st%dzb_dtau(j,i) = tmr%dtime_inv*(st%zb_neu(j,i)-st%zb(j,i))
            st%zs_neu(j,i)   = st%zb_neu(j,i)+st%H_neu(j,i)
         else
            st%maske(j,i)    = 0
            st%zb_neu(j,i)   = st%zl_neu(j,i)
            st%dzb_dtau(j,i) = st%dzl_dtau(j,i)
            st%zs_neu(j,i)   = st%zb_neu(j,i)+st%H_neu(j,i)
         end if

      else   ! if (H_neu(j,i) <= eps_H) then

         if (st%zl_neu(j,i)>st%z_sl(j,i)) then
            st%maske(j,i)    = 1
            st%zb_neu(j,i)   = st%zl_neu(j,i)
            st%dzb_dtau(j,i) = tmr%dtime_inv*(st%zb_neu(j,i)-st%zb(j,i))
            st%zs_neu(j,i)   = st%zb_neu(j,i)
         else
            st%maske(j,i)    = 2
            st%zb_neu(j,i)   = st%z_sl(j,i)
            st%dzb_dtau(j,i) = 0.0_wp
            st%zs_neu(j,i)   = st%z_sl(j,i)
         end if

      end if

   end if

end do
end do

if (par%ice_shelf_calving==2) then

do i=0, grd%IMAX
do j=0, grd%JMAX
   if(st%zl_fil(j,i).le.par%zl_deep) then 
     st%H_calv(j,i) = par%h_calv_deep
   else if(st%zl_fil(j,i).ge.par%zl_shallow) then
     st%H_calv(j,i) = par%h_calv_shallow
   else
     st%H_calv(j,i) = par%h_calv_deep + (st%zl_fil(j,i)-par%zl_deep)* & 
              (par%h_calv_deep-par%h_calv_shallow)/(par%zl_deep-par%zl_shallow)
   end if
end do
end do

do

   st%flag_calving_front_1 = .false.
   flag_calving_event   = .false.

   do i=1, grd%IMAX-1
   do j=1, grd%JMAX-1

      if ( (st%maske(j,i)==3) &   ! floating ice
           .and. &
             (    (st%maske(j,i+1)==2)   &   ! with
              .or.(st%maske(j,i-1)==2)   &   ! one
              .or.(st%maske(j+1,i)==2)   &   ! neighbouring
              .or.(st%maske(j-1,i)==2) ) &   ! sea point
         ) &
         st%flag_calving_front_1(j,i) = .true.   ! preliminary detection
                                              ! of the calving front

      ! new calving parameterisation by Reinhard
      if ( st%flag_calving_front_1(j,i).and.st%h(j,i)<st%H_calv(j,i) &  !fixme: st%h_neu(j,i), not st%h(j,i)  
           .and. &
              (    (st%h_neu(j,i+1)<st%H_calv(j,i))   &   ! with
              .and.(st%h_neu(j,i-1)<st%H_calv(j,i))   &   ! all
              .and.(st%h_neu(j+1,i)<st%H_calv(j,i))   &   ! neighbouring
              .and.(st%h_neu(j-1,i)<st%H_calv(j,i)) ) &   ! too thin ice points 
         ) then
                                              
!      if ( st%flag_calving_front_1(j,i).and.(st%h_neu(j,i) < st%H_calv(j,i)) ) then
         flag_calving_event = .true.  ! calving event,
         st%maske(j,i)      = 2   ! floating ice point changes to sea point
      end if

   end do
   end do

   if (.not.flag_calving_event) exit

end do

else if (par%ice_shelf_calving==3) then

do i=0, grd%IMAX
do j=0, grd%JMAX
   if(st%zl_fil(j,i).le.par%zl_deep) then 
     st%H_calv(j,i) = par%h_calv_deep
   else if(st%zl_fil(j,i).ge.par%zl_shallow) then
     st%H_calv(j,i) = par%h_calv_shallow
   else
     st%H_calv(j,i) = par%h_calv_deep + (st%zl_fil(j,i)-par%zl_deep)* & 
              (par%h_calv_deep-par%h_calv_shallow)/(par%zl_deep-par%zl_shallow)
   end if
end do
end do

   st%flag_calving_front_1 = .false.
   do i=1, grd%IMAX-1
   do j=1, grd%JMAX-1

      if ( (st%maske(j,i)==3) &   ! floating ice
           .and. &
             (    (st%maske(j,i+1)==2)   &   ! with
              .or.(st%maske(j,i-1)==2)   &   ! one
              .or.(st%maske(j+1,i)==2)   &   ! neighbouring
              .or.(st%maske(j-1,i)==2) ) &   ! sea point
         ) &
         st%flag_calving_front_1(j,i) = .true.   ! detection calving front

      ! new calving parameterisation by Reinhard
      if ( st%flag_calving_front_1(j,i).and.st%h_neu(j,i)<st%H_calv(j,i) &
           .and. &
              (    (st%h_neu(j,i+1)<st%H_calv(j,i))   &   ! with
              .and.(st%h_neu(j,i-1)<st%H_calv(j,i))   &   ! all
              .and.(st%h_neu(j+1,i)<st%H_calv(j,i))   &   ! neighbouring
              .and.(st%h_neu(j-1,i)<st%H_calv(j,i)) ) &   ! too thin ice points 
         ) then
         st%h_neu(j,i) = max(st%h_neu(j,i) - (st%H_calv(j,i)-st%h_neu(j,i))*dtime_tau_calv, 0.0_wp)
         if(st%H_neu(j,i) <= eps_H) st%maske(j,i) = 2
         ! this can get thinner only, i.e. maske= 3 or 2
      end if

   end do
   end do

else if (par%ice_shelf_calving==4) then

do i=0, grd%IMAX
do j=0, grd%JMAX

   if (st%maske(j,i)==3) st%maske(j,i) = 2
      ! float-kill: all floating ice points changed to sea points

end do
end do

else if (par%ice_shelf_calving==9) then

#if (defined(MISMIPP))

do i=0, grd%IMAX
do j=0, grd%JMAX

   if ((st%maske(j,i)==3).and.(xi(i) > Lx)) then
      st%maske(j,i) = 2   ! floating ice point changes to sea point
   end if

end do
end do

#else
write(6, fmt='(a)') ' >>> calc_thk_mask_update_aux3:'
write(6, fmt='(a)') '          Option ICE_SHELF_CALVING==9'
write(6, fmt='(a)') '          only defined for MISMIP+!'
stop
#endif

endif

endif

!  ------ Adjustment due to prescribed target topography

if (par%thk_evol==2) then
if (tmr%time >= par%time_target_topo_final) st%maske = st%maske_target
endif

!-------- Correction of zs_neu and zb_neu
!         for ice-free land, sea and floating ice --------

do i=0, grd%IMAX
do j=0, grd%JMAX

   if (st%maske(j,i) == 1) then   ! ice-free land

      st%zs_neu(j,i) = st%zb_neu(j,i)   ! this prevents zs_neu(j,i)
      st%H_neu(j,i)  = 0.0_wp        ! from being below zb_neu(j,i)

   else if (st%maske(j,i) == 2) then   ! sea

      st%zs_neu(j,i) = st%z_sl(j,i)          ! both zs_neu(j,i) and zb_neu(j,i)
      st%zb_neu(j,i) = st%z_sl(j,i)          ! set to the current sea level stand
      st%H_neu(j,i)  = 0.0_wp

   else if (st%maske(j,i) == 3) then   ! floating ice

      st%H_neu(j,i) = st%zs_neu(j,i)-st%zb_neu(j,i)   ! ice thickness

      st%zb_neu(j,i) = st%z_sl(j,i) - rho_rhosw_ratio*st%H_neu(j,i)   ! floating condition
      st%zs_neu(j,i) = st%zb_neu(j,i) + st%H_neu(j,i)

   end if

end do
end do

!-------- Special scheme for alpine glaciation --------

if (par%l_alpine) then
  call alpine_ice(st,grd,tmr,par)
endif

!-------- Limit thickness of isolated ice points --------

call limit_thickness_isolated_ice(st,grd,par)

!-------- Computation of further quantities --------

do i=0, grd%IMAX
do j=0, grd%JMAX

   if ( (st%maske(j,i) == 0).or.(st%maske(j,i) == 3) ) then
                                        ! grounded or floating ice

      if (st%n_cts(j,i) == 1) then
         if (st%H(j,i) > 0.0_wp) then
            H_inv        = 1.0_wp/st%H(j,i)
            st%H_c_neu(j,i) = st%H_c(j,i) * st%H_neu(j,i)*H_inv
            st%H_t_neu(j,i) = st%H_t(j,i) * st%H_neu(j,i)*H_inv
         else
            st%H_c_neu(j,i) = 0.99_wp * st%H_neu(j,i)   ! this case should not occur,
            st%H_t_neu(j,i) = 0.01_wp * st%H_neu(j,i)   ! just for safety
         end if
         st%zm_neu(j,i) = st%zb_neu(j,i)+st%H_t_neu(j,i)
      else
         st%H_c_neu(j,i) = st%H_neu(j,i)
         st%H_t_neu(j,i) = 0.0_wp
         st%zm_neu(j,i)  = st%zb_neu(j,i)
      end if

   else   ! maske(j,i) == 1 or 2, ice-free land or sea

      st%H_c_neu(j,i) = 0.0_wp
      st%H_t_neu(j,i) = 0.0_wp
      st%H_neu(j,i)   = 0.0_wp
      st%zs_neu(j,i) = st%zb_neu(j,i)
      st%zm_neu(j,i) = st%zb_neu(j,i)

   end if

end do
end do

end subroutine calc_thk_mask_update_aux3


!-------------------------------------------------------------------------------
!> Limit thickness of isolated ice points.
!<------------------------------------------------------------------------------
subroutine limit_thickness_isolated_ice(st,grd,par)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_par_class), intent(in) :: par

integer :: i, j
real(dp)     :: rho_rhosw_ratio
real(dp)     :: H_isol_max
logical      :: flag_kill_isolated_ice

rho_rhosw_ratio = RHO/RHO_SW

H_isol_max = par%H_isol_max
if (par%H_isol_max < eps_H) then
   H_isol_max = 0.0_wp
   flag_kill_isolated_ice = .true.
else
   flag_kill_isolated_ice = .false.
end if

do i=1, grd%IMAX-1
do j=1, grd%JMAX-1

   if (st%maske(j,i) == 0) then   ! grounded ice

      if ( ((st%maske(j,i+1) == 1).or.(st%maske(j,i+1) == 2)) &
           .and. ((st%maske(j,i-1) == 1).or.(st%maske(j,i-1) == 2)) &
           .and. ((st%maske(j+1,i) == 1).or.(st%maske(j+1,i) == 2)) &
           .and. ((st%maske(j-1,i) == 1).or.(st%maske(j-1,i) == 2)) &
         ) then   ! all nearest neighbours ice-free

         st%H_neu(j,i)  = min(st%H_neu(j,i), H_isol_max)
         st%zs_neu(j,i) = st%zb_neu(j,i)+st%H_neu(j,i)

         if (flag_kill_isolated_ice) then
            if (st%zb_neu(j,i) >= st%z_sl(j,i)) then
               st%maske(j,i) = 1
            else
               st%maske(j,i) = 2
            end if
         end if

      end if

   else if (st%maske(j,i) == 3) then   ! floating ice

      if ( ((st%maske(j,i+1) == 1).or.(st%maske(j,i+1) == 2)) &
           .and. ((st%maske(j,i-1) == 1).or.(st%maske(j,i-1) == 2)) &
           .and. ((st%maske(j+1,i) == 1).or.(st%maske(j+1,i) == 2)) &
           .and. ((st%maske(j-1,i) == 1).or.(st%maske(j-1,i) == 2)) &
         ) then   ! all nearest neighbours ice-free

         st%H_neu(j,i)  = min(st%H_neu(j,i), H_isol_max)
         st%zb_neu(j,i) = st%z_sl(j,i)-rho_rhosw_ratio*st%H_neu(j,i)
         st%zs_neu(j,i) = st%zb_neu(j,i)+st%H_neu(j,i)

         if (flag_kill_isolated_ice) then
            st%maske(j,i) = 2
         end if

      end if

   end if

end do
end do

end subroutine limit_thickness_isolated_ice

!-------------------------------------------------------------------------------
!> Special scheme for alpine glaciation 
!<------------------------------------------------------------------------------
subroutine alpine_ice(st,grd,tmr,par)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_timer_class), intent(in) :: tmr
type(sico_par_class), intent(in) :: par

integer :: i, j, in_0, jn_0
!real(wp), parameter :: H_thresh = 10._wp
real(wp) :: H_thresh
real(wp), parameter :: zs_std_min = 200._wp
real(wp), parameter :: zs_std_max = 1000._wp
real(wp), parameter :: H_thresh_min = 10._wp
real(wp), parameter :: H_thresh_max = 1000._wp
logical :: is_neighbor, is_neighbor_im, is_neighbor_ip, is_neighbor_jm, is_neighbor_jp
real(wp) :: zs_n_min, f_redi
real(wp), dimension(:,:), allocatable :: zs_neu_2


  allocate(zs_neu_2(0:grd%JMAX,0:grd%IMAX))
  zs_neu_2(:,:) = st%zs_neu(:,:)

  do i=1, grd%IMAX-1
    do j=1, grd%JMAX-1

      if( (st%maske(j,i).eq.0 .and. st%as_perp(j,i).gt.0._wp) .and. st%zs_std(j,i).gt.zs_std_min ) then

        H_thresh = H_thresh_min*(((zs_std_max-zs_std_min)/(st%zs_std(j,i)-zs_std_min)))
        H_thresh = min(H_thresh_max,H_thresh)
        H_thresh = max(H_thresh_min,H_thresh)

        is_neighbor=.false.
        is_neighbor_im=.false.
        is_neighbor_ip=.false.
        is_neighbor_jm=.false.
        is_neighbor_jp=.false.

        if(st%zs_neu(j,i-1).lt.st%zs_neu(j,i) .and. st%zs_neu(j,i-1)-st%zb_neu(j,i-1).lt.H_thresh) then
          is_neighbor_im=.true.
          is_neighbor=.true.
        end if
        if(st%zs_neu(j,i+1).lt.st%zs_neu(j,i) .and. st%zs_neu(j,i+1)-st%zb_neu(j,i+1).lt.H_thresh) then
          is_neighbor_ip=.true.
          is_neighbor=.true.
        end if
        if(st%zs_neu(j-1,i).lt.st%zs_neu(j,i) .and. st%zs_neu(j-1,i)-st%zb_neu(j-1,i).lt.H_thresh) then
          is_neighbor_jm=.true.
          is_neighbor=.true.
        end if
        if(st%zs_neu(j+1,i).lt.st%zs_neu(j,i) .and. st%zs_neu(j+1,i)-st%zb_neu(j+1,i).lt.H_thresh) then
          is_neighbor_jp=.true.
          is_neighbor=.true.
        end if

        zs_n_min=1.e10_wp
        if(is_neighbor_im .and. st%zs_neu(j,i-1).le.zs_n_min) then
          zs_n_min=st%zs_neu(j,i-1)
          in_0=i-1
          jn_0=j
        end if
        if(is_neighbor_ip .and. st%zs_neu(j,i+1).le.zs_n_min) then
          zs_n_min=st%zs_neu(j,i+1)
          in_0=i+1
          jn_0=j
        end if
        if(is_neighbor_jm .and. st%zs_neu(j-1,i).le.zs_n_min) then
          zs_n_min=st%zs_neu(j-1,i)
          in_0=i
          jn_0=j-1
        end if
        if(is_neighbor_jp .and. st%zs_neu(j+1,i).le.zs_n_min) then
          zs_n_min=st%zs_neu(j+1,i)
          in_0=i
          jn_0=j+1
        end if

        if(is_neighbor) then
          f_redi=min((st%zs_neu(j,i)-st%zb_neu(j,i))/H_thresh,1.1_wp)
          ! It will be taken maximal as much ice away as there is
          f_redi=min(f_redi,(st%zs_neu(j,i)-st%zb_neu(j,i))/(st%as_perp(j,i)*tmr%dtime))
          zs_neu_2(j,i)       = st%zs_neu(j,i)       - f_redi*st%as_perp(j,i)*tmr%dtime
          zs_neu_2(jn_0,in_0) = st%zs_neu(jn_0,in_0) + f_redi*st%as_perp(j,i)*tmr%dtime
        end if
      end if
    end do
  end do

  ! updating of the ice surface and forbid negative thickness 

  do i=0, grd%IMAX
    do j=0, grd%JMAX
      st%zs_neu(j,i) = zs_neu_2(j,i)
      if(st%zs_neu(j,i).lt.st%zb_neu(j,i)) st%zs_neu(j,i) = st%zb_neu(j,i)
      st%H_neu(j,i) = st%zs_neu(j,i)-st%zb_neu(j,i)
    end do
  end do

  deallocate(zs_neu_2) 


end subroutine alpine_ice 


!-------------------------------------------------------------------------------
!> Enforce connectivity of the ocean.
!<------------------------------------------------------------------------------
subroutine ocean_connect(st,grd)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd

integer                           :: i, j
integer, dimension(:,:), allocatable :: mask_connect
integer, dimension(:,:), allocatable :: mask_connect_save, mask_connect_diff
logical                                :: flag_change


allocate(mask_connect(0:grd%JMAX,0:grd%IMAX))
allocate(mask_connect_save(0:grd%JMAX,0:grd%IMAX))
allocate(mask_connect_diff(0:grd%JMAX,0:grd%IMAX))

!-------- Determine connected area allowed to be ocean --------

mask_connect = 0

mask_connect(0:1                 , :                  ) = 1
mask_connect(grd%JMAX-1:grd%JMAX , :                  ) = 1
mask_connect(:                   , 0:1                ) = 1
mask_connect(:                   , grd%IMAX-1:grd%IMAX) = 1

flag_change = .true.

do while (flag_change)

   mask_connect_save = mask_connect

   do i=1, grd%IMAX-1
   do j=1, grd%JMAX-1

      if (mask_connect_save(j,i) == 1) then
         if (st%maske(j  ,i+1) >= 2) mask_connect(j  ,i+1) = 1
         if (st%maske(j  ,i-1) >= 2) mask_connect(j  ,i-1) = 1
         if (st%maske(j+1,i  ) >= 2) mask_connect(j+1,i  ) = 1
         if (st%maske(j-1,i  ) >= 2) mask_connect(j-1,i  ) = 1
         if (st%maske(j+1,i+1) >= 2) mask_connect(j+1,i+1) = 1
         if (st%maske(j+1,i-1) >= 2) mask_connect(j+1,i-1) = 1
         if (st%maske(j-1,i+1) >= 2) mask_connect(j-1,i+1) = 1
         if (st%maske(j-1,i-1) >= 2) mask_connect(j-1,i-1) = 1
      end if

   end do
   end do

   mask_connect_diff = abs(mask_connect-mask_connect_save)

   if (maxval(mask_connect_diff) > 0) then
      flag_change = .true.
   else
      flag_change = .false.
   end if

end do

!-------- Reset disconnected "ocean islands" to ice-free land --------

where ((st%maske == 2).and.(mask_connect == 0)) st%maske = 1

deallocate(mask_connect)
deallocate(mask_connect_save)
deallocate(mask_connect_diff)


end subroutine ocean_connect


!-------------------------------------------------------------------------------
!> Determination of the several components of the mass balance:
!! Runoff, calving, basal melt. 
!<------------------------------------------------------------------------------
subroutine account_mb_source(st,grd,tmr)

  implicit none

  type(sico_state_class), intent(inout) :: st
  type(sico_grid_class), intent(in) :: grd
  type(sico_timer_class), intent(in) :: tmr

  integer :: i, j


!-------- Compute mass balance components --------

  st%mb_source_apl      = 0.0_wp
  st%as_perp_apl        = 0.0_wp
  st%accum_apl          = 0.0_wp
  st%runoff_apl         = 0.0_wp
  st%Q_b_apl            = 0.0_wp
  st%calving_apl        = 0.0_wp
  st%mask_ablation_type = 0

  do i=0, grd%IMAX
  do j=0, grd%JMAX

     if (grd%flag_inner_point(j,i)) then   ! inner point

        if (st%H_neu(j,i) >= eps_H) then   ! glaciated point

           st%mb_source_apl(j,i) = (st%H_neu(j,i)-st%H_neu_flow(j,i))*tmr%dtime_inv
                                ! applied MB, based on ice thickness change

           ! Store melt quantities here for later accounting
           st%accum_apl(j,i)   = st%accum(j,i)
           st%runoff_apl(j,i)  = st%runoff(j,i)
           st%Q_b_apl(j,i)     = st%Q_b_tot(j,i)
           st%calving_apl(j,i) = st%calving(j,i)
           st%as_perp_apl(j,i) = st%accum_apl(j,i) - st%runoff_apl(j,i)

           ! Store melting mask value. This became an ice point,
           ! should be either grounded or floating ice.
           if (st%maske(j,i)==0) then
              st%mask_ablation_type(j,i) = 1 ! grounded ice
           else if (st%maske(j,i)==3) then
              st%mask_ablation_type(j,i) = 3 ! floating ice
           else
              st%mask_ablation_type(j,i) = 9 ! misaccounted (if appears)
           end if

        else if (st%H_neu(j,i) < eps_H .and. st%H_neu_flow(j,i) >= eps_H) then
                                           ! ice disappeared

           st%mb_source_apl(j,i) = (st%H_neu(j,i)-st%H_neu_flow(j,i))*tmr%dtime_inv
                                ! applied MB, based on ice thickness change

           st%accum_apl(j,i) = st%accum(j,i)

           if (st%maske(j,i)==2) then
              ! all mass that flows into the ocean counted as
              ! (large-scale) calving
              st%calving_apl(j,i) = -(st%mb_source_apl(j,i)-st%accum_apl(j,i))

           else if (st%runoff(j,i) > 0.0_wp .or. st%calving(j,i) > 0.0_wp &
                                         .or. st%Q_b_tot(j,i) > 0.0_wp) then
              ! land or shelf where ice disappered,
              ! three melting contingents shared proportionally
              st%runoff_apl(j,i)  = -(st%mb_source_apl(j,i)-st%accum_apl(j,i)) &
                                       *st%runoff(j,i) &
                                       /(st%runoff(j,i)+st%calving(j,i) &
                                                    +st%Q_b_tot(j,i))
              st%calving_apl(j,i) = -(st%mb_source_apl(j,i)-st%accum_apl(j,i)) &
                                       *st%calving(j,i) &
                                       /(st%runoff(j,i)+st%calving(j,i) &
                                                    +st%Q_b_tot(j,i))
              st%Q_b_apl(j,i)     = -(st%mb_source_apl(j,i)-st%accum_apl(j,i)) &
                                       *st%Q_b_tot(j,i) &
                                       /(st%runoff(j,i)+st%calving(j,i) &
                                                    +st%Q_b_tot(j,i))

           !!! else st%runoff(j,i)=0.0_wp .and. st%calving(j,i) = 0.0_wp &
           !!!                         .and. st%Q_b_tot(j,i) = 0.0_wp

           end if

           st%as_perp_apl(j,i) = st%accum_apl(j,i) - st%runoff_apl(j,i)

           ! Store melting mask value.
           ! This grid point is ice free now, therefore it
           ! can be either ocean or land (using actual maske value).

           if (st%maske(j,i) == 2) then   ! ocean point (hidden melt)
              st%mask_ablation_type(j,i) = -2
           else                            ! land (hidden melt)
              st%mask_ablation_type(j,i) = -1
           end if

        end if

        if (st%calving_apl(j,i).ne.st%calving_apl(j,i)) then
          print *,'calving_apl is NaN',i,j
          print *,'H_neu, H_neu_flow',st%H_neu(j,i),st%H_neu_flow(j,i)
          print *,'accum',st%accum(j,i)
          print *,'runoff',st%runoff(j,i)
          print *,'Q_b',st%Q_b_tot(j,i)
          print *,'as_perp',st%as_perp(j,i)
          print *,'calving',st%calving(j,i)
          print *,'accum_apl',st%accum_apl(j,i)
          print *,'runoff_apl',st%runoff_apl(j,i)
          print *,'Q_b_apl',st%Q_b_apl(j,i)
          print *,'as_perp_apl',st%as_perp_apl(j,i)
          print *,'calving_apl',st%calving_apl(j,i)
          print *,'mb_source_apl',st%mb_source_apl(j,i)
        endif
        if (abs(st%calving_apl(j,i)).gt.1.e-5_wp) then
          print *,'calving_apl is large: ',i,j,st%calving_apl(j,i)
          print *,'H_neu, H_neu_flow',st%H_neu(j,i),st%H_neu_flow(j,i)
          print *,'mb_source',st%mb_source(j,i)
          print *,'accum',st%accum(j,i)
          print *,'runoff',st%runoff(j,i)
          print *,'Q_b_tot',st%Q_b_tot(j,i)
          print *,'Q_bm',st%Q_bm(j,i)
          print *,'Q_bm_float',st%Q_bm_float(j,i)
          print *,'Q_tld',st%Q_tld(j,i)
          print *,'as_perp',st%as_perp(j,i)
          print *,'calving',st%calving(j,i)
          print *,'accum_apl',st%accum_apl(j,i)
          print *,'runoff_apl',st%runoff_apl(j,i)
          print *,'Q_b_apl',st%Q_b_apl(j,i)
          print *,'as_perp_apl',st%as_perp_apl(j,i)
          print *,'calving_apl',st%calving_apl(j,i)
          print *,'mb_source_apl',st%mb_source_apl(j,i)
        endif

     end if

  end do
  end do

end subroutine account_mb_source

!-------------------------------------------------------------------------------

end module calc_thk_m
!
