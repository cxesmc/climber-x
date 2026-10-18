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
  use climber_grid, only : lon
  use ocn_params, only : drag, drag_bcl, rtv, rtv3, fcor, fcorv, fcormin, fcormin_ref, dt
  use ocn_params, only : drag_par
  use ocn_params, only : bering_tf_fix, drake_tf_fix, indo_tf_fix, davis_tf_fix, rho0
  use ocn_params, only : i_cor_form, nbw, qcor
  use ocn_grid, only : maxi, maxj, maxk, mpxi, mpxj 
  use ocn_grid, only : maxisles, n_isles, psiles, map_isles
  use ocn_grid, only : lm_eur, lm_ame, lm_ant, lm_aus, lm_grl, lm_name, i_lm, j_lm
  use ocn_grid, only : dx, dy, dz, zro, s, sv, cv, rh, k1, h

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
  ! contiguous band-storage copies of the LU factors used by ubarsolv (built once per mask
  ! update from ratm/gap), so the per-step substitutions are unit-stride
  real(wp), allocatable :: Lband(:,:)   ! lower factors,  Lband(off,i) = ratm(i+off,off)
  real(wp), allocatable :: Uband(:,:)   ! upper factors,  Uband(off,i) = gap(i,nbw+1+off)
  real(wp), allocatable :: Udiag(:)     ! diagonal,        Udiag(i)    = gap(i,nbw+1)
  real(wp), allocatable :: psisl(:,:,:)
  real(wp), allocatable :: ubisl(:,:,:,:)
  real(wp), allocatable :: erisl(:,:)
  integer, allocatable :: erisl_ipiv(:)   ! row-pivot vector from matinv, replayed by matmult each step
  real(wp), allocatable :: psibc(:)
  ! prescribed barotropic throughflows between landmasses (bering/drake/indo/davis_tf_fix):
  ! psibc(isl) - psibc(isl_fix_ref(isl)) = psi_fix(isl), with psibc(0) = 0 for the mainland
  integer, parameter :: n_strait = 4
  character(len=10), parameter :: strait_name(n_strait) = ["Bering    ","Drake     ","Indonesian","Davis     "]
  ! landmass on the positive side (target) and on the negative side (reference) of each strait,
  ! so that the throughflow is (psi(target)-psi(reference))/rho0 in the sign convention of ocn_ts.nc:
  ! Bering/Davis positive north (target to the east), Drake/Indonesian positive east (target to the south)
  integer, parameter :: strait_lm_tgt(n_strait) = [lm_ame, lm_ant, lm_aus, lm_grl]
  integer, parameter :: strait_lm_ref(n_strait) = [lm_eur, lm_ame, lm_eur, lm_ame]
  integer :: strait_status(n_strait) = -1   ! 0 = off, 1 = active, 2 = reference cell not land, 3 = same landmass on both sides
  integer :: strait_isl(n_strait) = 0       ! solver island index whose streamfunction is prescribed for each active strait
  integer, allocatable :: isl_fix_ref(:)    ! solver index of the reference landmass (0 = mainland), -1 = island not prescribed
  real(wp), allocatable :: psi_fix(:)       ! prescribed streamfunction difference to the reference landmass [kg/s]
  real(wp), allocatable :: tmpdrg(:,:)

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

    ! wind stress term for barotropic streamfunction
    !$ time1 = omp_get_wtime()
    call wind(tau, &    ! in 
              ubar_wind)       ! out
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'momentum: wind',time2-time1

    ! jbar term for barotropic streamfunction
    !$ time1 = omp_get_wtime()
    call jbar(rho, &    ! in
              bp, ubar_jbar)     ! out
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
    call ubarsolv(Lband,Uband,Udiag, &
                  gb, &
                  ub,psi)
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'momentum: ubarsolv',time2-time1

    ! find island path integral due to wind and jbar terms 
    !$ time1 = omp_get_wtime()
    do isl=1,n_isles
       call island(ub,tau,bp, &
                   isl,1, &
                   erisl(isl,n_isles+1), psi)
    enddo
    ! prescribed throughflows: the row of a prescribed island reads psibc(isl) - psibc(ref) = -rhs (momentum_update_grid)
    do isl=1,n_isles
      if (isl_fix_ref(isl).ge.0) erisl(isl,n_isles+1) = - psi_fix(isl)
    enddo
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'momentum: islands',time2-time1

    ! solve system of simultaneous equations
    !$ time1 = omp_get_wtime()
    if (n_isles > 1) then
       call matmult(n_isles,erisl,erisl_ipiv,erisl(:,n_isles+1))
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
    call velc(rho,dtau_dz2,dtav_dz2,ub, &
              u)
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
    integer :: i1, i1p, j1, ii
    integer :: n_band, nm_band, off
    real(wp) :: f_wet, h_wet


    ! barotropic drag at psi points (grid cell corners!): background adrag plus two independent enhancements
    ! from the wet cells around the psi point,
    !  - partial-cell drag drag_frac*(min(1/f_ocn,drag_frac_max)-1) of the least-wet cell: the flow through a cell
    !    that is only partly ocean is squeezed into a fraction f_ocn of its width
    !  - linearised bottom friction cd_bot*u_bot/H of the shallowest cell: damps shallow cells, negligible in the deep ocean
    do j=0,maxj
      do i=0,maxi
        ii = i
        if (ii.eq.0) ii=maxi
        f_wet = 1._wp
        h_wet = -1._wp
        do j1=max(1,j),min(maxj,j+1)
          do i1=ii,ii+1
            i1p = 1 + mod(i1-1,maxi)
            if (k1(i1p,j1).le.maxk) then
              f_wet = min(f_wet,f_ocn(i1p,j1))
              if (h_wet.lt.0._wp .or. h(3,i1p,j1).lt.h_wet) h_wet = h(3,i1p,j1)
            endif
          enddo
        enddo
        tmpdrg(i,j) = drag_par%adrag + drag_par%drag_frac*(min(1._wp/max(f_wet,1.e-3_wp),drag_par%drag_frac_max)-1._wp)
        if (h_wet.gt.0._wp) tmpdrg(i,j) = tmpdrg(i,j) + drag_par%cd_bot*drag_par%u_bot/h_wet
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

    ! baroclinic velocities: background drag only
    drag_bcl(:,:,:) = drag_par%adrag

    ! on rows where the Coriolis parameter is set to zero (|f| < fcormin, see momentum_init) add drag instead,
    ! r_eff = r + fcormin^2/r, so that r_eff/(0+r_eff^2) = r/(fcormin^2+r^2)
    do j=1,maxj
      if (abs(2._wp*omega*s(j)).lt.fcormin) then     ! u-grid rows (tracer latitudes)
        do i=1,maxi
          drag(1,i,j)     = drag(1,i,j)     + fcormin**2/drag(1,i,j)
          drag_bcl(1,i,j) = drag_bcl(1,i,j) + fcormin**2/drag_bcl(1,i,j)
        enddo
      endif
      if (abs(2._wp*omega*sv(j)).lt.fcormin) then    ! v-grid rows (e.g. the equator line)
        do i=1,maxi+1
          drag(2,i,j)     = drag(2,i,j)     + fcormin**2/drag(2,i,j)
          drag_bcl(2,i,j) = drag_bcl(2,i,j) + fcormin**2/drag_bcl(2,i,j)
        enddo
      endif
    enddo


    ! arrays for efficiency
    do j=1,maxj
      do i=1,maxi
        rtv(i,j) = 1._wp/(fcor(j)**2 + drag_bcl(1,i,j)*drag_bcl(1,i,j))
        rtv3(i,j) = 1._wp/(fcorv(j)**2 + drag_bcl(2,i,j)*drag_bcl(2,i,j))
      enddo
    enddo


    ! for barotropic velocity

    if (i_cor_form.eq.1) then
      ! f/H at tracer points for the J3 Coriolis term (zero on land, so dry cells drop out of the stencil)
      qcor(:,:) = 0._wp
      do j=1,maxj
        do i=0,maxi+1
          if (k1(i,j).le.maxk) qcor(i,j) = fcor(j)*rh(3,i,j)
        enddo
      enddo
    endif

    call invert(gap,ratm)

    ! repack the banded LU factors into contiguous band storage so the per-step ubarsolv
    ! substitutions are unit-stride (cache- and SIMD-friendly). done once per mask update;
    ! the arithmetic is unchanged, only the memory layout differs
    n_band  = nbw
    nm_band = mpxi*mpxj
    do i=1,nm_band
      Udiag(i) = gap(i,n_band+1)
      do off=1,min(n_band,nm_band-i)
        Lband(off,i) = ratm(i+off,off)        ! lower factor (forward substitution)
        Uband(off,i) = gap(i,n_band+1+off)     ! upper factor (back substitution)
      enddo
    enddo

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
       call ubarsolv(Lband,Uband,Udiag, &
                     gb,ubisl(:,:,:,isol),psisl(:,:,isol))

       ! find island path integral due to unit source on boundary
       do isl=1,n_isles
          call island(ubisl(:,:,:,isol),tau,bp, &  ! tau and bp not used in this call (indj==0)!
                      isl,0, &
                      erisl(isl,isol), psisl(:,:,isol))  ! erisl in m/s2
       enddo
    enddo

    !print*,'island path integrals due to unit sources'
    !do isl=1,n_isles
    !  do isol=1,n_isles
    !    print *,isl,isol,erisl(isl,isol)
    !  enddo
    !enddo

    ! prescribed throughflows: replace the path-integral equations of the prescribed islands
    call strait_constraints

    ! partially invert inland integral error matrix for psi bc calc
    call matinv(n_isles,erisl,erisl_ipiv)  ! LU factors of island matrix (+ pivots)

    !print*,'island path integrals due to unit sources after matinv'
    !do isl=1,n_isles
    !  do isol=1,n_isles
    !    print *,isl,isol,erisl(isl,isol)
    !  enddo
    !enddo


    return

  end subroutine momentum_update_grid


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s t r a i t _ c o n s t r a i n t s
  !   Purpose    :  prescribe the barotropic throughflow through straits
  !                 between two landmasses
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine strait_constraints

    ! The net barotropic transport between two landmasses equals the difference of their streamfunction values
    ! divided by rho0, whatever the path of the section. A throughflow is prescribed by replacing the path-integral
    ! equation of one of the two landmasses (target) by psibc(target) - psibc(ref) = rho0*tf_fix, with psibc = 0 on
    ! the mainland. The island whose equation is replaced is the designated target (America for Bering, Antarctica
    ! for Drake, Australia for the Indonesian throughflow, Greenland for Davis), or the other side if the target is
    ! the mainland or is already prescribed. The landmasses are identified by their reference cells (ocn_grid), so the
    ! constraints follow the island configuration (i_isl, l_isl_*): a strait is ignored when both sides belong to the
    ! same landmass. The constraints must not form a closed loop of landmasses, which would overdetermine the system.

    implicit none

    integer :: n, it, ir, rt, rr, isl, lm_fix
    integer :: status(n_strait), isl_strait(n_strait)
    integer :: lm_strait(n_strait)
    integer :: uf(0:maxisles)     ! union-find parent of the landmass graph spanned by the constraints (0 = mainland)
    real(wp) :: tf_fix(n_strait), sgn
    logical :: l_print


    tf_fix = [bering_tf_fix, drake_tf_fix, indo_tf_fix, davis_tf_fix]

    isl_fix_ref(:) = -1
    psi_fix(:) = 0._wp
    status(:) = 0
    isl_strait(:) = 0
    lm_strait(:) = 0
    do n=0,maxisles
      uf(n) = n
    enddo

    do n=1,n_strait
      if (tf_fix(n).lt.-900._wp) cycle

      ! solver indices of the landmasses on both sides (landmass label - 1, mainland = 0, ocean < 0)
      it = map_isles(i_lm(strait_lm_tgt(n)),j_lm(strait_lm_tgt(n))) - 1
      ir = map_isles(i_lm(strait_lm_ref(n)),j_lm(strait_lm_ref(n))) - 1
      if (it.lt.0 .or. ir.lt.0) then
        status(n) = 2
        cycle
      endif
      if (it.eq.ir) then
        status(n) = 3
        cycle
      endif

      ! a second constraint between landmasses that are already connected by constraints closes a loop
      rt = uf_root(it)
      rr = uf_root(ir)
      if (rt.eq.rr) then
        print *,'ERROR: prescribed ',trim(strait_name(n)),' throughflow closes a loop of prescribed throughflows'
        print *,'between ',trim(lm_name(strait_lm_tgt(n))),' and ',trim(lm_name(strait_lm_ref(n))),': disable one of them'
        stop 'ERROR in strait_constraints: overdetermined prescribed throughflows'
      endif

      ! island whose equation is replaced: designated target, else the reference side with reversed sign
      sgn = 1._wp
      lm_fix = strait_lm_tgt(n)
      if (.not.is_free(it)) then
        isl = it
        it = ir
        ir = isl
        sgn = -1._wp
        lm_fix = strait_lm_ref(n)
      endif
      if (.not.is_free(it)) then
        print *,'ERROR: no free island left to prescribe the ',trim(strait_name(n)),' throughflow:'
        print *,'both ',trim(lm_name(strait_lm_tgt(n))),' and ',trim(lm_name(strait_lm_ref(n))), &
          ' are the mainland or already prescribed by another strait'
        stop 'ERROR in strait_constraints: prescribed throughflows cannot be assigned to islands'
      endif

      isl_fix_ref(it) = ir
      psi_fix(it) = sgn*rho0*tf_fix(n)*1.e6_wp  ! kg/s
      uf(rt) = rr
      status(n) = 1
      isl_strait(n) = it
      lm_strait(n) = lm_fix
    enddo

    ! replace the island rows
    do isl=1,n_isles
      if (isl_fix_ref(isl).ge.0) then
        erisl(isl,1:n_isles) = 0._wp
        erisl(isl,isl) = 1._wp
        if (isl_fix_ref(isl).gt.0) erisl(isl,isl_fix_ref(isl)) = -1._wp
      endif
    enddo

    ! report on the first call and whenever the configuration changes (e.g. a strait closing with the land-sea mask)
    l_print = any(status.ne.strait_status) .or. any(isl_strait.ne.strait_isl)
    if (l_print) then
      do n=1,n_strait
        select case (status(n))
        case (1)
          print '(a,a,a,f8.3,a,a,a,i0,a)','momentum: ',trim(strait_name(n)),' throughflow prescribed to ',tf_fix(n), &
            ' Sv (streamfunction of ',trim(lm_name(lm_strait(n))),', island ',isl_strait(n),', fixed)'
        case (2)
          print *,'WARNING: ',trim(strait_name(n)),' throughflow prescribed, but the reference cell of ', &
            trim(lm_name(strait_lm_tgt(n))),' or ',trim(lm_name(strait_lm_ref(n))),' is not land: ignored'
        case (3)
          print *,'WARNING: ',trim(strait_name(n)),' throughflow prescribed, but ',trim(lm_name(strait_lm_tgt(n))), &
            ' and ',trim(lm_name(strait_lm_ref(n))),' belong to the same landmass (not an island or strait closed): ignored'
        end select
      enddo
    endif
    strait_status = status
    strait_isl = isl_strait

    return

  contains

    ! island (not the mainland) whose equation is not yet replaced by a constraint
    logical function is_free(k)
      integer, intent(in) :: k
      is_free = .false.
      if (k.gt.0) is_free = isl_fix_ref(k).lt.0
    end function is_free

    integer function uf_root(k)
      integer, intent(in) :: k
      uf_root = k
      do while (uf(uf_root).ne.uf_root)
        uf_root = uf(uf_root)
      enddo
    end function uf_root

  end subroutine strait_constraints


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
    allocate(qcor(0:maxi+1,0:maxj+1))
    qcor(:,:) = 0._wp
    ! half band width of the streamfunction matrix: 5-point stencil for both i_cor_form, 9-point with the
    ! checkerboard friction, whose diagonal neighbours across the periodic boundary are 2*maxi-1 rows away
    if (drag_par%drag_cb.gt.0._wp) then
      nbw = 2*mpxi
    else
      nbw = mpxi+1
    endif
    allocate(ratm(mpxi*mpxj,nbw))
    allocate(gb(mpxi*mpxj))
    allocate(gap(mpxi*mpxj,2*nbw+1))
    allocate(Lband(nbw,mpxi*mpxj))
    allocate(Uband(nbw,mpxi*mpxj))
    allocate(Udiag(mpxi*mpxj))
    allocate(bp(maxi+1,maxj,maxk))
    allocate(psisl(0:maxi,0:maxj,maxisles))
    allocate(ubisl(2,0:maxi+1,0:maxj,maxisles))
    allocate(erisl(maxisles,maxisles+1))
    allocate(erisl_ipiv(maxisles))
    allocate(psibc(maxisles))
    allocate(isl_fix_ref(maxisles))
    allocate(psi_fix(maxisles))
    isl_fix_ref(:) = -1
    psi_fix(:) = 0._wp
    allocate(drag(2,maxi+1,maxj))
    allocate(drag_bcl(2,maxi+1,maxj))
    allocate(tmpdrg(0:maxi,0:maxj))
    allocate(fcor(0:maxj+1))
    allocate(fcorv(0:maxj+1))
    allocate(rtv(maxi,maxj))
    allocate(rtv3(maxi,maxj))


    drag_par%adrag = 1._wp/(drag_par%adrag*86400._wp) ! 1/s
    if (drag_par%drag_frac.gt.0._wp) drag_par%drag_frac = 1._wp/(drag_par%drag_frac*86400._wp) ! 1/s
    if (drag_par%drag_cb.gt.0._wp) drag_par%drag_cb = 1._wp/(drag_par%drag_cb*86400._wp) ! 1/s

    ! Coriolis parameter
    ! Rows with |f| < fcormin (only the equator v-line on the 5 deg grid) get f = 0 and instead a larger drag,
    ! r_eff = r + fcormin^2/r (see momentum_update_grid), which keeps the frictional-geostrophic mobility
    ! r/(f^2+r^2) at the value a Coriolis floor of fcormin would give - the same numerical damping of the
    ! cross-equatorial pressure gradient, which is needed for stability (fcormin must stay >~ the drag rate) -
    ! but is symmetric about the equator instead of deflecting the equatorial flow like a NH row.
    do j=1,maxj
      if (abs(2._wp*omega*s(j)).lt.fcormin) then
        fcor(j) = 0._wp
      else
        fcor(j) = 2._wp*omega*s(j)
      endif
    enddo
    do j=0,maxj
      if (abs(2._wp*omega*sv(j)).lt.fcormin) then
        fcorv(j) = 0._wp
      else
        fcorv(j) = 2._wp*omega*sv(j)
      endif
    enddo

    return

  end subroutine momentum_init

end module momentum_mod
 
