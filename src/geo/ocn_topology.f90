!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module : o c n _ t o p o l o g y _ m o d
!
!  Purpose : topology-preserving coarse ocean mask
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
!
! CLIMBER-X is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CLIMBER-X.  If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Thresholding the ocean fraction with f_crit decides connectivity by accident: a strait cell below the threshold
! disappears (Hormuz, Bab-el-Mandeb, Torres), and two coarse ocean cells on either side of a high resolution land or
! ice barrier get an open face (Panama, Suez, Nares at the LGM). Here the coarse mask is made locally faithful to the
! connectivity of the high resolution ocean:
!  - high resolution graph H: nodes are the ocean components within each coarse cell (a cell holding Pacific and
!    Caribbean water has two nodes), edges are high resolution ocean pixels touching across a coarse cell face
!  - coarse graph C: nodes are the coarse ocean cells connected to the ocean origin, edges are all faces between
!    ocean cells, as opened by the ocean model
!  shortcut : two nodes within k_ocn_topo cells in C, but more than n_excess_ocn_topo hops further apart in H
!  missing  : two nodes within k_ocn_topo hops in H, but more than n_excess_ocn_topo hops further apart in C
! Only nodes larger than f_sub_ocn_topo of a cell are checked, and a second ocean body inside a cell is checked only
! if it connects to another ocean cell (a sliver of another ocean in a dead-end corner is harmless).
! Re-open the high resolution path of a missing connection through the least land, close the
! cell with the smallest ocean fraction on the coarse path of a shortcut, iterate, and finally undo the changes that
! turn out not to be needed. 
!
module ocn_topology_mod

  use precision, only : wp
  use geo_params, only : z_sill_ocn_topo, f_sub_ocn_topo, k_ocn_topo, n_excess_ocn_topo
  use geo_grid, only : ni, nj, n_topo_sur, i0_topo, i1_topo, j0_topo, j1_topo
  use fill_ocean_mod, only : fill_ocean_lowres

  implicit none

  integer, parameter :: d_inf = huge(1)
  integer, parameter :: kind_shortcut = 1
  integer, parameter :: kind_missing = 2
  integer, parameter :: maxit = 100

  ! high resolution graph, the nodes of a coarse cell are stored contiguously
  integer :: nnode
  integer, allocatable :: node_i(:), node_j(:), node_npix(:)
  integer, allocatable :: cell_n0(:,:), cell_n1(:,:)    ! first and last node of each coarse cell
  integer, allocatable :: cell_main(:,:)                ! largest node of each coarse cell
  integer, allocatable :: adj0(:), adj(:)               ! neighbours of node p are adj(adj0(p):adj0(p+1)-1)
  logical, allocatable :: sig(:)                        ! node large enough to be checked

  ! breadth-first search work arrays (-1 = not visited), reset after each search
  integer, allocatable :: dist_h(:), queue_h(:)
  integer :: nqueue_h
  integer, allocatable :: dist_c(:,:), queue_c(:,:)
  integer :: nqueue_c

  ! flags found by the last check
  integer :: nflag
  integer, allocatable :: flag_kind(:), flag_p(:), flag_q(:), flag_dc(:), flag_dh(:)

  ! flags that could not be resolved
  integer :: nignored
  integer :: ignored_p(maxit), ignored_q(maxit)

  real(wp), allocatable :: lon_c(:), lat_c(:)
  real(wp) :: lon_origin, lat_origin

  ! mask from the previous call, to report changes only when the result changes
  integer, allocatable, save :: mask_prev(:,:)

  private
  public :: ocn_topology

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  o c n _ t o p o l o g y
  ! Purpose    :  repair the coarse ocean mask so that it preserves the connectivity of the high resolution ocean
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine ocn_topology(hires_mask, hires_z_bed, lon, lat, lon_ocn_origin, lat_ocn_origin, f_ocn, mask_fix, &
      mask_ocn)

    implicit none

    integer, intent(in) :: hires_mask(:,:)       !! high resolution mask
    real(wp), intent(in) :: hires_z_bed(:,:)     !! high resolution bedrock elevation relative to sea level [m]
    real(wp), intent(in) :: lon(:), lat(:)       !! coarse grid coordinates
    real(wp), intent(in) :: lon_ocn_origin, lat_ocn_origin
    real(wp), intent(in) :: f_ocn(:,:)           !! ocean fraction before applying any threshold
    logical, intent(in) :: mask_fix(:,:)         !! cells that must remain land (explicit strait closures)
    integer, intent(inout) :: mask_ocn(:,:)      !! coarse ocean mask (1 = ocean, 0 = land), after thresholding

    integer :: it, n, i, j, nmiss0, nmiss_all, nshort_ini, nmiss_ini, nlive0, worst, wp_, wq_, nchange, npath, ncand
    real(wp) :: score, score_max
    logical :: done
    logical, dimension(ni,nj) :: keep, C, locked
    integer, dimension(4*maxit) :: change_i, change_j, change_a
    integer, allocatable :: path(:), cand_i(:), cand_j(:)


    allocate(lon_c(ni), lat_c(nj))
    lon_c = lon
    lat_c = lat
    lon_origin = lon_ocn_origin
    lat_origin = lat_ocn_origin

    call build_graph(hires_mask, hires_z_bed)

    allocate(path(nnode), cand_i(ni*nj), cand_j(ni*nj))
    nchange = 0
    nignored = 0

    keep = mask_ocn.eq.1 .and. .not.mask_fix
    locked = .false.    ! cells re-opened for a missing connection, never closed again

    do it=1,maxit

      call fill(keep, C)
      call check(keep, C)
      if (it.eq.1) then
        nshort_ini = nlive(kind_shortcut)
        nmiss_ini = nlive(kind_missing)
      endif
      if (nlive(0).eq.0) exit

      nmiss0 = nlive(kind_missing)
      nmiss_all = nall(kind_missing)

      if (nmiss0.gt.0) then

        ! worst missing connection first: re-open the high resolution path that crosses the least land
        worst = 0
        score_max = -1._wp
        do n=1,nflag
          if (flag_kind(n).ne.kind_missing .or. is_ignored(flag_p(n),flag_q(n))) cycle
          if (flag_dc(n).eq.d_inf) then
            score = 1.e6_wp
          else
            score = real(flag_dc(n)-flag_dh(n),wp)
          endif
          if (score.gt.score_max) then
            score_max = score
            worst = n
          endif
        enddo
        wp_ = flag_p(worst)
        wq_ = flag_q(worst)
        call path_min_land(keep, C, mask_fix, f_ocn, wp_, wq_, path, npath)
        done = .false.
        do n=1,npath
          i = node_i(path(n))
          j = node_j(path(n))
          if (.not.keep(i,j)) then
            keep(i,j) = .true.
            locked(i,j) = .true.
            call add_change(i, j, 1)
            done = .true.
          endif
        enddo
        if (.not.done) call ignore(wp_, wq_)

      else

        ! worst shortcut first: close the cell with the smallest ocean fraction on the coarse path
        worst = 0
        score_max = -1._wp
        do n=1,nflag
          if (flag_kind(n).ne.kind_shortcut .or. is_ignored(flag_p(n),flag_q(n))) cycle
          if (flag_dh(n).eq.d_inf) then
            score = 1.e6_wp
          else
            score = real(flag_dh(n)-flag_dc(n),wp)
          endif
          if (score.gt.score_max) then
            score_max = score
            worst = n
          endif
        enddo
        wp_ = flag_p(worst)
        wq_ = flag_q(worst)
        call path_c(C, node_i(wp_), node_j(wp_), node_i(wq_), node_j(wq_), cand_i, cand_j, ncand)
        call sort_cells(f_ocn, cand_i, cand_j, ncand)
        done = .false.
        do n=1,ncand
          i = cand_i(n)
          j = cand_j(n)
          if (locked(i,j)) cycle
          keep(i,j) = .false.
          call fill(keep, C)
          call check(keep, C)
          ! all missing connections count here, also the ignored ones, so that a closure can not re-create them
          if (nall(kind_missing).le.nmiss_all) then
            call add_change(i, j, -1)
            done = .true.
            exit
          endif
          keep(i,j) = .true.
        enddo
        if (.not.done) call ignore(wp_, wq_)

      endif

      if (nchange.ge.size(change_i)-ni) exit
      if (nignored.ge.maxit) exit

    enddo

    if (it.gt.maxit) print *,'WARNING: ocn_topology: repair did not converge in',maxit,' iterations'

    ! undo changes that are not needed, closures of the largest ocean fractions first
    ! (all flags count here, also the ignored ones, so that pruning can not bring back an unresolved flag)
    call fill(keep, C)
    call check(keep, C)
    nlive0 = nflag
    call sort_changes(f_ocn, change_i, change_j, change_a, nchange)
    do n=1,nchange
      i = change_i(n)
      j = change_j(n)
      keep(i,j) = .not.keep(i,j)
      call fill(keep, C)
      call check(keep, C)
      if (nflag.le.nlive0) then
        change_a(n) = 0   ! not needed
      else
        keep(i,j) = .not.keep(i,j)
      endif
    enddo

    where (keep)
      mask_ocn = 1
    elsewhere
      mask_ocn = 0
    endwhere

    ! flags of the final mask, to report the ones that remained unresolved
    call fill(keep, C)
    call check(keep, C)

    ! report, only when the mask differs from the previous call
    if (.not.allocated(mask_prev)) then
      allocate(mask_prev(ni,nj))
      mask_prev = -1
    endif
    if (any(mask_prev.ne.mask_ocn)) then
      print '(a,i5,a,i5,a)',' ocn_topology: thresholded mask has',nshort_ini,' shortcuts and',nmiss_ini,' missing connections'
      do n=1,nchange
        if (change_a(n).eq.0) cycle
        i = change_i(n)
        j = change_j(n)
        print '(a,a,2i4,a,f7.1,a,f6.1,a,f5.2)',' ocn_topology: ',merge('open  cell','close cell',change_a(n).eq.1), &
          i,j,' lon',lon_c(i),' lat',lat_c(j),' f_ocn',f_ocn(i,j)
      enddo
      do n=1,nignored
        if (.not.is_flagged(ignored_p(n),ignored_q(n))) cycle
        print '(a,2(a,f7.1,a,f6.1,a))',' WARNING: ocn_topology: unresolved', &
          ' (',lon_c(node_i(ignored_p(n))),',',lat_c(node_j(ignored_p(n))),')', &
          ' (',lon_c(node_i(ignored_q(n))),',',lat_c(node_j(ignored_q(n))),')'
      enddo
      mask_prev = mask_ocn
    endif

    deallocate(path, cand_i, cand_j)
    call free_graph

  contains

    subroutine add_change(i, j, a)
      integer, intent(in) :: i, j, a
      nchange = nchange+1
      change_i(nchange) = i
      change_j(nchange) = j
      change_a(nchange) = a
    end subroutine add_change

  end subroutine ocn_topology


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  b u i l d _ g r a p h
  ! Purpose    :  ocean components within each coarse cell and their connections across coarse cell faces
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine build_graph(hires_mask, hires_z_bed)

    implicit none

    integer, intent(in) :: hires_mask(:,:)
    real(wp), intent(in) :: hires_z_bed(:,:)

    integer :: i, j, ii, jj, a, b, n, m, p, ip1, ns, nmax, ne, nface, npix
    integer, allocatable :: lab(:,:), si(:), sj(:), ea(:), eb(:), fa(:), fb(:), deg(:), pos(:)
    logical, allocatable :: ocn(:,:), pass(:,:)


    allocate(ocn(size(hires_mask,1),size(hires_mask,2)))
    allocate(pass(size(hires_mask,1),size(hires_mask,2)))
    allocate(lab(size(hires_mask,1),size(hires_mask,2)))
    ocn = hires_mask.eq.2
    ! a crossing counts only where both pixels are deeper than the sill depth
    pass = ocn .and. hires_z_bed.lt.-z_sill_ocn_topo
    lab = 0

    ! label the ocean components within each coarse cell
    nmax = 0
    do j=1,nj
      do i=1,ni
        nmax = max(nmax, (i1_topo(i)-i0_topo(i)+1)*(j1_topo(j)-j0_topo(j)+1))
      enddo
    enddo
    allocate(si(nmax), sj(nmax))
    allocate(node_i(4*ni*nj), node_j(4*ni*nj), node_npix(4*ni*nj))
    allocate(cell_n0(ni,nj), cell_n1(ni,nj), cell_main(ni,nj))
    nnode = 0
    do j=1,nj
      do i=1,ni
        cell_n0(i,j) = nnode+1
        cell_main(i,j) = 0
        do jj=j0_topo(j),j1_topo(j)
          do ii=i0_topo(i),i1_topo(i)
            if (.not.ocn(ii,jj) .or. lab(ii,jj).ne.0) cycle
            ! new node, flood fill within the cell
            if (nnode.eq.size(node_i)) then
              call grow(node_i) ; call grow(node_j) ; call grow(node_npix)
            endif
            nnode = nnode+1
            node_i(nnode) = i
            node_j(nnode) = j
            npix = 0
            ns = 1
            si(1) = ii
            sj(1) = jj
            lab(ii,jj) = nnode
            do while (ns.gt.0)
              a = si(ns)
              b = sj(ns)
              ns = ns-1
              npix = npix+1
              call push(a-1, b) ; call push(a+1, b) ; call push(a, b-1) ; call push(a, b+1)
            enddo
            node_npix(nnode) = npix
            if (cell_main(i,j).eq.0) then
              cell_main(i,j) = nnode
            else if (npix.gt.node_npix(cell_main(i,j))) then
              cell_main(i,j) = nnode
            endif
          enddo
        enddo
        cell_n1(i,j) = nnode
      enddo
    enddo

    ! connections across coarse cell faces, unique per face
    allocate(ea(8*ni*nj), eb(8*ni*nj))
    allocate(fa(nmax), fb(nmax))
    ne = 0
    do j=1,nj
      do i=1,ni
        ! eastern face, periodic
        ip1 = i+1
        if (ip1.gt.ni) ip1 = 1
        nface = 0
        do jj=j0_topo(j),j1_topo(j)
          call cross(i1_topo(i), jj, i0_topo(ip1), jj)
        enddo
        call add_face
        ! northern face
        if (j.lt.nj) then
          nface = 0
          do ii=i0_topo(i),i1_topo(i)
            call cross(ii, j1_topo(j), ii, j0_topo(j+1))
          enddo
          call add_face
        endif
      enddo
    enddo

    ! adjacency in compressed row storage
    allocate(deg(nnode), pos(nnode), adj0(nnode+1), adj(2*ne))
    deg = 0
    do n=1,ne
      deg(ea(n)) = deg(ea(n))+1
      deg(eb(n)) = deg(eb(n))+1
    enddo
    adj0(1) = 1
    do p=1,nnode
      adj0(p+1) = adj0(p)+deg(p)
    enddo
    pos = adj0(1:nnode)
    do n=1,ne
      adj(pos(ea(n))) = eb(n) ; pos(ea(n)) = pos(ea(n))+1
      adj(pos(eb(n))) = ea(n) ; pos(eb(n)) = pos(eb(n))+1
    enddo

    allocate(sig(nnode))
    sig = node_npix(1:nnode).ge.f_sub_ocn_topo*real(n_topo_sur,wp)

    allocate(dist_h(nnode), queue_h(nnode))
    dist_h = -1
    nqueue_h = 0
    allocate(dist_c(ni,nj), queue_c(2,ni*nj))
    dist_c = -1
    nqueue_c = 0
    allocate(flag_kind(256), flag_p(256), flag_q(256), flag_dc(256), flag_dh(256))
    nflag = 0

    deallocate(ocn, pass, lab, si, sj, ea, eb, fa, fb, deg, pos)

  contains

    subroutine push(a, b)
      integer, intent(in) :: a, b
      if (a.lt.i0_topo(i) .or. a.gt.i1_topo(i) .or. b.lt.j0_topo(j) .or. b.gt.j1_topo(j)) return
      if (.not.ocn(a,b) .or. lab(a,b).ne.0) return
      lab(a,b) = nnode
      ns = ns+1
      si(ns) = a
      sj(ns) = b
    end subroutine push

    subroutine cross(ia, ja, ib, jb)
      integer, intent(in) :: ia, ja, ib, jb
      integer :: k, la, lb
      if (.not.(pass(ia,ja) .and. pass(ib,jb))) return
      la = lab(ia,ja)
      lb = lab(ib,jb)
      do k=1,nface
        if (fa(k).eq.la .and. fb(k).eq.lb) return
      enddo
      nface = nface+1
      fa(nface) = la
      fb(nface) = lb
    end subroutine cross

    subroutine add_face
      integer :: k
      do k=1,nface
        if (ne.eq.size(ea)) then
          call grow(ea) ; call grow(eb)
        endif
        ne = ne+1
        ea(ne) = fa(k)
        eb(ne) = fb(k)
      enddo
    end subroutine add_face

  end subroutine build_graph


  subroutine free_graph

    deallocate(node_i, node_j, node_npix, cell_n0, cell_n1, cell_main, adj0, adj, sig)
    deallocate(dist_h, queue_h, dist_c, queue_c)
    deallocate(flag_kind, flag_p, flag_q, flag_dc, flag_dh)
    deallocate(lon_c, lat_c)

  end subroutine free_graph


  subroutine grow(x)

    integer, allocatable, intent(inout) :: x(:)
    integer, allocatable :: tmp(:)

    allocate(tmp(2*size(x)))
    tmp(1:size(x)) = x
    call move_alloc(tmp, x)

  end subroutine grow


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  c h e c k
  ! Purpose    :  find shortcuts and missing connections of the coarse mask
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine check(C0, C)

    implicit none

    logical, intent(in) :: C0(:,:)   !! ocean cells after thresholding, including cells cut off from the ocean
    logical, intent(in) :: C(:,:)    !! ocean cells connected to the ocean origin

    integer :: p, q, n, ic, jc, dc, dh, dlim
    logical, allocatable :: act(:), act0(:)


    allocate(act(nnode), act0(nnode))
    call active(C, act)
    call active(C0, act0)
    dlim = k_ocn_topo+n_excess_ocn_topo
    nflag = 0

    do p=1,nnode
      if (.not.(act(p) .or. act0(p))) cycle

      call bfs_h(p, dlim)

      ! shortcut: q within k cells in C, but further than dc+excess in H
      if (act(p)) then
        call bfs_c(C, node_i(p), node_j(p), k_ocn_topo)
        do n=1,nqueue_c
          ic = queue_c(1,n)
          jc = queue_c(2,n)
          dc = dist_c(ic,jc)
          do q=cell_n0(ic,jc),cell_n1(ic,jc)
            if (q.le.p .or. .not.act(q)) cycle
            dh = dist_h(q)
            if (dh.lt.0) dh = d_inf
            if (dh.eq.d_inf .or. dh.gt.dc+n_excess_ocn_topo) call add_flag(kind_shortcut, p, q, dc, dh)
          enddo
        enddo
        call reset_c
      endif

      ! missing: q within k hops in H, but further than dh+excess in C or not connected at all
      if (act0(p)) then
        if (C(node_i(p),node_j(p))) call bfs_c(C, node_i(p), node_j(p), dlim)
        do n=1,nqueue_h
          q = queue_h(n)
          dh = dist_h(q)
          if (dh.gt.k_ocn_topo .or. q.le.p .or. .not.act0(q)) cycle
          dc = d_inf
          if (C(node_i(p),node_j(p))) then
            if (dist_c(node_i(q),node_j(q)).ge.0) dc = dist_c(node_i(q),node_j(q))
          endif
          if (dc.eq.d_inf .or. dc.gt.dh+n_excess_ocn_topo) call add_flag(kind_missing, p, q, dc, dh)
        enddo
        call reset_c
      endif

      call reset_h
    enddo

    deallocate(act, act0)

  end subroutine check


  ! a node takes part in the checks if it is large enough, its cell is ocean, and it is either the main ocean body of
  ! its cell or connects to another ocean cell
  subroutine active(Cm, act)

    implicit none

    logical, intent(in) :: Cm(:,:)
    logical, intent(out) :: act(:)

    integer :: p, n, w

    act = .false.
    do p=1,nnode
      if (.not.(sig(p) .and. Cm(node_i(p),node_j(p)))) cycle
      if (p.eq.cell_main(node_i(p),node_j(p))) then
        act(p) = .true.
        cycle
      endif
      do n=adj0(p),adj0(p+1)-1
        w = adj(n)
        if ((node_i(w).ne.node_i(p) .or. node_j(w).ne.node_j(p)) .and. Cm(node_i(w),node_j(w))) then
          act(p) = .true.
          exit
        endif
      enddo
    enddo

  end subroutine active


  ! one flag per pair of cells and kind
  subroutine add_flag(kind, p, q, dc, dh)

    implicit none

    integer, intent(in) :: kind, p, q, dc, dh
    integer :: n

    do n=1,nflag
      if (flag_kind(n).ne.kind) cycle
      if (same_cells(flag_p(n), flag_q(n), p, q)) return
    enddo
    if (nflag.eq.size(flag_kind)) then
      call grow(flag_kind) ; call grow(flag_p) ; call grow(flag_q) ; call grow(flag_dc) ; call grow(flag_dh)
    endif
    nflag = nflag+1
    flag_kind(nflag) = kind
    flag_p(nflag) = p
    flag_q(nflag) = q
    flag_dc(nflag) = dc
    flag_dh(nflag) = dh

  end subroutine add_flag


  logical function same_cells(p1, q1, p2, q2)

    integer, intent(in) :: p1, q1, p2, q2

    same_cells = (same_cell(p1,p2) .and. same_cell(q1,q2)) .or. (same_cell(p1,q2) .and. same_cell(q1,p2))

  end function same_cells


  logical function same_cell(p, q)

    integer, intent(in) :: p, q

    same_cell = node_i(p).eq.node_i(q) .and. node_j(p).eq.node_j(q)

  end function same_cell


  ! number of flags not ignored, of one kind or of all kinds (kind=0)
  integer function nlive(kind)

    integer, intent(in) :: kind
    integer :: n

    nlive = 0
    do n=1,nflag
      if (kind.ne.0 .and. flag_kind(n).ne.kind) cycle
      if (is_ignored(flag_p(n),flag_q(n))) cycle
      nlive = nlive+1
    enddo

  end function nlive


  ! number of flags of one kind, including the ignored ones
  integer function nall(kind)

    integer, intent(in) :: kind

    nall = count(flag_kind(1:nflag).eq.kind)

  end function nall


  logical function is_flagged(p, q)

    integer, intent(in) :: p, q
    integer :: n

    is_flagged = .false.
    do n=1,nflag
      if (flag_p(n).eq.p .and. flag_q(n).eq.q) then
        is_flagged = .true.
        return
      endif
    enddo

  end function is_flagged


  logical function is_ignored(p, q)

    integer, intent(in) :: p, q
    integer :: n

    is_ignored = .false.
    do n=1,nignored
      if (ignored_p(n).eq.p .and. ignored_q(n).eq.q) then
        is_ignored = .true.
        return
      endif
    enddo

  end function is_ignored


  subroutine ignore(p, q)

    integer, intent(in) :: p, q

    if (nignored.ge.maxit) return
    nignored = nignored+1
    ignored_p(nignored) = p
    ignored_q(nignored) = q

  end subroutine ignore


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! breadth-first searches
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine bfs_h(p, dmax)

    integer, intent(in) :: p, dmax
    integer :: h, v, n, w

    dist_h(p) = 0
    nqueue_h = 1
    queue_h(1) = p
    h = 1
    do while (h.le.nqueue_h)
      v = queue_h(h)
      h = h+1
      if (dist_h(v).ge.dmax) cycle
      do n=adj0(v),adj0(v+1)-1
        w = adj(n)
        if (dist_h(w).ge.0) cycle
        dist_h(w) = dist_h(v)+1
        nqueue_h = nqueue_h+1
        queue_h(nqueue_h) = w
      enddo
    enddo

  end subroutine bfs_h


  subroutine reset_h

    integer :: n

    do n=1,nqueue_h
      dist_h(queue_h(n)) = -1
    enddo
    nqueue_h = 0

  end subroutine reset_h


  subroutine bfs_c(C, i0, j0, dmax)

    logical, intent(in) :: C(:,:)
    integer, intent(in) :: i0, j0, dmax
    integer :: h, i, j, n, ii, jj

    nqueue_c = 0
    if (.not.C(i0,j0)) return
    dist_c(i0,j0) = 0
    nqueue_c = 1
    queue_c(:,1) = [i0, j0]
    h = 1
    do while (h.le.nqueue_c)
      i = queue_c(1,h)
      j = queue_c(2,h)
      h = h+1
      if (dist_c(i,j).ge.dmax) cycle
      do n=1,4
        call nbr(i, j, n, ii, jj)
        if (jj.lt.1 .or. jj.gt.nj) cycle
        if (.not.C(ii,jj) .or. dist_c(ii,jj).ge.0) cycle
        dist_c(ii,jj) = dist_c(i,j)+1
        nqueue_c = nqueue_c+1
        queue_c(:,nqueue_c) = [ii, jj]
      enddo
    enddo

  end subroutine bfs_c


  subroutine reset_c

    integer :: n

    do n=1,nqueue_c
      dist_c(queue_c(1,n),queue_c(2,n)) = -1
    enddo
    nqueue_c = 0

  end subroutine reset_c


  ! neighbour n (west, east, south, north) of cell (i,j), periodic in longitude
  subroutine nbr(i, j, n, ii, jj)

    integer, intent(in) :: i, j, n
    integer, intent(out) :: ii, jj

    ii = i
    jj = j
    select case (n)
    case (1)
      ii = i-1
      if (ii.lt.1) ii = ni
    case (2)
      ii = i+1
      if (ii.gt.ni) ii = 1
    case (3)
      jj = j-1
    case (4)
      jj = j+1
    end select

  end subroutine nbr


  ! ocean cells connected to the ocean origin
  subroutine fill(keep, C)

    logical, intent(in) :: keep(:,:)
    logical, intent(out) :: C(:,:)

    real(wp) :: f(ni,nj)
    integer :: mask_connect(ni,nj)

    f = merge(1._wp, 0._wp, keep)
    call fill_ocean_lowres(f, lon_c, lat_c, lon_origin, lat_origin, mask_connect)
    C = mask_connect.eq.1

  end subroutine fill


  ! shortest coarse path from cell a to cell b
  subroutine path_c(C, ia, ja, ib, jb, ci, cj, ncell)

    logical, intent(in) :: C(:,:)
    integer, intent(in) :: ia, ja, ib, jb
    integer, intent(out) :: ci(:), cj(:)
    integer, intent(out) :: ncell

    integer :: i, j, n, ii, jj

    ncell = 0
    call bfs_c(C, ib, jb, d_inf)
    if (dist_c(ia,ja).ge.0) then
      i = ia
      j = ja
      ncell = 1
      ci(1) = i
      cj(1) = j
      do while (dist_c(i,j).gt.0)
        do n=1,4
          call nbr(i, j, n, ii, jj)
          if (jj.lt.1 .or. jj.gt.nj) cycle
          if (dist_c(ii,jj).eq.dist_c(i,j)-1) exit
        enddo
        i = ii
        j = jj
        ncell = ncell+1
        ci(ncell) = i
        cj(ncell) = j
      enddo
    endif
    call reset_c

  end subroutine path_c


  ! high resolution path from node p to node q (or, if the cell of p is not connected to the ocean, to the nearest
  ! node in a connected cell) that re-opens the least land: entering a node costs one minus the ocean fraction of its
  ! cell if that cell is closed, plus a small cost per hop
  subroutine path_min_land(keep, C, mask_fix, f_ocn, p, q, path, npath)

    logical, intent(in) :: keep(:,:), C(:,:), mask_fix(:,:)
    real(wp), intent(in) :: f_ocn(:,:)
    integer, intent(in) :: p, q
    integer, intent(out) :: path(:)
    integer, intent(out) :: npath

    real(wp), allocatable :: dist(:)
    integer, allocatable :: par(:)
    logical, allocatable :: done(:)
    logical :: to_ocean
    integer :: v, w, n, hit
    real(wp) :: d, dmin

    allocate(dist(nnode), par(nnode), done(nnode))
    dist = huge(1._wp)
    par = 0
    done = .false.
    dist(p) = 0._wp
    par(p) = p
    to_ocean = .not.C(node_i(p),node_j(p))
    hit = 0

    do
      ! next node, linear search (the graph has a few thousand nodes and this is called a few times)
      v = 0
      dmin = huge(1._wp)
      do n=1,nnode
        if (.not.done(n) .and. dist(n).lt.dmin) then
          dmin = dist(n)
          v = n
        endif
      enddo
      if (v.eq.0) exit
      done(v) = .true.
      if (to_ocean) then
        if (C(node_i(v),node_j(v))) then
          hit = v
          exit
        endif
      else if (v.eq.q) then
        hit = v
        exit
      endif
      do n=adj0(v),adj0(v+1)-1
        w = adj(n)
        if (done(w) .or. mask_fix(node_i(w),node_j(w))) cycle
        d = dist(v) + 0.01_wp
        if (.not.keep(node_i(w),node_j(w))) d = d + 1._wp-f_ocn(node_i(w),node_j(w))
        if (d.lt.dist(w)) then
          dist(w) = d
          par(w) = v
        endif
      enddo
    enddo

    npath = 0
    if (hit.gt.0) then
      v = hit
      do
        npath = npath+1
        path(npath) = v
        if (v.eq.p) exit
        v = par(v)
      enddo
    endif

    deallocate(dist, par, done)

  end subroutine path_min_land


  ! sort cells by increasing ocean fraction (insertion sort, stable)
  subroutine sort_cells(f_ocn, ci, cj, ncell)

    real(wp), intent(in) :: f_ocn(:,:)
    integer, intent(inout) :: ci(:), cj(:)
    integer, intent(in) :: ncell
    integer :: n, m, ti, tj

    do n=2,ncell
      ti = ci(n)
      tj = cj(n)
      m = n-1
      do while (m.ge.1)
        if (f_ocn(ci(m),cj(m)).le.f_ocn(ti,tj)) exit
        ci(m+1) = ci(m)
        cj(m+1) = cj(m)
        m = m-1
      enddo
      ci(m+1) = ti
      cj(m+1) = tj
    enddo

  end subroutine sort_cells


  ! order of the pruning: closures by decreasing, openings by increasing ocean fraction (key -f or f)
  subroutine sort_changes(f_ocn, ci, cj, ca, nchange)

    real(wp), intent(in) :: f_ocn(:,:)
    integer, intent(inout) :: ci(:), cj(:), ca(:)
    integer, intent(in) :: nchange
    integer :: n, m, ti, tj, ta

    do n=2,nchange
      ti = ci(n)
      tj = cj(n)
      ta = ca(n)
      m = n-1
      do while (m.ge.1)
        if (key(ci(m),cj(m),ca(m)).le.key(ti,tj,ta)) exit
        ci(m+1) = ci(m)
        cj(m+1) = cj(m)
        ca(m+1) = ca(m)
        m = m-1
      enddo
      ci(m+1) = ti
      cj(m+1) = tj
      ca(m+1) = ta
    enddo

  contains

    real(wp) function key(i, j, a)
      integer, intent(in) :: i, j, a
      if (a.eq.-1) then
        key = -f_ocn(i,j)
      else
        key = f_ocn(i,j)
      endif
    end function key

  end subroutine sort_changes

end module ocn_topology_mod
