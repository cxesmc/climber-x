!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : l n d _ g r i d
!
!  Purpose : land model grid
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit and Andrey Ganopolski
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
module lnd_grid

  use precision, only : wp
  use nml
  use control, only : out_dir
  use climber_grid, only : ni, nj, lat

  implicit none

     integer,  parameter :: nx = ni
     integer,  parameter :: ny = nj

     ! soil levels, read from lnd_par.nml
     ! the layer interface depths z_int are prescribed, the nodes z are the layer midpoints
     ! note: src/lndvc/lndvc_grid.f90, src/smb/smb_grid.f90 and src/bmb/bmb_grid.f90 keep their
     ! own hardcoded vertical columns and are not affected by the settings here
     integer :: nl
     real(wp), allocatable, dimension(:) :: z          ! z(0) is overwritten with -0.5*h_snow in the local copies
     real(wp), allocatable, dimension(:) :: z_int, dz, rdz, rdz_pos, rdz_neg
     ! levels for lakes, bottom layer adjusted for lake depth
     integer :: nl_l
     real(wp), allocatable, dimension(:) :: z_l        ! z_l(0) is overwritten with -0.5*h_snow in the local copies
     real(wp), allocatable, dimension(:) :: z_int_l, dz_l, rdz_l, rdz_pos_l, rdz_neg_l

     ! levels for carbon (including burial layer)
     integer :: nlc
     real(wp), allocatable, dimension(:) :: z_c
     real(wp), allocatable, dimension(:) :: z_int_c, dz_c, rdz_c, rdz_pos_c, rdz_neg_c

     ! surface types
     integer,  parameter :: nsurf = 8 ! 5 pfts + bare soil + lake + ice
     integer,  parameter, dimension(nsurf) :: flag_pft = (/1,1,1,1,1,0,0,0/)
     integer,  parameter, dimension(nsurf) :: flag_veg = (/1,1,1,1,1,1,0,0/)
     integer,  parameter, dimension(nsurf) :: i_surf = (/1,2,3,4,5,6,7,8/)
     integer,  parameter :: i_bare = 6
     integer,  parameter :: i_lake = 7
     integer,  parameter :: i_ice = 8
     ! PFTs
     integer,  parameter :: npft = 5
     integer,  parameter :: ntrees = 2
     integer,  parameter :: ngrass = 2
     integer,  parameter :: nshrub = 1
     integer,  parameter, dimension(npft) :: i_pft = (/1,2,3,4,5/)
     integer,  parameter, dimension(ntrees) :: i_trees = (/1,2/)
     integer,  parameter, dimension(ngrass) :: i_grass = (/3,4/)
     integer,  parameter, dimension(nshrub) :: i_shrub = (/5/)
     integer,  parameter, dimension(npft) :: flag_tree  = (/1,1,0,0,0/)
     integer,  parameter, dimension(npft) :: flag_grass = (/0,0,1,1,0/)
     integer,  parameter, dimension(npft) :: flag_shrub = (/0,0,0,0,1/)
     integer,  parameter :: nveg = npft+1
     ! carbon
     integer,  parameter :: ncarb = 5
     integer,  parameter :: ic_min = 1
     integer,  parameter :: ic_peat = 2
     integer,  parameter :: ic_shelf = 3
     integer,  parameter :: ic_ice = 4
     integer,  parameter :: ic_lake = 5
     ! ground/snow types
     integer,  parameter :: nsoil = 3  ! vegetation, ice, lake
     integer,  parameter, dimension(nsoil) :: i_soil = (/1,2,3/)
     integer,  parameter :: is_veg = 1
     integer,  parameter :: is_ice = 2
     integer,  parameter :: is_lake = 3


contains

  subroutine lnd_grid_init

  implicit none

  integer :: k
  character (len=256) :: fnm


    ! read the vertical grid from the namelist.
    ! note that lnd_grid_init is called before lnd_params_init, so the namelist is read
    ! here directly, following the same approach used by atm_grid_init
    fnm = trim(out_dir)//"/lnd_par.nml"
    call nml_read(fnm,"lnd_par","nl",nl)
    call nml_read(fnm,"lnd_par","nl_lake",nl_l)
    nlc = nl+1

    ! allocate and initialise the vertical grid arrays.
    ! zero initialisation matters: rdz_pos(0), rdz_neg(0) and their lake/carbon counterparts
    ! are never assigned below but are used by some of the local grid copies
    allocate(z(0:nl), z_int(0:nl), dz(0:nl), rdz(0:nl), rdz_pos(0:nl), rdz_neg(0:nl))
    allocate(z_l(0:nl_l), z_int_l(0:nl_l), dz_l(0:nl_l), rdz_l(0:nl_l), rdz_pos_l(0:nl_l), rdz_neg_l(0:nl_l))
    allocate(z_c(0:nlc), z_int_c(0:nlc), dz_c(0:nlc), rdz_c(0:nlc), rdz_pos_c(0:nlc), rdz_neg_c(0:nlc))
    z     = 0._wp; z_int     = 0._wp; dz     = 0._wp; rdz     = 0._wp; rdz_pos     = 0._wp; rdz_neg     = 0._wp
    z_l   = 0._wp; z_int_l   = 0._wp; dz_l   = 0._wp; rdz_l   = 0._wp; rdz_pos_l   = 0._wp; rdz_neg_l   = 0._wp
    z_c   = 0._wp; z_int_c   = 0._wp; dz_c   = 0._wp; rdz_c   = 0._wp; rdz_pos_c   = 0._wp; rdz_neg_c   = 0._wp

    ! depth of the vertical layer interfaces, nl+1 (nl_lake+1) values starting with 0
    call nml_read(fnm,"lnd_par","z_int_soil",z_int(0:nl))
    call nml_read(fnm,"lnd_par","z_int_lake",z_int_l(0:nl_l))

    call check_levels("z_int_soil",nl,z_int)
    call check_levels("z_int_lake",nl_l,z_int_l)

    ! for soil
    call make_levels(nl,z_int,z,dz,rdz,rdz_pos,rdz_neg)

    print *,'z',z
    print *,'z_int',z_int
    print *,'dz',dz

    ! for lakes
    call make_levels(nl_l,z_int_l,z_l,dz_l,rdz_l,rdz_pos_l,rdz_neg_l)

    print *,'z_l',z_l
    print *,'z_int_l',z_int_l
    print *,'dz_l',dz_l

    ! for carbon, the soil column plus one burial layer of thickness dz(nl) at the bottom
    z_int_c(0:nl) = z_int(0:nl)
    z_int_c(nlc)  = z_int(nl) + dz(nl)
    call make_levels(nlc,z_int_c,z_c,dz_c,rdz_c,rdz_pos_c,rdz_neg_c)

    print *,'z_c',z_c
    print *,'z_int_c',z_int_c
    print *,'dz_c',dz_c

  return

  end subroutine lnd_grid_init


  ! check that the layer interface depths read from the namelist are usable
  subroutine check_levels(name,n,zi)

  implicit none

  character (len=*), intent(in) :: name
  integer, intent(in) :: n
  real(wp), dimension(0:), intent(in) :: zi

  integer :: k


    if (n.lt.2) then
      print *,'ERROR: at least 2 layers are required, got ',n,' for ',trim(name)
      stop
    endif
    if (zi(0).ne.0._wp) then
      print *,'ERROR: the first value of ',trim(name),' must be 0., got ',zi(0)
      stop
    endif
    do k=1,n
      if (zi(k).le.zi(k-1)) then
        print *,'ERROR: ',trim(name),' must be strictly increasing, but level ',k,' = ',zi(k), &
                ' is not larger than level ',k-1,' = ',zi(k-1)
        print *,'       ',trim(name),' needs ',n+1,' increasing values, the first one being 0.'
        stop
      endif
    enddo

  return

  end subroutine check_levels


  ! derive layer thicknesses, node depths and reciprocals from the layer interface depths
  subroutine make_levels(n,zi,zz,dzz,rdzz,rdzz_pos,rdzz_neg)

  implicit none

  integer, intent(in) :: n
  real(wp), dimension(0:), intent(in)  :: zi
  real(wp), dimension(0:), intent(out) :: zz, dzz, rdzz, rdzz_pos, rdzz_neg

  integer :: k


    zz  = 0._wp
    dzz = 0._wp
    rdzz = 0._wp
    rdzz_pos = 0._wp
    rdzz_neg = 0._wp

    ! layer thicknesses and node depths (layer midpoints)
    do k=1,n
     dzz(k) = zi(k) - zi(k-1)
     zz(k)  = 0.5_wp * ( zi(k-1) + zi(k) )
    enddo

    ! reciprocals to speed up fortran
    rdzz(1:n) = 1._wp/dzz(1:n)
    do k=1,n-1
     rdzz_pos(k) = 1._wp/(zz(k+1)-zz(k))
    enddo
    rdzz_pos(n) = 0._wp
    do k=1,n
     rdzz_neg(k) = 1._wp/(zz(k)-zz(k-1))
    enddo

  return

  end subroutine make_levels


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s o i l _ a t _ d e p t h
  ! Purpose  :  linearly interpolate a soil profile given at the layer nodes
  !             z(1:nl) to an arbitrary depth, so that parameterisations can refer
  !             to a fixed physical depth instead of to a layer index
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function soil_at_depth(prof,depth) result(val)

  implicit none

  real(wp), dimension(:), intent(in) :: prof    ! soil profile, prof(k) is the value at node z(k), k=1..nl
  real(wp), intent(in) :: depth                 ! m, depth to interpolate to
  real(wp) :: val

  integer :: k
  real(wp) :: w


    if (depth.le.z(1)) then
      ! above the first node, no extrapolation towards the surface
      val = prof(1)
    else if (depth.ge.z(nl)) then
      ! below the last node, no extrapolation towards the bottom
      val = prof(nl)
    else
      val = prof(nl)
      do k=1,nl-1
        if (depth.le.z(k+1)) then
          w = (depth-z(k)) / (z(k+1)-z(k))
          val = (1._wp-w)*prof(k) + w*prof(k+1)
          exit
        endif
      enddo
    endif

  return

  end function soil_at_depth

end module lnd_grid



