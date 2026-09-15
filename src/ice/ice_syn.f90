!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module : i c e _ s y n
!
!  Purpose : Synthetic ice-sheet "model". Builds an ice-sheet geometry
!            on a regional ice grid from a prescribed (optionally
!            transient) target-extent mask: signed distance to the mask
!            boundary -> surface-elevation profile (perfect-plastic or
!            linear, see ice_syn_topo) -> grounded ice thickness on the
!            current bedrock. Stateless: the geometry is a function of
!            (target mask, bedrock, sea level) at each call.
!
!            Selected with ice_model_name = 'syn' in control.nml.
!            Parameters in ice_syn_par.nml (group &ice_syn).
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2024 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit and Reinhard Calov
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
module ice_syn

    use precision, only : wp, dp
    use constants, only : rho_i, rho_sw, map_gen
    use coords, only : grid_class, grid_init, map_class, map_init, map_field
    use ice_syn_topo, only : compute_signed_distance, compute_z_syn_linear, compute_z_syn_plastic
    use nml
    use ncio

    implicit none
    private

    public :: ice_syn_class
    public :: ice_syn_init
    public :: ice_syn_update
    public :: ice_syn_write_step
    public :: ice_syn_end

    !-----------------------------------------------------------------
    ! Parameters (namelist group &ice_syn)
    !-----------------------------------------------------------------
    type :: ice_syn_par_type
        ! target-extent forcing
        character(len=512) :: mask_file = ""      ! lat/lon file with variable mask_var(lon,lat[,time])
        character(len=56)  :: mask_var  = "mask"  ! variable name
        real(wp)           :: var_thresh = 0.5_wp ! cells with mask_var > var_thresh are target ice
                                                  ! (0.5 for 0/1 masks, e.g. 10. for thickness fields)
        ! synthetic-elevation profile
        logical  :: use_plastic = .true.          ! T => plastic Nye/Vialov; F => linear wedge
        real(wp) :: slope       = 6.09e-3_wp      ! m/m, linear inside (use_plastic=F)
        real(wp) :: tau0        = 36.0e3_wp       ! Pa, plastic basal yield stress
        real(wp) :: slope_out   = 6.09e-3_wp      ! m/m, linear outside ramp
        real(wp) :: z_max_in    = 2500.0_wp       ! m, inside cap
        real(wp) :: z_max_out   = 200.0_wp        ! m, outside floor
        real(wp) :: rho_ice     = 910.0_wp        ! kg/m3
        real(wp) :: g           = 9.81_wp         ! m/s2
        ! grounded-only constraint: inside the target mask ice is at least
        ! thick enough to be grounded (H >= (z_sl-z_bed)*rho_sw/rho_i + h_grd_min)
        real(wp) :: h_grd_min   = 10.0_wp         ! m
    end type ice_syn_par_type

    !-----------------------------------------------------------------
    ! Per-domain state
    !-----------------------------------------------------------------
    type :: ice_syn_class
        type(grid_class)       :: grid            ! ice grid
        type(ice_syn_par_type) :: par

        ! target-mask forcing (file grid, lat/lon)
        type(grid_class)       :: grid_file
        type(map_class)        :: maps_file_to_ice
        integer                :: ni_file, nj_file
        logical                :: l_transient      ! file has a time dimension
        integer                :: ntime
        real(wp), allocatable  :: time_file(:)
        integer                :: i0, i1
        real(wp), allocatable  :: f_ice_0(:,:)     ! target ice fraction on ice grid, slice i0
        real(wp), allocatable  :: f_ice_1(:,:)     ! target ice fraction on ice grid, slice i1

        ! geometry on the ice grid
        real(wp), allocatable  :: x_m(:,:), y_m(:,:)   ! coordinates [m]
        real(wp), allocatable  :: f_ice_target(:,:)    ! time-interpolated target ice fraction [1]
        logical,  allocatable  :: mask_target(:,:)     ! target ice mask
        real(wp), allocatable  :: d_m(:,:)             ! signed distance to mask boundary [m]
        real(wp), allocatable  :: z_syn(:,:)           ! synthetic surface elevation [m]
        real(wp), allocatable  :: H_ice(:,:)           ! ice thickness [m]
        real(wp), allocatable  :: z_sur(:,:)           ! surface elevation [m]

        ! output
        character(len=1024)    :: file2D
        integer                :: nout
    end type ice_syn_class

contains

    !=================================================================
    ! init
    !=================================================================
    subroutine ice_syn_init(syn, grid, par_file, time, z_bed, z_sl, out_dir, file_prefix)

        type(ice_syn_class), intent(inout) :: syn
        type(grid_class),    intent(in)    :: grid
        character(len=*),    intent(in)    :: par_file
        real(wp),            intent(in)    :: time
        real(wp),            intent(in)    :: z_bed(:,:)
        real(wp),            intent(in)    :: z_sl(:,:)
        character(len=*),    intent(in)    :: out_dir
        character(len=*),    intent(in), optional :: file_prefix   ! default "ice_"

        integer :: nx, ny
        character(len=64) :: prefix

        syn%grid = grid
        nx = grid%G%nx
        ny = grid%G%ny

        call ice_syn_par_load(syn%par, par_file)

        allocate(syn%x_m(nx,ny), syn%y_m(nx,ny))
        allocate(syn%f_ice_target(nx,ny), syn%mask_target(nx,ny))
        allocate(syn%d_m(nx,ny), syn%z_syn(nx,ny), syn%H_ice(nx,ny), syn%z_sur(nx,ny))
        allocate(syn%f_ice_0(nx,ny), syn%f_ice_1(nx,ny))

        ! grid coordinates in metres (grid axes are in grid%cs%units, e.g. km)
        syn%x_m = real(grid%x,wp) * real(grid%cs%xy_conv,wp)
        syn%y_m = real(grid%y,wp) * real(grid%cs%xy_conv,wp)

        call ice_syn_mask_init(syn, time)

        call ice_syn_update(syn, time, z_bed, z_sl)

        ! 2D output file
        prefix = "ice_"
        if (present(file_prefix)) prefix = file_prefix
        syn%file2D = trim(out_dir)//"/"//trim(prefix)//trim(grid%name)//".nc"
        syn%nout = 0
        call ice_syn_write_init(syn)

    end subroutine ice_syn_init

    !=================================================================
    ! update: recompute geometry for current time / bedrock / sea level
    !=================================================================
    subroutine ice_syn_update(syn, time, z_bed, z_sl)

        type(ice_syn_class), intent(inout) :: syn
        real(wp),            intent(in)    :: time
        real(wp),            intent(in)    :: z_bed(:,:)
        real(wp),            intent(in)    :: z_sl(:,:)

        integer  :: i, j
        real(wp) :: h_grd

        ! target mask at current time
        call ice_syn_mask_update(syn, time)

        ! signed distance to mask boundary [m]
        call compute_signed_distance(syn%d_m, syn%mask_target, syn%x_m, syn%y_m, "m")

        ! synthetic surface elevation on top of the bedrock
        if (syn%par%use_plastic) then
            call compute_z_syn_plastic(syn%z_syn, syn%d_m, z_bed, syn%mask_target, &
                                       syn%par%tau0, syn%par%slope_out, syn%par%z_max_in, syn%par%z_max_out, &
                                       rho_ice=syn%par%rho_ice, g=syn%par%g)
        else
            call compute_z_syn_linear(syn%z_syn, syn%d_m, z_bed, syn%mask_target, &
                                      syn%par%slope, syn%par%z_max_in, syn%par%z_max_out)
        end if

        ! grounded ice thickness inside the target mask, no ice outside
        do j = 1, syn%grid%G%ny
            do i = 1, syn%grid%G%nx
                if (syn%mask_target(i,j)) then
                    h_grd = max(z_sl(i,j)-z_bed(i,j), 0.0_wp)*rho_sw/rho_i + syn%par%h_grd_min
                    syn%H_ice(i,j) = max(syn%z_syn(i,j)-z_bed(i,j), h_grd)
                else
                    syn%H_ice(i,j) = 0.0_wp
                end if
                syn%z_sur(i,j) = z_bed(i,j) + syn%H_ice(i,j)
            end do
        end do

    end subroutine ice_syn_update

    !=================================================================
    ! target-mask forcing
    !=================================================================

    subroutine ice_syn_mask_init(syn, time)

        type(ice_syn_class), intent(inout) :: syn
        real(wp),            intent(in)    :: time

        real(wp), allocatable :: lon(:), lat(:)
        character(len=256) :: fnm
        integer :: ppos, spos

        fnm = syn%par%mask_file
        if (len_trim(fnm) == 0) then
            error stop "ice_syn_mask_init: mask_file must be set in &ice_syn"
        end if

        syn%ni_file = nc_size(trim(fnm),"lon")
        syn%nj_file = nc_size(trim(fnm),"lat")
        allocate(lon(syn%ni_file), lat(syn%nj_file))
        call nc_read(trim(fnm),"lon",lon)
        call nc_read(trim(fnm),"lat",lat)

        spos = scan(trim(fnm),"/", BACK=.true.)+1
        ppos = scan(trim(fnm),".", BACK=.true.)-1
        call grid_init(syn%grid_file, name=trim(fnm(spos:ppos)), mtype="latlon", units="degrees", &
            x0=real(lon(1),dp), dx=real(lon(2)-lon(1),dp), nx=syn%ni_file, &
            y0=real(lat(1),dp), dy=real(lat(2)-lat(1),dp), ny=syn%nj_file)
        deallocate(lon, lat)

        call map_init(syn%maps_file_to_ice, syn%grid_file, syn%grid, method="con", &
            gen=map_gen, fldr="maps", load=.TRUE., clean=.FALSE.)

        syn%l_transient = nc_exists_var(trim(fnm),"time")
        if (syn%l_transient) then
            syn%ntime = nc_size(trim(fnm),"time")
            allocate(syn%time_file(syn%ntime))
            call nc_read(trim(fnm),"time",syn%time_file)
        else
            syn%ntime = 1
            allocate(syn%time_file(1))
            syn%time_file(1) = time
        end if

        ! force reading of both bracketing slices on first update
        syn%i0 = -1
        syn%i1 = -1

    end subroutine ice_syn_mask_init

    !-----------------------------------------------------------------
    !> Interpolate the target ice fraction to `time` (linear between
    !> bracketing slices; clamped outside the file range) and threshold
    !> it to the logical target mask.
    !-----------------------------------------------------------------
    subroutine ice_syn_mask_update(syn, time)

        type(ice_syn_class), intent(inout) :: syn
        real(wp),            intent(in)    :: time

        integer  :: i0, i1, imin, i0_old, i1_old
        real(wp) :: w0, w1

        if (.not. syn%l_transient) then
            i0 = 1; i1 = 1; w0 = 1.0_wp; w1 = 0.0_wp
        else if (time <= syn%time_file(1)) then
            i0 = 1; i1 = 1; w0 = 1.0_wp; w1 = 0.0_wp
        else if (time >= syn%time_file(syn%ntime)) then
            i0 = syn%ntime; i1 = syn%ntime; w0 = 1.0_wp; w1 = 0.0_wp
        else
            imin = minloc(abs(syn%time_file-time),1)
            if (syn%time_file(imin) <= time) then
                i0 = imin; i1 = imin+1
            else
                i0 = imin-1; i1 = imin
            end if
            w0 = 1.0_wp - (time-syn%time_file(i0))/(syn%time_file(i1)-syn%time_file(i0))
            w1 = 1.0_wp - w0
        end if

        ! (re)load slices only when the bracket changes
        i0_old = syn%i0
        i1_old = syn%i1
        if (i0 /= i0_old) then
            if (i0 == i1_old) then
                syn%f_ice_0 = syn%f_ice_1
            else
                call ice_syn_read_slice(syn, i0, syn%f_ice_0)
            end if
        end if
        if (i1 /= i1_old) then
            if (i1 == i0) then
                syn%f_ice_1 = syn%f_ice_0
            else
                call ice_syn_read_slice(syn, i1, syn%f_ice_1)
            end if
        end if
        syn%i0 = i0
        syn%i1 = i1

        syn%f_ice_target = w0*syn%f_ice_0 + w1*syn%f_ice_1
        syn%mask_target  = (syn%f_ice_target >= 0.5_wp)

    end subroutine ice_syn_mask_update

    !-----------------------------------------------------------------
    !> Read slice `idx` of mask_var, threshold to 0/1 on the file grid,
    !> conservatively map to the ice grid as an ice fraction.
    !-----------------------------------------------------------------
    subroutine ice_syn_read_slice(syn, idx, f_ice)

        type(ice_syn_class), intent(inout) :: syn
        integer,             intent(in)    :: idx
        real(wp),            intent(out)   :: f_ice(:,:)

        real(wp), allocatable :: tmp(:,:)

        allocate(tmp(syn%ni_file, syn%nj_file))
        if (syn%l_transient) then
            call nc_read(trim(syn%par%mask_file), trim(syn%par%mask_var), tmp, &
                         start=[1,1,idx], count=[syn%ni_file,syn%nj_file,1])
        else
            call nc_read(trim(syn%par%mask_file), trim(syn%par%mask_var), tmp)
        end if
        where (tmp > syn%par%var_thresh)
            tmp = 1.0_wp
        elsewhere
            tmp = 0.0_wp
        end where

        f_ice = 0.0_wp
        call map_field(syn%maps_file_to_ice, "f_ice", tmp, f_ice, stat="mean", missing_value=-9999._dp)
        where (f_ice < 0.0_wp) f_ice = 0.0_wp

        deallocate(tmp)

    end subroutine ice_syn_read_slice

    !=================================================================
    ! output
    !=================================================================

    subroutine ice_syn_write_init(syn)

        type(ice_syn_class), intent(inout) :: syn

        real(wp), parameter :: empty_time(1) = [-9999.0_wp]   ! placeholder, unlimited

        call nc_create(syn%file2D)
        call nc_write_dim(syn%file2D, "xc", x=syn%grid%G%x, axis="x", units=trim(syn%grid%cs%units))
        call nc_write_dim(syn%file2D, "yc", x=syn%grid%G%y, axis="y", units=trim(syn%grid%cs%units))
        call nc_write_dim(syn%file2D, "time", x=empty_time, axis="t", units="years", unlimited=.TRUE.)
        call nc_write(syn%file2D, "lon", real(syn%grid%lon,wp), dims=["xc","yc"], long_name="longitude", units="degrees_east")
        call nc_write(syn%file2D, "lat", real(syn%grid%lat,wp), dims=["xc","yc"], long_name="latitude", units="degrees_north")

    end subroutine ice_syn_write_init

    subroutine ice_syn_write_step(syn, time, z_bed, z_sl, smb)

        type(ice_syn_class), intent(inout) :: syn
        real(wp),            intent(in)    :: time
        real(wp),            intent(in)    :: z_bed(:,:)
        real(wp),            intent(in)    :: z_sl(:,:)
        real(wp),            intent(in)    :: smb(:,:)

        integer :: nx, ny, n, ncid
        character(len=1024) :: fnm

        nx = syn%grid%G%nx
        ny = syn%grid%G%ny
        syn%nout = syn%nout + 1
        n = syn%nout
        fnm = syn%file2D

        call nc_open(fnm, ncid)
        call nc_write(fnm, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)
        call nc_write(fnm, "f_ice_target", real(syn%f_ice_target,wp), dims=["xc  ","yc  ","time"], start=[1,1,n], count=[nx,ny,1], &
            long_name="target ice fraction (time-interpolated)", units="1", ncid=ncid)
        call nc_write(fnm, "mask_target", merge(1,0,syn%mask_target), dims=["xc  ","yc  ","time"], start=[1,1,n], count=[nx,ny,1], &
            long_name="target ice mask", units="1", ncid=ncid)
        call nc_write(fnm, "dist", real(syn%d_m,wp), dims=["xc  ","yc  ","time"], start=[1,1,n], count=[nx,ny,1], &
            long_name="signed distance to target-mask boundary", units="m", ncid=ncid)
        call nc_write(fnm, "z_syn", real(syn%z_syn,wp), dims=["xc  ","yc  ","time"], start=[1,1,n], count=[nx,ny,1], &
            long_name="synthetic surface elevation (profile)", units="m", ncid=ncid)
        call nc_write(fnm, "H_ice", real(syn%H_ice,wp), dims=["xc  ","yc  ","time"], start=[1,1,n], count=[nx,ny,1], &
            long_name="ice thickness", units="m", ncid=ncid)
        call nc_write(fnm, "z_srf", real(syn%z_sur,wp), dims=["xc  ","yc  ","time"], start=[1,1,n], count=[nx,ny,1], &
            long_name="surface elevation", units="m", ncid=ncid)
        call nc_write(fnm, "z_bed", real(z_bed,wp), dims=["xc  ","yc  ","time"], start=[1,1,n], count=[nx,ny,1], &
            long_name="bedrock elevation", units="m", ncid=ncid)
        call nc_write(fnm, "z_sl", real(z_sl,wp), dims=["xc  ","yc  ","time"], start=[1,1,n], count=[nx,ny,1], &
            long_name="sea/lake level", units="m", ncid=ncid)
        call nc_write(fnm, "smb", real(smb,wp), dims=["xc  ","yc  ","time"], start=[1,1,n], count=[nx,ny,1], &
            long_name="surface mass balance received from smb model (diagnostic, not applied)", units="m(ice)/s", ncid=ncid)
        call nc_close(ncid)

    end subroutine ice_syn_write_step

    !=================================================================
    ! end
    !=================================================================
    subroutine ice_syn_end(syn)

        type(ice_syn_class), intent(inout) :: syn

        if (allocated(syn%x_m))          deallocate(syn%x_m)
        if (allocated(syn%y_m))          deallocate(syn%y_m)
        if (allocated(syn%f_ice_target)) deallocate(syn%f_ice_target)
        if (allocated(syn%mask_target))  deallocate(syn%mask_target)
        if (allocated(syn%d_m))          deallocate(syn%d_m)
        if (allocated(syn%z_syn))        deallocate(syn%z_syn)
        if (allocated(syn%H_ice))        deallocate(syn%H_ice)
        if (allocated(syn%z_sur))        deallocate(syn%z_sur)
        if (allocated(syn%f_ice_0))      deallocate(syn%f_ice_0)
        if (allocated(syn%f_ice_1))      deallocate(syn%f_ice_1)
        if (allocated(syn%time_file))    deallocate(syn%time_file)

    end subroutine ice_syn_end

    !=================================================================
    ! namelist
    !=================================================================
    subroutine ice_syn_par_load(par, filename)

        type(ice_syn_par_type), intent(inout) :: par
        character(len=*),       intent(in)    :: filename

        character(len=8), parameter :: group = "ice_syn"

        call nml_read(filename, group, "mask_file",   par%mask_file)
        call nml_read(filename, group, "mask_var",    par%mask_var)
        call nml_read(filename, group, "var_thresh",  par%var_thresh)
        call nml_read(filename, group, "use_plastic", par%use_plastic)
        call nml_read(filename, group, "slope",       par%slope)
        call nml_read(filename, group, "tau0",        par%tau0)
        call nml_read(filename, group, "slope_out",   par%slope_out)
        call nml_read(filename, group, "z_max_in",    par%z_max_in)
        call nml_read(filename, group, "z_max_out",   par%z_max_out)
        call nml_read(filename, group, "rho_ice",     par%rho_ice)
        call nml_read(filename, group, "g",           par%g)
        call nml_read(filename, group, "h_grd_min",   par%h_grd_min)

    end subroutine ice_syn_par_load

end module ice_syn
