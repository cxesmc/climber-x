module lndvc_model
    ! Driver / orchestrator for the virtual-cell (lndvc) framework.
    !
    ! Owns the per-domain control (init/update/end), the loop over leaf virtual
    ! cells, dispatch to single-column physics by surface class, and the
    ! reduction back to coarse-cell state for coupling. Physics is ported from
    ! lnd/smb as single-column, single-class modules and called from here.
    ! See docs/design/virtual-cells.md.
    !
    ! Phase 1: allocation + decomposition + aggregation plumbing, wired into
    ! climber.x behind flag_lndvc (default off). Physics dispatch is a no-op
    ! until the port (Phase 2); numerical identity with the reference land model
    ! is therefore gated on Phase 2.

    use ncio
    use precision, only : sp, dp, wp

    use lndvc_def
    use lndvc_grid
    use lndvc_decomp,    only : lndvc_decompose, lndvc_remap_state
    use lndvc_downscale, only : lndvc_downscale_forcing
    use lndvc_aggregate, only : lndvc_aggregate_weights, lndvc_aggregate_cell, lndvc_conservation_check

    ! ported single-column physics
    use smb_snow_m,      only : snow_update

    implicit none

    private
    public :: lndvc_init
    public :: lndvc_update
    public :: lndvc_end

contains

    ! === Whole-domain control =================================================

    subroutine lndvc_init(lnd, nx, ny)
        implicit none
        type(lndvc_class), intent(inout) :: lnd
        integer,           intent(in)    :: nx, ny

        call lndvc_grid_init()
        call lndvc_alloc(lnd, nx, ny)

        return
    end subroutine lndvc_init

    subroutine lndvc_update(lnd)
        ! One coupling step: update every leaf virtual cell, then aggregate to
        ! coarse-cell state for coupling. The decomposition is assumed current
        ! (refreshed by the coupler via lndvc_decompose before this call).

        implicit none

        type(lndvc_class), intent(inout) :: lnd

        integer :: n, i, j, k
        real(wp), allocatable :: wt(:)

        allocate(wt(lnd%n_vc))

        ! TODO: OMP over the compressed active-leaf list (lnd%leaf_1d).
        do n = 1, lnd%ncells
            i = lnd%ij_1d(1,n)
            j = lnd%ij_1d(2,n)
            do k = 1, lnd%n_vc
                if (lnd%vc(i,j,k)%desc%class == 0) cycle
                call lndvc_update_vc(lnd%vc(i,j,k))
            end do
            call lndvc_aggregate_weights(lnd%vc(i,j,:), wt)
            call lndvc_aggregate_cell(lnd%cell(i,j), lnd%vc(i,j,:), wt)
        end do

        deallocate(wt)

        call lndvc_conservation_check(lnd)

        return

    end subroutine lndvc_update

    subroutine lndvc_end(lnd)
        implicit none
        type(lndvc_class), intent(inout) :: lnd
        call lndvc_dealloc(lnd)
        return
    end subroutine lndvc_end

    ! === Allocation ===========================================================

    subroutine lndvc_alloc(lnd, nx, ny)
        implicit none
        type(lndvc_class), intent(inout) :: lnd
        integer,           intent(in)    :: nx, ny

        lnd%n_vc   = 4          ! identity baseline: land, lake, ice (+shelf)
        lnd%ncells = 0

        if (allocated(lnd%vc))     deallocate(lnd%vc)
        if (allocated(lnd%cell))   deallocate(lnd%cell)
        if (allocated(lnd%id_map)) deallocate(lnd%id_map)
        if (allocated(lnd%ij_1d))  deallocate(lnd%ij_1d)

        allocate(lnd%vc(nx,ny,lnd%n_vc))
        allocate(lnd%cell(nx,ny))
        allocate(lnd%id_map(nx,ny))
        allocate(lnd%ij_1d(2,nx*ny))

        lnd%id_map = 0
        lnd%ij_1d  = 0

        return
    end subroutine lndvc_alloc

    subroutine lndvc_dealloc(lnd)
        implicit none
        type(lndvc_class), intent(inout) :: lnd

        if (allocated(lnd%vc))      deallocate(lnd%vc)
        if (allocated(lnd%cell))    deallocate(lnd%cell)
        if (allocated(lnd%id_map))  deallocate(lnd%id_map)
        if (allocated(lnd%ij_1d))   deallocate(lnd%ij_1d)
        if (allocated(lnd%z_edges)) deallocate(lnd%z_edges)
        if (allocated(lnd%leaf_1d)) deallocate(lnd%leaf_1d)
        if (allocated(lnd%z_lake))  deallocate(lnd%z_lake)

        return
    end subroutine lndvc_dealloc

    ! === Per-vc dispatch ======================================================

    subroutine lndvc_update_vc(vc)
        ! Update a single virtual cell: downscale forcing to its elevation, then
        ! dispatch to the single-column physics for its surface class.

        implicit none

        type(vc_t), intent(inout) :: vc

        call lndvc_downscale_forcing(vc%forc, vc%desc)

        select case(vc%desc%class)
            case(1)   ! land
                call lndvc_update_land(vc)
            case(2)   ! lake
                call lndvc_update_lake(vc)
            case(3)   ! ice
                call lndvc_update_ice(vc)
        end select

        return

    end subroutine lndvc_update_vc

    subroutine lndvc_update_land(vc)
        implicit none
        type(vc_t), intent(inout) :: vc
        ! TODO: ported single-column land physics (veg energy balance, soil
        !       thermal+permafrost, hydrology, snow, soil carbon).
        return
    end subroutine lndvc_update_land

    subroutine lndvc_update_lake(vc)
        implicit none
        type(vc_t), intent(inout) :: vc
        ! TODO: ported single-column lake physics.
        return
    end subroutine lndvc_update_lake

    subroutine lndvc_update_ice(vc)
        ! Single-column ice-surface update. Blocks are unpacked into the ported
        ! scalar-signature physics routines (pattern A). Being filled in as the
        ! ice/SMB path is ported file-by-file.

        implicit none

        type(vc_t), intent(inout) :: vc

        real(wp) :: evp

        ! TODO: evp (sublimation) from the ported ice-surface energy balance
        evp = 0._wp

        if (allocated(vc%snow)) then
            call snow_update(vc%snow%mask_snow, evp, &
                             vc%snow%w_snow, vc%snow%w_snow_old, vc%snow%w_snow_max, &
                             vc%snow%h_snow)
        end if

        ! TODO: ice/firn temperature, surface energy balance, SMB mass-budget diagnostic.

        return
    end subroutine lndvc_update_ice

end module lndvc_model
