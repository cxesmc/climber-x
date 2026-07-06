module lndvc_model
    ! Driver / orchestrator for the virtual-cell (lndvc) framework.
    !
    ! Owns the per-domain control (init/update/end), the loop over leaf virtual
    ! cells, dispatch to single-column physics by surface class, and the
    ! reduction back to coarse-cell state for coupling. Physics is ported from
    ! lnd/smb as single-column, single-class modules and called from here.
    ! See docs/design/virtual-cells.md.

    use ncio
    use precision, only : sp, dp, wp

    use lndvc_def
    use lndvc_grid
    use lndvc_decomp,    only : lndvc_decompose, lndvc_remap_state
    use lndvc_downscale, only : lndvc_downscale_forcing
    use lndvc_aggregate, only : lndvc_aggregate_weights, lndvc_aggregate_cell, lndvc_conservation_check

    implicit none

    private
    public :: lndvc_init
    public :: lndvc_update
    public :: lndvc_end

contains

    ! === Whole-domain control =================================================

    subroutine lndvc_init(lnd)
        implicit none
        type(lndvc_class), intent(inout) :: lnd
        ! TODO: allocate container, build fixed band edges, initial decomposition.
        call lndvc_decompose(lnd)
        return
    end subroutine lndvc_init

    subroutine lndvc_update(lnd)
        ! One coupling step: refresh decomposition weights, update every leaf
        ! virtual cell, then aggregate to coarse-cell state for coupling.

        implicit none

        type(lndvc_class), intent(inout) :: lnd

        integer :: n, i, j, k
        real(wp), allocatable :: wt(:)

        ! TODO: OMP over the compressed active-leaf list (lnd%leaf_1d).
        do n = 1, lnd%ncells
            i = lnd%ij_1d(1,n)
            j = lnd%ij_1d(2,n)
            do k = 1, lnd%n_vc
                call lndvc_update_vc(lnd%vc(i,j,k))
            end do
            allocate(wt(lnd%n_vc))
            call lndvc_aggregate_weights(lnd%vc(i,j,:), wt)
            call lndvc_aggregate_cell(lnd%cell(i,j), lnd%vc(i,j,:), wt)
            deallocate(wt)
        end do

        call lndvc_conservation_check(lnd)

        return

    end subroutine lndvc_update

    subroutine lndvc_end(lnd)
        implicit none
        type(lndvc_class), intent(inout) :: lnd
        ! TODO: deallocate container.
        return
    end subroutine lndvc_end

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
        implicit none
        type(vc_t), intent(inout) :: vc
        ! TODO: ported single-column ice-surface physics; SMB as mass-budget diagnostic.
        return
    end subroutine lndvc_update_ice

end module lndvc_model
