module lndvc_decomp
    ! Decomposition service for the virtual-cell framework.
    !
    ! Turns each coarse coupler cell into a list of leaf virtual cells via a
    ! composable pipeline of refine stages (see docs/design/virtual-cells.md, sec 4.3):
    !
    !   coarse cell --[optional bilinear regrid]--> mid-res cell
    !               --[hypsometric split]--------> virtual cell {z,dz,w,slope,class}
    !
    ! Also owns the conservative state-remap that runs when the decomposition
    ! changes (ice advance/retreat, topography shift).
    !
    ! Physics is grid-agnostic: all horizontal coupling lives here, never in the
    ! column physics.

    use precision, only : wp
    use lndvc_def
    use lndvc_grid

    implicit none

    private
    public :: lndvc_decompose        ! build the leaf virtual-cell list for a coarse cell
    public :: lndvc_remap_state       ! conservative state transfer when decomposition changes

contains

    subroutine lndvc_decompose(vc)
        ! Build the leaf virtual-cell descriptors for the domain from the
        ! high-res reference topography (hypsometry) and the optional mid-res
        ! horizontal refine stage. Sets z, dz, area weight w (Sum w = 1 per
        ! coarse cell), slope, and surface class for each leaf.

        implicit none

        type(lndvc_class), intent(inout) :: vc

        ! TODO: (1) optional bilinear regrid coarse->mid-res grid (coords/map_field)
        !       (2) per-(mid-res) cell hypsometry from high-res reference field
        !       (3) bin into elevation bands -> leaf vc descriptors + area weights
        !       (4) build compressed active-leaf list for OMP (generalizes ij_1d)

        return

    end subroutine lndvc_decompose

    subroutine lndvc_remap_state(vc)
        ! Conservatively transfer prognostic state (carbon, heat, snow, water)
        ! between leaf virtual cells when band membership or class changes.
        ! Replaces the ad-hoc "initialize newly vegetated/ice/lake cell" branches
        ! scattered through lnd_update today.

        implicit none

        type(lndvc_class), intent(inout) :: vc

        ! TODO: detect decomposition change; conservative redistribution.

        return

    end subroutine lndvc_remap_state

end module lndvc_decomp
