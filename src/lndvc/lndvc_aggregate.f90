module lndvc_aggregate
    ! Aggregation / conservation service for the virtual-cell framework.
    !
    ! Reduces the leaf virtual-cell ensemble back to a coarse-cell mean
    ! (vc_cell_t) for coupling with CLIMBER, using per-variable reduction rules
    ! (area-weighted mean / sum / flux-conserving). Enforces energy/water/carbon
    ! conservation at the coarse-cell boundary. See docs/design/virtual-cells.md, sec 5.
    !
    ! Conservation notes:
    !  - conservative fields (precip above all) must be renormalized after any
    !    non-conservative regrid so the coarse-cell mean is preserved;
    !  - the leaf reduction must be strictly area-conservative;
    !  - where a high-res backend is active, the coarse aggregate is the
    !    conservative reduction of the high-res field.

    use precision, only : wp
    use lndvc_def, only : vc_t, vc_cell_t, lndvc_class

    implicit none

    private
    public :: lndvc_aggregate_weights    ! area weights for the leaf reduction
    public :: lndvc_aggregate_cell       ! reduce leaf vc(:) -> coarse cell
    public :: lndvc_conservation_check   ! assert energy/water/carbon closure

contains

    subroutine lndvc_aggregate_weights(vc, wt)
        ! Return normalized area weights for the leaf reduction (Sum wt = 1).

        implicit none

        type(vc_t), intent(in)  :: vc(:)
        real(wp),   intent(out) :: wt(:)

        integer :: k

        do k = 1, size(vc)
            wt(k) = vc(k)%desc%w
        end do
        if (sum(wt) > 0._wp) wt = wt / sum(wt)

        return

    end subroutine lndvc_aggregate_weights

    subroutine lndvc_aggregate_cell(cell, vc, wt)
        ! Area-weighted reduction of the leaf ensemble to the coarse-cell mean.
        ! Class fractions from Sum of weights by class; coupling currency
        ! (t_skin, albedo, fluxes, runoff, evap) as area-weighted means; SMB as
        ! a conserving sum of ice-vc mass budgets.

        implicit none

        type(vc_cell_t), intent(out) :: cell
        type(vc_t),      intent(in)  :: vc(:)
        real(wp),        intent(in)  :: wt(:)

        ! TODO: per-field reduction rule (mean/sum/flux-conserving).

        return

    end subroutine lndvc_aggregate_cell

    subroutine lndvc_conservation_check(lnd)
        ! Assert energy/water/carbon closure across the aggregation boundary.

        implicit none

        type(lndvc_class), intent(in) :: lnd

        ! TODO: compare summed leaf budgets against coarse-cell aggregate.

        return

    end subroutine lndvc_conservation_check

end module lndvc_aggregate
