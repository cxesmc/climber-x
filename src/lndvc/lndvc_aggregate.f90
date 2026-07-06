module lndvc_aggregate
    ! Aggregation / conservation service for the virtual-cell framework.
    !
    ! Reduces the leaf virtual-cell ensemble back to coarse-cell means for
    ! coupling with CLIMBER, using per-variable reduction rules
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
    use lndvc_def

    implicit none

    private
    public :: lndvc_aggregate_weights    ! area weights for the leaf reduction
    public :: lndvc_aggregate_veg        ! reduce veg leaves -> coarse cell
    public :: lndvc_aggregate_smb        ! reduce smb/ice leaves -> coarse cell
    public :: lndvc_conservation_check   ! assert energy/water/carbon closure

contains

    subroutine lndvc_aggregate_weights(vc, wt)
        ! Return normalized area weights for the leaf reduction (Sum wt = 1).

        implicit none

        type(lndvc_vc_class), intent(in)  :: vc(:)
        real(wp),             intent(out) :: wt(:)

        wt = 0._wp
        ! TODO: wt(k) = area weight of leaf k relative to the coarse cell.
        if (sum(wt) > 0._wp) wt = wt / sum(wt)

        return

    end subroutine lndvc_aggregate_weights

    subroutine lndvc_aggregate_veg(veg, veg_vc, wt)
        ! Area-weighted reduction of vegetation/land leaf state to the cell mean.

        implicit none

        type(lndvc_veg_class), intent(out) :: veg
        type(lndvc_veg_class), intent(in)  :: veg_vc(:)
        real(wp),              intent(in)  :: wt(:)

        ! TODO: per-field reduction rule (mean/sum/flux-conserving).

        return

    end subroutine lndvc_aggregate_veg

    subroutine lndvc_aggregate_smb(smb, smb_vc, wt)
        ! Area-weighted reduction of smb/ice leaf state to the cell mean.

        implicit none

        type(lndvc_smb_class), intent(out) :: smb
        type(lndvc_smb_class), intent(in)  :: smb_vc(:)
        real(wp),              intent(in)  :: wt(:)

        ! TODO: per-field reduction rule; smb is a mass-budget diagnostic.

        return

    end subroutine lndvc_aggregate_smb

    subroutine lndvc_conservation_check(vc)
        ! Assert energy/water/carbon closure across the aggregation boundary.

        implicit none

        type(lndvc_class), intent(in) :: vc

        ! TODO: compare summed leaf budgets against coarse-cell aggregate.

        return

    end subroutine lndvc_conservation_check

end module lndvc_aggregate
