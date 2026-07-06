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
        !
        ! Phase 1: class fractions only (identity roundtrip of the geometry).

        implicit none

        type(vc_cell_t), intent(out) :: cell
        type(vc_t),      intent(in)  :: vc(:)
        real(wp),        intent(in)  :: wt(:)

        integer  :: k
        real(wp) :: w, wsum

        cell%f_veg  = 0._wp
        cell%f_lake = 0._wp
        cell%f_ice  = 0._wp

        do k = 1, size(vc)
            select case(vc(k)%desc%class)
                case(1); cell%f_veg  = cell%f_veg  + vc(k)%desc%w
                case(2); cell%f_lake = cell%f_lake + vc(k)%desc%w
                case(3); cell%f_ice  = cell%f_ice  + vc(k)%desc%w
            end select
        end do

        cell%f_land    = cell%f_veg + cell%f_lake + cell%f_ice
        cell%mask_lnd  = merge(1, 0, cell%f_land > 0._wp)

        ! Reduce the coupling currency from the leaf virtual cells whose surface
        ! physics is live. Only the ice/SMB path is ported so far, so surface
        ! fluxes are the ice-surface area-weighted mean and the mass budget is
        ! the conserving cell-area-weighted sum. Extends to land/lake vcs as
        ! those single-column paths come online.
        cell%t_skin = 0._wp
        cell%albedo = 0._wp
        cell%flx_sh = 0._wp
        cell%flx_lh = 0._wp
        cell%flx_g  = 0._wp
        cell%evap   = 0._wp
        cell%runoff = 0._wp
        cell%smb    = 0._wp
        cell%melt   = 0._wp
        cell%et     = 0._wp
        wsum = 0._wp

        do k = 1, size(vc)
            if (vc(k)%desc%class == 3 .and. allocated(vc(k)%ice)) then
                w = vc(k)%desc%w
                ! extensive mass budget (per unit cell area, conserving sum)
                cell%smb    = cell%smb    + w * vc(k)%ice%smb
                cell%melt   = cell%melt   + w * vc(k)%ice%melt
                cell%runoff = cell%runoff + w * vc(k)%ice%runoff
                ! intensive surface fluxes (accumulate weighted; normalized below)
                cell%t_skin = cell%t_skin + w * vc(k)%flx%t_skin(1)
                cell%albedo = cell%albedo + w * vc(k)%flx%albedo(1)
                cell%flx_sh = cell%flx_sh + w * vc(k)%flx%flx_sh(1)
                cell%flx_lh = cell%flx_lh + w * vc(k)%flx%flx_lh(1)
                cell%flx_g  = cell%flx_g  + w * vc(k)%flx%flx_g(1)
                cell%evap   = cell%evap   + w * vc(k)%flx%evap_surface(1)
                wsum = wsum + w
            end if
        end do

        if (wsum > 0._wp) then
            cell%t_skin = cell%t_skin / wsum
            cell%albedo = cell%albedo / wsum
            cell%flx_sh = cell%flx_sh / wsum
            cell%flx_lh = cell%flx_lh / wsum
            cell%flx_g  = cell%flx_g  / wsum
            cell%evap   = cell%evap   / wsum
        end if

        return

    end subroutine lndvc_aggregate_cell

    subroutine lndvc_conservation_check(lnd)
        ! Sanity guard on the aggregation boundary. Phase 1: aggregated class
        ! fractions must be non-negative and f_land must lie in [0,1].
        ! (Phase 2 extends this to energy/water/carbon closure.)

        implicit none

        type(lndvc_class), intent(in) :: lnd

        integer :: n, i, j
        real(wp), parameter :: eps = 1.e-6_wp

        do n = 1, lnd%ncells
            i = lnd%ij_1d(1,n)
            j = lnd%ij_1d(2,n)
            associate(c => lnd%cell(i,j))
            if (c%f_veg < -eps .or. c%f_lake < -eps .or. c%f_ice < -eps .or. &
                c%f_land < -eps .or. c%f_land > 1._wp+eps) then
                write(*,'(a,2i5,4f10.5)') 'lndvc warning: fraction out of range at i,j:', &
                    i, j, c%f_veg, c%f_lake, c%f_ice, c%f_land
            end if
            end associate
        end do

        return

    end subroutine lndvc_conservation_check

end module lndvc_aggregate
