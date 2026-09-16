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
    use lndvc_grid, only : nsurf, flag_veg

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

        integer  :: k, m
        real(wp) :: w, wsum, f_veg_vc, fw
        real(wp) :: tskin_vc, alb_vc, sh_vc, lh_vc, g_vc, ev_vc, et_vc

        cell%f_veg  = 0._wp
        cell%f_lake = 0._wp
        cell%f_ice  = 0._wp
        cell%f_peat = 0._wp

        do k = 1, size(vc)
            select case(vc(k)%desc%class)
                case(1); cell%f_veg  = cell%f_veg  + vc(k)%desc%w
                case(2); cell%f_lake = cell%f_lake + vc(k)%desc%w
                case(3); cell%f_ice  = cell%f_ice  + vc(k)%desc%w
            end select
            ! peatland fraction (cell fraction, on the land vc's carbon block)
            if (vc(k)%desc%class == 1 .and. allocated(vc(k)%carb)) &
                cell%f_peat = cell%f_peat + vc(k)%carb%f_peat
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
            else if (vc(k)%desc%class == 2 .and. allocated(vc(k)%lake)) then
                w = vc(k)%desc%w
                ! extensive lake surface water balance (P + M - E), conserving sum
                cell%runoff = cell%runoff + w * vc(k)%lake%lake_water_tendency
                ! intensive surface fluxes (accumulate weighted; normalized below)
                cell%t_skin = cell%t_skin + w * vc(k)%flx%t_skin(1)
                cell%albedo = cell%albedo + w * vc(k)%flx%albedo(1)
                cell%flx_sh = cell%flx_sh + w * vc(k)%flx%flx_sh(1)
                cell%flx_lh = cell%flx_lh + w * vc(k)%flx%flx_lh(1)
                cell%flx_g  = cell%flx_g  + w * vc(k)%flx%flx_g(1)
                cell%evap   = cell%evap   + w * vc(k)%flx%evap_surface(1)
                cell%et     = cell%et     + w * vc(k)%flx%et(1)
                wsum = wsum + w
            else if (vc(k)%desc%class == 1 .and. allocated(vc(k)%veg)) then
                ! land vc: collapse the veg sub-tiles (frac-weighted mean over
                ! the bare+PFT tiles) to a vc-level intensive value, then
                ! area-weight by the vc weight. Land runoff (surface + drainage,
                ! C.2) is an extensive flux and is added as a conserving sum.
                f_veg_vc = sum(vc(k)%flx%frac_surf, mask=flag_veg.eq.1)
                if (f_veg_vc > 0._wp) then
                    tskin_vc = 0._wp; alb_vc = 0._wp; sh_vc = 0._wp; lh_vc = 0._wp
                    g_vc = 0._wp; ev_vc = 0._wp; et_vc = 0._wp
                    do m = 1, nsurf
                        if (flag_veg(m).eq.1) then
                            fw = vc(k)%flx%frac_surf(m)/f_veg_vc
                            tskin_vc = tskin_vc + vc(k)%flx%t_skin(m)*fw
                            alb_vc   = alb_vc   + vc(k)%flx%albedo(m)*fw
                            sh_vc    = sh_vc    + vc(k)%flx%flx_sh(m)*fw
                            lh_vc    = lh_vc    + vc(k)%flx%flx_lh(m)*fw
                            g_vc     = g_vc     + vc(k)%flx%flx_g(m)*fw
                            ev_vc    = ev_vc    + vc(k)%flx%evap_surface(m)*fw
                            et_vc    = et_vc    + vc(k)%flx%et(m)*fw
                        end if
                    end do
                    w = vc(k)%desc%w
                    cell%t_skin = cell%t_skin + w * tskin_vc
                    cell%albedo = cell%albedo + w * alb_vc
                    cell%flx_sh = cell%flx_sh + w * sh_vc
                    cell%flx_lh = cell%flx_lh + w * lh_vc
                    cell%flx_g  = cell%flx_g  + w * g_vc
                    cell%evap   = cell%evap   + w * ev_vc
                    cell%et     = cell%et     + w * et_vc
                    ! extensive land runoff (surface + subsurface drainage)
                    cell%runoff = cell%runoff + w * vc(k)%soil%runoff(1)
                    wsum = wsum + w
                end if
            end if
        end do

        if (wsum > 0._wp) then
            cell%t_skin = cell%t_skin / wsum
            cell%albedo = cell%albedo / wsum
            cell%flx_sh = cell%flx_sh / wsum
            cell%flx_lh = cell%flx_lh / wsum
            cell%flx_g  = cell%flx_g  / wsum
            cell%evap   = cell%evap   / wsum
            cell%et     = cell%et     / wsum
        end if

        ! --- rich diagnostic aggregation (mirrors reference lnd_surf / lnd_ts coverage) ----
        ! Every ice/lake vc holds an identical forcing broadcast from its cell,
        ! so we take it from the first vc we find with a live block. Snowpack,
        ! water balance, vegetation, soil carbon are aggregated area-weighted
        ! over the vcs that contribute (ice-vc snow + veg-vc veg, etc.).
        cell%t2m           = 0._wp
        cell%q2m           = 0._wp
        cell%lwnet         = 0._wp
        cell%swnet         = 0._wp
        cell%rain          = 0._wp
        cell%snow_flx      = 0._wp
        cell%h_snow        = 0._wp
        cell%w_snow        = 0._wp
        cell%w_snow_max    = 0._wp
        cell%snowmelt      = 0._wp
        cell%icemelt       = 0._wp
        cell%f_snow        = 0._wp
        cell%transpiration = 0._wp
        cell%evap_surface  = 0._wp
        cell%evap_can      = 0._wp
        cell%runoff_sur    = 0._wp
        cell%runoff_gw     = 0._wp
        cell%drainage      = 0._wp
        cell%w_table       = 0._wp
        cell%w_table_peat  = 0._wp
        cell%f_wet         = 0._wp
        cell%f_wetland     = 0._wp
        cell%f_peat_pot    = 0._wp
        cell%lai           = 0._wp
        cell%sai           = 0._wp
        cell%gpp           = 0._wp
        cell%npp           = 0._wp
        cell%veg_c         = 0._wp
        cell%gdd5          = 0._wp
        cell%t2m_min_mon   = 0._wp
        cell%t_soil_top    = 0._wp
        cell%theta_w_top   = 0._wp
        cell%alt           = 0._wp
        cell%soil_c        = 0._wp
        cell%peat_c        = 0._wp
        cell%soil_resp     = 0._wp

        ! forcing scalars: one representative vc suffices (broadcast is uniform)
        do k = 1, size(vc)
            if (vc(k)%desc%class /= 0) then
                cell%t2m      = vc(k)%forc%t2m
                cell%q2m      = vc(k)%forc%q2m
                cell%swnet    = vc(k)%forc%swnet
                cell%rain     = vc(k)%forc%rain
                cell%snow_flx = vc(k)%forc%snow
                cell%lwnet    = vc(k)%forc%lwdown ! downwelling; net lw computed in ebal per vc
                exit
            end if
        end do

        ! snowpack + water/veg/soil aggregation over active vcs
        do k = 1, size(vc)
            w = vc(k)%desc%w
            if (vc(k)%desc%class == 3 .and. allocated(vc(k)%snow)) then
                ! ice-vc snow contribution
                cell%h_snow     = cell%h_snow     + w * vc(k)%snow%h_snow
                cell%w_snow     = cell%w_snow     + w * vc(k)%snow%w_snow
                cell%w_snow_max = cell%w_snow_max + w * vc(k)%snow%w_snow_max
                cell%snowmelt   = cell%snowmelt   + w * vc(k)%snow%snowmelt
                cell%icemelt    = cell%icemelt    + w * vc(k)%snow%icemelt
                cell%f_snow     = cell%f_snow     + w * vc(k)%snow%f_snow
            else if (vc(k)%desc%class == 2 .and. allocated(vc(k)%snow)) then
                ! lake-vc snow contribution
                cell%h_snow     = cell%h_snow     + w * vc(k)%snow%h_snow
                cell%w_snow     = cell%w_snow     + w * vc(k)%snow%w_snow
                cell%snowmelt   = cell%snowmelt   + w * vc(k)%snow%snowmelt
            else if (vc(k)%desc%class == 1) then
                ! land-vc: rich aggregation
                if (allocated(vc(k)%snow)) then
                    cell%h_snow   = cell%h_snow   + w * vc(k)%snow%h_snow
                    cell%w_snow   = cell%w_snow   + w * vc(k)%snow%w_snow
                    cell%snowmelt = cell%snowmelt + w * vc(k)%snow%snowmelt
                end if
                if (allocated(vc(k)%flx%transpiration)) then
                    cell%transpiration = cell%transpiration + w * sum(vc(k)%flx%transpiration)
                    cell%evap_surface  = cell%evap_surface  + w * sum(vc(k)%flx%evap_surface)
                    cell%evap_can      = cell%evap_can      + w * sum(vc(k)%flx%evap_can)
                end if
                if (allocated(vc(k)%soil)) then
                    cell%runoff_sur   = cell%runoff_sur + w * (sum(vc(k)%soil%runoff_sur) + vc(k)%soil%runoff_exc)
                    cell%runoff_gw    = cell%runoff_gw  + w * vc(k)%soil%runoff_gw
                    cell%drainage     = cell%drainage   + w * sum(vc(k)%soil%drainage)
                    cell%w_table      = cell%w_table    + w * vc(k)%soil%w_table
                    cell%w_table_peat = cell%w_table_peat + w * vc(k)%soil%w_table_peat
                    cell%f_wet        = cell%f_wet      + w * vc(k)%soil%f_wet
                    cell%f_wetland    = cell%f_wetland  + w * vc(k)%soil%f_wetland
                    cell%alt          = cell%alt        + w * vc(k)%soil%alt
                    if (allocated(vc(k)%soil%t_soil)) cell%t_soil_top = cell%t_soil_top + w * vc(k)%soil%t_soil(1)
                    if (allocated(vc(k)%soil%theta_w)) cell%theta_w_top = cell%theta_w_top + w * vc(k)%soil%theta_w(1)
                end if
                if (allocated(vc(k)%veg)) then
                    cell%lai         = cell%lai         + w * sum(vc(k)%veg%lai * vc(k)%veg%pft_frac)
                    cell%sai         = cell%sai         + w * sum(vc(k)%veg%sai * vc(k)%veg%pft_frac)
                    cell%gpp         = cell%gpp         + w * sum(vc(k)%veg%gpp * vc(k)%veg%pft_frac)
                    cell%npp         = cell%npp         + w * sum(vc(k)%veg%npp * vc(k)%veg%pft_frac)
                    cell%veg_c       = cell%veg_c       + w * sum(vc(k)%veg%veg_c * vc(k)%veg%pft_frac)
                    cell%gdd5        = cell%gdd5        + w * vc(k)%veg%gdd5
                    cell%t2m_min_mon = cell%t2m_min_mon + w * vc(k)%veg%t2m_min_mon
                end if
                if (allocated(vc(k)%carb)) then
                    cell%f_peat_pot  = cell%f_peat_pot + vc(k)%carb%f_peat_pot
                    cell%soil_c      = cell%soil_c     + w * sum(vc(k)%carb%soil_c_tot)
                    cell%soil_resp   = cell%soil_resp  + w * sum(vc(k)%carb%soil_resp)
                    cell%peat_c      = cell%peat_c     + vc(k)%carb%f_peat &
                        * (vc(k)%carb%litter_c_peat + vc(k)%carb%acro_c + sum(vc(k)%carb%cato_c))
                end if
            end if
        end do

        ! --- emission aggregation (port W.1) ----------------------------------
        ! per-cell CH4/N2O emission flux [kgC/m2 cell/s], mirroring the reference
        ! lnd_update_wrapper weighting. Land vc: wetland CH4 over the non-peat
        ! wetland fraction (f_wetland*w - f_peat) + peatland CH4 over f_peat;
        ! N2O over the veg fraction. Lake vc: lake CH4 over its area weight.
        ! Shelf CH4 is out of scope (ice/shelf carbon not ported). The per-vc
        ! sum is exact at n_vc=1; a multi-band split of f_peat is a Phase-3 item.
        cell%ch4_emis = 0._wp
        cell%n2o_emis = 0._wp
        do k = 1, size(vc)
            w = vc(k)%desc%w
            if (vc(k)%desc%class == 1 .and. allocated(vc(k)%carb)) then
                cell%ch4_emis = cell%ch4_emis &
                    + vc(k)%carb%ch4_emis_wetland * max(0._wp, vc(k)%soil%f_wetland*w - vc(k)%carb%f_peat) &
                    + vc(k)%carb%ch4_emis_peat * vc(k)%carb%f_peat
                cell%n2o_emis = cell%n2o_emis + vc(k)%carb%n2o_emis * w
            else if (vc(k)%desc%class == 2 .and. allocated(vc(k)%lake)) then
                cell%ch4_emis = cell%ch4_emis + vc(k)%lake%ch4_emis_lake * w
            end if
        end do

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
