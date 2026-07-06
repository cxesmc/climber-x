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
    use semi_m,          only : semi
    use smb_par_m,       only : p0, h_atm, gamma, prc_par, surf_par

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
        ! Single-column ice-surface update via the ported SEMI driver (pattern A):
        ! the vc's per-class blocks are unpacked into SEMI's scalar signature.
        ! Genuine carry-over state binds to block fields (snow/ice/flx); pure
        ! per-call diagnostics are call-local (kept lean, see design sec 12).
        ! SEMI does its own elevation downscaling from forc%z_sur_i to desc%z.

        implicit none

        type(vc_t), intent(inout) :: vc

        integer  :: ii, jj
        real(wp) :: dz, dT
        real(wp) :: f_ele, pressure

        ! downscaled/diagnostic scalars SEMI writes but we do not persist
        real(wp) :: tam, t2m, q2m, qsat, dqsatdT, r_a
        real(wp) :: u700, v700, wind, snow_rate, rain_rate, prc, f_wind
        real(wp) :: alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif
        real(wp) :: alb_bg, cld, swnet, swnet_min, swdown, lwdown
        real(wp) :: dflxg_dT, flx_melt, flx_lwu
        real(wp) :: num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw
        real(wp) :: f_sh, f_e, f_lh, f_lw

        ii = 0; jj = 0   ! SEMI uses these only for debug prints

        ! --- geometry-derived forcing (framework side) ------------------------
        ! barometric surface pressure at the vc elevation
        pressure = p0 * exp(-vc%desc%z / h_atm)
        vc%forc%pressure = pressure
        ! precipitation elevation-correction factor (Clausius-Clapeyron), cf. topo_factors
        if (prc_par%l_elevation_corr .and. vc%desc%z >= prc_par%z_sur_crit_fele) then
            dz = vc%desc%z - vc%forc%z_sur_i
            dT = -gamma*dz
            f_ele = exp(prc_par%dP_dT*dT)
        else
            f_ele = 1._wp
        end if
        if (vc%desc%z >= prc_par%z_sur_high_fele) f_ele = 0.1_wp*f_ele
        vc%desc%f_ele = f_ele
        ! pure ice vc: fully ice-covered patch, firn background albedo
        vc%forc%f_ice   = 1._wp
        vc%forc%alb_ice = surf_par%alb_firn

        ! --- surface energy + mass balance (ported SEMI, pattern A) ------------
        call semi(ii, jj, vc%forc%f_ice, vc%forc%alb_ice, &
            vc%desc%z, vc%forc%z_sur_i, vc%desc%z_sur_std, &
            vc%desc%dz_dx, vc%desc%dz_dy, vc%desc%dz_sur, vc%desc%f_ele, &
            vc%forc%tam_i, vc%forc%t2m_bias_i, vc%forc%dTvar, vc%forc%gam_i, vc%forc%tstd_i, vc%forc%ram_i, &
            pressure, vc%forc%u700_i, vc%forc%v700_i, vc%forc%wind_i, vc%forc%prc_i, vc%forc%prc_bias_i, &
            vc%forc%alb_vis_dir_i, vc%forc%alb_nir_dir_i, vc%forc%alb_vis_dif_i, vc%forc%alb_nir_dif_i, &
            vc%forc%swd_sur_vis_dir_i, vc%forc%swd_sur_nir_dir_i, vc%forc%swd_sur_vis_dif_i, vc%forc%swd_sur_nir_dif_i, &
            vc%forc%dswd_dalb_vis_dir_i, vc%forc%dswd_dalb_nir_dir_i, vc%forc%dswd_dalb_vis_dif_i, vc%forc%dswd_dalb_nir_dif_i, &
            vc%forc%dswd_dz_nir_dir_i, vc%forc%dswd_dz_nir_dif_i, vc%forc%dust_i, vc%forc%coszm_i, &
            vc%forc%swd_toa_i, vc%forc%swd_toa_min_i, vc%forc%cld_i, vc%forc%lwdown_i, vc%forc%gam_lw_i, &
            tam, t2m, vc%flx%t_skin(1), vc%flx%t_skin_old(1), vc%flx%t_skin_amp(1), &
            vc%ice%t_prof, vc%ice%t_prof_old, &
            q2m, qsat, dqsatdT, r_a, &
            u700, v700, wind, snow_rate, rain_rate, prc, f_wind, &
            vc%snow%mask_snow, vc%snow%f_snow, vc%snow%h_snow, vc%snow%w_snow, vc%snow%w_snow_old, vc%snow%w_snow_max, &
            vc%snow%snow_grain, vc%snow%dust_con, &
            alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif, &
            vc%snow%alb_snow_vis_dir, vc%snow%alb_snow_nir_dir, vc%snow%alb_snow_vis_dif, vc%snow%alb_snow_nir_dif, &
            vc%snow%dt_snowfree, alb_bg, vc%flx%albedo(1), cld, swnet, swnet_min, swdown, lwdown, &
            vc%flx%flx_g(1), dflxg_dT, flx_melt, vc%flx%flx_sh(1), flx_lwu, vc%flx%flx_lh(1), vc%flx%evap_surface(1), &
            num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw, &
            f_sh, f_e, f_lh, f_lw, &
            vc%snow%snowmelt, vc%snow%icemelt, vc%snow%refreezing, vc%snow%refreezing_sum, vc%ice%f_rfz_to_snow)

        ! --- SMB mass-budget diagnostic (kg/m2/s) -----------------------------
        vc%ice%melt   = vc%snow%snowmelt + vc%snow%icemelt
        vc%ice%runoff = vc%snow%snowmelt + vc%snow%icemelt + rain_rate - vc%snow%refreezing
        vc%ice%smb    = snow_rate - vc%flx%evap_surface(1) - vc%snow%snowmelt - vc%snow%icemelt + vc%snow%refreezing

        return
    end subroutine lndvc_update_ice

end module lndvc_model
