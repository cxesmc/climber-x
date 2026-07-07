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
    ! ported single-column lake physics
    use lndvc_surface_par_lnd, only : resist_aer_lake, snow_albedo_lake, surface_albedo_lake, resist_sur_lake
    use lndvc_lake_par_mod,    only : lake_par_thermal, sublake_par_thermal
    use lndvc_ebal_lake_mod,   only : ebal_lake, update_tskin_lake
    use lndvc_lake_temp_mod,   only : lake_temp
    use lndvc_sublake_temp_mod, only : sublake_temp
    use lndvc_lake_carbon_par_mod, only : lake_carbon_par
    use lndvc_lake_carbon_mod, only : lake_carbon
    use lndvc_hydrology_mod,   only : surface_hydrology_lake
    ! ported vegetated-land physics (port C)
    use lndvc_surface_par_lnd, only : resist_aer_veg, resist_sur_veg, surface_albedo_veg
    use lndvc_veg_par_mod,     only : phenology, root_frac_update, litter_in_frac_update, dynveg_par
    use lndvc_photosynthesis_mod, only : photosynthesis
    use lndvc_hydrology_mod,   only : canopy_water, surface_hydrology_veg
    use lndvc_soil_par_mod,    only : soil_par_thermal, soil_par_hydro, soil_par_update
    use lndvc_soil_hydro_mod,  only : soil_hydro
    use lndvc_water_deficit_mod, only : calculate_pet, calculate_cwd
    use lndvc_dyn_veg_mod,     only : dyn_veg
    use lndvc_soil_carbon_par_mod, only : soil_carbon_par
    use lndvc_soil_carbon_mod, only : soil_carbon
    use lndvc_peat_carbon_mod, only : peat_carbon
    use lndvc_n2o_emis_mod,    only : n2o_emission
    use lndvc_dust_emis_mod,   only : dust_emission
    use lndvc_carbon_flux_atm_lnd_mod, only : carbon_flux_atm_lnd
    use lndvc_weathering_mod,  only : weathering_gemco2, weathering_uhh
    use lndvc_carbon_export_mod, only : carbon_export
    use lndvc_carbon_inventory_mod, only : carbon_inventory
    use lndvc_ebal_veg_mod,    only : ebal_veg, update_tskin_veg
    use lndvc_soil_temp_mod,   only : soil_temp
    use lndvc_init_cell_mod,   only : lndvc_init_cell_veg
    use lnd_params,            only : lnd_surf_par => surf_par
    use lnd_params,            only : soil_par, hydro_par, veg_par, peat_par
    use lnd_params,            only : time_call_veg, time_call_carb, time_call_carb_p
    use lnd_params,            only : dt                       ! land timestep (global emission accumulation)
    use lnd_params,            only : i_weathering, l_river_export
    use lnd_params,            only : soilc_par                 ! l_burial (carbon burial branch)
    use climber_grid,          only : area                     ! coarse-cell area [m2]
    use timer,                 only : time_soy_lnd, time_eoy_lnd, time_eom_lnd
    use timer,                 only : year                     ! year index (running averages)
    use constants,             only : T0
    use wiso_params,           only : l_wiso, nwiso, i_o18, Rstd

    implicit none

    private
    public :: lndvc_init
    public :: lndvc_update
    public :: lndvc_init_land
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
        ! cell-level carbon-flux inputs assembled from the vcs (port W.2)
        real(wp) :: sresp(ncarb), sresp13(ncarb), sresp14(ncarb)
        real(wp) :: npp_real, npp13_real, npp14_real
        real(wp) :: fire_c_flux, fire_c13_flux, fire_c14_flux
        ! annual global carbon accounting (port W.5)
        real(wp) :: landc, landc13, landc14, burc, burc13, burc14
        real(wp) :: glandc, glandc13, glandc14, gburc, gburc13, gburc14
        real(wp) :: dburc, dburc13, dburc14, atot, wcarb, wsil

        allocate(wt(lnd%n_vc))

        ! reset the annual global accumulators at start-of-year (mirror
        ! lnd_update_wrapper); per-cell emission fluxes + the monthly
        ! atmosphere-land carbon flux accumulate below
        if (time_soy_lnd) then
            lnd%glob%ch4_emis = 0._wp
            lnd%glob%n2o_emis = 0._wp
            lnd%glob%Cflx_atm_lnd   = 0._wp
            lnd%glob%C13flx_atm_lnd = 0._wp
            lnd%glob%C14flx_atm_lnd = 0._wp
        end if

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
            ! accumulate the annual global CH4/N2O emission (kgC), mirroring the
            ! reference lnd_update_wrapper: per-cell flux [kgC/m2/s] * dt * area
            lnd%glob%ch4_emis = lnd%glob%ch4_emis + lnd%cell(i,j)%ch4_emis * dt * area(i,j)
            lnd%glob%n2o_emis = lnd%glob%n2o_emis + lnd%cell(i,j)%n2o_emis * dt * area(i,j)

            ! net atmosphere-land carbon flux (port W.2), monthly. Assemble the
            ! cell carbon-flux inputs from the vcs (n_vc=1: direct per-class
            ! values; a multi-band frac-weighted reduction is a Phase-3 item),
            ! then accumulate to the cell (kgC/s) + annual global (kgC/mon sum).
            ! ic_ice / ic_shelf respiration is zero (ice/shelf carbon out of
            ! scope) — the f_ice_grd / f_shelf weights below multiply zero.
            if (time_eom_lnd) then
                sresp(:) = 0._wp; sresp13(:) = 0._wp; sresp14(:) = 0._wp
                npp_real = 0._wp; npp13_real = 0._wp; npp14_real = 0._wp
                fire_c_flux = 0._wp; fire_c13_flux = 0._wp; fire_c14_flux = 0._wp
                do k = 1, lnd%n_vc
                    associate(vc => lnd%vc(i,j,k))
                    if (vc%desc%class == 1 .and. allocated(vc%carb)) then
                        sresp(ic_min)    = vc%carb%soil_resp(ic_min)
                        sresp(ic_peat)   = vc%carb%soil_resp(ic_peat)
                        sresp13(ic_min)  = vc%carb%soil_resp13(ic_min)
                        sresp13(ic_peat) = vc%carb%soil_resp13(ic_peat)
                        sresp14(ic_min)  = vc%carb%soil_resp14(ic_min)
                        sresp14(ic_peat) = vc%carb%soil_resp14(ic_peat)
                        npp_real   = vc%veg%npp_real;   npp13_real = vc%veg%npp13_real;   npp14_real = vc%veg%npp14_real
                        fire_c_flux = vc%veg%fire_c_flux; fire_c13_flux = vc%veg%fire_c13_flux; fire_c14_flux = vc%veg%fire_c14_flux
                    else if (vc%desc%class == 2 .and. allocated(vc%lake)) then
                        sresp(ic_lake)   = vc%lake%soil_resp_lake
                        sresp13(ic_lake) = vc%lake%soil_resp13_lake
                        sresp14(ic_lake) = vc%lake%soil_resp14_lake
                    end if
                    end associate
                end do
                call carbon_flux_atm_lnd(lnd%cell(i,j)%f_veg, lnd%cell(i,j)%f_peat, &
                    lnd%cell(i,j)%f_ice_grd, lnd%cell(i,j)%f_shelf, lnd%cell(i,j)%f_lake, area(i,j), &
                    npp_real, npp13_real, npp14_real, sresp, sresp13, sresp14, &
                    fire_c_flux, fire_c13_flux, fire_c14_flux, &
                    lnd%cell(i,j)%Cflx_atm_lnd, lnd%cell(i,j)%C13flx_atm_lnd, lnd%cell(i,j)%C14flx_atm_lnd, &
                    lnd%glob%Cflx_atm_lnd, lnd%glob%C13flx_atm_lnd, lnd%glob%C14flx_atm_lnd)
            end if
        end do

        ! ===== annual global carbon accounting (port W.5), end-of-year ========
        ! total land carbon + burial (carbon_inventory), the burial flux, and the
        ! running-average diagnostics — mirror lnd_update_wrapper's post-loop.
        ! ice/shelf sediment carbon is out of scope (zero pools); the total land
        ! carbon therefore omits it (documented gap, option A).
        if (time_eoy_lnd) then
            glandc = 0._wp; glandc13 = 0._wp; glandc14 = 0._wp
            gburc  = 0._wp; gburc13  = 0._wp; gburc14  = 0._wp
            do n = 1, lnd%ncells
                i = lnd%ij_1d(1,n); j = lnd%ij_1d(2,n)
                call lndvc_carbon_inventory_cell(lnd%vc(i,j,:), lnd%cell(i,j), &
                    landc, landc13, landc14, burc, burc13, burc14)
                glandc = glandc + landc*area(i,j); glandc13 = glandc13 + landc13*area(i,j); glandc14 = glandc14 + landc14*area(i,j)
                gburc  = gburc  + burc*area(i,j);  gburc13  = gburc13  + burc13*area(i,j);  gburc14  = gburc14  + burc14*area(i,j)
            end do

            ! burial flux: with burial the buried carbon leaves the system; without
            ! it, the carbon reaching the burial layer returns to the atmosphere
            dburc = gburc - lnd%glob%burc; dburc13 = gburc13 - lnd%glob%burc13; dburc14 = gburc14 - lnd%glob%burc14
            if (soilc_par%l_burial) then
                lnd%glob%Cflx_burial   = dburc
                lnd%glob%C13flx_burial = dburc13
                lnd%glob%C14flx_burial = dburc14
            else
                lnd%glob%Cflx_burial = 0._wp; lnd%glob%C13flx_burial = 0._wp; lnd%glob%C14flx_burial = 0._wp
                lnd%glob%Cflx_atm_lnd   = lnd%glob%Cflx_atm_lnd   - dburc
                lnd%glob%C13flx_atm_lnd = lnd%glob%C13flx_atm_lnd - dburc13
                lnd%glob%C14flx_atm_lnd = lnd%glob%C14flx_atm_lnd - dburc14
                atot = sum(area)
                do n = 1, lnd%ncells
                    i = lnd%ij_1d(1,n); j = lnd%ij_1d(2,n)
                    lnd%cell(i,j)%Cflx_atm_lnd    = lnd%cell(i,j)%Cflx_atm_lnd    - dburc  /dt/atot
                    lnd%cell(i,j)%C13flx_atm_lnd  = lnd%cell(i,j)%C13flx_atm_lnd  - dburc13/dt/atot
                    lnd%cell(i,j)%C14flx_atm_lnd  = lnd%cell(i,j)%C14flx_atm_lnd  - dburc14/dt/atot
                end do
            end if

            ! save totals for next year
            lnd%glob%landc = glandc; lnd%glob%landc13 = glandc13; lnd%glob%landc14 = glandc14
            lnd%glob%burc  = gburc;  lnd%glob%burc13  = gburc13;  lnd%glob%burc14  = gburc14

            ! running-average land-atmosphere carbon flux (from year 2)
            if (year > 1) then
                if (year == 2) lnd%glob%Cflx_avg = lnd%glob%Cflx_atm_lnd
                lnd%glob%Cflx_avg = 0.99_wp*lnd%glob%Cflx_avg + 0.01_wp*lnd%glob%Cflx_atm_lnd
            end if

            ! running-average global weathering fluxes (kgC/yr; from year 2)
            if (year > 1) then
                wcarb = 0._wp; wsil = 0._wp
                do n = 1, lnd%ncells
                    i = lnd%ij_1d(1,n); j = lnd%ij_1d(2,n)
                    do k = 1, lnd%n_vc
                        if (lnd%vc(i,j,k)%desc%class == 1 .and. allocated(lnd%vc(i,j,k)%carb) &
                            .and. lnd%vc(i,j,k)%forc%f_veg_cell > 0._wp) then
                            ! mol C/m2/yr * (cell veg area) * 12 g/mol * 1e-3 kg/g = kgC/yr
                            wcarb = wcarb + lnd%vc(i,j,k)%carb%weath_carb * area(i,j)*lnd%cell(i,j)%f_veg * 12._wp*1e-3_wp
                            wsil  = wsil  + lnd%vc(i,j,k)%carb%weath_sil  * area(i,j)*lnd%cell(i,j)%f_veg * 12._wp*1e-3_wp
                        end if
                    end do
                end do
                if (year == 2) lnd%glob%weath_carb_avg = wcarb
                lnd%glob%weath_carb_avg = 0.99_wp*lnd%glob%weath_carb_avg + 0.01_wp*wcarb
                if (year == 2) lnd%glob%weath_sil_avg = wsil
                lnd%glob%weath_sil_avg = 0.99_wp*lnd%glob%weath_sil_avg + 0.01_wp*wsil
            end if
        end if

        deallocate(wt)

        call lndvc_conservation_check(lnd)

        return

    end subroutine lndvc_update

    subroutine lndvc_carbon_inventory_cell(vc, cell, landc, landc13, landc14, burc, burc13, burc14)
        ! Assemble a coarse cell's carbon pools from its virtual cells (n_vc=1:
        ! direct per-class values) and evaluate the reference carbon_inventory.
        ! Vegetation + mineral + peat come from the land vc, lake sediment from
        ! the lake vc; ice/shelf sediment carbon is out of scope (zero pools, so
        ! the total land carbon omits it). A multi-band frac-weighted reduction
        ! is a Phase-3 item.

        implicit none

        type(vc_t),      intent(in)  :: vc(:)
        type(vc_cell_t), intent(in)  :: cell
        real(wp),        intent(out) :: landc, landc13, landc14, burc, burc13, burc14

        integer  :: k
        real(wp) :: fsurf(nsurf)
        real(wp) :: vgc(npft), vgc13(npft), vgc14(npft)
        real(wp) :: lit(nlc), fst(nlc), slw(nlc), lit13(nlc), fst13(nlc), slw13(nlc), lit14(nlc), fst14(nlc), slw14(nlc)
        real(wp) :: cato(nlc), cato13(nlc), cato14(nlc)
        real(wp) :: litp, acro, litp13, acro13, litp14, acro14
        real(wp) :: lkl(nlc), lkf(nlc), lks(nlc), lkl13(nlc), lkf13(nlc), lks13(nlc), lkl14(nlc), lkf14(nlc), lks14(nlc)
        real(wp) :: zc(nlc)   ! zero ice/shelf sediment pools (out of scope)

        fsurf = 0._wp; vgc = 0._wp; vgc13 = 0._wp; vgc14 = 0._wp
        lit = 0._wp; fst = 0._wp; slw = 0._wp; lit13 = 0._wp; fst13 = 0._wp; slw13 = 0._wp
        lit14 = 0._wp; fst14 = 0._wp; slw14 = 0._wp
        cato = 0._wp; cato13 = 0._wp; cato14 = 0._wp
        litp = 0._wp; acro = 0._wp; litp13 = 0._wp; acro13 = 0._wp; litp14 = 0._wp; acro14 = 0._wp
        lkl = 0._wp; lkf = 0._wp; lks = 0._wp; lkl13 = 0._wp; lkf13 = 0._wp; lks13 = 0._wp
        lkl14 = 0._wp; lkf14 = 0._wp; lks14 = 0._wp
        zc = 0._wp

        do k = 1, size(vc)
            if (vc(k)%desc%class == 1 .and. allocated(vc(k)%carb)) then
                fsurf = vc(k)%flx%frac_surf
                vgc = vc(k)%veg%veg_c; vgc13 = vc(k)%veg%veg_c13; vgc14 = vc(k)%veg%veg_c14
                lit = vc(k)%carb%litter_c; fst = vc(k)%carb%fast_c; slw = vc(k)%carb%slow_c
                lit13 = vc(k)%carb%litter_c13; fst13 = vc(k)%carb%fast_c13; slw13 = vc(k)%carb%slow_c13
                lit14 = vc(k)%carb%litter_c14; fst14 = vc(k)%carb%fast_c14; slw14 = vc(k)%carb%slow_c14
                cato = vc(k)%carb%cato_c; cato13 = vc(k)%carb%cato_c13; cato14 = vc(k)%carb%cato_c14
                litp = vc(k)%carb%litter_c_peat; acro = vc(k)%carb%acro_c
                litp13 = vc(k)%carb%litter_c13_peat; acro13 = vc(k)%carb%acro_c13
                litp14 = vc(k)%carb%litter_c14_peat; acro14 = vc(k)%carb%acro_c14
            else if (vc(k)%desc%class == 2 .and. allocated(vc(k)%lake)) then
                lkl = vc(k)%lake%litter_c_lake; lkf = vc(k)%lake%fast_c_lake; lks = vc(k)%lake%slow_c_lake
                lkl13 = vc(k)%lake%litter_c13_lake; lkf13 = vc(k)%lake%fast_c13_lake; lks13 = vc(k)%lake%slow_c13_lake
                lkl14 = vc(k)%lake%litter_c14_lake; lkf14 = vc(k)%lake%fast_c14_lake; lks14 = vc(k)%lake%slow_c14_lake
            end if
        end do

        call carbon_inventory(fsurf, cell%f_veg, cell%f_peat, cell%f_ice_grd, cell%f_shelf, cell%f_lake, &
            vgc, vgc13, vgc14, &
            lit, fst, slw, lit13, fst13, slw13, lit14, fst14, slw14, &
            litp, acro, cato, litp13, acro13, cato13, litp14, acro14, cato14, &
            zc, zc, zc, zc, zc, zc, zc, zc, zc, &   ! ice sediment (out of scope)
            zc, zc, zc, zc, zc, zc, zc, zc, zc, &   ! shelf sediment (out of scope)
            lkl, lkf, lks, lkl13, lkf13, lks13, lkl14, lkf14, lks14, &
            landc, landc13, landc14, burc, burc13, burc14)

        return

    end subroutine lndvc_carbon_inventory_cell

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

    subroutine lndvc_init_land(vc, m_theta_sat, m_k_sat, m_psi_sat, m_Bi, m_lambda_s, m_lambda_dry, &
                               cti_mean, cti_cdf, dyptop_k, dyptop_v, dyptop_xm, dyptop_fmax, &
                               c13_c12_atm, c14_c_atm)
        ! One-time physical init of a vegetated land vc (port C decision B):
        ! seed the soil parameters from the coarse-cell mineral texture (the
        ! reference lnd-init block, src/lnd/lnd_model.f90 ~1206) then run the
        ! ported init_cell_veg state init. frac_surf is set from pft_frac (bare
        ! is the remainder) since the multi-class surface_frac_up does not apply
        ! to a single-class veg vc. Called once per land vc from cmn_to_lndvc.

        implicit none

        type(vc_t), intent(inout) :: vc
        real(wp), intent(in) :: m_theta_sat, m_k_sat, m_psi_sat, m_Bi, m_lambda_s, m_lambda_dry
        real(wp), intent(in) :: cti_mean, cti_cdf(:)
        real(wp), intent(in) :: dyptop_k, dyptop_v, dyptop_xm, dyptop_fmax
        real(wp), intent(in) :: c13_c12_atm, c14_c_atm

        integer  :: k
        real(wp) :: sum_pft

        ! --- soil parameters from mineral texture -----------------------------
        do k = 1, nl
            vc%soil%theta_sat(k) = m_theta_sat
            vc%soil%psi_sat(k)   = m_psi_sat
            vc%soil%k_sat(k)     = m_k_sat
            vc%soil%k_exp(k)     = 2*nint(m_Bi)+3
            vc%soil%psi_exp(k)   = -nint(m_Bi)
            vc%soil%theta_field(k) = vc%soil%theta_sat(k) * (0.1_wp / (86400._wp * vc%soil%k_sat(k)))**(1._wp/vc%soil%k_exp(k))
            vc%soil%theta_wilt(k)  = vc%soil%theta_sat(k) * (hydro_par%p_psi_min / vc%soil%psi_sat(k))**(1._wp/vc%soil%psi_exp(k))
            vc%soil%lambda_s(k)   = m_lambda_s
            vc%soil%lambda_dry(k) = m_lambda_dry
        enddo
        if (soil_par%uniform_porosity) vc%soil%theta_sat(:) = soil_par%theta_sat_u
        if (soil_par%uniform_soil_par_therm) then
            vc%soil%lambda_s(:)   = soil_par%lambda_s_u
            vc%soil%lambda_dry(:) = soil_par%lambda_dry_u
        endif

        ! --- persist mineral-soil baseline for soil_par_update carbon feedback -
        vc%soil%mineral_theta_sat  = m_theta_sat
        vc%soil%mineral_k_sat      = m_k_sat
        vc%soil%mineral_psi_sat    = m_psi_sat
        vc%soil%mineral_Bi         = m_Bi
        vc%soil%mineral_lambda_s   = m_lambda_s
        vc%soil%mineral_lambda_dry = m_lambda_dry

        ! --- static wetland parameters (TOPMODEL cti + DYPTOP) ----------------
        vc%soil%cti_mean    = cti_mean
        vc%soil%cti_cdf(:)  = cti_cdf(:)
        vc%soil%dyptop_k    = dyptop_k
        vc%soil%dyptop_v    = dyptop_v
        vc%soil%dyptop_xm   = dyptop_xm
        vc%soil%dyptop_fmax = dyptop_fmax

        ! --- physical veg/soil state (ported init_cell_veg) -------------------
        call lndvc_init_cell_veg(c13_c12_atm, c14_c_atm, vc%desc%w, &
            vc%veg%pft_frac, vc%flx%t_skin, vc%soil%t_soil, vc%flx%w_can, vc%flx%s_can, &
            vc%snow%w_snow, vc%snow%w_snow_max, vc%snow%h_snow, vc%snow%mask_snow, &
            vc%soil%theta, vc%soil%theta_w, vc%soil%theta_i, vc%soil%w_w, vc%soil%w_i, vc%soil%theta_sat, &
            vc%soil%alt, vc%veg%gdd5, vc%veg%gdd, vc%veg%phen, vc%veg%phen_acc, vc%veg%lai_bal, &
            vc%veg%lai, vc%veg%sai, vc%soil%root_frac, vc%carb%litter_in_frac, &
            vc%veg%gamma_dist, vc%veg%gamma_fire, vc%veg%npp_ann, vc%veg%npp13_ann, vc%veg%npp14_ann, &
            vc%veg%leaf_c, vc%veg%root_c, vc%veg%stem_c, vc%veg%veg_h, vc%veg%veg_c, vc%veg%veg_c13, vc%veg%veg_c14, &
            vc%carb%f_peat, vc%carb%f_peat_pot, vc%soil%w_table_min, vc%soil%w_table_peat, vc%carb%dCpeat_dt, &
            vc%flx%z0m)

        ! --- sub-tile fractions within the vc (bare = remainder) --------------
        call set_frac_surf_veg(vc)

        return

    end subroutine lndvc_init_land

    subroutine set_frac_surf_veg(vc)
        ! Within-vc sub-tile fractions from the PFT fractions: each PFT tile takes
        ! its pft_frac directly, bare soil is the remainder (frac_surf sums to 1
        ! inside the veg vc). This is the vc-local analogue of the reference
        ! surface_frac_up, which instead builds a cell-global frac_surf folding in
        ! the ice/lake tiles that separate vcs own (decision 2). Used at init and
        ! after every pft_frac update (start of year + post dyn_veg).
        implicit none
        type(vc_t), intent(inout) :: vc
        integer  :: k
        real(wp) :: sum_pft
        sum_pft = 0._wp
        do k = 1, npft
            vc%flx%frac_surf(k) = vc%veg%pft_frac(k)
            sum_pft = sum_pft + vc%veg%pft_frac(k)
        enddo
        vc%flx%frac_surf(i_bare) = max(0._wp, 1._wp - sum_pft)
    end subroutine set_frac_surf_veg

    subroutine lndvc_update_land(vc)
        ! Single-column vegetated-land surface update (port C.1): pattern-A
        ! unpack of the vc blocks into the ported veg chain, mirroring the
        ! reference lnd veg path (phenology -> aerodynamic/surface resistance ->
        ! albedo -> photosynthesis -> canopy water -> soil thermal params ->
        ! surface energy balance -> soil temperature -> skin update). The vc
        ! holds a single scalar forcing; the pattern-A routines expect nsurf-tile
        ! forcing arrays, so the scalar is broadcast to local nsurf arrays here
        ! (identity baseline: one elevation, uniform forcing across sub-tiles).
        ! Soil hydrology (C.2), veg dynamics (C.3) and soil carbon (C.4) are
        ! separate stages; the water-isotope freeze/thaw redistribution after
        ! soil_temp is deferred (like the SMB/lake iso deferrals).

        implicit none

        type(vc_t), intent(inout) :: vc

        integer, parameter :: ii = 0, jj = 0   ! debug-print indices only
        integer  :: n
        real(wp) :: f_veg, t_skin_veg
        real(wp) :: flx_g_veg, dflxg_dT_veg, flx_melt_veg
        real(wp) :: soil_resp   ! combined mineral+peat soil respiration (n2o input)
        ! scalar forcing broadcast to the nsurf sub-tiles
        real(wp) :: tatm(nsurf), qatm(nsurf), t2m(nsurf), q2m(nsurf), wind(nsurf)
        real(wp) :: pressure(nsurf), swnet(nsurf), swnet_min(nsurf), lwdown(nsurf)
        real(wp) :: rain(nsurf), snow(nsurf)
        ! per-tile conservation diagnostics (vc keeps only scalar residuals)
        real(wp) :: encons1(nsurf), encons2(nsurf)
        ! potential evapotranspiration (per-tile local; only mcwd/cwd_mon persist)
        real(wp) :: pet(nsurf)

        ! --- broadcast the vc's scalar forcing to the sub-tiles ---------------
        tatm(:)      = vc%forc%tatm
        qatm(:)      = vc%forc%qatm
        t2m(:)       = vc%forc%t2m
        q2m(:)       = vc%forc%q2m
        wind(:)      = vc%forc%wind
        pressure(:)  = vc%forc%pressure
        swnet(:)     = vc%forc%swnet
        swnet_min(:) = vc%forc%swnet_min
        lwdown(:)    = vc%forc%lwdown
        rain(:)      = vc%forc%rain
        snow(:)      = vc%forc%snow

        ! ===== start-of-year vegetation update (port C.3) =====================
        ! sync sub-tile fractions to the PFT distribution (updated by last year's
        ! dyn_veg), then refresh soil parameters / root + litter profiles. Guarded
        ! on the cell veg fraction, mirroring the reference lnd start-of-year block.
        if (time_soy_lnd) then
            call set_frac_surf_veg(vc)
            if (vc%forc%f_veg_cell .gt. 0._wp) then
                if (.not.soil_par%constant_porosity .or. .not.soil_par%constant_soil_par_therm &
                    .or. .not.soil_par%constant_soil_par_hydro) then
                    call soil_par_update(vc%carb%f_peat, vc%carb%litter_c, vc%carb%fast_c, vc%carb%slow_c, &
                        vc%carb%litter_c_peat, vc%carb%acro_c, vc%carb%cato_c, &
                        vc%soil%mineral_theta_sat, vc%soil%mineral_k_sat, vc%soil%mineral_psi_sat, &
                        vc%soil%mineral_Bi, vc%soil%mineral_lambda_s, vc%soil%mineral_lambda_dry, &
                        vc%carb%frac_soc, vc%soil%theta_sat, vc%soil%psi_sat, vc%soil%k_sat, &
                        vc%soil%k_exp, vc%soil%psi_exp, vc%soil%theta_field, vc%soil%theta_wilt, &
                        vc%soil%lambda_s, vc%soil%lambda_dry)
                endif
                if (veg_par%lroot_frac) call root_frac_update(vc%soil%alt, vc%soil%root_frac)
                call litter_in_frac_update(vc%soil%alt, vc%carb%litter_in_frac)
            endif
        endif

        ! vegetated fraction within the vc (=1 for the identity baseline) and
        ! the veg-mean skin temperature from the previous step (for snow albedo)
        f_veg = sum(vc%flx%frac_surf, mask=flag_veg.eq.1)
        t_skin_veg = 0._wp
        if (f_veg .gt. 0._wp) then
            do n = 1, nsurf
                if (flag_veg(n).eq.1) t_skin_veg = t_skin_veg + vc%flx%t_skin(n)*vc%flx%frac_surf(n)/f_veg
            end do
        else
            t_skin_veg = T0
        end if

        ! carry-over snapshots for the step
        vc%flx%t_skin_old  = vc%flx%t_skin
        vc%snow%w_snow_old = vc%snow%w_snow

        ! isotope tagging of throughfall input (VSMOW) when isotopes are active
        if (l_wiso) then
            vc%flx%rain_iso(:,i_o18) = Rstd(i_o18) * rain(:)
            vc%flx%snow_iso(:,i_o18) = Rstd(i_o18) * snow(:)
        endif

        ! --- phenology --------------------------------------------------------
        call phenology(vc%flx%frac_surf, f_veg, t2m, vc%veg%t2m_min_mon, vc%veg%gdd5_temp, vc%veg%gdd5, &
            vc%veg%gdd, vc%veg%phen_acc, vc%veg%phen, vc%veg%gamma_leaf, vc%veg%lai, vc%veg%lai_bal)

        ! --- aerodynamic resistance -------------------------------------------
        call resist_aer_veg(vc%flx%frac_surf, vc%veg%veg_h, vc%veg%lai, vc%veg%sai, vc%snow%h_snow, &
            tatm, vc%flx%t_skin, wind, &
            vc%flx%z0m, vc%flx%rough_m, vc%flx%rough_h, vc%flx%Ch, vc%flx%r_a, vc%flx%r_a_can, vc%flx%Ri)

        ! --- snow + surface albedo --------------------------------------------
        call snow_albedo_lake(t_skin_veg, vc%forc%snow, vc%snow%w_snow, vc%snow%w_snow_max, &
            vc%forc%dust, vc%forc%coszm, &
            vc%snow%alb_snow_vis_dir, vc%snow%alb_snow_vis_dif, vc%snow%alb_snow_nir_dir, vc%snow%alb_snow_nir_dif, &
            vc%snow%snow_grain, vc%snow%dust_con)

        call surface_albedo_veg(vc%flx%frac_surf, vc%desc%z_sur_std, vc%snow%h_snow, vc%forc%coszm, &
            vc%veg%lai, vc%veg%sai, vc%flx%z0m, vc%flx%f_snow_can, &
            vc%snow%alb_snow_vis_dir, vc%snow%alb_snow_vis_dif, vc%snow%alb_snow_nir_dir, vc%snow%alb_snow_nir_dif, &
            vc%forc%alb_bare_vis, vc%forc%alb_bare_nir, vc%flx%f_snow, &
            vc%flx%alb_vis_dir, vc%flx%alb_vis_dif, vc%flx%alb_nir_dir, vc%flx%alb_nir_dif, vc%flx%albedo)

        ! --- photosynthesis (canopy conductance) ------------------------------
        call photosynthesis(vc%forc%co2, vc%forc%c13_c12_atm, vc%forc%c14_c_atm, &
            vc%flx%frac_surf, t2m, vc%veg%t2m_min_mon, vc%veg%gdd5, vc%soil%t_soil(1:nl), &
            q2m, pressure, swnet, vc%flx%albedo, vc%forc%daylength, &
            vc%soil%theta_w, vc%soil%theta_field, vc%soil%theta_wilt, vc%soil%wilt, vc%soil%root_frac, &
            vc%veg%lai, vc%veg%phen, vc%veg%leaf_c, vc%veg%stem_c, vc%veg%root_c, &
            vc%veg%discrimination, vc%veg%ci, vc%veg%g_can, vc%veg%gpp, vc%veg%npp, vc%veg%npp13, vc%veg%npp14, &
            vc%veg%npp_cum, vc%veg%npp13_cum, vc%veg%npp14_cum, vc%veg%npp_ann, vc%veg%npp13_ann, vc%veg%npp14_ann, &
            vc%veg%aresp)

        ! --- surface resistance for evapotranspiration ------------------------
        call resist_sur_veg(vc%flx%frac_surf, vc%snow%mask_snow, vc%soil%theta_w, vc%veg%g_can, &
            vc%flx%beta_s, vc%flx%r_s, vc%flx%beta_s_can, vc%flx%r_s_can)

        ! --- canopy interception (throughfall + canopy evap) ------------------
        call canopy_water(vc%flx%frac_surf, vc%veg%lai, vc%veg%sai, vc%flx%r_a, vc%flx%t_skin, &
            pressure, qatm, rain, snow, &
            vc%flx%w_can, vc%flx%w_can_old, vc%flx%s_can, vc%flx%s_can_old, &
            vc%flx%rain_ground, vc%flx%snow_ground, vc%flx%evap_can, vc%flx%subl_can, vc%flx%f_wat_can, vc%flx%f_snow_can, &
            vc%flx%rain_iso, vc%flx%snow_iso, &
            vc%flx%w_can_iso, vc%flx%w_can_iso_old, vc%flx%s_can_iso, vc%flx%s_can_iso_old, &
            vc%flx%rain_ground_iso, vc%flx%snow_ground_iso, vc%flx%evap_can_iso, vc%flx%subl_can_iso)

        ! --- soil thermal properties ------------------------------------------
        call soil_par_thermal(vc%soil%t_soil, vc%soil%theta_w, vc%soil%theta_i, vc%soil%theta, vc%soil%theta_sat, &
            vc%soil%lambda_s, vc%soil%lambda_dry, vc%snow%h_snow, &
            vc%soil%cap_soil, vc%soil%lambda_soil, vc%soil%lambda_int_soil)

        ! --- surface energy balance -------------------------------------------
        call ebal_veg(vc%flx%frac_surf, vc%snow%mask_snow, vc%snow%h_snow, vc%snow%w_snow, vc%soil%lambda_soil, &
            vc%flx%evap_can, vc%flx%subl_can, vc%flx%t_skin, vc%flx%t_skin_old, vc%soil%t_soil, &
            tatm, qatm, pressure, swnet, swnet_min, lwdown, &
            vc%flx%beta_s, vc%flx%r_s, vc%flx%beta_s_can, vc%flx%r_s_can, vc%flx%r_a, vc%flx%r_a_can, &
            vc%flx%flx_g, vc%flx%dflxg_dT, vc%flx%flx_melt, flx_g_veg, dflxg_dT_veg, flx_melt_veg, vc%flx%t_skin_amp, &
            vc%flx%num_lh, vc%flx%num_sh, vc%flx%num_sw, vc%flx%num_lw, vc%flx%denom_lh, vc%flx%denom_sh, vc%flx%denom_lw, &
            vc%flx%f_sh, vc%flx%f_e, vc%flx%f_t, vc%flx%f_le, vc%flx%f_lt, vc%flx%f_lw, vc%flx%lh_ecan, &
            vc%flx%qsat_e, vc%flx%dqsatdT_e, vc%flx%qsat_t, vc%flx%dqsatdT_t, &
            encons1, ii, jj)

        ! isotope snapshot before the phase-change step
        if (l_wiso) then
            vc%soil%w_w_iso_old = vc%soil%w_w_iso
            vc%soil%w_i_iso_old = vc%soil%w_i_iso
            vc%snow%w_snow_iso_old = vc%snow%w_snow_iso
        endif

        ! --- soil temperature (phase change, snowmelt) ------------------------
        call soil_temp(vc%snow%mask_snow, vc%snow%h_snow, vc%soil%cap_soil, vc%soil%lambda_int_soil, &
            vc%soil%psi_sat, vc%soil%psi_exp, vc%soil%theta_sat, vc%carb%soil_resp_l, vc%carb%f_peat, &
            flx_g_veg, dflxg_dT_veg, flx_melt_veg, &
            vc%soil%t_soil, vc%soil%t_soil_cum, vc%snow%w_snow, vc%soil%w_w, vc%soil%w_i, &
            vc%soil%theta_w, vc%soil%theta_i, &
            vc%snow%snowmelt, vc%soil%t_soil_old, vc%snow%w_snow_old, &
            vc%soil%w_w_old, vc%soil%w_i_old, vc%soil%w_w_phase, vc%soil%w_i_phase, &
            vc%energy_cons_soil, ii, jj)

        ! --- skin temperature + surface fluxes --------------------------------
        call update_tskin_veg(vc%flx%frac_surf, vc%snow%mask_snow, vc%flx%t_skin_old, vc%flx%dflxg_dT, &
            tatm, qatm, swnet, lwdown, vc%soil%t_soil, vc%soil%t_soil_old, vc%flx%evap_can, vc%flx%subl_can, &
            vc%flx%flx_g, vc%flx%flx_melt, vc%flx%t_skin, t_skin_veg, &
            vc%flx%flx_sh, vc%flx%flx_lwu, vc%flx%lwnet, vc%flx%flx_lh, &
            vc%flx%evap_surface, vc%flx%transpiration, vc%flx%et, &
            vc%flx%num_lh, vc%flx%num_sh, vc%flx%num_sw, vc%flx%num_lw, vc%flx%denom_lh, vc%flx%denom_sh, vc%flx%denom_lw, &
            vc%flx%f_sh, vc%flx%f_e, vc%flx%f_t, vc%flx%f_le, vc%flx%f_lt, vc%flx%f_lw, vc%flx%lh_ecan, &
            vc%flx%qsat_e, vc%flx%dqsatdT_e, vc%flx%qsat_t, vc%flx%dqsatdT_t, &
            encons2, &
            vc%soil%w_w, vc%soil%w_w_iso, vc%soil%wilt, vc%flx%evap_can_iso, vc%flx%subl_can_iso, &
            vc%flx%evap_surface_iso, vc%flx%transpiration_iso, vc%flx%et_iso)

        ! --- surface hydrology (snow layer, wetland, runoff, infiltration) ----
        call surface_hydrology_veg(vc%flx%frac_surf, vc%snow%mask_snow, &
            vc%flx%evap_surface, vc%flx%rain_ground, vc%flx%snow_ground, vc%snow%snowmelt, &
            vc%soil%theta, vc%soil%theta_sat, vc%soil%theta_field, vc%soil%k_sat, vc%soil%cap_soil(1), &
            vc%soil%cti_mean, vc%soil%cti_cdf, &
            vc%soil%dyptop_k, vc%soil%dyptop_v, vc%soil%dyptop_xm, vc%soil%dyptop_fmax, &
            vc%snow%w_snow_old, vc%snow%w_snow, vc%snow%w_snow_max, vc%soil%w_w, vc%soil%w_i, &
            vc%soil%w_table_cum, vc%soil%f_wet_cum, vc%soil%t_soil, &
            vc%snow%h_snow, vc%soil%calving(1), vc%soil%runoff_sur(1), vc%soil%infiltration, &
            vc%soil%w_table, vc%soil%f_wet, vc%soil%f_wet_max, vc%soil%cti_lim, &
            vc%flx%evap_surface_iso, vc%flx%rain_ground_iso, vc%flx%snow_ground_iso, vc%snow%snowmelt_iso, &
            vc%snow%w_snow_iso, vc%soil%w_w_iso, vc%soil%w_i_iso, &
            vc%soil%calving_iso, vc%soil%runoff_sur_iso, vc%soil%infiltration_iso)

        ! --- soil hydraulic properties + soil water update (only if any soil) --
        if (f_veg .gt. 0._wp) then
            call soil_par_hydro(vc%soil%theta_w, vc%soil%theta_sat, &
                vc%soil%w_w, vc%soil%w_i, vc%soil%psi_sat, vc%soil%k_sat, vc%soil%k_exp, vc%soil%psi_exp, &
                vc%soil%theta, vc%soil%psi, vc%soil%kappa_int)

            call soil_hydro(vc%flx%frac_surf, vc%snow%mask_snow, vc%soil%theta_sat, &
                vc%soil%k_sat, vc%soil%k_exp, vc%soil%psi_exp, vc%soil%kappa_int, vc%soil%psi, &
                vc%soil%w_table, &
                vc%flx%transpiration, vc%flx%evap_surface, vc%soil%infiltration, vc%soil%wilt, &
                vc%snow%w_snow, vc%soil%w_w, vc%soil%w_i, &
                vc%soil%theta_w, vc%soil%theta_i, vc%soil%theta, vc%soil%theta_w_cum, vc%soil%theta_i_cum, vc%veg%theta_fire_cum, &
                vc%soil%drainage(1), &
                vc%flx%transpiration_iso, vc%flx%evap_surface_iso, vc%soil%infiltration_iso, &
                vc%soil%w_w_iso, vc%soil%drainage_iso)
        else
            vc%soil%drainage(1) = 0._wp
        endif

        ! --- potential evapotranspiration + cumulative water deficit ----------
        call calculate_pet(vc%flx%frac_surf, t2m, q2m, pressure, swnet, lwdown, vc%flx%flx_lwu, vc%flx%r_a, &
            pet)
        if (time_soy_lnd) vc%soil%cwd_mon(:) = 0._wp
        call calculate_cwd(vc%flx%frac_surf, rain, snow, pet, vc%soil%cwd_mon)
        if (time_eoy_lnd) vc%soil%mcwd = maxval(vc%soil%cwd_mon(:))
        ! veg-mean pet diagnostic
        vc%soil%pet = 0._wp
        if (f_veg .gt. 0._wp) then
            do n = 1, nsurf
                if (flag_veg(n).eq.1) vc%soil%pet = vc%soil%pet + pet(n)*vc%flx%frac_surf(n)/f_veg
            end do
        endif

        ! --- total surface + subsurface runoff over the veg column ------------
        vc%soil%runoff(1) = vc%soil%runoff_sur(1) + vc%soil%drainage(1)
        if (l_wiso) vc%soil%runoff_iso(:) = vc%soil%runoff_sur_iso(:) + vc%soil%drainage_iso(:)

        ! ===== soil-carbon decomposition parameters (port C.4) ================
        ! runs before dyn_veg (needs the year's cumulative soil T/moisture); sets
        ! the k_* decomposition rates soil_carbon/peat_carbon consume below.
        if (time_call_carb_p .and. vc%forc%f_veg_cell .gt. 0._wp) then
            call soil_carbon_par(vc%forc%f_veg_cell, vc%soil%theta_field, vc%soil%theta_sat, &
                vc%carb%litter_c_peat, vc%carb%acro_c, vc%carb%cato_c, vc%forc%dust, &
                vc%soil%t_soil_cum, vc%soil%theta_w_cum, vc%soil%theta_i_cum, vc%soil%psi, &
                vc%soil%f_wet_cum, vc%soil%w_table_cum, vc%soil%f_wet_mon, vc%soil%f_wet_long, vc%soil%w_table_mon, &
                vc%soil%t_soil_max, &
                vc%soil%frozen_years, vc%soil%thaw_timer, vc%carb%k_slow_to_fast, &
                vc%carb%ftemp, vc%carb%fmoist, vc%carb%fdepth, &
                vc%carb%k_litter, vc%carb%k_fast, vc%carb%k_slow, vc%carb%diff_soilc, vc%carb%adv_soilc, &
                vc%carb%k_litter_wet, vc%carb%k_fast_wet, vc%carb%k_slow_wet, &
                vc%carb%k_litter_peat, vc%carb%k_acro, vc%carb%k_cato, vc%carb%k_litter_peat_anox, vc%carb%k_acro_anox, &
                vc%carb%ch4_frac_wet, vc%carb%ch4_frac_peat, &
                vc%carb%f_peat_pot, vc%carb%f_oxic_peat, vc%soil%f_wetland, &
                vc%soil%w_table_min, vc%soil%w_table_peat, vc%soil%alt, &
                vc%carb%acro_h, vc%carb%cato_h, vc%carb%peat_c_ini_year)
        endif

        ! ===== end-of-month disturbance / fire parameters (port C.3) ==========
        if (time_eom_lnd) then
            call dynveg_par(vc%veg%disturbance, vc%veg%t2m_min_mon, vc%veg%gdd5, vc%veg%veg_c_above, &
                vc%forc%f_veg_cell, vc%carb%f_peat, vc%veg%pft_frac, &
                vc%carb%litter_c(1), vc%carb%litter_c_peat, vc%veg%lai_bal, vc%soil%mcwd, vc%soil%mcwd_clim, &
                vc%veg%theta_fire_cum, vc%veg%gamma_dist_cum, vc%veg%gamma_fire_cum, &
                vc%veg%fuel, vc%veg%f_fire_fuel, vc%veg%f_fire_cwd)
        endif

        ! ===== vegetation carbon + distribution (port C.3) ====================
        if (time_call_veg) then
            call dyn_veg(vc%forc%co2, vc%forc%f_veg_cell, vc%forc%f_veg_old_cell, &
                vc%forc%f_ice_grd_cell, vc%forc%f_ice_grd_old_cell, vc%forc%f_ice_nbr_cell, &
                vc%forc%f_lake_cell, vc%forc%f_lake_old_cell, vc%forc%f_shelf_cell, vc%forc%f_shelf_old_cell, &
                vc%veg%f_crop, vc%veg%f_pasture, vc%veg%gamma_luc, &
                vc%desc%z_sur_std, vc%veg%gamma_ice, &
                vc%veg%gamma_dist, vc%veg%gamma_dist_cum, vc%veg%gamma_fire, vc%veg%gamma_fire_cum, &
                vc%veg%npp_ann, vc%veg%npp13_ann, vc%veg%npp14_ann, &
                vc%veg%gamma_leaf, vc%veg%lambda, vc%veg%lai_bal, vc%veg%sai, &
                vc%soil%root_frac, vc%carb%litter_in_frac, &
                vc%veg%veg_c, vc%veg%veg_c13, vc%veg%veg_c14, vc%veg%leaf_c, vc%veg%stem_c, vc%veg%root_c, &
                vc%veg%seed_frac, vc%veg%pft_frac, &
                vc%veg%veg_c_above, vc%veg%veg_c13_above, vc%veg%veg_c14_above, &
                vc%veg%veg_c_below, vc%veg%veg_c13_below, vc%veg%veg_c14_below, &
                vc%veg%veg_h, vc%carb%litterfall, vc%carb%litterfall13, vc%carb%litterfall14, &
                vc%veg%npp_real, vc%veg%npp13_real, vc%veg%npp14_real, &
                vc%veg%fire_c_flux, vc%veg%fire_c13_flux, vc%veg%fire_c14_flux, &
                vc%veg%fire_c_flux_pft, vc%veg%fire_c13_flux_pft, vc%veg%fire_c14_flux_pft, &
                vc%veg%carbon_bal_veg, vc%veg%carbon13_bal_veg, vc%veg%carbon14_bal_veg, ii, jj)
            ! sync sub-tile fractions to the updated PFT fractions
            call set_frac_surf_veg(vc)
        endif

        ! ===== soil + peat carbon (port C.4) ==================================
        ! consumes the dyn_veg litterfall + soil_carbon_par rates; mineral pool
        ! (ic_min) then peatland pool (ic_peat).
        if (time_call_carb) then
            if (vc%forc%f_veg_cell .gt. 0._wp) then
                call soil_carbon(vc%forc%f_veg_cell, vc%soil%f_wetland, vc%carb%f_peat, vc%veg%f_crop, vc%veg%f_pasture, &
                    vc%carb%litterfall(:,ic_min), vc%carb%litterfall13(:,ic_min), vc%carb%litterfall14(:,ic_min), &
                    vc%carb%litter_c, vc%carb%fast_c, vc%carb%slow_c, &
                    vc%carb%litter_c13, vc%carb%fast_c13, vc%carb%slow_c13, &
                    vc%carb%litter_c14, vc%carb%fast_c14, vc%carb%slow_c14, &
                    vc%carb%k_litter, vc%carb%k_fast, vc%carb%k_slow, vc%carb%k_slow_to_fast, &
                    vc%carb%k_litter_wet, vc%carb%k_fast_wet, vc%carb%k_slow_wet, vc%carb%diff_soilc, vc%carb%adv_soilc, vc%carb%ch4_frac_wet, &
                    vc%carb%soil_resp(ic_min), vc%carb%soil_resp13(ic_min), vc%carb%soil_resp14(ic_min), vc%carb%soil_resp_l(:,ic_min), &
                    vc%carb%soil_c_tot(ic_min), vc%carb%soil_c13_tot(ic_min), vc%carb%soil_c14_tot(ic_min), &
                    vc%carb%ch4_emis_wetland, vc%carb%c13h4_emis_wetland, &
                    vc%carb%carbon_cons_soil(ic_min), vc%carb%carbon13_cons_soil(ic_min), vc%carb%carbon14_cons_soil(ic_min))
            endif
            if (vc%forc%f_veg_cell .gt. 0._wp .and. peat_par%peat_carb) then
                call peat_carbon(vc%carb%f_oxic_peat, &
                    vc%carb%litterfall(:,ic_peat), vc%carb%litterfall13(:,ic_peat), vc%carb%litterfall14(:,ic_peat), &
                    vc%carb%litter_c_peat, vc%carb%acro_c, vc%carb%cato_c, &
                    vc%carb%litter_c13_peat, vc%carb%acro_c13, vc%carb%cato_c13, &
                    vc%carb%litter_c14_peat, vc%carb%acro_c14, vc%carb%cato_c14, &
                    vc%carb%k_litter_peat, vc%carb%k_acro, vc%carb%k_cato, &
                    vc%carb%k_litter_peat_anox, vc%carb%k_acro_anox, vc%carb%ch4_frac_peat, &
                    vc%carb%soil_resp(ic_peat), vc%carb%soil_resp13(ic_peat), vc%carb%soil_resp14(ic_peat), vc%carb%soil_resp_l(:,ic_peat), &
                    vc%carb%soil_c_tot(ic_peat), vc%carb%soil_c13_tot(ic_peat), vc%carb%soil_c14_tot(ic_peat), &
                    vc%carb%ch4_emis_peat, vc%carb%c13h4_emis_peat, &
                    vc%carb%carbon_cons_soil(ic_peat), vc%carb%carbon13_cons_soil(ic_peat), vc%carb%carbon14_cons_soil(ic_peat), &
                    vc%carb%peat_c_ini_year, vc%carb%dCpeat_dt)
            endif

            ! N2O emissions (port C.5) — combined mineral+peat soil respiration
            ! drives nitrification/denitrification (reference lnd_model order).
            if (vc%forc%f_veg_cell .gt. 0._wp) then
                soil_resp = vc%carb%soil_resp(ic_min)*(vc%forc%f_veg_cell-vc%carb%f_peat) &
                          + vc%carb%soil_resp(ic_peat)*vc%carb%f_peat
                call n2o_emission(soil_resp, vc%soil%t_soil(1), vc%soil%theta_w(1), &
                    vc%soil%theta_field(1), vc%soil%theta_sat(1), vc%carb%n2o_emis)
            else
                vc%carb%n2o_emis = 0._wp
            endif
        endif

        ! ===== dust emissions (port C.5) =====================================
        ! runs every land step (not time-gated), guarded f_veg>0 (reference lnd
        ! order, after the carbon block). z_veg_std uses desc%z_sur_std (as
        ! dyn_veg); tatm is the sub-tile broadcast; wind is the bare-tile scalar.
        if (vc%forc%f_veg_cell .gt. 0._wp) then
            call dust_emission(vc%flx%frac_surf, vc%desc%z_sur_std, &
                vc%forc%z_veg, vc%forc%z_veg_min, vc%forc%z_veg_max, vc%flx%f_snow(i_bare), &
                vc%veg%lai, vc%veg%sai, vc%snow%h_snow, vc%flx%t_skin, tatm, &
                vc%soil%theta_w(1), vc%soil%theta_i(1), vc%forc%wind, &
                vc%carb%dust_emis_d, vc%carb%dust_emis_g, vc%carb%dust_emis_s, vc%carb%dust_emis)
        else
            vc%carb%dust_emis_d = 0._wp
            vc%carb%dust_emis_g = 0._wp
            vc%carb%dust_emis_s = 0._wp
            vc%carb%dust_emis   = 0._wp
        endif

        ! ===== weathering + river carbon export (port W.4), end-of-year =======
        ! land-cell operations (reference lnd_model order, after dust). Weathering
        ! scales with f_veg internally; carbon_export leaches the veg+peat pools.
        ! Outputs land on the carb block (fed to the ocean carbon cycle later).
        if (time_eoy_lnd) then
            if (i_weathering.eq.1) then
                call weathering_gemco2(vc%forc%c13_c12_atm, vc%forc%c14_c_atm, vc%forc%weath_scale, &
                    vc%forc%f_veg_cell, vc%carb%lithology_gemco2, vc%soil%runoff_ann, &
                    vc%carb%weath_carb, vc%carb%weath13_carb, vc%carb%weath14_carb, &
                    vc%carb%weath_sil, vc%carb%weath13_sil, vc%carb%weath14_sil, vc%carb%weath_loess)
            else if (i_weathering.eq.2) then
                call weathering_uhh(vc%forc%c13_c12_atm, vc%forc%c14_c_atm, vc%forc%weath_scale, &
                    vc%forc%f_veg_cell, vc%carb%lithology_uhh, vc%soil%runoff_ann, vc%veg%t2m_ann_mean, &
                    vc%carb%weath_carb, vc%carb%weath13_carb, vc%carb%weath14_carb, &
                    vc%carb%weath_sil, vc%carb%weath13_sil, vc%carb%weath14_sil, vc%carb%weath_loess)
            endif
            if (l_river_export) then
                call carbon_export(vc%soil%runoff_ann, vc%forc%f_veg_cell, vc%carb%f_peat, &
                    vc%carb%litter_c, vc%carb%litter_c13, vc%carb%litter_c14, &
                    vc%carb%fast_c, vc%carb%fast_c13, vc%carb%fast_c14, &
                    vc%carb%slow_c, vc%carb%slow_c13, vc%carb%slow_c14, &
                    vc%carb%litter_c_peat, vc%carb%litter_c13_peat, vc%carb%litter_c14_peat, &
                    vc%carb%acro_c, vc%carb%acro_c13, vc%carb%acro_c14, &
                    vc%carb%cato_c, vc%carb%cato_c13, vc%carb%cato_c14, &
                    vc%carb%poc_export, vc%carb%poc13_export, vc%carb%poc14_export, &
                    vc%carb%doc_export, vc%carb%doc13_export, vc%carb%doc14_export)
            endif
        endif

        return

    end subroutine lndvc_update_land

    subroutine lndvc_update_lake(vc)
        ! Single-column lake-surface update (pattern A): unpack the vc's lake /
        ! snow / flux blocks into the ported single-tile lake chain, mirroring the
        ! reference lnd lake path (surface params -> lake thermal params -> energy
        ! balance -> lake temperature -> sublake soil thermal -> skin update ->
        ! surface hydrology). The sublake soil column (sublake_par_thermal/
        ! sublake_temp, port L.1) runs on mineral soil params seeded from the
        ! reference cell; its soil_resp_l heat source is zero until lake carbon
        ! is wired (L.3). t_skin_old is snapshotted inside ebal_lake; the melt-iso redistribution
        ! that lives in the reference lnd_model wrapper is deferred (like SMB iso).

        implicit none

        type(vc_t), intent(inout) :: vc

        integer, parameter :: ii = 0, jj = 0   ! debug-print indices only
        real(wp) :: z0m_lake
        real(wp) :: calving, runoff_sur              ! lake calving/runoff not yet aggregated
        real(wp) :: calving_iso(nwiso), runoff_sur_iso(nwiso)
        real(wp) :: litter_lake_zero(nlc)            ! litter input to lake (zero: cross-class, deferred to wrapper)

        ! carry-over snapshot for the step
        vc%snow%w_snow_old = vc%snow%w_snow

        ! no canopy over lake: throughfall = precipitation (tag iso at VSMOW)
        vc%flx%rain_ground(1) = vc%forc%rain
        vc%flx%snow_ground(1) = vc%forc%snow
        if (l_wiso) then
            vc%flx%rain_ground_iso(1,i_o18) = Rstd(i_o18) * vc%forc%rain
            vc%flx%snow_ground_iso(1,i_o18) = Rstd(i_o18) * vc%forc%snow
        endif

        ! lake surface momentum roughness (fixed lake/ice value)
        z0m_lake = lnd_surf_par%z0m_lake_ice
        vc%flx%z0m(1) = z0m_lake

        ! --- surface parameters -----------------------------------------------
        call resist_aer_lake(vc%snow%h_snow, vc%forc%tatm, vc%flx%t_skin(1), vc%forc%wind, &
            z0m_lake, vc%flx%rough_m(1), vc%flx%rough_h(1), vc%flx%Ch(1), vc%flx%r_a(1), vc%flx%Ri(1))

        call snow_albedo_lake(vc%flx%t_skin(1), vc%forc%snow, vc%snow%w_snow, vc%snow%w_snow_max, &
            vc%forc%dust, vc%forc%coszm, &
            vc%snow%alb_snow_vis_dir, vc%snow%alb_snow_vis_dif, vc%snow%alb_snow_nir_dir, vc%snow%alb_snow_nir_dif, &
            vc%snow%snow_grain, vc%snow%dust_con)

        call surface_albedo_lake(vc%snow%h_snow, vc%forc%coszm, vc%lake%f_lake_ice, &
            vc%snow%alb_snow_vis_dir, vc%snow%alb_snow_vis_dif, vc%snow%alb_snow_nir_dir, vc%snow%alb_snow_nir_dif, &
            vc%snow%f_snow, vc%flx%alb_vis_dir(1), vc%flx%alb_vis_dif(1), vc%flx%alb_nir_dir(1), vc%flx%alb_nir_dif(1), &
            vc%flx%albedo(1))

        call resist_sur_lake(vc%flx%beta_s(1), vc%flx%r_s(1))

        ! --- lake thermal parameters ------------------------------------------
        call lake_par_thermal(vc%snow%h_snow, vc%lake%h_lake, vc%lake%t_lake(1:nl_l), vc%lake%f_i_lake, &
            vc%forc%wind, vc%desc%lat, &
            vc%lake%cap_lake, vc%lake%lambda_lake, vc%lake%lambda_int_lake)

        ! --- surface energy balance -------------------------------------------
        call ebal_lake(vc%snow%mask_snow, vc%snow%h_snow, vc%lake%lambda_lake, &
            vc%flx%t_skin(1), vc%flx%t_skin_old(1), vc%lake%t_lake, vc%forc%tatm, vc%forc%qatm, vc%forc%pressure, &
            vc%forc%swnet, vc%forc%lwdown, &
            vc%flx%beta_s(1), vc%flx%r_s(1), vc%flx%r_a(1), &
            vc%flx%flx_g(1), vc%flx%dflxg_dT(1), vc%flx%flx_melt(1), vc%flx%t_skin_amp(1), &
            vc%flx%num_lh(1), vc%flx%num_sh(1), vc%flx%num_sw(1), vc%flx%num_lw(1), &
            vc%flx%denom_lh(1), vc%flx%denom_sh(1), vc%flx%denom_lw(1), &
            vc%flx%f_sh(1), vc%flx%f_e(1), vc%flx%f_le(1), vc%flx%f_lw(1), vc%flx%qsat_e(1), vc%flx%dqsatdT_e(1), &
            vc%energy_cons_surf1, ii, jj)

        ! --- lake temperature (convection handled inside lake_temp) ------------
        call lake_temp(vc%snow%mask_snow, vc%snow%h_snow, vc%lake%h_lake, vc%lake%cap_lake, vc%lake%lambda_int_lake, &
            vc%flx%flx_g(1), vc%flx%dflxg_dT(1), vc%flx%flx_melt(1), &
            vc%lake%t_lake, vc%snow%w_snow, vc%lake%w_w_lake, vc%lake%w_i_lake, vc%lake%f_i_lake, vc%lake%f_lake_ice, &
            vc%snow%snowmelt, vc%lake%t_lake_old, vc%snow%w_snow_old, &
            vc%lake%h_lake_conv, vc%lake%h_lake_mix, &
            vc%lake%energy_cons_lake, ii, jj)

        ! --- sublake soil column (thermal, port L.1) --------------------------
        ! bottom lake-layer temperature is the top boundary for the sublake soil
        ! profile; soil params are the seeded mineral values, soil_resp_l is the
        ! lake-carbon respiration heat (zero until lake carbon is wired in L.3).
        vc%lake%t_sublake(0) = vc%lake%t_lake(nl_l)
        call sublake_par_thermal(vc%lake%theta_w_sublake, vc%lake%theta_i_sublake, &
            vc%lake%theta_sat, vc%lake%lambda_s, vc%lake%cap_sublake, vc%lake%lambda_int_sublake)
        call sublake_temp(vc%lake%cap_sublake, vc%lake%lambda_int_sublake, &
            vc%lake%psi_sat, vc%lake%psi_exp, vc%lake%theta_sat, vc%lake%soil_resp_l, &
            vc%lake%t_sublake, vc%lake%w_w_sublake, vc%lake%w_i_sublake, &
            vc%lake%theta_w_sublake, vc%lake%theta_i_sublake, &
            vc%lake%t_sublake_cum, vc%lake%theta_w_sublake_cum, vc%lake%theta_i_sublake_cum, &
            vc%lake%energy_cons_lake, ii, jj)

        ! --- skin temperature + surface fluxes --------------------------------
        call update_tskin_lake(vc%snow%mask_snow, vc%flx%t_skin_old(1), vc%flx%dflxg_dT(1), &
            vc%forc%tatm, vc%forc%qatm, vc%forc%swnet, vc%forc%lwdown, &
            vc%lake%t_lake, vc%lake%t_lake_old, vc%flx%flx_g(1), vc%flx%flx_melt(1), &
            vc%flx%t_skin(1), vc%flx%flx_sh(1), vc%flx%flx_lwu(1), vc%flx%flx_lh(1), vc%flx%evap_surface(1), vc%flx%et(1), &
            vc%flx%num_lh(1), vc%flx%num_sh(1), vc%flx%num_sw(1), vc%flx%num_lw(1), &
            vc%flx%denom_lh(1), vc%flx%denom_sh(1), vc%flx%denom_lw(1), &
            vc%flx%f_sh(1), vc%flx%f_e(1), vc%flx%f_le(1), vc%flx%f_lw(1), vc%flx%qsat_e(1), vc%flx%dqsatdT_e(1), &
            vc%energy_cons_surf2, ii, jj, &
            vc%snow%w_snow, vc%snow%w_snow_iso, vc%lake%w_w_lake(1), vc%lake%w_w_lake_iso(1,:), &
            vc%flx%evap_surface_iso(1,:), vc%flx%et_iso(1,:))

        ! --- surface hydrology (snow budget + lake water balance -> runoff) ----
        call surface_hydrology_lake(vc%snow%mask_snow, vc%flx%evap_surface(1), vc%flx%snow_ground(1), vc%flx%rain_ground(1), &
            vc%snow%snowmelt, vc%lake%cap_lake(1), vc%lake%t_lake(1), vc%snow%w_snow_old, vc%snow%w_snow, vc%snow%w_snow_max, &
            vc%snow%h_snow, calving, runoff_sur, vc%lake%lake_water_tendency, &
            vc%flx%evap_surface_iso(1,:), vc%flx%snow_ground_iso(1,:), vc%flx%rain_ground_iso(1,:), &
            vc%snow%snowmelt_iso, vc%snow%w_snow_iso, calving_iso, runoff_sur_iso)

        ! --- lake sediment carbon (port L.3) ----------------------------------
        ! decomposition rates from the cumulative sublake T/moisture, then the
        ! carbon pool update. litterfall into the lake is a cross-class quantity
        ! (carbon_trans), zero here — deferred to the cell-level wrapper; at
        ! static n_vc=1 geometry carbon_trans is a no-op, so identity holds.
        ! No f_lake guard: the lake vc's existence is the gate (as the surface
        ! chain above). ch4_emis_lake stays per-vc, aggregation deferred.
        if (time_call_carb_p) then
            call lake_carbon_par(vc%lake%t_sublake_cum, vc%lake%theta_w_sublake_cum, vc%lake%theta_i_sublake_cum, &
                vc%lake%k_litter_lake, vc%lake%k_fast_lake, vc%lake%k_slow_lake, &
                vc%lake%diff_lakec, vc%lake%adv_lakec, vc%lake%ch4_frac_lake)
        endif
        if (time_call_carb) then
            litter_lake_zero(:) = 0._wp
            call lake_carbon(litter_lake_zero, litter_lake_zero, litter_lake_zero, vc%lake%ch4_frac_lake, &
                vc%lake%litter_c_lake, vc%lake%fast_c_lake, vc%lake%slow_c_lake, &
                vc%lake%litter_c13_lake, vc%lake%fast_c13_lake, vc%lake%slow_c13_lake, &
                vc%lake%litter_c14_lake, vc%lake%fast_c14_lake, vc%lake%slow_c14_lake, &
                vc%lake%k_litter_lake, vc%lake%k_fast_lake, vc%lake%k_slow_lake, vc%lake%diff_lakec, vc%lake%adv_lakec, &
                vc%lake%soil_resp_lake, vc%lake%soil_resp13_lake, vc%lake%soil_resp14_lake, vc%lake%soil_resp_l(:,ic_lake), &
                vc%lake%soil_c_tot_lake, vc%lake%soil_c13_tot_lake, vc%lake%soil_c14_tot_lake, &
                vc%lake%ch4_emis_lake, vc%lake%c13h4_emis_lake, &
                vc%lake%carbon_cons_lake, vc%lake%carbon13_cons_lake, vc%lake%carbon14_cons_lake)
        endif

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
