module lndvc_grid
    ! Re-export shim for the reference lnd_grid symbols.
    !
    ! Historically this module carried its own hard-coded soil grid
    ! (parameter nl = 5, z = 0, 0.1, 0.3, 0.7, 1.5, 3.1) alongside all the
    ! surface-type / PFT / carbon-class parameters. That created a fatal
    ! size mismatch with the reference lnd grid (nl variable, read from
    ! lnd_par.nml): the lndvc-framework top-level (lndvc_model, lndvc_decomp)
    ! saw nl = 5 while every file in src/lndvc/lnd/ pulled nl from lnd_grid
    ! (currently 7), so soil-column arrays were allocated at size 6 but the
    ! physics tridiag operated on 8 layers, corrupting memory and producing
    ! NaN in the first soil_temp call.
    !
    ! Rather than duplicate the grid, this module now re-exports every symbol
    ! from lnd_grid that the lndvc code uses, so the framework has a single
    ! source of truth. lnd_grid_init (called by lnd_init before lndvc_init in
    ! climber.f90) sets everything up; lndvc_grid_init is now a no-op kept
    ! only so the existing `call lndvc_grid_init()` in lndvc_init stays valid.

    use lnd_grid, only : nl, z, z_int, dz, rdz, rdz_pos, rdz_neg
    use lnd_grid, only : nl_l, z_l, z_int_l, dz_l, rdz_l, rdz_pos_l, rdz_neg_l
    use lnd_grid, only : nlc, z_c, z_int_c, dz_c, rdz_c, rdz_pos_c, rdz_neg_c
    use lnd_grid, only : nsurf, flag_pft, flag_veg, i_surf, i_bare, i_lake, i_ice
    use lnd_grid, only : npft, ntrees, ngrass, nshrub
    use lnd_grid, only : i_pft, i_trees, i_grass, i_shrub
    use lnd_grid, only : flag_tree, flag_grass, flag_shrub, nveg
    use lnd_grid, only : ncarb, ic_min, ic_peat, ic_shelf, ic_ice, ic_lake
    use lnd_grid, only : nsoil, i_soil, is_veg, is_ice, is_lake

    implicit none
    public

contains

    subroutine lndvc_grid_init
        ! No-op: lnd_grid_init (called from lnd_init, which runs before
        ! lndvc_init in climber.f90) already sets up every symbol re-exported
        ! above. Kept as an entry point so lndvc_init can retain its call.
        implicit none
        return
    end subroutine lndvc_grid_init

end module lndvc_grid
