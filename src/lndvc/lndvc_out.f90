!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module : l n d v c _ o u t
!
!  Purpose : diagnostics + netCDF output for the land-virtual-cell
!            framework. Parallel to src/lnd/lnd_out.f90 for the
!            reference land model, but consumes lnd%cell(i,j)%* — the
!            per-cell aggregate that lndvc_aggregate_cell populates
!            each timestep — instead of the reference lnd%l2d(i,j)
!            fields. Produces three files under out_dir:
!
!               lndvc_ts.nc     annual global-integrated scalars
!               lndvc_surf.nc   monthly-climatology + annual 2D fields
!                               on the coarse cmn grid  (write_surf)
!               lndvc_hires.nc  annual hires 2D fields painted from
!                               the coarse cmn grid    (write_hires)
!
!            At the identity baseline (n_vc=1 per class, no hypsometry)
!            the hires fields are cell-uniform within each cmn cell —
!            i.e. blocky. The plumbing is being set up here for the
!            eventual hypsometric split.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module lndvc_out

    use precision,     only : wp, sp
    use ncio
    use dim_name,      only : dim_time, dim_lon, dim_lat, dim_month
    use timer,         only : year, mon, year_now, year_clim, &
                              sec_year, nmon_year
    use timer,         only : time_soy_lndvc, time_eoy_lndvc, &
                              time_out_lndvc, time_out_ts_clim, &
                              nstep_mon_lnd, ny_out_ts, y_out_ts_clim, &
                              n_accel
    use climber_grid,  only : ni, nj, area, lon, lat
    use lnd_params,    only : dt
    use lndvc_def,     only : lndvc_class, vc_cell_t
    use lndvc_params,  only : lndvc_par_init, &
                              write_surf, write_hires, &
                              hires_write_h_snow, hires_write_w_snow, &
                              hires_write_t_skin, hires_write_smb, &
                              hires_write_melt, hires_write_f_ice, &
                              hires_write_veg_c, hires_write_soil_c
    use control,       only : out_dir
    use constants,     only : T0
    use geo_def,       only : geo_class

    implicit none

    private
    public :: lndvc_diag_init, lndvc_diag

    real(sp), parameter :: missing_value = -9999.0_sp

    integer :: nout   ! surf/hires write counter (grows by 1 per output frame)

    !--- annual-scalar accumulator (mirror of the useful subset of lnd_ts) ---
    type ts_out
        real(wp) :: land, veg, ice, lake, peat         ! areas [mln km2]
        real(wp) :: awet, aperm                        ! diagnostic areas
        real(wp) :: prc, evp, esur, trans, run         ! water fluxes [10^15 kg/yr]
        real(wp) :: runsur, drain                      ! runoff components
        real(wp) :: gpp, npp                           ! carbon fluxes [PgC/yr]
        real(wp) :: veg_c, soil_c, peat_c              ! carbon stocks [PgC]
        real(wp) :: Cflx_atm_lnd                       ! net atm-lnd C flux [PgC/yr]
        real(wp) :: ch4_emis                           ! TgCH4/yr
        real(wp) :: n2o_emis                           ! TgN2O-N/yr
        real(wp) :: snow_vol                           ! total snow volume [10^12 m^3]
        real(wp) :: t2m                                ! land-mean 2 m temperature [K]
        real(wp) :: h_snow                             ! area-mean snow depth [m]
        real(wp) :: lai                                ! area-mean LAI
    end type ts_out

    !--- 2D surface accumulator on the cmn (ni,nj) grid ------------------
    type surf_out
        real(wp), allocatable, dimension(:,:) :: f_land, f_veg, f_ice, f_lake, f_peat
        real(wp), allocatable, dimension(:,:) :: t_skin, albedo, t2m, q2m
        real(wp), allocatable, dimension(:,:) :: swnet, lwnet, flx_sh, flx_lh, flx_g
        real(wp), allocatable, dimension(:,:) :: h_snow, w_snow, w_snow_max
        real(wp), allocatable, dimension(:,:) :: snowmelt, icemelt, f_snow
        real(wp), allocatable, dimension(:,:) :: rain, snow_flx
        real(wp), allocatable, dimension(:,:) :: transpiration, evap_surface, evap_can, et
        real(wp), allocatable, dimension(:,:) :: runoff, runoff_sur, runoff_gw, drainage
        real(wp), allocatable, dimension(:,:) :: w_table, w_table_peat
        real(wp), allocatable, dimension(:,:) :: f_wet, f_wetland
        real(wp), allocatable, dimension(:,:) :: lai, sai, gpp, npp, veg_c
        real(wp), allocatable, dimension(:,:) :: gdd5, t2m_min_mon
        real(wp), allocatable, dimension(:,:) :: t_soil_top, theta_w_top, alt
        real(wp), allocatable, dimension(:,:) :: soil_c, soil_resp
        real(wp), allocatable, dimension(:,:) :: smb, melt
        real(wp), allocatable, dimension(:,:) :: Cflx_atm_lnd, ch4_emis, n2o_emis
    end type surf_out

    !--- hires accumulator (single annual snapshot painted from cmn) -----
    type hires_out
        real(wp), allocatable, dimension(:,:) :: h_snow, w_snow, t_skin
        real(wp), allocatable, dimension(:,:) :: smb, melt, f_ice
        real(wp), allocatable, dimension(:,:) :: veg_c, soil_c
    end type hires_out

    ! module-level accumulators
    type(ts_out),   allocatable :: ann_ts(:)         ! (ny_out_ts)
    type(surf_out), allocatable :: mon_su(:)         ! (nmon_year)
    type(surf_out)              :: ann_su
    type(hires_out)             :: hires_su

    ! hires grid metadata (captured at init from geo%hires%grid)
    integer :: nx_hi = 0, ny_hi = 0
    integer, allocatable :: ij_hi_i(:,:), ij_hi_j(:,:) ! (nx_hi,ny_hi) coarse index

    ! per-year running totals (rebuilt each year_soy)
    real(wp) :: acc_prc, acc_evp, acc_esur, acc_trans
    real(wp) :: acc_run, acc_runsur, acc_drain
    real(wp) :: acc_gpp, acc_npp
    real(wp) :: acc_ch4, acc_n2o
    real(wp) :: acc_cflx
    real(wp) :: acc_t2m_num, acc_t2m_den
    real(wp) :: acc_snowvol
    real(wp) :: acc_hsnow_num, acc_hsnow_den
    real(wp) :: acc_lai_num, acc_lai_den

contains

    !======================================================================
    ! l n d v c _ d i a g _ i n i t
    !----------------------------------------------------------------------
    subroutine lndvc_diag_init(lnd, geo)

        implicit none

        type(lndvc_class), intent(in) :: lnd
        type(geo_class),   intent(in) :: geo

        integer :: k

        ! parse lndvc_par.nml (mirrors smb_par_init idiom).
        call lndvc_par_init()

        nout = 0

        ! allocate ts buffer
        allocate(ann_ts(ny_out_ts))

        ! allocate monthly + annual surf accumulators
        allocate(mon_su(nmon_year))
        do k = 1, nmon_year
            call alloc_surf(mon_su(k))
        end do
        call alloc_surf(ann_su)

        ! allocate hires accumulator + build coarse-index map
        if (write_hires) then
            nx_hi = geo%hires%grid%G%nx
            ny_hi = geo%hires%grid%G%ny
            call alloc_hires(hires_su)
            call build_hires_map(geo)
        end if

        ! create NC files with dimension defs
        call ts_nc(trim(out_dir)//"/lndvc_ts.nc")
        if (write_surf ) call surf_nc(trim(out_dir)//"/lndvc_surf.nc")
        if (write_hires) call hires_nc(trim(out_dir)//"/lndvc_hires.nc", geo)

        ! zero per-year running totals so year 1 starts clean before the
        ! first time_soy_lndvc reset arrives
        call reset_year_accum()

        ! zero the mon_su + ann_su accumulators
        do k = 1, nmon_year
            call zero_surf(mon_su(k))
        end do
        call zero_surf(ann_su)

        return

    end subroutine lndvc_diag_init


    !======================================================================
    ! l n d v c _ d i a g
    !----------------------------------------------------------------------
    ! Called each timestep inside the flag_lndvc block, after lndvc_update
    ! has populated lnd%cell(:,:) for the current step. Accumulates monthly
    ! surface averages, annual scalars, and writes to disk at eoy when the
    ! output cadence hits.
    !======================================================================
    subroutine lndvc_diag(lnd, geo)

        implicit none

        type(lndvc_class), intent(in) :: lnd
        type(geo_class),   intent(in) :: geo

        integer  :: i, j, y
        real(wp) :: mon_avg
        real(wp) :: t_land, area_land
        real(wp) :: h_num, h_den
        real(wp) :: lai_num, lai_den
        real(wp) :: snow_vol_step

        y       = y_out_ts_clim
        mon_avg = 1._wp / real(nstep_mon_lnd, wp)

        ! -------------------------------------------------------------
        ! Year-boundary reset. Zero the monthly accumulators and the
        ! per-year running totals.
        ! -------------------------------------------------------------
        if (time_soy_lndvc) then
            do i = 1, nmon_year
                call zero_surf(mon_su(i))
            end do
            call zero_surf(ann_su)
            call reset_year_accum()
        end if

        ! -------------------------------------------------------------
        ! Per-step 2D accumulation into the current month's mon_su. Each
        ! contribution is scaled by mon_avg = 1/nstep_mon_lnd so mon_su
        ! ends the month as a proper time-mean.
        ! -------------------------------------------------------------
        !$omp parallel do collapse(2) private(i,j) schedule(static)
        do j = 1, nj
            do i = 1, ni
                if (lnd%cell(i,j)%mask_lnd .eq. 1) then
                    call accum_surf_cell(mon_su(mon), lnd%cell(i,j), i, j, mon_avg)
                end if
            end do
        end do
        !$omp end parallel do

        ! -------------------------------------------------------------
        ! Per-step contribution to the annual ts running totals. Same
        ! pattern as lnd_out: rate fluxes get *dt (accumulated over the
        ! whole year for a total), state variables get /nstep_year_lnd
        ! via mon_avg / nmon_year.
        ! -------------------------------------------------------------
        block
        real(wp) :: l_prc, l_evp, l_esur, l_trans, l_run, l_runsur, l_drain
        real(wp) :: l_gpp, l_npp, l_ch4, l_n2o

        t_land     = 0._wp
        area_land  = 0._wp
        h_num      = 0._wp
        h_den      = 0._wp
        lai_num    = 0._wp
        lai_den    = 0._wp
        snow_vol_step = 0._wp
        l_prc = 0._wp; l_evp = 0._wp; l_esur = 0._wp; l_trans = 0._wp
        l_run = 0._wp; l_runsur = 0._wp; l_drain = 0._wp
        l_gpp = 0._wp; l_npp = 0._wp; l_ch4 = 0._wp; l_n2o = 0._wp

        !$omp parallel do collapse(2) private(i,j) &
        !$omp reduction(+:t_land,area_land,h_num,h_den,lai_num,lai_den,snow_vol_step) &
        !$omp reduction(+:l_prc,l_evp,l_esur,l_trans,l_run,l_runsur,l_drain) &
        !$omp reduction(+:l_gpp,l_npp,l_ch4,l_n2o) schedule(static)
        do j = 1, nj
            do i = 1, ni
                if (lnd%cell(i,j)%mask_lnd .eq. 1) then
                    ! precipitation [10^15 kg/yr]
                    l_prc    = l_prc + (lnd%cell(i,j)%rain + lnd%cell(i,j)%snow_flx) &
                                     * area(i,j) * 1.e-15_wp * dt
                    ! total evapotranspiration
                    l_evp    = l_evp + lnd%cell(i,j)%et * area(i,j) * 1.e-15_wp * dt
                    ! surface evaporation
                    l_esur   = l_esur + lnd%cell(i,j)%evap_surface * area(i,j) * 1.e-15_wp * dt
                    ! transpiration
                    l_trans  = l_trans + lnd%cell(i,j)%transpiration * area(i,j) * 1.e-15_wp * dt
                    ! total runoff
                    l_run    = l_run + lnd%cell(i,j)%runoff * area(i,j) * 1.e-15_wp * dt
                    ! surface runoff
                    l_runsur = l_runsur + lnd%cell(i,j)%runoff_sur * area(i,j) * 1.e-15_wp * dt
                    ! drainage
                    l_drain  = l_drain + lnd%cell(i,j)%drainage * area(i,j) * 1.e-15_wp * dt
                    ! gpp / npp [PgC/yr] — kgC/m2/s * m2 * s (dt) * 1e-12
                    l_gpp    = l_gpp + lnd%cell(i,j)%gpp * area(i,j) * 1.e-12_wp * dt
                    l_npp    = l_npp + lnd%cell(i,j)%npp * area(i,j) * 1.e-12_wp * dt
                    ! ch4/n2o emissions — kgC/m2/s * m2 * s (dt) * 1e-9 [Tg/yr]
                    l_ch4    = l_ch4 + lnd%cell(i,j)%ch4_emis * area(i,j) * 1.e-9_wp * dt
                    l_n2o    = l_n2o + lnd%cell(i,j)%n2o_emis * area(i,j) * 1.e-9_wp * dt

                    ! 2 m land temperature: area-weighted, land only. Divisor
                    ! is applied after year sum.
                    t_land    = t_land    + lnd%cell(i,j)%t2m * lnd%cell(i,j)%f_veg * area(i,j)
                    area_land = area_land + lnd%cell(i,j)%f_veg * area(i,j)

                    ! area-mean snow depth over land
                    h_num = h_num + lnd%cell(i,j)%h_snow * lnd%cell(i,j)%f_land * area(i,j)
                    h_den = h_den + lnd%cell(i,j)%f_land * area(i,j)

                    ! area-mean LAI over the veg fraction
                    lai_num = lai_num + lnd%cell(i,j)%lai * lnd%cell(i,j)%f_veg * area(i,j)
                    lai_den = lai_den + lnd%cell(i,j)%f_veg * area(i,j)

                    ! total snow volume [10^12 m^3] — w_snow [kg/m2] / rho_w
                    ! * area / 1e12
                    snow_vol_step = snow_vol_step + lnd%cell(i,j)%w_snow * area(i,j) * 1.e-15_wp
                end if
            end do
        end do
        !$omp end parallel do

        ! roll thread-reduced locals into the module-level year accumulators
        acc_prc    = acc_prc    + l_prc
        acc_evp    = acc_evp    + l_evp
        acc_esur   = acc_esur   + l_esur
        acc_trans  = acc_trans  + l_trans
        acc_run    = acc_run    + l_run
        acc_runsur = acc_runsur + l_runsur
        acc_drain  = acc_drain  + l_drain
        acc_gpp    = acc_gpp    + l_gpp
        acc_npp    = acc_npp    + l_npp
        acc_ch4    = acc_ch4    + l_ch4
        acc_n2o    = acc_n2o    + l_n2o
        end block

        ! guarded means: divide by land-weighted area at the moment (per
        ! step), accumulate a running mean via mon_avg/nmon_year weighting
        if (area_land .gt. 0._wp) then
            acc_t2m_num = acc_t2m_num + t_land   * mon_avg / real(nmon_year, wp)
            acc_t2m_den = acc_t2m_den + area_land* mon_avg / real(nmon_year, wp)
        end if
        if (h_den .gt. 0._wp) then
            acc_hsnow_num = acc_hsnow_num + h_num * mon_avg / real(nmon_year, wp)
            acc_hsnow_den = acc_hsnow_den + h_den * mon_avg / real(nmon_year, wp)
        end if
        if (lai_den .gt. 0._wp) then
            acc_lai_num = acc_lai_num + lai_num * mon_avg / real(nmon_year, wp)
            acc_lai_den = acc_lai_den + lai_den * mon_avg / real(nmon_year, wp)
        end if
        acc_snowvol = acc_snowvol + snow_vol_step * mon_avg / real(nmon_year, wp)

        ! -------------------------------------------------------------
        ! End-of-year: finalise ts entry + build annual 2D fields, then
        ! write to disk when the ts cadence hits (ts always) and when
        ! time_out_lndvc hits (surf/hires).
        ! -------------------------------------------------------------
        if (time_eoy_lndvc) then

            ! Build the annual ts entry for this year.
            call build_ann_ts(lnd, ann_ts(y))

            ! Annual surface field = mean over the 12 months.
            call aggregate_mon_to_ann(mon_su, ann_su, lnd)

            ! Apply land mask (missing_value where mask_lnd==0) to the
            ! annual surface arrays for a friendlier NC file.
            call apply_mask_surf(ann_su, lnd)

            ! Paint annual surf onto the hires grid.
            if (write_hires) call paint_hires(ann_su, hires_su)

            ! ts writing: match lnd_out — writes past y years at once with
            ! the running index. Every year at eoy.
            if (time_out_ts_clim) then
                call ts_nc_write(trim(out_dir)//"/lndvc_ts.nc", &
                                 ann_ts(1:y), year_clim - y + 1, y)
            end if

            ! surf / hires writing: only at time_out_lndvc.
            if (time_out_lndvc) then

                nout = nout + 1

                if (write_surf) then
                    call surf_nc_write_all(trim(out_dir)//"/lndvc_surf.nc")
                end if

                if (write_hires) then
                    call hires_nc_write(trim(out_dir)//"/lndvc_hires.nc")
                end if

            end if

            ! Console line, mirror the reference lnd summary (single line
            ! per decade with a header refresh).
            if (mod(year, 10) .eq. 1) then
                print '(a7,a9,10a8)', 'lndvc', 'year', 'Cflx', 'NPP', 'CH4e', &
                    'N2Oe', 'prc', 'evp', 'run', 'Awet', 'Cveg', 'Csoil'
            end if
            print '(a7,i9,F8.3,3F8.2,3F8.1,3F8.1)', 'lndvc', year_now, &
                ann_ts(y)%Cflx_atm_lnd, ann_ts(y)%npp, ann_ts(y)%ch4_emis, &
                ann_ts(y)%n2o_emis, ann_ts(y)%prc, ann_ts(y)%evp, &
                ann_ts(y)%run, ann_ts(y)%awet, ann_ts(y)%veg_c, ann_ts(y)%soil_c

        end if

        return

    end subroutine lndvc_diag


    !======================================================================
    ! ts + surf + hires nc create (dimension definitions)
    !======================================================================
    subroutine ts_nc(fnm)

        implicit none

        character(len=*) :: fnm
        real(wp) :: empty_time(0)

        call nc_create(fnm)
        call nc_write_dim(fnm, dim_time, x=empty_time, axis="t", &
             units="years BP", unlimited=.TRUE.)
        call nc_write_dim(fnm, dim_lat, x=1, axis="y", units="1")
        call nc_write_dim(fnm, dim_lon, x=1, axis="x", units="1")

        return

    end subroutine ts_nc


    subroutine surf_nc(fnm)

        implicit none

        character(len=*) :: fnm
        integer :: ncid
        real(wp) :: empty_time(0)

        call nc_create(fnm)
        call nc_open(fnm, ncid)
        call nc_write_dim(fnm, dim_time, x=empty_time, axis="t", &
             units="years BP", unlimited=.TRUE., ncid=ncid)
        call nc_write_dim(fnm, dim_month, x=1._wp, dx=1._wp, nx=13, axis="e", &
             units="months", ncid=ncid)
        call nc_write_dim(fnm, dim_lat, x=lat, axis="y", &
             units="degrees_north", ncid=ncid)
        call nc_write_dim(fnm, dim_lon, x=lon, axis="x", &
             units="degrees_east", ncid=ncid)
        call nc_close(ncid)

        return

    end subroutine surf_nc


    subroutine hires_nc(fnm, geo)

        implicit none

        character(len=*) :: fnm
        type(geo_class), intent(in) :: geo

        integer :: ncid
        real(wp) :: empty_time(0)

        call nc_create(fnm)
        call nc_open(fnm, ncid)
        call nc_write_dim(fnm, dim_time, x=empty_time, axis="t", &
             units="years BP", unlimited=.TRUE., ncid=ncid)
        call nc_write_dim(fnm, dim_lat, x=geo%hires%grid%G%y, axis="y", &
             units="degrees_north", ncid=ncid)
        call nc_write_dim(fnm, dim_lon, x=geo%hires%grid%G%x, axis="x", &
             units="degrees_east", ncid=ncid)
        call nc_close(ncid)

        return

    end subroutine hires_nc


    !======================================================================
    ! ts write — mirror lnd_out::ts_nc_write shape (past y years).
    !======================================================================
    subroutine ts_nc_write(fnm, vars, ndat, y)

        implicit none

        type(ts_out), intent(in) :: vars(:)
        character(len=*)         :: fnm
        integer, intent(in)      :: ndat, y

        integer :: i, ncid

        call nc_open(fnm, ncid)
        call nc_write(fnm, "time", &
             real([(i, i = (year_now - (y - 1)*n_accel), year_now, n_accel)], wp), &
             dim1=dim_time, start=[ndat], count=[y], ncid=ncid)

        call nc_write(fnm, "temp",         sngl(vars%t2m - T0), &
             dim1=dim_time, start=[ndat], count=[y], &
             long_name="land 2 m temperature (area-weighted over veg)", &
             units="degC", missing_value=missing_value, ncid=ncid)

        call nc_write(fnm, "Aland",   sngl(vars%land), dim1=dim_time, &
             start=[ndat], count=[y], long_name="global land area", &
             units="mln km^2", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "Aveg",    sngl(vars%veg), dim1=dim_time, &
             start=[ndat], count=[y], long_name="global vegetated area", &
             units="mln km^2", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "Aice",    sngl(vars%ice), dim1=dim_time, &
             start=[ndat], count=[y], long_name="global ice area", &
             units="mln km^2", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "Alake",   sngl(vars%lake), dim1=dim_time, &
             start=[ndat], count=[y], long_name="global lake area", &
             units="mln km^2", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "Apeat",   sngl(vars%peat), dim1=dim_time, &
             start=[ndat], count=[y], long_name="global peatland area", &
             units="mln km^2", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "Awet",    sngl(vars%awet), dim1=dim_time, &
             start=[ndat], count=[y], long_name="global wetland area", &
             units="mln km^2", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "Aperm",   sngl(vars%aperm), dim1=dim_time, &
             start=[ndat], count=[y], long_name="permafrost area", &
             units="mln km^2", missing_value=missing_value, ncid=ncid)

        call nc_write(fnm, "prc",     sngl(vars%prc), dim1=dim_time, &
             start=[ndat], count=[y], long_name="precipitation", &
             units="10^15 kg/yr", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "evp",     sngl(vars%evp), dim1=dim_time, &
             start=[ndat], count=[y], long_name="evapotranspiration", &
             units="10^15 kg/yr", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "esur",    sngl(vars%esur), dim1=dim_time, &
             start=[ndat], count=[y], long_name="surface evaporation", &
             units="10^15 kg/yr", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "trans",   sngl(vars%trans), dim1=dim_time, &
             start=[ndat], count=[y], long_name="transpiration", &
             units="10^15 kg/yr", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "runoff",  sngl(vars%run), dim1=dim_time, &
             start=[ndat], count=[y], long_name="total runoff", &
             units="10^15 kg/yr", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "runsur",  sngl(vars%runsur), dim1=dim_time, &
             start=[ndat], count=[y], long_name="surface runoff", &
             units="10^15 kg/yr", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "drain",   sngl(vars%drain), dim1=dim_time, &
             start=[ndat], count=[y], long_name="drainage out of soil column", &
             units="10^15 kg/yr", missing_value=missing_value, ncid=ncid)

        call nc_write(fnm, "gpp",     sngl(vars%gpp), dim1=dim_time, &
             start=[ndat], count=[y], long_name="gross primary productivity", &
             units="PgC/yr", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "npp",     sngl(vars%npp), dim1=dim_time, &
             start=[ndat], count=[y], long_name="net primary productivity", &
             units="PgC/yr", missing_value=missing_value, ncid=ncid)

        call nc_write(fnm, "Cveg",    sngl(vars%veg_c), dim1=dim_time, &
             start=[ndat], count=[y], long_name="vegetation carbon", &
             units="PgC", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "Csoil",   sngl(vars%soil_c), dim1=dim_time, &
             start=[ndat], count=[y], long_name="total soil carbon", &
             units="PgC", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "Cpeat",   sngl(vars%peat_c), dim1=dim_time, &
             start=[ndat], count=[y], long_name="peatland carbon", &
             units="PgC", missing_value=missing_value, ncid=ncid)

        call nc_write(fnm, "Cflx_atm_lnd", sngl(vars%Cflx_atm_lnd), &
             dim1=dim_time, start=[ndat], count=[y], &
             long_name="net atmosphere-land carbon flux", &
             units="PgC/yr", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "ch4",     sngl(vars%ch4_emis), dim1=dim_time, &
             start=[ndat], count=[y], long_name="total CH4 emissions", &
             units="TgCH4/yr", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "n2o",     sngl(vars%n2o_emis), dim1=dim_time, &
             start=[ndat], count=[y], long_name="total N2O-N emissions", &
             units="TgN2O-N/yr", missing_value=missing_value, ncid=ncid)

        call nc_write(fnm, "Vsnow",   sngl(vars%snow_vol), dim1=dim_time, &
             start=[ndat], count=[y], long_name="total snow volume", &
             units="10^12 m^3", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "hsnow",   sngl(vars%h_snow), dim1=dim_time, &
             start=[ndat], count=[y], long_name="area-mean snow depth", &
             units="m", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "lai",     sngl(vars%lai), dim1=dim_time, &
             start=[ndat], count=[y], long_name="area-mean LAI", &
             units="/", missing_value=missing_value, ncid=ncid)

        call nc_close(ncid)

        return

    end subroutine ts_nc_write


    !======================================================================
    ! surf write — one output frame writes 12 monthly slots + annual (=13).
    !======================================================================
    subroutine surf_nc_write_all(fnm)

        implicit none

        character(len=*) :: fnm

        integer :: k, ncid

        call nc_open(fnm, ncid)
        call nc_write(fnm, dim_time, real(year_now, wp), dim1=dim_time, &
             start=[nout], count=[1], ncid=ncid)

        do k = 1, nmon_year
            call surf_nc_write(fnm, ncid, mon_su(k), k)
        end do
        call surf_nc_write(fnm, ncid, ann_su, nmon_year + 1)

        call nc_close(ncid)

        return

    end subroutine surf_nc_write_all


    subroutine surf_nc_write(fnm, ncid, vars, ndat)

        implicit none

        character(len=*)             :: fnm
        integer,      intent(in)     :: ncid, ndat
        type(surf_out), intent(in)   :: vars

        ! seasonally-invariant class fractions: only written on the annual
        ! slot (matches lnd_out convention for fsurf).
        if (ndat .eq. nmon_year + 1) then
            call nc_write(fnm, "f_land", sngl(vars%f_land), &
                 dims=[dim_lon,dim_lat,dim_time], start=[1,1,nout], count=[ni,nj,1], &
                 long_name="land fraction", units="/", missing_value=missing_value, ncid=ncid)
            call nc_write(fnm, "f_veg",  sngl(vars%f_veg), &
                 dims=[dim_lon,dim_lat,dim_time], start=[1,1,nout], count=[ni,nj,1], &
                 long_name="vegetated fraction", units="/", missing_value=missing_value, ncid=ncid)
            call nc_write(fnm, "f_ice",  sngl(vars%f_ice), &
                 dims=[dim_lon,dim_lat,dim_time], start=[1,1,nout], count=[ni,nj,1], &
                 long_name="ice-sheet fraction", units="/", missing_value=missing_value, ncid=ncid)
            call nc_write(fnm, "f_lake", sngl(vars%f_lake), &
                 dims=[dim_lon,dim_lat,dim_time], start=[1,1,nout], count=[ni,nj,1], &
                 long_name="lake fraction", units="/", missing_value=missing_value, ncid=ncid)
            call nc_write(fnm, "f_peat", sngl(vars%f_peat), &
                 dims=[dim_lon,dim_lat,dim_time], start=[1,1,nout], count=[ni,nj,1], &
                 long_name="peatland fraction", units="/", missing_value=missing_value, ncid=ncid)
            call nc_write(fnm, "gdd5",   sngl(vars%gdd5), &
                 dims=[dim_lon,dim_lat,dim_time], start=[1,1,nout], count=[ni,nj,1], &
                 long_name="growing degree days above 5 degC", units="K", &
                 missing_value=missing_value, ncid=ncid)
            call nc_write(fnm, "t2m_min_mon", sngl(vars%t2m_min_mon), &
                 dims=[dim_lon,dim_lat,dim_time], start=[1,1,nout], count=[ni,nj,1], &
                 long_name="minimum monthly 2 m air temperature", units="K", &
                 missing_value=missing_value, ncid=ncid)
            call nc_write(fnm, "alt",    sngl(vars%alt), &
                 dims=[dim_lon,dim_lat,dim_time], start=[1,1,nout], count=[ni,nj,1], &
                 long_name="active layer thickness", units="m", &
                 missing_value=missing_value, ncid=ncid)
        end if

        ! energy / near-surface state (monthly climatology)
        call nc_write(fnm, "t_skin",     sngl(vars%t_skin), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="skin temperature", units="K", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "albedo",     sngl(vars%albedo), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="surface albedo", units="/", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "t2m",        sngl(vars%t2m), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="2 m air temperature", units="K", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "q2m",        sngl(vars%q2m), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="2 m specific humidity", units="kg/kg", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "swnet",      sngl(vars%swnet), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="net shortwave radiation", units="W/m2", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "lwnet",      sngl(vars%lwnet), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="net longwave radiation", units="W/m2", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "flx_sh",     sngl(vars%flx_sh), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="sensible heat flux", units="W/m2", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "flx_lh",     sngl(vars%flx_lh), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="latent heat flux", units="W/m2", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "flx_g",      sngl(vars%flx_g), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="ground heat flux", units="W/m2", &
             missing_value=missing_value, ncid=ncid)

        ! snow
        call nc_write(fnm, "h_snow",     sngl(vars%h_snow), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="snow depth", units="m", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "w_snow",     sngl(vars%w_snow), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="snow water equivalent", units="kg/m2", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "w_snow_max", sngl(vars%w_snow_max), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="seasonal maximum SWE", units="kg/m2", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "snowmelt",   sngl(vars%snowmelt), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="snowmelt rate", units="kg/m2/s", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "icemelt",    sngl(vars%icemelt), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="icemelt rate", units="kg/m2/s", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "f_snow",     sngl(vars%f_snow), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="snow-cover fraction", units="/", &
             missing_value=missing_value, ncid=ncid)

        ! precip
        call nc_write(fnm, "rain",       sngl(vars%rain), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="rain flux", units="kg/m2/s", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "snow_flx",   sngl(vars%snow_flx), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="snowfall flux", units="kg/m2/s", &
             missing_value=missing_value, ncid=ncid)

        ! water balance
        call nc_write(fnm, "trans",      sngl(vars%transpiration), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="transpiration", units="kg/m2/s", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "esur",       sngl(vars%evap_surface), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="surface evaporation", units="kg/m2/s", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "ecan",       sngl(vars%evap_can), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="canopy evaporation", units="kg/m2/s", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "et",         sngl(vars%et), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="total evapotranspiration", units="kg/m2/s", &
             missing_value=missing_value, ncid=ncid)

        call nc_write(fnm, "runoff",     sngl(vars%runoff), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="total runoff", units="kg/m2/s", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "runsur",     sngl(vars%runoff_sur), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="surface runoff", units="kg/m2/s", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "runsub",     sngl(vars%runoff_gw), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="groundwater baseflow", units="kg/m2/s", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "drain",      sngl(vars%drainage), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="drainage out of soil column", units="kg/m2/s", &
             missing_value=missing_value, ncid=ncid)

        call nc_write(fnm, "wtab",       sngl(vars%w_table), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="water table depth", units="m", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "wtab_peat",  sngl(vars%w_table_peat), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="water table depth (peat)", units="m", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "fwet",       sngl(vars%f_wet), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="saturated fraction", units="/", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "fwetland",   sngl(vars%f_wetland), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="TOPMODEL wetland fraction", units="/", &
             missing_value=missing_value, ncid=ncid)

        ! vegetation + soil
        call nc_write(fnm, "lai",        sngl(vars%lai), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="leaf area index", units="/", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "sai",        sngl(vars%sai), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="stem area index", units="/", &
             missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "gpp",        sngl(vars%gpp), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="gross primary productivity", &
             units="kgC/m2/s", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "npp",        sngl(vars%npp), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="net primary productivity", &
             units="kgC/m2/s", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "veg_c",      sngl(vars%veg_c), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="vegetation carbon", &
             units="kgC/m2", missing_value=missing_value, ncid=ncid)

        call nc_write(fnm, "t_soil_top", sngl(vars%t_soil_top), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="top-layer soil temperature", &
             units="K", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "theta_w_top",sngl(vars%theta_w_top), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="top-layer volumetric soil moisture", &
             units="/", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "soil_c",     sngl(vars%soil_c), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="soil carbon (mineral+peat+lake+shelf+ice)", &
             units="kgC/m2", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "soil_resp",  sngl(vars%soil_resp), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="soil respiration", &
             units="kgC/m2/s", missing_value=missing_value, ncid=ncid)

        ! surface mass balance (ice + snow)
        call nc_write(fnm, "smb",        sngl(vars%smb), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="surface mass balance", &
             units="kg/m2/s", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "melt",       sngl(vars%melt), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="ice + snow melt", &
             units="kg/m2/s", missing_value=missing_value, ncid=ncid)

        ! cell-scale carbon / trace-gas fluxes
        call nc_write(fnm, "Cflx_atm_lnd", sngl(vars%Cflx_atm_lnd), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="net atmosphere-land carbon flux", &
             units="kgC/m2/s", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "ch4",        sngl(vars%ch4_emis), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="CH4 emissions", &
             units="kgC/m2/s", missing_value=missing_value, ncid=ncid)
        call nc_write(fnm, "n2o",        sngl(vars%n2o_emis), &
             dims=[dim_lon,dim_lat,dim_month,dim_time], start=[1,1,ndat,nout], &
             count=[ni,nj,1,1], long_name="N2O emissions", &
             units="kgN/m2/s", missing_value=missing_value, ncid=ncid)

        return

    end subroutine surf_nc_write


    !======================================================================
    ! hires write — one annual painted frame per output.
    !======================================================================
    subroutine hires_nc_write(fnm)

        implicit none

        character(len=*) :: fnm

        integer :: ncid

        call nc_open(fnm, ncid)
        call nc_write(fnm, dim_time, real(year_now, wp), dim1=dim_time, &
             start=[nout], count=[1], ncid=ncid)

        if (hires_write_h_snow) then
            call nc_write(fnm, "h_snow", sngl(hires_su%h_snow), &
                 dims=[dim_lon,dim_lat,dim_time], start=[1,1,nout], &
                 count=[nx_hi,ny_hi,1], long_name="snow depth", units="m", &
                 missing_value=missing_value, ncid=ncid)
        end if
        if (hires_write_w_snow) then
            call nc_write(fnm, "w_snow", sngl(hires_su%w_snow), &
                 dims=[dim_lon,dim_lat,dim_time], start=[1,1,nout], &
                 count=[nx_hi,ny_hi,1], long_name="snow water equivalent", &
                 units="kg/m2", missing_value=missing_value, ncid=ncid)
        end if
        if (hires_write_t_skin) then
            call nc_write(fnm, "t_skin", sngl(hires_su%t_skin), &
                 dims=[dim_lon,dim_lat,dim_time], start=[1,1,nout], &
                 count=[nx_hi,ny_hi,1], long_name="skin temperature", units="K", &
                 missing_value=missing_value, ncid=ncid)
        end if
        if (hires_write_smb) then
            call nc_write(fnm, "smb", sngl(hires_su%smb), &
                 dims=[dim_lon,dim_lat,dim_time], start=[1,1,nout], &
                 count=[nx_hi,ny_hi,1], long_name="surface mass balance", &
                 units="kg/m2/s", missing_value=missing_value, ncid=ncid)
        end if
        if (hires_write_melt) then
            call nc_write(fnm, "melt", sngl(hires_su%melt), &
                 dims=[dim_lon,dim_lat,dim_time], start=[1,1,nout], &
                 count=[nx_hi,ny_hi,1], long_name="ice + snow melt", &
                 units="kg/m2/s", missing_value=missing_value, ncid=ncid)
        end if
        if (hires_write_f_ice) then
            call nc_write(fnm, "f_ice", sngl(hires_su%f_ice), &
                 dims=[dim_lon,dim_lat,dim_time], start=[1,1,nout], &
                 count=[nx_hi,ny_hi,1], long_name="ice-sheet fraction (painted)", &
                 units="/", missing_value=missing_value, ncid=ncid)
        end if
        if (hires_write_veg_c) then
            call nc_write(fnm, "veg_c", sngl(hires_su%veg_c), &
                 dims=[dim_lon,dim_lat,dim_time], start=[1,1,nout], &
                 count=[nx_hi,ny_hi,1], long_name="vegetation carbon (painted)", &
                 units="kgC/m2", missing_value=missing_value, ncid=ncid)
        end if
        if (hires_write_soil_c) then
            call nc_write(fnm, "soil_c", sngl(hires_su%soil_c), &
                 dims=[dim_lon,dim_lat,dim_time], start=[1,1,nout], &
                 count=[nx_hi,ny_hi,1], long_name="soil carbon (painted)", &
                 units="kgC/m2", missing_value=missing_value, ncid=ncid)
        end if

        call nc_close(ncid)

        return

    end subroutine hires_nc_write


    !======================================================================
    ! surface allocation + zeroing + accumulation helpers
    !======================================================================
    subroutine alloc_surf(v)

        implicit none

        type(surf_out), intent(inout) :: v

        allocate(v%f_land(ni,nj), v%f_veg(ni,nj), v%f_ice(ni,nj), &
                 v%f_lake(ni,nj), v%f_peat(ni,nj))
        allocate(v%t_skin(ni,nj), v%albedo(ni,nj), v%t2m(ni,nj), v%q2m(ni,nj))
        allocate(v%swnet(ni,nj), v%lwnet(ni,nj), v%flx_sh(ni,nj), &
                 v%flx_lh(ni,nj), v%flx_g(ni,nj))
        allocate(v%h_snow(ni,nj), v%w_snow(ni,nj), v%w_snow_max(ni,nj), &
                 v%snowmelt(ni,nj), v%icemelt(ni,nj), v%f_snow(ni,nj))
        allocate(v%rain(ni,nj), v%snow_flx(ni,nj))
        allocate(v%transpiration(ni,nj), v%evap_surface(ni,nj), &
                 v%evap_can(ni,nj), v%et(ni,nj))
        allocate(v%runoff(ni,nj), v%runoff_sur(ni,nj), v%runoff_gw(ni,nj), &
                 v%drainage(ni,nj))
        allocate(v%w_table(ni,nj), v%w_table_peat(ni,nj))
        allocate(v%f_wet(ni,nj), v%f_wetland(ni,nj))
        allocate(v%lai(ni,nj), v%sai(ni,nj), v%gpp(ni,nj), v%npp(ni,nj), &
                 v%veg_c(ni,nj))
        allocate(v%gdd5(ni,nj), v%t2m_min_mon(ni,nj))
        allocate(v%t_soil_top(ni,nj), v%theta_w_top(ni,nj), v%alt(ni,nj))
        allocate(v%soil_c(ni,nj), v%soil_resp(ni,nj))
        allocate(v%smb(ni,nj), v%melt(ni,nj))
        allocate(v%Cflx_atm_lnd(ni,nj), v%ch4_emis(ni,nj), v%n2o_emis(ni,nj))

        return

    end subroutine alloc_surf


    subroutine zero_surf(v)

        implicit none

        type(surf_out), intent(inout) :: v

        v%f_land = 0._wp; v%f_veg = 0._wp; v%f_ice = 0._wp
        v%f_lake = 0._wp; v%f_peat = 0._wp
        v%t_skin = 0._wp; v%albedo = 0._wp; v%t2m = 0._wp; v%q2m = 0._wp
        v%swnet = 0._wp; v%lwnet = 0._wp
        v%flx_sh = 0._wp; v%flx_lh = 0._wp; v%flx_g = 0._wp
        v%h_snow = 0._wp; v%w_snow = 0._wp; v%w_snow_max = 0._wp
        v%snowmelt = 0._wp; v%icemelt = 0._wp; v%f_snow = 0._wp
        v%rain = 0._wp; v%snow_flx = 0._wp
        v%transpiration = 0._wp; v%evap_surface = 0._wp
        v%evap_can = 0._wp; v%et = 0._wp
        v%runoff = 0._wp; v%runoff_sur = 0._wp
        v%runoff_gw = 0._wp; v%drainage = 0._wp
        v%w_table = 0._wp; v%w_table_peat = 0._wp
        v%f_wet = 0._wp; v%f_wetland = 0._wp
        v%lai = 0._wp; v%sai = 0._wp; v%gpp = 0._wp; v%npp = 0._wp
        v%veg_c = 0._wp
        v%gdd5 = 0._wp; v%t2m_min_mon = 0._wp
        v%t_soil_top = 0._wp; v%theta_w_top = 0._wp; v%alt = 0._wp
        v%soil_c = 0._wp; v%soil_resp = 0._wp
        v%smb = 0._wp; v%melt = 0._wp
        v%Cflx_atm_lnd = 0._wp; v%ch4_emis = 0._wp; v%n2o_emis = 0._wp

        return

    end subroutine zero_surf


    subroutine accum_surf_cell(v, c, i, j, w)

        implicit none

        type(surf_out), intent(inout) :: v
        type(vc_cell_t), intent(in)   :: c
        integer,  intent(in)          :: i, j
        real(wp), intent(in)          :: w

        v%f_land(i,j)   = v%f_land(i,j)   + w * c%f_land
        v%f_veg(i,j)    = v%f_veg(i,j)    + w * c%f_veg
        v%f_ice(i,j)    = v%f_ice(i,j)    + w * c%f_ice
        v%f_lake(i,j)   = v%f_lake(i,j)   + w * c%f_lake
        v%f_peat(i,j)   = v%f_peat(i,j)   + w * c%f_peat

        v%t_skin(i,j)   = v%t_skin(i,j)   + w * c%t_skin
        v%albedo(i,j)   = v%albedo(i,j)   + w * c%albedo
        v%t2m(i,j)      = v%t2m(i,j)      + w * c%t2m
        v%q2m(i,j)      = v%q2m(i,j)      + w * c%q2m
        v%swnet(i,j)    = v%swnet(i,j)    + w * c%swnet
        v%lwnet(i,j)    = v%lwnet(i,j)    + w * c%lwnet
        v%flx_sh(i,j)   = v%flx_sh(i,j)   + w * c%flx_sh
        v%flx_lh(i,j)   = v%flx_lh(i,j)   + w * c%flx_lh
        v%flx_g(i,j)    = v%flx_g(i,j)    + w * c%flx_g

        v%h_snow(i,j)   = v%h_snow(i,j)     + w * c%h_snow
        v%w_snow(i,j)   = v%w_snow(i,j)     + w * c%w_snow
        v%w_snow_max(i,j)=v%w_snow_max(i,j) + w * c%w_snow_max
        v%snowmelt(i,j) = v%snowmelt(i,j)   + w * c%snowmelt
        v%icemelt(i,j)  = v%icemelt(i,j)    + w * c%icemelt
        v%f_snow(i,j)   = v%f_snow(i,j)     + w * c%f_snow

        v%rain(i,j)     = v%rain(i,j)     + w * c%rain
        v%snow_flx(i,j) = v%snow_flx(i,j) + w * c%snow_flx

        v%transpiration(i,j) = v%transpiration(i,j) + w * c%transpiration
        v%evap_surface(i,j)  = v%evap_surface(i,j)  + w * c%evap_surface
        v%evap_can(i,j)      = v%evap_can(i,j)      + w * c%evap_can
        v%et(i,j)            = v%et(i,j)            + w * c%et

        v%runoff(i,j)     = v%runoff(i,j)     + w * c%runoff
        v%runoff_sur(i,j) = v%runoff_sur(i,j) + w * c%runoff_sur
        v%runoff_gw(i,j)  = v%runoff_gw(i,j)  + w * c%runoff_gw
        v%drainage(i,j)   = v%drainage(i,j)   + w * c%drainage

        v%w_table(i,j)      = v%w_table(i,j)      + w * c%w_table
        v%w_table_peat(i,j) = v%w_table_peat(i,j) + w * c%w_table_peat
        v%f_wet(i,j)        = v%f_wet(i,j)        + w * c%f_wet
        v%f_wetland(i,j)    = v%f_wetland(i,j)    + w * c%f_wetland

        v%lai(i,j)     = v%lai(i,j)     + w * c%lai
        v%sai(i,j)     = v%sai(i,j)     + w * c%sai
        v%gpp(i,j)     = v%gpp(i,j)     + w * c%gpp
        v%npp(i,j)     = v%npp(i,j)     + w * c%npp
        v%veg_c(i,j)   = v%veg_c(i,j)   + w * c%veg_c
        v%gdd5(i,j)    = v%gdd5(i,j)    + w * c%gdd5
        v%t2m_min_mon(i,j) = v%t2m_min_mon(i,j) + w * c%t2m_min_mon

        v%t_soil_top(i,j)  = v%t_soil_top(i,j)  + w * c%t_soil_top
        v%theta_w_top(i,j) = v%theta_w_top(i,j) + w * c%theta_w_top
        v%alt(i,j)         = v%alt(i,j)         + w * c%alt
        v%soil_c(i,j)      = v%soil_c(i,j)      + w * c%soil_c
        v%soil_resp(i,j)   = v%soil_resp(i,j)   + w * c%soil_resp

        v%smb(i,j)  = v%smb(i,j)  + w * c%smb
        v%melt(i,j) = v%melt(i,j) + w * c%melt

        v%Cflx_atm_lnd(i,j) = v%Cflx_atm_lnd(i,j) + w * c%Cflx_atm_lnd
        v%ch4_emis(i,j)     = v%ch4_emis(i,j)     + w * c%ch4_emis
        v%n2o_emis(i,j)     = v%n2o_emis(i,j)     + w * c%n2o_emis

        return

    end subroutine accum_surf_cell


    subroutine aggregate_mon_to_ann(mvars, avars, lnd)

        implicit none

        type(surf_out), intent(in)    :: mvars(:)
        type(surf_out), intent(inout) :: avars
        type(lndvc_class), intent(in) :: lnd

        integer :: k, i, j
        real(wp) :: div

        div = real(size(mvars), wp)

        call zero_surf(avars)

        do k = 1, size(mvars)
            do j = 1, nj
                do i = 1, ni
                    avars%f_land(i,j)   = avars%f_land(i,j)   + mvars(k)%f_land(i,j)   / div
                    avars%f_veg(i,j)    = avars%f_veg(i,j)    + mvars(k)%f_veg(i,j)    / div
                    avars%f_ice(i,j)    = avars%f_ice(i,j)    + mvars(k)%f_ice(i,j)    / div
                    avars%f_lake(i,j)   = avars%f_lake(i,j)   + mvars(k)%f_lake(i,j)   / div
                    avars%f_peat(i,j)   = avars%f_peat(i,j)   + mvars(k)%f_peat(i,j)   / div

                    avars%t_skin(i,j)   = avars%t_skin(i,j)   + mvars(k)%t_skin(i,j)   / div
                    avars%albedo(i,j)   = avars%albedo(i,j)   + mvars(k)%albedo(i,j)   / div
                    avars%t2m(i,j)      = avars%t2m(i,j)      + mvars(k)%t2m(i,j)      / div
                    avars%q2m(i,j)      = avars%q2m(i,j)      + mvars(k)%q2m(i,j)      / div
                    avars%swnet(i,j)    = avars%swnet(i,j)    + mvars(k)%swnet(i,j)    / div
                    avars%lwnet(i,j)    = avars%lwnet(i,j)    + mvars(k)%lwnet(i,j)    / div
                    avars%flx_sh(i,j)   = avars%flx_sh(i,j)   + mvars(k)%flx_sh(i,j)   / div
                    avars%flx_lh(i,j)   = avars%flx_lh(i,j)   + mvars(k)%flx_lh(i,j)   / div
                    avars%flx_g(i,j)    = avars%flx_g(i,j)    + mvars(k)%flx_g(i,j)    / div

                    avars%h_snow(i,j)   = avars%h_snow(i,j)   + mvars(k)%h_snow(i,j)   / div
                    avars%w_snow(i,j)   = avars%w_snow(i,j)   + mvars(k)%w_snow(i,j)   / div
                    avars%w_snow_max(i,j)= max(avars%w_snow_max(i,j), mvars(k)%w_snow_max(i,j))
                    avars%snowmelt(i,j) = avars%snowmelt(i,j) + mvars(k)%snowmelt(i,j) / div
                    avars%icemelt(i,j)  = avars%icemelt(i,j)  + mvars(k)%icemelt(i,j)  / div
                    avars%f_snow(i,j)   = avars%f_snow(i,j)   + mvars(k)%f_snow(i,j)   / div

                    avars%rain(i,j)     = avars%rain(i,j)     + mvars(k)%rain(i,j)     / div
                    avars%snow_flx(i,j) = avars%snow_flx(i,j) + mvars(k)%snow_flx(i,j) / div

                    avars%transpiration(i,j) = avars%transpiration(i,j) + mvars(k)%transpiration(i,j) / div
                    avars%evap_surface(i,j)  = avars%evap_surface(i,j)  + mvars(k)%evap_surface(i,j)  / div
                    avars%evap_can(i,j)      = avars%evap_can(i,j)      + mvars(k)%evap_can(i,j)      / div
                    avars%et(i,j)            = avars%et(i,j)            + mvars(k)%et(i,j)            / div

                    avars%runoff(i,j)     = avars%runoff(i,j)     + mvars(k)%runoff(i,j)     / div
                    avars%runoff_sur(i,j) = avars%runoff_sur(i,j) + mvars(k)%runoff_sur(i,j) / div
                    avars%runoff_gw(i,j)  = avars%runoff_gw(i,j)  + mvars(k)%runoff_gw(i,j)  / div
                    avars%drainage(i,j)   = avars%drainage(i,j)   + mvars(k)%drainage(i,j)   / div

                    avars%w_table(i,j)      = avars%w_table(i,j)      + mvars(k)%w_table(i,j)      / div
                    avars%w_table_peat(i,j) = avars%w_table_peat(i,j) + mvars(k)%w_table_peat(i,j) / div
                    avars%f_wet(i,j)        = avars%f_wet(i,j)        + mvars(k)%f_wet(i,j)        / div
                    avars%f_wetland(i,j)    = avars%f_wetland(i,j)    + mvars(k)%f_wetland(i,j)    / div

                    avars%lai(i,j)   = avars%lai(i,j)   + mvars(k)%lai(i,j)   / div
                    avars%sai(i,j)   = avars%sai(i,j)   + mvars(k)%sai(i,j)   / div
                    avars%gpp(i,j)   = avars%gpp(i,j)   + mvars(k)%gpp(i,j)   / div
                    avars%npp(i,j)   = avars%npp(i,j)   + mvars(k)%npp(i,j)   / div
                    avars%veg_c(i,j) = avars%veg_c(i,j) + mvars(k)%veg_c(i,j) / div

                    ! phenology drivers: take annual max of the monthly
                    ! climatology instead of a mean (matches lnd_out).
                    avars%gdd5(i,j)        = max(avars%gdd5(i,j), mvars(k)%gdd5(i,j))
                    avars%t2m_min_mon(i,j) = avars%t2m_min_mon(i,j) + mvars(k)%t2m_min_mon(i,j) / div

                    avars%t_soil_top(i,j)  = avars%t_soil_top(i,j)  + mvars(k)%t_soil_top(i,j)  / div
                    avars%theta_w_top(i,j) = avars%theta_w_top(i,j) + mvars(k)%theta_w_top(i,j) / div
                    avars%alt(i,j)         = max(avars%alt(i,j), mvars(k)%alt(i,j))
                    avars%soil_c(i,j)      = avars%soil_c(i,j)      + mvars(k)%soil_c(i,j)      / div
                    avars%soil_resp(i,j)   = avars%soil_resp(i,j)   + mvars(k)%soil_resp(i,j)   / div

                    avars%smb(i,j)  = avars%smb(i,j)  + mvars(k)%smb(i,j)  / div
                    avars%melt(i,j) = avars%melt(i,j) + mvars(k)%melt(i,j) / div

                    avars%Cflx_atm_lnd(i,j) = avars%Cflx_atm_lnd(i,j) + mvars(k)%Cflx_atm_lnd(i,j) / div
                    avars%ch4_emis(i,j)     = avars%ch4_emis(i,j)     + mvars(k)%ch4_emis(i,j)     / div
                    avars%n2o_emis(i,j)     = avars%n2o_emis(i,j)     + mvars(k)%n2o_emis(i,j)     / div
                end do
            end do
        end do

        return

    end subroutine aggregate_mon_to_ann


    subroutine apply_mask_surf(v, lnd)

        implicit none

        type(surf_out), intent(inout) :: v
        type(lndvc_class), intent(in) :: lnd

        integer :: i, j

        do j = 1, nj
            do i = 1, ni
                if (lnd%cell(i,j)%mask_lnd .ne. 1) then
                    v%f_land(i,j)   = 0._wp
                    v%f_veg(i,j)    = 0._wp
                    v%f_ice(i,j)    = 0._wp
                    v%f_lake(i,j)   = 0._wp
                    v%f_peat(i,j)   = 0._wp
                    v%t_skin(i,j)   = real(missing_value, wp)
                    v%albedo(i,j)   = real(missing_value, wp)
                    v%t2m(i,j)      = real(missing_value, wp)
                    v%q2m(i,j)      = real(missing_value, wp)
                    v%swnet(i,j)    = real(missing_value, wp)
                    v%lwnet(i,j)    = real(missing_value, wp)
                    v%flx_sh(i,j)   = real(missing_value, wp)
                    v%flx_lh(i,j)   = real(missing_value, wp)
                    v%flx_g(i,j)    = real(missing_value, wp)
                    v%h_snow(i,j)   = real(missing_value, wp)
                    v%w_snow(i,j)   = real(missing_value, wp)
                    v%w_snow_max(i,j)=real(missing_value, wp)
                    v%snowmelt(i,j) = real(missing_value, wp)
                    v%icemelt(i,j)  = real(missing_value, wp)
                    v%f_snow(i,j)   = real(missing_value, wp)
                    v%rain(i,j)     = real(missing_value, wp)
                    v%snow_flx(i,j) = real(missing_value, wp)
                    v%transpiration(i,j)  = real(missing_value, wp)
                    v%evap_surface(i,j)   = real(missing_value, wp)
                    v%evap_can(i,j)       = real(missing_value, wp)
                    v%et(i,j)             = real(missing_value, wp)
                    v%runoff(i,j)         = real(missing_value, wp)
                    v%runoff_sur(i,j)     = real(missing_value, wp)
                    v%runoff_gw(i,j)      = real(missing_value, wp)
                    v%drainage(i,j)       = real(missing_value, wp)
                    v%w_table(i,j)        = real(missing_value, wp)
                    v%w_table_peat(i,j)   = real(missing_value, wp)
                    v%f_wet(i,j)          = real(missing_value, wp)
                    v%f_wetland(i,j)      = real(missing_value, wp)
                    v%lai(i,j)            = real(missing_value, wp)
                    v%sai(i,j)            = real(missing_value, wp)
                    v%gpp(i,j)            = real(missing_value, wp)
                    v%npp(i,j)            = real(missing_value, wp)
                    v%veg_c(i,j)          = real(missing_value, wp)
                    v%gdd5(i,j)           = real(missing_value, wp)
                    v%t2m_min_mon(i,j)    = real(missing_value, wp)
                    v%t_soil_top(i,j)     = real(missing_value, wp)
                    v%theta_w_top(i,j)    = real(missing_value, wp)
                    v%alt(i,j)            = real(missing_value, wp)
                    v%soil_c(i,j)         = real(missing_value, wp)
                    v%soil_resp(i,j)      = real(missing_value, wp)
                    v%smb(i,j)            = real(missing_value, wp)
                    v%melt(i,j)           = real(missing_value, wp)
                    v%Cflx_atm_lnd(i,j)   = real(missing_value, wp)
                    v%ch4_emis(i,j)       = real(missing_value, wp)
                    v%n2o_emis(i,j)       = real(missing_value, wp)
                end if
            end do
        end do

        return

    end subroutine apply_mask_surf


    !======================================================================
    ! ts (annual scalar) build helper.
    !======================================================================
    subroutine build_ann_ts(lnd, tsv)

        implicit none

        type(lndvc_class), intent(in) :: lnd
        type(ts_out), intent(inout)   :: tsv

        real(wp) :: fac_area
        integer  :: i, j

        fac_area = 1.e-12_wp   ! m2 -> mln km2

        ! Areas — snapshot at end of year.
        tsv%land = 0._wp; tsv%veg = 0._wp; tsv%ice = 0._wp
        tsv%lake = 0._wp; tsv%peat = 0._wp
        tsv%awet = 0._wp; tsv%aperm = 0._wp
        tsv%veg_c = 0._wp; tsv%soil_c = 0._wp; tsv%peat_c = 0._wp
        tsv%Cflx_atm_lnd = 0._wp

        do j = 1, nj
            do i = 1, ni
                if (lnd%cell(i,j)%mask_lnd .eq. 1) then
                    tsv%land = tsv%land + lnd%cell(i,j)%f_land * area(i,j) * fac_area
                    tsv%veg  = tsv%veg  + lnd%cell(i,j)%f_veg  * area(i,j) * fac_area
                    tsv%ice  = tsv%ice  + lnd%cell(i,j)%f_ice  * area(i,j) * fac_area
                    tsv%lake = tsv%lake + lnd%cell(i,j)%f_lake * area(i,j) * fac_area
                    tsv%peat = tsv%peat + lnd%cell(i,j)%f_peat * area(i,j) * fac_area

                    tsv%awet = tsv%awet &
                        + max(0._wp, lnd%cell(i,j)%f_wetland * lnd%cell(i,j)%f_veg &
                                     - lnd%cell(i,j)%f_peat) &
                        * area(i,j) * fac_area

                    if (lnd%cell(i,j)%alt .gt. -1._wp .and. lnd%cell(i,j)%f_veg .gt. 0._wp) then
                        tsv%aperm = tsv%aperm + lnd%cell(i,j)%f_veg * area(i,j) * fac_area
                    end if

                    ! carbon stocks [PgC] — kgC/m2 * m2 * 1e-12
                    tsv%veg_c  = tsv%veg_c  + lnd%cell(i,j)%veg_c  * lnd%cell(i,j)%f_veg  * area(i,j) * 1.e-12_wp
                    tsv%soil_c = tsv%soil_c + lnd%cell(i,j)%soil_c * lnd%cell(i,j)%f_land * area(i,j) * 1.e-12_wp
                    tsv%peat_c = tsv%peat_c + lnd%cell(i,j)%peat_c * lnd%cell(i,j)%f_veg  * area(i,j) * 1.e-12_wp

                    ! net atm-land C flux: kgC/m2/s * m2 * sec_year * 1e-12 [PgC/yr]
                    tsv%Cflx_atm_lnd = tsv%Cflx_atm_lnd &
                        + lnd%cell(i,j)%Cflx_atm_lnd * area(i,j) * sec_year * 1.e-12_wp
                end if
            end do
        end do

        ! per-year running totals
        tsv%prc    = acc_prc
        tsv%evp    = acc_evp
        tsv%esur   = acc_esur
        tsv%trans  = acc_trans
        tsv%run    = acc_run
        tsv%runsur = acc_runsur
        tsv%drain  = acc_drain
        tsv%gpp    = acc_gpp
        tsv%npp    = acc_npp
        tsv%ch4_emis = acc_ch4
        tsv%n2o_emis = acc_n2o

        ! per-year running means (guarded)
        if (acc_t2m_den .gt. 0._wp) then
            tsv%t2m = acc_t2m_num / acc_t2m_den
        else
            tsv%t2m = real(missing_value, wp)
        end if
        if (acc_hsnow_den .gt. 0._wp) then
            tsv%h_snow = acc_hsnow_num / acc_hsnow_den
        else
            tsv%h_snow = 0._wp
        end if
        if (acc_lai_den .gt. 0._wp) then
            tsv%lai = acc_lai_num / acc_lai_den
        else
            tsv%lai = 0._wp
        end if
        tsv%snow_vol = acc_snowvol

        return

    end subroutine build_ann_ts


    subroutine reset_year_accum()
        implicit none
        acc_prc    = 0._wp; acc_evp = 0._wp; acc_esur = 0._wp; acc_trans = 0._wp
        acc_run    = 0._wp; acc_runsur = 0._wp; acc_drain = 0._wp
        acc_gpp    = 0._wp; acc_npp = 0._wp
        acc_ch4    = 0._wp; acc_n2o = 0._wp
        acc_cflx   = 0._wp
        acc_t2m_num  = 0._wp; acc_t2m_den = 0._wp
        acc_hsnow_num= 0._wp; acc_hsnow_den = 0._wp
        acc_lai_num  = 0._wp; acc_lai_den = 0._wp
        acc_snowvol  = 0._wp
    end subroutine reset_year_accum


    !======================================================================
    ! hires allocation + coarse-index map + painting helpers.
    !======================================================================
    subroutine alloc_hires(h)

        implicit none

        type(hires_out), intent(inout) :: h

        allocate(h%h_snow(nx_hi, ny_hi))
        allocate(h%w_snow(nx_hi, ny_hi))
        allocate(h%t_skin(nx_hi, ny_hi))
        allocate(h%smb(nx_hi, ny_hi))
        allocate(h%melt(nx_hi, ny_hi))
        allocate(h%f_ice(nx_hi, ny_hi))
        allocate(h%veg_c(nx_hi, ny_hi))
        allocate(h%soil_c(nx_hi, ny_hi))

        h%h_snow = 0._wp; h%w_snow = 0._wp; h%t_skin = 0._wp
        h%smb    = 0._wp; h%melt   = 0._wp; h%f_ice  = 0._wp
        h%veg_c  = 0._wp; h%soil_c = 0._wp

        return

    end subroutine alloc_hires


    subroutine build_hires_map(geo)
        ! Build the (nx_hi, ny_hi) -> (ic, jc) nearest-neighbor mapping
        ! from the hires longitude/latitude arrays. Both grids are regular
        ! lat/lon so the index arithmetic is trivial; this materialises it
        ! once so paint_hires is a pure lookup.

        implicit none

        type(geo_class), intent(in) :: geo

        integer  :: ii, jj, ic, jc
        real(wp) :: lon_ii, lat_jj
        real(wp) :: dlon_c, dlat_c, lon0_c, lat0_c

        allocate(ij_hi_i(nx_hi, ny_hi))
        allocate(ij_hi_j(nx_hi, ny_hi))

        ! coarse (cmn) grid — 72x36 regular, lon = [-180, 180), lat = [-90, 90]
        dlon_c = 360._wp / real(ni, wp)
        dlat_c = 180._wp / real(nj, wp)
        lon0_c = -180._wp
        lat0_c =  -90._wp

        do jj = 1, ny_hi
            do ii = 1, nx_hi
                lon_ii = geo%hires%grid%lon(ii, jj)
                lat_jj = geo%hires%grid%lat(ii, jj)

                ic = int((lon_ii - lon0_c) / dlon_c) + 1
                jc = int((lat_jj - lat0_c) / dlat_c) + 1

                if (ic .lt. 1)  ic = 1
                if (ic .gt. ni) ic = ni
                if (jc .lt. 1)  jc = 1
                if (jc .gt. nj) jc = nj

                ij_hi_i(ii, jj) = ic
                ij_hi_j(ii, jj) = jc
            end do
        end do

        return

    end subroutine build_hires_map


    subroutine paint_hires(av, h)
        ! Nearest-neighbor paint from the coarse annual surf field onto the
        ! hires grid. Cell-uniform within each cmn cell at the identity
        ! baseline — a blocky picture — which is exactly what we want:
        ! it makes the framework visually obvious as it stands, and the
        ! same paint call becomes the hypsometric downscale later.

        implicit none

        type(surf_out),  intent(in)    :: av
        type(hires_out), intent(inout) :: h

        integer :: ii, jj, ic, jc

        !$omp parallel do collapse(2) private(ii,jj,ic,jc) schedule(static)
        do jj = 1, ny_hi
            do ii = 1, nx_hi
                ic = ij_hi_i(ii, jj)
                jc = ij_hi_j(ii, jj)
                h%h_snow(ii, jj) = av%h_snow(ic, jc)
                h%w_snow(ii, jj) = av%w_snow(ic, jc)
                h%t_skin(ii, jj) = av%t_skin(ic, jc)
                h%smb(ii, jj)    = av%smb(ic, jc)
                h%melt(ii, jj)   = av%melt(ic, jc)
                h%f_ice(ii, jj)  = av%f_ice(ic, jc)
                h%veg_c(ii, jj)  = av%veg_c(ic, jc)
                h%soil_c(ii, jj) = av%soil_c(ic, jc)
            end do
        end do
        !$omp end parallel do

        return

    end subroutine paint_hires

end module lndvc_out
