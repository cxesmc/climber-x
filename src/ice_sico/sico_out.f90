module sico_out

  use sico_types_m, only : wp
  use dim_name, only: dim_x, dim_y, dim_time, dim_month, dim_kc, dim_kt, dim_kr
  use timer, only : year, year_now, sec_year
  use timer, only : time_out_ice, time_out_ts, ny_out_ts, y_out_ts
  use control, only : out_dir
  use sico_grid_mod
  use sicopolis, only : sico_class, ts_out_class
  use sico_params, only : pi_180_inv, RHO_SW, RHO, RHO_W, G
  use enth_temp_omega_m, only : enth_fct_temp_omega
  use ncio

  implicit none

  real(wp), parameter :: no_value_neg_2 = -9.999e+03_wp

  private
  public :: sico_diag, sico_diag_init


contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  i c e _ d i a g _ i n i t
  ! Purpose  :  Initialize netcdf output for ice sheet
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sico_diag_init(ice)

    implicit none

    type(sico_class) :: ice
    character (len=256) :: fnm
 
    ! allocate
    allocate(ice%ts(ny_out_ts))

    fnm = trim(out_dir)//"/ice_"//trim(ice%grid%grid1%name)//"_ts.nc"
    call ts_nc(fnm)

    fnm = trim(out_dir)//"/ice_"//trim(ice%grid%grid1%name)//".nc"
    call ice_nc(ice%grid, fnm)

    ice%nout = 0

   return

  end subroutine sico_diag_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i c e _ d i a g
  !   Purpose    :  ice sheet diagnostics
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sico_diag(ice)

    implicit none

    type(sico_class) :: ice

    integer :: i, j, y
    integer :: ppos
    character (len=256) :: fnm
    character (len=256) :: dom
    real(wp) :: vs_help, Tbh_help
    real(wp) :: MB, LMT, GIMB, SIMB, LMH, LQH, OMH, PAT, PAH, MBMIS
    real(wp), parameter :: rhosw_rho_ratio = RHO_SW/RHO
    real(wp), parameter :: A_surf = 3.61132e+14_wp              ! global ocean area, in m^2
    real(wp), dimension(:,:), allocatable :: H, H_cold, H_temp


    ! current index
    y = y_out_ts

    ! Computation of the ice volume, the volume above flotation,
    ! the area covered by grounded ice, the area covered by floating ice,
    ! the total surface mass balance, the total basal mass balance
    ! and the total calving flux --------

    ! Thickness of the cold and temperate layers --------

    allocate(H(0:ice%grid%JMAX,0:ice%grid%IMAX))
    allocate(H_cold(0:ice%grid%JMAX,0:ice%grid%IMAX))
    allocate(H_temp(0:ice%grid%JMAX,0:ice%grid%IMAX))

    H = ice%state%H_c + ice%state%H_t
    H_cold = 0.0_wp
    H_temp = 0.0_wp

    if (ice%par%calcmod==1) then
      do i=1, ice%grid%IMAX-1
        do j=1, ice%grid%JMAX-1
          H_temp(j,i) = ice%state%H_t(j,i)
        end do
      end do
    else if (ice%par%calcmod==0 .or. ice%par%calcmod==2 .or. ice%par%calcmod==3 .or. ice%par%calcmod==-1) then
      do i=1, ice%grid%IMAX-1
        do j=1, ice%grid%JMAX-1
          H_temp(j,i) = ice%state%H_c(j,i)*ice%grid%eaz_c_quotient(ice%state%kc_cts(j,i))
        end do
      end do
    endif

    H_cold = H - H_temp


    ice%ts(y)%V_grounded = 0.0_wp
    ice%ts(y)%V_floating = 0.0_wp
    ice%ts(y)%V_gr_redu  = 0.0_wp
    ice%ts(y)%A_grounded = 0.0_wp
    ice%ts(y)%A_floating = 0.0_wp
    ice%ts(y)%V_temp     = 0.0_wp
    ice%ts(y)%A_temp     = 0.0_wp
    ice%ts(y)%Q_s        = 0.0_wp
    ice%ts(y)%Q_b        = 0.0_wp
    ice%ts(y)%Q_temp     = 0.0_wp
    ice%ts(y)%H_max      = 0.0_wp
    ice%ts(y)%H_t_max    = 0.0_wp
    ice%ts(y)%zs_max     = no_value_neg_2
    ice%ts(y)%vs_max     = 0.0_wp
    ice%ts(y)%Tbh_max    = no_value_neg_2

    do i=0, ice%grid%IMAX
      do j=0, ice%grid%JMAX

        if (ice%state%maske(j,i)==0) then   ! grounded ice

          if (ice%state%zs(j,i)       > ice%ts(y)%zs_max)  ice%ts(y)%zs_max  = ice%state%zs(j,i)
          if (H(j,i)                  > ice%ts(y)%H_max  ) ice%ts(y)%H_max   = H(j,i)
          if (H_temp(j,i)             > ice%ts(y)%H_t_max) ice%ts(y)%H_t_max = H_temp(j,i)

          ice%ts(y)%V_grounded = ice%ts(y)%V_grounded + H(j,i)*ice%grid%area(j,i)*1.e-15_wp ! mln km3
          ice%ts(y)%V_temp     = ice%ts(y)%V_temp     + H_temp(j,i)*ice%grid%area(j,i)*1.e-15_wp ! mln km3
          ice%ts(y)%A_grounded = ice%ts(y)%A_grounded + ice%grid%area(j,i)*1.e-12_wp ! mln km2
          ice%ts(y)%V_gr_redu  = ice%ts(y)%V_gr_redu &
            + rhosw_rho_ratio*max((ice%state%z_sl(j,i)-ice%state%zl(j,i)),0.0_wp)*ice%grid%area(j,i)*1.e-15_wp ! mln km3

          if (ice%state%n_cts(j,i) /= -1) ice%ts(y)%A_temp = ice%ts(y)%A_temp + ice%grid%area(j,i)*1.e-12_wp ! mln km2

          vs_help = sqrt( &
            0.25_wp*(ice%state%vx_c(ice%grid%KCMAX,j,i)+ice%state%vx_c(ice%grid%KCMAX,j,i-1))**2 &
            +0.25_wp*(ice%state%vy_c(ice%grid%KCMAX,j,i)+ice%state%vy_c(ice%grid%KCMAX,j-1,i))**2 )
          if (vs_help > ice%ts(y)%vs_max) ice%ts(y)%vs_max = vs_help

          if (ice%state%n_cts(j,i) >= 0) then   ! temperate base
            ice%ts(y)%Tbh_max = 0.0_wp
          else   ! cold base
            Tbh_help = min((ice%state%temp_c(0,j,i)-ice%state%temp_c_m(0,j,i)), 0.0_wp)
            if (Tbh_help > ice%ts(y)%Tbh_max) ice%ts(y)%Tbh_max = Tbh_help
          end if

        else if (ice%state%maske(j,i)==3) then   ! floating ice

          if (ice%state%zs(j,i)    > ice%ts(y)%zs_max) ice%ts(y)%zs_max = ice%state%zs(j,i)
          if (H(j,i)               > ice%ts(y)%H_max)  ice%ts(y)%H_max  = H(j,i)

          ice%ts(y)%V_floating = ice%ts(y)%V_floating + (ice%state%H_c(j,i)+ice%state%H_t(j,i))*ice%grid%area(j,i)*1.e-15_wp ! mln km3
          ice%ts(y)%A_floating = ice%ts(y)%A_floating + ice%grid%area(j,i)*1.e-12_wp ! mln km2

          vs_help = sqrt( &
            0.25_wp*(ice%state%vx_c(ice%grid%KCMAX,j,i)+ice%state%vx_c(ice%grid%KCMAX,j,i-1))**2 &
            +0.25_wp*(ice%state%vy_c(ice%grid%KCMAX,j,i)+ice%state%vy_c(ice%grid%KCMAX,j-1,i))**2 )
          if (vs_help > ice%ts(y)%vs_max) ice%ts(y)%vs_max = vs_help

          Tbh_help = min((ice%state%temp_c(0,j,i)-ice%state%temp_c_m(0,j,i)), 0.0_wp)
          if (Tbh_help > ice%ts(y)%Tbh_max) ice%ts(y)%Tbh_max = Tbh_help

        end if

        ice%ts(y)%Q_s   = ice%ts(y)%Q_s  + ice%state%as_perp_apl(j,i) * ice%grid%area(j,i)

        if (     (ice%state%maske(j,i)==0).or.(ice%state%maske_old(j,i)==0) &
          .or.(ice%state%maske(j,i)==3).or.(ice%state%maske_old(j,i)==3) ) then
          ! grounded or floating ice before or after the time step

          ice%ts(y)%Q_b   = ice%ts(y)%Q_b  + ice%state%Q_bm(j,i)        * ice%grid%area(j,i)

        end if

        if ( (ice%state%maske(j,i)==0).or.(ice%state%maske_old(j,i)==0) ) then
          ! grounded ice before or after the time step

          ice%ts(y)%Q_temp = ice%ts(y)%Q_temp + ice%state%Q_tld(j,i) * ice%grid%area(j,i)

        end if
      end do
    end do

    ice%ts(y)%A_tot   = ice%ts(y)%A_grounded + ice%ts(y)%A_floating

    ice%ts(y)%V_tot   = ice%ts(y)%V_grounded + ice%ts(y)%V_floating
    ice%ts(y)%V_af    = ice%ts(y)%V_grounded - ice%ts(y)%V_gr_redu

    ice%ts(y)%V_sle   = ice%ts(y)%V_af*1.e15_wp*(RHO/RHO_W)/A_surf     ! m^3 ice equiv./m^2 -> m water equiv.

    ice%ts(y)%Q_s       = ice%ts(y)%Q_s    *sec_year *RHO *1.e-12_wp ! m^3/s ice equiv. -> Gt/a
    ice%ts(y)%Q_b       = ice%ts(y)%Q_b    *sec_year *RHO *1.e-12_wp ! m^3/s ice equiv. -> Gt/a
    ice%ts(y)%Q_temp    = ice%ts(y)%Q_temp *sec_year *RHO *1.e-12_wp ! m^3/s ice equiv. -> Gt/a

    ice%ts(y)%vs_max    = ice%ts(y)%vs_max *sec_year  ! m/s -> m/a

    !  ------ Computation of further data

    ! Computation of volume fluxes

    ice%ts(y)%dV_dt      = 0.0_wp
    ice%ts(y)%precip_tot = 0.0_wp

    do i=1, ice%grid%IMAX-1
      do j=1, ice%grid%JMAX-1

        if (     (ice%state%maske(j,i)==0).or.(ice%state%maske_old(j,i)==0) &
          .or.(ice%state%maske(j,i)==3).or.(ice%state%maske_old(j,i)==3) ) then
          ! grounded or floating ice before or after the time step

          ! Volume change
          ice%ts(y)%dV_dt = ice%ts(y)%dV_dt + (ice%state%dzs_dtau(j,i)-ice%state%dzb_dtau(j,i))*ice%grid%area(j,i)
          ice%ts(y)%precip_tot    = ice%ts(y)%precip_tot + ice%state%accum_apl(j,i)*ice%grid%area(j,i)

        end if

      end do
    end do

    ice%ts(y)%precip_tot = ice%ts(y)%precip_tot *sec_year *RHO *1.e-12_wp ! m^3/s ice equiv. -> Gt/a 
    ice%ts(y)%dV_dt      = ice%ts(y)%dV_dt      *sec_year ! m^3/s ice equiv. -> m^3/a ice equiv.


    ! MB:   total mass balance as computed in subroutine apply_smb
    ! LMT:  total lost on land at the top of the ice sheet including shelf ice
    ! GIMB: total basal lost grounded ice
    ! SIMB: total basal lost shelf ice
    ! LMH:  total hidden lost through runoff on land
    ! OMH:  total hidden lost through flow into the ocean
    ! PAT:  total lost through ice discharge/calving parameterizations
    ! PAH:  total hidden lost through ice discharge/calving parameterization 
    ! LQH:  total hidden lost through melt at the base on land
    ! MBMIS: misaccounted

    MB    = 0.0_dp
    LMT   = 0.0_dp;
    GIMB  = 0.0_dp;   SIMB = 0.0_dp
    LMH   = 0.0_dp;   LQH  = 0.0_dp
    OMH   = 0.0_dp
    PAT   = 0.0_dp;   PAH  = 0.0_dp
    MBMIS = 0.0_dp;

    do i=1, ice%grid%IMAX-1
      do j=1, ice%grid%JMAX-1

        if ( ice%state%mask_ablation_type(j,i) /= 0 ) then
          ! glaciated land and ocean (including hidden melt points)

          ! Quantify what types of melt occurred
          select case ( ice%state%mask_ablation_type(j,i) )
          case( 3 )
            LMT  = LMT + ice%state%runoff_apl(j,i)  * ice%grid%area(j,i)
            PAT  = PAT + ice%state%calving_apl(j,i) * ice%grid%area(j,i)  ! could cause problems for Greenland
            SIMB = SIMB + ice%state%Q_b_apl(j,i)    * ice%grid%area(j,i)
          case( 1 )
            LMT  = LMT + ice%state%runoff_apl(j,i)  * ice%grid%area(j,i)
            PAT  = PAT + ice%state%calving_apl(j,i) * ice%grid%area(j,i)  ! ok
            GIMB = GIMB + ice%state%Q_b_apl(j,i)    * ice%grid%area(j,i)
          case( 9 )
            MBMIS = MBMIS + ice%state%mb_source_apl(j,i)*ice%grid%area(j,i)
          case( -1 )
            LMH = LMH + ice%state%runoff_apl(j,i)  * ice%grid%area(j,i)
            PAH = PAH + ice%state%calving_apl(j,i) * ice%grid%area(j,i)
            LQH = LQH + ice%state%Q_b_apl(j,i)     * ice%grid%area(j,i)
          case( -2 )
            OMH = OMH + ice%state%calving_apl(j,i)  * ice%grid%area(j,i)   ! only one contribution
          end select

        end if

        ! Actual ice mass balance (from top melt, bottom melt and calving)
        MB = MB + ice%state%mb_source_apl(j,i)*ice%grid%area(j,i)

      end do
    end do

    ! total mass balance 
    ice%ts(y)%mb_tot     = MB
    ! Runoff on land (excluding basal melt)
    ice%ts(y)%runoff_tot = LMT + LMH
    ! Ice discharge (excluding basal melt)
    ice%ts(y)%calv_tot   = OMH + PAT + PAH
    ! Ice discharge from ice flow, large scale
    ice%ts(y)%disc_lsc   = OMH
    ! Ice discharge from parameterization, small scale
    ice%ts(y)%disc_ssc   = PAT + PAH
    ! Basal mass balance
    ice%ts(y)%bmb_tot    = -GIMB-SIMB-LQH
    ! grounded ice
    ice%ts(y)%bmb_gr_tot = -GIMB-LQH
    ! shelf ice
    ice%ts(y)%bmb_fl_tot = -SIMB  ! hidden is counted as large scale calving

    ice%ts(y)%mb_tot     = ice%ts(y)%mb_tot     * sec_year *RHO *1.e-12_wp ! m^3/s ice equiv. -> Gt/a 
    ice%ts(y)%disc_lsc   = ice%ts(y)%disc_lsc   * sec_year *RHO *1.e-12_wp ! m^3/s ice equiv. -> Gt/a 
    ice%ts(y)%disc_ssc   = ice%ts(y)%disc_ssc   * sec_year *RHO *1.e-12_wp ! m^3/s ice equiv. -> Gt/a 
    ice%ts(y)%bmb_tot    = ice%ts(y)%bmb_tot    * sec_year *RHO *1.e-12_wp ! m^3/s ice equiv. -> Gt/a 
    ice%ts(y)%bmb_fl_tot = ice%ts(y)%bmb_fl_tot * sec_year *RHO *1.e-12_wp ! m^3/s ice equiv. -> Gt/a 
    ice%ts(y)%bmb_gr_tot = ice%ts(y)%bmb_gr_tot * sec_year *RHO *1.e-12_wp ! m^3/s ice equiv. -> Gt/a 
    ice%ts(y)%runoff_tot = ice%ts(y)%runoff_tot * sec_year *RHO *1.e-12_wp ! m^3/s ice equiv. -> Gt/a 
    ice%ts(y)%calv_tot   = ice%ts(y)%calv_tot   * sec_year *RHO *1.e-12_wp ! m^3/s ice equiv. -> Gt/a 

    if (ice%ts(y)%precip_tot /= 0.0_wp) then
      ice%ts(y)%mbp = ice%ts(y)%calv_tot/ice%ts(y)%precip_tot
    else
      ice%ts(y)%mbp = 0.0_wp
    end if

    ice%ts(y)%mb_resid = ice%ts(y)%Q_s + ice%ts(y)%bmb_tot - ice%ts(y)%calv_tot - ice%ts(y)%dV_dt*RHO*1.e-12_wp  ! Gt/a

    deallocate(H)
    deallocate(H_cold)
    deallocate(H_temp)


    if (time_out_ts) then
      fnm = trim(out_dir)//"/ice_"//trim(ice%grid%grid1%name)//"_ts.nc"
      call ts_nc_write(fnm,ice%ts(1:y),year-y+1)
    endif

    ppos = scan(trim(ice%grid%grid1%name),"-")-1
    dom = trim(ice%grid%grid1%name(1:ppos)) 
    ! print header
    if (mod(year,10).eq.1) then 
      print '(a7,a9,10a7)','ice_'//dom,'year','Vtot','Vgrd','Vflt','Atot','Agrd','Aflt','SMB','BMB','CLV','RUN'
    endif

    ! print values
    print '(a7,i9,6F7.1,4F7.0)', &
      'ice_'//dom,year_now, ice%ts(y)%V_tot, ice%ts(y)%V_grounded, ice%ts(y)%V_floating, &
      ice%ts(y)%A_tot, ice%ts(y)%A_grounded, ice%ts(y)%A_floating, ice%ts(y)%Q_s, ice%ts(y)%bmb_tot, ice%ts(y)%calv_tot, ice%ts(y)%runoff_tot


    ! spatially explicit output
    if (time_out_ice) then

      ice%nout = ice%nout+1
      call sico_diag_out(ice)

    endif


   return

  end subroutine sico_diag


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  i c e _ d i a g _ o u t
  ! Purpose  :  write ice sheet netcdf output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sico_diag_out(ice)

    implicit none

    type(sico_class) :: ice
    integer :: ncid
    character (len=256) :: fnm


    ! write to file
    fnm = trim(out_dir)//"/ice_"//trim(ice%grid%grid1%name)//".nc"
    call nc_open(fnm,ncid,writable=.TRUE.)
    call nc_write(fnm,dim_time,dble(year_now), dim1=dim_time, start=[ice%nout], count=[1],ncid=ncid)    
    call ice_nc_write(ice,fnm,ncid,ice%nout)
    call nc_close(ncid)


   return

  end subroutine sico_diag_out
  

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t s _ n c
  ! Purpose  :  initialize netcdf file for time series output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ts_nc(fnm)

    implicit none

    character (len=*) :: fnm
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_write_dim(fnm,dim_time, x=empty_time, axis="t", units="years BP", unlimited=.TRUE.)
    call nc_write_dim(fnm,dim_y, x=1, axis="y", units="1")
    call nc_write_dim(fnm,dim_x, x=1, axis="x", units="1")
    return

  end subroutine ts_nc
 

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t s _ n c _ w r i t e
  ! Purpose  :  write time series to netcdf
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ts_nc_write(fnm,ts,nout)

    implicit none

    type(ts_out_class) :: ts(:)

    character (len=*) :: fnm
    integer :: nout, ncid, i


    call nc_open(fnm,ncid,writable=.TRUE.)
    call nc_write(fnm,"time",  dble([(i,i=(year_now-y_out_ts+1),(year_now))]), dim1=dim_time,start=[nout],count=[y_out_ts],ncid=ncid)    
    call nc_write(fnm,"V_tot", ts%V_tot,           dims=[dim_time],start=[nout],count=[y_out_ts],long_name="total ice volume",units="mln km3",ncid=ncid)
    call nc_write(fnm,"V_sle", ts%V_sle,           dims=[dim_time],start=[nout],count=[y_out_ts],long_name="ice volume in SLE",units="m sle",ncid=ncid)
    call nc_write(fnm,"V_temp",ts%V_temp,          dims=[dim_time],start=[nout],count=[y_out_ts],long_name="volume of temperate ice",units="mln km3",ncid=ncid)
    call nc_write(fnm,"V_grd", ts%V_grounded,      dims=[dim_time],start=[nout],count=[y_out_ts],long_name="grounded ice volume",units="mln km3",ncid=ncid)
    call nc_write(fnm,"V_flt", ts%V_floating,      dims=[dim_time],start=[nout],count=[y_out_ts],long_name="floating ice volume",units="mln km3",ncid=ncid)
    call nc_write(fnm,"V_af",  ts%V_af,            dims=[dim_time],start=[nout],count=[y_out_ts],long_name="ice volume above floatation",units="mln km3",ncid=ncid)
    call nc_write(fnm,"dV_dt", ts%dV_dt,           dims=[dim_time],start=[nout],count=[y_out_ts],long_name="rate of ice volume change",units="m3/yr",ncid=ncid)
    call nc_write(fnm,"A_tot", ts%A_tot,           dims=[dim_time],start=[nout],count=[y_out_ts],long_name="total ice area",units="mln km2",ncid=ncid)
    call nc_write(fnm,"A_grd", ts%A_grounded,      dims=[dim_time],start=[nout],count=[y_out_ts],long_name="grounded ice area",units="mln km2",ncid=ncid)
    call nc_write(fnm,"A_flt", ts%A_floating,      dims=[dim_time],start=[nout],count=[y_out_ts],long_name="floating ice area",units="mln km2",ncid=ncid)
    call nc_write(fnm,"A_temp",ts%A_temp,          dims=[dim_time],start=[nout],count=[y_out_ts],long_name="area covered by temperate ice",units="mln km2",ncid=ncid)
    call nc_write(fnm,"f_temp",ts%A_temp/max(1e-20,ts%A_grounded), dims=[dim_time],start=[nout],count=[y_out_ts],long_name="fraction of basal area covered by temperate ice",units="1",ncid=ncid)
    call nc_write(fnm,"mb_tot", ts%mb_tot,         dims=[dim_time],start=[nout],count=[y_out_ts],long_name="total mass balance",units="Gt/a",ncid=ncid)
    call nc_write(fnm,"smb_tot", ts%Q_s,           dims=[dim_time],start=[nout],count=[y_out_ts],long_name="total surface mass balance",units="Gt/a",ncid=ncid)
    call nc_write(fnm,"bmb_tot", ts%bmb_tot,       dims=[dim_time],start=[nout],count=[y_out_ts],long_name="total basal mass balance",units="Gt/a",ncid=ncid)
    call nc_write(fnm,"bmb_gr_tot", ts%bmb_gr_tot, dims=[dim_time],start=[nout],count=[y_out_ts],long_name="total basal mass balance of grounded ice",units="Gt/a",ncid=ncid)
    call nc_write(fnm,"bmb_fl_tot", ts%bmb_fl_tot, dims=[dim_time],start=[nout],count=[y_out_ts],long_name="total basal mass balance of floating ice",units="Gt/a",ncid=ncid)
    call nc_write(fnm,"runoff_tot", ts%runoff_tot, dims=[dim_time],start=[nout],count=[y_out_ts],long_name="total runoff",units="Gt/a",ncid=ncid)
    call nc_write(fnm,"calv_tot", ts%calv_tot,     dims=[dim_time],start=[nout],count=[y_out_ts],long_name="total calving",units="Gt/a",ncid=ncid)
    call nc_write(fnm,"disc_lsc", ts%disc_lsc,     dims=[dim_time],start=[nout],count=[y_out_ts],long_name="Ice discharge from ice flow, large scale",units="Gt/a",ncid=ncid)
    call nc_write(fnm,"disc_ssc", ts%disc_ssc,     dims=[dim_time],start=[nout],count=[y_out_ts],long_name="Ice discharge from parameterization, small scale",units="Gt/a",ncid=ncid)
    call nc_write(fnm,"Q_b", ts%Q_b,               dims=[dim_time],start=[nout],count=[y_out_ts],long_name="basal melting rate",units="Gt/a",ncid=ncid)
    call nc_write(fnm,"Q_temp", ts%Q_temp,         dims=[dim_time],start=[nout],count=[y_out_ts],long_name="drainage rate from the temperate ice layer",units="Gt/a",ncid=ncid)
    call nc_write(fnm,"mb_resid", ts%mb_resid,     dims=[dim_time],start=[nout],count=[y_out_ts],long_name="residual of the total mass balance",units="Gt/a",ncid=ncid)
    call nc_write(fnm,"mbp", ts%mbp,               dims=[dim_time],start=[nout],count=[y_out_ts],long_name="mass balance partition",units="1",ncid=ncid)
    call nc_write(fnm,"H_max", ts%H_max,           dims=[dim_time],start=[nout],count=[y_out_ts],long_name="maximum ice thickness",units="m",ncid=ncid)
    call nc_write(fnm,"H_t_max", ts%H_t_max,       dims=[dim_time],start=[nout],count=[y_out_ts],long_name="maximum thickness of temperate ice",units="m",ncid=ncid)
    call nc_write(fnm,"zs_max", ts%zs_max,         dims=[dim_time],start=[nout],count=[y_out_ts],long_name="maximum surface elevation",units="m",ncid=ncid)
    call nc_write(fnm,"vs_max", ts%vs_max,         dims=[dim_time],start=[nout],count=[y_out_ts],long_name="maximum surface speed",units="m/yr",ncid=ncid)
    call nc_write(fnm,"Tbh_max", ts%Tbh_max,       dims=[dim_time],start=[nout],count=[y_out_ts],long_name="maximum basal temperature relative to pmp",units="decG",ncid=ncid)
    call nc_close(ncid)

!#if (DISC>0)
!
!!    ---- disc_lsc
!
!   call check( nf90_inq_dimid(ncid, 't', nc1d), thisroutine )
!   call check( nf90_def_var(ncid, 'disc_lsc', NF90_FLOAT, nc1d, ncv), &
!               thisroutine )
!   buffer = 'm3 ice equiv. a-1'
!   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
!               thisroutine )
!   buffer = 'large_scale_ice_lost_into_the_ocean'
!   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
!               thisroutine )
!   buffer = 'Large scale ice lost into the ocean'
!   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
!               thisroutine )
!
!!    ---- disc_ssc
!
!   call check( nf90_inq_dimid(ncid, 't', nc1d), thisroutine )
!   call check( nf90_def_var(ncid, 'disc_ssc', NF90_FLOAT, nc1d, ncv), &
!               thisroutine )
!   buffer = 'm3 ice equiv. a-1'
!   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
!               thisroutine )
!   buffer = 'small_scale_ice_lost_into_the_ocean'
!   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
!               thisroutine )
!   buffer = 'Small scale ice lost into the ocean'
!   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
!               thisroutine )
!
!!    ---- dT_glann
!
!   call check( nf90_inq_dimid(ncid, 't', nc1d), thisroutine )
!   call check( nf90_def_var(ncid, 'dT_glann', NF90_FLOAT, nc1d, ncv), &
!               thisroutine )
!   buffer = 'degC'
!   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
!               thisroutine )
!   buffer = 'global_annual_temperature_anomaly'
!   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
!               thisroutine )
!   buffer = 'Global annual temperature anomaly'
!   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
!               thisroutine )
!
!!    ---- dT_sub
!
!   call check( nf90_inq_dimid(ncid, 't', nc1d), thisroutine )
!   call check( nf90_def_var(ncid, 'dT_sub', NF90_FLOAT, nc1d, ncv), &
!               thisroutine )
!   buffer = 'degC'
!   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
!               thisroutine )
!   buffer = 'subsurface_ocean_temperature_anomaly'
!   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
!               thisroutine )
!   buffer = 'Subsurface ocean temperature anomaly'
!   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
!               thisroutine )
!
!#endif


   return

  end subroutine ts_nc_write


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  i c e _ n c
  ! Purpose  :  Initialize ice sheet netcdf output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ice_nc(grd, fnm)

    implicit none

    type(sico_grid_class) :: grd
    character (len=*) :: fnm

    integer :: ncid
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid,writable=.TRUE.)
    call nc_write_dim(fnm, dim_time, x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE.,ncid=ncid)
    call nc_write_dim(fnm,dim_kt,x=1,dx=1,nx=grd%KTMAX+1,axis="e", units="1",ncid=ncid)
    call nc_write_dim(fnm,dim_kc,x=1,dx=1,nx=grd%KCMAX+1,axis="z", units="1",ncid=ncid)
    call nc_write_dim(fnm,dim_y,x=grd%y0,dx=grd%dx,nx=grd%JMAX+1,axis="y",units="km",ncid=ncid)
    call nc_write_dim(fnm,dim_x,x=grd%x0,dx=grd%dx,nx=grd%IMAX+1,axis="x",units="km",ncid=ncid)
    call nc_write_dim(fnm,dim_kr,x=1,dx=1,nx=grd%KRMAX+1,units="1",ncid=ncid)

    call nc_write(fnm,"xi", grd%xi, &
          dims=[dim_x],start=[1],count=[grd%IMAX+1], &
          long_name="x-coordinate",grid_mapping="polar_stereographic",units="m",ncid=ncid)    
    call nc_write(fnm,"eta", grd%eta, &
          dims=[dim_y],start=[1],count=[grd%JMAX+1], &
          long_name="y-coordinate",grid_mapping="polar_stereographic",units="m",ncid=ncid)    
    call nc_write(fnm,"phi", transpose(grd%phi)*pi_180_inv, &
          dims=[dim_x,dim_y],start=[1,1],count=[grd%IMAX+1,grd%JMAX+1], &
          long_name="geographic latitude",grid_mapping="polar_stereographic",units="degN",ncid=ncid)    
    call nc_write(fnm,"lambda", transpose(grd%lambda)*pi_180_inv, &
          dims=[dim_x,dim_y],start=[1,1],count=[grd%IMAX+1,grd%JMAX+1], &
          long_name="geographic longitude",grid_mapping="polar_stereographic",units="degE",ncid=ncid)    
    call nc_write(fnm,"sigma_level_c", grd%eaz_c_quotient, &
          dims=[dim_kc],start=[1],count=[grd%KCMAX+1], &
          long_name="sigma level in cold ice",grid_mapping="polar_stereographic",units="1",ncid=ncid)    
    call nc_write(fnm,"sigma_level_t", grd%zeta_t, &
          dims=[dim_kt],start=[1],count=[grd%KTMAX+1], &
          long_name="sigma level in temperate ice",grid_mapping="polar_stereographic",units="1",ncid=ncid)    
    call nc_write(fnm,"sigma_level_r", grd%zeta_r, &
          dims=[dim_kr],start=[1],count=[grd%KRMAX+1], &
          long_name="sigma level in litosphere",grid_mapping="polar_stereographic",units="1",ncid=ncid)    
    call nc_close(ncid)

   return

  end subroutine ice_nc


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  i c e _ n c _ w r i t e
  ! Purpose  :  Output of ice sheet netcdf
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ice_nc_write(ice,fnm,ncid,nout)

    implicit none

    type(sico_class) :: ice

    character (len=*) :: fnm
    integer :: nout, ncid, i, j, kc, kt
    real(wp), dimension(:,:), allocatable :: H, H_cold, H_temp
    real(wp), dimension(:,:), allocatable :: vx_m_g, vy_m_g 
    real(wp), dimension(:,:), allocatable :: tau_b_driving, tau_b_drag


  call nc_write(fnm,"temp_s", real(transpose(ice%state%temp_s),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="surface temperature",grid_mapping="polar_stereographic",units="degC",ncid=ncid)    

  call nc_write(fnm,"temp_g", real(transpose(ice%state%temp_g),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="ground temperature",grid_mapping="polar_stereographic",units="degC",ncid=ncid)    

  call nc_write(fnm,"mb_source", real(transpose(ice%state%mb_source)*sec_year,sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="total mass balance",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"mb_source_apl", real(transpose(ice%state%mb_source_apl)*sec_year,sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="applied total mass balance",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"as_perp", real(transpose(ice%state%as_perp)*sec_year,sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="surface mass balance",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"as_perp_apl", real(transpose(ice%state%as_perp_apl)*sec_year,sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="applied surface mass balance",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"accum", real(transpose(ice%state%accum)*sec_year,sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="accumulation flux",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"accum_apl", real(transpose(ice%state%accum_apl)*sec_year,sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="applied accumulation flux",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"calving", real(transpose(ice%state%calving)*sec_year,sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="calving flux",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"calving_apl", real(transpose(ice%state%calving_apl)*sec_year,sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="applied calving flux",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"dis_perp", real(transpose(ice%state%dis_perp)*sec_year,sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="applied calving flux from discharge parameterisation",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"cst_dist", real(transpose(ice%state%cst_dist)*1.e-3,sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="coastal distance",grid_mapping="polar_stereographic",units="km",ncid=ncid)    

  call nc_write(fnm,"cos_grad_tc", real(transpose(ice%state%cos_grad_tc),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Cosine of angle between surface gradient and cst dist gradient",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

  call nc_write(fnm,"mask_mar", transpose(ice%state%mask_mar), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="margina ring mask",grid_mapping="polar_stereographic",units="1",ncid=ncid)    

  call nc_write(fnm,"runoff", real(transpose(ice%state%runoff)*sec_year,sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="runoff flux",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"runoff_apl", real(transpose(ice%state%runoff_apl)*sec_year,sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="applied runoff flux",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"Q_b", real(transpose(ice%state%Q_b_tot*sec_year),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Basal melting rate + Water drainage from the temperate layer",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"Q_b_apl", real(transpose(ice%state%Q_b_apl*sec_year),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="applied Basal melting rate + Water drainage from the temperate layer",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"Q_bm", real(transpose(ice%state%Q_bm*sec_year),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Basal melting rate",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"Q_tld", real(transpose(ice%state%Q_tld*sec_year),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Water drainage from the temperate layer",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"q_geo", real(transpose(ice%state%q_geo),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="geothermal heat flux",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid)    

  call nc_write(fnm,"maske", int(transpose(ice%state%maske),4), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="ice-land-sea mask",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

  call nc_write(fnm,"maske_old", int(transpose(ice%state%maske_old),4), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="ice-land-sea mask (old)",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

  call nc_write(fnm,"f_sed", transpose(ice%state%f_sed), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="sediment fraction",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

  call nc_write(fnm,"H_sed", transpose(ice%state%h_sed), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="sediment thickness",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

  call nc_write(fnm,"c_slide", transpose(ice%state%c_slide), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="basal sliding coefficient",grid_mapping="polar_stereographic",units="m/(a*Pa^(p-q))",ncid=ncid)    

  call nc_write(fnm,"c_fric", transpose(ice%state%c_fric), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="basal friction coefficient",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

  call nc_write(fnm,"vb_t", transpose(ice%state%vb_t)*sec_year, &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="threshold basal velocity for regularized Coulomb law",grid_mapping="polar_stereographic",units="m/a",ncid=ncid)    

  call nc_write(fnm,"delta", transpose(ice%state%delta), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Fraction of overburden pressure at saturation",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

  if (ice%par%thk_evol==4) then
  call nc_write(fnm,"mask_maxextent", int(transpose(ice%state%mask_maxextent),4), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="mask of maximum allowed ice sheet extent",grid_mapping="polar_stereographic",units="/",ncid=ncid)    
  endif

  call nc_write(fnm,"n_cts", int(transpose(ice%state%n_cts),4), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="mask for polythermal conditions",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

  call nc_write(fnm,"kc_cts", int(transpose(ice%state%kc_cts),4), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Grid index of the CTS position",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

  call nc_write(fnm,"mask_ablation_type", int(transpose(ice%state%mask_ablation_type),4), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Mask indicating ablation type",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

  call nc_write(fnm,"id_mask", int(transpose(ice%state%id_mask),4), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Mask indicating ice IDs",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

  call nc_write(fnm,"z_sl", real(transpose(ice%state%z_sl),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="sea level",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

  call nc_write(fnm,"zs", real(transpose(ice%state%zs),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Topography of the free surface",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

  call nc_write(fnm,"zm", real(transpose(ice%state%zm),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Topography of z=zm interface",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

  call nc_write(fnm,"zb", real(transpose(ice%state%zb),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Topography of the ice base",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

  call nc_write(fnm,"zl", real(transpose(ice%state%zl),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Topography of the lithosphere surface",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

  call nc_write(fnm,"zl_std", real(transpose(ice%state%zl_std),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="sub-grid standard deviation of topography of the lithosphere surface",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

  call nc_write(fnm,"zl_fil", real(transpose(ice%state%zl_fil),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Topography of the filtered lithosphere surface",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

  call nc_write(fnm,"zl0", real(transpose(ice%state%zl0),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Topography of isostatically relaxed lithosphere surface",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

  allocate(H(0:ice%grid%JMAX,0:ice%grid%IMAX))
  allocate(H_cold(0:ice%grid%JMAX,0:ice%grid%IMAX))
  allocate(H_temp(0:ice%grid%JMAX,0:ice%grid%IMAX))

  !-------- Ice thickness and time derivative --------
  H = ice%state%H_c + ice%state%H_t
  !-------- Thickness of the cold and temperate layers --------
  H_cold = 0.0_wp
  H_temp = 0.0_wp
  if (ice%par%calcmod==1) then
    do i=1, ice%grid%IMAX-1
      do j=1, ice%grid%JMAX-1
        H_temp(j,i) = ice%state%H_t(j,i)
      end do
    end do
  else if (ice%par%calcmod==0 .or. ice%par%calcmod==2 .or. ice%par%calcmod==3 .or. ice%par%calcmod==-1) then
    do i=1, ice%grid%IMAX-1
      do j=1, ice%grid%JMAX-1
        H_temp(j,i) = ice%state%H_c(j,i)*ice%grid%eaz_c_quotient(ice%state%kc_cts(j,i))
      end do
    end do
  endif
  H_cold = H - H_temp

  call nc_write(fnm,"H", real(transpose(H),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Thickness of ice",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

  call nc_write(fnm,"H_cold", real(transpose(H_cold),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Thickness of the cold ice layer",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

  call nc_write(fnm,"H_temp", real(transpose(H_temp),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="thickness of the temperate ice layer",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

  call nc_write(fnm,"H_flow", real(transpose(ice%state%H_neu_flow),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Thickness of ice after dynamics (without source term)",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

  deallocate(H_cold,H_temp)

  call nc_write(fnm,"H_eff", real(transpose(ice%state%H_eff),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Effective thickness at shelf ice front", grid_mapping="polar_stereographic", &
        units="m",ncid=ncid)

  call nc_write(fnm,"H_calv", real(transpose(ice%state%H_calv),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Calvin Threshold", grid_mapping="polar_stereographic", &
        units="m",ncid=ncid)

  call nc_write(fnm,"am_perp", real(transpose(ice%state%am_perp*sec_year),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Volume flux across the z=zm interface",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"qx", real(transpose(ice%state%qx*sec_year),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Horizontal volume flux qx",grid_mapping="polar_stereographic",units="m2/yr",ncid=ncid)    

  call nc_write(fnm,"qy", real(transpose(ice%state%qy*sec_year),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Horizontal volume flux qy",grid_mapping="polar_stereographic",units="m2/yr",ncid=ncid)    

  call nc_write(fnm,"dzs_dt", real(transpose(ice%state%dzs_dtau*sec_year),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Rate of change of the topography of the free surface",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"dzm_dt", real(transpose(ice%state%dzm_dtau*sec_year),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Rate of change of the topography of the z=zm interface",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"dzb_dt", real(transpose(ice%state%dzb_dtau*sec_year),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Rate of change of the topography of the ice base",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"dzl_dt", real(transpose(ice%state%dzl_dtau*sec_year),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Rate of change of the topography of the lithosphere surface",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"dH_c_dt", real(transpose(ice%state%dH_c_dtau*sec_year),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Rate of change of the thickness of the upper (kc) ice layer",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

  call nc_write(fnm,"dH_t_dt", real(transpose(ice%state%dH_t_dtau*sec_year),sp), &
        dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
        long_name="Rate of change of the thickness of the lower (kt) ice layer",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"dH_dt", real(transpose((ice%state%dH_c_dtau+ice%state%dH_t_dtau)*sec_year),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Rate of change of the ice thickness",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"vx_m_sia", real(transpose(ice%state%vx_m_sia*sec_year),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="x depth-averaged horizontal velocity from SIA ",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"vy_m_sia", real(transpose(ice%state%vy_m_sia*sec_year),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="y depth-averaged horizontal velocity from SIA ",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"vx_m_ssa", real(transpose(ice%state%vx_m_ssa*sec_year),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="x depth-averaged horizontal velocity from SSA ",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"vy_m_ssa", real(transpose(ice%state%vy_m_ssa*sec_year),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="y depth-averaged horizontal velocity from SSA ",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"vx_b_g", real(transpose(ice%state%vx_b_g*sec_year),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Horizontal velocity vx at the ice base",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"vy_b_g", real(transpose(ice%state%vy_b_g*sec_year),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Horizontal velocity vy at the ice base",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"vz_b", real(transpose(ice%state%vz_b*sec_year),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Vertical velocity vz at the ice base",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"vh_b", real(transpose(sqrt(ice%state%vx_b_g**2+ice%state%vy_b_g**2))*sec_year,sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Horizontal velocity vh at the ice base",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"vx_s_g", real(transpose(ice%state%vx_s_g*sec_year),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Horizontal velocity vx at the ice surface",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"vy_s_g", real(transpose(ice%state%vy_s_g*sec_year),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Horizontal velocity vy at the ice surface",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"vz_s", real(transpose(ice%state%vz_s*sec_year),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Vertical velocity vz at the ice surface",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"vh_s", real(transpose(sqrt(ice%state%vx_s_g**2+ice%state%vy_s_g**2))*sec_year,sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Horizontal velocity vh at the ice surface",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    


    !-------- Vertical mean of horizontal velocities --------
    allocate(vx_m_g(0:ice%grid%JMAX,0:ice%grid%IMAX))
    allocate(vy_m_g(0:ice%grid%JMAX,0:ice%grid%IMAX))

    vx_m_g = 0.0_wp
    vy_m_g = 0.0_wp

    do i=1, ice%grid%IMAX-1
      do j=1, ice%grid%JMAX-1
        if ( (ice%state%maske(j,i)==0).or.(ice%state%maske(j,i)==3) ) then
          vx_m_g(j,i) = 0.5_wp*(ice%state%vx_m(j,i)+ice%state%vx_m(j,i-1))
          vy_m_g(j,i) = 0.5_wp*(ice%state%vy_m(j,i)+ice%state%vy_m(j-1,i))
        end if
      end do
    end do

    call nc_write(fnm,"vx_m_g", real(transpose(vx_m_g*sec_year),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Vertical mean of horizontal velocity vx",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"vy_m_g", real(transpose(vy_m_g*sec_year),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Vertical mean of horizontal velocity vy",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"vh_m", real(transpose(sqrt(vx_m_g**2+vy_m_g**2))*sec_year,sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Vertical mean of horizontal velocity vh",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    ! SIA
    vx_m_g = 0.0_wp
    vy_m_g = 0.0_wp

    do i=1, ice%grid%IMAX-1
      do j=1, ice%grid%JMAX-1
        if ( (ice%state%maske(j,i)==0).or.(ice%state%maske(j,i)==3) ) then
          vx_m_g(j,i) = 0.5_wp*(ice%state%vx_m_sia(j,i)+ice%state%vx_m_sia(j,i-1))
          vy_m_g(j,i) = 0.5_wp*(ice%state%vy_m_sia(j,i)+ice%state%vy_m_sia(j-1,i))
        end if
      end do
    end do

    call nc_write(fnm,"vh_m_sia", real(transpose(sqrt(vx_m_g**2+vy_m_g**2))*sec_year,sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Vertical mean of SIA horizontal velocity",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    ! SSA
    vx_m_g = 0.0_wp
    vy_m_g = 0.0_wp

    do i=1, ice%grid%IMAX-1
      do j=1, ice%grid%JMAX-1
        if ( (ice%state%maske(j,i)==0).or.(ice%state%maske(j,i)==3) ) then
          vx_m_g(j,i) = 0.5_wp*(ice%state%vx_m_ssa(j,i)+ice%state%vx_m_ssa(j,i-1))
          vy_m_g(j,i) = 0.5_wp*(ice%state%vy_m_ssa(j,i)+ice%state%vy_m_ssa(j-1,i))
        end if
      end do
    end do

    call nc_write(fnm,"vh_m_ssa", real(transpose(sqrt(vx_m_g**2+vy_m_g**2))*sec_year,sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Vertical mean of SSA horizontal velocity",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    deallocate(vx_m_g,vy_m_g)

    call nc_write(fnm,"temp_b", real(transpose(ice%state%temp_b),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Temperature at the ice base",grid_mapping="polar_stereographic",units="degC",ncid=ncid)    

    call nc_write(fnm,"temph_b", real(transpose(ice%state%temph_b),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Temperature at the ice base relative to the pressure melting point",grid_mapping="polar_stereographic",units="degC",ncid=ncid)    

    !-------- Driving stress and basal drag --------
    allocate(tau_b_driving(0:ice%grid%JMAX,0:ice%grid%IMAX))
    allocate(tau_b_drag(0:ice%grid%JMAX,0:ice%grid%IMAX))

    tau_b_driving = 0.0_wp
    tau_b_drag    = 0.0_wp

    do i=0, ice%grid%IMAX
      do j=0, ice%grid%JMAX
        if (ice%state%maske(j,i)==0) then   ! grounded ice
          tau_b_driving(j,i) = RHO*G*H(j,i) * sqrt( ice%state%dzs_dxi_g(j,i)**2 + ice%state%dzs_deta_g(j,i)**2 )
          if (.not.ice%state%flag_shelfy_stream(j,i)) then
            tau_b_drag(j,i) = tau_b_driving(j,i)
          else
            tau_b_drag(j,i) = ice%state%beta_drag(j,i) * sqrt(ice%state%vx_b_g(j,i)**2+ice%state%vy_b_g(j,i)**2)
          end if
        else if (ice%state%maske(j,i)==3) then   ! floating ice
          tau_b_driving(j,i) = RHO*G*H(j,i) * sqrt( ice%state%dzs_dxi_g(j,i)**2 + ice%state%dzs_deta_g(j,i)**2 )
          tau_b_drag(j,i)    = 0.0_wp
        end if
      end do
    end do

    call nc_write(fnm,"tau_b_driving", real(transpose(tau_b_driving),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Driving stress",grid_mapping="polar_stereographic",units="Pa",ncid=ncid)    

    call nc_write(fnm,"tau_b_drag", real(transpose(tau_b_drag),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Basal drag",grid_mapping="polar_stereographic",units="Pa",ncid=ncid)    

    deallocate(tau_b_driving,tau_b_drag)
    deallocate(H)

    call nc_write(fnm,"p_b", real(transpose(ice%state%p_b),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Basal pressure",grid_mapping="polar_stereographic",units="Pa",ncid=ncid)    

    call nc_write(fnm,"p_b_w", real(transpose(ice%state%p_b_w),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Basal water pressure",grid_mapping="polar_stereographic",units="Pa",ncid=ncid)    

    call nc_write(fnm,"p_b_red_lim", real(transpose(ice%state%p_b_red_lim),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Reduced basal pressure",grid_mapping="polar_stereographic",units="Pa",ncid=ncid)    

    call nc_write(fnm,"H_w", real(transpose(ice%state%H_w),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Thickness of the water column under the ice base",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

    call nc_write(fnm,"q_gl_g", real(transpose(ice%state%q_gl_g)*sec_year,sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Horizontal volume flux across the grounding line",grid_mapping="polar_stereographic",units="m2/yr",ncid=ncid)    

    call nc_write(fnm,"ratio_sl_x", real(transpose(ice%state%ratio_sl_x),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Ratio of basal to surface velocity (slip ratio) in x-direction at (i+1/2,j)",grid_mapping="polar_stereographic",units="/",ncid=ncid)

    call nc_write(fnm,"ratio_sl_y", real(transpose(ice%state%ratio_sl_y),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Ratio of basal to surface velocity (slip ratio) in y-direction at (i+1/2,j)?",grid_mapping="polar_stereographic",units="/",ncid=ncid)

    call nc_write(fnm,"weigh_ssta_sia_x", real(transpose(ice%state%weigh_ssta_sia_x),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="weight of SStA vs SIA in x-direction at (i+1/2,j)",grid_mapping="polar_stereographic",units="/",ncid=ncid)

    call nc_write(fnm,"weigh_ssta_sia_y", real(transpose(ice%state%weigh_ssta_sia_y),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="weight of SStA vs SIA in y-direction at (i,j+1/2)",grid_mapping="polar_stereographic",units="/",ncid=ncid)

    call nc_write(fnm,"flag_shelfy_stream_x", transpose(ice%state%flag_shelfy_stream_x), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Shelfy stream flag in x-direction, at (i+1/2,j)",grid_mapping="polar_stereographic",units="/",ncid=ncid)

    call nc_write(fnm,"flag_shelfy_stream_y", transpose(ice%state%flag_shelfy_stream_y), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Shelfy stream flag in y-direction, at (i,j+1/2)",grid_mapping="polar_stereographic",units="/",ncid=ncid)

    call nc_write(fnm,"flag_shelfy_stream", transpose(ice%state%flag_shelfy_stream), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Shelfy stream flag",grid_mapping="polar_stereographic",units="/",ncid=ncid)

    call nc_write(fnm,"vis_int_g", real(transpose(ice%state%vis_int_g),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="Depth-integrated viscosity",grid_mapping="polar_stereographic",units="Pa s m",ncid=ncid)

    call nc_write(fnm,"beta_drag", real(transpose(ice%state%beta_drag),sp), &
          dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,1], &
          long_name="basal drag parameter for shelfy stream",grid_mapping="polar_stereographic",units="/",ncid=ncid)

if (ice%par%l_output_3d) then

    call nc_write(fnm,"vx_c", real(reshape(ice%state%vx_c,[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KCMAX+1],order=[3,2,1])*sec_year,sp), &
          dims=[dim_x,dim_y,dim_kc,dim_time],start=[1,1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KCMAX+1,1], &
          long_name="Horizontal velocity vx in the upper (kc) ice layer",grid_mapping="polar_stereographic",units="m/a",ncid=ncid)

    call nc_write(fnm,"vy_c", real(reshape(ice%state%vy_c,[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KCMAX+1],order=[3,2,1])*sec_year,sp), &
          dims=[dim_x,dim_y,dim_kc,dim_time],start=[1,1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KCMAX+1,1], &
          long_name="Horizontal velocity vy in the upper (kc) ice layer",grid_mapping="polar_stereographic",units="m/a",ncid=ncid)

    call nc_write(fnm,"vz_c", real(reshape(ice%state%vz_c,[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KCMAX+1],order=[3,2,1])*sec_year,sp), &
          dims=[dim_x,dim_y,dim_kc,dim_time],start=[1,1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KCMAX+1,1], &
          long_name="Vertical velocity vz in the upper (kc) ice layer",grid_mapping="polar_stereographic",units="m/a",ncid=ncid)

    call nc_write(fnm,"vx_t", real(reshape(ice%state%vx_t,[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KTMAX+1],order=[3,2,1])*sec_year,sp), &
          dims=[dim_x,dim_y,dim_kt,dim_time],start=[1,1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KTMAX+1,1], &
          long_name="Horizontal velocity vx in the lower (kt) ice layer",grid_mapping="polar_stereographic",units="m/a",ncid=ncid)

    call nc_write(fnm,"vy_t", real(reshape(ice%state%vy_t,[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KTMAX+1],order=[3,2,1])*sec_year,sp), &
          dims=[dim_x,dim_y,dim_kt,dim_time],start=[1,1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KTMAX+1,1], &
          long_name="Horizontal velocity vy in the lower (kt) ice layer",grid_mapping="polar_stereographic",units="m/a",ncid=ncid)

    call nc_write(fnm,"vz_t", real(reshape(ice%state%vz_t,[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KTMAX+1],order=[3,2,1])*sec_year,sp), &
          dims=[dim_x,dim_y,dim_kt,dim_time],start=[1,1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KTMAX+1,1], &
          long_name="Vertical velocity vz in the lower (kt) ice layer",grid_mapping="polar_stereographic",units="m/a",ncid=ncid)

    call nc_write(fnm,"temp_c", real(reshape(ice%state%temp_c,[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KCMAX+1],order=[3,2,1]),sp), &
          dims=[dim_x,dim_y,dim_kc,dim_time],start=[1,1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KCMAX+1,1], &
          long_name="Temperature in the upper (kc) ice layer",grid_mapping="polar_stereographic",units="degC",ncid=ncid)

    call nc_write(fnm,"omega_t", real(reshape(ice%state%omega_t,[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KTMAX+1],order=[3,2,1]),sp), &
          dims=[dim_x,dim_y,dim_kt,dim_time],start=[1,1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KTMAX+1,1], &
          long_name="Water content in the lower (kt) ice layer",grid_mapping="polar_stereographic",units="1",ncid=ncid)

    call nc_write(fnm,"omega_c", real(reshape(ice%state%omega_c,[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KCMAX+1],order=[3,2,1]),sp), &
          dims=[dim_x,dim_y,dim_kc,dim_time],start=[1,1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KCMAX+1,1], &
          long_name="Water content in the upper (kc) ice layer",grid_mapping="polar_stereographic",units="1",ncid=ncid)

    call nc_write(fnm,"temp_r", real(reshape(ice%state%temp_r,[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KRMAX+1],order=[3,2,1]),sp), &
          dims=[dim_x,dim_y,dim_kr,dim_time],start=[1,1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KRMAX+1,1], &
          long_name="Temperature in the lithosphere layer",grid_mapping="polar_stereographic",units="degC",ncid=ncid)

    !-------- Enthalpies for the non-enthalpy methods (POLY, COLD, ISOT) --------

    if (ice%par%calcmod==1) then

      do i=0, ice%grid%IMAX
        do j=0, ice%grid%JMAX
          do kc=0, ice%grid%KCMAX
            ice%state%enth_c(kc,j,i) = enth_fct_temp_omega(ice%state%temp_c(kc,j,i), 0.0_wp)
          end do
          if ( (ice%state%maske(j,i)==0).and.(ice%state%n_cts(j,i)==1) ) then
            do kt=0, ice%grid%KTMAX
              ice%state%enth_t(kt,j,i) = enth_fct_temp_omega(ice%state%temp_t_m(kt,j,i), ice%state%omega_t(kt,j,i))
            end do
          else
            do kt=0, ice%grid%KTMAX
              ice%state%enth_t(kt,j,i) = ice%state%enth_c(0,j,i)
            end do
          end if
        end do
      end do

    else if (ice%par%calcmod==0 .or. ice%par%calcmod==-1) then

      do i=0, ice%grid%IMAX
        do j=0, ice%grid%JMAX
          do kc=0, ice%grid%KCMAX
            ice%state%enth_c(kc,j,i) = enth_fct_temp_omega(ice%state%temp_c(kc,j,i), 0.0_wp)
          end do
          do kt=0, ice%grid%KTMAX
            ice%state%enth_t(kt,j,i) = ice%state%enth_c(0,j,i)
          end do
        end do
      end do

    endif

    call nc_write(fnm,"enth_c", real(reshape(ice%state%enth_c,[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KCMAX+1],order=[3,2,1]),sp), &
          dims=[dim_x,dim_y,dim_kc,dim_time],start=[1,1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KCMAX+1,1], &
          long_name="Enthalpy in the upper (kc) ice layer",grid_mapping="polar_stereographic",units="J7kg",ncid=ncid)

    call nc_write(fnm,"enth_t", real(reshape(ice%state%enth_t,[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KTMAX+1],order=[3,2,1]),sp), &
          dims=[dim_x,dim_y,dim_kt,dim_time],start=[1,1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KTMAX+1,1], &
          long_name="Enthalpy in the lower (kt) ice layer",grid_mapping="polar_stereographic",units="J7kg",ncid=ncid)

    call nc_write(fnm,"enh_c", real(reshape(ice%state%enh_c,[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KCMAX+1],order=[3,2,1]),sp), &
          dims=[dim_x,dim_y,dim_kc,dim_time],start=[1,1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KCMAX+1,1], &
          long_name="Flow enhancement factor in the upper (kc) ice layer",grid_mapping="polar_stereographic",units="1",ncid=ncid)

    call nc_write(fnm,"enh_t", real(reshape(ice%state%enh_t,[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KTMAX+1],order=[3,2,1]),sp), &
          dims=[dim_x,dim_y,dim_kt,dim_time],start=[1,1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KTMAX+1,1], &
          long_name="Flow enhancement factor in the lower (kt) ice layer",grid_mapping="polar_stereographic",units="1",ncid=ncid)

    call nc_write(fnm,"age_c", real(reshape(ice%state%enh_c,[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KCMAX+1],order=[3,2,1])*sec_year,sp), &
          dims=[dim_x,dim_y,dim_kc,dim_time],start=[1,1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KCMAX+1,1], &
          long_name="Age in the upper (kc) ice layer",grid_mapping="polar_stereographic",units="a",ncid=ncid)

    call nc_write(fnm,"age_t", real(reshape(ice%state%age_t,[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KTMAX+1],order=[3,2,1])*sec_year,sp), &
          dims=[dim_x,dim_y,dim_kt,dim_time],start=[1,1,1,nout],count=[ice%grid%IMAX+1,ice%grid%JMAX+1,ice%grid%KTMAX+1,1], &
          long_name="Age in the lower (kt) ice layer",grid_mapping="polar_stereographic",units="a",ncid=ncid)

endif

!    ---- End of the definitions

!#if (DISC>0)   /* Ice discharge parameterisation */
!
!!    ---- dis_perp
!
!call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
!            thisroutine )
!call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
!            thisroutine )
!call check( nf90_def_var(ncid, 'dis_perp', NF90_FLOAT, nc2d, ncv), &
!            thisroutine )
!buffer = 'm a-1'
!call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
!            thisroutine )
!buffer = 'ice_discharge'
!call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
!            thisroutine )
!buffer = 'Ice discharge'
!call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
!            thisroutine )
!call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'crs'), &
!            thisroutine )
!call check( nf90_put_att(ncid, ncv, 'coordinates', 'lat lon'), &
!            thisroutine )
!
!!    ---- cst_dist
!
!call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
!            thisroutine )
!call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
!            thisroutine )
!call check( nf90_def_var(ncid, 'cst_dist', NF90_FLOAT, nc2d, ncv), &
!            thisroutine )
!buffer = 'km'
!call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
!            thisroutine )
!buffer = 'coastal_distance'
!call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
!            thisroutine )
!buffer = 'Coastal distance'
!call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
!            thisroutine )
!call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'crs'), &
!            thisroutine )
!call check( nf90_put_att(ncid, ncv, 'coordinates', 'lat lon'), &
!            thisroutine )
!
!!    ---- cos_grad_tc
!
!call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
!            thisroutine )
!call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
!            thisroutine )
!call check( nf90_def_var(ncid, 'cos_grad_tc', NF90_FLOAT, nc2d, ncv), &
!            thisroutine )
!buffer = '1'
!call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
!            thisroutine )
!buffer = 'cos_alpha'
!call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
!            thisroutine )
!buffer = 'Cosine of angle between surface gradient and cst dist gradient'
!call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
!            thisroutine )
!call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'crs'), &
!            thisroutine )
!call check( nf90_put_att(ncid, ncv, 'coordinates', 'lat lon'), &
!            thisroutine )
!
!!    ---- mask_mar
!
!call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
!            thisroutine )
!call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
!            thisroutine )
!call check( nf90_def_var(ncid, 'mask_mar', NF90_SHORT, nc2d, ncv), &
!            thisroutine )
!buffer = 'marginal_ring_mask'
!call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
!            thisroutine )
!buffer = 'Marginal ring mask'
!call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
!            thisroutine )
!nc2flag = (/ 0, 1 /)
!call check( nf90_put_att(ncid, ncv, 'flag_values', nc2flag), &
!            thisroutine )
!buffer = 'no_ring '// &
!         'ring'
!call check( nf90_put_att(ncid, ncv, 'flag_meanings', trim(buffer)), &
!            thisroutine )
!call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'crs'), &
!            thisroutine )
!call check( nf90_put_att(ncid, ncv, 'coordinates', 'lat lon'), &
!            thisroutine )
!
!#endif





! fixme, is this for restart??
!
!!-------- Convert data to real*4 and seconds to years --------
!
!#if (!defined(OUT_TIMES) || OUT_TIMES==1)
!time_conv = real(time*sec_to_year,sp)
!#elif (OUT_TIMES==2)
!time_conv = real((time+year_zero)*sec_to_year,sp)
!#else
!stop ' >>> output1: OUT_TIMES must be either 1 or 2!'
!#endif
!
!delta_ts_conv   = real(delta_ts,sp)
!glac_index_conv = real(glac_index,sp)
!V_tot_conv      = real(V_tot,sp)
!V_af_conv       = real(V_af,sp)
!A_grounded_conv = real(A_grounded,sp)
!A_floating_conv = real(A_floating,sp)
!H_R_conv        = real(H_R,sp)
!
!do i=0, ice%grid%IMAX
!   xi_conv(i) = real(xi(i),sp)
!end do
!
!do j=0, JMAX
!   eta_conv(j) = real(eta(j),sp)
!end do
!
!do kc=0, KCMAX
!   sigma_level_c_conv(kc) = real(eaz_c_quotient(kc),sp)
!end do
!
!do kt=0, KTMAX
!   sigma_level_t_conv(kt) = real(zeta_t(kt),sp)
!end do
!
!do kr=0, KRMAX
!   sigma_level_r_conv(kr) = real(kr,sp)/real(KRMAX,sp)
!end do
!
!do i=0, ice%grid%IMAX
!do j=0, JMAX
!
!   maske_conv(i,j)     = maske(j,i)
!   maske_old_conv(i,j) = maske_old(j,i)
!   n_cts_conv(i,j)     = n_cts(j,i)
!   kc_cts_conv(i,j)    = kc_cts(j,i)
!
!   lambda_conv(i,j)    = real(lambda(j,i),sp) 
!   phi_conv(i,j)       = real(phi(j,i),sp) 
!   lon_conv(i,j)       = real(lambda(j,i)*pi_180_inv,sp)   ! longitude in deg
!   lon_conv(i,j)       = modulo(lon_conv(i,j)+180.0_sp, 360.0_sp)-180.0_sp
!                                    ! mapping to interval [-180 deg, +180 deg)
!   lat_conv(i,j)       = real(phi(j,i)   *pi_180_inv,sp)   ! latitute  in deg
!   if (lat_conv(i,j) >  90.0_sp) lat_conv(i,j) =  90.0_sp
!   if (lat_conv(i,j) < -90.0_sp) lat_conv(i,j) = -90.0_sp
!                                 ! constraining to interval [-90 deg, +90 deg]
!   temp_s_conv(i,j)      = real(temp_s(j,i),sp)
!   accum_conv(i,j)       = real(accum(j,i)*sec_year,sp)
!   as_perp_conv(i,j)     = real(as_perp(j,i)*sec_year,sp)
!   as_perp_apl_conv(i,j) = real(as_perp_apl(j,i)*sec_year,sp)
!
!#if (DISC>0)   /* Ice discharge parameterisation */
!   dis_perp_conv(i,j)  = real(dis_perp(j,i)*sec_year,sp)
!   cst_dist_conv(i,j)  = real(cst_dist(j,i)*0.001_wp,sp)
!   cos_grad_tc_conv(i,j) = real(cos_grad_tc(j,i),sp)
!   mask_mar_conv(i,j)  = mask_mar(j,i)
!#endif
!
!   q_geo_conv(i,j)     = real(q_geo(j,i),sp)
!   zs_conv(i,j)        = real(zs(j,i),sp)
!   zm_conv(i,j)        = real(zm(j,i),sp)
!   zb_conv(i,j)        = real(zb(j,i),sp)
!   zl_conv(i,j)        = real(zl(j,i),sp)
!   zl0_conv(i,j)       = real(zl0(j,i),sp)
!   H_cold_conv(i,j)    = real(H_cold(j,i),sp)
!   H_temp_conv(i,j)    = real(H_temp(j,i),sp)
!   H_conv(i,j)         = real(H(j,i),sp)
!   Q_bm_conv(i,j)      = real(Q_bm(j,i)*sec_year,sp)
!   Q_tld_conv(i,j)     = real(Q_tld(j,i)*sec_year,sp)
!   am_perp_conv(i,j)   = real(am_perp(j,i)*sec_year,sp)
!   qx_conv(i,j)        = real(qx(j,i)*sec_year,sp)
!   qy_conv(i,j)        = real(qy(j,i)*sec_year,sp)
!   dzs_dtau_conv(i,j)  = real(dzs_dtau(j,i)*sec_year,sp)
!   dzm_dtau_conv(i,j)  = real(dzm_dtau(j,i)*sec_year,sp)
!   dzb_dtau_conv(i,j)  = real(dzb_dtau(j,i)*sec_year,sp)
!   dzl_dtau_conv(i,j)  = real(dzl_dtau(j,i)*sec_year,sp)
!   dH_c_dtau_conv(i,j) = real(dH_c_dtau(j,i)*sec_year,sp)
!   dH_t_dtau_conv(i,j) = real(dH_t_dtau(j,i)*sec_year,sp)
!   dH_dtau_conv(i,j)   = real(dH_dtau(j,i)*sec_year,sp)
!   vx_b_g_conv(i,j)    = real(vx_b_g(j,i)*sec_year,sp)
!   vy_b_g_conv(i,j)    = real(vy_b_g(j,i)*sec_year,sp)
!   vz_b_conv(i,j)      = real(vz_b(j,i)*sec_year,sp)
!   vh_b_conv(i,j)      = sqrt( vx_b_g_conv(i,j)**2 + vy_b_g_conv(i,j)**2 )
!   vx_s_g_conv(i,j)    = real(vx_s_g(j,i)*sec_year,sp)
!   vy_s_g_conv(i,j)    = real(vy_s_g(j,i)*sec_year,sp)
!   vz_s_conv(i,j)      = real(vz_s(j,i)*sec_year,sp)
!   vh_s_conv(i,j)      = sqrt( vx_s_g_conv(i,j)**2 + vy_s_g_conv(i,j)**2 )
!   vx_m_g_conv(i,j)    = real(vx_m_g(j,i)*sec_year,sp)
!   vy_m_g_conv(i,j)    = real(vy_m_g(j,i)*sec_year,sp)
!   vh_m_conv(i,j)      = sqrt( vx_m_g_conv(i,j)**2 + vy_m_g_conv(i,j)**2 )
!   temp_b_conv(i,j)    = real(temp_b(j,i),sp)
!   temph_b_conv(i,j)   = real(temph_b(j,i),sp)
!   tau_b_driving_conv(i,j) = real(tau_b_driving(j,i),sp)
!   tau_b_drag_conv(i,j)    = real(tau_b_drag(j,i),sp)
!   p_b_w_conv(i,j)     = real(p_b_w(j,i),sp)
!   H_w_conv(i,j)       = real(H_w(j,i),sp)
!   q_gl_g_conv(i,j)    = real(q_gl_g(j,i)*sec_year,sp)
!   q_cf_g_conv(i,j)    = real(calving(j,i)*sec_year,sp)
!   ratio_sl_x_conv(i,j) = real(ratio_sl_x(j,i),sp)
!   ratio_sl_y_conv(i,j) = real(ratio_sl_y(j,i),sp)
!
!   if (flag_shelfy_stream_x(j,i)) then
!      flag_shelfy_stream_x_conv(i,j) = 1
!   else
!      flag_shelfy_stream_x_conv(i,j) = 0
!   end if
!
!   if (flag_shelfy_stream_y(j,i)) then
!      flag_shelfy_stream_y_conv(i,j) = 1
!   else
!      flag_shelfy_stream_y_conv(i,j) = 0
!   end if
!
!   if (flag_shelfy_stream(j,i)) then
!      flag_shelfy_stream_conv(i,j) = 1
!   else
!      flag_shelfy_stream_conv(i,j) = 0
!   end if
!
!   vis_int_g_conv(i,j) = real(vis_int_g(j,i),sp)
!
!   do kr=0, KRMAX
!      temp_r_conv(i,j,kr) = real(temp_r(kr,j,i),sp)
!   end do
!
!   do kt=0, KTMAX
!      vx_t_conv(i,j,kt)    = real(vx_t(kt,j,i)*sec_year,sp)
!      vy_t_conv(i,j,kt)    = real(vy_t(kt,j,i)*sec_year,sp)
!      vz_t_conv(i,j,kt)    = real(vz_t(kt,j,i)*sec_year,sp)
!      omega_t_conv(i,j,kt) = real(omega_t(kt,j,i),sp)
!      age_t_conv(i,j,kt)   = real(age_t(kt,j,i)*sec_to_year,sp)
!      enth_t_conv(i,j,kt)  = real(enth_t(kt,j,i),sp)
!      enh_t_conv(i,j,kt)   = real(enh_t(kt,j,i),sp)
!   end do
!
!   do kc=0, KCMAX
!      vx_c_conv(i,j,kc)    = real(vx_c(kc,j,i)*sec_year,sp)
!      vy_c_conv(i,j,kc)    = real(vy_c(kc,j,i)*sec_year,sp)
!      vz_c_conv(i,j,kc)    = real(vz_c(kc,j,i)*sec_year,sp)
!      temp_c_conv(i,j,kc)  = real(temp_c(kc,j,i),sp)
!      age_c_conv(i,j,kc)   = real(age_c(kc,j,i)*sec_to_year,sp)
!      enth_c_conv(i,j,kc)  = real(enth_c(kc,j,i),sp)
!      omega_c_conv(i,j,kc) = real(omega_c(kc,j,i),sp)
!      enh_c_conv(i,j,kc)   = real(enh_c(kc,j,i),sp)
!   end do
!
!end do
!end do

   return

  end subroutine ice_nc_write


end module sico_out
