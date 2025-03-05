module sicopolis

  use nml
  use ncio
  use coord, only : grid_class, grid_init
  use coord, only : map_scrip_class, map_scrip_init, map_scrip_field

  use timer, only : sec_year
  use control, only : restart_in_dir

  use sico_types_m, only : wp, sp

  use sico_state, only : sico_state_class
  use sico_params, only : sico_par_class, sico_params_init
  use sico_params, only : C, L, rho_w, rho_sw, rho
  use sico_grid_mod, only : sico_grid_class, sico_grid_init
  use sico_timer, only : sico_timer_class, sico_timer_init
  use sico_check, only : check_vel

  use boundary_m, only : boundary

  use init_temp_water_age_m, only : init_temp_water_age_1_1, init_temp_water_age_1_2, init_temp_water_age_1_3, init_temp_water_age_1_4
  use flag_update_gf_gl_cf_m, only : flag_update_gf_gl_cf
  use topograd_m, only : topograd_1, topograd_2
  use enth_temp_omega_m, only : enth_fct_temp_omega, calc_c_int_table, calc_c_int_inv_table
  use calc_temp_m, only : calc_temp_const, calc_temp_cold, calc_temp_poly
  use calc_temp_enth_m, only : calc_temp_enth
  use calc_enhance_m, only : calc_enhance
  use calc_vxy_m, only : calc_dzs_dxy_aux, calc_vxy_static, calc_vxy_b_sia, calc_vxy_sia, calc_vxy_ssa
  use calc_vz_m, only : calc_vz_static, calc_vz_grounded, calc_vz_floating
  use calc_dxyz_m, only : calc_dxyz
  use calc_thk_m, only : calc_thk_init, calc_thk_sia_expl, calc_thk_sia_impl_sor, calc_thk_sia_impl_lis
  use calc_thk_m, only : calc_thk_expl, calc_thk_impl_sor, calc_thk_impl_lis
  use calc_thk_m, only : account_mb_source, calc_thk_mask_update
  use calc_temp_melt_bas_m, only : calc_temp_melt, calc_temp_bas
  use calc_bas_melt_m, only : calc_qbm
  use calc_thk_water_bas_m, only : calc_thk_water_bas

  !$  use omp_lib

  implicit none

  type ts_out_class
    real(wp) :: V_tot, V_af, V_sle, V_grounded, V_floating, V_gr_redu  
    real(wp) :: A_tot, A_grounded, A_floating 
    real(wp) :: V_temp    
    real(wp) :: A_temp    
    real(wp) :: Q_s       
    real(wp) :: Q_b       
    real(wp) :: Q_temp    
    real(wp) :: H_max     
    real(wp) :: H_t_max   
    real(wp) :: zs_max   
    real(wp) :: vs_max   
    real(wp) :: Tbh_max  
    real(wp) :: dV_dt      
    real(wp) :: precip_tot 
    real(wp) :: runoff_tot 
    real(wp) :: calv_tot  
    real(wp) :: disc_lsc 
    real(wp) :: disc_ssc   
    real(wp) :: mb_tot   
    real(wp) :: bmb_tot   
    real(wp) :: bmb_gr_tot   
    real(wp) :: bmb_fl_tot   
    real(wp) :: mbp
    real(wp) :: mb_resid
  end type

  type sico_class
    logical :: error
    ! grid definition
    type(sico_grid_class) :: grid
    ! parameter definition
    type(sico_par_class) :: par
    ! timer definition
    type(sico_timer_class) :: timer
    ! state definition 
    type(sico_state_class) :: state
    ! time series output 
    type(ts_out_class), allocatable :: ts(:)
    ! for output
    integer :: nout
  end type


  private
  public :: sico_class, sico_init, sico_update, sico_end, sico_write_restart, ts_out_class

contains


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s i c o _ u p d a t e
  !   Purpose    :  update sicopolis
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sico_update(ice)

  implicit none

  type(sico_class) :: ice

  !$ real(wp) :: time1,time2


  ice%timer%itercount = ice%timer%itercount + 1

  if (.not.ice%par%l_restart .and. ice%timer%itercount==1) then   
    ! first sico_update call, initialize temperature profiles if surface temperature is needed!
    if (ice%par%i_temp_init==2) then
      ! In each ice column equal to the local ice surface temperature
      call init_temp_water_age_1_2(ice%state,ice%grid)
    else if (ice%par%i_temp_init==3) then
      ! Ice temperature linearly increasing with depth
      call init_temp_water_age_1_3(ice%state,ice%grid)
    else if (ice%par%i_temp_init==4) then
      ! Ice-temperature profiles from analytical solution for 1-d steady-state advection-diffusion equation
      ! under the assumption of linearly varying vz [Robin (1955) solution]
      call init_temp_water_age_1_4(ice%state,ice%grid)
    endif
  endif

!-------- Save old mask --------

  ice%state%maske_old = ice%state%maske

!-------- Boundary conditions --------

  call boundary(ice%state,ice%grid,ice%timer,ice%par)

!-------- Temperature, water content, age, flow enhancement factor --------

  if ( mod(ice%timer%itercount, ice%timer%iter_temp) == 0 ) then

!  ------ Temperature, water content, age

    !$ time1 = omp_get_wtime()
    if (ice%par%calcmod==1) then
     call calc_temp_poly(ice%state,ice%grid,ice%timer,ice%par) ! polythermal scheme (POLY)
    else if (ice%par%calcmod==0) then
     call calc_temp_cold(ice%state,ice%grid,ice%timer,ice%par) ! cold-ice scheme (COLD)
    else if (ice%par%calcmod==2 .or. ice%par%calcmod==3) then
     call calc_temp_enth(ice%state,ice%grid,ice%timer,ice%par) ! enthalpy scheme (ENTC or ENTM)
    else if (ice%par%calcmod==-1) then
     call calc_temp_const(ice%state,ice%grid,ice%timer,ice%par)   ! isothermal scheme (ISOT)
    else
     stop ' >>> sico_main_loop: Parameter CALCMOD must be between -1 and 3!'
    endif
    !$ time2 = omp_get_wtime()
    !$ if(ice%par%l_write_timer) print *,'calc_temp',time2-time1

!  ------ Time derivative of H_t (further time derivatives are
!         computed in subroutine calc_thk_xxx)

    !$ time1 = omp_get_wtime()
     ice%state%dH_t_dtau      = (ice%state%H_t_neu - ice%state%H_t)*ice%timer%dtime_temp_inv

!  ------ New values --> old values

     ice%state%n_cts   = ice%state%n_cts_neu
     ice%state%kc_cts  = ice%state%kc_cts_neu
     ice%state%zm      = ice%state%zm_neu
     ice%state%H_c     = ice%state%H_c_neu
     ice%state%H_t     = ice%state%H_t_neu
     ice%state%temp_c  = ice%state%temp_c_neu
     ice%state%age_c   = ice%state%age_c_neu
     ice%state%omega_t = ice%state%omega_t_neu
     ice%state%age_t   = ice%state%age_t_neu
     ice%state%temp_r  = ice%state%temp_r_neu

     if (ice%par%calcmod==2 .or. ice%par%calcmod==3) then
       ice%state%enth_c  = ice%state%enth_c_neu
       ice%state%enth_t  = ice%state%enth_t_neu
       ice%state%omega_c = ice%state%omega_c_neu
     endif
    !$ time2 = omp_get_wtime()
    !$ if(ice%par%l_write_timer) print *,'new->old',time2-time1

!  ------ Flow enhancement factor
     call calc_enhance(ice%state,ice%grid,ice%timer,ice%par)

  end if
  !     End of computation of temperature, water content, age and
  !     enhancement factor

!-------- Velocity --------

    !$ time1 = omp_get_wtime()
  call flag_update_gf_gl_cf(ice%state,ice%grid,ice%par)
    !$ time2 = omp_get_wtime()
    !$ if(ice%par%l_write_timer) print *,'updata_gf_gl',time2-time1

    !$ time1 = omp_get_wtime()
  if (ice%par%margin==3 .and. ice%par%gl_surf_grad==2) then
    ! one-sided gradients at the grounding line
    call calc_dzs_dxy_aux(ice%state,ice%grid)
  else
    ! use centered differences
    ice%state%dzs_dx_aux = ice%state%dzs_dxi
    ice%state%dzs_dy_aux = ice%state%dzs_deta
  endif
    !$ time2 = omp_get_wtime()
    !$ if(ice%par%l_write_timer) print *,'dzs_dxy_aux',time2-time1

  if (ice%par%dynamics==1 .or. ice%par%dynamics==2) then

    !$ time1 = omp_get_wtime()
    call calc_vxy_b_sia(ice%state,ice%grid,ice%par)
    !$ time2 = omp_get_wtime()
    !$ if(ice%par%l_write_timer) print *,'vxy_b_sia',time2-time1
    !$ time1 = omp_get_wtime()
    call calc_vxy_sia(ice%state,ice%grid,ice%par)
    !$ time2 = omp_get_wtime()
    !$ if(ice%par%l_write_timer) print *,'vxy_sia',time2-time1

    !$ time1 = omp_get_wtime()
    if (ice%par%margin==3 .or. ice%par%dynamics==2) then
      call calc_vxy_ssa(ice%state,ice%grid,ice%par)
    endif
    !$ time2 = omp_get_wtime()
    !$ if(ice%par%l_write_timer) print *,'vxy_ssa',time2-time1

    !$ time1 = omp_get_wtime()
    call calc_vz_grounded(ice%state,ice%grid)
    !$ time2 = omp_get_wtime()
    !$ if(ice%par%l_write_timer) print *,'vz_grounded',time2-time1

    !$ time1 = omp_get_wtime()
    if (ice%par%margin==3) then
      call calc_vz_floating(ice%state,ice%grid)
    endif
    !$ time2 = omp_get_wtime()
    !$ if(ice%par%l_write_timer) print *,'vz_floating',time2-time1

  else if (ice%par%dynamics==0) then

    call calc_vxy_static(ice%state,ice%par)
    call calc_vz_static(ice%state,ice%grid)

  else
    stop ' >>> sico_main_loop: DYNAMICS must be either 0, 1 or 2!'
  endif

  ! check velocities
  call check_vel(ice%state,ice%grid,ice%error)

    !$ time1 = omp_get_wtime()
  call calc_dxyz(ice%state,ice%grid,ice%par)
    !$ time2 = omp_get_wtime()
    !$ if(ice%par%l_write_timer) print *,'dxyz',time2-time1

!-------- ice topography --------

    !$ time1 = omp_get_wtime()
  call calc_thk_init(ice%state,ice%grid,ice%par)
    !$ time2 = omp_get_wtime()
    !$ if(ice%par%l_write_timer) print *,'thk_init',time2-time1

    !$ time1 = omp_get_wtime()
  if (ice%par%calcthk==1) then
    call calc_thk_sia_expl(ice%state,ice%grid,ice%timer,ice%par)
  else if (ice%par%calcthk==2) then
    call calc_thk_sia_impl_sor(ice%state,ice%grid,ice%timer,ice%par)
  else if (ice%par%calcthk==3) then
    call calc_thk_sia_impl_lis(ice%state,ice%grid,ice%timer,ice%par)
  else if (ice%par%calcthk==4) then
    call calc_thk_expl(ice%state,ice%grid,ice%timer,ice%par)
  else if (ice%par%calcthk==5) then
    call calc_thk_impl_sor(ice%state,ice%grid,ice%timer,ice%par)
  else if (ice%par%calcthk==6) then
    call calc_thk_impl_lis(ice%state,ice%grid,ice%timer,ice%par)
  endif
    !$ time2 = omp_get_wtime()
    !$ if(ice%par%l_write_timer) print *,'thk',time2-time1

    !$ time1 = omp_get_wtime()
  if (ice%par%margin==3) then ! /* coupled SIA/SSA or SIA/SStA/SSA dynamics */
    call calc_thk_mask_update(ice%state,ice%grid,ice%timer,ice%par, 3)
  else if (ice%par%dynamics==2) then ! /* hybrid SIA/SStA dynamics */
    call calc_thk_mask_update(ice%state,ice%grid,ice%timer,ice%par, 2)
  else                ! /* SIA-only dynamics */
    if (ice%par%calcthk==1 .or. ice%par%calcthk==2 .or. ice%par%calcthk==3) then
      call calc_thk_mask_update(ice%state,ice%grid,ice%timer,ice%par, 1)
    else if (ice%par%calcthk==4 .or. ice%par%calcthk==5 .or. ice%par%calcthk==6) then
      call calc_thk_mask_update(ice%state,ice%grid,ice%timer,ice%par, 2) 
    endif
  endif
    !$ time2 = omp_get_wtime()
    !$ if(ice%par%l_write_timer) print *,'thk_mask',time2-time1

  call account_mb_source(ice%state,ice%grid,ice%timer)

  call flag_update_gf_gl_cf(ice%state,ice%grid,ice%par)

!  ------ New values --> old values

  ice%state%zs = ice%state%zs_neu
  ice%state%zm = ice%state%zm_neu
  ice%state%zb = ice%state%zb_neu
  ice%state%zl = ice%state%zl_neu
  ice%state%H_c= ice%state%H_c_neu
  ice%state%H_t= ice%state%H_t_neu
  ice%state%dzl_dtau = 0._wp   ! fixme, needed?

!-------- Melting temperature --------

  call calc_temp_melt(ice%state,ice%grid)

!-------- Basal temperature --------

  call calc_temp_bas(ice%state,ice%grid)

!-------- Basal melting rate --------

  call calc_qbm(ice%state,ice%grid,ice%timer,ice%par)

!-------- Effective thickness of subglacial water  --------

  call calc_thk_water_bas(ice%state,ice%par,ice%timer)


  end subroutine sico_update
      

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s i c o _ i n i t
  !   Purpose    :  initialize sicopolis
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sico_init(ice, grid, l_restart, &
      z_bed, z_bed_fil, z_bed_rel, h_ice, q_geo, h_sed, id_mask)

! Include header for lis solver fortran interface
#include "lisf.h"

  implicit none

  type(sico_class) :: ice
  type(grid_class), intent(in) :: grid
  logical, intent(in) :: l_restart
  real(wp), dimension(:,:), intent(in) :: z_bed    ! bedrock elevation, used only if not from restart
  real(wp), dimension(:,:), intent(in) :: z_bed_fil    ! bedrock elevation, filtered, used only if not from restart
  real(wp), dimension(:,:), intent(in) :: z_bed_rel    ! icefree relaxed bedrock elevation
  real(wp), dimension(:,:), intent(in) :: h_ice    ! ice thickness, used only if not from restart
  real(wp), dimension(:,:), intent(in) :: q_geo    ! geothermal heat flux, used only if not from restart
  real(wp), dimension(:,:), intent(in) :: h_sed    ! sediment thickness, used only if not from restart
  integer, dimension(:,:), intent(in) :: id_mask   ! ice id mask

  integer :: i, j, ii, jj, kc, kt, i_f, j_f, n_filter, ni, nj
  integer :: ppos, spos
  integer :: ierr
  real(wp) :: sigma_filter, sigma_ref
  real(wp) :: dist, weigh, sum_weigh
  real(wp) :: Hice, freeboard_ratio
  type(grid_class) :: mask_maxextent_grid
  type(map_scrip_class) :: maps_maxextent_to_ice
  real(wp), dimension(:,:), allocatable :: tmp1
  integer, dimension(:,:), allocatable :: tmpi
  integer, dimension(:,:), allocatable :: maxi
  real(wp), dimension(:), allocatable :: lon_maxi, lat_maxi


  ! initialize parameters
  call sico_params_init(ice%par)

  ! initialize timer
  call sico_timer_init(ice%par%dtime, ice%par%dtime_temp, ice%par%dtime_mar_coa, ice%timer)

  ! initialize grid
  call sico_grid_init(grid, ice%grid,ice%timer,ice%par)

  ! allocate variables
  call sico_alloc(ice)

  !-------- Initialisation of the Library of Iterative Solvers Lis, if required -------- 
  if (ice%par%calcthk==3 .or. ice%par%calcthk==6 .or. ice%par%margin==3 .or. ice%par%dynamics==2) then
    call lis_initialize(ierr)
  endif

  !  ------ Some auxiliary quantities required for the enthalpy method

  call calc_c_int_table(C, -190, 10, real(L,wp))
  call calc_c_int_inv_table()

  !-------- Mean accumulation --------
  ice%state%mean_accum = ice%par%mean_accum*(1.0e-03_wp/sec_year)*(RHO_W/RHO) ! mm/a water equiv. -> m/s ice equiv.

  ! read maximum ice extent mask
  if (ice%par%thk_evol==4) then

    ! read from file
    ni = nc_size(trim(ice%par%mask_maxextent_file),"lon")
    nj = nc_size(trim(ice%par%mask_maxextent_file),"lat")
    allocate( maxi(ni,nj) )
    allocate( lon_maxi(ni) )
    allocate( lat_maxi(nj) )
    call nc_read(trim(ice%par%mask_maxextent_file),"lat",lat_maxi)
    call nc_read(trim(ice%par%mask_maxextent_file),"lon",lon_maxi)
    call nc_read(trim(ice%par%mask_maxextent_file),"mask",maxi)

    ! grid definition
    spos = scan(trim(ice%par%mask_maxextent_file),"/", BACK= .true.)+1
    ppos = scan(trim(ice%par%mask_maxextent_file),".", BACK= .true.)-1
    call grid_init(mask_maxextent_grid,name=trim(ice%par%mask_maxextent_file(spos:ppos)),mtype="latlon",units="degrees",x=lon_maxi,y=lat_maxi)
    ! map to ice grid
    allocate(tmpi(1:ice%grid%IMAX+1,1:ice%grid%JMAX+1))
    call map_scrip_init(maps_maxextent_to_ice,mask_maxextent_grid,ice%grid%grid1,method="nn",fldr="maps",load=.TRUE.,clean=.FALSE.)
    call map_scrip_field(maps_maxextent_to_ice,"mask",maxi,tmpi,method="mean",missing_value=-9999._wp, &
      filt_method="none",filt_par=[5.0*ice%grid%grid1%G%dx,ice%grid%grid1%G%dx])
    ice%state%mask_maxextent = transpose(tmpi)

    deallocate(maxi, lon_maxi, lat_maxi, tmpi)

  else

    ice%state%mask_maxextent(:,:) = 1

  endif

  ! set ice id mask 
  ice%state%id_mask = transpose(id_mask)


  !-------- Definition of initial values --------

  ice%par%l_restart = l_restart     ! save restart flag

  if (l_restart) then
    ! Read initial state from output of previous simulation

    call sico_read_restart(trim(restart_in_dir)//"/ice_sico_"//trim(ice%grid%grid1%name)//"_restart.nc",ice)
    print *,'read restart file ',trim(restart_in_dir)//"/ice_sico_"//trim(ice%grid%grid1%name)//"_restart.nc"

    if (ice%par%topograd==0) then
      call topograd_1(ice%state,ice%grid, 1)
    else if (ice%par%topograd==1) then
      call topograd_2(ice%state,ice%grid, 1)
    endif

    call boundary(ice%state, ice%grid, ice%timer, ice%par)

    where ((ice%state%maske==0).or.(ice%state%maske==3))
      ! grounded or floating ice
      ice%state%as_perp_apl = ice%state%as_perp
    elsewhere        ! maske==1 or 2, ice-free land or sea
      ice%state%as_perp_apl = 0.0_wp
    end where

    ice%state%Q_b_tot = ice%state%Q_bm + ice%state%Q_tld

  else 
    !  new start 

    allocate(tmp1(0:ice%grid%JMAX,0:ice%grid%IMAX))

    ! map input lithosphere elevation to ice grid
    ice%state%zl = transpose(z_bed)
    ice%state%zl_fil = transpose(z_bed_fil)
    ! filter lithosphere topography 
    if(ice%par%l_smooth_zl == 1) then
      tmp1 = ice%state%zl
      sigma_filter = ice%par%sigma_filter_zl/ice%grid%DX   ! half span of filtered area, in grid points
      n_filter     = ceiling(2.0_wp*sigma_filter)
      do i=0,ice%grid%IMAX
      do j=0,ice%grid%JMAX
        sum_weigh = 0.0_wp
        ice%state%zl(j,i) = 0._wp
        do ii=-n_filter, n_filter
        do jj=-n_filter, n_filter
          i_f = i+ii
          j_f = j+jj
          if (i_f <  0) i_f = 0
          if (i_f > ice%grid%IMAX) i_f = ice%grid%IMAX
          if (j_f <  0) j_f = 0
          if (j_f > ice%grid%JMAX) j_f = ice%grid%JMAX
          dist      = sqrt(real(ii,wp)**2+real(jj,wp)**2)
          weigh     = exp(-(dist/sigma_filter)**2)
          sum_weigh = sum_weigh + weigh
          ice%state%zl(j,i) = ice%state%zl(j,i) + weigh*tmp1(j_f,i_f)
        end do
        end do
        ice%state%zl(j,i) = ice%state%zl(j,i)/sum_weigh
      end do
      end do
    else if(ice%par%l_smooth_zl == 2) then
      tmp1 = ice%state%zl
      sigma_ref = ice%par%sigma_filter_zl/ice%grid%DX
      do i=0,ice%grid%IMAX
      do j=0,ice%grid%JMAX
        ! half span of filtered area, in grid points
        sigma_filter = min(max(sigma_ref, sigma_ref-ice%par%s_filter_zl* &
                               (tmp1(j,i)-ice%par%zl_cont)), 500.0_wp/ice%grid%DX)
        n_filter     = ceiling(2.0_wp*sigma_filter)
        sum_weigh = 0.0_wp
        ice%state%zl(j,i) = 0._wp
        do ii=-n_filter, n_filter
        do jj=-n_filter, n_filter
          i_f = i+ii
          j_f = j+jj
          if (i_f <  0) i_f = 0
          if (i_f > ice%grid%IMAX) i_f = ice%grid%IMAX
          if (j_f <  0) j_f = 0
          if (j_f > ice%grid%JMAX) j_f = ice%grid%JMAX
          dist      = sqrt(real(ii,wp)**2+real(jj,wp)**2)
          weigh     = exp(-(dist/sigma_filter)**2)
          sum_weigh = sum_weigh + weigh
          ice%state%zl(j,i) = ice%state%zl(j,i) + weigh*tmp1(j_f,i_f)
        end do
        end do
        ice%state%zl(j,i) = ice%state%zl(j,i)/sum_weigh
      end do
      end do
    end if 

    ! map input ice thickness to ice grid
    ice%state%H_c = transpose(h_ice)
    ! no ice at borders of domain
    ice%state%H_c(0:1,:) = 0._wp
    ice%state%H_c(:,0:1) = 0._wp
    ice%state%H_c(ice%grid%JMAX-1:ice%grid%JMAX,:) = 0._wp
    ice%state%H_c(:,ice%grid%IMAX-1:ice%grid%IMAX) = 0._wp

    ! map input relaxed lithosphere elevation to ice grid
    ice%state%zl0 = transpose(z_bed_rel)
    ! filter lithosphere topography
    if(ice%par%l_smooth_zl0 == 1) then
      tmp1 = ice%state%zl0
      sigma_filter = ice%par%sigma_filter_zl0/ice%grid%DX   ! half span of filtered area, in grid points
      n_filter     = ceiling(2.0_wp*sigma_filter)
      do i=0,ice%grid%IMAX
      do j=0,ice%grid%JMAX
        sum_weigh = 0.0_wp
        ice%state%zl0(j,i) = 0._wp
        do ii=-n_filter, n_filter
        do jj=-n_filter, n_filter
          i_f = i+ii
          j_f = j+jj
          if (i_f <  0) i_f = 0
          if (i_f > ice%grid%IMAX) i_f = ice%grid%IMAX
          if (j_f <  0) j_f = 0
          if (j_f > ice%grid%JMAX) j_f = ice%grid%JMAX
          dist      = sqrt(real(ii,wp)**2+real(jj,wp)**2)
          weigh     = exp(-(dist/sigma_filter)**2)
          sum_weigh = sum_weigh + weigh
          ice%state%zl0(j,i) = ice%state%zl0(j,i) + weigh*tmp1(j_f,i_f)
        end do
        end do
        ice%state%zl0(j,i) = ice%state%zl0(j,i)/sum_weigh
      end do
      end do
    else if(ice%par%l_smooth_zl0 == 2) then
      tmp1 = ice%state%zl0
      sigma_ref = ice%par%sigma_filter_zl0/ice%grid%DX
      do i=0,ice%grid%IMAX
      do j=0,ice%grid%JMAX
        ! half span of filtered area, in grid points
        sigma_filter = min(max(sigma_ref, sigma_ref-ice%par%s_filter_zl0* &
                               (tmp1(j,i)-ice%par%zl0_cont)), 500.0_wp/ice%grid%DX)
        n_filter     = ceiling(2.0_wp*sigma_filter)
        sum_weigh = 0.0_wp
        ice%state%zl0(j,i) = 0._wp
        do ii=-n_filter, n_filter
        do jj=-n_filter, n_filter
          i_f = i+ii
          j_f = j+jj
          if (i_f <  0) i_f = 0
          if (i_f > ice%grid%IMAX) i_f = ice%grid%IMAX
          if (j_f <  0) j_f = 0
          if (j_f > ice%grid%JMAX) j_f = ice%grid%JMAX
          dist      = sqrt(real(ii,wp)**2+real(jj,wp)**2)
          weigh     = exp(-(dist/sigma_filter)**2)
          sum_weigh = sum_weigh + weigh
          ice%state%zl0(j,i) = ice%state%zl0(j,i) + weigh*tmp1(j_f,i_f)
        end do
        end do
        ice%state%zl0(j,i) = ice%state%zl0(j,i)/sum_weigh
      end do
      end do
    endif

    deallocate(tmp1)

    ! geothermal heat flux
    ice%state%q_geo = transpose(q_geo)

    ! sediment thickness 
    ice%state%h_sed = transpose(h_sed)

    ! ensure consistency between elevations and mask

    freeboard_ratio = (RHO_SW-RHO)/RHO_SW
    do i=0, ice%grid%IMAX
      do j=0, ice%grid%JMAX

        if (ice%state%H_c(j,i)>0._wp .and. ice%state%H_c(j,i)>RHO_SW/RHO*(-ice%state%zl(j,i))) then  ! grounded ice

          ice%state%maske(j,i) = 0
          ice%state%zb(j,i) = ice%state%zl(j,i)
          ice%state%zs(j,i) = ice%state%zl(j,i)+ice%state%H_c(j,i)

        else if (ice%state%H_c(j,i)>0._wp) then ! floating ice

          ice%state%maske(j,i) = 3

          if (ice%par%margin==1 .or. (ice%par%margin==2) .and. ice%par%marine_ice_formation==1) then
            ice%state%maske(j,i) = 2   ! floating ice cut off
            ice%state%zs(j,i) = ice%state%zl(j,i)
            ice%state%zb(j,i) = ice%state%zl(j,i)
          else if (ice%par%margin==2 .and. ice%par%marine_ice_formation==2) then
            ice%state%maske(j,i) = 0   ! floating ice becomes "underwater ice"
            Hice   = ice%state%H_c(j,i)   ! ice thickness
            if (Hice < 0.0_wp) then
              Hice = 0.0_wp
              ice%state%maske(j,i) = 2
            end if
            ice%state%zs(j,i) = ice%state%zl(j,i)+Hice
            ice%state%zb(j,i) = ice%state%zl(j,i)
          else if (ice%par%margin==3) then
            Hice = ice%state%H_c(j,i)   ! ice thickness
            if (Hice < 0.0_wp) then
              Hice = 0.0_wp
              ice%state%maske(j,i) = 2
            end if
            ice%state%zs(j,i) = freeboard_ratio*Hice     ! ensure properly
            ice%state%zb(j,i) = ice%state%zs(j,i)-Hice   ! floating ice
          endif

        else if (ice%state%zl(j,i)>0._wp) then ! ice-free land

          ice%state%maske(j,i) = 1
          ice%state%zb(j,i) = ice%state%zl(j,i)  
          ice%state%zs(j,i) = ice%state%zl(j,i)  

        else if (ice%state%zl(j,i)<=0._wp) then ! ocean

          ice%state%maske(j,i) = 2
          if (ice%par%margin==1 .or. ice%par%margin==2) then
            ice%state%zs(j,i) = ice%state%zl(j,i)
            ice%state%zb(j,i) = ice%state%zl(j,i)
          else if (ice%par%margin==3) then
            ice%state%zs(j,i) = 0.0_wp   ! sea level
            ice%state%zb(j,i) = 0.0_wp   ! sea level
          endif

        else

          stop 'error in initial ice'

        endif

        ice%state%zm(j,i) = ice%state%zb(j,i)
        ice%state%n_cts(j,i) = -1
        ice%state%kc_cts(j,i) = 0

        ice%state%H_c(j,i) = ice%state%zs(j,i)-ice%state%zm(j,i)
        ice%state%H_t(j,i) = 0.0_wp

        ice%state%dzs_dtau(j,i)  = 0.0_wp
        ice%state%dzm_dtau(j,i)  = 0.0_wp
        ice%state%dzb_dtau(j,i)  = 0.0_wp
        ice%state%dzl_dtau(j,i)  = 0.0_wp
        ice%state%dH_c_dtau(j,i) = 0.0_wp
        ice%state%dH_t_dtau(j,i) = 0.0_wp

      end do
    end do

    ice%state%maske_old = ice%state%maske

    if (ice%par%topograd==0) then
      call topograd_1(ice%state,ice%grid, 1)
    else if (ice%par%topograd==1) then
      call topograd_2(ice%state,ice%grid, 1)
    endif

    call boundary(ice%state,ice%grid,ice%timer,ice%par)

    ice%state%as_perp_apl = 0.0_wp
    ice%state%am_perp = 0.0_wp

    ice%state%smb_corr = 0.0_wp

    ice%state%calving_apl = 0._wp

    ice%state%Q_bm    = 0.0_wp
    ice%state%Q_tld   = 0.0_wp
    ice%state%Q_b_tot = 0.0_wp 
    ice%state%Q_b_apl = 0.0_wp 

    ice%state%p_b_w   = 0.0_wp
    ice%state%H_w     = 0.0_wp

    ! initialize temperature, water and age profiles
    ! constant temperature in the entire ice sheet
    call init_temp_water_age_1_1(ice%state,ice%grid,ice%par)

    call calc_enhance(ice%state,ice%grid,ice%timer,ice%par)

    ice%state%vx_c = 0._wp
    ice%state%vy_c = 0._wp
    ice%state%vx_t = 0._wp
    ice%state%vy_t = 0._wp

    ice%state%vx_m_ssa = 0.0_wp
    ice%state%vy_m_ssa = 0.0_wp

    ice%state%qx = 0._wp
    ice%state%qy = 0._wp

    ice%state%vis_ave_g = ice%par%visc_init_ssa

  endif

  !-------- Grounding line and calving front flags --------

  ice%state%flag_grounding_line_1 = .false.
  ice%state%flag_grounding_line_2 = .false.

  ice%state%flag_calving_front_1  = .false.
  ice%state%flag_calving_front_2  = .false.

  if (ice%par%margin==3) then

    do i=1, ice%grid%IMAX-1
      do j=1, ice%grid%JMAX-1

        if ( (ice%state%maske(j,i)==0) &   ! grounded ice
          .and. &
          (    (ice%state%maske(j,i+1)==3)   &   ! with
          .or.(ice%state%maske(j,i-1)==3)   &   ! one
          .or.(ice%state%maske(j+1,i)==3)   &   ! neighbouring
          .or.(ice%state%maske(j-1,i)==3) ) &   ! floating ice point
          ) &
          ice%state%flag_grounding_line_1(j,i) = .true.

        if ( (ice%state%maske(j,i)==3) &   ! floating ice
          .and. &
          (    (ice%state%maske(j,i+1)==0)   &   ! with
          .or.(ice%state%maske(j,i-1)==0)   &   ! one
          .or.(ice%state%maske(j+1,i)==0)   &   ! neighbouring
          .or.(ice%state%maske(j-1,i)==0) ) &   ! grounded ice point
          ) &
          ice%state%flag_grounding_line_2(j,i) = .true.

        if ( (ice%state%maske(j,i)==3) &   ! floating ice
          .and. &
          (    (ice%state%maske(j,i+1)==2)   &   ! with
          .or.(ice%state%maske(j,i-1)==2)   &   ! one
          .or.(ice%state%maske(j+1,i)==2)   &   ! neighbouring
          .or.(ice%state%maske(j-1,i)==2) ) &   ! ocean point
          ) &
          ice%state%flag_calving_front_1(j,i) = .true.

        if ( (ice%state%maske(j,i)==2) &   ! ocean
          .and. &
          (    (ice%state%maske(j,i+1)==3)   &   ! with
          .or.(ice%state%maske(j,i-1)==3)   &   ! one
          .or.(ice%state%maske(j+1,i)==3)   &   ! neighbouring
          .or.(ice%state%maske(j-1,i)==3) ) &   ! floating ice point
          ) &
          ice%state%flag_calving_front_2(j,i) = .true.

      end do
    end do

    do i=1, ice%grid%IMAX-1

      j=0
      if ( (ice%state%maske(j,i)==2) &   ! ocean
        .and. (ice%state%maske(j+1,i)==3) &   ! with one neighbouring
        ) &                               ! floating ice point
        ice%state%flag_calving_front_2(j,i) = .true.

      j=ice%grid%JMAX
      if ( (ice%state%maske(j,i)==2) &   ! ocean
        .and. (ice%state%maske(j-1,i)==3) &   ! with one neighbouring
        ) &                               ! floating ice point
        ice%state%flag_calving_front_2(j,i) = .true.

    end do

    do j=1, ice%grid%JMAX-1

      i=0
      if ( (ice%state%maske(j,i)==2) &   ! ocean
        .and. (ice%state%maske(j,i+1)==3) &   ! with one neighbouring
        ) &                               ! floating ice point
        ice%state%flag_calving_front_2(j,i) = .true.

      i=ice%grid%IMAX
      if ( (ice%state%maske(j,i)==2) &   ! ocean
        .and. (ice%state%maske(j,i-1)==3) &   ! with one neighbouring
        ) &                               ! floating ice point
        ice%state%flag_calving_front_2(j,i) = .true.

    end do

  endif

  !-------- Initial velocities --------

  call calc_temp_melt(ice%state,ice%grid)

  ice%state%vx_m = 0._wp
  ice%state%vy_m = 0._wp

  if (ice%par%margin==3 .and. ice%par%gl_surf_grad==2) then
    ! one-sided gradients at the grounding line
    call calc_dzs_dxy_aux(ice%state,ice%grid)
  else
    ! use centered differences
    ice%state%dzs_dx_aux = ice%state%dzs_dxi
    ice%state%dzs_dy_aux = ice%state%dzs_deta
  endif

  if (ice%par%dynamics==1 .or. ice%par%dynamics==2) then

    call calc_vxy_b_sia(ice%state,ice%grid,ice%par)
    call calc_vxy_sia(ice%state,ice%grid,ice%par)

    if (ice%par%margin==3 .or. ice%par%dynamics==2) then
      call calc_vxy_ssa(ice%state,ice%grid,ice%par)
    endif

    call calc_vz_grounded(ice%state,ice%grid)

    if (ice%par%margin==3) then
      call calc_vz_floating(ice%state,ice%grid)
    endif

  else if (ice%par%dynamics==0) then

    call calc_vxy_static(ice%state,ice%par)
    call calc_vz_static(ice%state,ice%grid)

  else
    stop ' >>> sico_init: DYNAMICS must be either 0, 1 or 2!'
  endif

  call calc_dxyz(ice%state,ice%grid,ice%par)

  !-------- Initial enthalpies --------

  if (ice%par%calcmod==0 .or. ice%par%calcmod==-1) then

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

  else if (ice%par%calcmod==1) then

    do i=0, ice%grid%IMAX
      do j=0, ice%grid%JMAX

        do kc=0, ice%grid%KCMAX
          ice%state%enth_c(kc,j,i) = enth_fct_temp_omega(ice%state%temp_c(kc,j,i), 0.0_wp)
        end do

        if ( (ice%state%maske(j,i) == 0).and.(ice%state%n_cts(j,i) == 1) ) then
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

  else if (ice%par%calcmod==2 .or. ice%par%calcmod==3) then

    do i=0, ice%grid%IMAX
      do j=0, ice%grid%JMAX

        do kc=0, ice%grid%KCMAX
          ice%state%enth_c(kc,j,i) = enth_fct_temp_omega(ice%state%temp_c(kc,j,i), ice%state%omega_c(kc,j,i))
        end do

        do kt=0, ice%grid%KTMAX
          ice%state%enth_t(kt,j,i) = ice%state%enth_c(0,j,i)
        end do

      end do
    end do

  else

    stop ' >>> sico_init: Parameter CALCMOD must be either -1, 0, 1, 2 or 3!'

  endif

  print*
  print*,'======================================================='
  print*,' Initialisation of ICE ', trim(ice%grid%grid1%name),' complete'
  print*,'======================================================='
  print*

  return

  end subroutine sico_init
      

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i c e _ a l l o c
  !   Purpose    :  allocate ice model
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sico_alloc(ice)

    implicit none 

    type(sico_class) :: ice


    ! allocate all state variables
   allocate(ice%state%maske                 (0:ice%grid%JMAX,0:ice%grid%IMAX))    
   allocate(ice%state%maske_old             (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%maske_neu             (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%id_mask               (0:ice%grid%JMAX,0:ice%grid%IMAX))    
   allocate(ice%state%H_sed                 (0:ice%grid%JMAX,0:ice%grid%IMAX))    
   allocate(ice%state%f_sed                 (0:ice%grid%JMAX,0:ice%grid%IMAX))    
   allocate(ice%state%mask_ablation_type    (0:ice%grid%JMAX,0:ice%grid%IMAX))    
   allocate(ice%state%n_cts                 (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%n_cts_neu             (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%kc_cts                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%kc_cts_neu            (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%mask_maxextent        (0:ice%grid%JMAX,0:ice%grid%IMAX))    
   allocate(ice%state%maske_target          (0:ice%grid%JMAX,0:ice%grid%IMAX))    
   allocate(ice%state%H_target              (0:ice%grid%JMAX,0:ice%grid%IMAX))     
   allocate(ice%state%flag_ocean_point      (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%flag_grounding_line_1 (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%flag_grounding_line_2 (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%flag_calving_front_1  (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%flag_calving_front_2  (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%flag_shelfy_stream_x  (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%flag_shelfy_stream_y  (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%flag_shelfy_stream    (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%flag_calc_vxy_ssa_x   (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%flag_calc_vxy_ssa_y   (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%flag_grounded_front_a_1(0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%flag_grounded_front_a_2(0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%flag_grounded_front_b_1(0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%flag_grounded_front_b_2(0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%z_sl                  (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%zs                    (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%zs_std                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%zl_std                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%zm                    (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%zb                    (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%zl                    (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%zl_fil                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%zl0                   (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%H_c                   (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%H_t                   (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%H                     (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%H_neu                 (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%H_neu_flow            (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%H_neu_tmp             (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%H_smb                 (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%H_adv                 (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%H_eff                 (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%H_calv                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%H_sea_neu             (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%H_balance             (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dzs_dxi               (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dzs_dx_aux            (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dzm_dxi               (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dzb_dxi               (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dH_c_dxi              (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dH_t_dxi              (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dzs_deta              (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dzs_dy_aux            (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dzm_deta              (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dzb_deta              (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dH_c_deta             (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dH_t_deta             (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dzs_dxi_g             (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dzm_dxi_g             (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dzb_dxi_g             (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dH_c_dxi_g            (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dH_t_dxi_g            (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dzs_deta_g            (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dzm_deta_g            (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dzb_deta_g            (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dH_c_deta_g           (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dH_t_deta_g           (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dzs_dtau              (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dzm_dtau              (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dzb_dtau              (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dzl_dtau              (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dH_c_dtau             (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dH_t_dtau             (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dH_t_smooth           (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%p_weert               (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%q_weert               (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%p_weert_inv           (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%c_slide               (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%c_drag                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%c_fric                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vb_t                  (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%delta                 (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%p_b_w                 (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vx_b                  (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vy_b                  (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vx_m                  (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vy_m                  (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vx_m_1                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vy_m_1                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vx_m_2                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vy_m_2                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vx_m_sia              (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vy_m_sia              (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vx_m_ssa              (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vy_m_ssa              (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vx_m_prev             (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vy_m_prev             (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%d_help_b              (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%ratio_sl              (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%ratio_sl_x            (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%ratio_sl_y            (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vx_b_g                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vy_b_g                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vz_b                  (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vz_m                  (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vx_s_g                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vy_s_g                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vz_s                  (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vz_sl                 (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%flui_ave_sia          (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%h_diff                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%qx                    (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%qy                    (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%q_gl_g                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%q_geo                 (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%temp_b                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%temph_b               (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%Q_bm                  (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%Q_bm_float            (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%Q_tld                 (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%Q_b_tot               (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%Q_b_apl               (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%H_w                   (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%q_w                   (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%q_w_x                 (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%q_w_y                 (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%runoff                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%runoff_apl            (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%as_perp               (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%as_perp_apl           (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%mb_source_apl         (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%mb_source             (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%smb_corr              (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%calving               (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%calving_apl           (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%mask_mar              (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%cst_dist              (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%i_cst                 (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%j_cst                 (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%cos_grad_tc           (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%dis_perp              (0:ice%grid%JMAX,0:ice%grid%IMAX)) 
   allocate(ice%state%t_ocn                 (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%s_ocn                 (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%accum                 (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%accum_apl             (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%temp_s                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%temp_g                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%am_perp               (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%am_perp_st            (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%zs_neu                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%zm_neu                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%zb_neu                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%zl_neu                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%zs_aux                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%zm_aux                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%zb_aux                (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%H_c_neu               (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%H_t_neu               (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vx_c      (0:ice%grid%KCMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vy_c      (0:ice%grid%KCMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vz_c      (0:ice%grid%KCMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%temp_c    (0:ice%grid%KCMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%temp_c_neu(0:ice%grid%KCMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%temp_c_m  (0:ice%grid%KCMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%temp_c_help(0:ice%grid%KCMAX))
   allocate(ice%state%age_c     (0:ice%grid%KCMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%age_c_neu (0:ice%grid%KCMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%sigma_c   (0:ice%grid%KCMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%txz_c   (0:ice%grid%KCMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%tyz_c   (0:ice%grid%KCMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%enh_c     (0:ice%grid%KCMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%de_ssa    (0:ice%grid%JMAX,0:ice%grid%IMAX)) 
   allocate(ice%state%vis_int_g (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vis_ave_g (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vis_ave_g_smooth (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vis_int_sgxy (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%beta_drag (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vx_g      (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vy_g      (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vx_t           (0:ice%grid%KTMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vy_t           (0:ice%grid%KTMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%vz_t           (0:ice%grid%KTMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%omega_t        (0:ice%grid%KTMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%omega_t_neu    (0:ice%grid%KTMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%temp_t_m       (0:ice%grid%KTMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%age_t          (0:ice%grid%KTMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%age_t_neu      (0:ice%grid%KTMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%sigma_t        (0:ice%grid%KTMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%txz_t        (0:ice%grid%KTMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%tyz_t        (0:ice%grid%KTMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%enh_t          (0:ice%grid%KTMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%temp_r         (0:ice%grid%KRMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%temp_r_neu     (0:ice%grid%KRMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%enth_c         (0:ice%grid%KCMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%enth_c_neu     (0:ice%grid%KCMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%omega_c        (0:ice%grid%KCMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%omega_c_neu    (0:ice%grid%KCMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%enth_t         (0:ice%grid%KTMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%enth_t_neu     (0:ice%grid%KTMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%de_c           (0:ice%grid%KCMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%lambda_shear_c (0:ice%grid%KCMAX,0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%de_t           (0:ice%grid%KTMAX,0:ice%grid%JMAX,0:ice%grid%IMAX)) 
   allocate(ice%state%lambda_shear_t (0:ice%grid%KTMAX,0:ice%grid%JMAX,0:ice%grid%IMAX)) 
   allocate(ice%state%check_point    (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%weigh_ssta_sia_x(0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%weigh_ssta_sia_y(0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%p_b            (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%p_b_red        (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%p_b_red_lim    (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%tau_b          (0:ice%grid%JMAX,0:ice%grid%IMAX))

   allocate(ice%state%gamma_slide_inv    (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%sub_melt_flag      (0:ice%grid%JMAX,0:ice%grid%IMAX))

   allocate(ice%state%upH_x_1            (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%upH_x_2            (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%upH_y_1            (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%upH_y_2            (0:ice%grid%JMAX,0:ice%grid%IMAX))

   allocate(ice%state%czs2            (0:ice%grid%JMAX,0:ice%grid%IMAX))
   allocate(ice%state%czs3            (0:ice%grid%JMAX,0:ice%grid%IMAX))

   return

  end subroutine sico_alloc 

     
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s i c o _ end
  !   Purpose    :  end sicopolis
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sico_end(ice)

    implicit none

    type(sico_class) :: ice


    ! Deallocate all state variables to free memory
    call sico_dealloc(ice)
 

   return

  end subroutine sico_end


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i c e _ d e a l l o c
  !   Purpose    :  deallocate ice model
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sico_dealloc(ice)

    implicit none 

    type(sico_class) :: ice


    ! Deallocate all state variables to free memory
   deallocate(ice%state%maske                 )    
   deallocate(ice%state%maske_old             )
   deallocate(ice%state%maske_neu             )
   deallocate(ice%state%id_mask               )    
   deallocate(ice%state%H_sed                 )    
   deallocate(ice%state%f_sed                 )    
   deallocate(ice%state%mask_ablation_type    )    
   deallocate(ice%state%n_cts                 )
   deallocate(ice%state%n_cts_neu             )
   deallocate(ice%state%kc_cts                )
   deallocate(ice%state%kc_cts_neu            )
   deallocate(ice%state%mask_maxextent        )    
   deallocate(ice%state%maske_target          )    
   deallocate(ice%state%H_target              )     
   deallocate(ice%state%flag_ocean_point      )
   deallocate(ice%state%flag_grounding_line_1 )
   deallocate(ice%state%flag_grounding_line_2 )
   deallocate(ice%state%flag_calving_front_1  )
   deallocate(ice%state%flag_calving_front_2  )
   deallocate(ice%state%flag_shelfy_stream_x  )
   deallocate(ice%state%flag_shelfy_stream_y  )
   deallocate(ice%state%flag_shelfy_stream    )
   deallocate(ice%state%flag_calc_vxy_ssa_x   )
   deallocate(ice%state%flag_calc_vxy_ssa_y   )
   deallocate(ice%state%flag_grounded_front_a_1)
   deallocate(ice%state%flag_grounded_front_a_2)
   deallocate(ice%state%flag_grounded_front_b_1)
   deallocate(ice%state%flag_grounded_front_b_2)
   deallocate(ice%state%z_sl                  )
   deallocate(ice%state%zs                    )
   deallocate(ice%state%zs_std                )
   deallocate(ice%state%zl_std                )
   deallocate(ice%state%zm                    )
   deallocate(ice%state%zb                    )
   deallocate(ice%state%zl                    )
   deallocate(ice%state%zl_fil                )
   deallocate(ice%state%zl0                   )
   deallocate(ice%state%H_c                   )
   deallocate(ice%state%H_t                   )
   deallocate(ice%state%H                     )
   deallocate(ice%state%H_neu                 )
   deallocate(ice%state%H_neu_flow            )
   deallocate(ice%state%H_neu_tmp             )
   deallocate(ice%state%H_smb                 )
   deallocate(ice%state%H_adv                 )
   deallocate(ice%state%H_eff                 )
   deallocate(ice%state%H_calv                )
   deallocate(ice%state%H_sea_neu             )
   deallocate(ice%state%H_balance             )
   deallocate(ice%state%dzs_dxi               )
   deallocate(ice%state%dzs_dx_aux            )
   deallocate(ice%state%dzm_dxi               )
   deallocate(ice%state%dzb_dxi               )
   deallocate(ice%state%dH_c_dxi              )
   deallocate(ice%state%dH_t_dxi              )
   deallocate(ice%state%dzs_deta              )
   deallocate(ice%state%dzs_dy_aux            )
   deallocate(ice%state%dzm_deta              )
   deallocate(ice%state%dzb_deta              )
   deallocate(ice%state%dH_c_deta             )
   deallocate(ice%state%dH_t_deta             )
   deallocate(ice%state%dzs_dxi_g             )
   deallocate(ice%state%dzm_dxi_g             )
   deallocate(ice%state%dzb_dxi_g             )
   deallocate(ice%state%dH_c_dxi_g            )
   deallocate(ice%state%dH_t_dxi_g            )
   deallocate(ice%state%dzs_deta_g            )
   deallocate(ice%state%dzm_deta_g            )
   deallocate(ice%state%dzb_deta_g            )
   deallocate(ice%state%dH_c_deta_g           )
   deallocate(ice%state%dH_t_deta_g           )
   deallocate(ice%state%dzs_dtau              )
   deallocate(ice%state%dzm_dtau              )
   deallocate(ice%state%dzb_dtau              )
   deallocate(ice%state%dzl_dtau              )
   deallocate(ice%state%dH_c_dtau             )
   deallocate(ice%state%dH_t_dtau             )
   deallocate(ice%state%dH_t_smooth           )
   deallocate(ice%state%p_weert               )
   deallocate(ice%state%q_weert               )
   deallocate(ice%state%p_weert_inv           )
   deallocate(ice%state%c_slide               )
   deallocate(ice%state%c_drag                )
   deallocate(ice%state%c_fric                )
   deallocate(ice%state%vb_t                  )
   deallocate(ice%state%delta                 )
   deallocate(ice%state%p_b_w                 )
   deallocate(ice%state%vx_b                  )
   deallocate(ice%state%vy_b                  )
   deallocate(ice%state%vx_m                  )
   deallocate(ice%state%vy_m                  )
   deallocate(ice%state%vx_m_1                )
   deallocate(ice%state%vy_m_1                )
   deallocate(ice%state%vx_m_2                )
   deallocate(ice%state%vy_m_2                )
   deallocate(ice%state%vx_m_sia              )
   deallocate(ice%state%vy_m_sia              )
   deallocate(ice%state%vx_m_ssa              )
   deallocate(ice%state%vy_m_ssa              )
   deallocate(ice%state%vx_m_prev             )
   deallocate(ice%state%vy_m_prev             )
   deallocate(ice%state%d_help_b              )
   deallocate(ice%state%ratio_sl              )
   deallocate(ice%state%ratio_sl_x            )
   deallocate(ice%state%ratio_sl_y            )
   deallocate(ice%state%vx_b_g                )
   deallocate(ice%state%vy_b_g                )
   deallocate(ice%state%vz_b                  )
   deallocate(ice%state%vz_m                  )
   deallocate(ice%state%vx_s_g                )
   deallocate(ice%state%vy_s_g                )
   deallocate(ice%state%vz_s                  )
   deallocate(ice%state%vz_sl                 )
   deallocate(ice%state%flui_ave_sia          )
   deallocate(ice%state%h_diff                )
   deallocate(ice%state%qx                    )
   deallocate(ice%state%qy                    )
   deallocate(ice%state%q_gl_g                )
   deallocate(ice%state%q_geo                 )
   deallocate(ice%state%temp_b                )
   deallocate(ice%state%temph_b               )
   deallocate(ice%state%Q_bm                  )
   deallocate(ice%state%Q_bm_float            )
   deallocate(ice%state%Q_tld                 )
   deallocate(ice%state%Q_b_tot               )
   deallocate(ice%state%Q_b_apl               )
   deallocate(ice%state%H_w                   )
   deallocate(ice%state%q_w                   )
   deallocate(ice%state%q_w_x                 )
   deallocate(ice%state%q_w_y                 )
   deallocate(ice%state%runoff                )
   deallocate(ice%state%runoff_apl            )
   deallocate(ice%state%as_perp               )
   deallocate(ice%state%as_perp_apl           )
   deallocate(ice%state%mb_source_apl         )
   deallocate(ice%state%mb_source             )
   deallocate(ice%state%smb_corr              )
   deallocate(ice%state%calving               )
   deallocate(ice%state%calving_apl           )
   deallocate(ice%state%mask_mar              )
   deallocate(ice%state%cst_dist              )
   deallocate(ice%state%i_cst                 )
   deallocate(ice%state%j_cst                 )
   deallocate(ice%state%cos_grad_tc           )
   deallocate(ice%state%dis_perp              ) 
   deallocate(ice%state%t_ocn                 )
   deallocate(ice%state%s_ocn                 )
   deallocate(ice%state%accum                 )
   deallocate(ice%state%accum_apl             )
   deallocate(ice%state%temp_s                )
   deallocate(ice%state%temp_g                )
   deallocate(ice%state%am_perp               )
   deallocate(ice%state%am_perp_st            )
   deallocate(ice%state%zs_neu                )
   deallocate(ice%state%zm_neu                )
   deallocate(ice%state%zb_neu                )
   deallocate(ice%state%zl_neu                )
   deallocate(ice%state%zs_aux                )
   deallocate(ice%state%zm_aux                )
   deallocate(ice%state%zb_aux                )
   deallocate(ice%state%H_c_neu               )
   deallocate(ice%state%H_t_neu               )
   deallocate(ice%state%vx_c      )
   deallocate(ice%state%vy_c      )
   deallocate(ice%state%vz_c      )
   deallocate(ice%state%temp_c    )
   deallocate(ice%state%temp_c_neu)
   deallocate(ice%state%temp_c_m  )
   deallocate(ice%state%temp_c_help)
   deallocate(ice%state%age_c     )
   deallocate(ice%state%age_c_neu )
   deallocate(ice%state%sigma_c   )
   deallocate(ice%state%txz_c   )
   deallocate(ice%state%tyz_c   )
   deallocate(ice%state%enh_c     )
   deallocate(ice%state%de_ssa    ) 
   deallocate(ice%state%vis_int_g )
   deallocate(ice%state%vis_ave_g )
   deallocate(ice%state%vis_ave_g_smooth )
   deallocate(ice%state%vis_int_sgxy )
   deallocate(ice%state%beta_drag )
   deallocate(ice%state%vx_g      )
   deallocate(ice%state%vy_g      )
   deallocate(ice%state%vx_t           )
   deallocate(ice%state%vy_t           )
   deallocate(ice%state%vz_t           )
   deallocate(ice%state%omega_t        )
   deallocate(ice%state%omega_t_neu    )
   deallocate(ice%state%temp_t_m       )
   deallocate(ice%state%age_t          )
   deallocate(ice%state%age_t_neu      )
   deallocate(ice%state%sigma_t        )
   deallocate(ice%state%txz_t          )
   deallocate(ice%state%tyz_t          )
   deallocate(ice%state%enh_t          )
   deallocate(ice%state%temp_r         )
   deallocate(ice%state%temp_r_neu     )
   deallocate(ice%state%enth_c         )
   deallocate(ice%state%enth_c_neu     )
   deallocate(ice%state%omega_c        )
   deallocate(ice%state%omega_c_neu    )
   deallocate(ice%state%enth_t         )
   deallocate(ice%state%enth_t_neu     )
   deallocate(ice%state%de_c           )
   deallocate(ice%state%lambda_shear_c )
   deallocate(ice%state%de_t           ) 
   deallocate(ice%state%lambda_shear_t ) 
   deallocate(ice%state%check_point    )
   deallocate(ice%state%weigh_ssta_sia_x)
   deallocate(ice%state%weigh_ssta_sia_y)
   deallocate(ice%state%p_b            )
   deallocate(ice%state%p_b_red        )
   deallocate(ice%state%p_b_red_lim    )
   deallocate(ice%state%tau_b          )

   deallocate(ice%state%gamma_slide_inv    )
   deallocate(ice%state%sub_melt_flag      )

   deallocate(ice%state%upH_x_1            )
   deallocate(ice%state%upH_x_2            )
   deallocate(ice%state%upH_y_1            )
   deallocate(ice%state%upH_y_2            )

   deallocate(ice%state%czs2            )
   deallocate(ice%state%czs3            )

   return

  end subroutine sico_dealloc 


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s i c o _ w r i t e _ r e s t a r t
  ! Purpose  :  Write sicopolis restart file
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sico_write_restart(fnm,ice)

    use dim_name, only: dim_x, dim_y, dim_kc, dim_kt, dim_kr

    implicit none

    type(sico_class) :: ice

    character (len=*) :: fnm
    integer :: ncid, i, j, kc, kt
    real(wp), dimension(:,:), allocatable :: H, H_cold, H_temp


    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,dim_x,x=ice%grid%x0,dx=ice%grid%dx,nx=ice%grid%IMAX+1,axis="x",units="m",ncid=ncid)
    call nc_write_dim(fnm,dim_y,x=ice%grid%y0,dx=ice%grid%dx,nx=ice%grid%JMAX+1,axis="y",units="m",ncid=ncid)
    call nc_write_dim(fnm,dim_kc,x=1,dx=1,nx=ice%grid%KCMAX+1,axis="z",units="#",ncid=ncid)
    call nc_write_dim(fnm,dim_kt,x=1,dx=1,nx=ice%grid%KTMAX+1,axis="z",units="#",ncid=ncid)
    call nc_write_dim(fnm,dim_kr,x=1,dx=1,nx=ice%grid%KRMAX+1,axis="z",units="#",ncid=ncid)

    call nc_write(fnm,"x", ice%grid%xi, &
      dims=[dim_x],start=[1],count=[ice%grid%IMAX+1], &
      long_name="x-coordinate",grid_mapping="polar_stereographic",units="m",ncid=ncid)    
    call nc_write(fnm,"y", ice%grid%eta, &
      dims=[dim_y],start=[1],count=[ice%grid%JMAX+1], &
      long_name="y-coordinate",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

    call nc_write(fnm,"temp_s", real(ice%state%temp_s,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="surface temperature",grid_mapping="polar_stereographic",units="degC",ncid=ncid)    

    call nc_write(fnm,"temp_g", real(ice%state%temp_g,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="ground temperature",grid_mapping="polar_stereographic",units="degC",ncid=ncid)    

    call nc_write(fnm,"as_perp", real(ice%state%as_perp,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="surface mass balance",grid_mapping="polar_stereographic",units="m/s",ncid=ncid)    

    call nc_write(fnm,"as_perp_apl", real(ice%state%as_perp_apl,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="applied surface mass balance",grid_mapping="polar_stereographic",units="m/s",ncid=ncid)    

    call nc_write(fnm,"calving_apl", real(ice%state%calving_apl,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="applied calving",grid_mapping="applied calving",units="m/s",ncid=ncid)    

    call nc_write(fnm,"Q_b_apl", real(ice%state%Q_b_apl,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="applied basal melting",grid_mapping="polar_stereographic",units="m/s",ncid=ncid)    

    call nc_write(fnm,"maske", int(ice%state%maske,4), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="ice-land-sea mask",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

    call nc_write(fnm,"maske_old", int(ice%state%maske_old,4), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="ice-land-sea mask (old)",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

    call nc_write(fnm,"n_cts", int(ice%state%n_cts,4), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="mask for polythermal conditions",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

    call nc_write(fnm,"kc_cts", int(ice%state%kc_cts,4), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Grid index of the CTS position",grid_mapping="polar_stereographic",units="/",ncid=ncid)    

    call nc_write(fnm,"zs", real(ice%state%zs,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Topography of the free surface",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

    call nc_write(fnm,"zs_std", real(ice%state%zs_std,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="standard deviation of the free surface elevation",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

    call nc_write(fnm,"zm", real(ice%state%zm,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Topography of z=zm interface",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

    call nc_write(fnm,"zb", real(ice%state%zb,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Topography of the ice base",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

    call nc_write(fnm,"zl", real(ice%state%zl,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Topography of the lithosphere surface",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

    call nc_write(fnm,"zl_fil", real(ice%state%zl_fil,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Topography of the smoothed lithosphere surface",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

    call nc_write(fnm,"zl0", real(ice%state%zl0,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Topography of the relaxed lithosphere surface",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

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

    call nc_write(fnm,"H", real(H,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Thickness of the temperate ice layer",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

    call nc_write(fnm,"H_cold", real(H_cold,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Thickness of the cold ice layer",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

    call nc_write(fnm,"H_temp", real(H_temp,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="thickness of the temperate ice layer",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

    deallocate(H_cold,H_temp)

    call nc_write(fnm,"Q_bm", real(ice%state%Q_bm,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Basal melting rate",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"Q_bm_float", real(ice%state%Q_bm_float,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Basal melting rate below ice shelf",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"Q_tld", real(ice%state%Q_tld,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Water drainage from the temperate layer",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"am_perp", real(ice%state%am_perp,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Volume flux across the z=zm interface",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"qx", real(ice%state%qx,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Horizontal volume flux qx",grid_mapping="polar_stereographic",units="m2/yr",ncid=ncid)    

    call nc_write(fnm,"qy", real(ice%state%qy,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Horizontal volume flux qy",grid_mapping="polar_stereographic",units="m2/yr",ncid=ncid)    

    call nc_write(fnm,"dzs_dt", real(ice%state%dzs_dtau,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Rate of change of the topography of the free surface",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"dzm_dt", real(ice%state%dzm_dtau,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Rate of change of the topography of the z=zm interface",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"dzb_dt", real(ice%state%dzb_dtau,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Rate of change of the topography of the ice base",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"dzl_dt", real(ice%state%dzb_dtau,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Rate of change of the topography of the ice base",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"dH_c_dt", real(ice%state%dH_c_dtau,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Rate of change of the thickness of the upper (kc) ice layer",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"dH_t_dt", real(ice%state%dH_t_dtau,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Rate of change of the thickness of the lower (kt) ice layer",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"vx_b_g", real(ice%state%vx_b_g,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Horizontal velocity vx at the ice base",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"vy_b_g", real(ice%state%vy_b_g,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Horizontal velocity vy at the ice base",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"vz_b", real(ice%state%vz_b,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Vertical velocity vz at the ice base",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"vx_s_g", real(ice%state%vx_s_g,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Horizontal velocity vx at the ice surface",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"vy_s_g", real(ice%state%vy_s_g,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Horizontal velocity vy at the ice surface",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"vz_s", real(ice%state%vz_s,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Vertical velocity vz at the ice surface",grid_mapping="polar_stereographic",units="m/yr",ncid=ncid)    

    call nc_write(fnm,"temp_b", real(ice%state%temp_b,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Temperature at the ice base",grid_mapping="polar_stereographic",units="degC",ncid=ncid)    

    call nc_write(fnm,"temph_b", real(ice%state%temph_b,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Temperature at the ice base relative to the pressure melting point",grid_mapping="polar_stereographic",units="degC",ncid=ncid)    

    call nc_write(fnm,"p_b_w", real(ice%state%p_b_w,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Basal water pressure",grid_mapping="polar_stereographic",units="Pa",ncid=ncid)    

    call nc_write(fnm,"H_w", real(ice%state%H_w,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Thickness of the water column under the ice base",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

    call nc_write(fnm,"vx_c", real(ice%state%vx_c,wp), &
      dims=[dim_kc,dim_y,dim_x],start=[1,1,1],count=[ice%grid%KCMAX+1,ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Horizontal velocity vx in the upper (kc) ice layer",grid_mapping="polar_stereographic",units="m/a",ncid=ncid)

    call nc_write(fnm,"vy_c", real(ice%state%vy_c,wp), &
      dims=[dim_kc,dim_y,dim_x],start=[1,1,1],count=[ice%grid%KCMAX+1,ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Horizontal velocity vy in the upper (kc) ice layer",grid_mapping="polar_stereographic",units="m/a",ncid=ncid)

    call nc_write(fnm,"vz_c", real(ice%state%vz_c,wp), &
      dims=[dim_kc,dim_y,dim_x],start=[1,1,1],count=[ice%grid%KCMAX+1,ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Horizontal velocity vy in the upper (kc) ice layer",grid_mapping="polar_stereographic",units="m/a",ncid=ncid)

    call nc_write(fnm,"vx_t", real(ice%state%vx_t,wp), &
      dims=[dim_kt,dim_y,dim_x],start=[1,1,1],count=[ice%grid%KTMAX+1,ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Horizontal velocity vx in the lower (kt) ice layer",grid_mapping="polar_stereographic",units="m/a",ncid=ncid)

    call nc_write(fnm,"vy_t", real(ice%state%vy_t,wp), &
      dims=[dim_kt,dim_y,dim_x],start=[1,1,1],count=[ice%grid%KTMAX+1,ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Horizontal velocity vy in the lower (kt) ice layer",grid_mapping="polar_stereographic",units="m/a",ncid=ncid)

    call nc_write(fnm,"vz_t", real(ice%state%vz_t,wp), &
      dims=[dim_kt,dim_y,dim_x],start=[1,1,1],count=[ice%grid%KTMAX+1,ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Horizontal velocity vy in the lower (kt) ice layer",grid_mapping="polar_stereographic",units="m/a",ncid=ncid)

    call nc_write(fnm,"vx_m_ssa", real(ice%state%vx_m_ssa,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="SSA x velocity",grid_mapping="polar_stereographic",units="m/a",ncid=ncid)    

    call nc_write(fnm,"vy_m_ssa", real(ice%state%vy_m_ssa,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="SSA y velocity",grid_mapping="polar_stereographic",units="m/a",ncid=ncid)    

    call nc_write(fnm,"vis_ave_g", real(ice%state%vis_ave_g,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Depth-averaged viscosity of the SSA",grid_mapping="polar_stereographic",units="",ncid=ncid)    

    call nc_write(fnm,"temp_c", real(ice%state%temp_c,wp), &
      dims=[dim_kc,dim_y,dim_x],start=[1,1,1],count=[ice%grid%KCMAX+1,ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Temperature in the upper (kc) ice layer",grid_mapping="polar_stereographic",units="degC",ncid=ncid)

    call nc_write(fnm,"omega_t", real(ice%state%omega_t,wp), &
      dims=[dim_kt,dim_y,dim_x],start=[1,1,1],count=[ice%grid%KTMAX+1,ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Water content in the lower (kt) ice layer",grid_mapping="polar_stereographic",units="1",ncid=ncid)

    call nc_write(fnm,"omega_c", real(ice%state%omega_c,wp), &
      dims=[dim_kc,dim_y,dim_x],start=[1,1,1],count=[ice%grid%KCMAX+1,ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Water content in the upper (kc) ice layer",grid_mapping="polar_stereographic",units="1",ncid=ncid)

    call nc_write(fnm,"temp_r", real(ice%state%temp_r,wp), &
      dims=[dim_kr,dim_y,dim_x],start=[1,1,1],count=[ice%grid%KRMAX+1,ice%grid%JMAX+1,ice%grid%IMAX+1], &
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

    call nc_write(fnm,"enth_c", real(ice%state%enth_c,wp), &
      dims=[dim_kc,dim_y,dim_x],start=[1,1,1],count=[ice%grid%KCMAX+1,ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Enthalpy in the upper (kc) ice layer",grid_mapping="polar_stereographic",units="J7kg",ncid=ncid)

    call nc_write(fnm,"enth_t", real(ice%state%enth_t,wp), &
      dims=[dim_kt,dim_y,dim_x],start=[1,1,1],count=[ice%grid%KTMAX+1,ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Enthalpy in the lower (kt) ice layer",grid_mapping="polar_stereographic",units="J7kg",ncid=ncid)

    call nc_write(fnm,"enh_c", real(ice%state%enh_c,wp), &
      dims=[dim_kc,dim_y,dim_x],start=[1,1,1],count=[ice%grid%KCMAX+1,ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Flow enhancement factor in the upper (kc) ice layer",grid_mapping="polar_stereographic",units="1",ncid=ncid)

    call nc_write(fnm,"enh_t", real(ice%state%enh_t,wp), &
      dims=[dim_kt,dim_y,dim_x],start=[1,1,1],count=[ice%grid%KTMAX+1,ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Flow enhancement factor in the lower (kt) ice layer",grid_mapping="polar_stereographic",units="1",ncid=ncid)

    call nc_write(fnm,"age_c", real(ice%state%age_c,wp), &
      dims=[dim_kc,dim_y,dim_x],start=[1,1,1],count=[ice%grid%KCMAX+1,ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Age in the upper (kc) ice layer",grid_mapping="polar_stereographic",units="a",ncid=ncid)

    call nc_write(fnm,"age_t", real(ice%state%age_t,wp), &
      dims=[dim_kt,dim_y,dim_x],start=[1,1,1],count=[ice%grid%KTMAX+1,ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="Age in the lower (kt) ice layer",grid_mapping="polar_stereographic",units="a",ncid=ncid)

    call nc_write(fnm,"q_geo", real(ice%state%q_geo,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="geothermal heat flux",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid)    

    call nc_write(fnm,"h_sed", real(ice%state%h_sed,wp), &
      dims=[dim_y,dim_x],start=[1,1],count=[ice%grid%JMAX+1,ice%grid%IMAX+1], &
      long_name="sediment thickness",grid_mapping="polar_stereographic",units="m",ncid=ncid)    

    call nc_close(ncid)


  end subroutine sico_write_restart


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  i c e _ r e a d _ r e s t a r t
  ! Purpose  :  read ice restart file
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sico_read_restart(fnm,ice)

    implicit none

    type(sico_class), intent(inout) :: ice

    character (len=*) :: fnm
    real(wp), dimension(:,:), allocatable :: H, H_cold, H_temp


    call nc_read(fnm,"x", ice%grid%xi)

    call nc_read(fnm,"y", ice%grid%eta)

    call nc_read(fnm,"temp_s", ice%state%temp_s)

    call nc_read(fnm,"temp_g", ice%state%temp_g)

    call nc_read(fnm,"as_perp", ice%state%as_perp)

    call nc_read(fnm,"as_perp_apl", ice%state%as_perp_apl)

    call nc_read(fnm,"calving_apl", ice%state%calving_apl)

    call nc_read(fnm,"Q_b_apl", ice%state%Q_b_apl)

    call nc_read(fnm,"maske", ice%state%maske)

    call nc_read(fnm,"maske_old", ice%state%maske_old)

    call nc_read(fnm,"n_cts", ice%state%n_cts)

    call nc_read(fnm,"kc_cts", ice%state%kc_cts)

    call nc_read(fnm,"zs", ice%state%zs)

    call nc_read(fnm,"zs_std", ice%state%zs_std)

    call nc_read(fnm,"zm", ice%state%zm)

    call nc_read(fnm,"zb", ice%state%zb)

    call nc_read(fnm,"zl", ice%state%zl)

    call nc_read(fnm,"zl_fil", ice%state%zl_fil)

    call nc_read(fnm,"zl0", ice%state%zl0)

    allocate(H(0:ice%grid%JMAX,0:ice%grid%IMAX))
    allocate(H_cold(0:ice%grid%JMAX,0:ice%grid%IMAX))
    allocate(H_temp(0:ice%grid%JMAX,0:ice%grid%IMAX))

    call nc_read(fnm,"H", H)

    call nc_read(fnm,"H_cold", H_cold)

    call nc_read(fnm,"H_temp", H_temp)

    if (ice%par%calcmod==1) then
      ice%state%H_c     = H_cold
      ice%state%H_t     = H_temp
    else if (ice%par%calcmod==0 .or. ice%par%calcmod==2 .or. ice%par%calcmod==3 .or. ice%par%calcmod==-1) then
      ice%state%H_c     = H
      ice%state%H_t     = 0.0_wp
    endif

    deallocate(H_cold,H_temp)

    call nc_read(fnm,"Q_bm", ice%state%Q_bm)

    call nc_read(fnm,"Q_bm_float", ice%state%Q_bm_float)

    call nc_read(fnm,"Q_tld", ice%state%Q_tld)

    call nc_read(fnm,"am_perp", ice%state%am_perp)

    call nc_read(fnm,"qx", ice%state%qx)

    call nc_read(fnm,"qy", ice%state%qy)

    call nc_read(fnm,"dzs_dt", ice%state%dzs_dtau)

    call nc_read(fnm,"dzm_dt", ice%state%dzm_dtau)

    call nc_read(fnm,"dzb_dt", ice%state%dzb_dtau)

    call nc_read(fnm,"dzl_dt", ice%state%dzb_dtau)

    call nc_read(fnm,"dH_c_dt", ice%state%dH_c_dtau)

    call nc_read(fnm,"dH_t_dt", ice%state%dH_t_dtau)

    call nc_read(fnm,"vx_b_g", ice%state%vx_b_g)

    call nc_read(fnm,"vy_b_g", ice%state%vy_b_g)

    call nc_read(fnm,"vz_b", ice%state%vz_b)

    call nc_read(fnm,"vx_s_g", ice%state%vx_s_g)

    call nc_read(fnm,"vy_s_g", ice%state%vy_s_g)

    call nc_read(fnm,"vz_s", ice%state%vz_s)

    call nc_read(fnm,"temp_b", ice%state%temp_b)

    call nc_read(fnm,"temph_b", ice%state%temph_b)

    call nc_read(fnm,"p_b_w", ice%state%p_b_w)

    call nc_read(fnm,"H_w", ice%state%H_w)

    call nc_read(fnm,"vx_c", ice%state%vx_c)

    call nc_read(fnm,"vy_c", ice%state%vy_c)

    call nc_read(fnm,"vz_c", ice%state%vz_c)

    call nc_read(fnm,"vx_t", ice%state%vx_t)

    call nc_read(fnm,"vy_t", ice%state%vy_t)

    call nc_read(fnm,"vz_t", ice%state%vz_t)

    call nc_read(fnm,"vx_m_ssa", ice%state%vx_m_ssa)

    call nc_read(fnm,"vy_m_ssa", ice%state%vy_m_ssa)

    call nc_read(fnm,"vis_ave_g", ice%state%vis_ave_g)

    call nc_read(fnm,"temp_c", ice%state%temp_c)

    call nc_read(fnm,"omega_t", ice%state%omega_t)

    call nc_read(fnm,"omega_c", ice%state%omega_c)

    call nc_read(fnm,"temp_r", ice%state%temp_r)

    call nc_read(fnm,"enth_c", ice%state%enth_c)

    call nc_read(fnm,"enth_t", ice%state%enth_t)

    call nc_read(fnm,"enh_c", ice%state%enh_c)

    call nc_read(fnm,"enh_t", ice%state%enh_t)

    call nc_read(fnm,"age_c", ice%state%age_c)

    call nc_read(fnm,"age_t", ice%state%age_t)

    call nc_read(fnm,"q_geo", ice%state%q_geo)

    call nc_read(fnm,"h_sed", ice%state%h_sed)

  end subroutine sico_read_restart

end module sicopolis
