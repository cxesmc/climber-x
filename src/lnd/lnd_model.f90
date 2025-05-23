!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : l n d _ m o d e l
!
!  Purpose : main land model
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit and Andrey Ganopolski
!
! This file is part of CLIMBER-X.
!
! CLIMBER-X is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! CLIMBER-X is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with CLIMBER-X.  If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module lnd_model 

    use precision, only : wp
    use constants, only : T0, sigma, rho_w
    use timer, only : doy, year, year_now, sec_day, sec_year, dt_lnd, nday_year, nmon_year
    use timer, only : time_soy_lnd, time_eoy_lnd, time_eom_lnd
    use climber_grid, only : ni, nj, lon, lat, area
    use control, only : out_dir, lnd_restart, restart_in_dir, check_water, in_dir
    use lnd_grid, only : lnd_grid_init
    use lnd_grid, only : z, z_l, z_c, z_int, dz
    use lnd_grid, only : nx, ny, npft, nsurf, nsoil, ncarb, nl, nl_l, nlc, rdz_pos, rdz_neg
    use lnd_grid, only : ic_min, ic_peat, ic_shelf, ic_ice, ic_lake
    use lnd_grid, only : i_bare, i_ice, i_lake, i_surf, i_pft, i_soil, is_lake
    use lnd_params, only : lnd_params_init, i_init_veg
    use lnd_params, only : dt, rdt, dt_day_c, dt_day_v, dt_day_carb, dt_day_veg
    use lnd_params, only : time_call_veg, time_call_carb, time_call_carb_p
    use lnd_params, only : l_dynveg, pft_fix_file, l_fixlai, lai_fix_file
    use lnd_params, only : l_co2_fert_lim, co2_fert_lim_min, co2_fert_lim_max, i_weathering, l_river_export
    use lnd_params, only : veg_par, pft_par, snow_par, surf_par, hydro_par, soil_par, soilc_par, peat_par
    use lnd_params, only : mineral, topmodel, dyptop, nmonwet
    use lnd_params, only : weath_gemco2_par, weath_uhh_par
    use lnd_def, only : lnd_class, lnd_0d_class, lnd_2d_class
    use ncio


    implicit none 

    private
    public :: lnd_init, lnd_update_wrapper, lnd_end
    public :: lnd_read_restart, lnd_write_restart

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l n d _ u p d a t e _ w r a p p e r
  !   Purpose    :  wrapper for land update
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lnd_update_wrapper(lnd)

    !$  use omp_lib
    use carbon_inventory_mod
    use carbon_flux_atm_lnd_mod

    implicit none

    type(lnd_class) :: lnd

    integer :: i, j, ii, jj, n, ip1, im1, jp1, jm1
    real(wp) :: pft_max
    real(wp) :: global_landc, global_landc13, global_landc14
    real(wp) :: landc, landc13, landc14
    real(wp) :: global_burc, global_burc13, global_burc14
    real(wp) :: burc, burc13, burc14
    real(wp) :: lnd_cell_area, weath_carb_glob, weath_sil_glob


    time_call_veg    = (mod(doy,dt_day_v) .eq. 0)
    time_call_carb   = (mod(doy,dt_day_c) .eq. 0)
    time_call_carb_p = (mod(doy,dt_day_carb) .eq. 0)

    if (time_soy_lnd .and. veg_par%iseed.eq.2) then
      ! update seed distribution
      ! needs to be done here because it requires knowledge of 2D fields
      do j=1,ny
        do i=1,nx
          if( lnd%l2d(i,j)%mask_lnd.eq.1 ) then
            jm1 = max(1,j-1)
            jp1 = min(ny,j+1)
            im1 = i-1
            if (im1.eq.0) im1 = nx
            ip1 = i+1
            if (ip1.eq.nx+1) ip1 = 1
            do n=1,npft
              pft_max = 0._wp
              do jj=jm1,jp1
                do ii=im1,ip1
                  if( lnd%l2d(ii,jj)%mask_lnd.eq.1 ) then
                    pft_max = max(pft_max,lnd%l2d(ii,jj)%pft_frac(n))
                  endif
                enddo
              enddo
              if (pft_max.gt.veg_par%seed_pft_min) then
                lnd%l2d(i,j)%seed_frac(n) = veg_par%seed_fraction
              else
                lnd%l2d(i,j)%seed_frac(n) = 0._wp
              endif
            enddo
          endif
        enddo
      enddo
    endif

    !$omp parallel do &
    !$omp private ( i, j, n)
    do n=1,lnd%ncells
      i = lnd%ij_1d(1,n)
      j = lnd%ij_1d(2,n)
      !!$ print *,n,i,j,omp_get_thread_num(),'/',omp_get_num_threads()
      ! update land model
      call lnd_update(lnd%l2d(i,j),lnd%l0d,i,j)
    enddo
    !$omp end parallel do

    ! compute atmosphere-land carbon flux
    if (time_soy_lnd) then
      lnd%l0d%Cflx_atm_lnd   = 0._wp  
      lnd%l0d%C13flx_atm_lnd = 0._wp
      lnd%l0d%C14flx_atm_lnd = 0._wp
      lnd%l0d%ch4_emis = 0._wp  
    endif
    if (time_eom_lnd) then
       do j=1,ny
         do i=1,nx
           if( lnd%l2d(i,j)%mask_lnd.eq.1 ) then
             call carbon_flux_atm_lnd(lnd%l2d(i,j)%f_veg,lnd%l2d(i,j)%f_peat,lnd%l2d(i,j)%f_ice_grd,lnd%l2d(i,j)%f_shelf,lnd%l2d(i,j)%f_lake,area(i,j), &
               lnd%l2d(i,j)%npp_real,lnd%l2d(i,j)%npp13_real,lnd%l2d(i,j)%npp14_real,lnd%l2d(i,j)%soil_resp,lnd%l2d(i,j)%soil_resp13,lnd%l2d(i,j)%soil_resp14, &
               lnd%l2d(i,j)%Cflx_atm_lnd,lnd%l2d(i,j)%C13flx_atm_lnd,lnd%l2d(i,j)%C14flx_atm_lnd, &
               lnd%l0d%Cflx_atm_lnd,lnd%l0d%C13flx_atm_lnd,lnd%l0d%C14flx_atm_lnd) 
           endif
         enddo
       enddo
     endif

     ! cumulate annual global methane emissions
     do j=1,ny
       do i=1,nx
         lnd%l0d%ch4_emis = lnd%l0d%ch4_emis  &
           + (lnd%l2d(i,j)%ch4_emis_wetland*max(0._wp,lnd%l2d(i,j)%f_wetland*lnd%l2d(i,j)%f_veg-lnd%l2d(i,j)%f_peat) &  ! wetland
           + lnd%l2d(i,j)%ch4_emis_shelf*lnd%l2d(i,j)%f_shelf & ! shelf
           + lnd%l2d(i,j)%ch4_emis_lake*lnd%l2d(i,j)%f_lake & ! lakes
           + lnd%l2d(i,j)%ch4_emis_peat*lnd%l2d(i,j)%f_peat) & ! peatland
           * dt * area(i,j)   ! kgCH4
       enddo
     enddo

     ! carbon inventory 
     if (time_eoy_lnd) then
       global_landc = 0._wp
       global_landc13 = 0._wp
       global_landc14 = 0._wp
       global_burc = 0._wp
       global_burc13 = 0._wp
       global_burc14 = 0._wp
       do j=1,ny
         do i=1,nx
           if( lnd%l2d(i,j)%mask_lnd.eq.1 ) then
             call carbon_inventory(lnd%l2d(i,j)%frac_surf,lnd%l2d(i,j)%f_veg,lnd%l2d(i,j)%f_peat, &
               lnd%l2d(i,j)%f_ice_grd,lnd%l2d(i,j)%f_shelf,lnd%l2d(i,j)%f_lake, &
               lnd%l2d(i,j)%veg_c,lnd%l2d(i,j)%veg_c13,lnd%l2d(i,j)%veg_c14, &
               lnd%l2d(i,j)%litter_c,lnd%l2d(i,j)%fast_c,lnd%l2d(i,j)%slow_c, &
               lnd%l2d(i,j)%litter_c13,lnd%l2d(i,j)%fast_c13,lnd%l2d(i,j)%slow_c13, &
               lnd%l2d(i,j)%litter_c14,lnd%l2d(i,j)%fast_c14,lnd%l2d(i,j)%slow_c14, &
               lnd%l2d(i,j)%litter_c_peat,lnd%l2d(i,j)%acro_c,lnd%l2d(i,j)%cato_c, &
               lnd%l2d(i,j)%litter_c13_peat,lnd%l2d(i,j)%acro_c13,lnd%l2d(i,j)%cato_c13, &
               lnd%l2d(i,j)%litter_c14_peat,lnd%l2d(i,j)%acro_c14,lnd%l2d(i,j)%cato_c14, &
               lnd%l2d(i,j)%litter_c_ice,lnd%l2d(i,j)%fast_c_ice,lnd%l2d(i,j)%slow_c_ice,&
               lnd%l2d(i,j)%litter_c13_ice,lnd%l2d(i,j)%fast_c13_ice,lnd%l2d(i,j)%slow_c13_ice, &
               lnd%l2d(i,j)%litter_c14_ice,lnd%l2d(i,j)%fast_c14_ice,lnd%l2d(i,j)%slow_c14_ice, &
               lnd%l2d(i,j)%litter_c_shelf,lnd%l2d(i,j)%fast_c_shelf,lnd%l2d(i,j)%slow_c_shelf, &
               lnd%l2d(i,j)%litter_c13_shelf,lnd%l2d(i,j)%fast_c13_shelf,lnd%l2d(i,j)%slow_c13_shelf, &
               lnd%l2d(i,j)%litter_c14_shelf,lnd%l2d(i,j)%fast_c14_shelf,lnd%l2d(i,j)%slow_c14_shelf, &
               lnd%l2d(i,j)%litter_c_lake,lnd%l2d(i,j)%fast_c_lake,lnd%l2d(i,j)%slow_c_lake, &
               lnd%l2d(i,j)%litter_c13_lake,lnd%l2d(i,j)%fast_c13_lake,lnd%l2d(i,j)%slow_c13_lake, &
               lnd%l2d(i,j)%litter_c14_lake,lnd%l2d(i,j)%fast_c14_lake,lnd%l2d(i,j)%slow_c14_lake, &
               landc,landc13,landc14,burc,burc13,burc14)
             global_landc   = global_landc   + landc * area(i,j) ! KgC
             global_landc13 = global_landc13 + landc13 * area(i,j) ! KgC
             global_landc14 = global_landc14 + landc14 * area(i,j) ! KgC
             global_burc   = global_burc   + burc * area(i,j) ! KgC
             global_burc13 = global_burc13 + burc13 * area(i,j) ! KgC
             global_burc14 = global_burc14 + burc14 * area(i,j) ! KgC
           endif
         enddo
       enddo

       if (soilc_par%l_burial) then
         ! carbon is effectively buried and lost from the system (has to return through geological fluxes to avoid drifting)
         ! carbon flux to burial
         lnd%l0d%Cflx_burial    = (global_burc   - lnd%l0d%burc)   ! kgC/yr
         lnd%l0d%C13flx_burial  = (global_burc13 - lnd%l0d%burc13) ! kgC/yr
         lnd%l0d%C14flx_burial  = (global_burc14 - lnd%l0d%burc14) ! kgC/yr
       else
         ! no burial
         lnd%l0d%Cflx_burial    = 0._wp  
         lnd%l0d%C13flx_burial  = 0._wp
         lnd%l0d%C14flx_burial  = 0._wp
         ! carbon reaching the burial layer goes directly into the atmosphere
         lnd%l0d%Cflx_atm_lnd    = lnd%l0d%Cflx_atm_lnd   - (global_burc   - lnd%l0d%burc)   ! kgC/yr
         lnd%l0d%C13flx_atm_lnd  = lnd%l0d%C13flx_atm_lnd - (global_burc13 - lnd%l0d%burc13) ! kgC/yr
         lnd%l0d%C14flx_atm_lnd  = lnd%l0d%C14flx_atm_lnd - (global_burc14 - lnd%l0d%burc14) ! kgC/yr
         lnd%l2d(:,:)%Cflx_atm_lnd    = lnd%l2d(:,:)%Cflx_atm_lnd   - (global_burc   - lnd%l0d%burc)   / dt_lnd / sum(area(:,:)) ! kgC/s
         lnd%l2d(:,:)%C13flx_atm_lnd  = lnd%l2d(:,:)%C13flx_atm_lnd - (global_burc13 - lnd%l0d%burc13) / dt_lnd / sum(area(:,:)) ! kgC/s
         lnd%l2d(:,:)%C14flx_atm_lnd  = lnd%l2d(:,:)%C14flx_atm_lnd - (global_burc14 - lnd%l0d%burc14) / dt_lnd / sum(area(:,:)) ! kgC/s
       endif

       ! save total land carbon for next year
       lnd%l0d%landc   = global_landc  
       lnd%l0d%landc13 = global_landc13  
       lnd%l0d%landc14 = global_landc14  
       ! save buried land carbon for next year
       lnd%l0d%burc   = global_burc    
       lnd%l0d%burc13 = global_burc13  
       lnd%l0d%burc14 = global_burc14  
     endif

     ! average land-atmosphere carbon flux
     if (time_eoy_lnd .and. year>1) then
       if (year==2) lnd%l0d%Cflx_avg = lnd%l0d%Cflx_atm_lnd ! kgC/yr
       lnd%l0d%Cflx_avg = 0.99_wp*lnd%l0d%Cflx_avg + 0.01_wp*lnd%l0d%Cflx_atm_lnd
     endif

     ! average global weathering fluxes
     if (time_eoy_lnd .and. year>1) then
       ! compute global weathering fluxes
       weath_carb_glob = 0._wp
       weath_sil_glob  = 0._wp
       do j=1,nj
         do i=1,ni
           if (lnd%l2d(i,j)%f_veg.gt.0._wp) then
             lnd_cell_area = area(i,j)*lnd%l2d(i,j)%f_veg 
             ! carbonate and silicate weathering
             weath_carb_glob = weath_carb_glob + lnd%l2d(i,j)%weath_carb * lnd_cell_area * 12._wp *1e-3_wp ! mol C/m2/yr * m2 *g/mol * kg/g = kgC yr
             weath_sil_glob  = weath_sil_glob  + lnd%l2d(i,j)%weath_sil  * lnd_cell_area * 12._wp *1e-3_wp ! mol C/m2/yr * m2 *g/mol * kg/g = kgC yr
           endif
         enddo
       enddo
       ! average global weathering fluxes
       if (year==2) lnd%l0d%weath_carb_avg = weath_carb_glob ! kgC/yr
       lnd%l0d%weath_carb_avg = 0.99_wp*lnd%l0d%weath_carb_avg + 0.01_wp*weath_carb_glob
       if (year==2) lnd%l0d%weath_sil_avg = weath_sil_glob ! kgC/yr
       lnd%l0d%weath_sil_avg = 0.99_wp*lnd%l0d%weath_sil_avg + 0.01_wp*weath_sil_glob
     endif

  return

  end subroutine lnd_update_wrapper


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  lnd_update
  !   Purpose    :  update land column
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lnd_update(lnd,lndp,i,j)

    use surface_par_lnd
    use soil_par_mod
    use ice_par_mod
    use lake_par_mod
    use shelf_par_mod
    use veg_par_mod
    use soil_carbon_par_mod
    use shelf_carbon_par_mod
    use ice_carbon_par_mod
    use lake_carbon_par_mod
    use photosynthesis_mod
    use ebal_veg_mod
    use ebal_ice_mod
    use ebal_lake_mod
    use soil_temp_mod
    use ice_temp_mod
    use lake_temp_mod
    use sublake_temp_mod
    use shelf_temp_mod
    use hydrology_mod
    use soil_hydro_mod
    use water_check_mod
    use dyn_veg_mod
    use soil_carbon_mod
    use peat_carbon_mod
    use shelf_carbon_mod
    use ice_carbon_mod
    use lake_carbon_mod
    use carbon_trans_mod
    use dust_emis_mod
    use weathering_mod
    use carbon_export_mod
    use init_cell
    use end_cell

    implicit none 

    type(lnd_2d_class) :: lnd
    type(lnd_0d_class) :: lndp

    integer i, j, n
    real(wp) :: w_lgm


      if( time_soy_lnd ) then

       ! initialize soil if newly formed vegetated grid cell fraction
       if (lnd%f_veg.gt.0._wp .and. lnd%f_veg_old.eq.0._wp) then
        call init_cell_veg(lndp%c13_c12_atm,lndp%c14_c_atm, &
                          lnd%f_veg,lnd%f_shelf_old,lnd%frac_surf,lnd%pft_frac,lnd%t_shelf,lnd%theta_i_shelf,lnd%theta_w_shelf,lnd%theta_sat, &
                          lnd%t_skin,lnd%t_skin_veg,lnd%t_soil,lnd%w_can,lnd%s_can, &
                          lnd%w_snow(is_veg),lnd%w_snow_max(is_veg),lnd%h_snow(is_veg),lnd%mask_snow(is_veg), &
                          lnd%theta,lnd%theta_w,lnd%theta_i,lnd%w_w,lnd%w_i, &
                          lnd%alt,lnd%gdd5,lnd%gdd,lnd%phen,lnd%phen_acc,lnd%lai_bal,lnd%lai,lnd%sai,lnd%root_frac, &
                          lnd%litter_in_frac,lnd%gamma_dist, &
                          lnd%npp_ann,lnd%npp13_ann,lnd%npp14_ann, &
                          lnd%leaf_c,lnd%root_c, &
                          lnd%stem_c,lnd%veg_h,lnd%veg_c,lnd%veg_c13,lnd%veg_c14, &
                          lnd%litter_c,lnd%litter_c13,lnd%litter_c14,lnd%fast_c,lnd%fast_c13,lnd%fast_c14, &
                          lnd%slow_c,lnd%slow_c13,lnd%slow_c14, &
                          lnd%cato_c,lnd%cato_c13,lnd%cato_c14,lnd%litter_c_peat,lnd%litter_c13_peat,lnd%litter_c14_peat, &
                          lnd%acro_c,lnd%acro_c13,lnd%acro_c14,lnd%f_peat,lnd%f_peat_pot,lnd%w_table_min,lnd%w_table_peat,lnd%dCpeat_dt, &
                          lnd%z0m)
       endif

       ! initialize ice if newly formed ice grid cell fraction
       if (lnd%f_ice.gt.0._wp .and. lnd%f_ice_old.eq.0._wp) then
        call init_cell_ice(lnd%frac_surf,lnd%t_skin,lnd%t_ice, &
                          lnd%w_snow(is_ice),lnd%w_snow_max(is_ice),lnd%h_snow(is_ice),lnd%mask_snow(is_ice))
       endif

       ! initialize soil carbon below ice if newly formed grounded ice 
       if (lnd%f_ice_grd.gt.0._wp .and. lnd%f_ice_grd_old.eq.0._wp) then
        call init_cell_ice_grd(lnd%litter_c_ice,lnd%litter_c13_ice,lnd%litter_c14_ice,lnd%fast_c_ice, &
                              lnd%fast_c13_ice,lnd%fast_c14_ice,lnd%slow_c_ice,lnd%slow_c13_ice,lnd%slow_c14_ice)
       endif

       ! initialize shelf if newly formed shelf grid cell fraction
       if (lnd%f_shelf.gt.0._wp .and. lnd%f_shelf_old.eq.0._wp) then
        call init_cell_shelf(lnd%f_veg_old,lnd%t_soil,lnd%theta_i,lnd%theta_sat, &
                          lnd%t_shelf,lnd%w_w_shelf,lnd%w_i_shelf,lnd%theta_w_shelf,lnd%theta_i_shelf, &
                          lnd%litter_c_shelf,lnd%litter_c13_shelf,lnd%litter_c14_shelf, &
                          lnd%fast_c_shelf,lnd%fast_c13_shelf,lnd%fast_c14_shelf, &
                          lnd%slow_c_shelf,lnd%slow_c13_shelf,lnd%slow_c14_shelf)
       endif
     
       ! initialize lake if newly formed lake grid cell fraction
       if (lnd%f_lake.gt.0._wp .and. lnd%f_lake_old.eq.0._wp) then
        call init_cell_lake(lnd%frac_surf,lnd%t_skin, &
                          lnd%w_snow(is_lake),lnd%w_snow_max(is_lake),lnd%h_snow(is_lake),lnd%mask_snow(is_lake), &
                          lnd%f_veg_old,lnd%t_soil,lnd%theta_i,lnd%theta_sat, &
                          lnd%t_sublake,lnd%w_w_sublake,lnd%w_i_sublake,lnd%theta_w_sublake,lnd%theta_i_sublake, &
                          lnd%t_lake,lnd%f_i_lake,lnd%f_lake_ice, &
                          lnd%litter_c_lake,lnd%litter_c13_lake,lnd%litter_c14_lake, &
                          lnd%fast_c_lake,lnd%fast_c13_lake,lnd%fast_c14_lake, &
                          lnd%slow_c_lake,lnd%slow_c13_lake,lnd%slow_c14_lake)
       endif

       ! transfer soil carbon between different pools to account for changes in ice/shelf/lake area
       call carbon_trans(lnd%f_veg,lnd%f_veg_old,lnd%f_ice_grd,lnd%f_ice_grd_old,lnd%f_lake,lnd%f_lake_old,lnd%f_shelf,lnd%f_shelf_old, &
                         lnd%dCpeat_dt,lnd%f_peat_pot,lnd%f_peat,&
                         lnd%litter_c,lnd%fast_c,lnd%slow_c,lnd%litter_c13,lnd%fast_c13,lnd%slow_c13, &
                         lnd%litter_c14,lnd%fast_c14,lnd%slow_c14, &
                         lnd%litter_c_shelf,lnd%fast_c_shelf,lnd%slow_c_shelf,lnd%litter_c13_shelf, &
                         lnd%fast_c13_shelf,lnd%slow_c13_shelf, &
                         lnd%litter_c14_shelf,lnd%fast_c14_shelf,lnd%slow_c14_shelf, &
                         lnd%litter_c_lake,lnd%fast_c_lake,lnd%slow_c_lake,lnd%litter_c13_lake, &
                         lnd%fast_c13_lake,lnd%slow_c13_lake, &
                         lnd%litter_c14_lake,lnd%fast_c14_lake,lnd%slow_c14_lake, &
                         lnd%litter_c_ice,lnd%fast_c_ice,lnd%slow_c_ice,lnd%litter_c13_ice,lnd%fast_c13_ice,lnd%slow_c13_ice, &
                         lnd%litter_c14_ice,lnd%fast_c14_ice,lnd%slow_c14_ice, &
                         lnd%litter_c_peat,lnd%acro_c,lnd%cato_c,lnd%litter_c13_peat,lnd%acro_c13, &
                         lnd%cato_c13,lnd%litter_c14_peat,lnd%acro_c14,lnd%cato_c14)

       ! reset values if vegetated grid cell fraction disappeard
       if( lnd%f_veg .eq. 0._wp .and. lnd%f_veg_old .gt. 0._wp ) then
        call end_cell_veg(lndp%c13_c12_atm,lndp%c14_c_atm, &
                          lnd%pft_frac,lnd%alt, &
                          lnd%gdd5,lnd%gdd,lnd%phen,lnd%phen_acc,lnd%lai_bal,lnd%lai,lnd%sai,&
                          lnd%npp_ann,lnd%npp13_ann,lnd%npp14_ann, &
                          lnd%leaf_c,lnd%root_c, &
                          lnd%stem_c,lnd%veg_h,lnd%veg_c,lnd%veg_c13,lnd%veg_c14, &
                          lnd%w_can,lnd%s_can,lnd%t_skin,lnd%t_skin_veg,lnd%t_soil, &
                          lnd%w_snow(is_veg),lnd%h_snow(is_veg),lnd%mask_snow(is_veg), &
                          lnd%theta,lnd%theta_w,lnd%theta_i,lnd%theta_sat,lnd%w_w,lnd%w_i, &
                          lnd%litter_c,lnd%litter_c13,lnd%litter_c14,lnd%fast_c,lnd%fast_c13,lnd%fast_c14, &
                          lnd%slow_c,lnd%slow_c13,lnd%slow_c14, &
                          lnd%cato_c,lnd%cato_c13,lnd%cato_c14,lnd%litter_c_peat,lnd%litter_c13_peat,lnd%litter_c14_peat, &
                          lnd%acro_c,lnd%acro_c13,lnd%acro_c14,lnd%f_peat,lnd%f_peat_pot,lnd%w_table_min,lnd%w_table_peat, &
                          lnd%soil_resp(ic_min:ic_peat),lnd%soil_resp13(ic_min:ic_peat),lnd%soil_resp14(ic_min:ic_peat), &
                          lnd%soil_c_tot(ic_min:ic_peat),lnd%soil_c13_tot(ic_min:ic_peat),lnd%soil_c14_tot(ic_min:ic_peat), &
                          lnd%ch4_emis_wetland,lnd%ch4_emis_peat,lnd%c13h4_emis_wetland,lnd%c13h4_emis_peat, &
                          lnd%z0m)
       endif

       ! reset ice if ice grid cell fraction disappeard
       if( lnd%f_ice .eq. 0._wp .and. lnd%f_ice_old .gt. 0._wp ) then
        call end_cell_ice(lnd%t_skin,lnd%t_ice, &
                          lnd%w_snow(is_ice),lnd%h_snow(is_ice),lnd%mask_snow(is_ice))
       endif

       ! reset ice if ice grid cell fraction disappeard
       if( lnd%f_ice_grd .eq. 0._wp .and. lnd%f_ice_grd_old .gt. 0._wp ) then
        call end_cell_ice_grd(lnd%litter_c_ice,lnd%litter_c13_ice,lnd%litter_c14_ice,lnd%fast_c_ice, &
                          lnd%fast_c13_ice,lnd%fast_c14_ice,lnd%slow_c_ice,lnd%slow_c13_ice,lnd%slow_c14_ice, &
                          lnd%soil_resp(ic_ice),lnd%soil_resp13(ic_ice),lnd%soil_resp14(ic_ice), &
                          lnd%soil_c_tot(ic_ice),lnd%soil_c13_tot(ic_ice),lnd%soil_c14_tot(ic_ice))
       endif

       ! reset shelf if shelf grid cell fraction disappeard
       if( lnd%f_shelf .eq. 0._wp .and. lnd%f_shelf_old .gt. 0._wp ) then
        call end_cell_shelf(lnd%litter_c_shelf,lnd%litter_c13_shelf,lnd%litter_c14_shelf, &
                          lnd%fast_c_shelf,lnd%fast_c13_shelf,lnd%fast_c14_shelf, &
                          lnd%slow_c_shelf,lnd%slow_c13_shelf,lnd%slow_c14_shelf, &
                          lnd%soil_resp(ic_shelf),lnd%soil_resp13(ic_shelf),lnd%soil_resp14(ic_shelf), &
                          lnd%soil_c_tot(ic_shelf),lnd%soil_c13_tot(ic_shelf),lnd%soil_c14_tot(ic_shelf), &
                          lnd%ch4_emis_shelf,lnd%c13h4_emis_shelf)
       endif

       ! reset lake if lake grid cell fraction disappeard
       if( lnd%f_lake .eq. 0._wp .and. lnd%f_lake_old .gt. 0._wp ) then
        call end_cell_lake(lnd%t_skin,lnd%t_lake, &
                          lnd%w_snow(is_lake),lnd%h_snow(is_lake),lnd%mask_snow(is_lake), &
                          lnd%litter_c_lake,lnd%litter_c13_lake,lnd%litter_c14_lake, &
                          lnd%fast_c_lake,lnd%fast_c13_lake,lnd%fast_c14_lake, &
                          lnd%slow_c_lake,lnd%slow_c13_lake,lnd%slow_c14_lake, &
                          lnd%soil_resp(ic_lake),lnd%soil_resp13(ic_lake),lnd%soil_resp14(ic_lake), &
                          lnd%soil_c_tot(ic_lake),lnd%soil_c13_tot(ic_lake),lnd%soil_c14_tot(ic_lake), &
                          lnd%ch4_emis_lake,lnd%c13h4_emis_lake)
       endif

       ! update surface fractions
       call surface_frac_up(lnd%f_ice,lnd%f_ice_grd,lnd%f_shelf,lnd%f_lake,lnd%f_veg,lnd%pft_frac,lnd%frac_surf)

       ! update soil properties given the soil carbon content
       if( (.not.soil_par%constant_porosity .or. .not.soil_par%constant_soil_par_therm .or. .not.soil_par%constant_soil_par_hydro) .and. lnd%f_veg .gt. 0._wp ) then 
        call soil_par_update(lnd%f_peat,lnd%litter_c,lnd%fast_c,lnd%slow_c,lnd%litter_c_peat,lnd%acro_c,lnd%cato_c, &
                            mineral(i,j)%theta_sat,mineral(i,j)%k_sat,mineral(i,j)%psi_sat,mineral(i,j)%Bi,mineral(i,j)%lambda_s,mineral(i,j)%lambda_dry, &
                            lnd%frac_soc,lnd%theta_sat,lnd%psi_sat,lnd%k_sat,lnd%k_exp,lnd%psi_exp, &
                            lnd%theta_field,lnd%theta_wilt,lnd%lambda_s, &
                            lnd%lambda_dry)
       endif

       ! update vertical root distribution to account for active layer thickness
       if (veg_par%lroot_frac .and. lnd%f_veg.gt.0._wp) then
        call root_frac_update(lnd%alt,lnd%root_frac)
       endif

       ! update vertical litter input distribution to account for active layer thickness
       if (lnd%f_veg.gt.0._wp) then
        call litter_in_frac_update(lnd%alt,lnd%litter_in_frac)
       endif

      endif

      ! update phenological status of the vegetation 
      call phenology(lnd%frac_surf,lnd%f_veg, lnd%t2m, lnd%t2m_min_mon, lnd%gdd5_temp, lnd%gdd5, lnd%gdd, &
                    lnd%phen_acc, lnd%phen, lnd%gamma_leaf, lnd%lai, lnd%lai_bal, i,j)

      ! aerodynamic resistance
      call resist_aer(lnd%frac_surf,lnd%veg_h,lnd%lai,lnd%sai,lnd%h_snow,&
                     lnd%tatm,lnd%t_skin,lnd%wind, &
                     lnd%z0m,lnd%rough_m,lnd%rough_h, &
                     lnd%Ch,lnd%r_a,lnd%r_a_can,lnd%Ri)
 
      ! compute snow albedo, to be used for surface albedo calculation later
      call snow_albedo(lnd%f_veg,lnd%f_ice,lnd%f_lake,lnd%w_snow,lnd%w_snow_max, &
                      lnd%t_skin_veg,lnd%t_skin(i_ice),lnd%t_skin(i_lake), &
                      lnd%snow,lnd%dust_dep,lnd%coszm(doy), &
                      lnd%alb_snow_vis_dir,lnd%alb_snow_vis_dif, lnd%alb_snow_nir_dir, lnd%alb_snow_nir_dif, &
                      lnd%snow_grain, lnd%dust_con)

      ! compute surface albedo
      call surface_albedo(lnd%frac_surf,lnd%z_veg_std,lnd%h_snow,lnd%coszm(doy), &
                         lnd%lai,lnd%sai,lnd%z0m,lnd%f_snow_can,lnd%f_lake_ice, &
                         lnd%alb_snow_vis_dir,lnd%alb_snow_vis_dif,lnd%alb_snow_nir_dir,lnd%alb_snow_nir_dif, &
                         surf_par%alb_bare_vis(i,j),surf_par%alb_bare_nir(i,j), &
                         lnd%f_snow, &
                         lnd%alb_vis_dir,lnd%alb_vis_dif, &
                         lnd%alb_nir_dir,lnd%alb_nir_dif,lnd%albedo)

      if( lnd%f_veg .gt. 0._wp ) then
       ! do photosynthesis
       if (l_co2_fert_lim) then
         lndp%co2 = min(lndp%co2,co2_fert_lim_max)
         lndp%co2 = max(lndp%co2,co2_fert_lim_min)
       endif
       call photosynthesis(lndp%co2,lndp%c13_c12_atm,lndp%c14_c_atm, &
                          lnd%frac_surf,lnd%t2m,lnd%t2m_min_mon,lnd%gdd5,lnd%t_soil(1:nl), &
                          lnd%q2m,lnd%pressure,lnd%swnet,lnd%albedo,lnd%daylength(doy), &
                          lnd%theta_w,lnd%theta_field,lnd%theta_wilt,lnd%wilt,lnd%root_frac, &
                          lnd%lai,lnd%phen, &
                          lnd%leaf_c, lnd%stem_c, lnd%root_c, &
                          lnd%discrimination,lnd%ci,lnd%g_can, &
                          lnd%gpp,lnd%npp,lnd%npp13,lnd%npp14, &
                          lnd%npp_cum,lnd%npp13_cum,lnd%npp14_cum,lnd%npp_ann,lnd%npp13_ann,lnd%npp14_ann, lnd%aresp,i,j)
      endif

      ! surface resistance for evapotranspiration
      call resist_sur(lnd%frac_surf,lnd%mask_snow,lnd%w_snow, &
                     lnd%theta,lnd%theta_w,lnd%theta_field,lnd%theta_sat,lnd%psi_sat,lnd%psi_exp,lnd%k_exp,lnd%g_can, &
                     lnd%beta_s,lnd%r_s,lnd%beta_s_can,lnd%r_s_can)

      ! compute evaporation of intercepted water
      call canopy_water(lnd%frac_surf,lnd%lai,lnd%sai,lnd%r_a,lnd%t_skin, &
                       lnd%pressure,lnd%qatm,lnd%rain,lnd%snow, &
                       lnd%w_can,lnd%w_can_old,lnd%s_can,lnd%s_can_old, &
                       lnd%rain_ground,lnd%snow_ground,lnd%evap_can,lnd%subl_can,lnd%f_snow_can,lnd%f_wat_can)

      ! update soil thermal properties 
      if( lnd%f_veg .gt. 0._wp ) then
       call soil_par_thermal(lnd%t_soil,lnd%theta_w,lnd%theta_i,lnd%theta,lnd%theta_sat, &
                            lnd%lambda_s,lnd%lambda_dry, &
                            lnd%h_snow(is_veg), &
                            lnd%cap_soil,lnd%lambda_soil,lnd%lambda_int_soil)
      endif

      ! solve surface energy balance equation to get among others the ground heat flux 
      if( lnd%f_veg .gt. 0._wp ) then
       call ebal_veg(lnd%frac_surf,lnd%mask_snow(is_veg),lnd%h_snow(is_veg),lnd%w_snow(is_veg),lnd%lambda_soil, &
                    lnd%evap_can,lnd%subl_can, &
                    lnd%t_skin,lnd%t_skin_old, &
                    lnd%t_soil,lnd%tatm,lnd%qatm,lnd%pressure, &
                    lnd%swnet,lnd%swnet_min,lnd%lwdown, &
                    lnd%beta_s,lnd%r_s,lnd%beta_s_can,lnd%r_s_can,lnd%r_a,lnd%r_a_can, &
                    lnd%flx_g,lnd%dflxg_dT,lnd%flx_melt, &
                    lnd%flx_g_veg,lnd%dflxg_dT_veg,lnd%flx_melt_veg,lnd%t_skin_amp, &
                    lnd%num_lh,lnd%num_sh,lnd%num_sw,lnd%num_lw,lnd%denom_lh,lnd%denom_sh,lnd%denom_lw, &
                    lnd%f_sh,lnd%f_e,lnd%f_t,lnd%f_le,lnd%f_lt,lnd%f_lw,lnd%lh_ecan,lnd%qsat_e,lnd%dqsatdT_e,lnd%qsat_t,lnd%dqsatdT_t, &
                    lnd%energy_cons_surf1,i,j)
      endif

      ! update ice thermal properties 
      if( lnd%f_ice .gt. 0._wp ) then
       call ice_par_thermal(lnd%h_snow(is_ice), &
                            lnd%cap_ice,lnd%lambda_ice,lnd%lambda_int_ice)
      endif

      if( lnd%f_ice .gt. 0._wp ) then
       call ebal_ice(lnd%mask_snow(is_ice),lnd%h_snow(is_ice),lnd%w_snow(is_ice),lnd%lambda_ice, &
                    lnd%t_skin(i_ice),lnd%t_skin_old(i_ice), &
                    lnd%t_ice,lnd%tatm(i_ice),lnd%qatm(i_ice),lnd%pressure(i_ice), &
                    lnd%swnet(i_ice),lnd%swnet_min(i_ice),lnd%lwdown(i_ice), &
                    lnd%beta_s(i_ice),lnd%r_s(i_ice),lnd%r_a(i_ice), &
                    lnd%flx_g(i_ice),lnd%dflxg_dT(i_ice),lnd%flx_melt(i_ice),lnd%t_skin_amp(i_ice), &
                    lnd%num_lh(i_ice),lnd%num_sh(i_ice),lnd%num_sw(i_ice),lnd%num_lw(i_ice),lnd%denom_lh(i_ice),lnd%denom_sh(i_ice),lnd%denom_lw(i_ice), &
                    lnd%f_sh(i_ice),lnd%f_e(i_ice),lnd%f_le(i_ice),lnd%f_lw(i_ice),lnd%qsat_e(i_ice),lnd%dqsatdT_e(i_ice), &
                    lnd%energy_cons_surf1(i_ice),i,j)
      endif

      ! update lake thermal properties 
      if( lnd%f_lake .gt. 0._wp ) then
       call lake_par_thermal(lnd%h_snow(is_lake),lnd%h_lake,lnd%t_lake(1:nl_l),lnd%f_i_lake,lnd%wind(i_lake),lat(j), &
                            lnd%cap_lake,lnd%lambda_lake,lnd%lambda_int_lake)
      endif

      if( lnd%f_lake .gt. 0._wp ) then
       call ebal_lake(lnd%mask_snow(is_lake),lnd%h_snow(is_lake),lnd%w_snow(is_lake),lnd%lambda_lake, &
                    lnd%t_skin(i_lake),lnd%t_skin_old(i_lake), &
                    lnd%t_lake,lnd%tatm(i_lake),lnd%qatm(i_lake),lnd%pressure(i_lake), &
                    lnd%swnet(i_lake),lnd%swnet_min(i_lake),lnd%lwdown(i_lake), &
                    lnd%beta_s(i_lake),lnd%r_s(i_lake),lnd%r_a(i_lake), &
                    lnd%flx_g(i_lake),lnd%dflxg_dT(i_lake),lnd%flx_melt(i_lake),lnd%t_skin_amp(i_lake), &
                    lnd%num_lh(i_lake),lnd%num_sh(i_lake),lnd%num_sw(i_lake),lnd%num_lw(i_lake),lnd%denom_lh(i_lake),lnd%denom_sh(i_lake),lnd%denom_lw(i_lake), &
                    lnd%f_sh(i_lake),lnd%f_e(i_lake),lnd%f_le(i_lake),lnd%f_lw(i_lake),lnd%qsat_e(i_lake),lnd%dqsatdT_e(i_lake), &
                    lnd%energy_cons_surf1(i_lake),i,j)
      endif

      ! compute new soil temperature profile, including latent heat effects of water phase changes
      if( lnd%f_veg .gt. 0._wp ) then
       call soil_temp(lnd%mask_snow(is_veg),lnd%h_snow(is_veg), &
                     lnd%cap_soil,lnd%lambda_int_soil, &
                     lnd%psi_sat,lnd%psi_exp,lnd%theta_sat,lnd%soil_resp_l,lnd%f_peat, &
                     lnd%flx_g_veg,lnd%dflxg_dT_veg,lnd%flx_melt_veg, &
                     lnd%t_soil,lnd%t_soil_cum, &
                     lnd%w_snow(is_veg),lnd%w_w,lnd%w_i, &
                     lnd%theta_w,lnd%theta_i, &
                     lnd%snowmelt(is_veg),lnd%t_soil_old,lnd%w_snow_old(is_veg), &
                     lnd%w_w_old,lnd%w_i_old, &
                     lnd%w_w_phase,lnd%w_i_phase, &
                     lnd%energy_cons_soil,i,j)
      endif

      if( lnd%f_ice .gt. 0._wp ) then
       call ice_temp(lnd%mask_snow(is_ice),lnd%h_snow(is_ice),lnd%cap_ice,lnd%lambda_int_ice,lnd%flx_g(i_ice),lnd%dflxg_dT(i_ice),lnd%flx_melt(i_ice), &
                    lnd%t_ice,lnd%w_snow(is_ice),lnd%snowmelt(is_ice),lnd%icemelt(is_ice),lnd%t_ice_old,lnd%w_snow_old(is_ice), &
                    lnd%energy_cons_ice,i,j)
      endif

      if( lnd%f_lake .gt. 0._wp ) then
       call lake_temp(lnd%mask_snow(is_lake),lnd%h_snow(is_lake),lnd%h_lake, &
                     lnd%cap_lake,lnd%lambda_int_lake, &
                     lnd%flx_g(i_lake),lnd%dflxg_dT(i_lake),lnd%flx_melt(i_lake), &
                     lnd%t_lake,lnd%w_snow(is_lake),lnd%w_w_lake,lnd%w_i_lake,lnd%f_i_lake,lnd%f_lake_ice,lnd%snowmelt(is_lake), &
                     lnd%t_lake_old,lnd%w_snow_old(is_lake), &
                     lnd%h_lake_conv, lnd%h_lake_mix, &
                     lnd%energy_cons_lake,i,j)

       ! use bottom layer lake temperature as top boundary condition for sublake temperature profile
       lnd%t_sublake(0) = lnd%t_lake(nl_l)

       call sublake_par_thermal(lnd%theta_w_sublake,lnd%theta_i_sublake,lnd%theta_sat,lnd%lambda_s, &
                             lnd%cap_sublake,lnd%lambda_int_sublake)

       call sublake_temp(lnd%cap_sublake,lnd%lambda_int_sublake, &
                      lnd%psi_sat,lnd%psi_exp,lnd%theta_sat,lnd%soil_resp_l, &
                      lnd%t_sublake,lnd%w_w_sublake,lnd%w_i_sublake,lnd%theta_w_sublake,lnd%theta_i_sublake, &
                      lnd%t_sublake_cum,lnd%theta_w_sublake_cum,lnd%theta_i_sublake_cum, &
                      lnd%energy_cons_lake,i,j)
      endif

      if( lnd%f_shelf .gt. 0._wp ) then

       call shelf_par_thermal(lnd%theta_w_shelf,lnd%theta_i_shelf,lnd%theta_sat,lnd%lambda_s, &
                             lnd%cap_shelf,lnd%lambda_int_shelf)

       call shelf_temp(lnd%cap_shelf,lnd%lambda_int_shelf, &
                      lnd%psi_sat,lnd%psi_exp,lnd%theta_sat,lnd%soil_resp_l, &
                      lnd%t_shelf,lnd%w_w_shelf,lnd%w_i_shelf,lnd%theta_w_shelf,lnd%theta_i_shelf, &
                      lnd%t_shelf_cum,lnd%theta_w_shelf_cum,lnd%theta_i_shelf_cum, &
                      lnd%energy_cons_shelf,i,j)
      endif

      ! update skin temperature given the new top soil temperature and diagnose all surface energy fluxes
      if( lnd%f_veg .gt. 0._wp ) then
       call update_tskin_veg(lnd%frac_surf,lnd%mask_snow(is_veg),lnd%t_skin_old,lnd%dflxg_dT, &
                       lnd%tatm,lnd%qatm,lnd%swnet,lnd%lwdown, &
                       lnd%t_soil,lnd%t_soil_old,lnd%evap_can,lnd%subl_can, &
                       lnd%flx_g,lnd%flx_melt,lnd%t_skin,lnd%t_skin_veg, &
                       lnd%flx_sh,lnd%flx_lwu,lnd%lwnet,lnd%flx_lh, &
                       lnd%evap_surface,lnd%transpiration,lnd%et, &
                       lnd%num_lh,lnd%num_sh,lnd%num_sw,lnd%num_lw,lnd%denom_lh,lnd%denom_sh,lnd%denom_lw, &
                       lnd%f_sh,lnd%f_e,lnd%f_t,lnd%f_le,lnd%f_lt,lnd%f_lw,lnd%lh_ecan,lnd%qsat_e,lnd%dqsatdT_e,lnd%qsat_t,lnd%dqsatdT_t, &
                       lnd%energy_cons_surf2)

      endif

      if( lnd%f_ice .gt. 0._wp ) then
       call update_tskin_ice(lnd%mask_snow(is_ice),lnd%t_skin_old(i_ice), &
                       lnd%dflxg_dT(i_ice), &
                       lnd%tatm(i_ice),lnd%qatm(i_ice),lnd%swnet(i_ice),lnd%lwdown(i_ice), &
                       lnd%t_ice,lnd%t_ice_old, &
                       lnd%flx_g(i_ice),lnd%flx_melt(i_ice), &
                       lnd%t_skin(i_ice),lnd%flx_sh(i_ice), &
                       lnd%flx_lwu(i_ice),lnd%flx_lh(i_ice), &
                       lnd%evap_surface(i_ice),lnd%et(i_ice), &
                       lnd%num_lh(i_ice),lnd%num_sh(i_ice),lnd%num_sw(i_ice),lnd%num_lw(i_ice), &
                       lnd%denom_lh(i_ice),lnd%denom_sh(i_ice),lnd%denom_lw(i_ice), &
                       lnd%f_sh(i_ice),lnd%f_e(i_ice),lnd%f_le(i_ice),lnd%f_lw(i_ice),lnd%qsat_e(i_ice),lnd%dqsatdT_e(i_ice), &
                       lnd%energy_cons_surf2(i_ice),i,j)
      endif

      if( lnd%f_lake .gt. 0._wp ) then
       call update_tskin_lake(lnd%mask_snow(is_lake),lnd%t_skin_old(i_lake), &
                       lnd%dflxg_dT(i_lake), &
                       lnd%tatm(i_lake),lnd%qatm(i_lake),lnd%swnet(i_lake),lnd%lwdown(i_lake), &
                       lnd%t_lake,lnd%t_lake_old, &
                       lnd%flx_g(i_lake),lnd%flx_melt(i_lake), &
                       lnd%t_skin(i_lake),lnd%flx_sh(i_lake), &
                       lnd%flx_lwu(i_lake),lnd%flx_lh(i_lake), &
                       lnd%evap_surface(i_lake),lnd%et(i_lake), &
                       lnd%num_lh(i_lake),lnd%num_sh(i_lake),lnd%num_sw(i_lake),lnd%num_lw(i_lake), &
                       lnd%denom_lh(i_lake),lnd%denom_sh(i_lake),lnd%denom_lw(i_lake), &
                       lnd%f_sh(i_lake),lnd%f_e(i_lake),lnd%f_le(i_lake),lnd%f_lw(i_lake),lnd%qsat_e(i_lake),lnd%dqsatdT_e(i_lake), &
                       lnd%energy_cons_surf2(i_lake),i,j)
      endif

      ! surface hydrology
      call surface_hydrology(lnd%frac_surf,lnd%mask_snow, &
                            lnd%evap_surface,lnd%rain_ground,lnd%snow_ground, &
                            lnd%snowmelt, lnd%icemelt, &
                            lnd%theta, &
                            lnd%theta_sat,lnd%theta_field,lnd%k_sat, &
                            lnd%cap_soil(1),lnd%cap_ice(1),lnd%cap_lake(1),lnd%f_peat,lnd%w_table_peat, &
                            topmodel(i,j)%cti_mean,topmodel(i,j)%cti_cdf, &
                            dyptop(i,j)%k,dyptop(i,j)%v,dyptop(i,j)%xm,dyptop(i,j)%fmax, &
                            lnd%w_snow_old,lnd%w_snow,lnd%w_snow_max,lnd%w_w,lnd%w_i, &
                            lnd%w_table_cum,lnd%f_wet_cum,lnd%tatm,lnd%t_soil,lnd%t_ice,lnd%t_lake, &
                            lnd%h_snow,lnd%calving,lnd%runoff_sur, &
                            lnd%infiltration,lnd%w_table,lnd%f_wet,lnd%f_wet_max,lnd%cti_lim,lnd%lake_water_tendency)
 
      ! if any soil
      if( lnd%f_veg .gt. 0._wp ) then
       ! update hydrological soil properties if any soil
       call soil_par_hydro(lnd%theta_w,lnd%theta_sat, &
                          lnd%w_w,lnd%w_i,lnd%psi_sat,lnd%k_sat,lnd%k_exp,lnd%psi_exp, &
                          lnd%theta, &
                          lnd%psi,lnd%kappa_int)

       ! update water content of the soil profile
       call soil_hydro(lnd%frac_surf,lnd%mask_snow(is_veg),lnd%theta_sat, &
                      lnd%k_sat,lnd%k_exp,lnd%psi_exp,lnd%kappa_int,lnd%psi, &
                      lnd%w_table, &
                      lnd%transpiration,lnd%evap_surface,lnd%infiltration,lnd%wilt, &
                      lnd%w_snow(is_veg),lnd%w_w,lnd%w_i, &
                      lnd%theta_w,lnd%theta_i,lnd%theta,lnd%theta_w_cum,lnd%theta_i_cum,lnd%theta_fire_cum, &
                      lnd%drainage)

      else
       lnd%drainage = 0._wp
      endif

      ! total runoff over gridcell
      lnd%runoff(is_veg) = lnd%runoff_sur(is_veg)+lnd%drainage(is_veg)   ! runoff from vegetated grid cell fraction
      lnd%runoff(is_ice) = lnd%runoff_sur(is_ice)+lnd%drainage(is_ice)     ! runoff from glacier grid cell fraction
      lnd%runoff(is_lake)= lnd%runoff_sur(is_lake)+lnd%drainage(is_lake)  ! runoff from lake cell fraction, fixme

      ! cumulate runoff over the year for weathering, vegetated part only fixme?
      lnd%runoff_ann = lnd%runoff_ann + lnd%runoff(is_veg)*dt ! kg/m2

      if( check_water ) then
      if( (lnd%f_veg .gt. 0._wp .and. lnd%f_veg_old .gt. 0._wp) .or. (lnd%f_ice .gt. 0._wp .and. lnd%f_ice_old .gt. 0._wp) ) then
       ! water conservation check
       call water_check(i,j,lnd%frac_surf,lnd%f_veg,lnd%rain,lnd%snow, &
                         lnd%et,lnd%runoff_sur,lnd%calving,lnd%drainage,lnd%icemelt, &
                         lnd%w_w,lnd%w_i,lnd%w_snow,lnd%w_can,lnd%s_can, &
                         lnd%w_w_old,lnd%w_i_old,lnd%w_snow_old,lnd%w_can_old,lnd%s_can_old,lnd%lake_water_tendency, &
                         lnd%water_cons)
      endif
      endif

      if( time_call_carb_p ) then

       if (lnd%f_veg.gt.0._wp) then
        ! update soil carbon decomposition properties
        call soil_carbon_par(lnd%f_veg,lnd%theta_field,lnd%theta_sat, &
                            lnd%litter_c_peat,lnd%acro_c,lnd%cato_c,lnd%dust_dep, &
                            lnd%t_soil_cum,lnd%theta_w_cum,lnd%theta_i_cum,lnd%psi, &
                            lnd%f_wet_cum,lnd%w_table_cum, &
                            lnd%f_wet_mon,lnd%f_wet_long,lnd%w_table_mon, &
                            lnd%t_soil_max, &
                            lnd%ftemp,lnd%fmoist,lnd%fdepth, &
                            lnd%k_litter,lnd%k_fast,lnd%k_slow,lnd%diff_soilc,lnd%adv_soilc, &
                            lnd%k_litter_wet,lnd%k_fast_wet,lnd%k_slow_wet, &
                            lnd%k_litter_peat,lnd%k_acro,lnd%k_cato, &
                            lnd%k_litter_peat_anox,lnd%k_acro_anox,lnd%ch4_frac_wet,lnd%ch4_frac_peat, &
                            lnd%f_peat_pot,lnd%f_oxic_peat,lnd%f_wetland, &
                            lnd%w_table_min,lnd%w_table_peat,lnd%alt, &
                            lnd%acro_h,lnd%cato_h,lnd%peat_c_ini_year)
       endif

       if (lnd%f_shelf.gt.0._wp) then
        ! update shelf carbon decomposition properties
        call shelf_carbon_par(lnd%t_shelf_cum,lnd%theta_w_shelf_cum,lnd%theta_i_shelf_cum, &
                              lnd%k_litter_shelf,lnd%k_fast_shelf,lnd%k_slow_shelf,lnd%diff_shelfc,lnd%adv_shelfc, &
                              lnd%ch4_frac_shelf)
       endif

       if (lnd%f_ice_grd.gt.0._wp) then
        call ice_carbon_par(lnd%k_litter_ice,lnd%k_fast_ice,lnd%k_slow_ice,lnd%diff_icec,lnd%adv_icec)
       endif

       if (lnd%f_lake.gt.0._wp) then
        ! update lake carbon decomposition properties
        call lake_carbon_par(lnd%t_sublake_cum,lnd%theta_w_sublake_cum,lnd%theta_i_sublake_cum, &
                              lnd%k_litter_lake,lnd%k_fast_lake,lnd%k_slow_lake,lnd%diff_lakec,lnd%adv_lakec, &
                              lnd%ch4_frac_lake)
       endif

      endif


      if( time_eom_lnd ) then
       call dynveg_par(lnd%disturbance,lnd%t2m_min_mon,lnd%gdd5,lnd%veg_c_above,lnd%theta_fire_cum,lnd%gamma_dist_cum)
      endif

      if( time_call_veg ) then
       ! update vegetation carbon and vegetation distribution
       call dyn_veg(lndp%co2,lnd%f_veg,lnd%f_veg_old,lnd%f_ice_grd,lnd%f_ice_grd_old,lnd%f_ice_nbr, &
                   lnd%f_lake,lnd%f_lake_old,lnd%f_shelf,lnd%f_shelf_old, lnd%f_crop, lnd%f_pasture, &
                   lnd%gamma_luc, lnd%z_veg_std, lnd%gamma_ice, lnd%gamma_dist,lnd%gamma_dist_cum, &
                   lnd%npp_ann,lnd%npp13_ann,lnd%npp14_ann, &
                   lnd%gamma_leaf, lnd%lambda, &
                   lnd%lai_bal,lnd%sai,lnd%root_frac,lnd%litter_in_frac, &
                   lnd%veg_c,lnd%veg_c13,lnd%veg_c14, &
                   lnd%leaf_c, &
                   lnd%stem_c, &
                   lnd%root_c,lnd%seed_frac,lnd%pft_frac, &
                   lnd%veg_c_above,lnd%veg_c13_above,lnd%veg_c14_above,lnd%veg_c_below,lnd%veg_c13_below,lnd%veg_c14_below, &
                   lnd%veg_h, &
                   lnd%litterfall,lnd%litterfall13,lnd%litterfall14, &
                   lnd%npp_real,lnd%npp13_real,lnd%npp14_real, &
                   lnd%carbon_cons_veg,lnd%carbon13_cons_veg,lnd%carbon14_cons_veg,i,j)
        ! update surface fractions
        call surface_frac_up(lnd%f_ice,lnd%f_ice_grd,lnd%f_shelf,lnd%f_lake,lnd%f_veg,lnd%pft_frac,lnd%frac_surf)
       endif

       if( time_call_carb ) then 
        ! update soil carbon content
        if( lnd%f_veg .gt. 0._wp ) then
         call soil_carbon(lnd%f_veg,lnd%f_wetland,lnd%f_peat, &
                        lnd%litterfall(:,ic_min),lnd%litterfall13(:,ic_min),lnd%litterfall14(:,ic_min), &
                        lnd%litter_c,lnd%fast_c,lnd%slow_c, &
                        lnd%litter_c13,lnd%fast_c13,lnd%slow_c13, &
                        lnd%litter_c14,lnd%fast_c14,lnd%slow_c14, &
                        lnd%k_litter,lnd%k_fast,lnd%k_slow, &
                        lnd%k_litter_wet,lnd%k_fast_wet,lnd%k_slow_wet,lnd%diff_soilc,lnd%adv_soilc,lnd%ch4_frac_wet, &
                        lnd%soil_resp(ic_min),lnd%soil_resp13(ic_min),lnd%soil_resp14(ic_min),lnd%soil_resp_l(:,ic_min), &
                        lnd%soil_c_tot(ic_min),lnd%soil_c13_tot(ic_min),lnd%soil_c14_tot(ic_min), &
                        lnd%ch4_emis_wetland,lnd%c13h4_emis_wetland, &
                        lnd%carbon_cons_soil(ic_min),lnd%carbon13_cons_soil(ic_min), &
                        lnd%carbon14_cons_soil(ic_min))
        endif

        if( lnd%f_veg .gt. 0._wp .and. peat_par%peat_carb ) then
         ! update peatland carbon content
         call peat_carbon(lnd%f_oxic_peat, &
                        lnd%litterfall(:,ic_peat),lnd%litterfall13(:,ic_peat),lnd%litterfall14(:,ic_peat), &
                        lnd%litter_c_peat,lnd%acro_c,lnd%cato_c, &
                        lnd%litter_c13_peat,lnd%acro_c13,lnd%cato_c13, &
                        lnd%litter_c14_peat,lnd%acro_c14,lnd%cato_c14, &
                        lnd%k_litter_peat,lnd%k_acro,lnd%k_cato, &
                        lnd%k_litter_peat_anox,lnd%k_acro_anox,lnd%ch4_frac_peat, &
                        lnd%soil_resp(ic_peat),lnd%soil_resp13(ic_peat),lnd%soil_resp14(ic_peat),lnd%soil_resp_l(:,ic_peat), &
                        lnd%soil_c_tot(ic_peat),lnd%soil_c13_tot(ic_peat),lnd%soil_c14_tot(ic_peat), &
                        lnd%ch4_emis_peat,lnd%c13h4_emis_peat, &
                        lnd%carbon_cons_soil(ic_peat),lnd%carbon13_cons_soil(ic_peat),lnd%carbon14_cons_soil(ic_peat), &
                        lnd%peat_c_ini_year,lnd%dCpeat_dt)
        endif

        if( lnd%f_shelf .gt. 0._wp ) then
         call shelf_carbon(lnd%litterfall(:,ic_shelf),lnd%litterfall13(:,ic_shelf),lnd%litterfall14(:,ic_shelf), &
                         lnd%ch4_frac_shelf, &
                         lnd%litter_c_shelf,lnd%fast_c_shelf,lnd%slow_c_shelf, &
                         lnd%litter_c13_shelf,lnd%fast_c13_shelf,lnd%slow_c13_shelf, &
                         lnd%litter_c14_shelf,lnd%fast_c14_shelf,lnd%slow_c14_shelf, &
                         lnd%k_litter_shelf,lnd%k_fast_shelf,lnd%k_slow_shelf,lnd%diff_shelfc,lnd%adv_shelfc, &
                         lnd%soil_resp(ic_shelf),lnd%soil_resp13(ic_shelf),lnd%soil_resp14(ic_shelf),lnd%soil_resp_l(:,ic_shelf), &
                         lnd%soil_c_tot(ic_shelf),lnd%soil_c13_tot(ic_shelf),lnd%soil_c14_tot(ic_shelf), &
                         lnd%ch4_emis_shelf,lnd%c13h4_emis_shelf, &
                         lnd%carbon_cons_soil(ic_shelf),lnd%carbon13_cons_soil(ic_shelf),lnd%carbon14_cons_soil(ic_shelf))
        endif

        if( lnd%f_ice_grd .gt. 0._wp ) then
         call ice_carbon(lnd%litterfall(:,ic_ice),lnd%litterfall13(:,ic_ice),lnd%litterfall14(:,ic_ice), &
                          lnd%litter_c_ice,lnd%fast_c_ice,lnd%slow_c_ice, &
                          lnd%litter_c13_ice,lnd%fast_c13_ice,lnd%slow_c13_ice, &
                          lnd%litter_c14_ice,lnd%fast_c14_ice,lnd%slow_c14_ice, &
                          lnd%k_litter_ice,lnd%k_fast_ice,lnd%k_slow_ice,lnd%diff_icec,lnd%adv_icec, &
                          lnd%soil_resp(ic_ice),lnd%soil_resp13(ic_ice),lnd%soil_resp14(ic_ice), &
                          lnd%soil_c_tot(ic_ice),lnd%soil_c13_tot(ic_ice),lnd%soil_c14_tot(ic_ice), &
                          lnd%carbon_cons_soil(ic_ice),lnd%carbon13_cons_soil(ic_ice),lnd%carbon14_cons_soil(ic_ice))
        endif

        if( lnd%f_lake .gt. 0._wp ) then
         call lake_carbon(lnd%litterfall(:,ic_lake),lnd%litterfall13(:,ic_lake),lnd%litterfall14(:,ic_lake), &
                         lnd%ch4_frac_lake, &
                         lnd%litter_c_lake,lnd%fast_c_lake,lnd%slow_c_lake, &
                         lnd%litter_c13_lake,lnd%fast_c13_lake,lnd%slow_c13_lake, &
                         lnd%litter_c14_lake,lnd%fast_c14_lake,lnd%slow_c14_lake, &
                         lnd%k_litter_lake,lnd%k_fast_lake,lnd%k_slow_lake,lnd%diff_lakec,lnd%adv_lakec, &
                         lnd%soil_resp(ic_lake),lnd%soil_resp13(ic_lake),lnd%soil_resp14(ic_lake),lnd%soil_resp_l(:,ic_lake), &
                         lnd%soil_c_tot(ic_lake),lnd%soil_c13_tot(ic_lake),lnd%soil_c14_tot(ic_lake), &
                         lnd%ch4_emis_lake,lnd%c13h4_emis_lake, &
                         lnd%carbon_cons_soil(ic_lake),lnd%carbon13_cons_soil(ic_lake),lnd%carbon14_cons_soil(ic_lake))
        endif

        if (minval(lnd%litter_c).lt.0._wp) then
          print *,'litter<0',doy,year,i,j
          print *,lnd%litter_c
          print *,lnd%litterfall(:,ic_min)*86400.*30.
          print *,lnd%f_veg,lnd%f_veg_old
          print *,lnd%f_ice_grd,lnd%f_ice_grd_old
          print *,lnd%f_shelf,lnd%f_shelf_old
          print *,lnd%alt
          !if (minval(lnd%litter_c).lt.-1e-10) stop
        endif
        if (lnd%litter_c_peat.lt.0._wp) then
          print *,lnd%litter_c_peat
          print *,lnd%litterfall(:,ic_peat)
          !stop
        endif
        if (minval(lnd%fast_c).lt.0._wp) then
          print *,'fast_c<0',doy,year,i,j
          print *,lnd%fast_c
        endif
        if (minval(lnd%slow_c).lt.0._wp) then
          print *,'slow_c<0',doy,year,i,j
          print *,lnd%slow_c
        endif

       endif


       ! dust emissions
       if (lnd%f_veg.gt.0._wp) then
         call dust_emission(lnd%frac_surf,lnd%z_veg_std,lnd%z_veg,lnd%z_veg_min,lnd%z_veg_max,lnd%f_snow(i_bare), &
                        lnd%lai,lnd%sai,lnd%h_snow(is_veg),lnd%t_skin,lnd%tatm,lnd%theta_w(1),lnd%theta_i(1),lnd%wind(i_bare), &
                        lnd%dust_emis_d,lnd%dust_emis_g,lnd%dust_emis_s,lnd%dust_emis)
       else
         lnd%dust_emis_d = 0._wp
         lnd%dust_emis_g = 0._wp
         lnd%dust_emis_s = 0._wp
         lnd%dust_emis = 0._wp
       endif


       if (time_eoy_lnd) then
       
         ! weathering
         if (i_weathering.eq.1) then
           ! GEM-CO2 model

           ! update lithological map over shelf areas
           if (lnd%f_land.gt.0._wp) then
             if (lnd%f_land.gt.lnd%f_land0) then
               ! more land than at present sea level
               do n=1,weath_gemco2_par%nlit
                 if (n==1) then
                   ! carbonates
                   lnd%lithology_gemco2(n) = (lnd%f_land0*weath_gemco2_par%lit_map(i,j,n) &
                     + (lnd%f_land-lnd%f_land0)*lnd%f_carb) / lnd%f_land
                 else
                   ! equally share other lithologies
                   lnd%lithology_gemco2(n) = (lnd%f_land0*weath_gemco2_par%lit_map(i,j,n) &
                     + (lnd%f_land-lnd%f_land0)*(1._wp-lnd%f_carb)/real(weath_gemco2_par%nlit-1,wp)) / lnd%f_land
                 endif
               enddo
             else
               lnd%lithology_gemco2 = weath_gemco2_par%lit_map(i,j,:)
             endif
           else
             lnd%lithology_gemco2(:) = 0._wp
           endif

           call weathering_gemco2(lndp%c13_c12_atm,lndp%c14_c_atm,lndp%weath_scale,lnd%f_veg,lnd%lithology_gemco2,lnd%runoff_ann, &
             lnd%weath_carb,lnd%weath13_carb,lnd%weath14_carb, &
             lnd%weath_sil,lnd%weath13_sil,lnd%weath14_sil, lnd%weath_loess)

         else if (i_weathering.eq.2) then
           ! UHH (PalMod) model

           w_lgm = max(0._wp,min(1._wp,year_now/(-21000._wp)))

           ! update lithological map over shelf areas
           if (lnd%f_land.gt.0._wp) then
             if (lnd%f_land.gt.lnd%f_land0) then
               ! more land than at present sea level
               do n=1,weath_uhh_par%nlit
                 ! average between fractions for present day land and shelf
                 lnd%lithology_uhh(n) = (lnd%f_land0*weath_uhh_par%lit_map(i,j,n) &
                 !lnd%lithology_uhh(n) = (lnd%f_land0*((1._wp-w_lgm)*weath_uhh_par%lit_map(i,j,n) + w_lgm*weath_uhh_par%lit_map_lgm(i,j,n)) &
                   + (lnd%f_land-lnd%f_land0)*weath_uhh_par%lit_map_shelf(i,j,n)) / lnd%f_land
               enddo
             else
               lnd%lithology_uhh = weath_uhh_par%lit_map(i,j,:)
               !lnd%lithology_uhh = (1._wp-w_lgm)*weath_uhh_par%lit_map(i,j,:) + w_lgm*weath_uhh_par%lit_map_lgm(i,j,:)
             endif
           else
             lnd%lithology_uhh(:) = 0._wp
           endif

           call weathering_uhh(lndp%c13_c12_atm,lndp%c14_c_atm,lndp%weath_scale,lnd%f_veg,lnd%lithology_uhh,lnd%runoff_ann,lnd%t2m_ann_mean, &
             lnd%weath_carb,lnd%weath13_carb,lnd%weath14_carb, &
             lnd%weath_sil,lnd%weath13_sil,lnd%weath14_sil, lnd%weath_loess)

         endif

         ! POC and DOC carbon export from land through rivers
         if (l_river_export) then

           call carbon_export(lnd%runoff_ann, &
             lnd%f_veg,lnd%f_peat, &
             lnd%litter_c,lnd%litter_c13,lnd%litter_c14,lnd%fast_c,lnd%fast_c13,lnd%fast_c14,lnd%slow_c,lnd%slow_c13,lnd%slow_c14, &
             lnd%litter_c_peat,lnd%litter_c13_peat,lnd%litter_c14_peat,lnd%acro_c,lnd%acro_c13,lnd%acro_c14,lnd%cato_c,lnd%cato_c13,lnd%cato_c14, &
             lnd%poc_export,lnd%poc13_export,lnd%poc14_export, &
             lnd%doc_export,lnd%doc13_export,lnd%doc14_export)

         else

           lnd%poc_export   = 0._wp  
           lnd%poc13_export = 0._wp  
           lnd%poc14_export = 0._wp  
           lnd%doc_export   = 0._wp  
           lnd%doc13_export = 0._wp  
           lnd%doc14_export = 0._wp  

         endif

         ! reset cumulated runoff
         lnd%runoff_ann = 0._wp

       endif


    return

end subroutine lnd_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l n d _ i n i t 
  !   Purpose    :  land initialisation
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lnd_init(lnd,f_lnd,f_ocn,f_ice,f_ice_grd,f_lake)

    use surface_par_lnd

    implicit none 

    type(lnd_class) :: lnd

    real(wp), intent(in) :: f_lnd(:,:)
    real(wp), intent(in) :: f_ocn(:,:)
    real(wp), intent(in) :: f_ice(:,:)
    real(wp), intent(in) :: f_ice_grd(:,:)
    real(wp), intent(in) :: f_lake(:,:)

    integer :: i, j, n, k
    real(wp), allocatable :: f_veg(:,:)


    ! initialize grid
    call lnd_grid_init

    lnd%z_lake(:) = z_l(1:nl_l)

    ! initialize parameters
    call lnd_params_init

    ! allocate
    allocate(lnd%l2d(ni,nj))
    call lnd_alloc(lnd) 

    ! initialise soil properties
    do i = 1,ni
      do j = 1,nj
        do k = 1,nl
          ! mineral + organic soil characteristics
          lnd%l2d(i,j)%theta_sat(k) = mineral(i,j)%theta_sat
          ! hydraulic properties
          lnd%l2d(i,j)%psi_sat(k) = mineral(i,j)%psi_sat
          lnd%l2d(i,j)%k_sat(k) = mineral(i,j)%k_sat
          lnd%l2d(i,j)%k_exp(k) = 2*nint(mineral(i,j)%Bi)+3
          lnd%l2d(i,j)%psi_exp(k) = -nint(mineral(i,j)%Bi)
          ! field capacity
          lnd%l2d(i,j)%theta_field(k) = lnd%l2d(i,j)%theta_sat(k) * (0.1_wp / (86400._wp * lnd%l2d(i,j)%k_sat(k)))**(1._wp/lnd%l2d(i,j)%k_exp(k))
          ! wilting point
          lnd%l2d(i,j)%theta_wilt(k) = lnd%l2d(i,j)%theta_sat(k) * (hydro_par%p_psi_min / lnd%l2d(i,j)%psi_sat(k))**(1._wp/lnd%l2d(i,j)%psi_exp(k))
          ! thermal properties
          lnd%l2d(i,j)%lambda_s(k) = mineral(i,j)%lambda_s
          lnd%l2d(i,j)%lambda_dry(k) = mineral(i,j)%lambda_dry
        enddo
      enddo
    enddo

    if (soil_par%uniform_porosity) then
      do i = 1,ni
        do j = 1,nj
          do k=1,nl
            lnd%l2d(i,j)%theta_sat(k) = soil_par%theta_sat_u
          enddo
        enddo
      enddo
    endif

    if (soil_par%uniform_soil_par_therm) then
      do i = 1,ni
        do j = 1,nj
          do k=1,nl
            lnd%l2d(i,j)%lambda_s(k) = soil_par%lambda_s_u
            lnd%l2d(i,j)%lambda_dry(k) = soil_par%lambda_dry_u
          enddo
        enddo
      enddo
    endif

    if (soil_par%uniform_soil_par_hydro) then
      do i = 1,ni
        do j = 1,nj
          do k=1,nl
            lnd%l2d(i,j)%psi_sat(k) = soil_par%psi_sat_u
            lnd%l2d(i,j)%k_sat(k) = soil_par%k_sat_u
            lnd%l2d(i,j)%k_exp(k) = 2*soil_par%b_u+3
            lnd%l2d(i,j)%psi_exp(k) = -soil_par%b_u
            lnd%l2d(i,j)%theta_field(k) = soil_par%theta_field_u
            lnd%l2d(i,j)%theta_wilt(k) = soil_par%theta_wilt_u
          enddo
        enddo
      enddo
    endif


    allocate(f_veg(ni,nj))

    if( lnd_restart ) then

      call lnd_read_restart(trim(restart_in_dir)//"/lnd_restart.nc",lnd%l2d,lnd%l0d)

      ! update surface fractions given the prescribed input fractions of land/ocean/ice/lake
      ! potentially vegetated fraction
      where (f_lnd.gt.0._wp) 
        f_veg = f_lnd - f_ice_grd - f_lake
      elsewhere
        f_veg = 0._wp
      endwhere
      ! f_veg can be negative for numerical reasons
      where (f_veg.lt.0._wp)
        f_veg = 0._wp
      endwhere
      ! if prescribed vegetation read PFT fractions from file
      if (.not.l_dynveg) then
        do i=1,nx
          do j=1,ny
            call nc_read(trim(pft_fix_file),"pft_frac",lnd%l2d(i,j)%pft_frac(:),start=[1,i,j],count=[npft,1,1])
          enddo
        enddo
      endif
      ! if fixed LAI, read balanced LAI from file
      if (l_fixlai) then
        do i=1,nx
          do j=1,ny
            call nc_read(trim(lai_fix_file),"lai_bal",lnd%l2d(i,j)%lai_bal(:),start=[1,i,j],count=[npft,1,1])
          enddo
        enddo
      endif
      ! get surface fractions
      do i = 1,ni
        do j = 1,nj
          call surface_frac_up(f_ice(i,j),f_ice_grd(i,j),f_ocn(i,j),f_lake(i,j),f_veg(i,j),lnd%l2d(i,j)%pft_frac,lnd%l2d(i,j)%frac_surf)
        enddo
      enddo

    else ! initialize

      f_veg = f_lnd - f_ice_grd - f_lake
      f_veg = max(0._wp,f_veg)

      !lnd%l2d%f_veg = 0._wp
      !lnd%l2d%f_ice = 0._wp
      !lnd%l2d%f_ice_grd = 0._wp
      !lnd%l2d%f_shelf = 0._wp
      !lnd%l2d%f_lake = 0._wp
      lnd%l2d%f_veg = f_veg 
      lnd%l2d%f_ice = f_ice 
      lnd%l2d%f_ice_grd = f_ice_grd
      lnd%l2d%f_shelf = f_ocn
      lnd%l2d%f_lake = f_lake

      lnd%l2d%mask_lnd = 1  ! land everywhere, fixme, find smart way

      lnd%l2d%alt = -1._wp

      lnd%l2d%t2m_min_mon = T0+20._wp

      lnd%l2d(:,:)%Cflx_atm_lnd = 0._wp
      lnd%l2d(:,:)%C13flx_atm_lnd = 0._wp
      lnd%l2d(:,:)%C14flx_atm_lnd = 0._wp

      lnd%l0d%landc = 0._wp
      lnd%l0d%landc13 = 0._wp
      lnd%l0d%landc14 = 0._wp

      lnd%l0d%burc = 0._wp
      lnd%l0d%burc13 = 0._wp
      lnd%l0d%burc14 = 0._wp

      lnd%l0d%weath_scale = 1._wp

      ! seeds everywhere, overwrite nml value
      veg_par%iseed = 1

      ! initialize prognostic variables
      do i = 1,ni
        do j = 1,nj
          lnd%l2d(i,j)%frac_surf(:) = 0._wp
          if (i_init_veg.eq.1) then
            ! bare soil everywhere
            lnd%l2d(i,j)%frac_surf(i_bare) = f_lnd(i,j)
            lnd%l2d(i,j)%pft_frac        = veg_par%seed_fraction
            lnd%l2d(i,j)%lai_bal         = pft_par%lai_min
          else if (i_init_veg.eq.2) then
            ! forest everywhere
            !lnd%l2d(i,j)%frac_surf(1) = 0.5_wp*f_lnd(i,j)   ! tropical
            !lnd%l2d(i,j)%frac_surf(2) = 0.5_wp*f_lnd(i,j)   ! boreal
            lnd%l2d(i,j)%pft_frac        = veg_par%seed_fraction
            lnd%l2d(i,j)%lai_bal         = pft_par%lai_min
            lnd%l2d(i,j)%pft_frac(1:2)  = 0.5_wp
            lnd%l2d(i,j)%lai_bal(1:2)   = pft_par%lai_max(1:2)  
          endif
          ! if fixed LAI, read balanced LAI from file
          if (l_fixlai) then
            call nc_read(trim(lai_fix_file),"lai_bal",lnd%l2d(i,j)%lai_bal(:),start=[1,i,j],count=[npft,1,1])
          endif
          ! if prescribed vegetation read PFT fractions from file
          if (.not.l_dynveg) then
            call nc_read(trim(pft_fix_file),"pft_frac",lnd%l2d(i,j)%pft_frac(:),start=[1,i,j],count=[npft,1,1])
          endif
          !lnd%l2d(i,j)%frac_surf(i_ice)  = f_ice(i,j)-f_ice_grd(i,j)
          call surface_frac_up(f_ice(i,j),f_ice_grd(i,j),f_ocn(i,j),f_lake(i,j),f_veg(i,j),lnd%l2d(i,j)%pft_frac,lnd%l2d(i,j)%frac_surf)

          lnd%l2d(i,j)%albedo = 0.1_wp   ! needed for offline simulations
          lnd%l2d(i,j)%rough_m = 1.e-3_wp
          lnd%l2d(i,j)%rough_h = 1.e-3_wp
          lnd%l2d(i,j)%Ch      = 1.e-3_wp 
          lnd%l2d(i,j)%t_skin = T0
          lnd%l2d(i,j)%et             = 0._wp 
          lnd%l2d(i,j)%flx_sh         = 1._wp 
          lnd%l2d(i,j)%flx_lh         = 1._wp 
          lnd%l2d(i,j)%flx_lwu        = sigma*T0**4 
          lnd%l2d(i,j)%alb_vis_dir    = 0._wp 
          lnd%l2d(i,j)%alb_vis_dif    = 0._wp 
          lnd%l2d(i,j)%alb_nir_dir    = 0._wp 
          lnd%l2d(i,j)%alb_nir_dif    = 0._wp 

          lnd%l2d(i,j)%w_can          = 0._wp 
          lnd%l2d(i,j)%s_can          = 0._wp 
          lnd%l2d(i,j)%mask_snow      = 0._wp  
          lnd%l2d(i,j)%h_snow         = 0._wp  
          lnd%l2d(i,j)%w_snow         = 0._wp  
          lnd%l2d(i,j)%w_snow_max     = 0._wp  
          lnd%l2d(i,j)%npp_ann         = 0._wp 
          lnd%l2d(i,j)%npp13_ann       = 0._wp 
          lnd%l2d(i,j)%npp14_ann       = 0._wp 
          lnd%l2d(i,j)%lai             = lnd%l2d(i,j)%lai_bal
          lnd%l2d(i,j)%sai             = lnd%l2d(i,j)%lai_bal * veg_par%sai_scale 
          lnd%l2d(i,j)%leaf_c          = lnd%l2d(i,j)%lai_bal / pft_par%sla
          lnd%l2d(i,j)%stem_c          = pft_par%awl*lnd%l2d(i,j)%lai_bal**pft_par%bwl
          lnd%l2d(i,j)%root_c          = lnd%l2d(i,j)%leaf_c
          lnd%l2d(i,j)%veg_c           = lnd%l2d(i,j)%leaf_c + lnd%l2d(i,j)%stem_c + lnd%l2d(i,j)%root_c
          lnd%l2d(i,j)%veg_c13         = lnd%l2d(i,j)%veg_c * lnd%l0d%c13_c12_atm 
          lnd%l2d(i,j)%veg_c14         = lnd%l2d(i,j)%veg_c * lnd%l0d%c14_c_atm 
          lnd%l2d(i,j)%phen            = 0._wp 
          lnd%l2d(i,j)%phen_acc        = 0._wp 
          lnd%l2d(i,j)%gdd             = 0._wp 
          lnd%l2d(i,j)%veg_h           = pft_par%awh * lnd%l2d(i,j)%lai_bal
          lnd%l2d(i,j)%seed_frac       = 0._wp 
          lnd%l2d(i,j)%gamma_dist      = 0.001_wp 
          lnd%l2d(i,j)%gamma_luc       = 0._wp 
          lnd%l2d(i,j)%gamma_ice       = 0._wp 
          lnd%l2d(i,j)%root_frac       = pft_par%root_frac
          lnd%l2d(i,j)%litter_in_frac  = pft_par%litter_in_frac
          lnd%l2d(i,j)%t_soil          = T0 
          lnd%l2d(i,j)%t_ice           = T0 
          lnd%l2d(i,j)%t_shelf         = T0 
          lnd%l2d(i,j)%t_lake          = T0 
          lnd%l2d(i,j)%t_sublake       = T0 
          lnd%l2d(i,j)%theta_w(:)      = lnd%l2d(i,j)%theta_sat(:) 
          lnd%l2d(i,j)%theta_i         = 0._wp 
          lnd%l2d(i,j)%w_w             = lnd%l2d(i,j)%theta_w * dz(1:nl) * rho_w 
          lnd%l2d(i,j)%w_i             = 0._wp 
          lnd%l2d(i,j)%theta           = lnd%l2d(i,j)%w_w/(rho_w*dz(1:nl)) 
          lnd%l2d(i,j)%theta_w_shelf   = 0._wp 
          lnd%l2d(i,j)%theta_i_shelf   = 0._wp 
          lnd%l2d(i,j)%w_w_shelf       = 0._wp 
          lnd%l2d(i,j)%w_i_shelf       = 0._wp 
          lnd%l2d(i,j)%theta_w_sublake = 0._wp 
          lnd%l2d(i,j)%theta_i_sublake = 0._wp 
          lnd%l2d(i,j)%w_w_sublake     = 0._wp 
          lnd%l2d(i,j)%w_i_sublake     = 0._wp 
          lnd%l2d(i,j)%f_i_lake        = 0._wp 
          lnd%l2d(i,j)%f_lake_ice      = 0._wp 
          lnd%l2d(i,j)%litter_c        = 0._wp 
          lnd%l2d(i,j)%fast_c          = 0._wp 
          lnd%l2d(i,j)%slow_c          = 0._wp 
          lnd%l2d(i,j)%litter_c13      = 0._wp 
          lnd%l2d(i,j)%fast_c13        = 0._wp 
          lnd%l2d(i,j)%slow_c13        = 0._wp 
          lnd%l2d(i,j)%litter_c14      = 0._wp 
          lnd%l2d(i,j)%fast_c14        = 0._wp 
          lnd%l2d(i,j)%slow_c14        = 0._wp 
          lnd%l2d(i,j)%litter_c_shelf  = 0._wp 
          lnd%l2d(i,j)%fast_c_shelf    = 0._wp 
          lnd%l2d(i,j)%slow_c_shelf    = 0._wp 
          lnd%l2d(i,j)%litter_c13_shelf= 0._wp 
          lnd%l2d(i,j)%fast_c13_shelf  = 0._wp 
          lnd%l2d(i,j)%slow_c13_shelf  = 0._wp 
          lnd%l2d(i,j)%litter_c14_shelf= 0._wp 
          lnd%l2d(i,j)%fast_c14_shelf  = 0._wp 
          lnd%l2d(i,j)%slow_c14_shelf  = 0._wp 
          lnd%l2d(i,j)%litter_c_ice    = 0._wp 
          lnd%l2d(i,j)%fast_c_ice      = 0._wp 
          lnd%l2d(i,j)%slow_c_ice      = 0._wp 
          lnd%l2d(i,j)%litter_c13_ice  = 0._wp 
          lnd%l2d(i,j)%fast_c13_ice    = 0._wp 
          lnd%l2d(i,j)%slow_c13_ice    = 0._wp 
          lnd%l2d(i,j)%litter_c14_ice  = 0._wp 
          lnd%l2d(i,j)%fast_c14_ice    = 0._wp 
          lnd%l2d(i,j)%slow_c14_ice    = 0._wp 
          lnd%l2d(i,j)%litter_c_lake  = 0._wp 
          lnd%l2d(i,j)%fast_c_lake    = 0._wp 
          lnd%l2d(i,j)%slow_c_lake    = 0._wp 
          lnd%l2d(i,j)%litter_c13_lake= 0._wp 
          lnd%l2d(i,j)%fast_c13_lake  = 0._wp 
          lnd%l2d(i,j)%slow_c13_lake  = 0._wp 
          lnd%l2d(i,j)%litter_c14_lake= 0._wp 
          lnd%l2d(i,j)%fast_c14_lake  = 0._wp 
          lnd%l2d(i,j)%slow_c14_lake  = 0._wp 
          lnd%l2d(i,j)%cato_c          = 0._wp 
          lnd%l2d(i,j)%cato_c13        = 0._wp 
          lnd%l2d(i,j)%cato_c14           = 0._wp  
          lnd%l2d(i,j)%f_carb    = 1._wp  
          lnd%l2d(i,j)%weath_carb    = 0._wp  
          lnd%l2d(i,j)%weath13_carb    = 0._wp  
          lnd%l2d(i,j)%weath14_carb    = 0._wp  
          lnd%l2d(i,j)%weath_sil     = 0._wp  
          lnd%l2d(i,j)%weath13_sil     = 0._wp  
          lnd%l2d(i,j)%weath14_sil     = 0._wp  
          lnd%l2d(i,j)%weath_loess    = 0._wp  
          lnd%l2d(i,j)%poc_export     = 0._wp  
          lnd%l2d(i,j)%poc13_export     = 0._wp  
          lnd%l2d(i,j)%poc14_export     = 0._wp  
          lnd%l2d(i,j)%doc_export     = 0._wp  
          lnd%l2d(i,j)%doc13_export     = 0._wp  
          lnd%l2d(i,j)%doc14_export     = 0._wp  
        enddo
      enddo

    endif

    deallocate(f_veg)

    ! initialize
    do i = 1,ni
      do j = 1,nj
        lnd%l2d(i,j)%f_wet_cum      = 0._wp  
        lnd%l2d(i,j)%runoff_ann     = 0._wp  

        lnd%l2d(i,j)%coszm          = 0._wp 
        lnd%l2d(i,j)%daylength      = 0._wp 
        lnd%l2d(i,j)%f_wet_mon      = 0._wp 
        lnd%l2d(i,j)%w_table_mon    = 0._wp 
        lnd%l2d(i,j)%f_wet_long     = 0._wp 
        lnd%l2d(i,j)%tatm           = T0
        lnd%l2d(i,j)%t2m            = T0
        lnd%l2d(i,j)%t_skin_amp      = 0._wp
        lnd%l2d(i,j)%qatm           = 0._wp 
        lnd%l2d(i,j)%q2m            = 0._wp 
        lnd%l2d(i,j)%swnet          = 0._wp 
        lnd%l2d(i,j)%swnet_min      = 0._wp 
        lnd%l2d(i,j)%lwdown         = 0._wp 
        lnd%l2d(i,j)%Ri             = 0._wp 
        lnd%l2d(i,j)%r_a            = 0._wp 
        lnd%l2d(i,j)%r_s            = 0._wp 
        lnd%l2d(i,j)%beta_s         = 0._wp 
        lnd%l2d(i,j)%rain_ground    = 0._wp 
        lnd%l2d(i,j)%evap_can       = 0._wp 
        lnd%l2d(i,j)%snow_ground    = 0._wp 
        lnd%l2d(i,j)%subl_can       = 0._wp 
        lnd%l2d(i,j)%w_can_old      = 0._wp 
        lnd%l2d(i,j)%s_can_old      = 0._wp 
        lnd%l2d(i,j)%f_wat_can      = 0._wp 
        lnd%l2d(i,j)%f_snow_can     = 0._wp 
        lnd%l2d(i,j)%transpiration  = 0._wp 
        lnd%l2d(i,j)%evap_surface   = 0._wp
        lnd%l2d(i,j)%flx_g          = 0._wp 
        lnd%l2d(i,j)%dflxg_dT       = 0._wp 
        lnd%l2d(i,j)%flx_melt       = 0._wp 
        lnd%l2d(i,j)%lwnet          = 0._wp 

        lnd%l2d(i,j)%t_skin_old     = T0
        lnd%l2d(i,j)%num_lh         = 0._wp 
        lnd%l2d(i,j)%num_sh         = 0._wp 
        lnd%l2d(i,j)%num_sw         = 0._wp 
        lnd%l2d(i,j)%num_lw         = 0._wp 
        lnd%l2d(i,j)%denom_lh       = 0._wp 
        lnd%l2d(i,j)%denom_sh       = 0._wp 
        lnd%l2d(i,j)%denom_lw       = 0._wp 
        lnd%l2d(i,j)%f_sh           = 0._wp 
        lnd%l2d(i,j)%f_e            = 0._wp 
        lnd%l2d(i,j)%f_t            = 0._wp 
        lnd%l2d(i,j)%f_le           = 0._wp 
        lnd%l2d(i,j)%f_lt           = 0._wp 
        lnd%l2d(i,j)%f_lw           = 0._wp 
        lnd%l2d(i,j)%lh_ecan        = 0._wp 
        lnd%l2d(i,j)%qsat_e         = 0._wp 
        lnd%l2d(i,j)%dqsatdT_e      = 0._wp 
        lnd%l2d(i,j)%qsat_t         = 0._wp 
        lnd%l2d(i,j)%dqsatdT_t      = 0._wp 
        lnd%l2d(i,j)%energy_cons_surf1 = 0._wp   
        lnd%l2d(i,j)%energy_cons_surf2 = 0._wp
        lnd%l2d(i,j)%alb_snow_vis_dir  = 0._wp
        lnd%l2d(i,j)%alb_snow_vis_dif  = 0._wp
        lnd%l2d(i,j)%alb_snow_nir_dir  = 0._wp
        lnd%l2d(i,j)%alb_snow_nir_dif  = 0._wp
        lnd%l2d(i,j)%w_snow_old     = 0._wp  
        lnd%l2d(i,j)%snowmelt       = 0._wp  
        lnd%l2d(i,j)%icemelt        = 0._wp  

        lnd%l2d(i,j)%snow_grain     = 0._wp  
        lnd%l2d(i,j)%dust_con       = 0._wp  
        lnd%l2d(i,j)%runoff         = 0._wp  
        lnd%l2d(i,j)%runoff_sur     = 0._wp  
        lnd%l2d(i,j)%calving        = 0._wp  
        lnd%l2d(i,j)%drainage       = 0._wp  
        lnd%l2d(i,j)%water_cons     = 0._wp  
        lnd%l2d(i,j)%disturbance    = 0._wp 
        lnd%l2d(i,j)%r_a_can        = 0._wp 
        lnd%l2d(i,j)%r_s_can        = 0._wp 
        lnd%l2d(i,j)%beta_s_can     = 0._wp 
        lnd%l2d(i,j)%ci             = 0._wp 
        lnd%l2d(i,j)%g_can          = 0._wp 
        lnd%l2d(i,j)%gpp            = 0._wp 
        lnd%l2d(i,j)%npp            = 0._wp 
        lnd%l2d(i,j)%npp13          = 0._wp 
        lnd%l2d(i,j)%npp14          = 0._wp 
        lnd%l2d(i,j)%aresp          = 0._wp 
        lnd%l2d(i,j)%discrimination = 0._wp 
        lnd%l2d(i,j)%npp_cum         = 0._wp 
        lnd%l2d(i,j)%npp13_cum       = 0._wp 
        lnd%l2d(i,j)%npp14_cum       = 0._wp 
        lnd%l2d(i,j)%gamma_dist_cum  = 0._wp 
        lnd%l2d(i,j)%wilt      = 0._wp 
        lnd%l2d(i,j)%litterfall    = 0._wp
        lnd%l2d(i,j)%litterfall13  = 0._wp 
        lnd%l2d(i,j)%litterfall14  = 0._wp 
        lnd%l2d(i,j)%lambda_soil     = 0._wp 
        lnd%l2d(i,j)%lambda_int_soil = 0._wp 
        lnd%l2d(i,j)%cap_soil        = 0._wp 
        lnd%l2d(i,j)%lambda_ice      = 0._wp 
        lnd%l2d(i,j)%lambda_int_ice  = 0._wp 
        lnd%l2d(i,j)%cap_ice         = 0._wp 
        lnd%l2d(i,j)%lambda_int_shelf= 0._wp 
        lnd%l2d(i,j)%lambda_lake     = 0._wp 
        lnd%l2d(i,j)%lambda_int_lake = 0._wp 
        lnd%l2d(i,j)%cap_lake        = 0._wp 
        lnd%l2d(i,j)%lambda_int_sublake = 0._wp 
        lnd%l2d(i,j)%cap_sublake        = 0._wp 
        lnd%l2d(i,j)%t_soil_old      = T0
        lnd%l2d(i,j)%t_soil_max      = T0 
        lnd%l2d(i,j)%t_ice_old       = T0 
        lnd%l2d(i,j)%t_shelf_old     = T0 
        lnd%l2d(i,j)%t_shelf_max     = T0 
        lnd%l2d(i,j)%t_lake_old      = T0 
        lnd%l2d(i,j)%veg_c_below     = 0._wp 
        lnd%l2d(i,j)%veg_c13_below   = 0._wp 
        lnd%l2d(i,j)%veg_c14_below   = 0._wp 
        lnd%l2d(i,j)%cap_shelf       = 0._wp 
        lnd%l2d(i,j)%kappa_int       = 0._wp 

        lnd%l2d(i,j)%w_w_old         = 0._wp 
        lnd%l2d(i,j)%w_i_old         = 0._wp 
        lnd%l2d(i,j)%w_w_phase       = 0._wp 
        lnd%l2d(i,j)%w_i_phase       = 0._wp 

        lnd%l2d(i,j)%t_soil_cum      = 0._wp 
        lnd%l2d(i,j)%theta_w_cum     = 0._wp 
        lnd%l2d(i,j)%theta_i_cum     = 0._wp 
        lnd%l2d(i,j)%t_shelf_cum     = 0._wp 
        lnd%l2d(i,j)%theta_w_shelf_cum= 0._wp
        lnd%l2d(i,j)%theta_i_shelf_cum= 0._wp
        lnd%l2d(i,j)%t_sublake_cum= 0._wp 
        lnd%l2d(i,j)%theta_w_sublake_cum= 0._wp
        lnd%l2d(i,j)%theta_i_sublake_cum= 0._wp
        lnd%l2d(i,j)%psi             = 0._wp 
        lnd%l2d(i,j)%ftemp           = 0._wp 
        lnd%l2d(i,j)%fmoist          = 0._wp 
        lnd%l2d(i,j)%fdepth          = 0._wp 
        lnd%l2d(i,j)%k_litter        = 0._wp 
        lnd%l2d(i,j)%k_fast          = 0._wp 
        lnd%l2d(i,j)%k_slow          = 0._wp 
        lnd%l2d(i,j)%k_litter_wet    = 0._wp 
        lnd%l2d(i,j)%k_fast_wet      = 0._wp 
        lnd%l2d(i,j)%k_slow_wet      = 0._wp 
        lnd%l2d(i,j)%diff_soilc      = 0._wp 
        lnd%l2d(i,j)%adv_soilc       = 0._wp 
        lnd%l2d(i,j)%k_litter_shelf  = 0._wp 
        lnd%l2d(i,j)%k_fast_shelf    = 0._wp 
        lnd%l2d(i,j)%k_slow_shelf    = 0._wp 
        lnd%l2d(i,j)%diff_shelfc     = 0._wp 
        lnd%l2d(i,j)%adv_shelfc      = 0._wp 
        lnd%l2d(i,j)%k_litter_ice    = 0._wp 
        lnd%l2d(i,j)%k_fast_ice      = 0._wp 
        lnd%l2d(i,j)%k_slow_ice      = 0._wp 
        lnd%l2d(i,j)%diff_icec       = 0._wp 
        lnd%l2d(i,j)%adv_icec        = 0._wp 
        lnd%l2d(i,j)%k_litter_lake  = 0._wp 
        lnd%l2d(i,j)%k_fast_lake    = 0._wp 
        lnd%l2d(i,j)%k_slow_lake    = 0._wp 
        lnd%l2d(i,j)%diff_lakec     = 0._wp 
        lnd%l2d(i,j)%adv_lakec      = 0._wp 
        lnd%l2d(i,j)%k_cato          = 0._wp 
        lnd%l2d(i,j)%ch4_frac_wet    = 0._wp 
        lnd%l2d(i,j)%ch4_frac_peat   = 0._wp 
        lnd%l2d(i,j)%ch4_frac_shelf  = 0._wp 
        lnd%l2d(i,j)%ch4_frac_lake   = 0._wp 
        lnd%l2d(i,j)%soil_resp_l = 0._wp
        lnd%l2d(i,j)%soil_c_tot         = 0._wp
        lnd%l2d(i,j)%soil_resp          = 0._wp
        lnd%l2d(i,j)%soil_c13_tot       = 0._wp
        lnd%l2d(i,j)%soil_resp13        = 0._wp
        lnd%l2d(i,j)%soil_c14_tot       = 0._wp
        lnd%l2d(i,j)%soil_resp14        = 0._wp
        lnd%l2d(i,j)%frac_soc        = 0._wp 
        lnd%l2d(i,j)%carbon_cons_soil   = 0._wp
        lnd%l2d(i,j)%carbon13_cons_soil = 0._wp
        lnd%l2d(i,j)%carbon14_cons_soil = 0._wp
      enddo
    enddo

    ! seeds everywhere
    if (veg_par%iseed.eq.1) then
      do i = 1,ni
        do j = 1,nj
          if (lnd%l2d(i,j)%mask_lnd.eq.1) then
            lnd%l2d(i,j)%seed_frac = veg_par%seed_fraction
          else
            lnd%l2d(i,j)%seed_frac = 0._wp
          endif
        enddo
      enddo
    endif

    ! roughness length
    do i = 1,ni
      do j = 1,nj
        lnd%l2d(i,j)%z0m(1:npft)  = 0.1_wp * lnd%l2d(i,j)%veg_h(1:npft)  
        lnd%l2d(i,j)%z0m(i_bare)  = surf_par%z0m_bare
        lnd%l2d(i,j)%z0m(i_lake)  = surf_par%z0m_lake
        lnd%l2d(i,j)%z0m(i_ice)   = surf_par%z0m_ice
      enddo
    enddo

    ! initialize mapping between 2-D grid and 1-D array of column models
    lnd%ncells = count(lnd%l2d%mask_lnd==1)    
    allocate(lnd%ij_1d(2,lnd%ncells))
    allocate(lnd%id_map(ni, nj))
    n = 0  ! counter
    do i = 1,ni
      do j = 1,nj
        if (lnd%l2d(i,j)%mask_lnd==1) then
          n = n + 1
          lnd%id_map(i, j) = n
          lnd%ij_1d(1, n) = i
          lnd%ij_1d(2, n) = j
        else
          lnd%id_map(i, j) = 0    ! unassigned points
        endif
      enddo
    enddo

    print*
    print*,'======================================================='
    print*,' Initialisation of LAND complete'
    print*,'======================================================='
    print*

    return

  end subroutine lnd_init 


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l n d _ a l l o c
  !   Purpose    :  land allocate
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lnd_alloc(lnd)

    implicit none 

    type(lnd_class) :: lnd

    integer :: i,j


    do i = 1,ni
      do j = 1,nj
        allocate(lnd%l2d(i,j)%coszm       (nday_year)) 
        allocate(lnd%l2d(i,j)%daylength   (nday_year)) 
        allocate(lnd%l2d(i,j)%f_wet_mon   (nmon_year)) 
        allocate(lnd%l2d(i,j)%w_table_mon (nmon_year)) 
        allocate(lnd%l2d(i,j)%f_wet_long  (nmonwet)) 
        allocate(lnd%l2d(i,j)%frac_surf         (nsurf))
        allocate(lnd%l2d(i,j)%pressure          (nsurf))
        allocate(lnd%l2d(i,j)%tatm              (nsurf))
        allocate(lnd%l2d(i,j)%t2m               (nsurf))
        allocate(lnd%l2d(i,j)%qatm              (nsurf))
        allocate(lnd%l2d(i,j)%q2m               (nsurf))
        allocate(lnd%l2d(i,j)%swnet             (nsurf))
        allocate(lnd%l2d(i,j)%swnet_min         (nsurf))
        allocate(lnd%l2d(i,j)%lwdown            (nsurf))
        allocate(lnd%l2d(i,j)%wind              (nsurf))
        allocate(lnd%l2d(i,j)%rain              (nsurf))
        allocate(lnd%l2d(i,j)%snow              (nsurf))
        allocate(lnd%l2d(i,j)%rough_m           (nsurf))
        allocate(lnd%l2d(i,j)%rough_h           (nsurf))
        allocate(lnd%l2d(i,j)%Ch                (nsurf))
        allocate(lnd%l2d(i,j)%z0m               (nsurf))
        allocate(lnd%l2d(i,j)%Ri                (nsurf))
        allocate(lnd%l2d(i,j)%r_a               (nsurf))
        allocate(lnd%l2d(i,j)%r_s               (nsurf))
        allocate(lnd%l2d(i,j)%beta_s            (nsurf))
        allocate(lnd%l2d(i,j)%albedo            (nsurf))
        allocate(lnd%l2d(i,j)%alb_vis_dir       (nsurf))
        allocate(lnd%l2d(i,j)%alb_vis_dif       (nsurf))
        allocate(lnd%l2d(i,j)%alb_nir_dir       (nsurf))
        allocate(lnd%l2d(i,j)%alb_nir_dif       (nsurf))
        allocate(lnd%l2d(i,j)%rain_ground       (nsurf))
        allocate(lnd%l2d(i,j)%evap_can          (nsurf))
        allocate(lnd%l2d(i,j)%snow_ground       (nsurf))
        allocate(lnd%l2d(i,j)%subl_can          (nsurf))
        allocate(lnd%l2d(i,j)%w_can             (nsurf))
        allocate(lnd%l2d(i,j)%w_can_old         (nsurf))
        allocate(lnd%l2d(i,j)%s_can             (nsurf))
        allocate(lnd%l2d(i,j)%s_can_old         (nsurf))
        allocate(lnd%l2d(i,j)%f_snow            (nsurf))
        allocate(lnd%l2d(i,j)%f_wat_can         (nsurf))
        allocate(lnd%l2d(i,j)%f_snow_can        (nsurf))
        allocate(lnd%l2d(i,j)%transpiration     (nsurf))
        allocate(lnd%l2d(i,j)%evap_surface      (nsurf))
        allocate(lnd%l2d(i,j)%et                (nsurf))
        allocate(lnd%l2d(i,j)%flx_sh            (nsurf))
        allocate(lnd%l2d(i,j)%flx_lh            (nsurf))
        allocate(lnd%l2d(i,j)%flx_g             (nsurf))
        allocate(lnd%l2d(i,j)%dflxg_dT          (nsurf))
        allocate(lnd%l2d(i,j)%flx_melt          (nsurf))
        allocate(lnd%l2d(i,j)%flx_lwu           (nsurf))
        allocate(lnd%l2d(i,j)%lwnet             (nsurf))
        allocate(lnd%l2d(i,j)%t_skin            (nsurf))
        allocate(lnd%l2d(i,j)%t_skin_old        (nsurf))
        allocate(lnd%l2d(i,j)%t_skin_amp        (nsurf))
        allocate(lnd%l2d(i,j)%num_lh            (nsurf))
        allocate(lnd%l2d(i,j)%num_sh            (nsurf))
        allocate(lnd%l2d(i,j)%num_sw            (nsurf))
        allocate(lnd%l2d(i,j)%num_lw            (nsurf))
        allocate(lnd%l2d(i,j)%denom_lh          (nsurf))
        allocate(lnd%l2d(i,j)%denom_sh          (nsurf))
        allocate(lnd%l2d(i,j)%denom_lw          (nsurf))
        allocate(lnd%l2d(i,j)%f_sh              (nsurf))
        allocate(lnd%l2d(i,j)%f_e               (nsurf))
        allocate(lnd%l2d(i,j)%f_t               (nsurf))
        allocate(lnd%l2d(i,j)%f_le              (nsurf))
        allocate(lnd%l2d(i,j)%f_lt              (nsurf))
        allocate(lnd%l2d(i,j)%f_lw              (nsurf))
        allocate(lnd%l2d(i,j)%lh_ecan           (nsurf))
        allocate(lnd%l2d(i,j)%qsat_e            (nsurf))
        allocate(lnd%l2d(i,j)%dqsatdT_e         (nsurf))
        allocate(lnd%l2d(i,j)%qsat_t            (nsurf))
        allocate(lnd%l2d(i,j)%dqsatdT_t         (nsurf))
        allocate(lnd%l2d(i,j)%energy_cons_surf1 (nsurf))
        allocate(lnd%l2d(i,j)%energy_cons_surf2 (nsurf))
        allocate(lnd%l2d(i,j)%alb_snow_vis_dir (nsoil))
        allocate(lnd%l2d(i,j)%alb_snow_vis_dif (nsoil))
        allocate(lnd%l2d(i,j)%alb_snow_nir_dir (nsoil))
        allocate(lnd%l2d(i,j)%alb_snow_nir_dif (nsoil))
        allocate(lnd%l2d(i,j)%mask_snow        (nsoil))
        allocate(lnd%l2d(i,j)%h_snow           (nsoil))
        allocate(lnd%l2d(i,j)%w_snow           (nsoil))
        allocate(lnd%l2d(i,j)%w_snow_max       (nsoil))
        allocate(lnd%l2d(i,j)%w_snow_old       (nsoil))
        allocate(lnd%l2d(i,j)%snowmelt         (nsoil))
        allocate(lnd%l2d(i,j)%icemelt          (nsoil))
        allocate(lnd%l2d(i,j)%snow_grain       (nsoil))
        allocate(lnd%l2d(i,j)%dust_con         (nsoil))
        allocate(lnd%l2d(i,j)%runoff           (nsoil))
        allocate(lnd%l2d(i,j)%runoff_sur       (nsoil))
        allocate(lnd%l2d(i,j)%calving          (nsoil))
        allocate(lnd%l2d(i,j)%drainage         (nsoil))
        allocate(lnd%l2d(i,j)%water_cons       (nsoil))
        allocate(lnd%l2d(i,j)%disturbance     (npft))
        allocate(lnd%l2d(i,j)%r_a_can         (npft))
        allocate(lnd%l2d(i,j)%r_s_can         (npft))
        allocate(lnd%l2d(i,j)%beta_s_can      (npft))
        allocate(lnd%l2d(i,j)%ci              (npft))
        allocate(lnd%l2d(i,j)%g_can           (npft))
        allocate(lnd%l2d(i,j)%gpp             (npft))
        allocate(lnd%l2d(i,j)%npp             (npft))
        allocate(lnd%l2d(i,j)%npp13           (npft))
        allocate(lnd%l2d(i,j)%npp14           (npft))
        allocate(lnd%l2d(i,j)%aresp           (npft))
        allocate(lnd%l2d(i,j)%discrimination  (npft))
        allocate(lnd%l2d(i,j)%lai             (npft))
        allocate(lnd%l2d(i,j)%sai             (npft))
        allocate(lnd%l2d(i,j)%phen            (npft))
        allocate(lnd%l2d(i,j)%phen_acc        (npft))
        allocate(lnd%l2d(i,j)%gdd             (npft))
        allocate(lnd%l2d(i,j)%gamma_leaf      (npft))
        allocate(lnd%l2d(i,j)%lambda          (npft))
        allocate(lnd%l2d(i,j)%lai_bal         (npft))
        allocate(lnd%l2d(i,j)%npp_cum         (npft))
        allocate(lnd%l2d(i,j)%npp13_cum       (npft))
        allocate(lnd%l2d(i,j)%npp14_cum       (npft))
        allocate(lnd%l2d(i,j)%npp_ann         (npft))
        allocate(lnd%l2d(i,j)%npp13_ann       (npft))
        allocate(lnd%l2d(i,j)%npp14_ann       (npft))
        allocate(lnd%l2d(i,j)%veg_c           (npft))
        allocate(lnd%l2d(i,j)%veg_h           (npft))
        allocate(lnd%l2d(i,j)%pft_frac        (npft))
        allocate(lnd%l2d(i,j)%seed_frac       (npft))
        allocate(lnd%l2d(i,j)%gamma_luc       (npft))
        allocate(lnd%l2d(i,j)%gamma_ice       (npft))
        allocate(lnd%l2d(i,j)%gamma_dist      (npft))
        allocate(lnd%l2d(i,j)%gamma_dist_cum  (npft))
        allocate(lnd%l2d(i,j)%leaf_c          (npft))
        allocate(lnd%l2d(i,j)%stem_c          (npft))
        allocate(lnd%l2d(i,j)%root_c          (npft))
        allocate(lnd%l2d(i,j)%veg_c13         (npft))
        allocate(lnd%l2d(i,j)%veg_c14         (npft))
        allocate(lnd%l2d(i,j)%root_frac (nl,npft))  
        allocate(lnd%l2d(i,j)%litter_in_frac (nl))  
        allocate(lnd%l2d(i,j)%wilt      (nl,npft))  
        allocate(lnd%l2d(i,j)%litterfall  (nlc,ncarb)) 
        allocate(lnd%l2d(i,j)%litterfall13(nlc,ncarb)) 
        allocate(lnd%l2d(i,j)%litterfall14(nlc,ncarb)) 
        allocate(lnd%l2d(i,j)%lambda_soil     (0:nl)) 
        allocate(lnd%l2d(i,j)%lambda_int_soil (0:nl)) 
        allocate(lnd%l2d(i,j)%cap_soil        (0:nl)) 
        allocate(lnd%l2d(i,j)%lambda_ice      (0:nl)) 
        allocate(lnd%l2d(i,j)%lambda_int_ice  (0:nl)) 
        allocate(lnd%l2d(i,j)%cap_ice         (0:nl)) 
        allocate(lnd%l2d(i,j)%lambda_int_shelf(0:nl)) 
        allocate(lnd%l2d(i,j)%lambda_lake     (0:nl_l)) 
        allocate(lnd%l2d(i,j)%lambda_int_lake (0:nl_l)) 
        allocate(lnd%l2d(i,j)%cap_lake        (0:nl_l)) 
        allocate(lnd%l2d(i,j)%lambda_int_sublake (0:nl)) 
        allocate(lnd%l2d(i,j)%cap_sublake        (nl)) 
        allocate(lnd%l2d(i,j)%t_soil          (0:nl)) 
        allocate(lnd%l2d(i,j)%t_soil_old      (0:nl)) 
        allocate(lnd%l2d(i,j)%t_soil_max      (0:nl)) 
        allocate(lnd%l2d(i,j)%t_ice           (0:nl)) 
        allocate(lnd%l2d(i,j)%t_ice_old       (0:nl)) 
        allocate(lnd%l2d(i,j)%t_lake          (0:nl_l)) 
        allocate(lnd%l2d(i,j)%t_lake_old      (0:nl_l)) 
        allocate(lnd%l2d(i,j)%t_sublake       (0:nl)) 
        allocate(lnd%l2d(i,j)%t_shelf         (0:nl)) 
        allocate(lnd%l2d(i,j)%t_shelf_old     (0:nl)) 
        allocate(lnd%l2d(i,j)%t_shelf_max     (0:nl)) 
        allocate(lnd%l2d(i,j)%veg_c_below      (nl)) 
        allocate(lnd%l2d(i,j)%veg_c13_below    (nl)) 
        allocate(lnd%l2d(i,j)%veg_c14_below    (nl)) 
        allocate(lnd%l2d(i,j)%cap_shelf        (nl)) 
        allocate(lnd%l2d(i,j)%lambda_s         (nl)) 
        allocate(lnd%l2d(i,j)%lambda_dry       (nl)) 
        allocate(lnd%l2d(i,j)%kappa_int        (nl)) 
        allocate(lnd%l2d(i,j)%theta_sat        (nl)) 
        allocate(lnd%l2d(i,j)%k_sat            (nl)) 
        allocate(lnd%l2d(i,j)%psi_sat          (nl)) 
        allocate(lnd%l2d(i,j)%theta_field      (nl)) 
        allocate(lnd%l2d(i,j)%theta_wilt       (nl)) 
        allocate(lnd%l2d(i,j)%theta_w          (nl)) 
        allocate(lnd%l2d(i,j)%theta_i          (nl)) 
        allocate(lnd%l2d(i,j)%theta            (nl)) 
        allocate(lnd%l2d(i,j)%w_w              (nl)) 
        allocate(lnd%l2d(i,j)%w_i              (nl)) 
        allocate(lnd%l2d(i,j)%w_w_old          (nl)) 
        allocate(lnd%l2d(i,j)%w_i_old          (nl)) 
        allocate(lnd%l2d(i,j)%w_w_phase        (nl)) 
        allocate(lnd%l2d(i,j)%w_i_phase        (nl)) 
        allocate(lnd%l2d(i,j)%theta_w_shelf    (nl)) 
        allocate(lnd%l2d(i,j)%theta_i_shelf    (nl)) 
        allocate(lnd%l2d(i,j)%w_w_shelf        (nl)) 
        allocate(lnd%l2d(i,j)%w_i_shelf        (nl)) 
        allocate(lnd%l2d(i,j)%theta_w_sublake  (nl)) 
        allocate(lnd%l2d(i,j)%theta_i_sublake  (nl)) 
        allocate(lnd%l2d(i,j)%w_w_sublake      (nl)) 
        allocate(lnd%l2d(i,j)%w_i_sublake      (nl)) 
        allocate(lnd%l2d(i,j)%t_soil_cum       (nl)) 
        allocate(lnd%l2d(i,j)%theta_w_cum      (nl)) 
        allocate(lnd%l2d(i,j)%theta_i_cum      (nl)) 
        allocate(lnd%l2d(i,j)%t_shelf_cum      (nl)) 
        allocate(lnd%l2d(i,j)%theta_w_shelf_cum(nl)) 
        allocate(lnd%l2d(i,j)%theta_i_shelf_cum(nl)) 
        allocate(lnd%l2d(i,j)%t_sublake_cum (nl)) 
        allocate(lnd%l2d(i,j)%theta_w_sublake_cum (nl)) 
        allocate(lnd%l2d(i,j)%theta_i_sublake_cum (nl)) 
        allocate(lnd%l2d(i,j)%w_w_lake         (nl_l)) 
        allocate(lnd%l2d(i,j)%w_i_lake         (nl_l)) 
        allocate(lnd%l2d(i,j)%f_i_lake         (nl_l)) 
        allocate(lnd%l2d(i,j)%psi              (nl)) 
        allocate(lnd%l2d(i,j)%k_exp            (nl)) 
        allocate(lnd%l2d(i,j)%psi_exp          (nl)) 
        allocate(lnd%l2d(i,j)%ftemp            (nl)) 
        allocate(lnd%l2d(i,j)%fmoist           (nl)) 
        allocate(lnd%l2d(i,j)%fdepth           (nl)) 
        allocate(lnd%l2d(i,j)%k_litter         (nlc)) 
        allocate(lnd%l2d(i,j)%k_fast           (nlc)) 
        allocate(lnd%l2d(i,j)%k_slow           (nlc)) 
        allocate(lnd%l2d(i,j)%k_litter_wet     (nlc)) 
        allocate(lnd%l2d(i,j)%k_fast_wet       (nlc)) 
        allocate(lnd%l2d(i,j)%k_slow_wet       (nlc)) 
        allocate(lnd%l2d(i,j)%diff_soilc       (nlc)) 
        allocate(lnd%l2d(i,j)%adv_soilc        (nlc)) 
        allocate(lnd%l2d(i,j)%k_litter_shelf   (nlc)) 
        allocate(lnd%l2d(i,j)%k_fast_shelf     (nlc)) 
        allocate(lnd%l2d(i,j)%k_slow_shelf     (nlc)) 
        allocate(lnd%l2d(i,j)%diff_shelfc      (nlc)) 
        allocate(lnd%l2d(i,j)%adv_shelfc       (nlc)) 
        allocate(lnd%l2d(i,j)%k_litter_ice     (nlc)) 
        allocate(lnd%l2d(i,j)%k_fast_ice       (nlc)) 
        allocate(lnd%l2d(i,j)%k_slow_ice       (nlc)) 
        allocate(lnd%l2d(i,j)%diff_icec        (nlc)) 
        allocate(lnd%l2d(i,j)%adv_icec         (nlc)) 
        allocate(lnd%l2d(i,j)%k_litter_lake    (nlc)) 
        allocate(lnd%l2d(i,j)%k_fast_lake      (nlc)) 
        allocate(lnd%l2d(i,j)%k_slow_lake      (nlc)) 
        allocate(lnd%l2d(i,j)%diff_lakec       (nlc)) 
        allocate(lnd%l2d(i,j)%adv_lakec        (nlc)) 
        allocate(lnd%l2d(i,j)%k_cato           (nlc)) 
        allocate(lnd%l2d(i,j)%ch4_frac_wet     (nl)) 
        allocate(lnd%l2d(i,j)%ch4_frac_peat    (nl)) 
        allocate(lnd%l2d(i,j)%ch4_frac_shelf   (nl)) 
        allocate(lnd%l2d(i,j)%ch4_frac_lake    (nl)) 
        allocate(lnd%l2d(i,j)%frac_soc         (nl)) 
        allocate(lnd%l2d(i,j)%litter_c         (nlc)) 
        allocate(lnd%l2d(i,j)%fast_c           (nlc)) 
        allocate(lnd%l2d(i,j)%slow_c           (nlc)) 
        allocate(lnd%l2d(i,j)%litter_c13       (nlc)) 
        allocate(lnd%l2d(i,j)%fast_c13         (nlc)) 
        allocate(lnd%l2d(i,j)%slow_c13         (nlc)) 
        allocate(lnd%l2d(i,j)%litter_c14       (nlc)) 
        allocate(lnd%l2d(i,j)%fast_c14         (nlc)) 
        allocate(lnd%l2d(i,j)%slow_c14         (nlc)) 
        allocate(lnd%l2d(i,j)%litter_c_shelf   (nlc)) 
        allocate(lnd%l2d(i,j)%fast_c_shelf     (nlc)) 
        allocate(lnd%l2d(i,j)%slow_c_shelf     (nlc)) 
        allocate(lnd%l2d(i,j)%litter_c13_shelf (nlc)) 
        allocate(lnd%l2d(i,j)%fast_c13_shelf   (nlc)) 
        allocate(lnd%l2d(i,j)%slow_c13_shelf   (nlc)) 
        allocate(lnd%l2d(i,j)%litter_c14_shelf (nlc)) 
        allocate(lnd%l2d(i,j)%fast_c14_shelf   (nlc)) 
        allocate(lnd%l2d(i,j)%slow_c14_shelf   (nlc)) 
        allocate(lnd%l2d(i,j)%litter_c_ice     (nlc)) 
        allocate(lnd%l2d(i,j)%fast_c_ice       (nlc)) 
        allocate(lnd%l2d(i,j)%slow_c_ice       (nlc)) 
        allocate(lnd%l2d(i,j)%litter_c13_ice   (nlc)) 
        allocate(lnd%l2d(i,j)%fast_c13_ice     (nlc)) 
        allocate(lnd%l2d(i,j)%slow_c13_ice     (nlc)) 
        allocate(lnd%l2d(i,j)%litter_c14_ice   (nlc)) 
        allocate(lnd%l2d(i,j)%fast_c14_ice     (nlc)) 
        allocate(lnd%l2d(i,j)%slow_c14_ice     (nlc)) 
        allocate(lnd%l2d(i,j)%litter_c_lake    (nlc)) 
        allocate(lnd%l2d(i,j)%fast_c_lake      (nlc)) 
        allocate(lnd%l2d(i,j)%slow_c_lake      (nlc)) 
        allocate(lnd%l2d(i,j)%litter_c13_lake  (nlc)) 
        allocate(lnd%l2d(i,j)%fast_c13_lake    (nlc)) 
        allocate(lnd%l2d(i,j)%slow_c13_lake    (nlc)) 
        allocate(lnd%l2d(i,j)%litter_c14_lake  (nlc)) 
        allocate(lnd%l2d(i,j)%fast_c14_lake    (nlc)) 
        allocate(lnd%l2d(i,j)%slow_c14_lake    (nlc)) 
        allocate(lnd%l2d(i,j)%cato_c           (nlc)) 
        allocate(lnd%l2d(i,j)%cato_c13         (nlc)) 
        allocate(lnd%l2d(i,j)%cato_c14         (nlc))     
        allocate(lnd%l2d(i,j)%soil_resp_l(nl,ncarb))
        allocate(lnd%l2d(i,j)%soil_c_tot        (ncarb)) 
        allocate(lnd%l2d(i,j)%soil_resp         (ncarb)) 
        allocate(lnd%l2d(i,j)%soil_c13_tot      (ncarb)) 
        allocate(lnd%l2d(i,j)%soil_resp13       (ncarb)) 
        allocate(lnd%l2d(i,j)%soil_c14_tot      (ncarb)) 
        allocate(lnd%l2d(i,j)%soil_resp14       (ncarb)) 
        allocate(lnd%l2d(i,j)%carbon_cons_soil  (ncarb)) 
        allocate(lnd%l2d(i,j)%carbon13_cons_soil(ncarb)) 
        allocate(lnd%l2d(i,j)%carbon14_cons_soil(ncarb)) 
        allocate(lnd%l2d(i,j)%lithology_gemco2(weath_gemco2_par%nlit)) 
        allocate(lnd%l2d(i,j)%lithology_uhh(weath_uhh_par%nlit)) 
      enddo
    enddo


    return

  end subroutine lnd_alloc


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l n d _ e n d
  !   Purpose    :  end land
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lnd_end(lnd)

    implicit none

    type(lnd_2d_class), allocatable :: lnd(:,:)


    deallocate(lnd)
 

   return

  end subroutine lnd_end


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  w r i t e _ r e s t a r t
  ! Purpose  :  Write restart netcdf file 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lnd_write_restart(fnm,lnd,l0d)

    use ncio
    use dim_name, only: dim_lon, dim_lat, dim_depth, dim_depth0, dim_depth1, dim_depthl, dim_depth0l, dim_npft, dim_nsurf, dim_nsoil

    implicit none

    integer :: ncid
    character (len=*) :: fnm
    type(lnd_2d_class) :: lnd(:,:)
    type(lnd_0d_class) :: l0d

    integer :: i, j
    integer, dimension(:,:,:), allocatable :: vari_n
    real(wp), dimension(:,:,:), allocatable :: var_n

   call nc_create(fnm)
   call nc_write_dim(fnm,"p",x=1)
   call nc_write_dim(fnm,dim_lon,x=lon,axis="x")
   call nc_write_dim(fnm,dim_lat,x=lat,axis="y")
   call nc_write_dim(fnm,dim_depth,x=z(1:nl),units="m",axis="z")
   call nc_write_dim(fnm,dim_depth0,x=z(0:nl),units="m",axis="z")
   call nc_write_dim(fnm,dim_depthl,x=z_l(1:nl_l),units="m",axis="z")
   call nc_write_dim(fnm,dim_depth0l,x=z_l(0:nl_l),units="m",axis="z")
   call nc_write_dim(fnm,dim_depth1,x=z_c(1:nlc),units="m",axis="z")
   call nc_write_dim(fnm,dim_npft,x=i_pft,units="n/a")
   call nc_write_dim(fnm,dim_nsurf,x=i_surf,units="n/a")
   call nc_write_dim(fnm,dim_nsoil,x=i_soil,units="n/a")


   call nc_open(fnm,ncid)

   call nc_write(fnm,"Cflx_avg",  l0d%Cflx_avg,dim1="p",long_name="average land-atm carbon flux",ncid=ncid)
   call nc_write(fnm,"weath_carb_avg",  l0d%weath_carb_avg,dim1="p",long_name="average carbonate weathering flux",ncid=ncid)
   call nc_write(fnm,"weath_sil_avg",   l0d%weath_sil_avg, dim1="p",long_name="average silicate weathering flux",ncid=ncid)

   call nc_write(fnm,"landc",  l0d%landc,dim1="p",long_name="total land carbon",ncid=ncid)
   call nc_write(fnm,"landc13",l0d%landc13,dim1="p",long_name="total land carbon 13",ncid=ncid)
   call nc_write(fnm,"landc14",l0d%landc14,dim1="p",long_name="total land carbon 14",ncid=ncid)

   call nc_write(fnm,"burc",  l0d%burc,dim1="p",long_name="total buried carbon",ncid=ncid)
   call nc_write(fnm,"burc13",l0d%burc13,dim1="p",long_name="total buried carbon 13",ncid=ncid)
   call nc_write(fnm,"burc14",l0d%burc14,dim1="p",long_name="total buried carbon 14",ncid=ncid)

   call nc_write(fnm,"weath_scale",  l0d%weath_scale,dim1="p",long_name="weathering scaling factor to match alkalinity input into the ocean",ncid=ncid)

   call nc_write(fnm,"mask_lnd",     lnd%mask_lnd,    dims=[dim_lon,dim_lat],long_name="land mask",units="/",ncid=ncid)

   call nc_write(fnm,"f_land0",    lnd%f_land0,   dims=[dim_lon,dim_lat],long_name="land fraction at present day sea level",units="/",ncid=ncid)
   call nc_write(fnm,"f_land",     lnd%f_land,    dims=[dim_lon,dim_lat],long_name="land fraction",units="/",ncid=ncid)
   call nc_write(fnm,"f_ice",      lnd%f_ice,     dims=[dim_lon,dim_lat],long_name="ice fraction",units="/",ncid=ncid)
   call nc_write(fnm,"f_ice_grd",      lnd%f_ice_grd,     dims=[dim_lon,dim_lat],long_name="grounded ice fraction",units="/",ncid=ncid)
   call nc_write(fnm,"f_shelf",    lnd%f_shelf,   dims=[dim_lon,dim_lat],long_name="ocean shelf fraction",units="/",ncid=ncid)
   call nc_write(fnm,"f_lake",     lnd%f_lake,    dims=[dim_lon,dim_lat],long_name="lake fraction",units="/",ncid=ncid)
   call nc_write(fnm,"f_veg",      lnd%f_veg,     dims=[dim_lon,dim_lat],long_name="vegetation fraction",units="/",ncid=ncid)

   call nc_write(fnm,"alt",        lnd%alt,    dims=[dim_lon,dim_lat],long_name="active layer thickness",units="m",ncid=ncid)
   call nc_write(fnm,"gdd5",       lnd%gdd5,           dims=[dim_lon,dim_lat],long_name="growing degree days above 5 degC",units="K",ncid=ncid)
   call nc_write(fnm,"t2m_min_mon",lnd%t2m_min_mon,  dims=[dim_lon,dim_lat],long_name="coldest month temperature",units="K",ncid=ncid)
   call nc_write(fnm,"f_lake_ice", lnd%f_lake_ice,     dims=[dim_lon,dim_lat],long_name="lake ice fraction",units="/",ncid=ncid)

   call nc_write(fnm,"t_skin_veg",      lnd%t_skin_veg,     dims=[dim_lon,dim_lat],long_name="mean cell skin temperature",units="K",ncid=ncid)

   allocate(var_n(npft,nx,ny))
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%gdd(:)
     enddo
   enddo
   call nc_write(fnm,"gdd",              var_n,    dims=[dim_npft,dim_lon,dim_lat],start=[1,1,1],count=[npft,ni,nj],long_name="growing degree days",units="K",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%phen(:)
     enddo
   enddo
   call nc_write(fnm,"phen",             var_n,      dims=[dim_npft,dim_lon,dim_lat],start=[1,1,1],count=[npft,ni,nj],long_name="phenology",units="/",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%phen_acc(:)
     enddo
   enddo
   call nc_write(fnm,"phen_acc",         var_n,  dims=[dim_npft,dim_lon,dim_lat],start=[1,1,1],count=[npft,ni,nj],long_name="phenology accumulated",units="/",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%lai_bal(:)
     enddo
   enddo
   call nc_write(fnm,"lai_bal",          var_n,   dims=[dim_npft,dim_lon,dim_lat],start=[1,1,1],count=[npft,ni,nj],long_name="balanced leaf area index",units="m2/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%lai(:)
     enddo
   enddo
   call nc_write(fnm,"lai",              var_n,       dims=[dim_npft,dim_lon,dim_lat],start=[1,1,1],count=[npft,ni,nj],long_name="leaf area index",units="m2/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%sai(:)
     enddo
   enddo
   call nc_write(fnm,"sai",              var_n,       dims=[dim_npft,dim_lon,dim_lat],start=[1,1,1],count=[npft,ni,nj],long_name="stem area index",units="m2/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%seed_frac(:)
     enddo
   enddo
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%npp_ann(:)
     enddo
   enddo
   call nc_write(fnm,"npp_ann",              var_n,       dims=[dim_npft,dim_lon,dim_lat],start=[1,1,1],count=[npft,ni,nj],long_name="leaf area index",units="kgC/m2/s",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%npp13_ann(:)
     enddo
   enddo
   call nc_write(fnm,"npp13_ann",              var_n,       dims=[dim_npft,dim_lon,dim_lat],start=[1,1,1],count=[npft,ni,nj],long_name="leaf area index",units="kgC/m2/s",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%npp14_ann(:)
     enddo
   enddo
   call nc_write(fnm,"npp14_ann",              var_n,       dims=[dim_npft,dim_lon,dim_lat],start=[1,1,1],count=[npft,ni,nj],long_name="leaf area index",units="kgC/m2/s",ncid=ncid)
   call nc_write(fnm,"seed_frac",        var_n, dims=[dim_npft,dim_lon,dim_lat],start=[1,1,1],count=[npft,ni,nj],long_name="seed fraction",units="/",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%pft_frac(:)
     enddo
   enddo
   call nc_write(fnm,"pft_frac",         var_n,  dims=[dim_npft,dim_lon,dim_lat],start=[1,1,1],count=[npft,ni,nj],long_name="vegetation fraction",units="/",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%veg_h(:)
     enddo
   enddo
   call nc_write(fnm,"veg_h",            var_n,     dims=[dim_npft,dim_lon,dim_lat],start=[1,1,1],count=[npft,ni,nj],long_name="vegetation height",units="m",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%veg_c(:)
     enddo
   enddo
   call nc_write(fnm,"veg_c",            var_n,     dims=[dim_npft,dim_lon,dim_lat],start=[1,1,1],count=[npft,ni,nj],long_name="vegetation carbon",units="kgC/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%veg_c13(:)
     enddo
   enddo
   call nc_write(fnm,"veg_c13",          var_n,   dims=[dim_npft,dim_lon,dim_lat],start=[1,1,1],count=[npft,ni,nj],long_name="vegetation carbon 13",units="kgC/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%veg_c14(:)
     enddo
   enddo
   call nc_write(fnm,"veg_c14",          var_n,   dims=[dim_npft,dim_lon,dim_lat],start=[1,1,1],count=[npft,ni,nj],long_name="vegetation carbon 14",units="kgC/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%leaf_c(:)
     enddo
   enddo
   call nc_write(fnm,"leaf_c",           var_n,    dims=[dim_npft,dim_lon,dim_lat],start=[1,1,1],count=[npft,ni,nj],long_name="leaf carbon",units="kgC/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%root_c(:)
     enddo
   enddo
   call nc_write(fnm,"root_c",           var_n,    dims=[dim_npft,dim_lon,dim_lat],start=[1,1,1],count=[npft,ni,nj],long_name="root carbon",units="kgC/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%stem_c(:)
     enddo
   enddo
   call nc_write(fnm,"stem_c",           var_n,    dims=[dim_npft,dim_lon,dim_lat],start=[1,1,1],count=[npft,ni,nj],long_name="stem carbon",units="kgC/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%gamma_dist(:)
     enddo
   enddo
   call nc_write(fnm,"gamma_dist",       var_n,dims=[dim_npft,dim_lon,dim_lat],start=[1,1,1],count=[npft,ni,nj],long_name="vegetation disturbance rate",units="1/s",ncid=ncid)
   deallocate(var_n)

   do i=1,nx
     do j=1,ny
       call nc_write(fnm,"root_frac",    lnd(i,j)%root_frac, dims=[dim_depth,dim_npft,dim_lon,dim_lat],start=[1,1,i,j],count=[nl,npft,1,1],long_name="root fraction in layers",units="/",ncid=ncid)
     enddo
   enddo

   allocate(var_n(nsurf,nx,ny))
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%frac_surf(:)
     enddo
   enddo
   call nc_write(fnm,"frac_surf",     var_n, dims=[dim_nsurf,dim_lon,dim_lat],start=[1,1,1],count=[nsurf,ni,nj],long_name="surface types fraction",units="/",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%rough_m(:)
     enddo
   enddo
   call nc_write(fnm,"rough_m",       var_n, dims=[dim_nsurf,dim_lon,dim_lat],start=[1,1,1],count=[nsurf,ni,nj],long_name="surface roughness for momentum",units="m",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%rough_h(:)
     enddo
   enddo
   call nc_write(fnm,"rough_h",       var_n, dims=[dim_nsurf,dim_lon,dim_lat],start=[1,1,1],count=[nsurf,ni,nj],long_name="surface roughness for heat",units="m",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%Ch(:)
     enddo
   enddo
   call nc_write(fnm,"Cd",       var_n, dims=[dim_nsurf,dim_lon,dim_lat],start=[1,1,1],count=[nsurf,ni,nj],long_name="exchange coefficient for heat",units="1",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%albedo(:)
     enddo
   enddo
   call nc_write(fnm,"albedo",        var_n,    dims=[dim_nsurf,dim_lon,dim_lat],start=[1,1,1],count=[nsurf,ni,nj],long_name="surface albedo",units="/",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%w_can(:)
     enddo
   enddo
   call nc_write(fnm,"w_can",         var_n,     dims=[dim_nsurf,dim_lon,dim_lat],start=[1,1,1],count=[nsurf,ni,nj],long_name="canopy water",units="kg/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%s_can(:)
     enddo
   enddo
   call nc_write(fnm,"s_can",         var_n,     dims=[dim_nsurf,dim_lon,dim_lat],start=[1,1,1],count=[nsurf,ni,nj],long_name="canopy snow",units="kg/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%t_skin(:)
     enddo
   enddo
   call nc_write(fnm,"t_skin",        var_n,    dims=[dim_nsurf,dim_lon,dim_lat],start=[1,1,1],count=[nsurf,ni,nj],long_name="skin temperature",units="K",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%flx_sh(:)
     enddo
   enddo
   call nc_write(fnm,"flx_sh",        var_n,    dims=[dim_nsurf,dim_lon,dim_lat],start=[1,1,1],count=[nsurf,ni,nj],long_name="sensible heat flux",units="W/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%flx_lh(:)
     enddo
   enddo
   call nc_write(fnm,"flx_lh",        var_n,    dims=[dim_nsurf,dim_lon,dim_lat],start=[1,1,1],count=[nsurf,ni,nj],long_name="latent heat flux",units="W/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%flx_lwu(:)
     enddo
   enddo
   call nc_write(fnm,"flx_lwu",       var_n,    dims=[dim_nsurf,dim_lon,dim_lat],start=[1,1,1],count=[nsurf,ni,nj],long_name="longwave up at the surface",units="W/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%et(:)
     enddo
   enddo
   call nc_write(fnm,"et",            var_n,    dims=[dim_nsurf,dim_lon,dim_lat],start=[1,1,1],count=[nsurf,ni,nj],long_name="evapotranspiration",units="kg/m2/s",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%alb_vis_dir(:)
     enddo
   enddo
   call nc_write(fnm,"alb_vis_dir",   var_n,    dims=[dim_nsurf,dim_lon,dim_lat],start=[1,1,1],count=[nsurf,ni,nj],long_name="alb_vis_dir",units="1",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%alb_vis_dif(:)
     enddo
   enddo
   call nc_write(fnm,"alb_vis_dif",   var_n,    dims=[dim_nsurf,dim_lon,dim_lat],start=[1,1,1],count=[nsurf,ni,nj],long_name="alb_vis_dif",units="1",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%alb_nir_dir(:)
     enddo
   enddo
   call nc_write(fnm,"alb_nir_dir",   var_n,    dims=[dim_nsurf,dim_lon,dim_lat],start=[1,1,1],count=[nsurf,ni,nj],long_name="alb_nir_dir",units="1",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%alb_nir_dif(:)
     enddo
   enddo
   call nc_write(fnm,"alb_nir_dif",   var_n,    dims=[dim_nsurf,dim_lon,dim_lat],start=[1,1,1],count=[nsurf,ni,nj],long_name="alb_nir_dif",units="1",ncid=ncid)
   deallocate(var_n)

   allocate(var_n(nsoil,nx,ny))
   allocate(vari_n(nsoil,nx,ny))
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%w_snow(:)
     enddo
   enddo
   call nc_write(fnm,"w_snow",          var_n,          dims=[dim_nsoil,dim_lon,dim_lat],start=[1,1,1],count=[nsoil,ni,nj],long_name="snow water equivalent",units="kg/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%w_snow_max(:)
     enddo
   enddo
   call nc_write(fnm,"w_snow_max",      var_n, dims=[dim_nsoil,dim_lon,dim_lat],start=[1,1,1],count=[nsoil,ni,nj],long_name="maximum snow water equivalent",units="kg/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%h_snow(:)
     enddo
   enddo
   call nc_write(fnm,"h_snow",          var_n,          dims=[dim_nsoil,dim_lon,dim_lat],start=[1,1,1],count=[nsoil,ni,nj],long_name="snow thickness",units="m",ncid=ncid)
   do i=1,nx
     do j=1,ny
       vari_n(:,i,j) = lnd(i,j)%mask_snow(:)
     enddo
   enddo
   call nc_write(fnm,"mask_snow",       vari_n,       dims=[dim_nsoil,dim_lon,dim_lat],start=[1,1,1],count=[nsoil,ni,nj],long_name="snow mask",units="/",ncid=ncid)
   deallocate(var_n)
   deallocate(vari_n)

   allocate(var_n(0:nl,nx,ny))
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%t_soil(:)
     enddo
   enddo
   call nc_write(fnm,"t_soil",           var_n,          dims=[dim_depth0,dim_lon,dim_lat],start=[1,1,1],count=[nl+1,ni,nj],long_name="soil temperature",units="K",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%t_ice(:)
     enddo
   enddo
   call nc_write(fnm,"t_ice",            var_n,           dims=[dim_depth0,dim_lon,dim_lat],start=[1,1,1],count=[nl+1,ni,nj],long_name="ice temperature",units="K",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%t_shelf(:)
     enddo
   enddo
   call nc_write(fnm,"t_shelf",          var_n,         dims=[dim_depth0,dim_lon,dim_lat],start=[1,1,1],count=[nl+1,ni,nj],long_name="shelf soil temperature",units="K",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%t_sublake(:)
     enddo
   enddo
   call nc_write(fnm,"t_sublake",     var_n,    dims=[dim_depth0,dim_lon,dim_lat],start=[1,1,1],count=[nl+1,ni,nj],long_name="soil temperature below_lake",units="K",ncid=ncid)
   deallocate(var_n)

   allocate(var_n(0:nl_l,nx,ny))
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%t_lake(:)
     enddo
   enddo
   call nc_write(fnm,"t_lake",     var_n,    dims=[dim_depth0l,dim_lon,dim_lat],start=[1,1,1],count=[nl_l+1,ni,nj],long_name="lake temperature",units="K",ncid=ncid)
   deallocate(var_n)

   allocate(var_n(nl,nx,ny))
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%theta(:)
     enddo
   enddo
   call nc_write(fnm,"theta",            var_n,                  dims=[dim_depth,dim_lon,dim_lat],start=[1,1,1],count=[nl,ni,nj],long_name="soil total volumetric water content",units="m3/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%theta_i(:)
     enddo
   enddo
   call nc_write(fnm,"theta_i",          var_n,         dims=[dim_depth,dim_lon,dim_lat],start=[1,1,1],count=[nl,ni,nj],long_name="soil frozen volumetric water content",units="m3/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%theta_w(:)
     enddo
   enddo
   call nc_write(fnm,"theta_w",          var_n,         dims=[dim_depth,dim_lon,dim_lat],start=[1,1,1],count=[nl,ni,nj],long_name="soil liquid volumetric water content",units="m3/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%w_w(:)
     enddo
   enddo
   call nc_write(fnm,"w_w",              var_n,             dims=[dim_depth,dim_lon,dim_lat],start=[1,1,1],count=[nl,ni,nj],long_name="soil liquid water",units="kg/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%w_i(:)
     enddo
   enddo
   call nc_write(fnm,"w_i",              var_n,             dims=[dim_depth,dim_lon,dim_lat],start=[1,1,1],count=[nl,ni,nj],long_name="soil frozen water equivalent",units="kg/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%theta_i_shelf(:)
     enddo
   enddo
   call nc_write(fnm,"theta_i_shelf",    var_n,   dims=[dim_depth,dim_lon,dim_lat],start=[1,1,1],count=[nl,ni,nj],long_name="shelf frozen volumetric water content",units="m3/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%theta_w_shelf(:)
     enddo
   enddo
   call nc_write(fnm,"theta_w_shelf",    var_n,   dims=[dim_depth,dim_lon,dim_lat],start=[1,1,1],count=[nl,ni,nj],long_name="shelf liquid volumetric water content",units="m3/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%w_w_shelf(:)
     enddo
   enddo
   call nc_write(fnm,"w_w_shelf",        var_n,       dims=[dim_depth,dim_lon,dim_lat],start=[1,1,1],count=[nl,ni,nj],long_name="shelf liquid water",units="kg/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%w_i_shelf(:)
     enddo
   enddo
   call nc_write(fnm,"w_i_shelf",        var_n,       dims=[dim_depth,dim_lon,dim_lat],start=[1,1,1],count=[nl,ni,nj],long_name="shelf frozen water equivalent",units="kg/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%theta_i_sublake(:)
     enddo
   enddo
   call nc_write(fnm,"theta_i_sublake",    var_n,     dims=[dim_depth,dim_lon,dim_lat],start=[1,1,1],count=[nl,ni,nj],long_name="sublake frozen volumetric water content",units="m3/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%theta_w_sublake(:)
     enddo
   enddo
   call nc_write(fnm,"theta_w_sublake",    var_n,     dims=[dim_depth,dim_lon,dim_lat],start=[1,1,1],count=[nl,ni,nj],long_name="sublake liquid volumetric water content",units="m3/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%w_w_sublake(:)
     enddo
   enddo
   call nc_write(fnm,"w_w_sublake",        var_n,         dims=[dim_depth,dim_lon,dim_lat],start=[1,1,1],count=[nl,ni,nj],long_name="sublake liquid water",units="kg/m2",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%w_i_sublake(:)
     enddo
   enddo
   call nc_write(fnm,"w_i_sublake",        var_n,       dims=[dim_depth,dim_lon,dim_lat],start=[1,1,1],count=[nl,ni,nj],long_name="sublake frozen water equivalent",units="kg/m2",ncid=ncid)
   deallocate(var_n)

   allocate(var_n(nl_l,nx,ny))
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%f_i_lake(:)
     enddo
   enddo
   call nc_write(fnm,"f_i_lake",        var_n,       dims=[dim_depthl,dim_lon,dim_lat],start=[1,1,1],count=[nl_l,ni,nj],long_name="lake frozen water fraction",units="/",ncid=ncid)
   deallocate(var_n)

   allocate(var_n(nlc,nx,ny))
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%litter_c(:)
     enddo
   enddo
   call nc_write(fnm,"litter_c",         var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="litter carbon",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%fast_c(:)
     enddo
   enddo
   call nc_write(fnm,"fast_c",           var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="fast soil carbon",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%slow_c(:)
     enddo
   enddo
   call nc_write(fnm,"slow_c",           var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="slow soil carbon",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%litter_c13(:)
     enddo
   enddo
   call nc_write(fnm,"litter_c13",       var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="litter carbon 13",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%fast_c13(:)
     enddo
   enddo
   call nc_write(fnm,"fast_c13",         var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="fast soil carbon 13",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%slow_c13(:)
     enddo
   enddo
   call nc_write(fnm,"slow_c13",         var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="slow soil carbon 13",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%litter_c14(:)
     enddo
   enddo
   call nc_write(fnm,"litter_c14",       var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="litter carbon 14",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%fast_c14(:)
     enddo
   enddo
   call nc_write(fnm,"fast_c14",         var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="fast soil carbon 14",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%slow_c14(:)
     enddo
   enddo
   call nc_write(fnm,"slow_c14",         var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="slow soil carbon 14",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%cato_c(:)
     enddo
   enddo
   call nc_write(fnm,"cato_c",           var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="catotelm carbon",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%cato_c13(:)
     enddo
   enddo
   call nc_write(fnm,"cato_c13",         var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="catotelm carbon 13",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%cato_c14(:)
     enddo
   enddo
   call nc_write(fnm,"cato_c14",         var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="catotelm carbon 14",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%litter_c_shelf(:)
     enddo
   enddo
   call nc_write(fnm,"litter_c_shelf",   var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="litter shelf carbon",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%fast_c_shelf(:)
     enddo
   enddo
   call nc_write(fnm,"fast_c_shelf",     var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="fast shelf carbon",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%slow_c_shelf(:)
     enddo
   enddo
   call nc_write(fnm,"slow_c_shelf",     var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="slow sheöf carbon",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%litter_c13_shelf(:)
     enddo
   enddo
   call nc_write(fnm,"litter_c13_shelf", var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="litter shelf carbon 13",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%fast_c13_shelf(:)
     enddo
   enddo
   call nc_write(fnm,"fast_c13_shelf",   var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="fast shelf carbon 13",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%slow_c13_shelf(:)
     enddo
   enddo
   call nc_write(fnm,"slow_c13_shelf",   var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="slow shelf carbon 13",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%litter_c14_shelf(:)
     enddo
   enddo
   call nc_write(fnm,"litter_c14_shelf", var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="litter shelf carbon 14",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%fast_c14_shelf(:)
     enddo
   enddo
   call nc_write(fnm,"fast_c14_shelf",   var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="fast shelf carbon 14",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%slow_c14_shelf(:)
     enddo
   enddo
   call nc_write(fnm,"slow_c14_shelf",   var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="slow shelf carbon 14",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%litter_c_ice(:)
     enddo
   enddo
   call nc_write(fnm,"litter_c_ice",     var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="litter ice carbon",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%fast_c_ice(:)
     enddo
   enddo
   call nc_write(fnm,"fast_c_ice",       var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="fast ice carbon",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%slow_c_ice(:)
     enddo
   enddo
   call nc_write(fnm,"slow_c_ice",       var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="slow ice carbon",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%litter_c13_ice(:)
     enddo
   enddo
   call nc_write(fnm,"litter_c13_ice",   var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="litter ice carbon 13",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%fast_c13_ice(:)
     enddo
   enddo
   call nc_write(fnm,"fast_c13_ice",     var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="fast ice carbon 13",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%slow_c13_ice(:)
     enddo
   enddo
   call nc_write(fnm,"slow_c13_ice",     var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="slow ice carbon 13",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%litter_c14_ice(:)
     enddo
   enddo
   call nc_write(fnm,"litter_c14_ice",   var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="litter ice carbon 14",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%fast_c14_ice(:)
     enddo
   enddo
   call nc_write(fnm,"fast_c14_ice",     var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="fast ice carbon 14",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%slow_c14_ice(:)
     enddo
   enddo
   call nc_write(fnm,"slow_c14_ice",     var_n,      dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="slow ice carbon 14",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%litter_c_lake(:)
     enddo
   enddo
   call nc_write(fnm,"litter_c_lake",   var_n,       dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="litter lake carbon",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%fast_c_lake(:)
     enddo
   enddo
   call nc_write(fnm,"fast_c_lake",     var_n,       dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="fast lake carbon",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%slow_c_lake(:)
     enddo
   enddo
   call nc_write(fnm,"slow_c_lake",     var_n,       dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="slow sheöf carbon",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%litter_c13_lake(:)
     enddo
   enddo
   call nc_write(fnm,"litter_c13_lake", var_n,       dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="litter lake carbon 13",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%fast_c13_lake(:)
     enddo
   enddo
   call nc_write(fnm,"fast_c13_lake",   var_n,       dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="fast lake carbon 13",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%slow_c13_lake(:)
     enddo
   enddo
   call nc_write(fnm,"slow_c13_lake",   var_n,       dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="slow lake carbon 13",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%litter_c14_lake(:)
     enddo
   enddo
   call nc_write(fnm,"litter_c14_lake", var_n,       dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="litter lake carbon 14",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%fast_c14_lake(:)
     enddo
   enddo
   call nc_write(fnm,"fast_c14_lake",   var_n,       dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="fast lake carbon 14",units="kgC/m3",ncid=ncid)
   do i=1,nx
     do j=1,ny
       var_n(:,i,j) = lnd(i,j)%slow_c14_lake(:)
     enddo
   enddo
   call nc_write(fnm,"slow_c14_lake",   var_n,       dims=[dim_depth1,dim_lon,dim_lat],start=[1,1,1],count=[nlc,ni,nj],long_name="slow lake carbon 14",units="kgC/m3",ncid=ncid)
   deallocate(var_n)


   call nc_write(fnm,"litter_c_peat",     lnd%litter_c_peat,    dims=[dim_lon,dim_lat],long_name="peatland litter carbon",units="kgC/m2",ncid=ncid)
   call nc_write(fnm,"acro_c",           lnd%acro_c,          dims=[dim_lon,dim_lat],long_name="acrotelm carbon",units="kgC/m2",ncid=ncid)
   call nc_write(fnm,"litter_c13_peat",   lnd%litter_c13_peat,  dims=[dim_lon,dim_lat],long_name="peatland litter carbon 13",units="kgC/m2",ncid=ncid)
   call nc_write(fnm,"acro_c13",         lnd%acro_c13,        dims=[dim_lon,dim_lat],long_name="acrotelm carbon 13",units="kgC/m2",ncid=ncid)
   call nc_write(fnm,"litter_c14_peat",   lnd%litter_c14_peat,  dims=[dim_lon,dim_lat],long_name="peatland litter carbon 14",units="kgC/m2",ncid=ncid)
   call nc_write(fnm,"acro_c14",         lnd%acro_c14,        dims=[dim_lon,dim_lat],long_name="acrotelm carbon 14",units="kgC/m2",ncid=ncid)
   call nc_write(fnm,"f_peat",           lnd%f_peat,          dims=[dim_lon,dim_lat],long_name="peatland fraction",units="/",ncid=ncid)
   call nc_write(fnm,"f_peat_pot",       lnd%f_peat_pot,      dims=[dim_lon,dim_lat],long_name="potential peatland fraction",units="/",ncid=ncid)
   call nc_write(fnm,"dCpeat_dt",       lnd%dCpeat_dt,        dims=[dim_lon,dim_lat],long_name="peat accumulation rate",units="kgC/m2/s",ncid=ncid)
   call nc_write(fnm,"w_table_min",      lnd%w_table_min,     dims=[dim_lon,dim_lat],long_name="minumum yearly water table depth",units="m",ncid=ncid)
   call nc_write(fnm,"w_table_peat",     lnd%w_table_peat,    dims=[dim_lon,dim_lat],long_name="peatland water table depth",units="m",ncid=ncid)

   call nc_write(fnm,"f_carb",  lnd%f_carb, dims=[dim_lon,dim_lat],long_name="carbonate weathering rate",units="mol C/m2/yr",ncid=ncid)
   call nc_write(fnm,"weath_carb",  lnd%weath_carb, dims=[dim_lon,dim_lat],long_name="carbonate weathering rate",units="mol C/m2/yr",ncid=ncid)
   call nc_write(fnm,"weath13_carb",  lnd%weath13_carb, dims=[dim_lon,dim_lat],long_name="carbonate 13 weathering rate",units="mol C/m2/yr",ncid=ncid)
   call nc_write(fnm,"weath14_carb",  lnd%weath14_carb, dims=[dim_lon,dim_lat],long_name="carbonate 14 weathering rate",units="mol C/m2/yr",ncid=ncid)
   call nc_write(fnm,"weath_sil",   lnd%weath_sil,  dims=[dim_lon,dim_lat],long_name="silicate weathering rate",units="mol C/m2/yr",ncid=ncid)
   call nc_write(fnm,"weath13_sil",   lnd%weath13_sil,  dims=[dim_lon,dim_lat],long_name="silicate 13 weathering rate",units="mol C/m2/yr",ncid=ncid)
   call nc_write(fnm,"weath14_sil",   lnd%weath14_sil,  dims=[dim_lon,dim_lat],long_name="silicate 14 weathering rate",units="mol C/m2/yr",ncid=ncid)

   call nc_write(fnm,"poc_export",     lnd%poc_export,    dims=[dim_lon,dim_lat],long_name="POC export rate",units="kgC/m2/yr",ncid=ncid)
   call nc_write(fnm,"poc13_export",   lnd%poc13_export,  dims=[dim_lon,dim_lat],long_name="POC 13 export rate",units="kgC/m2/yr",ncid=ncid)
   call nc_write(fnm,"poc14_export",   lnd%poc14_export,  dims=[dim_lon,dim_lat],long_name="POC 14 export rate",units="kgC/m2/yr",ncid=ncid)
   call nc_write(fnm,"doc_export",     lnd%doc_export,    dims=[dim_lon,dim_lat],long_name="DOC export rate",units="kgC/m2/yr",ncid=ncid)
   call nc_write(fnm,"doc13_export",   lnd%doc13_export,  dims=[dim_lon,dim_lat],long_name="DOC 13 export rate",units="kgC/m2/yr",ncid=ncid)
   call nc_write(fnm,"doc14_export",   lnd%doc14_export,  dims=[dim_lon,dim_lat],long_name="DOC 14 export rate",units="kgC/m2/yr",ncid=ncid)

   call nc_close(ncid)

   return

  end subroutine lnd_write_restart


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  r e a d _ r e s t a r t
  ! Purpose  :  Read restart netcdf file 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lnd_read_restart(fnm,lnd,l0d)

    use ncio

    implicit none

    character (len=*) :: fnm
    type(lnd_2d_class) :: lnd(:,:)
    type(lnd_0d_class) :: l0d

    integer :: i, j, ncid


    call nc_open(fnm,ncid)
    call nc_read(fnm,"Cflx_avg",l0d%Cflx_avg,ncid=ncid)
    call nc_read(fnm,"weath_carb_avg",l0d%weath_carb_avg,ncid=ncid)
    call nc_read(fnm,"weath_sil_avg",l0d%weath_sil_avg,ncid=ncid)

    call nc_read(fnm,"landc",l0d%landc,ncid=ncid)
    call nc_read(fnm,"landc13",l0d%landc13,ncid=ncid)
    call nc_read(fnm,"landc14",l0d%landc14,ncid=ncid)

    call nc_read(fnm,"burc",  l0d%burc,ncid=ncid)
    call nc_read(fnm,"burc13",l0d%burc13,ncid=ncid)
    call nc_read(fnm,"burc14",l0d%burc14,ncid=ncid)

    call nc_read(fnm,"weath_scale",  l0d%weath_scale,ncid=ncid)

    call nc_read(fnm,"mask_lnd",      lnd%mask_lnd,ncid=ncid)

    call nc_read(fnm,"f_land0",     lnd%f_land0,ncid=ncid)
    call nc_read(fnm,"f_land",      lnd%f_land,ncid=ncid)
    call nc_read(fnm,"f_ice",      lnd%f_ice,ncid=ncid)
    call nc_read(fnm,"f_ice_grd",      lnd%f_ice_grd,ncid=ncid)
    call nc_read(fnm,"f_shelf",     lnd%f_shelf,ncid=ncid)
    call nc_read(fnm,"f_lake",     lnd%f_lake,ncid=ncid)
    call nc_read(fnm,"f_veg",      lnd%f_veg,ncid=ncid)

    call nc_read(fnm,"alt",       lnd%alt,ncid=ncid)
    call nc_read(fnm,"gdd5",             lnd%gdd5,ncid=ncid)
    call nc_read(fnm,"t2m_min_mon",    lnd%t2m_min_mon,ncid=ncid)
    call nc_read(fnm,"f_lake_ice",    lnd%f_lake_ice,ncid=ncid)

    call nc_read(fnm,"t_skin_veg",      lnd%t_skin_veg,ncid=ncid)

    do i=1,nx
      do j=1,ny
        call nc_read(fnm,"gdd",lnd(i,j)%gdd,start=[1,i,j],count=[npft,1,1],ncid=ncid)
        call nc_read(fnm,"phen",lnd(i,j)%phen,start=[1,i,j],count=[npft,1,1],ncid=ncid)
        call nc_read(fnm,"phen_acc",lnd(i,j)%phen_acc,start=[1,i,j],count=[npft,1,1],ncid=ncid)
        call nc_read(fnm,"lai_bal",lnd(i,j)%lai_bal,start=[1,i,j],count=[npft,1,1],ncid=ncid)
        call nc_read(fnm,"lai",lnd(i,j)%lai,start=[1,i,j],count=[npft,1,1],ncid=ncid)
        call nc_read(fnm,"sai",lnd(i,j)%sai,start=[1,i,j],count=[npft,1,1],ncid=ncid)
        call nc_read(fnm,"npp_ann",lnd(i,j)%npp_ann,start=[1,i,j],count=[npft,1,1],ncid=ncid)
        call nc_read(fnm,"npp13_ann",lnd(i,j)%npp13_ann,start=[1,i,j],count=[npft,1,1],ncid=ncid)
        call nc_read(fnm,"npp14_ann",lnd(i,j)%npp14_ann,start=[1,i,j],count=[npft,1,1],ncid=ncid)
        call nc_read(fnm,"seed_frac",lnd(i,j)%seed_frac,start=[1,i,j],count=[npft,1,1],ncid=ncid)
        call nc_read(fnm,"pft_frac",lnd(i,j)%pft_frac,start=[1,i,j],count=[npft,1,1],ncid=ncid)
        call nc_read(fnm,"veg_h",lnd(i,j)%veg_h,start=[1,i,j],count=[npft,1,1],ncid=ncid)
        call nc_read(fnm,"veg_c",lnd(i,j)%veg_c,start=[1,i,j],count=[npft,1,1],ncid=ncid)
        call nc_read(fnm,"veg_c13",lnd(i,j)%veg_c13,start=[1,i,j],count=[npft,1,1],ncid=ncid)
        call nc_read(fnm,"veg_c14",lnd(i,j)%veg_c14,start=[1,i,j],count=[npft,1,1],ncid=ncid)
        call nc_read(fnm,"leaf_c",lnd(i,j)%leaf_c,start=[1,i,j],count=[npft,1,1],ncid=ncid)
        call nc_read(fnm,"root_c",lnd(i,j)%root_c,start=[1,i,j],count=[npft,1,1],ncid=ncid)
        call nc_read(fnm,"stem_c",lnd(i,j)%stem_c,start=[1,i,j],count=[npft,1,1],ncid=ncid)
        call nc_read(fnm,"gamma_dist",lnd(i,j)%gamma_dist,start=[1,i,j],count=[npft,1,1],ncid=ncid)
        call nc_read(fnm,"root_frac",lnd(i,j)%root_frac,start=[1,1,i,j],count=[nl,npft,1,1],ncid=ncid)

        call nc_read(fnm,"frac_surf",lnd(i,j)%frac_surf,start=[1,i,j],count=[nsurf,1,1],ncid=ncid)
        call nc_read(fnm,"rough_m",lnd(i,j)%rough_m,start=[1,i,j],count=[nsurf,1,1],ncid=ncid)
        call nc_read(fnm,"rough_h",lnd(i,j)%rough_h,start=[1,i,j],count=[nsurf,1,1],ncid=ncid)
        call nc_read(fnm,"Cd",lnd(i,j)%Ch,start=[1,i,j],count=[nsurf,1,1],ncid=ncid)
        call nc_read(fnm,"albedo",lnd(i,j)%albedo,start=[1,i,j],count=[nsurf,1,1],ncid=ncid)
        call nc_read(fnm,"w_can",lnd(i,j)%w_can,start=[1,i,j],count=[nsurf,1,1],ncid=ncid)
        call nc_read(fnm,"s_can",lnd(i,j)%s_can,start=[1,i,j],count=[nsurf,1,1],ncid=ncid)
        call nc_read(fnm,"t_skin",lnd(i,j)%t_skin,start=[1,i,j],count=[nsurf,1,1],ncid=ncid)
        call nc_read(fnm,"flx_sh",lnd(i,j)%flx_sh,start=[1,i,j],count=[nsurf,1,1],ncid=ncid)
        call nc_read(fnm,"flx_lh",lnd(i,j)%flx_lh,start=[1,i,j],count=[nsurf,1,1],ncid=ncid)
        call nc_read(fnm,"flx_lwu",lnd(i,j)%flx_lwu,start=[1,i,j],count=[nsurf,1,1],ncid=ncid)
        call nc_read(fnm,"et",lnd(i,j)%et,start=[1,i,j],count=[nsurf,1,1],ncid=ncid)
        call nc_read(fnm,"alb_vis_dir",lnd(i,j)%alb_vis_dir,start=[1,i,j],count=[nsurf,1,1],ncid=ncid)
        call nc_read(fnm,"alb_vis_dif",lnd(i,j)%alb_vis_dif,start=[1,i,j],count=[nsurf,1,1],ncid=ncid)
        call nc_read(fnm,"alb_nir_dir",lnd(i,j)%alb_nir_dir,start=[1,i,j],count=[nsurf,1,1],ncid=ncid)
        call nc_read(fnm,"alb_nir_dif",lnd(i,j)%alb_nir_dif,start=[1,i,j],count=[nsurf,1,1],ncid=ncid)

        call nc_read(fnm,"w_snow",lnd(i,j)%w_snow,start=[1,i,j],count=[nsoil,1,1],ncid=ncid)
        call nc_read(fnm,"w_snow_max",lnd(i,j)%w_snow_max,start=[1,i,j],count=[nsoil,1,1],ncid=ncid)
        call nc_read(fnm,"h_snow",lnd(i,j)%h_snow,start=[1,i,j],count=[nsoil,1,1],ncid=ncid)
        call nc_read(fnm,"mask_snow",lnd(i,j)%mask_snow,start=[1,i,j],count=[nsoil,1,1],ncid=ncid)

        call nc_read(fnm,"t_soil",lnd(i,j)%t_soil,start=[1,i,j],count=[nl+1,1,1],ncid=ncid)
        call nc_read(fnm,"t_ice",lnd(i,j)%t_ice,start=[1,i,j],count=[nl+1,1,1],ncid=ncid)
        call nc_read(fnm,"t_shelf",lnd(i,j)%t_shelf,start=[1,i,j],count=[nl+1,1,1],ncid=ncid)
        call nc_read(fnm,"t_sublake",lnd(i,j)%t_sublake,start=[1,i,j],count=[nl+1,1,1],ncid=ncid)
        call nc_read(fnm,"t_lake",lnd(i,j)%t_lake,start=[1,i,j],count=[nl_l+1,1,1],ncid=ncid)

        call nc_read(fnm,"theta",lnd(i,j)%theta,start=[1,i,j],count=[nl,1,1],ncid=ncid)
        call nc_read(fnm,"theta_i",lnd(i,j)%theta_i,start=[1,i,j],count=[nl,1,1],ncid=ncid)
        call nc_read(fnm,"theta_w",lnd(i,j)%theta_w,start=[1,i,j],count=[nl,1,1],ncid=ncid)
        call nc_read(fnm,"w_w",lnd(i,j)%w_w,start=[1,i,j],count=[nl,1,1],ncid=ncid)
        call nc_read(fnm,"w_i",lnd(i,j)%w_i,start=[1,i,j],count=[nl,1,1],ncid=ncid)
        call nc_read(fnm,"theta_i_shelf",lnd(i,j)%theta_i_shelf,start=[1,i,j],count=[nl,1,1],ncid=ncid)
        call nc_read(fnm,"theta_w_shelf",lnd(i,j)%theta_w_shelf,start=[1,i,j],count=[nl,1,1],ncid=ncid)
        call nc_read(fnm,"w_w_shelf",lnd(i,j)%w_w_shelf,start=[1,i,j],count=[nl,1,1],ncid=ncid)
        call nc_read(fnm,"w_i_shelf",lnd(i,j)%w_i_shelf,start=[1,i,j],count=[nl,1,1],ncid=ncid)
        call nc_read(fnm,"theta_i_sublake",lnd(i,j)%theta_i_sublake,start=[1,i,j],count=[nl,1,1],ncid=ncid)
        call nc_read(fnm,"theta_w_sublake",lnd(i,j)%theta_w_sublake,start=[1,i,j],count=[nl,1,1],ncid=ncid)
        call nc_read(fnm,"w_w_sublake",lnd(i,j)%w_w_sublake,start=[1,i,j],count=[nl,1,1],ncid=ncid)
        call nc_read(fnm,"w_i_sublake",lnd(i,j)%w_i_sublake,start=[1,i,j],count=[nl,1,1],ncid=ncid)
        call nc_read(fnm,"f_i_lake",lnd(i,j)%f_i_lake,start=[1,i,j],count=[nl_l,1,1],ncid=ncid)

        call nc_read(fnm,"litter_c",lnd(i,j)%litter_c,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"fast_c",lnd(i,j)%fast_c,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"slow_c",lnd(i,j)%slow_c,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"litter_c13",lnd(i,j)%litter_c13,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"fast_c13",lnd(i,j)%fast_c13,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"slow_c13",lnd(i,j)%slow_c13,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"litter_c14",lnd(i,j)%litter_c14,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"fast_c14",lnd(i,j)%fast_c14,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"slow_c14",lnd(i,j)%slow_c14,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"litter_c_shelf",lnd(i,j)%litter_c_shelf,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"fast_c_shelf",lnd(i,j)%fast_c_shelf,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"slow_c_shelf",lnd(i,j)%slow_c_shelf,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"litter_c13_shelf",lnd(i,j)%litter_c13_shelf,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"fast_c13_shelf",lnd(i,j)%fast_c13_shelf,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"slow_c13_shelf",lnd(i,j)%slow_c13_shelf,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"litter_c14_shelf",lnd(i,j)%litter_c14_shelf,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"fast_c14_shelf",lnd(i,j)%fast_c14_shelf,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"slow_c14_shelf",lnd(i,j)%slow_c14_shelf,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"litter_c_ice",lnd(i,j)%litter_c_ice,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"fast_c_ice",lnd(i,j)%fast_c_ice,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"slow_c_ice",lnd(i,j)%slow_c_ice,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"litter_c13_ice",lnd(i,j)%litter_c13_ice,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"fast_c13_ice",lnd(i,j)%fast_c13_ice,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"slow_c13_ice",lnd(i,j)%slow_c13_ice,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"litter_c14_ice",lnd(i,j)%litter_c14_ice,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"fast_c14_ice",lnd(i,j)%fast_c14_ice,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"slow_c14_ice",lnd(i,j)%slow_c14_ice,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"litter_c_lake",lnd(i,j)%litter_c_lake,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"fast_c_lake",lnd(i,j)%fast_c_lake,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"slow_c_lake",lnd(i,j)%slow_c_lake,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"litter_c13_lake",lnd(i,j)%litter_c13_lake,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"fast_c13_lake",lnd(i,j)%fast_c13_lake,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"slow_c13_lake",lnd(i,j)%slow_c13_lake,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"litter_c14_lake",lnd(i,j)%litter_c14_lake,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"fast_c14_lake",lnd(i,j)%fast_c14_lake,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"slow_c14_lake",lnd(i,j)%slow_c14_lake,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"cato_c",lnd(i,j)%cato_c,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"cato_c13",lnd(i,j)%cato_c13,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
        call nc_read(fnm,"cato_c14",lnd(i,j)%cato_c14,start=[1,i,j],count=[nlc,1,1],ncid=ncid)
      enddo
    enddo

    call nc_read(fnm,"litter_c_peat",    lnd%litter_c_peat,ncid=ncid)
    call nc_read(fnm,"acro_c",           lnd%acro_c,ncid=ncid)
    call nc_read(fnm,"litter_c13_peat",  lnd%litter_c13_peat,ncid=ncid)
    call nc_read(fnm,"acro_c13",         lnd%acro_c13,ncid=ncid)
    call nc_read(fnm,"litter_c14_peat",  lnd%litter_c14_peat,ncid=ncid)
    call nc_read(fnm,"acro_c14",         lnd%acro_c14,ncid=ncid)
    call nc_read(fnm,"f_peat",           lnd%f_peat,ncid=ncid)
    call nc_read(fnm,"f_peat_pot",       lnd%f_peat_pot,ncid=ncid)
    call nc_read(fnm,"dCpeat_dt",        lnd%dCpeat_dt,ncid=ncid)
    call nc_read(fnm,"w_table_min",      lnd%w_table_min,ncid=ncid)
    call nc_read(fnm,"w_table_peat",     lnd%w_table_peat,ncid=ncid)

    call nc_read(fnm,"f_carb",  lnd%f_carb,ncid=ncid)
    call nc_read(fnm,"weath_carb",  lnd%weath_carb,ncid=ncid)
    call nc_read(fnm,"weath13_carb",  lnd%weath13_carb,ncid=ncid)
    call nc_read(fnm,"weath14_carb",  lnd%weath14_carb,ncid=ncid)
    call nc_read(fnm,"weath_sil",   lnd%weath_sil,ncid=ncid)
    call nc_read(fnm,"weath13_sil",   lnd%weath13_sil,ncid=ncid)
    call nc_read(fnm,"weath14_sil",   lnd%weath14_sil,ncid=ncid)

    call nc_read(fnm,"poc_export",     lnd%poc_export,ncid=ncid)
    call nc_read(fnm,"poc13_export",   lnd%poc13_export,ncid=ncid)
    call nc_read(fnm,"poc14_export",   lnd%poc14_export,ncid=ncid)
    call nc_read(fnm,"doc_export",     lnd%doc_export,ncid=ncid)
    call nc_read(fnm,"doc13_export",   lnd%doc13_export,ncid=ncid)
    call nc_read(fnm,"doc14_export",   lnd%doc14_export,ncid=ncid)

    call nc_close(ncid)

    print *,'read restart file ',fnm

   return

  end subroutine lnd_read_restart

end module lnd_model 
