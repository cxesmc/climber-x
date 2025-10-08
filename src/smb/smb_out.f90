!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s m b _ o u t
!
!  Purpose : SMB model diagnostics and output
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2024 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit
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
module smb_out

  use precision, only : wp
  use constants, only : rho_i, frac_vu
  use dim_name, only: dim_x, dim_y, dim_depth, dim_time, dim_month, dim_day
  use timer, only : doy, nyears, year, year_smb, year_now, sec_day, sec_mon, sec_year, &
  mon, nmon_year, nday_year 
  use timer, only : n_year_smb, nstep_year_smb, nstep_mon_smb, nyout_smb, &
  ny_out_ts, y_out_ts_smb, time_out_ts_smb
  use timer, only : time_soy_smb, time_eoy_smb, time_eom_smb, time_out_smb, year_out_start
  use control, only : out_dir
  use smb_def, only : smb_class, ts_out, s_out 
  use smb_grid, only : z, nl
  use smb_params, only : dt, i_smb, l_daily_output, l_monthly_output
  use ncio
  use coord, only : grid_class

  implicit none

  private
  public :: smb_diag, smb_diag_init


contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s m b _ d i a g _ i n i t
  ! Purpose  :  Initialize netcdf output for smb
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_diag_init(smb)

    implicit none

    type(smb_class) :: smb

    integer :: k
    character (len=256) :: fnm

    smb%nout = 0

    ! allocate
    allocate(smb%ann_ts(ny_out_ts))

    ! initialize netcdf output
    fnm = trim(out_dir)//"/smb_"//trim(smb%grid%name)//"_ts.nc"
    call ts_nc(fnm)

    fnm = trim(out_dir)//"/smb_"//trim(smb%grid%name)//".nc"
    call smb_nc(fnm,smb%grid)

    if (l_daily_output) then
      fnm = trim(out_dir)//"/smb_"//trim(smb%grid%name)//"_daily.nc"
      call smb_daily_nc(fnm,smb%grid)
    endif

    allocate(smb%ann_s%z_sur_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%t2m_i(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%t2m_bias_i(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%tam_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%gam_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%ram_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%prc_i(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%prc_bias_i(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%u700_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%v700_i(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%wind_i(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%cld_i(smb%grid%G%nx,smb%grid%G%ny))    
    allocate(smb%ann_s%dust_i(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%swdown_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%swd_toa_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%swd_sur_vis_dir_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%swd_sur_nir_dir_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%swd_sur_vis_dif_i(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%swd_sur_nir_dif_i(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%dswd_dalb_vis_dir_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%dswd_dalb_nir_dir_i(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%ann_s%dswd_dalb_vis_dif_i(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%ann_s%dswd_dalb_nir_dif_i(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%ann_s%dswd_dz_nir_dir_i(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%ann_s%dswd_dz_nir_dif_i(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%ann_s%coszm_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%lwdown_i(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%gam_lw_i(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%t_ground_i(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%area(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%dz_dx_sur(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%dz_dy_sur(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%dz_sur(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%mask_smb(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%mask_ice(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%mask_maxice(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%mask_margin(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%z_sur(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%z_sur_eff(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%z_sur_fil(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%z_sur_std(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%tam(smb%grid%G%nx,smb%grid%G%ny))    
    allocate(smb%ann_s%t2m(smb%grid%G%nx,smb%grid%G%ny))    
    allocate(smb%ann_s%t_skin(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%tskin_tam(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%t_skin_amp(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%t_prof(smb%grid%G%nx,smb%grid%G%ny,0:nl)) 
    allocate(smb%ann_s%q2m(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%pressure(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%mask_snow(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%ann_s%f_snow(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%h_snow(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%w_snow(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%ann_s%w_snow_max(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%snowmelt(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%icemelt(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%refreezing(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%f_rfz_to_snow(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%runoff(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%albedo(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%f_ice(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%h_ice(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%alb_ice(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%alb_bg(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%alb_vis_dir(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%alb_nir_dir(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%alb_vis_dif(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%alb_nir_dif(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%alb_snow_vis_dir(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%alb_snow_nir_dir(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%alb_snow_vis_dif(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%alb_snow_nir_dif(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%snow_grain(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%dust_con(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%cld(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%swnet(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%swnet_min(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%swdown(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%lwdown(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%ann_s%r_a(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%flx_g(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%ann_s%dflxg_dT(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%flx_melt(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%flx_sh(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%flx_lwu(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%flx_lh(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%evp(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%prc(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%f_ele(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%f_wind(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%rain(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%ann_s%snow(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%u700(smb%grid%G%nx,smb%grid%G%ny))    
    allocate(smb%ann_s%v700(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%wind(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%ann_s%t_ice(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%t_ground(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%pdd(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%ann_s%ann_smb_pdd(smb%grid%G%nx,smb%grid%G%ny))        
    allocate(smb%ann_s%ann_smb(smb%grid%G%nx,smb%grid%G%ny))        
    allocate(smb%ann_s%ann_smb_ice(smb%grid%G%nx,smb%grid%G%ny))        
    allocate(smb%ann_s%ann_smb_noice(smb%grid%G%nx,smb%grid%G%ny))        
    allocate(smb%ann_s%ann_prc(smb%grid%G%nx,smb%grid%G%ny))        
    allocate(smb%ann_s%ann_prc_ice(smb%grid%G%nx,smb%grid%G%ny))        
    allocate(smb%ann_s%ann_snow(smb%grid%G%nx,smb%grid%G%ny))       
    allocate(smb%ann_s%ann_snow_ice(smb%grid%G%nx,smb%grid%G%ny))       
    allocate(smb%ann_s%ann_ablation(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%ann_ablation_ice(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%ann_s%ann_melt(smb%grid%G%nx,smb%grid%G%ny))       
    allocate(smb%ann_s%ann_melt_ice(smb%grid%G%nx,smb%grid%G%ny))       
    allocate(smb%ann_s%ann_evp(smb%grid%G%nx,smb%grid%G%ny))        
    allocate(smb%ann_s%ann_evp_ice(smb%grid%G%nx,smb%grid%G%ny))        
    allocate(smb%ann_s%ann_runoff(smb%grid%G%nx,smb%grid%G%ny))     
    allocate(smb%ann_s%ann_runoff_ice(smb%grid%G%nx,smb%grid%G%ny))     
    allocate(smb%ann_s%ann_refreezing(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%ann_s%ann_refreezing_ice(smb%grid%G%nx,smb%grid%G%ny)) 

    do k=1,nmon_year
    allocate(smb%mon_s(k)%t2m_i(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%mon_s(k)%t2m_bias_i(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%mon_s(k)%tam_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%gam_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%ram_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%prc_i(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%prc_bias_i(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%mon_s(k)%u700_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%v700_i(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%wind_i(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%mon_s(k)%cld_i(smb%grid%G%nx,smb%grid%G%ny))    
    allocate(smb%mon_s(k)%dust_i(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%mon_s(k)%swdown_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%swd_toa_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%swd_sur_vis_dir_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%swd_sur_nir_dir_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%swd_sur_vis_dif_i(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%swd_sur_nir_dif_i(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%dswd_dalb_vis_dir_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%dswd_dalb_nir_dir_i(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%mon_s(k)%dswd_dalb_vis_dif_i(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%mon_s(k)%dswd_dalb_nir_dif_i(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%mon_s(k)%dswd_dz_nir_dir_i(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%mon_s(k)%dswd_dz_nir_dif_i(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%mon_s(k)%coszm_i(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%lwdown_i(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%mon_s(k)%gam_lw_i(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%mon_s(k)%tam(smb%grid%G%nx,smb%grid%G%ny))    
    allocate(smb%mon_s(k)%t2m(smb%grid%G%nx,smb%grid%G%ny))    
    allocate(smb%mon_s(k)%t_skin(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%tskin_tam(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%t_skin_amp(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%t_prof(smb%grid%G%nx,smb%grid%G%ny,0:nl)) 
    allocate(smb%mon_s(k)%q2m(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%pressure(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%mask_snow(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%mon_s(k)%f_snow(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%h_snow(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%w_snow(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%mon_s(k)%w_snow_max(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%snowmelt(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%icemelt(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%runoff(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%refreezing(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%f_rfz_to_snow(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%albedo(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%alb_vis_dir(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%alb_nir_dir(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%alb_vis_dif(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%alb_nir_dif(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%alb_snow_vis_dir(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%alb_snow_nir_dir(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%alb_snow_vis_dif(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%alb_snow_nir_dif(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%snow_grain(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%dust_con(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%cld(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%swnet(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%swnet_min(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%swdown(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%lwdown(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%mon_s(k)%r_a(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%flx_g(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%mon_s(k)%dflxg_dT(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%flx_melt(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%flx_sh(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%flx_lwu(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%flx_lh(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%evp(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%mon_s(k)%prc(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%f_wind(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%mon_s(k)%rain(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%mon_s(k)%snow(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%mon_s(k)%u700(smb%grid%G%nx,smb%grid%G%ny))    
    allocate(smb%mon_s(k)%v700(smb%grid%G%nx,smb%grid%G%ny))   
    allocate(smb%mon_s(k)%wind(smb%grid%G%nx,smb%grid%G%ny))
    enddo

    do k=1,nday_year
    allocate(smb%day_s(k)%t2m(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%day_s(k)%t_skin(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%day_s(k)%t_skin_amp(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%day_s(k)%t_prof(smb%grid%G%nx,smb%grid%G%ny,0:nl)) 
    allocate(smb%day_s(k)%mask_snow(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%day_s(k)%f_snow(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%day_s(k)%h_snow(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%day_s(k)%w_snow(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%day_s(k)%w_snow_max(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%day_s(k)%snowmelt(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%day_s(k)%icemelt(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%day_s(k)%runoff(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%day_s(k)%refreezing(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%day_s(k)%flx_g(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%day_s(k)%dflxg_dT(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%day_s(k)%flx_melt(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%day_s(k)%flx_sh(smb%grid%G%nx,smb%grid%G%ny))  
    allocate(smb%day_s(k)%flx_lwu(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%day_s(k)%flx_lh(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%day_s(k)%rain(smb%grid%G%nx,smb%grid%G%ny))
    allocate(smb%day_s(k)%snow(smb%grid%G%nx,smb%grid%G%ny)) 
    allocate(smb%day_s(k)%dust_con(smb%grid%G%nx,smb%grid%G%ny)) 
    enddo


   return

  end subroutine smb_diag_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s m b _ d i a g
  !   Purpose    :  sea ice diagnostics
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_diag(smb)

    implicit none

    type(smb_class) :: smb

    integer :: m, k, y
    integer :: ppos
    character (len=256) :: fnm
    character (len=256) :: dom
    real(wp) :: mon_avg, ann_avg
    real(wp) :: area_ice


    ! current index
    y = y_out_ts_smb

    mon_avg = 1._wp/nstep_mon_smb
    ann_avg = 1._wp/nstep_year_smb

    if( time_eoy_smb ) then

     ! integrated surface mass balance of ice sheet
     smb%ann_ts(y)%smb = sum(smb%ann_smb*smb%grid%area,mask=smb%mask_ice==1) * 1.e-6_wp ! Gt/yr 
     ! integrated precipitation on ice sheet
     smb%ann_ts(y)%prc = sum(smb%ann_prc*smb%grid%area,mask=smb%mask_ice==1) * 1.e-6_wp ! Gt/yr 
     ! integrated snowfall on ice sheet
     smb%ann_ts(y)%snow= sum(smb%ann_snow*smb%grid%area,mask=smb%mask_ice==1) * 1.e-6_wp ! Gt/yr 
     ! integrated runoff from ice sheet
     smb%ann_ts(y)%run = sum(smb%ann_runoff*smb%grid%area,mask=smb%mask_ice==1) * 1.e-6_wp ! Gt/yr 
     ! integrated melt of ice sheet
     smb%ann_ts(y)%melt = sum(smb%ann_melt*smb%grid%area,mask=smb%mask_ice==1) * 1.e-6_wp ! Gt/yr 
     ! integrated refreezing of ice sheet
     smb%ann_ts(y)%refreezing= sum(smb%ann_refreezing*smb%grid%area,mask=smb%mask_ice==1) * 1.e-6_wp ! Gt/yr 
     ! integrated sublimation of ice sheet
     smb%ann_ts(y)%evp = sum(smb%ann_evp*smb%grid%area,mask=smb%mask_ice==1) * 1.e-6_wp ! Gt/yr 

     area_ice = sum(smb%grid%area,mask=smb%mask_ice==1)
     ! ice area
     smb%ann_ts(y)%Aice = area_ice*1.e-6_wp  ! 10^6 km^2
     if (area_ice.gt.0._wp) then
       ! averaged surface mass balance of ice sheet
       smb%ann_ts(y)%smb_avg = sum(smb%ann_smb*smb%grid%area,mask=smb%mask_ice==1) / area_ice  ! kg/m2/yr 
       ! averaged precipitation on ice sheet
       smb%ann_ts(y)%prc_avg = sum(smb%ann_prc*smb%grid%area,mask=smb%mask_ice==1) / area_ice  ! kg/m2/yr 
       ! averaged snowfall on ice sheet
       smb%ann_ts(y)%snow_avg= sum(smb%ann_snow*smb%grid%area,mask=smb%mask_ice==1) / area_ice  ! kg/m2/yr 
       ! averaged runoff from ice sheet
       smb%ann_ts(y)%run_avg = sum(smb%ann_runoff*smb%grid%area,mask=smb%mask_ice==1) / area_ice  ! kg/m2/yr 
       ! averaged melt of ice sheet
       smb%ann_ts(y)%melt_avg = sum(smb%ann_melt*smb%grid%area,mask=smb%mask_ice==1) / area_ice  ! kg/m2/yr 
       ! averaged refreezing of ice sheet
       smb%ann_ts(y)%refreezing_avg= sum(smb%ann_refreezing*smb%grid%area,mask=smb%mask_ice==1) / area_ice  ! kg/m2/yr 
       ! averaged sublimation of ice sheet
       smb%ann_ts(y)%evp_avg = sum(smb%ann_evp*smb%grid%area,mask=smb%mask_ice==1) / area_ice  ! kg/m2/yr 
     else
       ! averaged surface mass balance of ice sheet
       smb%ann_ts(y)%smb_avg = 0._wp 
       ! averaged precipitation on ice sheet
       smb%ann_ts(y)%prc_avg = 0._wp 
       ! averaged snowfall on ice sheet
       smb%ann_ts(y)%snow_avg= 0._wp 
       ! averaged runoff from ice sheet
       smb%ann_ts(y)%run_avg = 0._wp 
       ! averaged melt of ice sheet
       smb%ann_ts(y)%melt_avg = 0._wp 
       ! averaged refreezing of ice sheet
       smb%ann_ts(y)%refreezing_avg= 0._wp 
       ! averaged sublimation of ice sheet
       smb%ann_ts(y)%evp_avg = 0._wp 
     endif

     ! write to netcdf file 
     if (time_out_ts_smb) then
       fnm = trim(out_dir)//"/smb_"//trim(smb%grid%name)//"_ts.nc"
       call ts_nc_write(fnm,smb%ann_ts(1:y),year_smb-y+1,y)
     endif
     ppos = scan(trim(smb%grid%name),"-")-1
     dom = trim(smb%grid%name(1:ppos)) 
     ! print header
     if (mod(year,10).eq.1) then
       print '(a7,a9,7a7)','smb_'//dom,'year','SMB','prc','snow','melt','runoff','refrz','sub'
     endif

     ! print values
     print '(a7,i9,7F7.0)', &
       'smb_'//dom,year_now,smb%ann_ts(y)%smb,smb%ann_ts(y)%prc,smb%ann_ts(y)%snow,smb%ann_ts(y)%melt, &
     smb%ann_ts(y)%run,smb%ann_ts(y)%refreezing,smb%ann_ts(y)%evp

    endif


    ! spatially explicit output
    if ( time_out_smb ) then

      if( time_soy_smb ) then
        do m=1,nmon_year
          smb%mon_s(m)%t2m_i            = 0._wp 
          smb%mon_s(m)%t2m_bias_i       = 0._wp 
          smb%mon_s(m)%tam_i            = 0._wp 
          smb%mon_s(m)%gam_i            = 0._wp 
          smb%mon_s(m)%ram_i            = 0._wp 
          smb%mon_s(m)%prc_i            = 0._wp 
          smb%mon_s(m)%prc_bias_i       = 0._wp 
          smb%mon_s(m)%u700_i           = 0._wp 
          smb%mon_s(m)%v700_i           = 0._wp 
          smb%mon_s(m)%wind_i           = 0._wp 
          smb%mon_s(m)%cld_i            = 0._wp 
          smb%mon_s(m)%dust_i           = 0._wp 
          smb%mon_s(m)%swdown_i         = 0._wp 
          smb%mon_s(m)%swd_toa_i        = 0._wp 
          smb%mon_s(m)%swd_sur_vis_dir_i = 0._wp 
          smb%mon_s(m)%swd_sur_nir_dir_i = 0._wp 
          smb%mon_s(m)%swd_sur_vis_dif_i = 0._wp 
          smb%mon_s(m)%swd_sur_nir_dif_i = 0._wp 
          smb%mon_s(m)%dswd_dalb_vis_dir_i = 0._wp
          smb%mon_s(m)%dswd_dalb_nir_dir_i = 0._wp
          smb%mon_s(m)%dswd_dalb_vis_dif_i = 0._wp
          smb%mon_s(m)%dswd_dalb_nir_dif_i = 0._wp
          smb%mon_s(m)%dswd_dz_nir_dir_i = 0._wp
          smb%mon_s(m)%dswd_dz_nir_dif_i = 0._wp
          smb%mon_s(m)%coszm_i          = 0._wp 
          smb%mon_s(m)%lwdown_i         = 0._wp 
          smb%mon_s(m)%gam_lw_i         = 0._wp 
          smb%mon_s(m)%tam              = 0._wp 
          smb%mon_s(m)%t2m              = 0._wp 
          smb%mon_s(m)%t_skin           = 0._wp 
          smb%mon_s(m)%tskin_tam        = 0._wp 
          smb%mon_s(m)%t_skin_amp       = 0._wp 
          smb%mon_s(m)%t_prof           = 0._wp 
          smb%mon_s(m)%q2m              = 0._wp 
          smb%mon_s(m)%pressure         = 0._wp 
          smb%mon_s(m)%mask_snow        = 0._wp
          smb%mon_s(m)%f_snow           = 0._wp 
          smb%mon_s(m)%h_snow           = 0._wp 
          smb%mon_s(m)%w_snow           = 0._wp 
          smb%mon_s(m)%w_snow_max       = 0._wp 
          smb%mon_s(m)%snowmelt         = 0._wp 
          smb%mon_s(m)%icemelt          = 0._wp 
          smb%mon_s(m)%refreezing       = 0._wp 
          smb%mon_s(m)%f_rfz_to_snow= 0._wp 
          smb%mon_s(m)%runoff           = 0._wp 
          smb%mon_s(m)%albedo           = 0._wp 
          smb%mon_s(m)%alb_vis_dir      = 0._wp 
          smb%mon_s(m)%alb_nir_dir      = 0._wp 
          smb%mon_s(m)%alb_vis_dif      = 0._wp 
          smb%mon_s(m)%alb_nir_dif      = 0._wp 
          smb%mon_s(m)%alb_snow_vis_dir = 0._wp 
          smb%mon_s(m)%alb_snow_nir_dir = 0._wp 
          smb%mon_s(m)%alb_snow_vis_dif = 0._wp 
          smb%mon_s(m)%alb_snow_nir_dif = 0._wp 
          smb%mon_s(m)%snow_grain       = 0._wp 
          smb%mon_s(m)%dust_con         = 0._wp 
          smb%mon_s(m)%cld              = 0._wp 
          smb%mon_s(m)%swnet            = 0._wp 
          smb%mon_s(m)%swnet_min        = 0._wp 
          smb%mon_s(m)%swdown           = 0._wp 
          smb%mon_s(m)%lwdown           = 0._wp 
          smb%mon_s(m)%r_a              = 0._wp 
          smb%mon_s(m)%flx_g            = 0._wp 
          smb%mon_s(m)%dflxg_dT         = 0._wp 
          smb%mon_s(m)%flx_melt         = 0._wp 
          smb%mon_s(m)%flx_sh           = 0._wp 
          smb%mon_s(m)%flx_lwu          = 0._wp 
          smb%mon_s(m)%flx_lh           = 0._wp 
          smb%mon_s(m)%evp              = 0._wp 
          smb%mon_s(m)%prc              = 0._wp 
          smb%mon_s(m)%f_wind           = 0._wp 
          smb%mon_s(m)%rain             = 0._wp 
          smb%mon_s(m)%snow             = 0._wp 
          smb%mon_s(m)%u700             = 0._wp 
          smb%mon_s(m)%v700             = 0._wp 
          smb%mon_s(m)%wind             = 0._wp 
        enddo
      endif

      where (smb%mask_smb.eq.1) 
        smb%mon_s(mon)%t2m_i            = smb%mon_s(mon)%t2m_i            + smb%t2m_i          * mon_avg
        smb%mon_s(mon)%t2m_bias_i       = smb%mon_s(mon)%t2m_bias_i       + smb%t2m_bias_i(:,:,doy) * mon_avg
        smb%mon_s(mon)%tam_i            = smb%mon_s(mon)%tam_i            + smb%tam_i          * mon_avg
        smb%mon_s(mon)%gam_i            = smb%mon_s(mon)%gam_i            + smb%gam_i          * mon_avg
        smb%mon_s(mon)%ram_i            = smb%mon_s(mon)%ram_i            + smb%ram_i          * mon_avg
        smb%mon_s(mon)%prc_i            = smb%mon_s(mon)%prc_i            + smb%prc_i*sec_day          * mon_avg
        smb%mon_s(mon)%prc_bias_i       = smb%mon_s(mon)%prc_bias_i       + smb%prc_bias_i(:,:,doy) * mon_avg
        smb%mon_s(mon)%u700_i           = smb%mon_s(mon)%u700_i           + smb%u700_i         * mon_avg
        smb%mon_s(mon)%v700_i           = smb%mon_s(mon)%v700_i           + smb%v700_i         * mon_avg
        smb%mon_s(mon)%wind_i           = smb%mon_s(mon)%wind_i           + smb%wind_i          * mon_avg
        smb%mon_s(mon)%cld_i            = smb%mon_s(mon)%cld_i            + smb%cld_i          * mon_avg
        smb%mon_s(mon)%dust_i           = smb%mon_s(mon)%dust_i           + smb%dust_i          * mon_avg
        smb%mon_s(mon)%swdown_i         = smb%mon_s(mon)%swdown_i + ((1._wp-smb%cld_i)*(frac_vu*smb%swd_sur_vis_dir_i+(1._wp-frac_vu)*smb%swd_sur_nir_dir_i) &
           + smb%cld_i * (frac_vu*smb%swd_sur_vis_dif_i+(1._wp-frac_vu)*smb%swd_sur_nir_dif_i)) * mon_avg
        smb%mon_s(mon)%swd_toa_i     = smb%mon_s(mon)%swd_toa_i     + smb%swd_toa_i          * mon_avg
        smb%mon_s(mon)%swd_sur_vis_dir_i = smb%mon_s(mon)%swd_sur_vis_dir_i + smb%swd_sur_vis_dir_i          * mon_avg
        smb%mon_s(mon)%swd_sur_nir_dir_i = smb%mon_s(mon)%swd_sur_nir_dir_i + smb%swd_sur_nir_dir_i          * mon_avg
        smb%mon_s(mon)%swd_sur_vis_dif_i = smb%mon_s(mon)%swd_sur_vis_dif_i + smb%swd_sur_vis_dif_i          * mon_avg
        smb%mon_s(mon)%swd_sur_nir_dif_i = smb%mon_s(mon)%swd_sur_nir_dif_i + smb%swd_sur_nir_dif_i          * mon_avg
        smb%mon_s(mon)%dswd_dalb_vis_dir_i = smb%mon_s(mon)%dswd_dalb_vis_dir_i + smb%dswd_dalb_vis_dir_i * mon_avg   
        smb%mon_s(mon)%dswd_dalb_nir_dir_i = smb%mon_s(mon)%dswd_dalb_nir_dir_i + smb%dswd_dalb_nir_dir_i * mon_avg
        smb%mon_s(mon)%dswd_dalb_vis_dif_i = smb%mon_s(mon)%dswd_dalb_vis_dif_i + smb%dswd_dalb_vis_dif_i * mon_avg
        smb%mon_s(mon)%dswd_dalb_nir_dif_i = smb%mon_s(mon)%dswd_dalb_nir_dif_i + smb%dswd_dalb_nir_dif_i * mon_avg
        smb%mon_s(mon)%dswd_dz_nir_dir_i   = smb%mon_s(mon)%dswd_dz_nir_dir_i   + smb%dswd_dz_nir_dir_i   * mon_avg
        smb%mon_s(mon)%dswd_dz_nir_dif_i   = smb%mon_s(mon)%dswd_dz_nir_dif_i   + smb%dswd_dz_nir_dif_i   * mon_avg
        smb%mon_s(mon)%coszm_i          = smb%mon_s(mon)%coszm_i          + smb%coszm_i          * mon_avg
        smb%mon_s(mon)%lwdown_i         = smb%mon_s(mon)%lwdown_i         + smb%lwdown_i      * mon_avg
        smb%mon_s(mon)%gam_lw_i         = smb%mon_s(mon)%gam_lw_i         + smb%gam_lw_i         * mon_avg
        smb%mon_s(mon)%tam              = smb%mon_s(mon)%tam              + smb%tam          * mon_avg
        smb%mon_s(mon)%t2m              = smb%mon_s(mon)%t2m              + smb%t2m          * mon_avg
        smb%mon_s(mon)%t_skin           = smb%mon_s(mon)%t_skin           + smb%t_skin          * mon_avg
        smb%mon_s(mon)%tskin_tam        = smb%mon_s(mon)%tskin_tam        + (smb%t_skin-smb%tam)       * mon_avg
        smb%mon_s(mon)%t_skin_amp       = smb%mon_s(mon)%t_skin_amp       + smb%t_skin_amp      * mon_avg
        smb%mon_s(mon)%q2m              = smb%mon_s(mon)%q2m              + smb%q2m          * mon_avg
        smb%mon_s(mon)%pressure         = smb%mon_s(mon)%pressure         + smb%pressure          * mon_avg
        smb%mon_s(mon)%mask_snow        = smb%mon_s(mon)%mask_snow        + smb%mask_snow          * mon_avg
        smb%mon_s(mon)%f_snow           = smb%mon_s(mon)%f_snow           + smb%f_snow          * mon_avg
        smb%mon_s(mon)%h_snow           = smb%mon_s(mon)%h_snow           + smb%h_snow          * mon_avg
        smb%mon_s(mon)%w_snow           = smb%mon_s(mon)%w_snow           + smb%w_snow          * mon_avg
        smb%mon_s(mon)%w_snow_max       = smb%mon_s(mon)%w_snow_max       + smb%w_snow_max          * mon_avg
        smb%mon_s(mon)%snowmelt         = smb%mon_s(mon)%snowmelt         + smb%snowmelt*sec_day          * mon_avg
        smb%mon_s(mon)%icemelt          = smb%mon_s(mon)%icemelt          + smb%icemelt*sec_day          * mon_avg
        smb%mon_s(mon)%refreezing       = smb%mon_s(mon)%refreezing       + smb%refreezing*sec_day          * mon_avg
        smb%mon_s(mon)%f_rfz_to_snow       = smb%mon_s(mon)%f_rfz_to_snow       + smb%f_rfz_to_snow          * mon_avg
        smb%mon_s(mon)%runoff           = smb%mon_s(mon)%runoff           + smb%runoff*sec_day          * mon_avg
        smb%mon_s(mon)%albedo           = smb%mon_s(mon)%albedo           + smb%albedo          * mon_avg
        smb%mon_s(mon)%alb_vis_dir      = smb%mon_s(mon)%alb_vis_dir      + smb%alb_vis_dir          * mon_avg
        smb%mon_s(mon)%alb_nir_dir      = smb%mon_s(mon)%alb_nir_dir      + smb%alb_nir_dir          * mon_avg
        smb%mon_s(mon)%alb_vis_dif      = smb%mon_s(mon)%alb_vis_dif      + smb%alb_vis_dif          * mon_avg
        smb%mon_s(mon)%alb_nir_dif      = smb%mon_s(mon)%alb_nir_dif      + smb%alb_nir_dif          * mon_avg
        smb%mon_s(mon)%alb_snow_vis_dir = smb%mon_s(mon)%alb_snow_vis_dir + smb%alb_snow_vis_dir* mon_avg
        smb%mon_s(mon)%alb_snow_nir_dir = smb%mon_s(mon)%alb_snow_nir_dir + smb%alb_snow_nir_dir* mon_avg
        smb%mon_s(mon)%alb_snow_vis_dif = smb%mon_s(mon)%alb_snow_vis_dif + smb%alb_snow_vis_dif* mon_avg
        smb%mon_s(mon)%alb_snow_nir_dif = smb%mon_s(mon)%alb_snow_nir_dif + smb%alb_snow_nir_dif* mon_avg
        smb%mon_s(mon)%snow_grain       = smb%mon_s(mon)%snow_grain       + smb%snow_grain* mon_avg
        smb%mon_s(mon)%dust_con         = smb%mon_s(mon)%dust_con         + smb%dust_con*1.d6* mon_avg
        smb%mon_s(mon)%cld              = smb%mon_s(mon)%cld              + smb%cld          * mon_avg
        smb%mon_s(mon)%swnet            = smb%mon_s(mon)%swnet            + smb%swnet          * mon_avg
        smb%mon_s(mon)%swnet_min        = smb%mon_s(mon)%swnet_min        + smb%swnet_min      * mon_avg
        smb%mon_s(mon)%swdown           = smb%mon_s(mon)%swdown           + smb%swdown         * mon_avg
        smb%mon_s(mon)%lwdown           = smb%mon_s(mon)%lwdown           + smb%lwdown          * mon_avg
        smb%mon_s(mon)%r_a              = smb%mon_s(mon)%r_a              + smb%r_a          * mon_avg
        smb%mon_s(mon)%flx_g            = smb%mon_s(mon)%flx_g            + smb%flx_g          * mon_avg
        smb%mon_s(mon)%dflxg_dT         = smb%mon_s(mon)%dflxg_dT         + smb%dflxg_dT          * mon_avg
        smb%mon_s(mon)%flx_melt         = smb%mon_s(mon)%flx_melt         + smb%flx_melt          * mon_avg
        smb%mon_s(mon)%flx_sh           = smb%mon_s(mon)%flx_sh           + smb%flx_sh          * mon_avg
        smb%mon_s(mon)%flx_lwu          = smb%mon_s(mon)%flx_lwu          + smb%flx_lwu          * mon_avg
        smb%mon_s(mon)%flx_lh           = smb%mon_s(mon)%flx_lh           + smb%flx_lh          * mon_avg
        smb%mon_s(mon)%evp              = smb%mon_s(mon)%evp              + smb%evp*sec_day          * mon_avg
        smb%mon_s(mon)%prc              = smb%mon_s(mon)%prc              + smb%prc*sec_day          * mon_avg
        smb%mon_s(mon)%f_wind           = smb%mon_s(mon)%f_wind           + smb%f_wind          * mon_avg
        smb%mon_s(mon)%rain             = smb%mon_s(mon)%rain             + smb%rain*sec_day          * mon_avg
        smb%mon_s(mon)%snow             = smb%mon_s(mon)%snow             + smb%snow*sec_day          * mon_avg
        smb%mon_s(mon)%u700             = smb%mon_s(mon)%u700             + smb%u700          * mon_avg
        smb%mon_s(mon)%v700             = smb%mon_s(mon)%v700             + smb%v700          * mon_avg
        smb%mon_s(mon)%wind             = smb%mon_s(mon)%wind             + smb%wind          * mon_avg
      endwhere

      if (l_daily_output) then
        where (smb%mask_smb.eq.1) 
          smb%day_s(doy)%t2m              = smb%t2m
          smb%day_s(doy)%t_skin           = smb%t_skin         
          smb%day_s(doy)%t_skin_amp       = smb%t_skin_amp     
          smb%day_s(doy)%mask_snow        = smb%mask_snow      
          smb%day_s(doy)%h_snow           = smb%h_snow         
          smb%day_s(doy)%w_snow           = smb%w_snow         
          smb%day_s(doy)%w_snow_max       = smb%w_snow_max     
          smb%day_s(doy)%snowmelt         = smb%snowmelt*sec_day  
          smb%day_s(doy)%icemelt          = smb%icemelt*sec_day  
          smb%day_s(doy)%refreezing       = smb%refreezing*sec_day 
          smb%day_s(doy)%flx_g            = smb%flx_g         
          smb%day_s(doy)%dflxg_dT         = smb%dflxg_dT      
          smb%day_s(doy)%flx_melt         = smb%flx_melt      
          smb%day_s(doy)%flx_sh           = smb%flx_sh        
          smb%day_s(doy)%flx_lwu          = smb%flx_lwu       
          smb%day_s(doy)%flx_lh           = smb%flx_lh        
          smb%day_s(doy)%rain             = smb%rain*sec_day  
          smb%day_s(doy)%snow             = smb%snow*sec_day  
          smb%day_s(doy)%dust_con         = smb%dust_con  
        endwhere
        do k=0,nl
          where (smb%mask_smb.eq.1)
            smb%day_s(doy)%t_prof(:,:,k)           =  smb%t_prof(:,:,k)
          endwhere
        enddo
      endif

      do k=0,nl
        where (smb%mask_smb.eq.1)
          smb%mon_s(mon)%t_prof(:,:,k)           = smb%mon_s(mon)%t_prof(:,:,k)           + smb%t_prof(:,:,k)          * mon_avg
        endwhere
      enddo

    endif

    if (time_out_smb .and. time_eoy_smb) then

      smb%ann_s%area             = smb%grid%area
      smb%ann_s%dz_dx_sur        = smb%dz_dx_sur 
      smb%ann_s%dz_dy_sur        = smb%dz_dy_sur
      smb%ann_s%dz_sur           = smb%dz_sur  
      smb%ann_s%mask_smb         = smb%mask_smb
      smb%ann_s%mask_ice         = smb%mask_ice  
      smb%ann_s%f_ice            = smb%f_ice
      smb%ann_s%h_ice            = smb%h_ice
      smb%ann_s%mask_maxice      = smb%mask_maxice
      smb%ann_s%mask_margin      = smb%mask_margin
      smb%ann_s%z_sur            = smb%z_sur  
      smb%ann_s%z_sur_eff        = smb%z_sur_eff
      smb%ann_s%z_sur_fil        = smb%z_sur_fil
      smb%ann_s%z_sur_std        = smb%z_sur_std
      smb%ann_s%z_sur_i          = smb%z_sur_i  
      smb%ann_s%f_ele            = smb%f_ele
      smb%ann_s%alb_ice          = smb%alb_ice
      smb%ann_s%alb_bg           = smb%alb_bg

      smb%ann_s%t_ice            = smb%t_ice    
      smb%ann_s%t_ground_i       = smb%t_ground_i
      smb%ann_s%t_ground         = smb%t_ground
      smb%ann_s%pdd              = smb%simple%pdd     ! degC

      smb%ann_s%ann_smb_pdd      = smb%simple%smb*sec_year
      smb%ann_s%ann_smb          = smb%ann_smb            
      smb%ann_s%ann_smb_ice      = smb%ann_smb
      where (smb%mask_ice.eq.0) smb%ann_s%ann_smb_ice = 0._wp
      smb%ann_s%ann_smb_noice    = smb%ann_smb
      where (smb%mask_ice.eq.1) smb%ann_s%ann_smb_noice = 0._wp
      smb%ann_s%ann_prc          = smb%ann_prc
      smb%ann_s%ann_prc_ice      = smb%ann_prc
      where (smb%mask_ice.eq.0) smb%ann_s%ann_prc_ice = 0._wp
      smb%ann_s%ann_snow         = smb%ann_snow
      smb%ann_s%ann_snow_ice     = smb%ann_snow
      where (smb%mask_ice.eq.0) smb%ann_s%ann_snow_ice = 0._wp
      smb%ann_s%ann_ablation     = smb%ann_ablation
      smb%ann_s%ann_ablation_ice = smb%ann_ablation
      where (smb%mask_ice.eq.0) smb%ann_s%ann_ablation_ice = 0._wp
      smb%ann_s%ann_melt         = smb%ann_melt
      smb%ann_s%ann_melt_ice     = smb%ann_melt
      where (smb%mask_ice.eq.0) smb%ann_s%ann_melt_ice = 0._wp
      smb%ann_s%ann_evp          = smb%ann_evp
      smb%ann_s%ann_evp_ice      = smb%ann_evp
      where (smb%mask_ice.eq.0) smb%ann_s%ann_evp_ice = 0._wp
      smb%ann_s%ann_runoff       = smb%ann_runoff
      smb%ann_s%ann_runoff_ice   = smb%ann_runoff
      where (smb%mask_ice.eq.0) smb%ann_s%ann_runoff_ice = 0._wp
      smb%ann_s%ann_refreezing   = smb%ann_refreezing
      smb%ann_s%ann_refreezing_ice = smb%ann_refreezing
      where (smb%mask_ice.eq.0) smb%ann_s%ann_refreezing_ice = 0._wp

      smb%nout = smb%nout+1
      call smb_diag_out(smb,smb%grid)
    endif


   return

  end subroutine smb_diag


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s m b _ d i a g _ o u t
  ! Purpose  :  write sea ice netcdf output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_diag_out(smb,grid)

    implicit none

    type(smb_class) :: smb
    type(grid_class) :: grid
    integer :: k, ncid
    character (len=256) :: fnm


    ! Get annual values
    call smb_ave( smb%mon_s,smb%ann_s )

    ! write to file
    fnm = trim(out_dir)//"/smb_"//trim(smb%grid%name)//".nc"
    call nc_open(fnm,ncid)
    call nc_write(fnm,dim_time,dble(year_now), dim1=dim_time, start=[smb%nout], count=[1],ncid=ncid)    
    do k = 1, nmon_year
       call smb_nc_write(fnm,ncid,smb%mon_s(k),grid,k,smb%nout)
    end do
    call smb_nc_write(fnm,ncid,smb%ann_s,grid,nmon_year+1,smb%nout)
    call nc_close(ncid)

    if (l_daily_output) then
      ! write to file
      fnm = trim(out_dir)//"/smb_"//trim(smb%grid%name)//"_daily.nc"
      call nc_open(fnm,ncid)
      call nc_write(fnm,dim_time,dble(year_now), dim1=dim_time, start=[smb%nout], count=[1],ncid=ncid)    
      do k = 1, nday_year
        call smb_daily_nc_write(fnm,ncid,smb%day_s(k),grid,k,smb%nout)
      end do
      call nc_close(ncid)
    endif


   return

  end subroutine smb_diag_out
  

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
    call nc_write_dim(fnm, dim_x, x=1, axis="x", units="1")
    call nc_write_dim(fnm, dim_y, x=1, axis="y", units="1")

    return

  end subroutine ts_nc
 

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t s _ n c _ w r i t e
  ! Purpose  :  write time series to netcdf
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ts_nc_write(fnm,vars,nout,y)

    implicit none

    type(ts_out) :: vars(:)

    character (len=*) :: fnm
    integer :: nout, ncid, y, i

    call nc_open(fnm,ncid)
    call nc_write(fnm,"time", dble([(i,i=(year_now-(y-1)*n_year_smb),(year_now),n_year_smb)]), &
    dim1=dim_time,start=[nout],count=[y],ncid=ncid)
    call nc_write(fnm,"Aice", vars%Aice, dims=[dim_time],start=[nout],count=[y],&
    long_name="total ice sheet area",units="mln km2",ncid=ncid)
    call nc_write(fnm,"smb", vars%smb, dims=[dim_time],start=[nout],count=[y],&
    long_name="integrated surface mass balance of ice sheets",units="Gt/yr",ncid=ncid)
    call nc_write(fnm,"prc", vars%prc, dims=[dim_time],start=[nout],count=[y],&
    long_name="integrated precipitation over ice sheets",units="Gt/yr",ncid=ncid)
    call nc_write(fnm,"snow", vars%snow, dims=[dim_time],start=[nout],count=[y],&
    long_name="integrated snowfall on ice sheets",units="Gt/yr",ncid=ncid)
    call nc_write(fnm,"melt", vars%melt, dims=[dim_time],start=[nout],count=[y],&
    long_name="total melt (ice+snow) of ice sheets",units="Gt/yr",ncid=ncid)
    call nc_write(fnm,"run", vars%run, dims=[dim_time],start=[nout],count=[y],&
    long_name="integrated runoff from ice sheets",units="Gt/yr",ncid=ncid)
    call nc_write(fnm,"rfz", vars%refreezing, dims=[dim_time],start=[nout],count=[y],&
    long_name="integrated refreezing on ice sheets",units="Gt/yr",ncid=ncid)
    call nc_write(fnm,"subl", vars%evp, dims=[dim_time],start=[nout],count=[y],&
    long_name="integrated sublimation ice sheets",units="Gt/yr",ncid=ncid)
    call nc_write(fnm,"smb_avg", vars%smb_avg, dims=[dim_time],start=[nout],count=[y],&
    long_name="averaged surface mass balance of ice sheets",units="kg/m2/yr",ncid=ncid)
    call nc_write(fnm,"prc_avg", vars%prc_avg, dims=[dim_time],start=[nout],count=[y],&
    long_name="averaged precipitation over ice sheets",units="kg/m2/yr",ncid=ncid)
    call nc_write(fnm,"snow_avg", vars%snow_avg, dims=[dim_time],start=[nout],count=[y],&
    long_name="averaged snowfall on ice sheets",units="kg/m2/yr",ncid=ncid)
    call nc_write(fnm,"melt_avg", vars%melt_avg, dims=[dim_time],start=[nout],count=[y],&
    long_name="total melt (ice+snow) of ice sheets",units="kg/m2/yr",ncid=ncid)
    call nc_write(fnm,"run_avg", vars%run_avg, dims=[dim_time],start=[nout],count=[y],&
    long_name="averaged runoff from ice sheets",units="kg/m2/yr",ncid=ncid)
    call nc_write(fnm,"rfz_avg", vars%refreezing_avg, dims=[dim_time],start=[nout],count=[y],&
    long_name="averaged refreezing on ice sheets",units="kg/m2/yr",ncid=ncid)
    call nc_write(fnm,"subl_avg", vars%evp_avg, dims=[dim_time],start=[nout],count=[y],&
    long_name="averaged sublimation ice sheets",units="kg/m2/yr",ncid=ncid)
    call nc_close(ncid)


   return

  end subroutine ts_nc_write


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s m b _ n c
  ! Purpose  :  Initialize smb netcdf output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_nc(fnm,grid)

    implicit none

    character (len=*) :: fnm
    type(grid_class) :: grid
    integer :: ncid
    real(wp) :: empty_time(0)


    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,dim_time, x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE.,ncid=ncid)
    call nc_write_dim(fnm, dim_month, x=1._wp, dx=1._wp, nx=13, axis="e", &
    units="months", ncid=ncid)
    call nc_write_dim(fnm, dim_y, x=grid%G%y0, dx=grid%G%dy, nx=grid%G%ny,&
    axis="y", units="km", ncid=ncid)
    call nc_write_dim(fnm, dim_x, x=grid%G%x0, dx=grid%G%dx, nx=grid%G%nx,&
    axis="x", units="km", ncid=ncid)
    call nc_write_dim(fnm,dim_depth,x=z,units="m",ncid=ncid)

    call nc_write(fnm,"phi",sngl(grid%lat), dims=[dim_x,dim_y],start=[1,1],count=[grid%G%nx,grid%G%ny],ncid=ncid)
    call nc_write(fnm,"lambda",sngl(grid%lon), dims=[dim_x,dim_y],start=[1,1],count=[grid%G%nx,grid%G%ny],ncid=ncid)
    call nc_close(ncid)

   return

  end subroutine smb_nc


  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s m b _ n c _ w r i t e
  ! Purpose  :  Output of smb netcdf
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_nc_write(fnm,ncid,vars,grid,ndat,nout)

    implicit none

    type(s_out) :: vars
    type(grid_class) :: grid

    character (len=*) :: fnm
    integer :: ndat, nout, ncid

if (l_monthly_output) then

if (i_smb.ne.1) then
call nc_write(fnm,"t2m_i", sngl(vars%t2m_i), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated 2m temperature",grid_mapping="polar_stereographic",units="K",ncid=ncid) 
endif
call nc_write(fnm,"t2m", sngl(vars%t2m             ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="2m temperature",grid_mapping="polar_stereographic",units="K",ncid=ncid) 
if (i_smb.ne.3) then
call nc_write(fnm,"prc_i", sngl(vars%prc_i           ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated precipitation",grid_mapping="polar_stereographic",units="kg/m2/day",ncid=ncid) 
call nc_write(fnm,"u700_i", sngl(vars%u700_i             ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated zonal wind velocity at 700 hPa",grid_mapping="polar_stereographic",units="m/s",ncid=ncid) 
call nc_write(fnm,"v700_i", sngl(vars%v700_i             ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated meridional wind velocity at 700 hPa",grid_mapping="polar_stereographic",units="m/s",ncid=ncid) 
call nc_write(fnm,"wind_i", sngl(vars%wind_i          ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated wind speed",grid_mapping="polar_stereographic",units="m/s",ncid=ncid) 
call nc_write(fnm,"u700", sngl(vars%u700               ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="zonal wind velocity at 700 hPa",grid_mapping="polar_stereographic",units="m/s",ncid=ncid) 
call nc_write(fnm,"v700", sngl(vars%v700               ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="meridional wind velocity at 700 hPa",grid_mapping="polar_stereographic",units="m/s",ncid=ncid) 
call nc_write(fnm,"wind", sngl(vars%wind            ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="wind speed",grid_mapping="polar_stereographic",units="m/s",ncid=ncid) 
call nc_write(fnm,"prc", sngl(vars%prc             ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="precipitation",grid_mapping="polar_stereographic",units="kg/m2/day",ncid=ncid) 
call nc_write(fnm,"f_wind", sngl(vars%f_wind             ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="wind factor for precipitation downscaling",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
call nc_write(fnm,"rain", sngl(vars%rain            ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="rainfall rate",grid_mapping="polar_stereographic",units="kg/m2/day",ncid=ncid) 
call nc_write(fnm,"snow", sngl(vars%snow            ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="snowfall rate",grid_mapping="polar_stereographic",units="kg/m2/day",ncid=ncid) 
endif

if (i_smb==1) then

call nc_write(fnm,"cld_i", sngl(vars%cld_i           ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated total cloud cover fraction",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
call nc_write(fnm,"t2m_bias_i", sngl(vars%t2m_bias_i            ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated 2m temperature bias for present-day",grid_mapping="polar_stereographic",units="K",ncid=ncid) 
call nc_write(fnm,"tam_i", sngl(vars%tam_i            ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated atmospheric temperature",grid_mapping="polar_stereographic",units="K",ncid=ncid) 
call nc_write(fnm,"gam_i", sngl(vars%gam_i            ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated atmospheric temperature lapse rate",grid_mapping="polar_stereographic",units="K/m",ncid=ncid) 
call nc_write(fnm,"ram_i", sngl(vars%ram_i            ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated atmospheric relative humidity",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
call nc_write(fnm,"prc_bias_i", sngl(vars%prc_bias_i            ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated precipitation bias for present-day (ratio of model/obs)",grid_mapping="polar_stereographic",units="1",ncid=ncid) 
call nc_write(fnm,"dust_i", sngl(vars%dust_i          ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated dust deposition rate",grid_mapping="polar_stereographic",units="kg/m2/s",ncid=ncid) 
call nc_write(fnm,"swdown_i", sngl(vars%swdown_i    ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated downward shortwave radiation at the surface",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"swd_toa_i", sngl(vars%swd_toa_i    ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated downward shortwave radiation at TOA",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"swd_sur_vis_dir_i", sngl(vars%swd_sur_vis_dir_i), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated downward direct shortwave radiation at the surface",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"swd_sur_nir_dir_i", sngl(vars%swd_sur_nir_dir_i), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated downward direct shortwave radiation at the surface",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"swd_sur_vis_dif_i", sngl(vars%swd_sur_vis_dif_i), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated downward diffude shortwave radiation at the surface",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"swd_sur_nir_dif_i", sngl(vars%swd_sur_nir_dif_i), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated downward diffude shortwave radiation at the surface",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 

call nc_write(fnm,"dswd_dalb_vis_dir_i ", sngl(vars%dswd_dalb_vis_dir_i ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"dswd_dalb_nir_dir_i ", sngl(vars%dswd_dalb_nir_dir_i ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"dswd_dalb_vis_dif_i ", sngl(vars%dswd_dalb_vis_dif_i ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"dswd_dalb_nir_dif_i ", sngl(vars%dswd_dalb_nir_dif_i ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"dswd_dz_nir_dir_i   ", sngl(vars%dswd_dz_nir_dir_i   ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"dswd_dz_nir_dif_i   ", sngl(vars%dswd_dz_nir_dif_i   ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 

call nc_write(fnm,"coszm_i", sngl(vars%coszm_i         ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated cosine of the solar zenith angle",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
call nc_write(fnm,"lwdown_i", sngl(vars%lwdown_i            ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated downward longwave radiation",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"gam_lw_i", sngl(vars%gam_lw_i            ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="interpolated downward longwave radiation",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"tam", sngl(vars%tam           ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="atmospheric temperature",grid_mapping="polar_stereographic",units="K",ncid=ncid) 
call nc_write(fnm,"t2m", sngl(vars%t2m           ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="surface air temperature",grid_mapping="polar_stereographic",units="K",ncid=ncid) 
call nc_write(fnm,"t_skin", sngl(vars%t_skin          ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="skin temperature",grid_mapping="polar_stereographic",units="K",ncid=ncid) 
call nc_write(fnm,"tskin_tam", sngl(vars%tskin_tam          ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="skin - atmospheric temperature",grid_mapping="polar_stereographic",units="K",ncid=ncid) 
call nc_write(fnm,"t_skin_amp", sngl(vars%t_skin_amp          ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="amplitude of diurnal cycle of skin temperature",grid_mapping="polar_stereographic",units="K",ncid=ncid) 
call nc_write(fnm,"t_prof", sngl(vars%t_prof          ), dims=[dim_x,dim_y,dim_depth,dim_month,dim_time],start=[1,1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,nl+1,1,1],long_name="temperature profile in snow+ice/soil",grid_mapping="polar_stereographic",units="K",ncid=ncid) 
call nc_write(fnm,"q2m", sngl(vars%q2m            ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="surface air humidity",grid_mapping="polar_stereographic",units="kg/kg",ncid=ncid) 
call nc_write(fnm,"pressure", sngl(vars%pressure        ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="surface pressure",grid_mapping="polar_stereographic",units="Pa",ncid=ncid) 
call nc_write(fnm,"mask_snow", sngl(vars%mask_snow       ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="snow mask",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
call nc_write(fnm,"f_snow", sngl(vars%f_snow          ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="snow cover fraction",grid_mapping="polar_stereographic",units="m",ncid=ncid) 
call nc_write(fnm,"h_snow", sngl(vars%h_snow          ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="snow thickness",grid_mapping="polar_stereographic",units="m",ncid=ncid) 
call nc_write(fnm,"w_snow", sngl(vars%w_snow          ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="snow water equivalent",grid_mapping="polar_stereographic",units="kg/m2",ncid=ncid) 
call nc_write(fnm,"w_snow_max", sngl(vars%w_snow_max      ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="seasonal max of snow water equivalent",grid_mapping="polar_stereographic",units="kg/m2",ncid=ncid) 
call nc_write(fnm,"snowmelt", sngl(vars%snowmelt        ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="snow melt rate",grid_mapping="polar_stereographic",units="kg/m2/day",ncid=ncid) 
call nc_write(fnm,"icemelt", sngl(vars%icemelt         ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="ice melt rate",grid_mapping="polar_stereographic",units="kg/m2/day",ncid=ncid) 
call nc_write(fnm,"runoff", sngl(vars%runoff         ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="runoff",grid_mapping="polar_stereographic",units="kg/m2/day",ncid=ncid) 
call nc_write(fnm,"refreezing", sngl(vars%refreezing        ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="refreezing",grid_mapping="polar_stereographic",units="kg/m2/day",ncid=ncid) 
call nc_write(fnm,"f_rfz_to_snow", sngl(vars%f_rfz_to_snow        ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="fraction of refreezing going into snow",grid_mapping="polar_stereographic",units="1",ncid=ncid) 
call nc_write(fnm,"albedo", sngl(vars%albedo          ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="surface albedo",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
call nc_write(fnm,"alb_vis_dir", sngl(vars%alb_vis_dir     ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="clear-sky visible surface albedo",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
call nc_write(fnm,"alb_nir_dir", sngl(vars%alb_nir_dir     ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="clear-sky near-infrared surface albedo",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
call nc_write(fnm,"alb_vis_dif", sngl(vars%alb_vis_dif     ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="diffuse visible surface albedo",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
call nc_write(fnm,"alb_nir_dif", sngl(vars%alb_nir_dif     ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="diffuse near-infrared surface albedo",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
call nc_write(fnm,"alb_snow_vis_dir", sngl(vars%alb_snow_vis_dir), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="clear-sky visible snow albedo",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
call nc_write(fnm,"alb_snow_nir_dir", sngl(vars%alb_snow_nir_dir), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="clear-sky near-infrared snow albedo",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
call nc_write(fnm,"alb_snow_vis_dif", sngl(vars%alb_snow_vis_dif), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="diffuse visible snow albedo",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
call nc_write(fnm,"alb_snow_nir_dif", sngl(vars%alb_snow_nir_dif), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="diffuse near-infrared snow albedo",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
call nc_write(fnm,"snow_grain", sngl(vars%snow_grain        ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="snow grain size",grid_mapping="polar_stereographic",units="um",ncid=ncid) 
call nc_write(fnm,"dust_con", sngl(vars%dust_con        ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="dust concentration in surface snow",grid_mapping="polar_stereographic",units="mg/kg",ncid=ncid) 
call nc_write(fnm,"cld", sngl(vars%cld           ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="cloud cover fraction",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"swnet", sngl(vars%swnet           ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="net surface shortwave radiation",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"swnet_min", sngl(vars%swnet_min          ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="minimum diurnal net surface shortwave radiation",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"swdown", sngl(vars%swdown          ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="downward surface shortwave radiation",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"lwdown", sngl(vars%lwdown          ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="downward longwave radiation at the surface",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"r_a", sngl(vars%r_a             ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="aerodynamic resistance",grid_mapping="polar_stereographic",units="s/m?",ncid=ncid) 
call nc_write(fnm,"flx_g", sngl(vars%flx_g           ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="ground heat flux",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"dflxg_dT", sngl(vars%dflxg_dT        ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="derivative of ground heat flux wrt top soil temperature",grid_mapping="polar_stereographic",units="W/m2/K",ncid=ncid) 
call nc_write(fnm,"flx_melt", sngl(vars%flx_melt        ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="heat flux going into the melt of snow",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"flx_sh", sngl(vars%flx_sh          ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="sensible heat flux",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"flx_lwu", sngl(vars%flx_lwu         ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="upward longwave radiation",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"flx_lh", sngl(vars%flx_lh          ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="latent_heat_flux",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"evp", sngl(vars%evp             ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="evaporation/sublimation",grid_mapping="polar_stereographic",units="kg/m2/day",ncid=ncid) 
endif

endif

if (ndat.eq.13) then

if (i_smb.eq.1) then
call nc_write(fnm,"mask_smb", real(vars%mask_smb       ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="mask where semi is applied",grid_mapping="polar_stereographic",units="/",ncid=ncid)
call nc_write(fnm,"smb_pdd", sngl(vars%ann_smb_pdd        ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="annual surface mass balance using PDD scheme",grid_mapping="polar_stereographic",units="kg/m2/yr",ncid=ncid) 
call nc_write(fnm,"alb_bg", sngl(vars%alb_bg          ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="backround (ice/bare soil) albedo",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
call nc_write(fnm,"alb_ice", sngl(vars%alb_ice          ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="ice albedo",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
call nc_write(fnm,"f_ice", sngl(vars%f_ice          ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="sub-grid ice sheet fraction due to orography",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
call nc_write(fnm,"h_ice", sngl(vars%h_ice          ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="ice thickness",grid_mapping="polar_stereographic",units="m",ncid=ncid) 
call nc_write(fnm,"mask_margin", vars%mask_margin         , dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="ice margin mask",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
endif
call nc_write(fnm,"area", sngl(vars%area       ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="grid cell area",grid_mapping="polar_stereographic",units="km^2",ncid=ncid) 
call nc_write(fnm,"mask_ice", vars%mask_ice          , dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="ice mask",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
call nc_write(fnm,"mask_maxice", vars%mask_maxice          , dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="mask of maximum allowed ice sheet extent",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
    call nc_write(fnm,"z_sur_i", sngl(vars%z_sur_i         ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="interpolated surface elevation",grid_mapping="polar_stereographic",units="m",ncid=ncid)    
call nc_write(fnm,"dz_dx_sur", sngl(vars%dz_dx_sur       ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="zonal gradient of surface elevation",grid_mapping="polar_stereographic",units="m/m",ncid=ncid) 
call nc_write(fnm,"dz_dy_sur", sngl(vars%dz_dy_sur       ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="meridional gradient of surface elevation",grid_mapping="polar_stereographic",units="m/m",ncid=ncid) 
call nc_write(fnm,"dz_sur", sngl(vars%dz_sur          ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="surface elevation gradient",grid_mapping="polar_stereographic",units="m/m",ncid=ncid) 
call nc_write(fnm,"z_sur", sngl(vars%z_sur          ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="surface elevation",grid_mapping="polar_stereographic",units="m",ncid=ncid) 
call nc_write(fnm,"z_sur_eff", sngl(vars%z_sur_eff          ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="effective surface elevation",grid_mapping="polar_stereographic",units="m",ncid=ncid) 
call nc_write(fnm,"z_sur_fil", sngl(vars%z_sur_fil ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="smoothed surface elevation",grid_mapping="polar_stereographic",units="m",ncid=ncid) 
call nc_write(fnm,"z_sur_std", sngl(vars%z_sur_std ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="sub-grid standard deviation of surface elevation",grid_mapping="polar_stereographic",units="m",ncid=ncid) 
call nc_write(fnm,"f_ele", sngl(vars%f_ele             ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="elevation factor for precipitation downscaling",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
call nc_write(fnm,"t_ice", sngl(vars%t_ice            ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="surface ice temperature (10 m firn)",grid_mapping="polar_stereographic",units="degC",ncid=ncid) 
call nc_write(fnm,"t_ground_i", sngl(vars%t_ground_i        ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="interpolated ground temperature",grid_mapping="polar_stereographic",units="degC",ncid=ncid) 
call nc_write(fnm,"t_ground", sngl(vars%t_ground        ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="ground temperature (soil temperature over land and bottom water temperatature over ocean)",grid_mapping="polar_stereographic",units="degC",ncid=ncid) 
call nc_write(fnm,"pdd", sngl(vars%pdd            ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="positive degree days",grid_mapping="polar_stereographic",units="degC",ncid=ncid) 
call nc_write(fnm,"smb", sngl(vars%ann_smb        ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="annual surface mass balance",grid_mapping="polar_stereographic",units="kg/m2/yr",ncid=ncid) 
call nc_write(fnm,"smb_ice", sngl(vars%ann_smb_ice        ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="annual surface mass balance, for ice cells only",grid_mapping="polar_stereographic",units="kg/m2/yr",ncid=ncid) 
call nc_write(fnm,"smb_noice", sngl(vars%ann_smb_noice        ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="annual surface mass balance, for ice-free cells only",grid_mapping="polar_stereographic",units="kg/m2/yr",ncid=ncid) 
call nc_write(fnm,"prc_ann", sngl(vars%ann_prc        ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="annual precipitation",grid_mapping="polar_stereographic",units="kg/m2/yr",ncid=ncid) 
call nc_write(fnm,"prc_ann_ice", sngl(vars%ann_prc_ice        ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="annual precipitation, for ice cells only",grid_mapping="polar_stereographic",units="kg/m2/yr",ncid=ncid) 
call nc_write(fnm,"snow_ann", sngl(vars%ann_snow       ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="annual snowfall",grid_mapping="polar_stereographic",units="kg/m2/yr",ncid=ncid) 
call nc_write(fnm,"snow_ann_ice", sngl(vars%ann_snow_ice        ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="annual snowfall, for ice cells only",grid_mapping="polar_stereographic",units="kg/m2/yr",ncid=ncid) 
call nc_write(fnm,"melt_ann", sngl(vars%ann_melt        ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="annual surface melt",grid_mapping="polar_stereographic",units="kg/m2/yr",ncid=ncid) 
call nc_write(fnm,"melt_ann_ice", sngl(vars%ann_melt_ice        ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="annual melt, for ice cells only",grid_mapping="polar_stereographic",units="kg/m2/yr",ncid=ncid) 
call nc_write(fnm,"runoff_ann", sngl(vars%ann_runoff       ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="annual runoff",grid_mapping="polar_stereographic",units="kg/m2/yr",ncid=ncid) 
call nc_write(fnm,"runoff_ann_ice", sngl(vars%ann_runoff_ice        ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="annual runoff, for ice cells only",grid_mapping="polar_stereographic",units="kg/m2/yr",ncid=ncid) 
call nc_write(fnm,"evp_ann", sngl(vars%ann_evp        ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="annual sublimation",grid_mapping="polar_stereographic",units="kg/m2/yr",ncid=ncid) 
call nc_write(fnm,"evp_ann_ice", sngl(vars%ann_evp_ice        ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="annual sublimation, for ice cells only",grid_mapping="polar_stereographic",units="kg/m2/yr",ncid=ncid) 
call nc_write(fnm,"refreezing_ann", sngl(vars%ann_refreezing        ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="refreezing",grid_mapping="polar_stereographic",units="kg/m2/yr",ncid=ncid) 
call nc_write(fnm,"refreezing_ann_ice", sngl(vars%ann_refreezing_ice        ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1],long_name="refreezing, for ice cells only",grid_mapping="polar_stereographic",units="kg/m2/yr",ncid=ncid) 
    endif

   return

  end subroutine smb_nc_write


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s m b _ d a i l y_ n c
  ! Purpose  :  Initialize sea ice netcdf output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_daily_nc(fnm,grid)

    implicit none

    character (len=*) :: fnm
    type(grid_class) :: grid
    integer :: ncid
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,dim_time, x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE., ncid=ncid)
    call nc_write_dim(fnm, dim_day, x=1._wp, dx=1._wp, nx=nday_year, axis="e", units="days", ncid=ncid)
    call nc_write_dim(fnm, dim_depth, x=z, axis="z", units="m", ncid=ncid)
    call nc_write_dim(fnm, dim_y, x=grid%G%y0, dx=grid%G%dy, nx=grid%G%ny, axis="y", units="km", ncid=ncid)
    call nc_write_dim(fnm, dim_x, x=grid%G%x0, dx=grid%G%dx, nx=grid%G%nx, axis="x", units="km", ncid=ncid)

    call nc_write(fnm,"phi",sngl(grid%lat), dims=[dim_x,dim_y], start=[1,1], &
    count=[grid%G%nx,grid%G%ny],ncid=ncid)
    call nc_write(fnm,"lambda",sngl(grid%lon), dims=[dim_x,dim_y], start=[1,1], &
    count=[grid%G%nx,grid%G%ny], ncid=ncid)
    call nc_close(ncid)

   return

  end subroutine smb_daily_nc


  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s m b _ d a i l y _ n c _ w r i t e
  ! Purpose  :  Output of sea ice netcdf
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_daily_nc_write(fnm,ncid,vars,grid,ndat,nout)

    implicit none

    type(s_out) :: vars
    type(grid_class) :: grid

    character (len=*) :: fnm
    integer :: ndat, nout, ncid


call nc_write(fnm,"rain", sngl(vars%rain            ), dims=[dim_x,dim_y,dim_day,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="rainfall rate",grid_mapping="polar_stereographic",units="kg/m2/day",ncid=ncid) 
call nc_write(fnm,"snow", sngl(vars%snow            ), dims=[dim_x,dim_y,dim_day,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="snowfall rate",grid_mapping="polar_stereographic",units="kg/m2/day",ncid=ncid) 
call nc_write(fnm,"t2m", sngl(vars%t2m           ), dims=[dim_x,dim_y,dim_day,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="surface air temperature",grid_mapping="polar_stereographic",units="K",ncid=ncid) 
call nc_write(fnm,"t_skin", sngl(vars%t_skin          ), dims=[dim_x,dim_y,dim_day,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="skin temperature",grid_mapping="polar_stereographic",units="K",ncid=ncid) 
call nc_write(fnm,"t_skin_amp", sngl(vars%t_skin_amp          ), dims=[dim_x,dim_y,dim_day,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="amplitude of diurnal cycle of skin temperature",grid_mapping="polar_stereographic",units="K",ncid=ncid) 
call nc_write(fnm,"t_prof", sngl(vars%t_prof          ), dims=[dim_x,dim_y,dim_depth,dim_day,dim_time],start=[1,1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,nl+1,1,1],long_name="temperature profile in snow+ice/soil",grid_mapping="polar_stereographic",units="K",ncid=ncid) 
call nc_write(fnm,"mask_snow", sngl(vars%mask_snow       ), dims=[dim_x,dim_y,dim_day,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="snow mask",grid_mapping="polar_stereographic",units="/",ncid=ncid) 
call nc_write(fnm,"h_snow", sngl(vars%h_snow          ), dims=[dim_x,dim_y,dim_day,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="snow thickness",grid_mapping="polar_stereographic",units="m",ncid=ncid) 
call nc_write(fnm,"w_snow", sngl(vars%w_snow          ), dims=[dim_x,dim_y,dim_day,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="snow water equivalent",grid_mapping="polar_stereographic",units="kg/m2",ncid=ncid) 
call nc_write(fnm,"w_snow_max", sngl(vars%w_snow_max      ), dims=[dim_x,dim_y,dim_day,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="seasonal max of snow water equivalent",grid_mapping="polar_stereographic",units="kg/m2",ncid=ncid) 
call nc_write(fnm,"snowmelt", sngl(vars%snowmelt        ), dims=[dim_x,dim_y,dim_day,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="snow melt rate",grid_mapping="polar_stereographic",units="kg/m2/day",ncid=ncid) 
call nc_write(fnm,"icemelt", sngl(vars%icemelt         ), dims=[dim_x,dim_y,dim_day,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="ice melt rate",grid_mapping="polar_stereographic",units="kg/m2/day",ncid=ncid) 
call nc_write(fnm,"runoff", sngl(vars%runoff         ), dims=[dim_x,dim_y,dim_day,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="runoff",grid_mapping="polar_stereographic",units="kg/m2/day",ncid=ncid) 
call nc_write(fnm,"refreezing", sngl(vars%refreezing        ), dims=[dim_x,dim_y,dim_day,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="refreezing",grid_mapping="polar_stereographic",units="kg/m2/day",ncid=ncid) 
call nc_write(fnm,"flx_g", sngl(vars%flx_g           ), dims=[dim_x,dim_y,dim_day,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="ground heat flux",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"dflxg_dT", sngl(vars%dflxg_dT        ), dims=[dim_x,dim_y,dim_day,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="derivative of ground heat flux wrt top soil temperature",grid_mapping="polar_stereographic",units="W/m2/K",ncid=ncid) 
call nc_write(fnm,"flx_melt", sngl(vars%flx_melt        ), dims=[dim_x,dim_y,dim_day,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="heat flux going into the melt of snow",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"flx_sh", sngl(vars%flx_sh          ), dims=[dim_x,dim_y,dim_day,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="sensible heat flux",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"flx_lwu", sngl(vars%flx_lwu         ), dims=[dim_x,dim_y,dim_day,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="upward longwave radiation",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"flx_lh", sngl(vars%flx_lh          ), dims=[dim_x,dim_y,dim_day,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="latent_heat_flux",grid_mapping="polar_stereographic",units="W/m2",ncid=ncid) 
call nc_write(fnm,"dust_con", sngl(vars%dust_con      ), dims=[dim_x,dim_y,dim_day,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1],long_name="dust concentration in snow",grid_mapping="polar_stereographic",units="ppm",ncid=ncid) 


   return

  end subroutine smb_daily_nc_write


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s m b _ a v e
  ! Purpose  :  Average (or sum) the sea ice fields
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_ave(d,ave)

    implicit none

    type(s_out) :: d(:), ave

    integer :: k, n
    real(wp) :: div

    n = size(d)
    div = dble(n)

    ! Set all values to zero
    ave%t2m_i            = 0._wp    
    ave%t2m_bias_i       = 0._wp    
    ave%tam_i            = 0._wp    
    ave%gam_i            = 0._wp    
    ave%ram_i            = 0._wp    
    ave%prc_i            = 0._wp    
    ave%prc_bias_i       = 0._wp    
    ave%u700_i           = 0._wp    
    ave%v700_i           = 0._wp    
    ave%wind_i           = 0._wp    
    ave%cld_i            = 0._wp    
    ave%dust_i           = 0._wp    
    ave%swdown_i      = 0._wp    
    ave%swd_toa_i     = 0._wp    
    ave%swd_sur_vis_dir_i = 0._wp    
    ave%swd_sur_nir_dir_i = 0._wp    
    ave%swd_sur_vis_dif_i = 0._wp    
    ave%swd_sur_nir_dif_i = 0._wp    
    ave%dswd_dalb_vis_dir_i = 0._wp 
    ave%dswd_dalb_nir_dir_i = 0._wp 
    ave%dswd_dalb_vis_dif_i = 0._wp 
    ave%dswd_dalb_nir_dif_i = 0._wp 
    ave%dswd_dz_nir_dir_i   = 0._wp 
    ave%dswd_dz_nir_dif_i   = 0._wp 
    ave%coszm_i          = 0._wp    
    ave%lwdown_i         = 0._wp    
    ave%gam_lw_i         = 0._wp    
    ave%tam              = 0._wp    
    ave%t2m              = 0._wp    
    ave%t_skin           = 0._wp    
    ave%tskin_tam        = 0._wp    
    ave%t_skin_amp       = 0._wp    
    ave%t_prof           = 0._wp    
    ave%q2m             = 0._wp    
    ave%pressure         = 0._wp    
    ave%mask_snow        = 0._wp    
    ave%f_snow           = 0._wp    
    ave%h_snow           = 0._wp    
    ave%w_snow           = 0._wp    
    ave%w_snow_max       = 0._wp    
    ave%snowmelt         = 0._wp    
    ave%icemelt          = 0._wp    
    ave%runoff           = 0._wp    
    ave%refreezing       = 0._wp    
    ave%f_rfz_to_snow       = 0._wp    
    ave%albedo           = 0._wp    
    ave%alb_vis_dir      = 0._wp    
    ave%alb_nir_dir      = 0._wp    
    ave%alb_vis_dif      = 0._wp    
    ave%alb_nir_dif      = 0._wp    
    ave%alb_snow_vis_dir = 0._wp    
    ave%alb_snow_nir_dir = 0._wp    
    ave%alb_snow_vis_dif = 0._wp    
    ave%alb_snow_nir_dif = 0._wp    
    ave%snow_grain       = 0._wp    
    ave%dust_con         = 0._wp    
    ave%cld              = 0._wp    
    ave%swnet            = 0._wp    
    ave%swnet_min        = 0._wp    
    ave%swdown           = 0._wp    
    ave%lwdown           = 0._wp    
    ave%r_a              = 0._wp    
    ave%flx_g            = 0._wp    
    ave%dflxg_dT         = 0._wp    
    ave%flx_melt         = 0._wp    
    ave%flx_sh           = 0._wp    
    ave%flx_lwu          = 0._wp    
    ave%flx_lh           = 0._wp    
    ave%evp              = 0._wp    
    ave%prc              = 0._wp    
    ave%f_wind           = 0._wp    
    ave%rain             = 0._wp    
    ave%snow             = 0._wp    
    ave%u700             = 0._wp    
    ave%v700             = 0._wp    
    ave%wind             = 0._wp    


    ! Loop over the time indices to sum up and average (if necessary)
    do k = 1, n
      ave%t2m_i            = ave%t2m_i            + d(k)%t2m_i            / div       
      ave%t2m_bias_i       = ave%t2m_bias_i       + d(k)%t2m_bias_i       / div       
      ave%tam_i            = ave%tam_i            + d(k)%tam_i            / div       
      ave%gam_i            = ave%gam_i            + d(k)%gam_i            / div       
      ave%ram_i            = ave%ram_i            + d(k)%ram_i            / div       
      ave%prc_i            = ave%prc_i            + d(k)%prc_i            / div       
      ave%prc_bias_i       = ave%prc_bias_i       + d(k)%prc_bias_i       / div       
      ave%u700_i           = ave%u700_i           + d(k)%u700_i           / div       
      ave%v700_i           = ave%v700_i           + d(k)%v700_i           / div       
      ave%wind_i           = ave%wind_i           + d(k)%wind_i           / div       
      ave%cld_i            = ave%cld_i            + d(k)%cld_i            / div       
      ave%dust_i           = ave%dust_i           + d(k)%dust_i           / div       
      ave%swdown_i      = ave%swdown_i      + d(k)%swdown_i     / div       
      ave%swd_toa_i     = ave%swd_toa_i     + d(k)%swd_toa_i     / div       
      ave%swd_sur_vis_dir_i = ave%swd_sur_vis_dir_i + d(k)%swd_sur_vis_dir_i / div       
      ave%swd_sur_nir_dir_i = ave%swd_sur_nir_dir_i + d(k)%swd_sur_nir_dir_i / div       
      ave%swd_sur_vis_dif_i = ave%swd_sur_vis_dif_i + d(k)%swd_sur_vis_dif_i / div       
      ave%swd_sur_nir_dif_i = ave%swd_sur_nir_dif_i + d(k)%swd_sur_nir_dif_i / div       
      ave%dswd_dalb_vis_dir_i = ave%dswd_dalb_vis_dir_i + d(k)%dswd_dalb_vis_dir_i / div   
      ave%dswd_dalb_nir_dir_i = ave%dswd_dalb_nir_dir_i + d(k)%dswd_dalb_nir_dir_i / div
      ave%dswd_dalb_vis_dif_i = ave%dswd_dalb_vis_dif_i + d(k)%dswd_dalb_vis_dif_i / div
      ave%dswd_dalb_nir_dif_i = ave%dswd_dalb_nir_dif_i + d(k)%dswd_dalb_nir_dif_i / div
      ave%dswd_dz_nir_dir_i   = ave%dswd_dz_nir_dir_i   + d(k)%dswd_dz_nir_dir_i   / div
      ave%dswd_dz_nir_dif_i   = ave%dswd_dz_nir_dif_i   + d(k)%dswd_dz_nir_dif_i   / div
      ave%coszm_i          = ave%coszm_i          + d(k)%coszm_i          / div       
      ave%lwdown_i         = ave%lwdown_i         + d(k)%lwdown_i         / div       
      ave%gam_lw_i         = ave%gam_lw_i         + d(k)%gam_lw_i         / div       
      ave%tam              = ave%tam              + d(k)%tam              / div       
      ave%t2m              = ave%t2m              + d(k)%t2m              / div       
      ave%t_skin           = ave%t_skin           + d(k)%t_skin           / div       
      ave%tskin_tam        = ave%tskin_tam        + d(k)%tskin_tam        / div       
      ave%t_skin_amp       = ave%t_skin_amp       + d(k)%t_skin_amp       / div       
      ave%t_prof           = ave%t_prof           + d(k)%t_prof           / div       
      ave%q2m             = ave%q2m             + d(k)%q2m             / div       
      ave%pressure         = ave%pressure         + d(k)%pressure         / div       
      ave%mask_snow        = ave%mask_snow        + d(k)%mask_snow        / div       
      ave%f_snow           = ave%f_snow           + d(k)%f_snow           / div       
      ave%h_snow           = ave%h_snow           + d(k)%h_snow           / div       
      ave%w_snow           = ave%w_snow           + d(k)%w_snow           / div       
      ave%w_snow_max       = ave%w_snow_max       + d(k)%w_snow_max       / div       
      ave%snowmelt         = ave%snowmelt         + d(k)%snowmelt         / div       
      ave%icemelt          = ave%icemelt          + d(k)%icemelt          / div       
      ave%runoff           = ave%runoff           + d(k)%runoff           / div       
      ave%refreezing       = ave%refreezing       + d(k)%refreezing       / div       
      ave%f_rfz_to_snow       = ave%f_rfz_to_snow       + d(k)%f_rfz_to_snow       / div       
      ave%albedo           = ave%albedo           + d(k)%albedo           / div       
      ave%alb_vis_dir      = ave%alb_vis_dir      + d(k)%alb_vis_dir      / div       
      ave%alb_nir_dir      = ave%alb_nir_dir      + d(k)%alb_nir_dir      / div       
      ave%alb_vis_dif      = ave%alb_vis_dif      + d(k)%alb_vis_dif      / div       
      ave%alb_nir_dif      = ave%alb_nir_dif      + d(k)%alb_nir_dif      / div       
      ave%alb_snow_vis_dir = ave%alb_snow_vis_dir + d(k)%alb_snow_vis_dir / div       
      ave%alb_snow_nir_dir = ave%alb_snow_nir_dir + d(k)%alb_snow_nir_dir / div       
      ave%alb_snow_vis_dif = ave%alb_snow_vis_dif + d(k)%alb_snow_vis_dif / div       
      ave%alb_snow_nir_dif = ave%alb_snow_nir_dif + d(k)%alb_snow_nir_dif / div       
      ave%snow_grain       = ave%snow_grain       + d(k)%snow_grain       / div       
      ave%dust_con         = ave%dust_con         + d(k)%dust_con         / div       
      ave%cld              = ave%cld              + d(k)%cld              / div       
      ave%swnet            = ave%swnet            + d(k)%swnet            / div       
      ave%swnet_min        = ave%swnet_min        + d(k)%swnet_min        / div       
      ave%swdown           = ave%swdown           + d(k)%swdown           / div       
      ave%lwdown           = ave%lwdown           + d(k)%lwdown           / div       
      ave%r_a              = ave%r_a              + d(k)%r_a              / div       
      ave%flx_g            = ave%flx_g            + d(k)%flx_g            / div       
      ave%dflxg_dT         = ave%dflxg_dT         + d(k)%dflxg_dT         / div       
      ave%flx_melt         = ave%flx_melt         + d(k)%flx_melt         / div       
      ave%flx_sh           = ave%flx_sh           + d(k)%flx_sh           / div       
      ave%flx_lwu          = ave%flx_lwu          + d(k)%flx_lwu          / div       
      ave%flx_lh           = ave%flx_lh           + d(k)%flx_lh           / div       
      ave%evp              = ave%evp              + d(k)%evp              / div       
      ave%prc              = ave%prc              + d(k)%prc              / div       
      ave%f_wind           = ave%f_wind           + d(k)%f_wind           / div       
      ave%rain             = ave%rain             + d(k)%rain             / div       
      ave%snow             = ave%snow             + d(k)%snow             / div       
      ave%u700             = ave%u700             + d(k)%u700             / div       
      ave%v700             = ave%v700             + d(k)%v700             / div       
      ave%wind             = ave%wind             + d(k)%wind             / div       

    end do

   return

  end subroutine smb_ave


end module smb_out
