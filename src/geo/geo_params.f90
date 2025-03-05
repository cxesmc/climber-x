!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : g e o _ p a r a m s
!
!  Purpose : geography model parameters
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
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
module geo_params

  use precision, only : wp
  use nml
  use control, only : out_dir

  implicit none

  integer :: i_geo

  logical :: l_visc_3d
  character (len=256) :: visc_1d_file
  character (len=256) :: visc_3d_file
  character (len=256) :: vilma_grid_file

  character (len=256) :: geo_ref_file

  character (len=256) :: z_bed_rel_file

  character (len=256) :: z_bed_1min_file

  real(wp) :: f_crit
  real(wp) :: f_crit_eq
  logical :: l_ocn_below_shelf
  logical :: l_connect_ocn
  logical :: l_close_panama
  logical :: l_close_bering
  logical :: l_close_hudson
  logical :: l_close_baltic
  integer :: i_fix_cell_grl

  real(wp) :: h_ice_min

  integer :: i_z_min_max
  real(wp) :: sigma_filter
  real(wp) :: z_ocn_max_quant

  integer :: i_lakes
  real(wp) :: lake_area_crit 
  real(wp) :: lake_sea_z_crit 

  real(wp) :: lon_ocn_origin
  real(wp) :: lat_ocn_origin

  integer :: n_coast_cells

  real(wp) :: gia_time_lag

  integer :: i_q_geo
  real(wp) :: q_geo_const
  character (len=256) :: q_geo_file

  character (len=256) :: sed_file

  logical :: l_write_timer

  logical :: l_output_hires

contains

    subroutine geo_params_init

    implicit none

    call geo_par_load(trim(out_dir)//"/geo_par.nml")


    return

    end subroutine geo_params_init


subroutine geo_par_load(filename)

    implicit none

    character (len=*) :: filename


    ! Read parameters from file
    write(*,*) "geo parameters ==========="
    call nml_read(filename,"geo_par","i_geo",i_geo)
    call nml_read(filename,"geo_par","l_visc_3d",l_visc_3d)
    call nml_read(filename,"geo_par","visc_1d_file",visc_1d_file)
    call nml_read(filename,"geo_par","visc_3d_file",visc_3d_file)
    call nml_read(filename,"geo_par","vilma_grid_file",vilma_grid_file)
    call nml_read(filename,"geo_par","f_crit",f_crit)
    call nml_read(filename,"geo_par","f_crit_eq",f_crit_eq)
    call nml_read(filename,"geo_par","l_ocn_below_shelf",l_ocn_below_shelf)
    call nml_read(filename,"geo_par","l_connect_ocn",l_connect_ocn)
    call nml_read(filename,"geo_par","l_close_panama",l_close_panama)
    call nml_read(filename,"geo_par","l_close_bering",l_close_bering)
    call nml_read(filename,"geo_par","l_close_hudson",l_close_hudson)
    call nml_read(filename,"geo_par","l_close_baltic",l_close_baltic)
    call nml_read(filename,"geo_par","i_fix_cell_grl",i_fix_cell_grl)
    call nml_read(filename,"geo_par","h_ice_min",h_ice_min)
    call nml_read(filename,"geo_par","i_z_min_max",i_z_min_max)
    call nml_read(filename,"geo_par","sigma_filter",sigma_filter)
    call nml_read(filename,"geo_par","z_ocn_max_quant",z_ocn_max_quant)
    call nml_read(filename,"geo_par","geo_ref_file",geo_ref_file)
    call nml_read(filename,"geo_par","z_bed_rel_file",z_bed_rel_file)
    call nml_read(filename,"geo_par","z_bed_1min_file",z_bed_1min_file)
    call nml_read(filename,"geo_par","i_lakes",i_lakes)
    call nml_read(filename,"geo_par","lake_area_crit",lake_area_crit)
    call nml_read(filename,"geo_par","lake_sea_z_crit",lake_sea_z_crit)
    call nml_read(filename,"geo_par","lon_ocn_origin",lon_ocn_origin)
    call nml_read(filename,"geo_par","lat_ocn_origin",lat_ocn_origin)
    call nml_read(filename,"geo_par","n_coast_cells",n_coast_cells)
    call nml_read(filename,"geo_par","gia_time_lag",gia_time_lag)
    call nml_read(filename,"geo_par","i_q_geo",i_q_geo)
    call nml_read(filename,"geo_par","q_geo_const",q_geo_const)
    call nml_read(filename,"geo_par","q_geo_file",q_geo_file)
    call nml_read(filename,"geo_par","sed_file",sed_file)
    call nml_read(filename,"geo_par","l_write_timer",l_write_timer)
    call nml_read(filename,"geo_par","l_output_hires",l_output_hires)

 
   return

end subroutine geo_par_load


end module geo_params
