!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : d i m _ n a m e
!
!  Purpose : define dimensions for netcdf output 
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
module dim_name
  ! dimensions to build arrays with same-length characters (gfortran)

  implicit none

  integer, parameter :: len = 10

  character(len=len), parameter :: dim_lon = "lon"
  character(len=len), parameter :: dim_lonu = "lonu"
  character(len=len), parameter :: dim_lat = "lat"
  character(len=len), parameter :: dim_latv = "latv"
  character(len=len), parameter :: dim_latv1 = "latv1"
  character(len=len), parameter :: dim_lev = "lev"
  character(len=len), parameter :: dim_levw = "levw"
  character(len=len), parameter :: dim_x = "x"
  character(len=len), parameter :: dim_y = "y"
  character(len=len), parameter :: dim_day = "day"
  character(len=len), parameter :: dim_month = "month"
  character(len=len), parameter :: dim_time = "time"
  character(len=len), parameter :: dim_lon1 = "lon1"
  character(len=len), parameter :: dim_lat1 = "lat1"
  character(len=len), parameter :: dim_type = "type"
  character(len=len), parameter :: dim_depth = "depth"
  character(len=len), parameter :: dim_depth0 = "depth0"
  character(len=len), parameter :: dim_depthl = "depthl"
  character(len=len), parameter :: dim_depth0l = "depth0l"
  character(len=len), parameter :: dim_depth1 = "depth1"
  character(len=len), parameter :: dim_depth_sed = "depth_sed"
  character(len=len), parameter :: dim_depth_eu = "depth_eu"
  character(len=len), parameter :: dim_isles = "isles"
  character(len=len), parameter :: dim_kc = "kc"
  character(len=len), parameter :: dim_kt = "kt"
  character(len=len), parameter :: dim_kr = "kr"
  character(len=len), parameter :: dim_npft = "npft"
  character(len=len), parameter :: dim_nsurf = "nsurf"
  character(len=len), parameter :: dim_nsoil = "nsoil"
  character(len=len), parameter :: dim_ncarb = "ncarb"
  character(len=len), parameter :: dim_nlit  = "nlit"
  character(len=len), parameter :: dim_dir = "dir"
  character(len=len), parameter :: dim_nocetra = "nocetra" 
  character(len=len), parameter :: dim_npowtra = "npowtra"
  character(len=len), parameter :: dim_nsedtra = "nsedtra"

end module


