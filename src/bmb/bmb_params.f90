!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : b m b _ p a r a m s
!
!  Purpose : BMB model parameters
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
module bmb_params

  use precision, only : wp
  use nml
  use timer, only : dt_bmb, sec_year
  use control, only : out_dir

  implicit none

  real(wp) :: dt, rdt

  integer :: i_bmb
  real(wp) :: bmb_const
  real(wp) :: k_1
  real(wp) :: k_2

  logical :: l_fix_depth
  real(wp) :: fix_depth
  real(wp) :: depth_disc

  integer :: i_bmb_lake

  logical :: l_depth_scale
  real(wp) :: zl_ref

  logical :: l_monthly_output

contains

    subroutine bmb_params_init

    implicit none


    ! time step
    dt = dt_bmb
    rdt = 1._wp/dt

    ! read bmb parameter file
    call bmb_par_load(trim(out_dir)//"/bmb_par.nml")


    return

    end subroutine bmb_params_init


subroutine bmb_par_load(filename)

    implicit none

    character (len=*) :: filename


    ! Read parameters from file
    write(*,*) "bmb parameters ==========="
    call nml_read(filename,"bmb_par","i_bmb",i_bmb)
    call nml_read(filename,"bmb_par","bmb_const",bmb_const)
    bmb_const = bmb_const/sec_year  ! kg/m2/a -> kg/m2/s
    call nml_read(filename,"bmb_par","k_1",k_1)
    call nml_read(filename,"bmb_par","k_2",k_2)
    call nml_read(filename,"bmb_par","l_fix_depth",l_fix_depth)
    call nml_read(filename,"bmb_par","fix_depth",fix_depth)
    call nml_read(filename,"bmb_par","depth_disc",depth_disc)
    call nml_read(filename,"bmb_par","i_bmb_lake",i_bmb_lake)
    call nml_read(filename,"bmb_par","l_depth_scale",l_depth_scale)
    call nml_read(filename,"bmb_par","zl_ref",zl_ref)
    call nml_read(filename,"bmb_par","l_monthly_output",l_monthly_output)

   return

end subroutine bmb_par_load


end module bmb_params
