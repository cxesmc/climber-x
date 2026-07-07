!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c a r b o n _ i n v e n t o r y _ m o d
!
!  Purpose : carbon pools inventory
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
module lndvc_carbon_inventory_mod

  use precision, only : wp
  use lnd_grid, only : npft, nl, nlc, dz_c

  implicit none

  private
  public :: carbon_inventory

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  c a r b o n _ i n v e n t o r y
  !   Purpose    :  carbon inventory
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine carbon_inventory(frac_surf, f_veg, f_peat, f_ice_grd, f_shelf, f_lake, &
      veg_c, veg_c13, veg_c14, &
      litter_c, fast_c, slow_c, litter_c13, fast_c13, slow_c13, litter_c14, fast_c14, slow_c14, &
      litter_c_peat, acro_c, cato_c, litter_c13_peat, acro_c13, cato_c13, litter_c14_peat, acro_c14, cato_c14, &
      litter_c_ice,   fast_c_ice,   slow_c_ice,   litter_c13_ice,   fast_c13_ice,   slow_c13_ice,   litter_c14_ice,   fast_c14_ice,   slow_c14_ice, &
      litter_c_shelf, fast_c_shelf, slow_c_shelf, litter_c13_shelf, fast_c13_shelf, slow_c13_shelf, litter_c14_shelf, fast_c14_shelf, slow_c14_shelf, &
      litter_c_lake,  fast_c_lake,  slow_c_lake,  litter_c13_lake,  fast_c13_lake,  slow_c13_lake,  litter_c14_lake,  fast_c14_lake,  slow_c14_lake, &
      landc, landc13, landc14, burc, burc13, burc14)

    implicit none

    real(wp), dimension(:), intent(in) :: frac_surf
    real(wp), intent(in) :: f_veg, f_ice_grd, f_shelf, f_lake, f_peat
    real(wp), dimension(:), intent(in) :: veg_c, veg_c13, veg_c14
    real(wp), dimension(:), intent(in) :: litter_c, fast_c, slow_c, litter_c13, fast_c13, slow_c13, litter_c14, fast_c14, slow_c14
    real(wp), dimension(:), intent(in) :: litter_c_shelf, fast_c_shelf, slow_c_shelf, litter_c13_shelf, fast_c13_shelf, slow_c13_shelf, litter_c14_shelf, fast_c14_shelf, slow_c14_shelf
    real(wp), dimension(:), intent(in) :: litter_c_lake, fast_c_lake, slow_c_lake, litter_c13_lake, fast_c13_lake, slow_c13_lake, litter_c14_lake, fast_c14_lake, slow_c14_lake
    real(wp), dimension(:), intent(in) :: litter_c_ice, fast_c_ice, slow_c_ice, litter_c13_ice, fast_c13_ice, slow_c13_ice, litter_c14_ice, fast_c14_ice, slow_c14_ice
    real(wp), dimension(:), intent(in) :: cato_c, cato_c13, cato_c14
    real(wp), intent(in) :: litter_c_peat, acro_c, litter_c13_peat, acro_c13, litter_c14_peat, acro_c14
    real(wp), intent(out) :: landc, landc13, landc14
    real(wp), intent(out) :: burc, burc13, burc14

    integer :: n, k
    real(wp) :: vegc, soilc, minc, peatc, icec, shelfc, lakec
    real(wp) :: vegc13, soilc13, minc13, peatc13, icec13, shelfc13, lakec13
    real(wp) :: vegc14, soilc14, minc14, peatc14, icec14, shelfc14, lakec14


    ! vegetation carbon
    vegc = 0._wp
    do k = 1, npft
      vegc = vegc + veg_c(k)*frac_surf(k) ! kgC/m2
    enddo
    ! soil carbon
    soilc = 0._wp
    burc = 0._wp
    do n=1,nlc
      minc  = (litter_c(n) + fast_c(n) + slow_c(n)) * dz_c(n) * (f_veg-f_peat) ! kgC/m2
      if (n.eq.1) then
        peatc = (litter_c_peat/dz_c(n) + acro_c/dz_c(n) + cato_c(n)) * dz_c(n) * f_peat  ! kgC/m2
      else
        peatc = cato_c(n) * dz_c(n) * f_peat ! kgC/m2
      endif
      icec   = (litter_c_ice(n) + fast_c_ice(n) + slow_c_ice(n)) * dz_c(n) * f_ice_grd  ! kgC/m2
      shelfc = (litter_c_shelf(n) + fast_c_shelf(n) + slow_c_shelf(n)) * dz_c(n) * f_shelf  ! kgC/m2
      lakec = (litter_c_lake(n) + fast_c_lake(n) + slow_c_lake(n)) * dz_c(n) * f_lake  ! kgC/m2
      ! depth integrated soil carbon
      if (n.le.nl) then
        soilc = soilc + minc + peatc + icec + shelfc + lakec ! kgC/m2
      else if (n.eq.nlc) then
        burc = burc + minc + peatc + icec + shelfc + lakec ! kgC/m2
      endif
    enddo
    ! total land carbon in grid cell
    landc = soilc + vegc ! kgC/m2

    ! vegetation carbon 13
    vegc13 = 0._wp
    do k = 1, npft
      vegc13 = vegc13 + veg_c13(k)*frac_surf(k) ! kgC/m2
    enddo
    ! soil carbon
    soilc13 = 0._wp
    burc13 = 0._wp
    do n=1,nlc
      minc13  = (litter_c13(n) + fast_c13(n) + slow_c13(n)) * dz_c(n) * (f_veg-f_peat)  ! kgC/m2
      if (n.eq.1) then
        peatc13 = (litter_c13_peat/dz_c(n) + acro_c13/dz_c(n) + cato_c13(n)) * dz_c(n) * f_peat  ! kgC/m2
      else
        peatc13 = cato_c13(n) * dz_c(n) * f_peat  ! kgC/m2 
      endif
      icec13   = (litter_c13_ice(n) + fast_c13_ice(n) + slow_c13_ice(n)) * dz_c(n) * f_ice_grd  ! kgC/m2
      shelfc13 = (litter_c13_shelf(n) + fast_c13_shelf(n) + slow_c13_shelf(n)) * dz_c(n) * f_shelf  ! kgC/m2
      lakec13 = (litter_c13_lake(n) + fast_c13_lake(n) + slow_c13_lake(n)) * dz_c(n) * f_lake  ! kgC/m2
      ! depth integrated soil carbon
      if (n.le.nl) then
        soilc13  = soilc13 + minc13 + peatc13 + icec13 + shelfc13 + lakec13 ! kgC/m2
      else if (n.eq.nlc) then
        burc13  = burc13 + minc13 + peatc13 + icec13 + shelfc13 + lakec13 ! kgC/m2
      endif
    enddo
    ! total land carbon in grid cell
    landc13 = soilc13 + vegc13 ! kgC/m2

    ! vegetation carbon
    vegc14 = 0._wp
    do k = 1, npft
      vegc14 = vegc14 + veg_c14(k)*frac_surf(k) ! kgC/m2
    enddo
    ! soil carbon
    soilc14 = 0._wp
    burc14 = 0._wp
    do n=1,nlc
      minc14  = (litter_c14(n) + fast_c14(n) + slow_c14(n)) * dz_c(n) * (f_veg-f_peat)  ! kgC/m2
      if (n.eq.1) then
        peatc14 = (litter_c14_peat/dz_c(n) + acro_c14/dz_c(n) + cato_c14(n)) * dz_c(n) * f_peat  ! kgC/m2
      else
        peatc14 = cato_c14(n) * dz_c(n) * f_peat  ! kgC/m2 
      endif
      icec14   = (litter_c14_ice(n) + fast_c14_ice(n) + slow_c14_ice(n)) * dz_c(n) * f_ice_grd  ! kgC/m2
      shelfc14 = (litter_c14_shelf(n) + fast_c14_shelf(n) + slow_c14_shelf(n)) * dz_c(n) * f_shelf  ! kgC/m2
      lakec14 = (litter_c14_lake(n) + fast_c14_lake(n) + slow_c14_lake(n)) * dz_c(n) * f_lake  ! kgC/m2
      ! depth integrated soil carbon
      if (n.le.nl) then
        soilc14  = soilc14 + minc14 + peatc14 + icec14 + shelfc14 + lakec14 ! kgC/m2
      else if (n.eq.nlc) then
        burc14  = burc14 + minc14 + peatc14 + icec14 + shelfc14 + lakec14 ! kgC/m2
      endif
      
    enddo
    ! total land carbon in grid cell
    landc14 = soilc14 + vegc14 ! kgC/m2


    return

  end subroutine carbon_inventory

end module lndvc_carbon_inventory_mod
