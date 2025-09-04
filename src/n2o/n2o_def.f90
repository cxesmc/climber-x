!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : n 2 o _ d e f
!
!  Purpose : definition of N2O model class 
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Andrey Ganopolski and Matteo Willeit
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
module n2o_def

    use precision, only : wp

    implicit none
    
    type n2o_class   
      real(wp) :: n2o             !! atmospheric n2o concentration [ppb]
      real(wp) :: n2om            !! atmospheric n2o mass [kgN2O]
      real(wp) :: n2o_slow        !! slow (relaxed) atmospheric ch4 concentration [ppb]
      real(wp) :: dn2oocn_dt      !! ocn methane flux [kgN2O]
      real(wp) :: dn2olnd_dt      !! lnd methane flux [kgN2O]
      real(wp) :: dn2oemis_dt     !! methane emissions [kgN2O]
      real(wp) :: dn2oox_dt       !! methane oxidation [kgN2O]
      real(wp) :: tau             !! atmospheric methane lifetime [years]
    end type


    private
    public :: n2o_class

end module
