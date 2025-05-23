!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c h 4 _ d e f
!
!  Purpose : definition of CH4 model class 
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
module ch4_def

    use precision, only : wp

    implicit none
    
    type ch4_class   
      real(wp) :: ch4             !! atmospheric ch4 concentration [ppb]
      real(wp) :: ch4m            !! atmospheric ch4 mass [kgCH4]
      real(wp) :: ch4_slow        !! slow (relaxed) atmospheric ch4 concentration [ppb]
      real(wp) :: dch4ocn_dt      !! ocn methane flux [kgCH4]
      real(wp) :: dch4lnd_dt      !! lnd methane flux [kgCH4]
      real(wp) :: dch4emis_dt     !! methane emissions [kgCH4]
      real(wp) :: f_ch4emis_agro  !! fraction of anthropogenic methane emissions originating from agriculture
      real(wp) :: dch4ox_dt       !! methane oxidation [kgCH4]
      real(wp) :: tau             !! atmospheric methane lifetime [years]
      real(wp) :: t2m_glob_ann    !! global mean annual near-surface air temperature [K]
      real(wp) :: demis_other_dt  !! other methane flux [kg CH4]
    end type


    private
    public :: ch4_class

end module
