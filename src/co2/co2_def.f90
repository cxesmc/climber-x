!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c o 2 _ d e f
!
!  Purpose : definition of CO2 model class 
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
module co2_def

    use precision, only : wp

    implicit none
    
    type co2_class   
      real(wp) :: co2           !! atmospheric CO2 concentration [ppm]
      real(wp) :: Catm          !! atmospheric carbon content [GtC]
      real(wp) :: C13atm        !! atmospheric carbon 13 content [GtC]
      real(wp) :: C14atm        !! atmospheric carbon 14 content [GtC]
      real(wp) :: c13_c12       !! carbon 13 to carbon 12 ratio []
      real(wp) :: c14_c         !! carbon 14 to carbon 12 ratio []
      real(wp) :: dCocn_dt      !! atm-ocn carbon flux [kgC/yr]
      real(wp) :: dC13ocn_dt    !! atm-ocn carbon 13 flux [kgC/yr]
      real(wp) :: dC14ocn_dt    !! atm-ocn carbon 14 flux [kgC/yr]
      real(wp) :: dClnd_dt      !! atm-lnd carbon flux [kgC/yr]
      real(wp) :: dC13lnd_dt    !! atm-lnd carbon 13 flux [kgC/yr]
      real(wp) :: dC14lnd_dt    !! atm-lnd carbon 14 flux [kgC/yr]
      real(wp) :: dCweath_dt    !! weathering carbon flux [kgC/yr]
      real(wp) :: weath_sil_avg !! initial average silicate weathering rate [kgC/yr]
      real(wp) :: dC13weath_dt  !! weathering carbon 13 flux [kgC/yr]
      real(wp) :: dC14weath_dt  !! weathering carbon 14 flux [kgC/yr]
      real(wp) :: dCocn_dt_avg     !! average atm-ocn carbon flux [kgC/yr]
      real(wp) :: dClnd_dt_avg     !! average atm-lnd carbon flux [kgC/yr]
      real(wp) :: dCvolc_dt     !! volcanic degassing carbon [kgC]/yr
      real(wp) :: dCvolc_dt_eq  !! volcanic degassing carbon needed to keep CO2 in equilibrium [kgC]/yr
      real(wp) :: dC13volc_dt   !! volcanic degassing carbon 13 [kgC/yr]
      real(wp) :: dCemis_dt     !! carbon emissions [kgC/yr]
      real(wp) :: dC13emis_dt   !! carbon 13 emissions [kgC/yr]
      real(wp) :: d13C_emis     !! d13C of CO2 emissions [permil]
      real(wp) :: Cemis_cum     !! cumulated carbon emissions [GtC]
      real(wp) :: dC14dec_dt    !! radiocarbon decay [kgC/yr]
      real(wp) :: dC14prod_dt   !! radiocarbon production [kgC/yr]
      real(wp) :: dCH4_dt      !! carbon flux from oxidation of anthropogenic and ocean (hydrate) CH4 to CO2 [kgC/yr]
      real(wp) :: dC13H4_dt    !! 13C of it (ocean hydrate CH4 only, anthropogenic CH4 carries no 13C) [kgC13/yr]
      real(wp) :: dCemis_extra_dt     !! carbon emissions from extra feedbacks [kgC/yr]
      real(wp) :: Cemis_extra_cum     !! cumulated carbon emissions from extra feedbacks [GtC]
      real(wp) :: T_glob        !! global mean annual surface air temperature [K]
      real(wp) :: dT_glob_cum   !! cumulative global mean annual surface air temperature change [K]
    end type

    private
    public :: co2_class

end module
