!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : b g c _ d e f
!
!  Purpose : definition of ocean biogeochemistry model class 
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
module bgc_def

    use precision, only : wp
    USE bgc_params, only: natm, nocetra, nsedtra, npowtra

    implicit none
    
    integer, parameter :: bgc_ntra = nocetra   ! name for the public interface (to intialize the ocean model)

    INTEGER, PARAMETER :: kt = 23              ! tracer layers (edit as necessary: must be large enough)
    INTEGER, PARAMETER :: ks = 12              ! sediment layers

    !! column type for the carbon pool and related variables
    ! ------------------------------------------------------
    TYPE carb_t

      REAL(wp) :: ocetra(kt, nocetra)  !! ocean tracers (column)
      REAL(wp) :: ocetra_cons(kt, nocetra)  !! ocean tracers compensation for negative values (column)

      REAL(wp) :: solco2

      REAL(wp) :: aksp(kt)

      REAL(wp) :: AKS3(Kt) 
      REAL(wp) :: AKF3(Kt)   
      REAL(wp) :: AK1P3(Kt)
      REAL(wp) :: AK2P3(Kt)
      REAL(wp) :: AK3P3(Kt)
      REAL(wp) :: AKSI3(Kt)
      REAL(wp) :: akw3(kt)
      REAL(wp) :: akb3(kt)
      REAL(wp) :: ak13(kt)
      REAL(wp) :: ak23(kt)
      REAL(wp) :: satoxy(kt)
      REAL(wp) :: satn2
      REAL(wp) :: satn2o
      REAL(wp) :: atdifv
      REAL(wp) :: suppco2
      REAL(wp) :: co2aq
      REAL(wp) :: aksurf(10)    !! dissociation constants for DIC,bor and water at the sea surface, no pressure dependency 
      ! n2budget closes the mass balance for the alkalinity for biogenic induced changes in N2
      REAL(wp) :: n2budget(kt)
      ! h2obudget closes the mass balance for the oxygen for biogenic induced changes in N2
      REAL(wp) :: h2obudget(kt)

      real(wp) :: wpoc(kt)
      real(wp) :: wcal(kt)
      real(wp) :: wopal(kt)
      real(wp) :: wdust(kt)
      real(wp) :: visco(kt)

    END TYPE

    ! input fluxes for one column
    TYPE flux_t
      ! for input of dust fields
      REAL(wp) :: dustdep 
      REAL(wp) :: dustdep_volcano  = 0 ! (default to zero)

      REAL(wp) :: inp_doc      !! dissolved organic carbon input from river flux
      REAL(wp) :: inp_doc13    !! dissolved organic carbon 13 input from river flux
      REAL(wp) :: inp_doc14    !! dissolved organic carbon 14 input from river flux
      REAL(wp) :: inp_poc      !! particulate organic carbon input from river flux
      REAL(wp) :: inp_poc13    !! particulate organic carbon 13 input from river flux
      REAL(wp) :: inp_poc14    !! particulate organic carbon 14 input from river flux
      REAL(wp) :: inp_dic      !! dissolved inorganic carbon input from river flux
      REAL(wp) :: inp_dic13    !! dissolved inorganic carbon 13 input from river flux
      REAL(wp) :: inp_dic14    !! dissolved inorganic carbon 14 input from river flux
      REAL(wp) :: inp_alk      !! alkalinity input from river flux
      REAL(wp) :: inp_sil      !! silica input from river flux
    END TYPE

    ! climate input type
    type clim_t
      real(wp), dimension(:,:,:), allocatable :: temp       ! potential temperature [deg C], i=1 at surface
      real(wp), dimension(:,:,:), allocatable :: sal        ! salinity [psu]
      real(wp), dimension(:,:),   allocatable :: fw             ! net freshwater flux at the ocean surface (+ive into the ocean) [kg/m2/s]
      real(wp), dimension(:,:),   allocatable :: swr_toa_24h     ! diurnal cycle of SW radiation at TOA [W/m2]
      real(wp), dimension(:),     allocatable :: daylength        ! daylength [hours]
      real(wp), dimension(:,:),   allocatable :: swr_sur_dw     ! short wave flux over ocean [W/m2]
      real(wp), dimension(:,:),   allocatable :: ice_frac      ! sea ice fraction [1]
      real(wp), dimension(:,:),   allocatable :: pressure      ! surface pressure [Pa]
      real(wp), dimension(:,:),   allocatable :: wind10           ! wind speed at 10m [m/s]
      real(wp), dimension(:,:),   allocatable :: sst_min    ! minimum annual surface layer temperature [degC]
      real(wp), dimension(:,:),   allocatable :: sst_max    ! maximum annual surface layer temperature [degC]
      real(wp) :: t2m_glob_ann                              ! global mean annual temperature [degC]
    end type

    ! daily climate input type
    type daily_input_save_t
      real(wp), dimension(:,:,:,:), allocatable :: temp       ! potential temperature [deg C], i=1 at surface
      real(wp), dimension(:,:,:,:), allocatable :: sal        ! salinity [psu]
      real(wp), dimension(:,:,:),   allocatable :: fw             ! net freshwater flux at the ocean surface (+ive into the ocean) [kg/m2/s]
      real(wp), dimension(:,:,:),   allocatable :: swr_sur_dw     ! short wave flux over ocean [W/m2]
      real(wp), dimension(:,:,:),   allocatable :: ice_frac      ! sea ice fraction [1]
      real(wp), dimension(:,:,:),   allocatable :: pressure      ! surface pressure [Pa]
      real(wp), dimension(:,:,:),   allocatable :: wind10           ! wind speed at 10m [m/s]
      real(wp), dimension(:,:,:),   allocatable :: dust_dep          ! dust deposition [kg/m2/s]
    end type

    ! sediment flux type
    type flx_sed_t
      real(wp), dimension(:,:),   allocatable  :: poc
      real(wp), dimension(:,:),   allocatable  :: poc13
      real(wp), dimension(:,:),   allocatable  :: poc14
      real(wp), dimension(:,:),   allocatable  :: caco3
      real(wp), dimension(:,:),   allocatable  :: caco313
      real(wp), dimension(:,:),   allocatable  :: caco314
      real(wp), dimension(:,:),   allocatable  :: opal
      real(wp) :: poc_tot
      real(wp) :: poc13_tot
      real(wp) :: poc14_tot
      real(wp) :: caco3_tot
      real(wp) :: caco313_tot
      real(wp) :: caco314_tot
      real(wp) :: opal_tot
    end type
    
    ! Sediment type, spatially variable (one instance per water column)
    TYPE sedmnt_t

      REAL(wp) :: sedlay (ks, nsedtra)    !! solid sediment tracers
      REAL(wp) :: powtra (ks, npowtra)    !! sediment pore water tracers
      REAL(wp) :: burial(nsedtra)     !! burial layer
      REAL(wp) :: burdiag(nsedtra)     !! burial fluxes diagnostic
      ! pown2bud closes the mass balance for the alkalinity for biogenic induced changes in N2 in sediment
      REAL(wp) :: pown2bud(ks)
      ! powh2obud closes the mass balance for oxygen
      REAL(wp) :: powh2obud(ks)
      REAL(wp) :: sedhpl(ks)  !! sediment pH
      REAL(wp) :: powcar(ks)  !! sediment carbonate ion concentration

      REAL(wp) :: silpro
      REAL(wp) :: prorca  
      REAL(wp) :: prcaca 
      REAL(wp) :: pror13    
      REAL(wp) :: prca13
      REAL(wp) :: pror14
      REAL(wp) :: prca14
      REAL(wp) :: produs

    END TYPE

    !! BGC type defined for one water column
    !  -------------------------------------
    type bgc_1d_t

      ! time control parameters
      integer  :: ldtrunbgc               !!  actual time steps of run.

      integer :: kbo                !! index of bottom cell
      real(wp) :: bolay             !! thickness of bottom cell (--> todo: remove)

      type(carb_t) :: tra      !! carbon pool
      type(flux_t) :: flux     !! carbon flux
      type(sedmnt_t) :: sed    !! bottom sediment

      real(wp) :: delta_c !! cumulated air-sea carbon flux over one year [kgc/s]
      real(wp) :: delta_c13 !! cumulated air-sea carbon 13 flux over one year [kgc13/s]
      real(wp) :: delta_c14 !! cumulated air-sea carbon 14 flux over one year [kgc14/s]

    end type

    ! Vertical grid variables
    type grid_t
      integer :: ni, nj, nk                                    !! grid dimensions
      integer :: ncells                                        !! total number of active ocean columns
      integer, dimension(:), allocatable :: idx_cell_active    !! index of active ocean cells
      integer, dimension(:, :), allocatable :: ij_1d           !! n --> i,j  (2 x ncol)
      integer, dimension(:, :), allocatable :: id_map          !! i,j --> n  (ni, nj)
      integer, dimension(:,:), allocatable :: mask_ocn        !! ocean mask
      integer, dimension(:,:), allocatable :: kbo     !! index of bottom layer
      real(wp), dimension(:), allocatable :: layer_thk(:)          !! layer thickness [m]
      real(wp), dimension(:), allocatable :: level_depth(:)        !! level depth (at layer interface), downward [m] - size + 1
      real(wp), dimension(:), allocatable :: layer_depth(:)        !! layer depth (at layer center), downward [m]
      real(wp), dimension(:), allocatable :: lat
      real(wp), dimension(:), allocatable :: lon
      real(wp) :: ocn_area_tot
      real(wp), dimension(:,:), allocatable :: ocn_area
      real(wp), dimension(:,:), allocatable :: ocn_area_old
      real(wp) :: ocn_vol_tot
      real(wp), dimension(:,:,:), allocatable :: ocn_vol
      real(wp), dimension(:,:,:), allocatable :: mask3d
      real(wp), dimension(:,:,:), allocatable :: coral_f_area
      real(wp), dimension(:,:,:), allocatable :: coral_f_topo
    end type

    type bgc_class
      type(grid_t) :: grid
      type(bgc_1d_t), dimension(:), allocatable :: bgc_1d      !! n column models (ncol)
      type(clim_t) :: clim                                 !! climate input fields
      type(daily_input_save_t) :: daily_input_save          !! daily input fields over one year
      type(flx_sed_t) :: flx_sed                           !! fluxes to sediment for spinup
      logical, dimension(nocetra) :: l_trans_tracers           !! flag for tracer transport to be passed to ocean model
      logical, dimension(nocetra) :: l_tracer_dic           !! flag for tracer transport to be passed to ocean model
      logical, dimension(nocetra) :: l_isodiff_tracers           !! flag for isopycnal diffusion to be passed to ocean model
      real(wp), dimension(natm) :: atm                         !! global mean atmosphere
      real(wp) :: delta_c !! cumulated air-sea carbon flux over one year [kgc]
      real(wp) :: delta_c13 !! cumulated air-sea carbon 13 flux over one year [kgc13]
      real(wp) :: delta_c14 !! cumulated air-sea carbon 14 flux over one year [kgc14]
      real(wp) :: n2o_flx  !! cumulated air-sea N2O flux over one year [kgN2O]
      real(wp) :: Cflx_avg !! average air-sea carbon flux [kgc/yr]
      real(wp) :: A_shelf !! flooded shelf area with a water depth between 0 and 100 m [m2]
      real(wp) :: A_shelf_old !! old flooded shelf area with a water depth between 0 and 100 m [m2]
      real(wp) :: dic_tot_ini
      real(wp) :: alk_tot_ini
      real(wp) :: po4_tot_ini
      real(wp) :: sil_tot_ini
      real(wp) :: dic_tot
      real(wp) :: alk_tot
      real(wp) :: po4_tot
      real(wp) :: sil_tot
    end type

    private
    public :: bgc_class, bgc_ntra

end module
