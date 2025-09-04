!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : b n d
!
!  Purpose : boundary conditions
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
module bnd_mod

    use precision, only : wp
    use dim_name, only: dim_lon, dim_lat
    use timer, only : year, year_ini, year_now, nyears, time_soy_bnd, n_year_geo
    use climber_grid, only : ni, nj, lat
    use control, only : iorbit, orbit_file
    use control, only : isol, sol_const, sol_file
    use control, only : ivolc, volc_const, volc_file, volc_scale
    use control, only : isea_level, sea_level_const, sea_level_file
    use control, only : ico2, co2_const, co2_file
    use control, only : ico2_rad, co2_rad_const, co2_rad_file
    use control, only : id13c, d13c_atm_const, d13c_atm_file
    use control, only : iD14c, D14c_atm_const, D14c_atm_file
    use control, only : ich4, ch4_const, ch4_file
    use control, only : ich4_rad, ch4_rad_const, ch4_rad_file
    use control, only : in2o, n2o_const, n2o_file
    use control, only : in2o_rad, n2o_rad_const, n2o_rad_file
    use control, only : iso4, so4_const, so4_file
    use control, only : io3, o3_const, o3_file_const, o3_file_var
    use control, only : icfc, cfc11_const, cfc12_const, cfc_file
    use control, only : iluc, luc_file
    use control, only : idist, dist_file
    use control, only : ocn_restore_temp, ocn_restore_sal
    use control, only : atm_fix_tau
    use control, only : flag_atm, flag_ocn, flag_sic, flag_lnd, flag_dust, flag_ice, ifake_ice, flag_geo, flag_smb
    use control, only : l_feedbacks, l_spinup_cc
    use constants, only : c13_c12_std, c14_c_std
    use sol_mod
    use volc_mod
    use insolation
    use sea_level_mod
    use co2_mod
    use co2_rad_mod
    use d13c_atm_mod
    use D14c_atm_mod
    use ch4_mod
    use ch4_rad_mod
    use n2o_mod
    use n2o_rad_mod
    use so4_mod
    use o3_mod
    use cfc_mod
    use luc_mod
    use dist_mod
    use fake_atm_mod
    use fake_dust_mod
    use fake_lnd_mod
    use fake_ocn_mod
    use fake_sic_mod
    use fake_ice_mod
    use fake_geo_mod
    use ncio
    use coord, only : grid_class

    implicit none 

    type bnd_class 
     real(wp) :: eccentricity
     real(wp) :: precession
     real(wp) :: obliquity
     real(wp), dimension(:,:,:), allocatable :: solar
     real(wp), dimension(:,:), allocatable :: solarm
     real(wp), dimension(:,:), allocatable :: solarmin
     real(wp), dimension(:,:), allocatable :: solarmax
     real(wp), dimension(:,:,:), allocatable :: cosz
     real(wp), dimension(:,:), allocatable :: coszm
     real(wp), dimension(:,:), allocatable :: daylength
     type(fake_atm_type) :: atm
     type(fake_dust_type) :: dust
     type(fake_lnd_type) :: lnd
     type(fake_ocn_type) :: ocn
     type(fake_sic_type) :: sic
     type(fake_ice_type) :: ice
     type(fake_geo_type) :: geo
     real(wp) :: sol
     real(wp) :: volc
     real(wp) :: co2, co2_rad, d13c_atm, D14c_atm, c13_c12_atm, c14_c_atm
     real(wp) :: ch4, ch4_rad
     real(wp) :: n2o, n2o_rad
     real(wp), dimension(:,:), allocatable :: so4
     real(wp), dimension(:), allocatable :: o3_pl
     real(wp), dimension(:,:,:), allocatable :: o3
     real(wp), dimension(:,:), allocatable :: f_crop
     real(wp), dimension(:,:), allocatable :: f_pasture
     real(wp), dimension(:,:,:), allocatable :: disturbance
     real(wp) :: sea_level
     real(wp) :: cfc11, cfc12
    end type 


    private
    public :: bnd_class
    public :: bnd_init, bnd_update, bnd_end


contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  b n d _ u p d a t e
  ! Purpose  :  update bnd
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bnd_update(bnd)

    implicit none 

    type(bnd_class) :: bnd

    real(wp) :: sol
    real(wp) :: year_orb


    ! updates at start of new year
    if (time_soy_bnd) then

      ! update sea level if required
      if( isea_level.eq.1 .or. isea_level.eq.2 ) then
        call sea_level_update(real(year_now,kind=wp),bnd%sea_level)     
        print *,'sea level:',bnd%sea_level
      endif

      ! update solar 'constant' if required 
      if( isol.eq.1 ) then
        call sol_update(isol, real(year_now,kind=wp), bnd%sol)     
      endif

      ! update volcanic forcing if required 
      if( ivolc.eq.1 ) then
        call volc_update(ivolc, real(year_now,kind=wp), bnd%volc)     
      endif

      ! apply volcanic radiative forcing to solar constant (convert radiative forcing to solar constant anomalies!)
      sol = bnd%sol + volc_scale * 4._wp/0.7_wp*bnd%volc

      ! get current orbital forcing if required
      if( iorbit.eq.2 .or. isol.eq.1 .or. ivolc.eq.1 ) then
        if (iorbit.eq.2) then
          year_orb = real(year_now,kind=wp)
        else
          year_orb = real(year_ini,kind=wp)
        endif
        call sinsol(sol, year_orb, lat, &
          bnd%eccentricity, bnd%precession, bnd%obliquity, bnd%solar, bnd%solarm, bnd%solarmin, bnd%solarmax, bnd%cosz, bnd%coszm)
        call day_length(lat,bnd%daylength)
      endif

      ! update atmospheric CO2 concentration if required
      if( ico2.ge.1) then
        call co2_update(ico2, real(year_now,kind=wp), bnd%co2)     
      endif
      ! if feedback analysis
      if (l_feedbacks) then
        if (year.gt.nyears/2)   bnd%co2 = co2_const
      endif

      ! update atmospheric CO2 concentration for radiation if required
      if( ico2_rad.eq.2 .or. ico2_rad.ge.3) then
        call co2_rad_update(ico2_rad,real(year_now,kind=wp), bnd%co2_rad)
      endif

      ! update atmospheric d13C if required
      if( id13c .eq. 1 ) then
        call d13c_atm_update(real(year_now,kind=wp),bnd%d13c_atm)
        bnd%c13_c12_atm = (bnd%d13c_atm/1000._wp+1._wp)*c13_c12_std
      endif

      ! update atmospheric D14C if required
      if( iD14c .eq. 1 ) then
        call D14c_atm_update(real(year_now,kind=wp),bnd%D14c_atm)
        bnd%c14_c_atm = (bnd%D14c_atm/1000._wp+1._wp)*c14_c_std
      endif

      ! update atmospheric CH4 concentration if required
      if( ich4.eq.1 .or. ich4.eq.2) then
        call ch4_update(ich4, real(year_now,kind=wp), bnd%ch4)     
      endif

      ! update atmospheric CH4 concentration for radiation if required
      if( ich4_rad.eq.2) then
        call ch4_rad_update(ich4_rad, real(year_now,kind=wp), bnd%ch4_rad)     
      endif

      ! update atmospheric N2O concentration if required
      if( in2o.eq.1 .or. in2o.eq.2) then
        call n2o_update(in2o, real(year_now,kind=wp), bnd%n2o)     
      endif

      ! update atmospheric N2O concentration for radiation if required
      if( in2o_rad.eq.2) then
        call n2o_rad_update(in2o_rad, real(year_now,kind=wp), bnd%n2o_rad)     
      endif

      ! update atmospheric SO4 load if required
      if( iso4.eq.1) then
        call so4_update(iso4, real(year_now,kind=wp), bnd%so4)     
      endif

      ! update atmospheric O3 concentration if required
      if( io3.eq.1) then
        call o3_update(io3, real(year_now,kind=wp), bnd%o3)     
      endif

      ! update atmospheric CFCs concentration if required
      if( icfc.eq.1 ) then
        call cfc_update(real(year_now,kind=wp),bnd%cfc11,bnd%cfc12)
      endif

      ! update land use state if required
      if( iluc.eq.2 ) then
        call luc_update(iluc, real(year_now,kind=wp), bnd%f_crop, bnd%f_pasture)     
      endif

      ! update vegetation disturbance rate if required
      if( idist.ge.2 ) then
        call dist_update(idist, real(year_now,kind=wp), bnd%disturbance)     
      endif

      ! update fake ice (ice sheet thickness)
      if( .not.flag_ice .or. ifake_ice.eq.1 ) then
        if (mod(year,n_year_geo).eq.1) then
          call fake_ice_update(real(year_now,kind=wp),bnd%ice)
        endif
      endif

      ! update fake geo (bedrock topography) 
      if( .not.flag_geo ) then
        if (mod(year,n_year_geo).eq.1) then
          call fake_geo_update(real(year_now,kind=wp),bnd%geo)
        endif
      endif

    endif

    ! update fake atmosphere
    if( (.not.flag_atm .and. .not.l_spinup_cc) .or. flag_smb .or. atm_fix_tau) then
      call fake_atm_update(real(year_now,kind=wp),bnd%atm)
    endif

    ! update fake dust
    if( (.not.flag_dust) ) then
      call fake_dust_update(real(year_now,kind=wp),bnd%dust)
    endif

    ! update fake land
    if( .not.flag_lnd .and. .not.l_spinup_cc) then
      call fake_lnd_update(bnd%lnd)
    endif

    ! update fake ocean
    if( (.not.flag_ocn .and. .not.l_spinup_cc) .or. ocn_restore_sal .or. ocn_restore_temp) then
      call fake_ocn_update(real(year_now,kind=wp),bnd%ocn)
    endif

    ! update fake sea ice
    if( .not.flag_sic.and. .not.l_spinup_cc ) then
      call fake_sic_update(real(year_now,kind=wp),bnd%sic)
    endif


    return

  end subroutine bnd_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  b n d _ i n i t
  ! Purpose  :  initialize bnd
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bnd_init(cmn_grid,bnd)

    implicit none 

    type(grid_class), intent(in) :: cmn_grid
    type(bnd_class) :: bnd 


    ! allocate variables
    call bnd_alloc(bnd)

    ! initialize solar 'constant'
    if (isol.eq.0) then
      bnd%sol = sol_const
    elseif (isol.eq.1) then
      ! read solar constant file
      call sol_init(sol_file)
      ! get initial solar constant value
      call sol_update(isol,real(year_ini,kind=wp),bnd%sol)
    endif

    ! initialize volcanic forcing
    if (ivolc.eq.0) then
      bnd%volc = volc_const
    elseif (ivolc.eq.1) then
      ! read volcainc forcing file
      call volc_init(volc_file)
      ! get initial volcanic forcing value
      call volc_update(ivolc,real(year_ini,kind=wp),bnd%volc)
    endif

    ! apply volcanic radiative forcing to solar constant (convert radiative forcing to solar constant anomalies!)
    bnd%sol = bnd%sol + volc_scale * 4._wp/0.7_wp*bnd%volc

    ! initialize insolation quantities
    call ini_sinsol(orbit_file)
    ! get initial insolation quantities
    call sinsol(bnd%sol, real(year_ini,kind=wp),lat, &
        bnd%eccentricity, bnd%precession, bnd%obliquity, bnd%solar, bnd%solarm, bnd%solarmin, bnd%solarmax, bnd%cosz, bnd%coszm)
    call day_length(lat,bnd%daylength)

    ! initialize sea level
    if (isea_level.eq.0 .or. isea_level.eq.2) then
      bnd%sea_level = sea_level_const
    elseif (isea_level.eq.1) then
      ! read sea level file
      call sea_level_init(sea_level_file)
      ! get initial sea level value
      call sea_level_update(real(year_ini,kind=wp),bnd%sea_level)
    endif

    ! initialize atmospheric CO2
    if (ico2.eq.0) then
      bnd%co2 = co2_const
    elseif (ico2.eq.1 .or. ico2.eq.-1) then
      ! read co2 file
      call co2_init(co2_file)
      ! get initial co2 value
      call co2_update(ico2,real(year_ini,kind=wp),bnd%co2)
    elseif (ico2.eq.2) then
      bnd%co2 = co2_const
    elseif (ico2.eq.3 .or. ico2.eq.4 .or. ico2.eq.5 .or. ico2.eq.6 .or. ico2.eq.7 .or. ico2.eq.8) then
      bnd%co2 = co2_const
    elseif (ico2.eq.9) then
      bnd%co2 = 280._wp
    endif
    ! if feedback analysis
    if (l_feedbacks) then
      bnd%co2 = co2_const*2._wp
    endif

    ! initialize atmospheric CO2 for radiation
    if (ico2_rad.eq.1) then
      bnd%co2_rad = co2_rad_const
    elseif (ico2_rad.eq.2) then
      ! read co2_rad file
      call co2_rad_init(co2_rad_file)
      ! get initial co2_rad value
      call co2_rad_update(ico2_rad,real(year_ini,kind=wp),bnd%co2_rad)
    elseif (ico2_rad.ge.3) then
      bnd%co2_rad = co2_rad_const
    endif

    ! initialize atmospheric d13C
    if (id13c.eq.0) then
      bnd%d13c_atm = d13c_atm_const
    elseif (id13c.eq.1) then
      ! read d13c file
      call d13c_atm_init(d13c_atm_file)
      ! get initial d13c value
      call d13c_atm_update(real(year_ini,kind=wp),bnd%d13c_atm)
    endif

    ! initialize atmospheric D14C
    if (iD14c.eq.0) then
      bnd%D14c_atm = D14c_atm_const
    elseif (iD14c.eq.1) then
      ! read D14c file
      call D14c_atm_init(D14c_atm_file)
      ! get initial D14c value
      call D14c_atm_update(real(year_ini,kind=wp),bnd%D14c_atm)
    endif

    ! initialize atmospheric CH4
    if (ich4.eq.0) then
      bnd%ch4 = ch4_const
    elseif (ich4.eq.1 .or. ich4.eq.-1) then
      ! read ch4 file
      call ch4_init(ch4_file)
      ! get initial ch4 value
      call ch4_update(ich4,real(year_ini,kind=wp),bnd%ch4)
    endif

    ! initialize atmospheric CH4 for radiation
    if (ich4_rad.eq.1) then
      bnd%ch4_rad = ch4_rad_const
    elseif (ich4_rad.eq.2) then
      ! read ch4_rad file
      call ch4_rad_init(ch4_rad_file)
      ! get initial ch4_rad value
      call ch4_rad_update(ich4_rad,real(year_ini,kind=wp),bnd%ch4_rad)
    endif

    ! initialize atmospheric N2O
    if (in2o.eq.0) then
      bnd%n2o = n2o_const
    elseif (in2o.eq.1 .or. in2o.eq.-1) then
      ! read n2o file
      call n2o_init(n2o_file)
      ! get initial n2o value
      call n2o_update(in2o,real(year_ini,kind=wp),bnd%n2o)
    endif

    ! initialize atmospheric N2O for radiation
    if (in2o_rad.eq.1) then
      bnd%n2o_rad = n2o_rad_const
    elseif (in2o_rad.eq.2) then
      ! read n2o_rad file
      call n2o_rad_init(n2o_rad_file)
      ! get initial n2o_rad value
      call n2o_rad_update(in2o_rad,real(year_ini,kind=wp),bnd%n2o_rad)
    endif

    ! initialize atmospheric SO4
    if (iso4.eq.0) then
      bnd%so4(:,:) = so4_const
    elseif (iso4.eq.1) then
      ! read so4 file
      call so4_init(so4_file)
      ! get initial so4 value
      call so4_update(iso4,real(year_ini,kind=wp),bnd%so4)
    endif

    ! initialize atmospheric O3
    bnd%o3_pl = (/1._wp,0.9_wp,0.8_wp,0.7_wp,0.6_wp,0.5_wp,0.4_wp,0.3_wp,0.2_wp,0.1_wp,0._wp/)  ! pressure levels for 3D O3 field
    if (io3.eq.0) then
      bnd%o3(:,:,:) = o3_const
    elseif (io3.eq.1) then
      ! read o3 file
      call o3_init(io3, o3_file_const)
      ! get initial o3 value
      call o3_update(io3,real(year_ini,kind=wp),bnd%o3)
    elseif (io3.eq.2) then
      ! read o3 file
      call o3_init(io3, o3_file_var)
      ! get initial o3 value
      call o3_update(io3,real(year_ini,kind=wp),bnd%o3)
    endif

    if (icfc.eq.0) then
      bnd%cfc11 = cfc11_const
      bnd%cfc12 = cfc12_const
    else if (icfc.eq.1) then
      ! read cfc file
      call cfc_init(cfc_file)
      ! get initial cfc value
      call cfc_update(real(year_ini,kind=wp),bnd%cfc11,bnd%cfc12)
    endif

    ! initialize land use state
    if (iluc.eq.0) then
      bnd%f_crop(:,:)    = 0._wp
      bnd%f_pasture(:,:) = 0._wp
    elseif (iluc.eq.1 .or. iluc.eq.2) then
      ! read luc file
      call luc_init(iluc, luc_file)
      ! get initial luc value
      call luc_update(iluc,real(year_ini,kind=wp),bnd%f_crop,bnd%f_pasture)
    endif

    ! initialize vegetation disturbance rate
    if (idist.eq.0) then
      bnd%disturbance(:,:,:)    = 0._wp
    elseif (idist.eq.1 .or. idist.eq.2) then
      ! read disturbance file
      call dist_init(idist, dist_file)
      ! get initial disturbance value
      call dist_update(idist,real(year_ini,kind=wp),bnd%disturbance)
    else if (idist.ge.3) then
      bnd%disturbance(:,:,:)    = 0._wp
    endif

    ! initialize fake atmosphere
    if (.not.flag_atm .or. flag_smb .or. atm_fix_tau) then
      call fake_atm_init(real(year_ini,kind=wp),bnd%atm)
    endif

    ! initialize fake dust
    if (.not.flag_dust) then
      call fake_dust_init(real(year_ini,kind=wp),bnd%dust)
    endif

    ! initialize fake land
    if (.not.flag_lnd) then
      call fake_lnd_init(bnd%lnd)
    endif
    
    ! initialize fake ocean
    if (.not.flag_ocn .or. ocn_restore_sal .or. ocn_restore_temp) then
      call fake_ocn_init(real(year_ini,kind=wp),bnd%ocn)
    endif

    ! initialize fake sea ice
    if (.not.flag_sic) then
      call fake_sic_init(real(year_ini,kind=wp),bnd%sic)
    endif
    
    ! initialize fake ice (ice sheet thickness)
    call fake_ice_init(real(year_ini,kind=wp),cmn_grid,bnd%ice)

    ! initialize fake geo (bedrock topography)
    call fake_geo_init(real(year_ini,kind=wp),cmn_grid,bnd%sea_level,bnd%geo)

    print*
    print*,'======================================================='
    print*,' Initialisation of BOUNDARY CONDITIONS complete'
    print*,'======================================================='
    print*

    return

  end subroutine bnd_init 


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  b n d _ e n d
  ! Purpose  :  end bnd
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bnd_end(bnd)

    implicit none 

    type(bnd_class) :: bnd

    ! Finalize model variables on domain 'dom'
    ! Call `model_dealloc` to deallocate state variables 
    call bnd_dealloc(bnd)

    return

  end subroutine bnd_end 


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  b n d _ a l l o c
  ! Purpose  :  allocate bnd
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bnd_alloc(bnd)

    implicit none 

    type(bnd_class) :: bnd

    ! Allocate all state variables with dimensions, eg, (ni,nj)

    allocate(bnd%solar(nday_year,24,nj))
    allocate(bnd%solarm(nday_year,nj))
    allocate(bnd%solarmin(nday_year,nj))
    allocate(bnd%solarmax(nday_year,nj))
    allocate(bnd%cosz(nday_year,24,nj))
    allocate(bnd%coszm(nday_year,nj))
    allocate(bnd%daylength(nday_year,nj))

    allocate(bnd%so4(ni,nj))
    allocate(bnd%o3_pl(11))
    allocate(bnd%o3(ni,nj,11))

    allocate(bnd%f_crop(ni,nj))
    allocate(bnd%f_pasture(ni,nj))
    allocate(bnd%disturbance(5,ni,nj))

    return

  end subroutine bnd_alloc 


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  b n d _ d e a l l o c
  ! Purpose  :  deallocate bnd
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bnd_dealloc(bnd)

    implicit none 

    type(bnd_class) :: bnd

    ! Deallocate all state variables to free memory
    deallocate(bnd%solar)
    deallocate(bnd%solarm)
    deallocate(bnd%solarmin)
    deallocate(bnd%solarmax)
    deallocate(bnd%cosz)
    deallocate(bnd%coszm)
    deallocate(bnd%daylength)

    deallocate(bnd%so4)
    deallocate(bnd%o3_pl)
    deallocate(bnd%o3)

    deallocate(bnd%f_crop)
    deallocate(bnd%f_pasture)
    deallocate(bnd%disturbance)

    deallocate(bnd%atm%taux)
    deallocate(bnd%atm%tauy)
    deallocate(bnd%atm%tair)
    deallocate(bnd%atm%tair_min_mon)
    deallocate(bnd%atm%qair)
    deallocate(bnd%atm%prc)
    deallocate(bnd%atm%rain)
    deallocate(bnd%atm%snow)
    deallocate(bnd%atm%pressure)
    deallocate(bnd%atm%swdown)
    deallocate(bnd%atm%lwdown)
    deallocate(bnd%atm%wind)

    deallocate(bnd%ocn%sst)
    deallocate(bnd%ocn%sss)
    deallocate(bnd%ocn%t1l)
    deallocate(bnd%ocn%s1l)
    deallocate(bnd%ocn%t_shelf)
    deallocate(bnd%ocn%t)
    deallocate(bnd%ocn%s)
    deallocate(bnd%ocn%uo1)
    deallocate(bnd%ocn%vo1)

    deallocate(bnd%sic%f_sic)

    deallocate(bnd%lnd%runoff)
    deallocate(bnd%lnd%discharge)


    return

  end subroutine bnd_dealloc 

end module bnd_mod
