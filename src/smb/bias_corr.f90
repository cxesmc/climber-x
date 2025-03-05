!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : b i a s _ c o r r _ m o d 
!
!  Purpose : climate bias correction
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
module bias_corr_mod

   use precision, only : wp, dp
   use timer, only: doy, nday_year, monthly2daily
   use control, only : out_dir
   use smb_params, only : smb_ref_file, smb_cx_ref_file 
   use coord, only : grid_init, grid_class
   use coord, only : map_scrip_init, map_scrip_class, map_scrip_field
   use ncio

   implicit none

   real(wp), dimension(:,:), allocatable :: t2m_bias_ann
   real(wp), dimension(:,:,:), allocatable :: t2m_bias_mon
   real(wp), dimension(:,:,:), allocatable :: prc_bias_mon
   integer, dimension(nday_year) :: m0, m1
   real(wp), dimension(nday_year) :: wtm0, wtm1

   private
   public :: smb_ref_write, smb_bias_corr, t2m_bias_corr, prc_bias_corr

contains

   subroutine smb_ref_write(grid, year_ini_smb_ref, year_end_smb_ref, &
       ann_smb_avg_ref, ann_prc_avg_ref, ann_evp_avg_ref)

   implicit none

   type(grid_class), intent(in) :: grid
   integer, intent(in) :: year_ini_smb_ref
   integer, intent(in) :: year_end_smb_ref
   real(wp), intent(in) :: ann_smb_avg_ref(:,:)
   real(wp), intent(in) :: ann_prc_avg_ref(:,:)
   real(wp), intent(in) :: ann_evp_avg_ref(:,:)

   integer :: ncid
   character (len=256) :: fnm
   character (len=256) :: i1, i2

   ! Create the netcdf file and the dimension variables
   write (i1,'(I4)') year_ini_smb_ref+2000
   write (i2,'(I4)') year_end_smb_ref+2000
   fnm = trim(out_dir)//"/CLIMBER-X_SMB_"//trim(grid%name)//"_clim_"//trim(i1)//"_"//trim(i2)//".nc"
   call nc_create(fnm)
   call nc_open(fnm,ncid)
   call nc_write_dim(fnm, "y", x=grid%G%y0, dx=grid%G%dy, nx=grid%G%ny,&
     axis="y", units="km", ncid=ncid)
   call nc_write_dim(fnm, "x", x=grid%G%x0, dx=grid%G%dx, nx=grid%G%nx,&
     axis="x", units="km", ncid=ncid)
   call nc_write(fnm,"phi",sngl(grid%lat), dims=["x","y"],start=[1,1],count=[grid%G%nx,grid%G%ny],ncid=ncid)
   call nc_write(fnm,"lambda",sngl(grid%lon), dims=["x","y"],start=[1,1],count=[grid%G%nx,grid%G%ny],ncid=ncid)
   call nc_write(fnm,"smb", sngl(ann_smb_avg_ref), dims=["x","y"],start=[1,1],count=[grid%G%nx,grid%G%ny],long_name="Climatology of annual surface mass balance",grid_mapping="polar_stereographic",units="kg/m2/yr",ncid=ncid) 
   call nc_write(fnm,"prc", sngl(ann_prc_avg_ref), dims=["x","y"],start=[1,1],count=[grid%G%nx,grid%G%ny],long_name="Climatology of annual precipitation",grid_mapping="polar_stereographic",units="kg/m2/yr",ncid=ncid) 
   call nc_write(fnm,"evp", sngl(ann_evp_avg_ref), dims=["x","y"],start=[1,1],count=[grid%G%nx,grid%G%ny],long_name="Climatology of annual sublimation",grid_mapping="polar_stereographic",units="kg/m2/yr",ncid=ncid) 
   call nc_close(ncid)

   return

   end subroutine smb_ref_write


   subroutine smb_bias_corr(grid, &
        ann_smb_ref, ann_smb_cx_ref, &
        ann_prc_ref, ann_prc_cx_ref, &
        ann_evp_ref, ann_evp_cx_ref)

   implicit none

   type(grid_class), intent(in) :: grid
   real(wp), intent(out) :: ann_smb_ref(:,:)
   real(wp), intent(out) :: ann_smb_cx_ref(:,:)
   real(wp), intent(out) :: ann_prc_ref(:,:)
   real(wp), intent(out) :: ann_prc_cx_ref(:,:)
   real(wp), intent(out) :: ann_evp_ref(:,:)
   real(wp), intent(out) :: ann_evp_cx_ref(:,:)

   integer :: ni, nj, nlon, nlat
   integer :: ppos, spos
   real(wp) :: dlon, dlat
   real(wp), dimension(:), allocatable :: lon, lat
   real(wp), dimension(:,:), allocatable :: smb, prc, evp
   character(len=256) :: fnm
   type(grid_class) :: grid_ref
   type(map_scrip_class) :: maps_to_smb


   ! read CLIMBER-X reference SMB file, has to be on the same grid as the current domain
   fnm = smb_cx_ref_file 
   call nc_read(fnm,"smb",ann_smb_cx_ref) 
   call nc_read(fnm,"prc",ann_prc_cx_ref) 
   call nc_read(fnm,"evp",ann_evp_cx_ref) 

   ! read reference SMB file, is NOT in the same grid as the current domain, mapping needed!
   fnm = smb_ref_file 
   ! read dimensions
   ni = nc_size(trim(fnm),"lon")
   nj = nc_size(trim(fnm),"lat")
   ! read lat lon
   allocate(lon(ni))
   allocate(lat(nj))
   call nc_read(trim(fnm),"lon",lon)
   call nc_read(trim(fnm),"lat",lat)
   dlon = lon(2)-lon(1)
   dlat = lat(2)-lat(1)

   allocate(smb(ni,nj))
   allocate(prc(ni,nj))
   allocate(evp(ni,nj))
   call nc_read(fnm,"smb",smb(:,:))
   call nc_read(fnm,"prc",prc(:,:))
   call nc_read(fnm,"evp",evp(:,:))

   ! mapping to current smb domain

   ! create grid object
   spos = scan(trim(fnm),"/", BACK= .true.)+1
   ppos = scan(trim(fnm),".", BACK= .true.)-1
   call grid_init(grid_ref,name=trim(fnm(spos:ppos)),mtype="latlon",units="degrees", &
     x0=real(lon(1),dp),dx=real(dlon,dp),nx=ni,y0=real(lat(1),dp),dy=real(dlat,dp),ny=nj)

   ! initialize map
   call map_scrip_init(maps_to_smb,grid_ref,grid,method="con",fldr="maps",load=.TRUE.,clean=.FALSE.)

   ! mapping
   call map_scrip_field(maps_to_smb,"smb",smb,ann_smb_ref,method="mean") 
   call map_scrip_field(maps_to_smb,"prc",prc,ann_prc_ref,method="mean") 
   call map_scrip_field(maps_to_smb,"evp",evp,ann_evp_ref,method="mean") 

   deallocate(smb)
   deallocate(prc)
   deallocate(evp)


   return

   end subroutine smb_bias_corr


   subroutine t2m_bias_corr(i_t2m_bias_corr, bias_corr_file, t2m_bias)

   implicit none

   integer, intent(in) :: i_t2m_bias_corr
   character (len=*), intent(in) :: bias_corr_file

   real(wp), intent(out) :: t2m_bias(:,:,:)

   integer :: d, nmon, nlon, nlat


   if (i_t2m_bias_corr.eq.1) then

     ! read temperature bias file
     nlon = nc_size(trim(bias_corr_file),"lon")
     nlat = nc_size(trim(bias_corr_file),"lat")
     allocate( t2m_bias_ann(nlon,nlat) )
     call nc_read(trim(bias_corr_file),"t2m_bias",t2m_bias_ann) 

     do d=1,nday_year
       t2m_bias(:,:,d) = t2m_bias_ann 
     enddo

   else if (i_t2m_bias_corr.eq.2) then

     ! get weights for interpolation from monthly to daily
     call monthly2daily(m0,m1,wtm0,wtm1)

     ! read temperature bias file
     nmon = nc_size(trim(bias_corr_file),"mon")
     if (nmon.gt.13) stop 'bias_corr ERROR: i_bias_corr==2 requires input file with monthly values' 
     nlon = nc_size(trim(bias_corr_file),"lon")
     nlat = nc_size(trim(bias_corr_file),"lat")
     allocate( t2m_bias_mon(nlon,nlat,0:nmon+1) )
     call nc_read(trim(bias_corr_file),"t2m_bias",t2m_bias_mon(:,:,1:12),start=[1,1,1],count=[nlon,nlat,12]) 

     ! expand month for easier interpolation at borders
     t2m_bias_mon(:,:,0)  = t2m_bias_mon(:,:,12)
     t2m_bias_mon(:,:,13) = t2m_bias_mon(:,:,1)

     do d=1,nday_year
       t2m_bias(:,:,d) = wtm0(d)*t2m_bias_mon(:,:,m0(d)) + wtm1(d)*t2m_bias_mon(:,:,m1(d))
     enddo

   endif

   return

   end subroutine t2m_bias_corr


   subroutine prc_bias_corr(i_prc_bias_corr, bias_corr_file, prc_bias)

   implicit none

   integer, intent(in) :: i_prc_bias_corr
   character (len=*), intent(in) :: bias_corr_file

   real(wp), intent(out) :: prc_bias(:,:,:)

   integer :: d, nmon, nlon, nlat


   if (i_prc_bias_corr.eq.1) then

     ! read temperature bias file
     nlon = nc_size(trim(bias_corr_file),"lon")
     nlat = nc_size(trim(bias_corr_file),"lat")
     nmon = nc_size(trim(bias_corr_file),"mon")
     allocate( prc_bias_mon(nlon,nlat,nmon) )
     call nc_read(trim(bias_corr_file),"prc_bias",prc_bias_mon) 

     do d=1,nday_year
       prc_bias(:,:,d) = sum(prc_bias_mon(:,:,5:9),3)/5._wp 
     enddo

   else if (i_prc_bias_corr.eq.2) then

     ! get weights for interpolation from monthly to daily
     call monthly2daily(m0,m1,wtm0,wtm1)

     ! read temperature bias file
     nmon = nc_size(trim(bias_corr_file),"mon")
     if (nmon.gt.13) stop 'bias_corr ERROR: i_bias_corr==2 requires input file with monthly values' 
     nlon = nc_size(trim(bias_corr_file),"lon")
     nlat = nc_size(trim(bias_corr_file),"lat")
     allocate( prc_bias_mon(nlon,nlat,0:nmon+1) )
     call nc_read(trim(bias_corr_file),"prc_bias",prc_bias_mon(:,:,1:12),start=[1,1,1],count=[nlon,nlat,12]) 

     ! expand month for easier interpolation at borders
     prc_bias_mon(:,:,0)  = prc_bias_mon(:,:,12)
     prc_bias_mon(:,:,13) = prc_bias_mon(:,:,1)

     do d=1,nday_year
       prc_bias(:,:,d) = wtm0(d)*prc_bias_mon(:,:,m0(d)) + wtm1(d)*prc_bias_mon(:,:,m1(d))
     enddo

   endif

   return

   end subroutine prc_bias_corr

end module bias_corr_mod
