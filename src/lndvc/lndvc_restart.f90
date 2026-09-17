module lndvc_restart_m
    !
    ! Restart I/O for the virtual-cell framework: per-vc prognostic state to a
    ! NetCDF file, mirroring src/lnd/lnd_model.f90::lnd_write_restart /
    ! lnd_read_restart. Layout is the reference (nl,ni,nj)-style buffers with an
    ! extra slowest-varying "vc" dimension. Inactive vc slots (desc%class==0)
    ! and unallocated per-class blocks are written as zeros so the file is
    ! rectangular; the reader guards symmetrically before landing values back
    ! into the allocated blocks.
    !

    use precision,     only : wp, sp
    use ncio
    use dim_name,      only : dim_lon, dim_lat, dim_depth, dim_depth0, dim_depth1, &
                              dim_depthl, dim_depth0l, dim_npft, dim_nsurf, dim_nsoil
    use climber_grid,  only : ni, nj, lon, lat
    use lnd_grid,      only : nl, z, z_int, nl_l, z_l, z_int_l, nlc, z_c, &
                              npft, i_pft, nsurf, i_surf, nsoil, i_soil
    use smb_grid_m,    only : nl_smb => nl
    use lndvc_def,     only : lndvc_class

    implicit none

    private
    public :: lndvc_write_restart, lndvc_read_restart

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  l n d v c _ w r i t e _ r e s t a r t
  ! Purpose    :  Persist per-vc prognostic state to a NetCDF restart file
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lndvc_write_restart(fnm, lndvc)

    implicit none

    character(len=*),  intent(in) :: fnm
    type(lndvc_class), intent(in) :: lndvc

    integer :: ncid
    integer :: i, j, k, n_vc
    integer,  dimension(:,:,:),   allocatable :: vari
    real(wp), dimension(:,:,:),   allocatable :: var
    real(wp), dimension(:,:,:,:), allocatable :: var_n
    real(wp), dimension(:,:,:,:,:), allocatable :: var_np

    n_vc = lndvc%n_vc

    call nc_create(fnm)
    call nc_write_dim(fnm,"p",x=1)
    call nc_write_dim(fnm,dim_lon,x=lon,axis="x")
    call nc_write_dim(fnm,dim_lat,x=lat,axis="y")
    call nc_write_dim(fnm,"vc",x=[(k,k=1,n_vc)],units="n/a")
    call nc_write_dim(fnm,"class",x=[1,2,3,4],units="n/a")
    call nc_write_dim(fnm,dim_depth,x=z(1:nl),units="m",axis="z")
    call nc_write_dim(fnm,dim_depth0,x=z(0:nl),units="m",axis="z")
    call nc_write_dim(fnm,dim_depthl,x=z_l(1:nl_l),units="m",axis="z")
    call nc_write_dim(fnm,dim_depth0l,x=z_l(0:nl_l),units="m",axis="z")
    call nc_write_dim(fnm,dim_depth1,x=z_c(1:nlc),units="m",axis="z")
    call nc_write_dim(fnm,"depth_smb0",x=[(k,k=0,nl_smb)],units="n/a",axis="z")
    call nc_write_dim(fnm,dim_npft,x=i_pft,units="n/a")
    call nc_write_dim(fnm,dim_nsurf,x=i_surf,units="n/a")
    call nc_write_dim(fnm,dim_nsoil,x=i_soil,units="n/a")

    call nc_open(fnm,ncid)

    ! ---- A. Global 0-D scalars ----
    call nc_write(fnm,"Cflx_avg",       lndvc%glob%Cflx_avg,       dim1="p",long_name="average land-atm carbon flux",ncid=ncid)
    call nc_write(fnm,"weath_carb_avg", lndvc%glob%weath_carb_avg, dim1="p",long_name="average carbonate weathering flux",ncid=ncid)
    call nc_write(fnm,"weath_sil_avg",  lndvc%glob%weath_sil_avg,  dim1="p",long_name="average silicate weathering flux",ncid=ncid)
    call nc_write(fnm,"landc",          lndvc%glob%landc,          dim1="p",long_name="total land carbon",ncid=ncid)
    call nc_write(fnm,"landc13",        lndvc%glob%landc13,        dim1="p",long_name="total land carbon 13",ncid=ncid)
    call nc_write(fnm,"landc14",        lndvc%glob%landc14,        dim1="p",long_name="total land carbon 14",ncid=ncid)
    call nc_write(fnm,"burc",           lndvc%glob%burc,           dim1="p",long_name="total buried carbon",ncid=ncid)
    call nc_write(fnm,"burc13",         lndvc%glob%burc13,         dim1="p",long_name="total buried carbon 13",ncid=ncid)
    call nc_write(fnm,"burc14",         lndvc%glob%burc14,         dim1="p",long_name="total buried carbon 14",ncid=ncid)
    call nc_write(fnm,"weath_scale",    lndvc%glob%weath_scale,    dim1="p",long_name="weathering scaling factor to match alkalinity input into the ocean",ncid=ncid)
    call nc_write(fnm,"co2",            lndvc%glob%co2,            dim1="p",long_name="atmospheric CO2",ncid=ncid)
    call nc_write(fnm,"c13_c12_atm",    lndvc%glob%c13_c12_atm,    dim1="p",long_name="atmospheric C13/C12 ratio",ncid=ncid)
    call nc_write(fnm,"c14_c_atm",      lndvc%glob%c14_c_atm,      dim1="p",long_name="atmospheric C14/C ratio",ncid=ncid)

    ! ---- B. Cell aggregation ----
    call nc_write(fnm,"mask_lnd",     lndvc%cell%mask_lnd,     dims=[dim_lon,dim_lat],long_name="land mask",units="/",ncid=ncid)
    call nc_write(fnm,"f_land",       lndvc%cell%f_land,       dims=[dim_lon,dim_lat],long_name="land fraction",units="/",ncid=ncid)
    call nc_write(fnm,"f_land0",      lndvc%cell%f_land0,      dims=[dim_lon,dim_lat],long_name="land fraction at present day sea level",units="/",ncid=ncid)
    call nc_write(fnm,"f_ice",        lndvc%cell%f_ice,        dims=[dim_lon,dim_lat],long_name="ice fraction",units="/",ncid=ncid)
    call nc_write(fnm,"f_ice_grd",    lndvc%cell%f_ice_grd,    dims=[dim_lon,dim_lat],long_name="grounded ice fraction",units="/",ncid=ncid)
    call nc_write(fnm,"f_shelf",      lndvc%cell%f_shelf,      dims=[dim_lon,dim_lat],long_name="ocean shelf fraction",units="/",ncid=ncid)
    call nc_write(fnm,"f_lake",       lndvc%cell%f_lake,       dims=[dim_lon,dim_lat],long_name="lake fraction",units="/",ncid=ncid)
    call nc_write(fnm,"f_veg",        lndvc%cell%f_veg,        dims=[dim_lon,dim_lat],long_name="vegetation fraction",units="/",ncid=ncid)
    call nc_write(fnm,"f_ice_old",    lndvc%cell%f_ice_old,    dims=[dim_lon,dim_lat],long_name="previous-step ice fraction",units="/",ncid=ncid)
    call nc_write(fnm,"f_ice_grd_old",lndvc%cell%f_ice_grd_old,dims=[dim_lon,dim_lat],long_name="previous-step grounded ice fraction",units="/",ncid=ncid)
    call nc_write(fnm,"f_ice_nbr",    lndvc%cell%f_ice_nbr,    dims=[dim_lon,dim_lat],long_name="neighbourhood ice fraction",units="/",ncid=ncid)
    call nc_write(fnm,"f_shelf_old",  lndvc%cell%f_shelf_old,  dims=[dim_lon,dim_lat],long_name="previous-step shelf fraction",units="/",ncid=ncid)
    call nc_write(fnm,"f_lake_old",   lndvc%cell%f_lake_old,   dims=[dim_lon,dim_lat],long_name="previous-step lake fraction",units="/",ncid=ncid)
    call nc_write(fnm,"f_veg_old",    lndvc%cell%f_veg_old,    dims=[dim_lon,dim_lat],long_name="previous-step vegetation fraction",units="/",ncid=ncid)
    call nc_write(fnm,"f_peat",       lndvc%cell%f_peat,       dims=[dim_lon,dim_lat],long_name="peatland fraction",units="/",ncid=ncid)

    ! ---- C. Per-vc descriptor + phys_init ----
    allocate(vari(ni,nj,n_vc))
    do k=1,n_vc; do j=1,nj; do i=1,ni
      vari(i,j,k) = lndvc%vc(i,j,k)%desc%class
    enddo; enddo; enddo
    call nc_write(fnm,"desc_class", vari, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="vc surface class",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      vari(i,j,k) = lndvc%vc(i,j,k)%desc%i
    enddo; enddo; enddo
    call nc_write(fnm,"desc_i",     vari, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="vc coarse-cell i index",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      vari(i,j,k) = lndvc%vc(i,j,k)%desc%j
    enddo; enddo; enddo
    call nc_write(fnm,"desc_j",     vari, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="vc coarse-cell j index",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      vari(i,j,k) = merge(1,0,lndvc%vc(i,j,k)%phys_init)
    enddo; enddo; enddo
    call nc_write(fnm,"phys_init",  vari, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="one-time physical-init guard",units="/",ncid=ncid)
    deallocate(vari)

    allocate(var(ni,nj,n_vc))
    do k=1,n_vc; do j=1,nj; do i=1,ni
      var(i,j,k) = lndvc%vc(i,j,k)%desc%z
    enddo; enddo; enddo
    call nc_write(fnm,"desc_z",         var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="vc band elevation",units="m",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      var(i,j,k) = lndvc%vc(i,j,k)%desc%dz
    enddo; enddo; enddo
    call nc_write(fnm,"desc_dz",        var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="vc band width",units="m",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      var(i,j,k) = lndvc%vc(i,j,k)%desc%w
    enddo; enddo; enddo
    call nc_write(fnm,"desc_w",         var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="vc area weight",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      var(i,j,k) = lndvc%vc(i,j,k)%desc%lat
    enddo; enddo; enddo
    call nc_write(fnm,"desc_lat",       var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="vc latitude",units="deg",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      var(i,j,k) = lndvc%vc(i,j,k)%desc%dz_dx
    enddo; enddo; enddo
    call nc_write(fnm,"desc_dz_dx",     var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="vc surface slope x",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      var(i,j,k) = lndvc%vc(i,j,k)%desc%dz_dy
    enddo; enddo; enddo
    call nc_write(fnm,"desc_dz_dy",     var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="vc surface slope y",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      var(i,j,k) = lndvc%vc(i,j,k)%desc%grad
    enddo; enddo; enddo
    call nc_write(fnm,"desc_grad",      var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="vc |grad z|",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      var(i,j,k) = lndvc%vc(i,j,k)%desc%z_sur_std
    enddo; enddo; enddo
    call nc_write(fnm,"desc_z_sur_std", var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="vc sub-grid elevation std",units="m",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      var(i,j,k) = lndvc%vc(i,j,k)%desc%dz_sur
    enddo; enddo; enddo
    call nc_write(fnm,"desc_dz_sur",    var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="vc sub-grid elevation range for precip",units="m",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      var(i,j,k) = lndvc%vc(i,j,k)%desc%f_ele
    enddo; enddo; enddo
    call nc_write(fnm,"desc_f_ele",     var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="vc elevation-desertification factor",units="/",ncid=ncid)

    ! ---- D. Flux carry-over (per-vc scalars at tile index 1) ----
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%flx%t_skin)) then
        var(i,j,k) = lndvc%vc(i,j,k)%flx%t_skin(1)
      else
        var(i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"t_skin",       var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="skin temperature",units="K",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%flx%t_skin_old)) then
        var(i,j,k) = lndvc%vc(i,j,k)%flx%t_skin_old(1)
      else
        var(i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"t_skin_old",   var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="previous-step skin temperature",units="K",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%flx%albedo)) then
        var(i,j,k) = lndvc%vc(i,j,k)%flx%albedo(1)
      else
        var(i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"albedo",       var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="surface albedo",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%flx%alb_vis_dir)) then
        var(i,j,k) = lndvc%vc(i,j,k)%flx%alb_vis_dir(1)
      else
        var(i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"alb_vis_dir",  var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="alb_vis_dir",units="1",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%flx%alb_vis_dif)) then
        var(i,j,k) = lndvc%vc(i,j,k)%flx%alb_vis_dif(1)
      else
        var(i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"alb_vis_dif",  var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="alb_vis_dif",units="1",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%flx%alb_nir_dir)) then
        var(i,j,k) = lndvc%vc(i,j,k)%flx%alb_nir_dir(1)
      else
        var(i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"alb_nir_dir",  var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="alb_nir_dir",units="1",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%flx%alb_nir_dif)) then
        var(i,j,k) = lndvc%vc(i,j,k)%flx%alb_nir_dif(1)
      else
        var(i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"alb_nir_dif",  var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="alb_nir_dif",units="1",ncid=ncid)
    deallocate(var)

    ! ---- E. Snow block (classes 1/2/3) ----
    allocate(vari(ni,nj,n_vc))
    allocate(var(ni,nj,n_vc))
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; vari(i,j,k) = 0; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%snow)) then
        vari(i,j,k) = lndvc%vc(i,j,k)%snow%mask_snow
      else
        vari(i,j,k) = 0
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"mask_snow",         vari, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="snow mask",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%snow)) then; var(i,j,k) = lndvc%vc(i,j,k)%snow%w_snow; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"w_snow",            var,  dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="snow water equivalent",units="kg/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%snow)) then; var(i,j,k) = lndvc%vc(i,j,k)%snow%w_snow_max; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"w_snow_max",        var,  dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="maximum snow water equivalent",units="kg/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%snow)) then; var(i,j,k) = lndvc%vc(i,j,k)%snow%h_snow; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"h_snow",            var,  dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="snow thickness",units="m",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%snow)) then; var(i,j,k) = lndvc%vc(i,j,k)%snow%f_snow; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"f_snow",            var,  dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="snow-covered fraction",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%snow)) then; var(i,j,k) = lndvc%vc(i,j,k)%snow%snow_grain; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"snow_grain",        var,  dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="snow grain size",units="m",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%snow)) then; var(i,j,k) = lndvc%vc(i,j,k)%snow%dust_con; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"dust_con",          var,  dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="snow dust concentration",units="kg/kg",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%snow)) then; var(i,j,k) = lndvc%vc(i,j,k)%snow%refreezing_sum; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"refreezing_sum",    var,  dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="refreezing-capacity accumulator",units="kg/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%snow)) then; var(i,j,k) = lndvc%vc(i,j,k)%snow%dt_snowfree; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"dt_snowfree",       var,  dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="time since snow-free",units="s",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%snow)) then; var(i,j,k) = lndvc%vc(i,j,k)%snow%alb_snow_vis_dir; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"alb_snow_vis_dir",  var,  dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="snow albedo vis dir",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%snow)) then; var(i,j,k) = lndvc%vc(i,j,k)%snow%alb_snow_vis_dif; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"alb_snow_vis_dif",  var,  dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="snow albedo vis dif",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%snow)) then; var(i,j,k) = lndvc%vc(i,j,k)%snow%alb_snow_nir_dir; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"alb_snow_nir_dir",  var,  dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="snow albedo nir dir",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%snow)) then; var(i,j,k) = lndvc%vc(i,j,k)%snow%alb_snow_nir_dif; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"alb_snow_nir_dif",  var,  dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="snow albedo nir dif",units="/",ncid=ncid)
    deallocate(vari, var)

    ! ---- F. Soil block (class 1: land) ----
    ! F.1 temperature profiles on 0:nl
    allocate(var_n(0:nl,ni,nj,n_vc))
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%soil%t_soil; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"t_soil",     var_n, dims=[dim_depth0,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl+1,ni,nj,n_vc],long_name="soil temperature",units="K",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%soil%t_soil_max; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"t_soil_max", var_n, dims=[dim_depth0,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl+1,ni,nj,n_vc],long_name="max soil temperature",units="K",ncid=ncid)
    deallocate(var_n)

    ! F.2 nl-layer soil hydrology + permafrost
    allocate(var_n(nl,ni,nj,n_vc))
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%soil%theta; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"theta",         var_n, dims=[dim_depth,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl,ni,nj,n_vc],long_name="soil total volumetric water content",units="m3/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%soil%theta_w; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"theta_w",       var_n, dims=[dim_depth,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl,ni,nj,n_vc],long_name="soil liquid volumetric water content",units="m3/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%soil%theta_i; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"theta_i",       var_n, dims=[dim_depth,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl,ni,nj,n_vc],long_name="soil frozen volumetric water content",units="m3/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%soil%w_w; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"w_w",           var_n, dims=[dim_depth,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl,ni,nj,n_vc],long_name="soil liquid water",units="kg/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%soil%w_i; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"w_i",           var_n, dims=[dim_depth,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl,ni,nj,n_vc],long_name="soil frozen water equivalent",units="kg/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%soil%frozen_years; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"frozen_years",  var_n, dims=[dim_depth,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl,ni,nj,n_vc],long_name="permafrost-thaw: consecutive perennially-frozen years",units="yr",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%soil%thaw_timer; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"thaw_timer",    var_n, dims=[dim_depth,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl,ni,nj,n_vc],long_name="permafrost-thaw: priming years remaining",units="yr",ncid=ncid)
    deallocate(var_n)

    ! F.3 soil scalars per-vc
    allocate(var(ni,nj,n_vc))
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var(i,j,k) = lndvc%vc(i,j,k)%soil%alt; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"alt",           var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="active layer thickness",units="m",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var(i,j,k) = lndvc%vc(i,j,k)%soil%w_table; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"w_table",       var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="water table depth",units="m",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var(i,j,k) = lndvc%vc(i,j,k)%soil%w_table_min; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"w_table_min",   var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="minimum yearly water table depth",units="m",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var(i,j,k) = lndvc%vc(i,j,k)%soil%w_table_peat; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"w_table_peat",  var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="peatland water table depth",units="m",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var(i,j,k) = lndvc%vc(i,j,k)%soil%w_table_perch; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"w_table_perch", var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="perched water table depth",units="m",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var(i,j,k) = lndvc%vc(i,j,k)%soil%w_table_eff; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"w_table_eff",   var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="effective water table depth",units="m",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var(i,j,k) = lndvc%vc(i,j,k)%soil%fz_eff; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"fz_eff",        var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="effective frozen depth",units="m",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var(i,j,k) = lndvc%vc(i,j,k)%soil%f_wet; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"f_wet",         var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="saturated fraction",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var(i,j,k) = lndvc%vc(i,j,k)%soil%f_wet_max; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"f_wet_max",     var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="maximum saturated fraction",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var(i,j,k) = lndvc%vc(i,j,k)%soil%f_wetland; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"f_wetland",     var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="TOPMODEL wetland fraction",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var(i,j,k) = lndvc%vc(i,j,k)%soil%cti_lim; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"cti_lim",       var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="wetland CTI threshold",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var(i,j,k) = lndvc%vc(i,j,k)%soil%cti_mean; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"cti_mean",      var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="mean CTI",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var(i,j,k) = lndvc%vc(i,j,k)%soil%mcwd; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"mcwd",          var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="maximum cumulative water deficit",units="mm",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var(i,j,k) = lndvc%vc(i,j,k)%soil%mcwd_clim; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"mcwd_clim",     var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="climatological MCWD",units="mm",ncid=ncid)
    deallocate(var)

    ! ---- G. Vegetation block (class 1: land) ----
    ! G.1 (npft) arrays
    allocate(var_n(npft,ni,nj,n_vc))
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%veg%gdd; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"gdd",       var_n, dims=[dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[npft,ni,nj,n_vc],long_name="growing degree days",units="K",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%veg%phen; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"phen",      var_n, dims=[dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[npft,ni,nj,n_vc],long_name="phenology",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%veg%phen_acc; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"phen_acc",  var_n, dims=[dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[npft,ni,nj,n_vc],long_name="phenology accumulated",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%veg%lai_bal; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"lai_bal",   var_n, dims=[dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[npft,ni,nj,n_vc],long_name="balanced leaf area index",units="m2/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%veg%lai; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"lai",       var_n, dims=[dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[npft,ni,nj,n_vc],long_name="leaf area index",units="m2/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%veg%sai; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"sai",       var_n, dims=[dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[npft,ni,nj,n_vc],long_name="stem area index",units="m2/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%veg%seed_frac; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"seed_frac", var_n, dims=[dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[npft,ni,nj,n_vc],long_name="seed fraction",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%veg%pft_frac; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"pft_frac",  var_n, dims=[dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[npft,ni,nj,n_vc],long_name="vegetation fraction",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%veg%veg_h; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"veg_h",     var_n, dims=[dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[npft,ni,nj,n_vc],long_name="vegetation height",units="m",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%veg%veg_c; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"veg_c",     var_n, dims=[dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[npft,ni,nj,n_vc],long_name="vegetation carbon",units="kgC/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%veg%veg_c13; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"veg_c13",   var_n, dims=[dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[npft,ni,nj,n_vc],long_name="vegetation carbon 13",units="kgC/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%veg%veg_c14; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"veg_c14",   var_n, dims=[dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[npft,ni,nj,n_vc],long_name="vegetation carbon 14",units="kgC/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%veg%leaf_c; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"leaf_c",    var_n, dims=[dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[npft,ni,nj,n_vc],long_name="leaf carbon",units="kgC/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%veg%stem_c; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"stem_c",    var_n, dims=[dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[npft,ni,nj,n_vc],long_name="stem carbon",units="kgC/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%veg%root_c; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"root_c",    var_n, dims=[dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[npft,ni,nj,n_vc],long_name="root carbon",units="kgC/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%veg%npp_ann; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"npp_ann",   var_n, dims=[dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[npft,ni,nj,n_vc],long_name="annual NPP",units="kgC/m2/s",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%veg%npp13_ann; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"npp13_ann", var_n, dims=[dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[npft,ni,nj,n_vc],long_name="annual NPP 13",units="kgC/m2/s",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%veg%npp14_ann; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"npp14_ann", var_n, dims=[dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[npft,ni,nj,n_vc],long_name="annual NPP 14",units="kgC/m2/s",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%veg%gamma_dist; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"gamma_dist",var_n, dims=[dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[npft,ni,nj,n_vc],long_name="vegetation disturbance rate",units="1/s",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%veg%gamma_fire; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"gamma_fire",var_n, dims=[dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[npft,ni,nj,n_vc],long_name="fire disturbance rate",units="1/s",ncid=ncid)
    deallocate(var_n)

    ! G.2 vegetation scalars
    allocate(var(ni,nj,n_vc))
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var(i,j,k) = lndvc%vc(i,j,k)%veg%gdd5; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"gdd5",         var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="growing degree days above 5 degC",units="K",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var(i,j,k) = lndvc%vc(i,j,k)%veg%t2m_min_mon; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"t2m_min_mon",  var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="coldest month temperature",units="K",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%veg)) then; var(i,j,k) = lndvc%vc(i,j,k)%veg%t2m_ann_mean; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"t2m_ann_mean", var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="annual mean 2m temperature",units="K",ncid=ncid)
    deallocate(var)

    ! G.3 root_frac (nl,npft)
    allocate(var_np(nl,npft,ni,nj,n_vc))
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_np(:,:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%soil)) then; var_np(:,:,i,j,k) = lndvc%vc(i,j,k)%soil%root_frac; else; var_np(:,:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"root_frac", var_np, dims=[dim_depth,dim_npft,dim_lon,dim_lat,"vc"],start=[1,1,1,1,1],count=[nl,npft,ni,nj,n_vc],long_name="root fraction in layers",units="/",ncid=ncid)
    deallocate(var_np)

    ! ---- H. Carbon block (class 1: land) ----
    ! H.1 (nlc) arrays
    allocate(var_n(nlc,ni,nj,n_vc))
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%carb%litter_c; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"litter_c",   var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="litter carbon",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%carb%fast_c; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"fast_c",     var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="fast soil carbon",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%carb%slow_c; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"slow_c",     var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="slow soil carbon",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%carb%litter_c13; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"litter_c13", var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="litter carbon 13",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%carb%fast_c13; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"fast_c13",   var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="fast soil carbon 13",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%carb%slow_c13; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"slow_c13",   var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="slow soil carbon 13",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%carb%litter_c14; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"litter_c14", var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="litter carbon 14",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%carb%fast_c14; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"fast_c14",   var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="fast soil carbon 14",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%carb%slow_c14; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"slow_c14",   var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="slow soil carbon 14",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%carb%cato_c; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"cato_c",     var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="catotelm carbon",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%carb%cato_c13; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"cato_c13",   var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="catotelm carbon 13",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%carb%cato_c14; else; var_n(:,i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"cato_c14",   var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="catotelm carbon 14",units="kgC/m3",ncid=ncid)
    deallocate(var_n)

    ! H.2 peat scalars per-vc
    allocate(var(ni,nj,n_vc))
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var(i,j,k) = lndvc%vc(i,j,k)%carb%litter_c_peat; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"litter_c_peat",   var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="peatland litter carbon",units="kgC/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var(i,j,k) = lndvc%vc(i,j,k)%carb%acro_c; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"acro_c",          var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="acrotelm carbon",units="kgC/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var(i,j,k) = lndvc%vc(i,j,k)%carb%litter_c13_peat; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"litter_c13_peat", var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="peatland litter carbon 13",units="kgC/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var(i,j,k) = lndvc%vc(i,j,k)%carb%acro_c13; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"acro_c13",        var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="acrotelm carbon 13",units="kgC/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var(i,j,k) = lndvc%vc(i,j,k)%carb%litter_c14_peat; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"litter_c14_peat", var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="peatland litter carbon 14",units="kgC/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var(i,j,k) = lndvc%vc(i,j,k)%carb%acro_c14; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"acro_c14",        var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="acrotelm carbon 14",units="kgC/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var(i,j,k) = lndvc%vc(i,j,k)%carb%f_peat; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"f_peat_vc",       var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="peatland fraction (per-vc)",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var(i,j,k) = lndvc%vc(i,j,k)%carb%f_peat_pot; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"f_peat_pot",      var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="potential peatland fraction",units="/",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var(i,j,k) = lndvc%vc(i,j,k)%carb%dCpeat_dt; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"dCpeat_dt",       var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="peat accumulation rate",units="kgC/m2/s",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var(i,j,k) = lndvc%vc(i,j,k)%carb%acro_h; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"acro_h",          var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="acrotelm thickness",units="m",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%carb)) then; var(i,j,k) = lndvc%vc(i,j,k)%carb%cato_h; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"cato_h",          var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="catotelm thickness",units="m",ncid=ncid)
    deallocate(var)

    ! ---- I. Ice block (class 3: ice) ----
    ! I.1 firn/skin thermal profile (0:nl_smb)
    allocate(var_n(0:nl_smb,ni,nj,n_vc))
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%t_prof)) then
          var_n(:,i,j,k) = lndvc%vc(i,j,k)%ice%t_prof
        else
          var_n(:,i,j,k) = 0._wp
        endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"t_prof", var_n, dims=["depth_smb0",dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl_smb+1,ni,nj,n_vc],long_name="snow+ice firn thermal profile",units="K",ncid=ncid)
    deallocate(var_n)

    ! I.2 t_ice, t_shelf, t_shelf_max (0:nl)
    allocate(var_n(0:nl,ni,nj,n_vc))
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%t_ice)) then
          var_n(:,i,j,k) = lndvc%vc(i,j,k)%ice%t_ice
        else
          var_n(:,i,j,k) = 0._wp
        endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"t_ice",       var_n, dims=[dim_depth0,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl+1,ni,nj,n_vc],long_name="ice temperature",units="K",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%t_shelf)) then
          var_n(:,i,j,k) = lndvc%vc(i,j,k)%ice%t_shelf
        else
          var_n(:,i,j,k) = 0._wp
        endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"t_shelf",     var_n, dims=[dim_depth0,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl+1,ni,nj,n_vc],long_name="shelf soil temperature",units="K",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%t_shelf_max)) then
          var_n(:,i,j,k) = lndvc%vc(i,j,k)%ice%t_shelf_max
        else
          var_n(:,i,j,k) = 0._wp
        endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"t_shelf_max", var_n, dims=[dim_depth0,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl+1,ni,nj,n_vc],long_name="shelf max soil temperature",units="K",ncid=ncid)
    deallocate(var_n)

    ! I.3 shelf hydrology (nl)
    allocate(var_n(nl,ni,nj,n_vc))
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%theta_w_shelf)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%ice%theta_w_shelf; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"theta_w_shelf", var_n, dims=[dim_depth,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl,ni,nj,n_vc],long_name="shelf liquid volumetric water content",units="m3/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%theta_i_shelf)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%ice%theta_i_shelf; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"theta_i_shelf", var_n, dims=[dim_depth,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl,ni,nj,n_vc],long_name="shelf frozen volumetric water content",units="m3/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%w_w_shelf)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%ice%w_w_shelf; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"w_w_shelf",     var_n, dims=[dim_depth,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl,ni,nj,n_vc],long_name="shelf liquid water",units="kg/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%w_i_shelf)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%ice%w_i_shelf; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"w_i_shelf",     var_n, dims=[dim_depth,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl,ni,nj,n_vc],long_name="shelf frozen water equivalent",units="kg/m2",ncid=ncid)
    deallocate(var_n)

    ! I.4 subglacial + subshelf carbon pools (nlc each) — 9 ice + 9 shelf pools
    allocate(var_n(nlc,ni,nj,n_vc))
    ! ice pools
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%litter_c_ice)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%ice%litter_c_ice; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"litter_c_ice",   var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="litter ice carbon",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%fast_c_ice)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%ice%fast_c_ice; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"fast_c_ice",     var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="fast ice carbon",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%slow_c_ice)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%ice%slow_c_ice; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"slow_c_ice",     var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="slow ice carbon",units="kgC/m3",ncid=ncid)
    ! shelf pools
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%litter_c_shelf)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%ice%litter_c_shelf; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"litter_c_shelf", var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="litter shelf carbon",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%fast_c_shelf)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%ice%fast_c_shelf; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"fast_c_shelf",   var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="fast shelf carbon",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%slow_c_shelf)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%ice%slow_c_shelf; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"slow_c_shelf",   var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="slow shelf carbon",units="kgC/m3",ncid=ncid)
    deallocate(var_n)
    ! NOTE: ice_col_t only carries bulk (litter/fast/slow)_c_(ice|shelf) — no c13/c14
    ! isotope tracers for subglacial/subshelf carbon exist yet in lndvc_def.f90.
    ! The 12 isotope fields requested by the reference pattern are therefore not
    ! written here; restarting an isotope-enabled run will not carry these two
    ! pools' 13C/14C inventories across the restart boundary until those fields
    ! are added to ice_col_t.

    ! ---- J. Lake block (class 2: lake) ----
    ! J.1 temperature profiles
    allocate(var_n(0:nl_l,ni,nj,n_vc))
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%t_lake)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%lake%t_lake; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"t_lake",        var_n, dims=[dim_depth0l,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl_l+1,ni,nj,n_vc],long_name="lake temperature",units="K",ncid=ncid)
    deallocate(var_n)

    allocate(var_n(0:nl,ni,nj,n_vc))
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%t_sublake)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%lake%t_sublake; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"t_sublake",     var_n, dims=[dim_depth0,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl+1,ni,nj,n_vc],long_name="soil temperature below lake",units="K",ncid=ncid)
    deallocate(var_n)

    ! J.2 sublake hydrology (nl)
    allocate(var_n(nl,ni,nj,n_vc))
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%theta_w_sublake)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%lake%theta_w_sublake; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"theta_w_sublake", var_n, dims=[dim_depth,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl,ni,nj,n_vc],long_name="sublake liquid volumetric water content",units="m3/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%theta_i_sublake)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%lake%theta_i_sublake; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"theta_i_sublake", var_n, dims=[dim_depth,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl,ni,nj,n_vc],long_name="sublake frozen volumetric water content",units="m3/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%w_w_sublake)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%lake%w_w_sublake; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"w_w_sublake",     var_n, dims=[dim_depth,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl,ni,nj,n_vc],long_name="sublake liquid water",units="kg/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%w_i_sublake)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%lake%w_i_sublake; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"w_i_sublake",     var_n, dims=[dim_depth,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl,ni,nj,n_vc],long_name="sublake frozen water equivalent",units="kg/m2",ncid=ncid)
    deallocate(var_n)

    ! J.3 lake-ice water fraction (nl_l)
    allocate(var_n(nl_l,ni,nj,n_vc))
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%w_w_lake)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%lake%w_w_lake; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"w_w_lake",  var_n, dims=[dim_depthl,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl_l,ni,nj,n_vc],long_name="lake liquid water",units="kg/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%w_i_lake)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%lake%w_i_lake; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"w_i_lake",  var_n, dims=[dim_depthl,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl_l,ni,nj,n_vc],long_name="lake frozen water",units="kg/m2",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%f_i_lake)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%lake%f_i_lake; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"f_i_lake",  var_n, dims=[dim_depthl,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nl_l,ni,nj,n_vc],long_name="lake frozen water fraction",units="/",ncid=ncid)
    deallocate(var_n)

    ! J.4 lake scalars
    allocate(var(ni,nj,n_vc))
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then; var(i,j,k) = lndvc%vc(i,j,k)%lake%h_lake; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"h_lake",      var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="lake depth",units="m",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then; var(i,j,k) = lndvc%vc(i,j,k)%lake%h_lake_conv; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"h_lake_conv", var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="lake convective mixing depth",units="m",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then; var(i,j,k) = lndvc%vc(i,j,k)%lake%h_lake_mix; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"h_lake_mix",  var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="lake surface mixed-layer depth",units="m",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var(i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then; var(i,j,k) = lndvc%vc(i,j,k)%lake%f_lake_ice; else; var(i,j,k) = 0._wp; endif
    enddo; enddo; enddo
    call nc_write(fnm,"f_lake_ice",  var, dims=[dim_lon,dim_lat,"vc"],start=[1,1,1],count=[ni,nj,n_vc],long_name="lake ice fraction",units="/",ncid=ncid)
    deallocate(var)

    ! J.5 lake sediment carbon (nlc each)
    allocate(var_n(nlc,ni,nj,n_vc))
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%litter_c_lake)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%lake%litter_c_lake; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"litter_c_lake",   var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="litter lake carbon",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%fast_c_lake)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%lake%fast_c_lake; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"fast_c_lake",     var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="fast lake carbon",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%slow_c_lake)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%lake%slow_c_lake; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"slow_c_lake",     var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="slow lake carbon",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%litter_c13_lake)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%lake%litter_c13_lake; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"litter_c13_lake", var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="litter lake carbon 13",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%fast_c13_lake)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%lake%fast_c13_lake; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"fast_c13_lake",   var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="fast lake carbon 13",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%slow_c13_lake)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%lake%slow_c13_lake; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"slow_c13_lake",   var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="slow lake carbon 13",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%litter_c14_lake)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%lake%litter_c14_lake; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"litter_c14_lake", var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="litter lake carbon 14",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%fast_c14_lake)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%lake%fast_c14_lake; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"fast_c14_lake",   var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="fast lake carbon 14",units="kgC/m3",ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class == 0) then; var_n(:,i,j,k) = 0._wp; cycle; endif
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%slow_c14_lake)) then; var_n(:,i,j,k) = lndvc%vc(i,j,k)%lake%slow_c14_lake; else; var_n(:,i,j,k) = 0._wp; endif
      else
        var_n(:,i,j,k) = 0._wp
      endif
    enddo; enddo; enddo
    call nc_write(fnm,"slow_c14_lake",   var_n, dims=[dim_depth1,dim_lon,dim_lat,"vc"],start=[1,1,1,1],count=[nlc,ni,nj,n_vc],long_name="slow lake carbon 14",units="kgC/m3",ncid=ncid)
    deallocate(var_n)

    call nc_close(ncid)

    print *,'wrote lndvc restart file ',fnm

    return

  end subroutine lndvc_write_restart


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  l n d v c _ r e a d _ r e s t a r t
  ! Purpose    :  Read per-vc prognostic state from a NetCDF restart file
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lndvc_read_restart(fnm, lndvc)

    implicit none

    character(len=*),  intent(in)    :: fnm
    type(lndvc_class), intent(inout) :: lndvc

    integer :: ncid
    integer :: i, j, k, n_vc
    integer :: n_vc_restart
    integer,  dimension(:,:,:),   allocatable :: vari
    real(wp), dimension(:,:,:),   allocatable :: var
    real(wp), dimension(:,:,:,:), allocatable :: var_n
    real(wp), dimension(:,:,:,:,:), allocatable :: var_np

    n_vc = lndvc%n_vc

    ! the restart file has to be on the same vertical grid as the model
    call check_restart_levels(fnm, dim_depth , nl , "nl")
    call check_restart_levels(fnm, dim_depth1, nlc, "nl (nlc=nl+1)")
    call check_restart_levels(fnm, dim_depthl, nl_l, "nl_lake")

    n_vc_restart = nc_size(fnm,"vc")
    if (n_vc_restart /= n_vc) then
      print *, 'ERROR: n_vc mismatch in lndvc restart: file=', n_vc_restart, ' model=', n_vc
      stop
    endif

    call nc_open(fnm,ncid,writable=.false.)

    ! ---- A. Global 0-D scalars ----
    call nc_read(fnm,"Cflx_avg", lndvc%glob%Cflx_avg, ncid=ncid)
    call nc_read(fnm,"weath_carb_avg", lndvc%glob%weath_carb_avg, ncid=ncid)
    call nc_read(fnm,"weath_sil_avg", lndvc%glob%weath_sil_avg, ncid=ncid)
    call nc_read(fnm,"landc", lndvc%glob%landc, ncid=ncid)
    call nc_read(fnm,"landc13", lndvc%glob%landc13, ncid=ncid)
    call nc_read(fnm,"landc14", lndvc%glob%landc14, ncid=ncid)
    call nc_read(fnm,"burc", lndvc%glob%burc, ncid=ncid)
    call nc_read(fnm,"burc13", lndvc%glob%burc13, ncid=ncid)
    call nc_read(fnm,"burc14", lndvc%glob%burc14, ncid=ncid)
    call nc_read(fnm,"weath_scale", lndvc%glob%weath_scale, ncid=ncid)
    call nc_read(fnm,"co2", lndvc%glob%co2, ncid=ncid)
    call nc_read(fnm,"c13_c12_atm", lndvc%glob%c13_c12_atm, ncid=ncid)
    call nc_read(fnm,"c14_c_atm", lndvc%glob%c14_c_atm, ncid=ncid)

    ! ---- B. Cell aggregation ----
    call nc_read(fnm,"mask_lnd", lndvc%cell%mask_lnd, ncid=ncid)
    call nc_read(fnm,"f_land", lndvc%cell%f_land, ncid=ncid)
    call nc_read(fnm,"f_land0", lndvc%cell%f_land0, ncid=ncid)
    call nc_read(fnm,"f_ice", lndvc%cell%f_ice, ncid=ncid)
    call nc_read(fnm,"f_ice_grd", lndvc%cell%f_ice_grd, ncid=ncid)
    call nc_read(fnm,"f_shelf", lndvc%cell%f_shelf, ncid=ncid)
    call nc_read(fnm,"f_lake", lndvc%cell%f_lake, ncid=ncid)
    call nc_read(fnm,"f_veg", lndvc%cell%f_veg, ncid=ncid)
    call nc_read(fnm,"f_ice_old", lndvc%cell%f_ice_old, ncid=ncid)
    call nc_read(fnm,"f_ice_grd_old", lndvc%cell%f_ice_grd_old, ncid=ncid)
    call nc_read(fnm,"f_ice_nbr", lndvc%cell%f_ice_nbr, ncid=ncid)
    call nc_read(fnm,"f_shelf_old", lndvc%cell%f_shelf_old, ncid=ncid)
    call nc_read(fnm,"f_lake_old", lndvc%cell%f_lake_old, ncid=ncid)
    call nc_read(fnm,"f_veg_old", lndvc%cell%f_veg_old, ncid=ncid)
    call nc_read(fnm,"f_peat", lndvc%cell%f_peat, ncid=ncid)

    ! ---- C. Descriptor (inline; no allocation guard) ----
    allocate(vari(ni,nj,n_vc))
    call nc_read(fnm,"desc_class", vari, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      lndvc%vc(i,j,k)%desc%class = vari(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"desc_i", vari, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      lndvc%vc(i,j,k)%desc%i = vari(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"desc_j", vari, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      lndvc%vc(i,j,k)%desc%j = vari(i,j,k)
    enddo; enddo; enddo
    deallocate(vari)

    allocate(var(ni,nj,n_vc))
    call nc_read(fnm,"desc_z", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      lndvc%vc(i,j,k)%desc%z = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"desc_dz", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      lndvc%vc(i,j,k)%desc%dz = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"desc_w", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      lndvc%vc(i,j,k)%desc%w = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"desc_lat", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      lndvc%vc(i,j,k)%desc%lat = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"desc_dz_dx", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      lndvc%vc(i,j,k)%desc%dz_dx = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"desc_dz_dy", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      lndvc%vc(i,j,k)%desc%dz_dy = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"desc_grad", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      lndvc%vc(i,j,k)%desc%grad = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"desc_z_sur_std", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      lndvc%vc(i,j,k)%desc%z_sur_std = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"desc_dz_sur", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      lndvc%vc(i,j,k)%desc%dz_sur = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"desc_f_ele", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      lndvc%vc(i,j,k)%desc%f_ele = var(i,j,k)
    enddo; enddo; enddo
    deallocate(var)

    ! ---- D. Flux carry-over (per-vc scalars at tile index 1) ----
    allocate(var(ni,nj,n_vc))
    call nc_read(fnm,"t_skin", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%flx%t_skin)) lndvc%vc(i,j,k)%flx%t_skin(1) = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"t_skin_old", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%flx%t_skin_old)) lndvc%vc(i,j,k)%flx%t_skin_old(1) = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"albedo", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%flx%albedo)) lndvc%vc(i,j,k)%flx%albedo(1) = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"alb_vis_dir", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%flx%alb_vis_dir)) lndvc%vc(i,j,k)%flx%alb_vis_dir(1) = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"alb_vis_dif", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%flx%alb_vis_dif)) lndvc%vc(i,j,k)%flx%alb_vis_dif(1) = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"alb_nir_dir", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%flx%alb_nir_dir)) lndvc%vc(i,j,k)%flx%alb_nir_dir(1) = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"alb_nir_dif", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%flx%alb_nir_dif)) lndvc%vc(i,j,k)%flx%alb_nir_dif(1) = var(i,j,k)
    enddo; enddo; enddo
    deallocate(var)

    ! ---- E. Snow block (classes 1/2/3) ----
    allocate(vari(ni,nj,n_vc))
    call nc_read(fnm,"mask_snow", vari, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%snow)) lndvc%vc(i,j,k)%snow%mask_snow = vari(i,j,k)
    enddo; enddo; enddo
    deallocate(vari)

    allocate(var(ni,nj,n_vc))
    call nc_read(fnm,"w_snow", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%snow)) lndvc%vc(i,j,k)%snow%w_snow = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"w_snow_max", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%snow)) lndvc%vc(i,j,k)%snow%w_snow_max = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"h_snow", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%snow)) lndvc%vc(i,j,k)%snow%h_snow = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"f_snow", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%snow)) lndvc%vc(i,j,k)%snow%f_snow = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"snow_grain", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%snow)) lndvc%vc(i,j,k)%snow%snow_grain = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"dust_con", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%snow)) lndvc%vc(i,j,k)%snow%dust_con = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"refreezing_sum", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%snow)) lndvc%vc(i,j,k)%snow%refreezing_sum = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"dt_snowfree", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%snow)) lndvc%vc(i,j,k)%snow%dt_snowfree = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"alb_snow_vis_dir", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%snow)) lndvc%vc(i,j,k)%snow%alb_snow_vis_dir = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"alb_snow_vis_dif", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%snow)) lndvc%vc(i,j,k)%snow%alb_snow_vis_dif = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"alb_snow_nir_dir", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%snow)) lndvc%vc(i,j,k)%snow%alb_snow_nir_dir = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"alb_snow_nir_dif", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%snow)) lndvc%vc(i,j,k)%snow%alb_snow_nir_dif = var(i,j,k)
    enddo; enddo; enddo
    deallocate(var)

    ! ---- F. Soil block (class 1: land) ----
    ! F.1 temperature profiles on 0:nl
    allocate(var_n(0:nl,ni,nj,n_vc))
    call nc_read(fnm,"t_soil", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%t_soil = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"t_soil_max", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%t_soil_max = var_n(:,i,j,k)
    enddo; enddo; enddo
    deallocate(var_n)

    ! F.2 nl-layer soil hydrology + permafrost
    allocate(var_n(nl,ni,nj,n_vc))
    call nc_read(fnm,"theta", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%theta = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"theta_w", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%theta_w = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"theta_i", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%theta_i = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"w_w", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%w_w = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"w_i", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%w_i = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"frozen_years", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%frozen_years = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"thaw_timer", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%thaw_timer = var_n(:,i,j,k)
    enddo; enddo; enddo
    deallocate(var_n)

    ! F.3 soil scalars per-vc
    allocate(var(ni,nj,n_vc))
    call nc_read(fnm,"alt", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%alt = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"w_table", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%w_table = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"w_table_min", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%w_table_min = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"w_table_peat", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%w_table_peat = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"w_table_perch", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%w_table_perch = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"w_table_eff", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%w_table_eff = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"fz_eff", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%fz_eff = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"f_wet", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%f_wet = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"f_wet_max", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%f_wet_max = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"f_wetland", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%f_wetland = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"cti_lim", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%cti_lim = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"cti_mean", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%cti_mean = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"mcwd", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%mcwd = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"mcwd_clim", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%mcwd_clim = var(i,j,k)
    enddo; enddo; enddo
    deallocate(var)

    ! ---- G. Vegetation block (class 1: land) ----
    ! G.1 (npft) arrays
    allocate(var_n(npft,ni,nj,n_vc))
    call nc_read(fnm,"gdd", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%gdd = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"phen", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%phen = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"phen_acc", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%phen_acc = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"lai_bal", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%lai_bal = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"lai", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%lai = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"sai", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%sai = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"seed_frac", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%seed_frac = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"pft_frac", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%pft_frac = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"veg_h", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%veg_h = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"veg_c", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%veg_c = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"veg_c13", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%veg_c13 = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"veg_c14", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%veg_c14 = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"leaf_c", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%leaf_c = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"stem_c", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%stem_c = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"root_c", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%root_c = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"npp_ann", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%npp_ann = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"npp13_ann", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%npp13_ann = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"npp14_ann", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%npp14_ann = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"gamma_dist", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%gamma_dist = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"gamma_fire", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%gamma_fire = var_n(:,i,j,k)
    enddo; enddo; enddo
    deallocate(var_n)

    ! G.2 vegetation scalars
    allocate(var(ni,nj,n_vc))
    call nc_read(fnm,"gdd5", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%gdd5 = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"t2m_min_mon", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%t2m_min_mon = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"t2m_ann_mean", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%veg)) lndvc%vc(i,j,k)%veg%t2m_ann_mean = var(i,j,k)
    enddo; enddo; enddo
    deallocate(var)

    ! G.3 root_frac (nl,npft)
    allocate(var_np(nl,npft,ni,nj,n_vc))
    call nc_read(fnm,"root_frac", var_np, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%soil)) lndvc%vc(i,j,k)%soil%root_frac = var_np(:,:,i,j,k)
    enddo; enddo; enddo
    deallocate(var_np)

    ! ---- H. Carbon block (class 1: land) ----
    ! H.1 (nlc) arrays
    allocate(var_n(nlc,ni,nj,n_vc))
    call nc_read(fnm,"litter_c", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%litter_c = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"fast_c", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%fast_c = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"slow_c", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%slow_c = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"litter_c13", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%litter_c13 = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"fast_c13", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%fast_c13 = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"slow_c13", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%slow_c13 = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"litter_c14", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%litter_c14 = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"fast_c14", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%fast_c14 = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"slow_c14", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%slow_c14 = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"cato_c", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%cato_c = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"cato_c13", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%cato_c13 = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"cato_c14", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%cato_c14 = var_n(:,i,j,k)
    enddo; enddo; enddo
    deallocate(var_n)

    ! H.2 peat scalars per-vc
    allocate(var(ni,nj,n_vc))
    call nc_read(fnm,"litter_c_peat", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%litter_c_peat = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"acro_c", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%acro_c = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"litter_c13_peat", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%litter_c13_peat = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"acro_c13", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%acro_c13 = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"litter_c14_peat", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%litter_c14_peat = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"acro_c14", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%acro_c14 = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"f_peat_vc", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%f_peat = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"f_peat_pot", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%f_peat_pot = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"dCpeat_dt", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%dCpeat_dt = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"acro_h", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%acro_h = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"cato_h", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%carb)) lndvc%vc(i,j,k)%carb%cato_h = var(i,j,k)
    enddo; enddo; enddo
    deallocate(var)

    ! ---- I. Ice block (class 3: ice) ----
    ! I.1 firn/skin thermal profile (0:nl_smb)
    allocate(var_n(0:nl_smb,ni,nj,n_vc))
    call nc_read(fnm,"t_prof", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%t_prof)) lndvc%vc(i,j,k)%ice%t_prof = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    deallocate(var_n)

    ! I.2 t_ice, t_shelf, t_shelf_max (0:nl)
    allocate(var_n(0:nl,ni,nj,n_vc))
    call nc_read(fnm,"t_ice", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%t_ice)) lndvc%vc(i,j,k)%ice%t_ice = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    call nc_read(fnm,"t_shelf", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%t_shelf)) lndvc%vc(i,j,k)%ice%t_shelf = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    call nc_read(fnm,"t_shelf_max", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%t_shelf_max)) lndvc%vc(i,j,k)%ice%t_shelf_max = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    deallocate(var_n)

    ! I.3 shelf hydrology (nl)
    allocate(var_n(nl,ni,nj,n_vc))
    call nc_read(fnm,"theta_w_shelf", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%theta_w_shelf)) lndvc%vc(i,j,k)%ice%theta_w_shelf = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    call nc_read(fnm,"theta_i_shelf", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%theta_i_shelf)) lndvc%vc(i,j,k)%ice%theta_i_shelf = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    call nc_read(fnm,"w_w_shelf", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%w_w_shelf)) lndvc%vc(i,j,k)%ice%w_w_shelf = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    call nc_read(fnm,"w_i_shelf", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%w_i_shelf)) lndvc%vc(i,j,k)%ice%w_i_shelf = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    deallocate(var_n)

    ! I.4 subglacial + subshelf carbon pools (nlc each)
    allocate(var_n(nlc,ni,nj,n_vc))
    call nc_read(fnm,"litter_c_ice", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%litter_c_ice)) lndvc%vc(i,j,k)%ice%litter_c_ice = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    call nc_read(fnm,"fast_c_ice", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%fast_c_ice)) lndvc%vc(i,j,k)%ice%fast_c_ice = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    call nc_read(fnm,"slow_c_ice", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%slow_c_ice)) lndvc%vc(i,j,k)%ice%slow_c_ice = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    call nc_read(fnm,"litter_c_shelf", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%litter_c_shelf)) lndvc%vc(i,j,k)%ice%litter_c_shelf = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    call nc_read(fnm,"fast_c_shelf", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%fast_c_shelf)) lndvc%vc(i,j,k)%ice%fast_c_shelf = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    call nc_read(fnm,"slow_c_shelf", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%ice)) then
        if (allocated(lndvc%vc(i,j,k)%ice%slow_c_shelf)) lndvc%vc(i,j,k)%ice%slow_c_shelf = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    deallocate(var_n)
    ! NOTE: c13/c14 subglacial+subshelf carbon isotope tracers are not
    ! restored here — see the corresponding NOTE in lndvc_write_restart.

    ! ---- J. Lake block (class 2: lake) ----
    ! J.1 temperature profiles
    allocate(var_n(0:nl_l,ni,nj,n_vc))
    call nc_read(fnm,"t_lake", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%t_lake)) lndvc%vc(i,j,k)%lake%t_lake = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    deallocate(var_n)

    allocate(var_n(0:nl,ni,nj,n_vc))
    call nc_read(fnm,"t_sublake", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%t_sublake)) lndvc%vc(i,j,k)%lake%t_sublake = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    deallocate(var_n)

    ! J.2 sublake hydrology (nl)
    allocate(var_n(nl,ni,nj,n_vc))
    call nc_read(fnm,"theta_w_sublake", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%theta_w_sublake)) lndvc%vc(i,j,k)%lake%theta_w_sublake = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    call nc_read(fnm,"theta_i_sublake", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%theta_i_sublake)) lndvc%vc(i,j,k)%lake%theta_i_sublake = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    call nc_read(fnm,"w_w_sublake", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%w_w_sublake)) lndvc%vc(i,j,k)%lake%w_w_sublake = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    call nc_read(fnm,"w_i_sublake", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%w_i_sublake)) lndvc%vc(i,j,k)%lake%w_i_sublake = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    deallocate(var_n)

    ! J.3 lake-ice water fraction (nl_l)
    allocate(var_n(nl_l,ni,nj,n_vc))
    call nc_read(fnm,"w_w_lake", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%w_w_lake)) lndvc%vc(i,j,k)%lake%w_w_lake = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    call nc_read(fnm,"w_i_lake", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%w_i_lake)) lndvc%vc(i,j,k)%lake%w_i_lake = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    call nc_read(fnm,"f_i_lake", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) then
        if (allocated(lndvc%vc(i,j,k)%lake%f_i_lake)) lndvc%vc(i,j,k)%lake%f_i_lake = var_n(:,i,j,k)
      endif
    enddo; enddo; enddo
    deallocate(var_n)

    ! J.4 lake scalars
    allocate(var(ni,nj,n_vc))
    call nc_read(fnm,"h_lake", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) lndvc%vc(i,j,k)%lake%h_lake = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"h_lake_conv", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) lndvc%vc(i,j,k)%lake%h_lake_conv = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"h_lake_mix", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) lndvc%vc(i,j,k)%lake%h_lake_mix = var(i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"f_lake_ice", var, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) lndvc%vc(i,j,k)%lake%f_lake_ice = var(i,j,k)
    enddo; enddo; enddo
    deallocate(var)

    ! J.5 lake sediment carbon (nlc each)
    allocate(var_n(nlc,ni,nj,n_vc))
    call nc_read(fnm,"litter_c_lake", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) lndvc%vc(i,j,k)%lake%litter_c_lake = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"fast_c_lake", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) lndvc%vc(i,j,k)%lake%fast_c_lake = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"slow_c_lake", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) lndvc%vc(i,j,k)%lake%slow_c_lake = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"litter_c13_lake", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) lndvc%vc(i,j,k)%lake%litter_c13_lake = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"fast_c13_lake", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) lndvc%vc(i,j,k)%lake%fast_c13_lake = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"slow_c13_lake", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) lndvc%vc(i,j,k)%lake%slow_c13_lake = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"litter_c14_lake", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) lndvc%vc(i,j,k)%lake%litter_c14_lake = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"fast_c14_lake", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) lndvc%vc(i,j,k)%lake%fast_c14_lake = var_n(:,i,j,k)
    enddo; enddo; enddo
    call nc_read(fnm,"slow_c14_lake", var_n, ncid=ncid)
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (allocated(lndvc%vc(i,j,k)%lake)) lndvc%vc(i,j,k)%lake%slow_c14_lake = var_n(:,i,j,k)
    enddo; enddo; enddo
    deallocate(var_n)

    ! one-time physical-init guard: a restart is only ever taken once the
    ! per-class physics for every active vc has already run at least once
    do k=1,n_vc; do j=1,nj; do i=1,ni
      if (lndvc%vc(i,j,k)%desc%class /= 0) lndvc%vc(i,j,k)%phys_init = .true.
    enddo; enddo; enddo

    call nc_close(ncid)

    print *,'read lndvc restart file ',fnm

    return

  end subroutine lndvc_read_restart


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  c h e c k _ r e s t a r t _ l e v e l s
  ! Purpose    :  stop with an informative message if the number of vertical
  !               levels in the restart file does not match the model grid
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine check_restart_levels(fnm,dimname,n,parname)

    implicit none

    character (len=*), intent(in) :: fnm      ! restart file name
    character (len=*), intent(in) :: dimname  ! name of the level dimension in the restart file
    integer,           intent(in) :: n        ! number of levels expected by the model
    character (len=*), intent(in) :: parname  ! name of the namelist parameter setting n

    integer :: n_restart


    n_restart = nc_size(fnm,dimname)

    if (n_restart.ne.n) then
      print *,'ERROR: vertical grid mismatch in the lndvc restart file'
      print *,'       file        : ',trim(fnm)
      print *,'       dimension   : ',trim(dimname)
      print *,'       levels in file        : ',n_restart
      print *,'       levels requested by ',trim(parname),' in lndvc_par.nml : ',n
      print *,'       the restart file has to be regenerated for the new vertical grid,'
      print *,'       or the vertical grid in lndvc_par.nml has to be set back to match the restart file'
      stop
    endif

  return

  end subroutine check_restart_levels


end module lndvc_restart_m
