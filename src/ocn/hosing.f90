module hosing_mod

  use precision, only : wp, sp
  use constants, only : pi
  use control, only : out_Dir
  use timer, only: sec_year, time_soy_ocn, year_ini
  use climber_grid, only : lat, lon, basin_mask, i_atlantic, i_pacific, i_southern
  use ocn_grid, only : maxi, maxj, dx, dy
  use ocn_params, only: dt, rho0
  use ocn_params, only : n_hosing_domain, n_hosing_domain_max, hosing_domain_name, hosing_comp_basin
  use ncio
  use nml

  implicit none

  real(wp), allocatable, dimension(:,:) :: hosing_comp_mask

  type hosing_type
    real(wp), allocatable, dimension(:) :: time
    real(wp), allocatable, dimension(:) :: fwf
    real(wp), allocatable, dimension(:,:) :: mask
    real(wp) :: area
  end type
  type(hosing_type) :: hosing(n_hosing_domain_max)

  private
  public :: hosing_init, hosing_update

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  h o s i n g _ u p d a t e
  !   Purpose    :  update freshwater forcing
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine hosing_update(time,f_ocn,fw_hosing,fw_hosing_comp) 

    implicit none

    real(wp), intent(in) :: time                          ! [yr] Current time, years before 2000 AD
    real(wp), dimension(:,:), intent(in) :: f_ocn
    real(wp), dimension(:,:), intent(inout) :: fw_hosing    ! kg/m2/s
    real(wp), dimension(:,:), intent(inout) :: fw_hosing_comp    ! kg/m2/s

    integer :: n, i0, i1, imin
    real(wp) :: w0, w1
    real(wp) :: fwf_now
    real(wp) :: fw_hosing_tot
    real(wp) :: hosing_comp_area


    fw_hosing(:,:) = 0._wp
    fw_hosing_comp(:,:) = 0._wp
    fw_hosing_tot = 0._wp

    ! combine freshwater hosing flux from the different hosing domains

    do n=1,n_hosing_domain

      ! compute area of hosing domain
      call hosing_domain_area(f_ocn, hosing(n)%mask, hosing(n)%area)

      ! interpolate hosing rate to current year
      if (time.le.hosing(n)%time(1)) then
        fwf_now = hosing(n)%fwf(1)
      elseif (time.ge.hosing(n)%time(ubound(hosing(n)%time,1))) then
        fwf_now = hosing(n)%fwf(ubound(hosing(n)%fwf,1))
      else
        imin = minloc(abs(hosing(n)%time-time),1) 
        if (hosing(n)%time(imin).lt.time) then
          i0 = imin
          i1 = imin+1
        else
          i0 = imin-1
          i1 = imin
        endif
        if (hosing(n)%time(i1)-hosing(n)%time(i0).eq.0._wp) then
          w0 = 1._wp
        else
          w0 = 1._wp - abs(hosing(n)%time(i0)-time)/(hosing(n)%time(i1)-hosing(n)%time(i0))
        endif
        w1 = 1._wp - w0
        fwf_now = w0*hosing(n)%fwf(i0) + w1*hosing(n)%fwf(i1)
      endif

      ! apply freshwater flux
      where (f_ocn.gt.0._wp .and. hosing(n)%mask.gt.0._wp)
        fw_hosing(:,:) = fw_hosing(:,:) + hosing(n)%mask(:,:)*fwf_now*1.e6_wp * rho0 / hosing(n)%area ! m3/s * kg/m3 / m2 = kg/m2/s
      endwhere

      fw_hosing_tot = fw_hosing_tot + fwf_now   ! Sv

    enddo

    ! apply compensation flux if needed

    if (hosing_comp_basin.ge.0) then 

      ! compute area of hosing compensation domain
      call hosing_domain_area(f_ocn, hosing_comp_mask, hosing_comp_area)

      where (f_ocn.gt.0._wp .and. hosing_comp_mask.gt.0._wp)
        fw_hosing_comp(:,:) = fw_hosing_comp(:,:) - fw_hosing_tot*1.e6_wp * rho0 / hosing_comp_area ! m3/s * kg/m3 / m2 = kg/m2/s
      endwhere

    endif

   return

  end subroutine hosing_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  h o s i n g _ i n i t
  !   Purpose    :  initialize freshwater forcing
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine hosing_init

    implicit none

    integer :: n
  

    ! initialize hosing domains
    do n=1,n_hosing_domain
      call hosing_domain_init(hosing_domain_name(n), hosing(n)%time, hosing(n)%fwf, hosing(n)%mask)
    enddo

    ! initialize hosing compensation domain
    if (hosing_comp_basin.ge.0) then 
      allocate(hosing_comp_mask(maxi,maxj))
      call hosing_comp_domain_init(hosing_comp_mask)
      if (hosing_comp_basin.eq.0) then
        ! global compensation, but exclude hosing domains
        do n=1,n_hosing_domain
          where (hosing(n)%mask.gt.0._wp)
            hosing_comp_mask = 0._wp
          endwhere
        enddo
      endif
    endif


    return

  end subroutine hosing_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  h o s i n g _ d o m a i n _ i n i t
  !   Purpose    :  initialize hosing domains
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine hosing_domain_init(domain,time,fwf,mask)

    implicit none

    character(len=*), intent(in)  :: domain
    real(wp), intent(inout), allocatable :: time(:)
    real(wp), intent(inout), allocatable :: fwf(:)
    real(wp), intent(inout), allocatable :: mask(:,:)

    integer :: i, j
    integer :: ntime
    character(len=256) :: fnm
    logical :: l_hosing_rate_from_file, l_hosing_mask_from_file

    real(wp) :: hosing_ini, hosing_trend
    integer :: year_hosing_ini
    integer :: year_hosing_end
    integer :: hosing_basin
    real(wp) :: lat_min_hosing, lat_max_hosing
    real(wp) :: lon_min_hosing, lon_max_hosing

    character(len=256) :: hosing_rate_file
    character(len=256) :: hosing_mask_file


    fnm = trim(out_dir)//"/hosing.nml"

    call nml_read(fnm,domain,"l_hosing_rate_from_file",l_hosing_rate_from_file)
    call nml_read(fnm,domain,"l_hosing_mask_from_file",l_hosing_mask_from_file)

    !---- Hosing rate ----

    if (l_hosing_rate_from_file) then
 
      call nml_read(fnm,domain,"hosing_rate_file",hosing_rate_file)
      ntime = nc_size(trim(hosing_rate_file),"time")
      allocate( time(ntime) )
      allocate( fwf(ntime) )
      call nc_read(trim(hosing_rate_file),"time",time)    
      call nc_read(trim(hosing_rate_file),"fwf",fwf) 

    else

      call nml_read(fnm,domain,"hosing_ini      ",hosing_ini      )
      call nml_read(fnm,domain,"hosing_trend    ",hosing_trend    )
      call nml_read(fnm,domain,"year_hosing_ini ",year_hosing_ini )
      call nml_read(fnm,domain,"year_hosing_end ",year_hosing_end )
      hosing_trend = hosing_trend / 1000._wp    ! convert from Sv/kyr to Sv/yr

      allocate( time(5) )
      allocate( fwf(5) )
      time(1) = year_ini - 1._wp
      fwf(1) = 0._wp
      time(2) = max(year_ini, year_hosing_ini-1)
      fwf(2) = 0._wp
      time(3) = max(year_ini, year_hosing_ini+1)
      fwf(3) = hosing_ini
      time(4) = year_hosing_end
      fwf(4) = hosing_ini + hosing_trend*(year_hosing_end-year_hosing_ini)
      time(5) = year_hosing_end + 1._wp
      fwf(5) = 0._wp

    endif

! TODO, allow for 3rd option to prescribed hosing rate from code?
!      if (time.lt.year_hosing_ini) then
!        hosing = 0._wp
!      else if (time.ge.year_hosing_ini .and. time.lt.year_hosing_end) then
!        hosing = hosing_ini
!      else if (time.ge.year_hosing_end .and. time.lt.(year_hosing_end+(year_hosing_end-year_hosing_ini))) then
!        hosing = -hosing_ini
!      else
!        hosing = 0._wp
!      endif
!

    !---- Hosing mask ----

    allocate(mask(maxi,maxj))

    if (l_hosing_mask_from_file) then
 
      call nml_read(fnm,domain,"hosing_mask_file",hosing_mask_file)
      call nc_read(trim(hosing_mask_file),"mask",mask) 

    else

      call nml_read(fnm,domain,"hosing_basin",hosing_basin)
      call nml_read(fnm,domain,"lat_min_hosing",lat_min_hosing)
      call nml_read(fnm,domain,"lat_max_hosing",lat_max_hosing)
      call nml_read(fnm,domain,"lon_min_hosing",lon_min_hosing)
      call nml_read(fnm,domain,"lon_max_hosing",lon_max_hosing)

      ! derive mask 
      do j=1,maxj
        do i=1,maxi
          if (lon(i).ge.lon_min_hosing .and. lon(i).le.lon_max_hosing .and. lat(j).ge.lat_min_hosing .and. lat(j).le.lat_max_hosing) then
            if ((hosing_basin.eq.1 .and. basin_mask(i,j).eq.i_atlantic) .or. &
              (hosing_basin.eq.2 .and. basin_mask(i,j).eq.i_pacific) .or. &
              (hosing_basin.eq.3 .and. basin_mask(i,j).eq.i_southern)) then
              mask(i,j) = 1._wp
            else
              mask(i,j) = 0._wp
            endif
          endif
        enddo
      enddo

    endif

    return

  end subroutine hosing_domain_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  h o s i n g _ c o m p _ d o m a i n _ i n i t
  !   Purpose    :  initialize hosing compensation domain
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine hosing_comp_domain_init(mask)

    implicit none

    real(wp), intent(out) :: mask(:,:)

    integer :: i, j
    character(len=256) :: fnm

    real(wp) :: lat_min_hosing_comp, lat_max_hosing_comp
    real(wp) :: lon_min_hosing_comp, lon_max_hosing_comp

    fnm = trim(out_dir)//"/ocn_par.nml"
    call nml_read(fnm,"ocn_par","lat_min_hosing_comp",lat_min_hosing_comp)
    call nml_read(fnm,"ocn_par","lat_max_hosing_comp",lat_max_hosing_comp)
    call nml_read(fnm,"ocn_par","lon_min_hosing_comp",lon_min_hosing_comp)
    call nml_read(fnm,"ocn_par","lon_max_hosing_comp",lon_max_hosing_comp)

    ! derive mask 
    if (hosing_comp_basin.eq.0) then
      ! global compensation

      mask(:,:) = 1._wp

    else
      ! compensation over specific region

      do j=1,maxj
        do i=1,maxi
          if (lon(i).ge.lon_min_hosing_comp .and. lon(i).le.lon_max_hosing_comp .and. lat(j).ge.lat_min_hosing_comp .and. lat(j).le.lat_max_hosing_comp) then
            if ((hosing_comp_basin.eq.1 .and. basin_mask(i,j).eq.i_atlantic) .or. &
              (hosing_comp_basin.eq.2 .and. basin_mask(i,j).eq.i_pacific) .or. &
              (hosing_comp_basin.eq.3 .and. basin_mask(i,j).eq.i_southern)) then
              mask(i,j) = 1._wp
            else
              mask(i,j) = 0._wp
            endif
          endif
        enddo
      enddo

    endif

    return

  end subroutine hosing_comp_domain_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  h o s i n g _ d o m a i n _ a r e a 
  !   Purpose    :  compute area of hosing domain
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine hosing_domain_area(f_ocn, mask, area)

    implicit none

    real(wp), intent(in) :: f_ocn(:,:)
    real(wp), intent(in) :: mask(:,:)
    real(wp), intent(out) :: area

    integer :: i, j

    area = 0._wp
    do j=1,maxj
      do i=1,maxi
        area = area + mask(i,j) * f_ocn(i,j)*dx(j)*dy   ! m2
      enddo
    enddo

    return

  end subroutine hosing_domain_area


end module hosing_mod
