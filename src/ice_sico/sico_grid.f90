module sico_grid_mod

  use nml
  use control, only : out_dir
  use coords, only : grid_class, grid_init, grid_write

  use sico_types_m
  use sico_params, only : sico_par_class
  use sico_params, only : pi, RHO, G, KAPPA_R, RHO_C_R, BETA, L, NUE
  use sico_params, only : eps, eps_dp, pi_180
  use ice_material_properties_m, only : kappa_val
  use sico_timer, only : sico_timer_class
  use stereo_proj_m

  implicit none


   real(wp) :: H_R = 2.0d+03 !       Thickness of the modelled lithosphere layer = 2000 m 

   ! WGS84 ellipsoid
   real(wp), parameter :: R = 6.371221E6_wp    !       Radius of the Earth = 6371221 m
   real(wp), parameter :: A = 6.378137d+06     !       Semi-major axis of the Earth = 6378137 m
   real(wp), parameter :: B = 6.3567523142d+06 !       Semi-minor axis of the Earth = 6356752.3142 m

   type sico_grid_class
     type(grid_class) :: grid1
     integer :: IMAX
     integer :: JMAX
     integer :: KCMAX
     integer :: KTMAX
     integer :: KRMAX
     real(wp) :: X0
     real(wp) :: Y0
     real(wp) :: DX
     real(wp) :: LAMBDA0
     real(wp) :: PHI0
     real(wp) :: H_R
     real(wp) :: dxi, deta
     real(wp) :: dxi_inv, deta_inv
     real(wp) :: dxi12_inv, deta12_inv
     real(wp) :: dzeta_c, dzeta_t, dzeta_r
     real(wp), dimension(:), allocatable :: xi !! xi(i): Coordinate xi (= x) of grid point i
     real(wp), dimension(:), allocatable :: eta !! eta(j): Coordinate eta (= y) of grid point j
     real(wp), dimension(:), allocatable :: zeta_c !! zeta_c(kc): Sigma coordinate zeta_c of grid point kc
     real(wp), dimension(:), allocatable :: zeta_t !! zeta_t(kt): Sigma coordinate zeta_t of grid point kt
     real(wp), dimension(:), allocatable :: zeta_r !! zeta_r(kr): Sigma coordinate zeta_r of grid point kr
     real(wp) :: aa !! aa: Exponential stretch parameter of the non-equidistant vertical grid in the upper (kc) ice domain
     logical :: flag_aa_nonzero !! flag_aa_nonzero: Flag for the exponential stretch parameter aa. .true.: aa greater than zero (non-equidistant grid) .false.: aa equal to zero (equidistant grid)
     real(wp) :: ea !! ea: Abbreviation for exp(aa)
     real(wp), dimension(:), allocatable :: eaz_c !! eaz_c(kc): Abbreviation for exp(aa*zeta(kc))
     real(wp), dimension(:), allocatable :: eaz_c_quotient !! eaz_c_quotient(kc): Abbreviation for (eaz_c(kc)-1.0)/(ea-1.0)

     real(wp), dimension(:,:), allocatable :: dist_dxdy !! dist_dxdy(jr,ir): Distance between grid points with delta_i=ir, delta_j=jr

     real(wp), dimension(:,:), allocatable :: lambda !! lambda(j,i): Geographic longitude of grid point (i,j)
     real(wp), dimension(:,:), allocatable :: phi !! phi(j,i): Geographic latitude of grid point (i,j)

     real(wp), dimension(:,:), allocatable :: area !! area(j,i): Area of grid cell associated with grid point (i,j)
     real(wp), dimension(:,:), allocatable :: dt_darea 
     real(wp), dimension(:,:), allocatable :: g11_g 
     real(wp), dimension(:,:), allocatable :: g22_g 
     real(wp), dimension(:,:), allocatable :: sq_g11_g !! sq_g11_g(j,i): Square root of the coefficient g11 of the metric tensor on grid point (i,j)
     real(wp), dimension(:,:), allocatable :: sq_g22_g !! sq_g22_g(j,i): Square root of the coefficient g22 of the metric tensor on grid point (i,j)
     real(wp), dimension(:,:), allocatable :: insq_g11_g !! insq_g11_g(j,i): Inverse square root of g11 on grid point (i,j)
     real(wp), dimension(:,:), allocatable :: insq_g22_g !! insq_g22_g(j,i): Inverse square root of g22 on grid point (i,j)
     real(wp), dimension(:,:), allocatable :: g11_sgx 
     real(wp), dimension(:,:), allocatable :: g11_sgy 
     real(wp), dimension(:,:), allocatable :: g22_sgx 
     real(wp), dimension(:,:), allocatable :: g22_sgy 
     real(wp), dimension(:,:), allocatable :: sq_g11_sgx !! sq_g11_sgx(j,i): Square root of g11, at (i+1/2,j)
     real(wp), dimension(:,:), allocatable :: sq_g11_sgy !! sq_g11_sgy(j,i): Square root of g11, at (i,j+1/2)
     real(wp), dimension(:,:), allocatable :: sq_g22_sgx !! sq_g22_sgx(j,i): Square root of g22, at (i+1/2,j)
     real(wp), dimension(:,:), allocatable :: sq_g22_sgy !! sq_g22_sgy(j,i): Square root of g22, at (i,j+1/2)
     real(wp), dimension(:,:), allocatable :: insq_g11_sgx !! insq_g11_sgx(j,i): Inverse square root of g11, at (i+1/2,j)
     real(wp), dimension(:,:), allocatable :: insq_g22_sgy !! insq_g22_sgy(j,i): Inverse square root of g22, at (i,j+1/2)
     logical, dimension(:,:), allocatable :: flag_inner_point !! flag_inner_point(j,i): Inner-point flag. .true.: inner point, .false.: margin point
     real(wp), dimension(:,:), allocatable :: sq_g22_x_1, sq_g22_x_2
     real(wp), dimension(:,:), allocatable :: sq_g11_y_1, sq_g11_y_2

     integer, dimension(:), allocatable :: ii
     integer, dimension(:), allocatable :: jj
     integer, dimension(:,:), allocatable :: nn

     integer, dimension(:), allocatable :: n2i
     integer, dimension(:), allocatable :: n2j
     integer, dimension(:,:), allocatable :: ij2n

     real(wp), dimension(:), allocatable :: zeta_c_sgz, eaz_c_sgz, eaz_c_quotient_sgz
     real(wp), dimension(:,:), allocatable :: fact_x, fact_y
     real(wp), dimension(:), allocatable       :: fact_z_c
     real(wp)                           :: fact_z_t

     real(wp), dimension(:), allocatable :: avxy3
     real(wp), dimension(:), allocatable :: aqxy1

     real(wp), dimension(:), allocatable :: avz3

     real(wp) :: azs2
     real(wp) :: azs3

     real(wp) :: aqbm1 
     real(wp) :: aqbm2 
     real(wp) :: aqbm3a
     real(wp) :: aqbm3b
     real(wp) :: aqbm4 


     real(wp), dimension(:), allocatable :: atm1
     real(wp), dimension(:), allocatable :: atm2

     real(wp) :: at7 
     real(wp) :: aw1
     real(wp) :: aw2
     real(wp) :: aw3
     real(wp) :: aw4
     real(wp) :: aw5
     real(wp) :: aw7
     real(wp) :: aw8
     real(wp) :: aw9 
     real(wp) :: ai3
     real(wp) :: atr1 
     real(wp) :: am1
     real(wp) :: am2
     real(wp) :: acb1 
     real(wp) :: acb2
     real(wp) :: acb3
     real(wp) :: acb4
     real(wp) :: alb1
     real(wp) :: aqtld
     real(wp) :: dtt_2dxi 
     real(wp) :: dtt_2deta

     real(wp), dimension(:), allocatable :: at1
     real(wp), dimension(:), allocatable :: at2_1
     real(wp), dimension(:), allocatable :: at2_2
     real(wp), dimension(:), allocatable :: at3_1
     real(wp), dimension(:), allocatable :: at3_2
     real(wp), dimension(:), allocatable :: at4_1
     real(wp), dimension(:), allocatable :: at4_2
     real(wp), dimension(:), allocatable :: at5
     real(wp), dimension(:), allocatable :: at6
     real(wp), dimension(:), allocatable :: ai1
     real(wp), dimension(:), allocatable :: ai2
     real(wp), dimension(:), allocatable :: aqtlde
     real(wp), dimension(:), allocatable :: am3

   end type

contains

  subroutine sico_grid_init(grid, grd,tmr,par)

    type(grid_class) :: grid
    type(sico_grid_class) :: grd
    type(sico_timer_class) :: tmr
    type(sico_par_class) :: par

    integer :: i, j, m, n, kc, kt, kr, ir, jr
    real(wp) :: K
    character (len=256) :: filename


    ! Read parameters from file
    write(*,*) "ice parameters ==========="
    filename = trim(out_dir)//"/ice_sico_par.nml"
    call nml_read(filename,"ice_sico_par","nkc",grd%KCMAX)
    call nml_read(filename,"ice_sico_par","nkt",grd%KTMAX)
    call nml_read(filename,"ice_sico_par","nkr",grd%KRMAX)

    ! check for oblique stereographic projection
    if (trim(grid%mtype)=="stereographic" .and. par%grid.ne.0) then
      stop 'oblique stereographic projection only supported for grid==0'
    endif

grd%DX = grid%G%dx ! km
grd%X0 = grid%G%x0 ! km
grd%Y0 = grid%G%y0 ! km
grd%IMAX = grid%G%nx-1 ! index starts from 0
grd%JMAX = grid%G%ny-1 ! index starts from 0
grd%PHI0 = grid%proj%phi*pi_180 ! deg -> rad
grd%LAMBDA0 = grid%proj%lambda*pi_180 ! deg -> rad

grd%H_R = H_R

!-------- Compatibility check of the thermodynamics mode
!         (cold vs. polythermal vs. enthalpy method)
!         and the number of grid points in the lower (kt) ice domain --------
if (par%calcmod==0 .or. par%calcmod==2 .or. par%calcmod==3 .or. par%calcmod==-1) then
  if (grd%KTMAX > 2) then
    write(6, fmt='(a)') ' >>> sico_init: for options calcmod==0, 2, 3 or -1,'
    write(6, fmt='(a)') '                the separate kt domain is redundant.'
    write(6, fmt='(a)') '                therefore, consider setting ktmax to 2.'
    write(6, fmt='(a)') ' '
  end if
endif

!-------- Check whether for the shallow shelf
!               or shelfy stream approximation
!                  the chosen grid is Cartesian coordinates
!                             without distortion correction (GRID==0) --------


call sico_grid_alloc(grd)

grd%dxi  = grd%DX *1000.0_wp   ! km -> m
grd%deta = grd%DX *1000.0_wp   ! km -> m

grd%dxi_inv  = 1.0_wp/grd%dxi
grd%deta_inv = 1.0_wp/grd%deta

grd%dxi12_inv  = 1.0_wp/(12.0_wp*grd%dxi)
grd%deta12_inv = 1.0_wp/(12.0_wp*grd%deta)

grd%xi(0:grd%IMAX) = grd%grid1%G%x(1:grd%IMAX+1)*1000._wp ! m
grd%eta(0:grd%JMAX) = grd%grid1%G%y(1:grd%JMAX+1)*1000._wp ! m

!-------- Geographic coordinates, metric tensor,
!                                 gradients of the topography --------

do i=0, grd%IMAX
  do j=0, grd%JMAX

    if (par%grid==0 .or. par%grid==1) then  ! Stereographic projection
      call stereo_inv_ellipsoid(grd%xi(i), grd%eta(j), A, B, &
        grd%LAMBDA0, grd%PHI0, grd%lambda(j,i), grd%phi(j,i))
      !print *,'old lambda,phi', grd%lambda(j,i), grd%phi(j,i)
      grd%lambda(j,i) = grd%grid1%lon(i+1,j+1)*pi_180 ! deg to rad
      grd%phi(j,i)    = grd%grid1%lat(i+1,j+1)*pi_180 ! deg to rad
      !print *,'new lambda,phi', grd%lambda(j,i), grd%phi(j,i)
    else if (par%grid==2) then ! Geographic coordinates
      grd%lambda(j,i) = grd%xi(i)
      grd%phi(j,i)    = grd%eta(j)
    endif

  end do
end do

!-------- Further initializations --------
grd%dzeta_c = 1.0_wp/real(grd%KCMAX,dp)
grd%dzeta_t = 1.0_wp/real(grd%KTMAX,dp)
grd%dzeta_r = 1.0_wp/real(grd%KRMAX,dp)

!-------- General abbreviations --------

!  ------ kc domain

if (par%deform >= eps) then

   grd%flag_aa_nonzero = .true.   ! non-equidistant grid

   grd%aa = par%deform
   grd%ea = exp(grd%aa)

   kc=0
   grd%zeta_c(kc)         = 0.0_wp
   grd%eaz_c(kc)          = 1.0_wp
   grd%eaz_c_quotient(kc) = 0.0_wp

   do kc=1, grd%KCMAX-1
      grd%zeta_c(kc) = kc*grd%dzeta_c
      grd%eaz_c(kc)  = exp(grd%aa*grd%zeta_c(kc))
      grd%eaz_c_quotient(kc) = (grd%eaz_c(kc)-1.0_wp)/(grd%ea-1.0_wp)
   end do

   kc=grd%KCMAX
   grd%zeta_c(kc)         = 1.0_wp
   grd%eaz_c(kc)          = exp(grd%aa)
   grd%eaz_c_quotient(kc) = 1.0_wp

else

   grd%flag_aa_nonzero = .false.   ! equidistant grid

   grd%aa = 0.0_wp
   grd%ea = 1.0_wp

   kc=0
   grd%zeta_c(kc)         = 0.0_wp
   grd%eaz_c(kc)          = 1.0_wp
   grd%eaz_c_quotient(kc) = 0.0_wp

   do kc=1, grd%KCMAX-1
      grd%zeta_c(kc) = kc*grd%dzeta_c
      grd%eaz_c(kc)  = 1.0_wp
      grd%eaz_c_quotient(kc) = grd%zeta_c(kc)
   end do

   kc=grd%KCMAX
   grd%zeta_c(kc)         = 1.0_wp
   grd%eaz_c(kc)          = 1.0_wp
   grd%eaz_c_quotient(kc) = 1.0_wp

end if

!  ------ kt domain

kt=0
grd%zeta_t(kt) = 0.0_wp

do kt=1, grd%KTMAX-1
   grd%zeta_t(kt) = kt*grd%dzeta_t
end do

kt=grd%KTMAX
grd%zeta_t(kt) = 1.0_wp

!  ------ kr domain

kr=0
grd%zeta_r(kr) = 0.0_wp

do kr=1, grd%KRMAX-1
   grd%zeta_r(kr) = kr*grd%dzeta_r
end do

kr=grd%KRMAX
grd%zeta_r(kr) = 1.0_wp

!-------- Construction of a vector (with index n) from a 2-d array
!         (with indices i, j) by diagonal numbering --------

n=1

do m=0, grd%IMAX+grd%JMAX
   do i=m, 0, -1
      j = m-i
      if ((i <= grd%IMAX).and.(j <= grd%JMAX)) then
         grd%ii(n)   = i
         grd%jj(n)   = j
         grd%nn(j,i) = n
         n=n+1
      end if
   end do
end do

!-------- Reshaping of a 2-d array (with indices i, j)
!                                  to a vector (with index n) --------

n=1

do i=0, grd%IMAX
  do j=0, grd%JMAX
    grd%n2i(n)    = i
    grd%n2j(n)    = j
    grd%ij2n(j,i) = n
    n=n+1
  end do
end do


!-------- Distance between grid points with delta_i=ir, delta_j=jr --------

if (par%grid==0 .or. par%grid==1) then  ! Stereographic projection

  do ir=-grd%IMAX, grd%IMAX
    do jr=-grd%JMAX, grd%JMAX
      grd%dist_dxdy(jr,ir) = sqrt( (real(ir,dp)*grd%dxi)**2 + (real(jr,dp)*grd%deta)**2 )
      ! distortion due to stereographic projection not accounted for
    end do
  end do

endif


if (par%grid==0) then  ! Stereographic projection (distortion neglected)

!-------- Components g11, g22 on the grid points (_g) and between
!         the grid points (_sg) --------

  grd%g11_g   = 1.0_wp
  grd%g22_g   = 1.0_wp
  grd%g11_sgx = 1.0_wp
  grd%g11_sgy = 1.0_wp
  grd%g22_sgx = 1.0_wp
  grd%g22_sgy = 1.0_wp

else if (par%grid==1) then  ! Stereographic projection

  if (grd%PHI0 > eps_dp) then   ! for northern hemisphere
     K = (cos(0.25_wp*pi-0.5_wp*grd%PHI0))**2
  else if (grd%PHI0 < (-eps_dp)) then   ! for southern hemisphere
     K = (cos(0.25_wp*pi+0.5_wp*grd%PHI0))**2
  else
     stop ' >>> metric: PHI0 must be different from zero!'
  end if

!-------- Components g11, g22 on the grid points (_g) --------

  do i=0, grd%IMAX
  do j=0, grd%JMAX
     call metric_stereo(grd%xi(i), grd%eta(j), K, grd%g11_g(j,i), grd%g22_g(j,i))
  end do
  end do

!-------- Components g11, g22 between the grid points (_sg) --------

  do i=0, grd%IMAX-1
  do j=0, grd%JMAX
     call metric_stereo(0.5_wp*(grd%xi(i)+grd%xi(i+1)), grd%eta(j), K, &
                        grd%g11_sgx(j,i), grd%g22_sgx(j,i))
  end do
  end do

  do i=0, grd%IMAX
  do j=0, grd%JMAX-1
     call metric_stereo(grd%xi(i), 0.5_wp*(grd%eta(j)+grd%eta(j+1)), K, &
                        grd%g11_sgy(j,i), grd%g22_sgy(j,i))
  end do
  end do

else if (par%grid==2) then !  /* Geographical coordinates */

!-------- Components g11, g22 on the grid points (_g) --------

  do i=0, grd%IMAX
  do j=0, grd%JMAX
     call metric_geogr(grd%eta(j), grd%g11_g(j,i), grd%g22_g(j,i))
  end do
  end do

!-------- Components g11, g22 between the grid points (_sg) --------

  do i=0, grd%IMAX-1
  do j=0, grd%JMAX
     call metric_geogr(grd%eta(j), grd%g11_sgx(j,i), grd%g22_sgx(j,i))
  end do
  end do

  do i=0, grd%IMAX
  do j=0, grd%JMAX-1
     call metric_geogr(0.5_wp*(grd%eta(j)+grd%eta(j+1)), grd%g11_sgy(j,i), grd%g22_sgy(j,i))
  end do
  end do

endif

!-------- Square roots (sq_) and inverse square roots (insq_) of
!         g11 and g22 --------

  do i=0, grd%IMAX
  do j=0, grd%JMAX
     grd%sq_g11_g(j,i)   = sqrt(grd%g11_g(j,i))
     grd%sq_g22_g(j,i)   = sqrt(grd%g22_g(j,i))
     grd%insq_g11_g(j,i) = 1.0_wp/grd%sq_g11_g(j,i)
     grd%insq_g22_g(j,i) = 1.0_wp/grd%sq_g22_g(j,i)
  end do
  end do

  do i=0, grd%IMAX-1
  do j=0, grd%JMAX
     grd%sq_g11_sgx(j,i)   = sqrt(grd%g11_sgx(j,i))
     grd%sq_g22_sgx(j,i)   = sqrt(grd%g22_sgx(j,i))
     grd%insq_g11_sgx(j,i) = 1.0_wp/grd%sq_g11_sgx(j,i)
  end do
  end do

  do i=0, grd%IMAX
  do j=0, grd%JMAX-1
     grd%sq_g11_sgy(j,i)   = sqrt(grd%g11_sgy(j,i))
     grd%sq_g22_sgy(j,i)   = sqrt(grd%g22_sgy(j,i))
     grd%insq_g22_sgy(j,i) = 1.0_wp/grd%sq_g22_sgy(j,i)
  end do
  end do

!-------- Corresponding area of grid points --------

do i=0, grd%IMAX
do j=0, grd%JMAX
   grd%area(j,i) = grd%sq_g11_g(j,i)*grd%sq_g22_g(j,i)*grd%dxi*grd%deta
   grd%dt_darea(j,i) = tmr%dtime/grd%area(j,i)
end do
end do

!-------- Inner-point flag --------

grd%flag_inner_point = .true.

grd%flag_inner_point(0,:)    = .false.
grd%flag_inner_point(grd%JMAX,:) = .false.

grd%flag_inner_point(:,0)    = .false.
grd%flag_inner_point(:,grd%IMAX) = .false.

do i=0, grd%IMAX
do j=0, grd%JMAX

   if (grd%flag_inner_point(j,i)) then
      grd%sq_g22_x_1(j,i) = grd%sq_g22_sgx(j,i-1)
      grd%sq_g22_x_2(j,i) = grd%sq_g22_sgx(j,i)
      grd%sq_g11_y_1(j,i) = grd%sq_g11_sgy(j-1,i)
      grd%sq_g11_y_2(j,i) = grd%sq_g11_sgy(j,i)

    else
      grd%sq_g22_x_1(j,i) = 0.0_wp
      grd%sq_g22_x_2(j,i) = 0.0_wp
      grd%sq_g11_y_1(j,i) = 0.0_wp
      grd%sq_g11_y_2(j,i) = 0.0_wp

   end if

 enddo
 enddo

!-------- Term abbreviations

grd%at7 = 2.0_wp/RHO*tmr%dtime_temp

grd%aw1 = tmr%dtime_temp/grd%dzeta_t
grd%aw2 = tmr%dtime_temp/grd%dzeta_t
grd%aw3 = tmr%dtime_temp/grd%dzeta_t
grd%aw4 = tmr%dtime_temp/grd%dzeta_t
grd%aw5 = NUE/RHO*tmr%dtime_temp/(grd%dzeta_t**2)
grd%aw7 = 2.0_wp/(RHO*L)*tmr%dtime_temp
grd%aw8 = BETA**2/(RHO*L) &
      *(kappa_val(0.0_wp)-kappa_val(-1.0_wp))*tmr%dtime_temp
grd%aw9 = BETA/L*tmr%dtime_temp

grd%ai3 = par%agediff*tmr%dtime_temp/(grd%dzeta_t**2)

grd%atr1 = KAPPA_R/(RHO_C_R*grd%H_R**2)*tmr%dtime_temp/(grd%dzeta_r**2)

if (grd%flag_aa_nonzero) then
   grd%am1 = grd%aa*BETA*grd%dzeta_c/(grd%ea-1.0_wp)
   grd%am2 = grd%aa*L*RHO*grd%dzeta_c/(grd%ea-1.0_wp)
else
   grd%am1 = BETA*grd%dzeta_c
   grd%am2 = L*RHO*grd%dzeta_c
end if

if (grd%flag_aa_nonzero) then
   grd%acb1 = (grd%ea-1.0_wp)/grd%aa/grd%dzeta_c
else
   grd%acb1 = 1.0_wp/grd%dzeta_c
end if

grd%acb2 = KAPPA_R/grd%H_R/grd%dzeta_r
grd%acb3 = RHO*G
grd%acb4 = RHO*G

grd%alb1 = grd%H_R/KAPPA_R*grd%dzeta_r

grd%aqtld = grd%dzeta_t/tmr%dtime_temp

grd%dtt_2dxi  = 0.5_wp*tmr%dtime_temp/grd%dxi
grd%dtt_2deta = 0.5_wp*tmr%dtime_temp/grd%deta

do kc=0, grd%KCMAX

   if (grd%flag_aa_nonzero) then

      grd%at1(kc)   = (grd%ea-1.0_wp)/(grd%aa*grd%eaz_c(kc))*tmr%dtime_temp/grd%dzeta_c
      grd%at2_1(kc) = (grd%ea-1.0_wp)/(grd%aa*grd%eaz_c(kc))*tmr%dtime_temp/grd%dzeta_c
      grd%at2_2(kc) = (grd%eaz_c(kc)-1.0_wp)/(grd%aa*grd%eaz_c(kc)) &
                  *tmr%dtime_temp/grd%dzeta_c
      grd%at3_1(kc) = (grd%ea-1.0_wp)/(grd%aa*grd%eaz_c(kc))*tmr%dtime_temp/grd%dzeta_c
      grd%at3_2(kc) = (grd%eaz_c(kc)-1.0_wp)/(grd%aa*grd%eaz_c(kc)) &
                  *tmr%dtime_temp/grd%dzeta_c
      grd%at4_1(kc) = (grd%ea-1.0_wp)/(grd%aa*grd%eaz_c(kc))*tmr%dtime_temp/grd%dzeta_c
      grd%at4_2(kc) = (grd%eaz_c(kc)-1.0_wp)/(grd%aa*grd%eaz_c(kc)) &
                  *tmr%dtime_temp/grd%dzeta_c
      grd%at5(kc)   = (grd%ea-1.0_wp)/(RHO*grd%aa*grd%eaz_c(kc)) &
                  *tmr%dtime_temp/grd%dzeta_c
      if (kc /= grd%KCMAX) then
         grd%at6(kc) = (grd%ea-1.0_wp) &
                   /(grd%aa*exp(grd%aa*0.5_wp*(grd%zeta_c(kc)+grd%zeta_c(kc+1)))) &
                   /grd%dzeta_c
      else
         grd%at6(kc) = 0.0_wp
      end if
      grd%ai1(kc) = par%agediff*(grd%ea-1.0_wp)/(grd%aa*grd%eaz_c(kc)) &
                *tmr%dtime_temp/grd%dzeta_c
      if (kc /= grd%KCMAX) then
         grd%ai2(kc) = (grd%ea-1.0_wp) &
                   /(grd%aa*exp(grd%aa*0.5_wp*(grd%zeta_c(kc)+grd%zeta_c(kc+1)))) &
                   /grd%dzeta_c
      else
         grd%ai2(kc) = 0.0_wp
      end if
      grd%aqtlde(kc) = (grd%aa*grd%eaz_c(kc))/(grd%ea-1.0_wp)*grd%dzeta_c/tmr%dtime_temp
      grd%am3(kc)    = (grd%aa*grd%eaz_c(kc))/(grd%ea-1.0_wp)*grd%dzeta_c*BETA

   else

      grd%at1(kc)   = tmr%dtime_temp/grd%dzeta_c
      grd%at2_1(kc) = tmr%dtime_temp/grd%dzeta_c
      grd%at2_2(kc) = grd%zeta_c(kc) &
                  *tmr%dtime_temp/grd%dzeta_c
      grd%at3_1(kc) = tmr%dtime_temp/grd%dzeta_c
      grd%at3_2(kc) = grd%zeta_c(kc) &
                  *tmr%dtime_temp/grd%dzeta_c
      grd%at4_1(kc) = tmr%dtime_temp/grd%dzeta_c
      grd%at4_2(kc) = grd%zeta_c(kc) &
                  *tmr%dtime_temp/grd%dzeta_c
      grd%at5(kc)   = 1.0_wp/RHO &
                  *tmr%dtime_temp/grd%dzeta_c
      if (kc /= grd%KCMAX) then
         grd%at6(kc) = 1.0_wp/grd%dzeta_c
      else
         grd%at6(kc) = 0.0_wp
      end if
      grd%ai1(kc) = par%agediff &
                *tmr%dtime_temp/grd%dzeta_c
      if (kc /= grd%KCMAX) then
         grd%ai2(kc) = 1.0_wp &
                   /grd%dzeta_c
      else
         grd%ai2(kc) = 0.0_wp
      end if
      grd%aqtlde(kc) = grd%dzeta_c/tmr%dtime_temp
      grd%am3(kc)    = grd%dzeta_c*BETA

   end if

end do

!-------- Term abbreviations --------

do kc=0, grd%KCMAX
   if (grd%flag_aa_nonzero) then
      grd%avxy3(kc) = grd%aa*grd%eaz_c(kc)/(grd%ea-1.0_wp)*grd%dzeta_c
      grd%aqxy1(kc) = grd%aa/(grd%ea-1.0_wp)*grd%eaz_c(kc)*grd%dzeta_c
   else
      grd%avxy3(kc) = grd%dzeta_c
      grd%aqxy1(kc) = grd%dzeta_c
   end if
end do

!-------- Abbreviations --------

grd%azs2 = tmr%dtime/(grd%dxi*grd%dxi)
grd%azs3 = tmr%dtime/(grd%deta*grd%deta)

!-------- Abbreviations --------

if (grd%flag_aa_nonzero) then
   grd%aqbm1 = (grd%ea-1.0_wp)/(RHO*L*grd%aa*grd%dzeta_c)
else
   grd%aqbm1 = 1.0_wp/(RHO*L*grd%dzeta_c)
end if

grd%aqbm2  = KAPPA_R/(RHO*L*grd%H_R*grd%dzeta_r)
grd%aqbm3a = G/L
grd%aqbm3b = 1.0_wp/(RHO*L)
grd%aqbm4  = BETA/(RHO*L)

!-------- Term abbreviations --------

do kc=0, grd%KCMAX
   if (grd%flag_aa_nonzero) then
      grd%avz3(kc)  = grd%aa*grd%eaz_c(kc)/(grd%ea-1.0_wp)*grd%dzeta_c
   else
      grd%avz3(kc)  = grd%dzeta_c
   end if
end do

do kc=0, grd%KCMAX-1

   grd%zeta_c_sgz(kc) = (kc+0.5_wp)*grd%dzeta_c

   if (grd%flag_aa_nonzero) then
      grd%eaz_c_sgz(kc)          = exp(grd%aa*grd%zeta_c_sgz(kc))
      grd%eaz_c_quotient_sgz(kc) = (grd%eaz_c_sgz(kc)-1.0_wp)/(grd%ea-1.0_wp)
   else
      grd%eaz_c_sgz(kc)          = 1.0_wp
      grd%eaz_c_quotient_sgz(kc) = grd%zeta_c_sgz(kc)
   end if

end do

!-------- Term abbreviations --------

  grd%fact_x   = grd%dxi_inv *grd%insq_g11_g
  grd%fact_y   = grd%deta_inv*grd%insq_g22_g

  do kc=0, grd%KCMAX
     if (grd%flag_aa_nonzero) then
        grd%fact_z_c(kc)  = (grd%ea-1.0_wp)/(grd%aa*grd%eaz_c(kc)*grd%dzeta_c)
     else
        grd%fact_z_c(kc)  = 1.0_wp/grd%dzeta_c
     end if
  end do

  grd%fact_z_t  = 1.0_wp/grd%dzeta_t


!-------- Term abbreviations --------

  grd%atm1 = BETA*(1.0_wp-grd%eaz_c_quotient)
  grd%atm2 = BETA*(1.0_wp-grd%zeta_t)


  end subroutine sico_grid_init


  subroutine sico_grid_alloc(grd)

    type(sico_grid_class) :: grd

   allocate(grd%xi                 (0:grd%IMAX) ) 
   allocate(grd%eta                (0:grd%JMAX) )
   allocate(grd%zeta_c             (0:grd%KCMAX))
   allocate(grd%zeta_t             (0:grd%KTMAX))
   allocate(grd%zeta_r             (0:grd%KRMAX))
   allocate(grd%eaz_c              (0:grd%KCMAX))
   allocate(grd%eaz_c_quotient     (0:grd%KCMAX))
   allocate(grd%zeta_c_sgz         (0:grd%KCMAX))
   allocate(grd%eaz_c_sgz          (0:grd%KCMAX))
   allocate(grd%eaz_c_quotient_sgz (0:grd%KCMAX))

   allocate(grd%dist_dxdy(-grd%JMAX:grd%JMAX,-grd%IMAX:grd%IMAX))

   allocate(grd%lambda           (0:grd%JMAX,0:grd%IMAX)) 
   allocate(grd%phi              (0:grd%JMAX,0:grd%IMAX)) 
   allocate(grd%area             (0:grd%JMAX,0:grd%IMAX)) 
   allocate(grd%dt_darea          (0:grd%JMAX,0:grd%IMAX)) 

   allocate(grd%g11_g            (0:grd%JMAX,0:grd%IMAX)) 
   allocate(grd%g22_g            (0:grd%JMAX,0:grd%IMAX)) 
   allocate(grd%sq_g11_g         (0:grd%JMAX,0:grd%IMAX)) 
   allocate(grd%sq_g22_g         (0:grd%JMAX,0:grd%IMAX)) 
   allocate(grd%insq_g11_g       (0:grd%JMAX,0:grd%IMAX)) 
   allocate(grd%insq_g22_g       (0:grd%JMAX,0:grd%IMAX)) 
   allocate(grd%g11_sgx          (0:grd%JMAX,0:grd%IMAX)) 
   allocate(grd%g11_sgy          (0:grd%JMAX,0:grd%IMAX)) 
   allocate(grd%g22_sgx          (0:grd%JMAX,0:grd%IMAX)) 
   allocate(grd%g22_sgy          (0:grd%JMAX,0:grd%IMAX)) 
   allocate(grd%sq_g11_sgx       (0:grd%JMAX,0:grd%IMAX)) 
   allocate(grd%sq_g11_sgy       (0:grd%JMAX,0:grd%IMAX)) 
   allocate(grd%sq_g22_sgx       (0:grd%JMAX,0:grd%IMAX)) 
   allocate(grd%sq_g22_sgy       (0:grd%JMAX,0:grd%IMAX)) 
   allocate(grd%insq_g11_sgx     (0:grd%JMAX,0:grd%IMAX)) 
   allocate(grd%insq_g22_sgy     (0:grd%JMAX,0:grd%IMAX)) 
   allocate(grd%flag_inner_point (0:grd%JMAX,0:grd%IMAX)) 
   allocate(grd%sq_g22_x_1       (0:grd%JMAX,0:grd%IMAX))
   allocate(grd%sq_g22_x_2       (0:grd%JMAX,0:grd%IMAX))
   allocate(grd%sq_g11_y_1       (0:grd%JMAX,0:grd%IMAX))
   allocate(grd%sq_g11_y_2       (0:grd%JMAX,0:grd%IMAX))

   allocate(grd%ii((grd%IMAX+1)*(grd%JMAX+1)))
   allocate(grd%jj((grd%IMAX+1)*(grd%JMAX+1)))
   allocate(grd%nn(0:grd%JMAX,0:grd%IMAX))

   allocate(grd%n2i((grd%IMAX+1)*(grd%JMAX+1)))
   allocate(grd%n2j((grd%IMAX+1)*(grd%JMAX+1)))
   allocate(grd%ij2n(0:grd%JMAX,0:grd%IMAX))

   allocate(grd%fact_x(0:grd%JMAX,0:grd%IMAX))
   allocate(grd%fact_y(0:grd%JMAX,0:grd%IMAX))
   allocate(grd%fact_z_c(0:grd%KCMAX))

   allocate(grd%avxy3 (0:grd%KCMAX)) 
   allocate(grd%aqxy1 (0:grd%KCMAX)) 
   allocate(grd%avz3  (0:grd%KCMAX)) 

   allocate(grd%atm1  (0:grd%KCMAX)) 
   allocate(grd%atm2  (0:grd%KTMAX)) 

   allocate(grd%at1   (0:grd%KCMAX)) 
   allocate(grd%at2_1 (0:grd%KCMAX))
   allocate(grd%at2_2 (0:grd%KCMAX))
   allocate(grd%at3_1 (0:grd%KCMAX))
   allocate(grd%at3_2 (0:grd%KCMAX))
   allocate(grd%at4_1 (0:grd%KCMAX))
   allocate(grd%at4_2 (0:grd%KCMAX))
   allocate(grd%at5   (0:grd%KCMAX))
   allocate(grd%at6   (0:grd%KCMAX))
   allocate(grd%ai1   (0:grd%KCMAX))
   allocate(grd%ai2   (0:grd%KCMAX))
   allocate(grd%aqtlde(0:grd%KCMAX))
   allocate(grd%am3   (0:grd%KCMAX)) 

  end subroutine sico_grid_alloc

!-------------------------------------------------------------------------------
!! Components g11 and g22 of the metric tensor for the
!! stereographical projection.
!!------------------------------------------------------------------------------
  subroutine metric_stereo(x_val, y_val, K, g11_r, g22_r)

  implicit none

  real(wp), intent(in)  :: x_val, y_val
  real(wp), intent(in)  :: K
  real(wp), intent(out) :: g11_r, g22_r

  g11_r = 1.0_wp / ( K**2*(1.0_wp+(x_val**2+y_val**2)/(2.0_wp*R*K)**2)**2 )

  g22_r = g11_r

  end subroutine metric_stereo

!-------------------------------------------------------------------------------
! Components g11 and g22 of the metric tensor for geographical coordinates.
!------------------------------------------------------------------------------
  subroutine metric_geogr(phi_val, g11_r, g22_r)

  implicit none

  real(wp), intent(in)  :: phi_val
  real(wp), intent(out) :: g11_r, g22_r

  g11_r = R**2*(cos(phi_val))**2

  g22_r = R**2

  end subroutine metric_geogr


end module sico_grid_mod
