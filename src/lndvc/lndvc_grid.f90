module lndvc_grid

    use precision, only : wp

    implicit none
    
    ! levels
    integer,  parameter :: nl = 5
    real(wp) :: z(0:nl) = (/0._wp,0.1_wp,0.3_wp,0.7_wp,1.5_wp,3.1_wp/)! z(0) is overwritten with 0.5*h_snow
    real(wp) :: z_int(0:nl)
    real(wp) :: dz(0:nl)
    real(wp) :: rdz(0:nl)
    real(wp) :: rdz_pos(0:nl)
    real(wp) :: rdz_neg(0:nl)
    ! levels for lakes, bottom layer adjusted for lake depth
    !integer,  parameter :: nl_l = 7
    !real(wp), dimension(0:nl_l) :: z_l = (/0._wp,0.2_wp,0.5_wp,1._wp,2._wp,10._wp,20._wp,50._wp/)! z_l(0) is overwritten with 0.5*h_snow
    integer,  parameter :: nl_l = 10
    real(wp) :: z_l(0:nl_l) = (/0._wp,0.2_wp,0.5_wp,1._wp,2._wp,5._wp,10._wp,20._wp,30._wp,50._wp,100._wp/)! z_l(0) is overwritten with 0.5*h_snow
    !real(wp) :: z_l(0:nl_l) = (/0._wp,0.2_wp,0.5_wp,1._wp,2._wp,5._wp,8._wp,11._wp,14._wp,17._wp,20._wp/)! z_l(0) is overwritten with 0.5*h_snow
    real(wp) :: z_int_l(0:nl_l)
    real(wp) :: dz_l(0:nl_l)
    real(wp) :: rdz_l(0:nl_l)
    real(wp) :: rdz_pos_l(0:nl_l)
    real(wp) :: rdz_neg_l(0:nl_l)

    ! levels for carbon (including burial layer) 
    integer,  parameter :: nlc = nl+1
    real(wp) :: z_c(0:nlc) 
    real(wp) :: z_int_c(0:nlc)
    real(wp) :: dz_c(0:nlc)
    real(wp) :: rdz_c(0:nlc)
    real(wp) :: rdz_pos_c(0:nlc)
    real(wp) :: rdz_neg_c(0:nlc)

    ! surface types
    integer,  parameter :: nsurf = 8 ! 5 pfts + bare soil + lake + ice
    integer,  parameter :: flag_pft(nsurf) = (/1,1,1,1,1,0,0,0/)
    integer,  parameter :: flag_veg(nsurf) = (/1,1,1,1,1,1,0,0/)
    integer,  parameter :: i_surf(nsurf) = (/1,2,3,4,5,6,7,8/)
    integer,  parameter :: i_bare = 6
    integer,  parameter :: i_lake = 7
    integer,  parameter :: i_ice = 8
    ! PFTs
    integer,  parameter :: npft = 5
    integer,  parameter :: ntrees = 2
    integer,  parameter :: ngrass = 2
    integer,  parameter :: nshrub = 1
    integer,  parameter :: i_pft(npft) = (/1,2,3,4,5/)
    integer,  parameter :: i_trees(ntrees) = (/1,2/)
    integer,  parameter :: i_grass(ngrass) = (/3,4/)
    integer,  parameter :: i_shrub(nshrub) = (/5/)
    integer,  parameter :: flag_tree(npft)  = (/1,1,0,0,0/)
    integer,  parameter :: flag_grass(npft) = (/0,0,1,1,0/)
    integer,  parameter :: flag_shrub(npft) = (/0,0,0,0,1/)
    integer,  parameter :: nveg = npft+1
    ! carbon
    integer,  parameter :: ncarb = 5
    integer,  parameter :: ic_min = 1
    integer,  parameter :: ic_peat = 2
    integer,  parameter :: ic_shelf = 3
    integer,  parameter :: ic_ice = 4
    integer,  parameter :: ic_lake = 5
    ! ground/snow types
    integer,  parameter :: nsoil = 3  ! vegetation, ice, lake
    integer,  parameter :: i_soil(nsoil) = (/1,2,3/)
    integer,  parameter :: is_veg = 1
    integer,  parameter :: is_ice = 2
    integer,  parameter :: is_lake = 3

contains

    subroutine lndvc_grid_init

        implicit none

        integer :: k

        z(0) = 0._wp 
        dz(0) = 0._wp
        ! vertical layers thickness
        dz(1) = 0.5_wp * ( z(1) + z(2) )
        do k=2,nl-1
            dz(k) = 0.5_wp * ( z(k+1) - z(k-1) )
        end do
        dz(nl) = z(nl) - z(nl-1)

        ! depth of vertical layer interfaces
        z_int(0) = 0._wp ! snow - soil interface
        do k=1,nl-1
            z_int(k) = 0.5_wp * ( z(k) + z(k+1) )
        end do
        z_int(nl) = z(nl) + 0.5_wp * dz(nl)

        write(*,*) 'z',z
        write(*,*) 'z_int',z_int
        write(*,*) 'dz',dz

        ! reciprocals to speed up fortran
        rdz(1:nl) = 1._wp/dz(1:nl)
        do k=1,nl-1
            rdz_pos(k) = 1._wp/(z(k+1)-z(k))
        end do
        rdz_pos(nl) = 0._wp
        do k=2,nl
            rdz_neg(k) = 1._wp/(z(k)-z(k-1))
        end do

        ! for lakes
        z_l(0) = 0._wp 
        dz_l(0) = 0._wp
        ! vertical layers thickness
        dz_l(1) = 0.5_wp * ( z_l(1) + z_l(2) )
        do k=2,nl_l-1
            dz_l(k) = 0.5_wp * ( z_l(k+1) - z_l(k-1) )
        end do
        dz_l(nl_l) = z_l(nl_l) - z_l(nl_l-1)

        ! depth of vertical layer interfaces
        z_int_l(0) = 0._wp ! snow - lake interface
        do k=1,nl_l-1
            z_int_l(k) = 0.5_wp * ( z_l(k) + z_l(k+1) )
        end do
        z_int_l(nl_l) = z_l(nl_l) + 0.5_wp * dz_l(nl_l)

        ! reciprocals to speed up fortran
        rdz_l(1:nl_l) = 1._wp/dz_l(1:nl_l)
        do k=1,nl_l-1
            rdz_pos_l(k) = 1._wp/(z_l(k+1)-z_l(k))
        end do
        rdz_pos_l(nl_l) = 0._wp
        do k=2,nl_l
            rdz_neg_l(k) = 1._wp/(z_l(k)-z_l(k-1))
        end do

        ! for carbon
        z_c(0:nl) = z
        z_c(nlc) = z(nl) + dz(nl)! 4.7_wp
        dz_c(0) = 0._wp
        ! vertical layers thickness
        dz_c(1) = 0.5_wp * ( z_c(1) + z_c(2) )
        do k=2,nlc-1
            dz_c(k) = 0.5_wp * ( z_c(k+1) - z_c(k-1) )
        end do
        dz_c(nlc) = z_c(nlc) - z_c(nlc-1)

        ! depth of vertical layer interfaces
        z_int_c(0) = 0._wp ! snow - soil interface
        do k=1,nlc-1
            z_int_c(k) = 0.5_wp * ( z_c(k) + z_c(k+1) )
        end do
        z_int_c(nlc) = z_c(nlc) + 0.5_wp * dz_c(nlc)

        write(*,*) 'z_c',z_c
        write(*,*) 'z_int_c',z_int_c
        write(*,*) 'dz_c',dz_c

        ! reciprocals to speed up fortran
        rdz_c(1:nlc) = 1._wp/dz_c(1:nlc)
        do k=1,nlc-1
            rdz_pos_c(k) = 1._wp/(z_c(k+1)-z_c(k))
        end do
        rdz_pos_c(nlc) = 0._wp
        do k=2,nlc
            rdz_neg_c(k) = 1._wp/(z_c(k)-z_c(k-1))
        end do

        return

    end subroutine lndvc_grid_init

end module lndvc_grid



