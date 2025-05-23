!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : d y n _ v e g _ m o d
!
!  Purpose : dynamic vegetation
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit and Andrey Ganopolski
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
module dyn_veg_mod

  use precision, only : wp
  use timer, only : doy, year, day_year, time_eoy_lnd 
  use constants, only : T0
  use control, only : check_carbon
  use lnd_grid, only : npft, ncarb, nl, nlc, flag_tree, flag_shrub, flag_grass
  use lnd_grid, only : ic_min, ic_peat, ic_shelf, ic_ice, ic_lake
  use lnd_params, only : dt_v, dt_day_v
  use lnd_params, only : l_dynveg, l_fixlai, veg_par, pft_par


   private
   public :: dyn_veg

contains


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  d y n _ v e g
  !   Purpose    :  implicitly solve vegetation carbon and vegetation fraction dynamics
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine dyn_veg(co2, f_veg ,f_veg_old, f_ice_grd, f_ice_grd_old, f_ice_nbr, f_lake, f_lake_old, f_shelf, f_shelf_old, &
                    f_crop, f_pasture, gamma_luc, &
                    z_veg_std, gamma_ice, &
                    gamma_dist, gamma_dist_cum, npp_ann, npp13_ann, npp14_ann, &
                    gamma_leaf, lambda, lai_bal, sai, &
                    root_frac, litter_in_frac, &
                    veg_c, veg_c13, veg_c14, leaf_c, stem_c, root_c, &
                    seed_frac, pft_frac, &
                    veg_c_above, veg_c13_above, veg_c14_above, veg_c_below, veg_c13_below, veg_c14_below, &
                    veg_h, litterfall, litterfall13, litterfall14, npp_real, npp13_real, npp14_real, &
                    carbon_bal_veg, carbon13_bal_veg, carbon14_bal_veg, i,j)

    implicit none

    real(wp), intent(in) :: co2
    real(wp), intent(in) :: f_veg, f_veg_old
    real(wp), intent(in) :: f_ice_grd, f_ice_grd_old, f_ice_nbr
    real(wp), intent(in) :: f_lake, f_lake_old
    real(wp), intent(in) :: f_shelf, f_shelf_old
    real(wp), intent(in) :: f_crop, f_pasture
    real(wp), intent(in) :: z_veg_std
    real(wp), dimension(:), intent(inout) :: gamma_luc, gamma_ice, gamma_dist, gamma_dist_cum
    real(wp), dimension(:), intent(inout) :: npp_ann, npp13_ann, npp14_ann
    real(wp), dimension(:), intent(inout) :: gamma_leaf, lambda, lai_bal, sai
    real(wp), dimension(:,:), intent(in) :: root_frac
    real(wp), dimension(:), intent(in) :: litter_in_frac
    real(wp), dimension(:), intent(inout) :: veg_c, veg_c13, veg_c14
    real(wp), dimension(:), intent(inout) :: leaf_c
    real(wp), dimension(:), intent(inout) :: stem_c
    real(wp), dimension(:), intent(inout) :: root_c
    real(wp), dimension(:), intent(inout) :: pft_frac, seed_frac
    real(wp), intent(inout) :: veg_c_above, veg_c13_above, veg_c14_above
    real(wp), dimension(:), intent(inout) :: veg_c_below, veg_c13_below, veg_c14_below
    real(wp), dimension(:), intent(out) :: veg_h
    real(wp), dimension(:,:), intent(out) :: litterfall, litterfall13, litterfall14
    real(wp), intent(out) :: npp_real, npp13_real, npp14_real
    real(wp), intent(out) :: carbon_bal_veg, carbon13_bal_veg, carbon14_bal_veg

    integer :: i, j
    integer :: n, m, k, n_iter
    real(wp) :: growth_sum
    real(wp) :: sum_v
    real(wp) :: d_f_veg, d_f_ice, d_f_shelf, d_f_lake, fac_to_buried
    real(wp) :: a, b, c, cluc, f_woody
    real(wp), dimension(npft) :: pft_frac_old, growth, sum_c, cc, sum_cij
    real(wp) :: lai_max
    real(wp) :: f_z_std
    real(wp), dimension(npft, npft) :: delta
    real(wp), dimension(npft) :: npp, npp13, npp14
    real(wp), dimension(npft) :: veg_c_old, veg_c13_old, veg_c14_old
    real(wp), dimension(npft) :: stem_c_old, leaf_c_old, root_c_old, lai_bal_old 
    real(wp), dimension(nl,npft) :: litter_pft, litter_dist_comp
    real(wp), dimension(nl,npft) :: litter13_pft, litter13_dist_comp
    real(wp), dimension(nl,npft) :: litter14_pft, litter14_dist_comp
    real(wp), dimension(nl) :: litter_shelf, litter13_shelf, litter14_shelf
    real(wp), dimension(nl) :: litter_lake, litter13_lake, litter14_lake
    real(wp), dimension(nl) :: litter_ice_to_veg, litter13_ice_to_veg, litter14_ice_to_veg
    real(wp), dimension(nl) :: litter_ice_to_buried, litter13_ice_to_buried, litter14_ice_to_buried
    real(wp), parameter :: log_2 = log(2._wp)

    real(wp) :: lai_n, lai_np1, f, fp


    f_z_std = exp(-z_veg_std/veg_par%z_veg_std_crit)

    litterfall = 0._wp
    litterfall13 = 0._wp
    litterfall14 = 0._wp

    ! when ice or water area retreating reduce pft fractions to increase bare soil
    if (f_veg.gt.f_veg_old .and. doy.eq.dt_day_v) then
      pft_frac = pft_frac * f_veg_old/f_veg 
      where (pft_frac.lt.veg_par%seed_fraction) pft_frac = veg_par%seed_fraction
    endif

    pft_frac_old = pft_frac

    if (f_veg.gt.0._wp) then

      ! annual npp 
      npp   = npp_ann
      npp13 = npp13_ann
      npp14 = npp14_ann

      ! update disturbance rate at end of year
      if (time_eoy_lnd) then
        gamma_dist = gamma_dist_cum 
        gamma_dist_cum = 0._wp
      endif

      ! save old values
      veg_c_old    = veg_c
      leaf_c_old   = leaf_c
      stem_c_old   = stem_c
      root_c_old   = root_c
      lai_bal_old  = lai_bal

      veg_c13_old  = veg_c13
      veg_c14_old  = veg_c14

      do n=1,npft

        if( pft_par%i_c4(n) .eq. 0 ) then
          ! C3 path
          ! adjust lai_max for CO2 concentration, to avoid changing lambda (partitioning between Cveg and pft expansion)
          lai_max = pft_par%lai_max(n) * (1._wp+veg_par%c_lambda_co2*log(co2/280._wp)/log_2*(1._wp-exp(-co2/40._wp)))
        else
          ! C4 path
          ! constant lai_max
          lai_max = pft_par%lai_max(n)
        endif
        if (lai_bal(n) .ge. lai_max) then
          lambda(n) = 1._wp
        else if (lai_bal(n) .lt. lai_max .and. lai_bal(n) .gt. (pft_par%lai_min(n)+0.1_wp)) then
          lambda(n) = ((lai_bal(n)-pft_par%lai_min(n))/(lai_max-pft_par%lai_min(n)))**veg_par%lambda_exp
        else
          lambda(n) = 0._wp
        endif

        lambda(n) = max(veg_par%lambda_min,lambda(n))
        lambda(n) = min(veg_par%lambda_max,lambda(n))

        ! solve iteratively for new balanced leaf area index using Newton-Rhapson method
        lai_n = lai_bal(n)
        n_iter = 0
        do 
          n_iter = n_iter + 1
          f = 2._wp*lai_n/(pft_par%sla(n)*dt_v) + pft_par%awl(n)/dt_v*lai_n**pft_par%bwl &
            + gamma_leaf(n)/pft_par%sla(n)*lai_n + pft_par%gamma_root(n)/pft_par%sla(n)*lai_n + pft_par%gamma_stem(n)*pft_par%awl(n)*lai_n**pft_par%bwl &
            - (1._wp-lambda(n))*npp(n) - 2._wp*lai_bal(n)/(pft_par%sla(n)*dt_v) - pft_par%awl(n)/dt_v*lai_bal(n)**pft_par%bwl
          fp = 2._wp/(pft_par%sla(n)*dt_v) + pft_par%awl(n)*pft_par%bwl/dt_v*lai_n**(pft_par%bwl-1._wp) &
            + gamma_leaf(n)/pft_par%sla(n) + pft_par%gamma_root(n)/pft_par%sla(n) + pft_par%gamma_stem(n)*pft_par%awl(n)*pft_par%bwl*lai_n**(pft_par%bwl-1._wp) 
          lai_np1 = lai_n - f/fp
          if (abs(lai_np1-lai_n).lt.veg_par%delta_lai_conv .or. n_iter.ge.10) exit
          lai_n = lai_np1
          if (lai_n.le.0._wp) exit
        enddo
        if (.not.l_fixlai) then
          lai_bal(n) = lai_n
        endif
        
        ! limit balanced LAI
        if (lai_bal(n).lt.pft_par%lai_min(n)) then
          lai_bal(n) = pft_par%lai_min(n)
          ! no NPP if balanced LAI equal to minimum value
          npp(n) = 0._wp
          npp13(n) = 0._wp
          npp14(n) = 0._wp
        endif

        if (lai_bal(n).eq.pft_par%lai_min(n) .and. lai_bal_old(n).eq.pft_par%lai_min(n)) then
          ! no NPP if balanced LAI equal to minimum value
          npp(n) = 0._wp
          npp13(n) = 0._wp
          npp14(n) = 0._wp
        endif

        ! update vegetation carbon pools based on new balanced LAI
        leaf_c(n) = lai_bal(n)/pft_par%sla(n)
        root_c(n) = leaf_c(n)
        stem_c(n) = pft_par%awl(n)*lai_bal(n)**pft_par%bwl
        veg_c(n)  = leaf_c(n)+root_c(n)+stem_c(n) 

        ! stem area index simply scaled to balanced leaf area index
        sai(n) = lai_bal(n) * veg_par%sai_scale

        ! update height of PFTs
        veg_h(n) = pft_par%awh(n) * lai_bal(n) 
        veg_h(n) = max(veg_par%veg_h_min, veg_h(n))

      enddo ! PFT loop


      ! ----------------- vegetation dynamics ---------------------------

      ! competition and colonisation rate
      do n=1,npft
        cc(n) = lambda(n)*npp(n)/veg_c(n)
      enddo

      ! competition/dominance factors
      do n=1,npft
        do m=1,npft
          if (n.eq.m) then
            delta(n,m) = 1._wp
          else if (flag_tree(n).eq.1 .and. flag_tree(m).eq.0) then
            ! tree dominance over shrubs and grass
            delta(n,m) = 0._wp
            delta(m,n) = 1._wp
          else if (flag_shrub(n).eq.1 .and. flag_grass(m).eq.1) then ! shrub dominance over grass
            ! shrubs dominance over grass
            delta(n,m) = 0._wp
            delta(m,n) = 1._wp
          else if ((flag_tree(n).eq.1 .and. flag_tree(m).eq.1) .or. (flag_grass(n).eq.1 .and. flag_grass(m).eq.1)) then ! tree-tree or grass-grass competition
            if (cc(n).gt.cc(m)) then
              delta(n,m) = 0._wp
              delta(m,n) = 1._wp
            else
              delta(n,m) = 1._wp
              delta(m,n) = 0._wp
            endif
          endif
        enddo
      enddo

      ! fraction of woody PFTs in grid cell
      f_woody = 0._wp
      do n=1,npft
        if (flag_tree(n).eq.1 .or. flag_shrub(n).eq.1) then
          f_woody = f_woody + pft_frac_old(n)
        endif
      enddo
      ! 'disturbance rate' from land use change if fraction of woody PFTs in grid cell is
      ! larger than fraction not affected by land use 
      do n=1,npft
        if (flag_tree(n).eq.1 .or. flag_shrub(n).eq.1) then
          gamma_luc(n) = veg_par%gamma_luc*max(0._wp, (f_woody-(1._wp-f_crop-f_pasture)))
        else
          gamma_luc(n) = 0._wp
        endif
      enddo

      ! additional 'disturbance rate' for forests from presence of ice sheets in grid cell
      do n=1,npft
        if (flag_tree(n).eq.1) then
          ! linear with ice sheet fraction, reduced over rough terrain to account for possibility that ice forms on mointain tops 
          ! while valleys can still be warm and support vegetation
          gamma_ice(n) = veg_par%gamma_ice*f_ice_nbr*f_z_std 
        else
          gamma_ice(n) = 0._wp
        endif
      enddo

      do n=1,npft

        sum_c(n) = 0._wp
        sum_cij(n) = 0._wp
        do k=1,npft
          if( k .ne. n) sum_c(n) = sum_c(n) + delta(n,k)*pft_frac_old(k)
          if( k .ne. n) sum_cij(n) = sum_cij(n) + delta(n,k)
        enddo

        if(sum(delta(n,:)).gt.5._wp) print *,'sum_c(n)ij > 5',n,sum(delta(n,:))
        if(sum_c(n).gt.1._wp) then
          print *,'sum_c(n)ij*v > 1',n,sum_c(n)
          print *,sum(pft_frac_old),sum(pft_frac_old).gt.1._wp
          print *,'pft_Frac_old',pft_frac_old
          print *,'seed_frac',seed_frac
          !stop
        endif

        ! land use state term, following Burton 2019, eq. 2
        if (flag_tree(n).eq.1 .or. flag_shrub(n).eq.1) then
          cluc = (f_crop + f_pasture)
        else
          cluc = 0._wp
        endif

        if (l_dynveg) then
          !update vegetation fraction

          ! semi-implicit scheme, implicit in nth-pft
          ! solve quadratic equation
          a = lambda(n)*npp(n)/veg_c_old(n) * delta(n,n)
          b = 1._wp/dt_v &
            + gamma_dist(n) + gamma_luc(n) + gamma_ice(n) &
            - lambda(n)*npp(n)/veg_c_old(n) * (1._wp - sum_c(n) - cluc)
          c = -pft_frac_old(n)/dt_v
          !if (a.eq.0._wp) then
          if (a.lt.epsilon(1._wp)) then
            pft_frac(n) = -c / b
          else
            pft_frac(n) = ( -b + sqrt(b**2-4._wp*a*c) ) / ( 2._wp*a )  ! positive solution
          endif

          ! if no seeds, block pft fraction at minimum fraction
          if (seed_frac(n).eq.0._wp .and. pft_frac_old(n).eq.veg_par%seed_fraction) then
            pft_frac(n) = veg_par%seed_fraction
          endif

        endif

      enddo  ! end npft

      ! options to apply deforestation
      if (veg_par%i_deforest.eq.0) then
      else if (veg_par%i_deforest.eq.1) then
        ! bare soil everywhere
        pft_frac(1:npft) = 0._wp
      else if (veg_par%i_deforest.eq.2) then
        ! convert forest to grassland
        if (pft_frac(3).gt.pft_frac(4)) then
          pft_frac(3) = pft_frac(3) + pft_frac(1) + pft_frac(2)
        else
          pft_frac(4) = pft_frac(4) + pft_frac(1) + pft_frac(2)
        endif
        pft_frac(1) = 0._wp
        pft_frac(2) = 0._wp
      else if (veg_par%i_deforest.eq.21) then
        ! convert tropical forest to grassland
        if (pft_frac(3).gt.pft_frac(4)) then
          pft_frac(3) = pft_frac(3) + pft_frac(1) 
        else
          pft_frac(4) = pft_frac(4) + pft_frac(1)
        endif
        pft_frac(1) = 0._wp
      else if (veg_par%i_deforest.eq.22) then
        ! convert boreal forest to grassland
        if (pft_frac(3).gt.pft_frac(4)) then
          pft_frac(3) = pft_frac(3) + pft_frac(2) 
        else
          pft_frac(4) = pft_frac(4) + pft_frac(2)
        endif
        pft_frac(2) = 0._wp
      else if (veg_par%i_deforest.eq.3) then
        ! convert forest to shrubs
        pft_frac(5) = pft_frac(5) + pft_frac(1) + pft_frac(2)
        pft_frac(1) = 0._wp
        pft_frac(2) = 0._wp
      else if (veg_par%i_deforest.eq.4) then
        ! convert grassland to bare soil
        pft_frac(3) = 0._wp
        pft_frac(4) = 0._wp
      endif

      ! ensure that seed fraction is maintained
      do n=1,npft
        if (pft_frac(n).lt.veg_par%seed_fraction) then
          pft_frac(n) = veg_par%seed_fraction
        endif
      enddo

      ! check that vegetation fraction not greater than 1 and renormalize if necessary
      norm: do
        sum_v = 0._wp
        do n=1,npft
          sum_v = sum_v + max(veg_par%seed_fraction, pft_frac(n))
        enddo
        !sum_v = sum(pft_frac)
        if( sum_v .gt. 1._wp ) then
          growth_sum = 0._wp
          do n=1,npft
            growth(n) = lambda(n)*npp(n)/veg_c(n)*sum_cij(n)*pft_frac(n) / gamma_dist(n)
            growth(n) = max(1.e-10_wp,growth(n))
          enddo
          do n=1,npft
            if( pft_frac(n) .lt. veg_par%seed_lim ) then
              if( pft_frac(n) .lt. veg_par%seed_fraction) pft_frac(n) = veg_par%seed_fraction
            else
              growth_sum = growth_sum + growth(n)
            endif
          enddo
          do n=1,npft
            if( pft_frac(n) .ge. veg_par%seed_lim ) pft_frac(n) = pft_frac(n) - growth(n)/growth_sum * (sum_v-1._wp)
          enddo
        endif
        if(minval(pft_frac).lt.veg_par%seed_fraction) then
          cycle norm
        endif
        exit norm
      enddo norm

      ! compute litterfall
      do n=1,npft

        if (lai_bal(n).eq.pft_par%lai_min(n)) then
          litter_pft(:,n) = -((leaf_c(n)-leaf_c_old(n) + stem_c(n)-stem_c_old(n))*litter_in_frac(:) &
            + (root_c(n)-root_c_old(n))*root_frac(:,n))/dt_v
          if (sum(litter_pft(:,n)).lt.1e-20_wp) then
            litter_pft(:,n)   = 0._wp
          endif
        else
          ! distribute litterfall to soil layers 
          litter_pft(:,n) = (gamma_leaf(n)*leaf_c(n) &  ! kgC/m2/s, leaf litter
            + pft_par%gamma_stem(n)*stem_c(n))*litter_in_frac(:) &  ! kgC/m2/s, stem litter
            + pft_par%gamma_root(n)*root_c(n)*root_frac(:,n)  ! kgC/m2/s, root litter
        endif

        ! litter from disturbance and competition
        litter_dist_comp(:,n)   = - (veg_c_old(n) * (pft_frac(n)-pft_frac_old(n))/dt_v - lambda(n)*npp(n)*pft_frac(n)) &
          * ((leaf_c_old(n)+stem_c_old(n))*litter_in_frac(:)+root_frac(:,n)*root_c_old(n))/veg_c_old(n)
        litter13_dist_comp(:,n) = veg_c13_old(n)/veg_c_old(n)*litter_dist_comp(:,n)
        litter14_dist_comp(:,n) = veg_c14_old(n)/veg_c_old(n)*litter_dist_comp(:,n)

        litter_dist_comp(:,n)   = max(0._wp,litter_dist_comp(:,n))
        litter13_dist_comp(:,n) = max(0._wp,litter13_dist_comp(:,n))
        litter14_dist_comp(:,n) = max(0._wp,litter14_dist_comp(:,n))

        ! update C13 and C14
        ! implicit scheme
        veg_c13(n) = (veg_c13(n) + (1._wp-lambda(n))*npp13(n)*dt_v) / (1._wp + 1._wp/veg_c(n)*sum(litter_pft(:,n))*dt_v)
        veg_c14(n) = (veg_c14(n) + (1._wp-lambda(n))*npp14(n)*dt_v) / (1._wp + 1._wp/veg_c(n)*sum(litter_pft(:,n))*dt_v)

        ! isotope litter
        litter13_pft(:,n)  = veg_c13(n)/veg_c(n)*litter_pft(:,n)
        litter14_pft(:,n)  = veg_c14(n)/veg_c(n)*litter_pft(:,n)

      enddo

      do k=1,nl
        do n=1,npft
          litterfall(k,ic_min)   = litterfall(k,ic_min)   + pft_frac(n)*litter_pft(k,n)   + litter_dist_comp(k,n)
          litterfall13(k,ic_min) = litterfall13(k,ic_min) + pft_frac(n)*litter13_pft(k,n) + litter13_dist_comp(k,n)
          litterfall14(k,ic_min) = litterfall14(k,ic_min) + pft_frac(n)*litter14_pft(k,n) + litter14_dist_comp(k,n)
        enddo
      enddo

      ! local carbon balance
      if (check_carbon) then

        do n=1,npft

          if(abs((veg_c(n)-veg_c_old(n)) - (1._wp-lambda(n))*npp(n)*dt_v + sum(litter_pft(:,n))*dt_v) .gt. 1.e-5_wp) then
            print *,''
            print *,'local carbon balance',n,(veg_c(n)-veg_c_old(n)) - (1._wp-lambda(n))*npp(n)*dt_v + sum(litter_pft(:,n))*dt_v
            print *,'lai',lai_bal(n)
            print *,'lai_old',lai_bal_old(n)
            print *,'lambda',lambda(n)
            print *,'npp',npp(n)
            print *,'vegc,vegc_old',veg_c(n),veg_c_old(n)
            print *,'litter',sum(litter_pft(:,n))
            print *,'dvegc,lambda*npp,litter',veg_c(n)-veg_c_old(n),(1._wp-lambda(n))*npp(n)*dt_v,sum(litter_pft(:,n))*dt_v
            stop
          endif

          if(abs((veg_c13(n)-veg_c13_old(n)) - (1._wp-lambda(n))*npp13(n)*dt_v + sum(litter13_pft(:,n))*dt_v) .gt. 1.e-5_wp) then
            print *,''
            print *,year,doy
            print *,n
            print *,'local carbon 13 balance',n,veg_c13(n)-veg_c13_old(n) - (1._wp-lambda(n))*npp13(n)*dt_v + sum(litter13_pft(:,n))*dt_v
            print *,'lai',lai_bal(n)
            print *,'lai_old',lai_bal_old(n)
            print *,'dvegc13,lambda*npp13,litter13',veg_c13(n)-veg_c13_old(n),(1._wp-lambda(n))*npp13(n)*dt_v,sum(litter13_pft(:,n))*dt_v
            print *,'lambda',lambda(n)
            print *,'npp',npp(n)
            print *,'npp13',npp13(n)
            print *,'veg_c13,veg_c13_old',veg_c13(n),veg_c13_old(n)
            print *,'veg_c,veg_c_old',veg_c(n),veg_c_old(n)
            print *,'leaf_c,leaf_c_old',leaf_c(n),leaf_c_old(n)
            print *,'root_c,root_c_old',root_c(n),root_c_old(n)
            print *,'stem_c,stem_c_old',stem_c(n),stem_c_old(n)
            !stop
          endif

          if(abs((veg_c14(n)-veg_c14_old(n)) - (1._wp-lambda(n))*npp14(n)*dt_v + sum(litter14_pft(:,n))*dt_v) .gt. 1.e-10_wp) then
            print *,''
            print *,year,doy
            print *,'local carbon 14 balance',n,veg_c14(n)-veg_c14_old(n) - (1._wp-lambda(n))*npp14(n)*dt_v + sum(litter14_pft(:,n))*dt_v
            print *,'lai',lai_bal(n)
            print *,'lai_old',lai_bal_old(n)
            print *,'dvegc14,lambda*npp14,litter14',veg_c14(n)-veg_c14_old(n),(1._wp-lambda(n))*npp14(n)*dt_v,sum(litter14_pft(:,n))*dt_v
            print *,'lambda',lambda(n)
            print *,'npp',npp(n)
            print *,'npp14',npp14(n)
            print *,'veg_c14,veg_c14_old',veg_c14(n),veg_c14_old(n)
            print *,'veg_c,veg_c_old',veg_c(n),veg_c_old(n)
            print *,'leaf_c,leaf_c_old',leaf_c(n),leaf_c_old(n)
            print *,'root_c,root_c_old',root_c(n),root_c_old(n)
            print *,'stem_c,stem_c_old',stem_c(n),stem_c_old(n)
            !stop
          endif

        enddo

      endif

      ! check for negative litterfall, should not happen
      if(minval(litterfall(:,ic_min)) .lt. 0._wp ) then
        !print *
        !print *,' WARNING litterfall < 0',sum(litterfall(:,ic_min))*dt_v,i,j
        !print *,'litter',litterfall(:,ic_min)*dt_v
        !print *,'litter_dist_comp',litter_dist_comp*dt_v
        !print *,'litter_pft',litter_pft*dt_v
        !print *,'lai',lai_bal
        !print *,'lai_old',lai_bal_old
        !print *,'dveg_c',veg_c-veg_c_old
        !print *,'veg_c',veg_c
        !print *,'veg_c_old',veg_c_old
        !print *,'leaf_c_old',leaf_c_old
        !print *,'stem_c_old',stem_c_old
        !print *,'root_c_old',root_c_old
        !print *,'pft_frac',pft_frac
        !print *,'pft_frac_old',pft_frac_old
        !print *,'npp',npp*dt_v
        !print *,'lambda',lambda
        !print *,'root_frac',root_frac
        litterfall(:,ic_min) = 0._wp
      endif
      if( minval(litterfall13(:,ic_min)) .lt. 0._wp ) then
        !print *
        !print *,' WARNING litterfall13 < 0',sum(litterfall13(:,ic_min))*dt_v,i,j
        !print *,'litter13_dist_comp',litter13_dist_comp
        !print *,'litter13_pft',litter13_pft
        !print *,'lai',lai_bal
        !print *,'lai_old',lai_bal_old
        !print *,'dveg_c13',veg_c13-veg_c13_old
        !print *,'veg_c13',veg_c13
        !print *,'veg_c13_old',veg_c13_old
        !print *,'veg_f',pft_frac
        !print *,'veg_f_old',pft_frac_old
        !print *,'npp13',npp13*dt_v
        litterfall13(:,ic_min) = 0._wp
      endif

      if( minval(litterfall14(:,ic_min)) .lt. 0._wp ) then
        !print *,' WARNING litterfall14 < 0',sum(litterfall14(:,ic_min))*dt_v,i,j
        litterfall14(:,ic_min) = 0._wp
      endif

      ! actual net carbon uptake trough photosynthesis, kgC/m2/s, used to compute net land carbon balance
      npp_real   = sum(npp*pft_frac)
      npp13_real = sum(npp13*pft_frac)
      npp14_real = sum(npp14*pft_frac)

      if (check_carbon) then

        carbon_bal_veg = sum(veg_c*pft_frac)-sum(veg_c_old*pft_frac_old) &
          - sum(npp*pft_frac)*dt_v &
          + sum(litterfall(:,ic_min))*dt_v
        if (abs(carbon_bal_veg).gt.1.e-5_wp) then 
          print *,' '
          print *,'vegetation carbon balance',carbon_bal_veg,i,j
          print *,'year,doy',year,doy
          print *,'lai',lai_bal
          print *,'lai_old',lai_bal_old
          print *,'lambda',lambda
          print *,'npp',npp*dt_v
          print *,'gamma_dist',gamma_dist
          print *,'litter_pft',pft_frac*sum(litter_pft,1)*dt_v
          print *,'litter_dist',sum(litter_dist_comp,1)*dt_v
          print *,'dveg_c',veg_c-veg_c_old
          print *,'veg_c',veg_c
          print *,'veg_c_old',veg_c_old
          print *,'dveg_f',pft_frac-pft_frac_old
          print *,'veg_f',pft_frac
          print *,'veg_f_old',pft_frac_old
          print *,'veg_c*pft_frac - veg_c_old*pft_frac_old',veg_c*pft_frac-veg_c_old*pft_frac_old
          print *,'npp*pft_frac',-npp*pft_frac*dt_v
          print *,'litterfall pfts',pft_frac*sum(litter_pft,1)*dt_v+sum(litter_dist_comp,1)*dt_v
          print *,'sum(dveg_c*dveg_f)',sum(veg_c*pft_frac)-sum(veg_c_old*pft_frac_old)
          print *,'sum(npp*pft_frac)',-sum(npp*pft_frac)*dt_v
          print *,'sum(litter)',sum(litterfall(:,ic_min))*dt_v
          print *,'cc',cc
          stop
        endif

        carbon13_bal_veg = sum(veg_c13*pft_frac)-sum(veg_c13_old*pft_frac_old) &
          - sum(npp13*pft_frac)*dt_v &
          + sum(litterfall13(:,ic_min))*dt_v
        if( abs(carbon13_bal_veg) .gt. 1.e-5_wp) then
          print *,' '
          print *,'vegetation carbon 13 balance',carbon13_bal_veg,i,j
          print *,'lai',lai_bal
          print *,'lai_old',lai_bal_old
          print *,'dveg_c',veg_c-veg_c_old
          print *,'dveg_c13',veg_c13-veg_c13_old
          print *,'veg_c13',veg_c13
          print *,'veg_c13_old',veg_c13_old
          print *,'veg_f',pft_frac
          print *,'veg_f_old',pft_frac_old
          print *,'npp13',npp13*dt_v
          print *,'litter13_pft',sum(litter13_pft,1)*dt_v
          print *,'litter13_dist',sum(litter13_dist_comp,1)*dt_v
          print *,'litter13',sum(litterfall13(:,ic_min))*dt_v
          print *,'dveg13_c*dveg_f',sum(veg_c13*pft_frac)-sum(veg_c13_old*pft_frac_old)
          print *,'lambda npp13',- sum((1._wp-lambda)*npp13*pft_frac_old + lambda*npp13*pft_frac)*dt_v
          print *,'lambda npp13',lambda*npp13*pft_frac*dt_v
          !stop
        endif

        carbon14_bal_veg = sum(veg_c14*pft_frac)-sum(veg_c14_old*pft_frac_old) &
          - sum(npp14*pft_frac)*dt_v &
          + sum(litterfall14(:,ic_min))*dt_v
        if( abs(carbon14_bal_veg) .gt. 1.e-10_wp) then !print *,'vegetation carbon 14 balance',carbon14_bal_veg(i,j),i,j
          print *,' '
          print *,'vegetation carbon 14 balance',carbon14_bal_veg,i,j
          print *,'lai',lai_bal
          print *,'lai_old',lai_bal_old
          print *,'dveg_c',veg_c-veg_c_old
          print *,'dveg_c14',veg_c14-veg_c14_old
          print *,'veg_c14',veg_c14
          print *,'veg_c14_old',veg_c14_old
          print *,'veg_f',pft_frac
          print *,'veg_f_old',pft_frac_old
          print *,'npp14',npp14*dt_v
          print *,'litter14_pft',sum(litter14_pft,1)*dt_v
          print *,'litter14_dist',sum(litter14_dist_comp,1)*dt_v
          print *,'litter14',sum(litterfall14(:,ic_min))*dt_v
          print *,'dveg_c*dveg_f',sum(veg_c14*pft_frac)-sum(veg_c14_old*pft_frac_old)
          print *,'lambda npp',- sum((1._wp-lambda)*npp14*pft_frac_old + lambda*npp14*pft_frac)*dt_v
          !stop
        endif

      endif


    else ! only ice and/or shelf and/or lake

      npp_real = 0._wp
      npp13_real = 0._wp
      npp14_real = 0._wp

      lai_bal = pft_par%lai_min

      do n=1,npft
        ! update vegetation carbon based on new lai
        veg_c(n) = 2._wp * lai_bal(n)/pft_par%sla(n) + pft_par%awl(n)*lai_bal(n)**pft_par%bwl
        ! update carbon pools
        leaf_c(n) = lai_bal(n)/pft_par%sla(n)
        root_c(n) = leaf_c(n)
        stem_c(n) = pft_par%awl(n)*lai_bal(n)**pft_par%bwl
        ! update height of PFTs
        veg_h(n) = stem_c(n) / (pft_par%aws(n)*0.01_wp) * (pft_par%awl(n)/stem_c(n))**(1._wp/pft_par%bwl)
      enddo

      pft_frac = seed_frac

      carbon_bal_veg = 0._wp
      carbon13_bal_veg = 0._wp
      carbon14_bal_veg = 0._wp

    endif


    ! add litter from ice/lake/shelf expansion, at first carbon model call of the year only
    ! generate litter from present vegetation carbon (in kg/m2) 
    ! assume no prior knowledge about the distribution of pfts within the grid cell
    if( doy .eq. dt_day_v ) then
      ! save above- and below ground vegetation carbon to use for litter production for ice, lake and shelf
      if( f_veg .gt. 0._wp ) then
        ! above ground 
        veg_c_above = sum(leaf_c*pft_frac) + sum(stem_c*pft_frac) ! kgC/m2
        veg_c13_above = veg_c_above*sum(veg_c13*pft_frac)/sum(veg_c*pft_frac) ! kgC/m2
        veg_c14_above = veg_c_above*sum(veg_c14*pft_frac)/sum(veg_c*pft_frac) ! kgC/m2
        ! below ground, in each soil layer
        do n=1,nl
          veg_c_below(n) = sum(root_c*root_frac(n,:)*pft_frac) ! kgC/m2
          veg_c13_below(n) = veg_c_below(n)*sum(veg_c13*pft_frac)/sum(veg_c*pft_frac) ! kgC/m2
          veg_c14_below(n) = veg_c_below(n)*sum(veg_c14*pft_frac)/sum(veg_c*pft_frac) ! kgC/m2
        enddo
      endif

      d_f_veg = f_veg - f_veg_old
      d_f_ice = f_ice_grd - f_ice_grd_old
      d_f_shelf = f_shelf - f_shelf_old
      d_f_lake = f_lake - f_lake_old

      litter_ice_to_buried = 0._wp
      litter_ice_to_veg = 0._wp
      litter_shelf = 0._wp
      litter_lake = 0._wp
      litter13_ice_to_buried = 0._wp
      litter13_ice_to_veg = 0._wp
      litter13_shelf = 0._wp
      litter13_lake = 0._wp
      litter14_ice_to_buried = 0._wp
      litter14_ice_to_veg = 0._wp
      litter14_shelf = 0._wp
      litter14_lake = 0._wp

      ! give priority to ice expansion relative to lake and shelf expansion
      if( d_f_ice .gt. 0._wp .and. d_f_veg .lt. 0._wp  ) then
        if( f_veg .gt. veg_par%f_veg_crit ) then
          fac_to_buried = veg_par%f_lit_to_ice
        else
          fac_to_buried = 1._wp - (1._wp - veg_par%f_lit_to_ice)/veg_par%f_veg_crit * f_veg
        endif
        ! divide aboveground litter (leaf + stem) from ice expansion into buried and vegetated grid cell part
        litter_ice_to_buried(1) = min(d_f_ice,-d_f_veg)/f_ice_grd * fac_to_buried        * veg_c_above  ! kgC/m2ice
        litter_ice_to_buried(1) = litter_ice_to_buried(1) + min(d_f_ice,-d_f_veg)/f_ice_grd * veg_c_below(1) ! kgC/m2ice
        litter13_ice_to_buried(1) = min(d_f_ice,-d_f_veg)/f_ice_grd * fac_to_buried        * veg_c13_above  ! kgC/m2ice
        litter13_ice_to_buried(1) = litter13_ice_to_buried(1) + min(d_f_ice,-d_f_veg)/f_ice_grd * veg_c13_below(1) ! kgC/m2ice
        litter14_ice_to_buried(1) = min(d_f_ice,-d_f_veg)/f_ice_grd * fac_to_buried        * veg_c14_above  ! kgC/m2ice
        litter14_ice_to_buried(1) = litter14_ice_to_buried(1) + min(d_f_ice,-d_f_veg)/f_ice_grd * veg_c14_below(1) ! kgC/m2ice
        ! litter from roots is all buried below the ice sheet
        do n=2,nl 
          litter_ice_to_buried(n) = min(d_f_ice,-d_f_veg)/f_ice_grd * veg_c_below(n) ! kgC/m2ice
          litter13_ice_to_buried(n) = min(d_f_ice,-d_f_veg)/f_ice_grd * veg_c13_below(n) ! kgC/m2ice
          litter14_ice_to_buried(n) = min(d_f_ice,-d_f_veg)/f_ice_grd * veg_c14_below(n) ! kgC/m2ice
        enddo
        if( f_veg .gt. 0._wp ) then
          litter_ice_to_veg(1) = min(d_f_ice,-d_f_veg)/f_veg * (1._wp-fac_to_buried) * veg_c_above ! kgC/m2veg
          litter13_ice_to_veg(1) = min(d_f_ice,-d_f_veg)/f_veg * (1._wp-fac_to_buried) * veg_c13_above ! kgC/m2veg
          litter14_ice_to_veg(1) = min(d_f_ice,-d_f_veg)/f_veg * (1._wp-fac_to_buried) * veg_c14_above ! kgC/m2veg
        endif
        ! recompute remaining vegetation and ice fraction change after accounting for ice
        d_f_veg = d_f_veg + min(d_f_ice,-d_f_veg)
        d_f_ice = d_f_ice - min(d_f_ice,-d_f_veg)
      endif
      ! next check for lake expansion into (remaining) vegetation fraction
      if( d_f_lake .gt. 0._wp .and. d_f_veg .lt. 0._wp  ) then
        ! put all litter from lake expansion into lake
        litter_lake(1) = min(d_f_lake,-d_f_veg)/f_lake * veg_c_above
        litter_lake(1) = litter_lake(1) + min(d_f_lake,-d_f_veg)/f_lake * veg_c_below(1)
        litter13_lake(1) = min(d_f_lake,-d_f_veg)/f_lake * veg_c13_above
        litter13_lake(1) = litter13_lake(1) + min(d_f_lake,-d_f_veg)/f_lake * veg_c13_below(1)
        litter14_lake(1) = min(d_f_lake,-d_f_veg)/f_lake * veg_c14_above
        litter14_lake(1) = litter14_lake(1) + min(d_f_lake,-d_f_veg)/f_lake * veg_c14_below(1)
        do n=2,nl
          litter_lake(n) = min(d_f_lake,-d_f_veg)/f_lake * veg_c_below(n)
          litter13_lake(n) = min(d_f_lake,-d_f_veg)/f_lake * veg_c13_below(n)
          litter14_lake(n) = min(d_f_lake,-d_f_veg)/f_lake * veg_c14_below(n)
        enddo
        ! recompute remaining vegetation and lake fraction change after accounting for lake
        d_f_veg = d_f_veg + min(d_f_lake,-d_f_veg)
        d_f_lake = d_f_lake - min(d_f_lake,-d_f_veg)
      endif
      ! next check for shelf expansion into (remaining) vegetation fraction
      if( d_f_shelf .gt. 0._wp .and. d_f_veg .lt. 0._wp  ) then
        ! put all litter from shelf expansion into shelf
        litter_shelf(1) = min(d_f_shelf,-d_f_veg)/f_shelf * veg_c_above
        litter_shelf(1) = litter_shelf(1) + min(d_f_shelf,-d_f_veg)/f_shelf * veg_c_below(1)
        litter13_shelf(1) = min(d_f_shelf,-d_f_veg)/f_shelf * veg_c13_above
        litter13_shelf(1) = litter13_shelf(1) + min(d_f_shelf,-d_f_veg)/f_shelf * veg_c13_below(1)
        litter14_shelf(1) = min(d_f_shelf,-d_f_veg)/f_shelf * veg_c14_above
        litter14_shelf(1) = litter14_shelf(1) + min(d_f_shelf,-d_f_veg)/f_shelf * veg_c14_below(1)
        do n=2,nl
          litter_shelf(n) = min(d_f_shelf,-d_f_veg)/f_shelf * veg_c_below(n)
          litter13_shelf(n) = min(d_f_shelf,-d_f_veg)/f_shelf * veg_c13_below(n)
          litter14_shelf(n) = min(d_f_shelf,-d_f_veg)/f_shelf * veg_c14_below(n)
        enddo
      endif

      d_f_veg = f_veg - f_veg_old
      d_f_ice = f_ice_grd - f_ice_grd_old
      d_f_shelf = f_shelf - f_shelf_old
      d_f_lake = f_lake - f_lake_old

      litterfall(1:nl,ic_min) = litterfall(1:nl,ic_min) + litter_ice_to_veg/dt_v  ! kg/m2/s
      litterfall(1:nl,ic_ice)  = litter_ice_to_buried/dt_v ! kg/m2/s
      litterfall(1:nl,ic_shelf) = litter_shelf/dt_v ! kg/m2/s
      litterfall(1:nl,ic_lake) = litter_lake/dt_v ! kg/m2/s

      litterfall13(1:nl,ic_min) = litterfall13(1:nl,ic_min) + litter13_ice_to_veg/dt_v  ! kg/m2/s
      litterfall13(1:nl,ic_ice)  = litter13_ice_to_buried/dt_v ! kg/m2/s
      litterfall13(1:nl,ic_shelf) = litter13_shelf/dt_v ! kg/m2/s
      litterfall13(1:nl,ic_lake) = litter13_lake/dt_v ! kg/m2/s

      litterfall14(1:nl,ic_min) = litterfall14(1:nl,ic_min) + litter14_ice_to_veg/dt_v  ! kg/m2/s
      litterfall14(1:nl,ic_ice)  = litter14_ice_to_buried/dt_v ! kg/m2/s
      litterfall14(1:nl,ic_shelf) = litter14_shelf/dt_v ! kg/m2/s
      litterfall14(1:nl,ic_lake) = litter14_lake/dt_v ! kg/m2/s


    endif

    litterfall(1:nl,ic_peat) = litterfall(1:nl,ic_min)
    litterfall13(1:nl,ic_peat) = litterfall13(1:nl,ic_min)
    litterfall14(1:nl,ic_peat) = litterfall14(1:nl,ic_min)

    return

  end subroutine dyn_veg

end module dyn_veg_mod
