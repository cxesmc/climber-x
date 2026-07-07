MODULE bgc_params

  USE precision, ONLY : wp

  IMPLICIT NONE

  PUBLIC

  logical :: l_sediments
  logical :: l_spinup_bgc 
  logical :: l_spinup_sed
  integer :: i_compensate
  logical :: l_conserve_phos
  logical :: l_conserve_sil 
  logical :: l_conserve_alk
  integer :: i_bgc_fw
  logical :: l_sed_hypsometry = .false.   !! sediment depth-class hypsometry (mirrors real bgc; off in dummy)
  real(wp) :: pk490 = 0.1_wp              !! light extinction coefficient (used by coupler coral hypsometry)
  real(wp), parameter :: rcar  = 122._wp

  ! advected tracers
  INTEGER, PARAMETER :: i_base_adv = 21,              &
       &                isco212    = 1,               &
       &                ialkali    = 2,               &
       &                iphosph    = 3,               &
       &                ioxygen    = 4,               &
       &                igasnit    = 5,               &
       &                iano3      = 6,               &
       &                isilica    = 7,               &
       &                idoc       = 8,               &
       &                iphy       = 9,               &
       &                izoo       = 10,              &
       &                ian2o      = 11,              &
       &                idms       = 12,              &
       &                iiron      = 13,              &
       &                ifdust     = 14,              &
       &                idet       = 15,              &         
       &                idetf      = 16,              &         
       &                idets      = 17,              &         
       &                icya       = 18,              &               
       &                ih2s       = 19,              &
       &                ihi        = 20,              &
       &                ico3       = 21

  INTEGER, PARAMETER :: isco213   = i_base_adv+1,     &
       &                isco214   = i_base_adv+2,     &
       &                idet13    = i_base_adv+3,     &
       &                idetf13   = i_base_adv+4,     &
       &                idets13   = i_base_adv+5,     &
       &                idet14    = i_base_adv+6,     &
       &                idetf14   = i_base_adv+7,     &
       &                idets14   = i_base_adv+8,     &
       &                iphy13    = i_base_adv+9,     &
       &                iphy14    = i_base_adv+10,     &       
       &                izoo13    = i_base_adv+11,     &
       &                izoo14    = i_base_adv+12,     &       
       &                idoc13    = i_base_adv+13,     &
       &                idoc14    = i_base_adv+14,     &       
       &                icya13    = i_base_adv+15,     &
       &                icya14    = i_base_adv+16,     &       
       &                i_iso_adv = 16

  ! total number of advected tracers
  INTEGER, PARAMETER :: ntraad=i_base_adv+i_iso_adv

  ! non-advected (fast sinking) tracers
  INTEGER, PARAMETER ::  icalc    = ntraad+1,         &
       &                 iopal    = ntraad+2,         &
       &                 i_base   = 2

  INTEGER, PARAMETER ::                               &
       &                 icalc13  = ntraad+i_base+1,  &
       &                 icalc14  = ntraad+i_base+2,  &
       &                 i_iso    = 2

  INTEGER, PARAMETER :: nocetra = ntraad + i_base + i_iso

  ! atmosphere
  INTEGER, PARAMETER :: iatmco2    = 1,               &
       &                iatmo2     = 2,               &
       &                iatmn2     = 3,               &
       &                i_base_atm = 3
  INTEGER, PARAMETER ::                               &
       &                iatmc13    = i_base_atm+1,    &
       &                iatmc14    = i_base_atm+2,    &
       &                i_iso_atm  = 2

  INTEGER, PARAMETER :: natm = i_base_atm + i_iso_atm

  ! sediment
  INTEGER, PARAMETER :: issso12   = 1,                &
       &                isssc12   = 2,                &
       &                issssil   = 3,                &
       &                issster   = 4,                &
       &                nsss_base = 4
  INTEGER, PARAMETER ::                               &
       &                issso13  = nsss_base + 1,     &
       &                issso14  = nsss_base + 2,     &
       &                isssc13  = nsss_base + 3,     &
       &                isssc14  = nsss_base + 4,     &
       &                nsss_iso = 4
  INTEGER, PARAMETER :: nsedtra = nsss_base + nsss_iso

  ! pore water tracers
  INTEGER, PARAMETER :: ipowaic    = 1,               &
       &                ipowaal    = 2,               &
       &                ipowaph    = 3,               &
       &                ipowaox    = 4,               &
       &                ipown2     = 5,               &
       &                ipowno3    = 6,               &
       &                ipowasi    = 7,               &
       &                ipowafe    = 8,               &
       &                ipowah2s   = 9,               &
       &                npowa_base = 9
  INTEGER, PARAMETER ::                               &
       &                ipowc13    = npowa_base + 1,  &
       &                ipowc14    = npowa_base + 2,  &
       &                npowa_iso  = 2
  INTEGER, PARAMETER :: npowtra=npowa_base+npowa_iso

END MODULE bgc_params
