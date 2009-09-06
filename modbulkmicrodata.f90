!> \file modbulkmicrodata.f90
!!  Variables necessary for the bulk microphysics

!>
!!  Variables necessary for the bulk microphysics
!>
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!

  module modbulkmicrodata

  use modglobal, only : ih,i1,jh,j1,k1, rhow
  implicit none
  save
  logical :: l_sb        = .true. , &!< SB scheme (.true.) / KK00 scheme (.false.)   (in namelist NAMMICROPHYSICS)
             l_sedc      = .true. , & !<  cloud droplet sedimentation flag             (in namelist NAMMICROPHYSICS)
             l_rain      = .true. , & !<  rain formation / evolution flag              (in namelist NAMMICROPHYSICS)
             l_mur_cst   = .false.    ! false = no constant value of mur (mur=f(Dv)) (in namelist NAMMICROPHYSICS)
  real    :: mur_cst     = 5        & !<  mur value if l_mur_cst=T                     (in namelist NAMMICROPHYSICS)
                 ,Nc_0 = 70e6       & !<  initial cloud droplet number
                 ,sig_g = 1.34      & !<  geom. std dev of cloud droplet DSD
                 ,sig_gr = 1.5        !<  geometric std dev of rain drop DSD

  logical :: l_lognormal = .false.    !<  log param of rain terminal velocities for rain sedim

  integer :: inr = 1, iqr=2

  real, parameter ::  D0_kk = 50e-6     & !<  diameter sep. cloud and prec. in KK00 scheme
                     ,qcmin = 1.e-7     & !<  Cloud mixing ratio treshold for calculations
                     ,qrmin = 1.e-13    & !<  Rain  mixing ratio treshold for calculations
!                     ,nuc = 0           & !< width parameter of cloud DSD
                     ,pluseps = 1e-25   &
                     ,mineps = -1e-50   &
                     ,epscloud= 0.01e-3 &
                     ,epsprec = 3.65e-5 & !<  RICO threshold
                     ,epsqr = 1.e-8     &
!  values picked by Verica Savic-Jovcic to optimize for Sc, note x and D have to be chosen consistently
!                ,xcmin = 4.2e-15       & !<  min mean mass of cw
!                ,xcmax = 6.5e-11       & !<  max mean mass of cw
!                ,xrmin = xcmax         & !<  min mean mass of pw
                     ,xrmax = 5.0e-6       & !<  max mean mass of pw
                     ,xrmaxkk = 5.2e-7     & !<  max mean mass of pw in KK00 scheme
!                ,Dvcmin = 2.0e-6       & !<  min mean diam. of cw
!                     ,Dvcmax = 49.8e-6      & !<  max mean diam. of cw
!                     ,Dvrmin = Dvcmax       & !<  min mean diam. of pw
!                     ,Dvrmax = 1000.0e-6    & !<  max mean diam. of pw
!  values given by SB2001
              ,xcmin = 4.2e-15     & !< \param xcmin  min mean mass of cw (D = 2.0e-6m)
              ,xcmax = 2.6e-10     & !<  max mean mass of cw (D = 80e-6m)
              ,xrmin = xcmax       & !<  min mean mass of pw
!               ,xrmax = 6.0e-07      & !<  max mean mass of pw
               ,Dvcmin = 2.0e-6     & !<  min mean diam. of cw
               ,Dvcmax = 79.2e-6    & !<  max mean diam. of cw
               ,Dvrmin = Dvcmax     & !<  min mean diam. of pw
               ,Dvrmax = 3000.0e-6  & !<  max mean diam. of pw
! NB in table1 in SB2006 komen weer andere getallen voor
! NB x_s is 'scheidingsdrop massa' en die mag dus best groter zijn dan bovengrens
! xcmax omdat die voor mean droplet mass staat!<  -> in gedachten houden
! NB de microphysica is heel gevoelig voor de waarde van de scheidingsdrop massa!<
        ,x_s = xcmax   & !<  drop mass sep. cloud and prec. part of DSD
        ,D_s = Dvcmax  & !<  diameter sep. cloud and prec. part of DSD
!          ,x_s = 2.6e-10  &
!          ,D_s = 79.2e-6  &
!          ,k_c = 9.44e9   & !<  Long Kernel coef. SB2001 [m^3 kg^-2 s^-1]
!          ,k_1 = 6.0e2    & !<  k_1 + k_2: coef. for phi function
!          ,k_2 = 0.68     & !<  in autoconversion rate SB2001
         ,k_c = 10.58e9 & !<  Long Kernel coef. SB2006 (k'cc)
!          ,k_c = 4.44e9   &  !<  Long Kernel coef. SB2006 (k'cc) for test
        ,k_1 = 4.0e2   & !<  k_1 + k_2: coef. for phi function
        ,k_2 = 0.70    & !<  in autoconversion rate SB2006
!
!          ,k_r = 5.78     & !<  Kernel coef. SB2001 [m^3 kg^-1 s^-1]
!          ,k_l = 5.e-4    & !<  coef for phi function in accr. rate
!          ,kappa_r = 0,   &
!          ,k_rr= k_r      &
         ,k_r = 5.25     & !<  Kernel SB2006
         ,k_l = 5.e-5    & !<  coef. for phi function in accr. rate
         ,kappa_r = 60.7 & !<  see eq. 11 SB2006
         ,k_rr = 7.12    & !<  idem dito

         ,Kt    = 2.5e-2  & !<  conductivity of heat [J/(sKm)]
         ,Dv    = 2.4e-5   & !<  diffusivity of water vapor [m2/s]
!  NB (see table 7.1 in Rogers: given Kt is for ~15 C while Dv is for > 30 C  2.4e-5
!  is value for ~ 15C How sensitive is G for this Aug 2006, ~5% -> Dv changed to 15 C value?
        ,c_St  = 1.19e8  & !<  Stokes fall vel. coef. [m^-1 s^-1]
!          ,pirhow = (pi*rhow)/6. & !< used in conversion of mass to diameter
                     ,pirhow = 3.14159*rhow/6.        &
         ,Rv = 461.5       & !<  specific gas constant for water vapor
         ,avf = 0.78       & !<  constants in vent. factor fv   (fv = 1. --> av=1,
         ,bvf = 0.308      & !<                                             bv=0 )
         ,nu_a = 1.41e-5   & !<  kin. viscosity of air [m2s-1]
         ,c_Nevap = 0.7    & !<  coeff for evap
         ,c_evapkk = 0.87  & !<  coeff for evap in KK00 scheme
         ,Sc_num = 0.71    &    !<  Schmidt number
         ,a_tvsb = 9.65    & !<  coeff in terminal velocity param
         ,b_tvsb = 9.8     & !<  coeff in terminal velocity param
         ,c_tvsb = 600.      !<  coeff in terminal velocity param


  real,allocatable, dimension(:,:,:) :: qc  & !<  cloud droplets mixing ratio [kg_w/kg_a]
                                       ,Nc  & !<  cloud droplets number conc.  [#/m^3]
                                       ,nuc & !<  width parameter of cloud DSD
                                       ,rhoz  !< slab averaged density in 3 dimensions

  real,allocatable, dimension(:,:,:) :: qr_spl, Nr_spl
                             !< prec. liq. water and conc. for sedim. time splitting
  real,allocatable, dimension(:,:,:) :: sedc,   & !<  sedimentation cloud droplets mix. ratio
                                        sed_qr, & !<  sedimentation rain drops mix. ratio
                                        sed_Nr    !<  sedimentation rain drop number conc.
  real ::  rho_c             &      !<  term to correct for density dep. of fall vel.
    ,k_au                     !<  coeff. for autoconversion rate
  real,allocatable, dimension(:,:,:) ::  &
    exnz               &      !<  3D exner function
    ,presz             &      !<  3D pressure
      ,Dvc               &      !<  cloud water mean diameter
          ,xc                &      !<  mean mass of cloud water droplets
    ,Dvr               &      !<  prec water mean diameter
    ,xr                &      !<  mean mass of prec. water drops
          ,mur               &      !<  mu parameter in rain gamma distribution
          ,lbdr              &      !<  slope parameter (lambda) in rain gamma distribution
    ,au                &      !<  autoconversion rate
    ,phi               &      !<  correction function (see SB2001)
    ,tau               &      !<  internal time scale
    ,ac                &      !<  accretion rate
    ,sc                &      !<  self collection rate
    ,br                &      !<  break-up rate
    ,evap              &      !<  mass tendency due to rain evap/cond
    ,Nevap             &      !<  concentration tendency due to rain evap/cond
    ,wfall_qr          &      !<  fall velocity for qr
    ,wfall_Nr                 !<  fall velocity for Nr
  real :: csed                      !<  parameter in cloud water grav. settling formula


  real, parameter ::  D_eq = 1.1E-3,  & !<  Parameters for break-up
            k_br = 1000       !<

   real,allocatable,dimension(:,:,:) :: Nr,Nrp,qltot,qr,qrp,thlpmcr,qtpmcr
   real,allocatable,dimension(:,:,:) :: precep

  real :: delt

  logical ,allocatable,dimension(:,:,:):: qcmask,qrmask
  end module modbulkmicrodata

