!!
!!  \author Jasper Tomas,TU Delft, 31 March 2014
!!  \par Revision list
!!  \todo Documentation

module modibmdata
  implicit none
  save

    integer, allocatable :: xwallsglobal(:,:)
    integer, allocatable :: ywallsglobal(:,:)
    integer, allocatable :: zwallsglobal(:,:)
!    integer, allocatable :: block(:,:)
    integer, allocatable :: xwallsshear(:,:)
    integer, allocatable :: ywallsp(:,:)
    integer, allocatable :: ywallsm(:,:)
    integer, allocatable :: zwallsshear(:,:)
    integer, allocatable :: xwallsnorm(:,:)
    integer, allocatable :: ywallsnorm(:,:)
    integer, allocatable :: zwallsnorm(:,:)
    integer              :: nxwall                   !number of local xwalls
    integer, allocatable :: ixwall(:)                !index of block that is on local processor, used to determine local xwalls, ils13, 16.02.2017
    integer              :: nywall                   !number of local xwalls
    integer, allocatable :: iyminwall(:,:)                !index of block that is on local processor, used to determine local xwalls, ils13, 16.02.2017
    integer              :: nyminwall                   !number of local xwalls                                                                                  
    integer, allocatable :: iywall(:)                !index of block that is on local processor, used to determine local xwalls, ils13, 16.02.2017
    integer              :: nypluswall                   !number of local xwalls                                                                                                    
    integer, allocatable :: iypluswall(:,:)                !index of block that is on local processor, used to determine local xwalls, ils13, 16.02.2017

    real, allocatable    :: ibmxforce(:,:)           ! spanwise- and time-averaged force by ibm method.
    real, allocatable    :: ibmxforcevol(:,:)        ! spanwise- and time-averaged force by ibm method (complete volume)
    real, allocatable    :: ibmxforcevolp(:,:)        ! spanwise- and time-averaged force by ibm method (complete volume minus dp/dx)

    !
    real                 :: sumctm=0.
    real                 :: bcTfluxA=0. 
    real                 :: bcqfluxA=0.
    !fluxes for temperature and humidity at immersed boundaries
    real                 :: bctfxm=0.
    real                 :: bctfxp=0.
    real                 :: bctfym=0.
    real                 :: bctfyp=0.
    real                 :: bctfz=0.
    real                 :: bcqfxm=0.
    real                 :: bcqfxp=0.
    real                 :: bcqfym=0.
    real                 :: bcqfyp=0.
    real                 :: bcqfz=0.


    integer              :: nxwallsnorm
    integer              :: nywallsnorm
    integer              :: nzwallsnorm
    integer              :: nxwallsshear
    integer              :: nywallsp
    integer              :: nywallsm
    integer              :: nzwallsshear
       
     integer :: offset=1  !why do we need offset in modglobal? just use same value here to get some indeces right, ils13 20/03/2017                                                                               
 
end module
