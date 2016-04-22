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
    integer, allocatable :: block(:,:)
    integer, allocatable :: xwallsshear(:,:)
    integer, allocatable :: ywallsp(:,:)
    integer, allocatable :: ywallsm(:,:)
    integer, allocatable :: zwallsshear(:,:)
    integer, allocatable :: xwallsnorm(:,:)
    integer, allocatable :: ywallsnorm(:,:)
    integer, allocatable :: zwallsnorm(:,:)
    integer, allocatable :: bcthltypex(:)            ! for 1:nxwallsshear set type of bc (0: flux, 1: fixed value)
    integer, allocatable :: bcthltypeyp(:)
    integer, allocatable :: bcthltypeym(:)
    integer, allocatable :: bcthltypez(:)
    integer, allocatable :: bcqttypex(:)            ! for 1:nxwallsshear set type of bc (0: flux, 1: fixed value)
    integer, allocatable :: bcqttypeyp(:)
    integer, allocatable :: bcqttypeym(:)
    integer, allocatable :: bcqttypez(:)
    real, allocatable    :: ibmxforce(:,:)           ! spanwise- and time-averaged force by ibm method.
    real, allocatable    :: ibmxforcevol(:,:)        ! spanwise- and time-averaged force by ibm method (complete volume)
    real, allocatable    :: ibmxforcevolp(:,:)        ! spanwise- and time-averaged force by ibm method (complete volume minus dp/dx)
    real, allocatable    :: bcthlvaluex(:)            ! for 1:nxwallsshear set value of bc (flux or value)
    real, allocatable    :: bcthlvalueyp(:)
    real, allocatable    :: bcthlvalueym(:)
    real, allocatable    :: bcthlvaluez(:)
    real, allocatable    :: bcqtvaluex(:)            ! for 1:nxwallsshear set value of bc (flux or value)
    real, allocatable    :: bcqtvalueyp(:)
    real, allocatable    :: bcqtvalueym(:)
    real, allocatable    :: bcqtvaluez(:)

    integer              :: nxwallsnorm
    integer              :: nywallsnorm
    integer              :: nzwallsnorm
    integer              :: nxwallsshear
    integer              :: nywallsp
    integer              :: nywallsm
    integer              :: nzwallsshear
        
end module
