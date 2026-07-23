!> \file modbasestate.f90
!! Hydrostatic base state derived from the initial profiles and ps (issue #302).
!! Fixed per start; recomputed identically on cold, warm and strat starts.
module modbasestate
   implicit none
   save
   public :: initbasestate, exitbasestate, ps
   real              :: ps = 101325.  !<  Surface pressure [Pa]
   real, allocatable :: thl_b(:)  !< base-state liquid water potential temperature [K]
   real, allocatable :: qt_b(:)   !< base-state total specific humidity [kg/kg]
   real, allocatable :: thv_b(:)  !< base-state virtual potential temperature [K]
   real, allocatable :: pf_b(:)   !< base-state hydrostatic pressure, full levels [Pa]
   real, allocatable :: ph_b(:)   !< base-state hydrostatic pressure, half levels [Pa]
   real, allocatable :: exnf_b(:) !< base-state Exner function, full levels
   real, allocatable :: exnh_b(:) !< base-state Exner function, half levels

contains

   subroutine initbasestate(thlprof, qtprof)
      use modglobal,   only : kb, ke, kh, zf, dzf, dzh, grav, cp, rd, rv, pref0
      use modmpi,      only : myid
      real, intent(in) :: thlprof(kb:ke), qtprof(kb:ke)
      real :: rdocp, thvh_b
      integer :: k

      rdocp = rd/cp

      if (.not. allocated(thl_b)) then
         allocate (thl_b(kb:ke+kh), qt_b(kb:ke+kh), thv_b(kb:ke+kh))
         allocate (pf_b(kb:ke+kh), ph_b(kb:ke+kh))
         allocate (exnf_b(kb:ke+kh), exnh_b(kb:ke+kh))
      end if

      thl_b(kb:ke) = thlprof
      qt_b(kb:ke)  = qtprof
      thl_b(ke+kh) = thlprof(ke)
      qt_b(ke+kh)  = qtprof(ke)
      thv_b = thl_b*(1.+(rv/rd - 1.)*qt_b) ! ql = 0 in the base state

      ! hydrostatic integration from ps at the domain bottom (z = 0), same
      ! discrete scheme as fromztop; defined over the full column, including
      ! levels inside terrain (reference-column continuation, backlog section 1.5)
      ph_b(kb) = ps
      pf_b(kb) = (ps**rdocp - grav*(pref0**rdocp)*zf(kb)/(cp*thv_b(kb)))**(1./rdocp)
      do k = kb + 1, ke + kh
         thvh_b  = (thv_b(k)*dzf(k-1) + thv_b(k-1)*dzf(k))/(2.*dzh(k))
         pf_b(k) = (pf_b(k-1)**rdocp - grav*(pref0**rdocp)*dzh(k)/(cp*thvh_b))**(1./rdocp)
         ph_b(k) = (ph_b(k-1)**rdocp - grav*(pref0**rdocp)*dzf(k-1)/(cp*thv_b(k-1)))**(1./rdocp)
      end do

      exnf_b = (pf_b/pref0)**rdocp
      exnh_b = (ph_b/pref0)**rdocp

      if (myid == 0) then
         write (*, '(A,F8.3,A,F10.1,A)') ' Base state: thv_b(kb) = ', thv_b(kb), &
            ' K, ps = ', ps, ' Pa (derived from prof.inp, #302)'
      end if
   end subroutine initbasestate

   subroutine exitbasestate
      if (allocated(thl_b)) deallocate (thl_b, qt_b, thv_b, pf_b, ph_b, exnf_b, exnh_b)
   end subroutine exitbasestate
end module modbasestate
