module modheatpump
  use mpi
  implicit none
  private

  public :: init_heatpump, heatpump, exit_heatpump
  public :: nhppoints_local, thl_dot_hp, w_hp_exhaust
#if defined(_GPU)
  public :: idhppts_local_d
#endif
  save

  integer, allocatable :: idhppts_global(:,:)  ! Global heat pump points array (indices with respect to cell centered coordinates)
  logical, allocatable :: lhpptsrank(:)        ! Logical array indicating if heat pump points are on this rank
  
  integer              :: nhppoints_local = 0  ! count of this rank's heat pump points
#if defined(_GPU)
  integer, allocatable, pinned :: idhppts_local(:,:)
  integer, device, allocatable :: idhppts_local_d(:,:)
#else
  integer, allocatable         :: idhppts_local(:,:) ! rank-local (i,j,k) indices
#endif

  real :: thl_dot_hp   ! Rate of change of temperature due to heat loss from the heat pump in [Km^3/s]
  real :: w_hp_exhaust ! Exhaust velocity at heat pump points in [m/s]

contains

  subroutine init_heatpump
    use modglobal, only : lheatpump, nhppoints, Q_dot_hp, QH_dot_hp, rhoa, cp, ifinput, cexpnr, ltempeq, ib, jb, dxi, dyi
    use modmpi,    only : myid, comm3d, mpierr
    use decomp_2d, only : zstart, zend
    implicit none

    integer :: n, m
    character(80) :: chmess

    if (.not.(lheatpump) .or. .not.(ltempeq) .or. (nhppoints<1)) return

    allocate(idhppts_global(nhppoints,3))  ! Allocate global heat pump points array
    allocate(lhpptsrank(nhppoints))        ! Allocate logical array for heat pump points on this rank

    ! read global heat pump points
    if(myid==0) then
      open (ifinput,file='heatpump.inp.'//cexpnr)
        read (ifinput,'(a80)') chmess
        read (ifinput,'(a80)') chmess
        do n = 1, nhppoints
          read (ifinput,*) idhppts_global(n,1), idhppts_global(n,2), idhppts_global(n,3)
        end do
      close (ifinput)
    end if
    ! Broadcast the heat pump points to all processes
    call MPI_BCAST(idhppts_global, nhppoints*3, MPI_INTEGER, 0, comm3d, mpierr)

    ! Determine whether points are on this rank
    do n = 1, nhppoints
      if ((idhppts_global(n,1) >= zstart(1) .and. idhppts_global(n,1) <= zend(1)) .and. &
          (idhppts_global(n,2) >= zstart(2) .and. idhppts_global(n,2) <= zend(2))) then
        lhpptsrank(n) = .true.
      else
        lhpptsrank(n) = .false.
      end if
    end do
    
    nhppoints_local = count(lhpptsrank)
    if (nhppoints_local > 0) then
      allocate(idhppts_local(nhppoints_local, 3))
      m = 0
      do n = 1, nhppoints
        if (lhpptsrank(n)) then
          m = m + 1
          idhppts_local(m, 1) = idhppts_global(n, 1) - zstart(1) + ib
          idhppts_local(m, 2) = idhppts_global(n, 2) - zstart(2) + jb
          idhppts_local(m, 3) = idhppts_global(n, 3)
        end if
      end do
#if defined(_GPU)
      allocate(idhppts_local_d(nhppoints_local, 3))
      idhppts_local_d = idhppts_local
#endif
    end if

    thl_dot_hp = QH_dot_hp / (nhppoints*rhoa*cp) ! Calculate temperature change rate from heat loss [Km^3/s]
    
    w_hp_exhaust = (Q_dot_hp/nhppoints)*dxi*dyi  ! Calculate exhaust velocity at heat pump points [m/s]

  end subroutine init_heatpump

  subroutine heatpump
    use modglobal,  only : lheatpump, lfan_hp, nhppoints, ltempeq, dxi, dyi
#if defined(_GPU)
    use modcuda,    only : wm_d, w0_d, wp_d, thlp_d, dzfi_d
#else
    use modglobal,  only : dzfi
    use modfields,  only : wm, w0, wp, thlp
#endif
    implicit none

    integer :: n, i, j, k

    if (.not.(lheatpump) .or. .not.(ltempeq) .or. (nhppoints<1)) return
    if (nhppoints_local == 0) return

#if defined(_GPU)
    !$acc parallel loop default(present) private(i, j, k)
    do n = 1, nhppoints_local
      i = idhppts_local_d(n, 1)
      j = idhppts_local_d(n, 2)
      k = idhppts_local_d(n, 3)

      if (lfan_hp) then
        wm_d(i, j, k+1) = w_hp_exhaust
        w0_d(i, j, k+1) = w_hp_exhaust
        wp_d(i, j, k+1) = 0.
      end if

      thlp_d(i, j, k) = thlp_d(i, j, k) - thl_dot_hp * dxi * dyi * dzfi_d(k)
    end do
    !$acc end parallel loop
#else
    if (lfan_hp) then ! Heat pump fan is on
      do n = 1, nhppoints_local
        i = idhppts_local(n,1)
        j = idhppts_local(n,2)
        k = idhppts_local(n,3)
        wm(i,j,k+1) = w_hp_exhaust ! Set exhaust velocity at heat pump point [m/s], at input 'w' cell face k+1
        w0(i,j,k+1) = w_hp_exhaust
        wp(i,j,k+1) = 0.
      end do
    end if
    do n = 1, nhppoints_local
      i = idhppts_local(n,1)
      j = idhppts_local(n,2)
      k = idhppts_local(n,3)
      thlp(i,j,k) = thlp(i,j,k) - thl_dot_hp * dxi * dyi * dzfi(k) ! [K/s], at cell center k
    end do
#endif
  end subroutine heatpump

  subroutine exit_heatpump
    use modglobal, only : lheatpump, nhppoints, ltempeq
    implicit none

    if (.not.(lheatpump) .or. .not.(ltempeq) .or. (nhppoints<1)) return

    deallocate(idhppts_global, lhpptsrank)
    if (allocated(idhppts_local)) deallocate(idhppts_local)
#if defined(_GPU)
    if (allocated(idhppts_local_d)) deallocate(idhppts_local_d)
#endif
  end subroutine exit_heatpump

end module modheatpump
