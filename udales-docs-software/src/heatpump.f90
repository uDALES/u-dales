module modheatpump
  use mpi
  implicit none
  private

  public :: init_heatpump, heatpump, exit_heatpump
  save

  integer, allocatable :: idhppts_global(:,:) ! Global heat pump points array (indices with respect to cell centered coordinates)
  logical, allocatable :: lhpptsrank(:)  ! Logical array indicating if heat pump points are on this rank

  real :: thl_dot_hp  ! Rate of change of temperature due to heat loss from the heat pump in [Km^3/s]
  real :: w_hp_exhaust ! Exhaust velocity at heat pump points in [m/s]

contains

  subroutine init_heatpump
    use modglobal, only : lheatpump, nhppoints, Q_dot_hp, QH_dot_hp, rhoa, cp, ifinput, cexpnr, ltempeq, dxi, dyi
    use modmpi,    only : myid, comm3d, mpierr
    use decomp_2d, only : zstart, zend
    implicit none

    integer :: n
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
    
    thl_dot_hp = QH_dot_hp / (nhppoints*rhoa*cp) ! Calculate temperature change rate from heat loss [Km^3/s]
    
    w_hp_exhaust = (Q_dot_hp/nhppoints)*dxi*dyi  ! Calculate exhaust velocity at heat pump points [m/s]

  end subroutine init_heatpump

  subroutine heatpump
    use modglobal,  only : lheatpump, lfan_hp, nhppoints, dxi, dyi, dzfi, ltempeq
    use modfields,  only : wm, w0, wp, thlp
    use modmpi,     only : myidx, myidy
    use decomp_2d,  only : zsize
    implicit none

    integer :: n, i, j, k

    if (.not.(lheatpump) .or. .not.(ltempeq) .or. (nhppoints<1)) return
    
    do n = 1, nhppoints
      if (lhpptsrank(n)) then
        i = idhppts_global(n,1) - myidx*zsize(1)
        j = idhppts_global(n,2) - myidy*zsize(2)
        k = idhppts_global(n,3)

        if (lfan_hp) then ! Heat pump fan is on
          wm(i,j,k+1) = w_hp_exhaust ! Set exhaust velocity at heat pump point [m/s], at input 'w' cell face k+1
          w0(i,j,k+1) = w_hp_exhaust
          wp(i,j,k+1) = 0.
        end if
        
        thlp(i,j,k) = thlp(i,j,k) - thl_dot_hp * dxi * dyi * dzfi(k)  ! [K/s], at cell center k
      end if
    end do
  end subroutine heatpump

  subroutine exit_heatpump
    use modglobal, only : lheatpump, nhppoints, ltempeq
    implicit none

    if (.not.(lheatpump) .or. .not.(ltempeq) .or. (nhppoints<1)) return

    deallocate(idhppts_global,lhpptsrank) ! Deallocate global heat pump points array and logical array
  end subroutine exit_heatpump

end module modheatpump
