!> Parametrised longwave radiation, ported from DALES dev 075ac4e
!> src/modradiation.f90 (radiation, radpar longwave branch). Case-specific
!> terms (DYCOMS free troposphere) are case hooks, see modcases.
module modradiation
  use mpi   ! explicit MPI interfaces: gfortran rejects MPI_BCAST calls of different types in one file otherwise
  use modraddata, only : irad_none, irad_par, iradiation, timerad, rad_longw, rad_shortw, rad_ls, &
                         rad_absorber, rad_absorber_sv, rka, dlwtop, dlwbot, thlprad, lwd, lwu
  implicit none
  private
  public :: initradiation, radiation, exitradiation

  real :: tnextrad = 0.   !< time of the next radiation update [s]

contains

  subroutine initradiation
    use modglobal, only : ifnamopt, fname_options, ib, ie, jb, je, kb, ke
    use modmpi,    only : myid, comm3d, mpierr, my_real, mpi_logical, mpi_integer, mpi_character
    implicit none
    integer :: ierr
    character(160) :: msg
    namelist/NAMRADIATION/ iradiation, timerad, rad_longw, rad_shortw, rad_ls, &
         rad_absorber, rad_absorber_sv, rka, dlwtop, dlwbot

    if (myid == 0) then
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read (ifnamopt, NAMRADIATION, iostat=ierr)
      if (ierr > 0) then
        write(msg, '(a,i0)') 'Problem in namoptions NAMRADIATION, iostat = ', ierr
      end if
      close(ifnamopt)
    end if
    call MPI_BCAST(ierr, 1, MPI_INTEGER, 0, comm3d, mpierr)   ! all ranks abort together
    if (ierr > 0) call fail(msg)
    call MPI_BCAST(iradiation,  1, MPI_INTEGER, 0, comm3d, mpierr)
    call MPI_BCAST(timerad,     1, MY_REAL,     0, comm3d, mpierr)
    call MPI_BCAST(rad_longw,   1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(rad_shortw,  1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(rad_ls,      1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(rad_absorber, len(rad_absorber), MPI_CHARACTER, 0, comm3d, mpierr)
    call MPI_BCAST(rad_absorber_sv, 1, MPI_INTEGER, 0, comm3d, mpierr)
    call MPI_BCAST(rka,         1, MY_REAL,     0, comm3d, mpierr)
    call MPI_BCAST(dlwtop,      1, MY_REAL,     0, comm3d, mpierr)
    call MPI_BCAST(dlwbot,      1, MY_REAL,     0, comm3d, mpierr)

    rad_absorber = adjustl(rad_absorber)

    if (iradiation == irad_none) then
      if (myid == 0 .and. (rka /= 0. .or. dlwtop /= 0. .or. dlwbot /= 0. .or. rad_absorber_sv /= 0 &
                           .or. trim(rad_absorber) /= 'ql')) &
        write(0, '(a)') 'WARNING: &NAMRADIATION settings ignored: iradiation = 0'
      return
    end if
    ! allocated only when the scheme is on: with iradiation = 0 nothing reads these arrays
    ! (radiation and bldiag are guarded by iradiation; exitradiation tests allocated)
    allocate(thlprad(ib:ie,jb:je,kb:ke)); thlprad = 0.
    allocate(lwd(ib:ie,jb:je,kb:ke+1));   lwd = 0.
    allocate(lwu(ib:ie,jb:je,kb:ke+1));   lwu = 0.
    call check_radiation_options
  end subroutine initradiation

  !> Configuration checks for iradiation /= 0 (all ranks; rank 0 reports).
  subroutine check_radiation_options
    use modglobal, only : nsv, lmoist, ltempeq
    use modmpi,    only : myid
    implicit none
    character(160) :: msg

    if (iradiation /= irad_par) then
      write(msg, '(a,i0)') 'iradiation must be 0 (none) or 2 (parametrised); got ', iradiation
      call fail(msg)
    end if
    if (rad_shortw) call fail('rad_shortw is not supported in uDALES')
    if (.not. ltempeq) call fail('iradiation = 2 requires ltempeq = .true.')
    if (timerad < 0.) then
      write(msg, '(a,es12.4)') '&NAMRADIATION timerad must be >= 0; got ', timerad
      call fail(msg)
    end if
    if (dlwtop < 0. .or. dlwbot < 0.) then
      write(msg, '(a,2es12.4)') '&NAMRADIATION dlwtop and dlwbot must be >= 0; got ', dlwtop, dlwbot
      call fail(msg)
    end if
    if (rka <= 0. .or. (dlwtop <= 0. .and. dlwbot <= 0.)) then
      write(msg, '(a,3es12.4)') '&NAMRADIATION: iradiation = 2 requires rka > 0 and dlwtop > 0 or dlwbot > 0'// &
                                ' (no defaults); got rka, dlwtop, dlwbot =', rka, dlwtop, dlwbot
      call fail(msg)
    end if
    select case (trim(rad_absorber))
    case ('ql')
      if (.not. lmoist) call fail("rad_absorber = 'ql' (liquid water) requires lmoist = .true.")
      if (rad_absorber_sv /= 0) call fail("rad_absorber_sv is set but rad_absorber = 'ql'; set rad_absorber = 'sv'")
    case ('sv')
      if (rad_absorber_sv < 1 .or. rad_absorber_sv > nsv) then
        write(msg, '(a,i0)') "rad_absorber = 'sv' requires 1 <= rad_absorber_sv <= nsv; got ", rad_absorber_sv
        call fail(msg)
      end if
    case default
      call fail("rad_absorber must be 'ql' or 'sv'; got '" // trim(rad_absorber) // "'")
    end select
    if (rad_ls .and. myid == 0) then
      write(0, *) 'WARNING: rad_ls ignored; the thlpcar column of lscale.inp is always applied in forces'
    end if
  end subroutine check_radiation_options

  !> Abort on all ranks with one message from rank 0.
  subroutine fail(msg)
    use modmpi, only : myid, abort_all_ranks
    implicit none
    character(*), intent(in) :: msg
    if (myid == 0) write(0, '(2a)') 'ERROR: ', trim(msg)
    call abort_all_ranks(1)
  end subroutine fail

  !> Recompute thlprad at the first RK substep when due (generic radpar, then the
  !> case hooks); add it to thlp every substep.
  subroutine radiation
    use modglobal, only : ib, ie, jb, je, kb, ke, timee, rk3step
    use modfields, only : thlp
    use modcases,  only : call_hooks_radiation
    implicit none
    if (iradiation == irad_none) return
    ! timee accumulates dt, so it can fall an ulp short of a grid time k timerad: both the test
    ! and the grid index allow a relative tolerance (1e-6 of timerad; 0 for timerad = 0, every step)
    if (rk3step == 1 .and. timee >= tnextrad - 1.e-6 * timerad) then
      ! the next update on the timerad grid (multiples of timerad) after timee: a warm start
      ! (tnextrad = 0) updates at its first step, thlprad not being in the restart files, and then
      ! rejoins the grid of the run it continues instead of updating at every step until tnextrad
      ! has caught up with timee
      if (timerad > 0.) tnextrad = timerad * (1. + floor(timee / timerad + 1.e-6))
      thlprad = 0.
      if (rad_longw) call radpar
      call call_hooks_radiation
    end if
    thlp(ib:ie,jb:je,kb:ke) = thlp(ib:ie,jb:je,kb:ke) + thlprad
  end subroutine radiation

  !> DALES dev 075ac4e modradiation.f90:radpar, longwave branch.
  subroutine radpar
    use modglobal, only : ib, ie, jb, je, kb, ke, cp, dzf
    use modfields, only : ql0, sv0, rhof, exnf, IIc
    implicit none
    real :: lwpt(kb:ke+1), lwpb(kb:ke+1), absorber(kb:ke)
    integer :: i, j, k

    do j = jb, je
      do i = ib, ie
        if (trim(rad_absorber) == 'sv') then
          absorber = sv0(i,j,kb:ke,rad_absorber_sv) * IIc(i,j,kb:ke)
        else
          absorber = ql0(i,j,kb:ke) * IIc(i,j,kb:ke)
        end if
        lwpt(ke+1) = 0.
        do k = ke, kb, -1
          lwpt(k) = lwpt(k+1) + rhof(k) * absorber(k) * dzf(k)
        end do
        lwpb(kb) = 0.
        ! DALES writes rhof(k)*absorber(k-1)*dzf(k-1); rhof(k-1) is the consistent
        ! density of the cell whose path is added (deliberate, tiny difference).
        do k = kb + 1, ke + 1
          lwpb(k) = lwpb(k-1) + rhof(k-1) * absorber(k-1) * dzf(k-1)
        end do
        do k = kb, ke + 1
          lwd(i,j,k) = dlwtop * exp(-rka * lwpt(k))
          lwu(i,j,k) = dlwbot * exp(-rka * lwpb(k))
        end do
        do k = kb, ke
          thlprad(i,j,k) = - ((lwd(i,j,k+1) + lwu(i,j,k+1)) - (lwd(i,j,k) + lwu(i,j,k))) &
                           / (rhof(k) * cp * exnf(k) * dzf(k))
        end do

        thlprad(i,j,kb:ke) = thlprad(i,j,kb:ke) * IIc(i,j,kb:ke)
      end do
    end do
  end subroutine radpar

  subroutine exitradiation
    implicit none
    if (allocated(thlprad)) deallocate(thlprad, lwd, lwu)
  end subroutine exitradiation

end module modradiation
