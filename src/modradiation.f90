!> Parametrised longwave radiation, ported from DALES dev 075ac4e
!> src/modradiation.f90 (radiation, radpar longwave branch). Case-specific
!> terms (DYCOMS free troposphere) are case hooks, see modcases.
module modradiation
  use modraddata
  implicit none
  private
  public :: initradiation, radiation, exitradiation

contains

  subroutine initradiation
    use modglobal, only : ifnamopt, fname_options, ib, ie, jb, je, kb, ke, nsv, lmoist, ltempeq
    use modmpi,    only : myid, comm3d, mpierr, my_real, mpi_logical, mpi_integer, mpi_character
    implicit none
    integer :: ierr
    namelist/NAMRADIATION/ iradiation, timerad, rad_longw, rad_shortw, rad_ls, &
         rad_absorber, rad_absorber_sv, rka, dlwtop, dlwbot

    if (myid == 0) then
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read (ifnamopt, NAMRADIATION, iostat=ierr)
      if (ierr > 0) then
        write(0, *) 'ERROR: Problem in namoptions NAMRADIATION'
        write(0, *) 'iostat error: ', ierr
        stop 1
      end if
      close(ifnamopt)
    end if
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

    allocate(thlprad(ib:ie,jb:je,kb:ke)); thlprad = 0.
    allocate(lwd(ib:ie,jb:je,kb:ke+1));   lwd = 0.
    allocate(lwu(ib:ie,jb:je,kb:ke+1));   lwu = 0.
    if (iradiation == irad_none) return

    if (iradiation /= irad_par) then
      write(0, *) 'ERROR: iradiation must be 0 (none) or 2 (parametrised); got ', iradiation
      stop 1
    end if
    if (rad_shortw) then
      write(0, *) 'ERROR: rad_shortw is not supported in uDALES'
      stop 1
    end if
    if (.not. ltempeq) then
      write(0, *) 'ERROR: iradiation = 2 requires ltempeq = .true.'
      stop 1
    end if
    select case (rad_absorber)
    case ('ql')
      if (.not. lmoist) then
        write(0, *) "ERROR: rad_absorber = 'ql' (liquid water) requires lmoist = .true."
        stop 1
      end if
    case ('sv')
      if (rad_absorber_sv < 1 .or. rad_absorber_sv > nsv) then
        write(0, *) "ERROR: rad_absorber = 'sv' requires 1 <= rad_absorber_sv <= nsv; got ", rad_absorber_sv
        stop 1
      end if
    case default
      write(0, *) "ERROR: rad_absorber must be 'ql' or 'sv'; got '", rad_absorber, "'"
      stop 1
    end select
    if (rad_ls .and. myid == 0) then
      write(0, *) 'WARNING: rad_ls ignored; the thlpcar column of lscale.inp is always applied in forces'
    end if
  end subroutine initradiation

  !> Recompute thlprad at the first RK substep when due (generic radpar, then the
  !> case hooks); add it to thlp every substep.
  subroutine radiation
    use modglobal, only : ib, ie, jb, je, kb, ke, timee, rk3step
    use modfields, only : thlp
    use modcases,  only : call_hooks_radiation
    implicit none
    if (iradiation == irad_none) return
    if (rk3step == 1 .and. timee >= tnextrad) then
      tnextrad = tnextrad + timerad
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
        if (rad_absorber == 'sv') then
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
