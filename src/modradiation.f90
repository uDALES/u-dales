!> Parametrised longwave radiation, ported from DALES dev 075ac4e
!> src/modradiation.f90 (radiation, radpar longwave branch) and the
!> DYCOMS rad_user of cases/dycoms_rf02/moduser.f90, which here is the
!> built-in option lrad_dycoms.
module modradiation
  use modraddata
  implicit none
  private
  public :: initradiation, radiation, exitradiation

contains

  subroutine initradiation
    use modglobal, only : ifnamopt, fname_options, ib, ie, jb, je, kb, ke, nsv, lmoist, ltempeq
    use modmpi,    only : myid, comm3d, mpierr, my_real, mpi_logical, mpi_integer
    implicit none
    integer :: ierr
    namelist/NAMRADIATION/ iradiation, timerad, rad_longw, rad_shortw, rad_smoke, rad_ls, &
         lrad_dycoms, isvsmoke, rka, dlwtop, dlwbot, rad_dycoms_D, rad_dycoms_qtthres

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
    call MPI_BCAST(rad_smoke,   1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(rad_ls,      1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(lrad_dycoms, 1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(isvsmoke,    1, MPI_INTEGER, 0, comm3d, mpierr)
    call MPI_BCAST(rka,         1, MY_REAL,     0, comm3d, mpierr)
    call MPI_BCAST(dlwtop,      1, MY_REAL,     0, comm3d, mpierr)
    call MPI_BCAST(dlwbot,      1, MY_REAL,     0, comm3d, mpierr)
    call MPI_BCAST(rad_dycoms_D,       1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(rad_dycoms_qtthres, 1, MY_REAL, 0, comm3d, mpierr)

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
    if (rad_smoke .and. (isvsmoke < 1 .or. isvsmoke > nsv)) then
      write(0, *) 'ERROR: rad_smoke requires 1 <= isvsmoke <= nsv'
      stop 1
    end if
    if ((.not. rad_smoke) .and. (.not. lmoist)) then
      write(0, *) 'ERROR: liquid-water radiation (rad_smoke = .false.) requires lmoist = .true.'
      stop 1
    end if
    if (lrad_dycoms .and. (.not. lmoist)) then
      write(0, *) 'ERROR: lrad_dycoms requires lmoist = .true.'
      stop 1
    end if
    if (rad_ls .and. myid == 0) then
      write(0, *) 'WARNING: rad_ls ignored; the thlpcar column of lscale.inp is always applied in forces'
    end if
  end subroutine initradiation

  !> Recompute thlprad at the first RK substep when due; add it to thlp every substep.
  subroutine radiation
    use modglobal, only : ib, ie, jb, je, kb, ke, timee, rk3step
    use modfields, only : thlp
    implicit none
    if (iradiation == irad_none) return
    if (rk3step == 1 .and. timee >= tnextrad) then
      tnextrad = tnextrad + timerad
      thlprad = 0.
      if (rad_longw) call radpar
    end if
    thlp(ib:ie,jb:je,kb:ke) = thlp(ib:ie,jb:je,kb:ke) + thlprad
  end subroutine radiation

  !> DALES dev 075ac4e modradiation.f90:radpar, longwave branch, plus the
  !> DYCOMS free-troposphere term (moduser rad_user of dycoms_rf02).
  subroutine radpar
    use modglobal, only : ib, ie, jb, je, kb, ke, cp, dzf, zf
    use modfields, only : ql0, qt0, sv0, rhof, exnf, IIc
    implicit none
    real :: lwpt(kb:ke+1), lwpb(kb:ke+1), absorber(kb:ke), zi
    integer :: i, j, k

    do j = jb, je
      do i = ib, ie
        if (rad_smoke) then
          absorber = sv0(i,j,kb:ke,isvsmoke) * IIc(i,j,kb:ke)
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

        if (lrad_dycoms) then
          ! Local inversion height: highest crossing of the qt threshold. Two fixes
          ! relative to DALES rad_user: zi is reset per column (DALES reuses the
          ! previous column's value) and the i loop covers the full range (DALES ran
          ! i = 2, j1). Sign: COOLING, -D/3[...]; DALES adds +D/3[...], which doubles
          ! subsidence warming instead of balancing it. With thl = 297.5 + (z-zi)^(1/3),
          ! -D/3[(z-zi)^(1/3) + zi (z-zi)^(-2/3)] = -(D/3) z (z-zi)^(-2/3)
          ! = w dthl/dz for w = -D z, i.e. exactly minus the subsidence warming.
          ! rho_i/rho ~ 1 and no 1/exner, as in DALES.
          zi = -1.
          do k = ke - 1, kb, -1
            if (qt0(i,j,k) > rad_dycoms_qtthres) then
              zi = zf(k) + (rad_dycoms_qtthres - qt0(i,j,k)) * (zf(k+1) - zf(k)) / (qt0(i,j,k+1) - qt0(i,j,k))
              exit
            end if
          end do
          if (zi > 0.) then
            do k = ke, kb, -1
              if (zf(k) <= zi) exit
              thlprad(i,j,k) = thlprad(i,j,k) - rad_dycoms_D / 3. * &
                               ((zf(k) - zi)**(1./3.) + zi * (zf(k) - zi)**(-2./3.))
            end do
          end if
        end if

        thlprad(i,j,kb:ke) = thlprad(i,j,kb:ke) * IIc(i,j,kb:ke)
      end do
    end do
  end subroutine radpar

  subroutine exitradiation
    implicit none
    if (allocated(thlprad)) deallocate(thlprad, lwd, lwu)
  end subroutine exitradiation

end module modradiation
