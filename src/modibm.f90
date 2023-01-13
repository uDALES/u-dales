!!> \file modibm.f90
!!!  adds forcing terms for immersed boundaries
!
!>
!!  \author Jasper Thomas TU Delft / Ivo Suter Imperial College London
!
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
module modibm
   use modibmdata
   !use wf_uno
   implicit none
   save
   public :: initibm, ibmnorm, ibmwallfun, bottom, &
             nsolpts_u, nsolpts_v, nsolpts_w, nsolpts_c, &
             nbndpts_u, nbndpts_v, nbndpts_w, nbndpts_c, &
             nfctsecs_u, nfctsecs_v, nfctsecs_w, nfctsecs_c

   ! u
   integer :: nsolpts_u
   integer, allocatable :: solpts_u(:,:)
   logical, allocatable :: lsolptsrank_u(:) !
   integer, allocatable :: solptsrank_u(:) ! indices of points on current rank
   integer :: nsolptsrank_u

   integer :: nbndpts_u
   integer, allocatable :: bndpts_u(:,:) ! ijk location of fluid boundary point
   real, allocatable    :: recpts_u(:,:) ! xyz location of reconstruction point
   real, allocatable    :: bndvec_u(:,:) ! vector from bound point to rec point (normalised)
   real, allocatable    :: bnddst_u(:) ! distance between surface & bound point
   real, allocatable    :: recdst_u(:) ! distance between bound point & rec point
   !logical, allocatable :: lbndptsrank_u(:) !
   !integer, allocatable :: bndptsrank_u(:) ! indices of points on current rank
   logical, allocatable :: lrec_u(:) ! Switch whether reconstruction point is a computational point
   !integer :: nbndptsrank_u

   integer :: nfctsecs_u
   integer, allocatable :: secbndptids_u(:)
   integer, allocatable :: secfacids_u(:)
   real,    allocatable :: secareas_u(:)

   !logical, allocatable :: lfctsecsrank_u(:)
   integer, allocatable :: fctsecsrank_u(:)
   integer :: nfctsecsrank_u

   ! v
   integer :: nsolpts_v
   integer, allocatable :: solpts_v(:,:)
   logical, allocatable :: lsolptsrank_v(:) !
   integer, allocatable :: solptsrank_v(:) ! indices of points on current rank
   integer :: nsolptsrank_v

   integer :: nbndpts_v
   integer, allocatable :: bndpts_v(:,:) ! ijk location of fluid boundary point
   real, allocatable    :: recpts_v(:,:) ! xyz location of reconstruction point
   real, allocatable    :: bndvec_v(:,:) ! vector from bound point to rec point (normalised)
   real, allocatable    :: bnddst_v(:) ! distance between surface & bound point
   real, allocatable    :: recdst_v(:) ! distance between bound point & rec point
   !logical, allocatable :: lbndptsrank_v(:) !
   !integer, allocatable :: bndptsrank_v(:) ! indices of points on current rank
   logical, allocatable :: lrec_v(:) ! Switch whether reconstruction point is a computational point
   !integer :: nbndptsrank_v

   integer :: nfctsecs_v
   integer, allocatable :: secbndptids_v(:)
   integer, allocatable :: secfacids_v(:)
   real,    allocatable :: secareas_v(:)

   !logical, allocatable :: lfctsecsrank_v(:)
   integer, allocatable :: fctsecsrank_v(:)
   integer :: nfctsecsrank_v

   ! w
   integer :: nsolpts_w
   integer, allocatable :: solpts_w(:,:)
   logical, allocatable :: lsolptsrank_w(:) !
   integer, allocatable :: solptsrank_w(:) !
   integer :: nsolptsrank_w

   integer :: nbndpts_w
   integer, allocatable :: bndpts_w(:,:) ! ijk location of fluid boundary point
   real, allocatable    :: recpts_w(:,:) ! xyz location of reconstruction point
   real, allocatable    :: bndvec_w(:,:) ! vector from bound point to rec point (normalised)
   real, allocatable    :: bnddst_w(:) ! distance between surface & bound point
   real, allocatable    :: recdst_w(:) ! distance between bound point & rec point
   !logical, allocatable :: lbndptsrank_w(:) !
   !integer, allocatable :: bndptsrank_w(:) ! indices of points on current rank
   logical, allocatable :: lrec_w(:) ! Switch whether reconstruction point is a computational point
   !integer :: nbndptsrank_w

   integer :: nfctsecs_w
   integer, allocatable :: secbndptids_w(:)
   integer, allocatable :: secfacids_w(:)
   real,    allocatable :: secareas_w(:)

   !logical, allocatable :: lfctsecsrank_w(:)
   integer, allocatable :: fctsecsrank_w(:)
   integer :: nfctsecsrank_w

   ! c
   integer :: nsolpts_c
   integer, allocatable :: solpts_c(:,:)
   logical, allocatable :: lsolptsrank_c(:) !
   integer, allocatable :: solptsrank_c(:) !
   integer :: nsolptsrank_c

   integer :: nbndpts_c
   integer, allocatable :: bndpts_c(:,:) ! ijk location of fluid boundary point
   real, allocatable    :: recpts_c(:,:) ! xyz location of reconstruction point
   real, allocatable    :: bndvec_c(:,:) ! vector from bound point to rec point (normalised)
   real, allocatable    :: bnddst_c(:) ! distance between surface & bound point
   real, allocatable    :: recdst_c(:) ! distance between bound point & rec point
   !logical, allocatable :: lbndptsrank_w(:) !
   !integer, allocatable :: bndptsrank_w(:) ! indices of points on current rank
   logical, allocatable :: lrec_c(:) ! Switch whether reconstruction point is a computational point
   !integer :: nbndptsrank_w

   integer :: nfctsecs_c
   integer, allocatable :: secbndptids_c(:)
   integer, allocatable :: secfacids_c(:)
   real,    allocatable :: secareas_c(:)

   !logical, allocatable :: lfctsecsrank_w(:)
   integer, allocatable :: fctsecsrank_c(:)
   integer :: nfctsecsrank_c

   contains

   subroutine initibm
     use modglobal, only : libm, xh, xf, yh, yf, zh, zf, xhat, yhat, zhat

     if (.not. libm) return

     call initibmnorm('solid_u.txt', nsolpts_u, solpts_u, nsolptsrank_u, solptsrank_u)
     call initibmnorm('solid_v.txt', nsolpts_v, solpts_v, nsolptsrank_v, solptsrank_v)
     call initibmnorm('solid_w.txt', nsolpts_w, solpts_w, nsolptsrank_w, solptsrank_w)
     call initibmnorm('solid_c.txt', nsolpts_c, solpts_c, nsolptsrank_c, solptsrank_c)

     call initibmwallfun('fluid_boundary_u.txt', 'facet_sections_u.txt', xh, yf, zf, &
                         nbndpts_u, bndpts_u, bnddst_u, recpts_u, recdst_u, bndvec_u, lrec_u, &
                         nfctsecs_u, secfacids_u, secareas_u, secbndptids_u, nfctsecsrank_u, fctsecsrank_u)

    call initibmwallfun('fluid_boundary_v.txt', 'facet_sections_v.txt', xf, yh, zf, &
                        nbndpts_v, bndpts_v, bnddst_v, recpts_v, recdst_v, bndvec_v, lrec_v, &
                        nfctsecs_v, secfacids_v, secareas_v, secbndptids_v, nfctsecsrank_v, fctsecsrank_v)

    call initibmwallfun('fluid_boundary_w.txt', 'facet_sections_w.txt', xf, yf, zh, &
                        nbndpts_w, bndpts_w, bnddst_w, recpts_w, recdst_w, bndvec_w, lrec_w, &
                        nfctsecs_w, secfacids_w, secareas_w, secbndptids_w, nfctsecsrank_w, fctsecsrank_w)

    call initibmwallfun('fluid_boundary_c.txt', 'facet_sections_c.txt', xf, yf, zf, &
                        nbndpts_c, bndpts_c, bnddst_c, recpts_c, recdst_c, bndvec_c, lrec_c, &
                        nfctsecs_c, secfacids_c, secareas_c, secbndptids_c, nfctsecsrank_c, fctsecsrank_c)

   end subroutine initibm


   subroutine initibmnorm(fname, nsolpts, solpts, nsolptsrank, solptsrank)
     use modglobal, only : ifinput
     use modmpi,    only : myid, comm3d, MPI_INTEGER, mpierr
     use decomp_2d, only : zstart, zend

     integer, intent(in)  :: nsolpts
     integer, intent(out) :: nsolptsrank
     integer, intent(out), dimension(:,:), allocatable :: solpts
     integer, intent(out), dimension(:), allocatable :: solptsrank

     logical :: lsolptsrank(nsolpts)
     integer n, m
     character(11) fname
     character(80) chmess

     allocate(solpts(nsolpts,3))

     ! read u points
     if (myid == 0) then
       open (ifinput, file=fname)
       read (ifinput, '(a80)') chmess
       do n = 1, nsolpts
         read (ifinput, *) solpts(n,1), solpts(n,2), solpts(n,3)
       end do
       close (ifinput)
     end if

     call MPI_BCAST(solpts, nsolpts*3, MPI_INTEGER, 0, comm3d, mpierr)

     ! Determine whether points are on this rank
     nsolptsrank = 0
     do n = 1, nsolpts
       if ((solpts(n,1) >= zstart(1) .and. solpts(n,1) <= zend(1)) .and. &
          (solpts(n,2) >= zstart(2) .and. solpts(n,2) <= zend(2))) then
          lsolptsrank(n) = .true.
          nsolptsrank = nsolptsrank + 1
        else
          lsolptsrank(n) = .false.
       end if
     end do

     ! Store indices of points on current rank - only loop through these points
     allocate(solptsrank(nsolptsrank))
     m = 0
     do n = 1, nsolpts
       if (lsolptsrank(n)) then
          m = m + 1
          solptsrank(m) = n
       end if
     end do

     write(*,*) "rank ", myid, " has ", nsolptsrank, " solid points from ", fname

   end subroutine initibmnorm


   subroutine initibmwallfun(fname_bnd, fname_sec, x, y, z, nbndpts, bndpts, bnddst, recpts, recdst, bndvec, lrec, &
                             nfctsecs, secfacids, secareas, secbndptids, nfctsecsrank, fctsecsrank)
     use modglobal, only : ifinput, ib, itot, ih, jb, jtot, jh, kb, ktot, kh
     use modmpi,    only : myid, comm3d, MPI_INTEGER, MY_REAL, mpierr
     use decomp_2d, only : zstart, zend

     character(20), intent(in) :: fname_bnd, fname_sec
     integer, intent(in) :: nbndpts, nfctsecs
     real,    intent(in), dimension(ib:itot+ih) :: x
     real,    intent(in), dimension(jb:jtot+jh) :: y
     real,    intent(in), dimension(kb:ktot+kh) :: z
     integer, intent(out) :: nfctsecsrank
     real,    intent(out), dimension(:),   allocatable :: bnddst, recdst, secareas
     real,    intent(out), dimension(:,:), allocatable :: recpts, bndvec
     integer, intent(out), dimension(:),   allocatable :: secfacids, secbndptids, fctsecsrank
     integer, intent(out), dimension(:,:), allocatable :: bndpts
     logical, intent(out), dimension(:),   allocatable :: lrec

     logical, dimension(nbndpts) :: lbndptsrank
     logical, dimension(nfctsecs) :: lfctsecsrank
     integer i, j, k, n, m, nbndptsrank
     character(80) chmess

     allocate(bndpts(nbndpts,3))
     allocate(bnddst(nbndpts))
     allocate(recpts(nbndpts,3))
     allocate(bndvec(nbndpts,3))
     allocate(recdst(nbndpts))
     allocate(lrec(nbndpts))

     ! read u points
     if (myid == 0) then
       open (ifinput, file=fname_bnd)
       read (ifinput, '(a80)') chmess
       do n = 1, nbndpts
         read (ifinput, *) bndpts(n,1), bndpts(n,2), bndpts(n,3), bnddst(n), recpts(n,1), recpts(n,2), recpts(n,3)
       end do
       close (ifinput)

       ! Calculate vector
       do n = 1,nbndpts
         bndvec(n,1) = recpts(n,1) - x(bndpts(n,1))
         bndvec(n,2) = recpts(n,2) - y(bndpts(n,2))
         bndvec(n,3) = recpts(n,3) - z(bndpts(n,3))
         recdst(n) = norm2(bndvec(n,:))
         bndvec(n,:) = bndvec(n,:) / recdst(n)
       end do
     end if

     call MPI_BCAST(bndpts, nbndpts*3, MPI_INTEGER, 0, comm3d, mpierr)
     call MPI_BCAST(recpts, nbndpts*3, MY_REAL,     0, comm3d, mpierr)
     call MPI_BCAST(bndvec, nbndpts*3, MY_REAL,     0, comm3d, mpierr)
     call MPI_BCAST(bnddst, nbndpts,   MY_REAL,     0, comm3d, mpierr)
     call MPI_BCAST(recdst, nbndpts,   MY_REAL,     0, comm3d, mpierr)

     ! Determine whether points are on this rank
     nbndptsrank = 0
     do n = 1, nbndpts
       if ((bndpts(n,1) >= zstart(1) .and. bndpts(n,1) <= zend(1)) .and. &
          (bndpts(n,2) >= zstart(2) .and. bndpts(n,2) <= zend(2))) then
          lbndptsrank(n) = .true.
          nbndptsrank = nbndptsrank + 1
        else
          lbndptsrank(n) = .false.
       end if
     end do

     write(*,*) "rank ", myid, " has ", nbndptsrank, "points from ", fname_bnd


     ! Work out whether reconstruction point is a computational point
     ! This should only occur when the normal is in the x direction and the distance between them is dx
     ! This means the facet section is aligned
     do n = 1,nbndpts
       ! if (computational point)
          lrec(n) = .true.
       ! else
          !lrec(n) = .false.
       ! end if
     end do

     allocate(secfacids(nfctsecs))
     allocate(secareas(nfctsecs))
     allocate(secbndptids(nfctsecs))

     ! read u facet sections
     if (myid == 0) then
       open (ifinput, file=fname_sec)
       read (ifinput, '(a80)') chmess
       do n = 1, nfctsecs
         read (ifinput, *) secfacids(n), secareas(n), secbndptids(n)
       end do
       close (ifinput)
     end if

     call MPI_BCAST(secfacids,   nfctsecs, MPI_INTEGER, 0, comm3d, mpierr)
     call MPI_BCAST(secareas,    nfctsecs, MY_REAL,     0, comm3d, mpierr)
     call MPI_BCAST(secbndptids, nfctsecs, MPI_INTEGER, 0, comm3d, mpierr)

     ! Determine whether section needs to be updated by this rank
     nfctsecsrank = 0
     do n = 1, nfctsecs
       if (lbndptsrank(secbndptids(n))) then
          lfctsecsrank(n) = .true.
          nfctsecsrank = nfctsecsrank + 1
        else
          lfctsecsrank(n) = .false.
       end if
     end do

     ! Store indices of sections on current rank - only loop through these sections
     allocate(fctsecsrank(nfctsecsrank))
     m = 0
     do n = 1, nfctsecs
       if (lfctsecsrank(n)) then
          m = m + 1
          fctsecsrank(m) = n
       end if
     end do

     write(*,*) "rank ", myid, " has ", nfctsecsrank, " facet sections from ", fname_sec

   end subroutine initibmwallfun


   subroutine ibmnorm
     use modglobal,   only : libm, ltempeq
     use modfields,   only : um, vm, wm, thlm, up, vp, wp, thlp
     use modboundary, only : halos
     use decomp_2d,   only : zstart, zend
     use modmpi, only : myid

     integer i, j, k, n, m

     if (.not. libm) return

     ! Set internal velocities to zero
     call solid(nsolpts_u, solpts_u, nsolptsrank_u, solptsrank_u, um, up, 0.)
     call solid(nsolpts_v, solpts_v, nsolptsrank_v, solptsrank_v, vm, vp, 0.)
     call solid(nsolpts_w, solpts_w, nsolptsrank_w, solptsrank_w, wm, wp, 0.)

     ! scalars
     if (ltempeq) call solid(nsolpts_c, solpts_c, nsolptsrank_c, solptsrank_c, thlm, thlp, 288.)

   end subroutine ibmnorm


   subroutine solid(nsolpts, solpts, nsolptsrank, solptsrank, var, rhs, val)
     use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh
     use decomp_2d, only : zstart

     real, intent(in) :: val
     integer, intent(in) :: nsolpts, nsolptsrank
     integer, intent(in) :: solpts(nsolpts,3), solptsrank(nsolptsrank)
     real, intent(inout) :: var(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
     real, intent(inout) :: rhs(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)

     integer :: i, j, k, n, m

     do m=1,nsolptsrank
      n = solptsrank(m)
       !if (lsolptsrank_u(n)) then
         i = solpts(n,1) - zstart(1) + 1
         j = solpts(n,2) - zstart(2) + 1
         k = solpts(n,3) - zstart(3) + 1
         var(i,j,k) = val
         rhs(i,j,k) = 0.
       !end if
     end do

   end subroutine


   subroutine ibmwallfun
     use modglobal, only : libm, iwallmom, xhat, yhat, zhat
     use modfields, only : up, vp, wp, tau_x, tau_y, tau_z

      if (.not. libm) return

      if (iwallmom == 1) return

      ! momfluxb = 0.
      ! tfluxb = 0.
      ! qfluxb = 0.
      !call uwallfun
      tau_x = up
      call wallfunmom(xhat, up, nbndpts_u, bndpts_u, bnddst_u, recpts_u, recdst_u, bndvec_u, lrec_u, &
                       nfctsecs_u, secfacids_u, secareas_u, secbndptids_u, nfctsecsrank_u, fctsecsrank_u)
      tau_x = up - tau_x

      tau_y = vp
      call wallfunmom(yhat, vp, nbndpts_v, bndpts_v, bnddst_v, recpts_v, recdst_v, bndvec_v, lrec_v, &
                      nfctsecs_v, secfacids_v, secareas_v, secbndptids_v, nfctsecsrank_v, fctsecsrank_v)
      tau_y = vp - tau_y

      tau_z = wp
      call wallfunmom(zhat, wp, nbndpts_w, bndpts_w, bnddst_w, recpts_w, recdst_w, bndvec_w, lrec_w, &
                      nfctsecs_w, secfacids_w, secareas_w, secbndptids_w, nfctsecsrank_w, fctsecsrank_w)
      tau_z = wp - tau_z

    end subroutine ibmwallfun


   subroutine wallfunmom(dir, rhs, nbndpts, bndpts, bnddst, recpts, recdst, bndvec, lrec, &
                         nfctsecs, secfacids, secareas, secbndptids, nfctsecsrank, fctsecsrank)
     use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh, eps1, fkar, dx, dy, dzf, iwallmom, xhat, yhat, zhat
     use modfields, only : u0, v0, w0, thl0, momfluxb, tfluxb, qfluxb
     use initfac,   only : facT, facz0, facz0h
     use modsurfdata, only : z0, z0h
     use decomp_2d, only : zstart

     real, intent(in)    :: dir(3)
     real, intent(inout) :: rhs(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
     integer, intent(in) :: nbndpts, nfctsecs, nfctsecsrank
     real,    intent(in), dimension(nbndpts)   :: bnddst, recdst, secareas
     real,    intent(in), dimension(nbndpts,3) :: recpts, bndvec
     integer, intent(in), dimension(nfctsecs)  :: secfacids, secbndptids, fctsecsrank
     integer, intent(in), dimension(nbndpts,3) :: bndpts
     logical, intent(in), dimension(nbndpts)   :: lrec

     integer i, j, k, n, m, sec, pt, fac
     real dist, stress, area, vol, momvol, Tair, Tsurf
     real, dimension(3) :: uvec, norm

     do m = 1,nfctsecsrank
       sec = fctsecsrank(m) ! index of section
       n = secbndptids(sec) ! index of boundary point
       fac = secfacids(sec) ! index of facet

       ! if surface normal is aligned with coordinate direction, there cannot be any tangential component
       if ((abs(abs(dot_product(bndvec(n,:), dir)) - 1.)) < eps1) cycle

       i = bndpts(n,1) - zstart(1) + 1 ! should be on this rank!
       j = bndpts(n,2) - zstart(2) + 1 ! should be on this rank!
       k = bndpts(n,3) - zstart(3) + 1 ! should be on this rank!
       if ((i < ib) .or. (i > ie) .or. (j < jb) .or. (j > je)) write(*,*) "problem", i, j

       if (lrec(n)) then ! Section aligned with grid - don't interpolate, use this point's velocity
         if (all(abs(dir - xhat) < eps1)) then
           uvec(1) = u0(i,j,k)
           uvec(2) = 0.25 * (v0(i,j,k) + v0(i,j+1,k) + v0(i-1,j,k) + v0(i-1,j+1,k))
           uvec(3) = 0.25 * (w0(i,j,k) + w0(i,j,k+1) + w0(i-1,j,k) + w0(i-1,j,k+1)) !only for equidistant grid!
           if (iwallmom == 2) then
             Tair = 0.5 * (thl0(i,j,k) + thl0(i-1,j,k))
           end if
         else if (all(abs(dir - yhat) < eps1)) then
           uvec(1) = 0.25 * (u0(i,j,k) + u0(i+1,j,k) + u0(i,j-1,k) + u0(i+1,j-1,k))
           uvec(2) = v0(i,j,k)
           uvec(3) = 0.25 * (w0(i,j,k) + w0(i,j,k+1) + w0(i,j-1,k) + w0(i,j-1,k+1)) !only for equidistant grid!
           if (iwallmom == 2) then
             Tair = 0.5 * (thl0(i,j,k) + thl0(i,j-1,k))
           end if
         else if (all(abs(dir - zhat) < eps1)) then
           uvec(1) = 0.25 * (u0(i,j,k) + u0(i+1,j,k) + u0(i,j-1,k) + u0(i+1,j-1,k))
           uvec(2) = 0.25 * (v0(i,j,k) + v0(i,j+1,k) + v0(i-1,j,k) + v0(i-1,j+1,k))
           uvec(3) = w0(i,j,k)
           if (iwallmom == 2) then
             Tair = 0.5 * (thl0(i,j,k) + thl0(i,j,k-1))
           end if
         else
           write(*,*) 'ERROR: none of x,y, or z'
           stop 1
         end if
         dist = bnddst(n)
       else ! Interpolate velocities at reconstruction point
         write(0, *) 'ERROR: interpolation at reconstruction point not supported'
         stop 1
       end if

       if (all(abs(uvec) < eps1)) cycle

       ! local coordinate system
       norm(1) = bndvec(n,1)
       norm(2) = bndvec(n,2)
       norm(3) = bndvec(n,3)

       if (iwallmom == 2) then
         stress = calc_stress(uvec, dir, norm, dist, facz0(fac), facz0h(fac), Tair, facT(fac, 1))
       else if (iwallmom == 3) then
         stress = calc_stress(uvec, dir, norm, dist, facz0(fac))
       end if

       area = secareas(sec)
       vol = dx*dy*dzf(k)
       momvol = stress * area / vol
       rhs(i,j,k) = rhs(i,j,k) - momvol
       !stressout(i,j,k) = - momvol

     end do

  end subroutine wallfunmom


   real function calc_stress(uvec, dir, norm, dist, z0, z0h, Tair, Tsurf) ! make interface which excludes z0h and Tsurf for neutral cases
     use modglobal, only : fkar, iwallmom, xhat, yhat, zhat, eps1

     real, intent(in) :: dist, z0
     real, intent(in), optional :: z0h, Tair, Tsurf
     real, dimension(3), intent(in) :: uvec, dir, norm

     real :: utan, a_is, a_xn, a_yn, a_zn, ctm, stress_ix, stress_iy, stress_iz, stress
     real, dimension(3) :: strm, span, stressvec

     span = cross_product(norm, uvec)
     if (norm2(span) < eps1) then
       ! velocity is pointing into or outof the surface, so no tangential stress
       calc_stress = 0.
       return
     else
       span = span / norm2(span)
     end if
     strm = cross_product(span, norm)
     utan = dot_product(uvec, strm)

     ! Rotation from local (strm,span,norm) to global (xhat,yhat,zhat) basis
     ! \tau'_ij = a_ip a_jq \tau_pq
     ! \tau_pq in local coordinates is something like \tau \delta_13, because we only have \tau_{strm,norm})
     a_is = dot_product(dir, strm)
     a_xn = dot_product(xhat, norm)
     a_yn = dot_product(yhat, norm)
     a_zn = dot_product(zhat, norm)

     if (iwallmom == 2) then ! convective
       ctm = unom(dist, z0, z0h, Tsurf, Tair, utan)
     else if (iwallmom == 3) then ! neutral
       ctm = (fkar / log(dist / z0))**2 ! z0 depends on facet in general
     end if

     stress = ctm * utan**2

     stress_ix = a_is * a_xn * stress
     stress_iy = a_is * a_yn * stress
     stress_iz = a_is * a_zn * stress

     stressvec(1) = stress_ix
     stressvec(2) = stress_iy
     stressvec(3) = stress_iz
     calc_stress  = SIGN(norm2(stressvec), utan)

   end function calc_stress

   function cross_product(a,b)
     implicit none
     real, dimension(3) :: cross_product, a, b

     cross_product(1) = a(2)*b(3) - a(3)*b(2)
     cross_product(2) =-a(3)*b(1) + a(1)*b(3)
     cross_product(3) = a(1)*b(2) - a(2)*b(1)

   end function cross_product

   real function unom(dist, z0, z0h, Tsurf, Tair, utan)
     use modglobal, only : grav, prandtlmol, fkar

      IMPLICIT NONE
      REAL, INTENT(in) :: dist, z0, z0h, Tsurf, Tair, utan
      REAL :: Ribl1, Fm, Fh, cm, ch, Ctm, M
      REAL, PARAMETER :: b1 = 9.4 !parameters from Uno1995
      REAL, PARAMETER :: b2 = 4.7
      REAL, PARAMETER :: dm = 7.4
      REAL, PARAMETER :: dh = 5.3
      real :: dT, Ribl0, logdz, logdzh, logzh, sqdz, fkar2

      dT = Tair - Tsurf
      Ribl0 = grav * dist * dT / (Tsurf * utan) !Eq. 6, guess initial Ri

      logdz = LOG(dist/z0)
      logdzh = LOG(dist/z0h)
      logzh = LOG(z0/z0h)
      sqdz = SQRT(dist/z0)
      fkar2 = fkar**2

      IF (Ribl0 > 0.21) THEN !0.25 approx critical for bulk Richardson number  => stable
         Fm = 1./(1. + b2*Ribl0)**2 !Eq. 4
         Fh = Fm !Eq. 4
      ELSE ! => unstable
         cm = (dm*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
         ch = (dh*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
         Fm = 1. - (b1*Ribl0)/(1. + cm*SQRT(ABS(Ribl0))) !Eq. 3
         Fh = 1. - (b1*Ribl0)/(1. + ch*SQRT(ABS(Ribl0))) !Eq. 3
      END IF

      M = prandtlmol*logdz*SQRT(Fm)/Fh !Eq. 14

      Ribl1 = Ribl0 - Ribl0*prandtlmol*logzh/(prandtlmol*logzh + M) !Eq. 17

      !interate to get new Richardson number
      IF (Ribl1 > 0.21) THEN !0.25 approx critical for bulk Richardson number  => stable
         Fm = 1./(1. + b2*Ribl1)**2 !Eq. 4
      ELSE ! => unstable
         cm = (dm*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
         Fm = 1. - (b1*Ribl1)/(1. + cm*SQRT(ABS(Ribl1))) !Eq. 3
      END IF

      Ctm = fkar2/(logdz**2)*Fm !Eq. 7
      unom = Ctm !Eq. 2, Eq. 8
   END FUNCTION unom

   ! !for scalar
   ! !FUNCTION unoh(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !for heat, the bit that does not change no matter what wall
   ! SUBROUTINE unoh(otf, octh, logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !for heat, the bit that does not change no matter what wall
   ! !flux in Km/s
   !    IMPLICIT NONE
   !    REAL, INTENT(in) :: logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2
   !    REAL, INTENT(out) :: octh, otf
   !    REAL :: Ribl1, Fm, Fh, cm, ch, M, dTrough, cth
   !    REAL, PARAMETER :: b1 = 9.4 !parameters from Uno1995
   !    REAL, PARAMETER :: b2 = 4.7
   !    REAL, PARAMETER :: dm = 7.4
   !    REAL, PARAMETER :: dh = 5.3
   !    REAL, PARAMETER :: prandtlmol = 0.71
   !    REAL, PARAMETER :: prandtlmoli = 1/0.71
   !
   !    octh = 0.
   !    otf = 0.
   !    IF (Ribl0 > 0.21) THEN !0.25 approx critical for bulk Richardson number  => stable
   !       Fm = 1./(1. + b2*Ribl0)**2 !Eq. 4
   !       Fh = Fm !Eq. 4
   !    ELSE ! => unstable
   !       cm = (dm*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
   !       ch = (dh*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
   !       Fm = 1. - (b1*Ribl0)/(1. + cm*SQRT(ABS(Ribl0))) !Eq. 3
   !       Fh = 1. - (b1*Ribl0)/(1. + ch*SQRT(ABS(Ribl0))) !Eq. 3
   !    END IF
   !
   !    M = prandtlmol*logdz*SQRT(Fm)/Fh !Eq. 14
   !
   !    Ribl1 = Ribl0 - Ribl0*prandtlmol*logzh/(prandtlmol*logzh + M) !Eq. 17
   !
   !    !interate to get new Richardson number
   !    IF (Ribl1 > 0.21) THEN !0.25 approx critical for bulk Richardson number  => stable
   !       Fm = 1./(1. + b2*Ribl1)**2 !Eq. 4
   !       Fh = Fm !Eq. 4
   !    ELSE ! => unstable
   !       cm = (dm*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
   !       ch = (dh*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
   !       Fm = 1. - (b1*Ribl1)/(1. + cm*SQRT(ABS(Ribl1))) !Eq. 3
   !       Fh = 1. - (b1*Ribl1)/(1. + ch*SQRT(ABS(Ribl1))) !Eq. 3
   !    END IF
   !    M = prandtlmol*logdz*SQRT(Fm)/Fh !Eq. 14
   !
   !    dTrough = dT*1./(prandtlmol*logzh/M + 1.) !Eq. 13a
   !
   !    octh = SQRT(utangInt)*fkar2/(logdz*logdzh)*prandtlmoli*Fh !Eq. 8
   !    otf = octh*dTrough !Eq. 2, Eq. 8
   !
   ! END SUBROUTINE unoh

   subroutine bottom
      !kind of obsolete when road facets are being used
      !vegetated floor not added (could simply be copied from vegetated horizontal facets)
      use modglobal, only:ib, ie, ih, jh, kh, jb, je, kb, numol, prandtlmol, dzh, nsv, &
         dxf, dxhi, dzf, dzfi, numoli, ltempeq, khc, lmoist, BCbotT, BCbotq, BCbotm, BCbots, dzh2i
      use modfields, only : u0,v0,e120,um,vm,w0,wm,e12m,thl0,qt0,sv0,thlm,qtm,svm,up,vp,thlp,qtp,svp,shear,momfluxb,tfluxb,cth
      use modsurfdata, only:thlflux, qtflux, svflux, ustar, thvs, wtsurf, wqsurf, thls, z0, z0h
      use modsubgriddata, only:ekm, ekh
      use modmpi, only:myid
      implicit none
      integer :: i, j, jp, jm, m
      !momentum
      if (BCbotm.eq.2) then
      call wfuno(ih, jh, kh, up, vp, thlp, momfluxb, tfluxb, cth, bcTfluxA, u0, v0, thl0, thls, z0, z0h, 0, 1, 91)
      elseif (BCbotm.eq.3) then
      call wfmneutral(ih, jh, kh, up, vp, momfluxb, u0, v0, z0, 0, 1, 91)
      else
      write(0, *) "ERROR: bottom boundary type for momentum undefined"
      stop 1
      end if

      if (ltempeq) then
         if (BCbotT.eq.1) then !neumann/fixed flux bc for temperature
            do j = jb, je
               do i = ib, ie
                  thlp(i, j, kb) = thlp(i, j, kb) &
                                   + ( &
                                   0.5*(dzf(kb - 1)*ekh(i, j, kb) + dzf(kb)*ekh(i, j, kb - 1)) &
                                   *(thl0(i, j, kb) - thl0(i, j, kb - 1)) &
                                   *dzh2i(kb) &
                                   - wtsurf &
                                   )*dzfi(kb)
               end do
            end do
         else if (BCbotT.eq.2) then !wall function bc for temperature (fixed temperature)
            call wfuno(ih, jh, kh, up, vp, thlp, momfluxb, tfluxb, cth, bcTfluxA, u0, v0, thl0, thls, z0, z0h, 0, 1, 92)
         else
         write(0, *) "ERROR: bottom boundary type for temperature undefined"
         stop 1
         end if
      end if ! ltempeq

      if (lmoist) then
         if (BCbotq.eq.1) then !neumann/fixed flux bc for moisture
            do j = jb, je
               do i = ib, ie
                  qtp(i, j, kb) = qtp(i, j, kb) + ( &
                                  0.5*(dzf(kb - 1)*ekh(i, j, kb) + dzf(kb)*ekh(i, j, kb - 1)) &
                                  *(qt0(i, j, kb) - qt0(i, j, kb - 1)) &
                                  *dzh2i(kb) &
                                  + wqsurf &
                                  )*dzfi(kb)
               end do
            end do
         else
          write(0, *) "ERROR: bottom boundary type for moisture undefined"
          stop 1
         end if !
      end if !lmoist

      if (nsv>0) then
         if (BCbots.eq.1) then !neumann/fixed flux bc for moisture
            do j = jb, je
               do i = ib, ie
                  do m = 1, nsv
                      svp(i, j, kb, m) = svp(i, j, kb, m) + ( &
                                      0.5*(dzf(kb - 1)*ekh(i, j, kb) + dzf(kb)*ekh(i, j, kb - 1)) &
                                     *(sv0(i, j, kb, m) - sv0(i, j, kb - 1, m)) &
                                     *dzh2i(kb) &
                                     + 0. &
                                     )*dzfi(kb)
                  end do
               end do
            end do
         else
          write(0, *) "ERROR: bottom boundary type for scalars undefined"
          stop 1
         end if !
      end if

      e120(:, :, kb - 1) = e120(:, :, kb)
      e12m(:, :, kb - 1) = e12m(:, :, kb)
      ! wm(:, :, kb) = 0. ! SO moved to modboundary
      ! w0(:, :, kb) = 0.
      return
   end subroutine bottom

end module modibm
