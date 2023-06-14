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
   public :: initibm, ibmnorm, ibmwallfun, bottom, lbottom, createmasks, &
             nsolpts_u, nsolpts_v, nsolpts_w, nsolpts_c, &
             nbndpts_u, nbndpts_v, nbndpts_w, nbndpts_c, &
             nfctsecs_u, nfctsecs_v, nfctsecs_w, nfctsecs_c, &
             mask_u, mask_v, mask_w, mask_c

    abstract interface
      function interp_velocity(i,j,k)
        real :: interp_velocity(3)
        integer, intent(in) :: i, j, k
      end function interp_velocity
    end interface

    abstract interface
      real function interp_temperature(i,j,k)
        integer, intent(in) :: i, j, k
      end function interp_temperature
    end interface

   logical :: lbottom = .false.

   ! read from namoptions
   integer :: nsolpts_u, nsolpts_v, nsolpts_w, nsolpts_c, &
              nbndpts_u, nbndpts_v, nbndpts_w, nbndpts_c, &
              nfctsecs_u, nfctsecs_v, nfctsecs_w, nfctsecs_c

   real, allocatable, target, dimension(:,:,:) :: mask_u, mask_v, mask_w, mask_c

   TYPE solid_info_type
     integer :: nsolpts
     integer, allocatable :: solpts(:,:)
     logical, allocatable :: lsolptsrank(:) !
     integer, allocatable :: solptsrank(:) ! indices of points on current rank
     integer :: nsolptsrank
   end TYPE solid_info_type

   type(solid_info_type) :: solid_info_u, solid_info_v, solid_info_w, solid_info_c

   TYPE bound_info_type
     integer :: nbndpts
     integer, allocatable :: bndpts(:,:) ! ijk location of fluid boundary point
     real, allocatable    :: intpts(:,:) ! xyz location of boundary intercept point
     real, allocatable    :: bndvec(:,:) ! vector from boundary to fluid point (normalised)
     real, allocatable    :: recpts(:,:) ! xyz location of reconstruction point
     integer, allocatable :: recids_u(:,:) ! ijk location of u grid cell that rec point is in
     integer, allocatable :: recids_v(:,:) ! ijk location of u grid cell that rec point is in
     integer, allocatable :: recids_w(:,:) ! ijk location of u grid cell that rec point is in
     integer, allocatable :: recids_c(:,:) ! ijk location of u grid cell that rec point is in
     real, allocatable    :: bnddst(:) ! distance between surface & bound point
     integer, allocatable :: bndptsrank(:) ! indices of points on current rank
     logical, allocatable :: lcomprec(:) ! Switch whether reconstruction point is a computational point
     logical, allocatable :: lskipsec(:) ! Switch whether to skip finding the shear stress at this point
     integer :: nbndptsrank

     integer :: nfctsecs
     integer, allocatable :: secbndptids(:)
     integer, allocatable :: secfacids(:)
     real,    allocatable :: secareas(:)
     integer, allocatable :: fctsecsrank(:)
     integer :: nfctsecsrank
   end TYPE bound_info_type

   type(bound_info_type) :: bound_info_u, bound_info_v, bound_info_w, bound_info_c

   contains

   subroutine initibm
     use modglobal, only : libm, xh, xf, yh, yf, zh, zf, xhat, yhat, zhat, vec0, &
                           ib, ie, ih, ihc, jb, je, jh, jhc, kb, ke, kh, khc, nsv, &
                           iwallmom, lmoist, ltempeq
     use decomp_2d, only : exchange_halo_z

     real, allocatable :: rhs(:,:,:)

     if (.not. libm) return

     solid_info_u%nsolpts = nsolpts_u
     solid_info_v%nsolpts = nsolpts_v
     solid_info_w%nsolpts = nsolpts_w
     call initibmnorm('solid_u.txt', solid_info_u)
     call initibmnorm('solid_v.txt', solid_info_v)
     call initibmnorm('solid_w.txt', solid_info_w)

     ! Define (real) masks
     ! Hopefully this can be removed eventually if (integer) IIx halos can be communicated
     ! These are only used in modibm, to cancel subgrid term across solid boundaries
     allocate(mask_u(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)); mask_u = 1.
     allocate(mask_v(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)); mask_v = 1.
     allocate(mask_w(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)); mask_w = 1.
     mask_w(:,:,kb) = 0.     ! In future this shouldn't be needed?
     mask_u(:,:,kb-kh) = 0.
     mask_v(:,:,kb-kh) = 0.
     mask_w(:,:,kb-kh) = 0.

     allocate(rhs(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
     call solid(solid_info_u, mask_u, rhs, 0., ih, jh, kh)
     call solid(solid_info_v, mask_v, rhs, 0., ih, jh, kh)
     call solid(solid_info_w, mask_w, rhs, 0., ih, jh, kh)
     call exchange_halo_z(mask_u)!, opt_zlevel=(/ih,jh,0/))
     call exchange_halo_z(mask_v)!, opt_zlevel=(/ih,jh,0/))
     call exchange_halo_z(mask_w)!, opt_zlevel=(/ih,jh,0/))

     if (iwallmom > 1) then
       bound_info_u%nbndpts = nbndpts_u
       bound_info_v%nbndpts = nbndpts_v
       bound_info_w%nbndpts = nbndpts_w
       bound_info_u%nfctsecs = nfctsecs_u
       bound_info_v%nfctsecs = nfctsecs_v
       bound_info_w%nfctsecs = nfctsecs_w
       call initibmwallfun('fluid_boundary_u.txt', 'facet_sections_u.txt', xhat, bound_info_u)
       call initibmwallfun('fluid_boundary_v.txt', 'facet_sections_v.txt', yhat, bound_info_v)
       call initibmwallfun('fluid_boundary_w.txt', 'facet_sections_w.txt', zhat, bound_info_w)
     end if

     if (ltempeq .or. lmoist .or. nsv>0) then
       solid_info_c%nsolpts = nsolpts_c
       call initibmnorm('solid_c.txt', solid_info_c)

       bound_info_c%nbndpts = nbndpts_c
       bound_info_c%nfctsecs = nfctsecs_c
       call initibmwallfun('fluid_boundary_c.txt', 'facet_sections_c.txt', vec0, bound_info_c)

       allocate(mask_c(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)); mask_c = 1.
       mask_c(:,:,kb-kh) = 0.
       call solid(solid_info_c, mask_c, rhs, 0., ih, jh, kh)
       call exchange_halo_z(mask_c)!, opt_zlevel=(/ih,jh,0/))
     end if

     deallocate(rhs)

   end subroutine initibm


   subroutine initibmnorm(fname, solid_info)
     use modglobal, only : ifinput
     use modmpi,    only : myid, comm3d, MPI_INTEGER, mpierr
     use decomp_2d, only : zstart, zend

     character(11), intent(in) :: fname

     type(solid_info_type), intent(inout) :: solid_info

     logical :: lsolptsrank(solid_info%nsolpts)
     integer n, m

     character(80) chmess

     allocate(solid_info%solpts(solid_info%nsolpts,3))

     ! read u points
     if (myid == 0) then
       open (ifinput, file=fname)
       read (ifinput, '(a80)') chmess
       do n = 1, solid_info%nsolpts
         read (ifinput, *) solid_info%solpts(n,1), solid_info%solpts(n,2), solid_info%solpts(n,3)
       end do
       close (ifinput)
     end if

     call MPI_BCAST(solid_info%solpts, solid_info%nsolpts*3, MPI_INTEGER, 0, comm3d, mpierr)

     ! Determine whether points are on this rank
     solid_info%nsolptsrank = 0
     do n = 1, solid_info%nsolpts
       if ((solid_info%solpts(n,1) >= zstart(1) .and. solid_info%solpts(n,1) <= zend(1)) .and. &
           (solid_info%solpts(n,2) >= zstart(2) .and. solid_info%solpts(n,2) <= zend(2))) then
          lsolptsrank(n) = .true.
          solid_info%nsolptsrank = solid_info%nsolptsrank + 1
        else
          lsolptsrank(n) = .false.
       end if
     end do

     ! Store indices of points on current rank - only loop through these points
     allocate(solid_info%solptsrank(solid_info%nsolptsrank))
     m = 0
     do n = 1, solid_info%nsolpts
       if (lsolptsrank(n)) then
          m = m + 1
          solid_info%solptsrank(m) = n
       end if
     end do

     !write(*,*) "rank ", myid, " has ", solid_info%nsolptsrank, " solid points from ", fname

   end subroutine initibmnorm


   subroutine initibmwallfun(fname_bnd, fname_sec, dir, bound_info)
     use modglobal, only : ifinput, ib, ie, itot, ih, jb, je, jtot, jh, kb, ktot, kh, &
                           xf, yf, zf, xh, yh, zh, dx, dy, dzh, dzf, xhat, yhat, zhat
     use modmpi,    only : myid, comm3d, MPI_INTEGER, MY_REAL, MPI_LOGICAL, mpierr
     use initfac,   only : facnorm
     use decomp_2d, only : zstart, zend

     character(20), intent(in) :: fname_bnd, fname_sec
     type(bound_info_type) :: bound_info
     real, intent(in), dimension(3) :: dir
     real, dimension(ib:itot+ih) :: xgrid
     real, dimension(jb:jtot+jh) :: ygrid
     real, dimension(kb:ktot+kh) :: zgrid
     logical, dimension(bound_info%nbndpts)  :: lbndptsrank
     logical, dimension(bound_info%nfctsecs) :: lfctsecsrank
     real, dimension(3) :: norm, p0, p1, pxl, pxu, pyl, pyu, pzl, pzu
     integer, dimension(6) :: check
     integer, dimension(1) :: pos_min_dist
     real, dimension(6,3) :: inter
     real, dimension(6) :: inter_dists
     real :: xc, yc, zc, xl, yl, zl, xu, yu, zu, checkxl, checkxu, checkyl, checkyu, checkzl, checkzu, inter_dist
     integer i, j, k, n, m, norm_align, dir_align, pos, p
     real dst

     character(80) chmess

     allocate(bound_info%bndpts(bound_info%nbndpts,3))

     ! read u points
     if (myid == 0) then
       open (ifinput, file=fname_bnd)
       read (ifinput, '(a80)') chmess
       do n = 1, bound_info%nbndpts
         read (ifinput, *) bound_info%bndpts(n,1), bound_info%bndpts(n,2), bound_info%bndpts(n,3)
       end do
       close (ifinput)
     end if

     call MPI_BCAST(bound_info%bndpts, bound_info%nbndpts*3, MPI_INTEGER, 0, comm3d, mpierr)

     ! Determine whether points are on this rank
     bound_info%nbndptsrank = 0
     do n = 1, bound_info%nbndpts
       if ((bound_info%bndpts(n,1) >= zstart(1) .and. bound_info%bndpts(n,1) <= zend(1)) .and. &
           (bound_info%bndpts(n,2) >= zstart(2) .and. bound_info%bndpts(n,2) <= zend(2))) then
          lbndptsrank(n) = .true.
          bound_info%nbndptsrank = bound_info%nbndptsrank + 1
        else
          lbndptsrank(n) = .false.
       end if
     end do

     !write(*,*) "rank ", myid, " has ", bound_info%nbndptsrank, "points from ", fname_bnd

     ! Store indices of points on current rank - only loop through these points
     allocate(bound_info%bndptsrank(bound_info%nbndptsrank))
     m = 0
     do n = 1, bound_info%nbndpts
       if (lbndptsrank(n)) then
          i = bound_info%bndpts(n,1) - zstart(1) + 1
          j = bound_info%bndpts(n,2) - zstart(2) + 1
          k = bound_info%bndpts(n,3) - zstart(3) + 1
          if ((i < ib) .or. (i > ie) .or. (j < jb) .or. (j > je)) then
            write(*,*) "problem in initibmwallfun", i, j
            stop 1
          end if
          m = m + 1
          bound_info%bndptsrank(m) = n
       end if
     end do

     allocate(bound_info%secfacids(bound_info%nfctsecs))
     allocate(bound_info%secareas(bound_info%nfctsecs))
     allocate(bound_info%secbndptids(bound_info%nfctsecs))
     allocate(bound_info%intpts(bound_info%nfctsecs,3))
     allocate(bound_info%bnddst(bound_info%nfctsecs))
     allocate(bound_info%bndvec(bound_info%nfctsecs,3))
     allocate(bound_info%recpts(bound_info%nfctsecs,3))
     allocate(bound_info%recids_u(bound_info%nfctsecs,3))
     allocate(bound_info%recids_v(bound_info%nfctsecs,3))
     allocate(bound_info%recids_w(bound_info%nfctsecs,3))
     allocate(bound_info%recids_c(bound_info%nfctsecs,3))
     allocate(bound_info%lcomprec(bound_info%nfctsecs))
     allocate(bound_info%lskipsec(bound_info%nfctsecs))

     dir_align = alignment(dir)
     select case(dir_align)
     case(1)
       xgrid = xh
       ygrid = yf
       zgrid = zf
     case(2)
       xgrid = xf
       ygrid = yh
       zgrid = zf
     case(3)
       xgrid = xf
       ygrid = yf
       zgrid = zh
     case(0)
       xgrid = xf
       ygrid = yf
       zgrid = zf
     end select

     if (myid == 0) then
       open (ifinput, file=fname_sec)
       read (ifinput, '(a80)') chmess
       do n = 1, bound_info%nfctsecs
         read (ifinput, *) bound_info%secfacids(n), bound_info%secareas(n), bound_info%secbndptids(n), bound_info%bnddst(n)
                           !bound_info%intpts(n,1),  bound_info%intpts(n,2), bound_info%intpts(n,3)
       end do
       close (ifinput)

       ! Calculate vector
       do n = 1,bound_info%nfctsecs
         m = bound_info%secbndptids(n)
         !bound_info%bndvec(n,1) = xgrid(bound_info%bndpts(m,1)) - bound_info%intpts(n,1)
         !bound_info%bndvec(n,2) = ygrid(bound_info%bndpts(m,2)) - bound_info%intpts(n,2)
         !bound_info%bndvec(n,3) = zgrid(bound_info%bndpts(m,3)) - bound_info%intpts(n,3)
         !bound_info%bnddst(n) = norm2(bound_info%bndvec(n,:))
         !write(*,*) bound_info%bnddst(n)
         !bound_info%bndvec(n,:) = bound_info%bndvec(n,:) / bound_info%bnddst(n)

         norm = facnorm(bound_info%secfacids(n),:)
         norm_align = alignment(norm)

         if (norm_align /= 0) then ! this facet normal is in x, y, or z direction
           bound_info%lcomprec(n) = .true. ! simple reconstruction
           if (dir_align == norm_align) then ! the normal is in the same direction
             ! don't need to calculate shear stress as no tangential component
             bound_info%lskipsec(n) = .true.
           else
             bound_info%lskipsec(n) = .false.
           end if

         else ! need to reconstruct
           bound_info%lcomprec(n) = .false.
           bound_info%lskipsec(n) = .false.

           ! cell centre (of current grid)
           xc = xgrid(bound_info%bndpts(m,1))
           yc = ygrid(bound_info%bndpts(m,2))
           zc = zgrid(bound_info%bndpts(m,3))

           ! cell edges
           xl = xc - dx/2.
           xu = xc + dx/2.
           yl = yc - dy/2.
           yu = yc + dy/2.
           zl = zc - dzf(1)/2. ! assumes equidistant
           zu = zc + dzf(1)/2. ! assumes equidistant

           ! points on planes
           pxl = (/xl, yc, zc/)
           pxu = (/xu, yc, zc/)
           pyl = (/xc, yl, zc/)
           pyu = (/xc, yu, zc/)
           pzl = (/xc, yc, zl/)
           pzu = (/xc, yc, zu/)

           p0 = (/xc, yc, zc/)
           p1 = p0 + norm * sqrt(3.)*(dx*dy*dzf(1))**(1./3.)

           call plane_line_intersection(xhat, pxl, p0, p1, inter(1,:), check(1), inter_dists(1))
           call plane_line_intersection(xhat, pxu, p0, p1, inter(2,:), check(2), inter_dists(2))
           call plane_line_intersection(yhat, pyl, p0, p1, inter(3,:), check(3), inter_dists(3))
           call plane_line_intersection(yhat, pyu, p0, p1, inter(4,:), check(4), inter_dists(4))
           call plane_line_intersection(zhat, pzl, p0, p1, inter(5,:), check(5), inter_dists(5))
           call plane_line_intersection(zhat, pzu, p0, p1, inter(6,:), check(6), inter_dists(6))

           pos_min_dist = minloc(inter_dists, mask=check==1)
           pos = pos_min_dist(1)

           if (pos == 0) then
             write(*,*) "ERROR: no intersection found"
           else
             bound_info%recpts(n,:) = inter(pos,:) ! x y z
           end if

           ! find which cell the point lies in for
           ! findloc will return 0 when the point is on the lower domain edge(?)
           ! because e.g. xf = [0.5 0.5+dx 0.5+2*dx ... Lx]
           ! but the reconstruction point could be [0 y z]
           ! so when these ids are used we should include a point beyond this lower boundary,
           ! i.e. xf0 = [-0.5 0.5 0.5+dx ... Lx]
           bound_info%recids_u(n,1) = findloc(bound_info%recpts(n,1) >= xh, .true., 1, back=.true.)
           bound_info%recids_u(n,2) = findloc(bound_info%recpts(n,2) >= yf, .true., 1, back=.true.)
           bound_info%recids_u(n,3) = findloc(bound_info%recpts(n,3) >= zf, .true., 1, back=.true.)

           bound_info%recids_v(n,1) = findloc(bound_info%recpts(n,1) >= xf, .true., 1, back=.true.)
           bound_info%recids_v(n,2) = findloc(bound_info%recpts(n,2) >= yh, .true., 1, back=.true.)
           bound_info%recids_v(n,3) = findloc(bound_info%recpts(n,3) >= zf, .true., 1, back=.true.)

           bound_info%recids_w(n,1) = findloc(bound_info%recpts(n,1) >= xf, .true., 1, back=.true.)
           bound_info%recids_w(n,2) = findloc(bound_info%recpts(n,2) >= yf, .true., 1, back=.true.)
           bound_info%recids_w(n,3) = findloc(bound_info%recpts(n,3) >= zh, .true., 1, back=.true.)

           bound_info%recids_c(n,1) = findloc(bound_info%recpts(n,1) >= xf, .true., 1, back=.true.)
           bound_info%recids_c(n,2) = findloc(bound_info%recpts(n,2) >= yf, .true., 1, back=.true.)
           bound_info%recids_c(n,3) = findloc(bound_info%recpts(n,3) >= zf, .true., 1, back=.true.)

           ! Add check to see if recids are all available to current rank.

           !recpts should lie inside the box defined by these corners
           ! u
           if ((bound_info%recpts(n,1) < xh(bound_info%recids_u(n,1))) .or. &
               (bound_info%recpts(n,1) > xh(bound_info%recids_u(n,1)+1))) then
             write(*,*) "ERROR: x out of bounds"
             stop 1
           end if
           if ((bound_info%recpts(n,2) < yf(bound_info%recids_u(n,2))) .or. &
               (bound_info%recpts(n,2) > yf(bound_info%recids_u(n,2)+1))) then
             write(*,*) "ERROR: y out of bounds"
             stop 1
           end if
           if ((bound_info%recpts(n,3) < zf(bound_info%recids_u(n,3))) .or. &
               (bound_info%recpts(n,3) > zf(bound_info%recids_u(n,3)+1))) then
             write(*,*) "ERROR: z out of bounds"
             stop 1
           end if

           ! v
           if ((bound_info%recpts(n,1) < xf(bound_info%recids_v(n,1))) .or. &
               (bound_info%recpts(n,1) > xf(bound_info%recids_v(n,1)+1))) then
             write(*,*) "ERROR: x out of bounds"
             stop 1
           end if
           if ((bound_info%recpts(n,2) < yh(bound_info%recids_v(n,2))) .or. &
               (bound_info%recpts(n,2) > yh(bound_info%recids_v(n,2)+1))) then
             write(*,*) "ERROR: y out of bounds"
             stop 1
           end if
           if ((bound_info%recpts(n,3) < zf(bound_info%recids_v(n,3))) .or. &
               (bound_info%recpts(n,3) > zf(bound_info%recids_v(n,3)+1))) then
             write(*,*) "ERROR: z out of bounds"
             stop 1
           end if

           ! w
           if ((bound_info%recpts(n,1) < xf(bound_info%recids_w(n,1))) .or. &
               (bound_info%recpts(n,1) > xf(bound_info%recids_w(n,1)+1))) then
             write(*,*) "ERROR: x out of bounds"
             stop 1
           end if
           if ((bound_info%recpts(n,2) < yf(bound_info%recids_w(n,2))) .or. &
               (bound_info%recpts(n,2) > yf(bound_info%recids_w(n,2)+1))) then
             write(*,*) "ERROR: y out of bounds"
             stop 1
           end if
           if ((bound_info%recpts(n,3) < zh(bound_info%recids_w(n,3))) .or. &
               (bound_info%recpts(n,3) > zh(bound_info%recids_w(n,3)+1))) then
             write(*,*) "ERROR: z out of bounds"
             stop 1
           end if

         end if
       end do
     end if ! myid==0

     call MPI_BCAST(bound_info%secfacids,   bound_info%nfctsecs,   MPI_INTEGER, 0, comm3d, mpierr)
     call MPI_BCAST(bound_info%secareas,    bound_info%nfctsecs,   MY_REAL,     0, comm3d, mpierr)
     call MPI_BCAST(bound_info%secbndptids, bound_info%nfctsecs,   MPI_INTEGER, 0, comm3d, mpierr)
     call MPI_BCAST(bound_info%intpts,      bound_info%nfctsecs*3, MY_REAL,     0, comm3d, mpierr)
     call MPI_BCAST(bound_info%bndvec,      bound_info%nfctsecs*3, MY_REAL,     0, comm3d, mpierr)
     call MPI_BCAST(bound_info%bnddst,      bound_info%nfctsecs,   MY_REAL,     0, comm3d, mpierr)
     call MPI_BCAST(bound_info%recpts,      bound_info%nfctsecs*3, MY_REAL,     0, comm3d, mpierr)
     call MPI_BCAST(bound_info%recids_u,    bound_info%nfctsecs*3, MPI_INTEGER, 0, comm3d, mpierr)
     call MPI_BCAST(bound_info%recids_v,    bound_info%nfctsecs*3, MPI_INTEGER, 0, comm3d, mpierr)
     call MPI_BCAST(bound_info%recids_w,    bound_info%nfctsecs*3, MPI_INTEGER, 0, comm3d, mpierr)
     call MPI_BCAST(bound_info%recids_c,    bound_info%nfctsecs*3, MPI_INTEGER, 0, comm3d, mpierr)
     call MPI_BCAST(bound_info%lskipsec,    bound_info%nfctsecs,   MPI_LOGICAL, 0, comm3d, mpierr)
     call MPI_BCAST(bound_info%lcomprec,    bound_info%nfctsecs,   MPI_LOGICAL, 0, comm3d, mpierr)

     ! Determine whether section needs to be updated by this rank
     bound_info%nfctsecsrank = 0
     do n = 1, bound_info%nfctsecs
       if (lbndptsrank(bound_info%secbndptids(n))) then
          lfctsecsrank(n) = .true.
          bound_info%nfctsecsrank = bound_info%nfctsecsrank + 1
        else
          lfctsecsrank(n) = .false.
       end if
     end do

     ! Store indices of sections on current rank - only loop through these sections
     allocate(bound_info%fctsecsrank(bound_info%nfctsecsrank))
     m = 0
     do n = 1, bound_info%nfctsecs
       if (lfctsecsrank(n)) then
          m = m + 1
          bound_info%fctsecsrank(m) = n
       end if
     end do

   end subroutine initibmwallfun


   subroutine plane_line_intersection(norm, V0, P0, P1, I, check, dist)
     use modglobal, only : vec0, eps1
     implicit none
     ! determines the intersection of a plane and a line segment
     ! norm: plane normal
     ! V0: point on the plane
     ! P0: start of line segment
     ! P1: end of line segment
     ! I: intersection point
     ! dist: distance from P0 to intersection point
     ! check: 0 if no intersection
     !        1 if unique intersection
     !        2 if line segment is in the plane
     !        3 if intersection is outside line segment
     real, intent(in),  dimension(3) :: norm, V0, P0, P1
     real, intent(out), dimension(3) :: I
     integer, intent(out) :: check
     real, intent(out) :: dist
     real, dimension(3) :: u, w
     real :: D, N, sI

     I = vec0
     w = P0 - V0
     u = P1 - P0
     D = dot_product(norm, u)
     N =-dot_product(norm, w)

     if (abs(D) < eps1) then ! line orthogonal to plane normal -> segment parallel to plane
       if (abs(N) < eps1) then ! start point is on the plane -> segment lies in the plane
         check = 2
         return
       else
         check = 0
         return
       end if
     end if

     sI = N / D
     I = P0 + sI * u
     dist = norm2(I - P0)

     if ((sI < 0.) .or. (sI > 1.)) then
       check = 3
     else
       check = 1
     end if

   end subroutine plane_line_intersection


   subroutine ibmnorm
     use modglobal,   only : ih, jh, kh, ihc, jhc, khc, nsv, dzf, zh, kb, ke, kh, nsv, libm, ltempeq, lmoist
     use modfields,   only : um, vm, wm, thlm, qtm, svm, up, vp, wp, thlp, qtp, svp, thl0, qt0, sv0, thl0av
     use modboundary, only : halos
     use decomp_2d,   only : zstart, zend
     use modmpi, only : myid

     integer i, j, k, n, m

     if (.not. libm) return

     ! Set internal velocities to zero
     call solid(solid_info_u, um, up, 0., ih, jh, kh)
     call solid(solid_info_v, vm, vp, 0., ih, jh, kh)
     call solid(solid_info_w, wm, wp, 0., ih, jh, kh)

     ! scalars
     if (ltempeq) then
        ! Solid value should not matter - choose domain-average for visualization.
        call solid(solid_info_c, thlm, thlp, sum(thl0av(kb:ke)*dzf(kb:ke))/zh(ke+1), ih, jh, kh)
        call advecc2nd_corr_liberal(thl0, thlp)
     end if

     if (lmoist) then
       call solid(solid_info_c, qtm, qtp, 0., ih, jh, kh)
       call advecc2nd_corr_liberal(qt0, qtp)
    end if

    do n=1,nsv
      call solid(solid_info_c, svm(:,:,:,n), svp(:,:,:,n), 0., ihc, jhc, khc)
      call solid_boundary(bound_info_c, mask_c, svm(:,:,:,n), svp(:,:,:,n), ihc, jhc, khc)
      ! The above should be replaced with something like avecc2nd_corr,
      ! but using the kappa advection scheme rather than the second order one.
    end do

   end subroutine ibmnorm


   subroutine solid(solid_info, var, rhs, val, hi, hj, hk)
     use modglobal, only : ib, ie, jb, je, kb, ke
     use decomp_2d, only : zstart

     type(solid_info_type), intent(in) :: solid_info
     integer, intent(in) :: hi, hj, hk
     real, intent(inout) :: var(ib-hi:ie+hi,jb-hj:je+hj,kb-hk:ke+hk)
     real, intent(inout) :: rhs(ib-hi:ie+hi,jb-hj:je+hj,kb   :ke+hk)
     real, intent(in) :: val

     integer :: i, j, k, n, m

     do m=1,solid_info%nsolptsrank
      n = solid_info%solptsrank(m)
       !if (lsolptsrank_u(n)) then
         i = solid_info%solpts(n,1) - zstart(1) + 1
         j = solid_info%solpts(n,2) - zstart(2) + 1
         k = solid_info%solpts(n,3) - zstart(3) + 1
         var(i,j,k) = val
         rhs(i,j,k) = 0.
       !end if
     end do

   end subroutine solid


   subroutine solid_boundary(bound_info, mask, var, rhs, hi, hj, hk)
     use modglobal, only : eps1, ib, ie, ih, jb, je, jh, kb, ke, kh
     use decomp_2d, only : zstart
     ! uDALES 1 approach
     ! Not truly conservative
     type(bound_info_type), intent(in) :: bound_info
     integer, intent(in) :: hi, hj, hk
     real, intent(in)    :: mask(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
     real, intent(inout) :: var(ib-hi:ie+hi,jb-hj:je+hj,kb-hk:ke+hk)
     real, intent(inout) :: rhs(ib-hi:ie+hi,jb-hj:je+hj,kb   :ke+hk)

     integer i, j, k, n, m

     do m = 1,bound_info%nbndptsrank
       n = bound_info%bndptsrank(m)
       i = bound_info%bndpts(n,1) - zstart(1) + 1
       j = bound_info%bndpts(n,2) - zstart(2) + 1
       k = bound_info%bndpts(n,3) - zstart(3) + 1

       if (abs(mask(i,j+1,k)) < eps1) then
         ! rhs(i,j+1,k) = 0.
         rhs(i,j+1,k) = rhs(i,j,k)
         var(i,j+1,k) = var(i,j,k)
       end if

       if (abs(mask(i,j-1,k)) < eps1) then
         ! rhs(i,j-1,k) = 0.
         rhs(i,j-1,k) = rhs(i,j,k)
         var(i,j-1,k) = var(i,j,k)
       end if

       if (abs(mask(i,j,k+1)) < eps1) then
         ! rhs(i,j,k+1) = 0.
         rhs(i,j,k+1) = rhs(i,j,k)
         var(i,j,k+1) = var(i,j,k)
       end if

       if (abs(mask(i,j,k-1)) < eps1) then
         ! rhs(i,j,k-1) = 0.
         rhs(i,j,k-1) = rhs(i,j,k)
         var(i,j,k-1) = var(i,j,k)
       end if

       if (abs(mask(i+1,j,k)) < eps1) then
         ! rhs(i+1,j,k) = 0.
         rhs(i+1,j,k) = rhs(i,j,k)
         var(i+1,j,k) = var(i,j,k)
       end if

       if (abs(mask(i-1,j,k)) < eps1) then
         ! rhs(i-1,j,k) = 0.
         rhs(i-1,j,k) = rhs(i,j,k)
         var(i-1,j,k) = var(i,j,k)
       end if

     end do

   end subroutine


   subroutine advecc2nd_corr_conservative(var, rhs)
     ! Removes the advection contribution from solid velocities, which should be
     ! close to zero but are not necessarily due to pressure correction.
     ! Has a fairly drastic effect on the initial flow, but the scalar is
     ! conserved throughout the simulation.
     use modglobal,      only : eps1, ib, ie, ih, jb, je, jh, kb, ke, kh, &
                                dx2i, dxi5, dy2i, dyi5, dzf, dzh2i, dzfi, dzhi, dzfi5
     use modfields,      only : u0, v0, w0
     use modsubgriddata, only : ekh
     use decomp_2d,      only : zstart

     real, intent(in)    :: var(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
     real, intent(inout) :: rhs(ib-ih:ie+ih,jb-jh:je+jh,kb   :ke+kh)
     integer :: i, j, k, n, m

     do m = 1,bound_info_c%nbndptsrank
      n = bound_info_c%bndptsrank(m)
         i = bound_info_c%bndpts(n,1) - zstart(1) + 1
         j = bound_info_c%bndpts(n,2) - zstart(2) + 1
         k = bound_info_c%bndpts(n,3) - zstart(3) + 1

         if ((abs(mask_u(i+1,j,k)) < eps1) .or. (abs(mask_c(i+1,j,k)) < eps1)) then
           rhs(i,j,k) = rhs(i,j,k) + u0(i+1,j,k)*(var(i+1,j,k) + var(i,j,k))*dxi5
         end if

         if ((abs(mask_u(i,j,k)) < eps1) .or. (abs(mask_c(i-1,j,k)) < eps1)) then
           rhs(i,j,k) = rhs(i,j,k) - u0(i,j,k)*(var(i-1,j,k) + var(i,j,k))*dxi5
         end if

         if ((abs(mask_v(i,j+1,k)) < eps1) .or. (abs(mask_c(i,j+1,k)) < eps1)) then
           rhs(i,j,k) = rhs(i,j,k) + v0(i,j+1,k)*(var(i,j+1,k) + var(i,j,k))*dyi5
         end if

         if ((abs(mask_v(i,j,k)) < eps1) .or. (abs(mask_c(i,j-1,k)) < eps1)) then
           rhs(i,j,k) = rhs(i,j,k) - v0(i,j,k)*(var(i,j-1,k) + var(i,j,k))*dyi5
         end if

         if ((abs(mask_w(i,j,k+1)) < eps1) .or. (abs(mask_c(i,j,k+1)) < eps1)) then
           rhs(i,j,k) = rhs(i,j,k) + w0(i,j,k+1)*(var(i,j,k+1)*dzf(k) + var(i,j,k)*dzf(k+1))*dzhi(k+1)*dzfi5(k)
         end if

         if ((abs(mask_w(i,j,k)) < eps1) .or. (abs(mask_c(i,j,k-1)) < eps1)) then
           rhs(i,j,k) = rhs(i,j,k) - w0(i,j,k)*(var(i,j,k-1)*dzf(k) + var(i,j,k)*dzf(k-1))*dzhi(k)*dzfi5(k)
         end if

     end do

   end subroutine advecc2nd_corr_conservative


   subroutine advecc2nd_corr_liberal(var, rhs)
     ! Removes the advection contribution from solid scalar points as calculated
     ! by the 2nd order scheme, and replaces it with a contribution in which the
     ! value inside the solid is equal to the value outside, thereby modelling
     ! a zero (advective) flux condition.
     ! Due to potentially nonzero solid velocities due to the pressure correction,
     ! the IBM will not be conservative.
     use modglobal,      only : eps1, ib, ie, ih, jb, je, jh, kb, ke, kh, &
                                dx2i, dxi5, dy2i, dyi5, dzf, dzh2i, dzfi, dzhi, dzfi5
     use modfields,      only : u0, v0, w0
     use modsubgriddata, only : ekh
     use decomp_2d,      only : zstart

     real, intent(in)    :: var(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
     real, intent(inout) :: rhs(ib-ih:ie+ih,jb-jh:je+jh,kb   :ke+kh)
     integer :: i, j, k, n, m

     do m = 1,bound_info_c%nbndptsrank
      n = bound_info_c%bndptsrank(m)
         i = bound_info_c%bndpts(n,1) - zstart(1) + 1
         j = bound_info_c%bndpts(n,2) - zstart(2) + 1
         k = bound_info_c%bndpts(n,3) - zstart(3) + 1

         if (abs(mask_c(i+1,j,k)) < eps1) then ! var(i+1) is solid
           rhs(i,j,k) = rhs(i,j,k) + u0(i+1,j,k)*(var(i+1,j,k) + var(i,j,k))*dxi5 & ! negate contribution added in advection using var(i+1)
                                   - u0(i+1,j,k)*(var(i  ,j,k) + var(i,j,k))*dxi5   ! add corresponding contribution with var(i+1) = var(i)
         end if

         if (abs(mask_c(i-1,j,k)) < eps1) then ! var(i-1) is solid
           rhs(i,j,k) = rhs(i,j,k) - u0(i,j,k)*(var(i-1,j,k) + var(i,j,k))*dxi5 & ! negate contribution added in advection using var(i-1)
                                   + u0(i,j,k)*(var(i  ,j,k) + var(i,j,k))*dxi5   ! add corresponding contribution with var(i-1) = var(i)
         end if

         if (abs(mask_c(i,j+1,k)) < eps1) then ! var(j+1) is solid
           rhs(i,j,k) = rhs(i,j,k) + v0(i,j+1,k)*(var(i,j+1,k) + var(i,j,k))*dyi5 & ! negate contribution added in advection using var(j+1)
                                   - v0(i,j+1,k)*(var(i,j  ,k) + var(i,j,k))*dyi5   ! add corresponding contribution with var(j+1) = var(j)
         end if

         if (abs(mask_c(i,j-1,k)) < eps1) then ! var(j-1) is solid
           rhs(i,j,k) = rhs(i,j,k) - v0(i,j,k)*(var(i,j-1,k) + var(i,j,k))*dyi5 & ! negate contribution added in advection using var(j-1)
                                   + v0(i,j,k)*(var(i,j  ,k) + var(i,j,k))*dyi5   ! add corresponding contribution with var(j-1) = var(j)
         end if

         if (abs(mask_c(i,j,k+1)) < eps1) then ! var(k+1) is solid
           rhs(i,j,k) = rhs(i,j,k) + w0(i,j,k+1)*(var(i,j,k+1)*dzf(k) + var(i,j,k)*dzf(k+1))*dzhi(k+1)*dzfi5(k) & ! negate contribution added in advection using var(k+1)
                                   - w0(i,j,k+1)*(var(i,j,k  )*dzf(k) + var(i,j,k)*dzf(k+1))*dzhi(k+1)*dzfi5(k)   ! add corresponding contribution with var(k+1) = var(k)
         end if

         if (abs(mask_c(i,j,k-1)) < eps1) then ! var(k-1) is solid
           rhs(i,j,k) = rhs(i,j,k) - w0(i,j,k)*(var(i,j,k-1)*dzf(k) + var(i,j,k)*dzf(k-1))*dzhi(k)*dzfi5(k) & ! negate contribution added in advection using var(k-1)
                                   + w0(i,j,k)*(var(i,j,k  )*dzf(k) + var(i,j,k)*dzf(k-1))*dzhi(k)*dzfi5(k)   ! add corresponding contribution with var(k-1) = var(k)
         end if

     end do
   end subroutine advecc2nd_corr_liberal


   subroutine diffu_corr
     ! Negate subgrid rhs contributions from solid points (added by diffu in modsubgrid)
     use modglobal,      only : eps1, ib, ie, ih, jb, je, jh, kb, ke, kh, &
                                dx2i, dxi5, dy2i, dyi5, dzf, dzh2i, dzfi, dzhi, dzfi5, dzhiq
     use modfields,      only : u0, up
     use modsubgriddata, only : ekm
     use decomp_2d,      only : zstart

     real :: empo, emmo, emop, emom
     integer :: i, j, k, n, m

     do m = 1,bound_info_u%nbndptsrank
      n = bound_info_u%bndptsrank(m)
         i = bound_info_u%bndpts(n,1) - zstart(1) + 1
         j = bound_info_u%bndpts(n,2) - zstart(2) + 1
         k = bound_info_u%bndpts(n,3) - zstart(3) + 1

         if (abs(mask_u(i,j+1,k)) < eps1) then
           empo = 0.25 * ((ekm(i,j,k) + ekm(i,j+1,k)) + (ekm(i-1,j,k) + ekm(i-1,j+1,k)))
           up(i,j,k) = up(i,j,k) - empo * (u0(i,j+1,k) - u0(i,j,k))*dy2i
         end if

         if (abs(mask_u(i,j-1,k)) < eps1) then
           emmo = 0.25 * ((ekm(i,j,k) + ekm(i,j-1,k)) + (ekm(i-1,j-1,k) + ekm(i-1,j,k)))
           up(i,j,k) = up(i,j,k) + emmo * (u0(i,j,k) - u0(i,j-1,k))*dy2i
         end if

         if (abs(mask_u(i,j,k+1)) < eps1) then
           emop = (dzf(k+1) * ( ekm(i,j,k)   + ekm(i-1,j,k  ))  + &
                   dzf(k)   * ( ekm(i,j,k+1) + ekm(i-1,j,k+1))) * dzhiq(k+1)
           up(i,j,k) = up(i,j,k) - emop * (u0(i,j,k+1) - u0(i,j,k))*dzhi(k+1)*dzfi(k)
         end if

         if (abs(mask_u(i,j,k-1)) < eps1) then
           emom = (dzf(k-1) * (ekm(i,j,k  ) + ekm(i-1,j,k  ))  + &
                   dzf(k)   * (ekm(i,j,k-1) + ekm(i-1,j,k-1))) * dzhiq(k)
           up(i,j,k) = up(i,j,k) + emom * (u0(i,j,k) - u0(i,j,k-1))*dzhi(k)*dzfi(k)
         end if

     end do


   end subroutine diffu_corr


   subroutine diffv_corr
     ! Negate subgrid rhs contributions from solid points (added by diffv in modsubgrid)
     use modglobal,      only : eps1, ib, ie, ih, jb, je, jh, kb, ke, kh, &
                                dx2i, dxi5, dy2i, dyi5, dzf, dzh2i, dzfi, dzhi, dzfi5, dzhiq
     use modfields,      only : v0, vp
     use modsubgriddata, only : ekm
     use decomp_2d,      only : zstart

     real :: epmo, emmo, eomp, eomm
     integer :: i, j, k, n, m

     do m = 1,bound_info_v%nbndptsrank
      n = bound_info_v%bndptsrank(m)
         i = bound_info_v%bndpts(n,1) - zstart(1) + 1
         j = bound_info_v%bndpts(n,2) - zstart(2) + 1
         k = bound_info_v%bndpts(n,3) - zstart(3) + 1

         if (abs(mask_v(i+1,j,k)) < eps1) then
           epmo = 0.25 * (ekm(i,j,k) + ekm(i,j-1,k) + ekm(i+1,j-1,k) + ekm(i+1,j,k))
           vp(i,j,k) = vp(i,j,k) - epmo * (v0(i+1,j,k) - v0(i,j,k))*dx2i
         end if

         if (abs(mask_v(i-1,j,k)) < eps1) then
           emmo = 0.25 * (ekm(i,j,k) + ekm(i,j-1,k) + ekm(i-1,j-1,k) + ekm(i-1,j,k))
           vp(i,j,k) = vp(i,j,k) + emmo * (v0(i,j,k) - v0(i-1,j,k))*dx2i
         end if

         if (abs(mask_v(i,j,k+1)) < eps1) then
           eomp = ( dzf(k+1) * ( ekm(i,j,k)   + ekm(i,j-1,k)  )  + &
                    dzf(k  ) * ( ekm(i,j,k+1) + ekm(i,j-1,k+1))) * dzhiq(k+1)
           vp(i,j,k) = vp(i,j,k) - eomp * (v0(i,j,k+1) - v0(i,j,k))*dzhi(k+1)*dzfi(k)
         end if

         if (abs(mask_v(i,j,k-1)) < eps1) then
           eomm = ( dzf(k-1) * ( ekm(i,j,k  )  + ekm(i,j-1,k)   ) + &
                    dzf(k)   * ( ekm(i,j,k-1)  + ekm(i,j-1,k-1))) * dzhiq(k)
           vp(i,j,k) = vp(i,j,k) + eomm * (v0(i,j,k) - v0(i,j,k-1))*dzhi(k)*dzfi(k)
         end if

     end do

   end subroutine diffv_corr


   subroutine diffw_corr
     ! Negate subgrid rhs contributions from solid points (added by diffw in modsubgrid)
     use modglobal,      only : eps1, ib, ie, ih, jb, je, jh, kb, ke, kh, &
                                dx2i, dxi5, dy2i, dyi5, dzf, dzh2i, dzfi, dzhi, dzfi5, dzhiq
     use modfields,      only : w0, wp
     use modsubgriddata, only : ekm
     use decomp_2d,      only : zstart

     real :: epom, emom, eopm, eomm
     integer :: i, j, k, n, m

     do m = 1,bound_info_w%nbndptsrank
      n = bound_info_w%bndptsrank(m)
         i = bound_info_w%bndpts(n,1) - zstart(1) + 1
         j = bound_info_w%bndpts(n,2) - zstart(2) + 1
         k = bound_info_w%bndpts(n,3) - zstart(3) + 1

         ! Account for solid w points
         if (abs(mask_w(i+1,j,k)) < eps1) then
           epom = ( dzf(k-1) * ( ekm(i,j,k  ) + ekm(i+1,j,k  ))    + &
                    dzf(k  ) * ( ekm(i,j,k-1) + ekm(i+1,j,k-1))) * dzhiq(k)
           wp(i,j,k) = wp(i,j,k) - epom * (w0(i+1,j,k) - w0(i,j,k))*dx2i
         end if

         if (abs(mask_w(i-1,j,k)) < eps1) then
           emom = ( dzf(k-1) * ( ekm(i,j,k  ) + ekm(i-1,j,k  ))  + &
                    dzf(k  ) * ( ekm(i,j,k-1) + ekm(i-1,j,k-1))) * dzhiq(k)
           wp(i,j,k) = wp(i,j,k) + emom * (w0(i,j,k) - w0(i-1,j,k))*dx2i
         end if

         if (abs(mask_w(i,j+1,k)) < eps1) then
           eopm = ( dzf(k-1) * ( ekm(i,j,k  ) + ekm(i,j+1,k  ))  + &
                    dzf(k  ) * ( ekm(i,j,k-1) + ekm(i,j+1,k-1))) * dzhiq(k)
           wp(i,j,k) = wp(i,j,k) - eopm * (w0(i,j+1,k) - w0(i,j,k))*dy2i
         end if

         if (abs(mask_w(i,j-1,k)) < eps1) then
           eomm = ( dzf(k-1) * ( ekm(i,j,k  ) + ekm(i,j-1,k  ))  + &
                    dzf(k  ) * ( ekm(i,j,k-1) + ekm(i,j-1,k-1))) * dzhiq(k)
           wp(i,j,k) = wp(i,j,k) + eomm * (w0(i,j,k) - w0(i,j-1,k))*dy2i
         end if

     end do

   end subroutine diffw_corr


   subroutine diffc_corr(var, rhs, hi, hj, hk)
     ! Negate subgrid rhs contributions from solid points (added by diffc in modsubgrid)
     use modglobal,      only : eps1, ib, ie, jb, je, kb, ke, kh, &
                                dx2i, dxi5, dy2i, dyi5, dzf, dzh2i, dzfi, dzhi, dzfi5
     use modsubgriddata, only : ekh
     use decomp_2d,      only : zstart

     integer, intent(in) :: hi, hj, hk
     real, intent(in)    :: var(ib-hi:ie+hi,jb-hj:je+hj,kb-hk:ke+hk)
     real, intent(inout) :: rhs(ib-hi:ie+hi,jb-hj:je+hj,kb   :ke+hk)
     integer :: i, j, k, n, m

     do m = 1,bound_info_c%nbndptsrank
      n = bound_info_c%bndptsrank(m)
         i = bound_info_c%bndpts(n,1) - zstart(1) + 1
         j = bound_info_c%bndpts(n,2) - zstart(2) + 1
         k = bound_info_c%bndpts(n,3) - zstart(3) + 1

         if (abs(mask_c(i+1,j,k)) < eps1) then
           rhs(i,j,k) = rhs(i,j,k) - 0.5 * (ekh(i+1,j,k) + ekh(i,j,k)) * (var(i+1,j,k) - var(i,j,k))*dx2i
         end if

         if (abs(mask_c(i-1,j,k)) < eps1) then
           rhs(i,j,k) = rhs(i,j,k) + 0.5 * (ekh(i,j,k) + ekh(i-1,j,k)) * (var(i,j,k) - var(i-1,j,k))*dx2i
         end if

         if (abs(mask_c(i,j+1,k)) < eps1) then
           rhs(i,j,k) = rhs(i,j,k) - 0.5 * (ekh(i,j+1,k) + ekh(i,j,k)) * (var(i,j+1,k) - var(i,j,k))*dy2i
         end if

         if (abs(mask_c(i,j-1,k)) < eps1) then
           rhs(i,j,k) = rhs(i,j,k) + 0.5 * (ekh(i,j,k) + ekh(i,j-1,k)) * (var(i,j,k) - var(i,j-1,k))*dy2i
         end if

         if (abs(mask_c(i,j,k+1)) < eps1) then
           rhs(i,j,k) = rhs(i,j,k) - 0.5 * (dzf(k+1)*ekh(i,j,k) + dzf(k)*ekh(i,j,k+1)) &
                                         * (var(i,j,k+1) - var(i,j,k))*dzh2i(k+1)*dzfi(k)
         end if

         if (abs(mask_c(i,j,k-1)) < eps1) then
           rhs(i,j,k) = rhs(i,j,k) + 0.5 * (dzf(k-1)*ekh(i,j,k) + dzf(k)*ekh(i,j,k-1)) &
                                         * (var(i,j,k) - var(i,j,k-1))*dzh2i(k)*dzfi(k)
         end if

     end do

   end subroutine diffc_corr


   subroutine ibmwallfun
     use modglobal, only : libm, iwallmom, iwalltemp, xhat, yhat, zhat, ltempeq, lmoist, &
                           ib, ie, ih, ihc, jb, je, jh, jhc, kb, ke, kh, khc, nsv
     use modfields, only : u0, v0, w0, thl0, qt0, sv0, up, vp, wp, thlp, qtp, svp, &
                           tau_x, tau_y, tau_z, thl_flux
     use modsubgriddata, only : ekm, ekh

     real, allocatable :: rhs(:,:,:)
     integer n

      if (.not. libm) return

      allocate(rhs(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))

      if (iwallmom > 1) then
        rhs = up
        call wallfunmom(xhat, up, bound_info_u)
        tau_x(:,:,kb:ke+kh) = tau_x(:,:,kb:ke+kh) + (up - rhs)

        rhs = vp
        call wallfunmom(yhat, vp, bound_info_v)
        tau_y(:,:,kb:ke+kh) = tau_y(:,:,kb:ke+kh) + (vp - rhs)

        rhs = wp
        call wallfunmom(zhat, wp, bound_info_w)
        tau_z(:,:,kb:ke+kh) = tau_z(:,:,kb:ke+kh) + (wp - rhs)

        ! This replicates uDALES 1 behaviour, but probably should be done even if not using wall functions
        call diffu_corr
        call diffv_corr
        call diffw_corr
      end if

      if (ltempeq .or. lmoist) then
        rhs = thlp
        call wallfunheat
        thl_flux(:,:,kb:ke+kh) = thl_flux(:,:,kb:ke+kh) + (thlp - rhs)
        if (ltempeq) call diffc_corr(thl0, thlp, ih, jh, kh)
        if (lmoist)  call diffc_corr(qt0, qtp, ih, jh, kh)
      end if

      do n = 1,nsv
        call diffc_corr(sv0(:,:,:,n), svp(:,:,:,n), ihc, jhc, khc)
      end do

      deallocate(rhs)

    end subroutine ibmwallfun


   subroutine wallfunmom(dir, rhs, bound_info)
     use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh, xf, yf, zf, xh, yh, zh, &
                           eps1, fkar, dx, dy, dzf, iwallmom, xhat, yhat, zhat, vec0
     use modfields, only : u0, v0, w0, thl0, tau_x, tau_y, tau_z
     use initfac,   only : facT, facz0, facz0h, facnorm
     use decomp_2d, only : zstart

     real, intent(in)    :: dir(3)
     real, intent(inout) :: rhs(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
     type(bound_info_type) :: bound_info

     integer i, j, k, n, m, sec, pt, fac
     real dist, stress, stress_dir, stress_aligned, area, vol, momvol, Tair, Tsurf, x, y, z, &
          utan, udir, ctm, a, a_is, a_xn, a_yn, a_zn, stress_ix, stress_iy, stress_iz, xrec, yrec, zrec
     real, dimension(3) :: uvec, norm, strm, span, stressvec
     logical :: valid

     procedure(interp_velocity), pointer :: interp_velocity_ptr => null()
     procedure(interp_temperature), pointer :: interp_temperature_ptr => null()

     select case(alignment(dir))
     case(1)
       interp_velocity_ptr => interp_velocity_u
       interp_temperature_ptr => interp_temperature_u
     case(2)
       interp_velocity_ptr => interp_velocity_v
       interp_temperature_ptr => interp_temperature_v
     case(3)
       interp_velocity_ptr => interp_velocity_w
       interp_temperature_ptr => interp_temperature_w
     end select

     do m = 1,bound_info%nfctsecsrank
       sec = bound_info%fctsecsrank(m) ! index of section
       if (bound_info%lskipsec(sec)) cycle

       n = bound_info%secbndptids(sec) ! index of boundary point
       fac = bound_info%secfacids(sec) ! index of facet
       norm = facnorm(fac,:) ! facet normal

       i = bound_info%bndpts(n,1) - zstart(1) + 1
       j = bound_info%bndpts(n,2) - zstart(2) + 1
       k = bound_info%bndpts(n,3) - zstart(3) + 1
       if ((i < ib) .or. (i > ie) .or. (j < jb) .or. (j > je)) write(*,*) "problem", i, j

       if (bound_info%lcomprec(sec)) then
         uvec = interp_velocity_ptr(i, j, k)
         if (iwallmom == 2) then
           Tair = interp_temperature_ptr(i, j, k)
         end if
         dist = bound_info%bnddst(sec)
       else
         xrec = bound_info%recpts(sec,1)
         yrec = bound_info%recpts(sec,2)
         zrec = bound_info%recpts(sec,3)
         uvec(1) = trilinear_interp_var(u0, bound_info%recids_u(sec,:), xh, yf, zf, xrec, yrec, zrec)
         uvec(2) = trilinear_interp_var(v0, bound_info%recids_v(sec,:), xf, yh, zf, xrec, yrec, zrec)
         uvec(3) = trilinear_interp_var(w0, bound_info%recids_w(sec,:), xf, yf, zh, xrec, yrec, zrec)
         if (iwallmom == 2) Tair  = trilinear_interp_var(thl0, bound_info%recids_c(sec,:), xf, yf, zf, xrec, yrec, zrec)
         dist = bound_info%bnddst(sec) + norm2((/xrec - xf(bound_info%bndpts(n,1)), &
                                                 yrec - yf(bound_info%bndpts(n,2)), &
                                                 zrec - zf(bound_info%bndpts(n,3))/))
       end if

       if (is_equal(uvec, vec0)) cycle

       call local_coords(uvec, norm, span, strm, valid)
       if (.not. valid) cycle

       utan = dot_product(uvec, strm)
       !utan = max(0.01, utan) ! uDALES 1

       ! calcualate momentum transfer coefficient
       ! make into interface somehow? because iwallmom doesn't change in the loop
       if (iwallmom == 2) then ! stability included
         ctm = mom_transfer_coef_stability(utan, dist, facz0(fac), facz0h(fac), Tair, facT(fac,1))
       else if (iwallmom == 3) then ! neutral
         ctm = mom_transfer_coef_neutral(dist, facz0(fac))
       end if

       stress = ctm * utan**2

       if (bound_info%lcomprec(sec)) then
         a = dot_product(dir, strm)
         stress_dir = a * stress
       else
         ! Rotation from local (strm,span,norm) to global (xhat,yhat,zhat) basis
         ! \tau'_ij = a_ip a_jq \tau_pq
         ! \tau_pq in local coordinates is something like \tau \delta_13, because we only have \tau_{strm,norm})
         a_is = dot_product(dir, strm)
         a_xn = dot_product(xhat, norm)
         a_yn = dot_product(yhat, norm)
         a_zn = dot_product(zhat, norm)

         stress_ix = a_is * a_xn * stress
         stress_iy = a_is * a_yn * stress
         stress_iz = a_is * a_zn * stress

         stressvec(1) = stress_ix
         stressvec(2) = stress_iy
         stressvec(3) = stress_iz
         stress_dir = norm2(stressvec)
       end if

       stress_dir = sign(stress_dir, dot_product(uvec, dir))

       area = bound_info%secareas(sec)
       vol = dx*dy*dzf(k)
       momvol = stress_dir * area / vol
       rhs(i,j,k) = rhs(i,j,k) - momvol

     end do

   end subroutine wallfunmom


   subroutine wallfunheat
     use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh, xf, yf, zf, xh, yh, zh, dx, dy, dzh, eps1, &
                           xhat, yhat, zhat, vec0, fkar, ltempeq, lmoist, iwalltemp, iwallmoist, lEB
     use modfields, only : u0, v0, w0, thl0, thlp, qt0, qtp
     use initfac,   only : facT, facz0, facz0h, facnorm, faca, fachf, facef, facqsat, fachurel, facf, faclGR
     use modsurfdata, only : z0, z0h
     use modibmdata, only : bctfxm, bctfxp, bctfym, bctfyp, bctfz
     use decomp_2d, only : zstart

     integer i, j, k, n, m, sec, fac
     real :: dist, flux, area, vol, tempvol, Tair, Tsurf, utan, cth, cveg, hurel, qtair, qwall, resc, ress, xrec, yrec, zrec
     real, dimension(3) :: uvec, norm, span, strm
     logical :: valid

     do m = 1,bound_info_c%nfctsecsrank
       sec = bound_info_c%fctsecsrank(m) ! index of section
       n =   bound_info_c%secbndptids(sec) ! index of boundary point
       fac = bound_info_c%secfacids(sec) ! index of facet
       norm = facnorm(fac,:)

       i = bound_info_c%bndpts(n,1) - zstart(1) + 1 ! should be on this rank!
       j = bound_info_c%bndpts(n,2) - zstart(2) + 1 ! should be on this rank!
       k = bound_info_c%bndpts(n,3) - zstart(3) + 1 ! should be on this rank!
       if ((i < ib) .or. (i > ie) .or. (j < jb) .or. (j > je)) then
          write(*,*) "problem in wallfunheat", i, j
          stop 1
        end if

       if (bound_info_c%lcomprec(sec)) then ! section aligned with grid - use this cell's velocity
         uvec = interp_velocity_c(i, j, k)
         Tair = thl0(i,j,k)
         qtair = qt0(i,j,k)
         dist = bound_info_c%bnddst(sec)

       else ! use velocity at reconstruction point
         xrec = bound_info_c%recpts(sec,1)
         yrec = bound_info_c%recpts(sec,2)
         zrec = bound_info_c%recpts(sec,3)
         uvec(1) = trilinear_interp_var(u0, bound_info_c%recids_u(sec,:), xh, yf, zf, xrec, yrec, zrec)
         uvec(2) = trilinear_interp_var(v0, bound_info_c%recids_v(sec,:), xf, yh, zf, xrec, yrec, zrec)
         uvec(3) = trilinear_interp_var(w0, bound_info_c%recids_w(sec,:), xf, yf, zh, xrec, yrec, zrec)
         Tair  = trilinear_interp_var(thl0, bound_info_c%recids_c(sec,:), xf, yf, zf, xrec, yrec, zrec)
         qtair = trilinear_interp_var( qt0, bound_info_c%recids_c(sec,:), xf, yf, zf, xrec, yrec, zrec)
         dist = bound_info_c%bnddst(sec) + norm2((/xrec - xf(bound_info_c%bndpts(n,1)), &
                                                 yrec - yf(bound_info_c%bndpts(n,2)), &
                                                 zrec - zf(bound_info_c%bndpts(n,3))/))

       end if

       if (is_equal(uvec, vec0)) cycle

       call local_coords(uvec, norm, span, strm, valid)
       if (.not. valid) cycle
       utan = dot_product(uvec, strm)
       !utan = max(0.01, utan) ! uDALES 1

       ! Sensible heat
       if (ltempeq) then
         if (iwalltemp == 1) then ! probably remove this eventually, only relevant to grid-aligned facets
           !if     (all(abs(norm - xhat) < eps1)) then
           if     (is_equal(norm, xhat)) then
             flux = bctfxp
           !elseif (all(abs(norm + xhat) < eps1)) then
           elseif (is_equal(norm, -xhat)) then
             flux = bctfxm
           !elseif (all(abs(norm - yhat) < eps1)) then
           elseif (is_equal(norm, yhat)) then
             flux = bctfyp
           !elseif (all(abs(norm + yhat) < eps1)) then
           elseif (is_equal(norm, -yhat)) then
             flux = bctfxm
           !elseif (all(abs(norm - zhat) < eps1)) then
           elseif (is_equal(norm, zhat)) then
             flux = bctfz
           end if

         elseif (iwalltemp == 2) then
           call heat_transfer_coef_flux(utan, dist, facz0(fac), facz0h(fac), Tair, facT(fac, 1), cth, flux)
           ! Outputs 'cth' = heat transfer coefficient x utan = aerodynamic resistance as well as flux
         end if

         ! Heat transfer coefficient (cth) could be output here
         ! flux [Km/s]
         ! fluid volumetric sensible heat source/sink = flux * area / volume [K/s]
         ! facet sensible heat flux = volumetric heat capacity of air * flux * sectionarea / facetarea [W/m^2]
         thlp(i,j,k) = thlp(i,j,k) - flux * bound_info_c%secareas(sec) / (dx*dy*dzh(k))

         if (lEB) then
           fachf(fac) = fachf(fac) + flux * bound_info_c%secareas(sec) ! [Km^2/s] (will be divided by facetarea(fac) in modEB)
         end if
       end if

       ! Latent heat
       if (lmoist .and. faclGR(fac)) then
         if (iwallmoist == 1) then ! probably remove this eventually, only relevant to grid-aligned facets
           if     (is_equal(norm, xhat)) then
             flux = bcqfxp
           elseif (is_equal(norm, -xhat)) then
             flux = bcqfxm
           elseif (is_equal(norm, yhat)) then
             flux = bcqfyp
           elseif (is_equal(norm, -yhat)) then
             flux = bcqfym
           elseif (is_equal(norm, zhat)) then
             flux = bcqfz
           end if

         elseif (iwallmoist == 2) then
           qwall = facqsat(fac)
           hurel = fachurel(fac)
           resc = facf(fac,4)
           ress = facf(fac,5)
           cveg = 0.8
           ! flux = min(0., cveg * (qtair - qwall)         / (1/cth + resc) + &
           !           (1 - cveg)* (qtair - qwall * hurel) / (1/cth + ress))
           flux = moist_flux(cveg, cth, qtair, qwall, hurel, resc, ress)
         end if

         ! flux [kg/kg m/s]
         ! fluid volumetric latent heat source/sink = flux * area / volume [kg/kg / s]
         ! facet latent heat flux = volumetric heat capacity of air * flux * sectionarea / facetarea [W/m^2]
         qtp(i,j,k) = qtp(i,j,k) - flux * bound_info_c%secareas(sec) / (dx*dy*dzh(k))

         if (lEB) then
           facef(fac) = facef(fac) + flux * bound_info_c%secareas(sec) ! [Km^2/s] (will be divided by facetarea(fac) in modEB)
         end if
       end if

     end do

   end subroutine wallfunheat


   real function trilinear_interp_var(var, cell, xgrid, ygrid, zgrid, x, y, z)
     use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh, itot, jtot, ktot
     use decomp_2d, only : zstart
     implicit none
     real, intent(in)    :: var(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:kb+kh)
     integer, intent(in) :: cell(3) ! GLOBAL indices of cell containing the point
     real, intent(in), dimension(ib:itot+ih) :: xgrid
     real, intent(in), dimension(jb:jtot+jh) :: ygrid
     real, intent(in), dimension(kb:ktot+kh) :: zgrid
     real,    intent(in) :: x, y, z ! location of point to interpolate at
     real, dimension(8)  :: corners(8)
     real :: x0, y0, z0, x1, y1, z1
     integer :: i, j, k

     i = cell(1) - zstart(1) + 1
     j = cell(2) - zstart(2) + 1
     k = cell(3) - zstart(3) + 1
     if ((i < ib-1) .or. (i > ie+1) .or. (j < jb-1) .or. (j > je+1)) then
       write(*,*) "problem in trilinear_interp_var", i, j, k
       stop 1
     end if
     corners = eval_corners(var, i, j, k)

     x0 = xgrid(cell(1))
     y0 = ygrid(cell(2))
     z0 = zgrid(cell(3))
     x1 = xgrid(cell(1)+1)
     y1 = ygrid(cell(2)+1)
     z1 = zgrid(cell(3)+1)

     trilinear_interp_var = trilinear_interp(x, y, z, x0, y0, z0, x1, y1, z1, corners)

   end function trilinear_interp_var

   function eval_corners(var, i, j, k)
     use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh
     integer, intent(in) :: i, j, k ! LOCAL indices
     real, intent(in)    :: var(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:kb+kh)
     real, dimension(8)  :: eval_corners(8)

     eval_corners(1) = var(i  ,j  ,k  ) !c000
     eval_corners(2) = var(i+1,j  ,k  ) !c100
     eval_corners(3) = var(i  ,j+1,k  ) !c010
     eval_corners(4) = var(i+1,j+1,k  ) !c110
     eval_corners(5) = var(i  ,j  ,k+1) !c001
     eval_corners(6) = var(i+1,j  ,k+1) !c101
     eval_corners(7) = var(i  ,j+1,k+1) !c011
     eval_corners(8) = var(i+1,j+1,k+1) !c111

   end function eval_corners


   real function trilinear_interp(x, y, z, x0, y0, z0, x1, y1, z1, corners)
     real, intent(in) :: x, y, z, x0, y0, z0, x1, y1, z1, corners(8)
     real :: xd, yd, zd

     xd = (x - x0) / (x1 - x0)
     yd = (y - y0) / (y1 - y0)
     zd = (z - z0) / (z1 - z0)
     ! check all positive

     trilinear_interp = corners(1) * (1-xd)*(1-yd)*(1-zd) + & ! c000
                        corners(2) * (  xd)*(1-yd)*(1-zd) + & ! c100
                        corners(3) * (1-xd)*(  yd)*(1-zd) + & ! c010
                        corners(4) * (  xd)*(  yd)*(1-zd) + & ! c110
                        corners(5) * (1-xd)*(1-yd)*(  zd) + & ! c001
                        corners(6) * (  xd)*(1-yd)*(  zd) + & ! c101
                        corners(7) * (1-xd)*(  yd)*(  zd) + & ! c011
                        corners(8) * (  xd)*(  yd)*(  zd)     ! c111

   end function trilinear_interp


   integer function alignment(n)
     ! returns an integer determining whether a unit vector n is aligned with the
     ! coordinates axes.
     use modglobal, only : xhat, yhat, zhat
     implicit none
     real, dimension(3), intent(in) :: n ! must be unit vector

     if     (is_equal(n, xhat)) then
       alignment = 1
     elseif (is_equal(n, yhat)) then
       alignment = 2
     elseif (is_equal(n, zhat)) then
       alignment = 3
     elseif (is_equal(n, -xhat)) then
       alignment = -1
     elseif (is_equal(n, -yhat)) then
       alignment = -2
     elseif (is_equal(n, -zhat)) then
       alignment = -3
     else
       alignment = 0
     end if

   end function alignment

   ! overload this with real dimension(1)
   ! could put somewhere else because not specific to ibm
   logical function is_equal(a,b)
     ! determines whether two vectors are equal to each other within a tolerance of eps1
     use modglobal, only : eps1
     implicit none
     real, dimension(3), intent(in) :: a, b

     if (all(abs(a - b) < eps1)) then
       is_equal = .true.
     else
       is_equal = .false.
     end if

   end function is_equal


   function cross_product(a,b)
     ! Calculate the cross product (a x b)
     implicit none
     real, dimension(3) :: cross_product
     real, dimension(3), intent(in) :: a, b

     cross_product(1) = a(2)*b(3) - a(3)*b(2)
     cross_product(2) = a(3)*b(1) - a(1)*b(3)
     cross_product(3) = a(1)*b(2) - a(2)*b(1)

   end function cross_product


   function interp_velocity_u(i, j, k)
     ! interpolates the velocity at u-grid location i,j,k
     use modfields, only :  u0, v0, w0
     real ::  interp_velocity_u(3)
     integer, intent(in) :: i, j, k

     interp_velocity_u(1) = u0(i,j,k)
     interp_velocity_u(2) = 0.25 * (v0(i,j,k) + v0(i,j+1,k) + v0(i-1,j,k) + v0(i-1,j+1,k))
     interp_velocity_u(3) = 0.25 * (w0(i,j,k) + w0(i,j,k+1) + w0(i-1,j,k) + w0(i-1,j,k+1)) !only for equidistant grid!

     return
   end function interp_velocity_u


   function interp_velocity_v(i, j, k)
     ! interpolates the velocity at v-grid location i,j,k
     use modfields, only :  u0, v0, w0
     real ::  interp_velocity_v(3)
     integer, intent(in) :: i, j, k

     interp_velocity_v(1) = 0.25 * (u0(i,j,k) + u0(i+1,j,k) + u0(i,j-1,k) + u0(i+1,j-1,k))
     interp_velocity_v(2) = v0(i,j,k)
     interp_velocity_v(3) = 0.25 * (w0(i,j,k) + w0(i,j,k+1) + w0(i,j-1,k) + w0(i,j-1,k+1)) !only for equidistant grid!

     return
   end function interp_velocity_v


   function interp_velocity_w(i, j, k)
     ! interpolates the velocity at w-grid location i,j,k
     use modfields, only :  u0, v0, w0
     real ::  interp_velocity_w(3)
     integer, intent(in) :: i, j, k

     interp_velocity_w(1) = 0.25 * (u0(i,j,k) + u0(i+1,j,k) + u0(i,j-1,k) + u0(i+1,j-1,k))
     interp_velocity_w(2) = v0(i,j,k)
     interp_velocity_w(3) = 0.25 * (w0(i,j,k) + w0(i,j,k+1) + w0(i,j-1,k) + w0(i,j-1,k+1)) !only for equidistant grid!

     return
   end function interp_velocity_w


   function interp_velocity_c(i, j, k)
     ! interpolates the velocity at c-grid location i,j,k
     use modfields, only :  u0, v0, w0
     real ::  interp_velocity_c(3)
     integer, intent(in) :: i, j, k

     interp_velocity_c(1) = 0.5 * (u0(i,j,k) + u0(i+1,j,k))
     interp_velocity_c(2) = 0.5 * (v0(i,j,k) + v0(i,j+1,k))
     interp_velocity_c(3) = 0.5 * (w0(i,j,k) + w0(i,j,k+1))

     return
   end function interp_velocity_c


   real function interp_temperature_u(i, j, k)
     ! interpolates the temperature at u-grid location i,j,k
     use modfields, only :  thl0
     integer, intent(in) :: i, j, k

     !interp_temperature_u = 0.5 * (thl0(i,j,k) + thl0(i-1,j,k))
     interp_temperature_u = 0.5 * (thl0(i  ,j,k)*mask_c(i  ,j,k)*(2.-mask_c(i-1,j,k)) &
                                +  thl0(i-1,j,k)*mask_c(i-1,j,k)*(2.-mask_c(i  ,j,k)))

     return
   end function interp_temperature_u


   real function interp_temperature_v(i, j, k)
     ! interpolates the temperature at v-grid location i,j,k
     use modfields, only :  thl0
     integer, intent(in) :: i, j, k

     !interp_temperature_v = 0.5 * (thl0(i,j,k) + thl0(i,j-1,k))
     interp_temperature_v = 0.5 * (thl0(i,j  ,k)*mask_c(i,j  ,k)*(2.-mask_c(i,j-1,k)) &
                                 + thl0(i,j-1,k)*mask_c(i,j-1,k)*(2.-mask_c(i,j  ,k)))

     return
   end function interp_temperature_v


   real function interp_temperature_w(i, j, k)
     ! interpolates the temperature at w-grid location i,j,k
     use modfields, only :  thl0
     integer, intent(in) :: i, j, k

     !interp_temperature_w = 0.5 * (thl0(i,j,k) + thl0(i,j,k-1))
     interp_temperature_w = 0.5 * (thl0(i,j,k  )*mask_c(i,j,k  )*(2.-mask_c(i,j,k-1)) &
                                +  thl0(i,j,k-1)*mask_c(i,j,k-1)*(2.-mask_c(i,j,k  )))

     return
   end function interp_temperature_w


   subroutine local_coords(uvec, norm, span, strm, valid)
     ! returns the local streamwise (strm) and spanwise vectors (span) in the
     ! plane normal to (norm) containing the velocity vector (uvec)
     use modglobal, only : vec0
     real, intent(in),  dimension(3) :: uvec, norm
     real, intent(out), dimension(3) :: span, strm
     logical, intent(out) :: valid

     span = cross_product(norm, uvec)
     !if (is_equal(span, (/0.,0.,0./))) then
     ! velocity is pointing into or outof the surface
     if (is_equal(span, vec0)) then
       strm = 0.
       valid = .false.
     else
       span = span / norm2(span)
       valid = .true.
     end if
     strm = cross_product(span, norm)

   end subroutine local_coords


   real function mom_transfer_coef_stability(utan, dist, z0, z0h, Tair, Tsurf)
     ! By Ivo Suter. calculates the momentum transfer coefficient based on the
     ! surface tangential velocity 'utan' at a distance 'dist' from the surface,
     ! for a surface with momentum roughness length z0 and heat roughness length z0h.
     ! Stability are included using the air temperature Tair and surface temperature Tsurf.
     use modglobal, only : grav, fkar, prandtlturb

      implicit none
      real, intent(in) :: dist, z0, z0h, Tsurf, Tair, utan
      real, parameter :: b1 = 9.4 !parameters from uno1995
      real, parameter :: b2 = 4.7
      real, parameter :: dm = 7.4
      real, parameter :: dh = 5.3
      real :: dT, Ribl0, logdz, logdzh, logzh, sqdz, fkar2, Ribl1, Fm, Fh, cm, ch, Ctm, M

      dT = Tair - Tsurf
      Ribl0 = grav * dist * dT / (Tsurf * utan**2) !Eq. 6, guess initial Ri

      logdz = LOG(dist/z0)
      logdzh = LOG(dist/z0h)
      logzh = LOG(z0/z0h)
      sqdz = SQRT(dist/z0)
      fkar2 = fkar**2

      IF (Ribl0 > 0.) THEN !0.25 approx critical for bulk Richardson number  => stable
         Fm = 1./(1. + b2*Ribl0)**2 !Eq. 4
         Fh = Fm !Eq. 4
      ELSE ! => unstable
         cm = (dm*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
         ch = (dh*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
         Fm = 1. - (b1*Ribl0)/(1. + cm*SQRT(ABS(Ribl0))) !Eq. 3
         Fh = 1. - (b1*Ribl0)/(1. + ch*SQRT(ABS(Ribl0))) !Eq. 3
      END IF

      M = prandtlturb*logdz*SQRT(Fm)/Fh !Eq. 14

      Ribl1 = Ribl0 - Ribl0*prandtlturb*logzh/(prandtlturb*logzh + M) !Eq. 17

      !interate to get new Richardson number
      IF (Ribl1 > 0.) THEN !0.25 approx critical for bulk Richardson number  => stable
         Fm = 1./(1. + b2*Ribl1)**2 !Eq. 4
      ELSE ! => unstable
         cm = (dm*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
         Fm = 1. - (b1*Ribl1)/(1. + cm*SQRT(ABS(Ribl1))) !Eq. 3
      END IF

      mom_transfer_coef_stability = fkar2/(logdz**2)*Fm !Eq. 7

   end function mom_transfer_coef_stability


   real function mom_transfer_coef_neutral(dist, z0)
     ! calculates the heat transfer coefficient based on the (neutral) log law,
     ! for a distance 'dist' and a momentum roughness length 'z0'.
     use modglobal, only : fkar

     implicit none
     real, intent(in) :: dist, z0

     mom_transfer_coef_neutral = (fkar / log(dist / z0))**2

   end function mom_transfer_coef_neutral


   subroutine heat_transfer_coef_flux(utan, dist, z0, z0h, Tair, Tsurf, cth, flux)
     use modglobal, only : grav, fkar, prandtlturb

      implicit none
      real, intent(in)  :: dist, z0, z0h, Tsurf, Tair, utan
      real, intent(out) :: cth, flux
      real, parameter :: b1 = 9.4 !parameters from Uno1995
      real, parameter :: b2 = 4.7
      real, parameter :: dm = 7.4
      real, parameter :: dh = 5.3
      !real :: Pr
      real :: dT, Ribl0, logdz, logdzh, logzh, sqdz, fkar2, Ribl1, Fm, Fh, cm, ch, M, dTrough

      !Pr = 1.
      !Pr = prandtlmol
      dT = Tair - Tsurf
      Ribl0 = grav * dist * dT / (Tsurf * utan**2) !Eq. 6, guess initial Ri

      logdz = log(dist/z0)
      logdzh = log(dist/z0h)
      logzh = log(z0/z0h)
      sqdz = sqrt(dist/z0)
      fkar2 = fkar**2

      cth = 0.
      flux = 0.
      if (Ribl0 > 0.) then
         Fm = 1./(1. + b2*Ribl0)**2 !Eq. 4
         Fh = Fm !Eq. 4
      else ! => unstable
         cm = (dm*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
         ch = (dh*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
         Fm = 1. - (b1*Ribl0)/(1. + cm*sqrt(abs(Ribl0))) !Eq. 3
         Fh = 1. - (b1*Ribl0)/(1. + ch*sqrt(abs(Ribl0))) !Eq. 3
      end if

      M = prandtlturb*logdz*sqrt(Fm)/Fh !Eq. 14
      Ribl1 = Ribl0 - Ribl0*prandtlturb*logzh/(prandtlturb*logzh + M) !Eq. 17

      !interate to get new Richardson number
      if (Ribl1 > 0.) then
         Fm = 1./(1. + b2*Ribl1)**2 !Eq. 4
         Fh = Fm !Eq. 4
      else ! => unstable
         cm = (dm*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
         ch = (dh*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
         Fm = 1. - (b1*Ribl1)/(1. + cm*sqrt(abs(Ribl1))) !Eq. 3
         Fh = 1. - (b1*Ribl1)/(1. + ch*sqrt(abs(Ribl1))) !Eq. 3
      end if

      ! Uno (2)
      M = prandtlturb*logdz*sqrt(Fm)/Fh !Eq. 14
      dTrough = dT*1./(prandtlturb*logzh/M + 1.) !Eq. 13a
      cth = abs(utan)*fkar2/(logdz*logdz)*Fh/prandtlturb !Eq. 8
      flux = cth*dTrough !Eq. 2, Eq. 8

      ! ! Uno (8)
      ! cth = abs(utan)*fkar2/(logdz*logdzh)*Fh/prandtlturb !Eq. 8
      ! flux = cth*dT !Eq. 2, Eq. 8

   end subroutine heat_transfer_coef_flux


   real function moist_flux(cveg, cth, qtair, qwall, hurel, resc, ress)
     real, intent(in) :: cveg, cth, qtair, qwall, hurel, resc, ress

     moist_flux = min(0., cveg * (qtair - qwall)         / (1/cth + resc) + &
                     (1 - cveg)* (qtair - qwall * hurel) / (1/cth + ress))

   end function moist_flux


   subroutine bottom
     ! By Ivo Suter.
      !kind of obsolete when road facets are being used
      !vegetated floor not added (could simply be copied from vegetated horizontal facets)
      use modglobal, only:ib, ie, ih, jh, kb,ke,kh, jb, je, kb, numol, prandtlmol, dzh, nsv, &
         dxf, dxhi, dzf, dzfi, numoli, ltempeq, khc, lmoist, BCbotT, BCbotq, BCbotm, BCbots, dzh2i, libm
      use modfields, only : u0,v0,e120,um,vm,w0,wm,e12m,thl0,qt0,sv0,thlm,qtm,svm,up,vp,wp,thlp,qtp,svp,shear,momfluxb,tfluxb,cth,tau_x,tau_y,tau_z,thl_flux
      use modsurfdata, only:thlflux, qtflux, svflux, ustar, thvs, wtsurf, wqsurf, thls, z0, z0h
      use modsubgriddata, only:ekm, ekh
      use modmpi, only:myid
      implicit none
      integer :: i, j, jp, jm, m

      e120(:, :, kb - 1) = e120(:, :, kb)
      e12m(:, :, kb - 1) = e12m(:, :, kb)
      ! wm(:, :, kb) = 0. ! SO moved to modboundary
      ! w0(:, :, kb) = 0.
      tau_x(:,:,kb:ke+kh) = up
      tau_y(:,:,kb:ke+kh) = vp
      tau_z(:,:,kb:ke+kh) = wp
      thl_flux(:,:,kb:ke+kh) = thlp

      !if (.not.(libm)) then
      if (lbottom) then
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

      end if

      tau_x(:,:,kb:ke+kh) = up - tau_x(:,:,kb:ke+kh)
      tau_y(:,:,kb:ke+kh) = vp - tau_y(:,:,kb:ke+kh)
      tau_z(:,:,kb:ke+kh) = wp - tau_z(:,:,kb:ke+kh)
      thl_flux(:,:,kb:ke+kh) = thlp - thl_flux(:,:,kb:ke+kh)

      return
   end subroutine bottom


   subroutine createmasks
      use modglobal, only : libm, ib, ie, ih, ihc, jb, je, jh, jhc, kb, ke, kh, khc, itot, jtot, rslabs
      use modfields, only : IIc,  IIu,  IIv,  IIw,  IIuw,  IIvw,  IIuv,  &
                            IIcs, IIus, IIvs, IIws, IIuws, IIvws, IIuvs, &
                            IIct, IIut, IIvt, IIwt, IIuwt, um, u0, vm, v0, wm, w0
      use modmpi,    only : myid, comm3d, mpierr, MPI_INTEGER, MPI_DOUBLE_PRECISION, MY_REAL, nprocs, MPI_SUM
      use decomp_2d, only : zstart, exchange_halo_z

      integer :: IIcl(kb:ke + khc), IIul(kb:ke + khc), IIvl(kb:ke + khc), IIwl(kb:ke + khc), IIuwl(kb:ke + khc), IIvwl(kb:ke + khc), IIuvl(kb:ke + khc)
      integer :: IIcd(ib:ie, kb:ke)
      integer :: IIwd(ib:ie, kb:ke)
      integer :: IIuwd(ib:ie, kb:ke)
      integer :: IIud(ib:ie, kb:ke)
      integer :: IIvd(ib:ie, kb:ke)
      integer :: i, j, k, n, m

      ! II*l needn't be defined up to ke_khc, but for now would require large scale changes in modstatsdump so if works leave as is ! tg3315 04/07/18

      if (.not. libm) then
         IIc(:, :, :) = 1
         IIu(:, :, :) = 1
         IIv(:, :, :) = 1
         IIw(:, :, :) = 1
         IIuw(:, :, :) = 1
         IIvw(:, :, :) = 1
         IIuv(:, :, :) = 1
         IIcs(:) = nint(rslabs)
         IIus(:) = nint(rslabs)
         IIvs(:) = nint(rslabs)
         IIws(:) = nint(rslabs)
         IIuws(:) = nint(rslabs)
         IIvws(:) = nint(rslabs)
         IIuvs(:) = nint(rslabs)
         IIct(:, :) = jtot
         IIut(:, :) = jtot
         IIvt(:, :) = jtot
         IIwt(:, :) = jtot
         IIuwt(:, :) = jtot
         return
      end if
      ! Create masking matrices
      IIc = 1; IIu = 1; IIv = 1; IIct = 1; IIw = 1; IIuw = 1; IIvw = 1; IIuv = 1; IIwt = 1; IIut = 1; IIvt = 1; IIuwt = 1; IIcs = 1; IIus = 1; IIvs = 1; IIws = 1; IIuws = 1; IIvws = 1; IIuvs = 1

      do m = 1,solid_info_u%nsolptsrank
       n = solid_info_u%solptsrank(m)
          i = solid_info_u%solpts(n,1) - zstart(1) + 1
          j = solid_info_u%solpts(n,2) - zstart(2) + 1
          k = solid_info_u%solpts(n,3) - zstart(3) + 1
          IIu(i,j,k) = 0
      end do

      do m = 1,solid_info_v%nsolptsrank
       n = solid_info_v%solptsrank(m)
          i = solid_info_v%solpts(n,1) - zstart(1) + 1
          j = solid_info_v%solpts(n,2) - zstart(2) + 1
          k = solid_info_v%solpts(n,3) - zstart(3) + 1
          IIv(i,j,k) = 0
      end do

      do m = 1,solid_info_w%nsolptsrank
       n = solid_info_w%solptsrank(m)
          i = solid_info_w%solpts(n,1) - zstart(1) + 1
          j = solid_info_w%solpts(n,2) - zstart(2) + 1
          k = solid_info_w%solpts(n,3) - zstart(3) + 1
          IIw(i,j,k) = 0
      end do

      do m = 1,solid_info_c%nsolptsrank
       n = solid_info_c%solptsrank(m)
          i = solid_info_c%solpts(n,1) - zstart(1) + 1
          j = solid_info_c%solpts(n,2) - zstart(2) + 1
          k = solid_info_c%solpts(n,3) - zstart(3) + 1
          IIc(i,j,k) = 0
      end do

      IIw(:, :, kb) = 0; IIuw(:, :, kb) = 0; IIvw(:, :, kb) = 0

      do i=ib,ie
        do j=jb,je
          IIuv(i,j,kb) = IIu(i,j,kb) * IIu(i,j-1,kb) * IIv(i,j,kb) * IIv(i-1,j,kb)
          do k=kb+1,ke
            ! Classed as solid (set to zero) unless ALL points in the stencil are fluid
            IIuv(i,j,k) = IIu(i,j,k) * IIu(i,j-1,k) * IIv(i,j,k) * IIv(i-1,j,k)
            IIuw(i,j,k) = IIu(i,j,k) * IIu(i,j,k-1) * IIw(i,j,k) * IIw(i-1,j,k)
            IIvw(i,j,k) = IIv(i,j,k) * IIv(i,j,k-1) * IIw(i,j,k) * IIw(i,j-1,k)
          end do
        end do
      end do

      ! Can't do this because no interface for integers
      ! call exchange_halo_z(IIuv, opt_zlevel=(/ihc,jhc,0/))
      ! call exchange_halo_z(IIuv, opt_zlevel=(/ihc,jhc,0/))
      ! call exchange_halo_z(IIvw, opt_zlevel=(/ihc,jhc,0/))

      do k = kb, ke + khc
         IIcl(k) = sum(IIc(ib:ie, jb:je, k))
         IIul(k) = sum(IIu(ib:ie, jb:je, k))
         IIvl(k) = sum(IIv(ib:ie, jb:je, k))
         IIwl(k) = sum(IIw(ib:ie, jb:je, k))
         IIuwl(k) = sum(IIuw(ib:ie, jb:je, k))
         IIvwl(k) = sum(IIvw(ib:ie, jb:je, k))
         IIuvl(k) = sum(IIuv(ib:ie, jb:je, k))
      enddo

      call MPI_ALLREDUCE(IIcl, IIcs, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIul, IIus, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIvl, IIvs, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIwl, IIws, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIuwl, IIuws, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIvwl, IIvws, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIuvl, IIuvs, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)

      IIcd(ib:ie, kb:ke) = sum(IIc(ib:ie, jb:je, kb:ke), DIM=2)
      IIwd(ib:ie, kb:ke) = sum(IIw(ib:ie, jb:je, kb:ke), DIM=2)
      IIuwd(ib:ie, kb:ke) = sum(IIuw(ib:ie, jb:je, kb:ke), DIM=2)
      IIud(ib:ie, kb:ke) = sum(IIu(ib:ie, jb:je, kb:ke), DIM=2)
      IIvd(ib:ie, kb:ke) = sum(IIv(ib:ie, jb:je, kb:ke), DIM=2)

      call MPI_ALLREDUCE(IIwd(ib:ie, kb:ke), IIwt(ib:ie, kb:ke), (ke - kb + 1)*(ie - ib + 1), MPI_INTEGER, MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIcd(ib:ie, kb:ke), IIct(ib:ie, kb:ke), (ke - kb + 1)*(ie - ib + 1), MPI_INTEGER, MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIuwd(ib:ie, kb:ke), IIuwt(ib:ie, kb:ke), (ke - kb + 1)*(ie - ib + 1), MPI_INTEGER, MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIud(ib:ie, kb:ke), IIut(ib:ie, kb:ke), (ke - kb + 1)*(ie - ib + 1), MPI_INTEGER, MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIvd(ib:ie, kb:ke), IIvt(ib:ie, kb:ke), (ke - kb + 1)*(ie - ib + 1), MPI_INTEGER, MPI_SUM, comm3d, mpierr)

   end subroutine createmasks

end module modibm
