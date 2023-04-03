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
   public :: createwalls, ibmwallfun, xwallfun, ywallfunplus, ywallfunmin, &
             zwallfun, ibmnorm, nearwall, bottom

   contains
   subroutine createwalls
      use modglobal, only:ib, ie, jb, je, jgb, jge, kb, ke, jmax, dzf, zh, nblocks, &
         nsv, cexpnr, ifinput, libm, ih, kh, iwallmom, iwalltemp, iwallmoist, rslabs, bldT
      use modsurfdata, only:thls, z0h, z0, thvs
      use modfields, only:sv0, svm, thl0, thlm, qtp, qt0, IIc, IIu, IIv, IIw, IIct, IIwt, IIcs, IIus, IIvs, IIws, thlprof
      use modmpi, only:myid, comm3d, mpierr, MPI_INTEGER, MPI_DOUBLE_PRECISION, MY_REAL, nprocs, cmyid, &
         MPI_REAL8, MPI_REAL4, MPI_SUM
      use initfac, only:block
      integer n, nn, pn, mn, jbeg, jend, nxn, nxs, nyn, nzn, nzs, iu, il, ju, jl, ku, kl, sc, &
         i, k, dbi, dbj, dbk
      integer :: IIcl(kb:ke + kh), IIul(kb:ke + kh), IIvl(kb:ke + kh), IIwl(kb:ke + kh)
      integer :: IIcd(ib - ih:ie + ih, kb:ke + kh)
      integer :: IIwd(ib - ih:ie + ih, kb:ke + kh)
      real    :: thlprofvolav
      character(80) chmess, name2

      if (.not. libm) return

      nxwall = 0
      do n = 1, nblocks ! first x and z walls
         if (block(n, 4) < jb + myid*jmax) then ! no x-wall/z-wall in the range of this proc
            cycle
         elseif (block(n, 3) > je + myid*jmax) then ! no x-wall/z-wall in the range of this proc
            cycle
         else ! x-wall/z-wall found on this proc
            nxwall = nxwall + 1
         end if
      end do
      allocate (ixwall(nxwall)) !allocate the list that stores the indeces of the blocks on this cpu

      k = 1
      do n = 1, nblocks ! save indeces of the found x/z-walls by re-iterating
         if (block(n, 4) < jb + myid*jmax) then ! no x-wall/z-wall in the range of this proc
            cycle
         elseif (block(n, 3) > je + myid*jmax) then ! no x-wall/z-wall in the range of this proc
            cycle
         else
            ixwall(k) = n ! save index of block which is on this processor
            k = k + 1
         end if
      end do

      !!new approach both y walls#################################################
      !!store index of block and index of the wall (since block might not be on this cpu, but is needed for x and z coords)
      !!check if wall is on next cpu but not on this
      !!check if wall is on last cpu but not on first (periodicity in y)
      !check if wall is on first cpu but not on last (periodicity in y)
      !!check if wall is on this cpu and another one on the next (i.e. both blocks end at cpu boundary, but touch each other)

      nyminwall = 0
      nypluswall = 0

      do n = 1, nblocks

         jl = block(n, 3) - myid*jmax
         ju = block(n, 4) - myid*jmax

         !IMPORTANT: THESE LINES OF CODE SHOULD BE HERE BUT CAUSE TROUBLE! MAKE SURE BLOCK DOES NOT TOUCH BOUNDARY (EXCECPT ALL FLOORS)
         !SEE ALSO BELOW! (Like 40 lines or so)
         if ((myid == 0) .and. (block(n, 4) == jge)) then ! periodicity!
            nypluswall = nypluswall + 1
         else if ((block(n, 3) == jgb) .and. (myid == (nprocs - 1))) then ! periodicity!
            nyminwall = nyminwall + 1
         end if

         if ((ju < (jb - 1)) .or. (jl > (je + 1))) then
            cycle
         end if

         if (ju == (jb - 1)) then !block on previous cpu, north wall on this
            nypluswall = nypluswall + 1 !
            cycle
         end if

         if (jl == (je + 1)) then
            nyminwall = nyminwall + 1 !block on next cpu, southwall on this
            cycle
         end if

         if ((ju < je) .and. (ju >= jb)) then !block & northwall on this cpu
            nypluswall = nypluswall + 1
         end if

         if ((jl > jb) .and. (jl <= je)) then !block & southwall on this cpu
            nyminwall = nyminwall + 1
         end if
      end do

      allocate (iyminwall(1:nyminwall, 1:2)) !two indeces to store wall index and block index
      allocate (iypluswall(1:nypluswall, 1:2))
      iyminwall(:, 1) = 0
      iyminwall(:, 2) = 0
      iypluswall(:, 1) = 0
      iypluswall(:, 2) = 0
      pn = 1
      mn = 1

      do n = 1, nblocks

         jl = block(n, 3) - myid*jmax
         ju = block(n, 4) - myid*jmax

         !IMPORTANT: THESE LINES OF CODE SHOULD BE HERE BUT CAUSE TROUBLE! MAKE SURE BLOCK DOES NOT TOUCH BOUNDARY (EXCECPT ALL FLOORS)
         if ((myid == 0) .and. (block(n, 4) == jge)) then ! periodicity!
            iypluswall(pn, 1) = n
            iypluswall(pn, 2) = jb
            pn = pn + 1
         else if ((block(n, 3) == jgb) .and. (myid == (nprocs - 1))) then ! periodicity!
            iyminwall(mn, 1) = n
            iyminwall(mn, 2) = je
            mn = mn + 1
         end if

         if ((ju < (jb - 1)) .or. (jl > (je + 1))) then
            cycle
         end if

         if (ju == (jb - 1)) then !block on previous cpu, north wall on this
            iypluswall(pn, 1) = n
            iypluswall(pn, 2) = jb
            pn = pn + 1
            cycle
         end if

         if (jl == (je + 1)) then !block on next cpu, south wall on this
            iyminwall(mn, 1) = n
            iyminwall(mn, 2) = je
            mn = mn + 1
            cycle
         end if

         if ((ju < je) .and. (ju >= jb)) then !block & northwall on this cpu   !ILS13, 5.12.17 following Tom
            iypluswall(pn, 1) = n
            iypluswall(pn, 2) = ju + 1
            pn = pn + 1
         end if

         if ((jl > jb) .and. (jl <= je)) then !block & southwall on this cpu
            iyminwall(mn, 1) = n
            iyminwall(mn, 2) = jl - 1
            mn = mn + 1
         end if
      end do

   end subroutine createwalls

   subroutine ibmwallfun
      use modglobal, only:libm
      use modfields, only:momfluxb, tfluxb, qfluxb, thl_flux

      if (libm) then
         ! compute fluxes at IBM
         momfluxb = 0.
         tfluxb = 0.
         qfluxb = 0.

         call xwallfun
         call ywallfunplus ! due to parallellisation differentiation between + and - side
         call ywallfunmin ! due to parallellisation differentiation between + and - side
         call zwallfun

      end if
   end subroutine ibmwallfun

   subroutine xwallfun
      use modglobal, only:dzf, dzhiq, dzhi, dxf, dxfi, dxhi, dyi, lles, nsv, numol, ltempeq, lmoist, &
         ih, jh, kh, ihc, jhc, khc, dxh, dy, dt, totavtime, rk3step, ib, ie, kb, ke, iwallmom, iwalltemp, iwallmoist, iwallscal, nblocks
      use modfields, only:um, up, v0, w0, vp, wp, shear, thl0, thlp, qt0, qtp, sv0, svp, momfluxb, tfluxb, exnf, cth, qfluxb
      use initfac, only:fachf, block, faclGR, facef, facqsat, fachurel, facf, facT, facz0,facz0h
      integer i, j, k, n, nc, jl, ju, kl, ku, im, jm, jp, km, m

      if (iwallmom == 1) then !fixed flux
      !not implemented
      else if (iwallmom == 2) then !wall function
         do n = 1, nxwall
            k = block(ixwall(n), 8) !west side
            if (k /= 0) call wfuno(ih, jh, kh, vp, wp, thlp, momfluxb, tfluxb, cth, bcTfluxA, v0, w0, thl0, facT(k, 1), facz0(k), facz0h(k), ixwall(n), 1, 11)
            k = block(ixwall(n), 9) !east side
            if (k /= 0) call wfuno(ih, jh, kh, vp, wp, thlp, momfluxb, tfluxb, cth, bcTfluxA, v0, w0, thl0, facT(k, 1), facz0(k), facz0h(k), ixwall(n), 1, 21)
         end do
      else if (iwallmom == 3) then
         do n = 1, nxwall
            k = block(ixwall(n), 8) !west side
            if (k /= 0) call wfmneutral(ih, jh, kh, vp, wp, momfluxb, v0, w0, facz0(k), ixwall(n), 1, 11)
            k = block(ixwall(n), 9) !east side
            if (k /= 0) call wfmneutral(ih, jh, kh, vp, wp, momfluxb, v0, w0, facz0(k), ixwall(n), 1, 21)
        end do
      end if

      if (ltempeq) then
         if (iwalltemp == 1) then !fixed flux
            do n = 1, nxwall
               k = block(ixwall(n), 8) !west side
               if (k /= 0) call xwallscalarmin_advecc2nd_corr(ih, jh, kh, thl0, thlp, bctfxm, ixwall(n))
               k = block(ixwall(n), 9) !east side
               if (k /= 0) call xwallscalarplus_advecc2nd_corr(ih, jh, kh, thl0, thlp, bctfxp, ixwall(n))
            end do
         else if (iwalltemp == 2) then
            do n = 1, nxwall
               k = block(ixwall(n), 8) !west side
               if (k /= 0) then
                  call wfuno(ih, jh, kh, vp, wp, thlp, momfluxb, tfluxb, cth, bcTfluxA, v0, w0, thl0, facT(k, 1), facz0(k), facz0h(k), ixwall(n), 1, 12) !left wall
                  fachf(k) = fachf(k) + bcTfluxA !accumulate flux from that facet (can be on multiple processors, will be MPI_ALLREDUCEd in modEB)
               end if
               k = block(ixwall(n), 9) !east side
               if (k /= 0) then
                  call wfuno(ih, jh, kh, vp, wp, thlp, momfluxb, tfluxb, cth, bcTfluxA, v0, w0, thl0, facT(k, 1), facz0(k), facz0h(k), ixwall(n), 1, 22) !right wall
                  fachf(k) = fachf(k) + bcTfluxA
               end if
            end do
         end if
      end if

      if (lmoist) then
         if (iwallmoist == 1) then !fixed flux
            do n = 1, nxwall
               k = block(ixwall(n), 8) !west side
               if (k /= 0) call xwallscalarmin_advecc2nd_corr(ih, jh, kh, qt0, qtp, bcqfxm, ixwall(n))
               k = block(ixwall(n), 9) !east side
               if (k /= 0) call xwallscalarplus_advecc2nd_corr(ih, jh, kh, qt0, qtp, bcqfxp, ixwall(n))
            end do
         end if
         if ((ltempeq) .and. (iwallmoist == 2)) then
            do n = 1, nxwall
               k = block(ixwall(n), 8)
               if (k /= 0) then
                  call wfGR(ih, jh, kh, qtp, qfluxb, cth, bcqfluxA, qt0(:, :, :), facqsat(k), fachurel(k), facf(k, 4), facf(k, 5), ixwall(n), 1, 12) !left wall
                  facef(k) = facef(k) + bcqfluxA
               end if
               k = block(ixwall(n), 9)
               if (k /= 0) then
                  call wfGR(ih, jh, kh, qtp, qfluxb, cth, bcqfluxA, qt0(:, :, :), facqsat(k), fachurel(k), facf(k, 4), facf(k, 5), ixwall(n), 1, 22) !right wall
                  facef(k) = facef(k) + bcqfluxA
               end if
            end do
         end if
      end if

      if (nsv>0) then
         if (iwallscal == 1) then !fixed flux
            do n = 1, nxwall
               do m= 1, nsv
                  call xwallscalar(ihc, jhc, khc, sv0(:,:,:,m), svp(:,:,:,m), 0., 0., ixwall(n))
               end do
            end do
         end if
      end if

   end subroutine xwallfun

   subroutine xwallscalar(hi, hj, hk, putin, putout, bcvaluem, bcvaluep, n)
      use modglobal, only:jmax, dxf, dxfi, dxfi5, dxhi, dxh2i, nsv, ib, ie, jb, je, kb, ke, prandtlmoli, numol
      use modfields, only:u0
      use modmpi, only:myid
      use modsubgriddata, only:ekh
      use initfac, only:block
      integer i, j, k, jl, ju, kl, ku, iww, iee

      integer, intent(in) :: hi !<size of halo in i
      integer, intent(in) :: hj !<size of halo in j
      integer, intent(in) :: hk !<size of halo in k
      real, intent(in)    :: putin(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk)
      real, intent(inout) :: putout(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk)
      real, intent(in)    :: bcvaluem, bcvaluep
      integer, intent(in)    :: n

      iww = block(n, 1) - 1
      iee = block(n, 2) + 1

      jl = MAX(block(n, 3) - myid*jmax, 1) ! starting j-index
      ju = MIN(block(n, 4) - myid*jmax, jmax) ! ending j-index
      kl = block(n, 5) ! starting k-index
      ku = block(n, 6) ! ending k-index

      !fixed flux
      !remove standard diffusion term, add flux=bcvalue
      do k = kl, ku
         do j = jl, ju
            putout(iee, j, k) = putout(iee, j, k) + ( &
                                0.5*(ekh(iee, j, k)*dxf(iee - 1) + ekh(iee - 1, j, k)*dxf(iee))* &
                                (putin(iee, j, k) - putin(iee - 1, j, k))*dxh2i(iee) - &
                                bcvaluep)*dxfi(iee) !

            putout(iww, j, k) = putout(iww, j, k) + ( &
                                -0.5*(ekh(iww + 1, j, k)*dxf(iww) + ekh(iww, j, k)*dxf(iww + 1))* &
                                (putin(iww + 1, j, k) - putin(iww, j, k))*dxh2i(iww + 1) - &
                                bcvaluem)*dxfi(iww) !
         end do
      end do

   end subroutine xwallscalar


   subroutine xwallscalarmin_advecc2nd_corr(hi, hj, hk, putin, putout, bcvaluem, n)
      ! 12
      use modglobal, only:jmax, dxf, dxfi, dxfi5, dxhi, dxh2i, nsv, ib, ie, jb, je, kb, ke, prandtlmoli, numol
      use modfields, only:u0
      use modmpi, only:myid
      use modsubgriddata, only:ekh
      use initfac, only:block
      integer i, ip, j, k, jl, ju, kl, ku

      integer, intent(in) :: hi !<size of halo in i
      integer, intent(in) :: hj !<size of halo in j
      integer, intent(in) :: hk !<size of halo in k
      real, intent(in)    :: putin(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk)
      real, intent(inout) :: putout(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk)
      real, intent(in)    :: bcvaluem
      integer, intent(in)    :: n

      i = block(n, 1) - 1
      ip = i + 1
      jl = MAX(block(n, 3) - myid*jmax, 1) ! starting j-index
      ju = MIN(block(n, 4) - myid*jmax, jmax) ! ending j-index
      kl = block(n, 5) ! starting k-index
      ku = block(n, 6) ! ending k-index

      !fixed flux
      !remove standard diffusion term and 2nd order advection term, add flux=bcvalue
      do k = kl, ku
         do j = jl, ju
            putout(i, j, k) = putout(i, j, k) - bcvaluem*dxfi(i) &
                            - 0.5*(ekh(ip, j, k)*dxf(i) + ekh(i, j, k)*dxf(ip))*(putin(ip, j, k) - putin(i, j, k))*dxh2i(ip)*dxfi(i) &
                            + (u0(ip, j, k)*(putin(ip, j, k)*dxf(i) + putin(i, j, k)*dxf(ip))*dxhi(ip))*dxfi5(i) &
                            - (u0(ip, j, k)*(putin(i , j, k)*dxf(i) + putin(i, j, k)*dxf(ip))*dxhi(ip))*dxfi5(i)
         end do
      end do

   end subroutine xwallscalarmin_advecc2nd_corr


   subroutine xwallscalarplus_advecc2nd_corr(hi, hj, hk, putin, putout, bcvaluep, n)
      ! 22
      use modglobal, only:jmax, dxf, dxfi, dxfi5, dxhi, dxh2i, nsv, ib, ie, jb, je, kb, ke, prandtlmoli, numol
      use modfields, only:u0
      use modmpi, only:myid
      use modsubgriddata, only:ekh
      use initfac, only:block
      integer i, im, j, k, jl, ju, kl, ku

      integer, intent(in) :: hi !<size of halo in i
      integer, intent(in) :: hj !<size of halo in j
      integer, intent(in) :: hk !<size of halo in k
      real, intent(in)    :: putin(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk)
      real, intent(inout) :: putout(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk)
      real, intent(in)    :: bcvaluep
      integer, intent(in)    :: n

      i = block(n, 2) + 1
      im = i - 1
      jl = MAX(block(n, 3) - myid*jmax, 1) ! starting j-index
      ju = MIN(block(n, 4) - myid*jmax, jmax) ! ending j-index
      kl = block(n, 5) ! starting k-index
      ku = block(n, 6) ! ending k-index

      !fixed flux
      !remove standard diffusion term and 2nd order advection term, add flux=bcvalue
      do k = kl, ku
         do j = jl, ju
            putout(i, j, k) = putout(i, j, k) - bcvaluep*dxfi(i) & ! applied flux
                            + 0.5*(ekh(i, j, k)*dxf(im) + ekh(im,j,k)*dxf(i))*(putin(i, j, k) - putin(im, j, k))*dxh2i(i)*dxfi(i) &
                            + (- u0(i, j, k)*(putin(im, j, k)*dxf(i) + putin(i, j, k)*dxf(im))*dxhi(i))*dxfi5(i) &
                            - (- u0(i, j, k)*(putin(i , j, k)*dxf(i) + putin(i, j, k)*dxf(im))*dxhi(i))*dxfi5(i)
         end do
      end do

   end subroutine xwallscalarplus_advecc2nd_corr


   subroutine ywallfunplus
      use modglobal, only:dzf, dzhiq, dzhi, dxf, dxhi, dy, dyi, nsv, lles, numol, ltempeq, lmoist, &
         je, jb, ih, jh, kh, ihc, jhc, khc, iwallmom, iwallmoist, iwalltemp, iwallscal, nblocks
      use modfields, only:u0, w0, up, wp, shear, thlp, thl0, qtp, qt0, sv0, svp, tfluxb, momfluxb, exnf, cth, qfluxb
      use modsubgriddata, only:ekm
      use modmpi, only:myid
      use initfac, only:fachf, block, faclGR, facqsat, facef, fachurel, facf, facT, facz0, facz0h
      integer i, j, k, n, nc, il, iu, kl, ku, im, jm, km, m

      if (iwallmom == 1) then !fixed flux
      !not implemented
      else if (iwallmom == 2) then
         do n = 1, nypluswall
            k = block(iypluswall(n, 1), 10) !upper y wall = north wall
            if (k /= 0) call wfuno(ih, jh, kh, up, wp, thlp, momfluxb, tfluxb, cth, bcTfluxA, u0, w0, thl0, facT(k,1), facz0(k), facz0h(k), iypluswall(n, 1), iypluswall(n, 2), 31)
         end do
      else if (iwallmom == 3) then
         do n = 1, nypluswall
            k = block(iypluswall(n, 1), 10) !upper y wall = north wall
            if (k /= 0) call wfmneutral(ih, jh, kh, up, wp, momfluxb, u0, w0, facz0(k), iypluswall(n, 1), iypluswall(n, 2), 31)
         end do
      end if

      if (ltempeq) then
         if (iwalltemp == 1) then
            do n = 1, nypluswall ! loop over all shear x-walls
               k = block(iypluswall(n, 1), 10)
               if (k /= 0) call ywallscalarplus_advecc2nd_corr(ih, jh, kh, thl0, thlp, bctfyp, n)
            end do
         else if (iwalltemp == 2) then
            do n = 1, nypluswall
               k = block(iypluswall(n, 1), 10)
               if (k /= 0) then
                  call wfuno(ih, jh, kh, up, wp, thlp, momfluxb, tfluxb, cth, bcTfluxA, u0, w0, thl0, facT(k,1), facz0(k), facz0h(k), iypluswall(n, 1), iypluswall(n, 2), 32)
                  fachf(k) = fachf(k) + bcTfluxA
               end if
            end do
         end if
      end if

      if (lmoist) then
         if (iwallmoist == 1) then
            do n = 1, nypluswall ! loop over all shear x-walls
               k = block(iypluswall(n, 1), 10)
               if (k /= 0) call ywallscalarplus_advecc2nd_corr(ih, jh, kh, qt0, qtp, bcqfyp, n)
            end do
         end if
         if ((ltempeq) .and. (iwallmoist == 2)) then
            do n = 1, nypluswall
               k = block(iypluswall(n, 1), 10)
               if (k /= 0) then
                  call wfGR(ih,jh,kh,qtp,qfluxb,cth,bcqfluxA,qt0,facqsat(k),fachurel(k),facf(k,4),facf(k,5),iypluswall(n,1),iypluswall(n,2),32)
                  facef(k) = facef(k) + bcqfluxA
               end if
            end do
         end if
      end if

      if (nsv>0) then
         if (iwallscal == 1) then
            do n = 1, nypluswall ! loop over all shear x-walls
               do m = 1, nsv
                  call ywallscalarplus(ihc, jhc, khc, sv0(:,:,:,m), svp(:,:,:,m), 0., n)
               end do
            end do
         end if
      end if

   end subroutine ywallfunplus


   subroutine ywallscalarplus(hi, hj, hk, putin, putout, bcvaluep, n)
      use modglobal, only:dyi, ib, ie, jb, je, kb, ke, numol, prandtlmoli
      use modsubgriddata, only:ekh
      use modmpi, only:myid
      use initfac, only:block
      integer i, j, k, il, iu, kl, ku, m

      integer, intent(in) :: hi !<size of halo in i
      integer, intent(in) :: hj !<size of halo in j
      integer, intent(in) :: hk !<size of halo in k
      real, intent(in)    :: putin(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk)
      real, intent(inout) :: putout(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk)
      real, intent(in)    :: bcvaluep
      integer, intent(in)    :: n
      m = iypluswall(n, 1)
      j = iypluswall(n, 2)
      il = block(m, 1)
      iu = block(m, 2)
      kl = block(m, 5)
      ku = block(m, 6)

      !fixed flux
      do k = kl, ku
         do i = il, iu
            putout(i, j, k) = putout(i, j, k) + (0.5*(ekh(i, j, k) + ekh(i, j - 1, k))*(putin(i, j, k) - putin(i, j - 1, k))*dyi - bcvaluep)*dyi
         end do
      end do
   end subroutine ywallscalarplus


   subroutine ywallscalarplus_advecc2nd_corr(hi, hj, hk, putin, putout, bcvaluep, n)
      ! 32
      use modglobal, only : dyi, dy2i, dyi5, ib, ie, jb, je, kb, ke, numol, prandtlmoli
      use modfields, only : v0
      use modsubgriddata, only:ekh
      use modmpi, only:myid
      use initfac, only:block
      integer i, j, jm, k, il, iu, kl, ku, m

      integer, intent(in) :: hi !<size of halo in i
      integer, intent(in) :: hj !<size of halo in j
      integer, intent(in) :: hk !<size of halo in k
      real, intent(in)    :: putin(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk)
      real, intent(inout) :: putout(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk)
      real, intent(in)    :: bcvaluep
      integer, intent(in)    :: n

      m = iypluswall(n, 1)
      j = iypluswall(n, 2)
      jm = j - 1
      il = block(m, 1)
      iu = block(m, 2)
      kl = block(m, 5)
      ku = block(m, 6)

      !fixed flux
      do k = kl, ku
         do i = il, iu
            putout(i, j, k) = putout(i, j, k) - bcvaluep*dyi &
                            + (0.5*(ekh(i, j, k) + ekh(i, jm, k))*(putin(i, j, k) - putin(i, jm, k)))*dy2i &
                            + (- v0(i, j, k)*(putin(i, jm, k) + putin(i, j, k)))*dyi5 &
                            - (- v0(i, j, k)*(putin(i, j , k) + putin(i, j, k)))*dyi5
         end do
      end do

   end subroutine ywallscalarplus_advecc2nd_corr


   subroutine ywallfunmin
      use modglobal, only:dxf, dxhi, dy, dyi, dzhiq, dzf, dzhi, lles, nsv, numol, ltempeq, lmoist, &
         ih, jh, kh, ihc, jhc, khc, iwallmom, iwalltemp, iwallmoist, iwallscal, nblocks
      use modfields, only:u0, w0, up, wp, shear, thl0, thlp, qt0, qtp, sv0, svp, tfluxb, momfluxb, exnf, cth, qfluxb
      use initfac, only:fachf, block, faclGR, facqsat, facef, fachurel, facf, facT, facz0, facz0h
      !      use modsurfdata,     only : wtsurf
      integer i, j, k, n, nc, il, iu, kl, ku, im, jp, km, m

      if (iwallmom == 1) then
      !fixed flux, not implemented
      else if (iwallmom == 2) then
         do n = 1, nyminwall
         k = block(iyminwall(n, 1), 11)
         if (k /= 0) call wfuno(ih, jh, kh, up, wp, thlp, momfluxb, tfluxb, cth, bcTfluxA, u0, w0, thl0, facT(k,1), facz0(k), facz0h(k), iyminwall(n, 1), iyminwall(n, 2), 41)
         end do
      else if (iwallmom == 3) then
         do n = 1, nyminwall
         k = block(iyminwall(n, 1), 11)
         if (k /= 0) call wfmneutral(ih, jh, kh, up, wp, momfluxb, u0, w0, facz0(k), iyminwall(n, 1), iyminwall(n, 2), 41)
         end do
      end if !

      if (ltempeq) then
         if (iwalltemp == 1) then
            do n = 1, nyminwall !
               k = block(iyminwall(n, 1), 11)
               if (k /= 0) call ywallscalarmin_advecc2nd_corr(ih, jh, kh, thl0, thlp, bctfym, n)
            end do
         else if (iwalltemp == 2) then
            do n = 1, nyminwall
               k = block(iyminwall(n, 1), 11)
               if (k /= 0) then
                  call wfuno(ih, jh, kh, up, wp, thlp, momfluxb, tfluxb, cth, bcTfluxA, u0, w0, thl0, facT(k,1), facz0(k), facz0h(k), iyminwall(n, 1), iyminwall(n, 2), 42)
                  fachf(k) = fachf(k) + bcTfluxA
               end if
            end do
         end if
      end if

      if (lmoist) then
         if (iwallmoist == 1) then
            do n = 1, nyminwall !
               k = block(iyminwall(n, 1), 11)
               if (k /= 0) call ywallscalarmin_advecc2nd_corr(ih, jh, kh, qt0, qtp, bcqfym, n)
            end do
         end if
         if ((ltempeq) .and. (iwallmoist == 2)) then
            do n = 1, nyminwall
               k = block(iyminwall(n, 1), 11)
               if (k /= 0) then
               call wfGR(ih, jh, kh, qtp, qfluxb, cth, bcqfluxA, qt0, facqsat(k),fachurel(k),facf(k,4),facf(k,5),iyminwall(n, 1), iyminwall(n, 2), 42)
                  facef(k) = facef(k) + bcqfluxA
               end if
            end do
         end if
      end if

      if (nsv>0) then
         if (iwallscal == 1) then
            do n = 1, nyminwall !
               do m = 1, nsv
                  call ywallscalarmin(ihc, jhc, khc, sv0(:,:,:,m), svp(:,:,:,m), 0., n)
               end do
            end do
         end if
      end if

   end subroutine ywallfunmin


   subroutine ywallscalarmin(hi, hj, hk, putin, putout, bcvaluem, n)
      use modglobal, only : dyi, ib, ie, jb, je, kb, ke, prandtlmoli, numol
      use modsubgriddata, only:ekh
      use modmpi, only:myid
      use initfac, only:block
      integer i, j, k, il, iu, kl, ku, m

      integer, intent(in) :: hi !<size of halo in i
      integer, intent(in) :: hj !<size of halo in j
      integer, intent(in) :: hk !<size of halo in k
      real, intent(in)    :: putin(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk)
      real, intent(inout) :: putout(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk)
      real, intent(in)    :: bcvaluem
      integer, intent(in)    :: n

      j = iyminwall(n, 2)
      m = iyminwall(n, 1)
      il = block(m, 1)
      iu = block(m, 2)
      kl = block(m, 5)
      ku = block(m, 6)

      do k = kl, ku
         do i = il, iu
            putout(i, j, k) = putout(i, j, k) + ( &
                              -0.5*(ekh(i, j, k) + ekh(i, j + 1, k))*(putin(i, j + 1, k) - putin(i, j, k))*dyi &
                              - bcvaluem)*dyi
         end do
      end do

   end subroutine ywallscalarmin


   subroutine ywallscalarmin_advecc2nd_corr(hi, hj, hk, putin, putout, bcvaluem, n)
      ! 42
      use modglobal, only : dyi, dy2i, dyi5, ib, ie, jb, je, kb, ke, prandtlmoli, numol
      use modfields, only : v0
      use modsubgriddata, only : ekh
      use modmpi,  only : myid
      use initfac, only : block
      integer i, j, jp, k, il, iu, kl, ku, m

      integer, intent(in) :: hi !<size of halo in i
      integer, intent(in) :: hj !<size of halo in j
      integer, intent(in) :: hk !<size of halo in k
      real, intent(in)    :: putin(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk)
      real, intent(inout) :: putout(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk)
      real, intent(in)    :: bcvaluem
      integer, intent(in)    :: n

      j = iyminwall(n, 2)
      jp = j + 1
      m = iyminwall(n, 1)
      il = block(m, 1)
      iu = block(m, 2)
      kl = block(m, 5)
      ku = block(m, 6)

      do k = kl, ku
         do i = il, iu
            putout(i, j, k) = putout(i, j, k) - bcvaluem*dyi &
                            - 0.5*(ekh(i, j, k) + ekh(i, jp, k))*(putin(i, jp, k) - putin(i, j, k))*dy2i &
                            + v0(i, jp, k)*(putin(i, jp, k) + putin(i, j, k))*dyi5 &
                            - v0(i, jp, k)*(putin(i, j , k) + putin(i, j, k))*dyi5

         end do
      end do

   end subroutine ywallscalarmin_advecc2nd_corr


   subroutine zwallfun
      use modglobal, only:dzf, dzfi, dzhi, dzhiq, dxf, dxfi, dxhi, dyi, nsv, lles, numol, ltempeq, lmoist, &
         ib, ie, jb, je, kb, ke, ih, jh, kh, ihc, jhc, khc, iwallmom, iwalltemp, iwallmoist, iwallscal, fkar, grav
      use modfields, only:u0, v0, up, vp, shear, thl0, thlp, qt0, qtp, sv0, svp, tfluxb, momfluxb, exnf, cth, qfluxb, IIc, IIids
      use modmpi, only:myid
      use initfac, only:fachf, block, faclGR, facef, facqsat, fachurel, facf, facT, facz0, facz0h
      use modsubgriddata, only: ekm
      integer i, j, k, n, nc, il, iu, jl, ju, im, jm, km, m, kdum
      real :: adv_contq, utang1Intdum, utang2Intdum,utangIntdum, dTdum, Ribl0dum, deltadum,fkar2dum, &
              z0dum,z0hdum,logdzdum, logdzhdum, logzhdum, sqdzdum, ctmdum, dummydum, bcmomfluxdum, eommdum, eompdum

      fkar2dum = fkar**2

      if (iwallmom == 1) then
      !fixed flux
      else if (iwallmom == 2) then
         do n = 1, nxwall
            k = block(ixwall(n), 7)
            if (k /= 0) call wfuno(ih, jh, kh, up, vp, thlp, momfluxb, tfluxb, cth, bcTfluxA, u0, v0, thl0, facT(k, 1), facz0(k), facz0h(k), ixwall(n), 1, 51)
            ! k = block(ixwall(n), 12)
            ! if (k /= 0) call wfuno(ih, jh, kh, up, vp, thlp, momfluxb, tfluxb, cth, bcTfluxA, u0, v0, thl0, facT(k, 1), facz0(k), facz0h(k), ixwall(n), 1, 61)
         end do

         ! Case 51
         j = jb
         do i=ib,ie
           do k = kb+1,ke
             if(IIc(i,j-1,k-1)==0 .and. IIc(i,j,k)==1 .and. IIc(i,j-1,k)==1 .and. IIc(i,j,k-1)==1) THEN
               kdum = IIids(i,j-1,k-1)
               z0dum = facz0(kdum)
               z0hdum = facz0h(kdum)
               deltadum = 0.5*dzf(k)
               logdzdum = LOG(deltadum/z0dum)
               logdzhdum = LOG(deltadum/z0hdum)
               logzhdum = LOG(z0dum/z0hdum)
               sqdzdum = SQRT(deltadum/z0dum)
               utang1Intdum = (u0(i,j,k)+u0(i,j-1,k)+u0(i + 1, j - 1, k)+u0(i + 1, j, k))*0.25
               utang2Intdum = v0(i,j,k)
               utangIntdum = max(0.0001,(utang1Intdum**2 + utang2Intdum**2))
               dTdum = thl0(i,j-1,k) - facT(kdum,1)
               Ribl0dum = grav*deltadum*dTdum/(facT(kdum,1)*utangIntdum)
               ctmdum = unom(logdzdum, logdzhdum, logzhdum, sqdzdum, utangIntdum, dTdum, Ribl0dum, fkar2dum)
               !dummydum = (utang2Intdum**2)*ctmdum
               dummydum = abs(utang2Intdum)*sqrt(utangIntdum)*ctmdum
               bcmomfluxdum = SIGN(dummydum, utang2Intdum)
               eommdum = (dzf(k-1)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, k-1) + ekm(i, j - 1, k-1)))*dzhiq(k)
               vp(i, j, k) = vp(i, j, k) + (v0(i, j, k) - v0(i, j, k-1))*eommdum*dzhi(k)*dzfi(k) - bcmomfluxdum*dzfi(k)*0.5
             end if

             ! ! Case 61
             ! if(IIc(i,j-1,k+1)==0 .and. IIc(i,j,k)==1 .and. IIc(i,j-1,k)==1 .and. IIc(i,j,k+1)==1) THEN
             !   kdum = IIids(i,j-1,k+1)
             !   z0dum = facz0(kdum)
             !   z0hdum = facz0h(kdum)
             !   deltadum = 0.5*dzf(k)
             !   logdzdum = LOG(deltadum/z0dum)
             !   logdzhdum = LOG(deltadum/z0hdum)
             !   logzhdum = LOG(z0dum/z0hdum)
             !   sqdzdum = SQRT(deltadum/z0dum)
             !   utang1Intdum = (u0(i,j,k)+u0(i,j-1,k)+u0(i + 1, j - 1, k)+u0(i + 1, j, k))*0.25
             !   utang2Intdum = v0(i,j,k)
             !   utangIntdum = max(0.0001,(utang1Intdum**2 + utang2Intdum**2))
             !   dTdum = thl0(i,j-1,k) - facT(kdum,1)
             !   Ribl0dum = grav*deltadum*dTdum/(facT(kdum,1)*utangIntdum)
             !   ctmdum = unom(logdzdum, logdzhdum, logzhdum, sqdzdum, utangIntdum, dTdum, Ribl0dum, fkar2dum)
             !   dummydum = (utang2Intdum**2)*ctmdum
             !   bcmomfluxdum = SIGN(dummydum, utang2Intdum)
             !   eompdum = (dzf(k+1)*(ekm(i,j,k) + ekm(i,j-1,k)) + dzf(k)*(ekm(i,j,k+1) + ekm(i,j-1,k+1)))*dzhiq(k+1)
             !   vp(i, j, k) = vp(i, j, k) - (v0(i, j, k+1) - v0(i, j, k))*eompdum*dzhi(k)*dzfi(k) - bcmomfluxdum*dzfi(k)*0.5
             ! end if
           end do
        end do

      else if (iwallmom == 3) then
         do n = 1, nxwall
            k = block(ixwall(n), 7)
            if (k /= 0) call wfmneutral(ih, jh, kh, up, vp, momfluxb, u0, v0, facz0(k), ixwall(n), 1, 51)
            ! k = block(ixwall(n),12)
            ! if (k /= 0) call wfmneutral(ih, jh, kh, up, vp, momfluxb, u0, v0, facz0(k), ixwall(n), 1, 61)
         end do

         ! Case 51
         j = jb
         do i=ib,ie
           do k = kb+1,ke
             if(IIc(i,j-1,k-1)==0 .and. IIc(i,j,k)==1 .and. IIc(i,j-1,k)==1 .and. IIc(i,j,k-1)==1) THEN
               kdum = IIids(i,j-1,k-1)
               z0dum = facz0(kdum)
               deltadum = 0.5*dzf(k)
               logdzdum = LOG(deltadum/z0dum)
               utang2Intdum = v0(i,j,k)
               ctmdum = fkar**2/logdzdum**2
               !dummydum = (utang2Intdum**2)*ctmdum
               dummydum = abs(utang2Intdum)*sqrt(utangIntdum)*ctmdum
               bcmomfluxdum = SIGN(dummydum, utang2Intdum)
               eommdum = (dzf(k-1)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, k-1) + ekm(i, j - 1, k-1)))*dzhiq(k)
               vp(i, j, k) = vp(i, j, k) + (v0(i, j, k) - v0(i, j, k-1))*eommdum*dzhi(k)*dzfi(k) - bcmomfluxdum*dzfi(k)*0.5
             end if

             ! ! Case 61
             ! if(IIc(i,j-1,k+1)==0 .and. IIc(i,j,k)==1 .and. IIc(i,j-1,k)==1 .and. IIc(i,j,k+1)==1) THEN
             !   kdum = IIids(i,j-1,k+1)
             !   z0dum = facz0(kdum)
             !   deltadum = 0.5*dzf(k)
             !   logdzdum = LOG(deltadum/z0dum)
             !   utang2Intdum = v0(i,j,k)
             !   ctmdum = fkar**2/logdzdum**2
             !   dummydum = (utang2Intdum**2)*ctmdum
             !   bcmomfluxdum = SIGN(dummydum, utang2Intdum)
             !   eompdum = (dzf(k+1)*(ekm(i,j,k) + ekm(i,j-1,k)) + dzf(k)*(ekm(i,j,k+1) + ekm(i,j-1,k+1)))*dzhiq(k+1)
             !   vp(i, j, k) = vp(i, j, k) - (v0(i, j, k+1) - v0(i, j, k))*eompdum*dzhi(k)*dzfi(k) - bcmomfluxdum*dzfi(k)*0.5
             ! end if
           end do
        end do

      end if

      if (ltempeq) then
         if (iwalltemp == 1) then
            do n = 1, nxwall ! loop over all shear x-walls
               k = block(ixwall(n), 7)
               if (k /= 0) call zwallscalarplus_advecc2nd_corr(ih, jh, kh, thl0, thlp, bctfz, ixwall(n))
               ! k = block(ixwall(n), 12)
               ! if (k /= 0) call zwallscalarmin_advecc2nd_corr(ih, jh, kh, thl0, thlp, 0., ixwall(n))
            end do
         else if (iwalltemp == 2) then
            do n = 1, nxwall
               k = block(ixwall(n), 7)
               if (k /= 0) then
                  call wfuno(ih, jh, kh, up, vp, thlp, momfluxb, tfluxb, cth, bcTfluxA, u0, v0, thl0, facT(k, 1), facz0(k), facz0h(k), ixwall(n), 1, 52)
                  fachf(k) = fachf(k) + bcTfluxA
               end if
               ! k = block(ixwall(n), 12)
               ! if (k /= 0) then
               !    call wfuno(ih, jh, kh, up, vp, thlp, momfluxb, tfluxb, cth, bcTfluxA, u0, v0, thl0, facT(k, 1), facz0(k), facz0h(k), ixwall(n), 1, 62)
               !    fachf(k) = fachf(k) + bcTfluxA
               ! end if
            end do
         end if
      end if

      if (lmoist) then
         if (iwallmoist == 1) then
            do n = 1, nxwall ! loop over all shear x-walls
               k = block(ixwall(n), 7)
               if (k /= 0) call zwallscalarplus_advecc2nd_corr(ih, jh, kh, qt0, qtp, bcqfz, ixwall(n))
               ! k = block(ixwall(n), 12)
               ! if (k /= 0) call zwallscalarmin_advecc2nd_corr(ih, jh, kh, qt0, qtp, 0., ixwall(n))
            end do
         end if
         if ((ltempeq) .and. (iwallmoist == 2)) then
            do n = 1, nxwall
               k = block(ixwall(n), 7)
               if (k /= 0) then
                  call wfGR(ih, jh, kh, qtp, qfluxb, cth, bcqfluxA, qt0, facqsat(k), fachurel(k), facf(k, 4), facf(k, 5), ixwall(n), 1, 52)
                  facef(k) = facef(k) + bcqfluxA
               end if
               k = block(ixwall(n), 7)
               ! if (k /= 0) then
               !    call wfGR(ih, jh, kh, qtp, qfluxb, cth, bcqfluxA, qt0, facqsat(k), fachurel(k), facf(k, 4), facf(k, 5), ixwall(n), 1, 62)
               !    facef(k) = facef(k) + bcqfluxA
               ! end if
            end do
         end if
      end if

      if (nsv>0) then
         if (iwallscal == 1) then
            do n = 1, nxwall ! loop over all shear x-walls
               do m = 1, nsv
                  call zwallscalar(ihc, jhc, khc, sv0(:,:,:,m), svp(:,:,:,m), 0., ixwall(n))
               end do
            end do
         end if
      end if

   end subroutine zwallfun


   subroutine zwallscalar(hi, hj, hk, putin, putout, bcvalue, n)
      use modglobal, only:jmax, dzf, dzfi, dzhi, dzh2i, ib, ie, jb, je, kb, ke, prandtlmoli, numol
      use modsubgriddata, only:ekh
      use modmpi, only:myid
      use initfac, only:block
      integer i, j, k, il, iu, jl, ju, km
      integer, intent(in) :: hi !<size of halo in i
      integer, intent(in) :: hj !<size of halo in j
      integer, intent(in) :: hk !<size of halo in k
      real, intent(in)    :: putin(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk)
      real, intent(inout) :: putout(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk)
      real, intent(in)    :: bcvalue
      integer, intent(in)    :: n

      k = block(n, 6) + 1 !block location
      km = k - 1 !
      il = block(n, 1)
      iu = block(n, 2)
      jl = MAX(block(n, 3) - myid*jmax, 1)
      ju = MIN(block(n, 4) - myid*jmax, jmax)

      !  delta=putout(3,1,2)
      do j = jl, ju
         do i = il, iu
            putout(i, j, k) = putout(i, j, k) + ( &
                              0.5*(dzf(km)*ekh(i, j, k) + dzf(k)*ekh(i, j, km))* & ! zero flux
                              (putin(i, j, k) - putin(i, j, km))*dzh2i(k) - &
                              bcvalue)*dzfi(k)
         end do
      end do
   end subroutine zwallscalar


   subroutine zwallscalarplus_advecc2nd_corr(hi, hj, hk, putin, putout, bcvalue, n)
      ! 52
      use modglobal, only : jmax, dzf, dzfi, dzhi, dzh2i, dzfi5, ib, ie, jb, je, kb, ke, prandtlmoli, numol
      use modfields, only : w0
      use modsubgriddata, only : ekh
      use modmpi,  only : myid
      use initfac, only : block
      integer i, j, k, il, iu, jl, ju, km
      integer, intent(in) :: hi !<size of halo in i
      integer, intent(in) :: hj !<size of halo in j
      integer, intent(in) :: hk !<size of halo in k
      real, intent(in)    :: putin(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk)
      real, intent(inout) :: putout(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk)
      real, intent(in)    :: bcvalue
      integer, intent(in)    :: n

      k = block(n, 6) + 1 !block location
      if (k > ke) return
      km = k - 1 !
      il = block(n, 1)
      iu = block(n, 2)
      jl = MAX(block(n, 3) - myid*jmax, 1)
      ju = MIN(block(n, 4) - myid*jmax, jmax)

      do j = jl, ju
         do i = il, iu
            putout(i, j, k) = putout(i, j, k) - bcvalue*dzfi(k) &
                            + 0.5*(dzf(km)*ekh(i, j, k) + dzf(k)*ekh(i, j, km))*(putin(i, j, k) - putin(i, j, km))*dzh2i(k)*dzfi(k) &
                            + (- w0(i, j, k)*(putin(i, j, km)*dzf(k) + putin(i, j, k)*dzf(km))*dzhi(k))*dzfi5(k) &
                            - (- w0(i, j, k)*(putin(i, j, k )*dzf(k) + putin(i, j, k)*dzf(km))*dzhi(k))*dzfi5(k)
         end do
      end do

   end subroutine zwallscalarplus_advecc2nd_corr


   ! subroutine zwallscalarmin_advecc2nd_corr(hi, hj, hk, putin, putout, bcvalue, n)
   !    ! 62
   !    use modglobal, only : jmax, dzf, dzfi, dzhi, dzh2i, dzfi5, ib, ie, jb, je, kb, ke, prandtlmoli, numol
   !    use modfields, only : w0
   !    use modsubgriddata, only : ekh
   !    use modmpi,  only : myid
   !    use initfac, only : block
   !    integer i, j, k, kp, il, iu, jl, ju
   !    integer, intent(in) :: hi !<size of halo in i
   !    integer, intent(in) :: hj !<size of halo in j
   !    integer, intent(in) :: hk !<size of halo in k
   !    real, intent(in)    :: putin(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk)
   !    real, intent(inout) :: putout(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk)
   !    real, intent(in)    :: bcvalue
   !    integer, intent(in)    :: n
   !
   !    k = block(n, 5) - 1 !block location
   !    if (k < kb) return
   !    kp = k + 1 !
   !    il = block(n, 1)
   !    iu = block(n, 2)
   !    jl = MAX(block(n, 3) - myid*jmax, 1)
   !    ju = MIN(block(n, 4) - myid*jmax, jmax)
   !
   !    do j = jl, ju
   !       do i = il, iu
   !          putout(i, j, k) = putout(i, j, k) - bcvalue*dzfi(k) &
   !                         - 0.5*(dzf(kp)*ekh(i, j, k) + dzf(k)*ekh(i, j, kp))*(putin(i, j, kp) - putin(i, j, k))*dzh2i(k)*dzfi(k) &
   !                         + w0(i, j, kp)*(putin(i, j, kp)*dzf(k) + putin(i, j, k)*dzf(kp))*dzhi(kp)*dzfi5(k) &
   !                         - w0(i, j, kp)*(putin(i, j, k )*dzf(k) + putin(i, j, k)*dzf(kp))*dzhi(kp)*dzfi5(k)
   !       end do
   !    end do
   !
   ! end subroutine zwallscalarmin_advecc2nd_corr


   subroutine ibmnorm
      use modglobal, only:ib, ie, ih, jb, je, jh, kb, ke, kh, zh, rk3step, dt, libm, jmax, &
         nblocks, nsv, ltempeq, lmoist, rk3step, ih, kh, dt, totavtime, &
         dxh, dzf, dy, ih, kh, jh, jge, bldT
      use modfields, only:up, vp, wp, um, vm, wm, u0, v0, w0, thl0, thl0av, thlm, svp, svm, thlp, qtp, qt0, qtm, IIc, IIids
      use modmpi, only:myid, nprocs
      use initfac, only:block
      real, dimension(ib - ih:ie + ih, kb - kh:ke + kh)          ::  dummy
      real rk3coef, rk3coefi, thl0volav
      integer n, i, j, k, il, iu, jl, ju, kl, ku, sc

      if (libm) then
         rk3coef = dt/(4.-dble(rk3step))
         rk3coefi = 1./rk3coef

         do n = 1, nxwall
            il = block(ixwall(n), 1)
            iu = block(ixwall(n), 2) + 1
            jl = max(block(ixwall(n), 3) - myid*jmax, 1) !
            ju = min(block(ixwall(n), 4) - myid*jmax, jmax) !
            !kl = block(ixwall(n), 5)
            kl = kb
            ! tg3315 18.03.19 - use kb because for lEB buildings block starts at kb+1 but this leaves area underneath the buildings and horizontally between the roads where we have no block. Only leads to small velocities in these areas but this negates this issue. WARNING - for modelling overhangs this should be changed but this would also require another facade type etc. Similarly applied to y and z directions below.
            ku = block(ixwall(n), 6)

            !up(il:iu, jl:ju, kl:ku) = -um(il:iu, jl:ju, kl:ku)*rk3coefi
            up(iu, jl:ju, kl:ku) = -um(iu, jl:ju, kl:ku)*rk3coefi
            up(il, jl:ju, kl:ku) = -um(il, jl:ju, kl:ku)*rk3coefi
            up(il + 1:iu - 1, jl:ju, kl:ku) = 0. !internal velocity don't change or
            um(il + 1:iu - 1, jl:ju, kl:ku) = 0. !internal velocity = 0    or both?
         end do ! 1,nxwallsnorm

         do n = 1, nblocks
           jl = block(n,3)
           ju = block(n,4)
           il = block(n, 1)
           iu = block(n, 2)
           kl = kb ! tg3315 see comment for x-direction above
           !kl = block(n,5)
           ku = block(n, 6)
           if (jl>je+myid*je .or. ju<jb+myid*je) THEN
             cycle
           else if (ju<=je+myid*je .and. jl>=jb+myid*je) then
             jl = jl - myid*je
             ju = ju - myid*je
             do i = il,iu
               do k = kl,ku
                   if (IIc(i,jl-1,k) == 0) THEN
                     vp(i, jl, k) = 0.0
                     vm(i, jl, k) = 0.0

                   else
                     vp(i, jl, k) = -vm(i, jl,k)*rk3coefi
                   end if
                   if(IIc(i,ju+1,k) == 0) then
                     vp(i, ju+1, k) = 0.0
                     vm(i, ju+1, k) = 0.0
                   else
                     vp(i, ju+1, k) = -vm(i, ju+1,k)*rk3coefi
                   end if
               end do
             end do
             vp(il:iu, jl + 1:ju, kl:ku) = 0.0 !setting the internal cells
             vm(il:iu, jl + 1:ju, kl:ku) = 0.0
           else if  (ju>je+myid*je .and. jl<jb+myid*je) then
               jl = jb
               ju = je
               vp(il:iu, jl:ju, kl:ku) = 0.0 !setting the internal cells
               vm(il:iu, jl:ju, kl:ku) = 0.0
           else if  (ju<=je+myid*je .and. jl<jb+myid*je) then
             jl = jb
             ju = ju - myid*je
             do i = il,iu
               do k = kl,ku
                   if(IIc(i,ju+1,k) == 0) then
                     vp(i, ju+1, k) = 0.0
                     vm(i, ju+1, k) = 0.0
                   else
                     vp(i, ju+1, k) = -vm(i, ju+1,k)*rk3coefi
                   end if
               end do
             end do
             vp(il:iu, jl:ju, kl:ku) = 0.0 !setting the internal cells
             vm(il:iu, jl:ju, kl:ku) = 0.0
           else if  (ju>je+myid*je .and. jl>=jb+myid*je) then
             ju = je
             jl = jl - myid*je
             do i = il,iu
               do k = kl,ku
                   if(IIc(i,jl-1,k) == 0) then
                     vp(i, jl, k) = 0.0
                     vm(i, jl, k) = 0.0
                   else
                     vp(i, jl, k) = -vm(i, jl,k)*rk3coefi
                   end if
               end do
             end do
             vp(il:iu, jl + 1:ju, kl:ku) = 0.0 !setting the internal cells
             vm(il:iu, jl + 1:ju, kl:ku) = 0.0
           end if
         end do

         j = jb
         do i = ib,ie
           do k = kb,ke
             if (IIc(i,j,k) == 1 .and. IIc(i,j-1,k)==0 ) then
               vp(i, j, k) = -vm(i, j,k)*rk3coefi
             end if
           end do
         end do

         do n = 1, nxwall
            !kl = block(ixwall(n), 5)
            kl = kb ! tg3315 see comment for x-direction above
            ku = block(ixwall(n), 6) + 1
            il = block(ixwall(n), 1)
            iu = block(ixwall(n), 2)
            jl = max(block(ixwall(n), 3) - myid*jmax, 1)
            ju = min(block(ixwall(n), 4) - myid*jmax, jmax)

            !wp(il:iu, jl:ju, kl:ku) = -wm(il:iu, jl:ju, kl:ku)*rk3coefi
            wp(il:iu, jl:ju, kl) = -wm(il:iu, jl:ju, kl)*rk3coefi
            wp(il:iu, jl:ju, ku) = -wm(il:iu, jl:ju, ku)*rk3coefi
            wp(il:iu, jl:ju, kl + 1:ku - 1) = 0.
            wm(il:iu, jl:ju, kl + 1:ku - 1) = 0.
         end do !1,nxwall

         ! shouldn't be necessary if wp(kb) has not been altered so far,
         ! but just to be safe...
         wm(:, :, kb) = 0.
         w0(:, :, kb) = 0.

         if (ltempeq) then
            thl0volav = sum(thl0av(kb:ke)*dzf(kb:ke))/zh(ke+1)
            do n = 1, nblocks
               il = block(n, 1)
               iu = block(n, 2)
               !kl = block(n, 5)
               kl = kb ! tg3315 see comment for x-direction above
               ku = block(n, 6)
               jl = block(n, 3) - myid*jmax
               ju = block(n, 4) - myid*jmax
               if ((ju < jb) .or. (jl > je)) then
                  cycle
               else
                  if (ju > je) ju = je
                  if (jl < jb) jl = jb
                  thlm(il:iu, jl:ju, kl:ku) = thl0volav
                  thlp(il:iu, jl:ju, kl:ku) = 0.
               end if
            end do
         end if

         if (lmoist) then
            do n = 1, nblocks
               il = block(n, 1)
               iu = block(n, 2)
               !kl = block(n, 5)
               kl = kb ! tg3315 see comment for x-direction above
               ku = block(n, 6)
               jl = block(n, 3) - myid*jmax
               ju = block(n, 4) - myid*jmax
               if ((ju < jb) .or. (jl > je)) then
                  cycle
               else
                  if (ju > je) ju = je
                  if (jl < jb) jl = jb
                  qtm(il:iu, jl:ju, kl:ku) = 0.
                  qtp(il:iu, jl:ju, kl:ku) = 0.
               end if
            end do
         end if

         if (nsv > 0) then
            do n = 1, nblocks
               il = block(n, 1)
               iu = block(n, 2)
               !kl = block(n, 5)
               kl = kb ! tg3315 see comment for x-direction above
               ku = block(n, 6)
               jl = block(n, 3) - myid*jmax
               ju = block(n, 4) - myid*jmax
               if ((ju < jb) .or. (jl > je)) then
                  cycle
               else
                  if (ju > je) ju = je
                  if (jl < jb) jl = jb
                  svp(il:iu, jl:ju, kl:ku, :) = 0.
                  ! attempt to set zero advective flux
                  svm(il, jl:ju, kl:ku, :) = svm(il - 1, jl:ju, kl:ku,:) ! tg3315 swapped these around with jl, ju as was getting values in buildings as blovks are split along x in real topology
                  svm(iu, jl:ju, kl:ku, :) = svm(iu + 1, jl:ju, kl:ku,:)
                  svm(il:iu, jl, kl:ku, :) = svm(il:iu, jl - 1, kl:ku,:)
                  svm(il:iu, ju, kl:ku, :) = svm(il:iu, ju + 1, kl:ku,:)
                  svm(il:iu, jl:ju, ku, :) = svm(il:iu, jl:ju, ku + 1,:)
               end if
            end do
         end if
      end if ! libm

   end subroutine ibmnorm

   !> Determines the distance to the nearest wall for each cell center (used in v. Driest damping function)
   !> Output is a field with minimal distance to wall for each cell center
   !ILS13,10.07.17, not being called anymore
   !only for smagorinsky
   !indeces of walls are wrong (xwallsglobal etc don't exist anymore)
   subroutine nearwall

      use modglobal, only:ib, ie, jb, je, jgb, jge, jmax, kb, ke, xh, xf, dy, zh, zf, lwarmstart, nblocks, libm, lzerogradtop, lwalldist
      use modsubgriddata, only:lsmagorinsky, loneeqn
      use modfields, only:mindist, wall
      use modibmdata, only:xwallsglobal, ywallsglobal, zwallsglobal
      use modmpi, only:myid
      use initfac, only:block
      implicit none

      integer, allocatable :: ux0all(:, :, :), vy0all(:, :, :), wz0all(:, :, :)
      real, allocatable :: distxf(:, :), distxh(:, :), distyf(:, :), distyh(:, :), distzf(:, :), distzh(:, :), &
                           distxf2(:, :), distxh2(:, :), distyf2(:, :), distyh2(:, :), distzf2(:, :), distzh2(:, :), distance(:)
      real distx, disty, distz ! distx is the distance to nearest x-wall, etc.
      ! integer, allocatable :: optie(:)
      integer ic, jc, kc, i, j, k, optie, il, iu, jl, ju, kl, ku, n

      ! if (lwarmstart .or. lles.eqv..false. .or. lvreman) then
      if (((lsmagorinsky) .or. (loneeqn)) .and. (lwalldist)) then
      if (myid == 0) then
         write (6, *) 'Computing wall distances'
      end if
      ! if (lles.eqv..false. .or. lvreman) then
      mindist = 1.0e10

      allocate (ux0all(ib - 1:ie + 1, jgb - 1:jge + 1, kb - 1:ke + 1)) ! This contains ux0 + the lower and (possibly) the upper wall
      allocate (vy0all(ib - 1:ie + 1, jgb - 1:jge + 1, kb - 1:ke + 1)) ! This contains ux0 + the lower and (possibly) the upper wall
      allocate (wz0all(ib - 1:ie + 1, jgb - 1:jge + 1, kb - 1:ke + 1)) ! This contains ux0 + the lower and (possibly) the upper wall
      allocate (distxh(ib:ie, ib:ie + 1))
      allocate (distxf(ib:ie, ib:ie + 1))
      allocate (distyh(jb:je, jgb:jge + 1))
      allocate (distyf(jb:je, jgb:jge + 1))
      allocate (distzh(kb:ke, kb:ke + 1))
      allocate (distzf(kb:ke, kb:ke + 1))
      allocate (distxh2(ib:ie, ib:ie + 1))
      allocate (distxf2(ib:ie, ib:ie + 1))
      allocate (distyh2(jb:je, jgb:jge + 1))
      allocate (distyf2(jb:je, jgb:jge + 1))
      allocate (distzh2(kb:ke, kb:ke + 1))
      allocate (distzf2(kb:ke, kb:ke + 1))
      allocate (distance(4))

      ! initialize wall indicators
      ux0all = 0
      vy0all = 0
      wz0all = 0

      ! Determine for each cell face if an x/y/z-wall is present
      ! from immersed boundaries

      if (libm) then
         ! do loop over blocks
         do n = 1, nblocks
            il = block(n, 1)
            iu = block(n, 2)
            jl = block(n, 3)
            ju = block(n, 4)
            kl = block(n, 5)
            ku = block(n, 6)
            do k = kl, ku
               do j = jl, ju
                  ux0all(il, j, k) = 1 ! lower x-wall
                  ux0all(iu + 1, j, k) = 1 ! upper x-wall
               end do
            end do
            do k = kl, ku
               do i = il, iu
                  vy0all(i, jl, k) = 1 ! lower y-wall
                  vy0all(i, ju + 1, k) = 1 ! upper y-wall
               end do
            end do
            do j = jl, ju
               do i = il, iu
                  wz0all(i, j, kl) = 1 ! lower z-wall
                  wz0all(i, j, ku + 1) = 1 ! upper z-wall
               end do
            end do
         end do ! loop over nblocks

      end if ! libm = .true.

      ! add the global walls (probably upper and lower wall, or only lower wall)
      if (lzerogradtop) then
         do i = ib, ie
         do j = jgb, jge
            wz0all(i, j, kb) = 1 ! ground wall
         end do
         end do
      else
         do i = ib, ie
         do j = jgb, jge
            wz0all(i, j, kb) = 1 ! ground wall
            wz0all(i, j, ke + 1) = 1; ! top wall
         end do
         end do
      end if

      write (6, *) 'Determing distance matrices, proc=', myid
      ! Determine x-distance matrices:
      do ic = ib, ie ! cell-center index
      do i = ib, ie + 1 ! vertex-index (1 more than cell centers)
         distxh(ic, i) = xf(ic) - xh(i)
      end do
      end do

      do ic = ib, ie ! cell-center index
      do i = ib, ie + 1 ! center -index
         distxf(ic, i) = xf(ic) - xf(i)
      end do
      end do

      ! Determine y-distance matrices:
      do jc = jb, je ! cell-center index
      do j = jgb, jge + 1 ! vertex-index (1 more than cell centers) (global index to make sure distance to all cells is determined)
         distyh(jc, j) = (jc + myid*jmax - j)*dy + 0.5*dy
      end do
      end do

      do jc = jb, je ! cell-center index
      do j = jgb, jge + 1 ! center-index  (global index to make sure distance to all cells is determined)
         distyf(jc, j) = (jc + myid*jmax - j)*dy
      end do
      end do

      ! Determine z-distance matrices:
      do kc = kb, ke ! cell-center index
      do k = kb, ke + 1 ! vertex-index (1 more than cell centers)
         distzh(kc, k) = zf(kc) - zh(k)
      end do
      end do

      do kc = kb, ke ! cell-center index
      do k = kb, ke + 1 ! vertex-index (1 more than cell centers)
         distzf(kc, k) = zf(kc) - zf(k)
      end do
      end do

      distxh2 = distxh**2
      distyh2 = distyh**2
      distzh2 = distzh**2
      distxf2 = distxf**2
      distyf2 = distyf**2
      distzf2 = distzf**2

      write (6, *) 'Finished determing distance matrices, proc=', myid
      write (6, *) 'determing distance to nearest wall for each cell center, proc=', myid

      ! Loop over cells (ic,jc,kc) for which minimal wall-to-cell-center-distance needs to be determined
      !  do jc=jgb,jge
      do kc = kb, ke
      do jc = jb, je
      do ic = ib, ie
         ! Determine distance between cc of cell (ic,jc,kc) and faces of all cells (i,j,k)
         do k = kb, ke + 1 ! Level ke+1 is computed in a separate loop (only necessary with upper wall-> global approach=faster)
         do j = jgb, jge + 1 ! loop goes up to jge+1 because jge+1 contains the last vy0-wall
         do i = ib, ie + 1 ! loop goes up to ie+1 because ie+1 contains the last ux0-wall
            if (ux0all(i, j, k) == 1 .OR. vy0all(i, j, k) == 1 .OR. wz0all(i, j, k) == 1) then
               distx = 1.0e10 ! make sure distx is very large when no x-wall is present
               disty = 1.0e10 ! make sure disty is very large when no y-wall is present
               distz = 1.0e10 ! make sure distz is very large when no z-wall is present
               if (ux0all(i, j, k) == 1) then
                  distx = sqrt(distxh2(ic, i) + distyf2(jc, j) + distzf2(kc, k))
               end if
               if (vy0all(i, j, k) == 1) then
                  disty = sqrt(distxf2(ic, i) + distyh2(jc, j) + distzf2(kc, k))
               end if
               if (wz0all(i, j, k) == 1) then
                  distz = sqrt(distxf2(ic, i) + distyf2(jc, j) + distzh2(kc, k))
               end if
            else ! no walls are present in cell (i,j,k) -> distance does not need to be determined for this cell
               cycle ! go to next cell (i,j,k)
            end if
            ! determine minimal wall distance between cc of (ic,jc,kc) and faces of cell (i,j,k)
            distance = (/mindist(ic, jc, kc), distx, disty, distz/)
            optie = minloc(distance, 1)
            ! write(6,*) 'optie=', optie

            if (optie == 1) then
               cycle
            else if (optie == 2) then
               mindist(ic, jc, kc) = distx
               wall(ic, jc, kc, 2) = j
               wall(ic, jc, kc, 3) = k
               ! wall(ic,jc,kc,4) = 1     ! This means the wall closest to the cc of (ic,jc,kc) is at an x-wall at (i,j,k)
               if (ic >= i) then
                  wall(ic, jc, kc, 1) = i
                  wall(ic, jc, kc, 4) = 5 ! shear component index (variable: shear)
                  wall(ic, jc, kc, 5) = 9 ! shear component index (variable: shear)
               else
                  wall(ic, jc, kc, 1) = i - 1 ! in the subgrid this stress is computed in the cell i-1
                  wall(ic, jc, kc, 4) = 6 ! shear component index (variable: shear)
                  wall(ic, jc, kc, 5) = 10 ! shear component index (variable: shear)
               end if
            else if (optie == 3) then
               mindist(ic, jc, kc) = disty
               wall(ic, jc, kc, 1) = i
               wall(ic, jc, kc, 3) = k
               ! wall(ic,jc,kc,4) = 2     ! This means the wall closest to the cc of (ic,jc,kc) is at a y-wall at (i,j,k)
               if (jc + myid*jmax >= j) then
                  wall(ic, jc, kc, 2) = j
                  wall(ic, jc, kc, 4) = 1 ! shear component index (variable: shear)
                  wall(ic, jc, kc, 5) = 11 ! shear component index (variable: shear)
               else
                  wall(ic, jc, kc, 2) = j - 1 ! in the subgrid this stress is computed in the cell j-1
                  wall(ic, jc, kc, 4) = 2 ! shear component index (variable: shear)
                  wall(ic, jc, kc, 5) = 12 ! shear component index (variable: shear)
               end if
            else if (optie == 4) then
               mindist(ic, jc, kc) = distz
               wall(ic, jc, kc, 1) = i
               wall(ic, jc, kc, 2) = j
               ! wall(ic,jc,kc,4) = 3     ! This means the wall closest to the cc of (ic,jc,kc) is at a z-wall at (i,j,k)
               if (kc >= k) then
                  wall(ic, jc, kc, 3) = k
                  wall(ic, jc, kc, 4) = 3 ! shear component index (variable: shear)
                  wall(ic, jc, kc, 5) = 7 ! shear component index (variable: shear)
               else
                  wall(ic, jc, kc, 3) = k - 1 ! in the subgrid this stress is computed in the cel k-1
                  wall(ic, jc, kc, 4) = 4 ! shear component index (variable: shear)
                  wall(ic, jc, kc, 5) = 8 ! shear component index (variable: shear)
               end if
            end if
         ! mindist(ic,jc+myid*jmax,kc)=min(mindist(ic,jc+myid*jmax,kc),distx,disty,distz)   ! global j index
         end do
         end do
         end do
      ! if (myid==0) write(6,*) 'finished for cell (ic,jc,kc)=',ic,jc,kc
      end do
      end do
      end do

      write (6, *) 'Finished determing distance to nearest wall for each cell center, proc=', myid
      ! write(6,*) 'mindist(ib,jb+myid*jmax,kb),mindist(ib,je+myid*jmax,kb)',mindist(ib,jb+myid*jmax,kb),mindist(ib,je+myid*jmax,kb)

      else
      return
      end if !(lwarmstart)

      deallocate (ux0all, vy0all, wz0all)
      deallocate (xwallsglobal, ywallsglobal, zwallsglobal, block) ! used for determining boundaries
      return
   end subroutine nearwall

   subroutine bottom
      !kind of obsolete when road facets are being used
      !vegetated floor not added (could simply be copied from vegetated horizontal facets)
      use modglobal, only:ib, ie, ih, jh, kh, jb, je, kb, numol, prandtlmol, dzh, nsv, &
         dxf, dxhi, dzf, dzfi, numoli, ltempeq, khc, lmoist, BCbotT, BCbotq, BCbotm, BCbots, dzh2i, libm
      use modfields, only : u0,v0,e120,um,vm,w0,wm,e12m,thl0,qt0,sv0,thlm,qtm,svm,up,vp,thlp,qtp,svp,shear,momfluxb,tfluxb,cth
      use modsurfdata, only:thlflux, qtflux, svflux, ustar, thvs, wtsurf, wqsurf, thls, z0, z0h
      use modsubgriddata, only:ekm, ekh
      use modmpi, only:myid
      implicit none
      integer :: i, j, jp, jm, m
      !momentum

      if (.not.(libm)) then
         if (BCbotm.eq.1) then

         elseif (BCbotm.eq.2) then
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

      return
   end subroutine bottom


     REAL FUNCTION unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !for momentum, this bit is not depended on orientation etc
        !momentum flux in m2/s2
        !dT,utang and logdzh are unused and could be removed
           IMPLICIT NONE
           REAL, INTENT(in) :: logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2
           REAL :: Ribl1, Fm, Fh, cm, ch, Ctm, M
           REAL, PARAMETER :: b1 = 9.4 !parameters from Uno1995
           REAL, PARAMETER :: b2 = 4.7
           REAL, PARAMETER :: dm = 7.4
           REAL, PARAMETER :: dh = 5.3
           REAL, PARAMETER :: prandtlmol = 0.71
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

end module modibm
