!!> \file b2f.f90
!!  \author Ivo Suter
!
!   reads the necessary input files to deal with facets
!
!   WARNING: if walls with more than 3 layers (4 points) are to be considered, this file needs to be changed
!            e.g. walltypes needs to read 7+4*nlayers columns, offsets in reading facet properties also change accordingly
   module initfac
      use modglobal, only:ifinput, nblocks, nfcts, cexpnr, libm, bldT, rsmin, wsoil, wfc,&
                          nwalllayers, block,lEB
      use modmpi, only:myid, comm3d, mpierr, MPI_INTEGER, MPI_DOUBLE_PRECISION, MY_REAL, nprocs, cmyid, MPI_REAL8, MPI_REAL4, MPI_SUM, mpi_logical
      use netcdf
      implicit none
      public :: readfacetfiles,qsat,dqsatdT
      save

      !integer, allocatable :: block(:, :) !block coordinates and facet Nr corresponding to block faces
      !facet properties
      logical, allocatable :: faclGR(:) !logic array, is it a green (vegetated) facet?
      real, allocatable    :: facz0(:) !roughness for momentum on facets
      real, allocatable    :: facz0h(:) !roughness for heat and moisture on facets
      real, allocatable    :: facalb(:) !facet shortwave albedo
      real, allocatable    :: facem(:) !facet longwave emissivity of all 5 faces of the blocks
      real, allocatable    :: facdi(:, :) !inverse facet thickness
      real, allocatable    :: facd(:,:) !facet thickness
      real, allocatable    :: faccp(:, :) !facet specific heat capacity
      real, allocatable    :: faclami(:, :) !inverse facet heat conductivity
      real, allocatable    :: fackappa(:, :) !facet heat diffusivity (lambda/(rho cp))
      real, allocatable    :: faca(:) !facet area
      integer, allocatable :: facain(:) !facet area as sum of indeces
      integer, allocatable :: facets(:, :) !facet orientation, walltype, block and building Nr
      real, allocatable    :: walltypes(:, :) !the defined wall and rooftypes with properties
      !radiation
      real, allocatable    :: vf(:, :) !viewfactors between facets
      real, allocatable    :: svf(:) !sky-viewfactor of facets
      real, allocatable    :: netsw(:) !net shortwave radiation on facets
      real, allocatable    :: facLWin(:) !incoming longwave on facets [W/m2]
      !temperature
      real, allocatable    :: Tfacinit(:) !initial facet temperatures
      real, allocatable    :: facT(:, :) !wall temperatures on surfaces and between layers (1=outdoors,end=indoors)
      real, allocatable    :: facTdash(:, :)!temperature gradient dT/dz
      !fluxes
      real, allocatable    :: facef(:) !evaporative flux on facets [W/m2] (single processor)
      real, allocatable    :: facefi(:) !time integrated latent heat flux [J/m2] (used in modEB)
      real, allocatable    :: facefsum(:) !evaporative flux on facets [W/m2] (sum over all processors)
      real, allocatable    :: fachf(:) !heat flux on facets [Km/s] (single processor)
      real, allocatable    :: fachfi(:) !time integrated heat flux [Km]
      real, allocatable    :: fachfsum(:) !heat flux on facets [Km/s] (sum over all processors)
      !GR
      real, allocatable    :: facf(:, :) !dependence of stomatal resistance/soil resistance
      real, allocatable    :: fachurel(:) !relative humidity at ground surface
      real, allocatable    :: facwsoil(:) !soil moisture of facets
      real, allocatable    :: faccth(:) !sum of all transfer coefficients of the facet, used in Penman Moneith, unused
      real, allocatable    :: facqsat(:) !saturation absoulute humidity at facet temperature

      !misc
      integer, allocatable :: typeloc(:) !array to match the walltype to sequential integers for indexing
      integer              :: nwalltypes = 0 !number of different walltypes, will be determined automatically
      character(80)        :: chmess !dummy character string
      integer              :: nwallprops

   contains

      subroutine readfacetfiles
      use modglobal, only: block, cexpnr, iwalltemp
      implicit none

      !use modglobal, only:block
      !read initial&unchangeable facet values from files
      !read blocks and facet Nr corresponding to block faces (Order: Top, West, East, North, South)
      !define facets with properties and initial temperature
      !read facets.inp.xxx facetarea.inp.xxx vf.inp.xxx walltypes.inp.xxx
      !read netsw.inp.xxx (if sun is not constant, K needs to be calculated at every EB-timestep)
      !read tfacinit.inp.xxx
      !use modglobal, only : nblocks, nfcts, cexpnr, ifinput
      character (len = 13) :: FILE_VF = 'vf.nc.inp.xxx'
      integer :: ncid, varid
      integer :: n = 0, m = 0, i = 0, j = 0, k = 0, io = 0
      integer :: iret

      if (.not.(nfcts > 0)) return

      nwallprops = 6 + 4*nwalllayers+1
      !allocate (block(nblocks, 11))
      allocate (faclGR(0:nfcts))
      allocate (facz0(0:nfcts)) !0 is the default value (e.g. for internal walls)
      allocate (facz0h(0:nfcts))
      allocate (facalb(0:nfcts))
      allocate (facem(0:nfcts))
      allocate (facdi(0:nfcts, nwalllayers))
      allocate (facd(0:nfcts, nwalllayers))
      allocate (faccp(0:nfcts, nwalllayers))
      allocate (faclami(0:nfcts, nwalllayers))
      allocate (fackappa(0:nfcts, nwalllayers+1))
      allocate (faca(0:nfcts))
      allocate (facain(0:nfcts))
      allocate (facets(nfcts, 4))

      if (lEB .eqv. .true.) then
        allocate (vf(1:nfcts, 1:nfcts))
        allocate (svf(1:nfcts))
        allocate (netsw(1:nfcts))
        allocate (facLWin(1:nfcts))
      end if

      allocate (Tfacinit(1:nfcts))
      allocate (facT(0:nfcts, nwalllayers+1))
      allocate (facTdash(1:nfcts,nwalllayers+1))
      allocate (facef(1:nfcts))
      allocate (facefi(1:nfcts))
      allocate (facefsum(1:nfcts))
      allocate (fachf(0:nfcts))
      allocate (fachfi(0:nfcts))
      allocate (fachfsum(1:nfcts))
      allocate (facf(0:nfcts, 5))
      allocate (fachurel(0:nfcts))
      allocate (facwsoil(0:nfcts))
      allocate (faccth(0:nfcts))
      allocate (facqsat(0:nfcts))

      !block = 0;
      faclGR = .false.; facz0 = 0.; facz0h = 0.; facalb = 0.; facem = 0.; facd=0.; facdi = 0.; faccp = 0.
      faclami = 0.; fackappa = 0.; faca = 0.; facain = 0; facets = 0
      if (lEB .eqv. .true.) then
        vf = 0.; svf = 0.; netsw = 0.; facLWin = 0.
      end if
      Tfacinit = 0.; facT = 0.; facTdash = 0.
      facef = 0.; facefi = 0.; facefsum = 0.; fachf = 0.; fachfi = 0.; fachfsum = 0.
      facf = 0.; fachurel = 0.; facwsoil = 0.; faccth = 0.; facqsat = 0.;

      if (myid == 0 .and. libm) then ! read blocks corner coordinates and facet-Nr correspoinding to the block sides
      !   open (ifinput, file='blocks.inp.'//cexpnr)
      !   read (ifinput, '(a80)') chmess
      !   read (ifinput, '(a80)') chmess
      !   do n = 1, nblocks
      !     blockfile is :ilow iheigh jlow jheigh klow kheigh facTop facWest facEast facNorth facSouth
      !     read (ifinput, *) &
      !         block(n, 1), &
      !         block(n, 2), &
      !         block(n, 3), &
      !         block(n, 4), &
      !         block(n, 5), &
      !         block(n, 6), &
      !         block(n, 7), &
      !         block(n, 8), &
      !         block(n, 9), &
      !         block(n, 10), &
      !         block(n, 11)
      !   end do
      !   close (ifinput)

        if (lEB) then
        ! read facet areas
          open (ifinput, file='facetarea.inp.'//cexpnr)
          read (ifinput, '(a80)') chmess
          do n = 1, nfcts
            read (ifinput, *) &
              faca(n)
          end do
          close (ifinput)
          write (*, *) "faca", faca
        end if !lEB

        ! read wall & (green) roof types
        ! read once to determine number of types, allocate, read again
        nwalltypes = -3 !3 lines as headers
        open (ifinput, file='walltypes.inp.'//cexpnr)
        do
          read (ifinput, *, iostat=io)
          if (io /= 0) exit
          nwalltypes = nwalltypes + 1
        end do
        close (ifinput)
      end if !(myid == 0 .and. libm)

      call MPI_BCAST(nwalltypes, 1, MPI_Integer, 0, comm3d, mpierr)
      allocate (walltypes(1:nwalltypes, nwallprops))

      if (myid == 0) then
        walltypes = 0.
        open (ifinput, file='walltypes.inp.'//cexpnr)
        read (ifinput, '(a80)') chmess
        read (ifinput, '(a80)') chmess
        read (ifinput, '(a80)') chmess

        do n = 1, nwalltypes
          read (ifinput, *) (walltypes(n, m), m=1,nwallprops)
        end do
        close (ifinput)
      end if !myid==0

        call MPI_BCAST(nwalltypes, 1, MPI_Integer, 0, comm3d, mpierr)
        call MPI_BCAST(walltypes, nwallprops*nwalltypes, MY_REAL, 0, comm3d, mpierr)

        !create an array mapping walltypes to sequential integers for indexing
        !e.g. lets assume walltype -3,-1,1,2,3 and 5 are defined.
        !index: [-3,-2,-1,0,1,2,3,4,5]  -> [-3,-2,-1,0,1,2,3,4,5]
        !value: [ 0, 0, 0,0,0,0,0,0,0]  -> [ 1, 0, 2,0,3,4,5,0,6]
        allocate (typeloc(int(minval(walltypes(:, 1))):int(maxval(walltypes(:, 1)))))

        if (myid .eq. 0) then  !all the read processes just need to be done on one processor
          typeloc = 0
          do n = 1, nwalltypes
            typeloc(int(walltypes(n, 1))) = n
          end do

          ! read the facets
          open (ifinput, file='facets.inp.'//cexpnr)
          read (ifinput, '(a80)') chmess
          !orientation, walltype, block in block.inp, original block/building
          do n = 1, nfcts
            read (ifinput, *) &
              facets(n, 1), &
              facets(n, 2), &
              facets(n, 3), &
              facets(n, 4)
          end do
          close (ifinput)

          ! calculate the number of indeces per facet in both dimensions (~area)
          do n = 1, nfcts
            if (facets(n, 1) .eq. 1) then
              facain(n) = (block(facets(n, 3), 2) - block(facets(n, 3), 1) + 1)*(block(facets(n, 3), 4) - block(facets(n, 3), 3) + 1)
            else if (facets(n, 1) .eq. 2) then
              facain(n) = (block(facets(n, 3), 6) - block(facets(n, 3), 5) + 1)*(block(facets(n, 3), 4) - block(facets(n, 3), 3) + 1)
            else if (facets(n, 1) .eq. 3) then
              facain(n) = (block(facets(n, 3), 6) - block(facets(n, 3), 5) + 1)*(block(facets(n, 3), 4) - block(facets(n, 3), 3) + 1)
            else if (facets(n, 1) .eq. 4) then
              facain(n) = (block(facets(n, 3), 2) - block(facets(n, 3), 1) + 1)*(block(facets(n, 3), 6) - block(facets(n, 3), 5) + 1)
            else if (facets(n, 1) .eq. 5) then
              facain(n) = (block(facets(n, 3), 2) - block(facets(n, 3), 1) + 1)*(block(facets(n, 3), 6) - block(facets(n, 3), 5) + 1)
            end if
          end do

          ! assign the facet properties to their own arrays
          do n = 1, nfcts
            i = typeloc(facets(n, 2))
            faclGR(n) = (abs(walltypes(i, 2) - 1.00) < 1.0D-5) !logic for green surface, conversion from real to logical
            facz0(n) = walltypes(i, 3)  !surface momentum roughness
            facz0h(n) = walltypes(i, 4) !surface heat & moisture roughness
            facalb(n) = walltypes(i, 5) !surface shortwave albedo
            facem(n) = walltypes(i, 6)  !surface longwave emissivity

            if (facets(n, 2) .lt. -100) then !it's a bounding wall
              do j = 1, nwalllayers !bounding walls don't need properties
                facdi(n, j) = 0.   !bounding walls have no energy balance
                facd(n,j)=0.       !bounding walls have no thickness
                faclami(n, j) = 0. !bounding walls have no energy balance
                faccp(n, j) = 0.   !bounding walls have no heat capacity
              end do
            else
              do j = 1, nwalllayers !for all layers
                facdi(n, j) = 1/walltypes(i, j + 6) !inverse of facet thickness of layer j
                facd(n,j) = walltypes(i,j+6) !facet thickness of layer j
                faclami(n, j) = 1/walltypes(i, j+12) !inverse of heat conductivity of layer j
                faccp(n, j) = walltypes(i, j+9) !specific heat capacity of layer j
              end do
            end if

            do j= 1,nwalllayers+1
              fackappa(n, j) = walltypes(i, j+15) !heat duffusivity of layer 1
            end do
          end do

          !give some dummy values for block internal facets (i.e. facets with walltype 0)
          facz0(0) = 0.00999; facz0h(0) = 0.00999; facalb(0) = 0.999; facem(0) = 0.999; facd(0,1)=0.999; facd(0,2)=0.999; facdi(0, 1) = 0.999; facdi(0, 2) = 0.999; facdi(0, 3) = 0.999; faccp(0, 1) = 999.; faccp(0, 2) = 999.
          faccp(0, 3) = 999.; faclami(0, 1) = 0.999; faclami(0, 2) = 0.999; faclami(0, 3) = 0.999; fackappa(0, 1) = 0.00000999; fackappa(0, 2) = 0.00000999; fackappa(0, 3) = 0.00000999; faclGR(0) = .false.

          if (lEB) then
            ! read viewfactors between facets
            ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
            ! the file.
            FILE_VF = 'vf.nc.inp.'//cexpnr
            iret = nf90_open(FILE_VF, NF90_NOWRITE, ncid) ! Get the varid of the data variable, based on its name.
            iret = nf90_inq_varid(ncid, "view factor", varid)
            ! Read the data.
            iret = nf90_get_var(ncid, varid, vf)
            write(*,*) "vf(1,6),vf(1,7),vf(6,1),vf(7,1)",vf(1,6),vf(1,7),vf(6,1),vf(7,1)

            ! read skyviewfactors
            open (ifinput, file='svf.inp.'//cexpnr)
            read (ifinput, '(a80)') chmess
            do n = 1, nfcts
              read (ifinput, *) &
                svf(n)
              end do
              close (ifinput)

              ! read net shortwave radiation
              open (ifinput, file='netsw.inp.'//cexpnr)
              read (ifinput, '(a80)') chmess
              do n = 1, nfcts
                read (ifinput, *) &
                  netsw(n)
              end do
              close (ifinput)
            end if !lEB

            if ((lEB) .or. (iwalltemp == 2)) then
              ! read initial facet temepratures
              open (ifinput, file='Tfacinit.inp.'//cexpnr)
              read (ifinput, '(a80)') chmess
              do n = 1, nfcts
                read (ifinput, *) &
                  Tfacinit(n)
              end do
              close (ifinput)

              do n = 1, nfcts
                facT(n, 1) = Tfacinit(n) !building surfaces is given an initial temperature
                facT(n, nwalllayers+1) = bldT !inner most layer has the same temperature as the building interior
                do j = 2,nwalllayers
                  facT(n, j) = Tfacinit(n)-(Tfacinit(n)-bldT)/nwalllayers*(j-1) !scale linearly inside the wall
                end do
              end do

              do n = 1,nwalllayers
                 facT(0, n) = 288.
              end do
              facT(0, nwalllayers+1) = 299.
            end if !((lEB) .or. (iwalltemp == 2))

            ! assign initial soil moisture (only outermost layer)
            do n = 1, nfcts
              if (faclGR(n)) then
                facwsoil(n) = wsoil
                fachurel(n) = 0.5*(1. - cos(3.14159*wsoil/wfc))
              end if
            end do

         end if !(myid .eq. 0)

         !write (*, *) "starting broadcast of facet properties"
         call MPI_BCAST(block, 11*nblocks, MPI_INTEGER, 0, comm3d, mpierr)
         !many of these are actually only needed on processor 0...
         call MPI_BCAST(faclGR(0:nfcts), nfcts + 1, mpi_logical, 0, comm3d, mpierr)
         call MPI_BCAST(facz0(0:nfcts), nfcts + 1, MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(facz0h(0:nfcts), nfcts + 1, MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(facalb(0:nfcts), nfcts + 1, MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(facem(0:nfcts), nfcts + 1, MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(facd(0:nfcts,1:nwalllayers),(nfcts+1)*nwalllayers, MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(facdi(0:nfcts, 1:nwalllayers), (nfcts + 1)*nwalllayers, MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(faccp(0:nfcts, 1:nwalllayers), (nfcts + 1)*nwalllayers, MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(faclami(0:nfcts, 1:nwalllayers), (nfcts + 1)*nwalllayers, MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(fackappa(0:nfcts, 1:nwalllayers+1), (nfcts + 1)*(nwalllayers+1), MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(faca(0:nfcts), nfcts + 1, MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(facain(0:nfcts), nfcts + 1, MPI_Integer, 0, comm3d, mpierr)
         call MPI_BCAST(facets, 4*nfcts, MPI_Integer, 0, comm3d, mpierr)
         !walltypes is broadcast further up
         if (lEB .eqv. .true.) then
           call MPI_BCAST(svf(1:nfcts), nfcts, MY_REAL, 0, comm3d, mpierr)
           call MPI_BCAST(netsw(1:nfcts), nfcts, MY_REAL, 0, comm3d, mpierr)
         end if
         !facLWin currently not being broadcast..
         !vf currently not being broadcast
         call MPI_BCAST(Tfacinit(1:nfcts), nfcts, MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(facT(0:nfcts, 1:nwalllayers+1), (nfcts + 1)*(nwalllayers+1), MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(facTdash(1:nfcts, 1:nwalllayers+1), (nfcts)*(nwalllayers+1), MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(facef(1:nfcts), nfcts, MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(facefi(1:nfcts), nfcts, MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(facefsum(1:nfcts), nfcts, MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(fachf(0:nfcts), nfcts + 1, MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(fachfi(0:nfcts), nfcts + 1, MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(fachfsum(1:nfcts), nfcts, MY_REAL, 0, comm3d, mpierr)
         ! standard plant & soil resistance for grass (Manickathan2018) in s/m
         facf(:,4)=200.
         facf(:,5)=50.
         call MPI_BCAST(facf(0:nfcts, 1:5), (nfcts + 1)*5, MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(fachurel(0:nfcts), nfcts + 1, MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(facwsoil(0:nfcts), nfcts + 1, MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(faccth(0:nfcts), nfcts + 1, MY_REAL, 0, comm3d, mpierr)
         do n = 1, nfcts
         facqsat(n) = qsat(facT(n,1))
         end do
         call MPI_BCAST(facqsat(0:nfcts), nfcts + 1, MY_REAL, 0, comm3d, mpierr)


      end subroutine readfacetfiles

      real function qsat(T)
         implicit none
         real, intent(in) :: T
         real :: gres
         gres = 611.00*exp(17.27*(T - 273.15)/(T - 35.85)) ![Pa] Bolton 1980
         qsat = 0.62198*0.01*gres/(1000-0.01*gres) ![kg/kg] Murphy & Koop 2005 !1000 can be replaced with actual air pressure if desired
      end function qsat

      real function dqsatdT(T)
         implicit none
         real, intent(in) :: T
         dqsatdT = 0.1384832710e-2 + 0.7708409674e-4*(T - 300) + 0.2022064593e-5*(T - 300)**2 + 0.000000036561*(T - 300)**3 !expansion of qsat(T)
      end function dqsatdt

   end module initfac
