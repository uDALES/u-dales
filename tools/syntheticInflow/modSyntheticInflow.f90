!! This program generates the driver files for u-DALES target simulations.
!!
!! Main inputs :
!!
!!   1) namoptions to read: iexpnr, dtmax,runtime, nprocy, itot,jtot,ktot, xlen,ylen, ltempeq,lmoist
!!   2) prof.inp. to read z info
!!
!!   Note all the below variables should be at cell edges starting from z = 0 to z = zsize. 
!!   Hence, length(ut) is length(zh)+1. This is true for all the input files (3) to (18).
!!   Also, the mean shoud be computed along t, x, and y
!!
!!   3) length_time_scales_u.txt --> to read length and time scales of u as function of z.
!!   4) length_time_scales_v.txt --> to read length and time scales of u as function of z.
!!   5) length_time_scales_w.txt --> to read length and time scales of u as function of z.
!!
!!   6) ut.txt    --> to read mean u0 as function of z only.
!!   7) upup.txt  --> to read mean u'u' as function of z.
!!   8) upvp.txt  --> to read mean u'v' as function of z.
!!   9) upwp.txt  --> to read mean u'w' as function of z.
!!  10) vpvp.txt  --> to read mean v'v' as function of z.
!!  11) vpwp.txt  --> to read mean v'w' as function of z.
!!  12) wpwp.txt  --> to read mean w'w' as function of z.
!!
!!  IF (ltempeq .eqv. .TRUE.) THEN
!!      13) length_time_scales_temp.txt --> to read length and time scales of u as function of z.
!!      14) thlt.txt  --> to read mean temperature (thl0) as function of z.
!!      15) thlpthlpt.txt  --> to read mean thlpthlp as function of z.
!!      16) wpthlpt.txt  --> to read mean w'thlp as function of z.
!!  END IF
!!
!!  IF (lmoist .eqv. .TRUE.) THEN
!!      17) length_time_scales_qt.txt --> to read length and time scales of u as function of z.
!!      18) qtt.txt  --> to read mean moisture (qt0) as function of z.
!!  END IF

!! For creating the files (3) to (18), one can look at the 'write_Reynolds_stress.m' file inside the
!! 'u-dales/tools/syntheticInflow/' directory.


PROGRAM synInflowGen
    
    USE OMP_LIB

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(1.0d0)
    REAL(KIND=dp), PARAMETER :: pi = 4.0d0*ATAN(1.0d0)

    !! ################################################################################

    INTEGER :: iexpnr
    CHARACTER(3) :: cexpnr
    INTEGER :: ny, nz, nt
    INTEGER :: itot, jtot, ktot
    INTEGER :: nprocy
    REAL(KIND=dp) :: dt_sig
    REAL(KIND=dp) :: dx, dy, ddy
    REAL(KIND=dp) :: xlen, ylen, runtime, dtmax
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: yf, yh, zf, zh, ydriver, zdriver, dz, ddz
    REAL(KIND=dp) :: u_bulk, ztotal
    LOGICAL :: lcalc_time_and_length_scale = .FALSE.
    LOGICAL :: lcheck_driver_outputs = .TRUE.

    !! ################################################################################

    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: umean
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: R11, R21, R22, R31, R32, R33
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: a11, a21, a22, a31, a32, a33
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: tu_scale, tv_scale, tw_scale

    INTEGER, DIMENSION(:), ALLOCATABLE :: nluy, nluz, nlvy, nlvz, nlwy, nlwz
    INTEGER, DIMENSION(:), ALLOCATABLE :: NUY,  NUZ,  NVY,  NVZ,  NWY,  NWZ
    INTEGER, DIMENSION(:), ALLOCATABLE :: NUY2, NUZ2, NVY2, NVZ2, NWY2, NWZ2

    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: buy, buz, bvy, bvz, bwy, bwz
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: psiu, psiv, psiw
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: psiu_new, psiv_new, psiw_new

    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: udriver, vdriver, wdriver, u0, v0, w0

    !! ################################################################################

    LOGICAL :: ltempeq = .FALSE.
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: tempmean
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: R41, R42, R43, R44
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: a41, a42, a43, a44
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: ttemp_scale
    INTEGER, DIMENSION(:), ALLOCATABLE :: nltempy, nltempz
    INTEGER, DIMENSION(:), ALLOCATABLE :: NtempY,  NtempZ
    INTEGER, DIMENSION(:), ALLOCATABLE :: NtempY2, NtempZ2
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: btempy, btempz
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: psitemp
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: psitemp_new
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: tempdriver, thl0

    !! ################################################################################
    
    LOGICAL :: lmoist = .FALSE.
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: qtmean
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: R51, R52, R53, R54, R55
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: a51, a52, a53, a54, a55
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: tqt_scale
    INTEGER, DIMENSION(:), ALLOCATABLE :: nlqty, nlqtz
    INTEGER, DIMENSION(:), ALLOCATABLE :: NqtY,  NqtZ
    INTEGER, DIMENSION(:), ALLOCATABLE :: NqtY2, NqtZ2
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: bqty, bqtz
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: psiqt
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: psiqt_new
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: qtdriver, qt0

    !! ################################################################################

    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    !! ################################################################################

    !$ DOUBLE PRECISION :: start_main, end_main

    !! ################################################################################

    CALL read_namelist 
    
    CALL init_synInflowGen
    
    !$ start_main = OMP_GET_WTIME()
    CALL synInflowGen_main
    !$ end_main = OMP_GET_WTIME()
    !$ WRITE(*,*) "Main computation took ", end_main - start_main, " seconds"

    CALL finish_synInflowGen


    CONTAINS

    !! The below subroutine "read_namelist" reads the 'namoptions' file and computes the necessary inputs.
    !! Main important entries from 'namoptions' that are necessary for this program are:
    !! iexpnr  --> experiment number as 3 digits integer. This is the experiment number of the target simulation.
    !!             This is necessary for creating output file names.
    !! dtmax   --> time-step size for the driver data generation. There for be careful about this while creating namoptions file.
    !!             Do no put a high dtmax in namoptions file. In that case target simulation will need compute several 
    !!             time interpolations, and the accuracy of the simulation will be bad.
    !! runtime --> time length of the driver files. 
    !!             This program will generate driver data for t = 0.0 to t = runtime in stepse of dtmax
    !! nprocy  --> number of processor along y-direction to be used in the target simulation
    !! itot    --> number of division/grid cells along x-direction
    !! jtot    --> number of division/grid cells along y-direction
    !! ktot    --> number of division/grid cells along z-direction
    !! xlen    --> domain size along x-direction
    !! ylen    --> domain size along y-direction
    !! ltempeq --> logical switch to decide whether to compute/write driver files for temperature
    !! lmoist  --> logical switch to decide whether to compute/write driver files for moisture
    SUBROUTINE read_namelist
        IMPLICIT NONE
        
        CHARACTER(14) :: fname_options = 'namoptions.'
        CHARACTER(90) :: chmess
        INTEGER, PARAMETER :: ifinput = 1, ifnamopt = 3
        INTEGER :: ierr, j, k

        !! ################################################################################
        LOGICAL :: lwarmstart, lstratstart, ladaptive, libm, lles, lper2inout, lwalldist, lreadmean, lrandomize, &
                   lcoriol, lbuoyancy, lprofforc, lvinf, luoutflowr, lvoutflowr, luvolflowr, lvvolflowr, lnudge, &
                   ltimedepsurf, ltimedepnudge, ltimedeplw, ltimedepsw
        CHARACTER(90) :: startfile, author
        REAL :: trestart, randu, randthl, randqt, courant, diffnr, xlat, xlon, xday, xtime, ps, tscale, dpdx, &
                uflowrate, vflowrate, tnudge
        INTEGER :: irandom, krand, nprocx, ksp, igrw_damp, ifixuinf, nnudge, ntimedepsurf, ntimedepnudge, &
                   ntimedeplw, ntimedepsw
        !! ################################################################################
        
        NAMELIST /RUN/ &
            iexpnr, lwarmstart, lstratstart, startfile, &
            runtime, dtmax, trestart, ladaptive, &
            irandom, randu, randthl, randqt, krand, &
            courant, diffnr, author, &
            libm, lles, &
            lper2inout, lwalldist, &
            lreadmean, &
            nprocx, nprocy, &
            lrandomize

        NAMELIST /DOMAIN/ &
            itot, jtot, ktot, xlen, ylen, &
            xlat, xlon, xday, xtime, ksp

        NAMELIST /PHYSICS/ &
            ps, igrw_damp, lmoist, lcoriol, lbuoyancy, ltempeq, &
            lprofforc, ifixuinf, lvinf, tscale, dpdx, &
            luoutflowr, lvoutflowr, luvolflowr, lvvolflowr, &
            uflowrate, vflowrate, &
            lnudge, tnudge, nnudge, &
            ltimedepsurf, ntimedepsurf, ltimedepnudge, ntimedepnudge, &
            ltimedeplw, ntimedeplw, ltimedepsw, ntimedepsw

        IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
            CALL GET_COMMAND_ARGUMENT(1, fname_options(12:14))
        END IF

        OPEN(ifnamopt,FILE=fname_options,STATUS='OLD',IOSTAT=ierr)
            IF (ierr /= 0) THEN
                WRITE(0, *) 'ERROR: Namoptions does not exist'
                WRITE(0, *) 'iostat error: ', ierr
                STOP 1
            END IF

            READ(ifnamopt, DOMAIN, IOSTAT=ierr)
            IF (ierr > 0) THEN
                WRITE(0, *) 'ERROR: Problem in namoptions DOMAIN'
                WRITE(0, *) 'iostat error: ', ierr
                STOP 1
            END IF
            ! WRITE(6, DOMAIN)
            REWIND(ifnamopt)

            READ(ifnamopt, RUN, IOSTAT=ierr)
            IF (ierr > 0) THEN
                WRITE(0, *) 'ERROR: Problem in namoptions RUN'
                WRITE(0, *) 'iostat error: ', ierr
                STOP 1
            END IF
            ! WRITE(6, RUN)
            REWIND(ifnamopt)

            READ(ifnamopt, PHYSICS, iostat=ierr)
            IF (ierr > 0) THEN
                WRITE(0, *) 'ERROR: Problem in namoptions PHYSICS'
                WRITE(0, *) 'iostat error: ', ierr
                STOP 1
            END IF
            ! WRITE(6, PHYSICS)
            REWIND(ifnamopt)

        CLOSE(ifnamopt)
        !! ########## namoptions reading is done here. ###########

        WRITE(cexpnr,'(i3.3)') iexpnr
        ny = jtot
        nz = ktot
        nt = CEILING(runtime/dtmax)
        dt_sig = dtmax
        dx = xlen/REAL(itot, KIND=dp)
        dy = ylen/REAL(jtot, KIND=dp)
        ddy = 1.0d0/dy

        ALLOCATE( yf(1:jtot), yh(1:jtot), zf(1:ktot), zh(1:ktot), ydriver(0:ny), zdriver(0:nz) )
        ALLOCATE( dz(0:nz), ddz(0:nz) )
    
        yf(1:jtot) = [((j-1)*dy+(dy/2.0d0), j=1,jtot)]
        yh(1:jtot) = [((j-1)*dy, j=1,jtot)]
        
        OPEN(ifinput, FILE='prof.inp.'//cexpnr,STATUS='OLD',IOSTAT=ierr)
        IF (ierr /= 0) THEN
            WRITE(0, *) 'ERROR: prof.inp.', cexpnr, ' does not exist'
            STOP 1
        END IF
        READ(ifinput, '(A72)') chmess
        READ(ifinput, '(A72)') chmess
        DO k = 1,ktot
            READ(ifinput, *) zf(k)
        END DO
        CLOSE(ifinput)

        DO k = 1,ktot
            IF (k==1) THEN
                zh(k) = 0.0d0;
            ELSE
                zh(k) = (zf(k)+zf(k-1))/2.0d0
            END IF
        END DO

        ydriver(0:ny-1) = yh(1:jtot)
        ydriver(ny) = yf(jtot) + (yf(jtot)-yf(jtot-1))/2.0d0

        zdriver(0:nz-1) = zh(1:ktot)
        zdriver(nz) = zf(ktot) + (zf(ktot)-zf(ktot-1))/2.0d0

        DO k = 0,nz
            IF (k==nz) THEN
                dz(k) = zdriver(k)-zdriver(k-1)
            ELSE
                dz(k) = zdriver(k+1)-zdriver(k)
            END IF
        END DO

        ddz = 1.0d0/dz

        ! DO k = 0,nz
        !     WRITE(*,*) k , " ", zf(k), " ", zh(k), " ", zdriver(k), " ", dz(k)
        ! END DO

        IF (lcheck_driver_outputs) THEN
            OPEN(UNIT=20,FILE='z_driver.txt')
            DO k = 0,nz
                WRITE(UNIT=20,FMT='(f15.10)',ADVANCE='NO') zdriver(k)
                WRITE(UNIT=20,FMT='(A)',ADVANCE='NO') NEW_LINE('a')
            END DO
            CLOSE(UNIT=20)
        END IF

    END SUBROUTINE read_namelist


    SUBROUTINE init_synInflowGen
        IMPLICIT NONE

        INTEGER, DIMENSION(8) :: current_time
        INTEGER(kind=8) :: time_now
        INTEGER :: seed_size
        INTEGER :: k

        ALLOCATE( umean(0:nz), R11(0:nz), R21(0:nz), R22(0:nz), R31(0:nz), R32(0:nz), R33(0:nz) )
        ALLOCATE( a11(0:nz), a21(0:nz), a22(0:nz), a31(0:nz), a32(0:nz), a33(0:nz) )
        ALLOCATE( tu_scale(0:nz), tv_scale(0:nz), tw_scale(0:nz) )

        ALLOCATE( nluy(0:nz), nluz(0:nz), nlvy(0:nz), nlvz(0:nz), nlwy(0:nz), nlwz(0:nz) )
        ALLOCATE( NUY(0:nz),  NUZ(0:nz),  NVY(0:nz),  NVZ(0:nz),  NWY(0:nz),  NWZ(0:nz) )
        ALLOCATE( NUY2(0:nz), NUZ2(0:nz), NVY2(0:nz), NVZ2(0:nz), NWY2(0:nz), NWZ2(0:nz) )

        ALLOCATE( psiu(0:ny,0:nz), psiv(0:ny,0:nz), psiw(0:ny,0:nz) )
        ALLOCATE( psiu_new(0:ny,0:nz), psiv_new(0:ny,0:nz), psiw_new(0:ny,0:nz) )

        ALLOCATE( udriver(0:ny,0:nz), vdriver(0:ny,0:nz), wdriver(0:ny,0:nz) )
        ALLOCATE( u0(1:jtot,1:ktot), v0(1:jtot,1:ktot), w0(1:jtot,1:ktot) )
        
        IF (ltempeq) THEN
            ALLOCATE( tempmean(0:nz), R41(0:nz), R42(0:nz), R43(0:nz), R44(0:nz) )
            ALLOCATE( a41(0:nz), a42(0:nz), a43(0:nz), a44(0:nz) )
            ALLOCATE( ttemp_scale(0:nz) )
            ALLOCATE( nltempy(0:nz), nltempz(0:nz) )
            ALLOCATE( NtempY(0:nz),  NtempZ(0:nz) )
            ALLOCATE( NtempY2(0:nz), NtempZ2(0:nz) )
            ALLOCATE( psitemp(0:ny,0:nz) )
            ALLOCATE( psitemp_new(0:ny,0:nz) )
            ALLOCATE( tempdriver(0:ny,0:nz), thl0(1:jtot,1:ktot) )
        END IF
        
        IF (lmoist) THEN
            ALLOCATE( qtmean(0:nz), R51(0:nz), R52(0:nz), R53(0:nz), R54(0:nz), R55(0:nz) )
            ALLOCATE( a51(0:nz), a52(0:nz), a53(0:nz), a54(0:nz), a55(0:nz) )
            ALLOCATE( tqt_scale(0:nz) )
            ALLOCATE( nlqty(0:nz), nlqtz(0:nz) )
            ALLOCATE( NqtY(0:nz),  NqtZ(0:nz) )
            ALLOCATE( NqtY2(0:nz), NqtZ2(0:nz) )
            ALLOCATE( psiqt(0:ny,0:nz) )
            ALLOCATE( psiqt_new(0:ny,0:nz) )
            ALLOCATE( qtdriver(0:ny,0:nz), qt0(1:jtot,1:ktot) )
        END IF

        !! ################################################################################

        !! Define seed for Fortran inbuilt random number generator ########################
        CALL DATE_AND_TIME(VALUES=current_time)
        time_now = current_time(5) * 3600 + current_time(6) * 60 + current_time(7)
        CALL RANDOM_SEED()
        CALL RANDOM_SEED(SIZE = seed_size)
        ALLOCATE(seed(1:seed_size))
        DO  k = 1,seed_size
            seed(k) = time_now + k
        END DO
        CALL RANDOM_SEED(PUT = seed)
        
        !! ################################################################################
        
        CALL read_synInflow_inputs
        CALL calc_Lund_transformation_coeff

        NUY = 2*nluy
        NUZ = 2*nluz
        NVY = 2*nlvy
        NVZ = 2*nlvz
        NWY = 2*nlwy
        NWZ = 2*nlwz

        NUY2 = 2*NUY
        NUZ2 = 2*NUZ
        NVY2 = 2*NVY
        NVZ2 = 2*NVZ
        NWY2 = 2*NWY
        NWZ2 = 2*NWZ

        IF (ltempeq) THEN
            NtempY = 2*nltempy
            NtempZ = 2*nltempz
            NtempY2 = 2*NtempY
            NtempZ2 = 2*NtempZ
        END IF

        IF (lmoist) THEN
            NqtY = 2*nlqty
            NqtZ = 2*nlqtz
            NqtY2 = 2*NqtY
            NqtZ2 = 2*NqtZ
        END IF

        !! ################################################################################

        ALLOCATE( buy(0:MAXVAL(NUY2),0:nz), buz(0:MAXVAL(NUZ2),0:nz) )
        ALLOCATE( bvy(0:MAXVAL(NVY2),0:nz), bvz(0:MAXVAL(NVZ2),0:nz) )
        ALLOCATE( bwy(0:MAXVAL(NWY2),0:nz), bwz(0:MAXVAL(NWZ2),0:nz) )

        CALL calc_filterCoeff_b(nluy,NUY,NUY2,buy)
        CALL calc_filterCoeff_b(nluz,NUZ,NUZ2,buz)

        CALL calc_filterCoeff_b(nlvy,NVY,NVY2,bvy)
        CALL calc_filterCoeff_b(nlvz,NVZ,NVZ2,bvz)

        CALL calc_filterCoeff_b(nlwy,NWY,NWY2,bwy)
        CALL calc_filterCoeff_b(nlwz,NWZ,NWZ2,bwz)

        !! ################################################################################

        IF (ltempeq) THEN
            ALLOCATE( btempy(0:MAXVAL(NtempY2),0:nz), btempz(0:MAXVAL(NtempZ2),0:nz) )
            CALL calc_filterCoeff_b(nltempy,NtempY,NtempY2,btempy)
            CALL calc_filterCoeff_b(nltempz,NtempZ,NtempZ2,btempz)
        END IF

        !! ################################################################################

        IF (lmoist) THEN
            ALLOCATE( bqty(0:MAXVAL(NqtY2),0:nz), bqtz(0:MAXVAL(NqtZ2),0:nz) )
            CALL calc_filterCoeff_b(nlqty,NqtY,NqtY2,bqty)
            CALL calc_filterCoeff_b(nlqtz,NqtZ,NqtZ2,bqtz)
        END IF

        WRITE(*,*) "Initialization done."
    END SUBROUTINE init_synInflowGen


    SUBROUTINE read_synInflow_inputs
        IMPLICIT NONE

        INTEGER :: k
        LOGICAL :: lexist
        
        ! REAL(KIND=dp), PARAMETER :: zsize = 0.96d0, u_ABL = 0.348d0, z0 = 0.0009d0, von_K = 0.42d0
        ! REAL(KIND=dp) :: dz_t
        ! dz_t = zsize/REAL(nz,KIND=dp)
        ! !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k) SCHEDULE(DYNAMIC)
        ! DO k = 0,nz
        !     umean(k) = (u_ABL/von_K) * LOG( ( REAL(k*dz_t,KIND=dp) + z0 ) / z0 )
        ! END DO
        ! !$OMP END PARALLEL DO
        OPEN(80,FILE='ut.txt')
            READ(80,*)(umean(k),  k = 0, nz)
        CLOSE(80)

        ztotal = SUM( dz(1:nz-1) ) + (dz(0) + dz(nz)) / 2.0d0
        u_bulk = SUM( umean(1:nz-1)*dz(1:nz-1) ) + ( umean(0)*dz(0) + umean(nz)*dz(nz) ) / 2.0d0
        u_bulk = u_bulk/ztotal

        !! ################################################################################

        OPEN(81,FILE='upup.txt')
            READ(81,*)(R11(k),  k = 0, nz)
        CLOSE(81)
        OPEN(82,FILE='upvp.txt')
            READ(82,*)(R21(k),  k = 0, nz)
        CLOSE(82)
        OPEN(83,FILE='vpvp.txt')
            READ(83,*)(R22(k),  k = 0, nz)
        CLOSE(83)
        OPEN(84,FILE='upwp.txt')
            READ(84,*)(R31(k),  k = 0, nz)
        CLOSE(84)
        OPEN(85,FILE='vpwp.txt')
            READ(85,*)(R32(k),  k = 0, nz)
        CLOSE(85)
        OPEN(86,FILE='wpwp.txt')
            READ(86,*)(R33(k),  k = 0, nz)
        CLOSE(86)

        !! ################################################################################

        IF (ltempeq) THEN
            OPEN(87,FILE='thlt.txt')
                READ(87,*)(tempmean(k),  k = 0, nz)
            CLOSE(87)
            R41 = 0.0001d0
            R42 = 0.0001d0
            R43 = 0.0001d0
            R44 = 0.0001d0
            ! OPEN(88,FILE='thlpthlpt.txt')
            !     READ(88,*)(R43(k),  k = 0, nz)
            ! CLOSE(88)
            ! OPEN(89,FILE='wpthlpt.txt')
            !     READ(89,*)(R44(k),  k = 0, nz)
            ! CLOSE(89)
            R41(0) = 0.0d0
            R42(0) = 0.0d0
            R43(0) = 0.0d0
            R44(0) = 0.0d0
        END IF
        
        !! ################################################################################

        IF (lmoist) THEN
            OPEN(90,FILE='qtt.txt')
                READ(90,*)(qtmean(k),  k = 0, nz)
            CLOSE(90)
            R51 = 0.0001d0
            R52 = 0.0001d0
            R53 = 0.0001d0
            R54 = 0.0001d0
            R55 = 0.0001d0
            R51(0) = 0.0d0
            R52(0) = 0.0d0
            R53(0) = 0.0d0
            R54(0) = 0.0d0
            R55(0) = 0.0d0
        END IF

        !! ################################################################################

        IF (lcalc_time_and_length_scale) THEN

            CALL calc_time_and_length_scale

        ELSE

            INQUIRE(FILE='length_time_scales_u.txt',EXIST=lexist)
            IF (lexist) THEN
                OPEN(91, FILE='length_time_scales_u.txt')
                DO k = 0,nz
                    READ(91, *) nluy(k), nluz(k), tu_scale(k)
                END DO
                CLOSE(91)
            ELSE
                WRITE(*,*) " 'length_time_scales_u.txt' does not exist in the work folder."
                STOP 1
            END IF

            OPEN(92, FILE='length_time_scales_v.txt')
            DO k = 0,nz
                READ(92, *) nlvy(k), nlvz(k), tv_scale(k)
            END DO
            CLOSE(92)

            OPEN(93, FILE='length_time_scales_w.txt')
            DO k = 0,nz
                READ(93, *) nlwy(k), nlwz(k), tw_scale(k)
            END DO
            CLOSE(93)

            IF (ltempeq) THEN
                OPEN(94, FILE='length_time_scales_temp.txt')
                DO k = 0,nz
                    READ(94, *) nltempy(k), nltempz(k), ttemp_scale(k)
                END DO
                CLOSE(94)
            END IF

            IF (lmoist) THEN
                OPEN(95, FILE='length_time_scales_qt.txt')
                DO k = 0,nz
                    READ(95, *) nlqty(k), nlqtz(k), tqt_scale(k)
                END DO
                CLOSE(95)
            END IF

        END IF
        
    END SUBROUTINE read_synInflow_inputs



    !! Taken from PALM4u synthetic_turbulence_generator_mod.f90 --> SUBROUTINE calc_length_and_time_scale
    SUBROUTINE calc_time_and_length_scale
        IMPLICIT NONE

        INTEGER :: k
        REAL(KIND=dp) :: length_scale

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k) SCHEDULE(DYNAMIC)
        DO k = 1,nz

            length_scale = 8.0d0 * MIN( dx, dy, dz(k) )
            
            nluy(k) = MAX( INT( length_scale * ddy ), 1 )
            nluz(k) = MAX( INT( length_scale * ddz(k) ), 1 )
            
            nlvy(k) = MAX( INT( length_scale * ddy ), 1 )
            nlvz(k) = MAX( INT( length_scale * ddz(k) ), 1 )
            
            nlwy(k) = MAX( INT( length_scale * ddy ), 1 )
            nlwz(k) = MAX( INT( length_scale * ddz(k) ), 1 )

            tu_scale(k)  = length_scale / ( ABS( u_bulk ) + 1.0E-9 ) 
            tv_scale(k)  = length_scale / ( ABS( u_bulk ) + 1.0E-9 ) 
            tw_scale(k)  = SQRT( tu_scale(k)**2 + tv_scale(k)**2 )

        END DO
        !$OMP END PARALLEL DO
        
        nluy(0) = nluy(1)
        nluz(0) = nluz(1)
        nlvy(0) = nlvy(1)
        nlvz(0) = nlvz(1)
        nlwy(0) = nlwy(1)
        nlwz(0) = nlwz(1)
        tu_scale(0) = tu_scale(1)
        tv_scale(0) = tv_scale(1)
        tw_scale(0) = tw_scale(1)

        !! ################################################################################

        IF (ltempeq) THEN
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k) SCHEDULE(DYNAMIC)
            DO k = 1,nz
                length_scale = 8.0d0 * MIN( dx, dy, dz(k) )
                nltempy(k) = MAX( INT( length_scale * ddy ), 1 )
                nltempz(k) = MAX( INT( length_scale * ddz(k) ), 1 )
                ttemp_scale(k)  = length_scale / ( ABS( u_bulk ) + 1.0E-9 ) 
            END DO
            !$OMP END PARALLEL DO
            
            nltempy(0) = nltempy(1)
            nltempz(0) = nltempz(1)
            ttemp_scale(0) = ttemp_scale(1)
        END IF

        !! ################################################################################

        IF (lmoist) THEN
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k) SCHEDULE(DYNAMIC)
            DO k = 1,nz
                length_scale = 8.0d0 * MIN( dx, dy, dz(k) )
                nlqty(k) = MAX( INT( length_scale * ddy ), 1 )
                nlqtz(k) = MAX( INT( length_scale * ddz(k) ), 1 )
                tqt_scale(k)  = length_scale / ( ABS( u_bulk ) + 1.0E-9 ) 
            END DO
            !$OMP END PARALLEL DO
            
            nlqty(0) = nlqty(1)
            nlqtz(0) = nlqtz(1)
            tqt_scale(0) = tqt_scale(1)
        END IF

    END SUBROUTINE calc_time_and_length_scale



    SUBROUTINE calc_Lund_transformation_coeff   !!  Eq. (18) of Xie and Castro (2008)
        IMPLICIT NONE

        INTEGER :: k

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k) SCHEDULE(DYNAMIC)
        DO k = 0,nz
            a11(k) = SQRT( ABS(R11(k)) )

            a21(k) = R21(k)/(a11(k)+1.0E-9)
            a22(k) = SQRT( ABS( R22(k) - a21(k)*a21(k) ) )

            a31(k) = R31(k)/(a11(k)+1.0E-9)
            a32(k) = ( R32(k) - a21(k)*a31(k) ) / (a22(k)+1.0E-9)
            a33(k) = SQRT( ABS( R33(k) - a31(k)*a31(k) - a32(k)*a32(k) ) )
        END DO
        !$OMP END PARALLEL DO

        !! ################################################################################

        IF (ltempeq) THEN
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k) SCHEDULE(DYNAMIC)
            DO k = 0,nz
                a41(k) = R41(k)/(a11(k)+1.0E-9)
                a42(k) = ( R42(k) - a21(k)*a41(k) ) / (a22(k)+1.0E-9)
                a43(k) = ( R43(k) - a31(k)*a41(k) - a32(k)*a42(k) ) / (a33(k)+1.0E-9)
                a44(k) = SQRT( ABS( R44(k) - a41(k)*a41(k) - a42(k)*a42(k) - a43(k)*a43(k) ) )
            END DO
            !$OMP END PARALLEL DO
        END IF

        !! ################################################################################

        IF (lmoist) THEN
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k) SCHEDULE(DYNAMIC)
            DO k = 0,nz
                a51(k) = R51(k)/(a11(k)+1.0E-9)
                a52(k) = ( R52(k) - a21(k)*a51(k) ) / (a22(k)+1.0E-9)
                a53(k) = ( R53(k) - a31(k)*a51(k) - a32(k)*a52(k) ) / (a33(k)+1.0E-9)
                a54(k) = ( R54(k) - a41(k)*a51(k) - a42(k)*a52(k) - a43(k)*a53(k) ) / (a44(k)+1.0E-9)
                a55(k) = SQRT( ABS( R55(k) - a51(k)*a51(k) - a52(k)*a52(k) - a53(k)*a53(k) - a54(k)*a54(k) ) )
            END DO
            !$OMP END PARALLEL DO
        END IF

    END SUBROUTINE calc_Lund_transformation_coeff



    SUBROUTINE calc_filterCoeff_b(n,NL,NL2,b)
        IMPLICIT NONE

        INTEGER, DIMENSION(0:nz), INTENT(IN) :: n, NL, NL2
        REAL(KIND=dp), DIMENSION(0:MAXVAL(NL2),0:nz), INTENT(OUT) :: b
        
        REAL(KIND=dp), DIMENSION(0:nz) :: b_deno
        INTEGER :: j, k

        !!  Eq. (10) of Xie and Castro (2008)
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,k) SCHEDULE(DYNAMIC)
        DO k = 0,nz
            DO j = 0,NL2(k)
                b(j,k) = EXP( -pi * ( ABS( REAL( j-NL(k), KIND=dp ) ) / REAL( n(k), KIND=dp ) ) )
            END DO
        END DO
        !$OMP END PARALLEL DO
        
        ! !!  Denominator in Eq. (9) of Xie and Castro (2008)
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k) SCHEDULE(DYNAMIC)
        DO k = 0,nz
            b_deno(k) = SQRT( SUM( b(0:NL2(k),k)**2 ) )
        END DO
        !$OMP END PARALLEL DO

        !!  Eq. (9) of Xie and Castro (2008)
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,k) SCHEDULE(DYNAMIC)
        DO k = 0,nz
            DO j = 0,NL2(k)
                b(j,k) = b(j,k) / b_deno(k)
            END DO
        END DO
        !$OMP END PARALLEL DO

    END SUBROUTINE calc_filterCoeff_b



    SUBROUTINE synInflowGen_main
        IMPLICIT NONE

        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: expu1, expu2, expv1, expv2, expw1, expw2
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: exptemp1, exptemp2
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: expqt1, expqt2
        INTEGER :: j, k, it

        ! DO k = 0,nz
        !     WRITE(*,*) k, " ", nluy(k), " ", nluz(k), " ", nlvy(k), " ", nlvz(k), " ", nlwy(k), " ", nlwz(k)
        ! END DO
        ! DO k = 0,nz
        !     WRITE(*,*) k, " ", nltempy(k), " ", nltempz(k), " ", nlqty(k), " ", nlqtz(k)
        ! END DO
        ! DO k = 0,nz
        !     WRITE(*,*) k, " ", tu_scale(k)/dt_sig, " ", tv_scale(k), " ", tw_scale(k)/dt_sig, " ", ttemp_scale(k), " ", tqt_scale(k)
        ! END DO
        
        ALLOCATE( expu1(0:nz), expu2(0:nz), expv1(0:nz), expv2(0:nz), expw1(0:nz), expw2(0:nz) )

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k) SCHEDULE(DYNAMIC)
        DO k = 0,nz
            expu1(k) = EXP(-(pi/2.0d0)*(dt_sig/tu_scale(k)))          !! refer to Eq. (14) of Xie and Castro (2008)
            expu2(k) = SQRT(1.0d0 - EXP(-pi*(dt_sig/tu_scale(k))))    !! refer to Eq. (14) of Xie and Castro (2008)

            expv1(k) = EXP(-(pi/2.0d0)*(dt_sig/tv_scale(k)))
            expv2(k) = SQRT(1.0d0 - EXP(-pi*(dt_sig/tv_scale(k))))

            expw1(k) = EXP(-(pi/2.0d0)*(dt_sig/tw_scale(k)))
            expw2(k) = SQRT(1.0d0 - EXP(-pi*(dt_sig/tw_scale(k))))
        END DO
        !$OMP END PARALLEL DO


        !! ###################### For the first time step #################################

        CALL calc_psi(NUY2,NUZ2,buy,buz,psiu_new)
        CALL calc_psi(NVY2,NVZ2,bvy,bvz,psiv_new)
        CALL calc_psi(NWY2,NWZ2,bwy,bwz,psiw_new)

        psiu = psiu_new
        psiv = psiv_new
        psiw = psiw_new

        !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(j,k) SCHEDULE(DYNAMIC, 8)
        DO k = 0,nz
            DO j = 0,ny
                !!  Eq. (17) of Xie and Castro (2008)
                udriver(j,k) = umean(k) + a11(k) * psiu(j,k)
                vdriver(j,k) = a21(k) * psiu(j,k) + a22(k) * psiv(j,k)
                wdriver(j,k) = a31(k) * psiu(j,k) + a32(k) * psiv(j,k) + a33(k) * psiw(j,k)
            END DO
        END DO
        !$OMP END PARALLEL DO

        CALL massFluxCorrection
        CALL bilinear_interp

        IF (ltempeq) THEN

            ALLOCATE( exptemp1(0:nz), exptemp2(0:nz) )

            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k) SCHEDULE(DYNAMIC)
            DO k = 0,nz
                exptemp1(k) = EXP(-(pi/2.0d0)*(dt_sig/ttemp_scale(k)))
                exptemp2(k) = SQRT(1.0d0 - EXP(-pi*(dt_sig/ttemp_scale(k))))
            END DO
            !$OMP END PARALLEL DO

            CALL calc_psi(NtempY2,NtempZ2,btempy,btempz,psitemp_new)
            psitemp = psitemp_new

            !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(j,k) SCHEDULE(DYNAMIC, 8)
            DO k = 0,nz
                DO j = 0,ny
                    tempdriver(j,k) = tempmean(k) + a41(k) * psiu(j,k) + a42(k) * psiv(j,k) &
                                              + a43(k) * psiw(j,k) + a44(k) * psitemp(j,k)
                END DO
            END DO
            !$OMP END PARALLEL DO

        END IF

        IF (lmoist) THEN

            ALLOCATE( expqt1(0:nz), expqt2(0:nz) )

            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k) SCHEDULE(DYNAMIC)
            DO k = 0,nz
                expqt1(k) = EXP(-(pi/2.0d0)*(dt_sig/tqt_scale(k)))
                expqt2(k) = SQRT(1.0d0 - EXP(-pi*(dt_sig/tqt_scale(k))))
            END DO
            !$OMP END PARALLEL DO

            CALL calc_psi(NqtY2,NqtZ2,bqty,bqtz,psiqt_new)
            psiqt = psiqt_new

            !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(j,k) SCHEDULE(DYNAMIC, 8)
            DO k = 0,nz
                DO j = 0,ny
                    qtdriver(j,k) = qtmean(k) + a51(k) * psiu(j,k) + a52(k) * psiv(j,k) &
                                          + a53(k) * psiw(j,k) + a54(k) * psitemp(j,k) + a55(k) * psiqt(j,k)
                END DO
            END DO
            !$OMP END PARALLEL DO

        END IF
        
        CALL write_synInflow_driver_file(1)

        IF (lcheck_driver_outputs) THEN
            CALL write_driver_file_1(0,ny,0,nz,'u',udriver)
            CALL write_driver_file_1(0,ny,0,nz,'v',vdriver)
            CALL write_driver_file_1(0,ny,0,nz,'w',wdriver)
            IF(ltempeq) CALL write_driver_file_1(0,ny,0,nz,'t',tempdriver)
            IF(lmoist) CALL write_driver_file_1(0,ny,0,nz,'q',qtdriver)
        END IF


        !! ###################### Second time-step onwards ################################
        DO it= 1,nt

            CALL calc_psi(NUY2,NUZ2,buy,buz,psiu_new)
            CALL calc_psi(NVY2,NVZ2,bvy,bvz,psiv_new)
            CALL calc_psi(NWY2,NWZ2,bwy,bwz,psiw_new) 
            
            !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(j,k) SCHEDULE(DYNAMIC, 8)
            DO k = 0,nz
                DO j = 0,ny

                    !!  Eq. (14) of Xie and Castro (2008)
                    psiu(j,k) = psiu(j,k)*expu1(k) + psiu_new(j,k)*expu2(k)
                    psiv(j,k) = psiv(j,k)*expv1(k) + psiv_new(j,k)*expv2(k)
                    psiw(j,k) = psiw(j,k)*expw1(k) + psiw_new(j,k)*expw2(k)

                    !!  Eq. (17) of Xie and Castro (2008)
                    udriver(j,k) = umean(k) + a11(k) * psiu(j,k)
                    vdriver(j,k) = a21(k) * psiu(j,k) + a22(k) * psiv(j,k)
                    wdriver(j,k) = a31(k) * psiu(j,k) + a32(k) * psiv(j,k) + a33(k) * psiw(j,k)

                END DO
            END DO
            !$OMP END PARALLEL DO

            CALL massFluxCorrection

            CALL bilinear_interp
            
            IF (ltempeq) THEN
                CALL calc_psi(NtempY2,NtempZ2,btempy,btempz,psitemp_new)
                !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(j,k) SCHEDULE(DYNAMIC, 8)
                DO k = 0,nz
                    DO j = 0,ny
                        psitemp(j,k) = psitemp(j,k)*exptemp1(k) + psitemp_new(j,k)*exptemp2(k)
                        tempdriver(j,k) = tempmean(k) + a41(k) * psiu(j,k) + a42(k) * psiv(j,k) &
                                                  + a43(k) * psiw(j,k) + a44(k) * psitemp(j,k)
                    END DO
                END DO
                !$OMP END PARALLEL DO
            END IF

            IF (lmoist) THEN
                CALL calc_psi(NqtY2,NqtZ2,bqty,bqtz,psiqt_new)
                !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(j,k) SCHEDULE(DYNAMIC, 8)
                DO k = 0,nz
                    DO j = 0,ny
                        psiqt(j,k) = psiqt(j,k)*expqt1(k) + psiqt_new(j,k)*expqt2(k)
                        qtdriver(j,k) = qtmean(k) + a51(k) * psiu(j,k) + a52(k) * psiv(j,k) &
                                               + a53(k) * psiw(j,k) + a54(k) * psitemp(j,k) + a55(k) * psiqt(j,k)
                    END DO
                END DO
                !$OMP END PARALLEL DO
            END IF

            CALL write_synInflow_driver_file(it+1)

            IF (lcheck_driver_outputs) THEN
                CALL write_driver_file_1(0,ny,0,nz,'u',udriver)
                CALL write_driver_file_1(0,ny,0,nz,'v',vdriver)
                CALL write_driver_file_1(0,ny,0,nz,'w',wdriver)
                IF(ltempeq) CALL write_driver_file_1(0,ny,0,nz,'t',tempdriver)
                IF(lmoist) CALL write_driver_file_1(0,ny,0,nz,'q',qtdriver)
            END IF

        END DO

        DEALLOCATE(expu1, expu2, expv1, expv2, expw1, expw2)
        IF (ltempeq) DEALLOCATE( exptemp1, exptemp2 )
        IF (lmoist) DEALLOCATE( expqt1, expqt2 )

        WRITE(*,*) "Main computation completed successfully."

    END SUBROUTINE synInflowGen_main



    SUBROUTINE calc_psi(NLY2,NLZ2,by,bz,psi)
        IMPLICIT NONE

        INTEGER, DIMENSION(0:nz), INTENT(IN) :: NLY2, NLZ2
        REAL(KIND=dp), DIMENSION(0:MAXVAL(NLY2),0:nz), INTENT(IN) :: by
        REAL(KIND=dp), DIMENSION(0:MAXVAL(NLZ2),0:nz), INTENT(IN) :: bz
        REAL(KIND=dp), DIMENSION(0:ny,0:nz), INTENT(OUT) :: psi

        REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: randn
        
        INTEGER :: nry, nrz
        INTEGER :: iy, iz, j, k
        
        nry  = ny + MAX( MAXVAL(NLY2), MAXVAL(NLZ2) )*2
        nrz  = nz + MAX( MAXVAL(NLY2), MAXVAL(NLZ2) )*2

        ALLOCATE( randn(0:nry,0:nrz) )
        
        CALL generate_normal_rand(nry,nrz,randn)
        
        psi = 0.0d0

        !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(iy,iz,j,k) SCHEDULE(DYNAMIC, 8)
        DO iz = 0,nz
            DO iy = 0,ny 
                DO k = 0,NLZ2(iz)
                    DO j = 0,NLY2(iz)
                        psi(iy,iz) = psi(iy,iz) + by(j,iz)*bz(k,iz)*randn(iy+j,iz+k)
                    END DO
                END DO
            END DO
        END DO
        !$OMP END PARALLEL DO

        DEALLOCATE(randn)

    END SUBROUTINE calc_psi



    SUBROUTINE generate_normal_rand(nry,nrz,rndm)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: nry, nrz
        REAL(KIND=dp), DIMENSION(0:nry,0:nrz), INTENT(OUT) :: rndm
        
        REAL(KIND=dp) :: mean_rndm, size_rndm

        size_rndm = REAL((nry+1)*(nrz+1), KIND=dp)

        CALL RANDOM_NUMBER(rndm)
        
        rndm = rndm - SUM(rndm)/size_rndm
        mean_rndm = SUM(rndm)/size_rndm
        rndm = rndm / SQRT( SUM( (rndm-mean_rndm)**2 ) / size_rndm )
        ! WRITE(*,*) SUM(rndm)/size_rndm, "  ", SQRT( SUM( rndm**2 ) / size_rndm )

    END SUBROUTINE generate_normal_rand



    SUBROUTINE massFluxCorrection
        IMPLICIT NONE

        REAL(KIND=dp) :: u_bulk_compute
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: unew
        INTEGER :: iz

        ALLOCATE( unew(0:nz) )
        
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iz) SCHEDULE(DYNAMIC)
        DO iz = 0,nz
            unew(iz) = ( SUM( udriver(1:ny-1,iz) ) + ( (udriver(0,iz) + udriver(nz,iz)) / 2.0d0 ) ) &
                        / REAL(ny,KIND=dp)
        END DO
        !$OMP END PARALLEL DO
        
        u_bulk_compute = SUM( unew(1:nz-1)*dz(1:nz-1) ) + ( unew(0)*dz(0) + unew(nz)*dz(nz) ) / 2.0d0
        u_bulk_compute = u_bulk_compute/ztotal      

        udriver(:,:) = udriver(:,:) * ( u_bulk / u_bulk_compute )
        
        DEALLOCATE( unew )

    END SUBROUTINE massFluxCorrection



    SUBROUTINE bilinear_interp
        IMPLICIT NONE

        INTEGER :: j, k

        !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(j,k) SCHEDULE(DYNAMIC,8)
        DO k = 1,ktot
            DO j = 1,jtot
                u0(j,k) = ( ( ( udriver(j-1,k-1) + udriver(j,k-1) ) / 2.0d0 ) &
                          + ( ( udriver(j-1,k) + udriver(j,k) ) / 2.0d0 ) ) / 2.0d0
            END DO
        END DO
        !$OMP END PARALLEL DO

        !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(j,k) SCHEDULE(DYNAMIC,8)
        DO k = 1,ktot
            DO j = 1,jtot
                v0(j,k) = ( vdriver(j-1,k-1) + vdriver(j-1,k) ) / 2.0d0
            END DO
        END DO
        !$OMP END PARALLEL DO

        !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(j,k) SCHEDULE(DYNAMIC,8)
        DO k = 1,ktot
            DO j = 1,jtot
                w0(j,k) = ( wdriver(j-1,k-1) + wdriver(j,k-1) ) / 2.0d0
            END DO
        END DO
        !$OMP END PARALLEL DO

        IF (ltempeq) THEN
            !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(j,k) SCHEDULE(DYNAMIC,8)
            DO k = 1,ktot
                DO j = 1,jtot
                    thl0(j,k) = ( ( ( tempdriver(j-1,k-1) + tempdriver(j,k-1) ) / 2.0d0 ) &
                                   + ( ( tempdriver(j-1,k) + tempdriver(j,k) ) / 2.0d0 ) ) / 2.0d0
                END DO
            END DO
            !$OMP END PARALLEL DO
        END IF

        IF (lmoist) THEN
            !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(j,k) SCHEDULE(DYNAMIC,8)
            DO k = 1,ktot
                DO j = 1,jtot
                    qt0(j,k) = ( ( ( qtdriver(j-1,k-1) + qtdriver(j,k-1) ) / 2.0d0 ) &
                               + ( ( qtdriver(j-1,k) + qtdriver(j,k) ) / 2.0d0 ) ) / 2.0d0
                END DO
            END DO
            !$OMP END PARALLEL DO
        END IF

    END SUBROUTINE bilinear_interp



    SUBROUTINE write_driver_file_1(js,je,ks,ke,ch,vel)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: js,je,ks,ke
        CHARACTER, INTENT(IN) :: ch
        REAL(KIND=dp), DIMENSION(js:je,ks:ke), INTENT(IN) :: vel

        INTEGER, PARAMETER :: file_unit = 10
        INTEGER :: j, k, ios
        LOGICAL :: lexist
        CHARACTER*12 :: filename

        filename = ch // "_driver.txt"

        INQUIRE(FILE=filename,EXIST=lexist)
        IF (lexist) THEN
            OPEN(UNIT=file_unit,FILE=filename,FORM='FORMATTED',STATUS='OLD', &
                    POSITION='APPEND',ACTION='WRITE')
        ELSE
            OPEN(UNIT=file_unit,FILE=filename,FORM='FORMATTED',STATUS='REPLACE', &
                    POSITION='APPEND',ACTION='WRITE',IOSTAT=ios)
            IF (IOS > 0) WRITE(6,*) 'IOS = ',IOS,' while creating ',filename,' in write_driver_file_1'
        END IF
        DO k = ks,ke
            DO j = js,je
                WRITE(UNIT=file_unit,FMT='(f15.10,x)',ADVANCE='NO') vel(j,k)
            END DO
        END DO
        CLOSE(UNIT=file_unit)

    END SUBROUTINE write_driver_file_1



    SUBROUTINE write_synInflow_driver_file(nstepdriver)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: nstepdriver

        REAL(KIND=dp) :: tdriverstart, timee
        INTEGER :: y_driver_len

        INTEGER :: ios, driverid
        INTEGER :: filesizet, filesizev
        CHARACTER(15) :: name
        CHARACTER(3) :: cdriverid
        LOGICAL :: lexist
        ! use decomp_2d, only : zstart, zend
        
        IF (MOD(jtot,nprocy)/=0) THEN
            WRITE(*,*) "ERROR !! jtot must be divisible by nprocy; give correct input parameters."
            STOP
        END IF
        y_driver_len = jtot/nprocy
        tdriverstart = 0.0d0 * dt_sig


        IF (nstepdriver==1) THEN
            WRITE(6,*) '=================================================================='
            WRITE(6,*) '*** Starting to write data for driver simulation ***'
            WRITE(6,*) 'Driver recording variables:'
            WRITE(6,'(A,F9.2,A,I8,A,F12.9)') ' Starting time: ', tdriverstart, &
                                            ';    Total number of time steps to be stored: ', nt+1, &
                                            ';    Inlet record intervals: ', dt_sig
            WRITE(6,*) '=================================================================='
        END IF
        
        WRITE(6,*) '============ Writing driver files ============'
        WRITE(*,*) 'Driver timestep: ', nstepdriver

        timee = (nstepdriver-1)*dt_sig
        INQUIRE(IOLENGTH=filesizet)(timee-tdriverstart)
        
        name = 'tdriver_   .'
        name(9:11)= '000'
        name(13:15)= cexpnr
        INQUIRE(FILE=name,EXIST=lexist)
        IF (lexist) THEN
            ! WRITE(6,*) 'Writing Time stamp to file: ', name
            OPEN(UNIT=11,FILE=name,FORM='UNFORMATTED',STATUS='OLD',ACCESS='DIRECT', &
                    RECL=filesizet,ACTION='WRITE')
        ELSE
            ! WRITE(6,*) 'Creating Time stamp driver file: ', name
            OPEN(UNIT=11,FILE=name,FORM='UNFORMATTED',STATUS='REPLACE',ACCESS='DIRECT', &
                    RECL=filesizet,ACTION='WRITE',IOSTAT=ios)
            IF (IOS > 0) THEN
            WRITE(6,*) 'IOS = ',IOS,' while creating ',name,' file in modSynInflowGen'
            END IF
        END IF
        WRITE(11,REC=nstepdriver) (timee-tdriverstart)
        CLOSE(UNIT=11)
        WRITE(*,*) 'Driver time:' , timee-tdriverstart


        INQUIRE(IOLENGTH=filesizev)u0(:,:)
        
        DO driverid = 0,nprocy-1
            
            WRITE(cdriverid,'(I3.3)') driverid
            
            name = 'udriver_   .'
            name(9:11)= cdriverid
            name(13:15)= cexpnr
            INQUIRE(FILE=name,EXIST=lexist)
            IF (lexist) THEN
                ! WRITE(6,*) 'Writing Inlet u-velocity to file: ', name
                OPEN(UNIT=11,FILE=name,FORM='UNFORMATTED',STATUS='OLD',ACCESS='DIRECT', &
                        RECL=filesizev,ACTION='WRITE')
            ELSE
                ! WRITE(6,*) 'Creating Inlet u-velocity inlet file: ', name
                OPEN(UNIT=11,FILE=name,FORM='UNFORMATTED',STATUS='REPLACE',ACCESS='DIRECT', &
                        RECL=filesizev,ACTION='WRITE')
            END IF
            WRITE(11,REC=nstepdriver) (u0(y_driver_len*driverid+1:y_driver_len*(driverid+1),:))
            ! WRITE(11,REC=nstepdriver) (u0(zstart(2):zend(2),:))
            CLOSE(UNIT=11)

            name = 'vdriver_   .'
            name(9:11)= cdriverid
            name(13:15)= cexpnr
            INQUIRE(FILE=name,EXIST=lexist)
            IF (lexist) THEN
                ! WRITE(6,*) 'Writing Inlet v-velocity to file: ', name
                OPEN(UNIT=11,FILE=name,FORM='UNFORMATTED',STATUS='OLD',ACCESS='DIRECT', &
                        RECL=filesizev,ACTION='WRITE')
            ELSE
                ! WRITE(6,*) 'Creating v-velocity inlet file: ', name
                OPEN(UNIT=11,FILE=name,FORM='UNFORMATTED',STATUS='REPLACE',ACCESS='DIRECT', &
                        RECL=filesizev,ACTION='WRITE')
            END IF
            WRITE(11,REC=nstepdriver) (v0(y_driver_len*driverid+1:y_driver_len*(driverid+1),:))
            ! WRITE(11,REC=nstepdriver) (v0(zstart(2):zend(2),:))
            CLOSE (UNIT=11)

            name = 'wdriver_   .'
            name(9:11)= cdriverid
            name(13:15)= cexpnr
            INQUIRE(FILE=name,EXIST=lexist)
            IF (lexist) THEN
                ! WRITE(6,*) 'Writing Inlet w-velocity to file: ', name
                OPEN(UNIT=11,FILE=name,FORM='UNFORMATTED',STATUS='OLD',ACCESS='DIRECT', &
                        RECL=filesizev,ACTION='WRITE')
            ELSE
                ! WRITE(6,*) 'Creating w-velocity inlet file: ', name
                OPEN(UNIT=11,FILE=name,FORM='UNFORMATTED',STATUS='REPLACE',ACCESS='DIRECT', &
                        RECL=filesizev,ACTION='WRITE')
            END IF
            WRITE(11,REC=nstepdriver) (w0(y_driver_len*driverid+1:y_driver_len*(driverid+1),:))
            ! WRITE(11,REC=nstepdriver) (w0(zstart(2):zend(2),:))
            CLOSE(UNIT=11)

            IF (ltempeq) THEN
                name = 'hdriver_   .'
                name(9:11)= cdriverid
                name(13:15)= cexpnr
                INQUIRE(FILE=name,EXIST=lexist)
                IF (lexist) THEN
                    ! WRITE(6,*) 'Writing Inlet temperature to file: ', name
                    OPEN(UNIT=11,FILE=name,FORM='UNFORMATTED',STATUS='OLD',ACCESS='DIRECT', &
                        RECL=filesizev,ACTION='WRITE')
                ELSE
                    ! WRITE(6,*) 'Creating temperature inlet file: ', name
                    OPEN(UNIT=11,FILE=name,FORM='UNFORMATTED',STATUS='REPLACE',ACCESS='DIRECT', &
                            RECL=filesizev,ACTION='WRITE')
                END IF
                WRITE(11,REC=nstepdriver) (thl0(y_driver_len*driverid+1:y_driver_len*(driverid+1),:))
                ! WRITE(11,REC=nstepdriver) (thl0(zstart(2):zend(2),:))
                CLOSE (UNIT=11)
            END IF

            IF (lmoist) THEN
                name = 'qdriver_   .'
                name(9:11)= cdriverid
                name(13:15)= cexpnr
                INQUIRE(FILE=name,EXIST=lexist)
                IF (lexist) THEN
                    ! WRITE(6,*) 'Writing Inlet moisture to file: ', name
                    OPEN(UNIT=11,FILE=name,FORM='UNFORMATTED',STATUS='OLD',ACCESS='DIRECT', &
                            RECL=filesizev,ACTION='WRITE')
                ELSE
                    ! WRITE(6,*) 'Creating moisture inlet file: ', name
                    OPEN(UNIT=11,FILE=name,FORM='UNFORMATTED',STATUS='REPLACE',ACCESS='DIRECT', &
                            RECL=filesizev,ACTION='WRITE')
                END IF
                ! write(ifoutput)  (((storew0driver (j,k,n),j=jb,je),k=kb,ke+1),n=1,nstore)
                WRITE(11,REC=nstepdriver) (qt0(y_driver_len*driverid+1:y_driver_len*(driverid+1),:))
                ! WRITE(11,REC=nstepdriver) (qt0(zstart(2):zend(2),:))
                CLOSE (UNIT=11)
            END IF

        END DO

    END SUBROUTINE write_synInflow_driver_file



    SUBROUTINE read_synInflow_driver_file
        IMPLICIT NONE

        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: storetdriver
        REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: storeu0driver, storev0driver, storew0driver
        REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: storethl0driver, storeqt0driver
        REAL(KIND=dp) :: tdriverstart, timee
        CHARACTER(15) :: name
        CHARACTER(3) :: cdriverid
        INTEGER :: driverjobnr
        INTEGER :: ios, driverid, filesize
        INTEGER :: driverstore, nstepdriver, y_driver_len, j, k

        INTEGER :: BCxT, BCxT_driver

        y_driver_len = (ny+1)/nprocy
        driverstore = nt+1
        driverjobnr = iexpnr
        tdriverstart = 0.0d0 * dt_sig
        timee = 0

        IF (ltempeq) BCxT = 0
        IF (ltempeq) BCxT_driver = 0

        ALLOCATE( storetdriver(1:driverstore), storeu0driver(1:jtot,1:ktot,1:driverstore), &
                  storev0driver(1:jtot,1:ktot,1:driverstore), storew0driver(1:jtot,1:ktot,1:driverstore) )
        IF (ltempeq .AND. (BCxT == BCxT_driver)) ALLOCATE( storethl0driver(1:jtot,1:ktot,1:driverstore) )
        IF (lmoist) ALLOCATE( storeqt0driver(1:jtot,1:ktot,1:driverstore) )
        
        WRITE(*,*) '========================================================================'
        WRITE(*,*) '*** Reading precursor driver simulation ***'
        
        !##############################################################################################
        name = 'tdriver_   .'
        name(9:11) = '000'
        WRITE(name(13:15),'(i3.3)') driverjobnr

        INQUIRE(FILE=name,SIZE=filesize)

        WRITE(6,*) 'Reading time stamps: ', name
        WRITE(6,*) 'driverstore: ', driverstore
        WRITE(6,*) 'File size of time in bytes (/8) = ', filesize
        
        INQUIRE(IOLENGTH=filesize)(timee-tdriverstart)
        OPEN(UNIT=11,FILE=name,FORM='UNFORMATTED',STATUS='OLD',ACCESS='DIRECT', &
                    RECL=filesize,ACTION='READ',IOSTAT=ios)
        IF (IOS > 0) THEN
            WRITE(6,*) 'IOS = ',ios,' while openning ',name,' for reading driver data in modSynInflowGen'
        END IF
        DO nstepdriver = 1,driverstore
            READ(11, REC=nstepdriver, IOSTAT=ios) storetdriver(nstepdriver)
            IF (ios>0) THEN
                WRITE(6,*) 'IOS = ',ios,' while reading ',name,' in modSynInflowGen'
            ELSEIF (IOS<0) THEN
                WRITE(6,*) 'n =', nstepdriver
            END IF
            WRITE(6,'(A,e20.12)') ' Reading t:', storetdriver(nstepdriver)    
        END DO
        storetdriver = storetdriver + timee ! in case using a warmstart...
        CLOSE(UNIT=11)

        !##############################################################################################
        DO driverid = 0,nprocy-1
                
            WRITE(cdriverid,'(I3.3)') driverid
        
            name = 'udriver_   .'
            name(9:11) = cdriverid
            WRITE(name(13:15),'(i3.3)') driverjobnr
            WRITE(6,*) 'Reading Driver u-velocity: ', name
            INQUIRE(IOLENGTH=filesize)storeu0driver(:,:,0)
            IF (driverid==0) WRITE(6,*) 'record length ',filesize
            OPEN(UNIT=11,FILE=name,FORM='UNFORMATTED',STATUS='OLD',ACCESS='DIRECT', &
                         RECL=filesize,ACTION='READ')
            DO nstepdriver = 1,driverstore
                READ(11,REC=nstepdriver) &
                ((storeu0driver(j,k,nstepdriver),j=y_driver_len*driverid+1,y_driver_len*(driverid+1)),k=1,ktot)
            END DO
            CLOSE(UNIT=11)

            name = 'vdriver_   .'
            name(9:11) = cdriverid
            WRITE(name(13:15),'(i3.3)') driverjobnr
            WRITE(6,*) 'Reading Driver v-velocity: ', name
            OPEN(UNIT=11,FILE=name,FORM='UNFORMATTED',STATUS='OLD',ACCESS='DIRECT', &
                         RECL=filesize,ACTION='READ')
            DO nstepdriver = 1,driverstore
                READ(11,REC=nstepdriver) &  !((storev0driver (j,k,n),j=jb-jh,je+jh),k=kb-kh,ke+kh)
                ((storev0driver(j,k,nstepdriver),j=y_driver_len*driverid+1,y_driver_len*(driverid+1)),k=1,ktot)
            END DO
            CLOSE(UNIT=11)

            name = 'wdriver_   .'
            name(9:11) = cdriverid
            WRITE(name(13:15),'(i3.3)') driverjobnr
            WRITE(6,*) 'Reading Driver w-velocity: ', name
            OPEN(UNIT=11,FILE=name,FORM='UNFORMATTED',STATUS='OLD',ACCESS='DIRECT', &
                         RECL=filesize,ACTION='READ')
            DO nstepdriver = 1,driverstore
                READ(11,REC=nstepdriver) &  !((storew0driver (j,k,n),j=jb-jh,je+jh),k=kb-kh,ke+kh)
                ((storew0driver(j,k,nstepdriver),j=y_driver_len*driverid+1,y_driver_len*(driverid+1)),k=1,ktot)
            END DO
            CLOSE(UNIT=11)

            IF (ltempeq .AND. (BCxT == BCxT_driver)) THEN
                name = 'hdriver_   .'
                name(9:11) = cdriverid
                WRITE(name(13:15),'(i3.3)') driverjobnr
                WRITE(6,*) 'Reading Driver temperature: ', name
                OPEN(UNIT=11,FILE=name,FORM='UNFORMATTED',STATUS='OLD',ACCESS='DIRECT', &
                             RECL=filesize,ACTION='READ')
                DO nstepdriver = 1,driverstore
                    READ(11,REC=nstepdriver) &  !((storethl0driver (j,k,n),j=jb-jh,je+jh),k=kb-kh,ke+kh)
                    ((storethl0driver(j,k,nstepdriver),j=y_driver_len*driverid+1,y_driver_len*(driverid+1)),k=1,ktot)
                END DO
                CLOSE(UNIT=11)
            END IF
            
            IF (lmoist) THEN
                name = 'qdriver_   .'
                name(9:11)= cdriverid
                WRITE(name(13:15),'(i3.3)') driverjobnr
                WRITE(6,*) 'Reading Driver moisture: ', name
                OPEN(UNIT=11,FILE=name,FORM='UNFORMATTED',STATUS='OLD',ACCESS='DIRECT', &
                             RECL=filesize,ACTION='READ')
                DO nstepdriver = 1,driverstore
                    READ(11,REC=nstepdriver) &  !((storeqt0driver (j,k,n),j=jb-jh,je+jh),k=kb-kh,ke+kh)
                    ((storeqt0driver(j,k,nstepdriver),j=y_driver_len*driverid+1,y_driver_len*(driverid+1)),k=1,ktot)
                END DO
                CLOSE(UNIT=11)
            END IF

        END DO
        
        !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(5)
            !$OMP SECTIONS

                !$OMP SECTION
                !$    CALL write_driver_file(1,jtot,1,ktot,1,driverstore,'a',storeu0driver,OMP_GET_THREAD_NUM())

                !$OMP SECTION
                !$    CALL write_driver_file(1,jtot,1,ktot,1,driverstore,'b',storev0driver,OMP_GET_THREAD_NUM())

                !$OMP SECTION
                !$    CALL write_driver_file(1,jtot,1,ktot,1,driverstore,'c',storew0driver,OMP_GET_THREAD_NUM())

                !$OMP SECTION
                !$    IF (ltempeq .AND. (BCxT == BCxT_driver)) THEN
                !$        CALL write_driver_file(1,jtot,1,ktot,1,driverstore,'d',storethl0driver,OMP_GET_THREAD_NUM())
                !$    END IF

                !$OMP SECTION
                !$    IF (lmoist) CALL write_driver_file(1,jtot,1,ktot,1,driverstore,'e',storeqt0driver,OMP_GET_THREAD_NUM())

            !$OMP END SECTIONS
        !$OMP END PARALLEL
        
        ! CALL write_driver_file(1,jtot,1,ktot,1,driverstore,'a',storeu0driver,10)
        ! CALL write_driver_file(1,jtot,1,ktot,1,driverstore,'b',storev0driver,10)
        ! CALL write_driver_file(1,jtot,1,ktot,1,driverstore,'c',storew0driver,10)
        ! IF (ltempeq .AND. (BCxT == BCxT_driver)) THEN
        !     CALL write_driver_file(1,jtot,1,ktot,1,driverstore,'d',storethl0driver,10)
        ! END IF
        ! IF (lmoist) CALL write_driver_file(1,jtot,1,ktot,1,driverstore,'e',storeqt0driver,10)

        DEALLOCATE( storetdriver,storeu0driver,storev0driver,storew0driver )
        IF (ltempeq .AND. (BCxT == BCxT_driver)) DEALLOCATE( storethl0driver )
        IF (lmoist) DEALLOCATE( storeqt0driver )

    END SUBROUTINE read_synInflow_driver_file

    SUBROUTINE write_driver_file(js,je,ks,ke,ts,te,ch,vel,fileID)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: js,je,ks,ke,ts,te, fileID
        CHARACTER, INTENT(IN) :: ch
        REAL(KIND=dp), DIMENSION(js:je,ks:ke,ts:te), INTENT(IN) :: vel

        INTEGER :: j, k, it
        CHARACTER*12 :: filename

        filename = ch // "_driver.txt"
        
        OPEN(UNIT=fileID,FILE=filename)
        DO it = ts,te
            DO k = ks,ke
                DO j = js,je
                    WRITE(UNIT=fileID,FMT='(f15.10,x)',ADVANCE='NO') vel(j,k,it)
                END DO
            END DO
            WRITE(UNIT=fileID,FMT='(A)',ADVANCE='NO') NEW_LINE('a')
        END DO
        CLOSE(UNIT=fileID)

        WRITE(*,*) filename, " written successfully."
    END SUBROUTINE write_driver_file

    SUBROUTINE finish_synInflowGen
        IMPLICIT NONE

        ! CALL read_synInflow_driver_file

        DEALLOCATE( seed )
        DEALLOCATE( dz, ddz, yf, yh, zf, zh, ydriver, zdriver )
        DEALLOCATE( umean, R11, R21, R22, R31, R32, R33 )
        DEALLOCATE( a11, a21, a22, a31, a32, a33 )
        DEALLOCATE( tu_scale, tv_scale, tw_scale )
        DEALLOCATE( nluy, nluz, nlvy, nlvz, nlwy, nlwz )
        DEALLOCATE( NUY,  NUZ,  NVY,  NVZ,  NWY,  NWZ )
        DEALLOCATE( NUY2, NUZ2, NVY2, NVZ2, NWY2, NWZ2 )
        DEALLOCATE( buy, buz, bvy, bvz, bwy, bwz )
        DEALLOCATE( psiu, psiv, psiw )
        DEALLOCATE( psiu_new, psiv_new, psiw_new )
        DEALLOCATE( udriver, vdriver, wdriver, u0, v0, w0 )

        IF (ltempeq) THEN
            DEALLOCATE( tempmean, R41, R42, R43, R44 )
            DEALLOCATE( a41, a42, a43, a44 )
            DEALLOCATE( ttemp_scale )
            DEALLOCATE( nltempy, nltempz )
            DEALLOCATE( NtempY,  NtempZ )
            DEALLOCATE( NtempY2, NtempZ2 )
            DEALLOCATE( btempy, btempz )
            DEALLOCATE( psitemp )
            DEALLOCATE( psitemp_new )
            DEALLOCATE( tempdriver, thl0 )
        END IF

        IF (lmoist) THEN
            DEALLOCATE( qtmean, R51, R52, R53, R54 )
            DEALLOCATE( a51, a52, a53, a54, a55 )
            DEALLOCATE( tqt_scale )
            DEALLOCATE( nlqty, nlqtz )
            DEALLOCATE( NqtY,  NqtZ )
            DEALLOCATE( NqtY2, NqtZ2 )
            DEALLOCATE( bqty, bqtz )
            DEALLOCATE( psiqt )
            DEALLOCATE( psiqt_new )
            DEALLOCATE( qtdriver, qt0 )
        END IF

        WRITE(*,*) "Synthetic inflow generator completed successfully."
    END SUBROUTINE finish_synInflowGen

END PROGRAM synInflowGen


!! REFERENCES:
!! [1] Xie, Zheng-Tong, and Ian P. Castro. "Efficient generation of inflow conditions for large 
!!     eddy simulation of street-scale flows." Flow, turbulence and combustion 81 (2008): 449-470.
!!
!! [2] Kim, Yusik, Ian P. Castro, and Zheng-Tong Xie. "Divergence-free turbulence inflow conditions 
!!     for large-eddy simulations with incompressible flow solvers." Computers & Fluids 84 (2013): 56-68.
