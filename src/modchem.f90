module modchem
! version 1.0  3 chem_components on input and  4 on output but only interger coefficents.
! version 1.1  4 chem_components on input and 4 on output real on output allowed
! 2009_07_09 bug fixed wrong coefficients
! 2009_07_10  1916 and 1954       write(coef_str,'(f4.2)')PL_scheme(i)%PL(j)%coef     f4.2 ipv i2
! version 1.2 4 chem on input and NCCAA(zie modchem) on output on change recompile


  !____________________SETTINGS_AND_SWITCHES________________________
  !                     IN &NAMCHEM
  !
  !       &NAMCHEM
  !         lchem       = .false.  SWITCH: switches chemistry module on and off
  !         tnor        = 17       total number of reactions (at least >= number of reactions)
  !         firstchem   = 1        number of the colum in scalar.inp of the first chemical
  !         lastchem    = 15       number of the colum in scalar.inp of the last chemical
  !                                only important if there are other scalars available in scalar.inp
  !         ldiuvar     = .true.   SWITCH: diurnal photolysis reaction rates, If false it uses h_ref to calculate
  !                                a constant photolysis rate
  !         h_ref       = 12.      if above is FALSE use this hour to calculate the photolysis rates
  !         lchconst    = .false.  SWITCH: if TRUE then the calculation of the reaction rates
  !                                        uses t_ref, p_ref, q_ref and not the model temp, humidity and pressure
  !         lcloudKconst   = .false.  SWITCH: if TRUE then if there are clouds they will not change the K of the photolysis reactions
  !         t_ref	    = 298.
  !         p_ref	    = 100000.
  !         q_ref	    = 5.e-3
  !         lchmovie    = .false.  SWITCH: if TRUE gives extra output for movies(experimental)
  !         dtchmovie   = 60.      when to write extra output
  !         lsegr       = .false.  SWITCH: if TRUE gives information about segregation in a Mixed Layer approach
  !       /
  !-----------------------------------------------------------------

  !
  ! EXSAMPLE chem.imp.xxx file

! # FORMAT OF inputchem
! # EVERY ITEM SHOULD BE SEPERATED BY AT LEAST ONE SPACE
! # EXCEPT FOR THE COEFFICIENTS IN THE REACTIONS
! # PUT NONACTIVE CHEMICAL COMPONENTS WITH COEFFICIENTS IN ().
! # No empty lines are allowed. Comment line should start with a #.
! # The file should start with a line with a @ as the first character, the rest of the line is free.
! # and should be followed by 3 lines with no comments in between.
! # The first line should contain the chemical species in same order as  in scalar.inp, follwoed by
! # 2 lines with the atol and rtol value's
! # After the last reaction there should be a line with a $ as the first character.
! @  1       2       3       4       5      6       7       8       9       10      11       12      13      14     15
!   O3      NO      NO2     NO3     H2OP    HNO3    RH      R       HO      HO2     H2O2    CO       CO2     H2O    INERT
!   1e-5    1e-5    1e-5    1e-5    1e-5    1e-5    1e-5    1e-5    1e-5    1e-5    1e-5    1e-5    1e-5    1e-5    1e-5
!   0.01	  0.01 	  0.01    0.01	  0.01 	  0.01	  0.01 	  0.01	  0.01    0.01	  0.01     0.01    0.01    0.01   0.01
! #
! #--------------------------------------------------------------------------------------------------------------
! # kr2nd  |Name |Rad |func1|  A    |   B   |   C   |   D   |  E   |   F  |  G   | chemical reaction
! #        |     |Dep |     |       |       |       |       |      |      |      | with inactive species
! #        |5char|int |int  | real  | real  | real  | real  | real | real | real |	 in ( )
! #---------------------------------------------------------------------------------------------------------------
! # all reactions must have a name
! 2.396E-13 R_CO   0    1     1.0      1.0     1.0    1.0      1.0    1.0    1.0   HO + CO    -> HO2 + CO2
! 7.197E-11 R_RH   0    2    2.33e-11  444     1.0    1.0      1.0    1.0    1.0   RH + HO    -> HO2 + R
! #HxOy
! 0.2708E+1 R_23   0    2    4.8e-11   250     1.0    1.0      1.0    1.0    1.0   HO + HO2   -> H2OP + (O2)
! #in the above reaction we put H2OP (instead of H2O) to treat it as a separate chemical
! 0.167E-2  R_25   0    2    1.7e-12   -940    1.0    1.0      1.0    1.0    1.0   HO + O3    -> HO2 + (O2)
! 0.713E-1  R_26a  0    6    2.2e-13   600    1.9e-33 980     1.4e-21 2200   1.0   2HO2       -> H2O2 + (O2)
! 0.492E-4  R_28   0    3    2.03e-16  300     4.57   693      1.0    1.0    1.0   HO2 + O3   -> HO + (2O2)
! #NOx
! 0.271E+0  R_43   0    4    3.3e-30   -3.0    0.0    4.1e-11  0.0    0.0    0.4   HO + NO2   -> HNO3
! 0.217E+0  R_45   0    2    3.6e-12   270     1.0    1.0      1.0    1.0    1.0   HO2 + NO   -> NO2 + HO
! #Atkinson
! 0.443E-3  R_54A  0    2    1.4e-12  -1310    1.0    1.0      1.0    1.0    1.0   NO + O3    -> NO2
! #The 2 following reactions are radiation dependent
! 5.0667e-6 R_O3   1    4    6.29e-5   2.412   0.7  2.2e-10  2.9e-11  1.0    1.0   O3         -> 2HO + (O2)
! 0.0167    R_NO2  1    2    1.67e-2   -.575   1.2    1.0     1.0     1.0    1.0   NO2 + (O2) -> NO + O3
! $ end of chemical reactions specified by $ as first character on the line.
!
! # Reaction: R_CO   R_RH   R_23   R_25  R_26a  R_28  R_43  R_45  R_54A  R_O3  R_NO2
!# function:  1      2      2      2     6      3     4     2     2      4     2
! #
! #   function = 0 => Meaning K colum kn2rd in PPB*sec no temp dependence
! #   function = 1 => Meaning K colum kn2rd in cm3/molecule*sec no temp dependence
! #   function = 2 => K = K (cm3/molecule*sec)
! #            MEANING of A..G         K = A * exp (B/T)
! #   function = 3 => K = K (cm3/molecule*s)
! #            MEANING of A..G         K = A * (T/B)^C * exp (D/T)
! #   function = 4 => K = K'* K"/(K'+ K")* G (cm3/molecule*sec)
! #            MEANING of A..G         K'= A * (T/300)^B * exp(C/T)[M]     K" = D * (T/300)^E * exp(F/T)   Fc = G
! #   function = 5 => K = K'* K"/(K'+ K")* G (/sec)
! #            MEANING of A..G         K'= A * (T/300)^B * exp(C/T)[M]     K" = D * (T/300)^E * exp(F/T)   Fc = G
! #   function = 6 => K =  (K'+ K")*(1+K'") (cm3/molecule*s)
! #            special for 2HO2 ->H2O2 + (O2) because there are 2 reactions
! #            MEANING of A..G         K'= A * exp (B/T)   K"= C * exp(D/T)[M]  K'" = 1 + E * exp(F/T)[H20]
! #   function = 7 => K = K (cm6/molecule2*s)for above R63Ab reaction K = A * (T/1)^0 * exp(0/T)=A*1*1 = A
! #            MEANING of A..G         K = A * (T/B)^C * exp (D/T)
!
! # function for photolysis reactions
!#   function = 0 => ! constant independent of sza, K read from colum kn2rd
!#            Keff = kn2rd
!#   function = 1 => ! constant independent of sza
!#            Keff = A
!#   function = 2 =>) ! exponential function
!#            Keff = A * exp(B / coszen)
!#   function = 3 => ! powerfunction
!#            Keff = A * coszen ** B
!#   function = 4 => ! powerfunction but special for JO3 because of dependence on H2O and Air
!#            Keff = A * (coszen ** B) * D*[H2O] / (D*[H2O] + E*[M])
!#   function = any other number =>
!#            Keff = 1.
!#
! # reactions for nighttime chemistry
! #0.583E+0  R_57A  0    2    1.8e-11   110     1.0    1.0      1.0    1.0    1.0   NO + NO3    -> 2NO2
! #0.109E-5  R_58A  0    2    1.4e-13  -2470    1.0    1.0      1.0    1.0    1.0   NO2 + O3    -> NO3 + (O2)
! #0.567E-11 R_63Aa 0    2    2.5e-22   0.0     1.0    1.0      1.0    1.0    1.0   N2O5 + H2O  -> 2HNO3
! #1.0       R_63Ab 0    7    1.8e-39   1.0     0.0    0.0      1.0    1.0    1.0   N2O5 + 2H2O -> 2HNO3 + H2O
!
! #Photolysis
! #5.0667e-6 R_O3   1    1      1.0     1.0     1.0     1.0     1.0     1.0    1.0  O3         -> 2HO + (O2)
! #0.0167    R_NO2  1    1      1.0     1.0     1.0     1.0     1.0     1.0    1.0  NO2 + (O2) -> NO + O3
! #2.530e-6  R_O3   1    0      1.0     1.0     1.0     1.0     1.0     1.0    1.0  O3         -> 2HO + (O2)
! #8.330e-3  R_NO2  1    0      1.0     1.0     1.0     1.0     1.0     1.0    1.0  NO2 + (O2) -> NO + O3
!
!
implicit none
private
PUBLIC :: lchem, initchem,inputchem, twostep, PL_scheme, nchsp, firstchem, lastchem, RH, choffset
save

  ! namoptions
  integer tnor, firstchem, lastchem
  real dtchmovie,itermin
  real t_ref,q_ref,p_ref,h_ref
  logical lchem, ldiuvar,lchconst,lchmovie,lcloudKconst,lsegr

  real tnextwrite
  logical switch

  integer   mrpcc
  parameter (mrpcc = 20)

  integer,parameter :: NCCBA = 4   !Number Chemical Components Before Arrow !!!!!!!  4 is MAXIMUM !!!!!
  integer,parameter :: NCCAA = 8   !Number Chemical Components After Arrow  on change recompile !!!!!
  integer,parameter :: NNSPEC = 2*(NCCBA + NCCAA) - 1

  integer ,parameter :: numit = 2

  real, parameter :: MW_air = 28.97
  real, parameter :: MW_H2O = 18
  real, parameter :: Avogrado = 6.023e23
  real, parameter :: ppb  = 1.e9 ! convert to ppb
  integer PRODUCTION, LOSS
  parameter (PRODUCTION=1, LOSS=2)

  integer nr_raddep !number of photolysis reactions
  integer nchsp     !number of chemical species calculated from lastchem - firstchem
  integer choffset  !=firstchem -1 => offset for chemicals in sv0

  type RCdef
    character*6 rname
    integer raddep  ! 1 if reaction = radiation dependend
    real Kreact       !reaction konstant from input file
    real Keff         !to use with special circumstances cq with radiation and/or temperature depend reactions
    integer Kindex    !index to array with effective K due to clouds and temperature
    integer func1
    real A
    real B
    real C
    real D
    real E
    real F
    real G
  end type RCdef

  type (RCdef), allocatable :: RC(:) !tnor

  type,PUBLIC :: Form
    integer formula  !number of the formula to use
    integer r_nr    !reaction number, index to RC
    integer PorL    !0=> not active   1=>Production  2=>Loss
    real coef    !coefficient in formula
    integer comp1    !index to sv0
    integer exp1
    integer comp2    !index to sv0
    integer exp2
    integer comp3    !index to sv0
    integer exp3
    integer comp4    !index to sv0
    integer exp4
 end type Form

  type,PUBLIC :: Name_Number
    character (len=6) name  !name of chemical
    logical active      !active=1 else 0
    integer chem_number    !number (not really used)
    real atol
    real rtol
    integer nr_PL      !total number of reactions in which this chemical is used
    type (Form) PL(mrpcc)  !stucture holding the reaction components, reaction number etc
  end type Name_Number

  type (Name_Number), allocatable ::PL_scheme(:)   !(nchsp)

  integer, allocatable :: raddep_RCindex(:) !(nr_raddep)
  real, allocatable :: keff(:,:,:,:)    ! (i,j2,raddep_nr,kmax)
  real, allocatable :: keffT(:,:,:)    !(i2,j2,kefft_nr)
  real, allocatable :: keffT3D(:,:,:,:)  !(i2,j2,raddep_nr,kmax)
  real, allocatable :: atol(:),rtol(:)  ! nchsp
  real, allocatable :: rk1(:,:),rk2(:,:),rk(:,:),kreact(:,:)
  real, allocatable :: T_abs(:,:),convppb(:,:)

  real, allocatable :: ynew(:,:,:),yold(:,:,:),ysum(:,:,:)
  real, allocatable,target :: yl(:,:,:),yp(:,:,:)

  real, allocatable :: writearray(:,:)
  real*4, allocatable :: k3d(:,:,:,:)

  type, PUBLIC :: location
    character (len = 6) name
    integer   loc
  end type location

 type (location) :: INERT, PRODUC , O3, O1D, NO2, NO, NO3, N2O5, HNO3, RH, R, ISO, RO2, H2O2, HO2, HO, CO, CO2, H2O, NH3, H2SO4, CH2O, CH3O2, MVK


  real, allocatable :: kefftemp(:,:)  !(1:kmax,1:nr_raddep)

  real, allocatable :: segregation(:)          !segregation calculated (amount of reactions)
  real, allocatable :: segregation_vert(:,:)   !segregation calculated (amount of reactions,height)
  real, allocatable :: seg_conc_prod(:)        !average of multiplied concentrations (amount of reactions)
  real, allocatable :: seg_conc_prod_vert(:,:) !average of multiplied concentrations (amount of reactions,height)
  real, allocatable :: seg_conc(:,:,:)         !individual averaged concentration (max. input chemicals involved per reaction,height,amount of reactions)
  real, allocatable :: seg_conc_mult(:)        !multiplication of the relevant averaged concentrations (amount of reactions)
  real, allocatable :: seg_conc_mult_vert(:,:) !multiplication of the relevant averaged concentrations (amount of reactions,height)
  real, allocatable :: seg_concl(:,:)          !summation parameter for individual averaged concentration (max. amount of input chemicals involved,height)
  real, allocatable :: seg_conc_prodl(:)       !summation parameter for average of multiplied concentrations (height)
  logical, allocatable :: reaction_ev(:)       !switch to check whether the reaction is already evaluated (amount of reactions)

contains
!-----------------------------------------------------------------------------------------
SUBROUTINE initchem
  use modglobal, only : imax,jmax,i1,i2,ih, j1,j2,jh, k1, kmax, nsv, ifnamopt, fname_options, ifoutput, cexpnr,timeav_glob,btime
  use modmpi,    only : myid, mpi_logical, mpi_integer, my_real, comm3d, mpierr
  implicit none

  integer i, ierr

  namelist/NAMCHEM/ lchem, lcloudKconst, tnor, firstchem,lastchem,ldiuvar,h_ref,lchconst, t_ref, q_ref, p_ref,lchmovie, dtchmovie,lsegr

  itermin  = 1.e-6

  lchem    =.false.
  ldiuvar  = .false.
  h_ref    = 12.0
  lchconst = .false.
  t_ref    = 298.
  p_ref    = 100000.
  q_ref    = 5.e-3
  lchmovie = .false.
  dtchmovie= 60
  firstchem= 1
  lastchem = nsv
  nchsp    = nsv
  lcloudKconst  = .false.
  lsegr    = .false.

  if(myid==0) then
    open(ifnamopt,file=fname_options,status='old',iostat=ierr)
    read (ifnamopt,NAMCHEM,iostat=ierr)
    if (ierr > 0) then
      print *, 'Problem in namoptions NAMCHEM'
      print *, 'iostat error: ', ierr
      stop 'ERROR: Problem in namoptions NAMCHEM'
    endif
    write(6 ,NAMCHEM)
    close(ifnamopt)
  endif


  call MPI_BCAST(lchem     ,1,mpi_logical , 0,comm3d, mpierr)
  call MPI_BCAST(ldiuvar   ,1,mpi_logical , 0,comm3d, mpierr)
  call MPI_BCAST(lchconst  ,1,mpi_logical , 0,comm3d, mpierr)
  call MPI_BCAST(lchmovie  ,1,mpi_logical , 0,comm3d, mpierr)
  call MPI_BCAST(lsegr     ,1,mpi_logical , 0,comm3d, mpierr)
  call MPI_BCAST(lcloudKconst,1,mpi_logical ,0,comm3d,mpierr)
  call MPI_BCAST(tnor      ,1,mpi_integer , 0,comm3d, mpierr)
  call MPI_BCAST(firstchem ,1,mpi_integer , 0,comm3d, mpierr)
  call MPI_BCAST(lastchem  ,1,mpi_integer , 0,comm3d, mpierr)
  call MPI_BCAST(t_ref     ,1,MY_REAL     , 0,comm3d, mpierr)
  call MPI_BCAST(q_ref     ,1,MY_REAL     , 0,comm3d, mpierr)
  call MPI_BCAST(p_ref     ,1,MY_REAL     , 0,comm3d, mpierr)
  call MPI_BCAST(h_ref     ,1,MY_REAL     , 0,comm3d, mpierr)
  call MPI_BCAST(itermin   ,1,MY_REAL     , 0,comm3d, mpierr)
  call MPI_BCAST(dtchmovie ,1,MY_REAL     , 0,comm3d, mpierr)

  if (.not. (lchem)) return
  tnextwrite = timeav_glob-1e-3+btime
  switch = .false.
  nchsp = lastchem - firstchem  + 1
  choffset = firstchem - 1

  allocate(PL_scheme(nchsp))
  allocate (atol(nchsp), rtol(nchsp))

  allocate (ynew(2:i1,2:j1,nchsp),yold(2:i1,2:j1,nchsp),ysum(2:i1,2:j1,nchsp))
  allocate (yl(2:i1,2:j1,nchsp),yp(2:i1,2:j1,nchsp))
  allocate (rk1(2:i1,2:j1),rk2(2:i1,2:j1),rk(2:i1,2:j1),kreact(2:i1,2:j1))
  allocate (T_abs(2:i1,2:j1), convppb(2:i1,2:j1))

  PL_scheme(1)%name = '     '
  PL_scheme(1)%atol = 0.
  PL_scheme(1)%rtol = 0.
  PL_scheme(1)%chem_number = 0
  PL_scheme(1)%nr_PL = 0

  do i=1,mrpcc
    PL_scheme(1)%PL(i)%formula = 0
    PL_scheme(1)%PL(i)%r_nr = 0
    PL_scheme(1)%PL(i)%coef = 0.
    PL_scheme(1)%PL(i)%comp1 = 0
    PL_scheme(1)%PL(i)%exp1 = 0
    PL_scheme(1)%PL(i)%comp2 = 0
    PL_scheme(1)%PL(i)%exp2 = 0
    PL_scheme(1)%PL(i)%comp3 = 0
    PL_scheme(1)%PL(i)%exp3 = 0
    PL_scheme(1)%PL(i)%comp4 = 0
    PL_scheme(1)%PL(i)%exp4 = 0
  enddo

  do i=2, nchsp
    pl_scheme(i)=pl_scheme(1)
  enddo

  if(myid==0)then
    open (ifoutput,file='cloudstat.'//cexpnr,status='replace')
    close (ifoutput)
    open (ifoutput,file='keffs.'//cexpnr,status='replace')
    close (ifoutput)
    if (lsegr .EQV. .true.) then
      open (ifoutput,file='seg.'//cexpnr,status='replace') !File that describes segregation over whole Mixed Layer
      close(ifoutput)
      open (ifoutput,file='seg_h.'//cexpnr,status='replace') !File that describes segregation per height
      close(ifoutput)
    endif !Output for segregation
  endif

  call inputchem

end subroutine initchem

!----------------------------------------------------------
subroutine inputchem

!c***********************************************************************
!c
!c  Determine from the chemical equations in chem.inp.xxx
!c  the scheme to solve the chemical production and loss terms
!c
!c***********************************************************************
 use modglobal, only : imax,jmax,i1,i2,ih, j1,j2,jh, k1,kmax, nsv,cexpnr, ifoutput
 use modmpi,    only : myid
 implicit none

  integer i,j,k,l,react
  integer*2 number
  integer react_nr
  real reactconst,coefficient, fact(7)
  integer func1,raddep,nr_chemcomp,nr_active_chemicals
  integer keff_nr, keffT_nr
  character*11 spec(NNSPEC)
  character*255 line
  logical prod,found
  character (len=6) tempname
  character (len=6) name
  character (len=6) rname
  integer icoeff(4)
  character*30 formatstring

  type Chem
    real coeff
    character (len=6) name
    integer chem_nr
    integer index_sv0
  end type Chem

  type Reaction
    character*6  name
    real kr      !kn2rd
    integer RadDep   !reaction is radiation dependend
    integer keff_index
    integer Order  !orde of reaction
    integer nr_chem  !nr of chemicals in reaction (including non active species
    integer nr_chem_inp !nr of chem on input
    integer nr_chem_outp !nr of chem on output
    type (Chem) inp(NCCBA)
    type (Chem) outp(NCCAA)
  end type Reaction

  type (Name_Number),allocatable :: PL_temp(:)    !(nchsp)
  type (Reaction),allocatable    :: reactions(:)  !(totalnumberofreactions=tnor)
  character*6, allocatable       :: chem_name(:)  !(nchsp)


  if (.not. (lchem)) return

  react = 0
  number = 0
  keff_nr = 0
  keffT_nr = 0


  allocate (PL_temp(nchsp))
  allocate (reactions(tnor))
  allocate (RC(tnor))
  allocate (chem_name(nchsp))

  reactions(:)%kr = 0.0
  reactions(:)%RadDep = 0
  reactions(:)%keff_index = 0
  reactions(:)%order = 0
  reactions(:)%nr_chem = 0
  reactions(:)%nr_chem_inp = 0
  reactions(:)%nr_chem_outp= 0
  do i=1,4
    reactions(:)%inp(i)%Coeff = 0.
    reactions(:)%inp(i)%name = '     '
    reactions(:)%inp(i)%chem_nr = 0
    reactions(:)%inp(i)%index_sv0 = 999
    reactions(:)%outp(i)%Coeff = 0
    reactions(:)%outp(i)%name = '     '
    reactions(:)%outp(i)%chem_nr = 0
    reactions(:)%outp(i)%index_sv0 = 999
  enddo

  open (unit=10,file='chem.inp.'//cexpnr,err=100,status='old',form='formatted')

  do while(.true.)
    read(10,'(a)',err=100) line
    if (line(1:1)=='#') then
      if (myid == 0)   print *, trim(line)
    elseif (line(1:1) == '@') then
      call read_chem(chem_name)
    elseif (line(1:1) == '$') then
     !end of chemical reactions so we are done, jump out do while loop
      exit
    else
      number = number + 1
      react = react + 1
      if( react > tnor ) then
        write(6,*) 'Number of reactions is greater then tnor specified in namoptions',tnor
        STOP
      endif
      if (myid == 0) then
        write(*,'(i2,2x,a)') number,trim(line)
      endif

      do i=1,NNSPEC
        spec(i)='           '
      enddo

      read(line,*,end=300)reactconst,rname,raddep,func1,(fact(j),j=1,7),(spec(j),j=1,NNSPEC)
300   j=j-1

      if ((func1 == 6 .or. (raddep==1 .and. func1== 4)) .and. H2O%loc == 0) then
        write(*,*) 'Function 6 or 4 needs H2O and this is not specified as a chemical component'
        STOP
      endif

      !determine the number of chemical components
      i=1
      do while (len_trim(spec(i))>0 .and. i<NNSPEC)
        i=i+1
      end do
      !   if (myid==0) print *,'aantal componenten=', i,(i+1)/2,(spec(j),j=1,7)

      nr_chemcomp = (i+1)/2

      prod = .false.
      l=0
      reactions(react)%kr     = reactconst
      reactions(react)%name   = rname
      reactions(react)%nr_chem = nr_chemcomp

      if (raddep == 1) then
        keff_nr = keff_nr +1
        reactions(react)%RadDep = raddep
        reactions(react)%keff_index = keff_nr
        RC(react)%Kindex  = keff_nr
      else
        keffT_nr = keffT_nr +1
        reactions(react)%keff_index = keffT_nr
        RC(react)%Kindex  = keffT_nr
      end if

      RC(react)%Kreact  = reactconst
      RC(react)%Keff    = reactconst
      RC(react)%rname   = rname
      RC(react)%RadDep  = raddep
      RC(react)%func1   = func1
      RC(react)%A = fact(1)
      RC(react)%B = fact(2)
      RC(react)%C = fact(3)
      RC(react)%D = fact(4)
      RC(react)%E = fact(5)
      RC(react)%F = fact(6)
      RC(react)%G = fact(7)


!analyze reaction scheme,determine chem species and location in scalar.inp
!and store in reactions(react)

      do j=1, 2 * nr_chemcomp - 1
        select case (spec(j))
        case ('+          ')
!         print *,'found +'
        case ('->         ')
!         print *,'found ->'
          prod = .true.
          reactions(react)%nr_chem_inp = l
          l=0
        case default
          l=l+1
          if ( spec(j)(1:1) == '(' ) then
            !non active species forget it
            l=l-1
          else
            do i=1,len(spec(j))
              if( spec(j)(i:i) .GT. '@' ) then
                !starting name chemical component
                tempname = spec(j)(i:len(spec(j)))
                if( i == 1 ) then !nothing before chem comp
                  coefficient = 1.
                else !we have numbers before
                  read( spec(j)(1:i-1),*)coefficient
                endif
                if( (prod .eqv. .false.) .and. ((coefficient +.0005)< 1.))then
                  write(*,*) 'Sorry, coefficient on input should be a multiply of 1'
                  STOP
                endif
                if( (prod .eqv. .false.) .and. (coefficient/int(coefficient) > 1.005) ) then
                  write(*,*) 'Sorry, coefficient on input should be a multiply of 1 found:',spec(j)
                  STOP
                endif
                exit
              else
                if (i >= len(spec(j)) )then
                  write(*,*)'Probably space between coefficient and chemical component'
                  write(*,*) 'look between ',spec(j),spec(j+1)
                  STOP
                endif
              endif
            enddo

            !find index in sv0
            i=1
            do while(tempname /= chem_name(i) )
              i= i+1
              if (i > nchsp) then
                if (myid == 0) print *,'Name ',tempname, 'NOT FOUND in speciesline after @'
                STOP
              end if
            end do

            if (prod .EQV. .false.) then
              reactions(react)%inp(l)%name = tempname
              reactions(react)%inp(l)%coeff = coefficient
              reactions(react)%inp(l)%index_sv0 = i
            else
              reactions(react)%outp(l)%name = tempname
              reactions(react)%outp(l)%coeff = coefficient
              reactions(react)%outp(l)%index_sv0 = i
            endif
          endif
        end select
      enddo !j=1,2*nr_chemcomp - 1
      reactions(react)%nr_chem_outp = l
    endif
  end do ! while(1)

!total number of reactions is equal to react

  if (myid == 0)  print *, 'Total number of reactions is',react,'of which', keff_nr,'is/are radiation dependent'

  if (tnor > react ) tnor = react

  nr_raddep = keff_nr

  allocate (raddep_RCindex(keff_nr))
  allocate (keff(2:i1,2:j1,keff_nr,kmax))
  allocate (kefftemp(kmax,keff_nr))
  allocate (keffT(2:i1,2:j1,keffT_nr))
  allocate (writearray(k1,tnor+2))

!  if (lsegrk .EQV. .true.) then                    ! As of yet unused expansion of segregation
!    allocate (keffT3D(2:i1,2:j1,keffT_nr,kmax))
!  endif

  if (lsegr .EQV. .true.) then
    allocate (segregation(tnor))
    allocate (segregation_vert(tnor,kmax))
    allocate (seg_conc_prod(tnor))
    allocate (seg_conc_prod_vert(tnor,kmax))
    allocate (seg_conc(NCCBA,kmax,tnor))
    allocate (seg_conc_mult(tnor))
    allocate (seg_conc_mult_vert(tnor,kmax))
    allocate (seg_concl(NCCBA,kmax))
    allocate (seg_conc_prodl(kmax))
    allocate (reaction_ev(tnor))

    open (ifoutput,file='seg.'//cexpnr,position='append') !File that describes segregation over whole Mixed Layer
      write(formatstring,'(a,i3,a)') '(a9,1x,a7,',tnor,'(2x,a5,i3.3,a4))'
      write(ifoutput,formatstring) '#Time [s]','z_i [m]',('IsegR',i,' [-]',i=1,tnor)
  endif

  !fill raddep_RCindex and keffT
  do i=1,react
    if (reactions(i)%raddep == 1) then
      raddep_RCindex(reactions(i)%keff_index) = i
    else
      keffT(:,:,reactions(i)%keff_index) = reactions(i)%kr
    endif
  enddo

  !make a list of chemical species and in which reaction number it is formed and/or losst
  k=0
  do i=1,react
    do j=1,reactions(i)%nr_chem_inp ! look only on input side of reaction
      name = reactions(i)%inp(j)%name
      found = .false.
      do l=1,k
        if(name == PL_scheme(l)%name) then
          found = .true.
          reactions(i)%inp(j)%chem_nr = l    !put chem component number in reaction
          PL_scheme(l)%nr_PL = PL_scheme(l)%nr_PL +1 !count number of reactions
          if ( PL_scheme(l)%nr_PL > mrpcc ) then
            print *, 'mrpcc to low, increase mrpcc in modchem'
            stop
          end if
          PL_scheme(l)%PL(PL_scheme(l)%nr_PL)%r_nr = i   !store reaction number index to RC
          PL_scheme(l)%PL(PL_scheme(l)%nr_PL)%PorL = 2   !this is a loss reaction for this component
          exit
        end if
      enddo
      if (found .EQV. .false.) then
        k=k+1
        PL_scheme(k)%name=name
        PL_scheme(k)%chem_number = k
        reactions(i)%inp(j)%chem_nr = k
        PL_scheme(l)%nr_PL = PL_scheme(l)%nr_PL +1
        if ( PL_scheme(l)%nr_PL > mrpcc ) then
          print *, 'mrpcc to low, increase mrpcc in modchem'
          stop
        end if
        PL_scheme(l)%PL(PL_scheme(l)%nr_PL)%r_nr = i   !store reaction number
        PL_scheme(l)%PL(PL_scheme(l)%nr_PL)%PorL = 2   !this is a loss reaction for this component
      endif
    enddo
  enddo

  do i=1,react
    do j=1,reactions(i)%nr_chem_outp  !this is for after the ->
      name = reactions(i)%outp(j)%name
      found = .false.
      do l=1,k
        if(name == PL_scheme(l)%name) then
          found = .true.
          reactions(i)%outp(j)%chem_nr = l
          PL_scheme(l)%nr_PL = PL_scheme(l)%nr_PL +1
          if ( PL_scheme( l)%nr_PL > mrpcc ) then
            print *, 'mrpcc to low, increase mrpcc in modchem'
            stop
          end if
          PL_scheme(l)%PL(PL_scheme(l)%nr_PL)%r_nr = i   !store reaction number
          PL_scheme(l)%PL(PL_scheme(l)%nr_PL)%PorL = 1   !this is a production reaction for this component
          exit
        end if
      enddo
      if (found .EQV. .false.) then
        k=k+1
        PL_scheme(k)%name=name
        PL_scheme(k)%chem_number = k
        reactions(i)%outp(j)%chem_nr = k
        PL_scheme(l)%nr_PL = PL_scheme(l)%nr_PL +1
        if ( PL_scheme(l)%nr_PL > mrpcc ) then
          print *, 'mrpcc to low, increase mrpcc in modchem'
          stop
        end if
        PL_scheme(l)%PL(PL_scheme(l)%nr_PL)%r_nr = i   !store reaction number
        PL_scheme(l)%PL(PL_scheme(l)%nr_PL)%PorL = 1   !this is a production reaction for this component
      endif
    enddo !j=1,reactions(i)%nr_chem_outp
  enddo !i=1,react


  nr_active_chemicals = k

  if (nr_active_chemicals < nchsp ) then
    if (myid == 0)  print *, 'WARNING: More chemicals specified in @ line then actually used.', nr_active_chemicals,' <', nchsp
  endif

  !Determine from the reactions which formula to use for all the production and loss reactions
  !first do all reactions on the production side

  do i=1,nr_active_chemicals        !********************* misschien nchsp
    do j=1,PL_scheme(i)%nr_PL
      react_nr = PL_scheme(i)%PL(j)%r_nr
      if( PL_scheme(i)%PL(j)%PorL == PRODUCTION) then     !this is a PRODUCTION
        do k=1,reactions(react_nr)%nr_chem_outp
          if( reactions(react_nr)%outp(k)%name == PL_scheme(i)%name) then
            select case(reactions(react_nr)%nr_chem_inp)  !left of arrow 1 reactant
              case (1)
              icoeff(1) = int(reactions(react_nr)%inp(1)%coeff +0.05)
                 select case (icoeff(1))
                   case (1)
                     PL_scheme(i)%PL(j)%formula = 1
                     PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                   case (2)
                     PL_scheme(i)%PL(j)%formula = 2
                     PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                     PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(1)%index_sv0
                   case default
                     PL_scheme(i)%PL(j)%formula = 3
                     PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                     PL_scheme(i)%PL(j)%exp1  = icoeff(1)
                  end select
                  PL_scheme(i)%PL(j)%coef =  reactions(react_nr)%outp(k)%coeff
              case(2)  !there are 2 reacting species
              icoeff(1) = int(reactions(react_nr)%inp(1)%coeff +0.05)
              icoeff(2) = int(reactions(react_nr)%inp(2)%coeff +0.05)
                if((icoeff(1) == 1) .AND. (icoeff(2) == 1) ) then
                  PL_scheme(i)%PL(j)%formula = 2
                  PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                  PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
                else
                  PL_scheme(i)%PL(j)%formula = 4
                  PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                  PL_scheme(i)%PL(j)%exp1  = icoeff(1)
                  PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
                  PL_scheme(i)%PL(j)%exp2  = icoeff(2)
                endif
                PL_scheme(i)%PL(j)%coef =  reactions(react_nr)%outp(k)%coeff
              case (3) !! there are 3 reacting species
               icoeff(1) = int(reactions(react_nr)%inp(1)%coeff +0.05)
               icoeff(2) = int(reactions(react_nr)%inp(2)%coeff +0.05)
               icoeff(3) = int(reactions(react_nr)%inp(3)%coeff +0.05)
               if( (icoeff(1)==1) .and. (icoeff(2) == 1) .and. (icoeff(3) == 1)) then
                  ! we don't need exponents keep it simple
                  PL_scheme(i)%PL(j)%formula = 5
                  PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                  PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
                  PL_scheme(i)%PL(j)%comp3 = reactions(react_nr)%inp(3)%index_sv0
                else
                  PL_scheme(i)%PL(j)%formula = 6
                  PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                  PL_scheme(i)%PL(j)%exp1  = icoeff(1)
                  PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
                  PL_scheme(i)%PL(j)%exp2  = icoeff(2)
                  PL_scheme(i)%PL(j)%comp3 = reactions(react_nr)%inp(3)%index_sv0
                  PL_scheme(i)%PL(j)%exp3  = icoeff(3)
                endif
                PL_scheme(i)%PL(j)%coef =  reactions(react_nr)%outp(k)%coeff
              case (4)
               icoeff(1) = int(reactions(react_nr)%inp(1)%coeff +0.05)
               icoeff(2) = int(reactions(react_nr)%inp(2)%coeff +0.05)
               icoeff(3) = int(reactions(react_nr)%inp(3)%coeff +0.05)
               icoeff(4) = int(reactions(react_nr)%inp(4)%coeff +0.05)
               PL_scheme(i)%PL(j)%formula = 7
               PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
               PL_scheme(i)%PL(j)%exp1  = icoeff(1)
               PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
               PL_scheme(i)%PL(j)%exp2  = icoeff(2)
               PL_scheme(i)%PL(j)%comp3 = reactions(react_nr)%inp(3)%index_sv0
               PL_scheme(i)%PL(j)%exp3  = icoeff(3)
               PL_scheme(i)%PL(j)%comp4 = reactions(react_nr)%inp(4)%index_sv0
               PL_scheme(i)%PL(j)%exp4  = icoeff(4)
               PL_scheme(i)%PL(j)%coef =  reactions(react_nr)%outp(k)%coeff
            end select
          endif
        enddo ! k=1,reactions(react_nr)%nr_chem_outp
      endif !( PL_scheme(i)%PL(j)%PorL == 1)
    enddo !j=1,PL_scheme(i)%nr_PL
  enddo !1,nr_active_chemicals

  !do all reactions on the loss side
  do i=1,nr_active_chemicals
    do j=1,PL_scheme(i)%nr_PL
      react_nr = PL_scheme(i)%PL(j)%r_nr
      if( PL_scheme(i)%PL(j)%PorL == LOSS) then     !This is LOSS
      do k=1,reactions(react_nr)%nr_chem_inp
        icoeff(1) = int(reactions(react_nr)%inp(1)%coeff +0.05)
        icoeff(2) = int(reactions(react_nr)%inp(2)%coeff +0.05)
        icoeff(3) = int(reactions(react_nr)%inp(3)%coeff +0.05)
        icoeff(4) = int(reactions(react_nr)%inp(4)%coeff +0.05)
        select case(reactions(react_nr)%nr_chem_inp)
        case (1) !the loss comp is the only reactant
          select case (icoeff(1))
            case (1)
              PL_scheme(i)%PL(j)%formula = 0
            case (2)
              PL_scheme(i)%PL(j)%formula = 1
              PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
            case (3)
              PL_scheme(i)%PL(j)%formula = 3
              PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
              PL_scheme(i)%PL(j)%exp1  = icoeff(1) - 1
          end select
          PL_scheme(i)%PL(j)%coef = reactions(react_nr)%inp(1)%coeff
        case (2) ! we have 2 components in which one is current species
          if( reactions(react_nr)%inp(k)%name == PL_scheme(i)%name) then ! current selected
            select case (k)
            case(1)
              if( (icoeff(1) ==1) .and. (icoeff(2) == 1) ) then
                PL_scheme(i)%PL(j)%formula = 1
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(2)%index_sv0
              else if((icoeff(1) ==1) .and. (icoeff(2) > 1) ) then
                PL_scheme(i)%PL(j)%formula = 3
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(2)%index_sv0
                PL_scheme(i)%PL(j)%exp1  = icoeff(2)
              else if((icoeff(1) ==2) .and. (icoeff(2) == 1) ) then
                PL_scheme(i)%PL(j)%formula = 2
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
              else
                PL_scheme(i)%PL(j)%formula = 4
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                PL_scheme(i)%PL(j)%exp1  = icoeff(1) - 1
                PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
                PL_scheme(i)%PL(j)%exp2  = icoeff(2)
              endif
              PL_scheme(i)%PL(j)%coef = reactions(react_nr)%inp(1)%coeff
            case(2)
              if( (icoeff(1) ==1) .and. (icoeff(2) == 1) ) then
                PL_scheme(i)%PL(j)%formula = 1
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
              else if((icoeff(2) ==1) .and. (icoeff(1) > 1)) then
                PL_scheme(i)%PL(j)%formula = 3
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                PL_scheme(i)%PL(j)%exp1  = icoeff(1)
              else if((icoeff(2) ==2) .and. (icoeff(1) == 1) ) then
                PL_scheme(i)%PL(j)%formula = 2
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
              else
                PL_scheme(i)%PL(j)%formula = 4
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                PL_scheme(i)%PL(j)%exp1  = icoeff(1)
                PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
                PL_scheme(i)%PL(j)%exp2  = icoeff(2) - 1
              endif
              PL_scheme(i)%PL(j)%coef = reactions(react_nr)%inp(2)%coeff
            end select
          endif
        case (3) !we have 3 components on input
          if( reactions(react_nr)%inp(k)%name == PL_scheme(i)%name) then ! current selected
            select case (k)
              case (1)
                PL_scheme(i)%PL(j)%formula = 6
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(2)%index_sv0
                PL_scheme(i)%PL(j)%exp1  = icoeff(2)
                PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(3)%index_sv0
                PL_scheme(i)%PL(j)%exp2  = icoeff(3)
                PL_scheme(i)%PL(j)%coef = reactions(react_nr)%inp(1)%coeff
                if ( icoeff(1) == 1 ) then
                  PL_scheme(i)%PL(j)%formula = 4
                else
                  PL_scheme(i)%PL(j)%comp3 = reactions(react_nr)%inp(1)%index_sv0
                  PL_scheme(i)%PL(j)%exp3  = icoeff(1) - 1
                endif
              case (2)
                PL_scheme(i)%PL(j)%formula = 6
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                PL_scheme(i)%PL(j)%exp1  = icoeff(1)
                PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(3)%index_sv0
                PL_scheme(i)%PL(j)%exp2  = icoeff(3)
                PL_scheme(i)%PL(j)%coef = reactions(react_nr)%inp(2)%coeff
                if ( icoeff(2) == 1 ) then
                  PL_scheme(i)%PL(j)%formula = 4
                else
                  PL_scheme(i)%PL(j)%comp3 = reactions(react_nr)%inp(2)%index_sv0
                  PL_scheme(i)%PL(j)%exp3  = icoeff(2) -1
                endif
              case (3)
                PL_scheme(i)%PL(j)%formula = 6
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                PL_scheme(i)%PL(j)%exp1  = icoeff(1)
                PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
                PL_scheme(i)%PL(j)%exp2  = icoeff(2)
                PL_scheme(i)%PL(j)%coef = reactions(react_nr)%inp(3)%coeff
                if ( icoeff(3) == 1 ) then
                  PL_scheme(i)%PL(j)%formula = 4
                else
                  PL_scheme(i)%PL(j)%comp3 = reactions(react_nr)%inp(3)%index_sv0
                  PL_scheme(i)%PL(j)%exp3  = icoeff(3) -1
                endif
            end select
          endif
        case (4) !we have 4 components on input
          PL_scheme(i)%PL(j)%coef = 1  !????????????
          if( reactions(react_nr)%inp(k)%name == PL_scheme(i)%name) then ! current selected
            PL_scheme(i)%PL(j)%formula = 7
            PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
            PL_scheme(i)%PL(j)%exp1  = icoeff(1)
            PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
            PL_scheme(i)%PL(j)%exp2  = icoeff(2)
            PL_scheme(i)%PL(j)%comp3 = reactions(react_nr)%inp(3)%index_sv0
            PL_scheme(i)%PL(j)%exp3  = icoeff(3)
            PL_scheme(i)%PL(j)%comp4 = reactions(react_nr)%inp(4)%index_sv0
            PL_scheme(i)%PL(j)%exp4  = icoeff(4)
            select case(k)
              case (1)
                PL_scheme(i)%PL(j)%exp1  = icoeff(1) - 1
                PL_scheme(i)%PL(j)%coef = reactions(react_nr)%inp(1)%coeff
              case (2)
                PL_scheme(i)%PL(j)%exp2  = icoeff(2) - 1
                PL_scheme(i)%PL(j)%coef = reactions(react_nr)%inp(2)%coeff
              case (3)
                PL_scheme(i)%PL(j)%exp3  = icoeff(3) - 1
                PL_scheme(i)%PL(j)%coef = reactions(react_nr)%inp(3)%coeff
              case (4)
                PL_scheme(i)%PL(j)%exp4  = icoeff(4) - 1
                PL_scheme(i)%PL(j)%coef = reactions(react_nr)%inp(4)%coeff
            end select
          endif
        end select
      enddo
      endif
    enddo
  enddo

  !fill PL_scheme with rtol and atol
  do i=1,nchsp
    do j=1,nchsp
      if (PL_scheme(i)%name == chem_name(j)) then
        PL_scheme(i)%atol = atol(j)
        PL_scheme(i)%rtol = rtol(j)
        !we don't need the org chem_number now use it as index to sv0
        PL_scheme(i)%chem_number = j
        exit
      end if
    end do
  end do

  ! the order of chemicals is in the order in which they appear in the chem reactions but we want them
  ! in the order of the @ line which is stored in chem_name()
  do i=1,nchsp
    found = .false.
    do j=1,nchsp
      if (chem_name(i) == PL_scheme(j)%name ) then
        PL_temp(i)= PL_scheme(j)
        found = .true.
        exit
      endif
    enddo
    if (found .EQV. .false.) then
      PL_temp(i)%name = chem_name(i)
    endif
  enddo

  PL_scheme = PL_temp

  ! check which chemicals are really used
  do i=1, nchsp
    if (PL_scheme(i)%nr_PL == 0 ) then
      PL_scheme(i)%active = .FALSE.
    else
      PL_scheme(i)%active = .TRUE.
    endif
  enddo

  call write_chem_scheme(chem_name)

  deallocate(PL_temp, reactions)

RETURN

100 if (myid == 0)  then
      print *,'error in inputchem'
      STOP
    ENDIF


end subroutine inputchem

!-----------------------------------------------------------------------------------------

SUBROUTINE read_chem(chem_name)
  use modglobal, only : nsv
  use modmpi,    only : myid

  implicit none

  character*6, dimension(nchsp) ::chem_name
  character*255 scalarline
  integer i,j

  INERT%name  = 'INERT'
  PRODUC%name = 'PRODUC'
  O3%name     = 'O3'
  NO%name     = 'NO'
  NO2%name    = 'NO2'
  NO3%name    = 'NO3'
  N2O5%name   = 'N2O5'
  HNO3%name   = 'HNO3'
  HO2%name    = 'HO2'
  HO%name     = 'HO'
  H2O2%name   = 'H2O2'
  H2O%name    = 'H2O'
  CO%name     = 'CO'
  CO2%name    = 'CO2'
  RH%name     = 'RH'
  R%name      = 'R'
  NH3%name    = 'NH3'
  H2SO4%name  = 'H2SO4'
  ISO%name    = 'ISO'

  !chem species and atol and rtol
  read(10,'(a)',err=100)scalarline

  read(scalarline,*)(chem_name(j),j=1,nchsp)

  do i=1, nchsp
    if (O3%name    == chem_name(i)) then ; O3%loc   = i;  cycle; endif
    if (NO%name    == chem_name(i)) then ; NO%loc   = i;  cycle; endif
    if (NO2%name   == chem_name(i)) then ; NO2%loc  = i;  cycle; endif
    if (NO3%name   == chem_name(i)) then ; NO3%loc  = i;  cycle; endif
    if (N2O5%name  == chem_name(i)) then ; N2O5%loc = i;  cycle; endif
    if (HNO3%name  == chem_name(i)) then ; HNO3%loc = i;  cycle; endif
    if (HO%name    == chem_name(i)) then ; HO%loc   = i;  cycle; endif
    if (HO2%name   == chem_name(i)) then ; HO2%loc  = i;  cycle; endif
    if (H2O2%name  == chem_name(i)) then ; H2O2%loc = i;  cycle; endif
    if (H2O%name   == chem_name(i)) then ; H2O%loc  = i;  cycle; endif
    if (CO%name    == chem_name(i)) then ; CO%loc   = i;  cycle; endif
    if (CO2%name   == chem_name(i)) then ; CO2%loc  = i;  cycle; endif
    if (RH%name    == chem_name(i)) then ; RH%loc   = i;  cycle; endif
    if (R%name     == chem_name(i)) then ; R%loc    = i;  cycle; endif
    if (ISO%name   == chem_name(i)) then ; ISO%loc    = i;  cycle; endif
    if (NH3%name   == chem_name(i)) then ; NH3%loc  = i;  cycle; endif
    if (H2SO4%name == chem_name(i)) then ; H2SO4%loc  = i;  cycle; endif
    if (INERT%name == chem_name(i)) then ; INERT%loc  = i;  cycle; endif
    if (PRODUC%name== chem_name(i)) then ; PRODUC%loc  = i;  cycle; endif
  enddo
! the above loc gives the location of a chemical component relative to the first chemical position
! in SV0. The loc of the first chemical is always 1.  If you need the absolute
! position in SV0 you have to add CHOFFSET to the above loc position

  read(10,'(a)',err=100)scalarline
  read(scalarline,*)(atol(j),j=1,nchsp)

  read(10,'(a)',err=100)scalarline
  read(scalarline,*)(rtol(j),j=1,nchsp)

  do i=1,nchsp
    if (myid == 0)  print*, i,'  ',chem_name(i),atol(i),rtol(i)
  enddo

  print *,' '

RETURN

100 print *, 'error in reading chem species in inputchem'
    STOP


end subroutine read_chem


!-----------------------------------------------------------------------------------------


SUBROUTINE twostep()     !(t,te,y)   (timee, timee+dt, sv0)
use modglobal, only : rk3step,timee,timeav_glob
use modfields, only: svm
use modmpi, only: myid
implicit none

  if (.not. (lchem)) return

  if (rk3step/=3) return
  if(timee==0) return

  !!!! We only use the chemistry scalars in svm,
  !!!! in twostep2 we use them as y(:,:,:,1:nchsp)
  !!!! they may be starting at XX but we acces them with index 1 to nchsp
  call twostep2(svm(:,:,:,firstchem:lastchem))
  if (timee >= tnextwrite ) then
    tnextwrite = tnextwrite + timeav_glob
  endif
end subroutine twostep

!-----------------------------------------------------------------------------------------


SUBROUTINE twostep2(y)
!c
!c-----------------------------------------------------------------|
!c                                                                 |
!c*** *twostep*  chewmical solver                                  |
!c                                                                 |
!c     Jordi Vila    WUR           16/08/2004                      |
!c                                                                 |
!c     purpose.                                                    |
!c     --------                                                    |
!c                                                                 |
!c     It calculates numerically the chemical                      |
!c     term of the conservation equation of the scalars.           |
!c     The solver is based on the TWOSTEP (Verwer et al., 1996,    |
!c     Atmospheric Environment 30, 49-58). The method is based on  |
!c     the implicit, second-order, two-step backward differential  |
!c     formula.                                                    |
!c                                                                 |
!c**   interface.                                                  |
!c     ----------                                                  |
!c                                                                 |
!c    *twostep* is called from *program*                           |
!c                                                                 |
!c-----------------------------------------------------------------|
!c
use modglobal, only : ih,i1,jh,j1,i2,j2,k1,kmax, nsv, xtime, timee,dt,timeav_glob, xday,xlat,xlon,zf,dzf,ifoutput,cexpnr,dz,rslabs
use modfields, only : qt0
use modmpi, only: comm3d, mpierr,mpi_max,mpi_min,mpi_sum,my_real,mpi_real,myid,cmyid,nprocs
use modtimestat, only: we, zi, ziold, calcblheight

implicit none

!!!!! we acces the chemicals from 1 to nchsp so we are independend of other scalars in svm
  real    y(2-ih:i1+ih,2-jh:j1+jh,k1,1:nchsp)

  real, allocatable :: ybegin(:,:,:,:)
  real, allocatable :: writearrayg(:,:)
  real    t, te
  integer nfcn,naccpt,nrejec,nstart,i,j,n,pl,it_i,it_j,it_k,it_l
  real    tb, kdt, dtmin, kdtmax, ytol
  real    ratio,dtold,a1,a2,c,cp1,dtg,errlte,dyc
  logical accept,failer,restart
  real    kdtl, errltel
  character*40 formatstring
  real    epsilon
  parameter (epsilon=1.0e-10)
  !parameter (dtmin=.2)    !orgineel 1.e-6)
  integer k_zi,Rnr
  integer nsecs, nhrs, nminut
  real    rest_zi, zi_temp, we_temp

  !c Initialization of logical and counters
  t = timee
  te = timee + dt

  dtmin   = itermin
  naccpt  = 0
  nrejec  = 0
  nfcn    = 0
  nstart  = 0
  tb      = t
  kdtmax  = te-t

  if (H2O%loc /= 0) then
    if (lchconst .EQV. .true.) then
      y(2:i1,2:j1,:,H2O%loc) = q_ref * MW_air / MW_h2o * ppb
    else
      y(2:i1,2:j1,:,H2O%loc) = qt0(2:i1,2:j1,:) * MW_air / MW_h2o * ppb
    endif
  endif

  if(lchmovie .and. mod(timee,dtchmovie)==0 ) then
    allocate(k3d(2:i1,2:j1,kmax,tnor+2))  !2 extra for T_abs and convppb
    allocate(ybegin(2-ih:i1+ih,2-jh:j1+jh,k1,nchsp))
    ybegin = y
  endif

  !c Chemical solver
  CALL ratech

  do pl=1,kmax

    !c Initialization of logicals and counters
    failer  =.false.
    restart =.false.
    accept  =.true.
    nfcn = 0
    nstart = 0
    !c Initial stepsize computation.
    if (dtmin.eq.kdtmax) then
      write(6,*) 'dtmin.eq.kdtmax'
      stop
    endif

    call calc_K(pl)

    ynew(:,:,1:nchsp) = y(2:i1,2:j1, pl,1:nchsp)
    ysum(:,:,1:nchsp)  = y(2:i1,2:j1, pl,1:nchsp)
    call ITER (0.0,pl)
    nfcn=nfcn+1
    kdtl=te-t

    !  if n points to H2O or INERT skip calculations
    do n=1,nchsp
      if (n == H2O%loc .or. (PL_scheme(n)%active .eqv. .false.)) cycle
      do j=2,j1
        do i=2,i1
          dyc=yp(i,j,n)-y(i,j, pl, n)*yl(i,j, n)
          if (dyc >epsilon) then
            ytol=atol(n)+rtol(n)*abs(y(i,j, pl, n))
            kdtl=min(kdtl,ytol/abs(dyc))
          endif
        enddo
      enddo
    enddo

    call MPI_ALLREDUCE(kdtl, kdt, 1, MY_REAL,MPI_MIN,  comm3d, mpierr)

25  nstart=nstart+1

    if (restart) kdt=kdt/10.0
    restart=.true.
    kdt=max(dtmin,min(kdt,kdtmax))
    call FIT(t,te,kdt)
    kdt=min(kdt,(te-t)/2)
    ! sta=kdt

    !c The starting step is carried out, using the implicit Euler method.
    ynew(:,:,1:nchsp)  = y(2:i1,2:j1, pl,1:nchsp)
    yold(:,:,1:nchsp)  = y(2:i1,2:j1, pl,1:nchsp)
    ysum(:,:,1:nchsp)  = y(2:i1,2:j1, pl,1:nchsp)

    do i=1,numit
      call ITER(kdt,pl)
      nfcn=nfcn+1
    enddo

    naccpt=naccpt+1
    t=t+kdt
    y(2:i1,2:j1,pl,1:nchsp)=ynew(:,:,1:nchsp)

    !c Subsequent steps are carried out with the two-step BDF method.
    dtold = kdt
    ratio = 1.0

60  continue

    c   = 1.0/ratio
    cp1 = c+1.0
    a1  = ((c+1.0)**2)/(c*c+2.0*c)
    a2  =-1.0/(c*c+2.0*c)
    dtg = kdt*(1.0+c)/(2.0+c)

    ysum(:,:,1:nchsp) = a1*y(2:i1,2:j1,pl,1:nchsp)+a2*yold(:,:,1:nchsp)
    ynew(:,:,1:nchsp) = max(0.0,y(2:i1,2:j1,pl,1:nchsp)+ ratio*(y(2:i1,2:j1,pl,1:nchsp)-yold(:,:,1:nchsp)))

    do i=1,numit
      call ITER(dtg,pl)
      nfcn=nfcn+1
    enddo

    !c Stepsize control.
    errltel=0.0
    do n=1,nchsp
      if (n == H2O%loc .or. n==INERT%loc) cycle  !we aren't interested in the chemistry of water and is kept constant
      do j=2,j1
        do i=2,i1
          ytol=atol(n)+rtol(n)*abs(y(i,j,pl,n))
          errltel=max(errltel,abs(c*ynew(i,j,n)-cp1* y(i,j,pl,n)+yold(i,j,n))/ytol)
        enddo
      enddo
    enddo

    call MPI_ALLREDUCE(errltel, errlte, 1,MY_REAL,MPI_MAX, comm3d, mpierr)

    errlte=2.0*errlte/(c+c*c)
    call NEWDT(t,te,kdt,dtold,ratio,errlte,accept,dtmin,kdtmax)
    !c Check if step has been accepted.
    if (accept) then
      failer  = .false.
      restart = .false.
      naccpt  = naccpt+1

      yold(:,:,1:nchsp) = y(2:i1,2:j1, pl,1:nchsp)
      y(2:i1,2:j1, pl,1:nchsp) = ynew(:,:,1:nchsp)

      t = t+dtold

      if (t.ge.te) goto 120
        goto 60
    endif

    !c A restart check is carried out.kmax
    if (failer) then
      failer = .false.
      nrejec = nrejec+1
      naccpt = naccpt-1
      t=t-dtold
      y(2:i1,2:j1, pl,1:nchsp) = yold(:,:,1:nchsp)

      goto 25
    endif

    !c Here the step has been rejected.
    nrejec=nrejec+1
    failer=.true.
    goto 60

120 t=tb
  enddo  !pl=1,kmax


  if (lsegr .EQV. .true.) then
    if(timee>=tnextwrite) then
      we_temp = we
      zi_temp = ziold
      call calcblheight
      we      = we_temp
      ziold   = zi_temp
      k_zi    = int(zi/dz)
      rest_zi = (zi/dz)-k_zi

      reaction_ev = .FALSE.
      do n=1,nchsp
        if (PL_scheme(n)%active .EQV. .TRUE.) then
          do j=1,PL_scheme(n)%nr_PL
          if(PL_scheme(n)%PL(j)%PorL == PRODUCTION) then !Only look at the reaction formulas from the perspective of produced chemicals
            Rnr = PL_scheme(n)%PL(j)%r_nr  !Number of the reaction evaluated
            if (reaction_ev(Rnr) .EQV. .FALSE.) then !Only look at reactions that are not evaluated yet
              seg_concl = 0.0
              seg_conc_prodl = 0.0
              select case (PL_scheme(n)%PL(j)%formula) ! Determine the chemical equation type

              case (0) ! No input chemicals, so segregation is always 0 everywhere
                seg_conc_prod      = 1.0
                seg_conc_mult      = 1.0
                seg_conc_prod_vert = 1.0
                seg_conc_mult_vert = 1.0

              case (1) ! A -> ...
                do it_k=1,kmax
                  do it_j=2,j1
                    do it_i=2,i1
                      seg_concl(1,it_k)    = seg_concl(1,it_k) + y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp1)
                      seg_conc_prodl(it_k) = seg_conc_prodl(it_k) + y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp1)
                    enddo
                  enddo
                enddo
                do it_l=1,1
                  call MPI_ALLREDUCE(seg_concl(it_l,:) ,seg_conc(it_l,:,Rnr) ,kmax ,MY_REAL ,MPI_SUM ,comm3d ,mpierr)
                enddo
                call MPI_ALLREDUCE(seg_conc_prodl(:) ,seg_conc_prod_vert(Rnr,:) ,kmax ,MY_REAL ,MPI_SUM ,comm3d ,mpierr)
                seg_conc_mult_vert(Rnr,:) = seg_conc(1,:,Rnr)
                seg_conc_prod(Rnr) = sum(seg_conc_prod_vert(Rnr,1:k_zi))+rest_zi*seg_conc_prod_vert(Rnr,k_zi+1)
                seg_conc_mult(Rnr) = (sum(seg_conc(1,1:k_zi,Rnr)) + rest_zi*seg_conc(1,k_zi+1,Rnr))

              case(2) ! A + B -> ...
                do it_k=1,kmax
                  do it_j=2,j1
                    do it_i=2,i1
                      seg_concl(1,it_k)    = seg_concl(1,it_k) + y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp1)
                      seg_concl(2,it_k)    = seg_concl(2,it_k) + y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp2)
                      seg_conc_prodl(it_k) = seg_conc_prodl(it_k) + y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp1)*y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp2)
                    enddo
                  enddo
                enddo
                do it_l=1,2
                  call MPI_ALLREDUCE(seg_concl(it_l,:) ,seg_conc(it_l,:,Rnr) ,kmax ,MY_REAL ,MPI_SUM ,comm3d ,mpierr)
                enddo
                call MPI_ALLREDUCE(seg_conc_prodl(:) ,seg_conc_prod_vert(Rnr,:) ,kmax ,MY_REAL ,MPI_SUM ,comm3d ,mpierr)
                seg_conc_mult_vert(Rnr,:) = seg_conc(1,:,Rnr)*seg_conc(2,:,Rnr)/rslabs
                seg_conc_prod(Rnr) = sum(seg_conc_prod_vert(Rnr,1:k_zi))+rest_zi*seg_conc_prod_vert(Rnr,k_zi+1)
                seg_conc_mult(Rnr) = (sum(seg_conc(1,1:k_zi,Rnr)) + rest_zi*seg_conc(1,k_zi+1,Rnr)) &
                                   * (sum(seg_conc(2,1:k_zi,Rnr)) + rest_zi*seg_conc(2,k_zi+1,Rnr))/(rslabs*(k_zi+rest_zi))
              case(3) ! A^a -> ...
                do it_k=1,kmax
                  do it_j=2,j1
                    do it_i=2,i1
                      seg_concl(1,it_k)    = seg_concl(1,it_k) + y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp1)
                      seg_conc_prodl(it_k) = seg_conc_prodl(it_k) + (y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp1) ** PL_scheme(n)%PL(j)%exp1)
                    enddo
                  enddo
                enddo
                do it_l=1,1
                  call MPI_ALLREDUCE(seg_concl(it_l,:) ,seg_conc(it_l,:,Rnr) ,kmax ,MY_REAL ,MPI_SUM ,comm3d ,mpierr)
                enddo
                call MPI_ALLREDUCE(seg_conc_prodl(:) ,seg_conc_prod_vert(Rnr,:) ,kmax ,MY_REAL ,MPI_SUM ,comm3d ,mpierr)
                seg_conc_mult_vert(Rnr,:) = (seg_conc(1,:,Rnr) ** PL_scheme(n)%PL(j)%exp1) / (rslabs ** (PL_scheme(n)%PL(j)%exp1 - 1))
                seg_conc_prod(Rnr) = sum(seg_conc_prod_vert(Rnr,1:k_zi))+rest_zi*seg_conc_prod_vert(Rnr,k_zi+1)
                seg_conc_mult(Rnr) = ((sum(seg_conc(1,1:k_zi,Rnr)) + rest_zi*seg_conc(1,k_zi+1,Rnr)) ** PL_scheme(n)%PL(j)%exp1)&
                                   / ((rslabs*(k_zi+rest_zi)) ** (PL_scheme(n)%PL(j)%exp1 - 1))
                
              case(4) !A^a + B^b -> ...
                do it_k=1,kmax
                  do it_j=2,j1
                    do it_i=2,i1
                      seg_concl(1,it_k)    = seg_concl(1,it_k) + y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp1)
                      seg_concl(2,it_k)    = seg_concl(2,it_k) + y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp2)
                      seg_conc_prodl(it_k) = seg_conc_prodl(it_k) + (y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp1) ** PL_scheme(n)%PL(j)%exp1) &
                                           * (y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp2) ** PL_scheme(n)%PL(j)%exp2)
                    enddo
                  enddo
                enddo
                do it_l=1,2
                  call MPI_ALLREDUCE(seg_concl(it_l,:) ,seg_conc(it_l,:,Rnr) ,kmax ,MY_REAL ,MPI_SUM ,comm3d ,mpierr)
                enddo
                call MPI_ALLREDUCE(seg_conc_prodl(:) ,seg_conc_prod_vert(Rnr,:) ,kmax ,MY_REAL ,MPI_SUM ,comm3d ,mpierr)
                seg_conc_mult_vert(Rnr,:) = (seg_conc(1,:,Rnr) ** PL_scheme(n)%PL(j)%exp1) * (seg_conc(2,:,Rnr) ** PL_scheme(n)%PL(j)%exp2) &
                                          / (rslabs ** (PL_scheme(n)%PL(j)%exp1 + PL_scheme(n)%PL(j)%exp2 - 1))
                seg_conc_prod(Rnr) = sum(seg_conc_prod_vert(Rnr,1:k_zi))+rest_zi*seg_conc_prod_vert(Rnr,k_zi+1)
                seg_conc_mult(Rnr) = ((sum(seg_conc(1,1:k_zi,Rnr)) + rest_zi*seg_conc(1,k_zi+1,Rnr)) ** PL_scheme(n)%PL(j)%exp1)&
                                   * ((sum(seg_conc(2,1:k_zi,Rnr)) + rest_zi*seg_conc(2,k_zi+1,Rnr)) ** PL_scheme(n)%PL(j)%exp2)&
                                   / ((rslabs*(k_zi+rest_zi)) ** (PL_scheme(n)%PL(j)%exp1 + PL_scheme(n)%PL(j)%exp2 - 1))
                
              case(5) ! A + B + C -> ...   
                do it_k=1,kmax
                  do it_j=2,j1
                    do it_i=2,i1
                      seg_concl(1,it_k)    = seg_concl(1,it_k) + y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp1)
                      seg_concl(2,it_k)    = seg_concl(2,it_k) + y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp2)
                      seg_concl(3,it_k)    = seg_concl(3,it_k) + y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp3)
                      seg_conc_prodl(it_k) = seg_conc_prodl(it_k) + y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp1) &
                                           * y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp2) * y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp3)
                    enddo
                  enddo
                enddo
                do it_l=1,3
                  call MPI_ALLREDUCE(seg_concl(it_l,:) ,seg_conc(it_l,:,Rnr) ,kmax ,MY_REAL ,MPI_SUM ,comm3d ,mpierr)
                enddo
                call MPI_ALLREDUCE(seg_conc_prodl(:) ,seg_conc_prod_vert(Rnr,:) ,kmax ,MY_REAL ,MPI_SUM ,comm3d ,mpierr)
                seg_conc_mult_vert(Rnr,:) = seg_conc(1,:,Rnr)*seg_conc(2,:,Rnr)*seg_conc(3,:,Rnr)/(rslabs**2)
                seg_conc_prod(Rnr) = sum(seg_conc_prod_vert(Rnr,1:k_zi))+rest_zi*seg_conc_prod_vert(Rnr,k_zi+1)
                seg_conc_mult(Rnr) = (sum(seg_conc(1,1:k_zi,Rnr)) + rest_zi*seg_conc(1,k_zi+1,Rnr)) &
                                   * (sum(seg_conc(2,1:k_zi,Rnr)) + rest_zi*seg_conc(2,k_zi+1,Rnr)) &
                                   * (sum(seg_conc(3,1:k_zi,Rnr)) + rest_zi*seg_conc(3,k_zi+1,Rnr)) / ((rslabs*(k_zi+rest_zi))**2)
                
              case(6) ! A^a + B^b + C^c -> ...
                do it_k=1,kmax
                  do it_j=2,j1
                    do it_i=2,i1
                      seg_concl(1,it_k)    = seg_concl(1,it_k) + y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp1)
                      seg_concl(2,it_k)    = seg_concl(2,it_k) + y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp2)
                      seg_concl(3,it_k)    = seg_concl(3,it_k) + y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp3)
                      seg_conc_prodl(it_k) = seg_conc_prodl(it_k) + (y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp1) ** PL_scheme(n)%PL(j)%exp1) &
                                           * (y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp2) ** PL_scheme(n)%PL(j)%exp2) &
                                           * (y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp3) ** PL_scheme(n)%PL(j)%exp3)
                    enddo
                  enddo
                enddo
                do it_l=1,3
                  call MPI_ALLREDUCE(seg_concl(it_l,:) ,seg_conc(it_l,:,Rnr) ,kmax ,MY_REAL ,MPI_SUM ,comm3d ,mpierr)
                enddo
                call MPI_ALLREDUCE(seg_conc_prodl(:) ,seg_conc_prod_vert(Rnr,:) ,kmax ,MY_REAL ,MPI_SUM ,comm3d ,mpierr)
                seg_conc_mult_vert(Rnr,:) = (seg_conc(1,:,Rnr) ** PL_scheme(n)%PL(j)%exp1) &
                                          * (seg_conc(2,:,Rnr) ** PL_scheme(n)%PL(j)%exp2) &
                                          * (seg_conc(3,:,Rnr) ** PL_scheme(n)%PL(j)%exp3) &
                                          / (rslabs ** (PL_scheme(n)%PL(j)%exp1 + PL_scheme(n)%PL(j)%exp2 + PL_scheme(n)%PL(j)%exp3 - 1))
                seg_conc_prod(Rnr) = sum(seg_conc_prod_vert(Rnr,1:k_zi))+rest_zi*seg_conc_prod_vert(Rnr,k_zi+1)
                seg_conc_mult(Rnr) = ((sum(seg_conc(1,1:k_zi,Rnr)) + rest_zi*seg_conc(1,k_zi+1,Rnr)) ** PL_scheme(n)%PL(j)%exp1)&
                                   * ((sum(seg_conc(2,1:k_zi,Rnr)) + rest_zi*seg_conc(2,k_zi+1,Rnr)) ** PL_scheme(n)%PL(j)%exp2)&
                                   * ((sum(seg_conc(3,1:k_zi,Rnr)) + rest_zi*seg_conc(3,k_zi+1,Rnr)) ** PL_scheme(n)%PL(j)%exp3)&
                                   / ((rslabs*(k_zi+rest_zi)) ** (PL_scheme(n)%PL(j)%exp1 + PL_scheme(n)%PL(j)%exp2 + PL_scheme(n)%PL(j)%exp3 - 1))
                
              case(7) ! A^a + B^b + C^c + D^d -> ... 
                do it_k=1,kmax
                  do it_j=2,j1
                    do it_i=2,i1
                      seg_concl(1,it_k)    = seg_concl(1,it_k) + y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp1)
                      seg_concl(2,it_k)    = seg_concl(2,it_k) + y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp2)
                      seg_concl(3,it_k)    = seg_concl(3,it_k) + y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp3)
                      seg_concl(4,it_k)    = seg_concl(4,it_k) + y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp4)
                      seg_conc_prodl(it_k) = seg_conc_prodl(it_k) + (y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp1) ** PL_scheme(n)%PL(j)%exp1) &
                                           * (y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp2) ** PL_scheme(n)%PL(j)%exp2) &
                                           * (y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp3) ** PL_scheme(n)%PL(j)%exp3) &
                                           * (y(it_i,it_j,it_k,PL_scheme(n)%PL(j)%comp4) ** PL_scheme(n)%PL(j)%exp4)
                    enddo
                  enddo
                enddo
                do it_l=1,4
                  call MPI_ALLREDUCE(seg_concl(it_l,:) ,seg_conc(it_l,:,Rnr) ,kmax ,MY_REAL ,MPI_SUM ,comm3d ,mpierr)
                enddo
                call MPI_ALLREDUCE(seg_conc_prodl(:) ,seg_conc_prod_vert(Rnr,:) ,kmax ,MY_REAL ,MPI_SUM ,comm3d ,mpierr)
                seg_conc_mult_vert(Rnr,:) = (seg_conc(1,:,Rnr) ** PL_scheme(n)%PL(j)%exp1) &
                                          * (seg_conc(2,:,Rnr) ** PL_scheme(n)%PL(j)%exp2) &
                                          * (seg_conc(3,:,Rnr) ** PL_scheme(n)%PL(j)%exp3) &
                                          * (seg_conc(4,:,Rnr) ** PL_scheme(n)%PL(j)%exp4) &
                                          / (rslabs ** (PL_scheme(n)%PL(j)%exp1 + PL_scheme(n)%PL(j)%exp2 &
                                                      + PL_scheme(n)%PL(j)%exp3 + PL_scheme(n)%PL(j)%exp4 - 1))
                seg_conc_prod(Rnr) = sum(seg_conc_prod_vert(Rnr,1:k_zi))+rest_zi*seg_conc_prod_vert(Rnr,k_zi+1)
                seg_conc_mult(Rnr) = ((sum(seg_conc(1,1:k_zi,Rnr)) + rest_zi*seg_conc(1,k_zi+1,Rnr)) ** PL_scheme(n)%PL(j)%exp1)&
                                   * ((sum(seg_conc(2,1:k_zi,Rnr)) + rest_zi*seg_conc(2,k_zi+1,Rnr)) ** PL_scheme(n)%PL(j)%exp2)&
                                   * ((sum(seg_conc(3,1:k_zi,Rnr)) + rest_zi*seg_conc(3,k_zi+1,Rnr)) ** PL_scheme(n)%PL(j)%exp3)&
                                   * ((sum(seg_conc(4,1:k_zi,Rnr)) + rest_zi*seg_conc(4,k_zi+1,Rnr)) ** PL_scheme(n)%PL(j)%exp4)&
                                   / ((rslabs*(k_zi+rest_zi)) ** (PL_scheme(n)%PL(j)%exp1 + PL_scheme(n)%PL(j)%exp2 &
                                                                + PL_scheme(n)%PL(j)%exp3 + PL_scheme(n)%PL(j)%exp4 - 1))
              case default
                if(myid==0) then
                  Print *,'Reaction type of reaction #',Rnr,' is not supported for segregation calculations'
                endif
              end select !select case as function of the equation type
!
!              if(myid==0) then
!                Print *,'Finished segregation calculation for reaction #',Rnr
!              endif
!                
              reaction_ev(PL_scheme(n)%PL(j)%r_nr) = .TRUE. ! Reaction #<r_nr> has been evaluated
            endif !reaction evaluated yet?
          endif   !only PRODUCTION terms for chemical n
          enddo   !j=1,reactions for chemical n
        endif     !active
      enddo       !n=1,nchsp
      do j=1,tnor
        if(reaction_ev(j) .EQV. .FALSE.) then
          if(myid==0) then
            Print *,'Reaction ',j,' has no production terms and is therefore not evaluated'
          endif
        endif
      enddo

      where (seg_conc_mult .ne. 0)
        segregation = (seg_conc_prod-seg_conc_mult)/seg_conc_mult
      elsewhere
        segregation = -9999.0 !Error code, devision by 0
      endwhere
      where (seg_conc_mult_vert .ne. 0)
        segregation_vert = (seg_conc_prod_vert-seg_conc_mult_vert)/seg_conc_mult_vert
      elsewhere
        segregation_vert = -9999.0 !Error code, devision by 0
      endwhere

      if(myid==0) then
        open (ifoutput,file='seg.'//cexpnr,position='append') 
          write(formatstring,'(a,i3,a)') '(F9.2,1x,F7.2,',tnor,'(2x,e12.6))'
          write(ifoutput,formatstring) timee,zi,(segregation(i),i=1,tnor)
        close(ifoutput)

        nsecs   = nint(timee)
        nhrs    = int(nsecs/3600)
        nminut  = int(nsecs/60)-nhrs*60
        nsecs   = mod(nsecs,60)
        open (ifoutput,file='seg_h.'//cexpnr,position='append') 
          write(ifoutput,'(A,/A,F5.0,A,I4,A,I2.2,A,I2.2,A,/A)') &
          '#-----------------------------------------------------------------------------'&
          ,'#',(timeav_glob),'--- AVERAGING TIMESTEP --- '      &
          ,nhrs,':',nminut,':',nsecs      &
          ,'   HRS:MIN:SEC AFTER INITIALIZATION ','#-----------------------------------------------------------------------------'

          write(formatstring,'(a,i3,a)') '(a9,1x,a7,1x,a7,',tnor,'(2x,a5,i3.3,a4))'
          write(ifoutput,formatstring) '#Time [s]','z_i [m]','h [m]',('IsegR',i,' [-]',i=1,tnor)
 
          write(formatstring,'(a,i3,a)') '(F9.2,1x,F7.2,1x,F7.2,',tnor,'(2x,e12.6))'
          do it_k=1,kmax
            write(ifoutput,formatstring) timee,zi,(it_k-0.5)*dz,(segregation_vert(i,it_k),i=1,tnor)
          enddo !Loop over height
          write(ifoutput,'(//)') 
        close(ifoutput)
      endif !(myid == 0)
    endif   !timee>=tnextwrite
  endif     !lsegr .EQV. .true.


  if( timee>=tnextwrite) then
    allocate (writearrayg(k1,tnor+2))
    writearrayg = 0.
    do n=1,tnor+2
      call MPI_ALLREDUCE(writearray(:,n), writearrayg(:,n) ,k1, MY_REAL, MPI_SUM, comm3d, mpierr)
    enddo
    if(myid==0) then
      open(ifoutput,file='keffs.'//cexpnr,position='append')
      write(ifoutput,'(a,f8.1)') 'time =',timee
      write(formatstring,'(a,i3,a)') '(5x,',tnor+2,'(2x,a6,3x))'
      write(ifoutput,formatstring) (RC(i)%rname,i=1,tnor),'Tabs','convppb'
      write(formatstring,'(a,i3,a)') '(5x,',tnor+2,'(3x,i3,5x))'
      write(ifoutput,formatstring)( RC(i)%func1,i=1,tnor)
!      do pl=1,kmax
!        write(ifoutput,'(i3,x,19e11.4)') pl,(writearray(pl,i),i=1,tnor+2)
!      enddo
      do pl=1,kmax
        write(ifoutput,'(i3,1x,19e11.4)') pl,(writearrayg(pl,i)/nprocs,i=1,tnor+2)
      enddo

      close(ifoutput)
    endif !(myid == 0)
    deallocate(writearrayg)
  endif

  if (lchmovie .and. (mod(timee,dtchmovie)==0 )) then
    call chemmovie(ybegin)
    deallocate(k3d,ybegin)
  endif

  !close (100)
return

end subroutine twostep2


!c***********************************************************************

subroutine NEWDT(t,te,dt,dtold,ratio,errlte,accept,dtmin,dtmax)
implicit none

  real    t,te,dt,dtold,ratio,errlte,ts,dtmin,dtmax,eps
  logical accept
  parameter (eps=1.e-12)

  if ((errlte.gt.1.0).and.(dt.gt.dtmin)) then
    accept=.false.
    ts=t
  else
    accept=.true.
    dtold=dt
    ts=t+dtold
  endif

  dt=max(0.5,min(2.0,0.8/sqrt(errlte+eps)))*dt
  dt=max(dtmin,min(dt,dtmax))
  call FIT(ts,te,dt)
  ratio=dt/dtold

return
end subroutine NEWDT

!c***********************************************************************

subroutine FIT(t,te,dt)
implicit none

 real  t,te,dt,rns
 integer  ns

 rns=(te-t)/dt
 if (rns.gt.10.0)then
   ! do nothing

 else
   ns=int(rns)+1
   dt=(te-t)/ns
   dt=(dt+t)-t
 endif

return
end subroutine FIT

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine calc_K(k)
use modglobal, only :rlv, cp, i1,j1, imax,jmax,timeav_glob,timee
use modfields, only :thl0,exnf,qt0,ql0,presf,svm
use modmpi, only : myid
implicit none

  integer k
  integer i

  if (lchconst .EQV. .false.) then
    T_abs(:,:) = thl0(2:i1,2:j1,k) * exnf(k) + (rlv/cp) * ql0(2:i1,2:j1,k)
    ! concentrations are in ppb reactionconstant are1 cm3/moleculer.s =>convppb
    ! if K in
  !  convppb(:,:) = 6.023e8 * (presf(k)/100) / (8.314e-2 * T_abs(:,:))
    convppb(:,:) = Avogrado * 1.e-9 * 1e-6/8.314e-2/100 * (presf(k)) / ( T_abs(:,:))
     ! 2.46e10*1.e9 = 2.46e19     <= 0.7224*P/T
  else
    T_abs(:,:) = t_ref
    convppb(:,:) = Avogrado * 1.e-9 * 1e-6 * (p_ref/100) / (8.314e-2 * t_ref)
  endif

  if (lchmovie .and. (mod(timee,dtchmovie)==0)) then
    k3d(:,:,k,tnor+1) = T_abs
    k3d(:,:,k,tnor+2) = convppb
  endif

  if( timee>=tnextwrite) then
     writearray(k,tnor+1)=sum(T_abs)/(imax*jmax)
     writearray(k,tnor+2)=sum(convppb)/(imax*jmax)
  endif

  if( lchmovie .and. (mod(timee,dtchmovie)==0)) then
     writearray(k,tnor+1)=sum(T_abs)/(imax*jmax)
     writearray(k,tnor+2)=sum(convppb)/(imax*jmax)
  endif

  do i=1, tnor
    if (RC(i)%RadDep == 1 ) then
      if (timee>=tnextwrite) then
        writearray(k,i) = sum(keff(:,:,RC(i)%Kindex,k))/(imax*jmax)
      endif
      if( lchmovie .and. (mod(timee,dtchmovie)==0)) then
        k3d(:,:,k,i) = keff(:,:,RC(i)%Kindex,k)
      endif
      !do nothing this is done in ratech
    else   ! not hv reaction func1 can have different meaning
      select case ( RC(i)%func1 )
      case(0)
        !do nothing K is in PPB and constant
      case(1) ! K in cm3/molecules*sec and independent of temperature
        keffT(:,:,RC(i)%Kindex) = RC(i)%A * convppb(:,:)
      case(2) !temperature depence of K
        keffT(:,:,RC(i)%Kindex) = RC(i)%A * exp(RC(i)%B / T_abs(:,:)) * convppb(:,:)
      case (3) !more complex temperature dependence
        keffT(:,:,RC(i)%Kindex) = RC(i)%A * (T_abs(:,:)/RC(i)%B)**RC(i)%C * exp(RC(i)%D / T_abs(:,:))* convppb(:,:)
      case(4:5) !use Lindemann-Hinshelwood equation  4 = second order 5 = first order so no conv_XXX factor
          !first CBL
        rk1(:,:)= RC(i)%A * (T_abs(:,:)/300)**RC(i)%B * exp(RC(i)%C/(T_abs(:,:))) * convppb(:,:) * 1e9
        rk2(:,:)= RC(i)%D * (T_abs(:,:)/300)**RC(i)%E * exp(RC(i)%F/(T_abs(:,:)))
        rk(:,:) = rk1 * rk2 /( rk1+rk2) * RC(i)%G
        if (RC(i)%func1 == 4) then
          keffT(:,:,RC(i)%Kindex) = rk(:,:) * convppb(:,:)
        else
          keffT(:,:,RC(i)%Kindex) = rk(:,:)
        endif
      case(6)!special function of reaction which is dependend on [H2O] but it is not a reactant
        !first CBL
        rk1(:,:) = RC(i)%A * exp(RC(i)%B / T_abs(:,:)) * convppb(:,:)
        rk2(:,:) = RC(i)%C * exp(RC(i)%D / T_abs(:,:)) * convppb(:,:)**2 * 1e9
        rk(:,:) =  RC(i)%E * exp(RC(i)%F / T_abs(:,:)) * svm(2:i1,2:j1,k,(H2O%loc+choffset)) * convppb(:,:)
        keffT(:,:,RC(i)%Kindex) = (rk1(:,:) + rk2(:,:)) * (1. + rk(:,:))
      case(7) ! same as 3 but third order so conv_ppb to the power 2
        keffT(:,:,RC(i)%Kindex) = RC(i)%A * (T_abs(:,:)/RC(i)%B)**RC(i)%C * exp(RC(i)%D / T_abs(:,:))* (convppb(:,:)**2)
      case default !if someone put by mistake a different number
        write (*,*) 'Unknown function specified'
        STOP
      end select

      if( timee>=tnextwrite) then
        writearray(k,i)=sum(keffT(:,:,RC(i)%Kindex))/(imax*jmax)
!        if (lsegrk .EQV. .true.) then
!          keffT3D(:,:,:,k) = keffT(:,:,:)
!        endif
      endif

      if (lchmovie .and. (mod(timee,dtchmovie)==0)) then
        k3d(:,:,k,i) = keffT(:,:,RC(i)%Kindex)
      endif
    endif

  end do !tnor

end subroutine calc_K

subroutine ratech
!
!-----------------------------------------------------------------|
!                                                                 |
!*** *ratech*  calculate the photolyis rate perturbed by clouds   |
!                                                                 |
!     Jordi Vila   WUR          16/08/2004                        |
!                                                                 |
!    purpose.                                                     |
!     --------
!
!     It calculates the photolysis rate perturbed by the clouds
!     It takes the clear sky value prescribed at input_chem
!     and it modifyes according to the parameterization developed
!     by Chang et al.(eqs 13) (JGR, vol 92, 14,681).
!     The parametrization depends on:
!     - Cloud optical depth
!     - Solar zenith angles
!
!**   interface.
!     ----------
! KvdD
!     if zbase then from top down to ztop to catch double(broken) clouds
!     inside cloud scaled with ql0
!
!                 hv
!     O3         ----> O(1D) + O2   (jO3)
!     O(1D) + H2O ---> 2HO
!     O(1D) + AIR ---> DEACTIVATION
!overall:
!     O3 + H2O ---> 2HO + O2         (JO3)
!-----------------------------------------------------------------
!

  use modglobal, only : i1,i2,ih, j1,j2,jh, k1,kmax, timeav_glob,pi,xtime,timee,xday,xlat,xlon, &
                        zf,dzf, iexpnr,rslabs,ifoutput,cexpnr
  use modfields, only : sv0, qt0, ql0 ,rhof
  use modmpi,    only : myid, comm3d, mpierr, mpi_max, my_real, mpi_integer, mpi_sum
  implicit none

  real  sza
  real     timeav
  integer   i,j,k,l,m,n,r
  real  zbase, ztop, qlint, clddepth
  real  zbasesum, ztopmax, zbasecount
  real  zbasesuml, ztopmaxl, zbasecountl
  real  cloudheight, cloudheightmax
  real  cloudheightmaxl
  real  cloudheightsum, cloudheightsuml
  real  cloudcount,cloudcountl
  real  qlintsum,qlintsuml
  real  qlintmax,qlintmaxl
  real  qlintallsum,qlintallsuml
  real  qlintallmax,qlintallmaxl
  real  coszen, coszenmax
  real  rhow, re, tr, fba, fab, tauc, tau2
  real  ks,diff,qsum
  real  xhr
  logical lday

  lday = .false.

  timeav = timeav_glob
  xhr = xtime + timee/3600.

  re    = 1.e-05 ! mean cloud drop radius
  rhow  = 1000.  ! density of water
  tauc  = 5.     ! Optically thin clouds (critical value)

!   --------------------------
!   Calculating reaction rates
!   --------------------------

  sza = getth(xday,xlat,xlon,xhr)
  coszen = max(0.0,cos(sza)) !to avoid negative values

  if (coszen > 0.00) then !it is day
    lday = .true.
    if(myid==0 .and. (switch .eqv. .false.)) then
      switch = .true.
      write(*,*)'The SUN is UP at timee=',timee,xhr
    endif
    if (ldiuvar .eqv. .false.) then ! we have to fix the sza to the fixed hour h_ref
      sza = getth(xday,xlat,xlon,h_ref)
      coszen = max(0.0,cos(sza))
!write(*,*)sza,coszen
    endif
  else
    lday = .false.
    if(myid==0 .and. (switch .eqv. .true.) ) then
      switch = .false.
      write(*,*)'The SUN is DOWN at timee=',timee,xhr
    endif
  endif

  if (lday .eqv. .false.) then
    ! It is night, no photolysis reactions
    keff(:,:,:,:) = 0.0
    RETURN
  else
    do i=1, nr_raddep
      r = raddep_RCindex(i) !array to photolysis reactionnumbers
      select case (RC(r)%func1)
      case (1) ! constant independent of sza
        RC(r)%Keff = RC(r)%A
      case (2) ! exponential function
        RC(r)%Keff = RC(r)%A * exp(RC(r)%B / coszen)
!write(*,*)'func 2',RC(r)%Keff
      case (3) ! powerfunction
        RC(r)%Keff = RC(r)%A * coszen ** RC(r)%B
      case (4) ! powerfunction but special for JO3
        RC(r)%Keff = RC(r)%A * coszen ** RC(r)%B
!write(*,*)'func 4',RC(r)%Keff
      case default
        RC(r)%Keff = 1.
      end select
    enddo
  endif

  !fill the keff with values for no clouds
  do i=1, nr_raddep
    r = raddep_RCindex(i)
    select case (RC(r)%func1)
    case (4) !special for jO3  JO3 = jO3 * k1[H2O] / (k1[H2O] + k2[Air])
      if(lchconst .eqv. .true.) then
        keff(:,:,i,:) = RC(r)%Keff * RC(r)%D *  (q_ref * MW_air/MW_h2o ) / &
         (RC(r)%D * q_ref * MW_air/MW_h2o + RC(r)%E * (1.- q_ref* MW_air/MW_h2o))
!write(*,*)'func 4', keff(2,2,i,2),RC(r)%Keff,RC(r)%D *  (q_ref * MW_air/MW_h2o ),RC(r)%D * q_ref * MW_air/MW_h2o,RC(r)%E * (1.- q_ref* MW_air/MW_h2o)
      else
        keff(:,:,i,:) = RC(r)%Keff * RC(r)%D * ( qt0(2:i1,2:j1,1:kmax) * MW_air/MW_h2o ) / &
         (RC(r)%D * ( qt0(2:i1,2:j1,1:kmax) * MW_air/MW_h2o ) + RC(r)%E * (1.- qt0(2:i1,2:j1,1:kmax)* MW_air/MW_h2o))
      endif
    case default ! for now for all other photolysis reactions
      keff(:,:,i,:) = RC(r)%Keff
    end select
  enddo

  if( sum(ql0) == 0.0 .or. (lcloudKconst  .eqv. .true.)) then  ! maybe < 0.01
    !there is no liquid water in the domain so no clouds
    !or we like constamnt value's for K so we are finished here
  else
    !we have to look for the individual clouds

    zbasesum = 0.
    zbasecount = 0.
    clddepth = 0.
    ztopmax = 0.
    cloudheightmax = 0.0
    cloudheightsum = 0.0
    cloudcount = 0.0
    qlintsum = 0.0
    qlintallsum = 0.0
    qlintallmax = 0.0

    !for clouds the the max solar zenith angle is cutoff at 60 degrees
    coszenmax = min(60*pi/180,coszen)

    do j=2,j1
      do i=2,i1
        kefftemp = 1.0
        do k=2,kmax
          if (ql0(i,j,k) > 0.0) then !found bottom of cloud at level k
            do L=kmax,k,-1  !continue from top down to bottom
              if (ql0(i,j,L).gt.0.0) then !found top of cloud at level L
                qsum  = sum(ql0(i,j,k:L))
                qlint = sum(ql0(i,j,k:L) * rhof(k:L) * dzf(k:L)) ! ??? ###############rhof(one) or rhof(L)   ###########################
                exit !L loop
              endif
            enddo

            zbase = zf(k)
            ztop = zf(L)
            clddepth = ztop - zbase
            cloudcount = cloudcount + 1
            qlintallsum = qlintallsum + qlint
            qlintallmax = max(qlintallmax,qlint)

            !- Calculating transmission coefficient, cloud optical depth
            tau2 = (3./2.)*(qlint/(rhow*re))
            if (tau2 >= tauc ) then  ! 'dense' cloud
              ! smooting of cloud base and top
              zbase = zf(k) - (ql0(i,j,k)/(ql0(i,j,k) + ql0(i,j,k+1))) * dzf(k)  !!!!! of dzf(k+-?) with non equidistant grid
              ztop  = zf(L) + (ql0(i,j,L)/(ql0(i,j,L) + ql0(i,j,L-1))) * dzf(L)  !!!!! of dzf(l+-?)
              !for cloud statistics
              zbasesum = zbasesum + zbase
              zbasecount = zbasecount + 1
              ztopmax = max(ztopmax, ztop)
              cloudheight = ztop - zbase
              cloudheightmax = max(cloudheightmax, cloudheight)
              cloudheightsum = cloudheightsum + cloudheight
              qlintsum = qlintsum + qlint
              qlintmax = max(qlintmax,qlint)

              tr   = (5. - exp(-tau2))/(4.+ 3.*tau2*0.14)

              fba  = 1.6 * tr * coszenmax      !factor below clouds
              kefftemp(1:k-1,:)  = fba
              ! Only the factor above clouds is dependent on the chemical species
              do m = 1, nr_raddep
                r = raddep_RCindex(m)
                !RC(r)%C is the alfa factor of Chang
                fab  = 1. + RC(r)%C * (1-tr) * coszenmax   !factor above clouds
                ! above cloud
                kefftemp(l:kmax,m) = fab

                !inside cloud
                diff = fab - fba
                ks = fab
       !if(myid==0)write(*,*) raddep_RCindex(m),'fab=',fab,'fba=',fba
                ! scale inside cloud with ql0
                ! not sure this is correct with non equidistant grid
                do n=l-1,k,-1
                  ks = ks - ql0(i,j,n)/qsum * diff
                  kefftemp(n,m) = ks
                enddo
                keff(i,j,m,:) = keff(i,j,m,:) * kefftemp(:,m)
        !if(myid==0)write (*,*) myid,'xxx',timee,m,r,RC(r)%func1, RC(r)%Keff,keff(i,j,m,1)
              enddo !m = 1, nr_raddep

            endif !dense cloud
            EXIT ! k loop
          endif !ql0(i,j,k) >0
        enddo !k loop
      enddo !i loop
    enddo !j loop
  endif !ql0(:,:,:) = 0


  if(timee>=tnextwrite) then
     ztopmaxl = ztopmax
     zbasesuml = zbasesum
     zbasecountl = zbasecount
     cloudheightmaxl = cloudheightmax
     cloudheightsuml = cloudheightsum
     cloudcountl = cloudcount
     qlintsuml = qlintsum
     qlintmaxl = qlintmax
     qlintallsuml = qlintallsum
     qlintallmaxl = qlintallmax

     call MPI_ALLREDUCE(zbasecountl   , zbasecount   , 1,    MY_REAL,  MPI_SUM, comm3d, mpierr)
     call MPI_ALLREDUCE(cloudcountl   , cloudcount   , 1,    MY_REAL,  MPI_SUM, comm3d, mpierr)

     call MPI_ALLREDUCE(zbasesuml, zbasesum, 1,    MY_REAL,  MPI_SUM, comm3d, mpierr)
     call MPI_ALLREDUCE(ztopmaxl, ztopmax, 1,    MY_REAL,  MPI_MAX, comm3d, mpierr)

     call MPI_ALLREDUCE(cloudheightmaxl, cloudheightmax, 1,    MY_REAL,  MPI_MAX, comm3d, mpierr)
     call MPI_ALLREDUCE(cloudheightsuml, cloudheightsum, 1,    MY_REAL,  MPI_SUM, comm3d, mpierr)

     call MPI_ALLREDUCE(qlintmaxl, qlintmax, 1,    MY_REAL,  MPI_MAX, comm3d, mpierr)
     call MPI_ALLREDUCE(qlintsuml, qlintsum, 1,    MY_REAL,  MPI_SUM, comm3d, mpierr)

     call MPI_ALLREDUCE(qlintallmaxl, qlintallmax, 1,    MY_REAL,  MPI_MAX, comm3d, mpierr)
     call MPI_ALLREDUCE(qlintallsuml, qlintallsum, 1,    MY_REAL,  MPI_SUM, comm3d, mpierr)

     if ( myid ==0 ) then
       open (ifoutput,file='cloudstat.'//cexpnr,position='append')
       write(ifoutput,'(f7.0,f6.2,2f7.2,4f9.1,4f9.4)')timee, xhr, zbasecount/rslabs, cloudcount/rslabs, zbasesum / (zbasecount+1.0e-5), &
         ztopmax, cloudheightsum/(zbasecount+1.0e-5), cloudheightmax, qlintsum / (zbasecount+1.0e-5), &
         qlintmax, qlintallsum / (cloudcount + 1.0e-5), qlintallmax
       close (ifoutput)
     endif
  endif

  RETURN

  end  subroutine ratech
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine ITER(gdt,k)
use modglobal, only : i2, j2,nsv,timee
use modmpi, only : myid
implicit none

  real    gdt
  integer k
  real,pointer,dimension(:,:,:) :: YPL_pointer
  integer j,n

  YL(:,:,:) = 0.
  YP(:,:,:) = 0.

  do n=1,nchsp
  if (PL_scheme(n)%active .EQV. .TRUE.) then
    if (PL_scheme(n)%name == H2O%name )  cycle    !don't do calculation of H2O
!    if (PL_scheme(n)%name == PRODUC%name) cycle

    do j=1, PL_scheme(n)%nr_PL

      if(PL_scheme(n)%PL(j)%PorL == 1) then
        YPL_pointer => YP
      else
        YPL_pointer => YL
      endif

      if (RC(PL_scheme(n)%PL(j)%r_nr)%raddep == 1 ) then !we have to use effective Keff from Ratech
        select case (PL_scheme(n)%PL(j)%formula)
        case (0)
          YPL_pointer(:,:,n) = YPL_pointer(:,:,n) + PL_scheme(n)%PL(j)%coef * keff(:,:,RC(PL_scheme(n)%PL(j)%r_nr)%Kindex,k)
        case (1)
          YPL_pointer(:,:,n) = YPL_pointer(:,:,n) + PL_scheme(n)%PL(j)%coef * keff(:,:,RC(PL_scheme(n)%PL(j)%r_nr)%Kindex,k) * ynew(:,:,PL_scheme(n)%PL(j)%comp1)
        case (2)
          YPL_pointer(:,:,n) = YPL_pointer(:,:,n) + PL_scheme(n)%PL(j)%coef * keff(:,:,RC(PL_scheme(n)%PL(j)%r_nr)%Kindex,k) * ynew(:,:,PL_scheme(n)%PL(j)%comp1) * ynew(:,:,PL_scheme(n)%PL(j)%comp2)
        case (3)
          YPL_pointer(:,:,n) = YPL_pointer(:,:,n) + PL_scheme(n)%PL(j)%coef * keff(:,:,RC(PL_scheme(n)%PL(j)%r_nr)%Kindex,k) * ynew(:,:,PL_scheme(n)%PL(j)%comp1) ** PL_scheme(n)%PL(j)%exp1
        case (4)
          YPL_pointer(:,:,n) = YPL_pointer(:,:,n) + PL_scheme(n)%PL(j)%coef * keff(:,:,RC(PL_scheme(n)%PL(j)%r_nr)%Kindex,k) * ynew(:,:,PL_scheme(n)%PL(j)%comp1)  ** PL_scheme(n)%PL(j)%exp1 &
                   * ynew(:,:,PL_scheme(n)%PL(j)%comp2)  ** PL_scheme(n)%PL(j)%exp2
        case (5)
          YPL_pointer(:,:,n) = YPL_pointer(:,:,n) + PL_scheme(n)%PL(j)%coef * keff(:,:,RC(PL_scheme(n)%PL(j)%r_nr)%Kindex,k) &
                   * ynew(:,:,PL_scheme(n)%PL(j)%comp1) &
                   * ynew(:,:,PL_scheme(n)%PL(j)%comp2) &
                   * ynew(:,:,PL_scheme(n)%PL(j)%comp3)
        case (6)
          YPL_pointer(:,:,n) = YPL_pointer(:,:,n) + PL_scheme(n)%PL(j)%coef * keff(:,:,RC(PL_scheme(n)%PL(j)%r_nr)%Kindex,k) &
                   * ynew(:,:,PL_scheme(n)%PL(j)%comp1) ** PL_scheme(n)%PL(j)%exp1 &
                   * ynew(:,:,PL_scheme(n)%PL(j)%comp2) ** PL_scheme(n)%PL(j)%exp2 &
                   * ynew(:,:,PL_scheme(n)%PL(j)%comp3) ** PL_scheme(n)%PL(j)%exp3
        case (7)
          YPL_pointer(:,:,n) = YPL_pointer(:,:,n) + PL_scheme(n)%PL(j)%coef * keff(:,:,RC(PL_scheme(n)%PL(j)%r_nr)%Kindex,k) &
                   * ynew(:,:,PL_scheme(n)%PL(j)%comp1) ** PL_scheme(n)%PL(j)%exp1 &
                   * ynew(:,:,PL_scheme(n)%PL(j)%comp2) ** PL_scheme(n)%PL(j)%exp2 &
                   * ynew(:,:,PL_scheme(n)%PL(j)%comp3) ** PL_scheme(n)%PL(j)%exp3 &
                   * ynew(:,:,PL_scheme(n)%PL(j)%comp4) ** PL_scheme(n)%PL(j)%exp4
        end select
      else
        kreact(:,:) = PL_scheme(n)%PL(j)%coef * KeffT(:,:,RC(PL_scheme(n)%PL(j)%r_nr)%Kindex)
        select case (PL_scheme(n)%PL(j)%formula)
        case (0)
          YPL_pointer(:,:,n) = YPL_pointer(:,:,n) + kreact(:,:)
        case (1)
          YPL_pointer(:,:,n) = YPL_pointer(:,:,n) + kreact(:,:) * ynew(:,:,PL_scheme(n)%PL(j)%comp1)
        case (2)
          YPL_pointer(:,:,n) = YPL_pointer(:,:,n) + kreact(:,:) * ynew(:,:,PL_scheme(n)%PL(j)%comp1) * ynew(:,:,PL_scheme(n)%PL(j)%comp2)
        case (3)
          YPL_pointer(:,:,n) = YPL_pointer(:,:,n) + kreact(:,:) * ynew(:,:,PL_scheme(n)%PL(j)%comp1) ** PL_scheme(n)%PL(j)%exp1
        case (4)
          YPL_pointer(:,:,n) = YPL_pointer(:,:,n) + kreact(:,:) * ynew(:,:,PL_scheme(n)%PL(j)%comp1) ** PL_scheme(n)%PL(j)%exp1 &
                   * ynew(:,:,PL_scheme(n)%PL(j)%comp2)  ** PL_scheme(n)%PL(j)%exp2
        case (5)
          YPL_pointer(:,:,n) = YPL_pointer(:,:,n) + kreact(:,:) * ynew(:,:,PL_scheme(n)%PL(j)%comp1) &
                                                                * ynew(:,:,PL_scheme(n)%PL(j)%comp2) &
                                                                * ynew(:,:,PL_scheme(n)%PL(j)%comp3)
        case (6)
          YPL_pointer(:,:,n) = YPL_pointer(:,:,n) + kreact(:,:) * ynew(:,:,PL_scheme(n)%PL(j)%comp1) ** PL_scheme(n)%PL(j)%exp1 &
                                                                * ynew(:,:,PL_scheme(n)%PL(j)%comp2) ** PL_scheme(n)%PL(j)%exp2 &
                                                                * ynew(:,:,PL_scheme(n)%PL(j)%comp3) ** PL_scheme(n)%PL(j)%exp3
        case (7)
          YPL_pointer(:,:,n) = YPL_pointer(:,:,n) + kreact(:,:) * ynew(:,:,PL_scheme(n)%PL(j)%comp1) ** PL_scheme(n)%PL(j)%exp1 &
                                                                * ynew(:,:,PL_scheme(n)%PL(j)%comp2) ** PL_scheme(n)%PL(j)%exp2 &
                                                                * ynew(:,:,PL_scheme(n)%PL(j)%comp3) ** PL_scheme(n)%PL(j)%exp3 &
                                                                * ynew(:,:,PL_scheme(n)%PL(j)%comp4) ** PL_scheme(n)%PL(j)%exp4
        end select
      endif

    enddo
    ynew(:,:,n)=max(0.0,(ysum(:,:,n)+gdt*YP(:,:,n))/(1.0+gdt*YL(:,:,n)))
  else
    ynew(:,:,n) = ysum(:,:,n)
  endif !active

  enddo !n=1,nchsp

end subroutine ITER


!c
!c ---- Function to calculate solar zenith angle
!c
real function getth(daynr,lat,lon,xhr)
implicit none
  real daynr, lat, lon, xhr
  real  houra
  real  obliq,deday,delta,lonr,latr
  real  piby,pi

  pi = acos(-1.)
  piby = pi/ 180.
  lonr = lon*piby
  latr = lat*piby
  obliq = 23.45 * piby
  deday = 4.88 + 2*pi/365  * daynr
  delta = asin(sin(obliq)*sin(deday))
  houra = lonr - pi + xhr * (2.*pi/24.)
  getth = acos(sin(delta)*sin(latr) + cos(delta)*cos(latr)*cos(houra))

  return
end function getth

subroutine chemmovie(ybegin)
use modglobal, only : i1,j1,ih,jh,i2, j2, k1,kmax,nsv, timee, iexpnr
use modfields, only: svm,qt0,ql0
use modmpi, only: myid
implicit none

  real*4, allocatable ::dummy(:,:,:)
  real ,dimension(2-ih:i1+ih,2-jh:j1+jh,k1,nchsp)::ybegin
  integer fileout,i,j,k,n

  character (len=20) ::filenaam
  character (len=8)  ::id

  allocate(dummy(2:i1,2:j1,kmax))

  fileout=20
  write(id,'(a,i3.3,a,i3.3)')'.',myid,'.',iexpnr
  do n=1,nchsp
    write(filenaam,'(a,a)')trim(PL_scheme(n)%name),id
    open(fileout,file=filenaam,form='unformatted',position='append',action='write')
    dummy(:,:,:)=svm(2:i1,2:j1,1:kmax,n)
     write(fileout) (((dummy(i,j,k),i=2,i1),j=2,j1),k=1,kmax)
    close(fileout)
  enddo

  do n=1,nchsp
    write(filenaam,'(a,a,a)')'PL_',trim(PL_scheme(n)%name),id
    open(fileout,file=filenaam,form='unformatted',position='append',action='write')
    dummy(:,:,:) = svm(2:i1,2:j1,1:kmax,n) - ybegin(2:i1,2:j1,1:kmax,n)
     write(fileout) (((dummy(i,j,k),i=2,i1),j=2,j1),k=1,kmax)
    close(fileout)
  enddo

  do n=1,tnor
    write(filenaam,'(a,a)')trim(RC(n)%rname),id
    open(fileout,file=filenaam,form='unformatted',position='append',action='write')
    write(fileout)k3d(:,:,:,n)
    close(fileout)
  enddo

  write(filenaam,'(a,a)')'T_abs',id
  open(fileout,file=filenaam,form='unformatted',position='append',action='write')
  write(fileout)k3d(:,:,:,tnor+1)
  close(fileout)

  write(filenaam,'(a,a)')'conv',id
  open(fileout,file=filenaam,form='unformatted',position='append',action='write')
  write(fileout)k3d(:,:,:,tnor+2)
  close(fileout)

  write(filenaam,'(a,a)')'qt',id
  open(fileout,file=filenaam,form='unformatted',position='append',action='write')
  dummy(:,:,:) = qt0(2:i1,2:j1,1:kmax)
  write(fileout) dummy
  close(fileout)

  write(filenaam,'(a,a)')'ql',id
  open(fileout,file=filenaam,form='unformatted',position='append',action='write')
  dummy(:,:,:) = ql0(2:i1,2:j1,1:kmax)
  write(fileout) dummy
  close(fileout)

  deallocate(dummy)

end subroutine chemmovie

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine write_chem_scheme(chem_name)
use modglobal, only: nsv,ifoutput,cexpnr
use modmpi,    only: myid
implicit none

  character*6, dimension(nchsp) :: chem_name(:)
  integer i,j
  character (len=4) coef_str

  if ( myid == 0 ) then

  !Print out the reaction schemes
  open(ifoutput,file='reaction_scheme.'//CEXPNR,RECL=132)
  write(ifoutput,*)' '
  do i=1,nchsp
  write(ifoutput,*)'---------------------------------------'
  write(ifoutput,*)' '
 if (PL_scheme(i)%active .EQV. .FALSE. ) then
    write(ifoutput,'(a6,a2,i2,a,a)')PL_scheme(i)%name,'(',i,')','NO REACTION'
    write(ifoutput,*)' '
  else
    write(ifoutput,'(a6,a2,i2,a)')PL_scheme(i)%name,'(',i,')'
    write(ifoutput,'(a5)') 'YP = '
    do j=1,PL_scheme(i)%nr_PL
      write(coef_str,'(f4.2)')PL_scheme(i)%PL(j)%coef
      if (PL_scheme(i)%PL(j)%PorL == PRODUCTION ) then
        select case (PL_scheme(i)%PL(j)%formula)
        case (0)
          write(ifoutput,*) '    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),')  ( F=',PL_scheme(i)%PL(j)%formula,')'
        case (1)
          write(ifoutput,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),']','( F=',PL_scheme(i)%PL(j)%formula,')'
        case (2)
          write(ifoutput,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),']', &
                    '* Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp2)),'] ( F=',PL_scheme(i)%PL(j)%formula,')'
        case (3)
          write(ifoutput,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),'] ** (',PL_scheme(i)%PL(j)%exp1,')','( F=',PL_scheme(i)%PL(j)%formula,')'
        case(4)
          write(ifoutput,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),'] ** (',PL_scheme(i)%PL(j)%exp1,')', &
                    ' * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp2)),']** (',PL_scheme(i)%PL(j)%exp2,') ( F=',PL_scheme(i)%PL(j)%formula,')'
        case(5)
          write(ifoutput,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname), &
                    ') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),']', &
                    '* Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp2)), &
                    '] * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp3)),'] ( F=',PL_scheme(i)%PL(j)%formula,')'
        case(6)
          write(ifoutput,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname), &
                    ') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),']** (',PL_scheme(i)%PL(j)%exp1, &
                    ') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp2)),'] ** (',PL_scheme(i)%PL(j)%exp2, &
                    ') *  Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp3)),'] ** (',PL_scheme(i)%PL(j)%exp3,')( F=',PL_scheme(i)%PL(j)%formula,')'
        case(7)
          write(ifoutput,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname), &
                    ') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),']** (',PL_scheme(i)%PL(j)%exp1, &
                    ') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp2)),'] ** (',PL_scheme(i)%PL(j)%exp2, &
                    ') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp3)),'] ** (',PL_scheme(i)%PL(j)%exp3, &
                    ')* Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp4)),'] ** (',PL_scheme(i)%PL(j)%exp4,')( F=',PL_scheme(i)%PL(j)%formula,')'
        end select
      endif
    enddo

    write(ifoutput,'(a5)') 'YL = '
    do j=1,PL_scheme(i)%nr_PL
      if (PL_scheme(i)%PL(j)%PorL == LOSS ) then
        write(coef_str,'(f4.2)')PL_scheme(i)%PL(j)%coef
        select case (PL_scheme(i)%PL(j)%formula)
        case (0)
          write(ifoutput,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),')  ( F=',PL_scheme(i)%PL(j)%formula,')'
        case (1)
          write(ifoutput,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),']','( F=',PL_scheme(i)%PL(j)%formula,')'
        case (2)
          write(ifoutput,*) '    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),']', &
                    '* Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp2)),'] ( F=',PL_scheme(i)%PL(j)%formula,')'
        case (3)
          write(ifoutput,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),'] ** (',PL_scheme(i)%PL(j)%exp1,') ( F=',PL_scheme(i)%PL(j)%formula,')'
        case(4)
          write(ifoutput,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),'] ** (',PL_scheme(i)%PL(j)%exp1,')', &
                    ' * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp2)),']** (',PL_scheme(i)%PL(j)%exp2,') ( F=',PL_scheme(i)%PL(j)%formula,')'
        case(5)
          write(ifoutput,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname), &
                    ') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),']', &
                    '* Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp2)), &
                    '] * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp3)),'] ( F=',PL_scheme(i)%PL(j)%formula,')'
        case(6)
          write(ifoutput,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname), &
                    ') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),']** (',PL_scheme(i)%PL(j)%exp1, &
                    ') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp2)),'] ** (',PL_scheme(i)%PL(j)%exp2, &
                    ') *  Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp3)),'] ** (',PL_scheme(i)%PL(j)%exp3,')( F=',PL_scheme(i)%PL(j)%formula,')'
        case(7)
          write(ifoutput,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname), &
                    ') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),']** (',PL_scheme(i)%PL(j)%exp1, &
                    ') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp2)),'] ** (',PL_scheme(i)%PL(j)%exp2, &
                    ') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp3)),'] ** (',PL_scheme(i)%PL(j)%exp3, &
                    ')* Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp4)),'] ** (',PL_scheme(i)%PL(j)%exp4,')( F=',PL_scheme(i)%PL(j)%formula,')'
        end select
      endif
    enddo
    write(ifoutput,*) ' '
  endif
  enddo

!  do i=1,nchsp
!  write(ifoutput,*)'---------------------------------------'
!  write(ifoutput,*)' '
!  if (PL_scheme(i)%active .EQV. .FALSE. ) then
!    write(ifoutput,*)PL_scheme(i)%name,'(',i,')'
!      write(ifoutput,*)' '
!  else
!    write(ifoutput,*)PL_scheme(i)%name,'(',i,')'
!    do j=1,PL_scheme(i)%nr_PL
!      if (PL_scheme(i)%PL(j)%PorL  == 1 ) then
!        select case (PL_scheme(i)%PL(j)%formula)
!        case (0)
!          write(ifoutput,*) '   YPtemp = YPtemp +',PL_scheme(i)%PL(j)%coef,' * K(',PL_scheme(i)%PL(j)%r_nr,')'
!        case (1)
!          write(ifoutput,*) '   YPtemp = YPtemp +',PL_scheme(i)%PL(j)%coef,' * K(',RC(PL_scheme(i)%PL(j)%r_nr)%rname,') * Y[',chem_name(PL_scheme(i)%PL(j)%comp1),']'
!        case (2)
!          write(ifoutput,*) '   YPtemp = YPtemp +',PL_scheme(i)%PL(j)%coef,' * K(',RC(PL_scheme(i)%PL(j)%r_nr)%rname,') * Y[',chem_name(PL_scheme(i)%PL(j)%comp1),']', &
!                    '* Y[',chem_name(PL_scheme(i)%PL(j)%comp2),']'
!        case (3)
!          write(ifoutput,*) '   YPtemp = YPtemp +',PL_scheme(i)%PL(j)%coef,' * K(',RC(PL_scheme(i)%PL(j)%r_nr)%rname,') * Y[',chem_name(PL_scheme(i)%PL(j)%comp1),'] ** (',PL_scheme(i)%PL(j)%exp1,')'
!        case(4)
!          write(ifoutput,*) '   YPtemp = YPtemp +',PL_scheme(i)%PL(j)%coef,' * K(',RC(PL_scheme(i)%PL(j)%r_nr)%rname,') * Y[',chem_name(PL_scheme(i)%PL(j)%comp1),'] ** (',PL_scheme(i)%PL(j)%exp1,')', &
!                    ' * Y[',chem_name(PL_scheme(i)%PL(j)%comp2),']** (',PL_scheme(i)%PL(j)%exp2,')'
!        end select
!      endif
!    enddo
!    write(ifoutput,*)'YP(',PL_scheme(i)%name,') = YPtemp'
!
!    do j=1,PL_scheme(i)%nr_PL
!      if (PL_scheme(i)%PL(j)%PorL == 2 ) then
!        select case (PL_scheme(i)%PL(j)%formula)
!        case (0)
!          write(ifoutput,*) '   YLtemp = YLtemp +',PL_scheme(i)%PL(j)%coef,' * K(',RC(PL_scheme(i)%PL(j)%r_nr)%rname,')'
!        case (1)
!          write(ifoutput,*) '   YLtemp = YLtemp +',PL_scheme(i)%PL(j)%coef,' * K(',RC(PL_scheme(i)%PL(j)%r_nr)%rname,') * Y[',chem_name(PL_scheme(i)%PL(j)%comp1),']'
!        case (2)
!          write(ifoutput,*) '   YLtemp = YLtemp +',PL_scheme(i)%PL(j)%coef,' * K(',RC(PL_scheme(i)%PL(j)%r_nr)%rname,') * Y[',chem_name(PL_scheme(i)%PL(j)%comp1),']', &
!                    '* Y[',chem_name(PL_scheme(i)%PL(j)%comp2),']'
!        case (3)
!          write(ifoutput,*) '   YLtemp = YLtemp +',PL_scheme(i)%PL(j)%coef,' * K(',RC(PL_scheme(i)%PL(j)%r_nr)%rname,') * Y[',chem_name(PL_scheme(i)%PL(j)%comp1),'] ** (',PL_scheme(i)%PL(j)%exp1,')'
!        case(4)
!          write(ifoutput,*) '   YLtemp = YLtemp +',PL_scheme(i)%PL(j)%coef,' * K(',RC(PL_scheme(i)%PL(j)%r_nr)%rname,') * Y[',chem_name(PL_scheme(i)%PL(j)%comp1),'] ** (',PL_scheme(i)%PL(j)%exp1,')', &
!                    ' * Y[',chem_name(PL_scheme(i)%PL(j)%comp2),']** (',PL_scheme(i)%PL(j)%exp2,')'
!        end select
!      endif
!    enddo
!    write(ifoutput,*)'YL(',PL_scheme(i)%name,') = YLtemp'
!    write(ifoutput,*) ' '
!  endif
!  enddo
  close(ifoutput)
  endif !if myid==0

end subroutine write_chem_scheme

end module modchem
