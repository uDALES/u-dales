program parameter_extractor
    implicit none
    
    ! Parameter variables (subset for testing)
    integer :: iexpnr, nprocx, nprocy, itot, jtot, ktot
    real :: runtime, dtmax, courant, randu
    real :: xlen, ylen, ps
    real :: tfielddump, tsample, tstatsdump, tstatstart
    real :: u0, v0, thl0, qt0, facT
    logical :: lwarmstart, lstratstart, ladaptive, lrandomize
    logical :: libm, lles, ltempeq, lbuoyancy, lmoist, lcoriol
    logical :: lqlnr, lfielddump
    integer :: ipoiss, iadv_mom, iadv_tke, iadv_thl, iadv_qt, iadv_sv
    character(len=50) :: fieldvars
    
    ! Set some default values for demonstration
    ! In a real implementation, these would be read from the input files
    iexpnr = 1
    runtime = 3600.0
    dtmax = 10.0
    nprocx = 2
    nprocy = 2
    itot = 64
    jtot = 64
    ktot = 64
    xlen = 500.0
    ylen = 500.0
    ps = 101325.0
    courant = 1.2
    randu = 0.05
    lwarmstart = .false.
    lstratstart = .false.
    ladaptive = .true.
    lrandomize = .true.
    libm = .true.
    lles = .true.
    ltempeq = .true.
    lbuoyancy = .true.
    lmoist = .false.
    lcoriol = .false.
    lqlnr = .false.
    ipoiss = 0
    iadv_mom = 2
    iadv_tke = 2
    iadv_thl = 2
    iadv_qt = 2
    iadv_sv = 7
    lfielddump = .true.
    tfielddump = 600.0
    fieldvars = 'u0,v0,w0,th'
    tsample = 10.0
    tstatsdump = 300.0
    tstatstart = 600.0
    u0 = 5.0
    v0 = 0.0
    thl0 = 300.0
    qt0 = 0.01
    facT = 293.15
    
    ! Output parameters in a structured format
    write(*,*) '=== u-DALES Parameter Extraction ==='
    write(*,*) 'RUN.iexpnr =', iexpnr
    write(*,*) 'RUN.runtime =', runtime
    write(*,*) 'RUN.dtmax =', dtmax
    write(*,*) 'RUN.nprocx =', nprocx
    write(*,*) 'RUN.nprocy =', nprocy
    write(*,*) 'RUN.courant =', courant
    write(*,*) 'RUN.randu =', randu
    write(*,*) 'RUN.lwarmstart =', lwarmstart
    write(*,*) 'RUN.lstratstart =', lstratstart
    write(*,*) 'RUN.ladaptive =', ladaptive
    write(*,*) 'RUN.lrandomize =', lrandomize
    write(*,*) 'RUN.libm =', libm
    write(*,*) 'RUN.lles =', lles
    write(*,*) 'DOMAIN.itot =', itot
    write(*,*) 'DOMAIN.jtot =', jtot
    write(*,*) 'DOMAIN.ktot =', ktot
    write(*,*) 'DOMAIN.xlen =', xlen
    write(*,*) 'DOMAIN.ylen =', ylen
    write(*,*) 'PHYSICS.ps =', ps
    write(*,*) 'PHYSICS.ltempeq =', ltempeq
    write(*,*) 'PHYSICS.lbuoyancy =', lbuoyancy
    write(*,*) 'PHYSICS.lmoist =', lmoist
    write(*,*) 'PHYSICS.lcoriol =', lcoriol
    write(*,*) 'DYNAMICS.lqlnr =', lqlnr
    write(*,*) 'DYNAMICS.ipoiss =', ipoiss
    write(*,*) 'DYNAMICS.iadv_mom =', iadv_mom
    write(*,*) 'DYNAMICS.iadv_tke =', iadv_tke
    write(*,*) 'DYNAMICS.iadv_thl =', iadv_thl
    write(*,*) 'DYNAMICS.iadv_qt =', iadv_qt
    write(*,*) 'DYNAMICS.iadv_sv =', iadv_sv
    write(*,*) 'OUTPUT.lfielddump =', lfielddump
    write(*,*) 'OUTPUT.tfielddump =', tfielddump
    write(*,*) 'OUTPUT.tsample =', tsample
    write(*,*) 'OUTPUT.tstatsdump =', tstatsdump
    write(*,*) 'OUTPUT.tstatstart =', tstatstart
    write(*,*) 'INPS.u0 =', u0
    write(*,*) 'INPS.v0 =', v0
    write(*,*) 'INPS.thl0 =', thl0
    write(*,*) 'INPS.qt0 =', qt0
    write(*,*) 'INPS.facT =', facT
    write(*,*) '==================================='
    
end program parameter_extractor