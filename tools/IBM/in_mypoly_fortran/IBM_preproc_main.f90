program IBM_preproc
    use omp_lib
    use in_mypoly_functions
    use boundaryMasking
    use matchFacets2Cells
    use IBM_preproc_io
    implicit none

    integer                                :: n_threads
    real                                   :: tol
    
    real                                   :: dx, dy  !, dz
    real   , allocatable, dimension(:)     :: xf, xh, yf, yh, zf, zh

    real   , allocatable, dimension(:)     :: vertices, faceNormals, incenters
    integer, allocatable, dimension(:)     :: facets
    real   , allocatable, dimension(:,:)   :: verts, faceNormal
    integer, allocatable, dimension(:,:)   :: faces
    integer                                :: n_vert, n_fcts
    integer                                :: stl_ground_i

    real   ,              dimension(3)     :: Ray_dir_u, Ray_dir_v, Ray_dir_w, Ray_dir_c
    integer                                :: diag_neighbs_i, periodic_x_i, periodic_y_i
    
    logical, allocatable, dimension(:,:,:) :: solid_u, solid_v, solid_w, solid_c

    logical, allocatable, dimension(:,:,:) :: fluid_IB_u, fluid_IB_v, fluid_IB_w, fluid_IB_c, &
                                              solid_IB_u, solid_IB_v, solid_IB_w, solid_IB_c
    real   , allocatable, dimension(:,:)   :: fluid_IB_xyz_u, fluid_IB_xyz_v, fluid_IB_xyz_w, fluid_IB_xyz_c

    integer, allocatable, dimension(:)     :: secfacids_u, secbndptids_u, secfacids_v, secbndptids_v, &
                                              secfacids_w, secbndptids_w, secfacids_c, secbndptids_c
    real   , allocatable, dimension(:)     :: secareas_u, bnddst_u, secareas_v, bnddst_v, secareas_w, bnddst_w, secareas_c, bnddst_c

    integer                                :: nfluid_IB_u, nfluid_IB_v, nfluid_IB_w, nfluid_IB_c, &
                                              nfacsecs_u, nfacsecs_v, nfacsecs_w, nfacsecs_c

    logical :: stl_ground, diag_neighbs, periodic_x, periodic_y
    real    :: max_height, L_char
    integer :: i, itot, jtot, ktot
    real(8) :: start_time, end_time

    start_time = OMP_GET_WTIME()

    open(unit=50,file='inmypoly_inp_info.txt')
        read(unit=50,fmt='(f15.10,x,f15.10)') dx, dy  !, dz
        read(unit=50,fmt='(i5,x,i5,x,i5)') itot, jtot, ktot
        read(unit=50,fmt='(f15.10)') tol
        read(unit=50,fmt='(f15.10,x,f15.10,x,f15.10)') Ray_dir_u(1), Ray_dir_u(2), Ray_dir_u(3)
        read(unit=50,fmt='(f15.10,x,f15.10,x,f15.10)') Ray_dir_v(1), Ray_dir_v(2), Ray_dir_v(3)
        read(unit=50,fmt='(f15.10,x,f15.10,x,f15.10)') Ray_dir_w(1), Ray_dir_w(2), Ray_dir_w(3)
        read(unit=50,fmt='(f15.10,x,f15.10,x,f15.10)') Ray_dir_c(1), Ray_dir_c(2), Ray_dir_c(3)
        read(unit=50,fmt='(i8,x,i8)') n_vert, n_fcts
        read(unit=50,fmt='(i4)') n_threads
        read(unit=50,fmt='(i1,x,i1,x,i1,x,i1)') stl_ground_i, diag_neighbs_i, periodic_x_i, periodic_y_i
    close(unit=50)

    call num2logical(stl_ground  , stl_ground_i  , 'stl_ground'  )
    call num2logical(diag_neighbs, diag_neighbs_i, 'diag_neighbs')
    call num2logical(periodic_x  , periodic_x_i  , 'periodic_x'  )
    call num2logical(periodic_y  , periodic_y_i  , 'periodic_y'  )


    !!!!!!! Memory acclocation !!!!!!

    ! Mesh coordinates
    allocate(xf(itot))
    allocate(xh(itot))
    allocate(yf(jtot))
    allocate(yh(jtot))
    allocate(zf(ktot))
    allocate(zh(ktot))

    ! STL geometry information
    allocate(vertices(n_vert*3))
    allocate(facets(n_fcts*3))
    allocate(incenters(n_fcts*3))
    allocate(faceNormals(n_fcts*3))
    allocate(verts(n_vert, 3))
    allocate(faces(n_fcts, 3))
    allocate(faceNormal(n_fcts, 3))
    
    ! solid points flag
    allocate(solid_u(itot,jtot,ktot))
    allocate(solid_v(itot,jtot,ktot))
    allocate(solid_w(itot,jtot,ktot))
    allocate(solid_c(itot,jtot,ktot))

    ! fluid and solid boundary points flag
    allocate(fluid_IB_u(itot,jtot,ktot), solid_IB_u(itot,jtot,ktot))
    allocate(fluid_IB_v(itot,jtot,ktot), solid_IB_v(itot,jtot,ktot))
    allocate(fluid_IB_w(itot,jtot,ktot), solid_IB_w(itot,jtot,ktot))
    allocate(fluid_IB_c(itot,jtot,ktot), solid_IB_c(itot,jtot,ktot))


    !!!!!!! Cartesian fluid mesh generation and geometry information reading !!!!!!

    xf(1:itot) = [((i-1)*dx+(dx/2.0), i=1,itot)]
    xh(1:itot) = [((i-1)*dx, i=1,itot)]

    yf(1:jtot) = [((i-1)*dy+(dy/2.0), i=1,jtot)]
    yh(1:jtot) = [((i-1)*dy, i=1,jtot)]

    ! zf(1:ktot) = [((i-1)*dz+(dz/2.0), i=1,ktot)]
    ! zh(1:ktot) = [((i-1)*dz, i=1,ktot)]

    call read_data('vertices.txt',n_vert,'faces.txt',n_fcts,'zfgrid.txt','zhgrid.txt',ktot, &
                    vertices,facets,incenters,faceNormals,zf,zh)


    !!!!!!! Solid-Fluid identification !!!!!!
    
    max_height = MAXVAL(vertices(3:n_vert*3:3)) + tol
    L_char = max_facet_side(n_vert,vertices,n_fcts,facets) + tol

    ! u-grid
    solid_u = is_grid_in_mypoly_func(n_vert,vertices,n_fcts,facets,incenters,faceNormals, &
                                    itot,xh,jtot,yf,ktot,zf,Ray_dir_u,L_char,max_height,tol,n_threads)
    write(*,*) 'Determined solid points for u-grid.'

    ! v-grid
    solid_v = is_grid_in_mypoly_func(n_vert,vertices,n_fcts,facets,incenters,faceNormals, &
                                    itot,xf,jtot,yh,ktot,zf,Ray_dir_v,L_char,max_height,tol,n_threads)
    write(*,*) 'Determined solid points for v-grid.'

    ! w-grid
    solid_w = is_grid_in_mypoly_func(n_vert,vertices,n_fcts,facets,incenters,faceNormals, &
                                    itot,xf,jtot,yf,ktot,zh,Ray_dir_w,L_char,max_height,tol,n_threads)
    write(*,*) 'Determined solid points for w-grid.'

    ! c-grid
    solid_c = is_grid_in_mypoly_func(n_vert,vertices,n_fcts,facets,incenters,faceNormals, &
                                    itot,xf,jtot,yf,ktot,zf,Ray_dir_c,L_char,max_height,tol,n_threads)
    write(*,*) 'Determined solid points for c-grid.'

    !! Write solid_ files
    call print_solid_points_index(itot, jtot, ktot, solid_u, solid_v, solid_w, solid_c, &
                                  'solid_u.txt', 'solid_v.txt', 'solid_w.txt', 'solid_c.txt')
    write(*,*) 'Written solid_*.txt files.'


    !!!!!!! Boundary masking !!!!!!

    ! u-grid
    call boundaryMasks(fluid_IB_u, solid_IB_u, fluid_IB_xyz_u, nfluid_IB_u, &
                       'u', itot, jtot, ktot, xh, yf, zf, &
                       solid_u, diag_neighbs, stl_ground)

    ! v-grid
    call boundaryMasks(fluid_IB_v, solid_IB_v, fluid_IB_xyz_v, nfluid_IB_v, &
                       'v', itot, jtot, ktot, xf, yh, zf, &
                       solid_v, diag_neighbs, stl_ground)

    ! w-grid
    call boundaryMasks(fluid_IB_w, solid_IB_w, fluid_IB_xyz_w, nfluid_IB_w, &
                       'w', itot, jtot, ktot, xf, yf, zh, &
                       solid_w, diag_neighbs, stl_ground)

    ! c-grid
    call boundaryMasks(fluid_IB_c, solid_IB_c, fluid_IB_xyz_c, nfluid_IB_c, &
                       'c', itot, jtot, ktot, xf, yf, zf, &
                       solid_c, diag_neighbs, stl_ground)
    
    !! Write fluid_boundary_ files
    call print_IB_index(itot, jtot, ktot, fluid_IB_u, fluid_IB_v, fluid_IB_w, fluid_IB_c, &
                        'fluid_boundary_u.txt', 'fluid_boundary_v.txt', &
                        'fluid_boundary_w.txt', 'fluid_boundary_c.txt')
    write(*,*) 'Written fluid_boundary_*.txt files.'


    !!!!!!! MF2C: Matching Facets to Cells and Writing Facet Sections Files !!!!!!!

    verts      = reshape(vertices,    [n_vert, 3], order=[2,1])
    faces      = reshape(facets,      [n_fcts, 3], order=[2,1])
    faceNormal = reshape(faceNormals, [n_fcts, 3], order=[2,1])

    !! Computation of facet sections for each grid type
    ! u-grid
    call matchFacetsToCells(faces, faceNormal, n_fcts, verts, n_vert, &
         fluid_IB_u, solid_IB_u, fluid_IB_xyz_u, nfluid_IB_u, xh, yf, zf, itot, jtot, ktot, &
         diag_neighbs, periodic_x, periodic_y,  n_threads, &
         secfacids_u, secbndptids_u, secareas_u, bnddst_u, nfacsecs_u)
    write(*,*) 'Computation of facet_sections_u done.'

    ! v-grid
    call matchFacetsToCells(faces, faceNormal, n_fcts, verts, n_vert, &
         fluid_IB_v, solid_IB_v, fluid_IB_xyz_v, nfluid_IB_v, xf, yh, zf, itot, jtot, ktot, &
         diag_neighbs, periodic_x, periodic_y,  n_threads, &
         secfacids_v, secbndptids_v, secareas_v, bnddst_v, nfacsecs_v)
    write(*,*) 'Computation of facet_sections_v done.'

    ! w-grid
    call matchFacetsToCells(faces, faceNormal, n_fcts, verts, n_vert, &
         fluid_IB_w, solid_IB_w, fluid_IB_xyz_w, nfluid_IB_w, xf, yf, zh, itot, jtot, ktot, &
         diag_neighbs, periodic_x, periodic_y, n_threads, &
         secfacids_w, secbndptids_w, secareas_w, bnddst_w, nfacsecs_w)
    write(*,*) 'Computation of facet_sections_w done.'

    ! c-grid
    call matchFacetsToCells(faces, faceNormal, n_fcts, verts, n_vert, &
         fluid_IB_c, solid_IB_c, fluid_IB_xyz_c, nfluid_IB_c, xf, yf, zf, itot, jtot, ktot, &
         diag_neighbs, periodic_x, periodic_y, n_threads, &
         secfacids_c, secbndptids_c, secareas_c, bnddst_c, nfacsecs_c)
    write(*,*) 'Computation of facet_sections_c done.'

    !! Writing facet section files
    call print_facet_sections(nfacsecs_u, nfacsecs_v, nfacsecs_w, nfacsecs_c, &
                              secfacids_u, secareas_u, secbndptids_u, bnddst_u, &
                              secfacids_v, secareas_v, secbndptids_v, bnddst_v, &
                              secfacids_w, secareas_w, secbndptids_w, bnddst_w, &
                              secfacids_c, secareas_c, secbndptids_c, bnddst_c, &
                              'facet_sections_u.txt', 'facet_sections_v.txt', &
                              'facet_sections_w.txt', 'facet_sections_c.txt')
    write(*,*) 'Written facet_sections_*.txt files'


    !!!!!!! Writing number of line count file !!!!!!

    open(unit=100,file='info_fort.txt')
    write(unit=100,fmt='(A)') 'nfcts  nsolpts_u  nsolpts_v  nsolpts_w  nsolpts_c  ' // &
                              'nbndpts_u  nbndpts_v  nbndpts_w  nbndpts_c  ' // &
                              'nfctsecs_u  nfctsecs_v  nfctsecs_w  nfctsecs_c'
    write(unit=100,fmt='(i12,x,i12,x,i12,x,i12,x,i12,x,i12,x,i12,x,i12,x,i12,x,i12,x,i12,x,i12,x,i12)') &
         n_fcts, &
         count(solid_u), count(solid_v), count(solid_w), count(solid_c), &
         nfluid_IB_u, nfluid_IB_v, nfluid_IB_w, nfluid_IB_c, &
         nfacsecs_u, nfacsecs_v, nfacsecs_w, nfacsecs_c
    close(unit=100)


    !!!!!!! Memory clean up !!!!!!

    deallocate(xf,xh,yf,yh,zf,zh)
    deallocate(vertices,facets,incenters,faceNormals,verts,faces,faceNormal)
    deallocate(solid_u,solid_v,solid_w,solid_c)
    deallocate(fluid_IB_u, fluid_IB_v, fluid_IB_w, fluid_IB_c, &
               solid_IB_u, solid_IB_v, solid_IB_w, solid_IB_c)

    end_time = OMP_GET_WTIME()
    write(*,'(A,F10.3,A)') 'Elapsed time by IBM fortran routine: ', end_time - start_time, ' seconds.'

end program IBM_preproc
