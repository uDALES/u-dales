program IBM_flagging

    use omp_lib
    use in_mypoly_functions
    use ibm_necessary_functions
    implicit none

    real :: max_height, L_char, tol, Ray_dir_u(3), Ray_dir_v(3), Ray_dir_w(3), Ray_dir_c(3)
    integer :: n_vert, n_fcts, n_threads, stl_ground_i, diag_neighbs_i
    logical :: stl_ground, diag_neighbs
    real :: dx, dy  !, dz
    integer :: itot, jtot, ktot
    integer :: i
    real, allocatable, dimension(:) :: xf, xh, yf, yh, zf, zh, vertices, incenters, faceNormals
    integer, allocatable, dimension(:) :: facets
    logical, allocatable, dimension(:,:,:) :: solid_u, solid_v, solid_w, solid_c
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
    read(unit=50,fmt='(i1,x,i1)') stl_ground_i, diag_neighbs_i
    close(unit=50)

    if (stl_ground_i==1) then
        stl_ground = .true.
    elseif (stl_ground_i==0) then
        stl_ground = .false.
    else
        print *, "Incorrect input for 'stl_ground'."
    end if

    if (diag_neighbs_i==1) then
        diag_neighbs = .true.
    elseif (diag_neighbs_i==0) then
        diag_neighbs = .false.
    else
        print *, "Incorrect input for 'diag_neighbs'."
    end if

    allocate(vertices(n_vert*3))
    allocate(facets(n_fcts*3))
    allocate(incenters(n_fcts*3))
    allocate(faceNormals(n_fcts*3))
    allocate(xf(itot))
    allocate(xh(itot))
    allocate(yf(jtot))
    allocate(yh(jtot))
    allocate(zf(ktot))
    allocate(zh(ktot))
    allocate(solid_u(itot,jtot,ktot))
    allocate(solid_v(itot,jtot,ktot))
    allocate(solid_w(itot,jtot,ktot))
    allocate(solid_c(itot,jtot,ktot))


    xf(1:itot) = [((i-1)*dx+(dx/2.0), i=1,itot)]
    xh(1:itot) = [((i-1)*dx, i=1,itot)]

    yf(1:jtot) = [((i-1)*dy+(dy/2.0), i=1,jtot)]
    yh(1:jtot) = [((i-1)*dy, i=1,jtot)]

    ! zf(1:ktot) = [((i-1)*dz+(dz/2.0), i=1,ktot)]
    ! zh(1:ktot) = [((i-1)*dz, i=1,ktot)]

    call read_data('vertices.txt',n_vert,'faces.txt',n_fcts,'zfgrid.txt','zhgrid.txt',ktot, &
                    vertices,facets,incenters,faceNormals,zf,zh)

    max_height = MAXVAL(vertices(3:n_vert*3:3)) + tol
    L_char = max_facet_side(n_vert,vertices,n_fcts,facets) + tol


    solid_u = is_grid_in_mypoly_func(n_vert,vertices,n_fcts,facets,incenters,faceNormals, &
                                    itot,xh,jtot,yf,ktot,zf,Ray_dir_u,L_char,max_height,tol,n_threads)
    write(*,*) 'Determined solid points for u-grid.'
    solid_v = is_grid_in_mypoly_func(n_vert,vertices,n_fcts,facets,incenters,faceNormals, &
                                    itot,xf,jtot,yh,ktot,zf,Ray_dir_v,L_char,max_height,tol,n_threads)
    write(*,*) 'Determined solid points for v-grid.'
    solid_w = is_grid_in_mypoly_func(n_vert,vertices,n_fcts,facets,incenters,faceNormals, &
                                    itot,xf,jtot,yf,ktot,zh,Ray_dir_w,L_char,max_height,tol,n_threads)
    write(*,*) 'Determined solid points for w-grid.'
    solid_c = is_grid_in_mypoly_func(n_vert,vertices,n_fcts,facets,incenters,faceNormals, &
                                    itot,xf,jtot,yf,ktot,zf,Ray_dir_c,L_char,max_height,tol,n_threads)
    write(*,*) 'Determined solid points for c-grid.'


    ! call print_solid_flags(itot,jtot,ktot,solid_u,solid_v,solid_w,solid_c)
    call print_solid_points_index(itot,jtot,ktot,solid_u,solid_v,solid_w,solid_c)


    call boundaryMasks('u', itot, jtot, ktot, solid_u, diag_neighbs, stl_ground)
    call boundaryMasks('v', itot, jtot, ktot, solid_v, diag_neighbs, stl_ground)
    call boundaryMasks('w', itot, jtot, ktot, solid_w, diag_neighbs, stl_ground)
    call boundaryMasks('c', itot, jtot, ktot, solid_c, diag_neighbs, stl_ground)


    deallocate(xf,xh,yf,yh,zf,zh,vertices,facets,incenters,faceNormals,solid_u,solid_v,solid_w,solid_c)

    end_time = OMP_GET_WTIME()
    write(*,'(A,F10.6,A)') 'Elapsed time by IBM solid-fluid tagging: ', end_time - start_time, ' seconds.'

end program IBM_flagging
