program Xie
    
    use omp_lib
    use in_mypoly_functions
    
    implicit none
    
    character*3 :: expnr, dummy
    real :: max_height, L_char, tol, Ray_dir_u(3), Ray_dir_v(3), Ray_dir_w(3), Ray_dir_c(3)
    integer :: n_vert, n_fcts
    real :: dx, dy  !, dz
    integer :: itot, jtot, ktot
    integer :: i, ix, iy, iz
    real, allocatable, dimension(:) :: xf, xh, yf, yh, zf, zh, vertices, incenters, faceNormals
    integer, allocatable, dimension(:) :: facets
    logical, allocatable, dimension(:,:,:) :: solid_u, solid_v, solid_w, solid_c
    

    open(unit=50,file='inmypoly_inp_info.txt')
    read(unit=50,fmt='(3a)') expnr
    read(unit=50,fmt='(f15.10,x,f15.10)') dx, dy  !, dz
    read(unit=50,fmt='(i5,x,i5,x,i5)') itot, jtot, ktot
    read(unit=50,fmt='(f15.10)') tol
    read(unit=50,fmt='(f15.10,x,f15.10,x,f15.10)') Ray_dir_u(1), Ray_dir_u(2), Ray_dir_u(3)
    read(unit=50,fmt='(f15.10,x,f15.10,x,f15.10)') Ray_dir_v(1), Ray_dir_v(2), Ray_dir_v(3)
    read(unit=50,fmt='(f15.10,x,f15.10,x,f15.10)') Ray_dir_w(1), Ray_dir_w(2), Ray_dir_w(3)
    read(unit=50,fmt='(f15.10,x,f15.10,x,f15.10)') Ray_dir_c(1), Ray_dir_c(2), Ray_dir_c(3)
    read(unit=50,fmt='(i5,x,i5)') n_vert, n_fcts
    close(unit=50)


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

    call read_data('vertices.txt',n_vert,'Stl_data.txt',n_fcts,'zfgrid.txt','zhgrid.txt',ktot, &
                    vertices,facets,incenters,faceNormals,zf,zh)
    
    max_height = MAXVAL(vertices(3:n_vert*3:3)) + tol
    L_char = max_facet_side(n_vert,vertices,n_fcts,facets) + tol
    
    solid_u = is_grid_in_mypoly_func(n_vert,vertices,n_fcts,facets,incenters,faceNormals, &
                                    itot,xh,jtot,yf,ktot,zf,Ray_dir_u,L_char,max_height,tol)
    solid_v = is_grid_in_mypoly_func(n_vert,vertices,n_fcts,facets,incenters,faceNormals, &
                                    itot,xf,jtot,yh,ktot,zf,Ray_dir_v,L_char,max_height,tol)
    solid_w = is_grid_in_mypoly_func(n_vert,vertices,n_fcts,facets,incenters,faceNormals, &
                                    itot,xf,jtot,yf,ktot,zh,Ray_dir_w,L_char,max_height,tol)
    solid_c = is_grid_in_mypoly_func(n_vert,vertices,n_fcts,facets,incenters,faceNormals, &
                                    itot,xf,jtot,yf,ktot,zf,Ray_dir_c,L_char,max_height,tol)
    
    
    !$ call OMP_SET_NUM_THREADS(4)
    !$OMP parallel
        !$OMP sections private(ix,iy,iz)
            
            !$OMP section
            open(unit=1,file='flag_u.txt')
            do iy = 1,jtot
                do iz = 1,ktot
                    do ix = 1,itot
                        if (solid_u(ix,iy,iz) .eqv. .true.) then
                                write(unit=1,fmt='(i1,A)',advance='no') 1, NEW_LINE('a')
                        else
                                write(unit=1,fmt='(i1,A)',advance='no') 0, NEW_LINE('a')
                        end if
                    end do
                end do
            end do
            close(unit=1)
            
            !$OMP section
            open(unit=2,file='flag_v.txt')
            do iy = 1,jtot
                do iz = 1,ktot
                    do ix = 1,itot
                        if (solid_v(ix,iy,iz) .eqv. .true.) then
                                write(unit=2,fmt='(i1,A)',advance='no') 1, NEW_LINE('a')
                        else
                                write(unit=2,fmt='(i1,A)',advance='no') 0, NEW_LINE('a')
                        end if
                    end do
                end do
            end do
            close(unit=2)
            
            !$OMP section
            open(unit=3,file='flag_w.txt')
            do iy = 1,jtot
                do iz = 1,ktot
                    do ix = 1,itot
                        if (solid_w(ix,iy,iz) .eqv. .true.) then
                                write(unit=3,fmt='(i1,A)',advance='no') 1, NEW_LINE('a')
                        else
                                write(unit=3,fmt='(i1,A)',advance='no') 0, NEW_LINE('a')
                        end if
                    end do
                end do
            end do
            close(unit=3)
            
            !$OMP section
            open(unit=4,file='flag_c.txt')
            do iy = 1,jtot
                do iz = 1,ktot
                    do ix = 1,itot
                        if (solid_c(ix,iy,iz) .eqv. .true.) then
                                write(unit=4,fmt='(i1,A)',advance='no') 1, NEW_LINE('a')
                        else
                                write(unit=4,fmt='(i1,A)',advance='no') 0, NEW_LINE('a')
                        end if
                    end do
                end do
            end do
            close(unit=4)

        !$OMP end sections
    !$OMP end parallel

    !$ call OMP_SET_NUM_THREADS(4)
    !$OMP parallel
        !$OMP sections private(ix,iy,iz)

            !$OMP section
            open(unit=1,file='solid_u.txt')
            write(unit=1,fmt='(a18)') '# position (i,j,k)'
            do iy = 1,jtot
                do iz = 1,ktot
                    do ix = 1,itot
                        if (solid_u(ix,iy,iz)) then
                                write(unit=1,fmt='(i4,x,i4,x,i4,A)',advance='no') ix, iy, iz, NEW_LINE('a')
                        end if
                    end do
                end do
            end do
            close(unit=1)

            !$OMP section
            open(unit=2,file='solid_v.txt')
            write(unit=2,fmt='(a18)') '# position (i,j,k)'
            do iy = 1,jtot
                do iz = 1,ktot
                    do ix = 1,itot
                        if (solid_v(ix,iy,iz)) then
                                write(unit=2,fmt='(i4,x,i4,x,i4,A)',advance='no') ix, iy, iz, NEW_LINE('a')
                        end if
                    end do
                end do
            end do
            close(unit=2)

            !$OMP section
            open(unit=3,file='solid_w.txt')
            write(unit=3,fmt='(a18)') '# position (i,j,k)'
            do iy = 1,jtot
                do iz = 1,ktot
                    do ix = 1,itot
                        if (solid_w(ix,iy,iz)) then
                                write(unit=3,fmt='(i4,x,i4,x,i4,A)',advance='no') ix, iy, iz, NEW_LINE('a')
                        end if
                    end do
                end do
            end do
            close(unit=3)

            !$OMP section
            open(unit=4,file='solid_c.txt')
            write(unit=4,fmt='(a18)') '# position (i,j,k)'
            do iy = 1,jtot
                do iz = 1,ktot
                    do ix = 1,itot
                        if (solid_c(ix,iy,iz)) then
                                write(unit=4,fmt='(i4,x,i4,x,i4,A)',advance='no') ix, iy, iz, NEW_LINE('a')
                        end if
                    end do
                end do
            end do
            close(unit=4)

        !$OMP end sections
    !$OMP end parallel
 
    deallocate(xf,xh,yf,yh,zf,zh,vertices,facets,incenters,faceNormals,solid_u,solid_v,solid_w,solid_c)
    
end program Xie



!! Subroutine to read STL data from text files vertices_file and facets_file
subroutine read_data(vertices_file,n_vert,facets_file,n_fcts,zf_file,zh_file,ktot, &
                     vertices,facets,incenters,faceNormals,zf,zh)
    
    use omp_lib
    
    implicit none

    integer, intent(in) :: n_vert, n_fcts, ktot
    character*12, intent(in) :: vertices_file, facets_file
    character*10, intent(in) :: zf_file, zh_file
    real, dimension(n_vert*3), intent(out) :: vertices
    integer, dimension(n_fcts*3), intent(out) :: facets
    real, dimension(n_fcts*3), intent(out) :: incenters, faceNormals
    real, intent(out) :: zf(ktot), zh(ktot)

    integer :: i, kk

    !$ call OMP_SET_NUM_THREADS(4)
    !$OMP parallel
        !$OMP sections private(i,kk)

        !$OMP section
        open(unit=1,file=vertices_file)
        do i=1,n_vert
            read(unit=1,fmt='(f15.10,x,f15.10,x,f15.10)') vertices(3*i-2), vertices(3*i-1), vertices(3*i)
        end do
        close(unit=1)

        !$OMP section
        15 format(i5,x,i5,x,i5,x,f15.10,x,f15.10,x,f15.10,x,f15.10,x,f15.10,x,f15.10)
        open(unit=2,file=facets_file)
        kk=1
        do i=1,n_fcts
            read(unit=2,fmt=15) facets(kk), facets(kk+1), facets(kk+2), &
                                incenters(kk), incenters(kk+1), incenters(kk+2), &
                                faceNormals(kk), faceNormals(kk+1), faceNormals(kk+2)
            kk = kk+3
        end do
        close(unit=2)

        !$OMP section
        open(unit=3,file=zf_file)
        do i = 1,ktot    
            read(unit=3,fmt='(f15.10)') zf(i)
        end do
        close(unit=3)

        !$OMP section
        open(unit=4,file=zh_file)
            do i = 1,ktot   
                read(unit=4,fmt='(f15.10)') zh(i)
            end do
        close(unit=4)

        !$OMP end sections
    !$OMP end parallel

end subroutine read_data



!! Subroutine to print STL data
subroutine print_stl_data(n_vert,vertices,n_fcts,facets,incenters,faceNormals)
    implicit none

    integer, intent(in) :: n_vert, n_fcts
    real, dimension(n_vert*3), intent(in) :: vertices
    integer, dimension(n_fcts*3), intent(in) :: facets
    real, dimension(n_fcts*3), intent(in) :: incenters, faceNormals

    integer :: i, kk

    10 format(f15.10,x,f15.10,x,f15.10)

    kk=1
    print *,"VERTICES"
    do i=1,n_vert
        write(*,fmt=10,advance='no') vertices(kk), vertices(kk+1), vertices(kk+2)
        write(*,*)
        kk=kk+3
    end do

    kk=1
    print *,"FACETS"
    do i=1,n_fcts
        write(*,fmt='(i5,x,i5,x,i5)',advance='no') facets(kk), facets(kk+1), facets(kk+2)
        write(*,*)
        kk=kk+3
    end do

    kk=1
    print *,"IN-CENTERS"
    do i=1,n_fcts
        write(*,fmt=10,advance='no') incenters(kk), incenters(kk+1), incenters(kk+2)
        write(*,*)
        kk=kk+3
    end do

    kk=1
    print *,"FACE NORMALS"
    do i=1,n_fcts
        write(*,fmt=10,advance='no') faceNormals(kk), faceNormals(kk+1), faceNormals(kk+2)
        write(*,*)
        kk=kk+3
    end do

end subroutine print_stl_data
