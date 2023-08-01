module ibm_necessary_functions

    use omp_lib

    contains


    !! Subroutine to read STL data from text files vertices_file and facets_file
    subroutine read_data(vertices_file,n_vert,facets_file,n_fcts,zf_file,zh_file,ktot, &
        vertices,facets,incenters,faceNormals,zf,zh)
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
        15 format(i8,x,i8,x,i8,x,f15.10,x,f15.10,x,f15.10,x,f15.10,x,f15.10,x,f15.10)
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



    subroutine print_solid_flags(itot,jtot,ktot,solid_u,solid_v,solid_w,solid_c)
        implicit none

        integer, intent(in) :: itot, jtot, ktot
        logical, dimension(itot,jtot,ktot), intent(in) :: solid_u, solid_v, solid_w, solid_c

        integer :: ix, iy, iz

        !$ call OMP_SET_NUM_THREADS(4)
        !$OMP parallel
            !$OMP sections private(ix,iy,iz)
                
                !$OMP section
                open(unit=1,file='flag_u.txt')
                do iy = 1,jtot
                    do iz = 1,ktot
                        do ix = 1,itot
                            if (solid_u(ix,iy,iz)) then
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
                            if (solid_v(ix,iy,iz)) then
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
                            if (solid_w(ix,iy,iz)) then
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
                            if (solid_c(ix,iy,iz)) then
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

    end subroutine print_solid_flags



    subroutine print_solid_points_index(itot,jtot,ktot,solid_u,solid_v,solid_w,solid_c)
        implicit none

        integer, intent(in) :: itot, jtot, ktot
        logical, dimension(itot,jtot,ktot), intent(in) :: solid_u, solid_v, solid_w, solid_c

        integer :: ix, iy, iz

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
                write(*,*) 'Written solid_u.txt'

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
                write(*,*) 'Written solid_v.txt'

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
                write(*,*) 'Written solid_w.txt'

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
                write(*,*) 'Written solid_c.txt'

            !$OMP end sections
        !$OMP end parallel

    end subroutine print_solid_points_index



    subroutine boundaryMasks(uvwc, itot, jtot, ktot, solid_x, include_diagonals, stl_ground)
        implicit none

        integer, intent(in) :: itot, jtot, ktot
        logical, dimension(itot,jtot,ktot), intent(in) :: solid_x
        logical, intent(in) :: include_diagonals, stl_ground
        character, intent(in) :: uvwc

        logical, allocatable, dimension(:,:,:) :: fluid_IB, solid_IB, fluid_x, solid_x_w
        ! logical, dimension(itot,jtot,ktot) :: fluid_IB
        ! logical, dimension(itot,jtot,ktot) :: solid_IB
        integer :: ix,iy,iz
        character*20 :: filename
        character*14 :: filename2

        allocate(fluid_IB(itot,jtot,ktot))
        allocate(solid_IB(itot,jtot,ktot))
        allocate(fluid_x(itot,jtot,ktot))

        do iy = 1,jtot
            do iz = 1,ktot
                do ix = 1,itot
                    fluid_x(ix,iy,iz) = .not.(solid_x(ix,iy,iz))
                end do
            end do
        end do

        if (uvwc=='w') then

            allocate(solid_x_w(itot,jtot,ktot))
            solid_x_w = solid_x
            if (stl_ground) then
                do iy = 1,jtot
                    do ix = 1,itot
                        solid_x_w(ix,iy,1) = .true. !! % Bottom is always solid for w
                    end do
                end do
            end if

            call getBoundaryCells(fluid_IB, solid_IB, itot, jtot, ktot, fluid_x, solid_x_w, include_diagonals)

            do iy = 1,jtot
                do ix = 1,itot
                    fluid_IB(ix,iy,1) = .false. !! % Bottom is always solid for w
                end do
            end do
            deallocate(solid_x_w)

        else        !! For u, v and c

            call getBoundaryCells(fluid_IB, solid_IB, itot, jtot, ktot, fluid_x, solid_x, include_diagonals)

            if (stl_ground) then
                do iy = 1,jtot
                    do ix = 1,itot
                        if (.not.(solid_x(ix,iy,1))) then
                            fluid_IB(ix,iy,1) = .true.
                        end if
                    end do
                end do
            end if

        end if
        
        !$ call OMP_SET_NUM_THREADS(4)
        !$OMP parallel
            !$OMP sections private(ix,iy,iz,filename,filename2)

            !$OMP section
            filename = 'fluid_boundary_' // uvwc // '.txt'
            open(unit=1,file=filename)
                write(unit=1,fmt='(a18)') '# position (i,j,k)'
                do iy = 1,jtot
                    do iz = 1,ktot
                        do ix = 1,itot
                            if (fluid_IB(ix,iy,iz)) then
                                    write(unit=1,fmt='(i4,x,i4,x,i4,A)',advance='no') ix, iy, iz, NEW_LINE('a')
                            end if
                        end do
                    end do
                end do
            close(unit=1)
            write(*,*) 'Written ', filename

            !$OMP section
            filename = 'solid_boundary_' // uvwc // '.txt'
            open(unit=2,file=filename)
                write(unit=2,fmt='(a18)') '# position (i,j,k)'
                do iy = 1,jtot
                    do iz = 1,ktot
                        do ix = 1,itot
                            if (solid_IB(ix,iy,iz)) then
                                    write(unit=2,fmt='(i4,x,i4,x,i4,A)',advance='no') ix, iy, iz, NEW_LINE('a')
                            end if
                        end do
                    end do
                end do
            close(unit=2)
            
            !$OMP section
            filename2 = 'fluid_IB_' // uvwc // '.txt'
            open(unit=3,file=filename2)
            do iy = 1,jtot
                do iz = 1,ktot
                    do ix = 1,itot
                        if (fluid_IB(ix,iy,iz)) then
                                write(unit=3,fmt='(i1,A)',advance='no') 1, NEW_LINE('a')
                        else
                                write(unit=3,fmt='(i1,A)',advance='no') 0, NEW_LINE('a')
                        end if
                    end do
                end do
            end do
            close(unit=3)

            !$OMP section
            filename2 = 'solid_IB_' // uvwc // '.txt'
            open(unit=4,file=filename2)
            do iy = 1,jtot
                do iz = 1,ktot
                    do ix = 1,itot
                        if (solid_IB(ix,iy,iz)) then
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

        deallocate(fluid_IB,solid_IB,fluid_x)
    end subroutine boundaryMasks



    subroutine getBoundaryCells(fluid_IB, solid_IB, itot, jtot, ktot, fluid, solid, include_diagonals)
        implicit none

        integer, intent(in) :: itot, jtot, ktot
        logical, dimension(itot,jtot,ktot), intent(in) :: fluid, solid
        logical, intent(in) :: include_diagonals

        logical, dimension(itot,jtot,ktot), intent(out) :: fluid_IB
        logical, dimension(itot,jtot,ktot), intent(out) :: solid_IB

        integer :: ix, iy, iz

        ! Initializing with false
        do iy = 1,jtot
            do iz = 1,ktot
                do ix = 1,itot
                    fluid_IB(ix,iy,iz) = .false.
                    solid_IB(ix,iy,iz) = .false.
                end do
            end do
        end do

        do iy = 1,jtot
            do iz = 1,ktot
                do ix = 1,itot
                    if (fluid(ix,iy,iz)) then

                        ! Identify fluid IB points
                        if (ix/=1) then
                            if (solid(ix-1,iy,iz)) then
                                fluid_IB(ix,iy,iz) = .true.
                            end if
                        end if
                        if (ix/=itot) then
                            if (solid(ix+1,iy,iz)) then
                                fluid_IB(ix,iy,iz) = .true.
                            end if
                        end if
                        if (iy/=1) then
                            if (solid(ix,iy-1,iz)) then
                                fluid_IB(ix,iy,iz) = .true.
                            end if
                        end if
                        if (iy/=jtot) then
                            if (solid(ix,iy+1,iz)) then
                                fluid_IB(ix,iy,iz) = .true.
                            end if
                        end if
                        if (iz/=1) then
                            if (solid(ix,iy,iz-1)) then
                                fluid_IB(ix,iy,iz) = .true.
                            end if
                        end if
                        if (iz/=ktot) then
                            if (solid(ix,iy,iz+1)) then
                                fluid_IB(ix,iy,iz) = .true.
                            end if
                        end if

                        if (include_diagonals) then

                            if (ix/=1 .and. iy/=1) then
                                if (solid(ix-1,iy-1,iz)) then
                                    fluid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=1 .and. iy/=jtot) then
                                if (solid(ix-1,iy+1,iz)) then
                                    fluid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=itot .and. iy/=1) then
                                if (solid(ix+1,iy-1,iz)) then
                                    fluid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=itot .and. iy/=jtot) then
                                if (solid(ix+1,iy+1,iz)) then
                                    fluid_IB(ix,iy,iz) = .true.
                                end if
                            end if

                            if (ix/=1 .and. iz/=1) then
                                if (solid(ix-1,iy,iz-1)) then
                                    fluid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=1 .and. iz/=ktot) then
                                if (solid(ix-1,iy,iz+1)) then
                                    fluid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=itot .and. iz/=1) then
                                if (solid(ix+1,iy,iz-1)) then
                                    fluid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=itot .and. iz/=ktot) then
                                if (solid(ix+1,iy,iz+1)) then
                                    fluid_IB(ix,iy,iz) = .true.
                                end if
                            end if

                            if (iy/=1 .and. iz/=1) then
                                if (solid(ix,iy-1,iz-1)) then
                                    fluid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (iy/=1 .and. iz/=ktot) then
                                if (solid(ix,iy-1,iz+1)) then
                                    fluid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (iy/=jtot .and. iz/=1) then
                                if (solid(ix,iy+1,iz-1)) then
                                    fluid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (iy/=jtot .and. iz/=ktot) then
                                if (solid(ix,iy+1,iz+1)) then
                                    fluid_IB(ix,iy,iz) = .true.
                                end if
                            end if

                            if (ix/=1 .and. iy/=1 .and. iz/=1) then
                                if (solid(ix-1,iy-1,iz-1)) then
                                    fluid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=itot .and. iy/=1 .and. iz/=1) then
                                if (solid(ix+1,iy-1,iz-1)) then
                                    fluid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=1 .and. iy/=jtot .and. iz/=1) then
                                if (solid(ix-1,iy+1,iz-1)) then
                                    fluid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=itot .and. iy/=jtot .and. iz/=1) then
                                if (solid(ix+1,iy+1,iz-1)) then
                                    fluid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=1 .and. iy/=1 .and. iz/=ktot) then
                                if (solid(ix-1,iy-1,iz+1)) then
                                    fluid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=itot .and. iy/=1 .and. iz/=ktot) then
                                if (solid(ix+1,iy-1,iz+1)) then
                                    fluid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=1 .and. iy/=jtot .and. iz/=ktot) then
                                if (solid(ix-1,iy+1,iz+1)) then
                                    fluid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=itot .and. iy/=jtot .and. iz/=ktot) then
                                if (solid(ix+1,iy+1,iz+1)) then
                                    fluid_IB(ix,iy,iz) = .true.
                                end if
                            end if

                        end if

                    else
                        ! Identify solid IB points
                        if (ix/=1) then
                            if (fluid(ix-1,iy,iz)) then
                                solid_IB(ix,iy,iz) = .true.
                            end if
                        end if
                        if (ix/=itot) then
                            if (fluid(ix+1,iy,iz)) then
                                solid_IB(ix,iy,iz) = .true.
                            end if
                        end if
                        if (iy/=1) then
                            if (fluid(ix,iy-1,iz)) then
                                solid_IB(ix,iy,iz) = .true.
                            end if
                        end if
                        if (iy/=jtot) then
                            if (fluid(ix,iy+1,iz)) then
                                solid_IB(ix,iy,iz) = .true.
                            end if
                        end if
                        if (iz/=1) then
                            if (fluid(ix,iy,iz-1)) then
                                solid_IB(ix,iy,iz) = .true.
                            end if
                        end if
                        if (iz/=ktot) then
                            if (fluid(ix,iy,iz+1)) then
                                solid_IB(ix,iy,iz) = .true.
                            end if
                        end if

                        if (include_diagonals) then

                            if (ix/=1 .and. iy/=1) then
                                if (fluid(ix-1,iy-1,iz)) then
                                    solid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=1 .and. iy/=jtot) then
                                if (fluid(ix-1,iy+1,iz)) then
                                    solid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=itot .and. iy/=1) then
                                if (fluid(ix+1,iy-1,iz)) then
                                    solid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=itot .and. iy/=jtot) then
                                if (fluid(ix+1,iy+1,iz)) then
                                    solid_IB(ix,iy,iz) = .true.
                                end if
                            end if

                            if (ix/=1 .and. iz/=1) then
                                if (fluid(ix-1,iy,iz-1)) then
                                    solid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=1 .and. iz/=ktot) then
                                if (fluid(ix-1,iy,iz+1)) then
                                    solid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=itot .and. iz/=1) then
                                if (fluid(ix+1,iy,iz-1)) then
                                    solid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=itot .and. iz/=ktot) then
                                if (fluid(ix+1,iy,iz+1)) then
                                    solid_IB(ix,iy,iz) = .true.
                                end if
                            end if

                            if (iy/=1 .and. iz/=1) then
                                if (fluid(ix,iy-1,iz-1)) then
                                    solid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (iy/=1 .and. iz/=ktot) then
                                if (fluid(ix,iy-1,iz+1)) then
                                    solid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (iy/=jtot .and. iz/=1) then
                                if (fluid(ix,iy+1,iz-1)) then
                                    solid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (iy/=jtot .and. iz/=ktot) then
                                if (fluid(ix,iy+1,iz+1)) then
                                    solid_IB(ix,iy,iz) = .true.
                                end if
                            end if

                            if (ix/=1 .and. iy/=1 .and. iz/=1) then
                                if (fluid(ix-1,iy-1,iz-1)) then
                                    solid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=itot .and. iy/=1 .and. iz/=1) then
                                if (fluid(ix+1,iy-1,iz-1)) then
                                    solid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=1 .and. iy/=jtot .and. iz/=1) then
                                if (fluid(ix-1,iy+1,iz-1)) then
                                    solid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=itot .and. iy/=jtot .and. iz/=1) then
                                if (fluid(ix+1,iy+1,iz-1)) then
                                    solid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=1 .and. iy/=1 .and. iz/=ktot) then
                                if (fluid(ix-1,iy-1,iz+1)) then
                                    solid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=itot .and. iy/=1 .and. iz/=ktot) then
                                if (fluid(ix+1,iy-1,iz+1)) then
                                    solid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=1 .and. iy/=jtot .and. iz/=ktot) then
                                if (fluid(ix-1,iy+1,iz+1)) then
                                    solid_IB(ix,iy,iz) = .true.
                                end if
                            end if
                            if (ix/=itot .and. iy/=jtot .and. iz/=ktot) then
                                if (fluid(ix+1,iy+1,iz+1)) then
                                    solid_IB(ix,iy,iz) = .true.
                                end if
                            end if

                        end if

                    end if
                end do
            end do
        end do

    end subroutine getBoundaryCells

    
end module ibm_necessary_functions
