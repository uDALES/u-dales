subroutine run_ibm_preproc_f2py(n_vert, n_fcts, vertices, facets, incenters, faceNormals, &
                                dx, dy, itot, jtot, ktot, zt, zm, &
                                Ray_dir_u, Ray_dir_v, Ray_dir_w, Ray_dir_c, &
                                stl_ground, diag_neighbs, periodic_x, periodic_y, &
                                n_threads, tol, counts)
    use ibm_preproc_mod, only: run_ibm_preproc_core
    implicit none

    integer, intent(in)                      :: n_vert, n_fcts
    real   , intent(in), dimension(n_vert,3) :: vertices
    integer, intent(in), dimension(n_fcts,3) :: facets
    real   , intent(in), dimension(n_fcts,3) :: incenters, faceNormals
    real   , intent(in)                      :: dx, dy
    integer, intent(in)                      :: itot, jtot, ktot
    real   , intent(in), dimension(ktot)     :: zt, zm
    real   , intent(in), dimension(3)        :: Ray_dir_u, Ray_dir_v, Ray_dir_w, Ray_dir_c
    logical, intent(in)                      :: stl_ground, diag_neighbs, periodic_x, periodic_y
    integer, intent(in)                      :: n_threads
    real   , intent(in)                      :: tol

    integer, intent(out), dimension(13)      :: counts

    call run_ibm_preproc_core(n_vert, n_fcts, vertices, facets, incenters, faceNormals, &
                              dx, dy, itot, jtot, ktot, zt, zm, &
                              Ray_dir_u, Ray_dir_v, Ray_dir_w, Ray_dir_c, &
                              stl_ground, diag_neighbs, periodic_x, periodic_y, &
                              n_threads, tol, counts)
end subroutine run_ibm_preproc_f2py
