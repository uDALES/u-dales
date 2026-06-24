subroutine run_ibm_preproc_f2py(vertices, facets, incenters, faceNormals, zf, zh, &
                                dx, dy, itot, jtot, ktot, tol, &
                                Ray_dir_u, Ray_dir_v, Ray_dir_w, Ray_dir_c, &
                                n_threads, stl_ground_i, diag_neighbs_i, periodic_x_i, periodic_y_i, &
                                n_vert, n_fcts, counts)
    use ibm_preproc_mod, only: run_ibm_preproc_from_arrays
    implicit none

    integer, intent(in) :: n_vert, n_fcts, itot, jtot, ktot, n_threads
    integer, intent(in) :: stl_ground_i, diag_neighbs_i, periodic_x_i, periodic_y_i
    real   , intent(in) :: dx, dy, tol
    real   , intent(in), dimension(3) :: Ray_dir_u, Ray_dir_v, Ray_dir_w, Ray_dir_c
    real   , intent(in), dimension(n_vert, 3) :: vertices
    integer, intent(in), dimension(n_fcts, 3) :: facets
    real   , intent(in), dimension(n_fcts, 3) :: incenters, faceNormals
    real   , intent(in), dimension(ktot) :: zf, zh
    integer, intent(out), dimension(13) :: counts

    call run_ibm_preproc_from_arrays(vertices, facets, incenters, faceNormals, zf, zh, &
                                     dx, dy, itot, jtot, ktot, tol, &
                                     Ray_dir_u, Ray_dir_v, Ray_dir_w, Ray_dir_c, &
                                     n_threads, stl_ground_i, diag_neighbs_i, periodic_x_i, periodic_y_i, &
                                     n_vert, n_fcts, counts)
end subroutine run_ibm_preproc_f2py
