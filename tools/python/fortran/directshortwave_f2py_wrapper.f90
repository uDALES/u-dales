subroutine calculate_direct_shortwave_f2py(connectivityList, incenter, faceNormal, &
                                           vertices, nFaces, nVertices, nsun, &
                                           irradiance, resolution, Sdir)
    use directshortwave_mod, only: calculateDirectShortwave
    implicit none

    integer, intent(in) :: nFaces, nVertices
    integer, intent(in), dimension(nFaces, 3) :: connectivityList
    real   , intent(in), dimension(nFaces, 3) :: incenter, faceNormal
    real   , intent(in), dimension(nVertices, 3) :: vertices
    real   , intent(in) :: nsun(3), irradiance, resolution
    real   , intent(out), dimension(nFaces) :: Sdir

    call calculateDirectShortwave(connectivityList, incenter, faceNormal, &
                                  nFaces, vertices, nVertices, nsun, &
                                  irradiance, resolution, Sdir)
end subroutine calculate_direct_shortwave_f2py
