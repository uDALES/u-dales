!> \file modchem.f90
!!tg3315, 2 Nov 2017 

!> Input chemistry of the null cycle into DALES following Zhong, Cai 2015.
!> Assumes scalar fields are in order 1: NO, 2: NO2 and 3: O3.

module modchem
implicit none
save

contains
    subroutine chem
    use modglobal,  only : lchem,k1,JNO2,dt,rk3step,ib,ie,ihc,ih,jb,je,jhc,jh,kb,ke,khc,kh
    use modfields,  only : svp,svm,sv0,IIc
    implicit none
    real, dimension(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc)     :: dummyNO, dummyNO2, dummyO3

    if (lchem .eqv. .false.) return

    ! forward Euler

    ! [NO]
!    svp(:,:,:,1) = svp(:,:,:,1) + 30.006 * ( JNO2 * (svm(:,:,:,2)/46.005) - k1 * (svm(:,:,:,1) / 30.006) * (svm(:,:,:,3)/47.997) )

    ! [NO2]
!    svp(:,:,:,2) = svp(:,:,:,2) + 46.005 * (-JNO2 * (svm(:,:,:,2)/46.005) + k1 * (svm(:,:,:,1) / 30.006) * (svm(:,:,:,3)/47.997) )

    ! [O3]
!    svp(:,:,:,3) = svp(:,:,:,3) + 47.997 * ( JNO2 * (svm(:,:,:,2)/46.005) - k1 * (svm(:,:,:,1) / 30.006) * (svm(:,:,:,3)/47.997) )


    if (.not. rk3step==3)  return

    ! convert into mol/m^3
    dummyNO = IIc*sv0(:,:,kb:ke+khc,1)/30.006
    dummyNO2= IIc*sv0(:,:,kb:ke+khc,2)/46.005
    dummyO3 = IIc*sv0(:,:,kb:ke+khc,3)/47.997

    !backward Euler, semi-implicit
!    sv0(:,:,kb:ke+khc,1) = 30.006 * ( ( (sv0(:,:,kb:ke+khc,1)/30.006) + JNO2 * dummyNO2 * dt) / (1. + k1 * dummyO3 * dt) )

!    sv0(:,:,kb:ke+khc,2) = 46.005 * ( ( (sv0(:,:,kb:ke+khc,2)/46.005) + k1 * dummyNO * dummyO3 * dt )  / (1. + JNO2 * dt) )

!    sv0(:,:,kb:ke+khc,3) = 47.997 * ( ( (sv0(:,:,kb:ke+khc,3)/47.997) + JNO2 * dummyNO2 * dt) / (1. + k1 * dummyNO * dt) )

    !backward Euler, fully implicit. Derivation at (/projects/Chemistry/FullyImplicit.mw)
       sv0(:,:,kb:ke+khc,1) = 30.006 * ( (sv0(:,:,kb:ke+khc,1)/30.006) + ( dt * (- k1 * dummyNO * dummyO3 + JNO2 * dummyNO2) ) / &
                              ( 1. + ( ( dummyNO + dummyO3 ) * k1 + JNO2 ) * dt ) )

       sv0(:,:,kb:ke+khc,2) = 46.005 * ( (sv0(:,:,kb:ke+khc,2)/46.005) - ( dt * (- k1 * dummyNO * dummyO3 + JNO2 * dummyNO2) ) / &
                              ( 1. + ( ( dummyNO + dummyO3 ) * k1 + JNO2 ) * dt ) ) 

       sv0(:,:,kb:ke+khc,3) = 47.997 * ( (sv0(:,:,kb:ke+khc,3)/47.997) + ( dt * (- k1 * dummyNO * dummyO3 + JNO2 * dummyNO2) ) / &
                              ( 1. + ( ( dummyNO + dummyO3 ) * k1 + JNO2 ) * dt ) )
 
    ! alternative method in Zhong 2017 eqn. 7, analytical solution!

    end subroutine chem

  end module modchem
