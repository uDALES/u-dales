module definitions
  implicit none
  save

  ! Shared staggering/location selectors for domain-aware reductions.
  !
  ! Reduction naming rule:
  ! routines named reduce_* are assumed to return globally reduced results
  ! unless a *_local suffix is used explicitly.
  integer, parameter :: LOC_C   = 7, &
                     &  LOC_U   = 6, &
                     &  LOC_V   = 5, &
                     &  LOC_W   = 3, &
                     &  LOC_UV  = 4, &
                     &  LOC_WU  = 2, &
                     &  LOC_VW  = 1, &
                     &  LOC_UVW = 0
end module definitions
