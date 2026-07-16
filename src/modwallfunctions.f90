!  This file is part of uDALES.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 2006-2021 the uDALES Team.
!
! This module is an empty shell: the flat-surface wall functions it carried
! (the Uno1995/Cai2012 scheme and its stability helpers; see BCcleanup_backlog.md,
! Phase 1) were deleted after their last caller, the legacy flat-bottom scheme,
! was removed. The live IBM wall functions are wallfunmom/wallfunheat in modibm.
! The file itself can be dropped in a follow-up (the build globs src/*.f90).
module modwallfunctions

   implicit none

end module modwallfunctions
