# uDALES (https://github.com/uDALES/u-dales).
# Copyright (C) 2019 D. Meyer.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

include(ExternalProject)

ExternalProject_Add(fishpack-cmake-project
    GIT_REPOSITORY https://github.com/dmey/fishpack-cmake.git
    GIT_TAG dmey/cray-support
    UPDATE_DISCONNECTED ${SKIP_UPDATE_EXTERNAL_PROJECTS}
    CMAKE_ARGS
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DUSE_REAL8=ON
        -DTESTS=OFF
    INSTALL_COMMAND ""
)

ExternalProject_Add(vfftpack-cmake-project
    GIT_REPOSITORY https://github.com/dmey/vfftpack-cmake.git
    GIT_TAG dmey/cray-support
    UPDATE_DISCONNECTED ${SKIP_UPDATE_EXTERNAL_PROJECTS}
    CMAKE_ARGS
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DUSE_REAL8=ON
    INSTALL_COMMAND ""
)
