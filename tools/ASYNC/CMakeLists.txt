# Copyright (c) 2016, Technische Universitaet Muenchen
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

cmake_minimum_required( VERSION 2.8 )

project( Async )

option( USE_MPI "Compile with MPI" ON )

# Enable warnings
add_definitions( -Wall )

# MPI
if( USE_MPI )
	# Find MPI
	find_package( MPI REQUIRED )
	if( MPI_LIBRARY )
		include_directories( SYSTEM ${MPI_C_INCLUDE_PATH} )
	else( MPI_LIBRARY )
		#Workaround for stupid mpicc/mpicxx
		get_filename_component( MPI_BASE ${MPI_C_COMPILER} PATH )
		set( MPI_BASE ${MPI_BASE}/.. )
		include_directories( SYSTEM ${MPI_BASE}/include )
		link_directories( ${MPI_BASE}/lib )
	endif( MPI_LIBRARY )
	set( CMAKE_C_COMPILER ${MPI_C_COMPILER} )
	set( CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} )
	set( CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} )

	# Do not include OpenMPI C++ Bindings
	add_definitions( -DOMPI_SKIP_MPICXX )

	# Use MPI
	add_definitions( -DUSE_MPI )
endif( USE_MPI )

# Threads
find_package( Threads REQUIRED )

# Enable testing
enable_testing()

# Make sure we find the header files
include_directories( ${CMAKE_SOURCE_DIR} )

# Submodules
include_directories( ${CMAKE_SOURCE_DIR}/submodules )

add_subdirectory( tests )