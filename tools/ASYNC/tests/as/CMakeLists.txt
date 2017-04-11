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

# Add tests
add_cxx_test( TestSync ${CMAKE_CURRENT_SOURCE_DIR}/Sync.t.h )
add_cxx_test( TestThread ${CMAKE_CURRENT_SOURCE_DIR}/Thread.t.h )
add_cxx_test( TestThreadBuffer ${CMAKE_CURRENT_SOURCE_DIR}/Thread_buffer.t.h )
set_tests_properties( TestThreadBuffer PROPERTIES ENVIRONMENT "ASYNC_BUFFER_ALIGNMENT=65536" )
if( USE_MPI )
    add_mpi_cxx_test( TestMPIScheduler ${CMAKE_CURRENT_SOURCE_DIR}/MPIScheduler.t.h 5 )
    add_mpi_cxx_test( TestMPI ${CMAKE_CURRENT_SOURCE_DIR}/MPI.t.h 5 )
    set_tests_properties( TestMPI PROPERTIES ENVIRONMENT "ASYNC_MPI_MAX_SEND=1048576" )
    add_mpi_cxx_test( TestMPIAsync ${CMAKE_CURRENT_SOURCE_DIR}/MPI.t.h 5 )
    set_target_properties( TestMPIAsync PROPERTIES COMPILE_DEFINITIONS "ASYNC_MPI_COPY" )
    set_tests_properties( TestMPIAsync PROPERTIES ENVIRONMENT "ASYNC_MPI_MAX_SEND=1048576;" )
endif( USE_MPI )