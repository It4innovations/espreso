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

# Add cxxtest pathes
set( CXXTEST_PATH ${CMAKE_SOURCE_DIR}/submodules/cxxtest )
set( CMAKE_INCLUDE_PATH
    ${CMAKE_SOURCE_DIR}/submodules/cxxtest
    ${CMAKE_INCLUDE_PATH} )
set( CMAKE_MODULE_PATH
	${CMAKE_SOURCE_DIR}/submodules/cxxtest/build_tools/cmake
	${CMAKE_MODULE_PATH} )

# Get cxxtest
find_package( CxxTest REQUIRED )
include_directories( ${CXXTEST_INCLUDES} )

# Timout for tests
set( test_timeout 20 )

# Add test functions
function( add_cxx_test target source )
    cxx_test( ${target} ${source} )
    target_link_libraries( ${target} ${CMAKE_THREAD_LIBS_INIT} )
    set_tests_properties( ${target} PROPERTIES TIMEOUT ${test_timeout} ) # Detect deadlocks
endfunction( add_cxx_test )

set( cxx_template ${CMAKE_CURRENT_SOURCE_DIR}/mpirunner.tpl )
function( add_mpi_cxx_test target source ranks )
    get_filename_component(CPP_FILE_NAME ${source} NAME)
    string(REGEX REPLACE "h$|hpp$" "cpp" CPP_FILE_NAME ${CPP_FILE_NAME})
    set(CPP_FULL_NAME "${CMAKE_CURRENT_BINARY_DIR}/${CPP_FILE_NAME}")
    add_custom_command(
        OUTPUT "${CPP_FULL_NAME}"
        COMMAND ${CXXTESTGEN} --runner=ErrorPrinter --output "${CPP_FULL_NAME}"
         --template "${cxx_template}" "${source}"
        DEPENDS "${source}" "${cxx_template}"
    )
    add_executable(${target} ${CPP_FULL_NAME})
    set_target_properties(${target} PROPERTIES COMPILE_FLAGS "-Wno-effc++")
    target_link_libraries( ${target} ${CMAKE_THREAD_LIBS_INIT} )
    add_test( NAME ${target} COMMAND ${MPIEXEC} ${MPIEXEC_PREFLAGS}
        ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${CMAKE_CURRENT_BINARY_DIR}/${target} )
    set_tests_properties( ${target} PROPERTIES TIMEOUT ${test_timeout} ) # Detect deadlocks
endfunction( add_mpi_cxx_test )

function( add_serial_cxx_test target source )
    # Makes sure that MPI is initialized even when only one rank is required
    if( USE_MPI )
        add_mpi_cxx_test( ${target} ${source} 1 )
    else( USE_MPI )
        add_cxx_test( ${target} ${source} )
    endif( USE_MPI )
endfunction( add_serial_cxx_test )

# Add sub tests
add_subdirectory( as )

# Add tests
add_serial_cxx_test( TestConfig ${CMAKE_CURRENT_SOURCE_DIR}/Config.t.h )
add_serial_cxx_test( TestDispatcher ${CMAKE_CURRENT_SOURCE_DIR}/Dispatcher.t.h )
add_serial_cxx_test( TestModuleSync ${CMAKE_CURRENT_SOURCE_DIR}/Module.t.h )
set_tests_properties( TestModuleSync PROPERTIES ENVIRONMENT "ASYNC_MODE=SYNC" )
add_serial_cxx_test( TestModuleThread ${CMAKE_CURRENT_SOURCE_DIR}/Module.t.h )
set_tests_properties( TestModuleThread PROPERTIES ENVIRONMENT "ASYNC_MODE=THREAD" )
if( USE_MPI )
    add_mpi_cxx_test( TestModuleMPI ${CMAKE_CURRENT_SOURCE_DIR}/Module.t.h 3 )
    set_tests_properties( TestModuleMPI PROPERTIES ENVIRONMENT "ASYNC_MODE=MPI" )
    add_mpi_cxx_test( TestModuleMPIAsync ${CMAKE_CURRENT_SOURCE_DIR}/Module.t.h 3 )
    set_tests_properties( TestModuleMPIAsync PROPERTIES ENVIRONMENT "ASYNC_MODE=MPI;ASYNC_MPI_COPY=1" )
endif( USE_MPI )