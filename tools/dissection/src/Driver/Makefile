#! \file   Makefile
#  \brief  source files in Driver/ directory and depending libraries
#  \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
#  \date   Jul. 12th 2015
#  \date   Feb. 29th 2016

# This file is part of Dissection
# 
# Dissection is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Dissection is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Dissection.  If not, see <http://www.gnu.org/licenses/>.

DIRECTORY  = Driver
SOURCESCPP = C_threads_tasks.cpp \
             C_Dsub.cpp \
             C_DFullLDLt.cpp \
             DissectionMatrix.cpp \
             DissectionQueue.cpp \
             DissectionSolver.cpp \
             QueueRuntime.cpp \
             C_BlasRoutines.cpp \
             TridiagBlockMatrix.cpp \
             TridiagQueue.cpp \
             C_KernDetect.cpp \
             CopyMatrix.cpp
# PardisoSolver.cpp 

LIBDEP     = -L../../lib -lAlgebra -lSplitters -lUtilities -lCompiler -lBlock_tridiag_matrix_s -lModules_s -lFeti_Abort_s $(SHARED_BLAS) $(SHARED_INTEL)
include ../Defs.inc
include ../Make.inc
