/*! \file CopyMatrix.hpp
    \brief task mangemanet of dissection algorithm
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Aug. 09th 2015
    \date   Feb. 29th 2016
*/

// This file is part of Dissection
// 
// Dissection is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Dissection is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Dissection.  If not, see <http://www.gnu.org/licenses/>.

#ifndef _COPY_MATRIX_
#define _COPY_MATRIX_
#include <cassert>
#include <vector>

#include "Driver/DissectionMatrix.hpp"

#include <vector>

using std::vector;

void CopySparseMatrix(SparseMatrix<double> *b,
		      SparseMatrix<quadruple, double, double> *a);

void CopySparseMatrix(SparseMatrix<complex<double>, complex<double>, double> *b,
		      SparseMatrix<complex<quadruple>, complex<double>, double> *a);


template<typename T, typename W>
void CopySquareBlockMatrix(SquareBlockMatrix<W> &b,
			   SquareBlockMatrix<T> &a);

template<typename T, typename W>
void CopyRectBlockMatrix(RectBlockMatrix<W> &b,
			 RectBlockMatrix<T> &a);

template<typename T, typename U, typename W, typename Z>
void CopyTridiagBlockMatrix(TridiagBlockMatrix<W, Z> &b,
			    TridiagBlockMatrix<T, U> &a,
			    W *coef);

template<typename T, typename U, typename W, typename Z>
void CopyDissectionMatrix(DissectionMatrix<W, Z>* a,
			  DissectionMatrix<T, U>* b,
			  SquareBlockMatrix<W> *diag,           // pointers
			  RectBlockMatrix<W> *lower,
			  RectBlockMatrix<W> *upper);

template<typename T, typename W>
void CopySchurMatrix(SchurMatrix<W> &b,
		     SchurMatrix<T> &a);

template<typename T, typename W>
void CopyKernelMatrix(KernelMatrix<W> &b,
		      KernelMatrix<T> &a);

#endif
