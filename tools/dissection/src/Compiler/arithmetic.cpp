/*! \file   arithmetic.cpp
    \brief  higher precision arithmetic for Kernel Detection
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Jul. 17th 2015
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

#include <cstdio>
#include "Compiler/arithmetic.hpp"
#include "Compiler/blas.hpp"

template<>
void printscalar<double>(FILE *fp, double x)
{
  fprintf(fp, "%16.8e ", x);
}

template<>
void printscalar<quadruple>(FILE *fp, quadruple x)
{
  fprintf(fp, "%16.8e ", quad2double(x));
}
template<>
void printscalar<complex<double> >(FILE *fp, complex<double> x)
{
  fprintf(fp, "(%16.8e %16.8e) ", x.real(), x.imag());
}

template<>
void printscalar<complex<quadruple> >(FILE *fp, complex<quadruple> x)
{
  fprintf(fp, "(%16.8e %16.8e) ",
	  quad2double(x.real()), quad2double(x.imag()));
}

template<typename T>
void printscalar(FILE *fp, T x)
{
  fprintf(stderr, "%s %d : printscalar is not implented\n", __FILE__, __LINE__);
}
#ifndef NO_OCTRUPLE
template
void printscalar<octruple>(FILE *fp, octruple x);
template
void printscalar<complex<octruple> >(FILE *fp, complex<octruple> x);
#endif

