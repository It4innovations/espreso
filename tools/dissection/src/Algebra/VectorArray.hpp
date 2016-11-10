/*! \file   ColumnMatrix.hpp
    \brief  Rectangular matrix view as a set of column vectors
    \author Xavier Juvigny, ONERA
    \date   Jan. 19th 2005
    \modification allocation of array by STL vector class
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Jun. 11th 2013
    \date   Jul. 12th 2015
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

#ifndef _ALGEBRA_VECTORARRAY_HPP
# define _ALGEBRA_VECTORARRAY_HPP

# include "Algebra/PlainMatrix.hpp"

template<typename T>
class VectorArray : public PlainMatrix<T>
{
public:

  using PlainMatrix<T>::coefs;
  using PlainMatrix<T>::addrCoefs;
  using PlainMatrix<T>::addrCoefs_pt;
  using PlainMatrix<T>::size;
  using PlainMatrix<T>::ZeroClear;
  using PlainMatrix<T>::init;
  using PlainMatrix<T>::free;

  VectorArray() : PlainMatrix<T>() {}
  VectorArray(int n) : PlainMatrix<T>() {
    PlainMatrix<T>::init(n);
  }
 
  ~VectorArray() { }
  
  virtual T& operator () (int i, int j)
  {
    return coefs()[i];
  }

  virtual const T& operator () (int i, int j) const
  {
    return coefs()[i];
  }
  virtual VectorArray<T>* clone() const
  {
    VectorArray<T> *ret = new VectorArray<T>;
    ret->copy(*this);
    return(ret);
  }

};

#endif
