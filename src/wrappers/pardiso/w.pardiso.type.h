
#ifndef SRC_WRAPPERS_PARDISO_W_PARDISO_TYPE_H_
#define SRC_WRAPPERS_PARDISO_W_PARDISO_TYPE_H_

#include "math/math.h"

namespace espreso {

template<typename T>
int _pardisoType(const Matrix_CSR<T> &x)
{
	switch (x.type) {
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE:    return  2;
	case Matrix_Type::REAL_SYMMETRIC_INDEFINITE:           return -2;
	case Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC:         return  1;
	case Matrix_Type::REAL_NONSYMMETRIC:                   return 11;
	case Matrix_Type::COMPLEX_HERMITIAN_POSITIVE_DEFINITE: return  4;
	case Matrix_Type::COMPLEX_HERMITIAN_INDEFINITE:        return -4;
	case Matrix_Type::COMPLEX_SYMMETRIC:                   return  6;
	case Matrix_Type::COMPLEX_STRUCTURALLY_SYMMETRIC:      return  3;
	case Matrix_Type::COMPLEX_NONSYMMETRIC:                return 13;
	}
	return 0;
}

}



#endif /* SRC_WRAPPERS_PARDISO_W_PARDISO_TYPE_H_ */
