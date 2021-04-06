
#ifndef SRC_WRAPPERS_MATH_MKL_W_MKL_H_
#define SRC_WRAPPERS_MATH_MKL_W_MKL_H_

#include "wrappers/pardiso/w.pardiso.h"
#include "mkl_spblas.h"

namespace espreso {
namespace MATH {

struct CSRHandlerData: public PARDISOParameters {
	sparse_matrix_t inspector;

	CSRHandlerData(): inspector{}
	{

	}

	void info(esint &nrows, esint &ncols, esint &nnz);
	void info(esint &nrows, esint &ncols, esint &nnz, esint* &rows, esint* &cols, double* &vals);
};

struct IJVHandlerData {
	sparse_matrix_t inspector;

	IJVHandlerData(): inspector{} {};
};

}
}



#endif /* SRC_WRAPPERS_MATH_MKL_W_MKL_H_ */
