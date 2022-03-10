
#ifndef SRC_WRAPPERS_MATH_MKL_W_MKL_H_
#define SRC_WRAPPERS_MATH_MKL_W_MKL_H_

#include "wrappers/pardiso/w.pardiso.h"
#include "mkl_spblas.h"

namespace espreso {

struct Matrix_CSR_External_Representation: public PARDISOParameters
{
	sparse_matrix_t inspector;
};

struct Matrix_IJV_External_Representation
{
	sparse_matrix_t inspector;
};

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
