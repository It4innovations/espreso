
#ifndef SOLVER_SPECIFIC_DENSESOLVER_H_
#define SOLVER_SPECIFIC_DENSESOLVER_H_

// #include "esbasis.h"
// #include "../generic/utils.h"
#include "../generic/SparseMatrix.h"

namespace espreso {

class DenseSolver
{

public:
	virtual ~DenseSolver() {};

	virtual void ImportMatrix(SparseMatrix & A) = 0;
	virtual void ImportMatrix_fl(SparseMatrix & A) = 0;

	virtual void ImportMatrix_wo_Copy(SparseMatrix & A) = 0;

	virtual int Factorization(const std::string &str) = 0;
	virtual void Clear() = 0;
	virtual void SetThreaded() = 0;


	virtual void Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT rhs_start_index, MKL_INT sol_start_index) = 0;
	virtual void Solve( SEQ_VECTOR <double> & rhs, SEQ_VECTOR <double> & sol, MKL_INT n_rhs) = 0;
	virtual void Solve( SEQ_VECTOR <double> & rhs_sol) = 0;
};

}

#endif
