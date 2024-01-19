
#ifndef SRC_ANALYSIS_LINEARSYSTEM_DIRECTSOLVER_H_
#define SRC_ANALYSIS_LINEARSYSTEM_DIRECTSOLVER_H_

#include "linearsystem.h"

#include "analysis/linearsystem/matrices/matrix_distributed.h"
#include "analysis/linearsystem/matrices/vector_distributed.h"

namespace espreso {

template <typename T>
struct DirectLinearSystemSolver: LinearSystemSolver<T> {

	DirectLinearSystemSolver()
	{
		LinearSystemSolver<T>::A = &A;
		LinearSystemSolver<T>::x = &x;
		LinearSystemSolver<T>::b = &b;
		LinearSystemSolver<T>::dirichlet = &dirichlet;
	}

	~DirectLinearSystemSolver()
	{

	}

protected:
	void setDirichlet();

	Matrix_Distributed<Matrix_CSR, T> A;
	Vector_Distributed<Vector_Dense, T> x, b;
	Vector_Distributed<Vector_Sparse, T> dirichlet;
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_DIRECTSOLVER_H_ */
