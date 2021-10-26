
#ifndef SRC_AXFETI_FETI_H_
#define SRC_AXFETI_FETI_H_

#include "iterativesolver/iterativesolver.h"
#include "preconditioner/preconditioner.h"
#include "projector/projector.h"

#include "config/ecf/linearsolver/feti.h"
#include "math2/generalization/matrix_distributed.h"
#include "math2/generalization/matrix_feti.h"
#include "math2/generalization/vector_feti.h"

namespace espreso {

template<typename T>
class AX_FETI {

public:
	struct Regularization {
		Matrix_FETI<Matrix_Dense, T> R1, R2;
		Matrix_FETI<Matrix_CSR, T> RegMat;
	};

	struct EqualityConstraints {
		Matrix_FETI<Matrix_IJV, T> B1Dirichlet, B1Gluing;
		Vector_FETI<Vector_Dense, T> B1c, B1Duplication;
	};

	AX_FETI(FETIConfiguration &configuration): configuration(configuration) {}

	bool set(const Matrix_FETI<Matrix_CSR, T> &K, const Regularization &regularization, const EqualityConstraints &equalityConstraints);

	bool update(const Matrix_Distributed<Matrix_CSR, T> &K);
	bool solve(const Vector_Distributed<Vector_Dense, T> &f, Vector_Distributed<Vector_Dense, T> &x);

	FETIConfiguration &configuration;

	Matrix_FETI<Matrix_CSR, T> Kplus;
	Matrix_Distributed<Matrix_CSR, T> G;

	IterativeSolver<T> *iterativeSolver = nullptr;
	Preconditioner<T> *preconditioner = nullptr;
	Projector<T> *projector = nullptr;
};

}

#endif /* SRC_AXFETI_FETI_H_ */
