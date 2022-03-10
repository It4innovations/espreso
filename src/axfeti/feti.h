
#ifndef SRC_AXFETI_FETI_H_
#define SRC_AXFETI_FETI_H_

#include "config/ecf/linearsolver/feti.h"
#include "math2/generalization/matrix_distributed.h"
#include "math2/generalization/matrix_feti.h"
#include "math2/generalization/vector_feti.h"

namespace espreso {

template <typename T> class IterativeSolver;
template <typename T> class Preconditioner;
template <typename T> class Projector;
template <typename T> class DualOperator;

template<typename T>
class AX_FETI {

public:
	struct SystemInfo {
		esint domains, clusters;
		esint R1offset, R2offset;
		esint R1size, R2size;
	};

	struct Regularization {
		Matrix_FETI<Matrix_Dense, T> R1, R2;
		Matrix_FETI<Matrix_CSR, T> RegMat;
	};

	struct EqualityConstraints {
		std::vector<std::vector<esint> > D2C; // domains to cluster
		std::vector<esint> C2G; // cluster to global
		serializededata<esint, int> *L2MPI = nullptr; // lambda to MPI processes
		Matrix_FETI<Matrix_IJV, T> B1;
		Vector_FETI<Vector_Dense, T> B1c, B1Duplication;

		~EqualityConstraints()
		{
			if (L2MPI) delete L2MPI;
		}
	};

	AX_FETI(FETIConfiguration &configuration): configuration(configuration) {}

	void info() const;

	bool set(const step::Step &step, Matrix_FETI<Matrix_CSR, T> &K, const Regularization &regularization, const EqualityConstraints &equalityConstraints);
	bool update(const step::Step &step, Matrix_FETI<Matrix_CSR, T> &K);
	bool solve(const step::Step &step, const Vector_FETI<Vector_Dense, T> &f, Vector_FETI<Vector_Dense, T> &x);

	FETIConfiguration &configuration;
	SystemInfo sinfo;

	const step::Step *step = nullptr;
	const Matrix_FETI<Matrix_CSR, T> *K = nullptr;
	const Regularization *regularization = nullptr;
	const EqualityConstraints *equalityConstraints = nullptr;

	IterativeSolver<T> *iterativeSolver = nullptr;
	Preconditioner<T> *preconditioner = nullptr;
	Projector<T> *projector = nullptr;
	DualOperator<T> *dualOperator = nullptr;
};

}

#endif /* SRC_AXFETI_FETI_H_ */
