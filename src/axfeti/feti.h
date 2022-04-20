
#ifndef SRC_AXFETI_FETI_H_
#define SRC_AXFETI_FETI_H_

#include "config/ecf/linearsolver/feti.h"
#include "esinfo/stepinfo.h"
#include "math/primitives/vector_dense.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"
#include "math/physics/matrix_feti.h"
#include "math/physics/vector_feti.h"
#include "math/feti/lmap.h"

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
		esint R1totalSize, R2totalSize;
		esint lambdasLocal, lambdasTotal;
	};

	struct Regularization {
		Matrix_FETI<Matrix_Dense, T> R1, R2;
		Matrix_FETI<Matrix_CSR, T> RegMat;
	};

	struct EqualityConstraints {
		struct Domain {
			esint nhalo;
			std::vector<esint> D2C;

			Matrix_CSR<T> B1;
			Vector_Dense<T> duplication;
		};

		esint global, nhalo, paired, local, nn;
		std::vector<LMAP> lmap;
		std::vector<esint> ordered;
		std::vector<Domain> domain;
		Vector_Dense<T> c;
	};

	AX_FETI(FETIConfiguration &configuration): configuration(configuration) {}

	void info() const;

	bool set(const step::Step &step, Matrix_FETI<Matrix_CSR, T> &K, const Regularization &regularization, const EqualityConstraints &equalityConstraints);
	bool update(const step::Step &step, Matrix_FETI<Matrix_CSR, T> &K, Vector_FETI<Vector_Dense, T> &f);
	bool solve(const step::Step &step, Vector_FETI<Vector_Dense, T> &x);

	FETIConfiguration &configuration;
	SystemInfo sinfo;

	const step::Step *step = nullptr;
	const Matrix_FETI<Matrix_CSR, T> *K = nullptr;
	const Regularization *regularization = nullptr;
	const EqualityConstraints *equalityConstraints = nullptr;

	Vector_FETI<Vector_Dense, T> *f = nullptr;
	Vector_FETI<Vector_Dense, T> *x = nullptr;

	IterativeSolver<T> *iterativeSolver = nullptr;
	Preconditioner<T> *preconditioner = nullptr;
	Projector<T> *projector = nullptr;
	DualOperator<T> *dualOperator = nullptr;
};

}

#endif /* SRC_AXFETI_FETI_H_ */
