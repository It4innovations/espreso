
#ifndef SRC_FETI_FETI_H_
#define SRC_FETI_FETI_H_

#include "config/ecf/linearsolver/feti.h"
#include "esinfo/stepinfo.h"
#include "math/primitives/vector_dense.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/matrix_csr.h"
#include "math/physics/matrix_feti.h"
#include "math/physics/vector_feti.h"

#include <map>
#include <set>

namespace espreso {

template <typename T> class IterativeSolver;
template <typename T> struct Preconditioner;
template <typename T> struct Projector;
template <typename T> class DualOperator;

template<typename T>
struct FETI {
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
			std::vector<int> D2C;

			Matrix_CSR<T, int> B1;
			Vector_Dense<T, int> duplication;
		};

		std::vector<esint> cmap; // size, ndomains <d0, d1, ..., dn>; size, ndomains <>; ...;

		esint dirichlet, nhalo, size;
		std::vector<Domain> domain;
		Vector_Dense<T, int> c;
	};

	FETI(FETIConfiguration &configuration);
	~FETI();

	void info() const;

	bool set(const step::Step &step);
	bool update(const step::Step &step);
	bool solve(const step::Step &step);

	FETIConfiguration &configuration;
	SystemInfo sinfo;

	DOFsDecomposition *decomposition;
	Matrix_FETI<Matrix_CSR, T> K;
	Vector_FETI<Vector_Dense, T> f, x;
	Regularization regularization;
	EqualityConstraints equalityConstraints;

	IterativeSolver<T> *iterativeSolver = nullptr;
	Preconditioner<T> *preconditioner = nullptr;
	Projector<T> *projector = nullptr;
	DualOperator<T> *dualOperator = nullptr;
};

}

#endif /* SRC_FETI_FETI_H_ */
