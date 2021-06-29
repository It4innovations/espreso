
#ifndef SRC_ANALYSIS_SCHEME_HARMONIC_COMPLEX_H_
#define SRC_ANALYSIS_SCHEME_HARMONIC_COMPLEX_H_

#include "analysis/linearsolver/linearsolver.h"

#include "math2/generalization/vector_base.h"
#include "math2/generalization/matrix_base.h"

namespace espreso {

struct AX_HarmonicComplex {

	void init(DirectSolver<double> *solver);
	void init(FETISolver<double> *solver);
	void init(MultigridSolver<double> *solver);

	Matrix_Base<std::complex<double> > *K, *M;
	Vector_Base<std::complex<double> > *f;
};

}

#endif /* SRC_ANALYSIS_SCHEME_HARMONIC_COMPLEX_H_ */
