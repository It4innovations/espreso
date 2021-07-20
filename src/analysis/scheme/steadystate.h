
#ifndef SRC_ANALYSIS_SCHEME_STEADYSTATE_H_
#define SRC_ANALYSIS_SCHEME_STEADYSTATE_H_

#include "scheme.h"

#include "analysis/linearsystem/linearsystem.h"
#include "math2/generalization/matrix_base.h"

namespace espreso {

struct AX_SteadyState: public AX_Scheme {

	void init(AX_LinearSystem<double> *system, AX_HeatTransfer &assembler);

	void reassemble(AX_HeatTransfer &assembler, bool &updatedA, bool &updatedB);
	void solved();

	Matrix_Base<double> *K;
	Vector_Base<double> *f;

	AX_LinearSystem<double> *system;
};

}

#endif /* SRC_ANALYSIS_SCHEME_STEADYSTATE_H_ */
