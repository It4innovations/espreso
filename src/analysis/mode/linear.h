
#ifndef SRC_ANALYSIS_MODE_LINEAR_H_
#define SRC_ANALYSIS_MODE_LINEAR_H_

#include "analysis/assembler/module/heattransfer.h"
#include "analysis/linearsystem/linearsystem.h"
#include "analysis/scheme/scheme.h"

namespace espreso {

struct AX_Linear {

	void init(AX_LinearSystem<double> *system);

	bool solve(AX_Scheme &scheme, AX_HeatTransfer &assembler);

private:
	AX_LinearSystem<double> *system;
};

}

#endif /* SRC_ANALYSIS_MODE_LINEAR_H_ */
