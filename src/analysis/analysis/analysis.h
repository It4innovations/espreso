
#ifndef SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_
#define SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_

#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

namespace step { struct Step; }

struct Analysis {
	virtual ~Analysis() {}

	virtual void analyze() =0;
	virtual void run(step::Step &step) =0;
};

}

#endif /* SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_ */
