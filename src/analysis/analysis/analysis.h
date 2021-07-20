
#ifndef SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_
#define SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_

#include "analysis/linearsystem/linearsystem.h"

namespace espreso {

struct Analysis {
	virtual ~Analysis() {}

	virtual void init() =0;
	virtual void run() =0;

};
}

#endif /* SRC_ANALYSIS_ANALYSIS_ANALYSIS_H_ */
