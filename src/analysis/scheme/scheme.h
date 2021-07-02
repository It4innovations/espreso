
#ifndef SRC_ANALYSIS_SCHEME_SCHEME_H_
#define SRC_ANALYSIS_SCHEME_SCHEME_H_

#include "analysis/assembler/module/heattransfer.h"

namespace espreso {

struct AX_Scheme {
	virtual ~AX_Scheme() {}

	virtual void reassemble(AX_HeatTransfer &assembler, bool &A, bool &b) {}
};

}



#endif /* SRC_ANALYSIS_SCHEME_SCHEME_H_ */
