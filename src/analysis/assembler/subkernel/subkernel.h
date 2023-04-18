
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_SUBKERNEL_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_SUBKERNEL_H_

#include "analysis/assembler/module/assembler.h"
#include "math/simd/simd.h"

namespace espreso {

struct SubKernel {
	int isconst, isactive;
	Assembler::Action action;

	SubKernel(): isconst(1), isactive(0), action(Assembler::VOID) {}
	virtual ~SubKernel() {}

	void setActiveness(Assembler::Action action)
	{
		isactive = isactive && (this->action & action);
	}

	void setActiveness(Assembler::Action action, int guard)
	{
		isactive = guard && isactive && (this->action & action);
	}

	virtual const char* name() const =0;
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_SUBKERNEL_H_ */
