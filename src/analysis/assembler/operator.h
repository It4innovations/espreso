
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATOR_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATOR_H_

#include "esinfo/eslog.h"
#include "math/simd/simd.h"

namespace espreso {

class Assembler;

struct ActionOperator {
	int isconst, update;

	ActionOperator(): isconst(1), update(1) {}
	virtual ~ActionOperator() {}

	virtual void operator++() =0;
	virtual void operator()() =0;
	virtual void simd() { eslog::error("not implemented simd operator\n"); }
	virtual void peel(size_t size) { simd(); }

	virtual void move(int n)  =0;

};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATOR_H_ */
