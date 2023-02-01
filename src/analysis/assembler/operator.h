
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATOR_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATOR_H_

#include "esinfo/eslog.h"
#include "math/simd/simd.h"

namespace espreso {

struct ActionOperator {
	enum Action {
		VOID     = 1 << 0,
		ASSEMBLE = 1 << 1,
		FILL     = 1 << 2,
		SOLUTION = 1 << 3
	};

	int isconst, update;
	Action action;

	ActionOperator(): isconst(1), update(1), action(Action::ASSEMBLE) {}
	virtual ~ActionOperator() {}

	virtual void operator++() {};
	virtual void operator()() {};

	virtual void move(int n) {};
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATOR_H_ */
