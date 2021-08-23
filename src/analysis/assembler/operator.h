
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATOR_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATOR_H_

namespace espreso {

class Assembler;

struct ActionOperator {
	static const int print = 2;

	int isconst, update;

	ActionOperator(): isconst(0), update(1) {}
	virtual ~ActionOperator() {}

	virtual void operator++() =0;
	virtual void operator()() =0;

	virtual void move(int n) =0;
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATOR_H_ */
