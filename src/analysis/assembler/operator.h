
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATOR_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATOR_H_

#include <cstdio>

namespace espreso {

class Assembler;

struct Operator {
	static const int print = 2;

	const int interval, isconst, update;

	Operator(int interval, bool isconst, bool update): interval(interval), isconst(isconst), update(update) {}
	virtual ~Operator() {}

	virtual void operator++() =0;
};

struct ActionOperator: public Operator
{
	using Operator::Operator;

	virtual ~ActionOperator() {}

	virtual void operator()() =0;
	virtual void reset() =0;
};


struct OperatorBuilder {
	const char *name;

	OperatorBuilder(const char *name): name(name) {}
	virtual ~OperatorBuilder() {}

	virtual void now() =0;
};

struct ElementOperatorBuilder: public OperatorBuilder {
	ElementOperatorBuilder(const char *name): OperatorBuilder(name) {}
	virtual ~ElementOperatorBuilder() {}

	void now();

	virtual void apply(int interval) =0;
};

struct BoundaryOperatorBuilder: public OperatorBuilder {
	BoundaryOperatorBuilder(const char *name): OperatorBuilder(name) {}
	virtual ~BoundaryOperatorBuilder() {}

	void now();

	virtual void apply(int region, int interval) =0;
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATOR_H_ */
