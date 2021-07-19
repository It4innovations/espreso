
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_

#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"

namespace espreso {

struct KFiller: public Operator {
	KFiller(const ParameterData &stiffness, double *K, const esint *perm, PerElementSize size, int interval)
	: Operator(interval, false, true),
	  stiffness(stiffness, interval),
	  K(K), perm(perm),
	  size(stiffness.increment(size, interval)) {}

	InputParameterIterator stiffness;
	double *K;
	const esint *perm;
	int size;

	void operator++()
	{
		// increment by operator()
	}
};

struct KSymmFiller: public KFiller {
	using KFiller::KFiller;

	void operator()()
	{
		for (int r = 0; r < size; ++r) {
			for (int c = 0; c < size; ++c, ++stiffness.data) {
				if (r <= c) {
					K[*perm++] += *stiffness.data;
				}
			}
		}
	}
};

struct KFullFiller: public KFiller {
	using KFiller::KFiller;

	void operator()()
	{
		for (int r = 0; r < size; ++r) {
			for (int c = 0; c < size; ++c, ++stiffness.data) {
				K[*perm++] += *stiffness.data;
			}
		}
	}
};

struct MatricesFiller: public ElementOperatorBuilder {
	bool issymmetric;
	double *K;
	const esint *perm;

	MatricesFiller(double *K, const esint *perm): ElementOperatorBuilder("FILL K"), issymmetric(false), K(K), perm(perm)
	{

	}

	void apply(int interval)
	{
//		if (issymmetric) {
//			iterate_elements(KSymmFiller(assembler.elements.stiffness, K, perm, enodes, interval));
//		} else {
//			iterate_elements(KFullFiller(assembler.elements.stiffness, K, perm, enodes, interval));
//		}
	}
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_ */
