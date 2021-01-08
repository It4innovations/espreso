
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_

#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"

namespace espreso {

struct KFiller: public Operator {
	KFiller(const ParameterData &stiffness, Kernel::InstanceFiller &filler, PerElementSize size, int interval)
	: Operator(interval, false, true),
	  stiffness(stiffness, interval, stiffness.size),
	  filler(filler),
	  size(stiffness.increment(size, interval)) {}

	InputParameterIterator stiffness;
	Kernel::InstanceFiller &filler;
	int size;

	void operator++()
	{
		// increment by operator()
	}
};

struct KSymmFiller: public KFiller {
	GET_NAME(KSymmFiller)
	using KFiller::KFiller;

	void operator()()
	{
		for (int r = 0; r < size; ++r) {
			for (int c = 0; c < size; ++c, ++stiffness.data) {
				if (r <= c) {
					filler.K[*filler.offset++] += *stiffness.data;
				}
			}
		}
	}
};

struct KFullFiller: public KFiller {
	GET_NAME(KFullFiller)
	using KFiller::KFiller;

	void operator()()
	{
		for (int r = 0; r < size; ++r) {
			for (int c = 0; c < size; ++c, ++stiffness.data) {
				filler.K[*filler.offset++] += *stiffness.data;
			}
		}
	}
};

struct MatricesFiller: public ElementOperatorBuilder {
	GET_NAME(MatricesFiller)

	HeatTransferKernelOpt &kernel;
	Kernel::InstanceFiller &filler;

	MatricesFiller(HeatTransferKernelOpt &kernel, Kernel::InstanceFiller &filler): kernel(kernel), filler(filler)
	{

	}

	void apply(int interval)
	{
		if (kernel.translationMotions.gp.builder) {
			iterate_elements(KFullFiller(kernel.linearSystem.stiffness, filler, enodes, interval));
		} else {
			iterate_elements(KSymmFiller(kernel.linearSystem.stiffness, filler, enodes, interval));
		}
	}
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_ */
