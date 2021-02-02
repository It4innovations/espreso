
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_

#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"

namespace espreso {

struct KFiller: public Operator {
	KFiller(const ParameterData &stiffness, double *K, esint *perm, PerElementSize size, int interval)
	: Operator(interval, false, true),
	  stiffness(stiffness, interval, stiffness.size),
	  K(K), perm(perm),
	  size(stiffness.increment(size, interval)) {}

	InputParameterIterator stiffness;
	double *K;
	esint *perm;
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
					K[*perm++] += *stiffness.data;
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
				K[*perm++] += *stiffness.data;
			}
		}
	}
};

struct MatricesFiller: public ElementOperatorBuilder {
	GET_NAME(MatricesFiller)

	HeatTransferModuleOpt &kernel;
	double *K;
	esint *perm;

	MatricesFiller(HeatTransferModuleOpt &kernel, double *K, esint *perm): kernel(kernel), K(K), perm(perm)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		return true;
	}

	void apply(int interval)
	{
		if (kernel.translationMotions.gp.isset) {
			iterate_elements(KFullFiller(kernel.linearSystem.stiffness, K, perm, enodes, interval));
		} else {
			iterate_elements(KSymmFiller(kernel.linearSystem.stiffness, K, perm, enodes, interval));
		}
	}
};

struct VectorFiller: public Operator {
	GET_NAME(RHSFiller)
	VectorFiller(const ParameterData &rhs, double *RHS, esint *perm, PerElementSize size, int interval)
	: Operator(interval, false, true),
	  rhs(rhs, interval, enodes),
	  RHS(RHS), perm(perm),
	  size(rhs.increment(size, interval)) {}

	InputParameterIterator rhs;
	double *RHS;
	esint *perm;
	int size;

	void operator++()
	{
		// increment by operator()
	}

	void operator()()
	{
		for (int r = 0; r < size; ++r) {
			RHS[*perm++] += *rhs.data++;
		}
	}
};

struct RHSFiller: public BoundaryOperatorBuilder {
	GET_NAME(MatricesFiller)

	HeatTransferModuleOpt &kernel;
	double *RHS;
	esint *perm;

	RHSFiller(HeatTransferModuleOpt &kernel, double *RHS, esint *perm): kernel(kernel), RHS(RHS), perm(perm)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		return true;
	}

	void apply(int region, int interval)
	{
		iterate_boundary(VectorFiller(kernel.linearSystem.boundary.rhs.regions[region], RHS, perm, enodes, interval), region);
	}
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_ */
