
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_

#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"

namespace espreso {

struct KFiller: public Operator {
	KFiller(const ParameterData &stiffness, double *K, esint *perm, PerElementSize size, int interval)
	: Operator(interval, false, true),
	  stiffness(stiffness, interval),
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
	HeatTransferModuleOpt &kernel;
	bool issymmetric;
	double *K;
	esint *perm;

	MatricesFiller(HeatTransferModuleOpt &kernel, double *K, esint *perm): ElementOperatorBuilder("FILL K"), kernel(kernel), issymmetric(true), K(K), perm(perm)
	{
		// TODO: allow various symmetricity for each interval
		for (size_t i = 0; i < kernel.translationMotions.gp.isset.size(); ++i) {
			issymmetric &= !kernel.translationMotions.gp.isset[i];
		}
		for (size_t i = 0; i < kernel.material.model.anisotropic.isset.size(); ++i) {
			issymmetric &= !kernel.material.model.anisotropic.isset[i];
		}
	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		return true;
	}

	void apply(int interval)
	{
		if (issymmetric) {
			iterate_elements(KSymmFiller(kernel.elements.stiffness, K, perm, enodes, interval));
		} else {
			iterate_elements(KFullFiller(kernel.elements.stiffness, K, perm, enodes, interval));
		}
	}
};

struct VectorFiller: public Operator {
	VectorFiller(const ParameterData &rhs, double *RHS, esint *perm, PerElementSize size, int interval)
	: Operator(interval, false, true),
	  rhs(rhs, interval),
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
	HeatTransferModuleOpt &kernel;
	double *RHS;
	esint *perm;

	RHSFiller(HeatTransferModuleOpt &kernel, double *RHS, esint *perm): BoundaryOperatorBuilder("FILL RHS"), kernel(kernel), RHS(RHS), perm(perm)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		return true;
	}

	void apply(int region, int interval)
	{
		iterate_boundary(VectorFiller(kernel.elements.boundary.rhs.regions[region], RHS, perm, enodes, interval), region);
	}
};

struct KFillerSimd: public Operator {
	KFillerSimd(const ParameterData &stiffness, double *K, esint *perm, PerElementSize size, int interval)
	: Operator(interval, false, true),
	  stiffness(stiffness, interval),
	  K(K), perm(perm),
	  size(stiffness.increment(size, interval)),
	  elementInPack(0) {}

	InputParameterIterator stiffness;
	double *K;
	esint *perm;
	int size;
	int elementInPack;

	void operator++()
	{
		// Loop iterate_elements_simd cant be used -- overfill of K when elements not multiple of SIMD lenght
		// Compicated operator++ will move iterator by SIMD elements when processed
		// indexing within SIMD pack done manually. (see derived operators). Workaround for now
		elementInPack++;
		if(elementInPack == SIMD::size)
		{
			stiffness += SIMD::size;
			elementInPack = 0;	
		}
	}
};

struct KSymmFillerSimd: public KFillerSimd {
	using KFillerSimd::KFillerSimd;

	void operator()()
	{
		for (int r = 0; r < size; ++r) {
			for (int c = 0; c < size; ++c) {
				if (r <= c) {
					int item = r * size + c;
					K[*perm++] += stiffness.data[item * SIMD::size + elementInPack];
				}
			}
		}
	}
};

struct KFullFillerSimd: public KFillerSimd {
	using KFillerSimd::KFillerSimd;

	void operator()()
	{
		for (int r = 0; r < size; ++r) {
			for (int c = 0; c < size; ++c) {
				int item = r * size + c;
				K[*perm++] += stiffness.data[item * SIMD::size + elementInPack];
			}
		}
	}
};

struct MatricesFillerSimd: public ElementOperatorBuilder {
	HeatTransferModuleOpt &kernel;
	bool issymmetric;
	double *K;
	esint *perm;

	MatricesFillerSimd(HeatTransferModuleOpt &kernel, double *K, esint *perm): ElementOperatorBuilder("FILL K"), kernel(kernel), issymmetric(true), K(K), perm(perm)
	{
		// TODO: allow various symmetricity for each interval
		for (size_t i = 0; i < kernel.translationMotions.gp.isset.size(); ++i) {
			issymmetric &= !kernel.translationMotions.gp.isset[i];
		}
		for (size_t i = 0; i < kernel.material.model.anisotropic.isset.size(); ++i) {
			issymmetric &= !kernel.material.model.anisotropic.isset[i];
		}
	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		return true;
	}

	void apply(int interval)
	{
		if (issymmetric) {
			iterate_elements(KSymmFillerSimd(kernel.elementsSimd.stiffness, K, perm, enodes, interval));
		} else {
			iterate_elements(KFullFillerSimd(kernel.elementsSimd.stiffness, K, perm, enodes, interval));
		}
	}
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_ */
