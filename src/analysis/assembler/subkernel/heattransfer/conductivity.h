
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_CONDUCTIVITY_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_CONDUCTIVITY_H_

#include "subkernels.h"

namespace espreso {

struct Conductivity: SubKernel {
	const char* name() const { return "ConductivityKernel"; }

	Conductivity()
	: conductivity(nullptr), direct(true)
	{
		action = Assembler::ASSEMBLE | Assembler::REASSEMBLE;
	}

	void activate(const ThermalConductivityConfiguration *conductivity, bool direct)
	{
		this->conductivity = conductivity;
		this->direct = direct;
		this->isactive = 1;
	}

	const ThermalConductivityConfiguration *conductivity;
	bool direct;
};

template <size_t gps, size_t ndim, size_t etype, class Physics> struct ConductivityKernel;

template <size_t gps, size_t ndim, class Physics>
struct ConductivityKernel<gps, ndim, HeatTransferElementType::SYMMETRIC_ISOTROPIC, Physics>: Conductivity, Physics {
	ConductivityKernel(const Conductivity &base)
	: Conductivity(base),
	  kxx(this->conductivity->values.get(0, 0).evaluator)
	{
		isconst = kxx->isConst();
	}

	Evaluator *kxx;

	void simd(typename Physics::Element &element)
	{
		if (direct) {
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.conductivity[gp][s] = kxx->evaluate();
				}
			}
		} else {
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][s] = kxx->evaluate();
				}
			}
		}
	}
};

template <size_t gps, size_t ndim, class Physics>
struct ConductivityKernel<gps, ndim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC, Physics>: Conductivity, Physics {
	ConductivityKernel(const Conductivity &base)
	: Conductivity(base),
	  kxx(this->conductivity->values.get(0, 0).evaluator)
	{
		isconst = kxx->isConst();
	}

	Evaluator *kxx;

	void simd(typename Physics::Element &element)
	{
		if (direct) {
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.conductivity[gp][s] = kxx->evaluate();
				}
			}
		} else {
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][s] = kxx->evaluate();
				}
			}
		}
	}
};

template <size_t gps, class Physics>
struct ConductivityKernel<gps, 2, HeatTransferElementType::SYMMETRIC_GENERAL, Physics>: Conductivity, Physics {
	ConductivityKernel(const Conductivity &base)
	: Conductivity(base),
	  kxx(this->conductivity->values.get(0, 0).evaluator),
	  kxy(this->conductivity->values.get(0, 1).evaluator),
	  kyy(this->conductivity->values.get(1, 1).evaluator)
	{
		isconst = kxx->isConst() && kxy->isConst() && kyy->isConst();
	}

	Evaluator *kxx, *kxy, *kyy;

	void simd(typename Physics::Element &element)
	{
		if (direct) {
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.conductivity[gp][0][s] = kxx->evaluate();
					element.conductivity[gp][1][s] = kxy->evaluate();
					element.conductivity[gp][2][s] = kyy->evaluate();
				}
			}
		} else {
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][0][s] = kxx->evaluate();
					element.ecf.conductivity[gp][1][s] = kxy->evaluate();
					element.ecf.conductivity[gp][2][s] = kyy->evaluate();
				}
			}
		}
	}
};

template <size_t gps, class Physics>
struct ConductivityKernel<gps, 3, HeatTransferElementType::SYMMETRIC_GENERAL, Physics>: Conductivity, Physics {
	ConductivityKernel(const Conductivity &base)
	: Conductivity(base),
	  kxx(this->conductivity->values.get(0, 0).evaluator),
	  kxy(this->conductivity->values.get(0, 1).evaluator),
	  kxz(this->conductivity->values.get(0, 2).evaluator),
	  kyy(this->conductivity->values.get(1, 1).evaluator),
	  kyz(this->conductivity->values.get(1, 2).evaluator),
	  kzz(this->conductivity->values.get(2, 2).evaluator)
	{
		isconst =
				kxx->isConst() && kxy->isConst() && kxz->isConst() &&
				kyy->isConst() && kyz->isConst() && kzz->isConst();
	}

	Evaluator *kxx, *kxy, *kxz, *kyy, *kyz, *kzz;

	void simd(typename Physics::Element &element)
	{
		if (direct) {
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.conductivity[gp][0][s] = kxx->evaluate();
					element.conductivity[gp][1][s] = kxy->evaluate();
					element.conductivity[gp][2][s] = kxz->evaluate();
					element.conductivity[gp][3][s] = kyy->evaluate();
					element.conductivity[gp][4][s] = kyz->evaluate();
					element.conductivity[gp][5][s] = kzz->evaluate();
				}
			}
		} else {
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][0][s] = kxx->evaluate();
					element.ecf.conductivity[gp][1][s] = kxy->evaluate();
					element.ecf.conductivity[gp][2][s] = kxz->evaluate();
					element.ecf.conductivity[gp][3][s] = kyy->evaluate();
					element.ecf.conductivity[gp][4][s] = kyz->evaluate();
					element.ecf.conductivity[gp][5][s] = kzz->evaluate();
				}
			}
		}
	}
};

template <size_t gps, class Physics>
struct ConductivityKernel<gps, 2, HeatTransferElementType::ASYMMETRIC_GENERAL, Physics>: Conductivity, Physics {
	ConductivityKernel(const Conductivity &base)
	: Conductivity(base),
	  kxx(this->conductivity->values.get(0, 0).evaluator),
	  kxy(this->conductivity->values.get(0, 1).evaluator),
	  kyx(this->conductivity->values.get(1, 0).evaluator),
	  kyy(this->conductivity->values.get(1, 1).evaluator)
	{
		isconst =
				kxx->isConst() && kxy->isConst() &&
				kyx->isConst() && kyy->isConst();
	}

	Evaluator *kxx, *kxy;
	Evaluator *kyx, *kyy;

	void simd(typename Physics::Element &element)
	{
		if (direct) {
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.conductivity[gp][0][s] = kxx->evaluate();
					element.conductivity[gp][1][s] = kxy->evaluate();
					element.conductivity[gp][2][s] = kyx->evaluate();
					element.conductivity[gp][3][s] = kyy->evaluate();
				}
			}
		} else {
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][0][s] = kxx->evaluate();
					element.ecf.conductivity[gp][1][s] = kxy->evaluate();
					element.ecf.conductivity[gp][2][s] = kyx->evaluate();
					element.ecf.conductivity[gp][3][s] = kyy->evaluate();
				}
			}
		}
	}
};

template <size_t gps, class Physics>
struct ConductivityKernel<gps, 3, HeatTransferElementType::ASYMMETRIC_GENERAL, Physics>: Conductivity, Physics {
	ConductivityKernel(const Conductivity &base)
	: Conductivity(base),
	  kxx(this->conductivity->values.get(0, 0).evaluator),
	  kxy(this->conductivity->values.get(0, 1).evaluator),
	  kxz(this->conductivity->values.get(0, 2).evaluator),
	  kyx(this->conductivity->values.get(1, 0).evaluator),
	  kyy(this->conductivity->values.get(1, 1).evaluator),
	  kyz(this->conductivity->values.get(1, 2).evaluator),
	  kzx(this->conductivity->values.get(2, 0).evaluator),
	  kzy(this->conductivity->values.get(2, 1).evaluator),
	  kzz(this->conductivity->values.get(2, 2).evaluator)
	{
		isconst =
				kxx->isConst() && kxy->isConst() && kxz->isConst() &&
				kyx->isConst() && kyy->isConst() && kyz->isConst() &&
				kzx->isConst() && kzy->isConst() && kzz->isConst();
	}

	Evaluator *kxx, *kxy, *kxz;
	Evaluator *kyx, *kyy, *kyz;
	Evaluator *kzx, *kzy, *kzz;

	void simd(typename Physics::Element &element)
	{
		if (direct) {
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.conductivity[gp][0][s] = kxx->evaluate();
					element.conductivity[gp][1][s] = kxy->evaluate();
					element.conductivity[gp][2][s] = kxz->evaluate();
					element.conductivity[gp][3][s] = kyx->evaluate();
					element.conductivity[gp][4][s] = kyy->evaluate();
					element.conductivity[gp][5][s] = kyz->evaluate();
					element.conductivity[gp][6][s] = kzx->evaluate();
					element.conductivity[gp][7][s] = kzy->evaluate();
					element.conductivity[gp][8][s] = kzz->evaluate();
				}
			}
		} else {
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][0][s] = kxx->evaluate();
					element.ecf.conductivity[gp][1][s] = kxy->evaluate();
					element.ecf.conductivity[gp][2][s] = kxz->evaluate();
					element.ecf.conductivity[gp][3][s] = kyx->evaluate();
					element.ecf.conductivity[gp][4][s] = kyy->evaluate();
					element.ecf.conductivity[gp][5][s] = kyz->evaluate();
					element.ecf.conductivity[gp][6][s] = kzx->evaluate();
					element.ecf.conductivity[gp][7][s] = kzy->evaluate();
					element.ecf.conductivity[gp][8][s] = kzz->evaluate();
				}
			}
		}

	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_CONDUCTIVITY_H_ */
