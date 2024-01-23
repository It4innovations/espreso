
#ifndef SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_OP_MATRIX_MASS_H_
#define SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_OP_MATRIX_MASS_H_

#include "analysis/assembler/general/subkernel.h"
#include "math/primitives/matrix_info.h"
#include "wrappers/simd/simd.h"


namespace espreso {

struct MatrixMass: public SubKernel {
	const char* name() const { return "MatrixMass"; }

	MatrixMass()
	{
		isconst = false;
		action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE;
	}

	void activate()
	{
		this->isactive = 1;
	}
};

template <size_t nodes>
struct MatrixMassKernel: MatrixMass {
	MatrixMassKernel(const MatrixMass &base): MatrixMass(base) {}

	template <typename Element>
	void simd(Element &element, size_t gp)
	{
		SIMD scale = element.det * load1(element.w[gp]) * element.ecf.density;
		for (size_t n = 0; n < nodes; ++n) {
			element.M[(n * nodes + n)] = element.M[(n * nodes + n)] + scale * load1(element.N[gp][n] * element.N[gp][n]);
			for (size_t m = n + 1; m < nodes; ++m) {
				SIMD k = scale * load1(element.N[gp][n] * element.N[gp][m]);
				element.M[(n * nodes + m)] = element.M[(n * nodes + m)] + k;
				element.M[(m * nodes + n)] = element.M[(m * nodes + n)] + k;
			}
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_OP_MATRIX_MASS_H_ */
