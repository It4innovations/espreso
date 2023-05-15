
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_DISPLACEMENT_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_DISPLACEMENT_H_

#include "subkernel.h"

namespace espreso {

struct Displacement: SubKernel {
	serializededata<esint, esint>::const_iterator enodes, end;
	const double * source;
	bool smallStrainTensor;

	Displacement()
	: enodes(info::mesh->elements->nodes->cbegin()),
	  end(info::mesh->elements->nodes->cend()),
	  source(nullptr), smallStrainTensor(false)
	{
		isconst = false;
		action = Assembler::ASSEMBLE | Assembler::REASSEMBLE | Assembler::ITERATION | Assembler::SOLUTION;
	}

	void activate(serializededata<esint, esint>::const_iterator enodes, serializededata<esint, esint>::const_iterator end, const double * source, bool smallStrainTensor)
	{
		this->enodes = enodes;
		this->end = end;
		this->source = source;
		this->smallStrainTensor = smallStrainTensor;
		this->isactive = 1;
	}
};

template <size_t nodes, size_t gps, size_t ndim, class Physics>
struct DisplacementKernel: Displacement, Physics {
	DisplacementKernel(const Displacement &base): Displacement(base) {}

	// B = dX  0  0    disp = x
	//      0 dY  0           y
	//      0  0 dZ           z
	//     dY dX  0
	//      0 dZ dY
	//     dZ  0 dX
	void simd(typename Physics::Element &element)
	{
		for (size_t s = 0; s < SIMD::size; ++s, ++enodes) {
			if (enodes == end) break;
			for (size_t n = 0; n < nodes; ++n) {
				for (size_t d = 0; d < ndim; ++d) {
					element.displacement[n][d][s] = source[ndim * enodes->at(n) + d];
				}
			}
		}
		if (smallStrainTensor) {
			for (size_t gp = 0; gp < gps; ++gp) {
				element.smallStrainTensor[gp][0] = zeros();
				element.smallStrainTensor[gp][1] = zeros();
				element.smallStrainTensor[gp][2] = zeros();
				element.smallStrainTensor[gp][3] = zeros();
				element.smallStrainTensor[gp][4] = zeros();
				element.smallStrainTensor[gp][5] = zeros();
				for (size_t n = 0; n < nodes; ++n) {
					element.smallStrainTensor[gp][0] = element.smallStrainTensor[gp][0] + element.dND[gp][n][0] * element.displacement[n][0];
					element.smallStrainTensor[gp][1] = element.smallStrainTensor[gp][1] + element.dND[gp][n][1] * element.displacement[n][1];
					element.smallStrainTensor[gp][2] = element.smallStrainTensor[gp][2] + element.dND[gp][n][2] * element.displacement[n][2];

					element.smallStrainTensor[gp][3] = element.smallStrainTensor[gp][3] + element.dND[gp][n][1] * element.displacement[n][0];
					element.smallStrainTensor[gp][3] = element.smallStrainTensor[gp][3] + element.dND[gp][n][0] * element.displacement[n][1];

					element.smallStrainTensor[gp][4] = element.smallStrainTensor[gp][4] + element.dND[gp][n][2] * element.displacement[n][1];
					element.smallStrainTensor[gp][4] = element.smallStrainTensor[gp][4] + element.dND[gp][n][1] * element.displacement[n][2];

					element.smallStrainTensor[gp][5] = element.smallStrainTensor[gp][5] + element.dND[gp][n][2] * element.displacement[n][0];
					element.smallStrainTensor[gp][5] = element.smallStrainTensor[gp][5] + element.dND[gp][n][0] * element.displacement[n][2];
				}
			}
		}
	}
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_DISPLACEMENT_H_ */
