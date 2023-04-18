
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_COORDINATES_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_COORDINATES_H_

#include "subkernel.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"

namespace espreso {

struct Coordinates: SubKernel {
	const char* name() const { return "Coordinates"; }

	serializededata<esint, esint>::const_iterator enodes, end;
	bool toGPs;

	Coordinates()
	: enodes(info::mesh->elements->nodes->cbegin()),
	  end(info::mesh->elements->nodes->cend()),
	  toGPs(false)
	{
		isconst = false;
		action = Assembler::ASSEMBLE | Assembler::REASSEMBLE | Assembler::SOLUTION;
	}

	void activate(serializededata<esint, esint>::const_iterator enodes, serializededata<esint, esint>::const_iterator end, bool toGPs)
	{
		this->enodes = enodes;
		this->end = end;
		this->toGPs = toGPs;
	}
};

template <size_t nodes, size_t gps, size_t ndim, class Physics>
struct CoordinatesKernel: Coordinates, Physics {
	CoordinatesKernel(const Coordinates &base): Coordinates(base) {}

	void simd(typename Physics::Element &element)
	{
		for (size_t s = 0; s < SIMD::size; ++s, ++enodes) {
			if (enodes == end) break;
			for (size_t n = 0; n < nodes; ++n) {
				for (size_t d = 0; d < ndim; ++d) {
					element.coords[n][d][s] = info::mesh->nodes->coordinates->datatarray()[enodes->at(n)][d];
				}
			}
		}
		if (toGPs) {
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t d = 0; d < ndim; ++d) {
					element.gpcoords[gp][d] = zeros();
				}
				for (size_t n = 0; n < nodes; ++n) {
					for (size_t d = 0; d < ndim; ++d) {
						element.gpcoords[gp][d] = element.gpcoords[gp][d] + load1(element.N[gp][n]) * element.coords[n][d];
					}
				}
			}
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_COORDINATES_H_ */
