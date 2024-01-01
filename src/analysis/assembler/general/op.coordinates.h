
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_COORDINATES_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_COORDINATES_H_

#include "element.h"
#include "subkernel.h"
#include "basis/containers/serializededata.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
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
		action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE | SubKernel::ITERATION | SubKernel::SOLUTION;
	}

	void activate(serializededata<esint, esint>::const_iterator enodes, serializededata<esint, esint>::const_iterator end, bool toGPs)
	{
		this->enodes = enodes;
		this->end = end;
		this->toGPs = toGPs;
	}
};

template <size_t nodes, size_t ndim>
struct CoordinatesKernel: Coordinates {
	CoordinatesKernel(const Coordinates &base): Coordinates(base) {}

	template <typename Element>
	void simd(Element &element)
	{
		for (size_t s = 0; s < SIMD::size; ++s, ++enodes) {
			if (enodes == end) {
				for (size_t n = 0; n < nodes; ++n) {
					for (size_t d = 0; d < ndim; ++d) {
						element.coords.node[n][d][s] = 0;
					}
				}
				break;
			} else {
				for (size_t n = 0; n < nodes; ++n) {
					for (size_t d = 0; d < ndim; ++d) {
						element.coords.node[n][d][s] = info::mesh->nodes->coordinates->datatarray()[enodes->at(n)][d];
					}
				}
			}
		}
	}
};

template <size_t nodes, size_t ndim>
struct CoordinatesToGPsKernel: Coordinates {
	CoordinatesToGPsKernel(const Coordinates &base): Coordinates(base)
	{
		this->isactive = toGPs;
	}

	template <typename Element>
	void simd(Element &element, size_t gp)
	{
		for (size_t d = 0; d < ndim; ++d) {
			element.coords.gp[d] = zeros();
		}
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t d = 0; d < ndim; ++d) {
				element.coords.gp[d] = element.coords.gp[d] + load1(element.N[gp][n]) * element.coords.node[n][d];
			}
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_COORDINATES_H_ */
