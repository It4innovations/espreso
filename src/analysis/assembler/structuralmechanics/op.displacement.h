
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_DISPLACEMENT_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_DISPLACEMENT_H_

#include "basis/containers/serializededata.h"
#include "analysis/assembler/general/element.h"
#include "analysis/assembler/general/subkernel.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

struct Displacement: SubKernel {
	Displacement()
	: enodes(info::mesh->elements->nodes->cbegin()),
	  end(info::mesh->elements->nodes->cend()),
	  source(nullptr)
	{
		isconst = false;
		action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE | SubKernel::ITERATION | SubKernel::SOLUTION;
	}

	void activate(serializededata<esint, esint>::const_iterator enodes, serializededata<esint, esint>::const_iterator end, const double * source)
	{
		this->enodes = enodes;
		this->end = end;
		this->source = source;
		this->isactive = 1;
	}

	serializededata<esint, esint>::const_iterator enodes, end;
	const double * source;
};

template <size_t nodes, size_t ndim> struct DisplacementKernel;

template <size_t nodes>
struct DisplacementKernel<nodes, 2>: Displacement {
	DisplacementKernel(const Displacement &base): Displacement(base) {}

	template <typename Element>
	void simd(Element &element)
	{
		for (size_t s = 0; s < SIMD::size; ++s, ++enodes) {
			if (enodes == end) break;
			for (size_t n = 0; n < nodes; ++n) {
				element.displacement[n][0][s] = source[2 * enodes->at(n) + 0];
				element.displacement[n][1][s] = source[2 * enodes->at(n) + 1];
			}
		}
	}
};


template <size_t nodes>
struct DisplacementKernel<nodes, 3>: Displacement {
	DisplacementKernel(const Displacement &base): Displacement(base) {}

	template <typename Element>
	void simd(Element &element)
	{
		for (size_t s = 0; s < SIMD::size; ++s, ++enodes) {
			if (enodes == end) break;
			for (size_t n = 0; n < nodes; ++n) {
				element.displacement[n][0][s] = source[3 * enodes->at(n) + 0];
				element.displacement[n][1][s] = source[3 * enodes->at(n) + 1];
				element.displacement[n][2][s] = source[3 * enodes->at(n) + 2];
			}
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_DISPLACEMENT_H_ */
