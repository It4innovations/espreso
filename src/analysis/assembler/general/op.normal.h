
#ifndef SRC_ANALYSIS_ASSEMBLER_GENERAL_OP_NORMAL_H_
#define SRC_ANALYSIS_ASSEMBLER_GENERAL_OP_NORMAL_H_

#include "element.h"
#include "subkernel.h"
#include "basis/containers/serializededata.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

struct Normal: SubKernel {
	const char* name() const { return "Normal"; }

	serializededata<esint, esint>::const_iterator enodes, end;
	double *target, *multiplicity;

	Normal()
	: enodes(info::mesh->elements->nodes->cbegin()),
	  end(info::mesh->elements->nodes->cend()),
	  target(nullptr),
	  multiplicity(nullptr)
	{
		isconst = false;
		action = SubKernel::PREPROCESS;
	}

	void activate(serializededata<esint, esint>::const_iterator enodes, serializededata<esint, esint>::const_iterator end, double *target, double *multiplicity)
	{
		this->enodes = enodes;
		this->end = end;
		this->target = target;
		this->multiplicity = multiplicity;
		this->isactive = 1;
	}
};

template <size_t nodes, size_t ndim>
struct StoreNormalKernel: Normal {
	StoreNormalKernel(const Normal &base): Normal(base) {}

	template <typename Element>
	void simd(Element &element)
	{
		// TODO: compute normal in nodes
		for (size_t s = 0; s < SIMD::size; ++s, ++enodes) {
			if (enodes == end) break;
			for (size_t n = 0; n < nodes; ++n) {
				for (size_t d = 0; d < ndim; ++d) {
					target[ndim * enodes->at(n) + d] += element.normal[d][s] * multiplicity[enodes->at(n)];
				}
			}
		}
	}
};

}
#endif /* SRC_ANALYSIS_ASSEMBLER_GENERAL_OP_NORMAL_H_ */
