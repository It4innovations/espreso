
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_THICKNESS_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_THICKNESS_H_

#include "subkernel.h"
#include "basis/containers/serializededata.h"
#include "esinfo/meshinfo.h"

namespace espreso {

struct Thickness: SubKernel {
	const char* name() const { return "Thickness"; }

	ECFExpression *expression;
	serializededata<esint, esint>::const_iterator enodes, end;
	double *target;

	Thickness()
	: expression(nullptr),
	  enodes(info::mesh->elements->nodes->cbegin()),
	  end(info::mesh->elements->nodes->cend()),
	  target(nullptr)
	{
		isconst = false;
		action = Assembler::ASSEMBLE | Assembler::REASSEMBLE | Assembler::ITERATION | Assembler::SOLUTION;
	}

	void activate(ECFExpression *expression, serializededata<esint, esint>::const_iterator enodes, serializededata<esint, esint>::const_iterator end, double * target)
	{
		this->expression = expression;
		this->enodes = enodes;
		this->end = end;
		this->target = target;
		if (this->expression) {
			this->isactive = 1;
		}
	}

	void activate(serializededata<esint, esint>::const_iterator enodes, serializededata<esint, esint>::const_iterator end, double * target)
	{
		this->enodes = enodes;
		this->end = end;
		this->target = target;
		this->isactive = 1;
	}
};

template <size_t nodes, size_t ndim, class Physics> struct ThicknessToNodes;

template <size_t nodes, class Physics>
struct ThicknessToNodes<nodes, 2, Physics>: Thickness, Physics {
	ThicknessToNodes(const Thickness &base): Thickness(base) {}

	void simd(typename Physics::Element &element)
	{
		for (size_t s = 0; s < SIMD::size; ++s, ++enodes) {
			if (enodes == end) break;
			for (size_t n = 0; n < nodes; ++n) {
				target[enodes->at(n)] = element.ecf.thickness[n][s];
			}
		}
	}
};

template <size_t nodes, class Physics>
struct ThicknessToNodes<nodes, 3, Physics>: Thickness, Physics {
	ThicknessToNodes(const Thickness &base): Thickness(base) {}

	void simd(typename Physics::Element &element)
	{

	}
};

template <size_t nodes, size_t ndim, class Physics> struct ThicknessFromNodes;

template <size_t nodes, class Physics>
struct ThicknessFromNodes<nodes, 2, Physics>: Thickness, Physics {
	ThicknessFromNodes(const Thickness &base): Thickness(base) {}

	void simd(typename Physics::Element &element)
	{
		for (size_t s = 0; s < SIMD::size; ++s, ++enodes) {
			if (enodes == end) break;
			for (size_t n = 0; n < nodes; ++n) {
				element.ecf.thickness[n][s] = target[enodes->at(n)];
			}
		}
	}
};

template <size_t nodes, class Physics>
struct ThicknessFromNodes<nodes, 3, Physics>: Thickness, Physics {
	ThicknessFromNodes(const Thickness &base): Thickness(base) {}

	void simd(typename Physics::Element &element)
	{

	}
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_THICKNESS_H_ */
