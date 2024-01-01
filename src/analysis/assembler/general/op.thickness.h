
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_THICKNESS_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_THICKNESS_H_

#include "subkernel.h"
#include "basis/containers/serializededata.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

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
		action = SubKernel::PREPROCESS | SubKernel::ASSEMBLE | SubKernel::REASSEMBLE | SubKernel::ITERATION | SubKernel::SOLUTION;
	}

	void activate(ECFExpression *expression, serializededata<esint, esint>::const_iterator enodes, serializededata<esint, esint>::const_iterator end, double * target)
	{
		this->expression = expression;
		this->enodes = enodes;
		this->end = end;
		this->target = target;
		if (this->expression) {
			this->isconst = this->expression->evaluator->isConst();
			this->isactive = !this->isconst && this->expression->isset;
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

template <size_t nodes, size_t ndim> struct ThicknessToGp;

template <size_t nodes>
struct ThicknessToGp<nodes, 2>: Thickness {
	ThicknessToGp(const Thickness &base): Thickness(base) {}

	template <typename Element>
	void simd(Element &element, size_t gp)
	{
		element.thickness.gp = zeros();
		for (size_t n = 0; n < nodes; ++n) {
			element.thickness.gp = element.thickness.gp + load1(element.N[gp][n]) * element.thickness.node[n];
		}
	}
};

template <size_t nodes>
struct ThicknessToGp<nodes, 3>: Thickness {
	ThicknessToGp(const Thickness &base): Thickness(base) {}

	template <typename Element>
	void simd(Element &element, size_t gp)
	{

	}
};

template <size_t nodes, size_t ndim> struct ThicknessToNodes;

template <size_t nodes>
struct ThicknessToNodes<nodes, 2>: Thickness {
	ThicknessToNodes(const Thickness &base): Thickness(base) {}

	template <typename Element>
	void simd(Element &element)
	{
		for (size_t s = 0; s < SIMD::size; ++s, ++enodes) {
			if (enodes == end) break;
			for (size_t n = 0; n < nodes; ++n) {
				target[enodes->at(n)] = element.thickness.node[n][s];
			}
		}
	}
};

template <size_t nodes>
struct ThicknessToNodes<nodes, 3>: Thickness {
	ThicknessToNodes(const Thickness &base): Thickness(base) {}

	template <typename Element>
	void simd(Element &element)
	{

	}
};

template <size_t nodes, size_t ndim> struct ThicknessFromNodes;

template <size_t nodes>
struct ThicknessFromNodes<nodes, 2>: Thickness {
	ThicknessFromNodes(const Thickness &base): Thickness(base) {}

	template <typename Element>
	void simd(Element &element)
	{
		for (size_t s = 0; s < SIMD::size; ++s, ++enodes) {
			if (enodes == end) break;
			for (size_t n = 0; n < nodes; ++n) {
				element.thickness.node[n][s] = target[enodes->at(n)];
			}
		}
	}
};

template <size_t nodes>
struct ThicknessFromNodes<nodes, 3>: Thickness {
	ThicknessFromNodes(const Thickness &base): Thickness(base) {}

	template <typename Element>
	void simd(Element &element)
	{

	}
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_THICKNESS_H_ */
