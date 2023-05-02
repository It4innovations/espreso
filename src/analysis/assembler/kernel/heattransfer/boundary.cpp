
#include "analysis/assembler/module/heattransfer.h"
#include "analysis/assembler/module/assembler.hpp"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"

#include "analysis/scheme/steadystate.h"
#include "math/physics/matrix_distributed.h"

#include <numeric>
#include <algorithm>

namespace espreso {

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void preprocess(HeatTransfer::BoundarySubKernels &subkernels)
{

}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void compute(const HeatTransfer::BoundarySubKernels &subkernels, Assembler::Action action)
{

}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void fill(const HeatTransfer::BoundarySubKernels &subkernels)
{

}

template <size_t ndim>
void initDirichlet(HeatTransfer::BoundarySubKernels &subkernels)
{
	typedef HeatTransferBoundaryDescriptor<1, 1, ndim, 0> Physics;
	if (subkernels.temperature.expression) {
		subkernels.expressions.push_back(new ExternalNodeExpression<1, Physics>(
				subkernels.temperature.expression->evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.temp[0][s] = value; }));
	}
}

template <size_t ndim>
void dirichlet(const HeatTransfer::BoundarySubKernels &subkernels)
{
	typedef HeatTransferBoundaryDescriptor<1, 1, ndim, 0> Physics;
	typename Physics::Element element;

	CoordinatesKernel<1, 1, ndim, Physics> coordinates(subkernels.coordinates);
	VectorSetterKernel<1, Physics> set(subkernels.dirichlet, [] (auto &element, size_t &n, size_t &d, size_t &s) { return element.temp[n][s]; });

	std::vector<ExternalNodeExpression<1, Physics>*> nonconst;
	for (size_t i = 0; i < subkernels.expressions.size(); ++i) {
		if (subkernels.expressions[i]->evaluator->isConst()) {
			dynamic_cast<ExternalNodeExpression<1, Physics>*>(subkernels.expressions[i])->simd(element);
		} else {
			nonconst.push_back(dynamic_cast<ExternalNodeExpression<1, Physics>*>(subkernels.expressions[i]));
		}
	}

	for (esint c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
		for (size_t i = 0; i < nonconst.size(); ++i) {
			nonconst[i]->simd(element);
		}
		set.simd(element);
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void runAction(HeatTransfer::BoundarySubKernels &subkernels, Assembler::Action action)
{
	switch (action) {
	case Assembler::Action::PREPROCESS: preprocess<code, nodes, gps, ndim, edim>(subkernels); break;
	case Assembler::Action::FILL: fill<code, nodes, gps, ndim, edim>(subkernels); break;
	default: compute<code, nodes, gps, ndim, edim>(subkernels, action); break;
	}
}

template <size_t ndim, size_t edim> void addBC(HeatTransfer::BoundarySubKernels &subkernels, Assembler::Action action);

template <> void addBC<2, 1>(HeatTransfer::BoundarySubKernels &subkernels, Assembler::Action action)
{

}

template <> void addBC<2, 0>(HeatTransfer::BoundarySubKernels &subkernels, Assembler::Action action)
{
	switch (action) {
	case Assembler::Action::PREPROCESS: initDirichlet<2>(subkernels); break;
	case Assembler::Action::ASSEMBLE: case Assembler::Action::REASSEMBLE: dirichlet<2>(subkernels); break;
	}
}
template <> void addBC<3, 2>(HeatTransfer::BoundarySubKernels &subkernels, Assembler::Action action)
{

}

template <> void addBC<3, 1>(HeatTransfer::BoundarySubKernels &subkernels, Assembler::Action action)
{

}

template <> void addBC<3, 0>(HeatTransfer::BoundarySubKernels &subkernels, Assembler::Action action)
{
	switch (action) {
	case Assembler::Action::PREPROCESS: initDirichlet<3>(subkernels); break;
	case Assembler::Action::ASSEMBLE: case Assembler::Action::REASSEMBLE: dirichlet<3>(subkernels); break;
	}
}

void HeatTransfer::runBoundary(Action action, size_t region, size_t interval)
{
	switch (info::mesh->dimension) {
	case 2:
		switch (info::mesh->boundaryRegions[region]->dimension) {
		case 0: addBC<2, 0>(boundary[region][interval], action); break;
		case 1: addBC<2, 1>(boundary[region][interval], action); break;
		} break;
	case 3:
		switch (info::mesh->boundaryRegions[region]->dimension) {
		case 0: addBC<3, 0>(boundary[region][interval], action); break;
		case 1: addBC<3, 1>(boundary[region][interval], action); break;
		case 2: addBC<3, 2>(boundary[region][interval], action); break;
		} break;
	}
}

}


