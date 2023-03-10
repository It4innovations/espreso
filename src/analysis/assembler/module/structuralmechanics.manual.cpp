
#include "structuralmechanics.h"
#include "structuralmechanics.generator.h"
#include "assembler.hpp"

#include "analysis/assembler/operators/info.h"
#include "analysis/assembler/operators/basis.h"
#include "analysis/assembler/operators/coordinates.h"
#include "analysis/assembler/operators/expression.h"
#include "analysis/assembler/operators/integration.h"
#include "analysis/assembler/operators/structuralmechanics.f.h"
#include "analysis/assembler/operators/structuralmechanics.K.h"
#include "analysis/assembler/operators/filler.h"

#include "basis/expression/variable.h"
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

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
Assembler::measurements StructuralMechanics::manualloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	double initStart, initEnd;
	eslog::info("       = LOOP TYPE                                                            MANUAL = \n");
	if (elements == 0) return {0.0, 0.0};
	initStart = eslog::time();
	typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element element;
	std::vector<DataDescriptor<nodes, gps, ndim, edim, etype>*> active; active.reserve(ops.size());

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
		if ((*op)->action & action) {
			if (elements > SIMD::size) {
				if ((*op)->isconst) {
					dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, etype>*>(*op)->simd(element);
				} else {
					active.push_back(dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, etype>*>(*op));
					active.back()->simd(element);
				}
			} else {
				dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, etype>*>(*op)->peel(element, elements);
			}
		}
	}
	initEnd = eslog::time();
	double start = eslog::time();
	esint chunks = elements / SIMD::size;
	for (esint c = 1; c < chunks; ++c) {
		for (auto op = active.cbegin(); op != active.cend(); ++op) {
			(*op)->simd(element);
		}
	}
	double end = eslog::time();

	if (elements % SIMD::size) {
		for (auto op = active.cbegin(); op != active.cend(); ++op) {
			(*op)->peel(element, elements % SIMD::size);
		}
	}

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
		if ((*op)->action & action) {
			if ((*op)->isconst) {
				(*op)->move(-(int)std::min(elements, (esint)SIMD::size));
			} else {
				(*op)->move(-elements);
			}
		}
	}
	return {initEnd - initStart, end - start};
}

template <int etype>
Assembler::measurements StructuralMechanics::instantiateManual2D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return manualloop<StructuralMechanicsDataDescriptor, 3, StructuralMechanicsGPC::TRIANGLE3, 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return manualloop<StructuralMechanicsDataDescriptor, 6, StructuralMechanicsGPC::TRIANGLE6, 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE4):   return manualloop<StructuralMechanicsDataDescriptor, 4, StructuralMechanicsGPC::SQUARE4  , 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE8):   return manualloop<StructuralMechanicsDataDescriptor, 8, StructuralMechanicsGPC::SQUARE8  , 2, 2, etype>(action, ops, interval, elements); break;
	default: return {0.0, 0.0};
	}
}

template <int etype>
Assembler::measurements StructuralMechanics::instantiateManual3D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TETRA4):    return manualloop<StructuralMechanicsDataDescriptor,  4, StructuralMechanicsGPC::TETRA4    , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::TETRA10):   return manualloop<StructuralMechanicsDataDescriptor, 10, StructuralMechanicsGPC::TETRA10   , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID5):  return manualloop<StructuralMechanicsDataDescriptor,  5, StructuralMechanicsGPC::PYRAMID5  , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID13): return manualloop<StructuralMechanicsDataDescriptor, 13, StructuralMechanicsGPC::PYRAMID13 , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA6):   return manualloop<StructuralMechanicsDataDescriptor,  6, StructuralMechanicsGPC::PRISMA6   , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA15):  return manualloop<StructuralMechanicsDataDescriptor, 15, StructuralMechanicsGPC::PRISMA15  , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA8):     return manualloop<StructuralMechanicsDataDescriptor,  8, StructuralMechanicsGPC::HEXA8     , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA20):    return manualloop<StructuralMechanicsDataDescriptor, 20, StructuralMechanicsGPC::HEXA20    , 3, 3, etype>(action, ops, interval, elements); break;
	default: return {0.0, 0.0};
	}
}

Assembler::measurements StructuralMechanics::instantiateManual(ActionOperator::Action action, int code, int etype, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (info::mesh->dimension) {
	case 2:
		switch (etype) {
		// elements
		case StructuralMechanicsElementType::SYMMETRIC_PLANE             : return instantiateManual2D<StructuralMechanicsElementType::SYMMETRIC_PLANE             >(action, code, ops, interval, elements);
		case StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC: return instantiateManual2D<StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC>(action, code, ops, interval, elements);
		}
	case 3:
		switch (etype) {
		// elements
		case StructuralMechanicsElementType::SYMMETRIC_VOLUME: return instantiateManual3D<StructuralMechanicsElementType::SYMMETRIC_VOLUME>(action, code, ops, interval, elements);
		}
	}
	return {0.0, 0.0};
}

}

