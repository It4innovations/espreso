
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_STRUCTURALMECHANICS_GENERATOR_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_STRUCTURALMECHANICS_GENERATOR_H_

#include "structuralmechanics.element.h"
#include "analysis/assembler/operator.h"
#include "analysis/assembler/operators/basis.h"
#include "analysis/assembler/operators/expression.h"

#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"

namespace espreso {

template <size_t ndim, int etype>
static void generateBaseFunctions(size_t interval, std::vector<std::vector<ActionOperator*> > &ops)
{
	switch (static_cast<Element::CODE>(info::mesh->elements->eintervals[interval].code)) {
	case Element::CODE::LINE2    : ops[interval].push_back(new Basis<Element::CODE::LINE2    ,  2, StructuralMechanicsGPC::LINE2    , 1, etype, StructuralMechanicsDataDescriptor< 2, StructuralMechanicsGPC::LINE2    , ndim, 1, etype> >()); break;
	case Element::CODE::LINE3    : ops[interval].push_back(new Basis<Element::CODE::LINE3    ,  3, StructuralMechanicsGPC::LINE3    , 1, etype, StructuralMechanicsDataDescriptor< 3, StructuralMechanicsGPC::LINE3    , ndim, 1, etype> >()); break;
	case Element::CODE::TRIANGLE3: ops[interval].push_back(new Basis<Element::CODE::TRIANGLE3,  3, StructuralMechanicsGPC::TRIANGLE3, 2, etype, StructuralMechanicsDataDescriptor< 3, StructuralMechanicsGPC::TRIANGLE3, ndim, 2, etype> >()); break;
	case Element::CODE::TRIANGLE6: ops[interval].push_back(new Basis<Element::CODE::TRIANGLE6,  6, StructuralMechanicsGPC::TRIANGLE6, 2, etype, StructuralMechanicsDataDescriptor< 6, StructuralMechanicsGPC::TRIANGLE6, ndim, 2, etype> >()); break;
	case Element::CODE::SQUARE4  : ops[interval].push_back(new Basis<Element::CODE::SQUARE4  ,  4, StructuralMechanicsGPC::SQUARE4  , 2, etype, StructuralMechanicsDataDescriptor< 4, StructuralMechanicsGPC::SQUARE4  , ndim, 2, etype> >()); break;
	case Element::CODE::SQUARE8  : ops[interval].push_back(new Basis<Element::CODE::SQUARE8  ,  8, StructuralMechanicsGPC::SQUARE8  , 2, etype, StructuralMechanicsDataDescriptor< 8, StructuralMechanicsGPC::SQUARE8  , ndim, 2, etype> >()); break;
	case Element::CODE::TETRA4   : ops[interval].push_back(new Basis<Element::CODE::TETRA4   ,  4, StructuralMechanicsGPC::TETRA4   , 3, etype, StructuralMechanicsDataDescriptor< 4, StructuralMechanicsGPC::TETRA4   , ndim, 3, etype> >()); break;
	case Element::CODE::TETRA10  : ops[interval].push_back(new Basis<Element::CODE::TETRA10  , 10, StructuralMechanicsGPC::TETRA10  , 3, etype, StructuralMechanicsDataDescriptor<10, StructuralMechanicsGPC::TETRA10  , ndim, 3, etype> >()); break;
	case Element::CODE::PYRAMID5 : ops[interval].push_back(new Basis<Element::CODE::PYRAMID5 ,  5, StructuralMechanicsGPC::PYRAMID5 , 3, etype, StructuralMechanicsDataDescriptor< 5, StructuralMechanicsGPC::PYRAMID5 , ndim, 3, etype> >()); break;
	case Element::CODE::PYRAMID13: ops[interval].push_back(new Basis<Element::CODE::PYRAMID13, 13, StructuralMechanicsGPC::PYRAMID13, 3, etype, StructuralMechanicsDataDescriptor<13, StructuralMechanicsGPC::PYRAMID13, ndim, 3, etype> >()); break;
	case Element::CODE::PRISMA6  : ops[interval].push_back(new Basis<Element::CODE::PRISMA6  ,  6, StructuralMechanicsGPC::PRISMA6  , 3, etype, StructuralMechanicsDataDescriptor< 6, StructuralMechanicsGPC::PRISMA6  , ndim, 3, etype> >()); break;
	case Element::CODE::PRISMA15 : ops[interval].push_back(new Basis<Element::CODE::PRISMA15 , 15, StructuralMechanicsGPC::PRISMA15 , 3, etype, StructuralMechanicsDataDescriptor<15, StructuralMechanicsGPC::PRISMA15 , ndim, 3, etype> >()); break;
	case Element::CODE::HEXA8    : ops[interval].push_back(new Basis<Element::CODE::HEXA8    ,  8, StructuralMechanicsGPC::HEXA8    , 3, etype, StructuralMechanicsDataDescriptor< 8, StructuralMechanicsGPC::HEXA8    , ndim, 3, etype> >()); break;
	case Element::CODE::HEXA20   : ops[interval].push_back(new Basis<Element::CODE::HEXA20   , 20, StructuralMechanicsGPC::HEXA20   , 3, etype, StructuralMechanicsDataDescriptor<20, StructuralMechanicsGPC::HEXA20   , ndim, 3, etype> >()); break;
	}
}

static void generateBaseFunctions(const std::vector<int> &etype, std::vector<std::vector<ActionOperator*> > &ops)
{
	GaussPoints<Element::CODE::LINE2    ,  2, StructuralMechanicsGPC::LINE2    , 1>::set();
	GaussPoints<Element::CODE::TRIANGLE3,  3, StructuralMechanicsGPC::TRIANGLE3, 2>::set();
	GaussPoints<Element::CODE::SQUARE4  ,  4, StructuralMechanicsGPC::SQUARE4  , 2>::set();
	GaussPoints<Element::CODE::TETRA4   ,  4, StructuralMechanicsGPC::TETRA4   , 3>::set();
	GaussPoints<Element::CODE::PYRAMID5 ,  5, StructuralMechanicsGPC::PYRAMID5 , 3>::set();
	GaussPoints<Element::CODE::PRISMA6  ,  6, StructuralMechanicsGPC::PRISMA6  , 3>::set();
	GaussPoints<Element::CODE::HEXA8    ,  8, StructuralMechanicsGPC::HEXA8    , 3>::set();
	GaussPoints<Element::CODE::LINE3    ,  3, StructuralMechanicsGPC::LINE3    , 1>::set();
	GaussPoints<Element::CODE::TRIANGLE6,  6, StructuralMechanicsGPC::TRIANGLE6, 2>::set();
	GaussPoints<Element::CODE::SQUARE8  ,  8, StructuralMechanicsGPC::SQUARE8  , 2>::set();
	GaussPoints<Element::CODE::TETRA10  , 10, StructuralMechanicsGPC::TETRA10  , 3>::set();
	GaussPoints<Element::CODE::PYRAMID13, 13, StructuralMechanicsGPC::PYRAMID13, 3>::set();
	GaussPoints<Element::CODE::PRISMA15 , 15, StructuralMechanicsGPC::PRISMA15 , 3>::set();
	GaussPoints<Element::CODE::HEXA20   , 20, StructuralMechanicsGPC::HEXA20   , 3>::set();

	for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		switch (info::mesh->dimension) {
		case 2:
			switch (etype[i]) {
			case StructuralMechanicsElementType::SYMMETRIC_PLANE             : generateBaseFunctions<2, StructuralMechanicsElementType::SYMMETRIC_PLANE             >(i, ops); break;
			case StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC: generateBaseFunctions<2, StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC>(i, ops); break;
			} break;
		case 3:
			switch (etype[i]) {
			case StructuralMechanicsElementType::SYMMETRIC_VOLUME:             generateBaseFunctions<3, StructuralMechanicsElementType::SYMMETRIC_VOLUME            >(i, ops); break;
			} break;
		}
	}
}

static void generateBaseFunctions(int axisymmetric, const std::vector<int> &bfilter, std::vector<std::vector<std::vector<ActionOperator*> > > &ops)
{
	for(size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (bfilter[r]) {
			if (axisymmetric) {

			}
			for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
				switch (info::mesh->boundaryRegions[r]->dimension) {
				case 2:
					switch (static_cast<Element::CODE>(info::mesh->boundaryRegions[r]->eintervals[i].code)) {
					case Element::CODE::TRIANGLE3: ops[r][i].push_back(new Basis<Element::CODE::TRIANGLE3,  3, StructuralMechanicsGPC::TRIANGLE3, 2, StructuralMechanicsElementType::FACE, StructuralMechanicsDataDescriptor< 3, StructuralMechanicsGPC::TRIANGLE3, 3, 2, StructuralMechanicsElementType::FACE> >()); break;
					case Element::CODE::TRIANGLE6: ops[r][i].push_back(new Basis<Element::CODE::TRIANGLE6,  6, StructuralMechanicsGPC::TRIANGLE6, 2, StructuralMechanicsElementType::FACE, StructuralMechanicsDataDescriptor< 6, StructuralMechanicsGPC::TRIANGLE6, 3, 2, StructuralMechanicsElementType::FACE> >()); break;
					case Element::CODE::SQUARE4  : ops[r][i].push_back(new Basis<Element::CODE::SQUARE4  ,  4, StructuralMechanicsGPC::SQUARE4  , 2, StructuralMechanicsElementType::FACE, StructuralMechanicsDataDescriptor< 4, StructuralMechanicsGPC::SQUARE4  , 3, 2, StructuralMechanicsElementType::FACE> >()); break;
					case Element::CODE::SQUARE8  : ops[r][i].push_back(new Basis<Element::CODE::SQUARE8  ,  8, StructuralMechanicsGPC::SQUARE8  , 2, StructuralMechanicsElementType::FACE, StructuralMechanicsDataDescriptor< 8, StructuralMechanicsGPC::SQUARE8  , 3, 2, StructuralMechanicsElementType::FACE> >()); break;
					}
				break;
				case 1:
					switch (info::mesh->dimension) {
					case 2:
						if (axisymmetric) {
							switch (static_cast<Element::CODE>(info::mesh->boundaryRegions[r]->eintervals[i].code)) {
							case Element::CODE::LINE2: ops[r][i].push_back(new Basis<Element::CODE::LINE2, 2, StructuralMechanicsGPC::LINE2, 1, StructuralMechanicsElementType::EDGE_AXISYMMETRIC, StructuralMechanicsDataDescriptor< 2, StructuralMechanicsGPC::LINE2, 2, 1, StructuralMechanicsElementType::EDGE_AXISYMMETRIC> >()); break;
							case Element::CODE::LINE3: ops[r][i].push_back(new Basis<Element::CODE::LINE3, 3, StructuralMechanicsGPC::LINE3, 1, StructuralMechanicsElementType::EDGE_AXISYMMETRIC, StructuralMechanicsDataDescriptor< 3, StructuralMechanicsGPC::LINE3, 2, 1, StructuralMechanicsElementType::EDGE_AXISYMMETRIC> >()); break;
							}
						} else {
							switch (static_cast<Element::CODE>(info::mesh->boundaryRegions[r]->eintervals[i].code)) {
							case Element::CODE::LINE2: ops[r][i].push_back(new Basis<Element::CODE::LINE2, 2, StructuralMechanicsGPC::LINE2, 1, StructuralMechanicsElementType::EDGE, StructuralMechanicsDataDescriptor< 2, StructuralMechanicsGPC::LINE2, 2, 1, StructuralMechanicsElementType::EDGE> >()); break;
							case Element::CODE::LINE3: ops[r][i].push_back(new Basis<Element::CODE::LINE3, 3, StructuralMechanicsGPC::LINE3, 1, StructuralMechanicsElementType::EDGE, StructuralMechanicsDataDescriptor< 3, StructuralMechanicsGPC::LINE3, 2, 1, StructuralMechanicsElementType::EDGE> >()); break;
							}
						}
						break;
					case 3:
						switch (static_cast<Element::CODE>(info::mesh->boundaryRegions[r]->eintervals[i].code)) {
						case Element::CODE::LINE2: ops[r][i].push_back(new Basis<Element::CODE::LINE2, 2, StructuralMechanicsGPC::LINE2, 1, StructuralMechanicsElementType::EDGE, StructuralMechanicsDataDescriptor< 2, StructuralMechanicsGPC::LINE2, 3, 1, StructuralMechanicsElementType::EDGE> >()); break;
						case Element::CODE::LINE3: ops[r][i].push_back(new Basis<Element::CODE::LINE3, 3, StructuralMechanicsGPC::LINE3, 1, StructuralMechanicsElementType::EDGE, StructuralMechanicsDataDescriptor< 3, StructuralMechanicsGPC::LINE3, 3, 1, StructuralMechanicsElementType::EDGE> >()); break;
						} break;
					}
					break;
				}
			}
		}
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Filler, class Target, class Setter>
static ActionOperator* generateNodeSetter2D(size_t region, size_t interval, size_t dofs, Target *target, const Setter &setter)
{
	return new Filler< 1, StructuralMechanicsGPC::POINT1, 2, 0, StructuralMechanicsElementType::NODE, StructuralMechanicsDataDescriptor< 1, StructuralMechanicsGPC::POINT1, 2, 0, StructuralMechanicsElementType::NODE>, Setter>(region, interval, dofs, target, setter);
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Filler, class Target, class Setter>
static ActionOperator* generateNodeSetter3D(size_t region, size_t interval, size_t dofs, Target *target, const Setter &setter)
{
	return new Filler< 1, StructuralMechanicsGPC::POINT1, 3, 0, StructuralMechanicsElementType::NODE, StructuralMechanicsDataDescriptor< 1, StructuralMechanicsGPC::POINT1, 3, 0, StructuralMechanicsElementType::NODE>, Setter>(region, interval, dofs, target, setter);
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Filler, class Target, class Setter>
static ActionOperator* generateNodeSetter(size_t region, size_t interval, size_t dofs, Target *target, const Setter &setter)
{
	switch (info::mesh->dimension) {
	case 2: return generateNodeSetter2D<Filler>(region, interval, dofs, target, setter);
	case 3: return generateNodeSetter3D<Filler>(region, interval, dofs, target, setter);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Filler, class ... Args>
static ActionOperator* generateNodeFiller2D(size_t region, size_t interval, size_t dofs, Args&& ... args)
{
	return new Filler< 1, StructuralMechanicsGPC::POINT1, 2, 0, StructuralMechanicsElementType::NODE, StructuralMechanicsDataDescriptor< 1, StructuralMechanicsGPC::POINT1, 2, 0, StructuralMechanicsElementType::NODE> >(region, interval, dofs, std::forward<Args>(args)...);
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Filler, class ... Args>
static ActionOperator* generateNodeFiller3D(size_t region, size_t interval, size_t dofs, Args&& ... args)
{
	return new Filler< 1, StructuralMechanicsGPC::POINT1, 3, 0, StructuralMechanicsElementType::NODE, StructuralMechanicsDataDescriptor< 1, StructuralMechanicsGPC::POINT1, 3, 0, StructuralMechanicsElementType::NODE> >(region, interval, dofs, std::forward<Args>(args)...);
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Filler, class ... Args>
static ActionOperator* generateNodeFiller(size_t region, size_t interval, size_t dofs, Args&& ... args)
{
	switch (info::mesh->dimension) {
	case 2: return generateNodeFiller2D<Filler>(region, interval, dofs, std::forward<Args>(args)...);
	case 3: return generateNodeFiller3D<Filler>(region, interval, dofs, std::forward<Args>(args)...);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Filler, class ... Args>
static ActionOperator* generateEdgeFiller2D(size_t region, size_t interval, size_t dofs, Args&& ... args)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::LINE2): return new Filler< 2, StructuralMechanicsGPC::LINE2, 2, 1, StructuralMechanicsElementType::EDGE, StructuralMechanicsDataDescriptor< 2, StructuralMechanicsGPC::LINE2, 2, 1, StructuralMechanicsElementType::EDGE> >(region, interval, dofs, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::LINE3): return new Filler< 3, StructuralMechanicsGPC::LINE3, 2, 1, StructuralMechanicsElementType::EDGE, StructuralMechanicsDataDescriptor< 3, StructuralMechanicsGPC::LINE3, 2, 1, StructuralMechanicsElementType::EDGE> >(region, interval, dofs, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Filler, class ... Args>
static ActionOperator* generateEdgeFiller2DAxisymmetric(size_t region, size_t interval, size_t dofs, Args&& ... args)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::LINE2): return new Filler< 2, StructuralMechanicsGPC::LINE2, 2, 1, StructuralMechanicsElementType::EDGE_AXISYMMETRIC, StructuralMechanicsDataDescriptor< 2, StructuralMechanicsGPC::LINE2, 2, 1, StructuralMechanicsElementType::EDGE_AXISYMMETRIC> >(region, interval, dofs, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::LINE3): return new Filler< 3, StructuralMechanicsGPC::LINE3, 2, 1, StructuralMechanicsElementType::EDGE_AXISYMMETRIC, StructuralMechanicsDataDescriptor< 3, StructuralMechanicsGPC::LINE3, 2, 1, StructuralMechanicsElementType::EDGE_AXISYMMETRIC> >(region, interval, dofs, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Filler, class ... Args>
static ActionOperator* generateEdgeFiller3D(size_t region, size_t interval, size_t dofs, Args&& ... args)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::LINE2): return new Filler< 2, StructuralMechanicsGPC::LINE2, 3, 1, StructuralMechanicsElementType::EDGE, StructuralMechanicsDataDescriptor< 2, StructuralMechanicsGPC::LINE2, 3, 1, StructuralMechanicsElementType::EDGE> >(region, interval, dofs, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::LINE3): return new Filler< 3, StructuralMechanicsGPC::LINE3, 3, 1, StructuralMechanicsElementType::EDGE, StructuralMechanicsDataDescriptor< 3, StructuralMechanicsGPC::LINE3, 3, 1, StructuralMechanicsElementType::EDGE> >(region, interval, dofs, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Filler, class ... Args>
static ActionOperator* generateEdgeFiller(int axisymmetric, size_t region, size_t interval, size_t dofs, Args&& ... args)
{
	if (axisymmetric) {
		return generateEdgeFiller2DAxisymmetric<Filler>(region, interval, dofs, std::forward<Args>(args)...);
	}
	switch (info::mesh->dimension) {
	case 2: return generateEdgeFiller2D<Filler>(region, interval, dofs, std::forward<Args>(args)...);
	case 3: return generateEdgeFiller3D<Filler>(region, interval, dofs, std::forward<Args>(args)...);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Filler, class ... Args>
static ActionOperator* generateFaceFiller(size_t region, size_t interval, size_t dofs, Args&& ... args)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::TRIANGLE3): return new Filler< 3, StructuralMechanicsGPC::TRIANGLE3, 3, 2, StructuralMechanicsElementType::FACE, StructuralMechanicsDataDescriptor< 3, StructuralMechanicsGPC::TRIANGLE3, 3, 2, StructuralMechanicsElementType::FACE> >(region, interval, dofs, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::TRIANGLE6): return new Filler< 6, StructuralMechanicsGPC::TRIANGLE6, 3, 2, StructuralMechanicsElementType::FACE, StructuralMechanicsDataDescriptor< 6, StructuralMechanicsGPC::TRIANGLE6, 3, 2, StructuralMechanicsElementType::FACE> >(region, interval, dofs, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::SQUARE4  ): return new Filler< 4, StructuralMechanicsGPC::SQUARE4  , 3, 2, StructuralMechanicsElementType::FACE, StructuralMechanicsDataDescriptor< 4, StructuralMechanicsGPC::SQUARE4  , 3, 2, StructuralMechanicsElementType::FACE> >(region, interval, dofs, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::SQUARE8  ): return new Filler< 8, StructuralMechanicsGPC::SQUARE8  , 3, 2, StructuralMechanicsElementType::FACE, StructuralMechanicsDataDescriptor< 8, StructuralMechanicsGPC::SQUARE8  , 3, 2, StructuralMechanicsElementType::FACE> >(region, interval, dofs, std::forward<Args>(args)...);
	default: return nullptr;
	}
}


template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateTypedExpressionNode2D(size_t region, size_t interval, Evaluator *evaluator, const Setter &setter)
{
	return new Expression< 1, StructuralMechanicsGPC::POINT1, 2, 0, StructuralMechanicsElementType::NODE, StructuralMechanicsDataDescriptor< 1, StructuralMechanicsGPC::POINT1, 2, 0, StructuralMechanicsElementType::NODE>, Setter>(interval, evaluator, setter);
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateTypedExpressionEdge2D(size_t region, size_t interval, Evaluator *evaluator, const Setter &setter)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::LINE2): return new Expression< 2, StructuralMechanicsGPC::LINE2, 2, 1, StructuralMechanicsElementType::EDGE, StructuralMechanicsDataDescriptor< 2, StructuralMechanicsGPC::LINE2, 2, 1, StructuralMechanicsElementType::EDGE>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::LINE3): return new Expression< 3, StructuralMechanicsGPC::LINE3, 2, 1, StructuralMechanicsElementType::EDGE, StructuralMechanicsDataDescriptor< 3, StructuralMechanicsGPC::LINE3, 2, 1, StructuralMechanicsElementType::EDGE>, Setter>(interval, evaluator, setter);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateTypedExpressionEdge2DAxisymmetric(size_t region, size_t interval, Evaluator *evaluator, const Setter &setter)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::LINE2): return new Expression< 2, StructuralMechanicsGPC::LINE2, 2, 1, StructuralMechanicsElementType::EDGE_AXISYMMETRIC, StructuralMechanicsDataDescriptor< 2, StructuralMechanicsGPC::LINE2, 2, 1, StructuralMechanicsElementType::EDGE_AXISYMMETRIC>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::LINE3): return new Expression< 3, StructuralMechanicsGPC::LINE3, 2, 1, StructuralMechanicsElementType::EDGE_AXISYMMETRIC, StructuralMechanicsDataDescriptor< 3, StructuralMechanicsGPC::LINE3, 2, 1, StructuralMechanicsElementType::EDGE_AXISYMMETRIC>, Setter>(interval, evaluator, setter);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateTypedExpressionNode3D(size_t region, size_t interval, Evaluator *evaluator, const Setter &setter)
{
	return new Expression< 1, StructuralMechanicsGPC::POINT1, 3, 0, StructuralMechanicsElementType::NODE, StructuralMechanicsDataDescriptor< 1, StructuralMechanicsGPC::POINT1, 3, 0, StructuralMechanicsElementType::NODE>, Setter>(interval, evaluator, setter);
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateTypedExpressionEdge3D(size_t region, size_t interval, Evaluator *evaluator, const Setter &setter)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::LINE2): return new Expression< 2, StructuralMechanicsGPC::LINE2, 3, 1, StructuralMechanicsElementType::EDGE, StructuralMechanicsDataDescriptor< 2, StructuralMechanicsGPC::LINE2, 3, 1, StructuralMechanicsElementType::EDGE>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::LINE3): return new Expression< 3, StructuralMechanicsGPC::LINE3, 3, 1, StructuralMechanicsElementType::EDGE, StructuralMechanicsDataDescriptor< 3, StructuralMechanicsGPC::LINE3, 3, 1, StructuralMechanicsElementType::EDGE>, Setter>(interval, evaluator, setter);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateTypedExpressionFace3D(size_t region, size_t interval, Evaluator *evaluator, const Setter &setter)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::TRIANGLE3): return new Expression< 3, StructuralMechanicsGPC::TRIANGLE3, 3, 2, StructuralMechanicsElementType::FACE, StructuralMechanicsDataDescriptor< 3, StructuralMechanicsGPC::TRIANGLE3, 3, 2, StructuralMechanicsElementType::FACE>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::TRIANGLE6): return new Expression< 6, StructuralMechanicsGPC::TRIANGLE6, 3, 2, StructuralMechanicsElementType::FACE, StructuralMechanicsDataDescriptor< 6, StructuralMechanicsGPC::TRIANGLE6, 3, 2, StructuralMechanicsElementType::FACE>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::SQUARE4  ): return new Expression< 4, StructuralMechanicsGPC::SQUARE4  , 3, 2, StructuralMechanicsElementType::FACE, StructuralMechanicsDataDescriptor< 4, StructuralMechanicsGPC::SQUARE4  , 3, 2, StructuralMechanicsElementType::FACE>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::SQUARE8  ): return new Expression< 8, StructuralMechanicsGPC::SQUARE8  , 3, 2, StructuralMechanicsElementType::FACE, StructuralMechanicsDataDescriptor< 8, StructuralMechanicsGPC::SQUARE8  , 3, 2, StructuralMechanicsElementType::FACE>, Setter>(interval, evaluator, setter);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, int etype, class Setter>
static ActionOperator* generateTypedExpression2D(size_t interval, Evaluator *evaluator, const Setter &setter)
{
	switch (info::mesh->elements->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::TRIANGLE3): return new Expression< 3, StructuralMechanicsGPC::TRIANGLE3, 2, 2, etype, StructuralMechanicsDataDescriptor< 3, StructuralMechanicsGPC::TRIANGLE3, 2, 2, etype>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::TRIANGLE6): return new Expression< 6, StructuralMechanicsGPC::TRIANGLE6, 2, 2, etype, StructuralMechanicsDataDescriptor< 6, StructuralMechanicsGPC::TRIANGLE6, 2, 2, etype>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::SQUARE4  ): return new Expression< 4, StructuralMechanicsGPC::SQUARE4  , 2, 2, etype, StructuralMechanicsDataDescriptor< 4, StructuralMechanicsGPC::SQUARE4  , 2, 2, etype>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::SQUARE8  ): return new Expression< 8, StructuralMechanicsGPC::SQUARE8  , 2, 2, etype, StructuralMechanicsDataDescriptor< 8, StructuralMechanicsGPC::SQUARE8  , 2, 2, etype>, Setter>(interval, evaluator, setter);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, int etype, class Setter>
static ActionOperator* generateTypedExpression3D(size_t interval, Evaluator *evaluator, const Setter &setter)
{
	switch (info::mesh->elements->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::TETRA4   ): return new Expression< 4, StructuralMechanicsGPC::TETRA4   , 3, 3, etype, StructuralMechanicsDataDescriptor< 4, StructuralMechanicsGPC::TETRA4   , 3, 3, etype>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::TETRA10  ): return new Expression<10, StructuralMechanicsGPC::TETRA10  , 3, 3, etype, StructuralMechanicsDataDescriptor<10, StructuralMechanicsGPC::TETRA10  , 3, 3, etype>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::PYRAMID5 ): return new Expression< 5, StructuralMechanicsGPC::PYRAMID5 , 3, 3, etype, StructuralMechanicsDataDescriptor< 5, StructuralMechanicsGPC::PYRAMID5 , 3, 3, etype>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::PYRAMID13): return new Expression<13, StructuralMechanicsGPC::PYRAMID13, 3, 3, etype, StructuralMechanicsDataDescriptor<13, StructuralMechanicsGPC::PYRAMID13, 3, 3, etype>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::PRISMA6  ): return new Expression< 6, StructuralMechanicsGPC::PRISMA6  , 3, 3, etype, StructuralMechanicsDataDescriptor< 6, StructuralMechanicsGPC::PRISMA6  , 3, 3, etype>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::PRISMA15 ): return new Expression<15, StructuralMechanicsGPC::PRISMA15 , 3, 3, etype, StructuralMechanicsDataDescriptor<15, StructuralMechanicsGPC::PRISMA15 , 3, 3, etype>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::HEXA8    ): return new Expression< 8, StructuralMechanicsGPC::HEXA8    , 3, 3, etype, StructuralMechanicsDataDescriptor< 8, StructuralMechanicsGPC::HEXA8    , 3, 3, etype>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::HEXA20   ): return new Expression<20, StructuralMechanicsGPC::HEXA20   , 3, 3, etype, StructuralMechanicsDataDescriptor<20, StructuralMechanicsGPC::HEXA20   , 3, 3, etype>, Setter>(interval, evaluator, setter);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, int etype, class Setter>
static ActionOperator* generateTypedExpression(size_t interval, Evaluator *evaluator, const Setter &setter)
{
	switch (info::mesh->dimension) {
	case 2: return generateTypedExpression2D<Expression, etype, Setter>(interval, evaluator, setter);
	case 3: return generateTypedExpression3D<Expression, etype, Setter>(interval, evaluator, setter);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateExpression(size_t interval, int etype, Evaluator *evaluator, const Setter &setter)
{
	switch (info::mesh->dimension) {
	case 2:
		switch (etype) {
		case StructuralMechanicsElementType::SYMMETRIC_PLANE             : return generateTypedExpression2D<Expression, StructuralMechanicsElementType::SYMMETRIC_PLANE             >(interval, evaluator, setter);
		case StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC: return generateTypedExpression2D<Expression, StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC>(interval, evaluator, setter);
		}
		break;
	case 3:
		switch (etype) {
		case StructuralMechanicsElementType::SYMMETRIC_VOLUME            : return generateTypedExpression3D<Expression, StructuralMechanicsElementType::SYMMETRIC_VOLUME            >(interval, evaluator, setter);
		}
		break;
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateExpression2D(size_t interval, int etype, Evaluator *evaluator, const Setter &setter)
{
	switch (etype) {
	case StructuralMechanicsElementType::SYMMETRIC_PLANE             : return generateTypedExpression2D<Expression, StructuralMechanicsElementType::SYMMETRIC_PLANE             >(interval, evaluator, setter);
	case StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC: return generateTypedExpression2D<Expression, StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC>(interval, evaluator, setter);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateExpression3D(size_t interval, int etype, Evaluator *evaluator, const Setter &setter)
{
	switch (etype) {
	case StructuralMechanicsElementType::SYMMETRIC_VOLUME            : return generateTypedExpression3D<Expression, StructuralMechanicsElementType::SYMMETRIC_VOLUME            >(interval, evaluator, setter);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static void generateElementExpression2D(const std::vector<int> &etype, std::vector<std::vector<ActionOperator*> > &ops, std::map<std::string, ECFExpression> &settings, const Setter &setter)
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (Assembler::getEvaluator(interval, settings)) {
			ops[interval].push_back(generateExpression2D<Expression>(interval, etype[interval], Assembler::getEvaluator(interval, settings), setter));
		}
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static void generateElementExpression3D(const std::vector<int> &etype, std::vector<std::vector<ActionOperator*> > &ops, std::map<std::string, ECFExpression> &settings, const Setter &setter)
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (Assembler::getEvaluator(interval, settings)) {
			ops[interval].push_back(generateExpression3D<Expression>(interval, etype[interval], Assembler::getEvaluator(interval, settings), setter));
		}
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static void generateElementExpression(const std::vector<int> &etype, std::vector<std::vector<ActionOperator*> > &ops, std::map<std::string, ECFExpression> &settings, const Setter &setter)
{
	switch (info::mesh->dimension) {
	case 2: generateElementExpression2D<Expression>(etype, ops, settings, setter); break;
	case 3: generateElementExpression3D<Expression>(etype, ops, settings, setter); break;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static void generateElementExpression2D(const std::vector<int> &etype, std::vector<std::vector<ActionOperator*> > &ops, std::map<std::string, ECFExpressionVector> &settings, int dim, const Setter &setter)
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (Assembler::getEvaluator(interval, settings, dim)) {
			ops[interval].push_back(generateExpression2D<Expression>(interval, etype[interval], Assembler::getEvaluator(interval, settings, dim), setter));
		}
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static void generateElementExpression3D(const std::vector<int> &etype, std::vector<std::vector<ActionOperator*> > &ops, std::map<std::string, ECFExpressionVector> &settings, int dim, const Setter &setter)
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (Assembler::getEvaluator(interval, settings, dim)) {
			ops[interval].push_back(generateExpression3D<Expression>(interval, etype[interval], Assembler::getEvaluator(interval, settings, dim), setter));
		}
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static void generateElementExpression(const std::vector<int> &etype, std::vector<std::vector<ActionOperator*> > &ops, std::map<std::string, ECFExpressionVector> &settings, int dim, const Setter &setter)
{
	switch (info::mesh->dimension) {
	case 2: generateElementExpression2D<Expression>(etype, ops, settings, dim, setter); break;
	case 3: generateElementExpression3D<Expression>(etype, ops, settings, dim, setter); break;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static void generateBoundaryExpression(int axisymmetric, size_t region, std::vector<std::vector<std::vector<ActionOperator*> > > &ops, Evaluator* evaluator, const Setter &setter)
{
	switch (info::mesh->dimension) {
	case 2:
		switch (info::mesh->boundaryRegions[region]->dimension) {
		case 0:
			for (size_t t = 0; t < info::mesh->boundaryRegions[region]->nodes->threads(); ++t) {
				ops[region][t].push_back(generateTypedExpressionNode2D<Expression>(region, t, evaluator, setter));
			}
			break;
		case 1:
			if (axisymmetric) {
				for (size_t interval = 0; interval < info::mesh->boundaryRegions[region]->eintervals.size(); ++interval) {
					ops[region][interval].push_back(generateTypedExpressionEdge2DAxisymmetric<Expression>(region, interval, evaluator, setter));
				}
			} else {
				for (size_t interval = 0; interval < info::mesh->boundaryRegions[region]->eintervals.size(); ++interval) {
					ops[region][interval].push_back(generateTypedExpressionEdge2D<Expression>(region, interval, evaluator, setter));
				}
			}
			break;
		}
		break;
	case 3:
		switch (info::mesh->boundaryRegions[region]->dimension) {
		case 0:
			for (size_t t = 0; t < info::mesh->boundaryRegions[region]->nodes->threads(); ++t) {
				ops[region][t].push_back(generateTypedExpressionNode3D<Expression>(region, t, evaluator, setter));
			}
			break;
		case 1:
			for (size_t interval = 0; interval < info::mesh->boundaryRegions[region]->eintervals.size(); ++interval) {
				ops[region][interval].push_back(generateTypedExpressionEdge3D<Expression>(region, interval, evaluator, setter));
			}
			break;
		case 2:
			for (size_t interval = 0; interval < info::mesh->boundaryRegions[region]->eintervals.size(); ++interval) {
				ops[region][interval].push_back(generateTypedExpressionFace3D<Expression>(region, interval, evaluator, setter));
			}
			break;
		}
		break;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static void generateBoundaryExpression(int axisymmetric, std::vector<std::vector<std::vector<ActionOperator*> > > &ops, std::map<std::string, ECFExpression> &settings, const Setter &setter)
{
	for(size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		auto it = settings.find(info::mesh->boundaryRegions[r]->name);
		if (it != settings.end()) {
			generateBoundaryExpression<Expression>(axisymmetric, r, ops, it->second.evaluator, setter);
		}
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static void generateBoundaryExpression(int axisymmetric, std::vector<std::vector<std::vector<ActionOperator*> > > &ops, std::map<std::string, ECFExpressionVector> &settings, int dim, const Setter &setter)
{
	for(size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		auto it = settings.find(info::mesh->boundaryRegions[r]->name);
		if (it != settings.end()) {
			generateBoundaryExpression<Expression>(axisymmetric, r, ops, it->second.data[dim].evaluator, setter);
		}
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static void generateBoundaryExpression(int axisymmetric, std::vector<std::vector<std::vector<ActionOperator*> > > &ops, std::map<std::string, ECFExpressionOptionalVector> &settings, int dim, const Setter &setter)
{
	for(size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		auto it = settings.find(info::mesh->boundaryRegions[r]->name);
		if (it != settings.end() && it->second.data[dim].isset) {
			generateBoundaryExpression<Expression>(axisymmetric, r, ops, it->second.data[dim].evaluator, setter);
		}
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, int etype, class ... Args>
static ActionOperator* generateElementTypedOperator2D(size_t interval, Args&& ... args)
{
	switch (info::mesh->elements->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::TRIANGLE3): return new Operator< 3, StructuralMechanicsGPC::TRIANGLE3, 2, 2, etype, StructuralMechanicsDataDescriptor< 3, StructuralMechanicsGPC::TRIANGLE3, 2, 2, etype> >(interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::TRIANGLE6): return new Operator< 6, StructuralMechanicsGPC::TRIANGLE6, 2, 2, etype, StructuralMechanicsDataDescriptor< 6, StructuralMechanicsGPC::TRIANGLE6, 2, 2, etype> >(interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::SQUARE4  ): return new Operator< 4, StructuralMechanicsGPC::SQUARE4  , 2, 2, etype, StructuralMechanicsDataDescriptor< 4, StructuralMechanicsGPC::SQUARE4  , 2, 2, etype> >(interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::SQUARE8  ): return new Operator< 8, StructuralMechanicsGPC::SQUARE8  , 2, 2, etype, StructuralMechanicsDataDescriptor< 8, StructuralMechanicsGPC::SQUARE8  , 2, 2, etype> >(interval, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, int etype, class ... Args>
static ActionOperator* generateElementTypedOperator3D(size_t interval, Args&& ... args)
{
	switch (info::mesh->elements->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::TETRA4   ): return new Operator< 4, StructuralMechanicsGPC::TETRA4   , 3, 3, etype, StructuralMechanicsDataDescriptor< 4, StructuralMechanicsGPC::TETRA4   , 3, 3, etype> >(interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::TETRA10  ): return new Operator<10, StructuralMechanicsGPC::TETRA10  , 3, 3, etype, StructuralMechanicsDataDescriptor<10, StructuralMechanicsGPC::TETRA10  , 3, 3, etype> >(interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::PYRAMID5 ): return new Operator< 5, StructuralMechanicsGPC::PYRAMID5 , 3, 3, etype, StructuralMechanicsDataDescriptor< 5, StructuralMechanicsGPC::PYRAMID5 , 3, 3, etype> >(interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::PYRAMID13): return new Operator<13, StructuralMechanicsGPC::PYRAMID13, 3, 3, etype, StructuralMechanicsDataDescriptor<13, StructuralMechanicsGPC::PYRAMID13, 3, 3, etype> >(interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::PRISMA6  ): return new Operator< 6, StructuralMechanicsGPC::PRISMA6  , 3, 3, etype, StructuralMechanicsDataDescriptor< 6, StructuralMechanicsGPC::PRISMA6  , 3, 3, etype> >(interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::PRISMA15 ): return new Operator<15, StructuralMechanicsGPC::PRISMA15 , 3, 3, etype, StructuralMechanicsDataDescriptor<15, StructuralMechanicsGPC::PRISMA15 , 3, 3, etype> >(interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::HEXA8    ): return new Operator< 8, StructuralMechanicsGPC::HEXA8    , 3, 3, etype, StructuralMechanicsDataDescriptor< 8, StructuralMechanicsGPC::HEXA8    , 3, 3, etype> >(interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::HEXA20   ): return new Operator<20, StructuralMechanicsGPC::HEXA20   , 3, 3, etype, StructuralMechanicsDataDescriptor<20, StructuralMechanicsGPC::HEXA20   , 3, 3, etype> >(interval, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, int etype, class ... Args>
static ActionOperator* generateElementTypedOperator(size_t interval, Args&& ... args)
{
	switch (info::mesh->dimension) {
	case 2: return generateElementTypedOperator2D<Operator, etype, Args...>(interval, std::forward<Args>(args)...);
	case 3: return generateElementTypedOperator3D<Operator, etype, Args...>(interval, std::forward<Args>(args)...);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static ActionOperator* generateElementOperator2D(size_t interval, int etype, Args&& ... args)
{
	switch (etype) {
	case StructuralMechanicsElementType::SYMMETRIC_PLANE             : return generateElementTypedOperator2D<Operator, StructuralMechanicsElementType::SYMMETRIC_PLANE             , Args...>(interval, std::forward<Args>(args)...);
	case StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC: return generateElementTypedOperator2D<Operator, StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC, Args...>(interval, std::forward<Args>(args)...);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static ActionOperator* generateElementOperator3D(size_t interval, int etype, Args&& ... args)
{
	switch (etype) {
	case StructuralMechanicsElementType::SYMMETRIC_VOLUME            : return generateElementTypedOperator3D<Operator, StructuralMechanicsElementType::SYMMETRIC_VOLUME            , Args...>(interval, std::forward<Args>(args)...);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static ActionOperator* generateElementOperator(size_t interval, int etype, Args&& ... args)
{
	switch (info::mesh->dimension) {
	case 2: return generateElementOperator2D<Operator>(interval, etype, std::forward<Args>(args)...);
	case 3: return generateElementOperator3D<Operator>(interval, etype, std::forward<Args>(args)...);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static void generateElementOperators(const std::vector<int> &etype, std::vector<std::vector<ActionOperator*> > &ops, Args&& ... args)
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		ops[interval].push_back(generateElementOperator<Operator, Args...>(interval, etype[interval], std::forward<Args>(args)...));
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static ActionOperator* generateBoundaryEdge2DOperator(size_t region, size_t interval, Args&& ... args)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::LINE2): return new Operator< 2, StructuralMechanicsGPC::LINE2, 2, 1, StructuralMechanicsElementType::EDGE, StructuralMechanicsDataDescriptor< 2, StructuralMechanicsGPC::LINE2, 2, 1, StructuralMechanicsElementType::EDGE> >(region, interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::LINE3): return new Operator< 3, StructuralMechanicsGPC::LINE3, 2, 1, StructuralMechanicsElementType::EDGE, StructuralMechanicsDataDescriptor< 3, StructuralMechanicsGPC::LINE3, 2, 1, StructuralMechanicsElementType::EDGE> >(region, interval, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static ActionOperator* generateBoundaryEdge2DAxisymmetricOperator(size_t region, size_t interval, Args&& ... args)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::LINE2): return new Operator< 2, StructuralMechanicsGPC::LINE2, 2, 1, StructuralMechanicsElementType::EDGE_AXISYMMETRIC, StructuralMechanicsDataDescriptor< 2, StructuralMechanicsGPC::LINE2, 2, 1, StructuralMechanicsElementType::EDGE_AXISYMMETRIC> >(region, interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::LINE3): return new Operator< 3, StructuralMechanicsGPC::LINE3, 2, 1, StructuralMechanicsElementType::EDGE_AXISYMMETRIC, StructuralMechanicsDataDescriptor< 3, StructuralMechanicsGPC::LINE3, 2, 1, StructuralMechanicsElementType::EDGE_AXISYMMETRIC> >(region, interval, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static ActionOperator* generateBoundaryEdge3DOperator(size_t region, size_t interval, Args&& ... args)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::LINE2): return new Operator< 2, StructuralMechanicsGPC::LINE2, 3, 1, StructuralMechanicsElementType::EDGE, StructuralMechanicsDataDescriptor< 2, StructuralMechanicsGPC::LINE2, 3, 1, StructuralMechanicsElementType::EDGE> >(region, interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::LINE3): return new Operator< 3, StructuralMechanicsGPC::LINE3, 3, 1, StructuralMechanicsElementType::EDGE, StructuralMechanicsDataDescriptor< 3, StructuralMechanicsGPC::LINE3, 3, 1, StructuralMechanicsElementType::EDGE> >(region, interval, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static ActionOperator* generateBoundaryFaceOperator(size_t region, size_t interval, Args&& ... args)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::TRIANGLE3): return new Operator< 3, StructuralMechanicsGPC::TRIANGLE3, 3, 2, StructuralMechanicsElementType::FACE, StructuralMechanicsDataDescriptor< 3, StructuralMechanicsGPC::TRIANGLE3, 3, 2, StructuralMechanicsElementType::FACE> >(region, interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::TRIANGLE6): return new Operator< 6, StructuralMechanicsGPC::TRIANGLE6, 3, 2, StructuralMechanicsElementType::FACE, StructuralMechanicsDataDescriptor< 6, StructuralMechanicsGPC::TRIANGLE6, 3, 2, StructuralMechanicsElementType::FACE> >(region, interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::SQUARE4  ): return new Operator< 4, StructuralMechanicsGPC::SQUARE4  , 3, 2, StructuralMechanicsElementType::FACE, StructuralMechanicsDataDescriptor< 4, StructuralMechanicsGPC::SQUARE4  , 3, 2, StructuralMechanicsElementType::FACE> >(region, interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::SQUARE8  ): return new Operator< 8, StructuralMechanicsGPC::SQUARE8  , 3, 2, StructuralMechanicsElementType::FACE, StructuralMechanicsDataDescriptor< 8, StructuralMechanicsGPC::SQUARE8  , 3, 2, StructuralMechanicsElementType::FACE> >(region, interval, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static ActionOperator* generateBoundaryOperator(int axisymmetric, size_t region, size_t interval, Args&& ... args)
{
	if (axisymmetric) {
		return generateBoundaryEdge2DAxisymmetricOperator<Operator>(region, interval, std::forward<Args>(args)...);
	}
	if (info::mesh->boundaryRegions[region]->dimension == 1) {
		switch (info::mesh->dimension) {
		case 2: return generateBoundaryEdge2DOperator<Operator>(region, interval, std::forward<Args>(args)...);
		case 3: return generateBoundaryEdge3DOperator<Operator>(region, interval, std::forward<Args>(args)...);
		}
	}
	if (info::mesh->boundaryRegions[region]->dimension == 2) {
		return generateBoundaryFaceOperator<Operator>(region, interval, std::forward<Args>(args)...);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static void generateBoundaryOperators(int axisymmetric, const std::vector<int> &bfilter, std::vector<std::vector<std::vector<ActionOperator*> > > &ops, Args&& ... args)
{
	for(size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (bfilter[r]) {
			for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
				ops[r][i].push_back(generateBoundaryOperator<Operator>(axisymmetric, r, i, std::forward<Args>(args)...));
			}
		}
	}
}

static void dropLastOperators(std::vector<std::vector<ActionOperator*> > &ops)
{
	for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		delete ops[i].back();
		ops[i].pop_back();
	}
}

static void dropLastOperators(const std::vector<int> &bfilter, std::vector<std::vector<std::vector<ActionOperator*> > > &ops)
{
	for(size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (bfilter[r]) {
			for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
				delete ops[r][i].back();
				ops[r][i].pop_back();
			}
		}
	}
}

}



#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_STRUCTURALMECHANICS_GENERATOR_H_ */
