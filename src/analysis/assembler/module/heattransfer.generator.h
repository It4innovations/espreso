
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_GENERATOR_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_GENERATOR_H_

#include "heattransfer.element.h"
#include "analysis/assembler/operator.h"
#include "analysis/assembler/operators/basis.h"
#include "analysis/assembler/operators/expression.h"
#include "analysis/assembler/operators/storage.h"

#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"

namespace espreso {

template <size_t ndim, int etype>
static void generateBaseFunctions(size_t interval, std::vector<std::vector<ActionOperator*> > &ops)
{
	switch (static_cast<Element::CODE>(info::mesh->elements->eintervals[interval].code)) {
	case Element::CODE::LINE2    : ops[interval].push_back(new Basis<Element::CODE::LINE2    ,  2, HeatTransferGPC::LINE2    , 1, etype, HeatTransferDataDescriptor< 2, HeatTransferGPC::LINE2    , ndim, 1, etype> >()); break;
	case Element::CODE::LINE3    : ops[interval].push_back(new Basis<Element::CODE::LINE3    ,  3, HeatTransferGPC::LINE3    , 1, etype, HeatTransferDataDescriptor< 3, HeatTransferGPC::LINE3    , ndim, 1, etype> >()); break;
	case Element::CODE::TRIANGLE3: ops[interval].push_back(new Basis<Element::CODE::TRIANGLE3,  3, HeatTransferGPC::TRIANGLE3, 2, etype, HeatTransferDataDescriptor< 3, HeatTransferGPC::TRIANGLE3, ndim, 2, etype> >()); break;
	case Element::CODE::TRIANGLE6: ops[interval].push_back(new Basis<Element::CODE::TRIANGLE6,  6, HeatTransferGPC::TRIANGLE6, 2, etype, HeatTransferDataDescriptor< 6, HeatTransferGPC::TRIANGLE6, ndim, 2, etype> >()); break;
	case Element::CODE::SQUARE4  : ops[interval].push_back(new Basis<Element::CODE::SQUARE4  ,  4, HeatTransferGPC::SQUARE4  , 2, etype, HeatTransferDataDescriptor< 4, HeatTransferGPC::SQUARE4  , ndim, 2, etype> >()); break;
	case Element::CODE::SQUARE8  : ops[interval].push_back(new Basis<Element::CODE::SQUARE8  ,  8, HeatTransferGPC::SQUARE8  , 2, etype, HeatTransferDataDescriptor< 8, HeatTransferGPC::SQUARE8  , ndim, 2, etype> >()); break;
	case Element::CODE::TETRA4   : ops[interval].push_back(new Basis<Element::CODE::TETRA4   ,  4, HeatTransferGPC::TETRA4   , 3, etype, HeatTransferDataDescriptor< 4, HeatTransferGPC::TETRA4   , ndim, 3, etype> >()); break;
	case Element::CODE::TETRA10  : ops[interval].push_back(new Basis<Element::CODE::TETRA10  , 10, HeatTransferGPC::TETRA10  , 3, etype, HeatTransferDataDescriptor<10, HeatTransferGPC::TETRA10  , ndim, 3, etype> >()); break;
	case Element::CODE::PYRAMID5 : ops[interval].push_back(new Basis<Element::CODE::PYRAMID5 ,  5, HeatTransferGPC::PYRAMID5 , 3, etype, HeatTransferDataDescriptor< 5, HeatTransferGPC::PYRAMID5 , ndim, 3, etype> >()); break;
	case Element::CODE::PYRAMID13: ops[interval].push_back(new Basis<Element::CODE::PYRAMID13, 13, HeatTransferGPC::PYRAMID13, 3, etype, HeatTransferDataDescriptor<13, HeatTransferGPC::PYRAMID13, ndim, 3, etype> >()); break;
	case Element::CODE::PRISMA6  : ops[interval].push_back(new Basis<Element::CODE::PRISMA6  ,  6, HeatTransferGPC::PRISMA6  , 3, etype, HeatTransferDataDescriptor< 6, HeatTransferGPC::PRISMA6  , ndim, 3, etype> >()); break;
	case Element::CODE::PRISMA15 : ops[interval].push_back(new Basis<Element::CODE::PRISMA15 , 15, HeatTransferGPC::PRISMA15 , 3, etype, HeatTransferDataDescriptor<15, HeatTransferGPC::PRISMA15 , ndim, 3, etype> >()); break;
	case Element::CODE::HEXA8    : ops[interval].push_back(new Basis<Element::CODE::HEXA8    ,  8, HeatTransferGPC::HEXA8    , 3, etype, HeatTransferDataDescriptor< 8, HeatTransferGPC::HEXA8    , ndim, 3, etype> >()); break;
	case Element::CODE::HEXA20   : ops[interval].push_back(new Basis<Element::CODE::HEXA20   , 20, HeatTransferGPC::HEXA20   , 3, etype, HeatTransferDataDescriptor<20, HeatTransferGPC::HEXA20   , ndim, 3, etype> >()); break;
	}
}

static void generateBaseFunctions(const std::vector<int> &etype, std::vector<std::vector<ActionOperator*> > &ops)
{
	GaussPoints<Element::CODE::LINE2    ,  2, HeatTransferGPC::LINE2    , 1>::set();
	GaussPoints<Element::CODE::TRIANGLE3,  3, HeatTransferGPC::TRIANGLE3, 2>::set();
	GaussPoints<Element::CODE::SQUARE4  ,  4, HeatTransferGPC::SQUARE4  , 2>::set();
	GaussPoints<Element::CODE::TETRA4   ,  4, HeatTransferGPC::TETRA4   , 3>::set();
	GaussPoints<Element::CODE::PYRAMID5 ,  5, HeatTransferGPC::PYRAMID5 , 3>::set();
	GaussPoints<Element::CODE::PRISMA6  ,  6, HeatTransferGPC::PRISMA6  , 3>::set();
	GaussPoints<Element::CODE::HEXA8    ,  8, HeatTransferGPC::HEXA8    , 3>::set();
	GaussPoints<Element::CODE::LINE3    ,  3, HeatTransferGPC::LINE3    , 1>::set();
	GaussPoints<Element::CODE::TRIANGLE6,  6, HeatTransferGPC::TRIANGLE6, 2>::set();
	GaussPoints<Element::CODE::SQUARE8  ,  8, HeatTransferGPC::SQUARE8  , 2>::set();
	GaussPoints<Element::CODE::TETRA10  , 10, HeatTransferGPC::TETRA10  , 3>::set();
	GaussPoints<Element::CODE::PYRAMID13, 13, HeatTransferGPC::PYRAMID13, 3>::set();
	GaussPoints<Element::CODE::PRISMA15 , 15, HeatTransferGPC::PRISMA15 , 3>::set();
	GaussPoints<Element::CODE::HEXA20   , 20, HeatTransferGPC::HEXA20   , 3>::set();

	for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		switch (info::mesh->dimension) {
		case 2:
			switch (etype[i]) {
			case HeatTransferElementType::SYMMETRIC_ISOTROPIC : generateBaseFunctions<2, HeatTransferElementType::SYMMETRIC_ISOTROPIC >(i, ops); break;
			case HeatTransferElementType::SYMMETRIC_GENERAL   : generateBaseFunctions<2, HeatTransferElementType::SYMMETRIC_GENERAL   >(i, ops); break;
			case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: generateBaseFunctions<2, HeatTransferElementType::ASYMMETRIC_ISOTROPIC>(i, ops); break;
			case HeatTransferElementType::ASYMMETRIC_GENERAL  : generateBaseFunctions<2, HeatTransferElementType::ASYMMETRIC_GENERAL  >(i, ops); break;
			} break;
		case 3:
			switch (etype[i]) {
			case HeatTransferElementType::SYMMETRIC_ISOTROPIC : generateBaseFunctions<3, HeatTransferElementType::SYMMETRIC_ISOTROPIC >(i, ops); break;
			case HeatTransferElementType::SYMMETRIC_GENERAL   : generateBaseFunctions<3, HeatTransferElementType::SYMMETRIC_GENERAL   >(i, ops); break;
			case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: generateBaseFunctions<3, HeatTransferElementType::ASYMMETRIC_ISOTROPIC>(i, ops); break;
			case HeatTransferElementType::ASYMMETRIC_GENERAL  : generateBaseFunctions<3, HeatTransferElementType::ASYMMETRIC_GENERAL  >(i, ops); break;
			} break;
		}
	}
}

static void generateBaseFunctions(const std::vector<int> &bfilter, std::vector<std::vector<std::vector<ActionOperator*> > > &ops)
{
	for(size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (bfilter[r]) {
			for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
				switch (info::mesh->boundaryRegions[r]->dimension) {
				case 2:
					switch (static_cast<Element::CODE>(info::mesh->boundaryRegions[r]->eintervals[i].code)) {
					case Element::CODE::TRIANGLE3: ops[r][i].push_back(new Basis<Element::CODE::TRIANGLE3,  3, HeatTransferGPC::TRIANGLE3, 2, HeatTransferElementType::FACE, HeatTransferDataDescriptor< 3, HeatTransferGPC::TRIANGLE3, 3, 2, HeatTransferElementType::FACE> >()); break;
					case Element::CODE::TRIANGLE6: ops[r][i].push_back(new Basis<Element::CODE::TRIANGLE6,  6, HeatTransferGPC::TRIANGLE6, 2, HeatTransferElementType::FACE, HeatTransferDataDescriptor< 6, HeatTransferGPC::TRIANGLE6, 3, 2, HeatTransferElementType::FACE> >()); break;
					case Element::CODE::SQUARE4  : ops[r][i].push_back(new Basis<Element::CODE::SQUARE4  ,  4, HeatTransferGPC::SQUARE4  , 2, HeatTransferElementType::FACE, HeatTransferDataDescriptor< 4, HeatTransferGPC::SQUARE4  , 3, 2, HeatTransferElementType::FACE> >()); break;
					case Element::CODE::SQUARE8  : ops[r][i].push_back(new Basis<Element::CODE::SQUARE8  ,  8, HeatTransferGPC::SQUARE8  , 2, HeatTransferElementType::FACE, HeatTransferDataDescriptor< 8, HeatTransferGPC::SQUARE8  , 3, 2, HeatTransferElementType::FACE> >()); break;
					}
				break;
				case 1:
					switch (info::mesh->dimension) {
					case 2:
						switch (static_cast<Element::CODE>(info::mesh->boundaryRegions[r]->eintervals[i].code)) {
						case Element::CODE::LINE2: ops[r][i].push_back(new Basis<Element::CODE::LINE2, 2, HeatTransferGPC::LINE2, 1, HeatTransferElementType::EDGE, HeatTransferDataDescriptor< 2, HeatTransferGPC::LINE2, 2, 1, HeatTransferElementType::EDGE> >()); break;
						case Element::CODE::LINE3: ops[r][i].push_back(new Basis<Element::CODE::LINE3, 3, HeatTransferGPC::LINE3, 1, HeatTransferElementType::EDGE, HeatTransferDataDescriptor< 3, HeatTransferGPC::LINE3, 2, 1, HeatTransferElementType::EDGE> >()); break;
						} break;
					case 3:
						switch (static_cast<Element::CODE>(info::mesh->boundaryRegions[r]->eintervals[i].code)) {
						case Element::CODE::LINE2: ops[r][i].push_back(new Basis<Element::CODE::LINE2, 2, HeatTransferGPC::LINE2, 1, HeatTransferElementType::EDGE, HeatTransferDataDescriptor< 2, HeatTransferGPC::LINE2, 3, 1, HeatTransferElementType::EDGE> >()); break;
						case Element::CODE::LINE3: ops[r][i].push_back(new Basis<Element::CODE::LINE3, 3, HeatTransferGPC::LINE3, 1, HeatTransferElementType::EDGE, HeatTransferDataDescriptor< 3, HeatTransferGPC::LINE3, 3, 1, HeatTransferElementType::EDGE> >()); break;
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
	return new Filler< 1, HeatTransferGPC::POINT1, 2, 0, HeatTransferElementType::NODE, HeatTransferDataDescriptor< 1, HeatTransferGPC::POINT1, 2, 0, HeatTransferElementType::NODE>, Setter>(region, interval, dofs, target, setter);
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Filler, class Target, class Setter>
static ActionOperator* generateNodeSetter3D(size_t region, size_t interval, size_t dofs, Target *target, const Setter &setter)
{
	return new Filler< 1, HeatTransferGPC::POINT1, 3, 0, HeatTransferElementType::NODE, HeatTransferDataDescriptor< 1, HeatTransferGPC::POINT1, 3, 0, HeatTransferElementType::NODE>, Setter>(region, interval, dofs, target, setter);
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
	return new Filler< 1, HeatTransferGPC::POINT1, 2, 0, HeatTransferElementType::NODE, HeatTransferDataDescriptor< 1, HeatTransferGPC::POINT1, 2, 0, HeatTransferElementType::NODE>>(region, interval, dofs, std::forward<Args>(args)...);
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Filler, class ... Args>
static ActionOperator* generateNodeFiller3D(size_t region, size_t interval, size_t dofs, Args&& ... args)
{
	return new Filler< 1, HeatTransferGPC::POINT1, 3, 0, HeatTransferElementType::NODE, HeatTransferDataDescriptor< 1, HeatTransferGPC::POINT1, 3, 0, HeatTransferElementType::NODE>>(region, interval, dofs, std::forward<Args>(args)...);
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
	case static_cast<int>(Element::CODE::LINE2): return new Filler< 2, HeatTransferGPC::LINE2, 2, 1, HeatTransferElementType::EDGE, HeatTransferDataDescriptor< 2, HeatTransferGPC::LINE2, 2, 1, HeatTransferElementType::EDGE> >(region, interval, dofs, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::LINE3): return new Filler< 3, HeatTransferGPC::LINE3, 2, 1, HeatTransferElementType::EDGE, HeatTransferDataDescriptor< 3, HeatTransferGPC::LINE3, 2, 1, HeatTransferElementType::EDGE> >(region, interval, dofs, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Filler, class ... Args>
static ActionOperator* generateEdgeFiller3D(size_t region, size_t interval, size_t dofs, Args&& ... args)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::LINE2): return new Filler< 2, HeatTransferGPC::LINE2, 3, 1, HeatTransferElementType::EDGE, HeatTransferDataDescriptor< 2, HeatTransferGPC::LINE2, 3, 1, HeatTransferElementType::EDGE> >(region, interval, dofs, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::LINE3): return new Filler< 3, HeatTransferGPC::LINE3, 3, 1, HeatTransferElementType::EDGE, HeatTransferDataDescriptor< 3, HeatTransferGPC::LINE3, 3, 1, HeatTransferElementType::EDGE> >(region, interval, dofs, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Filler, class ... Args>
static ActionOperator* generateEdgeFiller(size_t region, size_t interval, size_t dofs, Args&& ... args)
{
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
	case static_cast<int>(Element::CODE::TRIANGLE3): return new Filler< 3, HeatTransferGPC::TRIANGLE3, 3, 2, HeatTransferElementType::FACE, HeatTransferDataDescriptor< 3, HeatTransferGPC::TRIANGLE3, 3, 2, HeatTransferElementType::FACE> >(region, interval, dofs, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::TRIANGLE6): return new Filler< 6, HeatTransferGPC::TRIANGLE6, 3, 2, HeatTransferElementType::FACE, HeatTransferDataDescriptor< 6, HeatTransferGPC::TRIANGLE6, 3, 2, HeatTransferElementType::FACE> >(region, interval, dofs, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::SQUARE4  ): return new Filler< 4, HeatTransferGPC::SQUARE4  , 3, 2, HeatTransferElementType::FACE, HeatTransferDataDescriptor< 4, HeatTransferGPC::SQUARE4  , 3, 2, HeatTransferElementType::FACE> >(region, interval, dofs, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::SQUARE8  ): return new Filler< 8, HeatTransferGPC::SQUARE8  , 3, 2, HeatTransferElementType::FACE, HeatTransferDataDescriptor< 8, HeatTransferGPC::SQUARE8  , 3, 2, HeatTransferElementType::FACE> >(region, interval, dofs, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

static int getThread(size_t interval)
{
	return info::mesh->elements->eintervals[interval].thread;
}

static int getThread(size_t region, size_t interval)
{
	if (info::mesh->boundaryRegions[region]->dimension) {
		return info::mesh->boundaryRegions[region]->eintervals[interval].thread;
	} else {
		return interval;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateTypedExpressionNode2D(size_t region, size_t interval, Evaluator *evaluator, const Setter &setter)
{
	return new Expression< 1, HeatTransferGPC::POINT1, 2, 0, HeatTransferElementType::NODE, HeatTransferDataDescriptor< 1, HeatTransferGPC::POINT1, 2, 0, HeatTransferElementType::NODE>, Setter>(getThread(region, interval), interval, evaluator, setter);
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateTypedExpressionEdge2D(size_t region, size_t interval, Evaluator *evaluator, const Setter &setter)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::LINE2): return new Expression< 2, HeatTransferGPC::LINE2, 2, 1, HeatTransferElementType::EDGE, HeatTransferDataDescriptor< 2, HeatTransferGPC::LINE2, 2, 1, HeatTransferElementType::EDGE>, Setter>(getThread(region, interval), interval, evaluator, setter);
	case static_cast<int>(Element::CODE::LINE3): return new Expression< 3, HeatTransferGPC::LINE3, 2, 1, HeatTransferElementType::EDGE, HeatTransferDataDescriptor< 3, HeatTransferGPC::LINE3, 2, 1, HeatTransferElementType::EDGE>, Setter>(getThread(region, interval), interval, evaluator, setter);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateTypedExpressionNode3D(size_t region, size_t interval, Evaluator *evaluator, const Setter &setter)
{
	return new Expression< 1, HeatTransferGPC::POINT1, 3, 0, HeatTransferElementType::NODE, HeatTransferDataDescriptor< 1, HeatTransferGPC::POINT1, 3, 0, HeatTransferElementType::NODE>, Setter>(getThread(region, interval), interval, evaluator, setter);
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateTypedExpressionEdge3D(size_t region, size_t interval, Evaluator *evaluator, const Setter &setter)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::LINE2): return new Expression< 2, HeatTransferGPC::LINE2, 3, 1, HeatTransferElementType::EDGE, HeatTransferDataDescriptor< 2, HeatTransferGPC::LINE2, 3, 1, HeatTransferElementType::EDGE>, Setter>(getThread(region, interval), interval, evaluator, setter);
	case static_cast<int>(Element::CODE::LINE3): return new Expression< 3, HeatTransferGPC::LINE3, 3, 1, HeatTransferElementType::EDGE, HeatTransferDataDescriptor< 3, HeatTransferGPC::LINE3, 3, 1, HeatTransferElementType::EDGE>, Setter>(getThread(region, interval), interval, evaluator, setter);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateTypedExpressionFace3D(size_t region, size_t interval, Evaluator *evaluator, const Setter &setter)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::TRIANGLE3): return new Expression< 3, HeatTransferGPC::TRIANGLE3, 3, 2, HeatTransferElementType::FACE, HeatTransferDataDescriptor< 3, HeatTransferGPC::TRIANGLE3, 3, 2, HeatTransferElementType::FACE>, Setter>(getThread(region, interval), interval, evaluator, setter);
	case static_cast<int>(Element::CODE::TRIANGLE6): return new Expression< 6, HeatTransferGPC::TRIANGLE6, 3, 2, HeatTransferElementType::FACE, HeatTransferDataDescriptor< 6, HeatTransferGPC::TRIANGLE6, 3, 2, HeatTransferElementType::FACE>, Setter>(getThread(region, interval), interval, evaluator, setter);
	case static_cast<int>(Element::CODE::SQUARE4  ): return new Expression< 4, HeatTransferGPC::SQUARE4  , 3, 2, HeatTransferElementType::FACE, HeatTransferDataDescriptor< 4, HeatTransferGPC::SQUARE4  , 3, 2, HeatTransferElementType::FACE>, Setter>(getThread(region, interval), interval, evaluator, setter);
	case static_cast<int>(Element::CODE::SQUARE8  ): return new Expression< 8, HeatTransferGPC::SQUARE8  , 3, 2, HeatTransferElementType::FACE, HeatTransferDataDescriptor< 8, HeatTransferGPC::SQUARE8  , 3, 2, HeatTransferElementType::FACE>, Setter>(getThread(region, interval), interval, evaluator, setter);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, int etype, class Setter>
static ActionOperator* generateTypedExpression2D(size_t interval, Evaluator *evaluator, const Setter &setter)
{
	switch (info::mesh->elements->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::TRIANGLE3): return new Expression< 3, HeatTransferGPC::TRIANGLE3, 2, 2, etype, HeatTransferDataDescriptor< 3, HeatTransferGPC::TRIANGLE3, 2, 2, etype>, Setter>(getThread(interval), interval, evaluator, setter);
	case static_cast<int>(Element::CODE::TRIANGLE6): return new Expression< 6, HeatTransferGPC::TRIANGLE6, 2, 2, etype, HeatTransferDataDescriptor< 6, HeatTransferGPC::TRIANGLE6, 2, 2, etype>, Setter>(getThread(interval), interval, evaluator, setter);
	case static_cast<int>(Element::CODE::SQUARE4  ): return new Expression< 4, HeatTransferGPC::SQUARE4  , 2, 2, etype, HeatTransferDataDescriptor< 4, HeatTransferGPC::SQUARE4  , 2, 2, etype>, Setter>(getThread(interval), interval, evaluator, setter);
	case static_cast<int>(Element::CODE::SQUARE8  ): return new Expression< 8, HeatTransferGPC::SQUARE8  , 2, 2, etype, HeatTransferDataDescriptor< 8, HeatTransferGPC::SQUARE8  , 2, 2, etype>, Setter>(getThread(interval), interval, evaluator, setter);
	default: return nullptr;
	}
}
template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, int etype, class Setter>
static ActionOperator* generateTypedExpression3D(size_t interval, Evaluator *evaluator, const Setter &setter)
{
	switch (info::mesh->elements->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::TETRA4   ): return new Expression< 4, HeatTransferGPC::TETRA4   , 3, 3, etype, HeatTransferDataDescriptor< 4, HeatTransferGPC::TETRA4   , 3, 3, etype>, Setter>(getThread(interval), interval, evaluator, setter);
	case static_cast<int>(Element::CODE::TETRA10  ): return new Expression<10, HeatTransferGPC::TETRA10  , 3, 3, etype, HeatTransferDataDescriptor<10, HeatTransferGPC::TETRA10  , 3, 3, etype>, Setter>(getThread(interval), interval, evaluator, setter);
	case static_cast<int>(Element::CODE::PYRAMID5 ): return new Expression< 5, HeatTransferGPC::PYRAMID5 , 3, 3, etype, HeatTransferDataDescriptor< 5, HeatTransferGPC::PYRAMID5 , 3, 3, etype>, Setter>(getThread(interval), interval, evaluator, setter);
	case static_cast<int>(Element::CODE::PYRAMID13): return new Expression<13, HeatTransferGPC::PYRAMID13, 3, 3, etype, HeatTransferDataDescriptor<13, HeatTransferGPC::PYRAMID13, 3, 3, etype>, Setter>(getThread(interval), interval, evaluator, setter);
	case static_cast<int>(Element::CODE::PRISMA6  ): return new Expression< 6, HeatTransferGPC::PRISMA6  , 3, 3, etype, HeatTransferDataDescriptor< 6, HeatTransferGPC::PRISMA6  , 3, 3, etype>, Setter>(getThread(interval), interval, evaluator, setter);
	case static_cast<int>(Element::CODE::PRISMA15 ): return new Expression<15, HeatTransferGPC::PRISMA15 , 3, 3, etype, HeatTransferDataDescriptor<15, HeatTransferGPC::PRISMA15 , 3, 3, etype>, Setter>(getThread(interval), interval, evaluator, setter);
	case static_cast<int>(Element::CODE::HEXA8    ): return new Expression< 8, HeatTransferGPC::HEXA8    , 3, 3, etype, HeatTransferDataDescriptor< 8, HeatTransferGPC::HEXA8    , 3, 3, etype>, Setter>(getThread(interval), interval, evaluator, setter);
	case static_cast<int>(Element::CODE::HEXA20   ): return new Expression<20, HeatTransferGPC::HEXA20   , 3, 3, etype, HeatTransferDataDescriptor<20, HeatTransferGPC::HEXA20   , 3, 3, etype>, Setter>(getThread(interval), interval, evaluator, setter);
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
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_ISOTROPIC : return generateTypedExpression<Expression, HeatTransferElementType::SYMMETRIC_ISOTROPIC >(interval, evaluator, setter);
	case HeatTransferElementType::SYMMETRIC_GENERAL   : return generateTypedExpression<Expression, HeatTransferElementType::SYMMETRIC_GENERAL   >(interval, evaluator, setter);
	case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: return generateTypedExpression<Expression, HeatTransferElementType::ASYMMETRIC_ISOTROPIC>(interval, evaluator, setter);
	case HeatTransferElementType::ASYMMETRIC_GENERAL  : return generateTypedExpression<Expression, HeatTransferElementType::ASYMMETRIC_GENERAL  >(interval, evaluator, setter);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateExpression2D(size_t interval, int etype, Evaluator *evaluator, const Setter &setter)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_ISOTROPIC : return generateTypedExpression2D<Expression, HeatTransferElementType::SYMMETRIC_ISOTROPIC >(interval, evaluator, setter);
	case HeatTransferElementType::SYMMETRIC_GENERAL   : return generateTypedExpression2D<Expression, HeatTransferElementType::SYMMETRIC_GENERAL   >(interval, evaluator, setter);
	case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: return generateTypedExpression2D<Expression, HeatTransferElementType::ASYMMETRIC_ISOTROPIC>(interval, evaluator, setter);
	case HeatTransferElementType::ASYMMETRIC_GENERAL  : return generateTypedExpression2D<Expression, HeatTransferElementType::ASYMMETRIC_GENERAL  >(interval, evaluator, setter);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateIsotropicTypeExpression2D(size_t interval, int etype, Evaluator *evaluator, const Setter &setter)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_ISOTROPIC : return generateTypedExpression2D<Expression, HeatTransferElementType::SYMMETRIC_ISOTROPIC >(interval, evaluator, setter);
	case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: return generateTypedExpression2D<Expression, HeatTransferElementType::ASYMMETRIC_ISOTROPIC>(interval, evaluator, setter);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateGeneralTypeExpression2D(size_t interval, int etype, Evaluator *evaluator, const Setter &setter)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_GENERAL   : return generateTypedExpression2D<Expression, HeatTransferElementType::SYMMETRIC_GENERAL   >(interval, evaluator, setter);
	case HeatTransferElementType::ASYMMETRIC_GENERAL  : return generateTypedExpression2D<Expression, HeatTransferElementType::ASYMMETRIC_GENERAL  >(interval, evaluator, setter);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateAsymmetricTypeExpression2D(size_t interval, int etype, Evaluator *evaluator, const Setter &setter)
{
	switch (etype) {
	case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: return generateTypedExpression2D<Expression, HeatTransferElementType::ASYMMETRIC_ISOTROPIC>(interval, evaluator, setter);
	case HeatTransferElementType::ASYMMETRIC_GENERAL  : return generateTypedExpression2D<Expression, HeatTransferElementType::ASYMMETRIC_GENERAL  >(interval, evaluator, setter);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateExpression3D(size_t interval, int etype, Evaluator *evaluator, const Setter &setter)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_ISOTROPIC : return generateTypedExpression3D<Expression, HeatTransferElementType::SYMMETRIC_ISOTROPIC >(interval, evaluator, setter);
	case HeatTransferElementType::SYMMETRIC_GENERAL   : return generateTypedExpression3D<Expression, HeatTransferElementType::SYMMETRIC_GENERAL   >(interval, evaluator, setter);
	case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: return generateTypedExpression3D<Expression, HeatTransferElementType::ASYMMETRIC_ISOTROPIC>(interval, evaluator, setter);
	case HeatTransferElementType::ASYMMETRIC_GENERAL  : return generateTypedExpression3D<Expression, HeatTransferElementType::ASYMMETRIC_GENERAL  >(interval, evaluator, setter);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateIsotropicTypeExpression3D(size_t interval, int etype, Evaluator *evaluator, const Setter &setter)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_ISOTROPIC : return generateTypedExpression3D<Expression, HeatTransferElementType::SYMMETRIC_ISOTROPIC >(interval, evaluator, setter);
	case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: return generateTypedExpression3D<Expression, HeatTransferElementType::ASYMMETRIC_ISOTROPIC>(interval, evaluator, setter);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateGeneralTypeExpression3D(size_t interval, int etype, Evaluator *evaluator, const Setter &setter)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_GENERAL   : return generateTypedExpression3D<Expression, HeatTransferElementType::SYMMETRIC_GENERAL   >(interval, evaluator, setter);
	case HeatTransferElementType::ASYMMETRIC_GENERAL  : return generateTypedExpression3D<Expression, HeatTransferElementType::ASYMMETRIC_GENERAL  >(interval, evaluator, setter);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateAsymmetricTypeExpression3D(size_t interval, int etype, Evaluator *evaluator, const Setter &setter)
{
	switch (etype) {
	case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: return generateTypedExpression3D<Expression, HeatTransferElementType::ASYMMETRIC_ISOTROPIC>(interval, evaluator, setter);
	case HeatTransferElementType::ASYMMETRIC_GENERAL  : return generateTypedExpression3D<Expression, HeatTransferElementType::ASYMMETRIC_GENERAL  >(interval, evaluator, setter);
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
static void generateElementExpression2D(const std::vector<int> &etype, std::vector<std::vector<ActionOperator*> > &ops, std::map<std::string, ECFExpressionVector> &settings, int dimension, const Setter &setter)
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (Assembler::getEvaluator(interval, settings, dimension)) {
			ops[interval].push_back(generateExpression2D<Expression>(interval, etype[interval], Assembler::getEvaluator(interval, settings, dimension), setter));
		}
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static void generateElementExpression3D(const std::vector<int> &etype, std::vector<std::vector<ActionOperator*> > &ops, std::map<std::string, ECFExpressionVector> &settings, int dimension, const Setter &setter)
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (Assembler::getEvaluator(interval, settings, dimension)) {
			ops[interval].push_back(generateExpression3D<Expression>(interval, etype[interval], Assembler::getEvaluator(interval, settings, dimension), setter));
		}
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static void generateElementExpression(const std::vector<int> &etype, std::vector<std::vector<ActionOperator*> > &ops, std::map<std::string, ECFExpressionVector> &settings, int dimension, const Setter &setter)
{
	switch (info::mesh->dimension) {
	case 2: generateElementExpression2D<Expression>(etype, ops, settings, dimension, setter); break;
	case 3: generateElementExpression3D<Expression>(etype, ops, settings, dimension, setter); break;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, int etype, class Setter>
static void generateElementTypedExpression(std::vector<std::vector<ActionOperator*> > &ops, std::map<std::string, ECFExpressionVector> &settings, int dimension, const Setter &setter)
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (Assembler::getEvaluator(interval, settings, dimension)) {
			ops[interval].push_back(generateTypedExpression<Expression, etype>(interval, Assembler::getEvaluator(interval, settings, dimension), setter));
		}
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static void generateElementAsymmetricTypeExpression(const std::vector<int> &etype, std::vector<std::vector<ActionOperator*> > &ops, std::map<std::string, ECFExpressionVector> &settings, int dimension, const Setter &setter)
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (Assembler::getEvaluator(interval, settings, dimension)) {
			switch (info::mesh->dimension) {
			case 2: ops[interval].push_back(generateAsymmetricTypeExpression2D<Expression>(interval, etype[interval], Assembler::getEvaluator(interval, settings, dimension), setter)); break;
			case 3: ops[interval].push_back(generateAsymmetricTypeExpression3D<Expression>(interval, etype[interval], Assembler::getEvaluator(interval, settings, dimension), setter)); break;
			}
		}
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static void generateBoundaryExpression(std::vector<std::vector<std::vector<ActionOperator*> > > &ops, std::map<std::string, ECFExpression> &settings, const Setter &setter)
{
	for(size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		auto it = settings.find(info::mesh->boundaryRegions[r]->name);
		if (it != settings.end()) {
			switch (info::mesh->dimension) {
			case 2:
				switch (info::mesh->boundaryRegions[r]->dimension) {
				case 0:
					for (size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
						ops[r][t].push_back(generateTypedExpressionNode2D<Expression>(r, t, it->second.evaluator, setter));
					}
					break;
				case 1:
					for (size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
						ops[r][interval].push_back(generateTypedExpressionEdge2D<Expression>(r, interval, it->second.evaluator, setter));
					}
					break;
				}
				break;
			case 3:
				switch (info::mesh->boundaryRegions[r]->dimension) {
				case 0:
					for (size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
						ops[r][t].push_back(generateTypedExpressionNode3D<Expression>(r, t, it->second.evaluator, setter));
					}
					break;
				case 1:
					for (size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
						ops[r][interval].push_back(generateTypedExpressionEdge3D<Expression>(r, interval, it->second.evaluator, setter));
					}
					break;
				case 2:
					for (size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
						ops[r][interval].push_back(generateTypedExpressionFace3D<Expression>(r, interval, it->second.evaluator, setter));
					}
					break;
				}
				break;
			}
		}
	}
}


template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, int etype, class ... Args>
static ActionOperator* generateElementTypedOperator2D(size_t interval, Args&& ... args)
{
	switch (info::mesh->elements->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::TRIANGLE3): return new Operator< 3, HeatTransferGPC::TRIANGLE3, 2, 2, etype, HeatTransferDataDescriptor< 3, HeatTransferGPC::TRIANGLE3, 2, 2, etype> >(interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::TRIANGLE6): return new Operator< 6, HeatTransferGPC::TRIANGLE6, 2, 2, etype, HeatTransferDataDescriptor< 6, HeatTransferGPC::TRIANGLE6, 2, 2, etype> >(interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::SQUARE4  ): return new Operator< 4, HeatTransferGPC::SQUARE4  , 2, 2, etype, HeatTransferDataDescriptor< 4, HeatTransferGPC::SQUARE4  , 2, 2, etype> >(interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::SQUARE8  ): return new Operator< 8, HeatTransferGPC::SQUARE8  , 2, 2, etype, HeatTransferDataDescriptor< 8, HeatTransferGPC::SQUARE8  , 2, 2, etype> >(interval, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, int etype, class ... Args>
static ActionOperator* generateElementTypedOperator3D(size_t interval, Args&& ... args)
{
	switch (info::mesh->elements->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::TETRA4   ): return new Operator< 4, HeatTransferGPC::TETRA4   , 3, 3, etype, HeatTransferDataDescriptor< 4, HeatTransferGPC::TETRA4   , 3, 3, etype> >(interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::TETRA10  ): return new Operator<10, HeatTransferGPC::TETRA10  , 3, 3, etype, HeatTransferDataDescriptor<10, HeatTransferGPC::TETRA10  , 3, 3, etype> >(interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::PYRAMID5 ): return new Operator< 5, HeatTransferGPC::PYRAMID5 , 3, 3, etype, HeatTransferDataDescriptor< 5, HeatTransferGPC::PYRAMID5 , 3, 3, etype> >(interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::PYRAMID13): return new Operator<13, HeatTransferGPC::PYRAMID13, 3, 3, etype, HeatTransferDataDescriptor<13, HeatTransferGPC::PYRAMID13, 3, 3, etype> >(interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::PRISMA6  ): return new Operator< 6, HeatTransferGPC::PRISMA6  , 3, 3, etype, HeatTransferDataDescriptor< 6, HeatTransferGPC::PRISMA6  , 3, 3, etype> >(interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::PRISMA15 ): return new Operator<15, HeatTransferGPC::PRISMA15 , 3, 3, etype, HeatTransferDataDescriptor<15, HeatTransferGPC::PRISMA15 , 3, 3, etype> >(interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::HEXA8    ): return new Operator< 8, HeatTransferGPC::HEXA8    , 3, 3, etype, HeatTransferDataDescriptor< 8, HeatTransferGPC::HEXA8    , 3, 3, etype> >(interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::HEXA20   ): return new Operator<20, HeatTransferGPC::HEXA20   , 3, 3, etype, HeatTransferDataDescriptor<20, HeatTransferGPC::HEXA20   , 3, 3, etype> >(interval, std::forward<Args>(args)...);
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
static ActionOperator* generateElementOperator(size_t interval, int etype, Args&& ... args)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_ISOTROPIC : return generateElementTypedOperator<Operator, HeatTransferElementType::SYMMETRIC_ISOTROPIC , Args...>(interval, std::forward<Args>(args)...);
	case HeatTransferElementType::SYMMETRIC_GENERAL   : return generateElementTypedOperator<Operator, HeatTransferElementType::SYMMETRIC_GENERAL   , Args...>(interval, std::forward<Args>(args)...);
	case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: return generateElementTypedOperator<Operator, HeatTransferElementType::ASYMMETRIC_ISOTROPIC, Args...>(interval, std::forward<Args>(args)...);
	case HeatTransferElementType::ASYMMETRIC_GENERAL  : return generateElementTypedOperator<Operator, HeatTransferElementType::ASYMMETRIC_GENERAL  , Args...>(interval, std::forward<Args>(args)...);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static ActionOperator* generateElementOperator2D(size_t interval, int etype, Args&& ... args)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_ISOTROPIC : return generateElementTypedOperator2D<Operator, HeatTransferElementType::SYMMETRIC_ISOTROPIC , Args...>(interval, std::forward<Args>(args)...);
	case HeatTransferElementType::SYMMETRIC_GENERAL   : return generateElementTypedOperator2D<Operator, HeatTransferElementType::SYMMETRIC_GENERAL   , Args...>(interval, std::forward<Args>(args)...);
	case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: return generateElementTypedOperator2D<Operator, HeatTransferElementType::ASYMMETRIC_ISOTROPIC, Args...>(interval, std::forward<Args>(args)...);
	case HeatTransferElementType::ASYMMETRIC_GENERAL  : return generateElementTypedOperator2D<Operator, HeatTransferElementType::ASYMMETRIC_GENERAL  , Args...>(interval, std::forward<Args>(args)...);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static ActionOperator* generateElementAsymmetricOperator(size_t interval, int etype, Args&& ... args)
{
	switch (etype) {
	case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: return generateElementTypedOperator<Operator, HeatTransferElementType::ASYMMETRIC_ISOTROPIC, Args...>(interval, std::forward<Args>(args)...);
	case HeatTransferElementType::ASYMMETRIC_GENERAL  : return generateElementTypedOperator<Operator, HeatTransferElementType::ASYMMETRIC_GENERAL  , Args...>(interval, std::forward<Args>(args)...);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static ActionOperator* generateElementIsotropicTypeOperator(size_t interval, int etype, Args&& ... args)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_ISOTROPIC : return generateElementTypedOperator<Operator, HeatTransferElementType::SYMMETRIC_ISOTROPIC , Args...>(interval, std::forward<Args>(args)...);
	case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: return generateElementTypedOperator<Operator, HeatTransferElementType::ASYMMETRIC_ISOTROPIC, Args...>(interval, std::forward<Args>(args)...);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static ActionOperator* generateElementGeneralTypeOperator(size_t interval, int etype, Args&& ... args)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_GENERAL   : return generateElementTypedOperator<Operator, HeatTransferElementType::SYMMETRIC_GENERAL   , Args...>(interval, std::forward<Args>(args)...);
	case HeatTransferElementType::ASYMMETRIC_GENERAL  : return generateElementTypedOperator<Operator, HeatTransferElementType::ASYMMETRIC_GENERAL  , Args...>(interval, std::forward<Args>(args)...);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static ActionOperator* generateElementGeneralTypeOperator3D(size_t interval, int etype, Args&& ... args)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_GENERAL   : return generateElementTypedOperator3D<Operator, HeatTransferElementType::SYMMETRIC_GENERAL   , Args...>(interval, std::forward<Args>(args)...);
	case HeatTransferElementType::ASYMMETRIC_GENERAL  : return generateElementTypedOperator3D<Operator, HeatTransferElementType::ASYMMETRIC_GENERAL  , Args...>(interval, std::forward<Args>(args)...);
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
static void generateElementAsymmetricOperators(const std::vector<int> &etype, std::vector<std::vector<ActionOperator*> > &ops, Args&& ... args)
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		ops[interval].push_back(generateElementAsymmetricOperator<Operator, Args...>(interval, etype[interval], std::forward<Args>(args)...));
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static ActionOperator* generateBoundaryNode2DOperator(size_t region, size_t interval, Args&& ... args)
{
	return new Operator< 1, HeatTransferGPC::POINT1, 2, 0, HeatTransferElementType::NODE, HeatTransferDataDescriptor< 1, HeatTransferGPC::POINT1, 2, 0, HeatTransferElementType::NODE> >(region, interval, std::forward<Args>(args)...);
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static ActionOperator* generateBoundaryNode3DOperator(size_t region, size_t interval, Args&& ... args)
{
	return new Operator< 1, HeatTransferGPC::POINT1, 3, 0, HeatTransferElementType::NODE, HeatTransferDataDescriptor< 1, HeatTransferGPC::POINT1, 3, 0, HeatTransferElementType::NODE> >(region, interval, std::forward<Args>(args)...);
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static ActionOperator* generateBoundaryEdge2DOperator(size_t region, size_t interval, Args&& ... args)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::LINE2): return new Operator< 2, HeatTransferGPC::LINE2, 2, 1, HeatTransferElementType::EDGE, HeatTransferDataDescriptor< 2, HeatTransferGPC::LINE2, 2, 1, HeatTransferElementType::EDGE> >(region, interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::LINE3): return new Operator< 3, HeatTransferGPC::LINE3, 2, 1, HeatTransferElementType::EDGE, HeatTransferDataDescriptor< 3, HeatTransferGPC::LINE3, 2, 1, HeatTransferElementType::EDGE> >(region, interval, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static ActionOperator* generateBoundaryEdge3DOperator(size_t region, size_t interval, Args&& ... args)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::LINE2): return new Operator< 2, HeatTransferGPC::LINE2, 3, 1, HeatTransferElementType::EDGE, HeatTransferDataDescriptor< 2, HeatTransferGPC::LINE2, 3, 1, HeatTransferElementType::EDGE> >(region, interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::LINE3): return new Operator< 3, HeatTransferGPC::LINE3, 3, 1, HeatTransferElementType::EDGE, HeatTransferDataDescriptor< 3, HeatTransferGPC::LINE3, 3, 1, HeatTransferElementType::EDGE> >(region, interval, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static ActionOperator* generateBoundaryFaceOperator(size_t region, size_t interval, Args&& ... args)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::TRIANGLE3): return new Operator< 3, HeatTransferGPC::TRIANGLE3, 3, 2, HeatTransferElementType::FACE, HeatTransferDataDescriptor< 3, HeatTransferGPC::TRIANGLE3, 3, 2, HeatTransferElementType::FACE> >(region, interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::TRIANGLE6): return new Operator< 6, HeatTransferGPC::TRIANGLE6, 3, 2, HeatTransferElementType::FACE, HeatTransferDataDescriptor< 6, HeatTransferGPC::TRIANGLE6, 3, 2, HeatTransferElementType::FACE> >(region, interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::SQUARE4  ): return new Operator< 4, HeatTransferGPC::SQUARE4  , 3, 2, HeatTransferElementType::FACE, HeatTransferDataDescriptor< 4, HeatTransferGPC::SQUARE4  , 3, 2, HeatTransferElementType::FACE> >(region, interval, std::forward<Args>(args)...);
	case static_cast<int>(Element::CODE::SQUARE8  ): return new Operator< 8, HeatTransferGPC::SQUARE8  , 3, 2, HeatTransferElementType::FACE, HeatTransferDataDescriptor< 8, HeatTransferGPC::SQUARE8  , 3, 2, HeatTransferElementType::FACE> >(region, interval, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static ActionOperator* generateBoundaryOperator(size_t region, size_t interval, Args&& ... args)
{
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
static ActionOperator* generateBoundaryNodeOperator(size_t region, size_t interval, Args&& ... args)
{
	if (info::mesh->boundaryRegions[region]->dimension == 0) {
		switch (info::mesh->dimension) {
		case 2: return generateBoundaryNode2DOperator<Operator>(region, interval, std::forward<Args>(args)...);
		case 3: return generateBoundaryNode3DOperator<Operator>(region, interval, std::forward<Args>(args)...);
		}
	}
	return nullptr;
}


template <template <size_t, size_t, size_t, size_t, size_t, class> class Operator, class ... Args>
static void generateBoundaryOperators(const std::vector<int> &bfilter, std::vector<std::vector<std::vector<ActionOperator*> > > &ops, Args&& ... args)
{
	for(size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (bfilter[r]) {
			for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
				ops[r][i].push_back(generateBoundaryOperator<Operator>(r, i, std::forward<Args>(args)...));
			}
		}
	}
}

template <int etype, class Size, class Parameter>
static void addTypedElementStorage2D(std::vector<ActionOperator*> &ops, size_t interval, const Size &size, const Parameter &parameter)
{
	constexpr size_t ndim = 2, edim = 2;
	size_t elements = info::mesh->elements->eintervals[interval].end - info::mesh->elements->eintervals[interval].begin;
	switch (info::mesh->elements->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::TRIANGLE3): {
		typename HeatTransferDataDescriptor< 3, HeatTransferGPC::TRIANGLE3, ndim, edim, etype>::Element element;
		auto store = new StorageStore< 3, HeatTransferGPC::TRIANGLE3, ndim, edim, etype, HeatTransferDataDescriptor< 3, HeatTransferGPC::TRIANGLE3, ndim, edim, etype>, Parameter>(interval, ops.back(), elements, size(element), parameter);
		auto load  = new StorageLoad < 3, HeatTransferGPC::TRIANGLE3, ndim, edim, etype, HeatTransferDataDescriptor< 3, HeatTransferGPC::TRIANGLE3, ndim, edim, etype>, Parameter>(interval, store, parameter);
		ops.push_back(store);
		ops.push_back(load);
	} break;
	case static_cast<int>(Element::CODE::TRIANGLE6): {
		typename HeatTransferDataDescriptor< 6, HeatTransferGPC::TRIANGLE6, ndim, edim, etype>::Element element;
		auto store = new StorageStore< 6, HeatTransferGPC::TRIANGLE6, ndim, edim, etype, HeatTransferDataDescriptor< 6, HeatTransferGPC::TRIANGLE6, ndim, edim, etype>, Parameter>(interval, ops.back(), elements, size(element), parameter);
		auto load  = new StorageLoad < 6, HeatTransferGPC::TRIANGLE6, ndim, edim, etype, HeatTransferDataDescriptor< 6, HeatTransferGPC::TRIANGLE6, ndim, edim, etype>, Parameter>(interval, store, parameter);
		ops.push_back(store);
		ops.push_back(load);
	} break;
	case static_cast<int>(Element::CODE::SQUARE4): {
		typename HeatTransferDataDescriptor< 4, HeatTransferGPC::SQUARE4,   ndim, edim, etype>::Element element;
		auto store = new StorageStore< 4, HeatTransferGPC::SQUARE4  , ndim, edim, etype, HeatTransferDataDescriptor< 4, HeatTransferGPC::SQUARE4,   ndim, edim, etype>, Parameter>(interval, ops.back(), elements, size(element), parameter);
		auto load  = new StorageLoad < 4, HeatTransferGPC::SQUARE4  , ndim, edim, etype, HeatTransferDataDescriptor< 4, HeatTransferGPC::SQUARE4,   ndim, edim, etype>, Parameter>(interval, store, parameter);
		ops.push_back(store);
		ops.push_back(load);
	} break;
	case static_cast<int>(Element::CODE::SQUARE8): {
		typename HeatTransferDataDescriptor< 8, HeatTransferGPC::SQUARE8,   ndim, edim, etype>::Element element;
		auto store = new StorageStore< 8, HeatTransferGPC::SQUARE8  , ndim, edim, etype, HeatTransferDataDescriptor< 8, HeatTransferGPC::SQUARE8,   ndim, edim, etype>, Parameter>(interval, ops.back(), elements, size(element), parameter);
		auto load  = new StorageLoad < 8, HeatTransferGPC::SQUARE8  , ndim, edim, etype, HeatTransferDataDescriptor< 8, HeatTransferGPC::SQUARE8,   ndim, edim, etype>, Parameter>(interval, store, parameter);
		ops.push_back(store);
		ops.push_back(load);
	} break;
	}
}

template <int etype, class Size, class Parameter>
static void addTypedElementStorage3D(std::vector<ActionOperator*> &ops, size_t interval, const Size &size, const Parameter &parameter)
{
	constexpr size_t ndim = 3, edim = 3;
	size_t elements = info::mesh->elements->eintervals[interval].end - info::mesh->elements->eintervals[interval].begin;
	switch (info::mesh->elements->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::TETRA4): {
		typename HeatTransferDataDescriptor< 4, HeatTransferGPC::TETRA4, ndim, edim, etype>::Element element;
		auto store = new StorageStore< 4, HeatTransferGPC::TETRA4, ndim, edim, etype, HeatTransferDataDescriptor< 4, HeatTransferGPC::TETRA4, ndim, edim, etype>, Parameter>(interval, ops.back(), elements, size(element), parameter);
		auto load  = new StorageLoad < 4, HeatTransferGPC::TETRA4, ndim, edim, etype, HeatTransferDataDescriptor< 4, HeatTransferGPC::TETRA4, ndim, edim, etype>, Parameter>(interval, store, parameter);
		ops.push_back(store);
		ops.push_back(load);
	} break;
	case static_cast<int>(Element::CODE::TETRA10): {
		typename HeatTransferDataDescriptor<10, HeatTransferGPC::TETRA10, ndim, edim, etype>::Element element;
		auto store = new StorageStore<10, HeatTransferGPC::TETRA10, ndim, edim, etype, HeatTransferDataDescriptor<10, HeatTransferGPC::TETRA10, ndim, edim, etype>, Parameter>(interval, ops.back(), elements, size(element), parameter);
		auto load  = new StorageLoad <10, HeatTransferGPC::TETRA10, ndim, edim, etype, HeatTransferDataDescriptor<10, HeatTransferGPC::TETRA10, ndim, edim, etype>, Parameter>(interval, store, parameter);
		ops.push_back(store);
		ops.push_back(load);
	} break;
	case static_cast<int>(Element::CODE::PYRAMID5): {
		typename HeatTransferDataDescriptor< 5, HeatTransferGPC::PYRAMID5, ndim, edim, etype>::Element element;
		auto store = new StorageStore< 5, HeatTransferGPC::PYRAMID5, ndim, edim, etype, HeatTransferDataDescriptor< 5, HeatTransferGPC::PYRAMID5, ndim, edim, etype>, Parameter>(interval, ops.back(), elements, size(element), parameter);
		auto load  = new StorageLoad < 5, HeatTransferGPC::PYRAMID5, ndim, edim, etype, HeatTransferDataDescriptor< 5, HeatTransferGPC::PYRAMID5, ndim, edim, etype>, Parameter>(interval, store, parameter);
		ops.push_back(store);
		ops.push_back(load);
	} break;
	case static_cast<int>(Element::CODE::PYRAMID13): {
		typename HeatTransferDataDescriptor<13, HeatTransferGPC::PYRAMID13, ndim, edim, etype>::Element element;
		auto store = new StorageStore<13, HeatTransferGPC::PYRAMID13, ndim, edim, etype, HeatTransferDataDescriptor<13, HeatTransferGPC::PYRAMID13, ndim, edim, etype>, Parameter>(interval, ops.back(), elements, size(element), parameter);
		auto load  = new StorageLoad <13, HeatTransferGPC::PYRAMID13, ndim, edim, etype, HeatTransferDataDescriptor<13, HeatTransferGPC::PYRAMID13, ndim, edim, etype>, Parameter>(interval, store, parameter);
		ops.push_back(store);
		ops.push_back(load);
	} break;
	case static_cast<int>(Element::CODE::PRISMA6): {
		typename HeatTransferDataDescriptor< 6, HeatTransferGPC::PRISMA6, ndim, edim, etype>::Element element;
		auto store = new StorageStore< 6, HeatTransferGPC::PRISMA6, ndim, edim, etype, HeatTransferDataDescriptor< 6, HeatTransferGPC::PRISMA6, ndim, edim, etype>, Parameter>(interval, ops.back(), elements, size(element), parameter);
		auto load  = new StorageLoad < 6, HeatTransferGPC::PRISMA6, ndim, edim, etype, HeatTransferDataDescriptor< 6, HeatTransferGPC::PRISMA6, ndim, edim, etype>, Parameter>(interval, store, parameter);
		ops.push_back(store);
		ops.push_back(load);
	} break;
	case static_cast<int>(Element::CODE::PRISMA15): {
		typename HeatTransferDataDescriptor<15, HeatTransferGPC::PRISMA15, ndim, edim, etype>::Element element;
		auto store = new StorageStore<15, HeatTransferGPC::PRISMA15, ndim, edim, etype, HeatTransferDataDescriptor<15, HeatTransferGPC::PRISMA15, ndim, edim, etype>, Parameter>(interval, ops.back(), elements, size(element), parameter);
		auto load  = new StorageLoad <15, HeatTransferGPC::PRISMA15, ndim, edim, etype, HeatTransferDataDescriptor<15, HeatTransferGPC::PRISMA15, ndim, edim, etype>, Parameter>(interval, store, parameter);
		ops.push_back(store);
		ops.push_back(load);
	} break;
	case static_cast<int>(Element::CODE::HEXA8): {
		typename HeatTransferDataDescriptor< 8, HeatTransferGPC::HEXA8, ndim, edim, etype>::Element element;
		auto store = new StorageStore< 8, HeatTransferGPC::HEXA8, ndim, edim, etype, HeatTransferDataDescriptor< 8, HeatTransferGPC::HEXA8, ndim, edim, etype>, Parameter>(interval, ops.back(), elements, size(element), parameter);
		auto load  = new StorageLoad < 8, HeatTransferGPC::HEXA8, ndim, edim, etype, HeatTransferDataDescriptor< 8, HeatTransferGPC::HEXA8, ndim, edim, etype>, Parameter>(interval, store, parameter);
		ops.push_back(store);
		ops.push_back(load);
	} break;
	case static_cast<int>(Element::CODE::HEXA20): {
		typename HeatTransferDataDescriptor<20, HeatTransferGPC::HEXA20, ndim, edim, etype>::Element element;
		auto store = new StorageStore<20, HeatTransferGPC::HEXA20, ndim, edim, etype, HeatTransferDataDescriptor<20, HeatTransferGPC::HEXA20, ndim, edim, etype>, Parameter>(interval, ops.back(), elements, size(element), parameter);
		auto load  = new StorageLoad <20, HeatTransferGPC::HEXA20, ndim, edim, etype, HeatTransferDataDescriptor<20, HeatTransferGPC::HEXA20, ndim, edim, etype>, Parameter>(interval, store, parameter);
		ops.push_back(store);
		ops.push_back(load);
	} break;
	}
}

template <class Size, class Parameter>
static void addGeneralTypeElementStorage2D(std::vector<ActionOperator*> &ops, size_t interval, int etype, const Size &size, const Parameter &parameter)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_GENERAL : addTypedElementStorage2D<HeatTransferElementType::SYMMETRIC_GENERAL , Size, Parameter>(ops, interval, size, parameter); break;
	case HeatTransferElementType::ASYMMETRIC_GENERAL: addTypedElementStorage2D<HeatTransferElementType::ASYMMETRIC_GENERAL, Size, Parameter>(ops, interval, size, parameter); break;
	}
}

template <class Size, class Parameter>
static void addGeneralTypeElementStorage3D(std::vector<ActionOperator*> &ops, size_t interval, int etype, const Size &size, const Parameter &parameter)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_GENERAL : addTypedElementStorage3D<HeatTransferElementType::SYMMETRIC_GENERAL , Size, Parameter>(ops, interval, size, parameter); break;
	case HeatTransferElementType::ASYMMETRIC_GENERAL: addTypedElementStorage3D<HeatTransferElementType::ASYMMETRIC_GENERAL, Size, Parameter>(ops, interval, size, parameter); break;
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
			for (size_t i = 0; i < ops[r].size(); ++i) {
				delete ops[r][i].back();
				ops[r][i].pop_back();
			}
		}
	}
}

}


#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_GENERATOR_H_ */
