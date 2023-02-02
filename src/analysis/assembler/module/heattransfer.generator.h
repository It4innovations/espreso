
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_GENERATOR_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_GENERATOR_H_

#include "heattransfer.element.h"
#include "analysis/assembler/operator.h"
#include "analysis/assembler/operators/basis.h"
#include "analysis/assembler/operators/expression.h"

#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"

namespace espreso {

template <int ndim, int etype>
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
	GaussPoints<Element::CODE::LINE2, 2, 2, 1>::set();
	GaussPoints<Element::CODE::TRIANGLE3, 3, 6, 2>::set();
	GaussPoints<Element::CODE::SQUARE4, 4, 4, 2>::set();
	GaussPoints<Element::CODE::TETRA4, 4, 4, 3>::set();
	GaussPoints<Element::CODE::PYRAMID5, 5, 8, 3>::set();
	GaussPoints<Element::CODE::PRISMA6, 6, 9, 3>::set();
	GaussPoints<Element::CODE::HEXA8, 8, 8, 3>::set();
	GaussPoints<Element::CODE::LINE3, 3, 3, 1>::set();
	GaussPoints<Element::CODE::TRIANGLE6, 6, 6, 2>::set();
	GaussPoints<Element::CODE::SQUARE8, 8, 9, 2>::set();
	GaussPoints<Element::CODE::TETRA10, 10, 15, 3>::set();
	GaussPoints<Element::CODE::PYRAMID13, 13, 14, 3>::set();
	GaussPoints<Element::CODE::PRISMA15, 15, 9, 3>::set();
	GaussPoints<Element::CODE::HEXA20, 20, 8, 3>::set();

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		switch (info::mesh->dimension) {
		case 2:
			switch (etype[interval]) {
			case HeatTransferElementType::SYMMETRIC_ISOTROPIC: generateBaseFunctions<2, HeatTransferElementType::SYMMETRIC_ISOTROPIC>(interval, ops); break;
			case HeatTransferElementType::SYMMETRIC_GENERAL:   generateBaseFunctions<2, HeatTransferElementType::SYMMETRIC_GENERAL  >(interval, ops); break;
			} break;
		case 3:
			switch (etype[interval]) {
			case HeatTransferElementType::SYMMETRIC_ISOTROPIC: generateBaseFunctions<3, HeatTransferElementType::SYMMETRIC_ISOTROPIC>(interval, ops); break;
			case HeatTransferElementType::SYMMETRIC_GENERAL:   generateBaseFunctions<3, HeatTransferElementType::SYMMETRIC_GENERAL  >(interval, ops); break;
			} break;
		}
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Target, class Setter>
static ActionOperator* generateNodeSetter2D(size_t region, size_t interval, size_t dofs, Target *target, const Setter &setter)
{
	return new Expression< 1, HeatTransferGPC::POINT1, 2, 0, HeatTransferElementType::NODE, HeatTransferDataDescriptor< 1, HeatTransferGPC::POINT1, 2, 0, HeatTransferElementType::NODE>, Setter>(region, interval, dofs, target, setter);
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Target, class Setter>
static ActionOperator* generateNodeSetter3D(size_t region, size_t interval, size_t dofs, Target *target, const Setter &setter)
{
	return new Expression< 1, HeatTransferGPC::POINT1, 3, 0, HeatTransferElementType::NODE, HeatTransferDataDescriptor< 1, HeatTransferGPC::POINT1, 3, 0, HeatTransferElementType::NODE>, Setter>(region, interval, dofs, target, setter);
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateTypedExpressionNode2D(size_t region, size_t interval, Evaluator *evaluator, const Setter &setter)
{
	return new Expression< 1, HeatTransferGPC::POINT1, 2, 0, HeatTransferElementType::NODE, HeatTransferDataDescriptor< 1, HeatTransferGPC::POINT1, 2, 0, HeatTransferElementType::NODE>, Setter>(interval, evaluator, setter);
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateTypedExpressionEdge2D(size_t region, size_t interval, Evaluator *evaluator, const Setter &setter)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::LINE2): return new Expression< 2, HeatTransferGPC::LINE2, 2, 1, HeatTransferElementType::EDGE, HeatTransferDataDescriptor< 2, HeatTransferGPC::LINE2, 2, 1, HeatTransferElementType::EDGE>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::LINE3): return new Expression< 3, HeatTransferGPC::LINE3, 2, 1, HeatTransferElementType::EDGE, HeatTransferDataDescriptor< 3, HeatTransferGPC::LINE3, 2, 1, HeatTransferElementType::EDGE>, Setter>(interval, evaluator, setter);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateTypedExpressionNode3D(size_t region, size_t interval, Evaluator *evaluator, const Setter &setter)
{
	return new Expression< 1, HeatTransferGPC::POINT1, 3, 0, HeatTransferElementType::NODE, HeatTransferDataDescriptor< 1, HeatTransferGPC::POINT1, 3, 0, HeatTransferElementType::NODE>, Setter>(interval, evaluator, setter);
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateTypedExpressionEdge3D(size_t region, size_t interval, Evaluator *evaluator, const Setter &setter)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::LINE2): return new Expression< 2, HeatTransferGPC::LINE2, 3, 1, HeatTransferElementType::EDGE, HeatTransferDataDescriptor< 2, HeatTransferGPC::LINE2, 3, 1, HeatTransferElementType::EDGE>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::LINE3): return new Expression< 3, HeatTransferGPC::LINE3, 3, 1, HeatTransferElementType::EDGE, HeatTransferDataDescriptor< 3, HeatTransferGPC::LINE3, 3, 1, HeatTransferElementType::EDGE>, Setter>(interval, evaluator, setter);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateTypedExpressionFace3D(size_t region, size_t interval, Evaluator *evaluator, const Setter &setter)
{
	switch (info::mesh->boundaryRegions[region]->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::TRIANGLE3): return new Expression< 3, HeatTransferGPC::TRIANGLE3, 3, 2, HeatTransferElementType::FACE, HeatTransferDataDescriptor< 3, HeatTransferGPC::TRIANGLE3, 3, 2, HeatTransferElementType::FACE>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::TRIANGLE6): return new Expression< 6, HeatTransferGPC::TRIANGLE6, 3, 2, HeatTransferElementType::FACE, HeatTransferDataDescriptor< 6, HeatTransferGPC::TRIANGLE6, 3, 2, HeatTransferElementType::FACE>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::SQUARE4  ): return new Expression< 4, HeatTransferGPC::SQUARE4  , 3, 2, HeatTransferElementType::FACE, HeatTransferDataDescriptor< 4, HeatTransferGPC::SQUARE4  , 3, 2, HeatTransferElementType::FACE>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::SQUARE8  ): return new Expression< 8, HeatTransferGPC::SQUARE8  , 3, 2, HeatTransferElementType::FACE, HeatTransferDataDescriptor< 8, HeatTransferGPC::SQUARE8  , 3, 2, HeatTransferElementType::FACE>, Setter>(interval, evaluator, setter);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, int etype, class Setter>
static ActionOperator* generateTypedExpression2D(size_t interval, Evaluator *evaluator, const Setter &setter)
{
	switch (info::mesh->elements->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::TRIANGLE3): return new Expression< 3, HeatTransferGPC::TRIANGLE3, 2, 2, etype, HeatTransferDataDescriptor< 3, HeatTransferGPC::TRIANGLE3, 2, 2, etype>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::TRIANGLE6): return new Expression< 6, HeatTransferGPC::TRIANGLE6, 2, 2, etype, HeatTransferDataDescriptor< 6, HeatTransferGPC::TRIANGLE6, 2, 2, etype>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::SQUARE4  ): return new Expression< 4, HeatTransferGPC::SQUARE4  , 2, 2, etype, HeatTransferDataDescriptor< 4, HeatTransferGPC::SQUARE4  , 2, 2, etype>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::SQUARE8  ): return new Expression< 8, HeatTransferGPC::SQUARE8  , 2, 2, etype, HeatTransferDataDescriptor< 8, HeatTransferGPC::SQUARE8  , 2, 2, etype>, Setter>(interval, evaluator, setter);
	default: return nullptr;
	}
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, int etype, class Setter>
static ActionOperator* generateTypedExpression3D(size_t interval, Evaluator *evaluator, const Setter &setter)
{
	switch (info::mesh->elements->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::TETRA4   ): return new Expression< 4, HeatTransferGPC::TETRA4   , 3, 3, etype, HeatTransferDataDescriptor< 4, HeatTransferGPC::TETRA4   , 3, 3, etype>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::TETRA10  ): return new Expression<10, HeatTransferGPC::TETRA10  , 3, 3, etype, HeatTransferDataDescriptor<10, HeatTransferGPC::TETRA10  , 3, 3, etype>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::PYRAMID5 ): return new Expression< 5, HeatTransferGPC::PYRAMID5 , 3, 3, etype, HeatTransferDataDescriptor< 5, HeatTransferGPC::PYRAMID5 , 3, 3, etype>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::PYRAMID13): return new Expression<13, HeatTransferGPC::PYRAMID13, 3, 3, etype, HeatTransferDataDescriptor<13, HeatTransferGPC::PYRAMID13, 3, 3, etype>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::PRISMA6  ): return new Expression< 6, HeatTransferGPC::PRISMA6  , 3, 3, etype, HeatTransferDataDescriptor< 6, HeatTransferGPC::PRISMA6  , 3, 3, etype>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::PRISMA15 ): return new Expression<15, HeatTransferGPC::PRISMA15 , 3, 3, etype, HeatTransferDataDescriptor<15, HeatTransferGPC::PRISMA15 , 3, 3, etype>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::HEXA8    ): return new Expression< 8, HeatTransferGPC::HEXA8    , 3, 3, etype, HeatTransferDataDescriptor< 8, HeatTransferGPC::HEXA8    , 3, 3, etype>, Setter>(interval, evaluator, setter);
	case static_cast<int>(Element::CODE::HEXA20   ): return new Expression<20, HeatTransferGPC::HEXA20   , 3, 3, etype, HeatTransferDataDescriptor<20, HeatTransferGPC::HEXA20   , 3, 3, etype>, Setter>(interval, evaluator, setter);
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
	case HeatTransferElementType::SYMMETRIC_ISOTROPIC: return generateTypedExpression<Expression, HeatTransferElementType::SYMMETRIC_ISOTROPIC>(interval, evaluator, setter);
	case HeatTransferElementType::SYMMETRIC_GENERAL:   return generateTypedExpression<Expression, HeatTransferElementType::SYMMETRIC_GENERAL  >(interval, evaluator, setter);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateExpression2D(size_t interval, int etype, Evaluator *evaluator, const Setter &setter)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_ISOTROPIC: return generateTypedExpression2D<Expression, HeatTransferElementType::SYMMETRIC_ISOTROPIC>(interval, evaluator, setter);
	case HeatTransferElementType::SYMMETRIC_GENERAL:   return generateTypedExpression2D<Expression, HeatTransferElementType::SYMMETRIC_GENERAL  >(interval, evaluator, setter);
	}
	return nullptr;
}

template <template <size_t, size_t, size_t, size_t, size_t, class, class> class Expression, class Setter>
static ActionOperator* generateExpression3D(size_t interval, int etype, Evaluator *evaluator, const Setter &setter)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_ISOTROPIC: return generateTypedExpression3D<Expression, HeatTransferElementType::SYMMETRIC_ISOTROPIC>(interval, evaluator, setter);
	case HeatTransferElementType::SYMMETRIC_GENERAL:   return generateTypedExpression3D<Expression, HeatTransferElementType::SYMMETRIC_GENERAL  >(interval, evaluator, setter);
	}
	return nullptr;
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
	case HeatTransferElementType::SYMMETRIC_ISOTROPIC: return generateElementTypedOperator<Operator, HeatTransferElementType::SYMMETRIC_ISOTROPIC, Args...>(interval, std::forward<Args>(args)...);
	case HeatTransferElementType::SYMMETRIC_GENERAL:   return generateElementTypedOperator<Operator, HeatTransferElementType::SYMMETRIC_GENERAL  , Args...>(interval, std::forward<Args>(args)...);
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

static void dropLastOperators(std::vector<std::vector<ActionOperator*> > &ops)
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		ops[interval].pop_back();
	}
}

}


#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_GENERATOR_H_ */
