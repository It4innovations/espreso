
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_GENERATOR_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_GENERATOR_H_

#include "heattransfer.element.h"
#include "analysis/assembler/operator.h"
#include "analysis/assembler/operators/basis.h"

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

static int generateBaseFunctions(const std::vector<int> &etype, std::vector<std::vector<ActionOperator*> > &ops)
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

	return ParameterError::OK;
}

template <template <int, int, int, int, int, class> class Operator, int etype, class ... Args>
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

template <template <int, int, int, int, int, class> class Operator, int etype, class ... Args>
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

template <template <int, int, int, int, int, class> class Operator, class ... Args>
static ActionOperator* generateElementOperator2D(size_t interval, int etype, Args&& ... args)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_ISOTROPIC: return generateElementTypedOperator2D<Operator, HeatTransferElementType::SYMMETRIC_ISOTROPIC, Args...>(interval, std::forward<Args>(args)...);
	case HeatTransferElementType::SYMMETRIC_GENERAL:   return generateElementTypedOperator2D<Operator, HeatTransferElementType::SYMMETRIC_GENERAL  , Args...>(interval, std::forward<Args>(args)...);
	}
	return nullptr;
}

template <template <int, int, int, int, int, class> class Operator, class ... Args>
static ActionOperator* generateElementOperator3D(size_t interval, int etype, Args&& ... args)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_ISOTROPIC: return generateElementTypedOperator3D<Operator, HeatTransferElementType::SYMMETRIC_ISOTROPIC, Args...>(interval, std::forward<Args>(args)...);
	case HeatTransferElementType::SYMMETRIC_GENERAL:   return generateElementTypedOperator3D<Operator, HeatTransferElementType::SYMMETRIC_GENERAL  , Args...>(interval, std::forward<Args>(args)...);
	}
	return nullptr;
}

template <template <int, int, int, int, int, class> class Operator, class ... Args>
static int generateElementOperators2D(const std::vector<int> &etype, std::vector<std::vector<ActionOperator*> > &ops, Args&& ... args)
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		ops[interval].push_back(generateElementOperator2D<Operator, Args...>(interval, etype[interval], std::forward<Args>(args)...));
	}
	return ParameterError::OK;
}

template <template <int, int, int, int, int, class> class Operator, class ... Args>
static int generateElementOperators3D(const std::vector<int> &etype, std::vector<std::vector<ActionOperator*> > &ops, Args&& ... args)
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		ops[interval].push_back(generateElementOperator3D<Operator, Args...>(interval, etype[interval], std::forward<Args>(args)...));
	}
	return ParameterError::OK;
}

template <template <int, int, int, int, int, class> class Operator, class ... Args>
static int generateElementOperators(const std::vector<int> &etype, std::vector<std::vector<ActionOperator*> > &ops, Args&& ... args)
{
	switch (info::mesh->dimension) {
	case 2: generateElementOperators2D<Operator, Args...>(etype, ops, std::forward<Args>(args)...); break;
	case 3: generateElementOperators3D<Operator, Args...>(etype, ops, std::forward<Args>(args)...); break;
	}
	return ParameterError::OK;
}

static int dropLastOperators(std::vector<std::vector<ActionOperator*> > &ops)
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		ops[interval].pop_back();
	}
	return ParameterError::OK;
}

}


#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_GENERATOR_H_ */
