
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATOR_HPP_
#define SRC_PHYSICS_ASSEMBLER_OPERATOR_HPP_

#include "controller.h"
#include "operator.h"

#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"

namespace espreso {

template <class NGP, template <size_t N, size_t GP> class Operator, class ... Args>
static inline ActionOperator* _instantiate(int code, size_t interval, Args&& ... args)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::LINE2):     return new Operator< 2, NGP::LINE2    >(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::LINE3):     return new Operator< 3, NGP::LINE3    >(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return new Operator< 3, NGP::TRIANGLE3>(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return new Operator< 6, NGP::TRIANGLE6>(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::SQUARE4):   return new Operator< 4, NGP::SQUARE4  >(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::SQUARE8):   return new Operator< 8, NGP::SQUARE8  >(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::TETRA4):    return new Operator< 4, NGP::TETRA4   >(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::TETRA10):   return new Operator<10, NGP::TETRA10  >(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::PYRAMID5):  return new Operator< 5, NGP::PYRAMID5 >(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::PYRAMID13): return new Operator<13, NGP::PYRAMID13>(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::PRISMA6):   return new Operator< 6, NGP::PRISMA6  >(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::PRISMA15):  return new Operator<15, NGP::PRISMA15 >(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::HEXA8):     return new Operator< 8, NGP::HEXA8    >(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::HEXA20):    return new Operator<20, NGP::HEXA20   >(interval, std::forward<Args>(args)...); break;
	default: return nullptr;
	}
}


template <class NGP, size_t DIM, template <size_t N, size_t GP, size_t dimension> class Operator, class ... Args>
static inline ActionOperator* _instantiate(int code, size_t interval, Args&& ... args)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::LINE2):     return new Operator< 2, NGP::LINE2    , DIM>(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::LINE3):     return new Operator< 3, NGP::LINE3    , DIM>(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return new Operator< 3, NGP::TRIANGLE3, DIM>(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return new Operator< 6, NGP::TRIANGLE6, DIM>(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::SQUARE4):   return new Operator< 4, NGP::SQUARE4  , DIM>(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::SQUARE8):   return new Operator< 8, NGP::SQUARE8  , DIM>(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::TETRA4):    return new Operator< 4, NGP::TETRA4   , DIM>(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::TETRA10):   return new Operator<10, NGP::TETRA10  , DIM>(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::PYRAMID5):  return new Operator< 5, NGP::PYRAMID5 , DIM>(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::PYRAMID13): return new Operator<13, NGP::PYRAMID13, DIM>(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::PRISMA6):   return new Operator< 6, NGP::PRISMA6  , DIM>(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::PRISMA15):  return new Operator<15, NGP::PRISMA15 , DIM>(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::HEXA8):     return new Operator< 8, NGP::HEXA8    , DIM>(interval, std::forward<Args>(args)...); break;
	case static_cast<size_t>(Element::CODE::HEXA20):    return new Operator<20, NGP::HEXA20   , DIM>(interval, std::forward<Args>(args)...); break;
	default: return nullptr;
	}
}

template <class Operator, class ... Args>
static inline ActionOperator* instantiate(size_t interval, ParameterController &controller, Args&& ... args)
{
	auto op = new Operator(interval, std::forward<Args>(args)...);
	op->isconst = controller.getConstness(interval, args...);
	return op;
}

template <class NGP, template <size_t N, size_t GP> class Operator, class ... Args>
static inline ActionOperator* instantiate(size_t interval, ParameterController &controller, Args&& ... args)
{
	auto op = _instantiate<NGP, Operator, Args...>(info::mesh->elements->eintervals[interval].code, interval, std::forward<Args>(args)...);
	op->isconst = controller.getConstness(interval, args...);
	return op;
}

template <class NGP, size_t DIM, template <size_t N, size_t GP, size_t dimension> class Operator, class ... Args>
static inline ActionOperator* instantiate(size_t interval, ParameterController &controller, Args&& ... args)
{
	auto op = _instantiate<NGP, DIM, Operator, Args...>(info::mesh->elements->eintervals[interval].code, interval, std::forward<Args>(args)...);
	op->isconst = controller.getConstness(interval, args...);
	return op;
}

template <class NGP, template <size_t N, size_t GP> class Operator, class ... Args>
static inline ActionOperator* instantiate(size_t region, size_t interval, ParameterController &controller, Args&& ... args)
{
	auto op = _instantiate<NGP, Operator, Args...>(info::mesh->boundaryRegions[region]->eintervals[interval].code, interval, std::forward<Args>(args)...);
	op->isconst = controller.getConstness(interval, args...);
	return op;
}

template <class NGP, size_t DIM, template <size_t N, size_t GP, size_t dimension> class Operator, class ... Args>
static inline ActionOperator* instantiate(size_t region, size_t interval, ParameterController &controller, Args&& ... args)
{
	auto op = _instantiate<NGP, DIM, Operator, Args...>(info::mesh->boundaryRegions[region]->eintervals[interval].code, interval, std::forward<Args>(args)...);
	op->isconst = controller.getConstness(interval, args...);
	return op;
}

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATOR_HPP_ */
