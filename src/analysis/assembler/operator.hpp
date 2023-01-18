
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATOR_HPP_
#define SRC_PHYSICS_ASSEMBLER_OPERATOR_HPP_

#include "controller.h"
#include "operator.h"

#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"

namespace espreso {

template <class NGP, template <size_t, size_t> class Operator, class ... Args>
static inline ActionOperator* _instantiate(int code, size_t interval, Args&& ... args)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::POINT1):    return new Operator< 1, NGP::POINT1   >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::LINE2):     return new Operator< 2, NGP::LINE2    >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::LINE3):     return new Operator< 3, NGP::LINE3    >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return new Operator< 3, NGP::TRIANGLE3>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return new Operator< 6, NGP::TRIANGLE6>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::SQUARE4):   return new Operator< 4, NGP::SQUARE4  >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::SQUARE8):   return new Operator< 8, NGP::SQUARE8  >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TETRA4):    return new Operator< 4, NGP::TETRA4   >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TETRA10):   return new Operator<10, NGP::TETRA10  >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::PYRAMID5):  return new Operator< 5, NGP::PYRAMID5 >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::PYRAMID13): return new Operator<13, NGP::PYRAMID13>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::PRISMA6):   return new Operator< 6, NGP::PRISMA6  >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::PRISMA15):  return new Operator<15, NGP::PRISMA15 >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::HEXA8):     return new Operator< 8, NGP::HEXA8    >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::HEXA20):    return new Operator<20, NGP::HEXA20   >(interval, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <class NGP, template <size_t, size_t> class Operator, class ... Args>
static inline ActionOperator* instantiate(size_t interval, ParameterController &controller, Args&& ... args)
{
	auto op = _instantiate<NGP, Operator, Args...>(info::mesh->elements->eintervals[interval].code, interval, std::forward<Args>(args)...);
	controller.addOperator(op, interval, args...);
	return op;
}

template <class NGP, template <size_t, size_t, size_t, size_t, class> class Operator, template <size_t, size_t, size_t, size_t> class EData, class ... Args>
static inline ActionOperator* _instantiate2D(int code, size_t interval, Args&& ... args)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::POINT1):    return new Operator< 1, NGP::POINT1   , 2, 0, EData< 1, NGP::POINT1   , 2, 0> >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::LINE2):     return new Operator< 2, NGP::LINE2    , 2, 1, EData< 2, NGP::LINE2    , 2, 1> >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::LINE3):     return new Operator< 3, NGP::LINE3    , 2, 1, EData< 3, NGP::LINE3    , 2, 1> >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return new Operator< 3, NGP::TRIANGLE3, 2, 2, EData< 3, NGP::TRIANGLE3, 2, 2> >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return new Operator< 6, NGP::TRIANGLE6, 2, 2, EData< 6, NGP::TRIANGLE6, 2, 2> >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::SQUARE4):   return new Operator< 4, NGP::SQUARE4  , 2, 2, EData< 4, NGP::SQUARE4  , 2, 2> >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::SQUARE8):   return new Operator< 8, NGP::SQUARE8  , 2, 2, EData< 8, NGP::SQUARE8  , 2, 2> >(interval, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <class NGP, template <size_t, size_t, size_t, size_t, class> class Operator, template <size_t, size_t, size_t, size_t> class EData, class ... Args>
static inline ActionOperator* _instantiate3D(int code, size_t interval, Args&& ... args)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::POINT1):    return new Operator< 1, NGP::POINT1   , 3, 0, EData< 1, NGP::POINT1   , 3, 0> >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::LINE2):     return new Operator< 2, NGP::LINE2    , 3, 1, EData< 2, NGP::LINE2    , 3, 1> >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::LINE3):     return new Operator< 3, NGP::LINE3    , 3, 1, EData< 3, NGP::LINE3    , 3, 1> >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return new Operator< 3, NGP::TRIANGLE3, 3, 2, EData< 3, NGP::TRIANGLE3, 3, 2> >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return new Operator< 6, NGP::TRIANGLE6, 3, 2, EData< 6, NGP::TRIANGLE6, 3, 2> >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::SQUARE4):   return new Operator< 4, NGP::SQUARE4  , 3, 2, EData< 4, NGP::SQUARE4  , 3, 2> >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::SQUARE8):   return new Operator< 8, NGP::SQUARE8  , 3, 2, EData< 8, NGP::SQUARE8  , 3, 2> >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TETRA4):    return new Operator< 4, NGP::TETRA4   , 3, 3, EData< 4, NGP::TETRA4   , 3, 3> >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TETRA10):   return new Operator<10, NGP::TETRA10  , 3, 3, EData<10, NGP::TETRA10  , 3, 3> >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::PYRAMID5):  return new Operator< 5, NGP::PYRAMID5 , 3, 3, EData< 5, NGP::PYRAMID5 , 3, 3> >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::PYRAMID13): return new Operator<13, NGP::PYRAMID13, 3, 3, EData<13, NGP::PYRAMID13, 3, 3> >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::PRISMA6):   return new Operator< 6, NGP::PRISMA6  , 3, 3, EData< 6, NGP::PRISMA6  , 3, 3> >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::PRISMA15):  return new Operator<15, NGP::PRISMA15 , 3, 3, EData<15, NGP::PRISMA15 , 3, 3> >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::HEXA8):     return new Operator< 8, NGP::HEXA8    , 3, 3, EData< 8, NGP::HEXA8    , 3, 3> >(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::HEXA20):    return new Operator<20, NGP::HEXA20   , 3, 3, EData<20, NGP::HEXA20   , 3, 3> >(interval, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <class NGP, template <size_t, size_t, size_t, size_t, class> class Operator, template <size_t, size_t, size_t, size_t> class EData, class ... Args>
static inline ActionOperator* instantiate2D(size_t interval, ParameterController &controller, Args&& ... args)
{
	auto op = _instantiate2D<NGP, Operator, EData, Args...>(info::mesh->elements->eintervals[interval].code, interval, std::forward<Args>(args)...);
	controller.addOperator(op, interval, args...);
	return op;
}

template <class NGP, template <size_t, size_t, size_t, size_t, class> class Operator, template <size_t, size_t, size_t, size_t> class EData, class ... Args>
static inline ActionOperator* instantiate3D(size_t interval, ParameterController &controller, Args&& ... args)
{
	auto op = _instantiate3D<NGP, Operator, EData, Args...>(info::mesh->elements->eintervals[interval].code, interval, std::forward<Args>(args)...);
	controller.addOperator(op, interval, args...);
	return op;
}

template <class NGP, template <size_t, size_t, size_t, size_t, class> class Operator, template <size_t, size_t, size_t, size_t> class EData, class ... Args>
static inline ActionOperator* instantiate(size_t interval, ParameterController &controller, Args&& ... args)
{
	switch (info::mesh->dimension) {
	case 2: return instantiate2D<NGP, Operator, EData, Args...>(interval, controller, args...);
	case 3: return instantiate3D<NGP, Operator, EData, Args...>(interval, controller, args...);
	}
	return nullptr;
}

template <class NGP, size_t DIM, template <size_t, size_t, size_t> class Operator, class ... Args>
static inline ActionOperator* _instantiate(int code, size_t interval, Args&& ... args)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::POINT1):    return new Operator< 1, NGP::POINT1   , DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::LINE2):     return new Operator< 2, NGP::LINE2    , DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::LINE3):     return new Operator< 3, NGP::LINE3    , DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return new Operator< 3, NGP::TRIANGLE3, DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return new Operator< 6, NGP::TRIANGLE6, DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::SQUARE4):   return new Operator< 4, NGP::SQUARE4  , DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::SQUARE8):   return new Operator< 8, NGP::SQUARE8  , DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TETRA4):    return new Operator< 4, NGP::TETRA4   , DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TETRA10):   return new Operator<10, NGP::TETRA10  , DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::PYRAMID5):  return new Operator< 5, NGP::PYRAMID5 , DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::PYRAMID13): return new Operator<13, NGP::PYRAMID13, DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::PRISMA6):   return new Operator< 6, NGP::PRISMA6  , DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::PRISMA15):  return new Operator<15, NGP::PRISMA15 , DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::HEXA8):     return new Operator< 8, NGP::HEXA8    , DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::HEXA20):    return new Operator<20, NGP::HEXA20   , DIM>(interval, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <class NGP, size_t DIM, template <size_t, size_t, size_t> class Operator, class ... Args>
static inline ActionOperator* instantiate(size_t interval, ParameterController &controller, Args&& ... args)
{
	auto op = _instantiate<NGP, DIM, Operator, Args...>(info::mesh->elements->eintervals[interval].code, interval, std::forward<Args>(args)...);
	controller.addOperator(op, interval, args...);
	return op;
}

template <class NGP, size_t ndim, template <size_t, size_t, size_t, size_t, class, size_t> class Operator, template <size_t, size_t, size_t, size_t> class EData, size_t DIM, class ... Args>
static inline ActionOperator* _instantiate(int code, size_t interval, Args&& ... args)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::POINT1):    return new Operator< 1, NGP::POINT1   , ndim, 0, EData< 1, NGP::POINT1   , ndim, 0>, DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::LINE2):     return new Operator< 2, NGP::LINE2    , ndim, 1, EData< 2, NGP::LINE2    , ndim, 1>, DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::LINE3):     return new Operator< 3, NGP::LINE3    , ndim, 1, EData< 3, NGP::LINE3    , ndim, 1>, DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return new Operator< 3, NGP::TRIANGLE3, ndim, 2, EData< 3, NGP::TRIANGLE3, ndim, 2>, DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return new Operator< 6, NGP::TRIANGLE6, ndim, 2, EData< 6, NGP::TRIANGLE6, ndim, 2>, DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::SQUARE4):   return new Operator< 4, NGP::SQUARE4  , ndim, 2, EData< 4, NGP::SQUARE4  , ndim, 2>, DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::SQUARE8):   return new Operator< 8, NGP::SQUARE8  , ndim, 2, EData< 8, NGP::SQUARE8  , ndim, 2>, DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TETRA4):    return new Operator< 4, NGP::TETRA4   , ndim, 3, EData< 4, NGP::TETRA4   , ndim, 3>, DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TETRA10):   return new Operator<10, NGP::TETRA10  , ndim, 3, EData<10, NGP::TETRA10  , ndim, 3>, DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::PYRAMID5):  return new Operator< 5, NGP::PYRAMID5 , ndim, 3, EData< 5, NGP::PYRAMID5 , ndim, 3>, DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::PYRAMID13): return new Operator<13, NGP::PYRAMID13, ndim, 3, EData<13, NGP::PYRAMID13, ndim, 3>, DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::PRISMA6):   return new Operator< 6, NGP::PRISMA6  , ndim, 3, EData< 6, NGP::PRISMA6  , ndim, 3>, DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::PRISMA15):  return new Operator<15, NGP::PRISMA15 , ndim, 3, EData<15, NGP::PRISMA15 , ndim, 3>, DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::HEXA8):     return new Operator< 8, NGP::HEXA8    , ndim, 3, EData< 8, NGP::HEXA8    , ndim, 3>, DIM>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::HEXA20):    return new Operator<20, NGP::HEXA20   , ndim, 3, EData<20, NGP::HEXA20   , ndim, 3>, DIM>(interval, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <class NGP, template <size_t, size_t, size_t, size_t, class, size_t> class Operator, template <size_t, size_t, size_t, size_t> class EData, size_t DIM, class ... Args>
static inline ActionOperator* instantiate(size_t interval, ParameterController &controller, Args&& ... args)
{
	ActionOperator* op;
	switch (info::mesh->dimension) {
	case 2:
		op = _instantiate<NGP, 2, Operator, EData, DIM, Args...>(info::mesh->elements->eintervals[interval].code, interval, std::forward<Args>(args)...);
		controller.addOperator(op, interval, args...);
		break;
	case 3:
		op = _instantiate<NGP, 3, Operator, EData, DIM, Args...>(info::mesh->elements->eintervals[interval].code, interval, std::forward<Args>(args)...);
		controller.addOperator(op, interval, args...);
		break;
	}
	return op;
}


template <class NGP, template <size_t N, size_t GP> class Operator, class ... Args>
static inline ActionOperator* instantiate(size_t region, size_t interval, ParameterController &controller, Args&& ... args)
{
	if (info::mesh->boundaryRegions[region]->dimension) {
		auto op = _instantiate<NGP, Operator, Args...>(info::mesh->boundaryRegions[region]->eintervals[interval].code, interval, std::forward<Args>(args)...);
		controller.addOperator(op, interval, args...);
		return op;
	} else {
		auto op = _instantiate<NGP, Operator, Args...>(static_cast<size_t>(Element::CODE::POINT1), interval, std::forward<Args>(args)...);
		controller.addOperator(op, interval, args...);
		return op;
	}
}

template <class NGP, size_t DIM, template <size_t N, size_t GP, size_t dimension> class Operator, class ... Args>
static inline ActionOperator* instantiate(size_t region, size_t interval, ParameterController &controller, Args&& ... args)
{
	if (info::mesh->boundaryRegions[region]->dimension) {
		auto op = _instantiate<NGP, DIM, Operator, Args...>(info::mesh->boundaryRegions[region]->eintervals[interval].code, interval, std::forward<Args>(args)...);
		controller.addOperator(op, interval, args...);
		return op;
	} else {
		auto op = _instantiate<NGP, DIM, Operator, Args...>(static_cast<size_t>(Element::CODE::POINT1), interval, std::forward<Args>(args)...);
		controller.addOperator(op, interval, args...);
		return op;
	}
}

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATOR_HPP_ */
