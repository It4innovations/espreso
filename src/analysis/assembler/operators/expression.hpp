
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_EXPRESSION_HPP_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_EXPRESSION_HPP_

#include "expression.h"

namespace espreso {

template <class NGP, template <size_t, size_t, size_t, size_t, class, class> class Operator, template <size_t, size_t, size_t, size_t> class EData, class Setter, class ... Args>
static inline ActionOperator* _instantiate3D(int code, size_t interval, Args&& ... args)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::POINT1):    return new Operator< 1, NGP::POINT1   , 3, 0, EData< 1, NGP::POINT1   , 3, 0>, Setter>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::LINE2):     return new Operator< 2, NGP::LINE2    , 3, 1, EData< 2, NGP::LINE2    , 3, 1>, Setter>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::LINE3):     return new Operator< 3, NGP::LINE3    , 3, 1, EData< 3, NGP::LINE3    , 3, 1>, Setter>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return new Operator< 3, NGP::TRIANGLE3, 3, 2, EData< 3, NGP::TRIANGLE3, 3, 2>, Setter>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return new Operator< 6, NGP::TRIANGLE6, 3, 2, EData< 6, NGP::TRIANGLE6, 3, 2>, Setter>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::SQUARE4):   return new Operator< 4, NGP::SQUARE4  , 3, 2, EData< 4, NGP::SQUARE4  , 3, 2>, Setter>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::SQUARE8):   return new Operator< 8, NGP::SQUARE8  , 3, 2, EData< 8, NGP::SQUARE8  , 3, 2>, Setter>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TETRA4):    return new Operator< 4, NGP::TETRA4   , 3, 3, EData< 4, NGP::TETRA4   , 3, 3>, Setter>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TETRA10):   return new Operator<10, NGP::TETRA10  , 3, 3, EData<10, NGP::TETRA10  , 3, 3>, Setter>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::PYRAMID5):  return new Operator< 5, NGP::PYRAMID5 , 3, 3, EData< 5, NGP::PYRAMID5 , 3, 3>, Setter>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::PYRAMID13): return new Operator<13, NGP::PYRAMID13, 3, 3, EData<13, NGP::PYRAMID13, 3, 3>, Setter>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::PRISMA6):   return new Operator< 6, NGP::PRISMA6  , 3, 3, EData< 6, NGP::PRISMA6  , 3, 3>, Setter>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::PRISMA15):  return new Operator<15, NGP::PRISMA15 , 3, 3, EData<15, NGP::PRISMA15 , 3, 3>, Setter>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::HEXA8):     return new Operator< 8, NGP::HEXA8    , 3, 3, EData< 8, NGP::HEXA8    , 3, 3>, Setter>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::HEXA20):    return new Operator<20, NGP::HEXA20   , 3, 3, EData<20, NGP::HEXA20   , 3, 3>, Setter>(interval, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <class NGP, template <size_t, size_t, size_t, size_t, class, class> class Operator, template <size_t, size_t, size_t, size_t> class EData, class Setter, class ... Args>
static inline ActionOperator* _instantiate2D(int code, size_t interval, Args&& ... args)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::POINT1):    return new Operator< 1, NGP::POINT1   , 2, 0, EData< 1, NGP::POINT1   , 2, 0>, Setter>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::LINE2):     return new Operator< 2, NGP::LINE2    , 2, 1, EData< 2, NGP::LINE2    , 2, 1>, Setter>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::LINE3):     return new Operator< 3, NGP::LINE3    , 2, 1, EData< 3, NGP::LINE3    , 2, 1>, Setter>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return new Operator< 3, NGP::TRIANGLE3, 2, 2, EData< 3, NGP::TRIANGLE3, 2, 2>, Setter>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return new Operator< 6, NGP::TRIANGLE6, 2, 2, EData< 6, NGP::TRIANGLE6, 2, 2>, Setter>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::SQUARE4):   return new Operator< 4, NGP::SQUARE4  , 2, 2, EData< 4, NGP::SQUARE4  , 2, 2>, Setter>(interval, std::forward<Args>(args)...);
	case static_cast<size_t>(Element::CODE::SQUARE8):   return new Operator< 8, NGP::SQUARE8  , 2, 2, EData< 8, NGP::SQUARE8  , 2, 2>, Setter>(interval, std::forward<Args>(args)...);
	default: return nullptr;
	}
}

template <class NGP, template <size_t, size_t, size_t, size_t, class, class> class Operator, template <size_t, size_t, size_t, size_t> class EData, class Setter, class ... Args>
static inline ActionOperator* instantiate2D(size_t interval, ParameterController &controller, Args&& ... args)
{
	auto op = _instantiate2D<NGP, Operator, EData, Setter, Args...>(info::mesh->elements->eintervals[interval].code, interval, std::forward<Args>(args)...);
	controller.addOperator(op, interval, args...);
	return op;
}

template <class NGP, template <size_t, size_t, size_t, size_t, class, class> class Operator, template <size_t, size_t, size_t, size_t> class EData, class Setter, class ... Args>
static inline ActionOperator* instantiate3D(size_t interval, ParameterController &controller, Args&& ... args)
{
	auto op = _instantiate3D<NGP, Operator, EData, Setter, Args...>(info::mesh->elements->eintervals[interval].code, interval, std::forward<Args>(args)...);
	controller.addOperator(op, interval, args...);
	return op;
}

template <template <size_t, size_t, size_t, size_t> class Operator, template<size_t, size_t, size_t, size_t, class, class> class Target, class Module, class Setter>
void _fromExpression2D(Module &module, Setter setter, ParameterData &parameter, ExternalElementValue &value)
{
	if (std::all_of(value.evaluator.begin(), value.evaluator.end(), [] (const Evaluator *ev) { return ev == NULL; })) {
		return;
	}

	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		for (int d = 0; d < value.dimension; ++d) {
			if (value.evaluator[i * value.dimension + d]) {
				parameter.update[i] = 1;
				for (size_t p = 0; p < value.evaluator[i * value.dimension + d]->params.general.size(); ++p) {
					module.controller.addInput(i, parameter, value.evaluator[i * value.dimension + d]->params.general[p].variable);
				}
			}
		}
	}

	module.controller.prepare(parameter);

	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		for (int d = 0; d < value.dimension; ++d) {
			if (value.evaluator[i * value.dimension + d]) {
				module.elementOps[i].emplace_back(instantiate2D<typename Module::NGP, Target, Operator, Setter>(i, module.controller, setter, parameter, value.evaluator[i * value.dimension + d], d, value.dimension));
			}
		}
	}
}

template <template <size_t, size_t, size_t, size_t> class Operator, template<size_t, size_t, size_t, size_t, class, class> class Target, class Module, class Setter>
void _fromExpression3D(Module &module, Setter setter, ParameterData &parameter, ExternalElementValue &value)
{
	if (std::all_of(value.evaluator.begin(), value.evaluator.end(), [] (const Evaluator *ev) { return ev == NULL; })) {
		return;
	}

	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		for (int d = 0; d < value.dimension; ++d) {
			if (value.evaluator[i * value.dimension + d]) {
				parameter.update[i] = 1;
				for (size_t p = 0; p < value.evaluator[i * value.dimension + d]->params.general.size(); ++p) {
					module.controller.addInput(i, parameter, value.evaluator[i * value.dimension + d]->params.general[p].variable);
				}
			}
		}
	}

	module.controller.prepare(parameter);

	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		for (int d = 0; d < value.dimension; ++d) {
			if (value.evaluator[i * value.dimension + d]) {
				module.elementOps[i].emplace_back(instantiate3D<typename Module::NGP, Target, Operator, Setter>(i, module.controller, setter, parameter, value.evaluator[i * value.dimension + d], d, value.dimension));
			}
		}
	}
}

template <class Setter>
void fromExpression(HeatTransfer &module, ParameterData &parameter, ExternalElementNodesValue &value, Setter setter)
{
	switch (info::mesh->dimension){
	case 2: _fromExpression2D<HeatTransferOperator, ExpressionsToNodes2>(module, setter, parameter, value); break;
	case 3: _fromExpression3D<HeatTransferOperator, ExpressionsToNodes2>(module, setter, parameter, value); break;
	}
}

template <class Setter>
void fromExpression(HeatTransfer &module, ParameterData &parameter, ExternalElementGPsValue &value, Setter setter)
{
	switch (info::mesh->dimension){
	case 2: _fromExpression2D<HeatTransferOperator, ExpressionsToGPs2>(module, setter, parameter, value); break;
	case 3: _fromExpression3D<HeatTransferOperator, ExpressionsToGPs2>(module, setter, parameter, value); break;
	}
}

template <class Setter>
void fromExpression2D(HeatTransfer &module, ParameterData &parameter, ExternalElementGPsValue &value, Setter setter)
{
	_fromExpression2D<HeatTransferOperator, ExpressionsToGPs2>(module, setter, parameter, value);
}

template <class Setter>
void fromExpression3D(HeatTransfer &module, ParameterData &parameter, ExternalElementGPsValue &value, Setter setter)
{
	_fromExpression3D<HeatTransferOperator, ExpressionsToGPs2>(module, setter, parameter, value);
}

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_EXPRESSION_HPP_ */
