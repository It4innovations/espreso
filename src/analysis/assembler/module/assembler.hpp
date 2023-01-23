
#include "assembler.h"
#include "analysis/assembler/operators/expression.h"

#include "basis/expression/variable.h"
#include "basis/utilities/parser.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/domainstore.h"

using namespace espreso;

template <template <size_t, size_t, size_t, size_t> class Operator, size_t nodes, size_t gps, size_t ndim, size_t edim>
double simdloop(const std::vector<ActionOperator*> &ops, esint elements)
{
	typename Operator<nodes, gps, ndim, edim>::Element data;
	std::vector<Operator<nodes, gps, ndim, edim>*> active; active.reserve(ops.size());

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
//		if ((*op)->update) { // we must always update since clear object on the stact is created
			if ((*op)->isconst) {
				dynamic_cast<Operator<nodes, gps, ndim, edim>*>(*op)->simd(data);
			} else {
				active.push_back(dynamic_cast<Operator<nodes, gps, ndim, edim>*>(*op));
				active.back()->simd(data);
			}
//		}
	}

	double start = eslog::time();
	esint chunks = elements / SIMD::size;
	for (esint c = 1; c < chunks; ++c) {
		for (auto op = active.cbegin(); op != active.cend(); ++op) {
			(*op)->simd(data);
		}
	}
	double end = eslog::time();

	if (elements % SIMD::size) {
		for (auto op = active.cbegin(); op != active.cend(); ++op) {
			(*op)->peel(data, elements % SIMD::size);
		}
	}

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
//		if ((*op)->update) {
			if ((*op)->isconst) {
				(*op)->move(-(int)SIMD::size);
			} else {
				(*op)->move(-elements);
			}
//		}
	}
	return end - start;
}

template <template <size_t, size_t, size_t, size_t> class Operator, size_t nodes, size_t gps, size_t ndim, size_t edim>
double sisdloop(const std::vector<ActionOperator*> &ops, esint elements)
{
	typename Operator<nodes, gps, ndim, edim>::Element data;
	std::vector<Operator<nodes, gps, ndim, edim>*> active; active.reserve(ops.size());

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
//		if ((*op)->update) { // we must always update since clear object on the stact is created
			if ((*op)->isconst) {
				dynamic_cast<Operator<nodes, gps, ndim, edim>*>(*op)->sisd(data);
			} else {
				active.push_back(dynamic_cast<Operator<nodes, gps, ndim, edim>*>(*op));
				active.back()->sisd(data);
			}
//		}
	}

	double start = eslog::time();
	for (esint c = 1; c < elements; ++c) {
		for (auto op = active.cbegin(); op != active.cend(); ++op) {
			(*op)->sisd(data);
		}
	}
	double end = eslog::time();

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
//		if ((*op)->update) {
			if ((*op)->isconst) {
				(*op)->move(-1);
			} else {
				(*op)->move(-elements);
			}
//		}
	}
	return end - start;
}

template <typename Physics, template <size_t, size_t, size_t, size_t> class Operator, size_t ndim>
double simd(int code, const std::vector<ActionOperator*> &ops, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::POINT1):    return simdloop<Operator,  1, Physics::NGP::POINT1    , ndim, 0>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::LINE2):     return simdloop<Operator,  2, Physics::NGP::LINE2     , ndim, 1>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::LINE3):     return simdloop<Operator,  3, Physics::NGP::LINE3     , ndim, 1>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return simdloop<Operator,  3, Physics::NGP::TRIANGLE3 , ndim, 2>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return simdloop<Operator,  6, Physics::NGP::TRIANGLE6 , ndim, 2>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE4):   return simdloop<Operator,  4, Physics::NGP::SQUARE4   , ndim, 2>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE8):   return simdloop<Operator,  8, Physics::NGP::SQUARE8   , ndim, 2>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::TETRA4):    return simdloop<Operator,  4, Physics::NGP::TETRA4    , ndim, 3>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::TETRA10):   return simdloop<Operator, 10, Physics::NGP::TETRA10   , ndim, 3>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID5):  return simdloop<Operator,  5, Physics::NGP::PYRAMID5  , ndim, 3>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID13): return simdloop<Operator, 13, Physics::NGP::PYRAMID13 , ndim, 3>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA6):   return simdloop<Operator,  6, Physics::NGP::PRISMA6   , ndim, 3>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA15):  return simdloop<Operator, 15, Physics::NGP::PRISMA15  , ndim, 3>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA8):     return simdloop<Operator,  8, Physics::NGP::HEXA8     , ndim, 3>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA20):    return simdloop<Operator, 20, Physics::NGP::HEXA20    , ndim, 3>(ops, elements); break;
	}
	return 0;
}

template <typename Physics, template <size_t, size_t, size_t, size_t> class Operator, size_t ndim>
double sisd(int code, const std::vector<ActionOperator*> &ops, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::POINT1):    return sisdloop<Operator,  1, Physics::NGP::POINT1    , ndim, 0>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::LINE2):     return sisdloop<Operator,  2, Physics::NGP::LINE2     , ndim, 1>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::LINE3):     return sisdloop<Operator,  3, Physics::NGP::LINE3     , ndim, 1>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return sisdloop<Operator,  3, Physics::NGP::TRIANGLE3 , ndim, 2>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return sisdloop<Operator,  6, Physics::NGP::TRIANGLE6 , ndim, 2>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE4):   return sisdloop<Operator,  4, Physics::NGP::SQUARE4   , ndim, 2>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE8):   return sisdloop<Operator,  8, Physics::NGP::SQUARE8   , ndim, 2>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::TETRA4):    return sisdloop<Operator,  4, Physics::NGP::TETRA4    , ndim, 3>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::TETRA10):   return sisdloop<Operator, 10, Physics::NGP::TETRA10   , ndim, 3>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID5):  return sisdloop<Operator,  5, Physics::NGP::PYRAMID5  , ndim, 3>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID13): return sisdloop<Operator, 13, Physics::NGP::PYRAMID13 , ndim, 3>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA6):   return sisdloop<Operator,  6, Physics::NGP::PRISMA6   , ndim, 3>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA15):  return sisdloop<Operator, 15, Physics::NGP::PRISMA15  , ndim, 3>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA8):     return sisdloop<Operator,  8, Physics::NGP::HEXA8     , ndim, 3>(ops, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA20):    return sisdloop<Operator, 20, Physics::NGP::HEXA20    , ndim, 3>(ops, elements); break;
	}
	return 0;
}

template <typename Physics, template <size_t, size_t, size_t, size_t> class Operator>
double simd(int code, const std::vector<ActionOperator*> &ops, esint elements)
{
	switch (info::mesh->dimension) {
	case 2: return simd<Physics, Operator, 2>(code, ops, elements); break;
	case 3: return simd<Physics, Operator, 3>(code, ops, elements); break;
	}
	return 0;
}

template <typename Physics, template <size_t, size_t, size_t, size_t> class Operator>
double sisd(int code, const std::vector<ActionOperator*> &ops, esint elements)
{
	switch (info::mesh->dimension) {
	case 2: return sisd<Physics, Operator, 2>(code, ops, elements); break;
	case 3: return sisd<Physics, Operator, 3>(code, ops, elements); break;
	}
	return 0;
}

template <typename Physics, template <size_t, size_t, size_t, size_t> class Operator>
double Assembler::assemble()
{
	double time = 0;
	#pragma omp parallel for reduction(+:time)
	for (int t = 0; t < info::env::threads; ++t) {
		for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
			for (esint i = info::mesh->elements->eintervalsDistribution[d]; i < info::mesh->elements->eintervalsDistribution[d + 1]; ++i) {
				size_t elements = info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin;
				if (settings.simd && elements >= SIMD::size) {
					time += simd<Physics, Operator>(info::mesh->elements->eintervals[i].code, elementOps[i], elements);
				} else {
					time += sisd<Physics, Operator>(info::mesh->elements->eintervals[i].code, elementOps[i], elements);
				}
			}

			for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
				if (info::mesh->boundaryRegions[r]->dimension) {
					for (esint i = info::mesh->boundaryRegions[r]->eintervalsDistribution[d]; i < info::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1]; ++i) {
						size_t elementsInInterval = info::mesh->boundaryRegions[r]->eintervals[i].end - info::mesh->boundaryRegions[r]->eintervals[i].begin;

						for(size_t element = 0; element < elementsInInterval; ++element) {
							for (auto op = boundaryOps[r][i].begin(); op != boundaryOps[r][i].end(); ++op) {
								if((*op)->update) {
									if(element == 0 || !(*op)->isconst) {
										(**op)();
										++(**op);
									}
								}
							}
						}

						for (auto op = boundaryOps[r][i].begin(); op != boundaryOps[r][i].end(); ++op) {
							if((*op)->update) {
								if((*op)->isconst) {
									(*op)->move(-1);
								} else {
									(*op)->move(-elementsInInterval);
								}
							}
						}
					}
				}
			}
		}
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			if (info::mesh->boundaryRegions[r]->dimension == 0) {
				for (auto n = info::mesh->boundaryRegions[r]->nodes->datatarray().begin(t); n != info::mesh->boundaryRegions[r]->nodes->datatarray().end(t); ++n) {
					for (auto op = boundaryOps[r][t].begin(); op != boundaryOps[r][t].end(); ++op) {
						if((*op)->update) {
							if(n == info::mesh->boundaryRegions[r]->nodes->datatarray().begin(t) || !(*op)->isconst) {
								(**op)();
								++(**op);
							}
						}
					}
				}
				for (auto op = boundaryOps[r][t].begin(); op != boundaryOps[r][t].end(); ++op) {
					if((*op)->update) {
						if((*op)->isconst) {
							(*op)->move(-1);
						} else {
							(*op)->move(-info::mesh->boundaryRegions[r]->nodes->datatarray().size(t));
						}
					}
				}
			}
		}
	}

	return time / info::mesh->elements->eintervals.size();
}

template <class Element> struct ElementSize { size_t operator()(const Element& e) { return sizeof(e); } };

template <template <size_t, size_t, size_t, size_t> class Operator, size_t nodes, size_t gps, size_t ndim, size_t edim, template <class> class callback>
size_t simdloop(const std::vector<ActionOperator*> &ops, esint elements)
{
	typename Operator<nodes, gps, ndim, edim>::Element data;
	callback<typename Operator<nodes, gps, ndim, edim>::Element> call;
	return call(data);
}

template <typename Physics, template <size_t, size_t, size_t, size_t> class Operator, size_t ndim, template <class> class callback>
size_t simd(int code, const std::vector<ActionOperator*> &ops, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::POINT1):    return simdloop<Operator,  1, Physics::NGP::POINT1    , ndim, 0, callback>(ops, elements);
	case static_cast<size_t>(Element::CODE::LINE2):     return simdloop<Operator,  2, Physics::NGP::LINE2     , ndim, 1, callback>(ops, elements);
	case static_cast<size_t>(Element::CODE::LINE3):     return simdloop<Operator,  3, Physics::NGP::LINE3     , ndim, 1, callback>(ops, elements);
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return simdloop<Operator,  3, Physics::NGP::TRIANGLE3 , ndim, 2, callback>(ops, elements);
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return simdloop<Operator,  6, Physics::NGP::TRIANGLE6 , ndim, 2, callback>(ops, elements);
	case static_cast<size_t>(Element::CODE::SQUARE4):   return simdloop<Operator,  4, Physics::NGP::SQUARE4   , ndim, 2, callback>(ops, elements);
	case static_cast<size_t>(Element::CODE::SQUARE8):   return simdloop<Operator,  8, Physics::NGP::SQUARE8   , ndim, 2, callback>(ops, elements);
	case static_cast<size_t>(Element::CODE::TETRA4):    return simdloop<Operator,  4, Physics::NGP::TETRA4    , ndim, 3, callback>(ops, elements);
	case static_cast<size_t>(Element::CODE::TETRA10):   return simdloop<Operator, 10, Physics::NGP::TETRA10   , ndim, 3, callback>(ops, elements);
	case static_cast<size_t>(Element::CODE::PYRAMID5):  return simdloop<Operator,  5, Physics::NGP::PYRAMID5  , ndim, 3, callback>(ops, elements);
	case static_cast<size_t>(Element::CODE::PYRAMID13): return simdloop<Operator, 13, Physics::NGP::PYRAMID13 , ndim, 3, callback>(ops, elements);
	case static_cast<size_t>(Element::CODE::PRISMA6):   return simdloop<Operator,  6, Physics::NGP::PRISMA6   , ndim, 3, callback>(ops, elements);
	case static_cast<size_t>(Element::CODE::PRISMA15):  return simdloop<Operator, 15, Physics::NGP::PRISMA15  , ndim, 3, callback>(ops, elements);
	case static_cast<size_t>(Element::CODE::HEXA8):     return simdloop<Operator,  8, Physics::NGP::HEXA8     , ndim, 3, callback>(ops, elements);
	case static_cast<size_t>(Element::CODE::HEXA20):    return simdloop<Operator, 20, Physics::NGP::HEXA20    , ndim, 3, callback>(ops, elements);
	}
	return 0;
}

template <typename Physics, template <size_t, size_t, size_t, size_t> class Operator, template <class> class callback>
size_t simd(int code, const std::vector<ActionOperator*> &ops, esint elements)
{
	switch (info::mesh->dimension) {
	case 2: return simd<Physics, Operator, 2, callback>(code, ops, elements);
	case 3: return simd<Physics, Operator, 3, callback>(code, ops, elements);
	}
	return 0;
}

template <typename Physics, template <size_t nodes, size_t gps, size_t ndim, size_t edim> class Operator>
size_t Assembler::esize()
{
	size_t esize = 0;
	#pragma omp parallel for reduction(max:esize)
	for (int t = 0; t < info::env::threads; ++t) {
		for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
			for (esint i = info::mesh->elements->eintervalsDistribution[d]; i < info::mesh->elements->eintervalsDistribution[d + 1]; ++i) {
				size_t elements = info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin;
				esize = std::max(esize, simd<Physics, Operator, ElementSize>(info::mesh->elements->eintervals[i].code, elementOps[i], elements));
			}
		}
	}
	return esize;
}

template<typename Ttype>
void Assembler::validateRegionSettings(const std::string &name, const std::map<std::string, Ttype> &settings)
{
	for (auto ei = info::mesh->elements->eintervals.begin(); ei != info::mesh->elements->eintervals.end(); ++ei) {
		if (ei->region == -1) { // intersected regions
			std::vector<int> set;
			for (auto rindex = ei->regions.begin(); rindex != ei->regions.end(); ++rindex) {
				if (settings.find(info::mesh->elementsRegions[*rindex]->name) != settings.end()) {
					set.push_back(*rindex);
				}
			}
			if (set.size() > 1) {
				std::string regions = info::mesh->elementsRegions[set[0]]->name + ", " + info::mesh->elementsRegions[set[1]]->name;
				eslog::info("  INVALID %s SETTINGS (INTERSECTED REGIONS): %*s\n", name.c_str(), 50 - name.size(), regions.c_str());
				eslog::info(" ============================================================================================= \n");
				eslog::error("ESPRESO error: invalid configuration.\n");
			}
		}
	}
}

template<class TSecond>
bool Assembler::examineElementParameter(const std::string &name, std::map<std::string, TSecond> &settings, ExternalElementValue &value, int dimension, std::function<ECFExpression*(TSecond &expr)> getExpr)
{
	if (settings.size() == 1 && StringCompare::caseInsensitiveEq(settings.begin()->first, "ALL_ELEMENTS")) {
		ECFExpression* expr = getExpr(settings.begin()->second);
		if (expr->evaluator == nullptr) {
			if (!Variable::create(*expr)) {
				eslog::warning("  %s:  %*s \n", name.c_str(), 88 - name.size(), "INVALID EXPRESSION");
				return false;
			}
		}
		for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			value.evaluator[value.dimension * i + dimension] = expr->evaluator;
		}
		if (expr->evaluator->variables.size()) {
			std::string params = Parser::join(", ", expr->evaluator->variables);
			eslog::info("  %s:  %*s       FNC( %s )\n", name.c_str(), 55 - params.size(), "", params.c_str());
		} else {
			eslog::info("  %s:  %*g \n", name.c_str(), 88 - name.size(), expr->evaluator->eval(Evaluator::Params()));
		}
		return true;
	} else {
		if (settings.size() == 0) {
			eslog::info("  %s:  %*s \n", name.c_str(), 88 - name.size(), "UNKNOWN");
			return false;
		} else {
			eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");
			int rindex = 1, rlast = info::mesh->elementsRegions.size() - 1;
			for (auto reg = info::mesh->elementsRegions.begin() + 1; reg != info::mesh->elementsRegions.end(); ++reg, ++rindex) {
				auto ms = settings.find((*reg)->name);
				if (ms == settings.end()) {
					ms = settings.find("ALL_ELEMENTS");
				}
				if (ms != settings.end()) {
					ECFExpression* expr = getExpr(ms->second);
					if (expr->evaluator == nullptr) {
						if (!Variable::create(*expr)) {
							eslog::warning("   %30s:  %*s \n", (*reg)->name.c_str(), 60 - (*reg)->name.size(), "INVALID EXPRESSION");
							return false;
						}
					}
					if (expr->evaluator->variables.size()) {
						std::string params = Parser::join(", ", expr->evaluator->variables);
						eslog::info("   %30s:  %*s       FNC( %s )\n", (*reg)->name.c_str(), 43 - params.size(), "", params.c_str());
					} else {
						eslog::info("   %30s:  %57g \n", (*reg)->name.c_str(), expr->evaluator->eval(Evaluator::Params()));
					}
					for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
						if (info::mesh->elements->eintervals[i].region == rindex || (info::mesh->elements->eintervals[i].region == 0 && rindex == rlast)) {
							value.evaluator[value.dimension * i + dimension] = expr->evaluator;
						}
						if (info::mesh->elements->eintervals[i].region == -1) {
							const std::vector<int> &regions = info::mesh->elements->eintervals[i].regions;
							bool other = false;
							for (auto it = regions.begin(); it != regions.end(); ++it) {
								if (*it != rindex && settings.find((info::mesh->elementsRegions[*it])->name) != settings.end()) {
									other = true;
								}
							}
							if (!other) {
								value.evaluator[value.dimension * i + dimension] = expr->evaluator;
							}
						}
					}
					return true;
				} else {
					eslog::info("  %s:  %*s \n", name.c_str(), 88 - name.size(), "UNKNOWN");
					return false;
				}
			}
			eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		}
	}
	return false;
}

template<class TSecond>
bool Assembler::examineBoundaryParameter(const std::string &name, std::map<std::string, TSecond> &settings, ExternalBoundaryValue &value, int dimension, std::function<ECFExpression*(TSecond &expr)> getExpr)
{
	if (settings.size()) {
		eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");

		int rindex = 0;
		for (auto reg = info::mesh->boundaryRegions.begin(); reg != info::mesh->boundaryRegions.end(); ++reg, ++rindex) {
			auto ms = settings.find((*reg)->name);
			if (ms != settings.end()) {
				ECFExpression* expr = getExpr(ms->second);
				if (expr->isset) {
					if (expr->evaluator == nullptr) {
						if (!Variable::create(*expr, rindex)) {
							eslog::warning("   %30s:  %*s \n", (*reg)->name.c_str(), 60 - (*reg)->name.size(), "INVALID EXPRESSION");
							return false;
						}
					}

					if (expr->evaluator->variables.size()) {
						std::string params = Parser::join(", ", expr->evaluator->variables);
						eslog::info("   %30s:  %*s       FNC( %s )\n", (*reg)->name.c_str(), 43 - params.size(), "", params.c_str());
					} else {
						eslog::info("   %30s:  %57g \n", (*reg)->name.c_str(), expr->evaluator->eval(Evaluator::Params()));
					}
					value.evaluator[value.dimension * rindex + dimension] = expr->evaluator;
				}
			}
		}
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	}
	return true;
}
