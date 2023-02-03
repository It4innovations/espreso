
#include "assembler.h"

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

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, int nodes, int gps, int ndim, int edim, int etype>
double simdloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, esint elements)
{
	if (elements == 0) return 0;
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
	return end - start;
}

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, int nodes, int gps, int ndim, int edim, int etype>
double sisdloop(const std::vector<ActionOperator*> &ops, esint elements)
{
	typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element element;
	std::vector<DataDescriptor<nodes, gps, ndim, edim, etype>*> active; active.reserve(ops.size());

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
//		if ((*op)->update) { // we must always update since clear object on the stact is created
			if ((*op)->isconst) {
				dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, etype>*>(*op)->sisd(element);
			} else {
				active.push_back(dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, etype>*>(*op));
				active.back()->sisd(element);
			}
//		}
	}

	double start = eslog::time();
	for (esint c = 1; c < elements; ++c) {
		for (auto op = active.cbegin(); op != active.cend(); ++op) {
			(*op)->sisd(element);
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

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
double Assembler::loop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, esint elements)
{
	if (settings.simd && elements >= SIMD::size) {
		return simdloop<DataDescriptor, nodes, gps, ndim, edim, etype>(action, ops, elements);
	} else {
		return simdloop<DataDescriptor, nodes, gps, ndim, edim, etype>(action, ops, elements);
	}
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
