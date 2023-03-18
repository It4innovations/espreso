
#include "assembler.hpp"

#include "analysis/assembler/operators/expression.h"
#include "basis/evaluator/evaluator.h"
#include "basis/utilities/parser.h"
#include "basis/utilities/utils.h"
#include "config/holders/expression.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"

#include "wrappers/mpi/communication.h"
#include "wrappers/exprtk/exprtk.h"

#include <algorithm>
#include <numeric>

#include <iostream>

using namespace espreso;

Assembler::Assembler(PhysicsConfiguration &settings)
: settings(settings),
  etype(info::mesh->elements->eintervals.size()), bfilter(info::mesh->boundaryRegions.size()), btype(info::mesh->boundaryRegions.size()),
  elementOps(info::mesh->elements->eintervals.size()),
  boundaryOps(info::mesh->boundaryRegions.size()),
  K(nullptr), f(nullptr)
{
	for (size_t i = 0; i < info::mesh->boundaryRegions.size(); ++i) {
		if (info::mesh->boundaryRegions[i]->dimension) {
			btype[i].resize(info::mesh->boundaryRegions[i]->eintervals.size());
			boundaryOps[i].resize(info::mesh->boundaryRegions[i]->eintervals.size());
		} else {
			btype[i].resize(info::mesh->boundaryRegions[i]->nodes->threads());
			boundaryOps[i].resize(info::mesh->boundaryRegions[i]->nodes->threads());
		}
	}
}

Assembler::~Assembler()
{
	for (size_t i = 0; i < elementOps.size(); ++i) {
		for (size_t j = 0; j < elementOps[i].size(); ++j) {
			delete elementOps[i][j];
		}
	}
	for (size_t r = 0; r < boundaryOps.size(); ++r) {
		for (size_t i = 0; i < boundaryOps[r].size(); ++i) {
			for (size_t j = 0; j < boundaryOps[r][i].size(); ++j) {
				delete boundaryOps[r][i][j];
			}
		}
	}
}

Assembler::measurements Assembler::assemble(ActionOperator::Action action)
{
	Assembler::measurements times = {0.0, 0.0};
	if (action & ActionOperator::Action::FILL) {
		for (int t = 0; t < info::env::threads; ++t) {
			for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
				for (esint i = info::mesh->elements->eintervalsDistribution[d]; i < info::mesh->elements->eintervalsDistribution[d + 1]; ++i) {
					esint elements = info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin;
					times += instantiate(action, info::mesh->elements->eintervals[i].code, etype[i], elementOps[i], i, elements);
				}
			}
		}

		for (int t = 0; t < info::env::threads; ++t) {
			for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
				for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
					if (info::mesh->boundaryRegions[r]->dimension) {
						for (esint i = info::mesh->boundaryRegions[r]->eintervalsDistribution[d]; i < info::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1]; ++i) {
							size_t elements = info::mesh->boundaryRegions[r]->eintervals[i].end - info::mesh->boundaryRegions[r]->eintervals[i].begin;
							instantiate(action, info::mesh->boundaryRegions[r]->eintervals[i].code, btype[r][i], boundaryOps[r][i], i, elements);
						}
					} else {
						size_t elements = info::mesh->boundaryRegions[r]->nodes->datatarray().size(t);
						instantiate(action, static_cast<int>(Element::CODE::POINT1), btype[r][t], boundaryOps[r][t], r, elements);
					}
				}
			}
		}
	} else {
		#pragma omp parallel for reduction(+:times)
		for (int t = 0; t < info::env::threads; ++t) {
			for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
				for (esint i = info::mesh->elements->eintervalsDistribution[d]; i < info::mesh->elements->eintervalsDistribution[d + 1]; ++i) {
					esint elements = info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin;
					switch (info::ecf->loop) {
					case ECF::LOOP::INHERITANCE: times += instantiate          (action, info::mesh->elements->eintervals[i].code, etype[i], elementOps[i], i, elements); break;
					case ECF::LOOP::CONDITIONS : times += instantiateConditions(action, info::mesh->elements->eintervals[i].code, etype[i], elementOps[i], i, elements); break;
					case ECF::LOOP::MANUAL     : times += instantiateManual    (action, info::mesh->elements->eintervals[i].code, etype[i], elementOps[i], i, elements); break;
					}
				}
			}
		}
		#pragma omp parallel for
		for (int t = 0; t < info::env::threads; ++t) {
			for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
				for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
					if (info::mesh->boundaryRegions[r]->dimension) {
						for (esint i = info::mesh->boundaryRegions[r]->eintervalsDistribution[d]; i < info::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1]; ++i) {
							size_t elements = info::mesh->boundaryRegions[r]->eintervals[i].end - info::mesh->boundaryRegions[r]->eintervals[i].begin;
							instantiate(action, info::mesh->boundaryRegions[r]->eintervals[i].code, btype[r][i], boundaryOps[r][i], i, elements);
						}
					} else {
						size_t elements = info::mesh->boundaryRegions[r]->nodes->datatarray().size(t);
						instantiate(action, static_cast<int>(Element::CODE::POINT1), btype[r][t], boundaryOps[r][t], t, elements);
					}
				}
			}
		}
	}

	return times;
}

void Assembler::printElementVolume(std::vector<double> &volume)
{
	std::vector<double> sum(info::mesh->elementsRegions.size());
	for (size_t r = 1; r < info::mesh->elementsRegions.size(); ++r) {
		for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			if (info::mesh->elements->eintervals[i].region == (esint)r || (info::mesh->elements->eintervals[i].region == 0 && r == info::mesh->elementsRegions.size() - 1)) {
				sum[0] += volume[i];
				sum[r] += volume[i];
			}
		}
	}
	Communication::allReduce(sum, Communication::OP::SUM);

	eslog::info("  ELEMENT REGION VOLUME                                                                        \n");
	eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - \n");
	for (size_t r = 0; r < info::mesh->elementsRegions.size(); ++r) {
		eslog::info("     %30s :                                            %e   \n", info::mesh->elementsRegions[r]->name.c_str(), sum[r]);
	}
}

void Assembler::printBoundarySurface(std::vector<double> &surface)
{
	Communication::allReduce(surface, Communication::OP::SUM);

	eslog::info("\n  BOUDNARY REGION SURFACE                                                                      \n");
	eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - \n");
	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			eslog::info("     %30s :                                            %e   \n", info::mesh->boundaryRegions[r]->name.c_str(), surface[r]);
			info::mesh->boundaryRegions[r]->area = surface[r];
		}
	}
}

void Assembler::printMaterials(const std::map<std::string, std::string> &settings)
{
	for (auto reg = info::mesh->elementsRegions.begin() + 1; reg != info::mesh->elementsRegions.end(); ++reg) {
		auto ms = settings.find((*reg)->name);
		if (ms != settings.end()) {
			eslog::info("  %55s: %34s\n", (*reg)->name.c_str(), ms->second.c_str());
		} else {
			ms = settings.find("ALL_ELEMENTS");
			if (ms != settings.end()) {
				eslog::info("  %55s: %34s\n", (*reg)->name.c_str(), ms->second.c_str());
			} else {
				eslog::info("  %55s: %34s\n", (*reg)->name.c_str(), info::mesh->materials.front()->name.c_str());
			}
		}
	}
}

bool Assembler::checkExpression(const std::string &name, ECFExpression &expression)
{
	if (expression.evaluator == nullptr) {
		if (Exprtk::check(expression.value)) {
			eslog::warning("   %25s:  %62s \n", name.c_str(), "INVALID EXPRESSION");
			return false;
		}
		expression.evaluator = Evaluator::create(expression.value);
	}
	if (expression.evaluator->parameters.size()) {
		eslog::info("   %25s:  %62s \n", name.c_str(), expression.value.c_str());
	} else {
		eslog::info("   %25s:  %62g \n", name.c_str(), expression.evaluator->evaluate());
	}
	return true;
}

bool Assembler::checkElementParameter(const std::string &name, std::map<std::string, ECFExpression> &settings)
{
	if (settings.size() == 1 && StringCompare::caseInsensitiveEq(settings.begin()->first, "ALL_ELEMENTS")) {
		return checkExpression(name, settings.begin()->second);
	} else {
		eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");
		for (auto region = info::mesh->elementsRegions.crbegin(); region != info::mesh->elementsRegions.crend(); ++region) {
			auto it = settings.find((*region)->name);
			if (it != settings.end()) {
				if (!checkExpression(it->first, it->second)) { return false; }
			}
		}
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	}
	return true;
}

bool Assembler::checkElementParameter(const std::string &name, std::map<std::string, ECFExpressionVector> &settings)
{
	if (settings.size() == 1 && StringCompare::caseInsensitiveEq(settings.begin()->first, "ALL_ELEMENTS")) {
		switch (info::mesh->dimension) {
		case 2: if (!checkExpression(name + ".X", settings.begin()->second.x) || !checkExpression(name + ".Y", settings.begin()->second.y)) { return false; } break;
		case 3: if (!checkExpression(name + ".X", settings.begin()->second.x) || !checkExpression(name + ".Y", settings.begin()->second.y) || !checkExpression(name + ".Z", settings.begin()->second.z)) { return false; } break;
		}
	} else {
		eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");
		for (auto region = info::mesh->elementsRegions.crbegin(); region != info::mesh->elementsRegions.crend(); ++region) {
			auto it = settings.find((*region)->name);
			if (it != settings.end()) {
				switch (info::mesh->dimension) {
				case 2: if (!checkExpression(it->first + ".X", it->second.x) || !checkExpression(it->first + ".Y", it->second.y)) { return false; } break;
				case 3: if (!checkExpression(it->first + ".X", it->second.x) || !checkExpression(it->first + ".Y", it->second.y) || !checkExpression(it->first + ".Z", it->second.z)) { return false; } break;
				}
			}
		}
	}
	eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	return true;
}

bool Assembler::checkElementParameter(const std::string &name, std::map<std::string, ECFExpressionVector> &settings, int dim)
{
	if (settings.size() == 1 && StringCompare::caseInsensitiveEq(settings.begin()->first, "ALL_ELEMENTS")) {
		return checkExpression(name, settings.begin()->second.data[dim]);
	} else {
		eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");
		for (auto region = info::mesh->elementsRegions.crbegin(); region != info::mesh->elementsRegions.crend(); ++region) {
			auto it = settings.find((*region)->name);
			if (it != settings.end()) {
				return checkExpression(name, it->second.data[dim]);
			}
		}
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	}
	return true;
}

bool Assembler::checkBoundaryParameter(const std::string &name, std::map<std::string, ECFExpression> &settings)
{
	eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");
	for (auto region = info::mesh->boundaryRegions.crbegin(); region != info::mesh->boundaryRegions.crend(); ++region) {
		auto it = settings.find((*region)->name);
		if (it != settings.end()) {
			if (!checkExpression(it->first, it->second)) { return false; }
		}
	}
	eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	return true;
}

bool Assembler::checkBoundaryParameter(const std::string &name, std::map<std::string, ECFExpressionVector> &settings)
{
	eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");
	for (auto region = info::mesh->boundaryRegions.crbegin(); region != info::mesh->boundaryRegions.crend(); ++region) {
		auto it = settings.find((*region)->name);
		if (it != settings.end()) {
			switch (info::mesh->dimension) {
			case 2: if (!checkExpression(it->first + ".X", it->second.x) || !checkExpression(it->first + ".Y", it->second.y)) { return false; } break;
			case 3: if (!checkExpression(it->first + ".X", it->second.x) || !checkExpression(it->first + ".Y", it->second.y) || !checkExpression(it->first + ".Z", it->second.z)) { return false; } break;
			}
		}
	}
	eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	return true;
}

bool Assembler::checkBoundaryParameter(const std::string &name, std::map<std::string, ECFExpressionOptionalVector> &settings)
{
	eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");
	for (auto region = info::mesh->boundaryRegions.crbegin(); region != info::mesh->boundaryRegions.crend(); ++region) {
		auto it = settings.find((*region)->name);
		if (it != settings.end()) {
			switch (info::mesh->dimension) {
			case 2: if ((it->second.x.isset && !checkExpression(it->first + ".X", it->second.x)) || (it->second.y.isset && !checkExpression(it->first + ".Y", it->second.y))) { return false; } break;
			case 3: if ((it->second.x.isset && !checkExpression(it->first + ".X", it->second.x)) || (it->second.y.isset && !checkExpression(it->first + ".Y", it->second.y)) || (it->second.z.isset && !checkExpression(it->first + ".Z", it->second.z))) { return false; } break;
			}
		}
	}
	eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	return true;
}

Evaluator* Assembler::getEvaluator(size_t interval, std::map<std::string, ECFExpression> &settings)
{
	int region = info::mesh->elements->eintervals[interval].region;
	auto it = settings.find(info::mesh->elementsRegions[region]->name);
	if (it == settings.end()) {
		it = settings.find(info::mesh->elementsRegions[0]->name);
	}
	if (it != settings.end()) {
		return it->second.evaluator;
	}
	return nullptr;
}

Evaluator* Assembler::getEvaluator(size_t interval, std::map<std::string, ECFExpressionVector> &settings, int dim)
{
	int region = info::mesh->elements->eintervals[interval].region;
	auto it = settings.find(info::mesh->elementsRegions[region]->name);
	if (it == settings.end()) {
		it = settings.find(info::mesh->elementsRegions[0]->name);
	}
	if (it != settings.end()) {
		return it->second.data[dim].evaluator;
	}
	return nullptr;
}

