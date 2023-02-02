
#include "assembler.hpp"

#include "analysis/assembler/operators/expression.h"
#include "basis/evaluator/expressionevaluator.h"
#include "basis/expression/expression.h"
#include "basis/expression/variable.h"
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

#include <algorithm>
#include <numeric>

using namespace espreso;

Assembler::Assembler(PhysicsConfiguration &settings)
: settings(settings),
  etype(info::mesh->elements->eintervals.size()), btype(info::mesh->boundaryRegions.size()),
  elementOps(info::mesh->elements->eintervals.size()), elementFiller(info::mesh->elements->eintervals.size()), elementRes(info::mesh->elements->eintervals.size()),
  boundaryOps(info::mesh->boundaryRegions.size()), boundaryFiller(info::mesh->boundaryRegions.size()), boundaryRes(info::mesh->boundaryRegions.size())
{
	for (size_t i = 0; i < info::mesh->boundaryRegions.size(); ++i) {
		if (info::mesh->boundaryRegions[i]->dimension) {
			btype[i].resize(info::mesh->boundaryRegions[i]->eintervals.size());
			boundaryOps[i].resize(info::mesh->boundaryRegions[i]->eintervals.size());
			boundaryFiller[i].resize(info::mesh->boundaryRegions[i]->eintervals.size());
			boundaryRes[i].resize(info::mesh->boundaryRegions[i]->eintervals.size());
		} else {
			btype[i].resize(info::mesh->boundaryRegions[i]->nodes->threads());
			boundaryOps[i].resize(info::mesh->boundaryRegions[i]->nodes->threads());
			boundaryFiller[i].resize(info::mesh->boundaryRegions[i]->nodes->threads());
			boundaryRes[i].resize(info::mesh->boundaryRegions[i]->nodes->threads());
		}
	}
}

double Assembler::assemble(ActionOperator::Action action)
{
	double time = 0;
	#pragma omp parallel for reduction(+:time)
	for (int t = 0; t < info::env::threads; ++t) {
		for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
			for (esint i = info::mesh->elements->eintervalsDistribution[d]; i < info::mesh->elements->eintervalsDistribution[d + 1]; ++i) {
				esint elements = info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin;
				time += instantiate(action, info::mesh->elements->eintervals[i].code, etype[i], elementOps[i], elements);
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
						instantiate(action, info::mesh->boundaryRegions[r]->eintervals[i].code, btype[r][i], boundaryOps[r][i], elements);
					}
				} else {
					size_t elements = info::mesh->boundaryRegions[r]->nodes->datatarray().size(t);
					instantiate(action, static_cast<int>(Element::CODE::POINT1), btype[r][t], boundaryOps[r][t], elements);
				}
			}
		}
	}

	return time / info::mesh->elements->eintervals.size();
}

void Assembler::iterate()
{
	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; ++t) {
		for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
			for (esint i = info::mesh->elements->eintervalsDistribution[d]; i < info::mesh->elements->eintervalsDistribution[d + 1]; ++i) {
				size_t elements = info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin;
				for (size_t element = 0; element < elements; ++element) {
					for (auto op = elementOps[i].begin(); op != elementOps[i].end(); ++op) {
						if((*op)->update) {
							if(element == 0 || !(*op)->isconst) {
								(**op)();
								++(**op);
							}
						}
					}
				}
				for (auto op = elementOps[i].begin(); op != elementOps[i].end(); ++op) {
					if((*op)->update) {
						if((*op)->isconst) {
							(*op)->move(-1);
						} else {
							(*op)->move(-(long)elements);
						}
					}
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
}

void Assembler::fill()
{
//	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
//		size_t elements = info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin;
//		std::vector<InsertOperator*> ops; ops.reserve(elementFiller[i].size());
//		for (auto op = elementFiller[i].cbegin(); op != elementFiller[i].cend(); ++op) {
//			ops.push_back(dynamic_cast<InsertOperator*>(*op));
//		}
//		if (settings.simd && elements >= SIMD::size) {
//			esint chunks = elements / SIMD::size;
//			for (esint c = 0; c < chunks; ++c) {
//				for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
//					(*op)->simd();
//				}
//			}
//			for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
//				(*op)->peel(elements % SIMD::size);
//			}
//		} else {
//			for (esint c = 0; c < elements; ++c) {
//				for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
//					(*op)->sisd();
//				}
//			}
//		}
//		for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
//			(*op)->move(-elements);
//		}
//	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
				size_t elementsInInterval = info::mesh->boundaryRegions[r]->eintervals[i].end - info::mesh->boundaryRegions[r]->eintervals[i].begin;

				for(size_t element = 0; element < elementsInInterval; ++element) {
					for (auto op = boundaryFiller[r][i].begin(); op != boundaryFiller[r][i].end(); ++op) {
						(**op)();
						++(**op);
					}
				}

				for (auto op = boundaryFiller[r][i].begin(); op != boundaryFiller[r][i].end(); ++op) {
					(*op)->move(-elementsInInterval);
				}
			}
		} else {
			for (int t = 0; t < info::env::threads; ++t) {
				for (auto n = info::mesh->boundaryRegions[r]->nodes->datatarray().begin(t); n != info::mesh->boundaryRegions[r]->nodes->datatarray().end(t); ++n) {
					for (auto op = boundaryFiller[r][t].begin(); op != boundaryFiller[r][t].end(); ++op) {
						(**op)();
						++(**op);
					}
				}
				for (auto op = boundaryFiller[r][t].begin(); op != boundaryFiller[r][t].end(); ++op) {
					(*op)->move(-info::mesh->boundaryRegions[r]->nodes->datatarray().size(t));
				}
			}
		}
	}
}

void Assembler::results()
{
	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; ++t) {
		for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
			for (esint i = info::mesh->elements->eintervalsDistribution[d]; i < info::mesh->elements->eintervalsDistribution[d + 1]; ++i) {
				size_t elementsInInterval = info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin;

				for (size_t element = 0; element < elementsInInterval; ++element) {
					for (auto op = elementRes[i].begin(); op != elementRes[i].end(); ++op) {
//						if((*op)->update) {
//							if(element == 0 || !(*op)->isconst) {
								(**op)();
								++(**op);
//							}
//						}
					}
				}
				for (auto op = elementRes[i].begin(); op != elementRes[i].end(); ++op) {
//					if((*op)->update) {
//						if((*op)->isconst) {
//							(*op)->move(-1);
//						} else {
							(*op)->move(-elementsInInterval);
//						}
//					}
				}
			}

			for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
				if (info::mesh->boundaryRegions[r]->dimension) {
					for (esint i = info::mesh->boundaryRegions[r]->eintervalsDistribution[d]; i < info::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1]; ++i) {
						size_t elementsInInterval = info::mesh->boundaryRegions[r]->eintervals[i].end - info::mesh->boundaryRegions[r]->eintervals[i].begin;

						for(size_t element = 0; element < elementsInInterval; ++element) {
							for (auto op = boundaryRes[r][i].begin(); op != boundaryRes[r][i].end(); ++op) {
//								if((*op)->update) {
//									if(element == 0 || !(*op)->isconst) {
										(**op)();
										++(**op);
//									}
//								}
							}
						}

						for (auto op = boundaryRes[r][i].begin(); op != boundaryRes[r][i].end(); ++op) {
//							if((*op)->update) {
//								if((*op)->isconst) {
//									(*op)->move(-1);
//								} else {
									(*op)->move(-elementsInInterval);
//								}
//							}
						}
					}
				}
			}
		}
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			if (info::mesh->boundaryRegions[r]->dimension == 0) {
				for (auto n = info::mesh->boundaryRegions[r]->nodes->datatarray().begin(t); n != info::mesh->boundaryRegions[r]->nodes->datatarray().end(t); ++n) {
					for (auto op = boundaryRes[r][t].begin(); op != boundaryRes[r][t].end(); ++op) {
//						if((*op)->update) {
//							if(n == info::mesh->boundaryRegions[r]->nodes->datatarray().begin(t) || !(*op)->isconst) {
								(**op)();
								++(**op);
//							}
//						}
					}
				}
				for (auto op = boundaryRes[r][t].begin(); op != boundaryRes[r][t].end(); ++op) {
//					if((*op)->update) {
//						if((*op)->isconst) {
//							(*op)->move(-1);
//						} else {
							(*op)->move(-info::mesh->boundaryRegions[r]->nodes->datatarray().size(t));
//						}
//					}
				}
			}
		}
	}
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

void Assembler::printParameterStats(ParameterData &parameter)
{
	printf("parameter [isconst/update]:  ");
	for (size_t i = 0; i < parameter.update.size(); ++i) {
		if (parameter.data) {
			printf(" [%c/%c]", parameter.isconst[i] ? 'C' : ' ', parameter.update[i] > 0 ? 'U' : ' ');
		} else {
			printf(" [-/-]");
		}
	}
	printf(" %s\n", parameter.name.c_str());
}

void Assembler::printParameterStats(NamedData *data)
{
	printf("nameddata [isconst/update]:   [ /%c] %s\n", data->updated ? 'U' : ' ', data->name.c_str());
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
		if (!Variable::create(expression)) {
			eslog::warning("   %18s:  %69s \n", name.c_str(), "INVALID EXPRESSION");
			return false;
		}
	}
	if (expression.evaluator->variables.size()) {
		eslog::info("   %18s:  %69s \n", name.c_str(), expression.evaluator->toString().c_str());
	} else {
		eslog::info("   %18s:  %69g \n", name.c_str(), expression.evaluator->eval(Evaluator::Params()));
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
		case 2: return checkExpression(name + ".X", settings.begin()->second.x) && checkExpression(name + ".Y", settings.begin()->second.y);
		case 3: return checkExpression(name + ".X", settings.begin()->second.x) && checkExpression(name + ".Y", settings.begin()->second.y) && checkExpression(name, settings.begin()->second.z);
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

bool Assembler::examineMaterialParameter(const std::string &material, const std::string &name, ECFExpression &settings, ExternalElementValue &externalValue, int dimension)
{
	if (settings.evaluator == nullptr) {
		if (!Variable::create(settings)) {
			eslog::warning("   %18s:  %69s \n", name.c_str(), "INVALID EXPRESSION");
			return false;
		}
	}
	if (settings.evaluator->variables.size()) {
		std::string params = Parser::join(", ", settings.evaluator->variables);
		eslog::info("   %18s:  %*s       FNC( %s )\n", name.c_str(), 55 - params.size(), "", params.c_str());
	} else {
		eslog::info("   %18s:  %69g \n", name.c_str(), settings.evaluator->eval(Evaluator::Params()));
	}

	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		if (StringCompare::caseInsensitiveEq(info::mesh->materials[info::mesh->elements->eintervals[i].material]->name, material)) {
			externalValue.evaluator[externalValue.dimension * i + dimension] = settings.evaluator;
		}
	}
	return true;
}

bool Assembler::examineElementParameter(const std::string &name, std::map<std::string, ECFExpression> &settings, ExternalElementValue &externalValue)
{
	return examineElementParameter<ECFExpression>(name, settings, externalValue, 0, [] (ECFExpression &expr) { return &expr; });
}
bool Assembler::examineElementParameter(const std::string &name, std::map<std::string, ECFExpressionVector> &settings, ExternalElementValue &externalValue, int dimension)
{
	return examineElementParameter<ECFExpressionVector>(name, settings, externalValue, dimension, [&] (ECFExpressionVector &expr) { return &expr.data[dimension]; });
}

bool Assembler::examineBoundaryParameter(const std::string &name, std::map<std::string, ECFExpression> &settings, ExternalBoundaryValue &externalValue)
{
	if (settings.size()) {
		eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");

		int rindex = 0;
		for (auto reg = info::mesh->boundaryRegions.begin(); reg != info::mesh->boundaryRegions.end(); ++reg, ++rindex) {
			auto ms = settings.find((*reg)->name);
			if (ms != settings.end()) {
				if (!Variable::create(ms->second, rindex)) {
					eslog::warning("   %30s:  %57s \n", (*reg)->name.c_str(), "INVALID EXPRESSION");
					return false;
				}
				Evaluator *evaluator = ms->second.evaluator;
				externalValue.evaluator[externalValue.dimension * rindex] = evaluator;
				if (evaluator->variables.size()) {
					std::string params = Parser::join(", ", evaluator->variables);
					eslog::info("   %30s:  %*s       FNC( %s )\n", (*reg)->name.c_str(), 43 - params.size(), "", params.c_str());
				} else {
					eslog::info("   %30s:  %57g \n", (*reg)->name.c_str(), evaluator->eval(Evaluator::Params()));
				}
			}
		}
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	}
	return true;
}

bool Assembler::examineBoundaryParameter(const std::string &name, std::map<std::string, ECFExpressionVector> &settings, ExternalBoundaryValue &value, int dimension)
{
	return examineBoundaryParameter<ECFExpressionVector>(name, settings, value, dimension, [&] (ECFExpressionVector &expr) { return &expr.data[dimension]; });
}

bool Assembler::examineBoundaryParameter(const std::string &name, std::map<std::string, ECFExpressionOptionalVector> &settings, ExternalBoundaryValue &value, int dimension)
{
	return examineBoundaryParameter<ECFExpressionOptionalVector>(name, settings, value, dimension, [&] (ECFExpressionOptionalVector &expr) { return &expr.data[dimension]; });
}


//bool Assembler::examineBoundaryParameter(const std::string &name, std::map<std::string, ConvectionConfiguration> &settings, ParametersConvection &convection)
//{
//	return false;
//	if (settings.size()) {
//		eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");
//
//		int rindex = 0;
//		for (auto reg = info::mesh->boundaryRegions.begin(); reg != info::mesh->boundaryRegions.end(); ++reg, ++rindex) {
//			auto ms = settings.find((*reg)->name);
//			if (ms != settings.end()) {
//
//				auto name = [] (const std::string &name) {
//					eslog::info("   %30s:  %*s -- %s \n", "", 53 - name.size(), "", name.c_str());
//				};
//
//				auto addParam = [&] (ParametersConvection::ExternalParameter &param, const Evaluator *evaluator, const std::string &name) {
//					param.gp.builder->insert(rindex, 0, evaluator);
//					if (evaluator->variables.size()) {
//						std::string params = Parser::join(", ", evaluator->variables);
//						eslog::info("   %30s: %s %*s       FNC( %s )\n", "", name.c_str(), 43 - params.size() - name.size(), "", params.c_str());
//					} else {
//						eslog::info("   %30s: %s %*g \n", "", name.c_str(), 57 - name.size(), evaluator->eval(Evaluator::Params()));
//					}
//				};
//
//				auto setMaterial = [&] () {
//					switch (ms->second.fluid) {
//					case ConvectionConfiguration::FLUID::AIR:
//					case ConvectionConfiguration::FLUID::STEAM:
//						name("AIR");
//						addParam(convection.absolutePressure, ms->second.absolute_pressure.evaluator, "ABSOLUTE PRESSURE");
//						break;
//
//					case ConvectionConfiguration::FLUID::WATER:
//					case ConvectionConfiguration::FLUID::ENGINE_OIL:
//					case ConvectionConfiguration::FLUID::TRANSFORMER_OIL:
//						break;
//					}
//				};
//
//				convection.configuration.regions[rindex].isset = true;
//				convection.configuration.regions[rindex].settings.front() = &ms->second;
//				switch (ms->second.type) {
//
//				case ConvectionConfiguration::TYPE::USER: {
//					eslog::info("   %30s:  %57s \n", ms->first.c_str(), "USER");
//					addParam(convection.heatTransferCoeficient, ms->second.heat_transfer_coefficient.evaluator, "HTC");
//					addParam(convection.externalTemperature, ms->second.external_temperature.evaluator, "EXTERNAL TEMPERATURE");
//				} break;
//
//				case ConvectionConfiguration::TYPE::EXTERNAL_NATURAL:
//					eslog::info("   %30s:  %57s \n", ms->first.c_str(), "EXTERNAL NATURAL");
//					addParam(convection.externalTemperature, ms->second.external_temperature.evaluator, "EXTERNAL TEMPERATURE");
//					setMaterial();
//
//					switch (ms->second.variant) {
//					case ConvectionConfiguration::VARIANT::INCLINED_WALL:
//						name("INCLINED WALL");
//						addParam(convection.wallHeight, ms->second.wall_height.evaluator, "WALL HEIGHT");
//						addParam(convection.wallHeight, ms->second.tilt_angle.evaluator, "TILT ANGLE");
//						break;
//					case ConvectionConfiguration::VARIANT::VERTICAL_WALL:
//						name("VERTICAL WALL");
//						addParam(convection.wallHeight, ms->second.wall_height.evaluator, "WALL HEIGHT");
//						break;
//					case ConvectionConfiguration::VARIANT::HORIZONTAL_PLATE_UP:
//						name("HORIZONTAL PLATE UP");
//						addParam(convection.length, ms->second.length.evaluator, "LENGTH");
//						break;
//					case ConvectionConfiguration::VARIANT::HORIZONTAL_PLATE_DOWN:
//						name("HORIZONTAL PLATE DOWN");
//						addParam(convection.length, ms->second.length.evaluator, "LENGTH");
//						break;
//					case ConvectionConfiguration::VARIANT::HORIZONTAL_CYLINDER:
//						name("HORIZONTAL CYLINDER");
//						addParam(convection.diameter, ms->second.diameter.evaluator, "DIAMETER");
//						break;
//					case ConvectionConfiguration::VARIANT::SPHERE:
//						name("SPHERE");
//						addParam(convection.diameter, ms->second.diameter.evaluator, "DIAMETER");
//						break;
//					default:
//						break;
//					}
//					break;
//
//				case ConvectionConfiguration::TYPE::INTERNAL_NATURAL:
//					eslog::info("   %30s:  %57s \n", ms->first.c_str(), "INTERNAL NATURAL");
//					addParam(convection.length, ms->second.length.evaluator, "LENGTH");
//					setMaterial();
//					switch (ms->second.variant) {
//					case ConvectionConfiguration::VARIANT::PARALLEL_PLATES:
//						name("PARALLEL PLATES");
//						addParam(convection.wallHeight, ms->second.wall_height.evaluator, "WALL HEIGHT");
//						addParam(convection.length, ms->second.length.evaluator, "LENGTH");
//						break;
//					case ConvectionConfiguration::VARIANT::CIRCULAR_TUBE:
//						name("CIRCULAR TUBE");
//						addParam(convection.wallHeight, ms->second.wall_height.evaluator, "WALL HEIGHT");
//						addParam(convection.diameter, ms->second.diameter.evaluator, "DIAMETER");
//						break;
//					default:
//						break;
//					}
//					break;
//
//				case ConvectionConfiguration::TYPE::EXTERNAL_FORCED:
//					setMaterial();
//					switch (ms->second.variant) {
//					case ConvectionConfiguration::VARIANT::AVERAGE_PLATE:
//						name("AVERAGE PLATE");
//						addParam(convection.externalTemperature, ms->second.external_temperature.evaluator, "EXTERNAL TEMPERATURE");
//						addParam(convection.length, ms->second.length.evaluator, "LENGTH");
//						break;
//					default:
//						break;
//					}
//					break;
//
//				case ConvectionConfiguration::TYPE::INTERNAL_FORCED:
//					switch (ms->second.variant) {
//					case ConvectionConfiguration::VARIANT::TUBE:
//						setMaterial();
//						name("TUBE");
//						addParam(convection.externalTemperature, ms->second.external_temperature.evaluator, "EXTERNAL TEMPERATURE");
//						addParam(convection.fluidVelocity, ms->second.fluid_velocity.evaluator, "FLUID VELOCITY");
//						addParam(convection.diameter, ms->second.diameter.evaluator, "DIAMETER");
//						break;
//					case ConvectionConfiguration::VARIANT::QUENCH:
//						name("QUENCH");
//						addParam(convection.absolutePressure, ms->second.absolute_pressure.evaluator, "ABSOLUTE PRESSURE");
//						addParam(convection.externalTemperature, ms->second.external_temperature.evaluator, "EXTERNAL TEMPERATURE");
//						addParam(convection.experimentalConstant, ms->second.experimental_constant.evaluator, "EXPERIMENTAL CONSTANT");
//						break;
//					case ConvectionConfiguration::VARIANT::QUENCH_PARALLEL:
//						name("QUENCH PARALLEL");
//						addParam(convection.absolutePressure, ms->second.absolute_pressure.evaluator, "ABSOLUTE PRESSURE");
//						addParam(convection.externalTemperature, ms->second.external_temperature.evaluator, "EXTERNAL TEMPERATURE");
//						addParam(convection.experimentalConstant, ms->second.experimental_constant.evaluator, "EXPERIMENTAL CONSTANT");
//						break;
//					default:
//						break;
//					}
//					break;
//				}
//			}
//		}
//		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
//	}
//}

bool Assembler::examineBoundaryParameter(const std::string &name, std::map<std::string, ImpedanceConfiguration> &settings, ExternalBoundaryValue &impedance)
{
	if (settings.size()) {
		eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");

		int rindex = 0;
		for (auto reg = info::mesh->boundaryRegions.begin(); reg != info::mesh->boundaryRegions.end(); ++reg, ++rindex) {
			auto ms = settings.find((*reg)->name);
			if (ms != settings.end()) {
				if (!Variable::create(ms->second.impedance, rindex)) {
					eslog::warning("   %30s:  %57s \n", (*reg)->name.c_str(), "INVALID EXPRESSION");
					return false;
				}
				Evaluator *evaluator = ms->second.impedance.evaluator;
				impedance.evaluator[impedance.dimension * rindex] = evaluator;
				if (evaluator->variables.size()) {
					std::string params = Parser::join(", ", evaluator->variables);
					eslog::info("   %30s:  %*s       FNC( %s )\n", (*reg)->name.c_str(), 43 - params.size(), "", params.c_str());
				} else {
					eslog::info("   %30s:  %57g \n", (*reg)->name.c_str(), evaluator->eval(Evaluator::Params()));
				}
			}
		}
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	}

	return true;
}

bool Assembler::examineBoundaryParameter(const std::string &name, std::map<std::string, PointSourceConfiguration> &settings, ExternalBoundaryValue &point_source)
{
	if (settings.size()) {
		eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");

		int rindex = 0;
		for (auto reg = info::mesh->boundaryRegions.begin(); reg != info::mesh->boundaryRegions.end(); ++reg, ++rindex) {
			auto ms = settings.find((*reg)->name);
			if (ms != settings.end()) {
				if (!Variable::create(ms->second.source_amplitude, rindex)) {
					eslog::warning("   %30s:  %57s \n", (*reg)->name.c_str(), "INVALID EXPRESSION");
					return false;
				}
				Evaluator *evaluator = ms->second.source_amplitude.evaluator;
				point_source.evaluator[point_source.dimension * rindex] = evaluator;
				if (evaluator->variables.size()) {
					std::string params = Parser::join(", ", evaluator->variables);
					eslog::info("   %30s:  %*s       FNC( %s )\n", (*reg)->name.c_str(), 43 - params.size(), "", params.c_str());
				} else {
					eslog::info("   %30s:  %57g \n", (*reg)->name.c_str(), evaluator->eval(Evaluator::Params()));
				}
			}
		}
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	}

	return true;
}
