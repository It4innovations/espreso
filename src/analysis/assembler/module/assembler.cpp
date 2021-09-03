
#include "assembler.hpp"

#include "basis/evaluator/expressionevaluator.h"
#include "basis/expression/expression.h"
#include "basis/expression/variable.h"
#include "basis/utilities/parser.h"
#include "basis/utilities/utils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "analysis/assembler/operators/expression.h"

#include <algorithm>
#include <numeric>

#include "basis/utilities/print.h"

using namespace espreso;

Assembler::Assembler()
: elementOps(info::mesh->elements->eintervals.size()), elementRes(info::mesh->elements->eintervals.size()),
  boundaryOps(info::mesh->boundaryRegions.size()), boundaryRes(info::mesh->boundaryRegions.size()),
  version(0)
{
	for (size_t i = 0; i < info::mesh->boundaryRegions.size(); ++i) {
		if (info::mesh->boundaryRegions[i]->dimension) {
			boundaryOps[i].resize(info::mesh->boundaryRegions[i]->eintervals.size());
			boundaryRes[i].resize(info::mesh->boundaryRegions[i]->eintervals.size());
		}
	}
}

void Assembler::initDirichlet(std::map<std::string, ECFExpression> &settings, Vector_Sparse<double> &dirichlet)
{
	size_t dsize = 0;
	for (auto it = settings.begin(); it != settings.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		dsize += region->nodes->datatarray().size();
	}
	dirichletIndices.reserve(dsize);
	for (auto it = settings.begin(); it != settings.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		dirichletIndices.insert(dirichletIndices.end(), region->nodes->datatarray().begin(), region->nodes->datatarray().end());
	}
	dirichletPermutation.resize(dsize);
	std::iota(dirichletPermutation.begin(), dirichletPermutation.end(), 0);
	std::sort(dirichletPermutation.begin(), dirichletPermutation.end(), [&] (const esint &i, const esint &j) { return dirichletIndices[i] < dirichletIndices[j]; });
	dsize = 0;
	for (auto i = dirichletPermutation.begin(); i != dirichletPermutation.end(); ++i) {
		if (i == dirichletPermutation.begin() || dirichletIndices[*i] != dirichletIndices[*(i - 1)]) {
			++dsize;
		}
	}
	dirichlet.resize(info::mesh->nodes->IDs->datatarray().size(), dsize);
	auto dir = dirichlet.indices;
	for (auto i = dirichletPermutation.begin(); i != dirichletPermutation.end(); ++i) {
		if (i == dirichletPermutation.begin() || dirichletIndices[*i] != dirichletIndices[*(i - 1)]) {
			*dir++ = dirichletIndices[*i];
		}
	}
	std::sort(dirichletIndices.begin(), dirichletIndices.end());
}

void Assembler::fillDirichlet(std::map<std::string, ECFExpression> &settings, Vector_Sparse<double> &dirichlet)
{
	size_t offset = 0;
	std::vector<double> values(dirichletPermutation.size());
	for (auto it = settings.begin(); it != settings.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		it->second.evaluator->evalSelectedSparse(
				region->nodes->datatarray().size(),
				region->nodes->datatarray().data(),
				it->second.evaluator->params,
				values.data() + offset);
		offset += region->nodes->datatarray().size();
	}

	for (size_t i = 0, j = 0, v = 0; i < dirichletIndices.size(); i = j, ++v) {
		dirichlet.vals[v] = 0;
		while (j < dirichletIndices.size() && dirichletIndices[j] == dirichletIndices[i]) {
			dirichlet.vals[v] += values[dirichletPermutation[j++]];
		}
		dirichlet.vals[v] /= j - i;
	}
	dirichlet.touched = true;
}

void Assembler::iterate()
{
	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; ++t) {
		for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
			for (esint i = info::mesh->elements->eintervalsDistribution[d]; i < info::mesh->elements->eintervalsDistribution[d + 1]; ++i) {
				size_t elementsInInterval = info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin;

				for (size_t element = 0; element < elementsInInterval; ++element) {
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
					(**op).move(-elementsInInterval);
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
							(**op).move(-elementsInInterval);
						}
					}
				}
			}
		}
	}
}

void Assembler::printParameterStats(const char* name, ParameterData &parameter)
{
	printf("parameter[isconst/update]:  ");
	for (size_t i = 0; i < parameter.update.size(); ++i) {
		if (parameter.data) {
			printf(" [%c/%c]", parameter.isconst[i] ? 'C' : ' ', parameter.update[i] > 0 ? 'U' : ' ');
		} else {
			printf(" [-/-]");
		}
	}
	printf(" %s\n", name);
}

void Assembler::printParameterStats(const char* name, NamedData *data)
{
	printf("nameddata[isconst/update]:   [ /%c] %s\n", data->updated ? 'U' : ' ', name);
}

void Assembler::setMaterials(const std::map<std::string, std::string> &settings)
{
	for (auto ei = info::mesh->elements->eintervals.begin(); ei != info::mesh->elements->eintervals.end(); ++ei) {
		int region = ei->region;
		if (region == -1) { // intersected regions
			for (auto rindex = ei->regions.begin(); rindex != ei->regions.end(); ++rindex) {
				if (settings.find(info::mesh->elementsRegions[*rindex]->name) != settings.end()) {
					region = *rindex;
				}
			}
		}

		if (region == -1 || settings.find(info::mesh->elementsRegions[region]->name) == settings.end()) {
			region = 0;
		}
		auto mat = settings.find(info::mesh->elementsRegions[region]->name);
		if (mat == settings.end()) {
			eslog::error("Invalid material configuration: a region without a material settings found.\n");
		}

		for (size_t i = 0; i < info::mesh->materials.size(); ++i) {
			if (StringCompare::caseInsensitiveEq(info::mesh->materials[i]->name, mat->second)) {
				ei->material = i;
			}
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

bool Assembler::examineBoundaryParameter(const std::string &name, std::map<std::string, ConvectionConfiguration> &settings, ParametersConvection &convection)
{
	return false;
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
}

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
