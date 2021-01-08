
#include "kernel.opt.h"
#include "physics/assembler/operators/expression.h"

#include "basis/utilities/parser.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/elementsregionstore.h"

using namespace espreso;

template<typename Ttype>
void KernelOpt::validateRegionSettings(const std::string &name, const std::map<std::string, Ttype> &settings)
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
void KernelOpt::examineElementParameter(const std::string &name, const std::map<std::string, TSecond> &settings, ExpressionsToElements &builder, int dimension, std::function<const Evaluator*(const TSecond &expr)> getevaluator)
{
	builder.ecfname = name;
	if (settings.size() == 1 && StringCompare::caseInsensitiveEq(settings.begin()->first, "ALL_ELEMENTS")) {
		const Evaluator *evaluator = getevaluator(settings.begin()->second);
		for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			builder.evaluators[builder.dimension * i + dimension] = evaluator;
		}
		if (evaluator->variables.size()) {
			std::string params = Parser::join(", ", evaluator->variables);
			eslog::info("  %s:  %*s       FNC( %s )\n", name.c_str(), 54 - params.size(), "", params.c_str());
		} else {
			eslog::info("  %s:  %*g \n", name.c_str(), 88 - name.size(), evaluator->eval(Evaluator::Params()));
		}
	} else {
		if (settings.size() == 0) {
			eslog::info("  %s:  %*g \n", name.c_str(), 88 - name.size(), builder.defaultValue);
		} else {
			eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");
			int rindex = 1, rlast = info::mesh->elementsRegions.size() - 1;
			for (auto reg = info::mesh->elementsRegions.begin() + 1; reg != info::mesh->elementsRegions.end(); ++reg, ++rindex) {
				auto ms = settings.find((*reg)->name);
				if (ms == settings.end()) {
					ms = settings.find("ALL_ELEMENTS");
				}
				if (ms != settings.end()) {
					const Evaluator *evaluator = getevaluator(ms->second);
					if (evaluator->variables.size()) {
						std::string params = Parser::join(", ", evaluator->variables);
						eslog::info("   %30s:  %*s       FNC( %s )\n", (*reg)->name.c_str(), 43 - params.size(), "", params.c_str());
					} else {
						eslog::info("   %30s:  %57g \n", (*reg)->name.c_str(), evaluator->eval(Evaluator::Params()));
					}
					for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
						if (info::mesh->elements->eintervals[i].region == rindex || (info::mesh->elements->eintervals[i].region == 0 && rindex == rlast)) {
							builder.evaluators[builder.dimension * i + dimension] = evaluator;
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
								builder.evaluators[builder.dimension * i + dimension] = evaluator;
							}
						}
					}
				} else {
					eslog::info("   %30s:  %57g \n", (*reg)->name.c_str(), builder.defaultValue);
				}
			}
			eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		}
	}
}
