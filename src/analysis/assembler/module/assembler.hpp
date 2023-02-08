
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

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
double Assembler::loop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, esint elements)
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

