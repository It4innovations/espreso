
#include "coordinates.h"
#include "temperature.h"
#include "analysis/assembler/module/heattransfer.h"
#include "analysis/assembler/module/heattransfer.generator.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"

using namespace espreso;

void HeatTransfer::gatherInputs()
{
	for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[i].begin;
		bool cooToGPs = false;

		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[i].material];
		if (mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN) {
			cooToGPs |= mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC;
		}

		if (cooToGPs) {
			elementOps[i].push_back(generateElementOperator<CoordinatesToElementNodesAndGPs>(i, etype[i], procNodes));
		} else {
			elementOps[i].push_back(generateElementOperator<CoordinatesToElementNodes>(i, etype[i], procNodes));
		}

		if (Results::gradient != nullptr || Results::flux != nullptr) {
			elementOps[i].push_back(generateElementOperator<TemperatureToElementNodes>(i, etype[i], procNodes, Results::temperature->data.data()));
		}
	}

	for(size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (bfilter[r]) {
			if (info::mesh->boundaryRegions[r]->dimension) {
				for(size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
					auto procNodes = info::mesh->boundaryRegions[r]->elements->cbegin() + info::mesh->boundaryRegions[r]->eintervals[i].begin;
					boundaryOps[r][i].push_back(generateBoundaryOperator<CoordinatesToElementNodes>(r, i, procNodes));
				}
			}
		}
	}
}
