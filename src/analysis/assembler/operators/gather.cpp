
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
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[interval].begin;
		bool cooToGPs = false;

		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
		if (mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN) {
			cooToGPs |= mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC;
		}

		if (cooToGPs) {
			elementOps[interval].push_back(generateElementOperator<CoordinatesToElementNodesAndGPs>(interval, etype[interval], procNodes));
		} else {
			elementOps[interval].push_back(generateElementOperator<CoordinatesToElementNodes>(interval, etype[interval], procNodes));
		}

		if (Results::gradient != nullptr || Results::flux != nullptr) {
			elementOps[interval].push_back(generateElementOperator<TemperatureToElementNodes>(interval, etype[interval], procNodes, Results::temperature->data.data()));
		}
	}
}
