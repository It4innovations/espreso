
#include "filler.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/module/acoustic.h"
#include "analysis/assembler/module/heattransfer.h"

#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

template <size_t dofs, class Module>
void _add(Module &module, Matrix_Base<double> *A, ParameterData &parameter)
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		double *data = A->mapping.elements[interval].data;
		const esint *position = A->mapping.elements[interval].position;
		switch (A->shape) {
		case Matrix_Shape::FULL:  module.elementOps[interval].emplace_back(instantiate<typename Module::NGP, dofs, MatrixFullFiller >(interval, module.controller, parameter, data, position)); break;
		case Matrix_Shape::LOWER: module.elementOps[interval].emplace_back(instantiate<typename Module::NGP, dofs, MatrixLowerFiller>(interval, module.controller, parameter, data, position)); break;
		case Matrix_Shape::UPPER: module.elementOps[interval].emplace_back(instantiate<typename Module::NGP, dofs, MatrixUpperFiller>(interval, module.controller, parameter, data, position)); break;
		}
	}
}

void addFiller(AX_HeatTransfer &module)
{
	if (module.K != nullptr) {
		_add<1>(module, module.K, module.elements.stiffness);
	}
	if (module.M != nullptr) {
		_add<1>(module, module.M, module.elements.mass);
	}

	if (module.rhs != nullptr) {
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			if (info::mesh->boundaryRegions[r]->dimension && module.elements.boundary.rhs.regions[r].data != NULL) {
				for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
					module.boundaryOps[r][interval].emplace_back(instantiate<AX_HeatTransfer::NGP, 1, VectorFiller>(r, interval, module.controller, module.elements.boundary.rhs.regions[r], module.rhs->mapping.boundary[r][interval].data, module.rhs->mapping.boundary[r][interval].position));
				}
			}
		}
	}
}

void addFiller(AX_Acoustic &module)
{
	if (module.K != nullptr) {
		_add<1>(module, module.K, module.elements.stiffness);
	}
	if (module.M != nullptr) {
		_add<1>(module, module.M, module.elements.mass);
	}
	if (module.C != nullptr) {
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			if (info::mesh->boundaryRegions[r]->dimension && module.elements.boundary.mass.regions[r].data != NULL) {
				for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
					double *data = module.C->mapping.boundary[r][interval].data;
					const esint *position = module.C->mapping.boundary[r][interval].position;
					module.boundaryOps[r][interval].emplace_back(instantiate<AX_Acoustic::NGP, 1, MatrixFullFiller>(r, interval, module.controller, module.elements.boundary.mass.regions[r], data, position));
				}
			}
		}
	}

	if (module.re.rhs != nullptr) {
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			if (info::mesh->boundaryRegions[r]->dimension && module.elements.boundary.rhs.regions[r].data != NULL) {
				for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
					module.boundaryOps[r][interval].emplace_back(instantiate<AX_Acoustic::NGP, 1, VectorFiller>(r, interval, module.controller, module.elements.boundary.rhs.regions[r], module.re.rhs->mapping.boundary[r][interval].data, module.re.rhs->mapping.boundary[r][interval].position));
				}
			}
		}
	}
}

}
