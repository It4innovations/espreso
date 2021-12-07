
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
		case Matrix_Shape::FULL:  module.elementFiller[interval].emplace_back(instantiate<typename Module::NGP, dofs, MatrixFullFiller >(interval, module.controller, parameter, data, position)); break;
		case Matrix_Shape::LOWER: module.elementFiller[interval].emplace_back(instantiate<typename Module::NGP, dofs, MatrixLowerFiller>(interval, module.controller, parameter, data, position)); break;
		case Matrix_Shape::UPPER: module.elementFiller[interval].emplace_back(instantiate<typename Module::NGP, dofs, MatrixUpperFiller>(interval, module.controller, parameter, data, position)); break;
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
			if (info::mesh->boundaryRegions[r]->dimension && module.elements.boundary.rhs.isSet(r)) {
				for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
					module.boundaryFiller[r][interval].emplace_back(instantiate<AX_HeatTransfer::NGP, 1, VectorFiller>(r, interval, module.controller, module.elements.boundary.rhs.regions[r], module.rhs->mapping.boundary[r][interval].data, module.rhs->mapping.boundary[r][interval].position));
				}
			}
		}
		for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
			if (module.heatSource.gp.isSet(interval)) {
				module.elementFiller[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, 1, VectorFiller>(interval, module.controller, module.elements.rhs, module.rhs->mapping.elements[interval].data, module.rhs->mapping.elements[interval].position));
			}
		}
	}

	if (module.dirichlet != nullptr) {
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			if (module.temperature.node.isSet(r)) {
				for(size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
					module.boundaryFiller[r][t].emplace_back(instantiate<AX_HeatTransfer::NGP, 1, VectorSetter>(r, t, module.controller, module.temperature.node.regions[r], module.dirichlet->mapping.boundary[r][t].data, module.dirichlet->mapping.boundary[r][t].position));
					module.boundaryFiller[r][t].back()->isconst = false;
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
			if (info::mesh->boundaryRegions[r]->dimension && module.elements.boundary.mass.isSet(r)) {
				for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
					double *data = module.C->mapping.boundary[r][interval].data;
					const esint *position = module.C->mapping.boundary[r][interval].position;
					module.boundaryFiller[r][interval].emplace_back(instantiate<AX_Acoustic::NGP, 1, MatrixFullFiller>(r, interval, module.controller, module.elements.boundary.mass.regions[r], data, position));
				}
			}
		}
	}

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (module.monopoleSource.gp.isSet(interval)) {
			module.elementFiller[interval].emplace_back(instantiate<AX_Acoustic::NGP, 1, VectorFiller>(interval, module.controller, module.elements.monopole, module.re.rhs->mapping.elements[interval].data, module.re.rhs->mapping.elements[interval].position));
		}
		if (module.dipoleSource.gp.isSet(interval)) {
			module.elementFiller[interval].emplace_back(instantiate<AX_Acoustic::NGP, 1, VectorFiller>(interval, module.controller, module.elements.dipole, module.re.rhs->mapping.elements[interval].data, module.re.rhs->mapping.elements[interval].position));
		}
	}

	if (module.re.rhs != nullptr) {
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			if (info::mesh->boundaryRegions[r]->dimension && module.elements.boundary.rhs.isSet(r)) {
				for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
					module.boundaryFiller[r][interval].emplace_back(instantiate<AX_Acoustic::NGP, 1, VectorFiller>(r, interval, module.controller, module.elements.boundary.rhs.regions[r], module.re.rhs->mapping.boundary[r][interval].data, module.re.rhs->mapping.boundary[r][interval].position));
				}
			}
		}
	}

	if (module.re.dirichlet != nullptr) {
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			if (module.pressure.node.isSet(r)) {
				for(size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
					module.boundaryFiller[r][t].emplace_back(instantiate<AX_HeatTransfer::NGP, 1, VectorSetter>(r, t, module.controller, module.pressure.node.regions[r], module.re.dirichlet->mapping.boundary[r][t].data, module.re.dirichlet->mapping.boundary[r][t].position));
					module.boundaryFiller[r][t].back()->isconst = false;
				}
			}
		}
	}
}

}
