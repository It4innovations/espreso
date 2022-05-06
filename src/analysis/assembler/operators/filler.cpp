
#include "filler.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/module/acoustic.h"
#include "analysis/assembler/module/heattransfer.h"
#include "analysis/scheme/steadystate.h"
#include "analysis/scheme/harmonic.h"

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

void addFiller(HeatTransfer &module, SteadyState &scheme)
{
	_add<1>(module, scheme.K, module.elements.stiffness);

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension && module.elements.boundary.rhs.isSet(r)) {
			for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
				module.boundaryFiller[r][interval].emplace_back(instantiate<HeatTransfer::NGP, 1, VectorFiller>(r, interval, module.controller, module.elements.boundary.rhs.regions[r], scheme.f->mapping.boundary[r][interval].data, scheme.f->mapping.boundary[r][interval].position));
			}
		}
	}
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (module.heatSource.gp.isSet(interval)) {
			module.elementFiller[interval].emplace_back(instantiate<HeatTransfer::NGP, 1, VectorFiller>(interval, module.controller, module.elements.rhs, scheme.f->mapping.elements[interval].data, scheme.f->mapping.elements[interval].position));
		}
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (module.temperature.node.isSet(r)) {
			for(size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
				module.boundaryFiller[r][t].emplace_back(instantiate<HeatTransfer::NGP, 1, VectorSetter>(r, t, module.controller, module.temperature.node.regions[r], scheme.dirichlet->mapping.boundary[r][t].data, scheme.dirichlet->mapping.boundary[r][t].position));
				module.boundaryFiller[r][t].back()->isconst = false;
			}
		}
	}
}

void addFiller(Acoustic &module, Harmonic &scheme)
{
	_add<1>(module, scheme.K, module.elements.stiffness);
	_add<1>(module, scheme.M, module.elements.mass);
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension && module.elements.boundary.mass.isSet(r)) {
			for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
				double *data = scheme.C->mapping.boundary[r][interval].data;
				const esint *position = scheme.C->mapping.boundary[r][interval].position;
				module.boundaryFiller[r][interval].emplace_back(instantiate<Acoustic::NGP, 1, MatrixFullFiller>(r, interval, module.controller, module.elements.boundary.mass.regions[r], data, position));
			}
		}
	}

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (module.monopoleSource.gp.isSet(interval)) {
			module.elementFiller[interval].emplace_back(instantiate<Acoustic::NGP, 1, VectorFiller>(interval, module.controller, module.elements.monopole, scheme.re.f->mapping.elements[interval].data, scheme.re.f->mapping.elements[interval].position));
		}
		if (module.dipoleSource.gp.isSet(interval)) {
			module.elementFiller[interval].emplace_back(instantiate<Acoustic::NGP, 1, VectorFiller>(interval, module.controller, module.elements.dipole, scheme.re.f->mapping.elements[interval].data, scheme.re.f->mapping.elements[interval].position));
		}
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension && module.elements.boundary.rhs.isSet(r)) {
			for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
				module.boundaryFiller[r][interval].emplace_back(instantiate<Acoustic::NGP, 1, VectorFiller>(r, interval, module.controller, module.elements.boundary.rhs.regions[r], scheme.re.f->mapping.boundary[r][interval].data, scheme.re.f->mapping.boundary[r][interval].position));
			}
		}

		if (info::mesh->boundaryRegions[r]->dimension && module.acceleration.gp.isSet(r)) {
			for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
				module.boundaryFiller[r][interval].emplace_back(instantiate<Acoustic::NGP, 1, VectorFiller>(r, interval, module.controller, module.proj_acceleration.gp.regions[r], scheme.re.f->mapping.boundary[r][interval].data, scheme.re.f->mapping.boundary[r][interval].position));
			}
		}
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (module.pressure.node.isSet(r)) {
			for(size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
				module.boundaryFiller[r][t].emplace_back(instantiate<HeatTransfer::NGP, 1, VectorSetter>(r, t, module.controller, module.pressure.node.regions[r], scheme.re.dirichlet->mapping.boundary[r][t].data, scheme.re.dirichlet->mapping.boundary[r][t].position));
				module.boundaryFiller[r][t].back()->isconst = false;
			}
		}
	}
}

}
