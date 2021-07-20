
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
		case Matrix_Shape::FULL:  module.actionOps[interval].emplace_back(instantiate<typename Module::NGP, dofs, MatrixFullFiller >(interval, parameter, data, position)); break;
		case Matrix_Shape::LOWER: module.actionOps[interval].emplace_back(instantiate<typename Module::NGP, dofs, MatrixLowerFiller>(interval, parameter, data, position)); break;
		case Matrix_Shape::UPPER: module.actionOps[interval].emplace_back(instantiate<typename Module::NGP, dofs, MatrixUpperFiller>(interval, parameter, data, position)); break;
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
		_add<1>(module, module.C, module.elements.damping);
	}
}

}
