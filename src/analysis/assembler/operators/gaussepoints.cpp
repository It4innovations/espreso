
#include "gausspoints.h"

#include "analysis/assembler/module/heattransfer.h"

namespace espreso {

template <class Module>
void _moveEnodesToGPs(Module &module, const ParameterData &from, ParameterData &to, int dimension)
{
	to.addInput(from);
	to.resize();
	module.addParameter(to);

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		switch (dimension) {
		case 1: module.elementOps[interval].emplace_back(instantiate<typename Module::NGP, 1, FromNodesToGaussPoints>(interval, module.integration.N, from, to)); break;
		case 2: module.elementOps[interval].emplace_back(instantiate<typename Module::NGP, 2, FromNodesToGaussPoints>(interval, module.integration.N, from, to)); break;
		case 3: module.elementOps[interval].emplace_back(instantiate<typename Module::NGP, 3, FromNodesToGaussPoints>(interval, module.integration.N, from, to)); break;
		}
	}
}

void moveEnodesToGPs(AX_HeatTransfer &module, const ParameterData &from, ParameterData &to, int dimension)
{
	_moveEnodesToGPs(module, from, to, dimension);
}

}