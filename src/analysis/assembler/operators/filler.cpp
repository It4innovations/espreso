
#include "filler.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/module/heattransfer.h"

#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

void addFiller(AX_HeatTransfer &module)
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (module.K != nullptr) {
			double *data = module.K->mapping.elements[interval].data;
			const esint *position = module.K->mapping.elements[interval].position;
			switch (module.K->shape) {
			case Matrix_Shape::FULL:  module.actionOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, 1, MatrixFullFiller >(interval, module.elements.stiffness, data, position)); break;
			case Matrix_Shape::LOWER: module.actionOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, 1, MatrixLowerFiller>(interval, module.elements.stiffness, data, position)); break;
			case Matrix_Shape::UPPER: module.actionOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, 1, MatrixUpperFiller>(interval, module.elements.stiffness, data, position)); break;
			}

		}
		if (module.M != nullptr) {
			double *data = module.M->mapping.elements[interval].data;
			const esint *position = module.M->mapping.elements[interval].position;
			switch (module.M->shape) {
			case Matrix_Shape::FULL:  module.actionOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, 1, MatrixFullFiller >(interval, module.elements.mass, data, position)); break;
			case Matrix_Shape::LOWER: module.actionOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, 1, MatrixLowerFiller>(interval, module.elements.mass, data, position)); break;
			case Matrix_Shape::UPPER: module.actionOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, 1, MatrixUpperFiller>(interval, module.elements.mass, data, position)); break;
			}
		}
	}
}

}
