
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_BASIS_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_BASIS_H_

#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"

#include "physics/assembler/operator.hpp"
#include "physics/assembler/modules/heattransfer.module.opt.h"

namespace espreso {

struct Basis: public ElementOperatorBuilder {
	Basis(): ElementOperatorBuilder("BASIC FUNCTIONS") {}

	bool build(HeatTransferModuleOpt &kernel) override;

	void apply(int interval)
	{

	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_BASIS_H_ */
