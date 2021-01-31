
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_BASIS_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_BASIS_H_

#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"

#include "physics/assembler/operator.hpp"
#include "physics/assembler/kernels/heattransfer.kernel.opt.h"

namespace espreso {

struct Basis: public ElementOperatorBuilder {
	GET_NAME(Basis)

	bool build(HeatTransferKernelOpt &kernel) override;

	void apply(int interval)
	{

	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_BASIS_H_ */
