
#ifndef SRC_ASSEMBLER_PHYSICS_HYPRE_ADVECTIONDIFFUSION2D_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_HYPRE_ADVECTIONDIFFUSION2D_ASSEMBLER_H_

#include "../assembler.h"
#include "LLNL_FEI_Impl.h"

namespace espreso {

struct HypreTemperature: public HyprePhysics
{
	enum class STABILIZATION {
		SUPG = 0,
		CAU = 1
	};

	bool uniformDOFs() const { return true; }

	HypreTemperature(Mesh &mesh, LLNL_FEI_Impl &feiPtr): HyprePhysics(mesh, SparseMatrix::MatrixType::REAL_UNSYMMETRIC), feiPtr(feiPtr) {};

	void init();

	static double sigma;
	static STABILIZATION stabilization;

protected:
	LLNL_FEI_Impl &feiPtr;

	void composeSubdomain(size_t subdomain);
};

}



#endif /* SRC_ASSEMBLER_PHYSICS_HYPRE_ADVECTIONDIFFUSION2D_ASSEMBLER_H_ */
