
#ifndef SRC_ASSEMBLER_PHYSICS_HYPRE_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_HYPRE_ASSEMBLER_H_

#include "../assembler.h"

//#include "hypre.h"

namespace espreso {

struct HyprePhysics: public Physics {

	virtual bool singular() const
	{
		return true;
	}

	virtual void assemble()
	{
		ESINFO(PROGRESS2) << "Assemble matrices K, kernels, and RHS.";
		cilk_for (size_t p = 0; p < _mesh.parts(); p++) {
			composeSubdomain(p);
			ESINFO(PROGRESS2) << Info::plain() << ".";
		}
		ESINFO(PROGRESS2);
	}

	HyprePhysics(Mesh &mesh, SparseMatrix::MatrixType mtype): Physics(mesh, mtype, {}, {}, {}, {}, {}) {};
	virtual ~HyprePhysics() {};



protected:
	virtual void composeSubdomain(size_t subdomain) =0;
};

}



#endif /* SRC_ASSEMBLER_PHYSICS_HYPRE_ASSEMBLER_H_ */
