
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct LinearPhysics: public Physics {

	virtual bool singular() const
	{
		return true;
	}

	virtual void assemble()
	{
		ESINFO(PROGRESS2) << "Assemble matrices K, kernels, and RHS";

		K.resize(_mesh.parts());
		R1.resize(_mesh.parts());
		R2.resize(_mesh.parts());
		RegMat.resize(_mesh.parts());
		f.resize(_mesh.parts());
		cilk_for (size_t p = 0; p < _mesh.parts(); p++) {
			composeSubdomain(p);
			K[p].mtype = mtype;

			const std::vector<BoundaryCondition*> &bc = _mesh.boundaryConditions();
			for (size_t i = 0; i < bc.size(); i++) {
				bc[i]->fillForces(DOFs, _mesh, f[p], p);
			}

			ESINFO(PROGRESS2) << Info::plain() << ".";
		}
		ESINFO(PROGRESS2);
	}

	LinearPhysics(const Mesh &mesh, std::vector<DOFType> DOFs, SparseMatrix::MatrixType mtype): Physics(mesh, DOFs, mtype) {};
	virtual ~LinearPhysics() {};

protected:
	virtual void composeSubdomain(size_t subdomain) =0;
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_LINEAR_ASSEMBLER_H_ */
