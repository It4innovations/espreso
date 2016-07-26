
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
		ESINFO(PROGRESS2) << "Assemble matrices K, kernels, and InitialCondition";

		K.resize(_mesh.parts());
		R1.resize(_mesh.parts());
		R2.resize(_mesh.parts());
		RegMat.resize(_mesh.parts());
		f.resize(_mesh.parts());
		cilk_for (size_t p = 0; p < _mesh.parts(); p++) {
			composeSubdomain(p);
			K[p].mtype = mtype;
			ESINFO(PROGRESS2) << Info::plain() << ".";
		}
		ESINFO(PROGRESS2);
	}

	LinearPhysics(const Mesh &mesh, const std::vector<Property> unknowns, const std::vector<Property> dirichlet, SparseMatrix::MatrixType mtype)
	: Physics(mesh, unknowns, dirichlet, mtype) {};
	virtual ~LinearPhysics() {};

protected:
	virtual void composeSubdomain(size_t subdomain) =0;
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_LINEAR_ASSEMBLER_H_ */
