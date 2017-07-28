
#ifndef SRC_ASSEMBLER_PHYSICS_LAMESTEKLOVPOINCARE3D_H_
#define SRC_ASSEMBLER_PHYSICS_LAMESTEKLOVPOINCARE3D_H_

#include "structuralmechanics3d.h"
#include "boundarybased3d.h"

namespace espreso {

struct LameSteklovPoincare3D: public BoundaryBased3D, public StructuralMechanics3D
{
	LameSteklovPoincare3D(Mesh *mesh, Instance *instance, const StructuralMechanics3DConfiguration &configuration);

	void prepareTotalFETI();
	void prepareHybridTotalFETIWithKernels();

	virtual void preprocessData(const Step &step);

	virtual void updateMatrix(const Step &step, Matrices matrices, size_t domain, const std::vector<Solution*> &solution);
	virtual void updateMatrix(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution);

	void processSolution(const Step &step);

protected:
	enum BEMSolutionIndex: size_t {
		BOUNDARY = 0,

		SIZE     = 1
	};

	static size_t BEMOffset;
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_LAMESTEKLOVPOINCARE3D_H_ */
