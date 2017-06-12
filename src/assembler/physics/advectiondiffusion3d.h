
#ifndef SRC_ASSEMBLER_PHYSICS_ADVECTIONDIFFUSION3D_H_
#define SRC_ASSEMBLER_PHYSICS_ADVECTIONDIFFUSION3D_H_

#include "physics3d.h"
#include "advectiondiffusion.h"

namespace espreso {

struct AdvectionDiffusion3DConfiguration;

struct NewAdvectionDiffusion3D: public AdvectionDiffusion, public Physics3D
{
	NewAdvectionDiffusion3D(Mesh *mesh, Instance *instance, const AdvectionDiffusion3DConfiguration &configuration);

	void prepareTotalFETI();

	void processElement(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processFace(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processEdge(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processNode(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processSolution(const Step &step);

protected:
	void assembleMaterialMatrix(const Step &step, const Element *e, eslocal node, double temp, DenseMatrix &K, DenseMatrix &CD) const;
	void postProcessElement(const Step &step, const Element *e, std::vector<Solution*> &solution);

	const AdvectionDiffusion3DConfiguration &_configuration;
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_ADVECTIONDIFFUSION3D_H_ */
