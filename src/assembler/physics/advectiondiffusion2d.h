
#ifndef SRC_ASSEMBLER_PHYSICS_ADVECTIONDIFFUSION2D_H_
#define SRC_ASSEMBLER_PHYSICS_ADVECTIONDIFFUSION2D_H_

#include "physics2d.h"
#include "advectiondiffusion.h"

namespace espreso {

enum class Property;
struct AdvectionDiffusion2DConfiguration;

struct NewAdvectionDiffusion2D: public AdvectionDiffusion, public Physics2D
{
	NewAdvectionDiffusion2D(Mesh *mesh, Instance *instance, const AdvectionDiffusion2DConfiguration &configuration);

	void prepareTotalFETI();

	virtual std::vector<std::pair<ElementType, Property> > properties() const;

	void processElement(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processFace(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processEdge(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processNode(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processSolution(const Step &step);

protected:
	void assembleMaterialMatrix(const Step &step, const Element *e, eslocal node, double temp, DenseMatrix &K, DenseMatrix &CD, bool tangentCorrection) const;
	void postProcessElement(const Step &step, const Element *e, std::vector<Solution*> &solution);

	const AdvectionDiffusion2DConfiguration &_configuration;
};

}



#endif /* SRC_ASSEMBLER_PHYSICS_ADVECTIONDIFFUSION2D_H_ */
