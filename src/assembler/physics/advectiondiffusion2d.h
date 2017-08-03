
#ifndef SRC_ASSEMBLER_PHYSICS_ADVECTIONDIFFUSION2D_H_
#define SRC_ASSEMBLER_PHYSICS_ADVECTIONDIFFUSION2D_H_

#include "physics2d.h"
#include "advectiondiffusion.h"

namespace espreso {

enum class Property;
struct AdvectionDiffusion2DConfiguration;

struct AdvectionDiffusion2D: public AdvectionDiffusion, public Physics2D
{
	AdvectionDiffusion2D(Mesh *mesh, Instance *instance, const AdvectionDiffusion2DConfiguration &configuration);

	void prepare();

	virtual std::vector<std::pair<ElementType, Property> > propertiesToStore() const;

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
