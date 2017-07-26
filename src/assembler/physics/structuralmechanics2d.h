
#ifndef SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS2D_H_
#define SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS2D_H_

#include "physics2d.h"
#include "structuralmechanics.h"

namespace espreso {

struct StructuralMechanics2DConfiguration;

struct StructuralMechanics2D: public StructuralMechanics, public Physics2D
{
	StructuralMechanics2D(Mesh *mesh, Instance *instance, const StructuralMechanics2DConfiguration &configuration);

	virtual std::vector<std::pair<ElementType, Property> > properties() const;

	void prepareTotalFETI();
	void analyticRegularization(size_t domain);

	void processElement(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processFace(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processEdge(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processNode(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processSolution(const Step &step);

	const std::vector<Property>& pointDOFs() const
	{
		static std::vector<Property> pointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y };
		return pointDOFs;
	}
	const std::vector<Property>& midPointDOFs() const
	{
		static std::vector<Property> midPointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y };
		return midPointDOFs;
	}

protected:
	void assembleMaterialMatrix(const Step &step, const Element *e, eslocal node, double temp, DenseMatrix &K) const;
	void postProcessElement(const Step &step, const Element *e, std::vector<Solution*> &solution);

	const StructuralMechanics2DConfiguration &_configuration;
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS2D_H_ */
