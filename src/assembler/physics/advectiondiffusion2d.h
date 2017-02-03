
#ifndef SRC_ASSEMBLER_PHYSICS_ADVECTIONDIFFUSION2D_H_
#define SRC_ASSEMBLER_PHYSICS_ADVECTIONDIFFUSION2D_H_

#include "physics2d.h"

namespace espreso {

enum class Property;
struct AdvectionDiffusion2DConfiguration;

struct NewAdvectionDiffusion2D: public Physics2D
{
	NewAdvectionDiffusion2D(Mesh *mesh, NewInstance *instance, const AdvectionDiffusion2DConfiguration &configuration);

	void prepareTotalFETI();

	void assembleStiffnessMatrix(size_t domain);
	void assembleStiffnessMatrix(Element *e, std::vector<eslocal> &DOFs, DenseMatrix &Ke, DenseMatrix &fe);

	void processElement(const Element *e, DenseMatrix &Ke, DenseMatrix &fe) const;
	void processFace(const Element *e, DenseMatrix &Ke, DenseMatrix &fe) const;
	void processEdge(const Element *e, DenseMatrix &Ke, DenseMatrix &fe) const;
	void processNode(const Element *e, DenseMatrix &Ke, DenseMatrix &fe) const;

	void analyticRegularization(size_t domain);

	void storeSolution(std::vector<std::vector<double> > &solution, store::ResultStore *store);

	const std::vector<Property>& pointDOFs() const
	{
		static std::vector<Property> pointDOFs = { Property::TEMPERATURE };
		return pointDOFs;
	}
	const std::vector<Property>& midPointDOFs() const
	{
		static std::vector<Property> midPointDOFs = { Property::TEMPERATURE };
		return midPointDOFs;
	}
	const std::vector<Property>& edgeDOFs() const
	{
		static std::vector<Property> edgeDOFs = { };
		return edgeDOFs;
	}
	const std::vector<Property>& faceDOFs() const
	{
		static std::vector<Property> faceDOFs = { };
		return faceDOFs;
	}
	const std::vector<Property>& elementDOFs() const
	{
		static std::vector<Property> elementDOFs = { };
		return elementDOFs;
	}

protected:
	const AdvectionDiffusion2DConfiguration &_configuration;
};

}



#endif /* SRC_ASSEMBLER_PHYSICS_ADVECTIONDIFFUSION2D_H_ */
