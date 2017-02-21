
#ifndef SRC_ASSEMBLER_PHYSICS_SHALLOWWATER2D_H_
#define SRC_ASSEMBLER_PHYSICS_SHALLOWWATER2D_H_

#include "physics2d.h"

namespace espreso {

enum class Property;
struct ShallowWater2DConfiguration;

struct ShallowWater2D: public Physics2D
{
	ShallowWater2D(Mesh *mesh, Instance *instance, const ShallowWater2DConfiguration &configuration);

	MatrixType getMatrixType(const Step &step, size_t domain) const;

	void prepareTotalFETI();

	void processElement(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processFace(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processEdge(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processNode(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processSolution(const Step &step);

	void analyticRegularization(size_t domain);

	const std::vector<Property>& pointDOFs() const
	{
		static std::vector<Property> pointDOFs = { Property::MOMENTUM_X, Property::MOMENTUN_Y, Property::PRESSURE };
		return pointDOFs;
	}
	const std::vector<Property>& midPointDOFs() const
	{
		static std::vector<Property> midPointDOFs = { Property::MOMENTUM_X, Property::MOMENTUN_Y, Property::PRESSURE };
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
	const ShallowWater2DConfiguration &_configuration;
};

}



#endif /* SRC_ASSEMBLER_PHYSICS_SHALLOWWATER2D_H_ */
