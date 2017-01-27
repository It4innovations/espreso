
#ifndef SRC_ASSEMBLER_PHYSICS_TRANSIENT_ELASTICITY_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_TRANSIENT_ELASTICITY_ASSEMBLER_H_

#include "../assembler.h"
#include "../../../../config/linearelasticity3d.h"

namespace espreso {

struct TransientElasticity: public TransientPhysics
{
	TransientElasticity(Mesh &mesh, Constraints &constraints, const ESPRESOSolver &configuration)
	: TransientPhysics(
			mesh, constraints, configuration,
			MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE,
			elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs)
	{
		beta = 0.25;
		gama = 0.5;
		timestep = 1e-3;
		deltaT = timestep;
		timeConstant = 1 / (beta * timestep * timestep);
		A = {
			  1.0 / (beta * deltaT * deltaT),
			 gama / (beta * deltaT),
			  1.0 / (beta * deltaT),

			  1.0 / (beta * 2) - 1.0,
			 gama /  beta      - 1.0,
			(gama / (beta * 2) - 1.0) * deltaT,

			(1.0 - gama) * deltaT,
			       gama  * deltaT
		};
	};

	void prepareMeshStructures();
	void assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const;
	void makeStiffnessMatricesRegular();
	void assembleB1() {};
	void assembleB0() {};

	void saveMeshProperties(store::Store &store);
	void saveMeshResults(store::Store &store, const std::vector<std::vector<double> > &results);

	static std::vector<Property> elementDOFs;
	static std::vector<Property> faceDOFs;
	static std::vector<Property> edgeDOFs;
	static std::vector<Property> pointDOFs;
	static std::vector<Property> midPointDOFs;

	double timeConstant;

	double beta;
	double gama;
	double timestep;
	double deltaT;

protected:
	void composeSubdomain(size_t subdomain);
};

}

#endif /* SRC_ASSEMBLER_PHYSICS_TRANSIENT_ELASTICITY_ASSEMBLER_H_ */
