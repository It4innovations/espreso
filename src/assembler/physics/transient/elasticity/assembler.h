
#ifndef SRC_ASSEMBLER_PHYSICS_TRANSIENT_ELASTICITY_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_TRANSIENT_ELASTICITY_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct TransientElasticity: public TransientPhysics
{
	TransientElasticity(const Mesh &mesh): TransientPhysics(mesh, { DOFType::DISPLACEMENT_X, DOFType::DISPLACEMENT_Y, DOFType::DISPLACEMENT_Z }, SparseMatrix::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE)
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
