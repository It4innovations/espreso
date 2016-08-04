
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_TEMPERATURE_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_TEMPERATURE_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct Temperature: public LinearPhysics
{
	Temperature(Mesh &mesh)
	: LinearPhysics(mesh, { Property::TEMPERATURE }, SparseMatrix::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {};

	void init();

	static std::vector<Property> elementDOFs;
	static std::vector<Property> faceDOFs;
	static std::vector<Property> edgeDOFs;
	static std::vector<Property> pointDOFs;
	static std::vector<Property> midPointDOFs;

protected:
	void composeSubdomain(size_t subdomain);
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_LINEAR_TEMPERATURE_ASSEMBLER_H_ */
