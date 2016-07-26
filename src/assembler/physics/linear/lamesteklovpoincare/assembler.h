
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_LAMESTEKLOVPOINCARE_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_LAMESTEKLOVPOINCARE_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct LameSteklovPoincare: public LinearPhysics
{
	LameSteklovPoincare(const Mesh &mesh)
	: LinearPhysics(
			mesh,
			{ Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z },
			{ Property::FIXED_DISPLACEMENT_X, Property::FIXED_DISPLACEMENT_Y, Property::FIXED_DISPLACEMENT_Z },
			SparseMatrix::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {};
protected:
	void composeSubdomain(size_t subdomain);
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_LINEAR_LAMESTEKLOVPOINCARE_ASSEMBLER_H_ */
