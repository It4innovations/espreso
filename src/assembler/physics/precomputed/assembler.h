
#ifndef SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct PrecomputedPhysics: public Physics {

	virtual void assemble()
	{
		ESINFO(PROGRESS2) << "Assemble matrices K and RHS";
		K.resize(_mesh.parts());
		f.resize(_mesh.parts());
		for (size_t p = 0; p < _mesh.parts(); p++) {
			composeSubdomain(p);
			ESINFO(PROGRESS2) << Info::plain() << ".";
		}
		ESINFO(PROGRESS2);
	}

	virtual void save()
	{
		ESINFO(PROGRESS2) << "Save matrices K and RHS";
		for (size_t p = 0; p < _mesh.parts(); p++) {
			std::ofstream osK(Logging::prepareFile(p, "K").c_str());
			osK << K[p];
			osK.close();

			std::ofstream osF(Logging::prepareFile(p, "f").c_str());
			osF << f[p];
			osF.close();
		}
	}

	PrecomputedPhysics(const APIMesh &mesh, size_t DOFs, double *rhs, size_t rhs_size)
	: Physics(mesh, DOFs), _apimesh(mesh), rhs(rhs), rhs_size(rhs_size) {};
	virtual ~PrecomputedPhysics() {};

	std::vector<SparseMatrix> K, T; // T will be deleted
	std::vector<std::vector<double> > f;

protected:
	virtual void composeSubdomain(size_t subdomain) =0;

	const APIMesh &_apimesh;
	double *rhs;
	size_t rhs_size;
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_ASSEMBLER_H_ */
