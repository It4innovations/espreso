
#ifndef SRC_ASSEMBLER_PHYSICS_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_ASSEMBLER_H_

#include "esmesh.h"
#include "essolver.h"

namespace espreso {

struct Physics {

	virtual bool singular() const =0;
	virtual void assemble() =0;
	virtual void save()
	{
		ESINFO(PROGRESS2) << "Save matrices K and RHS.";
		for (size_t p = 0; p < K.size(); p++) {
			std::ofstream osK(Logging::prepareFile(p, "Kreg").c_str());
			osK << K[p];
			osK.close();
		}

		for (size_t p = 0; p < f.size(); p++) {
			std::ofstream osF(Logging::prepareFile(p, "f").c_str());
			osF << f[p];
			osF.close();
		}

		for (size_t p = 0; p < R1.size(); p++) {
			std::ofstream osR(Logging::prepareFile(p, "R1_").c_str());
			SparseMatrix tmpR = R1[p];
			tmpR.ConvertDenseToCSR(0);
			osR << tmpR;
			osR.close();
		}

		for (size_t p = 0; p < R2.size(); p++) {
			std::ofstream osR(Logging::prepareFile(p, "R2_").c_str());
			SparseMatrix tmpR = R2[p];
			tmpR.ConvertDenseToCSR(0);
			osR << tmpR;
			osR.close();
		}
	}

	std::vector<Property> DOFs;

	SparseMatrix::MatrixType mtype;
	std::vector<SparseMatrix> K, R1, R2, RegMat;
	std::vector<std::vector<double> > f;

	Physics(const Mesh &mesh, const std::vector<Property> DOFs, SparseMatrix::MatrixType mtype)
	: _mesh(mesh), DOFs(DOFs), mtype(mtype) {};
	virtual ~Physics() {};

protected:
	const Mesh& _mesh;
};

}

#endif /* SRC_ASSEMBLER_PHYSICS_ASSEMBLER_H_ */
