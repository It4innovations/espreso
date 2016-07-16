
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
		ESINFO(PROGRESS2) << "Save matrices K and RHS";
		for (size_t p = 0; p < K.size(); p++) {
			std::ofstream osK(Logging::prepareFile(p, "K").c_str());
			osK << K[p];
			osK.close();
		}

		for (size_t p = 0; p < f.size(); p++) {
			std::ofstream osF(Logging::prepareFile(p, "f").c_str());
			osF << f[p];
			osF.close();
		}

		for (size_t p = 0; p < R1.size(); p++) {
			std::ofstream osR(Logging::prepareFile(p, "R1").c_str());
			SparseMatrix tmpR = R1[p];
			tmpR.ConvertDenseToCSR(0);
			osR << tmpR;
			osR.close();
		}

    // AM uncomment for unsymmetric case
    for (size_t p = 0; p < R2.size(); p++) {
        std::ofstream osR(Logging::prepareFile(p, "R2").c_str());
        SparseMatrix tmpR = R2[p];
        tmpR.ConvertDenseToCSR(0);
        osR << tmpR;
        osR.close();
    }
//
//		for (size_t p = 0; p < R1H.size(); p++) {
//			std::ofstream osR(Logging::prepareFile(p, "R1H").c_str());
//			SparseMatrix tmpR = R1H[p];
//			tmpR.ConvertDenseToCSR(0);
//			osR << tmpR;
//			osR.close();
//		}
//
//		for (size_t p = 0; p < R2H.size(); p++) {
//			std::ofstream osR(Logging::prepareFile(p, "R2H").c_str());
//			SparseMatrix tmpR = R2H[p];
//			tmpR.ConvertDenseToCSR(0);
//			osR << tmpR;
//			osR.close();
//		}
//
//		for (size_t p = 0; p < RegMat.size(); p++) {
//			std::ofstream osR(Logging::prepareFile(p, "RegMat").c_str());
//			SparseMatrix tmpR = RegMat[p];
//			tmpR.ConvertDenseToCSR(0);
//			osR << tmpR;
//			osR.close();
//		}
	}

	size_t DOFs;

	SparseMatrix::MatrixType mtype;
	std::vector<SparseMatrix> K, T, R1, R2, R1H, R2H, RegMat; // T will be deleted
	std::vector<std::vector<double> > f;

	Physics(const Mesh &mesh, size_t DOFs, SparseMatrix::MatrixType mtype): _mesh(mesh), DOFs(DOFs), mtype(mtype) {};
	virtual ~Physics() {};

protected:
	const Mesh& _mesh;
};

}

#endif /* SRC_ASSEMBLER_PHYSICS_ASSEMBLER_H_ */
