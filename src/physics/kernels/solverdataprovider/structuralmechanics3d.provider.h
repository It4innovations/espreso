
#ifndef SRC_PHYSICS_KERNELS_SOLVERDATAPROVIDER_STRUCTURALMECHANICS3D_PROVIDER_H_
#define SRC_PHYSICS_KERNELS_SOLVERDATAPROVIDER_STRUCTURALMECHANICS3D_PROVIDER_H_

#include "provider.h"
#include "basis/containers/point.h"

#include <cstddef>

namespace espreso {

class StructuralMechanicsLoadStepConfiguration;

struct StructuralMechanics3DSolverDataProvider: public SolverDataProvider {
	struct General: public SolverDataProvider::General {
		MatrixType getMatrixType();

		void dirichletIndices(std::vector<std::pair<esint, esint>> &indices);
		void dirichletValues(std::vector<double> &values);

		void inequalityIndices(std::vector<std::pair<esint, esint>> &indices);
		void inequalityNormals(std::vector<double> &values);
		void inequalityGaps(std::vector<double> &values);

		StructuralMechanicsLoadStepConfiguration &_configuration;
		General(StructuralMechanicsLoadStepConfiguration &configuration): _configuration(configuration) {}
	};

	struct FETI: public SolverDataProvider::FETI {
		MatrixType getMatrixType(esint domain);

		bool hasKernel(esint domain);
		int initKernels(MatrixCSRFETI &K, MatrixDenseFETI &N1, MatrixDenseFETI &N2, MatrixCSRFETI &RegMat, bool ortogonalizeCluster);
		void fillKernels(MatrixCSRFETI &K, MatrixCSRFETI &M, MatrixDenseFETI &N1, MatrixDenseFETI &N2, MatrixCSRFETI &RegMat, bool ortogonalizeCluster);

		StructuralMechanicsLoadStepConfiguration &_configuration;
		FETI(StructuralMechanicsLoadStepConfiguration &configuration): _configuration(configuration), _RegMat(NULL) {}
		~FETI();

		// Data for computation of analytic kernels
		std::vector<std::vector<esint> > dnodes;
		std::vector<Point> _cCenter, _cNorm;
		std::vector<double> _cr44, _cr45, _cr46, _cr55, _cr56;
		std::vector<size_t> _cNp;

		std::vector<Point> _dCenter, _dNorm;
		std::vector<double> _dr44, _dr45, _dr46, _dr55, _dr56;

		MatrixCSRFETI *_RegMat;
	};

	struct Hypre: public SolverDataProvider::Hypre {
		int numfnc();
		void initKernels(MatrixCSRDistributed &K, VectorsDenseDistributed &N);
		void fillKernels(MatrixCSRDistributed &K, VectorsDenseDistributed &N);

		StructuralMechanicsLoadStepConfiguration &_configuration;
		Hypre(StructuralMechanicsLoadStepConfiguration &configuration): _configuration(configuration) {}
	};

	StructuralMechanics3DSolverDataProvider(StructuralMechanicsLoadStepConfiguration &conf)
	: SolverDataProvider(new General(conf), new FETI(conf), new Hypre(conf)) { }
};

}

#endif /* SRC_PHYSICS_KERNELS_SOLVERDATAPROVIDER_STRUCTURALMECHANICS3D_PROVIDER_H_ */
