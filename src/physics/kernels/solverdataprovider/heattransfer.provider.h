
#ifndef SRC_PHYSICS_KERNELS_SOLVERDATAPROVIDER_HEATTRANSFER_PROVIDER_H_
#define SRC_PHYSICS_KERNELS_SOLVERDATAPROVIDER_HEATTRANSFER_PROVIDER_H_

#include "provider.h"

namespace espreso {

class HeatTransferLoadStepConfiguration;

struct HeatTransferSolverDataProvider: public SolverDataProvider {
	struct General: public SolverDataProvider::General {
		MatrixType getMatrixType();

		void dirichletIndices(std::vector<std::pair<esint, esint>> &indices);
		void dirichletValues(std::vector<double> &values);

		void inequalityIndices(std::vector<std::pair<esint, esint>> &indices);
		void inequalityNormals(std::vector<double> &values);
		void inequalityGaps(std::vector<double> &values);

		HeatTransferLoadStepConfiguration &_configuration;
		General(HeatTransferLoadStepConfiguration &configuration): _configuration(configuration) {}
	};

	struct FETI: public SolverDataProvider::FETI {
		MatrixType getMatrixType(esint domain);

		bool hasKernel(esint domain);
		void initKernels(MatrixCSRFETI &K, MatrixDenseFETI &N1, MatrixDenseFETI &N2, MatrixCSRFETI &RegMat, bool ortogonalizeCluster);
		void fillKernels(MatrixCSRFETI &K, MatrixDenseFETI &N1, MatrixDenseFETI &N2, MatrixCSRFETI &RegMat, bool ortogonalizeCluster);

		HeatTransferLoadStepConfiguration &_configuration;
		FETI(HeatTransferLoadStepConfiguration &configuration): _configuration(configuration) {}
	};

	struct Hypre: public SolverDataProvider::Hypre {
		int numfnc();
		void initKernels(MatrixCSRDistributed &K, VectorsDenseDistributed &N);
		void fillKernels(MatrixCSRDistributed &K, VectorsDenseDistributed &N);

		HeatTransferLoadStepConfiguration &_configuration;
		Hypre(HeatTransferLoadStepConfiguration &configuration): _configuration(configuration) {}
	};

	HeatTransferSolverDataProvider(HeatTransferLoadStepConfiguration &conf)
	: SolverDataProvider(new General(conf), new FETI(conf), new Hypre(conf)) { }
};

}

#endif /* SRC_PHYSICS_KERNELS_SOLVERDATAPROVIDER_HEATTRANSFER_PROVIDER_H_ */
