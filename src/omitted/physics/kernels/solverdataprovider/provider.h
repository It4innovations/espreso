
#ifndef SRC_PHYSICS_KERNELS_SOLVERDATAPROVIDER_PROVIDER_H_
#define SRC_PHYSICS_KERNELS_SOLVERDATAPROVIDER_PROVIDER_H_

#include <utility>
#include <vector>

namespace espreso {

class MatrixCSRDistributed;
class VectorsDenseDistributed;
class MatrixCSRFETI;
class MatrixDenseFETI;
class VectorsDenseDistributed;
enum class MatrixType: int;

struct SolverDataProvider {
	struct General {
		virtual MatrixType getMatrixType() =0;

		virtual void dirichletIndices(std::vector<std::pair<esint, esint>> &indices) =0;
		virtual void dirichletValues(std::vector<double> &values) =0;

		virtual void inequalityIndices(std::vector<std::pair<esint, esint>> &indices) =0;
		virtual void inequalityNormals(std::vector<double> &values) =0;
		virtual void inequalityGaps(std::vector<double> &values) =0;

		virtual ~General() {}
	};

	struct FETI {
		virtual MatrixType getMatrixType(esint domain) =0;

		virtual int initKernels(MatrixCSRFETI &K, MatrixDenseFETI &N1, MatrixDenseFETI &N2, MatrixCSRFETI &RegMat, bool ortogonalizeCluster) =0;
		virtual void fillKernels(MatrixCSRFETI &K, MatrixCSRFETI &M, MatrixDenseFETI &N1, MatrixDenseFETI &N2, MatrixCSRFETI &RegMat, bool ortogonalizeCluster) =0;

		virtual ~FETI() {}
	};

	struct Hypre {
		virtual int numfnc() =0;
		virtual void initKernels(MatrixCSRDistributed &K, VectorsDenseDistributed &N) =0;
		virtual void fillKernels(MatrixCSRDistributed &K, VectorsDenseDistributed &N) =0;

		virtual ~Hypre() {}
	};

	template<class TGeneral, class TFETI, class THypre>
	SolverDataProvider(TGeneral *general, TFETI *feti, THypre *hypre)
	: general(general), feti(feti), hypre(hypre) { }

	~SolverDataProvider()
	{
		delete general;
		delete feti;
		delete hypre;
	}

	General *general;
	FETI *feti;
	Hypre *hypre;
};

}

#endif /* SRC_PHYSICS_KERNELS_SOLVERDATAPROVIDER_PROVIDER_H_ */
