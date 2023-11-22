
#include "feti.h"
#include "dualoperator/dualoperator.h"
#include "iterativesolver/pcpg.h"
#include "projector/projector.h"
#include "preconditioner/preconditioner.h"

#include "esinfo/eslog.hpp"
#include "esinfo/systeminfo.h"
#include "wrappers/mpi/communication.h"

namespace espreso {

template <typename T>
FETI<T>::FETI(FETIConfiguration &configuration)
: configuration(configuration)
{

}

template <typename T>
FETI<T>::~FETI()
{
	delete iterativeSolver;
	delete projector;
	delete dualOperator;
	delete preconditioner;
}

template <typename T>
bool FETI<T>::set(const step::Step &step)
{
	double start = eslog::time();

	esint offset[2] = { 0, 0 };
	esint size[5] = { 0, 0, 0, 0, 0 };
	size[0] = K.domains.size();
	for (size_t d = 0; d < K.domains.size(); ++d) {
		sinfo.R1size = offset[0] = size[2] += regularization.R1.domains[d].ncols;
		sinfo.R2size = offset[1] = size[3] += regularization.R2.domains[d].ncols;
	}
	sinfo.lambdasLocal = equalityConstraints.global + equalityConstraints.paired + equalityConstraints.local + equalityConstraints.nn;
	size[4] = sinfo.lambdasLocal - equalityConstraints.nhalo;

	Communication::exscan(offset, NULL, 2, MPITools::getType<esint>().mpitype, MPI_SUM);
	Communication::allReduce(size, NULL, 5, MPITools::getType<esint>().mpitype, MPI_SUM);
	sinfo.domains = size[0];
	sinfo.R1totalSize = size[2];
	sinfo.R2totalSize = size[3];
	sinfo.lambdasTotal = size[4];
	sinfo.R1offset = info::mpi::rank ? offset[0] : 0;
	sinfo.R2offset = info::mpi::rank ? offset[1] : 0;

	Vector_Dual<T>::set(equalityConstraints.nhalo, sinfo.lambdasLocal, equalityConstraints.lmap, K.decomposition->neighbors);
	Vector_Kernel<T>::set(sinfo.R1offset, sinfo.R1size, sinfo.R1totalSize);

	eslog::checkpointln("FETI: SET INFO");

	iterativeSolver = IterativeSolver<T>::set(*this, step);
	projector = Projector<T>::set(*this, step);
	dualOperator = DualOperator<T>::set(*this, step);
	preconditioner = Preconditioner<T>::set(*this, step);

	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
	eslog::info(" = EXTERNAL LINEAR SOLVER %*s = \n", 66, DirectSolver<T, Matrix_CSR>::name());
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
	iterativeSolver->info();
	projector->info();
	dualOperator->info();
	preconditioner->info();

	eslog::info(" = FETI SOLVER SET                                                                %8.3f s = \n", eslog::time() - start);
	if (MPITools::node->rank == 0) {
		info::system::memory::solver = info::system::memoryAvail();
	}
	eslog::info(" = FETI SOLVER MEMORY FOOTPRINT [GB] %55.2f = \n", (info::system::memory::physics - info::system::memory::solver) / 1024. / 1024.);

	return true;
}

template <typename T>
bool FETI<T>::update(const step::Step &step)
{
	double start = eslog::time();
	projector->update(step);
	dualOperator->update(step);
	preconditioner->update(step);
	eslog::info("       = FETI SOLVER UPDATED                                                %8.3f s = \n", eslog::time() - start);
	return true;
}

template <typename T>
bool FETI<T>::solve(const step::Step &step)
{
	double start = eslog::time();
	IterativeSolverInfo info;
	iterativeSolver->solve(step, info);

	switch (info.error) {
	case IterativeSolverInfo::ERROR::OK: break;
	case IterativeSolverInfo::ERROR::STAGNATION:

		eslog::info("       = FETI SOLVER ERROR                                   NONDECREASING CONVERGENCE = \n");
		break;
	case IterativeSolverInfo::ERROR::MAX_ITERATIONS_REACHED:
		eslog::info("       = FETI SOLVER ERROR                                  MAXIMUM ITERATIONS REACHED = \n");
		break;
	case IterativeSolverInfo::ERROR::INVALID_DATA:
		eslog::info("       = FETI SOLVER ERROR                                          INVALID INPUT DATA = \n");
		break;
	case IterativeSolverInfo::ERROR::CONVERGENCE_ERROR:
		eslog::info("       = FETI SOLVER ERROR              SOLVER DOES NOT CONVERGE TO THE REQUESTED NORM = \n");
		break;
	}
	eslog::info("       = ITERATIONS TOTAL                                                    %9d = \n", info.iterations);
	eslog::info("       = FETI SOLVER TIME                                                   %8.3f s = \n", eslog::time() - start);
	eslog::info("       = ----------------------------------------------------------------------------- = \n");
	return true;
}

template struct FETI<double>;

}
