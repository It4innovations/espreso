
#include "feti.h"
#include "dualoperator/dualoperator.h"
#include "iterativesolver/pcpg.h"
#include "projector/projector.h"
#include "preconditioner/preconditioner.h"

#include "esinfo/eslog.hpp"

namespace espreso {

template <typename T>
static void setInfo(FETI<T> *feti)
{
	esint offset[2] = { 0, 0 };
	esint size[5] = { 0, 0, 0, 0, 0 };
	size[0] = feti->K->domains.size();
	for (size_t d = 0; d < feti->K->domains.size(); ++d) {
		feti->sinfo.R1size = offset[0] = size[2] += feti->regularization->R1.domains[d].ncols;
		feti->sinfo.R2size = offset[1] = size[3] += feti->regularization->R2.domains[d].ncols;
	}
	feti->sinfo.lambdasLocal = feti->equalityConstraints->global + feti->equalityConstraints->paired + feti->equalityConstraints->local + feti->equalityConstraints->nn;
	size[4] = feti->sinfo.lambdasLocal - feti->equalityConstraints->nhalo;

	Communication::exscan(offset, NULL, 2, MPITools::getType<esint>().mpitype, MPI_SUM);
	Communication::allReduce(size, NULL, 5, MPITools::getType<esint>().mpitype, MPI_SUM);
	feti->sinfo.domains = size[0];
	feti->sinfo.R1totalSize = size[2];
	feti->sinfo.R2totalSize = size[3];
	feti->sinfo.lambdasTotal = size[4];
	feti->sinfo.R1offset = info::mpi::rank ? offset[0] : 0;
	feti->sinfo.R2offset = info::mpi::rank ? offset[1] : 0;
}

template <typename T>
void _info(const FETI<T> *feti)
{
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
	eslog::info(" = EXTERNAL LINEAR SOLVER %*s = \n", 66, math::sparseSolver());
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
	feti->iterativeSolver->info();
	feti->projector->info();
	feti->dualOperator->info();
	feti->preconditioner->info();
}

template <typename T>
bool _set(FETI<T> *feti, const step::Step &step, Matrix_FETI<Matrix_CSR, T> &K, typename FETI<T>::Regularization &regularization, typename FETI<T>::EqualityConstraints &equalityConstraints)
{
	double start = eslog::time();
	feti->step = &step;
	feti->K = &K;
	feti->regularization = &regularization;
	feti->equalityConstraints = &equalityConstraints;

	setInfo<T>(feti);
	Vector_Dual<T>::set(feti->equalityConstraints->nhalo, feti->sinfo.lambdasLocal, equalityConstraints.lmap, K.decomposition->neighbors);
	Vector_Kernel<T>::set(feti->sinfo.R1offset, feti->sinfo.R1size, feti->sinfo.R1totalSize);

	eslog::checkpointln("FETI: SET INFO");

	feti->iterativeSolver = IterativeSolver<T>::set(feti);
	feti->projector = Projector<T>::set(feti);
	feti->dualOperator = DualOperator<T>::set(feti);
	feti->preconditioner = Preconditioner<T>::set(feti);

	_info(feti);

	eslog::info(" = FETI SOLVER SET                                                                %8.3f s = \n", eslog::time() - start);
	return true;
}

template <typename T>
bool _update(FETI<T> *feti, const step::Step &step, Matrix_FETI<Matrix_CSR, T> &K, Vector_FETI<Vector_Dense, T> &f)
{
	double start = eslog::time();
	feti->f = &f;

	feti->projector->update();
	feti->dualOperator->update();
	feti->preconditioner->update();

	eslog::info("       = FETI SOLVER UPDATED                                                %8.3f s = \n", eslog::time() - start);
	return true;
}

template <typename T>
bool _solve(FETI<T> *feti, const step::Step &step, Vector_FETI<Vector_Dense, T> &x)
{
	double start = eslog::time();
	feti->x = &x;

	IterativeSolverInfo info;
	feti->iterativeSolver->solve(info);

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

template <typename T>
static void _free(FETI<T> *feti)
{
	delete feti->iterativeSolver;
	delete feti->projector;
	delete feti->dualOperator;
	delete feti->preconditioner;
}

template <> bool FETI<double>::set(const step::Step &step, Matrix_FETI<Matrix_CSR, double> &K, Regularization &regularization, EqualityConstraints &equalityConstraints) { return _set(this, step, K, regularization, equalityConstraints); }
template <> bool FETI<std::complex<double> >::set(const step::Step &step, Matrix_FETI<Matrix_CSR, std::complex<double> > &K, Regularization &regularization, EqualityConstraints &equalityConstraints) { return _set(this, step, K, regularization, equalityConstraints); }

template <> bool FETI<double>::update(const step::Step &step, Matrix_FETI<Matrix_CSR, double> &K, Vector_FETI<Vector_Dense, double> &f) { return _update(this, step, K, f); }
template <> bool FETI<std::complex<double> >::update(const step::Step &step, Matrix_FETI<Matrix_CSR, std::complex<double> > &K, Vector_FETI<Vector_Dense, std::complex<double> > &f) { return _update(this, step, K, f); }

template <> bool FETI<double>::solve(const step::Step &step, Vector_FETI<Vector_Dense, double> &x) { return _solve(this, step, x); }
template <> bool FETI<std::complex<double> >::solve(const step::Step &step, Vector_FETI<Vector_Dense, std::complex<double> > &x) { return _solve(this, step, x); }

template <> FETI<double>::~FETI() { _free(this); }
template <> FETI<std::complex<double> >::~FETI() { _free(this); }

}
