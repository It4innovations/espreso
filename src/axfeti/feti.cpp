
#include "feti.h"
#include "iterativesolver/pcpg.h"
#include "preconditioner/emptypreconditioner.h"
#include "projector/orthogonal.h"
#include "dualoperator/totalfeti.h"

#include "math2/generalization/matrix_feti.h"
#include "math2/generalization/vector_feti.h"

namespace espreso {

template <typename T>
static void setInfo(AX_FETI<T> *feti)
{
	esint offset[2] = { 0, 0 };
	esint size[4] = { 0, 0, 0, 0 };
	size[0] = feti->K->domains.size();
	for (size_t d = 0; d < feti->K->domains.size(); ++d) {
		offset[0] = size[2] += feti->regularization->R1.domains[d].ncols;
		offset[1] = size[3] += feti->regularization->R2.domains[d].ncols;
	}

	Communication::exscan(offset, NULL, 2, MPITools::getType<esint>().mpitype, MPI_SUM);
	Communication::allReduce(size, NULL, 4, MPITools::getType<esint>().mpitype, MPI_SUM);
	feti->sinfo.domains = size[0];
	feti->sinfo.R1size = size[2];
	feti->sinfo.R2size = size[3];
	feti->sinfo.R1offset = info::mpi::rank ? offset[0] : 0;
	feti->sinfo.R2offset = info::mpi::rank ? offset[1] : 0;
}

template <typename T>
void _info(const AX_FETI<T> *feti)
{
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
	feti->dualOperator->info();
	feti->iterativeSolver->info();
	feti->preconditioner->info();
	feti->projector->info();
}

template <typename T>
bool _set(AX_FETI<T> *feti, const step::Step &step, const Matrix_FETI<Matrix_CSR, T> &K, const typename AX_FETI<T>::Regularization &regularization, const typename AX_FETI<T>::EqualityConstraints &equalityConstraints)
{
	feti->step = &step;
	feti->K = &K;
	feti->regularization = &regularization;
	feti->equalityConstraints = &equalityConstraints;

	setInfo<T>(feti);

	feti->dualOperator = DualOperator<T>::set(feti);
	feti->iterativeSolver = IterativeSolver<T>::set(feti);
	feti->preconditioner = Preconditioner<T>::set(feti);
	feti->projector = Projector<T>::set(feti);

	_info(feti);

	return true;
}

template <typename T>
bool _update(AX_FETI<T> *feti, const step::Step &step, const Matrix_FETI<Matrix_CSR, T> &K)
{
	feti->dualOperator->update();
	feti->iterativeSolver->update();
	feti->preconditioner->update();
	feti->projector->update();
	return true;
}

template <typename T>
bool _solve(AX_FETI<T> *feti, const step::Step &step, const Vector_FETI<Vector_Dense, T> &f, Vector_FETI<Vector_Dense, T> &x)
{
	printf("FETI SOLVE\n");
	for (size_t d = 0; d < x.domains.size(); ++d) {
		x.set(0);
	}
	return true;
}

template <> bool AX_FETI<double>::set(const step::Step &step, Matrix_FETI<Matrix_CSR, double> &K, const Regularization &regularization, const EqualityConstraints &equalityConstraints) { return _set(this, step, K, regularization, equalityConstraints); }
template <> bool AX_FETI<std::complex<double> >::set(const step::Step &step, Matrix_FETI<Matrix_CSR, std::complex<double> > &K, const Regularization &regularization, const EqualityConstraints &equalityConstraints) { return _set(this, step, K, regularization, equalityConstraints); }

template <> bool AX_FETI<double>::update(const step::Step &step, Matrix_FETI<Matrix_CSR, double> &K) { return _update(this, step, K); }
template <> bool AX_FETI<std::complex<double> >::update(const step::Step &step, Matrix_FETI<Matrix_CSR, std::complex<double> > &K) { return _update(this, step, K); }

template <> bool AX_FETI<double>::solve(const step::Step &step, const Vector_FETI<Vector_Dense, double> &f, Vector_FETI<Vector_Dense, double> &x) { return _solve(this, step, f, x); }
template <> bool AX_FETI<std::complex<double> >::solve(const step::Step &step, const Vector_FETI<Vector_Dense, std::complex<double> > &f, Vector_FETI<Vector_Dense, std::complex<double> > &x) { return _solve(this, step, f, x); }

}
