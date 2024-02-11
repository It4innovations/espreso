
#include "dualoperator.h"
#include "totalfeti.implicit.h"
#include "totalfeti.explicit.h"
#include "totalfeti.explicit.acc.h"

#include "esinfo/eslog.hpp"

#include <climits>

namespace espreso {

template<typename T>
DualOperator<T>* DualOperator<T>::set(FETI<T> &feti, const step::Step &step)
{
	DualOperator<T>* dual = nullptr;
	switch (feti.configuration.dual_operator) {
	case FETIConfiguration::DUAL_OPERATOR::IMPLICIT:
		eslog::info(" = DUAL OPERATOR                                                         IMPLICIT TOTAL FETI = \n");
		dual = new TotalFETIImplicit<T>(feti);
		break;
	case FETIConfiguration::DUAL_OPERATOR::EXPLICIT:
		eslog::info(" = DUAL OPERATOR                                                         EXPLICIT TOTAL FETI = \n");
		dual = new TotalFETIExplicit<T>(feti);
		break;
	case FETIConfiguration::DUAL_OPERATOR::EXPLICIT_GPU:
		if (DirectSparseSolver<T>::provideFactors()) {
			eslog::info(" = DUAL OPERATOR                                                  EXPLICIT TOTAL FETI ON GPU = \n");
			dual = new TotalFETIExplicitAcc<T,esint>(feti);
		} else {
			eslog::globalerror("Third party software problem: solver does not provide factors that are required for EXPLICIT TOTAL FETI ON GPU.\n");
		}
		break;
	}
	dual->set(step);
	return dual;
}


// y = pm.v0/norm(pm.v0);
// v = H(y);
// lambda = norm(v); if lambda<eps; error('Incorrect initial eigenvector v0!'); end;
// nit = 1; err = inf;

// while err>pm.tol && nit<pm.maxit
//   y = v/lambda;
//   v = H(y);
//   lambda_prev = lambda; lambda = norm(v);
//   err = abs(lambda-lambda_prev)/abs(lambda);
//   nit = nit+1;
// end

template<typename T>
void DualOperator<T>::estimateMaxEigenValue(double epsilon, int maxIterations, double &lambda, int &iterations)
{
	DualOperator<T> *F = feti.dualOperator;

	// {-1, 1} / norma
	Vector_Dual<T> y, v;
	double norm = 1 / std::sqrt(feti.sinfo.lambdasTotal);
	for (esint i = feti.lambdas.nhalo; i < y.size; ++i) {
		// y.vals[i] = (1 - 2 * ((feti.decomposition->halo + i) % 2)) * norm;
	}
	y.synchronize();
	F->apply(y, v);

	// double lambda = std::sqrt(v.dot());
	// double err = std::numeric_limits<T>::max();
	// for (int i = 0; i < maxIterations && epsilon < err; ++i) {
	// 	// y = v->scale(1 / lambda);
	// 	F->apply(y, v);
	// 	double _lambda = lambda;
	// 	// lambda = v->norm();
	// 	err = std::fabs(lambda - _lambda) / std::fabs(lambda);
	// }
}

template<typename T>
void DualOperator<T>::estimateMaxProjectedEigenValue(double epsilon, int maxIterations, double &lambda, int &iterations)
{
	DualOperator<double> *F = feti.dualOperator;
	Projector<double> *P = feti.projector;

	// {-1, 1} / norma
	Vector_Dual<T> y, v;
	// for (esint i = y.nhalo; i < y.size; ++i) {
	// 	y.vals[i] = 1; 
	// }
	// P->apply()
	F->apply(y, v);
	// P->apply()

	// double lambda = 0; //v->norm();
	// double err = 111; // DOUBLE_MAX;
	// for (int i = 0; i < maxIterations && epsilon < err; ++i) {
	// 	// y = v->scale(1 / lambda);
	// 	F->apply(y, v);
	// 	double _lambda = lambda;
	// 	// lambda = v->norm();
	// 	err = std::fabs(lambda - _lambda) / std::fabs(lambda);
	// }
}

template class DualOperator<double>;

}
