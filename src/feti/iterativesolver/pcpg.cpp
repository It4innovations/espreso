
#include "pcpg.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "feti/projector/projector.h"
#include "feti/dualoperator/dualoperator.h"
#include "feti/preconditioner/preconditioner.h"

namespace espreso {

template class PCPG<double>;
template class PCPG<std::complex<double> >;

// initialization
// l_0: L => Gt * inv(GGt) * e
// r_0: L => d - F * lambda_0
// w_0: L => r_0 - Gt * inv(GGt) * G * r_0 :: (I - Q) * r_0
// y_0: L -> r_0 - Gt * inv(GGt) * G * S * w_0 :: (I - Q) * S * w_0
// p_0: L => w_0

// loop
// gamma_k: 1 => (y_k,w_k) / (p_k, F * p_k)
//   x_k+1: L => x_k + gama_k * p_k
//   r_k+1: L => r_k - gama_k * F * p_k
//   w_k+1: L => r_k+1 - Gt * inv(GGt) * G * r_k+1 :: (I - Q) * r_k+1
//   y_k+1: L => (I - Q) * S * w_k+1
//  beta_k: 1 => (y_k+1,w_k+1) / (p_k,w_k)
//   p_k+1: L => y_k+1 + beta_k * p_k

template <typename T>
void PCPG<T>::info()
{
	eslog::info(" = PRECONDITIONED CONJUGATE PROJECTED GRADIENT SETTINGS                                      = \n");
	switch (feti.configuration.stopping_criterion) {
	case FETIConfiguration::STOPPING_CRITERION::RELATIVE:
		eslog::info(" =   STOPPING CRITERION                                                             RELATIVE = \n");
		break;
	case FETIConfiguration::STOPPING_CRITERION::ABSOLUTE:
		eslog::info(" =   STOPPING CRITERION                                                             ABSOLUTE = \n");
		break;
	case FETIConfiguration::STOPPING_CRITERION::ARIOLI:
		eslog::info(" =   STOPPING CRITERION                                                               ARIOLI = \n");
		break;
	}
	eslog::info(" =   PRECISION                                                                      %.2e = \n", feti.configuration.precision);
	if (feti.configuration.max_iterations == 0) {
		eslog::info(" =   MAX_ITERATIONS                                                                     AUTO = \n");
	} else {
		eslog::info(" =   MAX_ITERATIONS                                                                  %7d = \n", feti.configuration.max_iterations);
	}
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

template <typename T>
PCPG<T>::PCPG(FETI<T> &feti)
: IterativeSolver<T>(feti)
{
	l.resize();
	r.resize();
	w.resize();
	y.resize();
	z.resize();
	p.resize();

	x.resize();
	Fp.resize();
}

template <> void PCPG<double>::solve(const step::Step &step, IterativeSolverInfo &info)
{
	DualOperator<double> *F = feti.dualOperator;
	Projector<double> *P = feti.projector;
	Preconditioner<double> *S = feti.preconditioner;

	P->applyGtInvGGt(P->e, l);             // l = Gt * inv(GGt) * e

	F->apply(l, r);                        // r = d - F * l
	r.scale(-1);                           //
	r.add(1, F->d);                        //

	P->apply(r, w);                        // w = P * r
	S->apply(w, z);                        // z = S * w
	P->apply(z, y);                        // y = P * z (y = P * S * w)

	y.copyTo(p);                           // p = w
	l.copyTo(x);                           // x = l

	double yw = y.dot(w);
	setInfo(info, feti.configuration, yw);

	eslog::checkpointln("FETI: CPG INITIALIZATION");
	eslog::startln("PCPG: ITERATIONS STARTED", "pcpg");
	while (!info.converged) {
		// gamma = (y, w) / (p, F * p)
		F->apply(p, Fp);
		eslog::accumulatedln("pcpg: apply F");
		double pFp = p.dot(Fp), gamma = yw / pFp;
		eslog::accumulatedln("pcpg: dot(p, Fp)");

		// x = x + gamma * p
		// r = r - gamma * F * p
		x.add(gamma, p);
		r.add(-gamma, Fp);
		eslog::accumulatedln("pcpg: update x, r");

		// w = P * r
		P->apply(r, w);
		eslog::accumulatedln("pcpg: apply P * r");

		// z = S * w
		S->apply(w, z);
		eslog::accumulatedln("pcpg: apply S * w");

		// y = P * z
		P->apply(z, y);
		eslog::accumulatedln("pcpg: apply P * z");

		// beta = (y+1, w+1) / (y, w)
		double _yw = y.dot(w), beta = _yw / yw;
		eslog::accumulatedln("pcpg: dot(y, w)");

		// p = y + beta * p  (y is not used anymore)
		y.add(beta, p); y.swap(p);
		eslog::accumulatedln("pcpg: update p");

		updateInfo(info, feti.configuration, yw, 0, 0);
		yw = _yw; // keep yw for the next iteration
		eslog::accumulatedln("pcpg: check criteria");
	}
	eslog::endln("pcpg: finished");
	eslog::checkpointln("FETI: PCPG ITERATIONS");
	reconstructSolution(x, r);
	eslog::checkpointln("FETI: SOLUTION RECONSTRUCTION");
	eslog::info("       = ----------------------------------------------------------------------------- = \n");
}

template <> void PCPG<std::complex<double> >::solve(const step::Step &step, IterativeSolverInfo &info)
{

}

}


