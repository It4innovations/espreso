
#include "pcpg.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "feti/projector/projector.h"
#include "feti/dualoperator/dualoperator.h"

namespace espreso {

// initialization
// l_0: L => Gt * inv(GGt) * e
// r_0: L => d - F * lambda_0
// w_0: L => r_0 - Gt * inv(GGt) * G * r_0 :: (I - Q) * r_0
// y_0: L -> r_0 - Gt * inv(GGt) * G * S * w_0 :: (I - Q) * S * w_0
// p_0: L => w_0

// loop
// gamma_k: 1 => (p_k,w_k) / (p_k, F * p_k)
//   x_k+1: L => x_k + gama_k * p_k
//   r_k+1: L => r_k - gama_k * F * p_k
//   w_k+1: L => r_k+1 - Gt * inv(GGt) * G * r_k+1 :: (I - Q) * r_k+1
//   y_k+1: L => (I - Q) * S * w_k+1
//  beta_k: 1 => (y_k+1,w_k+1) / (p_k,w_k)
//   p_k+1: L => y_k+1 + beta_k * p_k

template <typename T>
static void _info(PCPG<T> *solver)
{
	eslog::info(" = PRECONDITIONED CONJUGATE PROJECTED GRADIENT SETTINGS                                      = \n");
	switch (solver->feti->configuration.stopping_criterion) {
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
	eslog::info(" =   PRECISION                                                                      %.2e = \n", solver->feti->configuration.precision);
	eslog::info(" =   MAX_ITERATIONS                                                                  %7d = \n", solver->feti->configuration.max_iterations);
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

template <typename T>
static void _set(PCPG<T> *solver)
{
	solver->l.resize();
	solver->r.resize();
	solver->w.resize();
	solver->y.resize();
	solver->z.resize();
	solver->p.resize();

	solver->x.resize();
	solver->Fp.resize();
}


template <> PCPG<double>::PCPG(FETI<double> *feti): IterativeSolver<double>(feti) { _set<double>(this); }
template <> PCPG<std::complex<double> >::PCPG(FETI<std::complex<double> > *feti): IterativeSolver<std::complex<double> >(feti) { _set<std::complex<double> >(this); }

template <> void PCPG<double>::info() { _info<double>(this); }
template <> void PCPG<std::complex<double> >::info() { _info<std::complex<double> >(this); }

template <> void PCPG<double>::solve(IterativeSolverInfo &info)
{
//	DualOperator<double> *F = feti->dualOperator;
//	Projector<double> *P = feti->projector;
//
//	P->applyGtInvGGt(P->e, l);             // l = Gt * inv(GGt) * e
//
//	F->apply(l, r);                        // r = d - F * l
//	r.scale(-1);                           //
//	r.add(1, F->d);                        //
//
//	P->apply(r, w);                        // w = P * r
//	w.copyTo(p);                           // p = w
//	l.copyTo(x);                           // x = l
//
////	r.copyTo(r0);
////	double rho = F->d.dot(l), rr, r0x;
//
//	eslog::checkpointln("FETI: CPG INITIALIZATION");
//	eslog::startln("CPG: ITERATIONS STARTED", "cpg");
//	double pw = w.dot(y);
//	setInfo(info, feti->configuration, pw);
//	while (!info.converged) {
//		F->apply(p, Fp);                            //
//		eslog::accumulatedln("cpg: apply F");       // gamma = (w, w) / (p, F * p)
//		double pFp = p.dot(Fp), gamma = pw / pFp;   //
//		eslog::accumulatedln("cpg: dot(p, Fp)");    //
//
//
//		x.add(gamma, p);                            // x = x + gamma * p
//		r.add(-gamma, Fp);                          // r = r - gamma * F * p
//		eslog::accumulatedln("cpg: update x, r");   //
//
//		P->apply(r, w);                             // w = P * r
//		eslog::accumulatedln("cpg: apply P");       //
//
//		double _yw = w.dot(y), beta = _yw / pw;      // beta = (w+1, w+1) / (w, w)
//		eslog::accumulatedln("cpg: dot(w, w)");     //
//		w.add(beta, p); w.swap(p);                  // p = w + beta * p  (w is not used anymore)
//		eslog::accumulatedln("cpg: update p");      //
//
////		rr = r.dot(), r0x = r0.dot(x);
////		eslog::accumulatedln("cpg: dot(r, r), dot(r0, x)");
//
//		updateInfo(info, feti->configuration, pw, 0, 0);
//		pw = _yw; // keep ww for the next iteration
//		eslog::accumulatedln("cpg: check criteria");
//	}
//	eslog::endln("cpg: finished");
//	eslog::checkpointln("FETI: CPG ITERATIONS");
//	reconstructSolution(x, r);
//	eslog::checkpointln("FETI: SOLUTION RECONSTRUCTION");
//	eslog::info("       = ----------------------------------------------------------------------------- = \n");
}

template <> void PCPG<std::complex<double> >::solve(IterativeSolverInfo &info)
{

}

}


