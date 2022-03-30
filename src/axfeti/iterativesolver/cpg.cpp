
#include "cpg.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "axfeti/projector/projector.h"
#include "axfeti/dualoperator/dualoperator.h"

namespace espreso {

// https://digital.library.unt.edu/ark:/67531/metadc739671/m2/1/high_res_d/792775.pdf
// page 12

// initialization
// l_0: L => Gt * inv(GGt) * e
// r_0: L => d - F * lambda_0
// w_0: L => r_0 - Gt * inv(GGt) * G * r_0 :: (I - Q) * r_0
// p_0: L => w_0

// loop
// gamma_k: 1 => (w_k,w_k) / (p_k, F * p_k)
//   x_k+1: L => x_k + gama_k * p_k
//   r_k+1: L => r_k - gama_k * F * p_k
//   w_k+1: L => r_k+1 - Gt * inv(GGt) * G * r_k+1 :: (I - Q) * r_k+1
//  beta_k: 1 => (w_k+1,w_k+1) / (w_k,w_k)
//   p_k+1: L => w_k+1 + beta_k * p_k

template <typename T>
static void _info(CPG<T> *solver)
{
	eslog::info(" = CONJUGATE PROJECTED GRADIENT SETTINGS                                                     = \n");
	switch (solver->feti->configuration.stopping_criterion) {
	case FETIConfiguration::STOPPING_CRITERION::RELATIVE:
		eslog::info(" =   STOPPING CRITERION                                                             RELATIVE = \n");
		break;
	case FETIConfiguration::STOPPING_CRITERION::ABSOLUTE:
		eslog::info(" =   STOPPING CRITERION                                                             ABSOLUTE = \n");
		break;
	}
	eslog::info(" =   PRECISION                                                                      %.2e = \n", solver->feti->configuration.precision);
	eslog::info(" =   MAX_ITERATIONS                                                                  %7d = \n", solver->feti->configuration.max_iterations);
	eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

template <typename T>
static void _set(CPG<T> *solver)
{
	esint lambdas = solver->feti->sinfo.lambdasLocal;

	solver->l.resize(lambdas);
	solver->r.resize(lambdas);
	solver->w.resize(lambdas);
	solver->p.resize(lambdas);

	solver->x.resize(lambdas);
	solver->Fp.resize(lambdas);
}

template <> CPG<double>::CPG(AX_FETI<double> *feti): IterativeSolver<double>(feti) { _set<double>(this); }
template <> CPG<std::complex<double> >::CPG(AX_FETI<std::complex<double> > *feti): IterativeSolver<std::complex<double> >(feti) { _set<std::complex<double> >(this); }

template <> void CPG<double>::info() { _info<double>(this); }
template <> void CPG<std::complex<double> >::info() { _info<std::complex<double> >(this); }

template <> void CPG<double>::solve(IterativeSolverInfo &info)
{
	DualOperator<double> *F = feti->dualOperator;
	Projector<double> *P = feti->projector;

	P->applyGtInvGGt(P->e, l);             // l = Gt * inv(GGt) * e
	F->apply(l, r);                        // r = d - F * l
	r.scale(-1);                           //
	r.add(1, F->d);                        //
	P->apply(r, w);                        // w = P * r
	w.copyTo(p);                           // p = w
	l.copyTo(x);                           // x = l

	// r0 = r;
	// rho = dt * l

	// iter:  (rho + r0 * x)

	eslog::info("   iteration       r        e    time[s]\n");
	double ww = w.dot();
	for (info.iterations = 0; info.iterations < feti->configuration.max_iterations; ++info.iterations) {
		double start = eslog::time();

		F->apply(p, Fp); double pFp = p.dot(Fp), gamma = ww / pFp; // gamma = (w, w) / (p, F * p)
		x.add(gamma, p);                                           // x = x + gamma * p
		r.add(-gamma, Fp);                                         // r = r - gamma * F * p
		P->apply(r, w);                                            // w = P * r
		double _ww = w.dot(), beta = _ww / ww;                     // beta = (w+1, w+1) / (w, w)
		w.add(beta, p); w.swap(p);                                 // p = w + beta * p  (w is not used anymore)

//		eslog::info("%6d  %.4e  %.0e  %7.5f\n", i + 1, std::sqrt(ww), feti->configuration.precision, eslog::time() - start);
		if (std::sqrt(ww) < feti->configuration.precision) {
			break;
		}
		ww = _ww; // keep ww for the next iteration
	}
	reconstructSolution(x, r);
}

template <> void CPG<std::complex<double> >::solve(IterativeSolverInfo &info)
{

}

}


