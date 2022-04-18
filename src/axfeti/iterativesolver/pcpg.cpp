
#include "pcpg.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "axfeti/projector/projector.h"
#include "axfeti/dualoperator/dualoperator.h"

namespace espreso {

// https://digital.library.unt.edu/ark:/67531/metadc739671/m2/1/high_res_d/792775.pdf
// page 12

// L - number of lambdas
// K - number of DOFs
// R - number of kernels

// K: KxK
// B: LxK
// f: K
// c: L

// G: RxL => Rt * Bt
// e: R => Rt * f
// inv(GGt): RxR

// F: LxL => B * K+ * Bt
// d: L => B * K+ * f
// Q: L => Gt * inv(GGt) * G

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

// inv(GGt): dense stripe -> need gathered vector
// Gt: LxR: L shared per NN -> need NN synchronization
// G: RxL: sparse stripe -> no communication needed when a vector is synchronized
// F: LxL:
//  :: Bt: KxL: no communication needed when a vector is synchronized : COPY TO DOMAIN VECTORS
//  :: K+: KxK: simple threading over all domains
//  :: B: LxK: L shared per NN -> need NN synchronization

/*
 * L: three types of lambdas
 * A: dirichlet -> local lambdas
 * B: gluing -> max. 2 MPI processes
 * C: mortars -> max. NN processes
 * D: rigid body modes -> global across each body
 *
 * communication buffer: 1 + 2
 * without mortars: each lambda for 2 processes max (comm. buffer is a simple copy)
 *
 * TODO: with mortars: lambdas across more processes (non-contiguous buffer!!)
 * try to implement MPI_Datatype and keep the buffer contiguous
 *
 */

template <typename T>
static void _info(PCPG<T> *solver)
{

}

template <typename T>
static void _set(PCPG<T> *solver)
{
	solver->l.resize();
	solver->r.resize();
	solver->w.resize();
	solver->p.resize();

	solver->x.resize();
	solver->Fp.resize();
}

template <typename T>
static void _update(PCPG<T> *solver)
{

}

template <> PCPG<double>::PCPG(AX_FETI<double> *feti): IterativeSolver<double>(feti) { _set<double>(this); }
template <> PCPG<std::complex<double> >::PCPG(AX_FETI<std::complex<double> > *feti): IterativeSolver<std::complex<double> >(feti) { _set<std::complex<double> >(this); }

template <> void PCPG<double>::info() { _info<double>(this); }
template <> void PCPG<std::complex<double> >::info() { _info<std::complex<double> >(this); }

template <> void PCPG<double>::solve(IterativeSolverInfo &info)
{
	DualOperator<double> *F = feti->dualOperator;
	Projector<double> *P = feti->projector;

	// l = Gt * inv(GGt) * e
	P->applyGtInvGGt(P->e, l);

	// r = d - F * l
	F->apply(l, r);
	r.scale(-1);
	r.add(1, F->d);

	// w = P * r
	P->apply(r, w);

	// p = w
	w.copyTo(p);

	// x = l
	l.copyTo(x);

	eslog::info("   iteration       r        e    time[s]\n");
	double ww = w.dot();
	for (size_t i = 0; i < feti->configuration.max_iterations; ++i) {
		double start = eslog::time();

		// gamma = (w, w) / (p, F * p)
		F->apply(p, Fp);
		double pFp = p.dot(Fp), gamma = ww / pFp;

		// x = x + gamma * p
		x.add(gamma, p);

		// r = r - gamma * F * p
		r.add(-gamma, Fp);

		// w = P * r
		P->apply(r, w);

		// beta = (w+1, w+1) / (w, w)
		double _ww = w.dot();
		double beta = _ww / ww;
		ww = _ww; // keep ww for the next iteration

		// p = w + beta * p
		w.add(beta, p); w.swap(p); // w is not used anymore

		eslog::info("%6d  %.4e  %.0e  %7.5f\n", i + 1, ww, feti->configuration.precision, eslog::time() - start);
		if (ww < feti->configuration.precision) {
			break;
		}
	}
}

template <> void PCPG<std::complex<double> >::solve(IterativeSolverInfo &info)
{

}

}


