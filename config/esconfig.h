
#ifndef ESCONFIG_H_
#define ESCONFIG_H_

#include <cstdlib>

namespace esconfig {

enum Discretization {
	FEM,
	BEM
};

extern int MPIrank;
extern int MPIsize;
extern Discretization discretization;

namespace solver {
	extern double 	epsilon;					// Solver requested precision
	extern size_t 	maxIterations;				//
	extern size_t 	FETI_METHOD;				// 0 - Total FETI; 1 - HFETI;
	extern size_t   USE_SCHUR_COMPLEMENT; 		// 1 - YES
	extern size_t	KEEP_FACTORS;				// 1 - YES; 0 - NO
	extern size_t   PRECONDITIONER;				// 0 - NO preconditioner; 1 - Lumped
	extern size_t	CG_SOLVER;					// 0 - Standard CG; 1 - Pipelined CG
	extern size_t	REGULARIZATION;				// 0 - from mesh; 1 - from stifness matrix

}

}


#endif /* ESCONFIG_H_ */
