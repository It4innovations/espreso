
#ifndef ESCONFIG_H_
#define ESCONFIG_H_

#include <cstdlib>
#include <string>

namespace esconfig {

extern int MPIrank;
extern int MPIsize;

namespace mesh {
	extern size_t subdomains;
	extern size_t fixPoints;

	extern size_t corners;
	extern bool vertexCorners;
	extern bool edgeCorners;
	extern bool faceCorners;

	extern bool averageEdges;
	extern bool averageFaces;

	enum Input { ANSYS_MATSOL, ANSYS_WORKBENCH, OPENFOAM, ESDATA_IN, GENERATOR };
	extern Input input;

	enum Output { VTK_FULL, VTK_SURFACE, ESDATA_OUT }; // only VTK_FULL is working
	extern Output output;
}

namespace assembler {
	enum Discretization { FEM, BEM, API };
	extern Discretization discretization;

	enum Assembler { LinearElasticity, Temperature };
	extern Assembler assembler;
}

namespace solver {
	extern double	epsilon;					// Solver requested precision
	extern size_t	maxIterations;				//
	extern size_t	FETI_METHOD;				// 0 - Total FETI; 1 - HFETI;
	extern size_t	USE_SCHUR_COMPLEMENT; 		// 1 - YES
	extern size_t	KEEP_FACTORS;				// 1 - YES; 0 - NO
	extern size_t	PRECONDITIONER;				// 0 - NO preconditioner; 1 - Lumped
	extern size_t	CG_SOLVER;					// 0 - Standard CG; 1 - Pipelined CG
	extern size_t	REGULARIZATION;				// 0 - from mesh; 1 - from stifness matrix
	extern size_t	KSOLVER;					// 0 - Direct DP, 1 - Iter, 2 - Direct SP, 3 - Direct MIXED Prec
	extern size_t   KSOLVER_SP_iter_steps;		// number of reiteration steps for SP direct solver
	extern size_t   F0_SOLVER;					// 0 - Direct DP if KSOLVER is DIRECT DP
												// 1 - DIRECT SP if KSOLVER is DIRECT SP
												// 1 - Direct DP if KSOLVER is DIRECT SP


}

namespace info {
	extern std::string output;

	extern bool printMatrices;
}

namespace tmp{
	extern size_t DOFS;
}



}


#endif /* ESCONFIG_H_ */
