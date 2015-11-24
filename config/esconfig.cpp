
#include "esconfig.h"

namespace esconfig {

int MPIrank = 0;
int MPIsize = 1;

namespace mesh {
	size_t subdomains = 1;
	size_t fixPoints = 4;
	Input input = GENERATOR;
	Output output = VTK_FULL;
}

namespace assembler {
	Discretization discretization = FEM; //API;
	Assembler assembler = LinearElasticity;
}
namespace solver {

	double 		epsilon 				= 0.0001;	// Solver requested precision
	size_t 		maxIterations			= 50;		//
	size_t   	FETI_METHOD				= 0; 		// 0 - Total FETI; 1 - HFETI;
	size_t   	USE_SCHUR_COMPLEMENT	= 0; 		// 1 - YES
	size_t		KEEP_FACTORS			= 1; 		// 1 - YES; 0 - NO
	size_t   	PRECONDITIONER			= 0;		// 0 - NO preconditioner; 1 - Lumped
	size_t		CG_SOLVER				= 0;		// 0 - Standard CG; 1 - Pipelined CG
	size_t		REGULARIZATION 			= 0;		// 0 - from mesh; 1 - from stifness matrix

}

}


