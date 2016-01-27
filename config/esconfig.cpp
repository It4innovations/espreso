
#include "esconfig.h"

namespace esconfig {

int MPIrank = 0;
int MPIsize = 1;

namespace mesh {
	size_t subdomains = 8;
	size_t fixPoints = 8;

	size_t corners = 3;
	bool vertexCorners = false;
	bool edgeCorners = false;
	bool faceCorners = false;

	bool averageEdges = true;
	bool averageFaces = false;

	Input input = GENERATOR;
	Output output = VTK_FULL;
}

namespace assembler {
	Discretization discretization = FEM;
	Assembler assembler = LinearElasticity;
}
namespace solver {

	double 		epsilon 				= 0.0001;	// Solver requested precision
	size_t 		maxIterations			= 200;		//
	size_t   	FETI_METHOD				= 1; 		// 0 - Total FETI; 1 - HFETI;
	size_t   	USE_SCHUR_COMPLEMENT	= 0; 		// 1 - YES
	size_t		KEEP_FACTORS			= 1; 		// 1 - YES; 0 - NO
	size_t   	PRECONDITIONER			= 0;		// 0 - NO preconditioner; 1 - Lumped; 2 - weight function;
	size_t		CG_SOLVER				= 0;		// 0 - Standard CG; 1 - Pipelined CG
	size_t		REGULARIZATION 			= 1;		// 0 - from mesh; 1 - from stifness matrix
	size_t		KSOLVER					= 0;		// 0 - Direct, 1 - Iter


}

}


