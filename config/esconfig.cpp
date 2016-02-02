
#include "esconfig.h"

namespace esconfig {

int MPIrank = 0;
int MPIsize = 1;

namespace mesh {
	size_t subdomains = 32;
	size_t fixPoints = 8;

	size_t corners = 1;
	bool vertexCorners = true;
	bool edgeCorners = true;
	bool faceCorners = false;

	bool averageEdges = false;
	bool averageFaces = false;

	Input input = GENERATOR;
	Output output = VTK_FULL;

	double materialDifference = 1e4;
}

namespace assembler {
	Discretization discretization = FEM;
	Assembler assembler = LinearElasticity;
}
namespace solver {

	double epsilon               = 1e-4;// Solver requested precision
	size_t maxIterations         = 500;
	size_t FETI_METHOD           = 1;   // 0 - Total FETI; 1 - HFETI;
	size_t USE_SCHUR_COMPLEMENT  = 0;   // 1 - YES
	size_t KEEP_FACTORS          = 1;   // 1 - YES; 0 - NO
	size_t PRECONDITIONER        = 1;   // 0 - NO preconditioner; 1 - Lumped; 2 - weight function;
	size_t CG_SOLVER             = 0;   // 0 - Standard CG; 1 - Pipelined CG
	size_t REGULARIZATION        = 0;   // 0 - from mesh; 1 - from stiffness matrix
	size_t KSOLVER               = 0;	// 0 - Direct DP, 1 - Iterative solver, 2 - Direct SP,  3 - Direct MIXED Prec
	size_t KSOLVER_SP_iter_steps = 0;   // number of reiteration steps for SP direct solver
	size_t F0_SOLVER             = 0;   // 0 - DIRECT DP if KSOLVER is DIRECT DP - the same precission as KSOLVER
										// 0 - DIRECT SP if KSOLVER is DIRECT SP - the same precission as KSOLVER
										// 1 - DIRECT DP if KSOLVER is DIRECT SP - F0 is in higher precision

}

namespace info {
	std::string output = "log";

	bool printMatrices = true;
}

}


