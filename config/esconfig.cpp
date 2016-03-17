
#include "esconfig.h"

namespace config {

int MPIrank = 0;
int MPIsize = 1;

namespace mesh {
	size_t subdomains = 8;
	size_t fixPoints = 8;

	size_t corners = 1;
	bool vertexCorners = true;
	bool edgeCorners = true;
	bool faceCorners = false;

	bool averageEdges = false;
	bool averageFaces = false;

	Input input = GENERATOR;


	double materialDifference = 1e-0;
}

namespace output {

	Output output = VTK;

	bool saveMesh = false;
	bool saveFixPoints = false;
	bool saveFaces = true;
	bool saveLines = false;
	bool saveCorners = false;
	bool saveDirichlet = false;
	bool saveAveraging = false;
	bool saveResults = true;

	double subdomainShrinkRatio = .95;
	double clusterShrinkRatio = .9;
}



namespace assembler {
	Discretization discretization = FEM;
	Assembler assembler = LinearElasticity;
}
namespace solver {

	double epsilon               = 1e-12;// Solver requested precision
	size_t maxIterations         = 1000;
	size_t FETI_METHOD           = 1;   // 0 - Total FETI; 1 - HFETI;
	bool   REDUNDANT_LAGRANGE    = 1;
	size_t USE_SCHUR_COMPLEMENT  = 1;   // 1 - YES
	size_t KEEP_FACTORS          = 0;   // 1 - YES; 0 - NO
	size_t PRECONDITIONER        = 1;   // 0 - NO preconditioner; 1 - Lumped; 2 - weight function;
	size_t CG_SOLVER             = 0;   // 0 - Standard CG; 1 - Pipelined CG
	size_t REGULARIZATION        = 0;   // 0 - from mesh; 1 - from stiffness matrix
	size_t KSOLVER               = 0;	// 0 - Direct DP, 1 - Iterative solver, 2 - Direct SP,  3 - Direct MIXED Prec
	size_t KSOLVER_SP_iter_steps = 0;   // number of reiteration steps for SP direct solver
	double KSOLVER_SP_iter_norm  = 1e-12;
	size_t F0_SOLVER             = 0;   // 0 - DIRECT DP if KSOLVER is DIRECT DP - the same precission as KSOLVER
										// 0 - DIRECT SP if KSOLVER is DIRECT SP - the same precission as KSOLVER
										// 1 - DIRECT DP if KSOLVER is DIRECT SP - F0 is in higher precision
  size_t N_MICS                = 2;

}

namespace info {
	std::string output = "log";

	size_t verboseLevel = 0;
	size_t testingLevel = 0;
	size_t measureLevel = 0;

	bool printMatrices = false;
}

}


