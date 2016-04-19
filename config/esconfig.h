
#ifndef ESCONFIG_H_
#define ESCONFIG_H_

#include <cstdlib>
#include <string>
#include <vector>

namespace espreso {

namespace input {
class Description;
}

namespace config {



enum FetiMethod {
	TOTAL_FETI,
	HYBRID_FETI
};

enum Preconditioner {
	NO_PRECONDITIONER,
	LUMPED,
	WEIGHT,
	DIRICHLET
};

enum Regularization {
	FIX_POINTS,
	NULL_PIVOTS
};

enum CGSolver {
	STANDARD,
	PIPELINED
};

enum KSolver {
	DIRECT_DOUBLE_PRECISION,
	ITERATIVE,
	DIRECT_SINGLE_PRECISION,
	DIRECT_MIXED_PREXISION
};

enum F0Solver {
	KSOLVER_PRECISION = 0,
	DOUBLE_PRECISION = 1
};

namespace env {
	extern int MPIrank;
	extern int MPIsize;

	extern size_t MKL_NUM_THREADS;
	extern size_t OMP_NUM_THREADS;
	extern size_t SOLVER_NUM_THREADS;
	extern size_t PAR_NUM_THREADS;
	extern size_t CILK_NWORKERS;

	extern std::string executable;
	extern std::string configurationFile;

	extern std::vector<input::Description> description;
};

namespace mesh {
	enum Input { ANSYS_MATSOL, ANSYS_WORKBENCH, OPENFOAM, ESDATA, GENERATOR };
	extern int input;
	extern std::string path;

	extern size_t subdomains;
	extern size_t fixPoints;

	extern size_t corners;
	extern bool vertexCorners;
	extern bool edgeCorners;
	extern bool faceCorners;

	extern bool averageEdges;
	extern bool averageFaces;

	extern double materialDifference;

	extern std::vector<input::Description> description;
};

namespace output {

	enum Output { VTK, ESDATA }; // only VTK_FULL is working
	extern Output format;

	extern bool saveMesh;
	extern bool saveFixPoints;
	extern bool saveFaces;
	extern bool saveLines;
	extern bool saveCorners;
	extern bool saveDirichlet;
	extern bool saveAveraging;
	extern bool saveResults;

	extern double subdomainShrinkRatio;
	extern double clusterShrinkRatio;

	extern std::vector<input::Description> description;
};

namespace assembler {
	enum Discretization { FEM, BEM, API };
	extern int discretization;

	enum Physics { LinearElasticity, Temperature };
	extern int physics;

	extern std::vector<input::Description> description;
};

namespace solver {
	extern double   epsilon;					// Solver requested precision
	extern size_t   maxIterations;				//
	extern size_t   FETI_METHOD;				// 0 - Total FETI; 1 - HFETI;
	extern bool     REDUNDANT_LAGRANGE;
	extern bool     USE_SCHUR_COMPLEMENT; 		// 1 - YES
	extern bool     KEEP_FACTORS;				// 1 - YES; 0 - NO
	extern size_t   PRECONDITIONER;				// 0 - NO preconditioner; 1 - Lumped
	extern size_t   CG_SOLVER;					// 0 - Standard CG; 1 - Pipelined CG
	extern size_t   REGULARIZATION;				// 0 - from mesh; 1 - from stifness matrix
	extern size_t   KSOLVER;					// 0 - Direct DP, 1 - Iter, 2 - Direct SP, 3 - Direct MIXED Prec
	extern size_t   KSOLVER_SP_iter_steps;		// number of reiteration steps for SP direct solver
	extern double   KSOLVER_SP_iter_norm;
	extern size_t   F0_SOLVER;					// 0 - Direct DP if KSOLVER is DIRECT DP
												// 1 - DIRECT SP if KSOLVER is DIRECT SP
												// 1 - Direct DP if KSOLVER is DIRECT SP
	extern size_t   N_MICS;

	extern std::vector<input::Description> description;
};

namespace info {
	extern std::string output;

	extern size_t verboseLevel;
	extern size_t testingLevel;
	extern size_t measureLevel;

	extern bool printMatrices;

	extern std::vector<input::Description> description;
};

}

}


#endif /* ESCONFIG_H_ */
