
#include "esconfig.h"
#include "esbasis.h"

namespace espreso {
using namespace input;

namespace config {

int MPIrank = 0;
int MPIsize = 1;
std::string executable;

namespace mesh {

	static std::vector<Description> description = {
			{INTEGER_PARAMETER, "SUBDOMAINS", "Number of subdomains in a cluster."},
			{INTEGER_PARAMETER, "FIXPOINTS" , "Number of fix points in a subdomain."},

			{INTEGER_PARAMETER, "CORNERS"        , "Number of corners on an edge or a face."},
			{BOOLEAN_PARAMETER, "VERTEX_CORNERS" , "Set corners to vertices."},
			{BOOLEAN_PARAMETER, "EDGE_CORNERS"   , "Set corners on edges. The number is defined by parameter CORNERS."},
			{BOOLEAN_PARAMETER, "FACE_CORNERS"   , "Set corners on faces. The number is defined by parameter CORNERS."},
			{BOOLEAN_PARAMETER, "AVERAGE_EDGES"  , "Average nodes on edges."},
			{BOOLEAN_PARAMETER, "AVERAGE_FACES"  , "Average nodes on faces."}
	};

	static Configuration configuration(description, "espreso.config");

	size_t subdomains  = configuration.value("SUBDOMAINS", 8);
	size_t fixPoints   = configuration.value("FIXPOINTS" , 8);

	size_t corners     = configuration.value("CORNERS"       , 1);
	bool vertexCorners = configuration.value("VERTEX_CORNERS", true);
	bool edgeCorners   = configuration.value("EDGE_CORNERS"  , true);
	bool faceCorners   = configuration.value("FACE_CORNERS"  , false);

	bool averageEdges  = configuration.value("AVERAGE_EDGES" , false);
	bool averageFaces  = configuration.value("AVERAGE_FACES" , false);

	Input input = GENERATOR; // set by command line options of in espreso.config file by parameter INPUT

	double materialDifference = 1e-0;
}

namespace output {

	static std::vector<Description> description = {
			{BOOLEAN_PARAMETER, "SAVE_MESH"      , "Save an input mesh."},
			{BOOLEAN_PARAMETER, "SAVE_FIXPOINTS" , "Save a mesh fix points."},
			{BOOLEAN_PARAMETER, "SAVE_FACES"     , "Save faces between subdomains."},
			{BOOLEAN_PARAMETER, "SAVE_EDGES"     , "Save edges among subdomains."},
			{BOOLEAN_PARAMETER, "SAVE_CORNERS"   , "Save corner nodes."},
			{BOOLEAN_PARAMETER, "SAVE_DIRICHLET" , "Save nodes with a dirichlet condition."},
			{BOOLEAN_PARAMETER, "SAVE_AVERAGING" , "Save averaged nodes."},
			{BOOLEAN_PARAMETER, "SAVE_RESULTS"   , "Save the results."},

			{DOUBLE_PARAMETER, "SUBDOMAIN_SHRINK_RATIO", "Shrink ratio for subdomains."},
			{DOUBLE_PARAMETER, "CLUSTER_SHRINK_RATIO"  , "Shrink ratio for clusters."}
	};

	static Configuration configuration(description, "espreso.config");

	Output output = VTK;

	bool saveMesh      = configuration.value("SAVE_MESH"     , false);
	bool saveFixPoints = configuration.value("SAVE_FIXPOINTS", false);
	bool saveFaces     = configuration.value("SAVE_FACES"    , false);
	bool saveLines     = configuration.value("SAVE_EDGES"    , false);
	bool saveCorners   = configuration.value("SAVE_CORNERS"  , false);
	bool saveDirichlet = configuration.value("SAVE_DIRICHLET", false);
	bool saveAveraging = configuration.value("SAVE_AVERAGING", false);
	bool saveResults   = configuration.value("SAVE_RESULTS"  , true);

	double subdomainShrinkRatio = configuration.value("SUBDOMAIN_SHRINK_RATIO", .95);
	double clusterShrinkRatio   = configuration.value("CLUSTER_SHRINK_RATIO"  , .9);
}



namespace assembler {
	Discretization discretization = FEM;
	Assembler assembler = LinearElasticity;
}
namespace solver {

	static std::vector<Description> description = {
			{DOUBLE_PARAMETER,  "EPSILON"              , "Solver requested precision."},
			{INTEGER_PARAMETER, "ITERATIONS"           , "Solver maximum interations."},
			{INTEGER_PARAMETER, "FETI_METHOD"          , "The method used by ESPRESO."},
			{BOOLEAN_PARAMETER, "REDUNDANT_LAGRANGE"   , "Set Lagrange multipliers also among HFETI corners."},
			{BOOLEAN_PARAMETER, "USE_SCHUR_COMPLEMENT" , "Use schur complement to compute ...?"},
			{BOOLEAN_PARAMETER, "KEEP_FACTORS"         , "Keep factors for whole iteration process."},
			{INTEGER_PARAMETER, "PRECONDITIONER"       , "Preconditioner: 0 - NO preconditioner, 1 - Lumped, 2 - weight function, 3 - Dirichlet"},
			{INTEGER_PARAMETER, "CGSOLVER"             , "Conjugate gradients solver: 0 - standard, 1 - pipelined"},
			{INTEGER_PARAMETER, "REGULARIZATION"       , "Regularization of stiffness matrix by: 0 - fix points, 1 - random detection of null pivots"},
			{INTEGER_PARAMETER, "KSOLVER"              , "K is solved: 0 - directly with DP, 1 - iteratively, 2 - directly with SP, 3 - directly with mixed precision"},
			{INTEGER_PARAMETER, "KSOLVER_SP_iter_steps", "Number of reiteration steps for SP direct solver."},
			{INTEGER_PARAMETER, "KSOLVER_SP_iter_norm" , "Number of reiteration steps for SP direct solver."},
			{INTEGER_PARAMETER, "F0SOLVER"             , "F0 is solved: 0 - with the same precision as KSOLVER, 1 - always with DP."},
			{INTEGER_PARAMETER, "N_MICS"               , "Number of MIC accelerators."}
	};

	static Configuration configuration(description, "espreso.config");

	double epsilon               = configuration.value("EPSILON"              , 1e-5);
	size_t maxIterations         = configuration.value("ITERATIONS"           , 1000);
	size_t FETI_METHOD           = configuration.value("FETI_METHOD"          , config::TOTAL_FETI);
	size_t PRECONDITIONER        = configuration.value("PRECONDITIONER"       , config::LUMPED);
	size_t REGULARIZATION        = configuration.value("REGULARIZATION"       , config::FIX_POINTS);
	size_t CG_SOLVER             = configuration.value("CGSOLVER"             , config::STANDARD);
	size_t KSOLVER               = configuration.value("KSOLVER"              , config::DIRECT_DOUBLE_PRECISION);

	bool   REDUNDANT_LAGRANGE    = configuration.value("REDUNDANT_LAGRANGE"   , true);
	bool   USE_SCHUR_COMPLEMENT  = configuration.value("USE_SCHUR_COMPLEMENT" , false);
	bool   KEEP_FACTORS          = configuration.value("KEEP_FACTORS"         , true);

	size_t KSOLVER_SP_iter_steps = configuration.value("KSOLVER_SP_iter_steps", 10);
	double KSOLVER_SP_iter_norm  = configuration.value("KSOLVER_SP_iter_norm" , 1e-12);;
	size_t F0_SOLVER             = configuration.value("F0SOLVER"             , config::KSOLVER_PRECISION);
	size_t N_MICS                = configuration.value("N_MICS"               , 2);

}

namespace info {

	static std::vector<Description> description = {
			{STRING_PARAMETER, "OUTPUT", "An output destination for ESPRESO logs."},

			{INTEGER_PARAMETER, "VERBOSE_LEVEL", "ESPRESO verbose level: <0, 3>"},
			{INTEGER_PARAMETER, "TESTING_LEVEL", "ESPRESO testing level: <0, 3>"},
			{INTEGER_PARAMETER, "MEASURE_LEVEL", "ESPRESO measure level: <0, 3>"},

			{BOOLEAN_PARAMETER, "PRINT_MATRICES", "ESPRESO print all preprocessed matrices."}
	};

	static Configuration configuration(description, "espreso.config");

	std::string output = configuration.value("OUTPUT", "log");

	size_t verboseLevel = configuration.value("VERBOSE_LEVEL", 0);
	size_t testingLevel = configuration.value("TESTING_LEVEL", 0);
	size_t measureLevel = configuration.value("MEASURE_LEVEL", 0);

	bool printMatrices = configuration.value("PRINT_MATRICES", false);
}

}
}


