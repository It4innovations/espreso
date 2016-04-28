
#include "esconfig.h"
#include "esbasis.h"

namespace espreso {
using namespace input;

namespace config {

//////////////////////////// ENVIRONMENT ///////////////////////////////////////

int env::MPIrank = 0;
int env::MPIsize = 1;

size_t env::MKL_NUM_THREADS    = Esutils::getEnv<size_t>("MKL_NUM_THREADS");
size_t env::OMP_NUM_THREADS    = Esutils::getEnv<size_t>("OMP_NUM_THREADS");
size_t env::SOLVER_NUM_THREADS = Esutils::getEnv<size_t>("SOLVER_NUM_THREADS");
size_t env::PAR_NUM_THREADS    = Esutils::getEnv<size_t>("PAR_NUM_THREADS");
size_t env::CILK_NWORKERS      = Esutils::getEnv<size_t>("CILK_NWORKERS");

std::string env::executable;
std::string env::configurationFile = "espreso.config";


///////////////////////////////// MESH /////////////////////////////////////////

std::string mesh::path;
int mesh::input = GENERATOR;

size_t mesh::subdomains = 8;
size_t mesh::fixPoints  = 8;

size_t mesh::corners       = 1;
bool   mesh::vertexCorners = true;
bool   mesh::edgeCorners   = true;
bool   mesh::faceCorners   = false;

bool   mesh::averageEdges  = false;
bool   mesh::averageFaces  = false;

/////////////////////////////// SOLVER /////////////////////////////////////////

double   solver::epsilon                = 1e-5;
size_t   solver::maxIterations          = 1000;
size_t   solver::FETI_METHOD            = TOTAL_FETI;
size_t   solver::PRECONDITIONER         = LUMPED;
size_t   solver::CG_SOLVER              = STANDARD;
size_t   solver::REGULARIZATION         = FIX_POINTS;

bool     solver::REDUNDANT_LAGRANGE     = true;
bool     solver::USE_SCHUR_COMPLEMENT   = false;
size_t   solver::SCHUR_COMPLEMENT_PREC  = SC_DOUBLE_PRECISION;
size_t   solver::SCHUR_COMPLEMENT_TYPE  = GENERAL;
bool     solver::COMBINE_SC_AND_SPDS    = true;
bool     solver::KEEP_FACTORS           = true;

size_t   solver::KSOLVER                = DIRECT_DOUBLE_PRECISION;
size_t   solver::KSOLVER_SP_iter_steps  = 1000;
double   solver::KSOLVER_SP_iter_norm   = 1e-12;
size_t   solver::F0_SOLVER              = KSOLVER_PRECISION;

size_t   solver::N_MICS                 = 2;

/////////////////////////////// ASSEMBLER //////////////////////////////////////

int    assembler::discretization = FEM;
int    assembler::physics        = LinearElasticity;
size_t assembler::timeSteps      = 1;

/////////////////////////////// OUTPUT /////////////////////////////////////////

bool output::saveMesh      = false;
bool output::saveFixPoints = false;
bool output::saveFaces     = false;
bool output::saveLines     = false;
bool output::saveCorners   = false;
bool output::saveDirichlet = false;
bool output::saveAveraging = false;
bool output::saveResults   = true;

double output::subdomainShrinkRatio = .95;
double output::clusterShrinkRatio   = .9;

//////////////////////////////// INFO //////////////////////////////////////////

std::string info::output = "log";

size_t info::verboseLevel = 0;
size_t info::testingLevel = 0;
size_t info::measureLevel = 0;

bool info::printMatrices = false;


/////////////////////////////// DESCRIPTION ////////////////////////////////////

std::vector<input::Description> env::description;

std::vector<Description> mesh::description = {
	{ "PATH", mesh::path, "A path to an example.", WRITE_TO_HELP},
	{ "INPUT", mesh::input, "A format of an input.", {
			"matsol",
			"workbench",
			"openfoam",
			"esdata",
			"generator" },  WRITE_TO_HELP},

	{ "SUBDOMAINS", mesh::subdomains, "Number of subdomains in a cluster.", WRITE_TO_HELP },
	{ "FIXPOINTS" , mesh::fixPoints , "Number of fix points in a subdomain." },

	{ "CORNERS"        , mesh::corners      , "Number of corners on an edge or a face." },
	{ "VERTEX_CORNERS" , mesh::vertexCorners, "Set corners to vertices." },
	{ "EDGE_CORNERS"   , mesh::edgeCorners  , "Set corners on edges. The number is defined by parameter CORNERS." },
	{ "FACE_CORNERS"   , mesh::faceCorners  , "Set corners on faces. The number is defined by parameter CORNERS." },

	{ "AVERAGE_EDGES"  , mesh::averageEdges, "Average nodes on edges." },
	{ "AVERAGE_FACES"  , mesh::averageFaces, "Average nodes on faces." }
};

std::vector<Description> assembler::description = {
	{ "DISCRETIZATION", config::assembler::discretization, "A used discretization.",
			{ "FEM", "BEM" }, WRITE_TO_HELP }
};

std::vector<Description> solver::description = {
	{ "EPSILON", solver::epsilon, "Solver requested precision.", WRITE_TO_HELP },
	{ "ITERATIONS", solver::maxIterations, "Solver maximum iterations.", WRITE_TO_HELP },
	{ "FETI_METHOD", solver::FETI_METHOD, "The FETI method used by ESPRESO.", {
			"Total FETI",
			"Hybrid Total FETI" }, WRITE_TO_HELP },

	{ "PRECONDITIONER", solver::PRECONDITIONER, "Preconditioner.", {
			"NO preconditioner",
			"Lumped",
			"weight function",
			"Dirichlet" }, WRITE_TO_HELP },

	{ "CGSOLVER", solver::CG_SOLVER, "Conjugate gradients solver", {
			"standard",
			"pipelined" }, WRITE_TO_HELP },


	{ "REGULARIZATION", solver::REGULARIZATION, "Regularization of stiffness matrix.", {
			"fix points",
			"random detection of null pivots" }},

	{ "KSOLVER", solver::KSOLVER, "K solver precision.", {
			"directly with double precision",
			"iteratively",
			"directly with single precision",
			"directly with mixed precision" }},

	{ "F0SOLVER", solver::F0_SOLVER, "F0 solver precision.", {
			"with the same precision as KSOLVER",
			"always with double precision." }},

	{ "REDUNDANT_LAGRANGE", solver::REDUNDANT_LAGRANGE, "Set Lagrange multipliers also among HFETI corners." },
	{ "USE_SCHUR_COMPLEMENT", solver::USE_SCHUR_COMPLEMENT, "Use schur complement for stiffness matrix processing" },
	{ "SCHUR_COMPLEMENT_PREC", solver::SCHUR_COMPLEMENT_PREC, "Schur complement precision." },
	{ "SCHUR_COMPLEMENT_TYPE", solver::SCHUR_COMPLEMENT_TYPE, "Schur complement matrix type.", {
			"general",
			"symmetric"}},

	{ "COMBINE_SC_AND_SPDS", solver::COMBINE_SC_AND_SPDS, "Combine Schur complement for GPU and sparse direct solver for CPU." },
	{ "KEEP_FACTORS", solver::KEEP_FACTORS, "Keep factors for whole iteration process." },

	{ "KSOLVER_SP_iter_steps", solver::KSOLVER_SP_iter_steps, "Number of reiteration steps for SP direct solver." },
	{ "KSOLVER_SP_iter_norm", solver::KSOLVER_SP_iter_norm , "Number of reiteration steps for SP direct solver." },

	{ "N_MICS", solver::N_MICS, "Number of MIC accelerators.", WRITE_TO_HELP }
};

namespace output {

std::vector<Description> description = {
	{ "SAVE_MESH"      , saveMesh     , "Save an input mesh.", WRITE_TO_HELP },
	{ "SAVE_FIXPOINTS" , saveFixPoints, "Save a mesh fix points." },
	{ "SAVE_FACES"     , saveFaces    , "Save faces between subdomains." },
	{ "SAVE_EDGES"     , saveLines    , "Save edges among subdomains." },
	{ "SAVE_CORNERS"   , saveCorners  , "Save corner nodes." },
	{ "SAVE_DIRICHLET" , saveDirichlet, "Save nodes with a dirichlet condition.", WRITE_TO_HELP },
	{ "SAVE_AVERAGING" , saveAveraging, "Save averaged nodes." },
	{ "SAVE_RESULTS"   , saveResults  , "Save the results.", WRITE_TO_HELP },

	{ "SUBDOMAIN_SHRINK_RATIO", subdomainShrinkRatio, "Shrink ratio for subdomains.", WRITE_TO_HELP },
	{ "CLUSTER_SHRINK_RATIO"  , clusterShrinkRatio  , "Shrink ratio for clusters.", WRITE_TO_HELP }
};

}

namespace info {

std::vector<Description> description = {
	{ "OUTPUT", output, "A location for saving output informations.", WRITE_TO_HELP },
	{ "VERBOSE_LEVEL", verboseLevel, "ESPRESO verbose level.", WRITE_TO_HELP },
	{ "TESTING_LEVEL", verboseLevel, "ESPRESO testing level.", WRITE_TO_HELP },
	{ "MEASURE_LEVEL", verboseLevel, "ESPRESO measure level.", WRITE_TO_HELP },
	{ "PRINT_MATRICES", printMatrices, "ESPRESO print solver input matrices." }
};

}

}
}


