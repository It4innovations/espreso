
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

double                                   solver::NORM                  = 0;
double                                   solver::EPSILON               = 1e-5;
size_t                                   solver::ITERATIONS            = 1000;
solver::FETI_METHODalternative           solver::FETI_METHOD           = solver::FETI_METHODalternative::TOTAL_FETI;
solver::PRECONDITIONERalternative        solver::PRECONDITIONER        = solver::PRECONDITIONERalternative::LUMPED;
solver::REGULARIZATIONalternative        solver::REGULARIZATION        = solver::REGULARIZATIONalternative::FIX_POINTS;

bool                                     solver::REDUNDANT_LAGRANGE    = true;
solver::B0_TYPEalternative               solver::B0_TYPE               = solver::B0_TYPEalternative::KERNELS;

bool                                     solver::USE_SCHUR_COMPLEMENT  = false;
solver::SCHUR_COMPLEMENT_PRECalternative solver::SCHUR_COMPLEMENT_PREC = solver::SCHUR_COMPLEMENT_PRECalternative::DOUBLE;
solver::SCHUR_COMPLEMENT_TYPEalternative solver::SCHUR_COMPLEMENT_TYPE = solver::SCHUR_COMPLEMENT_TYPEalternative::GENERAL;


bool                                     solver::COMBINE_SC_AND_SPDS   = true;
bool                                     solver::KEEP_FACTORS          = true;


solver::CGSOLVERalternative             solver::CGSOLVER             = solver::CGSOLVERalternative::STANDARD;

solver::KSOLVERalternative               solver::KSOLVER               = solver::KSOLVERalternative::DIRECT_DP;
double                                   solver::KSOLVER_SP_NORM       = 1e-12;
size_t                                   solver::KSOLVER_SP_STEPS      = 1000;

solver::F0SOLVERalternative              solver::F0SOLVER              = solver::F0SOLVERalternative::K_PRECISION;
solver::SASOLVERalternative              solver::SASOLVER              = solver::SASOLVERalternative::CPU_DENSE;

size_t                                   solver::N_MICS                = 2;
size_t                                   solver::TIME_STEPS            = 1;

/////////////////////////////// ASSEMBLER //////////////////////////////////////

assembler::DISCRETIZATIONalternative assembler::DISCRETIZATION = assembler::DISCRETIZATIONalternative::FEM;
assembler::DOFS_ORDERalternative     assembler::DOFS_ORDER     = assembler::DOFS_ORDERalternative::GROUP_ELEMENTS;

/////////////////////////////// OUTPUT /////////////////////////////////////////

bool output::SAVE_MESH       = false;
bool output::SAVE_FIX_POINTS = false;
bool output::SAVE_FACES      = false;
bool output::SAVE_EDGES      = false;
bool output::SAVE_CORNERS    = false;
bool output::SAVE_DIRICHLET  = false;
bool output::SAVE_AVERAGING  = false;
bool output::SAVE_RESULTS    = true;

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
	// ASSEMBLER DESCRIPTION
	{ "DISCRETIZATION", config::assembler::DISCRETIZATION, "A used discretization.", {
			{ "FEM", assembler::DISCRETIZATIONalternative::FEM, "Finite Element Method" },
			{ "BEM", assembler::DISCRETIZATIONalternative::BEM, "Boundary Element Method" } }},

	{ "DOFS_ORDER", config::assembler::DOFS_ORDER, "Order of DOFs in stiffness matrices.", {
			{ "GROUP_ELEMENTS", assembler::DOFS_ORDERalternative::GROUP_ELEMENTS, "Order: x1, y1, x2, y2, ..." },
			{ "GROUP_DOFS", assembler::DOFS_ORDERalternative::GROUP_DOFS, "Order: x1, x2, ..., y1, y2, ..." } }},

	// SOLVER DESCRIPTION
	{ "NORM", solver::NORM, "Solver requested norm.", WRITE_TO_HELP },
	{ "EPSILON", solver::EPSILON, "Solver requested precision.", WRITE_TO_HELP },
	{ "ITERATIONS", solver::ITERATIONS, "Solver maximum iterations.", WRITE_TO_HELP },
	{ "FETI_METHOD", solver::FETI_METHOD, "The FETI method used by ESPRESO.", {
			"Total FETI",
			"Hybrid Total FETI" }, WRITE_TO_HELP },

	{ "PRECONDITIONER", solver::PRECONDITIONER, "Preconditioner.", {
			{ "NONE", solver::PRECONDITIONERalternative::NONE, "Use no preconditioner" },
			{ "LUMPED", solver::PRECONDITIONERalternative::LUMPED, "Lumber preconditioner" },
			{ "WEIGHT_FUNCTION", solver::PRECONDITIONERalternative::WEIGHT_FUNCTION, "Use weight function" },
			{ "DIRICHLET", solver::PRECONDITIONERalternative::DIRICHLET, "Dirichlet preconditioner" },
			{ "SUPER_DIRICHLET", solver::PRECONDITIONERalternative::SUPER_DIRICHLET, "simplified Dirichlet preconditioner" },
			{ "MAGIC", solver::PRECONDITIONERalternative::MAGIC}}, WRITE_TO_HELP },

	{ "CGSOLVER", solver::CGSOLVER, "Conjugate gradients solver", {
			{ "STANDARD", solver::CGSOLVERalternative::STANDARD, "Standard" },
			{ "PIPELINED", solver::CGSOLVERalternative::PIPELINED, "Pipelined" },
			{ "FULL_ORTOGONAL", solver::CGSOLVERalternative::FULL_ORTOGONAL, "Full ortogonalization" },
			{ "GMRES", solver::CGSOLVERalternative::GMRES, "GMRES - allows non-symmetric systems" },
      { "BICGSTAB", solver::CGSOLVERalternative::BICGSTAB, "BICGSTAB - allows non-symmetric systems" }}, WRITE_TO_HELP },


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

	{ "SASOLVER", solver::SA_SOLVER, "SA solver type.", {
					"DENSE solver on CPU",
					"DENSE solver on ACC,"
					"SPARSE solver on CPU." }},


	{ "REDUNDANT_LAGRANGE", solver::REDUNDANT_LAGRANGE, "Set Lagrange multipliers also among HFETI corners." },
	{ "B0_TYPE", solver::B0_TYPE, "The source for B0 assembler." },
	{ "USE_SCHUR_COMPLEMENT", solver::USE_SCHUR_COMPLEMENT, "Use schur complement for stiffness matrix processing" },
	{ "SCHUR_COMPLEMENT_PREC", solver::SCHUR_COMPLEMENT_PREC, "Schur complement precision." },
	{ "SCHUR_COMPLEMENT_TYPE", solver::SCHUR_COMPLEMENT_TYPE, "Schur complement matrix type.", {
			"general",
			"symmetric"}},

	{ "COMBINE_SC_AND_SPDS", solver::COMBINE_SC_AND_SPDS, "Combine Schur complement for GPU and sparse direct solver for CPU." },
	{ "KEEP_FACTORS", solver::KEEP_FACTORS, "Keep factors for whole iteration process." },

	{ "KSOLVER_SP_iter_steps", solver::KSOLVER_SP_STEPS, "Number of reiteration steps for SP direct solver." },
	{ "KSOLVER_SP_iter_norm", solver::KSOLVER_SP_NORM , "Number of reiteration steps for SP direct solver." },

	{ "N_MICS", solver::N_MICS, "Number of MIC accelerators.", WRITE_TO_HELP },

	// OUTPUT DESCRIPTION
	{ "SAVE_MESH"       , output::SAVE_MESH      , "Save an input mesh.", WRITE_TO_HELP },
	{ "SAVE_FIX_POINTS" , output::SAVE_FIX_POINTS, "Save a mesh fix points." },
	{ "SAVE_FACES"      , output::SAVE_FACES      , "Save faces between subdomains." },
	{ "SAVE_EDGES"      , output::SAVE_EDGES      , "Save edges among subdomains." },
	{ "SAVE_CORNERS"    , output::SAVE_CORNERS    , "Save corner nodes." },
	{ "SAVE_DIRICHLET"  , output::SAVE_DIRICHLET  , "Save nodes with a dirichlet condition.", WRITE_TO_HELP },
	{ "SAVE_AVERAGING"  , output::SAVE_AVERAGING  , "Save averaged nodes." },
	{ "SAVE_RESULTS"    , output::SAVE_RESULTS    , "Save the results.", WRITE_TO_HELP },

	{ "SUBDOMAIN_SHRINK_RATIO", output::SUBDOMAINS_SHRINK_RATIO, "Shrink ratio for subdomains.", WRITE_TO_HELP },
	{ "CLUSTER_SHRINK_RATIO"  , output::CLUSTERS_SHRINK_RATIO  , "Shrink ratio for clusters.", WRITE_TO_HELP },

	// INFO DESCRIPTION
	{ "OUTPUT", info::OUTPUT, "A location for saving output informations.", WRITE_TO_HELP },
	{ "VERBOSE_LEVEL", info::VERBOSE_LEVEL, "ESPRESO verbose level.", WRITE_TO_HELP },
	{ "TESTING_LEVEL", info::TESTING_LEVEL, "ESPRESO testing level.", WRITE_TO_HELP },
	{ "MEASURE_LEVEL", info::MEASURE_LEVEL, "ESPRESO measure level.", WRITE_TO_HELP },
	{ "PRINT_MATRICES", info::PRINT_MATRICES, "ESPRESO print solver input matrices." }
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


