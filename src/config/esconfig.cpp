
#include "esconfig.h"
#include "esbasis.h"

namespace espreso {

namespace config {

//////////////////////////// ENVIRONMENT ///////////////////////////////////////

int env::MPIrank = 0;
int env::MPIsize = 1;

size_t env::MKL_NUM_THREADS    = Esutils::getEnv<size_t>("MKL_NUM_THREADS");
size_t env::OMP_NUM_THREADS    = Esutils::getEnv<size_t>("OMP_NUM_THREADS");
size_t env::SOLVER_NUM_THREADS = Esutils::getEnv<size_t>("SOLVER_NUM_THREADS");
size_t env::PAR_NUM_THREADS    = Esutils::getEnv<size_t>("PAR_NUM_THREADS");
size_t env::CILK_NWORKERS      = Esutils::getEnv<size_t>("CILK_NWORKERS");

std::string env::executable    = "espreso"; // Esutils::getEnv<std::string>("_") not works with DDT;
std::string env::configurationFile = "espreso.config";


///////////////////////////////// MESH /////////////////////////////////////////

std::string mesh::PATH;
mesh::INPUTalternative mesh::INPUT = mesh::INPUTalternative::GENERATOR;

size_t mesh::SUBDOMAINS = 8;
size_t mesh::FIX_POINTS  = 8;

size_t mesh::CORNERS       = 1;
bool   mesh::VERTEX_CORNERS = true;
bool   mesh::EDGE_CORNERS   = true;
bool   mesh::FACE_CORNERS   = false;

bool   mesh::AVERAGE_EDGES  = false;
bool   mesh::AVERAGE_FACES  = false;

/////////////////////////////// SOLVER /////////////////////////////////////////

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
solver::SASOLVERalternative 	         solver::SASOLVER              = solver::SASOLVERalternative::CPU_DENSE;

size_t                                   solver::N_MICS                = 2;

/////////////////////////////// ASSEMBLER //////////////////////////////////////

assembler::DISCRETIZATIONalternative assembler::DISCRETIZATION = assembler::DISCRETIZATIONalternative::FEM;
assembler::PHYSICSalternative        assembler::PHYSICS        = assembler::PHYSICSalternative::LINEAR_ELASTICITY;
size_t                               assembler::TIME_STEPS      = 1;

/////////////////////////////// OUTPUT /////////////////////////////////////////

bool output::SAVE_MESH       = false;
bool output::SAVE_FIX_POINTS = false;
bool output::SAVE_FACES      = false;
bool output::SAVE_LINES      = false;
bool output::SAVE_CORNERS    = false;
bool output::SAVE_DIRICHLET  = false;
bool output::SAVE_AVERAGING  = false;
bool output::SAVE_RESULTS    = true;

double output::SUBDOMAINS_SHRINK_RATIO = .95;
double output::CLUSTERS_SHRINK_RATIO   = .9;

//////////////////////////////// INFO //////////////////////////////////////////

std::string info::OUTPUT = "log";

size_t info::VERBOSE_LEVEL = 0;
size_t info::TESTING_LEVEL = 0;
size_t info::MEASURE_LEVEL = 0;

bool info::PRINT_MATRICES = false;


/////////////////////////////// DESCRIPTION ////////////////////////////////////


#define WRITE_TO_HELP  1

std::vector<espreso::Parameter> parameters = {

	// MESH DESCRIPTION
	{ "PATH", mesh::PATH, "A path to an example.", WRITE_TO_HELP },
	{ "INPUT", mesh::INPUT, "A format of an input.", {
			{ "MATSOL", mesh::INPUTalternative::MATSOL, "IT4I internal library" },
			{ "WORKBENCH", mesh::INPUTalternative::WORKBENCH, "Ansys Workbench input file" },
			{ "OPENFOAM", mesh::INPUTalternative::OPENFOAM, "OpenFOAM input format" },
			{ "ESDATA", mesh::INPUTalternative::ESDATA, "ESPRESO binary format" },
			{ "GENERATOR", mesh::INPUTalternative::GENERATOR, "ESPRESO internal generator" } },  WRITE_TO_HELP },

	{ "SUBDOMAINS", mesh::SUBDOMAINS, "Number of subdomains in a cluster.", WRITE_TO_HELP },
	{ "FIX_POINTS" , mesh::FIX_POINTS , "Number of fix points in a subdomain." },

	{ "CORNERS"        , mesh::CORNERS      , "Number of corners on an edge or a face." },
	{ "VERTEX_CORNERS" , mesh::VERTEX_CORNERS, "Set corners to vertices." },
	{ "EDGE_CORNERS"   , mesh::EDGE_CORNERS  , "Set corners on edges. The number is defined by parameter CORNERS." },
	{ "FACE_CORNERS"   , mesh::FACE_CORNERS  , "Set corners on faces. The number is defined by parameter CORNERS." },

	{ "AVERAGE_EDGES"  , mesh::AVERAGE_EDGES, "Average nodes on edges." },
	{ "AVERAGE_FACES"  , mesh::AVERAGE_FACES, "Average nodes on faces." },

	// ASSEMBLER DESCRIPTION
	{ "DISCRETIZATION", config::assembler::DISCRETIZATION, "A used discretization.",{
			{ "FEM", assembler::DISCRETIZATIONalternative::FEM, "Finite Element Method" },
			{ "BEM", assembler::DISCRETIZATIONalternative::BEM, "Boundary Element Method" } }},

	// SOLVER DESCRIPTION
	{ "EPSILON", solver::EPSILON, "Solver requested precision.", WRITE_TO_HELP },
	{ "ITERATIONS", solver::ITERATIONS, "Solver maximum iterations.", WRITE_TO_HELP },
	{ "FETI_METHOD", solver::FETI_METHOD, "The FETI method used by ESPRESO.", {
			{ "TOTAL_FETI", solver::FETI_METHODalternative::TOTAL_FETI, "Total FETI." },
			{ "HYBRID_FETI", solver::FETI_METHODalternative::HYBRID_FETI, "Hybrid Total FETI." } }, WRITE_TO_HELP },

	{ "PRECONDITIONER", solver::PRECONDITIONER, "Preconditioner.", {
			{ "NONE", solver::PRECONDITIONERalternative::NONE, "Use no preconditioner" },
			{ "LUMPED", solver::PRECONDITIONERalternative::LUMPED, "Lumber preconditioner" },
			{ "WEIGHT_FUNCTION", solver::PRECONDITIONERalternative::WEIGHT_FUNCTION, "Use weight function" },
			{ "DIRICHLET", solver::PRECONDITIONERalternative::DIRICHLET, "Dirichlet preconditioner" },
			{ "MAGIC", solver::PRECONDITIONERalternative::MAGIC}}, WRITE_TO_HELP },

	{ "CGSOLVER", solver::CGSOLVER, "Conjugate gradients solver", {
			{ "STANDARD", solver::CGSOLVERalternative::STANDARD, "Standard" },
			{ "PIPELINED", solver::CGSOLVERalternative::PIPELINED, "Pipelined" },
			{ "FULL_ORTOGONAL", solver::CGSOLVERalternative::FULL_ORTOGONAL, "Full ortogonalization" },}, WRITE_TO_HELP },


	{ "REGULARIZATION", solver::REGULARIZATION, "Regularization of stiffness matrix.", {
			{ "FIX_POINTS", solver::REGULARIZATIONalternative::FIX_POINTS, "From fix points" },
			{ "NULL_PIVOTS", solver::REGULARIZATIONalternative::NULL_PIVOTS, "Random null pivots" } }},

	{ "KSOLVER", solver::KSOLVER, "K solver type.", {
			{ "DIRECT_DP", solver::KSOLVERalternative::DIRECT_DP, "Directly with double precision" },
			{ "ITERATIVE", solver::KSOLVERalternative::ITERATIVE, "Iteratively" },
			{ "DIRECT_SP", solver::KSOLVERalternative::DIRECT_SP, "Directly with single precision" },
			{ "DIRECT_MX", solver::KSOLVERalternative::DIRECT_MP, "Directly with mixed precision" } }},

	{ "F0SOLVER", solver::F0SOLVER, "F0 solver precision.", {
			{ "K_PRECISION", solver::F0SOLVERalternative::K_PRECISION, "With the same precision as KSOLVER" },
			{ "DOUBLE", solver::F0SOLVERalternative::DOUBLE, "Always with double precision" } }},

	{ "SASOLVER", solver::SASOLVER, "SA solver type.", {
			{ "CPU_DENSE", solver::SASOLVERalternative::CPU_DENSE, "DENSE solver on CPU" },
			{ "ACC_DENSE", solver::SASOLVERalternative::ACC_DENSE, "DENSE solver on ACC" },
			{ "CPU_SPARSE", solver::SASOLVERalternative::CPU_SPARSE, "SPARSE solver on CPU." } }},


	{ "REDUNDANT_LAGRANGE", solver::REDUNDANT_LAGRANGE, "Set Lagrange multipliers also among HFETI corners." },
	{ "B0_TYPE", solver::B0_TYPE, "Type of cluster gluing matrix.", {
			{"CORNERS", solver::B0_TYPEalternative::CORNERS, "Gluing based on corners."},
			{"KERNELS", solver::B0_TYPEalternative::KERNELS, "Gluing based on kernels."} }},

	{ "USE_SCHUR_COMPLEMENT", solver::USE_SCHUR_COMPLEMENT, "Use schur complement for stiffness matrix processing" },
	{ "SCHUR_COMPLEMENT_PREC", solver::SCHUR_COMPLEMENT_PREC, "Schur complement precision.", {
			{"DOUBLE", solver::SCHUR_COMPLEMENT_PRECalternative::DOUBLE, "Double precision"},
			{"SINGLE", solver::SCHUR_COMPLEMENT_PRECalternative::SINGLE, "Single precision"}} },
	{ "SCHUR_COMPLEMENT_TYPE", solver::SCHUR_COMPLEMENT_TYPE, "Schur complement matrix type.", {
			{ "GENERAL", solver::SCHUR_COMPLEMENT_TYPEalternative::GENERAL, "Store a full matrix" },
			{ "SYMMETRIC", solver::SCHUR_COMPLEMENT_TYPEalternative::SYMMETRIC, "Store only a triangle" } }},

	{ "COMBINE_SC_AND_SPDS", solver::COMBINE_SC_AND_SPDS, "Combine Schur complement for GPU and sparse direct solver for CPU." },
	{ "KEEP_FACTORS", solver::KEEP_FACTORS, "Keep factors for whole iteration process." },

	{ "KSOLVER_SP_iter_steps", solver::KSOLVER_SP_STEPS, "Number of reiteration steps for SP direct solver." },
	{ "KSOLVER_SP_iter_norm", solver::KSOLVER_SP_NORM , "Number of reiteration steps for SP direct solver." },

	{ "N_MICS", solver::N_MICS, "Number of MIC accelerators.", WRITE_TO_HELP },

	// OUTPUT DESCRIPTION
	{ "SAVE_MESH"       , output::SAVE_MESH      , "Save an input mesh.", WRITE_TO_HELP },
	{ "SAVE_FIX_POINTS" , output::SAVE_FIX_POINTS, "Save a mesh fix points." },
	{ "SAVE_FACES"      , output::SAVE_FACES      , "Save faces between subdomains." },
	{ "SAVE_EDGES"      , output::SAVE_LINES      , "Save edges among subdomains." },
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

}


