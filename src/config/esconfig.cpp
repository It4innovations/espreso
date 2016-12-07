
#include "esconfig.h"
#include "esbasis.h"

namespace espreso {

namespace config {


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

double                                   solver::NORM                  = 0;
double                                   solver::EPSILON               = 1e-5;
size_t                                   solver::ITERATIONS            = 1000;
solver::FETI_METHODalternative           solver::FETI_METHOD           = solver::FETI_METHODalternative::TOTAL_FETI;
solver::PRECONDITIONERalternative        solver::PRECONDITIONER        = solver::PRECONDITIONERalternative::LUMPED;
solver::REGULARIZATIONalternative        solver::REGULARIZATION        = solver::REGULARIZATIONalternative::FIX_POINTS;

bool                                     solver::REDUNDANT_LAGRANGE    = true;
bool                                     solver::SCALING               = false;
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
bool                                     solver::LOAD_BALANCING        = true;
size_t                                   solver::TIME_STEPS            = 1;

/////////////////////////////// ASSEMBLER //////////////////////////////////////

assembler::DISCRETIZATIONalternative assembler::DISCRETIZATION = assembler::DISCRETIZATIONalternative::FEM;
assembler::DOFS_ORDERalternative     assembler::DOFS_ORDER     = assembler::DOFS_ORDERalternative::GROUP_ELEMENTS;

//////////////////////////////// HYPRE /////////////////////////////////////////

hypre::SOLVERalternative hypre::HYPRE_SOLVER = hypre::SOLVERalternative::GMRES;
hypre::PRECONDITIONERalternative hypre::HYPRE_PRECONDITIONER = hypre::PRECONDITIONERalternative::DIAGONAL;

}
}


