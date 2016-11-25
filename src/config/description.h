
#ifndef SRC_CONFIG_DESCRIPTION_H_
#define SRC_CONFIG_DESCRIPTION_H_

#include "configuration.h"

namespace espreso {

enum class INPUT {
	/// Ansys Workbench format
	WORKBENCH = 0,
	/// OpenFOAM format
	OPENFOAM = 1,
	/// ESPRESO binary format
	ESDATA = 2,
	/// ESPRESO internal problem generator
	GENERATOR = 3
};

enum class FETI_METHOD {
	/// Total FETI
	TOTAL_FETI = 0,
	/// Hybrid Total FETI
	HYBRID_FETI = 1,
};

enum class FETI_ITERATIVE_SOLVER {
	/// Projected conjugate gradients
	PCG = 0,
	/// Pipelined PCG
	pipePCG = 1,
	/// Full orthogonal PCG
	orthogonalPCG = 2,
	/// GMRES
	GMRES = 3,
	/// BICGSTAB
	BICGSTAB = 4,
	/// QPCE
	QPCE = 5
};

enum class FETI_PRECONDITIONER {
	/// No preconditioner is used
	NONE = 0,
	/// Lumped preconditioner     S = K_ss
	LUMPED = 1,
	/// Weight function
	WEIGHT_FUNCTION = 2,
	/// Dirichlet preconditioner  S = K_ss - K_sr * inv(K_rr) * K_sr
	DIRICHLET = 3,
	/// simplified Dirichlet      S = K_ss - K_sr * 1/diag(K_rr) * K_sr
	SUPER_DIRICHLET = 4,
	/// Lubos's preconditioner
	MAGIC = 5
};

enum class REGULARIZATION {
	/// Fix points
	FIX_POINTS = 0,
	/// Randomly found null pivots of stiffness matrix
	NULL_PIVOTS = 1
};

enum class B0_TYPE {
	/// Gluing based on corners
	CORNERS = 0,
	/// Gluing based on kernels of faces
	KERNELS = 1,
	/// Both corners and kernels
	COMBINED = 2
};

enum class FLOAT_PRECISION {
	/// Double precision
	DOUBLE = 0,
	/// Single precision
	SINGLE = 1
};

enum class FETI_KSOLVER {
	/// A direct solver with double precision
	DIRECT_DP = 0,
	/// An iterative solver
	ITERATIVE = 1,
	/// A direct solver with single precision
	DIRECT_SP = 2,
	/// A direct solver with mixed precision
	DIRECT_MP = 3
};

enum class FETI_F0SOLVER_PRECISION {
	/// The same precision as K solver
	K_PRECISION = 0,
	/// Double precision
	DOUBLE = 1
};

enum class FETI_SASOLVER {
	/// Dense on CPU
	CPU_DENSE = 0,
	/// Dense on accelerator
	ACC_DENSE = 1,
	/// Sparse on CPU
	CPU_SPARSE = 2
};

enum class MATRIX_STORAGE {
	/// A full matrix is stored
	GENERAL = 0,
	/// Store only triangle
	SYMMETRIC = 1
};

enum class MULTIGRID_SOLVER {
	CG = 0,
	GMRES = 1,
	FGMRES = 2,
	BICGS = 3,
	BICGSTAB = 4,
	TFQMR = 5,
	SYMQMR = 6,
	SUPERLU = 7,
	SUPERLUX = 8
};

enum class MULTIGRID_PRECONDITIONER {
	DIAGONAL = 0,
	PILUT = 1,
	EUCLID = 2,
	PARASAILS = 3,
	BOOMERAMG = 4,
	POLY = 5,
	MLI = 6
};

enum class OUTPUT_FORMAT {
	VTK_LEGACY_FORMAT = 0,
	VTK_BINARY_FORMAT = 1,
	VTK_MULTIBLOCK_FORMAT = 2,
	ENSIGHT_FORMAT = 3
};

enum class PHYSICS {
	LINEAR_ELASTICITY_2D,
	LINEAR_ELASTICITY_3D,
	TRANSIENT_ELASTICITY_2D,
	TRANSIENT_ELASTICITY_3D,
	ADVECTION_DIFFUSION_2D,
	ADVECTION_DIFFUSION_3D,
	STOKES
};

enum class MATERIAL_MODEL {
	LINEAR_ELASTIC_ISOTROPIC = 0,
	LINEAR_ELASTIC_ORTHOTROPIC = 1,
	LINEAR_ELASTIC_ANISOTROPIC = 2
};

enum class GENERATOR_SHAPE {
	GRID,
	SPHERE
};

enum class ELEMENT_TYPE {
	HEXA8,
	HEXA20,
	TETRA4,
	TETRA10,
	PRISMA6,
	PRISMA15,
	PYRAMID5,
	PYRAMID13,

	SQUARE4,
	SQUARE8,
	TRIANGLE3,
	TRIANGLE6
};

struct Environment: public Configuration {

	PARAMETER(int, MPIrank, "Rank of an MPI process.", 0);
	PARAMETER(int, MPIsize, "Size of an MPI communicator (MPI_COMM_WORLD).", 1);

	PARAMETER(size_t, MKL_NUM_THREADS, "Number of MKL threads.", Esutils::getEnv<size_t>("MKL_NUM_THREADS"));
	PARAMETER(size_t, OMP_NUM_THREADS, "Number of OMP threads.", Esutils::getEnv<size_t>("OMP_NUM_THREADS"));
	PARAMETER(size_t, SOLVER_NUM_THREADS, "Number of threads used in ESPRESO solver.", Esutils::getEnv<size_t>("SOLVER_NUM_THREADS"));
	PARAMETER(size_t, PAR_NUM_THREADS, "Number of parallel threads.", Esutils::getEnv<size_t>("PAR_NUM_THREADS"));
	PARAMETER(size_t, CILK_NWORKERS, "Number of cilk++ threads.", Esutils::getEnv<size_t>("CILK_NWORKERS"));

	PARAMETER(std::string, executable, "name", "espreso");

};

struct GridGenerator: Configuration {

	PARAMETER(double, start_x, "x-coordinate of grid starting point.", 0);
	PARAMETER(double, start_y, "y-coordinate of grid starting point.", 0);
	PARAMETER(double, start_z, "z-coordinate of grid starting point.", 0);
	PARAMETER(double, length_x, "x-length of generated grid.", 1);
	PARAMETER(double, length_y, "y-length of generated grid.", 1);
	PARAMETER(double, length_z, "z-length of generated grid.", 1);

	PARAMETER(double, blocks_x, "Number of blocks in x-direction of a grid.", 1);
	PARAMETER(double, blocks_y, "Number of blocks in y-direction of a grid.", 1);
	PARAMETER(double, blocks_z, "Number of blocks in z-direction of a grid.", 1);

	PARAMETER(double, clusters_x, "Number of clusters in x-direction of each grid square.", 1);
	PARAMETER(double, clusters_y, "Number of clusters in y-direction of each grid square.", 1);
	PARAMETER(double, clusters_z, "Number of clusters in z-direction of each grid square.", 1);

	PARAMETER(double, domains_x, "Number of domains in x-direction of each cluster.", 2);
	PARAMETER(double, domains_y, "Number of domains in y-direction of each cluster.", 2);
	PARAMETER(double, domains_z, "Number of domains in z-direction of each cluster.", 2);

	PARAMETER(double, elements_x, "Number of elements in x-direction of each domain.", 5);
	PARAMETER(double, elements_y, "Number of elements in y-direction of each domain.", 5);
	PARAMETER(double, elements_z, "Number of elements in z-direction of each domain.", 5);

	SUBMAP(size_t, bool, blocks, "List of grid blocks [<INDEX> <VALUE>]. Where value indicate if a block will be generated.", "<INDEX>", true);

	SUBMAP(std::string, std::string, nodes, "List of nodes regions.", "<REGION_NAME>", "<INTERVAL / PATTERN[ALL]>");
	SUBMAP(std::string, std::string, edges, "List of edges regions.", "<REGION_NAME>", "<INTERVAL / PATTERN[ALL]>");
	SUBMAP(std::string, std::string, faces, "List of faces regions.", "<REGION_NAME>", "<INTERVAL / PATTERN[ALL]>");
	SUBMAP(std::string, std::string, elements, "List of elements regions.", "<REGION_NAME>", "<INTERVAL / PATTERN[ALL]>");
};

struct ESPRESOGenerator: public Configuration {

	OPTION(GENERATOR_SHAPE, shape, "Generated shape", GENERATOR_SHAPE::GRID, OPTIONS({
		{ "GRID"  , GENERATOR_SHAPE::GRID  , "Rectangular grid with empty spaces." },
		{ "SPHERE", GENERATOR_SHAPE::SPHERE, "Hollow sphere." }
	}));

	OPTION(ELEMENT_TYPE, element, "Type of generated element", ELEMENT_TYPE::HEXA8, OPTIONS({
		{ "HEXA8"    , ELEMENT_TYPE::HEXA8    , "Hexahedron."},
		{ "HEXA20"   , ELEMENT_TYPE::HEXA20   , "Hexahedron with midpoints."},
		{ "TETRA4"   , ELEMENT_TYPE::TETRA4   , "Tetrahedron."},
		{ "TETRA10"  , ELEMENT_TYPE::TETRA10  , "Tetrahedron with midpoints."},
		{ "PRISMA6"  , ELEMENT_TYPE::PRISMA6  , "Prisma."},
		{ "PRISMA15" , ELEMENT_TYPE::PRISMA15 , "Prisma with midpoints."},
		{ "PYRAMID5" , ELEMENT_TYPE::PYRAMID5 , "Pyramid."},
		{ "PYRAMID13", ELEMENT_TYPE::PYRAMID13, "Pyramid with midpoints."},

		{ "SQUARE4"  , ELEMENT_TYPE::SQUARE4  , "Square."},
		{ "SQUARE8"  , ELEMENT_TYPE::SQUARE8  , "Square with midpoints."},
		{ "TRIANGLE3", ELEMENT_TYPE::TRIANGLE3, "Triangle."},
		{ "TRIANGLE6", ELEMENT_TYPE::TRIANGLE6, "Triangle with midpoints."},
	}));

	SUBCONFIG(GridGenerator, grid, "Detailed specification of grid shape.");
};

struct FETISolver: public Configuration {

	PARAMETER(double, epsilon, "Solver requested precision.", 1e-5);
	PARAMETER(size_t, iterations, "solver iterations", 100);

	OPTION(FETI_METHOD, method, "The FETI method used by ESPRESO.", FETI_METHOD::TOTAL_FETI, OPTIONS({
		{ "TOTAL_FETI" , FETI_METHOD::TOTAL_FETI , "Total FETI." },
		{ "HYBRID_FETI", FETI_METHOD::HYBRID_FETI, "Hybrid Total FETI." }
	}));

	OPTION(FETI_ITERATIVE_SOLVER, solver, "Used iterative solver", FETI_ITERATIVE_SOLVER::PCG, OPTIONS({
		{ "PCG"          , FETI_ITERATIVE_SOLVER::PCG          , "Standard Projected conjugate gradients." },
		{ "PIPEPCG"      , FETI_ITERATIVE_SOLVER::pipePCG      , "Pipelined PCG." },
		{ "ORTHOGONALPCG", FETI_ITERATIVE_SOLVER::orthogonalPCG, "Full ortogonalization PCG." },
		{ "GMRES"        , FETI_ITERATIVE_SOLVER::GMRES        , "GMRES - allows non-symmetric systems." },
		{ "BICGSTAB"     , FETI_ITERATIVE_SOLVER::BICGSTAB     , "BICGSTAB - allows non-symmetric systems." },
		{ "QPCE"         , FETI_ITERATIVE_SOLVER::QPCE         , "QPCE - allows contact." }
	}));

	OPTION(FETI_PRECONDITIONER, preconditioner, "Preconditioner", FETI_PRECONDITIONER::LUMPED, OPTIONS({
		{ "NONE"           , FETI_PRECONDITIONER::NONE           , "Use no preconditioner." },
		{ "LUMPED"         , FETI_PRECONDITIONER::LUMPED         , "Lumber preconditioner." },
		{ "WEIGHT_FUNCTION", FETI_PRECONDITIONER::WEIGHT_FUNCTION, "Use weight function." },
		{ "DIRICHLET"      , FETI_PRECONDITIONER::DIRICHLET      , "Dirichlet preconditioner." },
		{ "SUPER_DIRICHLET", FETI_PRECONDITIONER::SUPER_DIRICHLET, "simplified Dirichlet preconditioner." },
		{ "MAGIC"          , FETI_PRECONDITIONER::MAGIC          , "TODO." }
	}));

	OPTION(REGULARIZATION, regularization, "Type of regularization of stiffness matrix", REGULARIZATION::FIX_POINTS, OPTIONS({
		{ "FIX_POINTS" , REGULARIZATION::FIX_POINTS , "From fix points." },
		{ "NULL_PIVOTS", REGULARIZATION::NULL_PIVOTS, "Random null pivots." }
	}));

	PARAMETER(bool, redundant_lagrange, "If true, each pair of DOF are glued", true);
	PARAMETER(bool, scaling, "If true, Lagrange multiplicators are weighted according a value in the stiffness matrix", true);

	OPTION(B0_TYPE, B0_type, "Type of gluing matrix for Hybrid FETI", B0_TYPE::KERNELS, OPTIONS({
		{"CORNERS" , B0_TYPE::CORNERS , "Gluing based on corners."},
		{"KERNELS" , B0_TYPE::KERNELS , "Gluing based on kernels."},
		{"COMBINED", B0_TYPE::COMBINED, "Both corners and kernels."}
	}));

	PARAMETER(bool, use_schur_complement, "Use schur complement for stiffness matrix processing.", false);

	OPTION(FLOAT_PRECISION, schur_precision, "Precision of Schur complement", FLOAT_PRECISION::DOUBLE, OPTIONS({
		{"DOUBLE", FLOAT_PRECISION::DOUBLE, "Double precision."},
		{"SINGLE", FLOAT_PRECISION::SINGLE, "Single precision."}
	}));

	OPTION(FETI_KSOLVER, Ksolver, "Stiffness matrix solver type", FETI_KSOLVER::DIRECT_DP, OPTIONS({
		{ "DIRECT_DP", FETI_KSOLVER::DIRECT_DP, "Directly with double precision" },
		{ "ITERATIVE", FETI_KSOLVER::ITERATIVE, "Iteratively" },
		{ "DIRECT_SP", FETI_KSOLVER::DIRECT_SP, "Directly with single precision" },
		{ "DIRECT_MX", FETI_KSOLVER::DIRECT_MP, "Directly with mixed precision" }
	}));

	PARAMETER(size_t, Ksolver_iterations, "Number of reiteration steps for single precision direct solver", 1000);
	PARAMETER(double, Ksolver_epsilon   , "Reguested norm for single precision direct solver", 1e-12);

	OPTION(FETI_F0SOLVER_PRECISION, F0_precision, "Precision of F0", FETI_F0SOLVER_PRECISION::K_PRECISION, OPTIONS({
		{ "K_PRECISION", FETI_F0SOLVER_PRECISION::K_PRECISION, "With the same precision as KSOLVER" },
		{ "DOUBLE"     , FETI_F0SOLVER_PRECISION::DOUBLE     , "Always with double precision" }
	}));

	OPTION(FETI_SASOLVER, SAsolver, "S alfa solver type", FETI_SASOLVER::CPU_DENSE, OPTIONS({
		{ "CPU_DENSE" , FETI_SASOLVER::CPU_DENSE , "Dense solver on CPU" },
		{ "ACC_DENSE" , FETI_SASOLVER::ACC_DENSE , "Dense solver on ACC" },
		{ "CPU_SPARSE", FETI_SASOLVER::CPU_SPARSE, "Sparse solver on CPU." }
	}));

	OPTION(MATRIX_STORAGE, schur_type, "Schur complement matrix type.", MATRIX_STORAGE::GENERAL, OPTIONS({
		{ "GENERAL"  , MATRIX_STORAGE::GENERAL, "Store a full matrix." },
		{ "SYMMETRIC", MATRIX_STORAGE::SYMMETRIC, "Store only triangle." }
	}));

	PARAMETER(bool, combine_sc_and_spds, "Combine usage of SC for Accelerator and sparse direct solver for CPU.", true);
	PARAMETER(bool, keep_factors, "Keep factors between iterations.", true);

	PARAMETER(size_t, N_MICS, "Number of MIC accelerators", 2);

	PARAMETER(size_t, time_steps, "Number of time steps for transient problems", 1);
};

struct MULTIGRIDSolver: public Configuration {

	OPTION(MULTIGRID_SOLVER, solver, "Used solver", MULTIGRID_SOLVER::CG, OPTIONS({
		{"CG"      , MULTIGRID_SOLVER::CG      , "CG solver." },
		{"GMRES"   , MULTIGRID_SOLVER::GMRES   , "GMRES solver." },
		{"FGMRES"  , MULTIGRID_SOLVER::FGMRES  , "FGMRES solver." },
		{"BICGS"   , MULTIGRID_SOLVER::BICGS   , "BICGS solver." },
		{"BICGSTAB", MULTIGRID_SOLVER::BICGSTAB, "BICGSTAB solver." },
		{"TFQMR"   , MULTIGRID_SOLVER::TFQMR   , "TFQMR solver." },
		{"SYMQMR"  , MULTIGRID_SOLVER::SYMQMR  , "SYMQMR solver." },
		{"SUPERLU" , MULTIGRID_SOLVER::SUPERLU , "SUPERLU solver." },
		{"SUPERLUX", MULTIGRID_SOLVER::SUPERLUX, "SUPERLUX solver." }
	}));

	OPTION(MULTIGRID_PRECONDITIONER, preconditioner, "Used preconditioner", MULTIGRID_PRECONDITIONER::BOOMERAMG, OPTIONS({
		{"DIAGONAL" , MULTIGRID_PRECONDITIONER::DIAGONAL , "DIAGONAL preconditioner." },
		{"PILUT"    , MULTIGRID_PRECONDITIONER::PILUT    , "PILUT preconditioner." },
		{"EUCLID"   , MULTIGRID_PRECONDITIONER::EUCLID   , "EUCLID preconditioner." },
		{"PARASAILS", MULTIGRID_PRECONDITIONER::PARASAILS, "PARASAILS preconditioner." },
		{"BOOMERAMG", MULTIGRID_PRECONDITIONER::BOOMERAMG, "BOOMERAMG preconditioner." },
		{"POLY"     , MULTIGRID_PRECONDITIONER::POLY     , "POLY preconditioner." },
		{"MLI"      , MULTIGRID_PRECONDITIONER::MLI      , "MLI preconditioner." }
	}));
};

struct Output: public Configuration {

	OPTION(OUTPUT_FORMAT, format, "Format - only LEGACY format is supported without VTK library", OUTPUT_FORMAT::VTK_LEGACY_FORMAT, OPTIONS({
		{ "VTK_LEGACY"    , OUTPUT_FORMAT::VTK_LEGACY_FORMAT    , "*.vtk files" },
		{ "VTK_BINARY"    , OUTPUT_FORMAT::VTK_BINARY_FORMAT    , "*.vtu files" },
		{ "VTK_MULTIBLOCK", OUTPUT_FORMAT::VTK_MULTIBLOCK_FORMAT, "*.vtu + *.vtm files" },
		{ "ENSIGHT"       , OUTPUT_FORMAT::ENSIGHT_FORMAT       , "EnSight files" }
	}));

	PARAMETER(bool, compression, "Compression - needs VTK library", false);
	PARAMETER(double, decimation, "Decimation - needs VTK library", 0);

	PARAMETER(bool, result, "Save result", true);
	PARAMETER(bool, properties, "Save also input parameters", false);
	PARAMETER(bool, gluing, "Save lagrange multipliers", false);

	PARAMETER(double, subdomains_shring_ratio, "All subdomains are shrunk by this ratio", .95);
	PARAMETER(double, clusters_shring_ratio  , "All clusters are shrunk by this ratio"  , .9);

	PARAMETER(std::string, log_dir, "Log directory.", "log");

	PARAMETER(size_t, verbose_level, "Verbose level [0-3]", 0);
	PARAMETER(size_t, testing_level, "Testing level [0-3]", 0);
	PARAMETER(size_t, measure_level, "Measure level [0-3]", 0);

	PARAMETER(bool, print_matrices, "Print assembler matrices.", false);
};

struct MaterialParameters: public Configuration {

	PARAMETER(std::string, DENS, "Density"             , "7850");
	PARAMETER(std::string, MIXY, "Poisson ratio XY."   , "0.3");
	PARAMETER(std::string, MIXZ, "Poisson ratio XZ."   , "0.3");
	PARAMETER(std::string, MIYZ, "Poisson ratio YZ."   , "0.3");
	PARAMETER(std::string, EX  , "Young modulus X."    , "2.1e11");
	PARAMETER(std::string, EY  , "Young modulus Y."    , "2.1e11");
	PARAMETER(std::string, EZ  , "Young modulus Z."    , "2.1e11");
	PARAMETER(std::string, C   , "Termal capacity."    , "1");
	PARAMETER(std::string, KXX , "Termal conduction X.", "1");
	PARAMETER(std::string, KXY , "Termal conduction Y.", "1");
	PARAMETER(std::string, KXZ , "Termal conduction Z.", "1");
	PARAMETER(std::string, ALPX, "Termal expansion X." , "1");
	PARAMETER(std::string, ALPY, "Termal expansion Y." , "1");
	PARAMETER(std::string, ALPZ, "Termal expansion Z." , "1");

	OPTION(MATERIAL_MODEL, model, "Material model", MATERIAL_MODEL::LINEAR_ELASTIC_ISOTROPIC, OPTIONS({
		{ "LINEAR_ELASTIC_ISOTROPIC"  , MATERIAL_MODEL::LINEAR_ELASTIC_ISOTROPIC  , "Isotropic." },
		{ "LINEAR_ELASTIC_ORTHOTROPIC", MATERIAL_MODEL::LINEAR_ELASTIC_ORTHOTROPIC, "Orthotropic." },
		{ "LINEAR_ELASTIC_ANISOTROPIC", MATERIAL_MODEL::LINEAR_ELASTIC_ANISOTROPIC, "Anisotropic." }
	}));
};

struct GlobalConfiguration: public Configuration {

	OPTION(INPUT, input, "test input", INPUT::GENERATOR, OPTIONS({
			{ "WORKBENCH", INPUT::WORKBENCH, "Ansys Workbench input file" },
			{ "OPENFOAM", INPUT::OPENFOAM, "OpenFOAM input format" },
			{ "ESDATA", INPUT::ESDATA, "ESPRESO binary format" },
			{ "GENERATOR", INPUT::GENERATOR, "ESPRESO internal generator" }
	}));

	OPTION(PHYSICS, physics, "Used physics", PHYSICS::LINEAR_ELASTICITY_3D, OPTIONS({
		{ "LINEAR_ELASTICITY_2D"   , PHYSICS::LINEAR_ELASTICITY_2D   , "2D linear elasticity." },
		{ "LINEAR_ELASTICITY_3D"   , PHYSICS::LINEAR_ELASTICITY_3D   , "3D linear elasticity." },
		{ "TRANSIENT_ELASTICITY_2D", PHYSICS::TRANSIENT_ELASTICITY_2D, "2D transient elasticity." },
		{ "TRANSIENT_ELASTICITY_3D", PHYSICS::TRANSIENT_ELASTICITY_3D, "3D transient elasticity." },
		{ "ADVECTION_DIFFUSION_2D" , PHYSICS::ADVECTION_DIFFUSION_2D , "2D advection diffusion"},
		{ "ADVECTION_DIFFUSION_3D" , PHYSICS::ADVECTION_DIFFUSION_3D , "3D advection diffusion"},
		{ "STOKES"                 , PHYSICS::STOKES                 , "Stokes"}
	}));

	SUBCONFIG(Environment       , env         , "Environment dependent variables (set by ./env/scripts).");
	SUBCONFIG(ESPRESOGenerator  , generator   , "ESPRESO internal mesh generator.");
	SUBCONFIG(FETISolver        , feti        , "Internal FETI solver options.");
	SUBCONFIG(MULTIGRIDSolver   , multigrid   , "Multigrid setting.");
	SUBCONFIG(Output            , output      , "Output settings.");

	SUBVECTOR(MaterialParameters, materials   , "Vector of materials (counterd from 1).", "1", "Description of material with index 1");
	SUBMAP(size_t, std::string  , material_set, "Assign materials to regions", "<MATERIAL_INDEX>", "<REGION>");

	SUBMAP(std::string, std::string, displacement       , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, normal_presure     , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, initial_temperature, "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, temperature        , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, heat_source        , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, translation_motions, "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, acceleration       , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, thickness          , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, obstacle           , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
	SUBMAP(std::string, std::string, normal_direction   , "<REGION> <EXPRESSION>;", "<REGION>", "<EXPRESSION>");
};

extern GlobalConfiguration configuration;

}


#endif /* SRC_CONFIG_DESCRIPTION_H_ */
