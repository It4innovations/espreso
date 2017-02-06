
#ifndef SRC_CONFIGURATION_SOLVERESPRESO_H_
#define SRC_CONFIGURATION_SOLVERESPRESO_H_

#include "../configuration/configuration.h"
#include "../configuration/solverespresooptions.h"

namespace espreso {


struct ESPRESOSolver: public Configuration {

	PARAMETER(double, epsilon, "Solver requested precision.", 1e-5);
	PARAMETER(size_t, iterations, "solver iterations", 100);

	OPTION(ESPRESO_METHOD, method, "The FETI method used by ESPRESO.", ESPRESO_METHOD::TOTAL_FETI, OPTIONS({
		{ "TOTAL_FETI" , ESPRESO_METHOD::TOTAL_FETI , "Total FETI." },
		{ "HYBRID_FETI", ESPRESO_METHOD::HYBRID_FETI, "Hybrid Total FETI." }
	}));

	OPTION(ESPRESO_ITERATIVE_SOLVER, solver, "Used iterative solver", ESPRESO_ITERATIVE_SOLVER::PCG, OPTIONS({
		{ "PCG"          , ESPRESO_ITERATIVE_SOLVER::PCG          , "Standard Projected conjugate gradients." },
		{ "PIPEPCG"      , ESPRESO_ITERATIVE_SOLVER::pipePCG      , "Pipelined PCG." },
		{ "ORTHOGONALPCG", ESPRESO_ITERATIVE_SOLVER::orthogonalPCG, "Full ortogonalization PCG." },
		{ "GMRES"        , ESPRESO_ITERATIVE_SOLVER::GMRES        , "GMRES - allows non-symmetric systems." },
		{ "BICGSTAB"     , ESPRESO_ITERATIVE_SOLVER::BICGSTAB     , "BICGSTAB - allows non-symmetric systems." },
		{ "QPCE"         , ESPRESO_ITERATIVE_SOLVER::QPCE         , "QPCE - allows contact." }
	}));

	OPTION(ESPRESO_PRECONDITIONER, preconditioner, "Preconditioner", ESPRESO_PRECONDITIONER::LUMPED, OPTIONS({
		{ "NONE"           , ESPRESO_PRECONDITIONER::NONE           , "Use no preconditioner." },
		{ "LUMPED"         , ESPRESO_PRECONDITIONER::LUMPED         , "Lumber preconditioner." },
		{ "WEIGHT_FUNCTION", ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION, "Use weight function." },
		{ "DIRICHLET"      , ESPRESO_PRECONDITIONER::DIRICHLET      , "Dirichlet preconditioner." },
		{ "SUPER_DIRICHLET", ESPRESO_PRECONDITIONER::SUPER_DIRICHLET, "simplified Dirichlet preconditioner." },
		{ "MAGIC"          , ESPRESO_PRECONDITIONER::MAGIC          , "TODO." }
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

	OPTION(ESPRESO_KSOLVER, Ksolver, "Stiffness matrix solver type", ESPRESO_KSOLVER::DIRECT_DP, OPTIONS({
		{ "DIRECT_DP", ESPRESO_KSOLVER::DIRECT_DP, "Directly with double precision" },
		{ "ITERATIVE", ESPRESO_KSOLVER::ITERATIVE, "Iteratively" },
		{ "DIRECT_SP", ESPRESO_KSOLVER::DIRECT_SP, "Directly with single precision" },
		{ "DIRECT_MX", ESPRESO_KSOLVER::DIRECT_MP, "Directly with mixed precision" }
	}));

	PARAMETER(size_t, Ksolver_iterations, "Number of reiteration steps for single precision direct solver", 1000);
	PARAMETER(double, Ksolver_epsilon   , "Reguested norm for single precision direct solver", 1e-12);

	OPTION(ESPRESO_F0SOLVER_PRECISION, F0_precision, "Precision of F0", ESPRESO_F0SOLVER_PRECISION::K_PRECISION, OPTIONS({
		{ "K_PRECISION", ESPRESO_F0SOLVER_PRECISION::K_PRECISION, "With the same precision as KSOLVER" },
		{ "DOUBLE"     , ESPRESO_F0SOLVER_PRECISION::DOUBLE     , "Always with double precision" }
	}));

	OPTION(ESPRESO_SASOLVER, SAsolver, "S alfa solver type", ESPRESO_SASOLVER::CPU_DENSE, OPTIONS({
		{ "CPU_DENSE" , ESPRESO_SASOLVER::CPU_DENSE , "Dense solver on CPU" },
		{ "ACC_DENSE" , ESPRESO_SASOLVER::ACC_DENSE , "Dense solver on ACC" },
		{ "CPU_SPARSE", ESPRESO_SASOLVER::CPU_SPARSE, "Sparse solver on CPU." }
	}));

	OPTION(MATRIX_STORAGE, schur_type, "Schur complement matrix type.", MATRIX_STORAGE::GENERAL, OPTIONS({
		{ "GENERAL"  , MATRIX_STORAGE::GENERAL, "Store a full matrix." },
		{ "SYMMETRIC", MATRIX_STORAGE::SYMMETRIC, "Store only triangle." }
	}));

	PARAMETER(bool, combine_sc_and_spds, "Combine usage of SC for Accelerator and sparse direct solver for CPU.", true);
	PARAMETER(bool, keep_factors, "Keep factors between iterations.", true);

	PARAMETER(size_t, N_MICS, "Number of MIC accelerators", 2);
	PARAMETER(bool, load_balancing, "Load balancing of MIC accelerators", true);

	PARAMETER(size_t, time_steps, "Number of time steps for transient problems", 1);
};

}

#endif /* SRC_CONFIGURATION_SOLVERESPRESO_H_ */
