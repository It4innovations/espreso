
#include "feti.h"
#include "config/configuration.hpp"
#include "esinfo/envinfo.h"

espreso::FETIConfiguration::FETIConfiguration()
{
	method = METHOD::TOTAL_FETI;
	REGISTER(method, ECFMetaData()
			.setdescription({ "Type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("TOTAL_FETI").setdescription("FETI with Dirichlet in B1."))
			.addoption(ECFOption().setname("HYBRID_FETI").setdescription("HYBRID FETI with Dirichlet in B1.")));

	preconditioner = PRECONDITIONER::DIRICHLET;
	REGISTER(preconditioner, ECFMetaData()
			.setdescription({ "Preconditioner" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NONE").setdescription("Without precodition."))
			.addoption(ECFOption().setname("LUMPED").setdescription("Lumped precodition."))
			.addoption(ECFOption().setname("WEIGHT_FUNCTION").setdescription("Precondition by weight function."))
			.addoption(ECFOption().setname("DIRICHLET").setdescription("Dirichler precodition."))
			.addoption(ECFOption().setname("SUPER_DIRICHLET").setdescription("Diagonal Dirichlet precodition.")));

	precision = 1e-5;
	REGISTER(precision, ECFMetaData()
			.setdescription({ "Precision" })
			.setdatatype({ ECFDataType::FLOAT }));

	max_iterations = 200;
	REGISTER(max_iterations, ECFMetaData()
			.setdescription({ "Max iterations" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	num_directions = 3;
	REGISTER(num_directions, ECFMetaData()
			 .setdescription({ "Number of plane wave directions" })
			 .setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	iterative_solver = ITERATIVE_SOLVER::PCG;
	REGISTER(iterative_solver, ECFMetaData()
			.setdescription({ "Iterative solver" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("PCG").setdescription("Projected conjugate gradients."))
			.addoption(ECFOption().setname("pipePCG").setdescription("Pipelined PCG."))
			.addoption(ECFOption().setname("orthogonalPCG").setdescription("Orthogonal PCG."))
			.addoption(ECFOption().setname("GMRES").setdescription("GMRES - supports non-symmetric systems."))
			.addoption(ECFOption().setname("BICGSTAB").setdescription("BICGSTAB - supports non-symmetric systems."))
			.addoption(ECFOption().setname("QPCE").setdescription("QPCE"))
			.addoption(ECFOption().setname("orthogonalPCG_CP").setdescription("FETI Geneo with full ortogonalization CG"))
			.addoption(ECFOption().setname("PCG_CP").setdescription("FETI Geneo with regular CG"))
			);

	regularization = REGULARIZATION::ANALYTIC;
	REGISTER(regularization, ECFMetaData()
			.setdescription({ "Regularization" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("ANALYTIC").setdescription("Analytic regularization provided by a particular physics."))
			.addoption(ECFOption().setname("ALGEBRAIC").setdescription("Regularization based on NULL PIVOTS.")));

	stopping_criterion = STOPPING_CRITERION::RELATIVE;
	REGISTER(stopping_criterion, ECFMetaData()
			.setdescription({ "Stopping criterion" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("RELATIVE").setdescription("Relative stopping."))
			.addoption(ECFOption().setname("ABSOLUTE").setdescription("Absolute stopping.")));

	regularization_version = REGULARIZATION_VERSION::FIX_POINTS;
	REGISTER(regularization_version, ECFMetaData()
			.setdescription({ "Regularization version." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("FIX_POINTS").setdescription(""))
			.addoption(ECFOption().setname("EIGEN_VECTORS").setdescription("."))
			.addoption(ECFOption().setname("WAVE_DIRECTIONS").setdescription(".")));

	conjugate_projector = CONJ_PROJECTOR::NONE;
	REGISTER(conjugate_projector, ECFMetaData()
			.setdescription({ "Conjugate projector" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NONE").setdescription("No conjugate projector."))
			.addoption(ECFOption().setname("GENEO").setdescription("Randomly found null pivots of stiffness matrix."))
			.addoption(ECFOption().setname("CONJ_R").setdescription("Conjugate projector for transient problems from pseudo-kernel"))
			.addoption(ECFOption().setname("CONJ_K").setdescription("Conjugate projector for transient problems from stiffness matrix")));

	geneo_size = 6;
	restart_iteration = 10;
	num_restart = 8;

	REGISTER(geneo_size, ECFMetaData()
			.setdescription({ "GENEO vectors" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));
	REGISTER(restart_iteration, ECFMetaData()
			.setdescription({ "Max iterations in each cycle" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));
	REGISTER(num_restart, ECFMetaData()
			.setdescription({ "Number of restarts" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	exhaustive_info = false;
	REGISTER(exhaustive_info, ECFMetaData()
			.setdescription({ "Print all FETI solver input properties." })
			.setdatatype({ ECFDataType::BOOL }));

	orthogonal_K_kernels = false;
	REGISTER(orthogonal_K_kernels, ECFMetaData()
			.setdescription({ "Orthogonal kernels" })
			.setdatatype({ ECFDataType::BOOL }));

	redundant_lagrange = scaling = true;
	REGISTER(redundant_lagrange, ECFMetaData()
			.setdescription({ "Redundant Lagrange" })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(scaling, ECFMetaData()
			.setdescription({ "Scaling" })
			.setdatatype({ ECFDataType::BOOL }));

	B0_type = B0_TYPE::KERNELS;
	REGISTER(B0_type, ECFMetaData()
			.setdescription({ "HFETI averaging" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("CORNERS").setdescription("B0 composed from corners."))
			.addoption(ECFOption().setname("KERNELS").setdescription("B0 composed from K kernel.")));

	use_schur_complement = false;
	REGISTER(use_schur_complement, ECFMetaData()
			.setdescription({ "Use Schur complement for stiffness matrix processing" })
			.setdatatype({ ECFDataType::BOOL }));

	schur_precision = FLOAT_PRECISION::DOUBLE;
	REGISTER(schur_precision, ECFMetaData()
			.setdescription({ "Precision of Schur complement" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("DOUBLE").setdescription("Double precision."))
			.addoption(ECFOption().setname("SINGLE").setdescription("Single precision.")));

	Ksolver = KSOLVER::DIRECT_DP;
	REGISTER(Ksolver, ECFMetaData()
			.setdescription({ "K solver type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("DIRECT_DP").setdescription("Directly with double precision."))
			.addoption(ECFOption().setname("ITERATIVE").setdescription("Iteratively."))
			.addoption(ECFOption().setname("DIRECT_SP").setdescription("Directly with single precision."))
			.addoption(ECFOption().setname("DIRECT_MP").setdescription("Directly with mixed precision.")));

	Ksolver_max_iterations = 1000;
	REGISTER(Ksolver_max_iterations, ECFMetaData()
			.setdescription({ "Number of reiteration steps for single precision direct solver" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	Ksolver_precision = 1e-12;
	REGISTER(Ksolver_precision, ECFMetaData()
			.setdescription({ "Reguested norm for single precision direct solver" })
			.setdatatype({ ECFDataType::FLOAT }));

	F0_precision = F0SOLVER_PRECISION::K_PRECISION;
	REGISTER(F0_precision, ECFMetaData()
			.setdescription({ "Precision of F0 solver" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("K_PRECISION").setdescription("The same precision as K solver."))
			.addoption(ECFOption().setname("DOUBLE").setdescription("Always double precision.")));

	SAsolver = SASOLVER::CPU_DENSE;
	REGISTER(SAsolver, ECFMetaData()
			.setdescription({ "S alfa solver." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("CPU_DENSE").setdescription("Dense solver on CPU."))
			.addoption(ECFOption().setname("ACC_DENSE").setdescription("Dense solver on Accelerator."))
			.addoption(ECFOption().setname("CPU_SPARSE").setdescription("Sparse solver on CPU.")));

	schur_type = MATRIX_STORAGE::GENERAL;
	REGISTER(schur_type, ECFMetaData()
			.setdescription({ "Storage type for schur complement" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("GENERAL").setdescription("Store a full matrix."))
			.addoption(ECFOption().setname("SYMMETRIC").setdescription("Store only triangle.")));

	mp_pseudoinverse = false;
	combine_sc_and_spds = keep_factors = true;
	REGISTER(mp_pseudoinverse, ECFMetaData()
			.setdescription({ "Moore-Penrose Inverse for FETI Solvers" })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(combine_sc_and_spds, ECFMetaData()
			.setdescription({ "Combine usage of SC for Accelerator and sparse direct solver for CPU" })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(keep_factors, ECFMetaData()
			.setdescription({ "Keep factors between iterations." })
			.setdatatype({ ECFDataType::BOOL }));

	sc_size = 200;
	n_mics = 2;
	REGISTER(sc_size, ECFMetaData()
			.setdescription({ "The size of null pivots for analytics regularization" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));
	REGISTER(n_mics, ECFMetaData()
			.setdescription({ "Number of MIC accelerators." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

	load_balancing = load_balancing_preconditioner = true;
	REGISTER(load_balancing, ECFMetaData()
			.setdescription({ "Load balancing of MIC accelerators" })
			.setdatatype({ ECFDataType::BOOL }));
	REGISTER(load_balancing_preconditioner, ECFMetaData()
			.setdescription({ "Load balancing of Dirichlet preconditioner" })
			.setdatatype({ ECFDataType::BOOL }));

	REGISTER(auto_optimization, ECFMetaData().setdescription({ "Automatic optimization of FETI solver parameters." }));

	allowed_gpu_memory_mb = -1;
	REGISTER(allowed_gpu_memory_mb, ECFMetaData()
			.setdescription({ "The size of GPU memory in MB that is allowed for LSC assembly to decrease #LSC on GPU" })
			.setdatatype({ ECFDataType::INTEGER }));
	num_info_objects = 16;
	REGISTER(num_info_objects, ECFMetaData()
			.setdescription({ "The number of cuSparse CSRSM2 Info objects affects the size of GPU buffers for LSC assembly" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));
	gpu_fragmentation_ratio = 1.5;
	REGISTER(gpu_fragmentation_ratio, ECFMetaData()
			.setdescription({ "The ratio of GPU memory allocation limit to number of domains in cluster - temporal, will be removed" })
			.setdatatype({ ECFDataType::FLOAT }));
	num_streams = info::env::OMP_NUM_THREADS;
	REGISTER(num_streams, ECFMetaData()
			.setdescription({ "The number of CUDA streams for iterative solver" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));
}



