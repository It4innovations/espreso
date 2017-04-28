
#ifndef SRC_CONFIGURATION_SOLVER_ESPRESOOPTIONS_H_
#define SRC_CONFIGURATION_SOLVER_ESPRESOOPTIONS_H_

namespace espreso {

enum class ESPRESO_METHOD {
	/// Total FETI
	TOTAL_FETI = 0,
	/// Hybrid Total FETI
	HYBRID_FETI = 1,
};

enum class ESPRESO_ITERATIVE_SOLVER {
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
	QPCE = 5,
	/// FETI Geneo with full ortogonalization CG
	orthogonalPCG_CP = 6,
	/// FETI Geneo with regular CG
	PCG_CP = 7
};

enum class ESPRESO_PRECONDITIONER {
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

enum class CONJ_PROJECTOR {
	/// No conj projector
	NONE = 0,
	/// Randomly found null pivots of stiffness matrix
	GENEO = 1
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

enum class ESPRESO_KSOLVER {
	/// A direct solver with double precision
	DIRECT_DP = 0,
	/// An iterative solver
	ITERATIVE = 1,
	/// A direct solver with single precision
	DIRECT_SP = 2,
	/// A direct solver with mixed precision
	DIRECT_MP = 3
};

enum class ESPRESO_F0SOLVER_PRECISION {
	/// The same precision as K solver
	K_PRECISION = 0,
	/// Double precision
	DOUBLE = 1
};

enum class ESPRESO_SASOLVER {
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

}


#endif /* SRC_CONFIGURATION_SOLVER_ESPRESOOPTIONS_H_ */
