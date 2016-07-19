
#ifndef ESCONFIG_H_
#define ESCONFIG_H_

#include <cstdlib>
#include <string>
#include <vector>
#include <map>

namespace espreso {

namespace input {
class Description;
}

namespace config {



enum FetiMethod {
	TOTAL_FETI,
	HYBRID_FETI
};

enum B0Type {
	CORNERS,
	KERNELS
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

enum SAsolver {
	SA_DENSE_on_CPU = 0,
	SA_DENSE_on_ACC = 1,
	SA_SPARSE_on_CPU = 2
};

enum SCPrecision {
	SC_DOUBLE_PRECISION,
	SC_SINGLE_PRECISION
};

enum MatrixType {
	GENERAL,
	SYMMETRIC
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

	extern std::vector<input::Description> description;
};

namespace output {

	enum class OUTPUT_FORMATAlternatives {
		VTK_LEGACY = 0,
		VTK_BINARY = 1,
		VTK_MULTIBLOCK = 2,
		ENSIGHT = 3
	};
	/// Format of output data
	extern OUTPUT_FORMATAlternatives OUTPUT_FORMAT;

	/// All results are compressed by 'z' library
	extern bool OUTPUT_COMPRESSION;

	/// Mesh is decimated by this ratio
	extern double OUTPUT_DECIMATION;

	/// Save input to VTK files.
	extern bool SAVE_MESH;

	/// Save fix points to VTK files.
	extern bool SAVE_FIX_POINTS;

	/// Save sub-domains common edges to VTK files.
	extern bool SAVE_EDGES;

	extern std::vector<input::Description> description;
};

namespace assembler {
	enum Discretization { FEM, BEM, API };
	extern int discretization;

	enum Physics { LinearElasticity, Temperature, TransientElasticity };
	extern int physics;

	extern size_t timeSteps;

namespace assembler {
	enum class DISCRETIZATIONalternative {
		/// Finite Element Method
		FEM = 0,
		/// Boundary Element Method
		BEM = 1
	};
	/// Discretization of an example.
	/**
	 * Stiffness matrices are computed based on the discretization:
	 * - Finite Element Method is used
	 * - Boundary Element Method is used
	 */
	extern DISCRETIZATIONalternative DISCRETIZATION;

	enum class DOFS_ORDERalternative {
		/// Group elements - x1, y1, z1, x2, y2, z2, ....
		GROUP_ELEMENTS = 0,

		/// Group elements - x1, x2, ..., y1, y2, ..., z1, z2, ....
		GROUP_DOFS = 1
	};

	extern DOFS_ORDERalternative DOFS_ORDER;
};

namespace solver {
	/// ESPRESO checks the norm of the solution if NORM is not zero
	extern double NORM;

	/// The solver requested precision.
	extern double EPSILON;

	/// Maximum iterations for the solver.
	extern size_t ITERATIONS;


	enum class FETI_METHODalternative {
		/// Total FETI
		TOTAL_FETI = 0,
		/// Hybrid Total FETI
		HYBRID_FETI = 1,
		/// Multi-grid Hypre interface
		HYPRE = 2
	};
	/// A variant of FETI method used by the solver
	extern FETI_METHODalternative FETI_METHOD;


	enum class PRECONDITIONERalternative {
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
	/// Used preconditioner
	extern PRECONDITIONERalternative PRECONDITIONER;


	enum class REGULARIZATIONalternative {
		/// Fix points
		FIX_POINTS = 0,
		/// Randomly found null pivots of stiffness matrix
		NULL_PIVOTS = 1
	};

	/// A type of regularization of stiffness matrix.
	/**
	 * In the case of a singular stiffness matrix, regularization has to be applied.
	 * When an example is loaded with mesh, it is possible to use regularization
	 * from fix points. It is faster then regularization from null pivots.
	 * However, null pivots are more general and usable even without mesh
	 * (e.g. when API is used).
	 */
	extern REGULARIZATIONalternative REGULARIZATION;

	/// Use redundant Lagrange multipliers.
	/**
	 * In the case of Hybrid Total FETI, ESPRESO compose gluing matrix for
	 * each cluster. If this option is on, the multipliers from the cluster
	 * gluing matrix will be also in global gluing matrix.
	 */
	extern bool REDUNDANT_LAGRANGE;

	enum class B0_TYPEalternative {
		/// Gluing based on corners
		CORNERS = 0,
		/// Gluing based on kernels of faces
		KERNELS = 1
	};
	/// Type of cluster gluing matrix.
	/**
	 * Sub-domains in each cluster have to be glued together. Gluing can be based
	 * on corners (random nodes at interfaces between sub-domains) or kernels of
	 * faces between sub-domains.
	 */
	extern B0_TYPEalternative B0_TYPE;


	/// Schur complement will be used.
	extern bool USE_SCHUR_COMPLEMENT;

	enum class SCHUR_COMPLEMENT_PRECalternative {
		/// Double precision
		DOUBLE = 0,
		/// Single precision
		SINGLE = 1
	};
	/// Precision of Schur complement
	extern SCHUR_COMPLEMENT_PRECalternative SCHUR_COMPLEMENT_PREC;

	enum class SCHUR_COMPLEMENT_TYPEalternative {
		/// A full matrix is stored
		GENERAL = 0,
		/// Store only triangle
		SYMMETRIC = 1
	};
	extern SCHUR_COMPLEMENT_TYPEalternative SCHUR_COMPLEMENT_TYPE;

	/// Combine usage of SC for Accelerator and Sparse Direct Solver for CPU
	extern bool COMBINE_SC_AND_SPDS;

	/// Keep factors between iterations
	extern bool KEEP_FACTORS;


	enum class CGSOLVERalternative {
		STANDARD = 0,
		PIPELINED = 1,
		FULL_ORTOGONAL = 2,
		GMRES = 3,
		BICGSTAB = 4
	};

	/// A type of conjugate gradient solver
	extern CGSOLVERalternative CGSOLVER;


	enum class KSOLVERalternative {
		/// A direct solver with double precision
		DIRECT_DP = 0,
		/// An iterative solver
		ITERATIVE = 1,
		/// A direct solver with single precision
		DIRECT_SP = 2,
		/// A direct solver with mixed precision
		DIRECT_MP = 3
	};
	/// A type of stiffness matrix solver
	extern KSOLVERalternative KSOLVER;

	/// Number of reiteration steps for single precision direct solver
	extern size_t   KSOLVER_SP_STEPS;
	/// Iterative norm single precision direct solver
	extern double   KSOLVER_SP_NORM;


	enum class F0SOLVERalternative {
		/// The same precision as K solver
		K_PRECISION = 0,
		/// Double precision
		DOUBLE = 1
	};
	/// A type of F0 solver
	extern F0SOLVERalternative F0SOLVER;

	enum class SASOLVERalternative {
		CPU_DENSE = 0,
		ACC_DENSE = 1,
		CPU_SPARSE = 2
	};
	/// A type of S alfa solver
	extern SASOLVERalternative SASOLVER;

	/// Number of used MIC accelerators
	extern size_t N_MICS;

	/// The number of time steps for transient problems
	extern size_t TIME_STEPS;
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
