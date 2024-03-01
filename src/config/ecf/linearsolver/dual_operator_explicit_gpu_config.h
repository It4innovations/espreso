
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_DUAL_OPERATOR_EXPLICIT_GPU_CONFIG_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_DUAL_OPERATOR_EXPLICIT_GPU_CONFIG_H_

#include "config/description.h"

namespace espreso {

struct DualOperatorExplicitGpuConfig: public ECFDescription {

	enum class CONCURRENCY {
		DEFAULT,
		PARALLEL,
		SEQ_CONTINUE,
		SEQ_WAIT
	};

	enum class MATRIX_STORAGE {
		DEFAULT,
		SPARSE,
		DENSE
	};

	enum class TRSM1_SOLVE_TYPE {
		DEFAULT,
		L,
		LHH
	};

	enum class TRSM2_SOLVE_TYPE {
		DEFAULT,
		U,
		UHH
	};

	enum class MATRIX_ORDER {
		DEFAULT,
		ROW_MAJOR,
		COL_MAJOR
	};

	enum class PATH_IF_HERMITIAN {
		DEFAULT,
		TRSM,
		HERK
	};

	enum class TRIANGLE_MATRIX_SHARING {
		DEFAULT,
		PRIVATE,
		SHARED
	};

	enum class QUEUE_COUNT {
		DEFAULT,
		PER_THREAD,
		PER_DOMAIN
	};

	enum class DEVICE {
		DEFAULT,
		CPU,
		GPU
	};

	enum class TIMERS {
		NONE,
		BASIC,
		ALL
	};

	CONCURRENCY concurrency_set;
	CONCURRENCY concurrency_update;
	CONCURRENCY concurrency_apply;
	bool synchronize_update;
	MATRIX_STORAGE trsm1_factor_storage;
	MATRIX_STORAGE trsm2_factor_storage;
	TRSM1_SOLVE_TYPE trsm1_solve_type;
	TRSM2_SOLVE_TYPE trsm2_solve_type;
	MATRIX_ORDER trsm_rhs_sol_order;
	PATH_IF_HERMITIAN path_if_hermitian;
	TRIANGLE_MATRIX_SHARING f_sharing_if_hermitian;
	QUEUE_COUNT queue_count;
	DEVICE apply_scatter_gather_where;
	TIMERS timers;

	DualOperatorExplicitGpuConfig();

};

}



#endif /* SRC_CONFIG_ECF_LINEARSOLVER_DUAL_OPERATOR_EXPLICIT_GPU_CONFIG_H_ */
