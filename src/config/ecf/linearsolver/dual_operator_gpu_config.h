
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_DUAL_OPERATOR_GPU_CONFIG_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_DUAL_OPERATOR_GPU_CONFIG_H_

#include "config/description.h"

namespace espreso {

struct DualOperatorGpuConfig: public ECFDescription {

	enum class CONCURRENCY {
		AUTO,
		PARALLEL,
		SEQ_CONTINUE,
		SEQ_WAIT
	};

	enum class MATRIX_STORAGE {
		AUTO,
		SPARSE,
		DENSE
	};

	enum class TRS1_SOLVE_TYPE {
		AUTO,
		L,
		LHH
	};

	enum class TRS2_SOLVE_TYPE {
		AUTO,
		U,
		UHH
	};

	enum class MATRIX_ORDER {
		AUTO,
		ROW_MAJOR,
		COL_MAJOR
	};

	enum class PATH_IF_HERMITIAN {
		AUTO,
		TRSM,
		HERK
	};

	enum class TRIANGLE_MATRIX_SHARING {
		AUTO,
		PRIVATE,
		SHARED
	};

	enum class QUEUE_COUNT {
		AUTO,
		PER_THREAD,
		PER_DOMAIN
	};

	enum class DEVICE {
		AUTO,
		CPU,
		GPU
	};

	enum class TIMERS {
		NONE,
		BASIC,
		ALL
	};

	enum class MEMORY_INFO {
		NONE,
		BASIC,
		ALL
	};

	CONCURRENCY concurrency_set;
	CONCURRENCY concurrency_update;
	CONCURRENCY concurrency_apply;
	bool synchronize_update;
	MATRIX_STORAGE trs1_factor_storage;
	MATRIX_STORAGE trs2_factor_storage;
	TRS1_SOLVE_TYPE trs1_solve_type;
	TRS2_SOLVE_TYPE trs2_solve_type;
	MATRIX_ORDER trsm_rhs_sol_order;
	PATH_IF_HERMITIAN path_if_hermitian;
	TRIANGLE_MATRIX_SHARING f_sharing_if_hermitian;
	QUEUE_COUNT queue_count;
	DEVICE apply_scatter_gather_where;
	DEVICE transpose_where;
	TIMERS timers;
	MEMORY_INFO memory_info;

	DualOperatorGpuConfig();

};

}



#endif /* SRC_CONFIG_ECF_LINEARSOLVER_DUAL_OPERATOR_GPU_CONFIG_H_ */
