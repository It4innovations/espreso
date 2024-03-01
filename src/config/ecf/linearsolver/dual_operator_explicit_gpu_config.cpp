
#include "config/ecf/linearsolver/dual_operator_explicit_gpu_config.h"

#include "config/configuration.hpp"

using namespace espreso;

DualOperatorExplicitGpuConfig::DualOperatorExplicitGpuConfig()
{
	concurrency_set = CONCURRENCY::DEFAULT;
	REGISTER(concurrency_set, ECFMetaData()
			.setdescription({ "Concurrency during set." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("DEFAULT").setdescription("Default, depending on the specific GPU and libraries."))
			.addoption(ECFOption().setname("PARALLEL").setdescription("Parallel factorization and submitting, asynchronous execution on gpu."))
			.addoption(ECFOption().setname("SEQ_CONTINUE").setdescription("Sequential factorization and submitting, asynchronous execution on gpu."))
			.addoption(ECFOption().setname("SEQ_WAIT").setdescription("Sequential factorization and submitting, synchronous gpu execution.")));

	concurrency_update = CONCURRENCY::DEFAULT;
	REGISTER(concurrency_update, ECFMetaData()
			.setdescription({ "Concurrency during update." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("DEFAULT").setdescription("Default, depending on the specific GPU and libraries."))
			.addoption(ECFOption().setname("PARALLEL").setdescription("Parallel factorization and submitting, asynchronous execution on gpu."))
			.addoption(ECFOption().setname("SEQ_CONTINUE").setdescription("Sequential factorization and submitting, asynchronous execution on gpu."))
			.addoption(ECFOption().setname("SEQ_WAIT").setdescription("Sequential factorization and submitting, synchronous gpu execution.")));

	concurrency_apply = CONCURRENCY::DEFAULT;
	REGISTER(concurrency_apply, ECFMetaData()
			.setdescription({ "Concurrency during apply." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("DEFAULT").setdescription("Default, depending on the specific GPU and libraries."))
			.addoption(ECFOption().setname("PARALLEL").setdescription("Parallel submitting, asynchronous execution on gpu."))
			.addoption(ECFOption().setname("SEQ_CONTINUE").setdescription("Sequential submitting, asynchronous execution on gpu."))
			.addoption(ECFOption().setname("SEQ_WAIT").setdescription("Sequential submitting, synchronous gpu execution.")));
	
	synchronize_update = false;
	REGISTER(synchronize_update, ECFMetaData()
			.setdescription({ "Wait at the end of update for all the GPU kernels to finish." })
			.setdatatype({ ECFDataType::BOOL }));

	trsm1_factor_storage = MATRIX_STORAGE::DEFAULT;
	REGISTER(trsm1_factor_storage, ECFMetaData()
			.setdescription({ "Storage of the factor in the first trsm operation." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("DEFAULT").setdescription("Default, depending on the specific GPU and libraries."))
			.addoption(ECFOption().setname("SPARSE").setdescription("Sparse in CSR format."))
			.addoption(ECFOption().setname("DENSE").setdescription("Dense in row-major order.")));

	trsm2_factor_storage = MATRIX_STORAGE::DEFAULT;
	REGISTER(trsm2_factor_storage, ECFMetaData()
			.setdescription({ "Storage of the factor in the second trsm operation." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("DEFAULT").setdescription("Default, depending on the specific GPU and libraries."))
			.addoption(ECFOption().setname("SPARSE").setdescription("Sparse in CSR format."))
			.addoption(ECFOption().setname("DENSE").setdescription("Dense in row-major order.")));

	trsm1_solve_type = TRSM1_SOLVE_TYPE::DEFAULT;
	REGISTER(trsm1_solve_type, ECFMetaData()
			.setdescription({ "Type of the first solve." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("DEFAULT").setdescription("Default, depending on the specific GPU and libraries."))
			.addoption(ECFOption().setname("L").setdescription("Solve LY=X."))
			.addoption(ECFOption().setname("LHH").setdescription("Solve ((L*)*)Y=X.")));

	trsm2_solve_type = TRSM2_SOLVE_TYPE::DEFAULT;
	REGISTER(trsm2_solve_type, ECFMetaData()
			.setdescription({ "Type of the first solve." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("DEFAULT").setdescription("Default, depending on the specific GPU and libraries."))
			.addoption(ECFOption().setname("U").setdescription("Solve UZ=Y."))
			.addoption(ECFOption().setname("UHH").setdescription("Solve ((U*)*)Z=Y.")));

	trsm_rhs_sol_order = MATRIX_ORDER::DEFAULT;
	REGISTER(trsm_rhs_sol_order, ECFMetaData()
			.setdescription({ "Memory order of the rhs and sol matrices in trsm." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("DEFAULT").setdescription("Default, depending on the specific GPU and libraries."))
			.addoption(ECFOption().setname("ROW_MAJOR").setdescription("Row-major."))
			.addoption(ECFOption().setname("COL_MAJOR").setdescription("Col-major.")));

	path_if_hermitian = PATH_IF_HERMITIAN::DEFAULT;
	REGISTER(path_if_hermitian, ECFMetaData()
			.setdescription({ "Code path if the system is hermitian or symmetric." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("DEFAULT").setdescription("Default, depending on the specific GPU and libraries."))
			.addoption(ECFOption().setname("TRSM").setdescription("Solve LY=x, solve UZ=Y, sp-dn multiply F=B*Z."))
			.addoption(ECFOption().setname("HERK").setdescription("Solve LY=x, dn-dn multiply F=Yt*Y.")));
	
	f_sharing_if_hermitian = TRIANGLE_MATRIX_SHARING::DEFAULT;
	REGISTER(f_sharing_if_hermitian, ECFMetaData()
			.setdescription({ "Sharing of the F matrix if the system is hermitian or symmetric." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("DEFAULT").setdescription("Default, depending on the specific GPU and libraries."))
			.addoption(ECFOption().setname("PRIVATE").setdescription("Every domain has its own allocation to store F."))
			.addoption(ECFOption().setname("SHARED").setdescription("Two domains share the same allocation to store the triangles of F.")));

	queue_count = QUEUE_COUNT::DEFAULT;
	REGISTER(queue_count, ECFMetaData()
			.setdescription({ "Determines the number of gpu queues (streams) and handles." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("DEFAULT").setdescription("Default, depending on the specific GPU and libraries."))
			.addoption(ECFOption().setname("PER_THREAD").setdescription("One queue and handle per thread."))
			.addoption(ECFOption().setname("PER_DOMAIN").setdescription("One queue and handle per domain.")));

	apply_scatter_gather_where = DEVICE::DEFAULT;
	REGISTER(apply_scatter_gather_where, ECFMetaData()
			.setdescription({ "Where to perform the lambda vector C2Dscatter/D2Cgather operations during apply." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("DEFAULT").setdescription("Default, depending on the specific GPU and libraries."))
			.addoption(ECFOption().setname("CPU").setdescription("On the CPU."))
			.addoption(ECFOption().setname("GPU").setdescription("On the GPU.")));

	timers = TIMERS::NONE;
	REGISTER(timers, ECFMetaData()
			.setdescription({ "Verbosity of timers in the explicit gpu F functions." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NONE").setdescription("No timers are printed."))
			.addoption(ECFOption().setname("BASIC").setdescription("Only the basic timers are printed."))
			.addoption(ECFOption().setname("ALL").setdescription("All timers are printed.")));
}
