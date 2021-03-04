
#include "clusterGPU.h"
#include "basis/utilities/parser.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"
#include <sstream>
#include <iostream>
#include <algorithm>
#include "mkl.h"
#include "wrappers/mpi/communication.h"
#include "wrappers/cuda/w.cuda.h"
#include "wrappers/csparse/w.csparse.h"

// TODO: move to CUDA wrapper
#include "wrappers/cuda/helper_cuda.h"
#include "cudakernels.h"
#include "wrappers/nvtx/w.nvtx.h"

using namespace espreso;

//#define SHARE_SC

ClusterGPU::~ClusterGPU() {
	DestroyCudaStreamPool();

	for(size_t d = 0; d < domains_in_global_index.size(); d++) {
        if(domains[d].B1Kplus.is_on_acc) {
			domains[d].B1Kplus.ClearCUDA_Stream();
#ifdef SHARE_SC
			if(domains[d].B1Kplus.USE_FLOAT) {
				// Return the original pointer from backup
				domains[d].B1Kplus.d_dense_values_fl = SC_dense_val_orig_fl[d];

				domains[d].B1Kplus.FreeFromCUDA_Dev_fl();
			} else {
				// Return the original pointer from backup
				domains[d].B1Kplus.d_dense_values = SC_dense_val_orig[d];

				domains[d].B1Kplus.FreeFromCUDA_Dev();
			}
#else
			if(domains[d].B1Kplus.USE_FLOAT) {
				domains[d].B1Kplus.FreeFromCUDA_Dev_fl();
			} else {
				domains[d].B1Kplus.FreeFromCUDA_Dev();
			}
		}
#endif			
		if(domains[d].Prec.is_on_acc) {
			domains[d].Prec.ClearCUDA_Stream();
			
			if(domains[d].Prec.is_on_acc)
			{
				domains[d].Prec.FreeFromCUDA_Dev();
			}
		}

		if(domains[d].Prec.is_on_acc) {
			domains[d].Prec.ClearCUDA_Stream();
			if(domains[d].Prec.USE_FLOAT) {
				domains[d].Prec.FreeFromCUDA_Dev_fl();
			} else {
				domains[d].Prec.FreeFromCUDA_Dev();
			}
		}
	}
}

void ClusterGPU::GetGPU() {
	// If already set dont do anything
	if(device_id < 0) {
//		bool GPU_full = false;
		//GPU_full = true;
		int nDevices;
		checkCudaErrors(cudaGetDeviceCount(&nDevices));

		// TODO_GPU
		// - zde se rohoduje, na ktere GPU tento MPI proces pouziva
		// Faze 1 - 1 MPI process pouziva 1 GPU
		//		  - napsat kod, ktere si detekuje kolim MPI ranku je na uzlu a podle toho priradi min. 1 nebo vice MPI procesu na kazde GPU
		// Faze 2 - napsat podporu pro vice GPU na 1 MPI process

		// GPU memory management
		// Create new communicator within the node (OMPI_COMM_TYPE_NODE can be swapped out with MPI_COMM_TYPE_SHARED for portability)
//		MPI_Comm node_comm;
//		MPI_Comm_split_type(info::mpi::comm, MPI_COMM_TYPE_SHARED, info::mpi::rank, MPI_INFO_NULL, &node_comm);
//
//		// Get local size and id
//		int local_procs;
//		MPI_Comm_size(node_comm, &local_procs);
//
//		int local_id;
//		MPI_Comm_rank(node_comm, &local_id);

<<<<<<< HEAD
		size_t procs_per_gpu;
=======
	size_t procs_per_gpu;
>>>>>>> ENH #56: device_id in SparseMatrix

		if(MPITools::node->size > nDevices)
		{
			if ((MPITools::node->size % nDevices) != 0 )
			{
				eslog::error("Only integer multiply number of processes per GPU. Processes: %d, GPUs:  %d\n", MPITools::node->size, nDevices);
			}
			else
			{
				procs_per_gpu = MPITools::node->size / nDevices;
				device_id     = MPITools::node->rank    / procs_per_gpu;
			}
		}
		else
		{
			procs_per_gpu = 1;
			device_id     = MPITools::node->rank;
		}

		checkCudaErrors(cudaSetDevice(device_id));
		checkCudaErrors(cudaMemGetInfo(&GPU_free_mem, &GPU_total_mem));
		GPU_free_mem  /= procs_per_gpu;
		GPU_total_mem /= procs_per_gpu;

		/*OVERKILL PART 1
		// Get memory info for all devices
		std::vector <size_t>  GPU_free_mem(nDevices, 0);
		std::vector <size_t>  GPU_total_mem(nDevices, 0);

		if(local_id == 0)
		{
			for (int i = 0; i < nDevices; i++) {
				cudaSetDevice(i);
				size_t free, total;
				cudaMemGetInfo(&free, &total);
				GPU_free_mem[i] = free;
				GPU_total_mem[i] = total;
			}
		}

		// Assign device
		int device_id = local_id % nDevices;
		cudaSetDevice(device_id);

		// Get mapping from proc to device_id
		std::vector <int>  GPU_mapping(local_procs, 0);
		MPI_Gather(&device_id, 1, MPI_INT, &GPU_mapping[0], 1, MPI_INT, 0, node_comm);

		esint domains_on_GPU = 0;
		esint domains_on_CPU = 0;
		esint DOFs_GPU = 0;
		esint DOFs_CPU = 0;

		size_t local_SC_size_to_add = 0;
		std::vector <int> SC_size_to_add(local_procs, 0);
		std::vector <int> SC_total_size(nDevices, 0);
		std::vector <int> msg(local_procs, 0);
		int reply = 0;

		MPI_Request msg_request;
		MPI_Status msg_status;
		*/
	}
}


size_t ClusterGPU::CalculateGpuBufferSize(esint max_B1_nnz, esint max_B1_rows, esint max_B1_size, esint max_K_rows,
 esint max_L_nnz, esint max_U_nnz) {
	size_t total_buffers_size = 0;

	// 88 MB observed on Barbora with CUDA 11
	// TODO implement #92
	size_t cusparse_size = 100 * 1024 * 1024;
	
	size_t csr_B1_size = max_B1_nnz * sizeof(int) + (max_B1_rows + 1) * sizeof(int) + max_B1_nnz * sizeof(double);
	size_t csr_B1t_size = max_B1_nnz * sizeof(int) + (max_K_rows + 1) * sizeof(int) + max_B1_nnz * sizeof(double);
	size_t dense_X_size = max_B1_size * sizeof(double);

	// Forward and backward permutation vectors
	size_t perm_vectors = 2 * max_K_rows * sizeof(int);

	size_t max_LU_nnz = std::max(max_L_nnz, max_U_nnz);
	size_t csrsm2_buffer_size = ((max_K_rows + max_LU_nnz) * 4 + max_B1_size) * sizeof(double);
	// L and U factors
	size_t csr_L_size = max_L_nnz * sizeof(int) + (max_K_rows + 1) * sizeof(int) + max_L_nnz * sizeof(double);
	size_t csr_U_size = max_U_nnz * sizeof(int) + (max_K_rows + 1) * sizeof(int) + max_U_nnz * sizeof(double);
	// Info objects for L and U factors
	size_t csrsm2_info_size = 2 * max_K_rows * sizeof(int) + csr_L_size + csr_U_size;

	// TODO: Temporal solution - to be removed after proper cudaMalloc error handling implemented for cusparseAnalysis routines
	size_t gpu_fragmentation_limit = configuration.gpu_fragmentation_ratio * lsc_on_gpu_ids.size() * 1024 * 1024;

	total_buffers_size = n_streams_per_gpu * (csr_B1_size + csr_B1t_size + dense_X_size + csrsm2_buffer_size + perm_vectors)
	 + n_csrsm2_info_per_gpu * csrsm2_info_size + csr_L_size + csr_U_size + cusparse_size + gpu_fragmentation_limit;

	return total_buffers_size;
}


void ClusterGPU::Create_SC_perDomain(bool USE_FLOAT) {
<<<<<<< HEAD
	DEBUGOUT << "Creating B1*K+*B1t Schur Complements with CSparse and cuSparse on GPU\n";
	
	// Currently sets device_id for the cluster (all domains)
	GetGPU();

	void* d_blocked_memory;
	size_t blocked_memory_size = 0;
	if (configuration.allowed_gpu_memory_mb > 0)
		blocked_memory_size = GPU_free_mem - ((size_t) configuration.allowed_gpu_memory_mb * 1024 * 1024);
	cuda::Malloc((void **) &d_blocked_memory, blocked_memory_size);

	// TODO Initialize cusparse context
	
	// Update amount of free device memory
	checkCudaErrors(cudaMemGetInfo(&GPU_free_mem, &GPU_total_mem));

	// Read configuration
	n_csrsm2_info_per_gpu = configuration.num_info_objects;

	// Decide if LSC fits in GPU memory
=======
	// Currently sets device_id for the cluster (all domains)
	GetAvailableGPUmemory();
	std::cout << "Creating B1*K+*B1t Schur Complements with Pardiso SC and coping them to GPU";

>>>>>>> ENH #56: device_id in SparseMatrix
	esint status = 0;
	std::vector<size_t> local_SC_size_to_add(domains_in_global_index.size(), 0);

	esint domains_on_GPU = 0;
	esint domains_on_CPU = 0;
	esint DOFs_GPU = 0;
	esint DOFs_CPU = 0;

	SEQ_VECTOR<esint> lsc_to_get_factors_ids;
	// Reserve capacity to prevent reallocation (full utilization of the capacity is expected in the most cases)
	lsc_to_get_factors_ids.reserve(domains_in_global_index.size());

	// Calculate the size of LSC matrix and input and output vectors per domain
	size_t gpu_buffers_size = 0;
	// Calculate GPU buffers without factors
	esint max_B1_nnz = 0;
	esint max_B1_rows = 0;
	esint max_B1_size = 0;
	esint max_K_rows = 0;
	esint max_L_nnz = 0;
	esint max_U_nnz = 0;

	// First round - Select domains for factorization to get max factor sizes
	for(size_t d = 0; d < domains_in_global_index.size(); d++) {
		size_t vec_size = (domains[d].B1_comp_dom.rows > domains[d].B1t_Dir_perm_vec.size())
				? domains[d].B1_comp_dom.rows
				: domains[d].B1t_Dir_perm_vec.size();

		switch (configuration.schur_type) {
		case FETIConfiguration::MATRIX_STORAGE::GENERAL:

	// TODO_GPU - SHARE_SC asi nepatri sem pro general matrix, ale patri do symetrickych matic nize 
	// Radime potvrd.
	// RV: Ten switch se rozhoduje na zaklade promenne schur_type, ktera se nastavuje v ecf souboru.
	// Urcuje se tim, zda se pro symetricky systemu bude SC ukladat do ctevrcove matice, nebo do trojuhelniku.
	// V pripade cterce se pak na zaklade makra SHARE_SC urci, zda se do jedne matice ulozi 1 nebo 2 domeny.
	// Osetreni aby se tohle cele delalo jen v pripade symetrickeho systemu zatim bud neni, nebo je jinde - je potreba overit.
#ifdef SHARE_SC
			// SC_total_size will be halved in the case of 2 symmetric SCs in 1 full matrix 
			if (d%2 == 0) {
				esint sc1_rows = domains[d].B1_comp_dom.rows;
				esint sc2_rows = 0;
				esint SC_size = 0;
				esint vec_size = domains[d].B1_comp_dom.rows;

				if(d+1 < domains_in_global_index.size()) {
					sc2_rows = domains[d+1].B1_comp_dom.rows;
					vec_size += domains[d+1].B1_comp_dom.rows;
				}

				if(sc1_rows > sc2_rows) {
					SC_size = sc1_rows * sc1_rows;
//					vec_size = sc1_rows;
				} else if(sc1_rows == sc2_rows) {
					SC_size = sc1_rows * sc2_rows;
//					vec_size = sc1_rows;
				} else { // sc1_rows < sc2_rows
					SC_size = sc2_rows * sc2_rows;
//					vec_size = sc2_rows;
				}

				if (USE_FLOAT) {
					local_SC_size_to_add[d] = (SC_size + 2 * vec_size) * sizeof(float);
				} else {
					local_SC_size_to_add[d] = (SC_size + 2 * vec_size) * sizeof(double);
				}
			}
#else
			if (USE_FLOAT) {
				local_SC_size_to_add[d] =
						( domains[d].B1_comp_dom.rows * domains[d].B1_comp_dom.rows +
						  2 * vec_size
						) * sizeof(float);
			} else {
				local_SC_size_to_add[d] =
						( domains[d].B1_comp_dom.rows * domains[d].B1_comp_dom.rows +
						  2 * vec_size
						) * sizeof(double);
			}
#endif
			break;
		case FETIConfiguration::MATRIX_STORAGE::SYMMETRIC:
			if (USE_FLOAT) {
				local_SC_size_to_add[d] =
						(((domains[d].B1_comp_dom.rows + 1 ) * domains[d].B1_comp_dom.rows ) / 2
						 + 2 * vec_size
						) * sizeof(float);
			} else {
				local_SC_size_to_add[d] =
						(((domains[d].B1_comp_dom.rows + 1 ) * domains[d].B1_comp_dom.rows ) / 2
						 + 2 * vec_size
						) * sizeof(double);
			}
			break;
		default:
			break;
			eslog::error("ERROR - Not implemented type of Schur complement.");
		}

<<<<<<< HEAD
		max_B1_nnz = std::max(max_B1_nnz, domains[d].B1_comp_dom.nnz);
		max_B1_rows = std::max(max_B1_rows, domains[d].B1_comp_dom.rows);
		max_B1_size = std::max(max_B1_size, domains[d].B1_comp_dom.rows * domains[d].B1_comp_dom.cols);
		max_K_rows = std::max(max_K_rows, domains[d].K.rows);

		// max_L_nnz and max_U_nnz are 0 now
		gpu_buffers_size = CalculateGpuBufferSize(max_B1_nnz, max_B1_rows, max_B1_size, max_K_rows, max_L_nnz, max_U_nnz);
		
		// Possible optimization - overlap gpu_buffers_size with size of vectors, see #86
		if(local_SC_size_to_add[d] < (gpu_buffers_size > GPU_free_mem ? 0 : GPU_free_mem - gpu_buffers_size)) {
=======
		if(local_SC_size_to_add[d] < GPU_free_mem)
		{
			domains_on_GPU++;
			DOFs_GPU += domains[d].K.rows;
			domains[d].B1Kplus.is_on_acc = 1;
			// Prepared for the case of multi-GPU per cluster
			domains[d].B1Kplus.device_id = device_id;
>>>>>>> ENH #56: device_id in SparseMatrix
			GPU_free_mem -= local_SC_size_to_add[d];
			lsc_to_get_factors_ids.push_back(d);
		} else {
			domains_on_CPU++;
			DOFs_CPU += domains[d].K.rows;
			domains[d].B1Kplus.is_on_acc = 0;
			lsc_on_cpu_ids.push_back(d);
		}

		/* OVERKILL PART 2
		// Collect info on how much each process wants to add
		MPI_Gather(&local_SC_size_to_add, 1, MPI_INT, &SC_size_to_add[0], 1, MPI_INT, 0, node_comm);

		if(local_id == 0)
		{   // Proc with id 0 will decide if it fits
			for(int proc = 0; proc < local_procs; proc++)
			{
				if (SC_total_size[GPU_mapping[proc]] + SC_size_to_add[proc] < GPU_free_mem[GPU_mapping[proc]]) {
					msg[proc] = 1;
					SC_total_size[GPU_mapping[proc]] += SC_size_to_add[proc];
				} else {
					msg[proc] = 0;
				}
			}
		}
		// Get the decision from proc 0
		MPI_Scatter(&msg[0], 1, MPI_INT, &reply, 1, MPI_INT, 0, node_comm);

		if(reply)
		{
			domains_on_GPU++;
			DOFs_GPU += domains[d].K.rows;
                        domains[d].is_on_acc = 1;
		}
		else
		{
			domains_on_CPU++;
			DOFs_CPU += domains[d].K.rows;
                        domains[d].is_on_acc = 0;
		}
		*/
	}

	// Factorize and get factor sizes, process only domains that can fit GPU w/o factor-based buffers
	PUSH_RANGE("Fact+LSC mem alloc", 1)
	int order;
	SEQ_VECTOR<int> vec_L_nnz(lsc_to_get_factors_ids.size());
	SEQ_VECTOR<int*> vec_L_row_indexes(lsc_to_get_factors_ids.size());
	SEQ_VECTOR<int*> vec_L_col_pointers(lsc_to_get_factors_ids.size());
	SEQ_VECTOR<double*> vec_L_values(lsc_to_get_factors_ids.size());
	SEQ_VECTOR<int*> vec_perm(lsc_to_get_factors_ids.size());

	SEQ_VECTOR<int> vec_U_nnz;
	SEQ_VECTOR<int*> vec_U_row_indexes;
	SEQ_VECTOR<int*> vec_U_col_pointers;
	SEQ_VECTOR<double*> vec_U_values;
	SEQ_VECTOR<int*> vec_perm_2;

	SEQ_VECTOR<int> vec_K_rows(lsc_to_get_factors_ids.size());
	SEQ_VECTOR<int> vec_B1_rows(lsc_to_get_factors_ids.size());
	SEQ_VECTOR<int> vec_B1_nnz(lsc_to_get_factors_ids.size());
	SEQ_VECTOR<int> vec_B1_size(lsc_to_get_factors_ids.size());

	if(SYMMETRIC_SYSTEM) {
		order = 1; // 0 = natural, 1 = amd(A+A')

		#pragma omp parallel for
		for (size_t d = 0; d < lsc_to_get_factors_ids.size(); d++) {
			esint idx = lsc_to_get_factors_ids[d];

			PUSH_RANGE("Fact", 1)
			csparse::FactorizeChol(domains[idx].K, order, vec_L_nnz[d], vec_L_row_indexes[d], 
			vec_L_col_pointers[d], vec_L_values[d], vec_perm[d]);
			POP_RANGE

			vec_K_rows[d] = domains[idx].K.rows;
			vec_B1_rows[d] = domains[idx].B1_comp_dom.rows;
			vec_B1_nnz[d] = domains[idx].B1_comp_dom.CSR_V_values.size();
			vec_B1_size[d] = domains[idx].B1_comp_dom.rows * domains[idx].B1_comp_dom.cols;
		}
	} else {
		vec_U_nnz.resize(lsc_to_get_factors_ids.size());
		vec_U_row_indexes.resize(lsc_to_get_factors_ids.size());
		vec_U_col_pointers.resize(lsc_to_get_factors_ids.size());
		vec_U_values.resize(lsc_to_get_factors_ids.size());
		vec_perm_2.resize(lsc_to_get_factors_ids.size());

		order = 1; // 0 = natural, 1 = amd(A+A'), 2 = amd(S'*S), 3 = amd(A'*A)

		#pragma omp parallel for
		for (size_t d = 0; d < lsc_to_get_factors_ids.size(); d++) {
			esint idx = lsc_to_get_factors_ids[d];

			PUSH_RANGE("Fact", 1)
			csparse::FactorizeLu(domains[idx].K, order, vec_L_nnz[d], vec_L_row_indexes[d], 
			vec_L_col_pointers[d], vec_L_values[d], vec_U_nnz[d], vec_U_row_indexes[d], 
			vec_U_col_pointers[d], vec_U_values[d], vec_perm[d], vec_perm_2[d]);
			POP_RANGE

			vec_K_rows[d] = domains[idx].K.rows;
			vec_B1_rows[d] = domains[idx].B1_comp_dom.rows;
			vec_B1_nnz[d] = domains[idx].B1_comp_dom.CSR_V_values.size();
			vec_B1_size[d] = domains[idx].B1_comp_dom.rows * domains[idx].B1_comp_dom.cols;
		}
	}
	POP_RANGE

	// Reset max variables
	max_B1_nnz = 0;
	max_B1_rows = 0;
	max_B1_size = 0;
	max_K_rows = 0;
	checkCudaErrors(cudaMemGetInfo(&GPU_free_mem, &GPU_total_mem));

	// Second round - Assign domains to GPUs wrt max factor sizes
	// Reserve capacity to prevent reallocation (full utilization of the capacity is expected in the most cases)
	lsc_on_gpu_ids.reserve(lsc_to_get_factors_ids.size());
	for(size_t d = 0, deleted = 0; d < lsc_to_get_factors_ids.size(); d++) {
		esint idx = lsc_to_get_factors_ids[d];

		max_B1_nnz = std::max(max_B1_nnz, domains[idx].B1_comp_dom.nnz);
		max_B1_rows = std::max(max_B1_rows, domains[idx].B1_comp_dom.rows);
		max_B1_size = std::max(max_B1_size, domains[idx].B1_comp_dom.rows * domains[idx].B1_comp_dom.cols);
		max_K_rows = std::max(max_K_rows, domains[idx].K.rows);
		max_L_nnz = std::max(max_L_nnz, vec_L_nnz[d]);
		max_U_nnz = SYMMETRIC_SYSTEM ? max_L_nnz : std::max(max_U_nnz, vec_U_nnz[d]);

		// max_L_nnz and max_U_nnz are 0 now
		gpu_buffers_size = CalculateGpuBufferSize(max_B1_nnz, max_B1_rows, max_B1_size, max_K_rows, max_L_nnz, max_U_nnz);
		
		// Possible optimization - overlap gpu_buffers_size with size of vectors, see #86
		if(local_SC_size_to_add[idx] < (gpu_buffers_size > GPU_free_mem ? 0 : GPU_free_mem - gpu_buffers_size)) {
			domains_on_GPU++;
			DOFs_GPU += domains[idx].K.rows;
			domains[idx].B1Kplus.is_on_acc = 1;
			// Prepared for the case of multi-GPU per cluster
			domains[idx].B1Kplus.device_id = device_id;
			GPU_free_mem -= local_SC_size_to_add[idx];
			lsc_on_gpu_ids.push_back(idx);
		} else {
			domains_on_CPU++;
			DOFs_CPU += domains[idx].K.rows;
			domains[idx].B1Kplus.is_on_acc = 0;
			lsc_on_cpu_ids.push_back(idx);

			if (SYMMETRIC_SYSTEM) {
				csparse::FreeCholFactor(vec_L_row_indexes[d - deleted], vec_L_col_pointers[d - deleted], vec_L_values[d - deleted], vec_perm[d - deleted]);

				vec_L_nnz.erase(vec_L_nnz.begin() + (d - deleted));
				vec_L_row_indexes.erase(vec_L_row_indexes.begin() + (d - deleted));
				vec_L_col_pointers.erase(vec_L_col_pointers.begin() + (d - deleted));
				vec_L_values.erase(vec_L_values.begin() + (d - deleted));
				vec_perm.erase(vec_perm.begin() + (d - deleted));
			} else {
				csparse::FreeLuFactors(vec_L_row_indexes[d - deleted], vec_L_col_pointers[d - deleted], vec_L_values[d - deleted],
				vec_U_row_indexes[d - deleted], vec_U_col_pointers[d - deleted], vec_U_values[d - deleted], vec_perm[d - deleted], vec_perm_2[d - deleted]);

				vec_L_nnz.erase(vec_L_nnz.begin() + (d - deleted));
				vec_L_row_indexes.erase(vec_L_row_indexes.begin() + (d - deleted));
				vec_L_col_pointers.erase(vec_L_col_pointers.begin() + (d - deleted));
				vec_L_values.erase(vec_L_values.begin() + (d - deleted));
				vec_perm.erase(vec_perm.begin() + (d - deleted));

				vec_U_nnz.erase(vec_U_nnz.begin() + (d - deleted));
				vec_U_row_indexes.erase(vec_U_row_indexes.begin() + (d - deleted));
				vec_U_col_pointers.erase(vec_U_col_pointers.begin() + (d - deleted));
				vec_U_values.erase(vec_U_values.begin() + (d - deleted));
				vec_perm_2.erase(vec_perm_2.begin() + (d - deleted));

			}
			deleted++;
		}
	}

	// TODO_GPU - vsechny tyto std::cout se musi prepsat na logovani co ma Ondra M. 
	// Ondro nektere moje rutiny, napr. SpyText jsou napsane pro std::cout a ne printf. Jake je reseni ? 
//	std::vector <int> on_gpu (info::mpi::size, 0);
//	MPI_Gather(&domains_on_GPU,1,MPI_INT,&on_gpu[0],1,MPI_INT, 0, info::mpi::comm);
//
//	std::vector <int> on_cpu (info::mpi::size, 0);
//	MPI_Gather(&domains_on_CPU,1,MPI_INT,&on_cpu[0],1,MPI_INT, 0, info::mpi::comm);
//
//	std::vector <int> don_gpu (info::mpi::size, 0);
//	MPI_Gather(&DOFs_GPU,1,MPI_INT,&don_gpu[0],1,MPI_INT, 0, info::mpi::comm);
//
//	std::vector <int> don_cpu (info::mpi::size, 0);
//	MPI_Gather(&DOFs_CPU,1,MPI_INT,&don_cpu[0],1,MPI_INT, 0, info::mpi::comm);
//
//	if (info::mpi::rank == 0) {
//		DEBUGOUT << "Local Schur complement:" << "\n";
//		for (esint i = 0; i < info::mpi::size; i++) {
//			DEBUGOUT << " MPI rank " << i <<
//				"\t - GPU : domains = \t" << on_gpu[i] << "\t Total DOFs = \t" << don_gpu[i] <<
//				"\t - CPU : domains = \t" << on_cpu[i] << "\t Total DOFs = \t" << don_cpu[i] << "\n";
//		}
//	}

	esint info[2] = { domains_on_GPU, domains_on_CPU + domains_on_GPU };
	Communication::allReduce(info, NULL, 2, MPITools::getType<esint>().mpitype, MPI_SUM);
	std::string ratio = Parser::stringwithcommas(info[0]) + " / " + Parser::stringwithcommas(info[1]);
	eslog::solver("     - | ACCELERATED DOMAINS %57s | -\n", ratio.c_str());

// Assemble LSC matrices per domain
#ifdef SHARE_SC
	SEQ_VECTOR <esint> SC_dense_val_offsets(domains_in_global_index.size(), 0);

	// 2 domains per iteration processed
	#pragma omp parallel for
	for (size_t d = 0; d < domains_in_global_index.size(); d += 2 ) {

		if (domains[d].is_on_acc == 1 || !configuration.combine_sc_and_spds) {
			// Calculates SC on CPU and keeps it CPU memory
			GetSchurComplement(USE_FLOAT, d);
			std::cout << Info::plain() << ".";

			// Set if Upper or Lower part is referenced
			domains[d].B1Kplus.uplo = 'U';

			// Set the default lda
			domains[d].B1Kplus.extern_lda = domains[d].B1Kplus.rows;
		}

		if (d+1 < domains_in_global_index.size() && (domains[d+1].B1Kplus.is_on_acc == 1 || !configuration.combine_sc_and_spds)) {
			// Calculates SC on CPU and keeps it CPU memory
			GetSchurComplement(USE_FLOAT, d+1);
//			std::cout << Info::plain() << ".";

			esint sc1_rows = domains[d].B1Kplus.rows;
			esint sc2_rows = domains[d+1].B1Kplus.rows;

			// Set if Upper or Lower part is referenced
			domains[d+1].B1Kplus.uplo = 'L';

			// Both SCs stored in the first domain [d]
			if(sc1_rows > sc2_rows) {
				// First SC -> U
				if (USE_FLOAT) {
					for(esint r = 0; r < sc2_rows; r++) {
						std::copy(&domains[d+1].B1Kplus.dense_values_fl[r*sc2_rows + r], &domains[d+1].B1Kplus.dense_values_fl[(r+1) *sc2_rows],
						 &domains[d].B1Kplus.dense_values_fl[r*sc1_rows + r+1]);
					}
				} else {
					for(esint r = 0; r < sc2_rows; r++) {
						std::copy(&domains[d+1].B1Kplus.dense_values[r*sc2_rows + r], &domains[d+1].B1Kplus.dense_values[(r+1) *sc2_rows],
						 &domains[d].B1Kplus.dense_values[r*sc1_rows + r+1]);
					}
				}

				SC_dense_val_offsets[d] = 0;
				SC_dense_val_offsets[d+1] = 1;

				domains[d+1].B1Kplus.extern_lda = sc1_rows;

			} else if(sc1_rows == sc2_rows) {
				// First SC -> U
				if (USE_FLOAT) {
					SEQ_VECTOR <float> sc1_tmp_fl(sc1_rows * (sc1_rows + 1));

					for(esint r = 0; r < sc1_rows; r++) {
						std::copy(&domains[d+1].B1Kplus.dense_values_fl[r*sc2_rows + r], &domains[d+1].B1Kplus.dense_values_fl[(r+1) *sc2_rows],
						 &sc1_tmp_fl[r*sc1_rows + r]);
						std::copy(&domains[d].B1Kplus.dense_values_fl[r*sc1_rows], &domains[d].B1Kplus.dense_values_fl[r*sc1_rows + r+1],
						 &sc1_tmp_fl[(r+1)*sc1_rows]);
					}
					domains[d].B1Kplus.dense_values_fl.swap(sc1_tmp_fl);
				} else {
					SEQ_VECTOR <double> sc1_tmp(sc1_rows * (sc1_rows + 1));

					for(esint r = 0; r < sc1_rows; r++) {
						std::copy(&domains[d+1].B1Kplus.dense_values[r*sc2_rows + r], &domains[d+1].B1Kplus.dense_values[(r+1) *sc2_rows],
						 &sc1_tmp[r*sc1_rows + r]);
						std::copy(&domains[d].B1Kplus.dense_values[r*sc1_rows], &domains[d].B1Kplus.dense_values[r*sc1_rows + r+1],
						 &sc1_tmp[(r+1)*sc1_rows]);
					}
					domains[d].B1Kplus.dense_values.swap(sc1_tmp);
				}

				SC_dense_val_offsets[d] = sc1_rows;
				SC_dense_val_offsets[d+1] = 0;

				domains[d].B1Kplus.extern_lda = sc1_rows;
				domains[d+1].B1Kplus.extern_lda = sc2_rows;

			} else { // sc1_rows < sc2_rows
				// First SC -> U
				if (USE_FLOAT) {
					for(esint r = 0; r < sc1_rows; r++) {
						std::copy(&domains[d].B1Kplus.dense_values_fl[r*sc1_rows], &domains[d].B1Kplus.dense_values_fl[r*sc1_rows + r + 1],
						 &domains[d+1].B1Kplus.dense_values_fl[(r+1)*sc2_rows]);
					}
					domains[d].B1Kplus.dense_values_fl = std::move(domains[d+1].B1Kplus.dense_values_fl);
				} else {
					for(esint r = 0; r < sc1_rows; r++) {
						std::copy(&domains[d].B1Kplus.dense_values[r*sc1_rows], &domains[d].B1Kplus.dense_values[r*sc1_rows + r + 1],
						 &domains[d+1].B1Kplus.dense_values[(r+1)*sc2_rows]);
//						std::copy(&domains[d].B1Kplus.dense_values[r*sc1_rows + r], &domains[d].B1Kplus.dense_values[(r+1) *sc1_rows],
//						 &domains[d+1].B1Kplus.dense_values[r*sc2_rows + r+1]);
					}
					domains[d].B1Kplus.dense_values = std::move(domains[d+1].B1Kplus.dense_values);
				}

				SC_dense_val_offsets[d] = sc2_rows;
				SC_dense_val_offsets[d+1] = 0;

				domains[d].B1Kplus.extern_lda = sc2_rows;
				domains[d+1].B1Kplus.extern_lda = sc2_rows;
			}
		}
	}

	if (USE_FLOAT){
		SC_dense_val_orig_fl.resize(domains_in_global_index.size());
	} else {
		SC_dense_val_orig.resize(domains_in_global_index.size());
	}
#else
	// TODO refactoring: Possibly remove the last checking phase
	GetSchurComplementsGpu(USE_FLOAT, vec_L_nnz, vec_L_row_indexes, vec_L_col_pointers, vec_L_values,
 	 vec_U_nnz, vec_U_row_indexes, vec_U_col_pointers, vec_U_values, vec_perm, vec_perm_2, max_B1_nnz,
 	 max_B1_rows, max_B1_size, max_K_rows, max_L_nnz, max_U_nnz);

	GetSchurComplementsCpu(USE_FLOAT);
#endif
//	eslog::info("\n");

	// TODO: check if correct GPU is set in case of multi-GPU per cluster
	CreateCudaStreamPool();

	// TODO refactoring: device vectors should be allocated with memory for LSC and the following phase removed
	// Allocate vectors on GPU - now sequential
	// TODO_GPU - Faze 1: tohle je zpusob jak se budou kopirovat i Dirichlet predpodminovace na GPU - inspirovat se, pripadne se muze dat rovnou to teto smycky 
	// TODO_GPU - Faze 2: zjistit, jak je to pomale a pripadne optimalizovat. Ale vyresi se implementaci vypoctu LSC na GPU 
	for (size_t d = 0; d < domains_in_global_index.size(); d++ ) {
		esint status = 0;

		if (domains[d].B1Kplus.is_on_acc == 1) {
			// TODO_GPU: Test performance (same streams)

			// Assign CUDA stream from he pool
			// TODO: This is probably unnecessary for just vectors allocation -> check if can be removed once the USE_FLOAT branch is implemented	
			domains[d].B1Kplus.stream = cuda_stream_pool[d % configuration.num_streams];
			domains[d].B1Kplus.handle = cublas_handle_pool[d % configuration.num_streams];
                        //if(USE_PREC == FETIConfiguration::PRECONDITIONER::DIRICHLET)
//                        {
//                                domains[d].Prec.SetCUDA_Stream(cuda_stream_pool[d % configuration.num_streams]);
//                        }

#ifdef SHARE_SC
			if (USE_FLOAT){
				if(d%2==0) {
					status = domains[d].B1Kplus.CopyToCUDA_Dev_fl();

					// Backup original pointer for freeing
					SC_dense_val_orig_fl[d] = domains[d].B1Kplus.d_dense_values_fl;

					// Add the specific offset for the device memory
					domains[d].B1Kplus.d_dense_values_fl += SC_dense_val_offsets[d];
				} else {
					// Allocate GPU memory only for vectors
					status = domains[d].B1Kplus.MallocVecsOnCUDA_Dev_fl();

					// Backup original pointer for freeing
					SC_dense_val_orig_fl[d] = NULL;

					// Set the pointer to the device memory of the previous domain with the specific offset
					domains[d].B1Kplus.d_dense_values_fl = domains[d-1].B1Kplus.d_dense_values_fl - SC_dense_val_offsets[d-1] + SC_dense_val_offsets[d];
				}

				size_t memsize = (domains[d].B1Kplus.rows > domains[d].B1t_Dir_perm_vec.size()) ? domains[d].B1Kplus.rows : domains[d].B1t_Dir_perm_vec.size();
				// TODO_GPU - alokuji se pole pro vektory, kterymi se bude ve funkci Apply_A nasobit tato matice
				domains[d].ml cuda_pinned_buff_fl.resize(memsize);
				if (domains[d].cuda_pinned_buff_fl.size() != memsize)
				{
				      std::cout << "Error allocating pinned host memory";
				      status = -1;
				}
				cudaHostRegister(domains[d].cuda_pinned_buff_fl.data(), memsize * sizeof(float), cudaHostRegisterDefault);

			} else {
				if(d%2==0) {
					status = domains[d].B1Kplus.CopyToCUDA_Dev();

					// Backup original pointer for freeing
					SC_dense_val_orig[d] = domains[d].B1Kplus.d_dense_values;

					// Add the specific offset for the device memory;
					domains[d].B1Kplus.d_dense_values += SC_dense_val_offsets[d];

				} else {
					// Allocate GPU memory only for vectors
					status = domains[d].B1Kplus.MallocVecsOnCUDA_Dev();

					// Backup original pointer for freeing
					SC_dense_val_orig[d] = NULL;

					// Set the pointer to the device memory of the previous domain with the specific offset
					domains[d].B1Kplus.d_dense_values = domains[d-1].B1Kplus.d_dense_values - SC_dense_val_offsets[d-1] + SC_dense_val_offsets[d];
				}

				size_t memsize = (domains[d].B1Kplus.rows > domains[d].B1t_Dir_perm_vec.size()) ? domains[d].B1Kplus.rows : domains[d].B1t_Dir_perm_vec.size();

				domains[d].cuda_pinned_buff.resize(memsize);
				if (domains[d].cuda_pinned_buff.size() != memsize)
				{
				      std::cout << "Error allocating pinned host memory";
				      status = -1;
				}
				cudaHostRegister(domains[d].cuda_pinned_buff.data(), memsize * sizeof(double), cudaHostRegisterDefault);

			}
#else
			if (USE_FLOAT){
				 status      = domains[d].B1Kplus.CopyToCUDA_Dev_fl();

				size_t memsize = (domains[d].B1_comp_dom.rows > domains[d].B1t_Dir_perm_vec.size()) ? domains[d].B1_comp_dom.rows : domains[d].B1t_Dir_perm_vec.size();
				// TODO_GPU - alokuji se pole pro vektory, kterymi se bude ve funkci Apply_A nasobit tato matice 
				domains[d].cuda_pinned_buff_fl.resize(memsize);			
				if (domains[d].cuda_pinned_buff_fl.size() != memsize) {
				      eslog::error("Error allocating pinned host memory");
				      status = -1;
				}
                cudaHostRegister(domains[d].cuda_pinned_buff_fl.data(), memsize * sizeof(float), cudaHostRegisterDefault);
			} else {
				// LSC matrix is already in device memory at this point 
				// This B1Kplus device vectors will be shared with Prec matrix
				// TODO: Will be moved from SparseMatrix to Domain in #50 or #64
				// TODO allocate `memsize` for shared device vectors after move from SparseMatrix to Domain
				status = domains[d].B1Kplus.MallocVecsOnCUDA_Dev();

				size_t memsize = (domains[d].B1_comp_dom.rows > domains[d].B1t_Dir_perm_vec.size()) ? domains[d].B1_comp_dom.rows : domains[d].B1t_Dir_perm_vec.size();

                domains[d].cuda_pinned_buff.resize(memsize);
				if (domains[d].cuda_pinned_buff.size() != memsize) {
				      eslog::error("Error allocating pinned host memory");
				      status = -1;
				}
                cudaHostRegister(domains[d].cuda_pinned_buff.data(), memsize * sizeof(double), cudaHostRegisterDefault);
			}
#endif
			if (status == 0) {
				// TODO_GPU - LSCs ktere se uspesne prenesly do pameti GPU se smazou z pameti procesoru 
				domains[d].Kplus.keep_factors = false;

				if (USE_FLOAT) {
					SEQ_VECTOR <float>  ().swap (domains[d].B1Kplus.dense_values_fl);

					// eslog::info("g");
				} else {
					SEQ_VECTOR <double> ().swap (domains[d].B1Kplus.dense_values);

					// eslog::info("G");
				}
			} else {
				// pokud se domenu nepodar nahrat na GPU 
				domains[d].B1Kplus.is_on_acc = 0;
				domains_on_CPU++;
				domains_on_GPU--;
				DOFs_CPU += domains[d].K.rows;
				DOFs_GPU -= domains[d].K.rows;
				GPU_free_mem += local_SC_size_to_add[d];

				domains[d].B1Kplus.ClearCUDA_Stream();

				// TODO vyriesit vracanie Prec ak sa nepodary nahrat B1K+
//				domains[d].Prec.ClearCUDA_Stream();
//				if(domains[d].Prec.USE_FLOAT) {
//					domains[d].Prec.FreeFromCUDA_Dev_fl();
//				} else {
//					domains[d].Prec.FreeFromCUDA_Dev();
//				}

#ifdef SHARE_SC
				if(domains[d].B1Kplus.USE_FLOAT) {
					// Return the original pointer from backup
					domains[d].B1Kplus.d_dense_values_fl = SC_dense_val_orig_fl[d];

					domains[d].B1Kplus.FreeFromCUDA_Dev_fl();
				} else {
					// Return the original pointer from backup
					domains[d].B1Kplus.d_dense_values = SC_dense_val_orig[d];

					domains[d].B1Kplus.FreeFromCUDA_Dev();
				}

				if(d%2 == 0 && d+1 < domains_in_global_index.size()) {
					domains[d+1].B1Kplus.is_on_acc = 0;
					domains_on_CPU++;
					domains_on_GPU--;
					DOFs_CPU += domains[d].K.rows;
					DOFs_GPU -= domains[d].K.rows;
				}
#else
				if(domains[d].B1Kplus.USE_FLOAT) {
					domains[d].B1Kplus.FreeFromCUDA_Dev_fl();
				} else {
					domains[d].B1Kplus.FreeVecsFromCUDA_Dev();
				}
#endif
				if(domains[d].B1Kplus.USE_FLOAT) {
                    cudaHostUnregister(domains[d].cuda_pinned_buff_fl.data());
					SEQ_VECTOR <float>  ().swap (domains[d].cuda_pinned_buff_fl);
				} else {
                    cudaHostUnregister(domains[d].cuda_pinned_buff.data());
					SEQ_VECTOR <double>  ().swap (domains[d].cuda_pinned_buff);
				}

				if (configuration.combine_sc_and_spds) {
					SEQ_VECTOR <double> ().swap (domains[d].B1Kplus.dense_values);
					SEQ_VECTOR <float>  ().swap (domains[d].B1Kplus.dense_values_fl);

					// if (USE_FLOAT)
					// 	std::cout << "f";
					// else
					// 	std::cout << "F";

				} else {

					// if (USE_FLOAT)
					// 	std::cout << "c";
					// else
					// 	std::cout << "C";
				}
			}
		} else {
			if (configuration.combine_sc_and_spds) {

				// if (USE_FLOAT)
				// 	std::cout << "f";
				// else
				// 	std::cout << "F";

			} else {

				// if (USE_FLOAT)
				// 	std::cout << "c";
				// else
				// 	std::cout << "C";

			}
		}

		// std::cout << " Domain: " << d << " GPU : " << domains[d].B1Kplus.is_on_acc << "\n";
	}

//	eslog::info("\n LSCs assembled on GPU: %d\n", domains_on_GPU);
//	eslog::info(" LSCs remaining on CPU: %d\n", domains_on_CPU);

	checkCudaErrors(cudaFree(d_blocked_memory));

//	cilk_for (esint i = 0; i < domains_in_global_index.size(); i++ ) {
//
////		cudaSetDevice(1);
//
////		SparseSolverCPU tmpsps2;
////		if ( i == 0 && cluster_global_index == 1) tmpsps2.msglvl = 1;
////		tmpsps2.Create_non_sym_SC_w_Mat( domains[i].K, TmpB, domains[i].B0t_comp, domains[i].B0KplusB1_comp, false, 0 );
//
//		esint status = 0;
//		cudaError_t status_c;
//
//		if (!GPU_full || !configuration.combine_sc_and_spds) {
//
//			GetSchurComplement(USE_FLOAT, i);
//
//			if (!GPU_full) {
//				if (USE_FLOAT){
//					status = domains[i].B1Kplus.CopyToCUDA_Dev_fl();
//				} else {
//					status = domains[i].B1Kplus.CopyToCUDA_Dev();
//				}
//
//				if (USE_FLOAT){
//					size_t memsize = (domains[d].B1Kplus.rows > domains[d].B1t_Dir_perm_vec.size()) ? domains[d].B1Kplus.rows : domains[d].B1t_Dir_perm_vec.size();
//
//                                      domains[d].cuda_pinned_buff_fl.resize(memsize);
//                                      cudaHostRegister(domains[d].cuda_pinned_buff_fl.data(), memsize * sizeof(float), cudaHostRegisterDefault);
//					if (domains[d].cuda_pinned_buff_fl.size() != memsize)
//					{
//						std::cout << "Error allocating pinned host memory";
//						status = -1;
//                                      }
//				} else {
//					size_t memsize = (domains[d].B1Kplus.rows > domains[d].B1t_Dir_perm_vec.size()) ? domains[d].B1Kplus.rows : domains[d].B1t_Dir_perm_vec.size();
//
//                                      domains[d].cuda_pinned_buff.resize(memsize);
//					if (domains[d].cuda_pinned_buff.size() != memsize)
//					{
//						std::cout << "Error allocating pinned host memory";
//						status = -1;
//                                      }
//                                      cudaHostRegister(domains[d].cuda_pinned_buff.data(), memsize * sizeof(double), cudaHostRegisterDefault);
//				}
//			} else {
//				status = 1;
//			}
//
//			// if status == 0 - all buffers in GPU mem were sucesfuly allocated
//			if (status == 0) {
//				domains[i].B1Kplus.is_on_acc = 1;
//				SEQ_VECTOR <double> ().swap (domains[i].B1Kplus.dense_values);
//				SEQ_VECTOR <float>  ().swap (domains[i].B1Kplus.dense_values_fl);
//				domains[i].Kplus.keep_factors = false;
//				if (USE_FLOAT)
//					////ESINFO(PROGRESS3) << Info::plain() << "g";
//				else
//					////ESINFO(PROGRESS3) << Info::plain() << "G";
//			} else {
//				domains[i].B1Kplus.is_on_acc = 0;
//				GPU_full = true;
//				if (configuration.combine_sc_and_spds) {
//					SEQ_VECTOR <double> ().swap (domains[i].B1Kplus.dense_values);
//					SEQ_VECTOR <float>  ().swap (domains[i].B1Kplus.dense_values_fl);
//					if (USE_FLOAT)
//						////ESINFO(PROGRESS3) << Info::plain() << "p";
//					else
//						////ESINFO(PROGRESS3) << Info::plain() << "P";
//				} else {
//					if (USE_FLOAT)
//						////ESINFO(PROGRESS3) << Info::plain() << "c";
//					else
//						////ESINFO(PROGRESS3) << Info::plain() << "C";
//				}
//			}
//
//		} else {
//          domains[i].B1Kplus.is_on_acc = 0;
//			if (USE_FLOAT)
//				////ESINFO(PROGRESS3) << Info::plain() << "p";
//			else
//				////ESINFO(PROGRESS3) << Info::plain() << "P";
//		}
//
//		//GPU_full = true;
//
//	}

	////ESINFO(PROGRESS3);
}


void ClusterGPU::GetSchurComplement(bool USE_FLOAT, esint i) {
	// |      K      (B1_comp_dom)t |
	// | B1_comp_dom      0         |

	SparseMatrix TmpB;
	domains[i].B1_comp_dom.MatTranspose(TmpB);
	SparseSolverCPU tmpsps;
	//	if ( i == 0 && cluster_global_index == 1) {
	//		tmpsps.msglvl = Info::report(LIBRARIES) ? 1 : 0;
	//	}

	// int order = 0;
	if (domains[i].K.type =='S') {
		tmpsps.Create_SC_w_Mat(domains[i].K, TmpB, domains[i].B1Kplus, false, 0); // general
		// order = 0; // 0 = natural, 1 = amd(A+A')
	} else if (domains[i].K.type =='G') {
<<<<<<< HEAD
		tmpsps.Create_non_sym_SC_w_Mat(domains[i].K, TmpB, TmpB, domains[i].B1Kplus, false, 0);
		// order = 0; // 0 = natural, 1 = amd(A+A'), 2 = amd(S'*S), 3 = amd(A'*A)
=======
		// 1-based indexing must be converted to 0-based indexing
		SparseMatrix B1_comp_dom_copy = domains[i].B1_comp_dom;
		
		// Convert CSR to CSC
		// 1-based indexing must be converted to 0-based indexing
		// TODO: Do not copy if possible
		SparseMatrix K_csc;
		domains[i].K.MatTranspose(K_csc);

		// CSparse factorization - 0-based indexing, CSC-format
		int order = 0; // 0 = natural, 1 = amd(AaA'), 2 = amd(S'*S), 3 = amd(A'*A)
		int tol = 1; // partial pivoting tolerance
		int n_rhs = B1_comp_dom_copy.rows;
		int device_id = domains[i].B1Kplus.device_id;

		cuda::SetDevice(device_id);
		cuda::Malloc((void**)&domains[i].B1Kplus.d_dense_values, n_rhs * n_rhs * sizeof(double));
		csparse::CreateLscGpu(K_csc, B1_comp_dom_copy, order, tol, device_id, domains[i].B1Kplus);
>>>>>>> ENH #62: Set device_id from SparseMatrix to the CSparse call
	} else {
		eslog::error("Error - not defined type of K mat type.");
		exit(0);
	}

	// DEPRECATED - The old way to compute LSC on GPU using custom CSparse - tobe removed
	// // 1-based indexing must be converted to 0-based indexing
	// SparseMatrix B1_comp_dom_copy = domains[i].B1_comp_dom;
	
	// // Convert CSR to CSC
	// // 1-based indexing must be converted to 0-based indexing
	// // TODO: Do not copy if possible
	// SparseMatrix K_csc;
	// domains[i].K.MatTranspose(K_csc);

	// // CSparse factorization - 0-based indexing, CSC-format
	// int tol = 1; // partial pivoting tolerance
	// int n_rhs = B1_comp_dom_copy.rows;
	// int device_id = domains[i].B1Kplus.device_id;
	// int print_output = 0;

	// cuda::SetDevice(device_id);
	// cuda::Malloc((void**)&domains[i].B1Kplus.d_dense_values, n_rhs * n_rhs * sizeof(double));
	// csparse::CreateLscGpu(K_csc, B1_comp_dom_copy, order, tol, device_id, print_output, domains[i].B1Kplus);

	if (USE_FLOAT){
		domains[i].B1Kplus.ConvertDenseToDenseFloat( 1 );
		domains[i].B1Kplus.USE_FLOAT = true;
	}

	//////ESINFO(PROGRESS3) << Info::plain() << "s";

	// if Schur complement is symmetric - then remove lower part - slower for GPU but more mem. efficient
	if (configuration.schur_type == FETIConfiguration::MATRIX_STORAGE::SYMMETRIC) {
		domains[i].B1Kplus.RemoveLowerDense();
	}
}


void ClusterGPU::DistributeDomains(TGPU* gpus, int n_gpu, int n_lsc) {
	// Calculate number of LSCs processed by each GPU - should be replaced with mapping based on lsc_sizes
    for (int g = 0; g < n_gpu; g++) {
        gpus[g].n_lsc_gpu = n_lsc / n_gpu;
    }
    int n_lsc_remaining = n_lsc % n_gpu;
    for (int r = 0; r < n_lsc_remaining; r++) {
        gpus[r].n_lsc_gpu++;
    }

    int n = 0;
    for (int g = 0; g < n_gpu; g++) {
        // Allocate CPU memory for IDs for each block of domains
        gpus[g].h_array_lsc_id = (int *) malloc(gpus[g].n_lsc_gpu * sizeof(int));
		gpus[g].h_array_d_lsc = (double **) malloc(gpus[g].n_lsc_gpu * sizeof(double *));

        // Assign IDs and LSC device pointers of domains to the GPU
        for (int i = 0; i < gpus[g].n_lsc_gpu; i++) {
            gpus[g].h_array_lsc_id[i] = lsc_on_gpu_ids[n];
			gpus[g].h_array_d_lsc[i] = domains[lsc_on_gpu_ids[n]].B1Kplus.d_dense_values;
            n++;
        }        
    }

    // Assign ranges
    gpus[0].start = 0;
    gpus[0].end = gpus[0].start + gpus[0].n_lsc_gpu;

    for (int g = 1; g < n_gpu; g++) {
        gpus[g].start = gpus[g - 1].end;
        gpus[g].end = gpus[g].start + gpus[g].n_lsc_gpu;
    }
}


void ClusterGPU::GetSchurComplementsGpu(bool USE_FLOAT, SEQ_VECTOR<int>& vec_L_nnz,
 SEQ_VECTOR<int*>& vec_L_row_indexes, SEQ_VECTOR<int*>& vec_L_col_pointers, SEQ_VECTOR<double*>& vec_L_values,
 SEQ_VECTOR<int>& vec_U_nnz, SEQ_VECTOR<int*>& vec_U_row_indexes, SEQ_VECTOR<int*>& vec_U_col_pointers,
 SEQ_VECTOR<double*>& vec_U_values, SEQ_VECTOR<int*>& vec_perm, SEQ_VECTOR<int*>& vec_perm_2, esint max_B1_nnz,
 esint max_B1_rows, esint max_B1_size, esint max_K_rows, esint max_L_nnz, esint max_U_nnz) {
	if(USE_FLOAT) {
		eslog::warning("ESPRESO run-time warning: Support for GPU accelerated LSC in single precision not implemented yet. Assembling LSC on CPU using Pardiso.\n");
		// Temporary CPU backup:
		#pragma omp parallel for
		for (size_t d = 0; d < domains_in_global_index.size(); d++) {
			if (domains[d].B1Kplus.is_on_acc == 1 || !configuration.combine_sc_and_spds) {
				// Calculates SC on CPU and keeps it in CPU memory
				GetSchurComplement(USE_FLOAT, d);
//				eslog::info(".");
			}
		}
	} else {
	// Calculates SC on GPU and keeps it in GPU memory

	// TODO: DEPRECATED - to be removed (also needs uncomment calling CSparse routine inside GetSchurComplement())
	// #pragma omp parallel for
	// for (esint d = 0; d < domains_in_global_index.size(); d++) {
	// 	if (domains[d].B1Kplus.is_on_acc == 1 || !configuration.combine_sc_and_spds) {
	// 		// Calculates SC on CPU and keeps it in CPU memory
	// 		GetSchurComplement(USE_FLOAT, d);
	// 		eslog::info(".");
	// 	}
	// }
	// return;

		int n_lsc = lsc_on_gpu_ids.size();
		SEQ_VECTOR<SparseMatrix> vec_B1_comp_dom_copy(n_lsc);
		int order;

		if(SYMMETRIC_SYSTEM) {
			order = 1; // 0 = natural, 1 = amd(A+A')

			#pragma omp parallel for
			for (esint d = 0; d < n_lsc; d++) {
				// Allocate device memory for LSCs
				cuda::SetDevice(device_id);
				cuda::Malloc((void**)&domains[lsc_on_gpu_ids[d]].B1Kplus.d_dense_values,
				 domains[lsc_on_gpu_ids[d]].B1_comp_dom.rows * domains[lsc_on_gpu_ids[d]].B1_comp_dom.rows * sizeof(double));
				
				domains[lsc_on_gpu_ids[d]].B1Kplus.cols = domains[lsc_on_gpu_ids[d]].B1_comp_dom.rows;
				domains[lsc_on_gpu_ids[d]].B1Kplus.rows = domains[lsc_on_gpu_ids[d]].B1_comp_dom.rows;
				domains[lsc_on_gpu_ids[d]].B1Kplus.type = 'G';

				// TODO: This copy should be eliminated by setting 1-based indexing for B in cusparse routines and using original B1
				vec_B1_comp_dom_copy[d] = domains[lsc_on_gpu_ids[d]].B1_comp_dom;
				// Convert 1-based to 0-based indexing
				for(int& i : vec_B1_comp_dom_copy[d].CSR_I_row_indices) {
					i--;
				}
				for(int& i : vec_B1_comp_dom_copy[d].CSR_J_col_indices) {
					i--;
				}
			}
		} else {
			order = 1; // 0 = natural, 1 = amd(A+A'), 2 = amd(S'*S), 3 = amd(A'*A)

			#pragma omp parallel for
			for (esint d = 0; d < n_lsc; d++) {
				// Allocate device memory for LSCs
				cuda::SetDevice(device_id);
				cuda::Malloc((void**)&domains[lsc_on_gpu_ids[d]].B1Kplus.d_dense_values,
				 domains[lsc_on_gpu_ids[d]].B1_comp_dom.rows * domains[lsc_on_gpu_ids[d]].B1_comp_dom.rows * sizeof(double));

				domains[lsc_on_gpu_ids[d]].B1Kplus.cols = domains[lsc_on_gpu_ids[d]].B1_comp_dom.rows;
				domains[lsc_on_gpu_ids[d]].B1Kplus.rows = domains[lsc_on_gpu_ids[d]].B1_comp_dom.rows;
				domains[lsc_on_gpu_ids[d]].B1Kplus.type = 'G';

				// TODO: This should be eliminated by setting 1-based indexing for B in cusparse routines
				vec_B1_comp_dom_copy[d] = domains[lsc_on_gpu_ids[d]].B1_comp_dom;
				// Convert 1-based to 0-based indexing
				for(int& i : vec_B1_comp_dom_copy[d].CSR_I_row_indices) {
					i--;
				}
				for(int& i : vec_B1_comp_dom_copy[d].CSR_J_col_indices) {
					i--;
				}
			}
		}
		POP_RANGE

		PUSH_RANGE("Prepare host arrays", 2)

		// TODO - To be removed as of #88 is implemented
		// Get max dimensions
		int max_nnz = std::max(max_L_nnz, max_U_nnz);

		// General settings
		// Currently only 1 GPU with specified device_id per cluster
		int n_gpu = 1;
		int gpu_id = device_id;
		// Array of gpu structs
        TGPU *gpus = (TGPU*) malloc(n_gpu * sizeof(TGPU));

		// Distribute LSCs among GPUs in case of multi-GPU per cluster (sets gpus[g].n_lsc_gpu for the next step)
		// Assign LSC device pointers
		DistributeDomains(gpus, n_gpu, n_lsc);

		// Allocate host arrays
		// 1) For each GPU
		for (int g = 0; g < n_gpu; g++) {
            gpus[g].h_array_stream = (cudaStream_t *) malloc(n_streams_per_gpu * sizeof(cudaStream_t *));
            gpus[g].h_array_handle = (cusparseHandle_t *) malloc(n_streams_per_gpu * sizeof(cusparseHandle_t *));
            gpus[g].h_array_info_L = (csrsm2Info_t *) malloc(n_csrsm2_info_per_gpu * sizeof(csrsm2Info_t *));
            gpus[g].h_array_info_U = (csrsm2Info_t *) malloc(n_csrsm2_info_per_gpu * sizeof(csrsm2Info_t *));

            gpus[g].h_array_d_csr_B_row_ptr = (int **) malloc(n_streams_per_gpu * sizeof(int *));
            gpus[g].h_array_d_csr_B_col_ind = (int **) malloc(n_streams_per_gpu * sizeof(int *));
            gpus[g].h_array_d_csr_B_val = (double **) malloc(n_streams_per_gpu * sizeof(double *));
            gpus[g].h_array_d_buffer = (char **) malloc(n_streams_per_gpu * sizeof(char *));

            gpus[g].h_array_d_csr_Bt_row_ptr = (int **) malloc(n_streams_per_gpu * sizeof(int *));
            gpus[g].h_array_d_csr_Bt_col_ind = (int **) malloc(n_streams_per_gpu * sizeof(int *));
            gpus[g].h_array_d_csr_Bt_val = (double **) malloc(n_streams_per_gpu * sizeof(double *));

            gpus[g].h_array_d_X_reordered = (double **) malloc(n_streams_per_gpu * sizeof(double *));

			gpus[g].h_array_d_pinv = (int **) malloc(n_streams_per_gpu * sizeof(int *));
			gpus[g].h_array_d_q = (int **) malloc(n_streams_per_gpu * sizeof(int *));
        }
		// 2) cuSparse variables unique for each LSC
        // TODO check if one descriptor for all LSC is enough
        cusparseSpMatDescr_t *h_array_h_matB = (cusparseSpMatDescr_t *) malloc(n_lsc * sizeof(cusparseSpMatDescr_t *));
        cusparseDnMatDescr_t *h_array_h_matX = (cusparseDnMatDescr_t *) malloc(n_lsc * sizeof(cusparseDnMatDescr_t *));
        cusparseDnMatDescr_t *h_array_h_matLSC = (cusparseDnMatDescr_t *) malloc(n_lsc * sizeof(cusparseDnMatDescr_t *));
		POP_RANGE  // END Prepare host arrays

		PUSH_RANGE("Prepare device arrays", 3)
		// cuSparse CSRSM2 settings
		cusparseMatDescr_t h_descrL = NULL;
        cusparseMatDescr_t h_descrU = NULL;
        cusparseMatDescr_t h_descrBt = NULL;
        int solve_algo = 1;  // solve_algo = 0 is non-block version; solve_algo = 1 is block version
        cusparseOperation_t solve_transL = CUSPARSE_OPERATION_TRANSPOSE; // Transposed matrix require larger info object
        cusparseOperation_t solve_transU = SYMMETRIC_SYSTEM ? CUSPARSE_OPERATION_NON_TRANSPOSE : CUSPARSE_OPERATION_TRANSPOSE;  // L^T and U in CSC
        cusparseOperation_t solve_transB = CUSPARSE_OPERATION_NON_TRANSPOSE; // TODO LU probably change for U
        double solve_alpha = 1.0;
        cusparseSolvePolicy_t policy = CUSPARSE_SOLVE_POLICY_NO_LEVEL; // CUSPARSE_SOLVE_POLICY_USE_LEVEL
        // TODO: set according the largest domain #54
        size_t expected_buffer_size = ((max_K_rows + max_nnz) * 4 + max_B1_size) * sizeof(double);
		// This should be correct according NVIDIA but is too low in real
		// size_t expected_buffer_size = (max_K_rows + max_nnz) * sizeof(float) + max_B1_size * sizeof(double); 
        size_t solve_buffer_size = expected_buffer_size; // solve_buffer_size is overwritten during the computation
		// Configuration of matrix descriptor L, U, and B^T*/
        checkCudaErrors(cusparseCreateMatDescr(&h_descrL));
        /* Loaded A is base-1 but CSparse converts it to base-0 */
        checkCudaErrors(cusparseSetMatIndexBase(h_descrL, CUSPARSE_INDEX_BASE_ZERO));
        // A is actually symmetric (stored as lower triangular part only),
        // but csrsm2 routine supports only CUSPARSE_MATRIX_TYPE_GENERAL
        checkCudaErrors(cusparseSetMatType(h_descrL, CUSPARSE_MATRIX_TYPE_GENERAL));
        /* A is lower triangle in CSR but upper triangle in CSC*/
        checkCudaErrors(cusparseSetMatFillMode(h_descrL, CUSPARSE_FILL_MODE_UPPER));
        /* A has non unit diagonal */
        checkCudaErrors(cusparseSetMatDiagType(h_descrL, SYMMETRIC_SYSTEM ? CUSPARSE_DIAG_TYPE_NON_UNIT : CUSPARSE_DIAG_TYPE_UNIT));

        checkCudaErrors(cusparseCreateMatDescr(&h_descrU));
        /* Loaded A is base-1 but CSparse converts it to base-0 */
        checkCudaErrors(cusparseSetMatIndexBase(h_descrU, CUSPARSE_INDEX_BASE_ZERO));
        // A is actually symmetric (stored as lower triangular part only),
        // but csrsm2 routine supports only CUSPARSE_MATRIX_TYPE_GENERAL
        checkCudaErrors(cusparseSetMatType(h_descrU, CUSPARSE_MATRIX_TYPE_GENERAL));
        /* A is lower triangle in CSR but upper triangle in CSC*/
        checkCudaErrors(cusparseSetMatFillMode(h_descrU, SYMMETRIC_SYSTEM ? CUSPARSE_FILL_MODE_UPPER : CUSPARSE_FILL_MODE_LOWER));
        /* A has non unit diagonal */
        checkCudaErrors(cusparseSetMatDiagType(h_descrU, CUSPARSE_DIAG_TYPE_NON_UNIT));

        checkCudaErrors(cusparseCreateMatDescr(&h_descrBt));
        /* Loaded B is base-1 but CSparse converts it to base-0 */
        checkCudaErrors(cusparseSetMatIndexBase(h_descrBt, CUSPARSE_INDEX_BASE_ZERO));
        // A is actually symmetric (stored as lower triangular part only),
        // but csrsm2 routine supports only CUSPARSE_MATRIX_TYPE_GENERAL
        checkCudaErrors(cusparseSetMatType(h_descrBt, CUSPARSE_MATRIX_TYPE_GENERAL));
        /* A has non unit diagonal */
        checkCudaErrors(cusparseSetMatDiagType(h_descrBt, CUSPARSE_DIAG_TYPE_NON_UNIT));
		
		// cuSparse SpMM settings
		cusparseOperation_t spmm_transB = CUSPARSE_OPERATION_NON_TRANSPOSE;
        cusparseOperation_t spmm_transX = CUSPARSE_OPERATION_NON_TRANSPOSE;
        double spmm_alpha = 1.0;
        double spmm_beta = 0.0;
        cusparseSpMMAlg_t spmm_algo = CUSPARSE_CSRMM_ALG1; // CUSPARSE_MM_ALG_DEFAULT
        cudaDataType compute_type = cudaDataType::CUDA_R_64F;
        size_t spmm_buffer_size = 0;
        cusparseOrder_t spmm_order = CUSPARSE_ORDER_COL;
        cusparseIndexType_t spmm_csr_row_offsets_type = CUSPARSE_INDEX_32I;
        cusparseIndexType_t spmm_csr_col_ind_type = CUSPARSE_INDEX_32I;
        cusparseIndexBase_t spmm_idx_base = CUSPARSE_INDEX_BASE_ZERO;
		
		// Allocate device arrays and initialize CUDA streams on all assigned GPUs
		for (int g = 0; g < n_gpu; g++) {
            checkCudaErrors(cudaSetDevice(g + gpu_id));

            checkCudaErrors(cudaStreamCreateWithFlags(&gpus[g].data_stream, cudaStreamNonBlocking));
          
            for (int s = 0; s < n_streams_per_gpu; s++) {
                checkCudaErrors(cudaStreamCreateWithFlags(&gpus[g].h_array_stream[s], cudaStreamNonBlocking));
                checkCudaErrors(cusparseCreate(&gpus[g].h_array_handle[s]));
                checkCudaErrors(cusparseSetStream(gpus[g].h_array_handle[s], gpus[g].h_array_stream[s]));

                // Allocate device arrays shared among LSCs
                // Allocate device memory for CSR B, CSR B^T and dense B^T (aka X) matrices
                checkCudaErrors(cudaMalloc((void **)&gpus[g].h_array_d_csr_B_col_ind[s], max_B1_nnz * sizeof(int)));
                checkCudaErrors(cudaMalloc((void **)&gpus[g].h_array_d_csr_B_row_ptr[s], (max_B1_rows + 1) * sizeof(int)));
                checkCudaErrors(cudaMalloc((void **)&gpus[g].h_array_d_csr_B_val[s], max_B1_nnz * sizeof(double)));

                // Comment out in order to eliminate Bt_csr
                checkCudaErrors(cudaMalloc((void **)&gpus[g].h_array_d_csr_Bt_col_ind[s], max_B1_nnz * sizeof(int)));
                checkCudaErrors(cudaMalloc((void **)&gpus[g].h_array_d_csr_Bt_row_ptr[s], (max_K_rows + 1) * sizeof(int)));
                checkCudaErrors(cudaMalloc((void **)&gpus[g].h_array_d_csr_Bt_val[s], max_B1_nnz * sizeof(double)));

                checkCudaErrors(cudaMalloc((void **)&gpus[g].h_array_d_X_reordered[s], max_B1_size * sizeof(double)));
                checkCudaErrors(cudaMalloc((void **)&gpus[g].h_array_d_buffer[s], solve_buffer_size));

				if(SYMMETRIC_SYSTEM) {
					if(order) {
						checkCudaErrors(cudaMalloc((void **)&gpus[g].h_array_d_pinv[s], max_K_rows * sizeof(int)));
					} else {
						gpus[g].h_array_d_pinv[s] = NULL;
					}
					gpus[g].h_array_d_q[s] = NULL;
				} else {
					if(order) {
						checkCudaErrors(cudaMalloc((void **)&gpus[g].h_array_d_q[s], max_K_rows * sizeof(int)));
					} else {
						gpus[g].h_array_d_q[s] = NULL;
					}
					checkCudaErrors(cudaMalloc((void **)&gpus[g].h_array_d_pinv[s], max_K_rows * sizeof(int)));
				}
            }
            checkCudaErrors(cudaEventCreateWithFlags(&gpus[g].event_data_preload, cudaEventDisableTiming));
            checkCudaErrors(cudaEventCreateWithFlags(&gpus[g].event1, cudaEventDisableTiming));
            checkCudaErrors(cudaEventCreateWithFlags(&gpus[g].event2, cudaEventDisableTiming));

            for (int i = 0; i < n_csrsm2_info_per_gpu; i++) {
                checkCudaErrors(cusparseCreateCsrsm2Info(&gpus[g].h_array_info_L[i]));
                checkCudaErrors(cusparseCreateCsrsm2Info(&gpus[g].h_array_info_U[i]));
            }

            // Allocate only one set of device arrays for CSC L, U (L^T), and permutation vectors per GPU
            checkCudaErrors(cudaMalloc((void **)&gpus[g].d_csc_L_row_ind, max_L_nnz * sizeof(int)));
            checkCudaErrors(cudaMalloc((void **)&gpus[g].d_csc_L_col_ptr, (max_K_rows + 1) * sizeof(int)));
            checkCudaErrors(cudaMalloc((void **)&gpus[g].d_csc_L_val, max_L_nnz * sizeof(double)));

            checkCudaErrors(cudaMalloc((void **)&gpus[g].d_csc_U_row_ind, max_U_nnz * sizeof(int)));
            checkCudaErrors(cudaMalloc((void **)&gpus[g].d_csc_U_col_ptr, (max_K_rows + 1) * sizeof(int)));
            checkCudaErrors(cudaMalloc((void **)&gpus[g].d_csc_U_val, max_U_nnz * sizeof(double)));
        }
		POP_RANGE  // END Prepare device arrays

		// Assembly Schur complements
		for (int g = 0; g < n_gpu; g++) {
            checkCudaErrors(cudaSetDevice(g + gpu_id));
			int idx = gpus[g].h_array_lsc_id[0];

            // Copy CSR B and CSC L for the first domain from host to device
            PUSH_RANGE("B_csr 0 memcpy HtD", 4)
            checkCudaErrors(cudaMemcpyAsync(gpus[g].h_array_d_csr_B_val[0], vec_B1_comp_dom_copy[idx].CSR_V_values.data(),
			 domains[idx].B1_comp_dom.CSR_V_values.size() * sizeof(double), cudaMemcpyHostToDevice, gpus[g].data_stream));
			checkCudaErrors(cudaMemcpyAsync(gpus[g].h_array_d_csr_B_col_ind[0], vec_B1_comp_dom_copy[idx].CSR_J_col_indices.data(),
			 domains[idx].B1_comp_dom.CSR_V_values.size() * sizeof(int), cudaMemcpyHostToDevice, gpus[g].data_stream));
			checkCudaErrors(cudaMemcpyAsync(gpus[g].h_array_d_csr_B_row_ptr[0], vec_B1_comp_dom_copy[idx].CSR_I_row_indices.data(),
			 (domains[idx].B1_comp_dom.rows + 1) * sizeof(int), cudaMemcpyHostToDevice, gpus[g].data_stream));
            POP_RANGE

            PUSH_RANGE("L_csr 0 memcpy HtD", 4)
            checkCudaErrors(cudaMemcpyAsync(gpus[g].d_csc_L_val, vec_L_values[idx], vec_L_nnz[idx] * sizeof(double),
			 cudaMemcpyHostToDevice, gpus[g].data_stream));
            checkCudaErrors(cudaMemcpyAsync(gpus[g].d_csc_L_row_ind, vec_L_row_indexes[idx], vec_L_nnz[idx] * sizeof(int),
			 cudaMemcpyHostToDevice, gpus[g].data_stream));
            checkCudaErrors(cudaMemcpyAsync(gpus[g].d_csc_L_col_ptr, vec_L_col_pointers[idx], (domains[idx].K.rows + 1) * sizeof(int),
			 cudaMemcpyHostToDevice, gpus[g].data_stream));
            POP_RANGE

			PUSH_RANGE("perm 0 memcpy HtD", 4)
			if(SYMMETRIC_SYSTEM) {
				if(order) {
					checkCudaErrors(cudaMemcpyAsync(gpus[g].h_array_d_pinv[0], vec_perm[idx], domains[idx].K.rows * sizeof(int),
					 cudaMemcpyHostToDevice, gpus[g].data_stream));
				}
            } else {
				if(order) {
					checkCudaErrors(cudaMemcpyAsync(gpus[g].h_array_d_q[0], vec_perm_2[idx], domains[idx].K.rows * sizeof(int),
					 cudaMemcpyHostToDevice, gpus[g].data_stream));
				}
				checkCudaErrors(cudaMemcpyAsync(gpus[g].h_array_d_pinv[0], vec_perm[idx], domains[idx].K.rows * sizeof(int),
				 cudaMemcpyHostToDevice, gpus[g].data_stream));
            }
			POP_RANGE

            // Record all the remaining work of the data stream into the event_data_preload
            checkCudaErrors(cudaEventRecord(gpus[g].event_data_preload, gpus[g].data_stream));
        }

		// Prepared for multi-GPU per cluster
		//#pragma omp parallel private(tid) num_threads(n_gpu)
        {
            // OMP thread IDs are mapped on GPU IDs
            // int g = omp_get_thread_num();
            // int tid = omp_get_thread_num();
            int g = 0;
            checkCudaErrors(cudaSetDevice(g + gpu_id));
            
            // Iterations in a chunk must be consequent because of preloading data (L, B)
            for (int i = gpus[g].start; i < gpus[g].end; i++) {
                char range_name[32];
                snprintf(range_name, sizeof range_name, "LSC %d", i);
                PUSH_RANGE(range_name, 0) // LSC

                int i_gpu = i - gpus[g].start;
                int s_gpu = i_gpu % n_streams_per_gpu;
                int info_gpu = i_gpu % n_csrsm2_info_per_gpu;

                PUSH_RANGE("Create SpMM objects", 5)
                // Create B, X, LSC cusparse matrix objects
                checkCudaErrors(cusparseCreateCsr(&h_array_h_matB[i], (int64_t) domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.rows,
				 (int64_t) domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows, (int64_t) domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.CSR_V_values.size(), (void*) gpus[g].h_array_d_csr_B_row_ptr[s_gpu],
                (void*) gpus[g].h_array_d_csr_B_col_ind[s_gpu], (void*) gpus[g].h_array_d_csr_B_val[s_gpu],
                spmm_csr_row_offsets_type, spmm_csr_col_ind_type, spmm_idx_base, compute_type));
                // h_array_d_buffer used instead of h_array_d_X
                checkCudaErrors(cusparseCreateDnMat(&h_array_h_matX[i], (int64_t) domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows,
				 (int64_t) domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.rows, (int64_t) domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows,
                (void*) gpus[g].h_array_d_buffer[s_gpu], compute_type, spmm_order));
                checkCudaErrors(cusparseCreateDnMat(&h_array_h_matLSC[i], (int64_t) domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.rows,
				 (int64_t) domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.rows, (int64_t) domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.rows,
                (void*) gpus[g].h_array_d_lsc[i_gpu], compute_type, spmm_order));
                POP_RANGE // Create SpMM objects

                // Need to wait for finishing B_csr and L_csr data transfers from data preload (before loop)
                checkCudaErrors(cudaStreamWaitEvent(gpus[g].h_array_stream[0], gpus[g].event_data_preload, 0));
                // Need to wait for finishing L-Solve from the previous iteration
                checkCudaErrors(cudaStreamWaitEvent(gpus[g].h_array_stream[s_gpu], gpus[g].event1, 0));
                // Need to wait for finishing B_csr and L_csr data transfers started in the previous iteration - very unlikely - only for the case with data transfer longer than L-Solve
                checkCudaErrors(cudaStreamWaitEvent(gpus[g].h_array_stream[s_gpu], gpus[g].event2, 0));
                PUSH_RANGE("Transpose B_csr", 6)
                // Transpose B to B^T
                // CUSPARSE_CSR2CSC_ALG2 - faster, the same ordering not garanteed
                // CUSPARSE_ACTION_SYMBOLIC most probably only for symmetric matrices

                // Currently not necessary to get buffer size
// #ifdef DEBUG
//                 checkCudaErrors(cusparseCsr2cscEx2_bufferSize(gpus[g].h_array_handle[s_gpu], domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.rows, domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows, domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.CSR_V_values.size(), gpus[g].h_array_d_csr_B_val[s_gpu],
//                  gpus[g].h_array_d_csr_B_row_ptr[s_gpu], gpus[g].h_array_d_csr_B_col_ind[s_gpu], gpus[g].h_array_d_csr_Bt_val[s_gpu], gpus[g].h_array_d_csr_Bt_row_ptr[s_gpu],
//                   gpus[g].h_array_d_csr_Bt_col_ind[s_gpu], compute_type, CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO, CUSPARSE_CSR2CSC_ALG1,
//                    &solve_buffer_size));
//                 printf("Thr %d GPU %d: LSC %d: Only %d (# CUDA streams) shared %f MB buffers for all domains and both solves\n", tid, g, i, n_streams_per_gpu, (double)solve_buffer_size / 1024.0 / 1024.0);
//                 if(solve_buffer_size > CS_MAX(1024 * 1024, expected_buffer_size))
//                     printf("CSparse LSC warning: Thr %d GPU %d: LSC %d: The calculated Csr2Csc buffer size (%f MB) is larger than the size of the used buffer (%f MB).\n",
//                      tid, g, i, (double) solve_buffer_size / 1024.0 / 1024.0, (double) CS_MAX(1024 * 1024, expected_buffer_size) / 1024.0 / 1024.0);
//                 checkCudaErrors(cudaDeviceSynchronize());
// #endif

                checkCudaErrors(cusparseCsr2cscEx2(gpus[g].h_array_handle[s_gpu], domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.rows,
				 domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows, domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.CSR_V_values.size(), gpus[g].h_array_d_csr_B_val[s_gpu],
                gpus[g].h_array_d_csr_B_row_ptr[s_gpu], gpus[g].h_array_d_csr_B_col_ind[s_gpu], gpus[g].h_array_d_csr_Bt_val[s_gpu], 
                gpus[g].h_array_d_csr_Bt_row_ptr[s_gpu], gpus[g].h_array_d_csr_Bt_col_ind[s_gpu], compute_type, CUSPARSE_ACTION_NUMERIC,
                CUSPARSE_INDEX_BASE_ZERO, CUSPARSE_CSR2CSC_ALG1, gpus[g].h_array_d_buffer[s_gpu]));
                POP_RANGE // Transpose B_csr

                PUSH_RANGE("Convert Bt_csr to dense", 0)
                // Convert B^T from CSR to dense on device
                // m = Bt_cs->m = domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows
                // n = Bt_cs->n = domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.rows
                // h_array_d_buffer used instead of h_array_d_X
                checkCudaErrors(cusparseDcsr2dense(gpus[g].h_array_handle[s_gpu], domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows,
				 domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.rows, h_descrBt, gpus[g].h_array_d_csr_Bt_val[s_gpu],
                gpus[g].h_array_d_csr_Bt_row_ptr[s_gpu], gpus[g].h_array_d_csr_Bt_col_ind[s_gpu], (double *) gpus[g].h_array_d_buffer[s_gpu],
				 domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows));
                POP_RANGE //Convert Bt_csr to dense

                PUSH_RANGE("Input reordering", 1)
                // Perform reordering of dense values by CS on device
                // h_array_d_buffer used instead of h_array_d_X
                // Symm: h_array_d_X_reordered(S->pinv) = h_array_d_buffer
                // Unsym: h_array_d_X_reordered(N->pinv) = h_array_d_buffer
                cuda::IpvecReorderMrhs(gpus[g].h_array_d_pinv[s_gpu], (double *) gpus[g].h_array_d_buffer[s_gpu], gpus[g].h_array_d_X_reordered[s_gpu],
                domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows, domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.rows, gpus[g].h_array_stream[s_gpu]); /* b = P'*x */
                POP_RANGE // Input reordering

                PUSH_RANGE("Querry workspace lsolve", 2)
                /* step 3: query workspace */
                // requires sorted CSR
                checkCudaErrors(cusparseDcsrsm2_bufferSizeExt(gpus[g].h_array_handle[s_gpu], solve_algo, solve_transL, solve_transB,
				 domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows, domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.rows, vec_L_nnz[i],
                 &solve_alpha, h_descrL, gpus[g].d_csc_L_val, gpus[g].d_csc_L_col_ptr, gpus[g].d_csc_L_row_ind, gpus[g].h_array_d_X_reordered[s_gpu],
				  domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows, gpus[g].h_array_info_L[info_gpu], policy, &solve_buffer_size));
                
				if(solve_buffer_size > expected_buffer_size)
                    printf("CSparse LSC warning: GPU %d: LSC %d: The calculated solve buffer size (%f MB) is larger than the size of the used buffer (%f MB).\n",
                     g, i, (double) solve_buffer_size / 1024.0 / 1024.0, (double) expected_buffer_size / 1024.0 / 1024.0);
                POP_RANGE // Querry workspace lsolve
                
                /* step 4: analysis */
                PUSH_RANGE("Analysis lsolve", 3)
                checkCudaErrors(cusparseDcsrsm2_analysis(gpus[g].h_array_handle[s_gpu], solve_algo, solve_transL, solve_transB, domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows,
				 domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.rows, vec_L_nnz[i], &solve_alpha, h_descrL, gpus[g].d_csc_L_val, gpus[g].d_csc_L_col_ptr, gpus[g].d_csc_L_row_ind,
                 gpus[g].h_array_d_X_reordered[s_gpu], domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows, gpus[g].h_array_info_L[info_gpu], policy, gpus[g].h_array_d_buffer[s_gpu]));
                POP_RANGE // Analysis lsolve

                // Overlap CSC U (L^T) transfer with L-Solve
                PUSH_RANGE("U_csc memcpy HtD", 4)
                checkCudaErrors(cudaMemcpyAsync(gpus[g].d_csc_U_val, SYMMETRIC_SYSTEM ? vec_L_values[i] : vec_U_values[i],
				 SYMMETRIC_SYSTEM ? vec_L_nnz[i] : vec_U_nnz[i] * sizeof(double), cudaMemcpyHostToDevice, gpus[g].data_stream));
                checkCudaErrors(cudaMemcpyAsync(gpus[g].d_csc_U_row_ind, SYMMETRIC_SYSTEM ? vec_L_row_indexes[i] : vec_U_row_indexes[i],
				 SYMMETRIC_SYSTEM ? vec_L_nnz[i] : vec_U_nnz[i] * sizeof(int), cudaMemcpyHostToDevice, gpus[g].data_stream));
                checkCudaErrors(cudaMemcpyAsync(gpus[g].d_csc_U_col_ptr, SYMMETRIC_SYSTEM ? vec_L_col_pointers[i] : vec_U_col_pointers[i],
				 (domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows + 1) * sizeof(int), cudaMemcpyHostToDevice, gpus[g].data_stream));
                POP_RANGE // U_csc memcpy HtD

                /* step 5: solve X = L * X */
                PUSH_RANGE("Lsolve", 5)
                checkCudaErrors(cusparseDcsrsm2_solve(gpus[g].h_array_handle[s_gpu], solve_algo, solve_transL, solve_transB,
				 domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows, domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.rows, vec_L_nnz[i],
				 &solve_alpha, h_descrL, gpus[g].d_csc_L_val, gpus[g].d_csc_L_col_ptr, gpus[g].d_csc_L_row_ind, gpus[g].h_array_d_X_reordered[s_gpu],
                 domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows, gpus[g].h_array_info_L[info_gpu], policy, gpus[g].h_array_d_buffer[s_gpu]));
                POP_RANGE // Lsolve

                // Record all the remaining work of the compute stream into the event1
                checkCudaErrors(cudaEventRecord(gpus[g].event1, gpus[g].h_array_stream[s_gpu]));
                // Data stream waits until all current work in the compute stream is done
                checkCudaErrors(cudaStreamWaitEvent(gpus[g].data_stream, gpus[g].event1, 0));

                // Overlap CSR B and CSC L transfer for the following iteration with the L^T-Solve (U-Solve)
                if((i_gpu + 1) < gpus[g].n_lsc_gpu) {
                    PUSH_RANGE("B_csr+n_gpu memcpy HtD", 6)
                    int s_gpu_following = (i_gpu + 1) % n_streams_per_gpu;
					checkCudaErrors(cudaMemcpyAsync(gpus[g].h_array_d_csr_B_val[s_gpu_following], vec_B1_comp_dom_copy[i + 1].CSR_V_values.data(),
					 domains[gpus[g].h_array_lsc_id[i_gpu + 1]].B1_comp_dom.CSR_V_values.size() * sizeof(double), cudaMemcpyHostToDevice, gpus[g].data_stream));
                    checkCudaErrors(cudaMemcpyAsync(gpus[g].h_array_d_csr_B_col_ind[s_gpu_following], vec_B1_comp_dom_copy[i + 1].CSR_J_col_indices.data(),
					 domains[gpus[g].h_array_lsc_id[i_gpu + 1]].B1_comp_dom.CSR_V_values.size() * sizeof(int), cudaMemcpyHostToDevice, gpus[g].data_stream));
                    checkCudaErrors(cudaMemcpyAsync(gpus[g].h_array_d_csr_B_row_ptr[s_gpu_following], vec_B1_comp_dom_copy[i + 1].CSR_I_row_indices.data(),
					 (domains[gpus[g].h_array_lsc_id[i_gpu + 1]].B1_comp_dom.rows + 1) * sizeof(int), cudaMemcpyHostToDevice, gpus[g].data_stream));
                    POP_RANGE // B_csr+n_gpu memcpy HtD

                    PUSH_RANGE("L_csc+n_gpu memcpy HtD", 0)
                    checkCudaErrors(cudaMemcpyAsync(gpus[g].d_csc_L_val, vec_L_values[i + 1], vec_L_nnz[i + 1] * sizeof(double),
					 cudaMemcpyHostToDevice, gpus[g].data_stream));
                    checkCudaErrors(cudaMemcpyAsync(gpus[g].d_csc_L_row_ind, vec_L_row_indexes[i + 1], vec_L_nnz[i + 1] * sizeof(int),
					 cudaMemcpyHostToDevice, gpus[g].data_stream));
                    checkCudaErrors(cudaMemcpyAsync(gpus[g].d_csc_L_col_ptr, vec_L_col_pointers[i + 1], (domains[gpus[g].h_array_lsc_id[i_gpu + 1]].K.rows + 1) * sizeof(int),
					 cudaMemcpyHostToDevice, gpus[g].data_stream));
                    POP_RANGE // L_csc+n_gpu memcpy HtD

					if(SYMMETRIC_SYSTEM) {
						if(order) {
							checkCudaErrors(cudaMemcpyAsync(gpus[g].h_array_d_pinv[s_gpu_following], vec_perm[i + 1],
							domains[gpus[g].h_array_lsc_id[i_gpu + 1]].K.rows * sizeof(int), cudaMemcpyHostToDevice, gpus[g].data_stream));
						}
					} else {
						if(order) {
							checkCudaErrors(cudaMemcpyAsync(gpus[g].h_array_d_q[s_gpu_following], vec_perm_2[i + 1],
							domains[gpus[g].h_array_lsc_id[i_gpu + 1]].K.rows * sizeof(int), cudaMemcpyHostToDevice, gpus[g].data_stream));
						}
						checkCudaErrors(cudaMemcpyAsync(gpus[g].h_array_d_pinv[s_gpu_following], vec_perm[i + 1],
						 domains[gpus[g].h_array_lsc_id[i_gpu + 1]].K.rows * sizeof(int), cudaMemcpyHostToDevice, gpus[g].data_stream));
					}

                    // Record all the remaining work of the data stream into the event2
                    checkCudaErrors(cudaEventRecord(gpus[g].event2, gpus[g].data_stream));
                }

                /* step 6: query workspace */
                PUSH_RANGE("Querry workspace ltsolve", 1)
                // TODO LU check correctness in case of problems
                // requires sorted CSR
                checkCudaErrors(cusparseDcsrsm2_bufferSizeExt(gpus[g].h_array_handle[s_gpu], solve_algo, solve_transU, solve_transB,
				 domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows, domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.rows, SYMMETRIC_SYSTEM ? vec_L_nnz[i] : vec_U_nnz[i],
				 &solve_alpha, h_descrU, gpus[g].d_csc_U_val, gpus[g].d_csc_U_col_ptr, gpus[g].d_csc_U_row_ind, gpus[g].h_array_d_X_reordered[s_gpu],
				 domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows, gpus[g].h_array_info_U[info_gpu], policy, &solve_buffer_size));

                if(solve_buffer_size > expected_buffer_size)
                    printf("CSparse LSC warning: GPU %d: LSC %d: The calculated solve buffer size (%f MB) is larger than the size of the used buffer (%f MB).\n",
                     g, i, (double) solve_buffer_size / 1024.0 / 1024.0, (double) expected_buffer_size / 1024.0 / 1024.0);
                POP_RANGE // Querry workspace ltsolve

                /* step 7: analysis */
                PUSH_RANGE("Analysis ltsolve", 2)
                checkCudaErrors(cusparseDcsrsm2_analysis(gpus[g].h_array_handle[s_gpu], solve_algo, solve_transU, solve_transB,
				 domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows, domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.rows, SYMMETRIC_SYSTEM ? vec_L_nnz[i] : vec_U_nnz[i],
                 &solve_alpha, h_descrU, gpus[g].d_csc_U_val, gpus[g].d_csc_U_col_ptr, gpus[g].d_csc_U_row_ind, gpus[g].h_array_d_X_reordered[s_gpu],
				 domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows, gpus[g].h_array_info_U[info_gpu], policy, gpus[g].h_array_d_buffer[s_gpu]));
                POP_RANGE // Analysis ltsolve

                /* step 8: solve X = Lt\X */
                PUSH_RANGE("Ltsolve", 3)
                checkCudaErrors(cusparseDcsrsm2_solve(gpus[g].h_array_handle[s_gpu], solve_algo, solve_transU, solve_transB,
				 domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows, domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.rows, SYMMETRIC_SYSTEM ? vec_L_nnz[i] : vec_U_nnz[i],
				 &solve_alpha, h_descrU, gpus[g].d_csc_U_val, gpus[g].d_csc_U_col_ptr, gpus[g].d_csc_U_row_ind, gpus[g].h_array_d_X_reordered[s_gpu],
				domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows, gpus[g].h_array_info_U[info_gpu], policy, gpus[g].h_array_d_buffer[s_gpu]));
                POP_RANGE // Ltsolve

                PUSH_RANGE("Output reordering", 4)
                // Perform CS backward reordering on device
                // h_array_d_buffer used instead of h_array_d_X
                // Symm: h_array_d_buffer = h_array_d_X_reordered(S->pinv)
                // Unsym: h_array_d_buffer(S->q) = h_array_d_X_reordered
                if (SYMMETRIC_SYSTEM) {
                    cuda::PvecReorderMrhs(gpus[g].h_array_d_pinv[s_gpu], gpus[g].h_array_d_X_reordered[s_gpu], (double *) gpus[g].h_array_d_buffer[s_gpu],
					domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows, domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.rows, gpus[g].h_array_stream[s_gpu]); /* b = P'*x */
                } else {
                    cuda::IpvecReorderMrhs(gpus[g].h_array_d_q[s_gpu], gpus[g].h_array_d_X_reordered[s_gpu], (double *) gpus[g].h_array_d_buffer[s_gpu],
					domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows, domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.rows, gpus[g].h_array_stream[s_gpu]);
                }
                POP_RANGE // Output reordering
// #ifdef DEBUG
//                 // Currently not necessary to get the buffer size (returns 0 for CSR algo)
//                 checkCudaErrors(cusparseSpMM_bufferSize(gpus[g].h_array_handle[s_gpu], spmm_transB, spmm_transX, &spmm_alpha, h_array_h_matB[i],
//                 h_array_h_matX[i], &spmm_beta, h_array_h_matLSC[i], compute_type, spmm_algo, &spmm_buffer_size));
    
//                 printf("Thr %d GPU %d: LSC %d: Calculated SpMM spmm_buffer_size = %f MB\n", tid, g, i, (double)spmm_buffer_size / 1024.0 / 1024.0);
//                 if(spmm_buffer_size > (domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows * domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.rows * sizeof(double)))
//                     printf("CSparse LSC warning: Thr %d GPU %d: LSC %d: The calculated SpMM buffer size (%f MB) is larger than the size of the used buffer (%f MB).\n",
//                      tid, g, i, (double) spmm_buffer_size / 1024.0 / 1024.0, (double) (domains[gpus[g].h_array_lsc_id[i_gpu]].K.rows * domains[gpus[g].h_array_lsc_id[i_gpu]].B1_comp_dom.rows * sizeof(double)) / 1024.0 / 1024.0);
//                 checkCudaErrors(cudaDeviceSynchronize());
// #endif
                PUSH_RANGE("SpMM", 5)
                // h_array_d_buffer used instead of h_array_d_X
                // h_array_d_X_reordered used instead of h_array_d_buffer
                // LSC = B_cs * X
                checkCudaErrors(cusparseSpMM(gpus[g].h_array_handle[s_gpu], spmm_transB, spmm_transX, &spmm_alpha, h_array_h_matB[i],
                h_array_h_matX[i], &spmm_beta, h_array_h_matLSC[i], compute_type, spmm_algo, gpus[g].h_array_d_X_reordered[s_gpu]));
                POP_RANGE // SpMM

                PUSH_RANGE("Destroy CSRSM2 and SpMMobjects", 6);
                checkCudaErrors(cusparseDestroySpMat(h_array_h_matB[i]));
                checkCudaErrors(cusparseDestroyDnMat(h_array_h_matX[i]));
                checkCudaErrors(cusparseDestroyDnMat(h_array_h_matLSC[i]));                
                POP_RANGE // Destroy CSRSM2 and SpMMobjects

				// TODO: Check - The first iteration probably should not do this
                if(i_gpu % n_csrsm2_info_per_gpu == 0) {
                    for (int j = 0; j < n_csrsm2_info_per_gpu; j++) {
                        checkCudaErrors(cusparseDestroyCsrsm2Info(gpus[g].h_array_info_L[j]));
                        checkCudaErrors(cusparseDestroyCsrsm2Info(gpus[g].h_array_info_U[j]));
                        checkCudaErrors(cusparseCreateCsrsm2Info(&gpus[g].h_array_info_L[j]));
                        checkCudaErrors(cusparseCreateCsrsm2Info(&gpus[g].h_array_info_U[j]));
                    }
                }
                POP_RANGE // LSC
            } // LSC loop

			checkCudaErrors(cudaDeviceSynchronize());
		} // OMP parallel region

		// Free memory
		PUSH_RANGE("Free", 6)
		for (int g = 0; g < n_gpu; g++) {
            checkCudaErrors(cudaSetDevice(g + gpu_id));
            checkCudaErrors(cudaStreamDestroy(gpus[g].data_stream));

            for (int s = 0; s < n_streams_per_gpu; s++) {
                checkCudaErrors(cudaStreamDestroy(gpus[g].h_array_stream[s]));
                checkCudaErrors(cusparseDestroy(gpus[g].h_array_handle[s]));
                checkCudaErrors(cudaFree(gpus[g].h_array_d_buffer[s]));
                checkCudaErrors(cudaFree(gpus[g].h_array_d_X_reordered[s]));
                checkCudaErrors(cudaFree(gpus[g].h_array_d_csr_B_col_ind[s]));
                checkCudaErrors(cudaFree(gpus[g].h_array_d_csr_B_row_ptr[s]));
                checkCudaErrors(cudaFree(gpus[g].h_array_d_csr_B_val[s]));
                checkCudaErrors(cudaFree(gpus[g].h_array_d_csr_Bt_col_ind[s]));
                checkCudaErrors(cudaFree(gpus[g].h_array_d_csr_Bt_row_ptr[s]));
                checkCudaErrors(cudaFree(gpus[g].h_array_d_csr_Bt_val[s]));

				if(SYMMETRIC_SYSTEM) {
					if(order) {
						checkCudaErrors(cudaFree(gpus[g].h_array_d_pinv[s]));
					}
				} else {
					if(order) {
						checkCudaErrors(cudaFree(gpus[g].h_array_d_q[s]));
					}
					checkCudaErrors(cudaFree(gpus[g].h_array_d_pinv[s]));
				}
            }
            checkCudaErrors(cudaEventDestroy(gpus[g].event_data_preload));
            checkCudaErrors(cudaEventDestroy(gpus[g].event1));
            checkCudaErrors(cudaEventDestroy(gpus[g].event2));

            for (int i = 0; i < n_csrsm2_info_per_gpu; i++) {
                checkCudaErrors(cusparseDestroyCsrsm2Info(gpus[g].h_array_info_L[i]));
                checkCudaErrors(cusparseDestroyCsrsm2Info(gpus[g].h_array_info_U[i]));
            }
            checkCudaErrors(cudaFree(gpus[g].d_csc_L_col_ptr));
            checkCudaErrors(cudaFree(gpus[g].d_csc_L_row_ind));
            checkCudaErrors(cudaFree(gpus[g].d_csc_L_val));
            checkCudaErrors(cudaFree(gpus[g].d_csc_U_col_ptr));
            checkCudaErrors(cudaFree(gpus[g].d_csc_U_row_ind));
            checkCudaErrors(cudaFree(gpus[g].d_csc_U_val));

            free(gpus[g].h_array_stream);
            free(gpus[g].h_array_handle);
            free(gpus[g].h_array_info_L);
            free(gpus[g].h_array_info_U);
            free(gpus[g].h_array_d_csr_B_row_ptr);
            free(gpus[g].h_array_d_csr_B_col_ind);
            free(gpus[g].h_array_d_csr_B_val);
            free(gpus[g].h_array_d_csr_Bt_row_ptr);
            free(gpus[g].h_array_d_csr_Bt_col_ind);
            free(gpus[g].h_array_d_csr_Bt_val);
            free(gpus[g].h_array_d_buffer);
            free(gpus[g].h_array_d_X_reordered);
            free(gpus[g].h_array_d_lsc);
            free(gpus[g].h_array_lsc_id);
        }
        
        checkCudaErrors(cusparseDestroyMatDescr(h_descrL));
        checkCudaErrors(cusparseDestroyMatDescr(h_descrU));
        checkCudaErrors(cusparseDestroyMatDescr(h_descrBt));

        free(h_array_h_matB);
        free(h_array_h_matX);
        free(h_array_h_matLSC);
        free(gpus);

		// Free factors
		#pragma omp parallel for
		for (esint d = 0; d < n_lsc; d++) {
			if (SYMMETRIC_SYSTEM) {
				csparse::FreeCholFactor(vec_L_row_indexes[d], vec_L_col_pointers[d], vec_L_values[d], vec_perm[d]);
			} else {
				csparse::FreeLuFactors(vec_L_row_indexes[d], vec_L_col_pointers[d], vec_L_values[d],
				vec_U_row_indexes[d], vec_U_col_pointers[d], vec_U_values[d], vec_perm[d], vec_perm_2[d]);
			}
		}
		POP_RANGE // Free
	}
}


void ClusterGPU::GetSchurComplementsCpu(bool USE_FLOAT) {
	#pragma omp parallel for
	for (size_t d = 0; d < domains_in_global_index.size(); d++ ) {
		if (domains[d].B1Kplus.is_on_acc == 0 && !configuration.combine_sc_and_spds) {
			// Calculates SC on CPU and keeps it in CPU memory
			GetSchurComplement(USE_FLOAT, d);
			// ESINFO(PROGRESS3) << Info::plain() << ".";
		}
	}
}


void ClusterGPU::CreateDirichletPrec( DataHolder *instance) {
	DEBUGOUT << "Creating Dirichlet Preconditioner with Pardiso SC and copying them to GPU\n";

	GetGPU();

	esint status = 0;

	std::vector<size_t> local_Prec_size_to_add(domains_in_global_index.size(), 0);

	esint domains_on_GPU = 0;
	esint domains_on_CPU = 0;
	esint DOFs_GPU = 0;
	esint DOFs_CPU = 0;

	// Smycka napocitava velikost LSCs pres vsechny domeny
	for (size_t d = 0; d < domains_in_global_index.size(); d++) {
		if (domains[d].Prec.USE_FLOAT) {
			local_Prec_size_to_add[d] =
					( domains[d].B1t_Dir_perm_vec.size() * domains[d].B1t_Dir_perm_vec.size()
					) * sizeof(float);
		} else {
			local_Prec_size_to_add[d] =
					( domains[d].B1t_Dir_perm_vec.size() * domains[d].B1t_Dir_perm_vec.size()
					) * sizeof(double);
		}
		// If B1Kplus not on GPU we also need 2 vectors
		if (domains[d].B1Kplus.is_on_acc == 0) {
			if (domains[d].Prec.USE_FLOAT) {
				local_Prec_size_to_add[d] +=
						( 2 * domains[d].B1t_Dir_perm_vec.size()
						) * sizeof(float);
			} else {
				local_Prec_size_to_add[d] +=
						( 2 * domains[d].B1t_Dir_perm_vec.size()
						) * sizeof(double);
			}
		}

		if(local_Prec_size_to_add[d] < GPU_free_mem) {
			domains_on_GPU++;
			DOFs_GPU += domains[d].Prec.rows;
			domains[d].Prec.is_on_acc = 1;
			GPU_free_mem -= local_Prec_size_to_add[d];
		} else {
			domains_on_CPU++;
			DOFs_CPU += domains[d].Prec.rows;
			domains[d].Prec.is_on_acc = 0;
		}

		/* OVERKILL PART 2
		// Collect info on how much each process wants to add
		MPI_Gather(&local_SC_size_to_add, 1, MPI_INT, &SC_size_to_add[0], 1, MPI_INT, 0, node_comm);

		if(local_id == 0)
		{   // Proc with id 0 will decide if it fits
		      for(int proc = 0; proc < local_procs; proc++)
		      {
			      if (SC_total_size[GPU_mapping[proc]] + SC_size_to_add[proc] < GPU_free_mem[GPU_mapping[proc]]) {
				      msg[proc] = 1;
				      SC_total_size[GPU_mapping[proc]] += SC_size_to_add[proc];
			      } else {
				      msg[proc] = 0;
			      }
		      }
		}
		// Get the decision from proc 0
		MPI_Scatter(&msg[0], 1, MPI_INT, &reply, 1, MPI_INT, 0, node_comm);

		if(reply)
		{
		      domains_on_GPU++;
		      DOFs_GPU += domains[d].K.rows;
		      domains[d].B1Kplus.is_on_acc = 1;
		}
		else
		{
		      domains_on_CPU++;
		      DOFs_CPU += domains[d].K.rows;
		      domains[d].B1Kplus.is_on_acc = 0;
		}
	*/
	}

	// TODO_GPU - vsechny tyto std::cout se musi prepsat na logovani co ma Ondra M.
	// Ondro nektere moje rutiny, napr. SpyText jsou napsane pro std::cout a ne printf. Jake je reseni ?

//	std::vector <int> on_gpu (info::mpi::size, 0);
//	MPI_Gather(&domains_on_GPU,1,MPI_INT,&on_gpu[0],1,MPI_INT, 0, info::mpi::comm);
//
//	std::vector <int> on_cpu (info::mpi::size, 0);
//	MPI_Gather(&domains_on_CPU,1,MPI_INT,&on_cpu[0],1,MPI_INT, 0, info::mpi::comm);
//
//	std::vector <int> don_gpu (info::mpi::size, 0);
//	MPI_Gather(&DOFs_GPU,1,MPI_INT,&don_gpu[0],1,MPI_INT, 0, info::mpi::comm);
//
//	std::vector <int> don_cpu (info::mpi::size, 0);
//	MPI_Gather(&DOFs_CPU,1,MPI_INT,&don_cpu[0],1,MPI_INT, 0, info::mpi::comm);
//
//	for (esint i = 0; i < info::mpi::size; i++) {
//		eslog::info(" MPI rank %d\t - GPU: Precs = %d Total DOFs = %d\t - CPU: Precs = %d Total DOFs = %d\n",
//		 i, on_gpu[i], don_gpu[i], on_cpu[i], don_cpu[i]);
//	}

	esint info[2] = { domains_on_GPU, domains_on_CPU + domains_on_GPU };
	Communication::allReduce(info, NULL, 2, MPITools::getType<esint>().mpitype, MPI_SUM);
	std::string ratio = Parser::stringwithcommas(info[0]) + " / " + Parser::stringwithcommas(info[1]);
	eslog::solver("     - | ACCELERATED PRECONDITIONERS %49s | -\n", ratio.c_str());

	#pragma omp parallel for
	for (size_t d = 0; d < domains_in_global_index.size(); d++ ) {
		// Calculates Prec on CPU and keeps it CPU memory
		// Base solver method
		CreateDirichletPrec_perDomain(instance, d);
//		std::cout << ".";
	}
//	if (info::mpi::rank == 0) {
//		std::cout << "\n";
//	}

	CreateCudaStreamPool();

	for (size_t d = 0; d < domains_in_global_index.size(); d++) {
		esint status = 0;

		if (domains[d].Prec.is_on_acc == 1) {
			// TODO_GPU: Test performance (same streams)
			// Assign CUDA stream from he pool
			domains[d].Prec.stream = cuda_stream_pool[d % configuration.num_streams];
			domains[d].Prec.handle = cublas_handle_pool[d % configuration.num_streams];

			if (domains[d].Prec.USE_FLOAT){
				status = domains[d].Prec.CopyToCUDA_Dev_fl();
				// B1Kplus is not on GPU, threfore vectors arent alocated
				if (domains[d].B1Kplus.is_on_acc == 0) {
					size_t memsize = domains[d].Prec.rows;
					// TODO_GPU - alokuji se pole pro vektory, kterymi se bude ve funkci Apply_Prec nasobit tato matice
					domains[d].cuda_pinned_buff_fl.resize(memsize);
					cudaHostRegister(domains[d].cuda_pinned_buff_fl.data(), memsize * sizeof(float), cudaHostRegisterDefault);
				}
			} else {
				// Share already allocated device vectors with B1Kplus after move from SparseMatrix to Domain
				// Temporal solution
				if(domains[d].B1Kplus.rows >= domains[d].Prec.rows) {
					domains[d].Prec.d_x_in = domains[d].B1Kplus.d_x_in;
					domains[d].Prec.d_y_out = domains[d].B1Kplus.d_y_out;
				}
				status = domains[d].Prec.CopyToCUDA_Dev();
				// B1Kplus is not on GPU, threfore vectors arent alocated
				if (domains[d].B1Kplus.is_on_acc == 0) {
					size_t memsize = domains[d].Prec.rows;

					domains[d].cuda_pinned_buff.resize(memsize);
					cudaHostRegister(domains[d].cuda_pinned_buff.data(), memsize * sizeof(double), cudaHostRegisterDefault);
				}
			}

			if (status == 0) {
                                // TODO_GPU - Precs ktere se uspesne prenesly do pameti GPU se smazou z pameti procesoru

				if (domains[d].Prec.USE_FLOAT) {
					SEQ_VECTOR <float>  ().swap (domains[d].Prec.dense_values_fl);

					// std::cout << "g";
				} else {
					SEQ_VECTOR <double> ().swap (domains[d].Prec.dense_values);
				}

					// std::cout << "G";
			} else {
				// pokud se domenu nepodar nahrat na GPU
				domains[d].Prec.is_on_acc = 0;
				domains_on_CPU++;
				domains_on_GPU--;
				DOFs_CPU += domains[d].K.rows;
				DOFs_GPU -= domains[d].K.rows;
				GPU_free_mem += local_Prec_size_to_add[d];

				domains[d].Prec.ClearCUDA_Stream();

				if(domains[d].Prec.USE_FLOAT) {
					domains[d].Prec.FreeFromCUDA_Dev_fl();
				} else {
					domains[d].Prec.FreeFromCUDA_Dev();
				}

				if(domains[d].Prec.USE_FLOAT) {
					if (domains[d].B1Kplus.is_on_acc == 0) {
						cudaHostUnregister(domains[d].cuda_pinned_buff_fl.data());
						SEQ_VECTOR <float>  ().swap (domains[d].cuda_pinned_buff_fl);
					}
				} else {
					if (domains[d].B1Kplus.is_on_acc == 0) {
						cudaHostUnregister(domains[d].cuda_pinned_buff.data());
						SEQ_VECTOR <double>  ().swap (domains[d].cuda_pinned_buff);
					}
				}
			}
		}

		// std::cout << " Domain: " << d << " GPU : " << domains[d].Prec.is_on_acc << "\n";
	}
//	eslog::info("\n Precs moved to GPU: %d\n", domains_on_GPU);
//	eslog::info(" Precs remaining on CPU: %d\n", domains_on_CPU);
}


void ClusterGPU::SetupKsolvers ( ) {
	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++) {

		// Import of Regularized matrix K into Kplus (Sparse Solver)
		switch (configuration.Ksolver) {
		case FETIConfiguration::KSOLVER::DIRECT_DP:
			domains[d].Kplus.ImportMatrix_wo_Copy (domains[d].K);
			break;
		case FETIConfiguration::KSOLVER::ITERATIVE:
			domains[d].Kplus.ImportMatrix_wo_Copy (domains[d].K);
			break;
		case FETIConfiguration::KSOLVER::DIRECT_SP:
			domains[d].Kplus.ImportMatrix_wo_Copy_fl(domains[d].K);
			//domains[d].Kplus.ImportMatrix_fl(domains[d].K);
			break;
		case FETIConfiguration::KSOLVER::DIRECT_MP:
			domains[d].Kplus.ImportMatrix_fl(domains[d].K);
			break;
//		case 4:
//			domains[d].Kplus.ImportMatrix_fl(domains[d].K);
//			break;
		default:
			////ESINFO(ERROR) << "Invalid KSOLVER value.";
			exit(EXIT_FAILURE);
		}

		if (configuration.keep_factors) {

			if (!configuration.combine_sc_and_spds) { // if both CPU and GPU uses Schur Complement
				std::stringstream ss;
				ss << "init -> rank: " << info::mpi::rank << ", subdomain: " << d;
				domains[d].Kplus.keep_factors = true;
				if (configuration.Ksolver != FETIConfiguration::KSOLVER::ITERATIVE) {
					domains[d].Kplus.Factorization (ss.str());
				}
			} else {
				if ( domains[d].B1Kplus.is_on_acc == 0 ) {
					std::stringstream ss;
					ss << "init -> rank: " << info::mpi::rank << ", subdomain: " << d;
					domains[d].Kplus.keep_factors = true;
					if (configuration.Ksolver != FETIConfiguration::KSOLVER::ITERATIVE) {
						domains[d].Kplus.Factorization (ss.str());
					}
				}
			}

		} else {
			domains[d].Kplus.keep_factors = false;
			//domains[d].Kplus.rank = info::mpi::rank;
		}

		if ( d == 0 && info::mpi::rank == 0) {
			domains[d].Kplus.msglvl = 0; //Info::report(LIBRARIES) ? 1 : 0;
		}
	}
}


void ClusterGPU::multKplusGlobal_GPU(SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {
	mkl_set_num_threads(1);

	cluster_time.totalTime.start();

	vec_fill_time.start();
	fill(vec_g0.begin(), vec_g0.end(), 0); // reset entire vector to 0
	vec_fill_time.end();

	// loop over domains in the cluster
	loop_1_1_time.start();
	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++)
	{
		domains[d].B0Kplus_comp.DenseMatVec(x_in[d], tm2[d]);			// g0 - with comp B0Kplus
		domains[d].Kplus_R.DenseMatVec(x_in[d], tm3[d], 'T');			// e0
	}
	loop_1_1_time.end();

	loop_1_2_time.start();
	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++)
	{
		esint e0_start	=  d	* domains[d].Kplus_R.cols;
		esint e0_end		= (d+1) * domains[d].Kplus_R.cols;

		for (esint i = e0_start; i < e0_end; i++ )
			vec_e0[i] = - tm3[d][i - e0_start];
	}

	for (size_t d = 0; d < domains.size(); d++)
		for (esint i = 0; i < domains[d].B0Kplus_comp.rows; i++)
			vec_g0[ domains[d].B0_comp_map_vec[i] - 1 ] += tm2[d][i];

	// end loop over domains
	loop_1_2_time.end();

	mkl_set_num_threads(PAR_NUM_THREADS);
	clusCP_time.start();

	clus_F0_1_time.start();
	F0.Solve(vec_g0, tm1[0], 0, 0);
	clus_F0_1_time.end();

	clus_G0_time.start();
	G0.MatVec(tm1[0], tm2[0], 'N');
	clus_G0_time.end();

	#pragma omp parallel for
	for (size_t i = 0; i < vec_e0.size(); i++)
		tm2[0][i] = tm2[0][i] - vec_e0[i];
	//cblas_daxpy(vec_e0.size(), -1.0, &vec_e0[0], 1, &tm2[0][0], 1);

	 clus_Sa_time.start();
//#ifdef SPARSE_SA
//	 Sa.Solve(tm2[0], vec_alfa,0,0);
//#else
//	char U = 'U';
//	esint nrhs = 1;
//	esint info = 0;
//	vec_alfa = tm2[0];
//	dsptrs( &U, &SaMat.rows, &nrhs, &SaMat.dense_values[0], &SaMat.ipiv[0], &vec_alfa[0], &SaMat.rows, &info );
//#endif
//
////TODO: New version with DenseSolver routine
////#ifdef SPARSE_SA
////	Sa.Solve(tm2[0], vec_alfa,0,0);
////#else
////	esint nrhs = 1;
////	Sa_dense.Solve(tm2[0], vec_alfa, nrhs);
////#endif

	 if (configuration.SAsolver == FETIConfiguration::SASOLVER::CPU_SPARSE) {
		 Sa.Solve(tm2[0], vec_alfa,0,0);
	 }

	 if (configuration.SAsolver == FETIConfiguration::SASOLVER::CPU_DENSE) {
			esint nrhs = 1;
			Sa_dense_cpu.Solve(tm2[0], vec_alfa, nrhs);
	 }

	 if (configuration.SAsolver == FETIConfiguration::SASOLVER::ACC_DENSE) {
			esint nrhs = 1;
			Sa_dense_acc.Solve(tm2[0], vec_alfa, nrhs);
	 }

	 clus_Sa_time.end();

	 clus_G0t_time.start();
	G0.MatVec(vec_alfa, tm1[0], 'T'); 	// lambda
	 clus_G0t_time.end();

	#pragma omp parallel for
	for (size_t i = 0; i < vec_g0.size(); i++)
		tm1[0][i] = vec_g0[i] - tm1[0][i];

	clus_F0_2_time.start();
	F0.Solve(tm1[0], vec_lambda,0,0);
	clus_F0_2_time.end();

	clusCP_time.end();

	// Kplus_x
	mkl_set_num_threads(1);
	loop_2_1_time.start();
	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++)
	{
		esint domain_size = domains[d].domain_prim_size;
		SEQ_VECTOR < double > tmp_vec (domains[d].B0_comp_map_vec.size(), 0.0);

		bool MIXED_SC_FACT = configuration.combine_sc_and_spds;

		if (domains[d].B1Kplus.is_on_acc == 0 && MIXED_SC_FACT) {
		// 4. 3. 2021 Ineffective computation
		// if (domains[d].B0Kplus_comp.is_on_acc == 0 && MIXED_SC_FACT) {

		// 	for (esint i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
		// 		tmp_vec[i] = vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;

		// 	domains[d].B0_comp.MatVec(tmp_vec, tm1[d], 'T');

		// 	for (esint i = 0; i < domain_size; i++)
		// 		tm1[d][i] = x_in[d][i] - tm1[d][i];

		// 	domains[d].multKplusLocal(tm1[d] , tm2[d]);
		// } else {
		for (size_t i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
			tmp_vec[i] = -1.0 * vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;

		domains[d].B0Kplus_comp.DenseMatVec(tmp_vec, tm2[d], 'T' );
		// }

		esint e0_start	=  d	* domains[d].Kplus_R.cols;
		esint e0_end		= (d+1) * domains[d].Kplus_R.cols;
		domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);

		for (esint i = 0; i < domain_size; i++)
			x_in[d][i] = tm2[d][i] + tm3[d][i];
		}
	}

	loop_2_1_time.end();

	cluster_time.totalTime.end();
}


void ClusterGPU::multKplus_HF(SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {

	cluster_time.totalTime.start();

	multKplus_HF_Loop1(x_in);

	multKplus_HF_CP();

	multKplus_HF_Loop2_MIX(x_in);

	cluster_time.totalTime.end();
}


void ClusterGPU::multKplus_HF_SC(SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {

	cluster_time.totalTime.start();

	multKplus_HF_Loop1(x_in);

	multKplus_HF_CP();

	multKplus_HF_Loop2_SC(x_in, x_in);

	cluster_time.totalTime.end();
}


void ClusterGPU::multKplus_HF_SC(SEQ_VECTOR<SEQ_VECTOR<double> > & x_in, SEQ_VECTOR<SEQ_VECTOR<double> > & y_out) {

	cluster_time.totalTime.start();

	multKplus_HF_Loop1(x_in);

	multKplus_HF_CP();

	multKplus_HF_Loop2_SC(x_in, y_out);

	cluster_time.totalTime.end();
}


void ClusterGPU::multKplus_HF_SPDS(SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {

	cluster_time.totalTime.start();

	multKplus_HF_Loop1(x_in);

	multKplus_HF_CP();

	multKplus_HF_Loop2_SPDS(x_in);

	cluster_time.totalTime.end();
}


void ClusterGPU::multKplus_HF_Loop1(SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {

	//mkl_set_num_threads(1);

	vec_fill_time.start();
	fill(vec_g0.begin(), vec_g0.end(), 0); // reset entire vector to 0
	vec_fill_time.end();

	// loop over domains in the cluster
	loop_1_1_time.start();
	loop_1_1_time.end();
	loop_1_2_time.start();
	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++)
	{

		domains[d].B0Kplus_comp.DenseMatVec(x_in[d], tm2[d]);			// g0 - with comp B0Kplus
		domains[d].Kplus_R.     DenseMatVec(x_in[d], tm3[d], 'T');	    // e0

		esint e0_start	=  d	* domains[d].Kplus_R.cols;
		esint e0_end		= (d+1) * domains[d].Kplus_R.cols;

		for (esint i = e0_start; i < e0_end; i++ )
			vec_e0[i] = - tm3[d][i - e0_start];
	}

	for (size_t d = 0; d < domains.size(); d++)
		for (esint i = 0; i < domains[d].B0Kplus_comp.rows; i++)
			vec_g0[ domains[d].B0_comp_map_vec[i] - 1 ] += tm2[d][i];

	loop_1_2_time.end();

}

void ClusterGPU::multKplus_HF_CP( ) {

	//mkl_set_num_threads(PAR_NUM_THREADS);
	 clusCP_time.start();

	 clus_F0_1_time.start();
	F0.Solve(vec_g0, tm1[0], 0, 0);
	 clus_F0_1_time.end();

	 clus_G0_time.start();
	G0.MatVec(tm1[0], tm2[0], 'N');
	 clus_G0_time.end();

	#pragma omp parallel for
	for (size_t i = 0; i < vec_e0.size(); i++)
		tm2[0][i] = tm2[0][i] - vec_e0[i];
	//cblas_daxpy(vec_e0.size(), -1.0, &vec_e0[0], 1, &tm2[0][0], 1);


	 clus_Sa_time.start();
	if (configuration.SAsolver == FETIConfiguration::SASOLVER::CPU_SPARSE) {
		Sa.Solve(tm2[0], vec_alfa,0,0);
	}

	if (configuration.SAsolver == FETIConfiguration::SASOLVER::CPU_DENSE) {
		esint nrhs = 1;
		Sa_dense_cpu.Solve(tm2[0], vec_alfa, nrhs);
	}

	if (configuration.SAsolver == FETIConfiguration::SASOLVER::ACC_DENSE) {
		esint nrhs = 1;
		Sa_dense_acc.Solve(tm2[0], vec_alfa, nrhs);// lambda
	}
	 clus_Sa_time.end();


	 clus_G0t_time.start();
	G0.MatVec(vec_alfa, tm1[0], 'T'); 	// lambda
	 clus_G0t_time.end();

	#pragma omp parallel for
	for (size_t i = 0; i < vec_g0.size(); i++)
		tm1[0][i] = vec_g0[i] - tm1[0][i];


	 clus_F0_2_time.start();
	F0.Solve(tm1[0], vec_lambda,0,0);
	 clus_F0_2_time.end();

	 clusCP_time.end();

}

void ClusterGPU::multKplus_HF_Loop2_SC(SEQ_VECTOR<SEQ_VECTOR<double> > & x_in, SEQ_VECTOR<SEQ_VECTOR<double> > & y_out) {

	// Kplus_x
	//mkl_set_num_threads(1);
	 loop_2_1_time.start();

	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++)
	{
		esint domain_size = domains[d].domain_prim_size;
		SEQ_VECTOR < double > tmp_vec (domains[d].B0_comp_map_vec.size(), 0.0);

		for (size_t i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
			tmp_vec[i] = -1.0 * vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;

		domains[d].B0Kplus_comp.DenseMatVec(tmp_vec, tm2[d], 'T' );

		esint e0_start	=  d	* domains[d].Kplus_R.cols;
		esint e0_end		= (d+1) * domains[d].Kplus_R.cols;

		domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);

		for (esint i = 0; i < domain_size; i++)
			y_out[d][i] = tm2[d][i] + tm3[d][i];

	}
	 loop_2_1_time.end();

}

void ClusterGPU::multKplus_HF_Loop2_SPDS (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {

	// Kplus_x
	//mkl_set_num_threads(1);
	 loop_2_1_time.start();

	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++)
	{
		esint domain_size = domains[d].domain_prim_size;
		SEQ_VECTOR < double > tmp_vec (domains[d].B0_comp_map_vec.size(), 0.0);

		for (size_t i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
			tmp_vec[i] = vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;

		domains[d].B0_comp.MatVec(tmp_vec, tm1[d], 'T');

		for (esint i = 0; i < domain_size; i++)
			tm1[d][i] = x_in[d][i] - tm1[d][i];

		domains[d].multKplusLocal(tm1[d] , tm2[d]);

		esint e0_start	=  d	* domains[d].Kplus_R.cols;
		esint e0_end		= (d+1) * domains[d].Kplus_R.cols;

		domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);


		for (esint i = 0; i < domain_size; i++)
			x_in[d][i] = tm2[d][i] + tm3[d][i];

	}
	 loop_2_1_time.end();

}

void ClusterGPU::multKplus_HF_Loop2_MIX (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {

	// Kplus_x
	// mkl_set_num_threads(1);
	 loop_2_1_time.start();

	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++)
	{
		esint domain_size = domains[d].domain_prim_size;
		SEQ_VECTOR < double > tmp_vec (domains[d].B0_comp_map_vec.size(), 0.0);

		bool MIXED_SC_FACT = configuration.combine_sc_and_spds;

		if ( (domains[d].B1Kplus.is_on_acc == 0 && MIXED_SC_FACT ) || !configuration.use_schur_complement ) {

			for (size_t i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
				tmp_vec[i] = vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;

			domains[d].B0_comp.MatVec(tmp_vec, tm1[d], 'T');

			for (esint i = 0; i < domain_size; i++)
				tm1[d][i] = x_in[d][i] - tm1[d][i];

			domains[d].multKplusLocal(tm1[d] , tm2[d]);

		// 4. 3. 2021 Ineffective computation
        // if ( (domains[d].B0Kplus_comp.is_on_acc == 0 && MIXED_SC_FACT ) || !configuration.use_schur_complement ) {
			// for (esint i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
			// 	tmp_vec[i] = vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;

			// domains[d].B0_comp.MatVec(tmp_vec, tm1[d], 'T');

			// for (esint i = 0; i < domain_size; i++)
			// 	tm1[d][i] = x_in[d][i] - tm1[d][i];

			// domains[d].multKplusLocal(tm1[d] , tm2[d]);

			// esint e0_start	=  d	* domains[d].Kplus_R.cols;
			// esint e0_end		= (d+1) * domains[d].Kplus_R.cols;

			// domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);
		// } else {
		for (size_t i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
			tmp_vec[i] = -1.0 * vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;

		domains[d].B0Kplus_comp.DenseMatVec(tmp_vec, tm2[d], 'T' );

		esint e0_start	=  d	* domains[d].Kplus_R.cols;
		esint e0_end		= (d+1) * domains[d].Kplus_R.cols;

		domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);
		// }

		for (esint i = 0; i < domain_size; i++)
			x_in[d][i] = tm2[d][i] + tm3[d][i];

	}
	 loop_2_1_time.end();
	}
}


void ClusterGPU::CreateCudaStreamPool() {
	if (cuda_stream_pool.size() < configuration.num_streams) {
		DestroyCudaStreamPool();

		cuda_stream_pool.resize(configuration.num_streams);
		cublas_handle_pool.resize(configuration.num_streams);

		for(size_t i = 0; i < configuration.num_streams; i++) {
			checkCudaErrors(cudaStreamCreate(&cuda_stream_pool[i]));
			checkCudaErrors(cublasCreate(&cublas_handle_pool[i]));
			checkCudaErrors(cublasSetStream(cublas_handle_pool[i], cuda_stream_pool[i]));
		}
		////ESINFO(VERBOSE_LEVEL3) << "CUDA stream pool created with " << configuration.num_streams << " streams per cluster";
	}
}


void ClusterGPU::DestroyCudaStreamPool() {
	for(size_t i = 0; i < cuda_stream_pool.size(); i++) {
		checkCudaErrors(cudaStreamDestroy(cuda_stream_pool[i]));
		checkCudaErrors(cublasDestroy(cublas_handle_pool[i]));
	}

	SEQ_VECTOR <cudaStream_t>().swap(cuda_stream_pool);
	SEQ_VECTOR <cublasHandle_t>().swap(cublas_handle_pool);
}
