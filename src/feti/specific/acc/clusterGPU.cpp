#include "clusterGPU.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/mpiinfo.h"
#include <sstream>
#include <iostream>
#include <algorithm>
#include "mkl.h"

using namespace espreso;

//#define STREAM_NUM 16
//#define SHARE_SC

ClusterGPU::~ClusterGPU() {
	DestroyCudaStreamPool();

	for(esint d = 0; d < domains_in_global_index.size(); d++) {
                if(domains[d].isOnACC) {
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
#endif
		}
	}
}

void ClusterGPU::GetAvailableGPUmemory() {

	bool GPU_full = false;
	//GPU_full = true;
	int nDevices;
	cudaGetDeviceCount(&nDevices);

	std::cout << "\n*** GPU Accelerators available on the server " << "\n\n";
	for (int i = 0; i < nDevices; i++) {
		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, i);
		std::cout << " Device Number: " << i << "\n";
		std::cout << " Device name: " << prop.name << "\n";
		std::cout << " Memory Clock Rate (KHz): " << prop.memoryClockRate << "\n";
		std::cout << " Memory Bus Width (bits): " << prop.memoryBusWidth << "\n";
		std::cout << " Peak Memory Bandwidth (GB/s): " << 2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6 << "\n";

		cudaSetDevice(i);
		size_t free, total;
		cudaMemGetInfo(&free, &total);
		std::cout << " GPU Total Memory [MB]: " << total/1024/1024 << "\n";
		std::cout << " GPU Free Memory [MB]:  " << free/1024/1024 << "\n\n";

	}

	// TODO_GPU
	// - zde se rohoduje, na ktere GPU tento MPI proces pouziva
	// Faze 1 - 1 MPI process pouziva 1 GPU
	//		  - napsat kod, ktere si detekuje kolim MPI ranku je na uzlu a podle toho priradi min. 1 nebo vice MPI procesu na kazde GPU
	// Faze 2 - napsat podporu pro vice GPU na 1 MPI process

	// GPU memory management
	// Create new communicator within the node (OMPI_COMM_TYPE_NODE can be swapped out with MPI_COMM_TYPE_SHARED for portability)
	MPI_Comm node_comm;
	MPI_Comm_split_type(info::mpi::comm, MPI_COMM_TYPE_SHARED, info::mpi::rank, MPI_INFO_NULL, &node_comm);

	// Get local size and id
	int local_procs;
	MPI_Comm_size(node_comm, &local_procs);

	int local_id;
	MPI_Comm_rank(node_comm, &local_id);

	size_t procs_per_gpu;
	int device_id;
	if ((local_procs % nDevices) != 0)
	{
	  std::cout<<" Only integer multiply number of processes per GPU. Processes: "<< local_procs << " GPUs: "<< nDevices << "\n";
	  exit(0);
	}
	else
	{
	  procs_per_gpu = local_procs / nDevices;
	  device_id     = local_id    / procs_per_gpu;
	}

	cudaSetDevice(device_id);
	cudaMemGetInfo(&GPU_free_mem, &GPU_total_mem);
	GPU_free_mem  /= procs_per_gpu;
	GPU_total_mem /= procs_per_gpu;

	/* OVERKILL PART 1
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

void ClusterGPU::Create_SC_perDomain(bool USE_FLOAT) {

	GetAvailableGPUmemory();
	std::cout << "Creating B1*K+*B1t Schur Complements with Pardiso SC and coping them to GPU";

	esint status = 0;
	cudaError_t status_c;

	std::vector<size_t> local_SC_size_to_add (domains_in_global_index.size(), 0);

	esint domains_on_GPU = 0;
	esint domains_on_CPU = 0;
	esint DOFs_GPU = 0;
	esint DOFs_CPU = 0;


	// Smycka napocitava velikost LSCs pres vsechny domeny
	for (esint d = 0; d < domains_in_global_index.size(); d++ ) {

		switch (configuration.schur_type) {
		case FETIConfiguration::MATRIX_STORAGE::GENERAL:

	// TODO_GPU - SHARE_SC asi nepatri sem pro general matrix, ale patri do symetrickych matic nize 
	// Radime potvrd.
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
						  2 * domains[d].B1_comp_dom.rows
						) * sizeof(float);
			} else {
				local_SC_size_to_add[d] =
						( domains[d].B1_comp_dom.rows * domains[d].B1_comp_dom.rows +
						  2 * domains[d].B1_comp_dom.rows
						) * sizeof(double);
			}
#endif
			break;
		case FETIConfiguration::MATRIX_STORAGE::SYMMETRIC:
			if (USE_FLOAT) {
				local_SC_size_to_add[d] =
						(((domains[d].B1_comp_dom.rows + 1 ) * domains[d].B1_comp_dom.rows ) / 2
						 + 2 * domains[d].B1_comp_dom.rows
						) * sizeof(float);
			} else {
				local_SC_size_to_add[d] =
						(((domains[d].B1_comp_dom.rows + 1 ) * domains[d].B1_comp_dom.rows ) / 2
						 + 2 * domains[d].B1_comp_dom.rows
						) * sizeof(double);
			}
			break;
		default:
			break;
			std::cout << "ERROR - Not implemented type of Schur complement.";
		}

		if(local_SC_size_to_add[d] < GPU_free_mem)
		{
			domains_on_GPU++;
			DOFs_GPU += domains[d].K.rows;
			domains[d].B1Kplus.isOnACC = 1;
			GPU_free_mem -= local_SC_size_to_add[d];
		}
		else
		{
			domains_on_CPU++;
			DOFs_CPU += domains[d].K.rows;
			domains[d].B1Kplus.isOnACC = 0;
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
                        domains[d].isOnACC = 1;
		}
		else
		{
			domains_on_CPU++;
			DOFs_CPU += domains[d].K.rows;
                        domains[d].isOnACC = 0;
		}
		*/
	}

	// TODO_GPU - vsechny tyto std::cout se musi prepsat na logovani co ma Ondra M. 
	// Ondro nektere moje rutiny, napr. SpyText jsou napsane pro std::cout a ne printf. Jake je reseni ? 
	std::cout << "\n Domains on GPU : " << domains_on_GPU << "\n";
	std::cout << " Domains on CPU : " << domains_on_CPU << "\n";

	std::vector <int> on_gpu (info::mpi::size, 0);
	MPI_Gather(&domains_on_GPU,1,MPI_INT,&on_gpu[0],1,MPI_INT, 0, info::mpi::comm);

	std::vector <int> on_cpu (info::mpi::size, 0);
	MPI_Gather(&domains_on_CPU,1,MPI_INT,&on_cpu[0],1,MPI_INT, 0, info::mpi::comm);

	std::vector <int> don_gpu (info::mpi::size, 0);
	MPI_Gather(&DOFs_GPU,1,MPI_INT,&don_gpu[0],1,MPI_INT, 0, info::mpi::comm);

	std::vector <int> don_cpu (info::mpi::size, 0);
	MPI_Gather(&DOFs_CPU,1,MPI_INT,&don_cpu[0],1,MPI_INT, 0, info::mpi::comm);


	for (esint i = 0; i < info::mpi::size; i++) {
		std::cout << " MPI rank " << i <<
			"\t - GPU : domains = \t" << on_gpu[i] << "\t Total DOFs = \t" << don_gpu[i] <<
			"\t - CPU : domains = \t" << on_cpu[i] << "\t Total DOFs = \t" << don_cpu[i] << "\n";
	}

#ifdef SHARE_SC
	SEQ_VECTOR <esint> SC_dense_val_offsets(domains_in_global_index.size(), 0);

	// 2 domains per iteration processed
	#pragma omp parallel for
	for (esint d = 0; d < domains_in_global_index.size(); d += 2 ) {

		if (domains[d].isOnACC == 1 || !configuration.combine_sc_and_spds) {
			// Calculates SC on CPU and keeps it CPU memory
			GetSchurComplement(USE_FLOAT, d);
			std::cout << Info::plain() << ".";

			// Set if Upper or Lower part is referenced
			domains[d].B1Kplus.uplo = 'U';

			// Set the default lda
			domains[d].B1Kplus.extern_lda = domains[d].B1Kplus.rows;
		}

		if (d+1 < domains_in_global_index.size() && (domains[d+1].isOnACC == 1 || !configuration.combine_sc_and_spds)) {
			// Calculates SC on CPU and keeps it CPU memory
			GetSchurComplement(USE_FLOAT, d+1);
			std::cout << Info::plain() << ".";

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
	#pragma omp parallel for
	for (esint d = 0; d < domains_in_global_index.size(); d++ ) {
			if (domains[d].B1Kplus.isOnACC == 1 || !configuration.combine_sc_and_spds) {
				// Calculates SC on CPU and keeps it CPU memory
				GetSchurComplement(USE_FLOAT, d);
				std::cout << ".";
			}
		}
#endif

	std::cout << "\n";

	CreateCudaStreamPool();

	// SC transfer to GPU - now sequential
	// TODO_GPU - Faze 1: tohle je zpusob jak se budou kopirovat i Dirichlet predpodminovace na GPU - inspirovat se, pripadne se muze dat rovnou to teto smycky 
	// TODO_GPU - Faze 2: zjistit, jak je to pomale a pripadne optimalizovat. Ale vyresi se implementaci vypoctu LSC na GPU 
	for (esint d = 0; d < domains_in_global_index.size(); d++ ) {

		esint status = 0;

		if (domains[d].B1Kplus.isOnACC == 1) {

			// TODO_GPU: Test performance (same streams)

			// Assign CUDA stream from he pool
			domains[d].B1Kplus.SetCUDA_Stream(cuda_stream_pool[d % STREAM_NUM]);
                        //if(USE_PREC == FETIConfiguration::PRECONDITIONER::DIRICHLET)
//                        {
//                                domains[d].Prec.SetCUDA_Stream(cuda_stream_pool[d % STREAM_NUM]);
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

				status_c = cudaMallocHost((void**)&domains[d].cuda_pinned_buff_fl, domains[d].B1_comp_dom.rows * sizeof(float));
				if (status_c != cudaSuccess) {
					std::cout << "Error allocating pinned host memory";
					status = -1;
				}
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

				status_c = cudaMallocHost((void**)&domains[d].cuda_pinned_buff, domains[d].B1_comp_dom.rows * sizeof(double));
				if (status_c != cudaSuccess) {
					std::cout << "Error allocating pinned host memory";
					status = -1;
				}
			}
#else
			if (USE_FLOAT){
                                status      = domains[d].B1Kplus.CopyToCUDA_Dev_fl();

                                size_t memsize = (domains[d].B1Kplus.rows > domains[d].Prec.rows) ? domains[d].B1Kplus.rows : domains[d].Prec.rows;
				// TODO_GPU - alokuji se pole pro vektory, kterymi se bude ve funkci Apply_A nasobit tato matice 
                                domains[d].cuda_pinned_buff_fl.resize(memsize);
                                cudaHostRegister(domains[d].cuda_pinned_buff_fl.data(), memsize * sizeof(float), cudaHostRegisterDefault);

			} else {

				status      = domains[d].B1Kplus.CopyToCUDA_Dev();

                                size_t memsize = (domains[d].B1Kplus.rows > domains[d].Prec.rows) ? domains[d].B1Kplus.rows : domains[d].Prec.rows;

                                domains[d].cuda_pinned_buff.resize(memsize);
                                cudaHostRegister(domains[d].cuda_pinned_buff.data(), memsize * sizeof(double), cudaHostRegisterDefault);

			}
#endif

			if (status == 0) {
				// TODO_GPU - LSCs ktere se uspesne prenesly do pameti GPU se smazou z pameti procesoru 
				domains[d].Kplus.keep_factors = false;

				if (USE_FLOAT) {
					SEQ_VECTOR <float>  ().swap (domains[d].B1Kplus.dense_values_fl);

					std::cout << "g";
				} else {
					SEQ_VECTOR <double> ().swap (domains[d].B1Kplus.dense_values);

					std::cout << "G";
				}
			} else {
				// pokud se domenu nepodar nahrat na GPU 


				domains[d].B1Kplus.isOnACC = 0;
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
					domains[d+1].isOnACC = 0;
					domains_on_CPU++;
					domains_on_GPU--;
					DOFs_CPU += domains[d].K.rows;
					DOFs_GPU -= domains[d].K.rows;
				}
#else


				if(domains[d].B1Kplus.USE_FLOAT) {
					domains[d].B1Kplus.FreeFromCUDA_Dev_fl();
				} else {
					domains[d].B1Kplus.FreeFromCUDA_Dev();
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

					if (USE_FLOAT)
						std::cout << "f";
					else
						std::cout << "F";

				} else {

					if (USE_FLOAT)
						std::cout << "c";
					else
						std::cout << "C";
				}
			}

		} else {

			if (configuration.combine_sc_and_spds) {

				if (USE_FLOAT)
					std::cout << "f";
				else
					std::cout << "F";

			} else {

				if (USE_FLOAT)
					std::cout << "c";
				else
					std::cout << "C";

			}

		}

		std::cout << " Domain: " << d << " GPU : " << domains[d].B1Kplus.isOnACC << "\n";
	}

	std::cout << "\n Domains transfered to GPU : " << domains_on_GPU << "\n";
	std::cout << " Domains on CPU : " << domains_on_CPU << "\n";

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
//					status_c = cudaMallocHost((void**)&domains[i].cuda_pinned_buff_fl, domains[i].B1_comp_dom.rows * sizeof(float));
//					if (status_c != cudaSuccess) {
//						////ESINFO(ERROR) << "Error allocating pinned host memory";
//						status = 1;
//					}
//				} else {
//					status_c = cudaMallocHost((void**)&domains[i].cuda_pinned_buff, domains[i].B1_comp_dom.rows * sizeof(double));
//					if (status_c != cudaSuccess) {
//						////ESINFO(ERROR) << "Error allocating pinned host memory";
//						status = 1;
//					}
//				}
//			} else {
//				status = 1;
//			}
//
//			// if status == 0 - all buffers in GPU mem were sucesfuly allocated
//			if (status == 0) {
//				domains[i].isOnACC = 1;
//				SEQ_VECTOR <double> ().swap (domains[i].B1Kplus.dense_values);
//				SEQ_VECTOR <float>  ().swap (domains[i].B1Kplus.dense_values_fl);
//				domains[i].Kplus.keep_factors = false;
//				if (USE_FLOAT)
//					////ESINFO(PROGRESS3) << Info::plain() << "g";
//				else
//					////ESINFO(PROGRESS3) << Info::plain() << "G";
//			} else {
//				domains[i].isOnACC = 0;
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
//                        domains[i].isOnACC = 0;
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

void ClusterGPU::GetSchurComplement( bool USE_FLOAT, esint i ) {

	SparseMatrix TmpB;
	domains[i].B1_comp_dom.MatTranspose(TmpB);

	SparseSolverCPU tmpsps;
//	if ( i == 0 && cluster_global_index == 1) {
//		tmpsps.msglvl = Info::report(LIBRARIES) ? 1 : 0;
//	}

	if (domains[i].K.type =='S')
	{
		tmpsps.Create_SC_w_Mat        ( domains[i].K, TmpB, domains[i].B1Kplus, false, 0 ); // general
	} else {
		if (domains[i].K.type =='G')
		{
			tmpsps.Create_non_sym_SC_w_Mat( domains[i].K, TmpB, TmpB, domains[i].B1Kplus, false, 0 );
		}
		else
		{
			std::cout << "Error - not defined type of K mat type.";
			exit(0);
		}
	}

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

void ClusterGPU::CreateDirichletPrec( DataHolder *instance) {

	std::cout << "Creating Dirichlet Preconditioner with Pardiso SC and coping them to GPU";

	esint status = 0;
	cudaError_t status_c;

	std::vector<size_t> local_Prec_size_to_add(domains_in_global_index.size(), 0);

	esint domains_on_GPU = 0;
	esint domains_on_CPU = 0;
	esint DOFs_GPU = 0;
	esint DOFs_CPU = 0;


	// Smycka napocitava velikost LSCs pres vsechny domeny
	for (esint d = 0; d < domains_in_global_index.size(); d++ ) {

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
		if (domains[d].B1Kplus.isOnACC == 0)
		{
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


		if(local_Prec_size_to_add[d] < GPU_free_mem)
		{
		      domains_on_GPU++;
		      DOFs_GPU += domains[d].Prec.rows;
		      domains[d].Prec.isOnACC = 1;
		      GPU_free_mem -= local_Prec_size_to_add[d];
		}
		else
		{
		      domains_on_CPU++;
		      DOFs_CPU += domains[d].Prec.rows;
		      domains[d].Prec.isOnACC = 0;
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
		      domains[d].isOnACC = 1;
		}
		else
		{
		      domains_on_CPU++;
		      DOFs_CPU += domains[d].K.rows;
		      domains[d].isOnACC = 0;
		}
	*/
	}

	// TODO_GPU - vsechny tyto std::cout se musi prepsat na logovani co ma Ondra M.
	// Ondro nektere moje rutiny, napr. SpyText jsou napsane pro std::cout a ne printf. Jake je reseni ?
	std::cout << "\n Domains on GPU : " << domains_on_GPU << "\n";
	std::cout << " Domains on CPU : " << domains_on_CPU << "\n";

	std::vector <int> on_gpu (info::mpi::size, 0);
	MPI_Gather(&domains_on_GPU,1,MPI_INT,&on_gpu[0],1,MPI_INT, 0, info::mpi::comm);

	std::vector <int> on_cpu (info::mpi::size, 0);
	MPI_Gather(&domains_on_CPU,1,MPI_INT,&on_cpu[0],1,MPI_INT, 0, info::mpi::comm);

	std::vector <int> don_gpu (info::mpi::size, 0);
	MPI_Gather(&DOFs_GPU,1,MPI_INT,&don_gpu[0],1,MPI_INT, 0, info::mpi::comm);

	std::vector <int> don_cpu (info::mpi::size, 0);
	MPI_Gather(&DOFs_CPU,1,MPI_INT,&don_cpu[0],1,MPI_INT, 0, info::mpi::comm);


	for (esint i = 0; i < info::mpi::size; i++) {
	      std::cout << " MPI rank " << i <<
		      "\t - GPU : domains = \t" << on_gpu[i] << "\t Total DOFs = \t" << don_gpu[i] <<
		      "\t - CPU : domains = \t" << on_cpu[i] << "\t Total DOFs = \t" << don_cpu[i] << "\n";
	}

	#pragma omp parallel for
	for (esint d = 0; d < domains_in_global_index.size(); d++ ) {
		// Calculates Prec on CPU and keeps it CPU memory
		GetDirichletPrec(instance, d);
		std::cout << ".";
	}
	std::cout << "\n";

	for (esint d = 0; d < domains_in_global_index.size(); d++ ) {

		esint status = 0;

		if (domains[d].Prec.isOnACC == 1) {

			// TODO_GPU: Test performance (same streams)

			// Assign CUDA stream from he pool
			domains[d].Prec.SetCUDA_Stream(cuda_stream_pool[d % STREAM_NUM]);


			if (domains[d].Prec.USE_FLOAT){

				status = domains[d].Prec.CopyToCUDA_Dev_fl();
				// B1Kplus is not on GPU, threfore vectors arent alocated
				if (domains[d].B1Kplus.isOnACC == 0)
				{
					size_t memsize = domains[d].Prec.rows;
					// TODO_GPU - alokuji se pole pro vektory, kterymi se bude ve funkci Apply_Prec nasobit tato matice
					domains[d].cuda_pinned_buff_fl.resize(memsize);
					cudaHostRegister(domains[d].cuda_pinned_buff_fl.data(), memsize * sizeof(float), cudaHostRegisterDefault);
				}

			} else {

				status = domains[d].Prec.CopyToCUDA_Dev();
				// B1Kplus is not on GPU, threfore vectors arent alocated
				if (domains[d].B1Kplus.isOnACC == 0)
				{
					size_t memsize = domains[d].Prec.rows;

					domains[d].cuda_pinned_buff.resize(memsize);
					cudaHostRegister(domains[d].cuda_pinned_buff.data(), memsize * sizeof(double), cudaHostRegisterDefault);
				}

			}

			if (status == 0) {
				// TODO_GPU - LSCs ktere se uspesne prenesly do pameti GPU se smazou z pameti procesoru

				if (domains[d].Prec.USE_FLOAT) {
					SEQ_VECTOR <float>  ().swap (domains[d].Prec.dense_values_fl);

					std::cout << "g";
				} else {
					SEQ_VECTOR <double> ().swap (domains[d].Prec.dense_values);

					std::cout << "G";
				}
			} else {
				// pokud se domenu nepodar nahrat na GPU

				domains[d].Prec.isOnACC = 0;
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
					if (domains[d].B1Kplus.isOnACC == 0)
					{
						cudaHostUnregister(domains[d].cuda_pinned_buff_fl.data());
						SEQ_VECTOR <float>  ().swap (domains[d].cuda_pinned_buff_fl);
					}

				} else {
					if (domains[d].B1Kplus.isOnACC == 0)
					{
						cudaHostUnregister(domains[d].cuda_pinned_buff.data());
						SEQ_VECTOR <double>  ().swap (domains[d].cuda_pinned_buff);
					}
				}



				if (domains[d].Prec.USE_FLOAT)
					std::cout << "c";
				else
					std::cout << "C";

			}

		} else {

			if (domains[d].Prec.USE_FLOAT)
				std::cout << "c";
			else
				std::cout << "C";

		}

		std::cout << " Domain: " << d << " GPU : " << domains[d].Prec.isOnACC << "\n";
	}

	std::cout << "\n Domains transfered to GPU : " << domains_on_GPU << "\n";
	std::cout << " Domains on CPU : " << domains_on_CPU << "\n";

}

void ClusterGPU::GetDirichletPrec( DataHolder *instance, esint d) {

//for (size_t d = 0; d < instance->K.size(); d++) {

	SEQ_VECTOR<esint> perm_vec = domains[d].B1t_Dir_perm_vec;
	SEQ_VECTOR<esint> perm_vec_full(instance->K[domains[d].domain_global_index].rows);// (instance->K[d].rows);
	SEQ_VECTOR<esint> perm_vec_diff(instance->K[domains[d].domain_global_index].rows);// (instance->K[d].rows);

	SEQ_VECTOR<esint> I_row_indices_p(instance->K[domains[d].domain_global_index].nnz);// (instance->K[d].nnz);
	SEQ_VECTOR<esint> J_col_indices_p(instance->K[domains[d].domain_global_index].nnz);// (instance->K[d].nnz);

	for (size_t i = 0; i < perm_vec.size(); i++) {
		perm_vec[i] = perm_vec[i] - 1;
	}

	for (size_t i = 0; i < perm_vec_full.size(); i++) {
		perm_vec_full[i] = i;
	}

	auto it = std::set_difference(perm_vec_full.begin(), perm_vec_full.end(), perm_vec.begin(), perm_vec.end(), perm_vec_diff.begin());
	perm_vec_diff.resize(it - perm_vec_diff.begin());

	perm_vec_full = perm_vec_diff;
	perm_vec_full.insert(perm_vec_full.end(), perm_vec.begin(), perm_vec.end());

	SparseMatrix K_modif = instance->K[domains[d].domain_global_index]; //[d];
#ifdef BEM4I_TO_BE_REMOVED
//TODO: Alex - da se nejak spocist Dir prec z dense matic K ??
	K_modif.ConvertDenseToCSR(1);
#endif
	SparseMatrix RegMatCRS = instance->RegMat[domains[d].domain_global_index]; //[d];
	RegMatCRS.ConvertToCSRwithSort(0);
	K_modif.MatAddInPlace(RegMatCRS, 'N', -1);
	// K_modif.RemoveLower();

	SEQ_VECTOR<SEQ_VECTOR<esint >> vec_I1_i2(K_modif.rows, SEQ_VECTOR<esint >(2, 1));
	esint offset = K_modif.CSR_I_row_indices[0] ? 1 : 0;

	for (esint i = 0; i < K_modif.rows; i++) {
		vec_I1_i2[i][0] = perm_vec_full[i];
		vec_I1_i2[i][1] = i; // position to create reverse permutation
	}

	std::sort(vec_I1_i2.begin(), vec_I1_i2.end(), [](const SEQ_VECTOR <esint >& a, const SEQ_VECTOR<esint>& b) {return a[0] < b[0];});

	// permutations made on matrix in COO format
	K_modif.ConvertToCOO(0);
	esint I_index, J_index;
	bool unsymmetric = !SYMMETRIC_SYSTEM;
	for (esint i = 0; i < K_modif.nnz; i++) {
		I_index = vec_I1_i2[K_modif.I_row_indices[i] - offset][1] + offset;
		J_index = vec_I1_i2[K_modif.J_col_indices[i] - offset][1] + offset;
		if (unsymmetric || I_index <= J_index) {
			I_row_indices_p[i] = I_index;
			J_col_indices_p[i] = J_index;
		} else {
			I_row_indices_p[i] = J_index;
			J_col_indices_p[i] = I_index;
		}
	}
	for (esint i = 0; i < K_modif.nnz; i++) {
		K_modif.I_row_indices[i] = I_row_indices_p[i];
		K_modif.J_col_indices[i] = J_col_indices_p[i];
	}
	K_modif.ConvertToCSRwithSort(1);
//	{
//		if (info::ecf->output.print_matrices) {
//			std::ofstream osS(Logging::prepareFile(d, "K_modif"));
//			osS << K_modif;
//			osS.close();
//		}
//	}

	// ------------------------------------------------------------------------------------------------------------------
	bool diagonalized_K_rr = configuration.preconditioner == FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET;
	//        PRECONDITIONER==NONE              - 0
	//        PRECONDITIONER==LUMPED            - 1
	//        PRECONDITIONER==WEIGHT_FUNCTION   - 2
	//        PRECONDITIONER==DIRICHLET         - 3
	//        PRECONDITIONER==SUPER_DIRICHLET   - 4
	//
	//        When next line is uncomment, var. PRECONDITIONER==DIRICHLET and PRECONDITIONER==SUPER_DIRICHLET provide identical preconditioner.
	//        bool diagonalized_K_rr = false
	// ------------------------------------------------------------------------------------------------------------------

	esint sc_size = perm_vec.size();

	if (sc_size == instance->K[domains[d].domain_global_index].rows) {
		domains[d].Prec = instance->K[domains[d].domain_global_index];
		domains[d].Prec.ConvertCSRToDense(1);
		// if physics.K[d] does not contain inner DOF
	} else {

		if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::DIRICHLET) {
			SparseSolverCPU createSchur;
//          createSchur.msglvl=1;
			esint sc_size = perm_vec.size();
			createSchur.ImportMatrix_wo_Copy(K_modif);
			createSchur.Create_SC(domains[d].Prec, sc_size, false);
			domains[d].Prec.ConvertCSRToDense(1);
		} else {
			SparseMatrix K_rr;
			SparseMatrix K_rs;
			SparseMatrix K_sr;
			SparseMatrix KsrInvKrrKrs;

			esint i_start = 0;
			esint nonsing_size = K_modif.rows - sc_size - i_start;
			esint j_start = nonsing_size;

			K_rs.getSubBlockmatrix_rs(K_modif, K_rs, i_start, nonsing_size, j_start, sc_size);

			if (SYMMETRIC_SYSTEM) {
				K_rs.MatTranspose(K_sr);
			} else {
				K_sr.getSubBlockmatrix_rs(K_modif, K_sr, j_start, sc_size, i_start, nonsing_size);
			}

			domains[d].Prec.getSubDiagBlockmatrix(K_modif, domains[d].Prec, nonsing_size, sc_size);
			SEQ_VECTOR<double> diagonals;
			SparseSolverCPU K_rr_solver;

			// K_rs is replaced by:
			// a) K_rs = 1/diag(K_rr) * K_rs          (simplified Dirichlet precond.)
			// b) K_rs =    inv(K_rr) * K_rs          (classical Dirichlet precond. assembled by own - not via PardisoSC routine)
			if (diagonalized_K_rr) {
				diagonals = K_modif.getDiagonal();
				// diagonals is obtained directly from K_modif (not from K_rr to avoid assembling) thanks to its structure
				//      K_modif = [K_rr, K_rs]
				//                [K_sr, K_ss]
				//
				for (esint i = 0; i < K_rs.rows; i++) {
					for (esint j = K_rs.CSR_I_row_indices[i]; j < K_rs.CSR_I_row_indices[i + 1]; j++) {
						K_rs.CSR_V_values[j - offset] /= diagonals[i];
					}
				}
			} else {
				K_rr.getSubDiagBlockmatrix(K_modif, K_rr, i_start, nonsing_size);
				K_rr_solver.ImportMatrix_wo_Copy(K_rr);
//            K_rr_solver.msglvl = 1;
				K_rr_solver.SolveMat_Dense(K_rs);
			}

			KsrInvKrrKrs.MatMat(K_sr, 'N', K_rs);
			domains[d].Prec.MatAddInPlace(KsrInvKrrKrs, 'N', -1);
//          if (!diagonalized_K_rr){
//				    cluster.domains[d].Prec.ConvertCSRToDense(1);
//          }
		}

	}

//	if (info::ecf->output.print_matrices) {
//		std::ofstream osS(Logging::prepareFile(domains[d].domain_global_index, "S"));
//		SparseMatrix SC = domains[d].Prec;
//		if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::DIRICHLET) {
//			SC.ConvertDenseToCSR(1);
//		}
//		osS << SC;
//		osS.close();
//	}

//	//ESINFO(PROGRESS3) << Info::plain() << ".";

}

void ClusterGPU::SetupKsolvers ( ) {

	#pragma omp parallel for
for (esint d = 0; d < domains.size(); d++) {

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
                                if ( domains[d].isOnACC == 0 ) {
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
for (esint d = 0; d < domains.size(); d++)
	{
		domains[d].B0Kplus_comp.DenseMatVec(x_in[d], tm2[d]);			// g0 - with comp B0Kplus
		domains[d].Kplus_R.DenseMatVec(x_in[d], tm3[d], 'T');			// e0
	}
	loop_1_1_time.end();

	loop_1_2_time.start();
	#pragma omp parallel for
for (esint d = 0; d < domains.size(); d++)
	{
		esint e0_start	=  d	* domains[d].Kplus_R.cols;
		esint e0_end		= (d+1) * domains[d].Kplus_R.cols;

		for (esint i = e0_start; i < e0_end; i++ )
			vec_e0[i] = - tm3[d][i - e0_start];
	}


	for (esint d = 0; d < domains.size(); d++)
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
for (esint i = 0; i < vec_e0.size(); i++)
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
for (esint i = 0; i < vec_g0.size(); i++)
		tm1[0][i] = vec_g0[i] - tm1[0][i];


	clus_F0_2_time.start();
	F0.Solve(tm1[0], vec_lambda,0,0);
	clus_F0_2_time.end();

	clusCP_time.end();


	// Kplus_x
	mkl_set_num_threads(1);
	loop_2_1_time.start();
	#pragma omp parallel for
for (esint d = 0; d < domains.size(); d++)
	{
		esint domain_size = domains[d].domain_prim_size;
		SEQ_VECTOR < double > tmp_vec (domains[d].B0_comp_map_vec.size(), 0.0);

		bool MIXED_SC_FACT = configuration.combine_sc_and_spds;

                if (domains[d].isOnACC == 0 && MIXED_SC_FACT) {

			for (esint i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
				tmp_vec[i] = vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;

			domains[d].B0_comp.MatVec(tmp_vec, tm1[d], 'T');

			for (esint i = 0; i < domain_size; i++)
				tm1[d][i] = x_in[d][i] - tm1[d][i];

			domains[d].multKplusLocal(tm1[d] , tm2[d]);

			esint e0_start	=  d	* domains[d].Kplus_R.cols;
			esint e0_end		= (d+1) * domains[d].Kplus_R.cols;

			domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);

		} else {

			for (esint i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
				tmp_vec[i] = -1.0 * vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;

			domains[d].B0Kplus_comp.DenseMatVec(tmp_vec, tm2[d], 'T' );

			esint e0_start	=  d	* domains[d].Kplus_R.cols;
			esint e0_end		= (d+1) * domains[d].Kplus_R.cols;
			domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);

		}

		for (esint i = 0; i < domain_size; i++)
			x_in[d][i] = tm2[d][i] + tm3[d][i];

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
for (esint d = 0; d < domains.size(); d++)
	{

		domains[d].B0Kplus_comp.DenseMatVec(x_in[d], tm2[d]);			// g0 - with comp B0Kplus
		domains[d].Kplus_R.     DenseMatVec(x_in[d], tm3[d], 'T');	    // e0

		esint e0_start	=  d	* domains[d].Kplus_R.cols;
		esint e0_end		= (d+1) * domains[d].Kplus_R.cols;

		for (esint i = e0_start; i < e0_end; i++ )
			vec_e0[i] = - tm3[d][i - e0_start];
	}

	for (esint d = 0; d < domains.size(); d++)
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
for (esint i = 0; i < vec_e0.size(); i++)
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
for (esint i = 0; i < vec_g0.size(); i++)
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
for (esint d = 0; d < domains.size(); d++)
	{
		esint domain_size = domains[d].domain_prim_size;
		SEQ_VECTOR < double > tmp_vec (domains[d].B0_comp_map_vec.size(), 0.0);

		for (esint i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
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
for (esint d = 0; d < domains.size(); d++)
	{
		esint domain_size = domains[d].domain_prim_size;
		SEQ_VECTOR < double > tmp_vec (domains[d].B0_comp_map_vec.size(), 0.0);

		for (esint i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
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
for (esint d = 0; d < domains.size(); d++)
	{
		esint domain_size = domains[d].domain_prim_size;
		SEQ_VECTOR < double > tmp_vec (domains[d].B0_comp_map_vec.size(), 0.0);

		bool MIXED_SC_FACT = configuration.combine_sc_and_spds;

                if ( (domains[d].isOnACC == 0 && MIXED_SC_FACT ) || !configuration.use_schur_complement ) {

			for (esint i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
				tmp_vec[i] = vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;

			domains[d].B0_comp.MatVec(tmp_vec, tm1[d], 'T');

			for (esint i = 0; i < domain_size; i++)
				tm1[d][i] = x_in[d][i] - tm1[d][i];

			domains[d].multKplusLocal(tm1[d] , tm2[d]);

			esint e0_start	=  d	* domains[d].Kplus_R.cols;
			esint e0_end		= (d+1) * domains[d].Kplus_R.cols;

			domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);

		} else {

			for (esint i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
				tmp_vec[i] = -1.0 * vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;

			domains[d].B0Kplus_comp.DenseMatVec(tmp_vec, tm2[d], 'T' );

			esint e0_start	=  d	* domains[d].Kplus_R.cols;
			esint e0_end		= (d+1) * domains[d].Kplus_R.cols;

			domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);

		}

		for (esint i = 0; i < domain_size; i++)
			x_in[d][i] = tm2[d][i] + tm3[d][i];

	}
	 loop_2_1_time.end();

}



void ClusterGPU::CreateCudaStreamPool() {
#ifdef STREAM_NUM
	if (cuda_stream_pool.size() < STREAM_NUM) {

		DestroyCudaStreamPool();

		cuda_stream_pool.resize(STREAM_NUM);

		for(esint i = 0; i < STREAM_NUM; i++) {
			cudaStreamCreate(&cuda_stream_pool[i]);
		}
		////ESINFO(VERBOSE_LEVEL3) << "CUDA stream pool created with " << STREAM_NUM << " streams per cluster";
	}
#endif
}

void ClusterGPU::DestroyCudaStreamPool() {
	for(esint i = 0; i < cuda_stream_pool.size(); i++) {
		cudaStreamDestroy(cuda_stream_pool[i]);
	}

	SEQ_VECTOR <cudaStream_t>().swap(cuda_stream_pool);
}
