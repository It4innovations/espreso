#include "clusterGPU.h"

using namespace espreso;

//#define STREAM_NUM 16
//#define SHARE_SC

ClusterGPU::~ClusterGPU() {
	DestroyCudaStreamPool();

	for(eslocal d = 0; d < domains_in_global_index.size(); d++) {
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

void ClusterGPU::Create_SC_perDomain(bool USE_FLOAT) {

	ESINFO(PROGRESS2) << "Creating B1*K+*B1t Schur Complements with Pardiso SC and coping them to GPU";

	bool GPU_full = false;
	//GPU_full = true;

	eslocal status = 0;
	cudaError_t status_c;

	int nDevices;
	cudaGetDeviceCount(&nDevices);

	ESINFO(PROGRESS2) << Info::plain() << "\n*** GPU Accelerators available on the server " << "\n\n";
	for (int i = 0; i < nDevices; i++) {
		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, i);
		ESINFO(PROGRESS2) << Info::plain() << " Device Number: " << i << "\n";
		ESINFO(PROGRESS2) << Info::plain() << " Device name: " << prop.name << "\n";
		ESINFO(PROGRESS2) << Info::plain() << " Memory Clock Rate (KHz): " << prop.memoryClockRate << "\n";
		ESINFO(PROGRESS2) << Info::plain() << " Memory Bus Width (bits): " << prop.memoryBusWidth << "\n";
		ESINFO(PROGRESS2) << Info::plain() << " Peak Memory Bandwidth (GB/s): " << 2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6 << "\n";

		cudaSetDevice(i);
		size_t free, total;
		cudaMemGetInfo(&free, &total);
		ESINFO(PROGRESS2) << Info::plain() << " GPU Total Memory [MB]: " << total/1024/1024 << "\n";
		ESINFO(PROGRESS2) << Info::plain() << " GPU Free Memory [MB]:  " << free/1024/1024 << "\n\n";

	}

	// GPU memory management
	#if DEVICE_ID == 1
		ESINFO(VERBOSE_LEVEL3) << "Selected CUDA device 1";
		cudaSetDevice(1);
	#else
		cudaSetDevice(0);
	#endif

	size_t GPU_free_mem, GPU_total_meml;
	cudaMemGetInfo(&GPU_free_mem, &GPU_total_meml);

	eslocal domains_on_GPU = 0;
	eslocal domains_on_CPU = 0;
	eslocal DOFs_GPU = 0;
	eslocal DOFs_CPU = 0;

	size_t SC_total_size = 0;
	for (eslocal d = 0; d < domains_in_global_index.size(); d++ ) {

		switch (config::solver::SCHUR_COMPLEMENT_TYPE) {
		case config::solver::SCHUR_COMPLEMENT_TYPEalternative::GENERAL:
#ifdef SHARE_SC
			// SC_total_size will be halved in the case of 2 symmetric SCs in 1 full matrix
			if (d%2 == 0) {
				eslocal sc1_rows = domains[d].B1_comp_dom.rows;
				eslocal sc2_rows = 0;
				eslocal SC_size = 0;
				eslocal vec_size = domains[d].B1_comp_dom.rows;

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
					SC_total_size += (SC_size + 2 * vec_size) * sizeof(float);
				} else {
					SC_total_size += (SC_size + 2 * vec_size) * sizeof(double);
				}
			}
#else
			if (USE_FLOAT) {
				SC_total_size +=
						( domains[d].B1_comp_dom.rows * domains[d].B1_comp_dom.rows +
						  2 * domains[d].B1_comp_dom.rows
						) * sizeof(float);
			} else {
				SC_total_size +=
						( domains[d].B1_comp_dom.rows * domains[d].B1_comp_dom.rows +
						  2 * domains[d].B1_comp_dom.rows
						) * sizeof(double);
			}
#endif
			break;
		case config::solver::SCHUR_COMPLEMENT_TYPEalternative::SYMMETRIC:
			if (USE_FLOAT) {
				SC_total_size +=
						(((domains[d].B1_comp_dom.rows + 1 ) * domains[d].B1_comp_dom.rows ) / 2
						 + 2 * domains[d].B1_comp_dom.rows
						) * sizeof(float);
			} else {
				SC_total_size +=
						(((domains[d].B1_comp_dom.rows + 1 ) * domains[d].B1_comp_dom.rows ) / 2
						 + 2 * domains[d].B1_comp_dom.rows
						) * sizeof(double);
			}
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented type of Schur complement.";
		}


		if (SC_total_size < GPU_free_mem) {
			domains_on_GPU++;
			DOFs_GPU += domains[d].K.rows;
			domains[d].isOnACC = 1;
		} else {
			domains_on_CPU++;
			DOFs_CPU += domains[d].K.rows;
			domains[d].isOnACC = 0;
		}
	}

	ESINFO(PROGRESS2) << Info::plain() << "\n Domains on GPU : " << domains_on_GPU << "\n";
	ESINFO(PROGRESS2) << Info::plain() << " Domains on CPU : " << domains_on_CPU << "\n";

	std::vector <int> on_gpu (config::env::MPIsize, 0);
	MPI_Gather(&domains_on_GPU,1,MPI_INT,&on_gpu[0],1,MPI_INT, 0, MPI_COMM_WORLD);

	std::vector <int> on_cpu (config::env::MPIsize, 0);
	MPI_Gather(&domains_on_CPU,1,MPI_INT,&on_cpu[0],1,MPI_INT, 0, MPI_COMM_WORLD);

	std::vector <int> don_gpu (config::env::MPIsize, 0);
	MPI_Gather(&DOFs_GPU,1,MPI_INT,&don_gpu[0],1,MPI_INT, 0, MPI_COMM_WORLD);

	std::vector <int> don_cpu (config::env::MPIsize, 0);
	MPI_Gather(&DOFs_CPU,1,MPI_INT,&don_cpu[0],1,MPI_INT, 0, MPI_COMM_WORLD);


	for (eslocal i = 0; i < config::env::MPIsize; i++) {
		ESINFO(PROGRESS2) << Info::plain()
			<< " MPI rank " << i <<
			"\t - GPU : domains = \t" << on_gpu[i] << "\t Total DOFs = \t" << don_gpu[i] <<
			"\t - CPU : domains = \t" << on_cpu[i] << "\t Total DOFs = \t" << don_cpu[i] << "\n";
	}


	// Schur complement calculation on CPU
//	cilk_for (eslocal d = 0; d < domains_in_global_index.size(); d++ ) {
//		if (domains[d].isOnACC == 1 || !config::solver::COMBINE_SC_AND_SPDS) {
//			// Calculates SC on CPU and keeps it CPU memory
//			GetSchurComplement(USE_FLOAT, d);
//			ESINFO(PROGRESS2) << Info::plain() << ".";
//		}
//	}


#ifdef SHARE_SC
	SEQ_VECTOR <eslocal> SC_dense_val_offsets(domains_in_global_index.size(), 0);

	// 2 domains per iteration processed
	cilk_for (eslocal d = 0; d < domains_in_global_index.size(); d += 2 ) {

		if (domains[d].isOnACC == 1 || !config::solver::COMBINE_SC_AND_SPDS) {
			// Calculates SC on CPU and keeps it CPU memory
			GetSchurComplement(USE_FLOAT, d);
			ESINFO(PROGRESS2) << Info::plain() << ".";

			// Set if Upper or Lower part is referenced
			domains[d].B1Kplus.uplo = 'U';

			// Set the default lda
			domains[d].B1Kplus.extern_lda = domains[d].B1Kplus.rows;
		}

		if (d+1 < domains_in_global_index.size() && (domains[d+1].isOnACC == 1 || !config::solver::COMBINE_SC_AND_SPDS)) {
			// Calculates SC on CPU and keeps it CPU memory
			GetSchurComplement(USE_FLOAT, d+1);
			ESINFO(PROGRESS2) << Info::plain() << ".";

			eslocal sc1_rows = domains[d].B1Kplus.rows;
			eslocal sc2_rows = domains[d+1].B1Kplus.rows;

			// Set if Upper or Lower part is referenced
			domains[d+1].B1Kplus.uplo = 'L';

			// Both SCs stored in the first domain [d]
			if(sc1_rows > sc2_rows) {
				// First SC -> U
				if (USE_FLOAT) {
					for(eslocal r = 0; r < sc2_rows; r++) {
						std::copy(&domains[d+1].B1Kplus.dense_values_fl[r*sc2_rows + r], &domains[d+1].B1Kplus.dense_values_fl[(r+1) *sc2_rows],
						 &domains[d].B1Kplus.dense_values_fl[r*sc1_rows + r+1]);
					}
				} else {
					for(eslocal r = 0; r < sc2_rows; r++) {
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

					for(eslocal r = 0; r < sc1_rows; r++) {
						std::copy(&domains[d+1].B1Kplus.dense_values_fl[r*sc2_rows + r], &domains[d+1].B1Kplus.dense_values_fl[(r+1) *sc2_rows],
						 &sc1_tmp_fl[r*sc1_rows + r]);
						std::copy(&domains[d].B1Kplus.dense_values_fl[r*sc1_rows], &domains[d].B1Kplus.dense_values_fl[r*sc1_rows + r+1],
						 &sc1_tmp_fl[(r+1)*sc1_rows]);
					}
					domains[d].B1Kplus.dense_values_fl.swap(sc1_tmp_fl);
				} else {
					SEQ_VECTOR <double> sc1_tmp(sc1_rows * (sc1_rows + 1));

					for(eslocal r = 0; r < sc1_rows; r++) {
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
					for(eslocal r = 0; r < sc1_rows; r++) {
						std::copy(&domains[d].B1Kplus.dense_values_fl[r*sc1_rows], &domains[d].B1Kplus.dense_values_fl[r*sc1_rows + r + 1],
						 &domains[d+1].B1Kplus.dense_values_fl[(r+1)*sc2_rows]);
					}
					domains[d].B1Kplus.dense_values_fl = std::move(domains[d+1].B1Kplus.dense_values_fl);
				} else {
					for(eslocal r = 0; r < sc1_rows; r++) {
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
	cilk_for (eslocal d = 0; d < domains_in_global_index.size(); d++ ) {
			if (domains[d].isOnACC == 1 || !config::solver::COMBINE_SC_AND_SPDS) {
				// Calculates SC on CPU and keeps it CPU memory
				GetSchurComplement(USE_FLOAT, d);
				ESINFO(PROGRESS2) << Info::plain() << ".";
			}
		}
#endif

	ESINFO(PROGRESS2) << Info::plain() << "\n";

	CreateCudaStreamPool();

	// SC transfer to GPU - now sequential
	for (eslocal d = 0; d < domains_in_global_index.size(); d++ ) {

		eslocal status = 0;

		if (domains[d].isOnACC == 1) {

			// Calculates SC on CPU and keeps it CPU memory
			//GetSchurComplement(USE_FLOAT, d);


			// TODO Test performance (same streams)
			// Assign CUDA stream from he pool
			domains[d].B1Kplus.SetCUDA_Stream(cuda_stream_pool[d % STREAM_NUM]);

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
					ESINFO(ERROR) << "Error allocating pinned host memory";
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
					ESINFO(ERROR) << "Error allocating pinned host memory";
					status = -1;
				}
			}
#else
			if (USE_FLOAT){
				status = domains[d].B1Kplus.CopyToCUDA_Dev_fl();

				status_c = cudaMallocHost((void**)&domains[d].cuda_pinned_buff_fl, domains[d].B1_comp_dom.rows * sizeof(float));
				if (status_c != cudaSuccess) {
					ESINFO(ERROR) << "Error allocating pinned host memory";
					status = -1;
				}
			} else {
				status = domains[d].B1Kplus.CopyToCUDA_Dev();

				status_c = cudaMallocHost((void**)&domains[d].cuda_pinned_buff, domains[d].B1_comp_dom.rows * sizeof(double));
				if (status_c != cudaSuccess) {
					ESINFO(ERROR) << "Error allocating pinned host memory";
					status = -1;
				}
			}
#endif

			if (status == 0) {
				domains[d].Kplus.keep_factors = false;

				if (USE_FLOAT) {
					SEQ_VECTOR <float>  ().swap (domains[d].B1Kplus.dense_values_fl);

					ESINFO(PROGRESS2) << Info::plain() << "g";
				} else {
					SEQ_VECTOR <double> ().swap (domains[d].B1Kplus.dense_values);

					ESINFO(PROGRESS2) << Info::plain() << "G";
				}
			} else {

				domains[d].isOnACC = 0;
				domains_on_CPU++;
				domains_on_GPU--;
				DOFs_CPU += domains[d].K.rows;
				DOFs_GPU -= domains[d].K.rows;

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
					cudaFreeHost(domains[d].cuda_pinned_buff_fl);
				} else {
					cudaFreeHost(domains[d].cuda_pinned_buff);
				}

				if (config::solver::COMBINE_SC_AND_SPDS) {

					SEQ_VECTOR <double> ().swap (domains[d].B1Kplus.dense_values);
					SEQ_VECTOR <float>  ().swap (domains[d].B1Kplus.dense_values_fl);

					if (USE_FLOAT)
						ESINFO(PROGRESS2) << Info::plain() << "f";
					else
						ESINFO(PROGRESS2) << Info::plain() << "F";

				} else {

					if (USE_FLOAT)
						ESINFO(PROGRESS2) << Info::plain() << "c";
					else
						ESINFO(PROGRESS2) << Info::plain() << "C";
				}
			}

		} else {

			if (config::solver::COMBINE_SC_AND_SPDS) {

				if (USE_FLOAT)
					ESINFO(PROGRESS2) << Info::plain() << "f";
				else
					ESINFO(PROGRESS2) << Info::plain() << "F";

			} else {

				//GetSchurComplement(USE_FLOAT, d);

				if (USE_FLOAT)
					ESINFO(PROGRESS2) << Info::plain() << "c";
				else
					ESINFO(PROGRESS2) << Info::plain() << "C";

			}

		}

//		ESINFO(PROGRESS2) << Info::plain() << " Domain: " << d << " GPU : " << domains[d].isOnACC << "\n";
	}

	ESINFO(PROGRESS2) << Info::plain() << "\n Domains transfered to GPU : " << domains_on_GPU << "\n";
	ESINFO(PROGRESS2) << Info::plain() << " Domains on CPU : " << domains_on_CPU << "\n";

//	cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ ) {
//
////		cudaSetDevice(1);
//
////		SparseSolverCPU tmpsps2;
////		if ( i == 0 && cluster_global_index == 1) tmpsps2.msglvl = 1;
////		tmpsps2.Create_non_sym_SC_w_Mat( domains[i].K, TmpB, domains[i].B0t_comp, domains[i].B0KplusB1_comp, false, 0 );
//
//		eslocal status = 0;
//		cudaError_t status_c;
//
//		if (!GPU_full || !config::solver::COMBINE_SC_AND_SPDS) {
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
//						ESINFO(ERROR) << "Error allocating pinned host memory";
//						status = 1;
//					}
//				} else {
//					status_c = cudaMallocHost((void**)&domains[i].cuda_pinned_buff, domains[i].B1_comp_dom.rows * sizeof(double));
//					if (status_c != cudaSuccess) {
//						ESINFO(ERROR) << "Error allocating pinned host memory";
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
//					ESINFO(PROGRESS2) << Info::plain() << "g";
//				else
//					ESINFO(PROGRESS2) << Info::plain() << "G";
//			} else {
//				domains[i].isOnACC = 0;
//				GPU_full = true;
//				if (config::solver::COMBINE_SC_AND_SPDS) {
//					SEQ_VECTOR <double> ().swap (domains[i].B1Kplus.dense_values);
//					SEQ_VECTOR <float>  ().swap (domains[i].B1Kplus.dense_values_fl);
//					if (USE_FLOAT)
//						ESINFO(PROGRESS2) << Info::plain() << "p";
//					else
//						ESINFO(PROGRESS2) << Info::plain() << "P";
//				} else {
//					if (USE_FLOAT)
//						ESINFO(PROGRESS2) << Info::plain() << "c";
//					else
//						ESINFO(PROGRESS2) << Info::plain() << "C";
//				}
//			}
//
//		} else {
//                        domains[i].isOnACC = 0;
//			if (USE_FLOAT)
//				ESINFO(PROGRESS2) << Info::plain() << "p";
//			else
//				ESINFO(PROGRESS2) << Info::plain() << "P";
//		}
//
//		//GPU_full = true;
//
//	}

	ESINFO(PROGRESS2);

}

void ClusterGPU::GetSchurComplement( bool USE_FLOAT, eslocal i ) {

	SparseMatrix TmpB;
	domains[i].B1_comp_dom.MatTranspose(TmpB);

	SparseSolverCPU tmpsps;
//	if ( i == 0 && cluster_global_index == 1) {
//		tmpsps.msglvl = Info::report(LIBRARIES) ? 1 : 0;
//	}

	tmpsps.Create_SC_w_Mat( domains[i].K, TmpB, domains[i].B1Kplus, false, 0 ); // general

	if (USE_FLOAT){
		domains[i].B1Kplus.ConvertDenseToDenseFloat( 1 );
		domains[i].B1Kplus.USE_FLOAT = true;
	}

	//ESINFO(PROGRESS2) << Info::plain() << "s";

	// if Schur complement is symmetric - then remove lower part - slower for GPU but more mem. efficient
	if (config::solver::SCHUR_COMPLEMENT_TYPE == config::solver::SCHUR_COMPLEMENT_TYPEalternative::SYMMETRIC) {
		domains[i].B1Kplus.RemoveLowerDense();
	}

}

void ClusterGPU::SetupKsolvers ( ) {

	cilk_for (eslocal d = 0; d < domains.size(); d++) {

		// Import of Regularized matrix K into Kplus (Sparse Solver)
		switch (config::solver::KSOLVER) {
		case config::solver::KSOLVERalternative::DIRECT_DP:
			domains[d].Kplus.ImportMatrix_wo_Copy (domains[d].K);
			break;
		case config::solver::KSOLVERalternative::ITERATIVE:
			domains[d].Kplus.ImportMatrix_wo_Copy (domains[d].K);
			break;
		case config::solver::KSOLVERalternative::DIRECT_SP:
			domains[d].Kplus.ImportMatrix_wo_Copy_fl(domains[d].K);
			//domains[d].Kplus.ImportMatrix_fl(domains[d].K);
			break;
		case config::solver::KSOLVERalternative::DIRECT_MP:
			domains[d].Kplus.ImportMatrix_fl(domains[d].K);
			break;
//		case 4:
//			domains[d].Kplus.ImportMatrix_fl(domains[d].K);
//			break;
		default:
			ESINFO(ERROR) << "Invalid KSOLVER value.";
			exit(EXIT_FAILURE);
		}

		if (config::solver::KEEP_FACTORS) {

			if (!config::solver::COMBINE_SC_AND_SPDS) { // if both CPU and GPU uses Schur Complement
				std::stringstream ss;
				ss << "init -> rank: " << config::env::MPIrank << ", subdomain: " << d;
				domains[d].Kplus.keep_factors = true;
				if (config::solver::KSOLVER != config::solver::KSOLVERalternative::ITERATIVE) {
					domains[d].Kplus.Factorization (ss.str());
				}
			} else {
				if ( domains[d].isOnACC == 0 ) {
					std::stringstream ss;
					ss << "init -> rank: " << config::env::MPIrank << ", subdomain: " << d;
					domains[d].Kplus.keep_factors = true;
					if (config::solver::KSOLVER != config::solver::KSOLVERalternative::ITERATIVE) {
						domains[d].Kplus.Factorization (ss.str());
					}
				}
			}

		} else {
			domains[d].Kplus.keep_factors = false;
			domains[d].Kplus.MPIrank = config::env::MPIrank;
		}

		domains[d].domain_prim_size = domains[d].Kplus.cols;

		if ( d == 0 && config::env::MPIrank == 0) {
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
	cilk_for (eslocal d = 0; d < domains.size(); d++)
	{
		domains[d].B0Kplus_comp.DenseMatVec(x_in[d], tm2[d]);			// g0 - with comp B0Kplus
		domains[d].Kplus_R.DenseMatVec(x_in[d], tm3[d], 'T');			// e0
	}
	loop_1_1_time.end();

	loop_1_2_time.start();
	cilk_for (eslocal d = 0; d < domains.size(); d++)
	{
		eslocal e0_start	=  d	* domains[d].Kplus_R.cols;
		eslocal e0_end		= (d+1) * domains[d].Kplus_R.cols;

		for (eslocal i = e0_start; i < e0_end; i++ )
			vec_e0[i] = - tm3[d][i - e0_start];
	}


	for (eslocal d = 0; d < domains.size(); d++)
		for (eslocal i = 0; i < domains[d].B0Kplus_comp.rows; i++)
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

	cilk_for (eslocal i = 0; i < vec_e0.size(); i++)
		tm2[0][i] = tm2[0][i] - vec_e0[i];
	//cblas_daxpy(vec_e0.size(), -1.0, &vec_e0[0], 1, &tm2[0][0], 1);

	 clus_Sa_time.start();
//#ifdef SPARSE_SA
//	 Sa.Solve(tm2[0], vec_alfa,0,0);
//#else
//	char U = 'U';
//	eslocal nrhs = 1;
//	eslocal info = 0;
//	vec_alfa = tm2[0];
//	dsptrs( &U, &SaMat.rows, &nrhs, &SaMat.dense_values[0], &SaMat.ipiv[0], &vec_alfa[0], &SaMat.rows, &info );
//#endif
//
////TODO: New version with DenseSolver routine
////#ifdef SPARSE_SA
////	Sa.Solve(tm2[0], vec_alfa,0,0);
////#else
////	eslocal nrhs = 1;
////	Sa_dense.Solve(tm2[0], vec_alfa, nrhs);
////#endif

	 if (config::solver::SASOLVER == config::solver::SASOLVERalternative::CPU_SPARSE) {
		 Sa.Solve(tm2[0], vec_alfa,0,0);
	 }

	 if (config::solver::SASOLVER == config::solver::SASOLVERalternative::CPU_DENSE) {
			eslocal nrhs = 1;
			Sa_dense_cpu.Solve(tm2[0], vec_alfa, nrhs);
	 }

	 if (config::solver::SASOLVER == config::solver::SASOLVERalternative::ACC_DENSE) {
			eslocal nrhs = 1;
			Sa_dense_acc.Solve(tm2[0], vec_alfa, nrhs);
	 }

	 clus_Sa_time.end();



	 clus_G0t_time.start();
	G0.MatVec(vec_alfa, tm1[0], 'T'); 	// lambda
	 clus_G0t_time.end();

	cilk_for (eslocal i = 0; i < vec_g0.size(); i++)
		tm1[0][i] = vec_g0[i] - tm1[0][i];


	clus_F0_2_time.start();
	F0.Solve(tm1[0], vec_lambda,0,0);
	clus_F0_2_time.end();

	clusCP_time.end();


	// Kplus_x
	mkl_set_num_threads(1);
	loop_2_1_time.start();
	cilk_for (eslocal d = 0; d < domains.size(); d++)
	{
		eslocal domain_size = domains[d].domain_prim_size;
		SEQ_VECTOR < double > tmp_vec (domains[d].B0_comp_map_vec.size(), 0.0);

		bool MIXED_SC_FACT = config::solver::COMBINE_SC_AND_SPDS;

		if (domains[d].isOnACC == 0 && MIXED_SC_FACT) {

			for (eslocal i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
				tmp_vec[i] = vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;

			domains[d].B0_comp.MatVec(tmp_vec, tm1[d], 'T');

			for (eslocal i = 0; i < domain_size; i++)
				tm1[d][i] = x_in[d][i] - tm1[d][i];

			domains[d].multKplusLocal(tm1[d] , tm2[d]);

			eslocal e0_start	=  d	* domains[d].Kplus_R.cols;
			eslocal e0_end		= (d+1) * domains[d].Kplus_R.cols;

			domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);

		} else {

			for (eslocal i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
				tmp_vec[i] = -1.0 * vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;

			domains[d].B0Kplus_comp.DenseMatVec(tmp_vec, tm2[d], 'T' );

			eslocal e0_start	=  d	* domains[d].Kplus_R.cols;
			eslocal e0_end		= (d+1) * domains[d].Kplus_R.cols;
			domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);

		}

		for (eslocal i = 0; i < domain_size; i++)
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
	cilk_for (eslocal d = 0; d < domains.size(); d++)
	{

		domains[d].B0Kplus_comp.DenseMatVec(x_in[d], tm2[d]);			// g0 - with comp B0Kplus
		domains[d].Kplus_R.     DenseMatVec(x_in[d], tm3[d], 'T');	    // e0

		eslocal e0_start	=  d	* domains[d].Kplus_R.cols;
		eslocal e0_end		= (d+1) * domains[d].Kplus_R.cols;

		for (eslocal i = e0_start; i < e0_end; i++ )
			vec_e0[i] = - tm3[d][i - e0_start];
	}

	for (eslocal d = 0; d < domains.size(); d++)
		for (eslocal i = 0; i < domains[d].B0Kplus_comp.rows; i++)
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

	cilk_for (eslocal i = 0; i < vec_e0.size(); i++)
		tm2[0][i] = tm2[0][i] - vec_e0[i];
	//cblas_daxpy(vec_e0.size(), -1.0, &vec_e0[0], 1, &tm2[0][0], 1);


	 clus_Sa_time.start();
	if (config::solver::SASOLVER == config::solver::SASOLVERalternative::CPU_SPARSE) {
		Sa.Solve(tm2[0], vec_alfa,0,0);
	}

	if (config::solver::SASOLVER == config::solver::SASOLVERalternative::CPU_DENSE) {
		eslocal nrhs = 1;
		Sa_dense_cpu.Solve(tm2[0], vec_alfa, nrhs);
	}

	if (config::solver::SASOLVER == config::solver::SASOLVERalternative::ACC_DENSE) {
		eslocal nrhs = 1;
		Sa_dense_acc.Solve(tm2[0], vec_alfa, nrhs);// lambda
	}
	 clus_Sa_time.end();


	 clus_G0t_time.start();
	G0.MatVec(vec_alfa, tm1[0], 'T'); 	// lambda
	 clus_G0t_time.end();

	cilk_for (eslocal i = 0; i < vec_g0.size(); i++)
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

	cilk_for (eslocal d = 0; d < domains.size(); d++)
	{
		eslocal domain_size = domains[d].domain_prim_size;
		SEQ_VECTOR < double > tmp_vec (domains[d].B0_comp_map_vec.size(), 0.0);

		for (eslocal i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
			tmp_vec[i] = -1.0 * vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;

		domains[d].B0Kplus_comp.DenseMatVec(tmp_vec, tm2[d], 'T' );

		eslocal e0_start	=  d	* domains[d].Kplus_R.cols;
		eslocal e0_end		= (d+1) * domains[d].Kplus_R.cols;

		domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);

		for (eslocal i = 0; i < domain_size; i++)
			y_out[d][i] = tm2[d][i] + tm3[d][i];

	}
	 loop_2_1_time.end();

}

void ClusterGPU::multKplus_HF_Loop2_SPDS (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {

	// Kplus_x
	//mkl_set_num_threads(1);
	 loop_2_1_time.start();

	cilk_for (eslocal d = 0; d < domains.size(); d++)
	{
		eslocal domain_size = domains[d].domain_prim_size;
		SEQ_VECTOR < double > tmp_vec (domains[d].B0_comp_map_vec.size(), 0.0);

		for (eslocal i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
			tmp_vec[i] = vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;

		domains[d].B0_comp.MatVec(tmp_vec, tm1[d], 'T');

		for (eslocal i = 0; i < domain_size; i++)
			tm1[d][i] = x_in[d][i] - tm1[d][i];

		domains[d].multKplusLocal(tm1[d] , tm2[d]);

		eslocal e0_start	=  d	* domains[d].Kplus_R.cols;
		eslocal e0_end		= (d+1) * domains[d].Kplus_R.cols;

		domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);


		for (eslocal i = 0; i < domain_size; i++)
			x_in[d][i] = tm2[d][i] + tm3[d][i];

	}
	 loop_2_1_time.end();

}

void ClusterGPU::multKplus_HF_Loop2_MIX (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {

	// Kplus_x
	// mkl_set_num_threads(1);
	 loop_2_1_time.start();

	cilk_for (eslocal d = 0; d < domains.size(); d++)
	{
		eslocal domain_size = domains[d].domain_prim_size;
		SEQ_VECTOR < double > tmp_vec (domains[d].B0_comp_map_vec.size(), 0.0);

		bool MIXED_SC_FACT = config::solver::COMBINE_SC_AND_SPDS;

		if ( (domains[d].isOnACC == 0 && MIXED_SC_FACT ) || !config::solver::USE_SCHUR_COMPLEMENT ) {

			for (eslocal i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
				tmp_vec[i] = vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;

			domains[d].B0_comp.MatVec(tmp_vec, tm1[d], 'T');

			for (eslocal i = 0; i < domain_size; i++)
				tm1[d][i] = x_in[d][i] - tm1[d][i];

			domains[d].multKplusLocal(tm1[d] , tm2[d]);

			eslocal e0_start	=  d	* domains[d].Kplus_R.cols;
			eslocal e0_end		= (d+1) * domains[d].Kplus_R.cols;

			domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);

		} else {

			for (eslocal i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
				tmp_vec[i] = -1.0 * vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;

			domains[d].B0Kplus_comp.DenseMatVec(tmp_vec, tm2[d], 'T' );

			eslocal e0_start	=  d	* domains[d].Kplus_R.cols;
			eslocal e0_end		= (d+1) * domains[d].Kplus_R.cols;

			domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);

		}

		for (eslocal i = 0; i < domain_size; i++)
			x_in[d][i] = tm2[d][i] + tm3[d][i];

	}
	 loop_2_1_time.end();

}

void ClusterGPU::CreateCudaStreamPool() {
#ifdef STREAM_NUM
	if (cuda_stream_pool.size() < STREAM_NUM) {

		DestroyCudaStreamPool();

		cuda_stream_pool.resize(STREAM_NUM);

		for(eslocal i = 0; i < STREAM_NUM; i++) {
			cudaStreamCreate(&cuda_stream_pool[i]);
		}
		ESINFO(VERBOSE_LEVEL3) << "CUDA stream pool created with " << STREAM_NUM << " streams per cluster";
	}
#endif
}

void ClusterGPU::DestroyCudaStreamPool() {
	for(eslocal i = 0; i < cuda_stream_pool.size(); i++) {
		cudaStreamDestroy(cuda_stream_pool[i]);
	}

	SEQ_VECTOR <cudaStream_t>().swap(cuda_stream_pool);
}
