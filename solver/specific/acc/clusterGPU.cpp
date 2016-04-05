#include "clusterGPU.h"

using namespace espreso;

void ClusterGPU::Create_SC_perDomain(bool USE_FLOAT) {

	ESINFO(PROGRESS2) << "Creating B1*K+*B1t Schur Complements with Pardiso SC and coping them to GPU";

	bool GPU_full = false;

	cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ ) {

		//cudaSetDevice(1);

		SparseMatrix TmpB;
		domains[i].B1_comp_dom.MatTranspose(TmpB);

		SparseSolverCPU tmpsps;
		if ( i == 0 && cluster_global_index == 1) {
			tmpsps.msglvl = Info::report(LIBRARIES) ? 1 : 0;
		}

		tmpsps.Create_SC_w_Mat( domains[i].K, TmpB, domains[i].B1Kplus, false, 0 );

		if (USE_FLOAT){
			domains[i].B1Kplus.ConvertDenseToDenseFloat( 1 );
			domains[i].B1Kplus.USE_FLOAT = true;
		}

//		SparseSolverCPU tmpsps2;
//		if ( i == 0 && cluster_global_index == 1) tmpsps2.msglvl = 1;
//		tmpsps2.Create_non_sym_SC_w_Mat( domains[i].K, TmpB, domains[i].B0t_comp, domains[i].B0KplusB1_comp, false, 0 );

//#ifdef CUDA

		//domains[i].B1Kplus.RemoveLowerDense();
		eslocal status;
		if (!GPU_full) {
			if (USE_FLOAT){
				status = domains[i].B1Kplus.CopyToCUDA_Dev_fl();
			} else {
				status = domains[i].B1Kplus.CopyToCUDA_Dev();
			}
		} else {
			status = -1;
		}

		if (status == 0) {
			domains[i].isOnACC = 1;
		} else {
			domains[i].isOnACC = 0;
			GPU_full = true;
		}
//#endif


//#ifdef CUDA
		cudaError_t status_c = cudaMallocHost((void**)&domains[i].cuda_pinned_buff, domains[i].B1_comp_dom.rows * sizeof(double));
		if (status_c != cudaSuccess) {
			ESINFO(ERROR) << "Error allocating pinned host memory";
			domains[i].isOnACC = 0;
			GPU_full = true;
		}
//#endif
		if (domains[i].isOnACC == 1) {
			ESINFO(PROGRESS2) << Info::plain() << "G";

			if (USE_FLOAT){
				SEQ_VECTOR <double> ().swap (domains[i].B1Kplus.dense_values);
			} else {
				SEQ_VECTOR <float> ().swap (domains[i].B1Kplus.dense_values_fl);
			}

		} else {
			ESINFO(PROGRESS2) << Info::plain() << "c";
		}
	}

	ESINFO(PROGRESS2);

}


void ClusterGPU::SetupKsolvers ( ) {

	cilk_for (eslocal d = 0; d < domains.size(); d++) {

		// Import of Regularized matrix K into Kplus (Sparse Solver)
		switch (config::solver::KSOLVER) {
		case 0: {
			domains[d].Kplus.ImportMatrix_wo_Copy (domains[d].K);
			break;
		}
		case 1: {
			domains[d].Kplus.ImportMatrix_wo_Copy (domains[d].K);
			break;
		}
		case 2: {
			domains[d].Kplus.ImportMatrix_fl(domains[d].K);
			break;
		}
		case 3: {
			domains[d].Kplus.ImportMatrix_fl(domains[d].K);
			break;
		}
		case 4: {
			domains[d].Kplus.ImportMatrix_fl(domains[d].K);
			break;
		}
		default:
			ESINFO(ERROR) << "Invalid KSOLVER value.";
			exit(EXIT_FAILURE);
		}

		if (config::solver::KEEP_FACTORS) {
			std::stringstream ss;
			ss << "init -> rank: " << config::MPIrank << ", subdomain: " << d;
			domains[d].Kplus.keep_factors = true;
			if (config::solver::KSOLVER != 1) {
				domains[d].Kplus.Factorization (ss.str());
			}
		} else {
			domains[d].Kplus.keep_factors = false;
			domains[d].Kplus.MPIrank = config::MPIrank;
		}

		domains[d].domain_prim_size = domains[d].Kplus.cols;

		if ( d == 0 && config::MPIrank == 0) {
			domains[d].Kplus.msglvl = Info::report(LIBRARIES) ? 1 : 0;
		}
	}

}
