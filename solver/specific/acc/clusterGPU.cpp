#include "clusterGPU.h"

using namespace espreso;

void ClusterGPU::Create_SC_perDomain(bool USE_FLOAT) {

	ESINFO(PROGRESS2) << "Creating B1*K+*B1t Schur Complements with Pardiso SC and coping them to GPU";

	bool GPU_full = false;

	//GPU_full = true;


	cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ ) {

		cudaSetDevice(1);

//		SparseSolverCPU tmpsps2;
//		if ( i == 0 && cluster_global_index == 1) tmpsps2.msglvl = 1;
//		tmpsps2.Create_non_sym_SC_w_Mat( domains[i].K, TmpB, domains[i].B0t_comp, domains[i].B0KplusB1_comp, false, 0 );

		eslocal status = 0;
		cudaError_t status_c;

		if (!GPU_full || !config::solver::COMBINE_SC_AND_SPDS) {

			GetSchurComplement(USE_FLOAT, i);

			if (!GPU_full) {
				if (USE_FLOAT){
					status = domains[i].B1Kplus.CopyToCUDA_Dev_fl();
				} else {
					status = domains[i].B1Kplus.CopyToCUDA_Dev();
				}

				if (USE_FLOAT){
					status_c = cudaMallocHost((void**)&domains[i].cuda_pinned_buff_fl, domains[i].B1_comp_dom.rows * sizeof(float));
					if (status_c != cudaSuccess) {
						ESINFO(ERROR) << "Error allocating pinned host memory";
						status = 1;
					}
				} else {
					status_c = cudaMallocHost((void**)&domains[i].cuda_pinned_buff, domains[i].B1_comp_dom.rows * sizeof(double));
					if (status_c != cudaSuccess) {
						ESINFO(ERROR) << "Error allocating pinned host memory";
						status = 1;
					}
				}
			} else {
				status = 1;
			}

			// if status == 0 - all buffers in GPU mem were sucesfuly allocated
			if (status == 0) {
				domains[i].isOnACC = 1;
				SEQ_VECTOR <double> ().swap (domains[i].B1Kplus.dense_values);
				SEQ_VECTOR <float>  ().swap (domains[i].B1Kplus.dense_values_fl);
				domains[i].Kplus.keep_factors = false;
				if (USE_FLOAT)
					ESINFO(PROGRESS2) << Info::plain() << "g";
				else
					ESINFO(PROGRESS2) << Info::plain() << "G";
			} else {
				domains[i].isOnACC = 0;
				GPU_full = true;
				if (config::solver::COMBINE_SC_AND_SPDS) {
					SEQ_VECTOR <double> ().swap (domains[i].B1Kplus.dense_values);
					SEQ_VECTOR <float>  ().swap (domains[i].B1Kplus.dense_values_fl);
					if (USE_FLOAT)
						ESINFO(PROGRESS2) << Info::plain() << "p";
					else
						ESINFO(PROGRESS2) << Info::plain() << "P";
				} else {
					if (USE_FLOAT)
						ESINFO(PROGRESS2) << Info::plain() << "c";
					else
						ESINFO(PROGRESS2) << Info::plain() << "C";
				}
			}

		} else {
                        domains[i].isOnACC = 0;
			if (USE_FLOAT)
				ESINFO(PROGRESS2) << Info::plain() << "p";
			else
				ESINFO(PROGRESS2) << Info::plain() << "P";
		}

		//GPU_full = true;

	}

	ESINFO(PROGRESS2);

}

void ClusterGPU::GetSchurComplement( bool USE_FLOAT, eslocal i ) {

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

	//ESINFO(PROGRESS2) << Info::plain() << "s";

	// if Schur complement is symmetric - then remove lower part - slower for GPU but more mem. efficient
	if (config::solver::SCHUR_COMPLEMENT_TYPE == 1)
		domains[i].B1Kplus.RemoveLowerDense();

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


			if (!config::solver::COMBINE_SC_AND_SPDS) {
				std::stringstream ss;
				ss << "init -> rank: " << config::MPIrank << ", subdomain: " << d;
				domains[d].Kplus.keep_factors = true;
				if (config::solver::KSOLVER != 1) {
					domains[d].Kplus.Factorization (ss.str());
				}
			} else {
				if ( domains[d].isOnACC == 0 ) {
					std::stringstream ss;
					ss << "init -> rank: " << config::MPIrank << ", subdomain: " << d;
					domains[d].Kplus.keep_factors = true;
					if (config::solver::KSOLVER != 1) {
						domains[d].Kplus.Factorization (ss.str());
					}
				}
			}

		} else {
			domains[d].Kplus.keep_factors = false;
			domains[d].Kplus.MPIrank = config::MPIrank;
		}

		domains[d].domain_prim_size = domains[d].Kplus.cols;

		if ( d == 0 && config::MPIrank == 0) {
			domains[d].Kplus.msglvl = 0; //Info::report(LIBRARIES) ? 1 : 0;
		}
	}

}


void ClusterGPU::multKplusGlobal_GPU(SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {

	//ESINFO(PROGRESS2) << "K+ multiply HFETI";
	mkl_set_num_threads(1);

	cluster_time.totalTime.start();

	vec_fill_time.start();
	fill(vec_g0.begin(), vec_g0.end(), 0); // reset entire vector to 0
	//fill(vec_e0.begin(), vec_e0.end(), 0); // reset entire vector to 0
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


//	for (int i = 0; i < vec_g0.size(); i++)
//	printf (       "Test probe 1: %d norm = %1.30f \n", i, vec_g0[i] );

	clus_F0_1_time.start();
	F0.Solve(vec_g0, tm1[0], 0, 0);
	clus_F0_1_time.end();

//	for (int i = 0; i < tm1[0].size(); i++)
//	printf (       "Test probe 2: %d norm = %1.30f \n", i, tm1[0][i] );

	clus_G0_time.start();
	G0.MatVec(tm1[0], tm2[0], 'N');
	clus_G0_time.end();

//	for (int i = 0; i < tm1[0].size(); i++)
//	printf (       "Test probe 3: %d norm = %1.30f \n", i, tm1[0][i] );

	cilk_for (eslocal i = 0; i < vec_e0.size(); i++)
		tm2[0][i] = tm2[0][i] - vec_e0[i];
	//cblas_daxpy(vec_e0.size(), -1.0, &vec_e0[0], 1, &tm2[0][0], 1);

	 clus_Sa_time.start();
#ifdef SPARSE_SA
	 Sa.Solve(tm2[0], vec_alfa,0,0);
#else
	char U = 'U';
	eslocal nrhs = 1;
	eslocal info = 0;
	vec_alfa = tm2[0];
	dsptrs( &U, &SaMat.rows, &nrhs, &SaMat.dense_values[0], &SaMat.ipiv[0], &vec_alfa[0], &SaMat.rows, &info );
#endif
	 clus_Sa_time.end();

//		for (int i = 0; i < vec_alfa.size(); i++)
//		printf (       "Test probe 4: %d norm = %1.30f \n", i, vec_alfa[i] );

	 clus_G0t_time.start();
	G0.MatVec(vec_alfa, tm1[0], 'T'); 	// lambda
	 clus_G0t_time.end();

//		for (int i = 0; i < tm1[0].size(); i++)
//		printf (       "Test probe 5: %d norm = %1.30f \n", i, tm1[0][i] );

	cilk_for (eslocal i = 0; i < vec_g0.size(); i++)
		tm1[0][i] = vec_g0[i] - tm1[0][i];


	clus_F0_2_time.start();
	F0.Solve(tm1[0], vec_lambda,0,0);
	clus_F0_2_time.end();

	clusCP_time.end();

//	for (int i = 0; i < vec_lambda.size(); i++)
//	printf (       "Test probe 6: %d norm = %1.30f \n", i, vec_lambda[i] );

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

		//ESINFO(PROGRESS2) << Info::plain() << ".";
	}
	//ESINFO(PROGRESS2);
	loop_2_1_time.end();

	cluster_time.totalTime.end();
}
