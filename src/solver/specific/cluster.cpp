
#include "../specific/cluster.h"

//#define SPARSE_SA

// *******************************************************************
// **** CLUSTER CLASS ************************************************

using namespace espreso;

void ClusterBase::ShowTiming()  {

	cluster_time.addEvent(vec_fill_time);
	cluster_time.addEvent(loop_1_1_time);
	cluster_time.addEvent(loop_1_2_time);

	cluster_time.addEvent(clus_F0_1_time);
	cluster_time.addEvent(clus_G0_time);
	cluster_time.addEvent(clus_Sa_time);
	cluster_time.addEvent(clus_G0t_time);
	cluster_time.addEvent(clus_F0_2_time);

	cluster_time.addEvent(clusCP_time);
	cluster_time.addEvent(loop_2_1_time);

	cluster_time.printStatsMPI();
}

void ClusterBase::SetDynamicParameters(double set_dynamic_timestep, double set_dynamic_beta, double set_dynamic_gama) {

	dynamic_timestep = set_dynamic_timestep;
	dynamic_beta     = set_dynamic_beta;
	dynamic_gama     = set_dynamic_gama;

	for (eslocal d = 0; d < domains.size(); d++)
		domains[d].SetDynamicParameters(dynamic_timestep, dynamic_beta, dynamic_gama);

}


void ClusterBase::InitClusterPC( eslocal * subdomains_global_indices, eslocal number_of_subdomains ) {

	// *** Init the vector of domains *****************************************************
	//LoadBinVectorInt(domains_in_global_index, string(path) + string(filename_DOMAINS));
	domains_in_global_index.resize( number_of_subdomains ) ;
	// domains_in_global_index[0] = index_of_first_subdomain ;
	// *** END - Init the vector of domains ***********************************************

	domains.resize( number_of_subdomains );

	if (USE_HFETI == 1) {
		for (eslocal i = 0; i < number_of_subdomains; i++)									// HFETI
			domains_in_global_index[i] = subdomains_global_indices[i] + 1;				// HFETI; [+1] -> domain numbering in espreso is from 1
	} else {
		//domains_in_global_index[0]     = 1 + ( number_of_subdomains * ( domains_in_global_index[0] - 1 ) );	// MFETI
		//for (eslocal i = 1; i < number_of_subdomains; i++)														// MFETI
		//	domains_in_global_index[i] = domains_in_global_index[0] + i;									    // MFETI
		for (eslocal i = 0; i < number_of_subdomains; i++)									// MFETI
			domains_in_global_index[i] = subdomains_global_indices[i] + 1;				// MFETI; [+1] -> domain numbering in espreso is from 1
	}

	// *** Init all domains of the cluster ********************************************
	for (eslocal i = 0; i < number_of_subdomains; i++ ) {
		domains[i].domain_global_index = domains_in_global_index[i];
		domains[i].USE_KINV    	 = USE_KINV;
		domains[i].USE_HFETI   	 = USE_HFETI;
		domains[i].USE_DYNAMIC 	 = USE_DYNAMIC;
		domains[i].DOFS_PER_NODE = DOFS_PER_NODE;
		domains[i].domain_index  = i;

	}
	// *** END - Init all domains of the cluster ***************************************
}

void ClusterBase::SetClusterPC( SEQ_VECTOR <SEQ_VECTOR <eslocal> > & lambda_map_sub ) {


	//map <eslocal,eslocal> my_lamdas_map_indices;


	//// *** Set up the dual size ********************************************************
	int MPIrank; 	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
	dual_size = domains[0].B1.rows;

	if (USE_HFETI == 1) {

		// *** Alocate temporarly vectors for inter-cluster processing *********************
		// *** - based on uncompressed matrix B0
		tm1.resize(domains.size());
		tm2.resize(domains.size());
		tm3.resize(domains.size());

		cilk_for (eslocal d = 0; d < domains.size(); d++) {
			eslocal max_tmp_vec_size = domains[d].B0.cols;

			if (domains[d].B0.rows > domains[d].B0.cols)
				max_tmp_vec_size = domains[d].B0.rows;

			tm1[d].resize( max_tmp_vec_size );
			tm2[d].resize( max_tmp_vec_size );
			tm3[d].resize( max_tmp_vec_size );
		}
		// *** END - Alocate temporarly vectors for inter-cluster processing *****************
	}


	//// *** Detection of affinity of lag. multipliers to specific subdomains **************
	//// *** - will be used to compress vectors and matrices for higher efficiency
	//my_neighs = neigh_domains;


	ESLOG(MEMORY) << "Setting vectors for lambdas";
	ESLOG(MEMORY) << "process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";
	ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";


	my_lamdas_indices.resize( lambda_map_sub.size() );
	for (eslocal i = 0; i < lambda_map_sub.size(); i++)
		my_lamdas_indices[i] = lambda_map_sub[i][0];


	SEQ_VECTOR< SEQ_VECTOR <eslocal> > lambdas_per_subdomain ( domains.size() * NUMBER_OF_CLUSTERS );
	my_lamdas_ddot_filter.resize( lambda_map_sub.size(), 0.0 );
	for (eslocal i = 0; i < lambda_map_sub.size(); i++) {
		if ( lambda_map_sub[i].size() > 2 ) {
			if ( lambda_map_sub[i][1] < lambda_map_sub[i][2] )
				my_lamdas_ddot_filter[i] = 1.0;

			 lambdas_per_subdomain[lambda_map_sub[i][1]].push_back(lambda_map_sub[i][0]);
			 lambdas_per_subdomain[lambda_map_sub[i][2]].push_back(lambda_map_sub[i][0]);

		} else {
			my_lamdas_ddot_filter[i] = 1.0;
		}
	}

	ESLOG(MEMORY) << "Setting vectors for lambdas communicators";
	ESLOG(MEMORY) << "process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";
	ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

	my_comm_lambdas_indices .resize(my_neighs.size());
	my_comm_lambdas			.resize(my_neighs.size());
	my_recv_lambdas			.resize(my_neighs.size());

	//ESLOG(MEMORY) << "1 process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";

	cilk_for (eslocal i = 0; i < my_neighs.size(); i++) {
		my_comm_lambdas_indices[i] = lambdas_per_subdomain[my_neighs[i]];
		my_comm_lambdas[i].resize(my_comm_lambdas_indices[i].size());
		my_recv_lambdas[i].resize(my_comm_lambdas_indices[i].size());
	}

	//ESLOG(MEMORY) << "2 process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";

	compressed_tmp    .resize( my_lamdas_indices.size(), 0 );
	//compressed_tmp2   .resize( my_lamdas_indices.size(), 0 );

	//ESLOG(MEMORY) << "3 process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";

	cilk_for (eslocal d = 0; d < domains.size(); d++ )
		if (USE_KINV == 1 ) {
			domains[d].compressed_tmp.resize( domains[d].B1.I_row_indices.size(), 0);
			domains[d].compressed_tmp2.resize( domains[d].B1.I_row_indices.size(), 0);
		} else {
			domains[d].compressed_tmp.resize( 1, 0);
			domains[d].compressed_tmp2.resize( 1, 0);
		}

	// mapping/compression vector for cluster
	for (eslocal i = 0; i <my_lamdas_indices.size(); i++)
		_my_lamdas_map_indices.insert(make_pair(my_lamdas_indices[i],i));

	//ESLOG(MEMORY) << "5 process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";

	// mapping/compression vector for domains
	cilk_for (eslocal i = 0; i < domains.size(); i++) {
		for (eslocal j = 0; j < domains[i].lambda_map_sub.size(); j++) {
			domains[i].my_lamdas_map_indices.insert(make_pair(domains[i].lambda_map_sub[j] ,j));
		}
	}

	//ESLOG(MEMORY) << "6 process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";

	cilk_for (eslocal d = 0; d < domains.size(); d++) {

            if (domains[d].lambda_map_sub.size() > 0 ) {

                eslocal i = 0;
				eslocal j = 0;
				do
				{
					eslocal big_index   = my_lamdas_indices[i];
					eslocal small_index = domains[d].lambda_map_sub[j];

					if (big_index >  small_index) j++;

					if (big_index  < small_index) i++;

					if (big_index == small_index) {
						domains[d].lambda_map_sub_local.push_back(i);
						i++; j++;
					}


				} while ( i < my_lamdas_indices.size() && j < domains[d].lambda_map_sub.size() );
            }
        }
	//// *** END - Detection of affinity of lag. multipliers to specific subdomains ***************

	//ESLOG(MEMORY) << "7 process " << config::evn::MPIrank << " uses " << Measure::processMemory() << " MB";


	//// *** Create a vector of communication pattern needed for AllReduceLambdas function *******
	my_comm_lambdas_indices_comp.resize(my_neighs.size());
	cilk_for (eslocal i = 0; i < my_neighs.size(); i++) {
		my_comm_lambdas_indices_comp[i].resize( lambdas_per_subdomain[my_neighs[i]].size() );
		for (eslocal j = 0; j < lambdas_per_subdomain[my_neighs[i]].size(); j++ )
			my_comm_lambdas_indices_comp[i][j] = _my_lamdas_map_indices[lambdas_per_subdomain[my_neighs[i]][j]];
	}
	//// *** END - Create a vector of communication pattern needed for AllReduceLambdas function *

	//ESLOG(MEMORY) << "8 process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";


	//// *** Compression of Matrix B1 to work with compressed lambda vectors *****************
	ESLOG(MEMORY) << "B1 compression";
	ESLOG(MEMORY) << "process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";
	ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

	cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ ) {

		domains[i].B1_comp_dom.I_row_indices = domains[i].B1.I_row_indices;
		domains[i].B1_comp_dom.J_col_indices = domains[i].B1.J_col_indices;
		domains[i].B1_comp_dom.V_values      = domains[i].B1.V_values;

		domains[i].B1_comp_dom.rows = domains[i].B1.rows;
		domains[i].B1_comp_dom.cols = domains[i].B1.cols;
		domains[i].B1_comp_dom.nnz  = domains[i].B1.nnz;
		domains[i].B1_comp_dom.type = domains[i].B1.type;

		for (eslocal j = 0; j < domains[i].B1_comp_dom.I_row_indices.size(); j++ ) {
			eslocal tmp_new = domains[i].my_lamdas_map_indices[domains[i].B1_comp_dom.I_row_indices [j] - 1] + 1;  // numbering from 1 in matrix
			domains[i].B1_comp_dom.I_row_indices [j] = tmp_new;									               // j + 1; // numbering matrix from 1
		}

		domains[i].B1_comp_dom.rows = domains[i].lambda_map_sub.size();
		domains[i].B1_comp_dom.ConvertToCSRwithSort( 1 );

		if (config::solver::PRECONDITIONER == config::solver::PRECONDITIONERalternative::DIRICHLET || 
        config::solver::PRECONDITIONER == config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET) {
			domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_DirPr);
			auto last = std::unique(domains[i].B1t_DirPr.CSR_I_row_indices.begin(), domains[i].B1t_DirPr.CSR_I_row_indices.end());
			domains[i].B1t_DirPr.CSR_I_row_indices.erase(last, domains[i].B1t_DirPr.CSR_I_row_indices.end());
			domains[i].B1t_DirPr.rows = domains[i].B1t_DirPr.CSR_I_row_indices.size();

			domains[i].B1t_Dir_perm_vec = domains[i].B1_comp_dom.CSR_J_col_indices;
			std::sort(domains[i].B1t_Dir_perm_vec.begin(), domains[i].B1t_Dir_perm_vec.end());
			last = std::unique(domains[i].B1t_Dir_perm_vec.begin(), domains[i].B1t_Dir_perm_vec.end());
			domains[i].B1t_Dir_perm_vec.erase(last, domains[i].B1t_Dir_perm_vec.end() );
		}

		//************************
		if ( config::solver::REGULARIZATION == 0 )
                    domains[i].B1.Clear();

		domains[i].B1t.Clear();
//		domains[i].B1_comp.Clear();
//		domains[i].B1t_comp.Clear();

		domains[i].my_lamdas_map_indices.clear();

	}


	//// *** END - Compression of Matrix B1 to work with compressed lambda vectors *************



	//// *** Compression of Matrix G1 to work with compressed lambda vectors *******************
	ESLOG(MEMORY) << "G1 compression";
	ESLOG(MEMORY) << "process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";
	ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

	if (USE_DYNAMIC == 0) {
		if ( ! (USE_HFETI == 1 && config::solver::REGULARIZATION == 1 )) {
			Compress_G1();
		}
	}
	//// *** END - Compression of Matrix G1 to work with compressed lambda vectors ***************

	ESLOG(MEMORY) << "Lambdas end";
	ESLOG(MEMORY) << "process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";
	ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

}

void ClusterBase::ImportKmatrixAndRegularize ( SEQ_VECTOR <SparseMatrix> & K_in, SEQ_VECTOR <SparseMatrix> & RegMat ) {

	cilk_for (eslocal d = 0; d < domains.size(); d++) {
		if ( d == 0 && config::env::MPIrank == 0) {
			domains[d].Kplus.msglvl = Info::report(LIBRARIES) ? 1 : 0;
		}


		domains[d].K.swap(K_in[d]);
		domains[d]._RegMat.swap(RegMat[d]);


		if ( config::solver::PRECONDITIONER == config::solver::PRECONDITIONERalternative::MAGIC ) {
			//domains[d].Prec = domains[d].K;
			domains[d]._RegMat.ConvertToCSR(1);
			domains[d].Prec.MatAdd(domains[d].K, domains[d]._RegMat, 'N', -1);
			domains[d]._RegMat.ConvertToCOO(1);
		}

		domains[d].enable_SP_refinement = true;
	    //}
	    ESINFO(PROGRESS2) << Info::plain() << ".";
	}
	ESINFO(PROGRESS2);
}


void ClusterBase::SetClusterHFETI () {
	// *** Create Matrices and allocate vectors for Hybrid FETI **************************
	if (USE_HFETI == 1) {

		TimeEval HFETI_prec_timing (" HFETI - preprocessing timing");
		HFETI_prec_timing.totalTime.start();

		int MPIrank;
		MPI_Comm_rank (MPI_COMM_WORLD, &MPIrank);
		ESINFO(PROGRESS2) << "HFETI preprocessing start";

		TimeEvent B0_time("Compress B0 per cluster");
		B0_time.start();

		CompressB0();

		B0_time.end();
		B0_time.printStatMPI();
		HFETI_prec_timing.addEvent(B0_time);


		TimeEvent G0_time("Create G0 per cluster");
		G0_time.start();

		CreateG0();
		vec_g0.resize(G0.cols);
		vec_e0.resize(G0.rows);

		G0_time.end();
		G0_time.printStatMPI();
		HFETI_prec_timing.addEvent(G0_time);


		TimeEvent F0_time("Create F0 per cluster");
		F0_time.start();

		CreateF0();
		vec_lambda.resize(F0.m_Kplus_size);

		F0_time.end();
		HFETI_prec_timing.addEvent(F0_time);


		TimeEvent Sa_time("Create Salfa per cluster");
		Sa_time.start();

		CreateSa();
		vec_alfa.resize(Sa.m_Kplus_size);

		Sa_time.end();
		HFETI_prec_timing.addEvent(Sa_time);

		HFETI_prec_timing.totalTime.end();
		HFETI_prec_timing.printStatsMPI();

	}
	// *** END - Create Matrices for Hybrid FETI *****************************************
}

void ClusterBase::SetClusterPC_AfterKplus () {

	//// *** Alocate temporarly vectors for Temporary vectors for Apply_A function *********
	//// *** - temporary vectors for work primal domain size *******************************
	x_prim_cluster1.resize( domains.size() );
	x_prim_cluster2.resize( domains.size() );
	x_prim_cluster3.resize( domains.size() );

	for (eslocal d = 0; d < domains.size(); d++) {
		x_prim_cluster1[d].resize( domains[d].domain_prim_size );
		x_prim_cluster2[d].resize( domains[d].domain_prim_size );
		x_prim_cluster3[d].resize( domains[d].domain_prim_size );
	}
	//// *** END - Alocate temporarly vectors for Temporary vectors for Apply_A function ***

	//// *** Prepare the initial right hand side in dual *************************************
	//CreateVec_b_perCluster();

}


void ClusterBase::multKplusGlobal(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, SEQ_VECTOR<eslocal> & cluster_map_vec) {

//	vec_g0.resize(G0.cols);
//	fill(vec_g0.begin(), vec_g0.end(), 0); // reset entire vector to 0
//
//	vec_e0.resize(G0.rows);
//	fill(vec_e0.begin(), vec_e0.end(), 0); // reset entire vector to 0
//
//	//y_out.resize(x_clust_size);
//
//	// temp vectors
//	SEQ_VECTOR <double> t1,t2,t3;
//
//
//	// loop over domains in the cluster
//	for (eslocal d = 0; d < domains.size(); d++) {
//		eslocal x_in_vector_start_index = x_clust_domain_map_vec[d] + cluster_map_vec[this->cluster_global_index-1];
//		eslocal domain_size = domains[d].Kplus.m_Kplus_size;
//
//		t1.resize(domain_size);
//		t2.resize(vec_g0.size());
//		fill(t1.begin(), t1.end(), 0);
//		fill(t2.begin(), t2.end(), 0);
//
//
//		// g0
//		domains[d].multKplusLocal( x_in, t1, x_in_vector_start_index, 0 );
//
//		//					  in   out trans
//		domains[d].B0.MatVec( t1 , t2,  'N' );
//
//		for (eslocal i = 0; i < vec_g0.size(); i++ )
//			vec_g0[i] = vec_g0[i] + t2[i];
//
//		// e0
//		t2.resize(domains[d].Kplus_R.cols); // resize na pocet sloupcu matice Kplus_R - velikost jadra
//		fill(t2.begin(), t2.end(), 0); // reset t2 - migh not be necessary
//
//		domains[d].Kplus_R.MatVec(x_in, t2, 'T', x_in_vector_start_index, 0);
//
//		eslocal e0_start	=  d	* domains[d].Kplus_R.cols;
//		eslocal e0_end		= (d+1) * domains[d].Kplus_R.cols;
//		for (eslocal i = e0_start; i < e0_end; i++ )
//			vec_e0[i] = - t2[i - e0_start];
//
//	} // end loop over domains
//
//
//	// alfa
//	t1.resize(F0.m_Kplus_size);
//	fill(t1.begin(), t1.end(), 0);
//
//	F0.Solve(vec_g0, t1,0,0);
//
//	t2.resize(G0.rows);
//	fill(t2.begin(), t2.end(), 0);
//
//	G0.MatVec(t1, t2, 'N');
//	for (eslocal i = 0; i < vec_e0.size(); i++)
//		t2[i] = t2[i] - vec_e0[i];
//
//	vec_alfa.resize(Sa.m_Kplus_size);
//	Sa.Solve(t2, vec_alfa,0,0);
//
//	// lambda
//	t1.resize(G0.cols);
//	fill(t1.begin(), t1.end(), 0);
//
//	G0.MatVec(vec_alfa, t1, 'T');
//
//	for (eslocal i = 0; i < vec_g0.size(); i++)
//		t1[i] = vec_g0[i] - t1[i];
//
//	vec_lambda.resize(F0.m_Kplus_size);
//	F0.Solve(t1, vec_lambda,0,0);
//
//	// Kplus_x
//	for (eslocal d = 0; d < domains.size(); d++) {
//		//eslocal x_in_vector_start_index = x_clust_domain_map_vec[d];
//		eslocal x_in_vector_start_index = x_clust_domain_map_vec[d] + cluster_map_vec[this->cluster_global_index-1];
//		eslocal domain_size = domains[d].Kplus.m_Kplus_size;
//
//		t1.resize(domain_size);
//		fill(t1.begin(), t1.end(), 0);
//
//		domains[d].B0.MatVec(vec_lambda, t1, 'T');
//
//		for (eslocal i = 0; i < domain_size; i++)
//			t1[i] = x_in[x_in_vector_start_index + i] - t1[i];
//
//		t2.resize(domain_size);
//		fill(t2.begin(), t2.end(), 0);
//
//		domains[d].multKplusLocal(t1 , t2, 0, 0);
//
//		t3.resize(domain_size);
//		fill(t3.begin(), t3.end(), 0);
//
//		eslocal e0_start	=  d	* domains[d].Kplus_R.cols;
//		eslocal e0_end		= (d+1) * domains[d].Kplus_R.cols;
//		domains[d].Kplus_R.MatVec(vec_alfa, t3, 'N', e0_start,0);
//
//		for (eslocal i = 0; i < domain_size; i++)
//			y_out[x_in_vector_start_index + i] = t2[i] + t3[i];
//
//		double t_sum = 0;
//		for (eslocal i = 0; i < domain_size; i++)
//			t_sum +=  (t2[i] + t3[i]);
//
//		t_sum = t_sum;
//
//	}
//
//	t1.clear();
//	t2.clear();
//	t3.clear();
}

void ClusterBase::multKplusGlobal_l(SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {

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
//#ifdef SPARSE_SA
//	 Sa.Solve(tm2[0], vec_alfa,0,0);
//#else
//	eslocal nrhs = 1;
//	Sa_dense.Solve(tm2[0], vec_alfa, nrhs);
//#endif

	 if (config::solver::SA_SOLVER == config::SA_SPARSE_on_CPU) {
		 Sa.Solve(tm2[0], vec_alfa,0,0);
	 }

	 if (config::solver::SA_SOLVER == config::SA_DENSE_on_CPU) {
		eslocal nrhs = 1;
		Sa_dense_cpu.Solve(tm2[0], vec_alfa, nrhs);
	 }

	 if (config::solver::SA_SOLVER == config::SA_DENSE_on_ACC) {
		eslocal nrhs = 1;
		Sa_dense_acc.Solve(tm2[0], vec_alfa, nrhs);
	 }

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

		SEQ_VECTOR < double > tmp_vec (domains[d].B0_comp_map_vec.size(), 0);
		for (eslocal i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
			tmp_vec[i] = vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;
		//domains[d].B0t_comp.MatVec(tmp_vec, tm1[d], 'N');
		domains[d].B0_comp.MatVec(tmp_vec, tm1[d], 'T');

		for (eslocal i = 0; i < domain_size; i++)
			tm1[d][i] = x_in[d][i] - tm1[d][i];

		domains[d].multKplusLocal(tm1[d] , tm2[d]);

		eslocal e0_start	=  d	* domains[d].Kplus_R.cols;
		eslocal e0_end		= (d+1) * domains[d].Kplus_R.cols;

		domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);

		for (eslocal i = 0; i < domain_size; i++)
			x_in[d][i] = tm2[d][i] + tm3[d][i];

		//ESINFO(PROGRESS2) << Info::plain() << ".";
	}
	//ESINFO(PROGRESS2);
	loop_2_1_time.end();

	cluster_time.totalTime.end();
}

void ClusterBase::multKplusGlobal_Kinv( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in ) {

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
		domains[d].Kplus_R.     DenseMatVec(x_in[d], tm3[d], 'T');		// e0
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
// #ifdef SPARSE_SA
//    Sa.Solve(tm2[0], vec_alfa,0,0);
// #else
//    eslocal nrhs = 1;
//    Sa_dense.Solve(tm2[0], vec_alfa, nrhs);
// #endif

	 if (config::solver::SA_SOLVER == config::SA_SPARSE_on_CPU) {
		 Sa.Solve(tm2[0], vec_alfa,0,0);
	 }

	 if (config::solver::SA_SOLVER == config::SA_DENSE_on_CPU) {
		eslocal nrhs = 1;
		Sa_dense_cpu.Solve(tm2[0], vec_alfa, nrhs);
	 }

	 if (config::solver::SA_SOLVER == config::SA_DENSE_on_ACC) {
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
		for (eslocal i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
			tmp_vec[i] = -1.0 * vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;
		domains[d].B0Kplus_comp.DenseMatVec(tmp_vec, tm2[d], 'T' );

		eslocal e0_start	=  d	* domains[d].Kplus_R.cols;
		eslocal e0_end		= (d+1) * domains[d].Kplus_R.cols;
		domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);

		for (eslocal i = 0; i < domain_size; i++)
			x_in[d][i] = tm2[d][i] + tm3[d][i];

	}
	loop_2_1_time.end();

	cluster_time.totalTime.end();
}

void ClusterBase::multKplusGlobal_Kinv_2( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in ) {

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
   		//domains[d].B0Kplus_comp.DenseMatVec(x_in[d], tm2[d]);				// g0 - with comp B0Kplus
		//domains[d].Kplus_R.     DenseMatVec(x_in[d], tm3[d], 'T');		// e0
		domains[d].B0KplusB1_comp .DenseMatVec(x_in[d], tm2[d], 'N');		// g0 - with comp B0Kplus
		domains[d].Kplus_R_B1_comp.DenseMatVec(x_in[d], tm3[d], 'N');		// e0
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
//	Sa.Solve(tm2[0], vec_alfa,0,0);
	 if (config::solver::SA_SOLVER == config::SA_SPARSE_on_CPU) {
		 Sa.Solve(tm2[0], vec_alfa,0,0);
	 }

	 if (config::solver::SA_SOLVER == config::SA_DENSE_on_CPU) {
		eslocal nrhs = 1;
		Sa_dense_cpu.Solve(tm2[0], vec_alfa, nrhs);
	 }

	 if (config::solver::SA_SOLVER == config::SA_DENSE_on_ACC) {
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
		for (eslocal i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
			tmp_vec[i] = -1.0 * vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;

		//domains[d].B0Kplus_comp.DenseMatVec(tmp_vec, tm2[d], 'T' );
		domains[d].B0KplusB1_comp .DenseMatVec(tmp_vec, tm2[d], 'T');
		domains[d].B0KplusB1_comp .MatVec(tmp_vec, tm2[d], 'T');

		eslocal e0_start	=  d	* domains[d].Kplus_R.cols;
		eslocal e0_end		= (d+1) * domains[d].Kplus_R.cols;

		//domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);
		domains[d].Kplus_R_B1_comp.DenseMatVec(vec_alfa, tm3[d], 'T', e0_start);
		//domains[d].Kplus_R_B1_comp.MatVec(vec_alfa, tm3[d], 'T', e0_start);

		for (eslocal i = 0; i < domain_size; i++)
			x_in[d][i] = tm2[d][i] + tm3[d][i];

	}
	loop_2_1_time.end();

	cluster_time.totalTime.end();
}


////backup March 31 2015
//void ClusterBase::multKplusGlobal_l(SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {
//
//	//eslocal MPIrank;
//	//MPI_Comm_rank (MPI_COMM_WORLD, &MPIrank);
//	//if (MPIrank == 0 ) { cout << "MultKplusGlobal - Cilk workers = " << __cilkrts_get_nworkers()      << endl; }
//	//if (MPIrank == 0 ) { cout << "MultKplusGlobal - Cilk workers = " << __cilkrts_get_total_workers() << endl; }
//	//
//
//	eslocal num_threads = domains.size();
//
//	//__cilkrts_set_param("nworkers", num_threads);
//	mkl_set_num_threads(1);
//
//	cluster_time.totalTime.start();
//
//	vec_fill_time.start();
//	fill(vec_g0.begin(), vec_g0.end(), 0); // reset entire vector to 0
//	//fill(vec_e0.begin(), vec_e0.end(), 0); // reset entire vector to 0
//	vec_fill_time.end();
//
//	// loop over domains in the cluster
//	loop_1_1_time.start();
//	cilk_for (eslocal d = 0; d < domains.size(); d++)
//	{
//		//PROD - domains[d].B0Kplus.MatVec(x_in[d], tm2[d], 'N');			// g0 using B0Kplus
//		//cout << "B0Klus - sparsity - " << 100.0 * (double)domains[d].B0Kplus.nnz / (double)(domains[d].B0Kplus.cols * domains[d].B0Kplus.rows) << endl;
//
//		//OLD - domains[d].multKplusLocal( x_in[d], tm1[d], 0, 0 );		// g0
//		//OLD - domains[d].B0.MatVec( tm1[d] , tm2[d],  'N' );			// g0
//
//		domains[d].B0Kplus_comp.DenseMatVec(x_in[d], tm2[d]);			// g0 - with comp B0
//
//		//domains[d].Kplus_R.MatVec     (x_in[d], tm3[d], 'T');		        // e0
//		domains[d].Kplus_R.DenseMatVec(x_in[d], tm3[d], 'T');
//	}
//	loop_1_1_time.end();
//
//	loop_1_2_time.start();
//	cilk_for (eslocal d = 0; d < domains.size(); d++)
//	{
//		eslocal e0_start	=  d	* domains[d].Kplus_R.cols;
//		eslocal e0_end		= (d+1) * domains[d].Kplus_R.cols;
//
//		for (eslocal i = e0_start; i < e0_end; i++ )
//			vec_e0[i] = - tm3[d][i - e0_start];
//	}
//
//
//	for (eslocal d = 0; d < domains.size(); d++)
//		for (eslocal i = 0; i < domains[d].B0Kplus_comp.rows; i++)
//			vec_g0[ domains[d].B0_comp_map_vec[i] - 1 ] += tm2[d][i];
//
//
//	//PROD -
//	//cilk_for (eslocal i = 0; i < vec_g0.size(); i++ )
//	//{
//	//	vec_g0[i] = tm2[0][i];
//	//}
//	//for (eslocal d = 1; d < domains.size(); d++) {
//	//	cilk_for (eslocal i = 0; i < vec_g0.size(); i++ )
//	//	{
//	//		vec_g0[i] = vec_g0[i] + tm2[d][i];
//	//	}
//	//}
//
//
//	// end loop over domains
//	loop_1_2_time.end();
//
//	mkl_set_num_threads(24); //pozor
//	clusCP_time.start();
//
//	clus_F0_1_time.start();
//	F0.Solve(vec_g0, tm1[0], 0, 0);
//	clus_F0_1_time.end();
//
//	clus_G0_time.start();
//	G0.MatVec(tm1[0], tm2[0], 'N');
//	clus_G0_time.end();
//
//	cilk_for (eslocal i = 0; i < vec_e0.size(); i++)
//		tm2[0][i] = tm2[0][i] - vec_e0[i];
//	//cblas_daxpy(vec_e0.size(), -1.0, &vec_e0[0], 1, &tm2[0][0], 1);
//
//	clus_Sa_time.start();
//	Sa.Solve(tm2[0], vec_alfa,0,0);
//	clus_Sa_time.end();
//
//	clus_G0t_time.start();
//	G0.MatVec(vec_alfa, tm1[0], 'T'); 	// lambda
//	clus_G0t_time.end();
//
//	cilk_for (eslocal i = 0; i < vec_g0.size(); i++)
//		tm1[0][i] = vec_g0[i] - tm1[0][i];
//
//
//	clus_F0_2_time.start();
//	F0.Solve(tm1[0], vec_lambda,0,0);
//	clus_F0_2_time.end();
//
//	clusCP_time.end();
//
//
//	// Kplus_x
//	mkl_set_num_threads(1);
//	loop_2_1_time.start();
//	cilk_for (eslocal d = 0; d < domains.size(); d++)
//	{
//		eslocal domain_size = domains[d].domain_prim_size;
//
//
//		SEQ_VECTOR < double > tmp_vec (domains[d].B0_comp_map_vec.size(), 0);
//		for (eslocal i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
//			tmp_vec[i] = vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;
//		domains[d].B0t_comp.MatVec(tmp_vec, tm1[d], 'N');
//
//
//		//domains[d].B0.MatVec(vec_lambda, tm1[d], 'T');  // both about the same performance
//		//domains[d].B0t.MatVec(vec_lambda, tm1[d], 'N');   // both about the same performance
//
//
//		for (eslocal i = 0; i < domain_size; i++)
//			tm1[d][i] = x_in[d][i] - tm1[d][i];
//
//		domains[d].multKplusLocal(tm1[d] , tm2[d], 0, 0);
//
//		eslocal e0_start	=  d	* domains[d].Kplus_R.cols;
//		eslocal e0_end		= (d+1) * domains[d].Kplus_R.cols;
//
//		//domains[d].Kplus_R.MatVec(vec_alfa, tm3[d], 'N', e0_start,0);
//		domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);
//
//		for (eslocal i = 0; i < domain_size; i++)
//			x_in[d][i] = tm2[d][i] + tm3[d][i];
//
//	}
//	loop_2_1_time.end();
//
//	cluster_time.totalTime.end();
//}


void ClusterBase::CompressB0() {

	cilk_for (eslocal d = 0; d < domains.size(); d++) {
		domains[d].B0.MatTranspose(domains[d].B0t);
		domains[d].B0_comp = domains[d].B0;

		for (size_t i = 0; i < domains[d].B0_comp.rows; i++) {
			if (domains[d].B0_comp.CSR_I_row_indices[i] != domains[d].B0_comp.CSR_I_row_indices[i + 1]) {
				domains[d].B0_comp_map_vec.push_back(i + 1);
			}
		}

		auto it = std::unique(domains[d].B0_comp.CSR_I_row_indices.begin(), domains[d].B0_comp.CSR_I_row_indices.end());
		domains[d].B0_comp.rows = std::distance(domains[d].B0_comp.CSR_I_row_indices.begin(), it) - 1;
		domains[d].B0_comp.CSR_I_row_indices.resize(domains[d].B0_comp.rows + 1);

		domains[d].B0_comp.MatTranspose(domains[d].B0t_comp);
	}

}

void ClusterBase::CreateG0() {

	mkl_set_num_threads(1);

	SEQ_VECTOR <SparseMatrix> G0LocalTemp( domains.size() );

	cilk_for (eslocal i = 0; i<domains.size(); i++) {
		domains[i].Kplus_R.ConvertDenseToCSR(0);

		G0LocalTemp[i].MatMat(domains[i].B0, 'N', domains[i].Kplus_R );
		G0LocalTemp[i].MatTranspose(-1.0);

		SEQ_VECTOR<eslocal>().swap( domains[i].Kplus_R.CSR_I_row_indices );
		SEQ_VECTOR<eslocal>().swap( domains[i].Kplus_R.CSR_J_col_indices );
		SEQ_VECTOR<double> ().swap( domains[i].Kplus_R.CSR_V_values );
	}

	for (eslocal i = 0; i<domains.size(); i++) {
		G0.MatAppend(G0LocalTemp[i]);
		G0LocalTemp[i].Clear();
	}

}

void ClusterBase::CreateF0() {

	 TimeEval F0_timing (" HFETI - F0 preprocessing timing");
	 F0_timing.totalTime.start();

	mkl_set_num_threads(1);

	int MPIrank; MPI_Comm_rank (MPI_COMM_WORLD, &MPIrank);

	SEQ_VECTOR <SparseMatrix> tmpF0v (domains.size());

	ESINFO(PROGRESS2) << "HFETI - Create F0";

	 TimeEvent solve_F0_time("B0 compression; F0 multiple RHS solve");
	 solve_F0_time.start();

	cilk_for (eslocal d = 0; d < domains.size(); d++) {

		if (MPIrank == 0 && d == 0)
			domains[d].Kplus.msglvl=0;
		else
			domains[d].Kplus.msglvl=0;

		// FO solve in doble in case K is in single
		if (   (config::solver::KSOLVER == 2 || config::solver::KSOLVER == 3 )
			 && config::solver::F0_SOLVER == 1
			) {
			SparseSolverCPU Ktmp;
			Ktmp.ImportMatrix_wo_Copy(domains[d].K);
			std::stringstream ss;
			ss << "Create F0 -> rank: " << config::env::MPIrank << ", subdomain: " << d;
			Ktmp.Factorization(ss.str());
			Ktmp.SolveMat_Dense(domains[d].B0t_comp, domains[d].B0Kplus_comp);
		} else {
			// F0 uses same precision is K
			domains[d].Kplus.SolveMat_Dense(domains[d].B0t_comp, domains[d].B0Kplus_comp);
		}

		domains[d].B0t_comp.Clear();

		domains[d].B0Kplus = domains[d].B0Kplus_comp;
		domains[d].B0Kplus_comp.MatTranspose();
		domains[d].B0Kplus_comp.ConvertCSRToDense(1);

		for (eslocal i = 0; i < domains[d].B0Kplus.CSR_J_col_indices.size() - 1; i++)
			domains[d].B0Kplus.CSR_J_col_indices[i] = domains[d].B0_comp_map_vec [ domains[d].B0Kplus.CSR_J_col_indices[i] - 1 ];

		domains[d].B0Kplus.cols = domains[d].B0.rows;;

		// Reduces the work for HFETI iteration - reduces the
		// New multKlpusGlobal_Kinv2
		if ( 0 == 1 ) {
			domains[d].B0KplusB1_comp. MatMat(domains[d].B1_comp_dom, 'N', domains[d].B0Kplus);
			domains[d].B0KplusB1_comp.ConvertCSRToDense(0);

			domains[d].Kplus_R_B1_comp.MatMat(domains[d].B1_comp_dom, 'N', domains[d].Kplus_R);
			domains[d].Kplus_R_B1_comp.ConvertCSRToDense(0);
		}
		// END - New multKlpusGlobal

		tmpF0v[d].MatMat(domains[d].B0, 'N', domains[d].B0Kplus);
		domains[d].B0Kplus.Clear();

		domains[d].Kplus.msglvl=0;
		ESINFO(PROGRESS2) << Info::plain() << ".";
	}

	ESINFO(PROGRESS2);

	 solve_F0_time.end();
	 solve_F0_time.printStatMPI();
	 F0_timing.addEvent(solve_F0_time);

	 TimeEvent reduction_F0_time("F0 reduction time");
	 reduction_F0_time.start();

	for (eslocal j = 1; j < tmpF0v.size(); j = j * 2 ) {
		cilk_for (eslocal i = 0; i < tmpF0v.size(); i = i + 2*j) {
			if ( i+j < tmpF0v.size()) {
				tmpF0v[i    ].MatAddInPlace( tmpF0v[i + j], 'N', 1.0 );
				tmpF0v[i + j].Clear();
			}
		}
	}
	F0_Mat = tmpF0v[0];

	 reduction_F0_time.end(); reduction_F0_time.printStatMPI(); F0_timing.addEvent(reduction_F0_time);


	 TimeEvent fact_F0_time("B0 Kplus Factorization "); fact_F0_time.start();

	mkl_set_num_threads(PAR_NUM_THREADS);
	F0_Mat.RemoveLower();
	F0.ImportMatrix(F0_Mat);

	bool PARDISO_SC = true;
	if (!PARDISO_SC)
		F0_fast.ImportMatrix(F0_Mat);

	//F0_Mat.Clear();
	F0.SetThreaded();
	std::stringstream ss;
	ss << "F0 -> rank: " << config::env::MPIrank;
	F0.Factorization(ss.str());

	mkl_set_num_threads(1);

	if (MPIrank == 0) F0.msglvl = 0;

	 fact_F0_time.end(); fact_F0_time.printStatMPI(); F0_timing.addEvent(fact_F0_time);

	F0_timing.totalTime.end();
	F0_timing.printStatsMPI();

	// *** POZOR **************************************************************
	cilk_for (eslocal d = 0; d<domains.size(); d++) {
		domains[d].B0.Clear();
		domains[d].B0t.Clear();
	}
};

void ClusterBase::CreateSa() {

	bool PARDISO_SC = true;
	bool get_kernel_from_mesh;

	if ( config::solver::REGULARIZATION == 0 )
  		get_kernel_from_mesh = true	;
  	else
  		get_kernel_from_mesh = false	;


	MKL_Set_Num_Threads(PAR_NUM_THREADS);
	int MPIrank; MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

	SparseMatrix Salfa; SparseMatrix tmpM;

	TimeEval Sa_timing (" HFETI - Salfa preprocessing timing"); Sa_timing.totalTime.start();

	 TimeEvent G0trans_Sa_time("G0 transpose"); G0trans_Sa_time.start();
	SparseMatrix G0t;
	G0.MatTranspose(G0t);
	 G0trans_Sa_time.end(); G0trans_Sa_time.printStatMPI(); Sa_timing.addEvent(G0trans_Sa_time);

	 TimeEvent G0solve_Sa_time("SolveMatF with G0t as RHS"); G0solve_Sa_time.start();
	if (!PARDISO_SC) {
		if (MPIrank == 0) {
			F0_fast.msglvl = Info::report(LIBRARIES) ? 1 : 0;
		}
		//SolaveMatF is obsolete - use Schur Complement Instead
		F0_fast.SolveMatF(G0t,tmpM, true);
		if (MPIrank == 0) F0_fast.msglvl = 0;
	} else {
		SparseSolverCPU tmpsps;
		if (MPIrank == 0) {
			tmpsps.msglvl = Info::report(LIBRARIES) ? 1 : 0;
		}
		tmpsps.Create_SC_w_Mat( F0_Mat, G0t, Salfa, true, 0 );

        Salfa.ConvertDenseToCSR(1);
        Salfa.RemoveLower();
		if (MPIrank == 0) tmpsps.msglvl = 0;
	}
	F0_Mat.Clear();
	 G0solve_Sa_time.end(); G0solve_Sa_time.printStatMPI(); Sa_timing.addEvent(G0solve_Sa_time);

	 TimeEvent SaMatMat_Sa_time("Salfa = MatMat G0 * solution "); SaMatMat_Sa_time.start();
	if (!PARDISO_SC) {
		Salfa.MatMat(G0, 'N', tmpM);
		Salfa.RemoveLower();
		tmpM.Clear();
	}
	 SaMatMat_Sa_time.end(); SaMatMat_Sa_time.printStatMPI(); Sa_timing.addEvent(SaMatMat_Sa_time);

	 if (!get_kernel_from_mesh) {
		 SparseMatrix Kernel_Sa;
		 ESINFO(PROGRESS2) << "Salfa - regularization from matrix";

		 SparseMatrix GGt;

		 GGt.MatMat(G0,'N',G0t);
		 GGt.RemoveLower();
//		 GGt.get_kernel_from_K(GGt, Kernel_Sa);


		 double tmp_double;
		 eslocal tmp_int;
		 SparseMatrix TSak, _tmpSparseMat;
		 Salfa.get_kernel_from_K(Salfa,_tmpSparseMat,Kernel_Sa,tmp_double, tmp_int, -1);
		 TSak.Clear();


	   if (config::info::printMatrices) {
			//SparseMatrix RT = cluster.domains[d].Kplus_R;
			//RT.ConvertDenseToCSR(1);

			std::ofstream osSa(Logging::prepareFile("Salfa"));
			osSa << Salfa;
			osSa.close();
     }


		 //domains[0].get_kernel_from_K(Salfa, Kernel_Sa);

//		 Salfa.printMatCSR2("Salfa.txt");

//		 SparseMatrix T1;
//		 T1.MatMat(G0t,'N', Kernel_Sa);
//		 SpyText(T1);
//		 T1.printMatCSR("BlMat");
//		 Kernel_Sa.printMatCSR2("H.txt");

		 char str000[128];
		 for (int d = 0; d < domains.size(); d++) {
			 SparseMatrix tR;
			 SEQ_VECTOR < eslocal > rows_inds (Kernel_Sa.cols);
			 for (int i = 0; i < Kernel_Sa.cols; i++)
				 rows_inds[i] = 1 + d * Kernel_Sa.cols + i;
			 tR.CreateMatFromRowsFromMatrix_NewSize(Kernel_Sa,rows_inds);
			 sprintf(str000,"%s%d%s","tr",d,".txt");
//			 tR.printMatCSR2(str000);
			 SparseMatrix TmpR;
			 domains[d].Kplus_R.ConvertDenseToCSR(0);
			 TmpR.MatMat( domains[d].Kplus_R, 'N', tR );
			 domains[d].Kplus_Rb = TmpR;
			 domains[d].Kplus_Rb.ConvertCSRToDense(0);
			 //SparseMatrix T2;
			 //T2.MatMat(domains[d].Prec, 'N', domains[d].Kplus_R);
		 }


	 } else {
		// Regularization of Salfa
		 TimeEvent reg_Sa_time("Salfa regularization "); reg_Sa_time.start();

		SparseMatrix Eye, N, Nt, NNt;
		eslocal dtmp = 0;
		if (DOFS_PER_NODE == 3) {
			dtmp = 6;
		}

		if (DOFS_PER_NODE == 1) {
			dtmp = 1;
		}

		Eye.CreateEye(dtmp); N.CreateEye(dtmp); Nt.CreateEye(dtmp);

		for (int i=0; i < domains.size()-1; i++) {
			N.MatAppend(Eye);
			Nt.MatAppend(Eye);
		}

		Nt.MatTranspose();
		NNt.MatMat(N, 'N', Nt);
		NNt.RemoveLower();

		double ro = Salfa.GetMeanOfDiagonalOfSymmetricMatrix();
		ro = 0.5 * ro;

		Salfa.MatAddInPlace(NNt,'N', ro);
		// End regularization of Salfa

		 reg_Sa_time.end(); reg_Sa_time.printStatMPI(); Sa_timing.addEvent(reg_Sa_time);
	 }


	 if (config::solver::SA_SOLVER == config::SA_SPARSE_on_CPU) {
     //#ifdef SPARSE_SA
		 TimeEvent fact_Sa_time("Salfa factorization "); fact_Sa_time.start();
		if (MPIrank == 0) Sa.msglvl = 1;
		Sa.ImportMatrix(Salfa);
		Sa.Factorization("salfa");
		if (MPIrank == 0) Sa.msglvl = 0;
		 fact_Sa_time.end(); fact_Sa_time.printStatMPI(); Sa_timing.addEvent(fact_Sa_time);
     //#else
	 }

	 if (config::solver::SA_SOLVER == config::SA_DENSE_on_CPU) {
		 TimeEvent factd_Sa_time("Salfa factorization - dense "); factd_Sa_time.start();

		Salfa.ConvertCSRToDense(1);

		Sa_dense_cpu.ImportMatrix(Salfa);
		Sa_dense_cpu.Factorization("Salfa - dense ");

		 factd_Sa_time.end(); factd_Sa_time.printStatMPI(); Sa_timing.addEvent(factd_Sa_time);
	 }

	 if (config::solver::SA_SOLVER == config::SA_DENSE_on_ACC) {
		 TimeEvent factd_Sa_time("Salfa factorization - dense "); factd_Sa_time.start();

//#if defined(SOLVER_CUDA)
		 Salfa.type = 'G';
//#endif
//#if defined(SOLVER_CUDA_7)
//	 Salfa.type = 'G';
//#endif

		Salfa.ConvertCSRToDense(1);
		Salfa.type = 'S';

		Sa_dense_acc.ImportMatrix(Salfa);
		Sa_dense_acc.Factorization("Salfa - dense ");

		factd_Sa_time.end(); factd_Sa_time.printStatMPI(); Sa_timing.addEvent(factd_Sa_time);
	 }

	Sa.m_Kplus_size = Salfa.cols;

//#endif // //#ifdef SPARSE_SA


	Sa_timing.totalTime.end(); Sa_timing.printStatsMPI();
	MKL_Set_Num_Threads(1);
}


void ClusterBase::Create_G1_perCluster() {

	SparseMatrix tmpM;

	if (USE_HFETI == 1) {
		//G1 = G1 + trans(B1 * domains[d].Kplus_R) for all domains
		//for (eslocal d = 0; d < domains.size(); d++)
		//{
		//	tmpM.MatMat( domains[d].B1, 'N', domains[d].Kplus_R);
		//	G1.MatAddInPlace( tmpM, 'N', 1.0 );
		//	tmpM.Clear();
		//}
		//G1.MatTranspose();

		//// OK - but sequential
		//for (eslocal d = 0; d < domains.size(); d++)					//HFETI
		//{
		//	tmpM.MatMat( domains[d].B1t, 'T', domains[d].Kplus_R);
		//	G1.MatAddInPlace( tmpM, 'N', 1.0 );
		//	tmpM.Clear();
		//}
		//G1.MatTranspose();




		//eslocal threads = 24;
		//eslocal n_domains = domains.size();
		//vector < SparseMatrix > tmp_Mat (threads);
		//cilk_for (eslocal t = 0; t < threads; t++ ) {
		//	for (eslocal i = t*(n_domains/threads+1); i < (t+1)*(n_domains/threads+1); i++ ) {
		//		if (i < n_domains) {
		//			SparseMatrix tmpM_l;
		//			//tmpM_l.MatMat( domains[i].Kplus_R, 'T', domains[i].B1t);
		//			tmpM_l.MatMat( domains[i].B1t, 'T', domains[i].Kplus_R);
		//			tmpM_l.MatTranspose();
		//			tmp_Mat[t].MatAddInPlace( tmpM_l, 'N', 1.0 );
		//		}
		//	}
		//}

		TimeEvent G1_1_time ("Create G1 per clust t. : MatMat+MatTrans ");
		G1_1_time.start();
		TimeEvent G1_1_mem  ("Create G1 per clust mem: MatMat+MatTrans ");
		G1_1_mem.startWithoutBarrier(GetProcessMemory_u());


		//SparseMatrix Rt;
		//SparseMatrix B;
		//
		//domains[0].Kplus_R.MatTranspose(Rt);
		//B = domains[0].B1;
		//Rt.ConvertCSRToDense(0);
		//
		//vector < vector < double > > tmpG (B.nnz, vector <double> (Rt.rows,0));
		//vector <eslocal > G_I_row_indices;
		//G_I_row_indices.resize(B.nnz);

		//eslocal indx = 0;
		//for (eslocal r = 0; r < Rt.rows; r++)
		//	tmpG[indx][r] += B.V_values[0] * Rt.dense_values[Rt.rows * (B.J_col_indices[0]-1) + r];
		//
		//G_I_row_indices[indx] = B.I_row_indices[0];

		//for (eslocal i = 1; i < B.I_row_indices.size(); i++) {
		//
		//	if (B.I_row_indices[i-1] != B.I_row_indices[i])
		//		indx++;

		//	for (eslocal r = 0; r < Rt.rows; r++)
		//		tmpG[indx][r] += B.V_values[i] * Rt.dense_values[Rt.rows * (B.J_col_indices[i]-1) + r];
		//
		//	G_I_row_indices[indx] = B.I_row_indices[i];
		//}
		//
		//G_I_row_indices.resize(indx+1);
		//tmpG.resize(indx+1);

		//vector <int>    G_I; G_I.reserve( tmpG.size() * tmpG[0].size());
		//vector <int>    G_J; G_J.reserve( tmpG.size() * tmpG[0].size());
		//vector <double> G_V; G_V.reserve( tmpG.size() * tmpG[0].size());

		//for (eslocal i = 0; i < tmpG.size(); i++) {
		//	for (eslocal j = 0; j < tmpG[i].size(); j++){
		//		if (tmpG[i][j] != 0) {
		//			G_I.push_back(G_I_row_indices[i]);
		//			G_J.push_back(j+1);
		//			G_V.push_back(tmpG[i][j]);
		//		}
		//	}
		//}

		//SparseMatrix Gcoo;
		//Gcoo.I_row_indices = G_J; //G_I;
		//Gcoo.J_col_indices = G_I; //G_J;
		//Gcoo.V_values      = G_V;
		//Gcoo.cols = B.rows; //R.cols;
		//Gcoo.rows = Rt.rows; //B.rows;
		//Gcoo.nnz  = G_I.size();
		//Gcoo.type = 'G';
		//
		//Gcoo.ConvertToCSRwithSort(1);




		//SparseMatrix Gtmp;
		//SparseMatrix Gtmpt;
		//
		//Gtmp.MatMat( domains[0].B1t, 'T', domains[0].Kplus_R);
		//Gtmp.MatTranspose(Gtmpt);
		//
		//Gtmp.ConvertToCOO(0);
		//Gtmpt.ConvertToCOO(0);

		int MPIrank;
		MPI_Comm_rank (MPI_COMM_WORLD, &MPIrank);

		PAR_VECTOR < SparseMatrix > tmp_Mat (domains.size());
		cilk_for (eslocal j = 0; j < tmp_Mat.size(); j++) {
			// V1
			//tmp_Mat[j].MatMat( domains[j].B1t, 'T', domains[j].Kplus_R);
			//tmp_Mat[j].MatTranspose();

			// V2 - not cool
			//tmp_Mat[j].MatMatSorted( domains[j].Kplus_R, 'T', domains[j].B1t); // - pozor mooooc pomale a MKL rutina zere mooooooc pameti

			// V3
			SparseMatrix Gcoo;
			if (domains[j].B1.nnz > 0) {

				SparseMatrix Rt;
				SparseMatrix B;

				if ( config::solver::REGULARIZATION == 0 ) {
					Rt = domains[j].Kplus_R;
					Rt.ConvertDenseToCSR(1);
					Rt.MatTranspose();
					//domains[j].Kplus_R.MatTranspose(Rt);
				} else {
					Rt = domains[j].Kplus_Rb;
					Rt.ConvertDenseToCSR(1);
					Rt.MatTranspose();
					//domains[j].Kplus_Rb.MatTranspose(Rt);
				}

				//Rt = domains[j].Kplus_R;
				//Rt.MatTranspose();

				Rt.ConvertCSRToDense(1);
				B = domains[j].B1;

				SEQ_VECTOR < SEQ_VECTOR < double > > tmpG (B.nnz, SEQ_VECTOR <double> (Rt.rows,0));
				SEQ_VECTOR <eslocal > G_I_row_indices;
				G_I_row_indices.resize(B.nnz);

				eslocal indx = 0;
				for (eslocal r = 0; r < Rt.rows; r++)
					tmpG[indx][r] += B.V_values[0] * Rt.dense_values[Rt.rows * (B.J_col_indices[0]-1) + r];

				G_I_row_indices[indx] = B.I_row_indices[0];

				for (eslocal i = 1; i < B.I_row_indices.size(); i++) {

					if (B.I_row_indices[i-1] != B.I_row_indices[i])
						indx++;

					for (eslocal r = 0; r < Rt.rows; r++)
						tmpG[indx][r] += B.V_values[i] * Rt.dense_values[Rt.rows * (B.J_col_indices[i]-1) + r];

					G_I_row_indices[indx] = B.I_row_indices[i];
				}

				G_I_row_indices.resize(indx+1);
				tmpG.resize(indx+1);

				SEQ_VECTOR <eslocal>    G_I; G_I.reserve( tmpG.size() * tmpG[0].size());
				SEQ_VECTOR <eslocal>    G_J; G_J.reserve( tmpG.size() * tmpG[0].size());
				SEQ_VECTOR <double> G_V; G_V.reserve( tmpG.size() * tmpG[0].size());

				for (eslocal i = 0; i < tmpG.size(); i++) {
					for (eslocal j = 0; j < tmpG[i].size(); j++){
						if (tmpG[i][j] != 0) {
							G_I.push_back(G_I_row_indices[i]);
							G_J.push_back(j+1);
							G_V.push_back(tmpG[i][j]);
						}
					}
				}




				Gcoo.I_row_indices = G_J; //G_I;
				Gcoo.J_col_indices = G_I; //G_J;
				Gcoo.V_values      = G_V;
				Gcoo.cols = B.rows; //R.cols;
				Gcoo.rows = Rt.rows; //B.rows;
				Gcoo.nnz  = G_I.size();
				Gcoo.type = 'G';

				Gcoo.ConvertToCSRwithSort(1);

			} else {
				Gcoo.cols = domains[j].B1.rows;
				if ( config::solver::REGULARIZATION == 0 )
					Gcoo.rows = domains[j].Kplus_R.rows;
				else
					Gcoo.rows = domains[j].Kplus_Rb.rows;
				Gcoo.I_row_indices.resize(1,0);
				Gcoo.J_col_indices.resize(1,0);
				Gcoo.V_values.resize(0,0);

				Gcoo.nnz  = 0;
				Gcoo.type = 'G';

				Gcoo.ConvertToCSRwithSort(1);
			}

			tmp_Mat[j] = Gcoo;

			// END - V3

			//if (MPIrank == 0)
			//	cout << j << endl;

		}

		G1_1_time.end();
		G1_1_time.printStatMPI();
		//G1_1_time.printLastStatMPIPerNode();
		G1_1_mem.endWithoutBarrier(GetProcessMemory_u());
		G1_1_mem.printStatMPI();
		//G1_1_mem.printLastStatMPIPerNode();

		TimeEvent G1_2_time ("Create G1 per clust t. : Par.red.+MatAdd ");
		G1_2_time.start();
		TimeEvent G1_2_mem  ("Create G1 per clust mem: Par.red.+MatAdd ");
		G1_2_mem.startWithoutBarrier(GetProcessMemory_u());

		for (eslocal j = 1; j < tmp_Mat.size(); j = j * 2 ) {
			cilk_for (eslocal i = 0; i < tmp_Mat.size(); i = i + 2*j) {
				if ( i+j < tmp_Mat.size()) {
					tmp_Mat[i    ].MatAddInPlace( tmp_Mat[i + j], 'N', 1.0 ); //  MFETI - MatAppend(tmp_Mat[i + j]);
					tmp_Mat[i + j].Clear();
				}
			}
		}

		G1_2_time.end();
		G1_2_time.printStatMPI();
		//G1_2_time.printLastStatMPIPerNode();

		G1_2_mem.endWithoutBarrier(GetProcessMemory_u());
		G1_2_mem.printStatMPI();
		//G1_2_mem.printLastStatMPIPerNode();

		G1 = tmp_Mat[0];
		tmp_Mat[0].Clear();
		//G1.MatTranspose();


	} else {

		PAR_VECTOR < SparseMatrix > tmp_Mat (domains.size());
		cilk_for (eslocal j = 0; j < domains.size(); j++) {

			// V1
			//tmp_Mat[j].MatMat( domains[j].B1t, 'T', domains[j].Kplus_R);
			//tmp_Mat[j].MatTranspose();

			// V2
			//tmp_Mat[j].MatMat( domains[j].Kplus_R, 'T', domains[j].B1t);

			// V3
			SparseMatrix Rt;
			SparseMatrix B;

			if ( config::solver::REGULARIZATION == 0 ) {
				Rt = domains[j].Kplus_R;
				Rt.ConvertDenseToCSR(1);
				Rt.MatTranspose();
				//domains[j].Kplus_R.MatTranspose(Rt);
			} else {
				Rt = domains[j].Kplus_Rb;
				Rt.ConvertDenseToCSR(1);
				Rt.MatTranspose();
				//domains[j].Kplus_Rb.MatTranspose(Rt);
			}

			Rt.ConvertCSRToDense(1);
			B = domains[j].B1;

			SEQ_VECTOR < SEQ_VECTOR < double > > tmpG (B.nnz, SEQ_VECTOR <double> (Rt.rows,0));
			SEQ_VECTOR <eslocal > G_I_row_indices;
			G_I_row_indices.resize(B.nnz);

			eslocal indx = 0;
			for (eslocal r = 0; r < Rt.rows; r++)
				tmpG[indx][r] += B.V_values[0] * Rt.dense_values[Rt.rows * (B.J_col_indices[0]-1) + r];

			G_I_row_indices[indx] = B.I_row_indices[0];

			for (eslocal i = 1; i < B.I_row_indices.size(); i++) {

				if (B.I_row_indices[i-1] != B.I_row_indices[i])
					indx++;

				for (eslocal r = 0; r < Rt.rows; r++)
					tmpG[indx][r] += B.V_values[i] * Rt.dense_values[Rt.rows * (B.J_col_indices[i]-1) + r];

				G_I_row_indices[indx] = B.I_row_indices[i];
			}

			G_I_row_indices.resize(indx+1);
			tmpG.resize(indx+1);

			SEQ_VECTOR <eslocal>    G_I; G_I.reserve( tmpG.size() * tmpG[0].size());
			SEQ_VECTOR <eslocal>    G_J; G_J.reserve( tmpG.size() * tmpG[0].size());
			SEQ_VECTOR <double> G_V; G_V.reserve( tmpG.size() * tmpG[0].size());

			for (eslocal i = 0; i < tmpG.size(); i++) {
				for (eslocal j = 0; j < tmpG[i].size(); j++){
					if (tmpG[i][j] != 0) {
						G_I.push_back(G_I_row_indices[i]);
						G_J.push_back(j+1);
						G_V.push_back(tmpG[i][j]);
					}
				}
			}

			SparseMatrix Gcoo;
			Gcoo.I_row_indices = G_J; //G_I;
			Gcoo.J_col_indices = G_I; //G_J;
			Gcoo.V_values      = G_V;
			Gcoo.cols = B.rows; //R.cols;
			Gcoo.rows = Rt.rows; //B.rows;
			Gcoo.nnz  = G_I.size();
			Gcoo.type = 'G';

			Gcoo.ConvertToCSRwithSort(1);

			tmp_Mat[j] = Gcoo;

		}

		for (eslocal j = 1; j < tmp_Mat.size(); j = j * 2 ) {
			cilk_for (eslocal i = 0; i < tmp_Mat.size(); i = i + 2*j) {
				if ( i+j < tmp_Mat.size()) {
					tmp_Mat[i    ].MatAppend(tmp_Mat[i + j]);
					tmp_Mat[i + j].Clear();
				}
			}
		}
		G1.MatAppend(tmp_Mat[0]);
		//for (eslocal d = 0; d < domains.size(); d++)					//MFETI
		//{
		//	//tmpM.MatMat( domains[d].B1, 'N', domains[d].Kplus_R);
		//	tmpM.MatMat( domains[d].B1t, 'T', domains[d].Kplus_R);
		//	tmpM.MatTranspose();
		//	G1.MatAppend(tmpM);
		//	tmpM.Clear();
		//}
	}

	// for both MFETI and HFETI
	for (eslocal i = 0; i < G1.CSR_V_values.size(); i++)
		G1.CSR_V_values[i] = -1.0 * G1.CSR_V_values[i];


//	std::stringstream ss;
//	ss << "G" << ".txt";
//	std::ofstream os(ss.str().c_str());
//	os << G1;
//	os.close();

}

void ClusterBase::Compress_G1() {

	G1.ConvertToCOO( 1 );

	cilk_for (eslocal j = 0; j < G1.J_col_indices.size(); j++ )
		G1.J_col_indices[j] = _my_lamdas_map_indices[ G1.J_col_indices[j] -1 ] + 1;  // numbering from 1 in matrix

	G1.cols = my_lamdas_indices.size();
	G1.ConvertToCSRwithSort( 1 );

	G1_comp.CSR_I_row_indices.swap( G1.CSR_I_row_indices );
	G1_comp.CSR_J_col_indices.swap( G1.CSR_J_col_indices );
	G1_comp.CSR_V_values     .swap( G1.CSR_V_values		 );

	G1_comp.rows = G1.rows;
	G1_comp.cols = G1.cols;
	G1_comp.nnz  = G1.nnz;
	G1_comp.type = G1.type;

	G1.Clear();

}

void ClusterBase::CreateVec_d_perCluster( SEQ_VECTOR<SEQ_VECTOR <double> > & f ) {

	eslocal size_d = domains[0].Kplus_R.cols; // because transpose of R

	if ( USE_HFETI == 1 )
		vec_d.resize( size_d );
	else
		vec_d.resize( domains.size() * size_d );	// MFETI

	if ( USE_HFETI == 1) {
		for (eslocal d = 0; d < domains.size(); d++) {
			if ( config::solver::REGULARIZATION == 0 ) {
				domains[d].Kplus_R .DenseMatVec(f[d], vec_d, 'T', 0, 0         , 1.0 );
			} else {
				domains[d].Kplus_Rb.DenseMatVec(f[d], vec_d, 'T', 0, 0         , 1.0 );
			}
		}
	} else {
		for (eslocal d = 0; d < domains.size(); d++) {											// MFETI
			if ( config::solver::REGULARIZATION == 0 ) {
				domains[d].Kplus_R.DenseMatVec(f[d], vec_d, 'T', 0, d * size_d, 0.0 );				// MFETI
			} else {
				domains[d].Kplus_Rb.DenseMatVec(f[d], vec_d, 'T', 0, d * size_d, 0.0 );				// MFETI
			}
		}
	}

	for (eslocal i = 0; i < vec_d.size(); i++)
		vec_d[i] = (-1.0) *  vec_d[i];


//	std::stringstream ss;
//	ss.precision(40);
//	ss << "vec_d" << ".txt";
//	std::ofstream os(ss.str().c_str());
//	os.precision(40);
//	os << vec_d;
//	os.close();

}


void ClusterBase::CreateVec_b_perCluster( SEQ_VECTOR<SEQ_VECTOR <double> > & f )  {

	SEQ_VECTOR<SEQ_VECTOR<double> > x_prim_cluster ( domains.size() );
	cilk_for (eslocal d = 0; d < domains.size(); d++) {
		x_prim_cluster[d] = f[d];
	}

	if (USE_HFETI == 0) {
		cilk_for (eslocal d = 0; d < domains.size(); d++) {				// MFETI
			domains[d].multKplusLocal( x_prim_cluster[d] );
		}
	} else {
		multKplusGlobal_l( x_prim_cluster );						// HFETI
	}

//	// pro ukoncovani v primaru - vypocet up0
//	cilk_for (eslocal d = 0; d < domains.size(); d++) {
//		domains[d].up0 = x_prim_cluster[d];
//	}

	vec_b_compressed.resize(my_lamdas_indices.size(), 0.0);

	SEQ_VECTOR < double > y_out_tmp (domains[0].B1_comp_dom.rows);

	for (eslocal d = 0; d < domains.size(); d++) {
		y_out_tmp.resize( domains[d].B1_comp_dom.rows );
		domains[d].B1_comp_dom.MatVec (x_prim_cluster[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

		for (eslocal i = 0; i < domains[d].lambda_map_sub_local.size(); i++)
			vec_b_compressed[ domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i] - domains[d].vec_c[i];

		//for (eslocal i = 0; lambda_map_sub_local.size(); i++)
		//	vec_b_compressed[i] -= domains[d].vec_c[i];
	}

}







void ClusterBase::compress_lambda_vector  ( SEQ_VECTOR <double> & decompressed_vec_lambda )
{
	//compress vector for CG in main loop
	for (eslocal i = 0; i < my_lamdas_indices.size(); i++)
		decompressed_vec_lambda[i] = decompressed_vec_lambda[my_lamdas_indices[i]];

	decompressed_vec_lambda.resize(my_lamdas_indices.size());
}

void ClusterBase::decompress_lambda_vector( SEQ_VECTOR <double> &   compressed_vec_lambda )
{
	SEQ_VECTOR <double> decompressed_vec_lambda (domains[0].B1.rows,0);

	for (eslocal i = 0; i < my_lamdas_indices.size(); i++)
		decompressed_vec_lambda[my_lamdas_indices[i]] = compressed_vec_lambda[i];

	compressed_vec_lambda = decompressed_vec_lambda;
}


void ClusterBase::B1_comp_MatVecSum( SEQ_VECTOR < SEQ_VECTOR <double> > & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose ) {

	ESINFO(ERROR) << " B1_comp_MatVecSum - not implemented ";

	exit(0);

	//if (domains.size() <= 3) {

//		domains[0].B1_comp.MatVec (x_in[0], y_out, T_for_transpose_N_for_non_transpose, 0, 0, 0.0); // first vector overwrites temp vector
//		for (eslocal d = 1; d < domains.size(); d++) // reduction
//			domains[d].B1_comp.MatVec (x_in[d], y_out, T_for_transpose_N_for_non_transpose, 0, 0, 1.0); // will add (summation per elements) all partial results into y_out

	//} else {
	//
	//	cilk_for (eslocal d = 0; d < domains.size(); d++) {
	//		domains[d].B1_comp.MatVec (x_in[d], domains[d].compressed_tmp, T_for_transpose_N_for_non_transpose, 0, 0, 0.0);
	//	}

	//	for ( eslocal r = 2 * (domains.size() / 2);  r > 1;  r = r / 2 ) {
	//		cilk_for ( eslocal d = 0;  d < r;  d++ ) {
	//			if ( (r + d)  <  domains.size() ) {
	//				for (eslocal i = 0; i < domains[d].compressed_tmp.size(); i++) {
	//					domains[d].compressed_tmp[i] = domains[d].compressed_tmp[i] + domains[r + d].compressed_tmp[i];
	//				}
	//			}
	//		}
	//	}

	//	for (eslocal i = 0; i < domains[0].compressed_tmp.size(); i++)
	//		y_out[i] = domains[0].compressed_tmp[i] + domains[1].compressed_tmp[i];
	//
	//}

}

// **** END - CLUSTER CLASS ************************************************
// *******************************************************************







