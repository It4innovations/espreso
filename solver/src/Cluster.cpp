#include "Cluster.h"
#include <tbb/mutex.h>
//#define SPARSE_SA

// *******************************************************************
// **** CLUSTER CLASS ************************************************

Cluster::Cluster(eslocal cluster_index){

	cluster_global_index = cluster_index;

	cluster_time  = TimeEval ("Cluster Timing ");

	vec_fill_time = TimeEvent("Reseting vec_g0 and vec_e0");
	loop_1_1_time = TimeEvent("Loop 1: Kplus-sv, B0-mv, KpluR-mv");
	loop_1_2_time = TimeEvent("Loop 1: vec_e0 and vec_g0");

	clusCP_time   = TimeEvent("Cluster CP - F0,GO,Sa,G0t,F0 ");
	clus_F0_1_time= TimeEvent("F0 solve - 1st ");
	clus_F0_2_time= TimeEvent("F0 solve - 2nd ");
	clus_G0_time  = TimeEvent("G0  Mult ");
	clus_G0t_time = TimeEvent("G0t Mult ");
	clus_Sa_time  = TimeEvent("Sa solve ");

	loop_2_1_time = TimeEvent("Loop2: Kplus-sv, B0-mv, Kplus-mv");

	iter_cnt_comm = 0;
}



Cluster::Cluster() {

	cluster_time  = TimeEval("Cluster Timing ");

	vec_fill_time = TimeEvent("Reseting vec_g0 and vec_e0");
	loop_1_1_time = TimeEvent("Loop 1: Kplus-sv, B0-mv, KpluR-mv");
	loop_1_2_time = TimeEvent("Loop 1: vec_e0 and vec_g0");

	clusCP_time   = TimeEvent("Cluster CP - F0,GO,Sa,G0t,F0 ");
	clus_F0_1_time= TimeEvent("F0 solve - 1st ");
	clus_F0_2_time= TimeEvent("F0 solve - 2nd ");
	clus_G0_time  = TimeEvent("G0  Mult ");
	clus_G0t_time = TimeEvent("G0t Mult ");
	clus_Sa_time  = TimeEvent("Sa solve ");

	loop_2_1_time = TimeEvent("Loop2: Kplus-sv, B0-mv, Kplus-mv");

	iter_cnt_comm = 0;
}


void Cluster::ShowTiming()  {

	cluster_time.AddEvent(vec_fill_time);
	cluster_time.AddEvent(loop_1_1_time);
	cluster_time.AddEvent(loop_1_2_time);

	cluster_time.AddEvent(clus_F0_1_time);
	cluster_time.AddEvent(clus_G0_time);
	cluster_time.AddEvent(clus_Sa_time);
	cluster_time.AddEvent(clus_G0t_time);
	cluster_time.AddEvent(clus_F0_2_time);

	cluster_time.AddEvent(clusCP_time);
	cluster_time.AddEvent(loop_2_1_time);

	cluster_time.PrintStatsMPI();
}

void Cluster::SetDynamicParameters(double set_dynamic_timestep, double set_dynamic_beta, double set_dynamic_gama) {

	dynamic_timestep = set_dynamic_timestep;
	dynamic_beta     = set_dynamic_beta;
	dynamic_gama     = set_dynamic_gama;

	for (eslocal d = 0; d < domains.size(); d++)
		domains[d].SetDynamicParameters(dynamic_timestep, dynamic_beta, dynamic_gama);

}


void Cluster::InitClusterPC( eslocal * subdomains_global_indices, eslocal number_of_subdomains ) {

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

	}
	// *** END - Init all domains of the cluster ***************************************
}

void Cluster::SetClusterPC( SEQ_VECTOR <SEQ_VECTOR <eslocal> > & lambda_map_sub ) {


	map <eslocal,eslocal> my_lamdas_map_indices;


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


	if (MPIrank == 0) {
		cout << " ******************************************************************************************************************************* " << endl;
		cout << " *** setting vectors for lambda map ******************************************************************************************** " << endl;
		GetProcessMemoryStat_u ( ); GetMemoryStat_u( );
		cout << endl;
	}


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


	if (MPIrank == 0) {
		cout << " ******************************************************************************************************************************* " << endl;
		cout << " *** setting vectors for lambda communication  ********************************************************************************* " << endl;
		GetProcessMemoryStat_u ( ); GetMemoryStat_u( );
		cout << endl;
	}

	my_comm_lambdas_indices .resize(my_neighs.size());
	my_comm_lambdas			.resize(my_neighs.size());
	my_recv_lambdas			.resize(my_neighs.size());

	cilk_for (eslocal i = 0; i < my_neighs.size(); i++) {
		my_comm_lambdas_indices[i] = lambdas_per_subdomain[my_neighs[i]];
		my_comm_lambdas[i].resize(my_comm_lambdas_indices[i].size());
		my_recv_lambdas[i].resize(my_comm_lambdas_indices[i].size());
	}

	compressed_tmp    .resize( my_lamdas_indices.size(), 0 );
	//compressed_tmp2   .resize( my_lamdas_indices.size(), 0 );


	cilk_for (eslocal d = 0; d < domains.size(); d++ )
		if (USE_KINV == 1 )
			domains[d].compressed_tmp.resize( my_lamdas_indices.size(), 0);
		else
			domains[d].compressed_tmp.resize( 1, 0);


	// mapping/compression vector for cluster
	for (eslocal i = 0; i <my_lamdas_indices.size(); i++)
		my_lamdas_map_indices.insert(make_pair(my_lamdas_indices[i],i));

	// mapping/compression vector for domains
	cilk_for (eslocal i = 0; i < domains.size(); i++) {
		for (eslocal j = 0; j < domains[i].lambda_map_sub.size(); j++) {
			domains[i].my_lamdas_map_indices.insert(make_pair(domains[i].lambda_map_sub[j] ,j));
		}
	}

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



	//// *** Create a vector of communication pattern needed for AllReduceLambdas function *******
	my_comm_lambdas_indices_comp.resize(my_neighs.size());
	cilk_for (eslocal i = 0; i < my_neighs.size(); i++) {
		my_comm_lambdas_indices_comp[i].resize( lambdas_per_subdomain[my_neighs[i]].size() );
		for (eslocal j = 0; j < lambdas_per_subdomain[my_neighs[i]].size(); j++ )
			my_comm_lambdas_indices_comp[i][j] = my_lamdas_map_indices[lambdas_per_subdomain[my_neighs[i]][j]];
	}
	//// *** END - Create a vector of communication pattern needed for AllReduceLambdas function *


	//// *** Compression of Matrix B1 to work with compressed lambda vectors *****************
	if (MPIrank == 0) {
		cout << " ******************************************************************************************************************************* " << endl;
		cout << " *** B1 compression ************************************************************************************************************ " << endl;
		GetProcessMemoryStat_u ( ); GetMemoryStat_u( );
		cout << endl;
	}

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

		//domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);


		//************************

		//TODO: pozor vratit
		domains[i].B1.Clear();
		domains[i].B1t.Clear();
//		domains[i].B1_comp.Clear();
//		domains[i].B1t_comp.Clear();

		domains[i].my_lamdas_map_indices.clear();

	}


	//// *** END - Compression of Matrix B1 to work with compressed lambda vectors *************



	//// *** Compression of Matrix G1 to work with compressed lambda vectors *******************
	if (MPIrank == 0) {
		cout << " ******************************************************************************************************************************* " << endl;
		cout << " *** G1 compression ************************************************************************************************************ " << endl;
		GetProcessMemoryStat_u ( ); GetMemoryStat_u( );
		cout << endl;
	}

	if (USE_DYNAMIC == 0) {

		G1.ConvertToCOO( 1 );
		cilk_for (eslocal j = 0; j < G1.J_col_indices.size(); j++ )
			G1.J_col_indices[j] = my_lamdas_map_indices[ G1.J_col_indices[j] -1 ] + 1;  // numbering from 1 in matrix

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
	//// *** END - Compression of Matrix G1 to work with compressed lambda vectors ***************

	bool R_from_mesh = true;
	if ( esconfig::solver::REGULARIZATION == 0 )
  		R_from_mesh = true	;
  	else
  		R_from_mesh = false	;

	if (! R_from_mesh)
		_my_lamdas_map_indices = my_lamdas_map_indices;

	if (MPIrank == 0) {
		cout << " ******************************************************************************************************************************* " << endl;
		cout << " *** setting vectors end ******************************************************************************************************* " << endl;
		GetProcessMemoryStat_u ( ); GetMemoryStat_u( );
		cout << endl;
	}

}

void Cluster::SetClusterHFETI (bool R_from_mesh) {
	// *** Create Matrices and allocate vectors for Hybrid FETI **************************
	if (USE_HFETI == 1) {

		TimeEval HFETI_prec_timing (" HFETI - preprocessing timing");
		HFETI_prec_timing.totalTime.AddStart(omp_get_wtime());

		int MPIrank;
		MPI_Comm_rank (MPI_COMM_WORLD, &MPIrank);
		if (MPIrank == 0 ) { cout << endl << "*** HFETI - Preprocessing start ************************************************************************************************ " << endl << endl; }


		TimeEvent B0_time("Compress B0 per cluster");
		B0_time.AddStart(omp_get_wtime());

		CompressB0();

		B0_time.AddEnd(omp_get_wtime());
		B0_time.PrintStatMPI(0.0);
		HFETI_prec_timing.AddEvent(B0_time);


		TimeEvent G0_time("Create G0 per cluster");
		G0_time.AddStart(omp_get_wtime());

		CreateG0();
		vec_g0.resize(G0.cols);
		vec_e0.resize(G0.rows);

		G0_time.AddEnd(omp_get_wtime());
		G0_time.PrintStatMPI(0.0);
		HFETI_prec_timing.AddEvent(G0_time);


		TimeEvent F0_time("Create F0 per cluster");
		F0_time.AddStart(omp_get_wtime());

		CreateF0();
		vec_lambda.resize(F0.m_Kplus_size);

		F0_time.AddEnd(omp_get_wtime());
		HFETI_prec_timing.AddEvent(F0_time);


		TimeEvent Sa_time("Create Salfa per cluster");
		Sa_time.AddStart(omp_get_wtime());

		CreateSa();
		vec_alfa.resize(Sa.m_Kplus_size);

		Sa_time.AddEnd(omp_get_wtime());
		HFETI_prec_timing.AddEvent(Sa_time);

		HFETI_prec_timing.totalTime.AddEnd(omp_get_wtime());
		HFETI_prec_timing.PrintStatsMPI();

		if (! R_from_mesh) {

			Create_G1_perCluster();

			if (USE_DYNAMIC == 0) {

				G1.ConvertToCOO( 1 );
				cilk_for (int j = 0; j < G1.J_col_indices.size(); j++ )
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

			}
		}


	}
	// *** END - Create Matrices for Hybrid FETI *****************************************
}

void Cluster::SetClusterPC_AfterKplus () {

	//// *** Alocate temporarly vectors for Temporary vectors for Apply_A function *********
	//// *** - temporary vectors for work primal domain size *******************************
	x_prim_cluster1.resize( domains.size() );
	x_prim_cluster2.resize( domains.size() );

	for (eslocal d = 0; d < domains.size(); d++) {
		x_prim_cluster1[d].resize( domains[d].domain_prim_size );
		x_prim_cluster2[d].resize( domains[d].domain_prim_size );
	}
	//// *** END - Alocate temporarly vectors for Temporary vectors for Apply_A function ***

	//// *** Prepare the initial right hand side in dual *************************************
	//CreateVec_b_perCluster();

}


void Cluster::multKplusGlobal(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, SEQ_VECTOR<eslocal> & cluster_map_vec) {

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

void Cluster::multKplusGlobal_l(SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {

	mkl_set_num_threads(1);

	cluster_time.totalTime.AddStart();

	vec_fill_time.AddStart();
	fill(vec_g0.begin(), vec_g0.end(), 0); // reset entire vector to 0
	//fill(vec_e0.begin(), vec_e0.end(), 0); // reset entire vector to 0
	vec_fill_time.AddEnd();

	// loop over domains in the cluster
	loop_1_1_time.AddStart();
	cilk_for (eslocal d = 0; d < domains.size(); d++)
	{
		domains[d].B0Kplus_comp.DenseMatVec(x_in[d], tm2[d]);			// g0 - with comp B0Kplus
		domains[d].Kplus_R.DenseMatVec(x_in[d], tm3[d], 'T');			// e0
	}
	loop_1_1_time.AddEnd();

	loop_1_2_time.AddStart();
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
	loop_1_2_time.AddEnd();

	mkl_set_num_threads(PAR_NUM_THREADS);
	clusCP_time.AddStart();


//	for (int i = 0; i < vec_g0.size(); i++)
//	printf (       "Test probe 1: %d norm = %1.30f \n", i, vec_g0[i] );

	clus_F0_1_time.AddStart();
	F0.Solve(vec_g0, tm1[0], 0, 0);
	clus_F0_1_time.AddEnd();

//	for (int i = 0; i < tm1[0].size(); i++)
//	printf (       "Test probe 2: %d norm = %1.30f \n", i, tm1[0][i] );

	clus_G0_time.AddStart();
	G0.MatVec(tm1[0], tm2[0], 'N');
	clus_G0_time.AddEnd();

//	for (int i = 0; i < tm1[0].size(); i++)
//	printf (       "Test probe 3: %d norm = %1.30f \n", i, tm1[0][i] );

	cilk_for (eslocal i = 0; i < vec_e0.size(); i++)
		tm2[0][i] = tm2[0][i] - vec_e0[i];
	//cblas_daxpy(vec_e0.size(), -1.0, &vec_e0[0], 1, &tm2[0][0], 1);

	 clus_Sa_time.AddStart();
#ifdef SPARSE_SA
	 Sa.Solve(tm2[0], vec_alfa,0,0);
#else
	char U = 'U';
	eslocal nrhs = 1;
	eslocal info = 0;
	vec_alfa = tm2[0];
	dsptrs( &U, &SaMat.rows, &nrhs, &SaMat.dense_values[0], &SaMat.ipiv[0], &vec_alfa[0], &SaMat.rows, &info );
#endif
	 clus_Sa_time.AddEnd();

//		for (int i = 0; i < vec_alfa.size(); i++)
//		printf (       "Test probe 4: %d norm = %1.30f \n", i, vec_alfa[i] );

	 clus_G0t_time.AddStart();
	G0.MatVec(vec_alfa, tm1[0], 'T'); 	// lambda
	 clus_G0t_time.AddEnd();

//		for (int i = 0; i < tm1[0].size(); i++)
//		printf (       "Test probe 5: %d norm = %1.30f \n", i, tm1[0][i] );

	cilk_for (eslocal i = 0; i < vec_g0.size(); i++)
		tm1[0][i] = vec_g0[i] - tm1[0][i];


	clus_F0_2_time.AddStart();
	F0.Solve(tm1[0], vec_lambda,0,0);
	clus_F0_2_time.AddEnd();

	clusCP_time.AddEnd();

//	for (int i = 0; i < vec_lambda.size(); i++)
//	printf (       "Test probe 6: %d norm = %1.30f \n", i, vec_lambda[i] );

	// Kplus_x
	mkl_set_num_threads(1);
	loop_2_1_time.AddStart();
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

		domains[d].multKplusLocal(tm1[d] , tm2[d], 0, 0);

		eslocal e0_start	=  d	* domains[d].Kplus_R.cols;
		eslocal e0_end		= (d+1) * domains[d].Kplus_R.cols;

		domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);

		for (eslocal i = 0; i < domain_size; i++)
			x_in[d][i] = tm2[d][i] + tm3[d][i];

	}
	loop_2_1_time.AddEnd();

	cluster_time.totalTime.AddEnd();
}

void Cluster::multKplusGlobal_Kinv( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in ) {

	mkl_set_num_threads(1);
	cluster_time.totalTime.AddStart();

	vec_fill_time.AddStart();
	fill(vec_g0.begin(), vec_g0.end(), 0); // reset entire vector to 0
	//fill(vec_e0.begin(), vec_e0.end(), 0); // reset entire vector to 0
	vec_fill_time.AddEnd();

	// loop over domains in the cluster
	loop_1_1_time.AddStart();
	cilk_for (eslocal d = 0; d < domains.size(); d++)
	{
   		domains[d].B0Kplus_comp.DenseMatVec(x_in[d], tm2[d]);			// g0 - with comp B0Kplus
		domains[d].Kplus_R.     DenseMatVec(x_in[d], tm3[d], 'T');		// e0
	}
	loop_1_1_time.AddEnd();

	loop_1_2_time.AddStart();
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
	loop_1_2_time.AddEnd();

	mkl_set_num_threads(PAR_NUM_THREADS);
	clusCP_time.AddStart();

	clus_F0_1_time.AddStart();
	F0.Solve(vec_g0, tm1[0], 0, 0);
	clus_F0_1_time.AddEnd();

	clus_G0_time.AddStart();
	G0.MatVec(tm1[0], tm2[0], 'N');
	clus_G0_time.AddEnd();

	cilk_for (eslocal i = 0; i < vec_e0.size(); i++)
		tm2[0][i] = tm2[0][i] - vec_e0[i];
	//cblas_daxpy(vec_e0.size(), -1.0, &vec_e0[0], 1, &tm2[0][0], 1);

	 clus_Sa_time.AddStart();
 #ifdef SPARSE_SA
    Sa.Solve(tm2[0], vec_alfa,0,0);
 #else
    char U = 'U';
    eslocal nrhs = 1;
    eslocal info = 0;
    vec_alfa = tm2[0];
    dsptrs( &U, &SaMat.rows, &nrhs, &SaMat.dense_values[0], &SaMat.ipiv[0], &vec_alfa[0], &SaMat.rows, &info );
 #endif
     clus_Sa_time.AddEnd();

	 clus_G0t_time.AddStart();
	G0.MatVec(vec_alfa, tm1[0], 'T'); 	// lambda
	 clus_G0t_time.AddEnd();

	cilk_for (eslocal i = 0; i < vec_g0.size(); i++)
		tm1[0][i] = vec_g0[i] - tm1[0][i];


	clus_F0_2_time.AddStart();
	F0.Solve(tm1[0], vec_lambda,0,0);
	clus_F0_2_time.AddEnd();

	clusCP_time.AddEnd();


	// Kplus_x
	mkl_set_num_threads(1);
	loop_2_1_time.AddStart();
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
	loop_2_1_time.AddEnd();

	cluster_time.totalTime.AddEnd();
}

void Cluster::multKplusGlobal_Kinv_2( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in ) {

	mkl_set_num_threads(1);
	cluster_time.totalTime.AddStart();

	vec_fill_time.AddStart();
	fill(vec_g0.begin(), vec_g0.end(), 0); // reset entire vector to 0
	//fill(vec_e0.begin(), vec_e0.end(), 0); // reset entire vector to 0
	vec_fill_time.AddEnd();

	// loop over domains in the cluster
	loop_1_1_time.AddStart();
	cilk_for (eslocal d = 0; d < domains.size(); d++)
	{
   		//domains[d].B0Kplus_comp.DenseMatVec(x_in[d], tm2[d]);				// g0 - with comp B0Kplus
		//domains[d].Kplus_R.     DenseMatVec(x_in[d], tm3[d], 'T');		// e0
		domains[d].B0KplusB1_comp .DenseMatVec(x_in[d], tm2[d], 'N');		// g0 - with comp B0Kplus
		domains[d].Kplus_R_B1_comp.DenseMatVec(x_in[d], tm3[d], 'N');		// e0
	}
	loop_1_1_time.AddEnd();

	loop_1_2_time.AddStart();
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
	loop_1_2_time.AddEnd();

	mkl_set_num_threads(PAR_NUM_THREADS);
	clusCP_time.AddStart();

	clus_F0_1_time.AddStart();
	F0.Solve(vec_g0, tm1[0], 0, 0);
	clus_F0_1_time.AddEnd();

	clus_G0_time.AddStart();
	G0.MatVec(tm1[0], tm2[0], 'N');
	clus_G0_time.AddEnd();

	cilk_for (eslocal i = 0; i < vec_e0.size(); i++)
		tm2[0][i] = tm2[0][i] - vec_e0[i];
	//cblas_daxpy(vec_e0.size(), -1.0, &vec_e0[0], 1, &tm2[0][0], 1);

	clus_Sa_time.AddStart();
	Sa.Solve(tm2[0], vec_alfa,0,0);
	clus_Sa_time.AddEnd();

	clus_G0t_time.AddStart();
	G0.MatVec(vec_alfa, tm1[0], 'T'); 	// lambda
	clus_G0t_time.AddEnd();

	cilk_for (eslocal i = 0; i < vec_g0.size(); i++)
		tm1[0][i] = vec_g0[i] - tm1[0][i];


	clus_F0_2_time.AddStart();
	F0.Solve(tm1[0], vec_lambda,0,0);
	clus_F0_2_time.AddEnd();

	clusCP_time.AddEnd();


	// Kplus_x
	mkl_set_num_threads(1);
	loop_2_1_time.AddStart();
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
	loop_2_1_time.AddEnd();

	cluster_time.totalTime.AddEnd();
}


////backup March 31 2015
//void Cluster::multKplusGlobal_l(SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {
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
//	cluster_time.totalTime.AddStart();
//
//	vec_fill_time.AddStart();
//	fill(vec_g0.begin(), vec_g0.end(), 0); // reset entire vector to 0
//	//fill(vec_e0.begin(), vec_e0.end(), 0); // reset entire vector to 0
//	vec_fill_time.AddEnd();
//
//	// loop over domains in the cluster
//	loop_1_1_time.AddStart();
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
//	loop_1_1_time.AddEnd();
//
//	loop_1_2_time.AddStart();
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
//	loop_1_2_time.AddEnd();
//
//	mkl_set_num_threads(24); //pozor
//	clusCP_time.AddStart();
//
//	clus_F0_1_time.AddStart();
//	F0.Solve(vec_g0, tm1[0], 0, 0);
//	clus_F0_1_time.AddEnd();
//
//	clus_G0_time.AddStart();
//	G0.MatVec(tm1[0], tm2[0], 'N');
//	clus_G0_time.AddEnd();
//
//	cilk_for (eslocal i = 0; i < vec_e0.size(); i++)
//		tm2[0][i] = tm2[0][i] - vec_e0[i];
//	//cblas_daxpy(vec_e0.size(), -1.0, &vec_e0[0], 1, &tm2[0][0], 1);
//
//	clus_Sa_time.AddStart();
//	Sa.Solve(tm2[0], vec_alfa,0,0);
//	clus_Sa_time.AddEnd();
//
//	clus_G0t_time.AddStart();
//	G0.MatVec(vec_alfa, tm1[0], 'T'); 	// lambda
//	clus_G0t_time.AddEnd();
//
//	cilk_for (eslocal i = 0; i < vec_g0.size(); i++)
//		tm1[0][i] = vec_g0[i] - tm1[0][i];
//
//
//	clus_F0_2_time.AddStart();
//	F0.Solve(tm1[0], vec_lambda,0,0);
//	clus_F0_2_time.AddEnd();
//
//	clusCP_time.AddEnd();
//
//
//	// Kplus_x
//	mkl_set_num_threads(1);
//	loop_2_1_time.AddStart();
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
//	loop_2_1_time.AddEnd();
//
//	cluster_time.totalTime.AddEnd();
//}


void Cluster::CompressB0() {

	cilk_for (eslocal d = 0; d < domains.size(); d++) {

		domains[d].B0.MatTranspose(domains[d].B0t);

		domains[d].B0_comp = domains[d].B0;
		domains[d].B0_comp.ConvertToCOO(1);
		domains[d].B0_comp_map_vec = domains[d].B0_comp.I_row_indices;

		for (eslocal i = 0; i < domains[d].B0_comp.I_row_indices.size(); i++)
			domains[d].B0_comp.I_row_indices[i] = i + 1;

		domains[d].B0_comp.rows = domains[d].B0_comp.I_row_indices.size();
		domains[d].B0_comp.ConvertToCSR(1);
		domains[d].B0_comp.MatTranspose(domains[d].B0t_comp);

		//domains[d].Kplus_R.ConvertCSRToDense(0); // TODO: - keep CSR data

	}

}

void Cluster::CreateG0() {

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

void Cluster::CreateF0() {

	 TimeEval F0_timing (" HFETI - F0 preprocessing timing");
	 F0_timing.totalTime.AddStart(omp_get_wtime());

	mkl_set_num_threads(1);

	int MPIrank; MPI_Comm_rank (MPI_COMM_WORLD, &MPIrank);

	SEQ_VECTOR <SparseMatrix> tmpF0v (domains.size());  
		
	if (MPIrank == 0 ) {cout << "HFETI - Create F0 - domain : " << endl; };
		
	 TimeEvent solve_F0_time("B0 compression; F0 multiple RHS solve");
	 solve_F0_time.AddStart(omp_get_wtime());

	cilk_for (eslocal d = 0; d < domains.size(); d++) {

		if (MPIrank == 0 && d == 0)
			domains[d].Kplus.msglvl=0;
		else
			domains[d].Kplus.msglvl=0;

		domains[d].Kplus.SolveMat_Dense(domains[d].B0t_comp, domains[d].B0Kplus_comp);
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
		if (MPIrank == 0 ) cout << "."; //{cout << d << " "; };
	}	

	if (MPIrank == 0 ) {cout << endl; };

	 solve_F0_time.AddEnd(omp_get_wtime());
	 solve_F0_time.PrintStatMPI(0.0);
	 F0_timing.AddEvent(solve_F0_time);
	
	if (MPIrank == 0 ) {cout << endl; };

	 TimeEvent reduction_F0_time("F0 reduction time");
	 reduction_F0_time.AddStart(omp_get_wtime());

	for (eslocal j = 1; j < tmpF0v.size(); j = j * 2 ) {
		cilk_for (eslocal i = 0; i < tmpF0v.size(); i = i + 2*j) {
			if ( i+j < tmpF0v.size()) {
				tmpF0v[i    ].MatAddInPlace( tmpF0v[i + j], 'N', 1.0 ); 
				tmpF0v[i + j].Clear();
			}
		}
	} 
	F0_Mat = tmpF0v[0]; 

	 reduction_F0_time.AddEnd(omp_get_wtime()); reduction_F0_time.PrintStatMPI(0.0); F0_timing.AddEvent(reduction_F0_time);


	 TimeEvent fact_F0_time("B0 Kplus Factorization "); fact_F0_time.AddStart(omp_get_wtime());

	mkl_set_num_threads(PAR_NUM_THREADS);
	F0_Mat.RemoveLower();
	F0.ImportMatrix(F0_Mat); 

	F0_fast.ImportMatrix(F0_Mat);

	//F0_Mat.Clear();
	F0.SetThreaded();
	F0.Factorization();

	mkl_set_num_threads(1);

	if (MPIrank == 0) F0.msglvl = 0;

	 fact_F0_time.AddEnd(omp_get_wtime()); fact_F0_time.PrintStatMPI(0.0); F0_timing.AddEvent(fact_F0_time);

	F0_timing.totalTime.AddEnd(omp_get_wtime());
	F0_timing.PrintStatsMPI();

	// *** POZOR **************************************************************
	cilk_for (eslocal d = 0; d<domains.size(); d++) {
		domains[d].B0.Clear(); 
		domains[d].B0t.Clear(); 
	}
};

void Cluster::CreateSa() {

	bool PARDISO_SC = true;
	bool get_kernel_from_mesh;
	
	if ( esconfig::solver::REGULARIZATION == 0 )
  		get_kernel_from_mesh = true	;
  	else
  		get_kernel_from_mesh = false	;
	
	
	MKL_Set_Num_Threads(PAR_NUM_THREADS);
	int MPIrank; MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

	SparseMatrix Salfa; SparseMatrix tmpM;

	TimeEval Sa_timing (" HFETI - Salfa preprocessing timing"); Sa_timing.totalTime.AddStart(omp_get_wtime());

	 TimeEvent G0trans_Sa_time("G0 transpose"); G0trans_Sa_time.AddStart(omp_get_wtime());
	SparseMatrix G0t; 
	G0.MatTranspose(G0t);
	 G0trans_Sa_time.AddEnd(omp_get_wtime()); G0trans_Sa_time.PrintStatMPI(0.0); Sa_timing.AddEvent(G0trans_Sa_time); 

	 TimeEvent G0solve_Sa_time("SolveMatF with G0t as RHS"); G0solve_Sa_time.AddStart(omp_get_wtime());
	if (!PARDISO_SC) {
		if (MPIrank == 0) F0_fast.msglvl = 1;
		//SolaveMatF is obsolete - use Schur Complement Instead
		F0_fast.SolveMatF(G0t,tmpM, true);
		if (MPIrank == 0) F0_fast.msglvl = 0;
	} else {
		SparseSolver tmpsps;
		if (MPIrank == 0) tmpsps.msglvl = 1;
		tmpsps.Create_SC_w_Mat( F0_Mat, G0t, Salfa, true, 0 );
        Salfa.ConvertDenseToCSR(1);
        Salfa.RemoveLower();
		if (MPIrank == 0) tmpsps.msglvl = 0;
	}
	F0_Mat.Clear();
	 G0solve_Sa_time.AddEnd(omp_get_wtime()); G0solve_Sa_time.PrintStatMPI(0.0); Sa_timing.AddEvent(G0solve_Sa_time); 

	 TimeEvent SaMatMat_Sa_time("Salfa = MatMat G0 * solution "); SaMatMat_Sa_time.AddStart(omp_get_wtime());
	if (!PARDISO_SC) {
		Salfa.MatMat(G0, 'N', tmpM);
		Salfa.RemoveLower();
		tmpM.Clear();
	}
	 SaMatMat_Sa_time.AddEnd(omp_get_wtime()); SaMatMat_Sa_time.PrintStatMPI(0.0); Sa_timing.AddEvent(SaMatMat_Sa_time); 
	
	 if (!get_kernel_from_mesh) {
		 SparseMatrix Kernel_Sa;
		 printf("Salfa\n");

		 SparseMatrix GGt;

		 GGt.MatMat(G0,'N',G0t);
		 GGt.RemoveLower();
		 GGt.get_kernel_from_K(GGt, Kernel_Sa);

		 SparseMatrix TSak;
		 GGt.get_kernel_from_K(Salfa,TSak);
		 TSak.Clear();

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
			 TmpR.MatMat( domains[d].Kplus_R, 'N', tR );
			 domains[d].Kplus_Rb = TmpR;
			 domains[d].Kplus_Rb.ConvertCSRToDense(0);
			 //SparseMatrix T2;
			 //T2.MatMat(domains[d].Prec, 'N', domains[d].Kplus_R);
		 }


	 } else {
		// Regularization of Salfa
		 TimeEvent reg_Sa_time("Salfa regularization "); reg_Sa_time.AddStart(omp_get_wtime());

		SparseMatrix Eye, N, Nt, NNt;
		Eye.CreateEye(6); N.CreateEye(6); Nt.CreateEye(6);

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

		 reg_Sa_time.AddEnd(omp_get_wtime()); reg_Sa_time.PrintStatMPI(0.0); Sa_timing.AddEvent(reg_Sa_time);
	 }

#ifdef SPARSE_SA
	 TimeEvent fact_Sa_time("Salfa factorization "); fact_Sa_time.AddStart(omp_get_wtime());
	if (MPIrank == 0) Sa.msglvl = 1;
	Sa.ImportMatrix(Salfa);
	Sa.Factorization();
	if (MPIrank == 0) Sa.msglvl = 0;
	 fact_Sa_time.AddEnd(omp_get_wtime()); fact_Sa_time.PrintStatMPI(0.0); Sa_timing.AddEvent(fact_Sa_time);
#else
	 TimeEvent factd_Sa_time("Salfa factorization - dense "); factd_Sa_time.AddStart(omp_get_wtime());
	SaMat = Salfa;
	SaMat.ConvertCSRToDense(1);
	eslocal info;
	char U = 'U';
	SaMat.ipiv.resize(SaMat.cols);
	dsptrf( &U, &SaMat.cols, &SaMat.dense_values[0], &SaMat.ipiv[0], &info );
	 factd_Sa_time.AddEnd(omp_get_wtime()); factd_Sa_time.PrintStatMPI(0.0); Sa_timing.AddEvent(factd_Sa_time);
#endif


	Sa_timing.totalTime.AddEnd(omp_get_wtime()); Sa_timing.PrintStatsMPI(); 
	MKL_Set_Num_Threads(1);
}


void Cluster::Create_G1_perCluster() {

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
		G1_1_time.AddStart(omp_get_wtime());
		TimeEvent G1_1_mem  ("Create G1 per clust mem: MatMat+MatTrans ");
		G1_1_mem.AddStartWOBarrier(GetProcessMemory_u());


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

				if ( esconfig::solver::REGULARIZATION == 0 ) {
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
				if ( esconfig::solver::REGULARIZATION == 0 )
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

		G1_1_time.AddEnd(omp_get_wtime());
		G1_1_time.PrintStatMPI(0.0);
		//G1_1_time.PrintLastStatMPI_PerNode(0.0);
		G1_1_mem.AddEndWOBarrier(GetProcessMemory_u());
		G1_1_mem.PrintStatMPI(0.0);
		//G1_1_mem.PrintLastStatMPI_PerNode(0.0);

		TimeEvent G1_2_time ("Create G1 per clust t. : Par.red.+MatAdd ");
		G1_2_time.AddStart(omp_get_wtime());
		TimeEvent G1_2_mem  ("Create G1 per clust mem: Par.red.+MatAdd ");
		G1_2_mem.AddStartWOBarrier(GetProcessMemory_u());

		for (eslocal j = 1; j < tmp_Mat.size(); j = j * 2 ) {
			cilk_for (eslocal i = 0; i < tmp_Mat.size(); i = i + 2*j) {
				if ( i+j < tmp_Mat.size()) {
					tmp_Mat[i    ].MatAddInPlace( tmp_Mat[i + j], 'N', 1.0 ); //  MFETI - MatAppend(tmp_Mat[i + j]);
					tmp_Mat[i + j].Clear();
				}
			}
		}

		G1_2_time.AddEnd(omp_get_wtime());
		G1_2_time.PrintStatMPI(0.0);
		//G1_2_time.PrintLastStatMPI_PerNode(0.0);

		G1_2_mem.AddEndWOBarrier(GetProcessMemory_u());
		G1_2_mem.PrintStatMPI(0.0);
		//G1_2_mem.PrintLastStatMPI_PerNode(0.0);

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

			if ( esconfig::solver::REGULARIZATION == 0 ) {
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
}

void Cluster::CreateVec_d_perCluster( SEQ_VECTOR<SEQ_VECTOR <double> > & f ) {

	eslocal size_d = domains[0].Kplus_R.cols; // because transpose of R

	if ( USE_HFETI == 1 )
		vec_d.resize( size_d );
	else
		vec_d.resize( domains.size() * size_d );	// MFETI

	if ( USE_HFETI == 1) {
		for (eslocal d = 0; d < domains.size(); d++) {
			if ( esconfig::solver::REGULARIZATION == 0 ) {
				domains[d].Kplus_R .DenseMatVec(f[d], vec_d, 'T', 0, 0         , 1.0 );
			} else {
				domains[d].Kplus_Rb.DenseMatVec(f[d], vec_d, 'T', 0, 0         , 1.0 );
			}
		}
	} else {
		for (eslocal d = 0; d < domains.size(); d++) {											// MFETI
			if ( esconfig::solver::REGULARIZATION == 0 ) {
				domains[d].Kplus_R.DenseMatVec(f[d], vec_d, 'T', 0, d * size_d, 0.0 );				// MFETI
			} else {
				domains[d].Kplus_Rb.DenseMatVec(f[d], vec_d, 'T', 0, d * size_d, 0.0 );				// MFETI
			}
		}
	}

	for (eslocal i = 0; i < vec_d.size(); i++)
		vec_d[i] = (-1.0) *  vec_d[i];

}


void Cluster::CreateVec_b_perCluster( SEQ_VECTOR<SEQ_VECTOR <double> > & f )  {

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




void Cluster::Create_Kinv_perDomain() {

	cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ )
		domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);

	if (cluster_global_index == 1)
		cout << "Creating B1*K+*B1t : ";

this->NUM_MICS = 2;

#ifdef MIC

	// compute sizes of data to be offloaded to MIC
	eslocal maxDevNumber = this->NUM_MICS;
	if (this->NUM_MICS == 0) {
		maxDevNumber = 1;
	}
	eslocal matrixPerPack = domains.size() / maxDevNumber;
	eslocal offset = 0;
	bool symmetric = true;
	this->B1KplusPacks.resize( maxDevNumber );
	eslocal * dom2dev = new eslocal[ domains.size() ];
	eslocal * offsets = new eslocal[maxDevNumber];

	for ( eslocal i = 0; i < maxDevNumber; i++ ) {
		if ( i == maxDevNumber - 1 ) {
			matrixPerPack += domains.size() % maxDevNumber;
		}

		long dataSize = 0;
		offsets[i] = offset;

		for ( eslocal j = offset; j < offset + matrixPerPack; j++ ) {
			if (!symmetric) {
				dataSize += domains[j].B1t_comp_dom.cols * domains[j].B1t_comp_dom.cols;
			} else {
				// isPacked => is symmetric
				dataSize += ( ( 1.0 + ( double ) domains[j].B1t_comp_dom.cols ) *
					( ( double ) domains[j].B1t_comp_dom.cols ) / 2.0 );
			}
			dom2dev[ j ] = i;
		}

		this->B1KplusPacks[i].Resize( matrixPerPack, dataSize );

		for ( eslocal j = offset; j < offset + matrixPerPack; j++ ) {
			this->B1KplusPacks[ i ].PreparePack( j - offset, domains[j].B1t_comp_dom.cols,
				domains[j].B1t_comp_dom.cols,  symmetric );
		}
		offset += matrixPerPack;
	}
//	tbb::mutex m;
#endif

	cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ ) {

		domains[i].KplusF.msglvl = 0;

		if ( i == 0 && cluster_global_index == 1) domains[i].KplusF.msglvl=1;

		//SolveMatF is obsolete - use Schur complement instead
		domains[i].KplusF.SolveMatF(domains[i].B1t_comp_dom, domains[i].B1Kplus, false);
		domains[i].B1Kplus.MatTranspose();

		if (cluster_global_index == 1 && i == 0)
			cout << "Creating B1*K+*B1t : ";

		if (cluster_global_index == 1) {
			//SpyText(domains[i].B1Kplus);
			//cout << "B1Klus - sparsity - " << 100.0 * (double)domains[i].B1Kplus.nnz / (double)(domains[i].B1Kplus.cols * domains[i].B1Kplus.rows) << endl;
			cout << " " << i ;
		}

		SparseMatrix Btmp;
		Btmp.MatAddInPlace(domains[i].B1Kplus, 'N', 1.0);

		domains[i].B1Kplus.Clear ();
		domains[i].B1Kplus.MatMat(Btmp,'N', domains[i].B1t_comp_dom);
		domains[i].B1Kplus.ConvertCSRToDense(0);
		//domains[i].B1Kplus.ConvertDenseToDenseFloat(0);


#ifdef CUDA

	//domains[i].B1Kplus.RemoveLowerDense();
	domains[i].B1Kplus.CopyToCUDA_Dev();
	//domains[i].B1Kplus.CopyToCUDA_Dev_fl();

#endif

#ifdef MIC
	this->B1KplusPacks[ dom2dev[ i ] ].AddDenseMatrix( i - offsets[dom2dev[i]], &(domains[i].B1Kplus.dense_values[0]) );
	domains[i].B1Kplus.Clear();
	//domains[i].B1t_comp_dom.Clear();
	//if (numDevices > 0) {
	//	domains[i].B1Kplus.CopyToMIC_Dev();
	//}
#endif

#ifdef CUDA
		if ( USE_KINV == 1 ) {
			cilk_for (eslocal d = 0; d < domains.size(); d++) {
				cudaError_t status = cudaMallocHost((void**)&domains[d].cuda_pinned_buff, domains[d].B1_comp_dom.rows * sizeof(double));
				if (status != cudaSuccess)
					printf("Error allocating pinned host memory \n");

				//status = cudaMallocHost((void**)&domains[d].cuda_pinned_buff_fl, domains[d].B1_comp_dom.rows * sizeof(float));
				//if (status != cudaSuccess)
				//	printf("Error allocating pinned host memory \n");
			}
		}
#endif


	}

	cout << endl;

#ifdef MIC
	delete [] dom2dev;
	delete [] offsets;
	if (this->NUM_MICS == 0) {
		this->B1KplusPacks[0].AllocateVectors( );
	}
	for (eslocal i = 0; i < this->NUM_MICS ; i++ ) {
		this->B1KplusPacks[i].AllocateVectors( );
		this->B1KplusPacks[i].SetDevice( i );
		this->B1KplusPacks[i].CopyToMIC();
	}

#endif

	cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ )
		domains[i].B1t_comp_dom.Clear();

//	std::cout << "This function is obsolete - use Create_SC_perDomain" << std::endl;
//	return();

}


void Cluster::Create_SC_perDomain(bool USE_FLOAT) {

	cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ )
		domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);

	if (cluster_global_index == 1)
		cout << "Creating B1*K+*B1t : using Pardiso SC : ";

	this->NUM_MICS = 2;
	#ifdef MIC

		// compute sizes of data to be offloaded to MIC
		eslocal maxDevNumber = this->NUM_MICS;
		if (this->NUM_MICS == 0) {
			maxDevNumber = 1;
		}
		eslocal matrixPerPack = domains.size() / maxDevNumber;
		eslocal offset = 0;
		bool symmetric = true;
		this->B1KplusPacks.resize( maxDevNumber );
		eslocal * dom2dev = new eslocal[ domains.size() ];
		eslocal * offsets = new eslocal[maxDevNumber];

		for ( eslocal i = 0; i < maxDevNumber; i++ ) {
			if ( i == maxDevNumber - 1 ) {
				matrixPerPack += domains.size() % maxDevNumber;
			}

			long dataSize = 0;
			offsets[i] = offset;

			for ( eslocal j = offset; j < offset + matrixPerPack; j++ ) {
				if (!symmetric) {
					dataSize += domains[j].B1t_comp_dom.cols * domains[j].B1t_comp_dom.cols;
				} else {
					// isPacked => is symmetric
					dataSize += ( ( 1.0 + ( double ) domains[j].B1t_comp_dom.cols ) *
						( ( double ) domains[j].B1t_comp_dom.cols ) / 2.0 );
				}
				dom2dev[ j ] = i;
			}

			this->B1KplusPacks[i].Resize( matrixPerPack, dataSize );

			for ( eslocal j = offset; j < offset + matrixPerPack; j++ ) {
				this->B1KplusPacks[ i ].PreparePack( j - offset, domains[j].B1t_comp_dom.cols,
					domains[j].B1t_comp_dom.cols,  symmetric );
			}
			offset += matrixPerPack;
		}
	//	tbb::mutex m;
	#endif




	cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ ) {

		if (cluster_global_index == 1) cout << "."; // << i ;

		SparseSolver tmpsps;
		if ( i == 0 && cluster_global_index == 1) tmpsps.msglvl = 1;
		tmpsps.Create_SC_w_Mat( domains[i].K, domains[i].B1t_comp_dom, domains[i].B1Kplus, false, 1 );

		if (USE_FLOAT){
			domains[i].B1Kplus.ConvertDenseToDenseFloat( 1 );
			domains[i].B1Kplus.USE_FLOAT = true;
		}

//		SparseSolver tmpsps2;
//		if ( i == 0 && cluster_global_index == 1) tmpsps2.msglvl = 1;
//		tmpsps2.Create_non_sym_SC_w_Mat( domains[i].K, domains[i].B1t_comp_dom, domains[i].B0t_comp, domains[i].B0KplusB1_comp, false, 0 );

#ifdef CUDA
		//domains[i].B1Kplus.RemoveLowerDense();
		eslocal status;
		status = domains[i].B1Kplus.CopyToCUDA_Dev();
		//domains[i].B1Kplus.CopyToCUDA_Dev_fl();
		if (status == 0)
			domains[i].isOnACC = 1;
		else
			domains[i].isOnACC = 0;
#endif

#ifdef MIC
	this->B1KplusPacks[ dom2dev[ i ] ].AddDenseMatrix( i - offsets[dom2dev[i]], &(domains[i].B1Kplus.dense_values[0]) );
	domains[i].B1Kplus.Clear();
	//domains[i].B1t_comp_dom.Clear();
	//if (numDevices > 0) {
	//	domains[i].B1Kplus.CopyToMIC_Dev();
	//}
#endif

#ifdef CUDA
		if ( USE_KINV == 1 ) {
			cilk_for (eslocal d = 0; d < domains.size(); d++) {
				cudaError_t status = cudaMallocHost((void**)&domains[d].cuda_pinned_buff, domains[d].B1_comp_dom.rows * sizeof(double));
				if (status != cudaSuccess)
					printf("Error allocating pinned host memory \n");

				//status = cudaMallocHost((void**)&domains[d].cuda_pinned_buff_fl, domains[d].B1_comp_dom.rows * sizeof(float));
				//if (status != cudaSuccess)
				//	printf("Error allocating pinned host memory \n");
			}
		}
#endif

	}


#ifdef MIC
	delete [] dom2dev;
	delete [] offsets;
	if (this->NUM_MICS == 0) {
		this->B1KplusPacks[0].AllocateVectors( );
	}
	for (eslocal i = 0; i < this->NUM_MICS ; i++ ) {
		this->B1KplusPacks[i].AllocateVectors( );
		this->B1KplusPacks[i].SetDevice( i );
		this->B1KplusPacks[i].CopyToMIC();
	}

#endif

	cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ )
		domains[i].B1t_comp_dom.Clear();

	if (cluster_global_index == 1)
		cout << endl;

}









void Cluster::compress_lambda_vector  ( SEQ_VECTOR <double> & decompressed_vec_lambda )
{
	//compress vector for CG in main loop
	for (eslocal i = 0; i < my_lamdas_indices.size(); i++)
		decompressed_vec_lambda[i] = decompressed_vec_lambda[my_lamdas_indices[i]];

	decompressed_vec_lambda.resize(my_lamdas_indices.size());
}

void Cluster::decompress_lambda_vector( SEQ_VECTOR <double> &   compressed_vec_lambda )
{
	SEQ_VECTOR <double> decompressed_vec_lambda (domains[0].B1.rows,0);

	for (eslocal i = 0; i < my_lamdas_indices.size(); i++)
		decompressed_vec_lambda[my_lamdas_indices[i]] = compressed_vec_lambda[i];

	compressed_vec_lambda = decompressed_vec_lambda;
}


void Cluster::B1_comp_MatVecSum( SEQ_VECTOR < SEQ_VECTOR <double> > & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose ) {

	std::cout << " B1_comp_MatVecSum - not implemented " << std::endl;

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
