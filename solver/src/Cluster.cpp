#include "Cluster.h"

// *******************************************************************
// **** CLUSTER CLASS ************************************************

Cluster::Cluster(int cluster_index){

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

	for (int d = 0; d < domains.size(); d++)
		domains[d].SetDynamicParameters(dynamic_timestep, dynamic_beta, dynamic_gama);

}

void Cluster::LoadCluster(string directory_path, int use_dynamic_1_no_dynamic_0, int use_kinv_1_no_kinv_0 ) {

	USE_DYNAMIC = use_dynamic_1_no_dynamic_0;
	USE_KINV    = use_kinv_1_no_kinv_0;

	// Filenames of matrices to load
	string filename_DOMAINS;

	string filename_G0;
	string filename_F0;
	string filename_Sa;

	string filename_vector_lambda;
	string filename_vector_alfa;
	string filename_vector_g0;
	string filename_vector_e0;
	// ****************

	filename_DOMAINS		= "DOMAINS";
	filename_G0				= "G0";
	filename_F0				= "F0";
	filename_Sa				= "Sapa";

	// ****************

	data_directory = string(directory_path);

	char number[8];
	sprintf(number, "%06d", cluster_global_index);
	string path = string(directory_path) + "CLUSTER_" + number + "_";

	// *** Load vector of domains *****************************************************
	LoadBinVectorInt(domains_in_global_index, string(path) + string(filename_DOMAINS));
	// *** END - Load vector of domains ***********************************************

	if (USE_HFETI == 1) {
		domains.resize( domains_in_global_index.size() );
	} else {
		domains.resize( SUBDOM_PER_CLUSTER );					// MFETI
		domains_in_global_index.resize( SUBDOM_PER_CLUSTER );	// MFETI

		//domains_in_global_index[0]     = 1 + ( SUBDOM_PER_CLUSTER * ( domains_in_global_index[0] - 1 ) );	// MFETI
		//for (int i = 1; i < SUBDOM_PER_CLUSTER; i++)														// MFETI
		//	domains_in_global_index[i] = domains_in_global_index[0] + i;									// MFETI
	}

	// *** Load all domains of the cluster ********************************************
	mkl_set_num_threads(1);

#ifdef _OPENMP
	//#pragma omp parallel for
#endif
	cilk_for (int i = 0; i < domains_in_global_index.size(); i++ ) {

		domains[i].domain_global_index = domains_in_global_index[i];

		if (USE_DYNAMIC == 1)
			domains[i].SetDynamicParameters(dynamic_timestep, dynamic_beta, dynamic_gama);

		domains[i].LoadDomain(directory_path, USE_HFETI, USE_DYNAMIC);

	}
	// *** END - Load all domains of the cluster ***************************************



	// *** Set up the dual *************************************************************
	dual_size = domains[0].B1.rows;


	// *** Alocate temporarly vectors for inter-cluster processing *********************
	// *** - based on uncompressed matrix B0
	tm1.resize(domains.size());
	tm2.resize(domains.size());
	tm3.resize(domains.size());

	for (int d = 0; d < domains.size(); d++) {
		int max_tmp_vec_size = domains[d].B0.cols;

		if (domains[d].B0.rows > domains[d].B0.cols)
			max_tmp_vec_size = domains[d].B0.rows;

		tm1[d].resize( max_tmp_vec_size );
		tm2[d].resize( max_tmp_vec_size );
		tm3[d].resize( max_tmp_vec_size );
	}
	// *** END - Alocate temporarly vectors for inter-cluster processing *****************


	// *** Create Matrices and allocate vectors for Hybrid FETI **************************
	if (USE_HFETI == 1) {

		CompressB0();

		CreateG0();
		//G0.LoadMatrix( string(path) + string(filename_G0), 'G');
		vec_g0.resize(G0.cols);
		vec_e0.resize(G0.rows);

		CreateF0();
		//F0.LoadMatrix( string(path) + string(filename_F0) );
		vec_lambda.resize(F0.m_Kplus_size);

		CreateSa();
		vec_alfa.resize(Sa.m_Kplus_size);
	}
	// *** END - Create Matrices for Hybrid FETI *****************************************




	// *** Alocate temporarly vectors for Temporary vectors for Apply_A function *********
	// *** - temporary vectors for work primal domain size *******************************
	// *** - ...
	x_prim_cluster1.resize( domains.size() );
	x_prim_cluster2.resize( domains.size() );

	for (int d = 0; d < domains.size(); d++) {
		x_prim_cluster1[d].resize( domains[d].domain_prim_size );
		x_prim_cluster2[d].resize( domains[d].domain_prim_size );
	}
	// *** END - Alocate temporarly vectors for Temporary vectors for Apply_A function ***




	// *** Detection of affinity of lag. multipliers to specific subdomains **************
	// *** - will be used to compress vectors and matrices for higher efficiency

	SEQ_VECTOR <SEQ_VECTOR <int> > lambda_map;
	LoadBinVecVec(lambda_map, string(directory_path) + string("LMtoSD"));

	//if (SUBDOM_PER_CLUSTER > 1) {															// MFETI
	//	for (int i = 0; i < lambda_map.size(); i++) {										// MFETI
	//		for (int j = 1; j < lambda_map[i][0] + 1; j++) {								// MFETI
	//			lambda_map[i][j] = 1 + ( (lambda_map[i][j] - 1) / SUBDOM_PER_CLUSTER ) ;	// MFETI
	//		}																				// MFETI
	//	}																					// MFETI
	//}																						// MFETI

	int max_dom_index = lambda_map.size(); // pozor overshot

	for (int i = 0; i < lambda_map.size(); i++)
		lambdas_filter.push_back(lambda_map[i][1]);

	SEQ_VECTOR< SEQ_VECTOR <int> > lambdas_per_subdomain (max_dom_index);

	int my_dom_index = this->cluster_global_index; // pozor
	std::SEQ_VECTOR<int>::iterator it;

	for (int lam_i = 0; lam_i < lambda_map.size(); lam_i++) {
		SEQ_VECTOR <int> my_neighs_tmp;
		int is_neigh = 0;

		if (lambda_map[lam_i][0] > 1) {

			it = unique(lambda_map[lam_i].begin() + 1, lambda_map[lam_i].begin() + 1 +  lambda_map[lam_i][0]);
			lambda_map[lam_i][0] = std::distance(lambda_map[lam_i].begin() + 1, it);
			lambda_map[lam_i].resize( std::distance(lambda_map[lam_i].begin(), it) );

			for (int dom_i = 1; dom_i < lambda_map[lam_i][0] + 1; dom_i++) {
				if ( lambda_map[lam_i][dom_i] != my_dom_index) {
					my_neighs_tmp.push_back(lambda_map[lam_i][dom_i] - 1);
				} else {
					is_neigh = 1;
					my_lamdas_indices.push_back(lam_i);

					if (dom_i == 1)
						my_lamdas_ddot_filter.push_back(1.0);
					else
						my_lamdas_ddot_filter.push_back(0.0);
				}
			}
		} else {
			if( lambda_map[lam_i][1] == my_dom_index ) {
				my_lamdas_indices.push_back(lam_i);
				my_lamdas_ddot_filter.push_back(1.0);
			}
		}

		if (is_neigh == 1) {

			for (int i = 0; i < my_neighs_tmp.size(); i++)
				lambdas_per_subdomain[my_neighs_tmp[i]].push_back(lam_i);

			my_neighs.insert(my_neighs.end(), my_neighs_tmp.begin(), my_neighs_tmp.end());
			sort(my_neighs.begin(), my_neighs.end());
			it = unique(my_neighs.begin(), my_neighs.end());
			my_neighs.resize( std::distance(my_neighs.begin(),it) );
		}
	}

	my_comm_lambdas_indices .resize(my_neighs.size());
	my_comm_lambdas			.resize(my_neighs.size());
	my_recv_lambdas			.resize(my_neighs.size());

	for (int i = 0; i < my_neighs.size(); i++) {
		my_comm_lambdas_indices[i] = lambdas_per_subdomain[my_neighs[i]];
		my_comm_lambdas[i].resize(my_comm_lambdas_indices[i].size());
		my_recv_lambdas[i].resize(my_comm_lambdas_indices[i].size());
	}

	compressed_tmp    .resize( my_lamdas_indices.size(), 0 );
	compressed_tmp2   .resize( my_lamdas_indices.size(), 0 );

	for (int d = 0; d < domains.size(); d++ )
		domains[d].compressed_tmp.resize( my_lamdas_indices.size(), 0);

	for (int i = 0; i <my_lamdas_indices.size(); i++)
		my_lamdas_map_indices.insert(make_pair(my_lamdas_indices[i],i));

	// *** END - Detection of affinity of lag. multipliers to specific subdomains ***************



	// *** Create a vector of communication pattern needed for AllReduceLambdas function *******
	my_comm_lambdas_indices_comp.resize(my_neighs.size());
	for (int i = 0; i < my_neighs.size(); i++) {
		my_comm_lambdas_indices_comp[i].resize( lambdas_per_subdomain[my_neighs[i]].size() );
		for (int j = 0; j < lambdas_per_subdomain[my_neighs[i]].size(); j++ )
			my_comm_lambdas_indices_comp[i][j] = my_lamdas_map_indices[lambdas_per_subdomain[my_neighs[i]][j]];
	}
	// *** END - Create a vector of communication pattern needed for AllReduceLambdas function *



	// *** Prepare local parts of the global matrices and vectors necessary for projector **
	if (USE_DYNAMIC == 0)
		Create_G1_perCluster  ();

	// *** Prepare the initial RBMs ********************************************************
	if (USE_DYNAMIC == 0)
		CreateVec_d_perCluster();

	//// *** Prepare the initial right hand side in dual *************************************
	//if (USE_DYNAMIC == 0)
	//	CreateVec_b_perCluster();




	// *** Compression of Matrix B1 to work with compressed lambda vectors *****************
	for (int i = 0; i < domains_in_global_index.size(); i++ ) {
		for (int j = 0; j < domains[i].B1_comp.I_row_indices.size(); j++ ) {
			int tmp_new = my_lamdas_map_indices[domains[i].B1_comp.I_row_indices [j]-1]+1;  // numbering from 1 in matrix
			domains[i].B1_comp.I_row_indices [j] = tmp_new;									// j + 1; // numbering matrix from 1
		}
		domains[i].B1_comp.rows = my_lamdas_indices.size();
		domains[i].B1_comp.ConvertToCSR( 0 );

		domains[i].B1_comp.MatTranspose(domains[i].B1t_comp);
	}
	// *** END - Compression of Matrix B1 to work with compressed lambda vectors *************


	// *** Prepare the initial right hand side in dual *************************************
	if (USE_DYNAMIC == 0)
		CreateVec_b_perCluster();


	// *** Compression of Matrix G1 to work with compressed lambda vectors *******************
	if (USE_DYNAMIC == 0) {
		G1_comp = G1;
		G1_comp.ConvertToCOO(1);
		for (int j = 0; j < G1_comp.J_col_indices.size(); j++ ) {
			int tmp_new = my_lamdas_map_indices[G1_comp.J_col_indices[j]-1]+1;  // numbering from 1 in matrix
			G1_comp.J_col_indices[j] = tmp_new; // j + 1;						// numbering matrix from 1
		}
		G1_comp.cols = my_lamdas_indices.size();
		G1_comp.ConvertToCSR(0);

		G1_comp.MatTranspose(G1t_comp);
	}
	// *** END - Compression of Matrix G1 to work with compressed lambda vectors ***************




	// *** Create inverse matrix of K+ matrix **************************************************
	if (USE_KINV == 1)
		Create_Kinv_perDomain();
	// *** END - Create inverse matrix of K+ matrix ********************************************


}


void Cluster::InitClusterPC( int * subdomains_global_indices, int number_of_subdomains ) {

	// *** Init the vector of domains *****************************************************
	//LoadBinVectorInt(domains_in_global_index, string(path) + string(filename_DOMAINS));
	domains_in_global_index.resize( number_of_subdomains ) ;
	// domains_in_global_index[0] = index_of_first_subdomain ;
	// *** END - Init the vector of domains ***********************************************

	domains.resize( number_of_subdomains );

	if (USE_HFETI == 1) {
		for (int i = 0; i < number_of_subdomains; i++)									// HFETI
			domains_in_global_index[i] = subdomains_global_indices[i] + 1;				// HFETI; [+1] -> domain numbering in espreso is from 1
	} else {
		//domains_in_global_index[0]     = 1 + ( number_of_subdomains * ( domains_in_global_index[0] - 1 ) );	// MFETI
		//for (int i = 1; i < number_of_subdomains; i++)														// MFETI
		//	domains_in_global_index[i] = domains_in_global_index[0] + i;									    // MFETI
		for (int i = 0; i < number_of_subdomains; i++)									// MFETI
			domains_in_global_index[i] = subdomains_global_indices[i] + 1;				// MFETI; [+1] -> domain numbering in espreso is from 1
	}

	// *** Init all domains of the cluster ********************************************
	for (int i = 0; i < number_of_subdomains; i++ ) {
		domains[i].domain_global_index = domains_in_global_index[i];
		domains[i].USE_KINV    = USE_KINV;
		domains[i].USE_HFETI   = USE_HFETI;
		domains[i].USE_DYNAMIC = USE_DYNAMIC;
	}
	// *** END - Init all domains of the cluster ***************************************
}

void Cluster::SetClusterPC( SEQ_VECTOR <SEQ_VECTOR <int> > & lambda_map_sub, SEQ_VECTOR < int > & neigh_domains ) {

	//USE_DYNAMIC = use_dynamic_1_no_dynamic_0;
	//USE_KINV    = use_kinv_1_no_kinv_0;

	// *** Load all domains of the cluster ********************************************
	//mkl_set_num_threads(1);

	//for (int i = 0; i < domains_in_global_index.size(); i++ ) {

		//if (USE_DYNAMIC == 1)
		//	domains[i].SetDynamicParameters(dynamic_timestep, dynamic_beta, dynamic_gama);

		//domains[i].SetDomain(USE_HFETI, USE_DYNAMIC);
				//.LoadDomain(directory_path, USE_HFETI, USE_DYNAMIC);

	//}
	// *** END - Load all domains of the cluster *****************************************



	//// *** Set up the dual size ********************************************************
	int MPIrank; 	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
	dual_size = domains[0].B1.rows;

	if (USE_HFETI == 1) {

		// *** Alocate temporarly vectors for inter-cluster processing *********************
		// *** - based on uncompressed matrix B0
		tm1.resize(domains.size());
		tm2.resize(domains.size());
		tm3.resize(domains.size());

		for (int d = 0; d < domains.size(); d++) {
			int max_tmp_vec_size = domains[d].B0.cols;

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
	my_neighs = neigh_domains;


	if (MPIrank == 0) {
		cout << " ******************************************************************************************************************************* " << endl;
		cout << " *** setting vectors for lambda map ******************************************************************************************** " << endl;
		GetProcessMemoryStat_u ( ); GetMemoryStat_u( );
		cout << endl;
	}


	my_lamdas_indices.resize( lambda_map_sub.size() );
	for (int i = 0; i < lambda_map_sub.size(); i++)
		my_lamdas_indices[i] = lambda_map_sub[i][0];


	SEQ_VECTOR< SEQ_VECTOR <int> > lambdas_per_subdomain ( domains.size() * NUMBER_OF_CLUSTERS );
	my_lamdas_ddot_filter.resize( lambda_map_sub.size(), 0.0 );
	for (int i = 0; i < lambda_map_sub.size(); i++) {
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

	for (int i = 0; i < my_neighs.size(); i++) {
		my_comm_lambdas_indices[i] = lambdas_per_subdomain[my_neighs[i]];
		my_comm_lambdas[i].resize(my_comm_lambdas_indices[i].size());
		my_recv_lambdas[i].resize(my_comm_lambdas_indices[i].size());
	}

	compressed_tmp    .resize( my_lamdas_indices.size(), 0 );
	compressed_tmp2   .resize( my_lamdas_indices.size(), 0 );

#ifdef DEVEL
	for (int d = 0; d < domains.size(); d++ )
	if (USE_KINV == 1 )
		domains[d].compressed_tmp.resize( my_lamdas_indices.size(), 0);
	else
		domains[d].compressed_tmp.resize( 1, 0);
#else
	for (int d = 0; d < domains.size(); d++ )
		domains[d].compressed_tmp.resize( my_lamdas_indices.size(), 0);
#endif

	// mapping/compression vector for cluster
	for (int i = 0; i <my_lamdas_indices.size(); i++)
		my_lamdas_map_indices.insert(make_pair(my_lamdas_indices[i],i));

	// mapping/compression vector for domains
	for (int i = 0; i < domains.size(); i++) {
		for (int j = 0; j < domains[i].lambda_map_sub.size(); j++) {
			domains[i].my_lamdas_map_indices.insert(make_pair(domains[i].lambda_map_sub[j] ,j));
		}
	}

	for (int d = 0; d < domains.size(); d++) {
		int i = 0;
		int j = 0;
		do
		{
			int big_index   = my_lamdas_indices[i];
			int small_index = domains[d].lambda_map_sub[j];

			if (big_index >  small_index) j++;

			if (big_index  < small_index) i++;

			if (big_index == small_index) {
				domains[d].lambda_map_sub_local.push_back(i);
				i++; j++;
			}


		} while ( i < my_lamdas_indices.size() && j < domains[d].lambda_map_sub.size() );
	}
	//// *** END - Detection of affinity of lag. multipliers to specific subdomains ***************



	//// *** Create a vector of communication pattern needed for AllReduceLambdas function *******
	my_comm_lambdas_indices_comp.resize(my_neighs.size());
	for (int i = 0; i < my_neighs.size(); i++) {
		my_comm_lambdas_indices_comp[i].resize( lambdas_per_subdomain[my_neighs[i]].size() );
		for (int j = 0; j < lambdas_per_subdomain[my_neighs[i]].size(); j++ )
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

#ifdef DEVEL

	cilk_for (int i = 0; i < domains_in_global_index.size(); i++ ) {

		//domains[i].B1_comp.I_row_indices = domains[i].B1.I_row_indices;
		//domains[i].B1_comp.J_col_indices = domains[i].B1.J_col_indices;
		//domains[i].B1_comp.V_values      = domains[i].B1.V_values;

		//domains[i].B1_comp.rows = domains[i].B1.rows;
		//domains[i].B1_comp.cols = domains[i].B1.cols;
		//domains[i].B1_comp.nnz  = domains[i].B1.nnz;
		//domains[i].B1_comp.type = domains[i].B1.type;

		//for (int j = 0; j < domains[i].B1_comp.I_row_indices.size(); j++ ) {
		//	int tmp_new = my_lamdas_map_indices[domains[i].B1_comp.I_row_indices [j] - 1] + 1;  // numbering from 1 in matrix
		//	domains[i].B1_comp.I_row_indices [j] = tmp_new;									    // j + 1; // numbering matrix from 1
		//}

		//domains[i].B1_comp.rows = my_lamdas_indices.size();
		//domains[i].B1_comp.ConvertToCSRwithSort( 1 );

		//domains[i].B1_comp.MatTranspose(domains[i].B1t_comp);

		//******************

		domains[i].B1_comp_dom.I_row_indices = domains[i].B1.I_row_indices;
		domains[i].B1_comp_dom.J_col_indices = domains[i].B1.J_col_indices;
		domains[i].B1_comp_dom.V_values      = domains[i].B1.V_values;

		domains[i].B1_comp_dom.rows = domains[i].B1.rows;
		domains[i].B1_comp_dom.cols = domains[i].B1.cols;
		domains[i].B1_comp_dom.nnz  = domains[i].B1.nnz;
		domains[i].B1_comp_dom.type = domains[i].B1.type;

		for (int j = 0; j < domains[i].B1_comp_dom.I_row_indices.size(); j++ ) {
			int tmp_new = domains[i].my_lamdas_map_indices[domains[i].B1_comp_dom.I_row_indices [j] - 1] + 1;  // numbering from 1 in matrix
			domains[i].B1_comp_dom.I_row_indices [j] = tmp_new;									               // j + 1; // numbering matrix from 1
		}

		domains[i].B1_comp_dom.rows = domains[i].lambda_map_sub.size();
		domains[i].B1_comp_dom.ConvertToCSRwithSort( 1 );

		domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);


		//************************

		domains[i].B1.Clear();
		domains[i].B1t.Clear();
		//domains[i].B1_comp.Clear();
		//domains[i].B1t_comp.Clear();

	}

#else

	for (int i = 0; i < domains_in_global_index.size(); i++ ) {

		domains[i].B1_comp.I_row_indices = domains[i].B1.I_row_indices;
		domains[i].B1_comp.J_col_indices = domains[i].B1.J_col_indices;
		domains[i].B1_comp.V_values      = domains[i].B1.V_values;

		domains[i].B1_comp.rows = domains[i].B1.rows;
		domains[i].B1_comp.cols = domains[i].B1.cols;
		domains[i].B1_comp.nnz  = domains[i].B1.nnz;
		domains[i].B1_comp.type = domains[i].B1.type;

		for (int j = 0; j < domains[i].B1_comp.I_row_indices.size(); j++ ) {
			int tmp_new = my_lamdas_map_indices[domains[i].B1_comp.I_row_indices [j]-1]+1;  // numbering from 1 in matrix
			domains[i].B1_comp.I_row_indices [j] = tmp_new;									// j + 1; // numbering matrix from 1
		}

		domains[i].B1_comp.rows = my_lamdas_indices.size();
		domains[i].B1_comp.ConvertToCSR( 0 );

		domains[i].B1_comp.MatTranspose(domains[i].B1t_comp);

		//******************

		domains[i].B1_comp_dom.I_row_indices = domains[i].B1.I_row_indices;
		domains[i].B1_comp_dom.J_col_indices = domains[i].B1.J_col_indices;
		domains[i].B1_comp_dom.V_values      = domains[i].B1.V_values;

		domains[i].B1_comp_dom.rows = domains[i].B1.rows;
		domains[i].B1_comp_dom.cols = domains[i].B1.cols;
		domains[i].B1_comp_dom.nnz  = domains[i].B1.nnz;
		domains[i].B1_comp_dom.type = domains[i].B1.type;

		for (int j = 0; j < domains[i].B1_comp_dom.I_row_indices.size(); j++ ) {
			int tmp_new = domains[i].my_lamdas_map_indices[domains[i].B1_comp_dom.I_row_indices [j]-1]+1;  // numbering from 1 in matrix
			domains[i].B1_comp_dom.I_row_indices [j] = tmp_new;									// j + 1; // numbering matrix from 1
		}

		domains[i].B1_comp_dom.rows = domains[i].lambda_map_sub.size();// my_lamdas_indices.size();
		domains[i].B1_comp_dom.ConvertToCSR( 0 );

		domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);


		//************************

		domains[i].B1.Clear();

	}
#endif
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
		cilk_for (int j = 0; j < G1.J_col_indices.size(); j++ )
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


	if (MPIrank == 0) {
		cout << " ******************************************************************************************************************************* " << endl;
		cout << " *** setting vectors end ******************************************************************************************************* " << endl;
		GetProcessMemoryStat_u ( ); GetMemoryStat_u( );
		cout << endl;
	}

}

void Cluster::SetClusterHFETI () {
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



	}
	// *** END - Create Matrices for Hybrid FETI *****************************************
}

void Cluster::SetClusterPC_AfterKplus () {

	//// *** Alocate temporarly vectors for Temporary vectors for Apply_A function *********
	//// *** - temporary vectors for work primal domain size *******************************
	x_prim_cluster1.resize( domains.size() );
	x_prim_cluster2.resize( domains.size() );

	for (int d = 0; d < domains.size(); d++) {
		x_prim_cluster1[d].resize( domains[d].domain_prim_size );
		x_prim_cluster2[d].resize( domains[d].domain_prim_size );
	}
	//// *** END - Alocate temporarly vectors for Temporary vectors for Apply_A function ***

	//// *** Prepare the initial right hand side in dual *************************************
	CreateVec_b_perCluster();

}


void Cluster::multKplusGlobal(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, SEQ_VECTOR<int> & cluster_map_vec) {

	vec_g0.resize(G0.cols);
	fill(vec_g0.begin(), vec_g0.end(), 0); // reset entire vector to 0

	vec_e0.resize(G0.rows);
	fill(vec_e0.begin(), vec_e0.end(), 0); // reset entire vector to 0

	//y_out.resize(x_clust_size);

	// temp vectors
	SEQ_VECTOR <double> t1,t2,t3;


	// loop over domains in the cluster
	for (int d = 0; d < domains.size(); d++) {
		int x_in_vector_start_index = x_clust_domain_map_vec[d] + cluster_map_vec[this->cluster_global_index-1];
		int domain_size = domains[d].Kplus.m_Kplus_size;

		t1.resize(domain_size);
		t2.resize(vec_g0.size());
		fill(t1.begin(), t1.end(), 0);
		fill(t2.begin(), t2.end(), 0);


		// g0
		domains[d].multKplusLocal( x_in, t1, x_in_vector_start_index, 0 );

		//					  in   out trans
		domains[d].B0.MatVec( t1 , t2,  'N' );

		for (int i = 0; i < vec_g0.size(); i++ )
			vec_g0[i] = vec_g0[i] + t2[i];

		// e0
		t2.resize(domains[d].Kplus_R.cols); // resize na pocet sloupcu matice Kplus_R - velikost jadra
		fill(t2.begin(), t2.end(), 0); // reset t2 - migh not be necessary

		domains[d].Kplus_R.MatVec(x_in, t2, 'T', x_in_vector_start_index, 0);

		int e0_start	=  d	* domains[d].Kplus_R.cols;
		int e0_end		= (d+1) * domains[d].Kplus_R.cols;
		for (int i = e0_start; i < e0_end; i++ )
			vec_e0[i] = - t2[i - e0_start];

	} // end loop over domains


	// alfa
	t1.resize(F0.m_Kplus_size);
	fill(t1.begin(), t1.end(), 0);

	F0.Solve(vec_g0, t1,0,0);

	t2.resize(G0.rows);
	fill(t2.begin(), t2.end(), 0);

	G0.MatVec(t1, t2, 'N');
	for (int i = 0; i < vec_e0.size(); i++)
		t2[i] = t2[i] - vec_e0[i];

	vec_alfa.resize(Sa.m_Kplus_size);
	Sa.Solve(t2, vec_alfa,0,0);

	// lambda
	t1.resize(G0.cols);
	fill(t1.begin(), t1.end(), 0);

	G0.MatVec(vec_alfa, t1, 'T');

	for (int i = 0; i < vec_g0.size(); i++)
		t1[i] = vec_g0[i] - t1[i];

	vec_lambda.resize(F0.m_Kplus_size);
	F0.Solve(t1, vec_lambda,0,0);

	// Kplus_x
	for (int d = 0; d < domains.size(); d++) {
		//int x_in_vector_start_index = x_clust_domain_map_vec[d];
		int x_in_vector_start_index = x_clust_domain_map_vec[d] + cluster_map_vec[this->cluster_global_index-1];
		int domain_size = domains[d].Kplus.m_Kplus_size;

		t1.resize(domain_size);
		fill(t1.begin(), t1.end(), 0);

		domains[d].B0.MatVec(vec_lambda, t1, 'T');

		for (int i = 0; i < domain_size; i++)
			t1[i] = x_in[x_in_vector_start_index + i] - t1[i];

		t2.resize(domain_size);
		fill(t2.begin(), t2.end(), 0);

		domains[d].multKplusLocal(t1 , t2, 0, 0);

		t3.resize(domain_size);
		fill(t3.begin(), t3.end(), 0);

		int e0_start	=  d	* domains[d].Kplus_R.cols;
		int e0_end		= (d+1) * domains[d].Kplus_R.cols;
		domains[d].Kplus_R.MatVec(vec_alfa, t3, 'N', e0_start,0);

		for (int i = 0; i < domain_size; i++)
			y_out[x_in_vector_start_index + i] = t2[i] + t3[i];

		double t_sum = 0;
		for (int i = 0; i < domain_size; i++)
			t_sum +=  (t2[i] + t3[i]);

		t_sum = t_sum;

	}

	t1.clear();
	t2.clear();
	t3.clear();
}

void Cluster::multKplusGlobal_l(SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {

	int num_threads = domains.size();

	mkl_set_num_threads(1);

	cluster_time.totalTime.AddStart();

	vec_fill_time.AddStart();
	fill(vec_g0.begin(), vec_g0.end(), 0); // reset entire vector to 0
	//fill(vec_e0.begin(), vec_e0.end(), 0); // reset entire vector to 0
	vec_fill_time.AddEnd();

	// loop over domains in the cluster
	loop_1_1_time.AddStart();
	cilk_for (int d = 0; d < domains.size(); d++)
	{
		domains[d].B0Kplus_comp.DenseMatVec(x_in[d], tm2[d]);			// g0 - with comp B0Kplus
		domains[d].Kplus_R.DenseMatVec(x_in[d], tm3[d], 'T');			// e0
	}
	loop_1_1_time.AddEnd();

	loop_1_2_time.AddStart();
	cilk_for (int d = 0; d < domains.size(); d++)
	{
		int e0_start	=  d	* domains[d].Kplus_R.cols;
		int e0_end		= (d+1) * domains[d].Kplus_R.cols;

		for (int i = e0_start; i < e0_end; i++ )
			vec_e0[i] = - tm3[d][i - e0_start];
	}


	for (int d = 0; d < domains.size(); d++)
		for (int i = 0; i < domains[d].B0Kplus_comp.rows; i++)
			vec_g0[ domains[d].B0_comp_map_vec[i] - 1 ] += tm2[d][i];

	// end loop over domains
	loop_1_2_time.AddEnd();

	mkl_set_num_threads(24); //pozor
	clusCP_time.AddStart();

	clus_F0_1_time.AddStart();
	F0.Solve(vec_g0, tm1[0], 0, 0);
	clus_F0_1_time.AddEnd();

	clus_G0_time.AddStart();
	G0.MatVec(tm1[0], tm2[0], 'N');
	clus_G0_time.AddEnd();

	cilk_for (int i = 0; i < vec_e0.size(); i++)
		tm2[0][i] = tm2[0][i] - vec_e0[i];
	//cblas_daxpy(vec_e0.size(), -1.0, &vec_e0[0], 1, &tm2[0][0], 1);

	clus_Sa_time.AddStart();
	Sa.Solve(tm2[0], vec_alfa,0,0);
	clus_Sa_time.AddEnd();

	clus_G0t_time.AddStart();
	G0.MatVec(vec_alfa, tm1[0], 'T'); 	// lambda
	clus_G0t_time.AddEnd();

	cilk_for (int i = 0; i < vec_g0.size(); i++)
		tm1[0][i] = vec_g0[i] - tm1[0][i];


	clus_F0_2_time.AddStart();
	F0.Solve(tm1[0], vec_lambda,0,0);
	clus_F0_2_time.AddEnd();

	clusCP_time.AddEnd();


	// Kplus_x
	mkl_set_num_threads(1);
	loop_2_1_time.AddStart();
	cilk_for (int d = 0; d < domains.size(); d++)
	{
		int domain_size = domains[d].domain_prim_size;


		SEQ_VECTOR < double > tmp_vec (domains[d].B0_comp_map_vec.size(), 0);
		for (int i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
			tmp_vec[i] = vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;
		domains[d].B0t_comp.MatVec(tmp_vec, tm1[d], 'N');

		for (int i = 0; i < domain_size; i++)
			tm1[d][i] = x_in[d][i] - tm1[d][i];

		domains[d].multKplusLocal(tm1[d] , tm2[d], 0, 0);

		int e0_start	=  d	* domains[d].Kplus_R.cols;
		int e0_end		= (d+1) * domains[d].Kplus_R.cols;

		domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);

		for (int i = 0; i < domain_size; i++)
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
	cilk_for (int d = 0; d < domains.size(); d++)
	{
		domains[d].B0Kplus_comp.DenseMatVec(x_in[d], tm2[d]);			// g0 - with comp B0Kplus
		domains[d].Kplus_R.     DenseMatVec(x_in[d], tm3[d], 'T');		// e0
	}
	loop_1_1_time.AddEnd();

	loop_1_2_time.AddStart();
	cilk_for (int d = 0; d < domains.size(); d++)
	{
		int e0_start	=  d	* domains[d].Kplus_R.cols;
		int e0_end		= (d+1) * domains[d].Kplus_R.cols;

		for (int i = e0_start; i < e0_end; i++ )
			vec_e0[i] = - tm3[d][i - e0_start];
	}

	for (int d = 0; d < domains.size(); d++)
		for (int i = 0; i < domains[d].B0Kplus_comp.rows; i++)
			vec_g0[ domains[d].B0_comp_map_vec[i] - 1 ] += tm2[d][i];

	// end loop over domains
	loop_1_2_time.AddEnd();

	mkl_set_num_threads(24); //pozor
	clusCP_time.AddStart();

	clus_F0_1_time.AddStart();
	F0.Solve(vec_g0, tm1[0], 0, 0);
	clus_F0_1_time.AddEnd();

	clus_G0_time.AddStart();
	G0.MatVec(tm1[0], tm2[0], 'N');
	clus_G0_time.AddEnd();

	cilk_for (int i = 0; i < vec_e0.size(); i++)
		tm2[0][i] = tm2[0][i] - vec_e0[i];
	//cblas_daxpy(vec_e0.size(), -1.0, &vec_e0[0], 1, &tm2[0][0], 1);

	clus_Sa_time.AddStart();
	Sa.Solve(tm2[0], vec_alfa,0,0);
	clus_Sa_time.AddEnd();

	clus_G0t_time.AddStart();
	G0.MatVec(vec_alfa, tm1[0], 'T'); 	// lambda
	clus_G0t_time.AddEnd();

	cilk_for (int i = 0; i < vec_g0.size(); i++)
		tm1[0][i] = vec_g0[i] - tm1[0][i];


	clus_F0_2_time.AddStart();
	F0.Solve(tm1[0], vec_lambda,0,0);
	clus_F0_2_time.AddEnd();

	clusCP_time.AddEnd();


	// Kplus_x
	mkl_set_num_threads(1);
	loop_2_1_time.AddStart();
	cilk_for (int d = 0; d < domains.size(); d++)
	{
		int domain_size = domains[d].domain_prim_size;

		SEQ_VECTOR < double > tmp_vec (domains[d].B0_comp_map_vec.size(), 0.0);
		for (int i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
			tmp_vec[i] = -1.0 * vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;
		domains[d].B0Kplus_comp.DenseMatVec(tmp_vec, tm2[d], 'T' );

		int e0_start	=  d	* domains[d].Kplus_R.cols;
		int e0_end		= (d+1) * domains[d].Kplus_R.cols;
		domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);

		for (int i = 0; i < domain_size; i++)
			x_in[d][i] = tm2[d][i] + tm3[d][i];

	}
	loop_2_1_time.AddEnd();

	cluster_time.totalTime.AddEnd();
}


////backup March 31 2015
//void Cluster::multKplusGlobal_l(SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {
//
//	//int MPIrank;
//	//MPI_Comm_rank (MPI_COMM_WORLD, &MPIrank);
//	//if (MPIrank == 0 ) { cout << "MultKplusGlobal - Cilk workers = " << __cilkrts_get_nworkers()      << endl; }
//	//if (MPIrank == 0 ) { cout << "MultKplusGlobal - Cilk workers = " << __cilkrts_get_total_workers() << endl; }
//	//
//
//	int num_threads = domains.size();
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
//	cilk_for (int d = 0; d < domains.size(); d++)
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
//	cilk_for (int d = 0; d < domains.size(); d++)
//	{
//		int e0_start	=  d	* domains[d].Kplus_R.cols;
//		int e0_end		= (d+1) * domains[d].Kplus_R.cols;
//
//		for (int i = e0_start; i < e0_end; i++ )
//			vec_e0[i] = - tm3[d][i - e0_start];
//	}
//
//
//	for (int d = 0; d < domains.size(); d++)
//		for (int i = 0; i < domains[d].B0Kplus_comp.rows; i++)
//			vec_g0[ domains[d].B0_comp_map_vec[i] - 1 ] += tm2[d][i];
//
//
//	//PROD -
//	//cilk_for (int i = 0; i < vec_g0.size(); i++ )
//	//{
//	//	vec_g0[i] = tm2[0][i];
//	//}
//	//for (int d = 1; d < domains.size(); d++) {
//	//	cilk_for (int i = 0; i < vec_g0.size(); i++ )
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
//	cilk_for (int i = 0; i < vec_e0.size(); i++)
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
//	cilk_for (int i = 0; i < vec_g0.size(); i++)
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
//	cilk_for (int d = 0; d < domains.size(); d++)
//	{
//		int domain_size = domains[d].domain_prim_size;
//
//
//		SEQ_VECTOR < double > tmp_vec (domains[d].B0_comp_map_vec.size(), 0);
//		for (int i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
//			tmp_vec[i] = vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;
//		domains[d].B0t_comp.MatVec(tmp_vec, tm1[d], 'N');
//
//
//		//domains[d].B0.MatVec(vec_lambda, tm1[d], 'T');  // both about the same performance
//		//domains[d].B0t.MatVec(vec_lambda, tm1[d], 'N');   // both about the same performance
//
//
//		for (int i = 0; i < domain_size; i++)
//			tm1[d][i] = x_in[d][i] - tm1[d][i];
//
//		domains[d].multKplusLocal(tm1[d] , tm2[d], 0, 0);
//
//		int e0_start	=  d	* domains[d].Kplus_R.cols;
//		int e0_end		= (d+1) * domains[d].Kplus_R.cols;
//
//		//domains[d].Kplus_R.MatVec(vec_alfa, tm3[d], 'N', e0_start,0);
//		domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);
//
//		for (int i = 0; i < domain_size; i++)
//			x_in[d][i] = tm2[d][i] + tm3[d][i];
//
//	}
//	loop_2_1_time.AddEnd();
//
//	cluster_time.totalTime.AddEnd();
//}


void Cluster::CompressB0() {

	cilk_for (int d = 0; d < domains.size(); d++) {

		domains[d].B0.MatTranspose(domains[d].B0t);

		domains[d].B0_comp = domains[d].B0;
		domains[d].B0_comp.ConvertToCOO(1);
		domains[d].B0_comp_map_vec = domains[d].B0_comp.I_row_indices;

		for (int i = 0; i < domains[d].B0_comp.I_row_indices.size(); i++)
			domains[d].B0_comp.I_row_indices[i] = i + 1;

		domains[d].B0_comp.rows = domains[d].B0_comp.I_row_indices.size();
		domains[d].B0_comp.ConvertToCSR(1);
		domains[d].B0_comp.MatTranspose(domains[d].B0t_comp);

		domains[d].Kplus_R.ConvertCSRToDense(0); // pozor - keep CSR data

	}

}

void Cluster::CreateG0() {

	mkl_set_num_threads(1);

	SEQ_VECTOR <SparseMatrix> G0LocalTemp( domains.size() );

	cilk_for (int i = 0; i<domains.size(); i++) {
		G0LocalTemp[i].MatMat(domains[i].B0, 'N', domains[i].Kplus_R );
		G0LocalTemp[i].MatTranspose(-1.0);
	}

	for (int i = 0; i<domains.size(); i++) {
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
	SparseMatrix F0_Mat;

	if (MPIrank == 0 ) {cout << "HFETI - Create F0 - domain : " << endl; };

	TimeEvent solve_F0_time("B0 compression; F0 multiple RHS solve");
	solve_F0_time.AddStart(omp_get_wtime());


	cilk_for (int d = 0; d < domains.size(); d++) {

		domains[d].Kplus.SolveMat_Dense(domains[d].B0t_comp, domains[d].B0Kplus_comp);

		domains[d].B0Kplus = domains[d].B0Kplus_comp;
		domains[d].B0Kplus_comp.MatTranspose();
		domains[d].B0Kplus_comp.ConvertCSRToDense(1);

		for (int i = 0; i < domains[d].B0Kplus.CSR_J_col_indices.size() - 1; i++)
			domains[d].B0Kplus.CSR_J_col_indices[i] = domains[d].B0_comp_map_vec [ domains[d].B0Kplus.CSR_J_col_indices[i] - 1 ];

		domains[d].B0Kplus.cols = domains[d].B0.rows;;

		tmpF0v[d].MatMat(domains[d].B0, 'N', domains[d].B0Kplus);
		domains[d].B0Kplus.Clear();

		if (MPIrank == 0 ) {cout << d << " "; };

	}

	if (MPIrank == 0 ) {cout << endl; };

	solve_F0_time.AddEnd(omp_get_wtime());
	solve_F0_time.PrintStatMPI(0.0);
	F0_timing.AddEvent(solve_F0_time);

	if (MPIrank == 0 ) {cout << endl; };

	TimeEvent reduction_F0_time("F0 reduction time");
	reduction_F0_time.AddStart(omp_get_wtime());

	for (int j = 1; j < tmpF0v.size(); j = j * 2 ) {
		cilk_for (int i = 0; i < tmpF0v.size(); i = i + 2*j) {
			if ( i+j < tmpF0v.size()) {
				tmpF0v[i    ].MatAddInPlace( tmpF0v[i + j], 'N', 1.0 );
				tmpF0v[i + j].Clear();
			}
		}
	}
	F0_Mat = tmpF0v[0];

	reduction_F0_time.AddEnd(omp_get_wtime()); reduction_F0_time.PrintStatMPI(0.0); F0_timing.AddEvent(reduction_F0_time);


	TimeEvent fact_F0_time("B0 Kplus Factorization "); fact_F0_time.AddStart(omp_get_wtime());

	mkl_set_num_threads(24);
	F0_Mat.RemoveLower();
	F0.ImportMatrix(F0_Mat);

	F0_fast.ImportMatrix(F0_Mat);

	F0_Mat.Clear();
	F0.Factorization();
	mkl_set_num_threads(1);

	if (MPIrank == 0)
		F0.msglvl = 0;

	fact_F0_time.AddEnd(omp_get_wtime()); fact_F0_time.PrintStatMPI(0.0); F0_timing.AddEvent(fact_F0_time);


	F0_timing.totalTime.AddEnd(omp_get_wtime());
	F0_timing.PrintStatsMPI();

	// *** POZOR **************************************************************
	for (int d = 0; d<domains.size(); d++) {
		domains[d].B0.Clear();
		domains[d].B0t.Clear();
	}
};

void Cluster::CreateSa() {

	MKL_Set_Num_Threads(24);
	int MPIrank; MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

	SparseMatrix Salfa;
	SparseMatrix tmpM;

	TimeEval Sa_timing (" HFETI - Salfa preprocessing timing"); Sa_timing.totalTime.AddStart(omp_get_wtime());

	 TimeEvent G0trans_Sa_time("G0 transpose"); G0trans_Sa_time.AddStart(omp_get_wtime());
	SparseMatrix G0t;
	G0.MatTranspose(G0t);
	 G0trans_Sa_time.AddEnd(omp_get_wtime()); G0trans_Sa_time.PrintStatMPI(0.0); Sa_timing.AddEvent(G0trans_Sa_time);

	 TimeEvent G0solve_Sa_time("SolveMatF with G0t as RHS"); G0solve_Sa_time.AddStart(omp_get_wtime());
	if (MPIrank == 0) F0_fast.msglvl = 1;
	F0_fast.SolveMatF(G0t,tmpM, true);
	if (MPIrank == 0) F0_fast.msglvl = 0;
	 G0solve_Sa_time.AddEnd(omp_get_wtime()); G0solve_Sa_time.PrintStatMPI(0.0); Sa_timing.AddEvent(G0solve_Sa_time);


	 TimeEvent SaMatMat_Sa_time("Salfa = MatMat G0 * solution "); SaMatMat_Sa_time.AddStart(omp_get_wtime());
	Salfa.MatMat(G0, 'N', tmpM);
	Salfa.RemoveLower();
	tmpM.Clear();
	 SaMatMat_Sa_time.AddEnd(omp_get_wtime()); SaMatMat_Sa_time.PrintStatMPI(0.0); Sa_timing.AddEvent(SaMatMat_Sa_time);


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

	 TimeEvent fact_Sa_time("Salfa factorization "); fact_Sa_time.AddStart(omp_get_wtime());
	if (MPIrank == 0) Sa.msglvl = 0;
	//Sa.LoadMatrix( string(path) + string(filename_Sa) );
	Sa.ImportMatrix(Salfa);
	Sa.Factorization();
	if (MPIrank == 0) Sa.msglvl = 0;

	 fact_Sa_time.AddEnd(omp_get_wtime()); fact_Sa_time.PrintStatMPI(0.0); Sa_timing.AddEvent(fact_Sa_time);

	Sa_timing.totalTime.AddEnd(omp_get_wtime()); Sa_timing.PrintStatsMPI();
	MKL_Set_Num_Threads(1);
}


void Cluster::Create_G1_perCluster() {

	SparseMatrix tmpM;

	if (USE_HFETI == 1) {
		//G1 = G1 + trans(B1 * domains[d].Kplus_R) for all domains
		//for (int d = 0; d < domains.size(); d++)
		//{
		//	tmpM.MatMat( domains[d].B1, 'N', domains[d].Kplus_R);
		//	G1.MatAddInPlace( tmpM, 'N', 1.0 );
		//	tmpM.Clear();
		//}
		//G1.MatTranspose();

		//// OK - but sequential
		//for (int d = 0; d < domains.size(); d++)					//HFETI
		//{
		//	tmpM.MatMat( domains[d].B1t, 'T', domains[d].Kplus_R);
		//	G1.MatAddInPlace( tmpM, 'N', 1.0 );
		//	tmpM.Clear();
		//}
		//G1.MatTranspose();




		//int threads = 24;
		//int n_domains = domains.size();
		//vector < SparseMatrix > tmp_Mat (threads);
		//cilk_for (int t = 0; t < threads; t++ ) {
		//	for (int i = t*(n_domains/threads+1); i < (t+1)*(n_domains/threads+1); i++ ) {
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
		//vector <int > G_I_row_indices;
		//G_I_row_indices.resize(B.nnz);

		//int indx = 0;
		//for (int r = 0; r < Rt.rows; r++)
		//	tmpG[indx][r] += B.V_values[0] * Rt.dense_values[Rt.rows * (B.J_col_indices[0]-1) + r];
		//
		//G_I_row_indices[indx] = B.I_row_indices[0];

		//for (int i = 1; i < B.I_row_indices.size(); i++) {
		//
		//	if (B.I_row_indices[i-1] != B.I_row_indices[i])
		//		indx++;

		//	for (int r = 0; r < Rt.rows; r++)
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

		//for (int i = 0; i < tmpG.size(); i++) {
		//	for (int j = 0; j < tmpG[i].size(); j++){
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
		cilk_for (int j = 0; j < tmp_Mat.size(); j++) {
			// V1
			//tmp_Mat[j].MatMat( domains[j].B1t, 'T', domains[j].Kplus_R);
			//tmp_Mat[j].MatTranspose();

			// V2 - not cool
			//tmp_Mat[j].MatMatSorted( domains[j].Kplus_R, 'T', domains[j].B1t); // - pozor mooooc pomale a MKL rutina zere mooooooc pameti

			// V3

			SparseMatrix Rt;
			SparseMatrix B;

			domains[j].Kplus_R.MatTranspose(Rt);
			//Rt = domains[j].Kplus_R;
			//Rt.MatTranspose();

			Rt.ConvertCSRToDense(1);
			B = domains[j].B1;

			SEQ_VECTOR < SEQ_VECTOR < double > > tmpG (B.nnz, SEQ_VECTOR <double> (Rt.rows,0));
			SEQ_VECTOR <int > G_I_row_indices;
			G_I_row_indices.resize(B.nnz);

			int indx = 0;
			for (int r = 0; r < Rt.rows; r++)
				tmpG[indx][r] += B.V_values[0] * Rt.dense_values[Rt.rows * (B.J_col_indices[0]-1) + r];

			G_I_row_indices[indx] = B.I_row_indices[0];

			for (int i = 1; i < B.I_row_indices.size(); i++) {

				if (B.I_row_indices[i-1] != B.I_row_indices[i])
					indx++;

				for (int r = 0; r < Rt.rows; r++)
					tmpG[indx][r] += B.V_values[i] * Rt.dense_values[Rt.rows * (B.J_col_indices[i]-1) + r];

				G_I_row_indices[indx] = B.I_row_indices[i];
			}

			G_I_row_indices.resize(indx+1);
			tmpG.resize(indx+1);

			SEQ_VECTOR <int>    G_I; G_I.reserve( tmpG.size() * tmpG[0].size());
			SEQ_VECTOR <int>    G_J; G_J.reserve( tmpG.size() * tmpG[0].size());
			SEQ_VECTOR <double> G_V; G_V.reserve( tmpG.size() * tmpG[0].size());

			for (int i = 0; i < tmpG.size(); i++) {
				for (int j = 0; j < tmpG[i].size(); j++){
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

		for (int j = 1; j < tmp_Mat.size(); j = j * 2 ) {
			cilk_for (int i = 0; i < tmp_Mat.size(); i = i + 2*j) {
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
		cilk_for (int j = 0; j < domains.size(); j++) {

			// V1
			//tmp_Mat[j].MatMat( domains[j].B1t, 'T', domains[j].Kplus_R);
			//tmp_Mat[j].MatTranspose();

			// V2
			//tmp_Mat[j].MatMat( domains[j].Kplus_R, 'T', domains[j].B1t);

			// V3
			SparseMatrix Rt;
			SparseMatrix B;

			domains[j].Kplus_R.MatTranspose(Rt);
			Rt.ConvertCSRToDense(1);
			B = domains[j].B1;

			SEQ_VECTOR < SEQ_VECTOR < double > > tmpG (B.nnz, SEQ_VECTOR <double> (Rt.rows,0));
			SEQ_VECTOR <int > G_I_row_indices;
			G_I_row_indices.resize(B.nnz);

			int indx = 0;
			for (int r = 0; r < Rt.rows; r++)
				tmpG[indx][r] += B.V_values[0] * Rt.dense_values[Rt.rows * (B.J_col_indices[0]-1) + r];

			G_I_row_indices[indx] = B.I_row_indices[0];

			for (int i = 1; i < B.I_row_indices.size(); i++) {

				if (B.I_row_indices[i-1] != B.I_row_indices[i])
					indx++;

				for (int r = 0; r < Rt.rows; r++)
					tmpG[indx][r] += B.V_values[i] * Rt.dense_values[Rt.rows * (B.J_col_indices[i]-1) + r];

				G_I_row_indices[indx] = B.I_row_indices[i];
			}

			G_I_row_indices.resize(indx+1);
			tmpG.resize(indx+1);

			SEQ_VECTOR <int>    G_I; G_I.reserve( tmpG.size() * tmpG[0].size());
			SEQ_VECTOR <int>    G_J; G_J.reserve( tmpG.size() * tmpG[0].size());
			SEQ_VECTOR <double> G_V; G_V.reserve( tmpG.size() * tmpG[0].size());

			for (int i = 0; i < tmpG.size(); i++) {
				for (int j = 0; j < tmpG[i].size(); j++){
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

		for (int j = 1; j < tmp_Mat.size(); j = j * 2 ) {
			cilk_for (int i = 0; i < tmp_Mat.size(); i = i + 2*j) {
				if ( i+j < tmp_Mat.size()) {
					tmp_Mat[i    ].MatAppend(tmp_Mat[i + j]);
					tmp_Mat[i + j].Clear();
				}
			}
		}
		G1.MatAppend(tmp_Mat[0]);
		//for (int d = 0; d < domains.size(); d++)					//MFETI
		//{
		//	//tmpM.MatMat( domains[d].B1, 'N', domains[d].Kplus_R);
		//	tmpM.MatMat( domains[d].B1t, 'T', domains[d].Kplus_R);
		//	tmpM.MatTranspose();
		//	G1.MatAppend(tmpM);
		//	tmpM.Clear();
		//}
	}

	// for both MFETI and HFETI
	for (int i = 0; i < G1.CSR_V_values.size(); i++)
		G1.CSR_V_values[i] = -1.0 * G1.CSR_V_values[i];
}

void Cluster::CreateVec_d_perCluster() {
	int size_d = domains[0].Kplus_R.cols; // because transpose of R


	if ( USE_HFETI == 1 )
		vec_d.resize( size_d );
	else
		vec_d.resize( domains.size() * size_d );	// MFETI

	if ( USE_HFETI == 1) {
		for (int d = 0; d < domains.size(); d++)
			domains[d].Kplus_R.MatVec(domains[d].f, vec_d, 'T', 0, 0         , 1.0 );
	} else {
		for (int d = 0; d < domains.size(); d++)										// MFETI
			domains[d].Kplus_R.MatVec(domains[d].f, vec_d, 'T', 0, d * size_d, 0.0 );	// MFETI
	}

	for (int i = 0; i < vec_d.size(); i++)
		vec_d[i] = (-1.0) *  vec_d[i];
}

void Cluster::CreateVec_b_perCluster() {

	//int MPIrank;
	//MPI_Comm_rank (MPI_COMM_WORLD, &MPIrank);

	SEQ_VECTOR<SEQ_VECTOR<double> > x_prim_cluster;
	for (int d = 0; d < domains.size(); d++) {
		x_prim_cluster.push_back( domains[d].f );
	}

	if (USE_HFETI == 0) {
		cilk_for (int d = 0; d < domains.size(); d++) {				// MFETI
			domains[d].multKplusLocal( x_prim_cluster[d] );
		}
	} else {
		multKplusGlobal_l( x_prim_cluster );						// HFETI
	}

	// pro ukoncovani v primaru - vypocet up0
	cilk_for (int d = 0; d < domains.size(); d++) {
		domains[d].up0 = x_prim_cluster[d];
	}

#ifdef DEVEL
	vec_b_compressed.resize(my_lamdas_indices.size(), 0.0);

	SEQ_VECTOR < double > y_out_tmp (domains[0].B1_comp_dom.rows);

	for (int d = 0; d < domains.size(); d++) {
		y_out_tmp.resize( domains[d].B1_comp_dom.rows );
		domains[d].B1_comp_dom.MatVec (x_prim_cluster[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out
		for (int i = 0; i < domains[d].lambda_map_sub_local.size(); i++)
			vec_b_compressed[ domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i];
	}


#else
	vec_b_compressed.resize(my_lamdas_indices.size());
	for (int d = 0; d < domains.size(); d++)
		domains[d].B1_comp.MatVec( x_prim_cluster[d], vec_b_compressed, 'N', 0, 0, 1.0);
#endif

	int a = 0;
}


void Cluster::Create_Kinv_perDomain() {

	if (cluster_global_index == 1)
		cout << "Creating B1*K+*B1t : ";

#ifdef DEVEL
	cilk_for (int i = 0; i < domains_in_global_index.size(); i++ ) {

		domains[i].KplusF.msglvl = 0;

		if ( i == 0 && cluster_global_index == 1) domains[i].KplusF.msglvl=1;

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


#ifdef CUDA
		if ( USE_KINV == 1 ) {
			cilk_for (int d = 0; d < domains.size(); d++) {
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

	if (cluster_global_index == 1)
		cout << endl;

#else
	cilk_for (int i = 0; i < domains_in_global_index.size(); i++ ) {

		domains[i].KplusF.msglvl = 0;

		//if (cluster_global_index == 1)
		if ( i == 0 && cluster_global_index == 1)
			domains[i].KplusF.msglvl=1;

		domains[i].KplusF.SolveMatF(domains[i].B1t_comp, domains[i].B1Kplus);
		domains[i].B1Kplus.MatTranspose();

		if (cluster_global_index == 1) {
			//SpyText(domains[i].B1Kplus);
			cout << "B1Klus - sparsity - " << 100.0 * (double)domains[i].B1Kplus.nnz / (double)(domains[i].B1Kplus.cols * domains[i].B1Kplus.rows) << endl;
		}

		SparseMatrix Btmp;
		Btmp.MatAddInPlace(domains[i].B1Kplus, 'N', 1.0);

		domains[i].B1Kplus.Clear ();
		domains[i].B1Kplus.MatMat(Btmp,'N', domains[i].B1t_comp);
		domains[i].B1Kplus.ConvertCSRToDense(0);
	}
#endif // DEVEL

}



void Cluster::Create_SC_perDomain() {

	if (cluster_global_index == 1)
		cout << "Creating B1*K+*B1t : using Pardiso SC";

	cilk_for (int i = 0; i < domains_in_global_index.size(); i++ ) {

		if (cluster_global_index == 1)
			cout << " " << i ;

		domains[i].KplusF.msglvl = 0;

		if ( i == 0 && cluster_global_index == 1) domains[i].KplusF.msglvl=1;


		SparseMatrix K_sc1;
		SparseMatrix K_b_tmp;
		SparseMatrix Sc_eye;

		K_b_tmp = domains[i].B1t_comp_dom;
		K_b_tmp.MatTranspose();

		Sc_eye.CreateEye(K_b_tmp.rows, 0.0, 0, K_b_tmp.cols);

		K_sc1 = domains[i].K;
		K_sc1.MatTranspose();
		K_sc1.MatAppend(K_b_tmp);
		K_sc1.MatTranspose();
		K_sc1.MatAppend(Sc_eye);

		domains[i].KplusF.ImportMatrix(K_sc1);
		domains[i].KplusF.Create_SC(domains[i].B1Kplus, K_b_tmp.rows, false);
		domains[i].B1Kplus.type = 'G';

		SparseMatrix SC_tmp;
		SC_tmp = domains[i].B1Kplus;
		SC_tmp.SetDiagonalOfSymmetricMatrix(0.0);
		SC_tmp.MatTranspose();

		domains[i].B1Kplus.MatAddInPlace(SC_tmp,'N',1.0);

		domains[i].B1Kplus.MatScale(-1.0);

		domains[i].B1Kplus.ConvertCSRToDense(0);
		//domains[i].B1Kplus.ConvertDenseToDenseFloat(0);




#ifdef CUDA
		//domains[i].B1Kplus.RemoveLowerDense();
		domains[i].B1Kplus.CopyToCUDA_Dev();
		//domains[i].B1Kplus.CopyToCUDA_Dev_fl();
#endif


#ifdef CUDA
		if ( USE_KINV == 1 ) {
			cilk_for (int d = 0; d < domains.size(); d++) {
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

	if (cluster_global_index == 1)
	cout << endl;

}









void Cluster::compress_lambda_vector  ( SEQ_VECTOR <double> & decompressed_vec_lambda )
{
	//compress vector for CG in main loop
	for (int i = 0; i < my_lamdas_indices.size(); i++)
		decompressed_vec_lambda[i] = decompressed_vec_lambda[my_lamdas_indices[i]];

	decompressed_vec_lambda.resize(my_lamdas_indices.size());
}

void Cluster::decompress_lambda_vector( SEQ_VECTOR <double> &   compressed_vec_lambda )
{
	SEQ_VECTOR <double> decompressed_vec_lambda (domains[0].B1.rows,0);

	for (int i = 0; i < my_lamdas_indices.size(); i++)
		decompressed_vec_lambda[my_lamdas_indices[i]] = compressed_vec_lambda[i];

	compressed_vec_lambda = decompressed_vec_lambda;
}


void Cluster::B1_comp_MatVecSum( SEQ_VECTOR < SEQ_VECTOR <double> > & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose ) {

	//if (domains.size() <= 3) {

		domains[0].B1_comp.MatVec (x_in[0], y_out, T_for_transpose_N_for_non_transpose, 0, 0, 0.0); // first vector overwrites temp vector
		for (int d = 1; d < domains.size(); d++) // reduction
			domains[d].B1_comp.MatVec (x_in[d], y_out, T_for_transpose_N_for_non_transpose, 0, 0, 1.0); // will add (summation per elements) all partial results into y_out

	//} else {
	//
	//	cilk_for (int d = 0; d < domains.size(); d++) {
	//		domains[d].B1_comp.MatVec (x_in[d], domains[d].compressed_tmp, T_for_transpose_N_for_non_transpose, 0, 0, 0.0);
	//	}

	//	for ( int r = 2 * (domains.size() / 2);  r > 1;  r = r / 2 ) {
	//		cilk_for ( int d = 0;  d < r;  d++ ) {
	//			if ( (r + d)  <  domains.size() ) {
	//				for (int i = 0; i < domains[d].compressed_tmp.size(); i++) {
	//					domains[d].compressed_tmp[i] = domains[d].compressed_tmp[i] + domains[r + d].compressed_tmp[i];
	//				}
	//			}
	//		}
	//	}

	//	for (int i = 0; i < domains[0].compressed_tmp.size(); i++)
	//		y_out[i] = domains[0].compressed_tmp[i] + domains[1].compressed_tmp[i];
	//
	//}

}

// **** END - CLUSTER CLASS ************************************************
// *******************************************************************
