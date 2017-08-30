
#include "../specific/cluster.h"

#include "../../assembler/instance.h"
#include "../../basis/utilities/utils.h"

#ifdef READEX_LEVEL_1
#include "readex.h"
#include "readex_regions.h"
#endif

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


void ClusterBase::InitClusterPC( eslocal * subdomains_global_indices, eslocal number_of_subdomains ) {

	// *** Init the vector of domains *****************************************************
	domains_in_global_index.resize( number_of_subdomains ) ;

	domains.reserve(number_of_subdomains);
	for (eslocal d = 0; d < number_of_subdomains; d++) {
		domains.push_back( (Domain(configuration, instance, subdomains_global_indices[d], USE_HFETI)) );
		domains[d].domain_index = d;

		// Verbose level for K_plus
		if ( d == 0 && environment->MPIrank == 0) {
			domains[d].Kplus.msglvl = Info::report(LIBRARIES) ? 1 : 0;
		}
	}




	#pragma omp parallel for
	for (eslocal d = 0; d < domains.size(); d++ ) {
		domains_in_global_index[d] = subdomains_global_indices[d];
		domains[d].SetDomain();
	}

	//// *** Alocate temporarly vectors for Temporary vectors for Apply_A function *********
	//// *** - temporary vectors for work primal domain size *******************************
	x_prim_cluster1.resize( domains.size() );
	x_prim_cluster2.resize( domains.size() );
	x_prim_cluster3.resize( domains.size() );

	// *** Init all domains of the cluster ********************************************
	#pragma omp parallel for
	for (eslocal d = 0; d < number_of_subdomains; d++ ) {

		// *** Temporary vectors for work primal domain size *******************************
		x_prim_cluster1[d].resize( domains[d].domain_prim_size );
    	x_prim_cluster2[d].resize( domains[d].domain_prim_size );
    	x_prim_cluster3[d].resize( domains[d].domain_prim_size );


		domains[d].USE_KINV    	 		= USE_KINV;
		domains[d].USE_HFETI   	 		= USE_HFETI;
		domains[d].domain_index  		= d;
		domains[d].isOnACC  			= 0;

		// Verbose level for K_plus
		if ( d == 0 && environment->MPIrank == 0) {
			domains[d].Kplus.msglvl = 0;
		}

	}

	//// *** Set up the dual size ********************************************************
	dual_size = domains[0].B1.rows;

//	if (USE_HFETI == 1) {
//
//		// *** Alocate temporarly vectors for inter-cluster processing *********************
//		// *** - based on uncompressed matrix B0
//		tm1.resize(domains.size());
//		tm2.resize(domains.size());
//		tm3.resize(domains.size());
//
//		#pragma omp parallel for
//		for (size_t d = 0; d < domains.size(); d++) {
//			eslocal max_tmp_vec_size = domains[d].B0.cols;
//
//			if (domains[d].B0.rows > domains[d].B0.cols)
//				max_tmp_vec_size = domains[d].B0.rows;
//
//			tm1[d].resize( max_tmp_vec_size );
//			tm2[d].resize( max_tmp_vec_size );
//			tm3[d].resize( max_tmp_vec_size );
//		}
//		// *** END - Alocate temporarly vectors for inter-cluster processing *****************
//	}

	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++ ) {
		if (USE_KINV == 1 ) {
			domains[d].compressed_tmp.resize( domains[d].B1.I_row_indices.size(), 0);
			domains[d].compressed_tmp2.resize( domains[d].B1.I_row_indices.size(), 0);
		} else {
			domains[d].compressed_tmp.resize( 1, 0);
			domains[d].compressed_tmp2.resize( 1, 0);
		}
	}


	// local comm layer insode cluster

	SEQ_VECTOR <SEQ_VECTOR <eslocal> > & lambda_map_sub = instance->B1clustersMap;

	SEQ_VECTOR <eslocal> l_my_lamdas_indices ( lambda_map_sub.size() );
	for (size_t i = 0; i < lambda_map_sub.size(); i++)
		l_my_lamdas_indices[i] = lambda_map_sub[i][0];

	compressed_tmp    .resize( l_my_lamdas_indices.size(), 0 );	// TODO: might not be necessary

	// mapping/compression vector for domains
	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++) {
		for (size_t j = 0; j < domains[d].lambda_map_sub.size(); j++) {
			domains[d].my_lamdas_map_indices.insert(make_pair(domains[d].lambda_map_sub[j] ,j));
		}
	}

	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++) {
		if (domains[d].lambda_map_sub.size() > 0 ) {
			size_t i = 0;
			size_t j = 0;
			do
			{
				eslocal big_index   = l_my_lamdas_indices[i];
				eslocal small_index = domains[d].lambda_map_sub[j];

				if (big_index >  small_index) j++;

				if (big_index  < small_index) i++;

				if (big_index == small_index) {
					domains[d].lambda_map_sub_local.push_back(i);
					i++; j++;
				}


			} while ( i < l_my_lamdas_indices.size() && j < domains[d].lambda_map_sub.size() );
		}
	}
	//// *** END - Detection of affinity of lag. multipliers to specific subdomains ***************


	//// *** Compression of Matrix B1 to work with compressed lambda vectors *****************
//	ESLOG(MEMORY) << "B1 compression";
//	ESLOG(MEMORY) << "process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
//	ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++ ) {

		domains[d].B1_comp_dom.I_row_indices = domains[d].B1.I_row_indices;
		domains[d].B1_comp_dom.J_col_indices = domains[d].B1.J_col_indices;
		domains[d].B1_comp_dom.V_values      = domains[d].B1.V_values;

		domains[d].B1_comp_dom.rows = domains[d].B1.rows;
		domains[d].B1_comp_dom.cols = domains[d].B1.cols;
		domains[d].B1_comp_dom.nnz  = domains[d].B1.nnz;
		domains[d].B1_comp_dom.type = domains[d].B1.type;

		for (size_t j = 0; j < domains[d].B1_comp_dom.I_row_indices.size(); j++ ) {
			eslocal tmp_new = domains[d].my_lamdas_map_indices[domains[d].B1_comp_dom.I_row_indices [j] - 1] + 1;  // +1 means --> numbering from 1 in matrices
			domains[d].B1_comp_dom.I_row_indices [j] = tmp_new;
		}

		domains[d].B1_comp_dom.rows = domains[d].lambda_map_sub.size();
		domains[d].B1_comp_dom.ConvertToCSRwithSort( 1 );

		if (configuration.preconditioner == ESPRESO_PRECONDITIONER::DIRICHLET ||
        configuration.preconditioner == ESPRESO_PRECONDITIONER::SUPER_DIRICHLET) {

			domains[d].B1_comp_dom.MatTranspose(domains[d].B1t_DirPr);
			Esutils::removeDuplicity(domains[d].B1t_DirPr.CSR_I_row_indices);
			domains[d].B1t_DirPr.rows = domains[d].B1t_DirPr.CSR_I_row_indices.size() - 1;

			domains[d].B1t_Dir_perm_vec = domains[d].B1_comp_dom.CSR_J_col_indices;
			std::sort(domains[d].B1t_Dir_perm_vec.begin(), domains[d].B1t_Dir_perm_vec.end());
			Esutils::removeDuplicity(domains[d].B1t_Dir_perm_vec);

		}

//		if (configuration.regularization == REGULARIZATION::FIX_POINTS) {
//			domains[i].B1.Clear();
//		}

		domains[d].B1t.Clear();
		domains[d].my_lamdas_map_indices.clear();

	}
	//// *** END - Compression of Matrix B1 to work with compressed lambda vectors *************


	ESLOG(MEMORY) << "Lambdas end";
	ESLOG(MEMORY) << "process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
	ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

}

void ClusterBase::SetClusterPC( ) {

#ifdef READEX_LEVEL_1
	READEX_REGION_START(REG_Cluster_SetClusterPC, "Cluster--SetClusterPC", SCOREP_USER_REGION_TYPE_COMMON);
#endif

	SEQ_VECTOR <SEQ_VECTOR <eslocal> > & lambda_map_sub = instance->B1clustersMap;

	int MPIrank = environment->MPIrank;

	//// *** Detection of affinity of lag. multipliers to specific subdomains **************
	//// *** - will be used to compress vectors and matrices for higher efficiency

	ESLOG(MEMORY) << "Setting vectors for lambdas";
	ESLOG(MEMORY) << "process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
	ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

	my_lamdas_indices.resize( lambda_map_sub.size() );
	for (size_t i = 0; i < lambda_map_sub.size(); i++)
		my_lamdas_indices[i] = lambda_map_sub[i][0];


	SEQ_VECTOR< SEQ_VECTOR <eslocal> > lambdas_per_subdomain ( domains.size() * environment->MPIsize);
	my_lamdas_ddot_filter.resize( lambda_map_sub.size(), 0.0 );
	for (size_t i = 0; i < lambda_map_sub.size(); i++) {
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
	ESLOG(MEMORY) << "process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
	ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

	my_comm_lambdas_indices .resize(my_neighs.size());
	my_comm_lambdas			.resize(my_neighs.size());
	my_recv_lambdas			.resize(my_neighs.size());

	//ESLOG(MEMORY) << "1 process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";

	#pragma omp parallel for
	for (size_t i = 0; i < my_neighs.size(); i++) {
		my_comm_lambdas_indices[i] = lambdas_per_subdomain[my_neighs[i]];
		my_comm_lambdas[i].resize(my_comm_lambdas_indices[i].size());
		my_recv_lambdas[i].resize(my_comm_lambdas_indices[i].size());
	}

	//ESLOG(MEMORY) << "2 process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";

//	compressed_tmp    .resize( my_lamdas_indices.size(), 0 );
//	//compressed_tmp2   .resize( my_lamdas_indices.size(), 0 );

	//ESLOG(MEMORY) << "3 process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";


	// mapping/compression vector for cluster
	for (size_t i = 0; i <my_lamdas_indices.size(); i++)
		_my_lamdas_map_indices.insert(make_pair(my_lamdas_indices[i],i));

	//ESLOG(MEMORY) << "5 process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";

//	// mapping/compression vector for domains
//	#pragma omp parallel for
//	for (size_t i = 0; i < domains.size(); i++) {
//		for (size_t j = 0; j < domains[i].lambda_map_sub.size(); j++) {
//			domains[i].my_lamdas_map_indices.insert(make_pair(domains[i].lambda_map_sub[j] ,j));
//		}
//	}

	//ESLOG(MEMORY) << "6 process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";

//	#pragma omp parallel for
//	for (size_t d = 0; d < domains.size(); d++) {
//
//            if (domains[d].lambda_map_sub.size() > 0 ) {
//
//            	size_t i = 0;
//            	size_t j = 0;
//				do
//				{
//					eslocal big_index   = my_lamdas_indices[i];
//					eslocal small_index = domains[d].lambda_map_sub[j];
//
//					if (big_index >  small_index) j++;
//
//					if (big_index  < small_index) i++;
//
//					if (big_index == small_index) {
//						domains[d].lambda_map_sub_local.push_back(i);
//						i++; j++;
//					}
//
//
//				} while ( i < my_lamdas_indices.size() && j < domains[d].lambda_map_sub.size() );
//            }
//        }
//	//// *** END - Detection of affinity of lag. multipliers to specific subdomains ***************

	//ESLOG(MEMORY) << "7 process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";


	//// *** Create a vector of communication pattern needed for AllReduceLambdas function *******
	my_comm_lambdas_indices_comp.resize(my_neighs.size());
	#pragma omp parallel for
	for (size_t i = 0; i < my_neighs.size(); i++) {
		my_comm_lambdas_indices_comp[i].resize( lambdas_per_subdomain[my_neighs[i]].size() );
		for (size_t j = 0; j < lambdas_per_subdomain[my_neighs[i]].size(); j++ )
			my_comm_lambdas_indices_comp[i][j] = _my_lamdas_map_indices[lambdas_per_subdomain[my_neighs[i]][j]];
	}
	//// *** END - Create a vector of communication pattern needed for AllReduceLambdas function *

	//ESLOG(MEMORY) << "8 process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";


//	//// *** Compression of Matrix B1 to work with compressed lambda vectors *****************
//	ESLOG(MEMORY) << "B1 compression";
//	ESLOG(MEMORY) << "process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
//	ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";
//
//	#pragma omp parallel for
//	for (size_t i = 0; i < domains_in_global_index.size(); i++ ) {
//
//		domains[i].B1_comp_dom.I_row_indices = domains[i].B1.I_row_indices;
//		domains[i].B1_comp_dom.J_col_indices = domains[i].B1.J_col_indices;
//		domains[i].B1_comp_dom.V_values      = domains[i].B1.V_values;
//
//		domains[i].B1_comp_dom.rows = domains[i].B1.rows;
//		domains[i].B1_comp_dom.cols = domains[i].B1.cols;
//		domains[i].B1_comp_dom.nnz  = domains[i].B1.nnz;
//		domains[i].B1_comp_dom.type = domains[i].B1.type;
//
//		for (size_t j = 0; j < domains[i].B1_comp_dom.I_row_indices.size(); j++ ) {
//			eslocal tmp_new = domains[i].my_lamdas_map_indices[domains[i].B1_comp_dom.I_row_indices [j] - 1] + 1;  // numbering from 1 in matrix
//			domains[i].B1_comp_dom.I_row_indices [j] = tmp_new;									               // j + 1; // numbering matrix from 1
//		}
//
//		domains[i].B1_comp_dom.rows = domains[i].lambda_map_sub.size();
//		domains[i].B1_comp_dom.ConvertToCSRwithSort( 1 );
//
//		if (configuration.preconditioner == ESPRESO_PRECONDITIONER::DIRICHLET ||
//        configuration.preconditioner == ESPRESO_PRECONDITIONER::SUPER_DIRICHLET) {
//			domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_DirPr);
//			Esutils::removeDuplicity(domains[i].B1t_DirPr.CSR_I_row_indices);
////			auto last = std::unique(domains[i].B1t_DirPr.CSR_I_row_indices.begin(), domains[i].B1t_DirPr.CSR_I_row_indices.end());
////			domains[i].B1t_DirPr.CSR_I_row_indices.erase(last, domains[i].B1t_DirPr.CSR_I_row_indices.end());
//			domains[i].B1t_DirPr.rows = domains[i].B1t_DirPr.CSR_I_row_indices.size() - 1;
//
//			domains[i].B1t_Dir_perm_vec = domains[i].B1_comp_dom.CSR_J_col_indices;
//			std::sort(domains[i].B1t_Dir_perm_vec.begin(), domains[i].B1t_Dir_perm_vec.end());
//			Esutils::removeDuplicity(domains[i].B1t_Dir_perm_vec);
////			last = std::unique(domains[i].B1t_Dir_perm_vec.begin(), domains[i].B1t_Dir_perm_vec.end());
////			domains[i].B1t_Dir_perm_vec.erase(last, domains[i].B1t_Dir_perm_vec.end() );
//		}
//
////		if (configuration.regularization == REGULARIZATION::FIX_POINTS) {
////			domains[i].B1.Clear();
////		}
//
//		domains[i].B1t.Clear();
//		domains[i].my_lamdas_map_indices.clear();
//
//	}
//	//// *** END - Compression of Matrix B1 to work with compressed lambda vectors *************
//
//
//	ESLOG(MEMORY) << "Lambdas end";
//	ESLOG(MEMORY) << "process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
//	ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

#ifdef READEX_LEVEL_1
	READEX_REGION_STOP(REG_Cluster_SetClusterPC);
#endif

}

void ClusterBase::SetupPreconditioner ( ) {

	switch (configuration.preconditioner) {
	case ESPRESO_PRECONDITIONER::LUMPED:
		// nothing needs to be done
#ifdef BEM4I_TO_BE_REMOVED
		ESINFO(GLOBAL_ERROR) << "Memory efficient Lumped not possible for BEM, used fast Lumped --> MAGIC (or 5)";
#endif
		break;
	case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
		// nothing needs to be done
		break;
	case ESPRESO_PRECONDITIONER::DIRICHLET:
		CreateDirichletPrec(instance);
		break;
	case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
		CreateDirichletPrec(instance);
		break;
	case ESPRESO_PRECONDITIONER::MAGIC: // Fast Lumped
		#pragma omp parallel for
		for (size_t d = 0; d < domains.size(); d++) {
			if (domains[d]._RegMat.nnz > 0) {
#ifdef BEM4I_TO_BE_REMOVED
				domains[d]._RegMat.ConvertToCSR(1);
				domains[d].K.ConvertDenseToCSR(0);
				domains[d].Prec.MatAdd(domains[d].K, domains[d]._RegMat, 'N', -1);
				SEQ_VECTOR<eslocal>().swap( domains[d].K.CSR_I_row_indices );
				SEQ_VECTOR<eslocal>().swap( domains[d].K.CSR_J_col_indices );
				SEQ_VECTOR<double>().swap( domains[d].K.CSR_V_values );

				domains[d]._RegMat.ConvertToCOO(1);
#else
				domains[d]._RegMat.ConvertToCSR(1);
				domains[d].Prec.MatAdd(domains[d].K, domains[d]._RegMat, 'N', -1);
				domains[d]._RegMat.ConvertToCOO(1);
#endif
			} else {
				domains[d].Prec = domains[d].K;
			}
		}
		break;
	case ESPRESO_PRECONDITIONER::NONE:
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
	}

}


void ClusterBase::SetClusterHFETI () {
	// *** Create Matrices and allocate vectors for Hybrid FETI **************************

#ifdef READEX_LEVEL_1
	READEX_REGION_START(REG_Cluster_HFETIpreprocessing, "Cluster--HFETIpreprocessing", SCOREP_USER_REGION_TYPE_COMMON);
#endif

	if (USE_HFETI == 1) {

		if (domains.size() > 1) {

			TimeEval HFETI_prec_timing (" HFETI - preprocessing timing");
			if (Measure::report(CLUSTER)) { HFETI_prec_timing.totalTime.start(); };

			ESINFO(PROGRESS3) << "HFETI preprocessing start";

			TimeEvent B0_time("Compress B0 per cluster");
			if (Measure::report(CLUSTER)) { B0_time.start(); };

			CompressB0();

			if (Measure::report(CLUSTER)) { B0_time.end(); B0_time.printStatMPI(); HFETI_prec_timing.addEvent(B0_time); }

			TimeEvent G0_time("Create G0 per cluster");
			if (Measure::report(CLUSTER)) { G0_time.start(); }

			CreateG0();

			if (Measure::report(CLUSTER)) { G0_time.end(); 	G0_time.printStatMPI();	HFETI_prec_timing.addEvent(G0_time); }


			TimeEvent F0_time("Create F0 per cluster");
			if (Measure::report(CLUSTER)) { F0_time.start();}

			CreateF0();

			if (Measure::report(CLUSTER)) { F0_time.end(); HFETI_prec_timing.addEvent(F0_time); }


			TimeEvent Sa_time("Create Salfa per cluster");
			if (Measure::report(CLUSTER)) { Sa_time.start(); }

			CreateSa();

			if (Measure::report(CLUSTER)) { Sa_time.end(); HFETI_prec_timing.addEvent(Sa_time); }

			if (Measure::report(CLUSTER)) { HFETI_prec_timing.totalTime.end(); HFETI_prec_timing.printStatsMPI(); }

		} else {

			//ESINFO(ALWAYS) << Info::TextColor::YELLOW
			std::cout << "Cluster " << cluster_global_index << " on MPI rank " << environment->MPIrank << " has only one domain -> Using TFETI" << std::endl;
			USE_HFETI = 0;

		}

	}

#ifdef READEX_LEVEL_1
	READEX_REGION_STOP(REG_Cluster_HFETIpreprocessing);
#endif

	// *** END - Create Matrices for Hybrid FETI *****************************************
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

	//ESINFO(PROGRESS3) << "K+ multiply HFETI";
	mkl_set_num_threads(1);

	if (Measure::report(CLUSTER)) { cluster_time.totalTime.start(); }

	if (Measure::report(CLUSTER)) { vec_fill_time.start(); }
	fill(vec_g0.begin(), vec_g0.end(), 0); // reset entire vector to 0
	//fill(vec_e0.begin(), vec_e0.end(), 0); // reset entire vector to 0
	if (Measure::report(CLUSTER)) { vec_fill_time.end(); }

	// loop over domains in the cluster
	if (Measure::report(CLUSTER)) { loop_1_1_time.start(); }
	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++)
	{
		domains[d].B0Kplus_comp.DenseMatVec(x_in[d], tm2[d]);			// g0 - with comp B0Kplus
		if (SYMMETRIC_SYSTEM) {
			domains[d].Kplus_R.DenseMatVec(x_in[d], tm3[d], 'T');			// e0
		} else {
			domains[d].Kplus_R2.DenseMatVec(x_in[d], tm3[d], 'T');			// e0
		}
	}
	if (Measure::report(CLUSTER)) { loop_1_1_time.end();}

	if (Measure::report(CLUSTER)) { loop_1_2_time.start(); }
	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++)
	{

		SEQ_VECTOR <eslocal> kerindices (domains.size() + 1, 0);
		kerindices[0] = 0;
		for (size_t k = 1; k < kerindices.size(); k++) {
			kerindices[k] = kerindices[k-1] + domains[k-1].Kplus_R.cols;
		}

		eslocal e0_start	= kerindices[d];
		eslocal e0_end		= kerindices[d+1];

		//eslocal e0_start	=  d	* domains[d].Kplus_R.cols;
		//eslocal e0_end		= (d+1) * domains[d].Kplus_R.cols;

		for (eslocal i = e0_start; i < e0_end; i++ )
			vec_e0[i] = - tm3[d][i - e0_start];
	}


	for (size_t d = 0; d < domains.size(); d++)
		for (eslocal i = 0; i < domains[d].B0Kplus_comp.rows; i++)
			vec_g0[ domains[d].B0_comp_map_vec[i] - 1 ] += tm2[d][i];

	// end loop over domains
	if (Measure::report(CLUSTER)) { loop_1_2_time.end(); }

	mkl_set_num_threads(PAR_NUM_THREADS);
	if (Measure::report(CLUSTER)) { clusCP_time.start(); }


//	for (int i = 0; i < vec_g0.size(); i++)
//	printf (       "Test probe 1: %d norm = %1.30f \n", i, vec_g0[i] );

	if (Measure::report(CLUSTER)) { clus_F0_1_time.start(); }
	F0.Solve(vec_g0, tm1[0], 0, 0);
	if (Measure::report(CLUSTER)) { clus_F0_1_time.end(); }

//	for (int i = 0; i < tm1[0].size(); i++)
//	printf (       "Test probe 2: %d norm = %1.30f \n", i, tm1[0][i] );

	if (Measure::report(CLUSTER)) { clus_G0_time.start(); }
	if (SYMMETRIC_SYSTEM) {
		G0.MatVec(tm1[0], tm2[0], 'N');
	} else {
		G02.MatVec(tm1[0], tm2[0], 'N');
	}

	if (Measure::report(CLUSTER)) { clus_G0_time.end(); }

//	for (int i = 0; i < tm1[0].size(); i++)
//	printf (       "Test probe 3: %d norm = %1.30f \n", i, tm1[0][i] );

	#pragma omp parallel for
	for (size_t i = 0; i < vec_e0.size(); i++)
		tm2[0][i] = tm2[0][i] - vec_e0[i];
	//cblas_daxpy(vec_e0.size(), -1.0, &vec_e0[0], 1, &tm2[0][0], 1);

	if (Measure::report(CLUSTER)) { clus_Sa_time.start(); }
//#ifdef SPARSE_SA
//	 Sa.Solve(tm2[0], vec_alfa,0,0);
//#else
//	eslocal nrhs = 1;
//	Sa_dense.Solve(tm2[0], vec_alfa, nrhs);
//#endif

	switch (configuration.SAsolver) {
	case ESPRESO_SASOLVER::CPU_SPARSE:
		Sa.Solve(tm2[0], vec_alfa, 0, 0);
		break;
	case ESPRESO_SASOLVER::CPU_DENSE:
		Sa_dense_cpu.Solve(tm2[0], vec_alfa, 1);
		break;
	case ESPRESO_SASOLVER::ACC_DENSE:
		Sa_dense_acc.Solve(tm2[0], vec_alfa, 1);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented S alfa solver.";
	}


	if (Measure::report(CLUSTER)) { clus_Sa_time.end(); }

//		for (int i = 0; i < vec_alfa.size(); i++)
//		printf (       "Test probe 4: %d norm = %1.30f \n", i, vec_alfa[i] );

	if (Measure::report(CLUSTER)) { clus_G0t_time.start(); }
	G0.MatVec(vec_alfa, tm1[0], 'T'); 	// lambda
	if (Measure::report(CLUSTER)) {clus_G0t_time.end(); }

//		for (int i = 0; i < tm1[0].size(); i++)
//		printf (       "Test probe 5: %d norm = %1.30f \n", i, tm1[0][i] );

	#pragma omp parallel for
	for (size_t i = 0; i < vec_g0.size(); i++)
		tm1[0][i] = vec_g0[i] - tm1[0][i];


	if (Measure::report(CLUSTER)) { clus_F0_2_time.start(); }
	F0.Solve(tm1[0], vec_lambda,0,0);
	if (Measure::report(CLUSTER)) { clus_F0_2_time.end(); }

	if (Measure::report(CLUSTER)) { clusCP_time.end(); }

//	for (int i = 0; i < vec_lambda.size(); i++)
//	printf (       "Test probe 6: %d norm = %1.30f \n", i, vec_lambda[i] );

	// Kplus_x
	mkl_set_num_threads(1);
	if (Measure::report(CLUSTER)) { loop_2_1_time.start(); }
	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++)
	{
		eslocal domain_size = domains[d].domain_prim_size;

		SEQ_VECTOR < double > tmp_vec (domains[d].B0_comp_map_vec.size(), 0);
		for (size_t i = 0; i < domains[d].B0_comp_map_vec.size(); i++)
			tmp_vec[i] = vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;
		//domains[d].B0t_comp.MatVec(tmp_vec, tm1[d], 'N');
		domains[d].B0_comp.MatVec(tmp_vec, tm1[d], 'T');

		for (eslocal i = 0; i < domain_size; i++)
			tm1[d][i] = x_in[d][i] - tm1[d][i];

		domains[d].multKplusLocal(tm1[d] , tm2[d]);

		SEQ_VECTOR <eslocal> kerindices (domains.size() + 1, 0);
		kerindices[0] = 0;
		for (size_t k = 1; k < kerindices.size(); k++) {
			kerindices[k] = kerindices[k-1] + domains[k-1].Kplus_R.cols;
		}

		eslocal e0_start	= kerindices[d];
		//eslocal e0_start	=  d	* domains[d].Kplus_R.cols;

		domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);

		for (eslocal i = 0; i < domain_size; i++)
			x_in[d][i] = tm2[d][i] + tm3[d][i];

		//ESINFO(PROGRESS3) << Info::plain() << ".";
	}
	//ESINFO(PROGRESS3);
	if (Measure::report(CLUSTER)) { loop_2_1_time.end(); }

	if (Measure::report(CLUSTER)) { cluster_time.totalTime.end(); }
}

void ClusterBase::multKplusGlobal_Kinv( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in ) {

	mkl_set_num_threads(1);
	if (Measure::report(CLUSTER)) { cluster_time.totalTime.start(); }

	if (Measure::report(CLUSTER)) { vec_fill_time.start(); }
	fill(vec_g0.begin(), vec_g0.end(), 0); // reset entire vector to 0
	//fill(vec_e0.begin(), vec_e0.end(), 0); // reset entire vector to 0
	if (Measure::report(CLUSTER)) { vec_fill_time.end(); }

	// loop over domains in the cluster
	if (Measure::report(CLUSTER)) { loop_1_1_time.start(); }
	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++)
	{
   		domains[d].B0Kplus_comp.DenseMatVec(x_in[d], tm2[d]);			// g0 - with comp B0Kplus
		domains[d].Kplus_R.     DenseMatVec(x_in[d], tm3[d], 'T');		// e0
	}
	if (Measure::report(CLUSTER)) { loop_1_1_time.end(); }

	if (Measure::report(CLUSTER)) { loop_1_2_time.start(); }
	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++)
	{
		eslocal e0_start	=  d	* domains[d].Kplus_R.cols;
		eslocal e0_end		= (d+1) * domains[d].Kplus_R.cols;

		for (eslocal i = e0_start; i < e0_end; i++ )
			vec_e0[i] = - tm3[d][i - e0_start];
	}

	for (size_t d = 0; d < domains.size(); d++)
		for (eslocal i = 0; i < domains[d].B0Kplus_comp.rows; i++)
			vec_g0[ domains[d].B0_comp_map_vec[i] - 1 ] += tm2[d][i];

	// end loop over domains
	if (Measure::report(CLUSTER)) { loop_1_2_time.end(); }

	mkl_set_num_threads(PAR_NUM_THREADS);
	if (Measure::report(CLUSTER)) { clusCP_time.start(); }

	if (Measure::report(CLUSTER)) { clus_F0_1_time.start(); }
	F0.Solve(vec_g0, tm1[0], 0, 0);
	if (Measure::report(CLUSTER)) { clus_F0_1_time.end(); }

	if (Measure::report(CLUSTER)) { clus_G0_time.start(); }
	G0.MatVec(tm1[0], tm2[0], 'N');
	if (Measure::report(CLUSTER)) { clus_G0_time.end(); }

	#pragma omp parallel for
	for (size_t i = 0; i < vec_e0.size(); i++) {
		tm2[0][i] = tm2[0][i] - vec_e0[i];
	}
	//cblas_daxpy(vec_e0.size(), -1.0, &vec_e0[0], 1, &tm2[0][0], 1);

	if (Measure::report(CLUSTER)) { clus_Sa_time.start(); }
// #ifdef SPARSE_SA
//    Sa.Solve(tm2[0], vec_alfa,0,0);
// #else
//    eslocal nrhs = 1;
//    Sa_dense.Solve(tm2[0], vec_alfa, nrhs);
// #endif

	switch (configuration.SAsolver) {
	case ESPRESO_SASOLVER::CPU_SPARSE:
		Sa.Solve(tm2[0], vec_alfa, 0, 0);
		break;
	case ESPRESO_SASOLVER::CPU_DENSE:
		Sa_dense_cpu.Solve(tm2[0], vec_alfa, 1);
		break;
	case ESPRESO_SASOLVER::ACC_DENSE:
		Sa_dense_acc.Solve(tm2[0], vec_alfa, 1);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented S alfa solver.";
	}

	if (Measure::report(CLUSTER)) { clus_Sa_time.end(); }

	if (Measure::report(CLUSTER)) { clus_G0t_time.start(); }
	G0.MatVec(vec_alfa, tm1[0], 'T'); 	// lambda
	if (Measure::report(CLUSTER)) { clus_G0t_time.end(); }

	#pragma omp parallel for
	for (size_t i = 0; i < vec_g0.size(); i++) {
		tm1[0][i] = vec_g0[i] - tm1[0][i];
	}


	if (Measure::report(CLUSTER)) { clus_F0_2_time.start(); }
	F0.Solve(tm1[0], vec_lambda,0,0);
	if (Measure::report(CLUSTER)) { clus_F0_2_time.end(); }

	if (Measure::report(CLUSTER)) { clusCP_time.end(); }


	// Kplus_x
	mkl_set_num_threads(1);
	if (Measure::report(CLUSTER)) { loop_2_1_time.start(); }
	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++)
	{
		eslocal domain_size = domains[d].domain_prim_size;

		SEQ_VECTOR < double > tmp_vec (domains[d].B0_comp_map_vec.size(), 0.0);
		for (size_t i = 0; i < domains[d].B0_comp_map_vec.size(); i++) {
			tmp_vec[i] = -1.0 * vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;
		}
		domains[d].B0Kplus_comp.DenseMatVec(tmp_vec, tm2[d], 'T' );

		eslocal e0_start	=  d	* domains[d].Kplus_R.cols;
		domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);

		for (eslocal i = 0; i < domain_size; i++)
			x_in[d][i] = tm2[d][i] + tm3[d][i];

	}
	if (Measure::report(CLUSTER)) { loop_2_1_time.end(); }

	if (Measure::report(CLUSTER)) { cluster_time.totalTime.end(); }
}

void ClusterBase::multKplusGlobal_Kinv_2( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in ) {

	mkl_set_num_threads(1);
	if (Measure::report(CLUSTER)) { cluster_time.totalTime.start(); }

	if (Measure::report(CLUSTER)) { vec_fill_time.start(); }
	fill(vec_g0.begin(), vec_g0.end(), 0); // reset entire vector to 0
	//fill(vec_e0.begin(), vec_e0.end(), 0); // reset entire vector to 0
	if (Measure::report(CLUSTER)) { vec_fill_time.end(); }

	// loop over domains in the cluster
	if (Measure::report(CLUSTER)) { loop_1_1_time.start(); }
	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++)
	{
   		//domains[d].B0Kplus_comp.DenseMatVec(x_in[d], tm2[d]);				// g0 - with comp B0Kplus
		//domains[d].Kplus_R.     DenseMatVec(x_in[d], tm3[d], 'T');		// e0
		domains[d].B0KplusB1_comp .DenseMatVec(x_in[d], tm2[d], 'N');		// g0 - with comp B0Kplus
		domains[d].Kplus_R_B1_comp.DenseMatVec(x_in[d], tm3[d], 'N');		// e0
	}
	if (Measure::report(CLUSTER)) { loop_1_1_time.end(); }

	if (Measure::report(CLUSTER)) { loop_1_2_time.start(); }
	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++)
	{
		eslocal e0_start	=  d	* domains[d].Kplus_R.cols;
		eslocal e0_end		= (d+1) * domains[d].Kplus_R.cols;

		for (eslocal i = e0_start; i < e0_end; i++ )
			vec_e0[i] = - tm3[d][i - e0_start];
	}

	for (size_t d = 0; d < domains.size(); d++)
		for (eslocal i = 0; i < domains[d].B0Kplus_comp.rows; i++)
			vec_g0[ domains[d].B0_comp_map_vec[i] - 1 ] += tm2[d][i];

	// end loop over domains
	if (Measure::report(CLUSTER)) { loop_1_2_time.end(); }

	mkl_set_num_threads(PAR_NUM_THREADS);
	if (Measure::report(CLUSTER)) { clusCP_time.start(); }

	if (Measure::report(CLUSTER)) { clus_F0_1_time.start(); }
	F0.Solve(vec_g0, tm1[0], 0, 0);
	if (Measure::report(CLUSTER)) { clus_F0_1_time.end(); }

	if (Measure::report(CLUSTER)) { clus_G0_time.start(); }
	G0.MatVec(tm1[0], tm2[0], 'N');
	if (Measure::report(CLUSTER)) { clus_G0_time.end(); }

	#pragma omp parallel for
	for (size_t i = 0; i < vec_e0.size(); i++)
		tm2[0][i] = tm2[0][i] - vec_e0[i];
	//cblas_daxpy(vec_e0.size(), -1.0, &vec_e0[0], 1, &tm2[0][0], 1);

	if (Measure::report(CLUSTER)) { clus_Sa_time.start(); }
//	Sa.Solve(tm2[0], vec_alfa,0,0);

	switch (configuration.SAsolver) {
	case ESPRESO_SASOLVER::CPU_SPARSE:
		Sa.Solve(tm2[0], vec_alfa, 0, 0);
		break;
	case ESPRESO_SASOLVER::CPU_DENSE:
		Sa_dense_cpu.Solve(tm2[0], vec_alfa, 1);
		break;
	case ESPRESO_SASOLVER::ACC_DENSE:
		Sa_dense_acc.Solve(tm2[0], vec_alfa, 1);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented S alfa solver.";
	}

	if (Measure::report(CLUSTER)) { clus_Sa_time.end(); }

	if (Measure::report(CLUSTER)) { clus_G0t_time.start(); }
	G0.MatVec(vec_alfa, tm1[0], 'T'); 	// lambda
	if (Measure::report(CLUSTER)) { clus_G0t_time.end(); }

	#pragma omp parallel for
	for (size_t i = 0; i < vec_g0.size(); i++) {
		tm1[0][i] = vec_g0[i] - tm1[0][i];
	}


	if (Measure::report(CLUSTER)) { clus_F0_2_time.start(); }
	F0.Solve(tm1[0], vec_lambda,0,0);
	if (Measure::report(CLUSTER)) { clus_F0_2_time.end(); }

	if (Measure::report(CLUSTER)) { clusCP_time.end(); }


	// Kplus_x
	mkl_set_num_threads(1);
	if (Measure::report(CLUSTER)) { loop_2_1_time.start(); }
	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++)
	{
		eslocal domain_size = domains[d].domain_prim_size;

		SEQ_VECTOR < double > tmp_vec (domains[d].B0_comp_map_vec.size(), 0.0);
		for (size_t i = 0; i < domains[d].B0_comp_map_vec.size(); i++) {
			tmp_vec[i] = -1.0 * vec_lambda[domains[d].B0_comp_map_vec[i] - 1] ;
		}

		//domains[d].B0Kplus_comp.DenseMatVec(tmp_vec, tm2[d], 'T' );
		domains[d].B0KplusB1_comp .DenseMatVec(tmp_vec, tm2[d], 'T');
		domains[d].B0KplusB1_comp .MatVec(tmp_vec, tm2[d], 'T');

		eslocal e0_start	=  d	* domains[d].Kplus_R.cols;

		//domains[d].Kplus_R.DenseMatVec(vec_alfa, tm3[d],'N', e0_start);
		domains[d].Kplus_R_B1_comp.DenseMatVec(vec_alfa, tm3[d], 'T', e0_start);
		//domains[d].Kplus_R_B1_comp.MatVec(vec_alfa, tm3[d], 'T', e0_start);

		for (eslocal i = 0; i < domain_size; i++)
			x_in[d][i] = tm2[d][i] + tm3[d][i];

	}
	if (Measure::report(CLUSTER)) { loop_2_1_time.end(); }

	if (Measure::report(CLUSTER)) { cluster_time.totalTime.end(); }
}


////backup March 31 2015
//void ClusterBase::multKplusGlobal_l(SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {
//
//	//eslocal MPIrank;
//	//MPI_Comm_rank (environment->MPICommunicator, &MPIrank);
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

	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++) {

		//TODO: Need fix - dual copy - prepisuju data asembleru
		domains[d].B0 = instance->B0[domains[d].domain_global_index];
		domains[d].B0.type = 'G';
		domains[d].B0.ConvertToCSRwithSort(1);


		domains[d].B0.MatTranspose(domains[d].B0t);
		domains[d].B0_comp = domains[d].B0;

		for (eslocal i = 0; i < domains[d].B0_comp.rows; i++) {
			if (domains[d].B0_comp.CSR_I_row_indices[i] != domains[d].B0_comp.CSR_I_row_indices[i + 1]) {
				domains[d].B0_comp_map_vec.push_back(i + 1);
			}
		}

		size_t unique = 0;
		for (size_t i = 1; i < domains[d].B0_comp.CSR_I_row_indices.size(); i++) {
			if (domains[d].B0_comp.CSR_I_row_indices[unique] != domains[d].B0_comp.CSR_I_row_indices[i]) {
				domains[d].B0_comp.CSR_I_row_indices[++unique] = domains[d].B0_comp.CSR_I_row_indices[i];
			}
		}
		domains[d].B0_comp.rows = unique;
		domains[d].B0_comp.CSR_I_row_indices.resize(unique + 1);

		// WARNING: There is a problem with 'std::unique'.
		// WARNING: Combination of 'Cilk' and '-O2' results in memory error in 'std::unique'.
		//
		//auto it = std::unique(domains[d].B0_comp.CSR_I_row_indices.begin(), domains[d].B0_comp.CSR_I_row_indices.end());
		//domains[d].B0_comp.rows = std::distance(domains[d].B0_comp.CSR_I_row_indices.begin(), it) - 1;
		//domains[d].B0_comp.CSR_I_row_indices.resize(domains[d].B0_comp.rows + 1);

		domains[d].B0_comp.MatTranspose(domains[d].B0t_comp);
	}
}

void ClusterBase::CreateG0() {

	mkl_set_num_threads(1);

	SEQ_VECTOR <SparseMatrix> G0LocalTemp( domains.size() );

	#pragma omp parallel for
	for (size_t i = 0; i<domains.size(); i++) {
		domains[i].Kplus_R.ConvertDenseToCSR(0);

		G0LocalTemp[i].MatMat(domains[i].B0, 'N', domains[i].Kplus_R );
		G0LocalTemp[i].MatTranspose(-1.0);

		SEQ_VECTOR<eslocal>().swap( domains[i].Kplus_R.CSR_I_row_indices );
		SEQ_VECTOR<eslocal>().swap( domains[i].Kplus_R.CSR_J_col_indices );
		SEQ_VECTOR<double> ().swap( domains[i].Kplus_R.CSR_V_values );
	}

	for (size_t i = 0; i<domains.size(); i++) {
		G0.MatAppend(G0LocalTemp[i]);
		G0LocalTemp[i].Clear();
	}

	if (!SYMMETRIC_SYSTEM) {

		SEQ_VECTOR <SparseMatrix> G0LocalTemp2( domains.size() );

		#pragma omp parallel for
		for (size_t i = 0; i<domains.size(); i++) {
			domains[i].Kplus_R2.ConvertDenseToCSR(0);

			G0LocalTemp2[i].MatMat(domains[i].B0, 'N', domains[i].Kplus_R2 );
			G0LocalTemp2[i].MatTranspose(-1.0);

			SEQ_VECTOR<eslocal>().swap( domains[i].Kplus_R2.CSR_I_row_indices );
			SEQ_VECTOR<eslocal>().swap( domains[i].Kplus_R2.CSR_J_col_indices );
			SEQ_VECTOR<double> ().swap( domains[i].Kplus_R2.CSR_V_values );
		}

		for (size_t i = 0; i<domains.size(); i++) {
			G02.MatAppend(G0LocalTemp2[i]);
			G0LocalTemp2[i].Clear();
		}

		if (environment->print_matrices) {
			SparseMatrix tmpG02 = G02;
			std::ofstream osG0(Logging::prepareFile("G02"));
			osG0 <<  tmpG02;
			osG0.close();
		}


	}


	if (environment->print_matrices) {
		SparseMatrix tmpG01 = G0;
		std::ofstream osG0(Logging::prepareFile("G01"));
		osG0 <<  tmpG01;
		osG0.close();
	}

	vec_g0.resize(G0.cols);
	vec_e0.resize(G0.rows);

}

void ClusterBase::CreateF0() {

	 TimeEval F0_timing (" HFETI - F0 preprocessing timing");
	 if (Measure::report(CLUSTER)) { F0_timing.totalTime.start(); }

	mkl_set_num_threads(1);

	int MPIrank; MPI_Comm_rank (environment->MPICommunicator, &MPIrank);

	SEQ_VECTOR <SparseMatrix> tmpF0v (domains.size());

	ESINFO(PROGRESS3) << "HFETI - Create F0";

#ifdef READEX_LEVEL_1
	READEX_REGION_START(REG_Cluster_CreateF0_AssembleF0, "Cluster--CreateF0-AssembleF0", SCOREP_USER_REGION_TYPE_COMMON);
#endif

	 TimeEvent solve_F0_time("B0 compression; F0 multiple InitialCondition solve");
	 if (Measure::report(CLUSTER)) { solve_F0_time.start(); }

	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++) {

		if (MPIrank == 0 && d == 0)
			domains[d].Kplus.msglvl=0;
		else
			domains[d].Kplus.msglvl=0;

		// FO solve in double in case K is in single
		if (
				configuration.F0_precision == ESPRESO_F0SOLVER_PRECISION::DOUBLE
				&& (configuration.Ksolver == ESPRESO_KSOLVER::DIRECT_SP
				|| configuration.Ksolver == ESPRESO_KSOLVER::DIRECT_MP )
			) {
			SparseSolverMKL Ktmp;
			Ktmp.ImportMatrix_wo_Copy(domains[d].K);
			std::stringstream ss;
			ss << "Create F0 -> rank: " << environment->MPIrank << ", subdomain: " << d;
			Ktmp.Factorization(ss.str());
			Ktmp.SolveMat_Dense(domains[d].B0t_comp, domains[d].B0Kplus_comp);
			domains[d].B0Kplus = domains[d].B0Kplus_comp;
		} else {
			// F0 uses same precision as K
			if (SYMMETRIC_SYSTEM) {
				if(configuration.mp_pseudoinverse) {
					domains[d].Kplus_R.ConvertDenseToCSR(0);

					SparseMatrix tmR, tmR2;
					SparseMatrix B0tm;

					B0tm = domains[d].B0t_comp;
					tmR.MatMat(domains[d].Kplus_R,'T', domains[d].B0t_comp);
					tmR2.MatMat(domains[d].Kplus_R,'N',tmR);
					B0tm.MatAddInPlace(tmR2,'N', -1.0);

					domains[d].Kplus.SolveMat_Dense(B0tm, domains[d].B0Kplus_comp);

					tmR.Clear();
					tmR2.Clear();
					B0tm.Clear();
					B0tm = domains[d].B0Kplus_comp;
					tmR.MatMat(domains[d].Kplus_R,'T', domains[d].B0Kplus_comp);
					tmR2.MatMat(domains[d].Kplus_R,'N',tmR);
					B0tm.MatAddInPlace(tmR2,'N', -1.0);

					domains[d].B0Kplus = B0tm;
					domains[d].B0Kplus_comp = B0tm;



				} else {
					domains[d].Kplus.SolveMat_Dense(domains[d].B0t_comp, domains[d].B0Kplus_comp);
					domains[d].B0Kplus = domains[d].B0Kplus_comp;

//					ESINFO(PROGRESS1) << domains[d].B0t_comp.SpyText();
				}
			} else {
				//TODO: The Klus.Solve - does not have to be called twice here - can be done with Transpose
				//TODO: Alex
				eslocal set_bckp = domains[d].Kplus.iparm[11];
				domains[d].Kplus.iparm[11] = 2;
				domains[d].Kplus.SolveMat_Dense(domains[d].B0t_comp, domains[d].B0Kplus_comp);
				domains[d].Kplus.iparm[11] = set_bckp;
				domains[d].Kplus.SolveMat_Dense(domains[d].B0t_comp, domains[d].B0Kplus);
			}
		}

		domains[d].B0t_comp.Clear();


		domains[d].B0Kplus_comp.MatTranspose();
		domains[d].B0Kplus_comp.ConvertCSRToDense(1);

		for (size_t i = 0; i < domains[d].B0Kplus.CSR_J_col_indices.size() - 1; i++) {
			domains[d].B0Kplus.CSR_J_col_indices[i] = domains[d].B0_comp_map_vec [ domains[d].B0Kplus.CSR_J_col_indices[i] - 1 ];
		}

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
		ESINFO(PROGRESS3) << Info::plain() << ".";
	}

	ESINFO(PROGRESS3);

	if (Measure::report(CLUSTER)) { solve_F0_time.end(); }
	if (Measure::report(CLUSTER)) { solve_F0_time.printStatMPI(); }
	if (Measure::report(CLUSTER)) {F0_timing.addEvent(solve_F0_time);}

	TimeEvent reduction_F0_time("F0 reduction time");
	if (Measure::report(CLUSTER)) {reduction_F0_time.start();}

	for (size_t j = 1; j < tmpF0v.size(); j *= 2) {
		#pragma omp parallel for
		for (size_t i = 0; i <= tmpF0v.size() / (2 * j); i++) {
			if (i * 2 * j + j < tmpF0v.size()) {
				tmpF0v[i * 2 * j].MatAddInPlace( tmpF0v[i * 2 * j + j], 'N', 1.0 );
				tmpF0v[i * 2 * j + j].Clear();
			}
		}
	}
	F0_Mat = tmpF0v[0];

	if (environment->print_matrices) {
		SparseMatrix tmpF0 = F0_Mat;
		std::ofstream osF0(Logging::prepareFile("F0"));
		osF0 <<  tmpF0;
		osF0.close();
	}

	if (Measure::report(CLUSTER)) { reduction_F0_time.end(); reduction_F0_time.printStatMPI(); F0_timing.addEvent(reduction_F0_time); }


	TimeEvent fact_F0_time("B0 Kplus Factorization ");
	if (Measure::report(CLUSTER)) { fact_F0_time.start();}

	mkl_set_num_threads(PAR_NUM_THREADS);

#ifdef READEX_LEVEL_1
	READEX_REGION_STOP(REG_Cluster_CreateF0_AssembleF0);
#endif

#ifdef READEX_LEVEL_1
	READEX_REGION_START(REG_Cluster_CreateF0_FactF0, "Cluster--CreateF0-FactF0", SCOREP_USER_REGION_TYPE_COMMON);
#endif

	if (SYMMETRIC_SYSTEM) {
		F0_Mat.RemoveLower();
	}

	F0_Mat.mtype = mtype;
	F0.ImportMatrix(F0_Mat);

//	bool PARDISO_SC = true;
//	if (!PARDISO_SC)
//		F0_fast.ImportMatrix(F0_Mat);

	//F0_Mat.Clear();
	F0.SetThreaded();
	std::stringstream ss;
	ss << "F0 -> rank: " << environment->MPIrank;
	F0.Factorization(ss.str());

	mkl_set_num_threads(1);

	if (MPIrank == 0) F0.msglvl = 0;

	if (Measure::report(CLUSTER)) {fact_F0_time.end(); fact_F0_time.printStatMPI(); F0_timing.addEvent(fact_F0_time);}

#ifdef READEX_LEVEL_1
	READEX_REGION_STOP(REG_Cluster_CreateF0_FactF0);
#endif

	if (Measure::report(CLUSTER)) { F0_timing.totalTime.end(); }
	//TODO: Need fix - for detailed measurements
//	if (min_numClusters_per_MPI > cluster_local_index)
//		F0_timing.printStatsMPI();

	// *** POZOR **************************************************************
	#pragma omp parallel for
	for (size_t d = 0; d<domains.size(); d++) {
		domains[d].B0.Clear();
		domains[d].B0t.Clear();
	}

	vec_lambda.resize(F0.m_Kplus_size);

};

void ClusterBase::CreateSa() {

#ifdef READEX_LEVEL_1
	READEX_REGION_START(REG_Cluster_CreateSa, "Cluster--CreateSa", SCOREP_USER_REGION_TYPE_COMMON);
#endif

	 TimeEval Sa_timing (" HFETI - Salfa preprocessing timing"); if (Measure::report(CLUSTER)) { Sa_timing.totalTime.start(); }

//	bool PARDISO_SC = true;
	bool get_kernel_from_mesh = configuration.regularization == REGULARIZATION::FIX_POINTS;
	int  MPIrank = environment->MPIrank;

	MKL_Set_Num_Threads(PAR_NUM_THREADS);


	SparseMatrix Salfa, tmpM;


	 TimeEvent G0trans_Sa_time("G0 transpose"); if (Measure::report(CLUSTER)) { G0trans_Sa_time.start(); }
	SparseMatrix G0t;
	SparseMatrix G02t;
	G0.MatTranspose(G0t);
	if (!SYMMETRIC_SYSTEM) {
		G02.MatTranspose(G02t);
	}
	 if (Measure::report(CLUSTER)) { G0trans_Sa_time.end(); G0trans_Sa_time.printStatMPI(); Sa_timing.addEvent(G0trans_Sa_time); }

#ifdef READEX_LEVEL_1
	READEX_REGION_START(REG_Cluster_CreateSa_SolveF0vG0, "Cluster--CreateSa-SolveF0vG0", SCOREP_USER_REGION_TYPE_COMMON);
#endif

	 TimeEvent G0solve_Sa_time("SolveMatF with G0t as InitialCondition"); if (Measure::report(CLUSTER)) { G0solve_Sa_time.start(); }
	SparseSolverMKL Salfa_SC_solver;
	if (MPIrank == 0) Salfa_SC_solver.msglvl = Info::report(LIBRARIES) ? 1 : 0;
	if (SYMMETRIC_SYSTEM) {
		Salfa_SC_solver.Create_SC_w_Mat( F0_Mat, G0t, Salfa, true, 0 );
		Salfa.ConvertDenseToCSR(1);
		Salfa.RemoveLower();
	} else {
		Salfa_SC_solver.Create_non_sym_SC_w_Mat( F0_Mat, G0t, G02t, Salfa, true, 0 );
		//TODO: SC je dive transponopvany
		//Salfa.ConvertDenseToCSR(1);
		//Salfa.MatTranspose();
	}
	 if (Measure::report(CLUSTER)) { G0solve_Sa_time.end(); G0solve_Sa_time.printStatMPI(); Sa_timing.addEvent(G0solve_Sa_time); }

	Salfa.mtype = this->mtype;
	F0_Mat.Clear();

	if (environment->print_matrices) {
		std::ofstream osSa(Logging::prepareFile("Salfa"));
		osSa << Salfa;
		osSa.close();
	}

//	if (!PARDISO_SC) {
//		if (MPIrank == 0) { F0_fast.msglvl = Info::report(LIBRARIES) ? 1 : 0; }
//
//		//SolaveMatF is obsolete - use Schur Complement Instead
//		F0_fast.SolveMatF(G0t,tmpM, true);
//		if (MPIrank == 0) F0_fast.msglvl = 0;
//
//		Salfa.MatMat(G0, 'N', tmpM);
//		Salfa.RemoveLower();
//		tmpM.Clear();
//	}

#ifdef READEX_LEVEL_1
	READEX_REGION_STOP(REG_Cluster_CreateSa_SolveF0vG0);
#endif

#ifdef READEX_LEVEL_1
	READEX_REGION_START(REG_Cluster_CreateSa_SaReg, "Cluster--CreateSa-SaReg", SCOREP_USER_REGION_TYPE_COMMON);
#endif

	// *** Regularization of Sa from NULL PIVOTS
	if ( configuration.regularization == REGULARIZATION::NULL_PIVOTS ) {

		SparseMatrix Kernel_Sa;
		SparseMatrix Kernel_Sa2;
		ESINFO(PROGRESS3) << "Salfa - regularization from matrix";

		double tmp_double;
		eslocal tmp_int;
		SparseMatrix _tmpSparseMat; // Regularization matrix for Salfa

		if (SYMMETRIC_SYSTEM) {

			int regularization_type = 1;

			switch (regularization_type) { //(configuration.SAsolver) {

				case 0: { //ESPRESO_SASOLVER::CPU_DENSE: {
					// *** Get kernel from GGt - faster as GGt is very sparse matrix - works only for symmetric systems
					//  - kernel G0*G0t is used as regularization matrix for Salfa - needs to be scaled as values of G0*G0t are 10^10 smaller than values of Salfa
					SparseMatrix GGt;
					GGt.MatMat(G0,'N',G0t);
					GGt.RemoveLower();
					GGt.get_kernel_from_K(GGt, _tmpSparseMat,Kernel_Sa,tmp_double, tmp_int, 1000000, configuration.SC_SIZE);
					_tmpSparseMat.ConvertToCSR(1);
					// *** END - Get kernel from GGt - faster as GGt is very sparse matrix - works only for symmetric systems
					break;
				}

				case 1: { //ESPRESO_SASOLVER::CPU_SPARSE: {
					// *** Get kernel from Salfa - slower as Salfa is dense and takes more resources to find kernel
					Salfa.get_kernel_from_K(Salfa, _tmpSparseMat,Kernel_Sa,tmp_double, tmp_int, -1, configuration.SC_SIZE);
					// *** END - Get kernel from Salfa
					break;
				}

				case 2: { //ESPRESO_SASOLVER::ACC_DENSE: {
					//TODO: Musi se zkonzultovat s Alexem - nefunguje na 100%

					// *** Regularization of Salfa from kernels of G0*G0t matrix - it is designed for Dissection solver - it does not provide regularization matrix

					// *** Make regularization matrix from kernels - warning - output is dense matrix -- as it is added to the dense Salfa who gives a fuck

					SparseMatrix GGt;
					GGt.MatMat(G0,'N',G0t);
					GGt.RemoveLower();

//#if defined(SOLVER_DISSECTION)
//					SparseSolverCPU GGt_solver;
//					GGt_solver.ImportMatrix_wo_Copy(GGt);
//					GGt_solver.Factorization ("G0G0t matrix");
//					GGt_solver.GetKernel(Kernel_Sa); // TODO: Kplus.GetKernels(Kernel_Sa, Kernel_Sa2) - upravit na tuto funkci - v sym. pripade bude Kernel_Sa2 prazdna
//#else
					GGt.get_kernel_from_K(GGt, _tmpSparseMat,Kernel_Sa,tmp_double, tmp_int, -1, configuration.SC_SIZE);
//#endif

					SparseMatrix Kernel_Sat;
					Kernel_Sa.MatTranspose(Kernel_Sat);
					_tmpSparseMat.MatMat(Kernel_Sa, 'N',Kernel_Sat);
					_tmpSparseMat.RemoveLower();
					// *** END - Make regularization matrix from kernels
					break;
				}

				case 3: { //ESPRESO_SASOLVER::ACC_DENSE: {
					//TODO: Musi se zkonzultovat s Alexem - nefunguje na 100%

					// *** Regularization of Salfa from kernels of G0*G0t matrix - it is designed for Dissection solver - it does not provide regularization matrix

					// *** Make regularization matrix from kernels - warning - output is dense matrix -- as it is added to the dense Salfa who gives a fuck

//#if defined(SOLVER_DISSECTION)
//					SparseSolverCPU GGt_solver;
//					GGt_solver.ImportMatrix_wo_Copy(Salfa);
//					GGt_solver.Factorization ("Salfa matrix");
//					GGt_solver.GetKernel(Kernel_Sa); // TODO: Kplus.GetKernels(Kernel_Sa, Kernel_Sa2) - upravit na tuto funkci - v sym. pripade bude Kernel_Sa2 prazdna
//#else
					Salfa.get_kernel_from_K(Salfa, _tmpSparseMat,Kernel_Sa,tmp_double, tmp_int, -1, configuration.SC_SIZE);
//#endif

					SparseMatrix Kernel_Sat;
					Kernel_Sa.MatTranspose(Kernel_Sat);
					_tmpSparseMat.MatMat(Kernel_Sa, 'N',Kernel_Sat);
					_tmpSparseMat.RemoveLower();
					// *** END - Make regularization matrix from kernels
					break;
				}

				default:
					ESINFO(GLOBAL_ERROR) << "Unknown method for symmetric Salfa regularization";

			}


			// *** Actuall Regularization of Salfa using regularization matrix
			double ro = Salfa.GetMeanOfDiagonalOfSymmetricMatrix();
			Salfa.MatAddInPlace(_tmpSparseMat,'N', ro);
			// *** END - Actuall Regularization of Salfa using regularization matrix



//			if (1 == 1) {
//				// *** Get kernel from GGt - faster as GGt is very sparse matrix - works only for symmetric systems
//				//  - kernel G0*G0t is used as regularization matrix for Salfa - needs to be scaled as values of G0*G0t are 10^10 smaller than values of Salfa
//				SparseMatrix GGt;
//				GGt.MatMat(G0,'N',G0t);
//				GGt.RemoveLower();
//				GGt.get_kernel_from_K(GGt, _tmpSparseMat,Kernel_Sa,tmp_double, tmp_int, 1000000, configuration.SC_SIZE);
//				_tmpSparseMat.ConvertToCSR(1);
//				// *** END - Get kernel from GGt - faster as GGt is very sparse matrix - works only for symmetric systems
//			} else {
//				// *** Get kernel from Salfa - slower as Salfa is dense and takes more resources to find kernel
//				Salfa.get_kernel_from_K(Salfa, _tmpSparseMat,Kernel_Sa,tmp_double, tmp_int, -1, configuration.SC_SIZE);
//				// *** END - Get kernel from Salfa
//			}
//
//			// *** Regularization of Salfa from kernels of G0*G0t matrix - it is designed for Dissection solver - it does not provide regularization matrix
//			if (1 == 0) {
//				// *** Make regularization matrix from kernels - warning - output is dense matrix -- as it is added to the dense Salfa who gives a fuck
//				Kernel_Sa.Clear();
//				_tmpSparseMat.Clear();
//				SparseMatrix Kernel_Sat;
//				Kernel_Sa.MatTranspose(Kernel_Sat);
//				_tmpSparseMat.MatMat(Kernel_Sa, 'N',Kernel_Sat);
//				_tmpSparseMat.RemoveLower();
//				// *** END - Make regularization matrix from kernels
//			}
//
//			// *** Regularization of Salfa
//			double ro = Salfa.GetMeanOfDiagonalOfSymmetricMatrix();
//			Salfa.MatAddInPlace(_tmpSparseMat,'N', ro);
//			// *** END - Regularization of Salfa

		}

		if (!SYMMETRIC_SYSTEM) {


			int regularization_type = 1;

			switch (regularization_type) { //(configuration.SAsolver) {

				case 0: { //ESPRESO_SASOLVER::CPU_DENSE: {
					ESINFO(GLOBAL_ERROR) << "For non-symmetric systems G0*G0t cannot be used for regularization of Salfa - use direct regularization of Salfa";

					// TODO: Alex prozatim nevi jak spocist regularizacni matici pro nesym. system z GGt

					//			SparseMatrix GGt;
					//			GGt.MatMat(G0,'N',G0t);
					//			GGt.RemoveLower();
					//			GGt.get_kernels_from_nonsym_K(GGt, _tmpSparseMat, Kernel_Sa, Kernel_Sa2, tmp_double, tmp_int, -1, configuration.SC_SIZE);

					break;
				}

				case 1: { //ESPRESO_SASOLVER::CPU_SPARSE: {
					// *** Get kernel from Salfa - slower as Salfa is dense and takes more resources to find kernel
					Salfa.get_kernels_from_nonsym_K(Salfa, _tmpSparseMat, Kernel_Sa, Kernel_Sa2, tmp_double, tmp_int, -1, configuration.SC_SIZE);
					// *** END - Get kernel from Salfa
					break;
				}

				case 2: { //ESPRESO_SASOLVER::ACC_DENSE: {
					// *** Regularization of Salfa from kernels of G0*G0t matrix - it is designed for Dissection solver - it does not provide regularization matrix
					ESINFO(GLOBAL_ERROR) << "This type of Salfa regularization for nonsymmetric systems is not implemented yet";

					// *** Make regularization matrix from kernels - warning - output is dense matrix -- as it is added to the dense Salfa who gives a fuck

//					SparseMatrix GGt;
//					GGt.MatMat(G0,'N',G0t);
//					GGt.RemoveLower();
//
////#if defined(SOLVER_DISSECTION)
////					SparseSolverCPU GGt_solver;
////					GGt_solver.ImportMatrix_wo_Copy(GGt);
////					GGt_solver.Factorization ("G0G0t matrix");
////					GGt_solver.GetKernel(Kernel_Sa); // TODO: Kplus.GetKernels(Kernel_Sa, Kernel_Sa2) - upravit na tuto funkci - v sym. pripade bude Kernel_Sa2 prazdna
////#else
//					GGt.get_kernel_from_K(GGt, _tmpSparseMat,Kernel_Sa,tmp_double, tmp_int, -1, configuration.SC_SIZE);
////#endif
//
//					SparseMatrix Kernel_Sat;
//					Kernel_Sa.MatTranspose(Kernel_Sat);
//					_tmpSparseMat.MatMat(Kernel_Sa, 'N',Kernel_Sat);
//					_tmpSparseMat.RemoveLower();
					// *** END - Make regularization matrix from kernels

					break;
				}

				case 3: {//ESPRESO_SASOLVER::ACC_DENSE: {
					// *** Regularization of Salfa from kernels of G0*G0t matrix - it is designed for Dissection solver - it does not provide regularization matrix
					ESINFO(GLOBAL_ERROR) << "This type of Salfa regularization for nonsymmetric systems is not implemented yet";

//					// *** Regularization of Salfa from kernels of G0*G0t matrix - it is designed for Dissection solver - it does not provide regularization matrix
//
//					// *** Make regularization matrix from kernels - warning - output is dense matrix -- as it is added to the dense Salfa who gives a fuck
//
////#if defined(SOLVER_DISSECTION)
////					SparseSolverCPU GGt_solver;
////					GGt_solver.ImportMatrix_wo_Copy(Salfa);
////					GGt_solver.Factorization ("Salfa matrix");
////					GGt_solver.GetKernel(Kernel_Sa); // TODO: Kplus.GetKernels(Kernel_Sa, Kernel_Sa2) - upravit na tuto funkci - v sym. pripade bude Kernel_Sa2 prazdna
////#else
//					Salfa.get_kernel_from_K(Salfa, _tmpSparseMat,Kernel_Sa,tmp_double, tmp_int, -1, configuration.SC_SIZE);
////#endif
//
//					SparseMatrix Kernel_Sat;
//					Kernel_Sa.MatTranspose(Kernel_Sat);
//					_tmpSparseMat.MatMat(Kernel_Sa, 'N',Kernel_Sat);
//					_tmpSparseMat.RemoveLower();
//					// *** END - Make regularization matrix from kernels

					break;
				}

				default:
					ESINFO(GLOBAL_ERROR) << "Unknown method for symmetric Salfa regularization";

			}


		}
		// *** END - Regularization of Sa from NULL PIVOTS

		if (environment->print_matrices) {
			std::ofstream osSa(Logging::prepareFile("Kernel_Sa"));
			osSa << Kernel_Sa;
			osSa.close();
		}

		if (environment->print_matrices) {
			std::ofstream osSa(Logging::prepareFile("Kernel_Sa2"));
			osSa << Kernel_Sa2;
			osSa.close();
		}

		// *** Kernel (Kplus_R and Kplus_R2) correction for HTFETI method
		if (SYMMETRIC_SYSTEM) {
			for (size_t d = 0; d < domains.size(); d++) {
				SparseMatrix tR;

				SEQ_VECTOR < eslocal > rows_inds (Kernel_Sa.cols);

				//TODO: Tady to asi neni dobre pro singu/regu
				for (int i = 0; i < Kernel_Sa.cols; i++)
					rows_inds[i] = 1 + d * Kernel_Sa.cols + i;

				tR.CreateMatFromRowsFromMatrix_NewSize(Kernel_Sa,rows_inds);

				SparseMatrix TmpR;
				domains[d].Kplus_R.ConvertDenseToCSR(0);
				TmpR.MatMat( domains[d].Kplus_R, 'N', tR );

				if (TmpR.nnz == 0) {
					; //domains[d].Kplus_Rb = domains[d].Kplus_R;
				} else {
					domains[d].Kplus_Rb = TmpR;
					domains[d].Kplus_Rb.ConvertCSRToDense(0);
				}
			}
		}


		if (!SYMMETRIC_SYSTEM) {

			SparseMatrix LAMN_RHS;
			LAMN_RHS.MatMat(G02, 'T', Kernel_Sa2);

			SparseMatrix LAMN;

			eslocal set_bckp_F0 = F0.iparm[11];
			F0.iparm[11] = 2;
			F0.SolveMat_Sparse( LAMN_RHS, LAMN );
			F0.iparm[11] = set_bckp_F0;

			for (eslocal i = 0;i<LAMN.nnz;i++){
				LAMN.CSR_V_values[i] *= -1;
			}

			if (environment->print_matrices) {
				std::ofstream osSa(Logging::prepareFile("LAMN"));
				osSa << LAMN;
				osSa.close();
			}

			for (size_t d = 0; d < domains.size(); d++) {


				SparseMatrix LAMN_local;
				LAMN_local.CreateMatFromRowsFromMatrix_NewSize(LAMN, domains[d].B0_comp_map_vec);
				SparseMatrix B0t_LAMNlocal;
				B0t_LAMNlocal.MatMat(domains[d].B0_comp, 'T', LAMN_local);

				eslocal set_bckp = domains[d].Kplus.iparm[11];
				domains[d].Kplus.iparm[11] = 2;
				domains[d].Kplus.SolveMat_Sparse(B0t_LAMNlocal);
				domains[d].Kplus.iparm[11] = set_bckp;

				SparseMatrix Rb2;
				Rb2 = domains[d].Kplus_R2;
				Rb2.ConvertDenseToCSR(1);


				SparseMatrix tR; SparseMatrix tR2;

				SEQ_VECTOR < eslocal > rows_inds (Kernel_Sa.cols);
				for (int i = 0; i < Kernel_Sa.cols; i++)
					rows_inds[i] = 1 + d * Kernel_Sa.cols + i;

				tR. CreateMatFromRowsFromMatrix_NewSize(Kernel_Sa ,rows_inds);
				tR2.CreateMatFromRowsFromMatrix_NewSize(Kernel_Sa2,rows_inds);

				SparseMatrix TmpR;
				domains[d].Kplus_R.ConvertDenseToCSR(0);
				TmpR.MatMat( domains[d].Kplus_R, 'N', tR );

				if (TmpR.nnz == 0) {
					; //domains[d].Kplus_Rb = domains[d].Kplus_R;
				} else {
					domains[d].Kplus_Rb = TmpR;
					domains[d].Kplus_Rb.ConvertCSRToDense(0);
				}


				SparseMatrix TmpR2;
				domains[d].Kplus_R2.ConvertDenseToCSR(0);
				TmpR2.MatMat( domains[d].Kplus_R2, 'N', tR2 );
				TmpR2.MatAddInPlace(B0t_LAMNlocal,'N', 1.0);

				if (TmpR2.nnz == 0) {
					; //domains[d].Kplus_Rb = domains[d].Kplus_R;
				} else {
					domains[d].Kplus_Rb2 = TmpR2;
					domains[d].Kplus_Rb2.ConvertCSRToDense(0);
				}


				if (environment->print_matrices) {
					std::ofstream osR(Logging::prepareFile(d, "Rb_").c_str());
					SparseMatrix tmpR = domains[d].Kplus_Rb;
					tmpR.ConvertDenseToCSR(0);
					osR << tmpR;
					osR.close();

				}

				if (environment->print_matrices) {
					std::ofstream osR(Logging::prepareFile(d, "Rb2_").c_str());
					SparseMatrix tmpR = domains[d].Kplus_Rb2;
					tmpR.ConvertDenseToCSR(0);
					osR << tmpR;
					osR.close();

				}


			}

		}
		// *** Kernel (Kplu_R and Kplus_R2) correction for HTFETI method

	}

	// *** END - Regularization of Sa from NULL PIVOTS


	if ( configuration.regularization == REGULARIZATION::FIX_POINTS ) {

		// *** Regularization of Salfa from FIX points

		// TODO: GENERALIZE

		 TimeEvent reg_Sa_time("Salfa regularization "); if (Measure::report(CLUSTER)) { reg_Sa_time.start(); }

		SparseMatrix N, Nt, NNt;

		for (size_t i = 0; i < domains.size(); i++) {
			SparseMatrix Eye;
			Eye.CreateEye(domains[i].Kplus_R.cols);
			N.MatAppend(Eye);
			Nt.MatAppend(Eye);
		}

		Nt.MatTranspose();
		NNt.MatMat(N, 'N', Nt);
		NNt.RemoveLower();

		double ro = Salfa.GetMeanOfDiagonalOfSymmetricMatrix();
		ro = 0.5 * ro;

		Salfa.MatAddInPlace(NNt,'N', ro);
		Salfa.mtype = this->mtype;

		 if (Measure::report(CLUSTER)) { reg_Sa_time.end(); reg_Sa_time.printStatMPI(); Sa_timing.addEvent(reg_Sa_time); }

		// *** END - Regularization of Salfa from FIX points

	}
	// *** END - Sa regularization

#ifdef READEX_LEVEL_1
	READEX_REGION_STOP(REG_Cluster_CreateSa_SaReg);
#endif

#ifdef READEX_LEVEL_1
	READEX_REGION_START(REG_Cluster_CreateSa_SaFactorization, "Cluster--CreateSa-SaFactorization", SCOREP_USER_REGION_TYPE_COMMON);
#endif

	if (environment->print_matrices) {
		std::ofstream osSa(Logging::prepareFile("Salfa_reg"));
		osSa << Salfa;
		osSa.close();
	}

	// *** Salfa import into dense|sparse|acc solver and factorization
	switch (configuration.SAsolver) {

		// Default settings and fastest for CPU
		case ESPRESO_SASOLVER::CPU_DENSE: {
			TimeEvent factd_Sa_time("Salfa factorization - dense ");
			if (Measure::report(CLUSTER)) { factd_Sa_time.start(); }

			Salfa.ConvertCSRToDense(1);

			Sa_dense_cpu.ImportMatrix(Salfa);
			Sa_dense_cpu.Factorization("Salfa - dense ");

			if (Measure::report(CLUSTER)) { factd_Sa_time.end(); }//factd_Sa_time.printStatMPI();
			if (Measure::report(CLUSTER)) { Sa_timing.addEvent(factd_Sa_time); }
			break;
		}

		case ESPRESO_SASOLVER::CPU_SPARSE: {
			TimeEvent fact_Sa_time("Salfa factorization ");
			if (Measure::report(CLUSTER)) { fact_Sa_time.start(); }

			if (MPIrank == 0)  {
				Sa.msglvl = 1;
			}
			Sa.ImportMatrix(Salfa);
			Sa.Factorization("salfa");
			if (MPIrank == 0) {
				Sa.msglvl = 0;
			}
			if (Measure::report(CLUSTER)) { fact_Sa_time.end(); } //fact_Sa_time.printStatMPI(); }
			if (Measure::report(CLUSTER)) { Sa_timing.addEvent(fact_Sa_time); }
			break;
		}

		case ESPRESO_SASOLVER::ACC_DENSE: {
			TimeEvent factd_Sa_time("Salfa factorization - dense ");
			if (Measure::report(CLUSTER)) {factd_Sa_time.start();}

			//TODO: This works only for CuSolver
			//TODO: Radim Vavrik
			Salfa.type = 'G'; // for cuSolver only
			Salfa.ConvertCSRToDense(1);
			Salfa.type = 'S'; //for cuSolver only

			Sa_dense_acc.ImportMatrix(Salfa);
			Sa_dense_acc.Factorization("Salfa - dense ");

			if (Measure::report(CLUSTER)) { factd_Sa_time.end(); } //factd_Sa_time.printStatMPI();
			if (Measure::report(CLUSTER)) { Sa_timing.addEvent(factd_Sa_time); }
			break;
		}

		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented Salfa solver.";

	}
	// *** END - Salfa import into dense|sparse|acc solver and factorization


	if (Measure::report(CLUSTER)) { Sa_timing.totalTime.end(); }

	//TODO: Need fix - for detailed measurements
	//	Sa_timing.printStatsMPI();
	MKL_Set_Num_Threads(1);

	Sa.m_Kplus_size = Salfa.cols;
	vec_alfa.resize(Sa.m_Kplus_size);

#ifdef READEX_LEVEL_1
	READEX_REGION_STOP(REG_Cluster_CreateSa_SaFactorization);
#endif

#ifdef READEX_LEVEL_1
	READEX_REGION_STOP(REG_Cluster_CreateSa);
#endif

}

void ClusterBase::Create_G_perCluster() {

#ifdef READEX_LEVEL_1
	READEX_REGION_START(REG_Cluster_CreateG1_perCluster, "Cluster--CreateG1-perCluster", SCOREP_USER_REGION_TYPE_COMMON);
#endif

	SparseMatrix tmpM;

	TimeEvent G1_1_time ("Create G1 per clust t. : MatMat+MatTrans ");
	if (Measure::report(CLUSTER)) { G1_1_time.start(); }
	TimeEvent G1_1_mem  ("Create G1 per clust mem: MatMat+MatTrans ");
	if (Measure::report(CLUSTER)) { G1_1_mem.startWithoutBarrier(GetProcessMemory_u()); }

	int MPIrank;
	MPI_Comm_rank (environment->MPICommunicator, &MPIrank);

	PAR_VECTOR < SparseMatrix > tmp_Mat (domains.size());
	PAR_VECTOR < SparseMatrix > tmp_Mat2 (domains.size());

	#pragma omp parallel for
	for (size_t j = 0; j < domains.size(); j++) {

	 	if (domains[j].Kplus_R.nnz != 0) {

			SparseMatrix Rt;
			SparseMatrix Rt2;

			SparseMatrix B;
			B = domains[j].B1;

			switch (configuration.regularization) {
			case REGULARIZATION::FIX_POINTS:
				Rt = domains[j].Kplus_R;
				Rt.ConvertDenseToCSR(1);
				Rt.MatTranspose();

				if (!SYMMETRIC_SYSTEM) {
					Rt2 = domains[j].Kplus_R2;
					Rt2.ConvertDenseToCSR(1);
					Rt2.MatTranspose();
				}

				break;
			case REGULARIZATION::NULL_PIVOTS:
				Rt = domains[j].Kplus_Rb;
				Rt.ConvertDenseToCSR(1);
				Rt.MatTranspose();

				if (!SYMMETRIC_SYSTEM) {
					Rt2 = domains[j].Kplus_Rb2;
					Rt2.ConvertDenseToCSR(1);
					Rt2.MatTranspose();
				}

				break;
			default:
				ESINFO(GLOBAL_ERROR) << "Not implemented type of regularization.";
			}

			Rt.ConvertCSRToDense(1);
			//Create_G1_perSubdomain(Rt, domains[j].B1, tmp_Mat[j]);
			Create_G_perSubdomain(Rt, B, tmp_Mat[j]);

			if (!SYMMETRIC_SYSTEM) {
				Rt2.ConvertCSRToDense(1);
				//Create_G1_perSubdomain(Rt2, domains[j].B1, tmp_Mat2[j]);
				Create_G_perSubdomain(Rt2, B, tmp_Mat2[j]);
			}

	 	} // END: if

	} //end cilk for

	if (Measure::report(CLUSTER)) {
		G1_1_time.end();
		G1_1_time.printStatMPI();
		//G1_1_time.printLastStatMPIPerNode();
		G1_1_mem.endWithoutBarrier(GetProcessMemory_u());
		G1_1_mem.printStatMPI();
		//G1_1_mem.printLastStatMPIPerNode();
	}
	//G1_1_mem.printLastStatMPIPerNode();
	TimeEvent G1_2_time ("Create G1 per clust t. : Par.red.+MatAdd ");
	if (Measure::report(CLUSTER)) { G1_2_time.start(); }
	TimeEvent G1_2_mem  ("Create G1 per clust mem: Par.red.+MatAdd ");
	if (Measure::report(CLUSTER)) { G1_2_mem.startWithoutBarrier(GetProcessMemory_u()); }

	for (size_t j = 1; j < tmp_Mat.size(); j *= 2) {
		#pragma omp parallel for
		for (size_t i = 0; i <= tmp_Mat.size() / (2 * j); i++) {

			if (i * 2 * j + j < tmp_Mat.size()) {
				if (USE_HFETI == 1) {
					tmp_Mat[i * 2 * j].MatAddInPlace( tmp_Mat[i * 2 * j + j], 'N', 1.0 ); 	//Fixed for empty matrix in MatAddInPlace
				} else {
					tmp_Mat[i * 2 * j].MatAppend(tmp_Mat[i * 2 * j + j]); 					//Fixed for empty matrix in MatAddInPlace
				}
				tmp_Mat[i * 2 * j + j].Clear();

				if (!SYMMETRIC_SYSTEM) {
					if (USE_HFETI == 1) {
						tmp_Mat2[i * 2 * j].MatAddInPlace(tmp_Mat2[i * 2 * j + j], 'N', 1.0 );	//Fixed for empty matrix in MatAddInPlace
					} else {
						tmp_Mat2[i * 2 * j].MatAppend(tmp_Mat2[i * 2 * j + j]);					//Fixed for empty matrix in MatAddInPlace
					}
					tmp_Mat2[i * 2 * j + j].Clear();
				}

			}
		}
	}

	if (Measure::report(CLUSTER)) {
		G1_2_time.end();
		G1_2_time.printStatMPI();
		//G1_2_time.printLastStatMPIPerNode();
		G1_2_mem.endWithoutBarrier(GetProcessMemory_u());
		G1_2_mem.printStatMPI();
		//G1_2_mem.printLastStatMPIPerNode();
	}

	// Save resulting matrix G1
	G1.swap( tmp_Mat[0] );

	for (size_t i = 0; i < G1.CSR_V_values.size(); i++) {
		G1.CSR_V_values[i] = -1.0 * G1.CSR_V_values[i];
	}


	if (!SYMMETRIC_SYSTEM) {
		G2.swap( tmp_Mat2[0] );
		for (size_t i = 0; i < G2.CSR_V_values.size(); i++) {
			G2.CSR_V_values[i] = -1.0 * G2.CSR_V_values[i];
		}

	}

#ifdef READEX_LEVEL_1
	READEX_REGION_STOP(REG_Cluster_CreateG1_perCluster);
#endif

}

void ClusterBase::Create_G_perSubdomain (SparseMatrix &R_in, SparseMatrix &B_in, SparseMatrix &G_out) {

	if (B_in.nnz > 0) {

		SEQ_VECTOR < SEQ_VECTOR < double > > tmpG (B_in.nnz, SEQ_VECTOR <double> (R_in.rows,0));
		SEQ_VECTOR <eslocal > G_I_row_indices;

		G_I_row_indices.resize(B_in.nnz);

		eslocal indx = 0;
		for (eslocal r = 0; r < R_in.rows; r++)
			tmpG[indx][r] += B_in.V_values[0] * R_in.dense_values[R_in.rows * (B_in.J_col_indices[0]-1) + r];

		G_I_row_indices[indx] = B_in.I_row_indices[0];

		for (size_t i = 1; i < B_in.I_row_indices.size(); i++) {

			if (B_in.I_row_indices[i-1] != B_in.I_row_indices[i])
				indx++;

			for (eslocal r = 0; r < R_in.rows; r++)
				tmpG[indx][r] += B_in.V_values[i] * R_in.dense_values[R_in.rows * (B_in.J_col_indices[i]-1) + r];

			G_I_row_indices[indx] = B_in.I_row_indices[i];
		}

		G_I_row_indices.resize(indx+1);
		tmpG.resize(indx+1);

		SEQ_VECTOR <eslocal>    G_I; G_I.reserve( tmpG.size() * tmpG[0].size());
		SEQ_VECTOR <eslocal>    G_J; G_J.reserve( tmpG.size() * tmpG[0].size());
		SEQ_VECTOR <double> G_V; G_V.reserve( tmpG.size() * tmpG[0].size());

		for (size_t i = 0; i < tmpG.size(); i++) {
			for (size_t j = 0; j < tmpG[i].size(); j++){
				if (tmpG[i][j] != 0) {
					G_I.push_back(G_I_row_indices[i]);
					G_J.push_back(j+1);
					G_V.push_back(tmpG[i][j]);
				}
			}
		}

		G_out.I_row_indices = G_J;
		G_out.J_col_indices = G_I;
		G_out.V_values      = G_V;
		G_out.cols = B_in.rows;
		G_out.rows = R_in.rows;
		G_out.nnz  = G_I.size();
		G_out.type = 'G';

	} else {

		G_out.cols = B_in.rows;
		G_out.rows = R_in.rows;

		G_out.I_row_indices.resize(1,0);
		G_out.J_col_indices.resize(1,0);
		G_out.V_values.resize(0,0);

		G_out.nnz  = 0;
		G_out.type = 'G';
	}

	G_out.ConvertToCSRwithSort(1);

}

//void ClusterBase::Create_G1_perCluster() {
//
//	SparseMatrix tmpM;
//
//	if (USE_HFETI == 1) {
//		//G1 = G1 + trans(B1 * domains[d].Kplus_R) for all domains
//		//for (eslocal d = 0; d < domains.size(); d++)
//		//{
//		//	tmpM.MatMat( domains[d].B1, 'N', domains[d].Kplus_R);
//		//	G1.MatAddInPlace( tmpM, 'N', 1.0 );
//		//	tmpM.Clear();
//		//}
//		//G1.MatTranspose();
//
//		//// OK - but sequential
//		//for (eslocal d = 0; d < domains.size(); d++)					//HFETI
//		//{
//		//	tmpM.MatMat( domains[d].B1t, 'T', domains[d].Kplus_R);
//		//	G1.MatAddInPlace( tmpM, 'N', 1.0 );
//		//	tmpM.Clear();
//		//}
//		//G1.MatTranspose();
//
//
//
//
//		//eslocal threads = 24;
//		//eslocal n_domains = domains.size();
//		//vector < SparseMatrix > tmp_Mat (threads);
//		//cilk_for (eslocal t = 0; t < threads; t++ ) {
//		//	for (eslocal i = t*(n_domains/threads+1); i < (t+1)*(n_domains/threads+1); i++ ) {
//		//		if (i < n_domains) {
//		//			SparseMatrix tmpM_l;
//		//			//tmpM_l.MatMat( domains[i].Kplus_R, 'T', domains[i].B1t);
//		//			tmpM_l.MatMat( domains[i].B1t, 'T', domains[i].Kplus_R);
//		//			tmpM_l.MatTranspose();
//		//			tmp_Mat[t].MatAddInPlace( tmpM_l, 'N', 1.0 );
//		//		}
//		//	}
//		//}
//
//		TimeEvent G1_1_time ("Create G1 per clust t. : MatMat+MatTrans ");
//		G1_1_time.start();
//		TimeEvent G1_1_mem  ("Create G1 per clust mem: MatMat+MatTrans ");
//		G1_1_mem.startWithoutBarrier(GetProcessMemory_u());
//
//
//		//SparseMatrix Rt;
//		//SparseMatrix B;
//		//
//		//domains[0].Kplus_R.MatTranspose(Rt);
//		//B = domains[0].B1;
//		//Rt.ConvertCSRToDense(0);
//		//
//		//vector < vector < double > > tmpG (B.nnz, vector <double> (Rt.rows,0));
//		//vector <eslocal > G_I_row_indices;
//		//G_I_row_indices.resize(B.nnz);
//
//		//eslocal indx = 0;
//		//for (eslocal r = 0; r < Rt.rows; r++)
//		//	tmpG[indx][r] += B.V_values[0] * Rt.dense_values[Rt.rows * (B.J_col_indices[0]-1) + r];
//		//
//		//G_I_row_indices[indx] = B.I_row_indices[0];
//
//		//for (eslocal i = 1; i < B.I_row_indices.size(); i++) {
//		//
//		//	if (B.I_row_indices[i-1] != B.I_row_indices[i])
//		//		indx++;
//
//		//	for (eslocal r = 0; r < Rt.rows; r++)
//		//		tmpG[indx][r] += B.V_values[i] * Rt.dense_values[Rt.rows * (B.J_col_indices[i]-1) + r];
//		//
//		//	G_I_row_indices[indx] = B.I_row_indices[i];
//		//}
//		//
//		//G_I_row_indices.resize(indx+1);
//		//tmpG.resize(indx+1);
//
//		//vector <int>    G_I; G_I.reserve( tmpG.size() * tmpG[0].size());
//		//vector <int>    G_J; G_J.reserve( tmpG.size() * tmpG[0].size());
//		//vector <double> G_V; G_V.reserve( tmpG.size() * tmpG[0].size());
//
//		//for (eslocal i = 0; i < tmpG.size(); i++) {
//		//	for (eslocal j = 0; j < tmpG[i].size(); j++){
//		//		if (tmpG[i][j] != 0) {
//		//			G_I.push_back(G_I_row_indices[i]);
//		//			G_J.push_back(j+1);
//		//			G_V.push_back(tmpG[i][j]);
//		//		}
//		//	}
//		//}
//
//		//SparseMatrix Gcoo;
//		//Gcoo.I_row_indices = G_J; //G_I;
//		//Gcoo.J_col_indices = G_I; //G_J;
//		//Gcoo.V_values      = G_V;
//		//Gcoo.cols = B.rows; //R.cols;
//		//Gcoo.rows = Rt.rows; //B.rows;
//		//Gcoo.nnz  = G_I.size();
//		//Gcoo.type = 'G';
//		//
//		//Gcoo.ConvertToCSRwithSort(1);
//
//
//
//
//		//SparseMatrix Gtmp;
//		//SparseMatrix Gtmpt;
//		//
//		//Gtmp.MatMat( domains[0].B1t, 'T', domains[0].Kplus_R);
//		//Gtmp.MatTranspose(Gtmpt);
//		//
//		//Gtmp.ConvertToCOO(0);
//		//Gtmpt.ConvertToCOO(0);
//
//		int MPIrank;
//		MPI_Comm_rank (environment->MPICommunicator, &MPIrank);
//
//		PAR_VECTOR < SparseMatrix > tmp_Mat (domains.size());
//		#pragma omp parallel for
//for (size_t j = 0; j < tmp_Mat.size(); j++) {
//			// V1
//			//tmp_Mat[j].MatMat( domains[j].B1t, 'T', domains[j].Kplus_R);
//			//tmp_Mat[j].MatTranspose();
//
//			// V2 - not cool
//			//tmp_Mat[j].MatMatSorted( domains[j].Kplus_R, 'T', domains[j].B1t); // - pozor mooooc pomale a MKL rutina zere mooooooc pameti
//
//			// V3
//			SparseMatrix Gcoo;
//			if (domains[j].B1.nnz > 0) {
//
//				SparseMatrix Rt;
//				SparseMatrix B;
//
//				switch (configuration.regularization) {
//				case REGULARIZATION::FIX_POINTS:
//					Rt = domains[j].Kplus_R;
//					Rt.ConvertDenseToCSR(1);
//					Rt.MatTranspose();
//					//domains[j].Kplus_R.MatTranspose(Rt);
//					break;
//				case REGULARIZATION::NULL_PIVOTS:
//					Rt = domains[j].Kplus_Rb;
//					Rt.ConvertDenseToCSR(1);
//					Rt.MatTranspose();
//					//domains[j].Kplus_Rb.MatTranspose(Rt);
//					break;
//				default:
//					ESINFO(GLOBAL_ERROR) << "Not implemented regularization.";
//				}
//
//				//Rt = domains[j].Kplus_R;
//				//Rt.MatTranspose();
//
//				Rt.ConvertCSRToDense(1);
//				B = domains[j].B1;
//
//				SEQ_VECTOR < SEQ_VECTOR < double > > tmpG (B.nnz, SEQ_VECTOR <double> (Rt.rows,0));
//				SEQ_VECTOR <eslocal > G_I_row_indices;
//				G_I_row_indices.resize(B.nnz);
//
//				eslocal indx = 0;
//				for (eslocal r = 0; r < Rt.rows; r++)
//					tmpG[indx][r] += B.V_values[0] * Rt.dense_values[Rt.rows * (B.J_col_indices[0]-1) + r];
//
//				G_I_row_indices[indx] = B.I_row_indices[0];
//
//				for (size_t i = 1; i < B.I_row_indices.size(); i++) {
//
//					if (B.I_row_indices[i-1] != B.I_row_indices[i])
//						indx++;
//
//					for (eslocal r = 0; r < Rt.rows; r++)
//						tmpG[indx][r] += B.V_values[i] * Rt.dense_values[Rt.rows * (B.J_col_indices[i]-1) + r];
//
//					G_I_row_indices[indx] = B.I_row_indices[i];
//				}
//
//				G_I_row_indices.resize(indx+1);
//				tmpG.resize(indx+1);
//
//				SEQ_VECTOR <eslocal>    G_I; G_I.reserve( tmpG.size() * tmpG[0].size());
//				SEQ_VECTOR <eslocal>    G_J; G_J.reserve( tmpG.size() * tmpG[0].size());
//				SEQ_VECTOR <double> G_V; G_V.reserve( tmpG.size() * tmpG[0].size());
//
//				for (size_t i = 0; i < tmpG.size(); i++) {
//					for (size_t j = 0; j < tmpG[i].size(); j++){
//						if (tmpG[i][j] != 0) {
//							G_I.push_back(G_I_row_indices[i]);
//							G_J.push_back(j+1);
//							G_V.push_back(tmpG[i][j]);
//						}
//					}
//				}
//
//
//
//
//				Gcoo.I_row_indices = G_J; //G_I;
//				Gcoo.J_col_indices = G_I; //G_J;
//				Gcoo.V_values      = G_V;
//				Gcoo.cols = B.rows; //R.cols;
//				Gcoo.rows = Rt.rows; //B.rows;
//				Gcoo.nnz  = G_I.size();
//				Gcoo.type = 'G';
//
//				Gcoo.ConvertToCSRwithSort(1);
//
//			} else {
//				Gcoo.cols = domains[j].B1.rows;
//				if (configuration.regularization == REGULARIZATION::FIX_POINTS) {
//					Gcoo.rows = domains[j].Kplus_R.rows;
//				} else {
//					Gcoo.rows = domains[j].Kplus_Rb.rows;
//				}
//				Gcoo.I_row_indices.resize(1,0);
//				Gcoo.J_col_indices.resize(1,0);
//				Gcoo.V_values.resize(0,0);
//
//				Gcoo.nnz  = 0;
//				Gcoo.type = 'G';
//
//				Gcoo.ConvertToCSRwithSort(1);
//			}
//
//			tmp_Mat[j] = Gcoo;
//
//			// END - V3
//
//			//if (MPIrank == 0)
//			//	cout << j << endl;
//
//		}
//
//		G1_1_time.end();
//		G1_1_time.printStatMPI();
//		//G1_1_time.printLastStatMPIPerNode();
//		G1_1_mem.endWithoutBarrier(GetProcessMemory_u());
//		G1_1_mem.printStatMPI();
//		//G1_1_mem.printLastStatMPIPerNode();
//
//		TimeEvent G1_2_time ("Create G1 per clust t. : Par.red.+MatAdd ");
//		G1_2_time.start();
//		TimeEvent G1_2_mem  ("Create G1 per clust mem: Par.red.+MatAdd ");
//		G1_2_mem.startWithoutBarrier(GetProcessMemory_u());
//
//		for (size_t j = 1; j < tmp_Mat.size(); j *= 2) {
//			#pragma omp parallel for
//for (size_t i = 0; i <= tmp_Mat.size() / (2 * j); i++) {
//				if (i * 2 * j + j < tmp_Mat.size()) {
//					tmp_Mat[i * 2 * j].MatAddInPlace(tmp_Mat[i * 2 * j + j], 'N', 1.0 ); //  MFETI - MatAppend(tmp_Mat[i + j]);
//					tmp_Mat[i * 2 * j + j].Clear();
//				}
//			}
//		}
//
//		G1_2_time.end();
//		G1_2_time.printStatMPI();
//		//G1_2_time.printLastStatMPIPerNode();
//
//		G1_2_mem.endWithoutBarrier(GetProcessMemory_u());
//		G1_2_mem.printStatMPI();
//		//G1_2_mem.printLastStatMPIPerNode();
//
//		G1 = tmp_Mat[0];
//		tmp_Mat[0].Clear();
//		//G1.MatTranspose();
//
//
//	} else {
//
//		PAR_VECTOR < SparseMatrix > tmp_Mat (domains.size());
//		#pragma omp parallel for
//for (size_t j = 0; j < domains.size(); j++) {
//
//			// V1
//			//tmp_Mat[j].MatMat( domains[j].B1t, 'T', domains[j].Kplus_R);
//			//tmp_Mat[j].MatTranspose();
//
//			// V2
//			//tmp_Mat[j].MatMat( domains[j].Kplus_R, 'T', domains[j].B1t);
//
//			// V3
//			SparseMatrix Rt;
//			SparseMatrix B;
//
//			switch (configuration.regularization) {
//			case REGULARIZATION::FIX_POINTS:
//				Rt = domains[j].Kplus_R;
//				Rt.ConvertDenseToCSR(1);
//				Rt.MatTranspose();
//				//domains[j].Kplus_R.MatTranspose(Rt);
//				break;
//			case REGULARIZATION::NULL_PIVOTS:
//				Rt = domains[j].Kplus_Rb;
//				Rt.ConvertDenseToCSR(1);
//				Rt.MatTranspose();
//				//domains[j].Kplus_Rb.MatTranspose(Rt);
//				break;
//			default:
//				ESINFO(GLOBAL_ERROR) << "Not implemented regularization.";
//			}
//
//			Rt.ConvertCSRToDense(1);
//			B = domains[j].B1;
//
//			SEQ_VECTOR < SEQ_VECTOR < double > > tmpG (B.nnz, SEQ_VECTOR <double> (Rt.rows,0));
//			SEQ_VECTOR <eslocal > G_I_row_indices;
//			G_I_row_indices.resize(B.nnz);
//
//			eslocal indx = 0;
//			for (eslocal r = 0; r < Rt.rows; r++)
//				tmpG[indx][r] += B.V_values[0] * Rt.dense_values[Rt.rows * (B.J_col_indices[0]-1) + r];
//
//			G_I_row_indices[indx] = B.I_row_indices[0];
//
//			for (size_t i = 1; i < B.I_row_indices.size(); i++) {
//
//				if (B.I_row_indices[i-1] != B.I_row_indices[i])
//					indx++;
//
//				for (eslocal r = 0; r < Rt.rows; r++)
//					tmpG[indx][r] += B.V_values[i] * Rt.dense_values[Rt.rows * (B.J_col_indices[i]-1) + r];
//
//				G_I_row_indices[indx] = B.I_row_indices[i];
//			}
//
//			G_I_row_indices.resize(indx+1);
//			tmpG.resize(indx+1);
//
//			SEQ_VECTOR <eslocal>    G_I; G_I.reserve( tmpG.size() * tmpG[0].size());
//			SEQ_VECTOR <eslocal>    G_J; G_J.reserve( tmpG.size() * tmpG[0].size());
//			SEQ_VECTOR <double> G_V; G_V.reserve( tmpG.size() * tmpG[0].size());
//
//			for (size_t i = 0; i < tmpG.size(); i++) {
//				for (size_t j = 0; j < tmpG[i].size(); j++){
//					if (tmpG[i][j] != 0) {
//						G_I.push_back(G_I_row_indices[i]);
//						G_J.push_back(j+1);
//						G_V.push_back(tmpG[i][j]);
//					}
//				}
//			}
//
//			SparseMatrix Gcoo;
//			Gcoo.I_row_indices = G_J; //G_I;
//			Gcoo.J_col_indices = G_I; //G_J;
//			Gcoo.V_values      = G_V;
//			Gcoo.cols = B.rows; //R.cols;
//			Gcoo.rows = Rt.rows; //B.rows;
//			Gcoo.nnz  = G_I.size();
//			Gcoo.type = 'G';
//
//			Gcoo.ConvertToCSRwithSort(1);
//
//			tmp_Mat[j] = Gcoo;
//
//		}
//
//		for (size_t j = 1; j < tmp_Mat.size(); j *= 2) {
//			#pragma omp parallel for
//			for (size_t i = 0; i <= tmp_Mat.size() / (2 * j); i++) {
//				if (i * 2 * j + j < tmp_Mat.size()) {
//					tmp_Mat[i * 2 * j].MatAppend(tmp_Mat[i * 2 * j + j]);
//					tmp_Mat[i * 2 * j + j].Clear();
//				}
//			}
//		}
//		G1.MatAppend(tmp_Mat[0]);
//		//for (eslocal d = 0; d < domains.size(); d++)					//MFETI
//		//{
//		//	//tmpM.MatMat( domains[d].B1, 'N', domains[d].Kplus_R);
//		//	tmpM.MatMat( domains[d].B1t, 'T', domains[d].Kplus_R);
//		//	tmpM.MatTranspose();
//		//	G1.MatAppend(tmpM);
//		//	tmpM.Clear();
//		//}
//	}
//
//	// for both MFETI and HFETI
//	for (size_t i = 0; i < G1.CSR_V_values.size(); i++) {
//		G1.CSR_V_values[i] = -1.0 * G1.CSR_V_values[i];
//	}
//
//
////	std::stringstream ss;
////	ss << "G" << ".txt";
////	std::ofstream os(ss.str().c_str());
////	os << G1;
////	os.close();
//
//}

void ClusterBase::Compress_G1() {

	G1_comp.Clear();
	Compress_G(G1, G1_comp);
	G1.Clear();

	if (!SYMMETRIC_SYSTEM) {
		G2_comp.Clear();
		Compress_G(G2, G2_comp);
		G2.Clear();
	}

}

void ClusterBase::Compress_G( SparseMatrix &G_in, SparseMatrix &G_comp_out ) {

	G_in.ConvertToCOO( 1 );

	#pragma omp parallel for
	for (size_t j = 0; j < G_in.J_col_indices.size(); j++ ) {
		G_in.J_col_indices[j] = _my_lamdas_map_indices[ G_in.J_col_indices[j] -1 ] + 1;  // numbering from 1 in matrix
	}

	G_in.cols = my_lamdas_indices.size();
	G_in.ConvertToCSRwithSort( 1 );

	G_comp_out.CSR_I_row_indices.swap( G_in.CSR_I_row_indices );
	G_comp_out.CSR_J_col_indices.swap( G_in.CSR_J_col_indices );
	G_comp_out.CSR_V_values     .swap( G_in.CSR_V_values		 );

	G_comp_out.rows = G_in.rows;
	G_comp_out.cols = G_in.cols;
	G_comp_out.nnz  = G_in.nnz;
	G_comp_out.type = G_in.type;

}

void ClusterBase::CreateVec_d_perCluster( SEQ_VECTOR<SEQ_VECTOR <double> > & f ) {

	SEQ_VECTOR <eslocal> kerindices (domains.size() + 1, 0);
	kerindices[0] = 0;

	for (size_t k = 1; k < kerindices.size(); k++) {
		kerindices[k] = kerindices[k-1] + domains[k-1].Kplus_R.cols;
	}

	eslocal size_d = 0;
	if ( USE_HFETI == 1 ) {
		for (size_t d = 0; d < domains.size(); d++) {
			if (size_d < domains[d].Kplus_R.cols)
				size_d = domains[d].Kplus_R.cols;
		}
	} else {
		for (size_t d = 0; d < domains.size(); d++) {
				size_d += domains[d].Kplus_R.cols;
		}
	}
	vec_d.resize( size_d );


	if (SYMMETRIC_SYSTEM) {

		if ( USE_HFETI == 1) {
			for (size_t d = 0; d < domains.size(); d++) {
				if ( configuration.regularization == REGULARIZATION::FIX_POINTS ) {
					domains[d].Kplus_R .DenseMatVec(f[domains[d].domain_global_index], vec_d, 'T', 0, 0         , 1.0 );
				} else {
					domains[d].Kplus_Rb.DenseMatVec(f[domains[d].domain_global_index], vec_d, 'T', 0, 0         , 1.0 );
				}
			}
		} else {
			for (size_t d = 0; d < domains.size(); d++) {											// MFETI
				if ( configuration.regularization == REGULARIZATION::FIX_POINTS ) {
					domains[d].Kplus_R. DenseMatVec(f[domains[d].domain_global_index], vec_d, 'T', 0, kerindices[d], 0.0 );				// MFETI
				} else {
					domains[d].Kplus_Rb.DenseMatVec(f[domains[d].domain_global_index], vec_d, 'T', 0, kerindices[d], 0.0 );				// MFETI
				}
			}
		}

	} else {
		// NON-SYMMETRIC SYSTEM
		if ( USE_HFETI == 1) {
			for (size_t d = 0; d < domains.size(); d++) {
				if ( configuration.regularization == REGULARIZATION::FIX_POINTS ) {
					domains[d].Kplus_R2 .DenseMatVec(f[domains[d].domain_global_index], vec_d, 'T', 0, 0         , 1.0 );
				} else {
					domains[d].Kplus_Rb2.DenseMatVec(f[domains[d].domain_global_index], vec_d, 'T', 0, 0         , 1.0 );
				}
			}
		} else {
			for (size_t d = 0; d < domains.size(); d++) {											// MFETI
				if ( configuration.regularization == REGULARIZATION::FIX_POINTS ) {
					domains[d].Kplus_R2. DenseMatVec(f[domains[d].domain_global_index], vec_d, 'T', 0, kerindices[d], 0.0 );				// MFETI
				} else {
					domains[d].Kplus_Rb2.DenseMatVec(f[domains[d].domain_global_index], vec_d, 'T', 0, kerindices[d], 0.0 );				// MFETI
				}
			}
		}

	}

	for (size_t i = 0; i < vec_d.size(); i++) {
		vec_d[i] = (-1.0) *  vec_d[i];
	}

}


void ClusterBase::CreateVec_b_perCluster( SEQ_VECTOR<SEQ_VECTOR <double> > & f )  {

	SEQ_VECTOR<SEQ_VECTOR<double> > x_prim_cluster ( domains.size() );

	#pragma omp parallel for
	for (size_t d = 0; d < domains.size(); d++) {
		x_prim_cluster[d] = f[domains[d].domain_global_index];
	}

	if (USE_HFETI == 0) {
		#pragma omp parallel for
		for (size_t d = 0; d < domains.size(); d++) {				// MFETI
				domains[d].multKplusLocal( x_prim_cluster[d] );
		}
	} else {
		multKplusGlobal_l( x_prim_cluster );						// HFETI
	}

//	// pro ukoncovani v primaru - vypocet up0
//	cilk_for (eslocal d = 0; d < domains.size(); d++) {
//		domains[d].up0 = x_prim_cluster[d];
//	}

//	vec_b_compressed.clear();
//	vec_b_compressed.resize(my_lamdas_indices.size(), 0.0);

	SEQ_VECTOR < double > y_out_tmp (domains[0].B1_comp_dom.rows);

	for (size_t d = 0; d < domains.size(); d++) {
		y_out_tmp.resize( domains[d].B1_comp_dom.rows );
		domains[d].B1_comp_dom.MatVec (x_prim_cluster[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

		for (size_t i = 0; i < domains[d].lambda_map_sub_local.size(); i++)
			vec_b_compressed[ domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i] - domains[d].vec_c[i];

	}

}


void ClusterBase::CreateVec_c_perCluster( SEQ_VECTOR <double> & vec_c_out )  {

	//vec_c_out.resize(my_lamdas_indices.size(), 0.0);

	for (size_t d = 0; d < domains.size(); d++) {
		for (size_t i = 0; i < domains[d].lambda_map_sub_local.size(); i++) {
			vec_c_out[ domains[d].lambda_map_sub_local[i] ] = domains[d].vec_c[i];
		}
	}


}

void ClusterBase::CreateVec_lb_perCluster( SEQ_VECTOR <double> & vec_lb_out )  {

	//vec_lb_out.resize(my_lamdas_indices.size(), -std::numeric_limits<double>::infinity());

	for (size_t d = 0; d < domains.size(); d++) {
		for (size_t i = 0; i < domains[d].lambda_map_sub_local.size(); i++){
			vec_lb_out[ domains[d].lambda_map_sub_local[i] ] = domains[d].vec_lb[i];
		}
	}
}


void ClusterBase::compress_lambda_vector  ( SEQ_VECTOR <double> & decompressed_vec_lambda )
{
	//compress vector for CG in main loop
	for (size_t i = 0; i < my_lamdas_indices.size(); i++)
		decompressed_vec_lambda[i] = decompressed_vec_lambda[my_lamdas_indices[i]];

	decompressed_vec_lambda.resize(my_lamdas_indices.size());
}

void ClusterBase::decompress_lambda_vector( SEQ_VECTOR <double> &   compressed_vec_lambda )
{
	SEQ_VECTOR <double> decompressed_vec_lambda (domains[0].B1.rows,0);

	for (size_t i = 0; i < my_lamdas_indices.size(); i++)
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

void ClusterBase::CreateDirichletPrec(Instance *instance)
{
	#pragma omp parallel for
	//for (size_t d = 0; d < instance->K.size(); d++) {
	for (size_t d = 0; d < domains.size(); d++) {

		SEQ_VECTOR<eslocal> perm_vec = domains[d].B1t_Dir_perm_vec;
		SEQ_VECTOR<eslocal> perm_vec_full(instance->K[domains[d].domain_global_index].rows);// (instance->K[d].rows);
		SEQ_VECTOR<eslocal> perm_vec_diff(instance->K[domains[d].domain_global_index].rows);// (instance->K[d].rows);

		SEQ_VECTOR<eslocal> I_row_indices_p(instance->K[domains[d].domain_global_index].nnz);// (instance->K[d].nnz);
		SEQ_VECTOR<eslocal> J_col_indices_p(instance->K[domains[d].domain_global_index].nnz);// (instance->K[d].nnz);

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
		K_modif.ConvertToCSR(1);
#endif
		SparseMatrix RegMatCRS = instance->RegMat[domains[d].domain_global_index]; //[d];
		RegMatCRS.ConvertToCSRwithSort(0);
		K_modif.MatAddInPlace(RegMatCRS, 'N', -1);
		// K_modif.RemoveLower();

		SEQ_VECTOR<SEQ_VECTOR<eslocal >> vec_I1_i2(K_modif.rows, SEQ_VECTOR<eslocal >(2, 1));
		eslocal offset = K_modif.CSR_I_row_indices[0] ? 1 : 0;

		for (eslocal i = 0; i < K_modif.rows; i++) {
			vec_I1_i2[i][0] = perm_vec_full[i];
			vec_I1_i2[i][1] = i; // position to create reverse permutation
		}

		std::sort(vec_I1_i2.begin(), vec_I1_i2.end(), [](const SEQ_VECTOR <eslocal >& a, const SEQ_VECTOR<eslocal>& b) {return a[0] < b[0];});

		// permutations made on matrix in COO format
		K_modif.ConvertToCOO(0);
		eslocal I_index, J_index;
		bool unsymmetric = !SYMMETRIC_SYSTEM;
		for (eslocal i = 0; i < K_modif.nnz; i++) {
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
		for (eslocal i = 0; i < K_modif.nnz; i++) {
			K_modif.I_row_indices[i] = I_row_indices_p[i];
			K_modif.J_col_indices[i] = J_col_indices_p[i];
		}
		K_modif.ConvertToCSRwithSort(1);
		{
			if (environment->print_matrices) {
				std::ofstream osS(Logging::prepareFile(d, "K_modif"));
				osS << K_modif;
				osS.close();
			}
		}

		// ------------------------------------------------------------------------------------------------------------------
		bool diagonalized_K_rr = configuration.preconditioner == ESPRESO_PRECONDITIONER::SUPER_DIRICHLET;
		//        PRECONDITIONER==NONE              - 0
		//        PRECONDITIONER==LUMPED            - 1
		//        PRECONDITIONER==WEIGHT_FUNCTION   - 2
		//        PRECONDITIONER==DIRICHLET         - 3
		//        PRECONDITIONER==SUPER_DIRICHLET   - 4
		//
		//        When next line is uncomment, var. PRECONDITIONER==DIRICHLET and PRECONDITIONER==SUPER_DIRICHLET provide identical preconditioner.
		//        bool diagonalized_K_rr = false
		// ------------------------------------------------------------------------------------------------------------------

		eslocal sc_size = perm_vec.size();

		if (sc_size == instance->K[domains[d].domain_global_index].rows) {
			domains[d].Prec = instance->K[domains[d].domain_global_index];
			domains[d].Prec.ConvertCSRToDense(1);
			// if physics.K[d] does not contain inner DOF
		} else {

			if (configuration.preconditioner == ESPRESO_PRECONDITIONER::DIRICHLET) {
				SparseSolverCPU createSchur;
//          createSchur.msglvl=1;
				eslocal sc_size = perm_vec.size();
				createSchur.ImportMatrix_wo_Copy(K_modif);
				createSchur.Create_SC(domains[d].Prec, sc_size, false);
				domains[d].Prec.ConvertCSRToDense(1);
			} else {
				SparseMatrix K_rr;
				SparseMatrix K_rs;
				SparseMatrix K_sr;
				SparseMatrix KsrInvKrrKrs;

				eslocal i_start = 0;
				eslocal nonsing_size = K_modif.rows - sc_size - i_start;
				eslocal j_start = nonsing_size;

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
					for (eslocal i = 0; i < K_rs.rows; i++) {
						for (eslocal j = K_rs.CSR_I_row_indices[i]; j < K_rs.CSR_I_row_indices[i + 1]; j++) {
							K_rs.CSR_V_values[j - offset] /= diagonals[i];
						}
					}
				} else {
					K_rr.getSubDiagBlockmatrix(K_modif, K_rr, i_start, nonsing_size);
					K_rr_solver.ImportMatrix_wo_Copy(K_rr);
//            K_rr_solver.msglvl = 1;
					K_rr_solver.SolveMat_Dense(K_rs);
				}
			}
		}
	}

	// implemented in clustercpu.cpp

//	#pragma omp parallel for
//	for (size_t d = 0; d < instance->K.size(); d++) {
//		SEQ_VECTOR<eslocal> perm_vec = domains[d].B1t_Dir_perm_vec;
//		SEQ_VECTOR<eslocal> perm_vec_full(instance->K[d].rows);
//		SEQ_VECTOR<eslocal> perm_vec_diff(instance->K[d].rows);
//
//		SEQ_VECTOR<eslocal> I_row_indices_p(instance->K[d].nnz);
//		SEQ_VECTOR<eslocal> J_col_indices_p(instance->K[d].nnz);
//
//		for (size_t i = 0; i < perm_vec.size(); i++) {
//			perm_vec[i] = perm_vec[i] - 1;
//		}
//
//		for (size_t i = 0; i < perm_vec_full.size(); i++) {
//			perm_vec_full[i] = i;
//		}
//
//		auto it = std::set_difference(perm_vec_full.begin(), perm_vec_full.end(), perm_vec.begin(), perm_vec.end(), perm_vec_diff.begin());
//		perm_vec_diff.resize(it - perm_vec_diff.begin());
//
//		perm_vec_full = perm_vec_diff;
//		perm_vec_full.insert(perm_vec_full.end(), perm_vec.begin(), perm_vec.end());
//
//		SparseMatrix K_modif = instance->K[d];
//		SparseMatrix RegMatCRS = instance->RegMat[d];
//		RegMatCRS.ConvertToCSRwithSort(0);
//		K_modif.MatAddInPlace(RegMatCRS, 'N', -1);
//		// K_modif.RemoveLower();
//
//		SEQ_VECTOR<SEQ_VECTOR<eslocal >> vec_I1_i2(K_modif.rows, SEQ_VECTOR<eslocal >(2, 1));
//		eslocal offset = K_modif.CSR_I_row_indices[0] ? 1 : 0;
//
//		for (eslocal i = 0; i < K_modif.rows; i++) {
//			vec_I1_i2[i][0] = perm_vec_full[i];
//			vec_I1_i2[i][1] = i; // position to create reverse permutation
//		}
//
//		std::sort(vec_I1_i2.begin(), vec_I1_i2.end(), [](const SEQ_VECTOR <eslocal >& a, const SEQ_VECTOR<eslocal>& b) {return a[0] < b[0];});
//
//		// permutations made on matrix in COO format
//		K_modif.ConvertToCOO(0);
//		eslocal I_index, J_index;
//		bool unsymmetric = !SYMMETRIC_SYSTEM;
//		for (eslocal i = 0; i < K_modif.nnz; i++) {
//			I_index = vec_I1_i2[K_modif.I_row_indices[i] - offset][1] + offset;
//			J_index = vec_I1_i2[K_modif.J_col_indices[i] - offset][1] + offset;
//			if (unsymmetric || I_index <= J_index) {
//				I_row_indices_p[i] = I_index;
//				J_col_indices_p[i] = J_index;
//			} else {
//				I_row_indices_p[i] = J_index;
//				J_col_indices_p[i] = I_index;
//			}
//		}
//		for (eslocal i = 0; i < K_modif.nnz; i++) {
//			K_modif.I_row_indices[i] = I_row_indices_p[i];
//			K_modif.J_col_indices[i] = J_col_indices_p[i];
//		}
//		K_modif.ConvertToCSRwithSort(1);
//		{
//			if (environment->print_matrices) {
//				std::ofstream osS(Logging::prepareFile(d, "K_modif"));
//				osS << K_modif;
//				osS.close();
//			}
//		}
//
//		// ------------------------------------------------------------------------------------------------------------------
//		bool diagonalized_K_rr = configuration.preconditioner == ESPRESO_PRECONDITIONER::SUPER_DIRICHLET;
//		//        PRECONDITIONER==NONE              - 0
//		//        PRECONDITIONER==LUMPED            - 1
//		//        PRECONDITIONER==WEIGHT_FUNCTION   - 2
//		//        PRECONDITIONER==DIRICHLET         - 3
//		//        PRECONDITIONER==SUPER_DIRICHLET   - 4
//		//
//		//        When next line is uncomment, var. PRECONDITIONER==DIRICHLET and PRECONDITIONER==SUPER_DIRICHLET provide identical preconditioner.
//		//        bool diagonalized_K_rr = false
//		// ------------------------------------------------------------------------------------------------------------------
//
//		eslocal sc_size = perm_vec.size();
//
//		if (sc_size == instance->K[d].rows) {
//			domains[d].Prec = instance->K[d];
//			domains[d].Prec.ConvertCSRToDense(1);
//			// if physics.K[d] does not contain inner DOF
//		} else {
//
//			if (configuration.preconditioner == ESPRESO_PRECONDITIONER::DIRICHLET) {
//				SparseSolverCPU createSchur;
////          createSchur.msglvl=1;
//				eslocal sc_size = perm_vec.size();
//				createSchur.ImportMatrix_wo_Copy(K_modif);
//				createSchur.Create_SC(domains[d].Prec, sc_size, false);
//				domains[d].Prec.ConvertCSRToDense(1);
//			} else {
//				SparseMatrix K_rr;
//				SparseMatrix K_rs;
//				SparseMatrix K_sr;
//				SparseMatrix KsrInvKrrKrs;
//
//				eslocal i_start = 0;
//				eslocal nonsing_size = K_modif.rows - sc_size - i_start;
//				eslocal j_start = nonsing_size;
//
//				K_rs.getSubBlockmatrix_rs(K_modif, K_rs, i_start, nonsing_size, j_start, sc_size);
//
//				if (SYMMETRIC_SYSTEM) {
//					K_rs.MatTranspose(K_sr);
//				} else {
//					K_sr.getSubBlockmatrix_rs(K_modif, K_sr, j_start, sc_size, i_start, nonsing_size);
//				}
//
//				domains[d].Prec.getSubDiagBlockmatrix(K_modif, domains[d].Prec, nonsing_size, sc_size);
//				SEQ_VECTOR<double> diagonals;
//				SparseSolverCPU K_rr_solver;
//
//				// K_rs is replaced by:
//				// a) K_rs = 1/diag(K_rr) * K_rs          (simplified Dirichlet precond.)
//				// b) K_rs =    inv(K_rr) * K_rs          (classical Dirichlet precond. assembled by own - not via PardisoSC routine)
//				if (diagonalized_K_rr) {
//					diagonals = K_modif.getDiagonal();
//					// diagonals is obtained directly from K_modif (not from K_rr to avoid assembling) thanks to its structure
//					//      K_modif = [K_rr, K_rs]
//					//                [K_sr, K_ss]
//					//
//					for (eslocal i = 0; i < K_rs.rows; i++) {
//						for (eslocal j = K_rs.CSR_I_row_indices[i]; j < K_rs.CSR_I_row_indices[i + 1]; j++) {
//							K_rs.CSR_V_values[j - offset] /= diagonals[i];
//						}
//					}
//				} else {
//					K_rr.getSubDiagBlockmatrix(K_modif, K_rr, i_start, nonsing_size);
//					K_rr_solver.ImportMatrix_wo_Copy(K_rr);
////            K_rr_solver.msglvl = 1;
//					K_rr_solver.SolveMat_Dense(K_rs);
//				}
//
//				KsrInvKrrKrs.MatMat(K_sr, 'N', K_rs);
//				domains[d].Prec.MatAddInPlace(KsrInvKrrKrs, 'N', -1);
////          if (!diagonalized_K_rr){
////				    cluster.domains[d].Prec.ConvertCSRToDense(1);
////          }
//			}
//
//		}
//
//		if (environment->print_matrices) {
//			std::ofstream osS(Logging::prepareFile(d, "S"));
//			SparseMatrix SC = domains[d].Prec;
//			if (configuration.preconditioner == ESPRESO_PRECONDITIONER::DIRICHLET) {
//				SC.ConvertDenseToCSR(1);
//			}
//			osS << SC;
//			osS.close();
//		}
//
//		ESINFO(PROGRESS3) << Info::plain() << ".";
//	}
}

// **** END - CLUSTER CLASS ************************************************
// *******************************************************************







