/*
 * supercluster.h
 *
 *  Created on: Jun 8, 2017
 *      Author: lriha
 */

#ifndef SRC_SOLVER_SPECIFIC_SUPERCLUSTER_H_
#define SRC_SOLVER_SPECIFIC_SUPERCLUSTER_H_

#include <omp.h>
#include "mpi.h"
// #include "mkl.h"

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <iomanip>
#include <map>

using std::vector;
using std::map;
using std::make_pair;

#include "../generic/SparseMatrix.h"
#include "sparsesolvers.h"
#include "clusters.h"
#include "../generic/utils.h"

namespace espreso {

class SuperClusterBase
{
public:
	SuperClusterBase(const ESPRESOSolver &configuration, Instance *instance_in):
    	configuration(configuration),
		instance(instance_in)
	{
		if (instance_in != NULL) {
			init();
		}
	}

	~SuperClusterBase() {
		clusters.clear();
		domains.clear();
		x_prim_cluster1.clear();
		x_prim_cluster2.clear();
		x_prim_cluster3.clear();
	}

	const ESPRESOSolver 	&configuration;
	Instance 				*instance;
	SEQ_VECTOR <Cluster> 	clusters;

	SEQ_VECTOR <Domain*> 	domains;			// vector of domains in the original numbering - for legacy purposes - only as pointers

	eslocal numClusters;
	eslocal number_of_subdomains_per_supercluster;

	bool USE_HFETI;
	bool USE_KINV;
	bool SYMMETRIC_SYSTEM;
	bool USE_FLOAT;

	eslocal PAR_NUM_THREADS;
	eslocal SOLVER_NUM_THREADS;
	MatrixType mtype;
	int MPIrank;
	eslocal dual_size;

	eslocal min_numClusters_per_MPI;

	void init() {
		numClusters 							= 1 + *std::max_element(instance->clustersMap.begin(), instance->clustersMap.end());
		number_of_subdomains_per_supercluster 	= instance->K.size();

		x_prim_cluster1.resize( number_of_subdomains_per_supercluster );
		x_prim_cluster2.resize( number_of_subdomains_per_supercluster );
		x_prim_cluster3.resize( number_of_subdomains_per_supercluster );

		domains.resize(number_of_subdomains_per_supercluster);
		for (eslocal c = 0; c < numClusters; c++) {
			clusters.push_back( Cluster(configuration, instance));
		}

		switch (configuration.method) {
		case ESPRESO_METHOD::TOTAL_FETI:
			USE_HFETI = false;
			break;
		case ESPRESO_METHOD::HYBRID_FETI:
			USE_HFETI = true;
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Unsupported FETI METHOD";
		}

		USE_KINV 			= configuration.use_schur_complement ? 1 : 0;
		PAR_NUM_THREADS 	= environment->PAR_NUM_THREADS;
		SOLVER_NUM_THREADS  = environment->SOLVER_NUM_THREADS;

		for (eslocal c = 0; c < numClusters; c++) {
			clusters[c].USE_HFETI 			= USE_HFETI;
			clusters[c].USE_KINV 			= USE_KINV;
			clusters[c].PAR_NUM_THREADS 	= PAR_NUM_THREADS;
			clusters[c].SOLVER_NUM_THREADS  = SOLVER_NUM_THREADS;
		}


//			std::vector<eslocal> domain_list(number_of_subdomains_per_supercluster, 0);
//			for (int i = 0; i < number_of_subdomains_per_supercluster; i++) {
//				domain_list[i] = i;
//			}


		mtype 	  = instance->K[0].mtype; // TODO: WARNING - Needs to be fixed

		switch (mtype) {
		case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
			SYMMETRIC_SYSTEM = true;
			break;
		case MatrixType::REAL_SYMMETRIC_INDEFINITE:
			SYMMETRIC_SYSTEM = true;
			break;
		case MatrixType::REAL_UNSYMMETRIC:
			SYMMETRIC_SYSTEM = false;
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Unknown matrix type";
		}


		// MPI calls to get total number of compute clusters - not superclusters
		int global_numClusters = 0;
		MPI_Allreduce(&numClusters, &global_numClusters, 1, MPI_INT, MPI_SUM, environment->MPICommunicator);
		int glob_clust_index = 0;
		MPI_Exscan(&numClusters, &glob_clust_index, 1, MPI_INT, MPI_SUM, environment->MPICommunicator);

		//instance->computeKernels(configuration.regularization, configuration.SC_SIZE);

		for (eslocal c = 0; c < numClusters; c++) {

			eslocal number_of_subdomains_per_cluster = 0;
			std::vector<eslocal> domain_list;
			for (eslocal i = 0; i < instance->clustersMap.size(); i++) {
				if (instance->clustersMap[i] == c) {
					number_of_subdomains_per_cluster++;
					domain_list.push_back(i);
				}
			}

			clusters[c].cluster_global_index = glob_clust_index + c + 1;
			clusters[c].cluster_local_index  = c;
			//clusters[c]->my_neighs = std::vector<eslocal>(instance->neighbours.begin(), instance->neighbours.end());
			clusters[c].mtype = mtype; // TODO: WARNING - Needs to be fixed

			switch (clusters[c].mtype) {
			case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
				clusters[c].SYMMETRIC_SYSTEM = true;
				break;
			case MatrixType::REAL_SYMMETRIC_INDEFINITE:
				clusters[c].SYMMETRIC_SYSTEM = true;
				break;
			case MatrixType::REAL_UNSYMMETRIC:
				clusters[c].SYMMETRIC_SYSTEM = false;
				break;
			default:
				ESINFO(GLOBAL_ERROR) << "Unknown matrix type";
			}

			// Init all compute clusters - communication layer and respective buffers are not allocated
			clusters[c].InitClusterPC(&domain_list[0], number_of_subdomains_per_cluster);

			// Get an original mapping of the subdomains
			for (int d = 0; d < clusters[c].domains.size(); d++) {
				domains[clusters[c].domains[d].domain_global_index] = & clusters[c].domains[d];
				x_prim_cluster1[clusters[c].domains[d].domain_global_index] = & clusters[c].x_prim_cluster1[d];
				x_prim_cluster2[clusters[c].domains[d].domain_global_index] = & clusters[c].x_prim_cluster2[d];
				x_prim_cluster3[clusters[c].domains[d].domain_global_index] = & clusters[c].x_prim_cluster3[d];

			}

		}


		// Setup communication layer of the supercluster
		MPIrank   = environment->MPIrank;
		my_neighs = std::vector<eslocal>(instance->neighbours.begin(), instance->neighbours.end());
		SetupCommunicationLayer();


		//TODO - Fix and remove
	}

	void SetClusterHFETI() {

		ESINFO(PROGRESS3) << "HFETI preprocessing start";

//		for (eslocal c = 0; c < clusters.size(); c++) {
//			clusters[c].SetClusterHFETI();
//		}

		for (eslocal c = 0; c < clusters.size(); c++) {
			if (clusters[c].domains.size() == 1) {
				//ESINFO(ALWAYS) << Info::TextColor::YELLOW
				std::cout << "Cluster " << clusters[c].cluster_global_index << " on MPI rank " << environment->MPIrank << " has only one domain -> Using TFETI" << std::endl;
				clusters[c].USE_HFETI = 0;
			}
		}

		 TimeEval HFETI_prec_timing (" HFETI - preprocessing timing");
		 HFETI_prec_timing.totalTime.start();

		 TimeEvent B0_time("Create and Compress B0 per cluster"); B0_time.start();

//		 instance->computeKernels(configuration.regularization, configuration.SC_SIZE);

		 instance->assembleB0(configuration.B0_type, instance->N1);

//#pragma omp parallel for
//for (size_t d = 0; d < domains.size(); d++) {
//	domains[d]->multKplusLocal(*x_in[d]);
//}

		for (eslocal c = 0; c < clusters.size(); c++) {
			if (clusters[c].domains.size() > 1) clusters[c].CompressB0();
		}
		 B0_time.end(); B0_time.printStatMPI(); HFETI_prec_timing.addEvent(B0_time);

		 for (eslocal c = 0; c < clusters.size(); c++) {

			 // *** Alocate temporarly vectors for inter-cluster processing *********************
			 // *** - based on uncompressed matrix B0
			 clusters[c].tm1.resize(clusters[c].domains.size());
			 clusters[c].tm2.resize(clusters[c].domains.size());
			 clusters[c].tm3.resize(clusters[c].domains.size());

			#pragma omp parallel for
			 for (size_t d = 0; d < clusters[c].domains.size(); d++) {

				 eslocal max_tmp_vec_size = clusters[c].domains[d].B0.cols;
				 if (clusters[c].domains[d].B0.rows > clusters[c].domains[d].B0.cols)
					 max_tmp_vec_size = clusters[c].domains[d].B0.rows;

				 clusters[c].tm1[d].resize( max_tmp_vec_size );
				 clusters[c].tm2[d].resize( max_tmp_vec_size );
				 clusters[c].tm3[d].resize( max_tmp_vec_size );
			 }

		 }


		 TimeEvent G0_time("Create G0 per cluster"); G0_time.start();
		for (eslocal c = 0; c < clusters.size(); c++) {
			if (clusters[c].domains.size() > 1) clusters[c].CreateG0();
		}
		 G0_time.end(); G0_time.printStatMPI(); HFETI_prec_timing.addEvent(G0_time);


		 TimeEvent F0_time("Create F0 per cluster"); F0_time.start();
		for (eslocal c = 0; c < clusters.size(); c++) {
			if (clusters[c].domains.size() > 1) clusters[c].CreateF0();
		}
		 F0_time.end(); HFETI_prec_timing.addEvent(F0_time);

		 TimeEvent Sa_time("Create Salfa per cluster"); Sa_time.start();
		for (eslocal c = 0; c < clusters.size(); c++) {
			if (clusters[c].domains.size() > 1) clusters[c].CreateSa();
		}
		 Sa_time.end(); HFETI_prec_timing.addEvent(Sa_time);

		 HFETI_prec_timing.totalTime.end();
		 HFETI_prec_timing.printStatsMPI();



	}

	SparseMatrix G1, G1_comp;
	SparseMatrix G2, G2_comp;
	SparseMatrix GGtinvM;

	void Create_G_perCluster() {
		for (eslocal c = 0; c < numClusters; c++) {
			clusters[c].Create_G_perCluster();
			G1.MatAppend(clusters[c].G1);
			G2.MatAppend(clusters[c].G2);
		}
	}


	void Compress_G1() {

		G1_comp.Clear();
		Compress_G(G1, G1_comp);
		G1.Clear();

		if (!SYMMETRIC_SYSTEM) {
			G2_comp.Clear();
			Compress_G(G2, G2_comp);
			G2.Clear();
		}

	}

	void Compress_G( SparseMatrix &G_in, SparseMatrix &G_comp_out ) {

		G_in.ConvertToCOO( 1 );

		#pragma omp parallel for
		for (size_t j = 0; j < G_in.J_col_indices.size(); j++ ) {
			G_in.J_col_indices[j] = my_lamdas_map_indices[ G_in.J_col_indices[j] -1 ] + 1;  // numbering from 1 in matrix
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


	// *** Communication layer setup members and routines

	SEQ_VECTOR < SEQ_VECTOR <eslocal> >	my_comm_lambdas_indices;
	SEQ_VECTOR < SEQ_VECTOR <double> >  my_comm_lambdas;
	SEQ_VECTOR < SEQ_VECTOR <double> >  my_recv_lambdas;

	SEQ_VECTOR < SEQ_VECTOR <eslocal> >	my_comm_lambdas_indices_comp;
	SEQ_VECTOR < SEQ_VECTOR <double> >  my_comm_lambdas_comp;
	SEQ_VECTOR < SEQ_VECTOR <double> >  my_recv_lambdas_comp;

	SEQ_VECTOR <eslocal>  my_neighs;
	SEQ_VECTOR <eslocal>  my_lamdas_indices;
	map <eslocal,eslocal> my_lamdas_map_indices;

	SEQ_VECTOR <eslocal>  my_lamdas_ddot_filter;
	SEQ_VECTOR <eslocal>  lambdas_filter;

	SEQ_VECTOR <double> compressed_tmp;

	void SetupCommunicationLayer() {

		SEQ_VECTOR <SEQ_VECTOR <eslocal> > & lambda_map_sub = instance->B1clustersMap;

		//// *** Detection of affinity of lag. multipliers to specific subdomains **************
		//// *** - will be used to compress vectors and matrices for higher efficiency

		ESLOG(MEMORY) << "Setting vectors for lambdas";
		ESLOG(MEMORY) << "process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		my_lamdas_indices.resize( lambda_map_sub.size() );
		for (size_t i = 0; i < lambda_map_sub.size(); i++)
			my_lamdas_indices[i] = lambda_map_sub[i][0];


		SEQ_VECTOR< SEQ_VECTOR <eslocal> > lambdas_per_subdomain (instance->K.size() * environment->MPIsize ); // ( domains.size() * environment->MPIsize);

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

		#pragma omp parallel for
		for (size_t i = 0; i < my_neighs.size(); i++) {
			my_comm_lambdas_indices[i] = lambdas_per_subdomain[my_neighs[i]];
			my_comm_lambdas[i].resize(my_comm_lambdas_indices[i].size());
			my_recv_lambdas[i].resize(my_comm_lambdas_indices[i].size());
		}

		// mapping/compression vector for cluster
		for (size_t i = 0; i <my_lamdas_indices.size(); i++)
			my_lamdas_map_indices.insert(make_pair(my_lamdas_indices[i],i));

		//// *** Create a vector of communication pattern needed for AllReduceLambdas function *******
		my_comm_lambdas_indices_comp.resize(my_neighs.size());
		#pragma omp parallel for
		for (size_t i = 0; i < my_neighs.size(); i++) {
			my_comm_lambdas_indices_comp[i].resize( lambdas_per_subdomain[my_neighs[i]].size() );
			for (size_t j = 0; j < lambdas_per_subdomain[my_neighs[i]].size(); j++ )
				my_comm_lambdas_indices_comp[i][j] = my_lamdas_map_indices[lambdas_per_subdomain[my_neighs[i]][j]];
		}
		//// *** END - Create a vector of communication pattern needed for AllReduceLambdas function *

		// Temp buffer for dual buffers
		compressed_tmp.resize( lambda_map_sub.size(), 0 );
		dual_size = my_lamdas_indices.size();

	}

	void compress_lambda_vector  ( SEQ_VECTOR <double> & decompressed_vec_lambda )
	{
		//compress vector for CG in main loop
		for (size_t i = 0; i < my_lamdas_indices.size(); i++)
			decompressed_vec_lambda[i] = decompressed_vec_lambda[my_lamdas_indices[i]];

		decompressed_vec_lambda.resize(my_lamdas_indices.size());
	}

	void decompress_lambda_vector( SEQ_VECTOR <double> &   compressed_vec_lambda )
	{
		SEQ_VECTOR <double> decompressed_vec_lambda (domains[0]->B1.rows,0); //TODO - need fix

		for (size_t i = 0; i < my_lamdas_indices.size(); i++)
			decompressed_vec_lambda[my_lamdas_indices[i]] = compressed_vec_lambda[i];

		compressed_vec_lambda = decompressed_vec_lambda;
	}

    SEQ_VECTOR <double> vec_d;
	void CreateVec_d_perCluster( SEQ_VECTOR<SEQ_VECTOR <double> > & f ) {
		vec_d.clear();
		for (int c = 0; c < numClusters; c++) {
			clusters[c].CreateVec_d_perCluster ( f );
			vec_d.insert(vec_d.end(), clusters[c].vec_d.begin(), clusters[c].vec_d.end());
		}
	}

    SEQ_VECTOR <double> vec_b;
    SEQ_VECTOR <double> vec_b_compressed;
	void CreateVec_b_perCluster( SEQ_VECTOR<SEQ_VECTOR <double> > & f )  {
		vec_b_compressed.clear();
		vec_b_compressed.resize(dual_size, 0.0);


		for (int c = 0; c < numClusters; c++) {
			clusters[c].vec_b_compressed.clear();
			clusters[c].vec_b_compressed.resize(dual_size, 0.0);
			clusters[c].CreateVec_b_perCluster ( f );
			for (int i=0; i< dual_size; i++) {
				vec_b_compressed[i] += clusters[c].vec_b_compressed[i];	// Can be improved - no copy is needed
			}
		}

	}

    //SEQ_VECTOR <double> vec_lb_out;
	void CreateVec_lb_perCluster( SEQ_VECTOR <double> & vec_lb_out )  {
		//TODO: Need check
		vec_lb_out.clear();
		vec_lb_out.resize(my_lamdas_indices.size(), -std::numeric_limits<double>::infinity());
		for (int c = 0; c < numClusters; c++) {
			clusters[c].CreateVec_lb_perCluster ( vec_lb_out );
		}

	}

	void CreateVec_c_perCluster( SEQ_VECTOR <double> & vec_c_out )  {
		//TODO: Need check
		vec_c_out.clear();
		vec_c_out.resize(my_lamdas_indices.size(), 0.0);
		for (int c = 0; c < numClusters; c++) {
			clusters[c].CreateVec_c_perCluster ( vec_c_out );
		}
	}

	SEQ_VECTOR<SEQ_VECTOR<double> *> x_prim_cluster1;
	SEQ_VECTOR<SEQ_VECTOR<double> *> x_prim_cluster2;
	SEQ_VECTOR<SEQ_VECTOR<double> *> x_prim_cluster3;

	void multKplusHFETI   (SEQ_VECTOR<SEQ_VECTOR<double> *> & x_in) { multKplusGlobal_l(x_in); }
	void multKplusGlobal_l(SEQ_VECTOR<SEQ_VECTOR<double> *> & x_in) {

		SEQ_VECTOR<SEQ_VECTOR<double> > x_prim_cluster;
		x_prim_cluster.resize(number_of_subdomains_per_supercluster);

		for (int c = 0; c < numClusters; c++) {

			for (int d = 0; d < clusters[c].domains.size(); d++)
				x_prim_cluster[d].swap( *x_in[clusters[c].domains[d].domain_global_index] );

			if (clusters[c].USE_HFETI == 1) {
				clusters[c].multKplusGlobal_l( x_prim_cluster );
			} else {
				clusters[c].domains[0].multKplusLocal(x_prim_cluster[0]); // FIX for cases where HTFETI cluster has only one subdomain
			}

			for (int d = 0; d < clusters[c].domains.size(); d++)
				(*x_in[clusters[c].domains[d].domain_global_index]).swap( x_prim_cluster[d] );

		}
	}

	void multKplusHFETI   (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) { multKplusGlobal_l(x_in); }
	void multKplusGlobal_l(SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {

		SEQ_VECTOR<SEQ_VECTOR<double> > x_prim_cluster;
		x_prim_cluster.resize(number_of_subdomains_per_supercluster);

		for (int c = 0; c < numClusters; c++) {
			for (int d = 0; d < clusters[c].domains.size(); d++)
				x_prim_cluster[d].swap( x_in[clusters[c].domains[d].domain_global_index] );

			if (clusters[c].USE_HFETI == 1) {
				clusters[c].multKplusGlobal_l( x_prim_cluster );
			} else {
				clusters[c].domains[0].multKplusLocal(x_prim_cluster[0]); // FIX for cases where HTFETI cluster has only one subdomain
			}

			for (int d = 0; d < clusters[c].domains.size(); d++)
				x_in[clusters[c].domains[d].domain_global_index].swap ( x_prim_cluster[d] );
		}
	}


	void multKplusGlobal_Kinv (SEQ_VECTOR<SEQ_VECTOR<double> *> & x_in) {
		SEQ_VECTOR<SEQ_VECTOR<double> > x_prim_cluster;
		x_prim_cluster.resize(number_of_subdomains_per_supercluster);

		for (int c = 0; c < numClusters; c++) {

			for (int d = 0; d < clusters[c].domains.size(); d++)
				x_prim_cluster[d].swap( *x_in[clusters[c].domains[d].domain_global_index] );

			clusters[c].multKplusGlobal_Kinv( x_prim_cluster );

			for (int d = 0; d < clusters[c].domains.size(); d++)
				(*x_in[clusters[c].domains[d].domain_global_index]).swap( x_prim_cluster[d] );

		}

	}

	void multKplusHFETI_LSC   (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) { multKplusGlobal_Kinv(x_in); };
	void multKplusGlobal_Kinv (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {
		SEQ_VECTOR<SEQ_VECTOR<double> > x_prim_cluster;
		x_prim_cluster.resize(number_of_subdomains_per_supercluster);

		for (int c = 0; c < numClusters; c++) {
			for (int d = 0; d < clusters[c].domains.size(); d++)
				x_prim_cluster[d].swap( x_in[clusters[c].domains[d].domain_global_index] );

			clusters[c].multKplusGlobal_Kinv( x_prim_cluster );

			for (int d = 0; d < clusters[c].domains.size(); d++)
				x_in[clusters[c].domains[d].domain_global_index].swap ( x_prim_cluster[d] );
		}
	}


	void multKplusFETI(SEQ_VECTOR<SEQ_VECTOR<double> *> & x_in) {
		#pragma omp parallel for
		for (size_t d = 0; d < domains.size(); d++) {
			domains[d]->multKplusLocal(*x_in[d]);
		}
	}

	void multKplusFETI(SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {
		#pragma omp parallel for
		for (size_t d = 0; d < domains.size(); d++) {
			domains[d]->multKplusLocal(x_in[d]);
		}
	}

	void multKplus(SEQ_VECTOR<SEQ_VECTOR<double> > & x_in) {

		if (  USE_HFETI ) multKplusHFETI(x_in);
		if ( !USE_HFETI ) multKplusFETI (x_in);


//		if ( USE_HFETI && !USE_KINV) multKplusHFETI(x_in);
//		if (!USE_HFETI && !USE_KINV) multKplusFETI (x_in);
//
//		if ( USE_HFETI &&  USE_KINV) multKplusHFETI_LSC(x_in);
//		if (!USE_HFETI &&  USE_KINV) multKplusFETI_LSC (x_in);


	}



};

}


#endif /* SRC_SOLVER_SPECIFIC_SUPERCLUSTER_H_ */
