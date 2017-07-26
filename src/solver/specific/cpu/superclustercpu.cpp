#include "superclustercpu.h"

#include "../../../assembler/instance.h"

using namespace espreso;

void SuperClusterCPU::init() {
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


