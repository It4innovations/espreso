
#include "clustercpu.h"

void ClusterCPU::Create_SC_perDomain(bool USE_FLOAT) {

	cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ )
		domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);

	if (cluster_global_index == 1)
		cout << "Creating B1*K+*B1t : using Pardiso SC : ";

	cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ ) {

		if (cluster_global_index == 1) cout << "."; // << i ;

		SparseSolverCPU tmpsps;
		if ( i == 0 && cluster_global_index == 1) tmpsps.msglvl = 1;
		tmpsps.Create_SC_w_Mat( domains[i].K, domains[i].B1t_comp_dom, domains[i].B1Kplus, false, 1 );

		if (USE_FLOAT){
			domains[i].B1Kplus.ConvertDenseToDenseFloat( 1 );
			domains[i].B1Kplus.USE_FLOAT = true;
		}


	}


	cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ )
		domains[i].B1t_comp_dom.Clear();

	if (cluster_global_index == 1)
		cout << endl;

}

void ClusterCPU::SetupKsolvers ( ) {

	cilk_for (eslocal d = 0; d < domains.size(); d++) {

		// Import of Regularized matrix K into Kplus (Sparse Solver)
		switch (esconfig::solver::KSOLVER) {
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
			ESLOG(eslog::ERROR) << "Invalid KSOLVER value.";
			exit(EXIT_FAILURE);
		}

		if (esconfig::solver::KEEP_FACTORS == 1) {
			std::stringstream ss;
			ss << "init -> rank: " << esconfig::MPIrank << ", subdomain: " << d;
			domains[d].Kplus.keep_factors = true;
			if (esconfig::solver::KSOLVER != 1) {
				domains[d].Kplus.Factorization (ss.str());
			}
		} else {
			domains[d].Kplus.keep_factors = false;
			domains[d].Kplus.MPIrank = esconfig::MPIrank;
		}

		domains[d].domain_prim_size = domains[d].Kplus.cols;

		if ( d == 0 && esconfig::MPIrank == 0) domains[d].Kplus.msglvl=0;
		if (esconfig::MPIrank == 0) std::cout << ".";

	}

}
