
#include "clustercpu.h"

using namespace espreso;

void ClusterCPU::Create_SC_perDomain(bool USE_FLOAT) {

    cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ ) {
        domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);
    }

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

void ClusterCPU::Create_Kinv_perDomain() {
    cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ )
        domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);


    if (cluster_global_index == 1)
        cout << "Creating B1*K+*B1t : ";

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

    }

    cout << endl;
    cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ )
        domains[i].B1t_comp_dom.Clear();

    //	std::cout << "This function is obsolete - use Create_SC_perDomain" << std::endl;
    //	return();

}


void ClusterCPU::SetupKsolvers ( ) {

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

        if (config::solver::KEEP_FACTORS == 1) {
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

        if ( d == 0 && config::MPIrank == 0) domains[d].Kplus.msglvl=0;
        if (config::MPIrank == 0) std::cout << ".";

    }

}
