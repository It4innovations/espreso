
#include "clustercpu.h"

using namespace espreso;

void ClusterCPU::Create_SC_perDomain(bool USE_FLOAT) {

    cilk_for (size_t i = 0; i < domains_in_global_index.size(); i++ ) {
        domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);
    }

    ESINFO(PROGRESS2) << "Creating B1*K+*B1t : using Pardiso SC";

    cilk_for (size_t i = 0; i < domains_in_global_index.size(); i++ ) {
        SparseSolverCPU tmpsps;
        if ( i == 0 && cluster_global_index == 1) {
        	tmpsps.msglvl = Info::report(LIBRARIES) ? 1 : 0;
        }
        tmpsps.Create_SC_w_Mat( domains[i].K, domains[i].B1t_comp_dom, domains[i].B1Kplus, false, 1 );

        if (USE_FLOAT){
            domains[i].B1Kplus.ConvertDenseToDenseFloat( 1 );
            domains[i].B1Kplus.USE_FLOAT = true;
        }
        ESINFO(PROGRESS2) << Info::plain() << ".";
    }
    ESINFO(PROGRESS2);


    cilk_for (size_t i = 0; i < domains_in_global_index.size(); i++ )
        domains[i].B1t_comp_dom.Clear();
}

void ClusterCPU::Create_Kinv_perDomain() {
    cilk_for (size_t i = 0; i < domains_in_global_index.size(); i++ )
        domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);

    ESINFO(PROGRESS2) << "Creating B1*K+*B1t";

    cilk_for (size_t i = 0; i < domains_in_global_index.size(); i++ ) {

        domains[i].KplusF.msglvl = 0;

        if ( i == 0 && cluster_global_index == 1) {
        	domains[i].KplusF.msglvl = Info::report(LIBRARIES) ? 1 : 0;
        }

        //SolveMatF is obsolete - use Schur complement instead
        domains[i].KplusF.SolveMatF(domains[i].B1t_comp_dom, domains[i].B1Kplus, false);
        domains[i].B1Kplus.MatTranspose();

        SparseMatrix Btmp;
        Btmp.MatAddInPlace(domains[i].B1Kplus, 'N', 1.0);

        domains[i].B1Kplus.Clear ();
        domains[i].B1Kplus.MatMat(Btmp,'N', domains[i].B1t_comp_dom);
        domains[i].B1Kplus.ConvertCSRToDense(0);
        //domains[i].B1Kplus.ConvertDenseToDenseFloat(0);
        ESINFO(PROGRESS2) << Info::plain() << ".";
    }
    ESINFO(PROGRESS2);

    cilk_for (size_t i = 0; i < domains_in_global_index.size(); i++ )
        domains[i].B1t_comp_dom.Clear();
}


void ClusterCPU::SetupKsolvers ( ) {

    cilk_for (size_t d = 0; d < domains.size(); d++) {

        // Import of Regularized matrix K into Kplus (Sparse Solver)
    	switch (config::solver::KSOLVER) {
		case config::solver::KSOLVERalternative::DIRECT_DP:
			domains[d].Kplus.ImportMatrix_wo_Copy (domains[d].K);
			break;
		case config::solver::KSOLVERalternative::ITERATIVE:
			domains[d].Kplus.ImportMatrix_wo_Copy (domains[d].K);
			break;
		case config::solver::KSOLVERalternative::DIRECT_SP:
			domains[d].Kplus.ImportMatrix_wo_Copy_fl(domains[d].K);
			break;
		case config::solver::KSOLVERalternative::DIRECT_MP:
			domains[d].Kplus.ImportMatrix_fl(domains[d].K);
			break;
//		case 4:
//			domains[d].Kplus.ImportMatrix_fl(domains[d].K);
//			break;
		default:
			ESINFO(ERROR) << "Invalid KSOLVER value.";
		}

    	//domains[d].Kplus.mtype = -2;

        if (config::solver::KEEP_FACTORS) {
            std::stringstream ss;
            ss << "init -> rank: " << config::env::MPIrank << ", subdomain: " << d;
            domains[d].Kplus.keep_factors = true;
            if (config::solver::KSOLVER != config::solver::KSOLVERalternative::ITERATIVE) {
                domains[d].Kplus.Factorization (ss.str());
            }
        } else {
            domains[d].Kplus.keep_factors = false;
            domains[d].Kplus.MPIrank = config::env::MPIrank;
        }

        domains[d].domain_prim_size = domains[d].Kplus.cols;

        //TODO: Hot Fix - needs to be done better
        if ( !SYMMETRIC_SYSTEM ) {
        	// 11 = Real and unsymmetric matrix
        	domains[d].Kplus.mtype = 11;
        } else {
        	// 2 = Real and symmetric positive definite
        	domains[d].Kplus.mtype = 2;
        }
        //TODO: else stokes = -2 = Real and symmetric indefinite

        if ( d == 0 && config::env::MPIrank == 0) {
        	domains[d].Kplus.msglvl = 0;
        }
        ESINFO(PROGRESS2) << Info::plain() << ".";
    }
    ESINFO(PROGRESS2);

}
