// Just for testing to get the matrix kernel using dissection
// #include <Driver/DissectionSolver.hpp>

#include "clustercpu.h"

#include "../../../assembler/instance.h"

using namespace espreso;

void ClusterCPU::Create_SC_perDomain(bool USE_FLOAT) {

    #pragma omp parallel for
	for (size_t i = 0; i < domains_in_global_index.size(); i++ ) {
        domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);
    }

    ESINFO(PROGRESS3) << "Creating B1*K+*B1t : using Pardiso SC";

    #pragma omp parallel for
    for (size_t i = 0; i < domains_in_global_index.size(); i++ ) {
        SparseSolverMKL tmpsps;
        if ( i == 0 && cluster_global_index == 1) {
            tmpsps.msglvl = Info::report(LIBRARIES) ? 1 : 0;
        }
        tmpsps.Create_SC_w_Mat( domains[i].K, domains[i].B1t_comp_dom, domains[i].B1Kplus, false, 1 );

        if (USE_FLOAT){
            domains[i].B1Kplus.ConvertDenseToDenseFloat( 1 );
            domains[i].B1Kplus.USE_FLOAT = true;
        }
        ESINFO(PROGRESS3) << Info::plain() << ".";
    }
    ESINFO(PROGRESS3);


    #pragma omp parallel for
    for (size_t i = 0; i < domains_in_global_index.size(); i++ )
        domains[i].B1t_comp_dom.Clear();
	}

void ClusterCPU::Create_Kinv_perDomain() {

	#pragma omp parallel for
	for (size_t i = 0; i < domains_in_global_index.size(); i++ )
        domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);

    ESINFO(PROGRESS3) << "Creating B1*K+*B1t";

    #pragma omp parallel for
    for (size_t i = 0; i < domains_in_global_index.size(); i++ ) {

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
        ESINFO(PROGRESS3) << Info::plain() << ".";
    }
    ESINFO(PROGRESS3);

    #pragma omp parallel for
	for (size_t i = 0; i < domains_in_global_index.size(); i++ )
			domains[i].B1t_comp_dom.Clear();
	}


void ClusterCPU::SetupKsolvers ( ) {

//   	#pragma omp parallel for
//	for (size_t d = 0; d < domains.size(); d++) {
//
//		domains[d].enable_SP_refinement = true;
//#if 0
//		// Just for testing to get the matrix kernel using dissection
//		domains[d]._RegMat.ConvertToCSR(0);
//		domains[d].K.MatAddInPlace(domains[d]._RegMat, 'N', -1.0);
//
//		domains[d].Kplus.ImportMatrix_wo_Copy(domains[d].K);
//		domains[d].Kplus.Factorization ("Dissection - kernel");
//		SEQ_VECTOR <double> kern_vect;
//		eslocal kern_dim = 0;
//		domains[d].Kplus.GetKernelVectors(kern_vect, kern_dim);
//
//		domains[d].Kplus.Clear();
//		domains[d].K.MatAddInPlace(domains[d]._RegMat, 'N', 1.0);
//		domains[d]._RegMat.ConvertToCOO(1);
//#endif
//        // Import of Regularized matrix K into Kplus (Sparse Solver)
//    	switch (configuration.Ksolver) {
//		case ESPRESO_KSOLVER::DIRECT_DP:
//			domains[d].Kplus.ImportMatrix_wo_Copy (domains[d].K);
//			break;
//		case ESPRESO_KSOLVER::ITERATIVE:
//			domains[d].Kplus.ImportMatrix_wo_Copy (domains[d].K);
//			break;
//		case ESPRESO_KSOLVER::DIRECT_SP:
//			domains[d].Kplus.ImportMatrix_wo_Copy_fl(domains[d].K);
//			break;
//		case ESPRESO_KSOLVER::DIRECT_MP:
//			domains[d].Kplus.ImportMatrix_fl(domains[d].K);
//			break;
////		case 4:
////			domains[d].Kplus.ImportMatrix_fl(domains[d].K);
////			break;
//		default:
//			ESINFO(ERROR) << "Invalid KSOLVER value.";
//		}
//
//        //domains[d].Kplus.mtype = -2;
//
//        if (configuration.keep_factors) {
//            std::stringstream ss;
//            ss << "init -> rank: " << environment->MPIrank << ", subdomain: " << d;
//            domains[d].Kplus.keep_factors = true;
//            if (configuration.Ksolver != ESPRESO_KSOLVER::ITERATIVE) {
//                domains[d].Kplus.Factorization (ss.str());
//            }
//        } else {
//            domains[d].Kplus.keep_factors = false;
//            domains[d].Kplus.MPIrank = environment->MPIrank;
//        }
//
//        //TODO: Hot Fix - needs to be done better
//#ifdef ESBEM
//        if ( !SYMMETRIC_SYSTEM ) {
//            // 11 = Real and unsymmetric matrix
//            domains[d].Kplus.mtype = espreso::MatrixType::REAL_UNSYMMETRIC; //11;
//        } else {
//            // 2 = Real and symmetric positive definite
//            domains[d].Kplus.mtype = espreso::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE; //2;
//        }
//        //TODO: else stokes = -2 = Real and symmetric indefinite
//#else
//        if ( !SYMMETRIC_SYSTEM ) {
//            // 11 = Real and unsymmetric matrix
//            domains[d].Kplus.mtype = 11; //espreso::MatrixType::REAL_UNSYMMETRIC; //11;
//        } else {
//            // 2 = Real and symmetric positive definite
//            domains[d].Kplus.mtype = 2; //espreso::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE; //2;
//        }
//        //TODO: else stokes = -2 = Real and symmetric indefinite
//#endif
//
//        if ( d == 0 && environment->MPIrank == 0) {
//        	domains[d].Kplus.msglvl = 0;
//        }
//        ESINFO(PROGRESS3) << Info::plain() << ".";
//    }
//    ESINFO(PROGRESS3);

}


void ClusterCPU::CreateDirichletPrec( Instance *instance ) {

#pragma omp parallel for
for (size_t d = 0; d < domains.size(); d++) {
//for (size_t d = 0; d < instance->K.size(); d++) {

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
#ifdef ESBEM
//TODO: Alex - da se nejak spocist Dir prec z dense matic K ??
        K_modif.ConvertDenseToCSR(1);
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
			SparseSolverMKL createSchur;
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
			SparseSolverMKL K_rr_solver;

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

			KsrInvKrrKrs.MatMat(K_sr, 'N', K_rs);
			domains[d].Prec.MatAddInPlace(KsrInvKrrKrs, 'N', -1);
//          if (!diagonalized_K_rr){
//				    cluster.domains[d].Prec.ConvertCSRToDense(1);
//          }
		}

	}

	if (environment->print_matrices) {
		std::ofstream osS(Logging::prepareFile(d, "S"));
		SparseMatrix SC = domains[d].Prec;
		if (configuration.preconditioner == ESPRESO_PRECONDITIONER::DIRICHLET) {
			SC.ConvertDenseToCSR(1);
		}
		osS << SC;
		osS.close();
	}

	ESINFO(PROGRESS3) << Info::plain() << ".";

}
}
