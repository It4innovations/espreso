
#include "Domain.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"

#include <cmath>

// *******************************************************************
// **** DOMAIN CLASS ************************************************

using namespace espreso;


void storeData(SparseMatrix matrices, const std::string &name, const std::string &description, int d)
{
	if (info::ecf->output.print_matrices) {
		eslog::storedata("STORING '%s'.\n", description.c_str());
		//for (size_t d = 0; d < matrices.size(); d++) {
//          std::ofstream os(Logging::prepareFile(d, name));
//          os << matrices;
//          os.close();
		//}
	}
}

void storeData(vector<double> vectors, const std::string &name, const std::string &description, int d)
{
	if (info::ecf->output.print_matrices) {
		eslog::storedata("STORING '%s'.\n", description.c_str());
		//for (size_t d = 0; d < vectors.size(); d++) {
//          std::ofstream os(Logging::prepareFile(d, name));
//          for (size_t i = 0; i < vectors.size(); i++) {
//              os << vectors[i] << std::endl;
//          }
//          os.close();
		//}
	}
}

Domain::Domain(const FETIConfiguration &configuration, DataHolder *instance_in, esint domain_index_in, esint USE_HTFETI_in):
		configuration(configuration),
		instance(instance_in),

		origK(instance_in->origK[domain_index_in]),
		K(instance_in->K[domain_index_in]),

		Kplus_R(instance_in->N1[domain_index_in]),
		Kplus_R2(instance_in->N2[domain_index_in]),

		Kplus_origR(instance_in->origKN1[domain_index_in]),
		Kplus_origR2(instance_in->origKN2[domain_index_in]),

		_RegMat(instance_in->RegMat[domain_index_in]),

		f(instance_in->F[domain_index_in]),

		B0(instance_in->B0[domain_index_in]),
		vec_c(instance_in->B1c[domain_index_in]),
		vec_lb(instance_in->LB[domain_index_in])
{
		domain_prim_size    = K.cols;
		domain_global_index = domain_index_in;
		USE_HFETI           = USE_HTFETI_in;

		// TODO: this is broken (ask Lubos)
		norm_f = 0;
}

void Domain::SetDomain() {

#if defined(SOLVER_DISSECTION)

	if ( configuration.regularization == FETIConfiguration::REGULARIZATION::ANALYTIC ) {
		instance->computeKernel(configuration.regularization, configuration.sc_size, domain_global_index, configuration.method == FETIConfiguration::METHOD::HYBRID_FETI);
		Kplus.ImportMatrix_wo_Copy(K);
		Kplus.Factorization ("K matrix");
	} else {
		//TODO: Modify for non-symmetric systems with two kernels
		Kplus.ImportMatrix_wo_Copy(K);
		Kplus.Factorization ("K matrix");
		Kplus.GetKernel(Kplus_R, Kplus_R2);
	}

// *** Debug only code ***
//
//  storeData(Kplus_R,              "R_dom_",      "R_dom",       domain_index);
//  storeData(Kplus_R.dense_values, "R_dom_dense", "R_dom_dense", domain_index);
//
//  Kplus_R.Clear();
//
//  instance->computeKernel(configuration.regularization, configuration.sc_size, domain_global_index, configuration.method == FETIConfiguration::METHOD::HYBRID_FETI);
//
//  Kplus.ImportMatrix(K); //_wo_Copy(K);
//  Kplus.Factorization ("K matrix");
//
//  Kplus_R.GramSchmidtOrtho();
//  Kplus_R.ConvertDenseToCSR(0);
//
//  storeData(Kplus_R,              "R_dom_a",       "R_dom_a",       domain_index);
//  storeData(Kplus_R.dense_values, "R_dom_a_dense", "R_dom_a_dense", domain_index);

#else

	switch (configuration.regularization) {
	case FETIConfiguration::REGULARIZATION::ANALYTIC:
		Kplus.ImportMatrix_wo_Copy(K);
		break;
	case FETIConfiguration::REGULARIZATION::ALGEBRAIC:
		if (
				configuration.conjugate_projector == FETIConfiguration::CONJ_PROJECTOR::CONJ_R ||
				configuration.conjugate_projector == FETIConfiguration::CONJ_PROJECTOR::CONJ_K) {

			double norm;
			esint defect;
			switch (K.mtype) {
			case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
			case MatrixType::REAL_SYMMETRIC_INDEFINITE:
				K.get_kernel_from_K(origK, _RegMat, Kplus_R, norm, defect, domain_index, configuration.sc_size);
			break;
			case MatrixType::REAL_UNSYMMETRIC:
				origK.get_kernels_from_nonsym_K(origK, _RegMat, Kplus_R, Kplus_R2, norm, defect, domain_index, configuration.sc_size);

			break;
			}
			Kplus.ImportMatrix_wo_Copy(origK);
		} else {
			double norm;
			esint defect;
			switch (K.mtype) {
			case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
			case MatrixType::REAL_SYMMETRIC_INDEFINITE:
				K.get_kernel_from_K(K, _RegMat, Kplus_R, norm, defect, domain_index, configuration.sc_size);
			break;
			case MatrixType::REAL_UNSYMMETRIC:
				K.get_kernels_from_nonsym_K(K, _RegMat, Kplus_R, Kplus_R2, norm, defect, domain_index, configuration.sc_size);
			break;
			}
			Kplus.ImportMatrix_wo_Copy(K);
		}
		break;
	}

	Kplus.Factorization ("K matrix");

//  if ( configuration.regularization == FETIConfiguration::REGULARIZATION::ANALYTIC ) {
//      instance->computeKernel(configuration.regularization, configuration.sc_size, domain_global_index, configuration.method == FETIConfiguration::METHOD::HYBRID_FETI);
//      Kplus.ImportMatrix_wo_Copy(K);
//      Kplus.Factorization ("K matrix");
//
//      if (    configuration.conjugate_projector == FETIConfiguration::CONJ_PROJECTOR::CONJ_R ||
//              configuration.conjugate_projector == FETIConfiguration::CONJ_PROJECTOR::CONJ_K)
//          instance->computeKernelFromOrigK(configuration.regularization, configuration.sc_size, domain_global_index, configuration.method == FETIConfiguration::METHOD::HYBRID_FETI);
//
//  } else {
//
//      // TODO: Change MKL solver so that it can extract kernels inside the object
//      //Kplus.ImportMatrix_wo_Copy(K);
//      //Kplus.Factorization ("K matrix");
//      //Kplus.GetKernel(Kplus_R); // TODO: Kplus.GetKernels(Kplus_R, Kplus_R2) - upravit na tuto funkci - v sym. pripade bude Kplus_R2 prazdna
//
//      // TODO: Temporary solution before MKL solver is updated
//      instance->computeKernel(configuration.regularization, configuration.sc_size, domain_global_index, configuration.method == FETIConfiguration::METHOD::HYBRID_FETI);
//      Kplus.ImportMatrix_wo_Copy(K);
//      Kplus.Factorization ("K matrix");
//
//      if (    configuration.conjugate_projector == FETIConfiguration::CONJ_PROJECTOR::CONJ_R ||
//              configuration.conjugate_projector == FETIConfiguration::CONJ_PROJECTOR::CONJ_K)
//          instance->computeKernelFromOrigK(configuration.regularization, configuration.sc_size, domain_global_index, configuration.method == FETIConfiguration::METHOD::HYBRID_FETI);
//
//  }



#endif


//  Kplus_R.ConvertDenseToCSR(0);
//
//  std::cout << K.SpyText();
//  std::cout << Kplus_R.SpyText();
//
//  storeData(K,"K_dom_","K_dom",domain_global_index);
//
//  storeData(Kplus_R,              "R_dom_a",       "R_dom_a",       domain_index);
//  storeData(Kplus_R.dense_values, "R_dom_a_dense", "R_dom_a_dense", domain_index);
//
//  std::cout << "Norm K_R alex = " << K.getNorm_K_R(K,Kplus_R,'N') << "\n";
//
//  _RegMat.ConvertToCSR(0);
//  K.MatAddInPlace(_RegMat,'N',-1.0);
//
//  std::cout << "Norm K_R alex = " << K.getNorm_K_R(K,Kplus_R,'N') << "\n";
//
//  SparseMatrix Kfull;
//  Kfull = K;
//  Kfull.MatTranspose();
//  Kfull.SetDiagonalOfSymmetricMatrix(0.0);
//  Kfull.MatAddInPlace(K,'N',1.0);
//  Kfull.type='G';
//  std::cout << Kfull.SpyText();
//
//
//  SparseMatrix TMa;
//  TMa.MatMat(Kfull,'T',Kplus_R);
//  storeData(TMa,              "KtimesRa",       "KtimesRa",       domain_index);
//
//  Kplus.ImportMatrix(K); //  _wo_Copy(K);
//  Kplus.Factorization ("Dissection - kernel");
//
//
//
////    Kplus_R.Clear();
////    Kplus.GetKernel(Kplus_R);
////    Kplus_R.ConvertDenseToCSR(0);
//
//
//  std::cout << "Norm K_R diss = " << K.getNorm_K_R(K,Kplus_R,'N') << "\n";
//
//  K.MatAddInPlace(_RegMat,'N',+1.0);
//
//  SparseMatrix TM;
//  TM.MatMat(Kfull,'T',Kplus_R);
//  storeData(TM,              "KtimesR",       "KtimesR",       domain_index);
//
//
////    std::cout << K.SpyText();
////    std::cout << Kplus_R.SpyText();
//
//
//
//  storeData(Kplus_R,              "R_dom_",      "R_dom",       domain_index);
//  storeData(Kplus_R.dense_values, "R_dom_dense", "R_dom_dense", domain_index);


	// *** Kernel setup
	if ( configuration.orthogonal_K_kernels ) {
		Kplus_R.GramSchmidtOrtho();
		if (Kplus_R2.nnz > 0)
			Kplus_R2.GramSchmidtOrtho();
	}

	Kplus_Rb  = Kplus_R;
	Kplus_Rb2 = Kplus_R2;
	// *** END - Kernel setup


	// *** Constraints and Dirichlet boundary condition
	B1 = instance->B1[domain_global_index];
	B1.type = 'G';
	B1t = B1;
	B1t.MatTransposeCOO();
	B1t.ConvertToCSRwithSort(1);

	B1_scale_vec = instance->B1duplication[domain_global_index];
	for (size_t i = 0; i < instance->B1[domain_global_index].I_row_indices.size(); ++i) {
		if (i == 0 || instance->B1[domain_global_index].I_row_indices[i - 1] != instance->B1[domain_global_index].I_row_indices[i]) {
			lambda_map_sub.push_back(instance->B1[domain_global_index].I_row_indices[i] - 1);
		}
	}
	// *** END - Constraints and Dirichlet boundary condition
}


void Domain::multKplusLocal(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out) {

	if (configuration.mp_pseudoinverse) {

		//            // a = ( x )
		//            // b =  Kplus * ( (( a )) - R * (RT * ( a ) ) )
		//            // result = (  b - R * ( RT * ( b ) ) )
		//            if (configuration.mp_pseudoinverse) {
		//              cluster.domains[d].Kplus_R.DenseMatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d], 'T');
		//              cluster.domains[d].Kplus_R.DenseMatVec(cluster.x_prim_cluster2[d], cluster.x_prim_cluster1[d], 'N', 0, 0, 1.0, -1.0);
		//            }

		SEQ_VECTOR<double> mp_tmp;
		mp_tmp.resize(Kplus_R.cols);

		y_out = x_in;

		Kplus_R.DenseMatVec(y_out, mp_tmp, 'T');
		Kplus_R.DenseMatVec(mp_tmp, y_out, 'N', 0, 0, 1.0, -1.0);

		multKplusLocalCore(y_out);

		Kplus_R.DenseMatVec(y_out, mp_tmp, 'T');
		Kplus_R.DenseMatVec(mp_tmp, y_out, 'N', 0, 0, 1.0, -1.0);



	} else {
		multKplusLocalCore(x_in, y_out);
	}
}


void Domain::multKplusLocal(SEQ_VECTOR <double> & x_in_y_out) {

	if (configuration.mp_pseudoinverse) {

		SEQ_VECTOR<double> mp_tmp;
		mp_tmp.resize(Kplus_R.cols);

		Kplus_R.DenseMatVec(x_in_y_out, mp_tmp, 'T');
		Kplus_R.DenseMatVec(mp_tmp, x_in_y_out, 'N', 0, 0, 1.0, -1.0);

		multKplusLocalCore(x_in_y_out);

		Kplus_R.DenseMatVec(x_in_y_out, mp_tmp, 'T');
		Kplus_R.DenseMatVec(mp_tmp, x_in_y_out, 'N', 0, 0, 1.0, -1.0);

	} else {
		multKplusLocalCore(x_in_y_out);
	}
}

void Domain::multKplusLocalCore(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out) {


	switch (configuration.Ksolver) {
	case FETIConfiguration::KSOLVER::DIRECT_DP:
		Kplus.Solve(x_in, y_out, 0, 0);
		break;
	case FETIConfiguration::KSOLVER::ITERATIVE:
		Kplus.SolveCG(K, x_in, y_out);
		break;
	case FETIConfiguration::KSOLVER::DIRECT_SP:
		Kplus.Solve(x_in, y_out, 0, 0);
		break;
	case FETIConfiguration::KSOLVER::DIRECT_MP: {

		SEQ_VECTOR<double> x (Kplus.m_Kplus_size, 0.0);
		SEQ_VECTOR<double> r (Kplus.m_Kplus_size, 0.0);
		SEQ_VECTOR<double> z (Kplus.m_Kplus_size, 0.0);

		Kplus.Solve(x_in, x, 0, 0);
		if (enable_SP_refinement) {
			for (size_t step = 0; step <= configuration.Ksolver_max_iterations; step++) {
				K.MatVec(x,r,'N');
				for (size_t i = 0; i < r.size(); i++) {
					r[i] = x_in[i] - r[i];
				}
				Kplus.Solve(r, z, 0, 0);
				//Kplus.SolveCG(K, r, z);
				for (size_t i = 0; i < r.size(); i++) {
					x[i] = x[i] + z[i];
				}

				double norm = 0.0;
				for (size_t i = 0; i < r.size(); i++) {
					norm += r[i]*r[i];
				}

				norm = sqrt(norm);

				if (norm < configuration.Ksolver_precision) {
					break;
				}
			}
		}

		for (size_t i = 0; i < r.size(); i++) {
			y_out[i] = x[i];
		}
		break;
	}
//  case 4:
//      SEQ_VECTOR<double> x (Kplus.m_Kplus_size, 0.0);
//      Kplus.Solve(x_in, x, 0, 0);
//      Kplus.SolveCG(K, x_in, y_out, x);
//      break;
	default:
		eslog::globalerror("Invalid KSOLVER value.\n");
		exit(EXIT_FAILURE);
	}
}

void Domain::multKplusLocalCore(SEQ_VECTOR <double> & x_in_y_out) {

	switch (configuration.Ksolver) {
	case FETIConfiguration::KSOLVER::DIRECT_DP:
		Kplus.Solve(x_in_y_out);
		break;
	case FETIConfiguration::KSOLVER::DIRECT_SP:
		Kplus.Solve(x_in_y_out);
		break;
	case FETIConfiguration::KSOLVER::DIRECT_MP: {
		SEQ_VECTOR<double> x (Kplus.m_Kplus_size, 0.0);
		SEQ_VECTOR<double> r (Kplus.m_Kplus_size, 0.0);
		SEQ_VECTOR<double> z (Kplus.m_Kplus_size, 0.0);

		Kplus.Solve(x_in_y_out, x, 0, 0);

		if (enable_SP_refinement) {
			for (size_t step = 0; step <= configuration.Ksolver_max_iterations; step++) {
				K.MatVec(x,r,'N');
				for (size_t i = 0; i < r.size(); i++) {
					r[i] = x_in_y_out[i] - r[i];
				}
				Kplus.Solve(r, z, 0, 0);
				for (size_t i = 0; i < r.size(); i++) {
					x[i] = x[i] + z[i];
				}

				double norm = 0.0;
				for (size_t i = 0; i < r.size(); i++) {
					norm += r[i]*r[i];
				}

				norm = sqrt(norm);

				if (norm < configuration.Ksolver_precision) {
					break;
				}

			}
		}

		for (size_t i = 0; i < r.size(); i++) {
			x_in_y_out[i] = x[i];
		}

		break;
	}
//  case 4: { // DIRECT MIX - 2xSP
//
//      SEQ_VECTOR<double> x (Kplus.m_Kplus_size, 0.0);
//      SEQ_VECTOR<double> z (Kplus.m_Kplus_size, 0.0);
//
//      Kplus.Solve(x_in_y_out, x, 0, 0);
//      Kplus.SolveCG(K, x_in_y_out, z, x);
//
//      for (esint i = 0; i < z.size(); i++)
//          x_in_y_out[i] = z[i];
//
//      break;
//  }
	case FETIConfiguration::KSOLVER::ITERATIVE:
		Kplus.SolveCG(K, x_in_y_out);
		break;
	default:
		eslog::globalerror("Invalid KSOLVER value.\n");
		exit(EXIT_FAILURE);
	}
}




// TODO: Obsolete functions - to be removed


void Domain::multKplusLocal(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, esint x_in_vector_start_index, esint y_out_vector_start_index) {

	eslog::globalerror("MultKplus-local with 4 parameters - currently dissabled.\n");

	switch (configuration.Ksolver) {
	case FETIConfiguration::KSOLVER::DIRECT_DP:
		Kplus.Solve(x_in, y_out, x_in_vector_start_index, y_out_vector_start_index);
		break;
	case FETIConfiguration::KSOLVER::DIRECT_SP:
		Kplus.Solve(x_in, y_out, x_in_vector_start_index, y_out_vector_start_index);
		break;
	default:
		eslog::globalerror("Invalid KSOLVER value.\n");
		exit(EXIT_FAILURE);
	}

}

















//void Domain::multKplusLocal(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, esint x_in_vector_start_index, esint y_out_vector_start_index) {
//  switch (configuration.Ksolver) {
//  case FETIConfiguration::KSOLVER::DIRECT_DP:
//      Kplus.Solve(x_in, y_out, x_in_vector_start_index, y_out_vector_start_index);
//      break;
//  case FETIConfiguration::KSOLVER::DIRECT_SP:
//      Kplus.Solve(x_in, y_out, x_in_vector_start_index, y_out_vector_start_index);
//      break;
//  default:
//      //ESINFO(GLOBAL_ERROR) << "Invalid KSOLVER value.";
//      exit(EXIT_FAILURE);
//  }
//
//
//  //      SEQ_VECTOR<double> x (Kplus.m_Kplus_size, 0.0);
//  //      SEQ_VECTOR<double> r (Kplus.m_Kplus_size, 0.0);
//  //      SEQ_VECTOR<double> z (Kplus.m_Kplus_size, 0.0);
//  //
//  //      Kplus.Solve(x_in, x, x_in_vector_start_index, 0);
//  //      if (enable_SP_refinement) {
//  //          for (esint step = 0; step <= esconfiguration.Ksolver_SP_iter_steps; step++) {
//  //              K.MatVec(x,r,'N');
//  //              for (esint i = 0; i < r.size(); i++)
//  //                  r[i] = x_in[i + x_in_vector_start_index] - r[i];
//  //              Kplus.Solve(r, z, 0, 0);
//  //              for (esint i = 0; i < r.size(); i++)
//  //                  x[i] = x[i] + z[i];
//  //
//  //              double norm = 0.0;
//  //              for (esint i = 0; i < r.size(); i++)
//  //                  norm += r[i]*r[i];
//  //
//  //              norm = sqrt(norm);
//  //
//  //              if (norm < esconfiguration.Ksolver_SP_iter_norm) {
//  //                  std::cout.precision(20);
//  //                  std::cout << "Refinement steps: " << step << " | norm: " << norm << std::endl;
//  //                  break;
//  //              }
//  //
//  //          }
//  //      }
//  //
//  //      for (esint i = 0; i < r.size(); i++)
//  //          y_out[y_out_vector_start_index + i] = x[i];
//
////    case 1: {
////        Kplus.SolveCG(K, x_in_y_out);
////        break;
////    }
//
//}
//
//void Domain::multKplusLocal(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out) {
//  switch (configuration.Ksolver) {
//  case FETIConfiguration::KSOLVER::DIRECT_DP:
//      Kplus.Solve(x_in, y_out, 0, 0);
//      break;
//  case FETIConfiguration::KSOLVER::ITERATIVE:
//      Kplus.SolveCG(K, x_in, y_out);
//      break;
//  case FETIConfiguration::KSOLVER::DIRECT_SP:
//      Kplus.Solve(x_in, y_out, 0, 0);
//      break;
//  case FETIConfiguration::KSOLVER::DIRECT_MP: {
//
//      SEQ_VECTOR<double> x (Kplus.m_Kplus_size, 0.0);
//      SEQ_VECTOR<double> r (Kplus.m_Kplus_size, 0.0);
//      SEQ_VECTOR<double> z (Kplus.m_Kplus_size, 0.0);
//
//      bool success = false;
//
//      Kplus.Solve(x_in, x, 0, 0);
//      if (enable_SP_refinement) {
//          for (size_t step = 0; step <= configuration.Ksolver_max_iterations; step++) {
//              K.MatVec(x,r,'N');
//              for (size_t i = 0; i < r.size(); i++) {
//                  r[i] = x_in[i] - r[i];
//              }
//              Kplus.Solve(r, z, 0, 0);
//              //Kplus.SolveCG(K, r, z);
//              for (size_t i = 0; i < r.size(); i++) {
//                  x[i] = x[i] + z[i];
//              }
//
//              double norm = 0.0;
//              for (size_t i = 0; i < r.size(); i++) {
//                  norm += r[i]*r[i];
//              }
//
//              norm = sqrt(norm);
//
//              if (norm < configuration.Ksolver_precision) {
//                  //ESINFO(PROGRESS3) << " " << step;
//                  success = true;
//                  break;
//              }
//          }
//      }
//
//      if (!success) {
//          //ESINFO(PROGRESS3) << "FAILED";
//      }
//
//      for (size_t i = 0; i < r.size(); i++) {
//          y_out[i] = x[i];
//      }
//      break;
//  }
////    case 4:
////        SEQ_VECTOR<double> x (Kplus.m_Kplus_size, 0.0);
////        Kplus.Solve(x_in, x, 0, 0);
////        Kplus.SolveCG(K, x_in, y_out, x);
////        break;
//  default:
//      //ESINFO(GLOBAL_ERROR) << "Invalid KSOLVER value.";
//      exit(EXIT_FAILURE);
//  }
//}
//
//void Domain::multKplusLocal(SEQ_VECTOR <double> & x_in_y_out) {
//  switch (configuration.Ksolver) {
//  case FETIConfiguration::KSOLVER::DIRECT_DP:
//      Kplus.Solve(x_in_y_out);
//      break;
//  case FETIConfiguration::KSOLVER::DIRECT_SP:
//      Kplus.Solve(x_in_y_out);
//      break;
//  case FETIConfiguration::KSOLVER::DIRECT_MP: {
//      bool success = false;
//
//      SEQ_VECTOR<double> x (Kplus.m_Kplus_size, 0.0);
//      SEQ_VECTOR<double> r (Kplus.m_Kplus_size, 0.0);
//      SEQ_VECTOR<double> z (Kplus.m_Kplus_size, 0.0);
//
//      Kplus.Solve(x_in_y_out, x, 0, 0);
//
//      if (enable_SP_refinement) {
//          for (size_t step = 0; step <= configuration.Ksolver_max_iterations; step++) {
//              K.MatVec(x,r,'N');
//              for (size_t i = 0; i < r.size(); i++) {
//                  r[i] = x_in_y_out[i] - r[i];
//              }
//              Kplus.Solve(r, z, 0, 0);
//              for (size_t i = 0; i < r.size(); i++) {
//                  x[i] = x[i] + z[i];
//              }
//
//              double norm = 0.0;
//              for (size_t i = 0; i < r.size(); i++) {
//                  norm += r[i]*r[i];
//              }
//
//              norm = sqrt(norm);
//
//              if (norm < configuration.Ksolver_precision) {
//                  //ESINFO(PROGRESS3) << " " << step;
//                  break;
//              }
//
//          }
//      }
//
//      if (!success) {
//          //ESINFO(PROGRESS3) << "FAILED";
//      }
//
//      for (size_t i = 0; i < r.size(); i++) {
//          x_in_y_out[i] = x[i];
//      }
//
//      break;
//  }
////    case 4: { // DIRECT MIX - 2xSP
////
////        SEQ_VECTOR<double> x (Kplus.m_Kplus_size, 0.0);
////        SEQ_VECTOR<double> z (Kplus.m_Kplus_size, 0.0);
////
////        Kplus.Solve(x_in_y_out, x, 0, 0);
////        Kplus.SolveCG(K, x_in_y_out, z, x);
////
////        for (esint i = 0; i < z.size(); i++)
////            x_in_y_out[i] = z[i];
////
////        break;
////    }
//  case FETIConfiguration::KSOLVER::ITERATIVE:
//      Kplus.SolveCG(K, x_in_y_out);
//      break;
//  default:
//      //ESINFO(ERROR) << "Invalid KSOLVER value.";
//      exit(EXIT_FAILURE);
//  }
//}

// **** END - DOMAIN CLASS *******************************************
// *******************************************************************
