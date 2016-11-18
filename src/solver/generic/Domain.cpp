#include "../generic/Domain.h"

#include <cmath>

// *******************************************************************
// **** DOMAIN CLASS ************************************************

using namespace espreso;

Domain::Domain(){

}

Domain::Domain(eslocal domain_index, eslocal use_dynamic_1_no_dynamic_0){
	domain_global_index = domain_index;
	USE_DYNAMIC         = use_dynamic_1_no_dynamic_0;
	isOnACC             = 0;
}

void Domain::SetDynamicParameters(double set_dynamic_timestep, double set_dynamic_beta, double set_dynamic_gama) {
	dynamic_timestep = set_dynamic_timestep;
	dynamic_beta     = set_dynamic_beta;
	dynamic_gama     = set_dynamic_gama;
}

void Domain::SetDomain(eslocal USE_HFETI, eslocal use_dynamic_1_no_dynamic_0) {

	USE_DYNAMIC = use_dynamic_1_no_dynamic_0;
	//K_regularizationFromR( );
	domain_prim_size = Kplus.cols;

}

void Domain::CreateKplus_R ( std::vector < std::vector < double > > coordinates )
{
	ESINFO(GLOBAL_ERROR) << "Create Kplus_R";
//	if (DOFS_PER_NODE == 3) {
//
//		eslocal elem_index = 0;
//
//		SparseMatrix R_per_element;
//
//		R_per_element.rows = 3;
//		R_per_element.cols = 2;
//		R_per_element.nnz  = 2;
//		R_per_element.type = 'G';
//
//		R_per_element.CSR_I_row_indices.resize(4);
//		R_per_element.CSR_J_col_indices.resize(2);
//		R_per_element.CSR_V_values.resize(2);
//
//		R_per_element.CSR_V_values[0] = 1;
//		R_per_element.CSR_V_values[1] = 1;
//
//		R_per_element.CSR_J_col_indices[0] = 1;
//		R_per_element.CSR_J_col_indices[1] = 2;
//
//		R_per_element.CSR_I_row_indices[0] = 1;
//		R_per_element.CSR_I_row_indices[1] = 2;
//		R_per_element.CSR_I_row_indices[2] = 3;
//		R_per_element.CSR_I_row_indices[3] = 3;
//
//		for (elem_index = 0; elem_index < coordinates.size(); elem_index++) {
//			Kplus_R.MatAppend(R_per_element);
//		}
//
//		R_per_element.Clear();
//
//		Kplus_R.ConvertCSRToDense(1);

//		eslocal elem_index = 0;
//
//		SparseMatrix R_per_element;
//
//		R_per_element.rows = 3;
//		R_per_element.cols = 6;
//		R_per_element.nnz  = 9;
//		R_per_element.type = 'G';
//
//		R_per_element.CSR_I_row_indices.resize(4);
//		R_per_element.CSR_J_col_indices.resize(9);
//		R_per_element.CSR_V_values.resize(9);
//
//		R_per_element.CSR_V_values[0] =   1;
//		R_per_element.CSR_V_values[1] = - coordinates[elem_index][2];
//		R_per_element.CSR_V_values[2] =   coordinates[elem_index][1];
//		R_per_element.CSR_V_values[3] =   1;
//		R_per_element.CSR_V_values[4] =   coordinates[elem_index][2];
//		R_per_element.CSR_V_values[5] = - coordinates[elem_index][0];
//		R_per_element.CSR_V_values[6] =   1 ;
//		R_per_element.CSR_V_values[7] = - coordinates[elem_index][1];
//		R_per_element.CSR_V_values[8] =   coordinates[elem_index][0];
//
//		R_per_element.CSR_J_col_indices[0] = 1;
//		R_per_element.CSR_J_col_indices[1] = 5;
//		R_per_element.CSR_J_col_indices[2] = 6;
//		R_per_element.CSR_J_col_indices[3] = 2;
//		R_per_element.CSR_J_col_indices[4] = 4;
//		R_per_element.CSR_J_col_indices[5] = 6;
//		R_per_element.CSR_J_col_indices[6] = 3;
//		R_per_element.CSR_J_col_indices[7] = 4;
//		R_per_element.CSR_J_col_indices[8] = 5;
//
//		R_per_element.CSR_I_row_indices[0] = 1;
//		R_per_element.CSR_I_row_indices[1] = 4;
//		R_per_element.CSR_I_row_indices[2] = 7;
//		R_per_element.CSR_I_row_indices[3] = 10;
//
//		Kplus_R.MatAppend(R_per_element);
//
//		for (elem_index = 1; elem_index < coordinates.size(); elem_index++) {
//			R_per_element.CSR_V_values[1] = - coordinates[elem_index][2];
//			R_per_element.CSR_V_values[2] =   coordinates[elem_index][1];
//			R_per_element.CSR_V_values[4] =   coordinates[elem_index][2];
//			R_per_element.CSR_V_values[5] = - coordinates[elem_index][0];
//			R_per_element.CSR_V_values[7] = - coordinates[elem_index][1];
//			R_per_element.CSR_V_values[8] =   coordinates[elem_index][0];
//
//			Kplus_R.MatAppend(R_per_element);
//		}
//
//		R_per_element.Clear();
//
//		Kplus_R.ConvertCSRToDense(1);
//
//	}
//
//
//	if (DOFS_PER_NODE == 1) {
//
//		Kplus_R.dense_values.resize( coordinates.size() );
//		double nsqrt = 1.0 / sqrt( coordinates.size() );
//
//		for (eslocal elem_index = 0; elem_index < coordinates.size(); elem_index++) {
//			Kplus_R.dense_values[elem_index] = nsqrt;
//		}
//
//		Kplus_R.rows = coordinates.size();
//		Kplus_R.cols = 1;
//		Kplus_R.nnz  = coordinates.size();
//		Kplus_R.type = 'G';
//
//		//Kplus_R.ConvertDenseToCSR(1);
//	}
}

void Domain::multKplusLocal(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, eslocal x_in_vector_start_index, eslocal y_out_vector_start_index) {
	switch (config::solver::KSOLVER) {
	case config::solver::KSOLVERalternative::DIRECT_DP:
		Kplus.Solve(x_in, y_out, x_in_vector_start_index, y_out_vector_start_index);
		break;
	case config::solver::KSOLVERalternative::DIRECT_SP:
		Kplus.Solve(x_in, y_out, x_in_vector_start_index, y_out_vector_start_index);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Invalid KSOLVER value.";
		exit(EXIT_FAILURE);
	}


	//		SEQ_VECTOR<double> x (Kplus.m_Kplus_size, 0.0);
	//		SEQ_VECTOR<double> r (Kplus.m_Kplus_size, 0.0);
	//		SEQ_VECTOR<double> z (Kplus.m_Kplus_size, 0.0);
	//
	//		Kplus.Solve(x_in, x, x_in_vector_start_index, 0);
	//		if (enable_SP_refinement) {
	//			for (eslocal step = 0; step <= esconfig::solver::KSOLVER_SP_iter_steps; step++) {
	//				K.MatVec(x,r,'N');
	//				for (eslocal i = 0; i < r.size(); i++)
	//					r[i] = x_in[i + x_in_vector_start_index] - r[i];
	//				Kplus.Solve(r, z, 0, 0);
	//				for (eslocal i = 0; i < r.size(); i++)
	//					x[i] = x[i] + z[i];
	//
	//				double norm = 0.0;
	//				for (eslocal i = 0; i < r.size(); i++)
	//					norm += r[i]*r[i];
	//
	//				norm = sqrt(norm);
	//
	//				if (norm < esconfig::solver::KSOLVER_SP_iter_norm) {
	//					std::cout.precision(20);
	//					std::cout << "Refinement steps: " << step << " | norm: " << norm << std::endl;
	//					break;
	//				}
	//
	//			}
	//		}
	//
	//		for (eslocal i = 0; i < r.size(); i++)
	//			y_out[y_out_vector_start_index + i] = x[i];

//	case 1: {
//		Kplus.SolveCG(K, x_in_y_out);
//		break;
//	}

}

void Domain::multKplusLocal(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out) {
	switch (config::solver::KSOLVER) {
	case config::solver::KSOLVERalternative::DIRECT_DP:
		Kplus.Solve(x_in, y_out, 0, 0);
		break;
	case config::solver::KSOLVERalternative::ITERATIVE:
		Kplus.SolveCG(K, x_in, y_out);
		break;
	case config::solver::KSOLVERalternative::DIRECT_SP:
		Kplus.Solve(x_in, y_out, 0, 0);
		break;
	case config::solver::KSOLVERalternative::DIRECT_MP: {

		SEQ_VECTOR<double> x (Kplus.m_Kplus_size, 0.0);
		SEQ_VECTOR<double> r (Kplus.m_Kplus_size, 0.0);
		SEQ_VECTOR<double> z (Kplus.m_Kplus_size, 0.0);

		bool success = false;

		Kplus.Solve(x_in, x, 0, 0);
		if (enable_SP_refinement) {
			for (size_t step = 0; step <= config::solver::KSOLVER_SP_STEPS; step++) {
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

				if (norm < config::solver::KSOLVER_SP_NORM) {
					ESINFO(PROGRESS2) << " " << step;
					success = true;
					break;
				}
			}
		}

		if (!success) {
			ESINFO(PROGRESS2) << "FAILED";
		}

		for (size_t i = 0; i < r.size(); i++) {
			y_out[i] = x[i];
		}
		break;
	}
//	case 4:
//		SEQ_VECTOR<double> x (Kplus.m_Kplus_size, 0.0);
//		Kplus.Solve(x_in, x, 0, 0);
//		Kplus.SolveCG(K, x_in, y_out, x);
//		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Invalid KSOLVER value.";
		exit(EXIT_FAILURE);
	}
}

void Domain::multKplusLocal(SEQ_VECTOR <double> & x_in_y_out) {
	switch (config::solver::KSOLVER) {
	case config::solver::KSOLVERalternative::DIRECT_DP:
		Kplus.Solve(x_in_y_out);
		break;
	case config::solver::KSOLVERalternative::DIRECT_SP:
		Kplus.Solve(x_in_y_out);
		break;
	case config::solver::KSOLVERalternative::DIRECT_MP: {
		bool success = false;

		SEQ_VECTOR<double> x (Kplus.m_Kplus_size, 0.0);
		SEQ_VECTOR<double> r (Kplus.m_Kplus_size, 0.0);
		SEQ_VECTOR<double> z (Kplus.m_Kplus_size, 0.0);

		Kplus.Solve(x_in_y_out, x, 0, 0);

		if (enable_SP_refinement) {
			for (size_t step = 0; step <= config::solver::KSOLVER_SP_STEPS; step++) {
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

				if (norm < config::solver::KSOLVER_SP_NORM) {
					ESINFO(PROGRESS2) << " " << step;
					break;
				}

			}
		}

		if (!success) {
			ESINFO(PROGRESS2) << "FAILED";
		}

		for (size_t i = 0; i < r.size(); i++) {
			x_in_y_out[i] = x[i];
		}

		break;
	}
//	case 4: { // DIRECT MIX - 2xSP
//
//		SEQ_VECTOR<double> x (Kplus.m_Kplus_size, 0.0);
//		SEQ_VECTOR<double> z (Kplus.m_Kplus_size, 0.0);
//
//		Kplus.Solve(x_in_y_out, x, 0, 0);
//		Kplus.SolveCG(K, x_in_y_out, z, x);
//
//		for (eslocal i = 0; i < z.size(); i++)
//			x_in_y_out[i] = z[i];
//
//		break;
//	}
	case config::solver::KSOLVERalternative::ITERATIVE:
		Kplus.SolveCG(K, x_in_y_out);
		break;
	default:
		ESINFO(ERROR) << "Invalid KSOLVER value.";
		exit(EXIT_FAILURE);
	}
}

//void Domain::K_regularization ( )
//{
//
//	SparseMatrix N;
//	N.rows = 0;
//	N.cols = 6; // pozor
//	N.nnz  = 0;
//	N.type = 'G';
//	eslocal row_fill = 1;
//	N.CSR_I_row_indices.reserve(K.rows+1);
//
//	eslocal old_index  = 0;
//	eslocal next_index = 0;
//
//	SEQ_VECTOR <double> dense_values;
//	double dense_pattern[] =
//	   {1, 0, 0,    0, 8, 8,
//		0, 1, 0,    8, 0, 8,
//		0, 0, 1,    8, 8, 0}; // 8 is ust a random number
//	eslocal dense_pattern_size = 18;
//
//	dense_values.resize(fix_nodes.size() * dense_pattern_size);
//	for (eslocal fpi = 0; fpi < fix_nodes.size(); fpi++) {
//		eslocal elem_index = fix_nodes[fpi];
//		dense_pattern[4]  = - coordinates[elem_index][2];
//		dense_pattern[5]  =   coordinates[elem_index][1];
//		dense_pattern[9]  =   coordinates[elem_index][2];
//		dense_pattern[11] = - coordinates[elem_index][0];
//		dense_pattern[15] = - coordinates[elem_index][1];
//		dense_pattern[16] =   coordinates[elem_index][0];
//		copy(dense_pattern, dense_pattern+dense_pattern_size, dense_values.begin() + (fpi*dense_pattern_size));
//		//dense_values.assign(dense_pattern + (fpi*dense_pattern_size), dense_pattern + ((fpi+1)*dense_pattern_size));
//	}
//
//	for (eslocal fpi = 0; fpi < fix_nodes.size(); fpi++) {
//
//		old_index  = next_index;
//		next_index = 3 * fix_nodes[fpi] - 2;
//		N.CSR_I_row_indices.resize( next_index );
//		fill(N.CSR_I_row_indices.begin() + old_index, N.CSR_I_row_indices.begin() + next_index, row_fill);
//		N.rows = next_index;
//
//		eslocal elem_index = fix_nodes[fpi];
//		SparseMatrix R_per_element;
//
//
//		R_per_element.rows =  3;
//		R_per_element.cols =  6;
//		R_per_element.nnz  =  9;
//		R_per_element.type = 'G';
//
//		R_per_element.CSR_I_row_indices.resize(4);
//		R_per_element.CSR_J_col_indices.resize(9);
//		R_per_element.CSR_V_values     .resize(9);
//
//		R_per_element.CSR_V_values[0] =   1;
//		R_per_element.CSR_V_values[1] = - coordinates[elem_index][2];
//		R_per_element.CSR_V_values[2] =   coordinates[elem_index][1];
//		R_per_element.CSR_V_values[3] =   1;
//		R_per_element.CSR_V_values[4] =   coordinates[elem_index][2];
//		R_per_element.CSR_V_values[5] = - coordinates[elem_index][0];
//		R_per_element.CSR_V_values[6] =   1 ;
//		R_per_element.CSR_V_values[7] = - coordinates[elem_index][1];
//		R_per_element.CSR_V_values[8] =   coordinates[elem_index][0];
//
//		R_per_element.CSR_J_col_indices[0] = 1;
//		R_per_element.CSR_J_col_indices[1] = 5;
//		R_per_element.CSR_J_col_indices[2] = 6;
//		R_per_element.CSR_J_col_indices[3] = 2;
//		R_per_element.CSR_J_col_indices[4] = 4;
//		R_per_element.CSR_J_col_indices[5] = 6;
//		R_per_element.CSR_J_col_indices[6] = 3;
//		R_per_element.CSR_J_col_indices[7] = 4;
//		R_per_element.CSR_J_col_indices[8] = 5;
//
//		R_per_element.CSR_I_row_indices[0] = 1;
//		R_per_element.CSR_I_row_indices[1] = 4;
//		R_per_element.CSR_I_row_indices[2] = 7;
//		R_per_element.CSR_I_row_indices[3] = 10;
//
//		N.MatAppend(R_per_element);
//		next_index = next_index + 3;
//		row_fill = N.CSR_I_row_indices[N.CSR_I_row_indices.size()-1];
//
//	}
//
//	old_index  = next_index;
//	next_index = K.rows;
//	N.CSR_I_row_indices.resize( next_index + 1);
//	fill(N.CSR_I_row_indices.begin() + old_index, N.CSR_I_row_indices.begin() + next_index + 1, row_fill);
//	N.rows = next_index;
//
//	SparseMatrix Nt;
//	N.MatTranspose(Nt);
//
//	SparseMatrix NtN_Mat;
//	NtN_Mat.MatMat(Nt,'N',N);
//	NtN_Mat.MatTranspose();
//	NtN_Mat.RemoveLower();
//
//	SparseSolverCPU NtN;
//	NtN.ImportMatrix(NtN_Mat);
//	NtN_Mat.Clear();
//
//	NtN.Factorization();
//	NtN.SolveMat_Sparse(Nt);
//
//	//NtN = Nt*N
//	//N.Clear();
//	//Nt.MatTranspose(N);
//	NtN_Mat.MatMat(N,'N',Nt);
//	NtN_Mat.RemoveLower();
//
//	double ro = K.GetMaxOfDiagonalOfSymmetricMatrix();
//	ro = 1.0 * ro;
//
//	K.    MatAddInPlace (NtN_Mat,'N', ro);
//	Kplus.ImportMatrix  (K);
//	Kplus.Factorization ();
//
//	KplusF.ImportMatrix (K);
//	//SpyText(K);
//	K.Clear();
//
//}

void Domain::K_regularizationFromR ( SparseMatrix & K_in ) {

	ESINFO(GLOBAL_ERROR) << "K_regularizationFromR";
//    if (USE_DYNAMIC == 0) {

//    	if (DOFS_PER_NODE == 3) {
//
//			SparseMatrix N;
//
//			Kplus_R.ConvertDenseToCSR(0);
//
//			N.CreateMatFromRowsFromMatrix( Kplus_R, fix_dofs);
//
//			SEQ_VECTOR<eslocal>().swap( Kplus_R.CSR_I_row_indices );
//			SEQ_VECTOR<eslocal>().swap( Kplus_R.CSR_J_col_indices );
//			SEQ_VECTOR<double> ().swap( Kplus_R.CSR_V_values );
//
//			SparseMatrix Nt;
//			N.MatTranspose( Nt );
//
//			SparseMatrix NtN_Mat;
//			NtN_Mat.MatMat( Nt,'N',N );
//			NtN_Mat.MatTranspose();
//			NtN_Mat.RemoveLower();
//
//			SparseSolverCPU NtN;
//			NtN.ImportMatrix(NtN_Mat);
//			NtN_Mat.Clear();
//
//			std::stringstream ss;
//			ss << "K regularization from R -> rank: " << config::env::MPIrank;
//			NtN.Factorization(ss.str());
//			NtN.SolveMat_Sparse(Nt);
//			NtN.Clear();
//
//			NtN_Mat.MatMat(N,'N',Nt);
//			NtN_Mat.RemoveLower();
//
//			double ro = K_in.GetMaxOfDiagonalOfSymmetricMatrix();
//			ro = 1.0 * ro;
//
//			_RegMat = NtN_Mat;
//			_RegMat.MatScale(ro);
//
//			K_in.MatAddInPlace(_RegMat, 'N', 1.0);
//
//			_RegMat.ConvertToCOO(1);
//
////			K.    MatAddInPlace (NtN_Mat,'N', ro);
//    	}
//
//    	if (DOFS_PER_NODE == 3) {
//    		double ro = K_in.GetMaxOfDiagonalOfSymmetricMatrix();
//
//    		//K.CSR_V_values[ K.CSR_I_row_indices[100] - 1 ] = +ro;
//
//    		//K_in.CSR_V_values[0] += ro;
//
//    		_RegMat.type = K_in.type;
//    		_RegMat.cols = K_in.cols;
//    		_RegMat.rows = K_in.rows;
//    		_RegMat.nnz  = 2;
//
//    		_RegMat.I_row_indices.push_back(1);
//    		_RegMat.J_col_indices.push_back(1);
//    		_RegMat.I_row_indices.push_back(2);
//			_RegMat.J_col_indices.push_back(2);
//    		_RegMat.V_values.push_back(ro);
//    		_RegMat.V_values.push_back(ro);
//
//    		_RegMat.ConvertToCSR(0);
//    		K_in.MatAddInPlace(_RegMat, 'N', ro);
//    	}
//
//    	if (DOFS_PER_NODE == 1) {
//    		double ro = K_in.GetMaxOfDiagonalOfSymmetricMatrix();
//
//    		//K.CSR_V_values[ K.CSR_I_row_indices[100] - 1 ] = +ro;
//
//    		K_in.CSR_V_values[0] += ro;
//
//    		_RegMat.type = K_in.type;
//    		_RegMat.cols = K_in.cols;
//    		_RegMat.rows = K_in.rows;
//    		_RegMat.nnz  = 1;
//
//    		_RegMat.I_row_indices.push_back(1);
//    		_RegMat.J_col_indices.push_back(1);
//    		_RegMat.V_values.push_back(ro);
//
//    	}
//
//    }

}

// **** END - DOMAIN CLASS *******************************************
// *******************************************************************
