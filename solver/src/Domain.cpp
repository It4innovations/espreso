#include "Domain.h"
#include <cmath>

// *******************************************************************
// **** DOMAIN CLASS ************************************************


Domain::Domain(){

}

Domain::Domain(eslocal domain_index, eslocal use_dynamic_1_no_dynamic_0){
	domain_global_index = domain_index; 
	USE_DYNAMIC = use_dynamic_1_no_dynamic_0; 
}

void Domain::SetDynamicParameters(double set_dynamic_timestep, double set_dynamic_beta, double set_dynamic_gama) {
	dynamic_timestep = set_dynamic_timestep;
	dynamic_beta     = set_dynamic_beta;
	dynamic_gama     = set_dynamic_gama;
}

void Domain::SetDomain(eslocal USE_HFETI, eslocal use_dynamic_1_no_dynamic_0) {

	USE_DYNAMIC = use_dynamic_1_no_dynamic_0; 
	K_regularizationFromR( );
	domain_prim_size = Kplus.cols;

}

void Domain::CreateKplus_R ( ) 
{
	if (DOFS_PER_NODE == 3) {

		eslocal elem_index = 0;

		SparseMatrix R_per_element;

		R_per_element.rows = 3;
		R_per_element.cols = 6;
		R_per_element.nnz  = 9;
		R_per_element.type = 'G';

		R_per_element.CSR_I_row_indices.resize(4);
		R_per_element.CSR_J_col_indices.resize(9);
		R_per_element.CSR_V_values.resize(9);

		R_per_element.CSR_V_values[0] =   1;
		R_per_element.CSR_V_values[1] = - coordinates[elem_index][2];
		R_per_element.CSR_V_values[2] =   coordinates[elem_index][1];
		R_per_element.CSR_V_values[3] =   1;
		R_per_element.CSR_V_values[4] =   coordinates[elem_index][2];
		R_per_element.CSR_V_values[5] = - coordinates[elem_index][0];
		R_per_element.CSR_V_values[6] =   1 ;
		R_per_element.CSR_V_values[7] = - coordinates[elem_index][1];
		R_per_element.CSR_V_values[8] =   coordinates[elem_index][0];

		R_per_element.CSR_J_col_indices[0] = 1;
		R_per_element.CSR_J_col_indices[1] = 5;
		R_per_element.CSR_J_col_indices[2] = 6;
		R_per_element.CSR_J_col_indices[3] = 2;
		R_per_element.CSR_J_col_indices[4] = 4;
		R_per_element.CSR_J_col_indices[5] = 6;
		R_per_element.CSR_J_col_indices[6] = 3;
		R_per_element.CSR_J_col_indices[7] = 4;
		R_per_element.CSR_J_col_indices[8] = 5;

		R_per_element.CSR_I_row_indices[0] = 1;
		R_per_element.CSR_I_row_indices[1] = 4;
		R_per_element.CSR_I_row_indices[2] = 7;
		R_per_element.CSR_I_row_indices[3] = 10;

		Kplus_R.MatAppend(R_per_element);

		for (elem_index = 1; elem_index < coordinates.size(); elem_index++) {
			R_per_element.CSR_V_values[1] = - coordinates[elem_index][2];
			R_per_element.CSR_V_values[2] =   coordinates[elem_index][1];
			R_per_element.CSR_V_values[4] =   coordinates[elem_index][2];
			R_per_element.CSR_V_values[5] = - coordinates[elem_index][0];
			R_per_element.CSR_V_values[7] = - coordinates[elem_index][1];
			R_per_element.CSR_V_values[8] =   coordinates[elem_index][0];

			Kplus_R.MatAppend(R_per_element);
		}

		R_per_element.Clear();

	}


	if (DOFS_PER_NODE == 1) {

		Kplus_R.dense_values.resize( coordinates.size() );
		double nsqrt = 1.0 / sqrt( coordinates.size() );

		for (eslocal elem_index = 1; elem_index < coordinates.size(); elem_index++) {
			Kplus_R.dense_values[elem_index] = nsqrt;
		}

		Kplus_R.rows = coordinates.size();
		Kplus_R.cols = 1;
		Kplus_R.nnz  = coordinates.size();
		Kplus_R.type = 'G';

		Kplus_R.ConvertDenseToCSR(1);
	}
}

void Domain::multKplusLocal(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, eslocal x_in_vector_start_index, eslocal y_out_vector_start_index) {
	Kplus.Solve(x_in, y_out, x_in_vector_start_index, y_out_vector_start_index); 
}

void Domain::multKplusLocal(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out) {
	Kplus.Solve(x_in, y_out, 0, 0); 
}

void Domain::multKplusLocal(SEQ_VECTOR <double> & x_in_y_out) {
	Kplus.Solve(x_in_y_out); 
}

void Domain::K_regularization ( ) 
{

	SparseMatrix N; 
	N.rows = 0;
	N.cols = 6; // pozor 
	N.nnz  = 0; 
	N.type = 'G';
	eslocal row_fill = 1;
	N.CSR_I_row_indices.reserve(K.rows+1);

	eslocal old_index  = 0;
	eslocal next_index = 0;

	SEQ_VECTOR <double> dense_values;
	double dense_pattern[] = 
	   {1, 0, 0,    0, 8, 8,
		0, 1, 0,    8, 0, 8,
		0, 0, 1,    8, 8, 0}; // 8 is ust a random number 
	eslocal dense_pattern_size = 18;

	dense_values.resize(fix_nodes.size() * dense_pattern_size);
	for (eslocal fpi = 0; fpi < fix_nodes.size(); fpi++) {
		eslocal elem_index = fix_nodes[fpi];
		dense_pattern[4]  = - coordinates[elem_index][2];
		dense_pattern[5]  =   coordinates[elem_index][1];
		dense_pattern[9]  =   coordinates[elem_index][2];
		dense_pattern[11] = - coordinates[elem_index][0];
		dense_pattern[15] = - coordinates[elem_index][1];
		dense_pattern[16] =   coordinates[elem_index][0];
		copy(dense_pattern, dense_pattern+dense_pattern_size, dense_values.begin() + (fpi*dense_pattern_size));
		//dense_values.assign(dense_pattern + (fpi*dense_pattern_size), dense_pattern + ((fpi+1)*dense_pattern_size));
	}

	for (eslocal fpi = 0; fpi < fix_nodes.size(); fpi++) {

		old_index  = next_index;
		next_index = 3 * fix_nodes[fpi] - 2;
		N.CSR_I_row_indices.resize( next_index );
		fill(N.CSR_I_row_indices.begin() + old_index, N.CSR_I_row_indices.begin() + next_index, row_fill);
		N.rows = next_index; 

		eslocal elem_index = fix_nodes[fpi];
		SparseMatrix R_per_element; 


		R_per_element.rows =  3;
		R_per_element.cols =  6; 
		R_per_element.nnz  =  9; 
		R_per_element.type = 'G'; 

		R_per_element.CSR_I_row_indices.resize(4);
		R_per_element.CSR_J_col_indices.resize(9);
		R_per_element.CSR_V_values     .resize(9);

		R_per_element.CSR_V_values[0] =   1;
		R_per_element.CSR_V_values[1] = - coordinates[elem_index][2];
		R_per_element.CSR_V_values[2] =   coordinates[elem_index][1];
		R_per_element.CSR_V_values[3] =   1;
		R_per_element.CSR_V_values[4] =   coordinates[elem_index][2];
		R_per_element.CSR_V_values[5] = - coordinates[elem_index][0];
		R_per_element.CSR_V_values[6] =   1 ;
		R_per_element.CSR_V_values[7] = - coordinates[elem_index][1];
		R_per_element.CSR_V_values[8] =   coordinates[elem_index][0];

		R_per_element.CSR_J_col_indices[0] = 1;
		R_per_element.CSR_J_col_indices[1] = 5;
		R_per_element.CSR_J_col_indices[2] = 6;
		R_per_element.CSR_J_col_indices[3] = 2;
		R_per_element.CSR_J_col_indices[4] = 4;
		R_per_element.CSR_J_col_indices[5] = 6;
		R_per_element.CSR_J_col_indices[6] = 3;
		R_per_element.CSR_J_col_indices[7] = 4;
		R_per_element.CSR_J_col_indices[8] = 5;

		R_per_element.CSR_I_row_indices[0] = 1;
		R_per_element.CSR_I_row_indices[1] = 4;
		R_per_element.CSR_I_row_indices[2] = 7;
		R_per_element.CSR_I_row_indices[3] = 10;

		N.MatAppend(R_per_element);
		next_index = next_index + 3; 
		row_fill = N.CSR_I_row_indices[N.CSR_I_row_indices.size()-1]; 

	}

	old_index  = next_index;
	next_index = K.rows;
	N.CSR_I_row_indices.resize( next_index + 1);
	fill(N.CSR_I_row_indices.begin() + old_index, N.CSR_I_row_indices.begin() + next_index + 1, row_fill);
	N.rows = next_index; 

	SparseMatrix Nt; 
	N.MatTranspose(Nt); 

	SparseMatrix NtN_Mat; 
	NtN_Mat.MatMat(Nt,'N',N); 
	NtN_Mat.MatTranspose();
	NtN_Mat.RemoveLower();

	SparseSolver NtN; 
	NtN.ImportMatrix(NtN_Mat);
	NtN_Mat.Clear();

	NtN.Factorization();
	NtN.SolveMat_Sparse(Nt); 

	//NtN = Nt*N
	//N.Clear();
	//Nt.MatTranspose(N);
	NtN_Mat.MatMat(N,'N',Nt);
	NtN_Mat.RemoveLower();

	double ro = K.GetMaxOfDiagonalOfSymmetricMatrix(); 
	ro = 1.0 * ro; 

	K.    MatAddInPlace (NtN_Mat,'N', ro);
	Kplus.ImportMatrix  (K);
	Kplus.Factorization ();
	 
	KplusF.ImportMatrix (K);
	//SpyText(K);
	K.Clear();

}

void Domain::K_regularizationFromR ( ) {

    if (USE_DYNAMIC == 0) {

    	if (DOFS_PER_NODE == 3) {

			SparseMatrix N;

			N.CreateMatFromRowsFromMatrix( Kplus_R, fix_dofs);

			SparseMatrix Nt;
			N.MatTranspose( Nt );

			SparseMatrix NtN_Mat;
			NtN_Mat.MatMat( Nt,'N',N );
			NtN_Mat.MatTranspose();
			NtN_Mat.RemoveLower();

			SparseSolver NtN;
			NtN.ImportMatrix(NtN_Mat);
			NtN_Mat.Clear();

			NtN.Factorization();
			NtN.SolveMat_Sparse(Nt);
			NtN.Clear();

			//NtN = Nt*N
			//N.Clear();
			//Nt.MatTranspose(N);
			NtN_Mat.MatMat(N,'N',Nt);
			NtN_Mat.RemoveLower();

			double ro = K.GetMaxOfDiagonalOfSymmetricMatrix();
			ro = 1.0 * ro;

			K.    MatAddInPlace (NtN_Mat,'N', ro);
			//KplusF.ImportMatrix (K);
    	}


    	if (DOFS_PER_NODE == 1) {
    		double ro = K.GetMaxOfDiagonalOfSymmetricMatrix();

    		//K.CSR_V_values[ K.CSR_I_row_indices[100] - 1 ] = +ro;

    		K.CSR_V_values[0] += ro;
    	}

    }
	
} 


void Domain::get_kernel_from_K(){
//      
// Routine calculates kernel Kplus_R of K satisfied euqality K * Kplus_R = O,
// where O is zero matrix, and it makes the matrix K non-singular (K_reg) 
// utilizing spectral conditions of Schur complement. Then ||K-K*inv(K_reg)*K||=0.0
//
//
// rev. 2015-10-23 (A.M.)
//==============================================================================
//
//    1) diagonalScaling
//  reducing of big jump coefficient effect (TODO include diagonal scaling into whole ESPRESO)
  bool diagonalScaling = true;            

//    2) permutVectorActive
//  random selection of singular DOFs
  eslocal permutVectorActive = 1; // 0 - no permut., 1 - std::vector shuffle, 2 - generating own random sequence -

//    3) use_null_pivots_or_s_set
  // NtN_Mat from null pivots or fixing DOFs
  bool use_null_pivots_or_s_set=true;

//    4) diagonalRegularization
//  regularization only on diagonal elements (big advantage: patern of K and K_regular is the same !!!)
//  size of set 's' = defect(K)
//  It's is active, only if and only if use_null_pivots_or_s_set
  bool diagonalRegularization=true;

//    5) get_n_first_and_n_last_eigenvals_from_dense_K
// get and preslocal K eigenvalues (A is temporarily converted to dense);
  eslocal get_n_first_and_n_last_eigenvals_from_dense_K = 0;

//    6) get_n_first_and_n_last_eigenvals_from_dense_S
// get and preslocal S eigenvalues
  eslocal get_n_first_and_n_last_eigenvals_from_dense_S = 0;

//    7) plot_n_first_n_last_eigenvalues
// get of K eigenvalues (A is temporarily converted to dense matrix);
  eslocal plot_n_first_n_last_eigenvalues = 0;

  if (!use_null_pivots_or_s_set) diagonalRegularization=false;

  printf("\n\n");
  printf(" ###################################################################\n");
  printf(" #                 Get kernel of K and null pivots                 #\n");
  printf(" ###################################################################\n");
// 
//    1) COND_NUMB_FOR_SINGULAR_MATRIX
//  If cond(A) > COND_NUMB_FOR_SINGULAR_MATRIX, A is considered as singular matrix.
  double COND_NUMB_FOR_SINGULAR_MATRIX=1e13; 

//    2) CHECK_NONSING
// if CHECK_NONSING>0, checking of K_rr non-singularity is activated and it is repeated 
// (CHECK_NONSING) times.
  eslocal CHECK_NONSING=0;

//    3) MAX_SIZE_OF_DENSE_MATRIX_TO_GET_EIGS
// if size of K is less then CHECK_N..., K is converted to dense format to get eigenvalues.
  eslocal MAX_SIZE_OF_DENSE_MATRIX_TO_GET_EIGS=2500;

//    4) SC_SIZE
// specification of size of Schur complement used for detection of zero eigenvalues.
// SC_SIZE >= expected defect 'd' (e.g. in elasticity d=6).
  eslocal SC_SIZE=50;

//    5) TWENTY
// testing last TWENTY eigenvalues of S to distinguish, if d-last ones are zero or not.
  eslocal TWENTY=20;
  // TWENTY eigenvalues are ascendly ordered in d = d[0],d[1], ..., d[n-2],d[n-1]
  
//    6) JUMP_IN_EIGENVALUES_ALERTING_SINGULARITY
// if d[i]/d[i+1]< JUMP_IN..., d[i] is last nonzero eigenvalue
  double JUMP_IN_EIGENVALUES_ALERTING_SINGULARITY=1.0e-5;

  
  //TODO if K.rows<=SC_SIZE, use directly input K instead of S
  //
  //
  //##########################################################################################
  //##########################################################################################
  //##########################################################################################
  //##########################################################################################
  //
  SparseMatrix S;
  SparseMatrix K_rr;
  SparseMatrix K_rs;
  eslocal i_start = 0;
  eslocal NONSING_SIZE = K.rows - SC_SIZE - i_start;
  eslocal j_start = NONSING_SIZE;
  SEQ_VECTOR <eslocal > permVec;
  permVec.resize(K.rows);
  SEQ_VECTOR <SEQ_VECTOR<eslocal >> vec_I1_i2(K.rows, SEQ_VECTOR<eslocal >(2, 1));
  eslocal offset = K.CSR_I_row_indices[0] ? 1 : 0;
  //  
  
  double cond_of_regular_part=1e307;
  eslocal *I_row_indices_p = new eslocal [K.nnz] ;
  eslocal *J_col_indices_p = new eslocal [K.nnz] ;
  SEQ_VECTOR <eslocal > tmp_vec_s;
  tmp_vec_s.resize(SC_SIZE);
  eslocal v1, n_mv, cnt_permut_vec;
  SEQ_VECTOR <eslocal >::iterator it;
  SEQ_VECTOR <eslocal > fix_dofs;
  fix_dofs.resize(SC_SIZE);
  SparseMatrix K_modif;

  double di=1,dj=1;
  eslocal cnt_iter_check_nonsing=0;


  K_modif = K;
  
//  K.MatCondNumb(K,"K_singular",plot_n_first_n_last_eigenvalues);
  // diagonal scaling of K_modif:
  // K_modif[i,j] = K_modif[i,j]/sqrt(K_modif[i,i]*K_modif[j,j]);
  for (eslocal i = 0;i<K_modif.rows;i++){
    if (diagonalScaling) {
      di=K_modif.CSR_V_values[K_modif.CSR_I_row_indices[i]-offset];
    }
    for (eslocal j = K_modif.CSR_I_row_indices[i];j<K_modif.CSR_I_row_indices[i+1];j++){
      if (diagonalScaling) {
        dj=K_modif.CSR_V_values[
          K_modif.CSR_I_row_indices[K_modif.CSR_J_col_indices[j-offset]-offset]-offset];
      }
      K_modif.CSR_V_values[j-offset]/=sqrt(di*dj);
    }
  }

  //################################################################################# 
  if (get_n_first_and_n_last_eigenvals_from_dense_K && 
      K_modif.cols<MAX_SIZE_OF_DENSE_MATRIX_TO_GET_EIGS && cnt_iter_check_nonsing==0) {
    K_modif.ConvertCSRToDense(0);
  // EIGENVALUES AND EIGENVECTORS OF STIFFNESS MATRIX
  // Works only for matrix until size ~ 2500
    char JOBZ_ = 'V';
    char UPLO_ = 'U';
    double *WK_modif = new double[K_modif.cols]; 
    double *ZK_modif = new double[K_modif.cols*K_modif.cols]; 
    MKL_INT info;
    MKL_INT ldzA = K_modif.cols;
    info = LAPACKE_dspev (LAPACK_COL_MAJOR, JOBZ_, UPLO_,
            K_modif.cols, &(K_modif.dense_values[0]), WK_modif, ZK_modif, ldzA);
    printf("eigenvals of %s d{1:%d} and d{%d:%d}\n",
          "K",get_n_first_and_n_last_eigenvals_from_dense_K,
          K_modif.rows-get_n_first_and_n_last_eigenvals_from_dense_K+2,K_modif.rows);
    for (eslocal i = 0 ; i < K_modif.rows; i++){
      if (i < get_n_first_and_n_last_eigenvals_from_dense_K || 
            i > K_modif.rows-get_n_first_and_n_last_eigenvals_from_dense_K){
        printf("%5d:  %3.8e \n",i+1, WK_modif[i]);
      }
    }
    if (info){
      printf("info = %d\n, something wrong with Schur complement in SparseSolver::generalIinverse",info);
    }
    delete [] WK_modif;
    delete [] ZK_modif;
  }
  //################################################################################# 




  while ( cond_of_regular_part > COND_NUMB_FOR_SINGULAR_MATRIX && cnt_iter_check_nonsing<(CHECK_NONSING+1)) {
    // loop checking non-singularity of K_rr matrix
    if (cnt_iter_check_nonsing>0){
      K_modif.Clear();
      K_modif=K;
      //diagonal scaling
      for (eslocal i = 0;i<K_modif.rows;i++){
        if (diagonalScaling) {
          di=K_modif.CSR_V_values[K_modif.CSR_I_row_indices[i]-offset];
        }
        for (eslocal j = K_modif.CSR_I_row_indices[i];j<K_modif.CSR_I_row_indices[i+1];j++){
          if (diagonalScaling) {
            dj=K_modif.CSR_V_values[
              K_modif.CSR_I_row_indices[K_modif.CSR_J_col_indices[j-offset]-offset]-offset];
          }
          K_modif.CSR_V_values[j-offset]/=sqrt(di*dj);
        }
      }
    }

    if (permutVectorActive<2){
      // set row permVec = {0,1,2,3,4,...,K.rows};
      for (eslocal i=0; i<K.rows; ++i) { permVec[i]=i;} // 0 1 2 A.rows-1
    }
//
    if (permutVectorActive==1){
      srand(time(NULL));
      random_shuffle ( permVec.begin(), permVec.end() );
      sort(permVec.begin(),permVec.begin()+NONSING_SIZE);
      sort(permVec.begin()+NONSING_SIZE,permVec.end());
    }
    else if (permutVectorActive==2){
      // random permutation
      n_mv = 0;                     // n_mv = size(unique(tmp_vec_s)) has to be equal to SC_SIZE
      cnt_permut_vec=0;
      srand(time(NULL));
      // loop controls, if series 'tmp_vec_s' with unique integers has suffisciant dimension.
      // If not, missing numbers are added and checked again.  
      do {
        for (eslocal i = 0;i<(SC_SIZE-n_mv);i++){
          v1 = rand() % K_modif.rows;
          tmp_vec_s[n_mv+i]=v1; 
        }
        it=tmp_vec_s.begin();
        std::sort(tmp_vec_s.begin(), tmp_vec_s.end());
        it = std::unique(tmp_vec_s.begin(), tmp_vec_s.end());
        n_mv = distance(tmp_vec_s.begin(),it);
        cnt_permut_vec++;
     } while (n_mv != SC_SIZE && cnt_permut_vec < 100);
      //
      eslocal ik=0,cnt_i=0;
      for (eslocal i = 0;i<permVec.size();i++){
        if (i==tmp_vec_s[ik]){
          permVec[ik+NONSING_SIZE]=tmp_vec_s[ik];
          ik++;
        }
        else{
          permVec[cnt_i]=i;
          cnt_i++;
        }
      }
      printf("n_mv: %d, SC_SIZE: %d, it. for RAND: %d\n",n_mv,SC_SIZE,cnt_permut_vec);
    }
    //      r = permVec[0:NONSING_SIZE-1]     (singular DOFs)
    //      s = permVec[NONSING_SIZE:end-1]   (non-singular DOFs)
    if (permutVectorActive>0){
//
      for (eslocal i = 0; i < K.rows;i++){
        vec_I1_i2[i][0]=permVec[i];
        vec_I1_i2[i][1]=i; // position to create revers permutation
      }

      std::sort(vec_I1_i2.begin(), vec_I1_i2.end(),
                [](const SEQ_VECTOR <eslocal >& a, const SEQ_VECTOR <eslocal >& b) {
        return a[0] < b[0];
      });
      // permutations made on matrix in COO format 
      K_modif.ConvertToCOO(0);
      eslocal I_index,J_index;
      for (eslocal i = 0;i<K_modif.nnz;i++){
        I_index = vec_I1_i2[K_modif.I_row_indices[i]-offset][1]+offset;
        J_index = vec_I1_i2[K_modif.J_col_indices[i]-offset][1]+offset;
        if (I_index>J_index){
          I_row_indices_p[i]=J_index; J_col_indices_p[i]=I_index;
        }
        else{
          I_row_indices_p[i]=I_index; J_col_indices_p[i]=J_index;
        }
      } 
//
      for (eslocal i = 0; i<K_modif.nnz;i++){
        K_modif.I_row_indices[i] = I_row_indices_p[i];
        K_modif.J_col_indices[i] = J_col_indices_p[i];
      }                            
      K_modif.ConvertToCSRwithSort(0);
    }
//
    for (eslocal i = 0;i<SC_SIZE;i++) fix_dofs[i]=permVec[NONSING_SIZE + i] + offset;
    K_rr.getSubDiagBlockmatrix(K_modif,K_rr,i_start, NONSING_SIZE);
//    K_rr.printMatCSR("K_rr");
    if (CHECK_NONSING!=0){
      cond_of_regular_part = K_rr.MatCondNumb(K_rr,"K_rr",plot_n_first_n_last_eigenvalues);
      printf("cond_of_regular_part=%3.9f\n",cond_of_regular_part);
    }
//
    cnt_iter_check_nonsing++;
  } 

  delete [] I_row_indices_p;
  delete [] J_col_indices_p;

//  for (eslocal i = 0 ; i< SC_SIZE;i++){ printf("%d ",permVec[NONSING_SIZE+i]); } printf("\n");
//
  K_rs.getSubBlockmatrix_rs(K_modif,K_rs,i_start, NONSING_SIZE,j_start,SC_SIZE);
//
  SparseSolver createSchur;
  // TODO Routine Create_SC can provide also factorization of K_rr which can avoid to latter factorization.
  createSchur.ImportMatrix(K_modif);
  createSchur.Create_SC(S,SC_SIZE,false);
  K_modif.Clear();
  createSchur.Clear();
  S.type='S';
  S.ConvertCSRToDense(1);
// EIGENVALUES AND EIGENVECTORS OF SCHUR COMPLEMENT
  char JOBZ = 'V';
  char UPLO = 'U';
  double *W = new double[S.cols]; 
  double *Z = new double[S.cols*S.cols]; 
  MKL_INT info;
  MKL_INT ldz = S.cols;
  info = LAPACKE_dspev (LAPACK_COL_MAJOR, JOBZ, UPLO, S.cols, &(S.dense_values[0]), W, Z, ldz);
  if (info){
    printf("info = %d\n, something wrong with Schur complement in SparseSolver::generalIinverse",info);
  }
// IDENTIFICATIONS OF ZERO EIGENVALUES 
  eslocal defect_A_in;// R_s_cols;
  double ratio; 
  eslocal itMax = TWENTY < S.rows ? TWENTY : S.rows ;
  for (eslocal i = itMax-1; i > 0;i--){
    ratio = fabs(W[i-1]/W[i]);
//    printf("eig[%d],eig[%d]= [%3.15e/%3.15e]\n",i,i-1,W[i],W[i-1]);
    if (ratio < JUMP_IN_EIGENVALUES_ALERTING_SINGULARITY){
      defect_A_in=i;
 //     printf("\tratio = %3.15e, DEFECT = %d\n",ratio,defect_A_in);
    }
  }
// 
  if (get_n_first_and_n_last_eigenvals_from_dense_S!=0){
    printf("eigenvals of %s d{1:%d} and d{%d:%d}\n",
          "S",get_n_first_and_n_last_eigenvals_from_dense_S,
          S.rows-get_n_first_and_n_last_eigenvals_from_dense_S+2,S.rows);
    for (eslocal i = 0 ; i < S.rows; i++){
      if (i < get_n_first_and_n_last_eigenvals_from_dense_S || 
            i > S.rows-get_n_first_and_n_last_eigenvals_from_dense_S){
        printf("%5d:  %3.8e \n",i+1, W[i]);
      }
    }
  }
// --------------- CREATING KERNEL R_s FOR SINGULAR PART (SCHUR COMPLEMENT)
  SparseMatrix R_s;
  R_s.nnz  = defect_A_in*S.rows;
  R_s.dense_values.resize(R_s.nnz);
  R_s.rows = S.rows;
  R_s.cols = defect_A_in;
  R_s.type = 'G';
 // 
  eslocal cnt=0;
  for (eslocal j = 0; j < defect_A_in; j++){
    for (eslocal i = 0; i < R_s.rows; i++){
	    R_s.dense_values[cnt] = Z[j*R_s.rows + i];			
      cnt++;
    }
  }
  R_s.ConvertDenseToCSR(0);
// --------------- CREATING KERNEL R_r FOR NON-SINGULAR PART
	SparseMatrix R_r; 
	R_r.MatMat(K_rs,'N',R_s); 
  K_rs.Clear();
  SparseSolver K_rr_solver; 
  K_rr_solver.ImportMatrix(K_rr);
  K_rr.Clear();
  K_rr_solver.Factorization();
  K_rr_solver.SolveMat_Sparse(R_r); // inv(K_rr)*K_rs*R_s 
  K_rr_solver.Clear();
  
  R_r.ConvertCSRToDense(0);
  R_s.ConvertCSRToDense(0);
// --------------- CREATING WHOLE KERNEL Kplus_R = [ (R_r)^T (R_s)^T ]^T
  Kplus_R.rows = R_r.rows+R_s.rows;
  Kplus_R.cols = R_r.cols;
  Kplus_R.nnz  = Kplus_R.cols*Kplus_R.rows;
  Kplus_R.type = 'G';
	Kplus_R.dense_values.resize(Kplus_R.nnz);			
  cnt=0;
  for (eslocal j = 0; j < Kplus_R.cols; j++){
    for (eslocal i = 0; i < R_r.rows; i++){
      if (diagonalScaling){
        di=K.CSR_V_values[K.CSR_I_row_indices[permVec[i]]-offset];
      }
      Kplus_R.dense_values[j*Kplus_R.rows + permVec[i]] = R_r.dense_values[j*R_r.rows + i]/sqrt(di);			
      cnt++;
    }
    for (eslocal i = 0; i < R_s.rows; i++){
      if (diagonalScaling){
        di=K.CSR_V_values[K.CSR_I_row_indices[permVec[i+R_r.rows]]-offset];
      }
	    Kplus_R.dense_values[j*Kplus_R.rows + permVec[i+R_r.rows]] =-R_s.dense_values[j*R_s.rows + i]/sqrt(di);
      cnt++;
    }
  }
//
  Kplus_R.GramSchmidtOrtho();
  SEQ_VECTOR <eslocal > null_pivots;
  Kplus_R.getNullPivots(null_pivots);

//  printf("null pivots, eslocal ");
//  for (eslocal i = 0;i<null_pivots.size();i++)
//    printf("%d ",null_pivots[i]);
//  printf("\n");

//
  double * AR =  new double [K.rows];
  double norm_AR_row,norm_AR = 0.0;
//  printf("||A*Kplus_R[:,i]|| ...   \n");
  for (eslocal i = 0;i<Kplus_R.cols;i++){
    memset(AR,0,Kplus_R.rows * sizeof(double));
  	K.spmv_( K,&(Kplus_R.dense_values[i*Kplus_R.rows]),AR);
    norm_AR_row=0.0;
    for (eslocal j = 0; j < Kplus_R.rows;j++){
      norm_AR_row+=AR[j]*AR[j];
    }
 //   printf("%3.3e  ",sqrt(norm_AR_row));
    norm_AR+=norm_AR_row;
  }
  delete [] AR;
  norm_AR=sqrt(norm_AR);
  printf("\n||A*Kplus_R|| = %3.9e \n",norm_AR);

  Kplus_R.ConvertDenseToCSR(1);
//  Kplus_R.printMatCSR("Kplus_R");
//  K.printMatCSR("K");
// REGULARIZATION OF MATRIX K

  double rho = K.GetMaxOfDiagonalOfSymmetricMatrix();
  rho = 1.0 * rho;
//
  if (diagonalRegularization){
    for (eslocal i = 0; i < null_pivots.size(); i++){
      K.CSR_V_values[K.CSR_I_row_indices[null_pivots[i]-offset]-offset]+=rho;
    }
  }
  else{
    SparseMatrix N;
    if (use_null_pivots_or_s_set){
      N.CreateMatFromRowsFromMatrix( Kplus_R, null_pivots);
    }
    else
    {
//      printf("fix_dofs...\n");
//      for (eslocal i = 0;i<fix_dofs.size();i++)
//        printf("%d ",fix_dofs[i]);
//      printf("\n");

      N.CreateMatFromRowsFromMatrix( Kplus_R, fix_dofs);
    }
  //null_pivots
    SparseMatrix Nt;
    N.MatTranspose( Nt );
    SparseMatrix NtN_Mat;
    NtN_Mat.MatMat( Nt,'N',N );
    NtN_Mat.MatTranspose();
    NtN_Mat.RemoveLower();
    SparseSolver NtN;
    NtN.ImportMatrix(NtN_Mat);
    NtN_Mat.Clear();
    NtN.Factorization();
    NtN.SolveMat_Sparse(Nt);
    NtN.Clear();
    NtN_Mat.MatMat(N,'N',Nt);
    NtN_Mat.RemoveLower();
    //K.MatCondNumb(K,"K_singular",plot_n_first_n_last_eigenvalues);
    K.MatAddInPlace (NtN_Mat,'N', rho);
  }
//  K.printMatCSR("K_regularized");
//  K.MatCondNumb(K,"K_regularized",plot_n_first_n_last_eigenvalues);


  printf(" =============================================================\n\n");

  delete [] W;
  delete [] Z;
}



// **** END - DOMAIN CLASS *******************************************
// *******************************************************************
