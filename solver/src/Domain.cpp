#include "Domain.h"
#include <cmath>

// *******************************************************************
// **** DOMAIN CLASS ************************************************


Domain::Domain(){

}

Domain::Domain(int domain_index, int use_dynamic_1_no_dynamic_0){
	domain_global_index = domain_index; 
	USE_DYNAMIC = use_dynamic_1_no_dynamic_0; 
}

void Domain::SetDynamicParameters(double set_dynamic_timestep, double set_dynamic_beta, double set_dynamic_gama) {
	dynamic_timestep = set_dynamic_timestep;
	dynamic_beta     = set_dynamic_beta;
	dynamic_gama     = set_dynamic_gama;
}

void Domain::SetDomain(int USE_HFETI, int use_dynamic_1_no_dynamic_0) {

	USE_DYNAMIC = use_dynamic_1_no_dynamic_0; 
	K_regularizationFromR( );
	domain_prim_size = Kplus.cols;

}

void Domain::CreateKplus_R ( ) 
{
	if (DOFS_PER_NODE == 3) {

		int elem_index = 0;

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

		for (int elem_index = 1; elem_index < coordinates.size(); elem_index++) {
			Kplus_R.dense_values[elem_index] = nsqrt;
		}

		Kplus_R.rows = coordinates.size();
		Kplus_R.cols = 1;
		Kplus_R.nnz  = coordinates.size();
		Kplus_R.type = 'G';

		Kplus_R.ConvertDenseToCSR(1);
	}
}

void Domain::multKplusLocal(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, int x_in_vector_start_index, int y_out_vector_start_index) {
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
	int row_fill = 1; 
	N.CSR_I_row_indices.reserve(K.rows+1);

	int old_index  = 0;
	int next_index = 0; 

	SEQ_VECTOR <double> dense_values;
	double dense_pattern[] = 
	   {1, 0, 0,    0, 8, 8,
		0, 1, 0,    8, 0, 8,
		0, 0, 1,    8, 8, 0}; // 8 is ust a random number 
	int dense_pattern_size = 18; 

	dense_values.resize(fix_nodes.size() * dense_pattern_size);
	for (int fpi = 0; fpi < fix_nodes.size(); fpi++) {
		int elem_index = fix_nodes[fpi];
		dense_pattern[4]  = - coordinates[elem_index][2];
		dense_pattern[5]  =   coordinates[elem_index][1];
		dense_pattern[9]  =   coordinates[elem_index][2];
		dense_pattern[11] = - coordinates[elem_index][0];
		dense_pattern[15] = - coordinates[elem_index][1];
		dense_pattern[16] =   coordinates[elem_index][0];
		copy(dense_pattern, dense_pattern+dense_pattern_size, dense_values.begin() + (fpi*dense_pattern_size));
		//dense_values.assign(dense_pattern + (fpi*dense_pattern_size), dense_pattern + ((fpi+1)*dense_pattern_size));
	}

	for (int fpi = 0; fpi < fix_nodes.size(); fpi++) {

		old_index  = next_index;
		next_index = 3 * fix_nodes[fpi] - 2;
		N.CSR_I_row_indices.resize( next_index );
		fill(N.CSR_I_row_indices.begin() + old_index, N.CSR_I_row_indices.begin() + next_index, row_fill);
		N.rows = next_index; 

		int elem_index = fix_nodes[fpi];
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
	N.Clear();
	Nt.MatTranspose(N);
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

void Domain::get_kernel_from_A(){
//      
// Routine calculates kernel Kplus_R of K satisfied euqality
//      K * Kplus_R = O,
// where O is zero matrix,
// and makes the matrix K non-singular (A_reg) utilizing spectral conditions
// Schur complement. Then 
//      || K - K*inv(A_reg)*K || = 0.0
//
//
// rev. 2015-10-11 (A.M.)
  printf("\n\n");
  printf(" ###################################################################\n");
  printf(" #                 GET KERNEL OF A AND NULL PIVOTS                 #\n");
  printf(" ###################################################################\n");
//
// 
// TODO FOR MORE GENERAL CASE MATRIX K SHOULD BE RANDOMLY PERMUTED. 
  int SC_SIZE = 50;         // size of S to detect singular part  
  int TWENTY  = 20;         // constant set-up to DEFECT <= TWENTY <= SC_SIZE 
  //TODO if K.rows<=SC_SIZE, use directlyi input K instead of S
  SparseMatrix S;
  SparseMatrix A_rr;
  SparseMatrix A_rs;
  SparseMatrix A_modif = K;
  bool diagonalScaling = true;
  bool permutVectorActive = true;
  bool diagonalRegularization=false;
  bool useNullPivotsOr_s_set=false;

  // diagonal scaling of A_modif:
  // A_modif[i,j] = A_modif[i,j]/sqrt(A_modif[i,i]*A_modif[j,j]);
  double di=1,dj=1;
  for (int i = 0;i<A_modif.rows;i++){
    if (diagonalScaling) {
      di=A_modif.CSR_V_values[A_modif.CSR_I_row_indices[i]-1];
    }
    for (int j = A_modif.CSR_I_row_indices[i];j<A_modif.CSR_I_row_indices[i+1];j++){
      if (diagonalScaling) {
        dj=A_modif.CSR_V_values[A_modif.CSR_I_row_indices[A_modif.CSR_J_col_indices[j-1]-1]-1];
      }
      A_modif.CSR_V_values[j-1]/=sqrt(di*dj);
    }
  }

//  A_modif.printMatCSR("A_modif");


  vector<int> permVec;
  int i_start = 0;
  int NONSING_SIZE = K.rows - SC_SIZE - i_start;
  int j_start = NONSING_SIZE;

  // set row permVec = {0,1,2,3,4,...,K.rows};
  for (int i=0; i<K.rows; ++i) permVec.push_back(i); // 0 1 2 A.rows
  // random permutation
  
  if (permutVectorActive) random_shuffle ( permVec.begin(), permVec.end() );

  //      r = permVec[0:NONSING_SIZE-1]     (singular DOFs)
  //      s = permVec[NONSING_SIZE:end-1]   (non-singular DOFs)
  sort(permVec.begin(),permVec.begin()+NONSING_SIZE);
  sort(permVec.begin()+NONSING_SIZE,permVec.end());

  std::vector<std::vector<int>> vec_I1_i2(K.rows, std::vector<int>(2, 1));

  for (int i = 0; i < K.rows;i++){
    vec_I1_i2[i][0]=permVec[i];
    vec_I1_i2[i][1]=i; // position to create revers permutation
  }

  std::sort(vec_I1_i2.begin(), vec_I1_i2.end(),
            [](const std::vector<int>& a, const std::vector<int>& b) {
    return a[0] < b[0];
  });

  A_modif.ConvertToCOO(0);
  int *I_row_indices_p = new int[A_modif.nnz] ;
  int *J_col_indices_p = new int[A_modif.nnz] ;
//	V_values; J_col_indices; I_row_indices;	    
    int I_index,J_index;
    int offset = K.CSR_I_row_indices[0] ? 1 : 0;
    for (int i = 0;i<A_modif.nnz;i++){
      I_index = vec_I1_i2[A_modif.I_row_indices[i]-offset][1]+offset;
      J_index = vec_I1_i2[A_modif.J_col_indices[i]-offset][1]+offset;
      if (I_index>J_index){
        I_row_indices_p[i]=J_index;
        J_col_indices_p[i]=I_index;
      }
      else{
        I_row_indices_p[i]=I_index;
        J_col_indices_p[i]=J_index;
      }
   } 
//
  for (int i = 0; i<A_modif.nnz;i++){
    A_modif.I_row_indices[i] = I_row_indices_p[i];
    A_modif.J_col_indices[i] = J_col_indices_p[i];
  }                            
  A_modif.ConvertToCSRwithSort(0);
  delete [] I_row_indices_p;
  delete [] J_col_indices_p;

//
  vector <int> fix_dofs;
  for (int i = 0;i<SC_SIZE;i++) fix_dofs.push_back(permVec[NONSING_SIZE + i] + offset); 
  A_rr.getSubDiagBlockmatrix(A_modif,A_rr,i_start, NONSING_SIZE);
  A_rr.MatCondNumb(A_rr,"A_rr");
  A_rs.getSubBlockmatrix_rs(A_modif,A_rs,i_start, NONSING_SIZE,j_start,SC_SIZE);
//
  SparseSolver createSchur;
  createSchur.ImportMatrix(A_modif);
  createSchur.Create_SC(S,SC_SIZE,false);
  createSchur.Clear();
//  S.printMatCSR("Schur_complement");
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
  int defect_A_in;// R_s_cols;
  double ratio; 
  int itMax = TWENTY < S.rows ? TWENTY : S.rows ;
  for (int i = itMax-1; i > 0;i--){ 
    ratio = fabs(W[i-1]/W[i]);
//    printf("eig[%d],eig[%d]= [%3.15e/%3.15e]\n",i,i-1,W[i],W[i-1]);
    if (ratio < 1e-5){
      defect_A_in=i;
      printf("\tratio = %3.15e, DEFECT = %d\n",ratio,defect_A_in);
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
  int cnt=0;
  for (int j = 0; j < defect_A_in; j++){
    for (int i = 0; i < R_s.rows; i++){
	    R_s.dense_values[cnt] = Z[j*R_s.rows + i];			
      cnt++;
    }
  }
  R_s.ConvertDenseToCSR(0);
//  R_s.printMatCSR("R_s");
// --------------- CREATING KERNEL R_r FOR NON-SINGULAR PART
	SparseMatrix R_r; 
	R_r.MatMat(A_rs,'N',R_s); 
  //R_r.printMatCSR("A_rs*R_s");
  SparseSolver A_rr_solver; 
  A_rr_solver.ImportMatrix(A_rr);
  A_rr.Clear();
  A_rr_solver.Factorization();
  A_rr_solver.SolveMat_Sparse(R_r); // inv(A_rr)*A_rs*R_s 
  
//  R_r.printMatCSR("R_r");
  R_r.ConvertCSRToDense(0);
  R_s.ConvertCSRToDense(0);
// --------------- CREATING WHOLE KERNEL Kplus_R = [ (R_r)^T (R_s)^T ]^T
  Kplus_R.rows = R_r.rows+R_s.rows;
  Kplus_R.cols = R_r.cols;
  Kplus_R.nnz  = Kplus_R.cols*Kplus_R.rows;
  Kplus_R.type = 'G';
	Kplus_R.dense_values.resize(Kplus_R.nnz);			
  cnt=0;
  for (int j = 0; j < Kplus_R.cols; j++){
    for (int i = 0; i < R_r.rows; i++){
      if (diagonalScaling){
        di=K.CSR_V_values[K.CSR_I_row_indices[permVec[i]]-1];
      }
      Kplus_R.dense_values[j*Kplus_R.rows + permVec[i]] = R_r.dense_values[j*R_r.rows + i]/sqrt(di);			
      cnt++;
    }
    for (int i = 0; i < R_s.rows; i++){
      if (diagonalScaling){
        di=K.CSR_V_values[K.CSR_I_row_indices[permVec[i+R_r.rows]]-1];
      }
	    Kplus_R.dense_values[j*Kplus_R.rows + permVec[i+R_r.rows]] =-R_s.dense_values[j*R_s.rows + i]/sqrt(di);
      cnt++;
    }
  }
//
  Kplus_R.GramSchmidtOrtho();
  SEQ_VECTOR <int> null_pivots;
  Kplus_R.getNullPivots(null_pivots);

  printf("null_pivots ...\n");
  for (int i = 0;i<null_pivots.size();i++)
    printf("%d ",null_pivots[i]);
  printf("\n");

//
  double * AR =  new double [K.rows];
  double norm_AR_row,norm_AR = 0.0;
  printf("||A*Kplus_R[:,i]|| ...   \n");
  for (int i = 0;i<Kplus_R.cols;i++){
    memset(AR,0,Kplus_R.rows * sizeof(double));
  	K.spmv_( K,&(Kplus_R.dense_values[i*Kplus_R.rows]),AR);
    norm_AR_row=0.0;
    for (int j = 0; j < Kplus_R.rows;j++){
      norm_AR_row+=AR[j]*AR[j];
    }
    printf("%3.3e  ",sqrt(norm_AR_row));
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
//  K.printMatCSR("A_singular");
  K.MatCondNumb(K,"A_singular");
  if (diagonalRegularization){
    for (int i = 0; i < null_pivots.size(); i++){
      K.CSR_V_values[K.CSR_I_row_indices[null_pivots[i]-1]-1]+=rho;
    }
  }
  else{
    SparseMatrix N;
    if (useNullPivotsOr_s_set){
      N.CreateMatFromRowsFromMatrix( Kplus_R, null_pivots);
    }
    else
    {
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
    //K.MatCondNumb(K,"A_singular");
    K.MatAddInPlace (NtN_Mat,'N', rho);
  }
//  K.printMatCSR("A_regularized");
  K.MatCondNumb(K,"A_regularized");


  printf(" =============================================================\n\n");

  delete [] W;
  delete [] Z;
  //delete [] perm;

}




void Domain::K_regularizationFromR ( ) {

    if (USE_DYNAMIC == 0) {

//    	if (DOFS_PER_NODE == 3) {
//
//			SparseMatrix N;
//			N.printMatCSR(N,"Ntt");
//
//			N.CreateMatFromRowsFromMatrix( Kplus_R, fix_dofs);
//
//			SparseMatrix Nt;
//			N.MatTranspose( Nt );
//
//			SparseMatrix NtN_Mat;
//			NtN_Mat.MatMat( Nt,'N',N );
//        NtN_Mat.MatTranspose();
//        NtN_Mat.RemoveLower();
//
//        SparseSolver NtN;
//        NtN.ImportMatrix(NtN_Mat);
//        NtN_Mat.Clear();
//
//        NtN.Factorization();
//        NtN.SolveMat_Sparse(Nt);
//        NtN.Clear();
//
//        //NtN = Nt*N
//        N.Clear();
//        Nt.MatTranspose(N);
//        NtN_Mat.MatMat(N,'N',Nt);
//        NtN_Mat.RemoveLower();
//
//        double ro = K.GetMaxOfDiagonalOfSymmetricMatrix();
//        ro = 1.0 * ro;
//
//        K.MatAddInPlace (NtN_Mat,'N', ro);
//    	}


//     	if (DOFS_PER_NODE == 1) {
//        SparseMatrix forGInv;
//        SparseMatrix Kplus_R; // automatically detected kernl of K
          
//        K.CSR_V_values[0]*=1.1;
//   	  }
    }



    Kplus.ImportMatrix  (K);
    //Kplus.Factorization ();

	// POZOR - jen pro Kinv
//	if (USE_KINV == 1 || USE_HFETI == 1)
//		KplusF.ImportMatrix (K);

	//K.    MatAddInPlace (NtN_Mat,'N', -ro);
	//K.Clear(); //TODO: potom se nekde musi K smazat
	
} 

// **** END - DOMAIN CLASS *******************************************
// *******************************************************************
