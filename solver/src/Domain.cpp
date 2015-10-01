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

void Domain::LoadDomain(string directory_path, int USE_HFETI, int use_dynamic_1_no_dynamic_0) {

	USE_DYNAMIC = use_dynamic_1_no_dynamic_0; 
	
	// Filenames of matrices to load 
	string filename_B0; 
	string filename_B1; 
	string filename_Kplus;

	string filename_map_vector_e0; 
	string filename_map_vector; 
	string filename_f; 
	string filename_coordinates; 
	string filename_fix_nodes; 
	string filename_local2global0; 
	string filename_nodeMulti; 
	// ****************

	filename_B0      = "B0"; 
	filename_B1      = "B1"; 
	filename_Kplus   = "K"; 

	filename_map_vector_e0 = "mapVec_e0"; 
	filename_map_vector	   = "mapVec";  
	filename_f			   = "f";
	filename_coordinates   = "COORDINATES";
	filename_fix_nodes	   = "fixNodes";
	filename_local2global0 = "local2global0";
	filename_nodeMulti	   = "nodeMulti"; 
	// ****************

	char number[8];
	sprintf(number, "%06d", domain_global_index);
	string path = string(directory_path) +  "DOMAIN_" + number + "_"; 


	// *** Load Matrix B0 - gluing matrix for domains inside the cluster **************
	if (USE_HFETI == 1)
		B0.LoadMatrixBin(string(path) + string(filename_B0), 'G', 0);
	// *** END - Load Matrix B0 - gluing matrix for domains inside the cluster ********


	// *** Load Matrix B1 - gluing matrix for FETI ************************************
	B1.		LoadMatrixBin(string(path) + string(filename_B1), 'G', 0); // this should not be necessary 
	B1_comp.LoadMatrixBin(string(path) + string(filename_B1), 'G', 0);
	// *** END - Load Matrix B1 - gluing matrix for FETI ******************************


	// *** Load the coordinates *******************************************************
	int tres; 
	tres = LoadCoordinates( string(path) + string(filename_coordinates) ); 
	if (tres == -1) {
		// !!! POZOR - works only for CUBE !!!!
		tres = LoadCoordinates( string(directory_path) + string(filename_coordinates) ); 

		SEQ_VECTOR <SEQ_VECTOR <double> > correct_coordinates; 
		LoadBinVecVec(correct_coordinates, string(directory_path) + string("Correct_COORDINATES"));

		for (int i = 0; i < coordinates.size(); i++) {
			for (int j = 0; j < coordinates[i].size(); j++) {
				coordinates[i][j] = coordinates[i][j] + correct_coordinates[domain_global_index - 1][j]; // domains numbering is from 1 
			}
		}

	}
	// *** END - Load the coordinates *************************************************


	// *** Create the FEM assembler object ********************************************
	string domain_name = string("DOMAIN_") + number + "_"; 
	FEM_Assembler K_asm(string(directory_path), domain_name, domain_global_index, USE_DYNAMIC ); 
	// *** END - Create the FEM assembler object **************************************


	if (USE_DYNAMIC == 0) {

		// *** Create the stifness matrix *********************************************
		K_asm.assemble_matrix(K, f);	//LoadBinVectorDouble(f, string(path) + string(filename_f));
		// *** END - Create the stifness matrix ***************************************
		
		// *** Prepare LUMPED preconditioner ******************************************
		Prec = K;		//Prec.MatAddInPlace(K,'N', 1);
		Prec.type = 'G'; 
		// *** END - Prepare LUMPED preconditioner ************************************

		K.RemoveLower();

		// *** Regularize K ***********************************************************
		// Loads fix points for each subdomain - if not present, loads global fixNodes file - valid only for cube. 
		tres = LoadBinVectorInt(fix_nodes, string(path) + string(filename_fix_nodes));
		if (tres == -1) {
			// !!! POZOR - works only for CUBE !!!!
			tres = LoadBinVectorInt(fix_nodes, string(directory_path) + string(filename_fix_nodes)); 
		}
		std::sort(fix_nodes.begin(),fix_nodes.end());
		K_regularization();
		// *** END - Regularize K *****************************************************

		// *** matrix - R *************************************************************
		CreateKplus_R();	//Kplus_R.LoadMatrixBin(string(path) + string("R") , 'G'); 
		// *** END - matrix - R *******************************************************

		domain_prim_size = Kplus.cols; 

	} else { // USE_DYNAMIC == 0

		//LoadBinVectorDouble(f, string(path) + string(filename_f));
		K_asm.assemble_matrix(K, M, f); //tmpf); 

		double time_const = 1.0 / ( dynamic_beta * dynamic_timestep * dynamic_timestep);  
		K.MatAddInPlace(M,'N', time_const); 
		M.type='G';

		Prec = K; //Prec.MatAddInPlace(K,'N', 1);
		Prec.type='G'; 

		K.RemoveLower();
		Kplus.ImportMatrix(K);
		Kplus.Factorization();

		KplusF.ImportMatrix(K);

		domain_prim_size = Kplus.cols; 
	}
	// *** END - matrix K - load and regularize ***********************************


	// POZOR - pozor - These vectors are needed only for findSolution function 
	LoadBinVectorInt(map_vector_local2global0, string(path) + string( filename_local2global0 )); 
	LoadBinVectorInt(nodeMulti, string(path) + string( filename_nodeMulti )); 
	LoadBinVectorInt(number_of_nodes_in_global0, string(directory_path) + string("nNodes"));
	// END - These vectors are needed only for findSolution function 


}


void Domain::SetDomain(int USE_HFETI, int use_dynamic_1_no_dynamic_0) {

	USE_DYNAMIC = use_dynamic_1_no_dynamic_0; 


	// *** Load Matrix B0 - gluing matrix for domains inside the cluster **************
	//if (USE_HFETI == 1)
	//	 B0.LoadMatrixBin(string(path) + string(filename_B0), 'G', 0);
	// *** END - Load Matrix B0 - gluing matrix for domains inside the cluster ********


	// *** Load Matrix B1 - gluing matrix for FETI ************************************
	// B1.		LoadMatrixBin(string(path) + string(filename_B1), 'G', 0); // this should not be necessary 
	// B1_comp.LoadMatrixBin(string(path) + string(filename_B1), 'G', 0);
	// B1_comp = B1; 
	// *** END - Load Matrix B1 - gluing matrix for FETI ******************************


	// *** Load the coordinates *******************************************************
	//int tres; 
	//tres = LoadCoordinates( string(path) + string(filename_coordinates) ); 
	//if (tres == -1) {
	//	// !!! POZOR - works only for CUBE !!!!
	//	tres = LoadCoordinates( string(directory_path) + string(filename_coordinates) ); 

	//	vector <vector <double> > correct_coordinates; 
	//	LoadBinVecVec(correct_coordinates, string(directory_path) + string("Correct_COORDINATES"));

	//	for (int i = 0; i < coordinates.size(); i++) {
	//		for (int j = 0; j < coordinates[i].size(); j++) {
	//			coordinates[i][j] = coordinates[i][j] + correct_coordinates[domain_global_index - 1][j]; // domains numbering is from 1 
	//		}
	//	}

	//}
	// *** END - Load the coordinates *************************************************


	// *** Create the FEM assembler object ********************************************
	//string domain_name = string("DOMAIN_") + number + "_"; 
	//FEM_Assembler K_asm(string(directory_path), domain_name, domain_global_index, USE_DYNAMIC ); 
	// *** END - Create the FEM assembler object **************************************


	//if (USE_DYNAMIC == 0) {

		// *** Create the stifness matrix *********************************************
		//K_asm.assemble_matrix(K, f);	
		// *** END - Create the stifness matrix ***************************************

		// *** Prepare LUMPED preconditioner ******************************************
		//Prec = K;		
		//Prec.type = 'S'; 
		// *** END - Prepare LUMPED preconditioner ************************************

		//K.RemoveLower();

		// *** Regularize K ***********************************************************
		// Loads fix points for each subdomain - if not present, loads global fixNodes file - valid only for cube. 
		//tres = LoadBinVectorInt(fix_nodes, string(path) + string(filename_fix_nodes));
		//if (tres == -1) {
		//	// !!! POZOR - works only for CUBE !!!!
		//	tres = LoadBinVectorInt(fix_nodes, string(directory_path) + string(filename_fix_nodes)); 
		//}
		//std::sort(fix_nodes.begin(),fix_nodes.end());
		//K_regularization();
		K_regularizationFromR( );
		// *** END - Regularize K *****************************************************

		// *** matrix - R *************************************************************
		//CreateKplus_R();	//Kplus_R.LoadMatrixBin(string(path) + string("R") , 'G'); 
		// *** END - matrix - R *******************************************************

		domain_prim_size = Kplus.cols; 

	//} else { // USE_DYNAMIC == 0

	//	//LoadBinVectorDouble(f, string(path) + string(filename_f));
	//	K_asm.assemble_matrix(K, M, f); //tmpf); 

	//	double time_const = 1.0 / ( dynamic_beta * dynamic_timestep * dynamic_timestep);  
	//	K.MatAddInPlace(M,'N', time_const); 
	//	M.type='G';

	//	Prec = K; //Prec.MatAddInPlace(K,'N', 1);
	//	Prec.type='G'; 

	//	K.RemoveLower();
	//	Kplus.ImportMatrix(K);
	//	Kplus.Factorization();

	//	KplusF.ImportMatrix(K);

	//	domain_prim_size = Kplus.cols; 
	//}
	// *** END - matrix K - load and regularize ***********************************


	// POZOR - pozor - These vectors are needed only for findSolution function 
	//LoadBinVectorInt(map_vector_local2global0, string(path) + string( filename_local2global0 )); 
	//LoadBinVectorInt(nodeMulti, string(path) + string( filename_nodeMulti )); 
	//LoadBinVectorInt(number_of_nodes_in_global0, string(directory_path) + string("nNodes"));
	// END - These vectors are needed only for findSolution function 


}


int  Domain::LoadCoordinates(string filename) {

	ifstream in (filename.c_str(), std::ios::binary);
	char delim = ';'; 
	string line, field;

	if ( in.is_open() ) {

		// Get parameters 
		getline(in,line);
		stringstream paramss(line);

		getline(paramss,field,delim); 
		int rows = atoi(field.c_str());		// get num of rows 

		getline(paramss,field,delim); 
		int cols = atoi(field.c_str());		// get num of columns 

		// Get data 

		coordinates.reserve(rows);

		for (int r_i = 0; r_i < rows; r_i++) {
			SEQ_VECTOR <double> tmp_coor (cols); 
			in.read((char*) &tmp_coor [0], cols*sizeof(double));
			//coordinates[r_i] = tmp_coor; 
			coordinates.push_back(tmp_coor); 
		}

		in.close(); 

		return 0; 
	} else {
		cout << "File " << filename << " not found ! " << endl; 
		return -1; 
	}
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
			N.Clear();
			Nt.MatTranspose(N);
			NtN_Mat.MatMat(N,'N',Nt);
			NtN_Mat.RemoveLower();

			double ro = K.GetMaxOfDiagonalOfSymmetricMatrix();
			ro = 1.0 * ro;

			K.    MatAddInPlace (NtN_Mat,'N', ro);
    	}


    	if (DOFS_PER_NODE == 1) {
    		double ro = K.GetMaxOfDiagonalOfSymmetricMatrix();
    		K.CSR_V_values[0]+=ro;
    	}

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
