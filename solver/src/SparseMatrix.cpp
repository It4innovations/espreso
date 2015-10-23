#include "SparseMatrix.h"


void SpyText (SparseMatrix & A) {
	
	SEQ_VECTOR<char> tmp (60,'-');  

	for( std::SEQ_VECTOR<char>::const_iterator i = tmp.begin(); i != tmp.end(); ++i)
		std::cout << *i << ' ';
    
	cout << endl; 

	int rows_coef = 1 + A.rows / 60;
	int cols_coef = 1 + A.cols / 60;

	int col_index = 0; 
	int row_index = 0; 
	for (int r = 0; r < A.rows; r = r + rows_coef) {
		int row_length = 0; 
		if (( r + rows_coef) < A.rows)
			row_length = A.CSR_I_row_indices[r+rows_coef] - A.CSR_I_row_indices[r];
		else
 			row_length = A.CSR_I_row_indices[A.rows] - A.CSR_I_row_indices[r];

			SEQ_VECTOR<char> tmp (60,' ');  
			SEQ_VECTOR<int> tmp_c (60,0);  
		for (int c = 0; c < row_length; c++) {
			if (A.CSR_V_values[col_index] != 0.0) {
				tmp_c[A.CSR_J_col_indices[col_index] / cols_coef]++;  
			} else {
				if (tmp_c[A.CSR_J_col_indices[col_index] / cols_coef] == 0)
					tmp_c[A.CSR_J_col_indices[col_index] / cols_coef] = -1; 
			}
			col_index++;
		}

		for (int c = 0; c < tmp_c.size(); c++) {
			if (tmp_c[c] > 0) {
				tmp[c] = '0' + tmp_c[c] / (cols_coef * 26); 
				if (tmp[c] == '0') tmp[c] = '.'; 
			} else {
				if (tmp_c[c] == -1)
					tmp[c] = 'x'; 
			}
		}

		cout << "|";
		for( std::SEQ_VECTOR<char>::const_iterator i = tmp.begin(); i != tmp.end(); ++i)
			std::cout << *i << ' ';

		cout << "|" << endl; 
	}

	//SEQ_VECTOR<char> tmp (60,'-');  
	
	for( std::SEQ_VECTOR<char>::const_iterator i = tmp.begin(); i != tmp.end(); ++i)
		std::cout << *i << ' ';

	cout << endl << endl; 
}

void sortMatrixInCOO(SparseMatrix & Matrix) 
{
	q_sort(Matrix, 0, Matrix.nnz - 1);

	//int k = 0;
	///*  Remove duplicates */
	//for( int i = 1; i < Matrix.nnz; i++) {
	//	if ( (Matrix.I_row_indices[k] != Matrix.I_row_indices[i]) || (Matrix.J_col_indices[k] != Matrix.J_col_indices[i]) ){
	//		k++;
	//		Matrix.I_row_indices[k] = Matrix.I_row_indices[i];
	//		Matrix.J_col_indices[k] = Matrix.J_col_indices[i];
	//		Matrix.V_values[k] = Matrix.V_values[i];
	//	} else {
	//		Matrix.V_values[k] += Matrix.V_values[i];
	//	}
	//}
	//Matrix.nnz = k+1;
	//Matrix.I_row_indices.resize(k+1);
	//Matrix.J_col_indices.resize(k+1);	
	//Matrix.V_values.resize(k+1);
	 
}

void SparseMatrix::sortInCOO() 
{
	q_sort_in(I_row_indices, J_col_indices,V_values, 0, nnz - 1);

	//int k = 0;
	///*  Remove duplicates */
	//for( int i = 1; i < Matrix.nnz; i++) {
	//	if ( (Matrix.I_row_indices[k] != Matrix.I_row_indices[i]) || (Matrix.J_col_indices[k] != Matrix.J_col_indices[i]) ){
	//		k++;
	//		Matrix.I_row_indices[k] = Matrix.I_row_indices[i];
	//		Matrix.J_col_indices[k] = Matrix.J_col_indices[i];
	//		Matrix.V_values[k] = Matrix.V_values[i];
	//	} else {
	//		Matrix.V_values[k] += Matrix.V_values[i];
	//	}
	//}
	//Matrix.nnz = k+1;
	//Matrix.I_row_indices.resize(k+1);
	//Matrix.J_col_indices.resize(k+1);	
	//Matrix.V_values.resize(k+1);

}


#define comp_x(_a,_b,_x,_p) (_a[_x]!=_a[_p] ? _a[_x]-_a[_p] : _b[_x]-_b[_p])

static void q_sort_in(vector <int>    & I_row_indices, 
					  vector <int>    & J_col_indices, 
					  vector <double> & V_values,  
					  int lo, int hi ) {
	int h, l, p, p1, p2, t;
	double p3, td;

	if(lo>=hi)
		return;

	l = lo;
	h = hi;
	p = hi;

	p1 = J_col_indices[p];
	p2 = I_row_indices[p];
	p3 = V_values[p];

	do {
		while ((l < h) && (comp_x(I_row_indices,J_col_indices,l,p)<=0))
			l = l+1;
		while ((h > l) && (comp_x(I_row_indices,J_col_indices,h,p)>=0))
			h = h-1;
		if (l < h) {
			t = J_col_indices[l];
			J_col_indices[l] = J_col_indices[h];
			J_col_indices[h] = t;

			t = I_row_indices[l];
			I_row_indices[l] = I_row_indices[h];
			I_row_indices[h] = t;

			td = V_values[l];
			V_values[l] = V_values[h];
			V_values[h] = td;

		}
	} while (l < h);

	J_col_indices[p] = J_col_indices[l];
	J_col_indices[l] = p1;

	I_row_indices[p] = I_row_indices[l];
	I_row_indices[l] = p2;

	V_values[p] = V_values[l];
	V_values[l] = p3;

	/* Sort smaller array first for less stack usage */
	if (l-lo<hi-l) {
		q_sort_in( I_row_indices,J_col_indices,V_values, lo, l-1 );
		q_sort_in( I_row_indices,J_col_indices,V_values, l+1, hi );
	} else {
		q_sort_in( I_row_indices,J_col_indices,V_values, l+1, hi );
		q_sort_in( I_row_indices,J_col_indices,V_values, lo, l-1 );
	}
}


static void q_sort(SparseMatrix & Matrix, int lo, int hi ) {
	int h, l, p, p1, p2, t;
	double p3, td;

	if(lo>=hi)
		return;

	l = lo;
	h = hi;
	p = hi;

	p1 = Matrix.J_col_indices[p];
	p2 = Matrix.I_row_indices[p];
	p3 = Matrix.V_values[p];

	do {
		while ((l < h) && (comp_x(Matrix.I_row_indices,Matrix.J_col_indices,l,p)<=0))
			l = l+1;
		while ((h > l) && (comp_x(Matrix.I_row_indices,Matrix.J_col_indices,h,p)>=0))
			h = h-1;
		if (l < h) {
			t = Matrix.J_col_indices[l];
			Matrix.J_col_indices[l] = Matrix.J_col_indices[h];
			Matrix.J_col_indices[h] = t;

			t = Matrix.I_row_indices[l];
			Matrix.I_row_indices[l] = Matrix.I_row_indices[h];
			Matrix.I_row_indices[h] = t;

			td = Matrix.V_values[l];
			Matrix.V_values[l] = Matrix.V_values[h];
			Matrix.V_values[h] = td;

		}
	} while (l < h);

	Matrix.J_col_indices[p] = Matrix.J_col_indices[l];
	Matrix.J_col_indices[l] = p1;

	Matrix.I_row_indices[p] = Matrix.I_row_indices[l];
	Matrix.I_row_indices[l] = p2;

	Matrix.V_values[p] = Matrix.V_values[l];
	Matrix.V_values[l] = p3;

	/* Sort smaller array first for less stack usage */
	if (l-lo<hi-l) {
		q_sort( Matrix, lo, l-1 );
		q_sort( Matrix, l+1, hi );
	} else {
		q_sort( Matrix, l+1, hi );
		q_sort( Matrix, lo, l-1 );
	}
}


SparseMatrix::SparseMatrix() {
	
	nnz  = 0;
	cols = 0; 
	rows = 0; 
	type = 0; 

	d_dense_values = NULL; 
	d_x_in		   = NULL;
	d_y_out		   = NULL; 

	d_dense_values_fl = NULL; 
	d_x_in_fl		  = NULL;
	d_y_out_fl		  = NULL; 

#ifdef CUDA
	handle		    = NULL;
	stream          = NULL; 
#endif	

}

SparseMatrix::SparseMatrix(char matrix_type_G_for_general_S_for_symmetric, string filename) {
	SparseMatrix::LoadMatrixBin(filename, matrix_type_G_for_general_S_for_symmetric);
}

SparseMatrix::SparseMatrix( const SparseMatrix &A_in) {
	
	rows = A_in.rows; 
	cols = A_in.cols;
	nnz  = A_in.nnz;
	type = A_in.type;

	// Sparse COO data 
	I_row_indices = A_in.I_row_indices;
	J_col_indices = A_in.CSR_J_col_indices;
	V_values	  = A_in.V_values;

	// Sparse CSR data 
	CSR_I_row_indices = A_in.CSR_I_row_indices;
	CSR_J_col_indices = A_in.CSR_J_col_indices;
	CSR_V_values	  = A_in.CSR_V_values;

	// Dense data 
	dense_values	  = A_in.dense_values; 
	dense_values_fl   = A_in.dense_values_fl;

	// GPU
	d_dense_values    = A_in.d_dense_values; 
	d_x_in			  = A_in.d_x_in;
	d_y_out			  = A_in.d_y_out;

	d_dense_values_fl = A_in.d_dense_values_fl; 
	d_x_in_fl		  = A_in.d_x_in_fl;
	d_y_out_fl		  = A_in.d_y_out_fl;

#ifdef CUDA
	handle			  = A_in.handle;
	stream			  = A_in.stream;
#endif	

}

SparseMatrix& SparseMatrix::operator= ( const SparseCSRMatrix<eslocal> &A_in ) {

	rows = A_in.rows();
	cols = A_in.columns();
	nnz  = A_in.rowPtrs()[rows];
	type = 'G';

	int offset = (A_in.rowPtrs()[0]) ? 0 : 1;
	nnz -= A_in.rowPtrs()[0];

	CSR_I_row_indices.resize(rows+1);
	CSR_J_col_indices.resize(nnz);
	CSR_V_values	 .resize(nnz);

	// Sparse CSR data
	//copy(rows, rows + n_cols + 1, K.CSR_I_row_indices.begin());
	for (int i = 0; i < CSR_I_row_indices.size(); i++)
		CSR_I_row_indices[i] = A_in.rowPtrs()[i] + offset;

	//copy(cols, cols + nnz, K.CSR_J_col_indices.begin());
	for (int i = 0; i < CSR_J_col_indices.size(); i++)
		CSR_J_col_indices[i] = A_in.columnIndices()[i] + offset;

	copy(A_in.values(), A_in.values() + nnz, CSR_V_values.begin());

//	// Sparse COO data
//	I_row_indices = NULL;
//	J_col_indices = NULL;
//	V_values	  = NULL;
//
//	// Dense data
//	dense_values	  = NULL;
//	dense_values_fl   = NULL;

	// GPU
	d_dense_values = NULL;
	d_x_in		   = NULL;
	d_y_out		   = NULL;

	d_dense_values_fl = NULL;
	d_x_in_fl		  = NULL;
	d_y_out_fl		  = NULL;

#ifdef CUDA
	handle		    = NULL;
	stream          = NULL;
#endif

	// return the existing object
	return *this;
}

SparseMatrix::SparseMatrix( const SparseCSRMatrix<eslocal> &A_in, char type_in ) {

	rows = A_in.rows();
	cols = A_in.columns();
	nnz  = A_in.rowPtrs()[rows];
	type = type_in;

	int offset = (A_in.rowPtrs()[0]) ? 0 : 1;
	nnz -= A_in.rowPtrs()[0];

	CSR_I_row_indices.resize(rows+1);
	CSR_J_col_indices.resize(nnz);
	CSR_V_values	 .resize(nnz);

	// Sparse CSR data
	//copy(rows, rows + n_cols + 1, K.CSR_I_row_indices.begin());
	for (int i = 0; i < CSR_I_row_indices.size(); i++)
		CSR_I_row_indices[i] = A_in.rowPtrs()[i] + offset;

	//copy(cols, cols + nnz, K.CSR_J_col_indices.begin());
	for (int i = 0; i < CSR_J_col_indices.size(); i++)
		CSR_J_col_indices[i] = A_in.columnIndices()[i] + offset;

	copy(A_in.values(), A_in.values() + nnz, CSR_V_values.begin());

//	// Sparse COO data
//	I_row_indices = NULL;
//	J_col_indices = NULL;
//	V_values	  = NULL;
//
//	// Dense data
//	dense_values	  = NULL;
//	dense_values_fl   = NULL;

	// GPU
	d_dense_values = NULL;
	d_x_in		   = NULL;
	d_y_out		   = NULL;

	d_dense_values_fl = NULL;
	d_x_in_fl		  = NULL;
	d_y_out_fl		  = NULL;

#ifdef CUDA
	handle		    = NULL;
	stream          = NULL;
#endif

}

SparseMatrix& SparseMatrix::operator= ( const SparseIJVMatrix<eslocal> &A_in ) {

	rows = A_in.rows();
	cols = A_in.columns();
	nnz  = A_in.nonZeroValues();
	type = 'G';

	int offset = A_in.indexing() ? 0 : 1;

	I_row_indices.resize(nnz);
	J_col_indices.resize(nnz);
	V_values	 .resize(nnz);

	// Sparse CSR data
	//copy(rows, rows + n_cols + 1, K.CSR_I_row_indices.begin());
	for (int i = 0; i < I_row_indices.size(); i++)
		I_row_indices[i] = A_in.rowIndices()[i] + offset;

	//copy(cols, cols + nnz, K.CSR_J_col_indices.begin());
	for (int i = 0; i < J_col_indices.size(); i++)
		J_col_indices[i] = A_in.columnIndices()[i] + offset;

	copy(A_in.values(), A_in.values() + nnz, V_values.begin());

	// Sparse CSR data
//	CSR_I_row_indices = NULL;
//	CSR_J_col_indices = NULL;
//	CSR_V_values	  = NULL;

	// Dense data
	//dense_values.clear();
	//dense_values_fl.clear();

	// GPU
	d_dense_values = NULL;
	d_x_in		   = NULL;
	d_y_out		   = NULL;

	d_dense_values_fl = NULL;
	d_x_in_fl		  = NULL;
	d_y_out_fl		  = NULL;

#ifdef CUDA
	handle		    = NULL;
	stream          = NULL;
#endif

	// return the existing object
	return *this;

}

SparseMatrix::SparseMatrix( const SparseIJVMatrix<eslocal> &A_in, char type_in ) {

	rows = A_in.rows();
	cols = A_in.columns();
	nnz  = A_in.nonZeroValues();
	type = type_in;

	int offset = A_in.indexing() ? 0 : 1;

	I_row_indices.resize(nnz);
	J_col_indices.resize(nnz);
	V_values	 .resize(nnz);

	// Sparse CSR data
	//copy(rows, rows + n_cols + 1, K.CSR_I_row_indices.begin());
	for (int i = 0; i < I_row_indices.size(); i++)
		I_row_indices[i] = A_in.rowIndices()[i] + offset;

	//copy(cols, cols + nnz, K.CSR_J_col_indices.begin());
	for (int i = 0; i < J_col_indices.size(); i++)
		J_col_indices[i] = A_in.columnIndices()[i] + offset;

	copy(A_in.values(), A_in.values() + nnz, V_values.begin());

//	// Sparse CSR data
//	CSR_I_row_indices = NULL;
//	CSR_J_col_indices = NULL;
//	CSR_V_values	  = NULL;
//
//	// Dense data
//	dense_values	  = NULL;
//	dense_values_fl   = NULL;

	// GPU
	d_dense_values = NULL;
	d_x_in		   = NULL;
	d_y_out		   = NULL;

	d_dense_values_fl = NULL;
	d_x_in_fl		  = NULL;
	d_y_out_fl		  = NULL;

#ifdef CUDA
	handle		    = NULL;
	stream          = NULL;
#endif

}


SparseMatrix& SparseMatrix::operator= (const SparseMatrix &A_in) {

	// check for self-assignment by comparing the address of the
	// implicit object and the parameter
	if (this == &A_in)
		return *this;
	
	// do the copy
	
		rows = A_in.rows; 
		cols = A_in.cols;
		nnz  = A_in.nnz;
		type = A_in.type;

		// Sparse COO data 
		I_row_indices = A_in.I_row_indices;
		J_col_indices = A_in.J_col_indices;
		V_values	  = A_in.V_values;

		// Sparse CSR data 
		CSR_I_row_indices = A_in.CSR_I_row_indices;
		CSR_J_col_indices = A_in.CSR_J_col_indices;
		CSR_V_values	  = A_in.CSR_V_values;

		// Dense data 
		dense_values	  = A_in.dense_values; 
		dense_values_fl   = A_in.dense_values_fl;

		// GPU
		d_dense_values    = A_in.d_dense_values; 
		d_x_in			  = A_in.d_x_in;
		d_y_out			  = A_in.d_y_out;

		d_dense_values_fl = A_in.d_dense_values_fl; 
		d_x_in_fl		  = A_in.d_x_in_fl;
		d_y_out_fl		  = A_in.d_y_out_fl;

#ifdef CUDA
		handle			  = A_in.handle;
		stream			  = A_in.stream;
#endif	


	// return the existing object
	return *this;
}

SparseMatrix::~SparseMatrix() {
	Clear();
}

void SparseMatrix::Clear() {

	rows = 0; 
	cols = 0; 
	nnz  = 0;  
	type = 0; 

	I_row_indices.clear();		
	J_col_indices.clear();		
	V_values.clear();			

	CSR_I_row_indices.clear();	
	CSR_J_col_indices.clear();	
	CSR_V_values.clear();		

	dense_values.clear();		
	
	SEQ_VECTOR<int>().swap( I_row_indices ); 
	SEQ_VECTOR<int>().swap( J_col_indices ); 
	SEQ_VECTOR<double>().swap( V_values ); 

	SEQ_VECTOR<int>().swap( CSR_I_row_indices ); 
	SEQ_VECTOR<int>().swap( CSR_J_col_indices ); 
	SEQ_VECTOR<double>().swap( CSR_V_values );

	SEQ_VECTOR<double>().swap( dense_values );  
	SEQ_VECTOR<float>().swap( dense_values_fl );  
	
	// GPU
	//d_dense_values = NULL; 
	//d_y_out        = NULL; 
	//d_x_in		   = NULL: 
}


int SparseMatrix::LoadMatrixBinInCOO(string filename, char matrix_type_G_for_general_S_for_symmetric) {
	
	type = matrix_type_G_for_general_S_for_symmetric;

	ifstream in (filename.c_str(), std::ios::binary);

	if ( in.is_open() ) {

		char delim = ';'; 
		string line, field;

		// Get parameters 
		getline(in,line);
		stringstream paramss(line);

		getline(paramss,field,delim); 
		rows = atoi(field.c_str());		// get num of rows 

		getline(paramss,field,delim); 
		cols = atoi(field.c_str());		// get num of columns 

		getline(paramss,field,delim); 
		nnz  = atoi(field.c_str());		// get num of non zero elements 

		// Get data 
		I_row_indices.resize(nnz);
		in.read((char*) &I_row_indices[0], nnz*sizeof(int));

		J_col_indices.resize(nnz);
		in.read((char*) &J_col_indices[0], nnz*sizeof(int));

		V_values.resize(nnz);
		in.read((char*) &V_values[0], nnz*sizeof(double));

		in.close();
	
		return 0;

	} else {
	
		cout << "Matrix file " << filename << " not found ! " << endl; 
		return -1; 

	}
}

int SparseMatrix::LoadMatrixBin(string filename, char matrix_type_G_for_general_S_for_symmetric) {
	int tmp = LoadMatrixBinInCOO(filename, matrix_type_G_for_general_S_for_symmetric);
	
	if (tmp == 0)	
		ConvertToCSR( 1 ); 
	
	return tmp;
}

int SparseMatrix::LoadMatrixBin(string filename, char matrix_type_G_for_general_S_for_symmetric, int clearCOO_1_keep_COO_0 ) {
	int tmp = LoadMatrixBinInCOO(filename, matrix_type_G_for_general_S_for_symmetric);
	
	if (tmp == 0)
		ConvertToCSR( clearCOO_1_keep_COO_0 ); 
	
	return tmp; 
}

//void SparseMatrix::LoadMatrixInCOO(string filename, char matrix_type_G_for_general_S_for_symmetric) {
//
//	type = matrix_type_G_for_general_S_for_symmetric;
//
//	ifstream in (filename.c_str());
//	char delim = ';'; 
//	string line, field;
//
//	// Get line with matrix parameters 
//	getline(in,line);
//	stringstream paramss(line);
//
//	getline(paramss,field,delim); 
//	rows = atoi(field.c_str());		// get num of rows 
//
//	getline(paramss,field,delim); 
//	cols = atoi(field.c_str());		// get num of columns 
//
//	getline(paramss,field,delim); 
//	nnz  = atoi(field.c_str());		// get num of non zero elements 
//	
//	// Get line with I values - INT 
//	getline(in,line); 
//	stringstream ssI(line);
//	while (getline(ssI,field,delim))  // break line into comma delimitted fields
//	{
//		I_row_indices.push_back(atoi(field.c_str()));  // add each field to the 1D array
//	}
//
//
//	// Get line with J values - INT 
//	getline(in,line); 
//	stringstream ssJ(line);
//	while (getline(ssJ,field,delim))  // break line into comma delimitted fields
//	{
//		J_col_indices.push_back(atoi(field.c_str()));  // add each field to the 1D array
//	}
//
//
//	// Get line with V values - DOUBLE 
//	getline(in,line); 
//	stringstream ssV(line);
//	while (getline(ssV,field,delim))  // break line into comma delimitted fields
//	{
//		V_values.push_back(atof(field.c_str()));  // add each field to the 1D array
//	}
//
//	in.close();
//}
//
//void SparseMatrix::LoadMatrix(string filename, char matrix_type_G_for_general_S_for_symmetric) {
//	LoadMatrixInCOO(filename, matrix_type_G_for_general_S_for_symmetric);
//	ConvertToCSR( 1 ); 
//}





void SparseMatrix::ConvertToCSR( ) {
	ConvertToCSR( 1 ); 
}

void SparseMatrix::ConvertToCSRwithSort( int clearCOO_1_keep_COO_0 ) {

	int job[8];//  = {0,0,0,0, 0,0,0,0};
	job[0] = 2; // the matrix in the coordinate format is converted to the CSR format
	job[1] = 1; // one-based indexing for the matrix in CSR format is used.
	job[2] = 1; // one-based indexing for the matrix in coordinate format is used.
	job[3] = 0;

	job[4] = nnz; //nnz - sets number of the non-zero elements of the matrix A if job(1)=1.
	job[5] = 0;   // For conversion to the CSR format: If job(6) = 0, all arrays acsr, ja, ia are filled in for the output storage.
	job[6] = 0;
	job[7] = 0;

	CSR_V_values.resize(nnz);			// Array containing non-zero elements of the matrix A. Its length is equal to the number of non-zero elements in the matrix A.
	CSR_J_col_indices.resize(nnz);		// Array containing the column indices for each non-zero element of the matrix A. Its length is equal to the length of the array acsr.
	CSR_I_row_indices.resize(rows+1);	// Array of length n + 1, containing indices of elements in the array acsr, such that ia(I) is the index in the array acsr of the first non-zero element from the row I. The value of the last element ia(n + 1) is equal to the number of non-zeros plus one.

	int info; 

	//void mkl_dcsrcoo ( MKL_INT * job, MKL_INT * n, double *Acsr,      MKL_INT * AJR,          MKL_INT * AIR,          MKL_INT * nnz,  double *Acoo,  MKL_INT * ir,       MKL_INT * jc,       MKL_INT * info);
	mkl_dcsrcoo		   ( job,           &rows,       &CSR_V_values[0],  &CSR_J_col_indices[0],  &CSR_I_row_indices[0],  &nnz,           &V_values[0],  &I_row_indices[0],  &J_col_indices[0],  &info ); 

	if (clearCOO_1_keep_COO_0 == 1 ) {

		V_values.clear();
		I_row_indices.clear();
		J_col_indices.clear(); 

		SEQ_VECTOR<double>().swap( V_values ); 
		SEQ_VECTOR<int>().swap( I_row_indices ); 
		SEQ_VECTOR<int>().swap( J_col_indices ); 

	}

	// m	INTEGER. Number of rows of the matrix A.
	// n	INTEGER. Number of columns of the matrix A.

}


void SparseMatrix::ConvertToCSR( int clearCOO_1_keep_COO_0 ) {

	int job[8];//  = {0,0,0,0, 0,0,0,0};
	job[0] = 2; // the matrix in the coordinate format is converted to the CSR format, and the column indices in CSR representation are sorted in the increasing order within each row.
	job[1] = 1; // one-based indexing for the matrix in CSR format is used.
	job[2] = 1; // one-based indexing for the matrix in coordinate format is used.
	job[3] = 0;

	job[4] = nnz; //nnz - sets number of the non-zero elements of the matrix A if job(1)=1.
	job[5] = 0;   // For conversion to the CSR format: If job(6) = 0, all arrays acsr, ja, ia are filled in for the output storage.
	job[6] = 0;
	job[7] = 0;

	CSR_V_values.resize(nnz);			// Array containing non-zero elements of the matrix A. Its length is equal to the number of non-zero elements in the matrix A.
	CSR_J_col_indices.resize(nnz);		// Array containing the column indices for each non-zero element of the matrix A. Its length is equal to the length of the array acsr.
	CSR_I_row_indices.resize(rows+1);	// Array of length n + 1, containing indices of elements in the array acsr, such that ia(I) is the index in the array acsr of the first non-zero element from the row I. The value of the last element ia(n + 1) is equal to the number of non-zeros plus one.

	int info; 

	//void mkl_dcsrcoo ( MKL_INT * job, MKL_INT * n, double *Acsr,      MKL_INT * AJR,          MKL_INT * AIR,          MKL_INT * nnz,  double *Acoo,  MKL_INT * ir,       MKL_INT * jc,       MKL_INT * info);
	mkl_dcsrcoo		   ( job,           &rows,       &CSR_V_values[0],  &CSR_J_col_indices[0],  &CSR_I_row_indices[0],  &nnz,           &V_values[0],  &I_row_indices[0],  &J_col_indices[0],  &info ); 

	if (clearCOO_1_keep_COO_0 == 1 ) {
		
		V_values.clear();
		I_row_indices.clear();
		J_col_indices.clear(); 

		SEQ_VECTOR<double>().swap( V_values ); 
		SEQ_VECTOR<int>().swap( I_row_indices ); 
		SEQ_VECTOR<int>().swap( J_col_indices ); 

	}

	// m	INTEGER. Number of rows of the matrix A.
	// n	INTEGER. Number of columns of the matrix A.

}

void SparseMatrix::ConvertToCOO( int clearCSR_1_keep_CSR_0 ) {

	int job[8];//  = {0,0,0,0, 0,0,0,0};
	job[0] = 0; // job(1)=0, the matrix in the CSR format is converted to the coordinate format;	
	job[1] = 1; // one-based indexing for the matrix in CSR format is used.
	job[2] = 1; // one-based indexing for the matrix in coordinate format is used.
	job[3] = 0;

	job[4] = nnz; //nnz - sets number of the non-zero elements of the matrix A if job(1)=1.
	job[5] = 3;   // For conversion to the coordinate format: If job(6)=3, all arrays rowind, colind, acoo are filled in for the output storage.
	job[6] = 0;
	job[7] = 0;

	V_values.resize(nnz);			
	J_col_indices.resize(nnz);		
	I_row_indices.resize(nnz);	    

	int info; 

	//void mkl_dcsrcoo ( MKL_INT * job, MKL_INT * n, double *Acsr,      MKL_INT * AJR,          MKL_INT * AIR,          MKL_INT * nnz,  double *Acoo,  MKL_INT * ir,       MKL_INT * jc,       MKL_INT * info);
	mkl_dcsrcoo		   ( job,           &rows,       &CSR_V_values[0],  &CSR_J_col_indices[0],  &CSR_I_row_indices[0],  &nnz,           &V_values[0],  &I_row_indices[0],  &J_col_indices[0],  &info ); 

	if (clearCSR_1_keep_CSR_0 == 1 ) {

		CSR_V_values.clear();
		CSR_I_row_indices.clear();
		CSR_J_col_indices.clear(); 

		SEQ_VECTOR<int>().swap( CSR_I_row_indices ); 
		SEQ_VECTOR<int>().swap( CSR_J_col_indices ); 
		SEQ_VECTOR<double>().swap( CSR_V_values );
		
	}

	// m	INTEGER. Number of rows of the matrix A.
	// n	INTEGER. Number of columns of the matrix A.

}


void SparseMatrix::ConvertCSRToDense( int clearCSR_1_keep_CSR_0 ) {

	int job[8];
	job[0] = 1; // if job(1)=1, the rectangular matrix A is restored from the CSR format.
	job[1] = 1; // if job(2)=1, one-based indexing for the rectangular matrix A is used.
	job[2] = 1; // if job(3)=1, one-based indexing for the matrix in CSR format is used.
	job[3] = 2; // If job(4)=2, adns is a whole matrix A.
	job[4] = 0; // job(5)=nzmax - maximum number of the non-zero elements allowed if job(1)=0.
	job[5] = 0; // job(6) - job indicator for conversion to CSR format.
				// If job(6)=0, only array ia is generated for the output storage.
				// If job(6)>0, arrays acsr, ia, ja are generated for the output storage.
	job[6] = 0; //
	job[7] = 0; //

	int m		= rows;
	int n		= cols; 
	int lda     = m;
	int info    = 0; 

	dense_values.resize(m * n, 0); 

	// Convert matrix to dense format 

	//void mkl_ddnscsr (
	//	MKL_INT *job, 
	//	MKL_INT *m, MKL_INT *n, 
	//	double *Adns, MKL_INT *lda, 
	//	double *Acsr, MKL_INT *AJ, MKL_INT *AI, 
	//	MKL_INT *info);

	mkl_ddnscsr (
		job, 
		&m, &n, 
		&dense_values[0], &lda, 
		&CSR_V_values[0], &CSR_J_col_indices[0], &CSR_I_row_indices[0], 
		&info);

	if (clearCSR_1_keep_CSR_0 == 1) {
		CSR_I_row_indices.clear();
		CSR_J_col_indices.clear();
		CSR_V_values.clear(); 

		SEQ_VECTOR<int>().swap( CSR_I_row_indices ); 
		SEQ_VECTOR<int>().swap( CSR_J_col_indices ); 
		SEQ_VECTOR<double>().swap( CSR_V_values );
	}

	if (type == 'S')
		this->RemoveLowerDense();

}

void SparseMatrix::ConvertDenseToCSR( int clearDense_1_keep_Dense_0 ){

	int m		= rows;
	int n		= cols; 
	int lda     = m;
	int info    = 0; 

	int job[8];
	
	// Convert to sparse format - find nnz step
	job[0] = 0; // If job(1)=0, the rectangular matrix A is converted to the CSR format;
	job[1] = 1; // if job(2)=1, one-based indexing for the rectangular matrix A is used.
	job[2] = 1; // if job(3)=1, one-based indexing for the matrix in CSR format is used.
	job[3] = 2; // If job(4)=2, adns is a whole matrix A.

	job[4] = 1; // job(5)=nzmax - maximum number of the non-zero elements allowed if job(1)=0.
	job[5] = 0; // job(6) - job indicator for conversion to CSR format.
				// If job(6)=0, only array ia is generated for the output storage.
				// If job(6)>0, arrays acsr, ia, ja are generated for the output storage.
	job[6] = 0; //
	job[7] = 0; //

	CSR_I_row_indices.resize(m + 1); 
	CSR_J_col_indices.resize(1);
	CSR_V_values.resize(1);

	mkl_ddnscsr (
		job, 
		&m, &n, 
		&dense_values[0], &lda, 
		&CSR_V_values[0], &CSR_J_col_indices[0], &CSR_I_row_indices[0], 
		&info);

	// Convert to sparse format - convert step 
	int nnzmax = CSR_I_row_indices[m] - 1; //- 1;  POZOR bez -1 polud se to pouzije ve funkci SolveMatF
	CSR_J_col_indices.resize(nnzmax);  // 
	CSR_V_values.resize(nnzmax);       // 

	job[4] = nnzmax; // job(5) = nzmax - maximum number of the non-zero elements allowed if job(1)=0.
	job[5] = 1; // job(6) - job indicator for conversion to CSR format.
				// If job(6)=0, only array ia is generated for the output storage.
				// If job(6)>0, arrays acsr, ia, ja are generated for the output storage.

	mkl_ddnscsr (
		job, 
		&m, &n, 
		&dense_values[0], &lda,
		&CSR_V_values[0], &CSR_J_col_indices[0], &CSR_I_row_indices[0], 
		&info);
	
	// Setup parameters for output matrix 
	nnz	= nnzmax; //POZOR  CSR_V_values.size();

	if (clearDense_1_keep_Dense_0 == 1) {
		dense_values.clear();  
		SEQ_VECTOR<double>().swap( dense_values );
	}

}

void SparseMatrix::ConvertDenseToDenseFloat( int clear_DoubleDense_1_keep_DoubleDense_0 ) {

	dense_values_fl.resize( dense_values.size() );
	
	for (int i = 0; i < dense_values.size(); i++)
		dense_values_fl[i] = (float)dense_values[i];

	if ( clear_DoubleDense_1_keep_DoubleDense_0 == 1)
		SEQ_VECTOR<double>().swap( dense_values );  

}

//void SparseMatrix::DenseMatVec(vector <double> & x_in, vector <double> & y_out) {
//
//	// void cblas_dgemv 
//	//  (const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA, 
//	//  const MKL_INT M, const MKL_INT N, 
//	//  const double alpha, const double *A, const MKL_INT lda, 
//	//  const double *X, const MKL_INT incX, 
//	//  const double beta, double *Y, const MKL_INT incY);
//
//	// y := alpha*A*x + beta*y,
//
//	double alpha = 1.0; 
//	double beta  = 0.0; 
//	int lda = rows; 
//
//	//char trans = 'N'; // CblasNoTrans=111,     /* trans='N' */
//						// CblasTrans=112,       /* trans='T' */
//						// CblasConjTrans=113};  /* trans='C' */
//
//	cblas_dgemv 
//	(CblasColMajor, CblasNoTrans,  
//	rows, cols, 
//	alpha, &dense_values[0], lda, 
//	&x_in[0], 1, 
//	beta, &y_out[0], 1);
//
//}

void SparseMatrix::DenseMatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out) {
	DenseMatVec(x_in, y_out, 'N', 0); 
}


void SparseMatrix::DenseMatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose ) {
	DenseMatVec(x_in, y_out, T_for_transpose_N_for_not_transpose, 0); 
}

void SparseMatrix::DenseMatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose, int x_in_vector_start_index) {

	// void cblas_dgemv 
	//  (const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA, 
	//  const MKL_INT M, const MKL_INT N, 
	//  const double alpha, const double *A, const MKL_INT lda, 
	//  const double *X, const MKL_INT incX, 
	//  const double beta, double *Y, const MKL_INT incY);

	// y := alpha*A*x + beta*y,

	double alpha = 1.0; 
	double beta  = 0.0; 
	int lda = rows; 

	//char trans = 'N'; // CblasNoTrans=111,     /* trans='N' */
	// CblasTrans=112,       /* trans='T' */
	// CblasConjTrans=113};  /* trans='C' */

	if ( T_for_transpose_N_for_not_transpose == 'T' ) 
		cblas_dgemv 
			(CblasColMajor, CblasTrans,  
			rows, cols, 
			alpha, &dense_values[0], lda, 
			&x_in[x_in_vector_start_index], 1, 
			beta, &y_out[0], 1);
	else 
		cblas_dgemv 
			(CblasColMajor, CblasNoTrans,  
			rows, cols, 
			alpha, &dense_values[0], lda, 
			&x_in[x_in_vector_start_index], 1, 
			beta, &y_out[0], 1);

}


void SparseMatrix::DenseMatVecCUDA_w_Copy(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose, int x_in_vector_start_index) {
#ifdef CUDA 

	double *d_x_in_t, *d_y_out_t, *d_Mat;
	int mat_size = rows * cols; 
	int lda = rows;

	cudaMalloc((void**)&d_x_in_t,  x_in.size()  * sizeof(double));
	cudaMalloc((void**)&d_y_out_t, y_out.size() * sizeof(double));
	cudaMalloc((void**)&d_Mat,   mat_size     * sizeof(double));

	// Create cublas instance
	cublasHandle_t handle;
	cublasCreate(&handle);

	// Set input matrices on device
	// cublasSetVector(int n, int elemSize, const void *x, int incx, void *y, int incy);
	cublasSetVector(x_in.size() , sizeof(double), &x_in[0] , 1, d_x_in_t , 1);
	cublasSetVector(y_out.size(), sizeof(double), &y_out[0], 1, d_y_out_t, 1);
	cublasSetMatrix(rows, cols  , sizeof(double), &dense_values[0], lda, d_Mat, lda);

	// DGEMM: C = alpha*A*B + beta*C
	double alpha = 1.0;
	double beta  = 0.0;

	if ( T_for_transpose_N_for_not_transpose == 'T' ) {
		cublasDgemv(handle, 
			CUBLAS_OP_T, 
			rows, cols, 
			&alpha, d_Mat, lda, 
			d_x_in_t, 1, 
			&beta, d_y_out_t, 1);
	} else {
		cublasDgemv(handle, 
			CUBLAS_OP_N,  
			rows, cols, 
			&alpha, d_Mat, lda, 
			d_x_in_t, 1, 
			&beta, d_y_out_t, 1);
	}
	
	// Retrieve result vector from device
	//cublasGetVector(int n, int elemSize, const void *x, int incx, void *y, int incy)
	cublasGetVector(y_out.size(), sizeof(double), d_y_out_t, 1, &y_out[0], 1);

	cudaFree(d_x_in_t);
	cudaFree(d_y_out_t);
	cudaFree(d_Mat);
	cublasDestroy(handle);


#endif
}


void SparseMatrix::DenseMatVecCUDA_wo_Copy(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose, int x_in_vector_start_index) {
#ifdef CUDA 

	int mat_size = rows * cols; 
	int lda = rows;

	if ( d_dense_values == NULL ) {
		cudaMalloc((void**)&d_dense_values,   mat_size * sizeof(double));
		cudaMalloc((void**)&d_x_in,  rows * sizeof(double));
		cudaMalloc((void**)&d_y_out, rows * sizeof(double));
		
		cublasSetMatrix(rows, cols  , sizeof(double), &dense_values[0], lda, d_dense_values, lda);
		cublasSetVector(rows,         sizeof(double), &dense_values[0], 1,   d_y_out,        1);

		// Create cublas instance - cublasHandle_t handle;
		
		cublasCreate(&handle);

		cudaStreamCreate(&stream);
		cublasSetStream(handle, stream);

		cout << "Y"; 
	}
	 
	// Set input matrices on device
	// cublasSetVector(int n, int elemSize, const void *x, int incx, void *y, int incy);
	cublasSetVector(x_in.size() , sizeof(double), &x_in[0] , 1, d_x_in , 1);

	// DGEMM: C = alpha*A*B + beta*C
	double alpha = 1.0;
	double beta  = 0.0;

	if ( T_for_transpose_N_for_not_transpose == 'T' ) {
		cublasDgemv(handle, 
			CUBLAS_OP_T, 
			rows, cols, 
			&alpha, d_dense_values, lda, 
			d_x_in, 1, 
			&beta, d_y_out, 1);
	} else {
		cublasDgemv(handle, 
			CUBLAS_OP_N,  
			rows, cols, 
			&alpha, d_dense_values, lda, 
			d_x_in, 1, 
			&beta, d_y_out, 1);
	}

	// Retrieve result vector from device
	//cublasGetVector(int n, int elemSize, const void *x, int incx, void *y, int incy)
	cublasGetVector(rows , sizeof(double), d_y_out, 1, &y_out[0], 1);

	//cudaFree(d_x_in_t);
	//cudaFree(d_y_out_t);
	//cudaFree(d_dense_values);
	//cublasDestroy(handle);
	
	cudaStreamSynchronize(stream);

#endif
}




void SparseMatrix::DenseMatVecCUDA_wo_Copy_start( double * x_in, double * y_out, char T_for_transpose_N_for_not_transpose, int x_in_vector_start_index) {
#ifdef CUDA 

	int mat_size = dense_values.size(); //rows * cols; 
	int lda = rows;

	if ( d_dense_values == NULL ) {
		cudaMalloc((void**)&d_dense_values,   mat_size * sizeof(double));
		cudaMalloc((void**)&d_x_in,  rows * sizeof(double));
		cudaMalloc((void**)&d_y_out, rows * sizeof(double));

		//cublasSetMatrix(rows, cols  , sizeof(double), &dense_values[0], lda, d_dense_values, lda);
		cudaMemcpy(d_dense_values, &dense_values[0], dense_values.size() * sizeof(double), cudaMemcpyHostToDevice);

		cublasSetVector(rows,         sizeof(double), &dense_values[0], 1,   d_y_out,        1);

		// Create cublas instance - cublasHandle_t handle;

		cublasCreate(&handle);

		cudaStreamCreate(       &stream);
		cublasSetStream (handle, stream);

		cout << "Y"; 
	}

	// Set input matrices on device
	// cublasSetVector(int n, int elemSize, const void *x, int incx, void *y, int incy);
	//cublasSetVector(rows, sizeof(double), x_in , 1, d_x_in , 1);
	cudaMemcpyAsync(d_x_in, x_in, rows * sizeof(double), cudaMemcpyHostToDevice, stream);

	// DGEMM: C = alpha*A*B + beta*C
	double alpha = 1.0;
	double beta  = 0.0;

	if ( T_for_transpose_N_for_not_transpose == 'T' ) {
		cublasDgemv(handle, 
			CUBLAS_OP_T, 
			rows, cols, 
			&alpha, d_dense_values, lda, 
			d_x_in, 1, 
			&beta, d_y_out, 1);
	} else {
		cublasDgemv(handle, 
			CUBLAS_OP_N,  
			rows, cols, 
			&alpha, d_dense_values, lda, 
			d_x_in, 1, 
			&beta, d_y_out, 1);
	}

		// jen o neco malo pomalejsi - cca 5% porad dobre, ale usetri se 50% pameti 
		//cublasDsymv(handle, 
		//	CUBLAS_FILL_MODE_UPPER, //CUBLAS_FILL_MODE_LOWER
		//	rows, 
		//	&alpha, d_dense_values, lda, 
		//	d_x_in, 1, 
		//	&beta, d_y_out, 1);

		// POMALE 
		//cublasDspmv(handle, 
		//	CUBLAS_FILL_MODE_UPPER, 
		//	rows, 
		//	&alpha, d_dense_values, 
		//	d_x_in, 1, 
		//	&beta, d_y_out, 1);


	// Retrieve result vector from device
	//cublasGetVector(int n, int elemSize, const void *x, int incx, void *y, int incy)
	//cublasGetVector(rows , sizeof(double), d_y_out, 1, &y_out[0], 1);
	cudaMemcpyAsync(y_out, d_y_out, rows * sizeof(double), cudaMemcpyDeviceToHost, stream);

	//cudaStreamSynchronize(stream);

#endif
}

void SparseMatrix::DenseMatVecCUDA_wo_Copy_sync ( ) {
#ifdef CUDA 

	// Retrieve result vector from device
	//cublasGetVector(int n, int elemSize, const void *x, int incx, void *y, int incy)
	//cublasGetVector(rows , sizeof(double), d_y_out, 1, &y_out[0], 1);
	cudaStreamSynchronize(stream);

#endif
}

void SparseMatrix::CopyToCUDA_Dev( ) {
#ifdef CUDA	 

	int mat_size = dense_values.size();// rows * cols; 
	int lda = rows;

	if ( d_dense_values == NULL ) {

		cudaError_t status = cudaMalloc((void**)&d_dense_values,   mat_size * sizeof(double));
		if (status != cudaSuccess)   {
			printf("Error allocating GPU memory \n");
			MPI_Finalize();
			exit(0);
		}
		

		status = cudaMalloc((void**)&d_x_in,  rows * sizeof(double));
		if (status != cudaSuccess) {
			printf("Error allocating GPU memory  \n");
			MPI_Finalize();
			exit(0);
		}
		

		status = cudaMalloc((void**)&d_y_out, rows * sizeof(double));
		if (status != cudaSuccess) {
			printf("Error allocating GPU memory \n");
			MPI_Finalize();
			exit(0);
		}
		

		//cublasSetMatrix(rows, cols  , sizeof(double), &dense_values[0], lda, d_dense_values, lda);
		cudaMemcpy(d_dense_values, &dense_values[0], dense_values.size() * sizeof(double), cudaMemcpyHostToDevice);

		cublasSetVector(rows,         sizeof(double), &dense_values[0], 1,   d_y_out,        1);
		// Create cublas instance - cublasHandle_t handle;

		cublasCreate(&handle);

		cudaStreamCreate(&stream);
		cublasSetStream(handle, stream);
		 
		cout << "X"; 
			
	}

#endif
}



void SparseMatrix::DenseMatVecCUDA_wo_Copy_start_fl( float * x_in, float * y_out, char T_for_transpose_N_for_not_transpose, int x_in_vector_start_index) {
#ifdef CUDA 

	int mat_size = dense_values.size(); //rows * cols; 
	int lda = rows;

	if ( d_dense_values_fl == NULL ) {
		cudaMalloc((void**)&d_dense_values_fl,   mat_size * sizeof(float));
		cudaMalloc((void**)&d_x_in_fl,  rows * sizeof(float));
		cudaMalloc((void**)&d_y_out_fl, rows * sizeof(float));

		//cublasSetMatrix(rows, cols  , sizeof(double), &dense_values[0], lda, d_dense_values, lda);
		cudaMemcpy(d_dense_values_fl, &dense_values_fl[0], dense_values_fl.size() * sizeof(float), cudaMemcpyHostToDevice);

		cublasSetVector(rows,         sizeof(float), &dense_values_fl[0], 1,   d_y_out_fl,        1);

		// Create cublas instance - cublasHandle_t handle;

		cublasCreate(&handle);

		cudaStreamCreate(       &stream);
		cublasSetStream (handle, stream);

		cout << "Y"; 
	}

	// Set input matrices on device
	// cublasSetVector(int n, int elemSize, const void *x, int incx, void *y, int incy);
	//cublasSetVector(rows, sizeof(double), x_in , 1, d_x_in , 1);
	cudaMemcpyAsync(d_x_in_fl, x_in, rows * sizeof(float), cudaMemcpyHostToDevice, stream);

	// DGEMM: C = alpha*A*B + beta*C
	float alpha = 1.0;
	float beta  = 0.0;

	if ( T_for_transpose_N_for_not_transpose == 'T' ) {
		cublasSgemv(handle, 
			CUBLAS_OP_T, 
			rows, cols, 
			&alpha, d_dense_values_fl, lda, 
			d_x_in_fl, 1, 
			&beta, d_y_out_fl, 1);
	} else {
		cublasSgemv(handle, 
			CUBLAS_OP_N,  
			rows, cols, 
			&alpha, d_dense_values_fl, lda, 
			d_x_in_fl, 1, 
			&beta, d_y_out_fl, 1);
	}

	// jen o neco malo pomalejsi - cca 5% porad dobre, ale usetri se 50% pameti 
	//cublasDsymv(handle, 
	//	CUBLAS_FILL_MODE_UPPER, //CUBLAS_FILL_MODE_LOWER
	//	rows, 
	//	&alpha, d_dense_values, lda, 
	//	d_x_in, 1, 
	//	&beta, d_y_out, 1);

	// POMALE 
	//cublasDspmv(handle, 
	//	CUBLAS_FILL_MODE_UPPER, 
	//	rows, 
	//	&alpha, d_dense_values, 
	//	d_x_in, 1, 
	//	&beta, d_y_out, 1);


	// Retrieve result vector from device
	//cublasGetVector(int n, int elemSize, const void *x, int incx, void *y, int incy)
	//cublasGetVector(rows , sizeof(double), d_y_out, 1, &y_out[0], 1);
	cudaMemcpyAsync(y_out, d_y_out_fl, rows * sizeof(float), cudaMemcpyDeviceToHost, stream);

	//cudaStreamSynchronize(stream);

#endif
}

void SparseMatrix::CopyToCUDA_Dev_fl ( ) {
#ifdef CUDA	 

	int mat_size = dense_values.size();// rows * cols; 
	int lda = rows;

	if ( d_dense_values_fl == NULL ) {

		cudaError_t status = cudaMalloc((void**)&d_dense_values_fl,   mat_size * sizeof(float));
		if (status != cudaSuccess)   {
			printf("Error allocating GPU memory \n");
			MPI_Finalize();
			exit(0);
		}


		status = cudaMalloc((void**)&d_x_in_fl,  rows * sizeof(float));
		if (status != cudaSuccess) {
			printf("Error allocating GPU memory  \n");
			MPI_Finalize();
			exit(0);
		}


		status = cudaMalloc((void**)&d_y_out_fl, rows * sizeof(float));
		if (status != cudaSuccess) {
			printf("Error allocating GPU memory \n");
			MPI_Finalize();
			exit(0);
		}


		//cublasSetMatrix(rows, cols  , sizeof(double), &dense_values[0], lda, d_dense_values, lda);
		cudaMemcpy(d_dense_values_fl, &dense_values_fl[0], dense_values_fl.size() * sizeof(float), cudaMemcpyHostToDevice);

		cublasSetVector(rows,         sizeof(float), &dense_values_fl[0], 1,   d_y_out_fl,        1);
		// Create cublas instance - cublasHandle_t handle;

		cublasCreate(&handle);

		cudaStreamCreate(&stream);
		cublasSetStream(handle, stream);

		cout << "X"; 

	}

#endif
}





void SparseMatrix::CopyFromCUDA_Dev() {
#ifdef CCUDA
	cudaFree(d_dense_values);
	d_dense_values = NULL; 
#endif
}

void SparseMatrix::FreeFromCUDA_Dev() {
#ifdef CCUDA
	cudaFree(d_dense_values);
	d_dense_values = NULL; 
#endif
}

void SparseMatrix::RemoveLowerDense( ) {

//	        dtrttp(                         uplo,            n,                   a,            lda,             ap, info )
//	LAPACKE_dtrttp( int matrix_layout, char uplo, lapack_int n, const <datatype>* a, lapack_int lda, <datatype>* ap )
//                      LAPACK_COL_MAJOR
//						LAPACK_ROW_MAJOR
	
	vector <double> tmp_dense ( (rows *(rows +1))/2, 0.0 ) ;

	LAPACKE_dtrttp(  LAPACK_COL_MAJOR,       'U',         rows,    &dense_values[0],           rows,  &tmp_dense[0] );

	dense_values.swap(tmp_dense);  

}


void SparseMatrix::MatVecCOO(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out ) {
	char T_for_transpose_N_for_non_transpose = 'N';
	MatVecCOO(x_in, y_out, T_for_transpose_N_for_non_transpose); 
}

void SparseMatrix::MatVecCOO(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose) {
	double beta = 0.0;
	MatVecCOO(x_in, y_out, T_for_transpose_N_for_non_transpose, beta); 
}

void SparseMatrix::MatVecCOO(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose, double beta) {

	char trans		 = T_for_transpose_N_for_non_transpose; 
	double alpha	 = 1.0; 
	//double beta		 = 0.0;

	char matdescra[] = {0,0,0,0,0,0};

	if (type == 'G') {
		matdescra[0] = 'G'; // General matrix 
		matdescra[1] = '0';  
		matdescra[2] = '0';
		matdescra[3] = 'F'; // one based indexing 
	} 

	if (type == 'S') {
		matdescra[0] = 'T'; // Triangular Matrix  
		matdescra[1] = 'U'; // Triangular indicator: upper 
		matdescra[2] = 'N'; // Main diagonal type: non-unit
		matdescra[3] = 'F'; // One based indexing 
	}
	//y_out.resize(rows);

	// y := alpha*A*x + beta*y
	//void mkl_dcoomv	(char *transa, MKL_INT *m, MKL_INT *k, double *alpha, char *matdescra, double *val,    MKL_INT *rowind,    MKL_INT *colind,    MKL_INT *nnz, double *x,  double *beta, double *y);
	mkl_dcoomv			( &trans,	   &rows,      &cols,      &alpha,        matdescra,      &V_values[0],    &I_row_indices[0],  &J_col_indices[0],  &nnz,         &x_in[0],   &beta,        &y_out[0]);

}



void SparseMatrix::MatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose ) {
	MatVec(x_in, y_out, T_for_transpose_N_for_non_transpose, 0 , 0);
}

void SparseMatrix::MatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose, int x_in_vector_start_index, int y_out_vector_start_index) {
	double beta = 0; 
	MatVec(x_in, y_out, T_for_transpose_N_for_non_transpose, x_in_vector_start_index, y_out_vector_start_index, beta);
}

void SparseMatrix::MatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose, int x_in_vector_start_index, int y_out_vector_start_index, double beta) {
	// y := alpha*A*x + beta*y
	// or 
	// y := alpha*A'*x + beta*y

	char trans;
	trans			 =  T_for_transpose_N_for_non_transpose; 
	double alpha	 =  1; 
	//double beta		 =  0;

	char matdescra[] = {0,0,0,0,0,0};

	if (type == 'G') {
		matdescra[0] = 'G'; // General matrix 
		matdescra[1] = '0'; // Ignored 
		matdescra[2] = '0'; // Ignored
		matdescra[3] = 'F'; // One based indexing 
	}

	if (type == 'S') {
		matdescra[0] = 'S'; // Triangular Matrix  
		matdescra[1] = 'U'; // Triangular indicator: upper 
		matdescra[2] = 'N'; // Main diagonal type: non-unit
		matdescra[3] = 'F'; // One based indexing 
	}

	// y := alpha*A*x + beta*y
	// void mkl_dcsrmv (char *transa, MKL_INT *m, MKL_INT *k, double *alpha, char *matdescra, double *val,       MKL_INT *indx,          MKL_INT *pntrb,     MKL_INT *pntre,     double *x,						   double *beta,  double *y);
	// mkl_ccsrmv      (&transa,      &m,         &m,         &alpha,              matdescra, values,            columns,                rowIndex,           &(rowIndex[1]),     sol_vec,						   &beta,         rhs_vec); 
	mkl_dcsrmv         (&trans,       &rows,      &cols,      &alpha,              matdescra, &CSR_V_values[0],  &CSR_J_col_indices[0],  &CSR_I_row_indices[0],  &CSR_I_row_indices[1],  &x_in[x_in_vector_start_index],   &beta,         &y_out[y_out_vector_start_index]);

}


void SparseMatrix::MatMat(SparseMatrix & A_in, char MatA_T_for_transpose_N_for_non_transpose, SparseMatrix & B_in) {
	// THIS := op(A)*B

	char transa = MatA_T_for_transpose_N_for_non_transpose; 

	int job; 

	int sort = 3;	   // 3	yes	yes	yes

	int m = A_in.rows; // Number of rows of the matrix A.
	int n = A_in.cols; // Number of columns of the matrix A.
	int k = B_in.cols; // Number of columns of the matrix B.
	
	SEQ_VECTOR<int>().swap( CSR_I_row_indices ); 
	SEQ_VECTOR<int>().swap( CSR_J_col_indices ); 
	SEQ_VECTOR<double>().swap( CSR_V_values );

	if (transa == 'T')
		CSR_I_row_indices.resize( n + 1 );
	else 
		CSR_I_row_indices.resize( m + 1 );
	
	CSR_J_col_indices.resize(1);
	CSR_V_values.resize(1); 

	double * a  = &A_in.CSR_V_values[0]; 
	int    * ia = &A_in.CSR_I_row_indices[0]; 
	int    * ja = &A_in.CSR_J_col_indices[0];

	double * b  = &B_in.CSR_V_values[0]; 
	int    * ib = &B_in.CSR_I_row_indices[0]; 
	int    * jb = &B_in.CSR_J_col_indices[0]; 

	int nnzmax = 1;  

	int ierr; 

	//void mkl_dcsrmultcsr (
	//	char *transa, MKL_INT *job, MKL_INT *sort, 
	//	MKL_INT *m, MKL_INT *n, MKL_INT *k, 
	//	double *a, MKL_INT *ja, MKL_INT *ia, 
	//	double *b, MKL_INT *jb, MKL_INT *ib, 
	//	double *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, 
	//	MKL_INT *ierr);
	   

	job = 1; 
	mkl_dcsrmultcsr        (		
		&transa, &job, &sort, 		
		&m, &n, &k, 		
		a, ja, ia, 		
		b, jb, ib, 		
		&CSR_V_values[0], &CSR_J_col_indices[0], &CSR_I_row_indices[0], &nnzmax,
		&ierr);

	job = 2; 
	nnz = CSR_I_row_indices[CSR_I_row_indices.size() - 1] - 1; //nnz = CSR_I_row_indices[m] - 1;
	nnzmax = nnz; 
	CSR_V_values.resize(nnz); 
	CSR_J_col_indices.resize(nnz);

	mkl_dcsrmultcsr        (		
		&transa, &job, &sort, 		
		&m, &n, &k, 		
		a, ja, ia, 		
		b, jb, ib, 		
		&CSR_V_values[0], &CSR_J_col_indices[0], &CSR_I_row_indices[0], &nnzmax,
		&ierr);
	
	if (transa == 'T') {
		rows = A_in.cols; 
		cols = B_in.cols; 
	} else {
		rows = A_in.rows; 
		cols = B_in.cols; 
	}
	
	type ='G'; 

}


void SparseMatrix::MatMatSorted(SparseMatrix & A_in, char MatA_T_for_transpose_N_for_non_transpose, SparseMatrix & B_in) {
	// THIS := op(A)*B

	char transa = MatA_T_for_transpose_N_for_non_transpose; 

	int job; 

	int sort = 7;	   // 3	no no no

	int m = A_in.rows; // Number of rows of the matrix A.
	int n = A_in.cols; // Number of columns of the matrix A.
	int k = B_in.cols; // Number of columns of the matrix B.

	SEQ_VECTOR<int>().swap( CSR_I_row_indices ); 
	SEQ_VECTOR<int>().swap( CSR_J_col_indices ); 
	SEQ_VECTOR<double>().swap( CSR_V_values );

	if (transa == 'T')
		CSR_I_row_indices.resize( n + 1 );
	else 
		CSR_I_row_indices.resize( m + 1 );

	CSR_J_col_indices.resize(1);
	CSR_V_values.resize(1); 

	double * a  = &A_in.CSR_V_values[0]; 
	int    * ia = &A_in.CSR_I_row_indices[0]; 
	int    * ja = &A_in.CSR_J_col_indices[0];

	double * b  = &B_in.CSR_V_values[0]; 
	int    * ib = &B_in.CSR_I_row_indices[0]; 
	int    * jb = &B_in.CSR_J_col_indices[0]; 

	int nnzmax = 1;  

	int ierr; 

	//void mkl_dcsrmultcsr (
	//	char *transa, MKL_INT *job, MKL_INT *sort, 
	//	MKL_INT *m, MKL_INT *n, MKL_INT *k, 
	//	double *a, MKL_INT *ja, MKL_INT *ia, 
	//	double *b, MKL_INT *jb, MKL_INT *ib, 
	//	double *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, 
	//	MKL_INT *ierr);


	job = 1; 
	mkl_dcsrmultcsr        (		
		&transa, &job, &sort, 		
		&m, &n, &k, 		
		a, ja, ia, 		
		b, jb, ib, 		
		&CSR_V_values[0], &CSR_J_col_indices[0], &CSR_I_row_indices[0], &nnzmax,
		&ierr);

	job = 2; 
	nnz = CSR_I_row_indices[CSR_I_row_indices.size() - 1] - 1; //nnz = CSR_I_row_indices[m] - 1;
	nnzmax = nnz; 
	CSR_V_values.resize(nnz); 
	CSR_J_col_indices.resize(nnz);

	mkl_dcsrmultcsr        (		
		&transa, &job, &sort, 		
		&m, &n, &k, 		
		a, ja, ia, 		
		b, jb, ib, 		
		&CSR_V_values[0], &CSR_J_col_indices[0], &CSR_I_row_indices[0], &nnzmax,
		&ierr);

	if (transa == 'T') {
		rows = A_in.cols; 
		cols = B_in.cols; 
	} else {
		rows = A_in.rows; 
		cols = B_in.cols; 
	}

	type ='G'; 

}



void SparseMatrix::MatAdd(SparseMatrix & A_in, SparseMatrix & B_in, char MatB_T_for_transpose_N_for_non_transpose, double beta) {
	//C := A+beta*op(B)

	char transa = MatB_T_for_transpose_N_for_non_transpose; 

	int job; 
	int sort = 3; // 3	yes	yes	yes

	int m = A_in.rows; // Number of rows of the matrix A.
	int n = A_in.cols; // Number of columns of the matrix A.

	int nnzmax; 
	int ierr; 

	double * a  = &A_in.CSR_V_values[0]; 
	int    * ia = &A_in.CSR_I_row_indices[0]; 
	int    * ja = &A_in.CSR_J_col_indices[0];

	double * b  = &B_in.CSR_V_values[0]; 
	int    * ib = &B_in.CSR_I_row_indices[0]; 
	int    * jb = &B_in.CSR_J_col_indices[0]; 

	SEQ_VECTOR<int>().swap( CSR_I_row_indices ); 
	SEQ_VECTOR<int>().swap( CSR_J_col_indices ); 
	SEQ_VECTOR<double>().swap( CSR_V_values );

	CSR_I_row_indices.resize( m + 1 );
	CSR_J_col_indices.resize(1);
	CSR_V_values.resize(1);

	//void mkl_dcsradd (
	//	char *transa, MKL_INT *job, MKL_INT *sort, 
	//	MKL_INT *m, MKL_INT *n, 
	//	double *a, MKL_INT *ja, MKL_INT *ia, 
	//	double *beta, double *b, MKL_INT *jb, MKL_INT *ib, 
	//	double *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, 
	//	MKL_INT *ierr);

	job	= 1;
	mkl_dcsradd (
		&transa, &job, &sort, 
		&m, &n, 
		a, ja, ia, 
		&beta, b, jb, ib, 
		&CSR_V_values[0], &CSR_J_col_indices[0], &CSR_I_row_indices[0], &nnzmax, 
		&ierr);

	job = 2;
	nnz = CSR_I_row_indices[m] - 1; 
	nnzmax = nnz; 
	CSR_V_values.resize(nnz); 
	CSR_J_col_indices.resize(nnz);

	mkl_dcsradd (
		&transa, &job, &sort, 
		&m, &n, 
		a, ja, ia, 
		&beta, b, jb, ib, 
		&CSR_V_values[0], &CSR_J_col_indices[0], &CSR_I_row_indices[0], &nnzmax, 
		&ierr);

	rows = A_in.rows; 
	cols = A_in.cols; 
	type = 'G'; 
}

// AM -start: ---------------------------------------------------------------------------
//
void SparseMatrix::spmv_(SparseMatrix & A, double *x, double *Ax){
  int nA = A.cols;
  int offset = A.CSR_I_row_indices[0] ? 1 : 0;
  memset(Ax,0,nA * sizeof(double));
  for (int i = 0; i < nA ; i++) {
    for (int j = A.CSR_I_row_indices[i];j<A.CSR_I_row_indices[i+1];j++) {
      Ax[i] += CSR_V_values[j-offset] * x[CSR_J_col_indices[j-offset]-offset];
      if (j > CSR_I_row_indices[i]) {
        Ax[CSR_J_col_indices[j-offset]-offset] +=
            CSR_V_values[j-offset] * x[CSR_J_col_indices[CSR_I_row_indices[i]-offset]-offset];
      }
    }
  }
}


void SparseMatrix::getSubDiagBlockmatrix(SparseMatrix & A_in, SparseMatrix & A_out, int i_start, int size_rr){
// 
// Function 'getSubDiagBlockmatrix' returns the diagonal block A_in(r,r) from original A_in,
// where r = { i_start , i_start+1 , i_start+2 , ... , istart + size_rr - 1 }
//
//
// rev. 2015-10-10 (A.M.)
//
// step 1: getting nnz of submatrix
  int nnz_new=0;
  int offset = A_in.CSR_I_row_indices[0] ? 1 : 0;
  printf("\toffset = %d\n",offset);
  for (int i = 0;i<size_rr;i++){
    for (int j = A_in.CSR_I_row_indices[i+i_start];j<A_in.CSR_I_row_indices[i+i_start+1];j++){
      if ((A_in.CSR_J_col_indices[j-offset]-offset)>=i_start && 
                      (A_in.CSR_J_col_indices[j-offset]-offset)<(i_start+size_rr)){
        nnz_new++;
      }
    }
  }
// step 2: allocation 1d arrays
  A_out.CSR_V_values.resize(nnz_new);
  A_out.CSR_J_col_indices.resize(nnz_new);
  A_out.CSR_I_row_indices.resize(size_rr+1);
  A_out.rows=size_rr;
  A_out.cols=size_rr;
  A_out.nnz=nnz_new;
	A_out.type = 'S'; 
// step 3: filling 1d arrays
  int ijcnt=0;
  A_out.CSR_I_row_indices[0]=offset;
  for (int i = 0;i<size_rr;i++){
    for (int j = A_in.CSR_I_row_indices[i+i_start];j<A_in.CSR_I_row_indices[i+i_start+1];j++){
      if ((A_in.CSR_J_col_indices[j-offset]-offset)>=i_start && 
                    (A_in.CSR_J_col_indices[j-offset]-offset)<(i_start+size_rr)){
        A_out.CSR_J_col_indices[ijcnt] = (A_in.CSR_J_col_indices[j-offset]) - i_start;
        A_out.CSR_V_values[ijcnt]=A_in.CSR_V_values[j-offset];
        ijcnt++;
      }
    }
    A_out.CSR_I_row_indices[i+1]=offset+ijcnt;
  }
}


void SparseMatrix::getSubBlockmatrix_rs( SparseMatrix & A_in, SparseMatrix & A_out, 
                                          int i_start, int i_size,
                                          int j_start, int j_size){
//
// Original matrix A_in is assembled from 4 submatrices
//
//      A_in = [A_in(r,r)  A_in(r,s)]
//             [A_in(s,r)  A_in(s,s)].
//
// Function 'getSubBlockmatrix_rs' returns square matrix A_in(r,s) in CSR format.
//
// rev. 2015-10-10 (A.M.)
//
// step 1: getting nnz of submatrix
  int nnz_new=0;
  int offset = A_in.CSR_I_row_indices[0] ? 1 : 0;
  printf("\toffset = %d\n",offset);
  for (int i = 0;i<i_size;i++){
    for (int j = A_in.CSR_I_row_indices[i+i_start];j<A_in.CSR_I_row_indices[i+i_start+1];j++){
      if ((A_in.CSR_J_col_indices[j-offset]-offset)>=j_start && 
                      (A_in.CSR_J_col_indices[j-offset]-offset)<(j_start+j_size)){
        nnz_new++;
      }
    }
  }
//  printf("nnz_new A_in(r,s)=%d\n",nnz_new);
// step 2: allocation 1d arrays
  A_out.CSR_V_values.resize(nnz_new);
  A_out.CSR_J_col_indices.resize(nnz_new);
  A_out.CSR_I_row_indices.resize(i_size+1);
  A_out.rows=i_size;
  A_out.cols=j_size;
  A_out.nnz=nnz_new;
	A_out.type = 'G'; 
// step 3: filling 1d arrays
  int ijcnt=0;
  A_out.CSR_I_row_indices[0]=offset;
  for (int i = 0;i<i_size;i++){
    for (int j = A_in.CSR_I_row_indices[i+i_start];j<A_in.CSR_I_row_indices[i+i_start+1];j++){
      if ((A_in.CSR_J_col_indices[j-offset]-offset)>=j_start && 
                    (A_in.CSR_J_col_indices[j-offset]-offset)<(j_start+j_size)){
        A_out.CSR_J_col_indices[ijcnt] = (A_in.CSR_J_col_indices[j-offset]) - j_start;
        A_out.CSR_V_values[ijcnt]=A_in.CSR_V_values[j-offset];
        ijcnt++;
      }
    }
    A_out.CSR_I_row_indices[i+1]=offset+ijcnt;
  }
}

void SparseMatrix::printMatCSR(char *str0){
  int offset = CSR_I_row_indices[0] ? 1 : 0;
  printf("%s = [ ...\n",str0);
  for (int i = 0;i<rows;i++){
    for (int j = CSR_I_row_indices[i];j<CSR_I_row_indices[i+1];j++){
      printf("%d %d %3.9e \n",i+1,CSR_J_col_indices[j-offset],CSR_V_values[j-offset]);
    }
  }
  printf("];%s = full(sparse(%s(:,1),%s(:,2),%s(:,3),%d,%d));\n",
                str0,str0,str0,str0,rows,cols);
  if (type=='S'){
    printf("%s=%s+%s'-diag(diag(%s));\n",str0,str0,str0,str0);
  }
}

//

void SparseMatrix::GramSchmidtOrtho(){
//
  double *w = new double [rows];
  double *R = new double [cols*cols];
  memset(R,0,(cols*cols) * sizeof(double));

  for (int j = 0;j<cols;j++){
    memcpy( w, &(dense_values[j*rows]) , sizeof( double ) * rows);
    for (int i = 0;i<j;i++){
      R[j*cols+i] = dot_e(w, &(dense_values[i*rows]),rows);
      for (int k=0;k<rows;k++){
        w[k]-=dense_values[i*rows+k]*R[j*cols+i];
      }
    }
    R[j*cols+j] = sqrt(dot_e(w,w,rows));
    for (int k=0;k<rows;k++){
      dense_values[j*rows+k] = w[k]/R[j*cols+j];
    }
  }
  delete [] w;
  delete [] R;
}


bool myfn(double i, double j) { return fabs(i)<=fabs(j); }

void SparseMatrix::getNullPivots(SEQ_VECTOR <int> & null_pivots){
	SEQ_VECTOR <double> N(dense_values);
  int nEl = rows*cols;
  std::vector <double>::iterator  it;
  int I,J,K,colInd,rowInd;
  double *tmpV = new double[rows];
  double pivot;
  int tmp_int;
  int *_nul_piv = new int[rows];
  for (int i = 0;i<rows;i++) _nul_piv[i]=i;

//TODO Ask about to the efficiency of next 2 lines.
  auto ij= [&]( int ii, int jj ) -> int 
   { return ii + rows*jj; };
 // 
  for (int j=0;j<cols;j++){
    it = std::max_element(N.begin(),N.end()-j*rows,myfn);
    I = it - N.begin();
    colInd = I/rows;
    rowInd = I-colInd*rows;
    for (int k=0;k<cols-j;k++){
      tmpV[k] = N[ij(rows-1-j,k)];
      N[ij(rows-1-j,k)] = N[ij(rowInd,k)];
      N[ij(rowInd,k)]= tmpV[k];
    }
    tmp_int = _nul_piv[rowInd];
    _nul_piv[rowInd] = _nul_piv[rows-1-j];
    _nul_piv[rows-1-j] = tmp_int;
    memcpy( tmpV, &(N[ij(0,cols-1-j)]) , sizeof( double ) * rows);
    memcpy( &(N[ij(0,cols-1-j)]), &(N[ij(0,colInd)]) , sizeof( double ) * rows);
    memcpy( &(N[ij(0,colInd)]),tmpV , sizeof( double ) * rows);
    pivot = N[ij(rows-1-j,cols-1-j)];
    printf("pivot = %3.9e \n",pivot);
    for (int J=0;J<cols-j-1;J++){
      for (int I=0;I<rows-j;I++){
        N[ij(I,J)] -= N[ij(I,cols-1-j)]*N[ij(rows-1-j,J)]/pivot;
      }
    }
  }  
// 
  printf("\n");
  for (int i = 0;i<cols;i++){
    null_pivots.push_back(_nul_piv[rows-1-i]+1);
  }
  sort(null_pivots.begin(),null_pivots.end());
//
  delete [] _nul_piv;
  delete [] tmpV;
//
}
//
double SparseMatrix::MatCondNumb( SparseMatrix & A_in, char *str0, int plot_n_first_n_last_eigenvalues){
//
  bool plot_a_and_b_defines_tridiag=false;
  int nA = A_in.rows;
  int nMax = 200; // size of tridiagonal matrix 
  //int nEigToplot = 10;
  double *s = new double[nA];
  double *s_bef = new double[nA];
  double *As = new double[nA];
  double *r = new double[nA];
  double tmp_a,alpha, beta = 1.0, beta_bef;
  double estim_cond;
//
//  if (nMax>nA); nMax=nA;
  nMax = nMax > nA ? nA : nMax;
  double *alphaVec = new double[nMax];
  double *betaVec  = new double[nMax];
  int cnt = 0;
//
  memset(s,0,nA * sizeof(double));
  for (int i = 0 ; i < nA; i++){ r[i] = i ; }
  tmp_a = sqrt(dot_e(r,r,nA));
  for (int i = 0 ; i < nA; i++){ r[i] /=  tmp_a; }
//
  for (int i = 0; i < nMax ; i++){
    memcpy( s_bef, s , sizeof( double ) * nA);
    beta_bef=beta;
    memcpy( s, r , sizeof( double ) * nA);
    for (int j =  0;j < nA; j++){
      s[j]/=beta;
    }
    spmv_(A_in,s,As);
    alpha = dot_e(s,As,nA);
    for (int j =  0;j < nA; j++){
      r[j]=As[j] - s[j]*alpha - s_bef[j]*beta;
    }
//
    beta = sqrt(dot_e(r,r,nA));
    alphaVec[i] = alpha;
    betaVec[i]  = beta;
//
    cnt++;
    if ( fabs(beta/beta_bef) < 1e-4 ){
      break; 
    }
  }
//
  if (plot_a_and_b_defines_tridiag){
    printf("\n alpha beta \n");
    for (int i = 0 ; i < cnt; i++){
      printf("%3.8e %3.8e\n",alphaVec[i],betaVec[i]);
    }
  }
  char JOBZ = 'N';
  double *Z = new double[cnt]; 
  MKL_INT info;
  MKL_INT ldz = cnt;
  info = LAPACKE_dstev(LAPACK_ROW_MAJOR, JOBZ, cnt, alphaVec, betaVec, Z, ldz);
  estim_cond=fabs(alphaVec[cnt-1]/alphaVec[0]);
  if (plot_n_first_n_last_eigenvalues>0){
    printf("cond(%s) = %3.15e\tit: %d\n",str0,estim_cond,cnt);
  }

  if (plot_n_first_n_last_eigenvalues>0){
    printf("eigenvals of %s d{1:%d} and d{%d:%d}\n",
          str0,plot_n_first_n_last_eigenvalues,cnt-plot_n_first_n_last_eigenvalues+2,cnt);
    for (int i = 0 ; i < cnt; i++){
      if (i < plot_n_first_n_last_eigenvalues || i > cnt-plot_n_first_n_last_eigenvalues){
        printf("%5d:  %3.8e \n",i+1, alphaVec[i]);
      }
    }
  }
//
  delete [] s;
  delete [] s_bef;
  delete [] As;
  delete [] r;
  delete [] alphaVec;
  delete [] betaVec;
  delete [] Z;

  return estim_cond;
//
}

double SparseMatrix::dot_e(double *x, double *y, int n){
  double dot_xy = 0.0;
  for (int i = 0; i< n; i++){
    dot_xy+=x[i]*y[i];
  }
  return dot_xy;
}

// AM -end:   ---------------------------------------------------------------------------



void SparseMatrix::MatAddInPlace(SparseMatrix & B_in, char MatB_T_for_transpose_N_for_non_transpose, double beta) {
	//C := A+beta*op(B)

	char transa = MatB_T_for_transpose_N_for_non_transpose; 

	// if this matrix is empty then we copy the input matrix 
	if (nnz == 0 && transa == 'N' && beta == 1.0) { // POZOR - what if we need to copy a transpose of the matrix 
		
		cols = B_in.cols;
		rows = B_in.rows;
		nnz  = B_in.nnz;
		type = B_in.type; 

		CSR_I_row_indices = B_in.CSR_I_row_indices;
		CSR_J_col_indices = B_in.CSR_J_col_indices;
		CSR_V_values      = B_in.CSR_V_values; 
		
		I_row_indices	  = B_in.I_row_indices;
		J_col_indices     = B_in.J_col_indices;
		V_values		  = B_in.V_values; 

		dense_values	  = B_in.dense_values; 

		return; 
	}


	if (nnz == 0 && transa == 'T' && beta == 1.0) { 
		cout << "Error in 'SparseMatrix::MatAddInPlace' - not implemented - " << "beta = " << beta << " Trans = " << transa << endl; 
		return; 
	}

	if (nnz == 0 && beta != 1.0) { 
		cout << "Error in 'SparseMatrix::MatAddInPlace' - not implemented - " << "beta = " << beta << " Trans = " << transa << endl; 
		return; 
	}


	int job; 
	int sort = 3; // 3	yes	yes	yes

	int m = rows; // Number of rows of the matrix A.
	int n = cols; // Number of columns of the matrix A.

	int nnzmax; 
	int ierr; 

	double * a  = &CSR_V_values[0]; 
	int    * ia = &CSR_I_row_indices[0]; 
	int    * ja = &CSR_J_col_indices[0];

	double * b  = &B_in.CSR_V_values[0]; 
	int    * ib = &B_in.CSR_I_row_indices[0]; 
	int    * jb = &B_in.CSR_J_col_indices[0]; 

	SEQ_VECTOR<int>		t_CSR_I_row_indices;	t_CSR_I_row_indices.resize( m + 1 );
	SEQ_VECTOR<int>		t_CSR_J_col_indices;	t_CSR_J_col_indices.resize(1);
	SEQ_VECTOR<double>	t_CSR_V_values;			t_CSR_V_values.resize(1); 

	//void mkl_dcsradd (
	//	char *transa, MKL_INT *job, MKL_INT *sort, 
	//	MKL_INT *m, MKL_INT *n, 
	//	double *a, MKL_INT *ja, MKL_INT *ia, 
	//	double *beta, double *b, MKL_INT *jb, MKL_INT *ib, 
	//	double *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, 
	//	MKL_INT *ierr);

	job	= 1;
	mkl_dcsradd (
		&transa, &job, &sort, 
		&m, &n, 
		a, ja, ia, 
		&beta, b, jb, ib, 
		&t_CSR_V_values[0], &t_CSR_J_col_indices[0], &t_CSR_I_row_indices[0], &nnzmax, 
		&ierr);

	job = 2;
	nnz = t_CSR_I_row_indices[m] - 1; 
	nnzmax = nnz; 
	t_CSR_V_values.resize(nnz); 
	t_CSR_J_col_indices.resize(nnz);

	mkl_dcsradd (
		&transa, &job, &sort, 
		&m, &n, 
		a, ja, ia, 
		&beta, b, jb, ib, 
		&t_CSR_V_values[0], &t_CSR_J_col_indices[0], &t_CSR_I_row_indices[0], &nnzmax, 
		&ierr);
	
	CSR_I_row_indices.swap(t_CSR_I_row_indices);
	CSR_J_col_indices.swap(t_CSR_J_col_indices);
	CSR_V_values.     swap(t_CSR_V_values); 

}

void SparseMatrix::MatScale(double alpha) {
	for (int i = 0; i < CSR_V_values.size(); i++) {
		CSR_V_values[i] = alpha * CSR_V_values[i]; 
	}
}

//void SparseMatrix::MatTranspose() {
//	MatTranspose(1.0);
//}
//
//void SparseMatrix::MatTranspose(double beta) {
//
//	SparseMatrix Temp; 
//
//	MatTranspose(Temp,beta);
//
//	rows = Temp.rows;
//	cols = Temp.cols;
//	nnz  = Temp.nnz; 
//	type = Temp.type; 
//
//	CSR_I_row_indices.clear();
//	CSR_J_col_indices.clear();
//	CSR_V_values.clear(); 
//
//	CSR_I_row_indices	= Temp.CSR_I_row_indices;
//	CSR_J_col_indices	= Temp.CSR_J_col_indices;
//	CSR_V_values		= Temp.CSR_V_values; 
//
//	Temp.Clear(); 
// 
//}
//
//void SparseMatrix::MatTranspose(SparseMatrix & A_out) {
//	MatTranspose(A_out,1.0);
//}
//
//void SparseMatrix::MatTranspose(SparseMatrix & A_out, double beta) {
//
//	char transa = 'T'; 
//	int job		=  1;  
//	int sort	=  3; // 3	yes	yes	yes
//	int nnzmax	=  1; 
//
//	int m = this->cols; // Number of rows of the matrix A.
//	int n = this->rows; // Number of columns of the matrix A.
//
//	int ierr; 
//
//	// Create an empty matrix with 1 element equal to zero 
//	SparseMatrix T; 
//	T.cols = rows;
//	T.rows = cols; 
//	T.nnz = 1; 
//	T.I_row_indices.push_back(1);
//	T.J_col_indices.push_back(1);
//	T.V_values.push_back(0);
//	T.type='G';
//	T.ConvertToCSR();
//
//	// Output matrix 
//
//	//A_out.CSR_I_row_indices.reserve (m + 10);
//
//	A_out.CSR_I_row_indices.resize( m + 1 ); 
//	A_out.CSR_J_col_indices.resize( 1 );
//	A_out.CSR_V_values.resize( 1 ); 
//
//	//void mkl_dcsradd (
//	//	char *transa, MKL_INT *job, MKL_INT *sort, 
//	//	MKL_INT *m, MKL_INT *n, 
//	//	double *a, MKL_INT *ja, MKL_INT *ia, 
//	//	double *beta, double *b, MKL_INT *jb, MKL_INT *ib, 
//	//	double *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, 
//	//	MKL_INT *ierr);
//
//	job	= 1;
//	mkl_dcsradd (
//		&transa, &job, &sort, 
//		&m, &n, 
//		&T.CSR_V_values[0],     &T.CSR_J_col_indices[0],     &T.CSR_I_row_indices[0],  //a, ja, ia, 
//		&beta,  &CSR_V_values[0],       &CSR_J_col_indices[0],       &CSR_I_row_indices[0],  //&beta, b, jb, ib, 
//		&A_out.CSR_V_values[0], &A_out.CSR_J_col_indices[0], &A_out.CSR_I_row_indices[0], &nnzmax,  
//		&ierr);
//
//	nnzmax = A_out.CSR_I_row_indices[m] - 1;
//	A_out.CSR_J_col_indices.resize( nnzmax );
//	A_out.CSR_V_values.resize( nnzmax ); 
//
//	job = 2; 
//	mkl_dcsradd (
//		&transa, &job, &sort, 
//		&m, &n, 
//		&T.CSR_V_values[0],     &T.CSR_J_col_indices[0],     &T.CSR_I_row_indices[0],  //a, ja, ia, 
//		&beta,  &CSR_V_values[0],       &CSR_J_col_indices[0],       &CSR_I_row_indices[0],  //&beta, b, jb, ib, 
//		&A_out.CSR_V_values[0], &A_out.CSR_J_col_indices[0], &A_out.CSR_I_row_indices[0], &nnzmax,   
//		&ierr);
//
//	A_out.rows = cols;
//	A_out.cols = rows;
//	A_out.nnz  = nnzmax; 
//	A_out.type = type; 
//
//}


void SparseMatrix::MatTranspose(double beta) {
	MatTranspose();
	MatScale(beta);
}

void SparseMatrix::MatTranspose(SparseMatrix & A_out, double beta) {
	MatTranspose(A_out);
	A_out.MatScale(beta);
}

void SparseMatrix::MatTranspose(SparseMatrix & A_out) {

	int job [] = { 0, 1, 1, 0, 0, 1 };
	int info;
	int m;

	A_out.cols = rows; 
	A_out.rows = cols;

	int row_size_backup = CSR_I_row_indices.size();
	if (cols > rows) {
		CSR_I_row_indices.resize( cols+1,  CSR_I_row_indices[CSR_I_row_indices.size()-1] );
		m = cols; 
	} else {
		m = rows;
	}

	SEQ_VECTOR<int>().swap( A_out.CSR_I_row_indices ); 
	SEQ_VECTOR<int>().swap( A_out.CSR_J_col_indices ); 
	SEQ_VECTOR<double>().swap( A_out.CSR_V_values );

	A_out.CSR_I_row_indices.resize(m + 1);
	A_out.CSR_J_col_indices.resize(nnz);
	A_out.CSR_V_values.		resize(nnz);

	//void mkl_dcsrcsc(int *job, int *m, double *acsr, int *ja, int *ia, double *acsc, int *ja1, int *ia1, int *info);
	//mkl_dcsrcsc( &job[0], &m, &CSR_V_values[0], &CSR_J_col_indices[0], &CSR_I_row_indices[0], &V_values[0], &J_col_indices[0], &I_row_indices[0], &info);

	mkl_dcsrcsc( &job[0], &m, 
		&CSR_V_values[0], &CSR_J_col_indices[0], &CSR_I_row_indices[0], 
		&A_out.CSR_V_values[0], &A_out.CSR_J_col_indices[0], &A_out.CSR_I_row_indices[0], 
		&info);

	if (cols > rows) {
		CSR_I_row_indices.resize(row_size_backup);
		SEQ_VECTOR<int> tmp; 
		tmp = CSR_I_row_indices; 
		CSR_I_row_indices.swap(tmp); 
	} else {
		A_out.CSR_I_row_indices.resize(A_out.rows + 1);
		SEQ_VECTOR<int> tmp; 
		tmp = A_out.CSR_I_row_indices; 
		A_out.CSR_I_row_indices.swap(tmp);
	}

	A_out.nnz  = nnz;
	A_out.type = type;

}

void SparseMatrix::MatTranspose() {

	int job [] = { 0, 1, 1, 0, 0, 1 };
	int info;
	int m;
	int row_size_backup = CSR_I_row_indices.size();

	SEQ_VECTOR <int> tCSR_I_row_indices;
	SEQ_VECTOR <int> tCSR_J_col_indices;
	SEQ_VECTOR <double> tCSR_V_values;

	if (cols > rows) {
		CSR_I_row_indices.resize( cols+1, CSR_I_row_indices[CSR_I_row_indices.size()-1] );
		m = cols; cols = rows; rows = m;
	} else {
		m = rows; rows = cols; cols = m; 
	}

	tCSR_I_row_indices.resize(m + 1);
	tCSR_J_col_indices.resize(nnz);
	tCSR_V_values.resize(nnz);

	//void mkl_dcsrcsc(int *job, int *m, double *acsr, int *ja, int *ia, double *acsc, int *ja1, int *ia1, int *info);
	//mkl_dcsrcsc( &job[0], &m, &CSR_V_values[0], &CSR_J_col_indices[0], &CSR_I_row_indices[0], &V_values[0], &J_col_indices[0], &I_row_indices[0], &info);

	mkl_dcsrcsc( &job[0], &m, 
		&CSR_V_values[0], &CSR_J_col_indices[0], &CSR_I_row_indices[0], 
		&tCSR_V_values[0], &tCSR_J_col_indices[0], &tCSR_I_row_indices[0], 
		&info);

	if (cols > rows) {
		tCSR_I_row_indices.resize(rows + 1);
	} //else {
		//CSR_I_row_indices.resize(row_size_backup); // neni nutne, stejne se prepise
	//}

	SEQ_VECTOR<int>().swap( CSR_I_row_indices ); 
	CSR_I_row_indices = tCSR_I_row_indices; //.swap(tCSR_I_row_indices);
	CSR_J_col_indices.swap(tCSR_J_col_indices);
	CSR_V_values     .swap(tCSR_V_values);

}

void SparseMatrix::MatTransposeCOO() {

	I_row_indices.swap(J_col_indices);
	int tmp = rows;
	rows = cols; 
	cols = tmp;  

	//sortInCOO();  // TODO: musi se zpatky povolit hned jak se opravit funkce sort v COO

}


void SparseMatrix::RemoveLower() {

	SEQ_VECTOR <int> t_CSR_I_row_indices;
	SEQ_VECTOR <int> t_CSR_J_col_indices;
	SEQ_VECTOR <double> t_CSR_V_values; 
	int l_nnz = 0; 

	for (int row = 0; row < CSR_I_row_indices.size() - 1; row++) {
		t_CSR_I_row_indices.push_back(l_nnz+1);
		int cols_in_row = CSR_I_row_indices[row+1] - CSR_I_row_indices[row]; 

		for (int col = 0; col <cols_in_row; col++) {
			int i = CSR_I_row_indices[row] - 1 + col;
			if ( CSR_J_col_indices[i] > row) {
				t_CSR_J_col_indices.push_back(CSR_J_col_indices[i]);
				t_CSR_V_values.push_back(CSR_V_values[i]);
				l_nnz++; 
			}
		}

	}

	t_CSR_I_row_indices.push_back(l_nnz+1);

	nnz = l_nnz; 
	type = 'S';

	//CSR_I_row_indices = t_CSR_I_row_indices;
	//CSR_J_col_indices = t_CSR_J_col_indices;
	//CSR_V_values = t_CSR_V_values; 
	
	// musis se otestovat verze se .swap()
	CSR_I_row_indices.swap( t_CSR_I_row_indices );
	CSR_J_col_indices.swap( t_CSR_J_col_indices );
	CSR_V_values.swap( t_CSR_V_values ); 

}

double SparseMatrix::GetMeanOfDiagonalOfSymmetricMatrix() {

	double sum = 0; 
	int count = 0; 

	for (int i = 0; i < CSR_I_row_indices.size() - 1; i++) {
		double val = CSR_V_values[ CSR_I_row_indices[i] - 1 ];
		sum = sum + val; 
		count++;
	}

	return sum/count; 
}

double SparseMatrix::GetMaxOfDiagonalOfSymmetricMatrix() {
	double vmax = 0; 

	for (int i = 0; i < CSR_I_row_indices.size() - 1; i++) {

		if ( vmax < CSR_V_values[ CSR_I_row_indices[i] - 1 ] )
			vmax = CSR_V_values[ CSR_I_row_indices[i] - 1 ]; 

	}

	return vmax; 
}


void SparseMatrix::SetDiagonalOfSymmetricMatrix( double val ) {
	for (int i = 0; i < CSR_I_row_indices.size() - 1; i++) {
			CSR_V_values[ CSR_I_row_indices[i] - 1 ] = val;
	}
}


void SparseMatrix::MatAppend(SparseMatrix & A) {

	if (nnz == 0 && rows == 0 && cols == 0) { // this matrix is empty 
		rows = A.rows;
		cols = A.cols;
		nnz = A.nnz;
		type = A.type; 

		CSR_I_row_indices = A.CSR_I_row_indices;
		CSR_J_col_indices = A.CSR_J_col_indices;
		CSR_V_values	  = A.CSR_V_values; 

	}else {
		// Just append the arrays 	
		CSR_J_col_indices.insert(CSR_J_col_indices.end(), A.CSR_J_col_indices.begin(), A.CSR_J_col_indices.end() ); 
		CSR_V_values.insert(CSR_V_values.end(), A.CSR_V_values.begin(), A.CSR_V_values.end() ); 

		//copy(CSR_J_col_indices.begin(), A.CSR_J_col_indices.begin(), A.CSR_J_col_indices.end());
		//copy(CSR_V_values.begin(),      A.CSR_V_values.begin(),      A.CSR_V_values.end());

		int last_row = CSR_I_row_indices[CSR_I_row_indices.size()-1]; 

		for (int i = 1; i < A.CSR_I_row_indices.size(); i++)
			CSR_I_row_indices.push_back(last_row + A.CSR_I_row_indices[i]-1);	

		rows = rows + A.rows; 
		nnz  = nnz + A.nnz; 
		cols = cols; 
		type = 'G'; 
	}
}


void SparseMatrix::CreateMatFromRowsFromMatrix(SparseMatrix & A_in, SEQ_VECTOR <int> & rows_to_add) {

	int old_index  = 0;
	int next_index = 0; 
	int row_fill   = 1;

	rows = A_in.rows;
	cols = A_in.cols;
	type = A_in.type; 

	CSR_I_row_indices.resize( rows + 1 );

	for (int i = 0; i < rows_to_add.size(); i++) {
		
		old_index  = next_index; 
		next_index = rows_to_add[i]; 

		fill(CSR_I_row_indices.begin() + old_index, CSR_I_row_indices.begin() + next_index, row_fill);
		
		int A_in_start_index = A_in.CSR_I_row_indices[rows_to_add[i] - 1 ] - 1 ; 
		int A_in_end_index   = A_in.CSR_I_row_indices[rows_to_add[i] + 1 - 1] - 1 ; 

		CSR_J_col_indices.insert(CSR_J_col_indices.end(), A_in.CSR_J_col_indices.begin() + A_in_start_index, A_in.CSR_J_col_indices.begin() + A_in_end_index ); 
		CSR_V_values.     insert(CSR_V_values.end(),      A_in.CSR_V_values.     begin() + A_in_start_index, A_in.CSR_V_values.     begin() + A_in_end_index ); 
		row_fill = 1 + CSR_J_col_indices.size(); 		

	}

	
	fill(CSR_I_row_indices.begin() + next_index, CSR_I_row_indices.begin() + rows + 1, row_fill);

	nnz = CSR_V_values.size(); 

}



int SparseMatrix::MatCompare(SparseMatrix & A) {
	int res = 0; 

	if (this->cols == A.cols && this->rows==A.rows && this->nnz == A.nnz && this->type == A.type ) {

		int tmp1 = 0; 
		int tmp2 = 0; 

		for (int i = 0; i < CSR_I_row_indices.size(); i++)
			if (CSR_I_row_indices[i] != A.CSR_I_row_indices[i])
				tmp1=1; 

		for (int i = 0; i < CSR_J_col_indices.size(); i++)
			if (CSR_J_col_indices[i] != A.CSR_J_col_indices[i])
				tmp1=1; 

		for (int i = 0; i < CSR_V_values.size(); i++)
			if (CSR_V_values[i] != A.CSR_V_values[i])
				tmp2=1; 

		res = 1000 * tmp1 + tmp2; 

	} else {
		res = -1; 
	}

	return res; 
}

int SparseMatrix::MatCompareCOO(SparseMatrix & A) {
	int res = 0; 

	if (this->cols == A.cols && this->rows==A.rows && this->nnz == A.nnz && this->type == A.type ) {

		int tmp1 = 0; 
		int tmp2 = 0; 
		int tmp3 = 0; 

		for (int i = 0; i < I_row_indices.size(); i++)
			if (I_row_indices[i] != A.I_row_indices[i])
				tmp1=1; 

		for (int i = 0; i < J_col_indices.size(); i++)
			if (J_col_indices[i] != A.J_col_indices[i])
				tmp2=1; 

		for (int i = 0; i < V_values.size(); i++)
			if (V_values[i] != A.V_values[i])
				tmp3=1; 

		res = 100 * tmp1 + 10 * tmp2 + tmp3; 

	} else {
		res = -1; 
	}

	return res; 
}


void SparseMatrix::CreateEye(int size) {

	for (int i = 0; i< size; i++) {
		J_col_indices.push_back(i+1);
		I_row_indices.push_back(i+1);
		V_values.push_back(1.0);
	}

	rows = size;
	cols = size;
	nnz  = size;
	type = 'G';

	ConvertToCSR(); 

}


void SparseMatrix::CreateEye(int size, double value, int offset_row, int offset_col) {

	for (int i = 0; i< size; i++) {
		J_col_indices.push_back(offset_col + i+1);
		I_row_indices.push_back(offset_row + i+1);
		V_values.push_back( value );
	}

	rows = offset_row + size;
	cols = offset_col + size;
	nnz  = size;
	type = 'G';

	ConvertToCSR();

}



void SparseMatrix::TestEye(int size) {

	for (int i = 0; i< size; i++) {
		J_col_indices.push_back(i+1);
		I_row_indices.push_back(i+1);
		V_values.     push_back(i+1.0);
	}

	rows = size;
	cols = size;
	nnz  = size;
	type = 'G';

	ConvertToCSR(); 

}

void SparseMatrix::TestMatRow(int size, int row_index) {

	for (int i = 0; i< size; i++) {
		J_col_indices.push_back(i+1);
		I_row_indices.push_back(row_index);
		V_values.     push_back(i+1.0);
	}

	rows = size;
	cols = size;
	nnz  = size;
	type = 'G';

	ConvertToCSR(); 

}



// **** END - Sparse Matrix CLASS ************************************
// *******************************************************************

void SparseMatrix::MatMatT(SparseMatrix & A_in, SparseMatrix & B_in) {
	
	//SEQ_VECTOR < SEQ_VECTOR < int    > > GGt_J (A_in.rows, SEQ_VECTOR < int    > () ); 
	//SEQ_VECTOR < SEQ_VECTOR < double > > GGt_V (A_in.rows, SEQ_VECTOR < double > () ); 
	
	SEQ_VECTOR<int>()   .swap( CSR_I_row_indices ); 
	SEQ_VECTOR<int>()   .swap( CSR_J_col_indices ); 
	SEQ_VECTOR<double>().swap( CSR_V_values );

	int glob_row_index = 0 + 1; 
	CSR_I_row_indices.push_back(glob_row_index);

	for (int i = 0; i < A_in.CSR_I_row_indices.size() - 1; i++ ) {

		int A_row_start = A_in.CSR_I_row_indices[i  ] - 1; 
		int A_row_end   = A_in.CSR_I_row_indices[i+1] - 1;

		if (A_row_start != A_row_end ) { // this row in B is NOT empty 
	
			for (int ii = 0; ii < B_in.CSR_I_row_indices.size() - 1; ii++ ) {

				int B_row_start = B_in.CSR_I_row_indices[ii  ] - 1; 
				int B_row_end   = B_in.CSR_I_row_indices[ii+1] - 1;

				if (B_row_start != B_row_end) { // this row in B is NOT empty 
					int A_ind = A_row_start; 
					int B_ind = B_row_start; 
					double C_v = 0; 
					do {

						int A_j = A_in.CSR_J_col_indices[A_ind]; 
						int B_j = B_in.CSR_J_col_indices[B_ind]; 
						if (A_j < B_j) A_ind++; 
						if (A_j > B_j) B_ind++; 
						
						if (A_j == B_j) {
							C_v += A_in.CSR_V_values[A_ind] * B_in.CSR_V_values[B_ind]; 
							A_ind++;
							B_ind++; 
						}
					} while ( A_ind < A_row_end && B_ind < B_row_end ); 

					if (C_v != 0.0) {
						//GGt_J[i].push_back(ii + 1);
						//GGt_V[i].push_back(C_v); 

						CSR_J_col_indices.push_back(ii + 1);
						CSR_V_values.push_back(C_v);
						glob_row_index++; 
					}

				}
			}
	
		}
		CSR_I_row_indices.push_back(glob_row_index); 
		int a = 10;
	}

	rows = A_in.rows;
	cols = A_in.rows;
	nnz  = CSR_V_values.size(); 
	type = 'G'; 

}
