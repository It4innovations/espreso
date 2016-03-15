#include "../generic/SparseMatrix.h"

#include "../specific/sparsesolvers.h"

using namespace espreso;

std::ostream& espreso::operator<<(std::ostream& os, const SparseMatrix &m)
{
	os << m.rows << " " << m.cols << " " << m.nnz << "\n";

	os.precision(15);

	SparseMatrix s = m;
	if (s.CSR_J_col_indices.size()) {
		s.ConvertToCOO(1);
	}

	for (size_t i = 0; i < s.nnz; i++) {
		os << s.I_row_indices[i] << " ";
		os << s.J_col_indices[i] << " ";
		os << std::scientific << s.V_values[i] << "\n";
	}
	return os;
}

#define comp_x(_a,_b,_x,_p) (_a[_x]!=_a[_p] ? _a[_x]-_a[_p] : _b[_x]-_b[_p])

static void q_sort_in(vector <eslocal>    & I_row_indices,
					  vector <eslocal>    & J_col_indices,
					  vector <double> & V_values,
					  eslocal lo, eslocal hi ) {
	eslocal h, l, p, p1, p2, t;
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


std::string SparseMatrix::SpyText()
{
	std::stringstream ss;
	SEQ_VECTOR<char> tmp (60,'-');

	for( std::SEQ_VECTOR<char>::const_iterator i = tmp.begin(); i != tmp.end(); ++i) {
		ss << *i << ' ';
	}
	ss << "\n";

	eslocal rows_coef = 1 + rows / 60;
	eslocal cols_coef = 1 + cols / 60;

	eslocal col_index = 0;
	eslocal row_index = 0;
	for (eslocal r = 0; r < rows; r = r + rows_coef) {
		eslocal row_length = 0;
		if (( r + rows_coef) < rows)
			row_length = CSR_I_row_indices[r+rows_coef] - CSR_I_row_indices[r];
		else
 			row_length = CSR_I_row_indices[rows] - CSR_I_row_indices[r];

			SEQ_VECTOR<char> tmp (60,' ');
			SEQ_VECTOR<eslocal> tmp_c (60,0);
		for (eslocal c = 0; c < row_length; c++) {
			if (CSR_V_values[col_index] != 0.0) {
				tmp_c[CSR_J_col_indices[col_index] / cols_coef]++;
			} else {
				if (tmp_c[CSR_J_col_indices[col_index] / cols_coef] == 0)
					tmp_c[CSR_J_col_indices[col_index] / cols_coef] = -1;
			}
			col_index++;
		}

		for (eslocal c = 0; c < tmp_c.size(); c++) {
			if (tmp_c[c] > 0) {
				tmp[c] = '0' + tmp_c[c] / (cols_coef * 26);
				if (tmp[c] == '0') tmp[c] = '.';
			} else {
				if (tmp_c[c] == -1)
					tmp[c] = 'x';
			}
		}

		ss << "|";
		for( std::SEQ_VECTOR<char>::const_iterator i = tmp.begin(); i != tmp.end(); ++i) {
			ss << *i << ' ';
		}

		ss << "|" << endl;
	}

	//SEQ_VECTOR<char> tmp (60,'-');

	for( std::SEQ_VECTOR<char>::const_iterator i = tmp.begin(); i != tmp.end(); ++i) {
		ss << *i << ' ';
	}

	ss << "\n";

	return ss.str();
}

void SparseMatrix::sortInCOO()
{
	q_sort_in(I_row_indices, J_col_indices,V_values, 0, nnz - 1);
}


SparseMatrix::SparseMatrix() {

	nnz  = 0;
	cols = 0;
	rows = 0;
	type = 0;

	USE_FLOAT = false;

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

	USE_FLOAT = A_in.USE_FLOAT;

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

	USE_FLOAT = false;

	eslocal offset = (A_in.rowPtrs()[0]) ? 0 : 1;
	nnz -= A_in.rowPtrs()[0];

	CSR_I_row_indices.resize(rows+1);
	CSR_J_col_indices.resize(nnz);
	CSR_V_values	 .resize(nnz);

	// Sparse CSR data
	//copy(rows, rows + n_cols + 1, K.CSR_I_row_indices.begin());
	for (eslocal i = 0; i < CSR_I_row_indices.size(); i++)
		CSR_I_row_indices[i] = A_in.rowPtrs()[i] + offset;

	//copy(cols, cols + nnz, K.CSR_J_col_indices.begin());
	for (eslocal i = 0; i < CSR_J_col_indices.size(); i++)
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

void SparseMatrix::swap ( SparseMatrix &A_in) {

	eslocal tmp;
	char ttype;

	tmp = rows; rows = A_in.rows; A_in.rows = tmp;
	tmp = cols; cols = A_in.cols; A_in.cols = tmp;
	tmp = nnz;  nnz  = A_in.nnz;  A_in.nnz  = tmp;
	ttype = type; type = A_in.type; A_in.type = ttype;

	bool tb = USE_FLOAT; USE_FLOAT = A_in.USE_FLOAT; A_in.USE_FLOAT = tb;

	// Sparse COO data
	I_row_indices.swap( A_in.I_row_indices );
	J_col_indices.swap( A_in.J_col_indices );
	V_values.swap     ( A_in.V_values );

	// Sparse CSR data
	CSR_I_row_indices.swap( A_in.CSR_I_row_indices );
	CSR_J_col_indices.swap( A_in.CSR_J_col_indices );
	CSR_V_values.swap     ( A_in.CSR_V_values );

	// Dense data
	dense_values	  .swap( A_in.dense_values );
	dense_values_fl   .swap( A_in.dense_values_fl );

	// GPU
	double * tmpp;
	tmpp = d_dense_values; d_dense_values = A_in.d_dense_values; A_in.d_dense_values = tmpp;
	tmpp = d_x_in;         d_x_in  = A_in.d_x_in;                A_in.d_x_in  = tmpp;
	tmpp = d_y_out;		   d_y_out = A_in.d_y_out;		         A_in.d_y_out = tmpp;

	float * tmppf;
	tmppf = d_dense_values_fl; d_dense_values_fl = A_in.d_dense_values_fl; A_in.d_dense_values_fl = tmppf;
	tmppf = d_x_in_fl		; d_x_in_fl         = A_in.d_x_in_fl;         A_in.d_x_in_fl         = tmppf;
	tmppf = d_y_out_fl		; d_y_out_fl        = A_in.d_y_out_fl;        A_in.d_y_out_fl        = tmppf;

#ifdef CUDA
	cublasHandle_t thandle;
	cudaStream_t   tstream;

	thandle = handle; handle = A_in.handle; A_in.handle = thandle;
	tstream = stream; stream = A_in.stream; A_in.stream = tstream;
#endif

}

SparseMatrix::SparseMatrix( const SparseCSRMatrix<eslocal> &A_in, char type_in ) {

	rows = A_in.rows();
	cols = A_in.columns();
	nnz  = A_in.rowPtrs()[rows];
	type = type_in;

	USE_FLOAT = false;

	eslocal offset = (A_in.rowPtrs()[0]) ? 0 : 1;
	nnz -= A_in.rowPtrs()[0];

	CSR_I_row_indices.resize(rows+1);
	CSR_J_col_indices.resize(nnz);
	CSR_V_values	 .resize(nnz);

	// Sparse CSR data
	//copy(rows, rows + n_cols + 1, K.CSR_I_row_indices.begin());
	for (eslocal i = 0; i < CSR_I_row_indices.size(); i++)
		CSR_I_row_indices[i] = A_in.rowPtrs()[i] + offset;

	//copy(cols, cols + nnz, K.CSR_J_col_indices.begin());
	for (eslocal i = 0; i < CSR_J_col_indices.size(); i++)
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

	USE_FLOAT = false;

	eslocal offset = A_in.indexing() ? 0 : 1;

	I_row_indices.resize(nnz);
	J_col_indices.resize(nnz);
	V_values	 .resize(nnz);

	// Sparse CSR data
	//copy(rows, rows + n_cols + 1, K.CSR_I_row_indices.begin());
	for (eslocal i = 0; i < I_row_indices.size(); i++)
		I_row_indices[i] = A_in.rowIndices()[i] + offset;

	//copy(cols, cols + nnz, K.CSR_J_col_indices.begin());
	for (eslocal i = 0; i < J_col_indices.size(); i++)
		J_col_indices[i] = A_in.columnIndices()[i] + offset;

	copy(A_in.values().begin(), A_in.values().end(), V_values.begin());

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

	USE_FLOAT = false;

	eslocal offset = A_in.indexing() ? 0 : 1;

	I_row_indices.resize(nnz);
	J_col_indices.resize(nnz);
	V_values	 .resize(nnz);

	// Sparse CSR data
	//copy(rows, rows + n_cols + 1, K.CSR_I_row_indices.begin());
	for (eslocal i = 0; i < I_row_indices.size(); i++)
		I_row_indices[i] = A_in.rowIndices()[i] + offset;

	//copy(cols, cols + nnz, K.CSR_J_col_indices.begin());
	for (eslocal i = 0; i < J_col_indices.size(); i++)
		J_col_indices[i] = A_in.columnIndices()[i] + offset;

	copy(A_in.values().begin(), A_in.values().end(), V_values.begin());

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

		USE_FLOAT = A_in.USE_FLOAT;

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

	SEQ_VECTOR<eslocal>().swap( I_row_indices );
	SEQ_VECTOR<eslocal>().swap( J_col_indices );
	SEQ_VECTOR<double>().swap( V_values );

	SEQ_VECTOR<eslocal>().swap( CSR_I_row_indices );
	SEQ_VECTOR<eslocal>().swap( CSR_J_col_indices );
	SEQ_VECTOR<double>().swap( CSR_V_values );

	SEQ_VECTOR<double>().swap( dense_values );
	SEQ_VECTOR<float>().swap( dense_values_fl );

	// GPU
	//d_dense_values = NULL;
	//d_y_out        = NULL;
	//d_x_in		   = NULL:
}

eslocal  SparseMatrix::SaveMatrixBinInCOO(string filename) {

	//ConvertToCSR(0);

	std::ofstream out (filename.c_str(), std::ios::out | std::ios::binary);

	if ( out.is_open() ) {
		char delim = ';';

		//write parameters
		out << "%% rows;cols;nnz;type" << endl;
		out << rows << ";" << cols << ";" << nnz << ";" << type << endl;

		out.write((char*)&CSR_I_row_indices[0], CSR_I_row_indices.size() * sizeof(eslocal));
		cout << endl;

		out.write((char*)&CSR_J_col_indices[0], CSR_J_col_indices.size() * sizeof(eslocal));
		cout << endl;

		out.write((char*)&CSR_V_values[0], CSR_V_values.size() * sizeof(double));
		cout << endl;

		out.close();
		return 0;

	} else {
		cout << "Matrix file " << filename << " cannot be created ! " << endl;
		return -1;
	}

}

eslocal SparseMatrix::LoadMatrixBinInCOO(string filename, char matrix_type_G_for_general_S_for_symmetric) {

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
		in.read((char*) &I_row_indices[0], nnz*sizeof(eslocal));

		J_col_indices.resize(nnz);
		in.read((char*) &J_col_indices[0], nnz*sizeof(eslocal));

		V_values.resize(nnz);
		in.read((char*) &V_values[0], nnz*sizeof(double));

		in.close();

		return 0;

	} else {

		cout << "Matrix file " << filename << " not found ! " << endl;
		return -1;

	}
}

eslocal SparseMatrix::LoadMatrixBin(string filename, char matrix_type_G_for_general_S_for_symmetric) {
	eslocal tmp = LoadMatrixBinInCOO(filename, matrix_type_G_for_general_S_for_symmetric);

	if (tmp == 0)
		ConvertToCSR( 1 );

	return tmp;
}

eslocal SparseMatrix::LoadMatrixBin(string filename, char matrix_type_G_for_general_S_for_symmetric, eslocal clearCOO_1_keep_COO_0 ) {
	eslocal tmp = LoadMatrixBinInCOO(filename, matrix_type_G_for_general_S_for_symmetric);

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
//	// Get line with I values - eslocal
//	getline(in,line);
//	stringstream ssI(line);
//	while (getline(ssI,field,delim))  // break line into comma delimitted fields
//	{
//		I_row_indices.push_back(atoi(field.c_str()));  // add each field to the 1D array
//	}
//
//
//	// Get line with J values - eslocal
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



void SparseMatrix::PrintMatSize( string Matname ) {

//	// Sparse COO data
//	I_row_indices = A_in.I_row_indices;
//	J_col_indices = A_in.J_col_indices;
//	V_values	  = A_in.V_values;
//
//	// Sparse CSR data
//	CSR_I_row_indices = A_in.CSR_I_row_indices;
//	CSR_J_col_indices = A_in.CSR_J_col_indices;
//	CSR_V_values	  = A_in.CSR_V_values;
//
//	// Dense data
//	dense_values	  = A_in.dense_values;
//	dense_values_fl   = A_in.dense_values_fl;

	eslocal dense_size = dense_values.size() * sizeof(double);
	eslocal CSR_size   = CSR_I_row_indices.size() * sizeof(eslocal) + CSR_J_col_indices.size() * sizeof(eslocal) + CSR_V_values.size() * sizeof(double);
	eslocal IJV_size   = I_row_indices.size() * sizeof(eslocal) 	+ J_col_indices.size() * sizeof(eslocal) 	  + V_values.size() * sizeof(double);

	std::cout << std::endl << "Matrix " << Matname << " sizes: "<< std::endl;
	std::cout << "DNS size: " << dense_size << " B" << std::endl;
	std::cout << "CSR size: " << CSR_size << " B"  << std::endl;
	std::cout << "IJV size: " << IJV_size << " B" << std::endl <<std::endl;
}


void SparseMatrix::ConvertToCSR( ) {
	ConvertToCSR( 1 );
}

void SparseMatrix::ConvertToCSRwithSort( eslocal clearCOO_1_keep_COO_0 ) {

	eslocal job[8];//  = {0,0,0,0, 0,0,0,0};
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

	eslocal info;

	//void mkl_dcsrcoo ( eslocal * job, eslocal * n, double *Acsr,      eslocal * AJR,          eslocal * AIR,          eslocal * nnz,  double *Acoo,  eslocal * ir,       eslocal * jc,       eslocal * info);
	mkl_dcsrcoo		   ( job,           &rows,       &CSR_V_values[0],  &CSR_J_col_indices[0],  &CSR_I_row_indices[0],  &nnz,           &V_values[0],  &I_row_indices[0],  &J_col_indices[0],  &info );

	if (clearCOO_1_keep_COO_0 == 1 ) {

		V_values.clear();
		I_row_indices.clear();
		J_col_indices.clear();

		SEQ_VECTOR<double>().swap( V_values );
		SEQ_VECTOR<eslocal>().swap( I_row_indices );
		SEQ_VECTOR<eslocal>().swap( J_col_indices );

	}

	// m	INTEGER. Number of rows of the matrix A.
	// n	INTEGER. Number of columns of the matrix A.

}


void SparseMatrix::ConvertToCSR( eslocal clearCOO_1_keep_COO_0 ) {

	eslocal job[8];//  = {0,0,0,0, 0,0,0,0};
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

	eslocal info;

	//void mkl_dcsrcoo ( eslocal * job, eslocal * n, double *Acsr,      eslocal * AJR,          eslocal * AIR,          eslocal * nnz,  double *Acoo,  eslocal * ir,       eslocal * jc,       eslocal * info);
	mkl_dcsrcoo		   ( job,           &rows,       &CSR_V_values[0],  &CSR_J_col_indices[0],  &CSR_I_row_indices[0],  &nnz,           &V_values[0],  &I_row_indices[0],  &J_col_indices[0],  &info );

	if (clearCOO_1_keep_COO_0 == 1 ) {

		V_values.clear();
		I_row_indices.clear();
		J_col_indices.clear();

		SEQ_VECTOR<double>().swap( V_values );
		SEQ_VECTOR<eslocal>().swap( I_row_indices );
		SEQ_VECTOR<eslocal>().swap( J_col_indices );

	}

	// m	INTEGER. Number of rows of the matrix A.
	// n	INTEGER. Number of columns of the matrix A.

}

void SparseMatrix::ConvertToCOO( eslocal clearCSR_1_keep_CSR_0 ) {

	eslocal job[8];//  = {0,0,0,0, 0,0,0,0};
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

	eslocal info;

	//void mkl_dcsrcoo ( eslocal * job, eslocal * n, double *Acsr,      eslocal * AJR,          eslocal * AIR,          eslocal * nnz,  double *Acoo,  eslocal * ir,       eslocal * jc,       eslocal * info);
	mkl_dcsrcoo		   ( job,           &rows,       &CSR_V_values[0],  &CSR_J_col_indices[0],  &CSR_I_row_indices[0],  &nnz,           &V_values[0],  &I_row_indices[0],  &J_col_indices[0],  &info );

	if (clearCSR_1_keep_CSR_0 == 1 ) {

		CSR_V_values.clear();
		CSR_I_row_indices.clear();
		CSR_J_col_indices.clear();

		SEQ_VECTOR<eslocal>().swap( CSR_I_row_indices );
		SEQ_VECTOR<eslocal>().swap( CSR_J_col_indices );
		SEQ_VECTOR<double>().swap( CSR_V_values );

	}

	// m	INTEGER. Number of rows of the matrix A.
	// n	INTEGER. Number of columns of the matrix A.

}


void SparseMatrix::ConvertCSRToDense( eslocal clearCSR_1_keep_CSR_0 ) {

	eslocal job[8];
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

	eslocal m		= rows;
	eslocal n		= cols;
	eslocal lda     = m;
	eslocal info    = 0;

	dense_values.resize(m * n, 0);

	// Convert matrix to dense format

	//void mkl_ddnscsr (
	//	eslocal *job,
	//	eslocal *m, eslocal *n,
	//	double *Adns, eslocal *lda,
	//	double *Acsr, eslocal *AJ, eslocal *AI,
	//	eslocal *info);

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

		SEQ_VECTOR<eslocal>().swap( CSR_I_row_indices );
		SEQ_VECTOR<eslocal>().swap( CSR_J_col_indices );
		SEQ_VECTOR<double>().swap( CSR_V_values );
	}

	if (type == 'S')
		this->RemoveLowerDense();

}

void SparseMatrix::ConvertDenseToCSR( eslocal clearDense_1_keep_Dense_0 ){

	eslocal m		= rows;
	eslocal n		= cols;
	eslocal lda     = m;
	eslocal info    = 0;

	eslocal job[8];

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
	eslocal nnzmax = CSR_I_row_indices[m] - 1; //- 1;  POZOR bez -1 polud se to pouzije ve funkci SolveMatF
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

void SparseMatrix::ConvertDenseToDenseFloat( eslocal clear_DoubleDense_1_keep_DoubleDense_0 ) {

	dense_values_fl.resize( dense_values.size() );

	for (eslocal i = 0; i < dense_values.size(); i++)
		dense_values_fl[i] = (float)dense_values[i];

	if ( clear_DoubleDense_1_keep_DoubleDense_0 == 1)
		SEQ_VECTOR<double>().swap( dense_values );

}

//void SparseMatrix::DenseMatVec(vector <double> & x_in, vector <double> & y_out) {
//
//	// void cblas_dgemv
//	//  (const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA,
//	//  const eslocal M, const eslocal N,
//	//  const double alpha, const double *A, const eslocal lda,
//	//  const double *X, const eslocal incX,
//	//  const double beta, double *Y, const eslocal incY);
//
//	// y := alpha*A*x + beta*y,
//
//	double alpha = 1.0;
//	double beta  = 0.0;
//	eslocal lda = rows;
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
	DenseMatVec(x_in, y_out, 'N', 0, 0, 0.0);
}


void SparseMatrix::DenseMatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose ) {
	DenseMatVec(x_in, y_out, T_for_transpose_N_for_not_transpose, 0, 0, 0.0);
}

void SparseMatrix::DenseMatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose, eslocal x_in_vector_start_index) {
	DenseMatVec(x_in, y_out, T_for_transpose_N_for_not_transpose, x_in_vector_start_index, 0, 0.0);
}

void SparseMatrix::DenseMatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose, eslocal x_in_vector_start_index, eslocal y_out_vector_start_index) {
	DenseMatVec(x_in, y_out, T_for_transpose_N_for_not_transpose, x_in_vector_start_index, y_out_vector_start_index, 0.0);
}

void SparseMatrix::DenseMatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose, eslocal x_in_vector_start_index, eslocal y_out_vector_start_index, double beta) {


	// void cblas_dgemv
	//  (const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA,
	//  const eslocal M, const eslocal N,
	//  const double alpha, const double *A, const eslocal lda,
	//  const double *X, const eslocal incX,
	//  const double beta, double *Y, const eslocal incY);

	// y := alpha*A*x + beta*y,

	double alpha = 1.0;
	//double beta  = 0.0;
	eslocal lda = rows;

	//char trans = 'N'; // CblasNoTrans=111,     /* trans='N' */
	// CblasTrans=112,       /* trans='T' */
	// CblasConjTrans=113};  /* trans='C' */

	if (type == 'G') {
		if ( T_for_transpose_N_for_not_transpose == 'T' )
			cblas_dgemv
				(CblasColMajor, CblasTrans,
				rows, cols,
				alpha, &dense_values[0], lda,
				&x_in[x_in_vector_start_index], 1,
				beta, &y_out[y_out_vector_start_index], 1);
		else
			cblas_dgemv
				(CblasColMajor, CblasNoTrans,
				rows, cols,
				alpha, &dense_values[0], lda,
				&x_in[x_in_vector_start_index], 1,
				beta, &y_out[y_out_vector_start_index], 1);
	} else {
		if ( T_for_transpose_N_for_not_transpose == 'T' ) {
                    std::cout << "Transposition is not supported for packed symmetric matrices" << std::endl;
                    return;
		} else {

			if (!USE_FLOAT) {
				cblas_dspmv(
						CblasColMajor, CblasUpper,
						rows,
						alpha, &dense_values[0],
						&x_in[x_in_vector_start_index], 1,
						beta, &y_out[y_out_vector_start_index], 1);
			} else {

				if (vec_fl_in.size()  < rows) vec_fl_in. resize(rows);
				if (vec_fl_out.size() < rows) vec_fl_out.resize(rows);

				for (eslocal i = 0; i < rows; i++)
					vec_fl_in[i] = (float)x_in[i + x_in_vector_start_index];

				cblas_sspmv(
						CblasColMajor, CblasUpper,
						rows,
						alpha, &dense_values_fl[0],
						&vec_fl_in[0], 1,
						beta, &vec_fl_out[0], 1);

				for (eslocal i = 0; i < rows; i++)
					y_out[i + y_out_vector_start_index] = (double)vec_fl_out[i];

				//cout << "using float " << endl;

			}
		}
	}
}


void SparseMatrix::DenseMatVecCUDA_w_Copy(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose, eslocal x_in_vector_start_index) {
#ifdef CUDA

	double *d_x_in_t, *d_y_out_t, *d_Mat;
	eslocal mat_size = rows * cols;
	eslocal lda = rows;

	cudaMalloc((void**)&d_x_in_t,  x_in.size()  * sizeof(double));
	cudaMalloc((void**)&d_y_out_t, y_out.size() * sizeof(double));
	cudaMalloc((void**)&d_Mat,   mat_size     * sizeof(double));

	// Create cublas instance
	cublasHandle_t handle;
	cublasCreate(&handle);

	// Set input matrices on device
	// cublasSetVector(eslocal n, eslocal elemSize, const void *x, eslocal incx, void *y, eslocal incy);
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
	//cublasGetVector(eslocal n, eslocal elemSize, const void *x, eslocal incx, void *y, eslocal incy)
	cublasGetVector(y_out.size(), sizeof(double), d_y_out_t, 1, &y_out[0], 1);

	cudaFree(d_x_in_t);
	cudaFree(d_y_out_t);
	cudaFree(d_Mat);
	cublasDestroy(handle);


#endif
}


void SparseMatrix::DenseMatVecCUDA_wo_Copy(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_not_transpose, eslocal x_in_vector_start_index) {
#ifdef CUDA

	eslocal mat_size = rows * cols;
	eslocal lda = rows;

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
	// cublasSetVector(eslocal n, eslocal elemSize, const void *x, eslocal incx, void *y, eslocal incy);
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
	//cublasGetVector(eslocal n, eslocal elemSize, const void *x, eslocal incx, void *y, eslocal incy)
	cublasGetVector(rows , sizeof(double), d_y_out, 1, &y_out[0], 1);

	//cudaFree(d_x_in_t);
	//cudaFree(d_y_out_t);
	//cudaFree(d_dense_values);
	//cublasDestroy(handle);

	cudaStreamSynchronize(stream);

#endif
}




void SparseMatrix::DenseMatVecCUDA_wo_Copy_start( double * x_in, double * y_out, char T_for_transpose_N_for_not_transpose, eslocal x_in_vector_start_index) {
#ifdef CUDA

	eslocal mat_size = dense_values.size(); //rows * cols;
	eslocal lda = rows;

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
	// cublasSetVector(eslocal n, eslocal elemSize, const void *x, eslocal incx, void *y, eslocal incy);
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
	//cublasGetVector(eslocal n, eslocal elemSize, const void *x, eslocal incx, void *y, eslocal incy)
	//cublasGetVector(rows , sizeof(double), d_y_out, 1, &y_out[0], 1);
	cudaMemcpyAsync(y_out, d_y_out, rows * sizeof(double), cudaMemcpyDeviceToHost, stream);

	//cudaStreamSynchronize(stream);

#endif
}

void SparseMatrix::DenseMatVecCUDA_wo_Copy_sync ( ) {
#ifdef CUDA

	// Retrieve result vector from device
	//cublasGetVector(eslocal n, eslocal elemSize, const void *x, eslocal incx, void *y, eslocal incy)
	//cublasGetVector(rows , sizeof(double), d_y_out, 1, &y_out[0], 1);
	cudaStreamSynchronize(stream);

#endif
}

eslocal SparseMatrix::CopyToCUDA_Dev( ) {
	eslocal error = 0;

#ifdef CUDA

	eslocal mat_size = dense_values.size();// rows * cols;
	eslocal lda = rows;

	if ( d_dense_values == NULL ) {

		cudaError_t status = cudaMalloc((void**)&d_dense_values,   mat_size * sizeof(double));
		if (status != cudaSuccess)   {
      std::cout <<"Error allocating GPU memory \n";
			MPI_Finalize();
			exit(0);
		}


		status = cudaMalloc((void**)&d_x_in,  rows * sizeof(double));
		if (status != cudaSuccess) {
      std::cout <<"Error allocating GPU memory for Matrix \n";
			//MPI_Finalize();
			//exit(0);
			error = -1;
		}


		status = cudaMalloc((void**)&d_y_out, rows * sizeof(double));
		if (status != cudaSuccess) {
      std::cout <<"Error allocating GPU memory for Vector \n";
			//MPI_Finalize();
			//exit(0);
			error = -1;
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
	return error;

}



void SparseMatrix::DenseMatVecCUDA_wo_Copy_start_fl( float * x_in, float * y_out, char T_for_transpose_N_for_not_transpose, eslocal x_in_vector_start_index) {
#ifdef CUDA

	eslocal mat_size = dense_values.size(); //rows * cols;
	eslocal lda = rows;

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
	// cublasSetVector(eslocal n, eslocal elemSize, const void *x, eslocal incx, void *y, eslocal incy);
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
	//cublasGetVector(eslocal n, eslocal elemSize, const void *x, eslocal incx, void *y, eslocal incy)
	//cublasGetVector(rows , sizeof(double), d_y_out, 1, &y_out[0], 1);
	cudaMemcpyAsync(y_out, d_y_out_fl, rows * sizeof(float), cudaMemcpyDeviceToHost, stream);

	//cudaStreamSynchronize(stream);

#endif
}

void SparseMatrix::CopyToCUDA_Dev_fl ( ) {
#ifdef CUDA

	eslocal mat_size = dense_values.size();// rows * cols;
	eslocal lda = rows;

	if ( d_dense_values_fl == NULL ) {

		cudaError_t status = cudaMalloc((void**)&d_dense_values_fl,   mat_size * sizeof(float));
		if (status != cudaSuccess)   {
      std::cout<< "Error allocating GPU memory \n";
			MPI_Finalize();
			exit(0);
		}


		status = cudaMalloc((void**)&d_x_in_fl,  rows * sizeof(float));
		if (status != cudaSuccess) {
      std::cout<<"Error allocating GPU memory  \n";
			MPI_Finalize();
			exit(0);
		}


		status = cudaMalloc((void**)&d_y_out_fl, rows * sizeof(float));
		if (status != cudaSuccess) {
      std::cout<<"Error allocating GPU memory \n";
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
//	LAPACKE_dtrttp( eslocal matrix_layout, char uplo, lapack_eslocal n, const <datatype>* a, lapack_eslocal lda, <datatype>* ap )
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
		matdescra[0] = 'S'; // Triangular Matrix
		matdescra[1] = 'U'; // Triangular indicator: upper
		matdescra[2] = 'N'; // Main diagonal type: non-unit
		matdescra[3] = 'F'; // One based indexing
	}
	//y_out.resize(rows);

	// y := alpha*A*x + beta*y
	//void mkl_dcoomv	(char *transa, eslocal *m, eslocal *k, double *alpha, char *matdescra, double *val,    eslocal *rowind,    eslocal *colind,    eslocal *nnz, double *x,  double *beta, double *y);
	mkl_dcoomv			( &trans,	   &rows,      &cols,      &alpha,        matdescra,      &V_values[0],    &I_row_indices[0],  &J_col_indices[0],  &nnz,         &x_in[0],   &beta,        &y_out[0]);

}



void SparseMatrix::MatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose ) {
	MatVec(x_in, y_out, T_for_transpose_N_for_non_transpose, 0 , 0);
}

void SparseMatrix::MatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose, eslocal x_in_vector_start_index, eslocal y_out_vector_start_index) {
	double beta = 0;
	MatVec(x_in, y_out, T_for_transpose_N_for_non_transpose, x_in_vector_start_index, y_out_vector_start_index, beta);
}

void SparseMatrix::MatVec(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose, eslocal x_in_vector_start_index, eslocal y_out_vector_start_index, double beta) {
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
	// void mkl_dcsrmv (char *transa, eslocal *m, eslocal *k, double *alpha, char *matdescra, double *val,       eslocal *indx,          eslocal *pntrb,     eslocal *pntre,     double *x,						   double *beta,  double *y);
	// mkl_ccsrmv      (&transa,      &m,         &m,         &alpha,              matdescra, values,            columns,                rowIndex,           &(rowIndex[1]),     sol_vec,						   &beta,         rhs_vec);
	mkl_dcsrmv         (&trans,       &rows,      &cols,      &alpha,              matdescra, &CSR_V_values[0],  &CSR_J_col_indices[0],  &CSR_I_row_indices[0],  &CSR_I_row_indices[1],  &x_in[x_in_vector_start_index],   &beta,         &y_out[y_out_vector_start_index]);

}


void SparseMatrix::MatMat(SparseMatrix & A_in, char MatA_T_for_transpose_N_for_non_transpose, SparseMatrix & B_in) {
	// THIS := op(A)*B

	char transa = MatA_T_for_transpose_N_for_non_transpose;

	eslocal job;

	eslocal sort = 3;	   // 3	yes	yes	yes

	eslocal m = A_in.rows; // Number of rows of the matrix A.
	eslocal n = A_in.cols; // Number of columns of the matrix A.
	eslocal k = B_in.cols; // Number of columns of the matrix B.

	SEQ_VECTOR<eslocal>().swap( CSR_I_row_indices );
	SEQ_VECTOR<eslocal>().swap( CSR_J_col_indices );
	SEQ_VECTOR<double>().swap( CSR_V_values );

	if (transa == 'T')
		CSR_I_row_indices.resize( n + 1 );
	else
		CSR_I_row_indices.resize( m + 1 );

	CSR_J_col_indices.resize(1);
	CSR_V_values.resize(1);

	double * a  = &A_in.CSR_V_values[0];
	eslocal    * ia = &A_in.CSR_I_row_indices[0];
	eslocal    * ja = &A_in.CSR_J_col_indices[0];

	double * b  = &B_in.CSR_V_values[0];
	eslocal    * ib = &B_in.CSR_I_row_indices[0];
	eslocal    * jb = &B_in.CSR_J_col_indices[0];

	eslocal nnzmax = 1;

	eslocal ierr;

	//void mkl_dcsrmultcsr (
	//	char *transa, eslocal *job, eslocal *sort,
	//	eslocal *m, eslocal *n, eslocal *k,
	//	double *a, eslocal *ja, eslocal *ia,
	//	double *b, eslocal *jb, eslocal *ib,
	//	double *c, eslocal *jc, eslocal *ic, eslocal *nnzmax,
	//	eslocal *ierr);


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

	eslocal job;

	eslocal sort = 7;	   // 3	no no no

	eslocal m = A_in.rows; // Number of rows of the matrix A.
	eslocal n = A_in.cols; // Number of columns of the matrix A.
	eslocal k = B_in.cols; // Number of columns of the matrix B.

	SEQ_VECTOR<eslocal>().swap( CSR_I_row_indices );
	SEQ_VECTOR<eslocal>().swap( CSR_J_col_indices );
	SEQ_VECTOR<double>().swap( CSR_V_values );

	if (transa == 'T')
		CSR_I_row_indices.resize( n + 1 );
	else
		CSR_I_row_indices.resize( m + 1 );

	CSR_J_col_indices.resize(1);
	CSR_V_values.resize(1);

	double * a  = &A_in.CSR_V_values[0];
	eslocal    * ia = &A_in.CSR_I_row_indices[0];
	eslocal    * ja = &A_in.CSR_J_col_indices[0];

	double * b  = &B_in.CSR_V_values[0];
	eslocal    * ib = &B_in.CSR_I_row_indices[0];
	eslocal    * jb = &B_in.CSR_J_col_indices[0];

	eslocal nnzmax = 1;

	eslocal ierr;

	//void mkl_dcsrmultcsr (
	//	char *transa, eslocal *job, eslocal *sort,
	//	eslocal *m, eslocal *n, eslocal *k,
	//	double *a, eslocal *ja, eslocal *ia,
	//	double *b, eslocal *jb, eslocal *ib,
	//	double *c, eslocal *jc, eslocal *ic, eslocal *nnzmax,
	//	eslocal *ierr);


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

	eslocal job;
	eslocal sort = 3; // 3	yes	yes	yes

	eslocal m = A_in.rows; // Number of rows of the matrix A.
	eslocal n = A_in.cols; // Number of columns of the matrix A.

	eslocal nnzmax;
	eslocal ierr;

	double * a  = &A_in.CSR_V_values[0];
	eslocal    * ia = &A_in.CSR_I_row_indices[0];
	eslocal    * ja = &A_in.CSR_J_col_indices[0];

	double * b  = &B_in.CSR_V_values[0];
	eslocal    * ib = &B_in.CSR_I_row_indices[0];
	eslocal    * jb = &B_in.CSR_J_col_indices[0];

	SEQ_VECTOR<eslocal>().swap( CSR_I_row_indices );
	SEQ_VECTOR<eslocal>().swap( CSR_J_col_indices );
	SEQ_VECTOR<double>().swap( CSR_V_values );

	CSR_I_row_indices.resize( m + 1 );
	CSR_J_col_indices.resize(1);
	CSR_V_values.resize(1);

	//void mkl_dcsradd (
	//	char *transa, eslocal *job, eslocal *sort,
	//	eslocal *m, eslocal *n,
	//	double *a, eslocal *ja, eslocal *ia,
	//	double *beta, double *b, eslocal *jb, eslocal *ib,
	//	double *c, eslocal *jc, eslocal *ic, eslocal *nnzmax,
	//	eslocal *ierr);

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

//TODO: Pozor pro zobecnenou verzi musi byt odkomentovane!
#define GENINVtools

// AM -start: ---------------------------------------------------------------------------
//
void SparseMatrix::spmv_(SparseMatrix & A, double *x, double *Ax){
#ifdef GENINVtools
	eslocal nA = A.cols;
  eslocal offset = A.CSR_I_row_indices[0] ? 1 : 0;
  memset(Ax,0,nA * sizeof(double));
  for (eslocal i = 0; i < nA ; i++) {
    for (eslocal j = A.CSR_I_row_indices[i];j<A.CSR_I_row_indices[i+1];j++) {
      Ax[i] += CSR_V_values[j-offset] * x[CSR_J_col_indices[j-offset]-offset];
      if (j > CSR_I_row_indices[i]) {
        Ax[CSR_J_col_indices[j-offset]-offset] +=
            CSR_V_values[j-offset] * x[CSR_J_col_indices[CSR_I_row_indices[i]-offset]-offset];
      }
    }
  }
#endif
}


void SparseMatrix::getSubDiagBlockmatrix(SparseMatrix & A_in, SparseMatrix & A_out, eslocal i_start, eslocal size_rr){
#ifdef GENINVtools
//
// Function 'getSubDiagBlockmatrix' returns the diagonal block A_in(r,r) from original A_in,
// where r = { i_start , i_start+1 , i_start+2 , ... , istart + size_rr - 1 }
//
//
// rev. 2015-10-10 (A.M.)
//
// step 1: getting nnz of submatrix
  eslocal nnz_new=0;
  eslocal offset = A_in.CSR_I_row_indices[0] ? 1 : 0;
  for (eslocal i = 0;i<size_rr;i++){
    for (eslocal j = A_in.CSR_I_row_indices[i+i_start];j<A_in.CSR_I_row_indices[i+i_start+1];j++){
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
  eslocal ijcnt=0;
  A_out.CSR_I_row_indices[0]=offset;
  for (eslocal i = 0;i<size_rr;i++){
    for (eslocal j = A_in.CSR_I_row_indices[i+i_start];j<A_in.CSR_I_row_indices[i+i_start+1];j++){
      if ((A_in.CSR_J_col_indices[j-offset]-offset)>=i_start &&
                    (A_in.CSR_J_col_indices[j-offset]-offset)<(i_start+size_rr)){
        A_out.CSR_J_col_indices[ijcnt] = (A_in.CSR_J_col_indices[j-offset]) - i_start;
        A_out.CSR_V_values[ijcnt]=A_in.CSR_V_values[j-offset];
        ijcnt++;
      }
    }
    A_out.CSR_I_row_indices[i+1]=offset+ijcnt;
  }
#endif
}


void SparseMatrix::getSubBlockmatrix_rs( SparseMatrix & A_in, SparseMatrix & A_out,
                                          eslocal i_start, eslocal i_size,
                                          eslocal j_start, eslocal j_size){
#ifdef GENINVtools
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
  eslocal nnz_new=0;
  eslocal offset = A_in.CSR_I_row_indices[0] ? 1 : 0;
  for (eslocal i = 0;i<i_size;i++){
    for (eslocal j = A_in.CSR_I_row_indices[i+i_start];j<A_in.CSR_I_row_indices[i+i_start+1];j++){
      if ((A_in.CSR_J_col_indices[j-offset]-offset)>=j_start &&
                      (A_in.CSR_J_col_indices[j-offset]-offset)<(j_start+j_size)){
        nnz_new++;
      }
    }
  }
// step 2: allocation 1d arrays
  A_out.CSR_V_values.resize(nnz_new);
  A_out.CSR_J_col_indices.resize(nnz_new);
  A_out.CSR_I_row_indices.resize(i_size+1);
  A_out.rows=i_size;
  A_out.cols=j_size;
  A_out.nnz=nnz_new;
	A_out.type = 'G';
// step 3: filling 1d arrays
  eslocal ijcnt=0;
  A_out.CSR_I_row_indices[0]=offset;
  for (eslocal i = 0;i<i_size;i++){
    for (eslocal j = A_in.CSR_I_row_indices[i+i_start];j<A_in.CSR_I_row_indices[i+i_start+1];j++){
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
  eslocal offset = CSR_I_row_indices[0] ? 1 : 0;
  printf("%s = [ ...\n",str0);
  for (eslocal i = 0;i<rows;i++){
    for (eslocal j = CSR_I_row_indices[i];j<CSR_I_row_indices[i+1];j++){
      printf("%d %d %3.9e \n",i+1,CSR_J_col_indices[j-offset],CSR_V_values[j-offset]);
    }
  }
  printf("];%s = full(sparse(%s(:,1),%s(:,2),%s(:,3),%d,%d));\n",
                str0,str0,str0,str0,rows,cols);
  if (type=='S'){
    printf("%s=%s+%s'-diag(diag(%s));\n",str0,str0,str0,str0);
  }
}

void SparseMatrix::printMatCSR2(char *str0){
  eslocal offset = CSR_I_row_indices[0] ? 1 : 0;

  FILE *fid = fopen(str0,"w");
  int isGeneral=0;
  if (type=='G') isGeneral=1;
  fprintf(fid,"%d %d %d\n",rows,cols,isGeneral);

  for (eslocal i = 0;i<rows;i++){
    for (eslocal j = CSR_I_row_indices[i];j<CSR_I_row_indices[i+1];j++){
      fprintf(fid,"%d %d %3.9e \n",i+1,CSR_J_col_indices[j-offset],CSR_V_values[j-offset]);
    }
  }
#endif
}


double SparseMatrix::getNorm_K_R(SparseMatrix & K, SparseMatrix &R_in_dense_format){
#ifdef GENINVtools
  double * AR =  new double [K.rows];
  double norm_AR_row,norm_AR = 0.0;
  for (eslocal i = 0;i<R_in_dense_format.cols;i++){
    memset(AR,0,R_in_dense_format.rows * sizeof(double));
  	K.spmv_( K,&(R_in_dense_format.dense_values[i*R_in_dense_format.rows]),AR);
    norm_AR_row=0.0;
    for (eslocal j = 0; j < R_in_dense_format.rows;j++){
      norm_AR_row+=AR[j]*AR[j];
    }
    norm_AR+=norm_AR_row;
  }
  delete [] AR;
  norm_AR=sqrt(norm_AR);
  return norm_AR;
#endif
}

//

void SparseMatrix::GramSchmidtOrtho(){
#ifdef GENINVtools
  double *w = new double [rows];
  double *R = new double [cols*cols];
  memset(R,0,(cols*cols) * sizeof(double));

  for (eslocal j = 0;j<cols;j++){
    memcpy( w, &(dense_values[j*rows]) , sizeof( double ) * rows);
    for (eslocal i = 0;i<j;i++){
      R[j*cols+i] = dot_e(w, &(dense_values[i*rows]),rows);
      for (eslocal k=0;k<rows;k++){
        w[k]-=dense_values[i*rows+k]*R[j*cols+i];
      }
    }
    R[j*cols+j] = sqrt(dot_e(w,w,rows));
    for (eslocal k=0;k<rows;k++){
      dense_values[j*rows+k] = w[k]/R[j*cols+j];
    }
  }
  delete [] w;
  delete [] R;
#endif
}


bool myfn(double i, double j) { return fabs(i)<=fabs(j); }

void SparseMatrix::getNullPivots(SEQ_VECTOR <eslocal> & null_pivots){
#ifdef GENINVtools
	SEQ_VECTOR <double> N(dense_values);
  eslocal nEl = rows*cols;
  std::vector <double>::iterator  it;
  eslocal I,J,K,colInd,rowInd;
  double *tmpV = new double[rows];
  double pivot;
  eslocal tmp_int;
  eslocal *_nul_piv = new eslocal[rows];
  for (eslocal i = 0;i<rows;i++) _nul_piv[i]=i;

//TODO Ask about to the efficiency of next 2 lines.
  auto ij= [&]( eslocal ii, eslocal jj ) -> eslocal
   { return ii + rows*jj; };
  for (eslocal j=0;j<cols;j++){
    it = std::max_element(N.begin(),N.end()-j*rows,myfn);
    I = it - N.begin();
    colInd = I/rows;
    rowInd = I-colInd*rows;
    for (eslocal k=0;k<cols-j;k++){
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
    for (eslocal J=0;J<cols-j-1;J++){
      for (eslocal I=0;I<rows-j;I++){
        N[ij(I,J)] -= N[ij(I,cols-1-j)]*N[ij(rows-1-j,J)]/pivot;
      }
    }
  }
  for (eslocal i = 0;i<cols;i++){
    null_pivots.push_back(_nul_piv[rows-1-i]+1);
  }
  sort(null_pivots.begin(),null_pivots.end());
//
  delete [] _nul_piv;
  delete [] tmpV;
//
#endif
}
//
double SparseMatrix::MatCondNumb( SparseMatrix & A_in, char *str0, eslocal plot_n_first_n_last_eigenvalues,
                    double *maxEig,  int nMax_input){
#ifdef GENINVtools
  bool plot_a_and_b_defines_tridiag=false;
  eslocal nA = A_in.rows;
  eslocal nMax = 200; // size of tridiagonal matrix
  eslocal offset = CSR_I_row_indices[0] ? 1 : 0;
  //eslocal nEigToplot = 10;
  double *s = new double[nA];
  double *s_bef = new double[nA];
  double *As = new double[nA];
  double *r = new double[nA];
  double tmp_a,alpha, beta = 1.0, beta_bef;
  double estim_cond;

  if (nMax_input>0 && nMax_input<nMax){
    nMax = nMax_input;
  }

//
  nMax = nMax > nA ? nA : nMax;
  double *alphaVec = new double[nMax];
  double *betaVec  = new double[nMax];
  eslocal cnt = 0;
//
  memset(s,0,nA * sizeof(double));
  for (eslocal i = 0 ; i < nA; i++){ r[i] = i ; }
  tmp_a = sqrt(dot_e(r,r,nA));
  for (eslocal i = 0 ; i < nA; i++){ r[i] /=  tmp_a; }
//
  for (eslocal i = 0; i < nMax ; i++){
    memcpy( s_bef, s , sizeof( double ) * nA);
    beta_bef=beta;
    memcpy( s, r , sizeof( double ) * nA);
    for (eslocal j =  0;j < nA; j++){
      s[j]/=beta;
    }
    spmv_(A_in,s,As);
    alpha = dot_e(s,As,nA);
    for (eslocal j =  0;j < nA; j++){
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
    for (eslocal i = 0 ; i < cnt; i++){
      printf("%3.8e %3.8e\n",alphaVec[i],betaVec[i]);
    }
  }
  char JOBZ = 'N';
  double *Z = new double[cnt];
  eslocal info;
  eslocal ldz = cnt;
  info = LAPACKE_dstev(LAPACK_ROW_MAJOR, JOBZ, cnt, alphaVec, betaVec, Z, ldz);
  estim_cond=fabs(alphaVec[cnt-1]/alphaVec[0]);
  if (plot_n_first_n_last_eigenvalues>0){
    printf("cond(%s) = %3.15e\tit: %d\n",str0,estim_cond,cnt);
  }
  *maxEig = alphaVec[cnt-1];

  if (plot_n_first_n_last_eigenvalues>0){
    std::cout<<"eigenvals of "<<str0 <<" d{1:" << plot_n_first_n_last_eigenvalues << "} and d{" <<
         cnt-plot_n_first_n_last_eigenvalues+2 << ":\t"<< cnt<< "}\n";


    for (eslocal i = 0 ; i < cnt; i++){
      if (i < plot_n_first_n_last_eigenvalues || i > cnt-plot_n_first_n_last_eigenvalues){
        std::cout<< i+1 <<":"<< alphaVec[i] << "\n";
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
#endif
  return 0;
}

double SparseMatrix::dot_e(double *x, double *y, eslocal n){
  double dot_xy = 0.0;
  for (eslocal i = 0; i< n; i++){
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


	eslocal job;
	eslocal sort = 3; // 3	yes	yes	yes

	eslocal m = rows; // Number of rows of the matrix A.
	eslocal n = cols; // Number of columns of the matrix A.

	eslocal nnzmax;
	eslocal ierr;

	double 	   * a  = &CSR_V_values[0];
	eslocal    * ia = &CSR_I_row_indices[0];
	eslocal    * ja = &CSR_J_col_indices[0];

	double     * b  = &B_in.CSR_V_values[0];
	eslocal    * ib = &B_in.CSR_I_row_indices[0];
	eslocal    * jb = &B_in.CSR_J_col_indices[0];

	SEQ_VECTOR<eslocal>		t_CSR_I_row_indices;	t_CSR_I_row_indices.resize( m + 1 );
	SEQ_VECTOR<eslocal>		t_CSR_J_col_indices;	t_CSR_J_col_indices.resize(1);
	SEQ_VECTOR<double>		t_CSR_V_values;			t_CSR_V_values.resize(1);

	//void mkl_dcsradd (
	//	char *transa, eslocal *job, eslocal *sort,
	//	eslocal *m, eslocal *n,
	//	double *a, eslocal *ja, eslocal *ia,
	//	double *beta, double *b, eslocal *jb, eslocal *ib,
	//	double *c, eslocal *jc, eslocal *ic, eslocal *nnzmax,
	//	eslocal *ierr);

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
	for (eslocal i = 0; i < CSR_V_values.size(); i++) {
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
//	eslocal job		=  1;
//	eslocal sort	=  3; // 3	yes	yes	yes
//	eslocal nnzmax	=  1;
//
//	eslocal m = this->cols; // Number of rows of the matrix A.
//	eslocal n = this->rows; // Number of columns of the matrix A.
//
//	eslocal ierr;
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
//	//	char *transa, eslocal *job, eslocal *sort,
//	//	eslocal *m, eslocal *n,
//	//	double *a, eslocal *ja, eslocal *ia,
//	//	double *beta, double *b, eslocal *jb, eslocal *ib,
//	//	double *c, eslocal *jc, eslocal *ic, eslocal *nnzmax,
//	//	eslocal *ierr);
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

	eslocal job [] = { 0, 1, 1, 0, 0, 1 };
	eslocal info;
	eslocal m;

	A_out.cols = rows;
	A_out.rows = cols;

	eslocal row_size_backup = CSR_I_row_indices.size();
	if (cols > rows) {
		CSR_I_row_indices.resize( cols+1,  CSR_I_row_indices[CSR_I_row_indices.size()-1] );
		m = cols;
	} else {
		m = rows;
	}

	SEQ_VECTOR<eslocal>().swap( A_out.CSR_I_row_indices );
	SEQ_VECTOR<eslocal>().swap( A_out.CSR_J_col_indices );
	SEQ_VECTOR<double>().swap( A_out.CSR_V_values );

	A_out.CSR_I_row_indices.resize(m + 1);
	A_out.CSR_J_col_indices.resize(nnz);
	A_out.CSR_V_values.		resize(nnz);

	//void mkl_dcsrcsc(eslocal *job, eslocal *m, double *acsr, eslocal *ja, eslocal *ia, double *acsc, eslocal *ja1, eslocal *ia1, eslocal *info);
	//mkl_dcsrcsc( &job[0], &m, &CSR_V_values[0], &CSR_J_col_indices[0], &CSR_I_row_indices[0], &V_values[0], &J_col_indices[0], &I_row_indices[0], &info);

	mkl_dcsrcsc( &job[0], &m,
		&CSR_V_values[0], &CSR_J_col_indices[0], &CSR_I_row_indices[0],
		&A_out.CSR_V_values[0], &A_out.CSR_J_col_indices[0], &A_out.CSR_I_row_indices[0],
		&info);

	if (cols > rows) {
		CSR_I_row_indices.resize(row_size_backup);
		SEQ_VECTOR<eslocal> tmp;
		tmp = CSR_I_row_indices;
		CSR_I_row_indices.swap(tmp);
	} else {
		A_out.CSR_I_row_indices.resize(A_out.rows + 1);
		SEQ_VECTOR<eslocal> tmp;
		tmp = A_out.CSR_I_row_indices;
		A_out.CSR_I_row_indices.swap(tmp);
	}

	A_out.nnz  = nnz;
	A_out.type = type;

}

void SparseMatrix::MatTranspose() {

	eslocal job [] = { 0, 1, 1, 0, 0, 1 };
	eslocal info;
	eslocal m;
	eslocal row_size_backup = CSR_I_row_indices.size();

	SEQ_VECTOR <eslocal> tCSR_I_row_indices;
	SEQ_VECTOR <eslocal> tCSR_J_col_indices;
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

	//void mkl_dcsrcsc(eslocal *job, eslocal *m, double *acsr, eslocal *ja, eslocal *ia, double *acsc, eslocal *ja1, eslocal *ia1, eslocal *info);
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

	SEQ_VECTOR<eslocal>().swap( CSR_I_row_indices );
	CSR_I_row_indices = tCSR_I_row_indices; //.swap(tCSR_I_row_indices);
	CSR_J_col_indices.swap(tCSR_J_col_indices);
	CSR_V_values     .swap(tCSR_V_values);

}

void SparseMatrix::MatTransposeCOO() {

	I_row_indices.swap(J_col_indices);
	eslocal tmp = rows;
	rows = cols;
	cols = tmp;

	//sortInCOO();  // TODO: musi se zpatky povolit hned jak se opravit funkce sort v COO

}


void SparseMatrix::RemoveLower() {

	SEQ_VECTOR <eslocal> t_CSR_I_row_indices;
	SEQ_VECTOR <eslocal> t_CSR_J_col_indices;
	SEQ_VECTOR <double> t_CSR_V_values;
	eslocal l_nnz = 0;

	for (eslocal row = 0; row < CSR_I_row_indices.size() - 1; row++) {
		t_CSR_I_row_indices.push_back(l_nnz+1);
		eslocal cols_in_row = CSR_I_row_indices[row+1] - CSR_I_row_indices[row];

		for (eslocal col = 0; col <cols_in_row; col++) {
			eslocal i = CSR_I_row_indices[row] - 1 + col;
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
	eslocal count = 0;

	for (eslocal i = 0; i < CSR_I_row_indices.size() - 1; i++) {
		double val = CSR_V_values[ CSR_I_row_indices[i] - 1 ];
		sum = sum + val;
		count++;
	}

	return sum/count;
}

double SparseMatrix::GetMaxOfDiagonalOfSymmetricMatrix() {
	double vmax = 0;

	for (eslocal i = 0; i < CSR_I_row_indices.size() - 1; i++) {

		if ( vmax < CSR_V_values[ CSR_I_row_indices[i] - 1 ] )
			vmax = CSR_V_values[ CSR_I_row_indices[i] - 1 ];

	}

	return vmax;
}


void SparseMatrix::SetDiagonalOfSymmetricMatrix( double val ) {
	for (eslocal i = 0; i < CSR_I_row_indices.size() - 1; i++) {
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

		eslocal last_row = CSR_I_row_indices[CSR_I_row_indices.size()-1];

		for (eslocal i = 1; i < A.CSR_I_row_indices.size(); i++)
			CSR_I_row_indices.push_back(last_row + A.CSR_I_row_indices[i]-1);

		rows = rows + A.rows;
		nnz  = nnz + A.nnz;
		cols = cols;
		type = 'G';
	}
}


void SparseMatrix::CreateMatFromRowsFromMatrix(SparseMatrix & A_in, SEQ_VECTOR <eslocal> & rows_to_add) {

	eslocal old_index  = 0;
	eslocal next_index = 0;
	eslocal row_fill   = 1;

	rows = A_in.rows;
	cols = A_in.cols;
	type = A_in.type;

	CSR_I_row_indices.resize( rows + 1 );

	for (eslocal i = 0; i < rows_to_add.size(); i++) {

		old_index  = next_index;
		next_index = rows_to_add[i];

		fill(CSR_I_row_indices.begin() + old_index, CSR_I_row_indices.begin() + next_index, row_fill);

                eslocal A_in_start_index;
                if (rows_to_add[i] > 0)
		   A_in_start_index = A_in.CSR_I_row_indices[rows_to_add[i] - 1 ] - 1 ;
		else
                   A_in_start_index = 0;

                eslocal A_in_end_index   = A_in.CSR_I_row_indices[rows_to_add[i] + 1 - 1] - 1 ;

		CSR_J_col_indices.insert(CSR_J_col_indices.end(), A_in.CSR_J_col_indices.begin() + A_in_start_index, A_in.CSR_J_col_indices.begin() + A_in_end_index );
		CSR_V_values.     insert(CSR_V_values.end(),      A_in.CSR_V_values.     begin() + A_in_start_index, A_in.CSR_V_values.     begin() + A_in_end_index );
		row_fill = 1 + CSR_J_col_indices.size();

	}


	fill(CSR_I_row_indices.begin() + next_index, CSR_I_row_indices.begin() + rows + 1, row_fill);

	nnz = CSR_V_values.size();

}

void SparseMatrix::CreateMatFromRowsFromMatrix_NewSize(SparseMatrix & A_in, SEQ_VECTOR <eslocal> & rows_to_add) {

	int old_index  = 0;
	int next_index = 0;
	int row_fill   = 1;

	rows = rows_to_add.size();
	cols = A_in.cols;
	type = A_in.type;

	//CSR_I_row_indices.resize( rows + 1 );

	for (int i = 0; i < rows_to_add.size(); i++) {

		old_index  = next_index;
		next_index = rows_to_add[i];

		//fill(CSR_I_row_indices.begin() + old_index, CSR_I_row_indices.begin() + next_index, row_fill);
                CSR_I_row_indices.push_back(row_fill);

		int A_in_start_index = A_in.CSR_I_row_indices[rows_to_add[i] - 1 ] - 1 ;
		int A_in_end_index   = A_in.CSR_I_row_indices[rows_to_add[i] + 1 - 1] - 1 ;

		CSR_J_col_indices.insert(CSR_J_col_indices.end(), A_in.CSR_J_col_indices.begin() + A_in_start_index, A_in.CSR_J_col_indices.begin() + A_in_end_index );
		CSR_V_values.     insert(CSR_V_values.end(),      A_in.CSR_V_values.     begin() + A_in_start_index, A_in.CSR_V_values.     begin() + A_in_end_index );
		row_fill = 1 + CSR_J_col_indices.size();

	}


	//fill(CSR_I_row_indices.begin() + next_index, CSR_I_row_indices.begin() + rows + 1, row_fill);
        CSR_I_row_indices.push_back(row_fill);

	nnz = CSR_V_values.size();

}

eslocal SparseMatrix::MatCompare(SparseMatrix & A) {
	eslocal res = 0;

	if (this->cols == A.cols && this->rows==A.rows && this->nnz == A.nnz && this->type == A.type ) {

		eslocal tmp1 = 0;
		eslocal tmp2 = 0;

		for (eslocal i = 0; i < CSR_I_row_indices.size(); i++)
			if (CSR_I_row_indices[i] != A.CSR_I_row_indices[i])
				tmp1=1;

		for (eslocal i = 0; i < CSR_J_col_indices.size(); i++)
			if (CSR_J_col_indices[i] != A.CSR_J_col_indices[i])
				tmp1=1;

		for (eslocal i = 0; i < CSR_V_values.size(); i++)
			if (CSR_V_values[i] != A.CSR_V_values[i])
				tmp2=1;

		res = 1000 * tmp1 + tmp2;

	} else {
		res = -1;
	}

	return res;
}

eslocal SparseMatrix::MatCompareCOO(SparseMatrix & A) {
	eslocal res = 0;

	if (this->cols == A.cols && this->rows==A.rows && this->nnz == A.nnz && this->type == A.type ) {

		eslocal tmp1 = 0;
		eslocal tmp2 = 0;
		eslocal tmp3 = 0;

		for (eslocal i = 0; i < I_row_indices.size(); i++)
			if (I_row_indices[i] != A.I_row_indices[i])
				tmp1=1;

		for (eslocal i = 0; i < J_col_indices.size(); i++)
			if (J_col_indices[i] != A.J_col_indices[i])
				tmp2=1;

		for (eslocal i = 0; i < V_values.size(); i++)
			if (V_values[i] != A.V_values[i])
				tmp3=1;

		res = 100 * tmp1 + 10 * tmp2 + tmp3;

	} else {
		res = -1;
	}

	return res;
}


void SparseMatrix::CreateEye(eslocal size) {

	for (eslocal i = 0; i< size; i++) {
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


void SparseMatrix::CreateEye(eslocal size, double value, eslocal offset_row, eslocal offset_col) {

	for (eslocal i = 0; i< size; i++) {
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



void SparseMatrix::TestEye(eslocal size) {

	for (eslocal i = 0; i< size; i++) {
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

void SparseMatrix::TestMatRow(eslocal size, eslocal row_index) {

	for (eslocal i = 0; i< size; i++) {
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

	//SEQ_VECTOR < SEQ_VECTOR < eslocal    > > GGt_J (A_in.rows, SEQ_VECTOR < eslocal    > () );
	//SEQ_VECTOR < SEQ_VECTOR < double > > GGt_V (A_in.rows, SEQ_VECTOR < double > () );

	SEQ_VECTOR<eslocal>()   .swap( CSR_I_row_indices );
	SEQ_VECTOR<eslocal>()   .swap( CSR_J_col_indices );
	SEQ_VECTOR<double>().swap( CSR_V_values );

	eslocal glob_row_index = 0 + 1;
	CSR_I_row_indices.push_back(glob_row_index);

	for (eslocal i = 0; i < A_in.CSR_I_row_indices.size() - 1; i++ ) {

		eslocal A_row_start = A_in.CSR_I_row_indices[i  ] - 1;
		eslocal A_row_end   = A_in.CSR_I_row_indices[i+1] - 1;

		if (A_row_start != A_row_end ) { // this row in B is NOT empty

			for (eslocal ii = 0; ii < B_in.CSR_I_row_indices.size() - 1; ii++ ) {

				eslocal B_row_start = B_in.CSR_I_row_indices[ii  ] - 1;
				eslocal B_row_end   = B_in.CSR_I_row_indices[ii+1] - 1;

				if (B_row_start != B_row_end) { // this row in B is NOT empty
					eslocal A_ind = A_row_start;
					eslocal B_ind = B_row_start;
					double C_v = 0;
					do {

						eslocal A_j = A_in.CSR_J_col_indices[A_ind];
						eslocal B_j = B_in.CSR_J_col_indices[B_ind];
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
		eslocal a = 10;
	}

	rows = A_in.rows;
	cols = A_in.rows;
	nnz  = CSR_V_values.size();
	type = 'G';

}



//void SparseMatrix::get_kernel_from_K() {
//	get_kernel_from_K(K, Kplus_R);
//}

void SparseMatrix::get_kernel_from_K(SparseMatrix &K, SparseMatrix &regMat,SparseMatrix &Kplus_R,double *norm_KR_d_pow_2,eslocal *defect_d,eslocal d_sub){
//
// Routine calculates kernel Kplus_R of K satisfied euqality K * Kplus_R = O,
// where O is zero matrix, and it makes the matrix K non-singular (K_reg)
// utilizing spectral conditions of Schur complement. Then ||K-K*inv(K_reg)*K||=0.0
//
//
// rev. 2016-02-03 (A.M.)
//==============================================================================
//
#ifndef VERBOSE_LEVEL
#define VERBOSE_LEVEL 0
#endif
//
//    1) diagonalScaling
//  reducing of big jump coefficient effect (TODO include diagonal scaling into whole ESPRESO)
//BOOL DIAGONALSCALING                                  = true;
  bool diagonalScaling                                  = true;

//    2) permutVectorActive
//  random selection of singular DOFs
// 0 - no permut., 1 - std::vector shuffle, 2 - generating own random sequence -
//ESLOCAL PERMUTVECTORACTIVE                            = 1;
  eslocal permutVectorActive                            = 1;

//    3) use_null_pivots_or_s_set
  // NtN_Mat from null pivots or fixing DOFs
//BOOL USE_NULL_PIVOTS_OR_S_SET                         = TRUE;
  bool use_null_pivots_or_s_set                         = true;

//    4) diagonalRegularization
//  regularization only on diagonal elements (big advantage: patern of K and K_regular is the same !!!)
//  size of set 's' = defect(K)
//  It's is active, only if and only if 'use_null_pivots_or_s_set = true'
//BOOL DIAGONALREGULARIZATION                           = TRUE;
  bool diagonalRegularization                           = true;

//    5) get_n_first_and_n_last_eigenvals_from_dense_K
// get and print 2*n K eigenvalues (K is temporarily converted to dense);
//ESLOCAL GET_N_FIRST_AND_N_LAST_EIGENVALS_FROM_DENSE_K = 0;
  eslocal get_n_first_and_n_last_eigenvals_from_dense_K = 0;

//    6) get_n_first_and_n_last_eigenvals_from_dense_S
// get and print 2*n S eigenvalues
//ESLOCAL GET_N_FIRST_AND_N_LAST_EIGENVALS_FROM_DENSE_S = 0;
  eslocal get_n_first_and_n_last_eigenvals_from_dense_S = 10;

//    7) plot_n_first_n_last_eigenvalues
// get of K eigenvalues (K is temporarily converted to dense matrix);
//ESLOCAL PLOT_N_FIRST_N_LAST_EIGENVALUES               = 0;
  eslocal plot_n_first_n_last_eigenvalues               = 0;

//    8) fixing_nodes_or_dof
// non-singular part determined by fixing nodes (FN),
// min(fixing_nodes_or_dof)>=3; if variable is nonzero,
// parameter sc_size is set to fixing_nodes_or_dof*dofPerNode
//ESLOCAL FIXING_NODES_OR_DOF                           = 0;
  eslocal fixing_nodes_or_dof = 0;
//ESLOCAL DOFPERNODE                                    = 3;
  eslocal dofPerNode                                    = 3;
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    1) cond_numb_for_singular_matrix
//  If cond(K) > cond_numb_for_singular_matrix, K is considered as singular matrix.
//DOUBLE COND_NUMB_FOR_SINGULAR_MATRIX                  = 1E13;
  double cond_numb_for_singular_matrix                  = 1e13;

//    2) check_nonsing
// if check_nonsing>0, checking of K_rr non-singularity is activated and it is repeated
// (check_nonsing) times.
//ESLOCAL CHECK_NONSING                                 = 0;
  eslocal check_nonsing                                 = 0;

//    3) max_size_of_dense_matrix_to_get_eigs
// if size of K is less then CHECK_N..., K is converted to dense format to get eigenvalues.
//ESLOCAL MAX_SIZE_OF_DENSE_MATRIX_TO_GET_EIGS          = 2500;
  eslocal max_size_of_dense_matrix_to_get_eigs          = 2500;

//    4) sc_size
// specification of size of Schur complement used for detection of zero eigenvalues.
//eslocal  sc_size >= expected defect 'd' (e.g. in elasticity d=6).
//ESLOCAL SC_SIZE                                       = 50;
  eslocal sc_size                                       = 50;

//    5) twenty
// testing last twenty eigenvalues of S to distinguish, if d-last ones are zero or not.
//ESLOCAL TWENTY                                        = 20;
  eslocal twenty                                        = 20;
  // twenty eigenvalues are ascendly ordered in d = d[0],d[1], ..., d[n-2],d[n-1]

//    6) jump_in_eigenvalues_alerting_singularity
// if d[i]/d[i+1]< jump_in_eigenvalues_alerting_singularity, d[i] is last nonzero eigenvalue
//DOUBLE JUMP_IN_EIGENVALUES_ALERTING_SINGULARITY       = 1.0E-5;
  double jump_in_eigenvalues_alerting_singularity       = 1.0e-5;

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

  int MPI_rank=config::MPIrank;
  double begin_time = omp_get_wtime();

// DEFAULT SET-UP
#if VERBOSE_LEVEL < 4
    diagonalScaling                               = DIAGONALSCALING;
    permutVectorActive                            = PERMUTVECTORACTIVE;
    use_null_pivots_or_s_set                      = USE_NULL_PIVOTS_OR_S_SET;
    diagonalRegularization                        = DIAGONALREGULARIZATION;
    get_n_first_and_n_last_eigenvals_from_dense_K = GET_N_FIRST_AND_N_LAST_EIGENVALS_FROM_DENSE_K;
    get_n_first_and_n_last_eigenvals_from_dense_S = GET_N_FIRST_AND_N_LAST_EIGENVALS_FROM_DENSE_S;
    plot_n_first_n_last_eigenvalues               = PLOT_N_FIRST_N_LAST_EIGENVALUES;
    fixing_nodes_or_dof                           = FIXING_NODES_OR_DOF;
    dofPerNode                                    = DOFPERNODE;
    cond_numb_for_singular_matrix                 = COND_NUMB_FOR_SINGULAR_MATRIX;
    check_nonsing                                 = CHECK_NONSING;
    max_size_of_dense_matrix_to_get_eigs          = MAX_SIZE_OF_DENSE_MATRIX_TO_GET_EIGS;
    sc_size                                       = SC_SIZE;
    twenty                                        = TWENTY;
    jump_in_eigenvalues_alerting_singularity      = JUMP_IN_EIGENVALUES_ALERTING_SINGULARITY;
#endif

// STATISTICS MADE AND PRINTED TO FILE
//        'kernel_detct_cX_dY.txt  (X - clust. number, Y - subdomain. number)
//  -BRIEF
#if VERBOSE_LEVEL == 2
    // print ||K*R|| to file
#endif
//  - DETAILED
#if VERBOSE_LEVEL == 3
    get_n_first_and_n_last_eigenvals_from_dense_S = 10;
    check_nonsing                                 = 1;
    // max(eig(K))
#endif
//  - OWN
#if VERBOSE_LEVEL == 4
    std::cout << "debug set-up \n";
#endif


  if (!use_null_pivots_or_s_set) diagonalRegularization=false;

#if VERBOSE_LEVEL>0
    char filenameI[128];
    if (d_sub==-1){
      sprintf(filenameI, "kernel_detct_c_%d_GGt.txt", MPI_rank);
    }
    else{
      sprintf(filenameI, "kernel_detct_c_%d_s%d.txt", MPI_rank,d_sub);
    }
    std::ofstream os (filenameI);
    os.precision(15);


    os << "Verbose Level:       " ;
    if (VERBOSE_LEVEL==1){
      os << "1/4   (times)\n";
    }
    else if (VERBOSE_LEVEL==2){
      os << "2/4   (brief)\n";
    }
    else if (VERBOSE_LEVEL==3){
      os << "3/4   (detailed)\n";
    }
    else if (VERBOSE_LEVEL==4){
      os << "4/4   (own)\n";
    }
    os << "diagonalScaling:     " << int(diagonalScaling)<< "\n";
    os << "permutVectorActive:  " << diagonalScaling<< "\n";
    if (use_null_pivots_or_s_set){
      os << "K regularized by null pivots";
      if (diagonalRegularization){
        os << ", pattern of K is not changed \n\tcontributions only on diagonal elements";
      }
      os << "\n";
    }
    else{
      os << "K regularized by set of fixing nodes\n";
    }

    if (fixing_nodes_or_dof==0){
      os << "non-singular part chosen by set of DOF\n";
    }
    else{
      os << "non-singular part chosen by set of fixing nodes (FN), number of FN: "<<fixing_nodes_or_dof << "\n";
      os << "DOF per node: " <<  dofPerNode << "\n";
    }

    os << std::scientific;
    os.precision(4);
    os << "cond_numb_for_singular_matrix:             " << cond_numb_for_singular_matrix << "\n";
    os << "jump_in_eigenvalues_alerting_singularity   " << jump_in_eigenvalues_alerting_singularity << "\n";
    os << std::fixed;
    os.precision(15);
    os << "check_nonsing:                             " << check_nonsing << "\n";
    os << "max_size_of_dense_matrix_to_get_eigs       " << max_size_of_dense_matrix_to_get_eigs << "\n";

#endif

  //TODO if K.rows<=sc_size, use directly input K instead of S
  //
  int n_nodsSub = 0;
  double rho = K.GetMaxOfDiagonalOfSymmetricMatrix();
  if (fixing_nodes_or_dof>0){
    sc_size = fixing_nodes_or_dof*dofPerNode;
    n_nodsSub = round(K.rows/dofPerNode);
  }
  //
  //##########################################################################################
  //
  SparseMatrix S;
  SparseMatrix K_rr;
  SparseMatrix K_rs;
  eslocal i_start = 0;

  if (K.rows<sc_size){
	  sc_size = K.rows;
    fixing_nodes_or_dof=0;
  }


#if VERBOSE_LEVEL>1
    os << "dim(SchurComplement):                      " <<  sc_size << "\n";
    os << "===================================================================================\n";
#endif

  eslocal nonsing_size = K.rows - sc_size - i_start;
  eslocal j_start = nonsing_size;
  SEQ_VECTOR <eslocal > permVec;

  permVec.resize(K.rows);
  SEQ_VECTOR <SEQ_VECTOR<eslocal >> vec_I1_i2(K.rows, SEQ_VECTOR<eslocal >(2, 1));
  eslocal offset = K.CSR_I_row_indices[0] ? 1 : 0;
  //

  double cond_of_regular_part=1e307;
  eslocal *I_row_indices_p = new eslocal [K.nnz] ;
  eslocal *J_col_indices_p = new eslocal [K.nnz] ;
  SEQ_VECTOR <eslocal > tmp_vec_s;
  tmp_vec_s.resize(sc_size);
  eslocal v1, n_mv, cnt_permut_vec;
  SEQ_VECTOR <eslocal >::iterator it;
  SEQ_VECTOR <eslocal > fix_dofs;
  fix_dofs.resize(sc_size);
  SparseMatrix K_modif;

  double di=1,dj=1;
  eslocal cnt_iter_check_nonsing=0;


  K_modif = K;

  double elapsed_secs[15];

  double time1 = omp_get_wtime();
//0 - allocation of vectors, copying of matrix
  elapsed_secs[0] = (time1 - begin_time) ;



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

#if VERBOSE_LEVEL>0
//1 - diagonal scaling
  time1 = omp_get_wtime();
  elapsed_secs[1] = (time1 - begin_time) ;
#endif
               //                                               |

  //#################################################################################
  if (get_n_first_and_n_last_eigenvals_from_dense_K &&
      K_modif.cols<max_size_of_dense_matrix_to_get_eigs && cnt_iter_check_nonsing==0) {
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


#if VERBOSE_LEVEL>1
      os<<"eigenvals of K d{1:" << get_n_first_and_n_last_eigenvals_from_dense_K << "} and d{" <<
           K_modif.rows-get_n_first_and_n_last_eigenvals_from_dense_K+2 << ":"<< K_modif.rows<< "}\n";


      for (eslocal i = 0 ; i < K_modif.rows; i++){
        if (i < get_n_first_and_n_last_eigenvals_from_dense_K ||
              i > K_modif.rows-get_n_first_and_n_last_eigenvals_from_dense_K){
          os<< i+1 <<":"<< WK_modif[i] << "\n";
        }
      }
#endif
    if (info){
      std::cout <<"info = " << info << " something wrong with Schur complement in SparseSolver::generalIinverse\n";
    }
    delete [] WK_modif;
    delete [] ZK_modif;
  }
  //#################################################################################


//    std::cout<<"eigenvals of K d{1:" << get_n_first_and_n_last_eigenvals_from_dense_K << " and d{" <<
//         K_modif.rows-get_n_first_and_n_last_eigenvals_from_dense_K+2 << ":"<< K_modif.rows<< "}\n";
//
//        std::cout<< i+1 <<":"<< WK_modif[i] << "\n";

#if VERBOSE_LEVEL>0
//2 - before singular test
  time1 = omp_get_wtime();
  elapsed_secs[2] = (time1 - begin_time) ;
#endif
               //                                               |

  while ( cond_of_regular_part > cond_numb_for_singular_matrix && cnt_iter_check_nonsing<(check_nonsing+1)) {
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
      if (fixing_nodes_or_dof==0 || (permutVectorActive==0)){
        for (eslocal i=0; i<K.rows; ++i) { permVec[i]=i;} // 0 1 2 K.rows-1
      }
      else
      {
        for (eslocal i=0; i<n_nodsSub; ++i) { permVec[i]=i;} // 0 1 2 n_nodsSub-1
      }
    }
//
    if (permutVectorActive==1){
      srand(time(NULL));
//      srand(0); // random will be constant until next compiling

      if (fixing_nodes_or_dof==0){
        random_shuffle ( permVec.begin(), permVec.end() );
      }
      else
      {
        std::srand(std::time(0));
        std::random_shuffle ( permVec.begin(), permVec.begin()+n_nodsSub);
        for (eslocal i=n_nodsSub;i>0;i--){
          for (eslocal j=0;j<dofPerNode;j++){
            permVec[dofPerNode*i-1-j] = dofPerNode*permVec[i-1]+j;
          }
        }
      }

      sort(permVec.begin(),permVec.begin()+nonsing_size);
      sort(permVec.begin()+nonsing_size,permVec.end());
    }
    else if (permutVectorActive==2){
      // random permutation
      n_mv = 0;                     // n_mv = size(unique(tmp_vec_s)) has to be equal to sc_size
      cnt_permut_vec=0;
      srand(time(NULL));
      // loop controls, if series 'tmp_vec_s' with unique integers has suffisciant dimension.
      // If not, missing numbers are added and checked again.
      do {
        for (eslocal i = 0;i<(sc_size-n_mv);i++){
          v1 = rand() % K_modif.rows;
          tmp_vec_s[n_mv+i]=v1;
        }
        it=tmp_vec_s.begin();
        std::sort(tmp_vec_s.begin(), tmp_vec_s.end());
        it = std::unique(tmp_vec_s.begin(), tmp_vec_s.end());
        n_mv = distance(tmp_vec_s.begin(),it);
        cnt_permut_vec++;
     } while (n_mv != sc_size && cnt_permut_vec < 100);
      //
      eslocal ik=0,cnt_i=0;
      for (eslocal i = 0;i<permVec.size();i++){
        if (i==tmp_vec_s[ik]){
          permVec[ik+nonsing_size]=tmp_vec_s[ik];
          ik++;
        }
        else{
          permVec[cnt_i]=i;
          cnt_i++;
        }
      }
#if VERBOSE_LEVEL>1
        os << "n_mv: " << n_mv <<", sc_size: " << sc_size << ", it. for RAND: "<< cnt_permut_vec<<"\n";
#endif
    }
    //      r = permVec[0:nonsing_size-1]     (singular DOFs)
    //      s = permVec[nonsing_size:end-1]   (non-singular DOFs)
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
    for (eslocal i = 0;i<sc_size;i++) fix_dofs[i]=permVec[nonsing_size + i] + offset;
    K_rr.getSubDiagBlockmatrix(K_modif,K_rr,i_start, nonsing_size);
    if (check_nonsing!=0){
      double lmx_K_rr;
      cond_of_regular_part = K_rr.MatCondNumb(K_rr,"K_rr",plot_n_first_n_last_eigenvalues,&lmx_K_rr,100);

#if VERBOSE_LEVEL>2
        os << "cond of regular part = "<< cond_of_regular_part <<"\n" << std::flush ;
#endif
    }
//
    cnt_iter_check_nonsing++;
  }

  delete [] I_row_indices_p;
  delete [] J_col_indices_p;


#if VERBOSE_LEVEL>0
//3 - after singular test
  time1 = omp_get_wtime();
  elapsed_secs[3] = (time1 - begin_time) ;
#endif
               //                                               |

//
  K_rs.getSubBlockmatrix_rs(K_modif,K_rs,i_start, nonsing_size,j_start,sc_size);
#if VERBOSE_LEVEL>0
//4 - creation of block K_rs
  time1 = omp_get_wtime();
  elapsed_secs[4] = double(time1 - begin_time) ;
#endif
//
  SparseSolverCPU K_rr_solver;
  std::stringstream ss;
  bool SC_via_K_rr=true;
  if (SC_via_K_rr){
    K_rr_solver.ImportMatrix(K_rr);
    K_rr.Clear();
    ss << "get kerner from K -> rank: " << config::MPIrank;
    K_rr_solver.Factorization(ss.str());
    S.getSubDiagBlockmatrix(K_modif,S,nonsing_size,sc_size);
    SparseMatrix invKrrKrs = K_rs;
    K_rr_solver.SolveMat_Dense(invKrrKrs);
    SparseMatrix KsrInvKrrKrs;
	  KsrInvKrrKrs.MatMat(K_rs,'T',invKrrKrs);
    S.MatAddInPlace(KsrInvKrrKrs,'N',-1);
    S.RemoveLower();
  }
  else{
    SparseSolverCPU createSchur;
    // TODO PARDISO_SC provides factor K_rr.
    // if SC_via_K_rr=false,  factorization is made redundantly later.
    createSchur.ImportMatrix(K_modif);
    createSchur.Create_SC(S,sc_size,false);
    K_modif.Clear();
    createSchur.Clear();
  }
  S.type='S';
  S.ConvertCSRToDense(1);
#if VERBOSE_LEVEL>0
//5 - Schur complement created
  time1 = omp_get_wtime();
  elapsed_secs[5] = double(time1 - begin_time) ;
#endif
// EIGENVALUES AND EIGENVECTORS OF SCHUR COMPLEMENT
  char JOBZ = 'V';
  char UPLO = 'U';
  double *W = new double[S.cols];
  double *Z = new double[S.cols*S.cols];
  MKL_INT info;
  MKL_INT ldz = S.cols;
  info = LAPACKE_dspev (LAPACK_COL_MAJOR, JOBZ, UPLO, S.cols, &(S.dense_values[0]), W, Z, ldz);
  if (info){
    std::cout <<"info = " << info << " something wrong with Schur complement in SparseSolverCPU::generalIinverse\n";
  }
#if VERBOSE_LEVEL>0
//6 - Schur complement eigenvalues obtained
  time1 = omp_get_wtime();
  elapsed_secs[6] = double(time1 - begin_time) ;
#endif
// IDENTIFICATIONS OF ZERO EIGENVALUES
  eslocal defect_K_in;// R_s_cols;
  double ratio;
  eslocal itMax = twenty < S.rows ? twenty : S.rows ;
  for (eslocal i = itMax-1; i > 0;i--){
    ratio = fabs(W[i-1]/W[i]);
    if (ratio < jump_in_eigenvalues_alerting_singularity){
      defect_K_in=i;
    }
  }
#if VERBOSE_LEVEL>0
//7 - zero eigenvalues detection
  time1 = omp_get_wtime();
  elapsed_secs[7] = double(time1 - begin_time) ;
#endif
//
#if VERBOSE_LEVEL>1
  if (get_n_first_and_n_last_eigenvals_from_dense_S!=0){
    int i1i = get_n_first_and_n_last_eigenvals_from_dense_S;
    if (i1i>S.rows){i1i=S.rows;}
    os<<"eigenvals of S d{1:" << i1i << "} and d{" <<
         S.rows-get_n_first_and_n_last_eigenvals_from_dense_S+2 << ":"<< S.rows<< "}\n";

    for (eslocal i = 0 ; i < S.rows; i++){
      if (i < get_n_first_and_n_last_eigenvals_from_dense_S ||
            i > S.rows-get_n_first_and_n_last_eigenvals_from_dense_S){
        os<< i+1 <<":\t"<< W[i] << "\n";
      }
    }
  }
#endif
// --------------- CREATING KERNEL R_s FOR SINGULAR PART (SCHUR COMPLEMENT)
  SparseMatrix R_s;
  R_s.nnz  = defect_K_in*S.rows;
  R_s.dense_values.resize(R_s.nnz);
  R_s.rows = S.rows;
  R_s.cols = defect_K_in;
  R_s.type = 'G';
  eslocal cntR=0;
  for (eslocal j = 0; j < defect_K_in; j++){
    for (eslocal i = 0; i < R_s.rows; i++){
	    R_s.dense_values[cntR] = Z[j*R_s.rows + i];
      cntR++;
    }
  }
  R_s.ConvertDenseToCSR(0);
#if VERBOSE_LEVEL>0
//8 - R_s created
  time1 = omp_get_wtime();
  elapsed_secs[8] = double(time1 - begin_time) ;
#endif
// --------------- CREATING KERNEL R_r FOR NON-SINGULAR PART
	SparseMatrix R_r;
	R_r.MatMat(K_rs,'N',R_s);
  K_rs.Clear();
  if (!SC_via_K_rr){
    K_rr_solver.ImportMatrix(K_rr);
    K_rr.Clear();
//    std::stringstream ss;
    ss << "get kerner from K -> rank: " << config::MPIrank;
    K_rr_solver.Factorization(ss.str());
  }
  K_rr_solver.SolveMat_Sparse(R_r); // inv(K_rr)*K_rs*R_s
  K_rr_solver.Clear();

  R_r.ConvertCSRToDense(0);
  R_s.ConvertCSRToDense(0);
#if VERBOSE_LEVEL>0
//9 - R_r created (applied K_rr)
  time1 = omp_get_wtime();
  elapsed_secs[9] = double(time1 - begin_time) ;
#endif
               //                                               |
// --------------- CREATING WHOLE KERNEL Kplus_R = [ (R_r)^T (R_s)^T ]^T
  Kplus_R.rows = R_r.rows+R_s.rows;
  Kplus_R.cols = R_r.cols;
  Kplus_R.nnz  = Kplus_R.cols*Kplus_R.rows;
  Kplus_R.type = 'G';
	Kplus_R.dense_values.resize(Kplus_R.nnz);
  cntR=0;
  for (eslocal j = 0; j < Kplus_R.cols; j++){
    for (eslocal i = 0; i < R_r.rows; i++){
      if (diagonalScaling){
        di=K.CSR_V_values[K.CSR_I_row_indices[permVec[i]]-offset];
      }
      Kplus_R.dense_values[j*Kplus_R.rows + permVec[i]] = R_r.dense_values[j*R_r.rows + i]/sqrt(di);
      cntR++;
    }
    for (eslocal i = 0; i < R_s.rows; i++){
      if (diagonalScaling){
        di=K.CSR_V_values[K.CSR_I_row_indices[permVec[i+R_r.rows]]-offset];
      }
	    Kplus_R.dense_values[j*Kplus_R.rows + permVec[i+R_r.rows]] =-R_s.dense_values[j*R_s.rows + i]/sqrt(di);
      cntR++;
    }
  }
//
  Kplus_R.GramSchmidtOrtho();
  SEQ_VECTOR <eslocal > null_pivots;
  Kplus_R.getNullPivots(null_pivots);

#if VERBOSE_LEVEL>0
//10 - R - Gram Schmidt Orthogonalization
  time1 = omp_get_wtime();
  elapsed_secs[10] = double(time1 - begin_time) ;
#endif
               //                                               |

  // norm of product K*R: second matrix has to be in dense format!!!
//11 - max(eig(K))
#if VERBOSE_LEVEL>1
  double lmx_K;
  K.MatCondNumb(K,"K_singular",plot_n_first_n_last_eigenvalues,&lmx_K,100);
  os << std::scientific;
  os << "max(eig(K)):     " << lmx_K << "\n";
  os << "max(diag(K)):    " << rho << "\n";
  double tmp = K.getNorm_K_R(K,Kplus_R);
  *norm_KR_d_pow_2 = (tmp*tmp)/(lmx_K*lmx_K);
  *defect_d = Kplus_R.cols;
  os << "defect(K):       " << Kplus_R.cols <<"\n";
  os << "norm_KR:         " << sqrt(*norm_KR_d_pow_2) <<"\n";
#endif
               //                                               |
  time1 = omp_get_wtime();
  elapsed_secs[11] = double(time1 - begin_time) ;

//
  Kplus_R.ConvertDenseToCSR(0);
//
//
  if (diagonalRegularization){
    eslocal tmp_int0;
    if (d_sub!=-1) {
	    regMat.rows = K.rows;
	    regMat.cols = K.cols;
	    regMat.type = 'S';
	    regMat.nnz= null_pivots.size();

      regMat.I_row_indices.resize(regMat.nnz);
      regMat.J_col_indices.resize(regMat.nnz);
      regMat.V_values.resize(regMat.nnz);
    }
    for (eslocal i = 0; i < null_pivots.size(); i++){
      tmp_int0=K.CSR_I_row_indices[null_pivots[i]-offset]-offset;
      K.CSR_V_values[tmp_int0]+=rho;
      if (d_sub!=-1) {
        regMat.I_row_indices[i] = null_pivots[i];
        regMat.J_col_indices[i] = null_pivots[i];
        regMat.V_values[i]      = rho ;
      }
    }
  }
  else{
    SparseMatrix N;
    if (use_null_pivots_or_s_set){
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
    SparseSolverCPU NtN;
    NtN.ImportMatrix(NtN_Mat);
    NtN_Mat.Clear();
    std::stringstream sss;
    sss << "get kernel from K -> rank: " << config::MPIrank;
    NtN.Factorization(ss.str());
    NtN.SolveMat_Sparse(Nt);
    NtN.Clear();
    NtN_Mat.MatMat(N,'N',Nt);
    NtN_Mat.MatScale(rho);
    NtN_Mat.RemoveLower();
    K.MatAddInPlace (NtN_Mat,'N', 1);
    // IF d_sub == -1, it is GGt0 of cluster and regMat is no need
    if (d_sub!=-1)
    {
      regMat=NtN_Mat;
      regMat.ConvertToCOO(1);
    }
  }
//  K.printMatCSR("K_regularized");
//  K.MatCondNumb(K,"K_regularized",plot_n_first_n_last_eigenvalues);


  delete [] W;
  delete [] Z;

  double end_time = omp_get_wtime();
//12 - Total time in kernel detection
  elapsed_secs[12] = double(end_time - begin_time) ;
  std::cout << "Total time in kernel detection:                 " << elapsed_secs[12] << "[s] \n";

#if VERBOSE_LEVEL>0
  os << std::fixed;
  os << "allocation of vectors, copying of matrix:       " << elapsed_secs[0]                 << "[s] \n";
  os << "diagonal scaling:                               " << elapsed_secs[1]-elapsed_secs[0] << "[s] \n";
  os << "before singular test:                           " << elapsed_secs[2]-elapsed_secs[1] << "[s] \n";
  os << "after singular test:                            " << elapsed_secs[3]-elapsed_secs[2] << "[s] \n";
  os << "creation of block K_rs:                         " << elapsed_secs[4]-elapsed_secs[3] << "[s] \n";
  os << "Schur complement created:                       " << elapsed_secs[5]-elapsed_secs[4] << "[s] \n";
  os << "Schur complement eigenvalues obtained:          " << elapsed_secs[6]-elapsed_secs[5] << "[s] \n";
  os << "zero eigenvalues detection:                     " << elapsed_secs[7]-elapsed_secs[6] << "[s] \n";
  os << "R_s created:                                    " << elapsed_secs[8]-elapsed_secs[7] << "[s] \n";
  os << "R_r created (applied K_rr):                     " << elapsed_secs[9]-elapsed_secs[8] << "[s] \n";
  os << "R - Gram Schmidt Orthogonalization:             " << elapsed_secs[10]-elapsed_secs[9] << "[s] \n";
  os << "max(eig(K)):                                    " << elapsed_secs[11]-elapsed_secs[10] << "[s] \n";
  os << "Total time in kernel detection:                 " << elapsed_secs[12]<< "[s] \n";
  os.close();
#endif


}

