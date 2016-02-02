
#include "utils.h"

//#include "mkl.h"
//#include <string>
//#include <sstream>
//#include <iostream>
//#include <vector>
//#include <fstream>
//#include <algorithm>
//
//#include <math.h>
//#include <stack>
//#include <ctime>

using std::string; 
using std::endl; 
using std::left;
using std::streampos; 
 
#pragma once 

class FEM_Assembler
{

public:
	FEM_Assembler(string directory_path, string domain_name, int domain_global_index, int use_dynamic_1_no_dynamic_0);
	//FEM_Assembler();
	//~FEM_Assembler();

	int meshType;

	SEQ_VECTOR <SEQ_VECTOR <double> > all_coordinates;
	SEQ_VECTOR <SEQ_VECTOR <int> >	  all_elements;
	SEQ_VECTOR <int>			      elementType;
	int USE_DYNAMIC;

	void LoadData(string directory_path, string domain_name, int domain_global_index); 
	void Clear(); 
	void assemble_matrix(SparseMatrix & K, SEQ_VECTOR <double>  & f_global);
	void assemble_matrix(SparseMatrix & K, SparseMatrix & M, SEQ_VECTOR <double>  & f_global);

private: 

	double LoadCoordinates(string filename, SEQ_VECTOR <SEQ_VECTOR <double> > & coordinates );
	double LoadElements(string filename, SEQ_VECTOR <SEQ_VECTOR <int> > & elements );
	double LoadBinVectorInt(string filename, SEQ_VECTOR <int> & Vector);

	void brick8_basis       (SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N, SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <int> & mapVecB, SEQ_VECTOR <int> & mapVecN); 
	void tetrahedra10_basis (SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N, SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <int> & mapVecB, SEQ_VECTOR <int> & mapVecN);
	
	void prism6_basis (SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N , SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <int> & mapVecB);
	void pyramid5_basis (SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N ,SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <int> & mapVecB);
	void brick20_basis (SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N ,SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <int> & mapVecB); 
	void prisma15_basis (SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N ,SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <int> & mapVecB);
	void pyramid13_basis (SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N ,SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <int> & mapVecB);

	
	void brick8_Elasticity      (SEQ_VECTOR <double> & fe, SEQ_VECTOR <double> & Ke, SEQ_VECTOR <double> & Me, SEQ_VECTOR <int> & mapVecB, SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N, SEQ_VECTOR <int> & mapVecN, SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <double> & coordinates, SEQ_VECTOR <double> & inertia, double  ex, double  mi); 
	void tetrahedra10_Elasticity(SEQ_VECTOR <double> & fe, SEQ_VECTOR <double> & Ke, SEQ_VECTOR <double> & Me, SEQ_VECTOR <int> & mapVecB, SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N, SEQ_VECTOR <int> & mapVecN, SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <double> & coordinates, SEQ_VECTOR <double> & inertia, double  ex, double  mi);
	
	void pyramid5_Elasticity(SEQ_VECTOR <double> & fe, SEQ_VECTOR <double> & Ke, SEQ_VECTOR <int> & mapVecB, SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N ,SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <double> & coordinates, SEQ_VECTOR <double> & inertia, double  ex, double  mi);
	void prisma6_Elasticity(SEQ_VECTOR <double> & fe, SEQ_VECTOR <double> & Ke, SEQ_VECTOR <int> & mapVecB, SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N ,SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <double> & coordinates, SEQ_VECTOR <double> & inertia, double  ex, double  mi);
	void pyramid13_Elasticity(SEQ_VECTOR <double> & fe, SEQ_VECTOR <double> & Ke, SEQ_VECTOR <int> & mapVecB, SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N ,SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <double> & coordinates, SEQ_VECTOR <double> & inertia, double  ex, double  mi);
	void prisma15_Elasticity(SEQ_VECTOR <double> & fe, SEQ_VECTOR <double> & Ke, SEQ_VECTOR <int> & mapVecB, SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N ,SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <double> & coordinates, SEQ_VECTOR <double> & inertia, double  ex, double  mi);
	void brick20_Elasticity(SEQ_VECTOR <double> & fe, SEQ_VECTOR <double> & Ke, SEQ_VECTOR <int> & mapVecB, SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N ,SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <double> & coordinates, SEQ_VECTOR <double> & inertia, double  ex, double  mi);
};
