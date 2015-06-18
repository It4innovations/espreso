#include "SparseMatrix.h"
#include "FEM_Assembler.h"

typedef std::pair<int,double> IV_elem;	
bool unique_comparator (  const IV_elem& l,  const IV_elem& r) { 
	if (l.first == r.first) {
	}
	return l.first == r.first; 
}
bool comparator (  const IV_elem& l,  const IV_elem& r) { 
	return l.first < r.first; 
}

FEM_Assembler::FEM_Assembler(string directory_path, string domain_name, int domain_global_index, int use_dynamic_1_no_dynamic_0) {
	USE_DYNAMIC = use_dynamic_1_no_dynamic_0;
	FEM_Assembler::LoadData(directory_path, domain_name, domain_global_index); 
}

void FEM_Assembler::LoadData(string directory_path, string domain_name, int domain_global_index) {

	string path = directory_path + domain_name; 

	string coordinates_file;
	string elements_file;
	string elementsType_file;
	string meshType_file; 

	coordinates_file  = string(path) + string("COORDINATES"); 
	elements_file     = string(path) + string("ELEMENTS"); 
	elementsType_file = string(path) + string("elementType"); ///////OPRAVIT ELEMENTY NA MATSOL
	meshType_file     = string(path) + string("emeshType");


	int tres; 
	tres = LoadCoordinates(coordinates_file, all_coordinates); 
	if (tres == -1) {
		// !!! POZOR - works only for CUBE !!!!
		tres = LoadCoordinates( string(directory_path) + string("COORDINATES"), all_coordinates ); 

		SEQ_VECTOR <SEQ_VECTOR <double> > correct_coordinates; 
		tres = LoadCoordinates( string(directory_path) + string("Correct_COORDINATES"), correct_coordinates);

		for (int i = 0; i < all_coordinates.size(); i++) {
			for (int j = 0; j < all_coordinates[i].size(); j++) {
				all_coordinates[i][j] = all_coordinates[i][j] + correct_coordinates[domain_global_index - 1][j]; // domains numbering is from 1 
			}
		}


	}
	
	tres = LoadElements(elements_file, all_elements); 
	if (tres == -1) {
		// !!! POZOR - works only for CUBE !!!!
		tres = LoadElements(string(directory_path) + string("ELEMENTS"), all_elements); 
	}

	tres = LoadBinVectorInt(elementsType_file, elementType);
	if (tres == -1) {
		tres = LoadBinVectorInt(string(directory_path) + string("elementType"), elementType);
	}

	SEQ_VECTOR<int> meshType_v(1,0); 
	tres = LoadBinVectorInt(meshType_file, meshType_v);
	if (tres == -1) {
		tres = LoadBinVectorInt(string(directory_path) + string("emeshType"), meshType_v);
	}
	meshType = meshType_v[0]; 

	
}

double FEM_Assembler::LoadCoordinates(string filename, SEQ_VECTOR <SEQ_VECTOR <double> > & coordinates ) {

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

		coordinates.reserve(rows); // pozor
		coordinates.resize(rows);

		for (int r_i = 0; r_i < rows; r_i++) {
			SEQ_VECTOR <double> tmp_coor (cols); 
			in.read((char*) &tmp_coor [0], cols*sizeof(double));
			//coordinates.push_back(tmp_coor); 
			coordinates[r_i] = tmp_coor;
		}

		in.close(); 
		
		return 0; 
	} else {
		//cout << "File " << filename << " not found ! " << endl; 
		return -1; 
	}

}

double FEM_Assembler::LoadElements(string filename, SEQ_VECTOR <SEQ_VECTOR <int> > & elements ) {

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
		elements.resize(rows); 
		for (int r_i = 0; r_i < rows; r_i++) {
			elements[r_i].resize(cols); 
		}
		
		SEQ_VECTOR <int> tmp_coor (cols);
		for (int r_i = 0; r_i < rows; r_i++) {
			in.read((char*) &elements[r_i][0], cols*sizeof(int));
			//in.read((char*) &tmp_coor [0], cols*sizeof(int));
			//elements.push_back(tmp_coor); 
		}

		in.close(); 

		return 0; 
	} else {
		//cout << "File " << filename << " not found ! " << endl; 
		return -1; 
	}
}

double FEM_Assembler::LoadBinVectorInt(string filename, SEQ_VECTOR <int> & Vector) {

	// open the file:
	streampos fileSize;
	ifstream file (filename.c_str(), std::ios::binary);

	if ( file.is_open() ) {

		// get its size:
		file.seekg(0, std::ios::end);
		fileSize = file.tellg();
		file.seekg(0, std::ios::beg);

		// read the data:
		Vector.resize(fileSize / sizeof(int)); 

		file.read((char*) &Vector[0], fileSize);

		file.close();
		return 0; 
	} else {
		//cout << "File " << filename << " not found ! " << endl; 
		return -1; 
	}
}

void FEM_Assembler::brick8_basis (SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N , SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <int> & mapVecB, SEQ_VECTOR <int> & mapVecN) {

	int gp_length  = 8; 
	dN.resize(gp_length);
	N.resize(gp_length);

	int CsQ_length = 24;
	int CsQ_cols   = 3; 
	int CsQ_rows   = 8; 
	int dN_length  = 24;

	double CsQ_scale = 0.577350269189626;
	double CsQ[] = {-1,  -1,  -1,
		-1,  -1,   1,
		-1,   1,  -1,
		-1,   1,   1,
		1,  -1,  -1,
		1,  -1,   1,
		1,   1,  -1,
		1,   1,   1};


	for (int i = 0; i < CsQ_length; i++) {
		CsQ[i] = CsQ[i] * CsQ_scale; 
	}

	double WeighFactor_temp[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	WeighFactor.assign(WeighFactor_temp, WeighFactor_temp + 8);

	for (int gp_index = 0; gp_index < gp_length; gp_index++) {
		double r = CsQ[gp_index * CsQ_cols + 0];
		double s = CsQ[gp_index * CsQ_cols + 1];
		double t = CsQ[gp_index * CsQ_cols + 2];

		///dN contains [dNr,dNs,dNt]
		dN[gp_index].resize(dN_length);

		// dNr - derivation of basis function
		dN[gp_index][0] = 0.125*(-(1-s)*(1-t)); 
		dN[gp_index][1] = 0.125*( (1-s)*(1-t));
		dN[gp_index][2] = 0.125*( (1+s)*(1-t));
		dN[gp_index][3] = 0.125*(-(1+s)*(1-t));
		dN[gp_index][4] = 0.125*(-(1-s)*(1+t)); 
		dN[gp_index][5] = 0.125*( (1-s)*(1+t));
		dN[gp_index][6] = 0.125*( (1+s)*(1+t));
		dN[gp_index][7] = 0.125*(-(1+s)*(1+t));

		// dNs - derivation of basis function
		dN[gp_index][8]  = 0.125*(-(1-r)*(1-t)); 
		dN[gp_index][9]  = 0.125*(-(1+r)*(1-t));
		dN[gp_index][10] = 0.125*( (1+r)*(1-t));
		dN[gp_index][11] = 0.125*( (1-r)*(1-t));
		dN[gp_index][12] = 0.125*(-(1-r)*(1+t));
		dN[gp_index][13] = 0.125*(-(1+r)*(1+t));
		dN[gp_index][14] = 0.125*( (1+r)*(1+t));
		dN[gp_index][15] = 0.125*( (1-r)*(1+t));

		// dNt - derivation of basis function
		dN[gp_index][16] = 0.125*(-(1-r)*(1-s));
		dN[gp_index][17] = 0.125*(-(1+r)*(1-s)); 
		dN[gp_index][18] = 0.125*(-(1+r)*(1+s));
		dN[gp_index][19] = 0.125*(-(1-r)*(1+s));
		dN[gp_index][20] = 0.125*( (1-r)*(1-s));
		dN[gp_index][21] = 0.125*( (1+r)*(1-s));
		dN[gp_index][22] = 0.125*( (1+r)*(1+s));
		dN[gp_index][23] = 0.125*( (1-r)*(1+s));

		// basis function
		N[gp_index].resize(gp_length); 
		N[gp_index][0] = 0.125 * (1-r)*(1-s)*(1-t);
		N[gp_index][1] = 0.125 * (r+1)*(1-s)*(1-t);
		N[gp_index][2] = 0.125 * (r+1)*(s+1)*(1-t);
		N[gp_index][3] = 0.125 * (1-r)*(s+1)*(1-t);
		N[gp_index][4] = 0.125 * (1-r)*(1-s)*(t+1);
		N[gp_index][5] = 0.125 * (r+1)*(1-s)*(t+1);
		N[gp_index][6] = 0.125 * (r+1)*(s+1)*(t+1);
		N[gp_index][7] = 0.125 * (1-r)*(s+1)*(t+1);

	}

	int mapVec_temp[] = {1, 4, 7, 10, 13, 16, 19, 22, 74, 77, 80, 83, 86, 89, 92, 95, 123, 126, 129, 132, 135, 138,
		141, 144, 26, 29, 32, 35, 38, 41, 44, 47, 73, 76, 79, 82, 85, 88, 91, 94, 99, 102, 105, 108,
		111, 114, 117, 120, 51, 54, 57, 60, 63, 66, 69, 72, 98, 101, 104, 107, 110, 113, 116, 119, 121, 124,
		127, 130, 133, 136, 139, 142};

	mapVecB.assign(mapVec_temp,mapVec_temp+72);
	
	// dynamic 
	int mapVec_temp2[] = {1, 4, 7, 10, 13, 16, 19, 22, 26, 29, 32, 35, 38, 41, 44, 47, 51, 54, 57, 60, 63, 66, 69, 72};
	mapVecN.assign(mapVec_temp2,mapVec_temp2+24);


}

void FEM_Assembler::prism6_basis (SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N , SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <int> & mapVecB) {

	int gp_length  = 9; 
	dN.resize(gp_length);
	N.resize(gp_length);

	int CsQ_length = 27;
	int CsQ_cols   = 3; 
	int CsQ_rows   = 9; 
	int dN_length  = 18;


	double CsQ[] = { 1/6, 1/6, -0.774596669241483,
		4/6, 1/6, -0.774596669241483,
		1/6, 4/6, -0.774596669241483,
		1/6, 1/6, 0,
		4/6, 1/6, 0,
		1/6, 4/6, 0,
		1/6, 1/6, 0.774596669241483,
		4/6, 1/6, 0.774596669241483,
		1/6, 4/6, 0.774596669241483};

	double WeighFactor_temp[] = {5/54, 5/54, 5/54,  8/54, 8/54, 8/54, 5/54, 5/54, 5/54};
	WeighFactor.assign(WeighFactor_temp, WeighFactor_temp + 9);

	for (int gp_index = 0; gp_index < gp_length; gp_index++) {
		double r = CsQ[gp_index * CsQ_cols + 0];
		double s = CsQ[gp_index * CsQ_cols + 1];
		double t = CsQ[gp_index * CsQ_cols + 2];

		dN[gp_index].resize(dN_length);
		// dNr - derivation of basis function
		dN[gp_index][0] = (  t/2 - 1/2); 
		dN[gp_index][1] = (  1/2 - t/2);
		dN[gp_index][2] =  0;
		dN[gp_index][3] = (- t/2 - 1/2);
		dN[gp_index][4] = (  t/2 + 1/2);
		dN[gp_index][5] =  0;

		// dNs - derivation of basis function
		dN[gp_index][6] =  ( t/2 - 1/2);
		dN[gp_index][7] =   0;
		dN[gp_index][8] =  ( 1/2 - t/2);
		dN[gp_index][9] =  (- t/2 - 1/2);
		dN[gp_index][10] =  0;
		dN[gp_index][11] = (t/2 + 1/2);

		// dNt - derivation of basis function
		dN[gp_index][12] =  ((r/2 + s/2) - 1/2);
		dN[gp_index][13] =  (-r/2);
		dN[gp_index][14] =  (-s/2);
		dN[gp_index][15] =  (1/2 - s/2 - r/2);
		dN[gp_index][16] =  (r/2);
		dN[gp_index][17] =  (s/2);

		N[gp_index].resize(dN_length);
		// basis function
		N[gp_index][0] = 0.5 * ((1 - r - s) * (1-t));
		N[gp_index][1] = 0.5 * (r * (1 - t));
		N[gp_index][2] = 0.5 * (s * (1 - t));
		N[gp_index][3] = 0.5 * ((1 - r - s) * (1 + t));
		N[gp_index][4] = 0.5 * (r * (t + 1));
		N[gp_index][5] = 0.5 * (s * (t + 1));


		int mapVec_temp[] = {1, 4, 7, 10, 13, 16, 56, 59, 62, 65, 68, 71, 93, 96, 99, 102, 105, 108, 20, 23, 26, 29, 32, 35, 55, 58, 61,
			64, 67, 70, 75, 78, 81, 84, 87, 90, 39, 42, 45, 48, 51, 54, 74, 77, 80, 83, 86, 89, 91, 94, 97, 100, 103, 106};

		mapVecB.assign(mapVec_temp, mapVec_temp + 54);

	}

}

void FEM_Assembler::pyramid5_basis (SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N ,SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <int> & mapVecB) {

	int gp_length  = 8; 
	dN.resize(gp_length);
	N.resize(gp_length);

	int CsQ_length = 24;
	int CsQ_cols   = 3; 
	int CsQ_rows   = 8; 
	int dN_length  = 15;

	double CsQ_scale = 0.577350269189626;
	double CsQ[] = {-1,  -1,  -1,
		-1,  -1,   1,
		-1,   1,  -1,
		-1,   1,   1,
		1,  -1,  -1,
		1,  -1,   1,
		1,   1,  -1,
		1,   1,   1};

	for (int i = 0; i < CsQ_length; i++) {
		CsQ[i] = CsQ[i] * CsQ_scale; 
	}

	double WeighFactor_temp[] = {1, 1, 1, 1, 1, 1, 1, 1};
	WeighFactor.assign(WeighFactor_temp, WeighFactor_temp + 8);

	for (int gp_index = 0; gp_index < gp_length; gp_index++) {
		double r = CsQ[gp_index * CsQ_cols + 0];
		double s = CsQ[gp_index * CsQ_cols + 1];
		double t = CsQ[gp_index * CsQ_cols + 2];

		dN[gp_index].resize(dN_length);
		// dNr - derivation of basis function
		dN[gp_index][0] = 0.125 * (-(1-s)*(1-t)); 
		dN[gp_index][1] = 0.125 * ((1-s)*(1-t));
		dN[gp_index][2] = 0.125 * ((1+s)*(1-t));
		dN[gp_index][3] = 0.125 * (-(1+s)*(1-t));
		dN[gp_index][4] = 0;

		// dNs - derivation of basis function
		dN[gp_index][5] = 0.125 * (-(1-r)*(1-t));
		dN[gp_index][6] = 0.125 * (-(1+r)*(1-t));
		dN[gp_index][7] = 0.125 * ((1+r)*(1-t));
		dN[gp_index][8] = 0.125 * ((1-r)*(1-t));
		dN[gp_index][9] = 0;

		// dNt - derivation of basis function
		dN[gp_index][10] = 0.125 * (-(1-r)*(1-s));
		dN[gp_index][11] = 0.125 * (-(1+r)*(1-s));
		dN[gp_index][12] = 0.125 * (-(1+r)*(1+s));
		dN[gp_index][13] = 0.125 * (-(1-r)*(1+s));
		dN[gp_index][14] = 0.125 * (4);

		N[gp_index].resize(dN_length);
		// basis function
		N[gp_index][0] = 0.125 * ((1-r)*(1-s)*(1-t));
		N[gp_index][1] = 0.125 * ((r+1)*(1-s)*(1-t));
		N[gp_index][2] = 0.125 * ((r+1)*(s+1)*(1-t));
		N[gp_index][3] = 0.125 * ((1-r)*(s+1)*(1-t));
		N[gp_index][4] = 0.125 * ( 4*(1+t));


		int mapVec_temp[] = { 1, 4, 7, 10, 13, 47, 50, 53, 56, 59, 78, 81, 84, 87, 90, 17, 20, 23, 26, 29, 46, 49, 52, 55, 58, 63, 66,
			69, 72, 75, 33, 36, 39, 42, 45, 62, 65, 68, 71, 74, 76, 79, 82, 85, 88};

		mapVecB.assign(mapVec_temp, mapVec_temp + 45);

	}
}

void FEM_Assembler::brick20_basis (SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N ,SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <int> & mapVecB) {

	int gp_length  = 14; 
	dN.resize(gp_length);
	N.resize(gp_length);

	int CsQ_length = 42;
	int CsQ_cols   = 3; 
	int CsQ_rows   = 14; 
	int dN_length  = 60;

	double CsQ_scale_1 = 0.758786910639329;
	double CsQ_scale_2 = 0.795822425754222;
	double CsQ[] = {-1.0,    -1.0,    -1.0,
		1.0,    -1.0,    -1.0,
		1.0,     1.0,    -1.0,
		-1.0,     1.0,    -1.0,
		-1.0,    -1.0,     1.0,
		1.0,    -1.0,     1.0,
		1.0,     1.0,     1.0,
		-1.0,     1.0,     1.0,	
		0.0,		0.0,  -1.0,
		0.0,	   -1.0,   0.0,
		1.0,		0.0,   0.0,
		0.0,		1.0,   0.0,
		-1.0,		0.0,   0.0,
		0.0,		0.0,   1.0};

	for (int i = 0; i < 24; i++) {
		CsQ[i] = CsQ[i] * CsQ_scale_1; 
	}
	for (int i = 24; i < 42; i++) {
		CsQ[i] = CsQ[i] * CsQ_scale_2; 
	}

	double WF_scale_1 =0.335180055401662;
	double WF_scale_2 =0.886426592797784;

	double WeighFactor_temp[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	for (int i = 0; i < 8; i++) {
		WeighFactor_temp[i] = WeighFactor_temp[i] * WF_scale_1; 
	}
	for (int i = 8; i < 14; i++) {
		WeighFactor_temp[i] = WeighFactor_temp[i] * WF_scale_2; 
	}

	WeighFactor.assign(WeighFactor_temp, WeighFactor_temp + 14);

	for (int gp_index = 0; gp_index < gp_length; gp_index++) {
		double r = CsQ[gp_index * CsQ_cols + 0];
		double s = CsQ[gp_index * CsQ_cols + 1];
		double t = CsQ[gp_index * CsQ_cols + 2];


		dN[gp_index].resize(dN_length);
		// dNr - derivation of basis function
		dN[gp_index][0] =   ((s - 1.0)*(t - 1.0)*(r + s + t + 2.0))/8.0 + ((r - 1.0)*(s - 1.0)*(t - 1.0))/8.0;
		dN[gp_index][1] =   ((r + 1.0)*(s - 1.0)*(t - 1.0))/8.0 - ((s - 1.0)*(t - 1.0)*(s - r + t + 2.0))/8.0;
		dN[gp_index][2] = - ((s + 1.0)*(t - 1.0)*(r + s - t - 2.0))/8.0 - ((r + 1.0)*(s + 1.0)*(t - 1.0))/8.0;
		dN[gp_index][3] = - ((s + 1.0)*(t - 1.0)*(r - s + t + 2.0))/8.0 - ((r - 1.0)*(s + 1.0)*(t - 1.0))/8.0;
		dN[gp_index][4] = - ((s - 1.0)*(t + 1.0)*(r + s - t + 2.0))/8.0 - ((r - 1.0)*(s - 1.0)*(t + 1.0))/8.0;
		dN[gp_index][5] = - ((s - 1.0)*(t + 1.0)*(r - s + t - 2.0))/8.0 - ((r + 1.0)*(s - 1.0)*(t + 1.0))/8.0;
		dN[gp_index][6] =   ((s + 1.0)*(t + 1.0)*(r + s + t - 2.0))/8.0 + ((r + 1.0)*(s + 1.0)*(t + 1.0))/8.0;
		dN[gp_index][7] =   ((s + 1.0)*(t + 1.0)*(r - s - t + 2.0))/8.0 + ((r - 1.0)*(s + 1.0)*(t + 1.0))/8.0;
		dN[gp_index][8] = - (r*(s - 1.0)*(t - 1.0))/2.0;
		dN[gp_index][9] =   ((pow(s,2.0) - 1.0)*(t - 1.0))/4.0;
		dN[gp_index][10] =  (r*(s + 1.0)*(t - 1.0))/2.0;
		dN[gp_index][11] = -((pow(s,2.0) - 1.0)*(t - 1.0))/4.0;
		dN[gp_index][12] =  (r*(s - 1.0)*(t + 1.0))/2.0;
		dN[gp_index][13] = -((pow(s,2.0) - 1.0)*(t + 1.0))/4.0;
		dN[gp_index][14] = -(r*(s + 1.0)*(t + 1.0))/2.0;
		dN[gp_index][15] =  ((pow(s,2.0) - 1.0)*(t + 1.0))/4.0;
		dN[gp_index][16] = -((pow(t,2.0) - 1.0)*(s - 1.0))/4.0;
		dN[gp_index][17] =  ((pow(t,2.0) - 1.0)*(s - 1.0))/4.0;
		dN[gp_index][18] = -((pow(t,2.0) - 1.0)*(s + 1.0))/4.0;
		dN[gp_index][19] =  ((pow(t,2.0) - 1.0)*(s + 1.0))/4.0;


		// dNs - derivation of basis function
		dN[gp_index][20] = ((r - 1.0)*(t - 1.0)*(r + s + t + 2.0))/8.0 + ((r - 1.0)*(s - 1.0)*(t - 1.0))/8.0;
		dN[gp_index][21] =  - ((r + 1.0)*(t - 1.0)*(s - r + t + 2.0))/8.0 - ((r + 1.0)*(s - 1.0)*(t - 1.0))/8.0;
		dN[gp_index][22] = - ((r + 1.0)*(t - 1.0)*(r + s - t - 2.0))/8.0 - ((r + 1.0)*(s + 1.0)*(t - 1.0))/8.0;
		dN[gp_index][23] = ((r - 1.0)*(s + 1.0)*(t - 1.0))/8.0 - ((r - 1.0)*(t - 1.0)*(r - s + t + 2.0))/8.0;
		dN[gp_index][24] = - ((r - 1.0)*(t + 1.0)*(r + s - t + 2.0))/8.0 - ((r - 1.0)*(s - 1.0)*(t + 1.0))/8.0;
		dN[gp_index][25] = ((r + 1.0)*(s - 1.0)*(t + 1.0))/8.0 - ((r + 1.0)*(t + 1.0)*(r - s + t - 2.0))/8.0;
		dN[gp_index][26] =  ((r + 1.0)*(t + 1.0)*(r + s + t - 2.0))/8.0 + ((r + 1.0)*(s + 1.0)*(t + 1.0))/8.0;
		dN[gp_index][27] = ((r - 1.0)*(t + 1.0)*(r - s - t + 2.0))/8.0 - ((r - 1.0)*(s + 1.0)*(t + 1.0))/8.0;
		dN[gp_index][28] = -((pow(r,2.0) - 1.0)*(t - 1.0))/4.0;
		dN[gp_index][29] = (s*(r + 1.0)*(t - 1.0))/2.0;
		dN[gp_index][30] = ((pow(r,2.0) - 1.0)*(t - 1.0))/4.0;
		dN[gp_index][31] = -(s*(r - 1.0)*(t - 1.0))/2.0;
		dN[gp_index][32] = ((pow(r,2.0) - 1.0)*(t + 1.0))/4.0;
		dN[gp_index][33] = -(s*(r + 1.0)*(t + 1.0))/2.0;
		dN[gp_index][34] = -((pow(r,2.0) - 1.0)*(t + 1.0))/4.0;
		dN[gp_index][35] = (s*(r - 1.0)*(t + 1.0))/2.0;
		dN[gp_index][36] = -((pow(t,2.0) - 1.0)*(r - 1.0))/4.0;
		dN[gp_index][37] = ((pow(t,2.0) - 1.0)*(r + 1.0))/4.0;
		dN[gp_index][38] = -((pow(t,2.0) - 1.0)*(r + 1.0))/4.0;
		dN[gp_index][39] = ((pow(t,2.0) - 1.0)*(r - 1.0))/4.0;

		// dNt - derivation of basis function
		dN[gp_index][40] = ((r - 1.0)*(s - 1.0)*(r + s + t + 2.0))/8.0 + ((r - 1.0)*(s - 1.0)*(t - 1.0))/8.0;
		dN[gp_index][41] = - ((r + 1.0)*(s - 1.0)*(s - r + t + 2.0))/8.0 - ((r + 1.0)*(s - 1.0)*(t - 1.0))/8.0;
		dN[gp_index][42] = ((r + 1.0)*(s + 1.0)*(t - 1.0))/8.0 - ((r + 1.0)*(s + 1.0)*(r + s - t - 2.0))/8.0;
		dN[gp_index][43] = - ((r - 1.0)*(s + 1.0)*(r - s + t + 2.0))/8.0 - ((r - 1.0)*(s + 1.0)*(t - 1.0))/8.0;
		dN[gp_index][44] = ((r - 1.0)*(s - 1.0)*(t + 1.0))/8.0 - ((r - 1)*(s - 1.0)*(r + s - t + 2.0))/8.0;
		dN[gp_index][45] = - ((r + 1.0)*(s - 1.0)*(r - s + t - 2.0))/8.0 - ((r + 1.0)*(s - 1.0)*(t + 1.0))/8.0;
		dN[gp_index][46] = ((r + 1.0)*(s + 1.0)*(r + s + t - 2.0))/8.0 + ((r + 1.0)*(s + 1.0)*(t + 1.0))/8.0;
		dN[gp_index][47] = ((r - 1.0)*(s + 1.0)*(r - s - t + 2.0))/8.0 - ((r - 1.0)*(s + 1.0)*(t + 1.0))/8.0;
		dN[gp_index][48] = -((pow(r,2.0) - 1.0)*(s - 1.0))/4.0; 
		dN[gp_index][49] = ((pow(s,2.0) - 1.0)*(r + 1.0))/4.0;
		dN[gp_index][50] = ((pow(r,2.0) - 1.0)*(s + 1.0))/4.0;
		dN[gp_index][51] = -((pow(s,2.0) - 1.0)*(r - 1.0))/4.0;
		dN[gp_index][52] = ((pow(r,2.0) - 1.0)*(s - 1.0))/4.0;
		dN[gp_index][53] = -((pow(s,2.0) - 1.0)*(r + 1.0))/4.0;
		dN[gp_index][54] = -((pow(r,2.0) - 1.0)*(s + 1.0))/4.0;
		dN[gp_index][55] = ((pow(s,2.0) - 1.0)*(r - 1.0))/4.0;
		dN[gp_index][56] = -(t*(r - 1.0)*(s - 1.0))/2.0;
		dN[gp_index][57] = (t*(r + 1.0)*(s - 1.0))/2.0;
		dN[gp_index][58] = -(t*(r + 1.0)*(s + 1.0))/2.0;
		dN[gp_index][59] = (t*(r - 1.0)*(s + 1.0))/2.0;

		N[gp_index].resize(dN_length);
		// basis function
		N[gp_index][0] = 0.125 * ((1.0-r)*(1-s)*(1.0-t)*(-2.0-r-s-t));
		N[gp_index][1] = 0.125 * ((r+1.0)*(1-s)*(1.0-t)*(r-s-t-2.0));
		N[gp_index][2] = 0.125 * ((r+1.0)*(s+1.0)*(1.0-t)*(r+s-t-2.0));
		N[gp_index][3] = 0.125 * ((1.0-r)*(s+1.0)*(1.0-t)*(-r+s-t-2.0));
		N[gp_index][4] = 0.125 * ((1.0-r)*(1-s)*(t+1.0)*(-r-s+t-2.0));
		N[gp_index][5] = 0.125 * ((r+1.0)*(1-s)*(t+1.0)*(r-s+t-2.0));
		N[gp_index][6] = 0.125 * ((r+1.0)*(s+1.0)*(t+1.0)*(r+s+t-2.0));
		N[gp_index][7] = 0.125 * ((1.0-r)*(s+1.0)*(t+1.0)*(-r+s+t-2.0));
		N[gp_index][8] = 0.25 * ((1.0-pow(r,2.0))*(1.0-s)*(1.0-t));
		N[gp_index][9] = 0.25 * ((1.0+r)*(1.0-pow(s,2.0))*(1.0-t));
		N[gp_index][10] = 0.25 * ((1.0-pow(r,2.0))*(1.0+s)*(1.0-t));
		N[gp_index][11] = 0.25 * ((1.0-r)*(1.0-pow(s,2.0))*(1.0-t));
		N[gp_index][12] = 0.25 * ((1.0-pow(r,2.0))*(1.0-s)*(1.0+t));
		N[gp_index][13] = 0.25 * ((1.0+r)*(1.0-pow(s,2.0))*(1.0+t));
		N[gp_index][14] = 0.25 * ((1.0-pow(r,2.0))*(1.0+s)*(1.0+t));
		N[gp_index][15] = 0.25 * ((1.0-r)*(1.0-pow(s,2.0))*(1.0+t));
		N[gp_index][16] = 0.25 * ((1.0-r)*(1.0-s)*(1.0-pow(t,2.0)));
		N[gp_index][17] = 0.25 * ((1.0+r)*(1.0-s)*(1.0-pow(t,2.0)));
		N[gp_index][18] = 0.25 * ((1.0+r)*(1.0+s)*(1.0-pow(t,2.0)));
		N[gp_index][19] = 0.25 * ((1.0-r)*(1.0+s)*(1.0-pow(t,2.0)));



		int mapVec_temp[] = { 1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 52, 55, 58, 182, 185, 188, 191, 194, 197, 200,
			203, 206, 209, 212, 215, 218, 221, 224, 227, 230, 233, 236, 239, 303, 306, 309, 312, 315, 318, 321, 324, 327, 330, 333, 336, 339, 342,
			345, 348, 351, 354, 357, 360,  62,  65,  68,  71,  74,  77,  80,  83,  86,  89,  92,  95,  98, 101, 104, 107, 110, 113, 116, 119, 181,
			184, 187, 190, 193, 196, 199, 202, 205, 208, 211, 214, 217, 220, 223, 226, 229, 232, 235, 238, 243, 246, 249, 252, 255, 258, 261, 264,
			267, 270, 273, 276, 279, 282, 285, 288, 291, 294, 297, 300, 123, 126, 129, 132, 135, 138, 141, 144, 147, 150, 153, 156, 159, 162, 165,
			168, 171, 174, 177, 180, 242, 245, 248, 251, 254, 257, 260, 263, 266, 269, 272, 275, 278, 281, 284, 287, 290, 293, 296, 299, 301, 304,
			307, 310, 313, 316, 319, 322, 325, 328, 331, 334, 337, 340, 343, 346, 349, 352, 355, 358};

		mapVecB.assign(mapVec_temp, mapVec_temp + 180);

	}
}

void FEM_Assembler::prisma15_basis (SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N ,SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <int> & mapVecB) {

	int gp_length  = 9; 
	dN.resize(gp_length);
	N.resize(gp_length);

	int CsQ_length = 27;
	int CsQ_cols   = 3; 
	int CsQ_rows   = 9; 
	int dN_length  = 45;


	double CsQ[] = { 1.0/6.0, 1.0/6.0, -0.774596669241483,
		                4.0/6.0, 1.0/6.0, -0.774596669241483,
		                1.0/6.0, 4.0/6.0, -0.774596669241483,
		                1.0/6.0, 1.0/6.0,  0.774596669241483,
		                4.0/6.0, 1.0/6.0,  0.774596669241483,
		                1.0/6.0, 4.0/6.0,  0.774596669241483,
		                1.0/6.0, 1.0/6.0,  0.0,
		                4.0/6.0, 1.0/6.0,  0.0,
		                1.0/6.0, 4.0/6.0,  0.0};


	double WF_scale_1 =5.0/54.0;
	double WF_scale_2 =8.0/54.0;

	double WeighFactor_temp[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };


	for (int i = 0; i < 6; i++) {
		WeighFactor_temp[i] = WeighFactor_temp[i] * WF_scale_1; 
	}
	for (int i = 6; i < 9; i++) {
		WeighFactor_temp[i] = WeighFactor_temp[i] * WF_scale_2; 
	}



	WeighFactor.assign(WeighFactor_temp, WeighFactor_temp + 9);


	for (int gp_index = 0; gp_index < gp_length; gp_index++) {
		double r = CsQ[gp_index * CsQ_cols + 0];
		double s = CsQ[gp_index * CsQ_cols + 1];
		double t = CsQ[gp_index * CsQ_cols + 2];

		dN[gp_index].resize(dN_length);

		// dNr - derivation of basis function
		dN[gp_index][0] = - (t - 1.0)*(r + s - 1.0) - ((t - 1.0)*(2.0*r + 2.0*s + t))/2.0;
		dN[gp_index][1] = ((t - 1.0)*(t - 2.0*r + 2.0))/2.0 - r*(t - 1.0);
		dN[gp_index][2] = 0;
		dN[gp_index][3] = ((t + 1.0)*(2.0*r + 2.0*s - t))/2.0 + (t + 1.0)*(r + s - 1.0);
		dN[gp_index][4] = r*(t + 1.0) + ((t + 1.0)*(2.0*r + t - 2.0))/2.0;
		dN[gp_index][5] = 0.0;
		dN[gp_index][6] = 2.0*(t - 1.0)*(r + s - 1.0) + 2.0*r*(t - 1.0);
		dN[gp_index][7] = (-2.0)*s*(t - 1.0);
		dN[gp_index][8] = 2.0*s*(t - 1.0);
		dN[gp_index][9] = - 2.0*(t + 1.0)*(r + s - 1.0) - 2.0*r*(t + 1.0);
		dN[gp_index][10] = 2.0*s*(t + 1.0);
		dN[gp_index][11] = (-2.0)*s*(t + 1.0);
		dN[gp_index][12] = pow(t,2.0) - 1.0;
		dN[gp_index][13] = 1.0 - pow(t,2.0);
		dN[gp_index][14] = 0.0;


		// dNs - derivation of basis function
		dN[gp_index][15] =  - (t - 1.0)*(r + s - 1.0) - ((t - 1.0)*(2.0*r + 2.0*s + t))/2.0;
		dN[gp_index][16] =  0.0;
		dN[gp_index][17] =  ((t - 1.0)*(t - 2.0*s + 2.0))/2.0 - s*(t - 1.0);
		dN[gp_index][18] =  ((t + 1.0)*(2.0*r + 2.0*s - t))/2.0 + (t + 1.0)*(r + s - 1.0);
		dN[gp_index][19] =  0.0;
		dN[gp_index][20] =  s*(t + 1.0) + ((t + 1.0)*(2.0*s + t - 2.0))/2.0;
		dN[gp_index][21] =  2.0*r*(t - 1.0);
		dN[gp_index][22] =  (-2.0)*r*(t - 1.0);
		dN[gp_index][23] =  2.0*(t - 1.0)*(r + s - 1.0) + 2.0*s*(t - 1.0);
		dN[gp_index][24] =  (-2.0)*r*(t + 1.0);
		dN[gp_index][25] =  2.0*r*(t + 1.0);
		dN[gp_index][26] =  - 2.0*(t + 1.0)*(r + s - 1.0) - 2.0*s*(t + 1.0);
		dN[gp_index][27] =  pow(t,2.0) - 1.0;
		dN[gp_index][28] =  0.0;
		dN[gp_index][29] =  1.0 - pow(t,2.0);

		// dNt - derivation of basis function
		dN[gp_index][30] = - ((r + s - 1.0)*(2.0*r + 2.0*s + t))/2.0 - ((t - 1.0)*(r + s - 1.0))/2.0;
		dN[gp_index][31] =  (r*(t - 2.0*r + 2.0))/2.0 + (r*(t - 1.0))/2.0;
		dN[gp_index][32] =  (s*(t - 2.0*s + 2.0))/2.0 + (s*(t - 1.0))/2.0;
		dN[gp_index][33] = ((r + s - 1.0)*(2.0*r + 2.0*s - t))/2.0 - ((t + 1.0)*(r + s - 1.0))/2.0;
		dN[gp_index][34] = (r*(2.0*r + t - 2.0))/2.0 + (r*(t + 1.0))/2.0;
		dN[gp_index][35] = (s*(2.0*s + t - 2.0))/2.0 + (s*(t + 1.0))/2.0;
		dN[gp_index][36] =  2.0*r*(r + s - 1.0);
		dN[gp_index][37] =  (-2.0)*r*s;
		dN[gp_index][38] =  2.0*s*(r + s - 1.0);
		dN[gp_index][39] =  (-2.0)*r*(r + s - 1.0);
		dN[gp_index][40] =  2.0*r*s;
		dN[gp_index][41] = (-2.0)*s*(r + s - 1.0);
		dN[gp_index][42] = 2.0*t*(r + s - 1.0);
		dN[gp_index][43] = (-2.0)*r*t;
		dN[gp_index][44] = (-2.0)*s*t;

		N[gp_index].resize(dN_length);
		// basis function
		N[gp_index][0] = -(1.0-r-s)*(1.0-t)*(2.0*r+2.0*s+t)/2.0;
		N[gp_index][1] =  r*(1.0-t)*(2.0*r-t-2.0)/2.0;
		N[gp_index][2] =  s*(1.0-t)*(2.0*s-t-2.0)/2.0;
		N[gp_index][3] =  -(1.0-r-s)*(1.0+t)*(2.0*r+2.0*s-t)/2.0;
		N[gp_index][4] =  r*(t+1.0)*(2.0*r+t-2.0)/2.0;
		N[gp_index][5] =  s*(t+1.0)*(2.0*s+t-2.0)/2.0;
		N[gp_index][6] =  2.0*r*(1.0-r-s)*(1.0-t);
		N[gp_index][7] =  2.0*r*s*(1.0-t);
		N[gp_index][8] =  2.0*s*(1.0-r-s)*(1.0-t);
		N[gp_index][9] =  2.0*r*(1.0-r-s)*(1.0+t);
		N[gp_index][10] = 2.0*r*s*(1.0+t);
		N[gp_index][11] = 2.0*s*(1.0-r-s)*(1.0+t);
		N[gp_index][12] = (1.0-r-s)*(1.0-pow(t,2.0));
		N[gp_index][13] = r*(1.0-pow(t,2.0));
		N[gp_index][14] = s*(1.0-pow(t,2.0));

		int mapVec_temp[] = { 1,   4,   7,  10,  13,  16,  19,  22,  25,  28,  31,  34,  37,  40,  43, 137, 140, 143, 146, 149, 152, 155, 158, 161, 164, 167, 170,
			173, 176, 179, 228, 231, 234, 237, 240, 243, 246, 249, 252, 255, 258, 261, 264, 267, 270,  47,  50, 53, 56,  59,  62,  65, 68,  71,
			74,  77,  80, 83,  86,  89, 136, 139, 142, 145, 148, 151, 154, 157, 160, 163, 166, 169, 172, 175, 178, 183, 186, 189, 192, 195, 198,
			201, 204, 207, 210, 213, 216, 219, 222, 225,  93,  96,  99, 102, 105, 108, 111, 114, 117, 120, 123, 126, 129, 132, 135, 182, 185, 188,
			191, 194, 197, 200, 203, 206, 209, 212, 215, 218, 221, 224, 226, 229, 232, 235, 238, 241, 244, 247, 250, 253, 256, 259, 262, 265,   268};

		mapVecB.assign(mapVec_temp, mapVec_temp + 135);

	}
}

void FEM_Assembler::pyramid13_basis (SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N ,SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <int> & mapVecB) {

	int gp_length  = 8; 
	dN.resize(gp_length);
	N.resize(gp_length);

	int CsQ_length = 24;
	int CsQ_cols   = 3; 
	int CsQ_rows   = 8; 
	int dN_length  = 39;


	double CsQ_scale = 0.577350269189626;
	double CsQ[] = {-1.0,  -1.0,  -1.0,
              		-1.0,  -1.0,   1.0,
              		-1.0,   1.0,  -1.0,
              		-1.0,   1.0,   1.0,
                	 1.0,  -1.0,  -1.0,
                	 1.0,  -1.0,   1.0,
                	 1.0,   1.0,  -1.0,
                	 1.0,   1.0,   1.0};


	for (int i = 0; i < CsQ_length; i++) {
		CsQ[i] = CsQ[i] * CsQ_scale; 
	}

	double WeighFactor_temp[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,};
	WeighFactor.assign(WeighFactor_temp, WeighFactor_temp + 8);


	for (int gp_index = 0; gp_index < gp_length; gp_index++) {
		double r = CsQ[gp_index * CsQ_cols + 0];
		double s = CsQ[gp_index * CsQ_cols + 1];
		double t = CsQ[gp_index * CsQ_cols + 2];

		dN[gp_index].resize(dN_length);

		// dNr - derivation of basis function
		dN[gp_index][0] = - (t/8.0 - 1.0/8.0)*(s - 1.0)*(r*(t/2.0 - 0.5) + s*(t/2.0 - 0.5) - 1.0) - (t/2.0 - 0.5)*(t/8.0 - 1.0/8.0)*(r - 1.0)*(s - 1.0);
		dN[gp_index][1] =  - (t/8.0 - 1.0/8.0)*(s - 1.0)*(r*(t/2.0 - 0.5) - s*(t/2.0 - 0.5) + 1.0) - (t/2.0 - 0.5)*(t/8.0 - 1.0/8.0)*(r + 1.0)*(s - 1.0);
		dN[gp_index][2] =  (t/8.0 - 1.0/8.0)*(s + 1.0)*(r*(t/2.0 - 0.5) + s*(t/2.0 - 0.5) + 1.0) + (t/2.0 - 0.5)*(t/8.0 - 1.0/8.0)*(r + 1.0)*(s + 1.0);
		dN[gp_index][3] =  (t/2.0 - 0.5)*(t/8.0 - 1.0/8.0)*(r - 1.0)*(s + 1.0) - (t/8.0 - 1.0/8.0)*(s + 1.0)*(s*(t/2.0 - 0.5) - r*(t/2.0 - 0.5) + 1.0);
		dN[gp_index][4] =  0.0;
		dN[gp_index][5] =  r*((t/2.0 - 0.5)*(t/2.0 - 0.5))*(s - 1.0);
		dN[gp_index][6] =  -((pow(s,2.0) - 1.0)*(t/2.0 - 0.5)*(t/2.0 - 0.5))/2.0;
		dN[gp_index][7] =  -r*((t/2.0 - 0.5)*(t/2.0 - 0.5))*(s + 1.0);
		dN[gp_index][8] =  ((pow(s,2.0) - 1)*(t/2.0 - 0.5)*(t/2.0 - 0.5))/2.0;
		dN[gp_index][9] =  -(t/2.0 - 0.5)*(t/2.0 + 0.5)*(s - 1.0);
		dN[gp_index][10] = (t/2.0 - 0.5)*(t/2.0 + 0.5)*(s - 1.0);
		dN[gp_index][11] = -(t/2.0 - 0.5)*(t/2.0 + 0.5)*(s + 1.0);
		dN[gp_index][12] =  (t/2.0 - 0.5)*(t/2.0 + 0.5)*(s + 1.0);

		// dNs - derivation of basis function
		dN[gp_index][13] = - (t/8.0 - 1.0/8.0)*(r - 1.0)*(r*(t/2.0 - 1.0/2.0) + s*(t/2.0 - 1.0/2.0) - 1.0) - (t/2.0 - 1.0/2.0)*(t/8.0 - 1.0/8.0)*(r - 1.0)*(s - 1.0);
		dN[gp_index][14] = (t/2.0 - 1.0/2.0)*(t/8.0 - 1.0/8.0)*(r + 1.0)*(s - 1.0) - (t/8.0 - 1.0/8.0)*(r + 1.0)*(r*(t/2.0 - 1.0/2.0) - s*(t/2.0 - 1.0/2.0) + 1.0);
		dN[gp_index][15] =  (t/8.0 - 1.0/8.0)*(r + 1.0)*(r*(t/2.0 - 1.0/2.0) + s*(t/2.0 - 1.0/2.0) + 1.0) + (t/2.0 - 1.0/2.0)*(t/8.0 - 1.0/8.0)*(r + 1.0)*(s + 1.0);
		dN[gp_index][16] =  - (t/8.0 - 1.0/8.0)*(r - 1.0)*(s*(t/2.0 - 1.0/2.0) - r*(t/2.0 - 1.0/2.0) + 1.0) - (t/2.0 - 1.0/2.0)*(t/8.0 - 1.0/8.0)*(r - 1.0)*(s + 1.0);
		dN[gp_index][17] =  0.0;
		dN[gp_index][18] =   ((pow(r,2.0) - 1.0)*(t/2.0 - 0.5)*(t/2.0 - 1.0/2.0))/2.0;
		dN[gp_index][19] =    -s*(t/2.0 - 1.0/2.0)*(t/2.0 - 0.5)*(r + 1.0);
		dN[gp_index][20] =   -((pow(r,2.0) - 1.0)*(t/2.0 - 0.5)*(t/2.0 - 1.0/2.0))/2.0;
		dN[gp_index][21] =   s*(t/2.0 - 0.5)*(t/2.0 - 0.5)*(r - 1.0);
		dN[gp_index][22] =  -(t/2.0 - 0.5)*(t/2.0 + 0.5)*(r - 1.0);
		dN[gp_index][23] =  (t/2.0 - 0.5)*(t/2.0 + 0.5)*(r + 1.0);
		dN[gp_index][24] = -(t/2.0 - 0.5)*(t/2.0 + 0.5)*(r + 1.0);
		dN[gp_index][25] =  (t/2.0 - 0.5)*(t/2.0 + 0.5)*(r - 1.0);

		// dNt - derivation of basis function
		dN[gp_index][26] =  - ((r - 1.0)*(s - 1.0)*(r*(t/2.0 - 0.5) + s*(t/2.0 - 0.5) - 1.0))/8.0 - (r/2.0 + s/2.0)*(t/8.0 - 1.0/8.0)*(r - 1.0)*(s - 1.0);
		dN[gp_index][27] =   - ((r + 1.0)*(s - 1.0)*(r*(t/2.0 - 0.5) - s*(t/2.0 - 0.5) + 1.0))/8.0 - (r/2.0 - s/2.0)*(t/8.0 - 1.0/8.0)*(r + 1.0)*(s - 1.0);
		dN[gp_index][28] =  ((r + 1.0)*(s + 1.0)*(r*(t/2.0 - 0.5) + s*(t/2.0 - 0.5) + 1.0))/8.0 + (r/2.0 + s/2.0)*(t/8.0 - 1.0/8.0)*(r + 1.0)*(s + 1.0);
		dN[gp_index][29] =  (r/2.0 - s/2.0)*(t/8.0 - 1.0/8.0)*(r - 1.0)*(s + 1.0) - ((r - 1.0)*(s + 1.0)*(s*(t/2.0 - 0.5) - r*(t/2.0 - 0.5) + 1.0))/8.0;
		dN[gp_index][30] =   t + 0.5;
		dN[gp_index][31] =   ((pow(r,2.0) - 1.0)*(t/2.0 - 0.5)*(s - 1.0))/2.0;
		dN[gp_index][32] =   -((pow(s,2.0) - 1.0)*(t/2.0 - 0.5)*(r + 1.0))/2.0;
		dN[gp_index][33] =  -((pow(r,2.0) - 1.0)*(t/2.0 - 0.5)*(s + 1.0))/2.0;
		dN[gp_index][34] =  ((pow(s,2.0) - 1.0)*(t/2.0 - 0.5)*(r - 1.0))/2.0;
		dN[gp_index][35] =  ((t/2.0 - 0.5)*(r + s - r*s - 1.0))/2.0 + ((t/2.0 + 0.5)*(r + s - r*s - 1.0))/2.0;
		dN[gp_index][36] =   - ((t/2.0 - 0.5)*(r - s - r*s + 1.0))/2.0 - ((t/2.0 + 0.5)*(r - s - r*s + 1.0))/2.0;
		dN[gp_index][37] =   - ((t/2.0 - 0.5)*(r + s + r*s + 1.0))/2.0 - ((t/2.0 + 0.5)*(r + s + r*s + 1.0))/2.0;
		dN[gp_index][38] =   ((t/2.0 - 0.5)*(r - s + r*s - 1.0))/2.0 + ((t/2.0 + 0.5)*(r - s + r*s - 1.0))/2.0;


		N[gp_index].resize(dN_length);
		// basis function
		N[gp_index][0] = ((0.5*(1.0-t))/4.0)*( (1.0-r)*(1.0-s)*(-1.0-(0.5*(1.0-t))*r-(0.5*(1.0-t))*s));
		N[gp_index][1] =    ((0.5*(1.0-t))/4.0)*( (1.0+r)*(1.0-s)*(-1.0+(0.5*(1-t))*r-(0.5*(1.0-t))*s));
		N[gp_index][2] =   ((0.5*(1.0-t))/4.0)*( (1.0+r)*(1.0+s)*(-1.0+(0.5*(1.0-t))*r+(0.5*(1.0-t))*s));
		N[gp_index][3] =  ((0.5*(1.0-t))/4.0)*( (1.0-r)*(1.0+s)*(-1.0-(0.5*(1.0-t))*r+(0.5*(1.0-t))*s));
		N[gp_index][4] =   (1.0-(0.5*(1.0-t)))*(1.0-2.0*(0.5*(1.0-t)));
		N[gp_index][5] =  (((0.5*(1.0-t))*(0.5*(1.0-t)))/2.0)*(1.0-s)*(1.0-pow(r,2.0));
		N[gp_index][6] =  (((0.5*(1.0-t))*(0.5*(1.0-t)))/2.0)*(1.0+r)*(1.0-pow(s,2.0));
		N[gp_index][7] =   (((0.5*(1.0-t))*(0.5*(1.0-t)))/2.0)*(1.0+s)*(1.0-pow(r,2.0));
		N[gp_index][8] =   (((0.5*(1.0-t))*(0.5*(1.0-t)))/2.0)*(1.0-r)*(1.0-pow(s,2.0));
		N[gp_index][9] =  (0.5*(1.0-t))*(1.0-(0.5*(1.0-t)))*(1.0-r-s+r*s);
		N[gp_index][10] =  (0.5*(1.0-t))*(1.0-(0.5*(1.0-t)))*(1.0+r-s-r*s);
		N[gp_index][11] = (0.5*(1.0-t))*(1.0-(0.5*(1.0-t)))*(1.0+r+s+r*s);
		N[gp_index][12] =  (0.5*(1.0-t))*(1.0-(0.5*(1.0-t)))*(1.0-r+s-r*s);

		int mapVec_temp[] = {  1,   4,   7,  10,  13,  16,  19,  22,  25,  28,  31,  34, 37,  119, 122, 125, 128, 131, 134, 137, 140, 143, 146, 149, 152, 155, 198,
			201, 204, 207, 210, 213, 216, 219, 222, 225, 228, 231, 234, 41,   44,  47,  50,  53, 56, 59,  62,  65,  68,  71,  74,  77, 118, 121,
			124, 127, 130, 133, 136, 139, 142, 145, 148, 151, 154, 159, 162, 165, 168, 171, 174, 177, 180, 183, 186, 189, 192, 195,  81,  84,  87,
			90,  93,  96,  99, 102, 105, 108, 111, 114, 117, 158, 161, 164, 167, 170,173, 176, 179, 182, 185, 188, 191, 194, 196, 199, 202, 205,
			208, 211, 214, 217, 220, 223, 226, 229, 232};

		mapVecB.assign(mapVec_temp, mapVec_temp + 117);

	}
}

void FEM_Assembler::tetrahedra10_basis (SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N, SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <int> & mapVecB, SEQ_VECTOR <int> & mapVecN) {
//void FEM_Assembler::tetrahedra10_basis (PAR_VECTOR<PAR_VECTOR <double> > & dN , vector<vector <double> > & N ,vector <double> & WeighFactor, vector <int> & mapVecB) {

	int gp_length  = 4; 
	dN.resize(gp_length);
	N.resize(gp_length);

	int CsQ_length = 16;
	int CsQ_cols   = 4; 
	int CsQ_rows   = 4; 
	int dN_length  = 40;

	double CsQ_scale = 2.236067977499790;
	double CsQ[] = {(5.0-CsQ_scale)/20.0,   (5.0-CsQ_scale)/20.0,   (5.0-CsQ_scale)/20.0,   (5.0+3.0*CsQ_scale)/20.0,
		(5.0-CsQ_scale)/20.0,   (5.0-CsQ_scale)/20.0,   (5.0+3.0*CsQ_scale)/20.0, (5.0-CsQ_scale)/20.0,
		(5.0-CsQ_scale)/20.0,   (5.0+3.0*CsQ_scale)/20.0, (5.0-CsQ_scale)/20.0,   (5.0-CsQ_scale)/20.0,
		(5.0+3.0*CsQ_scale)/20.0, (5.0-CsQ_scale)/20.0,   (5.0-CsQ_scale)/20.0,   (5.0-CsQ_scale)/20.0};

	double WeighFactor_scale = 1.0/24.0;
	double WeighFactor_temp[] = {1.0, 1.0, 1.0, 1.0};

	for (int i = 0; i < 4; i++) {
		WeighFactor_temp[i] = WeighFactor_temp[i] * WeighFactor_scale; 
	}

	WeighFactor.assign(WeighFactor_temp, WeighFactor_temp + 4);


	for (int gp_index = 0; gp_index < gp_length; gp_index++) {
		double r = CsQ[gp_index * CsQ_cols + 0];
		double s = CsQ[gp_index * CsQ_cols + 1];
		double t = CsQ[gp_index * CsQ_cols + 2];
		double q = CsQ[gp_index * CsQ_cols + 3];

		dN[gp_index].resize(dN_length);

		// dNq - derivation of basis function    dNq=[ 4*q - 1, 0, 0, 0, 4*r, 0, 4*s, 4*t, 0, 0];
		dN[gp_index][0] = 4.0*q - 1;
		dN[gp_index][1] = 0.0;
		dN[gp_index][2] = 0.0;
		dN[gp_index][3] = 0.0;
		dN[gp_index][4] = 4.0*r;
		dN[gp_index][5] = 0.0;
		dN[gp_index][6] = 4*s;
		dN[gp_index][7] = 4*t;
		dN[gp_index][8] = 0.0;
		dN[gp_index][9] = 0.0;


		// dNr - derivation of basis function  dNr=[0,4*r - 1, 0, 0, 4*q, 4*s, 0, 0, 4*t, 0];
		dN[gp_index][10] = 0.0;
		dN[gp_index][11] = 4.0*r - 1.0;
		dN[gp_index][12] = 0.0;
		dN[gp_index][13] = 0.0;
		dN[gp_index][14] = 4.0*q;
		dN[gp_index][15] = 4.0*s;
		dN[gp_index][16] = 0.0;
		dN[gp_index][17] = 0.0;
		dN[gp_index][18] = 4.0*t;
		dN[gp_index][19] = 0.0;

		// dNs - derivation of basis function     dNs=[ 0, 0, 4*s - 1, 0, 0,4*r,  4*q,0,  0, 4*t];
		dN[gp_index][20] = 0.0;
		dN[gp_index][21] = 0.0;
		dN[gp_index][22] = 4.0*s-1.0;
		dN[gp_index][23] = 0.0;
		dN[gp_index][24] = 0.0;
		dN[gp_index][25] = 4.0*r;
		dN[gp_index][26] = 4.0*q;
		dN[gp_index][27] = 0.0;
		dN[gp_index][28] = 0.0;
		dN[gp_index][29] = 4.0*t;

		// dNt - derivation of basis function    dNt=[ 0, 0, 0,4*t - 1,  0, 0, 0, 4*q, 4*r, 4*s];  
		dN[gp_index][30] = 0.0;
		dN[gp_index][31] = 0.0;
		dN[gp_index][32] = 0.0;
		dN[gp_index][33] = 4.0*t - 1.0;
		dN[gp_index][34] = 0.0;
		dN[gp_index][35] = 0.0;
		dN[gp_index][36] = 0.0;
		dN[gp_index][37] = 4.0*q;
		dN[gp_index][38] = 4.0*r;
		dN[gp_index][39] = 4.0*s;


		N[gp_index].resize(dN_length);
		// basis function
		N[gp_index][0] = (2.0*q-1.0)*q;
		N[gp_index][1] = (2.0*r-1.0)*r;
		N[gp_index][2] = (2.0*s-1.0)*s;
		N[gp_index][3] = (2.0*t-1.0)*t;
		N[gp_index][4] = 4.0*q*r;
		N[gp_index][5] = 4.0*r*s;
		N[gp_index][6] = 4.0*q*s;
		N[gp_index][7] = 4.0*q*t;
		N[gp_index][8] = 4.0*r*t;
		N[gp_index][9] = 4.0*s*t;


		int mapVec_temp[] = { 1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 92, 95, 98, 101, 104, 107, 110, 113, 116, 119, 153, 156, 159, 162, 165, 168, 171,
			174, 177, 180,  32, 35, 38, 41, 44, 47, 50, 53, 56, 59, 91, 94, 97, 100, 103, 106, 109, 112, 115, 118, 123, 126, 129, 132,
			135, 138, 141, 144, 147, 150, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90, 122, 125, 128, 131, 134, 137, 140, 143, 146, 149, 151,
			154, 157, 160, 163, 166, 169, 172, 175, 178};

		mapVecB.assign(mapVec_temp, mapVec_temp + 90);

		// dynamic 
		int mapVec_temp2[] = {1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 32, 35, 38, 41, 44, 47, 50, 53, 56, 59, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90};
		mapVecN.assign(mapVec_temp2,mapVec_temp2+30);


	}
}

void FEM_Assembler::pyramid5_Elasticity(SEQ_VECTOR <double> & fe, SEQ_VECTOR <double> & Ke, SEQ_VECTOR <int> & mapVecB, SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N ,SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <double> & coordinates, SEQ_VECTOR <double> & inertia, double  ex, double  mi) {

	int gp_length  = 8; 

	fill(Ke.begin(),Ke.end(),0); 
	fill(fe.begin(),fe.end(),0); 

	SEQ_VECTOR <double> MatC (6*6, 0.0); 

	double E   = ex/((1+mi)*(1-2*mi));
	double mi2 = E*(1.0 - mi); 
	double mi3 = E*(0.5 - mi);
	mi = E * mi;

	MatC[0]  = mi2;		MatC[1]  = mi; 		MatC[2]  = mi;
	MatC[6]  = mi; 		MatC[7]  = mi2;		MatC[8]  = mi;
	MatC[12] = mi;		MatC[13] = mi;		MatC[14] = mi2;
	MatC[21] = mi3;		MatC[28] = mi3;		MatC[35] = mi3;

	char trans   = 'N'; 
	double alpha =  1; 
	double beta  =  0; 
	int m = 3; 
	int n = 3;
	int k = 8; 

	int m1 = 6; 
	int n1 = 6;
	int k1 = 15; 

	SEQ_VECTOR <double> MatJ (9, 0);
	SEQ_VECTOR <double> invJ (9, 0);
	SEQ_VECTOR <double> dND (15, 0);
	SEQ_VECTOR <double> CB (90, 0);
	SEQ_VECTOR <double> B (90, 0);
	Ke.resize(225);//Ke.resize(576);
	fe.resize(15);//fe.resize(24);




	for (int gp_index = 0; gp_index < gp_length; gp_index++) {

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, &dN[gp_index][0], k, &coordinates[0], n, beta, &MatJ[0], n);

		double detJ = fabs( MatJ[0] * MatJ[4] * MatJ[8] + 
			MatJ[1] * MatJ[5] * MatJ[6] + 
			MatJ[2] * MatJ[3] * MatJ[7] - 
			MatJ[2] * MatJ[4] * MatJ[6] -
			MatJ[1] * MatJ[3] * MatJ[8] -
			MatJ[0] * MatJ[5] * MatJ[7]);

		double detJx = 1/detJ;

		invJ[0] = detJx * (  MatJ[8] * MatJ[4] - MatJ[7] * MatJ[5] );
		invJ[1] = detJx * (- MatJ[8] * MatJ[1] + MatJ[7] * MatJ[2] );
		invJ[2] = detJx * (  MatJ[5] * MatJ[1] - MatJ[4] * MatJ[2] );
		invJ[3] = detJx * (- MatJ[8] * MatJ[3] + MatJ[6] * MatJ[5] );		
		invJ[4] = detJx * (  MatJ[8] * MatJ[0] - MatJ[6] * MatJ[2] );		
		invJ[5] = detJx * (- MatJ[5] * MatJ[0] + MatJ[3] * MatJ[2] );		
		invJ[6] = detJx * (  MatJ[7] * MatJ[3] - MatJ[6] * MatJ[4] );		
		invJ[7] = detJx * (- MatJ[7] * MatJ[0] + MatJ[6] * MatJ[1] );
		invJ[8] = detJx * (  MatJ[4] * MatJ[0] - MatJ[3] * MatJ[1] );

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, k, n, alpha, &invJ[0], n, &dN[gp_index][0], k, beta,&dND[0], k);

		for (int j = 0; j < 8; j++) {

			B[mapVecB[j]-1]    = dND[j];
			B[mapVecB[j+5]-1]  = dND[j];
			B[mapVecB[j+10]-1] = dND[j];


			B[mapVecB[j+15]-1] = dND[j+5];
			B[mapVecB[j+20]-1] = dND[j+5];
			B[mapVecB[j+25]-1] = dND[j+5];


			B[mapVecB[j+30]-1] = dND[j+10];
			B[mapVecB[j+35]-1] = dND[j+10];
			B[mapVecB[j+40]-1] = dND[j+10];

		}

		//C*B'
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m1, k1, n1, alpha, &MatC[0], m1, &B[0], k1, beta, &CB[0], m1);

		//Ke = Ke+(B*(C*B'))*dJ*WF;
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, k1, k1, n1, detJ*WeighFactor[gp_index], &B[0], k1, &CB[0], n1, alpha, &Ke[0], k1);

		for (int i = 0; i < 15; i++) {
			fe[i]= fe[i] + (detJ*WeighFactor[gp_index])*(N[gp_index][ i / 3 ] * inertia[ i % 3]);
		}

	}
}

void FEM_Assembler::prisma6_Elasticity(SEQ_VECTOR <double> & fe, SEQ_VECTOR <double> & Ke, SEQ_VECTOR <int> & mapVecB, SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N ,SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <double> & coordinates, SEQ_VECTOR <double> & inertia, double  ex, double  mi) {

	int gp_length  = 8; 

	fill(Ke.begin(),Ke.end(),0); 
	fill(fe.begin(),fe.end(),0); 

	SEQ_VECTOR <double> MatC (6*6, 0.0); 

	double E   = ex/((1+mi)*(1-2*mi));
	double mi2 = E*(1.0 - mi); 
	double mi3 = E*(0.5 - mi);
	mi = E * mi;

	MatC[0]  = mi2;		MatC[1]  = mi; 		MatC[2]  = mi;
	MatC[6]  = mi; 		MatC[7]  = mi2;		MatC[8]  = mi;
	MatC[12] = mi;		MatC[13] = mi;		MatC[14] = mi2;
	MatC[21] = mi3;		MatC[28] = mi3;		MatC[35] = mi3;

	char trans   = 'N'; 
	double alpha =  1; 
	double beta  =  0; 
	int m = 3; 
	int n = 3;
	int k = 8; 

	int m1 = 6; 
	int n1 = 6;
	int k1 = 18; 

	SEQ_VECTOR <double> MatJ (9, 0);
	SEQ_VECTOR <double> invJ (9, 0);
	SEQ_VECTOR <double> dND (18, 0);
	SEQ_VECTOR <double> CB (108, 0);
	SEQ_VECTOR <double> B (108, 0);
	Ke.resize(324);//Ke.resize(576);
	fe.resize(18);//fe.resize(24);




	for (int gp_index = 0; gp_index < gp_length; gp_index++) {

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, &dN[gp_index][0], k, &coordinates[0], n, beta, &MatJ[0], n);

		double detJ = fabs( MatJ[0] * MatJ[4] * MatJ[8] + 
			MatJ[1] * MatJ[5] * MatJ[6] + 
			MatJ[2] * MatJ[3] * MatJ[7] - 
			MatJ[2] * MatJ[4] * MatJ[6] -
			MatJ[1] * MatJ[3] * MatJ[8] -
			MatJ[0] * MatJ[5] * MatJ[7]);

		double detJx = 1/detJ;

		invJ[0] = detJx * (  MatJ[8] * MatJ[4] - MatJ[7] * MatJ[5] );
		invJ[1] = detJx * (- MatJ[8] * MatJ[1] + MatJ[7] * MatJ[2] );
		invJ[2] = detJx * (  MatJ[5] * MatJ[1] - MatJ[4] * MatJ[2] );
		invJ[3] = detJx * (- MatJ[8] * MatJ[3] + MatJ[6] * MatJ[5] );		
		invJ[4] = detJx * (  MatJ[8] * MatJ[0] - MatJ[6] * MatJ[2] );		
		invJ[5] = detJx * (- MatJ[5] * MatJ[0] + MatJ[3] * MatJ[2] );		
		invJ[6] = detJx * (  MatJ[7] * MatJ[3] - MatJ[6] * MatJ[4] );		
		invJ[7] = detJx * (- MatJ[7] * MatJ[0] + MatJ[6] * MatJ[1] );
		invJ[8] = detJx * (  MatJ[4] * MatJ[0] - MatJ[3] * MatJ[1] );

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, k, n, alpha, &invJ[0], n, &dN[gp_index][0], k, beta,&dND[0], k);

		for (int j = 0; j < 8; j++) {

			B[mapVecB[j]-1]    = dND[j];
			B[mapVecB[j+6]-1]  = dND[j];
			B[mapVecB[j+12]-1] = dND[j];


			B[mapVecB[j+18]-1] = dND[j+6];
			B[mapVecB[j+24]-1] = dND[j+6];
			B[mapVecB[j+30]-1] = dND[j+6];


			B[mapVecB[j+36]-1] = dND[j+12];
			B[mapVecB[j+42]-1] = dND[j+12];
			B[mapVecB[j+48]-1] = dND[j+12];

		}

		//C*B'
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m1, k1, n1, alpha, &MatC[0], m1, &B[0], k1, beta, &CB[0], m1);

		//Ke = Ke+(B*(C*B'))*dJ*WF;
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, k1, k1, n1, detJ*WeighFactor[gp_index], &B[0], k1, &CB[0], n1, alpha, &Ke[0], k1);

		for (int i = 0; i < 18; i++) {		
			fe[i]= fe[i] + (detJ*WeighFactor[gp_index])*(N[gp_index][ i / 3 ] * inertia[ i % 3]);
		}

	}
}

void FEM_Assembler::brick8_Elasticity(SEQ_VECTOR <double> & fe, SEQ_VECTOR <double> & Ke, SEQ_VECTOR <double> & Me, SEQ_VECTOR <int> & mapVecB, SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N , SEQ_VECTOR <int> & mapVecN, SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <double> & coordinates, SEQ_VECTOR <double> & inertia, double  ex, double  mi) {

	int gp_length  = 8; 


	SEQ_VECTOR <double> MatC (6*6, 0.0); 

	double E   = ex/((1+mi)*(1-2*mi));
	double mi2 = E*(1.0 - mi); 
	double mi3 = E*(0.5 - mi);
	mi = E * mi;

	MatC[0]  = mi2;		MatC[1]  = mi; 		MatC[2]  = mi;
	MatC[6]  = mi; 		MatC[7]  = mi2;		MatC[8]  = mi;
	MatC[12] = mi;		MatC[13] = mi;		MatC[14] = mi2;
	MatC[21] = mi3;		MatC[28] = mi3;		MatC[35] = mi3;

	char trans   = 'N'; 
	double alpha =  1; 
	double beta  =  0; 
	int m = 3;
	int n = 3;
	int k = 8; 

	int m1 = 6; 
	int n1 = 6;
	int k1 = 24; 

	SEQ_VECTOR <double> MatJ (9, 0);
	SEQ_VECTOR <double> invJ (9, 0);
	SEQ_VECTOR <double> dND (24, 0);
	SEQ_VECTOR <double> CB (144, 0);
	SEQ_VECTOR <double> B (144, 0);
	Ke.resize(576);//Ke.resize(576);
	Me.resize(576); 
	fe.resize(24);//fe.resize(24);


	fill(Ke.begin(),Ke.end(),0); 
	fill(Me.begin(),Me.end(),0); 
	fill(fe.begin(),fe.end(),0); 


	for (int gp_index = 0; gp_index < gp_length; gp_index++) {

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, &dN[gp_index][0], k, &coordinates[0], n, beta, &MatJ[0], n);

		double detJ = fabs( MatJ[0] * MatJ[4] * MatJ[8] + 
			MatJ[1] * MatJ[5] * MatJ[6] + 
			MatJ[2] * MatJ[3] * MatJ[7] - 
			MatJ[2] * MatJ[4] * MatJ[6] -
			MatJ[1] * MatJ[3] * MatJ[8] -
			MatJ[0] * MatJ[5] * MatJ[7]);

		double detJx = 1/detJ;

		invJ[0] = detJx * (  MatJ[8] * MatJ[4] - MatJ[7] * MatJ[5] );
		invJ[1] = detJx * (- MatJ[8] * MatJ[1] + MatJ[7] * MatJ[2] );
		invJ[2] = detJx * (  MatJ[5] * MatJ[1] - MatJ[4] * MatJ[2] );
		invJ[3] = detJx * (- MatJ[8] * MatJ[3] + MatJ[6] * MatJ[5] );		
		invJ[4] = detJx * (  MatJ[8] * MatJ[0] - MatJ[6] * MatJ[2] );		
		invJ[5] = detJx * (- MatJ[5] * MatJ[0] + MatJ[3] * MatJ[2] );		
		invJ[6] = detJx * (  MatJ[7] * MatJ[3] - MatJ[6] * MatJ[4] );		
		invJ[7] = detJx * (- MatJ[7] * MatJ[0] + MatJ[6] * MatJ[1] );
		invJ[8] = detJx * (  MatJ[4] * MatJ[0] - MatJ[3] * MatJ[1] );

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, k, n, alpha, &invJ[0], n, &dN[gp_index][0], k, beta,&dND[0], k);

		for (int j = 0; j < 8; j++) {

			B[mapVecB[j]-1]    = dND[j];
			B[mapVecB[j+8]-1]  = dND[j];
			B[mapVecB[j+16]-1] = dND[j];


			B[mapVecB[j+24]-1] = dND[j+8];
			B[mapVecB[j+32]-1] = dND[j+8];
			B[mapVecB[j+40]-1] = dND[j+8];


			B[mapVecB[j+48]-1] = dND[j+16];
			B[mapVecB[j+56]-1] = dND[j+16];
			B[mapVecB[j+64]-1] = dND[j+16];

		}

		//C*B'
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m1, k1, n1, alpha, &MatC[0], m1, &B[0], k1, beta, &CB[0], m1);

		//Ke = Ke+(B*(C*B'))*dJ*WF;
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, k1, k1, n1, detJ*WeighFactor[gp_index], &B[0], k1, &CB[0], n1, alpha, &Ke[0], k1);

		for (int i = 0; i < 24; i++) {		
			fe[i]= fe[i] + (detJ*WeighFactor[gp_index])*(N[gp_index][ i / 3 ] * inertia[ i % 3]);
		}


		// dynamic 

		SEQ_VECTOR <double> NN (72, 0.0);
	
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 8; j++) {
				NN[mapVecN[8*i+j]-1] = N[gp_index][j];
			}
		}

		// Me = Me + WF*(DENS*dJ)*(NN*NN');
		double dense = 7.85e-9;
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 24, 24, 3, dense*detJ*WeighFactor[gp_index], &NN[0], 24, &NN[0], 24, alpha, &Me[0], 24);
		//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 24, 24, 3, 1.0, &NN[0], 24, &NN[0], 24, alpha, &Me[0], 24);


		


	}
}

void FEM_Assembler::tetrahedra10_Elasticity(SEQ_VECTOR <double> & fe, SEQ_VECTOR <double> & Ke, SEQ_VECTOR <double> & Me, SEQ_VECTOR <int> & mapVecB, SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N, SEQ_VECTOR <int> & mapVecN, SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <double> & coordinates, SEQ_VECTOR <double> & inertia, double  ex, double  mi) {
//void FEM_Assembler::tetrahedra10_Elasticity(vector <double> & fe, vector <double> & Ke, vector <int> & mapVecB, vector<vector <double> > & dN , vector<vector <double> > & N ,vector <double> & WeighFactor ,vector <double> & coordinates, vector <double> & inertia, double  ex, double  mi) {

	int gp_length  = WeighFactor.size(); 

	int matrixSize = coordinates.size();


	SEQ_VECTOR <double> MatC (6*6, 0.0); 

	double E   = ex/((1.0+mi)*(1.0-2.0*mi));
	double mi2 = E*(1.0 - mi); 
	double mi3 = E*(0.5 - mi);
	mi = E * mi;

	MatC[0]  = mi2;		MatC[1]  = mi; 		MatC[2]  = mi;
	MatC[6]  = mi; 		MatC[7]  = mi2;		MatC[8]  = mi;
	MatC[12] = mi;		MatC[13] = mi;		MatC[14] = mi2;
	MatC[21] = mi3;		MatC[28] = mi3;		MatC[35] = mi3;

	char trans   = 'N'; 
	double alpha =  1.0; 
	double beta  =  0.0; 
	int m = 4; 
	int n = 3;
	int k = 10; 

	int m1 = 6; 
	int n1 = 6;
	int k1 = 30; 

	SEQ_VECTOR <double> MatJ (12, 0);
	SEQ_VECTOR <double> invJ (12, 0);
	SEQ_VECTOR <double> dND (30, 0);
	SEQ_VECTOR <double> CB (180, 0);
	SEQ_VECTOR <double> B (180, 0);
	Ke.resize(matrixSize*matrixSize);
	fe.resize(matrixSize);

	fill(Ke.begin(),Ke.end(),0); 
	fill(fe.begin(),fe.end(),0); 

	//dynamic
	Me.resize(matrixSize*matrixSize);
	fill(Me.begin(),Me.end(),0); 
	

	for (int gp_index = 0; gp_index < gp_length; gp_index++) {

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, &dN[gp_index][0], k, &coordinates[0],n, beta, &MatJ[0], n);



		invJ[0] =  (MatJ[4]*MatJ[8] - MatJ[5]*MatJ[7] - MatJ[4]*MatJ[11] + MatJ[5]*MatJ[10] + MatJ[7]*MatJ[11] - MatJ[8]*MatJ[10])/(MatJ[0]*MatJ[4]*MatJ[8] - MatJ[0]*MatJ[5]*MatJ[7] - MatJ[1]*MatJ[3]*MatJ[8] + MatJ[1]*MatJ[5]*MatJ[6] + MatJ[2]*MatJ[3]*MatJ[7] - MatJ[2]*MatJ[4]*MatJ[6] - MatJ[0]*MatJ[4]*MatJ[11] + MatJ[0]*MatJ[5]*MatJ[10] + MatJ[1]*MatJ[3]*MatJ[11] - MatJ[1]*MatJ[5]*MatJ[9] - MatJ[2]*MatJ[3]*MatJ[10] + MatJ[2]*MatJ[4]*MatJ[9] + MatJ[0]*MatJ[7]*MatJ[11] - MatJ[0]*MatJ[8]*MatJ[10] - MatJ[1]*MatJ[6]*MatJ[11] + MatJ[1]*MatJ[8]*MatJ[9] + MatJ[2]*MatJ[6]*MatJ[10] - MatJ[2]*MatJ[7]*MatJ[9] - MatJ[3]*MatJ[7]*MatJ[11] + MatJ[3]*MatJ[8]*MatJ[10] + MatJ[4]*MatJ[6]*MatJ[11] - MatJ[4]*MatJ[8]*MatJ[9] - MatJ[5]*MatJ[6]*MatJ[10] + MatJ[5]*MatJ[7]*MatJ[9]);
		invJ[1] = -(MatJ[1]*MatJ[8] - MatJ[2]*MatJ[7] - MatJ[1]*MatJ[11] + MatJ[2]*MatJ[10] + MatJ[7]*MatJ[11] - MatJ[8]*MatJ[10])/(MatJ[0]*MatJ[4]*MatJ[8] - MatJ[0]*MatJ[5]*MatJ[7] - MatJ[1]*MatJ[3]*MatJ[8] + MatJ[1]*MatJ[5]*MatJ[6] + MatJ[2]*MatJ[3]*MatJ[7] - MatJ[2]*MatJ[4]*MatJ[6] - MatJ[0]*MatJ[4]*MatJ[11] + MatJ[0]*MatJ[5]*MatJ[10] + MatJ[1]*MatJ[3]*MatJ[11] - MatJ[1]*MatJ[5]*MatJ[9] - MatJ[2]*MatJ[3]*MatJ[10] + MatJ[2]*MatJ[4]*MatJ[9] + MatJ[0]*MatJ[7]*MatJ[11] - MatJ[0]*MatJ[8]*MatJ[10] - MatJ[1]*MatJ[6]*MatJ[11] + MatJ[1]*MatJ[8]*MatJ[9] + MatJ[2]*MatJ[6]*MatJ[10] - MatJ[2]*MatJ[7]*MatJ[9] - MatJ[3]*MatJ[7]*MatJ[11] + MatJ[3]*MatJ[8]*MatJ[10] + MatJ[4]*MatJ[6]*MatJ[11] - MatJ[4]*MatJ[8]*MatJ[9] - MatJ[5]*MatJ[6]*MatJ[10] + MatJ[5]*MatJ[7]*MatJ[9]);
		invJ[2] =  (MatJ[1]*MatJ[5] - MatJ[2]*MatJ[4] - MatJ[1]*MatJ[11] + MatJ[2]*MatJ[10] + MatJ[4]*MatJ[11] - MatJ[5]*MatJ[10])/(MatJ[0]*MatJ[4]*MatJ[8] - MatJ[0]*MatJ[5]*MatJ[7] - MatJ[1]*MatJ[3]*MatJ[8] + MatJ[1]*MatJ[5]*MatJ[6] + MatJ[2]*MatJ[3]*MatJ[7] - MatJ[2]*MatJ[4]*MatJ[6] - MatJ[0]*MatJ[4]*MatJ[11] + MatJ[0]*MatJ[5]*MatJ[10] + MatJ[1]*MatJ[3]*MatJ[11] - MatJ[1]*MatJ[5]*MatJ[9] - MatJ[2]*MatJ[3]*MatJ[10] + MatJ[2]*MatJ[4]*MatJ[9] + MatJ[0]*MatJ[7]*MatJ[11] - MatJ[0]*MatJ[8]*MatJ[10] - MatJ[1]*MatJ[6]*MatJ[11] + MatJ[1]*MatJ[8]*MatJ[9] + MatJ[2]*MatJ[6]*MatJ[10] - MatJ[2]*MatJ[7]*MatJ[9] - MatJ[3]*MatJ[7]*MatJ[11] + MatJ[3]*MatJ[8]*MatJ[10] + MatJ[4]*MatJ[6]*MatJ[11] - MatJ[4]*MatJ[8]*MatJ[9] - MatJ[5]*MatJ[6]*MatJ[10] + MatJ[5]*MatJ[7]*MatJ[9]);
		invJ[3] = -(MatJ[1]*MatJ[5] - MatJ[2]*MatJ[4] - MatJ[1]*MatJ[8] + MatJ[2]*MatJ[7] + MatJ[4]*MatJ[8] - MatJ[5]*MatJ[7])/(MatJ[0]*MatJ[4]*MatJ[8] - MatJ[0]*MatJ[5]*MatJ[7] - MatJ[1]*MatJ[3]*MatJ[8] + MatJ[1]*MatJ[5]*MatJ[6] + MatJ[2]*MatJ[3]*MatJ[7] - MatJ[2]*MatJ[4]*MatJ[6] - MatJ[0]*MatJ[4]*MatJ[11] + MatJ[0]*MatJ[5]*MatJ[10] + MatJ[1]*MatJ[3]*MatJ[11] - MatJ[1]*MatJ[5]*MatJ[9] - MatJ[2]*MatJ[3]*MatJ[10] + MatJ[2]*MatJ[4]*MatJ[9] + MatJ[0]*MatJ[7]*MatJ[11] - MatJ[0]*MatJ[8]*MatJ[10] - MatJ[1]*MatJ[6]*MatJ[11] + MatJ[1]*MatJ[8]*MatJ[9] + MatJ[2]*MatJ[6]*MatJ[10] - MatJ[2]*MatJ[7]*MatJ[9] - MatJ[3]*MatJ[7]*MatJ[11] + MatJ[3]*MatJ[8]*MatJ[10] + MatJ[4]*MatJ[6]*MatJ[11] - MatJ[4]*MatJ[8]*MatJ[9] - MatJ[5]*MatJ[6]*MatJ[10] + MatJ[5]*MatJ[7]*MatJ[9]);
		invJ[4] = -(MatJ[3]*MatJ[8] - MatJ[5]*MatJ[6] - MatJ[3]*MatJ[11] + MatJ[5]*MatJ[9] + MatJ[6]*MatJ[11] - MatJ[8]*MatJ[9])/(MatJ[0]*MatJ[4]*MatJ[8] - MatJ[0]*MatJ[5]*MatJ[7] - MatJ[1]*MatJ[3]*MatJ[8] + MatJ[1]*MatJ[5]*MatJ[6] + MatJ[2]*MatJ[3]*MatJ[7] - MatJ[2]*MatJ[4]*MatJ[6] - MatJ[0]*MatJ[4]*MatJ[11] + MatJ[0]*MatJ[5]*MatJ[10] + MatJ[1]*MatJ[3]*MatJ[11] - MatJ[1]*MatJ[5]*MatJ[9] - MatJ[2]*MatJ[3]*MatJ[10] + MatJ[2]*MatJ[4]*MatJ[9] + MatJ[0]*MatJ[7]*MatJ[11] - MatJ[0]*MatJ[8]*MatJ[10] - MatJ[1]*MatJ[6]*MatJ[11] + MatJ[1]*MatJ[8]*MatJ[9] + MatJ[2]*MatJ[6]*MatJ[10] - MatJ[2]*MatJ[7]*MatJ[9] - MatJ[3]*MatJ[7]*MatJ[11] + MatJ[3]*MatJ[8]*MatJ[10] + MatJ[4]*MatJ[6]*MatJ[11] - MatJ[4]*MatJ[8]*MatJ[9] - MatJ[5]*MatJ[6]*MatJ[10] + MatJ[5]*MatJ[7]*MatJ[9]);
		invJ[5] =  (MatJ[0]*MatJ[8] - MatJ[2]*MatJ[6] - MatJ[0]*MatJ[11] + MatJ[2]*MatJ[9] + MatJ[6]*MatJ[11] - MatJ[8]*MatJ[9])/(MatJ[0]*MatJ[4]*MatJ[8] - MatJ[0]*MatJ[5]*MatJ[7] - MatJ[1]*MatJ[3]*MatJ[8] + MatJ[1]*MatJ[5]*MatJ[6] + MatJ[2]*MatJ[3]*MatJ[7] - MatJ[2]*MatJ[4]*MatJ[6] - MatJ[0]*MatJ[4]*MatJ[11] + MatJ[0]*MatJ[5]*MatJ[10] + MatJ[1]*MatJ[3]*MatJ[11] - MatJ[1]*MatJ[5]*MatJ[9] - MatJ[2]*MatJ[3]*MatJ[10] + MatJ[2]*MatJ[4]*MatJ[9] + MatJ[0]*MatJ[7]*MatJ[11] - MatJ[0]*MatJ[8]*MatJ[10] - MatJ[1]*MatJ[6]*MatJ[11] + MatJ[1]*MatJ[8]*MatJ[9] + MatJ[2]*MatJ[6]*MatJ[10] - MatJ[2]*MatJ[7]*MatJ[9] - MatJ[3]*MatJ[7]*MatJ[11] + MatJ[3]*MatJ[8]*MatJ[10] + MatJ[4]*MatJ[6]*MatJ[11] - MatJ[4]*MatJ[8]*MatJ[9] - MatJ[5]*MatJ[6]*MatJ[10] + MatJ[5]*MatJ[7]*MatJ[9]);
		invJ[6] = -(MatJ[0]*MatJ[5] - MatJ[2]*MatJ[3] - MatJ[0]*MatJ[11] + MatJ[2]*MatJ[9] + MatJ[3]*MatJ[11] - MatJ[5]*MatJ[9])/(MatJ[0]*MatJ[4]*MatJ[8] - MatJ[0]*MatJ[5]*MatJ[7] - MatJ[1]*MatJ[3]*MatJ[8] + MatJ[1]*MatJ[5]*MatJ[6] + MatJ[2]*MatJ[3]*MatJ[7] - MatJ[2]*MatJ[4]*MatJ[6] - MatJ[0]*MatJ[4]*MatJ[11] + MatJ[0]*MatJ[5]*MatJ[10] + MatJ[1]*MatJ[3]*MatJ[11] - MatJ[1]*MatJ[5]*MatJ[9] - MatJ[2]*MatJ[3]*MatJ[10] + MatJ[2]*MatJ[4]*MatJ[9] + MatJ[0]*MatJ[7]*MatJ[11] - MatJ[0]*MatJ[8]*MatJ[10] - MatJ[1]*MatJ[6]*MatJ[11] + MatJ[1]*MatJ[8]*MatJ[9] + MatJ[2]*MatJ[6]*MatJ[10] - MatJ[2]*MatJ[7]*MatJ[9] - MatJ[3]*MatJ[7]*MatJ[11] + MatJ[3]*MatJ[8]*MatJ[10] + MatJ[4]*MatJ[6]*MatJ[11] - MatJ[4]*MatJ[8]*MatJ[9] - MatJ[5]*MatJ[6]*MatJ[10] + MatJ[5]*MatJ[7]*MatJ[9]);
		invJ[7] =  (MatJ[0]*MatJ[5] - MatJ[2]*MatJ[3] - MatJ[0]*MatJ[8] + MatJ[2]*MatJ[6] + MatJ[3]*MatJ[8] - MatJ[5]*MatJ[6])/(MatJ[0]*MatJ[4]*MatJ[8] - MatJ[0]*MatJ[5]*MatJ[7] - MatJ[1]*MatJ[3]*MatJ[8] + MatJ[1]*MatJ[5]*MatJ[6] + MatJ[2]*MatJ[3]*MatJ[7] - MatJ[2]*MatJ[4]*MatJ[6] - MatJ[0]*MatJ[4]*MatJ[11] + MatJ[0]*MatJ[5]*MatJ[10] + MatJ[1]*MatJ[3]*MatJ[11] - MatJ[1]*MatJ[5]*MatJ[9] - MatJ[2]*MatJ[3]*MatJ[10] + MatJ[2]*MatJ[4]*MatJ[9] + MatJ[0]*MatJ[7]*MatJ[11] - MatJ[0]*MatJ[8]*MatJ[10] - MatJ[1]*MatJ[6]*MatJ[11] + MatJ[1]*MatJ[8]*MatJ[9] + MatJ[2]*MatJ[6]*MatJ[10] - MatJ[2]*MatJ[7]*MatJ[9] - MatJ[3]*MatJ[7]*MatJ[11] + MatJ[3]*MatJ[8]*MatJ[10] + MatJ[4]*MatJ[6]*MatJ[11] - MatJ[4]*MatJ[8]*MatJ[9] - MatJ[5]*MatJ[6]*MatJ[10] + MatJ[5]*MatJ[7]*MatJ[9]);
		invJ[8] =  (MatJ[3]*MatJ[7] - MatJ[4]*MatJ[6] - MatJ[3]*MatJ[10] + MatJ[4]*MatJ[9] + MatJ[6]*MatJ[10] - MatJ[7]*MatJ[9])/(MatJ[0]*MatJ[4]*MatJ[8] - MatJ[0]*MatJ[5]*MatJ[7] - MatJ[1]*MatJ[3]*MatJ[8] + MatJ[1]*MatJ[5]*MatJ[6] + MatJ[2]*MatJ[3]*MatJ[7] - MatJ[2]*MatJ[4]*MatJ[6] - MatJ[0]*MatJ[4]*MatJ[11] + MatJ[0]*MatJ[5]*MatJ[10] + MatJ[1]*MatJ[3]*MatJ[11] - MatJ[1]*MatJ[5]*MatJ[9] - MatJ[2]*MatJ[3]*MatJ[10] + MatJ[2]*MatJ[4]*MatJ[9] + MatJ[0]*MatJ[7]*MatJ[11] - MatJ[0]*MatJ[8]*MatJ[10] - MatJ[1]*MatJ[6]*MatJ[11] + MatJ[1]*MatJ[8]*MatJ[9] + MatJ[2]*MatJ[6]*MatJ[10] - MatJ[2]*MatJ[7]*MatJ[9] - MatJ[3]*MatJ[7]*MatJ[11] + MatJ[3]*MatJ[8]*MatJ[10] + MatJ[4]*MatJ[6]*MatJ[11] - MatJ[4]*MatJ[8]*MatJ[9] - MatJ[5]*MatJ[6]*MatJ[10] + MatJ[5]*MatJ[7]*MatJ[9]);
		invJ[9] =  -(MatJ[0]*MatJ[7] - MatJ[1]*MatJ[6] - MatJ[0]*MatJ[10] + MatJ[1]*MatJ[9] + MatJ[6]*MatJ[10] - MatJ[7]*MatJ[9])/(MatJ[0]*MatJ[4]*MatJ[8] - MatJ[0]*MatJ[5]*MatJ[7] - MatJ[1]*MatJ[3]*MatJ[8] + MatJ[1]*MatJ[5]*MatJ[6] + MatJ[2]*MatJ[3]*MatJ[7] - MatJ[2]*MatJ[4]*MatJ[6] - MatJ[0]*MatJ[4]*MatJ[11] + MatJ[0]*MatJ[5]*MatJ[10] + MatJ[1]*MatJ[3]*MatJ[11] - MatJ[1]*MatJ[5]*MatJ[9] - MatJ[2]*MatJ[3]*MatJ[10] + MatJ[2]*MatJ[4]*MatJ[9] + MatJ[0]*MatJ[7]*MatJ[11] - MatJ[0]*MatJ[8]*MatJ[10] - MatJ[1]*MatJ[6]*MatJ[11] + MatJ[1]*MatJ[8]*MatJ[9] + MatJ[2]*MatJ[6]*MatJ[10] - MatJ[2]*MatJ[7]*MatJ[9] - MatJ[3]*MatJ[7]*MatJ[11] + MatJ[3]*MatJ[8]*MatJ[10] + MatJ[4]*MatJ[6]*MatJ[11] - MatJ[4]*MatJ[8]*MatJ[9] - MatJ[5]*MatJ[6]*MatJ[10] + MatJ[5]*MatJ[7]*MatJ[9]);
		invJ[10] =  (MatJ[0]*MatJ[4] - MatJ[1]*MatJ[3] - MatJ[0]*MatJ[10] + MatJ[1]*MatJ[9] + MatJ[3]*MatJ[10] - MatJ[4]*MatJ[9])/(MatJ[0]*MatJ[4]*MatJ[8] - MatJ[0]*MatJ[5]*MatJ[7] - MatJ[1]*MatJ[3]*MatJ[8] + MatJ[1]*MatJ[5]*MatJ[6] + MatJ[2]*MatJ[3]*MatJ[7] - MatJ[2]*MatJ[4]*MatJ[6] - MatJ[0]*MatJ[4]*MatJ[11] + MatJ[0]*MatJ[5]*MatJ[10] + MatJ[1]*MatJ[3]*MatJ[11] - MatJ[1]*MatJ[5]*MatJ[9] - MatJ[2]*MatJ[3]*MatJ[10] + MatJ[2]*MatJ[4]*MatJ[9] + MatJ[0]*MatJ[7]*MatJ[11] - MatJ[0]*MatJ[8]*MatJ[10] - MatJ[1]*MatJ[6]*MatJ[11] + MatJ[1]*MatJ[8]*MatJ[9] + MatJ[2]*MatJ[6]*MatJ[10] - MatJ[2]*MatJ[7]*MatJ[9] - MatJ[3]*MatJ[7]*MatJ[11] + MatJ[3]*MatJ[8]*MatJ[10] + MatJ[4]*MatJ[6]*MatJ[11] - MatJ[4]*MatJ[8]*MatJ[9] - MatJ[5]*MatJ[6]*MatJ[10] + MatJ[5]*MatJ[7]*MatJ[9]);
		invJ[11] =  -(MatJ[0]*MatJ[4] - MatJ[1]*MatJ[3] - MatJ[0]*MatJ[7] + MatJ[1]*MatJ[6] + MatJ[3]*MatJ[7] - MatJ[4]*MatJ[6])/(MatJ[0]*MatJ[4]*MatJ[8] - MatJ[0]*MatJ[5]*MatJ[7] - MatJ[1]*MatJ[3]*MatJ[8] + MatJ[1]*MatJ[5]*MatJ[6] + MatJ[2]*MatJ[3]*MatJ[7] - MatJ[2]*MatJ[4]*MatJ[6] - MatJ[0]*MatJ[4]*MatJ[11] + MatJ[0]*MatJ[5]*MatJ[10] + MatJ[1]*MatJ[3]*MatJ[11] - MatJ[1]*MatJ[5]*MatJ[9] - MatJ[2]*MatJ[3]*MatJ[10] + MatJ[2]*MatJ[4]*MatJ[9] + MatJ[0]*MatJ[7]*MatJ[11] - MatJ[0]*MatJ[8]*MatJ[10] - MatJ[1]*MatJ[6]*MatJ[11] + MatJ[1]*MatJ[8]*MatJ[9] + MatJ[2]*MatJ[6]*MatJ[10] - MatJ[2]*MatJ[7]*MatJ[9] - MatJ[3]*MatJ[7]*MatJ[11] + MatJ[3]*MatJ[8]*MatJ[10] + MatJ[4]*MatJ[6]*MatJ[11] - MatJ[4]*MatJ[8]*MatJ[9] - MatJ[5]*MatJ[6]*MatJ[10] + MatJ[5]*MatJ[7]*MatJ[9]);

		double detJ =  fabs( (MatJ[0])* (MatJ[4])* (MatJ[8]) -  (MatJ[0])* (MatJ[5])* (MatJ[7]) -  (MatJ[1])* (MatJ[3])* (MatJ[8]) +  (MatJ[1])* (MatJ[5])* (MatJ[6]) +  (MatJ[2])* (MatJ[3])* (MatJ[7]) -  (MatJ[2])* (MatJ[4])* (MatJ[6]) -  (MatJ[0])* (MatJ[4])* (MatJ[11]) +  (MatJ[0])* (MatJ[5])* (MatJ[10]) +  (MatJ[1])* (MatJ[3])* (MatJ[11]) -  (MatJ[1])* (MatJ[5])* (MatJ[9]) -  (MatJ[2])* (MatJ[3])* (MatJ[10]) +  (MatJ[2])* (MatJ[4])* (MatJ[9]) +  (MatJ[0])* (MatJ[7])* (MatJ[11]) -  (MatJ[0])* (MatJ[8])* (MatJ[10]) -  (MatJ[1])* (MatJ[6])* (MatJ[11]) +  (MatJ[1])* (MatJ[8])* (MatJ[9]) +  (MatJ[2])* (MatJ[6])* (MatJ[10]) -  (MatJ[2])* (MatJ[7])* (MatJ[9]) -  (MatJ[3])* (MatJ[7])* (MatJ[11]) +  (MatJ[3])* (MatJ[8])* (MatJ[10]) +  (MatJ[4])* (MatJ[6])* (MatJ[11]) -  (MatJ[4])* (MatJ[8])* (MatJ[9]) -  (MatJ[5])* (MatJ[6])* (MatJ[10]) +  (MatJ[5])* (MatJ[7])* (MatJ[9]));

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, k, m, alpha, &invJ[0], m, &dN[gp_index][0], k, beta,&dND[0], k);

		for (int j = 0; j < 10; j++) {

			B[mapVecB[j]-1]    = dND[j];
			B[mapVecB[j+10]-1] = dND[j];
			B[mapVecB[j+20]-1] = dND[j];		

			B[mapVecB[j+30]-1] = dND[j+10];
			B[mapVecB[j+40]-1] = dND[j+10];
			B[mapVecB[j+50]-1] = dND[j+10];

			B[mapVecB[j+60]-1] = dND[j+20];
			B[mapVecB[j+70]-1] = dND[j+20];
			B[mapVecB[j+80]-1] = dND[j+20];

		}

		//C*B'
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m1, k1, n1, alpha, &MatC[0], m1, &B[0], k1, beta, &CB[0], m1);

		//Ke = Ke+(B*(C*B'))*dJ*WF;
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, k1, k1, n1, detJ*WeighFactor[gp_index], &B[0], k1, &CB[0], n1, alpha, &Ke[0], k1);

		for (int i = 0; i < 30; i++) {		
			fe[i]= fe[i] + (detJ*WeighFactor[gp_index])*(N[gp_index][ i / 3 ] * inertia[ i % 3]);
		}


		// dynamic 
		SEQ_VECTOR <double> NN (90, 0.0);

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 10; j++) {
				NN[mapVecN[10*i+j]-1] = N[gp_index][j];
			}
		}

		// Me = Me + WF*(DENS*dJ)*(NN*NN');
		double dense = 7.85e-9;
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 30, 30, 3, dense*detJ*WeighFactor[gp_index], &NN[0], 30, &NN[0], 30, alpha, &Me[0], 30);

	}
}

void FEM_Assembler::pyramid13_Elasticity(SEQ_VECTOR <double> & fe, SEQ_VECTOR <double> & Ke, SEQ_VECTOR <int> & mapVecB, SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N ,SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <double> & coordinates, SEQ_VECTOR <double> & inertia, double  ex, double  mi) {

	int gp_length  = WeighFactor.size(); 

	int matrixSize = coordinates.size();


	fill(Ke.begin(),Ke.end(),0); 
	fill(fe.begin(),fe.end(),0); 

	SEQ_VECTOR <double> MatC (6*6, 0.0); 

	double E   = ex/((1+mi)*(1-2*mi));
	double mi2 = E*(1.0 - mi); 
	double mi3 = E*(0.5 - mi);
	mi = E * mi;

	MatC[0]  = mi2;		MatC[1]  = mi; 		MatC[2]  = mi;
	MatC[6]  = mi; 		MatC[7]  = mi2;		MatC[8]  = mi;
	MatC[12] = mi;		MatC[13] = mi;		MatC[14] = mi2;
	MatC[21] = mi3;		MatC[28] = mi3;		MatC[35] = mi3;

	char trans   = 'N'; 
	double alpha =  1; 
	double beta  =  0; 
	int m = 3; 
	int n = 3;
	int k = 13; 

	int m1 = 6; 
	int n1 = 6;
	int k1 = 39; 

	SEQ_VECTOR <double> MatJ (9, 0);
	SEQ_VECTOR <double> invJ (9, 0);
	SEQ_VECTOR <double> dND (39, 0);
	SEQ_VECTOR <double> CB (234, 0);
	SEQ_VECTOR <double> B (234, 0);

	Ke.resize(matrixSize*matrixSize);
	fe.resize(matrixSize);

	fill(Ke.begin(),Ke.end(),0); 
	fill(fe.begin(),fe.end(),0); 


	for (int gp_index = 0; gp_index < gp_length; gp_index++) {

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, &dN[gp_index][0], k, &coordinates[0], n, beta, &MatJ[0], n);

		double detJ = fabs( MatJ[0] * MatJ[4] * MatJ[8] + 
			MatJ[1] * MatJ[5] * MatJ[6] + 
			MatJ[2] * MatJ[3] * MatJ[7] - 
			MatJ[2] * MatJ[4] * MatJ[6] -
			MatJ[1] * MatJ[3] * MatJ[8] -
			MatJ[0] * MatJ[5] * MatJ[7]);

		double detJx = 1/detJ;

		invJ[0] = detJx * (  MatJ[8] * MatJ[4] - MatJ[7] * MatJ[5] );
		invJ[1] = detJx * (- MatJ[8] * MatJ[1] + MatJ[7] * MatJ[2] );
		invJ[2] = detJx * (  MatJ[5] * MatJ[1] - MatJ[4] * MatJ[2] );
		invJ[3] = detJx * (- MatJ[8] * MatJ[3] + MatJ[6] * MatJ[5] );		
		invJ[4] = detJx * (  MatJ[8] * MatJ[0] - MatJ[6] * MatJ[2] );		
		invJ[5] = detJx * (- MatJ[5] * MatJ[0] + MatJ[3] * MatJ[2] );		
		invJ[6] = detJx * (  MatJ[7] * MatJ[3] - MatJ[6] * MatJ[4] );		
		invJ[7] = detJx * (- MatJ[7] * MatJ[0] + MatJ[6] * MatJ[1] );
		invJ[8] = detJx * (  MatJ[4] * MatJ[0] - MatJ[3] * MatJ[1] );

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, k, n, alpha, &invJ[0], n, &dN[gp_index][0], k, beta,&dND[0], k);

		for (int j = 0; j < 13; j++) {

			B[mapVecB[j]-1]    = dND[j];
			B[mapVecB[j+13]-1]  = dND[j];
			B[mapVecB[j+26]-1] = dND[j];


			B[mapVecB[j+39]-1] = dND[j+13];
			B[mapVecB[j+52]-1] = dND[j+13];
			B[mapVecB[j+65]-1] = dND[j+13];


			B[mapVecB[j+78]-1] = dND[j+26];
			B[mapVecB[j+91]-1] = dND[j+26];
			B[mapVecB[j+104]-1] = dND[j+26];

		}

		//C*B'
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m1, k1, n1, alpha, &MatC[0], m1, &B[0], k1, beta, &CB[0], m1);

		//Ke = Ke+(B*(C*B'))*dJ*WF;
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, k1, k1, n1, detJ*WeighFactor[gp_index], &B[0], k1, &CB[0], n1, alpha, &Ke[0], k1);

		for (int i = 0; i < 45; i++) {		
			fe[i]= fe[i] + (detJ*WeighFactor[gp_index])*(N[gp_index][ i / 3 ] * inertia[ i % 3]);
		}

	}
}

void FEM_Assembler::prisma15_Elasticity(SEQ_VECTOR <double> & fe, SEQ_VECTOR <double> & Ke, SEQ_VECTOR <int> & mapVecB, SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N ,SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <double> & coordinates, SEQ_VECTOR <double> & inertia, double  ex, double  mi) {

	int gp_length  = WeighFactor.size(); 

	int matrixSize = coordinates.size();


	fill(Ke.begin(),Ke.end(),0); 
	fill(fe.begin(),fe.end(),0); 

	SEQ_VECTOR <double> MatC (6*6, 0.0); 

	double E   = ex/((1+mi)*(1-2*mi));
	double mi2 = E*(1.0 - mi); 
	double mi3 = E*(0.5 - mi);
	mi = E * mi;

	MatC[0]  = mi2;		MatC[1]  = mi; 		MatC[2]  = mi;
	MatC[6]  = mi; 		MatC[7]  = mi2;		MatC[8]  = mi;
	MatC[12] = mi;		MatC[13] = mi;		MatC[14] = mi2;
	MatC[21] = mi3;		MatC[28] = mi3;		MatC[35] = mi3;

	char trans   = 'N'; 
	double alpha =  1;
	double beta  =  0; 
	int m = 3; 
	int n = 3;
	int k = 15; 

	int m1 = 6; 
	int n1 = 6;
	int k1 = 45; 

	SEQ_VECTOR <double> MatJ (9, 0);
	SEQ_VECTOR <double> invJ (9, 0);
	SEQ_VECTOR <double> dND (45, 0);
	SEQ_VECTOR <double> CB (270, 0);
	SEQ_VECTOR <double> B (270, 0);

	Ke.resize(matrixSize*matrixSize);
	fe.resize(matrixSize);

	fill(Ke.begin(),Ke.end(),0); 
	fill(fe.begin(),fe.end(),0); 


	for (int gp_index = 0; gp_index < gp_length; gp_index++) {

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, &dN[gp_index][0], k, &coordinates[0], n, beta, &MatJ[0], n);

		double detJ = fabs( MatJ[0] * MatJ[4] * MatJ[8] + 
			MatJ[1] * MatJ[5] * MatJ[6] + 
			MatJ[2] * MatJ[3] * MatJ[7] - 
			MatJ[2] * MatJ[4] * MatJ[6] -
			MatJ[1] * MatJ[3] * MatJ[8] -
			MatJ[0] * MatJ[5] * MatJ[7]);

		double detJx = 1/detJ;

		invJ[0] = detJx * (  MatJ[8] * MatJ[4] - MatJ[7] * MatJ[5] );
		invJ[1] = detJx * (- MatJ[8] * MatJ[1] + MatJ[7] * MatJ[2] );
		invJ[2] = detJx * (  MatJ[5] * MatJ[1] - MatJ[4] * MatJ[2] );
		invJ[3] = detJx * (- MatJ[8] * MatJ[3] + MatJ[6] * MatJ[5] );		
		invJ[4] = detJx * (  MatJ[8] * MatJ[0] - MatJ[6] * MatJ[2] );		
		invJ[5] = detJx * (- MatJ[5] * MatJ[0] + MatJ[3] * MatJ[2] );		
		invJ[6] = detJx * (  MatJ[7] * MatJ[3] - MatJ[6] * MatJ[4] );
		invJ[7] = detJx * (- MatJ[7] * MatJ[0] + MatJ[6] * MatJ[1] );
		invJ[8] = detJx * (  MatJ[4] * MatJ[0] - MatJ[3] * MatJ[1] );

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, k, n, alpha, &invJ[0], n, &dN[gp_index][0], k, beta,&dND[0], k);

		for (int j = 0; j < 15; j++) {

			B[mapVecB[j]-1]    = dND[j];
			B[mapVecB[j+15]-1]  = dND[j];
			B[mapVecB[j+30]-1] = dND[j];


			B[mapVecB[j+45]-1] = dND[j+15];
			B[mapVecB[j+60]-1] = dND[j+15];
			B[mapVecB[j+75]-1] = dND[j+15];


			B[mapVecB[j+90]-1] = dND[j+30];
			B[mapVecB[j+105]-1] = dND[j+30];
			B[mapVecB[j+120]-1] = dND[j+30];

		}

		//C*B'
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m1, k1, n1, alpha, &MatC[0], m1, &B[0], k1, beta, &CB[0], m1);

		//Ke = Ke+(B*(C*B'))*dJ*WF;
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, k1, k1, n1, detJ*WeighFactor[gp_index], &B[0], k1, &CB[0], n1, alpha, &Ke[0], k1);

		for (int i = 0; i < 45; i++) {		
			fe[i]= fe[i] + (detJ*WeighFactor[gp_index])*(N[gp_index][ i / 3 ] * inertia[ i % 3]);
		}

	}
}

void FEM_Assembler::brick20_Elasticity(SEQ_VECTOR <double> & fe, SEQ_VECTOR <double> & Ke, SEQ_VECTOR <int> & mapVecB, SEQ_VECTOR<SEQ_VECTOR <double> > & dN , SEQ_VECTOR<SEQ_VECTOR <double> > & N ,SEQ_VECTOR <double> & WeighFactor, SEQ_VECTOR <double> & coordinates, SEQ_VECTOR <double> & inertia, double  ex, double  mi) {

	int gp_length  = WeighFactor.size(); 

	int matrixSize = coordinates.size();


	fill(Ke.begin(),Ke.end(),0); 
	fill(fe.begin(),fe.end(),0); 

	SEQ_VECTOR <double> MatC (6*6, 0.0); 

	double E   = ex/((1+mi)*(1-2*mi));
	double mi2 = E*(1.0 - mi); 
	double mi3 = E*(0.5 - mi);
	mi = E * mi;

	MatC[0]  = mi2;		MatC[1]  = mi; 		MatC[2]  = mi;
	MatC[6]  = mi; 		MatC[7]  = mi2;		MatC[8]  = mi;
	MatC[12] = mi;		MatC[13] = mi;		MatC[14] = mi2;
	MatC[21] = mi3;		MatC[28] = mi3;		MatC[35] = mi3;

	char trans   = 'N'; 
	double alpha =  1; 
	double beta  =  0; 
	int m = 3; 
	int n = 3;
	int k = 20; 

	int m1 = 6; 
	int n1 = 6;
	int k1 = 60; 

	SEQ_VECTOR <double> MatJ (9, 0);
	SEQ_VECTOR <double> invJ (9, 0);
	SEQ_VECTOR <double> dND (60, 0);
	SEQ_VECTOR <double> CB (360, 0);
	SEQ_VECTOR <double> B (360, 0);

	Ke.resize(matrixSize*matrixSize);
	fe.resize(matrixSize);

	fill(Ke.begin(),Ke.end(),0); 
	fill(fe.begin(),fe.end(),0); 




	for (int gp_index = 0; gp_index < gp_length; gp_index++) {

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, &dN[gp_index][0], k, &coordinates[0], n, beta, &MatJ[0], n);

		double detJ = fabs( MatJ[0] * MatJ[4] * MatJ[8] + 
			MatJ[1] * MatJ[5] * MatJ[6] + 
			MatJ[2] * MatJ[3] * MatJ[7] - 
			MatJ[2] * MatJ[4] * MatJ[6] -
			MatJ[1] * MatJ[3] * MatJ[8] -
			MatJ[0] * MatJ[5] * MatJ[7]);

		double detJx = 1/detJ;

		invJ[0] = detJx * (  MatJ[8] * MatJ[4] - MatJ[7] * MatJ[5] );
		invJ[1] = detJx * (- MatJ[8] * MatJ[1] + MatJ[7] * MatJ[2] );
		invJ[2] = detJx * (  MatJ[5] * MatJ[1] - MatJ[4] * MatJ[2] );
		invJ[3] = detJx * (- MatJ[8] * MatJ[3] + MatJ[6] * MatJ[5] );		
		invJ[4] = detJx * (  MatJ[8] * MatJ[0] - MatJ[6] * MatJ[2] );		
		invJ[5] = detJx * (- MatJ[5] * MatJ[0] + MatJ[3] * MatJ[2] );		
		invJ[6] = detJx * (  MatJ[7] * MatJ[3] - MatJ[6] * MatJ[4] );		
		invJ[7] = detJx * (- MatJ[7] * MatJ[0] + MatJ[6] * MatJ[1] );
		invJ[8] = detJx * (  MatJ[4] * MatJ[0] - MatJ[3] * MatJ[1] );

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, k, n, alpha, &invJ[0], n, &dN[gp_index][0], k, beta,&dND[0], k);

		for (int j = 0; j < 20; j++) {

			B[mapVecB[j]-1]    = dND[j];
			B[mapVecB[j+20]-1]  = dND[j];
			B[mapVecB[j+40]-1] = dND[j];


			B[mapVecB[j+60]-1] = dND[j+20];
			B[mapVecB[j+80]-1] = dND[j+20];
			B[mapVecB[j+100]-1] = dND[j+20];


			B[mapVecB[j+120]-1] = dND[j+40];
			B[mapVecB[j+140]-1] = dND[j+40];
			B[mapVecB[j+160]-1] = dND[j+40];

		}

		//C*B'
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m1, k1, n1, alpha, &MatC[0], m1, &B[0], k1, beta, &CB[0], m1);

		//Ke = Ke+(B*(C*B'))*dJ*WF;
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, k1, k1, n1, detJ*WeighFactor[gp_index], &B[0], k1, &CB[0], n1, alpha, &Ke[0], k1);

		for (int i = 0; i < 60; i++) {		
			fe[i]= fe[i] + (detJ*WeighFactor[gp_index])*(N[gp_index][ i / 3 ] * inertia[ i % 3]);
		}

	}
}

void FEM_Assembler::assemble_matrix(SparseMatrix & K, SEQ_VECTOR <double>  & f_global ){
	SparseMatrix M; 
	assemble_matrix(K, M, f_global);
}

void FEM_Assembler::assemble_matrix(SparseMatrix & K, SparseMatrix & M, SEQ_VECTOR <double>  & f_global ){

	SEQ_VECTOR<int> iii_tetra4;

	SEQ_VECTOR<SEQ_VECTOR <double> > dN_pyramid5; 
	SEQ_VECTOR<SEQ_VECTOR <double> >  N_pyramid5;
	SEQ_VECTOR <double> WeighFactor_pyramid5;
	SEQ_VECTOR<int> mapVecB_pyramid5;	
	SEQ_VECTOR<int> iii_pyramid5;

	SEQ_VECTOR<SEQ_VECTOR <double> > dN_prism6; 
	SEQ_VECTOR<SEQ_VECTOR <double> >  N_prism6;
	SEQ_VECTOR <double> WeighFactor_prism6;
	SEQ_VECTOR<int> mapVecB_prism6;	
	SEQ_VECTOR<int> iii_prism6;;

	SEQ_VECTOR<SEQ_VECTOR <double> > dN_brick8; 
	SEQ_VECTOR<SEQ_VECTOR <double> >  N_brick8;
	SEQ_VECTOR <double> WeighFactor_brick8;
	SEQ_VECTOR<int> mapVecB_brick8;	
	SEQ_VECTOR<int> iii_brick8;
	// dynamic 
	SEQ_VECTOR<int> mapVecN_brick8;	

	// meshType = 0 -- BILINEAR MESH
	if ( meshType == 0 ){

		int iii_tetra4_temp[] = {0, 1, 2, 4};
		iii_tetra4.assign(iii_tetra4_temp, iii_tetra4_temp + 4);

		////////////////////////////
		//PYRAMID5 BASIS
		pyramid5_basis (dN_pyramid5 , N_pyramid5 , WeighFactor_pyramid5 , mapVecB_pyramid5);
		////////////////////////////
		int iii_pyramid5_temp[] = {0, 1, 2, 3, 4};
		iii_pyramid5.assign(iii_pyramid5_temp, iii_pyramid5_temp + 5);


		////////////////////////////
		//PRISMA6 BASIS		
		prism6_basis (dN_prism6 , N_prism6 , WeighFactor_prism6 , mapVecB_prism6);
		////////////////////////////
		int iii_prism6_temp[] = {0, 1, 2, 4, 5, 6};
		iii_prism6.assign(iii_prism6_temp, iii_prism6_temp + 6);

		////////////////////////////
		//BRICK8 BASIS		
		brick8_basis (dN_brick8 , N_brick8 ,WeighFactor_brick8,mapVecB_brick8, mapVecN_brick8);
		int iii_brick8_temp[] = {0, 1, 2, 3, 4, 5, 6, 7};
		iii_brick8.assign(iii_brick8_temp, iii_brick8_temp + 8);
		////////////////////////////
	}


	SEQ_VECTOR<SEQ_VECTOR <double> > dN_tetra10; 
	SEQ_VECTOR<SEQ_VECTOR <double> >  N_tetra10;
	SEQ_VECTOR <double> WeighFactor_tetra10;
	SEQ_VECTOR<int> mapVecB_tetra10;
	SEQ_VECTOR<int> iii_tetra10;
	// dynamic 
	SEQ_VECTOR<int> mapVecN_tetra10;



	SEQ_VECTOR<SEQ_VECTOR <double> > dN_pyramid13; 
	SEQ_VECTOR<SEQ_VECTOR <double> >  N_pyramid13;
	SEQ_VECTOR <double> WeighFactor_pyramid13;
	SEQ_VECTOR<int> mapVecB_pyramid13;	
	SEQ_VECTOR<int> iii_pyramid13;

	SEQ_VECTOR<SEQ_VECTOR <double> > dN_prisma15; 
	SEQ_VECTOR<SEQ_VECTOR <double> >  N_prisma15;
	SEQ_VECTOR <double> WeighFactor_prisma15;
	SEQ_VECTOR<int> mapVecB_prisma15;	
	SEQ_VECTOR<int> iii_prisma15;

	SEQ_VECTOR<SEQ_VECTOR <double> > dN_brick20; 
	SEQ_VECTOR<SEQ_VECTOR <double> >  N_brick20;
	SEQ_VECTOR <double> WeighFactor_brick20;
	SEQ_VECTOR<int> mapVecB_brick20;
	SEQ_VECTOR<int> iii_brick20;

	// meshType = 1 -- QUADRI MESH (midpoints)
	if ( meshType == 1 ){
		////////////////////////////
		//TETRA10 BASIS		
		tetrahedra10_basis (dN_tetra10 , N_tetra10 , WeighFactor_tetra10 , mapVecB_tetra10 , mapVecN_tetra10);
		////////////////////////////	
		int iii_tetra10_temp[] = { 0,     1,     2,     4,     8,    9,    11,    16,    17,    18};
		iii_tetra10.assign(iii_tetra10_temp, iii_tetra10_temp + 10);

		////////////////////////////
		//PYRAMID13 BASIS		
		pyramid13_basis (dN_pyramid13 , N_pyramid13 , WeighFactor_pyramid13 , mapVecB_pyramid13);
		////////////////////////////
		int iii_pyramid13_temp[] = {0, 1, 2, 3, 4, 8, 9, 10, 11, 16, 17, 18, 19};
		iii_pyramid13.assign(iii_pyramid13_temp, iii_pyramid13_temp + 13);
		////////////////////////////
		//PRISMA15 BASIS		
		prisma15_basis (dN_prisma15 , N_prisma15 , WeighFactor_prisma15 , mapVecB_prisma15);
		////////////////////////////
		int iii_prisma15_temp[] = {0, 1, 2, 4, 5, 6, 8, 9, 11, 12, 13, 15, 16, 17, 18};
		iii_prisma15.assign(iii_prisma15_temp, iii_prisma15_temp + 15);
		////////////////////////////
		//BRICK20 BASIS		
		brick20_basis (dN_brick20 , N_brick20 , WeighFactor_brick20 , mapVecB_brick20);
		////////////////////////////
		int iii_brick20_temp[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
		iii_brick20.assign(iii_brick20_temp, iii_brick20_temp + 20);
	}


	int nK; 

	SEQ_VECTOR <double> inertia (3, 0.0); 

	inertia[0]=0.0; inertia[1]=0.0; inertia[2]=9810.0*7.85e-9;

	double ex; 
	double mi; 

	SEQ_VECTOR <double> fe; 
	SEQ_VECTOR <double> Ke; 
	SEQ_VECTOR <double> Me; 
	SEQ_VECTOR <int> DOFs; 
	SEQ_VECTOR <double> coordinates;

	DOFs.reserve(60);
	coordinates.reserve(60);
	Ke.reserve(3600);
	Me.reserve(3600);
	fe.reserve(60);

	ex = 2.1e5;
	mi = 0.3; 

	nK = all_coordinates.size() * 3;

	f_global.resize(nK);

	SEQ_VECTOR <SEQ_VECTOR <IV_elem> > IND_IV (nK);
	SEQ_VECTOR <SEQ_VECTOR <int> > IND_I (nK);
	SEQ_VECTOR <SEQ_VECTOR <double> > IND_V (nK);

	// dynamic 
	SEQ_VECTOR <SEQ_VECTOR <IV_elem> > IND_IV_Me (nK);

	int INT_SIZE_K = 0;
	int INT_SIZE_N = 0; 
	int IND_elementType = 0;

	for (int nElement_Count = 0; nElement_Count < all_elements.size(); nElement_Count++) {

		IND_elementType = elementType[nElement_Count];



		///////////////////////////////////////////
		// LOCAL STIFFNES MATRIX tetra4
		if (IND_elementType == 4){

			INT_SIZE_K = 12;
			INT_SIZE_N = 4; 

			coordinates.resize(INT_SIZE_K);

			for (int i = 0; i < INT_SIZE_N; i++) {
				coordinates[i*3] = all_coordinates[all_elements[nElement_Count][iii_tetra4[i]] -1   ][0];
				coordinates[i*3+1] = all_coordinates[all_elements[nElement_Count][iii_tetra4[i]] -1 ][1];
				coordinates[i*3+2] = all_coordinates[all_elements[nElement_Count][iii_tetra4[i]] -1 ][2];
			}

			DOFs.resize(INT_SIZE_K*3); 

			for (int i=0; i< INT_SIZE_N; i++){
				DOFs[i*3] = all_elements[nElement_Count][iii_tetra4[i]]*3 -2;
				DOFs[i*3 +1] = all_elements[nElement_Count][iii_tetra4[i]]*3 -1;
				DOFs[i*3 +2] = all_elements[nElement_Count][iii_tetra4[i]]*3;
			}

			//tetra4_Elasticity(fe, Ke, mapVecB_pyramid5, dN_pyramid5, N_pyramid5, WeighFactor_pyramid5, coordinates, inertia, ex, mi);

		}
		///////////////////////////////////////////END LOCAL STIFFNES MATRIX tetra4 



		///////////////////////////////////////////
		// LOCAL STIFFNES MATRIX pyramid_5
		if (IND_elementType == 5){

			INT_SIZE_K = 15;
			INT_SIZE_N = 5; 

			coordinates.resize(INT_SIZE_K); 

			for (int i = 0; i < INT_SIZE_N; i++) {
				coordinates[i*3] = all_coordinates[all_elements[nElement_Count][iii_pyramid5[i]] -1   ][0];
				coordinates[i*3+1] = all_coordinates[all_elements[nElement_Count][iii_pyramid5[i]] -1 ][1];
				coordinates[i*3+2] = all_coordinates[all_elements[nElement_Count][iii_pyramid5[i]] -1 ][2];
			}

			DOFs.resize(INT_SIZE_K*3); 

			for (int i=0; i< INT_SIZE_N; i++){
				DOFs[i*3] = all_elements[nElement_Count][iii_pyramid5[i]]*3 -2;
				DOFs[i*3 +1] = all_elements[nElement_Count][iii_pyramid5[i]]*3 -1;
				DOFs[i*3 +2] = all_elements[nElement_Count][iii_pyramid5[i]]*3;
			}

			pyramid5_Elasticity(fe, Ke, mapVecB_pyramid5, dN_pyramid5, N_pyramid5, WeighFactor_pyramid5, coordinates, inertia, ex, mi);

		}
		///////////////////////////////////////////END LOCAL STIFFNES MATRIX pyramid_5 


		///////////////////////////////////////////
		// LOCAL STIFFNES MATRIX prism_6 
		if (IND_elementType == 6){

			INT_SIZE_K = 18;
			INT_SIZE_N = 6; 

			coordinates.resize(INT_SIZE_K); 

			for (int i = 0; i < INT_SIZE_N; i++) {
				coordinates[i*3] = all_coordinates[all_elements[nElement_Count][iii_prism6[i]] -1   ][0];
				coordinates[i*3+1] = all_coordinates[all_elements[nElement_Count][iii_prism6[i]] -1 ][1];
				coordinates[i*3+2] = all_coordinates[all_elements[nElement_Count][iii_prism6[i]] -1 ][2];
			}

			DOFs.resize(INT_SIZE_K*3); 

			for (int i=0; i< INT_SIZE_N; i++){
				DOFs[i*3] = all_elements[nElement_Count][iii_prism6[i]]*3 -2;
				DOFs[i*3 +1] = all_elements[nElement_Count][iii_prism6[i]]*3 -1;
				DOFs[i*3 +2] = all_elements[nElement_Count][iii_prism6[i]]*3;
			}

			prisma6_Elasticity(fe, Ke, mapVecB_prism6, dN_prism6, N_prism6, WeighFactor_prism6, coordinates, inertia, ex, mi);

		}
		///////////////////////////////////////////END LOCAL STIFFNES MATRIX prism6 





		///////////////////////////////////////////
		// LOCAL STIFFNES MATRIX brick_8 
		if (IND_elementType == 8){

			INT_SIZE_K = 24;
			INT_SIZE_N = 8; 

			coordinates.resize(INT_SIZE_K); 

			for (int i = 0; i < INT_SIZE_N; i++) {
				coordinates[i*3] = all_coordinates[all_elements[nElement_Count][iii_brick8[i]] -1   ][0];
				coordinates[i*3+1] = all_coordinates[all_elements[nElement_Count][iii_brick8[i]] -1 ][1];
				coordinates[i*3+2] = all_coordinates[all_elements[nElement_Count][iii_brick8[i]] -1 ][2];
			}

			DOFs.resize(INT_SIZE_K); 

			for (int i=0; i< INT_SIZE_N; i++){
				DOFs[i*3] = all_elements[nElement_Count][iii_brick8[i]]*3 -2;
				DOFs[i*3 +1] = all_elements[nElement_Count][iii_brick8[i]]*3 -1;
				DOFs[i*3 +2] = all_elements[nElement_Count][iii_brick8[i]]*3;
			}

			brick8_Elasticity(fe, Ke, Me, mapVecB_brick8, dN_brick8, N_brick8, mapVecN_brick8, WeighFactor_brick8, coordinates, inertia, ex, mi);

		}
		///////////////////////////////////////////END LOCAL STIFFNES MATRIX brick8 


		///////////////////////////////////////////
		// LOCAL STIFFNES MATRIX tetrahedra_10 
		if (IND_elementType == 10){

			INT_SIZE_K = 30;
			INT_SIZE_N = 10; 

			coordinates.resize(INT_SIZE_K); 

			for (int i = 0; i < INT_SIZE_N; i++) {
				coordinates[i*3] = all_coordinates[all_elements[nElement_Count][iii_tetra10[i]] -1   ][0];
				coordinates[i*3+1] = all_coordinates[all_elements[nElement_Count][iii_tetra10[i]] -1 ][1];
				coordinates[i*3+2] = all_coordinates[all_elements[nElement_Count][iii_tetra10[i]] -1 ][2];
			}

			DOFs.resize(INT_SIZE_K*3); 

			for (int i=0; i< INT_SIZE_N; i++){
				DOFs[i*3] = all_elements[nElement_Count][iii_tetra10[i]]*3 -2;
				DOFs[i*3 +1] = all_elements[nElement_Count][iii_tetra10[i]]*3 -1;
				DOFs[i*3 +2] = all_elements[nElement_Count][iii_tetra10[i]]*3;
			}

			tetrahedra10_Elasticity(fe, Ke, Me, mapVecB_tetra10, dN_tetra10, N_tetra10, mapVecN_tetra10, WeighFactor_tetra10, coordinates, inertia, ex, mi);

		}
		///////////////////////////////////////////END LOCAL STIFFNES MATRIX tetrahedra_10 


		///////////////////////////////////////////
		// LOCAL STIFFNES MATRIX pyramid13 
		if (IND_elementType == 13){

			INT_SIZE_K = 49;
			INT_SIZE_N = 13; 

			coordinates.resize(INT_SIZE_K); 

			for (int i = 0; i < INT_SIZE_N; i++) {
				coordinates[i*3] = all_coordinates[all_elements[nElement_Count][iii_pyramid13[i]] -1   ][0];
				coordinates[i*3+1] = all_coordinates[all_elements[nElement_Count][iii_pyramid13[i]] -1 ][1];
				coordinates[i*3+2] = all_coordinates[all_elements[nElement_Count][iii_pyramid13[i]] -1 ][2];
			}

			DOFs.resize(INT_SIZE_K*3); 

			for (int i=0; i< INT_SIZE_N; i++){
				DOFs[i*3] = all_elements[nElement_Count][iii_pyramid13[i]]*3 -2;
				DOFs[i*3 +1] = all_elements[nElement_Count][iii_pyramid13[i]]*3 -1;
				DOFs[i*3 +2] = all_elements[nElement_Count][iii_pyramid13[i]]*3;
			}

			pyramid13_Elasticity (fe, Ke, mapVecB_pyramid13, dN_pyramid13, N_pyramid13, WeighFactor_pyramid13, coordinates, inertia, ex, mi);
		}
		///////////////////////////////////////////END LOCAL STIFFNES MATRIX pyramid13


		///////////////////////////////////////////
		// LOCAL STIFFNES MATRIX prisma15 
		if (IND_elementType == 15){

			INT_SIZE_K = 45;
			INT_SIZE_N = 15; 

			coordinates.resize(INT_SIZE_K); 

			for (int i = 0; i < INT_SIZE_N; i++) {
				coordinates[i*3] = all_coordinates[all_elements[nElement_Count][iii_prisma15[i]] -1   ][0];
				coordinates[i*3+1] = all_coordinates[all_elements[nElement_Count][iii_prisma15[i]] -1 ][1];
				coordinates[i*3+2] = all_coordinates[all_elements[nElement_Count][iii_prisma15[i]] -1 ][2];
			}

			DOFs.resize(INT_SIZE_K*3); 

			for (int i=0; i< INT_SIZE_N; i++){
				DOFs[i*3] = all_elements[nElement_Count][iii_prisma15[i]]*3 -2;
				DOFs[i*3 +1] = all_elements[nElement_Count][iii_prisma15[i]]*3 -1;
				DOFs[i*3 +2] = all_elements[nElement_Count][iii_prisma15[i]]*3;
			}

			prisma15_Elasticity (fe, Ke, mapVecB_prisma15, dN_prisma15, N_prisma15, WeighFactor_prisma15, coordinates, inertia, ex, mi);
		}
		///////////////////////////////////////////END LOCAL STIFFNES MATRIX prisma15 

		///////////////////////////////////////////
		// LOCAL STIFFNES MATRIX brick_20 
		if (IND_elementType == 20){

			INT_SIZE_K = 60;
			INT_SIZE_N = 20; 

			coordinates.resize(INT_SIZE_K); 

			for (int i = 0; i < INT_SIZE_N; i++) {
				coordinates[i*3] = all_coordinates[all_elements[nElement_Count][iii_brick20[i]] -1   ][0];
				coordinates[i*3+1] = all_coordinates[all_elements[nElement_Count][iii_brick20[i]] -1 ][1];
				coordinates[i*3+2] = all_coordinates[all_elements[nElement_Count][iii_brick20[i]] -1 ][2];
			}

			DOFs.resize(INT_SIZE_K*3); 

			for (int i=0; i< INT_SIZE_N; i++){
				DOFs[i*3] = all_elements[nElement_Count][iii_brick20[i]]*3 -2;
				DOFs[i*3 +1] = all_elements[nElement_Count][iii_brick20[i]]*3 -1;
				DOFs[i*3 +2] = all_elements[nElement_Count][iii_brick20[i]]*3;
			}

			brick20_Elasticity(fe, Ke, mapVecB_brick20, dN_brick20, N_brick20, WeighFactor_brick20, coordinates, inertia, ex, mi);

		}
		///////////////////////////////////////////END LOCAL STIFFNES MATRIX brick_20 

		int a = 9;

		for (int i=0; i < INT_SIZE_K; i++){
			for (int j=0; j < INT_SIZE_K; j++){
				if (Ke[i*INT_SIZE_K+j] != 0.0) {
					IND_IV   [DOFs[i] - 1].push_back(IV_elem (DOFs[j],Ke[i*INT_SIZE_K+j]));
				}
				
				if (USE_DYNAMIC == 1) {
					if (Me[i*INT_SIZE_K+j] != 0.0) {
						IND_IV_Me[DOFs[i] - 1].push_back(IV_elem (DOFs[j],Me[i*INT_SIZE_K+j])); // dynamic 
					}
				}
			}
			f_global[DOFs[i] - 1] = f_global[DOFs[i] - 1] + fe[i];
		}

		//// dynamic 
		//for (int i=0; i < INT_SIZE_K; i++){
		//	for (int j=0; j < INT_SIZE_K; j++){
		//		IND_IV_Me[DOFs[i] - 1].push_back(IV_elem (DOFs[j],Me[i*INT_SIZE_K+j]));
		//	}
		//}

	}


	std::SEQ_VECTOR<IV_elem>::iterator it;
	for (int i = 0; i < IND_IV.size(); i++) {
		sort(IND_IV[i].begin(), IND_IV[i].end(), comparator);

		for (int r = IND_IV[i].size() - 1; r > 0; r-- ) {
			if ( IND_IV[i][r-1].first == IND_IV[i][r].first)
				IND_IV[i][r-1].second += IND_IV[i][r].second;
		}

		it = unique(IND_IV[i].begin(), IND_IV[i].end(), unique_comparator);
		IND_IV[i].resize( distance(IND_IV[i].begin(),it) );
	}
		
	int row_index = 1; 
	K.CSR_I_row_indices.push_back(row_index); 
	for (int i = 0; i < IND_IV.size(); i++) {
		for (int j = 0; j < IND_IV[i].size(); j++) {
			K.CSR_J_col_indices.push_back(IND_IV[i][j].first);
			K.CSR_V_values.push_back(     IND_IV[i][j].second);
		}
		row_index += IND_IV[i].size(); 
		K.CSR_I_row_indices.push_back(row_index); 
	}

	K.cols = nK;
	K.rows = nK;
	K.nnz  = K.CSR_V_values.size(); 
	K.type = 'S';


	// dynamic 

	if (USE_DYNAMIC == 1) {
		std::SEQ_VECTOR<IV_elem>::iterator it2;
		for (int i = 0; i < IND_IV_Me.size(); i++) {
			sort(IND_IV_Me[i].begin(), IND_IV_Me[i].end(), comparator);

			for (int r = IND_IV_Me[i].size() - 1; r > 0; r-- ) {
				if ( IND_IV_Me[i][r-1].first == IND_IV_Me[i][r].first)
					IND_IV_Me[i][r-1].second += IND_IV_Me[i][r].second;
			}

			it2 = unique(IND_IV_Me[i].begin(), IND_IV_Me[i].end(), unique_comparator);
			IND_IV_Me[i].resize( distance(IND_IV_Me[i].begin(),it2) );
		}

		row_index = 1; 
		M.CSR_I_row_indices.push_back(row_index); 
		for (int i = 0; i < IND_IV_Me.size(); i++) {
			for (int j = 0; j < IND_IV_Me[i].size(); j++) {
				M.CSR_J_col_indices.push_back(IND_IV_Me[i][j].first);
				M.CSR_V_values.push_back(     IND_IV_Me[i][j].second);
			}
			row_index += IND_IV_Me[i].size(); 
			M.CSR_I_row_indices.push_back(row_index); 
		}

		M.cols = nK;
		M.rows = nK;
		M.nnz  = M.CSR_V_values.size(); 
		M.type = 'S';

	}

}
